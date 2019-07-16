/*********************************************************************
 * nl_solver.c
 *
 * Copyright (C) 2002 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Takaaki NAGATANI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: nl_solver.c 842 2006-12-04 06:21:20Z sasaoka $ */

#include <stdio.h>
#include <math.h>
#include "fem_solver.h"
#include "fem_struct.h"
#include "math_utl.h"
#include "rc.h"
#include "log_printf.h"
#include "nonzero_cg.h"

static RC update_NLsolve_iccg3(NODE_ARRAY node, ELEMENT_ARRAY element,
                               MATERIAL_PROP_ARRAY material,
                               PHYSICAL_PROP_ARRAY physical,
                               BC_ARRAY rest, int iter_num, BC_ARRAY force,
                               DEFAULT_BC_ARRAY b_force, DISP_ARRAY *disp,
                               STRESS_ARRAY *stress, STRAIN_ARRAY *strain,
                               int sol_ss_flag);
static int cal_linesearch_factor(int total_dof, int ls_num, double d_vect[],
                                 double r_vect[], double r_vect_old[],
                                 double *ls_fact);
static RC set_G_matrix(int dim, int node_num,
                       double **dNxyz, double **G_matrix);
static RC set_A_matrix(int dim, double **duvw, double **A_matrix);
static RC set_subM_matrix(int dim, double stress, double **M_matrix, int x,
                          int y);
static RC disp_add_node_location_matrix(ELEMENT *element,
                                        NODE_ARRAY node,
                                        DISP_ARRAY disp, double **matrix);


/* Update Lagrange 法による荷重増分法付き非線形解析ソルバ
 * ただし，微小ひずみ大変形理論に基づく定式化に従う       */
/* inc_numパラメータには荷重増分を行う回数を入力する
 * 荷重増分を行わない場合，inc_numには 1 を入れる         */
RC
fem_NLsolve_iccg3 (NODE_ARRAY node, ELEMENT_ARRAY element,
                   MATERIAL_PROP_ARRAY material,
                   PHYSICAL_PROP_ARRAY physical,
                   BC_ARRAY rest, int inc_num, BC_ARRAY force,
                   DEFAULT_BC_ARRAY b_force,
                   DISP_ARRAY *disp, STRESS_ARRAY *stress,
                   STRAIN_ARRAY *strain, int sol_ss_flag)
{
	int ii1;
	BC_ARRAY inc_force;
	DEFAULT_BC_ARRAY inc_b_force;

	RC_TRY( allocate_disp_array(0, disp) );
	RC_TRY( allocate_bc_array(force.size, &inc_force) );
	RC_TRY( allocate_default_bc_array(b_force.size, &inc_b_force) );

	if(force.size != 0){
		RC_TRY( copy_bc_array(force, &inc_force) );
		for(ii1=0; ii1<force.size; ii1++){
			inc_force.array[ii1].v.x = force.array[ii1].v.x / (double)inc_num;
			inc_force.array[ii1].v.y = force.array[ii1].v.y / (double)inc_num;
			inc_force.array[ii1].v.z = force.array[ii1].v.z / (double)inc_num;
		}
	}
	if( b_force.size != 0){
		RC_TRY( copy_default_bc_array(b_force, &inc_b_force) );
		for(ii1=0; ii1<b_force.size; ii1++){
			inc_b_force.array[ii1].grav.x = b_force.array[ii1].grav.x
			                                / (double)inc_num;
			inc_b_force.array[ii1].grav.y = b_force.array[ii1].grav.y
			                                / (double)inc_num;
			inc_b_force.array[ii1].grav.z = b_force.array[ii1].grav.z
			                                / (double)inc_num;
		}
	}

	for(ii1=1; ii1<(inc_num+1); ii1++){
		RC_TRY( log_printf(4, "inc num : %d/%d\n", ii1, inc_num) );
		if(ii1 < inc_num){
			RC_TRY( update_NLsolve_iccg3(node, element, material, physical,
			                             rest, ii1, inc_force, inc_b_force,
			                             disp, NULL, NULL, 0) );
		}else{
			RC_TRY( update_NLsolve_iccg3(node, element, material, physical,
			                             rest, ii1, inc_force, inc_b_force,
			                             disp, stress, strain, sol_ss_flag) );
		}
		RC_TRY( log_printf(4, "\n") );
	}

	RC_TRY( free_bc_array(&inc_force) );
	RC_TRY( free_default_bc_array(&inc_b_force) );

	return(NORMAL_RC);
}


/* LineSearch 処理付き Update Lagrange 法による非線形解析ソルバ */
static RC
update_NLsolve_iccg3 (NODE_ARRAY node, ELEMENT_ARRAY element,
                      MATERIAL_PROP_ARRAY material,
                      PHYSICAL_PROP_ARRAY physical,
                      BC_ARRAY rest, int iter_num, BC_ARRAY force,
                      DEFAULT_BC_ARRAY b_force, DISP_ARRAY *disp,
                      STRESS_ARRAY *stress, STRAIN_ARRAY *strain,
                      int sol_ss_flag)
{
	int ii1, ii2, ii3;
	static int dim, total_dof;
	static double *r_vect;
	const int ls_num = 5;
	const double eps_r = 1.0E-3;
	const double eps_d = 1.0E-3;
	double ls_fact = 1.0;
	double init_r_norm = 0.0;
	double init_d_norm = 0.0;
	double r_norm = 0.0;
	double d_norm = 0.0;
	NONZERO_MATRIX3 matrix;
	ELEM_MATRIX elem_matrix;
	DISP_ARRAY delta_disp, add_disp;
	double *d_vect, *f_vect, *react_vect, *r_vect_old;

	static int init_flag = 0;
	static int nl_status_flag = 0;

	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );
	if(init_flag == 0){
		if( (dim = analysis_dim(element)) != 3) return(ARG_ERROR_RC);
		RC_TRY( fem_total_dof(node, element, &total_dof) );
		RC_TRY( fem_allocate_vector(dim, node, &r_vect) );

		init_flag = 1;
	}
	if(iter_num == 1){
		RC_TRY( allocate_disp_array(node.size, disp) );
		RC_TRY( fem_zero_vector(dim, node, r_vect) );
		nl_status_flag = 0;
	}
	RC_TRY( allocate_disp_array(node.size, &add_disp) );
	RC_TRY( fem_allocate_vector(dim, node, &f_vect) );
	RC_TRY( fem_allocate_vector(dim, node, &react_vect) );
	RC_TRY( fem_allocate_vector(dim, node, &r_vect_old) );

	RC_TRY( add_nodal_force(dim, node, force, f_vect) );
	RC_TRY( add_body_force(dim, node, element, material, physical,
	                       b_force, f_vect) );

	for(ii1=0; ii1<total_dof; ii1++){
		r_vect[ii1] += f_vect[ii1];
	}
	init_r_norm = inner_product(total_dof, r_vect, r_vect);
	RC_TRY( log_printf(4, "              P               U"
	                      "LS_FACTOR\n") );

	for(ii1=0; ii1<30; ii1++){
		RC_TRY( allocate_nonzero3(node, element, &matrix) );
		for(ii2=0; ii2<element.size; ii2++){
			if(element.array[ii2].label < 0) continue;
			if( nl_status_flag == 0 ){
				RC_TRY( make_elem_matrix(element.array[ii2], node,
				                         material, physical, &elem_matrix) );
			}else{
				RC_TRY( make_tangent_matrix(element.array[ii2],
				                            node, material, physical,
				                            *disp, &elem_matrix) );
			}
			RC_TRY( impose_nonzero3(matrix, elem_matrix, 1.0) );
			RC_TRY( free_elem_matrix(&elem_matrix) );
		}
		RC_TRY( modify_nonzero3(node, rest, matrix, r_vect) );
		/*RC_TRY( iccg_nonzero3(matrix, r_vect, &d_vect, NULL) );*/
		RC_TRY( scale_iccg_nonzero3(matrix, r_vect, &d_vect, NULL) );
		RC_TRY( vect2disp(3, node, d_vect, &delta_disp) );

		/* LineSearch 処理 */
		ls_fact = 1.0;
		for(ii2=0; ; ii2++){
			if(disp->size > 0) RC_TRY( copy_disp_array(*disp, &add_disp) );
			RC_TRY( add_disp_array(ls_fact, delta_disp, &add_disp) );
			RC_TRY( fem_zero_vector(dim, node, react_vect) );
			RC_TRY( add_NLreact_force(element, node, material, physical,
			                          add_disp, react_vect) );
			for(ii3=0; ii3<total_dof; ii3++){
				r_vect[ii3] = f_vect[ii3]*iter_num - react_vect[ii3];
			}

			if(ii1 == 0 || ii2 == ls_num) break; /* 反復1回目or指定回数終了 */
			if( cal_linesearch_factor(total_dof, ls_num, d_vect, r_vect,
			                          r_vect_old, &ls_fact) ) break;
		}
		RC_TRY( copy_disp_array(add_disp, disp) );
		for(ii2=0; ii2<total_dof; ii2++){
			r_vect_old[ii2] = r_vect[ii2];
		}

		for(ii2=0; ii2<rest.size; ii2++){
			int index;
			if(rest.array[ii2].node < 0) continue;
			RC_NEG_CHK( index = node_renum_index(node, rest.array[ii2].node) );
			if(rest.array[ii2].v_type.x == BC_FIX) r_vect[3*index]   = 0.0;
			if(rest.array[ii2].v_type.y == BC_FIX) r_vect[3*index+1] = 0.0;
			if(rest.array[ii2].v_type.z == BC_FIX) r_vect[3*index+2] = 0.0;
		}

		/* ノルム計算&初期ノルム記憶 */
		r_norm = inner_product(total_dof, r_vect, r_vect);
		d_norm = inner_product(total_dof, d_vect, d_vect);
		if(ii1 == 0){
			init_d_norm = d_norm;
			RC_TRY( log_printf(4, "init: %15.7e %15.7e\n",
			                   init_r_norm, init_d_norm) );
			nl_status_flag = 1;
		}
		RC_TRY( log_printf(4, "%4d: %15.7e %15.7e  %10.5e\n",
		                      ii1, r_norm, d_norm,ls_fact) );
		RC_TRY( free_nonzero_matrix3(&matrix) );
		RC_TRY( fem_free_vector(&d_vect) );
		RC_TRY( free_disp_array(&delta_disp) );

		/* 収束判定 */
		if( !( r_norm>(init_r_norm*eps_r) || d_norm>(init_d_norm*eps_d) ) ){
			break;
		}
	}

	RC_TRY( free_disp_array(&add_disp) );
	RC_TRY( fem_free_vector(&react_vect) );
	RC_TRY( fem_free_vector(&f_vect) );
	RC_TRY( fem_free_vector(&r_vect_old) );

	/* Green ひずみ、応力(第2Piola-Krchhoff応力)の計算 */
	RC_TRY( stress_Green_strain(sol_ss_flag, node, element, *disp,
	                            material, physical, stress, strain) );

	return(NORMAL_RC);
}


/* 要素の大変位による反力 */
RC
add_NLreact_force (ELEMENT_ARRAY element, NODE_ARRAY node,
                   MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
                   DISP_ARRAY disp, double *f_vect)
{
	int dim, ssnum;
	int ii1, ii2, ii3;
	double thickness;
	double react_v[3*MAX_NODE];
	double stress_v[6];
	static int init_flag = 0;
	static double **D_matrix = NULL;
	static double **node_location = NULL;
	static double **B_matrix = NULL;
	static double **g_points = NULL;
	double g_weight[MAX_GAUSS_POINT];
	STRESS_STRAIN strain;
	int g_num_point;
	double det_J;

	if(init_flag == 0){
		RC_TRY( allocate2D(6, 6, &D_matrix) );
		RC_TRY( allocate2D(MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(6, 3*MAX_NODE, &B_matrix) );
		RC_TRY( allocate2D(MAX_GAUSS_POINT, 4, &g_points) );

		init_flag = 1;
	}

	if((dim = analysis_dim(element)) <= 1){
		fprintf(stderr, "analysis_dim() < ");
		return(ARG_ERROR_RC);
	}

	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		if( (ssnum = element_ssnum(element.array[ii1].type)) <= 0){
			fprintf(stderr, "element_ssnum < ");
			return(ARG_ERROR_RC);
		}

		RC_TRY( Gauss_const_default(element.array[ii1].type, 
		                            g_points, g_weight, &g_num_point) );
		RC_TRY( disp_add_node_location_matrix(&(element.array[ii1]), node,
		                                      disp, node_location) );

		RC_TRY( get_D_matrix(element.array[ii1],
		                     material, physical, D_matrix) );
		RC_TRY( get_thickness(element.array[ii1], physical, &thickness) );

		for(ii2=0; ii2<g_num_point; ii2++){
			RC_TRY( make_B_matrix(element.array[ii1], g_points[ii2],
			                      node_location, B_matrix, &det_J) );
			RC_TRY( local_Green_strain(element.array[ii1], g_points[ii2],
			                           node, disp, &strain) );
			
			if(dim == 2){
				for(ii3=0; ii3<ssnum; ii3++){
					stress_v[ii3] = D_matrix[ii3][0]*strain.v.x
					              + D_matrix[ii3][1]*strain.v.y
					              + D_matrix[ii3][2]*strain.v.xy;
				}
			}else{  /* dim == 3 */
				for(ii3=0; ii3<ssnum; ii3++){
					stress_v[ii3] = D_matrix[ii3][0]*strain.v.x
					              + D_matrix[ii3][1]*strain.v.y
					              + D_matrix[ii3][2]*strain.v.z
					              + D_matrix[ii3][3]*strain.v.yz
					              + D_matrix[ii3][4]*strain.v.zx
					              + D_matrix[ii3][5]*strain.v.xy;
				}
			}
			mul_trans_matrix_vect(ssnum, dim*element.array[ii1].node_num,
			                      B_matrix, stress_v, react_v);
			for(ii3=0; ii3<element.array[ii1].node_num; ii3++){
				int index;
				index = node_renum_index(node, element.array[ii1].node[ii3]);
				if(index < 0) return(SEEK_ERROR_RC);
				f_vect[index*dim]   += g_weight[ii2] * det_J
				                     * thickness*react_v[ii3*dim];
				f_vect[index*dim+1] += g_weight[ii2] * det_J
				                     * thickness*react_v[ii3*dim+1];
				if(dim == 3){
					f_vect[index*dim+2] += g_weight[ii2] * det_J
					                     * thickness*react_v[ii3*dim+2];
				}
			}
		}
	}

	return(NORMAL_RC);
}


/* LineSearch 法による factor 計算              */
/* エネルギー誤差値が収束していた場合、1 を返す */
static int
cal_linesearch_factor (int total_dof, int ls_num, double d_vect[],
                       double r_vect[], double r_vect_old[], double *ls_fact)
{
	double e_res, e_res0, ls_fact_tmp;
	const double e_allow = 0.5;
	static double e_res_old, ls_fact_old;

	/* エネルギー誤差計算 */
	e_res  = inner_product(total_dof, d_vect, r_vect);
	e_res0 = inner_product(total_dof, d_vect, r_vect_old);
	if( (e_res*e_res0 > 0.0) && (fabs(e_res0) < fabs(e_res)) ) e_res = 0.0;
	e_res /= e_res0;

	if(e_res > e_allow){ /* 2倍化スキーム処理 */
		ls_fact_tmp = 2.0*(*ls_fact);
	}else if(e_res < -e_allow){
		if(nearly_eq(*ls_fact, 1.0)){ /* 初回: 2次内挿処理 */
			ls_fact_tmp = (1.0 - sqrt(1.0 - 4.0*e_res))/(2.0*e_res);
		}else{                        /* 初回以降: 1次内挿処理 */
			if(e_res*e_res_old > 0.0) e_res_old /= 2.0;
			ls_fact_tmp = *ls_fact - e_res*(*ls_fact - ls_fact_old)
			                              /(e_res - e_res_old);
		}
	}else{               /* 収束 */
		return(1);
	}

	/* 求めた factor が既定値より小さい -> 2分化処理を行う */
	if(ls_fact_tmp < 1.0/(10.0*ls_num)){
		ls_fact_tmp = (ls_fact_tmp+(*ls_fact))/2.0;
	}

	e_res_old = e_res;
	ls_fact_old = *ls_fact;
	*ls_fact = ls_fact_tmp;

	return(0);
}


/* 応力(第2Piola-Kirchhoff応力)、Green ひずみの計算              */
/* stress, strain が NULL なら代入しない */
RC
stress_Green_strain (int sol_ss_flag, NODE_ARRAY node,
                     ELEMENT_ARRAY element, DISP_ARRAY disp,
                     MATERIAL_PROP_ARRAY material,
                     PHYSICAL_PROP_ARRAY physical,
                     STRESS_ARRAY *stress, STRAIN_ARRAY *strain)
{
	int ii1, ii2;
	int node_index;
	int *count = NULL;
	static int init_flag = 0;
	static double **local_points;
	STRESS_STRAIN local_stress, local_strain;
	VECT3DR *stress_sum = NULL;
	VECT3DR *strain_sum = NULL;
	VECT3DR init_v = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

	if( ((stress == NULL)&&(strain == NULL)) || (sol_ss_flag == 0) ){
		return(NORMAL_RC);
	}
	if(init_flag == 0){
		RC_TRY( allocate2D(MAX_NODE, 4, &local_points) );
		init_flag = 1;
	}

	if(sol_ss_flag & SOL_SS_NODE_AV){
		stress_sum = (VECT3DR *)mm_alloc((node.size)*sizeof(VECT3DR));
		strain_sum = (VECT3DR *)mm_alloc((node.size)*sizeof(VECT3DR));
		count = (int *)mm_alloc((node.size)*sizeof(int));
		if((stress_sum == NULL)||(strain_sum == NULL)||(count == NULL))
			return(ALLOC_ERROR_RC);
		for(ii1=0; ii1<node.size; ii1++){
			count[ii1] = 0;
			stress_sum[ii1] = strain_sum[ii1] = init_v;
		}
	}

	if(stress != NULL){
		RC_TRY( allocate_stress_array(0, stress) );
		stress->source = FEM_K_SOL;
	}
	if(strain != NULL){
		RC_TRY( allocate_strain_array(0, strain) );
		strain->source = FEM_K_SOL;
	}

	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		/* 要素中心 */
		if(sol_ss_flag & SOL_SS_CENTER){
			RC_TRY( set_center_point(element.array[ii1].type,
			                         local_points[0]) );
			RC_TRY( local_Green_strain(element.array[ii1], local_points[0],
			                           node, disp, &local_strain) );
			RC_TRY( local_strain2local_stress(element.array[ii1], material,
			                                  physical, local_strain,
			                                  &local_stress) );
			if(stress != NULL){
				RC_TRY( realloc_stress_array(stress) );
				stress->array[stress->size] = local_stress;
				stress->array[stress->size].node = -1;  /* center */
				(stress->size)++;
			}
			if(strain != NULL){
				RC_TRY( realloc_strain_array(strain) );
				strain->array[strain->size] = local_strain;
				strain->array[strain->size].node = -1;  /* center */
				(strain->size)++;
			}
		}

		/* 各節点 */
		if( (sol_ss_flag & SOL_SS_NODE) ||(sol_ss_flag & SOL_SS_NODE_AV) ){
			RC_TRY( set_local_node_points(element.array[ii1].type,
			                              local_points) );

			for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
				RC_TRY( local_Green_strain(element.array[ii1],
				                           local_points[ii2], node, disp,
				                           &local_strain) );
				RC_TRY( local_strain2local_stress(element.array[ii1],
				                                  material, physical,
				                                  local_strain,
				                                  &local_stress) );

				if(sol_ss_flag & SOL_SS_NODE_AV){
					node_index = search_node_label(node,
					                   element.array[ii1].node[ii2]);
					if(node_index < 0) return(SEEK_ERROR_RC);
					count[node_index]++;
					stress_sum[node_index]
					    = add_vect3dr(stress_sum[node_index], local_stress.v);
					strain_sum[node_index]
					    = add_vect3dr(strain_sum[node_index], local_strain.v);
				}

				if(stress != NULL){
					RC_TRY( realloc_stress_array(stress) );
					stress->array[stress->size] = local_stress;
					stress->array[stress->size].node
					      = element.array[ii1].node[ii2];
					(stress->size)++;
				}

				if(strain != NULL){
					RC_TRY( realloc_strain_array(strain) );
					strain->array[strain->size] = local_strain;
					strain->array[strain->size].node
					      = element.array[ii1].node[ii2];
					(strain->size)++;
				}
			}
		}
	}

	if(stress != NULL) RC_TRY( clean_stress_array(stress) );
	if(strain != NULL) RC_TRY( clean_strain_array(strain) );

	if(sol_ss_flag & SOL_SS_NODE_AV){
		if(stress != NULL){
			for(ii1=0; ii1<stress->size; ii1++){
				if(stress->array[ii1].element < 0) continue;
				node_index = search_node_label(node, stress->array[ii1].node);
				if(node_index < 0) return(SEEK_ERROR_RC);
				if(count[node_index] <= 0) return(UNKNOWN_ERROR_RC);
				stress->array[ii1].v
				 = mul_scalar_vect3dr(1.0/(double)(count[node_index]),
				                      stress_sum[node_index]);
			}
		}
		if(strain != NULL){
			for(ii1=0; ii1<strain->size; ii1++){
				if(strain->array[ii1].element < 0) continue;
				node_index = search_node_label(node, strain->array[ii1].node);
				if(node_index < 0) return(SEEK_ERROR_RC);
				if(count[node_index] <= 0) return(UNKNOWN_ERROR_RC);
				strain->array[ii1].v
				 = mul_scalar_vect3dr(1.0/(double)(count[node_index]),
				                      strain_sum[node_index]);
			}
		}
		RC_TRY( mm_free(stress_sum) );
		RC_TRY( mm_free(strain_sum) );
		RC_TRY( mm_free(count) );
	}

	return(NORMAL_RC);
}


/* 要素の局所 Green ひずみ */
RC
local_Green_strain (ELEMENT element, const double *local_point,
                    NODE_ARRAY node, DISP_ARRAY disp, STRESS_STRAIN *strain)
{
	int dim, ii1, ii2;
	static int init_flag = 0;
	static double **dN = NULL;
	static double **dNxyz = NULL;
	static double **Jacobi = NULL;
	static double **inverse_J = NULL;
	static double **node_location = NULL;
	static double **new_node_location = NULL;
	static double **F_tensor = NULL;
	static double **Ft = NULL;
	static double **strain_mat = NULL;

	if(init_flag == 0){
		RC_TRY( allocate2D(3, MAX_NODE, &dN) );
		RC_TRY( allocate2D(3, MAX_NODE, &dNxyz) );
		RC_TRY( allocate2D(3, 3, &Jacobi) );
		RC_TRY( allocate2D(3, 3, &inverse_J) );
		RC_TRY( allocate2D(MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(MAX_NODE, 3, &new_node_location) );
		RC_TRY( allocate2D(3, 3, &F_tensor) );
		RC_TRY( allocate2D(3, 3, &Ft) );
		RC_TRY( allocate2D(3, 3, &strain_mat) );

		init_flag = 1;
	}

	RC_NEG_ZERO_CHK( dim = element_dim(element.type) );
	init_stress_strain(strain);
	strain->element = element.label;

	/* 初期配置での節点座標マトリクス */
	RC_TRY( node_location_matrix(&element, node, node_location) );

	/* 現配置での節点座標マトリクス */
	RC_TRY( disp_add_node_location_matrix(&element, node, disp,
	                                      new_node_location) );

	RC_TRY( set_dN_matrix(element.type, local_point, dN) );
	mul_matrix(dim, element.node_num, dim, dN, node_location, Jacobi);
	RC_TRY( inverse_matrix(dim, Jacobi, inverse_J) );
	mul_matrix(dim, dim, element.node_num, inverse_J, dN, dNxyz);

	/* 変形勾配テンソル F_tensor */
	mul_matrix(dim, element.node_num, dim, dNxyz, new_node_location, F_tensor);
	transpose_matrix(dim, dim, F_tensor, Ft);

	mul_matrix(dim, dim, dim, Ft, F_tensor, strain_mat);

	/* Green-Lagrange ひずみテンソル */
	for(ii1=0; ii1<dim; ii1++){
		strain_mat[ii1][ii1] -= 1.0;
		for(ii2=0; ii2<dim; ii2++){
			strain_mat[ii1][ii2] *= 0.5;
		}
	}

	/* テンソルひずみ -> 工学ひずみ */
	if(dim == 2){
		strain->v.x  = strain_mat[0][0];
		strain->v.y  = strain_mat[1][1];
		strain->v.xy = strain_mat[0][1]*2.0;
	}else{  /* dim == 3 */
		strain->v.x  = strain_mat[0][0];
		strain->v.y  = strain_mat[1][1];
		strain->v.z  = strain_mat[2][2];
		strain->v.yz = strain_mat[1][2]*2.0;
		strain->v.zx = strain_mat[2][0]*2.0;
		strain->v.xy = strain_mat[0][1]*2.0;
	}

	return(NORMAL_RC);
}


/* 大変形の A,G マトリックス(Zienkiewicz 参照) */
RC
make_AG_matrix (ELEMENT element, const double *g_point,
                double **node_location, DISP_ARRAY disp,
                double **A_matrix, double **G_matrix)
{
	static int allocate_flag = 0;
	static double **dN = NULL;
	static double **dNxyz = NULL;
	static double **Jacobi = NULL;
	static double **inverse_J = NULL;
	static double **disp_matrix = NULL;
	static double **duvw = NULL;
	int dim;

	if(allocate_flag == 0){
		RC_TRY( allocate2D(3, MAX_NODE, &dN) );
		RC_TRY( allocate2D(3, MAX_NODE, &dNxyz) );
		RC_TRY( allocate2D(3, 3, &Jacobi) );
		RC_TRY( allocate2D(3, 3, &inverse_J) );
		RC_TRY( allocate2D(MAX_NODE, 6, &disp_matrix) );
		RC_TRY( allocate2D(3, 3, &duvw) );

		allocate_flag = 1;
	}

	if((dim = element_dim(element.type)) <= 0) return(ARG_ERROR_RC);

	/* g_point -> dN                */
	/* dN * node_location => Jacobi */
	/* Jacobi^-1 => inverse_J       */
	/* inverse_J * dN => dNxyz      */
	RC_TRY( set_dN_matrix(element.type, g_point, dN) );
	mul_matrix(dim, element.node_num, dim, dN, node_location, Jacobi);
	RC_TRY( inverse_matrix(dim, Jacobi, inverse_J) );
	mul_matrix(dim, dim, element.node_num, inverse_J, dN, dNxyz);

	RC_TRY( set_G_matrix(dim, element.node_num, dNxyz, G_matrix) );

	RC_TRY( fill_disp_matrix(element, disp, disp_matrix) );
	mul_matrix(dim, element.node_num, dim, dNxyz, disp_matrix, duvw);

	RC_TRY( set_A_matrix(dim, duvw, A_matrix) );

	return(NORMAL_RC);
}

/* 大変形のGマトリックスのみ作成(Zienkiewicz 参照) */
RC
make_G_matrix (ELEMENT element, const double *g_point,
               double **node_location, double **G_matrix)
{
	static int allocate_flag = 0;
	static double **dN = NULL;
	static double **dNxyz = NULL;
	static double **Jacobi = NULL;
	static double **inverse_J = NULL;
	int dim;

	if(allocate_flag == 0){
		RC_TRY( allocate2D(3, MAX_NODE, &dN) );
		RC_TRY( allocate2D(3, MAX_NODE, &dNxyz) );
		RC_TRY( allocate2D(3, 3, &Jacobi) );
		RC_TRY( allocate2D(3, 3, &inverse_J) );

		allocate_flag = 1;
	}

	if((dim = element_dim(element.type)) <= 0) return(ARG_ERROR_RC);

	/* g_point -> dN                */
	/* dN * node_location => Jacobi */
	/* Jacobi^-1 => inverse_J       */
	/* inverse_J * dN => dNxyz      */
	RC_TRY( set_dN_matrix(element.type, g_point, dN) );
	mul_matrix(dim, element.node_num, dim, dN, node_location, Jacobi);
	RC_TRY( inverse_matrix(dim, Jacobi, inverse_J) );
	mul_matrix(dim, dim, element.node_num, inverse_J, dN, dNxyz);

	RC_TRY( set_G_matrix(dim, element.node_num, dNxyz, G_matrix) );

	return(NORMAL_RC);
}

static RC
set_G_matrix (int dim, int node_num, double **dNxyz, double **G_matrix)
{
	int ii1;

	switch(dim){
	case 2:
		for(ii1=0;ii1<node_num;ii1++){
			G_matrix[0][dim*ii1]   = dNxyz[0][ii1];
			G_matrix[0][dim*ii1+1] = 0.0;
			G_matrix[1][dim*ii1]   = 0.0;
			G_matrix[1][dim*ii1+1] = dNxyz[0][ii1];

			G_matrix[2][dim*ii1]   = dNxyz[1][ii1];
			G_matrix[2][dim*ii1+1] = 0.0;
			G_matrix[3][dim*ii1]   = 0.0;
			G_matrix[3][dim*ii1+1] = dNxyz[1][ii1];
		}

		break;
	case 3:
		for(ii1=0;ii1<node_num;ii1++){
			G_matrix[0][dim*ii1]   = dNxyz[0][ii1];
			G_matrix[0][dim*ii1+1] = 0.0;
			G_matrix[0][dim*ii1+2] = 0.0;
			G_matrix[1][dim*ii1]   = 0.0;
			G_matrix[1][dim*ii1+1] = dNxyz[0][ii1];
			G_matrix[1][dim*ii1+2] = 0.0;
			G_matrix[2][dim*ii1]   = 0.0;
			G_matrix[2][dim*ii1+1] = 0.0;
			G_matrix[2][dim*ii1+2] = dNxyz[0][ii1];

			G_matrix[3][dim*ii1]   = dNxyz[1][ii1];
			G_matrix[3][dim*ii1+1] = 0.0;
			G_matrix[3][dim*ii1+2] = 0.0;
			G_matrix[4][dim*ii1]   = 0.0;
			G_matrix[4][dim*ii1+1] = dNxyz[1][ii1];
			G_matrix[4][dim*ii1+2] = 0.0;
			G_matrix[5][dim*ii1]   = 0.0;
			G_matrix[5][dim*ii1+1] = 0.0;
			G_matrix[5][dim*ii1+2] = dNxyz[1][ii1];

			G_matrix[6][dim*ii1]   = dNxyz[2][ii1];
			G_matrix[6][dim*ii1+1] = 0.0;
			G_matrix[6][dim*ii1+2] = 0.0;
			G_matrix[7][dim*ii1]   = 0.0;
			G_matrix[7][dim*ii1+1] = dNxyz[2][ii1];
			G_matrix[7][dim*ii1+2] = 0.0;
			G_matrix[8][dim*ii1]   = 0.0;
			G_matrix[8][dim*ii1+1] = 0.0;
			G_matrix[8][dim*ii1+2] = dNxyz[2][ii1];
		}
		
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
set_A_matrix (int dim, double **duvw, double **A_matrix)
{
	switch(dim){
	case 2:
		A_matrix[0][0] = duvw[0][0];
		A_matrix[0][1] = duvw[0][1];
		A_matrix[0][2] = 0.0;
		A_matrix[0][3] = 0.0;

		A_matrix[1][0] = 0.0;
		A_matrix[1][1] = 0.0;
		A_matrix[1][2] = duvw[1][0];
		A_matrix[1][3] = duvw[1][1];

		A_matrix[2][0] = duvw[1][0];
		A_matrix[2][1] = duvw[1][0];
		A_matrix[2][2] = duvw[0][0];
		A_matrix[2][3] = duvw[0][1];

		break;
	case 3:
		A_matrix[0][0] = duvw[0][0];
		A_matrix[0][1] = duvw[0][1];
		A_matrix[0][2] = duvw[0][2];
		A_matrix[0][3] = 0.0;
		A_matrix[0][4] = 0.0;
		A_matrix[0][5] = 0.0;
		A_matrix[0][6] = 0.0;
		A_matrix[0][7] = 0.0;
		A_matrix[0][8] = 0.0;

		A_matrix[1][0] = 0.0;
		A_matrix[1][1] = 0.0;
		A_matrix[1][2] = 0.0;
		A_matrix[1][3] = duvw[1][0];
		A_matrix[1][4] = duvw[1][1];
		A_matrix[1][5] = duvw[1][2];
		A_matrix[1][6] = 0.0;
		A_matrix[1][7] = 0.0;
		A_matrix[1][8] = 0.0;

		A_matrix[2][0] = 0.0;
		A_matrix[2][1] = 0.0;
		A_matrix[2][2] = 0.0;
		A_matrix[2][3] = 0.0;
		A_matrix[2][4] = 0.0;
		A_matrix[2][5] = 0.0;
		A_matrix[2][6] = duvw[2][0];
		A_matrix[2][7] = duvw[2][1];
		A_matrix[2][8] = duvw[2][2];

		A_matrix[3][0] = 0.0;
		A_matrix[3][1] = 0.0;
		A_matrix[3][2] = 0.0;
		A_matrix[3][3] = duvw[2][0];
		A_matrix[3][4] = duvw[2][1];
		A_matrix[3][5] = duvw[2][2];
		A_matrix[3][6] = duvw[1][0];
		A_matrix[3][7] = duvw[1][1];
		A_matrix[3][8] = duvw[1][2];

		A_matrix[4][0] = duvw[2][0];
		A_matrix[4][1] = duvw[2][1];
		A_matrix[4][2] = duvw[2][2];
		A_matrix[4][3] = 0.0;
		A_matrix[4][4] = 0.0;
		A_matrix[4][5] = 0.0;
		A_matrix[4][6] = duvw[0][0];
		A_matrix[4][7] = duvw[0][1];
		A_matrix[4][8] = duvw[0][2];

		A_matrix[5][0] = duvw[1][0];
		A_matrix[5][1] = duvw[1][1];
		A_matrix[5][2] = duvw[1][2];
		A_matrix[5][3] = duvw[0][0];
		A_matrix[5][4] = duvw[0][1];
		A_matrix[5][5] = duvw[0][2];
		A_matrix[5][6] = 0.0;
		A_matrix[5][7] = 0.0;
		A_matrix[5][8] = 0.0;

		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


/* 接線剛性マトリックスの作成 */
RC
make_tangent_matrix (ELEMENT element, NODE_ARRAY node,
                     MATERIAL_PROP_ARRAY material,
                     PHYSICAL_PROP_ARRAY physical,
                     DISP_ARRAY disp, ELEM_MATRIX *elem_matrix)
{
	int dim, ssnum;
	int ii1, ii2, ii3;
	int g_num_point;
	double thickness;
	double det_J;
	double stress_v[6];
	double g_weight[MAX_GAUSS_POINT];
	static double **D_matrix = NULL;
	static double **node_location = NULL;
	static double **B_matrix = NULL;
	static double **G_matrix = NULL;
	static double **M_matrix = NULL;
	static double **BtDB = NULL;
	static double **GtMG = NULL;
	static double **g_points = NULL;
	static int init_flag = 0;
	STRESS_STRAIN strain;

	if(init_flag == 0){
		RC_TRY( allocate2D(6, 6, &D_matrix) );
		RC_TRY( allocate2D(MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(6, 3*MAX_NODE, &B_matrix) );
		RC_TRY( allocate2D(3*3, 3*MAX_NODE, &G_matrix) );
		RC_TRY( allocate2D(3*3, 3*3, &M_matrix) );
		RC_TRY( allocate2D(3*MAX_NODE, 3*MAX_NODE, &BtDB) );
		RC_TRY( allocate2D(3*MAX_NODE, 3*MAX_NODE, &GtMG) );
		RC_TRY( allocate2D(MAX_GAUSS_POINT, 4, &g_points) );

		init_flag = 1;
	}

	if(element.label < 0) return(ARG_ERROR_RC);

	if((elem_matrix->dim = element_dim(element.type)) <= 0){
		return(ARG_ERROR_RC);
	}
	dim = elem_matrix->dim;
	if( (ssnum = element_ssnum(element.type)) <= 0) return(ARG_ERROR_RC);
	elem_matrix->element = element.label;
	elem_matrix->size = (element.node_num) * (elem_matrix->dim);

	/* elem_matrix->array[][] の確保 */
	RC_TRY( allocate2D(elem_matrix->size,
	                   elem_matrix->size, &(elem_matrix->matrix)) );

	elem_matrix->index = mm_alloc(sizeof(int)*(elem_matrix->size));
	if(elem_matrix->index == NULL) return(ALLOC_ERROR_RC);

	/* elem_matirx->index[] の確保と代入 */
	for(ii1=0; ii1<element.node_num; ii1++){
		int index = node_renum_index(node, element.node[ii1]);
		if(index < 0) return(SEEK_ERROR_RC);
		for(ii2=0; ii2<(elem_matrix->dim); ii2++){
			elem_matrix->index[ii1*(elem_matrix->dim)+ii2]
			             = index*(elem_matrix->dim) + ii2;
		}
	}

	RC_TRY( Gauss_const_default(element.type,
	                            g_points, g_weight, &g_num_point) );
	RC_TRY( disp_add_node_location_matrix(&element, node, disp,
	                                      node_location) );
	RC_TRY( get_D_matrix(element, material, physical, D_matrix) );
	RC_TRY( get_thickness(element, physical, &thickness) );

	for(ii1=0; ii1<g_num_point; ii1++){
		RC_TRY( make_B_matrix(element, g_points[ii1],
		                      node_location, B_matrix, &det_J) );
		RC_TRY( make_G_matrix(element, g_points[ii1], node_location,
		                      G_matrix) );
		RC_TRY( local_Green_strain(element, g_points[ii1], node, disp,
		                           &strain) );

		if(dim == 2){
			for(ii3=0; ii3<ssnum; ii3++){
				stress_v[ii3] = D_matrix[ii3][0]*strain.v.x
				              + D_matrix[ii3][1]*strain.v.y
				              + D_matrix[ii3][2]*strain.v.xy;
			}
		}else{  /* dim == 3 */
			for(ii3=0; ii3<ssnum; ii3++){
				stress_v[ii3] = D_matrix[ii3][0]*strain.v.x
				              + D_matrix[ii3][1]*strain.v.y
				              + D_matrix[ii3][2]*strain.v.z
				              + D_matrix[ii3][3]*strain.v.yz
				              + D_matrix[ii3][4]*strain.v.zx
				              + D_matrix[ii3][5]*strain.v.xy;
			}
		}
		RC_TRY( set_M_matrix(dim, stress_v, M_matrix) );

		mul_matrix_AtBA(ssnum, dim*element.node_num, B_matrix, D_matrix, BtDB);
		mul_matrix_AtBA(dim*dim, dim*element.node_num,
		                G_matrix, M_matrix, GtMG);

		/* (BtDB+GtMG) * detJ * weight * thickness */
		for(ii2=0; ii2<(elem_matrix->size); ii2++){
			for(ii3=0; ii3<(elem_matrix->size); ii3++){
				elem_matrix->matrix[ii2][ii3]
				 += (BtDB[ii2][ii3]+GtMG[ii2][ii3])
				  * det_J * g_weight[ii1] * thickness;
			}
		}
	}

	return(NORMAL_RC);
}


RC
set_M_matrix (int dim, double *stress, double **M_matrix)
{
	switch(dim){
	case 2:
		RC_TRY( set_subM_matrix(dim, stress[0], M_matrix, 0, 0) );
		RC_TRY( set_subM_matrix(dim, stress[2], M_matrix, 0, 1) );
		RC_TRY( set_subM_matrix(dim, stress[2], M_matrix, 1, 0) );
		RC_TRY( set_subM_matrix(dim, stress[1], M_matrix, 1, 1) );

		break;
	case 3:
		RC_TRY( set_subM_matrix(dim, stress[0], M_matrix, 0, 0) );
		RC_TRY( set_subM_matrix(dim, stress[5], M_matrix, 0, 1) );
		RC_TRY( set_subM_matrix(dim, stress[4], M_matrix, 0, 2) );
		RC_TRY( set_subM_matrix(dim, stress[5], M_matrix, 1, 0) );
		RC_TRY( set_subM_matrix(dim, stress[1], M_matrix, 1, 1) );
		RC_TRY( set_subM_matrix(dim, stress[3], M_matrix, 1, 2) );
		RC_TRY( set_subM_matrix(dim, stress[4], M_matrix, 2, 0) );
		RC_TRY( set_subM_matrix(dim, stress[3], M_matrix, 2, 1) );
		RC_TRY( set_subM_matrix(dim, stress[2], M_matrix, 2, 2) );

		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}

static RC
set_subM_matrix (int dim, double stress, double **M_matrix, int x, int y)
{
	int ii1, ii2;

	for(ii1=0; ii1<dim; ii1++){
		for(ii2=0; ii2<dim; ii2++){
			M_matrix[y*dim+ii1][x*dim+ii2] = 0.0;
		}
		M_matrix[y*dim+ii1][x*dim+ii1] = stress;
	}

	return(NORMAL_RC);
}


/* 要素 element を構成する節点座標値+変位のマトリックスを作成 */
static RC
disp_add_node_location_matrix (ELEMENT *element, NODE_ARRAY node,
                               DISP_ARRAY disp, double **matrix)
{
	int ii1;
	int index1, index2;

	for(ii1=0; ii1<element->node_num; ii1++){
		if(element->label < 0) continue;

		index1 = search_node_label(node, element->node[ii1]);
		index2 = search_disp_node(disp, element->node[ii1]);
		if((index1 < 0) || (index2 < 0)) return(SEEK_ERROR_RC);

		matrix[ii1][0] = node.array[index1].p.x + disp.array[index2].v.x;
		matrix[ii1][1] = node.array[index1].p.y + disp.array[index2].v.y;
		matrix[ii1][2] = node.array[index1].p.z + disp.array[index2].v.z;
	}

	return(NORMAL_RC);
}

