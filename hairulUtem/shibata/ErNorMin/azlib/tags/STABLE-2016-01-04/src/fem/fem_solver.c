/*********************************************************************
 * fem_solver.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Takaaki NAGATANI>
 *
 * Thanks to following contributors.
 *  <Tomoyuki TSUBATA> <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: fem_solver.c 1022 2012-04-30 05:00:38Z aoyama $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fem_solver.h"
#include "fem_struct.h"
#include "math_utl.h"
#include "sky_cholesky.h"
#include "rc.h"
#include "nst_component.h"
#include "nonzero_cg.h"
#include "memory_manager.h"

#define DEF_ANA_2D    ANA_PLANE_STRAIN


/* 要素上のガウスポイント毎の情報        */
/* ????[MAX_GAUSS_POINT]は要素中央の情報 */
typedef struct{
	int num_g_point;
	double **g_points;
	double weight[MAX_GAUSS_POINT+1];
	int B_matrix_row;
	int B_matrix_col;
	double **B_matrix[MAX_GAUSS_POINT+1];
	double det_J[MAX_GAUSS_POINT+1];
} G_PROP;

static RC allocate_g_prop(G_PROP *g_prop);
static RC fill_row_nodes(ELEM_LOOKUP elem_lookup,
                         NODE_ARRAY node, ELEMENT_ARRAY element,
                         int row_nodes[], int *row_size);
static RC element_body_force(int gdim, ELEMENT element, NODE_ARRAY node,
                             MATERIAL_PROP_ARRAY material,
                             PHYSICAL_PROP_ARRAY physical,
                             DEFAULT_BC_ARRAY def_b_force,
                             int *index_array, int *vect_size,
                             double *body_force);
static RC element_temp_load(int gdim, ELEMENT element, NODE_ARRAY node,
                            MATERIAL_PROP_ARRAY material,
                            PHYSICAL_PROP_ARRAY physical,
                            BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                            int *index_array, int *vect_size,
                            double *temp_load);
static RC element_strain_load(int gdim, ELEMENT element,
                              NODE_ARRAY node,
                              MATERIAL_PROP_ARRAY material,
                              PHYSICAL_PROP_ARRAY physical,
                              STRAIN_ARRAY init_strain,
                              int *index_array, int *vect_size,
                              double *strain_load);
static RC fill_index1(int dim, NODE_ARRAY node, ELEMENT_ARRAY element,
                      long *index1);
static RC RI_B_matrix(ELEMENT elem, G_PROP g_prop);
static RC set_g_prop(ELEMENT element, double **node_location,
                     G_PROP *g_prop);
static RC set_Gauss_const(ELEM_TYPE type, G_PROP *g_prop);
static RC set_B_matrix(int dim, int node_num,
                       double **dNxyz, double **B_matrix);
static RC mod_g_matrix(SKYLINE_MATRIX g_matrix);
static RC Gauss_const4mass_matrix(ELEM_TYPE type, double **points,
                                  double *weight, int *num_point);
static void get_D_matrix_therm_ortho(MATERIAL_PROP array, double **D_matrix);
static RC element_temp_load1(int gdim, ELEMENT element, NODE_ARRAY node,
                             MATERIAL_PROP_ARRAY material,
                             PHYSICAL_PROP_ARRAY physical,
                             BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                             double basis_temp, int *index_array,
                             int *vect_size, double *temp_load);
static RC temp_vector(double *f_vect, BC_ARRAY temp, NODE_ARRAY node);

RC
react_force (NODE_ARRAY node, ELEMENT_ARRAY element,
             MATERIAL_PROP_ARRAY material,
             PHYSICAL_PROP_ARRAY physical,
             DISP_ARRAY disp, REACT_ARRAY *react)
{
	double *r_vect;
	double *d_vect;
	ELEM_MATRIX elem_matrix;
	int dim;
	int ii1, ii2, ii3;

	dim = analysis_dim(element);
	RC_TRY( disp2vect(dim, node, disp, &d_vect) );
	RC_TRY( fem_allocate_vector(dim, node, &r_vect) );

	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		RC_TRY( make_elem_matrix(element.array[ii1], node,
		                         material, physical, &elem_matrix) );
		for(ii2=0; ii2<elem_matrix.size; ii2++){
			for(ii3=0; ii3<elem_matrix.size; ii3++){
				r_vect[elem_matrix.index[ii2]]
				      += elem_matrix.matrix[ii2][ii3]
				        *d_vect[elem_matrix.index[ii3]];
			}
		}
		RC_TRY( free_elem_matrix(&elem_matrix) );
	}

	RC_TRY( vect2react(dim, node, r_vect, react) );
	RC_TRY( fem_free_vector(&r_vect) );
	RC_TRY( fem_free_vector(&d_vect) );

	return(NORMAL_RC);
}


RC
fem_solve_iccg3 (NODE_ARRAY node, ELEMENT_ARRAY element,
                 MATERIAL_PROP_ARRAY material,
                 PHYSICAL_PROP_ARRAY physical,
                 BC_ARRAY rest, BC_ARRAY force,
                 DEFAULT_BC_ARRAY b_force,
                 BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                 DISP_ARRAY *disp, STRESS_ARRAY *stress,
                 STRAIN_ARRAY *strain, int sol_ss_flag)
{
	int ii1;
	NONZERO_MATRIX3 matrix;
	double *f_vect;
	double *d_vect;
	ELEM_MATRIX elem_matrix;
	int dim;

	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

	/* 全体剛性マトリックスの作成 */
	RC_TRY( allocate_nonzero3(node, element, &matrix) );
	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		RC_TRY( make_elem_matrix(element.array[ii1], node,
		                         material, physical, &elem_matrix) );
		RC_TRY( impose_nonzero3(matrix, elem_matrix, 1.0) );
		RC_TRY( free_elem_matrix(&elem_matrix) );
	}

	/* 荷重ベクトルの作成 */
	dim = analysis_dim(element);
	if(dim != 3) return(ARG_ERROR_RC);
	RC_TRY( fem_allocate_vector(dim, node, &f_vect) );
	RC_TRY( add_nodal_force(dim, node, force, f_vect) );
	RC_TRY( add_body_force(dim, node, element, material, physical,
	                       b_force, f_vect) );
	RC_TRY( add_temp_force(dim, node, element, material, physical,
	                       temp, def_temp, f_vect) );

	/* 変位の計算 */
	RC_TRY( modify_nonzero3(node, rest, matrix, f_vect) );
	/*RC_TRY( fem_iccg_nonzero3(matrix, f_vect, &d_vect, NULL) );*/
	RC_TRY( scale_iccg_nonzero3(matrix, f_vect, &d_vect, NULL) );
	RC_TRY( free_nonzero_matrix3(&matrix) );
	RC_TRY( fem_free_vector(&f_vect) );
	RC_TRY( vect2disp(3, node, d_vect, disp) );
	RC_TRY( fem_free_vector(&d_vect) );

	/* ひずみ、応力の計算 */
	RC_TRY( stress_strain(sol_ss_flag, node, element, *disp,
	                      material, physical, stress, strain) );

	return(NORMAL_RC);
}


RC
fem_solve (NODE_ARRAY node, ELEMENT_ARRAY element,
           MATERIAL_PROP_ARRAY material,
           PHYSICAL_PROP_ARRAY physical,
           BC_ARRAY rest, BC_ARRAY force,
           BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
           DISP_ARRAY *disp, STRESS_ARRAY *stress,
           STRAIN_ARRAY *strain, int sol_ss_flag)
{
	SKYLINE_MATRIX K_matrix;
	double *f_vect;

	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );
	RC_TRY( make_g_matrix(node, element, material, physical, &K_matrix) );
	RC_TRY( fem_allocate_vector(K_matrix.dim, node, &f_vect) );
	RC_TRY( add_nodal_force(K_matrix.dim, node, force, f_vect) );
	RC_TRY( add_temp_force(K_matrix.dim, node, element, material, physical,
	                       temp, def_temp, f_vect) );
	RC_TRY( modify_g_matrix(node, rest, K_matrix, f_vect) );
	RC_TRY( decomp_g_matrix(K_matrix) );
	RC_TRY( subst_g_matrix(node, K_matrix, f_vect, disp) );
	RC_TRY( free_g_matrix(&K_matrix) );
	RC_TRY( stress_strain(sol_ss_flag, node, element, *disp,
	                      material, physical, stress, strain) );
	RC_TRY( mm_free(f_vect) );

	return(NORMAL_RC);
}

RC
fem_solve_add_density(NODE_ARRAY node, ELEMENT_ARRAY element,
             MATERIAL_PROP_ARRAY material,
             PHYSICAL_PROP_ARRAY physical,
             BC_ARRAY rest, BC_ARRAY force,
             BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
             DISP_ARRAY *disp, STRESS_ARRAY *stress,
             STRAIN_ARRAY *strain, int sol_ss_flag
			 ,double *density)
{
	SKYLINE_MATRIX K_matrix;
	double *f_vect;

	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );
	RC_TRY( make_g_matrix_add_density(node, element, material, physical, &K_matrix,density) );
	RC_TRY( fem_allocate_vector(K_matrix.dim, node, &f_vect) );
	RC_TRY( add_nodal_force(K_matrix.dim, node, force, f_vect) );
	RC_TRY( add_temp_force(K_matrix.dim, node, element, material, physical,
	                       temp, def_temp, f_vect) );
	RC_TRY( modify_g_matrix(node, rest, K_matrix, f_vect) );
	RC_TRY( decomp_g_matrix(K_matrix) );
	RC_TRY( subst_g_matrix(node, K_matrix, f_vect, disp) );
	RC_TRY( free_g_matrix(&K_matrix) );
	RC_TRY( stress_strain(sol_ss_flag, node, element, *disp,
	                      material, physical, stress, strain) );
	RC_TRY( mm_free(f_vect) );


	return(NORMAL_RC);
}


RC
fem_solve_therm (NODE_ARRAY node, ELEMENT_ARRAY element,
                 MATERIAL_PROP_ARRAY material,
                 PHYSICAL_PROP_ARRAY physical,
                 BC_ARRAY bc_temp, BC_ARRAY *result_temp)
{
	SKYLINE_MATRIX g_matrix;
	double *s_vect;
	ELEM_MATRIX elem_matrix;
	int ii1;

	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

	/* g_matrix の確保と初期化 */
	RC_TRY( allocate_g_matrix(1, node, element, &g_matrix) );

	for(ii1=0;ii1<element.size;ii1++){
		if(element.array[ii1].label < 0) continue;

		RC_TRY( make_elem_matrix_therm(element.array[ii1], node, material,
		                                   physical, &elem_matrix) );

		/* elem_matrix -> g_matrix */
		RC_TRY( impose_elem_matrix(g_matrix, elem_matrix, 1.0) );
		RC_TRY( free_elem_matrix(&elem_matrix) );
	}

	RC_TRY( fem_allocate_vector(1, node, &s_vect) );
	RC_TRY( modify_g_matrix_therm(node, bc_temp, g_matrix, s_vect) );

	RC_TRY( decomp_g_matrix(g_matrix) );
	RC_TRY( subst_g_matrix_therm(node, g_matrix, s_vect, result_temp) );
	RC_TRY( free_g_matrix(&g_matrix) );
	RC_TRY( fem_free_vector(&s_vect) );

	return(NORMAL_RC);
}


/* 複数荷重条件の解析                                          */
/* forces, disps, stresses, strains は、num_force 個以上の配列 */
RC
fem_solve_multi_force (NODE_ARRAY node, ELEMENT_ARRAY element,
                       MATERIAL_PROP_ARRAY material,
                       PHYSICAL_PROP_ARRAY physical,
                       BC_ARRAY rest,
                       int num_force, const BC_ARRAY *forces,
                       BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                       DISP_ARRAY *disps, STRESS_ARRAY *stresses,
                       STRAIN_ARRAY *strains, int sol_ss_flag)
{
	SKYLINE_MATRIX K_matrix;
	double **f_vect;
	int ii1;

	if(num_force <= 0) return(ARG_ERROR_RC);

	f_vect = (double **)mm_alloc(num_force*sizeof(double *));
	if(f_vect == NULL) return(ALLOC_ERROR_RC);

	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

	RC_TRY( make_g_matrix(node, element, material, physical, &K_matrix) );
	for(ii1=0; ii1<num_force; ii1++){
		RC_TRY( fem_allocate_vector(K_matrix.dim, node, &(f_vect[ii1])) );
		RC_TRY( add_nodal_force(K_matrix.dim, node,
		                            forces[ii1], f_vect[ii1]) );
		RC_TRY( add_temp_force(K_matrix.dim, node, element,
		                           material, physical,
		                           temp, def_temp, f_vect[ii1]) );
	}
	RC_TRY( modify_g_matrix_n(node, rest, K_matrix, f_vect, num_force) );
	RC_TRY( decomp_g_matrix(K_matrix) );

	for(ii1=0; ii1<num_force; ii1++){
		RC_TRY( subst_g_matrix(node, K_matrix,
		                           f_vect[ii1], &(disps[ii1])) );
		RC_TRY( mm_free(f_vect[ii1]) );
	}
	
	RC_TRY( free_g_matrix(&K_matrix) );

	for(ii1=0; ii1<num_force; ii1++){
		RC_TRY( stress_strain(sol_ss_flag, node, element, disps[ii1],
		                          material, physical,
		                          &(stresses[ii1]), &(strains[ii1])) );
	}
	RC_TRY( mm_free(f_vect) );

	return(NORMAL_RC);
}


RC
free_g_matrix (SKYLINE_MATRIX *g_matrix)
{
	RC_TRY( s_cholesky_free(g_matrix->index2, g_matrix->array) );
	RC_TRY( mm_free(g_matrix->index1) );

	g_matrix->dim = 0;
	g_matrix->dof = 0;
	g_matrix->index1 = NULL;
	g_matrix->index2 = NULL;
	g_matrix->array = NULL;

	return(NORMAL_RC);
}


/* 応力、ひずみの計算                    */
/* stress, strain が NULL なら代入しない */
RC
stress_strain (int sol_ss_flag, NODE_ARRAY node,
               ELEMENT_ARRAY element, DISP_ARRAY disp,
               MATERIAL_PROP_ARRAY material,
               PHYSICAL_PROP_ARRAY physical,
               STRESS_ARRAY *stress, STRAIN_ARRAY *strain)
{
	int ii1, ii2;
	static int init_flag = 0;
	static double **local_points;
	STRESS_STRAIN local_stress, local_strain;
	VECT3DR *stress_sum = NULL;
	VECT3DR *strain_sum = NULL;
	int *count = NULL;
	VECT3DR init_v = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	int node_index;

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
		if( (stress_sum == NULL) ||(strain_sum == NULL) ||(count == NULL) )
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

			RC_TRY( local_points_stress_strain(element.array[ii1],
			                                   local_points[0],
			                                   node, disp, material, physical,
			                                   &local_stress, &local_strain) );
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
		if( (sol_ss_flag & SOL_SS_NODE)
		  ||(sol_ss_flag & SOL_SS_NODE_AV) ){
			RC_TRY( set_local_node_points(element.array[ii1].type,
			                              local_points) );

			for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
				RC_TRY( local_points_stress_strain(element.array[ii1],
				          local_points[ii2], node, disp, material, physical,
				          &local_stress, &local_strain) );

				if(sol_ss_flag & SOL_SS_NODE_AV){
					node_index = search_node_label(node,
					                   element.array[ii1].node[ii2]);
					RC_NEG_CHK( node_index );
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
				node_index = search_node_label(node,
				                                   stress->array[ii1].node);
				if(count[node_index] <= 0) return(UNKNOWN_ERROR_RC);
				stress->array[ii1].v
				 = mul_scalar_vect3dr(1.0/(double)(count[node_index]),
				                      stress_sum[node_index]);
			}
		}
		if(strain != NULL){
			for(ii1=0; ii1<strain->size; ii1++){
				if(strain->array[ii1].element < 0) continue;
				node_index = search_node_label(node,
				                                   strain->array[ii1].node);
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


/* 要素の局所ひずみ */
RC
local_points_strain (ELEMENT element, const double *local_point,
                     NODE_ARRAY node, DISP_ARRAY disp, STRESS_STRAIN *strain)
{
	double disp_v[3*MAX_NODE];
	static double **node_location = NULL;
	static double **B_matrix = NULL;
	static int init_flag = 0;

	if(init_flag == 0){
		RC_TRY( allocate2D(MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(6, 3*MAX_NODE, &B_matrix) );

		init_flag = 1;
	}

	if(element.label < 0) return(ARG_ERROR_RC);
	if(element_dim(element.type) <= 1) return(ARG_ERROR_RC);
	if(element_ssnum(element.type) <= 0) return(ARG_ERROR_RC);

	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( fill_disp_vector(element, disp, disp_v) );

	RC_TRY( local_points_strain_sub(element, local_point, node_location,
	                                disp_v, B_matrix, strain) );

	return(NORMAL_RC);
}


RC
local_points_strain_sub (ELEMENT element, const double local_point[],
                         double **node_location, const double disp_v[],
                         double **ws_mat_6_3N, STRESS_STRAIN *strain)
{
	int dim, ssnum;
	double strain_v[6];
	double **B_matrix = ws_mat_6_3N;

	if(element.label < 0) return(ARG_ERROR_RC);

	if((dim = element_dim(element.type)) <= 1) return(ARG_ERROR_RC);
	if((ssnum = element_ssnum(element.type)) <= 0) return(ARG_ERROR_RC);

	/* B_matrix, disp_v */
	RC_TRY( make_B_matrix(element, local_point,
	                      node_location, B_matrix, NULL) );

	/* B_matrix, disp_v -> strain_v */
	mul_matrix_vect(ssnum, dim*element.node_num, B_matrix, disp_v, strain_v);

	/* strain_v -> strain */
	init_stress_strain(strain);
	strain->element = element.label;
	RC_TRY( vect2stress_strain(ssnum, strain_v, strain) );

	return(NORMAL_RC);
}


/* 応力、ひずみの計算（熱歪のfact倍を差し引いた値）  */
/* stress, strain が NULL なら代入しない             */
RC
stress_strain_therm (int sol_ss_flag, NODE_ARRAY node,
                     ELEMENT_ARRAY element, DISP_ARRAY disp,
                     MATERIAL_PROP_ARRAY material,
                     PHYSICAL_PROP_ARRAY physical,
                     BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                     double fact,
                     STRESS_ARRAY *stress, STRAIN_ARRAY *strain)
{
	int ii1, ii2;
	static int init_flag = 0;
	static double **local_points;
	STRESS_STRAIN local_stress, local_strain;
	VECT3DR *stress_sum = NULL;
	VECT3DR *strain_sum = NULL;
	int *count = NULL;
	VECT3DR init_v = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	int node_index;

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
		if( (stress_sum == NULL) ||(strain_sum == NULL) ||(count == NULL) )
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
			RC_TRY( local_points_stress_strain_therm(element.array[ii1],
			                                         local_points[0],
			                                         node, disp, material,
			                                         physical, temp, def_temp,
			                                         fact, &local_stress,
			                                         &local_strain) );
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
		if( (sol_ss_flag & SOL_SS_NODE)
		  ||(sol_ss_flag & SOL_SS_NODE_AV) ){
			RC_TRY( set_local_node_points(element.array[ii1].type,
			                              local_points) );

			for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
				RC_TRY( local_points_stress_strain_therm(element.array[ii1],
				                                         local_points[ii2],
				                                         node, disp, material,
				                                         physical, temp,
				                                         def_temp, fact,
				                                         &local_stress,
				                                         &local_strain) );

				if(sol_ss_flag & SOL_SS_NODE_AV){
					node_index = search_node_label(node,
					                   element.array[ii1].node[ii2]);
					RC_NEG_CHK( node_index );
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
				node_index = search_node_label(node,
				                                   stress->array[ii1].node);
				if(count[node_index] <= 0) return(UNKNOWN_ERROR_RC);
				stress->array[ii1].v
				 = mul_scalar_vect3dr(1.0/(double)(count[node_index]),
				                      stress_sum[node_index]);
			}
		}
		if(strain != NULL){
			for(ii1=0; ii1<strain->size; ii1++){
				if(strain->array[ii1].element < 0) continue;
				node_index = search_node_label(node,
				                                   strain->array[ii1].node);
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


/* 要素の局所ひずみ（熱歪のfact倍を差し引いた値） */
RC
local_points_strain_therm (ELEMENT element, const double *local_point,
                           NODE_ARRAY node, DISP_ARRAY disp,
                           MATERIAL_PROP_ARRAY material,
                           PHYSICAL_PROP_ARRAY physical,
                           BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                           double fact,
                           STRESS_STRAIN *strain)
{
	int ii1;
	int dim, ssnum;
	double disp_v[3*MAX_NODE];
	double strain_v[6];
	static double **node_location = NULL;
	static double **B_matrix = NULL;
	static int init_flag = 0;
	double node_temp_array[MAX_NODE];
	double N[MAX_NODE];
	double local_temp;
	int m_index;
	MATERIAL_PROP mat;

	if(element.label < 0) return(ARG_ERROR_RC);

	if(init_flag == 0){
		RC_TRY( allocate2D(MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(6, 3*MAX_NODE, &B_matrix) );

		init_flag = 1;
	}

	for(ii1=0; ii1<element.node_num; ii1++){
		node_temp_array[ii1] = node_temp(temp, def_temp,
		                                 element.node[ii1]);
	}
	RC_TRY( set_N_vector(element.type, local_point, N) );
	local_temp = 0.0;
	for(ii1=0; ii1<element.node_num; ii1++){
		local_temp += node_temp_array[ii1] * N[ii1];
	}
	RC_NEG_CHK( m_index = search_material_prop_element(material,
	                                            physical, &element) );
	mat = material.array[m_index];
	RC_TRY( fill_material(&mat) );

	if( (dim = element_dim(element.type)) <= 1) return(ARG_ERROR_RC);
	if( (ssnum = element_ssnum(element.type)) <= 0) return(ARG_ERROR_RC);

	init_stress_strain(strain);
	strain->element = element.label;

	/* B_matrix, disp_v */
	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( make_B_matrix(element, local_point,
	                      node_location, B_matrix, NULL) );
	RC_TRY( fill_disp_vector(element, disp, disp_v) );

	/* B_matrix, disp_v -> strain_v */
	mul_matrix_vect(ssnum, dim*element.node_num, B_matrix, disp_v, strain_v);
	if(dim == 2){
		strain_v[0] -= fact * local_temp*mat.alpha_vect.x;
		strain_v[1] -= fact * local_temp*mat.alpha_vect.y;
		strain_v[2] -= fact * local_temp*mat.alpha_vect.xy;
	}else{   /* dim == 3 */
		strain_v[0] -= fact * local_temp*mat.alpha_vect.x;
		strain_v[1] -= fact * local_temp*mat.alpha_vect.y;
		strain_v[2] -= fact * local_temp*mat.alpha_vect.z;
		strain_v[3] -= fact * local_temp*mat.alpha_vect.yz;
		strain_v[4] -= fact * local_temp*mat.alpha_vect.zx;
		strain_v[5] -= fact * local_temp*mat.alpha_vect.xy;
	}

	/* strain_v, stress_v -> strain, stress */
	RC_TRY( vect2stress_strain(ssnum, strain_v, strain) );

	return(NORMAL_RC);
}


/* 要素の局所ひずみ -> 局所応力を計算 */
RC
local_strain2local_stress (ELEMENT element,
                           MATERIAL_PROP_ARRAY material,
                           PHYSICAL_PROP_ARRAY physical,
                           STRESS_STRAIN strain, STRESS_STRAIN *stress)
{
	int ssnum;
	double strain_v[6];
	double stress_v[6];
	double **D_matrix = NULL;

	if(element.label < 0) return(ARG_ERROR_RC);

	RC_TRY( allocate2D(6, 6, &D_matrix) );

	if( (ssnum = element_ssnum(element.type)) <= 0) return(ARG_ERROR_RC);

	init_stress_strain(stress);
	stress->element = element.label;
	RC_TRY( stress_strain2vect(ssnum, strain, strain_v) );

	/* D_matrix */
	RC_TRY( get_D_matrix(element, material, physical, D_matrix) );

	/* D_matrix, strain_v -> stress_v */
	mul_matrix_vect(ssnum, ssnum, D_matrix, strain_v, stress_v);

	/* stress_v -> stress */
	RC_TRY( vect2stress_strain(ssnum, stress_v, stress) );

	RC_TRY( free2D(6, 6, &D_matrix) );

	return(NORMAL_RC);
}


/* 要素の局所応力、ひずみ（熱歪のfact倍を差し引いた値） */
RC
local_points_stress_strain_therm (ELEMENT element, const double *local_point,
                                  NODE_ARRAY node, DISP_ARRAY disp,
                                  MATERIAL_PROP_ARRAY material,
                                  PHYSICAL_PROP_ARRAY physical,
                                  BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                                  double fact,
                                  STRESS_STRAIN *stress,
                                  STRESS_STRAIN *strain)
{
	int ii1;
	int dim, ssnum;
	double disp_v[3*MAX_NODE];
	double strain_v[6];
	double stress_v[6];
	static double **node_location = NULL;
	static double **B_matrix = NULL;
	static double **D_matrix = NULL;
	static int init_flag = 0;
	double node_temp_array[MAX_NODE];
	double N[MAX_NODE];
	double local_temp;
	int m_index;
	MATERIAL_PROP mat;

	if(element.label < 0) return(ARG_ERROR_RC);

	if(init_flag == 0){
		RC_TRY( allocate2D(MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(6, 3*MAX_NODE, &B_matrix) );
		RC_TRY( allocate2D(6, 6, &D_matrix) );

		init_flag = 1;
	}

	for(ii1=0; ii1<element.node_num; ii1++){
		node_temp_array[ii1] = node_temp(temp, def_temp,
		                                 element.node[ii1]);
	}
	RC_TRY( set_N_vector(element.type, local_point, N) );
	local_temp = 0.0;
	for(ii1=0; ii1<element.node_num; ii1++){
		local_temp += node_temp_array[ii1] * N[ii1];
	}
	RC_NEG_CHK( m_index = search_material_prop_element(material,
	                                            physical, &element) );
	mat = material.array[m_index];
	RC_TRY( fill_material(&mat) );

	if( (dim = element_dim(element.type)) <= 1) return(ARG_ERROR_RC);
	if( (ssnum = element_ssnum(element.type)) <= 0) return(ARG_ERROR_RC);

	init_stress_strain(stress);
	stress->element = element.label;

	init_stress_strain(strain);
	strain->element = element.label;

	/* B_matrix, disp_v */
	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( make_B_matrix(element, local_point,
	                      node_location, B_matrix, NULL) );
	RC_TRY( fill_disp_vector(element, disp, disp_v) );

	/* B_matrix, disp_v -> strain_v */
	mul_matrix_vect(ssnum, dim*element.node_num, B_matrix, disp_v, strain_v);
	if(dim == 2){
		strain_v[0] -= fact * local_temp*mat.alpha_vect.x;
		strain_v[1] -= fact * local_temp*mat.alpha_vect.y;
		strain_v[2] -= fact * local_temp*mat.alpha_vect.xy;
	}else{   /* dim == 3 */
		strain_v[0] -= fact * local_temp*mat.alpha_vect.x;
		strain_v[1] -= fact * local_temp*mat.alpha_vect.y;
		strain_v[2] -= fact * local_temp*mat.alpha_vect.z;
		strain_v[3] -= fact * local_temp*mat.alpha_vect.yz;
		strain_v[4] -= fact * local_temp*mat.alpha_vect.zx;
		strain_v[5] -= fact * local_temp*mat.alpha_vect.xy;
	}

	/* D_matrix */
	RC_TRY( get_D_matrix(element, material, physical, D_matrix) );

	/* D_matrix, strain_v -> stress_v */
	mul_matrix_vect(ssnum, ssnum, D_matrix, strain_v, stress_v);

	/* strain_v, stress_v -> strain, stress */
	RC_TRY( vect2stress_strain(ssnum, strain_v, strain) );
	RC_TRY( vect2stress_strain(ssnum, stress_v, stress) );

	return(NORMAL_RC);
}


/* 要素の局所応力、ひずみ */
RC
local_points_stress_strain (ELEMENT element, const double *local_point,
                            NODE_ARRAY node, DISP_ARRAY disp,
                            MATERIAL_PROP_ARRAY material,
                            PHYSICAL_PROP_ARRAY physical,
                            STRESS_STRAIN *stress, STRESS_STRAIN *strain)
{
	int dim, ssnum;
	double disp_v[3*MAX_NODE];
	double strain_v[6];
	double stress_v[6];
	static double **node_location = NULL;
	static double **B_matrix = NULL;
	static double **D_matrix = NULL;
	static int init_flag = 0;

	if(element.label < 0) return(ARG_ERROR_RC);

	if(init_flag == 0){
		RC_TRY( allocate2D(MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(6, 3*MAX_NODE, &B_matrix) );
		RC_TRY( allocate2D(6, 6, &D_matrix) );

		init_flag = 1;
	}

	if( (dim = element_dim(element.type)) <= 1) return(ARG_ERROR_RC);
	if( (ssnum = element_ssnum(element.type)) <= 0) return(ARG_ERROR_RC);

	init_stress_strain(stress);
	stress->element = element.label;

	init_stress_strain(strain);
	strain->element = element.label;

	/* B_matrix, disp_v */
	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( make_B_matrix(element, local_point,
	                      node_location, B_matrix, NULL) );
	RC_TRY( fill_disp_vector(element, disp, disp_v) );

	/* B_matrix, disp_v -> strain_v */
	mul_matrix_vect(ssnum, dim*element.node_num, B_matrix, disp_v, strain_v);

	/* D_matrix */
	RC_TRY( get_D_matrix(element, material, physical, D_matrix) );

	/* D_matrix, strain_v -> stress_v */
	mul_matrix_vect(ssnum, ssnum, D_matrix, strain_v, stress_v);

	/* strain_v, stress_v -> strain, stress */
	RC_TRY( vect2stress_strain(ssnum, strain_v, strain) );
	RC_TRY( vect2stress_strain(ssnum, stress_v, stress) );

	return(NORMAL_RC);
}


VECT3DR
d_mises_stress (STRESS_STRAIN stress)
{
	double fact;
	VECT3DR ret;

	fact = 0.5*( (stress.v.x - stress.v.y)*(stress.v.x - stress.v.y)
	           + (stress.v.y - stress.v.z)*(stress.v.y - stress.v.z)
	           + (stress.v.z - stress.v.x)*(stress.v.z - stress.v.x)
	           + 6.0*( (stress.v.yz * stress.v.yz)
	                 + (stress.v.zx * stress.v.zx)
	                 + (stress.v.xy * stress.v.xy) ) );
	fact = 0.5*pow(fact+ABS_TOL, -0.5);
	ret.x = fact*(2.0*stress.v.x - stress.v.y - stress.v.z);
	ret.y = fact*(2.0*stress.v.y - stress.v.z - stress.v.x);
	ret.z = fact*(2.0*stress.v.z - stress.v.x - stress.v.y);
	ret.yz = fact*6.0*(stress.v.yz);
	ret.zx = fact*6.0*(stress.v.zx);
	ret.xy = fact*6.0*(stress.v.xy);

#if 0
	/* fact = 0.5*(1.0/(sqrt(0.5*fact+ABS_TOL)+ABS_TOL)); */
	fact = sqrt(2.0)*0.5*(1.0/(sqrt(fact+ABS_TOL)+ABS_TOL));
	/*fact = 0.25*sqrt(2.0)*(1.0/(sqrt(fact+ABS_TOL)+ABS_TOL));*/
	ret.x = fact*(2.0*stress.v.x - stress.v.y - stress.v.z);
	ret.y = fact*(2.0*stress.v.y - stress.v.z - stress.v.x);
	ret.z = fact*(2.0*stress.v.z - stress.v.x - stress.v.y);
	ret.yz = fact*6.0*(stress.v.yz);
	ret.zx = fact*6.0*(stress.v.zx);
	ret.xy = fact*6.0*(stress.v.xy);
#endif

	return(ret);
}

double
mises_stress (STRESS_STRAIN stress)
{
	double ret;

	ret = (stress.v.x - stress.v.y)*(stress.v.x - stress.v.y)
	    + (stress.v.y - stress.v.z)*(stress.v.y - stress.v.z)
	    + (stress.v.z - stress.v.x)*(stress.v.z - stress.v.x)
	    + 6.0*((stress.v.yz * stress.v.yz)
	    +(stress.v.zx * stress.v.zx)
	    +(stress.v.xy * stress.v.xy));
	ret = sqrt(0.5*ret+ABS_TOL);

	return(ret);
}


RC
vect2stress_strain (int vect_size, const double *vect, STRESS_STRAIN *ss)
{
	switch(vect_size){
	case 3:
		ss->v.x = vect[0];
		ss->v.y = vect[1];
		ss->v.z = 0.0;
		ss->v.xy = vect[2];
		ss->v.yz = 0.0;
		ss->v.zx = 0.0;
		break;
	case 6:
		ss->v.x = vect[0];
		ss->v.y = vect[1];
		ss->v.z = vect[2];
		ss->v.yz = vect[3];
		ss->v.zx = vect[4];
		ss->v.xy = vect[5];
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


RC
stress_strain2vect (int vect_size, STRESS_STRAIN ss, double *vect)
{
	switch(vect_size){
	case 3:
		vect[0] = ss.v.x;
		vect[1] = ss.v.y;
		vect[2] = ss.v.xy;
		break;
	case 6:
		vect[0] = ss.v.x;
		vect[1] = ss.v.y;
		vect[2] = ss.v.z;
		vect[3] = ss.v.yz;
		vect[4] = ss.v.zx;
		vect[5] = ss.v.xy;
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


RC
fill_disp_vector (ELEMENT element, DISP_ARRAY disp, double *disp_v)
{
	int dim;
	int ii1;

	if((dim = element_dim(element.type)) <= 1) return(ARG_ERROR_RC);

	for(ii1=0; ii1<element.node_num; ii1++){
		int index;

		RC_NEG_CHK( index = search_disp_node(disp, element.node[ii1]) );

		disp_v[dim*ii1] = disp.array[index].v.x;
		disp_v[dim*ii1 + 1] = disp.array[index].v.y;
		if(dim == 3){
			disp_v[dim*ii1 + 2] = disp.array[index].v.z;
		}
	}

	return(NORMAL_RC);
}


/* 剛性マトリックスをLU分解 */
RC
decomp_g_matrix (SKYLINE_MATRIX g_matrix)
{
	if(g_matrix.dim == 3){
		RC_TRY( s_cholesky_decomp3(g_matrix.dof, g_matrix.index1,
		                           g_matrix.index2, g_matrix.array, 1) );
	}else{
		RC_TRY( s_cholesky_decomp(g_matrix.dof, g_matrix.index1,
		                          g_matrix.index2, g_matrix.array, 1) );
	}

	return(NORMAL_RC);
}


/* g_matrix, f_vect -> f_vect, disp */
RC
subst_g_matrix (NODE_ARRAY node, SKYLINE_MATRIX g_matrix,
                double *f_vect, DISP_ARRAY *disp)
{
	int ii1;
	int index;

	RC_TRY( s_cholesky_subst(g_matrix.dof, g_matrix.index1, g_matrix.index2,
	                         g_matrix.array, f_vect, 1) );

	/* f_vect -> disp */
	RC_TRY( allocate_disp_array(node.size, disp) );
	disp->source = FEM_K_SOL;
	for(ii1=0; ii1<node.size; ii1++){
		if(node.array[ii1].label < 0){
			disp->array[ii1].node = -1;
			continue;
		}
		index = node.array[ii1].renum;
		if((index < 0)||(index >= node.size)) return(ARG_ERROR_RC);

		disp->array[ii1].node = node.array[ii1].label;
		disp->array[ii1].v.x = f_vect[g_matrix.dim*index];
		disp->array[ii1].v.y = f_vect[g_matrix.dim*index + 1];
		if(g_matrix.dim == 3){
			disp->array[ii1].v.z = f_vect[g_matrix.dim*index + 2];
		}
	}
	disp->size = node.size;
	RC_TRY( clean_disp_array(disp) );

	return(NORMAL_RC);
}


/* g_matrix, s_vect -> s_vect, temp */
RC
subst_g_matrix_therm (NODE_ARRAY node, SKYLINE_MATRIX g_matrix,
                      double *s_vect, BC_ARRAY *temp)
{
	int ii1;
	int index;

	RC_TRY( s_cholesky_subst(g_matrix.dof, g_matrix.index1, g_matrix.index2,
	                         g_matrix.array, s_vect, 1) );

	/* s_vect -> temp */
	RC_TRY( allocate_bc_array(node.size, temp) );
	temp->source = FEM_K_SOL;
	for(ii1=0; ii1<node.size; ii1++){
		if(node.array[ii1].label < 0){
			temp->array[ii1].node = -1;
			continue;
		}
		index = node.array[ii1].renum;
		if((index < 0)||(index >= node.size)) return(ARG_ERROR_RC);

		temp->array[ii1].node = node.array[ii1].label;
		temp->array[ii1].s_type[0] = BC_TEMP;
		temp->array[ii1].s[0] = s_vect[index];
	}
	temp->size = node.size;
	RC_TRY( clean_bc_array(temp) );

	return(NORMAL_RC);
}


/* 拘束条件を考慮して全体剛性マトリックス、荷重ベクトルを修正 */
RC
modify_g_matrix (NODE_ARRAY node, BC_ARRAY rest,
                 SKYLINE_MATRIX g_matrix, double *f_vect)
{
	int ii1;
	int index;

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	for(ii1=0; ii1<rest.size; ii1++){
		if(rest.array[ii1].node < 0) continue;

		index = node_renum_index(node, rest.array[ii1].node);
		if(index < 0) return(SEEK_ERROR_RC);

		if(rest.array[ii1].v_type.x == BC_FIX){
			RC_TRY( s_cholesky_mod(g_matrix.dof, g_matrix.index1,
			                       g_matrix.index2, g_matrix.array, f_vect,
			                       g_matrix.dim*index,
			                       rest.array[ii1].v.x) );
		}
		if(rest.array[ii1].v_type.y == BC_FIX){
			RC_TRY( s_cholesky_mod(g_matrix.dof, g_matrix.index1,
			                       g_matrix.index2, g_matrix.array, f_vect,
			                       g_matrix.dim*index + 1,
			                       rest.array[ii1].v.y) );
		}
		if((g_matrix.dim == 3)&&(rest.array[ii1].v_type.z == BC_FIX)){
			RC_TRY( s_cholesky_mod(g_matrix.dof, g_matrix.index1,
			                       g_matrix.index2, g_matrix.array, f_vect,
			                       g_matrix.dim*index + 2,
			                       rest.array[ii1].v.z) );
		}
	}

	RC_TRY( mod_g_matrix(g_matrix) );

	return(NORMAL_RC);
}


/* 温度境界条件を考慮して全体剛性マトリックス、ベクトルを修正 */
RC
modify_g_matrix_therm (NODE_ARRAY node, BC_ARRAY temp,
                       SKYLINE_MATRIX g_matrix, double *s_vect)
{
	int ii1, ii2;
	int index;

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	for(ii1=0; ii1<temp.size; ii1++){
		if(temp.array[ii1].node < 0) continue;

		index = node_renum_index(node, temp.array[ii1].node);
		if(index < 0) return(SEEK_ERROR_RC);
		
		for(ii2=0; ii2<BC_SCALAR; ii2++){
			if(temp.array[ii1].s_type[ii2] == BC_TEMP){
				RC_TRY( s_cholesky_mod(g_matrix.dof, g_matrix.index1,
				                      g_matrix.index2, g_matrix.array,
				                      s_vect, index, temp.array[ii1].s[ii2]) );
			}
		}
	}

	RC_TRY( mod_g_matrix(g_matrix) );

	return(NORMAL_RC);
}


/* 拘束条件を考慮して全体剛性マトリックス、荷重ベクトルを修正 */
/* 複数荷重条件版(modify_g_matrix() との共通化の余地あり)     */
RC
modify_g_matrix_n (NODE_ARRAY node, BC_ARRAY rest,
                   SKYLINE_MATRIX g_matrix, double **f_vect, int f_vec_size)
{
	int ii1;
	int index;

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	for(ii1=0; ii1<rest.size; ii1++){
		if(rest.array[ii1].node < 0) continue;

		index = node_renum_index(node, rest.array[ii1].node);
		if(index < 0) return(SEEK_ERROR_RC);

		if(rest.array[ii1].v_type.x == BC_FIX){
			RC_TRY( s_cholesky_mod_n(g_matrix.dof, g_matrix.index1,
			                         g_matrix.index2, g_matrix.array, f_vect,
			                         g_matrix.dim*index,
			                         rest.array[ii1].v.x, f_vec_size) );
		}
		if(rest.array[ii1].v_type.y == BC_FIX){
			RC_TRY( s_cholesky_mod_n(g_matrix.dof, g_matrix.index1,
			                         g_matrix.index2, g_matrix.array, f_vect,
			                         g_matrix.dim*index + 1,
			                         rest.array[ii1].v.y, f_vec_size) );
		}
		if((g_matrix.dim == 3)&&(rest.array[ii1].v_type.z == BC_FIX)){
			RC_TRY( s_cholesky_mod_n(g_matrix.dof, g_matrix.index1,
			                         g_matrix.index2, g_matrix.array, f_vect,
			                         g_matrix.dim*index + 2,
			                         rest.array[ii1].v.z, f_vec_size) );
		}
	}

	RC_TRY( mod_g_matrix(g_matrix) );

	return(NORMAL_RC);
}


/* 代入されてない行(列)の対角要素を 1.0 に */
static RC
mod_g_matrix (SKYLINE_MATRIX g_matrix)
{
	RC rc;
	double *ptr;
	int ii1;

	for(ii1=0; ii1<g_matrix.dof; ii1++){
		if(g_matrix.index1[ii1] == ii1){ /* 対角要素のみ */
			ptr = s_cholesky_ptr(g_matrix.dof, g_matrix.index1,
			                     g_matrix.index2, g_matrix.array,
			                     ii1, ii1, &rc);
			if(rc != NORMAL_RC) return(rc);

			if( nearly_eq(*ptr, 0.0) ){
				*ptr = 1.0;
			}
		}
	}

	return(NORMAL_RC);
}


/* G_PROP の確保 */
static RC
allocate_g_prop (G_PROP *g_prop)
{
	int ii1;

	RC_TRY( allocate2D(MAX_GAUSS_POINT+1, 4, &(g_prop->g_points)) );

	for(ii1=0; ii1<MAX_GAUSS_POINT+1; ii1++){
		RC_TRY( allocate2D(6, 3*MAX_NODE, &(g_prop->B_matrix[ii1])) );
	}

	return(NORMAL_RC);
}


RC
fem_total_dof (NODE_ARRAY node, ELEMENT_ARRAY element, int *total_dof)
{
	int dim;

	/* 問題の次元をチェック */
	dim = analysis_dim(element);
	if(dim <= 1) return(ARG_ERROR_RC);

	*total_dof = dim * node.size;

	return(NORMAL_RC);
}


/* ベクトルの確保と初期化 */
RC
fem_allocate_vector (int dim, NODE_ARRAY node, double **f_vect)
{
	int dof;
	int ii1;

	dof = dim * count_valid_node(node);
	if(dof < 0) return(ARG_ERROR_RC);

	(*f_vect) = (double *)mm_alloc(dof * sizeof(double));
	if((*f_vect) == NULL) return(ALLOC_ERROR_RC);

	for(ii1=0; ii1<dof; ii1++){
		(*f_vect)[ii1] = 0.0;
	}

	return(NORMAL_RC);
}


/* 荷重ベクトルの初期化 */
RC
fem_zero_vector (int dim, NODE_ARRAY node, double f_vect[])
{
	int dof;
	int ii1;

	dof = dim * count_valid_node(node);
	if(dof < 0) return(ARG_ERROR_RC);

	for(ii1=0; ii1<dof; ii1++){
		f_vect[ii1] = 0.0;
	}

	return(NORMAL_RC);
}


RC
fem_free_vector (double **vect)
{
	if((*vect) == NULL) return(ARG_ERROR_RC);

	RC_TRY( mm_free((void *)(*vect)) );
	(*vect) = NULL;

	return(NORMAL_RC);
}


RC
force_vector2bc (const double *f_vect, ELEMENT_ARRAY element,
                 NODE_ARRAY node, BC_ARRAY *bc)
{
	int ii1, ii2;
	int dim;
	int index;
	int nonzero_flag;

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	/* 問題の次元をチェック */
	dim = analysis_dim(element);
	if(dim <= 1) return(ARG_ERROR_RC);

	RC_TRY( allocate_bc_array(0, bc) );

	for(ii1=0; ii1<node.size; ii1++){
		if(node.array[ii1].label < 0) continue;

		RC_NEG_CHK( index = node.array[ii1].renum );
		nonzero_flag = 0;
		for(ii2=0; ii2<dim; ii2++){
			if( !nearly_eq(f_vect[index*dim+ii2], 0.0) ) nonzero_flag = 1;
		}
		if(nonzero_flag){
			RC_TRY( realloc_bc_array(bc) );
			(bc->array)[bc->size].node = node.array[ii1].label;
			(bc->array)[bc->size].set_id = 1;
			(bc->array)[bc->size].v_type.x = BC_FORCE;
			(bc->array)[bc->size].v_type.y = BC_FORCE;
			(bc->array)[bc->size].v_type.z = BC_FORCE;
			(bc->array)[bc->size].v.x = f_vect[index*dim];
			(bc->array)[bc->size].v.y = f_vect[index*dim+1];
			if(dim == 3){
				(bc->array)[bc->size].v.z = f_vect[index*dim+2];
			}
			(bc->size)++;
		}
	}

	(bc->source) = FEM_UNKNOWN;
	RC_TRY( clean_bc_array(bc) );

	return(NORMAL_RC);
}


/* 荷重ベクトルに節点荷重を加算 */
RC
add_nodal_force (int dim, NODE_ARRAY node, BC_ARRAY force, double *f_vect)
{
	int index;
	int ii1;

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	for(ii1=0; ii1<force.size; ii1++){
		if(force.array[ii1].node < 0) continue;

		index = node_renum_index(node, force.array[ii1].node);
		if(index < 0) return(SEEK_ERROR_RC);

		if(force.array[ii1].v_type.x == BC_FORCE){
			f_vect[dim*index] += force.array[ii1].v.x;
		}
		if(force.array[ii1].v_type.y == BC_FORCE){
			f_vect[dim*index + 1] += force.array[ii1].v.y;
		}
		if(dim >= 3){
			if(force.array[ii1].v_type.z == BC_FORCE){
				f_vect[dim*index + 2] += force.array[ii1].v.z;
			}
		}
		if(dim >= 6){
			if(force.array[ii1].v_type.yz == BC_FORCE){
				f_vect[dim*index + 3] += force.array[ii1].v.yz;
			}
			if(force.array[ii1].v_type.zx == BC_FORCE){
				f_vect[dim*index + 4] += force.array[ii1].v.zx;
			}
			if(force.array[ii1].v_type.xy == BC_FORCE){
				f_vect[dim*index + 5] += force.array[ii1].v.xy;
			}
		}
	}

	return(NORMAL_RC);
}


RC
add_init_strain_force (int dim, NODE_ARRAY node, ELEMENT_ARRAY element,
                       MATERIAL_PROP_ARRAY material,
                       PHYSICAL_PROP_ARRAY physical,
                       STRAIN_ARRAY init_strain, double *f_vect)
{
	double strain_load[3*MAX_NODE];
	int index_array[3*MAX_NODE];
	int vect_size;
	int ii1, ii2;

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		RC_TRY( element_strain_load(dim, element.array[ii1], node, material,
		                            physical, init_strain,
		                            index_array, &vect_size, strain_load) );

		for(ii2=0; ii2<vect_size; ii2++){
			f_vect[index_array[ii2]] += strain_load[ii2];
		}
	}

	return(NORMAL_RC);
}


/* 要素の初期ひずみによる荷重 */
static RC
element_strain_load (int gdim, ELEMENT element, NODE_ARRAY node,
                     MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
                     STRAIN_ARRAY init_strain,
                     int *index_array, int *vect_size, double *strain_load)
{
	int dim;
	int ii1, ii2;
	double thickness;
	int m_index;
	VECT3DR strain_v[MAX_NODE];
	const VECT3DR zero_v = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	static int init_flag = 0;
	static double **D_matrix = NULL;
	static double **node_location = NULL;
	static double **trans_B = NULL;
	static double **BtD = NULL;
	static G_PROP g_prop;

	if(init_flag == 0){
		RC_TRY( allocate2D(6, 6, &D_matrix) );
		RC_TRY( allocate2D(MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(3*MAX_NODE, 6, &trans_B) );
		RC_TRY( allocate2D(3*MAX_NODE, 6, &BtD) );
		RC_TRY( allocate_g_prop(&g_prop) );
		init_flag = 1;
	}

	if((dim = element_dim(element.type)) <= 1) return(ARG_ERROR_RC);

	if(gdim < dim) return(ARG_ERROR_RC);
	*vect_size = dim*element.node_num;

	RC_NEG_CHK( m_index = search_material_prop_element(material,
	                                                   physical, &element) );

	/* strain_load を初期化 */
	/* index_array を初期化 */
	/* strain_v を初期化    */
	for(ii1=0; ii1<element.node_num; ii1++){
		int index;
		int str_index;

		index = node_renum_index(node, element.node[ii1]);
		if(index < 0) return(SEEK_ERROR_RC); 

		for(ii2=0; ii2<dim; ii2++){
			strain_load[dim*ii1 + ii2] = 0.0;
			index_array[dim*ii1 + ii2] = gdim*index + ii2;
		}

		str_index = search_strain_element_node(init_strain, element.label,
		                                       element.node[ii1]);
		if(str_index < 0){
			str_index = search_strain_element_node(init_strain,
			                                       element.label, -1);
		}
		if(str_index < 0){
			strain_v[ii1] = zero_v;
		}else{
			strain_v[ii1] = init_strain.array[str_index].v;
		}
	}

	/* 要素のD マトリックス */
	RC_TRY( get_D_matrix(element, material, physical, D_matrix) );
	RC_TRY( get_thickness(element, physical, &thickness) );

	/* ガウスポイントごとの B マトリックス */
	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( set_Gauss_const(element.type, &g_prop) );
	RC_TRY( set_g_prop(element, node_location, &g_prop) );

	for(ii1=0; ii1<g_prop.num_g_point; ii1++){
		VECT3DR local_strain;
		double local_force[3*MAX_NODE];
		double N[MAX_NODE];
		double init_strain[6];

		/* ガウスポイントの初期ひずみを内装して求める */
		RC_TRY( set_N_vector(element.type, g_prop.g_points[ii1], N) );
		local_strain = zero_v;
		for(ii2=0; ii2<element.node_num; ii2++){
			VECT3DR tmp_strain = mul_scalar_vect3dr(N[ii2], strain_v[ii1]);
			local_strain = add_vect3dr(local_strain, tmp_strain);
		}

		/* init_strain の初期化 */
		if(dim == 2){
			init_strain[0] = local_strain.x;
			init_strain[1] = local_strain.y;
			init_strain[2] = local_strain.xy;
		}else{   /* dim == 3 */
			init_strain[0] = local_strain.x;
			init_strain[1] = local_strain.y;
			init_strain[2] = local_strain.z;
			init_strain[3] = local_strain.yz;
			init_strain[4] = local_strain.zx;
			init_strain[5] = local_strain.xy;
		}

		/* B^T => trans_B                  */
		/* trans_B * D => BD               */
		/* BD * init_strain => local_force */
		transpose_matrix(g_prop.B_matrix_row, g_prop.B_matrix_col,
		                 g_prop.B_matrix[ii1], trans_B);
		mul_matrix(g_prop.B_matrix_col, g_prop.B_matrix_row,
		           g_prop.B_matrix_row, trans_B, D_matrix, BtD);
		mul_matrix_vect(g_prop.B_matrix_col, g_prop.B_matrix_row, 
		                BtD, init_strain, local_force);

		/* local_force * detJ * weight * thickness */
		for(ii2=0; ii2<(*vect_size); ii2++){
			strain_load[ii2] += local_force[ii2] * g_prop.det_J[ii1]
			                  * g_prop.weight[ii1] * thickness;
		}
	}

	return(NORMAL_RC);
}


/* 熱ひずみによる荷重を加算 */
RC
add_temp_force (int dim, NODE_ARRAY node, ELEMENT_ARRAY element,
                MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
                BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp, double *f_vect)
{
	double temp_load[3*MAX_NODE];
	int index_array[3*MAX_NODE];
	int vect_size;
	int ii1, ii2;

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;
		if(element.array[ii1].type == ELEM_BEAM1) continue;

		RC_TRY( element_temp_load(dim, element.array[ii1], node, material,
		                          physical, temp, def_temp, index_array,
		                          &vect_size, temp_load) );
		for(ii2=0; ii2<vect_size; ii2++){
			f_vect[index_array[ii2]] += temp_load[ii2];
		}
	}

	return(NORMAL_RC);
}


/* 体積力による荷重を加算 */
RC
add_body_force (int dim, NODE_ARRAY node, ELEMENT_ARRAY element,
                MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
                DEFAULT_BC_ARRAY def_b_force, double *f_vect)
{
	double body_force[3*MAX_NODE];
	int index_array[3*MAX_NODE];
	int vect_size;
	int ii1, ii2;

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	/* 問題の次元をチェック */

	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		RC_TRY( element_body_force(dim, element.array[ii1], node, material,
		                           physical, def_b_force, index_array,
		                           &vect_size, body_force) );
		for(ii2=0; ii2<vect_size; ii2++){
			f_vect[index_array[ii2]] += body_force[ii2];
		}
	}

	return(NORMAL_RC);
}


static RC
element_body_force (int gdim, ELEMENT element, NODE_ARRAY node,
                    MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
                    DEFAULT_BC_ARRAY def_b_force,
                    int *index_array, int *vect_size, double *body_force)
{
	int dim;
	int ii1, ii2;
	double thickness;
	int m_index;
	static int allocate_flag = 0;
	static double **node_location = NULL;
	static G_PROP g_prop;
	VECT3D b_force_vect = {0.0, 0.0, 0.0};

	if(allocate_flag == 0){
		RC_TRY( allocate2D(MAX_NODE, 3, &node_location) );
		RC_TRY( allocate_g_prop(&g_prop) );

		allocate_flag = 1;
	}

	/* vect_size, body_force[], index_array[] をセット */
	if((dim = element_dim(element.type)) <= 1) return(ARG_ERROR_RC);

	if(gdim < dim) return(ARG_ERROR_RC);
	*vect_size = dim*element.node_num;
	for(ii1=0; ii1<element.node_num; ii1++){
		int index;

		index = node_renum_index(node, element.node[ii1]);
		if(index < 0) return(SEEK_ERROR_RC);

		for(ii2=0; ii2<dim; ii2++){
			body_force[dim*ii1 + ii2] = 0.0;
			index_array[dim*ii1 + ii2] = gdim*index + ii2;
		}
	}

	/* 重力(grav)を b_force_vect にセット */
	if(def_b_force.size <= 0) return(NORMAL_RC);  /* 重力が作用していない */
	for(ii1=0; ii1<def_b_force.size; ii1++){
		if(def_b_force.array[ii1].set_id < 0) continue;

		b_force_vect.x += def_b_force.array[ii1].grav.x;
		b_force_vect.y += def_b_force.array[ii1].grav.y;
		b_force_vect.z += def_b_force.array[ii1].grav.z;
	}
	RC_NEG_CHK( m_index = search_material_prop_element(material,
	                                            physical, &element) );
	b_force_vect.x *= material.array[m_index].rho;
	b_force_vect.y *= material.array[m_index].rho;
	b_force_vect.z *= material.array[m_index].rho;

	RC_TRY( get_thickness(element, physical, &thickness) );
	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( set_Gauss_const(element.type, &g_prop) );
	RC_TRY( set_g_prop(element, node_location, &g_prop) );

	for(ii1=0; ii1<g_prop.num_g_point; ii1++){
		double N[MAX_NODE];

		RC_TRY( set_N_vector(element.type, g_prop.g_points[ii1], N) );

		for(ii2=0; ii2<element.node_num; ii2++){
			/* b_force_vect * thickness * N * |J| * weight * thickness */
			body_force[ii2*dim]   += b_force_vect.x * thickness * N[ii2]
			                       * g_prop.det_J[ii1] * g_prop.weight[ii1];
			body_force[ii2*dim+1] += b_force_vect.y * thickness * N[ii2]
			                       * g_prop.det_J[ii1] * g_prop.weight[ii1];
			if(dim == 3){
				body_force[ii2*dim+2] += b_force_vect.z * thickness * N[ii2]
				                       * g_prop.det_J[ii1] * g_prop.weight[ii1];
			}
		}
	}

	return(NORMAL_RC);
}


/* 要素の熱ひずみによる荷重 */
static RC
element_temp_load (int gdim, ELEMENT element, NODE_ARRAY node,
                   MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
                   BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                   int *index_array, int *vect_size, double *temp_load)
{
	int dim;
	int ii1, ii2;
	double thickness;
	int m_index;
	double node_temp_array[MAX_NODE];
	static int init_flag = 0;
	static double **D_matrix = NULL;
	static double **node_location = NULL;
	static double **trans_B = NULL;
	static double **BtD = NULL;
	static G_PROP g_prop;
	MATERIAL_PROP mat;

	if(init_flag == 0){
		RC_TRY( allocate2D(6, 6, &D_matrix) );
		RC_TRY( allocate2D(MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(3*MAX_NODE, 6, &trans_B) );
		RC_TRY( allocate2D(3*MAX_NODE, 6, &BtD) );
		RC_TRY( allocate_g_prop(&g_prop) );
		init_flag = 1;
	}

	if((dim = element_dim(element.type)) <= 1) return(ARG_ERROR_RC);

	if(gdim < dim) return(ARG_ERROR_RC);
	*vect_size = dim*element.node_num;

	m_index = search_material_prop_element(material, physical, &element);
	if(m_index < 0) return(SEEK_ERROR_RC);

	mat = material.array[m_index];
	RC_TRY( fill_material(&mat) );

	/* temp_load を初期化       */
	/* index_array を初期化     */
	/* node_temp_array を初期化 */
	for(ii1=0; ii1<element.node_num; ii1++){
		int index = node_renum_index(node, element.node[ii1]);
		if(index < 0) return(SEEK_ERROR_RC);

		for(ii2=0; ii2<dim; ii2++){
			temp_load[dim*ii1 + ii2] = 0.0;
			index_array[dim*ii1 + ii2] = gdim*index + ii2;
		}
		node_temp_array[ii1] = node_temp(temp, def_temp, element.node[ii1]);
	}

	/* 要素のD マトリックス */
	RC_TRY( get_D_matrix(element, material, physical, D_matrix) );
	RC_TRY( get_thickness(element, physical, &thickness) );

	/* ガウスポイントごとの B マトリックス */
	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( set_Gauss_const(element.type, &g_prop) );
	RC_TRY( set_g_prop(element, node_location, &g_prop) );

	for(ii1=0; ii1<g_prop.num_g_point; ii1++){
		double local_temp;
		double local_force[3*MAX_NODE];
		double N[MAX_NODE];
		double init_strain[6];
		int ii2;

		/* ガウスポイントの温度を内装して求める */
		RC_TRY( set_N_vector(element.type, g_prop.g_points[ii1], N) );
		local_temp = 0.0;
		for(ii2=0; ii2<element.node_num; ii2++){
			local_temp += node_temp_array[ii2] * N[ii2];
		}

		/* init_strain の初期化 */
		if(dim == 2){
			init_strain[0] = local_temp*mat.alpha_vect.x;
			init_strain[1] = local_temp*mat.alpha_vect.y;
			init_strain[2] = local_temp*mat.alpha_vect.xy;
		}else{   /* dim == 3 */
			init_strain[0] = local_temp*mat.alpha_vect.x;
			init_strain[1] = local_temp*mat.alpha_vect.y;
			init_strain[2] = local_temp*mat.alpha_vect.z;
			init_strain[3] = local_temp*mat.alpha_vect.yz;
			init_strain[4] = local_temp*mat.alpha_vect.zx;
			init_strain[5] = local_temp*mat.alpha_vect.xy;
		}

		/* B^T => trans_B */
		/* trans_B * D => BD */
		/* BD * init_strain => local_force */
		transpose_matrix(g_prop.B_matrix_row, g_prop.B_matrix_col,
		                 g_prop.B_matrix[ii1], trans_B);
		mul_matrix(g_prop.B_matrix_col, g_prop.B_matrix_row,
		           g_prop.B_matrix_row, trans_B, D_matrix, BtD);
		mul_matrix_vect(g_prop.B_matrix_col, g_prop.B_matrix_row, 
		                BtD, init_strain, local_force);

		/* local_force * detJ * weight * thickness */
		for(ii2=0; ii2<(*vect_size); ii2++){
			temp_load[ii2] += local_force[ii2] * g_prop.det_J[ii1]
			                * g_prop.weight[ii1] * thickness;
		}
	}

	return(NORMAL_RC);
}
                         

/* 節点の温度 */
double
node_temp (BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp, int node_label)
{
	int index;
	int ii1;
	double d_temp;

	/* 温度のデフォルト値は def_temp の最初の値もしくは 0.0 */
	if(def_temp.size <= 0){
		d_temp = 0.0;
	}else{
		d_temp = def_temp.array[0].temp;
	}

	index = search_bc_node(temp, node_label);
	if(index < 0) return(d_temp);

	for(ii1=0; ii1<BC_SCALAR; ii1++){
		if(temp.array[index].s_type[ii1] == BC_TEMP){
			return(temp.array[index].s[ii1]);
		}
	}

	return(d_temp);
}


/* 全体剛性マトリックスを作成 */
RC
make_g_matrix (NODE_ARRAY node, ELEMENT_ARRAY element,
               MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
               SKYLINE_MATRIX *g_matrix)
{
	int ii1;
	int dim;
	ELEM_MATRIX elem_matrix;

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	/* 問題の次元(g_matrix->dim) */
	dim = analysis_dim(element);
	if(dim <= 1) return(ARG_ERROR_RC);

	/* g_matrix の確保と初期化 */
	RC_TRY( allocate_g_matrix(dim, node, element, g_matrix) );

	for(ii1=0;ii1<element.size;ii1++){
		if(element.array[ii1].label < 0) continue;

		/* element.array[ii1], node, material, physical -> elem_matrix */
		RC_TRY( make_elem_matrix(element.array[ii1], node, material,
		                         physical, &elem_matrix) );

		/* elem_matrix -> g_matrix */
		RC_TRY( impose_elem_matrix(*g_matrix, elem_matrix, 1.0) );
		RC_TRY( free_elem_matrix(&elem_matrix) );
	}

	return(NORMAL_RC);
}

RC
make_g_matrix_add_density (NODE_ARRAY node, ELEMENT_ARRAY element,
               MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
               SKYLINE_MATRIX *g_matrix,double *density)
{
	int ii1;
	int dim;
	ELEM_MATRIX elem_matrix;

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	/* 問題の次元(g_matrix->dim) */
	dim = analysis_dim(element);
	if(dim <= 1) return(ARG_ERROR_RC);

	/* g_matrix の確保と初期化 */
	RC_TRY( allocate_g_matrix(dim, node, element, g_matrix) );

	for(ii1=0;ii1<element.size;ii1++){
		if(element.array[ii1].label < 0) continue;

		/* element.array[ii1], node, material, physical -> elem_matrix */
		RC_TRY( make_elem_matrix_add_density(element.array[ii1], node, material,
		                         physical, &elem_matrix, density) );

		/* elem_matrix -> g_matrix */
		RC_TRY( impose_elem_matrix(*g_matrix, elem_matrix, 1.0) );
		RC_TRY( free_elem_matrix(&elem_matrix) );
	}

	return(NORMAL_RC);
}


/* 全体剛性マトリックスの確保と問題の次元をチェック */
RC
allocate_g_matrix (int dim, NODE_ARRAY node, ELEMENT_ARRAY element,
                   SKYLINE_MATRIX *g_matrix)
{
	g_matrix->dim = dim;

	g_matrix->dof = (g_matrix->dim)*node.size;
	g_matrix->index1 = (long *)mm_alloc(g_matrix->dof * sizeof(long));
	if((g_matrix->index1) == NULL) return(ALLOC_ERROR_RC);

	RC_TRY( fill_index1(g_matrix->dim, node, element, g_matrix->index1) );

	RC_TRY( s_cholesky_alloc(g_matrix->dof, g_matrix->index1,
	                         &(g_matrix->index2), &(g_matrix->array),
	                         &(g_matrix->array_size)) );

	return(NORMAL_RC);
}


/* index1 を適切な値で初期化する           */
/* index1[i] ...i 列目の最初の非ゼロ行番号 */
static RC
fill_index1 (int dim, NODE_ARRAY node, ELEMENT_ARRAY element, long *index1)
{
	int ii1, ii2, ii3, ii4;
	long pos1, pos2;

	for(ii1=0; ii1<dim*(node.size); ii1++){
		index1[ii1] = ii1;
	}

	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
			pos1 = dim*node_renum_index(node, element.array[ii1].node[ii2]);
			if(pos1 < 0) return(SEEK_ERROR_RC);

			for(ii3=0; ii3<element.array[ii1].node_num; ii3++){
				pos2 = dim*node_renum_index(node,element.array[ii1].node[ii3]);
				if(pos2 < 0) return(SEEK_ERROR_RC);
				for(ii4=0; ii4<dim; ii4++){
					if(pos1 < index1[pos2+ii4]) index1[pos2+ii4] = pos1;
				}
			}
		}
	}

	return(NORMAL_RC);
}


/* 要素剛性マトリックスの作成 */
RC
make_elem_matrix (ELEMENT element, NODE_ARRAY node,
                  MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
                  ELEM_MATRIX *elem_matrix)
{
	int ii1, ii2, ii3;
	double thickness;
	static double **D_matrix = NULL;
	static double **BtDB = NULL;
	static int init_flag = 0;
	static double **node_location = NULL;
	static G_PROP g_prop;

	if(init_flag == 0){
		RC_TRY( allocate2D(6, 6, &D_matrix) );
		RC_TRY( allocate2D(3*MAX_NODE, 3*MAX_NODE, &BtDB) );
		RC_TRY( allocate2D(MAX_NODE, 3, &node_location) );
		RC_TRY( allocate_g_prop(&g_prop) );
		init_flag = 1;
	}

	if((elem_matrix->dim = element_dim(element.type)) <= 0){
		return(ARG_ERROR_RC);
	}

	if(element.label < 0) return(ARG_ERROR_RC);
	elem_matrix->element = element.label;
	elem_matrix->size = (element.node_num) * (elem_matrix->dim);

	/* elem_matrix->array[][] の確保 */
	RC_TRY( allocate2D(elem_matrix->size,
	                   elem_matrix->size, &(elem_matrix->matrix)) );

	/* elem_matirx->index[] の確保と代入 */

	elem_matrix->index = mm_alloc(sizeof(int)*(elem_matrix->size));
	if(elem_matrix->index == NULL) return(ALLOC_ERROR_RC);

	for(ii1=0; ii1<element.node_num; ii1++){
		int index = node_renum_index(node, element.node[ii1]);
		if(index < 0) return(SEEK_ERROR_RC);

		for(ii2=0; ii2<(elem_matrix->dim); ii2++){
			elem_matrix->index[ii1*(elem_matrix->dim)+ii2]
			             = index*(elem_matrix->dim) + ii2;
		}
	}

	/* D マトリックス */
	RC_TRY( get_D_matrix(element, material, physical, D_matrix) );
	RC_TRY( get_thickness(element, physical, &thickness) );

	//printf("%lf\n",thickness);

	/* B マトリックス */
	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( set_Gauss_const(element.type, &g_prop) );
	RC_TRY( set_g_prop(element, node_location, &g_prop) );
	RC_TRY( RI_B_matrix(element, g_prop) );

	for(ii1=0; ii1<g_prop.num_g_point; ii1++){
		/* B^T * D * B => BtDB */
		mul_matrix_AtBA(g_prop.B_matrix_row, g_prop.B_matrix_col,
		                g_prop.B_matrix[ii1], D_matrix, BtDB);

		/* BDB * detJ * weight * thickness */
		for(ii2=0; ii2<(elem_matrix->size); ii2++){
			for(ii3=0; ii3<(elem_matrix->size); ii3++){
				elem_matrix->matrix[ii2][ii3]
				      += BtDB[ii2][ii3] * g_prop.det_J[ii1]
				       * g_prop.weight[ii1] * thickness;
			}
		}
	}

	return(NORMAL_RC);
}

RC
make_elem_matrix_add_density (ELEMENT element, NODE_ARRAY node,
                  MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
                  ELEM_MATRIX *elem_matrix,double *density)
{
	int ii1, ii2, ii3;
	int index1;
	double phi;
	double thickness;
	static double **D_matrix = NULL;
	static double **BtDB = NULL;
	static int init_flag = 0;
	static double **node_location = NULL;
	static G_PROP g_prop;
	double *N;

	if(init_flag == 0){
		RC_TRY( allocate2D(6, 6, &D_matrix) );
		RC_TRY( allocate2D(3*MAX_NODE, 3*MAX_NODE, &BtDB) );
		RC_TRY( allocate2D(MAX_NODE, 3, &node_location) );
		RC_TRY( allocate_g_prop(&g_prop) );
		init_flag = 1;
	}

	RC_TRY( allocate1D(MAX_NODE,&N));

	if((elem_matrix->dim = element_dim(element.type)) <= 0){
		return(ARG_ERROR_RC);
	}

	if(element.label < 0) return(ARG_ERROR_RC);
	elem_matrix->element = element.label;
	elem_matrix->size = (element.node_num) * (elem_matrix->dim);

	/* elem_matrix->array[][] の確保 */
	RC_TRY( allocate2D(elem_matrix->size,
	                   elem_matrix->size, &(elem_matrix->matrix)) );

	/* elem_matirx->index[] の確保と代入 */

	elem_matrix->index = mm_alloc(sizeof(int)*(elem_matrix->size));
	if(elem_matrix->index == NULL) return(ALLOC_ERROR_RC);

	for(ii1=0; ii1<element.node_num; ii1++){
		int index = node_renum_index(node, element.node[ii1]);
		if(index < 0) return(SEEK_ERROR_RC);

		for(ii2=0; ii2<(elem_matrix->dim); ii2++){
			elem_matrix->index[ii1*(elem_matrix->dim)+ii2]
			             = index*(elem_matrix->dim) + ii2;
		}
	}

	/* D マトリックス */
	RC_TRY( get_D_matrix(element, material, physical, D_matrix) );
	RC_TRY( get_thickness(element, physical, &thickness) );

	/* B マトリックス */
	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( set_Gauss_const(element.type, &g_prop) );
	RC_TRY( set_g_prop(element, node_location, &g_prop) );
	RC_TRY( RI_B_matrix(element, g_prop) );


	for(ii1=0; ii1<g_prop.num_g_point; ii1++){
		/* B^T * D * B => BtDB */
		mul_matrix_AtBA(g_prop.B_matrix_row, g_prop.B_matrix_col,
		                g_prop.B_matrix[ii1], D_matrix, BtDB);

		phi = 0;

		for(ii2=0; ii2<element.node_num; ii2++){
			index1 = search_node_label(node,element.node[ii2]);
			RC_TRY( set_N_vector(element.type, g_prop.g_points[ii1], N) );
			phi += density[index1]* N[ii2];
		}		

		//printf("%lf\n",phi);



		/* BDB * detJ * weight * thickness */
		for(ii2=0; ii2<(elem_matrix->size); ii2++){
			for(ii3=0; ii3<(elem_matrix->size); ii3++){
				elem_matrix->matrix[ii2][ii3]
				      += BtDB[ii2][ii3] * g_prop.det_J[ii1]
				       * g_prop.weight[ii1] * thickness
						   * phi * phi;
			}
		}
	}

	RC_TRY( free1D(MAX_NODE,&N));

	return(NORMAL_RC);
}


/* 要素剛性マトリックス(熱伝導解析用) */
RC
make_elem_matrix_therm (ELEMENT element, NODE_ARRAY node,
                        MATERIAL_PROP_ARRAY material,
                        PHYSICAL_PROP_ARRAY physical,
                        ELEM_MATRIX *elem_matrix)
{
	int ii1, ii2, ii3;
	static int init_flag = 0;
	static double **node_location = NULL;
	static double **gauss_points = NULL;
	static double weight[MAX_GAUSS_POINT];
	static double **dN = NULL;
	static double **dNxyz = NULL;
	static double **jacobi = NULL;
	static double **inverse_J = NULL;
	static double **D_matrix = NULL;
	static double **BtDB = NULL;
	static int num_point;
	int elem_dim;
	double det_J, thickness;

	if(init_flag == 0){
		RC_TRY( allocate2D(MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(MAX_GAUSS_POINT, 4, &gauss_points) );
		RC_TRY( allocate2D(3, MAX_NODE, &dN) );
		RC_TRY( allocate2D(3, MAX_NODE, &dNxyz) );
		RC_TRY( allocate2D(3, 3, &jacobi) );
		RC_TRY( allocate2D(3, 3, &inverse_J) );
		RC_TRY( allocate2D(3, 3, &D_matrix) );
		RC_TRY( allocate2D(MAX_NODE, MAX_NODE, &BtDB) );

		init_flag = 1;
	}

	if(element.label < 0) return(ARG_ERROR_RC);
	elem_matrix->element = element.label;
	elem_dim = element_dim(element.type);

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	RC_TRY( Gauss_const_default(element.type, gauss_points,
	                            weight, &num_point) );
	RC_TRY( node_location_matrix(&element, node, node_location) );

	elem_matrix->dim = 1;
	elem_matrix->size = element.node_num;

	RC_TRY( get_D_matrix_therm(element, material, physical, D_matrix) );
	RC_TRY( get_thickness(element, physical, &thickness) );

	/* elem_matrix->array[][] の確保 */
	RC_TRY( allocate2D(elem_matrix->size,
	                   elem_matrix->size, &(elem_matrix->matrix)) );

	/* elem_matirx->index[] の確保と代入 */
	elem_matrix->index = mm_alloc(sizeof(int)*(elem_matrix->size));
	if(elem_matrix->index == NULL) return(ALLOC_ERROR_RC);

	for(ii1=0; ii1<element.node_num; ii1++){
		int index = node_renum_index(node, element.node[ii1]);
		if(index < 0) return(SEEK_ERROR_RC);

		elem_matrix->index[ii1] = index;
	}

	for(ii1=0; ii1<num_point; ii1++){
		RC_TRY( set_dN_matrix(element.type, gauss_points[ii1], dN) );
		mul_matrix(elem_dim, element.node_num, elem_dim,
		           dN, node_location, jacobi);
		RC_TRY( inverse_matrix(elem_dim, jacobi, inverse_J) );
		mul_matrix(elem_dim, elem_dim, element.node_num,
		           inverse_J, dN, dNxyz);
		det_J = determinant(elem_dim, jacobi);

		mul_matrix_AtBA(elem_dim, element.node_num, dNxyz, D_matrix, BtDB);

		for(ii2=0; ii2<element.node_num; ii2++){
			for(ii3=0; ii3<element.node_num; ii3++){
				elem_matrix->matrix[ii2][ii3]
				     += weight[ii1] * det_J * BtDB[ii2][ii3] * thickness;
			}
		}
	}

	return(NORMAL_RC);
}


/* バネ定数 */
RC
get_spring_const (PHYSICAL_PROP_ARRAY physical,
                  int p_label, double *spring_const)
{
	int p_index;

	if(physical.source == FEM_NASTRAN){
		p_index = search_physical_prop_label(physical, p_label);
		if(p_index < 0) return(SEEK_ERROR_RC);

		if(physical.array[p_index].i_info[1] == NST_PROP_PELAS){
			*spring_const = physical.array[p_index].d_info[0];
		}else{
			return(IMPLEMENT_ERROR_RC);
		}
	}else{
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


/* 板の厚さ or 梁の断面積 */
RC
get_thickness (ELEMENT element, PHYSICAL_PROP_ARRAY physical,
               double *thickness)
{
	int dim;
	int p_index;

	if((dim = element_dim(element.type)) <= 0) return(ARG_ERROR_RC);

	/* 厚さ */
	*thickness = 1.0;
	if((dim == 2)&&(physical.source == FEM_NASTRAN)){
		p_index = search_physical_prop_label(physical, element.physical);
		if(p_index < 0) return(SEEK_ERROR_RC);

		if(physical.array[p_index].i_info[1] == NST_PROP_PSHELL){
			/* PSHELL of NASTRAN */
			if(physical.array[p_index].d_info[0] > 0.0){
				*thickness = physical.array[p_index].d_info[0];
			}
		}
	}

	return(NORMAL_RC);
}


/* D マトリックス */
RC
get_D_matrix (ELEMENT element, MATERIAL_PROP_ARRAY material,
              PHYSICAL_PROP_ARRAY physical, double **D_matrix)
{
	int ii1, ii2;
	int dim, ssnum;
	int m_index;
	MATERIAL_PROP tmp_mat;

	m_index = search_material_prop_element(material, physical, &element);
	if(m_index < 0) return(SEEK_ERROR_RC);

	if((dim = element_dim(element.type)) <= 0) return(ARG_ERROR_RC);
	if( (ssnum = element_ssnum(element.type)) <= 0) return(ARG_ERROR_RC);

	tmp_mat = material.array[m_index];
	if(dim == 2){
		tmp_mat.ana_type = DEF_ANA_2D;
	}else{   /* dim == 3 */
		tmp_mat.ana_type = ANA_3D;
	}
	RC_TRY( fill_material(&tmp_mat) );

	/* tmp_mat.D_matrix -> D_matrix */
	for(ii1=0;ii1<ssnum;ii1++){
		for(ii2=0;ii2<ssnum;ii2++){
			D_matrix[ii1][ii2] = tmp_mat.D_matrix[ii1][ii2];
		}
	}

	return(NORMAL_RC);
}


/* D マトリックス（熱伝導解析用） */
RC
get_D_matrix_therm (ELEMENT element, MATERIAL_PROP_ARRAY material,
                    PHYSICAL_PROP_ARRAY physical, double **D_matrix)
{
	int ii1, ii2;
	int dim;
	int m_index;

	m_index = search_material_prop_element(material, physical, &element);
	if(m_index < 0) return(SEEK_ERROR_RC);

	if((dim = element_dim(element.type)) <= 0) return(ARG_ERROR_RC);

	switch(material.array[m_index].mat_type){
	case MAT_ISOTROPIC:
		for(ii1=0;ii1<dim;ii1++){
			for(ii2=0;ii2<dim;ii2++){
				D_matrix[ii1][ii2] = 0.0;
			}
			D_matrix[ii1][ii1] = material.array[m_index].K;
		}
		break;
	case MAT_ORTHOTROPIC:
		get_D_matrix_therm_ortho(material.array[m_index], D_matrix);
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}

/* D マトリックス（熱伝導解析用）異方性 */
void get_D_matrix_therm_ortho(MATERIAL_PROP array, double **D_matrix)
{
	D_matrix[0][0] = array.K_vect.x;
	D_matrix[0][1] = array.K_vect.xy;
	D_matrix[0][2] = array.K_vect.zx;
	D_matrix[1][0] = array.K_vect.xy;
	D_matrix[1][1] = array.K_vect.y;
	D_matrix[1][2] = array.K_vect.yz;
	D_matrix[2][0] = array.K_vect.zx;
	D_matrix[2][1] = array.K_vect.yz;
	D_matrix[2][2] = array.K_vect.z;
}


/* 次数低減積分の処理 (B_matrix を上書き) */
static RC
RI_B_matrix (ELEMENT elem, G_PROP g_prop)
{
/*
	int ii1, ii2;
	int dim;

	switch(elem.type){
	case ELEM_QUAD1:
		dim = 2;
		for(ii1=0; ii1<g_prop.num_g_point; ii1++){
			for(ii2=0;ii2<elem.node_num;ii2++){
				g_prop.B_matrix[ii1][2][dim*ii2]
				  = g_prop.B_matrix[MAX_GAUSS_POINT][2][dim*ii2];
				g_prop.B_matrix[ii1][2][dim*ii2+1]
				  = g_prop.B_matrix[MAX_GAUSS_POINT][2][dim*ii2+1];
			}
		}
		break;
	case ELEM_PENTA1:
	case ELEM_HEXA1:
		dim = 3;
		for(ii1=0; ii1<g_prop.num_g_point; ii1++){
			for(ii2=0;ii2<elem.node_num;ii2++){
				g_prop.B_matrix[ii1][3][dim*ii2]
				  = g_prop.B_matrix[MAX_GAUSS_POINT][3][dim*ii2];
				g_prop.B_matrix[ii1][3][dim*ii2+1]
				  = g_prop.B_matrix[MAX_GAUSS_POINT][3][dim*ii2+1];
				g_prop.B_matrix[ii1][3][dim*ii2+2]
				  = g_prop.B_matrix[MAX_GAUSS_POINT][3][dim*ii2+2];
				g_prop.B_matrix[ii1][4][dim*ii2]
				  = g_prop.B_matrix[MAX_GAUSS_POINT][4][dim*ii2];
				g_prop.B_matrix[ii1][4][dim*ii2+1]
				  = g_prop.B_matrix[MAX_GAUSS_POINT][4][dim*ii2+1];
				g_prop.B_matrix[ii1][4][dim*ii2+2]
				  = g_prop.B_matrix[MAX_GAUSS_POINT][4][dim*ii2+2];
				g_prop.B_matrix[ii1][5][dim*ii2]
				  = g_prop.B_matrix[MAX_GAUSS_POINT][5][dim*ii2];
				g_prop.B_matrix[ii1][5][dim*ii2+1]
				  = g_prop.B_matrix[MAX_GAUSS_POINT][5][dim*ii2+1];
				g_prop.B_matrix[ii1][5][dim*ii2+2]
				  = g_prop.B_matrix[MAX_GAUSS_POINT][5][dim*ii2+2];
			}
		}
		break;
	default:
		break;
	}
*/
	return(NORMAL_RC);
}


/* 要素上のガウスポイント毎の情報 */
static RC
set_g_prop (ELEMENT element, double **node_location, G_PROP *g_prop)
{
	int ii1;
	int dim;

	if( (g_prop->B_matrix_row = element_ssnum(element.type)) <= 0){
		return(ARG_ERROR_RC);
	}

	if((dim = element_dim(element.type)) <= 0) return(ARG_ERROR_RC);

	g_prop->B_matrix_col = dim * element.node_num;

	/* 各ガウスポイント */
	for(ii1=0; ii1<(g_prop->num_g_point); ii1++){
		RC_TRY( make_B_matrix(element, g_prop->g_points[ii1], node_location,
		                      g_prop->B_matrix[ii1], &(g_prop->det_J[ii1])) );
	}

	/* 要素中心 */
	RC_TRY( make_B_matrix(element, g_prop->g_points[MAX_GAUSS_POINT],
	                      node_location, g_prop->B_matrix[MAX_GAUSS_POINT],
	                      &(g_prop->det_J[MAX_GAUSS_POINT])) );

	return(NORMAL_RC);
}


/* B マトリックスと |J| */
RC
make_B_matrix (ELEMENT element, const double *g_point,
               double **node_location, double **B_matrix, double *det_J)
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
	/* inverse_J * dN => dNxy       */
	RC_TRY( set_dN_matrix(element.type, g_point, dN) );
	mul_matrix(dim, element.node_num, dim, dN, node_location, Jacobi);
	RC_TRY( inverse_matrix(dim, Jacobi, inverse_J) );
	mul_matrix(dim, dim, element.node_num, inverse_J, dN, dNxyz);
	RC_TRY( set_B_matrix(dim, element.node_num, dNxyz, B_matrix) );

	if(det_J != NULL){
		*det_J = determinant(dim, Jacobi);
	}

	return(NORMAL_RC);
}


/* 積分点の局所座標と重み */
static RC
set_Gauss_const (ELEM_TYPE type, G_PROP *g_prop)
{
	switch(type){
	case ELEM_TRI1:
	case ELEM_TRI2:
		g_prop->num_g_point = 3;
		RC_TRY( Gauss_const_tri(g_prop->num_g_point,
		                        g_prop->g_points, g_prop->weight) );
		RC_TRY( Gauss_const_tri(1, &(g_prop->g_points[MAX_GAUSS_POINT]),
		                        &(g_prop->weight[MAX_GAUSS_POINT])) );
		break;
	case ELEM_QUAD1:
		g_prop->num_g_point = 4;
		RC_TRY( Gauss_const_quad(g_prop->num_g_point,
		                         g_prop->g_points, g_prop->weight) );
		RC_TRY( Gauss_const_quad(1, &(g_prop->g_points[MAX_GAUSS_POINT]),
		                         &(g_prop->weight[MAX_GAUSS_POINT])) );
		break;
	case ELEM_QUAD2:
	case ELEM_QUAD2L:
		g_prop->num_g_point = 9;
		RC_TRY( Gauss_const_quad(g_prop->num_g_point,
		                         g_prop->g_points, g_prop->weight) );
		RC_TRY( Gauss_const_quad(1, &(g_prop->g_points[MAX_GAUSS_POINT]),
		                         &(g_prop->weight[MAX_GAUSS_POINT])) );
		break;
	case ELEM_TETRA1:
		g_prop->num_g_point = 4;
		RC_TRY( Gauss_const_tetra(g_prop->num_g_point,
		                          g_prop->g_points, g_prop->weight) );
		RC_TRY( Gauss_const_tetra(1, &(g_prop->g_points[MAX_GAUSS_POINT]),
		                          &(g_prop->weight[MAX_GAUSS_POINT])) );
		break;
	case ELEM_TETRA2:
		g_prop->num_g_point = 5;
		RC_TRY( Gauss_const_tetra(g_prop->num_g_point,
		                          g_prop->g_points, g_prop->weight) );
		RC_TRY( Gauss_const_tetra(1, &(g_prop->g_points[MAX_GAUSS_POINT]),
		                          &(g_prop->weight[MAX_GAUSS_POINT])) );
		break;
	case ELEM_PENTA1:
		g_prop->num_g_point = 6;
		RC_TRY( Gauss_const_penta(g_prop->num_g_point,
		                          g_prop->g_points, g_prop->weight) );
		RC_TRY( Gauss_const_penta(1, &(g_prop->g_points[MAX_GAUSS_POINT]),
		                          &(g_prop->weight[MAX_GAUSS_POINT])) );
		break;
	case ELEM_PENTA2:
		g_prop->num_g_point = 9;
		RC_TRY( Gauss_const_penta(g_prop->num_g_point,
		                          g_prop->g_points, g_prop->weight) );
		RC_TRY( Gauss_const_penta(1, &(g_prop->g_points[MAX_GAUSS_POINT]),
		                          &(g_prop->weight[MAX_GAUSS_POINT])) );
		break;
	case ELEM_HEXA1:
		g_prop->num_g_point = 8;
		RC_TRY( Gauss_const_hexa(g_prop->num_g_point,
		                         g_prop->g_points, g_prop->weight) );
		RC_TRY( Gauss_const_hexa(1, &(g_prop->g_points[MAX_GAUSS_POINT]),
		                         &(g_prop->weight[MAX_GAUSS_POINT])) );
		break;
	case ELEM_HEXA2:
		g_prop->num_g_point = 14;
		RC_TRY( Gauss_const_hexa(g_prop->num_g_point,
		                         g_prop->g_points, g_prop->weight) );
		RC_TRY( Gauss_const_hexa(1, &(g_prop->g_points[MAX_GAUSS_POINT]),
		                         &(g_prop->weight[MAX_GAUSS_POINT])) );
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


/* dNxyz を用いて B マトリックスを作成 */
/* dNxyz[0][i]...Ni の x 方向微分      */
/* dNxyz[1][i]...Ni の y 方向微分      */
/* dNxyz[2][i]...Ni の z 方向微分      */
/* dim=2 ...2次元, dim=3 ...3次元      */
static RC
set_B_matrix (int dim, int node_num, double **dNxyz, double **B_matrix)
{
	int ii1;

	switch(dim){
	case 2:
		for(ii1=0;ii1<node_num;ii1++){
			B_matrix[0][dim*ii1]   = dNxyz[0][ii1];
			B_matrix[0][dim*ii1+1] = 0.0;
			B_matrix[1][dim*ii1]   = 0.0;
			B_matrix[1][dim*ii1+1] = dNxyz[1][ii1];
			B_matrix[2][dim*ii1]   = dNxyz[1][ii1];
			B_matrix[2][dim*ii1+1] = dNxyz[0][ii1];
		}

		break;
	case 3:
		for(ii1=0;ii1<node_num;ii1++){
			B_matrix[0][dim*ii1]   = dNxyz[0][ii1];
			B_matrix[0][dim*ii1+1] = 0.0;
			B_matrix[0][dim*ii1+2] = 0.0;
			B_matrix[1][dim*ii1]   = 0.0;
			B_matrix[1][dim*ii1+1] = dNxyz[1][ii1];
			B_matrix[1][dim*ii1+2] = 0.0;
			B_matrix[2][dim*ii1]   = 0.0;
			B_matrix[2][dim*ii1+1] = 0.0;
			B_matrix[2][dim*ii1+2] = dNxyz[2][ii1];
			B_matrix[3][dim*ii1]   = 0.0;
			B_matrix[3][dim*ii1+1] = dNxyz[2][ii1];
			B_matrix[3][dim*ii1+2] = dNxyz[1][ii1];
			B_matrix[4][dim*ii1]   = dNxyz[2][ii1];
			B_matrix[4][dim*ii1+1] = 0.0;
			B_matrix[4][dim*ii1+2] = dNxyz[0][ii1];
			B_matrix[5][dim*ii1]   = dNxyz[1][ii1];
			B_matrix[5][dim*ii1+1] = dNxyz[0][ii1];
			B_matrix[5][dim*ii1+2] = 0.0;
		}
		
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


/* 要素剛性マトリックスを全体剛性マトリックスに組み込む */
RC
impose_elem_matrix (SKYLINE_MATRIX g_matrix,
                    ELEM_MATRIX elem_matrix, double fact)
{
	RC rc;
	int ii1, ii2;
	int ix1, ix2;
	double *ptr;

	for(ii1=0; ii1<elem_matrix.size; ii1++){
		for(ii2=0; ii2<elem_matrix.size; ii2++){
			ix1 = elem_matrix.index[ii1];
			ix2 = elem_matrix.index[ii2];
			if(ix1 < ix2) continue;

			ptr = s_cholesky_ptr(g_matrix.dof, g_matrix.index1,
			                     g_matrix.index2, g_matrix.array,
			                     ix1, ix2, &rc);
			if(rc != NORMAL_RC) return(rc);
			(*ptr) += fact * elem_matrix.matrix[ii1][ii2];
		}
	}

	return(NORMAL_RC);
}


RC
impose_nonzero3 (NONZERO_MATRIX3 matrix, ELEM_MATRIX elem_matrix, double fact)
{
	int ii1, ii2, ii3, ii4;
	int ix1, ix2;
	int col_index;
	if(elem_matrix.dim != 3) return(ARG_ERROR_RC);

	for(ii1=0; ii1<elem_matrix.size; ii1+=3){
		for(ii2=0; ii2<elem_matrix.size; ii2+=3){
			ix1 = elem_matrix.index[ii1]/3;
			ix2 = elem_matrix.index[ii2]/3;
			if(ix1 < ix2) continue;

			if(ix1 >= matrix.size) return(ARG_ERROR_RC);
			col_index = -1;
			for(ii3=0; ii3<matrix.row[ix1].size; ii3++){
				if(matrix.row[ix1].elem[ii3].col == ix2){
					col_index = ii3;
					break;
				}
			}
			if(col_index < 0) return(ARG_ERROR_RC);

			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					matrix.row[ix1].elem[col_index].v[ii3][ii4]
					      += fact * elem_matrix.matrix[ii1+ii3][ii2+ii4];
				}
			}
		}
	}

	return(NORMAL_RC);
}


RC
iccg_nonzero3 (NONZERO_MATRIX3 matrix, const double *f_vect, 
               double **d_vect, const double *init_vect)
{
	int ii1;

	*d_vect = (double *)mm_alloc(sizeof(double) * 3 * matrix.size);
	if(*d_vect == NULL) return(ALLOC_ERROR_RC);

	for(ii1=0; ii1<3*matrix.size; ii1++){
		if(init_vect == NULL){
			(*d_vect)[ii1] = 0.0;
		}else{
			(*d_vect)[ii1] = init_vect[ii1];
		}
	}

	RC_TRY( iccg_nonzero3s(matrix, f_vect, *d_vect, 1, 1) );

	return(NORMAL_RC);
}


RC
scale_iccg_nonzero3 (NONZERO_MATRIX3 matrix, double *f_vect, 
                     double **d_vect, const double *init_vect)
{
	int ii1;
	double *scale_fact;

	*d_vect = (double *)mm_alloc(sizeof(double) * 3 * matrix.size);
	scale_fact = (double *)mm_alloc(sizeof(double) * 3 * matrix.size);
	if((*d_vect == NULL)||(scale_fact == NULL)) return(ALLOC_ERROR_RC);

	for(ii1=0; ii1<3*matrix.size; ii1++){
		if(init_vect == NULL){
			(*d_vect)[ii1] = 0.0;
		}else{
			(*d_vect)[ii1] = init_vect[ii1];
		}
	}

	RC_TRY( nonzero3s_scaling_pre(matrix, scale_fact, f_vect) );
	RC_TRY( iccg_nonzero3s(matrix, f_vect, *d_vect, 1, 1) );
	RC_TRY( nonzero3s_scaling_post(matrix.size, scale_fact, *d_vect) );

	RC_TRY( mm_free(scale_fact) );

	return(NORMAL_RC);
}


RC
rest2vect (int dim, NODE_ARRAY node, BC_ARRAY rest, double **d_vect)
{
	int ii1;
	int index;

	if((dim != 2)&&(dim != 3)&&(dim != 6)) return(ARG_ERROR_RC);

	d_vect = mm_alloc((node.size) * dim * sizeof(double));
	if((d_vect== NULL)) return(ALLOC_ERROR_RC);

	for(ii1=0; ii1<dim*(node.size); ii1++){
		(*d_vect)[ii1] = 0.0;
	}

	for(ii1=0; ii1<rest.size; ii1++){
		if(rest.array[ii1].node < 0) continue;

		index = node_renum_index(node, rest.array[ii1].node);
		if(index < 0) return(SEEK_ERROR_RC);

		if(rest.array[ii1].v_type.x == BC_FIX){
			(*d_vect)[index*dim]     = rest.array[ii1].v.x;
		}
		if(rest.array[ii1].v_type.y == BC_FIX){
			(*d_vect)[index*dim + 1] = rest.array[ii1].v.y;
		}
		if(rest.array[ii1].v_type.z == BC_FIX){
			if(dim >= 3) (*d_vect)[index*dim + 2] = rest.array[ii1].v.z;
		}
		if(rest.array[ii1].v_type.yz == BC_FIX){
			if(dim == 6) (*d_vect)[index*dim + 3] = rest.array[ii1].v.yz;
		}
		if(rest.array[ii1].v_type.zx == BC_FIX){
			if(dim == 6) (*d_vect)[index*dim + 4] = rest.array[ii1].v.zx;
		}
		if(rest.array[ii1].v_type.xy == BC_FIX){
			if(dim == 6) (*d_vect)[index*dim + 5] = rest.array[ii1].v.xy;
		}
	}

	return(NORMAL_RC);
}


RC
disp2vect (int dim, NODE_ARRAY node, DISP_ARRAY disp, double **d_vect)
{
	int ii1;
	int index;

	if((dim != 2)&&(dim != 3)&&(dim != 6)) return(ARG_ERROR_RC);
	*d_vect = mm_alloc((node.size) * dim * sizeof(double));
	if((*d_vect== NULL)) return(ALLOC_ERROR_RC);

	for(ii1=0; ii1<dim*(node.size); ii1++){
		(*d_vect)[ii1] = 0.0;
	}

	for(ii1=0; ii1<disp.size; ii1++){
		if(disp.array[ii1].node < 0) continue;

		index = node_renum_index(node, disp.array[ii1].node);
		if(index < 0) return(SEEK_ERROR_RC);

		(*d_vect)[index*dim]     = disp.array[ii1].v.x;
		(*d_vect)[index*dim + 1] = disp.array[ii1].v.y;
		if(dim >= 3) (*d_vect)[index*dim + 2] = disp.array[ii1].v.z;
		if(dim == 6){
			(*d_vect)[index*dim + 3] = disp.array[ii1].v.yz;
			(*d_vect)[index*dim + 4] = disp.array[ii1].v.zx;
			(*d_vect)[index*dim + 5] = disp.array[ii1].v.xy;
		}
	}

	return(NORMAL_RC);
}


RC
vect2disp (int dim, NODE_ARRAY node, const double *d_vect, DISP_ARRAY *disp)
{
	int ii1;
	int index;

	if((dim != 2)&&(dim != 3)&&(dim != 6)) return(ARG_ERROR_RC);

	RC_TRY( allocate_disp_array(node.size, disp) );
	disp->source = FEM_K_SOL;
	for(ii1=0; ii1<node.size; ii1++){
		if(node.array[ii1].label < 0){
			disp->array[ii1].node = -1;
			continue;
		}
		index = node.array[ii1].renum;
		if((index < 0)||(index >= node.size)) return(ARG_ERROR_RC);

		disp->array[ii1].node = node.array[ii1].label;
		disp->array[ii1].v.x = d_vect[dim*index];
		disp->array[ii1].v.y = d_vect[dim*index + 1];
		if(dim >= 3) disp->array[ii1].v.z = d_vect[dim*index + 2];
		if(dim == 6){
			disp->array[ii1].v.yz = d_vect[dim*index + 3];
			disp->array[ii1].v.zx = d_vect[dim*index + 4];
			disp->array[ii1].v.xy = d_vect[dim*index + 5];
		}

	}
	disp->size = node.size;

	RC_TRY( clean_disp_array(disp) );

	return(NORMAL_RC);
}


RC
vect2react (int dim, NODE_ARRAY node, const double *r_vect,
            REACT_ARRAY *react)
{
	int ii1;
	int index;

	if((dim != 2)&&(dim != 3)) return(ARG_ERROR_RC);

	RC_TRY( allocate_react_array(node.size, react) );
	react->source = FEM_K_SOL;
	for(ii1=0; ii1<node.size; ii1++){
		if(node.array[ii1].label < 0){
			react->array[ii1].node = -1;
			continue;
		}
		index = node.array[ii1].renum;
		if((index < 0)||(index >= node.size)) return(ARG_ERROR_RC);

		react->array[ii1].node = node.array[ii1].label;
		react->array[ii1].v.x = r_vect[dim*index];
		react->array[ii1].v.y = r_vect[dim*index + 1];
		if(dim == 3) react->array[ii1].v.z = r_vect[dim*index + 2];
	}
	react->size = node.size;

	RC_TRY( clean_react_array(react) );

	return(NORMAL_RC);
}


RC
modify_nonzero3 (NODE_ARRAY node, BC_ARRAY rest,
                 NONZERO_MATRIX3 matrix, double *f_vect)
{
	return( modify_nonzero3_n(node, rest, matrix, &f_vect, 1) );
}


RC
modify_nonzero3_n (NODE_ARRAY node, BC_ARRAY rest, NONZERO_MATRIX3 matrix,
                   double **f_vect, int vect_size)
{
	int ii1;
	int index;

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	for(ii1=0; ii1<rest.size; ii1++){
		if(rest.array[ii1].node < 0) continue;

		RC_NEG_CHK( index = node_renum_index(node, rest.array[ii1].node) );

		if(rest.array[ii1].v_type.x == BC_FIX){
			RC_TRY( modify_nonzero3s_n(matrix, index, 0, rest.array[ii1].v.x,
			                           f_vect, vect_size) );
		}
		if(rest.array[ii1].v_type.y == BC_FIX){
			RC_TRY( modify_nonzero3s_n(matrix, index, 1, rest.array[ii1].v.y,
			                           f_vect, vect_size) );
		}
		if(rest.array[ii1].v_type.z == BC_FIX){
			RC_TRY( modify_nonzero3s_n(matrix, index, 2, rest.array[ii1].v.z,
			                           f_vect, vect_size) );
		}
	}

	return(NORMAL_RC);
}


/* 全体マトリックスの確保 */
RC
allocate_nonzero3 (NODE_ARRAY node, ELEMENT_ARRAY element,
                   NONZERO_MATRIX3 *matrix)
{
	int ii1, ii2;
	ELEM_LOOKUP_ARRAY elem_lookup;
	int *row_nodes;
	int row_size;
	int syn_size;
	int index;

	if(analysis_dim(element) != 3) return(ARG_ERROR_RC);
	if(node.renum_flag != RENUMBERED){
		RC_TRY( dummy_renumber(&node) );
	}

	/* 作業用配列 */
	row_nodes = (int *)mm_alloc(node.size*sizeof(int));
	if(row_nodes == NULL) return(ALLOC_ERROR_RC);

	RC_TRY( make_elem_lookup(node, element, &elem_lookup) );

	matrix->size = node.size;
	matrix->row = (NONZERO_ROW3 *)mm_alloc(matrix->size*sizeof(NONZERO_ROW3));
	if(matrix->row == NULL) return(ALLOC_ERROR_RC);

	for(ii1=0; ii1<matrix->size; ii1++){
		matrix->row[ii1].size = 0;
		matrix->row[ii1].elem = NULL;
		matrix->row[ii1].syn_size = 0;
		matrix->row[ii1].syn = NULL;
	}

	for(ii1=0; ii1<elem_lookup.size; ii1++){
		if(elem_lookup.array[ii1].node < 0) continue;

		index = node_renum_index(node, elem_lookup.array[ii1].node);
		if(index < 0) return(SEEK_ERROR_RC);

		RC_TRY( fill_row_nodes(elem_lookup.array[ii1], node, element,
		                       row_nodes, &row_size) );
		syn_size = 0;
		for(ii2=0; ii2<row_size; ii2++){
			if(row_nodes[ii2] > index){
				syn_size++;
			}
		}

		matrix->row[index].syn_size = syn_size;
		if(syn_size > 0){
			matrix->row[index].syn = mm_alloc(syn_size*sizeof(int));
			if(matrix->row[index].syn == NULL) return(ALLOC_ERROR_RC);
		}
		for(ii2=0; ii2<syn_size; ii2++){
			matrix->row[index].syn[ii2] = row_nodes[row_size - syn_size + ii2];
		}

		row_size -= syn_size;
		matrix->row[index].size = row_size;
		matrix->row[index].elem = mm_alloc(row_size*sizeof(NONZERO_ELEM3));
		if(matrix->row[index].elem == NULL) return(ALLOC_ERROR_RC);

		for(ii2=0; ii2<row_size; ii2++){
			matrix->row[index].elem[ii2].col = row_nodes[ii2];
			matrix->row[index].elem[ii2].v[0][0] = 0.0;
			matrix->row[index].elem[ii2].v[0][1] = 0.0;
			matrix->row[index].elem[ii2].v[0][2] = 0.0;
			matrix->row[index].elem[ii2].v[1][0] = 0.0;
			matrix->row[index].elem[ii2].v[1][1] = 0.0;
			matrix->row[index].elem[ii2].v[1][2] = 0.0;
			matrix->row[index].elem[ii2].v[2][0] = 0.0;
			matrix->row[index].elem[ii2].v[2][1] = 0.0;
			matrix->row[index].elem[ii2].v[2][2] = 0.0;
		}
	}

	RC_TRY( mm_free(row_nodes) );
	RC_TRY( free_elem_lookup_array(&elem_lookup) );

	return(NORMAL_RC);
}


/* 該当行に関連する列のインデックスを row_nodes に登録 */
static RC
fill_row_nodes (ELEM_LOOKUP elem_lookup,
                NODE_ARRAY node, ELEMENT_ARRAY element,
                int row_nodes[], int *row_size)
{
	int ii1, ii2, ii3;
	int flag, col_index;
	int elem_index;

	(*row_size) = 0;
	for(ii1=0; ii1<elem_lookup.size; ii1++){
		elem_index = search_element_label(element, elem_lookup.element[ii1]);
		if( elem_index < 0) return(SEEK_ERROR_RC);

		for(ii2=0; ii2<element.array[elem_index].node_num; ii2++){
			col_index = node_renum_index(node,
			                             element.array[elem_index].node[ii2]);
			if( col_index < 0) return(SEEK_ERROR_RC);

			flag = 0;
			for(ii3=0; ii3<(*row_size); ii3++){
				if(row_nodes[ii3] == col_index){
					flag = 1;
					break;
				}
			}
			if(flag == 1) continue;

			/* バブルソートしておく */
			row_nodes[*row_size] = col_index;
			for(ii3=(*row_size); ii3>=1; ii3--){
				if(row_nodes[ii3] < row_nodes[ii3-1]){
					int tmp = row_nodes[ii3];
					row_nodes[ii3] = row_nodes[ii3-1];
					row_nodes[ii3-1] = tmp;
				}else{
					break;
				}
			}
			(*row_size)++;
		}
	}

	return(NORMAL_RC);
}


/* 質量マトリックス */
RC
make_mass_matrix (ELEMENT element, NODE_ARRAY node,
                  MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
                  ELEM_MATRIX *elem_matrix)
{
	int ii1, ii2, ii3, ii4, ii5;
	static int init_flag = 0;
	static double **node_location = NULL;
	static double **gauss_points = NULL;
	static double weight[MAX_GAUSS_POINT];
	static double **dN = NULL;
	static double **jacobi = NULL;
	static double **rho_matrix = NULL;
	static int num_point;
	double rho;
	int mat_index;

	if(init_flag == 0){
		RC_TRY( allocate2D(MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(MAX_GAUSS_POINT, 4, &gauss_points) );
		RC_TRY( allocate2D(3, MAX_NODE, &dN) );
		RC_TRY( allocate2D(3, 3, &jacobi) );
		RC_TRY( allocate2D(3, 3, &rho_matrix) );

		init_flag = 1;
	}

	if(element.label < 0) return(ARG_ERROR_RC);
	elem_matrix->element = element.label;

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	RC_TRY( Gauss_const4mass_matrix(element.type, gauss_points,
	                                weight, &num_point) );
	RC_TRY( node_location_matrix(&element, node, node_location) );

	elem_matrix->dim = element_dim(element.type);
	if(elem_matrix->dim <= 0) return(ARG_ERROR_RC);

	elem_matrix->size = (element.node_num) * (elem_matrix->dim);

	/* elem_matrix->array[][] の確保 */
	RC_TRY( allocate2D(elem_matrix->size,
	                   elem_matrix->size, &(elem_matrix->matrix)) );

	/* elem_matirx->index[] の確保と代入 */
	elem_matrix->index = mm_alloc(sizeof(int)*(elem_matrix->size) );
	if(elem_matrix->index == NULL) return(ALLOC_ERROR_RC);

	for(ii1=0; ii1<element.node_num; ii1++){
		int index = node_renum_index(node, element.node[ii1]);
		if(index < 0) return(SEEK_ERROR_RC);

		for(ii2=0; ii2<(elem_matrix->dim); ii2++){
			elem_matrix->index[ii1*(elem_matrix->dim)+ii2]
			             = index*(elem_matrix->dim) + ii2;
		}
	}

	/* rho, rho_matrix */
	mat_index = search_material_prop_element(material, physical, &element);
	if(mat_index < 0) return(SEEK_ERROR_RC);

	rho = material.array[mat_index].rho;
	for(ii1=0; ii1<(elem_matrix->dim); ii1++){
		for(ii2=0; ii2<(elem_matrix->dim); ii2++){
			rho_matrix[ii1][ii2] = 0.0;
		}
		rho_matrix[ii1][ii1] = rho;
	}

	for(ii1=0; ii1<num_point; ii1++){
		double N[MAX_NODE];
		double det_J = 0.0;

		/* det_J, N を計算 */
		RC_TRY( set_dN_matrix(element.type, gauss_points[ii1], dN) );
		RC_TRY( set_N_vector(element.type, gauss_points[ii1], N) );
		mul_matrix((elem_matrix->dim), element.node_num,
		           (elem_matrix->dim), dN, node_location, jacobi);
		det_J = determinant((elem_matrix->dim), jacobi);

		for(ii2=0; ii2<element.node_num; ii2++){
			for(ii3=0; ii3<element.node_num; ii3++){
				for(ii4=0; ii4<(elem_matrix->dim); ii4++){
					int index_24 = ii2*(elem_matrix->dim) + ii4;

					for(ii5=0; ii5<(elem_matrix->dim); ii5++){
						int index_35 = ii3*(elem_matrix->dim) + ii5;

						/* elem_matrix->matrix[index_24][index_35]
						     += rho * weight[ii1] * det_J * N[ii2] * N[ii3]; */
						elem_matrix->matrix[index_24][index_35]
						     += rho_matrix[ii4][ii5] * weight[ii1]
						      * det_J * N[ii2] * N[ii3];
					}
				}
			}
		}
	}

	return(NORMAL_RC);
}

/* 質量マトリックス */
RC
make_mass_matrix_add_density (ELEMENT element, NODE_ARRAY node,
                  MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
                  ELEM_MATRIX *elem_matrix,double *density)
{
	int ii1, ii2, ii3, ii4, ii5;
	static int init_flag = 0;
	static double **node_location = NULL;
	static double **gauss_points = NULL;
	static double weight[MAX_GAUSS_POINT];
	static double **dN = NULL;
	static double **jacobi = NULL;
	static double **rho_matrix = NULL;
	static int num_point;
	double rho;
	double phi;
	int mat_index;
	int den_index;

	if(init_flag == 0){
		RC_TRY( allocate2D(MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(MAX_GAUSS_POINT, 4, &gauss_points) );
		RC_TRY( allocate2D(3, MAX_NODE, &dN) );
		RC_TRY( allocate2D(3, 3, &jacobi) );
		RC_TRY( allocate2D(3, 3, &rho_matrix) );

		init_flag = 1;
	}

	if(element.label < 0) return(ARG_ERROR_RC);
	elem_matrix->element = element.label;

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	RC_TRY( Gauss_const4mass_matrix(element.type, gauss_points,
	                                weight, &num_point) );
	RC_TRY( node_location_matrix(&element, node, node_location) );

	elem_matrix->dim = element_dim(element.type);
	if(elem_matrix->dim <= 0) return(ARG_ERROR_RC);

	elem_matrix->size = (element.node_num) * (elem_matrix->dim);

	/* elem_matrix->array[][] の確保 */
	RC_TRY( allocate2D(elem_matrix->size,
	                   elem_matrix->size, &(elem_matrix->matrix)) );

	/* elem_matirx->index[] の確保と代入 */
	elem_matrix->index = mm_alloc(sizeof(int)*(elem_matrix->size) );
	if(elem_matrix->index == NULL) return(ALLOC_ERROR_RC);

	for(ii1=0; ii1<element.node_num; ii1++){
		int index = node_renum_index(node, element.node[ii1]);
		if(index < 0) return(SEEK_ERROR_RC);

		for(ii2=0; ii2<(elem_matrix->dim); ii2++){
			elem_matrix->index[ii1*(elem_matrix->dim)+ii2]
			             = index*(elem_matrix->dim) + ii2;
		}
	}



	/* rho, rho_matrix */
	mat_index = search_material_prop_element(material, physical, &element);
	if(mat_index < 0) return(SEEK_ERROR_RC);

	rho = material.array[mat_index].rho;
	for(ii1=0; ii1<(elem_matrix->dim); ii1++){
		for(ii2=0; ii2<(elem_matrix->dim); ii2++){
			rho_matrix[ii1][ii2] = 0.0;
		}
		rho_matrix[ii1][ii1] = rho;
	}

	for(ii1=0; ii1<num_point; ii1++){
		double N[MAX_NODE];
		double det_J = 0.0;

		/* det_J, N を計算 */
		RC_TRY( set_dN_matrix(element.type, gauss_points[ii1], dN) );
		RC_TRY( set_N_vector(element.type, gauss_points[ii1], N) );
		mul_matrix((elem_matrix->dim), element.node_num,
		           (elem_matrix->dim), dN, node_location, jacobi);
		det_J = determinant((elem_matrix->dim), jacobi);

		phi = 0;
		for(ii2=0; ii2<element.node_num; ii2++){
			den_index = search_node_label(node,element.node[ii2]);
			phi += density[den_index]* N[ii2];
		}

		//printf("%lf\n",phi);


		for(ii2=0; ii2<element.node_num; ii2++){
			for(ii3=0; ii3<element.node_num; ii3++){
				for(ii4=0; ii4<(elem_matrix->dim); ii4++){
					int index_24 = ii2*(elem_matrix->dim) + ii4;

					for(ii5=0; ii5<(elem_matrix->dim); ii5++){
						int index_35 = ii3*(elem_matrix->dim) + ii5;

						/* elem_matrix->matrix[index_24][index_35]
						     += rho * weight[ii1] * det_J * N[ii2] * N[ii3]; */
						elem_matrix->matrix[index_24][index_35]
						     += rho_matrix[ii4][ii5] * weight[ii1]
						      * det_J * N[ii2] * N[ii3] * phi;
					}
				}
			}
		}
	}

	return(NORMAL_RC);
}


static RC
Gauss_const4mass_matrix (ELEM_TYPE type, double **points,
                         double *weight, int *num_point)
{
	switch(type){
	case ELEM_RBAR:
	case ELEM_BEAM1:
	case ELEM_BEAM2:
	case ELEM_LINE1:
	case ELEM_LINE2:
		*num_point = 3;
		RC_TRY( Gauss_const_line(*num_point, points, weight) );
		break;
	case ELEM_TRI1:
	case ELEM_TRI2:
		*num_point = 3;
		RC_TRY( Gauss_const_tri(*num_point, points, weight) );
		break;
	case ELEM_QUAD1:
		*num_point = 4;
		RC_TRY( Gauss_const_quad(*num_point, points, weight) );
		break;
	case ELEM_QUAD2:
	case ELEM_QUAD2L:
		*num_point = 9;
		RC_TRY( Gauss_const_quad(*num_point, points, weight) );
		break;
	case ELEM_TETRA1:
	case ELEM_TETRA2:
		*num_point = 4;
		RC_TRY( Gauss_const_tetra(*num_point, points, weight) );
		break;
	case ELEM_PENTA1:
		*num_point = 6;
		RC_TRY( Gauss_const_penta(*num_point, points, weight) );
		break;
	case ELEM_PENTA2:
		*num_point = 9;
		RC_TRY( Gauss_const_penta(*num_point, points, weight) );
		break;
	case ELEM_HEXA1:
		*num_point = 8;
		RC_TRY( Gauss_const_hexa(*num_point, points, weight) );
		break;
	case ELEM_HEXA2:
		*num_point = 14;
		RC_TRY( Gauss_const_hexa(*num_point, points, weight) );
		break;
	default :
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


/*
 * 座標変換マトリクスL[3][3]
 * 各要素に対して局所座標系を定義
 */
RC
set_L_matrix (ELEMENT element, double **node_location, double L[][3])
{
	switch(element.type){
	case ELEM_BEAM1:
	case ELEM_BEAM2:
		RC_TRY( set_L_matrix_beam(element, node_location, L) );
		break;
	case ELEM_TRI1:
	case ELEM_TRI2:
	case ELEM_QUAD1:
	case ELEM_QUAD2:
	case ELEM_HEXA1:
	case ELEM_HEXA2:
	case ELEM_PENTA1:
	case ELEM_PENTA2:
		RC_TRY( set_L_matrix_shell_solid(element, node_location, L) );
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


/*
 * シェル要素，ソリッド要素 座標変換マトリクスL[3][3]
 * 各要素に対して局所座標系を定義
 */
RC
set_L_matrix_shell_solid (ELEMENT element, double **node_location,
                          double L[][3])
{
	VECT3D vtmp1, vtmp2, vtmp3, vtmp4, vtmp5, vtmp6, vtmp7;
	VECT3D ex, ey, ez;


	switch(element.type){
	case ELEM_TRI1:
	case ELEM_TRI2:
		vtmp1.x = node_location[1][0] - node_location[0][0];
		vtmp1.y = node_location[1][1] - node_location[0][1];
		vtmp1.z = node_location[1][2] - node_location[0][2];
		if(nearly_eq(abs_vect3d(vtmp1), 0.0)) return(CAL_ERROR_RC);
		ex = unit_vect3d(vtmp1);

		vtmp2.x = node_location[2][0] - node_location[0][0];
		vtmp2.y = node_location[2][1] - node_location[0][1];
		vtmp2.z = node_location[2][2] - node_location[0][2];
		vtmp2 = gs_ortho3d(vtmp2, ex);
		if(nearly_eq(abs_vect3d(vtmp2), 0.0)) return(CAL_ERROR_RC);
		ey = unit_vect3d(vtmp2);

		ez = unit_vect3d( outer_product3d(ex, ey) );
		break;
	case ELEM_QUAD1:
	case ELEM_QUAD2:
		vtmp1.x = node_location[1][0] - node_location[3][0];
		vtmp1.y = node_location[1][1] - node_location[3][1];
		vtmp1.z = node_location[1][2] - node_location[3][2];

		vtmp2.x = node_location[2][0] - node_location[0][0];
		vtmp2.y = node_location[2][1] - node_location[0][1];
		vtmp2.z = node_location[2][2] - node_location[0][2];

		if(nearly_eq(abs_vect3d(vtmp1), 0.0)) return(CAL_ERROR_RC);
		if(nearly_eq(abs_vect3d(vtmp2), 0.0)) return(CAL_ERROR_RC);

		ex = unit_vect3d( add_vect3d(unit_vect3d(vtmp2),
		                             unit_vect3d(vtmp1)) );
		ey = unit_vect3d( sub_vect3d(unit_vect3d(vtmp2),
		                             unit_vect3d(vtmp1)) );
		ey = unit_vect3d( gs_ortho3d(ey, ex) );
		ez = unit_vect3d( outer_product3d(ex, ey) );
		break;
	case ELEM_HEXA1:
	case ELEM_HEXA2:
		vtmp1 = sub_vect3d(center3d_p4array(node_location[2],node_location[1],
		                                    node_location[5],node_location[6]),
		                   center3d_p4array(node_location[3],node_location[0],
		                                    node_location[4],node_location[7]));
		vtmp2 = sub_vect3d(center3d_p4array(node_location[3],node_location[2],
		                                    node_location[6],node_location[7]),
		                   center3d_p4array(node_location[0],node_location[1],
		                                    node_location[5],node_location[4]));
		vtmp3 = sub_vect3d(center3d_p4array(node_location[4],node_location[5],
		                                    node_location[6],node_location[7]),
		                   center3d_p4array(node_location[0],node_location[1],
		                                    node_location[2],node_location[3]));

		vtmp4 = gs_ortho3d(vtmp1, vtmp3);
		vtmp5 = gs_ortho3d(vtmp2, vtmp3);
		vtmp6 = unit_vect3d(vtmp3);
		vtmp7 = add_vect3d(unit_vect3d(vtmp4), unit_vect3d(vtmp5));
		vtmp4 = gs_ortho3d(vtmp4, vtmp7);
		vtmp5 = gs_ortho3d(vtmp5, vtmp7);
		vtmp4 = unit_vect3d( add_vect3d( unit_vect3d(vtmp4),
		                                 unit_vect3d(vtmp7) ) );
		vtmp5 = unit_vect3d( add_vect3d( unit_vect3d(vtmp5),
		                                 unit_vect3d(vtmp7) ) );
		ex = vtmp4;
		ey = vtmp5;
		ez = vtmp6;

		vtmp4 = gs_ortho3d(vtmp1, vtmp2);
		vtmp5 = unit_vect3d(vtmp2);
		vtmp6 = gs_ortho3d(vtmp3, vtmp2);
		vtmp7 = add_vect3d(unit_vect3d(vtmp4), unit_vect3d(vtmp6));
		vtmp4 = gs_ortho3d(vtmp4, vtmp7);
		vtmp6 = gs_ortho3d(vtmp6, vtmp7);
		vtmp4 = unit_vect3d( add_vect3d( unit_vect3d(vtmp4),
		                                 unit_vect3d(vtmp7) ) );
		vtmp6 = unit_vect3d( add_vect3d( unit_vect3d(vtmp6),
		                                 unit_vect3d(vtmp7) ) );
		ex = add_vect3d(ex, vtmp4);
		ey = add_vect3d(ey, vtmp5);
		ez = add_vect3d(ez, vtmp6);

		vtmp4 = unit_vect3d(vtmp1);
		vtmp5 = gs_ortho3d(vtmp2, vtmp1);
		vtmp6 = gs_ortho3d(vtmp3, vtmp1);
		vtmp7 = add_vect3d(unit_vect3d(vtmp5), unit_vect3d(vtmp6));
		vtmp5 = gs_ortho3d(vtmp5, vtmp7);
		vtmp6 = gs_ortho3d(vtmp6, vtmp7);
		vtmp5 = unit_vect3d( add_vect3d( unit_vect3d(vtmp5),
		                                 unit_vect3d(vtmp7) ) );
		vtmp6 = unit_vect3d( add_vect3d( unit_vect3d(vtmp6),
		                                 unit_vect3d(vtmp7) ) );

		ex = add_vect3d(ex, vtmp4);
		ey = add_vect3d(ey, vtmp5);
		ez = add_vect3d(ez, vtmp6);

		ex = unit_vect3d(ex);
		ey = unit_vect3d(gs_ortho3d(ey, ex));
		ez = unit_vect3d(gs_ortho3d(ez, ex));
		ez = unit_vect3d(gs_ortho3d(ez, ey));
		break;
	case ELEM_PENTA1:
	case ELEM_PENTA2:
		/* ez:上下面の重心間ベクトル */
		vtmp1 = sub_vect3d(center3d_p3array(node_location[3],node_location[4],
		                                    node_location[5]),
		                   center3d_p3array(node_location[0],node_location[1],
		                                    node_location[2]));
		ez = unit_vect3d(vtmp1);

		vtmp2.x = node_location[1][0] - node_location[0][0];
		vtmp2.y = node_location[1][1] - node_location[0][1];
		vtmp2.z = node_location[1][2] - node_location[0][2];
		vtmp3.x = node_location[4][0] - node_location[3][0];
		vtmp3.y = node_location[4][1] - node_location[3][1];
		vtmp3.z = node_location[4][2] - node_location[3][2];
		vtmp4 = mid_point3d(vtmp2, vtmp3);

		ex = unit_vect3d( gs_ortho3d(vtmp4, ez) );

		ey = unit_vect3d( outer_product3d(ex, ez) );
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	L[0][0] = ex.x;
	L[0][1] = ex.y;
	L[0][2] = ex.z;
	L[1][0] = ey.x;
	L[1][1] = ey.y;
	L[1][2] = ey.z;
	L[2][0] = ez.x;
	L[2][1] = ez.y;
	L[2][2] = ez.z;

	return(NORMAL_RC);
}


/*
 * ビーム要素 座標変換マトリクスL[3][3]
 * 各要素に対して局所座標系を定義
 */
RC
set_L_matrix_beam (ELEMENT element, double **node_location, double L[][3])
{
	VECT3D vtmp1, vtmp2;
	VECT3D ex, ey, ez;


	switch(element.type){
	case ELEM_BEAM1:
	case ELEM_BEAM2:
		vtmp1.x = node_location[1][0] - node_location[0][0];
		vtmp1.y = node_location[1][1] - node_location[0][1];
		vtmp1.z = node_location[1][2] - node_location[0][2];
		if(nearly_eq(abs_vect3d(vtmp1), 0.0)) return(CAL_ERROR_RC);
		ex = unit_vect3d(vtmp1);

		vtmp2.x = element.d_info[0];
		vtmp2.y = element.d_info[1];
		vtmp2.z = element.d_info[2];
		vtmp2 = gs_ortho3d(vtmp2, ex);
		if(nearly_eq(abs_vect3d(vtmp2), 0.0)) return(CAL_ERROR_RC);
		ey = unit_vect3d(vtmp2);

		ez = unit_vect3d( outer_product3d(ex, ey) );
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	L[0][0] = ex.x;
	L[0][1] = ex.y;
	L[0][2] = ex.z;
	L[1][0] = ey.x;
	L[1][1] = ey.y;
	L[1][2] = ey.z;
	L[2][0] = ez.x;
	L[2][1] = ez.y;
	L[2][2] = ez.z;

	return(NORMAL_RC);
}


RC
fem_solve1 (NODE_ARRAY node, ELEMENT_ARRAY element,
           MATERIAL_PROP_ARRAY material,
           PHYSICAL_PROP_ARRAY physical,
           BC_ARRAY rest, BC_ARRAY force,
           BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,double basis_temp,
           DISP_ARRAY *disp, STRESS_ARRAY *stress,
           STRAIN_ARRAY *strain, int sol_ss_flag)
{
	SKYLINE_MATRIX K_matrix;
	double *f_vect;

	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );
	RC_TRY( make_g_matrix(node, element, material, physical, &K_matrix) );
	RC_TRY( fem_allocate_vector(K_matrix.dim, node, &f_vect) );
	RC_TRY( add_nodal_force(K_matrix.dim, node, force, f_vect) );
	RC_TRY( add_temp_force1(K_matrix.dim, node, element, material, physical,
	                       temp, def_temp, basis_temp, f_vect) );
	RC_TRY( modify_g_matrix(node, rest, K_matrix, f_vect) );
	RC_TRY( decomp_g_matrix(K_matrix) );
	RC_TRY( subst_g_matrix(node, K_matrix, f_vect, disp) );
	RC_TRY( free_g_matrix(&K_matrix) );
	RC_TRY( stress_strain(sol_ss_flag, node, element, *disp,
	                      material, physical, stress, strain) );
	RC_TRY( mm_free(f_vect) );

	return(NORMAL_RC);
}


/* 熱ひずみによる荷重を加算 */
RC
add_temp_force1 (int dim, NODE_ARRAY node, ELEMENT_ARRAY element,
                MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
                BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp, double basis_temp,
                double *f_vect)
{
	double temp_load[3*MAX_NODE];
	int index_array[3*MAX_NODE];
	int vect_size;
	int ii1, ii2;

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;
		if(element.array[ii1].type == ELEM_BEAM1) continue;

		RC_TRY( element_temp_load1(dim, element.array[ii1], node, material,
		                           physical, temp, def_temp, basis_temp, 
		                           index_array, &vect_size, temp_load) );
		for(ii2=0; ii2<vect_size; ii2++){
			f_vect[index_array[ii2]] += temp_load[ii2];
		}
	}

	return(NORMAL_RC);
}


/* 要素の熱ひずみによる荷重 基準温度の導入*/
static RC
element_temp_load1 (int gdim, ELEMENT element, NODE_ARRAY node,
                    MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
                    BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp, double basis_temp,
                    int *index_array, int *vect_size, double *temp_load)
{
	int dim;
	int ii1, ii2;
	double thickness;
	int m_index;
	double node_temp_array[MAX_NODE];
	static int init_flag = 0;
	static double **D_matrix = NULL;
	static double **node_location = NULL;
	static double **trans_B = NULL;
	static double **BtD = NULL;
	static G_PROP g_prop;
	MATERIAL_PROP mat;

	if(init_flag == 0){
		RC_TRY( allocate2D(6, 6, &D_matrix) );
		RC_TRY( allocate2D(MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(3*MAX_NODE, 6, &trans_B) );
		RC_TRY( allocate2D(3*MAX_NODE, 6, &BtD) );
		RC_TRY( allocate_g_prop(&g_prop) );
		init_flag = 1;
	}

	if((dim = element_dim(element.type)) <= 1) return(ARG_ERROR_RC);

	if(gdim < dim) return(ARG_ERROR_RC);
	*vect_size = dim*element.node_num;

	m_index = search_material_prop_element(material, physical, &element);
	if(m_index < 0) return(SEEK_ERROR_RC);

	mat = material.array[m_index];
	RC_TRY( fill_material(&mat) );

	/* temp_load を初期化       */
	/* index_array を初期化     */
	/* node_temp_array を初期化 */
	for(ii1=0; ii1<element.node_num; ii1++){
		int index = node_renum_index(node, element.node[ii1]);
		if(index < 0) return(SEEK_ERROR_RC);

		for(ii2=0; ii2<dim; ii2++){
			temp_load[dim*ii1 + ii2] = 0.0;
			index_array[dim*ii1 + ii2] = gdim*index + ii2;
		}
		node_temp_array[ii1] = node_temp(temp, def_temp, element.node[ii1]);
	}

	/* 要素のD マトリックス */
	RC_TRY( get_D_matrix(element, material, physical, D_matrix) );
	RC_TRY( get_thickness(element, physical, &thickness) );

	/* ガウスポイントごとの B マトリックス */
	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( set_Gauss_const(element.type, &g_prop) );
	RC_TRY( set_g_prop(element, node_location, &g_prop) );

	for(ii1=0; ii1<g_prop.num_g_point; ii1++){
		double local_temp;
		double local_force[3*MAX_NODE];
		double N[MAX_NODE];
		double init_strain[6];
		int ii2;

		/* ガウスポイントの温度を内装して求める */
		RC_TRY( set_N_vector(element.type, g_prop.g_points[ii1], N) );
		local_temp = 0.0;
		for(ii2=0; ii2<element.node_num; ii2++){
			local_temp += (node_temp_array[ii2] - basis_temp) * N[ii2];
		}

		/* init_strain の初期化 */
		if(dim == 2){
			init_strain[0] = local_temp*mat.alpha_vect.x;
			init_strain[1] = local_temp*mat.alpha_vect.y;
			init_strain[2] = local_temp*mat.alpha_vect.xy;
		}else{   /* dim == 3 */
			init_strain[0] = local_temp*mat.alpha_vect.x;
			init_strain[1] = local_temp*mat.alpha_vect.y;
			init_strain[2] = local_temp*mat.alpha_vect.z;
			init_strain[3] = local_temp*mat.alpha_vect.yz;
			init_strain[4] = local_temp*mat.alpha_vect.zx;
			init_strain[5] = local_temp*mat.alpha_vect.xy;
		}

		/* B^T => trans_B */
		/* trans_B * D => BD */
		/* BD * init_strain => local_force */
		transpose_matrix(g_prop.B_matrix_row, g_prop.B_matrix_col,
		                 g_prop.B_matrix[ii1], trans_B);
		mul_matrix(g_prop.B_matrix_col, g_prop.B_matrix_row,
		           g_prop.B_matrix_row, trans_B, D_matrix, BtD);
		mul_matrix_vect(g_prop.B_matrix_col, g_prop.B_matrix_row, 
		                BtD, init_strain, local_force);

		/* local_force * detJ * weight * thickness */
		for(ii2=0; ii2<(*vect_size); ii2++){
			temp_load[ii2] += local_force[ii2] * g_prop.det_J[ii1]
			                * g_prop.weight[ii1] * thickness;
		}
	}

	return(NORMAL_RC);
}


/* 熱伝導解析を行う */
RC
fem_solve_thermal (NODE_ARRAY node, ELEMENT_ARRAY element,
                   ELEMENT_ARRAY surf, DEFAULT_BC_ARRAY def_temp,
                   MATERIAL_PROP_ARRAY material,
                   PHYSICAL_PROP_ARRAY physical,
                   BC_ARRAY bc_temp, BC_ARRAY flux,
                   ELEMENT_ARRAY transfer, BC_ARRAY *result_temp)
{
	SKYLINE_MATRIX g_matrix1;
	double *s_vect;

	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

	/* Kマトリックスの初期化 */
	RC_TRY( allocate_g_matrix(1, node, element, &g_matrix1) );

	RC_TRY( make_K_matrix_thermal(g_matrix1, node, element, surf, transfer,
	                              material, physical) );

	/* 右辺ベクトルの初期化 */
	RC_TRY( fem_allocate_vector(1, node, &s_vect) );
	RC_TRY( add_therm_vector(node, flux, s_vect) );
	RC_TRY( add_transfer_vector(surf, node, transfer, s_vect) );

	/* マトリックスを解く */
	RC_TRY( modify_g_matrix_therm(node, bc_temp, g_matrix1, s_vect) );
	RC_TRY( decomp_g_matrix(g_matrix1) );
	RC_TRY( subst_g_matrix_therm(node, g_matrix1, s_vect, result_temp) );

	/* 確保した分を解法 */
	RC_TRY( free_g_matrix(&g_matrix1) );
	RC_TRY( fem_free_vector(&s_vect) );

	return(NORMAL_RC);
}


/* 熱伝導場に関する全体剛性マトリックスの作成 */
RC
make_K_matrix_thermal (SKYLINE_MATRIX g_matrix1,
                       NODE_ARRAY node, ELEMENT_ARRAY element,
                       ELEMENT_ARRAY surf, ELEMENT_ARRAY transfer,
                       MATERIAL_PROP_ARRAY material,
                       PHYSICAL_PROP_ARRAY physical)
{
	int ii1;
	ELEM_MATRIX elem_matrix;

	/* 熱伝導に関するマトリックスの作成 */
	for(ii1=0;ii1<element.size;ii1++){
		if(element.array[ii1].label < 0) continue;

		RC_TRY( make_elem_matrix_therm(element.array[ii1], node, material,
		                               physical, &elem_matrix) );

		RC_TRY( impose_elem_matrix(g_matrix1, elem_matrix, 1.0) );
		RC_TRY( free_elem_matrix(&elem_matrix) );
	}

	/* 熱伝達に関するマトリックスの作成 */
	for(ii1=0; ii1 < transfer.size; ii1++){
		if(transfer.array[ii1].label < 0) continue;

		RC_TRY( make_elem_matrix_transfer(transfer.array[ii1], node, surf,
		                                  &elem_matrix) );
		RC_TRY( impose_elem_matrix(g_matrix1, elem_matrix, 1.0) );
		RC_TRY( free_elem_matrix(&elem_matrix) );
	}
	return(NORMAL_RC);
}


/* 熱伝達に関する要素マトリックスの作成 */
RC
make_elem_matrix_transfer (ELEMENT transfer, NODE_ARRAY node,
                           ELEMENT_ARRAY surf, ELEM_MATRIX *elem_matrix)
{
	int ii1,ii2,ii3;
	WORK_SET ws;
	ELEMENT surf_elem;
	int index = 0;
	int num_points;
	double *weights;
	double **gauss_points;
	double **node_location;

	RC_TRY(allocate_work_set(&ws));

	index = search_element_label(surf, transfer.label);
	if(index < 0) return(SEEK_ERROR_RC);
	surf_elem = surf.array[index];

	elem_matrix->element = surf_elem.label;
	elem_matrix->dim = 1;
	elem_matrix->size = surf_elem.node_num;

	/* elem_matrix->array[][] の確保 */
	RC_TRY( allocate2D(elem_matrix->size,
	                   elem_matrix->size, &(elem_matrix->matrix)) );

	/* elem_matirx->index[] の確保と代入 */
	elem_matrix->index = mm_alloc(sizeof(int) * (elem_matrix->size) );
	if(elem_matrix->index == NULL) return(ALLOC_ERROR_RC);

	for(ii1=0; ii1 < surf_elem.node_num; ii1++){
		int index1 = node_renum_index(node, surf_elem.node[ii1]);
		if(index1 < 0) return(SEEK_ERROR_RC);

		elem_matrix->index[ii1] = index1;
	}

	weights = ws.vec_G[0];
	gauss_points = ws.mat_G_4[0];
	node_location = ws.mat_N_3[0];

	RC_TRY( Gauss_const_default(surf_elem.type, gauss_points, weights,
	                            &num_points) );
	RC_TRY( node_location_matrix(&(surf_elem), node, node_location) );

	for(ii1=0; ii1 < num_points; ii1++){
		double *N = ws.vec_N[0];
		VECT3D d_A;
		double det_J;

		/* 積分の準備 */
		RC_TRY( set_N_vector(surf_elem.type, gauss_points[ii1], N) );
		RC_TRY( cal_d_A(surf_elem, gauss_points[ii1], node_location,
		                &d_A, ws.mat_3_N[0], ws.mat_3_3[0]) );

		det_J = sqrt(d_A.x*d_A.x + d_A.y*d_A.y + d_A.z*d_A.z);
		for(ii2=0; ii2 < surf_elem.node_num; ii2++){
			for(ii3=0; ii3 < surf_elem.node_num; ii3++){
				elem_matrix->matrix[ii2][ii3] += transfer.d_info[0] * N[ii2]
				                               * N[ii3] * weights[ii1] * det_J;
			}
		}
	}

	RC_TRY(free_work_set(&ws));
	
	return(NORMAL_RC);
}


/* 熱流束ベクトルを計算 */
RC
add_therm_vector (NODE_ARRAY node, BC_ARRAY flux, double *f_vect)
{
	int index, ii1;

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);
	
	for(ii1=0; ii1<flux.size; ii1++){
		if(flux.array[ii1].node < 0) continue;

		index = node_renum_index(node, flux.array[ii1].node);

		if(index < 0){
			fprintf(stderr, "node_renum_index() < ");
			return(SEEK_ERROR_RC);
		}

		f_vect[index] += flux.array[ii1].s[0];
	}
	return(NORMAL_RC);
}


/* 熱伝達ベクトルを計算 */
RC
add_transfer_vector (ELEMENT_ARRAY surf, NODE_ARRAY node,
                     ELEMENT_ARRAY transfer, double *s_vect)
{
	int ii1,ii2,ii3,index;
	WORK_SET ws;
	double *weights;
	double **gauss_points;
	double **node_location;
	int num_points;
	ELEMENT surf_elem;
	
	if(transfer.size == 0) return(NORMAL_RC);
	
	RC_TRY(allocate_work_set(&ws));
	
	weights = ws.vec_G[0];
	gauss_points = ws.mat_G_4[0];
	node_location = ws.mat_N_3[0];
	
	for(ii1=0;ii1<transfer.size;ii1++){
		if(transfer.array[ii1].label < 0) continue;
		
		index = search_element_label(surf, transfer.array[ii1].label);

		if(index < 0) return(SEEK_ERROR_RC);

		surf_elem = surf.array[index];

		RC_TRY( Gauss_const_default(surf_elem.type, gauss_points, weights,
		                            &num_points) );
		RC_TRY( node_location_matrix(&(surf_elem), node, node_location) );

		for(ii2=0; ii2<num_points; ii2++){
			double *N = ws.vec_N[0];
			VECT3D d_A;
			double det_J;

			RC_TRY( set_N_vector(surf_elem.type, gauss_points[ii2], N) );
			RC_TRY( cal_d_A(surf_elem, gauss_points[ii2], node_location,
			                &d_A, ws.mat_3_N[0], ws.mat_3_3[0]) );
			
			det_J = sqrt(d_A.x*d_A.x + d_A.y*d_A.y + d_A.z*d_A.z);
			for(ii3=0;ii3<surf_elem.node_num;ii3++){
				index = node_renum_index(node, surf_elem.node[ii3]);
				
				if(index < 0) return(SEEK_ERROR_RC);
				
				s_vect[index] += transfer.array[ii1].d_info[0]
				               * transfer.array[ii1].d_info[1] * weights[ii2]
				               * N[ii3] * det_J;
			}
		}
	}
	RC_TRY(free_work_set(&ws));

	return(NORMAL_RC);
}


/* 非定常熱伝導解析を行う 時間方向には差分法を用いている
   step_numがステップ数、delta_tがステップ間隔
   thetaが0.5の時クランクニコルソン法になる */
RC
fem_solve_thermal_transient (NODE_ARRAY node, ELEMENT_ARRAY element,
                             ELEMENT_ARRAY surf, DEFAULT_BC_ARRAY def_temp,
                             MATERIAL_PROP_ARRAY material,
                             PHYSICAL_PROP_ARRAY physical,
                             BC_ARRAY bc_temp, BC_ARRAY flux,
                             ELEMENT_ARRAY transfer, BC_ARRAY **result_temp,
                             int step_num, double delta_t, double theta)
{
	int ii1, ii2;
	double *s_vect;
	BC_ARRAY *tmp;
	SKYLINE_MATRIX g_matrix1,g_matrix2,g_matrix3; /* 左辺用、右辺用、コピー用 */

	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

	/* Kマトリックスの初期化 */
	RC_TRY( allocate_g_matrix(1, node, element, &g_matrix1) );
	RC_TRY( allocate_g_matrix(1, node, element, &g_matrix2) );

	tmp = mm_alloc(sizeof(BC_ARRAY)*(step_num + 1));
	if(tmp == NULL) return(ALLOC_ERROR_RC);
	*result_temp = tmp;

	RC_TRY( init_cond_temp(&tmp[0], def_temp, bc_temp, node) );

	/* Kマトリックスの作成 */
	RC_TRY( make_K_matrix_transient(g_matrix1, g_matrix2, node, element, surf,
	                                transfer, material, delta_t, physical,
	                                theta) );

	for(ii1=0; ii1<step_num; ii1++){
		/* 右辺ベクトルの初期化、作成 */
		RC_TRY( fem_allocate_vector(1, node, &s_vect) );
		RC_TRY( add_transient_therm_vector(s_vect, g_matrix2, node, tmp[ii1]) );
		RC_TRY( add_therm_vector(node, flux, s_vect) );
		RC_TRY( add_transfer_vector(surf, node, transfer, s_vect) );

		/* マトリックスを解く */
		RC_TRY( allocate_g_matrix(1, node, element, &g_matrix3) );
		/* Kマトリックスのコピー */
		for(ii2 = 0; ii2 < g_matrix1.array_size; ii2++){
			g_matrix3.array[ii2] = g_matrix1.array[ii2];
		}
		RC_TRY( modify_g_matrix_therm(node, bc_temp, g_matrix3, s_vect) );
		RC_TRY( decomp_g_matrix(g_matrix3) );
		RC_TRY( subst_g_matrix_therm(node, g_matrix3, s_vect, &tmp[ii1 + 1]) );

		RC_TRY( free_g_matrix(&g_matrix3) );
		RC_TRY( fem_free_vector(&s_vect) );
	}

	/* 確保した分を解法 */
	RC_TRY( free_g_matrix(&g_matrix1) );
	RC_TRY( free_g_matrix(&g_matrix2) );

	return(NORMAL_RC);
}


/* 初期時刻での温度分布 */
RC
init_cond_temp (BC_ARRAY *init_temp, DEFAULT_BC_ARRAY def_temp,
                BC_ARRAY bc_temp, NODE_ARRAY node)
{
	int ii1;
	int index;

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	RC_TRY( allocate_bc_array(node.size, init_temp) );

	for(ii1 = 0; ii1 < node.size; ii1++){
		if(node.array[ii1].label < 0){
			init_temp->array[ii1].node = -1;
			continue;
		}
		index = node.array[ii1].renum;
		if((index < 0)||(index >= node.size)) return(ARG_ERROR_RC);

		init_temp->array[ii1].node = node.array[ii1].label;
		init_temp->array[ii1].s_type[0] = BC_TEMP;
		init_temp->array[ii1].s[0] = node_temp(bc_temp, def_temp,
		                                       node.array[ii1].label);
	}
	init_temp->size = node.size;

	return(NORMAL_RC);
}


/*Kマトリックスの作成*/
RC
make_K_matrix_transient (SKYLINE_MATRIX g_matrix1, SKYLINE_MATRIX g_matrix2,
                         NODE_ARRAY node, ELEMENT_ARRAY element,
                         ELEMENT_ARRAY surf, ELEMENT_ARRAY transfer,
                         MATERIAL_PROP_ARRAY material, double delta_t,
                         PHYSICAL_PROP_ARRAY physical, double theta)
{
	int ii1;
	ELEM_MATRIX elem_matrix;
	/*熱伝導に関するマトリックスの作成*/
	for(ii1=0;ii1<element.size;ii1++){
		if(element.array[ii1].label < 0) continue;

		RC_TRY( make_elem_matrix_therm(element.array[ii1], node, material,
		                               physical, &elem_matrix) );

		/* elem_matrix -> g_matrix */
		RC_TRY( impose_elem_matrix(g_matrix1, elem_matrix, theta) );
		RC_TRY( impose_elem_matrix(g_matrix2, elem_matrix, -(1.0 - theta) ));
		RC_TRY( free_elem_matrix(&elem_matrix) );
	}

	/*熱伝達に関するマトリックスの作成*/
	for(ii1=0; ii1 < transfer.size; ii1++){
		if(transfer.array[ii1].label < 0) continue;

		RC_TRY( make_elem_matrix_transfer(transfer.array[ii1], node, surf,
		                                  &elem_matrix) );
		RC_TRY( impose_elem_matrix(g_matrix1, elem_matrix, theta) );
		RC_TRY( impose_elem_matrix(g_matrix2, elem_matrix, -(1.0 - theta) ));
		RC_TRY( free_elem_matrix(&elem_matrix) );
	}

	/*熱容量に関するマトリックスの作成*/
	for(ii1=0;ii1<element.size;ii1++){
		if(element.array[ii1].label < 0) continue;

		RC_TRY(make_elem_matrix_heatcapacity(element.array[ii1], node,
		                                     material, physical,&elem_matrix));
		RC_TRY( impose_elem_matrix(g_matrix1, elem_matrix, 1.0/delta_t) );
		RC_TRY( impose_elem_matrix(g_matrix2, elem_matrix, 1.0/delta_t) );
		RC_TRY( free_elem_matrix(&elem_matrix) );
	}
	return(NORMAL_RC);
}


RC
make_elem_matrix_heatcapacity (ELEMENT element, NODE_ARRAY node,
                               MATERIAL_PROP_ARRAY material,
                               PHYSICAL_PROP_ARRAY physical,
                               ELEM_MATRIX *elem_matrix)
{
	int ii1,ii2,ii3,ii4,m_index;
	double tmp;
	WORK_SET ws;
	double *weights = NULL;
	double **gauss_points = NULL;
	double **node_location = NULL;
	int num_points;

	m_index = search_material_prop_element(material, physical, &element);
	if(m_index < 0) return(SEEK_ERROR_RC);
	tmp = material.array[m_index].rho * material.array[m_index].sheat;

	RC_TRY(allocate_work_set(&ws));

	elem_matrix->element = element.label;
	elem_matrix->dim = element_dim(element.type);
	elem_matrix->size = element.node_num;

	/* elem_matrix->array[][] の確保 */
	RC_TRY( allocate2D(elem_matrix->size,
	                   elem_matrix->size, &(elem_matrix->matrix)) );

	/* elem_matirx->index[] の確保と代入 */
	elem_matrix->index = mm_alloc(sizeof(int)*(elem_matrix->size));
	if(elem_matrix->index == NULL) return(ALLOC_ERROR_RC);

	for(ii1=0; ii1<element.node_num; ii1++){
		int index = node_renum_index(node, element.node[ii1]);
		if(index < 0) return(SEEK_ERROR_RC);

		elem_matrix->index[ii1] = index;
	}

	weights = ws.vec_G[0];
	gauss_points = ws.mat_G_4[0];
	node_location = ws.mat_N_3[0];

	RC_TRY( Gauss_const_default(element.type, gauss_points, weights,
	                            &num_points) );

	RC_TRY( node_location_matrix(&(element), node, node_location) );

	for(ii2=0; ii2<num_points; ii2++){
		double *N = ws.vec_N[0];
		double **dN = ws.mat_3_N[0];
		double **jacobi = ws.mat_3_3[0];
		int dim = 0;
		double det_J = 0.0;

		/* 積分の準備 */
		RC_TRY( set_N_vector(element.type, gauss_points[ii2], N) );
		RC_TRY( set_dN_matrix(element.type, gauss_points[ii2], dN) );
		dim = element_dim(element.type);
		mul_matrix(dim, element.node_num, dim, dN, node_location, jacobi);
		det_J = determinant(dim, jacobi);

		for(ii3=0;ii3<element.node_num;ii3++){
			for(ii4=0;ii4<element.node_num;ii4++){
				elem_matrix->matrix[ii4][ii3] += tmp * N[ii3] * N[ii4]
				                               * weights[ii2] * det_J;
			}
		}
	}
	
	RC_TRY(free_work_set(&ws));
	
	return(NORMAL_RC);
}


/* 非定常項のベクトルを計算 */
RC
add_transient_therm_vector (double *s_vect, SKYLINE_MATRIX g_matrix,
                            NODE_ARRAY node, BC_ARRAY temp)
{
	int ii1, dof;
	double *s_vect1,*s_vect2;
	
	RC_TRY( fem_allocate_vector(1, node, &s_vect1) );
	RC_TRY( fem_allocate_vector(1, node, &s_vect2) );

	RC_TRY( temp_vector(s_vect2, temp, node) );

	RC_TRY( s_cholesky_mul_vect(g_matrix.dof, g_matrix.index1,
	                            g_matrix.index2, g_matrix.array, s_vect2,
	                            s_vect1) );

	dof = count_valid_node(node);
	for(ii1=0;ii1<dof;ii1++) s_vect[ii1] += s_vect1[ii1];

	RC_TRY( fem_free_vector(&s_vect1) );
	RC_TRY( fem_free_vector(&s_vect2) );

	return(NORMAL_RC);
}


/* 温度分布をベクトルにする */
static RC
temp_vector (double *f_vect, BC_ARRAY temp, NODE_ARRAY node)
{
	int ii1,ii2;
	int index;
	
	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

	for(ii1 = 0; ii1 < temp.size; ii1++){
		if(temp.array[ii1].node < 0) continue;

		index = node_renum_index(node, temp.array[ii1].node);
		
		for(ii2=0; ii2<BC_SCALAR; ii2++){
			if(temp.array[ii1].s_type[ii2] == BC_TEMP){
				f_vect[index] = temp.array[ii1].s[ii2];
			}
		}
	}

	return(NORMAL_RC);
}
