/*********************************************************************
 * fem_solver.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Takaaki NAGATANI>
 *
 * Thanks to following contributors.
 *  <Tomoyuki TSUBATA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: fem_solver.c,v 1.27 2004/02/27 06:48:34 nagatani Exp $ */

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

#define DEF_ANA_2D    ANA_PLANE_STRAIN


/* 要素上のガウスポイント毎の情報 */
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
static RC fill_row_nodes(FEM_ELEM_LOOKUP elem_lookup,
                         FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                         int row_nodes[], int *row_size);
static RC element_body_force(int gdim, FEM_ELEMENT element, FEM_NODE_ARRAY node,
                             FEM_MATERIAL_PROP_ARRAY material,
                             FEM_PHYSICAL_PROP_ARRAY physical,
                             FEM_DEFAULT_BC_ARRAY def_b_force,
                             int *index_array, int *vect_size,
                             double *body_force);
static RC element_temp_load(int gdim, FEM_ELEMENT element, FEM_NODE_ARRAY node,
                            FEM_MATERIAL_PROP_ARRAY material,
                            FEM_PHYSICAL_PROP_ARRAY physical,
                            FEM_BC_ARRAY temp, FEM_DEFAULT_BC_ARRAY def_temp,
                            int *index_array, int *vect_size,
                            double *temp_load);
static RC element_strain_load(int gdim, FEM_ELEMENT element,
                              FEM_NODE_ARRAY node,
                              FEM_MATERIAL_PROP_ARRAY material,
                              FEM_PHYSICAL_PROP_ARRAY physical,
                              FEM_STRAIN_ARRAY init_strain,
                              int *index_array, int *vect_size,
                              double *strain_load);
static double node_temp(FEM_BC_ARRAY temp, FEM_DEFAULT_BC_ARRAY def_temp,
                                                      int node_label);
static RC fill_index1(int dim, FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                      long *index1);
static RC RI_B_matrix(FEM_ELEMENT elem, G_PROP g_prop);
static RC set_g_prop(FEM_ELEMENT element, double **node_location,
                     G_PROP *g_prop);
static RC set_Gauss_const(FEM_ELEM_TYPE type, G_PROP *g_prop);
static RC set_B_matrix(int dim, int node_num,
                       double **dNxyz, double **B_matrix);
static RC mod_g_matrix(FEM_SKYLINE_MATRIX g_matrix);
static RC Gauss_const4mass_matrix(FEM_ELEM_TYPE type, double **points,
                                  double *weight, int *num_point);


RC fem_react_force(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                   FEM_MATERIAL_PROP_ARRAY material,
                   FEM_PHYSICAL_PROP_ARRAY physical,
                   FEM_DISP_ARRAY disp, FEM_REACT_ARRAY *react)
{
	double *r_vect;
	double *d_vect;
	FEM_ELEM_MATRIX elem_matrix;
	int dim;
	int ii1, ii2, ii3;

	dim = fem_analysis_dim(element);
	RC_TRY( fem_disp2vect(dim, node, disp, &d_vect) );
	RC_TRY( fem_allocate_vector(dim, node, &r_vect) );

	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		RC_TRY( fem_make_elem_matrix(element.array[ii1], node,
		                             material, physical, &elem_matrix) );
		for(ii2=0; ii2<elem_matrix.size; ii2++){
			for(ii3=0; ii3<elem_matrix.size; ii3++){
				r_vect[elem_matrix.index[ii2]]
				      += elem_matrix.matrix[ii2][ii3]
				        *d_vect[elem_matrix.index[ii3]];
			}
		}
		RC_TRY( free_fem_elem_matrix(&elem_matrix) );
	}

	RC_TRY( fem_vect2react(dim, node, r_vect, react) );
	RC_TRY( fem_free_vector(&r_vect) );
	RC_TRY( fem_free_vector(&d_vect) );

	return(NORMAL_RC);
}

RC fem_solve_iccg3(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                   FEM_MATERIAL_PROP_ARRAY material,
                   FEM_PHYSICAL_PROP_ARRAY physical,
                   FEM_BC_ARRAY rest, FEM_BC_ARRAY force,
                   FEM_DEFAULT_BC_ARRAY b_force,
                   FEM_BC_ARRAY temp, FEM_DEFAULT_BC_ARRAY def_temp,
                   FEM_DISP_ARRAY *disp, FEM_STRESS_ARRAY *stress,
                   FEM_STRAIN_ARRAY *strain, int fem_sol_ss_flag)
{
	int ii1;
	NONZERO_MATRIX3 matrix;
	double *f_vect;
	double *d_vect;
	FEM_ELEM_MATRIX elem_matrix;
	int dim;

	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

	/* 全体剛性マトリックスの作成 */
	RC_TRY( fem_allocate_nonzero3(node, element, &matrix) );
	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		RC_TRY( fem_make_elem_matrix(element.array[ii1], node,
		                             material, physical, &elem_matrix) );
		RC_TRY( fem_impose_nonzero3(matrix, elem_matrix, 1.0) );
		RC_TRY( free_fem_elem_matrix(&elem_matrix) );
	}

	/* 荷重ベクトルの作成 */
	dim = fem_analysis_dim(element);
	if(dim != 3) return(ARG_ERROR_RC);
	RC_TRY( fem_allocate_vector(dim, node, &f_vect) );
	RC_TRY( fem_add_nodal_force(dim, node, force, f_vect) );
	RC_TRY( fem_add_body_force(dim, node, element, material, physical,
	                           b_force, f_vect) );
	RC_TRY( fem_add_temp_force(dim, node, element, material, physical,
	                           temp, def_temp, f_vect) );

	/* 変位の計算 */
	RC_TRY( fem_modify_nonzero3(node, rest, matrix, f_vect) );
	/*RC_TRY( fem_iccg_nonzero3(matrix, f_vect, &d_vect, NULL) );*/
	RC_TRY( fem_scale_iccg_nonzero3(matrix, f_vect, &d_vect, NULL) );
	RC_TRY( free_nonzero_matrix3(&matrix) );
	RC_TRY( fem_free_vector(&f_vect) );
	RC_TRY( fem_vect2disp(3, node, d_vect, disp) );
	RC_TRY( fem_free_vector(&d_vect) );

	/* ひずみ、応力の計算 */
	RC_TRY( fem_stress_strain(fem_sol_ss_flag, node, element, *disp,
	                          material, physical, stress, strain) );

	return(NORMAL_RC);
}


RC fem_solve(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
             FEM_MATERIAL_PROP_ARRAY material,
             FEM_PHYSICAL_PROP_ARRAY physical,
             FEM_BC_ARRAY rest, FEM_BC_ARRAY force,
             FEM_BC_ARRAY temp, FEM_DEFAULT_BC_ARRAY def_temp,
             FEM_DISP_ARRAY *disp, FEM_STRESS_ARRAY *stress,
             FEM_STRAIN_ARRAY *strain, int fem_sol_ss_flag)
{
	FEM_SKYLINE_MATRIX K_matrix;
	double *f_vect;

	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );
	RC_TRY( fem_make_g_matrix(node, element, material, physical, &K_matrix) );
	RC_TRY( fem_allocate_vector(K_matrix.dim, node, &f_vect) );
	RC_TRY( fem_add_nodal_force(K_matrix.dim, node, force, f_vect) );
	RC_TRY( fem_add_temp_force(K_matrix.dim, node, element, material, physical,
	                           temp, def_temp, f_vect) );
	RC_TRY( fem_modify_g_matrix(node, rest, K_matrix, f_vect) );
	RC_TRY( fem_decomp_g_matrix(K_matrix) );
	RC_TRY( fem_subst_g_matrix(node, K_matrix, f_vect, disp) );
	RC_TRY( fem_free_g_matrix(&K_matrix) );
	RC_TRY( fem_stress_strain(fem_sol_ss_flag, node, element, *disp,
	                          material, physical, stress, strain) );
	free(f_vect);

	return(NORMAL_RC);
}


RC fem_solve_therm(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                   FEM_MATERIAL_PROP_ARRAY material,
                   FEM_PHYSICAL_PROP_ARRAY physical,
                   FEM_BC_ARRAY bc_temp, FEM_BC_ARRAY *result_temp)
{
	FEM_SKYLINE_MATRIX g_matrix;
	double *s_vect;
	FEM_ELEM_MATRIX elem_matrix;
	int ii1;

	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

	/* g_matrix の確保と初期化 */
	RC_TRY( fem_allocate_g_matrix(1, node, element, &g_matrix) );

	for(ii1=0;ii1<element.size;ii1++){
		if(element.array[ii1].label < 0) continue;

		RC_TRY( fem_make_elem_matrix_therm(element.array[ii1], node, material,
		                                   physical, &elem_matrix) );

		/* elem_matrix -> g_matrix */
		RC_TRY( fem_impose_elem_matrix(g_matrix, elem_matrix, 1.0) );
		RC_TRY( free_fem_elem_matrix(&elem_matrix) );
	}

	RC_TRY( fem_allocate_vector(1, node, &s_vect) );
	RC_TRY( fem_modify_g_matrix_therm(node, bc_temp, g_matrix, s_vect) );

	RC_TRY( fem_decomp_g_matrix(g_matrix) );
	RC_TRY( fem_subst_g_matrix_therm(node, g_matrix, s_vect, result_temp) );
	RC_TRY( fem_free_g_matrix(&g_matrix) );
	RC_TRY( fem_free_vector(&s_vect) );

	return(NORMAL_RC);
}


/* 複数荷重条件の解析 */
/* forces, disps, stresses, strains は、num_force 個以上の配列 */
RC fem_solve_multi_force(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                         FEM_MATERIAL_PROP_ARRAY material,
                         FEM_PHYSICAL_PROP_ARRAY physical,
                         FEM_BC_ARRAY rest,
                         int num_force, const FEM_BC_ARRAY *forces,
                         FEM_BC_ARRAY temp, FEM_DEFAULT_BC_ARRAY def_temp,
                         FEM_DISP_ARRAY *disps, FEM_STRESS_ARRAY *stresses,
                         FEM_STRAIN_ARRAY *strains, int fem_sol_ss_flag)
{
	FEM_SKYLINE_MATRIX K_matrix;
	double **f_vect;
	int ii1;

	if(num_force <= 0) return(ARG_ERROR_RC);
	f_vect = (double **)malloc(num_force*sizeof(double *));
	if(f_vect == NULL){
		fprintf(stderr, "malloc() < \n");
		return(ALLOC_ERROR_RC);
	}

	if(node.renum_flag != RENUMBERED){
		RC_TRY( dummy_renumber(&node) );
	}
	RC_TRY( fem_make_g_matrix(node, element, material, physical, &K_matrix) );
	for(ii1=0; ii1<num_force; ii1++){
		RC_TRY( fem_allocate_vector(K_matrix.dim, node, &(f_vect[ii1])) );
		RC_TRY( fem_add_nodal_force(K_matrix.dim, node,
		                            forces[ii1], f_vect[ii1]) );
		RC_TRY( fem_add_temp_force(K_matrix.dim, node, element,
		                           material, physical,
		                           temp, def_temp, f_vect[ii1]) );
	}
	RC_TRY( fem_modify_g_matrix_n(node, rest, K_matrix, f_vect, num_force) );
	RC_TRY( fem_decomp_g_matrix(K_matrix) );

	for(ii1=0; ii1<num_force; ii1++){
		RC_TRY( fem_subst_g_matrix(node, K_matrix,
		                           f_vect[ii1], &(disps[ii1])) );
		free(f_vect[ii1]);
	}
	
	RC_TRY( fem_free_g_matrix(&K_matrix) );

	for(ii1=0; ii1<num_force; ii1++){
		RC_TRY( fem_stress_strain(fem_sol_ss_flag, node, element, disps[ii1],
		                          material, physical,
		                          &(stresses[ii1]), &(strains[ii1])) );
	}
	free(f_vect);

	return(NORMAL_RC);
}


RC fem_free_g_matrix(FEM_SKYLINE_MATRIX *g_matrix)
{
	RC_TRY( s_cholesky_free(g_matrix->index2, g_matrix->array) );
	free(g_matrix->index1);

	g_matrix->dim = 0;
	g_matrix->dof = 0;
	g_matrix->index1 = NULL;
	g_matrix->index2 = NULL;
	g_matrix->array = NULL;

	return(NORMAL_RC);
}


/* 応力、ひずみの計算 */
/* stress, strain が NULL なら代入しない */
RC fem_stress_strain(int fem_sol_ss_flag, FEM_NODE_ARRAY node,
                     FEM_ELEMENT_ARRAY element, FEM_DISP_ARRAY disp,
                     FEM_MATERIAL_PROP_ARRAY material,
                     FEM_PHYSICAL_PROP_ARRAY physical,
                     FEM_STRESS_ARRAY *stress, FEM_STRAIN_ARRAY *strain)
{
	int ii1, ii2;
	static int init_flag = 0;
	static double **local_points;
	FEM_STRESS_STRAIN local_stress, local_strain;
	TRANS_ROTATION3D *stress_sum = NULL;
	TRANS_ROTATION3D *strain_sum = NULL;
	int *count = NULL;
	TRANS_ROTATION3D init_v = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	int node_index;

	if( ((stress == NULL)&&(strain == NULL)) || (fem_sol_ss_flag == 0) ){
		return(NORMAL_RC);
	}
	if(init_flag == 0){
		RC_TRY( allocate2D(FEM_MAX_NODE, 4, &local_points) );
		init_flag = 1;
	}

	if(fem_sol_ss_flag & FEM_SOL_SS_NODE_AV){
		RC_NULL_CHK( stress_sum = (TRANS_ROTATION3D *)
				              malloc((node.size)*sizeof(TRANS_ROTATION3D)) );
		RC_NULL_CHK( strain_sum = (TRANS_ROTATION3D *)
		                      malloc((node.size)*sizeof(TRANS_ROTATION3D)) );
		RC_NULL_CHK( count = (int *) malloc((node.size)*sizeof(int)) );
		for(ii1=0; ii1<node.size; ii1++){
			count[ii1] = 0;
			stress_sum[ii1] = strain_sum[ii1] = init_v;
		}
	}

	if(stress != NULL){
		RC_TRY( allocate_fem_stress_array(0, stress) );
		stress->source = FEM_K_SOL;
	}
	if(strain != NULL){
		RC_TRY( allocate_fem_strain_array(0, strain) );
		strain->source = FEM_K_SOL;
	}

	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		/* 要素中心 */
		if(fem_sol_ss_flag & FEM_SOL_SS_CENTER){
			RC_TRY( set_center_point(element.array[ii1].type,
			                         local_points[0]) );
			RC_TRY( fem_local_stress_strain(element.array[ii1], local_points[0],
			                                node, disp, material, physical,
			                                &local_stress, &local_strain) );
			if(stress != NULL){
				RC_TRY( realloc_fem_stress_array(stress) );
				stress->array[stress->size] = local_stress;
				stress->array[stress->size].node = -1;  /* center */
				(stress->size)++;
			}
			if(strain != NULL){
				RC_TRY( realloc_fem_strain_array(strain) );
				strain->array[strain->size] = local_strain;
				strain->array[strain->size].node = -1;  /* center */
				(strain->size)++;
			}
		}

		/* 各節点 */
		if( (fem_sol_ss_flag & FEM_SOL_SS_NODE)
		  ||(fem_sol_ss_flag & FEM_SOL_SS_NODE_AV) ){
			RC_TRY( set_local_node_points(element.array[ii1].type,
			                              local_points) );

			for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
				RC_TRY( fem_local_stress_strain(element.array[ii1],
				            local_points[ii2], node, disp, material, physical,
				            &local_stress, &local_strain) );

				if(fem_sol_ss_flag & FEM_SOL_SS_NODE_AV){
					node_index = search_fem_node_label(node,
					                   element.array[ii1].node[ii2]);
					RC_NEG_CHK( node_index );
					count[node_index]++;
					stress_sum[node_index]
					    = add_trans_rotation3d(stress_sum[node_index],
					                           local_stress.v);
					strain_sum[node_index]
					    = add_trans_rotation3d(strain_sum[node_index],
					                           local_strain.v);
				}

				if(stress != NULL){
					RC_TRY( realloc_fem_stress_array(stress) );
					stress->array[stress->size] = local_stress;
					stress->array[stress->size].node
					      = element.array[ii1].node[ii2];
					(stress->size)++;
				}

				if(strain != NULL){
					RC_TRY( realloc_fem_strain_array(strain) );
					strain->array[strain->size] = local_strain;
					strain->array[strain->size].node
					      = element.array[ii1].node[ii2];
					(strain->size)++;
				}
			}
		}
	}

	if(stress != NULL) RC_TRY( clean_fem_stress_array(stress) );
	if(strain != NULL) RC_TRY( clean_fem_strain_array(strain) );

	if(fem_sol_ss_flag & FEM_SOL_SS_NODE_AV){
		if(stress != NULL){
			for(ii1=0; ii1<stress->size; ii1++){
				if(stress->array[ii1].element < 0) continue;
				node_index = search_fem_node_label(node,
				                                   stress->array[ii1].node);
				if(count[node_index] <= 0) return(UNKNOWN_ERROR_RC);
				stress->array[ii1].v
				 = mul_scalar_trans_rotation3d(1.0/(double)(count[node_index]),
				                               stress_sum[node_index]);
			}
		}
		if(strain != NULL){
			for(ii1=0; ii1<strain->size; ii1++){
				if(strain->array[ii1].element < 0) continue;
				node_index = search_fem_node_label(node,
				                                   strain->array[ii1].node);
				if(count[node_index] <= 0) return(UNKNOWN_ERROR_RC);
				strain->array[ii1].v
				 = mul_scalar_trans_rotation3d(1.0/(double)(count[node_index]),
				                               strain_sum[node_index]);
			}
		}
		free(stress_sum);
		free(strain_sum);
		free(count);
	}

	return(NORMAL_RC);
}


/* 要素の局所ひずみ */
RC fem_local_strain(FEM_ELEMENT element, const double *local_point,
                    FEM_NODE_ARRAY node, FEM_DISP_ARRAY disp,
                    FEM_STRESS_STRAIN *strain)
{
	double disp_v[3*FEM_MAX_NODE];
	static double **node_location = NULL;
	static double **B_matrix = NULL;
	static int init_flag = 0;

	if(init_flag == 0){
		RC_TRY( allocate2D(FEM_MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(6, 3*FEM_MAX_NODE, &B_matrix) );

		init_flag = 1;
	}

	if(element.label < 0) return(ARG_ERROR_RC);
	if(element_dim(element.type) <= 1) return(ARG_ERROR_RC);
	if(element_ssnum(element.type) <= 0) return(ARG_ERROR_RC);

	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( fill_disp_vector(element, disp, disp_v) );

	RC_TRY( fem_local_strain_sub(element, local_point, node_location,
	                             disp_v, B_matrix, strain) );

	return(NORMAL_RC);
}


RC fem_local_strain_sub(FEM_ELEMENT element, const double local_point[],
                        double **node_location, const double disp_v[],
                        double **ws_mat_6_3N, FEM_STRESS_STRAIN *strain)
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
	init_fem_stress_strain(strain);
	strain->element = element.label;
	RC_TRY( vect2stress_strain(ssnum, strain_v, strain) );

	return(NORMAL_RC);
}

/* 要素の局所ひずみ -> 局所応力を計算 */
RC fem_local_strain2local_stress(FEM_ELEMENT element,
                                 FEM_MATERIAL_PROP_ARRAY material,
                                 FEM_PHYSICAL_PROP_ARRAY physical,
                                 FEM_STRESS_STRAIN strain,
                                 FEM_STRESS_STRAIN *stress)
{
	int ssnum;
	double strain_v[6];
	double stress_v[6];
	double **D_matrix = NULL;

	if(element.label < 0) return(ARG_ERROR_RC);

	RC_TRY( allocate2D(6, 6, &D_matrix) );

	if( (ssnum = element_ssnum(element.type)) <= 0){
		fprintf(stderr, "element_ssnum < ");
		return(ARG_ERROR_RC);
	}

	init_fem_stress_strain(stress);
	stress->element = element.label;
	RC_TRY( stress_strain2vect(ssnum, strain, strain_v) );

	/* D_matrix */
	RC_TRY( get_D_matrix(element, material, physical, D_matrix) );

	/* D_matrix, strain_v -> stress_v */
	mul_matrix_vect(ssnum, ssnum, D_matrix, strain_v, stress_v);

	/* stress_v -> stress */
	RC_TRY( vect2stress_strain(ssnum, stress_v, stress) );

	RC_TRY( free2D(6, 6, D_matrix) );

	return(NORMAL_RC);
}


/* 要素の局所応力、ひずみ（熱歪のfact倍を差し引いた値） */
RC fem_local_stress_strain_therm(FEM_ELEMENT element,
                             const double *local_point,
                             FEM_NODE_ARRAY node, FEM_DISP_ARRAY disp,
                             FEM_MATERIAL_PROP_ARRAY material,
                             FEM_PHYSICAL_PROP_ARRAY physical,
                             FEM_BC_ARRAY temp, FEM_DEFAULT_BC_ARRAY def_temp,
                             double fact,
                             FEM_STRESS_STRAIN *stress,
                             FEM_STRESS_STRAIN *strain)
{
	int ii1;
	int dim, ssnum;
	double disp_v[3*FEM_MAX_NODE];
	double strain_v[6];
	double stress_v[6];
	static double **node_location = NULL;
	static double **B_matrix = NULL;
	static double **D_matrix = NULL;
	static int init_flag = 0;
	double node_temp_array[FEM_MAX_NODE];
	double N[FEM_MAX_NODE];
	double local_temp;
	int m_index;
	FEM_MATERIAL_PROP mat;

	if(element.label < 0) return(ARG_ERROR_RC);

	if(init_flag == 0){
		RC_TRY( allocate2D(FEM_MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(6, 3*FEM_MAX_NODE, &B_matrix) );
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
	RC_NEG_CHK( m_index = search_fem_material_prop_element(material,
	                                            physical, &element) );
	mat = material.array[m_index];
	RC_TRY( fem_fill_material(&mat) );

	if((dim = element_dim(element.type)) <= 1){
		fprintf(stderr, "element_dim < ");
		return(ARG_ERROR_RC);
	}
	if( (ssnum = element_ssnum(element.type)) <= 0){
		fprintf(stderr, "element_ssnum < ");
		return(ARG_ERROR_RC);
	}

	init_fem_stress_strain(stress);
	stress->element = element.label;

	init_fem_stress_strain(strain);
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
RC fem_local_stress_strain(FEM_ELEMENT element, const double *local_point,
                           FEM_NODE_ARRAY node, FEM_DISP_ARRAY disp,
                           FEM_MATERIAL_PROP_ARRAY material,
                           FEM_PHYSICAL_PROP_ARRAY physical,
                           FEM_STRESS_STRAIN *stress,
                           FEM_STRESS_STRAIN *strain)
{
	int dim, ssnum;
	double disp_v[3*FEM_MAX_NODE];
	double strain_v[6];
	double stress_v[6];
	static double **node_location = NULL;
	static double **B_matrix = NULL;
	static double **D_matrix = NULL;
	static int init_flag = 0;

	if(element.label < 0) return(ARG_ERROR_RC);

	if(init_flag == 0){
		RC_TRY( allocate2D(FEM_MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(6, 3*FEM_MAX_NODE, &B_matrix) );
		RC_TRY( allocate2D(6, 6, &D_matrix) );

		init_flag = 1;
	}

	if((dim = element_dim(element.type)) <= 1){
		fprintf(stderr, "element_dim < ");
		return(ARG_ERROR_RC);
	}
	if( (ssnum = element_ssnum(element.type)) <= 0){
		fprintf(stderr, "element_ssnum < ");
		return(ARG_ERROR_RC);
	}

	init_fem_stress_strain(stress);
	stress->element = element.label;

	init_fem_stress_strain(strain);
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


TRANS_ROTATION3D d_mises_stress(FEM_STRESS_STRAIN stress)
{
	double fact;
	TRANS_ROTATION3D ret;

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

double mises_stress(FEM_STRESS_STRAIN stress)
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


RC vect2stress_strain(int vect_size, const double *vect,
                      FEM_STRESS_STRAIN *ss)
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


RC stress_strain2vect(int vect_size, FEM_STRESS_STRAIN ss, double *vect)
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


RC fill_disp_vector(FEM_ELEMENT element,
                    FEM_DISP_ARRAY disp, double *disp_v)
{
	int dim;
	int ii1;

	if((dim = element_dim(element.type)) <= 1){
		fprintf(stderr, "element_dim < ");
		return(ARG_ERROR_RC);
	}

	for(ii1=0; ii1<element.node_num; ii1++){
		int index;

		RC_NEG_CHK( index = search_fem_disp_node(disp, element.node[ii1]) );

		disp_v[dim*ii1] = disp.array[index].v.x;
		disp_v[dim*ii1 + 1] = disp.array[index].v.y;
		if(dim == 3){
			disp_v[dim*ii1 + 2] = disp.array[index].v.z;
		}
	}

	return(NORMAL_RC);
}


/* 剛性マトリックスをLU分解 */
RC fem_decomp_g_matrix(FEM_SKYLINE_MATRIX g_matrix)
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
RC fem_subst_g_matrix(FEM_NODE_ARRAY node, FEM_SKYLINE_MATRIX g_matrix,
                      double *f_vect, FEM_DISP_ARRAY *disp)
{
	int ii1;
	int index;

	RC_TRY( s_cholesky_subst(g_matrix.dof, g_matrix.index1, g_matrix.index2,
	                         g_matrix.array, f_vect, 1) );

	/* f_vect -> disp */
	RC_TRY( allocate_fem_disp_array(node.size, disp) );
	disp->source = FEM_K_SOL;
	for(ii1=0; ii1<node.size; ii1++){
		if(node.array[ii1].label < 0){
			disp->array[ii1].node = -1;
			continue;
		}
		index = node.array[ii1].renum;
		if((index < 0)||(index >= node.size)){
			fprintf(stderr, "node.array[index].renum < ");
			return(ARG_ERROR_RC);
		}

		disp->array[ii1].node = node.array[ii1].label;
		disp->array[ii1].v.x = f_vect[g_matrix.dim*index];
		disp->array[ii1].v.y = f_vect[g_matrix.dim*index + 1];
		if(g_matrix.dim == 3){
			disp->array[ii1].v.z = f_vect[g_matrix.dim*index + 2];
		}
	}
	disp->size = node.size;
	RC_TRY( clean_fem_disp_array(disp) );

	return(NORMAL_RC);
}


/* g_matrix, s_vect -> s_vect, temp */
RC fem_subst_g_matrix_therm(FEM_NODE_ARRAY node, FEM_SKYLINE_MATRIX g_matrix,
                            double *s_vect, FEM_BC_ARRAY *temp)
{
	int ii1;
	int index;

	RC_TRY( s_cholesky_subst(g_matrix.dof, g_matrix.index1, g_matrix.index2,
	                         g_matrix.array, s_vect, 1) );

	/* s_vect -> temp */
	RC_TRY( allocate_fem_bc_array(node.size, temp) );
	temp->source = FEM_K_SOL;
	for(ii1=0; ii1<node.size; ii1++){
		if(node.array[ii1].label < 0){
			temp->array[ii1].node = -1;
			continue;
		}
		index = node.array[ii1].renum;
		if((index < 0)||(index >= node.size)){
			fprintf(stderr, "node.array[ii1].renum < ");
			return(ARG_ERROR_RC);
		}

		temp->array[ii1].node = node.array[ii1].label;
		temp->array[ii1].s_type[0] = BC_TEMP;
		temp->array[ii1].s[0] = s_vect[index];
	}
	temp->size = node.size;
	RC_TRY( clean_fem_bc_array(temp) );

	return(NORMAL_RC);
}


/* 拘束条件を考慮して全体剛性マトリックス、荷重ベクトルを修正 */
RC fem_modify_g_matrix(FEM_NODE_ARRAY node, FEM_BC_ARRAY rest,
                       FEM_SKYLINE_MATRIX g_matrix, double *f_vect)
{
	int ii1;
	int index;

	if(node.renum_flag != RENUMBERED){
		return(ARG_ERROR_RC);
	}

	for(ii1=0; ii1<rest.size; ii1++){
		if(rest.array[ii1].node < 0) continue;

		index = node_renum_index(node, rest.array[ii1].node);
		if(index < 0){
			fprintf(stderr, "node_renum_index() < ");
			return(SEEK_ERROR_RC);
		}

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
RC fem_modify_g_matrix_therm(FEM_NODE_ARRAY node, FEM_BC_ARRAY temp,
                             FEM_SKYLINE_MATRIX g_matrix, double *s_vect)
{
	int ii1, ii2;
	int index;

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	for(ii1=0; ii1<temp.size; ii1++){
		if(temp.array[ii1].node < 0) continue;

		index = node_renum_index(node, temp.array[ii1].node);
		if(index < 0){
			fprintf(stderr, "node_renum_index() < ");
			return(SEEK_ERROR_RC);
		}
		
		for(ii2=0; ii2<FEM_BC_SCALAR; ii2++){
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
/* 複数荷重条件版(fem_modify_g_matrix() との共通化の余地あり) */
RC fem_modify_g_matrix_n(FEM_NODE_ARRAY node, FEM_BC_ARRAY rest,
                       FEM_SKYLINE_MATRIX g_matrix, double **f_vect,
                       int f_vec_size)
{
	int ii1;
	int index;

	if(node.renum_flag != RENUMBERED){
		return(ARG_ERROR_RC);
	}

	for(ii1=0; ii1<rest.size; ii1++){
		if(rest.array[ii1].node < 0) continue;

		index = node_renum_index(node, rest.array[ii1].node);
		if(index < 0){
			fprintf(stderr, "node_renum_index() < ");
			return(SEEK_ERROR_RC);
		}

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
static RC mod_g_matrix(FEM_SKYLINE_MATRIX g_matrix)
{
	RC rc;
	double *ptr;
	int ii1;

	for(ii1=0; ii1<g_matrix.dof; ii1++){
		if(g_matrix.index1[ii1] == ii1){ /* 対角要素のみ */
			ptr = s_cholesky_ptr(g_matrix.dof, g_matrix.index1,
			                     g_matrix.index2, g_matrix.array,
			                     ii1, ii1, &rc);
			if(rc != NORMAL_RC){
				fprintf(stderr, "s_cholesky_ptr < ");
				return(rc);
			}

			if( nearly_eq(*ptr, 0.0) ){
				*ptr = 1.0;
			}
		}
	}

	return(NORMAL_RC);
}


/* G_PROP の確保 */
static RC allocate_g_prop(G_PROP *g_prop)
{
	int ii1;

	RC_TRY( allocate2D(MAX_GAUSS_POINT+1, 4, &(g_prop->g_points)) );

	for(ii1=0; ii1<MAX_GAUSS_POINT+1; ii1++){
		RC_TRY( allocate2D(6, 3*FEM_MAX_NODE, &(g_prop->B_matrix[ii1])) );
	}

	return(NORMAL_RC);
}


RC fem_total_dof(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                 int *total_dof)
{
	int dim;

	/* 問題の次元をチェック */
	dim = fem_analysis_dim(element);
	if(dim <= 1){
		fprintf(stderr, "fem_analysis_dim() < ");
		return(ARG_ERROR_RC);
	}

	*total_dof = dim * node.size;

	return(NORMAL_RC);
}


/* ベクトルの確保と初期化 */
RC fem_allocate_vector(int dim, FEM_NODE_ARRAY node, double **f_vect)
{
	int dof;
	int ii1;

	dof = dim * count_valid_node(node);
	if(dof < 0) return(ARG_ERROR_RC);

	(*f_vect) = (double *)malloc(dof * sizeof(double));
	if((*f_vect) == NULL){
		fprintf(stderr, "malloc() < ");
		return(ALLOC_ERROR_RC);
	}
	for(ii1=0; ii1<dof; ii1++){
		(*f_vect)[ii1] = 0.0;
	}

	return(NORMAL_RC);
}


/* 荷重ベクトルの初期化 */
RC fem_zero_vector(int dim, FEM_NODE_ARRAY node, double f_vect[])
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


RC fem_free_vector(double **vect)
{
	if((*vect) == NULL) return(ARG_ERROR_RC);

	free((void *)(*vect));
	(*vect) = NULL;

	return(NORMAL_RC);
}


RC fem_force_vector2bc(const double *f_vect, FEM_ELEMENT_ARRAY element,
                       FEM_NODE_ARRAY node, FEM_BC_ARRAY *bc)
{
	int ii1, ii2;
	int dim;
	int index;
	int nonzero_flag;

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	/* 問題の次元をチェック */
	dim = fem_analysis_dim(element);
	if(dim <= 1){
		fprintf(stderr, "fem_analysis_dim() < ");
		return(ARG_ERROR_RC);
	}

	RC_TRY( allocate_fem_bc_array(0, bc) );

	for(ii1=0; ii1<node.size; ii1++){
		if(node.array[ii1].label < 0) continue;

		RC_NEG_CHK( index = node.array[ii1].renum );
		nonzero_flag = 0;
		for(ii2=0; ii2<dim; ii2++){
			if( !nearly_eq(f_vect[index*dim+ii2], 0.0) ) nonzero_flag = 1;
		}
		if(nonzero_flag){
			RC_TRY( realloc_fem_bc_array(bc) );
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
	RC_TRY( clean_fem_bc_array(bc) );

	return(NORMAL_RC);
}


/* 荷重ベクトルに節点荷重を加算 */
RC fem_add_nodal_force(int dim, FEM_NODE_ARRAY node,
                       FEM_BC_ARRAY force, double *f_vect)
{
	int index;
	int ii1;

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	for(ii1=0; ii1<force.size; ii1++){
		if(force.array[ii1].node < 0) continue;

		index = node_renum_index(node, force.array[ii1].node);
		if(index < 0){
			fprintf(stderr, "node_renum_index() < ");
			return(SEEK_ERROR_RC);
		}
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


RC fem_add_init_strain_force(int dim, FEM_NODE_ARRAY node,
                             FEM_ELEMENT_ARRAY element,
                             FEM_MATERIAL_PROP_ARRAY material,
                             FEM_PHYSICAL_PROP_ARRAY physical,
                             FEM_STRAIN_ARRAY init_strain, double *f_vect)
{
	double strain_load[3*FEM_MAX_NODE];
	int index_array[3*FEM_MAX_NODE];
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
static RC element_strain_load(int gdim, FEM_ELEMENT element,
                              FEM_NODE_ARRAY node,
                              FEM_MATERIAL_PROP_ARRAY material,
                              FEM_PHYSICAL_PROP_ARRAY physical,
                              FEM_STRAIN_ARRAY init_strain,
                              int *index_array, int *vect_size,
                              double *strain_load)
{
	int dim;
	int ii1, ii2;
	double thickness;
	int m_index;
	TRANS_ROTATION3D strain_v[FEM_MAX_NODE];
	const TRANS_ROTATION3D zero_v = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	static int init_flag = 0;
	static double **D_matrix = NULL;
	static double **node_location = NULL;
	static double **trans_B = NULL;
	static double **BtD = NULL;
	static G_PROP g_prop;

	if(init_flag == 0){
		RC_TRY( allocate2D(6, 6, &D_matrix) );
		RC_TRY( allocate2D(FEM_MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(3*FEM_MAX_NODE, 6, &trans_B) );
		RC_TRY( allocate2D(3*FEM_MAX_NODE, 6, &BtD) );
		RC_TRY( allocate_g_prop(&g_prop) );
		init_flag = 1;
	}

	if((dim = element_dim(element.type)) <= 1){
		fprintf(stderr, "element_dim() < ");
		return(ARG_ERROR_RC);
	}
	if(gdim < dim) return(ARG_ERROR_RC);
	*vect_size = dim*element.node_num;

	RC_NEG_CHK( m_index = search_fem_material_prop_element(material,
	                                            physical, &element) );

	/* strain_load を初期化 */
	/* index_array を初期化 */
	/* strain_v を初期化 */
	for(ii1=0; ii1<element.node_num; ii1++){
		int index;
		int str_index;

		RC_NEG_CHK( index = node_renum_index(node, element.node[ii1]) );

		for(ii2=0; ii2<dim; ii2++){
			strain_load[dim*ii1 + ii2] = 0.0;
			index_array[dim*ii1 + ii2] = gdim*index + ii2;
		}

		str_index = search_fem_strain_element_node(init_strain,
		                                           element.label,
		                                           element.node[ii1]);
		if(str_index < 0){
			str_index = search_fem_strain_element_node(init_strain,
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
		TRANS_ROTATION3D local_strain;
		double local_force[3*FEM_MAX_NODE];
		double N[FEM_MAX_NODE];
		double init_strain[6];

		/* ガウスポイントの初期ひずみを内装して求める */
		RC_TRY( set_N_vector(element.type, g_prop.g_points[ii1], N) );
		local_strain = zero_v;
		for(ii2=0; ii2<element.node_num; ii2++){
			TRANS_ROTATION3D tmp_strain;

			tmp_strain = mul_scalar_trans_rotation3d(N[ii2], strain_v[ii1]);
			local_strain = add_trans_rotation3d(local_strain, tmp_strain);
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
			strain_load[ii2] += local_force[ii2] * g_prop.det_J[ii1]
			                  * g_prop.weight[ii1] * thickness;
		}
	}

	return(NORMAL_RC);
}


/* 熱ひずみによる荷重を加算 */
RC fem_add_temp_force(int dim, FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                      FEM_MATERIAL_PROP_ARRAY material,
                      FEM_PHYSICAL_PROP_ARRAY physical, FEM_BC_ARRAY temp,
                      FEM_DEFAULT_BC_ARRAY def_temp, double *f_vect)
{
	double temp_load[3*FEM_MAX_NODE];
	int index_array[3*FEM_MAX_NODE];
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
RC fem_add_body_force(int dim, FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                      FEM_MATERIAL_PROP_ARRAY material,
                      FEM_PHYSICAL_PROP_ARRAY physical,
                      FEM_DEFAULT_BC_ARRAY def_b_force, double *f_vect)
{
	double body_force[3*FEM_MAX_NODE];
	int index_array[3*FEM_MAX_NODE];
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


static RC element_body_force(int gdim, FEM_ELEMENT element, FEM_NODE_ARRAY node,
                             FEM_MATERIAL_PROP_ARRAY material,
                             FEM_PHYSICAL_PROP_ARRAY physical,
                             FEM_DEFAULT_BC_ARRAY def_b_force,
                             int *index_array, int *vect_size,
                             double *body_force)
{
	int dim;
	int ii1, ii2;
	double thickness;
	int m_index;
	static int allocate_flag = 0;
	static double **node_location = NULL;
	static G_PROP g_prop;
	TRANSLATION3D b_force_vect = {0.0, 0.0, 0.0};

	if(allocate_flag == 0){
		RC_TRY( allocate2D(FEM_MAX_NODE, 3, &node_location) );
		RC_TRY( allocate_g_prop(&g_prop) );

		allocate_flag = 1;
	}

	/* vect_size, body_force[], index_array[] をセット */
	if((dim = element_dim(element.type)) <= 1){
		fprintf(stderr, "element_dim() < ");
		return(ARG_ERROR_RC);
	}
	if(gdim < dim) return(ARG_ERROR_RC);
	*vect_size = dim*element.node_num;
	for(ii1=0; ii1<element.node_num; ii1++){
		int index;

		index = node_renum_index(node, element.node[ii1]);
		if(index < 0){
			fprintf(stderr, "node_renum_index() < ");
			return(SEEK_ERROR_RC);
		}

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
	RC_NEG_CHK( m_index = search_fem_material_prop_element(material,
	                                            physical, &element) );
	b_force_vect.x *= material.array[m_index].rho;
	b_force_vect.y *= material.array[m_index].rho;
	b_force_vect.z *= material.array[m_index].rho;

	RC_TRY( get_thickness(element, physical, &thickness) );
	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( set_Gauss_const(element.type, &g_prop) );
	RC_TRY( set_g_prop(element, node_location, &g_prop) );

	for(ii1=0; ii1<g_prop.num_g_point; ii1++){
		double N[FEM_MAX_NODE];

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
static RC element_temp_load(int gdim, FEM_ELEMENT element, FEM_NODE_ARRAY node,
                            FEM_MATERIAL_PROP_ARRAY material,
                            FEM_PHYSICAL_PROP_ARRAY physical,
                            FEM_BC_ARRAY temp, FEM_DEFAULT_BC_ARRAY def_temp,
                            int *index_array, int *vect_size,
                            double *temp_load)
{
	int dim;
	int ii1, ii2;
	double thickness;
	int m_index;
	double node_temp_array[FEM_MAX_NODE];
	static int init_flag = 0;
	static double **D_matrix = NULL;
	static double **node_location = NULL;
	static double **trans_B = NULL;
	static double **BtD = NULL;
	static G_PROP g_prop;
	FEM_MATERIAL_PROP mat;

	if(init_flag == 0){
		RC_TRY( allocate2D(6, 6, &D_matrix) );
		RC_TRY( allocate2D(FEM_MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(3*FEM_MAX_NODE, 6, &trans_B) );
		RC_TRY( allocate2D(3*FEM_MAX_NODE, 6, &BtD) );
		RC_TRY( allocate_g_prop(&g_prop) );
		init_flag = 1;
	}

	if((dim = element_dim(element.type)) <= 1){
		fprintf(stderr, "element_dim() < ");
		return(ARG_ERROR_RC);
	}
	if(gdim < dim) return(ARG_ERROR_RC);
	*vect_size = dim*element.node_num;

	RC_NEG_CHK( m_index = search_fem_material_prop_element(material,
	                                            physical, &element) );
	mat = material.array[m_index];
	RC_TRY( fem_fill_material(&mat) );

	/* temp_load を初期化 */
	/* index_array を初期化 */
	/* node_temp_array を初期化 */
	for(ii1=0; ii1<element.node_num; ii1++){
		int index;

		RC_NEG_CHK( index = node_renum_index(node, element.node[ii1]) );

		for(ii2=0; ii2<dim; ii2++){
			temp_load[dim*ii1 + ii2] = 0.0;
			index_array[dim*ii1 + ii2] = gdim*index + ii2;
		}
		node_temp_array[ii1] = node_temp(temp, def_temp,
		                                 element.node[ii1]);
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
		double local_force[3*FEM_MAX_NODE];
		double N[FEM_MAX_NODE];
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
static double node_temp(FEM_BC_ARRAY temp, FEM_DEFAULT_BC_ARRAY def_temp,
                                                      int node_label)
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

	index = search_fem_bc_node(temp, node_label);
	if(index < 0) return(d_temp);

	for(ii1=0; ii1<FEM_BC_SCALAR; ii1++){
		if(temp.array[index].s_type[ii1] == BC_TEMP){
			return(temp.array[ii1].s[ii1]);
		}
	}

	return(d_temp);
}


/* 全体剛性マトリックスを作成 */
RC fem_make_g_matrix(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                     FEM_MATERIAL_PROP_ARRAY material,
                     FEM_PHYSICAL_PROP_ARRAY physical,
                     FEM_SKYLINE_MATRIX *g_matrix)
{
	int ii1;
	int dim;
	FEM_ELEM_MATRIX elem_matrix;

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	/* 問題の次元(g_matrix->dim) */
	dim = fem_analysis_dim(element);
	if(dim <= 1){
		fprintf(stderr, "fem_analysis_dim < ");
		return(ARG_ERROR_RC);
	}

	/* g_matrix の確保と初期化 */
	RC_TRY( fem_allocate_g_matrix(dim, node, element, g_matrix) );

	for(ii1=0;ii1<element.size;ii1++){
		if(element.array[ii1].label < 0) continue;

		/* element.array[ii1], node, material, physical -> elem_matrix */
		RC_TRY( fem_make_elem_matrix(element.array[ii1], node, material,
		                             physical, &elem_matrix) );

		/* elem_matrix -> g_matrix */
		RC_TRY( fem_impose_elem_matrix(*g_matrix, elem_matrix, 1.0) );
		RC_TRY( free_fem_elem_matrix(&elem_matrix) );
	}

	return(NORMAL_RC);
}


/* 全体剛性マトリックスの確保と問題の次元をチェック */
RC fem_allocate_g_matrix(int dim, FEM_NODE_ARRAY node,
                         FEM_ELEMENT_ARRAY element,
                         FEM_SKYLINE_MATRIX *g_matrix)
{
	g_matrix->dim = dim;

	g_matrix->dof = (g_matrix->dim)*node.size;
	g_matrix->index1 = (long *)malloc(g_matrix->dof * sizeof(long));
	if((g_matrix->index1) == NULL){
		fprintf(stderr, "malloc < ");
		return(ALLOC_ERROR_RC);
	}
	RC_TRY( fill_index1(g_matrix->dim, node, element, g_matrix->index1) );

	RC_TRY( s_cholesky_alloc(g_matrix->dof, g_matrix->index1,
	                         &(g_matrix->index2), &(g_matrix->array),
	                         &(g_matrix->array_size)) );

	return(NORMAL_RC);
}


/* index1 を適切な値で初期化する */
/* index1[i] ...i 列目の最初の非ゼロ行番号 */
static RC fill_index1(int dim, FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                      long *index1)
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
			if(pos1 < 0){
				fprintf(stderr, "node_renum_index() < ");
				return(SEEK_ERROR_RC);
			}

			for(ii3=0; ii3<element.array[ii1].node_num; ii3++){
				pos2 = dim*node_renum_index(node, element.array[ii1].node[ii3]);
				if(pos2 < 0){
					fprintf(stderr, "node_renum_index() < ");
					return(SEEK_ERROR_RC);
				}

				for(ii4=0; ii4<dim; ii4++){
					if(pos1 < index1[pos2+ii4]) index1[pos2+ii4] = pos1;
				}
			}
		}
	}

	return(NORMAL_RC);
}


/* 要素剛性マトリックスの作成 */
RC fem_make_elem_matrix(FEM_ELEMENT element, FEM_NODE_ARRAY node,
                        FEM_MATERIAL_PROP_ARRAY material,
                        FEM_PHYSICAL_PROP_ARRAY physical,
                        FEM_ELEM_MATRIX *elem_matrix)
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
		RC_TRY( allocate2D(3*FEM_MAX_NODE, 3*FEM_MAX_NODE, &BtDB) );
		RC_TRY( allocate2D(FEM_MAX_NODE, 3, &node_location) );
		RC_TRY( allocate_g_prop(&g_prop) );
		init_flag = 1;
	}

	if((elem_matrix->dim = element_dim(element.type)) <= 0){
		fprintf(stderr, "element_dim < ");
		return(ARG_ERROR_RC);
	}

	if(element.label < 0) return(ARG_ERROR_RC);
	elem_matrix->element = element.label;
	elem_matrix->size = (element.node_num) * (elem_matrix->dim);

	/* elem_matrix->array[][] の確保 */
	RC_TRY( allocate2D(elem_matrix->size,
	                   elem_matrix->size, &(elem_matrix->matrix)) );

	/* elem_matirx->index[] の確保と代入 */
	RC_TRY( rc_malloc(sizeof(int)*(elem_matrix->size),
	                  (void**)(&elem_matrix->index)) );
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


/* 要素剛性マトリックス(熱伝導解析用) */
RC fem_make_elem_matrix_therm(FEM_ELEMENT element, FEM_NODE_ARRAY node,
                              FEM_MATERIAL_PROP_ARRAY material,
                              FEM_PHYSICAL_PROP_ARRAY physical,
                              FEM_ELEM_MATRIX *elem_matrix)
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
		RC_TRY( allocate2D(FEM_MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(MAX_GAUSS_POINT, 4, &gauss_points) );
		RC_TRY( allocate2D(3, FEM_MAX_NODE, &dN) );
		RC_TRY( allocate2D(3, FEM_MAX_NODE, &dNxyz) );
		RC_TRY( allocate2D(3, 3, &jacobi) );
		RC_TRY( allocate2D(3, 3, &inverse_J) );
		RC_TRY( allocate2D(3, 3, &D_matrix) );
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &BtDB) );

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
	RC_TRY( rc_malloc(sizeof(int)*(elem_matrix->size),
	                  (void**)(&elem_matrix->index)) );
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
RC get_spring_const(FEM_PHYSICAL_PROP_ARRAY physical, int p_label,
                    double *spring_const)
{
	int p_index;

	if(physical.source == FEM_NASTRAN){
		RC_NEG_CHK( p_index
		            = search_fem_physical_prop_label(physical, p_label) );
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
RC get_thickness(FEM_ELEMENT element,
                 FEM_PHYSICAL_PROP_ARRAY physical, double *thickness)
{
	int dim;
	int p_index;

	if((dim = element_dim(element.type)) <= 0){
		fprintf(stderr, "element_dim() < ");
		return(ARG_ERROR_RC);
	}

	/* 厚さ */
	*thickness = 1.0;
	if((dim == 2)&&(physical.source == FEM_NASTRAN)){
		RC_NEG_CHK( p_index = search_fem_physical_prop_label(physical,
		                                             element.physical) );
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
RC get_D_matrix(FEM_ELEMENT element, FEM_MATERIAL_PROP_ARRAY material,
                FEM_PHYSICAL_PROP_ARRAY physical, double **D_matrix)
{
	int ii1, ii2;
	int dim, ssnum;
	int m_index;
	FEM_MATERIAL_PROP tmp_mat;

	RC_NEG_CHK( m_index = search_fem_material_prop_element(material,
	                                            physical, &element) );

	if((dim = element_dim(element.type)) <= 0){
		fprintf(stderr, "element_dim < ");
		return(ARG_ERROR_RC);
	}

	if( (ssnum = element_ssnum(element.type)) <= 0){
		fprintf(stderr, "element_ssnum < ");
		return(ARG_ERROR_RC);
	}

	tmp_mat = material.array[m_index];
	if(dim == 2){
		tmp_mat.ana_type = DEF_ANA_2D;
	}else{   /* dim == 3 */
		tmp_mat.ana_type = ANA_3D;
	}
	RC_TRY( fem_fill_material(&tmp_mat) );

	/* tmp_mat.D_matrix -> D_matrix */
	for(ii1=0;ii1<ssnum;ii1++){
		for(ii2=0;ii2<ssnum;ii2++){
			D_matrix[ii1][ii2] = tmp_mat.D_matrix[ii1][ii2];
		}
	}

	return(NORMAL_RC);
}


/* D マトリックス（熱伝導解析用） */
RC get_D_matrix_therm(FEM_ELEMENT element, FEM_MATERIAL_PROP_ARRAY material,
                      FEM_PHYSICAL_PROP_ARRAY physical, double **D_matrix)
{
	int ii1, ii2;
	int dim;
	int m_index;

	RC_NEG_CHK( m_index = search_fem_material_prop_element(material,
	                                            physical, &element) );

	if((dim = element_dim(element.type)) <= 0){
		fprintf(stderr, "element_dim < ");
		return(ARG_ERROR_RC);
	}

	if(material.array[m_index].mat_type == MAT_ISOTROPIC){
		for(ii1=0;ii1<dim;ii1++){
			for(ii2=0;ii2<dim;ii2++){
				D_matrix[ii1][ii2] = 0.0;
			}
			D_matrix[ii1][ii1] = material.array[m_index].K;
		}
	}else{
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


/* 次数低減積分の処理 (B_matrix を上書き) */
static RC RI_B_matrix(FEM_ELEMENT elem, G_PROP g_prop)
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
static RC set_g_prop(FEM_ELEMENT element, double **node_location,
                     G_PROP *g_prop)
{
	int ii1;
	int dim;

	if( (g_prop->B_matrix_row = element_ssnum(element.type)) <= 0){
		fprintf(stderr, "element_ssnum < ");
		return(ARG_ERROR_RC);
	}

	if((dim = element_dim(element.type)) <= 0){
		fprintf(stderr, "element_dim < ");
		return(ARG_ERROR_RC);
	}

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
RC make_B_matrix(FEM_ELEMENT element, const double *g_point,
                 double **node_location, double **B_matrix, double *det_J)
{
	static int allocate_flag = 0;
	static double **dN = NULL;
	static double **dNxyz = NULL;
	static double **Jacobi = NULL;
	static double **inverse_J = NULL;
	int dim;

	if(allocate_flag == 0){
		RC_TRY( allocate2D(3, FEM_MAX_NODE, &dN) );
		RC_TRY( allocate2D(3, FEM_MAX_NODE, &dNxyz) );
		RC_TRY( allocate2D(3, 3, &Jacobi) );
		RC_TRY( allocate2D(3, 3, &inverse_J) );

		allocate_flag = 1;
	}

	if((dim = element_dim(element.type)) <= 0){
		fprintf(stderr, "element_dim < ");
		return(ARG_ERROR_RC);
	}

	/* g_point -> dN */
	/* dN * node_location => Jacobi */
	/* Jacobi^-1 => inverse_J */
	/* inverse_J * dN => dNxy */
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
static RC set_Gauss_const(FEM_ELEM_TYPE type, G_PROP *g_prop)
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
/* dNxyz[0][i]...Ni の x 方向微分 */
/* dNxyz[1][i]...Ni の y 方向微分 */
/* dNxyz[2][i]...Ni の z 方向微分 */
/* dim=2 ...2次元, dim=3 ...3次元 */
static RC set_B_matrix(int dim, int node_num,
                       double **dNxyz, double **B_matrix)
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
RC fem_impose_elem_matrix(FEM_SKYLINE_MATRIX g_matrix,
                          FEM_ELEM_MATRIX elem_matrix, double fact)
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
			if(rc != NORMAL_RC){
				fprintf(stderr, "s_cholesky_ptr < ");
				return(rc);
			}
			(*ptr) += fact * elem_matrix.matrix[ii1][ii2];
		}
	}

	return(NORMAL_RC);
}


RC fem_impose_nonzero3(NONZERO_MATRIX3 matrix,
                       FEM_ELEM_MATRIX elem_matrix, double fact)
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
			if(col_index < 0){
				fprintf(stderr, "col_index < ");
				return(ARG_ERROR_RC);
			}

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

RC fem_iccg_nonzero3(NONZERO_MATRIX3 matrix, const double *f_vect, 
                     double **d_vect, const double *init_vect)
{
	int ii1;

	*d_vect = (double *)malloc(sizeof(double) * 3 * matrix.size);
	if(*d_vect == NULL){
		fprintf(stderr, "malloc() <\n");
		return(ALLOC_ERROR_RC);
	}

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

RC fem_scale_iccg_nonzero3(NONZERO_MATRIX3 matrix, double *f_vect, 
                           double **d_vect, const double *init_vect)
{
	int ii1;
	double *scale_fact;

	*d_vect = (double *)malloc(sizeof(double) * 3 * matrix.size);
	scale_fact = (double *)malloc(sizeof(double) * 3 * matrix.size);
	if((*d_vect == NULL)||(scale_fact == NULL)){
		fprintf(stderr, "malloc() <\n");
		return(ALLOC_ERROR_RC);
	}

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

	free(scale_fact);

	return(NORMAL_RC);
}


RC fem_rest2vect(int dim, FEM_NODE_ARRAY node, FEM_BC_ARRAY rest,
                 double **d_vect)
{
	int ii1;
	int index;

	if((dim != 2)&&(dim != 3)&&(dim != 6)) return(ARG_ERROR_RC);

	RC_TRY( rc_malloc((node.size) * dim * sizeof(double), (void **)d_vect) );
	for(ii1=0; ii1<dim*(node.size); ii1++){
		(*d_vect)[ii1] = 0.0;
	}

	for(ii1=0; ii1<rest.size; ii1++){
		if(rest.array[ii1].node < 0) continue;

		RC_NEG_CHK( index = node_renum_index(node, rest.array[ii1].node) );

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


RC fem_disp2vect(int dim, FEM_NODE_ARRAY node, FEM_DISP_ARRAY disp,
                 double **d_vect)
{
	int ii1;
	int index;

	if((dim != 2)&&(dim != 3)&&(dim != 6)) return(ARG_ERROR_RC);

	RC_TRY( rc_malloc((node.size) * dim * sizeof(double), (void **)d_vect) );
	for(ii1=0; ii1<dim*(node.size); ii1++){
		(*d_vect)[ii1] = 0.0;
	}

	for(ii1=0; ii1<disp.size; ii1++){
		if(disp.array[ii1].node < 0) continue;

		RC_NEG_CHK( index = node_renum_index(node, disp.array[ii1].node) );

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


RC fem_vect2disp(int dim, FEM_NODE_ARRAY node, const double *d_vect,
                 FEM_DISP_ARRAY *disp)
{
	int ii1;
	int index;

	if((dim != 2)&&(dim != 3)&&(dim != 6)) return(ARG_ERROR_RC);

	RC_TRY( allocate_fem_disp_array(node.size, disp) );
	disp->source = FEM_K_SOL;
	for(ii1=0; ii1<node.size; ii1++){
		if(node.array[ii1].label < 0){
			disp->array[ii1].node = -1;
			continue;
		}
		index = node.array[ii1].renum;
		if((index < 0)||(index >= node.size)){
			fprintf(stderr, "node.array[index].renum < ");
			return(ARG_ERROR_RC);
		}

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

	RC_TRY( clean_fem_disp_array(disp) );

	return(NORMAL_RC);
}


RC fem_vect2react(int dim, FEM_NODE_ARRAY node, const double *r_vect,
                  FEM_REACT_ARRAY *react)
{
	int ii1;
	int index;

	if((dim != 2)&&(dim != 3)) return(ARG_ERROR_RC);

	RC_TRY( allocate_fem_react_array(node.size, react) );
	react->source = FEM_K_SOL;
	for(ii1=0; ii1<node.size; ii1++){
		if(node.array[ii1].label < 0){
			react->array[ii1].node = -1;
			continue;
		}
		index = node.array[ii1].renum;
		if((index < 0)||(index >= node.size)){
			fprintf(stderr, "node.array[index].renum < ");
			return(ARG_ERROR_RC);
		}

		react->array[ii1].node = node.array[ii1].label;
		react->array[ii1].v.x = r_vect[dim*index];
		react->array[ii1].v.y = r_vect[dim*index + 1];
		if(dim == 3) react->array[ii1].v.z = r_vect[dim*index + 2];
	}
	react->size = node.size;

	RC_TRY( clean_fem_react_array(react) );

	return(NORMAL_RC);
}


RC fem_modify_nonzero3(FEM_NODE_ARRAY node, FEM_BC_ARRAY rest,
                       NONZERO_MATRIX3 matrix, double *f_vect)
{
	return( fem_modify_nonzero3_n(node, rest, matrix, &f_vect, 1) );
}


RC fem_modify_nonzero3_n(FEM_NODE_ARRAY node, FEM_BC_ARRAY rest,
                         NONZERO_MATRIX3 matrix, double **f_vect,
                         int vect_size)
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
RC fem_allocate_nonzero3(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                         NONZERO_MATRIX3 *matrix)
{
	int ii1, ii2;
	FEM_ELEM_LOOKUP_ARRAY elem_lookup;
	int *row_nodes;
	int row_size;
	int syn_size;
	int index;

	if(fem_analysis_dim(element) != 3) return(ARG_ERROR_RC);
	if(node.renum_flag != RENUMBERED){
		RC_TRY( dummy_renumber(&node) );
	}

	/* 作業用配列 */
	RC_NULL_CHK( row_nodes = (int *)malloc(node.size*sizeof(int)) );

	RC_TRY( make_elem_lookup(node, element, &elem_lookup) );

	matrix->size = node.size;
	matrix->row = (NONZERO_ROW3 *)malloc(matrix->size*sizeof(NONZERO_ROW3));
	if(matrix->row == NULL){
		fprintf(stderr, "malloc() < ");
		return(ALLOC_ERROR_RC);
	}
	for(ii1=0; ii1<matrix->size; ii1++){
		matrix->row[ii1].size = 0;
		matrix->row[ii1].elem = NULL;
		matrix->row[ii1].syn_size = 0;
		matrix->row[ii1].syn = NULL;
	}

	for(ii1=0; ii1<elem_lookup.size; ii1++){
		if(elem_lookup.array[ii1].node < 0) continue;

		RC_NEG_CHK( index = node_renum_index(node,
			                 elem_lookup.array[ii1].node) );

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
			RC_NULL_CHK( matrix->row[index].syn
			             = malloc(syn_size*sizeof(int)) );
		}
		for(ii2=0; ii2<syn_size; ii2++){
			matrix->row[index].syn[ii2] = row_nodes[row_size - syn_size + ii2];
		}

		row_size -= syn_size;
		matrix->row[index].size = row_size;
		RC_NULL_CHK( matrix->row[index].elem
		             = malloc(row_size*sizeof(NONZERO_ELEM3)) );
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

	free(row_nodes);
	RC_TRY( free_fem_elem_lookup_array(&elem_lookup) );

	return(NORMAL_RC);
}


/* 該当行に関連する列のインデックスを row_nodes に登録 */
static RC fill_row_nodes(FEM_ELEM_LOOKUP elem_lookup,
                         FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                         int row_nodes[], int *row_size)
{
	int ii1, ii2, ii3;
	int flag, col_index;
	int elem_index;

	(*row_size) = 0;
	for(ii1=0; ii1<elem_lookup.size; ii1++){
		elem_index = search_fem_element_label(element,
		                                      elem_lookup.element[ii1]);
		RC_NEG_CHK(elem_index);

		for(ii2=0; ii2<element.array[elem_index].node_num; ii2++){
			RC_NEG_CHK( col_index = node_renum_index(node,
			                  element.array[elem_index].node[ii2]) );

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
RC fem_make_mass_matrix(FEM_ELEMENT element, FEM_NODE_ARRAY node,
                        FEM_MATERIAL_PROP_ARRAY material,
                        FEM_PHYSICAL_PROP_ARRAY physical,
                        FEM_ELEM_MATRIX *elem_matrix)
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
		RC_TRY( allocate2D(FEM_MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(MAX_GAUSS_POINT, 4, &gauss_points) );
		RC_TRY( allocate2D(3, FEM_MAX_NODE, &dN) );
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
	if(elem_matrix->dim <= 0){
		fprintf(stderr, "element_dim() < ");
		return(ARG_ERROR_RC);
	}
	elem_matrix->size = (element.node_num) * (elem_matrix->dim);

	/* elem_matrix->array[][] の確保 */
	RC_TRY( allocate2D(elem_matrix->size,
	                   elem_matrix->size, &(elem_matrix->matrix)) );

	/* elem_matirx->index[] の確保と代入 */
	RC_TRY( rc_malloc(sizeof(int)*(elem_matrix->size),
	                  (void**)(&elem_matrix->index)) );
	for(ii1=0; ii1<element.node_num; ii1++){
		int index = node_renum_index(node, element.node[ii1]);
		if(index < 0) return(SEEK_ERROR_RC);

		for(ii2=0; ii2<(elem_matrix->dim); ii2++){
			elem_matrix->index[ii1*(elem_matrix->dim)+ii2]
			             = index*(elem_matrix->dim) + ii2;
		}
	}

	/* rho, rho_matrix */
	RC_NEG_CHK( mat_index = search_fem_material_prop_element(material,
	                                              physical, &element) );
	rho = material.array[mat_index].rho;
	for(ii1=0; ii1<(elem_matrix->dim); ii1++){
		for(ii2=0; ii2<(elem_matrix->dim); ii2++){
			rho_matrix[ii1][ii2] = 0.0;
		}
		rho_matrix[ii1][ii1] = rho;
	}

	for(ii1=0; ii1<num_point; ii1++){
		double N[FEM_MAX_NODE];
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


static RC Gauss_const4mass_matrix(FEM_ELEM_TYPE type, double **points,
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
