/*********************************************************************
 * flow_solver.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Shoji ITO> <Tomoyuki TSUBATA>
 *
 * Thanks to following contributors.
 *  <Kenzen TAKEUCHI> <Takaaki NAGATANI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: flow_solver.c,v 1.12 2003/12/05 02:47:32 nagatani Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "avs_utl.h"
#include "fem_solver.h"
#include "fem_struct.h"
#include "math_utl.h"
#include "sky_cholesky.h"
#include "sky_crout.h"
#include "rc.h"

#define NR_REL_ERROR (DBL_EPSILON * 1.0e+2)  /* Newton-Raphson 用相対誤差 */
#define ABS_ERROR (10.0*DBL_MIN/DBL_EPSILON) /* 絶対誤差 */

typedef struct{
	int num_point;
	double **g_points;
	double N[MAX_GAUSS_POINT][FEM_MAX_NODE];
	double weight[MAX_GAUSS_POINT];
	double det_J[MAX_GAUSS_POINT];
	double **dNxy[MAX_GAUSS_POINT];
} G_POINT_PROP;


/* 流れ場解析------------------------------------------------------------ */
static RC solve_middle_velocity(double delta_t, double Re, 
                                FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element, 
                                FEM_BC_ARRAY bc_velo, FEM_VELO_ARRAY vp,
                                FEM_VELO_ARRAY *vp1);
static RC set_f_vector_for_pressure(int dim, double delta_t,
                                    FEM_NODE_ARRAY node,
                                    FEM_ELEMENT_ARRAY element,
                                    FEM_VELO_ARRAY vp1, double *f_vector);
static RC set_elem_matrix_for_pressure(FEM_ELEMENT element, FEM_NODE_ARRAY node,
                                       int *index_array, int *matrix_size,
                                       double **elem_matrix);
static RC solve_pressure(int dim, double delta_t, FEM_NODE_ARRAY node,
                         FEM_ELEMENT_ARRAY element, FEM_VELO_ARRAY *vp1, 
                         FEM_BC_ARRAY bc);
static RC solve_velocity(double delta_t, FEM_NODE_ARRAY node,
                         FEM_ELEMENT_ARRAY element, FEM_BC_ARRAY bc_velo,
                         FEM_VELO_ARRAY vp1, FEM_VELO_ARRAY *vp);


/* 随伴解析---------------------------------------------------------------*/
static RC set_f_vector_for_adjoint_analysis(int dim, FEM_NODE_ARRAY node, 
                                            FEM_ELEMENT_ARRAY element,
                                            FEM_BC_ARRAY bc, FEM_VELO_ARRAY vp,
                                            double *f_vector);
static RC set_elem_matrix_for_adjoint_analysis(int dim, double Re,
                                               FEM_ELEMENT v_element,
                                               FEM_ELEMENT p_element,
                                               FEM_NODE_ARRAY node,
                                               FEM_VELO_ARRAY vp,
                                               int *index_array,
                                               int *matrix_size,
                                               double **elem_matrix);
static RC adjoint_modify_matrix(int dim, int matrix_size,
                                FEM_NODE_ARRAY node, FEM_BC_ARRAY bc,
                                double **array, double *vector);


/* 感度導出用ツール------------------------------------------------------ */
static RC local_velocity_strain(FEM_ELEMENT element, const double *local_point,
                                FEM_NODE_ARRAY node, FEM_VELO_ARRAY velo,
                                FEM_STRESS_STRAIN *strain_velo);
static RC fill_velo_vector(FEM_ELEMENT element, FEM_VELO_ARRAY velo,
                           double *velo_v);


/* 解析ツール------------------------------------------------------------ */
static RC allocate_g_prop(G_POINT_PROP *g_prop);
static RC set_g_prop(FEM_ELEMENT element, FEM_NODE_ARRAY node,
                     G_POINT_PROP *g_prop);
static RC allocate_matrix(int mode, int dim, FEM_NODE_ARRAY node,
                          FEM_ELEMENT_ARRAY element,
                          FEM_SKYLINE_MATRIX *g_matrix);
static RC flow_modify_matrix(int mode, FEM_NODE_ARRAY node, 
                             FEM_BC_ARRAY bc_press,
                             FEM_SKYLINE_MATRIX g_matrix, double *f_vector);
static RC fill_index1(int dim, FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                      long *index1);
static RC imp_g_matrix(int mode, int matrix_size, const int *index_array,
                       double **elem_matrix, FEM_SKYLINE_MATRIX g_matrix);
static RC set_Mass_matrix(int node_num, G_POINT_PROP g_prop, double **M);
static RC set_K_matrix(int dN_xyz, int v_xyz, int dim, FEM_NODE_ARRAY node,
                       FEM_ELEMENT element, G_POINT_PROP g_prop, double **V,
                       double **K);
static RC set_S_matrix(int dN_xyz1, int dN_xyz2, int dim, FEM_ELEMENT element, 
                       G_POINT_PROP g_prop, double **S);
static RC set_H_matrix(int dN_xyz, FEM_ELEMENT element, G_POINT_PROP g_prop,
                       double **H);
static RC set_H1_matrix(int dN_xyz, FEM_ELEMENT v_elem, FEM_ELEMENT p_elem,
                        G_POINT_PROP v_g_prop, G_POINT_PROP p_g_prop,
                        double **H1);
static RC set_H2_matrix(int dN_xyz, FEM_ELEMENT v_elem, FEM_ELEMENT p_elem,
                        G_POINT_PROP v_g_prop, G_POINT_PROP p_g_prop,
                        double **H2);
static RC set_advection_term(int dim, FEM_NODE_ARRAY node, FEM_ELEMENT element,
                             G_POINT_PROP g_prop, double **V, double **ADV);
static RC set_viscosity_term(int dim,  FEM_ELEMENT element, G_POINT_PROP g_prop,
                             double **V, double **VIS);
static RC make_press_element_data(FEM_ELEMENT_ARRAY velo_elem,
                                  FEM_ELEMENT_ARRAY *press_elem);
static RC pivotting_check(FEM_NODE_ARRAY node, FEM_SKYLINE_MATRIX g_matrix,
                          double *vector);
static RC set_Gauss_const(FEM_ELEM_TYPE type, G_POINT_PROP *g_prop);
static RC free_g_matrix(int mode, FEM_SKYLINE_MATRIX *g_matrix);
static RC allocate_g_prop(G_POINT_PROP *g_prop);
static RC flow_property_print(int dim, int iteration, double Re,
                              double delta_t);


static RC modify_velo_array(FEM_BC_ARRAY bc_velo, FEM_BC_ARRAY bc_press,
                            FEM_VELO_ARRAY *flow);
static RC modify_nonzero_matrix1a(FEM_BC_ARRAY bc_velo, FEM_BC_ARRAY bc_press,
                                  FEM_NODE_ARRAY node, int valid_node_num,
                                  NONZERO_MATRIX *matrix, double *global_vect); 
static RC modify_nonzero_matrix1a_diagonal(int diagonal_index, double value,
                                           NONZERO_MATRIX *matrix,
                                           double *global_vect);
static RC fem_impose_nonzero1a(FEM_ELEM_MATRIX elem_matrix,
                               const double *local_vect,
                               NONZERO_MATRIX *matrix, double *global_vect);
static RC make_elem_matrix(int dim, int model_order, FEM_ELEMENT elem,
                           FEM_NODE_ARRAY node, int valid_node_num,
                           int velo_order, int press_order,
                           double density, double viscosity,
                           const double *flow_vect, 
                           FEM_ELEM_MATRIX *elem_matrix, double **local_vect,
                           FEM_WORK_SET ws);
static RC fem_velo2vect(int dim, FEM_NODE_ARRAY node, FEM_VELO_ARRAY flow,
                        int valid_node_num, double *flow_vect);
static RC fem_vect2velo(int dim, FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY elem,
                        int velo_order, int press_order, int valid_node_num,
                        const double *flow_vect, FEM_VELO_ARRAY *flow);
static RC parameter_check(FEM_ELEMENT elem, int model_order, int para_order,
                          FEM_ELEM_TYPE *para_type, int *para_node_num);

static RC make_adjoint_elem_matrix(int dim, int model_order, FEM_ELEMENT elem,
                                   FEM_NODE_ARRAY node, int valid_node_num,
                                   int velo_order, int press_order,
                                   double density, double viscosity,
                                   const double *flow_vect,
                                   FEM_ELEM_MATRIX *elem_matrix,
                                   double **local_vect, FEM_WORK_SET ws);

/*
 * 定常流れ場の解析
 * 
 * 連立方程式の解法 Newton-Raphson 
 */
RC steady_flow_analysis(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY elem,
                        FEM_BC_ARRAY bc_velo, FEM_BC_ARRAY bc_press,
                        int velo_order, int press_order,
                        double density, double viscosity,
                        FEM_VELO_ARRAY *flow)
{
	int ii1;
	int dim;
	int model_order;
	int valid_node_num;

	NONZERO_MATRIX matrix;
	double *global_vect;
	double *delta;

	FEM_WORK_SET ws;

	FEM_ELEM_MATRIX elem_matrix;
	double *local_vect;
	double *flow_vect;
	double delta_max;
	double delta_norm;
	double flow_norm;
	int init_flag = 0;

	RC_TRY( allocate_fem_work_set(&ws) );

	dim = fem_analysis_dim(elem);
	if( (dim != 2) && (dim != 3) ){
		fprintf(stderr, "fem_analysis_dim <\n");
		return(IMPLEMENT_ERROR_RC);
	}

	model_order = fem_analysis_order(elem);

	if(model_order != 2){
		fprintf(stderr, "fem_model_order <\n");
		return(IMPLEMENT_ERROR_RC);
	}

	if( (velo_order < 1)
	 || (model_order < 1)
	 || (velo_order > model_order)
	 || (press_order > model_order) ){
		return(OVERFLOW_ERROR_RC); 
	}

	RC_TRY( locate_middle_node(node, elem) );
	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );
	valid_node_num = count_valid_node(node);

	/* 以前にこの関数で流速が計算されていたら確保しない */
	if( (flow->size != node.size) && (flow->source != FEM_K_SOL) ){
		init_flag = 1;
		/* FEM_VELO_ARRAY の確保 */
		RC_TRY( allocate_fem_velo_array(node.size, flow) );
		flow->source = FEM_K_SOL;
		flow->size = node.size;
		for(ii1=0; ii1<node.size; ii1++){
			if(node.array[ii1].label < 0){
				flow->array[ii1].node = -1;
				continue;
			}
			flow->array[ii1].node = node.array[ii1].label;
		}
		/* FEM_VELO_ARRAY に、境界条件の値を代入する */
		RC_TRY( modify_velo_array(bc_velo, bc_press, flow) );
	}

	/* FEM_VELO_ARRAY の計算用ベクトルの確保 -> 変換*/
	RC_TRY( allocate1D(valid_node_num*(dim+1), &flow_vect) );
	RC_TRY( fem_velo2vect(dim, node, *flow, valid_node_num, flow_vect) );

	/* Newton-Raphson 法による増分ベクトルの確保 */
	RC_TRY( allocate1D(valid_node_num*(dim+1), &delta) );

	while(1){
	fprintf(stderr, "Making Matrix........");
		RC_TRY( allocate_nonzero_matrix1a(valid_node_num*(dim+1), &matrix) );
		RC_TRY( allocate1D(valid_node_num*(dim+1), &global_vect) );

		for(ii1=0; ii1<elem.size; ii1++){
			if(elem.array[ii1].label < 0) continue;

			if(init_flag == 1){
				RC_TRY( make_elem_matrix( dim, model_order, elem.array[ii1],
				  	                      node, valid_node_num,
				                          velo_order, press_order,
				                          0.0, viscosity,
				                          flow_vect, &elem_matrix, &local_vect,
										  ws) );
			}else{
				RC_TRY( make_elem_matrix( dim, model_order, elem.array[ii1],
				  	                      node, valid_node_num,
				                          velo_order, press_order,
				                          density, viscosity,
				                          flow_vect, &elem_matrix, &local_vect,
                                          ws) );
			}

			RC_TRY( fem_impose_nonzero1a(elem_matrix, local_vect, &matrix,
			                             global_vect) );

			RC_TRY( free_fem_elem_matrix(&elem_matrix) );
			free( (void *)local_vect );
		}
		init_flag = 0;

		/* 連立方程式に境界条件を代入する */
		RC_TRY( modify_nonzero_matrix1a(bc_velo, bc_press, node,
		                                valid_node_num, &matrix, global_vect) );

		RC_TRY( nonzero1a_scaling(&matrix, global_vect) );
		fprintf(stderr, "OK\n");

		RC_TRY( cg_nonzero1a_ver1(matrix, global_vect, delta) );

		/* 流速、圧力に増分ベクトルを加える */
		for(ii1=0; ii1<matrix.size; ii1++){
			flow_vect[ii1] += delta[ii1];
		}

		delta_norm = 0.0;
		flow_norm = 0.0;
		for(ii1=0; ii1<matrix.size; ii1++){
			delta_norm += delta[ii1] * delta[ii1];
			flow_norm += flow_vect[ii1] * flow_vect[ii1];
		}
		delta_max = sqrt( (delta_norm/flow_norm)/matrix.size);
		fprintf(stderr, "   norm = %e\n", delta_max);

		RC_TRY( free_nonzero_matrix1a(&matrix) );
		free( (void *)global_vect );

		/* Newton-Raphson の収束 */
		if(delta_max < NR_REL_ERROR) break;

		/* Stokes 流れは線形連立方程式 */
		if(density < ABS_ERROR) break;
	}

	RC_TRY( fem_vect2velo(dim, node, elem, velo_order, press_order,
	                      valid_node_num, flow_vect, flow) );

	RC_TRY( free_fem_work_set(&ws) );
	free( (void *)flow_vect );
	free( (void *)delta );

	return(NORMAL_RC);
}


static RC modify_velo_array(FEM_BC_ARRAY bc_velo, FEM_BC_ARRAY bc_press,
                            FEM_VELO_ARRAY *flow)
{
	int index;
	int ii1;

	for(ii1=0; ii1<bc_press.size; ii1++){
		if(bc_press.array[ii1].node < 0) continue;

		index = search_fem_velo_node(*flow, bc_press.array[ii1].node);
		RC_NEG_CHK(index);

		if(bc_press.array[ii1].s_type[0] == BC_PRESS){
			flow->array[index].s = bc_press.array[ii1].s[0];
		}
	}

	for(ii1=0; ii1<bc_velo.size; ii1++){
		if(bc_velo.array[ii1].node < 0) continue;

		index = search_fem_velo_node(*flow, bc_velo.array[ii1].node);
		RC_NEG_CHK(index);

		if(bc_velo.array[ii1].v_type.x == BC_VELO){
			flow->array[index].v.x = bc_velo.array[ii1].v.x;
		}
		if(bc_velo.array[ii1].v_type.y == BC_VELO){
			flow->array[index].v.y = bc_velo.array[ii1].v.y;
		}
		if(bc_velo.array[ii1].v_type.z == BC_VELO){
			flow->array[index].v.z = bc_velo.array[ii1].v.z;
		}
	}

	return(NORMAL_RC);
}


static RC modify_nonzero_matrix1a(FEM_BC_ARRAY bc_velo, FEM_BC_ARRAY bc_press,
                                  FEM_NODE_ARRAY node, int valid_node_num,
                                  NONZERO_MATRIX *matrix, double *global_vect)
{
	int tmp;
	int index;
	int ii1;

	for(ii1=0; ii1<bc_press.size; ii1++){
		if(bc_press.array[ii1].node < 0) continue;

		index = search_fem_node_label(node, bc_press.array[ii1].node);
		RC_NEG_CHK(index);
		tmp = node.array[index].renum;

		if(bc_press.array[ii1].s_type[0] == BC_PRESS){
			RC_TRY( modify_nonzero_matrix1a_diagonal(valid_node_num*0 + tmp,
			                                         0, matrix, global_vect) );
		}
	}

	for(ii1=0; ii1<bc_velo.size; ii1++){
		if(bc_velo.array[ii1].node < 0) continue;

		index = search_fem_node_label(node, bc_velo.array[ii1].node);
		RC_NEG_CHK(index);
		tmp = node.array[index].renum;

		if(bc_velo.array[ii1].v_type.x == BC_VELO){
			RC_TRY( modify_nonzero_matrix1a_diagonal(valid_node_num*1 + tmp,
			                                         0, matrix, global_vect) );
		}
		if(bc_velo.array[ii1].v_type.y == BC_VELO){
			RC_TRY( modify_nonzero_matrix1a_diagonal(valid_node_num*2 + tmp,
			                                         0, matrix, global_vect) );
		}
		if(bc_velo.array[ii1].v_type.z == BC_VELO){
			RC_TRY( modify_nonzero_matrix1a_diagonal(valid_node_num*3 + tmp,
			                                         0, matrix, global_vect) );
		}
	}	

	return(NORMAL_RC);
}


/* B.C. より全体マトリックス、右辺ベクトルを修正 */
static RC modify_nonzero_matrix1a_diagonal(int diagonal_index, double value, 
                                           NONZERO_MATRIX *matrix,
                                           double *global_vect)
{
	int ii1, ii2;

	for(ii1=0; ii1<matrix->size; ii1++){
		for(ii2=0; ii2<matrix->row[ii1].size; ii2++){
			if(diagonal_index < matrix->row[ii1].col[ii2]) break;
			if(diagonal_index == matrix->row[ii1].col[ii2]){
				global_vect[ii1] -= matrix->row[ii1].v[ii2] * value;
				matrix->row[ii1].v[ii2] = 0.0;
			}
		}
	}

	/* 列に対してメモリ領域を1つにして1.0を代入 */
	global_vect[diagonal_index] = value;

	matrix->row[diagonal_index].size = 1;
	matrix->row[diagonal_index].col
	    = (int *)realloc(matrix->row[diagonal_index].col,
				         matrix->row[diagonal_index].size*sizeof(int) );
	RC_NULL_CHK(matrix->row[diagonal_index].col);
	matrix->row[diagonal_index].v
	    = (double *)realloc(matrix->row[diagonal_index].v,
	                        matrix->row[diagonal_index].size*sizeof(double) );
	RC_NULL_CHK(matrix->row[diagonal_index].v);
	matrix->row[diagonal_index].col[matrix->row[diagonal_index].size - 1]
	    = diagonal_index;
	matrix->row[diagonal_index].v[matrix->row[diagonal_index].size - 1] = 1.0;


	return(NORMAL_RC);
}


static RC fem_impose_nonzero1a(FEM_ELEM_MATRIX elem_matrix,
                               const double *local_vect,
                               NONZERO_MATRIX *matrix, double *global_vect)
{
	double *ptr;
	int ii1, ii2;
	int tmp1, tmp2;

	for(ii1=0; ii1<elem_matrix.size; ii1++){
		tmp1 = elem_matrix.index[ii1];
		if(tmp1 >= matrix->size) return(ARG_ERROR_RC);

		for(ii2=0; ii2<elem_matrix.size; ii2++){
			tmp2 = elem_matrix.index[ii2];
			if(tmp2 >= matrix->size) return(ARG_ERROR_RC);

			/* 全体マトリックスへの足し合わせ */
			RC_NULL_CHK( ptr = nonzero_matrix1a_ptr(matrix, tmp1, tmp2) );
			*ptr += elem_matrix.matrix[ii1][ii2];
		}
		/* 右辺ベクトルへの足し合わせ */
		global_vect[tmp1] += local_vect[ii1];
	}

	return(NORMAL_RC);
}


static RC make_elem_matrix(int dim, int model_order, FEM_ELEMENT elem,
                           FEM_NODE_ARRAY node, int valid_node_num,
                           int velo_order, int press_order,
                           double density, double viscosity,
                           const double *flow_vect, 
                           FEM_ELEM_MATRIX *elem_matrix, double **local_vect,
                           FEM_WORK_SET ws)
{
	static int allocate_flag = 0;
	static double ***Kx_tensor = NULL;
	static double ***Ky_tensor = NULL;
	static double ***Kz_tensor = NULL;
	static double **Sxx_matrix = NULL;
	static double **Syy_matrix = NULL;
	static double **Szz_matrix = NULL;
	static double **Sxy_matrix = NULL;
	static double **Syz_matrix = NULL;
	static double **Szx_matrix = NULL;
	static double **Hx_matrix = NULL;
	static double **Hy_matrix = NULL;
	static double **Hz_matrix = NULL;

	double **node_location = ws.mat_N_3[0];
	double **velo_dN = ws.mat_3_N[0];
	double **velo_dNxyz = ws.mat_3_N[1];
	double **Jacobi = ws.mat_3_3[0];
	double **inverse_J = ws.mat_3_3[1];
	double **g_points = ws.mat_G_4[0];
	double *velo_N = ws.vec_N[0];
	double *press_N = ws.vec_N[1];
	double *weight = ws.vec_G[0];
	FEM_ELEM_TYPE velo_type;
	FEM_ELEM_TYPE press_type;
	FEM_ELEM_TYPE base_type;
	double *elem_flow_vect;
	double det_J;
	int velo_node_num;
	int press_node_num;
	int base_node_num;
	int g_point_num;
	int index;
	int ii1, ii2, ii3, ii4;

	double tmp_x_1_x;
	double tmp_y_1_y;
	double tmp_z_1_z;

	double tmp_x_2_x;
	double tmp_x_2_y;
	double tmp_x_2_z;

	double tmp_y_2_x;
	double tmp_y_2_y;
	double tmp_y_2_z;

	double tmp_z_2_x;
	double tmp_z_2_y;
	double tmp_z_2_z;

	if(allocate_flag == 0){
		RC_TRY( allocate3D(FEM_MAX_NODE, FEM_MAX_NODE, FEM_MAX_NODE,
		                   &Kx_tensor) );
		RC_TRY( allocate3D(FEM_MAX_NODE, FEM_MAX_NODE, FEM_MAX_NODE,
		                   &Ky_tensor) );

		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &Sxx_matrix) );
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &Sxy_matrix) );
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &Syy_matrix) );
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &Hx_matrix) );
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &Hy_matrix) );

		if(dim == 3){
			RC_TRY( allocate3D(FEM_MAX_NODE, FEM_MAX_NODE, FEM_MAX_NODE,
		                       &Kz_tensor) );
			RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &Szz_matrix) );
			RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &Syz_matrix) );
			RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &Szx_matrix) );
			RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &Hz_matrix) );
		}

		allocate_flag = 1;
	}

	/* 初期化 */
	for(ii1=0; ii1<FEM_MAX_NODE; ii1++){
		 init_matrix(FEM_MAX_NODE, FEM_MAX_NODE, Kx_tensor[ii1]);
		 init_matrix(FEM_MAX_NODE, FEM_MAX_NODE, Ky_tensor[ii1]);
	}
	init_matrix(FEM_MAX_NODE, FEM_MAX_NODE, Sxx_matrix);
	init_matrix(FEM_MAX_NODE, FEM_MAX_NODE, Syy_matrix);
	init_matrix(FEM_MAX_NODE, FEM_MAX_NODE, Sxy_matrix);
	init_matrix(FEM_MAX_NODE, FEM_MAX_NODE, Hx_matrix);
	init_matrix(FEM_MAX_NODE, FEM_MAX_NODE, Hy_matrix);
	if(dim == 3){
		init_matrix(FEM_MAX_NODE, FEM_MAX_NODE, Szz_matrix);
		init_matrix(FEM_MAX_NODE, FEM_MAX_NODE, Syz_matrix);
		init_matrix(FEM_MAX_NODE, FEM_MAX_NODE, Szx_matrix);
		init_matrix(FEM_MAX_NODE, FEM_MAX_NODE, Hz_matrix);
		for(ii1=0; ii1<FEM_MAX_NODE; ii1++){
			init_matrix(FEM_MAX_NODE, FEM_MAX_NODE, Kz_tensor[ii1]);
		}
	}

	base_type = elem.type;
	base_node_num = elem.node_num;
	RC_TRY( parameter_check(elem, model_order, velo_order, &velo_type,
	                                                       &velo_node_num) );
	RC_TRY( parameter_check(elem, model_order, press_order, &press_type,
	                                                        &press_node_num) );
	RC_TRY( Gauss_const_max(base_type, g_points, weight, &g_point_num) );
	RC_TRY( node_location_matrix(&elem, node, node_location) );


	for(ii1=0; ii1<g_point_num; ii1++){
			
		/* 流速 */
		RC_TRY( set_N_vector(velo_type, g_points[ii1], velo_N) );
		RC_TRY( set_dN_matrix(velo_type, g_points[ii1], velo_dN) );
		mul_matrix(dim, velo_node_num, dim, velo_dN, node_location, Jacobi);
		det_J = determinant(dim, Jacobi);

		if(det_J < 0){
			fprintf(stderr, "Warning! Jacobian is negative!!!\n");
			return(NEGATIVE_ERROR_RC);
		}
		RC_TRY( inverse_matrix(dim, Jacobi, inverse_J) );
		mul_matrix(dim, dim, velo_node_num, inverse_J, velo_dN, velo_dNxyz);

		/* 圧力 */
		RC_TRY( set_N_vector(press_type, g_points[ii1], press_N) );

		for(ii2=0; ii2<velo_node_num; ii2++){
			for(ii3=0; ii3<velo_node_num; ii3++){

				/* S_matrix[velo_node_num][velo_node_num] */
				Sxx_matrix[ii2][ii3] += velo_dNxyz[0][ii2] * velo_dNxyz[0][ii3]
				                      * det_J * weight[ii1] * viscosity; 
				Syy_matrix[ii2][ii3] += velo_dNxyz[1][ii2] * velo_dNxyz[1][ii3]
				                      * det_J * weight[ii1] * viscosity;
				Sxy_matrix[ii2][ii3] += velo_dNxyz[0][ii2] * velo_dNxyz[1][ii3]
				                      * det_J * weight[ii1] * viscosity;
				if(dim == 3){
					Szz_matrix[ii2][ii3] += velo_dNxyz[2][ii2]
					                      * velo_dNxyz[2][ii3]
					                      * det_J * weight[ii1] * viscosity;
					Syz_matrix[ii2][ii3] += velo_dNxyz[1][ii2]
					                      * velo_dNxyz[2][ii3]
				                          * det_J * weight[ii1] * viscosity;
					Szx_matrix[ii2][ii3] += velo_dNxyz[2][ii2]
					                      * velo_dNxyz[0][ii3]
				                          * det_J * weight[ii1] * viscosity;
				}

				for(ii4=0; ii4<velo_node_num; ii4++){
					/* K_tensor[velo_node_num][velo_node_num][velo_node_num] */
					Kx_tensor[ii2][ii3][ii4] += velo_N[ii2] * velo_N[ii3]
					                          * velo_dNxyz[0][ii4] * det_J
					                          * weight[ii1] * density;
					Ky_tensor[ii2][ii3][ii4] += velo_N[ii2] * velo_N[ii3]
					                          * velo_dNxyz[1][ii4] * det_J
					                          * weight[ii1] * density;
					if(dim == 3){
						Kz_tensor[ii2][ii3][ii4] += velo_N[ii2] * velo_N[ii3]
					                              * velo_dNxyz[2][ii4] * det_J
					                              * weight[ii1] * density;
					}
				}

			}

			for(ii3=0; ii3<press_node_num; ii3++){
				/* H_matrix[velo_node_num][press_node_num] */
				Hx_matrix[ii2][ii3] += velo_dNxyz[0][ii2] * press_N[ii3]
				                     * det_J * weight[ii1];
				Hy_matrix[ii2][ii3] += velo_dNxyz[1][ii2] * press_N[ii3]
				                     * det_J * weight[ii1];
				if(dim == 3){
					Hz_matrix[ii2][ii3] += velo_dNxyz[2][ii2] * press_N[ii3]
				                         * det_J * weight[ii1];
				}

			}
		}

	}

	/* elem_matrix の確保 */
	elem_matrix->dim = dim;
	elem_matrix->element = elem.label;
	elem_matrix->size = base_node_num * (elem_matrix->dim + 1);
	elem_matrix->index = malloc( elem_matrix->size * sizeof(int) );
	RC_NULL_CHK(elem_matrix->index);
    RC_TRY( allocate2D( elem_matrix->size, elem_matrix->size,
	                    &(elem_matrix->matrix) ) );

	RC_TRY( allocate1D(elem_matrix->size, &elem_flow_vect) );

	/* elem_matrix.index[] の作成 */
	for(ii1=0; ii1<elem.node_num; ii1++){
		index = node_renum_index(node, elem.node[ii1]);
		RC_NEG_CHK(index);
		for(ii2=0; ii2<(elem_matrix->dim)+1; ii2++){
			elem_matrix->index[ii2*base_node_num + ii1]
        	= ii2 * valid_node_num + index;
			elem_flow_vect[ii2*base_node_num + ii1]
			= flow_vect[ii2*valid_node_num + index];
		}
	}

	/* 要素マトリックス用の右辺ベクトル作成 */
	RC_TRY( allocate1D(elem_matrix->size, local_vect) );

	for(ii1=0; ii1<velo_node_num; ii1++){
		for(ii2=0; ii2<velo_node_num; ii2++){

			/* tmp_x_1_x += tensor_x[][loop][] * velo_x[loop] */
			tmp_x_1_x = 0.0;
			tmp_y_1_y = 0.0;
			tmp_z_1_z = 0.0;

			tmp_x_2_x = 0.0;
			tmp_x_2_y = 0.0;
			tmp_x_2_z = 0.0;

			tmp_y_2_x = 0.0;
			tmp_y_2_y = 0.0;
			tmp_y_2_z = 0.0;

			tmp_z_2_x = 0.0;
			tmp_z_2_y = 0.0;
			tmp_z_2_z = 0.0;
			for(ii3=0; ii3<velo_node_num; ii3++){
				tmp_x_1_x += Kx_tensor[ii1][ii3][ii2]
				           * elem_flow_vect[ii3+base_node_num*1];
				tmp_y_1_y += Ky_tensor[ii1][ii3][ii2]
				           * elem_flow_vect[ii3+base_node_num*2];

				tmp_x_2_x += Kx_tensor[ii1][ii2][ii3]
				           * elem_flow_vect[ii3+base_node_num*1];
				tmp_x_2_y += Kx_tensor[ii1][ii2][ii3]
				           * elem_flow_vect[ii3+base_node_num*2];

				tmp_y_2_x += Ky_tensor[ii1][ii2][ii3]
				           * elem_flow_vect[ii3+base_node_num*1];
				tmp_y_2_y += Ky_tensor[ii1][ii2][ii3]
				           * elem_flow_vect[ii3+base_node_num*2];

				if(dim == 3){
					tmp_z_1_z += Kz_tensor[ii1][ii3][ii2]
				   	           * elem_flow_vect[ii3+base_node_num*3];
					tmp_x_2_z += Kx_tensor[ii1][ii2][ii3]
				   	           * elem_flow_vect[ii3+base_node_num*3];
					tmp_y_2_z += Ky_tensor[ii1][ii2][ii3]
				   	           * elem_flow_vect[ii3+base_node_num*3];
					tmp_z_2_x += Kz_tensor[ii1][ii2][ii3]
				   	           * elem_flow_vect[ii3+base_node_num*1];
					tmp_z_2_y += Kz_tensor[ii1][ii2][ii3]
				   	           * elem_flow_vect[ii3+base_node_num*2];
					tmp_z_2_z += Kz_tensor[ii1][ii2][ii3]
				   	           * elem_flow_vect[ii3+base_node_num*3];
				}
			}

			/* Gx_matrix */
			elem_matrix->matrix[ii1+base_node_num*1][ii2+base_node_num*1]
			= tmp_x_2_x + tmp_x_1_x + tmp_y_1_y
			+ 2.0 * Sxx_matrix[ii1][ii2] + Syy_matrix[ii1][ii2];

			/* Lxy_matrix */
			elem_matrix->matrix[ii1+base_node_num*1][ii2+base_node_num*2]
			= tmp_y_2_x + Sxy_matrix[ii2][ii1];

			/* Lyx_matrix */
			elem_matrix->matrix[ii1+base_node_num*2][ii2+base_node_num*1]
			= tmp_x_2_y + Sxy_matrix[ii1][ii2];

			/* Gy_matrix */
			elem_matrix->matrix[ii1+base_node_num*2][ii2+base_node_num*2]
			= tmp_y_2_y + tmp_y_1_y + tmp_x_1_x
			+ 2.0 * Syy_matrix[ii1][ii2] + Sxx_matrix[ii1][ii2];

			/* I_vector */
			(*local_vect)[ii1+ base_node_num*0]
			+= ( Hx_matrix[ii2][ii1] * elem_flow_vect[ii2+base_node_num*1]
			   + Hy_matrix[ii2][ii1] * elem_flow_vect[ii2+base_node_num*2] );

			/* Fx_vector */
			(*local_vect)[ii1+ base_node_num*1]
			-= (tmp_x_1_x * elem_flow_vect[ii2+base_node_num*1]
			  + tmp_y_1_y * elem_flow_vect[ii2+base_node_num*1]
			  + (2.0 * Sxx_matrix[ii1][ii2] + Syy_matrix[ii1][ii2])
			   * elem_flow_vect[ii2+base_node_num*1]
			  + Sxy_matrix[ii2][ii1] * elem_flow_vect[ii2+base_node_num*2]);

			/* Fy_vector */
			(*local_vect)[ii1+ base_node_num*2]
			-= (tmp_x_1_x * elem_flow_vect[ii2+base_node_num*2]
			  + tmp_y_1_y * elem_flow_vect[ii2+base_node_num*2]
			  + Sxy_matrix[ii1][ii2] * elem_flow_vect[ii2+base_node_num*1]
			  + (2.0 * Syy_matrix[ii1][ii2] + Sxx_matrix[ii1][ii2])
			   * elem_flow_vect[ii2 + base_node_num*2]);
			
			if(dim == 3){
				/* Gx_matrix */
				elem_matrix->matrix[ii1+base_node_num*1][ii2+base_node_num*1]
				+= tmp_z_1_z + Szz_matrix[ii1][ii2];

				/* Gy_matrix */
				elem_matrix->matrix[ii1+base_node_num*2][ii2+base_node_num*2]
				+= tmp_z_1_z + Szz_matrix[ii1][ii2];

				/* Gz_matrix */
				elem_matrix->matrix[ii1+base_node_num*3][ii2+base_node_num*3]
				= tmp_x_1_x + tmp_y_1_y + tmp_z_1_z + tmp_z_2_z
				+       Sxx_matrix[ii1][ii2]
				+       Syy_matrix[ii1][ii2]
				+ 2.0 * Szz_matrix[ii1][ii2];

				/* Lxz_matrix */
				elem_matrix->matrix[ii1+base_node_num*1][ii2+base_node_num*3]
				= tmp_z_2_x + Szx_matrix[ii1][ii2];

				/* Lzx_matrix */
				elem_matrix->matrix[ii1+base_node_num*3][ii2+base_node_num*1]
				= tmp_x_2_z + Szx_matrix[ii2][ii1];

				/* Lyz_matrix */
				elem_matrix->matrix[ii1+base_node_num*2][ii2+base_node_num*3]
				= tmp_z_2_y + Syz_matrix[ii2][ii1];

				/* Lzy_matrix */
				elem_matrix->matrix[ii1+base_node_num*3][ii2+base_node_num*2]
				= tmp_y_2_z + Syz_matrix[ii1][ii2];

				/* I_vector */
				(*local_vect)[ii1+ base_node_num*0]
				+= Hz_matrix[ii2][ii1] * elem_flow_vect[ii2+base_node_num*3];

				/* Fx_vector */
				(*local_vect)[ii1+ base_node_num*1]
				-= (tmp_z_1_z * elem_flow_vect[ii2+base_node_num*1]
				  + Szz_matrix[ii1][ii2] * elem_flow_vect[ii2+base_node_num*1]
				  + Szx_matrix[ii1][ii2] * elem_flow_vect[ii2+base_node_num*3]);

				/* Fy_vector */
				(*local_vect)[ii1+ base_node_num*2]
				-= (tmp_z_1_z * elem_flow_vect[ii2+base_node_num*2]
				  + Szz_matrix[ii1][ii2] * elem_flow_vect[ii2+base_node_num*2]
				  + Syz_matrix[ii2][ii1] * elem_flow_vect[ii2+base_node_num*3]);

				/* Fz_vector */
				(*local_vect)[ii1+ base_node_num*3]
				-= (tmp_x_1_x * elem_flow_vect[ii2+base_node_num*3]
			 	  + tmp_y_1_y * elem_flow_vect[ii2+base_node_num*3]
			 	  + tmp_z_1_z * elem_flow_vect[ii2+base_node_num*3]
			 	  + Szx_matrix[ii2][ii1] * elem_flow_vect[ii2+base_node_num*1]
			 	  + Syz_matrix[ii1][ii2] * elem_flow_vect[ii2+base_node_num*2]
			 	  + (        Sxx_matrix[ii1][ii2]
			 	     +       Syy_matrix[ii1][ii2]
				     + 2.0 * Szz_matrix[ii1][ii2])
			 	   * elem_flow_vect[ii2+base_node_num*3]);
			}

		}

		for(ii2=0; ii2<press_node_num; ii2++){
			elem_matrix->matrix[ii2+base_node_num*0][ii1+base_node_num*1]
			= - Hx_matrix[ii1][ii2];

			elem_matrix->matrix[ii2+base_node_num*0][ii1+base_node_num*2]
			= - Hy_matrix[ii1][ii2];

			elem_matrix->matrix[ii1+base_node_num*1][ii2+base_node_num*0]
			= - Hx_matrix[ii1][ii2];

			elem_matrix->matrix[ii1+base_node_num*2][ii2+base_node_num*0]
			= - Hy_matrix[ii1][ii2];

			(*local_vect)[ii1+ base_node_num*1]
			+= (Hx_matrix[ii1][ii2] * elem_flow_vect[ii2+base_node_num*0]);

			(*local_vect)[ii1+ base_node_num*2]
			+= (Hy_matrix[ii1][ii2] * elem_flow_vect[ii2+base_node_num*0]);

			if(dim == 3){
				elem_matrix->matrix[ii2+base_node_num*0][ii1+base_node_num*3]
				= - Hz_matrix[ii1][ii2];
				
				elem_matrix->matrix[ii1+base_node_num*3][ii2+base_node_num*0]
				= - Hz_matrix[ii1][ii2];

				(*local_vect)[ii1+ base_node_num*3]
				+= (Hz_matrix[ii1][ii2] * elem_flow_vect[ii2+base_node_num*0]);
			}
		}

	}

	/* 一次のパラメータに関する中間節点の変化量を 0 にする */
	for(ii1=press_node_num; ii1<base_node_num; ii1++){
		elem_matrix->matrix[ii1+base_node_num*0][ii1+base_node_num*0] = 1.0;
		(*local_vect)[ii1+ base_node_num*0] = 0.0;
	}
	for(ii1=velo_node_num; ii1<base_node_num; ii1++){
		elem_matrix->matrix[ii1+base_node_num*1][ii1+base_node_num*1] = 1.0;
		(*local_vect)[ii1+base_node_num*1] = 0.0;
		elem_matrix->matrix[ii1+base_node_num*2][ii1+base_node_num*2] = 1.0;
		(*local_vect)[ii1+base_node_num*2] = 0.0;
		if(dim == 3){
			elem_matrix->matrix[ii1+base_node_num*3][ii1+base_node_num*3] = 1.0;
			(*local_vect)[ii1+base_node_num*3] = 0.0;
		}
	}

	free( (void *)elem_flow_vect);

	return(NORMAL_RC);
}


static RC fem_velo2vect(int dim, FEM_NODE_ARRAY node, FEM_VELO_ARRAY flow,
                        int valid_node_num, double *flow_vect)
{
	int index;
	int ii1;

	for(ii1=0; ii1<flow.size; ii1++){
		if(flow.array[ii1].node < 0) continue;
		index = node_renum_index(node, flow.array[ii1].node);
		RC_NEG_CHK(index);

		flow_vect[valid_node_num*0 + index] = flow.array[ii1].s;
		flow_vect[valid_node_num*1 + index] = flow.array[ii1].v.x;
		flow_vect[valid_node_num*2 + index] = flow.array[ii1].v.y;
		if(dim == 3){
			flow_vect[valid_node_num*3 + index] = flow.array[ii1].v.z;
		}
	}

	return(NORMAL_RC);
}


static RC fem_vect2velo(int dim, FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY elem,
                        int velo_order, int press_order,
                        int valid_node_num,
                        const double *flow_vect, FEM_VELO_ARRAY *flow)
{
	int ii1, ii2, ii3;
	int index;
	double **local_coords;
	double N[FEM_MAX_NODE];
	double value[3];
	FEM_ELEM_TYPE velo_type;
	FEM_ELEM_TYPE press_type;
	int velo_node_num;
	int press_node_num;

	RC_TRY( allocate2D(FEM_MAX_NODE, 4, &local_coords) );

	RC_TRY( free_fem_velo_array(flow) );
	RC_TRY( allocate_fem_velo_array(node.size, flow) );

	for(ii1=0; ii1<node.size; ii1++){
		if(node.array[ii1].label < 0){
			flow->array[ii1].node = -1;
			continue;
		}
		RC_NEG_CHK( index = node.array[ii1].renum );

		flow->array[ii1].node = node.array[ii1].label;
		flow->array[ii1].s   = flow_vect[valid_node_num*0 + index];
		flow->array[ii1].v.x = flow_vect[valid_node_num*1 + index];
		flow->array[ii1].v.y = flow_vect[valid_node_num*2 + index];
		if(dim == 3){
			flow->array[ii1].v.z = flow_vect[valid_node_num*3 + index];
		}
	}
	flow->size = node.size;
	flow->source = FEM_K_SOL;

	if(press_order == 1){
		for(ii1=0; ii1<elem.size; ii1++){
			if(elem.array[ii1].label < 0) continue;

			RC_TRY( parameter_check(elem.array[ii1], 2, press_order,
			                        &press_type, &press_node_num) );
			RC_TRY( set_local_node_points(elem.array[ii1].type, local_coords) );

			for(ii2=press_node_num; ii2<elem.array[ii1].node_num; ii2++){ 
				RC_TRY( set_N_vector(press_type, local_coords[ii2], N) );

				value[0] = 0.0;
				for(ii3=0; ii3<press_node_num; ii3++){
					index = search_fem_velo_node(*flow,
					                             elem.array[ii1].node[ii3]);
					RC_NEG_CHK(index);
					value[0] += N[ii3] * flow->array[index].s;
				}

				index = search_fem_velo_node(*flow, elem.array[ii1].node[ii2]);
				RC_NEG_CHK(index);
				flow->array[index].s = value[0];
			}
		}
	}

	if(velo_order == 1){
		for(ii1=0; ii1<elem.size; ii1++){
			if(elem.array[ii1].label < 0) continue;

			RC_TRY( parameter_check(elem.array[ii1], 2, velo_order,
			                        &velo_type, &velo_node_num) );
			RC_TRY( set_local_node_points(elem.array[ii1].type, local_coords) );

			for(ii2=velo_node_num; ii2<elem.array[ii1].node_num; ii2++){ 
				RC_TRY( set_N_vector(velo_type, local_coords[ii2], N) );

				value[0] = 0.0;
				value[1] = 0.0;
				value[2] = 0.0;
				for(ii3=0; ii3<velo_node_num; ii3++){
					index = search_fem_velo_node(*flow,
					                             elem.array[ii1].node[ii3]);
					RC_NEG_CHK(index);
					value[0] += N[ii3] * flow->array[index].v.x;
					value[1] += N[ii3] * flow->array[index].v.y;
					value[2] += N[ii3] * flow->array[index].v.z;
				}

				index = search_fem_velo_node(*flow, elem.array[ii1].node[ii2]);
				RC_NEG_CHK(index);
				flow->array[index].v.x = value[0];
				flow->array[index].v.y = value[1];
				if(dim == 3){
					flow->array[index].v.z = value[2];
				}
			}
		}
	}

	RC_TRY( free2D(FEM_MAX_NODE, 4, local_coords) );

	return(NORMAL_RC);
}


static RC parameter_check(FEM_ELEMENT elem, int model_order, int para_order,
                          FEM_ELEM_TYPE *para_type, int *para_node_num)
{

	if( (model_order == 2) && (para_order == 1) ){
		RC_TRY( reduced_element_type(elem.type, para_type, para_node_num) );
		return(NORMAL_RC);
	}

	*para_type = elem.type;
	*para_node_num = elem.node_num;

	return(NORMAL_RC);
}


RC steady_adjflow_analysis(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY elem,
                           FEM_BC_ARRAY bc_velo, FEM_BC_ARRAY bc_press,
                           int velo_order, int press_order,
                           double density, double viscosity,
                           FEM_VELO_ARRAY flow, FEM_VELO_ARRAY *adjoint_flow)
{
	int ii1;
	int dim;
	int model_order;
	int valid_node_num;

	NONZERO_MATRIX matrix;
	double *global_vect;
	double *flow_vect;
	double *adjoint_flow_vect;

	FEM_ELEM_MATRIX elem_matrix;
	double *local_vect;

	FEM_WORK_SET ws;

	RC_TRY( allocate_fem_work_set(&ws) );

	dim = fem_analysis_dim(elem);
	if( (dim != 2) && (dim != 3) ){
		fprintf(stderr, "fem_analysis_dim <\n");
		return(IMPLEMENT_ERROR_RC);
	}

	model_order = fem_analysis_order(elem);
	if(model_order != 2){
		fprintf(stderr, "fem_model_order <\n");
		return(IMPLEMENT_ERROR_RC);
	}
	
	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );
	
	RC_NEG_ZERO_CHK(velo_order);
	RC_NEG_ZERO_CHK(press_order);
	if( (velo_order > model_order) || (press_order > model_order) ){
		return(OVERFLOW_ERROR_RC); 
	}

	valid_node_num = count_valid_node(node);

	/* FEM_VELO_ARRAY の確保 */
	RC_TRY( allocate_fem_velo_array(node.size, adjoint_flow) );
	adjoint_flow->size = node.size;
	for(ii1=0; ii1<node.size; ii1++){
		if(node.array[ii1].label < 0){
			adjoint_flow->array[ii1].node = -1;
			continue;
		}
		adjoint_flow->array[ii1].node = node.array[ii1].label;
	}

	/* 計算専用ベクトルの確保 */
	RC_TRY( allocate1D(valid_node_num*(dim+1), &flow_vect) );

	/* FEM_VELO_ARRAY -> 計算専用ベクトル の変換 */
	RC_TRY( fem_velo2vect(dim, node, flow, valid_node_num, flow_vect) );

	/* adjoint_flow用ベクトルの確保 */
	RC_TRY( allocate1D(valid_node_num*(dim+1), &adjoint_flow_vect) );

	RC_TRY( allocate_nonzero_matrix1a(valid_node_num*(dim+1), &matrix) );
	RC_TRY( allocate1D(valid_node_num*(dim+1), &global_vect) );

	for(ii1=0; ii1<elem.size; ii1++){
		if(elem.array[ii1].label < 0) continue;

		RC_TRY( make_adjoint_elem_matrix( dim, model_order, elem.array[ii1],
		  	                              node, valid_node_num,
		                                  velo_order, press_order,
		                                  density, viscosity,
		                                  flow_vect, &elem_matrix,
		                                  &local_vect, ws) );

		RC_TRY( fem_impose_nonzero1a(elem_matrix, local_vect, &matrix,
		                             global_vect) );

		RC_TRY( free_fem_elem_matrix(&elem_matrix) );
		free( (void *)local_vect );
	}

	RC_TRY( modify_nonzero_matrix1a(bc_velo, bc_press, node, valid_node_num,
	                                &matrix, global_vect) );

	RC_TRY( nonzero1a_scaling(&matrix, global_vect) );
	RC_TRY( cg_nonzero1a_ver1(matrix, global_vect, adjoint_flow_vect) );

	RC_TRY( free_nonzero_matrix1a(&matrix) );
	free( (void *)global_vect );
	free( (void *)flow_vect );

	RC_TRY( fem_vect2velo(dim, node, elem, velo_order, press_order,
	                      valid_node_num, adjoint_flow_vect, adjoint_flow) );

	free( (void *)adjoint_flow_vect );
	RC_TRY( free_fem_work_set(&ws) );

	return(NORMAL_RC);
}


/* 前半部は make_elem_matrix() と共通部分が多いため，モジュール化可能 */
static RC make_adjoint_elem_matrix(int dim, int model_order, FEM_ELEMENT elem,
                                   FEM_NODE_ARRAY node, int valid_node_num,
                                   int velo_order, int press_order,
                                   double density, double viscosity,
                                   const double *flow_vect,
                                   FEM_ELEM_MATRIX *elem_matrix,
                                   double **local_vect, FEM_WORK_SET ws)
{
	static int allocate_flag = 0;
	static double ***Kx_tensor = NULL;
	static double ***Ky_tensor = NULL;
	static double ***Kz_tensor = NULL;
	static double **Sxx_matrix = NULL;
	static double **Syy_matrix = NULL;
	static double **Szz_matrix = NULL;
	static double **Sxy_matrix = NULL;
	static double **Syz_matrix = NULL;
	static double **Szx_matrix = NULL;
	static double **Hx_matrix = NULL;
	static double **Hy_matrix = NULL;
	static double **Hz_matrix = NULL;

	double **node_location = ws.mat_N_3[0];
	double **velo_dN = ws.mat_3_N[0];
	double **velo_dNxyz = ws.mat_3_N[1];
	double **Jacobi = ws.mat_3_3[0];
	double **inverse_J = ws.mat_3_3[1];
	double **g_points = ws.mat_G_4[0];
	double *velo_N = ws.vec_N[0];
	double *press_N = ws.vec_N[1];
	double *weight = ws.vec_G[0];
	FEM_ELEM_TYPE velo_type;
	FEM_ELEM_TYPE press_type;
	FEM_ELEM_TYPE base_type;
	double det_J;
	int velo_node_num;
	int press_node_num;
	int base_node_num;
	int g_point_num;
	int index;
	int ii1, ii2, ii3, ii4;
	double *elem_flow_vect;

	double tmp_x_1_x;
	double tmp_y_1_y;
	double tmp_z_1_z;

	double tmp_x_2_x;
	double tmp_x_2_y;
	double tmp_x_2_z;

	double tmp_y_2_x;
	double tmp_y_2_y;
	double tmp_y_2_z;

	double tmp_z_2_x;
	double tmp_z_2_y;
	double tmp_z_2_z;

	if(allocate_flag == 0){
		RC_TRY( allocate3D(FEM_MAX_NODE, FEM_MAX_NODE, FEM_MAX_NODE,
		                   &Kx_tensor) );
		RC_TRY( allocate3D(FEM_MAX_NODE, FEM_MAX_NODE, FEM_MAX_NODE,
		                   &Ky_tensor) );

		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &Sxx_matrix) );
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &Sxy_matrix) );
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &Syy_matrix) );
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &Hx_matrix) );
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &Hy_matrix) );

		if(dim == 3){
			RC_TRY( allocate3D(FEM_MAX_NODE, FEM_MAX_NODE, FEM_MAX_NODE,
			                   &Kz_tensor) );
			RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &Szz_matrix) );
			RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &Syz_matrix) );
			RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &Szx_matrix) );
			RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &Hz_matrix) );
		}

		allocate_flag = 1;
	}

	/* 初期化 */
	for(ii1=0; ii1<FEM_MAX_NODE; ii1++){
		 init_matrix(FEM_MAX_NODE, FEM_MAX_NODE, Kx_tensor[ii1]);
		 init_matrix(FEM_MAX_NODE, FEM_MAX_NODE, Ky_tensor[ii1]);
	}
	init_matrix(FEM_MAX_NODE, FEM_MAX_NODE, Sxx_matrix);
	init_matrix(FEM_MAX_NODE, FEM_MAX_NODE, Syy_matrix);
	init_matrix(FEM_MAX_NODE, FEM_MAX_NODE, Sxy_matrix);
	init_matrix(FEM_MAX_NODE, FEM_MAX_NODE, Hx_matrix);
	init_matrix(FEM_MAX_NODE, FEM_MAX_NODE, Hy_matrix);
	if(dim == 3){
		for(ii1=0; ii1<FEM_MAX_NODE; ii1++){
			init_matrix(FEM_MAX_NODE, FEM_MAX_NODE, Kz_tensor[ii1]);
		}
		init_matrix(FEM_MAX_NODE, FEM_MAX_NODE, Szz_matrix);
		init_matrix(FEM_MAX_NODE, FEM_MAX_NODE, Syz_matrix);
		init_matrix(FEM_MAX_NODE, FEM_MAX_NODE, Szx_matrix);
		init_matrix(FEM_MAX_NODE, FEM_MAX_NODE, Hz_matrix);
	}
		

	RC_TRY( parameter_check(elem, model_order, velo_order, &velo_type,
	                                                       &velo_node_num) );
	RC_TRY( parameter_check(elem, model_order, press_order, &press_type,
	                                                        &press_node_num) );

	base_type = elem.type;
	base_node_num = elem.node_num;

	RC_TRY( Gauss_const_max(base_type, g_points, weight, &g_point_num) );
	RC_TRY( node_location_matrix(&elem, node, node_location) );

	for(ii1=0; ii1<g_point_num; ii1++){
			
		/* 流速 */
		RC_TRY( set_N_vector(velo_type, g_points[ii1], velo_N) );
		RC_TRY( set_dN_matrix(velo_type, g_points[ii1], velo_dN) );
		mul_matrix(dim, velo_node_num, dim, velo_dN, node_location, Jacobi);
		det_J = determinant(dim, Jacobi);
		if(det_J < 0){
			fprintf(stderr, "Warning! Jacobian is negative!!!\n");
			return(NEGATIVE_ERROR_RC);
		}

		RC_TRY( inverse_matrix(dim, Jacobi, inverse_J) );
		mul_matrix(dim, dim, velo_node_num, inverse_J, velo_dN, velo_dNxyz);

		/* 圧力 */
		RC_TRY( set_N_vector(press_type, g_points[ii1], press_N) );


		for(ii2=0; ii2<velo_node_num; ii2++){
			for(ii3=0; ii3<velo_node_num; ii3++){

				/* S_matrix[velo_node_num][velo_node_num] */
				Sxx_matrix[ii2][ii3] += velo_dNxyz[0][ii2] * velo_dNxyz[0][ii3]
				                      * det_J * weight[ii1] * viscosity; 
				Syy_matrix[ii2][ii3] += velo_dNxyz[1][ii2] * velo_dNxyz[1][ii3]
				                      * det_J * weight[ii1] * viscosity;
				Sxy_matrix[ii2][ii3] += velo_dNxyz[0][ii2] * velo_dNxyz[1][ii3]
				                      * det_J * weight[ii1] * viscosity;
				if(dim == 3){	
					Szz_matrix[ii2][ii3] += velo_dNxyz[2][ii2]
					                      * velo_dNxyz[2][ii3]
										  * det_J * weight[ii1] * viscosity;
					Syz_matrix[ii2][ii3] += velo_dNxyz[1][ii2]
					                      * velo_dNxyz[2][ii3]
										  * det_J * weight[ii1] * viscosity;
					Szx_matrix[ii2][ii3] += velo_dNxyz[2][ii2]
					                      * velo_dNxyz[0][ii3]
					                      * det_J * weight[ii1] * viscosity;
				}

				for(ii4=0; ii4<velo_node_num; ii4++){
					/* K_tensor[velo_node_num][velo_node_num][velo_node_num] */
					Kx_tensor[ii2][ii3][ii4] += velo_N[ii2] * velo_N[ii3]
					                          * velo_dNxyz[0][ii4] * det_J
					                          * weight[ii1] * density;
					Ky_tensor[ii2][ii3][ii4] += velo_N[ii2] * velo_N[ii3]
					                          * velo_dNxyz[1][ii4] * det_J
					                          * weight[ii1] * density;
					if(dim == 3){
						Kz_tensor[ii2][ii3][ii4] += velo_N[ii2] * velo_N[ii3]
						                          * velo_dNxyz[2][ii4] * det_J
						                          * weight[ii1] * density;
					}
				}

			}

			for(ii3=0; ii3<press_node_num; ii3++){
				/* H_matrix[velo_node_num][press_node_num] */
				Hx_matrix[ii2][ii3] += velo_dNxyz[0][ii2] * press_N[ii3]
				                     * det_J * weight[ii1];
				Hy_matrix[ii2][ii3] += velo_dNxyz[1][ii2] * press_N[ii3]
				                     * det_J * weight[ii1];
				if(dim == 3){
					Hz_matrix[ii2][ii3] += velo_dNxyz[2][ii2] * press_N[ii3]
					                     * det_J * weight[ii1];
				}
			}
		}

	}

	/* elem_matrix の確保 */
	elem_matrix->dim = dim;
	elem_matrix->element = elem.label;
	elem_matrix->size = base_node_num * (elem_matrix->dim + 1);
	elem_matrix->index = malloc( elem_matrix->size * sizeof(int) );
	RC_NULL_CHK(elem_matrix->index);
    RC_TRY( allocate2D( elem_matrix->size, elem_matrix->size,
	                    &(elem_matrix->matrix) ) );

	RC_TRY( allocate1D(elem_matrix->size, &elem_flow_vect) );

	/* elem_matrix.index[] の作成 */
	for(ii1=0; ii1<elem.node_num; ii1++){
		index = node_renum_index(node, elem.node[ii1]);
		RC_NEG_CHK(index);
		for(ii2=0; ii2<(elem_matrix->dim)+1; ii2++){
			elem_matrix->index[ii2*base_node_num + ii1]
        	= ii2 * valid_node_num + index;
			elem_flow_vect[ii2*base_node_num + ii1]
			= flow_vect[ii2*valid_node_num + index];
		}
	}

	/* 要素マトリックス用の右辺ベクトル作成 */
	RC_TRY( allocate1D(elem_matrix->size, local_vect) );

	for(ii1=0; ii1<velo_node_num; ii1++){
		for(ii2=0; ii2<velo_node_num; ii2++){

			/* tmp_x_1_x += tensor_x[][loop][] * velo_x[loop] */
			tmp_x_1_x = 0.0;
			tmp_y_1_y = 0.0;
			tmp_z_1_z = 0.0;

			tmp_x_2_x = 0.0;
			tmp_x_2_y = 0.0;
			tmp_x_2_z = 0.0;

			tmp_y_2_x = 0.0;
			tmp_y_2_y = 0.0;
			tmp_y_2_z = 0.0;

			tmp_z_2_x = 0.0;
			tmp_z_2_y = 0.0;
			tmp_z_2_z = 0.0;

			for(ii3=0; ii3<velo_node_num; ii3++){
				tmp_x_1_x += Kx_tensor[ii1][ii3][ii2]
				           * elem_flow_vect[ii3 + base_node_num*1];
				tmp_x_2_x += Kx_tensor[ii1][ii2][ii3]
				           * elem_flow_vect[ii3 + base_node_num*1];
				tmp_x_2_y += Kx_tensor[ii1][ii2][ii3]
				           * elem_flow_vect[ii3 + base_node_num*2];
				tmp_y_1_y += Ky_tensor[ii1][ii3][ii2]
				           * elem_flow_vect[ii3 + base_node_num*2];
				tmp_y_2_y += Ky_tensor[ii1][ii2][ii3]
				           * elem_flow_vect[ii3 + base_node_num*2];
				tmp_y_2_x += Ky_tensor[ii1][ii2][ii3]
				           * elem_flow_vect[ii3 + base_node_num*1];
				if(dim == 3){
					tmp_z_1_z += Kz_tensor[ii1][ii3][ii2]
					           * elem_flow_vect[ii3+base_node_num*3];
					tmp_x_2_z += Kx_tensor[ii1][ii2][ii3]
					           * elem_flow_vect[ii3+base_node_num*3];
					tmp_y_2_z += Ky_tensor[ii1][ii2][ii3]
					           * elem_flow_vect[ii3+base_node_num*3];
					tmp_z_2_x += Kz_tensor[ii1][ii2][ii3]
					           * elem_flow_vect[ii3+base_node_num*1];
					tmp_z_2_y += Kz_tensor[ii1][ii2][ii3]
					           * elem_flow_vect[ii3+base_node_num*2];
					tmp_z_2_z += Kz_tensor[ii1][ii2][ii3]
					           * elem_flow_vect[ii3+base_node_num*3];
				}
			}

			/* Gx_matrix */
			elem_matrix->matrix[ii1 + base_node_num*1][ii2 + base_node_num*1]
			= tmp_x_2_x + tmp_x_1_x + tmp_y_1_y
			+ 2.0 * Sxx_matrix[ii1][ii2] + Syy_matrix[ii1][ii2];

			/* Lxy_matrix */
			elem_matrix->matrix[ii1 + base_node_num*1][ii2 + base_node_num*2]
			= tmp_x_2_y + Sxy_matrix[ii1][ii2];

			/* Lyx_matrix */
			elem_matrix->matrix[ii1 + base_node_num*2][ii2 + base_node_num*1]
			= tmp_y_2_x + Sxy_matrix[ii2][ii1];

			/* Gy_matrix */
			elem_matrix->matrix[ii1 + base_node_num*2][ii2 + base_node_num*2]
			= tmp_y_2_y + tmp_x_1_x + tmp_y_1_y
			+ 2.0 * Syy_matrix[ii1][ii2] + Sxx_matrix[ii1][ii2];

			/* Fx_vector */
			(*local_vect)[ii1+ base_node_num*1]
			+= 4.0 * Sxx_matrix[ii1][ii2] * elem_flow_vect[ii2+base_node_num*1]
			 + 2.0 * Syy_matrix[ii1][ii2] * elem_flow_vect[ii2+base_node_num*1]
			 + 2.0 * Sxy_matrix[ii1][ii2] * elem_flow_vect[ii2+base_node_num*2];

			/* Fy_vector */
			(*local_vect)[ii1+ base_node_num*2]
			+= 4.0 * Syy_matrix[ii1][ii2] * elem_flow_vect[ii2+base_node_num*2]
			 + 2.0 * Sxx_matrix[ii1][ii2] * elem_flow_vect[ii2+base_node_num*2]
			 + 2.0 * Sxy_matrix[ii2][ii1] * elem_flow_vect[ii2+base_node_num*1];

			if(dim == 3){
				/* Gx_matrix */
				elem_matrix->matrix[ii1+base_node_num*1][ii2+base_node_num*1]
				+= tmp_z_1_z + Szz_matrix[ii1][ii2];

				/* Gy_matrix */
				elem_matrix->matrix[ii1+base_node_num*2][ii2+base_node_num*2]
				+= tmp_z_1_z + Szz_matrix[ii1][ii2];

				/* Gz_matrix */
				elem_matrix->matrix[ii1+base_node_num*3][ii2+base_node_num*3]
				= tmp_z_2_z + tmp_x_1_x + tmp_y_1_y + tmp_z_1_z
				+       Sxx_matrix[ii1][ii2]
				+       Syy_matrix[ii1][ii2]
				+ 2.0 * Szz_matrix[ii1][ii2];

				/* Lxz_matrix */
				elem_matrix->matrix[ii1+base_node_num*1][ii2+base_node_num*3]
				= tmp_x_2_z + Szx_matrix[ii2][ii1];

				/* Lzx_matrix */
				elem_matrix->matrix[ii1+base_node_num*3][ii2+base_node_num*1]
				= tmp_z_2_x + Szx_matrix[ii1][ii2];

				/* Lyz_matrix */
				elem_matrix->matrix[ii1+base_node_num*2][ii2+base_node_num*3]
				= tmp_y_2_z + Syz_matrix[ii1][ii2];

				/* Lzy_matrix */
				elem_matrix->matrix[ii1+base_node_num*3][ii2+base_node_num*2]
				= tmp_z_2_y + Syz_matrix[ii2][ii1];

				/* Fx_vector */
				(*local_vect)[ii1+ base_node_num*1]
				+= 2.0 * Szz_matrix[ii1][ii2]
				       * elem_flow_vect[ii2+base_node_num*1]
				 + 2.0 * Szx_matrix[ii2][ii1]
				       * elem_flow_vect[ii2+base_node_num*3];

				/* Fy_vector */
				(*local_vect)[ii1+ base_node_num*2]
				+= 2.0 * Szz_matrix[ii1][ii2]
				       * elem_flow_vect[ii2+base_node_num*2]
				 + 2.0 * Syz_matrix[ii1][ii2]
				       * elem_flow_vect[ii2+base_node_num*3];

				/* Fz_vector */
				(*local_vect)[ii1+ base_node_num*3]
				+= 4.0 * Szz_matrix[ii1][ii2]
				       * elem_flow_vect[ii2+base_node_num*3]
				 + 2.0 * Sxx_matrix[ii1][ii2]
				       * elem_flow_vect[ii2+base_node_num*3]
				 + 2.0 * Syy_matrix[ii1][ii2]
				       * elem_flow_vect[ii2+base_node_num*3]
				 + 2.0 * Szx_matrix[ii1][ii2]
				       * elem_flow_vect[ii2+base_node_num*1]
				 + 2.0 * Syz_matrix[ii2][ii1]
				       * elem_flow_vect[ii2+base_node_num*2];
			}
		}

		for(ii2=0; ii2<press_node_num; ii2++){
			elem_matrix->matrix[ii2 + base_node_num*0][ii1 + base_node_num*1]
			= - Hx_matrix[ii1][ii2];

			elem_matrix->matrix[ii2 + base_node_num*0][ii1 + base_node_num*2]
			= - Hy_matrix[ii1][ii2];

			elem_matrix->matrix[ii1 + base_node_num*1][ii2 + base_node_num*0]
			= - Hx_matrix[ii1][ii2];

			elem_matrix->matrix[ii1 + base_node_num*2][ii2 + base_node_num*0]
			= - Hy_matrix[ii1][ii2];
		}

		if(dim == 3){
			elem_matrix->matrix[ii2+base_node_num*0][ii1+base_node_num*3]
			= - Hz_matrix[ii1][ii2];

			elem_matrix->matrix[ii1+base_node_num*3][ii2+base_node_num*0]
			= - Hz_matrix[ii1][ii2];
		}

	}

	/* 一次のパラメータに関する中間節点の変化量を 0 にする */
	for(ii1=press_node_num; ii1<base_node_num; ii1++){
		elem_matrix->matrix[ii1 + base_node_num*0][ii1 + base_node_num*0] = 1.0;
		(*local_vect)[ii1+ base_node_num*0] = 0.0;
	}
	for(ii1=velo_node_num; ii1<base_node_num; ii1++){
		elem_matrix->matrix[ii1 + base_node_num*1][ii1 + base_node_num*1] = 1.0;
		(*local_vect)[ii1+ base_node_num*1] = 0.0;
		elem_matrix->matrix[ii1 + base_node_num*2][ii1 + base_node_num*2] = 1.0;
		(*local_vect)[ii1+ base_node_num*2] = 0.0;
		if(dim == 3){
			elem_matrix->matrix[ii1+base_node_num*3][ii1+base_node_num*3] = 1.0;
			(*local_vect)[ii1+base_node_num*3] = 0.0;
		}
	}

	free( (void *)elem_flow_vect );

	return(NORMAL_RC);
}






/*
 * 非定常流れ場の解析
 * 
 * 連立方程式の解法 SMAC
 */

/* 最適化解析における流れ場解析の実行(余計なファイルの出力は無し) */
RC flow_solve_for_optimization(double delta_t, int iteration, double Re, 
                               FEM_NODE_ARRAY node, 
                               FEM_ELEMENT_ARRAY element, 
                               FEM_BC_ARRAY bc_velo, FEM_BC_ARRAY bc_press,
                               FEM_VELO_ARRAY *vp)
{
	int ii1, ii2;
	int dim ;
	FEM_VELO_ARRAY tmp_vp1;
	FEM_VELO_ARRAY tmp_vp2;

	RC_TRY( allocate_fem_velo_array(node.size, &tmp_vp1) );
	RC_TRY( allocate_fem_velo_array(node.size, &tmp_vp2) );
	RC_TRY( allocate_fem_velo_array(node.size, vp) );

	for(ii1=0;ii1<node.size;ii1++){
		tmp_vp1.array[ii1].node = node.array[ii1].label;
		tmp_vp2.array[ii1].node = node.array[ii1].label;
		vp->array[ii1].node = node.array[ii1].label;
	}

	tmp_vp1.size = node.size;
	tmp_vp2.size = node.size;
	vp->size = node.size;

	dim = fem_analysis_dim(element); 

	for(ii1=0; ii1<iteration; ii1++){
		
		/* 中間流速の計算 */
		RC_TRY( solve_middle_velocity(delta_t, Re, node, element, bc_velo,
		                              tmp_vp1, &tmp_vp2) );
		/* 圧力の計算 */
		RC_TRY( solve_pressure(dim, delta_t, node, element, &tmp_vp2,
		                       bc_press) ); 
		/* 流速の計算 */
		RC_TRY( solve_velocity(delta_t, node, element, bc_velo, tmp_vp2,
		                       &tmp_vp1) );	
		/* 結果の出力 */
		for(ii2=0;ii2<node.size;ii2++){
    		vp->array[ii2].v.x = tmp_vp1.array[ii2].v.x;
    		vp->array[ii2].v.y = tmp_vp1.array[ii2].v.y;
			vp->array[ii2].v.z = tmp_vp1.array[ii2].v.z;
			vp->array[ii2].s   = tmp_vp1.array[ii2].s;
		}
	}
	RC_TRY( free_fem_velo_array(&tmp_vp1) );
	RC_TRY( free_fem_velo_array(&tmp_vp2) );
	RC_TRY( clean_fem_velo_array(vp) );

	return(NORMAL_RC);
}


/* 流れ場解析の実行とともに動画ファイルを出力する */
RC flow_solve_non_steady_flow(int plot_step, double delta_t, int iteration,
                              double Re, FEM_NODE_ARRAY node, 
                              FEM_ELEMENT_ARRAY element, FEM_BC_ARRAY bc_velo,
                              FEM_BC_ARRAY bc_press, FEM_VELO_ARRAY *vp)
{
	int ii1, ii2;
	int dim ;
	int step = 0;
	FEM_VELO_ARRAY tmp_vp1;
	FEM_VELO_ARRAY tmp_vp2;
	TRANSLATION3D delta_v;
	double delta_s;
	double delta_norm;
	double norm;
	FILE *fp;

	RC_TRY( allocate_fem_velo_array(node.size, &tmp_vp1) );
	RC_TRY( allocate_fem_velo_array(node.size, &tmp_vp2) );
	RC_TRY( allocate_fem_velo_array(node.size, vp) );

	for(ii1=0;ii1<node.size;ii1++){
		tmp_vp1.array[ii1].node = node.array[ii1].label;
		tmp_vp2.array[ii1].node = node.array[ii1].label;
		vp->array[ii1].node = node.array[ii1].label;
	}
	tmp_vp1.size = node.size;
	tmp_vp2.size = node.size;
	vp->size = node.size;

	dim = fem_analysis_dim(element); 
	RC_TRY( flow_property_print(dim, iteration, Re, delta_t) );

	fp = fopen_avs_step_data("./result_movie.inp", node, element);
	RC_NULL_CHK(fp);

	for(ii1=0; ii1<iteration; ii1++){
		
		/* 中間流速の計算 */
		RC_TRY( solve_middle_velocity(delta_t, Re, node, element, bc_velo,
		                              tmp_vp1, &tmp_vp2) );
		/* 圧力の計算 */
		RC_TRY( solve_pressure(dim, delta_t, node, element, &tmp_vp2,
		                       bc_press) ); 
		/* 流速の計算 */
		RC_TRY( solve_velocity(delta_t, node, element, bc_velo, tmp_vp2,
		                       &tmp_vp1) );	

		delta_norm = norm = 0.0;
		for(ii2=0; ii2<node.size; ii2++){
			delta_v.x = vp->array[ii2].v.x - tmp_vp1.array[ii2].v.x;
			delta_v.y = vp->array[ii2].v.y - tmp_vp1.array[ii2].v.y;
			delta_v.z = vp->array[ii2].v.z - tmp_vp1.array[ii2].v.z;
			delta_s = vp->array[ii2].s - tmp_vp1.array[ii2].s;
			delta_norm += delta_v.x * delta_v.x
			            + delta_v.y * delta_v.y
			            + delta_v.z * delta_v.z
			            + delta_s * delta_s;
			norm += tmp_vp1.array[ii2].v.x * tmp_vp1.array[ii2].v.x
			      + tmp_vp1.array[ii2].v.y * tmp_vp1.array[ii2].v.y
			      + tmp_vp1.array[ii2].v.z * tmp_vp1.array[ii2].v.z
			      + tmp_vp1.array[ii2].s * tmp_vp1.array[ii2].s;
		}
		fprintf(stderr, "%d : norm = %e ",
		                ii1, ( sqrt(delta_norm)/sqrt(norm) )/node.size );

		/* 結果の出力 */
		for(ii2=0; ii2<node.size; ii2++){
			vp->array[ii2].v.x = tmp_vp1.array[ii2].v.x;
			vp->array[ii2].v.y = tmp_vp1.array[ii2].v.y;
			vp->array[ii2].v.z = tmp_vp1.array[ii2].v.z;
			vp->array[ii2].s   = tmp_vp1.array[ii2].s;
		}

		if( (ii1 == 0) || (step == plot_step) ){
			RC_TRY( output_avs_step_velo(fp, *vp) );
			step = 0;
		}	

		step++;		
		fprintf(stderr, ".");
	}
	fprintf(stderr,".....Done\n");
	RC_TRY( free_fem_velo_array(&tmp_vp1) );
	RC_TRY( free_fem_velo_array(&tmp_vp2) );
	RC_TRY( clean_fem_velo_array(vp) );
	fclose(fp);

	return(NORMAL_RC);
}


/* 中間流速の計算 */
static RC solve_middle_velocity(double delta_t, double Re, 
                                FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                                FEM_BC_ARRAY bc_velo, FEM_VELO_ARRAY vp,
                                FEM_VELO_ARRAY *vp1)
{
	int dim = 0;
	int index;
	int ii1, ii2;
	double SS;
	static int init_flag = 0;
	static G_POINT_PROP g_prop;
	static double *AMB = NULL;
	static double **ADV = NULL;
	static double **VIS = NULL;
	static double **V1 = NULL;
	static double **V2 = NULL;

	if(init_flag == 0){
		RC_TRY( allocate2D(3, FEM_MAX_NODE, &ADV) );
		RC_TRY( allocate2D(3, FEM_MAX_NODE, &VIS) );
		RC_TRY( allocate2D(3, FEM_MAX_NODE, &V1) );
		RC_TRY( allocate2D(3, node.size, &V2) );
		RC_TRY( allocate_g_prop(&g_prop) );
		AMB = (double *)malloc(node.size * sizeof(double));
		init_flag = 1;
	}

	init_matrix(3, node.size, V2);

	for(ii1=0; ii1<node.size; ii1++){
		AMB[ii1] = 0.0;
	}

	RC_TRY( total_element_volume(element, node, &SS) );

	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		dim = element_dim(element.array[ii1].type);
		if(dim <= 1){
			fprintf(stderr, "element_dim < ");
			return(ARG_ERROR_RC);
		}
		
		/* 構造体 g_prop に各種パラメータを入力 */
		RC_TRY( set_g_prop(element.array[ii1], node, &g_prop) );

		/* 各要素ごとの流速を代入 */
		init_matrix(dim, element.array[ii1].node_num, V1);
		for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
			index = search_fem_node_label(node, element.array[ii1].node[ii2]);
			RC_NEG_CHK(index);
			AMB[index] += 1.0;
			V1[0][ii2] = vp.array[index].v.x;
			V1[1][ii2] = vp.array[index].v.y;
			if(dim == 3) V1[2][ii2] = vp.array[index].v.z;
		}

		RC_TRY( set_advection_term(dim, node, element.array[ii1], g_prop, V1,
		                           ADV) );
		RC_TRY( set_viscosity_term(dim, element.array[ii1], g_prop, V1, VIS) );

		/* 全体の流速を代入 */
		for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
			index = search_fem_node_label(node, element.array[ii1].node[ii2]);
			RC_NEG_CHK(index);
			V2[0][index] += (ADV[0][ii2] + VIS[0][ii2]/Re);
			V2[1][index] += (ADV[1][ii2] + VIS[1][ii2]/Re);
			if(dim == 3){
				V2[2][index] += (ADV[2][ii2] + VIS[2][ii2]/Re);
			}
		}

	}

	for(ii1=0; ii1<node.size; ii1++){
		if(AMB[ii1] < ABS_ERROR) return(CAL_ERROR_RC);
		SS = delta_t / AMB[ii1];
		vp1->array[ii1].v.x = vp.array[ii1].v.x - V2[0][ii1] * SS;
		vp1->array[ii1].v.y = vp.array[ii1].v.y - V2[1][ii1] * SS;
		if(dim == 3) vp1->array[ii1].v.z = vp.array[ii1].v.z - V2[2][ii1] * SS;
	}

	/* 中間流速の境界条件 */
	for(ii1=0; ii1<bc_velo.size; ii1++){
		index = search_fem_node_label(node, bc_velo.array[ii1].node);
		RC_NEG_CHK(index);
		if(bc_velo.array[ii1].v_type.x == BC_VELO){
			vp1->array[index].v.x = bc_velo.array[ii1].v.x ;
		} 
		if(bc_velo.array[ii1].v_type.y == BC_VELO){
			vp1->array[index].v.y = bc_velo.array[ii1].v.y ;
		}
		if(bc_velo.array[ii1].v_type.z == BC_VELO){
			vp1->array[index].v.z = bc_velo.array[ii1].v.z ;
		}
	}
/*
	for(ii1=0; ii1<node.size; ii1++){
		init_fem_disp_velo( &(vp.array[ii1]) );
	}
*/

	return(NORMAL_RC);
}


/* 対流項の計算 */
static RC set_advection_term(int dim, FEM_NODE_ARRAY node, FEM_ELEMENT element, 
                             G_POINT_PROP g_prop, double **V, double **ADV)
{
	int ii1, ii2, ii3, ii4;
	static int init_flag = 0;
	static double **K = NULL;
	static double **M = NULL;
	static double **inverse_M = NULL;
	static double **inverse_MK = NULL;
	
	if(init_flag == 0){
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &K) );
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &M) );
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &inverse_M) );
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &inverse_MK) );
		init_flag = 1;
	}

	RC_TRY( set_Mass_matrix(element.node_num, g_prop, M) );
	RC_TRY( inverse_matrix(element.node_num, M, inverse_M) );

	init_matrix(dim, element.node_num, ADV);

	for(ii1=0; ii1<dim; ii1++){
		for(ii2=0; ii2<dim; ii2++){
			RC_TRY( set_K_matrix(ii2, ii1, dim, node, element, g_prop, V, K) );
			init_matrix(element.node_num, element.node_num, inverse_MK);
			mul_matrix(element.node_num, element.node_num, element.node_num,
			           inverse_M, K, inverse_MK);
			for(ii3=0; ii3<element.node_num; ii3++){
				for(ii4=0; ii4<element.node_num; ii4++){
					ADV[ii1][ii3] += inverse_MK[ii3][ii4] * V[ii2][ii4];
				}
			}
		}
	}

	return(NORMAL_RC);
}


/* 粘性項の計算 */
static RC set_viscosity_term(int dim, FEM_ELEMENT element, G_POINT_PROP g_prop, 
                             double **V, double **VIS)
{
	int ii1, ii2, ii3, ii4;
	static int init_flag = 0;
	static double **S = NULL;
	static double **M = NULL;
	static double **inverse_M = NULL;
	static double **inverse_MS = NULL;

	if(init_flag == 0){
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &S) );
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &M) );
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &inverse_M) );
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &inverse_MS) );
		init_flag = 1;
	}

	RC_TRY( set_Mass_matrix(element.node_num, g_prop, M) );
	RC_TRY( inverse_matrix(element.node_num, M, inverse_M) );

	init_matrix(dim, element.node_num, VIS);

	for(ii1=0; ii1<dim; ii1++){
		for(ii2=0; ii2<dim; ii2++){
			RC_TRY( set_S_matrix(ii1, ii2, dim, element, g_prop, S) );
			init_matrix(element.node_num, element.node_num, inverse_MS);
			mul_matrix(element.node_num, element.node_num, element.node_num,
			           inverse_M, S, inverse_MS);
			for(ii3=0; ii3<element.node_num; ii3++){
				for(ii4=0; ii4<element.node_num; ii4++){
					VIS[ii1][ii3] += inverse_MS[ii3][ii4] * V[ii2][ii4];
				}
			}
		}
	}

	return(NORMAL_RC);
}


/* ポアソン方程式右辺（荷重項）の計算 */
static RC set_f_vector_for_pressure(int dim, double delta_t, 
                                    FEM_NODE_ARRAY node,
                                    FEM_ELEMENT_ARRAY element, 
                                    FEM_VELO_ARRAY vp1, double *f_vector)
{
	int index1, index2;
	int ii1, ii2, ii3, ii4;
	static double **V = NULL;
	static double **H = NULL;
	static G_POINT_PROP g_prop;
	static int init_flag = 0;

	if(init_flag == 0){
		RC_TRY( allocate2D(3, FEM_MAX_NODE, &V) );
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &H) );
		RC_TRY( allocate_g_prop(&g_prop) );
		init_flag = 1;
	}

	for(ii1=0; ii1<node.size; ii1++){
		f_vector[ii1] = 0.0;
	}
	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;
		/* 構造体 g_prop に各種パラメータを入力 */
		RC_TRY( set_g_prop(element.array[ii1], node, &g_prop) );

		for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
			index1 = search_fem_node_label(node, element.array[ii1].node[ii2]);
			RC_NEG_CHK(index1);
			V[0][ii2] = vp1.array[index1].v.x;
			V[1][ii2] = vp1.array[index1].v.y;
			V[2][ii2] = vp1.array[index1].v.z;
		}

		for(ii2=0; ii2<dim; ii2++){
			RC_TRY( set_H_matrix(ii2, element.array[ii1], g_prop, H) );
			for(ii3=0; ii3<element.array[ii1].node_num; ii3++){
				index2 = search_fem_node_label(node,
				                               element.array[ii1].node[ii3]);
				RC_NEG_CHK(index2);
				for(ii4=0; ii4<element.array[ii1].node_num; ii4++){
					f_vector[index2] -= (H[ii3][ii4] * V[ii2][ii4]) / delta_t;
				}
			}
		}
	}

	return(NORMAL_RC);
}


/* 圧力解析用要素剛性マトリックスの作成 */
/* index_array : 全体剛性マトリックス組み込み用インデックス */
/* matrix_size : 要素剛性マトリックスのサイズ */
/* elem_matrix : 要素剛性マトリックス */
static RC set_elem_matrix_for_pressure(FEM_ELEMENT element, FEM_NODE_ARRAY node,
                                       int *index_array, int *matrix_size,
                                       double **elem_matrix)
{
	int dim = 0;
	int ii1, ii2, ii3, ii4;
	static double **dNxy = NULL;
	static int init_flag = 0;
	static G_POINT_PROP g_prop;

	dim = element_dim(element.type);
	if(dim <= 1){
		fprintf(stderr, "element_dim < ");
		return(ARG_ERROR_RC);
	}

	if(init_flag == 0){
		RC_TRY( allocate2D(dim, FEM_MAX_NODE, &dNxy) );
		RC_TRY( allocate_g_prop(&g_prop) );
		init_flag = 1;
	}

	/* elem_matrix を初期化 */
	/* index_array を初期化 */
	for(ii1=0; ii1<element.node_num; ii1++){
		int index;

		RC_NEG_CHK( index = search_fem_node_label(node, element.node[ii1]) );
		index_array[ii1] = index;

		for(ii2=0; ii2<element.node_num; ii2++){
			elem_matrix[ii1][ii2] = 0.0;
		}
	}

	*matrix_size = element.node_num;

	/* 構造体 g_prop に各種パラメータを入力 */
	RC_TRY( set_g_prop(element, node, &g_prop) );

	for(ii1=0; ii1<g_prop.num_point; ii1++){
		dNxy = g_prop.dNxy[ii1];
		for(ii2=0; ii2<dim; ii2++){
			for(ii3=0; ii3<element.node_num; ii3++){
				for(ii4=0; ii4<element.node_num; ii4++){
					elem_matrix[ii3][ii4] += dNxy[ii2][ii3] * dNxy[ii2][ii4]
					                       * g_prop.det_J[ii1]
					                       * g_prop.weight[ii1];
				}
			}
		}
	}

	return(NORMAL_RC);
}


/* 圧力の計算 */
static RC solve_pressure(int dim, double delta_t, FEM_NODE_ARRAY node,
                         FEM_ELEMENT_ARRAY element, FEM_VELO_ARRAY *vp1,
                         FEM_BC_ARRAY bc_press)
{
	int ii1;
	int matrix_size;
	int index_array[FEM_MAX_NODE];
	double *f_vector;
	static double **elem_matrix = NULL;
	int p_flag = 1;
	static int init_flag = 0;
	FEM_SKYLINE_MATRIX g_matrix;


	if(init_flag == 0){
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &elem_matrix) );
		init_flag = 1;
	}

	RC_NEG_ZERO_CHK(node.size);
	RC_NULL_CHK( f_vector = (double *)malloc( node.size*sizeof(double) ) );
	RC_TRY( set_f_vector_for_pressure(dim, delta_t, node, element, *vp1,
	                                  f_vector) );

	/* g_matrix の確保と初期化 */
	RC_TRY( allocate_matrix(1, 1, node, element, &g_matrix) );

	/* 要素マトリックスを全体マトリックスに代入*/
	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		RC_TRY( set_elem_matrix_for_pressure(element.array[ii1], node,
		                                     index_array, &matrix_size,
                                             elem_matrix) );
		RC_TRY( imp_g_matrix(1, matrix_size, index_array, elem_matrix,
		                     g_matrix) );
	}

	/* 境界条件の考慮 */	
	RC_TRY( flow_modify_matrix(1, node, bc_press, g_matrix, f_vector) );

	RC_TRY( s_cholesky_decomp(g_matrix.dof, g_matrix.index1, g_matrix.index2,
	                          g_matrix.array,p_flag) );

	RC_TRY( s_cholesky_subst(g_matrix.dof, g_matrix.index1, g_matrix.index2,
	                         g_matrix.array, f_vector,p_flag) );

	for(ii1=0; ii1<node.size; ii1++){
		vp1->array[ii1].s = f_vector[ii1];
	}
	free( (void *)f_vector );

	RC_TRY( free_g_matrix(1, &g_matrix) );

	return(NORMAL_RC);
}


/* 流速の計算 */
static RC solve_velocity(double delta_t, FEM_NODE_ARRAY node,
                         FEM_ELEMENT_ARRAY element, FEM_BC_ARRAY bc_velo,
                         FEM_VELO_ARRAY vp1, FEM_VELO_ARRAY *vp)
{
	int dim = 0;
	int index1, index2, index3;
	int ii1, ii2, ii3, ii4;
	double SS;
	static int init_flag = 0;
	static double *AMB = NULL;
	static double **V = NULL;
	static double **H = NULL;
	static double **M = NULL;
	static double **inverse_M = NULL;
	static double **inverse_MH = NULL;
	static G_POINT_PROP g_prop;

	if(init_flag == 0){
		RC_TRY( allocate2D(3, node.size, &V) );
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &H) );
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &M) );
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &inverse_M) );
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &inverse_MH) );
		RC_TRY( allocate_g_prop(&g_prop) );
		AMB = (double *)malloc(node.size * sizeof(double));
		init_flag = 1;
	}

	init_matrix(3, node.size, V);

	for(ii1=0; ii1<node.size; ii1++){
		AMB[ii1] = 0.0;
	}

	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		dim = element_dim(element.array[ii1].type);
		if(dim <= 1){
			fprintf(stderr, "element_dim < ");
			return(ARG_ERROR_RC);
		}

		/* 構造体 g_prop に各種パラメータを入力 */
		RC_TRY( set_g_prop(element.array[ii1], node, &g_prop) );
		RC_TRY( set_Mass_matrix(element.array[ii1].node_num, g_prop, M) );
		RC_TRY( inverse_matrix(element.array[ii1].node_num, M, inverse_M) );

		for(ii2=0; ii2<dim; ii2++){
			RC_TRY( set_H_matrix(ii2, element.array[ii1], g_prop, H) );
			mul_matrix(element.array[ii1].node_num, element.array[ii1].node_num,
			           element.array[ii1].node_num, inverse_M, H, inverse_MH);
			for(ii3=0; ii3<element.array[ii1].node_num; ii3++){
				index1 = search_fem_node_label(node,
				                               element.array[ii1].node[ii3]);
				RC_NEG_CHK(index1);
				for(ii4=0; ii4<element.array[ii1].node_num; ii4++){
					index2 = search_fem_node_label(node,
					                              element.array[ii1].node[ii4]);
					RC_NEG_CHK(index2);
					V[ii2][index1] += inverse_MH[ii3][ii4]*vp1.array[index2].s;
				}
			}
		}
		for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
			index3 = search_fem_node_label(node, element.array[ii1].node[ii2]);
			RC_NEG_CHK(index3);
			AMB[index3] += 1.0;
		}
	}

	/* 流速の代入計算 */
	for(ii1=0; ii1<node.size; ii1++){
		if(AMB[ii1] < ABS_ERROR) return(CAL_ERROR_RC);
		SS = delta_t / AMB[ii1];
		vp->array[ii1].v.x = vp1.array[ii1].v.x - V[0][ii1] * SS;
		vp->array[ii1].v.y = vp1.array[ii1].v.y - V[1][ii1] * SS;
		if(dim == 3) vp->array[ii1].v.z = vp1.array[ii1].v.z - V[2][ii1] * SS;
	}

	/* 流速の境界条件 */
	for(ii1=0; ii1<bc_velo.size; ii1++){
		index3 = search_fem_node_label(node, bc_velo.array[ii1].node);
		RC_NEG_CHK(index3);
		if(bc_velo.array[ii1].v_type.x == BC_VELO){
			vp->array[index3].v.x = bc_velo.array[ii1].v.x ;
		}
		if(bc_velo.array[ii1].v_type.y == BC_VELO){
			vp->array[index3].v.y = bc_velo.array[ii1].v.y ;
		}	
		if(bc_velo.array[ii1].v_type.z == BC_VELO){
			vp->array[index3].v.z = bc_velo.array[ii1].v.z ;
		}
	}

	/* 圧力の代入 */
	for(ii1=0; ii1<node.size; ii1++){
		vp->array[ii1].s = vp1.array[ii1].s;
	}
/*
	for(ii1=0; ii1<node.size; ii1++){
		init_fem_disp_velo( &(vp1.array[ii1]) );
	}
*/

	return(NORMAL_RC);
}



/* 随伴解析---------------------------------------------------------------*/
/* 流速・圧力パラメータより随伴解析を行い、随伴変数を渡す */
RC flow_adjoint_analysis(double Re, FEM_NODE_ARRAY node,
                         FEM_ELEMENT_ARRAY element, FEM_BC_ARRAY bc,
                         FEM_VELO_ARRAY vp, FEM_VELO_ARRAY *adjoint_vp)
{
	int ii1;
	int dim;
	int matrix_size;
	int index_array[4*FEM_MAX_NODE];
	double *f_vector;
	static double **elem_matrix = NULL;
	int p_flag = 1;
	static int init_flag = 0;
	FEM_SKYLINE_MATRIX g_matrix;
	FEM_ELEMENT_ARRAY press_element;
	FILE *fp;

	RC_TRY( allocate1D(4*node.size, &f_vector) );
	dim = fem_analysis_dim(element);

	if(init_flag == 0){
		RC_TRY( allocate2D(4*FEM_MAX_NODE, 4*FEM_MAX_NODE, &elem_matrix) );
		RC_TRY( allocate_fem_element_array(element.size, &press_element) );
		RC_TRY( allocate_fem_velo_array(node.size, adjoint_vp) );
		init_flag = 1;
	}

	if( (fp = fopen("./test2.dat","w") ) == NULL){
		fprintf(stderr, "open_file < ");
	}

	for(ii1=0;ii1<node.size;ii1++){
		adjoint_vp->array[ii1].node = node.array[ii1].label;
	}
	adjoint_vp->size = node.size;

	RC_TRY( make_press_element_data(element, &press_element) );

	RC_TRY( set_f_vector_for_adjoint_analysis(dim, node, element, bc, vp,
	                                          f_vector) );

	RC_TRY( allocate_matrix(2, dim+1, node, element, &g_matrix) );

	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;
		RC_TRY( set_elem_matrix_for_adjoint_analysis(dim, Re,
		                                             element.array[ii1], 
		                                             press_element.array[ii1],
		                                             node, vp, index_array, 
		                                             &matrix_size,
		                                             elem_matrix) );
		RC_TRY( imp_g_matrix(2, matrix_size, index_array, elem_matrix,
		                     g_matrix) );
	}

	RC_TRY( flow_modify_matrix(3, node, bc, g_matrix, f_vector) );

	RC_TRY( pivotting_check(node, g_matrix, f_vector) );

	RC_TRY( s_crout_decomp(g_matrix.dof, g_matrix.index1, g_matrix.index2,
	                       g_matrix.array, p_flag) );
	RC_TRY( s_crout_subst(g_matrix.dof, g_matrix.index1, g_matrix.index2,
	                      g_matrix.array, f_vector,p_flag) );

	for(ii1=0; ii1<node.size; ii1++){
		switch(dim){
		case 2:
			adjoint_vp->array[ii1].v.x = f_vector[ii1*g_matrix.dim];
			adjoint_vp->array[ii1].v.y = f_vector[ii1*g_matrix.dim + 1];
			adjoint_vp->array[ii1].s   = f_vector[ii1*g_matrix.dim + 2];
			break;
		case 3:
			adjoint_vp->array[ii1].v.x = f_vector[ii1*g_matrix.dim];
			adjoint_vp->array[ii1].v.y = f_vector[ii1*g_matrix.dim + 1];
			adjoint_vp->array[ii1].v.z = f_vector[ii1*g_matrix.dim + 2];
			adjoint_vp->array[ii1].s   = f_vector[ii1*g_matrix.dim + 3];
			break;
		default:
			return(IMPLEMENT_ERROR_RC);
		}
	}

	RC_TRY( free_g_matrix(2, &g_matrix) );
	fclose(fp);
	free(f_vector);

	return(NORMAL_RC);
}


/* 流速・圧力パラメータより随伴解析を行い、随伴変数を渡す */
RC flow_adjoint_analysis2(double Re, FEM_NODE_ARRAY node,
                          FEM_ELEMENT_ARRAY element, FEM_BC_ARRAY bc,
                          FEM_VELO_ARRAY vp, FEM_VELO_ARRAY *adjoint_vp)
{
	int ii1, ii2, ii3;
	int dim;
	int ix1, ix2;
	int matrix_size;
	static int *index;
	int index_array[4*FEM_MAX_NODE];
	double dum;
	double *f_vector;
	double **K_matrix = NULL;
	static double **elem_matrix = NULL;
	static int init_flag = 0;
	FEM_ELEMENT_ARRAY press_element;
	FILE *fp;

	dim = fem_analysis_dim(element);
	RC_TRY( allocate1D(4*node.size, &f_vector) );

	if(init_flag == 0){
		index = (int *)malloc(4*node.size * sizeof(int));
		RC_TRY( allocate2D(4*node.size, 4*node.size, &K_matrix) );
		RC_TRY( allocate2D(4*FEM_MAX_NODE, 4*FEM_MAX_NODE, &elem_matrix) );
		RC_TRY( allocate_fem_element_array(element.size, &press_element) );
		RC_TRY( allocate_fem_velo_array(node.size, adjoint_vp) );
		init_flag = 1;
	}

	if( (fp = fopen("./test2.dat","w") ) == NULL){
		fprintf(stderr, "open_file < ");
	}
	for(ii1=0;ii1<node.size;ii1++){
		adjoint_vp->array[ii1].node = node.array[ii1].label;
	}
	adjoint_vp->size = node.size;

	RC_TRY( make_press_element_data(element, &press_element) );

	RC_TRY( set_f_vector_for_adjoint_analysis(dim, node, element, bc, vp,
	                                          f_vector) );

	for(ii1=0; ii1<element.size ; ii1++){
		if(element.array[ii1].label < 0) continue;
		RC_TRY( set_elem_matrix_for_adjoint_analysis(dim, Re,
		                                             element.array[ii1], 
		                                             press_element.array[ii1],
		                                             node, vp, index_array, 
		                                             &matrix_size,
	                                                 elem_matrix) );
		for(ii2=0; ii2<matrix_size; ii2++){
			for(ii3=0; ii3<matrix_size; ii3++){
				ix1 = index_array[ii2];
				ix2 = index_array[ii3];
				K_matrix[ix1][ix2] += elem_matrix[ii2][ii3];
			}
		}
	}
	RC_TRY( adjoint_modify_matrix(dim, (dim+1)*node.size, node, bc, K_matrix,
	                              f_vector) );

	RC_TRY( lu_decomp(K_matrix, (dim+1)*node.size, index, &dum) );
	RC_TRY( lu_subst(K_matrix, (dim+1)*node.size, index, f_vector) );

	RC_TRY( free2D(4*node.size, 4*node.size, K_matrix) );
	
	for(ii1=0; ii1<node.size; ii1++){
		switch(dim){
		case 2:
			adjoint_vp->array[ii1].v.x = f_vector[ii1*(dim+1)];
			adjoint_vp->array[ii1].v.y = f_vector[ii1*(dim+1) + 1];
			adjoint_vp->array[ii1].s   = f_vector[ii1*(dim+1) + 2];
			break;
		case 3:
			adjoint_vp->array[ii1].v.x = f_vector[ii1*(dim+1)];
			adjoint_vp->array[ii1].v.y = f_vector[ii1*(dim+1) + 1];
			adjoint_vp->array[ii1].v.z = f_vector[ii1*(dim+1) + 2];
			adjoint_vp->array[ii1].s   = f_vector[ii1*(dim+1) + 3];
			break;
		default:
			return(IMPLEMENT_ERROR_RC);
		}
	}

	fclose(fp);
	free(f_vector);

	return(NORMAL_RC);
}


/* 2次元のみ */
static RC set_f_vector_for_adjoint_analysis(int dim, FEM_NODE_ARRAY node, 
                                            FEM_ELEMENT_ARRAY element,
                                            FEM_BC_ARRAY bc, FEM_VELO_ARRAY vp,
                                            double *f_vector)
{
	int ii1, ii2, ii3;
	int index1, index2;
	int bc_index1, bc_index2;
	int check_flag;
	int check_index[4];
	double AJ, BJ, Up, Vp;
	double check_param;

	for(ii1=0; ii1<(dim+1)*node.size; ii1++){
		f_vector[ii1] = 0.0;
	}

	for(ii1=0; ii1<element.size; ii1++){
		check_flag = 0;
		check_param = 0.0;
		for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
			for(ii3=0; ii3<bc.size; ii3++){
				if(element.array[ii1].node[ii2] == bc.array[ii3].node){
					check_index[check_flag] = bc.array[ii3].node;
					check_param += bc.array[ii3].v.x
					             + bc.array[ii3].v.y
					             + bc.array[ii3].v.z;
					check_flag++;
				}
			}
		}
		if( (check_flag > 1) && (check_param > ABS_ERROR) && (dim == 2) ){
			index1 = search_fem_node_label(node, check_index[0]);
			index2 = search_fem_node_label(node, check_index[1]);
			bc_index1 =  search_fem_bc_node(bc, check_index[0]);
			bc_index2 =  search_fem_bc_node(bc, check_index[1]);
			AJ = fabs(node.array[index1].p.x - node.array[index2].p.x);
			BJ = fabs(node.array[index1].p.y - node.array[index2].p.y);
			Up = (bc.array[bc_index1].v.x + bc.array[bc_index2].v.x) / 2.0;
			Vp = (bc.array[bc_index1].v.y + bc.array[bc_index2].v.y) / 2.0;
			f_vector[(dim+1)*index1 + (dim)] +=  (Up * BJ) - (Vp * AJ);
			f_vector[(dim+1)*index2 + (dim)] +=  (Up * BJ) - (Vp * AJ);
			fprintf(stderr, "%d %d : \n", check_index[0] ,check_index[1]);
			fprintf(stderr, "%f %f %f %f\n", AJ, BJ, Up, Vp);
		}
	}

	return(NORMAL_RC);
}


static RC set_elem_matrix_for_adjoint_analysis(int dim, double Re,
                                               FEM_ELEMENT v_element,
                                               FEM_ELEMENT p_element,
                                               FEM_NODE_ARRAY node,
                                               FEM_VELO_ARRAY vp,
                                               int *index_array,
                                               int *matrix_size,
                                               double **elem_matrix)
{
	int ii1, ii2, ii3, ii4;
	int g_dim = 0;
	int n_num = 0;
	static double **V = NULL;
	static double **K = NULL;
	static double **S = NULL;
	static double **H = NULL;
	static double **H1 = NULL;
	static double **H2 = NULL;
	static int init_flag = 0;
	static G_POINT_PROP v_g_prop;
	static G_POINT_PROP p_g_prop;

	if(init_flag == 0){
		RC_TRY( allocate2D(3, FEM_MAX_NODE, &V) );
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &K) );
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &S) );
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &H) );
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &H1) );
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &H2) );
		RC_TRY( allocate_g_prop(&v_g_prop) );
		RC_TRY( allocate_g_prop(&p_g_prop) );
		init_flag = 1;
	}

	g_dim = dim + 1;
	n_num = v_element.node_num;

	/* elem_matrix を初期化 */
	for(ii1=0; ii1 < v_element.node_num; ii1++){
		int index;
		RC_NEG_CHK( index = search_fem_node_label(node, v_element.node[ii1]) );
		V[0][ii1] = vp.array[index].v.x;
		V[1][ii1] = vp.array[index].v.y;
		V[2][ii1] = vp.array[index].v.z;

		for(ii2=0; ii2<g_dim; ii2++){
			index_array[g_dim*ii1 + ii2] = g_dim*index + ii2;
			for(ii3=0; ii3<g_dim*v_element.node_num; ii3++){
				elem_matrix[g_dim*ii1 + ii2][ii3] = 0.0;
			}
		}
	}

	*matrix_size = g_dim * v_element.node_num;
	RC_TRY( set_g_prop(v_element, node, &v_g_prop) );
	RC_TRY( set_g_prop(p_element, node, &p_g_prop) );

	for(ii1=0; ii1<dim; ii1++){
        RC_TRY( set_H1_matrix(ii1, v_element, p_element, v_g_prop, p_g_prop,
		                      H1) );
        RC_TRY( set_H2_matrix(ii1, v_element, p_element, v_g_prop, p_g_prop,
		                      H2) );

		transpose_matrix(p_element.node_num, v_element.node_num, H1, H2);

		for(ii2=0; ii2<p_element.node_num; ii2++){
			for(ii3 = 0; ii3<v_element.node_num; ii3++){
				elem_matrix[ii2*g_dim][ii3*g_dim + ii1] += H1[ii2][ii3];
        	}
		}
		for(ii2=0; ii2<v_element.node_num; ii2++){
			for(ii3=0; ii3<p_element.node_num; ii3++){
				elem_matrix[ii2*g_dim+1 + ii1][ii3*g_dim+dim] += H2[ii2][ii3];
			}
		}

		RC_TRY( set_H_matrix(ii1, v_element, v_g_prop, H) );
		for(ii2=0; ii2<v_element.node_num; ii2++){
       		for(ii3=0; ii3<v_element.node_num; ii3++){
				elem_matrix[ii2*g_dim][ii3*g_dim + ii1] += H[ii2][ii3];
				elem_matrix[ii2*g_dim+1 + ii1][ii3*g_dim+dim] += H[ii2][ii3];
			}
		}

		for(ii2=0; ii2<dim; ii2++){
			RC_TRY( set_K_matrix(ii2, ii1, dim, node, v_element, v_g_prop,
			                     V, K) );
			RC_TRY( set_S_matrix(ii1, ii2, dim, v_element, v_g_prop, S) );

			for(ii3 = 0; ii3<v_element.node_num; ii3++){
				for(ii4=0; ii4<v_element.node_num; ii4++){
					elem_matrix[ii3*g_dim + 1 + ii1][ii4*g_dim + ii2]
					           += 2*K[ii3][ii4] + (S[ii3][ii4])/Re;
				}
			}
		}
	}

	return(NORMAL_RC);
}


static RC adjoint_modify_matrix(int dim, int matrix_size, FEM_NODE_ARRAY node,
                                FEM_BC_ARRAY bc, double **array,
                                double *vector)
{
	int ii1, ii2;
	int index;

	for(ii1=0; ii1<bc.size; ii1++){
		if(bc.array[ii1].node < 0) continue;

		index = search_fem_node_label(node, bc.array[ii1].node);
		RC_NEG_CHK(index);

		if(bc.array[ii1].v_type.x == BC_VELO){
			for(ii2=0; ii2<matrix_size; ii2++){
				array[index][ii2] = 0.0;
				array[ii2][index] = 0.0;
			}
			array[index][index] = 1.0;
			vector[index] = 0.0;
		}
	
		if(bc.array[ii1].v_type.y == BC_VELO){
			for(ii2=0; ii2<matrix_size; ii2++){
				array[index+1][ii2] = 0.0;
				array[ii2][index+1] = 0.0;
			}
			array[index+1][index+1] = 1.0;
			vector[index+1] = 0.0;
		}

		if(dim == 3){
			if(bc.array[ii1].v_type.z == BC_VELO){
				for(ii2=0; ii2<matrix_size; ii2++){
						array[index+2][ii2] = 0.0;
					array[ii2][index+2] = 0.0;
				}
				array[index+2][index+2] = 1.0;
				vector[index+2] = 0.0;
			}
		}

		if(bc.array[ii1].s_type[0] == BC_PRESS){
			for(ii2=0; ii2<matrix_size; ii2++){
				array[index+dim][ii2] = 0.0;
				array[ii2][index+dim] = 0.0;
			}
			array[index+dim][index+dim] = 1.0;
			vector[index+dim] = 0.0;
		}
	}

	return(NORMAL_RC);
}


/* 感度導出用ツール-------------------------------------------------------*/

/* ひずみ速度の計算 */
RC strain_velocity(int fem_sol_ss_flag, FEM_NODE_ARRAY node,
                   FEM_ELEMENT_ARRAY element, FEM_VELO_ARRAY velo,
                   FEM_STRAIN_ARRAY *strain_velo)
{
	int ii1, ii2;
	static int init_flag = 0;
	static double **local_points;

	if(init_flag == 0){
		RC_TRY( allocate2D(FEM_MAX_NODE, 4, &local_points) );
		init_flag = 1;
	}

	RC_TRY( allocate_fem_strain_array(0, strain_velo) );
	strain_velo->source = FEM_K_SOL;

	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		/* 要素中心 */
		if(fem_sol_ss_flag & FEM_SOL_SS_CENTER){
			RC_TRY( realloc_fem_strain_array(strain_velo) );
			RC_TRY( set_center_point(element.array[ii1].type,
			                         local_points[0]) );
			RC_TRY( local_velocity_strain( element.array[ii1], local_points[0],
                                           node, velo,
                                   &(strain_velo->array[strain_velo->size]) ) );
			strain_velo->array[strain_velo->size].node = -1;  /* center */
			(strain_velo->size)++;
		}

		/* 各節点 */
		if(fem_sol_ss_flag & FEM_SOL_SS_NODE){
			RC_TRY( set_local_node_points(element.array[ii1].type,
			                              local_points) );

			for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
				RC_TRY( realloc_fem_strain_array(strain_velo) );
				RC_TRY( local_velocity_strain( element.array[ii1],
				                               local_points[ii2], node, velo,
				                   &(strain_velo->array[strain_velo->size]) ) );
				strain_velo->array[strain_velo->size].node
				    = element.array[ii1].node[ii2];
				(strain_velo->size)++;
			}
		}
	}

	RC_TRY( clean_fem_strain_array(strain_velo) );

	return(NORMAL_RC);
}


/* 要素の局所ひずみ速度 */
static RC local_velocity_strain(FEM_ELEMENT element, const double *local_point,
                                FEM_NODE_ARRAY node, FEM_VELO_ARRAY velo, 
                                FEM_STRESS_STRAIN *strain_velo)
{
	int dim, ssnum;
	double velo_v[3*FEM_MAX_NODE];
	double strain_velo_v[6];
	double det_J;
	static double **node_location = NULL;
	static double **B_matrix = NULL;
	static int init_flag = 0;

	if(element.label < 0) return(ARG_ERROR_RC);

	if(init_flag == 0){
		RC_TRY( allocate2D(FEM_MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(6, 3*FEM_MAX_NODE, &B_matrix) );
		init_flag = 1;
	}

	if( (dim = element_dim(element.type) ) <= 1){
		fprintf(stderr, "element_dim < ");
		return(ARG_ERROR_RC);
	}

	if( (ssnum = element_ssnum(element.type)) <= 0){
		fprintf(stderr, "element_ssnum < ");
		return(ARG_ERROR_RC);
	}

	init_fem_stress_strain(strain_velo);
	strain_velo->element = element.label;

	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( make_B_matrix(element, local_point, node_location, B_matrix,
	                      &det_J) );
	RC_TRY( fill_velo_vector(element, velo, velo_v) );	

	mul_matrix_vect(ssnum, dim*element.node_num, B_matrix, velo_v,
	                strain_velo_v); 
	RC_TRY( vect2stress_strain(ssnum, strain_velo_v, strain_velo) );

    return(NORMAL_RC);
}


static RC fill_velo_vector(FEM_ELEMENT element,
                           FEM_VELO_ARRAY velo, double *velo_v)
{
	int dim;
	int ii1;

	dim = element_dim(element.type);
	if(dim <= 1){
		fprintf(stderr, "element_dim < ");
		return(ARG_ERROR_RC);
	}

	for(ii1=0; ii1<element.node_num; ii1++){
		int index;

		index = search_fem_velo_node(velo, element.node[ii1]);
		RC_NEG_CHK(index);

		velo_v[dim*ii1 + 0] = velo.array[index].v.x;
		velo_v[dim*ii1 + 1] = velo.array[index].v.y;
		if(dim == 3){
			velo_v[dim*ii1 + 2] = velo.array[index].v.z;
		}
	}

	return(NORMAL_RC);
}


/* 解析用ツール --------------------------------------------------------*/

/* 全境界節点から壁面境界節点を取り出す */
RC pick_out_wall_bc_node(FEM_BC_ARRAY bc, FEM_NODE_ARRAY *wall_node)
{
	int ii1;
	static int velo_flag = 0;

	RC_TRY( allocate_fem_node_array(0, wall_node) );
	wall_node->size = 0;

	for(ii1=0;ii1<bc.size;ii1++){
		velo_flag = 0;
		if(bc.array[ii1].node < 0) continue;

		if((bc.array[ii1].v_type.x == BC_VELO)&&(bc.array[ii1].v.x>ABS_ERROR)){
			velo_flag = 1;
		}
		if((bc.array[ii1].v_type.y == BC_VELO)&&(bc.array[ii1].v.y>ABS_ERROR)){
			velo_flag = 1;
		}
		if((bc.array[ii1].v_type.z == BC_VELO)&&(bc.array[ii1].v.z>ABS_ERROR)){
			velo_flag = 1;
		}

		if( (velo_flag != 1) && (bc.array[ii1].s_type[0] == BC_FREE) ){
			RC_TRY( realloc_fem_node_array(wall_node) );
			wall_node->array[wall_node->size].label = bc.array[ii1].node;
			(wall_node->size)++;
		}
	}
	
	return(NORMAL_RC);
}


/* 境界条件のタイプを替える(2次要素モデルの解析に用いる) */
/* ANSYS-FLOTRANでは流れ場の2次要素が用意されていないため、
 * 構造解析用のモデルを用いて「流速」を「強制変位」、
 *                           「圧力」に「境界温度」を代わりに与えている。
 */
RC change_bc_type4flow(FEM_BC_ARRAY *bc)
{
	int ii1;

	for(ii1=0; ii1<bc->size; ii1++){
		if(bc->array[ii1].node < 0) continue;
		if(bc->array[ii1].v_type.x == BC_FIX) bc->array[ii1].v_type.x = BC_VELO;
		if(bc->array[ii1].v_type.y == BC_FIX) bc->array[ii1].v_type.y = BC_VELO;
		if(bc->array[ii1].v_type.z == BC_FIX) bc->array[ii1].v_type.z = BC_VELO;
		if(bc->array[ii1].s_type[0] == BC_TEMP){
			bc->array[ii1].s_type[0] = BC_PRESS;
		} 
	}

	return(NORMAL_RC);
}


/* g_matrixの解放 */
/* mode : 1. コレスキー　2. クラウト*/
static RC free_g_matrix(int mode, FEM_SKYLINE_MATRIX *g_matrix)
{
	switch(mode){
	case 1:
		RC_TRY( s_cholesky_free(g_matrix->index2, g_matrix->array) );
		break;
	case 2:
		RC_TRY( s_crout_free(g_matrix->index2, g_matrix->array) );
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	free(g_matrix->index1);
	g_matrix->dim = 0;
	g_matrix->dof = 0;
	g_matrix->index1 = NULL;
	g_matrix->index2 = NULL;
	g_matrix->array = NULL;

	return(NORMAL_RC);
}


/* G_PROPへ各種値を格納*/
static RC set_g_prop(FEM_ELEMENT element, FEM_NODE_ARRAY node,
                     G_POINT_PROP *g_prop)
{
	int ii1, ii2;
	static double **node_location = NULL;
	static double **dN = NULL;
	static double **jacobi = NULL;
	static double **inverse_J = NULL;
	static double N[FEM_MAX_NODE];
	static int array_flag = 0;
	int dim;

	/* 静的変数の初期化 */
	if(array_flag == 0){
		RC_TRY( allocate2D(FEM_MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(3, 3, &jacobi) );
		RC_TRY( allocate2D(3, FEM_MAX_NODE, &dN) );
		RC_TRY( allocate2D(3, 3, &inverse_J) );
		array_flag = 1;
	}

	if( (dim = element_dim(element.type) ) <= 1){
		fprintf(stderr, "element_dim < ");
		return(ARG_ERROR_RC);
	}

	RC_TRY( set_Gauss_const(element.type, g_prop) );
	RC_TRY( node_location_matrix(&element, node, node_location) );

	for(ii1=0; ii1<g_prop->num_point; ii1++){
		init_matrix(dim, FEM_MAX_NODE, dN);
		RC_TRY( set_N_vector(element.type, g_prop->g_points[ii1], N) );
		RC_TRY( set_dN_matrix(element.type, g_prop->g_points[ii1], dN) );

		/* dN * xy => jacobi */
		/* jacobi^-1 => inverse_J */
		/* inverse_J * dN => dNxy */
		mul_matrix(dim, element.node_num, dim, dN, node_location, jacobi);
		RC_TRY( inverse_matrix(dim, jacobi, inverse_J) );
		mul_matrix(dim, dim, element.node_num, inverse_J, dN,
		           g_prop->dNxy[ii1]);

		/* N, det[J] */
		for(ii2=0; ii2<element.node_num; ii2++){
        	g_prop->N[ii1][ii2] = N[ii2];
		}
		g_prop->det_J[ii1]  = determinant(dim, jacobi);
	}

	return(NORMAL_RC);
}


/* 拘束条件を考慮して全体剛性マトリックス、荷重ベクトルを修正 */
/* mode : 1. 圧力解析用, 2, 随伴解析用 */
static RC flow_modify_matrix(int mode, FEM_NODE_ARRAY node,
                             FEM_BC_ARRAY bc_press,
                             FEM_SKYLINE_MATRIX g_matrix, double *f_vector)
{
	int ii1;
	int index;

	for(ii1=0; ii1<bc_press.size; ii1++){
		if(bc_press.array[ii1].node < 0) continue;

		index = search_fem_node_label(node, bc_press.array[ii1].node);
		RC_NEG_CHK(index);
		switch(mode){
		case 1:
			if(bc_press.array[ii1].s_type[0] == BC_PRESS){
				RC_TRY( s_cholesky_mod( g_matrix.dof, g_matrix.index1,
				                        g_matrix.index2, g_matrix.array,
				                        f_vector, g_matrix.dim*index,
				                        bc_press.array[ii1].s[0]) );
			}
			break;
/*
		case 2:
			if(bc.array[ii1].v_type.x == BC_VELO){
				RC_TRY( s_crout_mod(g_matrix.dof, g_matrix.index1,
				                    g_matrix.index2, g_matrix.array, f_vector,
				                    g_matrix.dim*index, 0.0) );
			}
			if(bc.array[ii1].v_type.y == BC_VELO){
				RC_TRY( s_crout_mod(g_matrix.dof, g_matrix.index1,
				                    g_matrix.index2, g_matrix.array, f_vector,
				                    g_matrix.dim*index+1, 0.0) );
			}
			if( (g_matrix.dim-1) == 3 ){
				if(bc.array[ii1].v_type.z == BC_VELO){
					RC_TRY( s_crout_mod(g_matrix.dof, g_matrix.index1,
				                        g_matrix.index2, g_matrix.array,
					                    f_vector, g_matrix.dim*index+2, 0.0) );
				}
			}
			if(bc.array[ii1].s_type[0] == BC_PRESS){
				RC_TRY( s_crout_mod( g_matrix.dof, g_matrix.index1,
				                     g_matrix.index2, g_matrix.array,f_vector,
				                     g_matrix.dim*index + g_matrix.dim-1,
				                     bc.array[ii1].s[0]) );
			}
			break;
*/
		default:
			return(IMPLEMENT_ERROR_RC);
		}
	}

	return(NORMAL_RC);
}


/* 全体剛性マトリックスの確保と問題の次元をチェック */
/* mode : 1. コレスキー, 2. クラウト */
static RC allocate_matrix(int mode, int dim, FEM_NODE_ARRAY node,
                          FEM_ELEMENT_ARRAY element,
                          FEM_SKYLINE_MATRIX *g_matrix)
{

	/* 問題の次元(g_matrix->dim) */
	RC_NEG_ZERO_CHK(dim);

	g_matrix->dim = dim;
	g_matrix->dof = (g_matrix->dim)*node.size;
	g_matrix->index1 = (long *)malloc(g_matrix->dof * sizeof(long));
	RC_NULL_CHK(g_matrix->index1);
	RC_TRY( fill_index1(g_matrix->dim, node, element, g_matrix->index1) );

	switch(mode){
	case 1:
		RC_TRY( s_cholesky_alloc(g_matrix->dof, g_matrix->index1,
		                         &(g_matrix->index2),
		                         &(g_matrix->array), &(g_matrix->array_size)) );
		break;
	case 2:
		RC_TRY( s_crout_alloc(g_matrix->dof, g_matrix->index1,
		                      &(g_matrix->index2),
		                      &(g_matrix->array), &(g_matrix->array_size)) );
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC fill_index1(int dim, FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                      long *index1)
{
	int ii1, ii2, ii3, ii4;
	long pos1, pos2;

	for(ii1=0; ii1<(dim*node.size); ii1++){
		index1[ii1] = ii1;
	}

	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
			pos1 = dim * search_fem_node_label(node,
			                                   element.array[ii1].node[ii2]);
			RC_NEG_CHK(pos1);
			for(ii3=0; ii3<element.array[ii1].node_num; ii3++){
				pos2 = dim * search_fem_node_label(node,
				                                  element.array[ii1].node[ii3]);
				RC_NEG_CHK(pos2);
				for(ii4=0; ii4<dim; ii4++){
					if(pos1 < index1[pos2+ii4]) index1[pos2+ii4] = pos1;
				}
			}
		}
	}

	return(NORMAL_RC);
}


/* 要素剛性マトリックスを全体剛性マトリックスに組み込む */
static RC imp_g_matrix(int mode, int matrix_size, const int *index_array,
                       double **elem_matrix, FEM_SKYLINE_MATRIX g_matrix)
{
	RC rc;
	int ii1, ii2;
	int ix1, ix2;
	double *ptr;

	for(ii1=0; ii1<matrix_size; ii1++){
		for(ii2=0; ii2<matrix_size; ii2++){
			ix1 = index_array[ii1];
			ix2 = index_array[ii2];

			switch(mode){
			case 1:
				if(ix1 < ix2) continue;
				ptr = s_cholesky_ptr(g_matrix.dof, g_matrix.index1,
				                     g_matrix.index2, g_matrix.array,
				                     ix1, ix2, &rc);
				break;
			case 2:
	        	ptr = s_crout_ptr(g_matrix.dof, g_matrix.index1,
				                  g_matrix.index2, g_matrix.array,
				                  ix1, ix2, &rc);
				break;
			default:
				return(IMPLEMENT_ERROR_RC);
			}
	
			if(rc != NORMAL_RC){
				fprintf(stderr, "ptr < ");
				return(rc);
			}
			(*ptr) += elem_matrix[ii1][ii2];
        }
	}

	return(NORMAL_RC);
}


/* 圧力用要素データの作成 */
static RC make_press_element_data(FEM_ELEMENT_ARRAY velo_elem,
                                  FEM_ELEMENT_ARRAY *press_elem)
{
	int ii1, ii2;

	press_elem->size = velo_elem.size;

	for(ii1=0; ii1<velo_elem.size; ii1++){
		press_elem->array[ii1].label = velo_elem.array[ii1].label;

		if(velo_elem.array[ii1].type == ELEM_TRI2){
			press_elem->array[ii1].type  = ELEM_TRI1;
			press_elem->array[ii1].node_num = 3;
		}else if(velo_elem.array[ii1].type == ELEM_QUAD2){
			press_elem->array[ii1].type = ELEM_QUAD1;
			press_elem->array[ii1].node_num = 4;
		}else if(velo_elem.array[ii1].type == ELEM_TETRA2){
			press_elem->array[ii1].type = ELEM_TETRA1;
			press_elem->array[ii1].node_num = 4;
		}else if(velo_elem.array[ii1].type == ELEM_HEXA2){
			press_elem->array[ii1].type = ELEM_HEXA1;
			press_elem->array[ii1].node_num = 8;
		}

		for(ii2=0; ii2<press_elem->array[ii1].node_num; ii2++){
			press_elem->array[ii1].node[ii2] = velo_elem.array[ii1].node[ii2];
		}
	}

	return(NORMAL_RC);
}


static RC pivotting_check(FEM_NODE_ARRAY node, FEM_SKYLINE_MATRIX g_matrix,
                          double *vector)
{
	int ii1, ii2, ii3, ii4, ii5;
	int g_dim;
	double get_param;
	double *ptr;
	RC rc;

	g_dim = g_matrix.dim;

	/* エラーのチェック */
	for(ii1=0;ii1<node.size;ii1++){
		for(ii2=0;ii2<g_dim;ii2++){
			if(fabs(g_matrix.array[g_matrix.index2[ii1*g_dim+ii2]])<ABS_ERROR){
				get_param = 0;
				for(ii3=0;ii3<g_dim;ii3++){
					get_param += s_crout_get(g_matrix.dof, g_matrix.index1,
					                         g_matrix.index2, g_matrix.array,
					                         ii1*g_dim+ii3, ii1*g_dim+ii2, &rc);
				}
				if(fabs(get_param) < ABS_ERROR){
					fprintf(stderr,"%d %d\n ", ii1+1, ii2+1);

					for(ii4=0;ii4<g_dim;ii4++){
						for(ii5=0;ii5<g_dim;ii5++){
							fprintf(stderr,"%f ", s_crout_get(g_matrix.dof,
							                                  g_matrix.index1,
							                                  g_matrix.index2,
							                                  g_matrix.array,
							                                  ii1*g_dim+ii4,
							                                  ii1*g_dim+ii5,
							                                  &rc));
                		}
                		fprintf(stderr,"\n ");
        			}

					fprintf(stderr," pivotting error\n");
					/* return(CAL_ERROR_RC); */
				}
			}
		}
	}

	/* ゼロ部の消去 */
	for(ii1=0;ii1<node.size;ii1++){
		for(ii2=0;ii2<g_dim;ii2++){
			if(fabs(g_matrix.array[g_matrix.index2[ii1*g_dim+ii2]])<ABS_ERROR){
				for(ii3=0;ii3<g_dim;ii3++){
					get_param = s_crout_get(g_matrix.dof, g_matrix.index1,
					                        g_matrix.index2, g_matrix.array,
					                        ii1*g_dim+ii3, ii1*g_dim+ii2, &rc);
					if(fabs(get_param) > ABS_ERROR){
						for(ii4=0;ii4<g_matrix.dof;ii4++){
							ptr = s_crout_ptr(g_matrix.dof, g_matrix.index1,
							                  g_matrix.index2, g_matrix.array,
							                  ii1*g_dim+ii2, ii4, &rc);
							(*ptr) += s_crout_get(g_matrix.dof, g_matrix.index1,
							                      g_matrix.index2,
							                      g_matrix.array, ii1*g_dim+ii3,
							                      ii4, &rc);
						}
						vector[ii1*g_dim+ii2] += vector[ii1*g_dim+ii3];
						break;
					}
				}
			}
		}
	}

	return(NORMAL_RC);
}


/* 質量行列の計算 */
static RC set_Mass_matrix(int node_num, G_POINT_PROP g_prop, double **M)
{
	int ii1, ii2, ii3;

	init_matrix(node_num, node_num, M);

	for(ii1=0; ii1<g_prop.num_point; ii1++){
		for(ii2=0; ii2<node_num; ii2++){
			for(ii3=0; ii3<node_num; ii3++){
				M[ii2][ii3] += g_prop.N[ii1][ii2] * g_prop.N[ii1][ii3] 
				             * g_prop.det_J[ii1] * g_prop.weight[ii1];
			}
		}
	}

	return(NORMAL_RC);
}


/* 移流項行列の計算 */
static RC set_K_matrix(int dN_xyz, int v_xyz, int dim, FEM_NODE_ARRAY node,
                       FEM_ELEMENT element, G_POINT_PROP g_prop, double **V,
                       double **K)
{
	int ii1, ii2, ii3;
	double sum_dNV;
	static double **dN;
	static double **NN;
	static int init_flag = 0;

	if(init_flag == 0){
		RC_TRY( allocate2D(FEM_MAX_NODE, FEM_MAX_NODE, &NN) );
		RC_TRY( allocate2D(3, FEM_MAX_NODE, &dN) );
		init_flag = 1;
	}

	init_matrix(element.node_num, element.node_num, K);

	for(ii1=0; ii1<g_prop.num_point; ii1++){
		sum_dNV = 0.0;
		dN = g_prop.dNxy[ii1];
		for(ii2=0; ii2<element.node_num; ii2++){
			for(ii3=0; ii3<element.node_num; ii3++){
				NN[ii2][ii3] = g_prop.N[ii1][ii2] * g_prop.N[ii1][ii3];
			}
			sum_dNV += dN[dN_xyz][ii2] * V[v_xyz][ii2];
		}
		for(ii2=0; ii2<element.node_num; ii2++){
			for(ii3=0; ii3<element.node_num; ii3++){
     			K[ii2][ii3] += sum_dNV * NN[ii2][ii3] * g_prop.det_J[ii1]
				                                      * g_prop.weight[ii1];
			}
		}
	}

	return(NORMAL_RC);
}


/* 粘性項行列の計算 */
static RC set_S_matrix(int dN_xyz1, int dN_xyz2, int dim, FEM_ELEMENT element, 
                       G_POINT_PROP g_prop, double **S)
{
	int ii1, ii2, ii3;
	static double **dNxy;
	static int init_flag = 0;

	if(init_flag == 0){
		RC_TRY( allocate2D(dim, element.node_num, &dNxy) );
		init_flag = 1;
	}

	init_matrix(element.node_num, element.node_num, S);

	for(ii1=0; ii1<g_prop.num_point; ii1++){
		dNxy = g_prop.dNxy[ii1];
		for(ii2=0; ii2<element.node_num; ii2++){
			for(ii3=0; ii3<element.node_num; ii3++){
				if(dN_xyz1 == dN_xyz2){
					S[ii2][ii3] += 2 * dNxy[dN_xyz1][ii2] * dNxy[dN_xyz2][ii3] 
					             * g_prop.det_J[ii1] * g_prop.weight[ii1];
				}else{
					S[ii2][ii3] += dNxy[dN_xyz1][ii2] * dNxy[dN_xyz2][ii3] 
					             * g_prop.det_J[ii1] * g_prop.weight[ii1];
				}
			}
		}
	}

	return(NORMAL_RC);
}


/* ポアソン方程式右辺行列の計算 */
static RC set_H_matrix(int dN_xyz, FEM_ELEMENT element, G_POINT_PROP g_prop,
                       double **H)
{
	int ii1, ii2, ii3;
	static double **dNxy;
	static int init_flag = 0;

	if(init_flag == 0){
		RC_TRY( allocate2D(3, FEM_MAX_NODE, &dNxy) );
		init_flag = 1;
	}

	init_matrix(element.node_num, element.node_num, H);

	for(ii1=0; ii1<g_prop.num_point; ii1++){
		dNxy = g_prop.dNxy[ii1];
		for(ii2=0; ii2<element.node_num; ii2++){
	        for(ii3=0; ii3<element.node_num; ii3++){
				H[ii2][ii3] += g_prop.N[ii1][ii2] * dNxy[dN_xyz][ii3]
				             * g_prop.det_J[ii1] * g_prop.weight[ii1];
        	}
		}
	}

	return(NORMAL_RC);
}


static RC set_H1_matrix(int dN_xyz, FEM_ELEMENT v_elem, FEM_ELEMENT p_elem,
                        G_POINT_PROP v_g_prop, G_POINT_PROP p_g_prop,
                        double **H1)
{
	int ii1, ii2, ii3;
	static double **dNxy;
	static int init_flag = 0;

	if(init_flag == 0){
		RC_TRY( allocate2D(3, FEM_MAX_NODE, &dNxy) );
		init_flag = 1;
	}

	init_matrix(p_elem.node_num, v_elem.node_num, H1);

	for(ii1=0; ii1<v_g_prop.num_point; ii1++){
		dNxy = p_g_prop.dNxy[ii1];
		for(ii2=0; ii2<p_elem.node_num; ii2++){

			/*fprintf(stderr, "%f", p_g_prop.N[ii1][ii2]);*/
				for(ii3=0; ii3<v_elem.node_num; ii3++){
        			H1[ii2][ii3] += v_g_prop.N[ii1][ii3] * dNxy[dN_xyz][ii2]
					              * v_g_prop.det_J[ii1] * v_g_prop.weight[ii1];
				}

		}
		/*fprintf(stderr, "\n");*/
	} 

	return(NORMAL_RC);
}


static RC set_H2_matrix(int dN_xyz, FEM_ELEMENT v_elem, FEM_ELEMENT p_elem,
                        G_POINT_PROP v_g_prop, G_POINT_PROP p_g_prop,
                        double **H2)
{
	int ii1, ii2, ii3;
	static double **dNxy;
	static int init_flag = 0;

	if(init_flag == 0){
		RC_TRY( allocate2D(3, FEM_MAX_NODE, &dNxy) );
		init_flag = 1;
	}

	init_matrix(v_elem.node_num, p_elem.node_num, H2);

	for(ii1=0; ii1<v_g_prop.num_point; ii1++){
		dNxy = p_g_prop.dNxy[ii1];
		for(ii2=0; ii2<v_elem.node_num; ii2++){
			for(ii3=0; ii3<p_elem.node_num; ii3++){
				H2[ii2][ii3] += v_g_prop.N[ii1][ii2] * dNxy[dN_xyz][ii3]
			                  * v_g_prop.det_J[ii1] * v_g_prop.weight[ii1];
			}
		}
	}
							  
	return(NORMAL_RC);
}


static RC set_Gauss_const(FEM_ELEM_TYPE type, G_POINT_PROP *g_prop)
{
	switch(type){
	case ELEM_TRI1:
		g_prop->num_point = 4;
		RC_TRY( Gauss_const_tri(g_prop->num_point, g_prop->g_points,
		                        g_prop->weight) );
		break;
	case ELEM_TRI2:
		g_prop->num_point = 7;
		RC_TRY( Gauss_const_tri(g_prop->num_point, g_prop->g_points,
		                        g_prop->weight) );
		break;
	case ELEM_QUAD1:
		g_prop->num_point = 9;
		RC_TRY( Gauss_const_quad(g_prop->num_point, g_prop->g_points,
		                         g_prop->weight) );
		break;
	case ELEM_QUAD2:
		g_prop->num_point = 9;	
		RC_TRY( Gauss_const_quad(g_prop->num_point, g_prop->g_points,
		                         g_prop->weight) );
		break;
	case ELEM_TETRA1:
	case ELEM_TETRA2:
		g_prop->num_point = 4;
		RC_TRY( Gauss_const_tetra(g_prop->num_point, g_prop->g_points,
		                          g_prop->weight) );
		break;
	case ELEM_PENTA1:
		g_prop->num_point = 6;
		RC_TRY( Gauss_const_penta(g_prop->num_point, g_prop->g_points,
		                          g_prop->weight) );
		break;
	case ELEM_PENTA2:
		g_prop->num_point = 9;
		RC_TRY( Gauss_const_penta(g_prop->num_point, g_prop->g_points,
		                          g_prop->weight) );
		break;
	case ELEM_HEXA1:
		g_prop->num_point = 8;
		RC_TRY( Gauss_const_hexa(g_prop->num_point, g_prop->g_points,
		                         g_prop->weight) );
		break;
	case ELEM_HEXA2:
		g_prop->num_point = 14;
		RC_TRY( Gauss_const_hexa(g_prop->num_point, g_prop->g_points,
		                         g_prop->weight) );
		break;
	default :
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


/* G_PROP の確保 */
static RC allocate_g_prop(G_POINT_PROP *g_prop)
{
	int ii1;

	RC_TRY( allocate2D(MAX_GAUSS_POINT, 4, &(g_prop->g_points)) );

	for(ii1=0; ii1<MAX_GAUSS_POINT; ii1++){
		RC_TRY( allocate2D(3, FEM_MAX_NODE, &(g_prop->dNxy[ii1]) ) );
	}

	return(NORMAL_RC);
}


static RC flow_property_print(int dim, int iteration, double Re, double delta_t)
{
	fprintf(stdout," DIMENSION : %7d\n", dim);
	fprintf(stdout," ITERATION : %7d\n", iteration);
	fprintf(stdout," DELTA_T   : %7.3e\n", delta_t);
	fprintf(stdout," RE NUMBER : %7.3e\n", Re);

	return(NORMAL_RC);
}  
