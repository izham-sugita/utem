/*********************************************************************
 * fem_solver_plate_shell.c
 *
 * Copyright (C) 2012 AzLib Developers Group
 *
 * Written by
 *  <Taiki AOYAMA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: fem_solver_plate_shell.c 1022 2012-04-30 05:00:38Z aoyama $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "base.h"
#include "mathematics.h"
#include "fem_solver.h"
#include "fem_struct.h"
#include "sky_cholesky.h"
#include "nst_component.h"
  

static RC set_B_matrix_plate_membrane(int node_num, double **dNxyz,
                                      double **B_matrix);
static RC set_B_matrix_plate_bend(int node_num, double **dNxyz,
                                  double **B_matrix);
static RC set_B_matrix_plate_shear(int node_num, double *N, double **dNxyz,
                                   double **B_matrix);
static RC local_disp_vect_shell(double L[][3], int node_num,
                                double disp_vect[]);
static RC set_drilling_dof(int mitc_node_num, ELEMENT *element,
                           NODE_ARRAY node, MATERIAL_PROP_ARRAY material,
                           PHYSICAL_PROP_ARRAY physical, double thickness,
                           ELEM_MATRIX *elem_matrix);
static RC set_N_matrix_membrane(int node_num, double *N, double **N_matrix);
static RC set_N_matrix_bend(int node_num, double *N, double **N_matrix);

/*
 * 平面シェル要素 微小変形解析ソルバー
 * 面内変形と曲げ変形は非連成
 * 
 * Reissner-Mindlin板曲げ理論使用
 * ロッキング回避のため剪断剛性のみ次数低減積分実行
 *
 * 連立方程式の求解にはICCG使用するが，板厚が薄くなると正定値性が弱くなるので
 * fem_solve_shellの使用を推奨する．
 */
RC
fem_solve_plate_shell_iccg6 (NODE_ARRAY node, ELEMENT_ARRAY element,
                             RIGID_ELEMENT_ARRAY rigid,
                             MATERIAL_PROP_ARRAY material,
                             PHYSICAL_PROP_ARRAY physical,
                             BC_ARRAY rest, BC_ARRAY force,
                             BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                             DEFAULT_BC_ARRAY b_force,
                             DISP_ARRAY *disp, STRESS_ARRAY *stress,
                             STRAIN_ARRAY *strain,
                             int sol_ss_flag, int sol_ss_face)
{
	int ii1;
	int dim;
	double *f_vect;
	double *u_vect;
	NONZERO_MATRIX3 matrix;
	LOCAL_COORD_ARRAY dummy_coord;
	WORK_SET ws;


	if((dim = analysis_dim_shell(element)) != 6){
		return(ARG_ERROR_RC);
	}
	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

	RC_TRY( allocate_nonzero6(node, element, rigid, &matrix) );

	RC_TRY( allocate_work_set(&ws) );
	for(ii1=0; ii1<element.size; ii1++){
		ELEM_MATRIX elem_matrix;

		if(element.array[ii1].label < 0) continue;

		RC_TRY( make_plate_shell_elem_matrix(element.array[ii1], node,
		                                     material, physical,
		                                     &elem_matrix, ws) );
		RC_TRY( impose_nonzero6(matrix, elem_matrix, 1.0) );

		RC_TRY( free_elem_matrix(&elem_matrix) );
	}
	RC_TRY( free_work_set(&ws) );

	RC_TRY( fem_allocate_vector(dim, node, &f_vect) );
	/* 節点荷重・モーメント */
	RC_TRY( add_nodal_force(dim, node, force, f_vect) );
	/* 熱荷重 未対応 */
	/* 重力荷重 未対応 */

	RC_TRY( allocate_local_coord_array(0, &dummy_coord) );
	RC_TRY( modify_nonzero6(node, dummy_coord, rest, rigid, matrix,
	                        &f_vect, 1, NULL) );
	RC_TRY( free_local_coord_array(&dummy_coord) );
	RC_TRY( scale_iccg_nonzero3(matrix, f_vect, &u_vect, NULL) );
	RC_TRY( free_nonzero_matrix3(&matrix) );
	RC_TRY( fem_free_vector(&f_vect) );
	RC_TRY( rigid_recover6(node, rigid, u_vect) );
	RC_TRY( vect2disp(dim, node, u_vect, disp) );
	RC_TRY( fem_free_vector(&u_vect) );

	/* 応力，ひずみ */
	RC_TRY( stress_strain_plate_shell(sol_ss_flag, sol_ss_face, node, element,
	                                  material, physical, *disp,
	                                  stress, strain) );

	return(NORMAL_RC);
}


/*
 * 平面シェル要素 微小変形解析ソルバー
 * 面内変形と曲げ変形は非連成
 * 
 * Reissner-Mindlin板曲げ理論使用
 * ロッキング回避のため剪断剛性のみ次数低減積分実行
 *
 * 連立方程式の求解にはスカイライン法使用
 */
RC
fem_solve_plate_shell (NODE_ARRAY node, ELEMENT_ARRAY element,
                       RIGID_ELEMENT_ARRAY rigid,
                       MATERIAL_PROP_ARRAY material,
                       PHYSICAL_PROP_ARRAY physical,
                       BC_ARRAY rest, BC_ARRAY force,
                       BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                       DEFAULT_BC_ARRAY b_force,
                       DISP_ARRAY *disp,
                       STRESS_ARRAY *stress,
                       STRAIN_ARRAY *strain, int sol_ss_flag, int sol_ss_face)
{
	int ii1;
	int dim;
	double *f_vect;
	SKYLINE_MATRIX matrix;
	WORK_SET ws;


	if((dim = analysis_dim_shell(element)) != 6){
		return(ARG_ERROR_RC);
	}
	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

	RC_TRY( allocate_g_matrix(dim, node, element, &matrix) );

	RC_TRY( allocate_work_set(&ws) );
	for(ii1=0; ii1<element.size; ii1++){
		ELEM_MATRIX elem_matrix;

		if(element.array[ii1].label < 0) continue;

		RC_TRY( make_plate_shell_elem_matrix(element.array[ii1], node,
		                                     material, physical,
		                                     &elem_matrix, ws) );
		RC_TRY( impose_elem_matrix(matrix, elem_matrix, 1.0) );

		RC_TRY( free_elem_matrix(&elem_matrix) );
	}
	RC_TRY( free_work_set(&ws) );

	RC_TRY( fem_allocate_vector(dim, node, &f_vect) );
	/* 節点荷重・モーメント */
	RC_TRY( add_nodal_force(dim, node, force, f_vect) );
	/* 熱荷重 (未対応) */
	/* 重力荷重 (未対応) */

	RC_TRY( modify_g_matrix6(node, rest, matrix, f_vect) );
	RC_TRY( decomp_g_matrix(matrix) );
	RC_TRY( s_cholesky_subst(matrix.dof, matrix.index1, matrix.index2,
	                         matrix.array, f_vect, 1) );
	RC_TRY( vect2disp(dim, node, f_vect, disp) );
	RC_TRY( free_g_matrix(&matrix) );
	RC_TRY( fem_free_vector(&f_vect) );

	/* 応力，ひずみ */
	RC_TRY( stress_strain_plate_shell(sol_ss_flag, sol_ss_face, node, element,
	                                  material, physical, *disp,
	                                  stress, strain) );

	return(NORMAL_RC);
}


/*
 * 要素剛性マトリクス
 *
 * Reissner-Mindlin板曲げ理論使用
 * 剪断剛性は次数低減積分実行
 */
RC
make_plate_shell_elem_matrix (ELEMENT element, NODE_ARRAY node,
                              MATERIAL_PROP_ARRAY material,
                              PHYSICAL_PROP_ARRAY physical,
                              ELEM_MATRIX *elem_matrix,
                              WORK_SET ws)
{
	int ii1, ii2, ii3;
	int dim;
	int ssnum;
	int g_num_point;
	int prop_index;
	int m_index;
	double thickness;
	double det_J;
	double *g_weights = ws.vec_G[0];
	double **g_points = ws.mat_G_4[0];
	double **node_location = ws.mat_N_3[0];
	double **BtDB = ws.mat_3N_3N[0];
	double L[3][3];
		

	if(element.label < 0) return(ARG_ERROR_RC);
	prop_index = search_physical_prop_label(physical, element.physical);
	if(prop_index < 0) return(SEEK_ERROR_RC);
	if(NST_PROP_PSHELL != physical.array[prop_index].i_info[1]){
		return(ARG_ERROR_RC);
	}

	if( (dim = element_dim_shell(element.type)) <= 0
	|| (ssnum = element_ssnum_shell(element.type)) <= 0 ){
		return(ARG_ERROR_RC);
	}

	elem_matrix->dim = dim;
	elem_matrix->element = element.label;
	elem_matrix->size = (element.node_num) * dim;

	RC_TRY( allocate1I(elem_matrix->size, &(elem_matrix->index)) );
	for(ii1=0; ii1<element.node_num; ii1++){
		int index;

		index = node_renum_index(node, element.node[ii1]);
		if(index < 0 ) return(SEEK_ERROR_RC);

		for(ii2=0; ii2<dim; ii2++){
			elem_matrix->index[ii1*dim+ii2] = index*dim + ii2;
		}
	}

	RC_TRY( allocate2D(elem_matrix->size, elem_matrix->size,
	                   &(elem_matrix->matrix)) );

	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( set_L_matrix_shell_solid(element, node_location, L) );
	RC_TRY( local_node_location(L, element.node_num, node_location) );

	RC_TRY( get_thickness(element, physical, &thickness) );

	m_index = search_material_prop_element(material, physical, &element);
	if(m_index < 0) return(SEEK_ERROR_RC);

	/* 膜剛性(Membrane Stiffness) */
	{
		double **dN = ws.mat_3_N[0];
		double **dNxyz = ws.mat_3_N[1];
		double **Jacobi = ws.mat_3_3[0];
		double **inverse_J = ws.mat_3_3[1];
		double **Dm_matrix = ws.mat_6_6[0];
		double **Bm_matrix = ws.mat_6_3N[0];
		MATERIAL_PROP mat;

		/* Set membrane D matrix */
		mat = material.array[m_index];
		mat.ana_type = ANA_PLANE_STRESS;
		RC_TRY( fill_material(&mat) );
		for(ii1=0;ii1<3;ii1++){
			for(ii2=0;ii2<3;ii2++){
				Dm_matrix[ii1][ii2] = mat.D_matrix[ii1][ii2];
			}
		}

		RC_TRY( Gauss_const_default(element.type, g_points, g_weights,
		                            &g_num_point) );

		for(ii1=0; ii1<g_num_point; ii1++){
			double factor;

			RC_TRY( set_dN_matrix(element.type, g_points[ii1], dN) );
			mul_matrix(2, element.node_num, 2, dN, node_location, Jacobi);
			RC_TRY( inverse_matrix(2, Jacobi, inverse_J) );
			mul_matrix(2, 2, element.node_num, inverse_J, dN, dNxyz);

			det_J = determinant(2, Jacobi);
			if(det_J < 0.0) return(CAL_ERROR_RC);

			RC_TRY( set_B_matrix_plate_membrane(element.node_num, dNxyz,
			                                    Bm_matrix) );

			mul_matrix_AtBA(3,element.node_num*dim,Bm_matrix,Dm_matrix,BtDB);

			factor = g_weights[ii1] * det_J * thickness;
			for(ii2=0; ii2<elem_matrix->size; ii2++){
				for(ii3=0; ii3<elem_matrix->size; ii3++){
					elem_matrix->matrix[ii2][ii3] += BtDB[ii2][ii3] * factor;
				}
			}
		}
	}

	/* 曲げ剛性(Bending Stiffness) */
	{
		double **dN = ws.mat_3_N[0];
		double **dNxyz = ws.mat_3_N[1];
		double **Jacobi = ws.mat_3_3[0];
		double **inverse_J = ws.mat_3_3[1];
		double **Db_matrix = ws.mat_6_6[0];
		double **Bb_matrix = ws.mat_6_3N[0];
		MATERIAL_PROP mat;

		/* Set plate bending D matrix */
		mat = material.array[m_index];
		mat.ana_type = ANA_PLATE_BEND;
		RC_TRY( fill_material(&mat) );
		for(ii1=0;ii1<3;ii1++){
			for(ii2=0;ii2<3;ii2++){
				Db_matrix[ii1][ii2] = mat.D_matrix[ii1][ii2];
			}
		}

		RC_TRY( Gauss_const_default(element.type, g_points, g_weights,
		                            &g_num_point) );

		for(ii1=0; ii1<g_num_point; ii1++){
			double factor;

			RC_TRY( set_dN_matrix(element.type, g_points[ii1], dN) );
			mul_matrix(2, element.node_num, 2, dN, node_location, Jacobi);
			RC_TRY( inverse_matrix(2, Jacobi, inverse_J) );
			mul_matrix(2, 2, element.node_num, inverse_J, dN, dNxyz);
			
			det_J = determinant(2, Jacobi);
			if(det_J < 0.0) return(CAL_ERROR_RC);

			RC_TRY( set_B_matrix_plate_bend(element.node_num, dNxyz,
			                                Bb_matrix) );

			mul_matrix_AtBA(3,element.node_num*dim,Bb_matrix,Db_matrix,BtDB);

			factor = g_weights[ii1]*det_J*(thickness*thickness*thickness/12.0);
			for(ii2=0; ii2<elem_matrix->size; ii2++){
				for(ii3=0; ii3<elem_matrix->size; ii3++){
					elem_matrix->matrix[ii2][ii3] += BtDB[ii2][ii3] * factor;
				}
			}
		}
	}

	/* せん断剛性(Shearing Stiffness) */
	{
		double *N = ws.vec_N[0];
		double **dN = ws.mat_3_N[0];
		double **dNxyz = ws.mat_3_N[1];
		double **Jacobi = ws.mat_3_3[0];
		double **inverse_J = ws.mat_3_3[1];
		double **Ds_matrix = ws.mat_6_6[0];
		double **Bs_matrix = ws.mat_6_3N[0];
		MATERIAL_PROP mat;

		/* Set plate bending D matrix for Shearing */
		mat = material.array[m_index];
		mat.ana_type = ANA_PLATE_SHEAR;
		RC_TRY( fill_material(&mat) );
		for(ii1=0;ii1<2;ii1++){
			for(ii2=0;ii2<2;ii2++){
				Ds_matrix[ii1][ii2] = mat.D_matrix[ii1][ii2];
			}
		}

		/* Selective Reduced Integration */
		RC_TRY( Gauss_const_RI(element.type, g_points, g_weights,
		                       &g_num_point) );
#if 0
		RC_TRY( Gauss_const_default(element.type, g_points, g_weights,
		                            &g_num_point) );
#endif

		for(ii1=0; ii1<g_num_point; ii1++){
			double factor;

			RC_TRY( set_N_vector(element.type, g_points[ii1], N) );
			RC_TRY( set_dN_matrix(element.type, g_points[ii1], dN) );
			mul_matrix(2, element.node_num, 2, dN, node_location, Jacobi);
			RC_TRY( inverse_matrix(2, Jacobi, inverse_J) );
			mul_matrix(2, 2, element.node_num, inverse_J, dN, dNxyz);
			
			det_J = determinant(2, Jacobi);
			if(det_J < 0.0) return(CAL_ERROR_RC);

			RC_TRY( set_B_matrix_plate_shear(element.node_num, N, dNxyz,
			                                 Bs_matrix) );

			mul_matrix_AtBA(2,element.node_num*dim,Bs_matrix,Ds_matrix,BtDB);

			factor = g_weights[ii1] * det_J * thickness;
			for(ii2=0; ii2<elem_matrix->size; ii2++){
				for(ii3=0; ii3<elem_matrix->size; ii3++){
					elem_matrix->matrix[ii2][ii3] += BtDB[ii2][ii3] * factor;
				}
			}
		}
	}

	/* 虚偽ドリリング自由度 */
	RC_TRY( set_drilling_dof(element.node_num, &element, node,
	                         material, physical, thickness, elem_matrix) );

	/* Local coordinates -> Global coordinates */
	RC_TRY( mul_matrix33_LALt(elem_matrix->size, elem_matrix->matrix, L) );


	return(NORMAL_RC);
}


/*
 * Bマトリクス
 *
 * Bm : Membrane
 * Bb : Bending 
 * Bs : Shearing
 * 
 * Bm_matrix,Bb_matrix,Bs_matrix,det_JにNULLを渡すと計算無視
 */
RC
make_B_matrix_plate_shell (ELEMENT element, const double *local_point,
                           NODE_ARRAY node,
                           double **Bm_matrix, double **Bb_matrix,
                           double **Bs_matrix, double *det_J)
{
	static int allocate_flag = 0;
	static double *N = NULL;
	static double **node_location = NULL;
	static double **dN = NULL;
	static double **dNxyz = NULL;
	static double **Jacobi = NULL;
	static double **inverse_J = NULL;
	int dim;
	double L[3][3];


	if(allocate_flag == 0){
		RC_TRY( allocate1D(MAX_NODE, &N) );
		RC_TRY( allocate2D(MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(3, MAX_NODE, &dN) );
		RC_TRY( allocate2D(3, MAX_NODE, &dNxyz) );
		RC_TRY( allocate2D(3, 3, &Jacobi) );
		RC_TRY( allocate2D(3, 3, &inverse_J) );

		allocate_flag = 1;
	}

	if((dim = element_dim(element.type)) != 2) return(ARG_ERROR_RC);

	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( set_L_matrix_shell_solid(element, node_location, L) );
	RC_TRY( local_node_location(L, element.node_num, node_location) );

	RC_TRY( set_N_vector(element.type, local_point, N) );
	RC_TRY( set_dN_matrix(element.type, local_point, dN) );
	mul_matrix(dim, element.node_num, dim, dN, node_location, Jacobi);
	RC_TRY( inverse_matrix(dim, Jacobi, inverse_J) );
	mul_matrix(dim, dim, element.node_num, inverse_J, dN, dNxyz);

	if(det_J != NULL){
		*det_J = determinant(dim, Jacobi);
	}

	if(Bm_matrix != NULL){
		RC_TRY( set_B_matrix_plate_membrane(element.node_num, dNxyz,
		                                    Bm_matrix) );
	}
	if(Bb_matrix != NULL){
		RC_TRY( set_B_matrix_plate_bend(element.node_num, dNxyz, Bb_matrix) );
	}
	if(Bs_matrix != NULL){
		RC_TRY( set_B_matrix_plate_shear(element.node_num, N, dNxyz,
		                                 Bs_matrix) );
	}

	return(NORMAL_RC);
}


static RC
set_B_matrix_plate_membrane (int node_num, double **dNxyz, double **B_matrix)
{
	int ii1;
	int dim;


	dim = 6;
	for(ii1=0;ii1<node_num;ii1++){
		B_matrix[0][dim*ii1]   = dNxyz[0][ii1];
		B_matrix[0][dim*ii1+1] = 0.0;
		B_matrix[0][dim*ii1+2] = 0.0;
		B_matrix[0][dim*ii1+3] = 0.0;
		B_matrix[0][dim*ii1+4] = 0.0;
		B_matrix[0][dim*ii1+5] = 0.0;

		B_matrix[1][dim*ii1]   = 0.0;
		B_matrix[1][dim*ii1+1] = dNxyz[1][ii1];
		B_matrix[1][dim*ii1+2] = 0.0;
		B_matrix[1][dim*ii1+3] = 0.0;
		B_matrix[1][dim*ii1+4] = 0.0;
		B_matrix[1][dim*ii1+5] = 0.0;

		B_matrix[2][dim*ii1]   = dNxyz[1][ii1];
		B_matrix[2][dim*ii1+1] = dNxyz[0][ii1];
		B_matrix[2][dim*ii1+2] = 0.0;
		B_matrix[2][dim*ii1+3] = 0.0;
		B_matrix[2][dim*ii1+4] = 0.0;
		B_matrix[2][dim*ii1+5] = 0.0;
	}

	return(NORMAL_RC);
}


static RC
set_B_matrix_plate_bend (int node_num, double **dNxyz, double **B_matrix)
{
	int ii1;
	int dim;


	dim = 6;
	for(ii1=0; ii1<node_num; ii1++){
		B_matrix[0][dim*ii1+0] = 0.0;
		B_matrix[0][dim*ii1+1] = 0.0;
		B_matrix[0][dim*ii1+2] = 0.0;
		B_matrix[0][dim*ii1+3] = 0.0;
		B_matrix[0][dim*ii1+4] = -dNxyz[0][ii1];
		B_matrix[0][dim*ii1+5] = 0.0;

		B_matrix[1][dim*ii1+0] = 0.0;
		B_matrix[1][dim*ii1+1] = 0.0;
		B_matrix[1][dim*ii1+2] = 0.0;
		B_matrix[1][dim*ii1+3] = dNxyz[1][ii1];
		B_matrix[1][dim*ii1+4] = 0.0;
		B_matrix[1][dim*ii1+5] = 0.0;

		B_matrix[2][dim*ii1+0] = 0.0;
		B_matrix[2][dim*ii1+1] = 0.0;
		B_matrix[2][dim*ii1+2] = 0.0;
		B_matrix[2][dim*ii1+3] = dNxyz[0][ii1];
		B_matrix[2][dim*ii1+4] = -dNxyz[1][ii1];
		B_matrix[2][dim*ii1+5] = 0.0;
	}

	return(NORMAL_RC);
}


static RC
set_B_matrix_plate_shear (int node_num, double *N, double **dNxyz,
                          double **B_matrix)
{
	int ii1;
	int dim;


	dim = 6;
	for(ii1=0; ii1<node_num; ii1++){
		B_matrix[0][dim*ii1+0] = 0.0;
		B_matrix[0][dim*ii1+1] = 0.0;
		B_matrix[0][dim*ii1+2] = dNxyz[0][ii1];
		B_matrix[0][dim*ii1+3] = 0.0;
		B_matrix[0][dim*ii1+4] = N[ii1];
		B_matrix[0][dim*ii1+5] = 0.0;

		B_matrix[1][dim*ii1+0] = 0.0;
		B_matrix[1][dim*ii1+1] = 0.0;
		B_matrix[1][dim*ii1+2] = dNxyz[1][ii1];
		B_matrix[1][dim*ii1+3] = -N[ii1];
		B_matrix[1][dim*ii1+4] = 0.0;
		B_matrix[1][dim*ii1+5] = 0.0;
	}

	return(NORMAL_RC);
}


/* 応力、ひずみ計算 */
RC
stress_strain_plate_shell (int sol_ss_flag, int sol_ss_face,
                           NODE_ARRAY node, ELEMENT_ARRAY element,
                           MATERIAL_PROP_ARRAY material,
                           PHYSICAL_PROP_ARRAY physical,
                           DISP_ARRAY disp,
                           STRESS_ARRAY *stress, STRAIN_ARRAY *strain)
{
	int ii1, ii2;
	int node_index;
	int *count = NULL;
	double **local_points = NULL;
	STRESS_STRAIN local_stress, local_strain;
	VECT3DR *stress_sum = NULL;
	VECT3DR *strain_sum = NULL;
	VECT3DR init_v = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	WORK_SET ws;


	if( ((stress == NULL)&&(strain == NULL)) || (sol_ss_flag == 0) ){
		return(NORMAL_RC);
	}

	RC_TRY( allocate2D(MAX_NODE, 4, &local_points) );
	RC_TRY( allocate_work_set(&ws) );

	if(sol_ss_flag & SOL_SS_NODE_AV){
		stress_sum = (VECT3DR *)mm_alloc((node.size)*sizeof(VECT3DR));
		strain_sum = (VECT3DR *)mm_alloc((node.size)*sizeof(VECT3DR));
		count = (int *)mm_alloc((node.size)*sizeof(int));
		if( (stress_sum == NULL) ||(strain_sum == NULL) ||(count == NULL) ){
			return(ALLOC_ERROR_RC);
		}

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
			RC_TRY( set_center_point_shell(sol_ss_face, element.array[ii1].type,
			                               local_points[0]) );
			RC_TRY( local_stress_strain_plate_shell(element.array[ii1],
			                                        local_points[0], node,
			                                        material, physical, disp,
			                                        &local_stress,
			                                        &local_strain,
			                                        sol_ss_face, ws) );
			if(stress != NULL){
				RC_TRY( realloc_stress_array(stress) );
				stress->array[stress->size] = local_stress;
				stress->array[stress->size].node = -1;
				(stress->size)++;
			}
			if(strain != NULL){
				RC_TRY( realloc_strain_array(strain) );
				strain->array[strain->size] = local_strain;
				strain->array[strain->size].node = -1;
				(strain->size)++;
			}
		}

		/* 各節点 */
		if( (sol_ss_flag & SOL_SS_NODE) || (sol_ss_flag & SOL_SS_NODE_AV) ){
			RC_TRY( set_local_node_points_shell(sol_ss_face,
			                                    element.array[ii1].type,
			                                    local_points) );

			for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
				RC_TRY( local_stress_strain_plate_shell(element.array[ii1],
				                                        local_points[ii2], node,
				                                        material, physical,
				                                        disp, &local_stress,
				                                        &local_strain,
				                                        sol_ss_face, ws) );

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
	RC_TRY( free2D(MAX_NODE, 4, &local_points) );
	RC_TRY( free_work_set(&ws) );

	if(stress != NULL) RC_TRY( clean_stress_array(stress) );
	if(strain != NULL) RC_TRY( clean_strain_array(strain) );

	if(sol_ss_flag & SOL_SS_NODE_AV){
		if(stress != NULL){
			for(ii1=0; ii1<stress->size; ii1++){
				if(stress->array[ii1].element < 0) continue;
				node_index = search_node_label(node, stress->array[ii1].node);
				if(count[node_index] <= 0) return(UNKNOWN_ERROR_RC);
				stress->array[ii1].v
				       = mul_scalar_vect3dr(1.0/(double)(count[node_index]),
				                            stress_sum[node_index]);
			}
		}
		if(strain != NULL){
			for(ii1=0; ii1<strain->size; ii1++){
				if(strain->array[ii1].element < 0) continue;
				node_index = search_node_label(node,strain->array[ii1].node);
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


/* 局所座標系での応力、ひずみ計算 */
RC
local_stress_strain_plate_shell (ELEMENT element, const double *local_point,
                                 NODE_ARRAY node,
                                 MATERIAL_PROP_ARRAY material,
                                 PHYSICAL_PROP_ARRAY physical,
                                 DISP_ARRAY disp,
                                 STRESS_STRAIN *stress,
                                 STRESS_STRAIN *strain,
                                 int sol_ss_face, WORK_SET ws)
{
	int ii1, ii2;
	int dim;
	int ssnum;
	int m_index;
	double disp_v[6*MAX_SHELL_NODE];
	double strain_vect[6];
	double stress_vect[6];
	double global_strain_vect[6];
	double global_stress_vect[6];
	double L[3][3];
	double **node_location = ws.mat_N_3[0];
	double **Bm_matrix = ws.mat_6_3N[0];
	double **Bb_matrix = ws.mat_6_3N[1];
	double **Bs_matrix = ws.mat_6_3N[2];
	double **Dm_matrix = ws.mat_6_6[0];
	double **Db_matrix = ws.mat_6_6[1];
	double **Ds_matrix = ws.mat_6_6[2];
	double **Q_matrix = ws.mat_6_6[3];
	MATERIAL_PROP mat;


	if( (dim = element_dim_shell(element.type)) <= 0
	||(ssnum = element_ssnum_shell(element.type)) <= 0 ){
		return(ARG_ERROR_RC);
	}

	/* Set L[3][3] */
	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( set_L_matrix_shell_solid(element, node_location, L) );

	/* Local displacement */
	RC_TRY( local_node_location(L, element.node_num, node_location) );
	RC_TRY( fill_disp_vect_shell(element, disp, disp_v) );
	RC_TRY( local_disp_vect_shell(L, element.node_num, disp_v) );

	RC_TRY( make_B_matrix_plate_shell(element, local_point, node,
	                                  Bm_matrix, Bb_matrix, Bs_matrix, NULL) );

	m_index = search_material_prop_element(material, physical, &element);
	if(m_index < 0) return(SEEK_ERROR_RC);

	init_vect(6, strain_vect);
	init_vect(6, stress_vect);

	/* 膜 */
	{
		double strain_membrane[3];
		double stress_membrane[3];

		/* Set membrane D matrix */
		mat = material.array[m_index];
		mat.ana_type = ANA_PLANE_STRESS;
		RC_TRY( fill_material(&mat) );
		for(ii1=0;ii1<3;ii1++){
			for(ii2=0;ii2<3;ii2++){
				Dm_matrix[ii1][ii2] = mat.D_matrix[ii1][ii2];
			}
		}

		mul_matrix_vect(3, element.node_num*dim, Bm_matrix, disp_v,
		                strain_membrane); 
		mul_matrix_vect(3, 3, Dm_matrix, strain_membrane, stress_membrane); 

		strain_vect[0] += strain_membrane[0];
		strain_vect[1] += strain_membrane[1];
		strain_vect[5] += strain_membrane[2];
		stress_vect[0] += stress_membrane[0];
		stress_vect[1] += stress_membrane[1];
		stress_vect[5] += stress_membrane[2];
	}

	/* 曲げ */
	{
		double z;
		double thickness;
		double strain_bend[3];
		double stress_bend[3];

		/* Set plate bending D matrix */
		mat = material.array[m_index];
		mat.ana_type = ANA_PLATE_BEND;
		RC_TRY( fill_material(&mat) );
		for(ii1=0;ii1<3;ii1++){
			for(ii2=0;ii2<3;ii2++){
				Db_matrix[ii1][ii2] = mat.D_matrix[ii1][ii2];
			}
		}

		mul_matrix_vect(3, element.node_num*dim, Bb_matrix, disp_v,
		                strain_bend); 

		if( fabs(local_point[2]) > 1.0 ) return(ARG_ERROR_RC);
		RC_TRY( get_thickness(element, physical, &thickness) );
		z = local_point[2] * 0.5 * thickness;
		for(ii1=0; ii1<3; ii1++) strain_bend[ii1] *= z;

		mul_matrix_vect(3, 3, Db_matrix, strain_bend, stress_bend); 

		strain_vect[0] += strain_bend[0];
		strain_vect[1] += strain_bend[1];
		strain_vect[5] += strain_bend[2];
		stress_vect[0] += stress_bend[0];
		stress_vect[1] += stress_bend[1];
		stress_vect[5] += stress_bend[2];
	}

	/* せん断 */
	{
		double strain_shear[2];
		double stress_shear[2];

		/* Set plate bending D matrix for Shearing */
		mat = material.array[m_index];
		mat.ana_type = ANA_PLATE_SHEAR;
		RC_TRY( fill_material(&mat) );
		for(ii1=0;ii1<2;ii1++){
			for(ii2=0;ii2<2;ii2++){
				Ds_matrix[ii1][ii2] = mat.D_matrix[ii1][ii2];
			}
		}

		mul_matrix_vect(2, element.node_num*dim, Bs_matrix, disp_v,
		                strain_shear); 
		mul_matrix_vect(2, 2, Ds_matrix, strain_shear, stress_shear); 

		strain_vect[4] += strain_shear[0];
		strain_vect[3] += strain_shear[1];
		stress_vect[4] += stress_shear[0];
		stress_vect[3] += stress_shear[1];
	}

	/* stress strain in global coordinates */
	RC_TRY( set_trans_matrix4stress_strain(L, Q_matrix) );
	mul_trans_matrix_vect(dim, dim, Q_matrix, strain_vect, global_strain_vect);
	mul_trans_matrix_vect(dim, dim, Q_matrix, stress_vect, global_stress_vect);

	if(strain != NULL){
		init_stress_strain(strain);
		strain->element = element.label;
		RC_TRY( vect2stress_strain(6, global_strain_vect, strain) );
	}
	if(stress != NULL){
		init_stress_strain(stress);
		stress->element = element.label;
		RC_TRY( vect2stress_strain(6, global_stress_vect, stress) );
	}

	return(NORMAL_RC);
}


/* 変換マトリックス(応力，ひずみベクトルを全体座標系へ) */
RC
set_trans_matrix4stress_strain (double L[][3], double **Q_matrix)
{
	double l1, l2, l3;
	double m1, m2, m3;
	double n1, n2, n3;


	l1 = L[0][0]; l2 = L[1][0]; l3 = L[2][0];
	m1 = L[0][1]; m2 = L[1][1]; m3 = L[2][1];
	n1 = L[0][2]; n2 = L[1][2]; n3 = L[2][2];

	Q_matrix[0][0] = l1*l1;
	Q_matrix[0][1] = m1*m1;
	Q_matrix[0][2] = n1*n1;
	Q_matrix[0][3] = m1*n1;
	Q_matrix[0][4] = n1*l1;
	Q_matrix[0][5] = l1*m1;

	Q_matrix[1][0] = l2*l2;
	Q_matrix[1][1] = m2*m2;
	Q_matrix[1][2] = n2*n2;
	Q_matrix[1][3] = m2*n2;
	Q_matrix[1][4] = n2*l2;
	Q_matrix[1][5] = l2*m2;

	Q_matrix[2][0] = l3*l3;
	Q_matrix[2][1] = m3*m3;
	Q_matrix[2][2] = n3*n3;
	Q_matrix[2][3] = m3*n3;
	Q_matrix[2][4] = n3*l3;
	Q_matrix[2][5] = l3*m3;

	Q_matrix[3][0] = 2.0*l2*l3;
	Q_matrix[3][1] = 2.0*m2*m3;
	Q_matrix[3][2] = 2.0*n2*n3;
	Q_matrix[3][3] = m2*n3 + m3*n2;
	Q_matrix[3][4] = n2*l3 + n3*l2;
	Q_matrix[3][5] = l2*m3 + l3*m2;;

	Q_matrix[3][0] = 2.0*l3*l1;
	Q_matrix[3][1] = 2.0*m3*m1;
	Q_matrix[3][2] = 2.0*n3*n1;
	Q_matrix[3][3] = m3*n1 + m1*n3;
	Q_matrix[3][4] = n3*l1 + n1*l3;
	Q_matrix[3][5] = l3*m1 + l1*m3;;

	Q_matrix[5][0] = 2.0*l1*l2;
	Q_matrix[5][1] = 2.0*m1*m2;
	Q_matrix[5][2] = 2.0*n1*n2;
	Q_matrix[5][3] = m1*n2 + m2*n1;
	Q_matrix[5][4] = n1*l2 + n2*l1;
	Q_matrix[5][5] = l1*m2 + l2*m1;;

	return(NORMAL_RC);
}


/* 要素節点変位ベクトル */
RC
fill_disp_vect_shell (ELEMENT element, DISP_ARRAY disp, double disp_vect[])
{
	int ii1;
	int dim;


	if( (dim = element_dim_shell(element.type)) <= 0 ){
		return(ARG_ERROR_RC);
	}

	for(ii1=0; ii1<element.node_num; ii1++){
		int index;

		RC_NEG_CHK( index = search_disp_node(disp, element.node[ii1]) );

		disp_vect[dim*ii1]   = disp.array[index].v.x;
		disp_vect[dim*ii1+1] = disp.array[index].v.y;
		disp_vect[dim*ii1+2] = disp.array[index].v.z;
		disp_vect[dim*ii1+3] = disp.array[index].v.yz;
		disp_vect[dim*ii1+4] = disp.array[index].v.zx;
		disp_vect[dim*ii1+5] = disp.array[index].v.xy;
	}

	return(NORMAL_RC);
}


static RC
local_disp_vect_shell (double L[][3], int node_num, double disp_vect[])
{
	int ii1;
	int dim;

	dim = 6;
	for(ii1=0; ii1<node_num; ii1++){
		mul_matrix33_vect(L, &(disp_vect[ii1*dim]));
		mul_matrix33_vect(L, &(disp_vect[ii1*dim+3]));
	}

	return(NORMAL_RC);
}


/* 虚偽のドリリング自由度 */
static RC
set_drilling_dof (int mitc_node_num, ELEMENT *element,
                  NODE_ARRAY node, MATERIAL_PROP_ARRAY material,
                  PHYSICAL_PROP_ARRAY physical, double thickness,
                  ELEM_MATRIX *elem_matrix)
{
	int ii1, ii2;
	int dim = 6;
	int m_index;
	double factor;
	double alpha = 0.03;
	double diagonal_val;
	double off_diagonal_val;
	MATERIAL_PROP tmp_mat;


	m_index = search_material_prop_element(material, physical, element);
	if(m_index < 0) return(SEEK_ERROR_RC);

	tmp_mat = material.array[m_index];
	if(tmp_mat.mat_type != MAT_ISOTROPIC) return(ARG_ERROR_RC);

	RC_TRY( element_volume(element, node) );

	factor = alpha * tmp_mat.E * element->volume * thickness;
	diagonal_val = factor;
	off_diagonal_val = -(1.0 / (double)(mitc_node_num - 1)) * factor;
	for(ii1=0; ii1<mitc_node_num; ii1++){
		for(ii2=(ii1+1); ii2<mitc_node_num; ii2++){
			elem_matrix->matrix[dim*ii1+5][dim*ii2+5] =
				elem_matrix->matrix[dim*ii2+5][dim*ii1+5] = off_diagonal_val;
		}
		elem_matrix->matrix[dim*ii1+5][dim*ii1+5] = diagonal_val;
	}

	return(NORMAL_RC);
}


RC
make_plate_shell_mass_matrix (ELEMENT element, NODE_ARRAY node,
                              MATERIAL_PROP_ARRAY material,
                              PHYSICAL_PROP_ARRAY physical,
                              ELEM_MATRIX *elem_matrix,
                              WORK_SET ws)
{
	int ii1, ii2, ii3;
	int dim;
	int ssnum;
	int g_num_point;
	int prop_index;
	int m_index;
	double thickness;
	double rho;
	double det_J;
	double *N = ws.vec_N[0];
	double *g_weights = ws.vec_G[0];
	double **g_points = ws.mat_G_4[0];
	double **node_location = ws.mat_N_3[0];
	double **Jacobi = ws.mat_3_3[0];
	double **dN = ws.mat_3_N[0];
	double **N_matrix = ws.mat_6_3N[0];
	double L[3][3];


	if(element.label < 0) return(ARG_ERROR_RC);
	prop_index = search_physical_prop_label(physical, element.physical);
	if(prop_index < 0) return(SEEK_ERROR_RC);
	if(NST_PROP_PSHELL != physical.array[prop_index].i_info[1]){
		return(ARG_ERROR_RC);
	}

	if( (dim = element_dim_shell(element.type)) <= 0
	|| (ssnum = element_ssnum_shell(element.type)) <= 0 ){
		return(ARG_ERROR_RC);
	}

	elem_matrix->dim = dim;
	elem_matrix->element = element.label;
	elem_matrix->size = (element.node_num) * dim;

	RC_TRY( allocate1I(elem_matrix->size, &(elem_matrix->index)) );
	for(ii1=0; ii1<element.node_num; ii1++){
		int index;

		index = node_renum_index(node, element.node[ii1]);
		if(index < 0 ) return(SEEK_ERROR_RC);

		for(ii2=0; ii2<dim; ii2++){
			elem_matrix->index[ii1*dim+ii2] = index*dim + ii2;
		}
	}

	RC_TRY( allocate2D(elem_matrix->size, elem_matrix->size,
	                   &(elem_matrix->matrix)) );

	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( set_L_matrix_shell_solid(element, node_location, L) );
	RC_TRY( local_node_location(L, element.node_num, node_location) );

	RC_TRY( get_thickness(element, physical, &thickness) );

	m_index = search_material_prop_element(material, physical, &element);
	if(m_index < 0) return(SEEK_ERROR_RC);
	rho = material.array[m_index].rho;

	/* membrane */
	{
		double **NtN = ws.mat_3N_3N[0];

		RC_TRY( Gauss_const_default(element.type, g_points, g_weights,
		                            &g_num_point) );

		for(ii1=0; ii1<g_num_point; ii1++){
			double factor;

			RC_TRY( set_N_vector(element.type, g_points[ii1], N) );
			RC_TRY( set_dN_matrix(element.type, g_points[ii1], dN) );
			mul_matrix(2, element.node_num, 2, dN, node_location, Jacobi);

			det_J = determinant(2, Jacobi);
			if(det_J < 0.0) return(CAL_ERROR_RC);

			RC_TRY( set_N_matrix_membrane(element.node_num, N, N_matrix) );
			mul_matrix_AtB(2, element.node_num*dim, element.node_num*dim,
			               N_matrix, N_matrix, NtN);

			factor = rho * g_weights[ii1] * det_J * thickness;
			for(ii2=0; ii2<elem_matrix->size; ii2++){
				for(ii3=0; ii3<elem_matrix->size; ii3++){
					elem_matrix->matrix[ii2][ii3] += NtN[ii2][ii3] * factor;
				}
			}
		}
	}

	/* bending */
	{
		double **I_matrix = ws.mat_6_6[0];
		double **NtIN = ws.mat_3N_3N[0];

		/* I matrix */
		init_matrix(3, 3, I_matrix);
		I_matrix[0][0] = rho * thickness;
		I_matrix[1][1] = I_matrix[2][2]
		               = rho * (thickness*thickness*thickness/12.0);

		RC_TRY( Gauss_const_default(element.type, g_points, g_weights,
		                            &g_num_point) );

		for(ii1=0; ii1<g_num_point; ii1++){
			double factor;

			RC_TRY( set_N_vector(element.type, g_points[ii1], N) );
			RC_TRY( set_dN_matrix(element.type, g_points[ii1], dN) );
			mul_matrix(2, element.node_num, 2, dN, node_location, Jacobi);

			det_J = determinant(2, Jacobi);
			if(det_J < 0.0) return(CAL_ERROR_RC);

			RC_TRY( set_N_matrix_bend(element.node_num, N, N_matrix) );
			mul_matrix_AtBA(3, element.node_num*dim, N_matrix, I_matrix, NtIN);

			factor = g_weights[ii1] * det_J;
			for(ii2=0; ii2<elem_matrix->size; ii2++){
				for(ii3=0; ii3<elem_matrix->size; ii3++){
					elem_matrix->matrix[ii2][ii3] += NtIN[ii2][ii3] * factor;
				}
			}
		}
	}

	RC_TRY( mul_matrix33_LALt(elem_matrix->size, elem_matrix->matrix, L) );

	return(NORMAL_RC);
}


static RC
set_N_matrix_membrane (int node_num, double *N, double **N_matrix)
{
	int ii1;
	int dim;


	dim = 6;
	for(ii1=0; ii1<node_num; ii1++){
		N_matrix[0][dim*ii1+0] = N[ii1];
		N_matrix[0][dim*ii1+1] = 0.0;
		N_matrix[0][dim*ii1+2] = 0.0;
		N_matrix[0][dim*ii1+3] = 0.0;
		N_matrix[0][dim*ii1+4] = 0.0;
		N_matrix[0][dim*ii1+5] = 0.0;

		N_matrix[1][dim*ii1+0] = 0.0;
		N_matrix[1][dim*ii1+1] = N[ii1];
		N_matrix[1][dim*ii1+2] = 0.0;
		N_matrix[1][dim*ii1+3] = 0.0;
		N_matrix[1][dim*ii1+4] = 0.0;
		N_matrix[1][dim*ii1+5] = 0.0;
	}
	
	return(NORMAL_RC);
}


static RC
set_N_matrix_bend (int node_num, double *N, double **N_matrix)
{
	int ii1;
	int dim;


	dim = 6;
	for(ii1=0; ii1<node_num; ii1++){
		N_matrix[0][dim*ii1+0] = 0.0;
		N_matrix[0][dim*ii1+1] = 0.0;
		N_matrix[0][dim*ii1+2] = N[ii1];
		N_matrix[0][dim*ii1+3] = 0.0;
		N_matrix[0][dim*ii1+4] = 0.0;
		N_matrix[0][dim*ii1+5] = 0.0;

		N_matrix[1][dim*ii1+0] = 0.0;
		N_matrix[1][dim*ii1+1] = 0.0;
		N_matrix[1][dim*ii1+2] = 0.0;
		N_matrix[1][dim*ii1+3] = 0.0;
		N_matrix[1][dim*ii1+4] = -N[ii1];
		N_matrix[1][dim*ii1+5] = 0.0;

		N_matrix[2][dim*ii1+0] = 0.0;
		N_matrix[2][dim*ii1+1] = 0.0;
		N_matrix[2][dim*ii1+2] = 0.0;
		N_matrix[2][dim*ii1+3] = N[ii1];
		N_matrix[2][dim*ii1+4] = 0.0;
		N_matrix[2][dim*ii1+5] = 0.0;
	}
	
	return(NORMAL_RC);
}



