/*********************************************************************
 * fem_solver_shell.c
 *
 * Copyright (C) 2003 AzLib Developers Group
 *
 * Written by
 *  <Yasuo ASAGA> <Taiki AOYAMA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: fem_solver_shell.c 1021 2012-04-30 04:09:18Z aoyama $*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "base.h"
#include "mathematics.h"
#include "fem_solver.h"
#include "fem_struct.h"
#include "sky_cholesky.h"
#include "nst_component.h"

#define ORTHOGONAL_ERROR (1.0e-6)


static RC set_drilling_dof(int mitc_node_num, ELEMENT *element,
                           NODE_ARRAY node, MATERIAL_PROP_ARRAY material,
                           PHYSICAL_PROP_ARRAY physical, double thickness,
                           ELEM_MATRIX *elem_matrix);
static RC trans_global_rotation_dof(int matrix_size, double **matrix,
                                    DIRECTOR_VECT elem_vect[]);
static RC set_center_node_director_vect(DIRECTOR_VECT elem_vect[],
                                        WORK_SET ws);
static RC set_center_node_location(double **node_location, WORK_SET ws);
static RC reset_shell_element_type(ELEM_TYPE input_type, ELEM_TYPE *ouput_type,
                                   int *output_node_num);
static RC set_unit_normal_vect_shell(double *local_point, ELEMENT element,
                                     NODE_ARRAY node, VECT3D *vect,
                                     WORK_SET ws);
static RC set_Q_matrix4stress(VECT3D G[], double **Q_matrix);
static RC set_Q_matrix4strain(VECT3D g[], double **Q_matrix);
static RC fill_disp5dof_vect_shell(ELEMENT element, ELEM_TYPE elem_type,
                                   int elem_node_num, DIRECTOR_VECT elem_vect[],
                                   DISP_ARRAY disp, double disp5dof_v[]);
static RC fill_disp6dof_vect_shell(ELEMENT element, ELEM_TYPE elem_type,
                                   int elem_node_num, DIRECTOR_VECT elem_vect[],
                                   DISP_ARRAY disp, double disp5dof_v[]);
static RC make_B_matrix_mitc3(ELEM_TYPE elem_type, int elem_node_num,
                              DIRECTOR_VECT elem_vect[], double thickness,
                              double **node_location, const double *local_point,
                              double **B_matrix,
                              MITC_TYING_SET mitc, WORK_SET ws);
static RC make_B_matrix_mitc6(ELEM_TYPE elem_type, int elem_node_num,
                              DIRECTOR_VECT elem_vect[], double thickness,
                              double **node_location, const double *local_point,
                              double **B_matrix,
                              MITC_TYING_SET mitc, WORK_SET ws);
static RC make_B_matrix_mitc4(ELEM_TYPE elem_type, int elem_node_num,
                              DIRECTOR_VECT elem_vect[], double thickness,
                              double **node_location,
                              const double *local_point, double **B_matrix,
                              MITC_TYING_SET mitc, WORK_SET ws);
static RC make_B_matrix_mitc8(ELEM_TYPE elem_type, int elem_node_num,
                              DIRECTOR_VECT elem_vect[], double thickness,
                              double **node_location,
                              const double *local_point, double **B_matrix,
                              MITC_TYING_SET mitc, WORK_SET ws);
static RC make_B_matrix_mitc9(ELEM_TYPE elem_type, int elem_node_num,
                              DIRECTOR_VECT elem_vect[], double thickness,
                              double **node_location, const double *local_point,
                              double **B_matrix,
                              MITC_TYING_SET mitc, WORK_SET ws);
static RC set_disp_N_matrix_shell (ELEM_TYPE elem_type, int elem_node_num,
                                   DIRECTOR_VECT elem_vect[], double thickness,
                                   const double *local_point,
                                   double **N_matrix, WORK_SET ws);
static RC Gauss_const_mitc_tri(int num_point, double **rst, double *weight);
static RC Gauss_const_mitc_quad(int num_point, double **rst, double *weight);
static RC N_mitc3(const double *rst, double *N);
static RC N_mitc4(const double *rst, double *N);
static RC N_mitc6(const double *rst, double *N);
static RC N_mitc8(const double *rst, double *N);
static RC N_mitc9(const double *rst, double *N);
static RC dN_mitc3(const double *rst, double **dN);
static RC dN_mitc4(const double *rst, double **dN);
static RC dN_mitc6(const double *rst, double **dN);
static RC dN_mitc8(const double *rst, double **dN);
static RC dN_mitc9(const double *rst, double **dN);


/*
 * シェル要素用微小変形解析ソルバー
 * ロッキング回避のためにMITCシェル要素使用
 *
 * 連立方程式の求解にはICCG使用するが，板厚が薄くなると正定値性が弱くなるので
 * fem_solve_shellの使用を推奨する．
 */
RC fem_solve_shell_iccg6 (NODE_ARRAY node, ELEMENT_ARRAY element,
                          RIGID_ELEMENT_ARRAY rigid,
                          MATERIAL_PROP_ARRAY material,
                          PHYSICAL_PROP_ARRAY physical,
                          BC_ARRAY rest, BC_ARRAY force,
                          BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                          DEFAULT_BC_ARRAY def_body_force,
                          DISP_ARRAY *disp,
                          STRESS_ARRAY *stress, STRAIN_ARRAY *strain,
                          int sol_ss_flag, int sol_ss_face)
{
	int ii1;
	int dim;
	double *f_vect;
	double *u_vect;
	NONZERO_MATRIX3 matrix;
	DIRECTOR_VECT_ARRAY vect;
	LOCAL_COORD_ARRAY dummy_coord;
	WORK_SET ws1, ws2;
	MITC_TYING_SET mitc;


	if((dim = analysis_dim_shell(element)) != 6){
		return(ARG_ERROR_RC);
	}
	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

	RC_TRY( make_nodal_director_vect(node, element, &vect) );

	RC_TRY( allocate_nonzero6(node, element, rigid, &matrix) );

	RC_TRY( allocate_work_set(&ws1) );
	RC_TRY( allocate_work_set(&ws2) );
	RC_TRY( allocate_mitc_tying_set(&mitc) );
	for(ii1=0; ii1<element.size; ii1++){
		ELEM_MATRIX elem_matrix;

		if(element.array[ii1].label < 0) continue;

		RC_TRY( make_shell_elem_matrix(element.array[ii1], node,
		                               material, physical, vect,
		                               &elem_matrix, mitc, ws1, ws2) );
		RC_TRY( impose_nonzero6(matrix, elem_matrix, 1.0) );

		RC_TRY( free_elem_matrix(&elem_matrix) );
	}
	RC_TRY( free_work_set(&ws1) );
	RC_TRY( free_work_set(&ws2) );
	RC_TRY( free_mitc_tying_set(&mitc) );

	RC_TRY( fem_allocate_vector(dim, node, &f_vect) );
	/* 節点荷重・モーメント */
	RC_TRY( add_nodal_force(dim, node, force, f_vect) );
	/* 熱荷重 */
	RC_TRY( add_temp_force_shell(node, element, material, physical,
	                             temp, def_temp, vect, f_vect) );
	/* 重力荷重 */
	RC_TRY( add_body_force_shell(node, element, material, physical,
	                             def_body_force, vect, f_vect) );

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
	RC_TRY( stress_strain_shell(sol_ss_flag, sol_ss_face, node, element,
	                            material, physical, vect, *disp,
	                            stress, strain) );
	RC_TRY( free_director_vect_array(&vect) );

	return(NORMAL_RC);
}


/*
 * シェル要素用微小変形解析ソルバー
 * ロッキング回避のためにMITCシェル要素使用
 *
 * 連立方程式の求解にはスカイライン法使用
 */
RC fem_solve_shell (NODE_ARRAY node, ELEMENT_ARRAY element,
                    RIGID_ELEMENT_ARRAY rigid,
                    MATERIAL_PROP_ARRAY material,
                    PHYSICAL_PROP_ARRAY physical,
                    BC_ARRAY rest, BC_ARRAY force,
                    BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                    DEFAULT_BC_ARRAY def_body_force,
                    DISP_ARRAY *disp,
                    STRESS_ARRAY *stress, STRAIN_ARRAY *strain,
                    int sol_ss_flag, int sol_ss_face)
{
	int ii1;
	int dim;
	double *f_vect;
	SKYLINE_MATRIX matrix;
	DIRECTOR_VECT_ARRAY vect;
	WORK_SET ws1, ws2;
	MITC_TYING_SET mitc;


	if((dim = analysis_dim_shell(element)) != 6){
		return(ARG_ERROR_RC);
	}
	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

	RC_TRY( make_nodal_director_vect(node, element, &vect) );

	RC_TRY( allocate_g_matrix(dim, node, element, &matrix) );

	RC_TRY( allocate_work_set(&ws1) );
	RC_TRY( allocate_work_set(&ws2) );
	RC_TRY( allocate_mitc_tying_set(&mitc) );
	for(ii1=0; ii1<element.size; ii1++){
		ELEM_MATRIX elem_matrix;

		if(element.array[ii1].label < 0) continue;

		RC_TRY( make_shell_elem_matrix(element.array[ii1], node,
		                               material, physical, vect,
		                               &elem_matrix, mitc, ws1, ws2) );
		RC_TRY( impose_elem_matrix(matrix, elem_matrix, 1.0) );

		RC_TRY( free_elem_matrix(&elem_matrix) );
	}
	RC_TRY( free_work_set(&ws1) );
	RC_TRY( free_work_set(&ws2) );
	RC_TRY( free_mitc_tying_set(&mitc) );

	RC_TRY( fem_allocate_vector(dim, node, &f_vect) );
	/* 節点荷重・モーメント */
	RC_TRY( add_nodal_force(dim, node, force, f_vect) );
	/* 熱荷重 */
	RC_TRY( add_temp_force_shell(node, element, material, physical,
	                             temp, def_temp, vect, f_vect) );
	/* 重力荷重 */
	RC_TRY( add_body_force_shell(node, element, material, physical,
	                             def_body_force, vect, f_vect) );

	RC_TRY( modify_g_matrix6(node, rest, matrix, f_vect) );
	RC_TRY( decomp_g_matrix(matrix) );
	RC_TRY( s_cholesky_subst(matrix.dof, matrix.index1, matrix.index2,
	                         matrix.array, f_vect, 1) );
	RC_TRY( vect2disp(dim, node, f_vect, disp) );
	RC_TRY( free_g_matrix(&matrix) );
	RC_TRY( fem_free_vector(&f_vect) );
	RC_TRY( stress_strain_shell(sol_ss_flag, sol_ss_face, node, element,
	                            material, physical, vect, *disp,
	                            stress, strain) );
	RC_TRY( free_director_vect_array(&vect) );

	return(NORMAL_RC);
}


RC
modify_g_matrix6 (NODE_ARRAY node, BC_ARRAY rest, SKYLINE_MATRIX g_matrix,
                  double f_vect[])
{
	int ii1;
	int index;
	double *ptr;
	RC rc;


	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);
	if(g_matrix.dim != 6) return(ARG_ERROR_RC);

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
		if(rest.array[ii1].v_type.z == BC_FIX){
			RC_TRY( s_cholesky_mod(g_matrix.dof, g_matrix.index1,
			                       g_matrix.index2, g_matrix.array, f_vect,
			                       g_matrix.dim*index + 2,
			                       rest.array[ii1].v.z) );
		}
		if(rest.array[ii1].v_type.yz == BC_FIX){
			RC_TRY( s_cholesky_mod(g_matrix.dof, g_matrix.index1,
			                       g_matrix.index2, g_matrix.array, f_vect,
			                       g_matrix.dim*index + 3,
			                       rest.array[ii1].v.yz) );
		}
		if(rest.array[ii1].v_type.zx == BC_FIX){
			RC_TRY( s_cholesky_mod(g_matrix.dof, g_matrix.index1,
			                       g_matrix.index2, g_matrix.array, f_vect,
			                       g_matrix.dim*index + 4,
			                       rest.array[ii1].v.zx) );
		}
		if(rest.array[ii1].v_type.xy == BC_FIX){
			RC_TRY( s_cholesky_mod(g_matrix.dof, g_matrix.index1,
			                       g_matrix.index2, g_matrix.array, f_vect,
			                       g_matrix.dim*index + 5,
			                       rest.array[ii1].v.xy) );
		}
	}

	/* 代入されてない行(列)の対角要素を 1.0 に */
	for(ii1=0; ii1<g_matrix.dof; ii1++){
		if(g_matrix.index1[ii1] == ii1){
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


RC
modify_g_matrix6_n (NODE_ARRAY node, BC_ARRAY rest,
                    SKYLINE_MATRIX g_matrix, double **f_vect, int f_vect_size)
{
	int ii1;
	int index;


	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);
	if(g_matrix.dim != 6) return(ARG_ERROR_RC);

	for(ii1=0; ii1<rest.size; ii1++){
		if(rest.array[ii1].node < 0) continue;

		index = node_renum_index(node, rest.array[ii1].node);
		if(index < 0) return(SEEK_ERROR_RC);

		if(rest.array[ii1].v_type.x == BC_FIX){
			RC_TRY( s_cholesky_mod_n(g_matrix.dof, g_matrix.index1,
			                         g_matrix.index2, g_matrix.array, f_vect,
			                         g_matrix.dim*index,
			                         rest.array[ii1].v.x, f_vect_size) );
		}
		if(rest.array[ii1].v_type.y == BC_FIX){
			RC_TRY( s_cholesky_mod_n(g_matrix.dof, g_matrix.index1,
			                         g_matrix.index2, g_matrix.array, f_vect,
			                         g_matrix.dim*index + 1,
			                         rest.array[ii1].v.y, f_vect_size) );
		}
		if(rest.array[ii1].v_type.z == BC_FIX){
			RC_TRY( s_cholesky_mod_n(g_matrix.dof, g_matrix.index1,
			                         g_matrix.index2, g_matrix.array, f_vect,
			                         g_matrix.dim*index + 2,
			                         rest.array[ii1].v.z, f_vect_size) );
		}
		if(rest.array[ii1].v_type.yz == BC_FIX){
			RC_TRY( s_cholesky_mod_n(g_matrix.dof, g_matrix.index1,
			                         g_matrix.index2, g_matrix.array, f_vect,
			                         g_matrix.dim*index + 3,
			                         rest.array[ii1].v.yz, f_vect_size) );
		}
		if(rest.array[ii1].v_type.zx == BC_FIX){
			RC_TRY( s_cholesky_mod_n(g_matrix.dof, g_matrix.index1,
			                         g_matrix.index2, g_matrix.array, f_vect,
			                         g_matrix.dim*index + 4,
			                         rest.array[ii1].v.zx, f_vect_size) );
		}
		if(rest.array[ii1].v_type.xy == BC_FIX){
			RC_TRY( s_cholesky_mod_n(g_matrix.dof, g_matrix.index1,
			                         g_matrix.index2, g_matrix.array, f_vect,
			                         g_matrix.dim*index + 5,
			                         rest.array[ii1].v.xy, f_vect_size) );
		}
	}

	/* 代入されてない行(列)の対角要素を 1.0 に */
	for(ii1=0; ii1<g_matrix.dof; ii1++){
		if(g_matrix.index1[ii1] == ii1){
			RC rc;
			double *ptr = s_cholesky_ptr(g_matrix.dof, g_matrix.index1,
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


/* 要素剛性マトリクス */
RC
make_shell_elem_matrix (ELEMENT element, NODE_ARRAY node,
                        MATERIAL_PROP_ARRAY material,
                        PHYSICAL_PROP_ARRAY physical,
                        DIRECTOR_VECT_ARRAY d_vect,
                        ELEM_MATRIX *elem_matrix,
                        MITC_TYING_SET mitc,
                        WORK_SET ws1, WORK_SET ws2)
{
	int ii1, ii2, ii3, ii4, ii5;
	int dim;
	int ssnum;
	int elem_node_num;
	int elem_matrix_size;
	int g_num_point;
	int prop_index;
	double thickness;
	double det_J;
	double *g_weights = ws1.vec_G[0];
	double **node_location = ws1.mat_N_3[0];
	double **g_points = ws1.mat_G_4[0];
	double **D_matrix = ws1.mat_6_6[0];
	double **B_matrix = ws1.mat_6_3N[0];
	double **BtDB = ws1.mat_3N_3N[0];
	VECT3D G0[3], g0[3];
	DIRECTOR_VECT elem_vect[MAX_SHELL_NODE];
	ELEM_TYPE elem_type;


	if(element.label < 0) return(ARG_ERROR_RC);

	prop_index = search_physical_prop_label(physical, element.physical);
	if(prop_index < 0) return(SEEK_ERROR_RC);
	if(NST_PROP_PSHELL != physical.array[prop_index].i_info[1]){
		return(ARG_ERROR_RC);
	}

	if( (dim = element_dim_shell(element.type)) != 6
	|| (ssnum = element_ssnum_shell(element.type)) != 5 ){
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

	/* 要素タイプの再設定 ELEM_QUAD2 -> ELEM_QUAD2L */
	RC_TRY( reset_shell_element_type(element.type, &elem_type,
	                                 &elem_node_num) );

	elem_matrix_size = elem_node_num * dim;
	RC_TRY( allocate2D(elem_matrix_size, elem_matrix_size,
	                   &(elem_matrix->matrix)) );

	RC_TRY( set_element_direcor_vect(element, d_vect, elem_vect) );
	RC_TRY( node_location_matrix(&element, node, node_location) );
	if(elem_type == ELEM_QUAD2L){
		RC_TRY( set_center_node_location(node_location, ws2) );
		RC_TRY( set_center_node_director_vect(elem_vect, ws2) );
	}

	RC_TRY( Gauss_const_default_shell(elem_type, g_points, g_weights,
	                                  &g_num_point) );
	RC_TRY( get_thickness(element, physical, &thickness) );

	for(ii1=0; ii1<g_num_point; ii1++){
		RC_TRY( make_G1G2G3g1g2g3_basis(elem_type, elem_node_num,
		                                elem_vect, thickness,
		                                node_location, g_points[ii1],
		                                G0, g0, ws2) );

		RC_TRY( make_B_matrix_shell(elem_type, elem_node_num,
		                            elem_vect, thickness, node_location,
		                            g_points[ii1], B_matrix, mitc, ws2) );

		RC_TRY( make_D_matrix_shell(element, material, physical,
		                            g0, D_matrix) );

		mul_matrix_AtBA(ssnum, 5*elem_node_num, B_matrix, D_matrix, BtDB);

		det_J = scalar_triple_product3d(G0[0], G0[1], G0[2]);

		for(ii2=0; ii2<elem_node_num; ii2++){
			for(ii3=0; ii3<elem_node_num; ii3++){
				for(ii4=0; ii4<5; ii4++){
					for(ii5=0; ii5<5; ii5++){
						elem_matrix->matrix[dim*ii2+ii4][dim*ii3+ii5]
						+= BtDB[5*ii2+ii4][5*ii3+ii5]*det_J*g_weights[ii1];
					}
				}
			}
		}
	}

	/* 虚偽のドリリング剛性 */
	RC_TRY( set_drilling_dof(elem_node_num, &element, node, material, physical,
	                         thickness, elem_matrix) );
	/* 全体座標系成分への変換 */
	RC_TRY( trans_global_rotation_dof(elem_matrix_size, elem_matrix->matrix,
	                                  elem_vect) );

	/* 余分な自由度を静縮小(Guyan reduction) */
	if(elem_type == ELEM_QUAD2L){
		int elem_size = elem_matrix->size;
		int realloc_size = elem_matrix_size*elem_matrix_size
		                 - dim*dim*elem_node_num;

		for(ii2=(dim-1); ii2>=0; ii2--){
			RC_TRY( static_condens_matrix(dim*(elem_node_num-1) + ii2 + 1,
			                              elem_matrix->matrix, NULL,
			                              dim*(elem_node_num-1) + ii2) );
		}
		elem_matrix->matrix[0] = (double *)mm_realloc(elem_matrix->matrix[0],
		                                   realloc_size*sizeof(double));
		elem_matrix->matrix = (double **)mm_realloc(elem_matrix->matrix,
		                                 elem_size*sizeof(double *));
	}

	return(NORMAL_RC);
}


/* 虚偽のドリリング剛性 */
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


/* 要素剛性マトリクスの回転成分を全体座標系成分へ変換 */
static RC
trans_global_rotation_dof (int matrix_size, double **matrix,
                           DIRECTOR_VECT elem_vect[])
{
	int ii1, ii2, ii3, ii4;
	int node_size;
	int n;
	double tmp[3][3];
	double V_matrix[3][3];

	if(matrix_size % 6 != 0) return(ARG_ERROR_RC);

	node_size = matrix_size / 6;
	n = matrix_size / 3;

	for(ii1=0; ii1<node_size; ii1++){
		int ii1_6p3 = ii1*6 + 3;
		VECT3D V1 = elem_vect[ii1].v1;
		VECT3D V2 = elem_vect[ii1].v2;
		VECT3D V3 = elem_vect[ii1].vn;

		V_matrix[0][0] = V1.x;
		V_matrix[0][1] = V1.y;
		V_matrix[0][2] = V1.z;
		V_matrix[1][0] = V2.x;
		V_matrix[1][1] = V2.y;
		V_matrix[1][2] = V2.z;
		V_matrix[2][0] = V3.x;
		V_matrix[2][1] = V3.y;
		V_matrix[2][2] = V3.z;

		for(ii2=0; ii2<n; ii2++){
			int ii2_3 = ii2 * 3;

			for(ii3=0; ii3<3; ii3++){
				tmp[ii3][0] = matrix[ii2_3+ii3][ii1_6p3]  *V_matrix[0][0]
				            + matrix[ii2_3+ii3][ii1_6p3+1]*V_matrix[1][0]
				            + matrix[ii2_3+ii3][ii1_6p3+2]*V_matrix[2][0];
				tmp[ii3][1] = matrix[ii2_3+ii3][ii1_6p3]  *V_matrix[0][1]
				            + matrix[ii2_3+ii3][ii1_6p3+1]*V_matrix[1][1]
				            + matrix[ii2_3+ii3][ii1_6p3+2]*V_matrix[2][1];
				tmp[ii3][2] = matrix[ii2_3+ii3][ii1_6p3]  *V_matrix[0][2]
				            + matrix[ii2_3+ii3][ii1_6p3+1]*V_matrix[1][2]
				            + matrix[ii2_3+ii3][ii1_6p3+2]*V_matrix[2][2];
			}
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					matrix[ii2_3+ii3][ii1_6p3+ii4] = tmp[ii3][ii4];
				}
			}
		}
		for(ii2=0; ii2<n; ii2++){
			int ii2_3 = ii2 * 3;

			for(ii3=0; ii3<3; ii3++){
				tmp[ii3][0] = V_matrix[0][ii3]*matrix[ii1_6p3][ii2_3]
				            + V_matrix[1][ii3]*matrix[ii1_6p3+1][ii2_3]
				            + V_matrix[2][ii3]*matrix[ii1_6p3+2][ii2_3];
				tmp[ii3][1] = V_matrix[0][ii3]*matrix[ii1_6p3][ii2_3+1]
				            + V_matrix[1][ii3]*matrix[ii1_6p3+1][ii2_3+1]
				            + V_matrix[2][ii3]*matrix[ii1_6p3+2][ii2_3+1];
				tmp[ii3][2] = V_matrix[0][ii3]*matrix[ii1_6p3][ii2_3+2]
				            + V_matrix[1][ii3]*matrix[ii1_6p3+1][ii2_3+2]
				            + V_matrix[2][ii3]*matrix[ii1_6p3+2][ii2_3+2];
			}
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					matrix[ii1_6p3+ii3][ii2_3+ii4] = tmp[ii3][ii4];
				}
			}
		}
	}

	return(NORMAL_RC);
}


/* 要素中心でのディレクターベクトル(内挿による導出) */
static RC
set_center_node_director_vect (DIRECTOR_VECT elem_vect[], WORK_SET ws)
{
	int ii1;
	double abs_val;
	double *center_point = ws.vec_N[0];
	double *N = ws.vec_N[1];
	VECT3D ey = {0.0, 1.0, 0.0};
	VECT3D ez = {0.0, 0.0, 1.0};
	VECT3D tmp_v;


	center_point[0] = 0.0;
	center_point[1] = 0.0;
	center_point[2] = 0.0;

	RC_TRY( set_N_vector(ELEM_QUAD2, center_point, N) );

	RC_TRY( init_vect3d(&(elem_vect[8].vn)) );
	RC_TRY( init_vect3d(&(elem_vect[8].v1)) );
	RC_TRY( init_vect3d(&(elem_vect[8].v2)) );
	for(ii1=0; ii1<8; ii1++){
		elem_vect[8].vn.x += N[ii1]*elem_vect[ii1].vn.x;
		elem_vect[8].vn.y += N[ii1]*elem_vect[ii1].vn.y;
		elem_vect[8].vn.z += N[ii1]*elem_vect[ii1].vn.z;
	}

	tmp_v = outer_product3d(ey, elem_vect[8].vn);
	abs_val = abs_vect3d(tmp_v);
	if(abs_val < ORTHOGONAL_ERROR){
		elem_vect[8].v1 = ez;
	}else{
		elem_vect[8].v1.x = tmp_v.x / abs_val;
		elem_vect[8].v1.y = tmp_v.y / abs_val;
		elem_vect[8].v1.z = tmp_v.z / abs_val;
	}
	elem_vect[8].v2 = unit_vect3d(outer_product3d(elem_vect[8].vn,
	                              elem_vect[8].v1) );

	return(NORMAL_RC);
}


/* 要素中心での節点座標マトリクス成分(内挿による導出) */
static RC
set_center_node_location (double **node_location, WORK_SET ws)
{
	int ii1;
	double *center_point = ws.vec_N[0];
	double *N = ws.vec_N[1];


	center_point[0] = 0.0;
	center_point[1] = 0.0;
	center_point[2] = 0.0;

	RC_TRY( set_N_vector(ELEM_QUAD2, center_point, N) );

	node_location[8][0] = 0.0;
	node_location[8][1] = 0.0;
	node_location[8][2] = 0.0;
	for(ii1=0; ii1<8; ii1++){
		node_location[8][0] += N[ii1] * node_location[ii1][0];
		node_location[8][1] += N[ii1] * node_location[ii1][1];
		node_location[8][2] += N[ii1] * node_location[ii1][2];
	}

	return(NORMAL_RC);
}


/* 要素タイプの再設定 ELEM_QUAD2 -> ELEM_QUAD2L */
static RC
reset_shell_element_type (ELEM_TYPE input_type, ELEM_TYPE *ouput_type,
                          int *output_node_num)
{
	switch(input_type){
	case ELEM_TRI1:
		*ouput_type = ELEM_TRI1;
		*output_node_num = 3;
		break;
	case ELEM_QUAD1:
		*ouput_type = ELEM_QUAD1;
		*output_node_num = 4;
		break;
	case ELEM_TRI2:
		*ouput_type = ELEM_TRI2;
		*output_node_num = 6;
		break;
	case ELEM_QUAD2:
		*ouput_type = ELEM_QUAD2L;
		*output_node_num = 9;
		break;
	case ELEM_QUAD2L:
		*ouput_type = ELEM_QUAD2L;
		*output_node_num = 9;
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


RC
add_body_force_shell (NODE_ARRAY node, ELEMENT_ARRAY element,
                      MATERIAL_PROP_ARRAY material,
                      PHYSICAL_PROP_ARRAY physical,
                      DEFAULT_BC_ARRAY def_b_force,
                      DIRECTOR_VECT_ARRAY d_vect, double f_vect[])
{
	int ii1, ii2, ii3;
	int dim;
	int mat_index;
	int elem_node_num;
	int g_num_point;
	double rho;
	double det_J;
	double thickness;
	double grav_vect[3];
	double Nt_grav[5*MAX_SHELL_NODE];
	double elem_b_force[5*MAX_SHELL_NODE];
	double *g_weights = NULL;
	double **g_points = NULL;
	double **node_location = NULL;
	double **N_matrix = NULL;
	VECT3D accel_grav;
	VECT3D G0[3];
	DIRECTOR_VECT elem_vect[MAX_SHELL_NODE];
	ELEM_TYPE elem_type;
	WORK_SET ws1, ws2;


	if((dim = analysis_dim_shell(element)) != 6){
		return(ARG_ERROR_RC);
	}

	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

	if(def_b_force.size <= 0) return(NORMAL_RC);

	RC_TRY( allocate_work_set(&ws1) );
	RC_TRY( allocate_work_set(&ws2) );

	g_weights = ws1.vec_G[0];
	node_location = ws1.mat_N_3[0];
	g_points = ws1.mat_G_4[0];
	N_matrix = ws1.mat_6_3N[0];

	/* acceleration of gravity */
	accel_grav.x = accel_grav.y = accel_grav.z = 0.0;
	for(ii1=0; ii1<def_b_force.size; ii1++){
		if(def_b_force.array[ii1].set_id < 0) continue;

		accel_grav.x += def_b_force.array[ii1].grav.x;
		accel_grav.y += def_b_force.array[ii1].grav.y;
		accel_grav.z += def_b_force.array[ii1].grav.z;
	}

	for(ii1=0; ii1<element.size; ii1++){
		ELEMENT elem = element.array[ii1];

		if(elem.label < 0) continue;

		elem_type = elem.type;
		elem_node_num = elem.node_num;

		RC_TRY( set_element_direcor_vect(elem, d_vect, elem_vect) );
		RC_TRY( node_location_matrix(&elem, node, node_location) );
		RC_TRY( Gauss_const_default_shell(elem_type, g_points, g_weights,
		                                  &g_num_point) );
		RC_TRY( get_thickness(elem, physical, &thickness) );

		/* density */
		mat_index = search_material_prop_element(material, physical, &elem);
		if(mat_index < 0) return(SEEK_ERROR_RC);

		rho = material.array[mat_index].rho;

		/* gravity vector */
		grav_vect[0] = rho * accel_grav.x;
		grav_vect[1] = rho * accel_grav.y;
		grav_vect[2] = rho * accel_grav.z;

		for(ii2=0; ii2<elem_node_num*5; ii2++) elem_b_force[ii2] = 0.0;
		for(ii2=0; ii2<g_num_point; ii2++){
			RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num,
			                          elem_vect, thickness,
			                          node_location, g_points[ii2],
			                          G0, ws2) );
			RC_TRY( set_disp_N_matrix_shell(elem_type, elem_node_num, elem_vect,
			                                thickness, g_points[ii2],
			                                N_matrix, ws2) );

			mul_trans_matrix_vect(3, elem_node_num*5, N_matrix,
			                      grav_vect, Nt_grav);

			det_J = scalar_triple_product3d(G0[0], G0[1], G0[2]);

			for(ii3=0; ii3<elem_node_num*5; ii3++){
				elem_b_force[ii3] += Nt_grav[ii3] * det_J * g_weights[ii2];
			}
		}

		for(ii2=0; ii2<elem.node_num; ii2++){
			int index;
			VECT3D V1 = elem_vect[ii2].v1;
			VECT3D V2 = elem_vect[ii2].v2;

			index = node_renum_index(node, elem.node[ii2]);
			if(index < 0 ) return(SEEK_ERROR_RC);

			f_vect[index*dim] += elem_b_force[ii2*5];
			f_vect[index*dim+1] += elem_b_force[ii2*5+1];
			f_vect[index*dim+2] += elem_b_force[ii2*5+2];
			/* local component -> Global component */
			f_vect[index*dim+3] += V1.x*elem_b_force[ii2*5+3]
			                     + V2.x*elem_b_force[ii2*5+4];
			f_vect[index*dim+4] += V1.y*elem_b_force[ii2*5+3]
			                     + V2.y*elem_b_force[ii2*5+4];
			f_vect[index*dim+5] += V1.z*elem_b_force[ii2*5+3]
			                     + V2.z*elem_b_force[ii2*5+4];
		}
	}

	RC_TRY( free_work_set(&ws1) );
	RC_TRY( free_work_set(&ws2) );

	return(NORMAL_RC);
}


RC
add_temp_force_shell (NODE_ARRAY node, ELEMENT_ARRAY element,
                      MATERIAL_PROP_ARRAY material,
                      PHYSICAL_PROP_ARRAY physical,
                      BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                      DIRECTOR_VECT_ARRAY d_vect, double f_vect[])
{
	int ii1, ii2, ii3;
	int dim;
	int mat_index;
	int elem_node_num;
	int g_num_point;
	int ssnum;
	double det_J;
	double thickness;
	double local_temp;
	double strain_vect[5];
	double stress_vect[5];
	double react_vect[5*MAX_SHELL_NODE];
	double elem_temp_force[5*MAX_SHELL_NODE];
	double *elem_node_temp = NULL;
	double *N = NULL;
	double *g_weights = NULL;
	double **g_points = NULL;
	double **node_location = NULL;
	double **B_matrix = NULL;
	double **D_matrix = NULL;
	double **Jt = NULL;
	double **strain_tensor = NULL;
	double **cova_strain_tensor = NULL;
	VECT3D G0[3], g0[3];
	DIRECTOR_VECT elem_vect[MAX_SHELL_NODE];
	ELEM_TYPE elem_type;
	MATERIAL_PROP mat;
	WORK_SET ws1, ws2;
	MITC_TYING_SET mitc;


	if((dim = analysis_dim_shell(element)) != 6){
		return(ARG_ERROR_RC);
	}

	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

	if( (def_temp.size <= 0) && (temp.size <= 0) ) return(NORMAL_RC);

	RC_TRY( allocate_work_set(&ws1) );
	RC_TRY( allocate_work_set(&ws2) );
	RC_TRY( allocate_mitc_tying_set(&mitc) );

	elem_node_temp = ws1.vec_N[0];
	N = ws1.vec_N[1];
	g_weights = ws1.vec_G[0];
	node_location = ws1.mat_N_3[0];
	B_matrix = ws1.mat_6_3N[0];
	D_matrix = ws1.mat_6_6[0];
	g_points = ws1.mat_G_4[0];
	Jt = ws1.mat_3_3[0];
	strain_tensor = ws1.mat_3_3[1];
	cova_strain_tensor = ws1.mat_3_3[2];

	for(ii1=0; ii1<element.size; ii1++){
		ELEMENT elem = element.array[ii1];

		if(elem.label < 0) continue;
		if((ssnum = element_ssnum_shell(elem.type)) != 5 ){
			return(ARG_ERROR_RC);
		}

		elem_type = elem.type;
		elem_node_num = elem.node_num;

		RC_TRY( set_element_direcor_vect(elem, d_vect, elem_vect) );
		RC_TRY( node_location_matrix(&elem, node, node_location) );
		RC_TRY( Gauss_const_default_shell(elem_type, g_points, g_weights,
		                                  &g_num_point) );
		RC_TRY( get_thickness(elem, physical, &thickness) );

		/* nodal temperature */
		for(ii2=0; ii2<elem.node_num; ii2++){
			elem_node_temp[ii2] = node_temp(temp, def_temp, elem.node[ii2]);
		}

		/* set coefficient of thermal expansion */
		mat_index = search_material_prop_element(material, physical, &elem);
		if(mat_index < 0) return(SEEK_ERROR_RC);
		mat = material.array[mat_index];
		RC_TRY( fill_material(&mat) );

		for(ii2=0; ii2<elem_node_num*5; ii2++) elem_temp_force[ii2] = 0.0;
		for(ii2=0; ii2<g_num_point; ii2++){
			RC_TRY( make_G1G2G3g1g2g3_basis(elem_type, elem_node_num,
			                                elem_vect, thickness,
			                                node_location, g_points[ii2],
			                                G0, g0, ws2) );
			RC_TRY( make_B_matrix_shell(elem_type, elem_node_num,
			                            elem_vect, thickness, node_location,
			                            g_points[ii2], B_matrix, mitc, ws2) );
			RC_TRY( make_D_matrix_shell(elem, material, physical, g0,
			                            D_matrix) );
			
			/* set temperature on Gauss point */
			RC_TRY( set_N_vector_shell(elem_type, g_points[ii2], N) );
			local_temp = 0.0;
			for(ii3=0; ii3<elem.node_num; ii3++){
				local_temp += elem_node_temp[ii3] * N[ii3];
			}
			/* set temperature strain tensor */
			strain_tensor[0][0] = local_temp * mat.alpha_vect.x;
			strain_tensor[1][1] = local_temp * mat.alpha_vect.y;
			strain_tensor[0][1] = strain_tensor[1][0]
			                    = local_temp * mat.alpha_vect.xy;
			strain_tensor[0][2] = strain_tensor[2][0]
			                    = local_temp * mat.alpha_vect.zx;
			strain_tensor[1][2] = strain_tensor[2][1]
			                    = local_temp * mat.alpha_vect.yz;
			/* thickness is constant (ez = 0) */
			strain_tensor[2][2] = 0.0;

			/* set covariant temperature strain vector */
			Jt[0][0] = G0[0].x; Jt[1][0] = G0[0].y; Jt[2][0] = G0[0].z;
			Jt[0][1] = G0[1].x; Jt[1][1] = G0[1].y; Jt[2][1] = G0[1].z;
			Jt[0][2] = G0[2].x; Jt[1][2] = G0[2].y; Jt[2][2] = G0[2].z;
			mul_matrix_AtBA(3, 3, Jt, strain_tensor, cova_strain_tensor);
			strain_vect[0] = cova_strain_tensor[0][0];
			strain_vect[1] = cova_strain_tensor[1][1];
			strain_vect[2] = cova_strain_tensor[0][1];
			strain_vect[3] = cova_strain_tensor[0][2];
			strain_vect[4] = cova_strain_tensor[1][2];

			/* contravariant stress vector */
			mul_matrix_vect(ssnum, ssnum, D_matrix, strain_vect, stress_vect);
			/* temperatur force on Gauss point */
			mul_trans_matrix_vect(ssnum, 5*elem_node_num, B_matrix,
			                      stress_vect, react_vect);

			det_J = scalar_triple_product3d(G0[0], G0[1], G0[2]);

			for(ii3=0; ii3<elem_node_num*5; ii3++){
				elem_temp_force[ii3]
				+= react_vect[ii3] * det_J * g_weights[ii2];
			}
		}

		for(ii2=0; ii2<elem.node_num; ii2++){
			int index;
			VECT3D V1 = elem_vect[ii2].v1;
			VECT3D V2 = elem_vect[ii2].v2;

			index = node_renum_index(node, elem.node[ii2]);
			if(index < 0 ) return(SEEK_ERROR_RC);

			f_vect[index*dim] += elem_temp_force[ii2*5];
			f_vect[index*dim+1] += elem_temp_force[ii2*5+1];
			f_vect[index*dim+2] += elem_temp_force[ii2*5+2];
			/* local component -> Global component */
			f_vect[index*dim+3] += V1.x*elem_temp_force[ii2*5+3]
			                     + V2.x*elem_temp_force[ii2*5+4];
			f_vect[index*dim+4] += V1.y*elem_temp_force[ii2*5+3]
			                     + V2.y*elem_temp_force[ii2*5+4];
			f_vect[index*dim+5] += V1.z*elem_temp_force[ii2*5+3]
			                     + V2.z*elem_temp_force[ii2*5+4];
		}
	}

	RC_TRY( free_work_set(&ws1) );
	RC_TRY( free_work_set(&ws2) );
	RC_TRY( free_mitc_tying_set(&mitc) );

	return(NORMAL_RC);
}


/* シェルモデルの次元(1節点当たりの自由度) */
int
analysis_dim_shell (ELEMENT_ARRAY element)
{
	int dim, tmp_dim;
	int ii1;


	dim = 0;
	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		tmp_dim = element_dim_shell(element.array[ii1].type);
		if(dim == 0){
			dim = tmp_dim;
		}else if(dim != tmp_dim){
			return(-1);
		}
	}
	if(dim <= 0){
		return(-1);
	}

	return(dim);
}


/* シェル要素の次元(1節点当たりの自由度) */
int
element_dim_shell (ELEM_TYPE type)
{
	switch(type){
	case ELEM_TRI1:
	case ELEM_TRI2:
	case ELEM_QUAD1:
	case ELEM_QUAD2:
	case ELEM_QUAD2L:
		return(6);
	default:
		break;
	}

	return(-1);
}


/* シェル要素の応力 or ひずみ成分数(xx, yy, yz, zx, xy) */
int
element_ssnum_shell (ELEM_TYPE type)
{
	switch(type){
	case ELEM_TRI1:
	case ELEM_TRI2:
	case ELEM_QUAD1:
	case ELEM_QUAD2:
	case ELEM_QUAD2L:
		return(5);
	default:
		break;
	}

	return(-1);
}


/* 応力，ひずみ計算 */
RC
stress_strain_shell (int sol_ss_flag, int sol_ss_face,
                     NODE_ARRAY node, ELEMENT_ARRAY element,
                     MATERIAL_PROP_ARRAY material,
                     PHYSICAL_PROP_ARRAY physical,
                     DIRECTOR_VECT_ARRAY d_vect, DISP_ARRAY disp,
                     STRESS_ARRAY *stress, STRAIN_ARRAY *strain)
{
	int ii1, ii2;
	int node_index;
	int prop_index;
	int *count = NULL;
	double **local_points = NULL;
	STRESS_STRAIN local_stress, local_strain;
	VECT3DR init_v = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	VECT3DR *stress_sum = NULL;
	VECT3DR *strain_sum = NULL;
	WORK_SET ws1, ws2;
	MITC_TYING_SET mitc;


	if( ((stress == NULL)&&(strain == NULL)) || (sol_ss_flag == 0) ){
		return(NORMAL_RC);
	}

	RC_TRY( allocate_work_set(&ws1) );
	RC_TRY( allocate_work_set(&ws2) );
	RC_TRY( allocate_mitc_tying_set(&mitc) );
	RC_TRY( allocate2D(MAX_SHELL_NODE, 4, &local_points) );

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

		prop_index = search_physical_prop_label(physical,
		                                        element.array[ii1].physical);
		if(prop_index < 0) return(SEEK_ERROR_RC);
		if(NST_PROP_PSHELL != physical.array[prop_index].i_info[1]){
			return(ARG_ERROR_RC);
		}


		/* center of element */
		if(sol_ss_flag & SOL_SS_CENTER){
			RC_TRY( set_center_point_shell(sol_ss_face, element.array[ii1].type,
			                               local_points[0]) );

			RC_TRY( local_stress_strain_shell(element.array[ii1],
			                                  local_points[0], node, material,
			                                  physical, d_vect, disp,
			                                  &local_stress, &local_strain,
			                                  mitc, ws1, ws2) );

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
		/* each nodes */
		if( (sol_ss_flag & SOL_SS_NODE) || (sol_ss_flag & SOL_SS_NODE_AV) ){
			RC_TRY( set_local_node_points_shell(sol_ss_face,
			                                    element.array[ii1].type,
			                                    local_points));
			for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
				RC_TRY( local_stress_strain_shell(element.array[ii1],
				                                  local_points[ii2], node,
				                                  material, physical,
				                                  d_vect, disp,
				                                  &local_stress,
				                                  &local_strain,
				                                  mitc, ws1, ws2) );

				if(sol_ss_flag & SOL_SS_NODE_AV){
					node_index = search_node_label(node,
							element.array[ii1].node[ii2]);
					RC_NEG_CHK( node_index );
					count[node_index]++;
					stress_sum[node_index]=add_vect3dr(stress_sum[node_index],
					                                   local_stress.v);
					strain_sum[node_index]=add_vect3dr(strain_sum[node_index],
					                                   local_strain.v);
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

	RC_TRY( free_work_set(&ws1) );
	RC_TRY( free_work_set(&ws2) );
	RC_TRY( free_mitc_tying_set(&mitc) );
	RC_TRY( free2D(MAX_SHELL_NODE, 4, &local_points) );

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


/*
 * 局所座標系での応力，ひずみ計算
 * ひずみ不要の場合 : strain = NULL
 * 応力不要の場合   : stress = NULL
 */
RC
local_stress_strain_shell (ELEMENT element, const double *local_point,
                           NODE_ARRAY node,
                           MATERIAL_PROP_ARRAY material,
                           PHYSICAL_PROP_ARRAY physical,
                           DIRECTOR_VECT_ARRAY d_vect,
                           DISP_ARRAY disp,
                           STRESS_STRAIN *stress,
                           STRESS_STRAIN *strain,
                           MITC_TYING_SET mitc,
                           WORK_SET ws1, WORK_SET ws2)
{
	int dim, ssnum;
	int elem_node_num;
	int m_index;
	double thickness;
	double disp5dof_v[5*MAX_SHELL_NODE];
	double strain_vect[6];
	double stress_vect[6];
	double global_strain_vect[6];
	double global_stress_vect[6];
	double **node_location = ws1.mat_N_3[0];
	double **B_matrix = ws1.mat_6_3N[0];
	double **D_matrix = ws1.mat_6_6[0];
	double **Q_matrix = ws1.mat_6_6[1];
	DIRECTOR_VECT elem_vect[MAX_SHELL_NODE];
	ELEM_TYPE elem_type;
	MATERIAL_PROP tmp_mat;
	VECT3D G[3], g[3];


	if( (dim = element_dim_shell(element.type)) <= 0
	||(ssnum = element_ssnum_shell(element.type)) <= 0 ){
		return(ARG_ERROR_RC);
	}

	RC_TRY( get_thickness(element, physical, &thickness) );
	RC_TRY( reset_shell_element_type(element.type, &elem_type,
	                                 &elem_node_num) );

	RC_TRY( set_element_direcor_vect(element, d_vect, elem_vect) );
	RC_TRY( node_location_matrix(&element, node, node_location) );
	if(elem_type == ELEM_QUAD2L){
		RC_TRY( set_center_node_location(node_location, ws2) );
		RC_TRY( set_center_node_director_vect(elem_vect, ws2) );
	}

	RC_TRY( make_G1G2G3g1g2g3_basis(elem_type, elem_node_num,
	                                elem_vect, thickness, node_location,
	                                local_point, G, g, ws2) );

	/* strain */
	RC_TRY( make_B_matrix_shell(elem_type, elem_node_num,
	                            elem_vect, thickness, node_location,
	                            local_point, B_matrix, mitc, ws2) );
	RC_TRY( fill_disp5dof_vect_shell(element, elem_type, elem_node_num,
	                                 elem_vect, disp, disp5dof_v) );
	mul_matrix_vect(ssnum, 5*elem_node_num, B_matrix, disp5dof_v, strain_vect);
	/* 厚さ方向の垂直ひずみ */
	m_index = search_material_prop_element(material, physical, &element);
	if(m_index < 0) return(SEEK_ERROR_RC);
	tmp_mat = material.array[m_index];
	if(tmp_mat.mat_type != MAT_ISOTROPIC) return(ARG_ERROR_RC);
	strain_vect[5] = (-tmp_mat.nu / (1.0 + tmp_mat.nu))
	               * (inner_product3d(g[0], g[0])*strain_vect[0]
	               + inner_product3d(g[0], g[1])*strain_vect[2]
	               + inner_product3d(g[1], g[0])*strain_vect[2]
	               + inner_product3d(g[1], g[1])*strain_vect[1]); 

	RC_TRY( set_Q_matrix4strain(g, Q_matrix) );
	mul_matrix_vect(dim, dim, Q_matrix, strain_vect, global_strain_vect);

	/* stress */
	RC_TRY( make_D_matrix_shell(element, material, physical, g, D_matrix) );
	mul_matrix_vect(ssnum, ssnum, D_matrix, strain_vect, stress_vect);

	RC_TRY( set_Q_matrix4stress(G, Q_matrix) );
	mul_matrix_vect(dim, dim, Q_matrix, stress_vect, global_stress_vect);

	if(strain != NULL){
		init_stress_strain(strain);
		strain->element = element.label;
		strain->v.x =  global_strain_vect[0];
		strain->v.y =  global_strain_vect[1];
		strain->v.z =  global_strain_vect[2];
		strain->v.yz = global_strain_vect[3];
		strain->v.zx = global_strain_vect[4];
		strain->v.xy = global_strain_vect[5];
	}
	if(stress != NULL){
		init_stress_strain(stress);
		stress->element = element.label;
		stress->v.x =  global_stress_vect[0];
		stress->v.y =  global_stress_vect[1];
		stress->v.z =  global_stress_vect[2];
		stress->v.yz = global_stress_vect[3];
		stress->v.zx = global_stress_vect[4];
		stress->v.xy = global_stress_vect[5];
	}

	return(NORMAL_RC);
}


/* 反変基底での応力を全体座標系成分へ変換するマトリクス */
static RC
set_Q_matrix4stress (VECT3D G[], double **Q_matrix)
{
	double l1, l2, l3;
	double m1, m2, m3;
	double n1, n2, n3;


	l1 = G[0].x; l2 = G[0].y; l3 = G[0].z;
	m1 = G[1].x; m2 = G[1].y; m3 = G[1].z;
	n1 = G[2].x; n2 = G[2].y; n3 = G[2].z;

	/* xx */
	Q_matrix[0][0] = l1*l1;
	Q_matrix[0][1] = m1*m1;
	Q_matrix[0][2] = 2.0*l1*m1;
	Q_matrix[0][3] = 2.0*n1*l1;
	Q_matrix[0][4] = 2.0*m1*n1;
	Q_matrix[0][5] = n1*n1;

	/* yy */
	Q_matrix[1][0] = l2*l2;
	Q_matrix[1][1] = m2*m2;
	Q_matrix[1][2] = 2.0*l2*m2;
	Q_matrix[1][3] = 2.0*n2*l2;
	Q_matrix[1][4] = 2.0*m2*n2;
	Q_matrix[1][5] = n2*n2;

	/* zz */
	Q_matrix[2][0] = l3*l3;
	Q_matrix[2][1] = m3*m3;
	Q_matrix[2][2] = 2.0*l3*m3;
	Q_matrix[2][3] = 2.0*n3*l3;
	Q_matrix[2][4] = 2.0*m3*n3;
	Q_matrix[2][5] = n3*n3;

	/* yz */
	Q_matrix[3][0] = l2*l3;
	Q_matrix[3][1] = m2*m3;
	Q_matrix[3][2] = l2*m3 + l3*m2;
	Q_matrix[3][3] = n2*l3 + n3*l2;
	Q_matrix[3][4] = m2*n3 + m3*n2;
	Q_matrix[3][5] = n2*n3;

	/* zx */
	Q_matrix[4][0] = l3*l1;
	Q_matrix[4][1] = m3*m1;
	Q_matrix[4][2] = l3*m1 + l1*m3;
	Q_matrix[4][3] = n3*l1 + n1*l3;
	Q_matrix[4][4] = m3*n1 + m1*n3;
	Q_matrix[4][5] = n3*n1;

	/* xy */
	Q_matrix[5][0] = l1*l2;
	Q_matrix[5][1] = m1*m2;
	Q_matrix[5][2] = l1*m2 + l2*m1;
	Q_matrix[5][3] = n1*l2 + n2*l1;
	Q_matrix[5][4] = m1*n2 + m2*n1;
	Q_matrix[5][5] = n1*n2;

	return(NORMAL_RC);
}


/* 共変基底でのひずみを全体座標系成分へ変換するマトリクス */
static RC
set_Q_matrix4strain (VECT3D g[], double **Q_matrix)
{
	double l1, l2, l3;
	double m1, m2, m3;
	double n1, n2, n3;


	l1 = g[0].x; l2 = g[0].y; l3 = g[0].z;
	m1 = g[1].x; m2 = g[1].y; m3 = g[1].z;
	n1 = g[2].x; n2 = g[2].y; n3 = g[2].z;

	/* xx */
	Q_matrix[0][0] = l1*l1;
	Q_matrix[0][1] = m1*m1;
	Q_matrix[0][2] = l1*m1;
	Q_matrix[0][3] = n1*l1;
	Q_matrix[0][4] = m1*n1;
	Q_matrix[0][5] = n1*n1;

	/* yy */
	Q_matrix[1][0] = l2*l2;
	Q_matrix[1][1] = m2*m2;
	Q_matrix[1][2] = l2*m2;
	Q_matrix[1][3] = n2*l2;
	Q_matrix[1][4] = m2*n2;
	Q_matrix[1][5] = n2*n2;

	/* zz */
	Q_matrix[2][0] = l3*l3;
	Q_matrix[2][1] = m3*m3;
	Q_matrix[2][2] = l3*m3;
	Q_matrix[2][3] = n3*l3;
	Q_matrix[2][4] = m3*n3;
	Q_matrix[2][5] = n3*n3;

	/* yz */
	Q_matrix[3][0] = 2.0*l2*l3;
	Q_matrix[3][1] = 2.0*m2*m3;
	Q_matrix[3][2] = l2*m3 + l3*m2;
	Q_matrix[3][3] = n2*l3 + n3*l2;
	Q_matrix[3][4] = m2*n3 + m3*n2;
	Q_matrix[3][5] = 2.0*n2*n3;

	/* zx */
	Q_matrix[4][0] = 2.0*l3*l1;
	Q_matrix[4][1] = 2.0*m3*m1;
	Q_matrix[4][2] = l3*m1 + l1*m3;
	Q_matrix[4][3] = n3*l1 + n1*l3;
	Q_matrix[4][4] = m3*n1 + m1*n3;
	Q_matrix[4][5] = 2.0*n3*n1;

	/* xy */
	Q_matrix[5][0] = 2.0*l1*l2;
	Q_matrix[5][1] = 2.0*m1*m2;
	Q_matrix[5][2] = l1*m2 + l2*m1;
	Q_matrix[5][3] = n1*l2 + n2*l1;
	Q_matrix[5][4] = m1*n2 + m2*n1;
	Q_matrix[5][5] = 2.0*n1*n2;

	return(NORMAL_RC);
}


/* 1節点5自由度での要素節点変位ベクトル */
static RC
fill_disp5dof_vect_shell (ELEMENT element, ELEM_TYPE elem_type,
                          int elem_node_num, DIRECTOR_VECT elem_vect[],
                          DISP_ARRAY disp, double disp5dof_v[])
{
	int ii1;
	int dim;
	double disp_v[6*MAX_SHELL_NODE];
	VECT3D V1;
	VECT3D V2;


	if( (dim = element_dim_shell(element.type)) <= 0 ){
		return(ARG_ERROR_RC);
	}

	RC_TRY( fill_disp6dof_vect_shell(element, elem_type, elem_node_num,
	                                 elem_vect, disp, disp_v) );

	for(ii1=0; ii1<elem_node_num; ii1++){
		V1 = elem_vect[ii1].v1;
		V2 = elem_vect[ii1].v2;

		disp5dof_v[5*ii1]   = disp_v[dim*ii1+0];
		disp5dof_v[5*ii1+1] = disp_v[dim*ii1+1];
		disp5dof_v[5*ii1+2] = disp_v[dim*ii1+2];
		disp5dof_v[5*ii1+3] = V1.x*disp_v[dim*ii1+3]
		                    + V1.y*disp_v[dim*ii1+4]
		                    + V1.z*disp_v[dim*ii1+5];
		disp5dof_v[5*ii1+4] = V2.x*disp_v[dim*ii1+3]
		                    + V2.y*disp_v[dim*ii1+4]
		                    + V2.z*disp_v[dim*ii1+5];
	}

	return(NORMAL_RC);
}


/* 1節点6自由度での要素節点変位ベクトル */
static RC
fill_disp6dof_vect_shell (ELEMENT element, ELEM_TYPE elem_type,
                          int elem_node_num, DIRECTOR_VECT elem_vect[],
                          DISP_ARRAY disp, double disp6dof_v[])
{
	int ii1;
	int dim;
	double center_point[4];
	double N[MAX_SHELL_NODE];


	if( (dim = element_dim_shell(element.type)) <= 0 ){
		return(ARG_ERROR_RC);
	}

	for(ii1=0; ii1<element.node_num; ii1++){
		int index;

		RC_NEG_CHK( index = search_disp_node(disp, element.node[ii1]) );

		disp6dof_v[dim*ii1]   = disp.array[index].v.x;
		disp6dof_v[dim*ii1+1] = disp.array[index].v.y;
		disp6dof_v[dim*ii1+2] = disp.array[index].v.z;
		disp6dof_v[dim*ii1+3] = disp.array[index].v.yz;
		disp6dof_v[dim*ii1+4] = disp.array[index].v.zx;
		disp6dof_v[dim*ii1+5] = disp.array[index].v.xy;
	}

	if(elem_type == ELEM_QUAD2L){
		int bubble_index;

		bubble_index = dim * 8;
		disp6dof_v[bubble_index]   = 0.0;
		disp6dof_v[bubble_index+1] = 0.0;
		disp6dof_v[bubble_index+2] = 0.0;
		disp6dof_v[bubble_index+3] = 0.0;
		disp6dof_v[bubble_index+4] = 0.0;
		disp6dof_v[bubble_index+5] = 0.0;

		RC_TRY( set_center_point(ELEM_QUAD2, center_point) );
		RC_TRY( set_N_vector(ELEM_QUAD2, center_point, N) );
		for(ii1=0; ii1<8; ii1++){
			disp6dof_v[bubble_index]   += disp6dof_v[dim*ii1]   * N[ii1];
			disp6dof_v[bubble_index+1] += disp6dof_v[dim*ii1+1] * N[ii1];
			disp6dof_v[bubble_index+2] += disp6dof_v[dim*ii1+2] * N[ii1];
			disp6dof_v[bubble_index+3] += disp6dof_v[dim*ii1+3] * N[ii1];
			disp6dof_v[bubble_index+4] += disp6dof_v[dim*ii1+4] * N[ii1];
			disp6dof_v[bubble_index+5] += disp6dof_v[dim*ii1+5] * N[ii1];
		}
	}

	return(NORMAL_RC);
}


/*
 * 各節点での Director Vector を作成(vect のメモリも確保される)
 * 厚さ方向ベクトルは各要素の外向き単位法線を平均化したベクトル
 */
RC
make_nodal_director_vect (NODE_ARRAY node, ELEMENT_ARRAY element,
		                          DIRECTOR_VECT_ARRAY *vect)
{
	int ii1, ii2;
	int node_index;
	int *count = NULL;
	double **local_point = NULL;
	VECT3D local_v;
	VECT3D init_v = {0.0, 0.0, 0.0};
	VECT3D ey = {0.0, 1.0, 0.0};
	VECT3D ez = {0.0, 0.0, 1.0};
	VECT3D *vect_sum = NULL;
	WORK_SET ws;


	if(node.size <= 0 || element.size <= 0) return(ARG_ERROR_RC);
	/* DIRECTOR_VECT_ARRAY は NODE_ARRAY と1対1対応 */
	vect->size = node.size;
	vect->alloc_size = node.size;
	vect->sort_flag = node.sort_flag;
	vect->array = (DIRECTOR_VECT *)mm_alloc((node.size)*sizeof(DIRECTOR_VECT));
	if(vect->array == NULL) return(ALLOC_ERROR_RC);
	for(ii1=0; ii1<node.size; ii1++){
		vect->array[ii1].label = node.array[ii1].label;
		vect->array[ii1].v1 = init_v;
		vect->array[ii1].v2 = init_v;
		vect->array[ii1].vn = init_v;
	}

	RC_TRY( allocate_work_set(&ws) );
	local_point = ws.mat_N_4[0];
	init_matrix(MAX_SHELL_NODE, 4, local_point);

	vect_sum = (VECT3D *)mm_alloc((node.size)*sizeof(VECT3D));
	count = (int *)mm_alloc((node.size)*sizeof(int));
	if(vect_sum == NULL || count == NULL) return(ALLOC_ERROR_RC);
	for(ii1=0; ii1<node.size; ii1++){
		vect_sum[ii1] = init_v;
		count[ii1] = 0;
	}

	for(ii1=0; ii1<element.size; ii1++){
		ELEMENT elem = element.array[ii1];
		if(elem.label < 0) continue;

		RC_TRY( set_local_node_points(elem.type,local_point) );
		for(ii2=0; ii2<elem.node_num; ii2++){
			RC_TRY( set_unit_normal_vect_shell(local_point[ii2], elem, node,
			                                   &local_v, ws) );
			node_index = search_node_label(node, elem.node[ii2]);
			RC_NEG_CHK(node_index);
			count[node_index]++;
			vect_sum[node_index] = add_vect3d(vect_sum[node_index], local_v);
		}
	}

	for(ii1=0; ii1<vect->size; ii1++){
		double c;
		double abs_val;
		VECT3D v;

		if(count[ii1] == 0
		|| vect->array[ii1].label < 0){
			continue;
		}

		/* 厚さ方向ベクトル Vn */
		c = 1.0 / (double)count[ii1];
		v.x = vect_sum[ii1].x * c;
		v.y = vect_sum[ii1].y * c;
		v.z = vect_sum[ii1].z * c;
		vect->array[ii1].vn = unit_vect3d(v);

		/* V1 */
		v = outer_product3d(ey, vect->array[ii1].vn);
		abs_val = abs_vect3d(v);
		if(abs_val < ORTHOGONAL_ERROR){
			vect->array[ii1].v1 = ez;
		}else{
			c = 1.0 / abs_val;
			vect->array[ii1].v1.x = v.x * c;
			vect->array[ii1].v1.y = v.y * c;
			vect->array[ii1].v1.z = v.z * c;
		}

		/* V2 */
		vect->array[ii1].v2 = outer_product3d(vect->array[ii1].vn,
		                                      vect->array[ii1].v1);
	}

	RC_TRY( mm_free((void *)vect_sum) );
	RC_TRY( mm_free((void *)count) );

	RC_TRY( free_work_set(&ws) );

	return(NORMAL_RC);
}


/* 要素節点の厚み方向ベクトル elem_vect[] の設定 */
RC
set_element_direcor_vect(ELEMENT element, DIRECTOR_VECT_ARRAY vect,
                         DIRECTOR_VECT elem_vect[])
{
	int ii1;


	for(ii1=0; ii1<element.node_num; ii1++){
		int index = search_director_vect_label(vect, element.node[ii1]);
		if(index < 0) return(SEEK_ERROR_RC);

		elem_vect[ii1] = vect.array[index];
	}

	return(NORMAL_RC);
}


/* 節点番号から厚さ方向ベクトルの index 探査 */
int
search_director_vect_label (DIRECTOR_VECT_ARRAY vect, int label)
{
	int found_index = -1;


	if(vect.array == NULL || vect.size <= 0) return(-1);

	if(vect.sort_flag == ARRAY_SORTED){
		/* Binary search */
		int low = 0;
		int high = vect.size - 1;
		int middle;
		DIRECTOR_VECT *v = vect.array;

		while(low <= high){
			middle = (low + high) / 2;

			if(v[middle].label == label){
				found_index = middle;
				break;
			}else if(v[middle].label < label){
				low = middle + 1;
			}else{
				high = middle - 1;
			}
		}
	}else{
		/* Linear search */
		int ii1;
		DIRECTOR_VECT *v = vect.array;

		for(ii1=0; ii1<vect.size; ii1++){
			if(v[ii1].label == label){
				found_index = ii1;
				break;
			}
		}
	}

	return(found_index);
}


/* 中立面上の外向き単位法線ベクトル */
static RC
set_unit_normal_vect_shell (double *local_point, ELEMENT element,
                            NODE_ARRAY node, VECT3D *vect, WORK_SET ws)
{
	double abs_val;
	VECT3D d_xi, d_eta;
	double **dN = NULL;
	double **node_location = NULL;
	double **jacobi = NULL;


	dN = ws.mat_3_N[0];
	node_location = ws.mat_N_3[0];
	jacobi = ws.mat_3_3[0];

	init_matrix(3, MAX_SHELL_NODE, dN);
	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( set_dN_matrix(element.type, local_point, dN) );
	mul_matrix(3, element.node_num, 3, dN, node_location, jacobi);

	switch(element.type){
	case ELEM_TRI1:
	case ELEM_TRI2:
	case ELEM_QUAD1:
	case ELEM_QUAD2:
		d_xi.x = jacobi[0][0];
		d_xi.y = jacobi[0][1];
		d_xi.z = jacobi[0][2];
		d_eta.x = jacobi[1][0];
		d_eta.y = jacobi[1][1];
		d_eta.z = jacobi[1][2];
		*vect = outer_product3d(d_xi, d_eta);
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	abs_val = abs_vect3d(*vect);
	if(nearly_eq(abs_val, 0.0)) return(CAL_ERROR_RC);
	vect->x /= abs_val;
	vect->y /= abs_val;
	vect->z /= abs_val;

	return(NORMAL_RC);
}


/* ディレクターベクトルのメモリ確保 */
RC
allocate_director_vect_array (int num, DIRECTOR_VECT_ARRAY *vect)
{
	int ii1;
	VECT3D init_v = {0.0, 0.0, 0.0};


	if(num > 0){
		vect->array = mm_alloc(num * sizeof(DIRECTOR_VECT));
		if(vect->array == NULL) return(ALLOC_ERROR_RC);
	}else{
		vect->array = NULL;
	}

	vect->size = 0;
	vect->alloc_size = num;
	vect->sort_flag = ARRAY_UNSORTED;
	for(ii1=0; ii1<num; ii1++){
		vect->array[ii1].v1 = init_v;
		vect->array[ii1].v2 = init_v;
		vect->array[ii1].vn = init_v;
	}

	return(NORMAL_RC);
}


/* ディレクターベクトルのメモリ解放 */
RC
free_director_vect_array (DIRECTOR_VECT_ARRAY *vect)
{
	if(vect->alloc_size > 0) RC_TRY( mm_free(vect->array) );
	vect->array = NULL;
	vect->alloc_size = 0;
	vect->size = 0;

	return(NORMAL_RC);
}


/* ディレクターベクトルのコピー */
RC
copy_director_vect_array (DIRECTOR_VECT_ARRAY src,
                          DIRECTOR_VECT_ARRAY *dest)
{
	int ii1;
	int alloc_size_tmp;
	DIRECTOR_VECT *array_tmp;


	if(src.size > dest->alloc_size) return(OVERFLOW_ERROR_RC);
	if(dest->array == NULL) return(NULL_ERROR_RC);

	array_tmp = dest->array;
	alloc_size_tmp = dest->alloc_size;
	*dest = src;
	dest->array = array_tmp;
	dest->alloc_size = alloc_size_tmp;

	for(ii1=0; ii1<(dest->size); ii1++){
		dest->array[ii1] = src.array[ii1];
	}

	return(NORMAL_RC);
}


RC
allocate_mitc_tying_set (MITC_TYING_SET *mitc)
{
	int ii1;

	for(ii1=0; ii1<MITC_SIZE; ii1++){
		RC_TRY( allocate1D(MAX_SHELL_NODE, &(mitc->vec_N[ii1])) );
		RC_TRY( allocate1D(TYING_MAX, &(mitc->vec_T[ii1])) );
		RC_TRY( allocate2D(TYING_MAX, 3, &(mitc->mat_T_3[ii1])) );
		RC_TRY( allocate2D(3, MAX_SHELL_NODE, &(mitc->mat_3_N[ii1])) );
		RC_TRY( allocate2D(3, 5*MAX_SHELL_NODE, &(mitc->mat_3_5N[ii1])) );
		RC_TRY( allocate2D(5*MAX_SHELL_NODE, 5*MAX_SHELL_NODE,
		                   &(mitc->mat_5N_5N[ii1])) );
	}

	return(NORMAL_RC);
}


RC
free_mitc_tying_set (MITC_TYING_SET *mitc)
{
	int ii1;


	for(ii1=0; ii1<MITC_SIZE; ii1++){
		RC_TRY( free1D(MAX_SHELL_NODE, &(mitc->vec_N[ii1])) );
		RC_TRY( free1D(TYING_MAX, &(mitc->vec_T[ii1])) );
		RC_TRY( free2D(TYING_MAX, 3, &(mitc->mat_T_3[ii1])) );
		RC_TRY( free2D(3, MAX_SHELL_NODE, &(mitc->mat_3_N[ii1])) );
		RC_TRY( free2D(3, 5*MAX_SHELL_NODE, &(mitc->mat_3_5N[ii1])) );
		RC_TRY( free2D(5*MAX_SHELL_NODE, 5*MAX_SHELL_NODE,
		               &(mitc->mat_5N_5N[ii1])) );
	}

	return(NORMAL_RC);
}


/* 弾性テンソル */
RC
make_D_matrix_shell (ELEMENT element, MATERIAL_PROP_ARRAY material,
                     PHYSICAL_PROP_ARRAY physical, VECT3D g[],
                     double **D_matrix)
{
	int m_index;
	int ssnum;
	double lambda;
	double mu;
	double nu;
	double E;
	double C3311, C3333, C3322, C3312, C3323, C3331;
	double Cij11, Cij33, Cij22, Cij12, Cij23, Cij31;
	double div_C3333;
	double gij[3][3];
	MATERIAL_PROP tmp_mat;


	if( (ssnum = element_ssnum_shell(element.type)) <= 0) return(ARG_ERROR_RC);

	m_index = search_material_prop_element(material, physical, &element);
	if(m_index < 0) return(SEEK_ERROR_RC);

	tmp_mat = material.array[m_index];
	if(tmp_mat.mat_type != MAT_ISOTROPIC) return(ARG_ERROR_RC);

	E = tmp_mat.E;
	nu = tmp_mat.nu;
	lambda = E * nu / ((1.0+nu)*(1.0-2.0*nu));
	mu = E / (2.0 * (1 + nu));

	/* metric tensor gij */
	gij[0][0] = inner_product3d(g[0], g[0]);
	gij[1][1] = inner_product3d(g[1], g[1]);
	gij[2][2] = inner_product3d(g[2], g[2]);
	gij[0][1] = gij[1][0] = inner_product3d(g[0], g[1]);
	gij[0][2] = gij[2][0] = inner_product3d(g[0], g[2]);
	gij[1][2] = gij[2][1] = inner_product3d(g[1], g[2]);

	/* Initialization of D */
	init_matrix(ssnum, ssnum, D_matrix);

	C3311 = lambda*gij[2][2]*gij[0][0]
	      + mu*(gij[2][0]*gij[2][0] + gij[2][0]*gij[2][0]);
	C3333 = lambda*gij[2][2]*gij[2][2]
	      + mu*(gij[2][2]*gij[2][2] + gij[2][2]*gij[2][2]);
	C3322 = lambda*gij[2][2]*gij[1][1]
	      + mu*(gij[2][1]*gij[2][1] + gij[2][1]*gij[2][1]);
	C3312 = lambda*gij[2][2]*gij[0][1]
	      + mu*(gij[2][0]*gij[2][1] + gij[2][1]*gij[2][0]);
	C3323 = lambda*gij[2][2]*gij[1][2]
	      + mu*(gij[2][1]*gij[2][2] + gij[2][2]*gij[2][1]);
	C3331 = lambda*gij[2][2]*gij[2][0]
	      + mu*(gij[2][2]*gij[2][0] + gij[2][0]*gij[2][2]);
	div_C3333 = 1.0 / C3333;

	/* D[0]->rr, D[1]->ss. D[2]->rs, D[3]->rt, D[4]->st */
	/* rr 11 */
	Cij11 = lambda*gij[0][0]*gij[0][0]
	      + mu*(gij[0][0]*gij[0][0] + gij[0][0]*gij[0][0]);
	Cij33 = lambda*gij[0][0]*gij[2][2]
	      + mu*(gij[0][2]*gij[0][2] + gij[0][2]*gij[0][2]);
	Cij22 = lambda*gij[0][0]*gij[1][1]
	      + mu*(gij[0][1]*gij[0][1] + gij[0][1]*gij[0][1]);
	Cij12 = lambda*gij[0][0]*gij[0][1]
	      + mu*(gij[0][0]*gij[0][1] + gij[0][1]*gij[0][0]);
	Cij23 = lambda*gij[0][0]*gij[1][2]
	      + mu*(gij[0][1]*gij[0][2] + gij[0][2]*gij[0][1]);
	Cij31 = lambda*gij[0][0]*gij[2][0]
	      + mu*(gij[0][2]*gij[0][0] + gij[0][0]*gij[0][2]);
	D_matrix[0][0] = Cij11 - Cij33*C3311*div_C3333;
	D_matrix[0][1] = Cij22 - Cij33*C3322*div_C3333;
	D_matrix[0][2] = 0.5*((Cij12 + Cij12) - Cij33*(C3312 + C3312)*div_C3333);
	D_matrix[0][3] = 0.5*((Cij31 + Cij31) - Cij33*(C3331 + C3331)*div_C3333);
	D_matrix[0][4] = 0.5*((Cij23 + Cij23) - Cij33*(C3323 + C3323)*div_C3333);

	/* ss 22 */
	Cij11 = lambda*gij[1][1]*gij[0][0]
	      + mu*(gij[1][0]*gij[1][0] + gij[1][0]*gij[1][0]);
	Cij33 = lambda*gij[1][1]*gij[2][2]
	      + mu*(gij[1][2]*gij[1][2] + gij[1][2]*gij[1][2]);
	Cij22 = lambda*gij[1][1]*gij[1][1]
	      + mu*(gij[1][1]*gij[1][1] + gij[1][1]*gij[1][1]);
	Cij12 = lambda*gij[1][1]*gij[0][1]
	      + mu*(gij[1][0]*gij[1][1] + gij[1][1]*gij[1][0]);
	Cij23 = lambda*gij[1][1]*gij[1][2]
	      + mu*(gij[1][1]*gij[1][2] + gij[1][2]*gij[1][1]);
	Cij31 = lambda*gij[1][1]*gij[2][0]
	      + mu*(gij[1][2]*gij[1][0] + gij[1][0]*gij[1][2]);
	D_matrix[1][0] = Cij11 - Cij33*C3311*div_C3333;
	D_matrix[1][1] = Cij22 - Cij33*C3322*div_C3333;
	D_matrix[1][2] = 0.5*((Cij12 + Cij12) - Cij33*(C3312 + C3312)*div_C3333);
	D_matrix[1][3] = 0.5*((Cij31 + Cij31) - Cij33*(C3331 + C3331)*div_C3333);
	D_matrix[1][4] = 0.5*((Cij23 + Cij23) - Cij33*(C3323 + C3323)*div_C3333);

	/* rs 12 */
	Cij11 = lambda*gij[0][1]*gij[0][0]
	      + mu*(gij[0][0]*gij[1][0] + gij[0][0]*gij[1][0]);
	Cij33 = lambda*gij[0][1]*gij[2][2]
	      + mu*(gij[0][2]*gij[1][2] + gij[0][2]*gij[1][2]);
	Cij22 = lambda*gij[0][1]*gij[1][1]
	      + mu*(gij[0][1]*gij[1][1] + gij[0][1]*gij[1][1]);
	Cij12 = lambda*gij[0][1]*gij[0][1]
	      + mu*(gij[0][0]*gij[1][1] + gij[0][1]*gij[1][0]);
	Cij23 = lambda*gij[0][1]*gij[1][2]
	      + mu*(gij[0][1]*gij[1][2] + gij[0][2]*gij[1][1]);
	Cij31 = lambda*gij[0][1]*gij[2][0]
	      + mu*(gij[0][2]*gij[1][0] + gij[0][0]*gij[1][2]);
	D_matrix[2][0] = Cij11 - Cij33*C3311*div_C3333;
	D_matrix[2][1] = Cij22 - Cij33*C3322*div_C3333;
	D_matrix[2][2] = 0.5*((Cij12 + Cij12) - Cij33*(C3312 + C3312)*div_C3333);
	D_matrix[2][3] = 0.5*((Cij31 + Cij31) - Cij33*(C3331 + C3331)*div_C3333);
	D_matrix[2][4] = 0.5*((Cij23 + Cij23) - Cij33*(C3323 + C3323)*div_C3333);

	/* rt 31 */
	Cij11 = lambda*gij[2][0]*gij[0][0]
	      + mu*(gij[2][0]*gij[0][0] + gij[2][0]*gij[0][0]);
	Cij33 = lambda*gij[2][0]*gij[2][2]
	      + mu*(gij[2][2]*gij[0][2] + gij[2][2]*gij[0][2]);
	Cij22 = lambda*gij[2][0]*gij[1][1]
	      + mu*(gij[2][1]*gij[0][1] + gij[2][1]*gij[0][1]);
	Cij12 = lambda*gij[2][0]*gij[0][1]
	      + mu*(gij[2][0]*gij[0][1] + gij[2][1]*gij[0][0]);
	Cij23 = lambda*gij[2][0]*gij[1][2]
	      + mu*(gij[2][1]*gij[0][2] + gij[2][2]*gij[0][1]);
	Cij31 = lambda*gij[2][0]*gij[2][0]
	      + mu*(gij[2][2]*gij[0][0] + gij[2][0]*gij[0][2]);
	D_matrix[3][0] = Cij11 - Cij33*C3311*div_C3333;
	D_matrix[3][1] = Cij22 - Cij33*C3322*div_C3333;
	D_matrix[3][2] = 0.5*((Cij12 + Cij12) - Cij33*(C3312 + C3312)*div_C3333);
	D_matrix[3][3] = 0.5*((Cij31 + Cij31) - Cij33*(C3331 + C3331)*div_C3333);
	D_matrix[3][4] = 0.5*((Cij23 + Cij23) - Cij33*(C3323 + C3323)*div_C3333);

	/* st 23 */
	Cij11 = lambda*gij[1][2]*gij[0][0]
	      + mu*(gij[1][0]*gij[2][0] + gij[1][0]*gij[2][0]);
	Cij33 = lambda*gij[1][2]*gij[2][2]
	      + mu*(gij[1][2]*gij[2][2] + gij[1][2]*gij[2][2]);
	Cij22 = lambda*gij[1][2]*gij[1][1]
	      + mu*(gij[1][1]*gij[2][1] + gij[1][1]*gij[2][1]);
	Cij12 = lambda*gij[1][2]*gij[0][1]
	      + mu*(gij[1][0]*gij[2][1] + gij[1][1]*gij[2][0]);
	Cij23 = lambda*gij[1][2]*gij[1][2]
	      + mu*(gij[1][1]*gij[2][2] + gij[1][2]*gij[2][1]);
	Cij31 = lambda*gij[1][2]*gij[2][0]
	      + mu*(gij[1][2]*gij[2][0] + gij[1][0]*gij[2][2]);
	D_matrix[4][0] = Cij11 - Cij33*C3311*div_C3333;
	D_matrix[4][1] = Cij22 - Cij33*C3322*div_C3333;
	D_matrix[4][2] = 0.5*((Cij12 + Cij12) - Cij33*(C3312 + C3312)*div_C3333);
	D_matrix[4][3] = 0.5*((Cij31 + Cij31) - Cij33*(C3331 + C3331)*div_C3333);
	D_matrix[4][4] = 0.5*((Cij23 + Cij23) - Cij33*(C3323 + C3323)*div_C3333);

	return(NORMAL_RC);
}


/* 要素中心での局所座標 */
RC
set_center_point_shell (int sol_ss_face, ELEM_TYPE type, double *lcoord)
{
	switch(type){
	case ELEM_TRI1:
		lcoord[0] = 0.333333333333333;
		lcoord[1] = 0.333333333333333;
		break;
	case ELEM_TRI2:
		lcoord[0] = 0.333333333333333;
		lcoord[1] = 0.333333333333333;
		break;
	case ELEM_QUAD1:
		lcoord[0] = 0.0;
		lcoord[1] = 0.0;
		break;
	case ELEM_QUAD2:
	case ELEM_QUAD2L:
		lcoord[0] = 0.0;
		lcoord[1] = 0.0;
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}
	if(sol_ss_face == SOL_SS_BOTTOM){
		lcoord[2] = -1.0;
	}else if(sol_ss_face == SOL_SS_TOP){
		lcoord[2] = +1.0;
	}else{
		lcoord[2] = 0.0;
	}

	return(NORMAL_RC);
}


RC
set_local_node_points_shell (int sol_ss_face, ELEM_TYPE type, double **lcoord)
{
	double t;


	if(sol_ss_face == SOL_SS_BOTTOM){
		t = -1.0;
	}else if(sol_ss_face == SOL_SS_TOP){
		t = +1.0;
	}else{
		t = 0.0;
	}

	switch(type){
	case ELEM_TRI1:
		lcoord[0][0] = 0.0;
		lcoord[0][1] = 0.0;
		lcoord[1][0] = 1.0;
		lcoord[1][1] = 0.0;
		lcoord[2][0] = 0.0;
		lcoord[2][1] = 1.0;
		lcoord[0][2] = t;
		lcoord[1][2] = t;
		lcoord[2][2] = t;
		break;
	case ELEM_TRI2:
		lcoord[0][0] = 0.0;
		lcoord[0][1] = 0.0;
		lcoord[1][0] = 1.0;
		lcoord[1][1] = 0.0;
		lcoord[2][0] = 0.0;
		lcoord[2][1] = 1.0;
		lcoord[3][0] = 0.5;
		lcoord[3][1] = 0.5;
		lcoord[4][0] = 0.0;
		lcoord[4][1] = 0.5;
		lcoord[5][0] = 0.5;
		lcoord[5][1] = 0.0;
		lcoord[0][2] = t;
		lcoord[1][2] = t;
		lcoord[2][2] = t;
		lcoord[3][2] = t;
		lcoord[4][2] = t;
		lcoord[5][2] = t;
		break;
	case ELEM_QUAD1:
		lcoord[0][0] = -1.0;
		lcoord[0][1] = -1.0;
		lcoord[1][0] = +1.0;
		lcoord[1][1] = -1.0;
		lcoord[2][0] = +1.0;
		lcoord[2][1] = +1.0;
		lcoord[3][0] = -1.0;
		lcoord[3][1] = +1.0;
		lcoord[0][2] = t;
		lcoord[1][2] = t;
		lcoord[2][2] = t;
		lcoord[3][2] = t;
		break;
	case ELEM_QUAD2:
		lcoord[0][0] = -1.0;
		lcoord[0][1] = -1.0;
		lcoord[1][0] = +1.0;
		lcoord[1][1] = -1.0;
		lcoord[2][0] = +1.0;
		lcoord[2][1] = +1.0;
		lcoord[3][0] = -1.0;
		lcoord[3][1] = +1.0;
		lcoord[4][0] =  0.0;
		lcoord[4][1] = -1.0;
		lcoord[5][0] = +1.0;
		lcoord[5][1] =  0.0;
		lcoord[6][0] =  0.0;
		lcoord[6][1] = +1.0;
		lcoord[7][0] = -1.0;
		lcoord[7][1] =  0.0;
		lcoord[0][2] = t;
		lcoord[1][2] = t;
		lcoord[2][2] = t;
		lcoord[3][2] = t;
		lcoord[4][2] = t;
		lcoord[5][2] = t;
		lcoord[6][2] = t;
		lcoord[7][2] = t;
		break;
	case ELEM_QUAD2L:
		lcoord[0][0] = -1.0;
		lcoord[0][1] = -1.0;
		lcoord[1][0] = +1.0;
		lcoord[1][1] = -1.0;
		lcoord[2][0] = +1.0;
		lcoord[2][1] = +1.0;
		lcoord[3][0] = -1.0;
		lcoord[3][1] = +1.0;
		lcoord[4][0] =  0.0;
		lcoord[4][1] = -1.0;
		lcoord[5][0] = +1.0;
		lcoord[5][1] =  0.0;
		lcoord[6][0] =  0.0;
		lcoord[6][1] = +1.0;
		lcoord[7][0] = -1.0;
		lcoord[7][1] =  0.0;
		lcoord[8][0] =  0.0;
		lcoord[8][1] =  0.0;
		lcoord[0][2] = t;
		lcoord[1][2] = t;
		lcoord[2][2] = t;
		lcoord[3][2] = t;
		lcoord[4][2] = t;
		lcoord[5][2] = t;
		lcoord[6][2] = t;
		lcoord[7][2] = t;
		lcoord[8][2] = t;
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


/* 共変基底ベクトル */
RC
make_G1G2G3_basis (ELEM_TYPE elem_type, int elem_node_num,
                   DIRECTOR_VECT elem_vect[],
                   double thickness,
                   double **node_location,
                   const double *local_point,
                   VECT3D G[], WORK_SET ws)
{
	int ii1;
	double c1, c2;
	double t = local_point[2];
	double *N = ws.vec_N[0];
	double **dN = ws.mat_3_N[0];
	double **Jacobi = ws.mat_3_3[0];
	VECT3D V3;


	RC_TRY( set_N_vector_shell(elem_type, local_point, N) );
	RC_TRY( set_dN_matrix_shell(elem_type, local_point, dN) );

	init_matrix(3, 3, Jacobi);
	c1 = 0.5 * thickness * t;
	c2 = 0.5 * thickness;
	for(ii1=0; ii1<elem_node_num; ii1++){
		V3 = elem_vect[ii1].vn;

		Jacobi[0][0] += dN[0][ii1] * (node_location[ii1][0] + c1*V3.x);
		Jacobi[0][1] += dN[0][ii1] * (node_location[ii1][1] + c1*V3.y);
		Jacobi[0][2] += dN[0][ii1] * (node_location[ii1][2] + c1*V3.z);
		Jacobi[1][0] += dN[1][ii1] * (node_location[ii1][0] + c1*V3.x);
		Jacobi[1][1] += dN[1][ii1] * (node_location[ii1][1] + c1*V3.y);
		Jacobi[1][2] += dN[1][ii1] * (node_location[ii1][2] + c1*V3.z);
		Jacobi[2][0] += N[ii1] * c2 * V3.x;
		Jacobi[2][1] += N[ii1] * c2 * V3.y;
		Jacobi[2][2] += N[ii1] * c2 * V3.z;
	}

	G[0].x = Jacobi[0][0];
	G[0].y = Jacobi[0][1];
	G[0].z = Jacobi[0][2];
	G[1].x = Jacobi[1][0];
	G[1].y = Jacobi[1][1];
	G[1].z = Jacobi[1][2];
	G[2].x = Jacobi[2][0];
	G[2].y = Jacobi[2][1];
	G[2].z = Jacobi[2][2];

	return(NORMAL_RC);
}


/* 共変基底ベクトル，反変基底ベクトル */
RC
make_G1G2G3g1g2g3_basis (ELEM_TYPE elem_type, int elem_node_num,
                         DIRECTOR_VECT elem_vect[],
                         double thickness,
                         double **node_location,
                         const double *local_point,
                         VECT3D G[], VECT3D g[], WORK_SET ws)
{
	int ii1;
	double c1, c2;
	double t = local_point[2];
	double *N = ws.vec_N[0];
	double **dN = ws.mat_3_N[0];
	double **Jacobi = ws.mat_3_3[0];
	double **inverse_J = ws.mat_3_3[1];
	VECT3D V3;


	RC_TRY( set_N_vector_shell(elem_type, local_point, N) );
	RC_TRY( set_dN_matrix_shell(elem_type, local_point, dN) );

	init_matrix(3, 3, Jacobi);
	init_matrix(3, 3, inverse_J);
	c1 = 0.5 * thickness * t;
	c2 = 0.5 * thickness;
	for(ii1=0; ii1<elem_node_num; ii1++){
		V3 = elem_vect[ii1].vn;

		Jacobi[0][0] += dN[0][ii1] * (node_location[ii1][0] + c1*V3.x);
		Jacobi[0][1] += dN[0][ii1] * (node_location[ii1][1] + c1*V3.y);
		Jacobi[0][2] += dN[0][ii1] * (node_location[ii1][2] + c1*V3.z);
		Jacobi[1][0] += dN[1][ii1] * (node_location[ii1][0] + c1*V3.x);
		Jacobi[1][1] += dN[1][ii1] * (node_location[ii1][1] + c1*V3.y);
		Jacobi[1][2] += dN[1][ii1] * (node_location[ii1][2] + c1*V3.z);
		Jacobi[2][0] += N[ii1] * c2 * V3.x;
		Jacobi[2][1] += N[ii1] * c2 * V3.y;
		Jacobi[2][2] += N[ii1] * c2 * V3.z;
	}

	RC_TRY( inverse_matrix(3, Jacobi, inverse_J) );

	G[0].x = Jacobi[0][0];
	G[0].y = Jacobi[0][1];
	G[0].z = Jacobi[0][2];
	G[1].x = Jacobi[1][0];
	G[1].y = Jacobi[1][1];
	G[1].z = Jacobi[1][2];
	G[2].x = Jacobi[2][0];
	G[2].y = Jacobi[2][1];
	G[2].z = Jacobi[2][2];

	g[0].x = inverse_J[0][0];
	g[0].y = inverse_J[1][0];
	g[0].z = inverse_J[2][0];
	g[1].x = inverse_J[0][1];
	g[1].y = inverse_J[1][1];
	g[1].z = inverse_J[2][1];
	g[2].x = inverse_J[0][2];
	g[2].y = inverse_J[1][2];
	g[2].z = inverse_J[2][2];

	return(NORMAL_RC);
}


/* 反変基底ベクトル */
RC
make_g1g2g3_basis (ELEM_TYPE elem_type, int elem_node_num,
                   DIRECTOR_VECT elem_vect[],
                   double thickness,
                   double **node_location,
                   const double *local_point,
                   VECT3D g[], WORK_SET ws)
{
	int ii1;
	double c1, c2;
	double t = local_point[2];
	double *N = ws.vec_N[0];
	double **dN = ws.mat_3_N[0];
	double **Jacobi = ws.mat_3_3[0];
	double **inverse_J = ws.mat_3_3[1];
	VECT3D V3;


	RC_TRY( set_N_vector_shell(elem_type, local_point, N) );
	RC_TRY( set_dN_matrix_shell(elem_type, local_point, dN) );

	init_matrix(3, 3, Jacobi);
	init_matrix(3, 3, inverse_J);
	c1 = 0.5 * thickness * t;
	c2 = 0.5 * thickness;
	for(ii1=0; ii1<elem_node_num; ii1++){
		V3 = elem_vect[ii1].vn;

		Jacobi[0][0] += dN[0][ii1] * (node_location[ii1][0] + c1*V3.x);
		Jacobi[0][1] += dN[0][ii1] * (node_location[ii1][1] + c1*V3.y);
		Jacobi[0][2] += dN[0][ii1] * (node_location[ii1][2] + c1*V3.z);
		Jacobi[1][0] += dN[1][ii1] * (node_location[ii1][0] + c1*V3.x);
		Jacobi[1][1] += dN[1][ii1] * (node_location[ii1][1] + c1*V3.y);
		Jacobi[1][2] += dN[1][ii1] * (node_location[ii1][2] + c1*V3.z);
		Jacobi[2][0] += N[ii1] * c2 * V3.x;
		Jacobi[2][1] += N[ii1] * c2 * V3.y;
		Jacobi[2][2] += N[ii1] * c2 * V3.z;
	}

	RC_TRY( inverse_matrix(3, Jacobi, inverse_J) );

	g[0].x = inverse_J[0][0];
	g[0].y = inverse_J[1][0];
	g[0].z = inverse_J[2][0];
	g[1].x = inverse_J[0][1];
	g[1].y = inverse_J[1][1];
	g[1].z = inverse_J[2][1];
	g[2].x = inverse_J[0][2];
	g[2].y = inverse_J[1][2];
	g[2].z = inverse_J[2][2];

	return(NORMAL_RC);
}


/* シェル要素 変位-ひずみマトリクス(B matrix) */
RC
make_B_matrix_shell (ELEM_TYPE elem_type, int elem_node_num,
                     DIRECTOR_VECT elem_vect[], double thickness,
                     double **node_location, const double *local_point,
                     double **B_matrix, MITC_TYING_SET mitc, WORK_SET ws)
{
	switch(elem_type){
	case ELEM_TRI1:
		RC_TRY( make_B_matrix_mitc3(elem_type, elem_node_num, elem_vect,
		                            thickness, node_location, local_point,
		                            B_matrix, mitc, ws) );
		break;
	case ELEM_TRI2:
		RC_TRY( make_B_matrix_mitc6(elem_type, elem_node_num, elem_vect,
		                            thickness, node_location, local_point,
		                            B_matrix, mitc, ws) );
		break;
	case ELEM_QUAD1:
		RC_TRY( make_B_matrix_mitc4(elem_type, elem_node_num, elem_vect,
		                            thickness, node_location, local_point,
		                            B_matrix, mitc, ws) );
		break;
	case ELEM_QUAD2:
		RC_TRY( make_B_matrix_mitc8(elem_type, elem_node_num, elem_vect,
		                            thickness, node_location, local_point,
		                            B_matrix, mitc, ws) );
		break;
	case ELEM_QUAD2L:
		RC_TRY( make_B_matrix_mitc9(elem_type, elem_node_num, elem_vect,
		                            thickness, node_location, local_point,
		                            B_matrix, mitc, ws) );
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


/* シェル要素 変位-ひずみマトリクス(B matrix) MITC3 */
static RC
make_B_matrix_mitc3 (ELEM_TYPE elem_type, int elem_node_num,
                     DIRECTOR_VECT elem_vect[], double thickness,
                     double **node_location, const double *local_point,
                     double **B_matrix, MITC_TYING_SET mitc, WORK_SET ws)
{
	int ii1, ii2;
	int dim = 5;
	int ErtEst_tying_num;
	double t = local_point[2];
	double c1, c2;
	double V1G1, V1G2;
	double V2G1, V2G2;
	double V1G1rt_rt, V1G2rt_st, V1G3rt_rt, V1G3rt_st;
	double V2G1rt_rt, V2G2rt_st, V2G3rt_rt, V2G3rt_st;
	double V1G1st_rt, V1G2st_st, V1G3st_rt, V1G3st_st;
	double V2G1st_rt, V2G2st_st, V2G3st_rt, V2G3st_st;
	double *Ert_rt_h = mitc.vec_T[0];
	double *Ert_st_h = mitc.vec_T[1];
	double *Est_rt_h = mitc.vec_T[2];
	double *Est_st_h = mitc.vec_T[3];
	double *Nrt_rt = mitc.vec_N[0];
	double *Nrt_st = mitc.vec_N[1];
	double *Nst_rt = mitc.vec_N[2];
	double *Nst_st = mitc.vec_N[3];
	double **Ert_rt_tying_point = mitc.mat_T_3[0];
	double **Ert_st_tying_point = mitc.mat_T_3[1];
	double **Est_rt_tying_point = mitc.mat_T_3[2];
	double **Est_st_tying_point = mitc.mat_T_3[3];
	double **dN = mitc.mat_3_N[0];
	double **dNrt_rt = mitc.mat_3_N[1];
	double **dNrt_st = mitc.mat_3_N[2];
	double **dNst_rt = mitc.mat_3_N[3];
	double **dNst_st = mitc.mat_3_N[4];
	VECT3D G[3];
	VECT3D Grt_rt[3];
	VECT3D Grt_st[3];
	VECT3D Gst_rt[3];
	VECT3D Gst_st[3];
	VECT3D V1, V2;


	c1 = 0.5 * thickness * t;
	c2 = 0.5 * thickness;

	RC_TRY( set_dN_matrix_shell(elem_type, local_point, dN) );
	RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_vect, thickness,
                              node_location, local_point, G, ws) );
	for(ii1=0; ii1<elem_node_num; ii1++){
		V1 = elem_vect[ii1].v1;
		V2 = elem_vect[ii1].v2;
		V1G1 = inner_product3d(V1, G[0]);
		V1G2 = inner_product3d(V1, G[1]);
		V2G1 = inner_product3d(V2, G[0]);
		V2G2 = inner_product3d(V2, G[1]);

		B_matrix[0][dim*ii1]   = dN[0][ii1] * G[0].x;
		B_matrix[0][dim*ii1+1] = dN[0][ii1] * G[0].y;
		B_matrix[0][dim*ii1+2] = dN[0][ii1] * G[0].z;
		B_matrix[0][dim*ii1+3] = -dN[0][ii1] * c1 * V2G1;
		B_matrix[0][dim*ii1+4] = dN[0][ii1] * c1 * V1G1;

		B_matrix[1][dim*ii1]   = dN[1][ii1] * G[1].x;
		B_matrix[1][dim*ii1+1] = dN[1][ii1] * G[1].y;
		B_matrix[1][dim*ii1+2] = dN[1][ii1] * G[1].z;
		B_matrix[1][dim*ii1+3] = -dN[1][ii1] * c1 * V2G2;
		B_matrix[1][dim*ii1+4] = dN[1][ii1] * c1 * V1G2;

		B_matrix[2][dim*ii1]   = G[1].x*dN[0][ii1] + G[0].x*dN[1][ii1];
		B_matrix[2][dim*ii1+1] = G[1].y*dN[0][ii1] + G[0].y*dN[1][ii1];
		B_matrix[2][dim*ii1+2] = G[1].z*dN[0][ii1] + G[0].z*dN[1][ii1];
		B_matrix[2][dim*ii1+3] = -c1 * (dN[0][ii1]*V2G2 + dN[1][ii1]*V2G1);
		B_matrix[2][dim*ii1+4] = c1 * (dN[0][ii1]*V1G2 + dN[1][ii1]*V1G1);

		B_matrix[3][dim*ii1]   = 0.0;
		B_matrix[3][dim*ii1+1] = 0.0;
		B_matrix[3][dim*ii1+2] = 0.0;
		B_matrix[3][dim*ii1+3] = 0.0;
		B_matrix[3][dim*ii1+4] = 0.0;
		B_matrix[4][dim*ii1]   = 0.0;
		B_matrix[4][dim*ii1+1] = 0.0;
		B_matrix[4][dim*ii1+2] = 0.0;
		B_matrix[4][dim*ii1+3] = 0.0;
		B_matrix[4][dim*ii1+4] = 0.0;

	}

	RC_TRY( set_tying_data_mitc3(local_point, &ErtEst_tying_num,
	                             Ert_rt_tying_point, Ert_st_tying_point,
	                             Est_rt_tying_point, Est_st_tying_point,
	                             Ert_rt_h, Ert_st_h, Est_rt_h, Est_st_h) );

	for(ii1=0; ii1<ErtEst_tying_num; ii1++){
		RC_TRY(set_N_vector_shell(elem_type,Ert_rt_tying_point[ii1],Nrt_rt));
		RC_TRY(set_N_vector_shell(elem_type,Ert_st_tying_point[ii1],Nrt_st));
		RC_TRY(set_N_vector_shell(elem_type,Est_rt_tying_point[ii1],Nst_rt));
		RC_TRY(set_N_vector_shell(elem_type,Est_st_tying_point[ii1],Nst_st));
		RC_TRY(set_dN_matrix_shell(elem_type,Ert_rt_tying_point[ii1],dNrt_rt));
		RC_TRY(set_dN_matrix_shell(elem_type,Ert_st_tying_point[ii1],dNrt_st));
		RC_TRY(set_dN_matrix_shell(elem_type,Est_rt_tying_point[ii1],dNst_rt));
		RC_TRY(set_dN_matrix_shell(elem_type,Est_st_tying_point[ii1],dNst_st));
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_vect,
		                          thickness, node_location,
		                          Ert_rt_tying_point[ii1], Grt_rt, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_vect,
		                          thickness, node_location,
		                          Ert_st_tying_point[ii1], Grt_st, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_vect,
		                          thickness, node_location,
		                          Est_rt_tying_point[ii1], Gst_rt, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_vect,
		                          thickness, node_location,
		                          Est_st_tying_point[ii1], Gst_st, ws) );

		for(ii2=0; ii2<elem_node_num; ii2++){
			V1 = elem_vect[ii2].v1;
			V2 = elem_vect[ii2].v2;
			V1G1rt_rt = inner_product3d(V1, Grt_rt[0]);
			V1G2rt_st = inner_product3d(V1, Grt_st[1]);
			V1G3rt_rt = inner_product3d(V1, Grt_rt[2]);
			V1G3rt_st = inner_product3d(V1, Grt_st[2]);
			V2G1rt_rt = inner_product3d(V2, Grt_rt[0]);
			V2G2rt_st = inner_product3d(V2, Grt_st[1]);
			V2G3rt_rt = inner_product3d(V2, Grt_rt[2]);
			V2G3rt_st = inner_product3d(V2, Grt_st[2]);
			V1G1st_rt = inner_product3d(V1, Gst_rt[0]);
			V1G2st_st = inner_product3d(V1, Gst_st[1]);
			V1G3st_rt = inner_product3d(V1, Gst_rt[2]);
			V1G3st_st = inner_product3d(V1, Gst_st[2]);
			V2G1st_rt = inner_product3d(V2, Gst_rt[0]);
			V2G2st_st = inner_product3d(V2, Gst_st[1]);
			V2G3st_rt = inner_product3d(V2, Gst_rt[2]);
			V2G3st_st = inner_product3d(V2, Gst_st[2]);

			/* Brt */
			B_matrix[3][dim*ii2]  +=Ert_rt_h[ii1]*(dNrt_rt[0][ii2]*Grt_rt[2].x);
			B_matrix[3][dim*ii2+1]+=Ert_rt_h[ii1]*(dNrt_rt[0][ii2]*Grt_rt[2].y);
			B_matrix[3][dim*ii2+2]+=Ert_rt_h[ii1]*(dNrt_rt[0][ii2]*Grt_rt[2].z);
			B_matrix[3][dim*ii2+3]+=Ert_rt_h[ii1]*(-c2*Nrt_rt[ii2]*V2G1rt_rt
			                                     -c1*dNrt_rt[0][ii2]*V2G3rt_rt);
			B_matrix[3][dim*ii2+4]+=Ert_rt_h[ii1]*(c2*Nrt_rt[ii2]*V1G1rt_rt
			                                     +c1*dNrt_rt[0][ii2]*V1G3rt_rt);

			B_matrix[3][dim*ii2]  +=Ert_st_h[ii1]*(dNrt_st[1][ii2]*Grt_st[2].x);
			B_matrix[3][dim*ii2+1]+=Ert_st_h[ii1]*(dNrt_st[1][ii2]*Grt_st[2].y);
			B_matrix[3][dim*ii2+2]+=Ert_st_h[ii1]*(dNrt_st[1][ii2]*Grt_st[2].z);
			B_matrix[3][dim*ii2+3]+=Ert_st_h[ii1]*(-c2*Nrt_st[ii2]*V2G2rt_st
			                                     -c1*dNrt_st[1][ii2]*V2G3rt_st);
			B_matrix[3][dim*ii2+4]+=Ert_st_h[ii1]*(c2*Nrt_st[ii2]*V1G2rt_st
			                                     +c1*dNrt_st[1][ii2]*V1G3rt_st);

			/* Bst */
			B_matrix[4][dim*ii2]  +=Est_rt_h[ii1]*(dNst_rt[0][ii2]*Gst_rt[2].x);
			B_matrix[4][dim*ii2+1]+=Est_rt_h[ii1]*(dNst_rt[0][ii2]*Gst_rt[2].y);
			B_matrix[4][dim*ii2+2]+=Est_rt_h[ii1]*(dNst_rt[0][ii2]*Gst_rt[2].z);
			B_matrix[4][dim*ii2+3]+=Est_rt_h[ii1]*(-c2*Nst_rt[ii2]*V2G1st_rt
			                                     -c1*dNst_rt[0][ii2]*V2G3st_rt);
			B_matrix[4][dim*ii2+4]+=Est_rt_h[ii1]*(c2*Nst_rt[ii2]*V1G1st_rt
			                                     +c1*dNst_rt[0][ii2]*V1G3st_rt);

			B_matrix[4][dim*ii2]  +=Est_st_h[ii1]*(dNst_st[1][ii2]*Gst_st[2].x);
			B_matrix[4][dim*ii2+1]+=Est_st_h[ii1]*(dNst_st[1][ii2]*Gst_st[2].y);
			B_matrix[4][dim*ii2+2]+=Est_st_h[ii1]*(dNst_st[1][ii2]*Gst_st[2].z);
			B_matrix[4][dim*ii2+3]+=Est_st_h[ii1]*(-c2*Nst_st[ii2]*V2G2st_st
			                                     -c1*dNst_st[1][ii2]*V2G3st_st);
			B_matrix[4][dim*ii2+4]+=Est_st_h[ii1]*(c2*Nst_st[ii2]*V1G2st_st
			                                     +c1*dNst_st[1][ii2]*V1G3st_st);
		}
	}

	return(NORMAL_RC);
}


/* シェル要素 変位-ひずみマトリクス(B matrix) MITC6 */
static RC
make_B_matrix_mitc6 (ELEM_TYPE elem_type, int elem_node_num,
                     DIRECTOR_VECT elem_vect[], double thickness,
                     double **node_location, const double *local_point,
                     double **B_matrix, MITC_TYING_SET mitc, WORK_SET ws)
{
	int ii1, ii2;
	int dim = 5;
	int ssnum = 5;
	int ErrEss_tying_num;
	int Ers_tying_num;
	int ErtEst_tying_num;
	double t = local_point[2];
	double c1, c2;
	double V1G1rr, V2G1rr;
	double V1G2ss, V2G2ss;
	double V1G1rs_rr, V2G1rs_rr;
	double V1G2rs_ss, V2G2rs_ss;
	double V1G1rs_rs, V1G2rs_rs, V2G1rs_rs, V2G2rs_rs;
	double V1G1rt_rt, V1G2rt_st, V1G3rt_rt, V1G3rt_st;
	double V2G1rt_rt, V2G2rt_st, V2G3rt_rt, V2G3rt_st;
	double V1G1st_rt, V1G2st_st, V1G3st_rt, V1G3st_st;
	double V2G1st_rt, V2G2st_st, V2G3st_rt, V2G3st_st;
	double *Err_h = mitc.vec_T[0];
	double *Ess_h = mitc.vec_T[1];
	double *Ers_rr_h = mitc.vec_T[2];
	double *Ers_ss_h = mitc.vec_T[3];
	double *Ers_rs_h = mitc.vec_T[4];
	double *Ert_rt_h = mitc.vec_T[5];
	double *Ert_st_h = mitc.vec_T[6];
	double *Est_rt_h = mitc.vec_T[7];
	double *Est_st_h = mitc.vec_T[8];
	double *Nrt_rt = mitc.vec_N[0];
	double *Nrt_st = mitc.vec_N[1];
	double *Nst_rt = mitc.vec_N[2];
	double *Nst_st = mitc.vec_N[3];
	double **Err_tying_point = mitc.mat_T_3[0];
	double **Ess_tying_point = mitc.mat_T_3[1];
	double **Ers_rr_tying_point = mitc.mat_T_3[2];
	double **Ers_ss_tying_point = mitc.mat_T_3[3];
	double **Ers_rs_tying_point = mitc.mat_T_3[4];
	double **Ert_rt_tying_point = mitc.mat_T_3[5];
	double **Ert_st_tying_point = mitc.mat_T_3[6];
	double **Est_rt_tying_point = mitc.mat_T_3[7];
	double **Est_st_tying_point = mitc.mat_T_3[8];
	double **dNrr = mitc.mat_3_N[0];
	double **dNss = mitc.mat_3_N[1];
	double **dNrs_rr = mitc.mat_3_N[2];
	double **dNrs_ss = mitc.mat_3_N[3];
	double **dNrs_rs = mitc.mat_3_N[4];
	double **dNrt_rt = mitc.mat_3_N[5];
	double **dNrt_st = mitc.mat_3_N[6];
	double **dNst_rt = mitc.mat_3_N[7];
	double **dNst_st = mitc.mat_3_N[8];
	VECT3D Grr[3];
	VECT3D Gss[3];
	VECT3D Grs_rr[3];
	VECT3D Grs_ss[3];
	VECT3D Grs_rs[3];
	VECT3D Grt_rt[3];
	VECT3D Grt_st[3];
	VECT3D Gst_rt[3];
	VECT3D Gst_st[3];
	VECT3D V1, V2;


	c1 = 0.5 * thickness * t;
	c2 = 0.5 * thickness;

	for(ii1=0; ii1<ssnum; ii1++){
		for(ii2=0; ii2<elem_node_num; ii2++){
			B_matrix[ii1][dim*ii2]   = 0.0;
			B_matrix[ii1][dim*ii2+1] = 0.0;
			B_matrix[ii1][dim*ii2+2] = 0.0;
			B_matrix[ii1][dim*ii2+3] = 0.0;
			B_matrix[ii1][dim*ii2+4] = 0.0;
		}
	}

	RC_TRY( set_tying_data_mitc6(local_point, &ErrEss_tying_num,
	                             &Ers_tying_num, &ErtEst_tying_num,
	                             Err_tying_point, Ess_tying_point,
	                             Ers_rr_tying_point, Ers_ss_tying_point,
	                             Ers_rs_tying_point, Ert_rt_tying_point,
	                             Ert_st_tying_point, Est_rt_tying_point,
	                             Est_st_tying_point, Err_h, Ess_h, Ers_rr_h,
	                             Ers_ss_h, Ers_rs_h, Ert_rt_h, Ert_st_h,
	                             Est_rt_h, Est_st_h) );

	/* Err, Ess */
	for(ii1=0; ii1<ErrEss_tying_num; ii1++){
		RC_TRY( set_dN_matrix_shell(elem_type, Err_tying_point[ii1], dNrr) );
		RC_TRY( set_dN_matrix_shell(elem_type, Ess_tying_point[ii1], dNss) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_vect,
		                          thickness, node_location,
		                          Err_tying_point[ii1], Grr, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_vect,
		                          thickness, node_location,
		                          Ess_tying_point[ii1], Gss, ws) );

		for(ii2=0; ii2<elem_node_num; ii2++){
			V1 = elem_vect[ii2].v1;
			V2 = elem_vect[ii2].v2;
			V1G1rr = inner_product3d(V1, Grr[0]);
			V2G1rr = inner_product3d(V2, Grr[0]);
			V1G2ss = inner_product3d(V1, Gss[1]);
			V2G2ss = inner_product3d(V2, Gss[1]);

			B_matrix[0][dim*ii2]   += Err_h[ii1] * (dNrr[0][ii2] * Grr[0].x);
			B_matrix[0][dim*ii2+1] += Err_h[ii1] * (dNrr[0][ii2] * Grr[0].y);
			B_matrix[0][dim*ii2+2] += Err_h[ii1] * (dNrr[0][ii2] * Grr[0].z);
			B_matrix[0][dim*ii2+3] += Err_h[ii1] * (-dNrr[0][ii2]*c1*V2G1rr);
			B_matrix[0][dim*ii2+4] += Err_h[ii1] * (dNrr[0][ii2]*c1*V1G1rr);

			B_matrix[1][dim*ii2]   += Ess_h[ii1] * (dNss[1][ii2] * Gss[1].x);
			B_matrix[1][dim*ii2+1] += Ess_h[ii1] * (dNss[1][ii2] * Gss[1].y);
			B_matrix[1][dim*ii2+2] += Ess_h[ii1] * (dNss[1][ii2] * Gss[1].z);
			B_matrix[1][dim*ii2+3] += Ess_h[ii1] * (-dNss[1][ii2]*c1*V2G2ss);
			B_matrix[1][dim*ii2+4] += Ess_h[ii1] * (dNss[1][ii2]*c1*V1G2ss);
		}
	}

	/* Ers */
	for(ii1=0; ii1<elem_node_num; ii1++){
		B_matrix[2][dim*ii1]   = B_matrix[0][dim*ii1]
		                       + B_matrix[1][dim*ii1];
		B_matrix[2][dim*ii1+1] = B_matrix[0][dim*ii1+1]
		                       + B_matrix[1][dim*ii1+1];
		B_matrix[2][dim*ii1+2] = B_matrix[0][dim*ii1+2]
		                       + B_matrix[1][dim*ii1+2];
		B_matrix[2][dim*ii1+3] = B_matrix[0][dim*ii1+3]
		                       + B_matrix[1][dim*ii1+3];
		B_matrix[2][dim*ii1+4] = B_matrix[0][dim*ii1+4]
		                       + B_matrix[1][dim*ii1+4];
	}
	for(ii1=0; ii1<Ers_tying_num; ii1++){
		RC_TRY(set_dN_matrix_shell(elem_type,Ers_rr_tying_point[ii1],dNrs_rr));
		RC_TRY(set_dN_matrix_shell(elem_type,Ers_ss_tying_point[ii1],dNrs_ss));
		RC_TRY(set_dN_matrix_shell(elem_type,Ers_rs_tying_point[ii1],dNrs_rs));
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_vect,
		                          thickness, node_location,
		                          Ers_rr_tying_point[ii1], Grs_rr, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_vect,
		                          thickness, node_location,
		                          Ers_ss_tying_point[ii1], Grs_ss, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_vect,
		                          thickness, node_location,
		                          Ers_rs_tying_point[ii1], Grs_rs, ws) );
		for(ii2=0; ii2<elem_node_num; ii2++){
			V1 = elem_vect[ii2].v1;
			V2 = elem_vect[ii2].v2;
			V1G1rs_rr = inner_product3d(V1, Grs_rr[0]);
			V2G1rs_rr = inner_product3d(V2, Grs_rr[0]);
			V1G2rs_ss = inner_product3d(V1, Grs_ss[1]);
			V2G2rs_ss = inner_product3d(V2, Grs_ss[1]);
			V1G1rs_rs = inner_product3d(V1, Grs_rs[0]);
			V1G2rs_rs = inner_product3d(V1, Grs_rs[1]);
			V2G1rs_rs = inner_product3d(V2, Grs_rs[0]);
			V2G2rs_rs = inner_product3d(V2, Grs_rs[1]);

			B_matrix[2][dim*ii2]  +=Ers_rr_h[ii1]*(dNrs_rr[0][ii2]*Grs_rr[0].x);
			B_matrix[2][dim*ii2+1]+=Ers_rr_h[ii1]*(dNrs_rr[0][ii2]*Grs_rr[0].y);
			B_matrix[2][dim*ii2+2]+=Ers_rr_h[ii1]*(dNrs_rr[0][ii2]*Grs_rr[0].z);
			B_matrix[2][dim*ii2+3]
			+= Ers_rr_h[ii1]*(-dNrs_rr[0][ii2]*c1*V2G1rs_rr);
			B_matrix[2][dim*ii2+4]
			+= Ers_rr_h[ii1]*(dNrs_rr[0][ii2]*c1*V1G1rs_rr);

			B_matrix[2][dim*ii2]  +=Ers_ss_h[ii1]*(dNrs_ss[1][ii2]*Grs_ss[1].x);
			B_matrix[2][dim*ii2+1]+=Ers_ss_h[ii1]*(dNrs_ss[1][ii2]*Grs_ss[1].y);
			B_matrix[2][dim*ii2+2]+=Ers_ss_h[ii1]*(dNrs_ss[1][ii2]*Grs_ss[1].z);
			B_matrix[2][dim*ii2+3]
			+= Ers_ss_h[ii1]*(-dNrs_ss[1][ii2]*c1*V2G2rs_ss);
			B_matrix[2][dim*ii2+4]
			+= Ers_ss_h[ii1]*(dNrs_ss[1][ii2]*c1*V1G2rs_ss);

			B_matrix[2][dim*ii2]  += Ers_rs_h[ii1]*(Grs_rs[1].x*dNrs_rs[0][ii2]
			                       + Grs_rs[0].x*dNrs_rs[1][ii2]);
			B_matrix[2][dim*ii2+1]+= Ers_rs_h[ii1]*(Grs_rs[1].y*dNrs_rs[0][ii2]
			                       + Grs_rs[0].y*dNrs_rs[1][ii2]);
			B_matrix[2][dim*ii2+2]+= Ers_rs_h[ii1]*(Grs_rs[1].z*dNrs_rs[0][ii2]
			                       + Grs_rs[0].z*dNrs_rs[1][ii2]);
			B_matrix[2][dim*ii2+3]
			+= Ers_rs_h[ii1] * (-c1*(dNrs_rs[0][ii2]*V2G2rs_rs
			 + dNrs_rs[1][ii2]*V2G1rs_rs));
			B_matrix[2][dim*ii2+4]
			+= Ers_rs_h[ii1] * (c1*(dNrs_rs[0][ii2]*V1G2rs_rs
			 + dNrs_rs[1][ii2]*V1G1rs_rs));
		}
	}

	/* Ert, Est */
	for(ii1=0; ii1<ErtEst_tying_num; ii1++){
		RC_TRY(set_N_vector_shell(elem_type,Ert_rt_tying_point[ii1],Nrt_rt));
		RC_TRY(set_N_vector_shell(elem_type,Ert_st_tying_point[ii1],Nrt_st));
		RC_TRY(set_N_vector_shell(elem_type,Est_rt_tying_point[ii1],Nst_rt));
		RC_TRY(set_N_vector_shell(elem_type,Est_st_tying_point[ii1],Nst_st));
		RC_TRY(set_dN_matrix_shell(elem_type,Ert_rt_tying_point[ii1],dNrt_rt));
		RC_TRY(set_dN_matrix_shell(elem_type,Ert_st_tying_point[ii1],dNrt_st));
		RC_TRY(set_dN_matrix_shell(elem_type,Est_rt_tying_point[ii1],dNst_rt));
		RC_TRY(set_dN_matrix_shell(elem_type,Est_st_tying_point[ii1],dNst_st));
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_vect,
		                          thickness, node_location,
		                          Ert_rt_tying_point[ii1], Grt_rt, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_vect,
		                          thickness, node_location,
		                          Ert_st_tying_point[ii1], Grt_st, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_vect,
		                          thickness, node_location,
		                          Est_rt_tying_point[ii1], Gst_rt, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_vect,
		                          thickness, node_location,
		                          Est_st_tying_point[ii1], Gst_st, ws) );
		for(ii2=0; ii2<elem_node_num; ii2++){
			V1 = elem_vect[ii2].v1;
			V2 = elem_vect[ii2].v2;
			V1G1rt_rt = inner_product3d(V1, Grt_rt[0]);
			V1G2rt_st = inner_product3d(V1, Grt_st[1]);
			V1G3rt_rt = inner_product3d(V1, Grt_rt[2]);
			V1G3rt_st = inner_product3d(V1, Grt_st[2]);
			V2G1rt_rt = inner_product3d(V2, Grt_rt[0]);
			V2G2rt_st = inner_product3d(V2, Grt_st[1]);
			V2G3rt_rt = inner_product3d(V2, Grt_rt[2]);
			V2G3rt_st = inner_product3d(V2, Grt_st[2]);
			V1G1st_rt = inner_product3d(V1, Gst_rt[0]);
			V1G2st_st = inner_product3d(V1, Gst_st[1]);
			V1G3st_rt = inner_product3d(V1, Gst_rt[2]);
			V1G3st_st = inner_product3d(V1, Gst_st[2]);
			V2G1st_rt = inner_product3d(V2, Gst_rt[0]);
			V2G2st_st = inner_product3d(V2, Gst_st[1]);
			V2G3st_rt = inner_product3d(V2, Gst_rt[2]);
			V2G3st_st = inner_product3d(V2, Gst_st[2]);

			/* Brt */
			B_matrix[3][dim*ii2]  +=Ert_rt_h[ii1]*(dNrt_rt[0][ii2]*Grt_rt[2].x);
			B_matrix[3][dim*ii2+1]+=Ert_rt_h[ii1]*(dNrt_rt[0][ii2]*Grt_rt[2].y);
			B_matrix[3][dim*ii2+2]+=Ert_rt_h[ii1]*(dNrt_rt[0][ii2]*Grt_rt[2].z);
			B_matrix[3][dim*ii2+3]+=Ert_rt_h[ii1]*(-c2*Nrt_rt[ii2]*V2G1rt_rt
			                                     -c1*dNrt_rt[0][ii2]*V2G3rt_rt);
			B_matrix[3][dim*ii2+4]+=Ert_rt_h[ii1]*(c2*Nrt_rt[ii2]*V1G1rt_rt
			                                     +c1*dNrt_rt[0][ii2]*V1G3rt_rt);

			B_matrix[3][dim*ii2]  +=Ert_st_h[ii1]*(dNrt_st[1][ii2]*Grt_st[2].x);
			B_matrix[3][dim*ii2+1]+=Ert_st_h[ii1]*(dNrt_st[1][ii2]*Grt_st[2].y);
			B_matrix[3][dim*ii2+2]+=Ert_st_h[ii1]*(dNrt_st[1][ii2]*Grt_st[2].z);
			B_matrix[3][dim*ii2+3]+=Ert_st_h[ii1]*(-c2*Nrt_st[ii2]*V2G2rt_st
			                                     -c1*dNrt_st[1][ii2]*V2G3rt_st);
			B_matrix[3][dim*ii2+4]+=Ert_st_h[ii1]*(c2*Nrt_st[ii2]*V1G2rt_st
			                                     +c1*dNrt_st[1][ii2]*V1G3rt_st);

			/* Bst */
			B_matrix[4][dim*ii2]  +=Est_rt_h[ii1]*(dNst_rt[0][ii2]*Gst_rt[2].x);
			B_matrix[4][dim*ii2+1]+=Est_rt_h[ii1]*(dNst_rt[0][ii2]*Gst_rt[2].y);
			B_matrix[4][dim*ii2+2]+=Est_rt_h[ii1]*(dNst_rt[0][ii2]*Gst_rt[2].z);
			B_matrix[4][dim*ii2+3]+=Est_rt_h[ii1]*(-c2*Nst_rt[ii2]*V2G1st_rt
			                                     -c1*dNst_rt[0][ii2]*V2G3st_rt);
			B_matrix[4][dim*ii2+4]+=Est_rt_h[ii1]*(c2*Nst_rt[ii2]*V1G1st_rt
			                                     +c1*dNst_rt[0][ii2]*V1G3st_rt);

			B_matrix[4][dim*ii2]  +=Est_st_h[ii1]*(dNst_st[1][ii2]*Gst_st[2].x);
			B_matrix[4][dim*ii2+1]+=Est_st_h[ii1]*(dNst_st[1][ii2]*Gst_st[2].y);
			B_matrix[4][dim*ii2+2]+=Est_st_h[ii1]*(dNst_st[1][ii2]*Gst_st[2].z);
			B_matrix[4][dim*ii2+3]+=Est_st_h[ii1]*(-c2*Nst_st[ii2]*V2G2st_st
			                                     -c1*dNst_st[1][ii2]*V2G3st_st);
			B_matrix[4][dim*ii2+4]+=Est_st_h[ii1]*(c2*Nst_st[ii2]*V1G2st_st
			                                     +c1*dNst_st[1][ii2]*V1G3st_st);
		}
	}

	return(NORMAL_RC);
}


/* シェル要素 変位-ひずみマトリクス(B matrix) MITC4 */
static RC
make_B_matrix_mitc4 (ELEM_TYPE elem_type, int elem_node_num,
                     DIRECTOR_VECT elem_vect[], double thickness,
                     double **node_location, const double *local_point,
                     double **B_matrix, MITC_TYING_SET mitc, WORK_SET ws)
{
	int ii1, ii2;
	int dim = 5;
	int ErtEst_tying_num;
	double t = local_point[2];
	double c1, c2;
	double V1G1, V1G2;
	double V2G1, V2G2;
	double V1G1rt, V1G3rt;
	double V2G1rt, V2G3rt;
	double V1G2st, V1G3st;
	double V2G2st, V2G3st;
	double *Ert_h = mitc.vec_T[0];
	double *Est_h = mitc.vec_T[1];
	double *Nrt = mitc.vec_N[0];
	double *Nst = mitc.vec_N[1];
	double **Ert_tying_point = mitc.mat_T_3[0];
	double **Est_tying_point = mitc.mat_T_3[1];
	double **dN = mitc.mat_3_N[0];
	double **dNrt = mitc.mat_3_N[1];
	double **dNst = mitc.mat_3_N[2];
	VECT3D G[3];
	VECT3D Grt[3];
	VECT3D Gst[3];
	VECT3D V1, V2;


	c1 = 0.5 * thickness * t;
	c2 = 0.5 * thickness;

	RC_TRY( set_dN_matrix_shell(elem_type, local_point, dN) );
	RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_vect, thickness,
                              node_location, local_point, G, ws) );
	for(ii1=0; ii1<elem_node_num; ii1++){
		V1 = elem_vect[ii1].v1;
		V2 = elem_vect[ii1].v2;
		V1G1 = inner_product3d(V1, G[0]);
		V1G2 = inner_product3d(V1, G[1]);
		V2G1 = inner_product3d(V2, G[0]);
		V2G2 = inner_product3d(V2, G[1]);

		B_matrix[0][dim*ii1]   = dN[0][ii1] * G[0].x;
		B_matrix[0][dim*ii1+1] = dN[0][ii1] * G[0].y;
		B_matrix[0][dim*ii1+2] = dN[0][ii1] * G[0].z;
		B_matrix[0][dim*ii1+3] = -dN[0][ii1] * c1 * V2G1;
		B_matrix[0][dim*ii1+4] = dN[0][ii1] * c1 * V1G1;

		B_matrix[1][dim*ii1]   = dN[1][ii1] * G[1].x;
		B_matrix[1][dim*ii1+1] = dN[1][ii1] * G[1].y;
		B_matrix[1][dim*ii1+2] = dN[1][ii1] * G[1].z;
		B_matrix[1][dim*ii1+3] = -dN[1][ii1] * c1 * V2G2;
		B_matrix[1][dim*ii1+4] = dN[1][ii1] * c1 * V1G2;

		B_matrix[2][dim*ii1]   = G[1].x*dN[0][ii1] + G[0].x*dN[1][ii1];
		B_matrix[2][dim*ii1+1] = G[1].y*dN[0][ii1] + G[0].y*dN[1][ii1];
		B_matrix[2][dim*ii1+2] = G[1].z*dN[0][ii1] + G[0].z*dN[1][ii1];
		B_matrix[2][dim*ii1+3] = -c1 * (dN[0][ii1]*V2G2 + dN[1][ii1]*V2G1);
		B_matrix[2][dim*ii1+4] = c1 * (dN[0][ii1]*V1G2 + dN[1][ii1]*V1G1);

		B_matrix[3][dim*ii1]   = 0.0;
		B_matrix[3][dim*ii1+1] = 0.0;
		B_matrix[3][dim*ii1+2] = 0.0;
		B_matrix[3][dim*ii1+3] = 0.0;
		B_matrix[3][dim*ii1+4] = 0.0;
		B_matrix[4][dim*ii1]   = 0.0;
		B_matrix[4][dim*ii1+1] = 0.0;
		B_matrix[4][dim*ii1+2] = 0.0;
		B_matrix[4][dim*ii1+3] = 0.0;
		B_matrix[4][dim*ii1+4] = 0.0;

	}

	RC_TRY( set_tying_data_mitc4(local_point, &ErtEst_tying_num,
	                             Ert_tying_point, Est_tying_point,
	                             Ert_h, Est_h) );

	for(ii1=0; ii1<ErtEst_tying_num; ii1++){
		RC_TRY( set_N_vector_shell(elem_type, Ert_tying_point[ii1], Nrt) );
		RC_TRY( set_N_vector_shell(elem_type, Est_tying_point[ii1], Nst) );
		RC_TRY( set_dN_matrix_shell(elem_type, Ert_tying_point[ii1], dNrt) );
		RC_TRY( set_dN_matrix_shell(elem_type, Est_tying_point[ii1], dNst) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_vect,
		                          thickness, node_location,
		                          Ert_tying_point[ii1], Grt, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_vect,
		                          thickness, node_location,
		                          Est_tying_point[ii1], Gst, ws) );

		for(ii2=0; ii2<elem_node_num; ii2++){
			V1 = elem_vect[ii2].v1;
			V2 = elem_vect[ii2].v2;
			V1G1rt = inner_product3d(V1, Grt[0]);
			V1G3rt = inner_product3d(V1, Grt[2]);
			V2G1rt = inner_product3d(V2, Grt[0]);
			V2G3rt = inner_product3d(V2, Grt[2]);
			V1G2st = inner_product3d(V1, Gst[1]);
			V1G3st = inner_product3d(V1, Gst[2]);
			V2G2st = inner_product3d(V2, Gst[1]);
			V2G3st = inner_product3d(V2, Gst[2]);
			B_matrix[3][dim*ii2]   += Ert_h[ii1] * (dNrt[0][ii2] * Grt[2].x);
			B_matrix[3][dim*ii2+1] += Ert_h[ii1] * (dNrt[0][ii2] * Grt[2].y);
			B_matrix[3][dim*ii2+2] += Ert_h[ii1] * (dNrt[0][ii2] * Grt[2].z);
			B_matrix[3][dim*ii2+3] += Ert_h[ii1] * (-c2*Nrt[ii2]*V2G1rt
			                                        -c1*dNrt[0][ii2]*V2G3rt);
			B_matrix[3][dim*ii2+4] += Ert_h[ii1] * (c2*Nrt[ii2]*V1G1rt
			                                        + c1*dNrt[0][ii2]*V1G3rt);


			B_matrix[4][dim*ii2]   += Est_h[ii1] * (dNst[1][ii2] * Gst[2].x);
			B_matrix[4][dim*ii2+1] += Est_h[ii1] * (dNst[1][ii2] * Gst[2].y);
			B_matrix[4][dim*ii2+2] += Est_h[ii1] * (dNst[1][ii2] * Gst[2].z);
			B_matrix[4][dim*ii2+3] += Est_h[ii1] * (-c2*Nst[ii2]*V2G2st
			                                        -c1*dNst[1][ii2]*V2G3st);
			B_matrix[4][dim*ii2+4] += Est_h[ii1] * (c2*Nst[ii2]*V1G2st
			                                        + c1*dNst[1][ii2]*V1G3st);
		}
	}

	return(NORMAL_RC);
}


/* シェル要素 変位-ひずみマトリクス(B matrix) MITC8 */
static RC
make_B_matrix_mitc8 (ELEM_TYPE elem_type, int elem_node_num,
                     DIRECTOR_VECT elem_vect[], double thickness,
                     double **node_location, const double *local_point,
                     double **B_matrix, MITC_TYING_SET mitc, WORK_SET ws)
{
	int ii1, ii2;
	int dim = 5;
	int ssnum = 5;
	int Ers_tying_num;
	int ErrEssErtEst_tying_num;
	double t = local_point[2];
	double c1, c2;
	double V1G1rr, V2G1rr;
	double V1G2ss, V2G2ss;
	double V1G1rs, V1G2rs;
	double V2G1rs, V2G2rs;
	double V1G1rt, V1G3rt;
	double V2G1rt, V2G3rt;
	double V1G2st, V1G3st;
	double V2G2st, V2G3st;
	double *Err_h = mitc.vec_T[0];
	double *Ess_h = mitc.vec_T[1];
	double *Ers_h = mitc.vec_T[2];
	double *Ert_h = mitc.vec_T[3];
	double *Est_h = mitc.vec_T[4];
	double *Nrt = mitc.vec_N[0];
	double *Nst = mitc.vec_N[1];
	double **Err_tying_point = mitc.mat_T_3[0];
	double **Ess_tying_point = mitc.mat_T_3[1];
	double **Ers_tying_point = mitc.mat_T_3[2];
	double **Ert_tying_point = mitc.mat_T_3[3];
	double **Est_tying_point = mitc.mat_T_3[4];
	double **dNrr = mitc.mat_3_N[0];
	double **dNss = mitc.mat_3_N[1];
	double **dNrs = mitc.mat_3_N[2];
	double **dNrt = mitc.mat_3_N[3];
	double **dNst = mitc.mat_3_N[4];
	VECT3D Grr[3];
	VECT3D Gss[3];
	VECT3D Grs[3];
	VECT3D Grt[3];
	VECT3D Gst[3];
	VECT3D V1, V2;


	c1 = 0.5 * thickness * t;
	c2 = 0.5 * thickness;

	for(ii1=0; ii1<ssnum; ii1++){
		for(ii2=0; ii2<elem_node_num; ii2++){
			B_matrix[ii1][dim*ii2]   = 0.0;
			B_matrix[ii1][dim*ii2+1] = 0.0;
			B_matrix[ii1][dim*ii2+2] = 0.0;
			B_matrix[ii1][dim*ii2+3] = 0.0;
			B_matrix[ii1][dim*ii2+4] = 0.0;
		}
	}

	RC_TRY( set_tying_data_mitc8(local_point, &Ers_tying_num,
	                             &ErrEssErtEst_tying_num, Err_tying_point,
	                             Ess_tying_point, Ers_tying_point,
	                             Ert_tying_point, Est_tying_point, Err_h,
	                             Ess_h, Ers_h, Ert_h, Est_h) );

	for(ii1=0; ii1<Ers_tying_num; ii1++){
		RC_TRY( set_dN_matrix_shell(elem_type, Ers_tying_point[ii1], dNrs) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_vect,
		                          thickness, node_location,
		                          Ers_tying_point[ii1], Grs, ws) );
		for(ii2=0; ii2<elem_node_num; ii2++){
			V1 = elem_vect[ii2].v1;
			V2 = elem_vect[ii2].v2;
			V1G1rs = inner_product3d(V1, Grs[0]);
			V1G2rs = inner_product3d(V1, Grs[1]);
			V2G1rs = inner_product3d(V2, Grs[0]);
			V2G2rs = inner_product3d(V2, Grs[1]);

			B_matrix[2][dim*ii2]   += Ers_h[ii1] * (Grs[1].x*dNrs[0][ii2]
			                                      + Grs[0].x*dNrs[1][ii2]);
			B_matrix[2][dim*ii2+1] += Ers_h[ii1] * (Grs[1].y*dNrs[0][ii2]
			                                      + Grs[0].y*dNrs[1][ii2]);
			B_matrix[2][dim*ii2+2] += Ers_h[ii1] * (Grs[1].z*dNrs[0][ii2]
			                                      + Grs[0].z*dNrs[1][ii2]);
			B_matrix[2][dim*ii2+3] += Ers_h[ii1] * (-c1*(dNrs[0][ii2]*V2G2rs
			                                           + dNrs[1][ii2]*V2G1rs));
			B_matrix[2][dim*ii2+4] += Ers_h[ii1] * (c1*(dNrs[0][ii2]*V1G2rs
			                                          + dNrs[1][ii2]*V1G1rs));

		}
	}

	/* Err, Ert, Ess, Est */
	for(ii1=0; ii1<ErrEssErtEst_tying_num; ii1++){
		RC_TRY( set_N_vector_shell(elem_type, Ert_tying_point[ii1], Nrt) );
		RC_TRY( set_N_vector_shell(elem_type, Est_tying_point[ii1], Nst) );
		RC_TRY( set_dN_matrix_shell(elem_type, Err_tying_point[ii1], dNrr) );
		RC_TRY( set_dN_matrix_shell(elem_type, Ess_tying_point[ii1], dNss) );
		RC_TRY( set_dN_matrix_shell(elem_type, Ert_tying_point[ii1], dNrt) );
		RC_TRY( set_dN_matrix_shell(elem_type, Est_tying_point[ii1], dNst) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_vect,
		                          thickness, node_location,
		                          Err_tying_point[ii1], Grr, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_vect,
		                          thickness, node_location,
		                          Ess_tying_point[ii1], Gss, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_vect,
		                          thickness, node_location,
		                          Ert_tying_point[ii1], Grt, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_vect,
		                          thickness, node_location,
		                          Est_tying_point[ii1], Gst, ws) );

		for(ii2=0; ii2<elem_node_num; ii2++){
			V1 = elem_vect[ii2].v1;
			V2 = elem_vect[ii2].v2;
			V1G1rr = inner_product3d(V1, Grr[0]);
			V2G1rr = inner_product3d(V2, Grr[0]);
			V1G2ss = inner_product3d(V1, Gss[1]);
			V2G2ss = inner_product3d(V2, Gss[1]);
			V1G1rt = inner_product3d(V1, Grt[0]);
			V1G3rt = inner_product3d(V1, Grt[2]);
			V2G1rt = inner_product3d(V2, Grt[0]);
			V2G3rt = inner_product3d(V2, Grt[2]);
			V1G2st = inner_product3d(V1, Gst[1]);
			V1G3st = inner_product3d(V1, Gst[2]);
			V2G2st = inner_product3d(V2, Gst[1]);
			V2G3st = inner_product3d(V2, Gst[2]);

			B_matrix[0][dim*ii2]   += Err_h[ii1] * (dNrr[0][ii2] * Grr[0].x);
			B_matrix[0][dim*ii2+1] += Err_h[ii1] * (dNrr[0][ii2] * Grr[0].y);
			B_matrix[0][dim*ii2+2] += Err_h[ii1] * (dNrr[0][ii2] * Grr[0].z);
			B_matrix[0][dim*ii2+3] += Err_h[ii1] * (-dNrr[0][ii2]*c1*V2G1rr);
			B_matrix[0][dim*ii2+4] += Err_h[ii1] * (dNrr[0][ii2]*c1*V1G1rr);

			B_matrix[1][dim*ii2]   += Ess_h[ii1] * (dNss[1][ii2] * Gss[1].x);
			B_matrix[1][dim*ii2+1] += Ess_h[ii1] * (dNss[1][ii2] * Gss[1].y);
			B_matrix[1][dim*ii2+2] += Ess_h[ii1] * (dNss[1][ii2] * Gss[1].z);
			B_matrix[1][dim*ii2+3] += Ess_h[ii1] * (-dNss[1][ii2]*c1*V2G2ss);
			B_matrix[1][dim*ii2+4] += Ess_h[ii1] * (dNss[1][ii2]*c1*V1G2ss);

			B_matrix[3][dim*ii2]   += Ert_h[ii1] * (dNrt[0][ii2] * Grt[2].x);
			B_matrix[3][dim*ii2+1] += Ert_h[ii1] * (dNrt[0][ii2] * Grt[2].y);
			B_matrix[3][dim*ii2+2] += Ert_h[ii1] * (dNrt[0][ii2] * Grt[2].z);
			B_matrix[3][dim*ii2+3] += Ert_h[ii1] * (-c2*Nrt[ii2]*V2G1rt
			                                        -c1*dNrt[0][ii2]*V2G3rt);
			B_matrix[3][dim*ii2+4] += Ert_h[ii1] * (c2*Nrt[ii2]*V1G1rt
			                                        + c1*dNrt[0][ii2]*V1G3rt);

			B_matrix[4][dim*ii2]   += Est_h[ii1] * (dNst[1][ii2] * Gst[2].x);
			B_matrix[4][dim*ii2+1] += Est_h[ii1] * (dNst[1][ii2] * Gst[2].y);
			B_matrix[4][dim*ii2+2] += Est_h[ii1] * (dNst[1][ii2] * Gst[2].z);
			B_matrix[4][dim*ii2+3] += Est_h[ii1] * (-c2*Nst[ii2]*V2G2st
			                                        -c1*dNst[1][ii2]*V2G3st);
			B_matrix[4][dim*ii2+4] += Est_h[ii1] * (c2*Nst[ii2]*V1G2st
			                                        + c1*dNst[1][ii2]*V1G3st);
		}
	}

	return(NORMAL_RC);
}


/* シェル要素 変位-ひずみマトリクス(B matrix) MITC9 */
static RC
make_B_matrix_mitc9 (ELEM_TYPE elem_type, int elem_node_num,
                     DIRECTOR_VECT elem_vect[], double thickness,
                     double **node_location, const double *local_point,
                     double **B_matrix, MITC_TYING_SET mitc, WORK_SET ws)
{
	int ii1, ii2;
	int dim = 5;
	int ssnum = 5;
	int Ers_tying_num;
	int ErrEssErtEst_tying_num;
	double t = local_point[2];
	double c1, c2;
	double V1G1rr, V2G1rr;
	double V1G2ss, V2G2ss;
	double V1G1rs, V1G2rs;
	double V2G1rs, V2G2rs;
	double V1G1rt, V1G3rt;
	double V2G1rt, V2G3rt;
	double V1G2st, V1G3st;
	double V2G2st, V2G3st;
	double *Err_h = mitc.vec_T[0];
	double *Ess_h = mitc.vec_T[1];
	double *Ers_h = mitc.vec_T[2];
	double *Ert_h = mitc.vec_T[3];
	double *Est_h = mitc.vec_T[4];
	double *Nrt = mitc.vec_N[0];
	double *Nst = mitc.vec_N[1];
	double **Err_tying_point = mitc.mat_T_3[0];
	double **Ess_tying_point = mitc.mat_T_3[1];
	double **Ers_tying_point = mitc.mat_T_3[2];
	double **Ert_tying_point = mitc.mat_T_3[3];
	double **Est_tying_point = mitc.mat_T_3[4];
	double **dNrr = mitc.mat_3_N[0];
	double **dNss = mitc.mat_3_N[1];
	double **dNrs = mitc.mat_3_N[2];
	double **dNrt = mitc.mat_3_N[3];
	double **dNst = mitc.mat_3_N[4];
	VECT3D Grr[3];
	VECT3D Gss[3];
	VECT3D Grs[3];
	VECT3D Grt[3];
	VECT3D Gst[3];
	VECT3D V1, V2;


	c1 = 0.5 * thickness * t;
	c2 = 0.5 * thickness;

	for(ii1=0; ii1<ssnum; ii1++){
		for(ii2=0; ii2<elem_node_num; ii2++){
			B_matrix[ii1][dim*ii2]   = 0.0;
			B_matrix[ii1][dim*ii2+1] = 0.0;
			B_matrix[ii1][dim*ii2+2] = 0.0;
			B_matrix[ii1][dim*ii2+3] = 0.0;
			B_matrix[ii1][dim*ii2+4] = 0.0;
		}
	}

	RC_TRY( set_tying_data_mitc9(local_point, &Ers_tying_num,
	                             &ErrEssErtEst_tying_num, Err_tying_point,
	                             Ess_tying_point, Ers_tying_point,
	                             Ert_tying_point, Est_tying_point, Err_h,
	                             Ess_h, Ers_h, Ert_h, Est_h) );

	for(ii1=0; ii1<Ers_tying_num; ii1++){
		RC_TRY( set_dN_matrix_shell(elem_type, Ers_tying_point[ii1], dNrs) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_vect,
		                          thickness, node_location,
		                          Ers_tying_point[ii1], Grs, ws) );
		for(ii2=0; ii2<elem_node_num; ii2++){
			V1 = elem_vect[ii2].v1;
			V2 = elem_vect[ii2].v2;
			V1G1rs = inner_product3d(V1, Grs[0]);
			V1G2rs = inner_product3d(V1, Grs[1]);
			V2G1rs = inner_product3d(V2, Grs[0]);
			V2G2rs = inner_product3d(V2, Grs[1]);

			B_matrix[2][dim*ii2]   += Ers_h[ii1] * (Grs[1].x*dNrs[0][ii2]
			                                      + Grs[0].x*dNrs[1][ii2]);
			B_matrix[2][dim*ii2+1] += Ers_h[ii1] * (Grs[1].y*dNrs[0][ii2]
			                                      + Grs[0].y*dNrs[1][ii2]);
			B_matrix[2][dim*ii2+2] += Ers_h[ii1] * (Grs[1].z*dNrs[0][ii2]
			                                      + Grs[0].z*dNrs[1][ii2]);
			B_matrix[2][dim*ii2+3] += Ers_h[ii1] * (-c1*(dNrs[0][ii2]*V2G2rs
			                                           + dNrs[1][ii2]*V2G1rs));
			B_matrix[2][dim*ii2+4] += Ers_h[ii1] * (c1*(dNrs[0][ii2]*V1G2rs
			                                          + dNrs[1][ii2]*V1G1rs));

		}
	}

	/* Err, Ert, Ess, Est */
	for(ii1=0; ii1<ErrEssErtEst_tying_num; ii1++){
		RC_TRY( set_N_vector_shell(elem_type, Ert_tying_point[ii1], Nrt) );
		RC_TRY( set_N_vector_shell(elem_type, Est_tying_point[ii1], Nst) );
		RC_TRY( set_dN_matrix_shell(elem_type, Err_tying_point[ii1], dNrr) );
		RC_TRY( set_dN_matrix_shell(elem_type, Ess_tying_point[ii1], dNss) );
		RC_TRY( set_dN_matrix_shell(elem_type, Ert_tying_point[ii1], dNrt) );
		RC_TRY( set_dN_matrix_shell(elem_type, Est_tying_point[ii1], dNst) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_vect,
		                          thickness, node_location,
		                          Err_tying_point[ii1], Grr, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_vect,
		                          thickness, node_location,
		                          Ess_tying_point[ii1], Gss, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_vect,
		                          thickness, node_location,
		                          Ert_tying_point[ii1], Grt, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_vect,
		                          thickness, node_location,
		                          Est_tying_point[ii1], Gst, ws) );

		for(ii2=0; ii2<elem_node_num; ii2++){
			V1 = elem_vect[ii2].v1;
			V2 = elem_vect[ii2].v2;
			V1G1rr = inner_product3d(V1, Grr[0]);
			V2G1rr = inner_product3d(V2, Grr[0]);
			V1G2ss = inner_product3d(V1, Gss[1]);
			V2G2ss = inner_product3d(V2, Gss[1]);
			V1G1rt = inner_product3d(V1, Grt[0]);
			V1G3rt = inner_product3d(V1, Grt[2]);
			V2G1rt = inner_product3d(V2, Grt[0]);
			V2G3rt = inner_product3d(V2, Grt[2]);
			V1G2st = inner_product3d(V1, Gst[1]);
			V1G3st = inner_product3d(V1, Gst[2]);
			V2G2st = inner_product3d(V2, Gst[1]);
			V2G3st = inner_product3d(V2, Gst[2]);

			B_matrix[0][dim*ii2]   += Err_h[ii1] * (dNrr[0][ii2] * Grr[0].x);
			B_matrix[0][dim*ii2+1] += Err_h[ii1] * (dNrr[0][ii2] * Grr[0].y);
			B_matrix[0][dim*ii2+2] += Err_h[ii1] * (dNrr[0][ii2] * Grr[0].z);
			B_matrix[0][dim*ii2+3] += Err_h[ii1] * (-dNrr[0][ii2]*c1*V2G1rr);
			B_matrix[0][dim*ii2+4] += Err_h[ii1] * (dNrr[0][ii2]*c1*V1G1rr);

			B_matrix[1][dim*ii2]   += Ess_h[ii1] * (dNss[1][ii2] * Gss[1].x);
			B_matrix[1][dim*ii2+1] += Ess_h[ii1] * (dNss[1][ii2] * Gss[1].y);
			B_matrix[1][dim*ii2+2] += Ess_h[ii1] * (dNss[1][ii2] * Gss[1].z);
			B_matrix[1][dim*ii2+3] += Ess_h[ii1] * (-dNss[1][ii2]*c1*V2G2ss);
			B_matrix[1][dim*ii2+4] += Ess_h[ii1] * (dNss[1][ii2]*c1*V1G2ss);

			B_matrix[3][dim*ii2]   += Ert_h[ii1] * (dNrt[0][ii2] * Grt[2].x);
			B_matrix[3][dim*ii2+1] += Ert_h[ii1] * (dNrt[0][ii2] * Grt[2].y);
			B_matrix[3][dim*ii2+2] += Ert_h[ii1] * (dNrt[0][ii2] * Grt[2].z);
			B_matrix[3][dim*ii2+3] += Ert_h[ii1] * (-c2*Nrt[ii2]*V2G1rt
			                                        -c1*dNrt[0][ii2]*V2G3rt);
			B_matrix[3][dim*ii2+4] += Ert_h[ii1] * (c2*Nrt[ii2]*V1G1rt
			                                        + c1*dNrt[0][ii2]*V1G3rt);

			B_matrix[4][dim*ii2]   += Est_h[ii1] * (dNst[1][ii2] * Gst[2].x);
			B_matrix[4][dim*ii2+1] += Est_h[ii1] * (dNst[1][ii2] * Gst[2].y);
			B_matrix[4][dim*ii2+2] += Est_h[ii1] * (dNst[1][ii2] * Gst[2].z);
			B_matrix[4][dim*ii2+3] += Est_h[ii1] * (-c2*Nst[ii2]*V2G2st
			                                        -c1*dNst[1][ii2]*V2G3st);
			B_matrix[4][dim*ii2+4] += Est_h[ii1] * (c2*Nst[ii2]*V1G2st
			                                        + c1*dNst[1][ii2]*V1G3st);
		}
	}

	return(NORMAL_RC);
}


RC
set_tying_data_mitc3 (const double *local_point, int *ErtEst_tying_num,
                      double **Ert_rt_tying_point, double **Ert_st_tying_point,
                      double **Est_rt_tying_point, double **Est_st_tying_point,
                      double *Ert_rt_h, double *Ert_st_h,
                      double *Est_rt_h, double *Est_st_h)
{
	double r = local_point[0];
	double s = local_point[1];
	double t = local_point[2];


	*ErtEst_tying_num = 2;

	Ert_rt_tying_point[0][0] = Est_rt_tying_point[0][0] = 0.5;
	Ert_rt_tying_point[0][1] = Est_rt_tying_point[0][1] = 0.0;
	Ert_rt_tying_point[0][2] = Est_rt_tying_point[0][2] = t;
	Ert_rt_tying_point[1][0] = Est_rt_tying_point[1][0] = 0.5;
	Ert_rt_tying_point[1][1] = Est_rt_tying_point[1][1] = 0.5;
	Ert_rt_tying_point[1][2] = Est_rt_tying_point[1][2] = t;
	Ert_st_tying_point[0][0] = Est_st_tying_point[0][0] = 0.0;
	Ert_st_tying_point[0][1] = Est_st_tying_point[0][1] = 0.5;
	Ert_st_tying_point[0][2] = Est_st_tying_point[0][2] = t;
	Ert_st_tying_point[1][0] = Est_st_tying_point[1][0] = 0.5;
	Ert_st_tying_point[1][1] = Est_st_tying_point[1][1] = 0.5;
	Ert_st_tying_point[1][2] = Est_st_tying_point[1][2] = t;

	Ert_rt_h[0] = 1.0 - s;
	Ert_rt_h[1] = s;
	Ert_st_h[0] = s;
	Ert_st_h[1] = -s;
	Est_rt_h[0] = r;
	Est_rt_h[1] = -r;
	Est_st_h[0] = 1.0 - r;
	Est_st_h[1] = r;

	return(NORMAL_RC);
}


RC
set_tying_data_mitc6 (const double *local_point, int *ErrEss_tying_num,
                      int *Ers_tying_num, int *ErtEst_tying_num,
                      double **Err_tying_point, double **Ess_tying_point,
                      double **Ers_rr_tying_point, double **Ers_ss_tying_point,
                      double **Ers_rs_tying_point,
                      double **Ert_rt_tying_point, double **Ert_st_tying_point,
                      double **Est_rt_tying_point, double **Est_st_tying_point,
                      double *Err_h, double *Ess_h,
                      double *Ers_rr_h, double *Ers_ss_h, double *Ers_rs_h,
                      double *Ert_rt_h, double *Ert_st_h,
                      double *Est_rt_h, double *Est_st_h)
{
	double r = local_point[0];
	double s = local_point[1];
	double t = local_point[2];
	double r1, r2, r3;
	double s1, s2, s3;
	double root3 = sqrt(3.0);
	double root3p3 = root3 + 3.0;
	double mroot3p3 = -root3 + 3.0;
	double root3p1 = root3 + 1.0;
	double mroot3p1 = -root3 + 1.0;
	double root3p2 = root3 + 2.0;
	double mroot3p2 = -root3 + 2.0;
	double mrmsp1 = 1.0 - r - s;


	*ErrEss_tying_num = 3;
	/* Err, Ess */
	r1 = s1 = 0.5 - (1.0 / (2.0*root3));
	r2 = s2 = 0.5 + (1.0 / (2.0*root3));
	r3 = s3 = 1.0 / root3;
	Err_tying_point[0][0] = r1;
	Err_tying_point[0][1] = 0.0;
	Err_tying_point[0][2] = t;
	Err_tying_point[1][0] = r2;
	Err_tying_point[1][1] = 0.0;
	Err_tying_point[1][2] = t;
	Err_tying_point[2][0] = r1;
	Err_tying_point[2][1] = s3;
	Err_tying_point[2][2] = t;
	Ess_tying_point[0][0] = 0.0;
	Ess_tying_point[0][1] = s1;
	Ess_tying_point[0][2] = t;
	Ess_tying_point[1][0] = 0.0;
	Ess_tying_point[1][1] = s2;
	Ess_tying_point[1][2] = t;
	Ess_tying_point[2][0] = r3;
	Ess_tying_point[2][1] = s1;
	Ess_tying_point[2][2] = t;

	Err_h[0] = 0.5*(-root3p3*s + root3p1) + 3.0*r1*s - root3*r;
	Err_h[1] = 0.5*(mroot3p3*s + mroot3p1) - 3.0*r1*s + root3*r;
	Err_h[2] = root3*s;
	Ess_h[0] = 0.5*(-root3p3*r + root3p1) + 3.0*s1*r - root3*s;
	Ess_h[1] = 0.5*(mroot3p3*r + mroot3p1) - 3.0*s1*r + root3*s;
	Ess_h[2] = root3*r;

	*Ers_tying_num = 3;
	/* Ers */
	Ers_rr_tying_point[0][0] = Ers_ss_tying_point[0][0]
	                         = Ers_rs_tying_point[0][0] = r2;
	Ers_rr_tying_point[0][1] = Ers_ss_tying_point[0][1]
	                         = Ers_rs_tying_point[0][1] = s1;
	Ers_rr_tying_point[0][2] = Ers_ss_tying_point[0][2]
	                         = Ers_rs_tying_point[0][2] = t;
	Ers_rr_tying_point[1][0] = Ers_ss_tying_point[1][0]
	                         = Ers_rs_tying_point[1][0] = r1;
	Ers_rr_tying_point[1][1] = Ers_ss_tying_point[1][1]
	                         = Ers_rs_tying_point[1][1] = s2;
	Ers_rr_tying_point[1][2] = Ers_ss_tying_point[1][2]
	                         = Ers_rs_tying_point[1][2] = t;
	Ers_rr_tying_point[2][0] = Ers_ss_tying_point[2][0]
	                         = Ers_rs_tying_point[2][0] = r1;
	Ers_rr_tying_point[2][1] = Ers_ss_tying_point[2][1]
	                         = Ers_rs_tying_point[2][1] = s1;
	Ers_rr_tying_point[2][2] = Ers_ss_tying_point[2][2]
	                         = Ers_rs_tying_point[2][2] = t;

	Ers_rr_h[0] = Ers_ss_h[0] = -0.5*(mroot3p1 + mroot3p3*mrmsp1) - root3*r
	                          + 3.0*r1*mrmsp1;
	Ers_rr_h[1] = Ers_ss_h[1] = 0.5*(-root3p1 + root3p3*mrmsp1) + root3*r
	                          - 3.0*r1*mrmsp1;
	Ers_rr_h[2] = Ers_ss_h[2] = -root3 * mrmsp1;
	Ers_rs_h[0] = -Ers_rr_h[0];
	Ers_rs_h[1] = -Ers_rr_h[1];
	Ers_rs_h[2] = -Ers_rr_h[2];

	*ErtEst_tying_num = 5;
	/* Ert, Est */
	r1 = s1 = 0.5 - (1.0/(2.0*root3));
	r2 = s2 = 0.5 + (1.0/(2.0*root3));
	r3 = s3 = 1.0 / 3.0;
	Ert_rt_tying_point[0][0] = Est_rt_tying_point[0][0] = r1;
	Ert_rt_tying_point[0][1] = Est_rt_tying_point[0][1] = 0.0;
	Ert_rt_tying_point[0][2] = Est_rt_tying_point[0][2] = t;
	Ert_rt_tying_point[1][0] = Est_rt_tying_point[1][0] = r2;
	Ert_rt_tying_point[1][1] = Est_rt_tying_point[1][1] = s1;
	Ert_rt_tying_point[1][2] = Est_rt_tying_point[1][2] = t;
	Ert_rt_tying_point[2][0] = Est_rt_tying_point[2][0] = r2;
	Ert_rt_tying_point[2][1] = Est_rt_tying_point[2][1] = 0.0;
	Ert_rt_tying_point[2][2] = Est_rt_tying_point[2][2] = t;
	Ert_rt_tying_point[3][0] = Est_rt_tying_point[3][0] = r1;
	Ert_rt_tying_point[3][1] = Est_rt_tying_point[3][1] = s2;
	Ert_rt_tying_point[3][2] = Est_rt_tying_point[3][2] = t;
	Ert_rt_tying_point[4][0] = Est_rt_tying_point[4][0] = r3;
	Ert_rt_tying_point[4][1] = Est_rt_tying_point[4][1] = s3;
	Ert_rt_tying_point[4][2] = Est_rt_tying_point[4][2] = t;
	Ert_st_tying_point[0][0] = Est_st_tying_point[0][0] = 0.0;
	Ert_st_tying_point[0][1] = Est_st_tying_point[0][1] = s1;
	Ert_st_tying_point[0][2] = Est_st_tying_point[0][2] = t;
	Ert_st_tying_point[1][0] = Est_st_tying_point[1][0] = r2;
	Ert_st_tying_point[1][1] = Est_st_tying_point[1][1] = s1;
	Ert_st_tying_point[1][2] = Est_st_tying_point[1][2] = t;
	Ert_st_tying_point[2][0] = Est_st_tying_point[2][0] = 0.0;
	Ert_st_tying_point[2][1] = Est_st_tying_point[2][1] = s2;
	Ert_st_tying_point[2][2] = Est_st_tying_point[2][2] = t;
	Ert_st_tying_point[3][0] = Est_st_tying_point[3][0] = r1;
	Ert_st_tying_point[3][1] = Est_st_tying_point[3][1] = s2;
	Ert_st_tying_point[3][2] = Est_st_tying_point[3][2] = t;
	Ert_st_tying_point[4][0] = Est_st_tying_point[4][0] = r3;
	Ert_st_tying_point[4][1] = Est_st_tying_point[4][1] = s3;
	Ert_st_tying_point[4][2] = Est_st_tying_point[4][2] = t;

	Ert_rt_h[0] = 0.5*(root3p1 + root3p3*s*s) - root3*r - root3p2*s
	            + root3*r*s;
	Ert_rt_h[1] = -s + 0.5*(root3p3*r*s + mroot3p3*s*s);
	Ert_rt_h[2] = 0.5*(mroot3p1 + mroot3p3*s*s) + root3*r - mroot3p2*s
	            - root3*r*s;
	Ert_rt_h[3] = -s + 0.5*(mroot3p3*r*s + root3p3*s*s);
	Ert_rt_h[4] = 6.0*s - 3.0*r*s - 6.0*s*s;

	Ert_st_h[0] = 0.5*(root3p1*s - root3p3*r*s) - root3*s*s;
	Ert_st_h[1] = s - 0.5*(root3p3*r*s + mroot3p3*s*s);
	Ert_st_h[2] = 0.5*(mroot3p1*s - mroot3p3*r*s) + root3*s*s;
	Ert_st_h[3] = s - 0.5*(mroot3p3*r*s + root3p3*s*s);
	Ert_st_h[4] = -3.0*s + 6.0*r*s + 3.0*s*s;

	Est_rt_h[0] = 0.5*(root3p1*r - root3p3*r*s) - root3*r*r;
	Est_rt_h[1] = r - 0.5*(root3p3*r*r + mroot3p3*r*s);
	Est_rt_h[2] = 0.5*(mroot3p1*r - mroot3p3*r*s) + root3*r*r;
	Est_rt_h[3] = r - 0.5*(mroot3p3*r*r + root3p3*r*s);
	Est_rt_h[4] = -3.0*r + 3.0*r*r + 6.0*r*s;

	Est_st_h[0] = 0.5*(root3p1 + root3p3*r*r) - root3p2*r - root3*s
	            + root3*r*s;
	Est_st_h[1] = -r + 0.5*(root3p3*r*r + mroot3p3*r*s);
	Est_st_h[2] = 0.5*(mroot3p1 + mroot3p3*r*r) - mroot3p2*r + root3*s
	            - root3*r*s;
	Est_st_h[3] = -r + 0.5*(mroot3p3*r*r + root3p3*r*s);
	Est_st_h[4] = 6.0*r - 6.0*r*r - 3.0*r*s;

	return(NORMAL_RC);
}


RC
set_tying_data_mitc4 (const double *local_point, int *ErtEst_tying_num,
                      double **Ert_tying_point, double **Est_tying_point,
                      double *Ert_h, double *Est_h)
{
	double r = local_point[0];
	double s = local_point[1];
	double t = local_point[2];


	*ErtEst_tying_num = 2;

	Ert_tying_point[0][0] = 0.0;
	Ert_tying_point[0][1] = 1.0;
	Ert_tying_point[0][2] = t;
	Ert_tying_point[1][0] = 0.0;
	Ert_tying_point[1][1] = -1.0;
	Ert_tying_point[1][2] = t;
	Ert_h[0] = 0.5 * (1.0 + s);
	Ert_h[1] = 0.5 * (1.0 - s);

	Est_tying_point[0][0] = 1.0;
	Est_tying_point[0][1] = 0.0;
	Est_tying_point[0][2] = t;
	Est_tying_point[1][0] = -1.0;
	Est_tying_point[1][1] = 0.0;
	Est_tying_point[1][2] = t;
	Est_h[0] = 0.5 * (1.0 + r);
	Est_h[1] = 0.5 * (1.0 - r);

	return(NORMAL_RC);
}


RC
set_tying_data_mitc8 (const double *local_point, int *Ers_tying_num,
                      int *ErrEssErtEst_tying_num,
                      double **Err_tying_point, double **Ess_tying_point,
                      double **Ers_tying_point, double **Ert_tying_point,
                      double **Est_tying_point,
                      double *Err_h, double *Ess_h, double *Ers_h,
                      double *Ert_h, double *Est_h)
{
	double r = local_point[0];
	double s = local_point[1];
	double t = local_point[2];
	double r1, r2, r3, r4, r5, r6;
	double s1, s2, s3, s4, s5, s6;


	*Ers_tying_num = 4;
	/* Ers */
	Ers_tying_point[0][0] = r1 = -0.577350296189626;
	Ers_tying_point[0][1] = s1 = -0.577350296189626;
	Ers_tying_point[0][2] = t;
	Ers_tying_point[1][0] = r2 = +0.577350296189626;
	Ers_tying_point[1][1] = s2 = -0.577350296189626;
	Ers_tying_point[1][2] = t;
	Ers_tying_point[2][0] = r3 = +0.577350296189626;
	Ers_tying_point[2][1] = s3 = +0.577350296189626;
	Ers_tying_point[2][2] = t;
	Ers_tying_point[3][0] = r4 = -0.577350296189626;
	Ers_tying_point[3][1] = s4 = +0.577350296189626;
	Ers_tying_point[3][2] = t;
	Ers_h[0] = ((r - r2)/(r1 - r2))*((s - s4)/(s1 - s4));
	Ers_h[1] = ((r - r1)/(r2 - r1))*((s - s3)/(s2 - s3));
	Ers_h[2] = ((r - r4)/(r3 - r4))*((s - s2)/(s3 - s2));
	Ers_h[3] = ((r - r3)/(r4 - r3))*((s - s1)/(s4 - s1));


	*ErrEssErtEst_tying_num = 6;
	/* Err, Ert */
	Err_tying_point[0][0] = Ert_tying_point[0][0] = r1 = -0.577350296189626;
	Err_tying_point[0][1] = Ert_tying_point[0][1] = s1 = -0.774596669241483;
	Err_tying_point[0][2] = Ert_tying_point[0][2] = t;
	Err_tying_point[1][0] = Ert_tying_point[1][0] = r2 = +0.577350296189626;
	Err_tying_point[1][1] = Ert_tying_point[1][1] = s2 = -0.774596669241483;
	Err_tying_point[1][2] = Ert_tying_point[1][2] = t;
	Err_tying_point[2][0] = Ert_tying_point[2][0] = r3 = +0.577350296189626;
	Err_tying_point[2][1] = Ert_tying_point[2][1] = s3 = 0.0;
	Err_tying_point[2][2] = Ert_tying_point[2][2] = t;
	Err_tying_point[3][0] = Ert_tying_point[3][0] = r4 = +0.577350296189626;
	Err_tying_point[3][1] = Ert_tying_point[3][1] = s4 = +0.774596669241483;
	Err_tying_point[3][2] = Ert_tying_point[3][2] = t;
	Err_tying_point[4][0] = Ert_tying_point[4][0] = r5 = -0.577350296189626;
	Err_tying_point[4][1] = Ert_tying_point[4][1] = s5 = +0.774596669241483;
	Err_tying_point[4][2] = Ert_tying_point[4][2] = t;
	Err_tying_point[5][0] = Ert_tying_point[5][0] = r6 = -0.577350296189626;
	Err_tying_point[5][1] = Ert_tying_point[5][1] = s6 = 0.0;
	Err_tying_point[5][2] = Ert_tying_point[5][2] = t;
	Err_h[0] = Ert_h[0] = ((r - r2) / (r1 - r2)) * (((s - s6)*(s - s5))
	                    / ((s1 - s6)*(s1 - s5)));
	Err_h[1] = Ert_h[1] = ((r - r1) / (r2 - r1)) * (((s - s3)*(s - s4))
	                    / ((s2 - s3) * (s2 - s4)));
	Err_h[2] = Ert_h[2] = ((r - r6) / (r3 - r6)) * (((s - s2)*(s - s4))
	                    / ((s3 - s2) * (s3 - s4)));
	Err_h[3] = Ert_h[3] = ((r - r5) / (r4 - r5)) * (((s - s3)*(s - s2))
	                    / ((s4 - s3) * (s4 - s2)));
	Err_h[4] = Ert_h[4] = ((r - r4) / (r5 - r4)) * (((s - s6)*(s - s1))
	                    / ((s5- s6) * (s5 - s1)));
	Err_h[5] = Ert_h[5] = ((r - r3) / (r6 - r3)) * (((s - s5)*(s - s1))
	                    / ((s6 - s5) * (s6 - s1)));

	/* Ess, Est */
	Ess_tying_point[0][0] = Est_tying_point[0][0] = r1 = -0.774596669241483;
	Ess_tying_point[0][1] = Est_tying_point[0][1] = s1 = -0.577350296189626;
	Ess_tying_point[0][2] = Est_tying_point[0][2] = t;
	Ess_tying_point[1][0] = Est_tying_point[1][0] = r2 = 0.0;
	Ess_tying_point[1][1] = Est_tying_point[1][1] = s2 = -0.577350296189626;
	Ess_tying_point[1][2] = Est_tying_point[1][2] = t;
	Ess_tying_point[2][0] = Est_tying_point[2][0] = r3 = +0.774596669241483;
	Ess_tying_point[2][1] = Est_tying_point[2][1] = s3 = -0.577350296189626;
	Ess_tying_point[2][2] = Est_tying_point[2][2] = t;
	Ess_tying_point[3][0] = Est_tying_point[3][0] = r4 = +0.774596669241483;
	Ess_tying_point[3][1] = Est_tying_point[3][1] = s4 = +0.577350296189626;
	Ess_tying_point[3][2] = Est_tying_point[3][2] = t;
	Ess_tying_point[4][0] = Est_tying_point[4][0] = r5 = 0.0;
	Ess_tying_point[4][1] = Est_tying_point[4][1] = s5 = +0.577350296189626;
	Ess_tying_point[4][2] = Est_tying_point[4][2] = t;
	Ess_tying_point[5][0] = Est_tying_point[5][0] = r6 = -0.774596669241483;
	Ess_tying_point[5][1] = Est_tying_point[5][1] = s6 = +0.577350296189626;
	Ess_tying_point[5][2] = Est_tying_point[5][2] = t;
	Ess_h[0] = Est_h[0] = (((r - r2)*(r - r3))/((r1 - r2)*(r1 - r3)))
	                    * ((s - s6)/(s1 - s6));
	Ess_h[1] = Est_h[1] = (((r - r1)*(r - r3))/((r2 - r1)*(r2 - r3)))
	                    * ((s - s5)/(s2 - s5));
	Ess_h[2] = Est_h[2] = (((r - r1)*(r - r2))/((r3 - r1)*(r3 - r2)))
	                    * ((s - s4)/(s3 - s4));
	Ess_h[3] = Est_h[3] = (((r - r5)*(r - r6))/((r4 - r5)*(r4 - r6)))
	                    * ((s - s3)/(s4 - s3));
	Ess_h[4] = Est_h[4] = (((r - r4)*(r - r6))/((r5 - r4)*(r5 - r6)))
	                    * ((s - s2)/(s5 - s2));
	Ess_h[5] = Est_h[5] = (((r - r5)*(r - r4))/((r6 - r5)*(r6 - r4)))
	                    * ((s - s1)/(s6 - s1));

	return(NORMAL_RC);
}


RC
set_tying_data_mitc9 (const double *local_point, int *Ers_tying_num,
                      int *ErrEssErtEst_tying_num,
                      double **Err_tying_point, double **Ess_tying_point,
                      double **Ers_tying_point, double **Ert_tying_point,
                      double **Est_tying_point,
                      double *Err_h, double *Ess_h, double *Ers_h,
                      double *Ert_h, double *Est_h)
{
	double r = local_point[0];
	double s = local_point[1];
	double t = local_point[2];
	double r1, r2, r3, r4, r5, r6;
	double s1, s2, s3, s4, s5, s6;


	*Ers_tying_num = 4;
	/* Ers */
	Ers_tying_point[0][0] = r1 = -0.577350296189626;
	Ers_tying_point[0][1] = s1 = -0.577350296189626;
	Ers_tying_point[0][2] = t;
	Ers_tying_point[1][0] = r2 = +0.577350296189626;
	Ers_tying_point[1][1] = s2 = -0.577350296189626;
	Ers_tying_point[1][2] = t;
	Ers_tying_point[2][0] = r3 = +0.577350296189626;
	Ers_tying_point[2][1] = s3 = +0.577350296189626;
	Ers_tying_point[2][2] = t;
	Ers_tying_point[3][0] = r4 = -0.577350296189626;
	Ers_tying_point[3][1] = s4 = +0.577350296189626;
	Ers_tying_point[3][2] = t;
	Ers_h[0] = ((r - r2)/(r1 - r2))*((s - s4)/(s1 - s4));
	Ers_h[1] = ((r - r1)/(r2 - r1))*((s - s3)/(s2 - s3));
	Ers_h[2] = ((r - r4)/(r3 - r4))*((s - s2)/(s3 - s2));
	Ers_h[3] = ((r - r3)/(r4 - r3))*((s - s1)/(s4 - s1));


	*ErrEssErtEst_tying_num = 6;
	/* Err, Ert */
	Err_tying_point[0][0] = Ert_tying_point[0][0] = r1 = -0.577350296189626;
	Err_tying_point[0][1] = Ert_tying_point[0][1] = s1 = -0.774596669241483;
	Err_tying_point[0][2] = Ert_tying_point[0][2] = t;
	Err_tying_point[1][0] = Ert_tying_point[1][0] = r2 = +0.577350296189626;
	Err_tying_point[1][1] = Ert_tying_point[1][1] = s2 = -0.774596669241483;
	Err_tying_point[1][2] = Ert_tying_point[1][2] = t;
	Err_tying_point[2][0] = Ert_tying_point[2][0] = r3 = +0.577350296189626;
	Err_tying_point[2][1] = Ert_tying_point[2][1] = s3 = 0.0;
	Err_tying_point[2][2] = Ert_tying_point[2][2] = t;
	Err_tying_point[3][0] = Ert_tying_point[3][0] = r4 = +0.577350296189626;
	Err_tying_point[3][1] = Ert_tying_point[3][1] = s4 = +0.774596669241483;
	Err_tying_point[3][2] = Ert_tying_point[3][2] = t;
	Err_tying_point[4][0] = Ert_tying_point[4][0] = r5 = -0.577350296189626;
	Err_tying_point[4][1] = Ert_tying_point[4][1] = s5 = +0.774596669241483;
	Err_tying_point[4][2] = Ert_tying_point[4][2] = t;
	Err_tying_point[5][0] = Ert_tying_point[5][0] = r6 = -0.577350296189626;
	Err_tying_point[5][1] = Ert_tying_point[5][1] = s6 = 0.0;
	Err_tying_point[5][2] = Ert_tying_point[5][2] = t;
	Err_h[0] = Ert_h[0] = ((r - r2) / (r1 - r2)) * (((s - s6)*(s - s5))
	                    / ((s1 - s6)*(s1 - s5)));
	Err_h[1] = Ert_h[1] = ((r - r1) / (r2 - r1)) * (((s - s3)*(s - s4))
	                    / ((s2 - s3) * (s2 - s4)));
	Err_h[2] = Ert_h[2] = ((r - r6) / (r3 - r6)) * (((s - s2)*(s - s4))
	                    / ((s3 - s2) * (s3 - s4)));
	Err_h[3] = Ert_h[3] = ((r - r5) / (r4 - r5)) * (((s - s3)*(s - s2))
	                    / ((s4 - s3) * (s4 - s2)));
	Err_h[4] = Ert_h[4] = ((r - r4) / (r5 - r4)) * (((s - s6)*(s - s1))
	                    / ((s5- s6) * (s5 - s1)));
	Err_h[5] = Ert_h[5] = ((r - r3) / (r6 - r3)) * (((s - s5)*(s - s1))
	                    / ((s6 - s5) * (s6 - s1)));

	/* Ess, Est */
	Ess_tying_point[0][0] = Est_tying_point[0][0] = r1 = -0.774596669241483;
	Ess_tying_point[0][1] = Est_tying_point[0][1] = s1 = -0.577350296189626;
	Ess_tying_point[0][2] = Est_tying_point[0][2] = t;
	Ess_tying_point[1][0] = Est_tying_point[1][0] = r2 = 0.0;
	Ess_tying_point[1][1] = Est_tying_point[1][1] = s2 = -0.577350296189626;
	Ess_tying_point[1][2] = Est_tying_point[1][2] = t;
	Ess_tying_point[2][0] = Est_tying_point[2][0] = r3 = +0.774596669241483;
	Ess_tying_point[2][1] = Est_tying_point[2][1] = s3 = -0.577350296189626;
	Ess_tying_point[2][2] = Est_tying_point[2][2] = t;
	Ess_tying_point[3][0] = Est_tying_point[3][0] = r4 = +0.774596669241483;
	Ess_tying_point[3][1] = Est_tying_point[3][1] = s4 = +0.577350296189626;
	Ess_tying_point[3][2] = Est_tying_point[3][2] = t;
	Ess_tying_point[4][0] = Est_tying_point[4][0] = r5 = 0.0;
	Ess_tying_point[4][1] = Est_tying_point[4][1] = s5 = +0.577350296189626;
	Ess_tying_point[4][2] = Est_tying_point[4][2] = t;
	Ess_tying_point[5][0] = Est_tying_point[5][0] = r6 = -0.774596669241483;
	Ess_tying_point[5][1] = Est_tying_point[5][1] = s6 = +0.577350296189626;
	Ess_tying_point[5][2] = Est_tying_point[5][2] = t;
	Ess_h[0] = Est_h[0] = (((r - r2)*(r - r3))/((r1 - r2)*(r1 - r3)))
	                    * ((s - s6)/(s1 - s6));
	Ess_h[1] = Est_h[1] = (((r - r1)*(r - r3))/((r2 - r1)*(r2 - r3)))
	                    * ((s - s5)/(s2 - s5));
	Ess_h[2] = Est_h[2] = (((r - r1)*(r - r2))/((r3 - r1)*(r3 - r2)))
	                    * ((s - s4)/(s3 - s4));
	Ess_h[3] = Est_h[3] = (((r - r5)*(r - r6))/((r4 - r5)*(r4 - r6)))
	                    * ((s - s3)/(s4 - s3));
	Ess_h[4] = Est_h[4] = (((r - r4)*(r - r6))/((r5 - r4)*(r5 - r6)))
	                    * ((s - s2)/(s5 - s2));
	Ess_h[5] = Est_h[5] = (((r - r5)*(r - r4))/((r6 - r5)*(r6 - r4)))
	                    * ((s - s1)/(s6 - s1));

	return(NORMAL_RC);
}


/* 要素質量マトリクス */
RC
make_shell_mass_matrix (ELEMENT element, NODE_ARRAY node,
                        MATERIAL_PROP_ARRAY material,
                        PHYSICAL_PROP_ARRAY physical,
                        DIRECTOR_VECT_ARRAY d_vect,
                        ELEM_MATRIX *elem_matrix,
                        WORK_SET ws1, WORK_SET ws2)
{
	int ii1, ii2, ii3, ii4, ii5;
	int dim;
	int elem_node_num;
	int elem_matrix_size;
	int g_num_point;
	int mat_index;
	int prop_index;
	double rho;
	double thickness;
	double det_J;
	double *g_weights = ws1.vec_G[0];
	double **node_location = ws1.mat_N_3[0];
	double **g_points = ws1.mat_G_4[0];
	double **N_matrix = ws1.mat_6_3N[0];
	double **NtN = ws1.mat_3N_3N[0];
	VECT3D G0[3];
	DIRECTOR_VECT elem_vect[MAX_SHELL_NODE];
	ELEM_TYPE elem_type;


	if(element.label < 0) return(ARG_ERROR_RC);

	prop_index = search_physical_prop_label(physical, element.physical);
	if(prop_index < 0) return(SEEK_ERROR_RC);
	if(NST_PROP_PSHELL != physical.array[prop_index].i_info[1]){
		return(ARG_ERROR_RC);
	}

	if( (dim = element_dim_shell(element.type)) != 6){
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

	elem_type = element.type;
	elem_node_num = element.node_num;
	elem_matrix_size = elem_node_num * dim;
	RC_TRY( allocate2D(elem_matrix_size, elem_matrix_size,
	                   &(elem_matrix->matrix)) );

	RC_TRY( set_element_direcor_vect(element, d_vect, elem_vect) );
	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( Gauss_const_default_shell(elem_type, g_points, g_weights,
	                                  &g_num_point) );
	RC_TRY( get_thickness(element, physical, &thickness) );

	/* density */
	mat_index = search_material_prop_element(material, physical, &element);
	if(mat_index < 0) return(SEEK_ERROR_RC);

	rho = material.array[mat_index].rho;

	for(ii1=0; ii1<g_num_point; ii1++){
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num,
		                          elem_vect, thickness,
		                          node_location, g_points[ii1],
		                          G0, ws2) );
		RC_TRY( set_disp_N_matrix_shell(elem_type, elem_node_num, elem_vect,
		                                thickness, g_points[ii1],
		                                N_matrix, ws2) );
		mul_matrix_AtB(3, 5*elem_node_num, 5*elem_node_num, N_matrix,
		               N_matrix, NtN);

		det_J = scalar_triple_product3d(G0[0], G0[1], G0[2]);

		for(ii2=0; ii2<elem_node_num; ii2++){
			for(ii3=0; ii3<elem_node_num; ii3++){
				for(ii4=0; ii4<5; ii4++){
					for(ii5=0; ii5<5; ii5++){
						elem_matrix->matrix[dim*ii2+ii4][dim*ii3+ii5]
						+= rho*NtN[5*ii2+ii4][5*ii3+ii5]*det_J*g_weights[ii1];
					}
				}
			}
		}
	}

	/* 回転自由度を全体座標系成分への変換 */
	RC_TRY( trans_global_rotation_dof(elem_matrix_size, elem_matrix->matrix,
	                                  elem_vect) );

	return(NORMAL_RC);
}


static RC
set_disp_N_matrix_shell (ELEM_TYPE elem_type, int elem_node_num,
                         DIRECTOR_VECT elem_vect[], double thickness,
                         const double *local_point,
                         double **N_matrix, WORK_SET ws)
{
	int ii1;
	double a = thickness;
	double t = local_point[2];
	double c = a * t * 0.5;
	double *N = ws.vec_N[0];


	RC_TRY( set_N_vector_shell(elem_type, local_point, N) );
	for(ii1=0; ii1<elem_node_num; ii1++){
		VECT3D V1 = elem_vect[ii1].v1;
		VECT3D V2 = elem_vect[ii1].v2;

		N_matrix[0][ii1*5+0] = N[ii1];
		N_matrix[0][ii1*5+1] = 0.0;
		N_matrix[0][ii1*5+2] = 0.0;
		N_matrix[0][ii1*5+3] = -c*N[ii1]*V2.x;
		N_matrix[0][ii1*5+4] = c*N[ii1]*V1.x;

		N_matrix[1][ii1*5+0] = 0.0;
		N_matrix[1][ii1*5+1] = N[ii1];
		N_matrix[1][ii1*5+2] = 0.0;
		N_matrix[1][ii1*5+3] = -c*N[ii1]*V2.y;
		N_matrix[1][ii1*5+4] = c*N[ii1]*V1.y;

		N_matrix[2][ii1*5+0] = 0.0;
		N_matrix[2][ii1*5+1] = 0.0;
		N_matrix[2][ii1*5+2] = N[ii1];
		N_matrix[2][ii1*5+3] = -c*N[ii1]*V2.z;
		N_matrix[2][ii1*5+4] = c*N[ii1]*V1.z;
	}

	return(NORMAL_RC);
}


/* シェル要素のガウスポイント設定 */
RC
Gauss_const_default_shell (ELEM_TYPE type, double **points,
                           double *weight, int *num_point)
{
	switch(type){
	case ELEM_TRI1:
		*num_point = 6;
		RC_TRY( Gauss_const_mitc_tri(*num_point, points, weight) );
		break;
	case ELEM_TRI2:
		*num_point = 14;
		RC_TRY( Gauss_const_mitc_tri(*num_point, points, weight) );
		break;
	case ELEM_QUAD1:
		*num_point = 8;
		RC_TRY( Gauss_const_mitc_quad(*num_point, points, weight) );
		break;
	case ELEM_QUAD2:
	case ELEM_QUAD2L:
		*num_point = 18;
		RC_TRY( Gauss_const_mitc_quad(*num_point, points, weight) );
		break;
	default :
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
Gauss_const_mitc_tri (int num_point, double **rst, double *weight)
{
	double a, b;
	double w1, w2, w3;

	switch(num_point){
	case 6:
		rst[0][0] = 0.166666666666667;
		rst[0][1] = 0.166666666666667;
		rst[0][2] = -0.577350269189626;
		weight[0] = 0.333333333333333*0.5;
		rst[1][0] = 0.166666666666667;
		rst[1][1] = 0.166666666666667;
		rst[1][2] = +0.577350269189626;
		weight[1] = 0.333333333333333*0.5;
		rst[2][0] = 0.666666666666667;
		rst[2][1] = 0.666666666666667;
		rst[2][2] = -0.577350269189626;
		weight[2] = 0.333333333333333*0.5;
		rst[3][0] = 0.666666666666667;
		rst[3][1] = 0.666666666666667;
		rst[3][2] = +0.577350269189626;
		weight[3] = 0.333333333333333*0.5;
		rst[4][0] = 0.166666666666667;
		rst[4][1] = 0.666666666666667;
		rst[4][2] = -0.577350269189626;
		weight[4] = 0.333333333333333*0.5;
		rst[5][0] = 0.166666666666667;
		rst[5][1] = 0.666666666666667;
		rst[5][2] = +0.577350269189626;
		weight[5] = 0.333333333333333*0.5;
		break;
	case 14:
		a = (6.0+sqrt(15.0))/21.0;
		b = 4.0/7.0 - a;
		w1 = 9.0/80.0;
		w2 = (155.0 + sqrt(15.0))/2400.0;
		w3 = 31.0/240.0 - w2;
		rst[0][0] = b;
		rst[0][1] = b;
		rst[0][2] = -0.577350269189626;
		weight[0] = w3;
		rst[1][0] = b;
		rst[1][1] = b;
		rst[1][2] = +0.577350269189626;
		weight[1] = w3;
		rst[2][0] = 1.0 - 2.0*b;
		rst[2][1] = b;
		rst[2][2] = -0.577350269189626;
		weight[2] = w3;
		rst[3][0] = 1.0 - 2.0*b;
		rst[3][1] = b;
		rst[3][2] = +0.577350269189626;
		weight[3] = w3;
		rst[4][0] = b;
		rst[4][1] = 1.0 - 2.0*b;
		rst[4][2] = -0.577350269189626;
		weight[4] = w3;
		rst[5][0] = b;
		rst[5][1] = 1.0 - 2.0*b;
		rst[5][2] = +0.577350269189626;
		weight[5] = w3;
		rst[6][0] = a;
		rst[6][1] = 1.0 - 2.0*a;
		rst[6][2] = -0.577350269189626;
		weight[6] = w2;
		rst[7][0] = a;
		rst[7][1] = 1.0 - 2.0*a;
		rst[7][2] = +0.577350269189626;
		weight[7] = w2;
		rst[8][0] = a;
		rst[8][1] = a;
		rst[8][2] = -0.577350269189626;
		weight[8] = w2;
		rst[9][0] = a;
		rst[9][1] = a;
		rst[9][2] = +0.577350269189626;
		weight[9] = w2;
		rst[10][0] = 1.0 - 2.0*a;
		rst[10][1] = a;
		rst[10][2] = -0.577350269189626;
		weight[10] = w2;
		rst[11][0] = 1.0 - 2.0*a;
		rst[11][1] = a;
		rst[11][2] = +0.577350269189626;
		weight[11] = w2;
		rst[12][0] = 1.0/3.0;
		rst[12][1] = 1.0/3.0;
		rst[12][2] = -0.577350269189626;
		weight[12] = w1;
		rst[13][0] = 1.0/3.0;
		rst[13][1] = 1.0/3.0;
		rst[13][2] = +0.577350269189626;
		weight[13] = w1;
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
Gauss_const_mitc_quad (int num_point, double **rst, double *weight)
{
	int ii1, ii2, ii3;
	int inc;
	double p[3];
	double pt[2];
	double w[4];
	double wt[2];

	switch(num_point){
	case 8:
		rst[0][0] = -0.577350269189626;
		rst[0][1] = -0.577350269189626;
		rst[0][2] = -0.577350269189626;
		weight[0] = 1.000000000000000;
		rst[1][0] = +0.577350269189626;
		rst[1][1] = -0.577350269189626;
		rst[1][2] = -0.577350269189626;
		weight[1] = 1.000000000000000;
		rst[2][0] = +0.577350269189626;
		rst[2][1] = +0.577350269189626;
		rst[2][2] = -0.577350269189626;
		weight[2] = 1.000000000000000;
		rst[3][0] = -0.577350269189626;
		rst[3][1] = +0.577350269189626;
		rst[3][2] = -0.577350269189626;
		weight[3] = 1.000000000000000;
		rst[4][0] = -0.577350269189626;
		rst[4][1] = -0.577350269189626;
		rst[4][2] = +0.577350269189626;
		weight[4] = 1.000000000000000;
		rst[5][0] = +0.577350269189626;
		rst[5][1] = -0.577350269189626;
		rst[5][2] = +0.577350269189626;
		weight[5] = 1.000000000000000;
		rst[6][0] = +0.577350269189626;
		rst[6][1] = +0.577350269189626;
		rst[6][2] = +0.577350269189626;
		weight[6] = 1.000000000000000;
		rst[7][0] = -0.577350269189626;
		rst[7][1] = +0.577350269189626;
		rst[7][2] = +0.577350269189626;
		weight[7] = 1.000000000000000;
		break;
	case 18:
		p[0] = -0.774596669241483;
		p[1] = +0.774596669241483;
		p[2] =  0.000000000000000;
		pt[0] = -0.577350269189626;
		pt[1] = +0.577350269189626;
		w[0] = 0.555555555555556;
		w[1] = 0.555555555555556;
		w[2] = 2.0 - (w[0] + w[1]);
		wt[0] = 1.000000000000000;
		wt[1] = 1.000000000000000;
		inc = 0;
		for(ii1=0; ii1<3; ii1++){
			for(ii2=0; ii2<3; ii2++){
				for(ii3=0; ii3<2; ii3++){
					rst[inc][0] = p[ii1];
					rst[inc][1] = p[ii2];
					rst[inc][2] = pt[ii3];
					weight[inc] = w[ii1] * w[ii2] * wt[ii3];
					inc++;
				}
			}
		}
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


/* シェル要素の形状関数の設定 */
RC
set_N_vector_shell (ELEM_TYPE type, const double *coord, double *N)
{
	if( (coord == NULL)||(N == NULL) ) return(ARG_ERROR_RC);

	switch(type){
	case ELEM_TRI1:
		RC_TRY( N_mitc3(coord, N) );
		break;
	case ELEM_QUAD1:
		RC_TRY( N_mitc4(coord, N) );
		break;
	case ELEM_TRI2:
		RC_TRY( N_mitc6(coord, N) );
		break;
	case ELEM_QUAD2:
		RC_TRY( N_mitc8(coord, N) );
		break;
	case ELEM_QUAD2L:
		RC_TRY( N_mitc9(coord, N) );
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
N_mitc3 (const double *rst, double *N)
{
	double r = rst[0];
	double s = rst[1];


	N[0] = 1.0 - r - s;
	N[1] = r;
	N[2] = s;

	return(NORMAL_RC);
}


static RC
N_mitc4 (const double *rst, double *N)
{
	double r = rst[0];
	double s = rst[1];
	double m_r = 1.0 - r;
	double p_r = 1.0 + r;
	double m_s = 1.0 - s;
	double p_s = 1.0 + s;


	N[0] = 0.25 * m_r * m_s;
	N[1] = 0.25 * p_r * m_s;
	N[2] = 0.25 * p_r * p_s;
	N[3] = 0.25 * m_r * p_s;

	return(NORMAL_RC);
}


static RC
N_mitc6 (const double *rst, double *N)
{
	double r = rst[0];
	double s = rst[1];
	double N4 = 4.0*r*s;
	double N5 = 4.0*s*(1.0 - r - s);
	double N6 = 4.0*r*(1.0 - r - s);


	N[0] = 1.0 - r - s - 0.5*N5 - 0.5*N6;
	N[1] = r - 0.5*N6 - 0.5*N4;
	N[2] = s - 0.5*N4 - 0.5*N5;
	N[3] = N4;
	N[4] = N5;
	N[5] = N6;

	return(NORMAL_RC);
}


static RC
N_mitc8 (const double *rst, double *N)
{
	double r = rst[0];
	double s = rst[1];
	double m_r = 1.0 - r;
	double m_s = 1.0 - s;
	double m_rr = 1.0 - r*r;
	double m_ss = 1.0 - s*s;
	double p_r = 1.0 + r;
	double p_s = 1.0 + s;
	double N5 = 0.5*m_rr*m_s;
	double N6 = 0.5*m_ss*p_r;
	double N7 = 0.5*m_rr*p_s;
	double N8 = 0.5*m_ss*m_r;


	N[0] = 0.25*m_r*m_s - 0.5*N8 - 0.5*N5;
	N[1] = 0.25*p_r*m_s - 0.5*N5 - 0.5*N6;
	N[2] = 0.25*p_r*p_s - 0.5*N6 - 0.5*N7;
	N[3] = 0.25*m_r*p_s - 0.5*N7 - 0.5*N8;
	N[4] = N5;
	N[5] = N6;
	N[6] = N7;
	N[7] = N8;

	return(NORMAL_RC);
}


static RC
N_mitc9 (const double *rst, double *N)
{
	double r = rst[0];
	double s = rst[1];
	double m_r = 1.0 - r;
	double m_s = 1.0 - s;
	double m_rr = 1.0 - r*r;
	double m_ss = 1.0 - s*s;
	double p_r = 1.0 + r;
	double p_s = 1.0 + s;
	double N9 = m_rr*m_ss;
	double N5 = 0.5*m_rr*m_s - 0.5*N9;
	double N6 = 0.5*m_ss*p_r - 0.5*N9;
	double N7 = 0.5*m_rr*p_s - 0.5*N9;
	double N8 = 0.5*m_ss*m_r - 0.5*N9;


	N[0] = 0.25*m_r*m_s - 0.5*N8 - 0.5*N5 - 0.25*N9;
	N[1] = 0.25*p_r*m_s - 0.5*N5 - 0.5*N6 - 0.25*N9;
	N[2] = 0.25*p_r*p_s - 0.5*N6 - 0.5*N7 - 0.25*N9;
	N[3] = 0.25*m_r*p_s - 0.5*N7 - 0.5*N8 - 0.25*N9;
	N[4] = N5;
	N[5] = N6;
	N[6] = N7;
	N[7] = N8;
	N[8] = N9;

	return(NORMAL_RC);
}


/* シェル要素の形状関数に対する偏導関数の設定 */
RC
set_dN_matrix_shell (ELEM_TYPE type, const double *coord, double **dN)
{
	if( (coord == NULL)||(dN == NULL) ) return(ARG_ERROR_RC);

	switch(type){
	case ELEM_TRI1:
		RC_TRY( dN_mitc3(coord, dN) );
		break;
	case ELEM_QUAD1:
		RC_TRY( dN_mitc4(coord, dN) );
		break;
	case ELEM_TRI2:
		RC_TRY( dN_mitc6(coord, dN) );
		break;
	case ELEM_QUAD2:
		RC_TRY( dN_mitc8(coord, dN) );
		break;
	case ELEM_QUAD2L:
		RC_TRY( dN_mitc9(coord, dN) );
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
dN_mitc3 (const double *rst, double **dN)
{
	dN[0][0] = -1.0;
	dN[0][1] = +1.0;
	dN[0][2] =  0.0;

	dN[1][0] = -1.0;
	dN[1][1] =  0.0;
	dN[1][2] = +1.0;

	return(NORMAL_RC);
}


static RC
dN_mitc4 (const double *rst, double **dN)
{
	double r = rst[0];
	double s = rst[1];
	double m_r = 1.0 - r;
	double p_r = 1.0 + r;
	double m_s = 1.0 - s;
	double p_s = 1.0 + s;


	dN[0][0] = -0.25 * m_s;
	dN[0][1] = +0.25 * m_s;
	dN[0][2] = +0.25 * p_s;
	dN[0][3] = -0.25 * p_s;

	dN[1][0] = -0.25 * m_r;
	dN[1][1] = -0.25 * p_r;
	dN[1][2] = +0.25 * p_r;
	dN[1][3] = +0.25 * m_r;

	return(NORMAL_RC);
}


static RC
dN_mitc6 (const double *rst, double **dN)
{
	double r = rst[0];
	double s = rst[1];
	double dNr4 = 4.0*s;
	double dNr5 = - 4.0*s;
	double dNr6 = 4.0 - 8.0*r - 4.0*s;
	double dNs4 = 4.0*r;
	double dNs5 = 4.0 - 4.0*r - 8.0*s;
	double dNs6 = - 4.0*r;


	dN[0][0] = - 1.0 - 0.5*dNr5 - 0.5*dNr6;
	dN[0][1] = 1.0 - 0.5*dNr6 - 0.5*dNr4;
	dN[0][2] = - 0.5*dNr4 - 0.5*dNr5;
	dN[0][3] = dNr4;
	dN[0][4] = dNr5;
	dN[0][5] = dNr6;

	dN[1][0] = - 1.0 - 0.5*dNs5 - 0.5*dNs6;
	dN[1][1] = - 0.5*dNs6 - 0.5*dNs4;
	dN[1][2] = 1.0 - 0.5*dNs4 - 0.5*dNs5;
	dN[1][3] = dNs4;
	dN[1][4] = dNs5;
	dN[1][5] = dNs6;

	return(NORMAL_RC);
}


static RC
dN_mitc8 (const double *rst, double **dN)
{
	double r = rst[0];
	double s = rst[1];
	double m_r = 1.0 - r;
	double m_s = 1.0 - s;
	double m_rr = 1.0 - r*r;
	double m_ss = 1.0 - s*s;
	double p_r = 1.0 + r;
	double p_s = 1.0 + s;
	double dNr5 = -r*m_s;
	double dNr6 = 0.5*m_ss;
	double dNr7 = -r*p_s;
	double dNr8 = -0.5*m_ss;
	double dNs5 = -0.5*m_rr;
	double dNs6 = -s*p_r;
	double dNs7 = 0.5*m_rr;
	double dNs8 = -s*m_r;


	dN[0][0] = -0.25*m_s - 0.5*dNr8 - 0.5*dNr5;
	dN[0][1] =  0.25*m_s - 0.5*dNr5 - 0.5*dNr6;
	dN[0][2] =  0.25*p_s - 0.5*dNr6 - 0.5*dNr7;
	dN[0][3] = -0.25*p_s - 0.5*dNr7 - 0.5*dNr8;
	dN[0][4] = dNr5;
	dN[0][5] = dNr6;
	dN[0][6] = dNr7;
	dN[0][7] = dNr8;

	dN[1][0] = -0.25*m_r - 0.5*dNs8 - 0.5*dNs5;
	dN[1][1] = -0.25*p_r - 0.5*dNs5 - 0.5*dNs6;
	dN[1][2] =  0.25*p_r - 0.5*dNs6 - 0.5*dNs7;
	dN[1][3] =  0.25*m_r - 0.5*dNs7 - 0.5*dNs8;
	dN[1][4] = dNs5;
	dN[1][5] = dNs6;
	dN[1][6] = dNs7;
	dN[1][7] = dNs8;

	return(NORMAL_RC);
}


static RC
dN_mitc9 (const double *rst, double **dN)
{
	double r = rst[0];
	double s = rst[1];
	double m_r = 1.0 - r;
	double m_s = 1.0 - s;
	double m_rr = 1.0 - r*r;
	double m_ss = 1.0 - s*s;
	double p_r = 1.0 + r;
	double p_s = 1.0 + s;
	double dNr9 = -2.0*r*m_ss;
	double dNr5 = -r*m_s - 0.5*dNr9;
	double dNr6 = 0.5*m_ss - 0.5*dNr9;
	double dNr7 = -r*p_s - 0.5*dNr9;
	double dNr8 = -0.5*m_ss - 0.5*dNr9;
	double dNs9 = -2.0*s*m_rr;
	double dNs5 = -0.5*m_rr - 0.5*dNs9;
	double dNs6 = -s*p_r - 0.5*dNs9;
	double dNs7 = 0.5*m_rr - 0.5*dNs9;
	double dNs8 = -s*m_r - 0.5*dNs9;

	dN[0][0] = -0.25*m_s - 0.5*dNr8 - 0.5*dNr5 - 0.25*dNr9;
	dN[0][1] =  0.25*m_s - 0.5*dNr5 - 0.5*dNr6 - 0.25*dNr9;
	dN[0][2] =  0.25*p_s - 0.5*dNr6 - 0.5*dNr7 - 0.25*dNr9;
	dN[0][3] = -0.25*p_s - 0.5*dNr7 - 0.5*dNr8 - 0.25*dNr9;
	dN[0][4] = dNr5;
	dN[0][5] = dNr6;
	dN[0][6] = dNr7;
	dN[0][7] = dNr8;
	dN[0][8] = dNr9;

	dN[1][0] = -0.25*m_r - 0.5*dNs8 - 0.5*dNs5 - 0.25*dNs9;
	dN[1][1] = -0.25*p_r - 0.5*dNs5 - 0.5*dNs6 - 0.25*dNs9;
	dN[1][2] =  0.25*p_r - 0.5*dNs6 - 0.5*dNs7 - 0.25*dNs9;
	dN[1][3] =  0.25*m_r - 0.5*dNs7 - 0.5*dNs8 - 0.25*dNs9;
	dN[1][4] = dNs5;
	dN[1][5] = dNs6;
	dN[1][6] = dNs7;
	dN[1][7] = dNs8;
	dN[1][8] = dNs9;

	return(NORMAL_RC);
}


