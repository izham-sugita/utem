/*********************************************************************
 * nl_solver_shell.c
 *
 * Copyright (C) 2010 AzLib Developers Group
 *
 * Written by
 *  <Taiki AOYAMA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: nl_solver_shell.c 1000 2012-01-29 08:20:16Z aoyama $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "base.h"
#include "mathematics.h"
#include "fem_solver.h"
#include "fem_struct.h"
#include "sky_cholesky.h"
#include "nst_component.h"

#define NR_MAX (30)         /* Newton-Raphson法最大反復回数 */


static RC set_rest_zero6(NODE_ARRAY node, BC_ARRAY rest, double resid_vect[]);
static RC make_large_rotation_matrix_mitc3(ELEM_TYPE elem_type,
                                           int elem_node_num,
                                           DIRECTOR_VECT elem_curr_vect[], 
                                           double thickness,
                                           double **update_location,
                                           const double *local_point,
                                           double stress_vect[],
                                           double *Kex,
                                           MITC_TYING_SET mitc, WORK_SET ws);
static RC make_large_rotation_matrix_mitc6(ELEM_TYPE elem_type,
                                           int elem_node_num,
                                           DIRECTOR_VECT elem_curr_vect[], 
                                           double thickness,
                                           double **update_location,
                                           const double *local_point,
                                           double stress_vect[],
                                           double *Kex,
                                           MITC_TYING_SET mitc, WORK_SET ws);
static RC make_large_rotation_matrix_mitc4(ELEM_TYPE elem_type,
                                           int elem_node_num,
                                           DIRECTOR_VECT elem_curr_vect[], 
                                           double thickness,
                                           double **update_location,
                                           const double *local_point,
                                           double stress_vect[],
                                           double *Kex,
                                           MITC_TYING_SET mitc, WORK_SET ws);
static RC make_large_rotation_matrix_mitc8(ELEM_TYPE elem_type,
                                           int elem_node_num,
                                           DIRECTOR_VECT elem_curr_vect[],
                                           double thickness,
                                           double **update_location,
                                           const double *local_point,
                                           double stress_vect[],
                                           double *Kex,
                                           MITC_TYING_SET mitc, WORK_SET ws);
static RC make_SZtZ_matrix_mitc3(ELEM_TYPE elem_type, int elem_node_num,
                                 DIRECTOR_VECT elem_curr_vect[],
                                 double thickness,
                                 const double *local_point,
                                 double stress_vect[],
                                 double **SZtZ_matrix,
                                 MITC_TYING_SET mitc, WORK_SET ws);
static RC make_SZtZ_matrix_mitc6(ELEM_TYPE elem_type, int elem_node_num,
                                 DIRECTOR_VECT elem_curr_vect[],
                                 double thickness,
                                 const double *local_point,
                                 double stress_vect[],
                                 double **SZtZ_matrix,
                                 MITC_TYING_SET mitc, WORK_SET ws);
static RC make_SZtZ_matrix_mitc4(ELEM_TYPE elem_type, int elem_node_num,
                                 DIRECTOR_VECT elem_curr_vect[],
                                 double thickness,
                                 const double *local_point,
                                 double stress_vect[],
                                 double **SZtZ_matrix,
                                 MITC_TYING_SET mitc, WORK_SET ws);
static RC make_SZtZ_matrix_mitc8(ELEM_TYPE elem_type, int elem_node_num,
                                 DIRECTOR_VECT elem_curr_vect[],
                                 double thickness,
                                 const double *local_point,
                                 double stress_vect[],
                                 double **SZtZ_matrix,
                                 MITC_TYING_SET mitc, WORK_SET ws);
static RC set_Z1_matrix(ELEM_TYPE elem_type, int elem_node_num,
                        DIRECTOR_VECT elem_curr_vect[],
                        double thickness, const double *local_point,
                        double **Z1_matrix, WORK_SET ws);
static RC set_Z2_matrix(ELEM_TYPE elem_type, int elem_node_num,
                        DIRECTOR_VECT elem_curr_vect[],
                        double thickness, const double *local_point,
                        double **Z2_matrix, WORK_SET ws);
static RC set_Z3_matrix(ELEM_TYPE elem_type, int elem_node_num,
                        DIRECTOR_VECT elem_curr_vect[],
                        double thickness, const double *local_point,
                        double **Z3_matrix, WORK_SET ws);
static RC local_cova_Green_strain_vect_mitc3(ELEM_TYPE elem_type,
                                             int elem_node_num,
                                             DIRECTOR_VECT elem_init_vect[],
                                             DIRECTOR_VECT elem_curr_vect[],
                                             double thickness,
                                             double **node_location,
                                             double **update_location,
                                             const double *local_point,
                                             double strain_vect[],
                                             MITC_TYING_SET mitc, WORK_SET ws);
static RC local_cova_Green_strain_vect_mitc6(ELEM_TYPE elem_type,
                                             int elem_node_num,
                                             DIRECTOR_VECT elem_init_vect[],
                                             DIRECTOR_VECT elem_curr_vect[],
                                             double thickness,
                                             double **node_location,
                                             double **update_location,
                                             const double *local_point,
                                             double strain_vect[],
                                             MITC_TYING_SET mitc, WORK_SET ws);
static RC local_cova_Green_strain_vect_mitc4(ELEM_TYPE elem_type,
                                             int elem_node_num,
                                             DIRECTOR_VECT elem_init_vect[],
                                             DIRECTOR_VECT elem_curr_vect[],
                                             double thickness,
                                             double **node_location,
                                             double **update_location,
                                             const double *local_point,
                                             double strain_vect[],
                                             MITC_TYING_SET mitc, WORK_SET ws);
static RC local_cova_Green_strain_vect_mitc8(ELEM_TYPE elem_type,
                                             int elem_node_num,
                                             DIRECTOR_VECT elem_init_vect[],
                                             DIRECTOR_VECT elem_curr_vect[],
                                             double thickness,
                                             double **node_location,
                                             double **update_location,
                                             const double *local_point,
                                             double strain_vect[],
                                             MITC_TYING_SET mitc, WORK_SET ws);
static RC set_drilling_dof(int mitc_node_num, ELEMENT *element,
                           NODE_ARRAY node, MATERIAL_PROP_ARRAY material,
                           PHYSICAL_PROP_ARRAY physical, double thickness,
                           ELEM_MATRIX *elem_matrix);
static RC trans_global_rotation_dof(int matrix_size, double **matrix,
                                    DIRECTOR_VECT elem_vect[]);
static RC set_Q_matrix4stress(VECT3D G[], double **Q_matrix);
static RC set_Q_matrix4strain(VECT3D g[], double **Q_matrix);


/*
 * Reference book is
 * 非線形有限要素法の基礎と応用
 * 久田俊明，野口裕久
 *
 * シェル要素用微小ひずみ大変形解析 Total Lagrange法ソルバー
 * ただし，座屈後解析不可
 * 応力は第2Piola-Kirchhoff応力，ひずみはGreenひずみを計算
 * 変位や応力だけでなく初期配置および現配置でのディレクターベクトルも計算
 * ロッキング回避のためにMITCシェル要素使用
 * 荷重制御(増分)法を使用しており，lcm_numには増分回数を指定
 * LCMはLoad Control Method(荷重制御法)の略
 *
 * 連立方程式の求解にはICCG法を使用するが，板厚が薄くなると正定値性が
 * 弱くなるので fem_NL_solve_shell_LCM の使用を推奨する．
 */
RC
fem_NL_solve_shell_iccg6LCM (NODE_ARRAY node, ELEMENT_ARRAY element,
                             RIGID_ELEMENT_ARRAY rigid,
                             MATERIAL_PROP_ARRAY material,
                             PHYSICAL_PROP_ARRAY physical,
                             BC_ARRAY rest, BC_ARRAY force,
                             BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                             DEFAULT_BC_ARRAY def_body_force,
                             DIRECTOR_VECT_ARRAY *initial_vect,
                             DIRECTOR_VECT_ARRAY *current_vect,
                             DISP_ARRAY *disp,
                             STRESS_ARRAY *stress, STRAIN_ARRAY *strain,
                             int sol_ss_flag, int sol_ss_face, int lcm_num)
{
	int ii1, ii2, ii3;
	int dim;
	int total_dof;
	double beta1 = 1.0e-4;
	double beta2 = 1.0e-7;
	double *f_vect;
	double *df_vect;
	double *react_vect;
	double *resid_vect;
	double *du_vect;
	double *u0_vect;
	double *scale_fact;
	NONZERO_MATRIX3 matrix;
	LOCAL_COORD_ARRAY dummy_coord;
	DISP_ARRAY delta_disp;
	WORK_SET ws1, ws2;
	MITC_TYING_SET mitc;


	if((dim = analysis_dim_shell(element)) != 6){
		return(ARG_ERROR_RC);
	}
	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

	total_dof = node.size * dim;
	RC_TRY( allocate1D(total_dof, &f_vect) );
	RC_TRY( allocate1D(total_dof, &df_vect) );
	RC_TRY( allocate1D(total_dof, &react_vect) );
	RC_TRY( allocate1D(total_dof, &resid_vect) );
	RC_TRY( allocate1D(total_dof, &du_vect) );
	RC_TRY( allocate1D(total_dof, &u0_vect) );
	RC_TRY( allocate1D(total_dof, &scale_fact) );
	RC_TRY( allocate_local_coord_array(0, &dummy_coord) );
	RC_TRY( vect2disp(dim, node, du_vect, disp) );
	RC_TRY( vect2disp(dim, node, du_vect, &delta_disp) );

	/* 初期配置，現配置でのディレクターベクトルの初期化 */
	RC_TRY( make_nodal_director_vect(node, element, initial_vect) );
	RC_TRY( allocate_director_vect_array(node.size, current_vect) );
	RC_TRY( copy_director_vect_array(*initial_vect, current_vect) );

	/* 節点荷重・モーメント */
	RC_TRY( add_nodal_force(dim, node, force, f_vect) );
	/* 未対応 */
#if 0
	/* 熱荷重 */
	RC_TRY( add_temp_force(dim, node, element, material, physical,
	                       temp, def_temp, f_vect) );
#endif
	/* 重力荷重 */
	RC_TRY( add_body_force_shell(node, element, material, physical,
	                             def_body_force, *initial_vect, f_vect) );
	
	/* 増分荷重ベクトル */
	for(ii1=0; ii1<total_dof; ii1++){
		df_vect[ii1] = f_vect[ii1] / (double)lcm_num;
		f_vect[ii1] = 0.0;
	}

	RC_TRY( allocate_work_set(&ws1) );
	RC_TRY( allocate_work_set(&ws2) );
	RC_TRY( allocate_mitc_tying_set(&mitc) );
	for(ii1=0; ii1<lcm_num; ii1++){
		RC_TRY( log_printf(4, " _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/\n") );
		RC_TRY( log_printf(4, "   Load Incremental Number :%4d /%4d   \n",
		                   ii1+1, lcm_num) );
		RC_TRY( log_printf(4, " _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/\n") );

		for(ii2=0; ii2<total_dof; ii2++){
			f_vect[ii2] += df_vect[ii2];
			resid_vect[ii2] = f_vect[ii2] - react_vect[ii2];
		}

		for(ii2=0; ii2<NR_MAX; ii2++){
			double du_norm, du_r, u0_norm, u0_f;
			ELEM_MATRIX elem_matrix;

			RC_TRY( log_printf(5, " +-------------------------+\n") );
			RC_TRY( log_printf(5, " | N-R iteration NUM : %3d |\n", ii2+1) );
			RC_TRY( log_printf(5, " +-------------------------+\n") );

			/* 接線剛性マトリクス作成 */
			RC_TRY( allocate_nonzero6(node, element, rigid, &matrix) );
			for(ii3=0; ii3<element.size; ii3++){
				if(element.array[ii3].label < 0) continue;

				RC_TRY( make_NL_shell_tangent_elem_matrix(element.array[ii3],
				                                          node, material,
				                                          physical,
				                                          *initial_vect,
				                                          *current_vect, *disp,
				                                          &elem_matrix, mitc,
				                                          ws1, ws2) );
				RC_TRY( impose_nonzero6(matrix, elem_matrix, 1.0) );
				RC_TRY( free_elem_matrix(&elem_matrix) );
			}

			RC_TRY( modify_nonzero6(node, dummy_coord, rest, rigid, matrix,
			                        &resid_vect, 1, NULL) );
			for(ii3=0; ii3<total_dof; ii3++){
				du_vect[ii3] = 0.0;
				scale_fact[ii3] = 0.0;
			}
			RC_TRY( nonzero3s_scaling_pre(matrix, scale_fact, resid_vect) );
			RC_TRY( iccg_nonzero3s(matrix, resid_vect, du_vect, 1, 1) );
			RC_TRY( nonzero3s_scaling_post(matrix.size, scale_fact, du_vect) );
			RC_TRY( free_nonzero_matrix3(&matrix) );
			RC_TRY( rigid_recover6(node, rigid, du_vect) );

			RC_TRY( add_vect2disp(dim, node, du_vect, *disp, 1.0) );
			RC_TRY( copy_vect2disp(dim, node, du_vect, delta_disp, 1.0) );

			if(ii2 == 0){
				for(ii3=0; ii3<total_dof; ii3++){
					u0_vect[ii3] = du_vect[ii3];
				}
			}

			/* 現配置でのディレクターベクトル更新 */
			RC_TRY( update_current_director_vect(delta_disp, current_vect) );

			/* 内力(反力)計算 */
			RC_TRY( make_NL_shell_react_force(element, node, material,
			                                  physical, *initial_vect,
			                                  *current_vect, *disp,
			                                  react_vect) );

			/* 残差荷重ベクトル */
			for(ii3=0; ii3<total_dof; ii3++){
				resid_vect[ii3] = f_vect[ii3] - react_vect[ii3];
			}
			RC_TRY( set_rest_zero6(node, rest, resid_vect) );

			/* 変位，エネルギー収束判定 */
			u0_norm = sqrt( inner_product(total_dof, u0_vect, u0_vect) );
			du_norm = sqrt( inner_product(total_dof, du_vect, du_vect) );
			u0_f = 0.0;
			du_r = 0.0;
			for(ii3=0; ii3<total_dof; ii3++){
				u0_f += fabs(u0_vect[ii3] * f_vect[ii3]);
				du_r += fabs(du_vect[ii3] * resid_vect[ii3]);
			}

			RC_TRY( log_printf(4, "    Disp   error :%15.7e\n",
			                   du_norm/u0_norm) );
			RC_TRY( log_printf(4, "    Energy error :%15.7e\n", du_r/u0_f) );

			/* 収束判定 */
			if(du_norm < beta1*u0_norm + ABS_TOL
			&& du_r < beta2*u0_f + ABS_TOL){
				break;
			}
			RC_TRY( log_printf(4, "\n") );
		}
		if(ii2 == NR_MAX){
			log_printf(3, "N-R iteration NOT converged!!\n");
			return(CAL_ERROR_RC);
		}
		RC_TRY( log_printf(4, "\n") );
		RC_TRY( log_printf(4, "\n") );
	}
	RC_TRY( free_work_set(&ws1) );
	RC_TRY( free_work_set(&ws2) );
	RC_TRY( free_mitc_tying_set(&mitc) );

	/* メモリ開放 */
	RC_TRY( free1D(total_dof, &f_vect) );
	RC_TRY( free1D(total_dof, &df_vect) );
	RC_TRY( free1D(total_dof, &react_vect) );
	RC_TRY( free1D(total_dof, &resid_vect) );
	RC_TRY( free1D(total_dof, &du_vect) );
	RC_TRY( free1D(total_dof, &u0_vect) );
	RC_TRY( free1D(total_dof, &scale_fact) );
	RC_TRY( free_local_coord_array(&dummy_coord) );
	RC_TRY( free_disp_array(&delta_disp) );

	/* 第2Piola-Kirchhoff応力，Greenひずみ */
	RC_TRY( PK2_stress_Green_strain_shell(sol_ss_flag, sol_ss_face,
	                                      node, element, material, physical,
	                                      *initial_vect, *current_vect,
	                                      *disp, stress, strain) );


	return(NORMAL_RC);
}


/*
 * Reference book is
 * 非線形有限要素法の基礎と応用
 * 久田俊明，野口裕久
 *
 * シェル要素用微小ひずみ大変形解析 Total Lagrange法ソルバー
 * ただし，座屈後解析不可
 * 応力は第2Piola-Kirchhoff応力，ひずみはGreenひずみを計算
 * 変位や応力だけでなく初期配置および現配置でのディレクターベクトルも計算
 * ロッキング回避のためにMITCシェル要素使用
 * 荷重制御(増分)法を使用しており，lcm_numには増分回数を指定
 * LCMはLoad Control Method(荷重制御法)の略
 *
 * 連立方程式の求解にはスカイライン法使用
 */
RC
fem_NL_solve_shell_LCM (NODE_ARRAY node, ELEMENT_ARRAY element,
                        RIGID_ELEMENT_ARRAY rigid,
                        MATERIAL_PROP_ARRAY material,
                        PHYSICAL_PROP_ARRAY physical,
                        BC_ARRAY rest, BC_ARRAY force,
                        BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                        DEFAULT_BC_ARRAY def_body_force,
                        DIRECTOR_VECT_ARRAY *initial_vect,
                        DIRECTOR_VECT_ARRAY *current_vect,
                        DISP_ARRAY *disp,
                        STRESS_ARRAY *stress, STRAIN_ARRAY *strain,
                        int sol_ss_flag, int sol_ss_face, int lcm_num)
{
	int ii1, ii2, ii3;
	int dim;
	int total_dof;
	double beta1 = 1.0e-4;
	double beta2 = 1.0e-7;
	double *f_vect;
	double *df_vect;
	double *react_vect;
	double *resid_vect;
	double *du_vect;
	double *u0_vect;
	SKYLINE_MATRIX matrix;
	DISP_ARRAY delta_disp;
	WORK_SET ws1, ws2;
	MITC_TYING_SET mitc;


	if((dim = analysis_dim_shell(element)) != 6){
		return(ARG_ERROR_RC);
	}
	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

	total_dof = node.size * dim;
	RC_TRY( allocate1D(total_dof, &f_vect) );
	RC_TRY( allocate1D(total_dof, &df_vect) );
	RC_TRY( allocate1D(total_dof, &react_vect) );
	RC_TRY( allocate1D(total_dof, &resid_vect) );
	RC_TRY( allocate1D(total_dof, &du_vect) );
	RC_TRY( allocate1D(total_dof, &u0_vect) );
	RC_TRY( vect2disp(dim, node, du_vect, disp) );
	RC_TRY( vect2disp(dim, node, du_vect, &delta_disp) );

	/* 初期配置，現配置でのディレクターベクトルの初期化 */
	RC_TRY( make_nodal_director_vect(node, element, initial_vect) );
	RC_TRY( allocate_director_vect_array(node.size, current_vect) );
	RC_TRY( copy_director_vect_array(*initial_vect, current_vect) );

	/* 節点荷重・モーメント */
	RC_TRY( add_nodal_force(dim, node, force, f_vect) );
	/* 未対応 */
#if 0
	/* 熱荷重 */
	RC_TRY( add_temp_force(dim, node, element, material, physical,
	                       temp, def_temp, f_vect) );
#endif
	/* 重力荷重 */
	RC_TRY( add_body_force_shell(node, element, material, physical,
	                             def_body_force, *initial_vect, f_vect) );
	
	/* 増分荷重ベクトル */
	for(ii1=0; ii1<total_dof; ii1++){
		df_vect[ii1] = f_vect[ii1] / (double)lcm_num;
		f_vect[ii1] = 0.0;
	}

	RC_TRY( allocate_work_set(&ws1) );
	RC_TRY( allocate_work_set(&ws2) );
	RC_TRY( allocate_mitc_tying_set(&mitc) );
	for(ii1=0; ii1<lcm_num; ii1++){
		RC_TRY( log_printf(4, " _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/\n") );
		RC_TRY( log_printf(4, "   Load Incremental Number :%4d /%4d   \n",
		                   ii1+1, lcm_num) );
		RC_TRY( log_printf(4, " _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/\n") );

		for(ii2=0; ii2<total_dof; ii2++){
			f_vect[ii2] += df_vect[ii2];
			resid_vect[ii2] = f_vect[ii2] - react_vect[ii2];
		}

		for(ii2=0; ii2<NR_MAX; ii2++){
			double du_norm, du_r, u0_norm, u0_f;
			ELEM_MATRIX elem_matrix;

			RC_TRY( log_printf(5, " +-------------------------+\n") );
			RC_TRY( log_printf(5, " | N-R iteration NUM : %3d |\n", ii2+1) );
			RC_TRY( log_printf(5, " +-------------------------+\n") );

			/* 接線剛性マトリクス作成 */
			RC_TRY( allocate_g_matrix(dim, node, element, &matrix) );
			for(ii3=0; ii3<element.size; ii3++){
				if(element.array[ii3].label < 0) continue;

				RC_TRY( make_NL_shell_tangent_elem_matrix(element.array[ii3],
				                                          node, material,
				                                          physical,
				                                          *initial_vect,
				                                          *current_vect, *disp,
				                                          &elem_matrix, mitc,
				                                          ws1, ws2) );
				RC_TRY( impose_elem_matrix(matrix, elem_matrix, 1.0) );
				RC_TRY( free_elem_matrix(&elem_matrix) );
			}

			RC_TRY( modify_g_matrix6(node, rest, matrix, resid_vect) );
			RC_TRY( decomp_g_matrix(matrix) );
			for(ii3=0; ii3<total_dof; ii3++){
				du_vect[ii3] = resid_vect[ii3];
			}
			RC_TRY( s_cholesky_subst(matrix.dof, matrix.index1, matrix.index2,
			                         matrix.array, du_vect, 1) );
			RC_TRY( free_g_matrix(&matrix) );

			RC_TRY( add_vect2disp(dim, node, du_vect, *disp, 1.0) );
			RC_TRY( copy_vect2disp(dim, node, du_vect, delta_disp, 1.0) );

			if(ii2 == 0){
				for(ii3=0; ii3<total_dof; ii3++){
					u0_vect[ii3] = du_vect[ii3];
				}
			}

			/* 現配置でのディレクターベクトル更新 */
			RC_TRY( update_current_director_vect(delta_disp, current_vect) );

			/* 内力(反力)計算 */
			RC_TRY( make_NL_shell_react_force(element, node, material,
			                                  physical, *initial_vect,
			                                  *current_vect, *disp,
			                                  react_vect) );

			/* 残差荷重ベクトル */
			for(ii3=0; ii3<total_dof; ii3++){
				resid_vect[ii3] = f_vect[ii3] - react_vect[ii3];
			}
			RC_TRY( set_rest_zero6(node, rest, resid_vect) );

			/* 変位，エネルギー収束判定 */
			u0_norm = sqrt( inner_product(total_dof, u0_vect, u0_vect) );
			du_norm = sqrt( inner_product(total_dof, du_vect, du_vect) );
			u0_f = 0.0;
			du_r = 0.0;
			for(ii3=0; ii3<total_dof; ii3++){
				u0_f += fabs(u0_vect[ii3] * f_vect[ii3]);
				du_r += fabs(du_vect[ii3] * resid_vect[ii3]);
			}

			RC_TRY( log_printf(4, "    Disp   error :%15.7e\n",
			                   du_norm/u0_norm) );
			RC_TRY( log_printf(4, "    Energy error :%15.7e\n", du_r/u0_f) );

			/* 収束判定 */
			if(du_norm < beta1*u0_norm + ABS_TOL
			&& du_r < beta2*u0_f + ABS_TOL){
				break;
			}
			RC_TRY( log_printf(4, "\n") );
		}
		if(ii2 == NR_MAX){
			log_printf(3, "N-R iteration NOT converged!!\n");
			return(CAL_ERROR_RC);
		}
		RC_TRY( log_printf(4, "\n") );
		RC_TRY( log_printf(4, "\n") );
	}
	RC_TRY( free_work_set(&ws1) );
	RC_TRY( free_work_set(&ws2) );
	RC_TRY( free_mitc_tying_set(&mitc) );

	/* メモリ開放 */
	RC_TRY( free1D(total_dof, &f_vect) );
	RC_TRY( free1D(total_dof, &df_vect) );
	RC_TRY( free1D(total_dof, &react_vect) );
	RC_TRY( free1D(total_dof, &resid_vect) );
	RC_TRY( free1D(total_dof, &du_vect) );
	RC_TRY( free1D(total_dof, &u0_vect) );
	RC_TRY( free_disp_array(&delta_disp) );

	/* 第2Piola-Kirchhoff応力，Greenひずみ */
	RC_TRY( PK2_stress_Green_strain_shell(sol_ss_flag, sol_ss_face,
	                                      node, element, material, physical,
	                                      *initial_vect, *current_vect,
	                                      *disp, stress, strain) );


	return(NORMAL_RC);
}


static RC
set_rest_zero6 (NODE_ARRAY node, BC_ARRAY rest, double resid_vect[])
{
	int ii1;


	for(ii1=0; ii1<rest.size; ii1++){
		int index;

		if(rest.array[ii1].node < 0) continue;

		index = node_renum_index(node, rest.array[ii1].node);
		if(index < 0) return(SEEK_ERROR_RC);

		if(rest.array[ii1].v_type.x == BC_FIX) resid_vect[6*index] = 0.0;
		if(rest.array[ii1].v_type.y == BC_FIX) resid_vect[6*index+1] = 0.0;
		if(rest.array[ii1].v_type.z == BC_FIX) resid_vect[6*index+2] = 0.0;
		if(rest.array[ii1].v_type.yz == BC_FIX) resid_vect[6*index+3] = 0.0;
		if(rest.array[ii1].v_type.zx == BC_FIX) resid_vect[6*index+4] = 0.0;
		if(rest.array[ii1].v_type.xy == BC_FIX) resid_vect[6*index+5] = 0.0;
	}

	return(NORMAL_RC);
}


/* Total Lagrange法の内力ベクトル(等価節点力ベクトル) */
RC
make_NL_shell_react_force (ELEMENT_ARRAY element, NODE_ARRAY node,
                           MATERIAL_PROP_ARRAY material,
                           PHYSICAL_PROP_ARRAY physical,
                           DIRECTOR_VECT_ARRAY initial_vect,
                           DIRECTOR_VECT_ARRAY current_vect,
                           DISP_ARRAY disp, double react_vect[])
{
	int ii1, ii2;
	int dim;
	int total_dof;
	int elem_force_index[MAX_SHELL_NODE*6];
	double elem_force[MAX_SHELL_NODE*6];
	MITC_TYING_SET mitc;
	WORK_SET ws1, ws2;


	if((dim = analysis_dim_shell(element)) != 6){
		return(ARG_ERROR_RC);
	}
	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

	RC_TRY( allocate_mitc_tying_set(&mitc) );
	RC_TRY( allocate_work_set(&ws1) );
	RC_TRY( allocate_work_set(&ws2) );

	total_dof = dim * node.size;

	for(ii1=0; ii1<total_dof; ii1++) react_vect[ii1] = 0.0;

	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		RC_TRY( make_NL_shell_elem_react_force(element.array[ii1], node,
		                                       material, physical,
		                                       initial_vect, current_vect,
		                                       disp, elem_force_index,
		                                       elem_force, mitc, ws1, ws2) );

		for(ii2=0; ii2<element.array[ii1].node_num*dim; ii2++){
			int index = elem_force_index[ii2];

			react_vect[index] += elem_force[ii2];
		}
	}

	RC_TRY( free_mitc_tying_set(&mitc) );
	RC_TRY( free_work_set(&ws1) );
	RC_TRY( free_work_set(&ws2) );

	return(NORMAL_RC);
}



/* Total Lagrange法の要素内力ベクトル(等価節点力ベクトル) */
RC
make_NL_shell_elem_react_force (ELEMENT element, NODE_ARRAY node,
                                MATERIAL_PROP_ARRAY material,
                                PHYSICAL_PROP_ARRAY physical,
                                DIRECTOR_VECT_ARRAY initial_vect,
                                DIRECTOR_VECT_ARRAY current_vect,
                                DISP_ARRAY disp,
                                int elem_force_index[], double elem_force[],
                                MITC_TYING_SET mitc,
                                WORK_SET ws1, WORK_SET ws2)
{
	int ii1, ii2;
	int dim;
	int ssnum;
	int elem_node_num;
	int g_num_point;
	int prop_index;
	double thickness;
	double det_J;
	double strain_vect[5];
	double stress_vect[5];
	double react_v[MAX_SHELL_NODE*5];
	double local_force[MAX_SHELL_NODE*5];
	double *g_weights = ws1.vec_G[0];
	double **node_location = ws1.mat_N_3[0];
	double **update_location = ws1.mat_N_3[1];
	double **g_points = ws1.mat_G_4[0];
	double **D_matrix = ws1.mat_6_6[0];
	double **B_matrix = ws1.mat_6_3N[0];
	VECT3D G0[3], g0[3];
	DIRECTOR_VECT elem_init_vect[MAX_SHELL_NODE];
	DIRECTOR_VECT elem_curr_vect[MAX_SHELL_NODE];
	ELEM_TYPE elem_type;


	prop_index = search_physical_prop_label(physical, element.physical);
	if(prop_index < 0) return(SEEK_ERROR_RC);
	if(NST_PROP_PSHELL != physical.array[prop_index].i_info[1]){
		return(ARG_ERROR_RC);
	}

	if( (dim = element_dim_shell(element.type)) != 6
	|| (ssnum = element_ssnum_shell(element.type)) != 5 ){
		return(ARG_ERROR_RC);
	}
	
	for(ii1=0; ii1<element.node_num*dim; ii1++) elem_force[ii1] = 0.0;
	for(ii1=0; ii1<element.node_num; ii1++){
		int index;

		index = node_renum_index(node, element.node[ii1]);
		if(index < 0 ) return(SEEK_ERROR_RC);

		for(ii2=0; ii2<dim; ii2++){
			elem_force_index[ii1*dim+ii2] = index*dim + ii2;
		}
	}

	elem_node_num = element.node_num;
	elem_type = element.type;

	RC_TRY( set_element_direcor_vect(element, initial_vect, elem_init_vect) );
	RC_TRY( set_element_direcor_vect(element, current_vect, elem_curr_vect) );
	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( update_node_location_matrix(&element,node,disp,update_location) );

	RC_TRY( Gauss_const_default_shell(elem_type, g_points, g_weights,
	                                  &g_num_point) );
	RC_TRY( get_thickness(element, physical, &thickness) );

	init_vect(elem_node_num*5, local_force);
	for(ii1=0; ii1<g_num_point; ii1++){
		RC_TRY( make_G1G2G3g1g2g3_basis(elem_type, elem_node_num,
		                                elem_init_vect, thickness,
		                                node_location, g_points[ii1],
		                                G0, g0, ws2) );

		RC_TRY( make_B_matrix_shell(elem_type, elem_node_num,
		                            elem_curr_vect, thickness,
		                            update_location,
		                            g_points[ii1], B_matrix, mitc, ws2) );

		RC_TRY( make_D_matrix_shell(element, material, physical,
		                            g0, D_matrix) );

		RC_TRY( local_covariant_Green_strain_vect(elem_type, elem_node_num,
		                                          elem_init_vect,
		                                          elem_curr_vect, thickness,
		                                          node_location,
		                                          update_location,
		                                          g_points[ii1],
		                                          strain_vect, mitc, ws2) );

		mul_matrix_vect(ssnum, ssnum, D_matrix, strain_vect, stress_vect);

		mul_trans_matrix_vect(ssnum, 5*elem_node_num, B_matrix, stress_vect,
		                      react_v);

		det_J = scalar_triple_product3d(G0[0], G0[1], G0[2]);

		for(ii2=0; ii2<elem_node_num*5; ii2++){
			local_force[ii2] += g_weights[ii1] * det_J * react_v[ii2];
		}
	}

	/* 全体座標系成分へ変換 */
	for(ii1=0; ii1<elem_node_num; ii1++){
		VECT3D V1 = elem_curr_vect[ii1].v1;
		VECT3D V2 = elem_curr_vect[ii1].v2;

		elem_force[ii1*dim+0] = local_force[ii1*5+0];
		elem_force[ii1*dim+1] = local_force[ii1*5+1];
		elem_force[ii1*dim+2] = local_force[ii1*5+2];
		elem_force[ii1*dim+3] = V1.x*local_force[ii1*5+3]
		                      + V2.x*local_force[ii1*5+4];
		elem_force[ii1*dim+4] = V1.y*local_force[ii1*5+3]
		                      + V2.y*local_force[ii1*5+4];
		elem_force[ii1*dim+5] = V1.z*local_force[ii1*5+3]
		                      + V2.z*local_force[ii1*5+4];
	}

	return(NORMAL_RC);
}


/*
 * Reference book is
 * 非線形有限要素法の基礎と応用
 * 久田俊明，野口裕久
 *
 * 要素接線剛性マトリクス作成(Total Lagrange法)
 */
RC
make_NL_shell_tangent_elem_matrix (ELEMENT element, NODE_ARRAY node,
                                   MATERIAL_PROP_ARRAY material,
                                   PHYSICAL_PROP_ARRAY physical,
                                   DIRECTOR_VECT_ARRAY initial_vect,
                                   DIRECTOR_VECT_ARRAY current_vect,
                                   DISP_ARRAY disp,
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
	double strain_vect[5];
	double stress_vect[5];
	double *Kex = ws1.vec_N[0];
	double *g_weights = ws1.vec_G[0];
	double **node_location = ws1.mat_N_3[0];
	double **update_location = ws1.mat_N_3[1];
	double **g_points = ws1.mat_G_4[0];
	double **D_matrix = ws1.mat_6_6[0];
	double **B_matrix = ws1.mat_6_3N[0];
	double **BtDB = ws1.mat_3N_3N[0];
	double **SZtZ = ws1.mat_3N_3N[1];
	VECT3D G0[3], g0[3];
	DIRECTOR_VECT elem_init_vect[MAX_SHELL_NODE];
	DIRECTOR_VECT elem_curr_vect[MAX_SHELL_NODE];
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

	elem_node_num = element.node_num;
	elem_type = element.type;

	elem_matrix_size = elem_node_num * dim;
	RC_TRY( allocate2D(elem_matrix_size, elem_matrix_size,
	                   &(elem_matrix->matrix)) );

	RC_TRY( set_element_direcor_vect(element, initial_vect, elem_init_vect) );
	RC_TRY( set_element_direcor_vect(element, current_vect, elem_curr_vect) );
	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( update_node_location_matrix(&element,node,disp,update_location) );
	RC_TRY( Gauss_const_default_shell(elem_type, g_points, g_weights,
	                                  &g_num_point) );
	RC_TRY( get_thickness(element, physical, &thickness) );

	for(ii1=0; ii1<g_num_point; ii1++){
		RC_TRY( make_G1G2G3g1g2g3_basis(elem_type, elem_node_num,
		                                elem_init_vect, thickness,
		                                node_location, g_points[ii1],
		                                G0, g0, ws2) );

		RC_TRY( make_B_matrix_shell(elem_type, elem_node_num,
		                            elem_curr_vect, thickness,
		                            update_location,
		                            g_points[ii1], B_matrix, mitc, ws2) );

		RC_TRY( make_D_matrix_shell(element, material, physical,
		                            g0, D_matrix) );

		mul_matrix_AtBA(ssnum, 5*elem_node_num, B_matrix, D_matrix, BtDB);

		RC_TRY( local_covariant_Green_strain_vect(elem_type, elem_node_num,
		                                          elem_init_vect,
		                                          elem_curr_vect, thickness,
		                                          node_location,
		                                          update_location,
		                                          g_points[ii1],
		                                          strain_vect, mitc, ws2) );

		mul_matrix_vect(ssnum, ssnum, D_matrix, strain_vect, stress_vect);

		/* 初期応力 SZtZ */
		RC_TRY( make_SZtZ_matrix(elem_type, elem_node_num,
		                         elem_curr_vect, thickness, g_points[ii1],
		                         stress_vect, SZtZ, mitc, ws2) );

		/* 有限回転 */
		RC_TRY( make_large_rotation_matrix(elem_type, elem_node_num,
                                           elem_curr_vect, thickness,
		                                   update_location, g_points[ii1],
                                           stress_vect, Kex, mitc, ws2) );

		det_J = scalar_triple_product3d(G0[0], G0[1], G0[2]);

		for(ii2=0; ii2<elem_node_num; ii2++){
			for(ii3=0; ii3<elem_node_num; ii3++){
				for(ii4=0; ii4<5; ii4++){
					for(ii5=0; ii5<5; ii5++){
						elem_matrix->matrix[dim*ii2+ii4][dim*ii3+ii5]
						+= (BtDB[5*ii2+ii4][5*ii3+ii5]
						  + SZtZ[5*ii2+ii4][5*ii3+ii5])*det_J*g_weights[ii1];
					}
				}
			}
			elem_matrix->matrix[dim*ii2+3][dim*ii2+3]
			+= Kex[ii2] * det_J * g_weights[ii1];
			elem_matrix->matrix[dim*ii2+4][dim*ii2+4]
			+= Kex[ii2] * det_J * g_weights[ii1];
		}
	}

	/* 虚偽のドリリング剛性 */
	RC_TRY( set_drilling_dof(elem_node_num, &element, node, material, physical,
	                         thickness, elem_matrix) );

	/* 全体座標系成分へ変換 */
	RC_TRY( trans_global_rotation_dof(elem_matrix_size, elem_matrix->matrix,
	                                  elem_curr_vect) );

	return(NORMAL_RC);
}


/* 増分変位による現配置でのディレクターベクトルの更新 */
RC
update_current_director_vect (DISP_ARRAY delta_disp,
                              DIRECTOR_VECT_ARRAY *current_vect)
{
	int ii1, ii2, ii3, ii4;
	double omega, omega2;
	double I[3][3];
	double R[3][3];
	double phi[3][3];
	double phi2[3][3];


	for(ii1=0; ii1<3; ii1++){
		for(ii2=0; ii2<3; ii2++){
			I[ii1][ii2] = 0.0;
			R[ii1][ii2] = 0.0;
			phi[ii1][ii2] = 0.0;
		}
		I[ii1][ii1] = 1.0;
	}

	for(ii1=0; ii1<current_vect->size; ii1++){
		int index;
		double factor1, factor2;
		double alpha, beta, gamma;
		double theta1, theta2, theta3;
		VECT3D V1, V2, V3;
		VECT3D update_V1, update_V2, update_V3;


		index = search_disp_node(delta_disp, current_vect->array[ii1].label);
		if(index < 0) return(SEEK_ERROR_RC);

		V1 = current_vect->array[ii1].v1;
		V2 = current_vect->array[ii1].v2;
		V3 = current_vect->array[ii1].vn;

		alpha = V1.x*delta_disp.array[index].v.yz
		      + V1.y*delta_disp.array[index].v.zx
		      + V1.z*delta_disp.array[index].v.xy;
		beta = V2.x*delta_disp.array[index].v.yz
		     + V2.y*delta_disp.array[index].v.zx
		     + V2.z*delta_disp.array[index].v.xy;
		gamma = 0.0;

		theta1 = V1.x*alpha + V2.x*beta;
		theta2 = V1.y*alpha + V2.y*beta;
		theta3 = V1.z*alpha + V2.z*beta;
#if 0
		gamma = V3.x*delta_disp.array[index].v.yz
		      + V3.y*delta_disp.array[index].v.zx
		      + V3.z*delta_disp.array[index].v.xy;

		theta1 = V1.x*alpha + V2.x*beta + V3.x*gamma;
		theta2 = V1.y*alpha + V2.y*beta + V3.y*gamma;
		theta3 = V1.z*alpha + V2.z*beta + V3.z*gamma;
#endif

		phi[0][1] = -theta3;
		phi[1][0] = theta3;
		phi[0][2] = theta2;
		phi[2][0] = -theta2;
		phi[1][2] = -theta1;
		phi[2][1] = theta1;

		for(ii2=0; ii2<3; ii2++){
			for(ii3=0; ii3<3; ii3++){
				phi2[ii2][ii3] = 0.0;
			}
		}
		for(ii2=0; ii2<3; ii2++){
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					phi2[ii2][ii4] += phi[ii2][ii3] * phi[ii3][ii4];
				}
			}
		}

		omega = sqrt(theta1*theta1 + theta2*theta2 + theta3*theta3);
		omega2 = 0.5 * omega;

		if( nearly_eq(0.0, omega) ){
			for(ii2=0; ii2<3; ii2++){
				for(ii3=0; ii3<3; ii3++){
					R[ii2][ii3] = I[ii2][ii3];
				}
			}
		}else{
			factor1 = sin(omega) / omega;
			factor2 = 0.5 * ((sin(omega2)/omega2)*(sin(omega2)/omega2));

			for(ii2=0; ii2<3; ii2++){
				for(ii3=0; ii3<3; ii3++){
					R[ii2][ii3] = I[ii2][ii3] + factor1*phi[ii2][ii3]
					                          + factor2*phi2[ii2][ii3];
				}
			}
		}

		/* Update V1 */
		update_V1.x = R[0][0]*V1.x + R[0][1]*V1.y + R[0][2]*V1.z;
		update_V1.y = R[1][0]*V1.x + R[1][1]*V1.y + R[1][2]*V1.z;
		update_V1.z = R[2][0]*V1.x + R[2][1]*V1.y + R[2][2]*V1.z;
		/* Update V2 */
		update_V2.x = R[0][0]*V2.x + R[0][1]*V2.y + R[0][2]*V2.z;
		update_V2.y = R[1][0]*V2.x + R[1][1]*V2.y + R[1][2]*V2.z;
		update_V2.z = R[2][0]*V2.x + R[2][1]*V2.y + R[2][2]*V2.z;
		/* Update V3(director) */
		update_V3.x = R[0][0]*V3.x + R[0][1]*V3.y + R[0][2]*V3.z;
		update_V3.y = R[1][0]*V3.x + R[1][1]*V3.y + R[1][2]*V3.z;
		update_V3.z = R[2][0]*V3.x + R[2][1]*V3.y + R[2][2]*V3.z;

		current_vect->array[ii1].v1 = unit_vect3d(update_V1);
		current_vect->array[ii1].v2 = unit_vect3d(update_V2);
		current_vect->array[ii1].vn = unit_vect3d(update_V3);
	}

	return(NORMAL_RC);
}


/*
 * Reference book is
 * 非線形有限要素法の基礎と応用
 * 久田俊明，野口裕久
 *
 * 要素接線剛性マトリクスの有限回転項
 */
RC
make_large_rotation_matrix (ELEM_TYPE elem_type, int elem_node_num,
                            DIRECTOR_VECT elem_curr_vect[], 
                            double thickness,
                            double **update_location,
                            const double *local_point,
                            double stress_vect[],
                            double *Kex, MITC_TYING_SET mitc, WORK_SET ws)
{
	switch(elem_type){
	case ELEM_TRI1:
		RC_TRY( make_large_rotation_matrix_mitc3(elem_type, elem_node_num,
		                                         elem_curr_vect, thickness,
		                                         update_location, local_point,
		                                         stress_vect, Kex, mitc, ws) );
		break;
	case ELEM_TRI2:
		RC_TRY( make_large_rotation_matrix_mitc6(elem_type, elem_node_num,
		                                         elem_curr_vect, thickness,
		                                         update_location, local_point,
		                                         stress_vect, Kex, mitc, ws) );
		break;
	case ELEM_QUAD1:
		RC_TRY( make_large_rotation_matrix_mitc4(elem_type, elem_node_num,
		                                         elem_curr_vect, thickness,
		                                         update_location, local_point,
		                                         stress_vect, Kex, mitc, ws) );
		break;
	case ELEM_QUAD2:
		RC_TRY( make_large_rotation_matrix_mitc8(elem_type, elem_node_num,
		                                         elem_curr_vect, thickness,
		                                         update_location, local_point,
		                                         stress_vect, Kex, mitc, ws) );
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


/* 要素接線剛性マトリクスの有限回転項 MITC3 */
static RC
make_large_rotation_matrix_mitc3 (ELEM_TYPE elem_type, int elem_node_num,
                                  DIRECTOR_VECT elem_curr_vect[], 
                                  double thickness,
                                  double **update_location,
                                  const double *local_point,
                                  double stress_vect[],
                                  double *Kex,
                                  MITC_TYING_SET mitc, WORK_SET ws)
{
	int ii1, ii2;
	int ErtEst_tying_num;
	double t = local_point[2];
	double c1, c2;
	double V3G1, V3G2;
	double V3G3rt_rt, V3G1rt_rt, V3G3rt_st, V3G2rt_st;
	double V3G3st_rt, V3G1st_rt, V3G3st_st, V3G2st_st;
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
	VECT3D V3;


	init_vect(elem_node_num, Kex);

	c1 = 0.5 * thickness * t;
	c2 = 0.5 * thickness;

	RC_TRY( set_dN_matrix_shell(elem_type, local_point, dN) );
	RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
	                          thickness, update_location, local_point, G, ws) );

	RC_TRY( set_tying_data_mitc3(local_point, &ErtEst_tying_num,
	                             Ert_rt_tying_point, Ert_st_tying_point,
	                             Est_rt_tying_point, Est_st_tying_point,
	                             Ert_rt_h, Ert_st_h, Est_rt_h, Est_st_h) );

	for(ii1=0; ii1<elem_node_num; ii1++){
		double Ert_rt, Ert_st;
		double Est_rt, Est_st;
		double stress_rr = stress_vect[0];
		double stress_ss = stress_vect[1];
		double stress_rs = stress_vect[2];
		double stress_rt = stress_vect[3];
		double stress_st = stress_vect[4];

		V3 = elem_curr_vect[ii1].vn;

		V3G1 = inner_product3d(V3, G[0]);
		V3G2 = inner_product3d(V3, G[1]);

		Kex[ii1] = - stress_rr*c1*dN[0][ii1]*V3G1
		           - stress_ss*c1*dN[1][ii1]*V3G2
		           - stress_rs*c1*(dN[0][ii1]*V3G2 + dN[1][ii1]*V3G1);

		for(ii2=0; ii2<ErtEst_tying_num; ii2++){
			RC_TRY( set_N_vector_shell(elem_type, Ert_rt_tying_point[ii2],
			                           Nrt_rt) );
			RC_TRY( set_N_vector_shell(elem_type, Ert_st_tying_point[ii2],
			                           Nrt_st) );
			RC_TRY( set_N_vector_shell(elem_type, Est_rt_tying_point[ii2],
			                           Nst_rt) );
			RC_TRY( set_N_vector_shell(elem_type, Est_st_tying_point[ii2],
			                           Nst_st) );
			RC_TRY( set_dN_matrix_shell(elem_type, Ert_rt_tying_point[ii2],
			                            dNrt_rt) );
			RC_TRY( set_dN_matrix_shell(elem_type, Ert_st_tying_point[ii2],
			                            dNrt_st) );
			RC_TRY( set_dN_matrix_shell(elem_type, Est_rt_tying_point[ii2],
			                            dNst_rt) );
			RC_TRY( set_dN_matrix_shell(elem_type, Est_st_tying_point[ii2],
			                            dNst_st) );
			RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
			                          thickness, update_location,
			                          Ert_rt_tying_point[ii2], Grt_rt, ws) );
			RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
			                          thickness, update_location,
			                          Ert_st_tying_point[ii2], Grt_st, ws) );
			RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
			                          thickness, update_location,
			                          Est_rt_tying_point[ii2], Gst_rt, ws) );
			RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
			                          thickness, update_location,
			                          Est_st_tying_point[ii2], Gst_st, ws) );

			V3G3rt_rt = inner_product3d(V3, Grt_rt[2]);
			V3G1rt_rt = inner_product3d(V3, Grt_rt[0]);
			V3G3rt_st = inner_product3d(V3, Grt_st[2]);
			V3G2rt_st = inner_product3d(V3, Grt_st[1]);
			V3G3st_rt = inner_product3d(V3, Gst_rt[2]);
			V3G1st_rt = inner_product3d(V3, Gst_rt[0]);
			V3G3st_st = inner_product3d(V3, Gst_st[2]);
			V3G2st_st = inner_product3d(V3, Gst_st[1]);

			Ert_rt = -c1*dNrt_rt[0][ii1]*V3G3rt_rt -c2*Nrt_rt[ii1]*V3G1rt_rt;
			Ert_st = -c1*dNrt_st[1][ii1]*V3G3rt_st -c2*Nrt_st[ii1]*V3G2rt_st;
			Est_rt = -c1*dNst_rt[0][ii1]*V3G3st_rt -c2*Nst_rt[ii1]*V3G1st_rt;
			Est_st = -c1*dNst_st[1][ii1]*V3G3st_st -c2*Nst_st[ii1]*V3G2st_st;

			Kex[ii1] += Ert_rt_h[ii2] * stress_rt * Ert_rt;
			Kex[ii1] += Ert_st_h[ii2] * stress_rt * Ert_st;
			Kex[ii1] += Est_rt_h[ii2] * stress_st * Est_rt;
			Kex[ii1] += Est_st_h[ii2] * stress_st * Est_st;
		}
	}

	return(NORMAL_RC);
}


/* 要素接線剛性マトリクスの有限回転項 MITC6 */
static RC
make_large_rotation_matrix_mitc6 (ELEM_TYPE elem_type, int elem_node_num,
                                  DIRECTOR_VECT elem_curr_vect[], 
                                  double thickness,
                                  double **update_location,
                                  const double *local_point,
                                  double stress_vect[],
                                  double *Kex,
                                  MITC_TYING_SET mitc, WORK_SET ws)
{
	int ii1, ii2;
	int ErrEss_tying_num;
	int Ers_tying_num;
	int ErtEst_tying_num;
	double t = local_point[2];
	double c1, c2;
	double V3G1rr;
	double V3G2ss;
	double V3G1rs_rr;
	double V3G2rs_ss;
	double V3G1rs_rs;
	double V3G2rs_rs;
	double V3G3rt_rt, V3G1rt_rt, V3G3rt_st, V3G2rt_st;
	double V3G3st_rt, V3G1st_rt, V3G3st_st, V3G2st_st;
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
	VECT3D V3;


	init_vect(elem_node_num, Kex);

	c1 = 0.5 * thickness * t;
	c2 = 0.5 * thickness;

	RC_TRY( set_tying_data_mitc6(local_point, &ErrEss_tying_num,
	                             &Ers_tying_num, &ErtEst_tying_num,
	                             Err_tying_point, Ess_tying_point,
	                             Ers_rr_tying_point, Ers_ss_tying_point,
	                             Ers_rs_tying_point, Ert_rt_tying_point,
	                             Ert_st_tying_point, Est_rt_tying_point,
	                             Est_st_tying_point, Err_h, Ess_h, Ers_rr_h,
	                             Ers_ss_h, Ers_rs_h, Ert_rt_h, Ert_st_h,
	                             Est_rt_h, Est_st_h) );

	for(ii1=0; ii1<elem_node_num; ii1++){
		double Err;
		double Ess;
		double Err_tmp;
		double Ess_tmp;
		double Ers_rr;
		double Ers_ss;
		double Ers_rs;
		double Ert_rt;
		double Ert_st;
		double Est_rt;
		double Est_st;
		double stress_rr = stress_vect[0];
		double stress_ss = stress_vect[1];
		double stress_rs = stress_vect[2];
		double stress_rt = stress_vect[3];
		double stress_st = stress_vect[4];

		V3 = elem_curr_vect[ii1].vn;

		/* Err, Ess */
		Err_tmp = 0.0;
		Ess_tmp = 0.0;
		for(ii2=0; ii2<ErrEss_tying_num; ii2++){
			RC_TRY( set_dN_matrix_shell(elem_type, Err_tying_point[ii2],
			                            dNrr) );
			RC_TRY( set_dN_matrix_shell(elem_type, Ess_tying_point[ii2],
			                            dNss) );
			RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
			                          thickness, update_location,
			                          Err_tying_point[ii2], Grr, ws) );
			RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
			                          thickness, update_location,
			                          Ess_tying_point[ii2], Gss, ws) );

			V3G1rr = inner_product3d(V3, Grr[0]);
			V3G2ss = inner_product3d(V3, Gss[1]);

			Err = -c1 * dNrr[0][ii1] * V3G1rr;
			Ess = -c1 * dNss[1][ii1] * V3G2ss;

			Kex[ii1] += Err_h[ii2] * stress_rr * Err;
			Kex[ii1] += Ess_h[ii2] * stress_ss * Ess;
			Err_tmp += Err_h[ii2] * Err;
			Ess_tmp += Ess_h[ii2] * Ess;
		}

		/* Ers */
		Kex[ii1] += stress_rs * (Err_tmp + Ess_tmp);
		for(ii2=0; ii2<Ers_tying_num; ii2++){
			RC_TRY( set_dN_matrix_shell(elem_type, Ers_rr_tying_point[ii2],
			                            dNrs_rr) );
			RC_TRY( set_dN_matrix_shell(elem_type, Ers_ss_tying_point[ii2],
			                            dNrs_ss) );
			RC_TRY( set_dN_matrix_shell(elem_type, Ers_rs_tying_point[ii2],
			                            dNrs_rs) );
			RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
			                          thickness, update_location,
			                          Ers_rr_tying_point[ii2], Grs_rr, ws) );
			RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
			                          thickness, update_location,
			                          Ers_ss_tying_point[ii2], Grs_ss, ws) );
			RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
			                          thickness, update_location,
			                          Ers_rs_tying_point[ii2], Grs_rs, ws) );

			V3G1rs_rr = inner_product3d(V3, Grs_rr[0]);
			V3G2rs_ss = inner_product3d(V3, Grs_ss[1]);
			V3G1rs_rs = inner_product3d(V3, Grs_rs[0]);
			V3G2rs_rs = inner_product3d(V3, Grs_rs[1]);

			Ers_rr = -c1 * dNrs_rr[0][ii1] * V3G1rs_rr;
			Ers_ss = -c1 * dNrs_ss[1][ii1] * V3G2rs_ss;
			Ers_rs = -c1 * (dNrs_rs[0][ii1] * V3G2rs_rs
			              + dNrs_rs[1][ii1] * V3G1rs_rs);

			Kex[ii1] += Ers_rr_h[ii2] * stress_rs * Ers_rr;
			Kex[ii1] += Ers_ss_h[ii2] * stress_rs * Ers_ss;
			Kex[ii1] += Ers_rs_h[ii2] * stress_rs * Ers_rs;
		}

		/* Ert, Est */
		for(ii2=0; ii2<ErtEst_tying_num; ii2++){
			RC_TRY( set_N_vector_shell(elem_type, Ert_rt_tying_point[ii2],
			                           Nrt_rt) );
			RC_TRY( set_N_vector_shell(elem_type, Ert_st_tying_point[ii2],
			                           Nrt_st) );
			RC_TRY( set_N_vector_shell(elem_type, Est_rt_tying_point[ii2],
			                           Nst_rt) );
			RC_TRY( set_N_vector_shell(elem_type, Est_st_tying_point[ii2],
			                           Nst_st) );
			RC_TRY( set_dN_matrix_shell(elem_type, Ert_rt_tying_point[ii2],
			                            dNrt_rt) );
			RC_TRY( set_dN_matrix_shell(elem_type, Ert_st_tying_point[ii2],
			                            dNrt_st) );
			RC_TRY( set_dN_matrix_shell(elem_type, Est_rt_tying_point[ii2],
			                            dNst_rt) );
			RC_TRY( set_dN_matrix_shell(elem_type, Est_st_tying_point[ii2],
			                            dNst_st) );
			RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
			                          thickness, update_location,
			                          Ert_rt_tying_point[ii2], Grt_rt, ws) );
			RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
			                          thickness, update_location,
			                          Ert_st_tying_point[ii2], Grt_st, ws) );
			RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
			                          thickness, update_location,
			                          Est_rt_tying_point[ii2], Gst_rt, ws) );
			RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
			                          thickness, update_location,
			                          Est_st_tying_point[ii2], Gst_st, ws) );

			V3G3rt_rt = inner_product3d(V3, Grt_rt[2]);
			V3G1rt_rt = inner_product3d(V3, Grt_rt[0]);
			V3G3rt_st = inner_product3d(V3, Grt_st[2]);
			V3G2rt_st = inner_product3d(V3, Grt_st[1]);
			V3G3st_rt = inner_product3d(V3, Gst_rt[2]);
			V3G1st_rt = inner_product3d(V3, Gst_rt[0]);
			V3G3st_st = inner_product3d(V3, Gst_st[2]);
			V3G2st_st = inner_product3d(V3, Gst_st[1]);

			Ert_rt = -c1*dNrt_rt[0][ii1]*V3G3rt_rt -c2*Nrt_rt[ii1]*V3G1rt_rt;
			Ert_st = -c1*dNrt_st[1][ii1]*V3G3rt_st -c2*Nrt_st[ii1]*V3G2rt_st;
			Est_rt = -c1*dNst_rt[0][ii1]*V3G3st_rt -c2*Nst_rt[ii1]*V3G1st_rt;
			Est_st = -c1*dNst_st[1][ii1]*V3G3st_st -c2*Nst_st[ii1]*V3G2st_st;

			Kex[ii1] += Ert_rt_h[ii2] * stress_rt * Ert_rt;
			Kex[ii1] += Ert_st_h[ii2] * stress_rt * Ert_st;
			Kex[ii1] += Est_rt_h[ii2] * stress_st * Est_rt;
			Kex[ii1] += Est_st_h[ii2] * stress_st * Est_st;
		}
	}

	return(NORMAL_RC);
}


/* 要素接線剛性マトリクスの有限回転項 MITC4 */
static RC
make_large_rotation_matrix_mitc4 (ELEM_TYPE elem_type, int elem_node_num,
                                  DIRECTOR_VECT elem_curr_vect[], 
                                  double thickness,
                                  double **update_location,
                                  const double *local_point,
                                  double stress_vect[],
                                  double *Kex,
                                  MITC_TYING_SET mitc, WORK_SET ws)
{
	int ii1, ii2;
	int ErtEst_tying_num;
	double t = local_point[2];
	double c1, c2;
	double V3G1, V3G2;
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
	VECT3D V3;


	init_vect(elem_node_num, Kex);

	c1 = 0.5 * thickness * t;
	c2 = 0.5 * thickness;

	RC_TRY( set_dN_matrix_shell(elem_type, local_point, dN) );
	RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
	                          thickness, update_location, local_point, G, ws) );

	RC_TRY( set_tying_data_mitc4(local_point, &ErtEst_tying_num,
	                             Ert_tying_point, Est_tying_point,
	                             Ert_h, Est_h) );

	for(ii1=0; ii1<elem_node_num; ii1++){
		double V3G3rt, V3G1rt;
		double V3G3st, V3G2st;
		double Ert, Est;
		double stress_rr = stress_vect[0];
		double stress_ss = stress_vect[1];
		double stress_rs = stress_vect[2];
		double stress_rt = stress_vect[3];
		double stress_st = stress_vect[4];

		V3 = elem_curr_vect[ii1].vn;

		V3G1 = inner_product3d(V3, G[0]);
		V3G2 = inner_product3d(V3, G[1]);

		Kex[ii1] = - stress_rr*c1*dN[0][ii1]*V3G1
		           - stress_ss*c1*dN[1][ii1]*V3G2
		           - stress_rs*c1*(dN[0][ii1]*V3G2 + dN[1][ii1]*V3G1);

		for(ii2=0; ii2<ErtEst_tying_num; ii2++){
			RC_TRY( set_N_vector_shell(elem_type, Ert_tying_point[ii2], Nrt) );
			RC_TRY( set_N_vector_shell(elem_type, Est_tying_point[ii2], Nst) );
			RC_TRY( set_dN_matrix_shell(elem_type, Ert_tying_point[ii2],
			                            dNrt) );
			RC_TRY( set_dN_matrix_shell(elem_type, Est_tying_point[ii2],
			                            dNst) );
			RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
			                          thickness, update_location,
			                          Ert_tying_point[ii2], Grt, ws) );
			RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
			                          thickness, update_location,
			                          Est_tying_point[ii2], Gst, ws) );

			V3G3rt = inner_product3d(V3, Grt[2]);
			V3G1rt = inner_product3d(V3, Grt[0]);
			V3G3st = inner_product3d(V3, Gst[2]);
			V3G2st = inner_product3d(V3, Gst[1]);

			Ert = -c1*dNrt[0][ii1]*V3G3rt -c2*Nrt[ii1]*V3G1rt;
			Est = -c1*dNst[1][ii1]*V3G3st -c2*Nst[ii1]*V3G2st;

			Kex[ii1] += Ert_h[ii2] * stress_rt * Ert;
			Kex[ii1] += Est_h[ii2] * stress_st * Est;
		}
	}


	return(NORMAL_RC);
}


/* 要素接線剛性マトリクスの有限回転項 MITC8 */
static RC
make_large_rotation_matrix_mitc8 (ELEM_TYPE elem_type, int elem_node_num,
                                  DIRECTOR_VECT elem_curr_vect[],
                                  double thickness,
                                  double **update_location,
                                  const double *local_point,
                                  double stress_vect[],
                                  double *Kex,
                                  MITC_TYING_SET mitc, WORK_SET ws)
{
	int ii1, ii2;
	int Ers_tying_num;
	int ErrEssErtEst_tying_num;
	double t = local_point[2];
	double c1, c2;
	double V3G1rr;
	double V3G2ss;
	double V3G3rt;
	double V3G1rt;
	double V3G3st;
	double V3G2st;
	double V3G1rs;
	double V3G2rs;
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
	VECT3D V3;


	init_vect(elem_node_num, Kex);

	c1 = 0.5 * thickness * t;
	c2 = 0.5 * thickness;

	RC_TRY( set_tying_data_mitc8(local_point, &Ers_tying_num,
	                             &ErrEssErtEst_tying_num, Err_tying_point,
	                             Ess_tying_point, Ers_tying_point,
	                             Ert_tying_point, Est_tying_point, Err_h,
	                             Ess_h, Ers_h, Ert_h, Est_h) );

	for(ii1=0; ii1<elem_node_num; ii1++){
		double Err;
		double Ess;
		double Ers;
		double Ert;
		double Est;
		double stress_rr = stress_vect[0];
		double stress_ss = stress_vect[1];
		double stress_rs = stress_vect[2];
		double stress_rt = stress_vect[3];
		double stress_st = stress_vect[4];

		V3 = elem_curr_vect[ii1].vn;


		for(ii2=0; ii2<Ers_tying_num; ii2++){
			RC_TRY( set_dN_matrix_shell(elem_type, Ers_tying_point[ii2],
			                            dNrs) );
			RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
			                          thickness, update_location,
			                          Ers_tying_point[ii2], Grs, ws) );

			V3G1rs = inner_product3d(V3, Grs[0]);
			V3G2rs = inner_product3d(V3, Grs[1]);

			Ers = -c1 * (dNrs[0][ii1] * V3G2rs + dNrs[1][ii1] * V3G1rs);

			Kex[ii1] += Ers_h[ii2] * stress_rs * Ers;
		}

		for(ii2=0; ii2<ErrEssErtEst_tying_num; ii2++){
			RC_TRY( set_dN_matrix_shell(elem_type, Err_tying_point[ii2],
			                            dNrr) );
			RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
			                          thickness, update_location,
			                          Err_tying_point[ii2], Grr, ws) );
			RC_TRY( set_dN_matrix_shell(elem_type, Ess_tying_point[ii2],
			                            dNss) );
			RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
			                          thickness, update_location,
			                          Ess_tying_point[ii2], Gss, ws) );
			RC_TRY( set_N_vector_shell(elem_type, Ert_tying_point[ii2], Nrt) );
			RC_TRY( set_dN_matrix_shell(elem_type, Ert_tying_point[ii2],
			                            dNrt) );
			RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
			                          thickness, update_location,
			                          Ert_tying_point[ii2], Grt, ws) );
			RC_TRY( set_N_vector_shell(elem_type, Est_tying_point[ii2], Nst) );
			RC_TRY( set_dN_matrix_shell(elem_type, Est_tying_point[ii2],
			                            dNst) );
			RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
			                          thickness, update_location,
			                          Est_tying_point[ii2], Gst, ws) );

			V3G1rr = inner_product3d(V3, Grr[0]);
			V3G2ss = inner_product3d(V3, Gss[1]);
			V3G3rt = inner_product3d(V3, Grt[2]);
			V3G1rt = inner_product3d(V3, Grt[0]);
			V3G3st = inner_product3d(V3, Gst[2]);
			V3G2st = inner_product3d(V3, Gst[1]);

			Err = -c1 * dNrr[0][ii1] * V3G1rr;
			Ess = -c1 * dNss[1][ii1] * V3G2ss;
			Ert = -c1*dNrt[0][ii1]*V3G3rt -c2*Nrt[ii1]*V3G1rt;
			Est = -c1*dNst[1][ii1]*V3G3st -c2*Nst[ii1]*V3G2st;

			Kex[ii1] += Err_h[ii2] * stress_rr * Err;
			Kex[ii1] += Ess_h[ii2] * stress_ss * Ess;
			Kex[ii1] += Ert_h[ii2] * stress_rt * Ert;
			Kex[ii1] += Est_h[ii2] * stress_st * Est;
		}
	}

	return(NORMAL_RC);
}


/* 初期応力マトリクス */
RC
make_SZtZ_matrix (ELEM_TYPE elem_type, int elem_node_num,
                  DIRECTOR_VECT elem_curr_vect[], double thickness,
                  const double *local_point, double stress_vect[],
                  double **SZtZ_matrix, MITC_TYING_SET mitc, WORK_SET ws)
{
	switch(elem_type){
	case ELEM_TRI1:
		RC_TRY( make_SZtZ_matrix_mitc3(elem_type, elem_node_num,
		                               elem_curr_vect, thickness, local_point,
		                               stress_vect, SZtZ_matrix, mitc, ws) );
		break;
	case ELEM_TRI2:
		RC_TRY( make_SZtZ_matrix_mitc6(elem_type, elem_node_num,
		                               elem_curr_vect, thickness, local_point,
		                               stress_vect, SZtZ_matrix, mitc, ws) );
		break;
	case ELEM_QUAD1:
		RC_TRY( make_SZtZ_matrix_mitc4(elem_type, elem_node_num,
		                               elem_curr_vect, thickness, local_point,
		                               stress_vect, SZtZ_matrix, mitc, ws) );
		break;
	case ELEM_QUAD2:
		RC_TRY( make_SZtZ_matrix_mitc8(elem_type, elem_node_num,
		                               elem_curr_vect, thickness, local_point,
		                               stress_vect, SZtZ_matrix, mitc, ws) );
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


/* 初期応力マトリクス MITC3 */
static RC
make_SZtZ_matrix_mitc3 (ELEM_TYPE elem_type, int elem_node_num,
                        DIRECTOR_VECT elem_curr_vect[], double thickness,
                        const double *local_point, double stress_vect[],
                        double **SZtZ_matrix, MITC_TYING_SET mitc, WORK_SET ws)
{
	int ii1, ii2, ii3, ii4;
	int dim = 5;
	int elem_dof;
	int ErtEst_tying_num;
	double *Ert_rt_h = mitc.vec_T[0];
	double *Ert_st_h = mitc.vec_T[1];
	double *Est_rt_h = mitc.vec_T[2];
	double *Est_st_h = mitc.vec_T[3];
	double **Ert_rt_tying_point = mitc.mat_T_3[0];
	double **Ert_st_tying_point = mitc.mat_T_3[1];
	double **Est_rt_tying_point = mitc.mat_T_3[2];
	double **Est_st_tying_point = mitc.mat_T_3[3];
	double **Z1 = NULL;
	double **Z2 = NULL;
	double **Z1rt_rt = NULL;
	double **Z1st_rt = NULL;
	double **Z2rt_st = NULL;
	double **Z2st_st = NULL;
	double **Z3rt_rt = NULL;
	double **Z3rt_st = NULL;
	double **Z3st_rt = NULL;
	double **Z3st_st = NULL;
	double **ZtZ11 = mitc.mat_5N_5N[0];
	double **ZtZ22 = mitc.mat_5N_5N[1];
	double **ZtZ12 = mitc.mat_5N_5N[2];
	double **ZtZ13rt_rt = mitc.mat_5N_5N[3];
	double **ZtZ13rt_st = mitc.mat_5N_5N[4];
	double **ZtZ23st_rt = mitc.mat_5N_5N[5];
	double **ZtZ23st_st = mitc.mat_5N_5N[6];


	elem_dof = dim * elem_node_num;

	init_matrix(elem_dof, elem_dof, SZtZ_matrix);
	init_matrix(elem_dof, elem_dof, ZtZ11);
	init_matrix(elem_dof, elem_dof, ZtZ22);
	init_matrix(elem_dof, elem_dof, ZtZ12);

	Z1 = mitc.mat_3_5N[0];
	Z2 = mitc.mat_3_5N[1];

	RC_TRY( set_Z1_matrix(elem_type, elem_node_num, elem_curr_vect, thickness,
	                      local_point, Z1, ws) );
	RC_TRY( set_Z2_matrix(elem_type, elem_node_num, elem_curr_vect, thickness,
	                      local_point, Z2, ws) );

	for(ii1=0; ii1<3; ii1++){
		for(ii2=0; ii2<elem_dof; ii2++){
			for(ii3=0; ii3<elem_dof; ii3++){
				ZtZ11[ii2][ii3] += Z1[ii1][ii2] * Z1[ii1][ii3];
				ZtZ22[ii2][ii3] += Z2[ii1][ii2] * Z2[ii1][ii3];
				ZtZ12[ii2][ii3] += Z1[ii1][ii2] * Z2[ii1][ii3]
				                 + Z2[ii1][ii2] * Z1[ii1][ii3];
			}
		}
	}
	for(ii1=0; ii1<elem_dof; ii1++){
		for(ii2=0; ii2<elem_dof; ii2++){
			SZtZ_matrix[ii1][ii2] += stress_vect[0]*ZtZ11[ii1][ii2]
			                       + stress_vect[1]*ZtZ22[ii1][ii2]
			                       + stress_vect[2]*ZtZ12[ii1][ii2];
		}
	}

	RC_TRY( set_tying_data_mitc3(local_point, &ErtEst_tying_num,
	                             Ert_rt_tying_point, Ert_st_tying_point,
	                             Est_rt_tying_point, Est_st_tying_point,
	                             Ert_rt_h, Ert_st_h, Est_rt_h, Est_st_h) );

	Z1rt_rt = mitc.mat_3_5N[0];
	Z1st_rt = mitc.mat_3_5N[1];
	Z2rt_st = mitc.mat_3_5N[2];
	Z2st_st = mitc.mat_3_5N[3];
	Z3rt_rt = mitc.mat_3_5N[4];
	Z3rt_st = mitc.mat_3_5N[5];
	Z3st_rt = mitc.mat_3_5N[6];
	Z3st_st = mitc.mat_3_5N[7];
	for(ii1=0; ii1<ErtEst_tying_num; ii1++){
		double stress_rt = stress_vect[3];
		double stress_st = stress_vect[4];

		RC_TRY( set_Z1_matrix(elem_type, elem_node_num, elem_curr_vect,
		                     thickness, Ert_rt_tying_point[ii1], Z1rt_rt, ws) );
		RC_TRY( set_Z1_matrix(elem_type, elem_node_num, elem_curr_vect,
		                     thickness, Est_rt_tying_point[ii1], Z1st_rt, ws) );
		RC_TRY( set_Z2_matrix(elem_type, elem_node_num, elem_curr_vect,
		                     thickness, Ert_st_tying_point[ii1], Z2rt_st, ws) );
		RC_TRY( set_Z2_matrix(elem_type, elem_node_num, elem_curr_vect,
		                     thickness, Est_st_tying_point[ii1], Z2st_st, ws) );
		RC_TRY( set_Z3_matrix(elem_type, elem_node_num, elem_curr_vect,
		                     thickness, Ert_rt_tying_point[ii1], Z3rt_rt, ws) );
		RC_TRY( set_Z3_matrix(elem_type, elem_node_num, elem_curr_vect,
		                     thickness, Ert_st_tying_point[ii1], Z3rt_st, ws) );
		RC_TRY( set_Z3_matrix(elem_type, elem_node_num, elem_curr_vect,
		                     thickness, Est_rt_tying_point[ii1], Z3st_rt, ws) );
		RC_TRY( set_Z3_matrix(elem_type, elem_node_num, elem_curr_vect,
		                     thickness, Est_st_tying_point[ii1], Z3st_st, ws) );

		/* ZtZ13 = Z1tZ3 + Z3tZ1 */
		/* ZtZ23 = Z2tZ3 + Z3tZ2 */
		init_matrix(elem_dof, elem_dof, ZtZ13rt_rt);
		init_matrix(elem_dof, elem_dof, ZtZ13rt_st);
		init_matrix(elem_dof, elem_dof, ZtZ23st_rt);
		init_matrix(elem_dof, elem_dof, ZtZ23st_st);
		for(ii2=0; ii2<3; ii2++){
			for(ii3=0; ii3<elem_dof; ii3++){
				for(ii4=0; ii4<elem_dof; ii4++){
					ZtZ13rt_rt[ii3][ii4] += Z1rt_rt[ii2][ii3]*Z3rt_rt[ii2][ii4]
					                      + Z3rt_rt[ii2][ii3]*Z1rt_rt[ii2][ii4];
					ZtZ13rt_st[ii3][ii4] += Z2rt_st[ii2][ii3]*Z3rt_st[ii2][ii4]
					                      + Z3rt_st[ii2][ii3]*Z2rt_st[ii2][ii4];
					ZtZ23st_rt[ii3][ii4] += Z1st_rt[ii2][ii3]*Z3st_rt[ii2][ii4]
					                      + Z3st_rt[ii2][ii3]*Z1st_rt[ii2][ii4];
					ZtZ23st_st[ii3][ii4] += Z2st_st[ii2][ii3]*Z3st_st[ii2][ii4]
					                      + Z3st_st[ii2][ii3]*Z2st_st[ii2][ii4];
				}
			}
		}

		for(ii2=0; ii2<elem_dof; ii2++){
			for(ii3=0; ii3<elem_dof; ii3++){
				SZtZ_matrix[ii2][ii3]
				+= Ert_rt_h[ii1]*stress_rt*ZtZ13rt_rt[ii2][ii3];
				SZtZ_matrix[ii2][ii3]
				+= Ert_st_h[ii1]*stress_rt*ZtZ13rt_st[ii2][ii3];
				SZtZ_matrix[ii2][ii3]
				+= Est_rt_h[ii1]*stress_st*ZtZ23st_rt[ii2][ii3];
				SZtZ_matrix[ii2][ii3]
				+= Est_st_h[ii1]*stress_st*ZtZ23st_st[ii2][ii3];
			}
		}
	}

	return(NORMAL_RC);
}


/* 初期応力マトリクス MITC6 */
static RC
make_SZtZ_matrix_mitc6 (ELEM_TYPE elem_type, int elem_node_num,
                        DIRECTOR_VECT elem_curr_vect[], double thickness,
                        const double *local_point, double stress_vect[],
                        double **SZtZ_matrix, MITC_TYING_SET mitc, WORK_SET ws)
{
	int ii1, ii2, ii3, ii4;
	int dim = 5;
	int elem_dof;
	int ErrEss_tying_num;
	int Ers_tying_num;
	int ErtEst_tying_num;
	double *Err_h = mitc.vec_T[0];
	double *Ess_h = mitc.vec_T[1];
	double *Ers_rr_h = mitc.vec_T[2];
	double *Ers_ss_h = mitc.vec_T[3];
	double *Ers_rs_h = mitc.vec_T[4];
	double *Ert_rt_h = mitc.vec_T[5];
	double *Ert_st_h = mitc.vec_T[6];
	double *Est_rt_h = mitc.vec_T[7];
	double *Est_st_h = mitc.vec_T[8];
	double **Err_tying_point = mitc.mat_T_3[0];
	double **Ess_tying_point = mitc.mat_T_3[1];
	double **Ers_rr_tying_point = mitc.mat_T_3[2];
	double **Ers_ss_tying_point = mitc.mat_T_3[3];
	double **Ers_rs_tying_point = mitc.mat_T_3[4];
	double **Ert_rt_tying_point = mitc.mat_T_3[5];
	double **Ert_st_tying_point = mitc.mat_T_3[6];
	double **Est_rt_tying_point = mitc.mat_T_3[7];
	double **Est_st_tying_point = mitc.mat_T_3[8];
	double **Z1rr = NULL;
	double **Z2ss = NULL;
	double **Z1rs_rr = NULL;
	double **Z2rs_ss = NULL;
	double **Z1rs_rs = NULL;
	double **Z2rs_rs = NULL;
	double **Z1rt_rt = NULL;
	double **Z1st_rt = NULL;
	double **Z2rt_st = NULL;
	double **Z2st_st = NULL;
	double **Z3rt_rt = NULL;
	double **Z3rt_st = NULL;
	double **Z3st_rt = NULL;
	double **Z3st_st = NULL;
	double **ZtZ11rr = NULL;
	double **ZtZ22ss = NULL;
	double **ZtZ11rr_tmp = NULL;
	double **ZtZ22ss_tmp = NULL;
	double **ZtZ11rs_rr = NULL;
	double **ZtZ22rs_ss = NULL;
	double **ZtZ12rs_rs = NULL;
	double **ZtZ13rt_rt = NULL;
	double **ZtZ13rt_st = NULL;
	double **ZtZ23st_rt = NULL;
	double **ZtZ23st_st = NULL;


	elem_dof = dim * elem_node_num;

	init_matrix(elem_dof, elem_dof, SZtZ_matrix);

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
	Z1rr = mitc.mat_3_5N[0];
	Z2ss = mitc.mat_3_5N[1];
	ZtZ11rr = mitc.mat_5N_5N[0];
	ZtZ22ss = mitc.mat_5N_5N[1];
	ZtZ11rr_tmp = mitc.mat_5N_5N[2];
	ZtZ22ss_tmp = mitc.mat_5N_5N[3];
	init_matrix(elem_dof, elem_dof, ZtZ11rr_tmp);
	init_matrix(elem_dof, elem_dof, ZtZ22ss_tmp);
	for(ii1=0; ii1<ErrEss_tying_num; ii1++){
		double stress_rr = stress_vect[0];
		double stress_ss = stress_vect[1];

		RC_TRY( set_Z1_matrix(elem_type, elem_node_num, elem_curr_vect,
		                      thickness, Err_tying_point[ii1], Z1rr, ws) );
		RC_TRY( set_Z2_matrix(elem_type, elem_node_num, elem_curr_vect,
		                      thickness, Ess_tying_point[ii1], Z2ss, ws) );

		init_matrix(elem_dof, elem_dof, ZtZ11rr);
		init_matrix(elem_dof, elem_dof, ZtZ22ss);
		for(ii2=0; ii2<3; ii2++){
			for(ii3=0; ii3<elem_dof; ii3++){
				for(ii4=0; ii4<elem_dof; ii4++){
					ZtZ11rr[ii3][ii4] += Z1rr[ii2][ii3] * Z1rr[ii2][ii4];
					ZtZ22ss[ii3][ii4] += Z2ss[ii2][ii3] * Z2ss[ii2][ii4];
				}
			}
		}

		for(ii2=0; ii2<elem_dof; ii2++){
			for(ii3=0; ii3<elem_dof; ii3++){
				SZtZ_matrix[ii2][ii3]+=Err_h[ii1]*stress_rr*ZtZ11rr[ii2][ii3];
				SZtZ_matrix[ii2][ii3]+=Ess_h[ii1]*stress_ss*ZtZ22ss[ii2][ii3];
				ZtZ11rr_tmp[ii2][ii3] += Err_h[ii1] * ZtZ11rr[ii2][ii3];
				ZtZ22ss_tmp[ii2][ii3] += Ess_h[ii1] * ZtZ22ss[ii2][ii3];
			}
		}
	}

	/* Ers */
	for(ii1=0; ii1<elem_dof; ii1++){
		for(ii2=0; ii2<elem_dof; ii2++){
			SZtZ_matrix[ii1][ii2]
			+= stress_vect[2] * (ZtZ11rr_tmp[ii1][ii2] + ZtZ22ss_tmp[ii1][ii2]);
		}
	}
	Z1rs_rr = mitc.mat_3_5N[0];
	Z2rs_ss = mitc.mat_3_5N[1];
	Z1rs_rs = mitc.mat_3_5N[2];
	Z2rs_rs = mitc.mat_3_5N[3];
	ZtZ11rs_rr = mitc.mat_5N_5N[0];
	ZtZ22rs_ss = mitc.mat_5N_5N[1];
	ZtZ12rs_rs = mitc.mat_5N_5N[2];
	for(ii1=0; ii1<Ers_tying_num; ii1++){
		double stress_rs = stress_vect[2];

		RC_TRY( set_Z1_matrix(elem_type, elem_node_num, elem_curr_vect,
		                      thickness, Ers_rr_tying_point[ii1],
		                      Z1rs_rr, ws) );
		RC_TRY( set_Z2_matrix(elem_type, elem_node_num, elem_curr_vect,
		                      thickness, Ers_ss_tying_point[ii1],
		                      Z2rs_ss, ws) );
		RC_TRY( set_Z1_matrix(elem_type, elem_node_num, elem_curr_vect,
		                      thickness, Ers_rs_tying_point[ii1],
		                      Z1rs_rs, ws) );
		RC_TRY( set_Z2_matrix(elem_type, elem_node_num, elem_curr_vect,
		                      thickness, Ers_rs_tying_point[ii1],
		                      Z2rs_rs, ws) );

		init_matrix(elem_dof, elem_dof, ZtZ11rs_rr);
		init_matrix(elem_dof, elem_dof, ZtZ22rs_ss);
		init_matrix(elem_dof, elem_dof, ZtZ12rs_rs);
		for(ii2=0; ii2<3; ii2++){
			for(ii3=0; ii3<elem_dof; ii3++){
				for(ii4=0; ii4<elem_dof; ii4++){
					ZtZ11rs_rr[ii3][ii4]+=Z1rs_rr[ii2][ii3]*Z1rs_rr[ii2][ii4];
					ZtZ22rs_ss[ii3][ii4]+=Z2rs_ss[ii2][ii3]*Z2rs_ss[ii2][ii4];
					ZtZ12rs_rs[ii3][ii4]+=Z1rs_rs[ii2][ii3]*Z2rs_rs[ii2][ii4]
					                     +Z2rs_rs[ii2][ii3]*Z1rs_rs[ii2][ii4];
				}
			}
		}
		for(ii2=0; ii2<elem_dof; ii2++){
			for(ii3=0; ii3<elem_dof; ii3++){
				SZtZ_matrix[ii2][ii3]
				+= Ers_rr_h[ii1] * stress_rs * ZtZ11rs_rr[ii2][ii3];
				SZtZ_matrix[ii2][ii3]
				+= Ers_ss_h[ii1] * stress_rs * ZtZ22rs_ss[ii2][ii3];
				SZtZ_matrix[ii2][ii3]
				+= Ers_rs_h[ii1] * stress_rs * ZtZ12rs_rs[ii2][ii3];
			}
		}
	}

	/* Ert, Est */
	Z1rt_rt = mitc.mat_3_5N[0];
	Z1st_rt = mitc.mat_3_5N[1];
	Z2rt_st = mitc.mat_3_5N[2];
	Z2st_st = mitc.mat_3_5N[3];
	Z3rt_rt = mitc.mat_3_5N[4];
	Z3rt_st = mitc.mat_3_5N[5];
	Z3st_rt = mitc.mat_3_5N[6];
	Z3st_st = mitc.mat_3_5N[7];
	ZtZ13rt_rt = mitc.mat_5N_5N[0];
	ZtZ13rt_st = mitc.mat_5N_5N[1];
	ZtZ23st_rt = mitc.mat_5N_5N[2];
	ZtZ23st_st = mitc.mat_5N_5N[3];
	for(ii1=0; ii1<ErtEst_tying_num; ii1++){
		double stress_rt = stress_vect[3];
		double stress_st = stress_vect[4];

		RC_TRY( set_Z1_matrix(elem_type, elem_node_num, elem_curr_vect,
		                     thickness, Ert_rt_tying_point[ii1], Z1rt_rt, ws) );
		RC_TRY( set_Z1_matrix(elem_type, elem_node_num, elem_curr_vect,
		                     thickness, Est_rt_tying_point[ii1], Z1st_rt, ws) );
		RC_TRY( set_Z2_matrix(elem_type, elem_node_num, elem_curr_vect,
		                     thickness, Ert_st_tying_point[ii1], Z2rt_st, ws) );
		RC_TRY( set_Z2_matrix(elem_type, elem_node_num, elem_curr_vect,
		                     thickness, Est_st_tying_point[ii1], Z2st_st, ws) );
		RC_TRY( set_Z3_matrix(elem_type, elem_node_num, elem_curr_vect,
		                     thickness, Ert_rt_tying_point[ii1], Z3rt_rt, ws) );
		RC_TRY( set_Z3_matrix(elem_type, elem_node_num, elem_curr_vect,
		                     thickness, Ert_st_tying_point[ii1], Z3rt_st, ws) );
		RC_TRY( set_Z3_matrix(elem_type, elem_node_num, elem_curr_vect,
		                     thickness, Est_rt_tying_point[ii1], Z3st_rt, ws) );
		RC_TRY( set_Z3_matrix(elem_type, elem_node_num, elem_curr_vect,
		                     thickness, Est_st_tying_point[ii1], Z3st_st, ws) );

		/* ZtZ13 = Z1tZ3 + Z3tZ1 */
		/* ZtZ23 = Z2tZ3 + Z3tZ2 */
		init_matrix(elem_dof, elem_dof, ZtZ13rt_rt);
		init_matrix(elem_dof, elem_dof, ZtZ13rt_st);
		init_matrix(elem_dof, elem_dof, ZtZ23st_rt);
		init_matrix(elem_dof, elem_dof, ZtZ23st_st);
		for(ii2=0; ii2<3; ii2++){
			for(ii3=0; ii3<elem_dof; ii3++){
				for(ii4=0; ii4<elem_dof; ii4++){
					ZtZ13rt_rt[ii3][ii4] += Z1rt_rt[ii2][ii3]*Z3rt_rt[ii2][ii4]
					                      + Z3rt_rt[ii2][ii3]*Z1rt_rt[ii2][ii4];
					ZtZ13rt_st[ii3][ii4] += Z2rt_st[ii2][ii3]*Z3rt_st[ii2][ii4]
					                      + Z3rt_st[ii2][ii3]*Z2rt_st[ii2][ii4];
					ZtZ23st_rt[ii3][ii4] += Z1st_rt[ii2][ii3]*Z3st_rt[ii2][ii4]
					                      + Z3st_rt[ii2][ii3]*Z1st_rt[ii2][ii4];
					ZtZ23st_st[ii3][ii4] += Z2st_st[ii2][ii3]*Z3st_st[ii2][ii4]
					                      + Z3st_st[ii2][ii3]*Z2st_st[ii2][ii4];
				}
			}
		}

		for(ii2=0; ii2<elem_dof; ii2++){
			for(ii3=0; ii3<elem_dof; ii3++){
				SZtZ_matrix[ii2][ii3]
				+= Ert_rt_h[ii1]*stress_rt*ZtZ13rt_rt[ii2][ii3];
				SZtZ_matrix[ii2][ii3]
				+= Ert_st_h[ii1]*stress_rt*ZtZ13rt_st[ii2][ii3];
				SZtZ_matrix[ii2][ii3]
				+= Est_rt_h[ii1]*stress_st*ZtZ23st_rt[ii2][ii3];
				SZtZ_matrix[ii2][ii3]
				+= Est_st_h[ii1]*stress_st*ZtZ23st_st[ii2][ii3];
			}
		}
	}

	return(NORMAL_RC);
}


/* 初期応力マトリクス MITC4 */
static RC
make_SZtZ_matrix_mitc4 (ELEM_TYPE elem_type, int elem_node_num,
                        DIRECTOR_VECT elem_curr_vect[], double thickness,
                        const double *local_point, double stress_vect[],
                        double **SZtZ_matrix, MITC_TYING_SET mitc, WORK_SET ws)
{
	int ii1, ii2, ii3, ii4;
	int dim = 5;
	int elem_dof;
	int ErtEst_tying_num;
	double *Ert_h = mitc.vec_T[0];
	double *Est_h = mitc.vec_T[1];
	double **Ert_tying_point = mitc.mat_T_3[0];
	double **Est_tying_point = mitc.mat_T_3[1];
	double **Z1 = mitc.mat_3_5N[0];
	double **Z2 = mitc.mat_3_5N[1];
	double **Z1rt = mitc.mat_3_5N[2];
	double **Z2st = mitc.mat_3_5N[3];
	double **Z3rt = mitc.mat_3_5N[4];
	double **Z3st = mitc.mat_3_5N[5];
	double **ZtZ11 = mitc.mat_5N_5N[0];
	double **ZtZ22 = mitc.mat_5N_5N[1];
	double **ZtZ12 = mitc.mat_5N_5N[2];
	double **ZtZ13 = mitc.mat_5N_5N[3];
	double **ZtZ23 = mitc.mat_5N_5N[4];


	elem_dof = dim * elem_node_num;

	init_matrix(elem_dof, elem_dof, SZtZ_matrix);
	init_matrix(elem_dof, elem_dof, ZtZ11);
	init_matrix(elem_dof, elem_dof, ZtZ22);
	init_matrix(elem_dof, elem_dof, ZtZ12);

	RC_TRY( set_Z1_matrix(elem_type, elem_node_num, elem_curr_vect, thickness,
	                      local_point, Z1, ws) );
	RC_TRY( set_Z2_matrix(elem_type, elem_node_num, elem_curr_vect, thickness,
	                      local_point, Z2, ws) );

	for(ii1=0; ii1<3; ii1++){
		for(ii2=0; ii2<elem_dof; ii2++){
			for(ii3=0; ii3<elem_dof; ii3++){
				ZtZ11[ii2][ii3] += Z1[ii1][ii2] * Z1[ii1][ii3];
				ZtZ22[ii2][ii3] += Z2[ii1][ii2] * Z2[ii1][ii3];
				ZtZ12[ii2][ii3] += Z1[ii1][ii2] * Z2[ii1][ii3]
				                 + Z2[ii1][ii2] * Z1[ii1][ii3];
			}
		}
	}
	for(ii1=0; ii1<elem_dof; ii1++){
		for(ii2=0; ii2<elem_dof; ii2++){
			SZtZ_matrix[ii1][ii2] += stress_vect[0]*ZtZ11[ii1][ii2]
			                       + stress_vect[1]*ZtZ22[ii1][ii2]
			                       + stress_vect[2]*ZtZ12[ii1][ii2];
		}
	}

	RC_TRY( set_tying_data_mitc4(local_point, &ErtEst_tying_num,
	                             Ert_tying_point, Est_tying_point,
	                             Ert_h, Est_h) );

	for(ii1=0; ii1<ErtEst_tying_num; ii1++){
		double stress_rt = stress_vect[3];
		double stress_st = stress_vect[4];

		RC_TRY( set_Z1_matrix(elem_type, elem_node_num, elem_curr_vect,
		                      thickness, Ert_tying_point[ii1], Z1rt, ws) );
		RC_TRY( set_Z3_matrix(elem_type, elem_node_num, elem_curr_vect,
		                      thickness, Ert_tying_point[ii1], Z3rt, ws) );
		RC_TRY( set_Z2_matrix(elem_type, elem_node_num, elem_curr_vect,
		                      thickness, Est_tying_point[ii1], Z2st, ws) );
		RC_TRY( set_Z3_matrix(elem_type, elem_node_num, elem_curr_vect,
		                      thickness, Est_tying_point[ii1], Z3st, ws) );

		/* ZtZ13 = Z1tZ3 + Z3tZ1 */
		/* ZtZ23 = Z2tZ3 + Z3tZ2 */
		init_matrix(elem_dof, elem_dof, ZtZ13);
		init_matrix(elem_dof, elem_dof, ZtZ23);
		for(ii2=0; ii2<3; ii2++){
			for(ii3=0; ii3<elem_dof; ii3++){
				for(ii4=0; ii4<elem_dof; ii4++){
					ZtZ13[ii3][ii4] += Z1rt[ii2][ii3] * Z3rt[ii2][ii4]
					                 + Z3rt[ii2][ii3] * Z1rt[ii2][ii4];
					ZtZ23[ii3][ii4] += Z2st[ii2][ii3] * Z3st[ii2][ii4]
					                 + Z3st[ii2][ii3] * Z2st[ii2][ii4];
				}
			}
		}

		for(ii2=0; ii2<elem_dof; ii2++){
			for(ii3=0; ii3<elem_dof; ii3++){
				SZtZ_matrix[ii2][ii3] += Ert_h[ii1]*stress_rt*ZtZ13[ii2][ii3];
				SZtZ_matrix[ii2][ii3] += Est_h[ii1]*stress_st*ZtZ23[ii2][ii3];
			}
		}
	}

	return(NORMAL_RC);
}


/* 初期応力マトリクス MITC8 */
static RC
make_SZtZ_matrix_mitc8 (ELEM_TYPE elem_type, int elem_node_num,
                        DIRECTOR_VECT elem_curr_vect[], double thickness,
                        const double *local_point, double stress_vect[],
                        double **SZtZ_matrix,
                        MITC_TYING_SET mitc, WORK_SET ws)
{
	int ii1, ii2, ii3, ii4;
	int dim = 5;
	int elem_dof;
	int Ers_tying_num;
	int ErrEssErtEst_tying_num;
	double *Err_h = mitc.vec_T[0];
	double *Ess_h = mitc.vec_T[1];
	double *Ers_h = mitc.vec_T[2];
	double *Ert_h = mitc.vec_T[3];
	double *Est_h = mitc.vec_T[4];
	double **Err_tying_point = mitc.mat_T_3[0];
	double **Ess_tying_point = mitc.mat_T_3[1];
	double **Ers_tying_point = mitc.mat_T_3[2];
	double **Ert_tying_point = mitc.mat_T_3[3];
	double **Est_tying_point = mitc.mat_T_3[4];
	double **Z1rr = mitc.mat_3_5N[0];
	double **Z2ss = mitc.mat_3_5N[1];
	double **Z1rt = mitc.mat_3_5N[2];
	double **Z2st = mitc.mat_3_5N[3];
	double **Z3rt = mitc.mat_3_5N[4];
	double **Z3st = mitc.mat_3_5N[5];
	double **Z1rs = mitc.mat_3_5N[6];
	double **Z2rs = mitc.mat_3_5N[7];
	double **ZtZ11 = mitc.mat_5N_5N[0];
	double **ZtZ22 = mitc.mat_5N_5N[1];
	double **ZtZ12 = mitc.mat_5N_5N[2];
	double **ZtZ13 = mitc.mat_5N_5N[3];
	double **ZtZ23 = mitc.mat_5N_5N[4];


	elem_dof = dim * elem_node_num;

	init_matrix(elem_dof, elem_dof, SZtZ_matrix);

	RC_TRY( set_tying_data_mitc8(local_point, &Ers_tying_num,
	                             &ErrEssErtEst_tying_num, Err_tying_point,
	                             Ess_tying_point, Ers_tying_point,
	                             Ert_tying_point, Est_tying_point, Err_h,
	                             Ess_h, Ers_h, Ert_h, Est_h) );

	for(ii1=0; ii1<Ers_tying_num; ii1++){
		double stress_rs = stress_vect[2];

		RC_TRY( set_Z1_matrix(elem_type, elem_node_num, elem_curr_vect,
		                      thickness, Ers_tying_point[ii1], Z1rs, ws) );
		RC_TRY( set_Z2_matrix(elem_type, elem_node_num, elem_curr_vect,
		                      thickness, Ers_tying_point[ii1], Z2rs, ws) );

		init_matrix(elem_dof, elem_dof, ZtZ12);
		for(ii2=0; ii2<3; ii2++){
			for(ii3=0; ii3<elem_dof; ii3++){
				for(ii4=0; ii4<elem_dof; ii4++){
					ZtZ12[ii3][ii4] += Z1rs[ii2][ii3]*Z2rs[ii2][ii4]
					                 + Z2rs[ii2][ii3]*Z1rs[ii2][ii4];
				}
			}
		}

		for(ii2=0; ii2<elem_dof; ii2++){
			for(ii3=0; ii3<elem_dof; ii3++){
				SZtZ_matrix[ii2][ii3] += Ers_h[ii1]*stress_rs*ZtZ12[ii2][ii3];
			}
		}
	}

	for(ii1=0; ii1<ErrEssErtEst_tying_num; ii1++){
		double stress_rr = stress_vect[0];
		double stress_ss = stress_vect[1];
		double stress_rt = stress_vect[3];
		double stress_st = stress_vect[4];

		RC_TRY( set_Z1_matrix(elem_type, elem_node_num, elem_curr_vect,
		                      thickness, Err_tying_point[ii1], Z1rr, ws) );
		RC_TRY( set_Z2_matrix(elem_type, elem_node_num, elem_curr_vect,
		                      thickness, Ess_tying_point[ii1], Z2ss, ws) );
		RC_TRY( set_Z1_matrix(elem_type, elem_node_num, elem_curr_vect,
		                      thickness, Ert_tying_point[ii1], Z1rt, ws) );
		RC_TRY( set_Z3_matrix(elem_type, elem_node_num, elem_curr_vect,
		                      thickness, Ert_tying_point[ii1], Z3rt, ws) );
		RC_TRY( set_Z2_matrix(elem_type, elem_node_num, elem_curr_vect,
		                      thickness, Est_tying_point[ii1], Z2st, ws) );
		RC_TRY( set_Z3_matrix(elem_type, elem_node_num, elem_curr_vect,
		                      thickness, Est_tying_point[ii1], Z3st, ws) );

		init_matrix(elem_dof, elem_dof, ZtZ11);
		init_matrix(elem_dof, elem_dof, ZtZ22);
		init_matrix(elem_dof, elem_dof, ZtZ13);
		init_matrix(elem_dof, elem_dof, ZtZ23);
		for(ii2=0; ii2<3; ii2++){
			for(ii3=0; ii3<elem_dof; ii3++){
				for(ii4=0; ii4<elem_dof; ii4++){
					ZtZ11[ii3][ii4] += Z1rr[ii2][ii3] * Z1rr[ii2][ii4];
					ZtZ22[ii3][ii4] += Z2ss[ii2][ii3] * Z2ss[ii2][ii4];
					ZtZ13[ii3][ii4] += Z1rt[ii2][ii3]*Z3rt[ii2][ii4]
					                 + Z3rt[ii2][ii3]*Z1rt[ii2][ii4];
					ZtZ23[ii3][ii4] += Z2st[ii2][ii3]*Z3st[ii2][ii4]
					                 + Z3st[ii2][ii3]*Z2st[ii2][ii4];
				}
			}
		}

		for(ii2=0; ii2<elem_dof; ii2++){
			for(ii3=0; ii3<elem_dof; ii3++){
				SZtZ_matrix[ii2][ii3] += Err_h[ii1]*stress_rr*ZtZ11[ii2][ii3];
				SZtZ_matrix[ii2][ii3] += Ess_h[ii1]*stress_ss*ZtZ22[ii2][ii3];
				SZtZ_matrix[ii2][ii3] += Ert_h[ii1]*stress_rt*ZtZ13[ii2][ii3];
				SZtZ_matrix[ii2][ii3] += Est_h[ii1]*stress_st*ZtZ23[ii2][ii3];
			}
		}
	}

	return(NORMAL_RC);
}


/* {du/dr(or du/dr1)} = [Z1]{u} */
static RC
set_Z1_matrix (ELEM_TYPE elem_type, int elem_node_num,
               DIRECTOR_VECT elem_curr_vect[],
               double thickness, const double *local_point,
               double **Z1_matrix, WORK_SET ws)
{
	int ii1;
	int dim = 5;
	double c;
	double t = local_point[2];
	double **dN = ws.mat_3_N[0];
	VECT3D V1, V2;


	RC_TRY( set_dN_matrix_shell(elem_type, local_point, dN) );

	c = 0.5 * thickness * t;

	for(ii1=0; ii1<elem_node_num; ii1++){
		double dNc = c * dN[0][ii1];

		V1 = elem_curr_vect[ii1].v1;
		V2 = elem_curr_vect[ii1].v2;

		Z1_matrix[0][ii1*dim+0] = dN[0][ii1];
		Z1_matrix[0][ii1*dim+1] = 0.0;
		Z1_matrix[0][ii1*dim+2] = 0.0;
		Z1_matrix[0][ii1*dim+3] = -dNc * V2.x;
		Z1_matrix[0][ii1*dim+4] = dNc * V1.x;

		Z1_matrix[1][ii1*dim+0] = 0.0;
		Z1_matrix[1][ii1*dim+1] = dN[0][ii1];
		Z1_matrix[1][ii1*dim+2] = 0.0;
		Z1_matrix[1][ii1*dim+3] = -dNc * V2.y;
		Z1_matrix[1][ii1*dim+4] = dNc * V1.y;

		Z1_matrix[2][ii1*dim+0] = 0.0;
		Z1_matrix[2][ii1*dim+1] = 0.0;
		Z1_matrix[2][ii1*dim+2] = dN[0][ii1];
		Z1_matrix[2][ii1*dim+3] = -dNc * V2.z;
		Z1_matrix[2][ii1*dim+4] = dNc * V1.z;
	}

	return(NORMAL_RC);
}


/* {du/ds(or du/dr2)} = [Z2]{u} */
static RC
set_Z2_matrix (ELEM_TYPE elem_type, int elem_node_num,
               DIRECTOR_VECT elem_curr_vect[],
               double thickness, const double *local_point,
               double **Z2_matrix, WORK_SET ws)
{
	int ii1;
	int dim = 5;
	double c;
	double t = local_point[2];
	double **dN = ws.mat_3_N[0];
	VECT3D V1, V2;


	RC_TRY( set_dN_matrix_shell(elem_type, local_point, dN) );

	c = 0.5 * thickness * t;

	for(ii1=0; ii1<elem_node_num; ii1++){
		double dNc = c * dN[1][ii1];

		V1 = elem_curr_vect[ii1].v1;
		V2 = elem_curr_vect[ii1].v2;

		Z2_matrix[0][ii1*dim+0] = dN[1][ii1];
		Z2_matrix[0][ii1*dim+1] = 0.0;
		Z2_matrix[0][ii1*dim+2] = 0.0;
		Z2_matrix[0][ii1*dim+3] = -dNc * V2.x;
		Z2_matrix[0][ii1*dim+4] = dNc * V1.x;

		Z2_matrix[1][ii1*dim+0] = 0.0;
		Z2_matrix[1][ii1*dim+1] = dN[1][ii1];
		Z2_matrix[1][ii1*dim+2] = 0.0;
		Z2_matrix[1][ii1*dim+3] = -dNc * V2.y;
		Z2_matrix[1][ii1*dim+4] = dNc * V1.y;

		Z2_matrix[2][ii1*dim+0] = 0.0;
		Z2_matrix[2][ii1*dim+1] = 0.0;
		Z2_matrix[2][ii1*dim+2] = dN[1][ii1];
		Z2_matrix[2][ii1*dim+3] = -dNc * V2.z;
		Z2_matrix[2][ii1*dim+4] = dNc * V1.z;
	}

	return(NORMAL_RC);
}


/* {du/dt(or du/dr3)} = [Z3]{u} */
static RC
set_Z3_matrix (ELEM_TYPE elem_type, int elem_node_num,
               DIRECTOR_VECT elem_curr_vect[],
               double thickness, const double *local_point,
               double **Z3_matrix, WORK_SET ws)
{
	int ii1;
	int dim = 5;
	double c;
	double *N = ws.vec_N[0];
	VECT3D V1, V2;


	RC_TRY( set_N_vector_shell(elem_type, local_point, N) );

	c = 0.5 * thickness;

	for(ii1=0; ii1<elem_node_num; ii1++){
		double Nc = c * N[ii1];

		V1 = elem_curr_vect[ii1].v1;
		V2 = elem_curr_vect[ii1].v2;

		Z3_matrix[0][ii1*dim+0] = 0.0;
		Z3_matrix[0][ii1*dim+1] = 0.0;
		Z3_matrix[0][ii1*dim+2] = 0.0;
		Z3_matrix[0][ii1*dim+3] = -Nc * V2.x;
		Z3_matrix[0][ii1*dim+4] = Nc * V1.x;

		Z3_matrix[1][ii1*dim+0] = 0.0;
		Z3_matrix[1][ii1*dim+1] = 0.0;
		Z3_matrix[1][ii1*dim+2] = 0.0;
		Z3_matrix[1][ii1*dim+3] = -Nc * V2.y;
		Z3_matrix[1][ii1*dim+4] = Nc * V1.y;

		Z3_matrix[2][ii1*dim+0] = 0.0;
		Z3_matrix[2][ii1*dim+1] = 0.0;
		Z3_matrix[2][ii1*dim+2] = 0.0;
		Z3_matrix[2][ii1*dim+3] = -Nc * V2.z;
		Z3_matrix[2][ii1*dim+4] = Nc * V1.z;
	}

	return(NORMAL_RC);
}


/* 第2Piola-Kirchhoff応力，Greenひずみ */
RC
PK2_stress_Green_strain_shell (int sol_ss_flag, int sol_ss_face,
                               NODE_ARRAY node, ELEMENT_ARRAY element,
                               MATERIAL_PROP_ARRAY material,
                               PHYSICAL_PROP_ARRAY physical,
                               DIRECTOR_VECT_ARRAY initial_vect,
                               DIRECTOR_VECT_ARRAY current_vect,
                               DISP_ARRAY disp,
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

			RC_TRY( local_PK2_stress_Green_strain_shell(element.array[ii1],
			                                            local_points[0], node,
			                                            material, physical,
			                                            initial_vect,
			                                            current_vect, disp,
			                                            &local_stress,
			                                            &local_strain,
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
				RC_TRY( local_PK2_stress_Green_strain_shell(element.array[ii1],
				                                            local_points[ii2],
				                                            node, material,
				                                            physical,
				                                            initial_vect,
				                                            current_vect, disp,
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
 * 局所座標系での第2Piola-Kirchhoff応力，Greeひずみn
 * ひずみ不要の場合 : strain = NULL
 * 応力不要の場合   : stress = NULL
 */
RC
local_PK2_stress_Green_strain_shell (ELEMENT element, const double *local_point,
                                     NODE_ARRAY node,
                                     MATERIAL_PROP_ARRAY material,
                                     PHYSICAL_PROP_ARRAY physical,
                                     DIRECTOR_VECT_ARRAY initial_vect,
                                     DIRECTOR_VECT_ARRAY current_vect,
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
	double strain_vect[6];
	double stress_vect[6];
	double global_strain_vect[6];
	double global_stress_vect[6];
	double **node_location = ws1.mat_N_3[0];
	double **update_location = ws1.mat_N_3[1];
	double **D_matrix = ws1.mat_6_6[0];
	double **Q_matrix = ws1.mat_6_6[1];
	DIRECTOR_VECT elem_init_vect[MAX_SHELL_NODE];
	DIRECTOR_VECT elem_curr_vect[MAX_SHELL_NODE];
	ELEM_TYPE elem_type;
	MATERIAL_PROP tmp_mat;
	VECT3D G[3], g[3];


	if( (dim = element_dim_shell(element.type)) <= 0
	||(ssnum = element_ssnum_shell(element.type)) <= 0 ){
		return(ARG_ERROR_RC);
	}

	RC_TRY( get_thickness(element, physical, &thickness) );

	elem_node_num = element.node_num;
	elem_type = element.type;

	RC_TRY( set_element_direcor_vect(element, initial_vect, elem_init_vect) );
	RC_TRY( set_element_direcor_vect(element, current_vect, elem_curr_vect) );
	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( update_node_location_matrix(&element,node,disp,update_location) );
	RC_TRY( make_G1G2G3g1g2g3_basis(elem_type, elem_node_num, elem_init_vect,
	                                thickness, node_location, local_point,
	                                G, g, ws2) );
	/* strain */
	RC_TRY( local_covariant_Green_strain_vect(elem_type, elem_node_num,
	                                          elem_init_vect, elem_curr_vect,
	                                          thickness, node_location,
	                                          update_location, local_point,
	                                          strain_vect, mitc, ws2) );
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
	mul_matrix_vect(6, 6, Q_matrix, strain_vect, global_strain_vect);

	/* stress */
	RC_TRY( make_D_matrix_shell(element, material, physical, g, D_matrix) );
	mul_matrix_vect(ssnum, ssnum, D_matrix, strain_vect, stress_vect);
	stress_vect[5] = 0.0;

	RC_TRY( set_Q_matrix4stress(G, Q_matrix) );
	mul_matrix_vect(6, 6, Q_matrix, stress_vect, global_stress_vect);

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


/* 局所座標系でのGreenひずみベクトル(共変基底成分) */
RC
local_covariant_Green_strain_vect (ELEM_TYPE elem_type, int elem_node_num,
                                   DIRECTOR_VECT elem_init_vect[],
                                   DIRECTOR_VECT elem_curr_vect[],
                                   double thickness,
                                   double **node_location,
                                   double **update_location,
                                   const double *local_point,
                                   double strain_vect[],
                                   MITC_TYING_SET mitc, WORK_SET ws)
{
	switch(elem_type){
	case ELEM_TRI1:
		RC_TRY( local_cova_Green_strain_vect_mitc3(elem_type, elem_node_num,
		                                           elem_init_vect,
		                                           elem_curr_vect,
		                                           thickness, node_location,
		                                           update_location, local_point,
		                                           strain_vect, mitc, ws) );
		break;
	case ELEM_TRI2:
		RC_TRY( local_cova_Green_strain_vect_mitc6(elem_type, elem_node_num,
		                                           elem_init_vect,
		                                           elem_curr_vect,
		                                           thickness, node_location,
		                                           update_location, local_point,
		                                           strain_vect, mitc, ws) );
		break;
	case ELEM_QUAD1:
		RC_TRY( local_cova_Green_strain_vect_mitc4(elem_type, elem_node_num,
		                                           elem_init_vect,
		                                           elem_curr_vect,
		                                           thickness, node_location,
		                                           update_location, local_point,
		                                           strain_vect, mitc, ws) );

		break;
	case ELEM_QUAD2:
		RC_TRY( local_cova_Green_strain_vect_mitc8(elem_type, elem_node_num,
		                                           elem_init_vect,
		                                           elem_curr_vect,
		                                           thickness, node_location,
		                                           update_location, local_point,
		                                           strain_vect, mitc, ws) );

		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


/* 局所座標系でのGreenひずみベクトル(共変基底成分) MITC3 */
static RC
local_cova_Green_strain_vect_mitc3 (ELEM_TYPE elem_type, int elem_node_num,
                                    DIRECTOR_VECT elem_init_vect[],
                                    DIRECTOR_VECT elem_curr_vect[],
                                    double thickness,
                                    double **node_location,
                                    double **update_location,
                                    const double *local_point,
                                    double strain_vect[],
                                    MITC_TYING_SET mitc, WORK_SET ws)
{
	int ii1;
	int ErtEst_tying_num;
	double *Ert_rt_h = mitc.vec_T[0];
	double *Ert_st_h = mitc.vec_T[1];
	double *Est_rt_h = mitc.vec_T[2];
	double *Est_st_h = mitc.vec_T[3];
	double **Ert_rt_tying_point = mitc.mat_T_3[0];
	double **Ert_st_tying_point = mitc.mat_T_3[1];
	double **Est_rt_tying_point = mitc.mat_T_3[2];
	double **Est_st_tying_point = mitc.mat_T_3[3];
	VECT3D G0[3], Gt[3];
	VECT3D G0_rt_rt[3], Gt_rt_rt[3];
	VECT3D G0_rt_st[3], Gt_rt_st[3];
	VECT3D G0_st_rt[3], Gt_st_rt[3];
	VECT3D G0_st_st[3], Gt_st_st[3];


	RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_init_vect,
	                          thickness, node_location, local_point, G0, ws) );
	RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
	                          thickness, update_location, local_point,
	                          Gt, ws) );

	strain_vect[0] = ( inner_product3d(Gt[0], Gt[0])
	                 - inner_product3d(G0[0], G0[0]) ) * 0.5;
	strain_vect[1] = ( inner_product3d(Gt[1], Gt[1])
	                 - inner_product3d(G0[1], G0[1]) ) * 0.5;
	strain_vect[2] = ( inner_product3d(Gt[0], Gt[1])
	                 - inner_product3d(G0[0], G0[1]) );

	RC_TRY( set_tying_data_mitc3(local_point, &ErtEst_tying_num,
	                             Ert_rt_tying_point, Ert_st_tying_point,
	                             Est_rt_tying_point, Est_st_tying_point,
	                             Ert_rt_h, Ert_st_h, Est_rt_h, Est_st_h) );

	strain_vect[3] = 0.0;
	strain_vect[4] = 0.0;
	for(ii1=0; ii1<ErtEst_tying_num; ii1++){
		double Ert_rt, Ert_st;
		double Est_rt, Est_st;

		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_init_vect,
		                          thickness, node_location,
		                          Ert_rt_tying_point[ii1], G0_rt_rt, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
		                          thickness, update_location,
		                          Ert_rt_tying_point[ii1], Gt_rt_rt, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_init_vect,
		                          thickness, node_location,
		                          Ert_st_tying_point[ii1], G0_rt_st, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
		                          thickness, update_location,
		                          Ert_st_tying_point[ii1], Gt_rt_st, ws) );

		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_init_vect,
		                          thickness, node_location,
		                          Est_rt_tying_point[ii1], G0_st_rt, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
		                          thickness, update_location,
		                          Est_rt_tying_point[ii1], Gt_st_rt, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_init_vect,
		                          thickness, node_location,
		                          Est_st_tying_point[ii1], G0_st_st, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
		                          thickness, update_location,
		                          Est_st_tying_point[ii1], Gt_st_st, ws) );

		Ert_rt = ( inner_product3d(Gt_rt_rt[0], Gt_rt_rt[2])
		         - inner_product3d(G0_rt_rt[0], G0_rt_rt[2]) );
		Ert_st = ( inner_product3d(Gt_rt_st[1], Gt_rt_st[2])
		         - inner_product3d(G0_rt_st[1], G0_rt_st[2]) );

		Est_rt = ( inner_product3d(Gt_st_rt[0], Gt_st_rt[2])
		         - inner_product3d(G0_st_rt[0], G0_st_rt[2]) );
		Est_st = ( inner_product3d(Gt_st_st[1], Gt_st_st[2])
		         - inner_product3d(G0_st_st[1], G0_st_st[2]) );

		strain_vect[3] += Ert_rt_h[ii1] * Ert_rt;
		strain_vect[3] += Ert_st_h[ii1] * Ert_st;

		strain_vect[4] += Est_rt_h[ii1] * Est_rt;
		strain_vect[4] += Est_st_h[ii1] * Est_st;
	}

	return(NORMAL_RC);
}


/* 局所座標系でのGreenひずみベクトル(共変基底成分) MITC6 */
static RC
local_cova_Green_strain_vect_mitc6 (ELEM_TYPE elem_type, int elem_node_num,
                                    DIRECTOR_VECT elem_init_vect[],
                                    DIRECTOR_VECT elem_curr_vect[],
                                    double thickness,
                                    double **node_location,
                                    double **update_location,
                                    const double *local_point,
                                    double strain_vect[],
                                    MITC_TYING_SET mitc, WORK_SET ws)
{
	int ii1;
	int ErrEss_tying_num;
	int Ers_tying_num;
	int ErtEst_tying_num;
	double *Err_h = mitc.vec_T[0];
	double *Ess_h = mitc.vec_T[1];
	double *Ers_rr_h = mitc.vec_T[2];
	double *Ers_ss_h = mitc.vec_T[3];
	double *Ers_rs_h = mitc.vec_T[4];
	double *Ert_rt_h = mitc.vec_T[5];
	double *Ert_st_h = mitc.vec_T[6];
	double *Est_rt_h = mitc.vec_T[7];
	double *Est_st_h = mitc.vec_T[8];
	double **Err_tying_point = mitc.mat_T_3[0];
	double **Ess_tying_point = mitc.mat_T_3[1];
	double **Ers_rr_tying_point = mitc.mat_T_3[2];
	double **Ers_ss_tying_point = mitc.mat_T_3[3];
	double **Ers_rs_tying_point = mitc.mat_T_3[4];
	double **Ert_rt_tying_point = mitc.mat_T_3[5];
	double **Ert_st_tying_point = mitc.mat_T_3[6];
	double **Est_rt_tying_point = mitc.mat_T_3[7];
	double **Est_st_tying_point = mitc.mat_T_3[8];
	VECT3D G0_rr[3], Gt_rr[3];
	VECT3D G0_ss[3], Gt_ss[3];
	VECT3D G0_rs_rr[3], Gt_rs_rr[3];
	VECT3D G0_rs_ss[3], Gt_rs_ss[3];
	VECT3D G0_rs_rs[3], Gt_rs_rs[3];
	VECT3D G0_rt_rt[3], Gt_rt_rt[3];
	VECT3D G0_rt_st[3], Gt_rt_st[3];
	VECT3D G0_st_rt[3], Gt_st_rt[3];
	VECT3D G0_st_st[3], Gt_st_st[3];


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
	strain_vect[0] = 0.0;
	strain_vect[1] = 0.0;
	for(ii1=0; ii1<ErrEss_tying_num; ii1++){
		double Err, Ess;

		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_init_vect,
		                          thickness, node_location,
		                          Err_tying_point[ii1],
		                          G0_rr, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
		                          thickness, update_location,
		                          Err_tying_point[ii1],
		                          Gt_rr, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_init_vect,
		                          thickness, node_location,
		                          Ess_tying_point[ii1],
		                          G0_ss, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
		                          thickness, update_location,
		                          Ess_tying_point[ii1],
		                          Gt_ss, ws) );

		Err = ( inner_product3d(Gt_rr[0], Gt_rr[0])
		      - inner_product3d(G0_rr[0], G0_rr[0]) ) * 0.5;
		Ess = ( inner_product3d(Gt_ss[1], Gt_ss[1])
		      - inner_product3d(G0_ss[1], G0_ss[1]) ) * 0.5;

		strain_vect[0] += Err_h[ii1] * Err;
		strain_vect[1] += Ess_h[ii1] * Ess;
	}

	/* Ers */
	strain_vect[2] = strain_vect[0] + strain_vect[1];
	for(ii1=0; ii1<Ers_tying_num; ii1++){
		double Ers_rr, Ers_ss, Ers_rs;

		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_init_vect,
		                          thickness, node_location,
		                          Ers_rr_tying_point[ii1],
		                          G0_rs_rr, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
		                          thickness, update_location,
		                          Ers_rr_tying_point[ii1],
		                          Gt_rs_rr, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_init_vect,
		                          thickness, node_location,
		                          Ers_ss_tying_point[ii1],
		                          G0_rs_ss, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
		                          thickness, update_location,
		                          Ers_ss_tying_point[ii1],
		                          Gt_rs_ss, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_init_vect,
		                          thickness, node_location,
		                          Ers_rs_tying_point[ii1],
		                          G0_rs_rs, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
		                          thickness, update_location,
		                          Ers_rs_tying_point[ii1],
		                          Gt_rs_rs, ws) );

		Ers_rr = ( inner_product3d(Gt_rs_rr[0], Gt_rs_rr[0])
		         - inner_product3d(G0_rs_rr[0], G0_rs_rr[0]) ) * 0.5;
		Ers_ss = ( inner_product3d(Gt_rs_ss[1], Gt_rs_ss[1])
		         - inner_product3d(G0_rs_ss[1], G0_rs_ss[1]) ) * 0.5;
		Ers_rs = ( inner_product3d(Gt_rs_rs[0], Gt_rs_rs[1])
	             - inner_product3d(G0_rs_rs[0], G0_rs_rs[1]) );

		strain_vect[2] += Ers_rr_h[ii1] * Ers_rr;
		strain_vect[2] += Ers_ss_h[ii1] * Ers_ss;
		strain_vect[2] += Ers_rs_h[ii1] * Ers_rs;
	}

	/* Ert, Est */
	strain_vect[3] = 0.0;
	strain_vect[4] = 0.0;
	for(ii1=0; ii1<ErtEst_tying_num; ii1++){
		double Ert_rt, Ert_st;
		double Est_rt, Est_st;

		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_init_vect,
		                          thickness, node_location,
		                          Ert_rt_tying_point[ii1], G0_rt_rt, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
		                          thickness, update_location,
		                          Ert_rt_tying_point[ii1], Gt_rt_rt, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_init_vect,
		                          thickness, node_location,
		                          Ert_st_tying_point[ii1], G0_rt_st, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
		                          thickness, update_location,
		                          Ert_st_tying_point[ii1], Gt_rt_st, ws) );

		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_init_vect,
		                          thickness, node_location,
		                          Est_rt_tying_point[ii1], G0_st_rt, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
		                          thickness, update_location,
		                          Est_rt_tying_point[ii1], Gt_st_rt, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_init_vect,
		                          thickness, node_location,
		                          Est_st_tying_point[ii1], G0_st_st, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
		                          thickness, update_location,
		                          Est_st_tying_point[ii1], Gt_st_st, ws) );

		Ert_rt = ( inner_product3d(Gt_rt_rt[0], Gt_rt_rt[2])
		         - inner_product3d(G0_rt_rt[0], G0_rt_rt[2]) );
		Ert_st = ( inner_product3d(Gt_rt_st[1], Gt_rt_st[2])
		         - inner_product3d(G0_rt_st[1], G0_rt_st[2]) );

		Est_rt = ( inner_product3d(Gt_st_rt[0], Gt_st_rt[2])
		         - inner_product3d(G0_st_rt[0], G0_st_rt[2]) );
		Est_st = ( inner_product3d(Gt_st_st[1], Gt_st_st[2])
		         - inner_product3d(G0_st_st[1], G0_st_st[2]) );

		strain_vect[3] += Ert_rt_h[ii1] * Ert_rt;
		strain_vect[3] += Ert_st_h[ii1] * Ert_st;

		strain_vect[4] += Est_rt_h[ii1] * Est_rt;
		strain_vect[4] += Est_st_h[ii1] * Est_st;
	}

	return(NORMAL_RC);
}


/* 局所座標系でのGreenひずみベクトル(共変基底成分) MITC4 */
static RC
local_cova_Green_strain_vect_mitc4 (ELEM_TYPE elem_type, int elem_node_num,
                                    DIRECTOR_VECT elem_init_vect[],
                                    DIRECTOR_VECT elem_curr_vect[],
                                    double thickness,
                                    double **node_location,
                                    double **update_location,
                                    const double *local_point,
                                    double strain_vect[],
                                    MITC_TYING_SET mitc, WORK_SET ws)
{
	int ii1;
	int ErtEst_tying_num;
	double *Ert_h = mitc.vec_T[0];
	double *Est_h = mitc.vec_T[1];
	double **Ert_tying_point = mitc.mat_T_3[0];
	double **Est_tying_point = mitc.mat_T_3[1];
	VECT3D G0[3], Gt[3];
	VECT3D G0_rt[3], Gt_rt[3];
	VECT3D G0_st[3], Gt_st[3];


	RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_init_vect,
	                          thickness, node_location, local_point, G0, ws) );
	RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
	                          thickness, update_location, local_point,
	                          Gt, ws) );

	strain_vect[0] = ( inner_product3d(Gt[0], Gt[0])
	                 - inner_product3d(G0[0], G0[0]) ) * 0.5;
	strain_vect[1] = ( inner_product3d(Gt[1], Gt[1])
	                 - inner_product3d(G0[1], G0[1]) ) * 0.5;
	strain_vect[2] = ( inner_product3d(Gt[0], Gt[1])
	                 - inner_product3d(G0[0], G0[1]) );

	RC_TRY( set_tying_data_mitc4(local_point, &ErtEst_tying_num,
	                             Ert_tying_point, Est_tying_point,
	                             Ert_h, Est_h) );

	strain_vect[3] = 0.0;
	strain_vect[4] = 0.0;
	for(ii1=0; ii1<ErtEst_tying_num; ii1++){
		double Ert, Est;

		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_init_vect,
		                          thickness, node_location,
		                          Ert_tying_point[ii1], G0_rt, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
		                          thickness, update_location,
		                          Ert_tying_point[ii1], Gt_rt, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_init_vect,
		                          thickness, node_location,
		                          Est_tying_point[ii1], G0_st, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
		                          thickness, update_location,
		                          Est_tying_point[ii1], Gt_st, ws) );

		Ert = ( inner_product3d(Gt_rt[0], Gt_rt[2])
		      - inner_product3d(G0_rt[0], G0_rt[2]) );
		Est = ( inner_product3d(Gt_st[1], Gt_st[2])
		      - inner_product3d(G0_st[1], G0_st[2]) );

		strain_vect[3] += Ert_h[ii1] * Ert;
		strain_vect[4] += Est_h[ii1] * Est;
	}

	return(NORMAL_RC);
}


/* 局所座標系でのGreenひずみベクトル(共変基底成分) MITC8 */
static RC
local_cova_Green_strain_vect_mitc8 (ELEM_TYPE elem_type, int elem_node_num,
                                    DIRECTOR_VECT elem_init_vect[],
                                    DIRECTOR_VECT elem_curr_vect[],
                                    double thickness,
                                    double **node_location,
                                    double **update_location,
                                    const double *local_point,
                                    double strain_vect[],
                                    MITC_TYING_SET mitc, WORK_SET ws)
{
	int ii1;
	int Ers_tying_num;
	int ErrEssErtEst_tying_num;
	double *Err_h = mitc.vec_T[0];
	double *Ess_h = mitc.vec_T[1];
	double *Ers_h = mitc.vec_T[2];
	double *Ert_h = mitc.vec_T[3];
	double *Est_h = mitc.vec_T[4];
	double **Err_tying_point = mitc.mat_T_3[0];
	double **Ess_tying_point = mitc.mat_T_3[1];
	double **Ers_tying_point = mitc.mat_T_3[2];
	double **Ert_tying_point = mitc.mat_T_3[3];
	double **Est_tying_point = mitc.mat_T_3[4];
	VECT3D G0_rr[3], Gt_rr[3];
	VECT3D G0_ss[3], Gt_ss[3];
	VECT3D G0_rs[3], Gt_rs[3];
	VECT3D G0_rt[3], Gt_rt[3];
	VECT3D G0_st[3], Gt_st[3];


	RC_TRY( set_tying_data_mitc8(local_point, &Ers_tying_num,
	                             &ErrEssErtEst_tying_num, Err_tying_point,
	                             Ess_tying_point, Ers_tying_point,
	                             Ert_tying_point, Est_tying_point, Err_h,
	                             Ess_h, Ers_h, Ert_h, Est_h) );

	init_vect(5, strain_vect);

	for(ii1=0; ii1<Ers_tying_num; ii1++){
		double Ers;

		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_init_vect,
		                          thickness, node_location,
		                          Ers_tying_point[ii1], G0_rs, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
		                          thickness, update_location,
		                          Ers_tying_point[ii1], Gt_rs, ws) );

		Ers = ( inner_product3d(Gt_rs[0], Gt_rs[1])
		      - inner_product3d(G0_rs[0], G0_rs[1]) );

		strain_vect[2] += Ers_h[ii1] * Ers;
	}

	for(ii1=0; ii1<ErrEssErtEst_tying_num; ii1++){
		double Err;
		double Ess;
		double Ert;
		double Est;

		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_init_vect,
		                          thickness, node_location,
		                          Err_tying_point[ii1], G0_rr, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
		                          thickness, update_location,
		                          Err_tying_point[ii1], Gt_rr, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_init_vect,
		                          thickness, node_location,
		                          Ess_tying_point[ii1], G0_ss, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
		                          thickness, update_location,
		                          Ess_tying_point[ii1], Gt_ss, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_init_vect,
		                          thickness, node_location,
		                          Ert_tying_point[ii1], G0_rt, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
		                          thickness, update_location,
		                          Ert_tying_point[ii1], Gt_rt, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_init_vect,
		                          thickness, node_location,
		                          Est_tying_point[ii1], G0_st, ws) );
		RC_TRY( make_G1G2G3_basis(elem_type, elem_node_num, elem_curr_vect,
		                          thickness, update_location,
		                          Est_tying_point[ii1], Gt_st, ws) );

		Err = ( inner_product3d(Gt_rr[0], Gt_rr[0])
		      - inner_product3d(G0_rr[0], G0_rr[0]) ) * 0.5;
		Ess = ( inner_product3d(Gt_ss[1], Gt_ss[1])
		      - inner_product3d(G0_ss[1], G0_ss[1]) ) * 0.5;
		Ert = ( inner_product3d(Gt_rt[0], Gt_rt[2])
		      - inner_product3d(G0_rt[0], G0_rt[2]) );
		Est = ( inner_product3d(Gt_st[1], Gt_st[2])
		      - inner_product3d(G0_st[1], G0_st[2]) );

		strain_vect[0] += Err_h[ii1] * Err;
		strain_vect[1] += Ess_h[ii1] * Ess;
		strain_vect[3] += Ert_h[ii1] * Ert;
		strain_vect[4] += Est_h[ii1] * Est;
	}

	return(NORMAL_RC);
}


/* 虚偽のドリリング剛性の設定 */
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


/* 要素接線剛性マトリクスの回転自由度を全体座標系成分へ変換 */
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


/* 反変基底での第2Piola-Kirchhoff応力を全体座標系成分へ変換するマトリクス */
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


/* 共変基底でのGreenひずみを全体座標系成分へ変換するマトリクス */
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


