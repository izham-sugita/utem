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

/* $Id: nl_solver.c 553 2005-11-09 08:29:48Z aoyama $ */

#include <stdio.h>
#include <math.h>
#include "fem_solver.h"
#include "fem_struct.h"
#include "math_utl.h"
#include "sky_cholesky.h"
#include "rc.h"
#include "log_printf.h"
#include "nonzero_cg.h"
#include "memory_manager.h"

#define NR_MAX         (50)           /* Newton-Raphson法最大反復回数 */
/*
 * Configuration of arc length method
 */
#define SCALING_FACTOR (1.0e+0)       /* 弧長増分法スケーリング係数 */
#define AL_MAX         (10)           /* 最大反復回数 */
#define MAX_BISECT     (10)           /* 2分化処理最大回数 */ 
#define DESIRED_NUM    (3)            /* 設定反復回数 */
#define MIN_RATIO      (0.33)         /* 下界弧長比 */
#define MAX_RATIO      (3.0)          /* 上界弧長比 */


static RC make_AG_matrix(ELEMENT element, const double *g_point, WORK_SET ws,
                         double **node_location, DISP_ARRAY disp,
                         double **A_matrix, double **G_matrix);
static RC set_G_matrix(int dim, int node_num,
                       double **dNxyz, double **G_matrix);
static RC set_A_matrix(int dim, double **duvw, double **A_matrix);
static RC set_subM_matrix(int dim, double stress, double **M_matrix, int x,
                          int y);
static RC set_rest_zero(NODE_ARRAY node, BC_ARRAY rest, double resid_v[]);
static RC update_node_location_matrix(ELEMENT *element, NODE_ARRAY node,
                                      DISP_ARRAY disp, double **matrix);


/* 
 * Reference book is
 * Non-linear Finite Element Analysis of Solids and Structures
 * Volume 1:Essentials
 * M.A.Crisfield 
 * 超円筒型弧長増分法による非線形座屈後解析ソルバー(ベータ版)
 * 連立方程式の求解にはスカイライン使用
 * 指定したload_factorを越えた時点で解析終了
 */
RC
fem_cylindrical_AL_solve (NODE_ARRAY node, ELEMENT_ARRAY element,
                          MATERIAL_PROP_ARRAY material,
                          PHYSICAL_PROP_ARRAY physical,
                          BC_ARRAY rest, BC_ARRAY force,
                          DEFAULT_BC_ARRAY b_force,
                          BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                          double load_factor,
                          DISP_ARRAY *disp, STRESS_ARRAY *stress,
                          STRAIN_ARRAY *strain, int sol_ss_flag)
{
	SKYLINE_MATRIX tan_K;
	DISP_ARRAY old_disp;
	int ii1, ii2;
	int dim, total_dof;
	int bisec_count;
	int converge_flag;
	double beta1 = 1.0e-3;
	double beta2 = 1.0e-7;
	double delta_AL;
	double lambda, old_lambda;
	double abs_lambda, old_abs_lambda;
	double delta_lambda;
	double *f_vect;
	double *resid_v;
	double *d_vect;
	double *old_d_vect;
	double *sub_vect;
	double **delta_d;

	dim = analysis_dim(element);
	if(dim <= 1 || dim > 3) return(ARG_ERROR_RC);

	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

	RC_TRY( fem_total_dof(node, element, &total_dof) );

	RC_TRY( allocate1D(total_dof, &f_vect) );
	RC_TRY( allocate1D(total_dof, &resid_v) );
	RC_TRY( allocate1D(total_dof, &d_vect) );
	RC_TRY( allocate1D(total_dof, &old_d_vect) );
	RC_TRY( allocate1D(total_dof, &sub_vect) );
	RC_TRY( allocate2D(2, total_dof, &delta_d) );
	RC_TRY( vect2disp(dim, node, d_vect, disp) );
	RC_TRY( allocate_disp_array(node.size, &old_disp) );

	/* 外力入力 */
	RC_TRY( add_nodal_force(dim, node, force, f_vect) );
	RC_TRY( add_body_force(dim, node, element, material, physical, b_force,
	                       f_vect) );
	/* 熱ひずみ(初期ひずみ) */
	RC_TRY( add_temp_force(dim, node, element, material, physical, temp,
	                       def_temp, f_vect) );

	bisec_count = 0;
	lambda = 0.0;
	abs_lambda = 0.0;
	delta_AL = load_factor * 0.1;
	while(lambda < load_factor){
		int cut_AL_flag = 0;
		double arc_ratio;
		double old_error1 = 1.0;
		double old_error2 = 1.0;

		RC_TRY(log_printf(4, "_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/\n") );
		RC_TRY(log_printf(4, "Arc length     : %10.5e\n", delta_AL) );
		RC_TRY(log_printf(4, "Current factor : %10.5e\n", lambda) );
		RC_TRY(log_printf(4, "Load factor    : %10.5e\n", load_factor) );

		/* Predictor solution */
		RC_TRY( make_NL_g_matrix(node, element, material, physical,
		                         *disp, &tan_K) );
		for(ii1=0; ii1<total_dof; ii1++){
			delta_d[0][ii1] = f_vect[ii1];
		}
		RC_TRY( modify_g_matrix(node, rest, tan_K, delta_d[0]) );
		RC_TRY( decomp_g_matrix(tan_K) );
		RC_TRY( log_printf(4, "\n") );
		RC_TRY( s_cholesky_subst(tan_K.dof, tan_K.index1, tan_K.index2,
		                         tan_K.array, delta_d[0],0) );
		RC_TRY( free_g_matrix(&tan_K) );

		delta_lambda = delta_AL
		             / sqrt( inner_product(total_dof, delta_d[0], delta_d[0]) );

		/* 求解方向判定 */
		if( inner_product(total_dof, old_d_vect, delta_d[0]) < 0.0 ){
			delta_lambda *= -1.0;
		}

		for(ii1=0; ii1<total_dof; ii1++){
			d_vect[ii1] = delta_lambda * delta_d[0][ii1];
		}

		/* 前回の結果を記憶(複素解or反復回数超過時のリスタート処理のため) */
		old_lambda = lambda;
		RC_TRY( copy_disp_array(*disp, &old_disp) );

		/* 荷重係数，全体変位更新(方向解) */
		lambda += delta_lambda;
		RC_TRY( add_vect2disp(dim, node, d_vect, *disp, 1.0) );

		old_abs_lambda = abs_lambda;
		abs_lambda += fabs(delta_lambda);

		/* 内力，残差力計算(方向解) */
		RC_TRY( make_NL_react_force(element, node, material, physical,
		                            *disp, resid_v) );
		for(ii2=0; ii2<total_dof; ii2++){
			resid_v[ii2] = lambda*f_vect[ii2] - resid_v[ii2];
			delta_d[0][ii2] = f_vect[ii2];
			delta_d[1][ii2] = resid_v[ii2];
		}
		RC_TRY( set_rest_zero(node, rest, resid_v) );

		/* AL反復(Newton-Raphson法) */
		RC_TRY( log_printf(4, " +---------------------+\n") );
		RC_TRY( log_printf(4, " | N-R iteration START |\n") );
		RC_TRY( log_printf(4, " +---------------------+\n") );
		converge_flag = 0;
		for(ii1=0; ii1<AL_MAX; ii1++){
			double a1, a2, a3, a4, a5;
			double root_value;
			double lambda1, lambda2;
			double du_norm, du_r, u_norm, u_f;

			RC_TRY( make_NL_g_matrix(node, element, material, physical,
			                         *disp, &tan_K) );
			RC_TRY( modify_g_matrix_n(node, rest, tan_K, delta_d, 2) );
			RC_TRY( decomp_g_matrix(tan_K) );
			RC_TRY( log_printf(4, "\n") );
			for(ii2=0; ii2<2; ii2++){
				RC_TRY(s_cholesky_subst(tan_K.dof, tan_K.index1, tan_K.index2,
				                        tan_K.array, delta_d[ii2], 0) );
			}
			RC_TRY( free_g_matrix(&tan_K) );

			/* 2次方程式の係数(a1,a2,a3)，鋭角検査の係数(a4,a5) */
			for(ii2=0; ii2<total_dof; ii2++){
				sub_vect[ii2] = d_vect[ii2] + delta_d[1][ii2];
			}
			a1 = inner_product(total_dof, delta_d[0], delta_d[0]);
			a2 = 2.0 * inner_product(total_dof, delta_d[0], sub_vect);
			a3 = inner_product(total_dof, sub_vect, sub_vect)
			   - (delta_AL * delta_AL);
			a4 = inner_product(total_dof, d_vect, delta_d[1])
			   + inner_product(total_dof, d_vect, d_vect);
			a5 = inner_product(total_dof, d_vect, delta_d[0]);

			/* 線形解 */
			if( nearly_eq(a1, 0.0) ){
				delta_lambda = - a3 / a2;
			/* 2次方程式求解 */
			}else{
				root_value = a2*a2 - 4.0*a1*a3;
				/* 複素解の発生 */
				if(root_value < 0.0){
					cut_AL_flag = 1;
					log_printf(4, "Encounter the Complex number!!\n");
					break;
				}

				lambda1 = (-a2 + sqrt(root_value)) / (2.0*a1);
				lambda2 = (-a2 - sqrt(root_value)) / (2.0*a1);

				/* 鋭角判定 */
				if( (a4+a5*lambda1) < (a4+a5*lambda2) ){
					delta_lambda = lambda2;
				}else{
					delta_lambda = lambda1;
				}
			}

			for(ii2=0; ii2<total_dof; ii2++){
				sub_vect[ii2] = delta_d[1][ii2]
				              + delta_lambda*delta_d[0][ii2];
				d_vect[ii2] += sub_vect[ii2];
			}

			/* 荷重係数，全体変位更新(AL反復) */
			lambda += delta_lambda;
			RC_TRY( add_vect2disp(dim, node, sub_vect, *disp, 1.0) );

			abs_lambda += fabs(delta_lambda);

			/* 内力，残差力計算(AL反復) */
			RC_TRY( make_NL_react_force(element, node, material, physical,
			                            *disp, resid_v) );
			for(ii2=0; ii2<total_dof; ii2++){
				resid_v[ii2] = lambda*f_vect[ii2] - resid_v[ii2];
				delta_d[0][ii2] = f_vect[ii2];
				delta_d[1][ii2] = resid_v[ii2];
			}
			RC_TRY( set_rest_zero(node, rest, resid_v) );

			/* 変位，エネルギー収束判定 */
			du_norm = sqrt( inner_product(total_dof, sub_vect, sub_vect) );
			du_r = 0.0;
			for(ii2=0; ii2<total_dof; ii2++){
				du_r += fabs(sub_vect[ii2] * resid_v[ii2]);
			}
			RC_TRY( copy_disp2vect(dim, node, *disp, sub_vect, 1.0) );
			u_norm = sqrt( inner_product(total_dof, sub_vect, sub_vect) );
			u_f = 0.0;
			for(ii2=0; ii2<total_dof; ii2++){
				u_f += fabs(sub_vect[ii2] * abs_lambda * f_vect[ii2]);
			}
			/* 発散判定 */
			if( (du_norm/u_norm) > old_error1
			||  (du_r/u_f) > old_error2 ){
				break;
			}
			RC_TRY( log_printf(4, "    Disp error  :%15.7e\n",du_norm/u_norm) );
			RC_TRY( log_printf(4, "    Energy error:%15.7e\n", du_r/u_f) );
			old_error1 = du_norm / u_norm;
			old_error2 = du_r / u_f;
			/* 収束判定 */
			if(du_norm < beta1*u_norm + ABS_TOL
			|| du_r < beta2*u_f + ABS_TOL){
				converge_flag = 1;
				break;
			}
		}
		RC_TRY( log_printf(4, "\n") );


		/* 複素解発生による弧長の2分化処理．その後リスタート */
		if(cut_AL_flag == 1 || converge_flag == 0){
			/* 2分化回数判定，指定回数オーバー時は計算終了 */
			if(bisec_count > MAX_BISECT){
				log_printf(3, "Bisection count OVER!!\n");
				return(NORMAL_RC);
			}
			delta_AL *= 0.5;
			lambda = old_lambda;
			abs_lambda = old_abs_lambda;
			RC_TRY( copy_disp_array(old_disp, disp) );
			bisec_count++;
			continue;
		}

		/* 反復回数による弧長設定 */
		arc_ratio = sqrt((double)DESIRED_NUM / (double)(ii1+1));
		/* 弧長比(arc_ratio)の上・下界判定 */
		if(arc_ratio < MIN_RATIO){
			arc_ratio = MIN_RATIO;
		}else if(arc_ratio > MAX_RATIO){
			arc_ratio = MAX_RATIO;
		}
		delta_AL *= arc_ratio;

		if( nearly_eq(delta_AL, 0.0) ) break;

		/* 2分化回数初期化 */
		bisec_count = 0;

		for(ii1=0; ii1<total_dof; ii1++){
			old_d_vect[ii1] = d_vect[ii1];
		}
	}
	RC_TRY(log_printf(4, "_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/\n") );
	RC_TRY(log_printf(4, "Arc length     : %10.5e\n", delta_AL) );
	RC_TRY(log_printf(4, "Current factor : %10.5e\n", lambda) );
	RC_TRY(log_printf(4, "Load factor    : %10.5e\n", load_factor) );
	/* メモリ開放 */
	RC_TRY( free1D(total_dof, &f_vect) );
	RC_TRY( free1D(total_dof, &resid_v) );
	RC_TRY( free1D(total_dof, &d_vect) );
	RC_TRY( free1D(total_dof, &old_d_vect) );
	RC_TRY( free1D(total_dof, &sub_vect) );
	RC_TRY( free2D(2, total_dof, &delta_d) );
	RC_TRY( free_disp_array(&old_disp) );

	/* 第2Piola-Kirchhoff応力，Green歪み */
	RC_TRY( PK2_stress_Green_strain(sol_ss_flag, node, element, *disp,
	                                material, physical, stress, strain) );

	return(NORMAL_RC);
}


/*
 * 微小歪み大変形解析 Total Lagrange法ソルバー(3次元問題専用)
 * 連立方程式の求解にはICCG使用，全体剛性マトリクスはノンゼロタイプ(3次元専用)
 * 応力は第2Piola-Kirchhoff応力，歪みはGreen歪みを計算
 * ただし，座屈後解析不可
 */
RC
fem_NL_solve_iccg3 (NODE_ARRAY node, ELEMENT_ARRAY element,
                    MATERIAL_PROP_ARRAY material,
                    PHYSICAL_PROP_ARRAY physical,
                    BC_ARRAY rest, BC_ARRAY force,
                    DEFAULT_BC_ARRAY b_force,
                    BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                    DISP_ARRAY *disp, STRESS_ARRAY *stress,
                    STRAIN_ARRAY *strain, int sol_ss_flag)
{
	NONZERO_MATRIX3 tan_K;
	int ii1, ii2;
	int dim, total_dof;
	double beta1 = 1.0e-3;
	double beta2 = 1.0e-7;
	double *f_vect;
	double *resid_v;
	double *react_v;
	double *delta_d;


	dim = analysis_dim(element);
	if(dim != 3) return(ARG_ERROR_RC);

	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

	/* 全自由度 */
	RC_TRY( fem_total_dof(node, element, &total_dof) );

	RC_TRY( allocate1D(total_dof, &f_vect) );
	RC_TRY( allocate1D(total_dof, &resid_v) );
	RC_TRY( allocate1D(total_dof, &react_v) );
	RC_TRY( allocate1D(total_dof, &delta_d) );
	RC_TRY( vect2disp(dim, node, delta_d, disp) );

	/* 節点荷重，重力入力 */
	RC_TRY( add_nodal_force(dim, node, force, f_vect) );
	RC_TRY( add_body_force(dim, node, element, material, physical, b_force,
	                       f_vect) );
	/* 熱ひずみ(初期ひずみ)による荷重入力 */
	RC_TRY( add_temp_force(dim, node, element, material, physical, temp,
	                       def_temp, f_vect) );

	for(ii1=0; ii1<total_dof; ii1++){
		resid_v[ii1] = f_vect[ii1];
	}

	for(ii1=0; ii1<NR_MAX; ii1++){
		double du_norm, du_r, u_norm, u_f;

		RC_TRY( log_printf(5, " +-------------------------+\n") );
		RC_TRY( log_printf(5, " | N-R iteration NUM : %3d |\n", ii1+1) );
		RC_TRY( log_printf(5, " +-------------------------+\n") );

		if(ii1 == 0){
			/* 初期接線剛性マトリクス(通常の剛性マトリクス)作成 */
			RC_TRY( allocate_nonzero3(node, element, &tan_K) );
			for(ii2=0; ii2<element.size; ii2++){
				ELEM_MATRIX elem_matrix;
				if(element.array[ii2].label < 0) continue;

				RC_TRY( make_elem_matrix(element.array[ii2], node, material,
				                         physical, &elem_matrix) );
				RC_TRY( impose_nonzero3(tan_K, elem_matrix, 1.0) );
				RC_TRY( free_elem_matrix(&elem_matrix) );
			}
		}else{
			/* 接線剛性マトリクス作成 */
			RC_TRY( make_NL_nonzero_matrix3(node, element, material,
			                                physical, *disp, &tan_K) );
		}

		RC_TRY( modify_nonzero3(node, rest, tan_K, resid_v) );
		RC_TRY( iccg_nonzero3s(tan_K, resid_v, delta_d, 1, 1) );
		RC_TRY( free_nonzero_matrix3(&tan_K) );
		RC_TRY( add_vect2disp(dim, node, delta_d, *disp, 1.0) );

		/* 反力計算 */
		RC_TRY( make_NL_react_force(element, node, material, physical,
		                            *disp, react_v) );
		/* 残差荷重ベクトル */
		for(ii2=0; ii2<total_dof; ii2++){
			resid_v[ii2] = f_vect[ii2] - react_v[ii2];
		}
		RC_TRY( set_rest_zero(node, rest, resid_v) );

		/* 変位，エネルギー収束判定 */
		du_norm = sqrt( inner_product(total_dof, delta_d, delta_d) );
		du_r = 0.0;
		for(ii2=0; ii2<total_dof; ii2++){
			du_r += fabs(delta_d[ii2] * resid_v[ii2]);
		}
		RC_TRY( copy_disp2vect(dim, node, *disp, delta_d, 1.0) );
		u_norm = sqrt( inner_product(total_dof, delta_d, delta_d) );
		u_f = 0.0;
		for(ii2=0; ii2<total_dof; ii2++){
			u_f += fabs(delta_d[ii2] * f_vect[ii2]);
		}
		RC_TRY( log_printf(5, "    Disp   error :%15.7e\n", du_norm/u_norm) );
		RC_TRY( log_printf(5, "    Energy error :%15.7e\n", du_r/u_f) );
		/* 収束判定 */
		if(du_norm < beta1*u_norm + ABS_TOL
		|| du_r < beta2*u_f + ABS_TOL){
			break;
		}
		RC_TRY( log_printf(5, "\n") );
	}
	if(ii1 == NR_MAX){
		log_printf(3, "N-R iteration NOT converged!!\n");
	}
	/* メモリ開放 */
	RC_TRY( mm_free(f_vect) );
	RC_TRY( mm_free(resid_v) );
	RC_TRY( mm_free(react_v) );
	RC_TRY( mm_free(delta_d) );

	/* 第2Piola-Kirchhoff応力，Green歪み */
	RC_TRY( PK2_stress_Green_strain(sol_ss_flag, node, element, *disp,
	                                material, physical, stress, strain) );

	return(NORMAL_RC);
}


/*
 * 全体接線剛性マトリクス作成(Total Lagrange法)
 * マトリクスはスカイラインタイプを使用
 */
RC
make_NL_g_matrix (NODE_ARRAY node, ELEMENT_ARRAY element, 
                  MATERIAL_PROP_ARRAY material,
                  PHYSICAL_PROP_ARRAY physical,
                  DISP_ARRAY disp,
                  SKYLINE_MATRIX *K)
{
	int ii1;
	int dim;
	WORK_SET ws;

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	dim = analysis_dim(element);
	if(dim <= 1 || dim > 3) return(ARG_ERROR_RC);

	RC_TRY( allocate_g_matrix(dim, node, element, K) );
	RC_TRY( allocate_work_set(&ws) );

	for(ii1=0; ii1<element.size; ii1++){
		ELEM_MATRIX elem_matrix;

		if(element.array[ii1].label < 0) continue;

		RC_TRY( make_tangent_elem_matrix(element.array[ii1], node,
		                                 material, physical, ws, disp,
		                                 &elem_matrix) );
		RC_TRY( impose_elem_matrix(*K, elem_matrix, 1.0) );
		RC_TRY( free_elem_matrix(&elem_matrix) );
	}

	RC_TRY( free_work_set(&ws) );

	return(NORMAL_RC);
}


/*
 * 全体接線剛性マトリクス作成(Total Lagrange法)
 * マトリクスはノンゼロタイプ(3次元問題専用)を使用
*/
RC
make_NL_nonzero_matrix3 (NODE_ARRAY node, ELEMENT_ARRAY element, 
                         MATERIAL_PROP_ARRAY material,
                         PHYSICAL_PROP_ARRAY physical,
                         DISP_ARRAY disp,
                         NONZERO_MATRIX3 *K)
{
	int ii1;
	int dim;
	WORK_SET ws;

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	dim = analysis_dim(element);
	if(dim != 3) return(ARG_ERROR_RC);

	RC_TRY( allocate_nonzero3(node, element, K) );
	RC_TRY( allocate_work_set(&ws) );

	for(ii1=0; ii1<element.size; ii1++){
		ELEM_MATRIX elem_matrix;

		if(element.array[ii1].label < 0) continue;

		RC_TRY( make_tangent_elem_matrix(element.array[ii1], node,
		                                 material, physical, ws, disp,
		                                 &elem_matrix) );
		RC_TRY( impose_nonzero3(*K, elem_matrix, 1.0) );
		RC_TRY( free_elem_matrix(&elem_matrix) );
	}

	RC_TRY( free_work_set(&ws) );

	return(NORMAL_RC);
}


/*
 * Reference book is
 * THE FINITE ELEMENT METHOD : Solid mechanics
 * 5th edition
 * O.C.Zienkiewicz, R.L.Taylor
 * 要素接線剛性マトリクス作成(Total Lagrange法)
 */
RC
make_tangent_elem_matrix (ELEMENT element, NODE_ARRAY node,
                          MATERIAL_PROP_ARRAY material,
                          PHYSICAL_PROP_ARRAY physical,
                          WORK_SET ws, DISP_ARRAY disp,
                          ELEM_MATRIX *elem_matrix)
{
	int dim, ssnum;
	int ii1, ii2, ii3;
	int g_num_point;
	double thickness;
	double det_J;
	double PK2_stress_v[6];
	double Green_strain_v[6];
	double disp_v[3*MAX_NODE];
	double *g_weights = ws.vec_G[0];
	double **node_location = ws.mat_N_3[0];
	double **D_matrix = ws.mat_6_6[0];
	double **B_matrix = ws.mat_6_3N[0];
	double **A_matrix = ws.mat_9_9[0];
	double **G_matrix = ws.mat_9_3N[0];
	double **AG = ws.mat_6_3N[1];
	double **B_strain = ws.mat_6_3N[2];
	double **M_matrix = ws.mat_9_9[1];
	double **BtDB = ws.mat_3N_3N[0];
	double **GtMG = ws.mat_3N_3N[1];
	double **g_points = ws.mat_G_4[0];

	if(element.label < 0) return(ARG_ERROR_RC);

	if( (elem_matrix->dim = element_dim(element.type)) <= 0 ||
	    (ssnum = element_ssnum(element.type)) <= 0 ){
		return(ARG_ERROR_RC);
	}

	dim = elem_matrix->dim;
	elem_matrix->element = element.label;
	elem_matrix->size = (element.node_num) * dim;

	/* elem_matrix->array[][]の確保 */
	RC_TRY( allocate2D(elem_matrix->size, elem_matrix->size,
	                                    &(elem_matrix->matrix)) );
	RC_TRY( allocate1I(elem_matrix->size, &(elem_matrix->index)) );

	for(ii1=0; ii1<element.node_num; ii1++){
		int index = node_renum_index(node, element.node[ii1]);
		if(index < 0) return(SEEK_ERROR_RC);
		for(ii2=0; ii2<dim; ii2++){
			elem_matrix->index[ii1*dim+ii2] = index*dim + ii2;
		}
	}

	RC_TRY( Gauss_const_default(element.type, g_points, g_weights,
	                                                   &g_num_point) );
	RC_TRY( fill_disp_vector(element, disp, disp_v) );
	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( get_D_matrix(element, material, physical, D_matrix) );
	RC_TRY( get_thickness(element, physical, &thickness) );

	for(ii1=0; ii1<g_num_point; ii1++){
		STRESS_STRAIN Green_strain;

		RC_TRY( make_B_matrix(element, g_points[ii1], node_location, B_matrix,
		                      &det_J) );
		RC_TRY( make_AG_matrix(element, g_points[ii1], ws, node_location, disp, 
		                       A_matrix, G_matrix) );
		mul_matrix(ssnum,dim*dim,dim*element.node_num,A_matrix,G_matrix,AG);
		RC_TRY( add_matrix(ssnum, dim*element.node_num, 1.0, B_matrix, 0.5, AG,
		                   B_strain) );
		mul_matrix_vect(ssnum, dim*element.node_num, B_strain, disp_v,
		                Green_strain_v);
		RC_TRY( vect2stress_strain(ssnum, Green_strain_v, &Green_strain) );
		RC_TRY( add_matrix(ssnum, dim*element.node_num, 1.0, B_matrix, 1.0, AG,
		                   B_matrix) );

		if(dim == 2){
			for(ii2=0; ii2<ssnum; ii2++){
				PK2_stress_v[ii2] = D_matrix[ii2][0]*Green_strain.v.x
				                  + D_matrix[ii2][1]*Green_strain.v.y
				                  + D_matrix[ii2][2]*Green_strain.v.xy;
			}
		}else{ /* dim == 3 */
			for(ii2=0; ii2<ssnum; ii2++){
				PK2_stress_v[ii2] = D_matrix[ii2][0]*Green_strain.v.x
				                  + D_matrix[ii2][1]*Green_strain.v.y
				                  + D_matrix[ii2][2]*Green_strain.v.z
				                  + D_matrix[ii2][3]*Green_strain.v.yz
				                  + D_matrix[ii2][4]*Green_strain.v.zx
				                  + D_matrix[ii2][5]*Green_strain.v.xy;
			}
		}

		RC_TRY( set_M_matrix(dim, PK2_stress_v, M_matrix) );

		mul_matrix_AtBA(ssnum, dim*element.node_num, B_matrix, D_matrix, BtDB);
		mul_matrix_AtBA(dim*dim,dim*element.node_num,G_matrix,M_matrix,GtMG);

		for(ii2=0; ii2<(elem_matrix->size); ii2++){
			for(ii3=0; ii3<(elem_matrix->size); ii3++){
				elem_matrix->matrix[ii2][ii3] += (BtDB[ii2][ii3]+GtMG[ii2][ii3])
				                               * det_J*g_weights[ii1]*thickness;
			}
		}
	}

	return(NORMAL_RC);
}


/* Total Lagrange法の内力ベクトル(等価節点力ベクトル) */
RC
make_NL_react_force (ELEMENT_ARRAY element, NODE_ARRAY node,
                     MATERIAL_PROP_ARRAY material,
                     PHYSICAL_PROP_ARRAY physical,
                     DISP_ARRAY disp, double f_vect[])
{
	int dim, ssnum;
	int ii1, ii2, ii3;
	int total_dof;
	int g_num_point;
	double thickness;
	double det_J;
	double react_v[3*MAX_NODE];
	double PK2_stress_v[6];
	double Green_strain_v[6];
	double disp_v[3*MAX_NODE];
	double *g_weights = NULL;
	double **D_matrix = NULL;
	double **B_matrix = NULL;
	double **A_matrix = NULL;
	double **G_matrix = NULL;
	double **AG = NULL;
	double **B_strain = NULL;
	double **node_location = NULL;
	double **g_points = NULL;
	WORK_SET ws;

	if( (dim = analysis_dim(element)) <= 1 || dim > 3){
		return(ARG_ERROR_RC);
	}

	RC_TRY( allocate_work_set(&ws) );
	node_location = ws.mat_N_3[0];
	g_weights = ws.vec_G[0];
	D_matrix = ws.mat_6_6[0];
	B_matrix = ws.mat_6_3N[0];
	A_matrix = ws.mat_9_9[0];
	G_matrix = ws.mat_9_3N[0];
	AG = ws.mat_6_3N[1];
	B_strain = ws.mat_6_3N[2];
	g_points = ws.mat_G_4[0];

	RC_TRY( fem_total_dof(node, element, &total_dof) );

	for(ii1=0; ii1<total_dof; ii1++) f_vect[ii1] = 0.0; 

	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		if( (ssnum = element_ssnum(element.array[ii1].type)) <= 0 ){
			return(ARG_ERROR_RC);
		}

		RC_TRY( Gauss_const_default(element.array[ii1].type, g_points,
		                            g_weights, &g_num_point) );
		RC_TRY( fill_disp_vector(element.array[ii1], disp, disp_v) );
		RC_TRY( node_location_matrix(&(element.array[ii1]), node,
		                             node_location) );
		RC_TRY( get_D_matrix(element.array[ii1],material,physical,D_matrix) );
		RC_TRY( get_thickness(element.array[ii1], physical, &thickness) );

		for(ii2=0; ii2<g_num_point; ii2++){
			STRESS_STRAIN Green_strain;

			RC_TRY( make_B_matrix(element.array[ii1], g_points[ii2],
			                      node_location, B_matrix, &det_J) );
			RC_TRY( make_AG_matrix(element.array[ii1], g_points[ii2], ws,
			                       node_location, disp, A_matrix, G_matrix) );
			mul_matrix(ssnum, dim*dim, dim*element.array[ii1].node_num,
			           A_matrix, G_matrix, AG);
			RC_TRY( add_matrix(ssnum, dim*element.array[ii1].node_num, 1.0,
			                   B_matrix, 0.5, AG, B_strain) );
			mul_matrix_vect(ssnum, dim*element.array[ii1].node_num, B_strain,
			                disp_v, Green_strain_v);
			RC_TRY( vect2stress_strain(ssnum, Green_strain_v, &Green_strain) );
			RC_TRY( add_matrix(ssnum, dim*element.array[ii1].node_num, 1.0,
			                   B_matrix, 1.0, AG, B_matrix) );

			if(dim==2){
				for(ii3=0; ii3<ssnum; ii3++){
					PK2_stress_v[ii3] = D_matrix[ii3][0]*Green_strain.v.x
					                  + D_matrix[ii3][1]*Green_strain.v.y
					                  + D_matrix[ii3][2]*Green_strain.v.xy;
				}
			}else{  /* dim == 3 */
				for(ii3=0; ii3<ssnum; ii3++){
					PK2_stress_v[ii3] = D_matrix[ii3][0]*Green_strain.v.x
					                  + D_matrix[ii3][1]*Green_strain.v.y
					                  + D_matrix[ii3][2]*Green_strain.v.z
					                  + D_matrix[ii3][3]*Green_strain.v.yz
					                  + D_matrix[ii3][4]*Green_strain.v.zx
					                  + D_matrix[ii3][5]*Green_strain.v.xy;
				}
			}

			mul_trans_matrix_vect(ssnum, dim*element.array[ii1].node_num,
			                      B_matrix, PK2_stress_v, react_v);

			for(ii3=0; ii3<element.array[ii1].node_num; ii3++){
				int index;

				index = node_renum_index(node, element.array[ii1].node[ii3]);
				if(index < 0) return(SEEK_ERROR_RC);

				f_vect[index*dim] += g_weights[ii2] * det_J * thickness
				                   * react_v[ii3*dim];
				f_vect[index*dim+1] += g_weights[ii2] * det_J * thickness
				                     * react_v[ii3*dim+1];
				if(dim == 3){
					f_vect[index*dim+2] += g_weights[ii2] * det_J * thickness
					                     * react_v[ii3*dim+2];
				}
			}
		}
	}

	RC_TRY( free_work_set(&ws) );

	return(NORMAL_RC);
}


/* 拘束点にゼロ代入 */
static RC
set_rest_zero (NODE_ARRAY node, BC_ARRAY rest, double resid_v[])
{
	int ii1;

	for(ii1=0; ii1<rest.size; ii1++){
		int index;

		if(rest.array[ii1].node < 0) continue;

		index = node_renum_index(node, rest.array[ii1].node);
		if(index < 0) return(SEEK_ERROR_RC);

		if(rest.array[ii1].v_type.x == BC_FIX) resid_v[3*index] = 0.0;
		if(rest.array[ii1].v_type.y == BC_FIX) resid_v[3*index+1] = 0.0;
		if(rest.array[ii1].v_type.z == BC_FIX) resid_v[3*index+2] = 0.0;
	}

	return(NORMAL_RC);
}


/* Cauchy応力，Almansi歪み */
RC
Cauchy_stress_Almansi_strain (int sol_ss_flag, NODE_ARRAY node,
                              ELEMENT_ARRAY element, DISP_ARRAY disp,
                              MATERIAL_PROP_ARRAY material,
                              PHYSICAL_PROP_ARRAY physical,
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

	RC_TRY( allocate_work_set(&ws) );
	local_points = ws.mat_N_4[0];

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
			RC_TRY(set_center_point(element.array[ii1].type, local_points[0]));
			RC_TRY( local_Cauchy_stress_Almansi_strain(element.array[ii1],
			                                           local_points[0],
			                                           node, disp, material,
			                                           physical, ws,
			                                           &local_stress,
			                                           &local_strain) );
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
			RC_TRY(set_local_node_points(element.array[ii1].type,local_points));

			for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
				RC_TRY( local_Cauchy_stress_Almansi_strain(element.array[ii1],
				                                           local_points[ii2],
				                                           node, disp, material,
				                                           physical, ws,
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


/* 第2Piola-Kirchhoff応力，Green歪み */
RC
PK2_stress_Green_strain (int sol_ss_flag, NODE_ARRAY node,
                         ELEMENT_ARRAY element, DISP_ARRAY disp,
                         MATERIAL_PROP_ARRAY material,
                         PHYSICAL_PROP_ARRAY physical,
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

	RC_TRY( allocate_work_set(&ws) );
	local_points = ws.mat_N_4[0];

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
			RC_TRY(set_center_point(element.array[ii1].type, local_points[0]));
			RC_TRY( local_PK2_stress_Green_strain(element.array[ii1],
			                                      local_points[0],
			                                      node, disp, material,
			                                      physical, ws,
			                                      &local_stress,
			                                      &local_strain) );
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
			RC_TRY(set_local_node_points(element.array[ii1].type,local_points));

			for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
				RC_TRY( local_PK2_stress_Green_strain(element.array[ii1],
				                                      local_points[ii2],
				                                      node, disp, material,
				                                      physical, ws,
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


/* 局所座標系でのCauchy応力，Almansi歪み */
RC
local_Cauchy_stress_Almansi_strain (ELEMENT element, const double *local_points,
                                    NODE_ARRAY node, DISP_ARRAY disp,
                                    MATERIAL_PROP_ARRAY material,
                                    PHYSICAL_PROP_ARRAY physical,
                                    WORK_SET ws,
                                    STRESS_STRAIN *stress,
                                    STRESS_STRAIN *strain)
{
	int ii1;
	int dim, ssnum;
	double stress_v[6];
	double **D_matrix = ws.mat_6_6[0];

	if(element.label < 0) return(ARG_ERROR_RC);

	if( (dim = element_dim(element.type)) <= 1) return(ARG_ERROR_RC);
	if( (ssnum = element_ssnum(element.type)) <= 0) return(ARG_ERROR_RC);

	RC_TRY( local_Almansi_strain(element,local_points,node,ws,disp,strain) );
	RC_TRY( get_D_matrix(element, material, physical, D_matrix) );

	if(dim == 2){
		for(ii1=0; ii1<ssnum; ii1++){
			stress_v[ii1] = D_matrix[ii1][0]*strain->v.x
			              + D_matrix[ii1][1]*strain->v.y
			              + D_matrix[ii1][2]*strain->v.xy;
		}
	}else{  /* dim == 3 */
		for(ii1=0; ii1<ssnum; ii1++){
			stress_v[ii1] = D_matrix[ii1][0]*strain->v.x
			              + D_matrix[ii1][1]*strain->v.y
			              + D_matrix[ii1][2]*strain->v.z
			              + D_matrix[ii1][3]*strain->v.yz
			              + D_matrix[ii1][4]*strain->v.zx
			              + D_matrix[ii1][5]*strain->v.xy;
		}
	}

	init_stress_strain(stress);
	stress->element = element.label;
	RC_TRY( vect2stress_strain(ssnum, stress_v, stress) );

	return(NORMAL_RC);
}


/* 局所座標系での第2Piola-Kirchhoff応力，Green歪み */
RC
local_PK2_stress_Green_strain (ELEMENT element, const double *local_points,
                               NODE_ARRAY node, DISP_ARRAY disp,
                               MATERIAL_PROP_ARRAY material,
                               PHYSICAL_PROP_ARRAY physical,
                               WORK_SET ws,
                               STRESS_STRAIN *stress,
                               STRESS_STRAIN *strain)
{
	int ii1;
	int dim, ssnum;
	double stress_v[6];
	double **D_matrix = ws.mat_6_6[0];

	if(element.label < 0) return(ARG_ERROR_RC);

	if( (dim = element_dim(element.type)) <= 1) return(ARG_ERROR_RC);
	if( (ssnum = element_ssnum(element.type)) <= 0) return(ARG_ERROR_RC);

	RC_TRY( local_Green_strain(element,local_points,node,ws,disp,strain) );
	RC_TRY( get_D_matrix(element, material, physical, D_matrix) );

	if(dim == 2){
		for(ii1=0; ii1<ssnum; ii1++){
			stress_v[ii1] = D_matrix[ii1][0]*strain->v.x
			              + D_matrix[ii1][1]*strain->v.y
			              + D_matrix[ii1][2]*strain->v.xy;
		}
	}else{  /* dim == 3 */
		for(ii1=0; ii1<ssnum; ii1++){
			stress_v[ii1] = D_matrix[ii1][0]*strain->v.x
			              + D_matrix[ii1][1]*strain->v.y
			              + D_matrix[ii1][2]*strain->v.z
			              + D_matrix[ii1][3]*strain->v.yz
			              + D_matrix[ii1][4]*strain->v.zx
			              + D_matrix[ii1][5]*strain->v.xy;
		}
	}

	init_stress_strain(stress);
	stress->element = element.label;
	RC_TRY( vect2stress_strain(ssnum, stress_v, stress) );

	return(NORMAL_RC);
}


/* 局所座標系でのAlmansi歪み */
RC
local_Almansi_strain (ELEMENT element, const double *local_point,
                      NODE_ARRAY node, WORK_SET ws, DISP_ARRAY disp,
                      STRESS_STRAIN *Almansi_strain)
{
	int dim, ii1, ii2;
	double **dN = ws.mat_3_N[0];
	double **dNxyz = ws.mat_3_N[1];
	double **Jacobi = ws.mat_3_3[0];
	double **inverse_J = ws.mat_3_3[1];
	double **node_location = ws.mat_N_3[1];
	double **update_location = ws.mat_N_3[2];
	double **inverse_F = ws.mat_3_3[2];
	double **inverse_Ft = ws.mat_3_3[3];
	double **strain_mat = ws.mat_6_6[1];

	RC_NEG_ZERO_CHK( dim = element_dim(element.type) );
	init_stress_strain(Almansi_strain);
	Almansi_strain->element = element.label;

	/* 初期配置での節点座標マトリクス */
	RC_TRY( node_location_matrix(&element, node, node_location) );

	/* 現配置での節点座標マトリクス */
	RC_TRY( update_node_location_matrix(&element, node, disp,
	                                    update_location) );

	RC_TRY( set_dN_matrix(element.type, local_point, dN) );
	mul_matrix(dim, element.node_num, dim, dN, update_location, Jacobi);
	RC_TRY( inverse_matrix(dim, Jacobi, inverse_J) );
	mul_matrix(dim, dim, element.node_num, inverse_J, dN, dNxyz);

	/* 逆変形勾配テンソル inverse_F, inverse_Ft */
	mul_matrix(dim, element.node_num, dim, dNxyz, node_location, inverse_Ft);
	transpose_matrix(dim, dim, inverse_Ft, inverse_F);
	mul_matrix(dim, dim, dim, inverse_Ft, inverse_F, strain_mat);

	/* Almansi歪みテンソル */
	for(ii1=0; ii1<dim; ii1++){
		for(ii2=0; ii2<dim; ii2++){
			strain_mat[ii1][ii2] *= -1.0;
		}
		strain_mat[ii1][ii1] += 1.0;
		for(ii2=0; ii2<dim; ii2++){
			strain_mat[ii1][ii2] *= 0.5;
		}
	}

	/* テンソル歪み -> 工学歪み */
	if(dim == 2){
		Almansi_strain->v.x = strain_mat[0][0];
		Almansi_strain->v.y = strain_mat[1][1];
		Almansi_strain->v.xy = strain_mat[0][1]*2.0;
	}else{ /* dim == 3 */
		Almansi_strain->v.x = strain_mat[0][0];
		Almansi_strain->v.y = strain_mat[1][1];
		Almansi_strain->v.z = strain_mat[2][2];
		Almansi_strain->v.yz = strain_mat[1][2]*2.0;
		Almansi_strain->v.zx = strain_mat[2][0]*2.0;
		Almansi_strain->v.xy = strain_mat[0][1]*2.0;
	}

	return(NORMAL_RC);
}


/* 局所座標系でのGreen歪み */
RC
local_Green_strain (ELEMENT element, const double *local_point,
                    NODE_ARRAY node, WORK_SET ws, DISP_ARRAY disp,
                    STRESS_STRAIN *Green_strain)
{
	int dim, ii1, ii2;
	double **dN = ws.mat_3_N[0];
	double **dNxyz = ws.mat_3_N[1];
	double **Jacobi = ws.mat_3_3[0];
	double **inverse_J = ws.mat_3_3[1];
	double **node_location = ws.mat_N_3[1];
	double **update_location = ws.mat_N_3[2];
	double **F = ws.mat_3_3[2];
	double **Ft = ws.mat_3_3[3];
	double **strain_mat = ws.mat_6_6[1];

	RC_NEG_ZERO_CHK( dim = element_dim(element.type) );
	init_stress_strain(Green_strain);
	Green_strain->element = element.label;

	/* 初期配置での節点座標マトリクス */
	RC_TRY( node_location_matrix(&element, node, node_location) );

	/* 現配置での節点座標マトリクス */
	RC_TRY( update_node_location_matrix(&element,node,disp,update_location) );

	RC_TRY( set_dN_matrix(element.type, local_point, dN) );
	mul_matrix(dim, element.node_num, dim, dN, node_location, Jacobi);
	RC_TRY( inverse_matrix(dim, Jacobi, inverse_J) );
	mul_matrix(dim, dim, element.node_num, inverse_J, dN, dNxyz);

	/* 変形勾配テンソル F，転置Ft */
	mul_matrix(dim, element.node_num, dim, dNxyz, update_location, Ft);
	transpose_matrix(dim, dim, Ft, F);
	mul_matrix(dim, dim, dim, Ft, F, strain_mat);

	/* Green歪みテンソル */
	for(ii1=0; ii1<dim; ii1++){
		strain_mat[ii1][ii1] -= 1.0;
		for(ii2=0; ii2<dim; ii2++){
			strain_mat[ii1][ii2] *= 0.5;
		}
	}

	/* テンソル歪み -> 工学歪み */
	if(dim == 2){
		Green_strain->v.x = strain_mat[0][0];
		Green_strain->v.y = strain_mat[1][1];
		Green_strain->v.xy = strain_mat[0][1]*2.0;
	}else{ /* dim == 3 */
		Green_strain->v.x = strain_mat[0][0];
		Green_strain->v.y = strain_mat[1][1];
		Green_strain->v.z = strain_mat[2][2];
		Green_strain->v.yz = strain_mat[1][2]*2.0;
		Green_strain->v.zx = strain_mat[2][0]*2.0;
		Green_strain->v.xy = strain_mat[0][1]*2.0;
	}

	return(NORMAL_RC);
}


/* 現配置の節点座標マトリクス */
static RC
update_node_location_matrix (ELEMENT *element, NODE_ARRAY node, DISP_ARRAY disp,
                             double **matrix)
{
	int ii1;

	for(ii1=0; ii1<element->node_num; ii1++){
		int node_index, disp_index;

		if(element->label < 0) continue;

		node_index = search_node_label(node, element->node[ii1]);
		disp_index = search_disp_node(disp, element->node[ii1]);
		if( (node_index < 0) || (disp_index < 0) ) return(SEEK_ERROR_RC);

		matrix[ii1][0] = node.array[node_index].p.x
		               + disp.array[disp_index].v.x;
		matrix[ii1][1] = node.array[node_index].p.y
		               + disp.array[disp_index].v.y;
		matrix[ii1][2] = node.array[node_index].p.z
		               + disp.array[disp_index].v.z;
	}

	return(NORMAL_RC);
}


/*
 * Reference book is
 * THE FINITE ELEMENT METHOD : Solid mechanics
 * 5th edition
 * O.C.Zienkiewicz, R.L.Taylor
 * 大変形の A,G マトリックス作成
 */
static RC
make_AG_matrix (ELEMENT element, const double *g_point, WORK_SET ws,
                double **node_location, DISP_ARRAY disp,
                double **A_matrix, double **G_matrix)
{
	int dim;
	double **dN = ws.mat_3_N[0];
	double **dNxyz = ws.mat_3_N[1];
	double **Jacobi = ws.mat_3_3[0];
	double **inverse_J = ws.mat_3_3[1];
	double **disp_matrix = ws.mat_N_6[0];
	double **duvw = ws.mat_3_3[2];


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


/*
 * Reference book is
 * THE FINITE ELEMENT METHOD : Solid mechanics
 * 5th edition
 * O.C.Zienkiewicz, R.L.Taylor
 * 大変形のGマトリックスのみ作成
 */
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


/* d_vectにdispのfactor倍を代入，ただしd_vectはメモリ確保済み */
RC
copy_disp2vect (int dim, NODE_ARRAY node, DISP_ARRAY disp, double d_vect[],
                double factor)
{
	int ii1;

	if( (dim != 2) && (dim != 3) && (dim != 6) ) return(ARG_ERROR_RC);
	if(disp.size <= 0) return(ARG_ERROR_RC);
	if(disp.source != FEM_K_SOL && disp.source != FEM_NASTRAN){
		return(ARG_ERROR_RC);
	}

	for(ii1=0; ii1<disp.size; ii1++){
		int index;

		if(disp.array[ii1].node < 0) continue;

		index = node_renum_index(node, disp.array[ii1].node);
		if(index < 0) return(SEEK_ERROR_RC);

		d_vect[index*dim] = factor * disp.array[ii1].v.x;
		d_vect[index*dim+1] = factor * disp.array[ii1].v.x;
		if(dim == 3) d_vect[index*dim+2] = factor * disp.array[ii1].v.z;
		if(dim == 6){
			d_vect[index*dim+3] = factor * disp.array[ii1].v.yz;
			d_vect[index*dim+4] = factor * disp.array[ii1].v.zx;
			d_vect[index*dim+5] = factor * disp.array[ii1].v.xy;
		}
	}

	return(NORMAL_RC);
}


/* dispにd_vectのfactor倍を代入(=)，ただしdispはデータ入力・メモリ確保済み */
RC
copy_vect2disp (int dim, NODE_ARRAY node, const double d_vect[],
                DISP_ARRAY disp, double factor)
{
	int ii1;

	if( (dim != 2) && (dim != 3) && (dim != 6) ) return(ARG_ERROR_RC);
	if(disp.size <= 0) return(ARG_ERROR_RC);
	if(disp.source != FEM_K_SOL && disp.source != FEM_NASTRAN){
		return(ARG_ERROR_RC);
	}

	for(ii1=0; ii1<disp.size; ii1++){
		int index;

		if(disp.array[ii1].node < 0) continue;

		index = node_renum_index(node, disp.array[ii1].node);
		if( (index < 0) || (index >= node.size) ) return(ARG_ERROR_RC);

		disp.array[ii1].v.x = factor * d_vect[dim*index];
		disp.array[ii1].v.y = factor * d_vect[dim*index+1];
		if(dim == 3) disp.array[ii1].v.z = factor * d_vect[dim*index+2];
		if(dim == 6){
			disp.array[ii1].v.yz = factor * d_vect[dim*index+3];
			disp.array[ii1].v.zx = factor * d_vect[dim*index+4];
			disp.array[ii1].v.xy = factor * d_vect[dim*index+5];
		}
	}


	return(NORMAL_RC);
}


/* d_vectにdispのfactor倍を加算(+=)，ただしd_vectはメモリ確保済み */
RC
add_disp2vect (int dim, NODE_ARRAY node, DISP_ARRAY disp, double d_vect[],
               double factor)
{
	int ii1;

	if( (dim != 2) && (dim != 3) && (dim != 6) ) return(ARG_ERROR_RC);
	if(disp.size <= 0) return(ARG_ERROR_RC);
	if(disp.source != FEM_K_SOL && disp.source != FEM_NASTRAN){
		return(ARG_ERROR_RC);
	}

	for(ii1=0; ii1<disp.size; ii1++){
		int index;

		if(disp.array[ii1].node < 0) continue;

		index = node_renum_index(node, disp.array[ii1].node);
		if(index < 0) return(SEEK_ERROR_RC);

		d_vect[index*dim] += factor * disp.array[ii1].v.x;
		d_vect[index*dim+1] += factor * disp.array[ii1].v.x;
		if(dim == 3) d_vect[index*dim+2] += factor * disp.array[ii1].v.z;
		if(dim == 6){
			d_vect[index*dim+3] += factor * disp.array[ii1].v.yz;
			d_vect[index*dim+4] += factor * disp.array[ii1].v.zx;
			d_vect[index*dim+5] += factor * disp.array[ii1].v.xy;
		}
	}

	return(NORMAL_RC);
}


/* dispにd_vectのfactor倍を加算(+=)，ただしdispはデータ入力・メモリ確保済み */
RC
add_vect2disp (int dim, NODE_ARRAY node, const double d_vect[], DISP_ARRAY disp,
               double factor)
{
	int ii1;

	if( (dim != 2) && (dim != 3) && (dim != 6) ) return(ARG_ERROR_RC);
	if(disp.size <= 0) return(ARG_ERROR_RC);
	if(disp.source != FEM_K_SOL && disp.source != FEM_NASTRAN){
		return(ARG_ERROR_RC);
	}

	for(ii1=0; ii1<disp.size; ii1++){
		int index;

		if(disp.array[ii1].node < 0) continue;

		index = node_renum_index(node, disp.array[ii1].node);
		if( (index < 0) || (index >= node.size) ) return(ARG_ERROR_RC);

		disp.array[ii1].v.x += factor * d_vect[dim*index];
		disp.array[ii1].v.y += factor * d_vect[dim*index+1];
		if(dim == 3) disp.array[ii1].v.z += factor * d_vect[dim*index+2];
		if(dim == 6){
			disp.array[ii1].v.yz += factor * d_vect[dim*index+3];
			disp.array[ii1].v.zx += factor * d_vect[dim*index+4];
			disp.array[ii1].v.xy += factor * d_vect[dim*index+5];
		}
	}

	return(NORMAL_RC);
}

