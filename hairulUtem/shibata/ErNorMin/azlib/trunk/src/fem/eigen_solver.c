/*********************************************************************
 * eigen_solver.c
 *
 * Copyright (C) 2002 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Taiki AOYAMA>
 *
 * Thanks to following contributors.
 *  <Tomoyuki TSUBATA> <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: eigen_solver.c 1091 2016-12-01 17:37:37Z wares $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "math_utl.h"
#include "nonzero_cg.h"
#include "rc.h"
#include "memory_manager.h"
#include "fem_solver.h"
#include "fem_struct.h"
#include "sky_cholesky.h"

#define STACK_SIZE (1024)


/* ランチョス法 */
static RC eigen_sortSI(int e_num, double *e_vals, double **col_vects);
static RC free_lanczosSI3s(int step, int total_dof,
                           double *alpha, double *beta, double **w_matrix);
static RC lanczosSI3s(NONZERO_MATRIX3 subA, NONZERO_MATRIX3 B, int max_step,
                      int *step, int *total_dof, double **alpha, double **beta,
                      double ***w_matrix);
static RC tri_evals_range(int n, const double *alpha, const double *beta,
                          double *lambda_min, double *lambda_max);
static RC sturm_stack(int n, double *alpha, double *beta,
                      double lambda_min, double lambda_max,
                      int *num_val, double **e_vals);
/*static RC sturm(int n, double *alpha, double *beta,
                double lambda_min, double lambda_max,
                int *num_val, double **e_vals);*/
static RC lambda_number(int n, const double *alpha, const double *beta,
                        double lambda, int *zero_flag, int *p);
static RC tri_inverse_iteration(int step, int e_num,
                                const double *alpha, const double *beta,
                                double *e_vals, double **e_vects);
static RC tri_gauss(int n, const double *alpha, const double *beta,
                    const double *y, double *x);

/* ブロックランチョス法 */
static RC block_lanczosSI_sky3(SKYLINE_MATRIX subA, NONZERO_MATRIX3 B,
                               NODE_ARRAY node,
                               BC_ARRAY rest,
                               int b_size, int max_step, int *step,
                               int *total_dof, double **T_matrix, int **d_index,
                               double ***q_matrix);
static RC get_initial_index(int b_size, int n, NODE_ARRAY node,
                            BC_ARRAY rest, int init_index[]);
static RC qr_decomp_BO(int n, int b_size, NONZERO_MATRIX3 B,
                       double **R, double **q, double **w, double **B_mat);
static RC modify_qr_decomp_BO(int n, int b_size, NONZERO_MATRIX3 B,
                              double **R, double **q, double **w,
                              double **B_mat);
static RC tri_reduction(int n, int b_size, const double T_array[],
                        const int index[], double **alpha, double **beta,
                        double ***g_matrix);
static RC tri_ql_implicit(int n, const double alpha[], const double beta[],
                                       double e_vals[], double **ql_vects);
static double cal_shift_ql(double a1, double a2, double b);



/*
 * シフト付きLanczos法(反復解法)による自由振動固有値解析ソルバ
 * シフトパラメータ(shift)近傍の固有値を求める．ただし，3次元問題のみに対応
 * また，多重固有値の導出不可
 */
RC
fem_eigen_solve_iccg3 (NODE_ARRAY node, ELEMENT_ARRAY element,
                       MATERIAL_PROP_ARRAY material,
                       PHYSICAL_PROP_ARRAY physical,
                       double shift, int max_step,
                       int max_e_num, int *e_num,
                       double evals[], DISP_ARRAY mode_disps[])
{
	NONZERO_MATRIX3 subK, M;
	int step = 0;
	double *alpha;
	double *beta;
	double **w;
	double lambda_min, lambda_max;
	double *tri_evals;
	double **tri_evects;
	int tri_e_num;
	double *disp_vect;
	int ii1, ii2, ii3;
	int total_dof;

	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

	RC_TRY( allocate_nonzero3(node, element, &subK) );
	RC_TRY( allocate_nonzero3(node, element, &M) );



	for(ii1=0; ii1<element.size; ii1++){
		ELEM_MATRIX elem_matrix;

		if(element.array[ii1].label < 0) continue;

		RC_TRY( make_elem_matrix(element.array[ii1], node,
		                             material, physical, &elem_matrix) );
		RC_TRY( impose_nonzero3(subK, elem_matrix, 1.0) );
		RC_TRY( free_elem_matrix(&elem_matrix) );

		RC_TRY( make_mass_matrix(element.array[ii1], node, material,
		                             physical, &elem_matrix) );

		RC_TRY( impose_nonzero3(subK, elem_matrix, -shift) );
		RC_TRY( impose_nonzero3(M, elem_matrix, 1.0) );
		RC_TRY( free_elem_matrix(&elem_matrix) );
	}

	RC_TRY( lanczosSI3s(subK, M, max_step,
	                    &step, &total_dof, &alpha, &beta, &w) );
	RC_TRY( tri_evals_range(step, alpha, beta, &lambda_min, &lambda_max) );
	RC_TRY( sturm_stack(step, alpha, beta,
	              lambda_min, lambda_max, &tri_e_num, &tri_evals) );
	RC_TRY( eigen_sortSI(tri_e_num, tri_evals, NULL) );



	(*e_num) = max_e_num;
	if((*e_num) > tri_e_num) (*e_num) = tri_e_num;
	RC_TRY( allocate2D(*e_num, step, &tri_evects) );
	RC_TRY( tri_inverse_iteration(step, *e_num, alpha, beta,
	                              tri_evals, tri_evects) );
	RC_NULL_CHK( disp_vect = (double *)mm_alloc(sizeof(double) * total_dof) );
	for(ii1=0; ii1<(*e_num); ii1++){
		for(ii2=0; ii2<total_dof; ii2++){
			disp_vect[ii2] = 0.0;
		}
		for(ii2=0; ii2<total_dof; ii2++){
			for(ii3=0; ii3<step; ii3++){
				disp_vect[ii2] += w[ii3][ii2] * tri_evects[ii1][ii3];
			}
		}
		RC_TRY( vect2disp(3, node, disp_vect, &(mode_disps[ii1])) );
		evals[ii1] = 1.0/tri_evals[ii1] + shift;
	}



	RC_TRY( mm_free((void *)disp_vect) );
	RC_TRY( free_lanczosSI3s(step, total_dof, alpha, beta, w) );
	RC_TRY( free_nonzero_matrix3(&subK) );
	RC_TRY( free_nonzero_matrix3(&M) );
	RC_TRY( mm_free((void *)tri_evals) );
	RC_TRY( free2D(*e_num, step, &tri_evects) );

	return(NORMAL_RC);
}


/* 固有値，固有ベクトルを昇順にソーティング */
static RC
eigen_sortSI (int e_num, double *e_vals, double **col_vects)
{
	int ii1, ii2;
	
	for(ii1=0; ii1<e_num; ii1++){
		for(ii2=ii1+1; ii2<e_num; ii2++){
			double temp;
			double *temp_ptr;
			if( (1.0/e_vals[ii2]) < (1.0/e_vals[ii1]) ){
				temp = e_vals[ii1];
				e_vals[ii1] = e_vals[ii2];
				e_vals[ii2] = temp;
				if(col_vects != NULL){
					temp_ptr = col_vects[ii1];
					col_vects[ii1] = col_vects[ii2];
					col_vects[ii2] = temp_ptr;
				}
			}
		}
	}

	return(NORMAL_RC);
}


static RC
free_lanczosSI3s (int step, int total_dof,
                  double *alpha, double *beta, double **w_matrix)
{
	int ii1;

	for(ii1=0; ii1<step; ii1++){
		RC_TRY( mm_free((void *)w_matrix[ii1]) );
	}
	RC_TRY( mm_free((void *)w_matrix) );
	RC_TRY( mm_free((void *)alpha) );
	RC_TRY( mm_free((void *)beta) );

	return(NORMAL_RC);
}


static RC
lanczosSI3s (NONZERO_MATRIX3 subA, NONZERO_MATRIX3 B, int max_step,
             int *step, int *total_dof, double **alpha, double **beta,
             double ***w_matrix)
{
	int ii1, ii2, ii3;
	int n;
	double *r;
	double *q;
	double *tmp_alpha;
	double *tmp_beta;
	double **v;
	double **w;

	if( (subA.size != B.size)||(subA.size <= 0)||(B.size <= 0) ){
		return(ARG_ERROR_RC);
	}
	(*total_dof) = n = 3 * subA.size;

	RC_TRY( allocate1D(n, &r) );
	r[n-1] = 1.0;

	RC_TRY( allocate1D(n, &q) );
	RC_TRY( mul_nonzero_matrix_vect3s(B, r, q) );

	RC_TRY( allocate1D(max_step, &tmp_alpha) );
	RC_TRY( allocate1D(max_step+1, &tmp_beta) );
	tmp_beta++;
	tmp_beta[-1] = sqrt(inner_product(n, q, r));
	if(tmp_beta[-1] < ABS_TOL) return(CAL_ERROR_RC);

	RC_NULL_CHK( v = (double **)mm_alloc(max_step * sizeof(double *)) );
	RC_NULL_CHK( w = (double **)mm_alloc((max_step+1) * sizeof(double *)) );
	w++;
	RC_TRY( allocate1D(n, &(w[-1])) );


	for(ii1=0; ii1<max_step; ii1++){
		RC_TRY( allocate1D(n, &(w[ii1])) );
		RC_TRY( allocate1D(n, &(v[ii1])) );
		for(ii2=0; ii2<n; ii2++){
			w[ii1][ii2] = r[ii2] / tmp_beta[ii1-1];
			v[ii1][ii2] = q[ii2] / tmp_beta[ii1-1];
		}

		for(ii2=0; ii2<n; ii2++){
			r[ii2] = 0.0;
		}
		RC_TRY( iccg_nonzero3s(subA, v[ii1], r, 1, 1) );
		for(ii2=0; ii2<n; ii2++){
			r[ii2] -= tmp_beta[ii1-1] * w[ii1-1][ii2];
		}

		tmp_alpha[ii1] = inner_product(n, v[ii1], r);

		for(ii2=0; ii2<n; ii2++){
			r[ii2] -= tmp_alpha[ii1] * w[ii1][ii2];
		}

		/* reorthogonalize */
		for(ii2=0; ii2<=ii1; ii2++){
			double v_r = inner_product(n, v[ii2], r);
			for(ii3=0; ii3<n; ii3++){
				r[ii3] -= v_r * w[ii2][ii3];
			}
		}

		RC_TRY( mul_nonzero_matrix_vect3s(B, r, q) );

		tmp_beta[ii1] = sqrt(inner_product(n, q, r));
		*step = ii1+1;

		if(tmp_beta[ii1] < -ABS_TOL) return(CAL_ERROR_RC);
		if(tmp_beta[ii1] < ABS_TOL) break;
	}

	RC_TRY( mm_free((void *)r) );
	RC_TRY( mm_free((void *)q) );

	RC_TRY( allocate1D(*step, alpha) );
	RC_TRY( allocate1D(*step, beta) );
	RC_NULL_CHK( *w_matrix = (double **)mm_alloc((*step) * sizeof(double *)) );
	for(ii1=0; ii1<(*step); ii1++){
		(*alpha)[ii1] = tmp_alpha[ii1];
		(*beta)[ii1] = tmp_beta[ii1];
		(*w_matrix)[ii1] = w[ii1];
	}
	RC_TRY( mm_free((void *)tmp_alpha) );
	tmp_beta--;
	RC_TRY( mm_free((void *)tmp_beta) );
	RC_TRY( mm_free((void *)(w[-1])) );
	RC_TRY( mm_free((void *)v) );
	w--;
	RC_TRY( mm_free((void *)w) );

	return (NORMAL_RC);
}


static RC
tri_evals_range (int n, const double *alpha, const double *beta,
                  double *lambda_min, double *lambda_max)
{
	int ii1;
	double sum = 0.0;
	
	for(ii1=0; ii1<n-1; ii1++){
		sum +=  (alpha[ii1] * alpha[ii1]) + 2*(beta[ii1] * beta[ii1]);
	}
	sum += (alpha[n-1] * alpha[n-1]);
	sum = sqrt(sum);

	*lambda_max = sum;
	*lambda_min = -sum;

	return(NORMAL_RC);
}


static RC
sturm_stack (int n, double *alpha, double *beta, double lambda_min,
             double lambda_max, int *num_val, double **e_vals)
{
	int p1, p2;
	int e_count;
	double center;
	double max_stack[STACK_SIZE];
	double min_stack[STACK_SIZE];
	int stack_ptr = 0;
	int hit_flag1, hit_flag2;

	max_stack[stack_ptr] = lambda_max;
	min_stack[stack_ptr] = lambda_min;
	stack_ptr++;
	lambda_number(n, alpha, beta, lambda_min, &hit_flag1, &p1);
	lambda_number(n, alpha, beta, lambda_max, &hit_flag2, &p2);
	(*num_val) = p1 - p2;
	if(*num_val <= 0){
		fprintf(stderr, "num_val <= 0\n");
		fprintf(stderr, "num_val = %d\n", *num_val);
		return(CAL_ERROR_RC);
	}
	
	RC_TRY( allocate1D(*num_val, e_vals) );

	e_count = 0;
	while( (e_count < (*num_val)) && (stack_ptr >= 1) ){
		stack_ptr--;
		lambda_max = max_stack[stack_ptr];
		lambda_min = min_stack[stack_ptr];
		center = (lambda_min + lambda_max)/2.0;
		lambda_number(n, alpha, beta, lambda_min, &hit_flag1, &p1);
		lambda_number(n, alpha, beta, lambda_max, &hit_flag2, &p2);
		
		if(p1 - p2 > 1){ 
			if(fabs(lambda_max - lambda_min) < 
			   fabs(center)*REL_TOL/10.0 + ABS_TOL/10.0){
				fprintf(stderr, "!");
				(*e_vals)[e_count] = center;
				e_count++;
			}else{
				if(stack_ptr >= STACK_SIZE) return(OVERFLOW_ERROR_RC);
				max_stack[stack_ptr] = lambda_max;
				min_stack[stack_ptr] = center;
				stack_ptr++;
				if(stack_ptr >= STACK_SIZE) return(OVERFLOW_ERROR_RC);
				max_stack[stack_ptr] = center;
				min_stack[stack_ptr] = lambda_min;
				stack_ptr++;
			}
		}else if(p1 - p2 == 0){
			if(hit_flag2){
				(*e_vals)[e_count] = center;
				e_count++;
			}
		}else if(p1 - p2 == 1){
			if(fabs(lambda_max - lambda_min) < 
			   fabs(center)*REL_TOL + ABS_TOL){
				(*e_vals)[e_count] = center;
				e_count++;
			}else{
				if(stack_ptr >= STACK_SIZE) return(OVERFLOW_ERROR_RC);
				max_stack[stack_ptr] = lambda_max;
				min_stack[stack_ptr] = center;
				stack_ptr++;
				if(stack_ptr >= STACK_SIZE) return(OVERFLOW_ERROR_RC);
				max_stack[stack_ptr] = center;
				min_stack[stack_ptr] = lambda_min;
				stack_ptr++;
			}
		}else{
			fprintf(stderr, "(p1 - p2) < 0\n");
			return (CAL_ERROR_RC);
		}
	}

	return (NORMAL_RC);
}


#if 0
static RC
sturm (int n, double *alpha, double *beta, double lambda_min,
       double lambda_max, int *num_val, double **e_vals)
{
	int p1, p2;
	int e_count;
	double center;
	double e_max = lambda_max;
	int hit_flag1, hit_flag2;
int loop = 0;

	lambda_number(n, alpha, beta, lambda_min, &hit_flag1, &p1);
	lambda_number(n, alpha, beta, lambda_max, &hit_flag2, &p2);
	(*num_val) = p1 - p2;
	if(*num_val <= 0){
		fprintf(stderr, "num_val <= 0\n");
		fprintf(stderr, "num_val = %d\n", *num_val);
		return(CAL_ERROR_RC);
	}
	
	RC_TRY( allocate1D(*num_val, e_vals) );

	e_count = 0;
	while(e_count < (*num_val)){
loop++;
		center = (lambda_min + lambda_max) * 0.5;
		lambda_number(n, alpha, beta, lambda_min, &hit_flag1, &p1);
		lambda_number(n, alpha, beta, center, &hit_flag2, &p2);
fprintf(stderr, "%d:%d[%20.12e:%20.12e]:%d\n", loop, e_count,
lambda_min, lambda_max, p1-p2);
		
		if(p1 - p2 > 1){ 
			if(fabs(lambda_min - lambda_max) < 
			   fabs(center)*REL_ERROR/10.0 + ABS_ERROR/10.0){
				return(CAL_ERROR_RC);
			}
			lambda_max = center;
		}else if(p1 - p2 == 0){
			if(hit_flag2){
				(*e_vals)[e_count] = center;
				e_count++;
			}
			lambda_min = center;
			lambda_max = e_max;
		}else if(p1 - p2 == 1){
			if(fabs(lambda_min - lambda_max) < 
			   fabs(center)*REL_ERROR + ABS_ERROR){
				(*e_vals)[e_count] = center;
				e_count++;
				lambda_min = center;
				lambda_max = e_max;
			}else{
				lambda_max = center;
			}
		}else{
			fprintf(stderr, "(p1 - p2) < 0\n");
			return (CAL_ERROR_RC);
		}
	}

	return (NORMAL_RC);
}
#endif


static RC
lambda_number (int n, const double *alpha, const double *beta,
               double lambda, int *zero_flag, int *p)
{
	int ii1;
	double sturm_q0 = 0.0;
	double sturm_q1 = 0.0;
	double sturm_q2 = 0.0;

	if(n < 1) return(ARG_ERROR_RC);

	(*p) = 0;
	sturm_q2 = alpha[0] - lambda;
	if(sturm_q2 > 0.0) (*p)++;

	for(ii1=1; ii1<n; ii1++){
		sturm_q0 = sturm_q1;
		sturm_q1 = sturm_q2;
		if( nearly_eq(sturm_q1, 0.0) ){
			sturm_q2 = - 1.0;
		}else if( nearly_eq(sturm_q0, 0.0) && (ii1 > 1) ){
			sturm_q2 = alpha[ii1] - lambda;
		}else{
			sturm_q2 = alpha[ii1] - lambda
			          - ( (beta[ii1-1] * beta[ii1-1])/sturm_q1 );
		}
		if(sturm_q2 > 0.0) (*p)++;
	}

	if( nearly_eq(sturm_q2, 0.0) ){
		*zero_flag = 1;
	}else{
		*zero_flag = 0;
	}

	return(NORMAL_RC);
}


/* 3重対角行列用逆反復法 */
static RC
tri_inverse_iteration (int step, int e_num, const double *alpha,
                       const double *beta, double *e_vals, double **e_vects)
{
	int ii1, ii2;
	double theta, resid_norm, y_y;
	double *y;
	double *sub_alpha;

	RC_TRY( allocate1D(step, &y) );
	RC_TRY( allocate1D(step, &sub_alpha) );

	for(ii1=0; ii1<e_num; ii1++){
		for(ii2=0; ii2<step; ii2++){
			sub_alpha[ii2] = alpha[ii2] - e_vals[ii1];
			y[ii2] = 0.0;
		}
		y[step-1] = 1.0;

		resid_norm = 1.0;
		theta = 0.0;
		while(sqrt(resid_norm) >= (REL_TOL * fabs(theta) + ABS_TOL)){
			y_y = sqrt(inner_product(step, y, y));
			if(nearly_eq(y_y, 0.0)) return(CAL_ERROR_RC);

			for(ii2=0; ii2<step; ii2++){
				e_vects[ii1][ii2] = y[ii2] / y_y;
			}
			
			RC_TRY( tri_gauss(step, sub_alpha, beta, e_vects[ii1], y) );

			theta = inner_product(step, e_vects[ii1], y);

			resid_norm = 0.0;
			for(ii2=0; ii2<step; ii2++){
				resid_norm += (y[ii2] - theta * e_vects[ii1][ii2])
				             *(y[ii2] - theta * e_vects[ii1][ii2]);
			}
		}	
		e_vals[ii1] += 1.0 / theta; 
	}

	RC_TRY( mm_free((void *)y) );
	RC_TRY( mm_free((void *)sub_alpha) );

	return (NORMAL_RC);
}


/* 3重対角行列用Gaussの消去法 */
static RC
tri_gauss (int n, const double *alpha, const double *beta,
           const double *y, double *x)
{
	int ii1;
	double *temp;

	RC_TRY( allocate1D(n, &temp) );

	temp[0] = alpha[0];
	if(temp[0] >= 0.0){
		temp[0] += ABS_TOL;
	}else{
		temp[0] -= ABS_TOL;
	}

	x[0] = y[0];
	for(ii1=1; ii1<n; ii1++){
		temp[ii1] = alpha[ii1] - ((beta[ii1-1] * beta[ii1-1]) / temp[ii1-1]);
		if(fabs(temp[ii1]) >= 0.0){
			temp[ii1] += ABS_TOL;
		}else{
			temp[ii1] -= ABS_TOL;
		}

		x[ii1] = y[ii1] - ((beta[ii1-1] * x[ii1-1]) / temp[ii1-1]);
	}

	x[n-1] /= temp[n-1];
	for(ii1=n-2; ii1>=0; ii1--){
		x[ii1] = (x[ii1] - x[ii1+1] * beta[ii1]) / temp[ii1];
	}

	RC_TRY( mm_free((void *)temp) );

	return (NORMAL_RC);
}


/* 
 * ここからブロックランチョス(Block Lanczos)法プログラム
 */

/*
 * Block Lanczos法による固有値，固有モード解析(3次元solid)
 * 連立方程式の求解にはスカイコレスキー法使用．
 * ブロックサイズ(b_size)以下の多重度の導出可能
 */
RC
fem_eigen_block_solve_sky3 (NODE_ARRAY node, ELEMENT_ARRAY element,
                            MATERIAL_PROP_ARRAY material,
                            PHYSICAL_PROP_ARRAY physical,
                            BC_ARRAY rest,
                            double shift, int b_size, int max_step,
                            int max_e_num, int *e_num,
                            double evals[], DISP_ARRAY mode_disps[])
{
	SKYLINE_MATRIX subK;
	NONZERO_MATRIX3 M;
	int ii1, ii2, ii3;
	int dim;
	int total_dof;
	int step;
	int tri_size;
	int *d_index;
	double *T_array;
	double *alpha;
	double *beta;
	double *g0_ptr;
	double *tri_evals;
	double *disp_vect;
	double **q_matrix;
	double **g_matrix;


	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

	dim = analysis_dim(element);
	if(dim != 3){
		return(ARG_ERROR_RC);
	}

	RC_TRY( allocate_g_matrix(dim, node, element, &subK) );
	RC_TRY( allocate_nonzero3(node, element, &M) );

	/* K-sM and M */
	for(ii1=0; ii1<element.size; ii1++){
		ELEM_MATRIX elem_matrix;

		if(element.array[ii1].label < 0) continue;

		RC_TRY( make_elem_matrix(element.array[ii1], node,
		                             material, physical, &elem_matrix) );
		RC_TRY( impose_elem_matrix(subK, elem_matrix, 1.0) );
		RC_TRY( free_elem_matrix(&elem_matrix) );

		RC_TRY( make_mass_matrix(element.array[ii1], node, material,
		                             physical, &elem_matrix) );
		RC_TRY( impose_elem_matrix(subK, elem_matrix, -shift) );
		RC_TRY( impose_nonzero3(M, elem_matrix, 1.0) );
		RC_TRY( free_elem_matrix(&elem_matrix) );
	}

	/* Block Lanczos iteration */
	RC_TRY( block_lanczosSI_sky3(subK, M, node, rest, b_size, max_step, &step,
	                             &total_dof, &T_array, &d_index, &q_matrix) );
	RC_TRY( free_g_matrix(&subK) );
	RC_TRY( free_nonzero_matrix3(&M) );

	tri_size = step * b_size;

	/* 3重対角行列へ相似変換*/
	RC_TRY( tri_reduction(tri_size, b_size, T_array, d_index, &alpha, &beta,
	                      &g_matrix) );
	RC_TRY( mm_free(T_array) );
	RC_TRY( mm_free(d_index) );

	RC_TRY( allocate1D(tri_size, &tri_evals) );

	/* T_arrayの固有値，固有ベクトル計算 */
	RC_TRY( tri_ql_implicit(tri_size, alpha, beta, tri_evals, g_matrix) );
	RC_TRY( mm_free(alpha) );
	RC_TRY( mm_free(beta) );

	/* ソーティングに備えて先頭アドレス記憶 */
	g0_ptr = g_matrix[0];
	/* 昇順にソーティング */
	RC_TRY( eigen_sortSI(tri_size, tri_evals, g_matrix) );

	/* 実際に得られた固有値数(e_num) */
	(*e_num) = max_e_num;
	if((*e_num) > tri_size) (*e_num) = tri_size;

	RC_TRY( allocate1D(total_dof, &disp_vect) );
	for(ii1=0; ii1<(*e_num); ii1++){
		for(ii2=0; ii2<total_dof; ii2++){
			disp_vect[ii2] = 0.0;
		}
		/* Lanczosベクトルによる変換 */
		for(ii2=0; ii2<total_dof; ii2++){
			for(ii3=0; ii3<tri_size; ii3++){
				disp_vect[ii2] += q_matrix[ii3][ii2] * g_matrix[ii1][ii3];
			}
		}
		RC_TRY( vect2disp(3, node, disp_vect, &(mode_disps[ii1])) );

		/*変換前の固有値計算 */
		evals[ii1] = 1.0/tri_evals[ii1] + shift;
	}

	RC_TRY( mm_free(disp_vect) );
	RC_TRY( mm_free(tri_evals) );
	RC_TRY( free2D(tri_size, total_dof, &q_matrix) );
	/* 解放の為に先頭アドレスを戻す */
	g_matrix[0] = g0_ptr;
	RC_TRY( free2D(tri_size, tri_size, &g_matrix) );

	return(NORMAL_RC);
}


/* Block Lanczos法の反復計算 */
static RC
block_lanczosSI_sky3 (SKYLINE_MATRIX subA, NONZERO_MATRIX3 B,
                      NODE_ARRAY node, BC_ARRAY rest,
                      int b_size, int max_step, int *step,
                      int *total_dof, double **T_array, int **d_index, 
                      double ***q_matrix)
{
	int ii1, ii2, ii3, ii4, ii5;
	int n;
	int index;
	int *sub_index;
	int *init_index;
	double *init_v;
	double *sub_array;
	double **R;
	double **W;
	double **A_mat;
	double **B_mat;
	double ***Q;


	if((subA.dof/subA.dim != B.size)||(subA.dof/subA.dim <= 0)||(B.size <= 0)){
		return(ARG_ERROR_RC);
	}
	(*total_dof) = n = 3 * B.size;

	if(n < 1 * b_size) return(ARG_ERROR_RC);

	RC_TRY( allocate1I(max_step*b_size, &sub_index) );
	RC_TRY( allocate1D(((b_size+1)*b_size*(2*max_step - 1))/2, &sub_array) );
	RC_TRY( allocate2D(b_size, n, &R) );
	RC_TRY( allocate3D(max_step+2, b_size, n, &Q) );
	Q++;

	RC_TRY( allocate2D(b_size, n, &W) );
	RC_TRY( allocate2D(b_size, b_size, &(A_mat)) );
	RC_TRY( allocate2D(b_size, b_size, &(B_mat)) );

	/* Provisional block R[][] */
	RC_TRY( allocate1D(n, &init_v) );
	RC_TRY( allocate1I(b_size, &init_index) );
	RC_TRY( get_initial_index(b_size, n, node, rest, init_index) );
	for(ii1=0; ii1<b_size; ii1++){
		init_v[init_index[ii1]*3] = 1.0;
		RC_TRY( mul_nonzero_matrix_vect3s(B, init_v, R[ii1]) );
		init_v[init_index[ii1]*3] = 0.0;
	}
	/* Modify B */
	RC_TRY( modify_nonzero3(node, rest, B, init_v) );
	RC_TRY( mm_free(init_v) );
	RC_TRY( mm_free(init_index) );

	/* Modify subA */
	RC_TRY( modify_g_matrix_n(node, rest, subA, R, b_size) );

	RC_TRY( s_cholesky_decomp3(subA.dof,subA.index1,subA.index2,subA.array,1) );

	/* Initialization */
	/* Starting block R[][], Q[0] */
	for(ii1=0; ii1<b_size; ii1++){
		RC_TRY( s_cholesky_subst(subA.dof, subA.index1, subA.index2,
		                                   subA.array, R[ii1], 0) );
	}
	RC_TRY( modify_qr_decomp_BO(n, b_size, B, R, Q[0], W, B_mat) );

	/* Lanczos Loop */
	index = 0;
	*step = 0;
	for(ii1=0; ii1<max_step; ii1++){
		RC_TRY( progress_printf(5, "L[%d/%d]", ii1+1, max_step) );
		for(ii2=0; ii2<b_size; ii2++){
			for(ii3=0; ii3<n; ii3++){
				R[ii2][ii3] = W[ii2][ii3];
			}
			RC_TRY( s_cholesky_subst(subA.dof, subA.index1, subA.index2,
			                                   subA.array, R[ii2], 0) );
			for(ii3=0; ii3<b_size; ii3++){
				double *r_ii2 = R[ii2];
				double *q_ii1_1_ii3 = Q[ii1-1][ii3];
				double b_ii2_ii3 = B_mat[ii2][ii3];
				for(ii4=0; ii4<n; ii4++){
					r_ii2[ii4] -= q_ii1_1_ii3[ii4] * b_ii2_ii3;
				}
			}

			for(ii3=(ii2+1); ii3<b_size; ii3++){
				A_mat[ii2][ii3] = A_mat[ii3][ii2]
				                = inner_product(n, R[ii2], W[ii3]);
			}
			A_mat[ii2][ii2] = inner_product(n, R[ii2], W[ii2]);

			for(ii3=0; ii3<b_size; ii3++){
				double *r_ii2 = R[ii2];
				double *q_ii1_ii3 = Q[ii1][ii3];
				double a_ii2_ii3 = A_mat[ii2][ii3];
				for(ii4=0; ii4<n; ii4++){
					r_ii2[ii4] -= q_ii1_ii3[ii4] * a_ii2_ii3;	
				}
			}
		}

		RC_TRY( qr_decomp_BO(n, b_size, B, R, Q[ii1+1], W, B_mat) );

		/* Full reorthonalization */
		for(ii2=0; ii2<=ii1; ii2++){
			double **q_ii1_1 = Q[ii1+1];
			double **q_ii2 = Q[ii2];
			for(ii3=0; ii3<b_size; ii3++){
				double *w_ii3 = W[ii3];
				double *q_ii1_1_ii3 = q_ii1_1[ii3];
				for(ii4=0; ii4<b_size; ii4++){
					double *q_ii2_ii4 = q_ii2[ii4];
					double Qi_MQj1 = inner_product(n, q_ii2_ii4, w_ii3);
					for(ii5=0; ii5<n; ii5++){
						q_ii1_1_ii3[ii5] -= Qi_MQj1 * q_ii2_ii4[ii5];
					}
				}
			}
		}
		for(ii2=0; ii2<b_size; ii2++){
			RC_TRY(mul_nonzero_matrix_vect3s(B, Q[ii1+1][ii2], W[ii2]));
		}

		(*step)++;

		/* Rank check */ 
		if(-ABS_TOL > B_mat[b_size-1][b_size-1]) return(CAL_ERROR_RC);
		if( ABS_TOL > B_mat[b_size-1][b_size-1] || ii1 == (max_step-1) ){
			for(ii2=0; ii2<b_size; ii2++){
				sub_index[ii1*b_size+ii2] = index;
				for(ii3=ii2; ii3<b_size; ii3++){
					sub_array[index] = A_mat[ii2][ii3];
					index++;
				}
			}
			break;
		}

		/* Block tridiagonal -> T_array(1D array type) */
		for(ii2=0; ii2<b_size; ii2++){
			sub_index[ii1*b_size+ii2] = index;
			for(ii3=ii2; ii3<b_size; ii3++){
				sub_array[index] = A_mat[ii2][ii3];
				index++;
			}
			for(ii3=0; ii3<=ii2; ii3++){
				sub_array[index] = B_mat[ii3][ii2];
				index++;
			}
		}
	}
	RC_TRY( free2D(b_size, b_size, &A_mat) );
	RC_TRY( free2D(b_size, b_size, &B_mat) );
	RC_TRY( free2D(b_size, n, &R) );
	RC_TRY( free2D(b_size, n, &W) );

	RC_TRY( allocate1I((*step)*b_size, d_index) );
	RC_TRY( allocate1D(index, T_array) );
	for(ii1=0; ii1<index; ii1++){
		(*T_array)[ii1] = sub_array[ii1];
	}
	for(ii1=0; ii1<(*step); ii1++){
		for(ii2=0; ii2<b_size; ii2++){
			(*d_index)[ii1*b_size+ii2] = sub_index[ii1*b_size+ii2];
		}
	}
	RC_TRY( mm_free(sub_index) );
	RC_TRY( mm_free(sub_array) );

	RC_TRY( allocate2D((*step)*b_size, n, q_matrix) );
	RC_NULL_CHK( (*q_matrix)[0] = (double *)memcpy((*q_matrix)[0], Q[0][0],
	                              (*step)*b_size*n*sizeof(double)) );
	Q--;
	RC_TRY( free3D(max_step+2, b_size, n, &Q) );

	return(NORMAL_RC);
}


/* 初期値を代入するindexを求める */
static RC
get_initial_index (int b_size, int n, NODE_ARRAY node, BC_ARRAY rest,
                   int *init_index)
{
	int ii1, ii2, ii3;
	int *rest_index;

	if(rest.size > 0){
		RC_TRY( allocate1I(rest.size, &rest_index) );
		for(ii1=0; ii1<rest.size; ii1++){
			rest_index[ii1] = node_renum_index(node, rest.array[ii1].node);
		}

		ii2 = 0;
		for(ii1=0; ii1<b_size; ii1++){
			int index = ii2;
			for(ii2=index; ii2<n; ii2++){
				for(ii3=0; ii3<rest.size; ii3++){
					if(ii2 == rest_index[ii3]) break;
				}
				if(ii3 == rest.size){
					init_index[ii1] = ii2;
					ii2++;
					break;
				}
			}
		}
		RC_TRY( mm_free(rest_index) );
	}else if(rest.size == 0){
		for(ii1=0; ii1<b_size; ii1++){
			init_index[ii1]= ii1;
		}
	}else{
		return(NEGATIVE_ERROR_RC);
	}


	return(NORMAL_RC);
}


/* B-直交QR分解，Gram-Schmidt法使用 */
static RC
qr_decomp_BO (int n, int b_size, NONZERO_MATRIX3 B,
              double **R, double **q, double **w, double **B_mat)
{
	int ii1, ii2, ii3;

	for(ii1=0; ii1<b_size; ii1++){
		for(ii2=0; ii2<n; ii2++){
			q[ii1][ii2] = R[ii1][ii2];
		}
		RC_TRY( mul_nonzero_matrix_vect3s(B, q[ii1], w[ii1]) );
	}

	for(ii1=0; ii1<b_size; ii1++){
		double *q_ii1 = q[ii1];
		double *w_ii1 = w[ii1];
		double qw, b_ii1_ii1;

		qw = inner_product(n, q_ii1, w_ii1);
		if(qw < 0.0) return(NEGATIVE_ERROR_RC);
		b_ii1_ii1 = B_mat[ii1][ii1] = sqrt(qw);
		if( nearly_eq(b_ii1_ii1, 0.0) ) return(CAL_ERROR_RC);

		for(ii2=0; ii2<n; ii2++){
			q_ii1[ii2] /= b_ii1_ii1;
			w_ii1[ii2] /= b_ii1_ii1;
		}
		for(ii2=(ii1+1); ii2<b_size; ii2++){
			double *q_ii2 = q[ii2];
			double *w_ii2 = w[ii2];
			double b_ii1_ii2 = B_mat[ii1][ii2]
			                 = inner_product(n, q_ii1, w_ii2);
			for(ii3=0; ii3<n; ii3++){
				q_ii2[ii3] -= b_ii1_ii2 * q_ii1[ii3];
				w_ii2[ii3] -= b_ii1_ii2 * w_ii1[ii3];
			}
		}
	}

	return(NORMAL_RC);
}


/* 修正B-直交QR分解 */
static RC
modify_qr_decomp_BO (int n, int b_size, NONZERO_MATRIX3 B, double **R,
                     double **q, double **w, double **B_mat)
{
	int ii1, ii2, ii3, ii4;
	int i_num = 2 * b_size;
	double **old_B; 
	double **B_hat; 

	RC_TRY( allocate2D(b_size, b_size, &old_B) );
	RC_TRY( allocate2D(b_size, b_size, &B_hat) );

	/* 初期化 */
	for(ii1=0; ii1<b_size; ii1++){
		for(ii2=(ii1+1); ii2<b_size; ii2++){
			old_B[ii2][ii1] = old_B[ii1][ii2] = 0.0;
		}
		old_B[ii1][ii1] = 1.0;

		for(ii2=0; ii2<n; ii2++){
			q[ii1][ii2] = R[ii1][ii2];
		}
	}

	for(ii1=0; ii1<i_num; ii1++){
		for(ii2=0; ii2<b_size; ii2++){
			RC_TRY( mul_nonzero_matrix_vect3s(B, q[ii2], w[ii2]) );
		}
		for(ii2=0; ii2<b_size; ii2++){
			double qw, b_ii2_ii2, *q_ii2, *w_ii2;

			qw = inner_product(n, q[ii2], w[ii2]);
			if(qw < 0.0) return(NEGATIVE_ERROR_RC);

			b_ii2_ii2 = B_hat[ii2][ii2] = sqrt(qw);
			if( nearly_eq(b_ii2_ii2, 0.0) ) return(CAL_ERROR_RC);

			q_ii2 = q[ii2];
			w_ii2 = w[ii2];
			for(ii3=0; ii3<n; ii3++){
				q_ii2[ii3] /= b_ii2_ii2;
				w_ii2[ii3] /= b_ii2_ii2;
			}
			for(ii3=ii2+1; ii3<b_size; ii3++){
				double b_ii2_ii3 = B_hat[ii2][ii3]
				                 = inner_product(n, q[ii2], w[ii3]);
				double *q_ii3 = q[ii3];
				double *w_ii3 = w[ii3];
				for(ii4=0; ii4<n; ii4++){
					q_ii3[ii4] -= b_ii2_ii3 * q_ii2[ii4];
					w_ii3[ii4] -= b_ii2_ii3 * w_ii2[ii4];
				}
			}
		}
		
		mul_matrix(b_size, b_size, b_size, B_hat, old_B, B_mat);
		for(ii2=0; ii2<b_size; ii2++){
			old_B[ii2][ii2] = B_mat[ii2][ii2];
			for(ii3=ii2+1; ii3<b_size; ii3++){
				old_B[ii2][ii3] = B_mat[ii2][ii3];
			}
		}
	}

	RC_TRY( free2D(b_size, b_size, &old_B) );
	RC_TRY( free2D(b_size, b_size, &B_hat) );

	return(NORMAL_RC);
}


/*
 * Reduction of block tridiafonal matrix (by Givens rotation)
 * Reduced matrix is a scalar tridiagonal matrix
 * ブロック三重対角行列を三重対角行列へ変換, Givens法使用
 */
static RC
tri_reduction (int n, int b_size, const double T_array[],
               const int index[], double **alpha, double **beta,
               double ***g_matrix)
{
	int ii1;

	RC_TRY( allocate1D(n, alpha) );
	RC_TRY( allocate1D(n-1, beta) );
	RC_TRY( allocate2D(n, n, g_matrix) );
	for(ii1=0; ii1<n; ii1++) (*g_matrix)[ii1][ii1] = 1.0;
	
	if(b_size == 1){
		for(ii1=0; ii1<n-1; ii1++){
			(*alpha)[ii1] = T_array[index[ii1]];
			(*beta)[ii1] = T_array[index[ii1]+1];
		}
		(*alpha)[ii1] = T_array[index[ii1]];
	}else if(b_size > 1){
		int ii2, ii3;
		int size;
		double *array;

		size = (b_size+1)*n - (b_size+1)*b_size/2;
		RC_TRY( allocate1D(size, &array) );
		for(ii1=0; ii1<size; ii1++){
			array[ii1] = T_array[ii1];
		}
		for(ii1=0; ii1<n-2; ii1++){
			int end_index = index[ii1+1] - index[ii1] - 1 + ii1;
			for(ii2=end_index; ii2>ii1+1; ii2--){
				int I = ii2 - 1;
				int J = ii2;
				double term, c, s, temp1, temp2;

				if( nearly_eq(array[index[ii1]+J-ii1], 0.0) ) continue;

				term = sqrt(array[index[ii1]+I-ii1]*array[index[ii1]+I-ii1]
				          + array[index[ii1]+J-ii1]*array[index[ii1]+J-ii1]);

				c = array[index[ii1]+I-ii1] / term;
				s = array[index[ii1]+J-ii1] / term;

				temp1 = array[index[I]]*c*c + array[index[J]]*s*s
				             + 2.0*array[index[I]+J-I]*s*c;
				temp2 = array[index[I]]*s*s + array[index[J]]*c*c
				             - 2.0*array[index[I]+J-I]*s*c;

				array[index[I]+J-I] = (array[index[J]] - array[index[I]])*s*c
				                    +  array[index[I]+J-I]*(c*c - s*s);

				array[index[I]] = temp1;
				array[index[J]] = temp2;

				for(ii3=J+1; ii3<(index[I+1]-index[I]+I); ii3++){
					temp1 = array[index[I]+ii3-I]*c + array[index[J]+ii3-J]*s;

					array[index[J]+ii3-J] = - array[index[I]+ii3-I]*s
					                        + array[index[J]+ii3-J]*c;

					array[index[I]+ii3-I] = temp1;
				}
				if(J + b_size < n){
					temp2 = array[index[J]+ii3-J]*c;
				}

				for(ii3=I-1; ii3>ii1; ii3--){
					temp1 = - array[index[ii3]+I-ii3]*s
					        + array[index[ii3]+J-ii3]*c;

					array[index[ii3]+I-ii3] = array[index[ii3]+I-ii3]*c
					                        + array[index[ii3]+J-ii3]*s;

					array[index[ii3]+J-ii3] = temp1;
				}
				array[index[ii3]+I-ii3] = array[index[ii3]+I-ii3]*c
				                        + array[index[ii3]+J-ii3]*s;

				/* g_matrix */
				for(ii3=0; ii3<n; ii3++){
					double gI = (*g_matrix)[I][ii3];
					double gJ = (*g_matrix)[J][ii3];
					(*g_matrix)[I][ii3] = c*gI + s*gJ;
					(*g_matrix)[J][ii3] = c*gJ - s*gI;
				}

				/* fill-in 処理 */
				while(J + b_size < n){
					double nonzero_elem = array[index[J]+b_size]*s;

					array[index[J]+b_size] = temp2;

					I += b_size;
					J += b_size;

					term = sqrt(array[index[I-b_size]+b_size]
					          * array[index[I-b_size]+b_size]
					          + nonzero_elem*nonzero_elem);

					c = array[index[I-b_size]+b_size] / term;
					s = nonzero_elem / term;


					temp1 = array[index[I]]*c*c + array[index[J]]*s*s
					      + 2.0*array[index[I]+J-I]*s*c;
					temp2 = array[index[I]]*s*s + array[index[J]]*c*c
					      - 2.0*array[index[I]+J-I]*s*c;

					array[index[I]+J-I] = (array[index[J]]-array[index[I]])*s*c
					                    +  array[index[I]+J-I]*(c*c - s*s);

					array[index[I]] = temp1;
					array[index[J]] = temp2;

					for(ii3=J+1; ii3<(index[I+1]-index[I]+I); ii3++){
						temp1 = array[index[I]+ii3-I]*c
						      + array[index[J]+ii3-J]*s;

						array[index[J]+ii3-J] = - array[index[I]+ii3-I]*s
						                        + array[index[J]+ii3-J]*c;

						array[index[I]+ii3-I] = temp1;
					}
					if(J + b_size < n){
						temp2 = array[index[J]+ii3-J]*c;
					}

					for(ii3=I-1; ii3>I-b_size; ii3--){
						temp1 = - array[index[ii3]+I-ii3]*s
					            + array[index[ii3]+J-ii3]*c;

						array[index[ii3]+I-ii3] = array[index[ii3]+I-ii3]*c
						                        + array[index[ii3]+J-ii3]*s;

						array[index[ii3]+J-ii3] = temp1;
					}
					array[index[ii3]+I-ii3] = array[index[ii3]+I-ii3]*c
					                        + nonzero_elem*s;

					/* g_matrix */
					for(ii3=0; ii3<n; ii3++){
						double gI = (*g_matrix)[I][ii3];
						double gJ = (*g_matrix)[J][ii3];
						(*g_matrix)[I][ii3] = c*gI + s*gJ;
						(*g_matrix)[J][ii3] = c*gJ - s*gI;
					}
				}
			}
			(*alpha)[ii1] = array[index[ii1]];
			(*beta)[ii1] = array[index[ii1]+1];
		}

		(*alpha)[n-2] = array[index[n-2]];
		(*beta)[n-2] = array[index[n-2]+1];
		(*alpha)[n-1] = array[index[n-1]];
		RC_TRY( mm_free(array) );
	}else{
		return(ARG_ERROR_RC);
	}

	return(NORMAL_RC);
}


/*
 * The implicit QL method for the symmetric tridiagonal matrix
 * 陰的シフト付きQL法による対称3重対角行列の固有値，固有ベクトル計算
 * alpha[](対角要素)，beta[](副対角要素)は上書きされることはない
 * ql_vectsは単位行列に初期化しておく
 */
static RC
tri_ql_implicit (int n, const double alpha[], const double beta[],
                              double e_vals[], double **ql_vects)
{
	if(n == 1){
		e_vals[0] = alpha[0];
	}else if(n > 1){
		int ii1, ii2, ii2x, ii3;
		double *copy_alpha; 
		double *copy_beta; 

		RC_TRY( allocate1D(n, &copy_alpha) );
		RC_TRY( allocate1D(n, &copy_beta) );
		for(ii1=0; ii1<n-1; ii1++){
			copy_alpha[ii1] = alpha[ii1];
			copy_beta[ii1] = beta[ii1];
		}
		copy_alpha[n-1] = alpha[n-1];

		for(ii1=0; ii1<n-1; ii1++){
			int i_num = 0;
			do{
				double p = 0.0;
				double g = 0.0;
				double r = 0.0;
				for(ii2x=ii1; ii2x<n-1; ii2x++){
					double dd = fabs(copy_alpha[ii2x])+fabs(copy_alpha[ii2x+1]);
					if( nearly_eq(fabs(copy_beta[ii2x]+dd), dd) ){
						break;
					}
				}

				if(ii2x != ii1){
					double shift, c, s;
					int zero_flag;
					if(i_num++ > n) {
						log_printf(1, "Too many iterations in QL\n");
						return(CAL_ERROR_RC);
					}
					shift = cal_shift_ql(copy_alpha[ii1],
					                     copy_alpha[ii1+1],
					                     copy_beta[ii1]); 

					/* アンダーフローのチェックフラグ */
					zero_flag = 0;
					c = 1.0;
					s = 1.0;
					g = copy_alpha[ii2x] - shift;
					p = 0.0;
					for(ii2=ii2x-1; ii2>=ii1; ii2--){
						double f = s * copy_beta[ii2];
						double b = c * copy_beta[ii2];

						r = sqrt(f*f + g*g);
						copy_beta[ii2+1] = r;

						/* アンダーフロー回避 */
						if( nearly_eq(r, 0.0) ){
							copy_alpha[ii2+1] -= p;
							copy_beta[ii2x] = 0.0;
							zero_flag = 1;
							break;
						}

						s = f / r;
						c = g / r;
						g = copy_alpha[ii2+1] - p;
						r = (copy_alpha[ii2] - g)*s + 2.0*c*b;
						p = s * r;
						copy_alpha[ii2+1] = g + p;
						g = c*r - b;

						/* ql_vects (列ベクトルタイプ) */
						for(ii3=0; ii3<n; ii3++){
							double ql_v = ql_vects[ii2][ii3];
							double ql_v1 = ql_vects[ii2+1][ii3];
							ql_vects[ii2][ii3] = c*ql_v - s*ql_v1;
							ql_vects[ii2+1][ii3] = c*ql_v1 + s*ql_v;
						}
					}
					if( zero_flag && ii2 >= ii1) continue;
					copy_alpha[ii1] -= p;
					copy_beta[ii1] = g;
					copy_beta[ii2x] = 0.0;
				}
			}while(ii2x != ii1);
		}
		for(ii1=0; ii1<n; ii1++){
			e_vals[ii1] = copy_alpha[ii1];
		}
		RC_TRY( mm_free(copy_alpha) );
		RC_TRY( mm_free(copy_beta) );
	}else{
		return(ARG_ERROR_RC);
	}

	return(NORMAL_RC);
}


/* 2×2小行列のWilkinsonシフト計算 */
static double
cal_shift_ql (double a1, double a2, double b)
{
	double delta;
	double sign_delta;

	if( nearly_eq(a1, a2) ) return(a1 - fabs(b));

	delta = (a2 - a1) * 0.5;
	sign_delta = (delta >= 0.0) ? 1.0 : -1.0;

	/*
	return(a1 - (sign_delta*b*b)/( fabs(delta) + sqrt(delta*delta + b*b) ));
	*/
	return((a1 + a2)*0.5 - sign_delta*sqrt(delta*delta + b*b));
}

//MAC method added by WARES
RC
cal_mac(NODE_ARRAY node, ELEMENT_ARRAY element, DISP_ARRAY mode_disps[],
		                DISP_ARRAY init_disps[], int e_num, double evals[]) 
{
	int ii1, ii2, ii3;
	int total_dof;
	int dim;
	double **mac;
	int *order;
	DISP_ARRAY tmp_disp;

    dim = analysis_dim(element); 
    RC_TRY( fem_total_dof(node, element, &total_dof) ); 

    RC_TRY( allocate2D(e_num, e_num, &mac) ); 
    RC_TRY( allocate1I(e_num, &order) ); 


    // MAC値を計算する 
    for(ii1=0; ii1<e_num; ii1++){
        double *d_vect1;
        RC_TRY( allocate1D(total_dof, &d_vect1) );
        RC_TRY( copy_disp2vect(dim, node, mode_disps[ii1], d_vect1, 1.0) ); 

        for(ii2=0; ii2<e_num; ii2++){
            double *d_vect2;
            RC_TRY( allocate1D(total_dof, &d_vect2) ); 
            RC_TRY( copy_disp2vect(dim, node, init_disps[ii2], d_vect2, 1.0) ); 
	//FORMULAR OF MAC **
            mac[ii1][ii2] = 
                ( inner_product(total_dof, d_vect1, d_vect2) 
                  * inner_product(total_dof, d_vect1, d_vect2) )
                / ( inner_product(total_dof, d_vect1, d_vect1)
                        * inner_product(total_dof, d_vect2, d_vect2) );

            RC_TRY( free1D(total_dof, &d_vect2) );
        }
        RC_TRY( free1D(total_dof, &d_vect1) );
    }

    // 入れ替わり発生時にモードを入れ替える exchange mode 
    for(ii1=0; ii1<e_num; ii1++){
       int max = 0; 
        for(ii2=0; ii2<e_num; ii2++){
            if(mac[ii1][ii2] > mac[ii1][max]){ //TRUE condition of MAC checking
//mac[ii1][max] means first column of MAC matrix
                max = ii2; //those diagonal order
            }
        }
        order[ii1] = max; //diagonal order is order[ii1]
 
	}

    for(ii1=0; ii1<e_num; ii1++){
        int max;
        double tmp_e_val;
        if(order[ii1] != ii1){ //if diagonal order is not order, it means "FLASE condition"
#if 1
	printf("\n-----------------------MAC is ---------------------\n");
            for(ii2=0; ii2<e_num; ii2++){
                for(ii3=0; ii3<e_num; ii3++){
		 printf("|%f", mac[ii2][ii3]);
                }
                printf("|\n");
            }
#endif
            max = order[ii1];
            RC_TRY( allocate_disp_array(node.size, &tmp_disp) );
            RC_TRY( copy_disp_array(mode_disps[ii1], &tmp_disp) );
	    
	    RC_TRY( realloc_disp_array(&mode_disps[ii1]) );
            RC_TRY( copy_disp_array(mode_disps[max], &mode_disps[ii1]) );
	    
	    RC_TRY( realloc_disp_array(&mode_disps[max]) );
            RC_TRY( copy_disp_array(tmp_disp, &mode_disps[max]) );
            RC_TRY(free_disp_array(&tmp_disp)); 
			
            order[max] = order[ii1];

            tmp_e_val = evals[max]; 
            evals[max] = evals[ii1];
            evals[ii1] = tmp_e_val;   

            printf("\n---------------- chenge_mode ---------------------\n");
            RC_TRY( log_printf(5, "change mode %d not equal to mode %d \n", order[ii1], ii1) );
            break;
#if 0
            double *d_vect1;
            double *d_vect2;
            RC_TRY( allocate1D(total_dof, &d_vect1) );
            RC_TRY( allocate1D(total_dof, &d_vect2) );
            RC_TRY( copy_disp2vect(dim, node, mode_disps[ii1], d_vect1, 1.0) ); 
            RC_TRY( copy_disp2vect(dim, node, mode_disps[max], d_vect2, 1.0) ); 
            RC_TRY( copy_vect2disp(dim, node, d_vect1, mode_disps[max], 1.0) );
            RC_TRY( copy_vect2disp(dim, node, d_vect2, mode_disps[ii1], 1.0) );

            RC_TRY( free1D(total_dof, &d_vect1) );
            RC_TRY( free1D(total_dof, &d_vect2) );
            break;
#endif
        }
	else{
	printf("\n---------Mode shape did not have problems. MAC is----------------\n");
	 	for(ii2=0; ii2<e_num; ii2++){
                	for(ii3=0; ii3<e_num; ii3++){
                    	printf("|%f", mac[ii2][ii3]);
                	}
                	printf("|\n");
            	}
	break;	
	}

    }

	RC_TRY( free2D(e_num, e_num, &mac) );

	RC_TRY( free1I(e_num, &order) );

	return(NORMAL_RC);
}
