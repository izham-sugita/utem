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

/* $Id: eigen_solver.c,v 1.10 2003/07/22 09:15:27 nagatani Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "math_utl.h"
#include "nonzero_cg.h"
#include "rc.h"
#include "fem_solver.h"
  
#define ABS_ERROR (DBL_MIN / DBL_EPSILON)
#define REL_ERROR (DBL_EPSILON * 1.0e+3)
#define STACK_SIZE (1024)


static RC e_vals_sortSI(int e_num, double *e_vals);
static RC free_lanczosSI3s(int step, int total_dof,
                           double *alpha, double *beta, double **w_matrix);
static RC lanczosSI3s(NONZERO_MATRIX3 subA, NONZERO_MATRIX3 B, int max_step,
                      int *step, int *total_dof, double **alpha, double **beta,
                      double ***w_matrix);
static RC tri_e_vals_range(int n, const double *alpha, const double *beta,
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



RC fem_eigen_solve_iccg3(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                         FEM_MATERIAL_PROP_ARRAY material,
                         FEM_PHYSICAL_PROP_ARRAY physical,
                         double *evals, FEM_DISP_ARRAY *disps,
                         double sigma, int max_e_num, int *e_num, int max_step)
{
	NONZERO_MATRIX3 subK, M;
	int step;
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

	RC_TRY( fem_allocate_nonzero3(node, element, &subK) );
	RC_TRY( fem_allocate_nonzero3(node, element, &M) );

	for(ii1=0; ii1<element.size; ii1++){
		FEM_ELEM_MATRIX elem_matrix;

		if(element.array[ii1].label < 0) continue;

		RC_TRY( fem_make_elem_matrix(element.array[ii1], node,
		                             material, physical, &elem_matrix) );
		RC_TRY( fem_impose_nonzero3(subK, elem_matrix, 1.0) );
		RC_TRY( free_fem_elem_matrix(&elem_matrix) );

		RC_TRY( fem_make_mass_matrix(element.array[ii1], node, material,
		                             physical, &elem_matrix) );
		RC_TRY( fem_impose_nonzero3(subK, elem_matrix, -sigma) );
		RC_TRY( fem_impose_nonzero3(M, elem_matrix, 1.0) );
		RC_TRY( free_fem_elem_matrix(&elem_matrix) );
	}

	RC_TRY( lanczosSI3s(subK, M, max_step,
	                    &step, &total_dof, &alpha, &beta, &w) );
	RC_TRY( tri_e_vals_range(step, alpha, beta, &lambda_min, &lambda_max) );
	RC_TRY( sturm_stack(step, alpha, beta,
	              lambda_min, lambda_max, &tri_e_num, &tri_evals) );
	RC_TRY( e_vals_sortSI(tri_e_num, tri_evals) );

	(*e_num) = max_e_num;
	if((*e_num) > tri_e_num) (*e_num) = tri_e_num;
	RC_TRY( allocate2D(*e_num, step, &tri_evects) );
	RC_TRY( tri_inverse_iteration(step, *e_num, alpha, beta,
	                              tri_evals, tri_evects) );
	RC_NULL_CHK( disp_vect = (double *)malloc(sizeof(double) * total_dof) );
	for(ii1=0; ii1<(*e_num); ii1++){
		for(ii2=0; ii2<total_dof; ii2++){
			disp_vect[ii2] = 0.0;
		}
		for(ii2=0; ii2<total_dof; ii2++){
			for(ii3=0; ii3<step; ii3++){
				disp_vect[ii2] += w[ii3][ii2] * tri_evects[ii1][ii3];
			}
		}
		RC_TRY( fem_vect2disp(3, node, disp_vect, &(disps[ii1])) );
		evals[ii1] = 1.0/tri_evals[ii1] + sigma;
	}

	free((void *)disp_vect);
	RC_TRY( free_lanczosSI3s(step, total_dof, alpha, beta, w) );
	RC_TRY( free_nonzero_matrix3(&subK) );
	RC_TRY( free_nonzero_matrix3(&M) );
	free((void *)tri_evals);
	RC_TRY( free2D(*e_num, step, tri_evects) );

	return(NORMAL_RC);
}


static RC e_vals_sortSI(int e_num, double *e_vals)
{
	int ii1, ii2;
	double temp;
	
	for(ii1=0; ii1<e_num; ii1++){
		for(ii2=ii1+1; ii2<e_num; ii2++){
			if( (1.0/e_vals[ii2]) < (1.0/e_vals[ii1]) ){
				temp = e_vals[ii1];
				e_vals[ii1] = e_vals[ii2];
				e_vals[ii2] = temp;
			}
		}
	}

	return(NORMAL_RC);
}


static RC free_lanczosSI3s(int step, int total_dof,
                           double *alpha, double *beta, double **w_matrix)
{
	free((void *)alpha);
	free((void *)beta);
	RC_TRY( free2D(step, total_dof, w_matrix) );

	return(NORMAL_RC);
}


static RC lanczosSI3s(NONZERO_MATRIX3 subA, NONZERO_MATRIX3 B, int max_step,
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
	if(tmp_beta[-1] < ABS_ERROR) return(CAL_ERROR_RC);

	RC_NULL_CHK( v = (double **)malloc(max_step * sizeof(double *)) );
	RC_NULL_CHK( w = (double **)malloc((max_step+1) * sizeof(double *)) );
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

		if(tmp_beta[ii1] < -ABS_ERROR) return(CAL_ERROR_RC);
		if(tmp_beta[ii1] < ABS_ERROR) break;
	}

	free((void *)r);
	free((void *)q);

	RC_TRY( allocate1D(*step, alpha) );
	RC_TRY( allocate1D(*step, beta) );
	RC_NULL_CHK( (*w_matrix) = (double **)malloc((*step) * sizeof(double *)) );
	for(ii1=0; ii1<(*step); ii1++){
		(*alpha)[ii1] = tmp_alpha[ii1];
		(*beta)[ii1] = tmp_beta[ii1];
		(*w_matrix)[ii1] = w[ii1];
	}
	free((void *)tmp_alpha);
	tmp_beta--;
	free((void *)tmp_beta);
	free((void *)(w[-1]));
	free((void *)v);
	w--;
	free((void *)w);

	return (NORMAL_RC);
}


static RC tri_e_vals_range(int n, const double *alpha, const double *beta,
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


static RC sturm_stack(int n, double *alpha, double *beta,
                      double lambda_min, double lambda_max,
                      int *num_val, double **e_vals)
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
			   fabs(center)*REL_ERROR/10.0 + ABS_ERROR/10.0){
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
			   fabs(center)*REL_ERROR + ABS_ERROR){
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
static RC sturm(int n, double *alpha, double *beta,
                double lambda_min, double lambda_max,
                int *num_val, double **e_vals)
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


static RC lambda_number(int n, const double *alpha, const double *beta,
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


static RC tri_inverse_iteration(int step, int e_num,
                                const double *alpha, const double *beta,
                                double *e_vals, double **e_vects)
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
		while(sqrt(resid_norm) >= (REL_ERROR * fabs(theta) + ABS_ERROR)){
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

	free( (void *)y );
	free( (void *)sub_alpha );

	return (NORMAL_RC);
}


static RC tri_gauss(int n, const double *alpha, const double *beta,
                    const double *y, double *x)
{
	int ii1;
	double *temp;

	RC_TRY( allocate1D(n, &temp) );

	temp[0] = alpha[0];
	if(temp[0] >= 0.0){
		temp[0] += ABS_ERROR;
	}else{
		temp[0] -= ABS_ERROR;
	}

	x[0] = y[0];
	for(ii1=1; ii1<n; ii1++){
		temp[ii1] = alpha[ii1] - ((beta[ii1-1] * beta[ii1-1]) / temp[ii1-1]);
		if(fabs(temp[ii1]) >= 0.0){
			temp[ii1] += ABS_ERROR;
		}else{
			temp[ii1] -= ABS_ERROR;
		}

		x[ii1] = y[ii1] - ((beta[ii1-1] * x[ii1-1]) / temp[ii1-1]);
	}

	x[n-1] /= temp[n-1];
	for(ii1=n-2; ii1>=0; ii1--){
		x[ii1] = (x[ii1] - x[ii1+1] * beta[ii1]) / temp[ii1];
	}

	free( (void *)temp );

	return (NORMAL_RC);
}


