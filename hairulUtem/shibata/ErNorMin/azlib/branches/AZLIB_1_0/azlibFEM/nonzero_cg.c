/*********************************************************************
 * nonzero_cg.c
 *
 * Copyright (C) 2001 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Tomoyuki TSUBATA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: nonzero_cg.c,v 1.18 2003/07/22 10:30:40 aoyama Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "nonzero_cg.h"
#include "math_utl.h"
#include "rc.h"

#define ABS_ERROR     (DBL_MIN/DBL_EPSILON)

#define REL_ERROR1     (DBL_EPSILON * 1.0e+2)
#define REL_ERROR2     (DBL_EPSILON * 1.0e+4)


static RC clone_nonzero_matrix3(NONZERO_MATRIX3 src,
                                NONZERO_MATRIX3 *dest, int n);
static RC count_correlation3s(NONZERO_MATRIX3 mat, int row_num, int corre[],
                              int n, int corre_index[], int *index_size);
static RC inv_33(double mat[][3], double inv[][3], int *neg_count);
static void fill_zero33(double mat[][3]);
static RC sub_mul_333(double A[][3], double B[][3],
                      double C[][3], double D[][3]);
static RC dia_inv3s(NONZERO_MATRIX3 matrix, long double d[]);
static RC scale_solve3s(int matrix_size, const long double d[],
                        const long double r[], long double p[]);
static RC bic_decomp3s_mod(double c, NONZERO_MATRIX3 matrix,
                           NONZERO_ELEM3 d[], NONZERO_ELEM3 d_inv[]);
static RC ic_decomp3s(double c, NONZERO_MATRIX3 matrix, long double *d);
static RC ic_subst3s(NONZERO_MATRIX3 matrix, const long double *d,
                     const long double *r, long double *p);
static RC bic_decomp3s(double c, NONZERO_MATRIX3 matrix,
                       NONZERO_ELEM3 d[], NONZERO_ELEM3 d_inv[]);
static RC bic_subst3s(NONZERO_MATRIX3 matrix, const NONZERO_ELEM3 d[],
                      const NONZERO_ELEM3 d_inv[], const long double r[],
                      long double *p);
static RC mul_matrix_vect3sl(NONZERO_MATRIX3 matrix,
                             const long double *vect, long double *ans);
static void mul_nonzero_matrix_ATAb(NONZERO_MATRIX matrix, const double *vector,
                                    double *ans);
static int search_col3(NONZERO_MATRIX3 matrix, int index, int col);
static void mul_matrix_33(double A[][3], double B[][3]);
static void mul_matrix_t33(double A[][3], double B[][3]);


RC modal_const_nonzero3s(NONZERO_MATRIX3 *matrix,
                         int mode_size, const double mode[])
{
	int ii1;
	int old_size = (matrix->size);
	int dia_index;
	int elem_size;
	int elem_index;
	
	if((mode_size/3) > (matrix->size)) return(ARG_ERROR_RC);

	elem_size = 0;
	for(ii1=0; ii1<(mode_size/3); ii1++){
		if( (nearly_eq(mode[3*ii1 + 0], 0.0))
		  &&(nearly_eq(mode[3*ii1 + 1], 0.0))
		  &&(nearly_eq(mode[3*ii1 + 2], 0.0)) ) continue;
		elem_size++;
	}
	elem_size++;

	(matrix->size)++;
	(matrix->row) = realloc((matrix->row), (matrix->size)*sizeof(NONZERO_ROW3));
	RC_NULL_CHK(matrix->row);
	(matrix->row)[old_size].size = elem_size;
	(matrix->row)[old_size].elem
	       = malloc( ((matrix->row)[old_size].size) * sizeof(NONZERO_ELEM3) );
	RC_NULL_CHK((matrix->row)[old_size].elem);
	(matrix->row)[old_size].syn_size = 0;
	(matrix->row)[old_size].syn = NULL;

	dia_index = (matrix->row)[old_size].size - 1;
	elem_index = 0;
	for(ii1=0; ii1<(mode_size/3); ii1++){
		if(elem_index >= elem_size) return(UNKNOWN_ERROR_RC);
		if( (nearly_eq(mode[3*ii1 + 0], 0.0))
		  &&(nearly_eq(mode[3*ii1 + 1], 0.0))
		  &&(nearly_eq(mode[3*ii1 + 2], 0.0)) ){
			continue;
		}

		(matrix->row)[old_size].elem[elem_index].col = ii1;
		fill_zero33((matrix->row)[old_size].elem[elem_index].v);
		(matrix->row)[old_size].elem[elem_index].v[0][0] = mode[ii1*3 + 0];
		(matrix->row)[old_size].elem[elem_index].v[0][1] = mode[ii1*3 + 1];
		(matrix->row)[old_size].elem[elem_index].v[0][2] = mode[ii1*3 + 2];

		(matrix->row)[ii1].syn
		     = realloc( (matrix->row)[ii1].syn,
		                ((matrix->row)[ii1].syn_size + 1)*sizeof(int) );
		RC_NULL_CHK( (matrix->row)[ii1].syn );
		(matrix->row)[ii1].syn[(matrix->row)[ii1].syn_size] = old_size;
		((matrix->row)[ii1].syn_size)++;
		elem_index++;
	}
	(matrix->row)[old_size].elem[dia_index].col = old_size;
	fill_zero33((matrix->row)[old_size].elem[dia_index].v);
	(matrix->row)[old_size].elem[dia_index].v[0][0] = 2.0*ABS_TOL;
	(matrix->row)[old_size].elem[dia_index].v[1][1] = 2.0*ABS_TOL;
	(matrix->row)[old_size].elem[dia_index].v[2][2] = 2.0*ABS_TOL;

	return(NORMAL_RC);
}


/* matrix $B$K;HMQ$7$F$$$k%a%b%j!<$N%P%$%H?t$r(B size $B$KBeF~(B */
size_t count_size_nonzero_matrix3(NONZERO_MATRIX3 matrix)
{
	int ii1;
	size_t size;

	size = sizeof(NONZERO_MATRIX3) + matrix.size * sizeof(NONZERO_ROW3);

	for(ii1=0; ii1<matrix.size; ii1++){
		size += matrix.row[ii1].size * sizeof(NONZERO_ELEM3);
		size += matrix.row[ii1].syn_size * sizeof(int);
	}

	return(size);
}


static void fill_zero33(double mat[][3])
{
	mat[0][0] = 0.0;
	mat[0][1] = 0.0;
	mat[0][2] = 0.0;
	mat[1][0] = 0.0;
	mat[1][1] = 0.0;
	mat[1][2] = 0.0;
	mat[2][0] = 0.0;
	mat[2][1] = 0.0;
	mat[2][2] = 0.0;
}


/* $BIT40A4%3%l%9%-!<J,2rMQ$K(B src $B$r(B dest $B$K%3%T!<(B */
/* $B$?$@$7!"(Bsyn[] $B$O%3%T!<$7$J$$(B */
/* n: $B3F9T$K$D$-5vMF$9$k(B fill-in $B$N8D?t(B */
/* n = 0 : Level 0 */
/* n > 0 : Level 1 */
static RC clone_nonzero_matrix3(NONZERO_MATRIX3 src,
                                NONZERO_MATRIX3 *dest, int n)
{
	int ii1, ii2, ii3;
	int *corre = NULL;
	int *corre_index = NULL;
	int index_size = 0;

	if(n > 0){
		RC_NULL_CHK( corre = (int *)malloc((src.size + n) * sizeof(int)) );
		corre_index = &(corre[src.size]);
	}

	*dest = src;
	dest->row = (NONZERO_ROW3 *)malloc(dest->size * sizeof(NONZERO_ROW3));
	RC_NULL_CHK( dest->row );

	for(ii1=0; ii1<dest->size; ii1++){
		dest->row[ii1] = src.row[ii1];

		if(n > 0){
			RC_TRY( count_correlation3s(src, ii1, corre, n,
			                            corre_index, &index_size) );
		}else{
			index_size = 0;
		}
		dest->row[ii1].size += index_size;

		if(dest->row[ii1].size > 0){
			dest->row[ii1].elem = (NONZERO_ELEM3 *)malloc(dest->row[ii1].size
			                                          * sizeof(NONZERO_ELEM3));
			RC_NULL_CHK( dest->row[ii1].elem );
		}else{
			dest->row[ii1].elem = NULL;
		}
		for(ii2=0; ii2<(src.row[ii1].size); ii2++){
			dest->row[ii1].elem[ii2] = src.row[ii1].elem[ii2];
		}
		for(ii2=0; ii2<index_size; ii2++){
			int elem_index = (src.row[ii1].size) + ii2;

			dest->row[ii1].elem[elem_index].col = corre_index[ii2];
			fill_zero33(dest->row[ii1].elem[elem_index].v);
			/* $B%P%V%k%=!<%H(B */
			for(ii3=elem_index; ii3>0; ii3--){
				if( dest->row[ii1].elem[ii3].col
				  < dest->row[ii1].elem[ii3-1].col ){
					NONZERO_ELEM3 tmp = dest->row[ii1].elem[ii3-1];
					dest->row[ii1].elem[ii3-1] = dest->row[ii1].elem[ii3];
					dest->row[ii1].elem[ii3] = tmp;
				}else{
					break;
				}
			}
		}

		dest->row[ii1].syn_size = 0;
		dest->row[ii1].syn = NULL;
	}

	if(n > 0){
		free(corre);
	}

	return(NORMAL_RC);
}


/* $B%<%m%V%m%C%/$NAj4X$r5a$a(B corre[] $B$KBeF~(B */
/* corre[]$BFb$GBg$-$JMWAG$r:GBg(B n $B8DA*$s$G!"$=$N%$%s%G%C%/%9$r(B */
/* corre_index[] $B$KBeF~!"(Bindex_size $B$O<B:]$KBeF~$7$?8D?t(B */
static RC count_correlation3s(NONZERO_MATRIX3 mat, int row_num, int corre[],
                              int n, int corre_index[], int *index_size)
{
	int ii1, ii2;
	int index0, max_index0;
	int index1, max_index1;
	int elem_index;
	
	for(ii1=mat.row[row_num].elem[0].col+1; ii1<row_num; ii1++){
		corre[ii1] = 0;
	}

	elem_index = 1;
	for(ii1=mat.row[row_num].elem[0].col+1; ii1<row_num; ii1++){
		if(ii1 == mat.row[row_num].elem[elem_index].col){
			/* $BHs%<%m%V%m%C%/$O%Q%9(B */
			elem_index++;
			continue;
		}

		index0 = 0;
		max_index0 = elem_index;
		index1 = 0;
		max_index1 = mat.row[ii1].size - 1;

		while((index0 <= elem_index)&&(index1 < max_index1)){
			int col0 = mat.row[row_num].elem[index0].col;
			int col1 = mat.row[ii1].elem[index1].col;

			if(col0 < col1){
				index0++;
			}else if(col0 > col1){
				index1++;
			}else{
				corre[ii1]++;
				index0++;
				index1++;
			}
		}
	}

	for(ii1=0; ii1<n; ii1++){
		int max_corre = 0;
		int max_corre_index = 0;

		for(ii2=mat.row[row_num].elem[0].col+1; ii2<row_num; ii2++){
			if(corre[ii2] > max_corre){
				max_corre = corre[ii2];
				max_corre_index = ii2;
			}
		}
		if(max_corre == 0) break;
		corre_index[ii1] = max_corre_index;
		corre[max_corre_index] = 0;
	}
	*index_size = ii1;

	return(NORMAL_RC);
}


/* 3 * 3 $B$N5U9TNs7W;;(B */
static RC inv_33(double mat[][3], double inv[][3], int *neg_count)
{
	int ii1, ii2, ii3;
	double tmp[3][3];
	double div, mul;

	for(ii1=0; ii1<3; ii1++){
		tmp[ii1][0] = mat[ii1][0];
		tmp[ii1][1] = mat[ii1][1];
		tmp[ii1][2] = mat[ii1][2];
	}

	for(ii1=0; ii1<3; ii1++){
		inv[ii1][0] = inv[ii1][1] = inv[ii1][2] = 0.0;
		inv[ii1][ii1] = 1.0;
	}

	for(ii1=0; ii1<3; ii1++){
		if( fabs(tmp[ii1][ii1])
		    < ((1.0e+5)*(DBL_EPSILON*fabs(mat[ii1][ii1]))+ABS_TOL) ){
			return(CAL_ERROR_RC);
		}
		if( tmp[ii1][ii1] < 0.0 ){
			tmp[ii1][ii1] *= -10.0;   /* $BHs@5DjCMBP:v(B */
			(*neg_count)++;
		}

		div = 1.0/tmp[ii1][ii1];
		for(ii2=ii1+1; ii2<3; ii2++) tmp[ii1][ii2] *= div;
		for(ii2=0; ii2<3; ii2++)     inv[ii1][ii2] *= div;

		for(ii2=0; ii2<3; ii2++){
			if(ii2 == ii1) continue;

			mul = tmp[ii2][ii1];
			for(ii3=ii1+1; ii3<3; ii3++) tmp[ii2][ii3] -= mul*tmp[ii1][ii3];
			for(ii3=0; ii3<3; ii3++)     inv[ii2][ii3] -= mul*inv[ii1][ii3];
		}
	}

	return(NORMAL_RC);
}


/* D - A*B*C^t => D */
static RC sub_mul_333(double A[][3], double B[][3],
                      double C[][3], double D[][3])
{
	int ii1;
	long double tmp[3][3];

	/* A*B => tmp */
	for(ii1=0; ii1<3; ii1++){
		tmp[ii1][0] = (long double)(A[ii1][0]) * B[0][0]
		            + (long double)(A[ii1][1]) * B[1][0]
		            + (long double)(A[ii1][2]) * B[2][0];
		tmp[ii1][1] = (long double)(A[ii1][0]) * B[0][1]
		            + (long double)(A[ii1][1]) * B[1][1]
		            + (long double)(A[ii1][2]) * B[2][1];
		tmp[ii1][2] = (long double)(A[ii1][0]) * B[0][2]
		            + (long double)(A[ii1][1]) * B[1][2]
		            + (long double)(A[ii1][2]) * B[2][2];
	}

	/* tmp*C^t => D */
	for(ii1=0; ii1<3; ii1++){
		D[ii1][0] -= (double)(tmp[ii1][0] * C[0][0]
		                    + tmp[ii1][1] * C[0][1]
		                    + tmp[ii1][2] * C[0][2]);
		D[ii1][1] -= (double)(tmp[ii1][0] * C[1][0]
		                    + tmp[ii1][1] * C[1][1]
		                    + tmp[ii1][2] * C[1][2]);
		D[ii1][2] -= (double)(tmp[ii1][0] * C[2][0]
		                   + tmp[ii1][1] * C[2][1]
		                   + tmp[ii1][2] * C[2][2]);
	}

	return(NORMAL_RC);
}


/* X[X_size][3*B.size]^T A X $B$r7W;;$7!"(BXtAX[X_size][X_size] $B$KBeF~(B */
RC mul_nonzero_matrix_XtAX3s(NONZERO_MATRIX3 A, int X_size, double **X,
                             double **XtAX)
{
	int ii1, ii2;
	int size = 3*(A.size);
	double *Ax;

	RC_TRY( allocate1D(size, &Ax) );

	for(ii1=0; ii1<X_size; ii1++){
		RC_TRY( mul_nonzero_matrix_vect3s(A, X[ii1], Ax) );
		for(ii2=0; ii2<X_size; ii2++){
			XtAX[ii1][ii2] = inner_product(size, X[ii2], Ax);
		}
	}
	for(ii1=0; ii1<X_size; ii1++){
		for(ii2=ii1+1; ii2<X_size; ii2++){
			XtAX[ii1][ii2] = XtAX[ii2][ii1]
			               = 0.5*(XtAX[ii1][ii2] + XtAX[ii2][ii1]);
		}
	}

	free(Ax);

	return(NORMAL_RC);
}


/* vects[vect_size][3*B.size] $B$r(B B-$BD>8r2=$7$F>e=q$-(B */
/* Bvects != NULL $B$J$i(B B * vects[] $B$rBeF~(B */
RC nonzero_B_orthogonalize3s(NONZERO_MATRIX3 B, int vect_size,
                             double **vects, double **Bvects)
{
	int ii1, ii2, ii3;
	int size = 3*(B.size);
	double *Bv = NULL;
	double vtBv, isq_vtBv;
	double v2tBv;

	if(Bvects == NULL){
		RC_TRY( allocate1D(size, &Bv) );
	}

	for(ii1=0; ii1<vect_size; ii1++){
		/* ii1 $BHVL\$N%Y%/%H%k$r(B B $B$K$D$$$F@55,2=(B */
		if(Bvects != NULL) Bv = Bvects[ii1];
		RC_TRY( mul_nonzero_matrix_vect3s(B, vects[ii1], Bv) );
		vtBv = inner_product(size, vects[ii1], Bv);
		if(vtBv < ABS_ERROR) return(CAL_ERROR_RC);
		isq_vtBv = 1.0/(sqrt(vtBv));
		for(ii2=0; ii2<size; ii2++){
			vects[ii1][ii2] *= isq_vtBv;
			Bv[ii2] *= isq_vtBv;
		}

		/* ii1 $BHVL\@.J,$r0z$/(B */
		for(ii2=ii1+1; ii2<vect_size; ii2++){
			v2tBv = inner_product(size, vects[ii2], Bv);
			for(ii3=0; ii3<size; ii3++){
				vects[ii2][ii3] -= v2tBv * vects[ii1][ii3];
			}
		}
	}

	if(Bvects == NULL){
		free(Bv);
	}

	return(NORMAL_RC);
}


int count_nonzero_matrix3(NONZERO_MATRIX3 matrix)
{
	int ii1;
	int ret = 0;

	for(ii1=0; ii1<matrix.size; ii1++){
		ret += matrix.row[ii1].size;
	}
	return(ret);
}


/* $B%<%m$N9`$r:o=|(B */
RC delete_zero_matrix3s(NONZERO_MATRIX3 matrix)
{
	int ii1, ii2, ii3;

	for(ii1=0; ii1<matrix.size; ii1++){
		int col_count = 0;
		int dia_elem = matrix.row[ii1].size-1;
		for(ii2=0; ii2<dia_elem; ii2++){
			double (*val)[3] = matrix.row[ii1].elem[ii2].v;
			int col = matrix.row[ii1].elem[ii2].col;

			if( (!nearly_eq(val[0][0], 0.0))||(!nearly_eq(val[0][1], 0.0))
			  ||(!nearly_eq(val[0][2], 0.0))||(!nearly_eq(val[1][0], 0.0))
			  ||(!nearly_eq(val[1][1], 0.0))||(!nearly_eq(val[1][2], 0.0))
			  ||(!nearly_eq(val[2][0], 0.0))||(!nearly_eq(val[2][1], 0.0))
			  ||(!nearly_eq(val[2][2], 0.0)) ){
				/* $BHs%<%m9`$rG[NsA0J}$K5M$a$k(B */
				matrix.row[ii1].elem[col_count] = matrix.row[ii1].elem[ii2];
				col_count++;
				continue;
			}
			for(ii3=0; ii3<(matrix.row[col].syn_size); ii3++){
				if(matrix.row[col].syn[ii3] == ii1){
					/* syn[] $B$K:o=|%U%i%0$rN)$F$k(B */
					matrix.row[col].syn[ii3] = -1;
					break;
				}
			}
			if(ii3 >= (matrix.row[col].syn_size)) return(ARG_ERROR_RC);
		}
		matrix.row[ii1].elem[col_count] = matrix.row[ii1].elem[dia_elem];
		col_count++;
		matrix.row[ii1].elem = realloc(matrix.row[ii1].elem,
		                               col_count * sizeof(NONZERO_ELEM3));
		RC_NULL_CHK( matrix.row[ii1].elem );
		matrix.row[ii1].size = col_count;
	}

	for(ii1=0; ii1<matrix.size; ii1++){
		if(matrix.row[ii1].syn_size > 0){
			int syn_count = 0;
			for(ii2=0; ii2<matrix.row[ii1].syn_size;ii2++){
				if(matrix.row[ii1].syn[ii2] >= 0){
					matrix.row[ii1].syn[syn_count] = matrix.row[ii1].syn[ii2];
					syn_count++;
				}
			}
			if(syn_count > 0){
				matrix.row[ii1].syn = realloc(matrix.row[ii1].syn,
				                              syn_count * sizeof(int));
				RC_NULL_CHK( matrix.row[ii1].syn );
			}else{
				matrix.row[ii1].syn = NULL;
			}
			matrix.row[ii1].syn_size = syn_count;
		}
	}

	return(NORMAL_RC);
}


double nonzero_matrix3s_val(NONZERO_MATRIX3 matrix, int row, int col)
{
	int index;

	if(row < col){
		int tmp = row;
		row = col;
		col = tmp;
	}

	index = search_col3(matrix, row/3, col/3);
	if(index < 0) return(0.0);

	return( matrix.row[row/3].elem[index].v[row%3][col%3] );
}


RC print_nonzero_matrix3(FILE *fp, NONZERO_MATRIX3 matrix)
{
	int ii1, ii2;

	for(ii1=0; ii1<matrix.size; ii1++){
		fprintf(fp, "%4d-%2d:%2d", ii1, matrix.row[ii1].size,
		        matrix.row[ii1].syn_size);
		for(ii2=0; ii2<matrix.row[ii1].size; ii2++){
			fprintf(fp, " %d", matrix.row[ii1].elem[ii2].col);
		}
		fprintf(fp, "\n          ");
		for(ii2=0; ii2<matrix.row[ii1].syn_size; ii2++){
			fprintf(fp, " %d", matrix.row[ii1].syn[ii2]);
		}
		fprintf(fp, "\n");
	}

	return(NORMAL_RC);
}


RC free_nonzero_matrix3(NONZERO_MATRIX3 *matrix)
{
	int ii1;

	if(matrix->size <= 0) return(ARG_ERROR_RC);
	if(matrix->row == NULL) return(ARG_ERROR_RC);

	for(ii1=0; ii1<matrix->size; ii1++){
		if(matrix->row[ii1].size > 0){
			if(matrix->row[ii1].elem == NULL) return(ARG_ERROR_RC);
			free(matrix->row[ii1].elem);
		}
		matrix->row[ii1].size = 0;
		matrix->row[ii1].elem = NULL;

		if(matrix->row[ii1].syn_size > 0){
			if(matrix->row[ii1].syn == NULL) return(ARG_ERROR_RC);
			free(matrix->row[ii1].syn);
		}
		matrix->row[ii1].syn_size = 0;
		matrix->row[ii1].syn = NULL;
	}
	free(matrix->row);
	matrix->row = NULL;
	matrix->size = 0;

	return(NORMAL_RC);
}


NONZERO_ELEM3 *nonzero_matrix3s_ptr(NONZERO_MATRIX3 matrix, int row, int col)
{
	int ii1;
	NONZERO_ELEM3 *tmp_ptr;
	NONZERO_ELEM3 *ret;
	int *itmp_ptr;
	int old_size;

	if(row < col) return(NULL);

	if((col < 0)||(row >= matrix.size)){
		fprintf(stderr, "col | row\n");
		return(NULL);
	}

	/* row $B9T(B col $BNs$r8!:w(B */
	for(ii1=0; ii1<matrix.row[row].size; ii1++){
		if(col == matrix.row[row].elem[ii1].col){
			return(&(matrix.row[row].elem[ii1]));
		}
	}

	/* $B8+$D$+$i$J$+$C$?$N$GMWAG$rDI2C(B */
	old_size = matrix.row[row].size;
	tmp_ptr = realloc(matrix.row[row].elem,
	                  (old_size + 1)*sizeof(NONZERO_ELEM3));
	if(tmp_ptr == NULL) return(NULL);
	matrix.row[row].elem = tmp_ptr;
	matrix.row[row].elem[old_size].col = col;
	fill_zero33(matrix.row[row].elem[old_size].v);
	(matrix.row[row].size)++;
	ret = &(matrix.row[row].elem[old_size]);
	for(ii1=old_size; ii1>=1; ii1--){
		if(matrix.row[row].elem[ii1].col < matrix.row[row].elem[ii1-1].col){
			NONZERO_ELEM3 tmp_elem = matrix.row[row].elem[ii1-1];
			matrix.row[row].elem[ii1-1] = matrix.row[row].elem[ii1];
			matrix.row[row].elem[ii1] = tmp_elem;
			ret = &(matrix.row[row].elem[ii1-1]);
		}else{
			break;
		}
	}

	/* $BBP>NMWAG$N(B syn[] $B$K$bDI2C(B */
	old_size = matrix.row[col].syn_size;
	itmp_ptr = realloc(matrix.row[col].syn, (old_size + 1)*sizeof(int));
	if(itmp_ptr == NULL) return(NULL);
	matrix.row[col].syn = itmp_ptr;
	matrix.row[col].syn[old_size] = row;
	(matrix.row[col].syn_size)++;
	for(ii1=old_size; ii1>=1; ii1--){
		if(matrix.row[col].syn[ii1] < matrix.row[col].syn[ii1-1]){
			int tmp_syn = matrix.row[col].syn[ii1-1];
			matrix.row[col].syn[ii1-1] = matrix.row[col].syn[ii1];
			matrix.row[col].syn[ii1] = tmp_syn;
		}else{
			break;
		}
	}

	return(ret);
}


/* precondition < 0 : $BA0=hM}$J$7(B */
/*              = 0 : $BBP3Q%9%1!<%j%s%0(B */
/*              = 1 : $B%]%$%s%HIT40A4J,2r!JBP3Q9`$N$_!K(B */
/*              = 2 : $B%V%m%C%/IT40A4J,2r!JBP3Q%V%m%C%/$N$_!K(B */
/*              = 3 : $B%V%m%C%/IT40A4J,2r!J(Bfill lebel 0$B!K(B */
/*              > 3 : $B%V%m%C%/IT40A4J,2r!J(Bfill lebel 1$B!K(B */
RC iccg_nonzero3s(NONZERO_MATRIX3 matrix, const double *b, double *x,
                                       int precondition, int verb_sw)
{
	long double *d = NULL;
	NONZERO_ELEM3 *block_d = NULL;
	NONZERO_ELEM3 *block_d_inv = NULL;
	NONZERO_MATRIX3 matrix_decomp;
	long double *p;
	long double *q;
	long double *r;
	long double *lx;
	long double c1, c2, c3;
	long double alpha, beta;
	long double r_r;
	long double b_b;
	double c;
	int ii1, ii2;
	double rel_error_b;
	int it_max = 2*3*matrix.size + 100;
	int update = 4*(int)(sqrt((double)(matrix.size))) + 10;
	int error_flag = 0;
	int conv_count = 0;

	if(precondition >= 2){
		RC_NULL_CHK( block_d = (NONZERO_ELEM3 *)malloc(sizeof(NONZERO_ELEM3)
		                                               * matrix.size) );
		RC_NULL_CHK( block_d_inv = (NONZERO_ELEM3 *)malloc(sizeof(NONZERO_ELEM3)
		                                                   * matrix.size) );
	}else if(precondition >= 0){
		RC_NULL_CHK( d =  (long double *)malloc(3 * sizeof(long double)
		                                          * matrix.size) );
	}
	RC_NULL_CHK( p =  (long double *)malloc(3 * sizeof(long double)
	                                          * matrix.size) );
	RC_NULL_CHK( q =  (long double *)malloc(3 * sizeof(long double)
	                                          * matrix.size) );
	RC_NULL_CHK( r =  (long double *)malloc(3 * sizeof(long double)
	                                          * matrix.size) );
	RC_NULL_CHK( lx = (long double *)malloc(3 * sizeof(long double)
	                                          * matrix.size) );

	c = 1.0;
	for(ii1=0; ii1<15; ii1++){
		RC rc = NORMAL_RC;

		if(precondition == 0){
			rc = dia_inv3s(matrix, d);
		}else if(precondition == 1){
			rc = ic_decomp3s(c, matrix, d);
		}else if(precondition == 2){
			rc = bic_decomp3s(c, matrix, block_d, block_d_inv);
		}else if(precondition >= 3){
			RC_TRY( clone_nonzero_matrix3(matrix, &matrix_decomp,
			                                      precondition-3) );
			rc = bic_decomp3s_mod(c, matrix_decomp, block_d, block_d_inv);
		}
		if(rc == NORMAL_RC) break;
		if(rc == CAL_ERROR_RC){
			if(precondition >= 3){
				RC_TRY( free_nonzero_matrix3(&matrix_decomp) );
			}
			c *= 1.2;
			continue;
		}
		fprintf(stderr, "ic_decomp3s() < \n");
		return(rc);
	}
	if(ii1 >= 15){
		fprintf(stderr, "Preconditioning failed!!\n");
		return(CAL_ERROR_RC);
	}
	/* if(verb_sw) fprintf(stderr, "c = %15.7e\n", c); */

	for(ii1=0; ii1<matrix.size; ii1++){
		lx[ii1*3]   = x[ii1*3];
		lx[ii1*3+1] = x[ii1*3+1];
		lx[ii1*3+2] = x[ii1*3+2];
	}
	RC_TRY( mul_matrix_vect3sl(matrix, lx, q) );
	b_b = ABS_ERROR;
	for(ii1=0; ii1<matrix.size; ii1++){
		r[ii1*3]   = b[ii1*3]   - q[ii1*3];
		r[ii1*3+1] = b[ii1*3+1] - q[ii1*3+1];
		r[ii1*3+2] = b[ii1*3+2] - q[ii1*3+2];
		b_b += b[ii1*3] * b[ii1*3]
		     + b[ii1*3+1] * b[ii1*3+1]
		     + b[ii1*3+2] * b[ii1*3+2];
	}
	if(precondition < 0){
		RC_TRY( scale_solve3s(matrix.size, NULL, r, p) );
	}else if(precondition == 0){
		RC_TRY( scale_solve3s(matrix.size, d, r, p) );
	}else if(precondition == 1){
		RC_TRY( ic_subst3s(matrix, d, r, p) );
	}else if(precondition == 2){
		RC_TRY( bic_subst3s(matrix, block_d, block_d_inv, r, p) );
	}else{ /* precondition >= 3 */
		RC_TRY( bic_subst3s(matrix_decomp, block_d, block_d_inv, r, p) );
	}

	c1 = ABS_ERROR;
	for(ii1=0; ii1<matrix.size; ii1++){
		c1 += (r[ii1*3]   * p[ii1*3]
		     + r[ii1*3+1] * p[ii1*3+1]
		     + r[ii1*3+2] * p[ii1*3+2]);
	}

	for(ii1=0; ii1<it_max; ii1++){
		RC_TRY( mul_matrix_vect3sl(matrix, p, q) );
		c2 = ABS_ERROR;
		for(ii2=0; ii2<matrix.size; ii2++){
			c2 += (p[ii2*3]   * q[ii2*3]
			     + p[ii2*3+1] * q[ii2*3+1]
			     + p[ii2*3+2] * q[ii2*3+2]);
		}
		alpha = c1/c2;
		for(ii2=0; ii2<matrix.size; ii2++){
			lx[ii2*3]   += alpha * p[ii2*3];
			lx[ii2*3+1] += alpha * p[ii2*3+1];
			lx[ii2*3+2] += alpha * p[ii2*3+2];
		}

		r_r = 0.0;
		if(ii1%update == 0){
			RC_TRY( mul_matrix_vect3sl(matrix, lx, q) );
			for(ii2=0; ii2<matrix.size; ii2++){
				r[ii2*3]   = b[ii2*3]   - q[ii2*3];
				r[ii2*3+1] = b[ii2*3+1] - q[ii2*3+1];
				r[ii2*3+2] = b[ii2*3+2] - q[ii2*3+2];
				r_r += (r[ii2*3]   * r[ii2*3]
				      + r[ii2*3+1] * r[ii2*3+1]
				      + r[ii2*3+2] * r[ii2*3+2]);
			}
		}else{
			for(ii2=0; ii2<matrix.size; ii2++){
				r[ii2*3]   -= alpha * q[ii2*3];
				r[ii2*3+1] -= alpha * q[ii2*3+1];
				r[ii2*3+2] -= alpha * q[ii2*3+2];
				r_r += (r[ii2*3]   * r[ii2*3]
				      + r[ii2*3+1] * r[ii2*3+1]
				      + r[ii2*3+2] * r[ii2*3+2]);
			}
		}

		rel_error_b = sqrt( (r_r/b_b)/(3*matrix.size) );
		if( verb_sw && ((ii1%5 == 0)||(matrix.size > 10000)) ){
			fprintf(stderr, "%6d:%10.2e\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",
		                                               ii1, rel_error_b);
		}

		/* $B8m:9$,(B REL_ERROR1 $B0J2<$H$J$k$+(B 10 $B2sO"B3$G(B REL_ERROR2 $B0J2<(B */
		/* $B$H$J$C$?;~E@$G=*N;(B                                         */
		if((rel_error_b <= REL_ERROR1)||(conv_count >= 10)){
			if(verb_sw) fprintf(stderr, "%6d: converged          \n", ii1);
			break;
		}
		if(rel_error_b <= REL_ERROR2){
			conv_count++;
		}else{
			conv_count = 0;
		}

		if(precondition < 0){
			RC_TRY( scale_solve3s(matrix.size, NULL, r, q) );
		}else if(precondition == 0){
			RC_TRY( scale_solve3s(matrix.size, d, r, q) );
		}else if(precondition == 1){
			RC_TRY( ic_subst3s(matrix, d, r, q) );
		}else if(precondition == 2){
			RC_TRY( bic_subst3s(matrix, block_d, block_d_inv, r, q) );
		}else{ /* precondition >= 3 */
			RC_TRY( bic_subst3s(matrix_decomp, block_d, block_d_inv, r, q) );
		}
		
		c3 = ABS_ERROR;
		for(ii2=0; ii2<matrix.size; ii2++){
			c3 += (r[ii2*3]   * q[ii2*3]
			     + r[ii2*3+1] * q[ii2*3+1]
			     + r[ii2*3+2] * q[ii2*3+2]);
		}
		beta = c3/c1;
		c1 = c3;
		for(ii2=0; ii2<matrix.size; ii2++){
			p[ii2*3]   = q[ii2*3]   + beta * p[ii2*3];
			p[ii2*3+1] = q[ii2*3+1] + beta * p[ii2*3+1];
			p[ii2*3+2] = q[ii2*3+2] + beta * p[ii2*3+2];
		}
	}
	if(ii1 >= it_max) error_flag = 1;

	for(ii1=0; ii1<matrix.size; ii1++){
		x[ii1*3]   = lx[ii1*3];
		x[ii1*3+1] = lx[ii1*3+1];
		x[ii1*3+2] = lx[ii1*3+2];
	}

	if(precondition >= 2){
		free(block_d);
		free(block_d_inv);
		if(precondition >= 3){
			RC_TRY( free_nonzero_matrix3(&matrix_decomp) );
		}
	}else if(precondition >= 0){
		free(d);
	}
	free(p);
	free(q);
	free(r);
	free(lx);

	if(error_flag){
		fprintf(stderr, "WARNING: CG iteration was NOT converged!!\n");
		/* return(CAL_ERROR_RC); */
	}

	return(NORMAL_RC);
}


/* matrix $B$NBP3Q9`$N5U?t$r(B d[] $B$KBeF~(B */
static RC dia_inv3s(NONZERO_MATRIX3 matrix, long double d[])
{
	int ii1;

	for(ii1=0; ii1<matrix.size; ii1++){
		int dia_pos = matrix.row[ii1].size - 1;

		if(matrix.row[ii1].elem[dia_pos].col != ii1) return(SEEK_ERROR_RC);
/*		if( (matrix.row[ii1].elem[dia_pos].v[0][0] < 0.1*ABS_TOL)
		  ||(matrix.row[ii1].elem[dia_pos].v[1][1] < 0.1*ABS_TOL)
		  ||(matrix.row[ii1].elem[dia_pos].v[2][2] < 0.1*ABS_TOL) ){
			return(CAL_ERROR_RC);
		}*/
		if( (fabs(matrix.row[ii1].elem[dia_pos].v[0][0]) < 10.0*ABS_TOL)
		  ||(fabs(matrix.row[ii1].elem[dia_pos].v[1][1]) < 10.0*ABS_TOL)
		  ||(fabs(matrix.row[ii1].elem[dia_pos].v[2][2]) < 10.0*ABS_TOL) ){
			d[ii1*3 + 0] = 1.0;
			d[ii1*3 + 1] = 1.0;
			d[ii1*3 + 2] = 1.0;
			continue;
		}
		d[ii1*3 + 0] = 1.0/matrix.row[ii1].elem[dia_pos].v[0][0];
		d[ii1*3 + 1] = 1.0/matrix.row[ii1].elem[dia_pos].v[1][1];
		d[ii1*3 + 2] = 1.0/matrix.row[ii1].elem[dia_pos].v[2][2];
	}

	return(NORMAL_RC);
}


static RC scale_solve3s(int matrix_size, const long double d[],
                        const long double r[], long double p[])
{
	int ii1;

	if(d == NULL){
		/* $BC1$J$k%3%T!<(B */
		for(ii1=0; ii1<matrix_size; ii1++){
			p[ii1*3 + 0] = r[ii1*3 + 0];
			p[ii1*3 + 1] = r[ii1*3 + 1];
			p[ii1*3 + 2] = r[ii1*3 + 2];
		}
	}else{
		/* d[] $B$r3]$1$J$,$i%3%T!<(B */
		for(ii1=0; ii1<matrix_size; ii1++){
			p[ii1*3 + 0] = d[ii1*3 + 0] * r[ii1*3 + 0];
			p[ii1*3 + 1] = d[ii1*3 + 1] * r[ii1*3 + 1];
			p[ii1*3 + 2] = d[ii1*3 + 2] * r[ii1*3 + 2];
		}
	}

	return(NORMAL_RC);
}



/* $B%V%m%C%/IT40A4J,2r(B */
/* matrix $B$r=q$-49$((B(Level 0) */
/* d[]$B$KJ,2r8e$NBP3Q%V%m%C%/!"(Bd_inv[] $B$KBP3Q%V%m%C%/$N5U9TNs$rBeF~(B */
static RC bic_decomp3s_mod(double c, NONZERO_MATRIX3 matrix,
                           NONZERO_ELEM3 d[], NONZERO_ELEM3 d_inv[])
{
	int ii1, ii2, ii3;
	int neg_count = 0;
	int max_neg_count = matrix.size/100000 + 3;

	for(ii1=0; ii1<matrix.size; ii1++){
		int dia_pos = matrix.row[ii1].size - 1;
		RC rc;

		if(matrix.row[ii1].elem[dia_pos].col != ii1) return(SEEK_ERROR_RC);
		for(ii2=1; ii2<dia_pos; ii2++){
			int index;
			int col_ii2 = matrix.row[ii1].elem[ii2].col;
			int max_index = matrix.row[col_ii2].size - 1;
			if(col_ii2 >= ii1) return(SEEK_ERROR_RC);

			for(ii3=0, index=0; (ii3<ii2)&&(index<max_index); ){
				if( matrix.row[col_ii2].elem[index].col
				  < matrix.row[ii1].elem[ii3].col ){
					index++;
				}else if( matrix.row[col_ii2].elem[index].col
				        > matrix.row[ii1].elem[ii3].col ){
					ii3++;
				}else{
					int col_ii3 = matrix.row[ii1].elem[ii3].col;
					RC_TRY( sub_mul_333(matrix.row[ii1].elem[ii3].v,
					                    d_inv[col_ii3].v,
					                    matrix.row[col_ii2].elem[index].v,
					                    matrix.row[ii1].elem[ii2].v) );
					index++;
					ii3++;
				}
			}
		}
					
		d[ii1] = matrix.row[ii1].elem[dia_pos];
		d[ii1].v[0][0] *= c;
		d[ii1].v[0][1] *= c;
		d[ii1].v[0][2] *= c;
		d[ii1].v[1][0] *= c;
		d[ii1].v[1][1] *= c;
		d[ii1].v[1][2] *= c;
		d[ii1].v[2][0] *= c;
		d[ii1].v[2][1] *= c;
		d[ii1].v[2][2] *= c;

		for(ii2=0; ii2<dia_pos; ii2++){
			int col = matrix.row[ii1].elem[ii2].col;
			if(col >= ii1) return(SEEK_ERROR_RC);

			RC_TRY( sub_mul_333(matrix.row[ii1].elem[ii2].v, d_inv[col].v,
			                    matrix.row[ii1].elem[ii2].v, d[ii1].v) );
		}
		rc = inv_33(d[ii1].v, d_inv[ii1].v, &neg_count);
		if( (neg_count > max_neg_count)||(rc == CAL_ERROR_RC) ){
			return(CAL_ERROR_RC);
		}
		RC_TRY(rc);
	}

	return(NORMAL_RC);
}


/* $B%V%m%C%/IT40A4J,2r(B */
/* d[]$B$KJ,2r8e$NBP3Q%V%m%C%/!"(Bd_inv[] $B$KBP3Q%V%m%C%/$N5U9TNs$rBeF~(B */
static RC bic_decomp3s(double c, NONZERO_MATRIX3 matrix,
                       NONZERO_ELEM3 d[], NONZERO_ELEM3 d_inv[])
{
	int ii1, ii2;
	int neg_count = 0;
	int max_neg_count = matrix.size/100000 + 3;

	for(ii1=0; ii1<matrix.size; ii1++){
		int dia_pos = matrix.row[ii1].size - 1;
		RC rc;

		if(matrix.row[ii1].elem[dia_pos].col != ii1) return(SEEK_ERROR_RC);
		d[ii1] = matrix.row[ii1].elem[dia_pos];
		d[ii1].v[0][0] *= c;
		d[ii1].v[0][1] *= c;
		d[ii1].v[0][2] *= c;
		d[ii1].v[1][0] *= c;
		d[ii1].v[1][1] *= c;
		d[ii1].v[1][2] *= c;
		d[ii1].v[2][0] *= c;
		d[ii1].v[2][1] *= c;
		d[ii1].v[2][2] *= c;

		for(ii2=0; ii2<dia_pos; ii2++){
			int col = matrix.row[ii1].elem[ii2].col;
			if(col >= ii1) return(SEEK_ERROR_RC);

			RC_TRY( sub_mul_333(matrix.row[ii1].elem[ii2].v, d_inv[col].v,
			                    matrix.row[ii1].elem[ii2].v, d[ii1].v) );
		}
		rc = inv_33(d[ii1].v, d_inv[ii1].v, &neg_count);
		if( (neg_count > max_neg_count)||(rc == CAL_ERROR_RC) ){
			return(CAL_ERROR_RC);
		}
		RC_TRY(rc);
	}

	return(NORMAL_RC);
}


static RC bic_subst3s(NONZERO_MATRIX3 matrix, const NONZERO_ELEM3 d[],
                      const NONZERO_ELEM3 d_inv[], const long double r[],
                      long double *p)
{
	int ii1, ii2;

	/* $BA0?JBeF~(B */
	for(ii1=0; ii1<matrix.size; ii1++){
		int dia_pos = matrix.row[ii1].size - 1;
		double (*val)[3];
		long double p_ii1_3_0 = r[ii1*3];
		long double p_ii1_3_1 = r[ii1*3+1];
		long double p_ii1_3_2 = r[ii1*3+2];

		for(ii2=0; ii2<dia_pos; ii2++){
			int col = matrix.row[ii1].elem[ii2].col;

			val = matrix.row[ii1].elem[ii2].v;
			p_ii1_3_0 -= (p[col*3]   * val[0][0]
			            + p[col*3+1] * val[0][1]
			            + p[col*3+2] * val[0][2]);
			p_ii1_3_1 -= (p[col*3]   * val[1][0]
			            + p[col*3+1] * val[1][1]
			            + p[col*3+2] * val[1][2]);
			p_ii1_3_2 -= (p[col*3]   * val[2][0]
			            + p[col*3+1] * val[2][1]
			            + p[col*3+2] * val[2][2]);
		}
		p[ii1*3+0] = d_inv[ii1].v[0][0] * p_ii1_3_0
		           + d_inv[ii1].v[0][1] * p_ii1_3_1
		           + d_inv[ii1].v[0][2] * p_ii1_3_2;
		p[ii1*3+1] = d_inv[ii1].v[1][0] * p_ii1_3_0
		           + d_inv[ii1].v[1][1] * p_ii1_3_1
		           + d_inv[ii1].v[1][2] * p_ii1_3_2;
		p[ii1*3+2] = d_inv[ii1].v[2][0] * p_ii1_3_0
		           + d_inv[ii1].v[2][1] * p_ii1_3_1
		           + d_inv[ii1].v[2][2] * p_ii1_3_2;
	}

	/* $BBP3Q(B */
	for(ii1=0; ii1<matrix.size; ii1++){
		long double p_ii1_3_0 = p[ii1*3];
		long double p_ii1_3_1 = p[ii1*3+1];
		long double p_ii1_3_2 = p[ii1*3+2];
		p[ii1*3+0] = d[ii1].v[0][0] * p_ii1_3_0
		           + d[ii1].v[0][1] * p_ii1_3_1
		           + d[ii1].v[0][2] * p_ii1_3_2;
		p[ii1*3+1] = d[ii1].v[1][0] * p_ii1_3_0
		           + d[ii1].v[1][1] * p_ii1_3_1
		           + d[ii1].v[1][2] * p_ii1_3_2;
		p[ii1*3+2] = d[ii1].v[2][0] * p_ii1_3_0
		           + d[ii1].v[2][1] * p_ii1_3_1
		           + d[ii1].v[2][2] * p_ii1_3_2;
	}

	/* $B8eB`BeF~(B */
	for(ii1=matrix.size-1; ii1>=0; ii1--){
		int dia_pos = matrix.row[ii1].size - 1;
		long double p_ii1_3_0 = p[ii1*3+0];
		long double p_ii1_3_1 = p[ii1*3+1];
		long double p_ii1_3_2 = p[ii1*3+2];

		p[ii1*3+0] = d_inv[ii1].v[0][0] * p_ii1_3_0
		           + d_inv[ii1].v[0][1] * p_ii1_3_1
		           + d_inv[ii1].v[0][2] * p_ii1_3_2;
		p[ii1*3+1] = d_inv[ii1].v[1][0] * p_ii1_3_0
		           + d_inv[ii1].v[1][1] * p_ii1_3_1
		           + d_inv[ii1].v[1][2] * p_ii1_3_2;
		p[ii1*3+2] = d_inv[ii1].v[2][0] * p_ii1_3_0
		           + d_inv[ii1].v[2][1] * p_ii1_3_1
		           + d_inv[ii1].v[2][2] * p_ii1_3_2;

		p_ii1_3_0 = p[ii1*3+0];
		p_ii1_3_1 = p[ii1*3+1];
		p_ii1_3_2 = p[ii1*3+2];
		for(ii2=0; ii2<dia_pos; ii2++){
			int col = matrix.row[ii1].elem[ii2].col;
			double (*val)[3] = matrix.row[ii1].elem[ii2].v;

			p[col*3+0] -= (p_ii1_3_0 * val[0][0]
			             + p_ii1_3_1 * val[1][0]
			             + p_ii1_3_2 * val[2][0]);
			p[col*3+1] -= (p_ii1_3_0 * val[0][1]
			             + p_ii1_3_1 * val[1][1]
			             + p_ii1_3_2 * val[2][1]);
			p[col*3+2] -= (p_ii1_3_0 * val[0][2]
			             + p_ii1_3_1 * val[1][2]
			             + p_ii1_3_2 * val[2][2]);
		}
	}

	return(NORMAL_RC);
}


/* matrix $B$rIT40A4%3%l%9%-!<J,2r(B(U^tDU)$B$9$k(B */
/* $B$?$@$7!"(Bmatrix $B$OJQ99$;$:!"(Bd[] $B$K(B U $B$NBP3Q@.J,$rBeF~(B */
static RC ic_decomp3s(double c, NONZERO_MATRIX3 matrix, long double *d)
{
	int ii1, ii2;

	for(ii1=0; ii1<matrix.size; ii1++){
		int dia_pos = matrix.row[ii1].size - 1;
		int zero_flag0 = 0;
		int zero_flag1 = 0;
		int zero_flag2 = 0;
		double d_tol_0, d_tol_1, d_tol_2;

		if(matrix.row[ii1].elem[dia_pos].col != ii1) return(SEEK_ERROR_RC);

		d[ii1*3]     = c * matrix.row[ii1].elem[dia_pos].v[0][0];
		d[ii1*3 + 1] = c * matrix.row[ii1].elem[dia_pos].v[1][1];
		d[ii1*3 + 2] = c * matrix.row[ii1].elem[dia_pos].v[2][2];
		d_tol_0 = (1.0e+5)*(DBL_EPSILON)
		         *(matrix.row[ii1].elem[dia_pos].v[0][0] + ABS_TOL);
		d_tol_1 = (1.0e+5)*(DBL_EPSILON)
		         *(matrix.row[ii1].elem[dia_pos].v[1][1] + ABS_TOL);
		d_tol_2 = (1.0e+5)*(DBL_EPSILON)
		         *(matrix.row[ii1].elem[dia_pos].v[2][2] + ABS_TOL);

		if(fabs(d[ii1*3]) < 10.0*ABS_TOL){
			zero_flag0 = 1;
			/* fprintf(stderr, "matrix[%d-%d]\n", ii1, 0);
			return(CAL_ERROR_RC); */
		}
		if(fabs(d[ii1*3+1]) < 10.0*ABS_TOL){
			zero_flag1 = 1;
			/* fprintf(stderr, "matrix[%d-%d]\n", ii1, 1);
			return(CAL_ERROR_RC); */
		}
		if(fabs(d[ii1*3+2]) < 10.0*ABS_TOL){
			zero_flag2 = 1;
			/* fprintf(stderr, "matrix[%d-%d]\n", ii1, 2);
			return(CAL_ERROR_RC); */
		}

		for(ii2=0; ii2<dia_pos; ii2++){
			int col = matrix.row[ii1].elem[ii2].col;
			if(col >= ii1) return(SEEK_ERROR_RC);

			d[ii1*3]   -= (matrix.row[ii1].elem[ii2].v[0][0]
			              *matrix.row[ii1].elem[ii2].v[0][0]/d[col*3]
			             + matrix.row[ii1].elem[ii2].v[0][1]
			              *matrix.row[ii1].elem[ii2].v[0][1]/d[col*3+1]
			             + matrix.row[ii1].elem[ii2].v[0][2]
			              *matrix.row[ii1].elem[ii2].v[0][2]/d[col*3+2]);
			d[ii1*3+1] -= (matrix.row[ii1].elem[ii2].v[1][0]
			              *matrix.row[ii1].elem[ii2].v[1][0]/d[col*3]
			             + matrix.row[ii1].elem[ii2].v[1][1]
			              *matrix.row[ii1].elem[ii2].v[1][1]/d[col*3+1]
			             + matrix.row[ii1].elem[ii2].v[1][2]
			              *matrix.row[ii1].elem[ii2].v[1][2]/d[col*3+2]);
			d[ii1*3+2] -= (matrix.row[ii1].elem[ii2].v[2][0]
			              *matrix.row[ii1].elem[ii2].v[2][0]/d[col*3]
			             + matrix.row[ii1].elem[ii2].v[2][1]
			              *matrix.row[ii1].elem[ii2].v[2][1]/d[col*3+1]
			             + matrix.row[ii1].elem[ii2].v[2][2]
			              *matrix.row[ii1].elem[ii2].v[2][2]/d[col*3+2]);
		}
		if(d[ii1*3] < d_tol_0){
			if( (zero_flag0)&&(fabs(d[ii1*3]) >= d_tol_0) ){
				/* $BHs@5DjCMBP:v(B */
				/* d[ii1*3] *= -1.0; */
				d[ii1*3] = 1.0;
			}else{
				return(CAL_ERROR_RC);
			}
		}
		d[ii1*3+1] -= matrix.row[ii1].elem[dia_pos].v[1][0]
		             *matrix.row[ii1].elem[dia_pos].v[1][0]/d[ii1*3];
		if(d[ii1*3+1] < d_tol_1){
			if( (zero_flag1)&&(fabs(d[ii1*3+1]) >= d_tol_1) ){
				/* $BHs@5DjCMBP:v(B */
				/* d[ii1*3+1] *= -1.0; */
				d[ii1*3+1] = 1.0;
			}else{
				return(CAL_ERROR_RC);
			}
		}
		d[ii1*3+2] -= (matrix.row[ii1].elem[dia_pos].v[2][0]
		              *matrix.row[ii1].elem[dia_pos].v[2][0]/d[ii1*3]
		             + matrix.row[ii1].elem[dia_pos].v[2][1]
		              *matrix.row[ii1].elem[dia_pos].v[2][1]/d[ii1*3+1]);
		if(d[ii1*3+2] < d_tol_2){
			if( (zero_flag2)&&(fabs(d[ii1*3+2]) >= d_tol_2) ){
				/* $BHs@5DjCMBP:v(B */
				/* d[ii1*3+2] *= -1.0; */
				d[ii1*3+2] = 1.0;
			}else{
				return(CAL_ERROR_RC);
			}
		}
	}

	return(NORMAL_RC);
}


/* r = matrix * p $B$r6a;wE*$K2r$$$F(B p $B$r5a$a$k(B */
/* d $B$O(B matrix $B$NIT40A4%3%l%9%-!<J,2r8e$NBP3Q9`(B */
static RC ic_subst3s(NONZERO_MATRIX3 matrix, const long double *d,
                     const long double *r, long double *p)
{
	int ii1, ii2;

	/* $BA0?JBeF~(B */
	for(ii1=0; ii1<matrix.size; ii1++){
		int dia_pos = matrix.row[ii1].size - 1;
		double (*val)[3];
		long double p_ii1_3_0 = r[ii1*3];
		long double p_ii1_3_1 = r[ii1*3+1];
		long double p_ii1_3_2 = r[ii1*3+2];

		for(ii2=0; ii2<dia_pos; ii2++){
			int col = matrix.row[ii1].elem[ii2].col;

			val = matrix.row[ii1].elem[ii2].v;
			p_ii1_3_0 -= (p[col*3]   * val[0][0]
			            + p[col*3+1] * val[0][1]
			            + p[col*3+2] * val[0][2]);
			p_ii1_3_1 -= (p[col*3]   * val[1][0]
			            + p[col*3+1] * val[1][1]
			            + p[col*3+2] * val[1][2]);
			p_ii1_3_2 -= (p[col*3]   * val[2][0]
			            + p[col*3+1] * val[2][1]
			            + p[col*3+2] * val[2][2]);
		}
		val = matrix.row[ii1].elem[dia_pos].v;

		p[ii1*3] = p_ii1_3_0/d[ii1*3];

		p[ii1*3+1] = (p_ii1_3_1 - p[ii1*3] * val[1][0])/d[ii1*3+1];

		p[ii1*3+2] = (p_ii1_3_2 - (p[ii1*3]   * val[2][0]
		                         + p[ii1*3+1] * val[2][1]))/d[ii1*3+2];
	}

	/* $BBP3Q(B */
	for(ii1=0; ii1<matrix.size; ii1++){
		p[ii1*3]   *= d[ii1*3];
		p[ii1*3+1] *= d[ii1*3+1];
		p[ii1*3+2] *= d[ii1*3+2];
	}

	/* $B8eB`BeF~(B */
	for(ii1=matrix.size-1; ii1>=0; ii1--){
		int dia_pos = matrix.row[ii1].size - 1;
		long double p_ii1_3_0;
		long double p_ii1_3_1;
		long double p_ii1_3_2;

		p[ii1*3+2] /= d[ii1*3+2];

		p[ii1*3+1] -= p[ii1*3+2] * matrix.row[ii1].elem[dia_pos].v[2][1];
		p[ii1*3+1] /= d[ii1*3+1];

		p[ii1*3] -= (p[ii1*3+2] * matrix.row[ii1].elem[dia_pos].v[2][0]
		           + p[ii1*3+1] * matrix.row[ii1].elem[dia_pos].v[1][0]);
		p[ii1*3] /= d[ii1*3];

		p_ii1_3_0 = p[ii1*3];
		p_ii1_3_1 = p[ii1*3+1];
		p_ii1_3_2 = p[ii1*3+2];
		for(ii2=0; ii2<dia_pos; ii2++){
			int col = matrix.row[ii1].elem[ii2].col;
			double (*val)[3] = matrix.row[ii1].elem[ii2].v;

			p[col*3]   -= (p_ii1_3_0 * val[0][0]
			             + p_ii1_3_1 * val[1][0]
			             + p_ii1_3_2 * val[2][0]);
			p[col*3+1] -= (p_ii1_3_0 * val[0][1]
			             + p_ii1_3_1 * val[1][1]
			             + p_ii1_3_2 * val[2][1]);
			p[col*3+2] -= (p_ii1_3_0 * val[0][2]
			             + p_ii1_3_1 * val[1][2]
			             + p_ii1_3_2 * val[2][2]);
		}
	}

	return(NORMAL_RC);
}


RC nonzero3s_scaling_pre(NONZERO_MATRIX3 matrix,
                         double *scale_fact, double *vect)
{
	int ii1, ii2;

	for(ii1=0; ii1<matrix.size; ii1++){
		int dia_pos = matrix.row[ii1].size - 1;
		double (*val)[3] = matrix.row[ii1].elem[dia_pos].v;

		if(val[0][0] < ABS_ERROR){
			fprintf(stdout, "Warning: MATRIX[%d-%d]: %15.7e -> %15.7e\n",
			        ii1, 0, val[0][0], ABS_ERROR);
			val[0][0] = ABS_ERROR;
		}
		scale_fact[ii1*3]   = 1.0/sqrt(val[0][0]);

		if(val[1][1] < ABS_ERROR){
			fprintf(stdout, "Warning: MATRIX[%d-%d]: %15.7e -> %15.7e\n",
			        ii1, 1, val[1][1], ABS_ERROR);
			val[1][1] = ABS_ERROR;
		}
		scale_fact[ii1*3+1] = 1.0/sqrt(val[1][1]);

		if(val[2][2] < ABS_ERROR){
			fprintf(stdout, "Warning: MATRIX[%d-%d]: %15.7e -> %15.7e\n",
			        ii1, 2, val[2][2], ABS_ERROR);
			val[2][2] = ABS_ERROR;
		}
		scale_fact[ii1*3+2] = 1.0/sqrt(val[2][2]);

		for(ii2=0; ii2<=dia_pos; ii2++){
			int col = matrix.row[ii1].elem[ii2].col;

			val = matrix.row[ii1].elem[ii2].v;
			val[0][0] *= scale_fact[ii1*3]  *scale_fact[col*3];
			val[0][1] *= scale_fact[ii1*3]  *scale_fact[col*3+1];
			val[0][2] *= scale_fact[ii1*3]  *scale_fact[col*3+2];
			val[1][0] *= scale_fact[ii1*3+1]*scale_fact[col*3];
			val[1][1] *= scale_fact[ii1*3+1]*scale_fact[col*3+1];
			val[1][2] *= scale_fact[ii1*3+1]*scale_fact[col*3+2];
			val[2][0] *= scale_fact[ii1*3+2]*scale_fact[col*3];
			val[2][1] *= scale_fact[ii1*3+2]*scale_fact[col*3+1];
			val[2][2] *= scale_fact[ii1*3+2]*scale_fact[col*3+2];
		}
		vect[ii1*3]   *= scale_fact[ii1*3];
		vect[ii1*3+1] *= scale_fact[ii1*3+1];
		vect[ii1*3+2] *= scale_fact[ii1*3+2];
	}

	return(NORMAL_RC);
}


RC nonzero3s_scaling_post(int matrix_size,
                          const double *scale_fact, double *vect)
{
	int ii1;

	for(ii1=0; ii1<matrix_size; ii1++){
		vect[ii1*3]   *= scale_fact[ii1*3];
		vect[ii1*3+1] *= scale_fact[ii1*3+1];
		vect[ii1*3+2] *= scale_fact[ii1*3+2];
	}

	return(NORMAL_RC);
}


/* matrix $B$H(B vect $B$N@Q(B => ans (double $BHG(B) */
RC mul_nonzero_matrix_vect3s(NONZERO_MATRIX3 matrix,
                             const double *vect, double *ans)
{
	int ii1, ii2;

	for(ii1=0; ii1<matrix.size; ii1++){
		int dia_pos = matrix.row[ii1].size - 1;
		double (*val)[3];
		double ans_ii1_3_0 = 0.0;
		double ans_ii1_3_1 = 0.0;
		double ans_ii1_3_2 = 0.0;
		double vect_ii1_3_0 = vect[ii1*3];
		double vect_ii1_3_1 = vect[ii1*3+1];
		double vect_ii1_3_2 = vect[ii1*3+2];

		for(ii2=0; ii2<dia_pos; ii2++){
			int col = matrix.row[ii1].elem[ii2].col;

			val = matrix.row[ii1].elem[ii2].v;
			ans_ii1_3_0 += (vect[col*3]   * val[0][0]
			              + vect[col*3+1] * val[0][1]
			              + vect[col*3+2] * val[0][2]);
			ans_ii1_3_1 += (vect[col*3]   * val[1][0]
			              + vect[col*3+1] * val[1][1]
			              + vect[col*3+2] * val[1][2]);
			ans_ii1_3_2 += (vect[col*3]   * val[2][0]
			              + vect[col*3+1] * val[2][1]
			              + vect[col*3+2] * val[2][2]);

			ans[col*3]   += (vect_ii1_3_0 * val[0][0]
			               + vect_ii1_3_1 * val[1][0]
			               + vect_ii1_3_2 * val[2][0]);
			ans[col*3+1] += (vect_ii1_3_0 * val[0][1]
			               + vect_ii1_3_1 * val[1][1]
			               + vect_ii1_3_2 * val[2][1]);
			ans[col*3+2] += (vect_ii1_3_0 * val[0][2]
			               + vect_ii1_3_1 * val[1][2]
			               + vect_ii1_3_2 * val[2][2]);
		}
		val = matrix.row[ii1].elem[dia_pos].v;
		ans[ii1*3]   = ans_ii1_3_0 + (vect[ii1*3]   * val[0][0]
		                            + vect[ii1*3+1] * val[0][1]
		                            + vect[ii1*3+2] * val[0][2]);
		ans[ii1*3+1] = ans_ii1_3_1 + (vect[ii1*3]   * val[1][0]
		                            + vect[ii1*3+1] * val[1][1]
		                            + vect[ii1*3+2] * val[1][2]);
		ans[ii1*3+2] = ans_ii1_3_2 + (vect[ii1*3]   * val[2][0]
		                            + vect[ii1*3+1] * val[2][1]
		                            + vect[ii1*3+2] * val[2][2]);
	}

	return(NORMAL_RC);
}


/* matrix $B$H(B vect $B$N@Q(B => ans (long double $BHG(B) */
static RC mul_matrix_vect3sl(NONZERO_MATRIX3 matrix,
                             const long double *vect, long double *ans)
{
	int ii1, ii2;

	for(ii1=0; ii1<matrix.size; ii1++){
		int dia_pos = matrix.row[ii1].size - 1;
		double (*val)[3];
		long double ans_ii1_3_0 = 0.0;
		long double ans_ii1_3_1 = 0.0;
		long double ans_ii1_3_2 = 0.0;
		long double vect_ii1_3_0 = vect[ii1*3];
		long double vect_ii1_3_1 = vect[ii1*3+1];
		long double vect_ii1_3_2 = vect[ii1*3+2];

		for(ii2=0; ii2<dia_pos; ii2++){
			int col = matrix.row[ii1].elem[ii2].col;

			val = matrix.row[ii1].elem[ii2].v;
			ans_ii1_3_0 += (vect[col*3]   * val[0][0]
			              + vect[col*3+1] * val[0][1]
			              + vect[col*3+2] * val[0][2]);
			ans_ii1_3_1 += (vect[col*3]   * val[1][0]
			              + vect[col*3+1] * val[1][1]
			              + vect[col*3+2] * val[1][2]);
			ans_ii1_3_2 += (vect[col*3]   * val[2][0]
			              + vect[col*3+1] * val[2][1]
			              + vect[col*3+2] * val[2][2]);

			ans[col*3]   += (vect_ii1_3_0 * val[0][0]
			               + vect_ii1_3_1 * val[1][0]
			               + vect_ii1_3_2 * val[2][0]);
			ans[col*3+1] += (vect_ii1_3_0 * val[0][1]
			               + vect_ii1_3_1 * val[1][1]
			               + vect_ii1_3_2 * val[2][1]);
			ans[col*3+2] += (vect_ii1_3_0 * val[0][2]
			               + vect_ii1_3_1 * val[1][2]
			               + vect_ii1_3_2 * val[2][2]);
		}
		val = matrix.row[ii1].elem[dia_pos].v;
		ans[ii1*3]   = ans_ii1_3_0 + (vect[ii1*3]   * val[0][0]
		                            + vect[ii1*3+1] * val[0][1]
		                            + vect[ii1*3+2] * val[0][2]);
		ans[ii1*3+1] = ans_ii1_3_1 + (vect[ii1*3]   * val[1][0]
		                            + vect[ii1*3+1] * val[1][1]
		                            + vect[ii1*3+2] * val[1][2]);
		ans[ii1*3+2] = ans_ii1_3_2 + (vect[ii1*3]   * val[2][0]
		                            + vect[ii1*3+1] * val[2][1]
		                            + vect[ii1*3+2] * val[2][2]);
	}

	return(NORMAL_RC);
}


RC modify_nonzero3s(NONZERO_MATRIX3 matrix,
                    int index, int xyz, double v, double *vect)
{
	return( modify_nonzero3s_n(matrix, index, xyz, v, &vect, 1) );
}


RC modify_nonzero3s_n(NONZERO_MATRIX3 matrix,
                      int index, int xyz, double v, double **vect,
                      int vect_size)
{
	int ii1, ii2, ii3;
	int dia_index = matrix.row[index].size - 1;

	if((vect_size >= 1)&&(vect == NULL)) return(ARG_ERROR_RC);
	if((xyz < 0)||(xyz >= 3)) return(ARG_ERROR_RC);
	if((index < 0)||(index >= matrix.size)) return(ARG_ERROR_RC);
	
	for(ii1=0; ii1<matrix.row[index].size-1; ii1++){
		int col_index = matrix.row[index].elem[ii1].col;

		for(ii2=0; ii2<vect_size; ii2++){
			vect[ii2][col_index*3]   -= v*matrix.row[index].elem[ii1].v[xyz][0];
			vect[ii2][col_index*3+1] -= v*matrix.row[index].elem[ii1].v[xyz][1];
			vect[ii2][col_index*3+2] -= v*matrix.row[index].elem[ii1].v[xyz][2];
		}
		matrix.row[index].elem[ii1].v[xyz][0] = 0.0;
		matrix.row[index].elem[ii1].v[xyz][1] = 0.0;
		matrix.row[index].elem[ii1].v[xyz][2] = 0.0;
	}

	for(ii1=0; ii1<3; ii1++){
		for(ii2=0; ii2<vect_size; ii2++){
			vect[ii2][index*3+ii1]
			      -= v*matrix.row[index].elem[dia_index].v[xyz][ii1];
		}
		matrix.row[index].elem[dia_index].v[xyz][ii1] = 0.0;
		matrix.row[index].elem[dia_index].v[ii1][xyz] = 0.0;
	}
	matrix.row[index].elem[dia_index].v[xyz][xyz] = 1.0;
	/* matrix.row[index].elem[dia_index].v[xyz][xyz] = 10.0*ABS_TOL; */
	for(ii1=0; ii1<vect_size; ii1++){
		vect[ii1][index*3+xyz] = v;
	}

	for(ii1=0; ii1<matrix.row[index].syn_size; ii1++){
		int row_index = matrix.row[index].syn[ii1];

		for(ii2=0; ii2<matrix.row[row_index].size; ii2++){
			if(matrix.row[row_index].elem[ii2].col != index) continue;

			for(ii3=0; ii3<vect_size; ii3++){
				vect[ii3][row_index*3]
				     -= v*matrix.row[row_index].elem[ii2].v[0][xyz];
				vect[ii3][row_index*3+1]
				     -= v*matrix.row[row_index].elem[ii2].v[1][xyz];
				vect[ii3][row_index*3+2]
				     -= v*matrix.row[row_index].elem[ii2].v[2][xyz];
			}

			matrix.row[row_index].elem[ii2].v[0][xyz] = 0.0;
			matrix.row[row_index].elem[ii2].v[1][xyz] = 0.0;
			matrix.row[row_index].elem[ii2].v[2][xyz] = 0.0;
			break;
		}
	}

	return(NORMAL_RC);
}


/* index[0-2] $B9T!$Ns$K(B s_matrix[3][3], s_matrix[3][3]^T $B$r3]$1$k(B */
RC row_col_mul_nonzero3s(NONZERO_MATRIX3 matrix, int index,
                         double s_matrix[][3])
{
	int ii1;
	int elem_index;
	int syn_index;

	for(ii1=0; ii1<matrix.row[index].size; ii1++){
		mul_matrix_33(s_matrix, matrix.row[index].elem[ii1].v);
	}
	mul_matrix_t33(matrix.row[index].elem[matrix.row[index].size-1].v,
	               s_matrix);
	for(ii1=0; ii1<matrix.row[index].syn_size; ii1++){
		syn_index = matrix.row[index].syn[ii1];
		elem_index = search_col3(matrix, syn_index, index);
		RC_NEG_CHK( elem_index );
		mul_matrix_t33(matrix.row[syn_index].elem[elem_index].v, s_matrix);
	}

	return(NORMAL_RC);
}


/* A[3][3] * B[3][3] => B[3][3] */
static void mul_matrix_33(double A[][3], double B[][3])
{
	int ii1, ii2;
	double AB[3][3];

	for(ii1=0; ii1<3; ii1++){
		for(ii2=0; ii2<3; ii2++){
			AB[ii1][ii2] = A[ii1][0] * B[0][ii2]
			             + A[ii1][1] * B[1][ii2]
			             + A[ii1][2] * B[2][ii2];
		}
	}
	for(ii1=0; ii1<3; ii1++){
		for(ii2=0; ii2<3; ii2++){
			B[ii1][ii2] = AB[ii1][ii2];
		}
	}

}


/* A[3][3] * B^T[3][3] => A[3][3] */
static void mul_matrix_t33(double A[][3], double B[][3])
{
	int ii1, ii2;
	double AB[3][3];

	for(ii1=0; ii1<3; ii1++){
		for(ii2=0; ii2<3; ii2++){
			AB[ii1][ii2] = A[ii1][0] * B[ii2][0]
			             + A[ii1][1] * B[ii2][1]
			             + A[ii1][2] * B[ii2][2];
		}
	}
	for(ii1=0; ii1<3; ii1++){
		for(ii2=0; ii2<3; ii2++){
			A[ii1][ii2] = AB[ii1][ii2];
		}
	}
}


/* index1,xyz1 $B9T!$Ns$N(B weight $BG\$r(B index2,xyz2 $B9T!$Ns$K2C;;(B */
RC row_col_add_nonzero3s(NONZERO_MATRIX3 matrix, int index1, int xyz1,
                         double weight, int index2, int xyz2)
{
	int ii1, ii2;
	int dia_index1;
	int dia_index2;

	if((xyz1 < 0)||(xyz1 >= 3)) return(ARG_ERROR_RC);
	if((xyz2 < 0)||(xyz2 >= 3)) return(ARG_ERROR_RC);
	if((index1 < 0)||(index1 >= matrix.size)) return(ARG_ERROR_RC);
	if((index2 < 0)||(index2 >= matrix.size)) return(ARG_ERROR_RC);
	if((index1 == index2)&&(xyz1 == xyz2)) return(ARG_ERROR_RC);


	/* $BBP3Q9`$N%$%s%G%C%/%9(B  $B<B9TCf$KJQ2=$9$k$N$GC`<!:F7W;;$,I,MW(B */
	dia_index1 = matrix.row[index1].size - 1;
	dia_index2 = matrix.row[index2].size - 1;

	/* $BBP3Q9`F1;N$N7W;;(B */
	matrix.row[index2].elem[dia_index2].v[xyz2][xyz2]
	    += weight * weight * matrix.row[index1].elem[dia_index1].v[xyz1][xyz1];


	for(ii1=0; ii1<matrix.row[index1].size-1; ii1++){
		int col_index1 = matrix.row[index1].elem[ii1].col;

		if(col_index1 < index2){
			NONZERO_ELEM3 *elem_index2;

			RC_NULL_CHK( elem_index2
			             = nonzero_matrix3s_ptr(matrix, index2, col_index1) );
			for(ii2=0; ii2<3; ii2++){
				elem_index2->v[xyz2][ii2]
				    += weight * matrix.row[index1].elem[ii1].v[xyz1][ii2];
			}
		}else if(col_index1 == index2){
			dia_index2 = matrix.row[index2].size - 1;
			for(ii2=0; ii2<3; ii2++){
				matrix.row[index2].elem[dia_index2].v[xyz2][ii2]
				    += weight * matrix.row[index1].elem[ii1].v[xyz1][ii2];
				matrix.row[index2].elem[dia_index2].v[ii2][xyz2]
				    += weight * matrix.row[index1].elem[ii1].v[xyz1][ii2];
			}
		}else{
			NONZERO_ELEM3 *elem_col_index1;

			RC_NULL_CHK( elem_col_index1
			             = nonzero_matrix3s_ptr(matrix, col_index1, index2) );
			for(ii2=0; ii2<3; ii2++){
				elem_col_index1->v[ii2][xyz2]
				    += weight * matrix.row[index1].elem[ii1].v[xyz1][ii2];
			}
		}
	}

	if(index1 < index2){
		NONZERO_ELEM3 *elem_index2;

		RC_NULL_CHK( elem_index2
		             = nonzero_matrix3s_ptr(matrix, index2, index1) );

		for(ii2=0; ii2<xyz1; ii2++){
			elem_index2->v[xyz2][ii2]
			    += weight * matrix.row[index1].elem[dia_index1].v[xyz1][ii2];
		}
		for(ii2=xyz1+1; ii2<3; ii2++){
			elem_index2->v[xyz2][ii2]
			    += weight * matrix.row[index1].elem[dia_index1].v[ii2][xyz1];
		}
	}else if(index1 == index2){
		for(ii2=0; ii2<3; ii2++){
			if(ii2 == xyz1) continue;
			matrix.row[index1].elem[dia_index1].v[xyz2][ii2]
			    += weight * matrix.row[index1].elem[dia_index1].v[xyz1][ii2];
		}
		for(ii2=0; ii2<3; ii2++){
			if(ii2 == xyz1) continue;
			matrix.row[index1].elem[dia_index1].v[ii2][xyz2]
			    += weight * matrix.row[index1].elem[dia_index1].v[ii2][xyz1];
		}
	}else{
		NONZERO_ELEM3 *elem_index1;
		RC_NULL_CHK (elem_index1
		             = nonzero_matrix3s_ptr(matrix, index1, index2) );
		dia_index1 = matrix.row[index1].size - 1;

		for(ii2=0; ii2<xyz1; ii2++){
			elem_index1->v[ii2][xyz2]
			    += weight * matrix.row[index1].elem[dia_index1].v[xyz1][ii2];
		}
		for(ii2=xyz1+1; ii2<3; ii2++){
			elem_index1->v[ii2][xyz2]
			    += weight * matrix.row[index1].elem[dia_index1].v[ii2][xyz1];
		}
	}

	for(ii1=0; ii1<matrix.row[index1].syn_size; ii1++){
		int syn_index1 = matrix.row[index1].syn[ii1];
		int elem_syn_index1;

		if(syn_index1 < index2){
			NONZERO_ELEM3 *elem_index2;

			RC_NULL_CHK( elem_index2
			           = nonzero_matrix3s_ptr(matrix, index2, syn_index1) );
			RC_NEG_CHK( elem_index2 );
			elem_syn_index1 = search_col3(matrix, syn_index1, index1);
			RC_NEG_CHK( elem_syn_index1 );
			for(ii2=0; ii2<3; ii2++){
				elem_index2->v[xyz2][ii2]
				    += weight * matrix.row[syn_index1]
				                      .elem[elem_syn_index1].v[ii2][xyz1];
			}
		}else if(syn_index1 == index2){
			dia_index2 = matrix.row[index2].size - 1;
			elem_syn_index1 = search_col3(matrix, syn_index1, index1);
			RC_NEG_CHK( elem_syn_index1 );
			for(ii2=0; ii2<3; ii2++){
				matrix.row[index2].elem[dia_index2].v[xyz2][ii2]
				    += weight * matrix.row[syn_index1]
				                      .elem[elem_syn_index1].v[ii2][xyz1];
				matrix.row[index2].elem[dia_index2].v[ii2][xyz2]
				    += weight * matrix.row[syn_index1]
				                      .elem[elem_syn_index1].v[ii2][xyz1];
			}
		}else{
			NONZERO_ELEM3 *elem_syn_index12;

			RC_NULL_CHK( elem_syn_index12
			             = nonzero_matrix3s_ptr(matrix, syn_index1, index2) );
			elem_syn_index1 = search_col3(matrix, syn_index1, index1);
			RC_NEG_CHK( elem_syn_index1 );
			for(ii2=0; ii2<3; ii2++){
				elem_syn_index12->v[ii2][xyz2]
				    += weight * matrix.row[syn_index1]
				                      .elem[elem_syn_index1].v[ii2][xyz1];
			}
		}
	}

	return(NORMAL_RC);
}


/* matrix $B$N(B index $B9T$G(B col $BNsL\$NMWAG$N%$%s%G%C%/%9$rJV5Q(B */
static int search_col3(NONZERO_MATRIX3 matrix, int index, int col)
{
	int ii1;

	for(ii1=0; ii1<matrix.row[index].size; ii1++){
		if(matrix.row[index].elem[ii1].col == col){
			return(ii1);
		}
	}

	return(-1);
}


/* index,xyz $B9T(B($BNs(B)$B$rI=<((B */
RC row_col_print_nonzero3s(FILE *fp, NONZERO_MATRIX3 matrix, int index, int xyz)
{
	int ii1, ii2;

	if((xyz < 0)||(xyz >= 3)) return(ARG_ERROR_RC);
	if((index < 0)||(index >= matrix.size)) return(ARG_ERROR_RC);
	
	for(ii1=0; ii1<matrix.row[index].size; ii1++){
		fprintf(fp, "%10.2e", matrix.row[index].elem[ii1].v[xyz][0]);
		fprintf(fp, "%10.2e", matrix.row[index].elem[ii1].v[xyz][1]);
		fprintf(fp, "%10.2e", matrix.row[index].elem[ii1].v[xyz][2]);
	}

	for(ii1=0; ii1<matrix.row[index].syn_size; ii1++){
		int row_index = matrix.row[index].syn[ii1];
		if((row_index < 0)||(row_index >= matrix.size)) return(ARG_ERROR_RC);

		for(ii2=0; ii2<matrix.row[row_index].size; ii2++){
			if(matrix.row[row_index].elem[ii2].col != index) continue;

			fprintf(fp, "%10.2e", matrix.row[row_index].elem[ii2].v[0][xyz]);
			fprintf(fp, "%10.2e", matrix.row[row_index].elem[ii2].v[1][xyz]);
			fprintf(fp, "%10.2e", matrix.row[row_index].elem[ii2].v[2][xyz]);
			break;
		}
	}

	return(NORMAL_RC);
}


/* index,xyz $B9T$HNs$r%<%m$rBeF~!%BP3Q9`$O(B 1.0 $B$rBeF~(B */
RC row_col_zero_nonzero3s(NONZERO_MATRIX3 matrix, int index, int xyz)
{
	int ii1, ii2;
	int dia_index = matrix.row[index].size - 1;

	if((xyz < 0)||(xyz >= 3)) return(ARG_ERROR_RC);
	if((index < 0)||(index >= matrix.size)) return(ARG_ERROR_RC);
	
	for(ii1=0; ii1<matrix.row[index].size-1; ii1++){
		matrix.row[index].elem[ii1].v[xyz][0] = 0.0;
		matrix.row[index].elem[ii1].v[xyz][1] = 0.0;
		matrix.row[index].elem[ii1].v[xyz][2] = 0.0;
	}

	for(ii1=0; ii1<3; ii1++){
		matrix.row[index].elem[dia_index].v[xyz][ii1] = 0.0;
		matrix.row[index].elem[dia_index].v[ii1][xyz] = 0.0;
	}
	/*matrix.row[index].elem[dia_index].v[xyz][xyz] = 1.0;*/
	matrix.row[index].elem[dia_index].v[xyz][xyz] = 10.0*ABS_TOL;

	for(ii1=0; ii1<matrix.row[index].syn_size; ii1++){
		int row_index = matrix.row[index].syn[ii1];
		if((row_index < 0)||(row_index >= matrix.size)) return(ARG_ERROR_RC);

		for(ii2=0; ii2<matrix.row[row_index].size; ii2++){
			if(matrix.row[row_index].elem[ii2].col == index){
				matrix.row[row_index].elem[ii2].v[0][xyz] = 0.0;
				matrix.row[row_index].elem[ii2].v[1][xyz] = 0.0;
				matrix.row[row_index].elem[ii2].v[2][xyz] = 0.0;
				break;
			}
		}
	}

	return(NORMAL_RC);
}


/*
 * $B$3$3$+$i$OHsBP>N%^%H%j%C%/%9@lMQ(B
 */


/* A^T[][] * A[][] * vector[] = ans[] */
static void mul_nonzero_matrix_ATAb(NONZERO_MATRIX matrix, const double *vector,
                                    double *ans)
{
	int ii1, ii2;
	double tmp;

	for(ii1=0; ii1<matrix.size; ii1++){
		ans[ii1] = 0.0;
	}

	for(ii1=0; ii1<matrix.size; ii1++){
		tmp = 0.0;
		for(ii2=0; ii2<matrix.row[ii1].size; ii2++){
			tmp += matrix.row[ii1].v[ii2]
			     * vector[ matrix.row[ii1].col[ii2] ];
		}
		for(ii2=0; ii2<matrix.row[ii1].size; ii2++){
			ans[ matrix.row[ii1].col[ii2] ] += matrix.row[ii1].v[ii2] * tmp;
		}
	}
}


void print_nonzero_matrix1a(FILE *fp, NONZERO_MATRIX matrix)
{
	int ii1, ii2;

	for(ii1=0; ii1<matrix.size; ii1++){
		fprintf(fp, "size %d\n", matrix.row[ii1].size);
		fprintf(fp, "diagonal_index %d\n", matrix.row[ii1].diagonal_index);
		for(ii2=0; ii2<matrix.row[ii1].size; ii2++){
			fprintf(fp, "row %d:   col %d:   value %e\n",
			             ii1, matrix.row[ii1].col[ii2],
			             matrix.row[ii1].v[ii2]);
		}
		fprintf(fp, "\n");
	}
}


/* matrix[][] * vect[] = ans[] */
void mul_nonzero_matrix_vect1a(NONZERO_MATRIX matrix, const double *vect,
                               double *ans)
{
	int ii1, ii2;
	
	for(ii1=0; ii1<matrix.size; ii1++){
		ans[ii1] = 0.0;
		for(ii2=0; ii2<matrix.row[ii1].size; ii2++){
			ans[ii1] += matrix.row[ii1].v[ii2]
			          * vect[ matrix.row[ii1].col[ii2] ];
		}
	}
}


/* matrix^T[][] * vect[] = ans[] */
void mul_nonzero_matrixT_vect1a(NONZERO_MATRIX matrix, const double *vect,
                                double *ans)
{
	int ii1, ii2;
	
	for(ii1=0; ii1<matrix.size; ii1++){
		ans[ii1] = 0.0;
	}
	for(ii1=0; ii1<matrix.size; ii1++){
		for(ii2=0; ii2<matrix.row[ii1].size; ii2++){
			ans[matrix.row[ii1].col[ii2]] += matrix.row[ii1].v[ii2] * vect[ii1];
		}
	}
}


/* matrix.size = total_dof, matrix.row[].size = 0 $B$,BeF~$5$l$k(B */
RC allocate_nonzero_matrix1a(int total_dof, NONZERO_MATRIX *matrix)
{
	int ii1;

	matrix->size = total_dof;
	RC_NULL_CHK( matrix->row = (NONZERO_ROW *)malloc(matrix->size
	                                               * sizeof(NONZERO_ROW) ) );
	for(ii1=0; ii1<matrix->size; ii1++){
		matrix->row[ii1].size = 0;
		matrix->row[ii1].diagonal_index = -1;
		matrix->row[ii1].col = NULL;
		matrix->row[ii1].v = NULL;
	}

	return(NORMAL_RC);
}


/* a[row][col] $B$,$9$G$K3NJ]$7$F$"$k$J$i!$$=$NCM$N%]%$%s%?!<$r!$(B */
/* $B3NJ]$7$F$J$$$J$i!$?7$?$K3NJ]$7$F$=$N%]%$%s%?!<$rJV$9!%(B       */
/* $B$3$N4X?t;HMQ;~$O!$(Ballocate_nonzero_matrix() $B$r;vA0$KDL$9!%(B   */
/* $BITHw$,@8$8$?>l9g!$(BNULL $B$,La$jCM$H$J$k!%(B                      */
double *nonzero_matrix1a_ptr(NONZERO_MATRIX *matrix, int row, int col)
{
#if 0
	int ii1;
	
	if( (0 > matrix->size) || (matrix->size <= row) ) return(NULL);

	for(ii1=0; ii1<matrix->row[row].size; ii1++){
		if(col == matrix->row[row].col[ii1]){
			return( &(matrix->row[row].v[ii1]) );
		}
	}

	(matrix->row[row].size)++;	

	matrix->row[row].col = (int *)realloc(matrix->row[row].col,
	                               matrix->row[row].size * sizeof(int) );
	if(matrix->row[row].col == NULL) return(NULL);
	matrix->row[row].col[matrix->row[row].size-1] = col;

	matrix->row[row].v = (double *)realloc(matrix->row[row].v,
	                                matrix->row[row].size * sizeof(double) );
	if(matrix->row[row].v == NULL) return(NULL);
	matrix->row[row].v[matrix->row[row].size-1] = 0.0;
		
	return( &(matrix->row[row].v[matrix->row[row].size-1]) );
#endif

	int ii1;

	if( (0 > matrix->size) || (matrix->size <= row) ) return(NULL);

	/* $B$3$l0JA0$KCM$,BeF~$5$l$F$$$k(B */
	for(ii1=matrix->row[row].size-1; ii1>=0; ii1--){
		if(col > matrix->row[row].col[ii1]) break;
		if(col == matrix->row[row].col[ii1]){
			return( &(matrix->row[row].v[ii1]) );
		}
	}

	(matrix->row[row].size)++;
	matrix->row[row].col = (int *)realloc(matrix->row[row].col,
	                               matrix->row[row].size * sizeof(int) );
	if(matrix->row[row].col == NULL) return(NULL);
	matrix->row[row].v = (double *)realloc(matrix->row[row].v,
	                                matrix->row[row].size * sizeof(double) );
	if(matrix->row[row].v == NULL) return(NULL);
		
	/* $B=i$a$FCM$,BeF~$5$l$k(B */
	for(ii1=matrix->row[row].size-1; ii1>=1; ii1--){
		if(col > matrix->row[row].col[ii1-1]){
			matrix->row[row].col[ii1] = col;
			matrix->row[row].v[ii1] = 0.0;
			if(row == col){
				matrix->row[row].diagonal_index = ii1;
			}else if(matrix->row[row].diagonal_index >= ii1){
				matrix->row[row].diagonal_index++;
			}
			return( &(matrix->row[row].v[ii1]) );
		}
		matrix->row[row].col[ii1] = matrix->row[row].col[ii1-1];
		matrix->row[row].v[ii1] = matrix->row[row].v[ii1-1];
	}
	matrix->row[row].col[0] = col;
	matrix->row[row].v[0] = 0.0;
	if(row == col){
		matrix->row[row].diagonal_index = 0;
	}else if(matrix->row[row].diagonal_index >= 0){
		matrix->row[row].diagonal_index++;
	}
		
	return( &(matrix->row[row].v[0]) );
}


RC free_nonzero_matrix1a(NONZERO_MATRIX *matrix)
{
	int ii1;

	if(matrix->size <= 0) return(ARG_ERROR_RC);
	if(matrix->row == NULL) return(ARG_ERROR_RC);

	for(ii1=0; ii1<matrix->size; ii1++){
		if(matrix->row[ii1].size > 0){
			if(matrix->row[ii1].col == NULL) return(ARG_ERROR_RC);
			if(matrix->row[ii1].v == NULL) return(ARG_ERROR_RC);
			free( (void *)matrix->row[ii1].col );
			free( (void *)matrix->row[ii1].v );
		}
		matrix->row[ii1].size = 0;
		matrix->row[ii1].diagonal_index = -1;
		matrix->row[ii1].col = NULL;
		matrix->row[ii1].v = NULL;
	}
	free( (void *)matrix->row );
	matrix->size = 0;
	matrix->row = NULL;

	return(NORMAL_RC);
}


RC copy_nonzero_matrix1a(NONZERO_MATRIX source, NONZERO_MATRIX *copy)
{
	int ii1, ii2;
	double *ptr;

	RC_TRY( allocate_nonzero_matrix1a(source.size, copy) );
	for(ii1=0; ii1<source.size; ii1++){
		copy->row[ii1].diagonal_index = source.row[ii1].diagonal_index;
		for(ii2=0; ii2<source.row[ii1].size; ii2++){
			ptr = nonzero_matrix1a_ptr(copy, ii1, source.row[ii1].col[ii2]);
			RC_NULL_CHK(ptr);
			*ptr = source.row[ii1].v[ii2];
		}
	}

	return(NORMAL_RC);
}


/* $B%^%H%j%C%/%9$N3F9T$N:GBgCM$r8+$D$1!$$=$NCM$G%9%1!<%j%s%0(B */
/* $B%^%H%j%C%/%9$NCM!$%Y%/%H%k$NCM$O=$@5$5$l$k(B */
RC nonzero1a_scaling(NONZERO_MATRIX *matrix, double *vector)
{
	double max;
	double tmp_d;
	int tmp_i;
	int ii1, ii2;

	for(ii1=0; ii1<matrix->size; ii1++){
		tmp_i = 0;
		for(ii2=0; ii2<matrix->row[ii1].size; ii2++){
			/* $B5-21$5$l$F$$$kCM$,(B 0 $B$@$C$?$i(B */
			if(fabs(matrix->row[ii1].v[ii2]) < DBL_MIN){
				/* $B$=$N>l=j$,BP3Q9`$@$C$?$i(B */
				if(matrix->row[ii1].col[ii2] == ii1){
					matrix->row[ii1].diagonal_index = -1;
				}
				continue;
			}
			matrix->row[ii1].col[tmp_i] = matrix->row[ii1].col[ii2];
			matrix->row[ii1].v[tmp_i] = matrix->row[ii1].v[ii2];
			/* $BBP3Q9`$N5-21>l=j$,0\F0$7$?$i(B */
			if(matrix->row[ii1].col[ii2] == ii1){
				matrix->row[ii1].diagonal_index = tmp_i;
			}
			tmp_i++;
		}
		matrix->row[ii1].size = tmp_i;
		matrix->row[ii1].col = (int *)realloc(matrix->row[ii1].col,
	                               matrix->row[ii1].size * sizeof(int) );
		RC_NULL_CHK(matrix->row[ii1].col);
		matrix->row[ii1].v = (double *)realloc(matrix->row[ii1].v,
	                                matrix->row[ii1].size * sizeof(double) );
		RC_NULL_CHK(matrix->row[ii1].v);

		max = 0.0;
		for(ii2=0; ii2<matrix->row[ii1].size; ii2++){
			if( fabs(max) < fabs(matrix->row[ii1].v[ii2]) ){
				max = matrix->row[ii1].v[ii2];
			}
		}
		/* $B9T$NA4$F$NCM$,(B 0 $B$@$C$?$i(B */
		if(fabs(max) < ABS_ERROR) return(CAL_ERROR_RC);
		tmp_d = 1.0 / max;
		for(ii2=0; ii2<matrix->row[ii1].size; ii2++){
			matrix->row[ii1].v[ii2] *= tmp_d;
		}
		vector[ii1] *= tmp_d;
	}

	return(NORMAL_RC);
}


/* x[] $B$K=i4|CM$rBeF~$9$k$H!$<B9T8e$O(B( [matrix]{x} = {b} )$B$N2r$,BeF~$5$l$k(B */
RC cg_nonzero1a_ver1(NONZERO_MATRIX matrix, const double *b, double *x)
{
	double *r; 
	double *p;
	double *tmp;
	double alpha;
	double beta;
	double error;
	long double r_p;
	long double p_Ap;
	long double b_b;
	long double r_r;
	long double r0_r0;
	int counter = 0;
	int conv_count = 0;
	int ii1;
	
	RC_NULL_CHK( r = (double *)malloc( matrix.size * sizeof(double) ) );
	RC_NULL_CHK( p = (double *)malloc( matrix.size * sizeof(double) ) );
	RC_NULL_CHK( tmp = (double *)malloc( matrix.size * sizeof(double) ) );

	mul_nonzero_matrix_vect1a(matrix, x, tmp);
	for(ii1=0; ii1<matrix.size; ii1++){
		tmp[ii1] = b[ii1] - tmp[ii1];
	}
	mul_nonzero_matrixT_vect1a(matrix, tmp, r);

	b_b = 0.0;
	r0_r0 = 0.0;
	for(ii1=0; ii1<matrix.size; ii1++){
		b_b += b[ii1] * b[ii1];
		r0_r0 += r[ii1] * r[ii1];
		p[ii1] = r[ii1];
	}
	if(b_b < ABS_ERROR){
		fprintf(stderr, "\n");
		return(CAL_ERROR_RC);
	}
	error = sqrt( (r0_r0/b_b)/matrix.size );
	
	fprintf(stderr, "Computing CG...");

	/* $B8m:9$,(B REL_ERROR1 $B0J2<$H$J$k$+(B 10 $B2sO"B3$G(B REL_ERROR2 $B0J2<(B */
	/* $B$H$J$C$?;~E@$G=*N;(B                                         */
	while( (error > REL_ERROR1) && (conv_count < 10) ){

/*
fprintf(stderr, "\r%d : %7.3e", counter, error);
*/

		mul_nonzero_matrix_ATAb(matrix, p, tmp);

		r_p = p_Ap = 0.0;
		for(ii1=0; ii1<matrix.size; ii1++){
			r_p += r[ii1] * p[ii1];
			p_Ap += p[ii1] * tmp[ii1];
		}
		if(p_Ap < ABS_ERROR){
			fprintf(stderr, "\n");
			return(CAL_ERROR_RC);
		}
		alpha = r_p / p_Ap;

		r_r = 0.0;
		for(ii1=0; ii1<matrix.size; ii1++){
			x[ii1] += alpha * p[ii1];
			r[ii1] -= alpha * tmp[ii1];
			r_r += r[ii1] * r[ii1];
		}
		beta = r_r / r0_r0;
		r0_r0 = r_r;

		error = sqrt( (r_r/b_b)/matrix.size );
		if(error <= REL_ERROR2){
			conv_count++;
		}else{
			conv_count = 0;
		}

		for(ii1=0; ii1<matrix.size; ii1++){
			p[ii1] = r[ii1] + beta * p[ii1];
		}

		counter ++;

	}

	fprintf(stderr, "\rCG Convergence(loop = %d)\n", counter);

	free( (void *)r );
	free( (void *)p );
	free( (void *)tmp );

	return(NORMAL_RC);
}


/* |r| $B%N%k%`$,3N<B$K2<$,$C$F$$$/HsBP>NMQ(B cg */
RC cg_nonzero1a_ver2(NONZERO_MATRIX matrix, const double *vector, double *ans)
{
	double *r; 
	double *p;
	double *tmp;
	double alpha;
	double beta;
	double error;
	long double c;
	long double c_before;
	long double d;
	long double r_r;
	long double b_b;
	int counter = 0;
	int conv_count = 0;
	int ii1;
	
	RC_NULL_CHK( r = (double *)malloc( matrix.size * sizeof(double) ) );
	RC_NULL_CHK( p = (double *)malloc( matrix.size * sizeof(double) ) );
	RC_NULL_CHK( tmp = (double *)malloc( matrix.size * sizeof(double) ) );

	mul_nonzero_matrix_vect1a(matrix, ans, tmp);

	b_b = 0.0;
	r_r = 0.0;
	for(ii1=0; ii1<matrix.size; ii1++){
		r[ii1] = vector[ii1] - tmp[ii1];
		r_r += r[ii1] * r[ii1];
		b_b += vector[ii1] * vector[ii1];
		p[ii1] = 0.0;
	}

	error = sqrt( (r_r/b_b)/matrix.size );

	c_before = 1.0;

	fprintf(stderr, "Computing CG...");
	while( (error > REL_ERROR1) && (conv_count < 10) ){

/*
fprintf(stderr, "\r%d : %7.5e", counter, error);
*/

		mul_nonzero_matrixT_vect1a(matrix, r, tmp);
		c = 0.0;
		for(ii1=0; ii1<matrix.size; ii1++){
			c += tmp[ii1] * tmp[ii1];
		}

		if( c_before <= ABS_ERROR ){
			fprintf(stderr, "\n");
			return(CAL_ERROR_RC);
		}
		beta = c / c_before;

		for(ii1=0; ii1<matrix.size; ii1++){
			p[ii1] = tmp[ii1] + beta * p[ii1];
		}

		mul_nonzero_matrix_vect1a(matrix, p, tmp);
		d = 0.0;
		for(ii1=0; ii1<matrix.size; ii1++){
			d += tmp[ii1] * tmp[ii1];
		}

		if( d <= ABS_ERROR ){
			fprintf(stderr, "\n");
			return(CAL_ERROR_RC);
		}
		alpha = c / d;
		c_before = c;

		r_r = 0.0;
		for(ii1=0; ii1<matrix.size; ii1++){
			ans[ii1] += alpha * p[ii1];
			r[ii1] -= alpha * tmp[ii1];
			r_r += r[ii1] * r[ii1];
		}

		error = sqrt( (r_r/b_b)/matrix.size );
		if(error <= REL_ERROR2){
			conv_count++;
		/* $B$3$N(B CG $B$O?6F0$7$J$$(B
		}else{
			conv_count = 0;
		*/
		}

		counter ++;
	}

	fprintf(stderr, "\rCG Convergence(loop = %d)\n", counter);

	free( (void *)r );
	free( (void *)p );
	free( (void *)tmp );

	return(NORMAL_RC);
}
