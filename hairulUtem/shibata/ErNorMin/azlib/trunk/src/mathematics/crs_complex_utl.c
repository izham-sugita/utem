/*********************************************************************
 * crs_complex_utl.c
 *
 * Copyright (C) 2016 AzLib Developers Group
 *
 * Written by
 *  <Takuya HAYASHI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: crs_complex_utl.c 1095 2016-12-05 07:29:35Z hayashi $ */

#include <stdio.h>
#include <math.h>
#include "rc.h"
#include "memory_manager.h"
#include "math_utl.h"
#include "complex_utl.h"


/* CRS行列 matrix と複素ベクトル vect の積を計算する．
 *   ret[m] = matrix[m*n] * vect[n]
 */
RC mul_crs_cmatrix_cvect (CRS_CMATRIX matrix, COMPLEX *vect, 
        COMPLEX *ret)
{
	RC_NULL_CHK( vect ); 
	RC_NULL_CHK( ret );

	int ii1, ii2; 

	for(ii1=0; ii1<matrix.size; ii1++){
		ret[ii1] = 0.0; 
		for(ii2=0; ii2<matrix.row[ii1].size; ii2++){
			ret[ii1] += matrix.row[ii1].v[ii2] * vect[matrix.row[ii1].col[ii2]];
		}
	}	

	return(NORMAL_RC);
}

/* CRS行列 matrix と複素ベクトル vect の積を複素ベクトル ret に足す．
 *   ret[m] += matrix[m*n] * vect[n]
 */
RC mul_add_crs_cmatrix_cvect (CRS_CMATRIX matrix, COMPLEX *vect, COMPLEX *ret)
{
	RC_NULL_CHK( vect ); 
	RC_NULL_CHK( ret );

	int ii1, ii2; 

	for(ii1=0; ii1<matrix.size; ii1++){
		for(ii2=0; ii2<matrix.row[ii1].size; ii2++){
			ret[ii1] += matrix.row[ii1].v[ii2] * vect[matrix.row[ii1].col[ii2]];
		}
	}	

	return(NORMAL_RC);
}

/* CRS行列 source を dest にコピーする．
 * dest はあらかじめ確保しておく．
 */
RC copy_crs_cmatrix (CRS_CMATRIX source, CRS_CMATRIX *dest)
{
	int ii1, ii2;
	COMPLEX *ptr;

	init_crs_cmatrix(dest);

	dest->size = source.size;
	for(ii1=0; ii1<source.size; ii1++){
		dest->row[ii1].diagonal_index = source.row[ii1].diagonal_index;
		for(ii2=0; ii2<source.row[ii1].size; ii2++){
			ptr = crs_cmatrix_ptr(dest, ii1, source.row[ii1].col[ii2]);
			RC_NULL_CHK( ptr );
			*ptr = source.row[ii1].v[ii2];
		}
	}

	return(NORMAL_RC);
}


/* CRS行列 matrix を fp に出力する．
 */
RC print_crs_cmatrix (FILE *fp, CRS_CMATRIX matrix)
{
	RC_NULL_CHK( fp );

	int ii1;

	for(ii1=0; ii1<matrix.size; ii1++){
		fprintf(fp, "\nrow %d\n", ii1);
		RC_TRY( print_crs_crow(fp, matrix.row[ii1]) );
	}

	return(NORMAL_RC);
}

/* CRS行列 row を fp に出力する．
 */
RC print_crs_crow(FILE *fp, CRS_CROW row)
{
	RC_NULL_CHK( fp );

	int ii1;

	fprintf(fp, "size %d\n", row.size);
	fprintf(fp, "diagonal_index %d\n", row.diagonal_index);
	for(ii1=0; ii1<row.size; ii1++){
		fprintf(fp, "  col %d:   value %15.7e + %15.7e i\n",
			row.col[ii1], creal(row.v[ii1]), cimag(row.v[ii1]));
	}

	return(NORMAL_RC);
}


/* CRS行列を動的確保する．
 * 行方向は dof 要素分確保する．
 * 列方向はここでは確保しない．crs_cmatrix_ptr()で値を代入するときに
 * 確保する．
 */
RC allocate_crs_cmatrix(int dof, CRS_CMATRIX *matrix)
{
    RC_NEG_ZERO_CHK( dof );
    RC_NULL_CHK( matrix );

    int ii1;

    matrix->size = dof;
    matrix->row = (CRS_CROW *)mm_alloc(matrix->size*sizeof(CRS_CROW));
    if(matrix->row == NULL) return(ALLOC_ERROR_RC);

    for(ii1=0; ii1<(matrix->size); ii1++){
        matrix->row[ii1].size = 0;
        matrix->row[ii1].diagonal_index = -1;
		matrix->row[ii1].col = NULL;
		matrix->row[ii1].v = NULL;
    }

    return(NORMAL_RC);
}

/* 動的確保したCRS行列を解放する．
 */
RC free_crs_cmatrix(CRS_CMATRIX *matrix)
{
    RC_NULL_CHK( matrix );

    int ii1;

    RC_NEG_ZERO_CHK( matrix->size );
    RC_NULL_CHK( matrix->row );

    for(ii1=0; ii1<(matrix->size); ii1++){
        if(matrix->row[ii1].size > 0){
            RC_NULL_CHK( matrix->row[ii1].col );
            RC_NULL_CHK( matrix->row[ii1].v );
            RC_TRY( mm_free( (void *)matrix->row[ii1].col) );
            RC_TRY( mm_free( (void *)matrix->row[ii1].v) );
        }
        matrix->row[ii1].size = 0;
        matrix->row[ii1].diagonal_index = -1;
        matrix->row[ii1].col = NULL;
        matrix->row[ii1].v = NULL;
    }
    RC_TRY( mm_free((void *)matrix->row) );
    matrix->size = 0;
    matrix->row = NULL;

    return(NORMAL_RC);
}

/* CRS行列を初期化する．
 */
RC init_crs_cmatrix(CRS_CMATRIX *matrix)
{
    RC_NULL_CHK( matrix );

    int ii1;

    for(ii1=0; ii1<(matrix->size); ii1++){
		if(matrix->row[ii1].size > 0){
			RC_TRY( free1C(matrix->row[ii1].size, &(matrix->row[ii1].v)) );
			RC_TRY( free1I(matrix->row[ii1].size, &(matrix->row[ii1].col)) );
		}
        matrix->row[ii1].size = 0;
        matrix->row[ii1].diagonal_index = -1;
    }

    return(NORMAL_RC);
}

/* row[col] がすでに確保してあるなら，その値のポインタを返す．
 * 確保してないなら，新たに確保してそのポインタを返す．
 * この関数使用時は，allocate_crs_cmatrix() を事前に通す．
 * 不備が生じた場合は，NULL が戻り値となる．
 */
COMPLEX *crs_crow_ptr (CRS_CROW *n_row, int row, int col)
{
	int ii1; 

	if(n_row->size < 0) return(NULL);

	/* これ以前に値が代入されているかを調べる． 
	 * 代入されていれば，ポインタを返却する． */
	for(ii1=n_row->size-1; ii1>=0; ii1--){
		if(col > n_row->col[ii1]) break;
		if(col == n_row->col[ii1]){
			return( &(n_row->v[ii1]) );
		}
	}    

	/* 代入されていなかったため，要素を1つ増やす． */
	(n_row->size)++;
	n_row->col = (int *)mm_realloc(n_row->col, (n_row->size)*sizeof(int) );
	if(n_row->col == NULL) return(NULL);
	n_row->v = (COMPLEX *)mm_realloc(n_row->v, (n_row->size)*sizeof(COMPLEX) );
	if(n_row->v == NULL) return(NULL);

    /* 値を代入し，列番号と対角項インデックスを格納する． */
	for(ii1=(n_row->size-1); ii1>=1; ii1--){
		if(col > n_row->col[ii1-1]){
			n_row->col[ii1] = col; 
			n_row->v[ii1] = 0.0; 
			if(row == col){
				n_row->diagonal_index = ii1; 
			}else if(n_row->diagonal_index >= ii1){
				n_row->diagonal_index++;
			}
			return( &(n_row->v[ii1]) );
		}
		n_row->col[ii1] = n_row->col[ii1-1];
		n_row->v[ii1] = n_row->v[ii1-1];
	}    
	n_row->col[0] = col;
	n_row->v[0] = 0.0;
	if(row == col){
		n_row->diagonal_index = 0;
	}else if(n_row->diagonal_index >= 0){
		n_row->diagonal_index++;
	}

	return( &(n_row->v[0]) );
}

/* matrix[row][col] がすでに確保してあるなら，その値のポインタを返す．
 * 確保してないなら，新たに確保してそのポインタを返す．
 * この関数使用時は，allocate_crs_cmatrix() を事前に通す．
 * 不備が生じた場合は，NULL が戻り値となる．
 */
COMPLEX *crs_cmatrix_ptr (CRS_CMATRIX *matrix, int row, int col)
{
	return( crs_crow_ptr(&(matrix->row[row]), row, col) );
}


/* CRS行列 dst の実部あるいは虚部に，fact 倍した NONZERO_MATRIX src をたす．
 * 実部と虚部のどちらにたすかは，part で指定する．
 *
 * part
 *   - 1: 実部
 *   - 2: 虚部
 */
RC wadd_crs_dmatrix2cmatrix(int part, double fact, NONZERO_MATRIX src, 
		CRS_CMATRIX *dst)
{
	int ii1;

    for(ii1=0; ii1<src.size; ii1++){
        RC_TRY( wadd_crs_drow2crow(part, fact, ii1, src.row[ii1], 
                    &(dst->row[ii1])) );
    }

	return(NORMAL_RC);
}

/* CRS行 dst の実部あるいは虚部に，fact 倍した NONZERO_ROW src をたす．
 * 実部と虚部のどちらにたすかは，part で指定する．
 *
 * part
 *   - 1: 実部
 *   - 2: 虚部
 */
RC wadd_crs_drow2crow(int part, double fact, int row, NONZERO_ROW src,
		CRS_CROW *dst)
{
	int ii1;
	COMPLEX *ptr;

	/* 実部にたす． */
	if(part == 1){
		for(ii1=0; ii1<src.size; ii1++){
			ptr = crs_crow_ptr(dst, row, src.col[ii1]);
			RC_NULL_CHK( ptr );
			*ptr += fact*src.v[ii1];
		}
	}
	/* 虚部にたす． */
	else if(part == 2){
		for(ii1=0; ii1<src.size; ii1++){
			ptr = crs_crow_ptr(dst, row, src.col[ii1]);
			RC_NULL_CHK( ptr );
			*ptr += fact*src.v[ii1]*COMPLEX_UNIT;
		}
	}else{
		return(ARG_ERROR_RC);
	}

	return(NORMAL_RC);
}

/* ゼロの項を削除する．
 */
RC delete_zero_crs_cmatrix (CRS_CMATRIX *matrix)
{
	int ii1, ii2, ii3;

	for(ii1=0; ii1<(matrix->size); ii1++){
		for(ii2=0; ii2<(matrix->row[ii1].size); ii2++){
			if(cabs(matrix->row[ii1].v[ii2]) >= DBL_MIN) continue;

			for(ii3=0; ii2+ii3<(matrix->row[ii1].size); ii3++){
				matrix->row[ii1].v[ii2+ii3] = matrix->row[ii1].v[ii2+ii3+1];
				matrix->row[ii1].col[ii2+ii3] = matrix->row[ii1].col[ii2+ii3+1];
				if(ii1 == matrix->row[ii1].col[ii2+ii3]){
					matrix->row[ii1].diagonal_index = ii2+ii3;
				}
			}
			(matrix->row[ii1].size)--;
			ii2--;
		}
	}

	return(NORMAL_RC);
}





