#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "base.h"
#include "fem.h"
#include "mathematics.h"
#include "ccs_format.h"
  

static RC modify_ccs_format_dof (int index_xyz, double rest_value, CCS *ccs,
                                 double f_vect[]);

RC
allocate_nonzero_ccs (int n, CCS *matrix)
{
	int ii1;


	matrix->size = matrix->array_size = n;
	matrix->alloc_array_size = n + 100;

	RC_TRY( allocate1I(matrix->size + 1, &(matrix->col_offset)) );
	RC_TRY( allocate1I(matrix->alloc_array_size, &(matrix->row_index)) );
	RC_TRY( allocate1D(matrix->alloc_array_size, &(matrix->array)) );

	/* Initialization */
	/* Diadonal components is initialized by (0.0) */
	for(ii1=0; ii1<matrix->size; ii1++){
#if 0
		matrix->row_index[ii1] = -1;
#endif
		matrix->row_index[ii1] = ii1;
		matrix->col_offset[ii1] = ii1;
		matrix->array[ii1] = 0.0;
	}
	matrix->col_offset[matrix->size] = matrix->size;

	return(NORMAL_RC);
}


RC
clean_nonzero_ccs (CCS *matrix)
{
	/*  Extra memory  is freed */
	if(matrix->alloc_array_size > matrix->array_size){
		matrix->alloc_array_size = matrix->array_size;
		matrix->array = (double *)mm_realloc(matrix->array,
		                          matrix->alloc_array_size*sizeof(double));
		matrix->row_index = (int *)mm_realloc(matrix->row_index,
		                           matrix->alloc_array_size*sizeof(int));
		if(matrix->array == NULL || matrix->row_index == NULL){
			return(NULL_ERROR_RC);
		}
	}

	return(NORMAL_RC);
}


RC
free_nonzero_ccs (CCS *matrix)
{
	RC_TRY( free1I(matrix->size + 1, &(matrix->col_offset)) );
	RC_TRY( free1I(matrix->alloc_array_size, &(matrix->row_index)) );
	RC_TRY( free1I(matrix->alloc_array_size, &(matrix->col_index)) );
	RC_TRY( free1D(matrix->alloc_array_size, &(matrix->array)) );

	matrix->size = matrix->array_size = 0;
	matrix->alloc_array_size = 0;
	matrix->row_index = NULL;
	matrix->col_offset = NULL;
	matrix->array = NULL;

	return(NORMAL_RC);
}


double *
nonzero_ccs_ptr (CCS *matrix, int row, int col)
{
	int ii1;
	int left_index;
	int right_index;
	int mid_index;
	int move_size;
	int *i_ptr = NULL;
	double *d_ptr = NULL;

	if(row > matrix->size || col > matrix->size || row < 0 || col < 0){
		return(NULL);
	}

#if 0
	/* The first access */
	if(matrix->row_index[matrix->col_offset[col]] < 0){
		matrix->row_index[matrix->col_offset[col]] = row;
		return(&(matrix->array[matrix->col_offset[col]]));
	}
#endif

	/* Binary search */
	mid_index = -1;
	left_index = matrix->col_offset[col];
	right_index = matrix->col_offset[col+1] - 1;
	while(left_index <= right_index){
		mid_index = (left_index + right_index) / 2;

		if(matrix->row_index[mid_index] == row){
			return(&(matrix->array[mid_index]));
		}

		if(matrix->row_index[mid_index] < row){
			left_index = mid_index + 1;
		}else{
			right_index = mid_index - 1;
		}
	}
	if(mid_index < 0) return(NULL);
	if( matrix->row_index[mid_index] < row
	&&  mid_index < matrix->col_offset[col+1]){
		mid_index++;
	}

	/* New nonzero value */
	(matrix->array_size)++;
	if(matrix->array_size == matrix->alloc_array_size){
		matrix->alloc_array_size += matrix->size * 10;
		matrix->array = (double *)mm_realloc(matrix->array,
		                          matrix->alloc_array_size*sizeof(double));
		matrix->row_index = (int *)mm_realloc(matrix->row_index,
		                           matrix->alloc_array_size*sizeof(int));
		if(matrix->array == NULL || matrix->row_index == NULL){
			return(NULL);
		}
	}

	for(ii1=(col+1); ii1<=matrix->size; ii1++){
		matrix->col_offset[ii1]++;
	}

	move_size = matrix->array_size - (mid_index + 1);
#if 0
	for(ii1=matrix->array_size-1; ii1>mid_index; ii1--){
		matrix->array[ii1] = matrix->array[ii1-1];
		matrix->row_index[ii1] = matrix->row_index[ii1-1];
	}
#endif
	d_ptr = (double *)memmove(&(matrix->array[mid_index+1]),
	                          &(matrix->array[mid_index]),
	                          move_size*sizeof(double));
	i_ptr = (int *)memmove(&(matrix->row_index[mid_index+1]),
	                       &(matrix->row_index[mid_index]),
	                       move_size*sizeof(int));
	if(d_ptr == NULL || i_ptr == NULL) return(NULL);

	matrix->row_index[mid_index] = row;
	matrix->array[mid_index] = 0.0;
	return(&(matrix->array[mid_index]));
}


int
search_nonzero_ccs_array_index (CCS matrix, int row, int col)
{
	int mid_index;
	int left_index;
	int right_index;
	

	/* Binary search */
	left_index = matrix.col_offset[col];
	right_index = matrix.col_offset[col+1] - 1;
	while(left_index <= right_index){
		mid_index = (left_index + right_index) / 2;

		if(matrix.row_index[mid_index] == row){
			return(mid_index);
		}

		if(matrix.row_index[mid_index] < row){
			left_index = mid_index + 1;
		}else{
			right_index = mid_index - 1;
		}
	}

	return(-1);
}


int
search_nonzero_crs_array_index (CRS matrix, int row, int col)
{
	int mid_index;
	int left_index;
	int right_index;


	/* Binary search */
	left_index = matrix.row_offset[row];
	right_index = matrix.row_offset[row+1] - 1;
	while(left_index <= right_index){
		mid_index = (left_index + right_index) / 2;

		if(matrix.col_index[mid_index] == col){
			return(mid_index);
		}

		if(matrix.col_index[mid_index] < col){
			left_index = mid_index + 1;
		}else{
			right_index = mid_index - 1;
		}
	}

	return(-1);
}


double *
nonzero_matrix1a_ccs_ptr (NONZERO_MATRIX_CCS *matrix, int row, int col)
{
	int ii1;


	if( (0 > matrix->size) || (matrix->size <= row) ) return(NULL);

	for(ii1=matrix->col[col].size-1; ii1>=0; ii1--){
		if(row > matrix->col[col].row[ii1]) break;
		if(row == matrix->col[col].row[ii1]){
			return( &(matrix->col[col].v[ii1]) );
		}
	}

	(matrix->col[col].size)++;
	matrix->col[col].row = (int *)mm_realloc(matrix->col[col].row,
	                              matrix->col[col].size * sizeof(int) );
	if(matrix->col[col].row == NULL) return(NULL);
	matrix->col[col].v = (double *)mm_realloc(matrix->col[col].v,
	                               matrix->col[col].size * sizeof(double) );
	if(matrix->col[col].v == NULL) return(NULL);

	for(ii1=matrix->col[col].size-1; ii1>=1; ii1--){
		if(row > matrix->col[col].row[ii1-1]){
			matrix->col[col].row[ii1] = row;
			matrix->col[col].v[ii1] = 0.0;
			if(col == row){
				matrix->col[col].diagonal_index = ii1;
			}else if(matrix->col[col].diagonal_index >= ii1){
				matrix->col[col].diagonal_index++;
			}
			return( &(matrix->col[col].v[ii1]) );
		}
		matrix->col[col].row[ii1] = matrix->col[col].row[ii1-1];
		matrix->col[col].v[ii1] = matrix->col[col].v[ii1-1];
	}
	matrix->col[col].row[0] = row;
	matrix->col[col].v[0] = 0.0;
	if(col == row){
		matrix->col[col].diagonal_index = 0;
	}else if(matrix->col[col].diagonal_index >= 0){
		matrix->col[col].diagonal_index++;
	}

	return( &(matrix->col[col].v[0]) );
}


RC
allocate_nonzero_matrix1a_ccs (int total_dof, NONZERO_MATRIX_CCS *matrix)
{
	int ii1;


	matrix->size = total_dof;
	matrix->col = (NONZERO_COL *)mm_alloc(matrix->size*sizeof(NONZERO_COL) );
	if(matrix->col == NULL) return(ALLOC_ERROR_RC);

	for(ii1=0; ii1<matrix->size; ii1++){
		matrix->col[ii1].size = 0;
		matrix->col[ii1].diagonal_index = -1;
		matrix->col[ii1].row = NULL;
		matrix->col[ii1].v = NULL;
	}

	return(NORMAL_RC);
}


RC
free_nonzero_matrix1a_ccs (NONZERO_MATRIX_CCS *matrix)
{
	int ii1;


	if(matrix->size <= 0) return(ARG_ERROR_RC);
	if(matrix->col == NULL) return(ARG_ERROR_RC);

	for(ii1=0; ii1<matrix->size; ii1++){
		if(matrix->col[ii1].size > 0){
			if(matrix->col[ii1].row == NULL) return(ARG_ERROR_RC);
			if(matrix->col[ii1].v == NULL) return(ARG_ERROR_RC);
			RC_TRY( mm_free( (void *)matrix->col[ii1].row ) );
			RC_TRY( mm_free( (void *)matrix->col[ii1].v ) );
		}
		matrix->col[ii1].size = 0;
		matrix->col[ii1].diagonal_index = -1;
		matrix->col[ii1].row = NULL;
		matrix->col[ii1].v = NULL;
	}
	RC_TRY( mm_free( (void *)matrix->col ) );
	matrix->size = 0;
	matrix->col = NULL;

	return(NORMAL_RC);
}


RC
copy_nonzero_matrix1a_ccs2ccs_format (NONZERO_MATRIX_CCS *matrix, CCS *ccs)
{
	int ii1, ii2;
	int nz;
	int nz_ptr;


	if(matrix->size <= 0) return(ARG_ERROR_RC);
	if(matrix->col == NULL) return(ARG_ERROR_RC);

	ccs->size = matrix->size;
	/* Total number of nonzero elements */
	nz = 0;
	for(ii1=0; ii1<matrix->size; ii1++){
		nz += matrix->col[ii1].size;
	}
	ccs->array_size = nz;
	ccs->alloc_array_size = nz;

	RC_TRY( allocate1I(ccs->size + 1, &(ccs->col_offset)) );
	RC_TRY( allocate1I(nz, &(ccs->row_index)) );
	RC_TRY( allocate1I(nz, &(ccs->col_index)) );
	RC_TRY( allocate1D(nz, &(ccs->array)) );

	nz_ptr = 0;
	for(ii1=0; ii1<matrix->size; ii1++){
		NONZERO_COL col_elem = matrix->col[ii1];

		ccs->col_offset[ii1] = nz_ptr;
		for(ii2=0; ii2<col_elem.size; ii2++){
			ccs->row_index[nz_ptr] = col_elem.row[ii2];
			ccs->col_index[nz_ptr] = ii1;
			ccs->array[nz_ptr] = col_elem.v[ii2];
			nz_ptr++;
		}
	}
	ccs->col_offset[ccs->size] = nz;

	return(NORMAL_RC);
}


RC
copy_nonzero_matrix1a_crs2crs_format (NONZERO_MATRIX *matrix, CRS *crs)
{
	int ii1, ii2;
	int nz;
	int nz_ptr;


	if(matrix->size <= 0) return(ARG_ERROR_RC);
	if(matrix->row == NULL) return(ARG_ERROR_RC);

	crs->size = matrix->size;
	/* Total number of nonzero elements */
	nz = 0;
	for(ii1=0; ii1<matrix->size; ii1++){
		nz += matrix->row[ii1].size;
	}
	crs->array_size = nz;
	crs->alloc_array_size = nz;

	RC_TRY( allocate1I(crs->size + 1, &(crs->row_offset)) );
	RC_TRY( allocate1I(nz, &(crs->col_index)) );
	RC_TRY( allocate1D(nz, &(crs->array)) );

	nz_ptr = 0;
	for(ii1=0; ii1<matrix->size; ii1++){
		NONZERO_ROW row_elem = matrix->row[ii1];

		crs->row_offset[ii1] = nz_ptr;
		for(ii2=0; ii2<row_elem.size; ii2++){
			crs->col_index[nz_ptr] = row_elem.col[ii2];
			crs->array[nz_ptr] = row_elem.v[ii2];
			nz_ptr++;
		}
	}
	crs->row_offset[crs->size] = nz_ptr;

	return(NORMAL_RC);
}


RC
free_nonzero_crs (CRS *matrix)
{
	RC_TRY( free1I(matrix->size + 1, &(matrix->row_offset)) );
	RC_TRY( free1I(matrix->alloc_array_size, &(matrix->col_index)) );
	RC_TRY( free1D(matrix->alloc_array_size, &(matrix->array)) );

	matrix->size = matrix->array_size = 0;
	matrix->alloc_array_size = 0;
	matrix->col_index = NULL;
	matrix->row_offset = NULL;
	matrix->array = NULL;

	return(NORMAL_RC);
}


RC
modify_nonzero1a_crs (int dim, NODE_ARRAY node, BC_ARRAY rest,
                      NONZERO_MATRIX *matrix, double f_vect[])
{
	int ii1;

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	for(ii1=0; ii1<rest.size; ii1++){
		int index;
		if(rest.array[ii1].node < 0) continue;

		RC_NEG_CHK( index = node_renum_index(node, rest.array[ii1].node) );

		if(rest.array[ii1].v_type.x == BC_FIX){
			RC_TRY( modify_nonzero1a_crs_dof(dim*index + 0, rest.array[ii1].v.x,
			                                 matrix, f_vect) );
		}
		if(rest.array[ii1].v_type.y == BC_FIX){
			RC_TRY( modify_nonzero1a_crs_dof(dim*index + 1, rest.array[ii1].v.y,
			                                 matrix, f_vect) );
		}
		if(dim == 3){
			if(rest.array[ii1].v_type.z == BC_FIX){
				RC_TRY( modify_nonzero1a_crs_dof(dim*index + 2,
				                                 rest.array[ii1].v.z,
				                                 matrix, f_vect) );
			}
		}
	}

	return(NORMAL_RC);
}


RC
modify_nonzero1a_crs_dof (int index_xyz, double rest_value,
                          NONZERO_MATRIX *matrix, double f_vect[])
{
	int ii1, ii2;


	for(ii1=0; ii1<matrix->size; ii1++){
		for(ii2=0; ii2<matrix->row[ii1].size; ii2++){
			if(matrix->row[ii1].col[ii2] > index_xyz) break;
			if(matrix->row[ii1].col[ii2] == index_xyz){
				f_vect[ii1] -= matrix->row[ii1].v[ii2] * rest_value;
				matrix->row[ii1].v[ii2] = 0.0;
			}
		}
	}

	f_vect[index_xyz] = rest_value;
	for(ii1=0; ii1<matrix->row[index_xyz].size; ii1++){
		matrix->row[index_xyz].v[ii1] = 0.0;
	}
	matrix->row[index_xyz].v[matrix->row[index_xyz].diagonal_index] = 1.0;


	return(NORMAL_RC);
}
                          
RC
modify_ccs_format (int dim, NODE_ARRAY node, BC_ARRAY rest, CCS *ccs,
                   double f_vect[])
{
	int ii1;
	int index;

	for(ii1=0; ii1<rest.size; ii1++){
		if(rest.array[ii1].node < 0) continue;
		index = node_renum_index(node, rest.array[ii1].node);
		RC_NEG_CHK(index);

		if(rest.array[ii1].v_type.x == BC_FIX){
			RC_TRY( modify_ccs_format_dof(dim*index + 0, rest.array[ii1].v.x,
                                          ccs, f_vect) );
		}
		if(rest.array[ii1].v_type.y == BC_FIX){
			RC_TRY( modify_ccs_format_dof(dim*index + 1, rest.array[ii1].v.y,
                                          ccs, f_vect) );
		}
		if(dim == 3){
			if(rest.array[ii1].v_type.z == BC_FIX){
				RC_TRY( modify_ccs_format_dof(dim*index + 2,
                                              rest.array[ii1].v.z, ccs,
                                              f_vect) );
			}
		}
	}

	return(NORMAL_RC);
}

static RC
modify_ccs_format_dof (int index_xyz, double rest_value, CCS *ccs,
                       double f_vect[])
{
	int ii1;
	int index;


	for(ii1=ccs->col_offset[index_xyz];
        ii1<ccs->col_offset[index_xyz+1]; ii1++){
		f_vect[ccs->row_index[ii1]] -= ccs->array[ii1] * rest_value;
		ccs->array[ii1] = 0.0;
	}
	f_vect[index_xyz] = rest_value;

	for(ii1=0; ii1<ccs->size; ii1++){
		if(ccs->row_index[ccs->col_offset[ii1]] > index_xyz
        || ccs->row_index[ccs->col_offset[ii1+1]-1] < index_xyz) continue;

		index = search_nonzero_ccs_array_index(*ccs, index_xyz, ii1);
		if(index < 0) continue;

		ccs->array[index] = 0.0;
	}

	index = search_nonzero_ccs_array_index(*ccs, index_xyz, index_xyz);
	if(index < 0) return(SEEK_ERROR_RC);
	ccs->array[index] = 1.0;

	return(NORMAL_RC);
}

