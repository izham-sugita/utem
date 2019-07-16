#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "base.h"
#include "mathematics.h"
#include "matrix_db_plus.h"

typedef struct {
	int row;
	int col;
	double v;
} MATRIX;


RC test_matrix_db(MATRIX_DB_INFO *mdbi, int size, int random_flag);
RC test_input(MATRIX_DB_INFO *mdbi, int size, int random_flag, MATRIX **array);
RC test_verify(MATRIX_DB_INFO *mdbi, int size, MATRIX *array);
RC test_delete(MATRIX_DB_INFO *mdbi, int size, MATRIX *array);
RC test_min_col_search(MATRIX_DB_INFO *mdbi, int row_min, int row_max);
RC test_get_row(MATRIX_DB_INFO *mdbi, int row_min, int row_max);
RC test_add_row(MATRIX_DB_INFO *mdbi);

int
main (int argc, char **argv)
{
	MATRIX_DB_INFO mdbi;
	int ii1;
	WRAP64_INT cache[2];
	int input_size;


	if(argc < 2){
		fprintf(stderr, "Usage : %s [scratch directory(ex:/tmp)]\n", argv[0]);
		return(EXIT_FAILURE);
	}

	/* Log Print Level */
	RC_TRY_MAIN( set_log_level(5) );

	/* Memory Manager */
	RC_TRY_MAIN( mm_init((size_t)256*MEGA_BYTE) );

	/* Size Check */
	print_matrix_db_struct_size(stderr);

	/* キャッシュの大きさ */
	/* cache  = size_t * asso * div; */
	cache[0] = sizeof(MATRIX_DB_NODE)*6*5;
	cache[1] = sizeof(MATRIX_DB_LEAF)*8*8;

	input_size = 30;
	RC_TRY_MAIN( log_printf(5, "Input Size : %d\n\n", input_size) );

	for(ii1=0; ii1<1; ii1++){
		/* Open Matrix DB  */
		RC_TRY_MAIN( open_matrix_db(&mdbi, argv[1], cache[0], cache[1], 6, 8));

		RC_TRY_MAIN( log_printf(5, "Input Size : %d\n\n", input_size+ii1) );
		RC_TRY_MAIN( test_matrix_db(&mdbi, input_size+ii1, 1) );
		/*
		RC_TRY_MAIN( print_matrix_db_tree(stderr, &mdbi) );
		RC_TRY_MAIN( print_matrix_db_cache_node_info(stderr, &mdbi.cache_node));
		RC_TRY_MAIN( print_matrix_db_cache_leaf_info(stderr, &mdbi.cache_leaf));
		*/

		/* Close Matrix DB  */
		RC_TRY_MAIN( close_matrix_db(&mdbi) );
	}

	/* Memory Manager */
	RC_TRY( mm_terminate() );

	return(EXIT_SUCCESS);
}


RC
test_matrix_db(MATRIX_DB_INFO *mdbi, int size, int random_flag)
{
	MATRIX *array;

	RC_TRY( test_input(mdbi, size, random_flag, &array) );
	RC_TRY( print_matrix_db_tree(stderr, mdbi) );
	RC_TRY( test_add_row(mdbi) );
	RC_TRY( print_matrix_db_tree(stderr, mdbi) );
	/*
	RC_TRY( print_matrix_db_cache_info(stderr, mdbi) );
	*/

	/*
	RC_TRY( test_get_row(mdbi, 0, 17) )
	RC_TRY( test_min_col_search(mdbi, 0, 15) );
	RC_TRY( test_verify(mdbi, size, array) );
	RC_TRY( test_delete(mdbi, size, array) );
	*/

	RC_TRY( mm_free((void *)array) );
		
	return(NORMAL_RC);
}


RC
test_input (MATRIX_DB_INFO *mdbi, int size, int random_flag, MATRIX **array)
{
	int ii1;
	double v[3][3];
	long rand_s;
	int div = (size>100 ? (size/100) : size);

	RC_NULL_CHK( *array = (MATRIX *)mm_alloc(sizeof(MATRIX)*size) );

	/*
	rand_s = random_seed_time();
	*/
	rand_s = 2;
	log_printf(5, "Random Seed = %ld\n", rand_s);
	init_d_random(rand_s);

	log_printf(5, "Generate input value...\n");
	for(ii1=0; ii1<size; ii1++){
		if(random_flag == 1){ /* ランダム入力 */
			/*
			(*array)[ii1].row = size*d_random();
			*/
			(*array)[ii1].row = ii1%2;
			(*array)[ii1].col = size*d_random();
		}else{ /* 連番入力 */
			(*array)[ii1].row = ii1;
			(*array)[ii1].col = ii1;
		}
		(*array)[ii1].v = 1.0;
		/*
		(*array)[ii1].v = d_random();
		log_printf(5, "[%3d] %3d, %3d, %15.7f\n", ii1,
		           (*array)[ii1].row, (*array)[ii1].col, (*array)[ii1].v);
		*/
		log_printf(5, "%3d, %3d, %15.7f\n",
		           (*array)[ii1].row, (*array)[ii1].col, (*array)[ii1].v);
	}
	log_printf(4, " Done.\n");

	/* Input */
	log_printf(4, "Input Start\n ");
	for(ii1=0; ii1<size; ii1++){
		v[0][0] = v[0][1] = v[0][2] = (*array)[ii1].v;
		v[1][0] = v[1][1] = v[1][2] = (*array)[ii1].v;
		v[2][0] = v[2][1] = v[2][2] = (*array)[ii1].v;
		RC_TRY( add_v_matrix_db((*array)[ii1].row,
		                        (*array)[ii1].col, v, mdbi) );
		/*
		if((ii1+1)%div==0)
			RC_TRY( progress_bar_printf(4, (double)(ii1+1)/size, "") );
		*/
	}
	log_printf(4, "\n");

	return(NORMAL_RC);
}


RC
test_verify (MATRIX_DB_INFO *mdbi, int size, MATRIX *array)
{
	int ii1;
	int found_flag;
	double v[3][3];
	int div = (size>100 ? (size/100) : size);

	/* Verification */
	log_printf(4, "Verification Start\n ");
	for(ii1=0; ii1<size; ii1++){
		RC_TRY( get_v_matrix_db(array[ii1].row, array[ii1].col, v,
		                        &found_flag, mdbi) ); 
		if(found_flag){
			if(
			    (v[0][0] != array[ii1].v )||
			    (v[0][1] != array[ii1].v )||
			    (v[0][2] != array[ii1].v )||
			    (v[1][0] != array[ii1].v )||
			    (v[1][1] != array[ii1].v )||
			    (v[1][2] != array[ii1].v )||
			    (v[2][0] != array[ii1].v )||
			    (v[2][1] != array[ii1].v )||
			    (v[2][2] != array[ii1].v ) ){
				/*
				print_matrix33(stderr, v);
				*/
				log_printf(4, "OverWrite or illegal value!!\n");
				log_printf(4, "array[%d] row = %d, col = %d, v = %15.7f\n ",
				        ii1, array[ii1].row, array[ii1].col, array[ii1].v);
			}
		}else{
			RC_TRY( log_printf(1, "array[%d] row=%d col=%d is missing\n", ii1,
			                   array[ii1].row, array[ii1].col) );
			RC_TRY( print_matrix_db_tree(stderr, mdbi) );
			return(UNKNOWN_ERROR_RC);
		}
		if((ii1+1)%div == 0)
			RC_TRY( progress_bar_printf(4, (double)(ii1+1)/size, "") );
	}
	log_printf(4, "\n");

	return(NORMAL_RC);
}


RC
test_delete (MATRIX_DB_INFO *mdbi, int size, MATRIX *array)
{
	int ii1;
	int div = (size>100 ? (size/100) : size);

	/* Delete */
	log_printf(4, "Delete Start\n ");
	for(ii1=0; ii1<size; ii1++){
		log_printf(6, "[%3d] %3d, %3d, %15.7f\n", ii1, array[ii1].row,
		           array[ii1].col, array[ii1].v);
		RC_TRY( delete_v_matrix_db(array[ii1].row, array[ii1].col, mdbi) );
		if((ii1+1)%div==0)
			RC_TRY( progress_bar_printf(6, (double)(ii1+1)/size, "") );
	}
	log_printf(6, "\n");

	return(NORMAL_RC);
}


RC
test_min_col_search(MATRIX_DB_INFO *mdbi, int row_min, int row_max)
{
	int ii1;
	int min_col;

	for(ii1=row_min; ii1<=row_max; ii1++){
		RC_TRY( min_col_matrix_db(mdbi, ii1, 0, &min_col) );
		log_printf(5, "row[%d] : min_col = %d\n\n", ii1,  min_col);
	}

	return(NORMAL_RC);
}


RC
test_add_row(MATRIX_DB_INFO *mdbi)
{
	int ii1;
	double (*m_array)[3][3];
	int *col_array;
	int size = 20;

	RC_TRY( allocate1I(size, &col_array) );
	RC_TRY( allocate1D33(size, &m_array) );

	for(ii1=0; ii1<size; ii1++){
		col_array[ii1] = ii1+0;
		unit_matrix33(m_array[ii1]);
	}

	RC_TRY( add_row_matrix_db(mdbi, 0, size, m_array, col_array) );

	RC_TRY( free1I(size, &col_array) );
	RC_TRY( free1D33(size, &m_array) );

	return(NORMAL_RC);
}


RC
test_get_row(MATRIX_DB_INFO *mdbi, int row_min, int row_max)
{
	int ii1, ii2;
	int size;
	int *col_array;
	double (*array)[3][3];

	RC_TRY( allocate1I(1000, &col_array) );
	RC_NULL_CHK( array = mm_alloc(sizeof(double[3][3])*1000) );

	for(ii1=row_min; ii1<row_max; ii1++){
		log_printf(5, "\nrow = %d\n", ii1);
		RC_TRY( get_row_matrix_db(mdbi, ii1, 0, INT_MAX, &size, array,
		                          col_array) );
		log_printf(5, "size = %d\n", size);
		for(ii2=0; ii2<size; ii2++){
			log_printf(5, "col_array[%d] = %d\n", ii2, col_array[ii2]);
		}
	}

	RC_TRY( free1I(1000, &col_array) );
	RC_TRY( mm_free(array) );

	return(NORMAL_RC);
}


/* 保存用 */
#if 0 /* search test */
	double (*v)[][3];
	int ii1;
	for(ii1=0; ii1<30; ii1++){
		v = (double(*)[][3])search_v_matrix_db(ii1, ii1, mdb[0]);
		if(v == NULL){
			log_printf(5, "no hit!!\n");
		}else{
			RC_TRY_MAIN( print_matrix33(stderr, *v) );
		}
	}

	v = (double(*)[][3])search_v_matrix_db(0, 1, mdb[0]);
	if(v == NULL){
		log_printf(5, "no hit!!\n");
	}else{
		RC_TRY_MAIN( print_matrix33(stderr, *v) );
	}
#endif
