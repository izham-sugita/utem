#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "base.h"
#include "mathematics.h"
#include "matrix_db.h"

typedef struct {
	int row;
	int col;
	double v;
} MATRIX;

RC test_matrix_db(MATRIX_DB_INFO *mdbi, int size, int random_flag);
RC input(MATRIX_DB_INFO *mdbi, int size, int random_flag, MATRIX **array);
RC verify(MATRIX_DB_INFO *mdbi, int size, MATRIX *array);
RC delete(MATRIX_DB_INFO *mdbi, int size, MATRIX *array);

int
main (int argc, char **argv)
{
	int ii1;
	int cache;
	WRAP64_INT cache_size;
	int input_size;
	MATRIX_DB_INFO mdbi;

	if(argc < 2){
		fprintf(stderr, "Usage : %s [scratch directory(ex:/tmp)]\n", argv[0]);
		return(EXIT_FAILURE);
	}

	/* Log Print Level */
	RC_TRY_MAIN( set_log_level(5) );

	/* Memory Manager */
	RC_TRY_MAIN( mm_init((size_t)512*MEGA_BYTE) );

	/* データモデル チェック */
	print_data_model(stderr);

	/* MATRIX_DB構造体の大きさ */
	print_matrix_db_size(stderr);

	/* キャッシュの大きさ */
	cache = 4;
	cache_size = sizeof(MATRIX_DB)*MATRIX_DB_ASSOCIATE*cache;
	RC_TRY_MAIN( log_printf(5, "Cache Size :"
	                        " [associate][divisor] = [%d][%d] (%d MByte)\n\n",
	                        MATRIX_DB_ASSOCIATE, cache, cache_size/MEGA_BYTE) );
	/*
	input_size = 747;
	*/
	input_size = 100;
	RC_TRY_MAIN( log_printf(5, "Input Size : %d\n\n", input_size) );

	for(ii1=0; ii1<1; ii1++){
		/* Open Matrix DB  */
		RC_TRY_MAIN( open_matrix_db(argv[1], cache_size, &mdbi) );

		/* MatrixDB書き込み 4048 でstratum2 満杯 */
		RC_TRY_MAIN( log_printf(5, "Input Size : %d\n\n", input_size+ii1) );
		RC_TRY_MAIN( test_matrix_db(&mdbi, input_size+ii1, 1) );

		/* Close Matrix DB  */
		RC_TRY_MAIN( close_matrix_db(&mdbi) );
	}

	/*
	RC_TRY_MAIN( print_btree_stratum_info(stderr, &mdbi) );
	*/

	/* Memory Manager */
	RC_TRY( mm_terminate() );

	return(EXIT_SUCCESS);
}


RC
test_matrix_db(MATRIX_DB_INFO *mdbi, int size, int random_flag)
{
	MATRIX *array;

	RC_TRY( input(mdbi, size, random_flag, &array) );

	/*
	RC_TRY( print_matrix_db_cache_info(stderr, mdbi) );
	RC_TRY( print_btree_matrix_db(stderr, mdbi) );
	double v[3][3];
	RC_TRY( print_matrix_db_cache_info(stderr, mdbi) );
	int ii1;
	for(ii1=0; ii1<50; ii1++){
		int col = min_col_matrix_db(mdbi, ii1);
		log_printf(5, "col[%d] = %d\n", ii1,  col);
		if(col < 0) exit(1);
	}
	*/

	/* add_row_matrix_db() TEST */
	{
		int ii1;
		double (*m_array)[3][3];
		int *col_array;
		int size = 20;

		RC_TRY( print_btree_matrix_db(stderr, mdbi) );
		RC_TRY( allocate1I(size, &col_array) );
		RC_TRY( allocate1D33(size, &m_array) );
		for(ii1=0; ii1<size; ii1++){
			col_array[ii1] = ii1+7;
			unit_matrix33(m_array[ii1]);
		}

		RC_TRY( add_row_matrix_db(mdbi, 0, size, m_array, col_array) );

		RC_TRY( print_btree_matrix_db(stderr, mdbi) );
		RC_TRY( print_matrix_db_cache_info(stderr, mdbi) );

	}
#if 0
	{ /* row 削除 */
		int ii1, ii2;
		int array_size;
		double (*m_array)[3][3];
		int *col_array;
		int delete_row;
		delete_row = 10;
		RC_TRY( allocate1I(1000, &col_array) );
		RC_NULL_CHK( m_array = mm_alloc(sizeof(double[3][3])*1000) );
		for(ii1=0; ii1<size; ii1++){
			get_row_matrix_db(mdbi, ii1, &array_size, m_array, col_array);
			log_printf(5, "size = %d\n", array_size);

			for(ii2=0; ii2<array_size; ii2++){
				fprintf(stdout, "row[%d]=%d col[%d]=%d\n", ii2, ii1,
				                                         ii2, col_array[ii2]);
				RC_TRY( delete_v_matrix_db(ii1, col_array[ii2], mdbi) );
			}
		}
		/*
		RC_TRY( print_btree_matrix_db(stderr, mdbi) );
		*/

		/*
		for(ii1=0; ii1<16; ii1++){
			log_printf(5, "\nrow = %d\n", ii1);
			RC_TRY( get_row_col_matrix_db(mdbi, ii1, 0, &array_size, m_array,
			                          col_array) );
			log_printf(5, "array_size = %d\n", array_size);
			for(ii2=0; ii2<array_size; ii2++){
				log_printf(5, "col_array[%d] = %d\n", ii2, col_array[ii2]);
			}
		}
		*/

		RC_TRY( free1I(100, &col_array) );
		RC_TRY( mm_free(m_array) );

	}
#endif

	/*
	RC_TRY( print_btree_matrix_db(stderr, mdbi) );
	{
		int ii1;
		for(ii1=0; ii1<16; ii1++){
			log_printf(5, "\x1b[31mDELETE ROW\x1b[0m = %d\n", ii1);
			RC_TRY( delete_row_range_matrix_db(mdbi, ii1, 0, INT_MAX) );
			//DEBUG_PRINT;
			//RC_TRY( print_btree_matrix_db(stderr, mdbi) );
		}
	}
	RC_TRY( print_btree_matrix_db(stderr, mdbi) );
	RC_TRY( print_btree_matrix_db(stderr, mdbi) );
	RC_TRY( print_matrix_db_cache_info(stderr, mdbi) );

			RC_TRY( print_btree_matrix_db(stderr, mdbi) );
			log_printf(5, "===================================\n");
	RC_TRY( print_stratum_info_matrix_db(stderr, mdbi) );
	RC_TRY( delete_row_range_matrix_db(mdbi, 1, 0, INT_MAX) );
	*/

	/*
	RC_TRY( verify(mdbi, size, array) );
	RC_TRY( delete(mdbi, size, array) );
	RC_TRY( print_matrix_db_cache_info(stderr, mdbi) );
	*/

	RC_TRY( mm_free((void *)array) );
		
	return(NORMAL_RC);
}


RC
input (MATRIX_DB_INFO *mdbi, int size, int random_flag, MATRIX **array)
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

	log_printf(5, "Generate input value...");
	for(ii1=0; ii1<size; ii1++){
		if(random_flag == 1){
			/* ランダム入力 */
			/*
			(*array)[ii1].row = size*d_random();
			*/
			(*array)[ii1].row = ii1%16;
			(*array)[ii1].col = size*d_random();
		}else{
			/* 連番入力 */
			(*array)[ii1].row = ii1;
			(*array)[ii1].col = ii1;
		}
		(*array)[ii1].v = 1.0;
		/*
		(*array)[ii1].v = d_random();
		fprintf(stderr, "[%d] %3d, %3d, %15.7f\n", ii1, (*array)[ii1].row,
		        (*array)[ii1].col, (*array)[ii1].v);
		*/
	}
	log_printf(4, " Done.\n");

	/* Input */
	log_printf(4, "Input Start\n ");
	for(ii1=0; ii1<size; ii1++){
		v[0][0] = v[0][1] = v[0][2] = (*array)[ii1].v;
		v[1][0] = v[1][1] = v[1][2] = (*array)[ii1].v;
		v[2][0] = v[2][1] = v[2][2] = (*array)[ii1].v;
		RC_TRY( put_v_matrix_db((*array)[ii1].row,
		                        (*array)[ii1].col, v, mdbi) );
		if((ii1+1)%div==0)
			RC_TRY( progress_bar_printf(4, (double)(ii1+1)/size, "") );
	}
	log_printf(4, "\n");

#if 0
	/* 検証と削除の邪魔になるので 同一のrow-colを削除する */
	log_printf(5, "Delete Same row-col\n ");
	for(ii1=0; ii1<size; ii1++){
		int ii2;
		for(ii2=ii1+1; ii2<size; ii2++){
			if( (array[ii1].row == array[ii2].row)
			 && (array[ii1].col == array[ii2].col) ){
				array[ii1] = array[size-1];
				size--;
				RC_TRY( log_printf(5, "Same row=%d, col=%d\n\n",
				                   array[ii2].row, array[ii2].col) );
			}
		}
		if((ii1+1)%div==0)
			RC_TRY( progress_bar_printf(5, (double)(ii1+1)/size, "") );
	}
	log_printf(5, "\n");
#endif
	return(NORMAL_RC);
}


RC
verify (MATRIX_DB_INFO *mdbi, int size, MATRIX *array)
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
		}
		if((ii1+1)%div == 0)
			RC_TRY( progress_bar_printf(4, (double)(ii1+1)/size, "") );
	}
	log_printf(4, "\n");

	return(NORMAL_RC);
}


RC
delete (MATRIX_DB_INFO *mdbi, int size, MATRIX *array)
{
	int ii1;
	int div = (size>100 ? (size/100) : size);

	/* Delete */
	log_printf(4, "Delete Start\n ");
	for(ii1=0; ii1<size; ii1++){
		/*
		if(ii1 > 4910)
		fprintf(stderr, "array[%d] row=%3d, col=%3d\n", ii1,
		        array[ii1].row, array[ii1].col);
		*/
		RC_TRY( delete_v_matrix_db(array[ii1].row, array[ii1].col, mdbi) );
#if 0
		/* 一回削除する度に削除していないデータを全て検索する． */
		if(ii1 > 0){ /* バグった寸前から検索したほうがいい */
			int ii2;
			log_printf(5, "Verification Start\n ");
			for(ii2=ii1+1; ii2<size; ii2++){
				RC_TRY( get_v_matrix_db(array[ii2].row, array[ii2].col,
				                        v, mdbi) );
				if(
				    (v[0][0] != array[ii2].v )||
				    (v[0][1] != array[ii2].v )||
				    (v[0][2] != array[ii2].v )||
				    (v[1][0] != array[ii2].v )||
				    (v[1][1] != array[ii2].v )||
				    (v[1][2] != array[ii2].v )||
				    (v[2][0] != array[ii2].v )||
				    (v[2][1] != array[ii2].v )||
				    (v[2][2] != array[ii2].v ) ){
					/*
					print_matrix33(stderr, v);
					*/
					log_printf(5, "Illegal value or MAY BE OverWrite!!\n");
					log_printf(5, " array[%d] row=%d, col=%d, v=%15.7f\n ",
					       ii1, array[ii1].row, array[ii1].col, array[ii1].v);
					exit(0);
				}
				if((ii2+1)%div==0)
					RC_TRY( progress_bar_printf(5, (double)(ii2+1)/size, "") );
			}
			log_printf(5, "\n");
		}
#endif
		if((ii1+1)%div==0)
			RC_TRY( progress_bar_printf(4, (double)(ii1+1)/size, "") );
	}
	log_printf(4, "\n");

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
