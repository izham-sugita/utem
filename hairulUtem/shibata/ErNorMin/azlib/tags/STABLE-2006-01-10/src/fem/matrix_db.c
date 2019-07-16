/*********************************************************************
 * mtrix_db.c
 *
 * Copyright (C) 2005 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA> <Kenzen TAKEUCHI>
 *
 * Thanks to following contributors.
 *  <Taiki AOYAMA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: matrix_db.c 722 2006-01-05 10:37:07Z sasaoka $ */

#include <stdio.h>
#include <math.h>
#ifndef WIN32
#include <unistd.h>
#endif  /* WIN32 */
#include <string.h>
#include <limits.h>
#include "wrapper_64.h"
#include "scratch_io.h"
#include "matrix_db.h"
#include "mathematics.h"
#include "rc.h"
#include "log_printf.h"
#include "memory_manager.h"
#include "macros.h"

typedef struct{
	int size;
	int alloc_size;
	int index;
	int *row;
	int *col;
	double (*v)[3][3];
} MATRIX_DB_ARY;

#define MAX_ACCESS_COUNT   (INT_MAX/16 + 1)
#define CACHE_ROW(offset, divisor) ( ((offset)/sizeof(MATRIX_DB))%(divisor) )

/* MATRIX_DB 構造体を初期化する時に使用 */
static MATRIX_DB Initialized_MDB;    /* 予め初期化しておくこと */
#define INIT_MATRIX_DB(mdb)  {mdb = Initialized_MDB;}

static RC min_col_search_parent_matrix_db(MATRIX_DB_INFO *mdbi,
                                          int row, int col, int *ret);
static RC add_row_clean_buf(MATRIX_DB_ARY *buf);
static RC add_row_sub(MATRIX_DB_INFO *mdbi, MATRIX_DB *mdb,
                      MATRIX_DB_ARY *data, MATRIX_DB_ARY *buf, int row,
                      int col, int which, int next_row, int next_col);
static RC add_row_child(MATRIX_DB_INFO *mdbi, MATRIX_DB *mdb, int index,
                        MATRIX_DB_ARY *data, MATRIX_DB_ARY *buf, int row,
                        int col, int which, int next_row, int next_col);
static RC add_row_leaf(MATRIX_DB_INFO *mdbi, MATRIX_DB *mdb, int index,
                       MATRIX_DB_ARY *data, MATRIX_DB_ARY *buf, int row,
                       int col, int which, int next_row, int next_col);
static RC add_row_min_row_col(MATRIX_DB_ARY *data, MATRIX_DB_ARY *buf,
                              int *row, int *col, int *which);
static RC get_row_range_sub(MATRIX_DB_INFO *mdbi, MATRIX_DB *mdb, int row,
                            int col_min, int col_max, int *size,
                            double array[][3][3], int col_array[]);
static RC get_row_range_child(MATRIX_DB_INFO *mdbi, MATRIX_DB *mdb, int index,
                              int row, int col_min, int col_max, int *size,
                              double array[][3][3], int col_array[]);
static RC delete_row_range_sub(MATRIX_DB_INFO *mdbi, MATRIX_DB *mdb, int row,
                               int col_min, int col_max);
static RC delete_row_range_child(MATRIX_DB_INFO *mdbi, MATRIX_DB *mdb,
                                 int index, int row, int col_min, int col_max);
static RC delete_row_range_tree(MATRIX_DB_INFO *mdbi, int p_file,
                                WRAP64_INT p_offset);
static RC delete_row_range_leaf(MATRIX_DB_INFO *mdbi, MATRIX_DB *mdb,
                                int index, int row, int col_min, int col_max);
static RC input_v_btree_matrix_db(int row, int col, double v[][3], int index,
                                  MATRIX_DB *mdb, MATRIX_DB_INFO *mdbi);
static RC search_parent_matrix_db(MATRIX_DB_INFO *mdbi, MATRIX_DB *mdb,
                                  int *index);
static RC split_btree_matrix_db_sub(MATRIX_DB *child, MATRIX_DB *parent,
                                    MATRIX_DB *new_child);
static MATRIX_DB *search_cache_ptr(MATRIX_DB_INFO *mdbi, int db_file,
                                   WRAP64_INT db_offset, int inc_access);
static MATRIX_DB *unused_cache_ptr(MATRIX_DB_INFO *mdbi, int crow);
static int compar_row_col(int row1, int col1, int row2, int col2);


/* cache_size はキャッシュとして利用可能なバイト数で最低でも */
/* sizeof(MATRIX_DB)*MATRIX_DB_ASSOCIATEバイト分は必要       */
RC
open_matrix_db (const char scratch_dir[], WRAP64_INT size,
                MATRIX_DB_INFO *mdbi)
{
	int ii1, ii2;
	WRAP64_INT mdb_asso_size = sizeof(MATRIX_DB)*MATRIX_DB_ASSOCIATE;

	RC_NULL_CHK( mdbi );
	if(size < mdb_asso_size) return(ARG_ERROR_RC);

	RC_TRY( init_matrix_db(&Initialized_MDB) );

	for(ii1=0; ii1<MATRIX_DB_SF_SIZE; ii1++){
		mdbi->scratch_array[ii1] = open_scratch_file(scratch_dir);
		if(mdbi->scratch_array[ii1] == NULL) return(OPEN_ERROR_RC);
		mdbi->last_p_offset[ii1] = 0;
		mdbi->blank_p_offset[ii1] = -1;
		mdbi->current_p_offset[ii1] = 0;
	}

	/* キャッシュ */
	mdbi->cache_divisor = (int)(size/mdb_asso_size);
	mdbi->cache = (MATRIX_DB (*)[MATRIX_DB_ASSOCIATE])mm_alloc(
	                            (size_t)(mdb_asso_size *(mdbi->cache_divisor)));
	RC_NULL_CHK( mdbi->cache );
	mdbi->access_count = (int (*)[MATRIX_DB_ASSOCIATE])mm_alloc(
	           (size_t)(sizeof(int)*MATRIX_DB_ASSOCIATE*(mdbi->cache_divisor)));
	RC_NULL_CHK( mdbi->access_count );
	for(ii1=0; ii1<(mdbi->cache_divisor); ii1++){
		for(ii2=0; ii2<MATRIX_DB_ASSOCIATE; ii2++){
			INIT_MATRIX_DB(mdbi->cache[ii1][ii2]);
			mdbi->access_count[ii1][ii2] = -1;
		}
	}

	/* 経路情報 */
	mdbi->route_size = 0;
	mdbi->route_count = 0;
	RC_TRY( allocate1I(MATRIX_DB_ROUTE_SIZE, &(mdbi->route_file)) );
	mdbi->route_offset = (WRAP64_INT *)mm_alloc(
	                      (size_t)(MATRIX_DB_ROUTE_SIZE*sizeof(WRAP64_INT)) );
	RC_NULL_CHK( mdbi->route_offset );

	/* Root節点を確定 */
	/* RootNodeを消さないためにreadする */
	RC_TRY( new_matrix_db_start(&mdbi->root, mdbi) );
	RC_TRY( set_root_matrix_db(mdbi->root, mdbi) );
	RC_TRY( read_matrix_db_start(mdbi->root->p_file_self,
	                             mdbi->root->p_offset_self,
	                             &mdbi->root, mdbi) );
	RC_TRY( new_matrix_db_end(mdbi->root) );

	return(NORMAL_RC);
}


RC
close_matrix_db (MATRIX_DB_INFO *mdbi)
{
	int ii1;

	RC_NULL_CHK( mdbi );
	RC_TRY( read_matrix_db_end(mdbi->root) );

	RC_TRY( mm_free(mdbi->cache) );
	RC_TRY( mm_free(mdbi->access_count) );
	RC_TRY( mm_free(mdbi->route_file) );
	RC_TRY( mm_free(mdbi->route_offset) );

	for(ii1=0; ii1<MATRIX_DB_SF_SIZE; ii1++){
		RC_TRY( close_scratch_file(mdbi->scratch_array[ii1]) );
	}

	return(NORMAL_RC);
}


RC
get_v_matrix_db (int row, int col, double v[][3], int *found_flag,
                 MATRIX_DB_INFO *mdbi)
{
	int index;
	int p_file;
	WRAP64_INT p_offset;
	MATRIX_DB *mdb;

	RC_NULL_CHK( v );
	RC_NULL_CHK( mdbi );
	RC_NULL_CHK( found_flag );
	if((row < 0) || (col < 0)) return(ARG_ERROR_RC);

	/* 挿入すべき Leaf Nodeを Rootから探索 */
	mdbi->route_size = 0;
	RC_TRY( search_v_btree_matrix_db(mdbi, mdbi->root, row, col, found_flag,
	                                 &index, &p_file, &p_offset) );
	if(p_file <= 0) return(UNKNOWN_ERROR_RC);
	if(p_file > MATRIX_DB_SF_SIZE) return(UNKNOWN_ERROR_RC);
	if(p_offset < 0)  return(UNKNOWN_ERROR_RC);

	if(*found_flag){
		RC_TRY( read_matrix_db_start(p_file, p_offset, &mdb, mdbi) );
		if( (index < 0)||(index >= mdb->size) ) return(UNKNOWN_ERROR_RC);
		RC_TRY( copy_matrix33(mdb->v[index], v) );
		RC_TRY( read_matrix_db_end(mdb) );
	}

	return(NORMAL_RC);
}


RC
put_v_matrix_db (int row, int col, double v[][3], MATRIX_DB_INFO *mdbi)
{
	int index, found_flag;
	int p_file;
	WRAP64_INT p_offset;
	MATRIX_DB *mdb;

	if( (row < 0) || (col < 0) ) return(ARG_ERROR_RC);
	RC_NULL_CHK( v );
	RC_NULL_CHK( mdbi );

	mdbi->route_size = 0;
	/* 挿入すべき Leaf Nodeを Rootから探索 */
	RC_TRY( search_v_btree_matrix_db(mdbi, mdbi->root, row, col, &found_flag,
	                                 &index, &p_file, &p_offset) );
	if(p_file <= 0) return(UNKNOWN_ERROR_RC);
	if(p_file > MATRIX_DB_SF_SIZE) return(UNKNOWN_ERROR_RC);
	if(p_offset < 0)  return(UNKNOWN_ERROR_RC);

	RC_TRY( write_matrix_db_start(p_file, p_offset, &mdb, mdbi) );
	if( (index < 0)||(index > mdb->size) ) return(UNKNOWN_ERROR_RC);
	if(found_flag){
		/* 上書きして終了 */
		RC_TRY( log_printf(6, "OverWrite row=%d, col=%d\n", row, col) );
		RC_TRY( set_v_matrix_db(row, col, v, index, mdb) );
	}else{
		RC_TRY( input_v_btree_matrix_db(row, col, v, index, mdb, mdbi) );
	}
	RC_TRY( write_matrix_db_end(mdb) );

	return(NORMAL_RC);
}


/* row-colの vがあれば足して上書きする．無ければ新規追加する */
RC
add_v_matrix_db (int row, int col, double v[][3], MATRIX_DB_INFO *mdbi)
{
	int index, found_flag;
	MATRIX_DB *mdb;
	double tmp[3][3];
	int p_file;
	WRAP64_INT p_offset;

	if( (row < 0) || (col < 0) ) return(ARG_ERROR_RC);
	RC_NULL_CHK( v );
	RC_NULL_CHK( mdbi );

	mdbi->route_size = 0;
	/* 挿入すべき Leaf Nodeを Rootから探索 */
	RC_TRY( search_v_btree_matrix_db(mdbi, mdbi->root, row, col, &found_flag,
	                                 &index, &p_file, &p_offset) );
	if(p_file <= 0) return(UNKNOWN_ERROR_RC);
	if(p_file > MATRIX_DB_SF_SIZE) return(UNKNOWN_ERROR_RC);
	if(p_offset < 0)  return(UNKNOWN_ERROR_RC);

	RC_TRY( write_matrix_db_start(p_file, p_offset, &mdb, mdbi) );
	if(found_flag){
		if( (index < 0)||(index >= mdb->size) ) return(UNKNOWN_ERROR_RC);
		RC_TRY( add_matrix33a(v, mdb->v[index], tmp) );
		RC_TRY( set_v_matrix_db(row, col, tmp, index, mdb) );
	}else{
		if( (index < 0)||(index > mdb->size) ) return(UNKNOWN_ERROR_RC);
		RC_TRY( input_v_btree_matrix_db(row, col, v, index, mdb, mdbi) );
	}
	RC_TRY( write_matrix_db_end(mdb) );

	return(NORMAL_RC);
}


/* row-colの vがあれば掛けて上書きする．無ければ初期化したvが新規追加される */
RC
mul_v_matrix_db (int row, int col, double v[][3], MATRIX_DB_INFO *mdbi)
{
	int index, found_flag;
	MATRIX_DB *mdb;
	double tmp[3][3];
	int p_file;
	WRAP64_INT p_offset;

	if( (row < 0) || (col < 0) ) return(ARG_ERROR_RC);
	RC_NULL_CHK( v );
	RC_NULL_CHK( mdbi );

	mdbi->route_size = 0;
	/* 挿入すべき Leaf Nodeを Rootから探索 */
	RC_TRY( search_v_btree_matrix_db(mdbi, mdbi->root, row, col, &found_flag,
	                                 &index, &p_file, &p_offset) );
	if(p_file <= 0) return(UNKNOWN_ERROR_RC);
	if(p_file > MATRIX_DB_SF_SIZE) return(UNKNOWN_ERROR_RC);
	if(p_offset < 0)  return(UNKNOWN_ERROR_RC);

	RC_TRY( write_matrix_db_start(p_file, p_offset, &mdb, mdbi) );
	if(found_flag){
		if( (index < 0)||(index >= mdb->size) ) return(UNKNOWN_ERROR_RC);
		RC_TRY( mul_matrix33(v, mdb->v[index], tmp) );
		RC_TRY( set_v_matrix_db(row, col, tmp, index, mdb) );
	}else{
		if( (index < 0)||(index > mdb->size) ) return(UNKNOWN_ERROR_RC);
		init_matrix33(tmp);
		RC_TRY( input_v_btree_matrix_db(row, col, tmp, index, mdb, mdbi) );
	}
	RC_TRY( write_matrix_db_end(mdb) );

	return(NORMAL_RC);
}


RC
delete_v_matrix_db (int row, int col, MATRIX_DB_INFO *mdbi)
{
	RC_NULL_CHK( mdbi );
	if( (row < 0) || (col < 0) ) return(ARG_ERROR_RC);
	/* 該当するrow-colが無い時は，エラーを返さない */
	RC_TRY( delete_v_btree_matrix_db(row, col, mdbi) );
	return(NORMAL_RC);
}


/* row 行目で最も小さい列(一番左側)の数を返却 */
int
min_col_matrix_db (MATRIX_DB_INFO *mdbi, int row)
{
	return( min_col_sub_matrix_db(mdbi, row, 0) );
}


/* row 行目の col以上で最小の列の数を返却 */
int
min_col_sub_matrix_db (MATRIX_DB_INFO *mdbi, int row, int col)
{
	RC rc;
	int found_flag;
	int index;
	MATRIX_DB *mdb;
	int ret;
	int p_file;
	WRAP64_INT p_offset;

	if( (mdbi == NULL) || (row < 0) || (col < 0) ){
		rc_error_print(__FILE__, __LINE__, "min_col_sub_matrix_db()",
		               ARG_ERROR_RC);
		return(-1);
	}

	/* (row, col) を探す．結果(row, min_col)のindexが得られる */
	mdbi->route_size = 0;
	rc = search_v_btree_matrix_db(mdbi, mdbi->root, row, col, &found_flag,
	                              &index, &p_file, &p_offset);
	if(rc == NORMAL_RC){
		if( (p_file <= 0) || (p_file > MATRIX_DB_SF_SIZE) || (p_offset < 0) ){
			rc_error_print(__FILE__, __LINE__, "min_col_sub_matrix_db()",
			               UNKNOWN_ERROR_RC);
			return(-2);
		}
		rc = read_matrix_db_start(p_file, p_offset, &mdb, mdbi);
		if(rc != NORMAL_RC){
			rc_error_print(__FILE__, __LINE__, "min_col_sub_matrix_db()", rc);
			return(-2);
		}
		if( (index < 0) || (index > mdb->size) ){
			rc_error_print(__FILE__, __LINE__, "min_col_sub_matrix_db()",
			               UNKNOWN_ERROR_RC);
			return(-2);
		}
		if(found_flag){
			ret = mdb->col[index];
		}else{
			if(index == mdb->size){
				/* 親にcolの最小があるかもしれない */
				rc = min_col_search_parent_matrix_db(mdbi, row, col, &ret);
				if(rc != NORMAL_RC){
					rc_error_print(__FILE__, __LINE__,
					               "min_col_sub_matrix_db()", rc);
					ret = -2;
				}
			}else if(mdb->row[index] != row){
				/* 親にcolの最小があるかもしれない */
				rc = min_col_search_parent_matrix_db(mdbi, row, col, &ret);
				if(rc != NORMAL_RC){
					rc_error_print(__FILE__, __LINE__,
					               "min_col_sub_matrix_db()", rc);
					ret = -2;
				}
			}else{ /* mdb->row[index] == row */
				ret = mdb->col[index];
			}
		}
	}else{
		rc_error_print(__FILE__, __LINE__, "min_col_sub_matrix_db()", rc);
		ret = -2;
	}

	rc = read_matrix_db_end(mdb);
	if(rc != NORMAL_RC){
		rc_error_print(__FILE__, __LINE__, "min_col_sub_matrix_db()", rc);
		ret = -2;
	}

	return(ret);
}


static RC
min_col_search_parent_matrix_db (MATRIX_DB_INFO *mdbi, int row, int col,
                                 int *ret)
{
	int ii1;
	int p_file;
	WRAP64_INT p_offset;
	MATRIX_DB *parent;
	int found_flag;
	int index;

	RC_NULL_CHK(mdbi);
	RC_NULL_CHK(ret);
	if( (row < 0) || (col < 0) ) return(ARG_ERROR_RC);

	*ret = -1;
	for(ii1=mdbi->route_size-2; ii1>=0; ii1--){
		RC_TRY( read_matrix_db_start(mdbi->route_file[ii1],
		                             mdbi->route_offset[ii1], &parent, mdbi) );
		RC_TRY( search_matrix_db(row, col, parent, &found_flag,
		                         &index, &p_file, &p_offset) );
		if( (index < 0)||(index > parent->size) ) return(UNKNOWN_ERROR_RC);
		if(index == parent->size){
			RC_TRY( read_matrix_db_end(parent) );
			continue;
		}else if(parent->row[index] != row){
			RC_TRY( read_matrix_db_end(parent) );
			continue;
		}else{ /* parent->row[index] == row */
			*ret = parent->col[index];
		}
		RC_TRY( read_matrix_db_end(parent) );
		if(*ret >= 0) return(NORMAL_RC);
	}

	return(NORMAL_RC);
}


/* row行目のデータを連続で addする．      */
/* col_arrayはソートされている必要がある．*/
RC
add_row_matrix_db (MATRIX_DB_INFO *mdbi, int row, int size,
                   double array[][3][3], int col_array[])
{
	MATRIX_DB_ARY buf;
	MATRIX_DB_ARY data;
	RC rc;
	int col, which;

	RC_NULL_CHK(mdbi);
	RC_NULL_CHK(array);
	RC_NULL_CHK(col_array);
	if( (row < 0)||(size < 0) ) return(ARG_ERROR_RC);

	if(size == 0) return(NORMAL_RC);

	/* バッファ */
	buf.size = 0;
	buf.alloc_size = size+MATRIX_DB_SIZE; /* sizeでもいいかもしれない */
	buf.index = 0;
	RC_TRY( allocate1I(buf.alloc_size, &(buf.row)) );
	RC_TRY( allocate1I(buf.alloc_size, &(buf.col)) );
	RC_TRY( allocate1D33(buf.alloc_size, &(buf.v)) );

	/* 追加するデータ */
	RC_TRY( allocate1I(1, &(data.row)) );
	data.size = size;
	data.alloc_size = size;
	data.index = 0;
	data.row[0] = row; /* 全部同じなので 1つでいい */
	data.col = col_array;
	data.v = array;

	while(1){ 
		/* 入力する dataとbufの row,colの最小値を探す */
		RC_TRY( add_row_min_row_col(&data, &buf, &row, &col, &which) );
		if(row < 0){
			/* 入力するデータが無いので終了 */
			return(NORMAL_RC);
		}
		rc = add_row_sub(mdbi, mdbi->root, &data, &buf, row, col, which,
		                 INT_MAX, INT_MAX);
		if(rc == NORMAL_RC){
			break;
		}
		if(rc != SPECIAL_RC) return(rc);
	}

	RC_TRY( free1I(buf.alloc_size, &(buf.row)) );
	RC_TRY( free1I(buf.alloc_size, &(buf.col)) );
	RC_TRY( free1D33(buf.alloc_size, &(buf.v)) );

	RC_TRY( free1I(1, &(data.row)) );

	return(NORMAL_RC);
}


/* 左に詰める */
static RC
add_row_clean_buf (MATRIX_DB_ARY *buf)
{
	int ii1;
	int pos;

	RC_NULL_CHK(buf);
	if(buf->index == 0) return(NORMAL_RC);
	if(buf->index == buf->size){
		buf->index = 0;
		buf->size  = 0;
		return(NORMAL_RC);
	}

	pos = 0;
	for(ii1=buf->index; ii1<buf->size; ii1++){
		buf->row[pos] = buf->row[ii1];
		buf->col[pos] = buf->col[ii1];
		RC_TRY( copy_matrix33(buf->v[ii1], buf->v[pos]) );
		pos++;
	}
	buf->index = 0;
	buf->size  = pos;

	return(NORMAL_RC);
}


static RC
add_row_sub (MATRIX_DB_INFO *mdbi, MATRIX_DB *mdb, MATRIX_DB_ARY *data,
             MATRIX_DB_ARY *buf, int row, int col, int which, int next_row,
             int next_col)
{
	RC rc;
	int p_file;
	WRAP64_INT p_offset;
	int found_flag;
	int index;

	RC_NULL_CHK(mdbi);
	RC_NULL_CHK(mdb);
	RC_NULL_CHK(buf);
	RC_NULL_CHK(data);
	if( (row < 0) || (col < 0) ) return(ARG_ERROR_RC);
	if( (next_row < 0) || (next_col < 0) ) return(ARG_ERROR_RC);
	if( (which < -1) || (which > 1) ) return(ARG_ERROR_RC);

	RC_TRY( search_matrix_db(row, col, mdb, &found_flag, &index,
	                         &p_file, &p_offset) );
	if( (index < 0)||(index > mdb->size ) ) return(UNKNOWN_ERROR_RC);

	if( is_node_matrix_db(mdb) ){
		if(index == mdb->size){
			/* next_row, next_colは引き継ぎ */
			rc = add_row_child(mdbi, mdb, index, data, buf, row, col,
			                   which, next_row, next_col);
			if(rc != NORMAL_RC) return(rc);
		}else if(mdb->row[index] > row){
			rc = add_row_child(mdbi, mdb, index, data, buf, row, col,
			                   which, mdb->row[index], mdb->col[index]);
			if(rc != NORMAL_RC) return(rc);
		}else if(mdb->row[index] < row){
			if(index == mdb->size-1){
				rc = add_row_child(mdbi, mdb, index+1, data, buf, row, col,
				                   which, next_row, next_col);
			}else{
				rc = add_row_child(mdbi, mdb, index+1, data, buf, row, col,
				                   which, mdb->row[index+1],
				                   mdb->col[index+1]);
			}
			if(rc != NORMAL_RC) return(rc);
		}else{ /* mdb->row[index] == row */
			if(mdb->col[index] > col){
				rc = add_row_child(mdbi, mdb, index, data, buf, row, col,
				                   which, mdb->row[index], mdb->col[index]);
				if(rc != NORMAL_RC) return(rc);
			}else if(mdb->col[index] < col){
				if(index == mdb->size-1){
					rc = add_row_child(mdbi, mdb, index+1, data, buf, row, col,
					                   which, next_row, next_col);
				}else{
					rc = add_row_child(mdbi, mdb, index+1, data, buf, row, col,
					                   which, mdb->row[index+1],
					                   mdb->col[index+1]);
				}
				if(rc != NORMAL_RC) return(rc);
			}else{ /* mdb->col[index] == col */
				RC_TRY( write_matrix_db_start(mdb->p_file_self,
				                              mdb->p_offset_self, &mdb,
				                              mdbi) );
				if(which == -1){
					add_matrix33(mdb->v[index], data->v[data->index]);
					data->index++;
				}else if(which == 0){
					double v_tmp[3][3];
					RC_TRY( add_matrix33a(data->v[data->index],
					                      buf->v[buf->index], v_tmp) );
					add_matrix33(mdb->v[index], v_tmp);
					data->index++;
					buf->index++;
				}else{ /* (which = 1) */
					add_matrix33(mdb->v[index], buf->v[buf->index]);
					buf->index++;
				}
				RC_TRY( write_matrix_db_end(mdb) );
				RC_TRY( add_row_min_row_col(data, buf, &row, &col, &which) );
				if(row < 0){
					/* 被入力データは無い */
					return(NORMAL_RC);
				}else{
					/* まだ残っている */
					if(index == mdb->size-1){
						rc = add_row_child(mdbi, mdb, index+1, data, buf, row,
						                   col, which, next_row, next_col);
					}else{
						rc = add_row_child(mdbi, mdb, index+1, data, buf, row,
						                   col, which, mdb->row[index+1],
						                   mdb->col[index+1]);
					}
					if(rc != NORMAL_RC) return(rc);
				}
			}
		}
	}else{ /* Leaf */
		RC_TRY( write_matrix_db_start(mdb->p_file_self,
		                              mdb->p_offset_self, &mdb, mdbi) );
		rc = add_row_leaf(mdbi, mdb, index, data, buf, row, col, which,
		                  next_row, next_col);
		if(rc == SPECIAL_RC){
			if(mdb->size == MATRIX_DB_SIZE){
				MATRIX_DB *new;
				RC_TRY( new_matrix_db_start(&new, mdbi) );
				RC_TRY( split_btree_matrix_db(mdb, new, mdbi) );
				RC_TRY( new_matrix_db_end(new) );
				RC_TRY( write_matrix_db_end(mdb) );
				return(SPECIAL_RC);
			}else{
				RC_TRY( write_matrix_db_end(mdb) );
				return(rc);
			}
		}else{
			RC_TRY( write_matrix_db_end(mdb) );
			return(rc);
		}
	}

	return(NORMAL_RC);
}


static RC
add_row_child (MATRIX_DB_INFO *mdbi, MATRIX_DB *mdb, int index,
               MATRIX_DB_ARY *data, MATRIX_DB_ARY *buf, int row, int col,
               int which, int next_row, int next_col)
{
	RC rc;
	MATRIX_DB *child;

	RC_NULL_CHK( mdbi );
	RC_NULL_CHK( mdb );
	RC_NULL_CHK( data );
	RC_NULL_CHK( buf );
	if( (index < 0) || (index > mdb->size) ) return(ARG_ERROR_RC);
	if( (row < 0) || (col < 0) ) return(ARG_ERROR_RC);
	if( (next_row < 0) || (next_col < 0) ) return(ARG_ERROR_RC);
	if( (which < -1) || (which > 1) ) return(ARG_ERROR_RC);

	RC_TRY( read_matrix_db_start(mdb->p_file[index], mdb->p_offset[index],
	                             &child, mdbi) );
	rc = add_row_sub(mdbi, child, data, buf, row, col, which, next_row,
	                 next_col);
	RC_TRY( read_matrix_db_end(child) );

	return(rc);
}


static RC
add_row_leaf (MATRIX_DB_INFO *mdbi, MATRIX_DB *mdb, int index,
              MATRIX_DB_ARY *data, MATRIX_DB_ARY *buf, int row, int col,
              int which, int next_row, int next_col)
{
	int comp;
	double v_tmp[3][3];

	RC_NULL_CHK(mdbi);
	RC_NULL_CHK(mdb);
	RC_NULL_CHK(data);
	RC_NULL_CHK(buf);
	if( (index < 0) || (index > mdb->size) )return(ARG_ERROR_RC);
	if( (row < 0) || (col < 0) ) return(ARG_ERROR_RC);
	if( (next_row < 0) || (next_col < 0) ) return(ARG_ERROR_RC);
	if( (which < -1) || (which > 1) ) return(ARG_ERROR_RC);

	if(mdb->size == MATRIX_DB_SIZE) return(SPECIAL_RC);

	while(1){
		/* 親ノードのキーより大きい値になったらrootからやり直し */
		comp = compar_row_col(row, col, next_row, next_col);
		if(comp >= 0){ /* key1 >= key2 */
			return(SPECIAL_RC);
		}

		if(index == mdb->size){
			/* mdbの後半に入力 */
			comp = -1;
		}else{
			comp = compar_row_col(row, col, mdb->row[index], mdb->col[index]);
		}
#if 0
		int ii1;
		RC_TRY( print_matrix_db_debug(stderr, mdb) );
		log_printf(5, "[%d] comp = %d\n", __LINE__, comp);
		log_printf(5, "[%d] row = %d\n", __LINE__, row);
		log_printf(5, "[%d] col = %d\n", __LINE__, col);
		log_printf(5, "[%d] which = %d\n", __LINE__, which);
		log_printf(5, "[%d] next_row = %d\n", __LINE__, next_row);
		log_printf(5, "[%d] next_col = %d\n", __LINE__, next_col);
		log_printf(5, "[%d] index = %d\n", __LINE__, index);
		log_printf(5, "[%d] mdb->size = %d\n", __LINE__, mdb->size);
		log_printf(5, "[%d] buf->size = %d\n", __LINE__, buf->size);
		log_printf(5, "[%d] buf->index = %d\n", __LINE__, buf->index);
		log_printf(5, "[%d] data->size = %d\n", __LINE__, data->size);
		log_printf(5, "[%d] data->index = %d\n", __LINE__, data->index);
		for(ii1=buf->index; ii1<buf->size; ii1++){
			log_printf(5, "[%d] buf->row,col = %d %d %f\n", __LINE__,
			           buf->row[ii1], buf->col[ii1], buf->v[ii1][0][0]);
		}
#endif

		if(comp == -1){ /* key1 < key2 */
			if(index < mdb->size){
				/* mdb->*[index]をバッファへ退避 */
				buf->row[buf->size] = mdb->row[index];
				buf->col[buf->size] = mdb->col[index];
				RC_TRY( copy_matrix33(mdb->v[index], buf->v[buf->size]) );
				buf->size++;
			}
			if(which == -1){ /* data < buf */
				RC_TRY( set_v_matrix_db(row, col, data->v[data->index],
				                        index, mdb) );
				data->index++;
			}else if(which == 0){ /* data = buf */
				RC_TRY( add_matrix33a(data->v[data->index],
				                      buf->v[buf->index], v_tmp) );
				RC_TRY( set_v_matrix_db(row, col, v_tmp, index, mdb) );
				data->index++;
				buf->index++;
			}else{ /* (which == 1), data > buf */
				RC_TRY( set_v_matrix_db(row, col, buf->v[buf->index], index,
				                        mdb) );
				buf->index++;
			}
			if(index == mdb->size){
				mdb->size++;
				if(mdb->size == MATRIX_DB_SIZE) return(SPECIAL_RC);
			}
		}else if(comp == 0){ /* key1 = key2 */
			if(which == -1){ /* data < buf */
				add_matrix33(mdb->v[index], data->v[data->index]);
				data->index++;
			}else if(which == 0){ /* data = buf */
				RC_TRY( add_matrix33a(data->v[data->index],
				                      buf->v[buf->index], v_tmp) );
				add_matrix33(mdb->v[index], v_tmp);
				data->index++;
				buf->index++;
			}else{ /* (which == 1), data > buf */
				add_matrix33(mdb->v[index], buf->v[buf->index]);
				buf->index++;
			}
		}else{ /* (comp == 1) key1 > key2 */
			/* do nothing */
		}
		index++;
		if(buf->size == buf->alloc_size){
			RC_TRY( add_row_clean_buf(buf) );
			if(buf->size >= buf->alloc_size) return(OVERFLOW_ERROR_RC);
		}
		/* 次に入力する row colを求める */
		RC_TRY( add_row_min_row_col(data, buf, &row, &col, &which) );
		if(row < 0){
			/* 入力するデータがないので終了 */
			return(NORMAL_RC);
		}
	}

	return(NORMAL_RC);
}


/* data, buf の各index番目を比較する                          */
/* data のrow,colが最小なら which = -1, bufなら 1, 共通なら 0 */
static RC
add_row_min_row_col (MATRIX_DB_ARY *data, MATRIX_DB_ARY *buf,
                     int *row, int *col, int *which)
{

	RC_NULL_CHK(data);
	RC_NULL_CHK(buf);
	RC_NULL_CHK(row);
	RC_NULL_CHK(col);
	RC_NULL_CHK(which);

	*row = -1;
	*col = -1;
	*which = 0;

	if(data->index == data->size){ /* data には被入力データが無い */
		if(buf->index != buf->size){
			if(buf->index >= buf->size) return(UNKNOWN_ERROR_RC);
			/* bufにはある */
			*row = buf->row[buf->index];
			*col = buf->col[buf->index];
			*which = 1;
		}
	}else if(buf->index == buf->size){ /* buf には被入力データが無い */
		if(data->index != data->size){
			if(data->index >= data->size) return(UNKNOWN_ERROR_RC);
			/* dataにはある */
			*row = data->row[0];
			*col = data->col[data->index];
			*which = -1;
		}
	}else{ /* dataにもbufにも被入力データがある */
		if(buf->index  >= buf->size) return(UNKNOWN_ERROR_RC);
		if(data->index >= data->size) return(UNKNOWN_ERROR_RC);
		*which = compar_row_col(data->row[0], data->col[data->index],
		                        buf->row[buf->index], buf->col[buf->index]);
		if(*which == -1){ /* bufの方が大きい */
			*row = data->row[0];
			*col = data->col[data->index];
		}else if(*which == 0){ /* 同じ */
			*row = data->row[0];
			*col = data->col[data->index];
		}else{ /* (*which == 1) dataの方が大きい */
			*row = buf->row[buf->index];
			*col = buf->col[buf->index];
		}
	}

	return(NORMAL_RC);
}


/*
 * row 行目の col列から後半の登録されたデータ全部を array にコピー，
 * 列は col_arrayにコピー，コピーした数を size に代入 array, col_array
 * は呼び出し時に十分な大きさが確保されているものとする
 */
RC
get_row_matrix_db (MATRIX_DB_INFO *mdbi, int row, int *size,
                   double array[][3][3], int col_array[])
{
	RC_TRY( get_row_range_matrix_db(mdbi, row, 0, INT_MAX, size, array,
	                                col_array) );
	return(NORMAL_RC);
}


/* row 行目の col列から後半の登録されたデータ全部を array にコピー */
RC
get_row_col_matrix_db (MATRIX_DB_INFO *mdbi, int row, int col, int *size,
                       double array[][3][3], int col_array[])
{
	RC_TRY( get_row_range_matrix_db(mdbi, row, col, INT_MAX, size, array,
	                                col_array) );
	return(NORMAL_RC);
}


/* row 行目の col_min列(含む)からcol_max列(含まない)までの登録された */
/* データを array にコピー                                           */
RC
get_row_range_matrix_db (MATRIX_DB_INFO *mdbi, int row, int col_min,
                         int col_max, int *size, double array[][3][3],
                         int col_array[])
{
	*size = 0;
	RC_TRY( get_row_range_sub(mdbi, mdbi->root, row, col_min, col_max, size,
	                          array, col_array) );
	return(NORMAL_RC);
}


static RC
get_row_range_sub (MATRIX_DB_INFO *mdbi, MATRIX_DB *mdb, int row, int col_min,
                   int col_max, int *size, double array[][3][3],
                   int col_array[])
{
	int ii1;
	int tmp_file;
	WRAP64_INT tmp_offset;
	int found_flag; 
	int index;

	RC_NULL_CHK(mdbi);
	RC_NULL_CHK(mdb);
	RC_NULL_CHK(size);
	RC_NULL_CHK(array);
	RC_NULL_CHK(col_array);
	if( (row < 0)||(col_min < 0) ) return(ARG_ERROR_RC);
	if(col_min >= col_max) return(ARG_ERROR_RC);

	RC_TRY( search_matrix_db(row, col_min, mdb, &found_flag, &index, 
	                         &tmp_file, &tmp_offset) );
	if( (index < 0)||(index > mdb->size) ) return(UNKNOWN_ERROR_RC);

	if( is_node_matrix_db(mdb) ){
		if(index == mdb->size){
			/* mdbにrowが無い場合 */
			RC_TRY( get_row_range_child(mdbi, mdb, index, row, col_min,
			                            col_max, size, array, col_array) );
		}else if(mdb->row[index] > row){
			/* mdbにrowが無い場合 */
			RC_TRY( get_row_range_child(mdbi, mdb, index, row, col_min,
			                            col_max, size, array, col_array) );
		}else if(mdb->row[index] < row){
			/* mdbにrowが無い場合 */
			RC_TRY( get_row_range_child(mdbi, mdb, index+1, row, col_min,
			                            col_max, size, array, col_array) );
		}else{ /* mdb->row[index] == row */
			if(mdb->col[index] > col_min){
				/* mdb->row[index]の左下を全部取ってくる */
				RC_TRY( get_row_range_child(mdbi, mdb, index, row, col_min,
				                            col_max, size, array, col_array) );
			}
			for(ii1=index; ii1<mdb->size; ii1++){
				if( is_range_matrix_db(mdb, ii1, row, col_min, col_max) ){
					RC_TRY( copy_matrix33(mdb->v[ii1], array[*size]) );
					col_array[*size] = mdb->col[ii1];
					(*size)++;
					/* mdb->row[ii1]の右下を全部取ってくる */
					RC_TRY( get_row_range_child(mdbi, mdb, ii1+1, row,
					                            col_min, col_max, size, array,
					                            col_array) );
				}else{
					if(ii1==index){
						if(mdb->col[ii1] < col_min){
							/* mdb->row[ii1]の右下を全部取ってくる */
							RC_TRY( get_row_range_child(mdbi, mdb, ii1+1, row,
							                            col_min, col_max, size,
							                            array, col_array) );
						}
					}
					break;
				}
			}
		}
	}else{
		for(ii1=index; ii1<mdb->size; ii1++){
			if( is_range_matrix_db(mdb, ii1, row, col_min, col_max) ){
				RC_TRY( copy_matrix33(mdb->v[ii1], array[*size]) );
				col_array[*size] = mdb->col[ii1];
				(*size)++;
			}else{
				break;
			}
		}
	}

	return(NORMAL_RC);
}


static RC
get_row_range_child (MATRIX_DB_INFO *mdbi, MATRIX_DB *mdb, int index, int row,
                     int col_min, int col_max, int *size, double array[][3][3],
                     int col_array[])
{
	MATRIX_DB *child;

	RC_NULL_CHK(mdbi);
	RC_NULL_CHK(mdb);
	RC_NULL_CHK(size);
	RC_NULL_CHK(array);
	RC_NULL_CHK(col_array);
	if( (index < 0) || (index > mdb->size) )return(ARG_ERROR_RC);
	if( (row < 0)||(col_min < 0) ) return(ARG_ERROR_RC);
	if(col_min >= col_max) return(ARG_ERROR_RC);

	RC_TRY( read_matrix_db_start(mdb->p_file[index], mdb->p_offset[index],
	                             &child, mdbi) );
	RC_TRY( get_row_range_sub(mdbi, child, row, col_min, col_max, size, array,
	                          col_array) );
	RC_TRY( read_matrix_db_end(child) );

	return(NORMAL_RC);
}


/* row 行目を消去 */
RC
delete_row_matrix_db (MATRIX_DB_INFO *mdbi, int row)
{
	RC_TRY( delete_row_range_matrix_db(mdbi, row, 0, INT_MAX) )
	return(NORMAL_RC);
}


/* row 行目の col_min列(含む)から後半を削除 */
RC
delete_row_col_matrix_db (MATRIX_DB_INFO *mdbi, int row, int col_min)
{
	RC_TRY( delete_row_range_matrix_db(mdbi, row, col_min, INT_MAX) )
	return(NORMAL_RC);
}


/* row 行目の col_min列(含む)からcol_max列(含まない)までを消去 */
RC
delete_row_range_matrix_db (MATRIX_DB_INFO *mdbi, int row, int col_min,
                            int col_max)
{
	RC rc;

	RC_NULL_CHK(mdbi);
	if( (row < 0) || (col_min < 0) ) return(ARG_ERROR_RC);
	if(col_min >= col_max) return(ARG_ERROR_RC);

	while(1){ 
		if(mdbi->root->size == 0) return(NORMAL_RC);
		rc = delete_row_range_sub(mdbi, mdbi->root, row, col_min, col_max);
		if(rc == SPECIAL_RC){
			continue;
		}else{
			return(rc);
		}
	}

	return(NORMAL_RC);
}


static RC
delete_row_range_sub (MATRIX_DB_INFO *mdbi, MATRIX_DB *mdb, int row,
                      int col_min, int col_max)
{
	int ii1, ii2;
	int p_file;
	WRAP64_INT p_offset;
	int found_flag; 
	int index;
	RC rc;

	RC_NULL_CHK( mdbi );
	RC_NULL_CHK( mdb );
	if(mdb->size <= 0) return(ARG_ERROR_RC);
	if( (row < 0) || (col_min < 0) ) return(ARG_ERROR_RC);
	if(col_min >= col_max) return(ARG_ERROR_RC);

	RC_TRY( search_matrix_db(row, col_min, mdb, &found_flag, &index, 
	                         &p_file, &p_offset) );
	if( (index < 0)||(index > mdb->size) ) return(UNKNOWN_ERROR_RC);

	if( is_node_matrix_db(mdb) ){
		if(index == mdb->size){
			/* mdbにrowが無い場合 */
			rc = delete_row_range_child(mdbi, mdb, index, row, col_min,
			                            col_max);
			if(rc != NORMAL_RC) return(rc);
		}else if(mdb->row[index] > row){
			/* mdbにrowが無い場合 */
			rc = delete_row_range_child(mdbi, mdb, index, row, col_min,
			                            col_max);
			if(rc != NORMAL_RC) return(rc);
		}else if(mdb->row[index] < row){
			/* mdbにrowが無い場合 */
			rc = delete_row_range_child(mdbi, mdb, index+1, row, col_min,
			                            col_max);
			if(rc != NORMAL_RC) return(rc);
		}else{ /* mdb->row[index] == row */
			int index2;
			RC_TRY( write_matrix_db_start(mdb->p_file_self,
			                              mdb->p_offset_self, &mdb, mdbi) );
			/* 一斉削除できるTreeがあれば先に消す */
			for(ii1=index+1; ii1<mdb->size; ii1++){
				if( is_range_matrix_db(mdb, ii1, row, col_min, col_max) ){
					/* Tree全消去 */
					RC_TRY( delete_row_range_tree(mdbi, mdb->p_file[ii1],
					                                    mdb->p_offset[ii1]) );
				}else{
					break;
				}
			}
			if(ii1 == mdb->size){
				mdb->p_file[index+1]   = mdb->p_file[mdb->size];
				mdb->p_offset[index+1] = mdb->p_offset[mdb->size];
				mdb->size = index+1;
				RC_TRY( clean_matrix_db(mdb) );
			}else if(ii1 > index+1){
				for(ii2=ii1, index2=index+1; ii2<mdb->size; ii2++){
					RC_TRY( set_v_matrix_db(mdb->row[ii2], mdb->col[ii2],
					                        mdb->v[ii2], index2, mdb) );
					mdb->p_file[index2]   = mdb->p_file[ii2];
					mdb->p_offset[index2] = mdb->p_offset[ii2];
					index2++;
				}
				mdb->p_file[index2]   = mdb->p_file[ii2];
				mdb->p_offset[index2] = mdb->p_offset[ii2];
				mdb->size = index2;
				RC_TRY( clean_matrix_db(mdb) );
			}else{
				/* do nothing */
			}
			/* nodeが少なくなっていたらバランスを取ってからやり直し */
			if( !is_root_matrix_db(mdb) && (mdb->size < MATRIX_DB_SIZE/2) ){
				while( !is_root_matrix_db(mdb)
				    && (mdb->size < MATRIX_DB_SIZE/2) ){
					RC_TRY( balance_node_matrix_db(mdb, mdbi) );
					if( (mdb->size  == 0) || is_root_matrix_db(mdb) ){
						RC_TRY( write_matrix_db_end(mdb) );
						return(SPECIAL_RC);
					}
				}
				RC_TRY( write_matrix_db_end(mdb) );
				return(SPECIAL_RC);
			}
			/* 一斉消去終了 */
			/* index左のTree消去 */
			if(mdb->col[index] > col_min){
				rc = delete_row_range_child(mdbi, mdb, index, row,
				                            col_min, col_max);
				if(rc != NORMAL_RC){
					RC_TRY( write_matrix_db_end(mdb) );
					return(rc);
				}
			}
			/* index右のTree消去 */
			if(mdb->col[index] < col_max-1){
				rc = delete_row_range_child(mdbi, mdb, index+1, row,
				                            col_min, col_max);
				if(rc != NORMAL_RC){
					RC_TRY( write_matrix_db_end(mdb) );
					return(rc);
				}
			}
			/* node の index番目を消去 */
			if( is_range_matrix_db(mdb, index, row, col_min, col_max) ){
				RC_TRY( delete_v_node_matrix_db(index, mdb, mdbi) );
				RC_TRY( write_matrix_db_end(mdb) );
				return(SPECIAL_RC);
			}
			RC_TRY( write_matrix_db_end(mdb) );
		}
	}else{
		if(index == mdb->size) return(NORMAL_RC);
		if( is_range_matrix_db(mdb, index, row, col_min, col_max) ){
			RC_TRY( write_matrix_db_start(mdb->p_file_self,
			                              mdb->p_offset_self, &mdb, mdbi) );
			RC_TRY( delete_row_range_leaf(mdbi, mdb, index, row, col_min,
			                              col_max) );
			if( !is_root_matrix_db(mdb) && (mdb->size < MATRIX_DB_SIZE/2) ){
				while( !is_root_matrix_db(mdb)
				    && (mdb->size < MATRIX_DB_SIZE/2) ){
					RC_TRY( balance_leaf_matrix_db(mdb, mdbi) );
					if( (mdb->size  == 0) || is_root_matrix_db(mdb) ){
						RC_TRY( write_matrix_db_end(mdb) );
						return(SPECIAL_RC);
					}
				}
				RC_TRY( write_matrix_db_end(mdb) );
				return(SPECIAL_RC);
			}
			RC_TRY( write_matrix_db_end(mdb) );
		}
	}
	
	return(NORMAL_RC);
}


static RC
delete_row_range_child (MATRIX_DB_INFO *mdbi, MATRIX_DB *mdb, int index,
                        int row, int col_min, int col_max)
{
	int p_file;
	WRAP64_INT p_offset;
	MATRIX_DB *child;
	RC rc;

	RC_NULL_CHK( mdbi );
	RC_NULL_CHK( mdb );
	if(mdb->size <= 0) return(ARG_ERROR_RC);
	if( (index < 0) || (index > mdb->size) )return(ARG_ERROR_RC);
	if( (row < 0) || (col_min < 0) ) return(ARG_ERROR_RC);
	if(col_min >= col_max) return(ARG_ERROR_RC);

	RC_TRY( read_matrix_db_start(mdb->p_file[index],
	                             mdb->p_offset[index], &child, mdbi) );
	rc = delete_row_range_sub(mdbi, child, row, col_min, col_max);
	if(rc == NORMAL_RC){
		RC_TRY( read_matrix_db_end(child) );
		if(mdb->size == 0){
			return(SPECIAL_RC);
		}else{
			return(NORMAL_RC);
		}
	}
	if(rc != SPECIAL_RC) return(rc);

	if(child->size == 0){
		p_file   = child->p_file_self;
		p_offset = child->p_offset_self;
		if(child->use_count != 1){
			return(UNKNOWN_ERROR_RC);
		}
		RC_TRY( read_matrix_db_end(child) );
		RC_TRY( delete_matrix_db(mdbi, p_file, p_offset) );
	}else{
		RC_TRY( read_matrix_db_end(child) );
	}
	mdbi->route_size = 0; /* 必須 */

	return(SPECIAL_RC);
}


/* p_file, p_offset 以下のTreeを全消去する */
static RC
delete_row_range_tree (MATRIX_DB_INFO *mdbi, int p_file, WRAP64_INT p_offset)
{
	int  ii1;
	MATRIX_DB *mdb;
	MATRIX_DB *child;

	RC_NULL_CHK( mdbi );
	if(p_file <= 0) return(ARG_ERROR_RC);
	if(p_offset < 0) return(ARG_ERROR_RC);

	RC_TRY( read_matrix_db_start(p_file, p_offset, &mdb, mdbi) );
	if( is_node_matrix_db(mdb) ){
		RC_TRY( read_matrix_db_start(mdb->p_file[0], mdb->p_offset[0],
		                             &child, mdbi) );
		if( is_leaf_matrix_db(child) ){
			RC_TRY( read_matrix_db_end(child) );
			for(ii1=0; ii1<=mdb->size; ii1++){
				RC_TRY( delete_matrix_db(mdbi, mdb->p_file[ii1],
				                               mdb->p_offset[ii1]) );
			}
		}else{
			RC_TRY( read_matrix_db_end(child) );
			for(ii1=0; ii1<=mdb->size; ii1++){
				RC_TRY( delete_row_range_tree(mdbi, mdb->p_file[ii1],
				                                    mdb->p_offset[ii1]) );
			}
		}
	}
	RC_TRY( read_matrix_db_end(mdb) );
	RC_TRY( delete_matrix_db(mdbi, p_file, p_offset) );

	return(NORMAL_RC);
}


static RC
delete_row_range_leaf (MATRIX_DB_INFO *mdbi, MATRIX_DB *mdb, int index,
                       int row, int col_min, int col_max)
{
	int ii1, ii2;
	int row0, col0;

	RC_NULL_CHK(mdb);
	RC_NULL_CHK(mdbi);
	if(mdb->size <= 0) return(UNDERFLOW_ERROR_RC);
	if( (index < 0) || (index > mdb->size) ) return(ARG_ERROR_RC);
	if(index >= mdb->size) return(ARG_ERROR_RC);
	if(row < 0) return(ARG_ERROR_RC);
	if(col_min < 0) return(ARG_ERROR_RC);
	if(col_min >= col_max) return(ARG_ERROR_RC);

	row0 = mdb->row[0];
	col0 = mdb->col[0];
	for(ii1=index+1; ii1<mdb->size; ii1++){
		if( !is_range_matrix_db(mdb, ii1, row, col_min, col_max) ){
			for(ii2=ii1; ii2<mdb->size; ii2++){
				RC_TRY( set_v_matrix_db(mdb->row[ii2], mdb->col[ii2],
				                        mdb->v[ii2], index, mdb) );
				index++;
			}
			break;
		}
	}
	mdb->size = index;
	RC_TRY( clean_matrix_db(mdb) );

	if(mdb->size == 0){
		mdb->key_row = row0;
		mdb->key_col = col0;
	}

	return(NORMAL_RC);
}


/* B-Treeにvを入力 */
static RC
input_v_btree_matrix_db (int row, int col, double v[][3], int index,
                         MATRIX_DB *mdb, MATRIX_DB_INFO *mdbi)
{
	int found_flag;
	int p_file;
	WRAP64_INT p_offset;

	RC_NULL_CHK( mdb );
	RC_NULL_CHK( mdbi );
	RC_NULL_CHK( v );
	if( (row < 0) || (col < 0) ) return(ARG_ERROR_RC);
	if( (index < 0) || (index >= MATRIX_DB_SIZE) ) return(ARG_ERROR_RC);
	if( !is_leaf_matrix_db(mdb) ) return(ARG_ERROR_RC); /* 葉に到達してない */

	/* 葉が一杯なら分割する */
	if(mdb->size >= MATRIX_DB_SIZE){
		MATRIX_DB *new_mdb, *new_child;

		RC_TRY( new_matrix_db_start(&new_child, mdbi) );
		RC_TRY( split_btree_matrix_db(mdb, new_child, mdbi) );
		RC_TRY( new_matrix_db_end(new_child) );

		/* 分割後の被入力葉を再検索 (注: 親から検索してはいけない) */
		mdbi->route_size = 0;
		RC_TRY( search_v_btree_matrix_db(mdbi, mdbi->root, row, col,
		                                 &found_flag, &index, &p_file,
		                                 &p_offset) );
		if(p_file <= 0) return(UNKNOWN_ERROR_RC);
		if(p_file > MATRIX_DB_SF_SIZE) return(UNKNOWN_ERROR_RC);
		if(p_offset < 0)  return(UNKNOWN_ERROR_RC);

		RC_TRY( write_matrix_db_start(p_file, p_offset, &new_mdb, mdbi) );
		if( (index < 0)||(index > new_mdb->size) ) return(UNKNOWN_ERROR_RC);
		RC_TRY( insert_v_leaf_matrix_db(row, col, v, index, new_mdb) );
		RC_TRY( write_matrix_db_end(new_mdb) );
	}else{
		RC_TRY( insert_v_leaf_matrix_db(row, col, v, index, mdb) );
	}

	return(NORMAL_RC);
}


/*
 * mdb以下の Tree内から row,colを含む MATRIX_DBを探索する．
 *
 * 通常，mdb = mdbi->rootとして RootNodeから検索する．row,colが一致する
 * 場所に辿り付いた場合，found_flag = 1となり，該当する MATRIX_DBの p_file,
 * p_offsetがセットされ，該当個所の配列 indexがセットされる．
 * found_flag = 0 の場合，row,colを insertする時の MATRIX_DBの p_file,
 * p_offsetがセットされ，insertする時の配列 indexがセットされる．
 */
RC
search_v_btree_matrix_db (MATRIX_DB_INFO *mdbi, MATRIX_DB *mdb,
                          int row, int col, int *found_flag, int *index,
                          int *p_file, WRAP64_INT *p_offset)
{
	int next_p_file;
	WRAP64_INT next_p_offset;
	MATRIX_DB *child;

	RC_NULL_CHK( mdb );
	RC_NULL_CHK( mdbi );
	RC_NULL_CHK( index );
	RC_NULL_CHK( found_flag );
	RC_NULL_CHK( p_file );
	RC_NULL_CHK( p_offset );
	if( (row < 0) || (col < 0) ) return(ARG_ERROR_RC);

	/* 経路情報を記録 */
	mdbi->route_file[mdbi->route_size]   = mdb->p_file_self;
	mdbi->route_offset[mdbi->route_size] = mdb->p_offset_self;
	mdbi->route_size++;
	if(mdbi->route_size >= MATRIX_DB_ROUTE_SIZE) return(OVERFLOW_ERROR_RC);

	RC_TRY( search_matrix_db(row, col, mdb, found_flag, index, 
	                         &next_p_file, &next_p_offset) );
	if(*found_flag){
		/* 一致するMATRIX_DBが見つかった */
		*p_file    = mdb->p_file_self;
		*p_offset  = mdb->p_offset_self;
	}else{
		if( is_leaf_matrix_db(mdb) ){
			/* 葉に到達した．次に探索するものがない */
			*p_file   = mdb->p_file_self;
			*p_offset = mdb->p_offset_self;
		}else{
			/* mdbは内部節点．葉になるまで次を探索 */
			RC_TRY( read_matrix_db_start(next_p_file, next_p_offset,
			                             &child, mdbi) );
			RC_TRY( search_v_btree_matrix_db(mdbi, child, row, col, found_flag,
			                                 index, p_file, p_offset) );
			RC_TRY( read_matrix_db_end(child) );
		}
	}

	return(NORMAL_RC);
}


static RC
search_parent_matrix_db (MATRIX_DB_INFO *mdbi, MATRIX_DB *mdb, int *index)
{
	int ii1;
	int found_flag;
	int idx;
	MATRIX_DB *tmp;
	int p_file;
	WRAP64_INT p_offset;

	RC_NULL_CHK( mdbi );
	RC_NULL_CHK( mdb );
	RC_NULL_CHK( index );
	if(mdbi->route_size < 0) return(ARG_ERROR_RC);
	if(mdb->size < 0) return(ARG_ERROR_RC);

	if( is_root_matrix_db(mdb) ){
		*index = -1;
		return(NORMAL_RC);
	}

	if(mdbi->route_count > MATRIX_DB_ROUTE_SIZE) return(SEEK_ERROR_RC);

	/* 経路情報が無いのでmdbまでの経路を空サーチ */
	if(mdbi->route_size == 0){
		if(mdb->size == 0){
			if((mdb->key_row < 0) || (mdb->key_col < 0)) return(ARG_ERROR_RC);
			RC_TRY( search_v_btree_matrix_db(mdbi, mdbi->root, mdb->key_row,
			                                 mdb->key_col, &found_flag, &idx,
			                                 &p_file, &p_offset) );
		}else{
			RC_TRY( search_v_btree_matrix_db(mdbi, mdbi->root, mdb->row[0],
			                                 mdb->col[0], &found_flag, &idx,
			                                 &p_file, &p_offset) );
		}
		if(idx < 0) return(UNKNOWN_ERROR_RC);
		if(mdbi->route_size < 0) return(UNKNOWN_ERROR_RC);
	}

	*index = -1;
	for(ii1=mdbi->route_size-1; ii1>0; ii1--){
		if( (mdbi->route_file[ii1]   == mdb->p_file_self)
		 && (mdbi->route_offset[ii1] == mdb->p_offset_self) ){
			*index = ii1-1;
			mdbi->route_count = 0;
			break;
		}
	}

#if 1
	/* チェック */
	if(*index >= 0){
		RC_TRY( read_matrix_db_start(mdbi->route_file[*index],
		                             mdbi->route_offset[*index], &tmp, mdbi) );
		if( !is_parent_matrix_db(mdb, tmp) ){
			/*
			RC_TRY( print_matrix_db_debug(stderr, mdb) );
			RC_TRY( print_matrix_db_debug(stderr, tmp) );
			RC_TRY( print_btree_matrix_db(stderr, mdbi) );
			for(ii1=0; ii1<mdbi->route_size; ii1++){
				log_printf(5, "route_offset[%d] = %lld\n", ii1,
				           mdbi->route_offset[ii1]);
			}
			*/
			return(UNKNOWN_ERROR_RC);
		}
		RC_TRY( read_matrix_db_end(tmp) );
	}
#endif

	if(*index < 0){
		/* 経路情報を再取得するために空サーチをする */
		mdbi->route_size = 0;
		mdbi->route_count++;
		if(mdb->size == 0){
			if((mdb->key_row < 0) || (mdb->key_col < 0)) return(ARG_ERROR_RC);
			RC_TRY( search_v_btree_matrix_db(mdbi, mdbi->root, mdb->key_row,
			                                 mdb->key_col, &found_flag, &idx,
			                                 &p_file, &p_offset) );
		}else{
			RC_TRY( search_v_btree_matrix_db(mdbi, mdbi->root, mdb->row[0],
			                                  mdb->col[0], &found_flag, &idx,
			                                  &p_file, &p_offset) );
		}
		if(idx < 0) return(UNKNOWN_ERROR_RC);
		RC_TRY( search_parent_matrix_db(mdbi, mdb, index) );
		if(*index < 0) return(UNKNOWN_ERROR_RC);
	}

	return(NORMAL_RC);
}


/* childをを分割する */
RC
split_btree_matrix_db (MATRIX_DB *child, MATRIX_DB *new_child,
                       MATRIX_DB_INFO *mdbi)
{
	int parent_new_flag = 0;
	int parent_change_flag = 0;
	int parent_file;
	WRAP64_INT parent_offset;
	MATRIX_DB *parent;

	RC_NULL_CHK( child );
	RC_NULL_CHK( new_child );
	RC_NULL_CHK( mdbi );

	/* 親節点の取得 */
	if(is_root_matrix_db(child)){
		/* RootNodeだったchildを分離する時に作られる新 RootNode */
		RC_TRY( new_matrix_db_start(&parent, mdbi) );
		parent->p_file[0]   = child->p_file_self;
		parent->p_offset[0] = child->p_offset_self;
		RC_TRY( set_root_matrix_db(parent, mdbi) );
		child->p_file_parent = -1;
		child->p_offset_parent = -1;
		parent_new_flag = 1;
		parent_file   = parent->p_file_self;
		parent_offset = parent->p_offset_self;
	}else{
		int pindex;
		RC_TRY( search_parent_matrix_db(mdbi, child, &pindex) );
		if(pindex < 0) return(UNKNOWN_ERROR_RC);
		parent_file   = mdbi->route_file[pindex];
		parent_offset = mdbi->route_offset[pindex];

		RC_TRY( write_matrix_db_start(parent_file, parent_offset,
		                              &parent, mdbi) );
		/* 親が一杯なら親を分割 */
		if(parent->size >= MATRIX_DB_SIZE){
			MATRIX_DB *new_parent;

			RC_TRY( new_matrix_db_start(&new_parent, mdbi) );
			RC_TRY( split_btree_matrix_db(parent, new_parent, mdbi) );
			if( is_parent_matrix_db(child, new_parent) ){
				int index;
				mdbi->route_size = 0;
				RC_TRY( search_parent_matrix_db(mdbi, child, &index) );
				if(index < 0) return(UNKNOWN_ERROR_RC);
				mdbi->route_file[index]   = new_parent->p_file_self;
				mdbi->route_offset[index] = new_parent->p_offset_self;
				parent_file   = new_parent->p_file_self;
				parent_offset = new_parent->p_offset_self;
			}
			RC_TRY( new_matrix_db_end(new_parent) );
			parent_change_flag = 1;
		}
	}

	if(parent_new_flag == 1){
		RC_TRY( new_matrix_db_end(parent) );
		RC_TRY( write_matrix_db_start(parent_file, parent_offset,
		                              &parent, mdbi) );
	}else if(parent_change_flag == 1){
		RC_TRY( write_matrix_db_end(parent) );
		RC_TRY( write_matrix_db_start(parent_file, parent_offset,
		                              &parent, mdbi) );
	}

	/* 値をそれぞれに移動 */
	RC_TRY( split_btree_matrix_db_sub(child, parent, new_child) );

	RC_TRY( write_matrix_db_end(parent) );

	return(NORMAL_RC);
}


/*
 * FULLになった child を分割する時の移動作業をする．
 *
 * 例:
 * MATRIX_DB_SIZE = 89 => (44 + 1 + 44)  に分断
 * 0〜(MATRIX_DB_SIZE/2) = 0〜43 -> 左の子(こっちは数が減るだけ)
 * MATRIX_DB_SIZE/2 = 44 -> 親へUP
 * (MATRIX_DB_SIZE/2+1)〜(MATRIX_DB_SIZE-1) = 45〜88 -> 右の子へ(新規作成)
 */
static RC
split_btree_matrix_db_sub (MATRIX_DB *child, MATRIX_DB *parent,
                           MATRIX_DB *new_child)
{
	int ii1, ii2;
	int index, found_flag;
	int p_file;
	WRAP64_INT p_offset;

	RC_NULL_CHK( child );
	RC_NULL_CHK( parent );
	RC_NULL_CHK( new_child );

	/* childは最低でもMATRIX_DB_SIZE/2よりも多くなければいけない */
	if(child->size < MATRIX_DB_SIZE/2) return(ARG_ERROR_RC);

	/* 親が埋まっていてはいけない */
	if(parent->size >= MATRIX_DB_SIZE) return(ARG_ERROR_RC); /* SIZE FULL */

	/* childの後半を new_childへ移動(+1は親に移動する分) */
	for(ii1=(MATRIX_DB_SIZE/2+1), ii2=0; ii1<child->size; ii1++){
		RC_TRY( set_v_matrix_db(child->row[ii1], child->col[ii1],
		                        child->v[ii1], ii2, new_child) );
		new_child->p_file[ii2]   = child->p_file[ii1];
		new_child->p_offset[ii2] = child->p_offset[ii1];
		ii2++;
	}
	new_child->p_file[ii2]   = child->p_file[ii1];
	new_child->p_offset[ii2] = child->p_offset[ii1];
	new_child->size = ii2;

	/* 親に挿入 */
	RC_TRY( search_matrix_db(child->row[MATRIX_DB_SIZE/2],
	                         child->col[MATRIX_DB_SIZE/2], parent, &found_flag,
	                         &index, &p_file, &p_offset) );
	if(found_flag) return(UNKNOWN_ERROR_RC);
	if( (index < 0) || (index > parent->size) ) return(UNKNOWN_ERROR_RC);

	RC_TRY( insert_v_node_matrix_db(child->row[MATRIX_DB_SIZE/2],
	                                child->col[MATRIX_DB_SIZE/2],
	                                child->v[MATRIX_DB_SIZE/2],
	                                new_child->p_file_self,
	                                new_child->p_offset_self, index, parent));
	/* child の後半を初期化 */
	child->size = MATRIX_DB_SIZE/2;
	RC_TRY( clean_matrix_db(child) );

	return(NORMAL_RC);
}


RC
delete_v_btree_matrix_db (int row, int col, MATRIX_DB_INFO *mdbi)
{
	int index, found_flag;
	MATRIX_DB *mdb;
	int p_file;
	WRAP64_INT p_offset;

	RC_NULL_CHK( mdbi );
	if( (row < 0) || (col < 0) ) return(ARG_ERROR_RC);

	/* 挿入すべき Leaf Nodeを Rootから探索 */
	mdbi->route_size = 0;
	RC_TRY( search_v_btree_matrix_db(mdbi, mdbi->root, row, col, &found_flag,
	                                 &index, &p_file, &p_offset) );
	if(p_file <= 0) return(UNKNOWN_ERROR_RC);
	if(p_file > MATRIX_DB_SF_SIZE) return(UNKNOWN_ERROR_RC);
	if(p_offset < 0)  return(UNKNOWN_ERROR_RC);

	if(found_flag == 0){
		RC_TRY( log_printf(6, "row = %d col = %d is missing. \n", row, col) );
		return(NORMAL_RC);
	}

	RC_TRY( write_matrix_db_start(p_file, p_offset, &mdb, mdbi) );
	if( (index < 0) || (index > mdb->size) ) return(UNKNOWN_ERROR_RC);

	if( is_leaf_matrix_db(mdb) ){
		RC_TRY( delete_v_leaf_matrix_db(index, mdb, mdbi) )
	}else if(is_node_matrix_db(mdb) ){
		RC_TRY( delete_v_node_matrix_db(index, mdb, mdbi) )
	}else{
		return(UNKNOWN_ERROR_RC);
	}

	if(mdb->size == 0){
		if( is_root_matrix_db(mdb) ){
			/* root は消さない */
			RC_TRY( write_matrix_db_end(mdb) );
		}else{
			p_file   = mdb->p_file_self;
			p_offset = mdb->p_offset_self;
			RC_TRY( write_matrix_db_end(mdb) );
			RC_TRY( delete_matrix_db(mdbi, p_file, p_offset) );
		}
	}else{
		RC_TRY( write_matrix_db_end(mdb) );
	}

	return(NORMAL_RC);
}


/* この関数は主に split_btree_matrix_db_sub() で使う */
RC
insert_v_node_matrix_db (int row, int col, double v[][3], int p_file,
                         WRAP64_INT p_offset, int index, MATRIX_DB *mdb)
{
	int ii1;

	RC_NULL_CHK(mdb);
	RC_NULL_CHK(v);
	if( (index < 0) || (index > mdb->size) ) return(ARG_ERROR_RC);
	if( (row < 0) || (col < 0) ) return(ARG_ERROR_RC);
	if( (p_file <= 0) || (p_file > MATRIX_DB_SF_SIZE) ) return(ARG_ERROR_RC);
	if(p_offset < 0)  return(ARG_ERROR_RC);

	/* SIZE FULL */
	if(mdb->size >= MATRIX_DB_SIZE) return(OVERFLOW_ERROR_RC);

	/* row,col,vだけ先に挿入 */
	RC_TRY( insert_v_leaf_matrix_db(row, col, v, index, mdb) );

	/* (index+1)番目より右を一つずらす */
	for(ii1=mdb->size; ii1>index+1; ii1--){
		mdb->p_file[ii1] = mdb->p_file[ii1-1];
		mdb->p_offset[ii1] = mdb->p_offset[ii1-1];
	}

	mdb->p_file[index+1] = p_file;
	mdb->p_offset[index+1] = p_offset;

	return(NORMAL_RC);
}


/* mdbのindex番目に row,col,vを入力．row-colが同じものがあれば上書きする */
RC
insert_v_leaf_matrix_db (int row, int col, double v[][3],
                         int index, MATRIX_DB *mdb)
{
	int ii1;

	RC_NULL_CHK(mdb);
	RC_NULL_CHK(v);
	if( (index < 0) || (index > mdb->size) ) return(ARG_ERROR_RC);
	if( (row < 0) || (col < 0) ) return(ARG_ERROR_RC);

	/* SIZE FULL */
	if(mdb->size >= MATRIX_DB_SIZE) return(OVERFLOW_ERROR_RC);

	for(ii1=mdb->size; ii1>index; ii1--){
		RC_TRY( set_v_matrix_db(mdb->row[ii1-1], mdb->col[ii1-1],
		                        mdb->v[ii1-1], ii1, mdb) );
	}
	(mdb->size)++;

	RC_TRY( set_v_matrix_db(row, col, v, index, mdb) );

	return(NORMAL_RC);
}


RC
delete_v_node_matrix_db (int index, MATRIX_DB *mdb, MATRIX_DB_INFO *mdbi)
{
	int ii1;
	MATRIX_DB *child_prev, *child_next;
	int p_file;
	WRAP64_INT p_offset;

	RC_NULL_CHK(mdb);
	RC_NULL_CHK(mdbi);
	if( (index < 0) || (index >= mdb->size) ) return(ARG_ERROR_RC);
	if( !is_node_matrix_db(mdb) ) return(ARG_ERROR_RC);

	RC_TRY( read_matrix_db_start(mdb->p_file[index],
	                             mdb->p_offset[index], &child_prev, mdbi) );
	RC_TRY( read_matrix_db_start(mdb->p_file[index+1],
	                             mdb->p_offset[index+1], &child_next, mdbi) );

	if(child_next->size > MATRIX_DB_SIZE/2){
		/* child_nextを頂点とする木の中から mdb[index]の直後のvを持つ葉  */
		/* まで降りて行き，その値を消し，その値で mdb[index]を上書きする */
		MATRIX_DB *mdb_next = child_next;
		while( !is_leaf_matrix_db(mdb_next) ){
			MATRIX_DB *tmp;
			RC_TRY( read_matrix_db_start(mdb_next->p_file[0],
			                             mdb_next->p_offset[0], &tmp, mdbi) );
			RC_TRY( read_matrix_db_end(mdb_next) );
			mdb_next = tmp;
		}
		RC_TRY( write_matrix_db_start(mdb_next->p_file_self,
		                              mdb_next->p_offset_self,
		                              &mdb_next, mdbi) );
		RC_TRY( set_v_matrix_db(mdb_next->row[0], mdb_next->col[0],
		                        mdb_next->v[0], index, mdb) );
		RC_TRY( delete_v_leaf_matrix_db(0, mdb_next, mdbi) );
		RC_TRY( write_matrix_db_end(mdb_next) );

		if(mdb_next->size == 0){
			p_file   = mdb_next->p_file_self;
			p_offset = mdb_next->p_offset_self;
			RC_TRY( read_matrix_db_end(mdb_next) );
			RC_TRY( delete_matrix_db(mdbi, p_file, p_offset) );
		}else{
			RC_TRY( read_matrix_db_end(mdb_next) );
		}
		RC_TRY( read_matrix_db_end(child_prev) );
	}else if(child_prev->size > MATRIX_DB_SIZE/2){
		/* child_prevを頂点とする木の中から mdb[index]の直前のvを持つ葉  */
		/* まで降りて行き，その値を消し，その値で mdb[index]を上書きする */
		int max;
		MATRIX_DB *mdb_prev = child_prev;
		while( !is_leaf_matrix_db(mdb_prev) ){
			MATRIX_DB *tmp;
			RC_TRY( read_matrix_db_start(mdb_prev->p_file[mdb_prev->size],
			                             mdb_prev->p_offset[mdb_prev->size],
			                             &tmp, mdbi) );
			RC_TRY( read_matrix_db_end(mdb_prev) );
			mdb_prev = tmp;
		}
		RC_TRY( write_matrix_db_start(mdb_prev->p_file_self,
		                              mdb_prev->p_offset_self,
		                              &mdb_prev, mdbi) );
		max = mdb_prev->size - 1;
		RC_TRY( set_v_matrix_db(mdb_prev->row[max], mdb_prev->col[max],
		                        mdb_prev->v[max], index, mdb) );
		RC_TRY( delete_v_leaf_matrix_db(max, mdb_prev, mdbi) );
		RC_TRY( write_matrix_db_end(mdb_prev) );
		if(mdb_prev->size == 0){
			p_file   = mdb_prev->p_file_self;
			p_offset = mdb_prev->p_offset_self;
			RC_TRY( read_matrix_db_end(mdb_prev) );
			RC_TRY( delete_matrix_db(mdbi, p_file, p_offset) );
		}else{
			RC_TRY( read_matrix_db_end(mdb_prev) );
		}
		RC_TRY( read_matrix_db_end(child_next) );
	}else if( (child_prev->size + 1 + child_next->size) <= MATRIX_DB_SIZE){
		/* index の左右の子の数が少ないのでマージして，再度削除工程を行う */
		/* child_prev と child_nextの間にmdb[index]を落としてから削除する */
		int del_index = child_prev->size;
		RC_TRY( write_matrix_db_start(child_prev->p_file_self,
		                              child_prev->p_offset_self,
		                              &child_prev, mdbi) );
		RC_TRY( write_matrix_db_start(child_next->p_file_self,
		                              child_next->p_offset_self,
		                              &child_next, mdbi) );
		if( is_leaf_matrix_db(child_prev) ){
			if( (child_prev->size == 0) && (child_next->size == 0) ){
				for(ii1=index+1; ii1<mdb->size; ii1++){
					RC_TRY( set_v_matrix_db(mdb->row[ii1], mdb->col[ii1],
					                        mdb->v[ii1], index, mdb) );
					mdb->p_file[index+1]   = mdb->p_file[ii1+1];
					mdb->p_offset[index+1] = mdb->p_offset[ii1+1];
					index++;
				}
				mdb->size--;
				RC_TRY( clean_matrix_db(mdb) );
				RC_TRY( balance_leaf_matrix_db(child_prev, mdbi) );
			}else{
				RC_TRY( marge_leaf_matrix_db(child_prev, child_next,
				                             mdb, mdbi, index) );
				RC_TRY( delete_v_leaf_matrix_db(del_index, child_prev, mdbi) );
			}
		}else{
			RC_TRY( marge_node_matrix_db(child_prev, child_next,
			                             mdb, mdbi, index) );
			RC_TRY( delete_v_node_matrix_db(del_index, child_prev, mdbi) );
		}
		RC_TRY( write_matrix_db_end(child_prev) );
		RC_TRY( write_matrix_db_end(child_next) );
		if(child_prev->size == 0){
			return(UNKNOWN_ERROR_RC);
		}else{
			RC_TRY( read_matrix_db_end(child_prev) );
		}
		if(child_next->size == 0){
			p_file   = child_next->p_file_self;
			p_offset = child_next->p_offset_self;
			RC_TRY( read_matrix_db_end(child_next) );
			RC_TRY( delete_matrix_db(mdbi, p_file, p_offset) );
		}else{
			return(UNKNOWN_ERROR_RC);
		}
	}else{
		return(UNKNOWN_ERROR_RC);
	}

	/* (mdb->size > 0)必須 */
	if( !is_root_matrix_db(mdb) && (mdb->size > 0)
	                            && (mdb->size < MATRIX_DB_SIZE/2) ){
		RC_TRY( balance_node_matrix_db(mdb, mdbi) );
	}

	return(NORMAL_RC);
}


/* RootNodeではない内部節点のバランスをとる． */
RC
balance_node_matrix_db (MATRIX_DB *mdb, MATRIX_DB_INFO *mdbi)
{
	int ii1;
	int pindex; /* mdbを指す親の index */
	int p_file;
	WRAP64_INT p_offset;
	MATRIX_DB *parent, *right, *left;
	int p_index;

	RC_NULL_CHK(mdb);
	RC_NULL_CHK(mdbi);
	if( is_root_matrix_db(mdb) ) return(ARG_ERROR_RC);
	if( is_leaf_matrix_db(mdb) ) return(ARG_ERROR_RC);

	RC_TRY( search_parent_matrix_db(mdbi, mdb, &p_index) );
	if(p_index < 0) return(UNKNOWN_ERROR_RC);
	RC_TRY( write_matrix_db_start(mdbi->route_file[p_index],
	                              mdbi->route_offset[p_index],
	                              &parent, mdbi) );
	/* 念のためチェック: MATRIX_DB_SIZEが小さい時に稀に問題になる */
#if 0
	if( !is_parent_matrix_db(mdb, parent) ){
		RC_TRY( write_matrix_db_end(parent) );
		mdbi->route_size = 0;
		RC_TRY( search_parent_matrix_db(mdbi, mdb, &p_index) );
		RC_TRY( write_matrix_db_start(mdbi->route_file[p_index],
		                              mdbi->route_offset[p_index],
		                              &parent, mdbi) );
	}
#endif

	/* mdbを指す親の p_file[index] のindexを探す．総当り */
	for(ii1=0, pindex=-1; ii1<parent->size+1; ii1++){
		if( (mdb->p_file_self   == parent->p_file[ii1])
		 && (mdb->p_offset_self == parent->p_offset[ii1]) ){
			pindex = ii1;
			break;
		}
	}
	if(pindex < 0) return(UNKNOWN_ERROR_RC);

	if(parent->size > pindex){ /* 右側と調整 */
		RC_TRY( write_matrix_db_start(parent->p_file[pindex+1],
		                              parent->p_offset[pindex+1],
		                              &right, mdbi) );
		if((mdb->size + 1 + right->size) <= MATRIX_DB_SIZE){
			/* 右を mdbへマージする */
			RC_TRY( marge_node_matrix_db(mdb, right, parent, mdbi, pindex) );
			p_file   = right->p_file_self;
			p_offset = right->p_offset_self;
			RC_TRY( write_matrix_db_end(right) );
			RC_TRY( delete_matrix_db(mdbi, p_file, p_offset) );
		}else if(right->size > MATRIX_DB_SIZE/2){
			/* 右からmdbへ1個移動 */
			RC_TRY( move_from_right_node_matrix_db(mdb, right, parent,
			                                       mdbi, pindex) );
			RC_TRY( write_matrix_db_end(right) );
		}else{
			return(UNKNOWN_ERROR_RC);
		}
	}else if(pindex > 0){ /* 左側と調整 */
		RC_TRY( write_matrix_db_start(parent->p_file[pindex-1],
		                              parent->p_offset[pindex-1],
		                              &left, mdbi) );
		if((mdb->size + 1 + left->size) <= MATRIX_DB_SIZE){
			/* 左へ mdbをマージする */
			RC_TRY( marge_node_matrix_db(left, mdb, parent, mdbi, pindex-1) );
		}else if(left->size > MATRIX_DB_SIZE/2){
			/* mdbから左へ1個移動 */
			RC_TRY( move_from_left_node_matrix_db(left, mdb, parent,
			                                      mdbi, pindex-1) );
		}else{
			return(UNKNOWN_ERROR_RC);
		}
		RC_TRY( write_matrix_db_end(left) );
	}else{
		return(UNKNOWN_ERROR_RC);
	}

	if(parent->size == 0){
		if(parent->use_count > 1){
			/* 上位の関数で delete_matrix_db() する */
			RC_TRY( write_matrix_db_end(parent) );
		}else if(parent->use_count == 1){
			p_file   = parent->p_file_self;
			p_offset = parent->p_offset_self;
			RC_TRY( write_matrix_db_end(parent) );
			RC_TRY( delete_matrix_db(mdbi, p_file, p_offset) );
		}else{
			return(UNKNOWN_ERROR_RC);
		}
	}else{
		RC_TRY( write_matrix_db_end(parent) );
	}

	return(NORMAL_RC);
}


/* leftからrightへ1個移動する               */
/* indexは leftとrightを指している親のindex */
RC
move_from_left_node_matrix_db (MATRIX_DB *left, MATRIX_DB *right,
                               MATRIX_DB *parent, MATRIX_DB_INFO *mdbi,
                               int index)
{
	int ii1;

	RC_NULL_CHK(left);
	RC_NULL_CHK(right);
	RC_NULL_CHK(parent);
	if( (index < 0) || (index >= parent->size) ) return(ARG_ERROR_RC);

	/* rightを1個右にずらして，parentから1個移動 */
	for(ii1=right->size; ii1>0; ii1--){
		RC_TRY( set_v_matrix_db(right->row[ii1-1], right->col[ii1-1],
		                        right->v[ii1-1], ii1, right) );
		right->p_file[ii1+1]   = right->p_file[ii1];
		right->p_offset[ii1+1] = right->p_offset[ii1];
	}
	right->p_file[ii1+1]   = right->p_file[ii1];
	right->p_offset[ii1+1] = right->p_offset[ii1];
	right->size++;

	RC_TRY( set_v_matrix_db(parent->row[index], parent->col[index],
	                        parent->v[index], 0, right) );

	/* leftの右端を 1個 parentへ移動 */
	RC_TRY( set_v_matrix_db(left->row[left->size-1], left->col[left->size-1],
	                        left->v[left->size-1], index, parent) );
	/* rightの左端の p_file, p_offset */
	right->p_file[0]   = left->p_file[left->size];
	right->p_offset[0] = left->p_offset[left->size];
	left->size--;
	RC_TRY( clean_matrix_db(left) );

	return(NORMAL_RC);
}


/* rightからleftへ1個移動する               */
/* indexは leftとrightを指している親のindex */
RC
move_from_right_node_matrix_db (MATRIX_DB *left, MATRIX_DB *right,
                                MATRIX_DB *parent, MATRIX_DB_INFO *mdbi,
                                int index)
{
	int ii1;

	RC_NULL_CHK(left);
	RC_NULL_CHK(right);
	RC_NULL_CHK(parent);
	if( (index < 0) || (index >= parent->size) ) return(ARG_ERROR_RC);

	/* parentから1個 leftの最後に移動 */
	RC_TRY( set_v_matrix_db(parent->row[index], parent->col[index],
	                        parent->v[index], left->size, left) );
	left->size++;
	left->p_file[left->size]   = right->p_file[0];
	left->p_offset[left->size] = right->p_offset[0];

	/* rightの最初を parentに移動 */
	RC_TRY( set_v_matrix_db(right->row[0], right->col[0], right->v[0],
	                        index, parent) );

	/* rightを左にずらす */
	right->size--;
	for(ii1=0; ii1<right->size; ii1++){
		RC_TRY( set_v_matrix_db(right->row[ii1+1], right->col[ii1+1],
		                        right->v[ii1+1], ii1, right) );
		right->p_file[ii1]   = right->p_file[ii1+1];
		right->p_offset[ii1] = right->p_offset[ii1+1];
	}
	right->p_file[ii1]   = right->p_file[ii1+1];
	right->p_offset[ii1] = right->p_offset[ii1+1];
	RC_TRY( clean_matrix_db(right) );

	return(NORMAL_RC);
}


/* 親のindex番目を挟んで leftに rightをマージする */
RC
marge_node_matrix_db (MATRIX_DB *left, MATRIX_DB *right, MATRIX_DB *parent,
                      MATRIX_DB_INFO *mdbi, int index)
{
	int ii1;

	RC_NULL_CHK(left);
	RC_NULL_CHK(right);
	RC_NULL_CHK(parent);
	if( (index < 0) || (index >= parent->size) ) return(ARG_ERROR_RC);

	if( (right->size + 1 + left->size) > MATRIX_DB_SIZE)
		return(ARG_ERROR_RC);

	/* parent の index番目を leftに移動 */
	RC_TRY( set_v_matrix_db(parent->row[index], parent->col[index],
	                        parent->v[index], left->size, left) );
	left->size++;

	/* parentのindexより右を左に1ずらす */
	for(ii1=index+1; ii1<parent->size; ii1++){
		RC_TRY( set_v_matrix_db(parent->row[ii1], parent->col[ii1],
		                        parent->v[ii1], ii1-1, parent) );
		parent->p_file[ii1]   = parent->p_file[ii1+1];
		parent->p_offset[ii1] = parent->p_offset[ii1+1];
	}
	parent->size--;
	RC_TRY( clean_matrix_db(parent) );

	/* right をleftに移動 */
	for(ii1=0; ii1<right->size; ii1++){
		RC_TRY( set_v_matrix_db(right->row[ii1], right->col[ii1],
		                        right->v[ii1], left->size, left) );
		left->p_file[left->size]   = right->p_file[ii1];
		left->p_offset[left->size] = right->p_offset[ii1];
		left->size++;
	}
	left->p_file[left->size]   = right->p_file[ii1];
	left->p_offset[left->size] = right->p_offset[ii1];

	right->size = 0; /* right は上位の関数で放棄される */

	if( is_root_matrix_db(parent) ){
		if(parent->size == 0){
			/* leftが rootになる */
			RC_TRY( set_root_matrix_db(left, mdbi) );
			parent->p_file_parent = -1;
			parent->p_offset_parent = -1;
		}
	}else{
		if(parent->size < MATRIX_DB_SIZE/2){
			RC_TRY( balance_node_matrix_db(parent, mdbi) );
		}
	}

	return(NORMAL_RC);
}


/* mdb の index番目を消す */
RC
delete_v_leaf_matrix_db (int index, MATRIX_DB *mdb, MATRIX_DB_INFO *mdbi)
{
	int ii1;
	int row0, col0;

	RC_NULL_CHK(mdb);
	RC_NULL_CHK(mdbi);
	if(mdb->size  <= 0) return(UNDERFLOW_ERROR_RC);
	if( (index  < 0) || (index >= mdb->size) ) return(ARG_ERROR_RC);

	row0 = mdb->row[0];
	col0 = mdb->col[0];
	for(ii1=index+1; ii1<mdb->size; ii1++){
		RC_TRY( set_v_matrix_db(mdb->row[ii1], mdb->col[ii1],
		                        mdb->v[ii1], ii1-1, mdb) );
	}
	(mdb->size)--;
	RC_TRY( clean_matrix_db(mdb) );

	if(mdb->size == 0){
		mdb->key_row = row0;
		mdb->key_col = col0;
	}

	/* 葉のバランスを取る */
	if( !is_root_matrix_db(mdb) && (mdb->size < MATRIX_DB_SIZE/2) ){
		RC_TRY( balance_leaf_matrix_db(mdb, mdbi));
	}

	return(NORMAL_RC);
}


/* 葉のバランスを取る */
RC
balance_leaf_matrix_db (MATRIX_DB *mdb, MATRIX_DB_INFO *mdbi)
{
	int ii1;
	int pindex; /* mdbを指している親の index */
	int p_file;
	WRAP64_INT p_offset;
	MATRIX_DB *parent, *right, *left;
	int p_index;

	RC_NULL_CHK(mdb);
	RC_NULL_CHK(mdbi);

	if( is_root_matrix_db(mdb) ) return(ARG_ERROR_RC);
	if( is_node_matrix_db(mdb) ) return(ARG_ERROR_RC);

	RC_TRY( search_parent_matrix_db(mdbi, mdb, &p_index) );
	if( (p_index < 0) || (p_index >= MATRIX_DB_ROUTE_SIZE) )
		return(UNKNOWN_ERROR_RC);
	RC_TRY( write_matrix_db_start(mdbi->route_file[p_index],
	                              mdbi->route_offset[p_index], &parent,mdbi) );
	/* 念のためチェック: MATRIX_DB_SIZEが小さい時に稀にに問題になる */
#if 0
	if( !is_parent_matrix_db(mdb, parent) ){
		RC_TRY( write_matrix_db_end(parent) );
		mdbi->route_size = 0;
		RC_TRY( search_parent_matrix_db(mdbi, mdb, &p_index) );
		if( (p_index < 0) || (p_index >= MATRIX_DB_ROUTE_SIZE) )
			return(UNKNOWN_ERROR_RC);
		RC_TRY( write_matrix_db_start(mdbi->route_file[p_index],
		                              mdbi->route_offset[p_index],
		                              &parent, mdbi) );
	}
#endif

	/* mdbを指す親の p_file[index] のindexを探す．総当り */
	for(ii1=0, pindex=-1; ii1<parent->size+1; ii1++){
		if( (mdb->p_file_self   == parent->p_file[ii1])
		 && (mdb->p_offset_self == parent->p_offset[ii1]) ){
			pindex = ii1;
			break;
		}
	}
	if(pindex < 0) return(UNKNOWN_ERROR_RC);

	if(parent->size > pindex){
		RC_TRY( write_matrix_db_start(parent->p_file[pindex+1],
		                              parent->p_offset[pindex+1],
		                              &right, mdbi) );
		if((mdb->size + 1 + right->size) <= MATRIX_DB_SIZE){
			/* 右を mdbへマージする */
			RC_TRY( marge_leaf_matrix_db(mdb, right, parent, mdbi, pindex) );
			p_file   = right->p_file_self;
			p_offset = right->p_offset_self;
			RC_TRY( write_matrix_db_end(right) );
			RC_TRY( delete_matrix_db(mdbi, p_file, p_offset) );
		}else if(right->size > MATRIX_DB_SIZE/2){
			/* 右からmdbへ1個移動 */
			RC_TRY( move_from_right_leaf_matrix_db(mdb, right, parent,
			                                       pindex) );
			RC_TRY( write_matrix_db_end(right) );
		}else{
			return(UNKNOWN_ERROR_RC);
		}
	}else if(pindex > 0){
		RC_TRY( write_matrix_db_start(parent->p_file[pindex-1],
		                              parent->p_offset[pindex-1],
		                              &left, mdbi) );
		if((mdb->size + 1 + left->size) <= MATRIX_DB_SIZE){
			/* 左へ mdbをマージする */
			RC_TRY( marge_leaf_matrix_db(left, mdb, parent, mdbi, pindex-1) );
		}else if(left->size > MATRIX_DB_SIZE/2){
			/* 左からmdbへ1個移動 */
			RC_TRY( move_from_left_leaf_matrix_db(left, mdb, parent,
			                                      pindex-1) );
		}else{
			return(UNKNOWN_ERROR_RC);
		}
		RC_TRY( write_matrix_db_end(left) );
	}else{
		return(UNKNOWN_ERROR_RC);
	}

	if(parent->size == 0){
		if(parent->use_count > 1){
			/* 上位の関数で delete_matrix_db() する */
			RC_TRY( write_matrix_db_end(parent) );
		}else if(parent->use_count == 1){
			p_file   = parent->p_file_self;
			p_offset = parent->p_offset_self;
			RC_TRY( write_matrix_db_end(parent) );
			RC_TRY( delete_matrix_db(mdbi, p_file, p_offset) );
		}else{
			return(UNKNOWN_ERROR_RC);
		}
	}else{
		RC_TRY( write_matrix_db_end(parent) );
	}

	return(NORMAL_RC);
}


/* leftからrightへ1個移動する               */
/* indexは leftとrightを指している親のindex */
RC
move_from_left_leaf_matrix_db (MATRIX_DB *left, MATRIX_DB *right,
                               MATRIX_DB *parent, int index)
{
	int ii1;

	RC_NULL_CHK(right);
	RC_NULL_CHK(left);
	RC_NULL_CHK(parent);
	if( (index < 0) || (index >= parent->size) ) return(ARG_ERROR_RC);

	/* rightを1個右にずらして，parentから1個移動 */
	for(ii1=right->size; ii1>0; ii1--){
		RC_TRY( set_v_matrix_db(right->row[ii1-1], right->col[ii1-1],
		                        right->v[ii1-1], ii1, right) );
	}
	right->size++;
	RC_TRY( set_v_matrix_db(parent->row[index], parent->col[index],
	                        parent->v[index], 0, right) );

	/* leftの右端を 1個 parentへ移動 */
	RC_TRY( set_v_matrix_db(left->row[left->size-1], left->col[left->size-1],
	                        left->v[left->size-1], index, parent) );
	left->size--;
	RC_TRY( clean_matrix_db(left) );

	return(NORMAL_RC);
}


/* rightから leftへ1個移動する              */
/* indexは leftとrightを指している親のindex */
RC
move_from_right_leaf_matrix_db (MATRIX_DB *left, MATRIX_DB *right,
                                MATRIX_DB *parent, int index)
{
	int ii1;

	RC_NULL_CHK(right);
	RC_NULL_CHK(left);
	RC_NULL_CHK(parent);
	if( (index < 0) || (index >= parent->size) ) return(ARG_ERROR_RC);

	/* parentから1個 leftの最後に移動 */
	RC_TRY( set_v_matrix_db(parent->row[index], parent->col[index],
	                        parent->v[index], left->size, left) );
	left->size++;

	/* rightの最初を parentに移動 */
	RC_TRY( set_v_matrix_db(right->row[0], right->col[0], right->v[0],
	                        index, parent) );
	right->size--;

	/* rightを左にずらす */
	for(ii1=0; ii1<right->size; ii1++){
		RC_TRY( set_v_matrix_db(right->row[ii1+1], right->col[ii1+1],
		                        right->v[ii1+1], ii1, right) );
	}
	RC_TRY( clean_matrix_db(right) );

	return(NORMAL_RC);
}


/* 親のindex番目を挟んで leftにrightをマージする */
RC
marge_leaf_matrix_db (MATRIX_DB *left, MATRIX_DB *right, MATRIX_DB *parent,
                      MATRIX_DB_INFO *mdbi, int index)
{
	int ii1;

	RC_NULL_CHK(right);
	RC_NULL_CHK(left);
	RC_NULL_CHK(parent);
	RC_NULL_CHK(mdbi);
	if( (index < 0) || (index >= parent->size) ) return(ARG_ERROR_RC);
	if( (right->size + 1 + left->size) > MATRIX_DB_SIZE) return(ARG_ERROR_RC);

	/* parent の index番目を leftに移動 */
	RC_TRY( set_v_matrix_db(parent->row[index], parent->col[index],
	                        parent->v[index], left->size, left) );
	left->size++;

	/* parentを空白分をずらす */
	for(ii1=index+1; ii1<parent->size; ii1++){
		RC_TRY( set_v_matrix_db(parent->row[ii1], parent->col[ii1],
		                        parent->v[ii1], ii1-1, parent) );
		parent->p_file[ii1]   = parent->p_file[ii1+1];
		parent->p_offset[ii1] = parent->p_offset[ii1+1];
	}
	parent->size--;
	RC_TRY( clean_matrix_db(parent) );

	/* right をleftに移動 */
	for(ii1=0; ii1<right->size; ii1++){
		RC_TRY( set_v_matrix_db(right->row[ii1], right->col[ii1],
		                        right->v[ii1], left->size, left) );
		left->size++;
	}
	right->size = 0;
	
	if( is_root_matrix_db(parent) ){
		if(parent->size == 0){
			/* leftが rootになる */
			RC_TRY( set_root_matrix_db(left, mdbi) );
			parent->p_file_parent = -1;
			parent->p_offset_parent = -1;
		}
	}else{
		if(parent->size < MATRIX_DB_SIZE/2){
			RC_TRY( balance_node_matrix_db(parent, mdbi) );
		}
	}

	return(NORMAL_RC);
}


/*
 * row, col が一致する v が mdb の中にあるか探索する．
 * row, col が一致すれば，found_flagが ONになり，そのindexを返す．
 * row, col が一致するものが無ければ，found_flagは OFFで，indexにはrow, col
 * を挿入すればよい位置(row,colのよりも一つ大きな値のindex)が返る．
 * また,found_flagは OFFならば，次の探査候補となる MATRIX_DBの p_file,
 * p_offsetがセットされる．
 * 通常のバイナリサーチと異なり，一致するものを探すだけでなく，row-colの
 * を挿入する配列のindexを探す場合にも対応しているため少々複雑になっている．
 */
RC
search_matrix_db (int row, int col, MATRIX_DB *mdb, int *found_flag, 
                  int *index, int *p_file, WRAP64_INT *p_offset)
{
	int comp;
	int low;
	int high;
	int middle;

	RC_NULL_CHK(found_flag);
	RC_NULL_CHK(mdb);
	RC_NULL_CHK(index);
	RC_NULL_CHK(p_file);
	RC_NULL_CHK(p_offset);
	if( (row < 0) || (col < 0) ) return(ARG_ERROR_RC);

	low = 0;
	high = mdb->size-1;

	*found_flag = 0;
	*index = -1;
	*p_file = 0;
	*p_offset = 0;

	if(mdb->size < 1){
		*index = 0; /* 必須!! */
		return(NORMAL_RC);
	}

#if 0
    /* 先に境界(0,max)のindexチェックをしてもスピードは変わらない．*/
	/* min */
	comp = compar_row_col(mdb->row[0], mdb->col[0], row, col);
	if(comp > 0){
		*index = 0;
		*p_file = mdb->p_file[0];
		*p_offset = mdb->p_offset[0];
		return(NORMAL_RC);
	}else if(comp == 0){
		*found_flag = 1;
		*index = 0;
		return(NORMAL_RC);
	}

	/* max */
	comp = compar_row_col(mdb->row[mdb->size-1], mdb->col[mdb->size-1],
	                      row, col);
	if(comp < 0){
		/*
		log_printf(5, "mdb->row[0]=%d, mdb->col[0]=%d, row=%d, col=%d\n",
		           mdb->row[mdb->size-1], mdb->col[mdb->size-1],row, col);
		print_matrix_db_debug(stderr,  mdb);
		*/
		*index = mdb->size;
		*p_file = mdb->p_file[*index];
		*p_offset = mdb->p_offset[*index];
		return(NORMAL_RC);
	}else if(comp == 0){
		*found_flag = 1;
		*index = mdb->size-1;
		return(NORMAL_RC);
	}
#endif

	while(low <= high){
		middle = (low + high)/2;
		comp = compar_row_col(mdb->row[middle], mdb->col[middle], row, col);
		if(comp == 0){
			*found_flag = 1;
			*index = middle;
			*p_file = 0;
			*p_offset = 0;
			break;
		}else if(comp > 0){
			high = middle - 1;
			if(high <= low){
				if(high < 0) high = 0;
				comp = compar_row_col(mdb->row[high],
				                      mdb->col[high], row, col);
				if(comp == 0){
					*found_flag = 1;
					*index = high;
					break;
				}else if(comp > 0){
					*index = high;
				}else{
					if(mdb->size < high+1){
						return(UNKNOWN_ERROR_RC);
						/* *index = high; */
					}else{
						*index = high+1;
					}
				}
				*p_file = mdb->p_file[*index];
				*p_offset = mdb->p_offset[*index];
				break;
			}
		}else{
			low = middle + 1;
			if(high <= low){
				comp = compar_row_col(mdb->row[low], mdb->col[low], row, col);
				if(comp == 0){
					*found_flag = 1;
					*index = high;
					break;
				}else if(comp > 0){
					*index = low;
				}else{
					if( (mdb->size == MATRIX_DB_SIZE) &&(mdb->size <= low+1) ){
						*index = low;
						*p_file = mdb->p_file[*index+1];
						*p_offset = mdb->p_offset[*index+1];
						break;
					}else if(mdb->size < low+1){
						*index = low;
					}else{
						*index = low+1;
					}
				}
				*p_file = mdb->p_file[*index];
				*p_offset = mdb->p_offset[*index];
				break;
			}
		}
	}

	return(NORMAL_RC);
}


/* rowでsort して colでsort */
static int
compar_row_col (int row1, int col1, int row2, int col2)
{
	if(row1 < row2){
		return(-1);
	}else if(row1 == row2){
		if(col1 < col2){
			return(-1);
		}else if(col1 == col2){
			return(0);
		}else{ /* col1 > col2 */
			return(1);
		}
	}else{ /* row1 > row2 */
		return(1);
	}
}

#if 0
/* colでsort して rowでsort */
static int
compar_row_col (int row1, int col1, int row2, int col2)
{
	if(col1 < col2){
		return(-1);
	}else if(col1 == col2){
		if(row1 < row2){
			return(-1);
		}else if(row1 == row2){
			return(0);
		}else{ /* row1 > row2 */
			return(1);
		}
	}else{ /* col1 > col2 */
		return(1);
	}
}
#endif


/*
 * cache 内で db_file, db_offset のデータを探してそのポインタを返却
 * 無かったら NULL を返却
 * アクセスカウンタの調整もここで行う．
 * inc_access : 0 ならアクセスカウンタに小さな値をセット
 *              1 ならアクセスカウンタをリセット
 */
static MATRIX_DB *
search_cache_ptr (MATRIX_DB_INFO *mdbi, int db_file, WRAP64_INT db_offset,
                  int inc_access)
{
	int ii1, ii2;
	int crow = (int)(CACHE_ROW(db_offset, mdbi->cache_divisor));

	for(ii1=0; ii1<MATRIX_DB_ASSOCIATE; ii1++){
		if(mdbi->cache[crow][ii1].use_count < 0){
			return(NULL);
		}else if( (mdbi->cache[crow][ii1].p_file_self == db_file)
		       && (mdbi->cache[crow][ii1].p_offset_self == db_offset) ){
			if(inc_access){
				if(mdbi->access_count[crow][ii1] + MAX_ACCESS_COUNT/32
				   >= MAX_ACCESS_COUNT){
					for(ii2=0; ii2<MATRIX_DB_ASSOCIATE; ii2++){
						(mdbi->access_count[crow][ii2]) /= 2;
					}
				}
				mdbi->access_count[crow][ii1] += MAX_ACCESS_COUNT/32;
			}else{
				mdbi->access_count[crow][ii1] = 0;
			}

			return( &(mdbi->cache[crow][ii1]) );
		}
	}

	return(NULL);
}


/*
 * cache[crow][] 内で未使用領域を探してそのポインタを返却
 * 未使用領域が無ければ，できるだけアクセス数の少なくアクセス中でない
 * データのポインタを返却
 * その際，dirty_flag が ON ならディスクに書き出して OFF にする．
 * 全てのデータがアクセス中なら NULL を返却
 * アクセスカウンタのリセットもここで行う．
 */ 
static MATRIX_DB *
unused_cache_ptr (MATRIX_DB_INFO *mdbi, int crow)
{
	int ii1;
	int min_count = INT_MAX;
	int min_col = -1;
	MATRIX_DB *ret;

	for(ii1=0; ii1<MATRIX_DB_ASSOCIATE; ii1++){
		if(mdbi->cache[crow][ii1].use_count < 0){
			/* 未使用領域 */
			min_count = -1;
			min_col = ii1;
			break;
		}else if(mdbi->cache[crow][ii1].use_count == 0){
			/* 使用されているがアクセス中ではない */
			int access = mdbi->access_count[crow][ii1];
			if(mdbi->cache[crow][ii1].dirty_flag){
				/* dirty_flag が ON の場合は，ディスク書き込みが発生するので */
				/* アクセス数を見かけ上大きくして，破棄されにくくする */
				access *= 2;
			}
			if(access < min_count){
				min_count = access;
				min_col = ii1;
			}
		}
	}

	if(min_col < 0) return(NULL);  /* 全部アクセス中 */

	ret = &(mdbi->cache[crow][min_col]);
	if(ret->dirty_flag){
		if(put_matrix_db(mdbi, ret->p_file_self,
		                 ret->p_offset_self, ret) != NORMAL_RC) return(NULL);
		ret->dirty_flag = 0;
	}

	/* 新規登録されたデータがすぐに破棄されないように，*/
	/* アクセスカウンタは最大値にする */
	mdbi->access_count[crow][min_col] = MAX_ACCESS_COUNT;

	return(ret);
}


RC
resize_cache_matrix_db (MATRIX_DB_INFO *mdbi, WRAP64_INT size)
{
	int ii1, ii2;
	int p_file_root;
	WRAP64_INT  p_offset_root;
	WRAP64_INT mdb_asso_size = sizeof(MATRIX_DB)*MATRIX_DB_ASSOCIATE;

	RC_NULL_CHK( mdbi);
	if(size < mdb_asso_size) return(ARG_ERROR_RC);

	/* root以外の use_countが無いかチェック */
	for(ii1=0; ii1<(mdbi->cache_divisor); ii1++){
		for(ii2=0; ii2<MATRIX_DB_ASSOCIATE; ii2++){
			if(mdbi->cache[ii1][ii2].use_count > 0){
				if( is_root_matrix_db(&(mdbi->cache[ii1][ii2])) ) continue;
				else return(ARG_ERROR_RC);
			}
		}
	}

	RC_TRY( flush_cache_matrix_db(mdbi) );
	p_file_root   = mdbi->root->p_file_self;
	p_offset_root = mdbi->root->p_offset_self;

	RC_TRY( mm_free(mdbi->cache) );
	RC_TRY( mm_free(mdbi->access_count) );

	mdbi->cache_divisor = (int)(size/mdb_asso_size);
	mdbi->cache = (MATRIX_DB (*)[MATRIX_DB_ASSOCIATE])mm_alloc(
	                         (size_t)(mdb_asso_size *(mdbi->cache_divisor)));
	RC_NULL_CHK( mdbi->cache );
	mdbi->access_count = (int (*)[MATRIX_DB_ASSOCIATE])mm_alloc(
	                     sizeof(int)*MATRIX_DB_ASSOCIATE*(mdbi->cache_divisor));
	RC_NULL_CHK( mdbi->access_count );

	for(ii1=0; ii1<(mdbi->cache_divisor); ii1++){
		for(ii2=0; ii2<MATRIX_DB_ASSOCIATE; ii2++){
			INIT_MATRIX_DB(mdbi->cache[ii1][ii2]);
			mdbi->access_count[ii1][ii2] = -1;
		}
	}

	RC_TRY( read_matrix_db_start(p_file_root, p_offset_root,
	                             &(mdbi->root), mdbi) );
	mdbi->root->p_file_parent = 0;
	mdbi->root->p_offset_parent = 0;

	return(NORMAL_RC);
}


/* キャッシュの dirty_flag = ON の MATRIX_DBをディスクに書き出す */
RC
flush_cache_matrix_db (MATRIX_DB_INFO *mdbi)
{
	int ii1, ii2;

	RC_NULL_CHK( mdbi );

	for(ii1=0; ii1<mdbi->cache_divisor; ii1++){
		for(ii2=0; ii2<MATRIX_DB_ASSOCIATE; ii2++){
			MATRIX_DB *cache = &(mdbi->cache[ii1][ii2]);
			if(cache->use_count == -1) continue;
			if(cache->dirty_flag == 1){
				cache->dirty_flag = 0;
				RC_TRY( put_matrix_db(mdbi, cache->p_file_self,
				                      cache->p_offset_self, cache) );
			}
		}
	}

	return(NORMAL_RC);
}


/* キャッシュとファイル上のデータ(db_file/db_offset)を削除（再利用可能に） */
RC
delete_matrix_db (MATRIX_DB_INFO *mdbi, int db_file, WRAP64_INT db_offset)
{
	MATRIX_DB tmp_db;
	MATRIX_DB *ptr;

	RC_NULL_CHK( mdbi );
	if((db_file <= 0)||(db_file > MATRIX_DB_SF_SIZE)) return(ARG_ERROR_RC);

	/* キャッシュ内にデータがあれば初期化し，p_offset_parent */
	/* に直近に削除されたデータのオフセットを書き込み        */
	ptr = search_cache_ptr(mdbi, db_file, db_offset, 0);
	if(ptr != NULL){
		if(ptr->use_count > 0){
			log_printf(5, "[%d] use_count > 0\n", __LINE__);
			return(UNKNOWN_ERROR_RC);
		}
		INIT_MATRIX_DB(*ptr);
		ptr->use_count = 0;
		ptr->p_offset_parent = mdbi->blank_p_offset[db_file - 1];
		ptr->p_file_self = db_file;
		ptr->p_offset_self = db_offset;
		ptr->dirty_flag = 1;

		/* blkank_p_offset を更新 */
		mdbi->blank_p_offset[db_file - 1] = db_offset;

		return(NORMAL_RC);
	}

	/* p_offset_parent に直近に削除されたデータのオフセットを書き込み */
	INIT_MATRIX_DB(tmp_db);
	tmp_db.p_offset_parent = mdbi->blank_p_offset[db_file - 1];
	tmp_db.p_file_self = db_file;
	tmp_db.p_offset_self = db_offset;
	RC_TRY( put_matrix_db(mdbi, db_file, db_offset, &tmp_db) );

	/* blkank_p_offset を更新 */
	mdbi->blank_p_offset[db_file - 1] = db_offset;

	return(NORMAL_RC);
}


RC
new_matrix_db_start (MATRIX_DB **mdb, MATRIX_DB_INFO *mdbi)
{
	int ii1;
	int db_file = 0;
	WRAP64_INT db_offset = 0;
	int reuse_flag = 0;

	RC_NULL_CHK( mdb );
	RC_NULL_CHK( mdbi );
	*mdb = NULL;

	/* 直近に削除された領域があれば，db_file/db_offset に取得 */
	for(ii1=0; ii1<MATRIX_DB_SF_SIZE; ii1++){
		if(mdbi->blank_p_offset[ii1] < 0) continue;

		db_file = ii1 + 1;
		db_offset = mdbi->blank_p_offset[ii1];
		reuse_flag = 1;

		/* 当該データはキャッシュ上にあるかも */
		*mdb = search_cache_ptr(mdbi, db_file, db_offset, 1);
		if(*mdb != NULL){
			mdbi->blank_p_offset[ii1] = (*mdb)->p_offset_parent;
			break;
		}

		/* 先に削除された領域のオフセットを読み込んで blank_p_offset へ */
		*mdb = unused_cache_ptr(mdbi,
		                 (size_t)(CACHE_ROW(db_offset, mdbi->cache_divisor)) );
		RC_NULL_CHK( *mdb );
		RC_TRY( get_matrix_db(mdbi, db_file, db_offset, *mdb) );

		mdbi->blank_p_offset[ii1] = (*mdb)->p_offset_parent;
		break;
	}

	if(db_file <= 0){
		/* 新規に作成して，db_file/db_offset をセット */

		/* ここで，スクラッチファイルを選択するが，*/
		/* 今のところスクラッチファイルが一つなので db_file = 1 に決め打ち */
		db_file = 1;
		db_offset = mdbi->last_p_offset[0];
		mdbi->last_p_offset[0] += sizeof(MATRIX_DB);

		*mdb = unused_cache_ptr(mdbi,
		              (size_t)(CACHE_ROW(db_offset, mdbi->cache_divisor)));
		RC_NULL_CHK( *mdb );
	}

	/* キャッシュ内にデータを作成 */
	RC_NULL_CHK( *mdb );    /* 念のため */
	INIT_MATRIX_DB(**mdb);
	(*mdb)->use_count = 1;
	(*mdb)->dirty_flag = 1;
	(*mdb)->p_file_self = db_file;
	(*mdb)->p_offset_self = db_offset;
	if(reuse_flag) return(NORMAL_RC);

	/* キャッシュのデータをファイルに書き込み */
	RC_TRY( put_matrix_db(mdbi, db_file, db_offset, *mdb) );

	return(NORMAL_RC);
}


RC
read_matrix_db_start (int p_file, WRAP64_INT p_offset, MATRIX_DB **mdb,
                      MATRIX_DB_INFO *mdbi)
{
	RC_NULL_CHK( mdb );
	RC_NULL_CHK( mdbi );
	if( (p_file <= 0) || (p_file > MATRIX_DB_SF_SIZE) ) return(ARG_ERROR_RC);
	if(p_offset < 0)  return(ARG_ERROR_RC);

	/* キャッシュ内を探す */
	*mdb = search_cache_ptr(mdbi, p_file, p_offset, 1);
	if(*mdb != NULL){
		(*mdb)->use_count++;
		return(NORMAL_RC);
	}

	/* キャッシュの空き領域に読み込む */
	*mdb = unused_cache_ptr(mdbi,
	              (int)(CACHE_ROW(p_offset, mdbi->cache_divisor)) );
	RC_NULL_CHK( *mdb );
	RC_TRY( get_matrix_db(mdbi, p_file, p_offset, *mdb) );

	if((*mdb)->use_count < 0){
		return(UNKNOWN_ERROR_RC);   /* 使われていない？ */
	}
	(*mdb)->use_count = 1;
	(*mdb)->dirty_flag = 0;

	return(NORMAL_RC);
}


RC
write_matrix_db_start (int p_file, WRAP64_INT p_offset, MATRIX_DB **mdb,
                       MATRIX_DB_INFO *mdbi)
{
	RC_TRY( read_matrix_db_start(p_file, p_offset, mdb, mdbi) );
	(*mdb)->dirty_flag = 1;

	return(NORMAL_RC);
}


RC
new_matrix_db_end (MATRIX_DB *mdb)
{
	RC_NULL_CHK( mdb );
	if(mdb->use_count <= 0) return(ARG_ERROR_RC);
	mdb->use_count--;

	return(NORMAL_RC);
}


RC
read_matrix_db_end (MATRIX_DB *mdb)
{
	RC_NULL_CHK( mdb );
	if(mdb->use_count <= 0) return(ARG_ERROR_RC);
	mdb->use_count--;

	return(NORMAL_RC);
}


RC
write_matrix_db_end (MATRIX_DB *mdb)
{
	RC_NULL_CHK( mdb );
	if(mdb->use_count <= 0) return(ARG_ERROR_RC);
	mdb->use_count--;

	return(NORMAL_RC);
}


RC
put_matrix_db (MATRIX_DB_INFO *mdbi, int p_file, WRAP64_INT p_offset,
               MATRIX_DB *ptr)
{
	SCRATCH_FILE *sf;
	FILE *fp;

	RC_NULL_CHK( mdbi );
	RC_NULL_CHK( ptr );
	if( (p_file <= 0) || (p_file > MATRIX_DB_SF_SIZE) ) return(ARG_ERROR_RC);
	if(p_offset < 0)  return(ARG_ERROR_RC);

	ptr->p_file_self = p_file;
	ptr->p_offset_self = p_offset;

	RC_NULL_CHK( sf = mdbi->scratch_array[p_file - 1] );
	RC_NULL_CHK( fp = sf->fp );

	if(p_offset != mdbi->current_p_offset[p_file - 1]){
		/* 現在の読み書き位置と異なる場合だけシークする */
		/* wrap64_fseek() 内でバッファがフラッシュされるのを極力避けるため */
		RC_TRY( wrap64_fseek(fp, p_offset, SEEK_SET) );
	}
	if(fwrite(ptr, sizeof(MATRIX_DB), 1, fp) != 1){
		return(WRITE_ERROR_RC);    /* ディスク書き込み失敗 */
	}
	mdbi->current_p_offset[p_file - 1] = p_offset + sizeof(MATRIX_DB);

	return(NORMAL_RC);
}


RC
get_matrix_db (MATRIX_DB_INFO *mdbi, int p_file, WRAP64_INT p_offset,
               MATRIX_DB *ptr)
{
	SCRATCH_FILE *sf;
	FILE *fp;

	RC_NULL_CHK( mdbi );
	RC_NULL_CHK( ptr );
	if( (p_file <= 0) || (p_file > MATRIX_DB_SF_SIZE) ) return(ARG_ERROR_RC);
	if(p_offset < 0)  return(ARG_ERROR_RC);

	RC_NULL_CHK( sf = mdbi->scratch_array[p_file - 1] );
	RC_NULL_CHK( fp = sf->fp );
	if(p_offset != mdbi->current_p_offset[p_file - 1]){
		/* 現在の読み書き位置と異なる場合だけシークする */
		/* wrap64_fseek() 内でバッファがフラッシュされるのを極力避けるため */
		RC_TRY( wrap64_fseek(fp, p_offset, SEEK_SET) );
	}
	if(fread(ptr, sizeof(MATRIX_DB), 1, fp) != 1){
		return(READ_ERROR_RC);    /* ディスク読み込み失敗 */
	}
	mdbi->current_p_offset[p_file - 1] = p_offset + sizeof(MATRIX_DB);

	return(NORMAL_RC);
}


/* 初期化されたMATRIX_DBならTRUE */
int
is_init_matrix_db (MATRIX_DB *mdb)
{
	return(mdb->p_file_self == 0);
}

int
is_leaf_matrix_db (MATRIX_DB *mdb)
{
	/* 葉には子がいない */
	return( (!mdb->p_file[0]) && !(is_init_matrix_db(mdb)) );
}

int
is_root_matrix_db (MATRIX_DB *mdb)
{
	/* rootには親がいない．子がいるとは限らない． */
	return( (mdb->p_file_parent == 0) && !(is_init_matrix_db(mdb)) );
}

int
is_node_matrix_db (MATRIX_DB *mdb)
{
	/* 少なくとも葉ではない */
	return( !(is_leaf_matrix_db(mdb)) );
}

int
is_parent_matrix_db (MATRIX_DB *mdb, MATRIX_DB *parent)
{
	int ii1;

	if(mdb == NULL) return(0);
	if(parent == NULL) return(0);

	for(ii1=0; ii1<parent->size+1; ii1++){
		if(mdb->p_file_self == parent->p_file[ii1]) {
			if(mdb->p_offset_self == parent->p_offset[ii1]){
				return(1);
			}
		}
	}
	return(0);
}

int
is_range_matrix_db (MATRIX_DB *mdb, int index, int row, int col_min,
                    int col_max)
{
	if(mdb == NULL) return(0);
	if( (index < 0) || (index >= mdb->size) ) return(0);
	if( (row < 0) || (col_min < 0) || (col_max < 0) ) return(0);

	if( (mdb->row[index] == row) && (mdb->col[index] >= col_min)
	                             && (mdb->col[index] <  col_max) ){
		return(1);
	}

	return(0);
}


RC
set_v_matrix_db (int row, int col, double v[][3], int index, MATRIX_DB *mdb)
{
	RC_NULL_CHK(mdb);
	RC_NULL_CHK(v);
	if(index < 0) return(ARG_ERROR_RC);
	/* DO NOT SET
	if(index > mdb->size) return(ARG_ERROR_RC);
	 */
	if( (row < 0) || (col < 0) ) return(ARG_ERROR_RC);

	mdb->row[index] = row;
	mdb->col[index] = col;
	RC_TRY( copy_matrix33(v, mdb->v[index]) );

	return(NORMAL_RC);
}


RC
set_root_matrix_db (MATRIX_DB *root, MATRIX_DB_INFO *mdbi)
{
	RC_NULL_CHK(root);
	RC_NULL_CHK(mdbi);

	RC_TRY( read_matrix_db_end(mdbi->root) );
	root->p_file_parent = 0;
	root->p_offset_parent = 0;
	mdbi->root = root;
	RC_TRY( read_matrix_db_start(root->p_file_self, root->p_offset_self,
	                             &(mdbi->root), mdbi) );
	return(NORMAL_RC);
}


/* mdb->size 以上の配列を初期化 */
RC
clean_matrix_db (MATRIX_DB *mdb)
{
	int ii1;

	RC_NULL_CHK(mdb);

	for(ii1=mdb->size+1; ii1<MATRIX_DB_SIZE+1; ii1++){
		mdb->p_file[ii1] = 0;
		mdb->p_offset[ii1] = 0;
	}
	if(mdb->size == 0){
		mdb->p_file[0] = 0;
		mdb->p_offset[0] = 0;
	}

	for(ii1=mdb->size; ii1<MATRIX_DB_SIZE; ii1++){
		mdb->row[ii1] = 0;
		mdb->col[ii1] = 0;
		init_matrix33(mdb->v[ii1]);
	}

	return(NORMAL_RC);
}


RC
init_matrix_db (MATRIX_DB *mdb)
{
	int ii1;

	RC_NULL_CHK(mdb);

	mdb->size = 0;
	mdb->dirty_flag = 0;
	mdb->use_count = -1;
	mdb->key_row = -1;
	mdb->key_col = -1;
	mdb->p_file_parent = -1;
	mdb->p_offset_parent = -1;
	mdb->p_file_self = 0;
	mdb->p_offset_self = 0;

	for(ii1=0; ii1<MATRIX_DB_HEADER_SIZE; ii1++){
		mdb->header[ii1] = 0;
	}

	for(ii1=0; ii1<MATRIX_DB_SIZE+1; ii1++){
		mdb->p_file[ii1] = 0;
		mdb->p_offset[ii1] = 0; /* WRAP64_INTの初期化. "型"注意!! */
	}

	for(ii1=0; ii1<MATRIX_DB_SIZE; ii1++){
		mdb->row[ii1] = 0;
		mdb->col[ii1] = 0;
		init_matrix33(mdb->v[ii1]);
	}

	return(NORMAL_RC);
}


RC
print_matrix_db_cache_info (FILE *fp, MATRIX_DB_INFO *mdbi)
{
	int ii1, ii2;

	RC_NULL_CHK(mdbi);

	for(ii1=0; ii1<mdbi->cache_divisor; ii1++){
		fprintf(fp, "cache_divisor : %2d\n", ii1);
		fprintf(fp, "asso dirty use     access file    offset type\n");
		for(ii2=0; ii2<MATRIX_DB_ASSOCIATE; ii2++){
			MATRIX_DB *cache = &mdbi->cache[ii1][ii2];
			if(cache->use_count == -1) continue;
#ifdef WIN32
			fprintf(fp, "[%2d]     %d  %2d %10d    %d %9I64d",
			             ii2,
			             cache->dirty_flag, cache->use_count,
			             mdbi->access_count[ii1][ii2],
			             cache->p_file_self, cache->p_offset_self);
#else /* WIN32 */
  #ifdef _LARGEFILE_SOURCE
			fprintf(fp, "[%2d]     %d  %2d %10d    %d %9lld",
			             ii2,
			             cache->dirty_flag, cache->use_count,
			             mdbi->access_count[ii1][ii2],
			             cache->p_file_self, cache->p_offset_self);
  #else
			fprintf(fp, "[%2d]     %d  %2d %10d    %d %9ld",
			             ii2,
			             cache->dirty_flag, cache->use_count,
			             mdbi->access_count[ii1][ii2],
			             cache->p_file_self, cache->p_offset_self);
  #endif
#endif  /* WIN32 */
			if( is_root_matrix_db(cache) ){
				fprintf(fp, " root\n");
			}else if( is_node_matrix_db(cache) ){
				fprintf(fp, " node\n");
			}else if( is_leaf_matrix_db(cache) ){
				fprintf(fp, " leaf\n");
			}else{
				fprintf(fp, " unknown\n");
			}
		}
	}

	return(NORMAL_RC);
}


/* Treeが大きい場合，かなりのメモリが必要 */
RC
print_stratum_info_matrix_db (FILE *fp, MATRIX_DB_INFO *mdbi)
{
	int ii1, ii2, ii3;
	int end_flag = 0;
	MATRIX_DB *mdb;
	int size;
	int *p_file;
	WRAP64_INT *p_offset;
	int size_n; /* next */
	int *p_file_n = NULL;
	WRAP64_INT *p_offset_n = NULL;

	RC_NULL_CHK(mdbi);

	/* RootNode を階層1としてセット */
	size = 1;
	RC_NULL_CHK( p_file = (int *)mm_alloc(sizeof(int)*(size)) );
	RC_NULL_CHK( p_offset = (WRAP64_INT *)mm_alloc(
	                                      sizeof(WRAP64_INT)*(size)) );
	p_file[0]   = mdbi->root->p_file_self;
	p_offset[0] = mdbi->root->p_offset_self;

	for(ii1=size, size_n=size;; ii1++){
		/* 下の階層のデータを入れる配列 確保 */
		size_n = (MATRIX_DB_SIZE+1)*size_n;
		RC_NULL_CHK( p_file_n = (int *)mm_alloc(sizeof(int)*size_n) );
		p_offset_n = (WRAP64_INT *)mm_alloc(sizeof(WRAP64_INT)*(size_n));
		RC_NULL_CHK( p_offset_n );
		/* 初期化 */
		for(ii2=0; ii2<size_n; ii2++){
			p_file_n[ii2] = 0;
			p_offset_n[ii2] = 0;
		}

		fprintf(fp, "----- stratum %d -----[%d MATRIX_DB]\n", ii1, size);
		fprintf(fp, "stratum         "
		            "          self        parent\n");
		fprintf(fp, "     index  "
		            "size    offset/file   offset/file    "
		            "row(min-max)    col(min-max)\n");

		for(ii2=0, end_flag=0; ii2<size; ii2++){
			if(p_file[ii2] != 0){
				RC_TRY( read_matrix_db_start(p_file[ii2], p_offset[ii2],
			                             &mdb, mdbi) );
#ifdef WIN32
				fprintf(fp, "[%d][%5d]   %3d %12I64d/%d %10I64d/%d"
				            " %7d-%7d %7d-%7d\n", ii1, ii2, mdb->size,
				            mdb->p_offset_self,   mdb->p_file_self,
				            mdb->p_offset_parent, mdb->p_file_parent,
				            mdb->row[0], mdb->row[mdb->size-1],
				            mdb->col[0], mdb->col[mdb->size-1]);
#else /* WIN32 */
  #ifdef _LARGEFILE_SOURCE
				fprintf(fp, "[%d][%5d]   %3d %12lld/%d %10lld/%d"
				            " %7d-%7d %7d-%7d\n", ii1, ii2, mdb->size,
				            mdb->p_offset_self,   mdb->p_file_self,
				            mdb->p_offset_parent, mdb->p_file_parent,
				            mdb->row[0], mdb->row[mdb->size-1],
				            mdb->col[0], mdb->col[mdb->size-1]);
  #else
				fprintf(fp, "[%d][%5d]   %3d %12ld/%d %10ld/%d"
				            " %7d-%7d %7d-%7d\n", ii1, ii2, mdb->size,
				            mdb->p_offset_self,   mdb->p_file_self,
				            mdb->p_offset_parent, mdb->p_file_parent,
				            mdb->row[0], mdb->row[mdb->size-1],
				            mdb->col[0], mdb->col[mdb->size-1]);
  #endif
#endif  /* WIN32 */
				if( is_leaf_matrix_db(mdb) ){
					end_flag = 1;
				}else{
					for(ii3=0; ii3<MATRIX_DB_SIZE+1; ii3++){
						int ii4 = ii2*(MATRIX_DB_SIZE+1)+ii3;
						p_file_n[ii4]   = mdb->p_file[ii3];
						p_offset_n[ii4] = mdb->p_offset[ii3];
					}
				}
				RC_TRY( read_matrix_db_end(mdb) );
			}
		}

		size = size_n;

		RC_TRY( mm_free(p_file) );
		RC_TRY( mm_free(p_offset) );
		p_file   = p_file_n;
		p_offset = p_offset_n;
		if( end_flag == 1){
			RC_TRY( mm_free(p_file_n) );
			RC_TRY( mm_free(p_offset_n) );
			break;
		}
		p_file_n = NULL;
		p_offset_n = NULL;
	}

	return(NORMAL_RC);
}


RC
print_btree_matrix_db (FILE *fp, MATRIX_DB_INFO *mdbi)
{
	int ii1, ii2, ii3;
	int end_flag = 0;
	MATRIX_DB *mdb;
	int size;
	int *p_file;
	WRAP64_INT *p_offset;
	int size_n; /* next */
	int *p_file_n = NULL;
	WRAP64_INT *p_offset_n = NULL;

	RC_NULL_CHK(mdbi);

	/* RootNode を階層1としてセット */
	size = 1;
	RC_NULL_CHK( p_file = (int *)mm_alloc(sizeof(int)*(size)) );
	RC_NULL_CHK( p_offset = (WRAP64_INT *)mm_alloc(
	                                      sizeof(WRAP64_INT)*(size)) );
	p_file[0]   = mdbi->root->p_file_self;
	p_offset[0] = mdbi->root->p_offset_self;

	for(ii1=size, size_n=size;; ii1++){
		/* 下の階層のデータを入れる配列 確保 */
		size_n = (MATRIX_DB_SIZE+1)*size_n;
		RC_NULL_CHK( p_file_n = (int *)mm_alloc(sizeof(int)*size_n) );
		p_offset_n = (WRAP64_INT *)mm_alloc(sizeof(WRAP64_INT)*(size_n));
		RC_NULL_CHK( p_offset_n );
		/* 初期化 */
		for(ii2=0; ii2<size_n; ii2++){
			p_file_n[ii2] = 0;
			p_offset_n[ii2] = 0;
		}

		fprintf(fp, "----- stratum %d -----\n", ii1);
		for(ii2=0, end_flag=0; ii2<size; ii2++){
			if(p_file[ii2] != 0){
				RC_TRY( read_matrix_db_start(p_file[ii2], p_offset[ii2],
				                             &mdb, mdbi) );
				RC_TRY( print_matrix_db_debug(fp, mdb) );
				if( is_leaf_matrix_db(mdb) ){
					end_flag = 1;
				}else{
					for(ii3=0; ii3<MATRIX_DB_SIZE+1; ii3++){
						int ii4 = ii2*(MATRIX_DB_SIZE+1)+ii3;
						p_file_n[ii4]   = mdb->p_file[ii3];
						p_offset_n[ii4] = mdb->p_offset[ii3];
					}
				}
				RC_TRY( read_matrix_db_end(mdb) );
			}
		}

		size = size_n;

		RC_TRY( mm_free(p_file) );
		RC_TRY( mm_free(p_offset) );
		p_file   = p_file_n;
		p_offset = p_offset_n;
		if( end_flag == 1){
			RC_TRY( mm_free(p_file_n) );
			RC_TRY( mm_free(p_offset_n) );
			break;
		}
		p_file_n = NULL;
		p_offset_n = NULL;
	}

	return(NORMAL_RC);
}


/* MATRIX_DB 構造体を定義のまま出力 */
RC
print_matrix_db (FILE *fp, MATRIX_DB *mdb)
{
	int ii1;

	RC_NULL_CHK(mdb);

	RC_TRY( print_matrix_db_info(fp, mdb) );

#if 0
	for(ii1=0; ii1<MATRIX_DB_HEADER_SIZE; ii1++){
		fprintf(fp, "header[%2d] = %d\n", ii1, mdb->header[ii1]);
	}
#endif

#ifdef WIN32
	for(ii1=0; ii1<MATRIX_DB_SIZE+1; ii1++){
		fprintf(fp, "p_file[%3d] = %d, p_offset[%3d] = %I64d\n",
		             ii1, mdb->p_file[ii1], ii1, mdb->p_offset[ii1]);
	}
#else /* WIN32 */
  #ifdef _LARGEFILE_SOURCE
	for(ii1=0; ii1<MATRIX_DB_SIZE+1; ii1++){
		fprintf(fp, "p_file[%3d] = %d, p_offset[%3d] = %lld\n",
		             ii1, mdb->p_file[ii1], ii1, mdb->p_offset[ii1]);
	}
  #else
	for(ii1=0; ii1<MATRIX_DB_SIZE+1; ii1++){
		fprintf(fp, "p_file[%3d] = %d, p_offset[%3d] = %ld\n",
		             ii1, mdb->p_file[ii1], ii1, mdb->p_offset[ii1]);
	}
  #endif
#endif  /* WIN32 */

	for(ii1=0; ii1<MATRIX_DB_SIZE; ii1++){
		fprintf(fp, "row[%3d] = %d, col[%3d] = %d\n", ii1, mdb->row[ii1],
		                                              ii1, mdb->col[ii1]);
	}

	for(ii1=0; ii1<MATRIX_DB_SIZE; ii1++){
		fprintf(fp, "v[%3d]\n", ii1);
		print_matrix33 (fp, mdb->v[ii1]);
	}

	return(NORMAL_RC);
}


RC
print_matrix_db_debug (FILE *fp, MATRIX_DB *mdb)
{
	int ii1, ii2=0;
	int size;

	RC_NULL_CHK(mdb);

	RC_TRY( print_matrix_db_info(fp, mdb) );
	if(mdb->size == 0) return(NORMAL_RC);
	
#if 1
	if( is_node_matrix_db(mdb) ){
		int mdb_size = mdb->size+1;
		if(mdb_size%2){
			size = mdb_size/2+1;
		}else{
			size = mdb_size/2;
		}
		if(mdb->size == MATRIX_DB_SIZE) size = MATRIX_DB_SIZE/2+1;
		for(ii1=0; ii1<size; ii1++){
			ii2 = ii1+size;
#ifdef WIN32
			fprintf(fp, "file[%3d] = %d, offset[%3d] = %7I64d  "
			            "file[%3d] = %d, offset[%3d] = %7I64d\n",
			             ii1, mdb->p_file[ii1], ii1, mdb->p_offset[ii1],
			             ii2, mdb->p_file[ii2], ii2, mdb->p_offset[ii2]);
#else /* WIN32 */
  #ifdef _LARGEFILE_SOURCE
			fprintf(fp, "file[%3d] = %d, offset[%3d] = %7lld  "
			            "file[%3d] = %d, offset[%3d] = %7lld\n",
			             ii1, mdb->p_file[ii1], ii1, mdb->p_offset[ii1],
			             ii2, mdb->p_file[ii2], ii2, mdb->p_offset[ii2]);
  #else
			fprintf(fp, "file[%3d] = %d, offset[%3d] = %7ld  "
			            "file[%3d] = %d, offset[%3d] = %7ld\n",
			             ii1, mdb->p_file[ii1], ii1, mdb->p_offset[ii1],
			             ii2, mdb->p_file[ii2], ii2, mdb->p_offset[ii2]);
  #endif
#endif  /* WIN32 */
		}
	}
#endif

	if(mdb->size%2){
		size = mdb->size/2+1;
	}else{
		size = mdb->size/2;
	}
	if(mdb->size == MATRIX_DB_SIZE) size = MATRIX_DB_SIZE/2;
	for(ii1=0; ii1<size; ii1++){
		ii2 = ii1+size;
		fprintf(fp, "row[%3d] = %4d, col[%3d] = %4d   "
		            "row[%3d] = %4d, col[%3d] = %4d\n",
		        ii1, mdb->row[ii1], ii1, mdb->col[ii1],
		        ii2, mdb->row[ii2], ii2, mdb->col[ii2]);
	}
	if(mdb->size == MATRIX_DB_SIZE){
		int ii2 = MATRIX_DB_SIZE-1;
		fprintf(fp, "                                   "
		            "row[%3d] = %4d, col[%3d] = %4d\n",
		        ii2, mdb->row[ii2], ii2, mdb->col[ii2]);
	}

#if 0
	/*
	for(ii1=0; ii1<MATRIX_DB_SIZE; ii1++)
	*/
	for(ii1=0; ii1<mdb->size; ii1++){
		fprintf(fp, "v[%2d]\n", ii1);
		print_matrix33 (fp, mdb->v[ii1]);
	}
#endif

	return(NORMAL_RC);
}


RC
print_matrix_db_info (FILE *fp, MATRIX_DB *mdb)
{
	RC_NULL_CHK(mdb);

	fprintf(fp, "size = %d\n", mdb->size);
	fprintf(fp, "dirty_flag = %d\n", mdb->dirty_flag);
	fprintf(fp, "use_count = %d\n", mdb->use_count);
	fprintf(fp, "key_row = %d\n", mdb->key_row);
	fprintf(fp, "key_col = %d\n", mdb->key_col);

#ifdef WIN32
	fprintf(fp, "p_file_parent = %d\n", mdb->p_file_parent);
	fprintf(fp, "p_offset_parent = %I64d\n", mdb->p_offset_parent);
	fprintf(fp, "p_file_self = %d\n", mdb->p_file_self);
	fprintf(fp, "p_offset_self = %I64d\n", mdb->p_offset_self);
#else /* WIN32 */
  #ifdef _LARGEFILE_SOURCE
	fprintf(fp, "p_file_parent = %d\n", mdb->p_file_parent);
	fprintf(fp, "p_offset_parent = %lld\n", mdb->p_offset_parent);
	fprintf(fp, "p_file_self = %d\n", mdb->p_file_self);
	fprintf(fp, "p_offset_self = %lld\n", mdb->p_offset_self);
  #else
	fprintf(fp, "p_file_parent = %d\n", mdb->p_file_parent);
	fprintf(fp, "p_offset_parent = %ld\n", mdb->p_offset_parent);
	fprintf(fp, "p_file_self = %d\n", mdb->p_file_self);
	fprintf(fp, "p_offset_self = %ld\n", mdb->p_offset_self);
  #endif
#endif  /* WIN32 */

	return(NORMAL_RC);
}


void
print_matrix_db_size (FILE *fp)
{
	fprintf(fp, "\n<MATRIX DB SIZE CHECK>\n");
	fprintf(fp, " #define MATRIX_DB_SIZE        = %d\n",MATRIX_DB_SIZE);
	fprintf(fp, " #define MATRIX_DB_HEADER_SIZE = %d\n",MATRIX_DB_HEADER_SIZE);
	fprintf(fp, " typedef struct{\n");

#ifdef WIN32
	fprintf(fp,
	"   [%5d]  int size;\n"
	"   [%5d]  int dirty_flag;\n"
	"   [%5d]  int use_count;\n"
	"   [%5d]  int key_row;\n"
	"   [%5d]  int key_col;\n"
	"   [%5d]  int header[MATRIX_DB_HEADER_SIZE];\n"
	"   [%5d]  int p_file_parent;\n"
	"   [%5d]  int p_file_self;\n"
	"   [%5d]  int p_file[MATRIX_DB_SIZE+1];\n"
	"   [%5d]  WRAP64_INT p_offset_parent;\n"
	"   [%5d]  WRAP64_INT p_offset_self;\n"
	"   [%5d]  WRAP64_INT p_offset[MATRIX_DB_SIZE+1];\n"
	"   [%5d]  int row[MATRIX_DB_SIZE];\n"
	"   [%5d]  int col[MATRIX_DB_SIZE];\n"
	"   [%5d]  double v[MATRIX_DB_SIZE][3][3];\n"
	" } MATRIX_DB;\n"
	" structure size of MATRIX_DB  = %d\n",
	sizeof(int), sizeof(int), sizeof(int), sizeof(int), sizeof(int),
	(MATRIX_DB_HEADER_SIZE)*sizeof(int),
	sizeof(int), sizeof(int), (MATRIX_DB_SIZE+1)*sizeof(int),
	sizeof(WRAP64_INT), sizeof(WRAP64_INT),
	(MATRIX_DB_SIZE+1)*sizeof(WRAP64_INT),
	(MATRIX_DB_SIZE)*sizeof(int), (MATRIX_DB_SIZE)*sizeof(int),
	MATRIX_DB_SIZE*3*3*sizeof(double), sizeof(MATRIX_DB) );
#else /* WIN32 */
  #ifdef _LARGEFILE_SOURCE
	fprintf(fp,
	"   [%5d]  int size;\n" 
	"   [%5d]  int dirty_flag;\n"
	"   [%5d]  int use_count;\n"
	"   [%5d]  int key_row;\n"
	"   [%5d]  int key_col;\n"
	"   [%5d]  int header[MATRIX_DB_HEADER_SIZE];\n"
	"   [%5d]  int p_file_parent;\n"
	"   [%5d]  int p_file_self;\n"
	"   [%5d]  int p_file[MATRIX_DB_SIZE+1];\n"
	"   [%5d]  WRAP64_INT p_offset_parent;\n"
	"   [%5d]  WRAP64_INT p_offset_self;\n"
	"   [%5d]  WRAP64_INT p_offset[MATRIX_DB_SIZE+1];\n"
	"   [%5d]  int row[MATRIX_DB_SIZE];\n"
	"   [%5d]  int col[MATRIX_DB_SIZE];\n"
	"   [%5d]  double v[MATRIX_DB_SIZE][3][3];\n"
	" } MATRIX_DB;\n"
	" structure size of MATRIX_DB  = %d\n",
	sizeof(int), sizeof(int), sizeof(int), sizeof(int), sizeof(int),
	(MATRIX_DB_HEADER_SIZE)*sizeof(int),
	sizeof(int), sizeof(int), (MATRIX_DB_SIZE+1)*sizeof(int),
	sizeof(WRAP64_INT), sizeof(WRAP64_INT),
	(MATRIX_DB_SIZE+1)*sizeof(WRAP64_INT),
	(MATRIX_DB_SIZE)*sizeof(int), (MATRIX_DB_SIZE)*sizeof(int),
	MATRIX_DB_SIZE*3*3*sizeof(double), sizeof(MATRIX_DB) );
  #else /* _LARGEFILE_SOURCE */
	fprintf(fp,
	"   [%5ld]  int size;\n" 
	"   [%5ld]  int dirty_flag;\n"
	"   [%5ld]  int use_count;\n"
	"   [%5ld]  int key_row;\n"
	"   [%5ld]  int key_col;\n"
	"   [%5ld]  int header[MATRIX_DB_HEADER_SIZE];\n"
	"   [%5ld]  int p_file_parent;\n"
	"   [%5ld]  int p_file_self;\n"
	"   [%5ld]  int p_file[MATRIX_DB_SIZE+1];\n"
	"   [%5ld]  WRAP64_INT p_offset_parent;\n"
	"   [%5ld]  WRAP64_INT p_offset_self;\n"
	"   [%5ld]  WRAP64_INT p_offset[MATRIX_DB_SIZE+1];\n"
	"   [%5ld]  int row[MATRIX_DB_SIZE];\n"
	"   [%5ld]  int col[MATRIX_DB_SIZE];\n"
	"   [%5ld]  double v[MATRIX_DB_SIZE][3][3];\n"
	" } MATRIX_DB;\n"
	" structure size of MATRIX_DB  = %ld\n",
	sizeof(int), sizeof(int), sizeof(int), sizeof(int), sizeof(int),
	(MATRIX_DB_HEADER_SIZE)*sizeof(int),
	sizeof(int), sizeof(int), (MATRIX_DB_SIZE+1)*sizeof(int),
	sizeof(WRAP64_INT), sizeof(WRAP64_INT),
	(MATRIX_DB_SIZE+1)*sizeof(WRAP64_INT),
	(MATRIX_DB_SIZE)*sizeof(int), (MATRIX_DB_SIZE)*sizeof(int),
	MATRIX_DB_SIZE*3*3*sizeof(double), sizeof(MATRIX_DB) );
  #endif /* _LARGEFILE_SOURCE */
#endif  /* WIN32 */

	fprintf(fp, "<MATRIX DB SIZE CHECK>\n\n");
}


/* memo
 * このプログラムは
 * 「アルゴリズムイントロダクション 第2巻 アルゴリズムの設計と解析手法」
 * の第19章 「B木」をベースに作られている．
 *
 * 1段 88(x 1) = 88 ブロック
 * 2段 88(x 1+89) = 7,920 ブロック
 * 3段 88(x 1+89 + 89^2) = 704,968 (約 70万)ブロック
 * 4段 88(x 1+89 + 89^2 + 89^3) = 62,742,240 (約 6千万)ブロック
 * 5段 88(x 1+89 + 89^2 + 89^3 + 89^4) = 5,584,059,448 (約 55億)ブロック
 *
 * 2005-09-29 offset を printfで出力する時に %lld (C99) を使わないと，
 * 前後の表示がおかしくなる．
 *
 */


