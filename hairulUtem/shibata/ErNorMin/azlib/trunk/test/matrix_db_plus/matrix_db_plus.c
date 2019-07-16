/*********************************************************************
* mtrix_db_plus.c
 *
 * Copyright (C) 2006 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA>
 *
 * Thanks to following contributors.
 *   <Kenzen TAKEUCHI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id$ */

/*
 * TODO:
 *  - add_row(), delete_row()
 *  - balance_leaf() でマージだけでなく，balance_node()のように分配をする
 *  - splitの可変分割
 */

#include <stdio.h>
#include <math.h>
#ifndef WIN32
#include <unistd.h>
#endif  /* WIN32 */
#include <string.h>
#include <limits.h>
#include "wrapper_64.h"
#include "scratch_io.h"
#include "matrix_db_plus.h"
#include "mathematics.h"
#include "rc.h"
#include "memory_manager.h"
#include "macros.h"

typedef struct{
	int size;
	int alloc_size;
	int index;
	int *row;
	int *col;
	double (*v)[3][3];
} MATRIX_DB_ARRAY;


static RC create_cache_node(MATRIX_DB_CACHE_NODE *node, int cache_size,
                            int associate);
static RC create_cache_leaf(MATRIX_DB_CACHE_LEAF *leaf, int cache_size,
                            int associate);
static RC delete_cache_node(MATRIX_DB_CACHE_NODE *node);
static RC delete_cache_leaf(MATRIX_DB_CACHE_LEAF *leaf);

static RC min_col_search(MATRIX_DB_INFO *mdbi, MATRIX_DB_NODE *node, int row,
                         int col, int *min_col);
static RC min_col_search_leaf(MATRIX_DB_INFO *mdbi, int file,WRAP64_INT offset,
                              int row, int col, int *min_col);

static RC get_row(MATRIX_DB_INFO *mdbi, MATRIX_DB_NODE *node, int row,
                  int col_min, int col_max, int *size, double array[][3][3],
                  int col_array[]);
static RC get_row_leaf(MATRIX_DB_INFO *mdbi, int file, WRAP64_INT offset,
                       int row, int col_min, int col_max, int *size,
                       double array[][3][3], int col_array[]);

static RC add_row(MATRIX_DB_INFO *mdbi, MATRIX_DB_NODE *node,
                  MATRIX_DB_ARRAY *data, MATRIX_DB_ARRAY *buf, int row,
                  int col, int which, int next_row, int next_col, int *retry);
static RC add_row_leaf(MATRIX_DB_INFO *mdbi, int file, WRAP64_INT offset,
                       MATRIX_DB_ARRAY *data, MATRIX_DB_ARRAY *buf, int row,
                       int col, int which, int limit_row, int limit_col,
                       int *retry);
static RC add_row_min_row_col(MATRIX_DB_ARRAY *data, MATRIX_DB_ARRAY *buf,
                              int *row, int *col, int *which);
static RC add_row_clean_buf(MATRIX_DB_ARRAY *buf);

static RC input_tree_leaf(MATRIX_DB_INFO *mdbi, MATRIX_DB_LEAF *leaf,
                          int index, int row, int col, double v[][3]);
static RC split_node(MATRIX_DB_INFO *mdbi, MATRIX_DB_NODE *node,
                     MATRIX_DB_NODE *new_node);
static RC split_leaf(MATRIX_DB_INFO *mdbi, MATRIX_DB_LEAF *leaf,
                     MATRIX_DB_LEAF *new_leaf);
static RC search_parent_edge(MATRIX_DB_INFO *mdbi, MATRIX_DB_LEAF *leaf,
                             int *file, WRAP64_INT *offset);
static RC search_parent_node(MATRIX_DB_INFO *mdbi, MATRIX_DB_NODE *node,
                             int *file, WRAP64_INT *offset);
static RC search_tree_leaf(MATRIX_DB_INFO *mdbi, MATRIX_DB_NODE *node,
                           int row, int col, int *file, WRAP64_INT *offset);
static RC search_tree_node(MATRIX_DB_INFO *mdbi, MATRIX_DB_NODE *node,
                           int row, int col, int *file, WRAP64_INT *offset);
static RC insert_key_node(MATRIX_DB_NODE *node, int index, int row, int col,
                          int file, WRAP64_INT offset);
static RC insert_v_leaf(MATRIX_DB_LEAF *leaf, int index, int row, int col,
                        double v[][3]);
static RC delete_key_node(MATRIX_DB_NODE *node, int file, WRAP64_INT offset);
static RC balance_node(MATRIX_DB_INFO *mdbi, MATRIX_DB_NODE *node);
static RC adjust_node_right(MATRIX_DB_INFO *mdbi, MATRIX_DB_NODE *node,
                            MATRIX_DB_NODE *parent, int index);
static RC adjust_node_left(MATRIX_DB_INFO *mdbi, MATRIX_DB_NODE *node,
                           MATRIX_DB_NODE *parent, int index);
static RC balance_leaf(MATRIX_DB_INFO *mdbi, MATRIX_DB_NODE *parent,
                       MATRIX_DB_LEAF *leaf);

static MATRIX_DB_NODE *search_node_cache_ptr(MATRIX_DB_CACHE_NODE *ca_node,
                                             int file, WRAP64_INT offset,
                                             int inc_access);
static MATRIX_DB_LEAF *search_leaf_cache_ptr(MATRIX_DB_CACHE_LEAF *ca_leaf,
                                             int file, WRAP64_INT offset,
                                             int inc_access);
static MATRIX_DB_NODE *unused_node_cache_ptr(MATRIX_DB_INFO *mdbi, int crow);
static MATRIX_DB_LEAF *unused_leaf_cache_ptr(MATRIX_DB_INFO *mdbi, int crow);

static RC new_node_start(MATRIX_DB_INFO *mdbi, MATRIX_DB_NODE **node);
static RC new_leaf_start(MATRIX_DB_INFO *mdbi, MATRIX_DB_LEAF **leaf);
static RC read_node_start(MATRIX_DB_INFO *mdbi, int file,
                          WRAP64_INT offset, MATRIX_DB_NODE **node);
static RC read_leaf_start(MATRIX_DB_INFO *mdbi, int file,
                          WRAP64_INT offset, MATRIX_DB_LEAF **leaf);
static RC write_node_start(MATRIX_DB_INFO *mdbi, int file,
                           WRAP64_INT offset, MATRIX_DB_NODE **node);
static RC write_leaf_start(MATRIX_DB_INFO *mdbi, int file,
                           WRAP64_INT offset, MATRIX_DB_LEAF **leaf);
static RC new_node_end(MATRIX_DB_NODE *node);
static RC new_leaf_end(MATRIX_DB_LEAF *leaf);
static RC read_node_end(MATRIX_DB_NODE *node);
static RC read_leaf_end(MATRIX_DB_LEAF *leaf);
static RC write_node_end(MATRIX_DB_NODE *node);
static RC write_leaf_end(MATRIX_DB_LEAF *leaf);

static RC delete_node(MATRIX_DB_INFO *mdbi, int file, WRAP64_INT offset);
static RC delete_leaf(MATRIX_DB_INFO *mdbi, int file, WRAP64_INT offset);
static RC put_node(MATRIX_DB_SCRATCH *mdbsf, int file, WRAP64_INT offset,
                   MATRIX_DB_NODE *ptr);
static RC put_leaf(MATRIX_DB_SCRATCH *mdbsf, int file, WRAP64_INT offset,
                   MATRIX_DB_LEAF *ptr);
static RC get_node(MATRIX_DB_SCRATCH *mdbsf, int file, WRAP64_INT offset,
                   MATRIX_DB_NODE *ptr);
static RC get_leaf(MATRIX_DB_SCRATCH *mdbsf, int file, WRAP64_INT offset,
                   MATRIX_DB_LEAF *ptr);

static RC set_root_node(MATRIX_DB_INFO *mdbi, MATRIX_DB_NODE *root,
                        MATRIX_DB_NODE_TYPE type);
static RC init_node(MATRIX_DB_NODE *node);
static RC init_leaf(MATRIX_DB_LEAF *leaf);
static RC clean_node(MATRIX_DB_NODE *node);
static RC clean_leaf(MATRIX_DB_LEAF *leaf);
static RC search_node(MATRIX_DB_NODE *node, int row, int col, int *index,
                      int *found_flag, int *file, WRAP64_INT *offset);
static RC search_leaf(MATRIX_DB_LEAF *leaf, int row, int col,
                      int *index, int *found_flag);

static int is_root_node(MATRIX_DB_NODE *node);
static int is_parent_edge(MATRIX_DB_LEAF *leaf, MATRIX_DB_NODE *edge_node);
static int is_edge_node(MATRIX_DB_NODE *node);
static int is_parent_node(MATRIX_DB_NODE *node, MATRIX_DB_NODE *parent);
/*
static int is_range_node(MATRIX_DB_NODE *node, int index, int row, int col_min,
                         int col_max);
 */
static int is_range_leaf(MATRIX_DB_LEAF *leaf, int index, int row, int col_min,
                         int col_max);
static int compar_row_col(int row1, int col1, int row2, int col2);


#define MAX_ACCESS_COUNT   (INT_MAX/16 + 1)
#define CACHE_ROW(offset, divisor, st_size) ( ((offset)/st_size)%(divisor) )

/* MATRIX_DB 構造体を初期化する時に使用 */
static MATRIX_DB_NODE Initialized_MDB_NODE;
static MATRIX_DB_LEAF Initialized_MDB_LEAF;

#define INIT_MATRIX_DB_NODE(node)  {node = Initialized_MDB_NODE;}
#define INIT_MATRIX_DB_LEAF(leaf)  {leaf = Initialized_MDB_LEAF;}

/*
 * index_cache_size : B+-Treeのツリー(index, node)部のキャッシュサイズ [Byte]
 * data_cache_size  : B+-Treeのデータ(data,  leaf)部のキャッシュサイズ [Byte]
 * *_associate  : 連結されたcache数
 * cache は (sizeof(MATRIX_DB_*)*associate) Byte は必要
 */
RC
open_matrix_db (MATRIX_DB_INFO *mdbi, const char scratch_dir[],
                WRAP64_INT index_cache_size, WRAP64_INT data_cache_size,
                int index_associate, int data_associate)
{
	int ii1;

	RC_NULL_CHK( mdbi );

	if(data_cache_size  < (sizeof(MATRIX_DB_LEAF)*data_associate))
		return(ARG_ERROR_RC);

	RC_TRY( init_node(&Initialized_MDB_NODE) );
	RC_TRY( init_leaf(&Initialized_MDB_LEAF) );

	for(ii1=0; ii1<MATRIX_DB_SF_SIZE; ii1++){
		mdbi->sf_node.scratch_array[ii1] = open_scratch_file(scratch_dir);
		if(mdbi->sf_node.scratch_array[ii1] == NULL) return(OPEN_ERROR_RC);
		mdbi->sf_node.last_offset[ii1] = 0;
		mdbi->sf_node.blank_offset[ii1] = -1;
		mdbi->sf_node.current_offset[ii1] = 0;

		mdbi->sf_leaf.scratch_array[ii1] = open_scratch_file(scratch_dir);
		if(mdbi->sf_leaf.scratch_array[ii1] == NULL) return(OPEN_ERROR_RC);
		mdbi->sf_leaf.last_offset[ii1] = 0;
		mdbi->sf_leaf.blank_offset[ii1] = -1;
		mdbi->sf_leaf.current_offset[ii1] = 0;
	}

	/* キャッシュ */
	RC_TRY( create_cache_node(&(mdbi->cache_node), index_cache_size,
	                          index_associate) );
	RC_TRY( create_cache_leaf(&(mdbi->cache_leaf), data_cache_size,
	                          data_associate) );

	/* 経路情報 */
	mdbi->route_size = 0;
	mdbi->route_count = 0;
	RC_TRY( allocate1I(MATRIX_DB_ROUTE_SIZE, &(mdbi->route_file)) );
	mdbi->route_offset = (WRAP64_INT *)mm_alloc(
	                      (size_t)(MATRIX_DB_ROUTE_SIZE*sizeof(WRAP64_INT)) );
	RC_NULL_CHK( mdbi->route_offset );

	/* Root節点を確定 */
	/* RootNodeを消さないためにreadする */
	RC_TRY( new_node_start(mdbi, &mdbi->root) );
	RC_TRY( set_root_node(mdbi, mdbi->root, MATRIX_DB_ROOT_EDGE) );
	{
		/* 初期leafの作成 */
		MATRIX_DB_LEAF *first_leaf;
		RC_TRY( new_leaf_start(mdbi, &first_leaf) );
		mdbi->root->file[0]   = first_leaf->file_self;
		mdbi->root->offset[0] = first_leaf->offset_self;
		mdbi->root->size = 0;
		RC_TRY( new_leaf_end(first_leaf) );
	}
	RC_TRY( read_node_start(mdbi, mdbi->root->file_self,
	                        mdbi->root->offset_self, &mdbi->root) );
	RC_TRY( new_node_end(mdbi->root) );

	return(NORMAL_RC);
}


RC
close_matrix_db (MATRIX_DB_INFO *mdbi)
{
	int ii1;

	RC_NULL_CHK( mdbi );
	RC_TRY( read_node_end(mdbi->root) );

	RC_TRY( delete_cache_node(&(mdbi->cache_node)) );
	RC_TRY( delete_cache_leaf(&(mdbi->cache_leaf)) );

	RC_TRY( mm_free(mdbi->route_file) );
	RC_TRY( mm_free(mdbi->route_offset) );

	for(ii1=0; ii1<MATRIX_DB_SF_SIZE; ii1++){
		RC_TRY( close_scratch_file(mdbi->sf_node.scratch_array[ii1]) );
		RC_TRY( close_scratch_file(mdbi->sf_leaf.scratch_array[ii1]) );
	}

	return(NORMAL_RC);
}


static RC
create_cache_node (MATRIX_DB_CACHE_NODE *node, int cache_size, int associate)
{
	int ii1, ii2;
	MATRIX_DB_NODE *ptr;

	RC_NULL_CHK(node);
	if(associate <= 0) return(ARG_ERROR_RC);
	if(cache_size < (sizeof(MATRIX_DB_NODE)*associate)) return(ARG_ERROR_RC);

	node->divisor = (int)(cache_size / (sizeof(MATRIX_DB_NODE) * associate));
	RC_NEG_ZERO_CHK(node->divisor);
	node->associate = associate;
	node->cache = (MATRIX_DB_NODE **)mm_alloc(
	                    (size_t)(sizeof(MATRIX_DB_NODE *)*node->divisor) );
	RC_NULL_CHK( node->cache );
	ptr = (MATRIX_DB_NODE *)mm_alloc((size_t)(sizeof(MATRIX_DB_NODE)
	                                 * node->associate * node->divisor) );
	RC_NULL_CHK( ptr );
	for(ii1=0; ii1<node->divisor; ii1++){
		node->cache[ii1] = &ptr[ii1*node->associate];
	}

	RC_TRY( allocate2I(node->divisor, node->associate, &(node->access_count)));
	RC_NULL_CHK( node->access_count );

	for(ii1=0; ii1<node->divisor; ii1++){
		for(ii2=0; ii2<node->associate; ii2++){
			INIT_MATRIX_DB_NODE(node->cache[ii1][ii2]);
			node->access_count[ii1][ii2] = -1;
		}
	}

	return(NORMAL_RC);
}


static RC
create_cache_leaf (MATRIX_DB_CACHE_LEAF *leaf, int cache_size, int associate)
{
	int ii1, ii2;
	MATRIX_DB_LEAF *ptr;

	RC_NULL_CHK(leaf);
	if(associate <= 0) return(ARG_ERROR_RC);
	if(cache_size < (sizeof(MATRIX_DB_LEAF)*associate)) return(ARG_ERROR_RC);

	leaf->divisor = (int)(cache_size / (sizeof(MATRIX_DB_LEAF) * associate));
	RC_NEG_ZERO_CHK(leaf->divisor);
	leaf->associate = associate;
	leaf->cache = (MATRIX_DB_LEAF **)mm_alloc(
	                    (size_t)(sizeof(MATRIX_DB_LEAF *) * leaf->divisor) );
	RC_NULL_CHK( leaf->cache );
	ptr = (MATRIX_DB_LEAF *)mm_alloc((size_t)(sizeof(MATRIX_DB_LEAF)
	                                 * leaf->associate * leaf->divisor) );
	RC_NULL_CHK( ptr );
	for(ii1=0; ii1<leaf->divisor; ii1++){
		leaf->cache[ii1] = &ptr[ii1*leaf->associate];
	}

	RC_TRY( allocate2I(leaf->divisor, leaf->associate, &(leaf->access_count)));
	RC_NULL_CHK( leaf->access_count );

	for(ii1=0; ii1<leaf->divisor; ii1++){
		for(ii2=0; ii2<leaf->associate; ii2++){
			INIT_MATRIX_DB_LEAF(leaf->cache[ii1][ii2]);
			leaf->access_count[ii1][ii2] = -1;
		}
	}

	return(NORMAL_RC);
}


static RC
delete_cache_node (MATRIX_DB_CACHE_NODE *node)
{
	RC_NULL_CHK(node);

	RC_TRY( mm_free(node->cache[0]) );
	RC_TRY( mm_free(node->cache) );
	RC_TRY( free2I(node->divisor, node->associate, &(node->access_count)));
	node->associate = 0;
	node->divisor = 0;

	return(NORMAL_RC);
}


static RC
delete_cache_leaf (MATRIX_DB_CACHE_LEAF *leaf)
{
	RC_NULL_CHK(leaf);

	RC_TRY( mm_free(leaf->cache[0]) );
	RC_TRY( mm_free(leaf->cache) );
	RC_TRY( free2I(leaf->divisor, leaf->associate, &(leaf->access_count)));
	leaf->associate = 0;
	leaf->divisor = 0;

	return(NORMAL_RC);
}



/* row 行目の col以上(col含む)の最小値を探す         */
/* min_col == -1 ならrow行目がない，もしくは検索失敗 */
RC
min_col_matrix_db (MATRIX_DB_INFO *mdbi, int row, int col, int *min_col)
{
	RC_NULL_CHK( mdbi );
	RC_NULL_CHK( min_col );
	if((row < 0) || (col < 0)) return(ARG_ERROR_RC);

	*min_col = -1;
	RC_TRY( min_col_search(mdbi, mdbi->root, row, col, min_col) );

	return(NORMAL_RC);
}


static RC
min_col_search (MATRIX_DB_INFO *mdbi, MATRIX_DB_NODE *node, int row, int col,
                int *min_col)
{
	int index;
	int found_flag;
	int file;
	WRAP64_INT offset;

	RC_NULL_CHK( mdbi );
	RC_NULL_CHK( node );
	RC_NULL_CHK( min_col );
	if( (row < 0) || (col < 0) ) return(ARG_ERROR_RC);

	RC_TRY( search_node(node, row, col, &index, &found_flag, &file, &offset) );
	if( (file < 1) || (offset < 0) || (index < 0) ){
		RC_TRY( print_matrix_db_node(stderr, node) );
		return(UNKNOWN_ERROR_RC);
	}

	if( is_edge_node(node) ){
		if(index < node->size){
			int comp = compar_row_col(node->row[index], node->col[index],
			                          row, col);
			if(comp > 0){
				/* index の左右のleafを探す */
				RC_TRY( min_col_search_leaf(mdbi, node->file[index],
				                            node->offset[index], row, col,
				                            min_col) );
				if(*min_col >= 0) return(NORMAL_RC);
				if(node->row[index] > row) return(NORMAL_RC);
				RC_TRY( min_col_search_leaf(mdbi, node->file[index+1],
				                            node->offset[index+1], row, col,
				                            min_col) );
			}else{
				/* indexの右のleafを探す */
				RC_TRY( min_col_search_leaf(mdbi, node->file[index+1],
				                            node->offset[index+1], row, col,
				                            min_col) );
			}
		}else if(index == node->size){
			if(index > 0){
				if(node->row[index-1] > row) return(NORMAL_RC);
			}
			RC_TRY( min_col_search_leaf(mdbi, node->file[index],
			                            node->offset[index], row, col,
			                            min_col) );
		}else{
			/* do nothing */
		}
	}else{
		MATRIX_DB_NODE *child;
		if(index < node->size){
			int comp = compar_row_col(node->row[index], node->col[index],
			                          row, col);
			if(comp > 0){
				RC_TRY( read_node_start(mdbi, node->file[index],
				                        node->offset[index], &child) );
				RC_TRY( min_col_search(mdbi, child, row, col, min_col) );
				RC_TRY( read_node_end(child) );
				if(*min_col >= 0) return(NORMAL_RC);
				if(node->row[index] > row) return(NORMAL_RC);
				RC_TRY( read_node_start(mdbi, node->file[index+1],
				                        node->offset[index+1], &child) );
				RC_TRY( min_col_search(mdbi, child, row, col, min_col) );
				RC_TRY( read_node_end(child) );
			}else{
				RC_TRY( read_node_start(mdbi, node->file[index+1],
				                        node->offset[index+1], &child) );
				RC_TRY( min_col_search(mdbi, child, row, col, min_col) );
				RC_TRY( read_node_end(child) );
			}
		}else if(index == node->size){
			if(index > 0){
				if(node->row[index-1] > row) return(NORMAL_RC);
			}
			RC_TRY( read_node_start(mdbi, node->file[index],
			                              node->offset[index], &child) );
			RC_TRY( min_col_search(mdbi, child, row, col, min_col) );
			RC_TRY( read_node_end(child) );
		}else{
			/* do nothing */
		}
	}

	return(NORMAL_RC);
}


/* min_col_search sub 関数 */
static RC
min_col_search_leaf (MATRIX_DB_INFO *mdbi, int file, WRAP64_INT offset,
                     int row, int col, int *min_col)
{
	MATRIX_DB_LEAF *leaf;
	int index;
	int found_flag;

	RC_NULL_CHK( mdbi );
	RC_NULL_CHK( min_col );
	if( (row < 0) || (col < 0) ) return(ARG_ERROR_RC);
	if( (file < 1) || (offset < 0) ) return(ARG_ERROR_RC);

	RC_TRY( read_leaf_start(mdbi, file, offset, &leaf) );
	RC_TRY( search_leaf(leaf, row, col, &index, &found_flag) );
	if( (index < 0) || (index >= MATRIX_DB_LEAF_SIZE) )
		return(UNKNOWN_ERROR_RC);
	if(index < leaf->size){
		if(leaf->row[index] == row){
			*min_col = leaf->col[index];
		}
	}
	RC_TRY( read_leaf_end(leaf) );

	return(NORMAL_RC);
}


/*
 * row 行目の row 行目の col_min列(含む)からcol_max列(含まない)までの
 * 登録されたデータ全部を array及び col_array にコピーする．コピーした
 * 数が size に代入される．array, col_arrayは呼び出し側で十分な大きさが
 * 確保されているものとする．
 */
RC
get_row_matrix_db (MATRIX_DB_INFO *mdbi, int row, int col_min, int col_max,
                   int *size, double array[][3][3], int col_array[])
{
	*size = 0;
	RC_TRY( get_row(mdbi, mdbi->root, row, col_min, col_max,
	                size, array, col_array) );
	return(NORMAL_RC);
}


static RC
get_row (MATRIX_DB_INFO *mdbi, MATRIX_DB_NODE *node, int row, int col_min,
         int col_max, int *size, double array[][3][3], int col_array[])
{
	int ii1;
	int index;
	int found_flag; 
	int file;
	WRAP64_INT offset;

	RC_NULL_CHK( mdbi );
	RC_NULL_CHK( node );
	RC_NULL_CHK( size );
	RC_NULL_CHK( array );
	RC_NULL_CHK( col_array );
	if( (row < 0)||(col_min < 0) ) return(ARG_ERROR_RC);
	if(col_min >= col_max) return(ARG_ERROR_RC);

	RC_TRY( search_node(node, row, col_min, &index, &found_flag,
	                    &file, &offset) );
	if( (file < 1) || (offset < 0) || (index < 0) ){
		RC_TRY( print_matrix_db_node(stderr, node) );
		return(UNKNOWN_ERROR_RC);
	}

	if( is_edge_node(node) ){
		if(index < node->size){
			int comp = compar_row_col(node->row[index], node->col[index],
			                          row, col_min);
			if(comp > 0){
				RC_TRY( get_row_leaf(mdbi, node->file[index],
				                     node->offset[index], row, col_min,
				                     col_max, size, array, col_array) );
			}
			for(ii1=index; ii1<node->size; ii1++){
				comp = compar_row_col(node->row[ii1], node->col[ii1],
				                      row, col_max); /* not col_min !! */
				if(comp >= 0) break;
				RC_TRY( get_row_leaf(mdbi, node->file[ii1+1],
				                     node->offset[ii1+1], row, col_min,
				                     col_max, size, array, col_array) );
			}
		}else if(index == node->size){
			if(index > 0){
				if(node->row[index-1] > row) return(NORMAL_RC);
			}
			RC_TRY( get_row_leaf(mdbi, node->file[index], node->offset[index],
			                     row, col_min, col_max, size, array,
			                     col_array) );
		}else{
			/* do nothing */
		}
	}else{
		MATRIX_DB_NODE *child;
		if(index < node->size){
			int comp = compar_row_col(node->row[index], node->col[index],
			                          row, col_min);
			if(comp > 0){
				RC_TRY( read_node_start(mdbi, node->file[index],
				                              node->offset[index], &child) );
				RC_TRY( get_row(mdbi, child, row, col_min, col_max, size,
				                array, col_array) );
				RC_TRY( read_node_end(child) );
			}
			for(ii1=index; ii1<node->size; ii1++){
				comp = compar_row_col(node->row[ii1], node->col[ii1],
				                      row, col_max); /* not col_min !! */
				if(comp >= 0) break;
				RC_TRY( read_node_start(mdbi, node->file[ii1+1],
				                              node->offset[ii1+1], &child) );
				RC_TRY( get_row(mdbi, child, row, col_min, col_max, size,
				                array, col_array) );
				RC_TRY( read_node_end(child) );
			}
		}else if(index == node->size){
			if(index > 0){
				if(node->row[index-1] > row) return(NORMAL_RC);
			}
			RC_TRY( read_node_start(mdbi, node->file[index],
			                              node->offset[index], &child) );
			RC_TRY( get_row(mdbi, child, row, col_min, col_max, size,
			                array, col_array) );
			RC_TRY( read_node_end(child) );
		}else{
			/* do nothing */
		}
	}

	return(NORMAL_RC);
}


static RC
get_row_leaf (MATRIX_DB_INFO *mdbi, int file, WRAP64_INT offset, int row,
              int col_min, int col_max, int *size, double array[][3][3],
              int col_array[])
{
	int ii1, index, found_flag;
	MATRIX_DB_LEAF *leaf;

	RC_NULL_CHK( size );
	RC_NULL_CHK( array );
	RC_NULL_CHK( col_array );
	if( (file < 1) || (offset < 0) ) return(ARG_ERROR_RC);
	if( (row < 0)||(col_min < 0) ) return(ARG_ERROR_RC);
	if(col_min >= col_max) return(ARG_ERROR_RC);

	RC_TRY( read_leaf_start(mdbi, file, offset, &leaf) );
	RC_TRY( search_leaf(leaf, row, col_min, &index, &found_flag) );
	if(index < 0) return(UNKNOWN_ERROR_RC);
	for(ii1=index; ii1<leaf->size; ii1++){
		if( is_range_leaf(leaf, ii1, row, col_min, col_max) ){
			RC_TRY( copy_matrix33(leaf->v[ii1], array[*size]) );
			col_array[*size] = leaf->col[ii1];
			(*size)++;
		}else{
			break;
		}
	}
	RC_TRY( read_leaf_end(leaf) );

	return(NORMAL_RC);
}


/*
 * row行目のデータを連続で addする．
 * col_array & array はソートされている必要がある．
 *
 * 追加するデータ，Treeにあるデータ，バッファにあるデータの三つを比べ
 * 必要があれば Treeからバッファにデータを退避し，小さいものから順に
 * add, insertをしていく．
 */
RC
add_row_matrix_db (MATRIX_DB_INFO *mdbi, int row, int size,
                   double array[][3][3], int col_array[])
{
	MATRIX_DB_ARRAY buf;
	MATRIX_DB_ARRAY data;
	int col, which;
	int retry;

	RC_NULL_CHK(mdbi);
	RC_NULL_CHK(array);
	RC_NULL_CHK(col_array);
	if( (row < 0)||(size < 0) ) return(ARG_ERROR_RC);

	if(size == 0) return(NORMAL_RC);

	/* バッファ */
	buf.size = 0;
	buf.alloc_size = size + MATRIX_DB_LEAF_SIZE; /* sizeでもいいかもしれない */
	buf.index = 0;
	RC_TRY( allocate1I(buf.alloc_size, &(buf.row)) );
	RC_TRY( allocate1I(buf.alloc_size, &(buf.col)) );
	RC_TRY( allocate1D33(buf.alloc_size, &(buf.v)) );

	/* 追加するデータ */
	data.size = size;
	data.alloc_size = size;
	data.index = 0;
	RC_TRY( allocate1I(1, &(data.row)) ); /* 全部同じなので 1つでいい */
	data.row[0] = row;
	data.col = col_array;
	data.v = array;

	/* 入力する dataとbufの row,colの最小値を探す */
	RC_TRY( add_row_min_row_col(&data, &buf, &row, &col, &which) );

	while(1){
		retry = 0;
		RC_TRY( add_row(mdbi, mdbi->root, &data, &buf, row, col, which,
		                INT_MAX, INT_MAX, &retry) );
		RC_TRY( add_row_min_row_col(&data, &buf, &row, &col, &which) );
		if(row < 0) break;
		log_printf(5, "[%d] retry\n", __LINE__);
	}

	RC_TRY( free1I(buf.alloc_size, &(buf.row)) );
	RC_TRY( free1I(buf.alloc_size, &(buf.col)) );
	RC_TRY( free1D33(buf.alloc_size, &(buf.v)) );

	RC_TRY( free1I(1, &(data.row)) );

	return(NORMAL_RC);
}


static RC
add_row (MATRIX_DB_INFO *mdbi, MATRIX_DB_NODE *node, MATRIX_DB_ARRAY *data,
         MATRIX_DB_ARRAY *buf, int row, int col, int which, int limit_row,
         int limit_col, int *retry)
{
	int ii1;
	int index;
	int found_flag;
	int file;
	WRAP64_INT offset;

	RC_NULL_CHK( mdbi );
	RC_NULL_CHK( node );
	RC_NULL_CHK( buf );
	RC_NULL_CHK( data );
	if( (row < 0) || (col < 0) ) return(ARG_ERROR_RC);
	if( (limit_row < 0) || (limit_col < 0) ) return(ARG_ERROR_RC);
	if( (which < -1) || (which > 1) ) return(ARG_ERROR_RC);

	RC_TRY( search_node(node, row, col, &index, &found_flag, &file, &offset) );
	if( (file < 1) || (offset < 0) || (index < 0) ){
		RC_TRY( print_matrix_db_node(stderr, node) );
		return(UNKNOWN_ERROR_RC);
	}

	if( is_edge_node(node) ){
		log_printf(5, "[%d] \x1b[33mSplit Node\x1b[0m (Edge) \n", __LINE__);
		if(index < node->size){
			int comp = compar_row_col(node->row[index],
			                          node->col[index], row, col);
			if(comp > 0){
				RC_TRY( add_row_leaf(mdbi, node->file[index],
				                     node->offset[index], data, buf, row, col,
				                     which, node->row[index],
				                     node->col[index], retry) );
				if(retry) return(NORMAL_RC);
				RC_TRY( add_row_min_row_col(data, buf, &row, &col, &which) );
				if(row < 0) return(NORMAL_RC);
			}
			for(ii1=index; ii1<node->size; ii1++){
				comp = compar_row_col(node->row[ii1], node->col[ii1], row,col);
				if(comp > 0) return(UNKNOWN_ERROR_RC);
				if(ii1+1 == node->size){
					log_printf(5, "limit_row = %d\n", limit_row);
					log_printf(5, "limit_col = %d\n", limit_col);
					log_printf(5, "row = %d\n", row);
					log_printf(5, "col = %d\n", col);
					RC_TRY( print_matrix_db_node(stderr, node) );
					RC_TRY( add_row_leaf(mdbi, node->file[ii1+1],
					                     node->offset[ii1+1], data, buf, row,
					                     col, which, limit_row, limit_col,
					                     retry) );
				}else{
					comp = compar_row_col(node->row[ii1+1],
					                      node->col[ii1+1], row, col);
					if(comp > 0){
						RC_TRY( add_row_leaf(mdbi, node->file[ii1+1],
						                     node->offset[ii1+1], data, buf,
						                     row, col, which, node->row[ii1+1],
						                     node->col[ii1+1], retry) );
					}
				}
				if(retry) return(NORMAL_RC);
				RC_TRY( add_row_min_row_col(data, buf, &row, &col, &which) );
				if(row < 0) return(NORMAL_RC);
			}
		}else if(index == node->size){
			if(index > 0){
				if(node->row[index-1] > row) return(NORMAL_RC);
			}
			RC_TRY( add_row_leaf(mdbi, node->file[index], node->offset[index],
			                     data, buf, row, col, which, limit_row,
			                     limit_col, retry) );
			if(retry) return(NORMAL_RC);
		}else{
			/* do nothing */
		}
	}else{
		log_printf(5, "[%d] \x1b[33mSplit Node\x1b[0m (Node) \n", __LINE__);
		MATRIX_DB_NODE *child;
		if(index < node->size){
			int comp = compar_row_col(node->row[index], node->col[index],
			                          row, col);
			if(comp > 0){
				RC_TRY( read_node_start(mdbi, node->file[index],
				                              node->offset[index], &child) );
				RC_TRY( add_row(mdbi, child, data, buf, row, col, which,
				                limit_row, limit_col, retry) );
				RC_TRY( read_node_end(child) );
				if(retry) return(NORMAL_RC);
				RC_TRY( add_row_min_row_col(data, buf, &row, &col, &which) );
				if(row < 0) return(NORMAL_RC);
			}
			for(ii1=index; ii1<node->size; ii1++){
				comp = compar_row_col(node->row[ii1], node->col[ii1], row,col);
				if(comp > 0){
					log_printf(5, "node->row[%d] = %d\n", ii1, node->row[ii1]);
					log_printf(5, "node->col[%d] = %d\n", ii1, node->col[ii1]);
					log_printf(5, "row = %d\n", row);
					log_printf(5, "col = %d\n", col);
					return(UNKNOWN_ERROR_RC);
				}
				if(ii1+1 == node->size){
					RC_TRY( read_node_start(mdbi, node->file[ii1+1],
					                        node->offset[ii1+1], &child) );
					RC_TRY( add_row(mdbi, child, data, buf, row, col, which,
					                limit_row, limit_col, retry) );
					RC_TRY( read_node_end(child) );
				}else{
					comp = compar_row_col(node->row[ii1+1],
					                      node->col[ii1+1], row, col);
					if(comp > 0){
						RC_TRY( read_node_start(mdbi, node->file[ii1+1],
						                        node->offset[ii1+1], &child) );
						RC_TRY( add_row(mdbi, child, data, buf, row, col,
						                which, node->row[ii1+1],
						                node->col[ii1+1], retry) );
						RC_TRY( read_node_end(child) );
					}
				}
				if(retry) return(NORMAL_RC);
				RC_TRY( add_row_min_row_col(data, buf, &row, &col, &which) );
				if(row < 0) return(NORMAL_RC);
			}
		}else if(index == node->size){
			if(index > 0){
				if(node->row[index-1] > row) return(NORMAL_RC);
			}
			RC_TRY( read_node_start(mdbi, node->file[index],
			                              node->offset[index], &child) );
			RC_TRY( add_row(mdbi, child, data, buf, row, col, which,
			                limit_row, limit_col, retry) );
			RC_TRY( read_node_end(child) );
			if(retry) return(NORMAL_RC);
		}else{
			/* do nothing */
		}
	}

	return(NORMAL_RC);
}


static RC
add_row_leaf (MATRIX_DB_INFO *mdbi, int file, WRAP64_INT offset,
              MATRIX_DB_ARRAY *data, MATRIX_DB_ARRAY *buf, int row, int col,
              int which, int limit_row, int limit_col, int *retry)
{
	int index;
	int found_flag;
	int comp;
	MATRIX_DB_LEAF *leaf;
	MATRIX_DB_LEAF *new_leaf;

	RC_NULL_CHK( mdbi );
	RC_NULL_CHK( data );
	RC_NULL_CHK( buf );
	if( (file < 1) || (offset < 0) ) return(ARG_ERROR_RC);
	if( (row < 0) || (col < 0) ) return(ARG_ERROR_RC);
	if( (limit_row < 0) || (limit_col < 0) ) return(ARG_ERROR_RC);
	if( (which < -1) || (which > 1) ) return(ARG_ERROR_RC);

	*retry = 0;
	RC_TRY( write_leaf_start(mdbi, file, offset, &leaf) );

	DEBUG_PRINT;
	RC_TRY( print_matrix_db_leaf_info(stderr, leaf) );
	DEBUG_PRINT;

	if(leaf->size >= MATRIX_DB_LEAF_SIZE){
		RC_TRY( new_leaf_start(mdbi, &new_leaf) );
		RC_TRY( split_leaf(mdbi, leaf, new_leaf) );
		comp = compar_row_col(new_leaf->row[0], new_leaf->col[0], row, col);
		if(comp <= 0){ /* key1 <= key2 */
			RC_TRY( write_leaf_end(leaf) );
			leaf = new_leaf;
		}else{ /* (comp > 0) key1 > key2 */
			limit_row = new_leaf->row[0];
			limit_col = new_leaf->col[0];
			RC_TRY( new_leaf_end(new_leaf) );
		}
		*retry = 1;
		DEBUG_PRINT;
		RC_TRY( print_matrix_db_tree(stderr, mdbi) );
		DEBUG_PRINT;
	}

	RC_TRY( search_leaf(leaf, row, col, &index, &found_flag) );
	if(index < 0) return(UNKNOWN_ERROR_RC);

	while(1){
		DEBUG_PRINT;
		RC_TRY( print_matrix_db_leaf(stderr, leaf) );
		DEBUG_PRINT;
		/* 親ノードのキー値より大きい値になったら終わり */
		comp = compar_row_col(row, col, limit_row, limit_col);
		log_printf(5, "row = %d\n", row);
		log_printf(5, "col = %d\n", col);
		log_printf(5, "limit_row = %d\n", limit_row);
		log_printf(5, "limit_col = %d\n", limit_col);
		if(comp >= 0){ /* key1 >= key2 */
			DEBUG_PRINT;
			break;
		}

		if(index == leaf->size){
			comp = -1; /* leafの後半に入力 */
		}else{
			comp = compar_row_col(row, col, leaf->row[index],leaf->col[index]);
		}

		if(comp < 0){ /* key1 < key2 */
			if(index < leaf->size){
				/* leaf->*[index]をバッファへ退避 */
				buf->row[buf->size] = leaf->row[index];
				buf->col[buf->size] = leaf->col[index];
				RC_TRY( copy_matrix33(leaf->v[index], buf->v[buf->size]) );
				buf->size++;
			}
			leaf->row[index] = row;
			leaf->col[index] = col;
			if(which == -1){ /* data < buf */
				RC_TRY( copy_matrix33(data->v[data->index], leaf->v[index]) );
				data->index++;
			}else if(which == 0){ /* data = buf */
				RC_TRY( add_matrix33a(data->v[data->index],
				                      buf->v[buf->index], leaf->v[index]) );
				data->index++;
				buf->index++;
			}else{ /* (which == 1), data > buf */
				RC_TRY( copy_matrix33(buf->v[buf->index], leaf->v[index]) );
				buf->index++;
			}
			if(index == leaf->size) leaf->size++;
		}else if(comp == 0){ /* key1 = key2 */
			if(which == -1){ /* data < buf */
				add_matrix33(leaf->v[index], data->v[data->index]);
				data->index++;
			}else if(which == 0){ /* data = buf */
				add_matrix33(leaf->v[index], data->v[data->index]);
				add_matrix33(leaf->v[index], buf->v[buf->index]);
				data->index++;
				buf->index++;
			}else{ /* (which == 1), data > buf */
				add_matrix33(leaf->v[index], buf->v[buf->index]);
				buf->index++;
			}
		}else{ /* (comp > 0) key1 > key2 */
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
			return(NORMAL_RC); /* 入力するデータがないので終了 */
		}

		if(leaf->size == MATRIX_DB_LEAF_SIZE){
			RC_TRY( new_leaf_start(mdbi, &new_leaf) );
			RC_TRY( split_leaf(mdbi, leaf, new_leaf) );
			RC_TRY( write_leaf_end(leaf) );
			leaf = new_leaf;
			// RC_TRY( print_matrix_db_tree(stderr, mdbi) );
			// exit(0);
			RC_TRY( search_leaf(leaf, row, col, &index, &found_flag) );
			if(index < 0) return(UNKNOWN_ERROR_RC);
			*retry = 1;
			DEBUG_PRINT;
			RC_TRY( print_matrix_db_tree(stderr, mdbi) );
			DEBUG_PRINT;
		}
	}

	RC_TRY( write_leaf_end(leaf) );

	return(NORMAL_RC);
}


/* data, buf の各index番目を比較する                          */
/* data のrow,colが最小なら which = -1, bufなら 1, 共通なら 0 */
static RC
add_row_min_row_col (MATRIX_DB_ARRAY *data, MATRIX_DB_ARRAY *buf,
                     int *row, int *col, int *which)
{

	RC_NULL_CHK( data );
	RC_NULL_CHK( buf );
	RC_NULL_CHK( row );
	RC_NULL_CHK( col );
	RC_NULL_CHK( which );

	*row = -1;
	*col = -1;
	*which = 0;

	if(data->index >= data->size){ /* data には被入力データが無い */
		if(buf->index < buf->size){
			if(buf->index >= buf->size) return(UNKNOWN_ERROR_RC);
			/* bufにはある */
			*row = buf->row[buf->index];
			*col = buf->col[buf->index];
			*which = 1;
		}
	}else if(buf->index >= buf->size){ /* buf には被入力データが無い */
		if(data->index < data->size){
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
		}else{ /* (*which == 1), dataの方が大きい */
			*row = buf->row[buf->index];
			*col = buf->col[buf->index];
		}
	}

	return(NORMAL_RC);
}


/* 左に詰める */
static RC
add_row_clean_buf (MATRIX_DB_ARRAY *buf)
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


RC
add_v_matrix_db (int row, int col, double v[][3], MATRIX_DB_INFO *mdbi)
{
	int index, found_flag;
	MATRIX_DB_LEAF *leaf;
	int file;
	WRAP64_INT offset;

	RC_NULL_CHK( mdbi );
	if( (row < 0) || (col < 0) ) return(ARG_ERROR_RC);

	/* 該当leafを検索 */
	mdbi->route_size = 0;
	RC_TRY( search_tree_leaf(mdbi, mdbi->root, row, col, &file, &offset) );

	if((file <= 0) || (file > MATRIX_DB_SF_SIZE)) return(UNKNOWN_ERROR_RC);
	if(offset < 0)  return(UNKNOWN_ERROR_RC);

	log_printf(6, "[%d] row = %d, col = %d\n", __LINE__, row, col);

	RC_TRY( write_leaf_start(mdbi, file, offset, &leaf) );
	RC_TRY( search_leaf(leaf, row, col, &index, &found_flag) );
	if(found_flag){
		log_printf(5, "[%d] OverWrite row=%d, col=%d\n", __LINE__, row, col);
		if( (index < 0)||(index >= leaf->size) ) return(UNKNOWN_ERROR_RC);
		leaf->row[index] = row;
		leaf->col[index] = col;
		add_matrix33(leaf->v[index], v);
	}else{
		if( (index < 0)||(index > leaf->size) ) return(UNKNOWN_ERROR_RC);
		RC_TRY( input_tree_leaf(mdbi, leaf, index, row, col, v) );
	}
	RC_TRY( write_leaf_end(leaf) );

	/*
	if( (row==974) && (col==498) ){
		RC_TRY( print_matrix_db_tree(stderr, mdbi) );
		RC_TRY( print_matrix_db_cache_node_info(stderr, &mdbi->cache_node) );
		RC_TRY( print_matrix_db_cache_leaf_info(stderr, &mdbi->cache_leaf) );

	}
	*/

	return(NORMAL_RC);
}


RC
get_v_matrix_db (int row, int col, double v[][3], int *found_flag,
                 MATRIX_DB_INFO *mdbi)
{
	int index;
	int file;
	WRAP64_INT offset;
	MATRIX_DB_LEAF *leaf;

	RC_NULL_CHK( mdbi );
	RC_NULL_CHK( found_flag );
	if((row < 0) || (col < 0)) return(ARG_ERROR_RC);

	/* 挿入すべき Leaf Nodeを Rootから探索 */
	mdbi->route_size = 0;
	RC_TRY( search_tree_leaf(mdbi, mdbi->root, row, col, &file, &offset) );
	if((file <= 0) || (file > MATRIX_DB_SF_SIZE)) return(UNKNOWN_ERROR_RC);
	if(offset < 0)  return(UNKNOWN_ERROR_RC);

	RC_TRY( read_leaf_start(mdbi, file, offset, &leaf) );
	RC_TRY( search_leaf(leaf, row, col, &index, found_flag) );
	if(*found_flag){
		if( (index < 0)||(index >= leaf->size) ) return(UNKNOWN_ERROR_RC);
		RC_TRY( copy_matrix33(leaf->v[index], v) );
	}
	RC_TRY( read_leaf_end(leaf) );

	return(NORMAL_RC);
}


static RC
input_tree_leaf (MATRIX_DB_INFO *mdbi, MATRIX_DB_LEAF *leaf,
                 int index, int row, int col, double v[][3])
{
	int found_flag;

	RC_NULL_CHK( leaf );
	RC_NULL_CHK( mdbi );
	if( (row < 0) || (col < 0) ) return(ARG_ERROR_RC);
	if( (index < 0) || (index >= MATRIX_DB_LEAF_SIZE) ) return(ARG_ERROR_RC);

	if(leaf->size >= MATRIX_DB_LEAF_SIZE){
		MATRIX_DB_LEAF *new_leaf;
		int comp;

		RC_TRY( new_leaf_start(mdbi, &new_leaf) );
		RC_TRY( split_leaf(mdbi, leaf, new_leaf) );
		comp = compar_row_col(new_leaf->row[0], new_leaf->col[0], row, col);
		if(comp < 0){ /* key1 < key2 */
			RC_TRY( search_leaf(new_leaf, row, col, &index, &found_flag) );
			if( (index < 0)||(index > new_leaf->size) )
				return(UNKNOWN_ERROR_RC);
			RC_TRY( insert_v_leaf(new_leaf, index, row, col, v) );
		}else if(comp == 0){ /* key1 = key2 */
			/* これはありえない */
			return(UNKNOWN_ERROR_RC);
		}else{ /* (comp > 0) key1 > key2 */
			RC_TRY( search_leaf(leaf, row, col, &index, &found_flag) );
			if( (index < 0)||(index > leaf->size) ) return(UNKNOWN_ERROR_RC);
			RC_TRY( insert_v_leaf(leaf, index, row, col, v) );
		}
		RC_TRY( new_leaf_end(new_leaf) );
	}else{
		RC_TRY( insert_v_leaf(leaf, index, row, col, v) );
	}

	return(NORMAL_RC);
}


/* fullになったnodeを分割 */
static RC
split_node (MATRIX_DB_INFO *mdbi, MATRIX_DB_NODE *node,
            MATRIX_DB_NODE *new_node)
{
	int ii1, ii2;
	int index, found_flag, file;
	WRAP64_INT offset;
	int parent_new_flag = 0;
	MATRIX_DB_NODE *parent;
	int key_row, key_col;

	RC_NULL_CHK( mdbi );
	RC_NULL_CHK( node );
	RC_NULL_CHK( new_node );
	if(node->size < MATRIX_DB_NODE_SIZE/2+1) return(ARG_ERROR_RC);
	if(new_node->size != 0) return(ARG_ERROR_RC);

	/* parent node の選定 */
	if( is_root_node(node) ){
		log_printf(6, "[%d] \x1b[33mSplit Node\x1b[0m (Root) \n", __LINE__);
		RC_TRY( new_node_start(mdbi, &parent) );
		RC_TRY( set_root_node(mdbi, parent, MATRIX_DB_ROOT) );
		if(node->type == MATRIX_DB_ROOT_EDGE){
			node->type     = MATRIX_DB_EDGE;
			new_node->type = MATRIX_DB_EDGE;
		}else if(node->type == MATRIX_DB_ROOT){
			node->type     = MATRIX_DB_GRID;
			new_node->type = MATRIX_DB_GRID;
		}else{
			return(UNKNOWN_ERROR_RC);
		}
		parent->file[0]   = node->file_self;
		parent->offset[0] = node->offset_self;
		parent_new_flag = 1;
	}else{
		log_printf(6, "[%d] \x1b[33mSplit Node\x1b[0m (Node) \n", __LINE__);
		mdbi->route_size = 0; /* 必須 */
		RC_TRY( search_parent_node(mdbi, node, &file, &offset) );
		if((file <= 0) || (file > MATRIX_DB_SF_SIZE)) return(UNKNOWN_ERROR_RC);
		if(offset < 0) return(UNKNOWN_ERROR_RC);
		RC_TRY( write_node_start(mdbi,file, offset, &parent) );
		/* 親が一杯なら親を分割 */
		if(parent->size >= MATRIX_DB_NODE_SIZE){
			MATRIX_DB_NODE *new_parent;
			RC_TRY( new_node_start(mdbi, &new_parent) );
			new_parent->type = parent->type;
			RC_TRY( split_node(mdbi, parent, new_parent) );
			if( is_parent_node(node, parent) ){
				RC_TRY( new_node_end(new_parent) );
			}else if( is_parent_node(node, new_parent) ){
				RC_TRY( write_node_end(parent) );
				parent = new_parent;
				parent_new_flag = 1;
			}else{
				return(UNKNOWN_ERROR_RC);
			}
		}
	}

	/* nodeの後半を new_nodeへ移動 */
	for(ii1=(MATRIX_DB_NODE_SIZE/2+1), ii2=0; ii1<node->size; ii1++){
		new_node->row[ii2] = node->row[ii1];
		new_node->col[ii2] = node->col[ii1];
		new_node->file[ii2]   = node->file[ii1];
		new_node->offset[ii2] = node->offset[ii1];
		ii2++;
	}
	new_node->file[ii2]   = node->file[ii1];
	new_node->offset[ii2] = node->offset[ii1];
	new_node->size = ii2;

	key_row = node->row[MATRIX_DB_NODE_SIZE/2];
	key_col = node->col[MATRIX_DB_NODE_SIZE/2];
	node->size = MATRIX_DB_NODE_SIZE/2;
	RC_TRY( clean_node(node) );

	/* parentへ挿入  */
	RC_TRY( search_node(parent, key_row, key_col, &index, &found_flag,
	                    &file, &offset) );
	if(index < 0) return(UNKNOWN_ERROR_RC);
	RC_TRY( insert_key_node(parent, index, key_row, key_col,
	                        new_node->file_self, new_node->offset_self) );

	if(parent_new_flag){
		RC_TRY( new_node_end(parent) );
	}else{
		RC_TRY( write_node_end(parent) );
	}

	return(NORMAL_RC);
}


/* fullになったleafを分割 */
static RC
split_leaf (MATRIX_DB_INFO *mdbi, MATRIX_DB_LEAF *leaf,
            MATRIX_DB_LEAF *new_leaf)
{
	int ii1, ii2;
	int index;
	int found_flag;
	int file;
	int new_node_flag = 0;
	WRAP64_INT offset;
	MATRIX_DB_NODE *node;

	RC_NULL_CHK( mdbi );
	RC_NULL_CHK( leaf );
	RC_NULL_CHK( new_leaf );
	if(leaf->size < MATRIX_DB_LEAF_SIZE/2) return(ARG_ERROR_RC);
	if(new_leaf->size != 0) return(ARG_ERROR_RC);

	log_printf(6, "[%d] \x1b[32mSplit Leaf\x1b[0m\n", __LINE__);
	/* 新leafへleafの後半を移動 */
	for(ii1=(MATRIX_DB_LEAF_SIZE/2), ii2=0; ii1<leaf->size; ii1++){
		new_leaf->row[ii2] = leaf->row[ii1];
		new_leaf->col[ii2] = leaf->col[ii1];
		RC_TRY( copy_matrix33(leaf->v[ii1], new_leaf->v[ii2]) );
		ii2++;
	}
	new_leaf->size = ii2;

	/* leaf の後半を初期化 */
	leaf->size = MATRIX_DB_LEAF_SIZE/2;
	RC_TRY( clean_leaf(leaf) );

	/* 親node(edge)を探す */
	RC_TRY( search_parent_edge(mdbi, leaf, &file, &offset) );
	if((file <= 0) || (file > MATRIX_DB_SF_SIZE)) return(UNKNOWN_ERROR_RC);
	if(offset < 0) return(UNKNOWN_ERROR_RC);

	RC_TRY( write_node_start(mdbi, file, offset, &node) );

	if(node->size >= MATRIX_DB_NODE_SIZE){
		MATRIX_DB_NODE *new_node;
		RC_TRY( new_node_start(mdbi, &new_node) );
		new_node->type = node->type;
		RC_TRY( split_node(mdbi, node, new_node) );
		if( is_parent_edge(leaf, node) ){
			RC_TRY( new_node_end(new_node) );
		}else if( is_parent_edge(leaf, new_node) ){
			RC_TRY( write_node_end(node) );
			new_node_flag = 1;
			node = new_node;
		}else{
			return(UNKNOWN_ERROR_RC);
		}
	}

	RC_TRY( search_node(node, new_leaf->row[0], new_leaf->col[0],
	                    &index, &found_flag, &file, &offset) );
	if(index < 0) return(UNKNOWN_ERROR_RC);

	/* 親にnew_leafの情報を登録 */
	RC_TRY( insert_key_node(node, index, new_leaf->row[0], new_leaf->col[0],
	                        new_leaf->file_self, new_leaf->offset_self) );
	if(new_node_flag){
		RC_TRY( new_node_end(node) );
	}else{
		RC_TRY( write_node_end(node) );
	}

	return(NORMAL_RC);
}


/* leafの親にあたるindex部末端のnodeを探し，file/offsetを返す */
static RC
search_parent_edge (MATRIX_DB_INFO *mdbi, MATRIX_DB_LEAF *leaf,
                    int *file, WRAP64_INT *offset)
{
	int tmp_file;
	WRAP64_INT tmp_offset;
	MATRIX_DB_NODE *edge_node;

	RC_NULL_CHK( mdbi );
	RC_NULL_CHK( leaf );
	RC_NULL_CHK( file );
	RC_NULL_CHK( offset );
	if(mdbi->route_size < 0) return(ARG_ERROR_RC);
	if(mdbi->route_count > MATRIX_DB_ROUTE_SIZE) return(SEEK_ERROR_RC);

	*file   =  0;
	*offset = -1;

	/* 経路情報が無いのでleafまでの経路を空サーチ */
	if(mdbi->route_size == 0){
		if(leaf->size == 0){
			if((leaf->key_row < 0) || (leaf->key_col < 0))
				return(ARG_ERROR_RC);
			RC_TRY( search_tree_leaf(mdbi, mdbi->root, leaf->key_row,
			                         leaf->key_col, &tmp_file, &tmp_offset) );
		}else{
			RC_TRY( search_tree_leaf(mdbi, mdbi->root, leaf->row[0],
			                         leaf->col[0], &tmp_file, &tmp_offset) );
		}
		if((tmp_file <= 0) || (tmp_file > MATRIX_DB_SF_SIZE))
			return(UNKNOWN_ERROR_RC);
		if(tmp_offset < 0) return(UNKNOWN_ERROR_RC);
		if(mdbi->route_size <= 0) return(UNKNOWN_ERROR_RC);
	}

	*file = mdbi->route_file[mdbi->route_size-1];
	*offset = mdbi->route_offset[mdbi->route_size-1];
	if( (*file < 1)||(*offset < 0) ) return(UNKNOWN_ERROR_RC);

	/* Verification Check */
	RC_TRY( read_node_start(mdbi, *file, *offset, &edge_node) );
	if( !is_edge_node(edge_node) ) return(UNKNOWN_ERROR_RC);
	if( !is_parent_edge(leaf, edge_node) ){
		mdbi->route_count++;
		mdbi->route_size = 0;
		RC_TRY( search_parent_edge(mdbi, leaf, file, offset) );
	}
	RC_TRY( read_node_end(edge_node) );

	return(NORMAL_RC);
}


/* nodeの親にあたる親 nodeを探し，file/offsetを返す */
static RC
search_parent_node (MATRIX_DB_INFO *mdbi, MATRIX_DB_NODE *node,
                    int *file, WRAP64_INT *offset)
{
	int tmp_file;
	WRAP64_INT tmp_offset;
	MATRIX_DB_NODE *parent;

	RC_NULL_CHK( mdbi );
	RC_NULL_CHK( node );
	RC_NULL_CHK( file );
	RC_NULL_CHK( offset );
	if(mdbi->route_size < 0) return(ARG_ERROR_RC);
	if(mdbi->route_count > MATRIX_DB_ROUTE_SIZE) return(SEEK_ERROR_RC);
	if( is_root_node(node) ) return(ARG_ERROR_RC);

	*file   =  0;
	*offset = -1;

	/* 経路情報が無いのでleafまでの経路を空サーチ */
	if(mdbi->route_size == 0){
		if(node->size == 0){
			if((node->key_row < 0) || (node->key_col < 0))
				return(ARG_ERROR_RC);
			RC_TRY( search_tree_node(mdbi, mdbi->root, node->key_row,
			                         node->key_col, &tmp_file, &tmp_offset) );
		}else{
			RC_TRY( search_tree_node(mdbi, mdbi->root, node->row[0],
			                         node->col[0], &tmp_file, &tmp_offset) );
		}
		if((tmp_file <= 0) || (tmp_file > MATRIX_DB_SF_SIZE))
			return(UNKNOWN_ERROR_RC);
		if(tmp_offset < 0) return(UNKNOWN_ERROR_RC);
		if(mdbi->route_size <= 0) return(UNKNOWN_ERROR_RC);
		if( (node->file_self != tmp_file)||(node->offset_self != tmp_offset) ){
			return(UNKNOWN_ERROR_RC);
		}
	}

	*file = mdbi->route_file[mdbi->route_size-2];
	*offset = mdbi->route_offset[mdbi->route_size-2];
	if( (*file < 1)||(*offset < 0) ) return(UNKNOWN_ERROR_RC);

	/* Verification Check */
	RC_TRY( read_node_start(mdbi, *file, *offset, &parent) );
	if( !is_parent_node(node, parent) ){
		mdbi->route_count++;
		mdbi->route_size = 0;
		RC_TRY( search_parent_node(mdbi, node, file, offset) );
	}
	RC_TRY( read_node_end(parent) );

	return(NORMAL_RC);
}


/* index-Tree内から row/col を含むであろう MATRIX_DB_LEAF の
 * file/offsetを取得する．
 */
static RC
search_tree_leaf (MATRIX_DB_INFO *mdbi, MATRIX_DB_NODE *node,
                  int row, int col, int *file, WRAP64_INT *offset)
{
	int index;
	int found_flag;
	MATRIX_DB_NODE *child;

	RC_NULL_CHK( mdbi );
	RC_NULL_CHK( node );
	RC_NULL_CHK( file );
	RC_NULL_CHK( offset );
	if( (row < 0) || (col < 0) ) return(ARG_ERROR_RC);

	*file   =  0;
	*offset = -1;

	/* 経路情報を記録 */
	mdbi->route_file[mdbi->route_size]   = node->file_self;
	mdbi->route_offset[mdbi->route_size] = node->offset_self;
	mdbi->route_size++;
	if(mdbi->route_size >= MATRIX_DB_ROUTE_SIZE) return(OVERFLOW_ERROR_RC);

	log_printf(6, "[%d] row = %d, col = %d\n", __LINE__, row, col);
	RC_TRY( search_node(node, row, col, &index, &found_flag, file, offset) );
	if( (*file < 1) || (*offset < 0) ){
		RC_TRY( print_matrix_db_node(stderr, node) );
		return(UNKNOWN_ERROR_RC);
	}

	if( !is_edge_node(node) ){
		log_printf(6, "[%d] file = %d, offset = %lld\n", __LINE__,
		           *file, *offset);
		RC_TRY( read_node_start(mdbi, *file, *offset, &child) );
		RC_TRY( search_tree_leaf(mdbi, child, row, col, file, offset) );
		RC_TRY( read_node_end(child) );
	}

	return(NORMAL_RC);
}


/* index-Tree内から row/col を含む MATRIX_DB_NODEのfile/offsetを取得する */
static RC
search_tree_node (MATRIX_DB_INFO *mdbi, MATRIX_DB_NODE *node,
                  int row, int col, int *file, WRAP64_INT *offset)
{
	int index;
	int found_flag;
	MATRIX_DB_NODE *child;

	RC_NULL_CHK( mdbi );
	RC_NULL_CHK( node );
	RC_NULL_CHK( file );
	RC_NULL_CHK( offset );
	if( (row < 0) || (col < 0) ) return(ARG_ERROR_RC);

	*file   =  0;
	*offset = -1;

	/* 経路情報を記録 */
	mdbi->route_file[mdbi->route_size]   = node->file_self;
	mdbi->route_offset[mdbi->route_size] = node->offset_self;
	mdbi->route_size++;
	if(mdbi->route_size >= MATRIX_DB_ROUTE_SIZE) return(OVERFLOW_ERROR_RC);

	RC_TRY( search_node(node, row, col, &index, &found_flag, file, offset) );
	if(found_flag){
		*file   = node->file_self;
		*offset = node->offset_self;
		return(NORMAL_RC);
	}else if( is_edge_node(node) ){
		return(END_RC);
	}else{
		if( (*file < 1) || (*offset < 0) ) return(UNKNOWN_ERROR_RC);
		RC_TRY( read_node_start(mdbi, *file, *offset, &child) );
		RC_TRY( search_tree_node(mdbi, child, row, col, file, offset) );
		RC_TRY( read_node_end(child) );
	}

	return(NORMAL_RC);
}


static RC
insert_key_node (MATRIX_DB_NODE *node, int index, int row, int col, int file,
                 WRAP64_INT offset)
{
	int ii1;

	RC_NULL_CHK(node);

	if( (index < 0) || (index > node->size) ) return(ARG_ERROR_RC);
	if( (row < 0) || (col < 0) ) return(ARG_ERROR_RC);
	if( (file <= 0) || (file > MATRIX_DB_SF_SIZE) ) return(ARG_ERROR_RC);
	if(offset < 0)  return(ARG_ERROR_RC);
	/* if(node->size <= 0) return(ARG_ERROR_RC); */
	if(node->size >= MATRIX_DB_NODE_SIZE) return(OVERFLOW_ERROR_RC);

	if(index == node->size){
		node->row[index] = row;
		node->col[index] = col;
		node->file[index+1] = file;
		node->offset[index+1] = offset;
		(node->size)++;
		return(NORMAL_RC);
	}

	for(ii1=node->size; ii1>index; ii1--){
		node->row[ii1] = node->row[ii1-1];
		node->col[ii1] = node->col[ii1-1];
	}
	(node->size)++;
	node->row[index] = row;
	node->col[index] = col;

	for(ii1=node->size; ii1>index; ii1--){
		node->file[ii1]   = node->file[ii1-1];
		node->offset[ii1] = node->offset[ii1-1];
	}
	node->file[index+1]   = file;
	node->offset[index+1] = offset;

	return(NORMAL_RC);
}


static RC
insert_v_leaf (MATRIX_DB_LEAF *leaf, int index, int row, int col,
               double v[][3])
{
	int ii1;

	RC_NULL_CHK(leaf);

	if( (index < 0) || (index > leaf->size) ) return(ARG_ERROR_RC);
	if( (row < 0) || (col < 0) ) return(ARG_ERROR_RC);
	if(leaf->size >= MATRIX_DB_LEAF_SIZE) return(OVERFLOW_ERROR_RC);

	for(ii1=leaf->size; ii1>index; ii1--){
		leaf->row[ii1] = leaf->row[ii1-1];
		leaf->col[ii1] = leaf->col[ii1-1];
		RC_TRY( copy_matrix33(leaf->v[ii1-1], leaf->v[ii1]) );
	}
	leaf->size++;

	leaf->row[index] = row;
	leaf->col[index] = col;
	RC_TRY( copy_matrix33(v, leaf->v[index]) );

	return(NORMAL_RC);
}


RC
delete_v_matrix_db (int row, int col, MATRIX_DB_INFO *mdbi)
{
	int ii1;
	int index, found_flag;
	int file;
	WRAP64_INT offset;
	MATRIX_DB_LEAF *leaf;
	MATRIX_DB_NODE *parent;

	RC_NULL_CHK( mdbi );
	if((row < 0) || (col < 0)) return(ARG_ERROR_RC);

	/* 該当leafを検索 */
	mdbi->route_size = 0;
	RC_TRY( search_tree_leaf(mdbi, mdbi->root, row, col, &file, &offset) );
	if((file <= 0) || (file > MATRIX_DB_SF_SIZE)) return(UNKNOWN_ERROR_RC);
	if(offset < 0)  return(UNKNOWN_ERROR_RC);

	log_printf(6, "[%d] row = %d, col = %d\n", __LINE__, row, col);

	RC_TRY( write_leaf_start(mdbi, file, offset, &leaf) );
	RC_TRY( search_leaf(leaf, row, col, &index, &found_flag) );
	if(found_flag){
		/* index番目を削除 */
		int row0 = leaf->row[0];
		int col0 = leaf->col[0];
		if( (index < 0) || (index >= leaf->size) ){
			DEBUG_PRINT;
			RC_TRY( print_matrix_db_leaf(stderr, leaf) );
			RC_TRY( print_matrix_db_tree(stderr, mdbi) );
			return(UNKNOWN_ERROR_RC);
		}
		for(ii1=index+1; ii1<leaf->size; ii1++){
			leaf->row[ii1-1] = leaf->row[ii1];
			leaf->col[ii1-1] = leaf->col[ii1];
			RC_TRY( copy_matrix33(leaf->v[ii1], leaf->v[ii1-1]) );
		}
		leaf->size--;
		if(leaf->size == 0){
			leaf->key_row = row0;
			leaf->key_col = col0;
		}else{
			RC_TRY( clean_leaf(leaf) );
		}
	}else{
		RC_TRY( log_printf(5, "[%d] row = %d col = %d is missing.\n",
		                   __LINE__, row, col) );
	}

	if(!(leaf->delete_flag) && (leaf->size < MATRIX_DB_LEAF_SIZE/2)){
		/* leafのバランス取り */
		RC_TRY( search_parent_edge(mdbi, leaf, &file, &offset) );
		if((file <= 0) || (file > MATRIX_DB_SF_SIZE)) return(UNKNOWN_ERROR_RC);
		if(offset < 0) return(UNKNOWN_ERROR_RC);
		RC_TRY( write_node_start(mdbi, file, offset, &parent) );
		if(parent->size > 0){
			RC_TRY( balance_leaf(mdbi, parent, leaf) );
			// RC_TRY( print_matrix_db_tree(stderr, mdbi) );
		}
		if((parent->size < MATRIX_DB_NODE_SIZE/2) && (!is_root_node(parent))){
			RC_TRY( balance_node(mdbi, parent) );
			// RC_TRY( print_matrix_db_tree(stderr, mdbi) );
			if(parent->delete_flag){
				if(parent->use_count == 1){
					file   = parent->file_self;
					offset = parent->offset_self;
					RC_TRY( write_node_end(parent) );
					RC_TRY( delete_node(mdbi, file, offset) );
				}else{
					DEBUG_PRINT;
					return(UNKNOWN_ERROR_RC);
				}
			}else{
				RC_TRY( write_node_end(parent) );
			}
		}else{
			RC_TRY( write_node_end(parent) );
		}
	}

	if(leaf->delete_flag){
		if(leaf->use_count == 1){
			file   = leaf->file_self;
			offset = leaf->offset_self;
			RC_TRY( write_leaf_end(leaf) );
			RC_TRY( delete_leaf(mdbi, file, offset) );
		}else{
			RC_TRY( write_leaf_end(leaf) );
		}
	}else{
		RC_TRY( write_leaf_end(leaf) );
	}


#if 0
	if((row == 44) && (col == 49)){
		RC_TRY( print_matrix_db_tree(stderr, mdbi) );
	}
	if((row == 45) && (col == 11)){
		RC_TRY( print_matrix_db_tree(stderr, mdbi) );
		exit(1);
	}
#endif

#if 0
	log_printf(5, "=========================================\n");
	RC_TRY( print_matrix_db_tree(stderr, mdbi) );
	log_printf(5, "=========================================\n\n");
#endif

	return(NORMAL_RC);
}


/* node から file/offset の値を持つ部分を削除する */
static RC
delete_key_node (MATRIX_DB_NODE *node, int file, WRAP64_INT offset)
{
	int ii1;
	int index = -1;

	RC_NULL_CHK( node );
	if((file <= 0) || (file > MATRIX_DB_SF_SIZE)) return(ARG_ERROR_RC);
	if(offset < 0)  return(ARG_ERROR_RC);

	for(ii1=0; ii1<node->size+1; ii1++){
		if( (file == node->file[ii1]) && (offset == node->offset[ii1]) ){
			index = ii1;
			break;
		}
	}

	if(index < 0) return(UNKNOWN_ERROR_RC);
	if(index == 0){
		for(ii1=index+1; ii1<node->size; ii1++){
			node->row[ii1-1] = node->row[ii1];
			node->col[ii1-1] = node->col[ii1];
			node->file[ii1-1]   = node->file[ii1];
			node->offset[ii1-1] = node->offset[ii1];
		}
		node->file[ii1-1]   = node->file[ii1];
		node->offset[ii1-1] = node->offset[ii1];
	}else if(index < node->size){ /* index >= 1 */
		for(ii1=index; ii1<node->size; ii1++){
			node->row[ii1-1] = node->row[ii1];
			node->col[ii1-1] = node->col[ii1];
			node->file[ii1]   = node->file[ii1+1];
			node->offset[ii1] = node->offset[ii1+1];
		}
	}
	node->size--;
	RC_TRY( clean_node(node) );

	return(NORMAL_RC);
}


/* 葉のバランスを取る */
static RC
balance_node (MATRIX_DB_INFO *mdbi, MATRIX_DB_NODE *node)
{
	int ii1;
	int index;
	int file;
	WRAP64_INT offset;
	MATRIX_DB_NODE *parent;

	RC_NULL_CHK( mdbi );
	RC_NULL_CHK( node );
	if( is_root_node(node) ) return(ARG_ERROR_RC);

	log_printf(6, "[%d] \x1b[33mBalance Node\x1b[0m\n", __LINE__);

	mdbi->route_size = 0;
	RC_TRY( search_parent_node(mdbi, node, &file, &offset) );
	if((file <= 0) || (file > MATRIX_DB_SF_SIZE)) return(UNKNOWN_ERROR_RC);
	if(offset < 0) return(UNKNOWN_ERROR_RC);
	RC_TRY( write_node_start(mdbi, file, offset, &parent) );

	// RC_TRY( print_matrix_db_node(stderr, parent) );

	/* nodeを指す親の file[index] の indexを探す．総当り */
	for(ii1=0, index=-1; ii1<parent->size+1; ii1++){
		if( (node->file_self   == parent->file[ii1])
		 && (node->offset_self == parent->offset[ii1]) ){
			index = ii1;
			break;
		}
	}
	log_printf(6, "[%d] index = %d\n", __LINE__, index);
	if(index < 0) return(UNKNOWN_ERROR_RC);

	/* nodeの右とバランス取り */
	if(parent->size > index){
		RC_TRY( adjust_node_right(mdbi, node, parent, index) );
	}

	/* nodeの左とバランス取り */
	if( (node->size < MATRIX_DB_NODE_SIZE/2)
	 && (parent->size > 0) && (parent->delete_flag == 0)
	 && (!is_root_node(node)) && (index > 0) ){
		RC_TRY( adjust_node_left(mdbi, node, parent, index) );
	}

	if(parent->size < MATRIX_DB_NODE_SIZE/2){
		if(!is_root_node(parent)){
			RC_TRY( balance_node(mdbi, parent) );
		}
		if(parent->delete_flag){
			if(parent->use_count == 1){
				file   = parent->file_self;
				offset = parent->offset_self;
				RC_TRY( write_node_end(parent) );
				RC_TRY( delete_node(mdbi, file, offset) );
			}else{
				DEBUG_PRINT;
				return(UNKNOWN_ERROR_RC);
			}
		}else{
			RC_TRY( write_node_end(parent) );
		}
	}else{
		RC_TRY( write_node_end(parent) );
	}

	return(NORMAL_RC);
}


/* nodeの右側とバランス取り．indexはnodeを指しているparentのindex */
static RC
adjust_node_right (MATRIX_DB_INFO *mdbi, MATRIX_DB_NODE *node,
                   MATRIX_DB_NODE *parent, int index)
{
	int ii1, ii2;
	int size;
	int file;
	WRAP64_INT offset;
	MATRIX_DB_NODE *right;

	RC_NULL_CHK( mdbi );
	RC_NULL_CHK( parent );
	RC_NULL_CHK( node );
	if(index < 0) return(ARG_ERROR_RC);

	log_printf(6, "[%d] \x1b[33mBalance Node right\x1b[0m\n", __LINE__);

	RC_TRY( write_node_start(mdbi, parent->file[index+1],
	                               parent->offset[index+1], &right) );
	size = node->size + 1 + right->size;
	if(size <= MATRIX_DB_NODE_SIZE){
		/* 右とマージ */
	   	log_printf(6, "[%d] \x1b[33mBalance Node right merge\x1b[0m\n",
		           __LINE__);
		node->row[node->size] = parent->row[index];
		node->col[node->size] = parent->col[index];
		node->size++;
		for(ii1=node->size, ii2=0; ii1<size; ii1++){
			node->row[ii1] = right->row[ii2];
			node->col[ii1] = right->col[ii2];
			node->file[ii1]   = right->file[ii2];
			node->offset[ii1] = right->offset[ii2];
			ii2++;
		}
		node->file[ii1]   = right->file[ii2];
		node->offset[ii1] = right->offset[ii2];
		node->size = size;
		right->size = 0;
		right->delete_flag = 1;
		if(parent->size == 1){
			parent->key_row = parent->row[0];
			parent->key_col = parent->col[0];
		}
		RC_TRY( delete_key_node(parent, right->file_self, right->offset_self));
		if(right->use_count == 1){
			/* right削除 */
			file   = right->file_self;
			offset = right->offset_self;
			RC_TRY( write_node_end(right) );
			RC_TRY( delete_node(mdbi, file, offset) );
		}else{
			DEBUG_PRINT;
			return(UNKNOWN_ERROR_RC);
		}
		if(is_root_node(parent) && (parent->size == 0)){
			if(is_edge_node(node)){
				RC_TRY( set_root_node(mdbi, node, MATRIX_DB_ROOT_EDGE) );
			}else{
				RC_TRY( set_root_node(mdbi, node, MATRIX_DB_ROOT) );
			}
			parent->delete_flag = 1;
		}
	}else{
	   log_printf(6,"[%d] \x1b[33mBalance Node right move\x1b[0m\n",__LINE__);
		/* nodeとrightがほぼ同じサイズになるように分配           */
		/* right->parent->nodeへと一つずつ移動するのを一度にやる */
		int move_size = (size-1)/2 - node->size;
		if(move_size < 0){
			return(UNKNOWN_ERROR_RC);
		}
		node->row[node->size] = parent->row[index];
		node->col[node->size] = parent->col[index];
		node->size++;
		node->file[node->size]   = right->file[0];
		node->offset[node->size] = right->offset[0];
		for(ii1=0; ii1<move_size-1; ii1++){
			node->row[node->size] = right->row[ii1];
			node->col[node->size] = right->col[ii1];
			node->file[node->size+1]   = right->file[ii1+1];
			node->offset[node->size+1] = right->offset[ii1+1];
			node->size++;
		}
		parent->row[index] = right->row[ii1];
		parent->col[index] = right->col[ii1];
		for(ii1=move_size, ii2=0; ii1<right->size; ii1++){
			right->row[ii2] = right->row[ii1];
			right->col[ii2] = right->col[ii1];
			right->file[ii2]   = right->file[ii1];
			right->offset[ii2] = right->offset[ii1];
			ii2++;
		}
		right->file[ii2]   = right->file[ii1];
		right->offset[ii2] = right->offset[ii1];
		right->size -= move_size;
		RC_TRY( clean_node(right) );
		RC_TRY( write_node_end(right) );
	}

	return(NORMAL_RC);
}


/* nodeの左側とバランス取り．indexは nodeを指している parentの index */
static RC
adjust_node_left (MATRIX_DB_INFO *mdbi, MATRIX_DB_NODE *node,
                  MATRIX_DB_NODE *parent, int index)
{
	int ii1, ii2;
	MATRIX_DB_NODE *left;
	int size;

	RC_NULL_CHK( mdbi );
	RC_NULL_CHK( parent );
	RC_NULL_CHK( node );
	if(index < 0) return(ARG_ERROR_RC);

	log_printf(6, "[%d] \x1b[33mBalance Node left\x1b[0m\n", __LINE__);

	RC_TRY( write_node_start(mdbi, parent->file[index-1],
	                               parent->offset[index-1], &left) );
	size = left->size + 1 + node->size;
	if(size <= MATRIX_DB_NODE_SIZE){
		log_printf(6,"[%d] \x1b[33mBalance Node left merge\x1b[0m\n",__LINE__);
		left->row[left->size] = parent->row[index-1];
		left->col[left->size] = parent->col[index-1];
		left->size++;
		for(ii1=left->size, ii2=0; ii1<size; ii1++){
			left->row[ii1] = node->row[ii2];
			left->col[ii1] = node->col[ii2];
			left->file[ii1]   = node->file[ii2];
			left->offset[ii1] = node->offset[ii2];
			ii2++;
		}
		left->file[ii1]   = node->file[ii2];
		left->offset[ii1] = node->offset[ii2];
		left->size = size;
		node->size = 0;
		node->delete_flag = 1; /* nodeは上位の関数で消す */
		node->key_row = node->row[0];
		node->key_col = node->col[0];
		if(parent->size == 1){
			parent->key_row = parent->row[0];
			parent->key_col = parent->col[0];
		}
		RC_TRY( delete_key_node(parent, node->file_self,
		                                node->offset_self) );
		if(is_root_node(parent) && (parent->size == 0)){
			if(is_edge_node(left)){
				RC_TRY( set_root_node(mdbi, left, MATRIX_DB_ROOT_EDGE) );
			}else{
				RC_TRY( set_root_node(mdbi, left, MATRIX_DB_ROOT) );
			}
			parent->delete_flag = 1;
		}
		RC_TRY( write_node_end(left) );
	}else{
		log_printf(6,"[%d] \x1b[33mBalance Node left move\x1b[0m\n",__LINE__);
		/* nodeとleftがほぼ同じサイズになるように分配           */
		/* left->parent->nodeへと一つずつ移動するのを一度にやる */
		int move_size = (size-1)/2 - node->size;
		if(move_size < 0){
			return(UNKNOWN_ERROR_RC);
		}
		for(ii1=(node->size-1); ii1>=0; ii1--){
			node->row[move_size+ii1] = node->row[ii1];
			node->col[move_size+ii1] = node->col[ii1];
			node->file[move_size+ii1+1]   = node->file[ii1+1];
			node->offset[move_size+ii1+1] = node->offset[ii1+1];
		}
		node->file[move_size+ii1+1]   = node->file[ii1+1];
		node->offset[move_size+ii1+1] = node->offset[ii1+1];
		node->size += move_size;

		node->row[move_size-1] = parent->row[index-1];
		node->col[move_size-1] = parent->col[index-1];
		node->file[move_size-1]   = left->file[left->size];
		node->offset[move_size-1] = left->offset[left->size];

		for(ii1=0, ii2=(left->size-move_size+1); ii1<move_size-1; ii1++){
			node->row[ii1] = left->row[ii2+ii1];
			node->col[ii1] = left->col[ii2+ii1];
			node->file[ii1]   = left->file[ii2+ii1];
			node->offset[ii1] = left->offset[ii2+ii1];
		}
		parent->row[index-1] = left->row[ii2-1];
		parent->col[index-1] = left->col[ii2-1];
		left->size -= move_size;
		RC_TRY( clean_node(left) );

		RC_TRY( write_node_end(left) );
	}

	return(NORMAL_RC);
}


/* 葉のバランスを取る */
static RC
balance_leaf (MATRIX_DB_INFO *mdbi, MATRIX_DB_NODE *parent,
              MATRIX_DB_LEAF *leaf)
{
	int ii1, ii2;
	int index;
	int file;
	WRAP64_INT offset;

	RC_NULL_CHK( mdbi );
	RC_NULL_CHK( parent );
	RC_NULL_CHK( leaf );
	if(parent->size <= 0) return(ARG_ERROR_RC);

	log_printf(6, "[%d] \x1b[32mBalance Leaf\x1b[0m\n", __LINE__);

	/* leafを指す親の file[index] の indexを探す．総当り */
	for(ii1=0, index=-1; ii1<parent->size+1; ii1++){
		if( (leaf->file_self   == parent->file[ii1])
		 && (leaf->offset_self == parent->offset[ii1]) ){
			index = ii1;
			break;
		}
	}
	if(index < 0) return(UNKNOWN_ERROR_RC);

	if(parent->size > index){
		MATRIX_DB_LEAF *right;
		int size;
		RC_TRY( write_leaf_start(mdbi, parent->file[index+1],
		                               parent->offset[index+1], &right) );
		size = leaf->size + right->size;
		if(size <= MATRIX_DB_LEAF_SIZE){
			for(ii1=leaf->size, ii2=0; ii1<size; ii1++){
				leaf->row[ii1] = right->row[ii2];
				leaf->col[ii1] = right->col[ii2];
				RC_TRY( copy_matrix33(right->v[ii2], leaf->v[ii1]) );
				ii2++;
			}
			leaf->size = size;
			right->size = 0;
			right->delete_flag = 1;
			if(parent->size == 1){
				parent->key_row = parent->row[0];
				parent->key_col = parent->col[0];
			}
			RC_TRY( delete_key_node(parent, right->file_self,
			                                right->offset_self) ) ;
			if(right->use_count == 1){
				/* right削除 */
				file   = right->file_self;
				offset = right->offset_self;
				RC_TRY( write_leaf_end(right) );
				RC_TRY( delete_leaf(mdbi, file, offset) );
				return(NORMAL_RC);
			}
		}
		RC_TRY( write_leaf_end(right) );
	}else if(index > 0){
		MATRIX_DB_LEAF *Left;
		int size;
		RC_TRY( write_leaf_start(mdbi, parent->file[index-1],
		                               parent->offset[index-1], &Left) );
		size = Left->size + leaf->size;
		if(size <= MATRIX_DB_LEAF_SIZE){
			for(ii1=Left->size, ii2=0; ii1<size; ii1++){
				Left->row[ii1] = leaf->row[ii2];
				Left->col[ii1] = leaf->col[ii2];
				RC_TRY( copy_matrix33(leaf->v[ii2], Left->v[ii1]) );
				ii2++;
			}
			Left->size = size;
			leaf->size = 0;
			leaf->delete_flag = 1;
			if(parent->size == 1){
				parent->key_row = parent->row[0];
				parent->key_col = parent->col[0];
			}
			RC_TRY( delete_key_node(parent, leaf->file_self,
			                                leaf->offset_self) ) ;
		}
		RC_TRY( write_leaf_end(Left) );
	}else{
		DEBUG_PRINT;
		return(UNKNOWN_ERROR_RC);
	}

	return(NORMAL_RC);
}


/*
 * cache 内で db_file, db_offset のデータを探してそのポインタを返却
 * 無かったら NULL を返却
 * アクセスカウンタの調整もここで行う．
 * inc_access : 0 ならアクセスカウンタに小さな値をセット
 *              1 ならアクセスカウンタをリセット
 */
static MATRIX_DB_NODE *
search_node_cache_ptr (MATRIX_DB_CACHE_NODE *ca_node, int file,
                       WRAP64_INT offset, int inc_access)
{
	int ii1, ii2;
	int crow = (int)(CACHE_ROW(offset, ca_node->divisor,
	                           sizeof(MATRIX_DB_NODE)) );

	for(ii1=0; ii1<ca_node->associate; ii1++){
		if(ca_node->cache[crow][ii1].use_count < 0){
			return(NULL);
		}else if( (ca_node->cache[crow][ii1].file_self == file)
		       && (ca_node->cache[crow][ii1].offset_self == offset) ){
			if(inc_access){
				if(ca_node->access_count[crow][ii1] + MAX_ACCESS_COUNT/32
				   >= MAX_ACCESS_COUNT){
					for(ii2=0; ii2<ca_node->associate; ii2++){
						(ca_node->access_count[crow][ii2]) /= 2;
					}
				}
				ca_node->access_count[crow][ii1] += MAX_ACCESS_COUNT/32;
			}else{
				ca_node->access_count[crow][ii1] = 0;
			}
			return( &(ca_node->cache[crow][ii1]) );
		}
	}

	return(NULL);
}


/* search_node_cache_ptr() と同じ */
static MATRIX_DB_LEAF *
search_leaf_cache_ptr (MATRIX_DB_CACHE_LEAF *ca_leaf, int file,
                       WRAP64_INT offset, int inc_access)
{
	int ii1, ii2;
	int crow = (int)(CACHE_ROW(offset, ca_leaf->divisor,
	                           sizeof(MATRIX_DB_LEAF)) );

	for(ii1=0; ii1<ca_leaf->associate; ii1++){
		if(ca_leaf->cache[crow][ii1].use_count < 0){
			return(NULL);
		}else if( (ca_leaf->cache[crow][ii1].file_self == file)
		       && (ca_leaf->cache[crow][ii1].offset_self == offset) ){
			if(inc_access){
				if(ca_leaf->access_count[crow][ii1] + MAX_ACCESS_COUNT/32
				   >= MAX_ACCESS_COUNT){
					for(ii2=0; ii2<ca_leaf->associate; ii2++){
						(ca_leaf->access_count[crow][ii2]) /= 2;
					}
				}
				ca_leaf->access_count[crow][ii1] += MAX_ACCESS_COUNT/32;
			}else{
				ca_leaf->access_count[crow][ii1] = 0;
			}
			return( &(ca_leaf->cache[crow][ii1]) );
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
static MATRIX_DB_NODE *
unused_node_cache_ptr (MATRIX_DB_INFO *mdbi, int crow)
{
	int ii1;
	int min_count = INT_MAX;
	int min_col = -1;
	MATRIX_DB_NODE *ret;

	for(ii1=0; ii1<mdbi->cache_node.associate; ii1++){
		if(mdbi->cache_node.cache[crow][ii1].use_count < 0){
			/* 未使用領域 */
			min_count = -1;
			min_col = ii1;
			break;
		}else if(mdbi->cache_node.cache[crow][ii1].use_count == 0){
			/* 使用されているがアクセス中ではない */
			int access = mdbi->cache_node.access_count[crow][ii1];
			if(mdbi->cache_node.cache[crow][ii1].dirty_flag){
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

	if(min_col < 0){
		print_matrix_db_cache_node_info(stderr, &mdbi->cache_node);
		return(NULL);  /* 全部アクセス中 */
	}

	ret = &(mdbi->cache_node.cache[crow][min_col]);
	if(ret->dirty_flag){
		if( put_node(&(mdbi->sf_node), ret->file_self, ret->offset_self, ret)
		    != NORMAL_RC) return(NULL);
		ret->dirty_flag = 0;
	}

	/* 新規登録されたデータがすぐに破棄されないように，*/
	/* アクセスカウンタは最大値にする */
	mdbi->cache_node.access_count[crow][min_col] = MAX_ACCESS_COUNT;

	return(ret);
}


/* unused_node_cache_ptr() と同じ */
static MATRIX_DB_LEAF *
unused_leaf_cache_ptr (MATRIX_DB_INFO *mdbi, int crow)
{
	int ii1;
	int min_count = INT_MAX;
	int min_col = -1;
	MATRIX_DB_LEAF *ret;

	for(ii1=0; ii1<mdbi->cache_leaf.associate; ii1++){
		if(mdbi->cache_leaf.cache[crow][ii1].use_count < 0){
			/* 未使用領域 */
			min_count = -1;
			min_col = ii1;
			break;
		}else if(mdbi->cache_leaf.cache[crow][ii1].use_count == 0){
			/* 使用されているがアクセス中ではない */
			int access = mdbi->cache_leaf.access_count[crow][ii1];
			if(mdbi->cache_leaf.cache[crow][ii1].dirty_flag){
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

	if(min_col < 0){
		print_matrix_db_cache_leaf_info(stderr, &mdbi->cache_leaf);
		return(NULL);  /* 全部アクセス中 */
	}

	ret = &(mdbi->cache_leaf.cache[crow][min_col]);
	if(ret->dirty_flag){
		if(put_leaf(&(mdbi->sf_leaf), ret->file_self,
		                ret->offset_self, ret) != NORMAL_RC) return(NULL);
		ret->dirty_flag = 0;
	}

	/* 新規登録されたデータがすぐに破棄されないように，*/
	/* アクセスカウンタは最大値にする */
	mdbi->cache_leaf.access_count[crow][min_col] = MAX_ACCESS_COUNT;

	return(ret);
}


static RC
new_node_start (MATRIX_DB_INFO *mdbi, MATRIX_DB_NODE **node)
{
	int ii1;
	int file = 0;
	WRAP64_INT offset = 0;
	int reuse_flag = 0;
	int crow;

	RC_NULL_CHK( node );
	RC_NULL_CHK( mdbi );
	*node = NULL;

	/* 直近に削除された領域があれば，file/offset に取得 */
	for(ii1=0; ii1<MATRIX_DB_SF_SIZE; ii1++){
		if(mdbi->sf_node.blank_offset[ii1] < 0) continue;

		file = ii1 + 1;
		offset = mdbi->sf_node.blank_offset[ii1];
		reuse_flag = 1;

		/* 当該データはキャッシュ上にあるかも */
		*node = search_node_cache_ptr(&(mdbi->cache_node), file, offset, 1);
		if(*node != NULL){
			mdbi->sf_node.blank_offset[ii1] = (*node)->offset[0];
			break;
		}

		/* 先に削除された領域のオフセットを読み込んで blank_offset へ */
		crow = (int)(CACHE_ROW(offset, mdbi->cache_node.divisor,
		                       sizeof(MATRIX_DB_NODE)));
		RC_NULL_CHK( *node = unused_node_cache_ptr(mdbi, crow) );
		RC_TRY( get_node(&(mdbi->sf_node), file, offset, *node) );
		mdbi->sf_node.blank_offset[ii1] = (*node)->offset[0];
		break;
	}

	if(file <= 0){
		/* 新規に作成して，file/offset をセット */
		/* ここで，スクラッチファイルを選択するが，*/
		/* 今のところスクラッチファイルが一つなので file = 1 に決め打ち */
		file = 1;
		offset = mdbi->sf_node.last_offset[0];
		mdbi->sf_node.last_offset[0] += sizeof(MATRIX_DB_NODE);
		crow = (int)(CACHE_ROW(offset, mdbi->cache_node.divisor,
		                       sizeof(MATRIX_DB_NODE)));
		RC_NULL_CHK( *node = unused_node_cache_ptr(mdbi, crow) );
	}

	/* キャッシュ内にデータを作成 */
	RC_NULL_CHK( *node );    /* 念のため */
	INIT_MATRIX_DB_NODE(**node);
	(*node)->use_count = 1;
	(*node)->dirty_flag = 1;
	(*node)->file_self = file;
	(*node)->offset_self = offset;
	if(reuse_flag) return(NORMAL_RC);

	/* キャッシュのデータをファイルに書き込み */
	RC_TRY( put_node(&(mdbi->sf_node), file, offset, *node) );

	return(NORMAL_RC);
}


static RC
new_leaf_start (MATRIX_DB_INFO *mdbi, MATRIX_DB_LEAF **leaf)
{
	int ii1;
	int file = 0;
	WRAP64_INT offset = 0;
	int reuse_flag = 0;
	int crow;

	RC_NULL_CHK( leaf );
	RC_NULL_CHK( mdbi );
	*leaf = NULL;

	/* 直近に削除された領域があれば，file/offset に取得 */
	for(ii1=0; ii1<MATRIX_DB_SF_SIZE; ii1++){
		if(mdbi->sf_leaf.blank_offset[ii1] < 0) continue;

		file = ii1 + 1;
		offset = mdbi->sf_leaf.blank_offset[ii1];
		reuse_flag = 1;

		/* 当該データはキャッシュ上にあるかも */
		*leaf = search_leaf_cache_ptr(&(mdbi->cache_leaf), file, offset, 1);
		if(*leaf != NULL){
			mdbi->sf_leaf.blank_offset[ii1] = (*leaf)->deleted_offset;
			break;
		}

		/* 先に削除された領域のオフセットを読み込んで blank_offset へ */
		crow = (int)(CACHE_ROW(offset, mdbi->cache_leaf.divisor,
		                       sizeof(MATRIX_DB_LEAF)));
		RC_NULL_CHK( *leaf = unused_leaf_cache_ptr(mdbi, crow) );
		RC_TRY( get_leaf(&(mdbi->sf_leaf), file, offset, *leaf) );
		mdbi->sf_leaf.blank_offset[ii1] = (*leaf)->deleted_offset;
		break;
	}

	if(file <= 0){
		/* 新規に作成して，file/offset をセット */
		/* ここで，スクラッチファイルを選択するが，*/
		/* 今のところスクラッチファイルが一つなので file = 1 に決め打ち */
		file = 1;
		offset = mdbi->sf_leaf.last_offset[0];
		mdbi->sf_leaf.last_offset[0] += sizeof(MATRIX_DB_LEAF);
		crow = (int)(CACHE_ROW(offset, mdbi->cache_leaf.divisor,
		                       sizeof(MATRIX_DB_LEAF)));
		RC_NULL_CHK( *leaf = unused_leaf_cache_ptr(mdbi, crow) );
	}

	/* キャッシュ内にデータを作成 */
	RC_NULL_CHK( *leaf );    /* 念のため */
	INIT_MATRIX_DB_LEAF(**leaf);
	(*leaf)->use_count = 1;
	(*leaf)->dirty_flag = 1;
	(*leaf)->file_self = file;
	(*leaf)->offset_self = offset;
	if(reuse_flag) return(NORMAL_RC);

	/* キャッシュのデータをファイルに書き込み */
	RC_TRY( put_leaf(&(mdbi->sf_leaf), file, offset, *leaf) );

	return(NORMAL_RC);
}


static RC
read_node_start (MATRIX_DB_INFO *mdbi, int file, WRAP64_INT offset,
                 MATRIX_DB_NODE **node)
{
	int crow;

	RC_NULL_CHK( node );
	RC_NULL_CHK( mdbi );
	if( (file <= 0) || (file > MATRIX_DB_SF_SIZE) ) return(ARG_ERROR_RC);
	if(offset < 0)  return(ARG_ERROR_RC);

	/* キャッシュ内を探す */
	*node = search_node_cache_ptr(&(mdbi->cache_node), file, offset, 1);
	if(*node != NULL){
		(*node)->use_count++;
		return(NORMAL_RC);
	}

	/* キャッシュの空き領域に読み込む */
	crow  = (int)(CACHE_ROW(offset, mdbi->cache_node.divisor,
	                        sizeof(MATRIX_DB_NODE)) );
	RC_NULL_CHK( *node = unused_node_cache_ptr(mdbi, crow) );
	RC_TRY( get_node(&(mdbi->sf_node), file, offset, *node) );

	if((*node)->use_count < 0){
		return(UNKNOWN_ERROR_RC);   /* 使われていない？ */
	}
	(*node)->use_count = 1;
	(*node)->dirty_flag = 0;

	return(NORMAL_RC);
}


static RC
read_leaf_start (MATRIX_DB_INFO *mdbi, int file, WRAP64_INT offset,
                 MATRIX_DB_LEAF **leaf)
{
	int crow;

	RC_NULL_CHK( leaf );
	RC_NULL_CHK( mdbi );
	if( (file <= 0) || (file > MATRIX_DB_SF_SIZE) ) return(ARG_ERROR_RC);
	if(offset < 0)  return(ARG_ERROR_RC);

	/* キャッシュ内を探す */
	*leaf = search_leaf_cache_ptr(&(mdbi->cache_leaf), file, offset, 1);
	if(*leaf != NULL){
		(*leaf)->use_count++;
		return(NORMAL_RC);
	}

	/* キャッシュの空き領域に読み込む */
	crow  = (int)(CACHE_ROW(offset, mdbi->cache_leaf.divisor,
	                        sizeof(MATRIX_DB_LEAF)) );
	RC_NULL_CHK( *leaf = unused_leaf_cache_ptr(mdbi, crow) );
	RC_TRY( get_leaf(&(mdbi->sf_leaf), file, offset, *leaf) );

	if((*leaf)->use_count < 0){
		return(UNKNOWN_ERROR_RC);   /* 使われていない？ */
	}
	(*leaf)->use_count = 1;
	(*leaf)->dirty_flag = 0;

	return(NORMAL_RC);
}


static RC
write_node_start (MATRIX_DB_INFO *mdbi, int file, WRAP64_INT offset,
                  MATRIX_DB_NODE **node)
{
	RC_TRY( read_node_start(mdbi, file, offset, node) );
	(*node)->dirty_flag = 1;

	return(NORMAL_RC);
}
	
static RC
write_leaf_start (MATRIX_DB_INFO *mdbi, int file, WRAP64_INT offset,
                  MATRIX_DB_LEAF **leaf)
{
	RC_TRY( read_leaf_start(mdbi, file, offset, leaf) );
	(*leaf)->dirty_flag = 1;

	return(NORMAL_RC);
}


static RC
new_node_end (MATRIX_DB_NODE *node)
{
	RC_NULL_CHK( node );
	if(node->use_count <= 0) return(ARG_ERROR_RC);
	node->use_count--;

	return(NORMAL_RC);
}


static RC
new_leaf_end (MATRIX_DB_LEAF *leaf)
{
	RC_NULL_CHK( leaf );
	if(leaf->use_count <= 0) return(ARG_ERROR_RC);
	leaf->use_count--;

	return(NORMAL_RC);
}


static RC
read_node_end (MATRIX_DB_NODE *node)
{
	RC_NULL_CHK( node );
	if(node->use_count <= 0) return(ARG_ERROR_RC);
	node->use_count--;

	return(NORMAL_RC);
}


static RC
read_leaf_end (MATRIX_DB_LEAF *leaf)
{
	RC_NULL_CHK( leaf );
	if(leaf->use_count <= 0) return(ARG_ERROR_RC);
	leaf->use_count--;

	return(NORMAL_RC);
}


static RC
write_node_end (MATRIX_DB_NODE *node)
{
	RC_NULL_CHK( node );
	if(node->use_count <= 0) return(ARG_ERROR_RC);
	node->use_count--;

	return(NORMAL_RC);
}


static RC
write_leaf_end (MATRIX_DB_LEAF *leaf)
{
	RC_NULL_CHK( leaf );
	if(leaf->use_count <= 0) return(ARG_ERROR_RC);
	leaf->use_count--;

	return(NORMAL_RC);
}



/* キャッシュとファイル上のデータ(file/offset)を削除(再利用可能に) */
static RC
delete_node (MATRIX_DB_INFO *mdbi, int file, WRAP64_INT offset)
{
	MATRIX_DB_NODE tmp_db;
	MATRIX_DB_NODE *ptr;

	RC_NULL_CHK( mdbi );
	if((file <= 0)||(file > MATRIX_DB_SF_SIZE)) return(ARG_ERROR_RC);

	/* キャッシュ内にデータがあれば初期化し，offset[0] */
	/* に直近に削除されたデータのオフセットを書き込み  */
	ptr = search_node_cache_ptr(&(mdbi->cache_node), file, offset, 0);
	if(ptr != NULL){
		if(ptr->use_count > 0){
			log_printf(5, "[%d] use_count > 0\n", __LINE__);
			return(UNKNOWN_ERROR_RC);
		}
		INIT_MATRIX_DB_NODE(*ptr);
		ptr->use_count = 0;
		ptr->offset[0] = mdbi->sf_node.blank_offset[file - 1];
		ptr->file_self = file;
		ptr->offset_self = offset;
		ptr->dirty_flag = 1;

		/* blkank_offset を更新 */
		mdbi->sf_node.blank_offset[file - 1] = offset;

		return(NORMAL_RC);
	}

	/* offset[0] に直近に削除されたデータのオフセットを書き込み */
	INIT_MATRIX_DB_NODE(tmp_db);
	tmp_db.offset[0] = mdbi->sf_node.blank_offset[file - 1];
	tmp_db.file_self = file;
	tmp_db.offset_self = offset;
	RC_TRY( put_node(&(mdbi->sf_node), file, offset, &tmp_db) );

	/* blkank_offset を更新 */
	mdbi->sf_node.blank_offset[file - 1] = offset;

	return(NORMAL_RC);
}


/* キャッシュとファイル上のデータ(file/offset)を削除(再利用可能に) */
static RC
delete_leaf (MATRIX_DB_INFO *mdbi, int file, WRAP64_INT offset)
{
	MATRIX_DB_LEAF tmp_db;
	MATRIX_DB_LEAF *ptr;

	RC_NULL_CHK( mdbi );
	if((file <= 0)||(file > MATRIX_DB_SF_SIZE)) return(ARG_ERROR_RC);

	/* キャッシュ内にデータがあれば初期化し，offset[0] */
	/* に直近に削除されたデータのオフセットを書き込み  */
	ptr = search_leaf_cache_ptr(&(mdbi->cache_leaf), file, offset, 0);
	if(ptr != NULL){
		if(ptr->use_count > 0){
			log_printf(5, "[%d] use_count > 0\n", __LINE__);
			return(UNKNOWN_ERROR_RC);
		}
		INIT_MATRIX_DB_LEAF(*ptr);
		ptr->use_count = 0;
		ptr->deleted_offset = mdbi->sf_leaf.blank_offset[file - 1];
		ptr->file_self = file;
		ptr->offset_self = offset;
		ptr->dirty_flag = 1;

		/* blkank_offset を更新 */
		mdbi->sf_leaf.blank_offset[file - 1] = offset;

		return(NORMAL_RC);
	}

	/* offset[0] に直近に削除されたデータのオフセットを書き込み */
	INIT_MATRIX_DB_LEAF(tmp_db);
	tmp_db.deleted_offset = mdbi->sf_leaf.blank_offset[file - 1];
	tmp_db.file_self = file;
	tmp_db.offset_self = offset;
	RC_TRY( put_leaf(&(mdbi->sf_leaf), file, offset, &tmp_db) );

	/* blkank_offset を更新 */
	mdbi->sf_leaf.blank_offset[file - 1] = offset;

	return(NORMAL_RC);
}


static RC
put_node (MATRIX_DB_SCRATCH *mdbsf, int file, WRAP64_INT offset,
          MATRIX_DB_NODE *ptr)
{
	SCRATCH_FILE *sf;
	FILE *fp;

	RC_NULL_CHK( mdbsf );
	RC_NULL_CHK( ptr );
	if( (file <= 0) || (file > MATRIX_DB_SF_SIZE) ) return(ARG_ERROR_RC);
	if(offset < 0)  return(ARG_ERROR_RC);

	ptr->file_self = file;
	ptr->offset_self = offset;

	RC_NULL_CHK( sf = mdbsf->scratch_array[file - 1] );
	RC_NULL_CHK( fp = sf->fp );

	if(offset != mdbsf->current_offset[file - 1]){
		/* 現在の読み書き位置と異なる場合だけシークする */
		/* wrap64_fseek() 内でバッファがフラッシュされるのを極力避けるため */
		RC_TRY( wrap64_fseek(fp, offset, SEEK_SET) );
	}
	if(fwrite(ptr, sizeof(MATRIX_DB_NODE), 1, fp) != 1){
		return(WRITE_ERROR_RC);    /* ディスク書き込み失敗 */
	}
	mdbsf->current_offset[file - 1] = offset + sizeof(MATRIX_DB_NODE);

	return(NORMAL_RC);
}


static RC
put_leaf (MATRIX_DB_SCRATCH *mdbsf, int file, WRAP64_INT offset,
          MATRIX_DB_LEAF *ptr)
{
	SCRATCH_FILE *sf;
	FILE *fp;

	RC_NULL_CHK( mdbsf );
	RC_NULL_CHK( ptr );
	if( (file <= 0) || (file > MATRIX_DB_SF_SIZE) ) return(ARG_ERROR_RC);
	if(offset < 0)  return(ARG_ERROR_RC);

	ptr->file_self = file;
	ptr->offset_self = offset;

	RC_NULL_CHK( sf = mdbsf->scratch_array[file - 1] );
	RC_NULL_CHK( fp = sf->fp );

	if(offset != mdbsf->current_offset[file - 1]){
		/* 現在の読み書き位置と異なる場合だけシークする */
		/* wrap64_fseek() 内でバッファがフラッシュされるのを極力避けるため */
		RC_TRY( wrap64_fseek(fp, offset, SEEK_SET) );
	}
	if(fwrite(ptr, sizeof(MATRIX_DB_LEAF), 1, fp) != 1){
		return(WRITE_ERROR_RC);    /* ディスク書き込み失敗 */
	}
	mdbsf->current_offset[file - 1] = offset + sizeof(MATRIX_DB_LEAF);

	return(NORMAL_RC);
}


static RC
get_node (MATRIX_DB_SCRATCH *mdbsf, int file, WRAP64_INT offset,
          MATRIX_DB_NODE *ptr)
{
	SCRATCH_FILE *sf;
	FILE *fp;

	RC_NULL_CHK( mdbsf);
	RC_NULL_CHK( ptr );
	if( (file <= 0) || (file > MATRIX_DB_SF_SIZE) ) return(ARG_ERROR_RC);
	if(offset < 0)  return(ARG_ERROR_RC);

	RC_NULL_CHK( sf = mdbsf->scratch_array[file - 1] );
	RC_NULL_CHK( fp = sf->fp );
	if(offset != mdbsf->current_offset[file - 1]){
		/* 現在の読み書き位置と異なる場合だけシークする */
		/* wrap64_fseek() 内でバッファがフラッシュされるのを極力避けるため */
		RC_TRY( wrap64_fseek(fp, offset, SEEK_SET) );
	}
	if(fread(ptr, sizeof(MATRIX_DB_NODE), 1, fp) != 1){
		return(READ_ERROR_RC);    /* ディスク読み込み失敗 */
	}
	mdbsf->current_offset[file - 1] = offset + sizeof(MATRIX_DB_NODE);

	return(NORMAL_RC);
}


static RC
get_leaf (MATRIX_DB_SCRATCH *mdbsf, int file, WRAP64_INT offset,
          MATRIX_DB_LEAF *ptr)
{
	SCRATCH_FILE *sf;
	FILE *fp;

	RC_NULL_CHK( mdbsf);
	RC_NULL_CHK( ptr );
	if( (file <= 0) || (file > MATRIX_DB_SF_SIZE) ) return(ARG_ERROR_RC);
	if(offset < 0)  return(ARG_ERROR_RC);

	RC_NULL_CHK( sf = mdbsf->scratch_array[file - 1] );
	RC_NULL_CHK( fp = sf->fp );
	if(offset != mdbsf->current_offset[file - 1]){
		/* 現在の読み書き位置と異なる場合だけシークする */
		/* wrap64_fseek() 内でバッファがフラッシュされるのを極力避けるため */
		RC_TRY( wrap64_fseek(fp, offset, SEEK_SET) );
	}
	if(fread(ptr, sizeof(MATRIX_DB_LEAF), 1, fp) != 1){
		return(READ_ERROR_RC);    /* ディスク読み込み失敗 */
	}
	mdbsf->current_offset[file - 1] = offset + sizeof(MATRIX_DB_LEAF);

	return(NORMAL_RC);
}


static RC
set_root_node (MATRIX_DB_INFO *mdbi, MATRIX_DB_NODE *root,
               MATRIX_DB_NODE_TYPE type)
{
	RC_NULL_CHK(mdbi);
	RC_NULL_CHK(root);

	RC_TRY( read_node_end(mdbi->root) );
	root->type = type;
	mdbi->root = root;
	RC_TRY( read_node_start(mdbi, root->file_self, root->offset_self,
	                        &(mdbi->root)) );
	return(NORMAL_RC);
}


static RC
init_node (MATRIX_DB_NODE *node)
{
	int ii1;

	RC_NULL_CHK(node);

	node->size = 0;
	node->dirty_flag = 0;
	node->delete_flag = 0;
	node->use_count = -1;
	node->key_row = -1;
	node->key_col = -1;
	node->type = MATRIX_DB_INIT;
	node->file_self = 0;
	node->offset_self = 0;

	for(ii1=0; ii1<MATRIX_DB_NODE_HEADER; ii1++){
		node->header[ii1] = 0;
	}

	for(ii1=0; ii1<MATRIX_DB_NODE_SIZE+1; ii1++){
		node->file[ii1] = 0;
		node->offset[ii1] = 0;
	}

	for(ii1=0; ii1<MATRIX_DB_NODE_SIZE; ii1++){
		node->row[ii1] = 0;
		node->col[ii1] = 0;
	}

	return(NORMAL_RC);
}


static RC
init_leaf (MATRIX_DB_LEAF *leaf)
{
	int ii1;

	RC_NULL_CHK(leaf);

	leaf->size = 0;
	leaf->dirty_flag = 0;
	leaf->delete_flag = 0;
	leaf->use_count = -1;
	leaf->key_row = -1;
	leaf->key_col = -1;
	leaf->file_self = 0;
	leaf->offset_self = 0;
	leaf->deleted_offset = 0;

	for(ii1=0; ii1<MATRIX_DB_LEAF_HEADER; ii1++){
		leaf->header[ii1] = 0;
	}

	for(ii1=0; ii1<MATRIX_DB_LEAF_SIZE; ii1++){
		leaf->row[ii1] = 0;
		leaf->col[ii1] = 0;
		init_matrix33(leaf->v[ii1]);
	}

	return(NORMAL_RC);
}


/* size 以上の配列を初期化 */
static RC
clean_node (MATRIX_DB_NODE *node)
{
	int ii1;

	RC_NULL_CHK(node);

	for(ii1=node->size+1; ii1<=MATRIX_DB_NODE_SIZE; ii1++){
		node->file[ii1]   = 0;
		node->offset[ii1] = 0;
	}
	/*
	if(node->size == 0){
		node->file[0]   = 0;
		node->offset[0] = 0;
	}
	*/

	for(ii1=node->size; ii1<MATRIX_DB_NODE_SIZE; ii1++){
		node->row[ii1] = 0;
		node->col[ii1] = 0;
	}

	return(NORMAL_RC);
}


/* size 以上の配列を初期化 */
static RC
clean_leaf (MATRIX_DB_LEAF *leaf)
{
	int ii1;

	RC_NULL_CHK(leaf);

	for(ii1=leaf->size; ii1<MATRIX_DB_LEAF_SIZE; ii1++){
		leaf->row[ii1] = 0;
		leaf->col[ii1] = 0;
		init_matrix33(leaf->v[ii1]);
	}

	return(NORMAL_RC);
}


/* B-Treeの時の search_matrix_db() とは違うので注意!! */
static RC
search_node (MATRIX_DB_NODE *node, int row, int col, int *index,
             int *found_flag, int *file, WRAP64_INT *offset)
{
	int comp;
	int low;
	int high;
	int middle;

	RC_NULL_CHK(node);
	RC_NULL_CHK(found_flag);
	RC_NULL_CHK(index);
	RC_NULL_CHK(file);
	RC_NULL_CHK(offset);
	if( (row < 0) || (col < 0) ) return(ARG_ERROR_RC);

	low = 0;
	high = node->size-1;

	*found_flag = 0;
	*index = -1;
	*file = 0;
	*offset = 0;

	if(node->size < 1){
		*index = 0; /* 必須!! */
		if((node->key_row == row) &&(node->key_col == col) )
			*found_flag = 1;
		if(node->file[0] != 0){
			*file   = node->file[0];
			*offset = node->offset[0];
		}
		return(NORMAL_RC);
	}

	while(low <= high){
		middle = (low + high)/2;
		comp = compar_row_col(node->row[middle], node->col[middle], row, col);
		if(comp == 0){ /* key1 = key2 */
			*found_flag = 1;
			*index = middle;
			*file = node->file[middle+1];
			*offset = node->offset[middle+1];
			break;
		}else if(comp > 0){ /* key1 > key2 */
			high = middle - 1;
			if(high <= low){
				if(high < 0) high = 0;
				comp = compar_row_col(node->row[high],
				                      node->col[high], row, col);
				if(comp == 0){ /* key1 = key2 */
					*found_flag = 1;
					*index = high;
					*file = node->file[*index+1];
					*offset = node->offset[*index+1];
				}else if(comp > 0){ /* key1 > key2 */
					*index = high;
					*file = node->file[*index];
					*offset = node->offset[*index];
				}else{ /* (comp < 0) */ /* key1 < key2 */
					if(node->size < high+1){
						return(UNKNOWN_ERROR_RC);
						/* *index = high; */
					}else{
						*index = high+1;
						*file = node->file[*index];
						*offset = node->offset[*index];
					}
				}
				break;
			}
		}else{ /* (comp < 0) */ /* key1 < key2 */
			low = middle + 1;
			if(high <= low){
				comp = compar_row_col(node->row[low], node->col[low],
				                      row, col);
				if(comp == 0){ /* key1 = key2 */
					*found_flag = 1;
					*index = low;
					*file = node->file[*index+1];
					*offset = node->offset[*index+1];
				}else if(comp > 0){ /* key1 > key2 */
					*index = low;
					*file = node->file[*index];
					*offset = node->offset[*index];
				}else{ /* (comp < 0) */ /* key1 < key2 */
					if( (node->size == MATRIX_DB_NODE_SIZE)
					 && (node->size <= low+1) ){
						*index = low;
						*file = node->file[*index+1];
						*offset = node->offset[*index+1];
					}else if(node->size < low+1){
						*index = low;
						*file = node->file[*index];
						*offset = node->offset[*index];
					}else{
						*index = low+1;
						*file = node->file[*index];
						*offset = node->offset[*index];
					}
				}
				break;
			}
		}
	}

	return(NORMAL_RC);
}


static RC
search_leaf (MATRIX_DB_LEAF *leaf, int row, int col, int *index,
             int *found_flag)
{
	int comp;
	int low;
	int high;
	int middle;

	RC_NULL_CHK(leaf);
	RC_NULL_CHK(found_flag);
	RC_NULL_CHK(index);
	if( (row < 0) || (col < 0) ) return(ARG_ERROR_RC);

	low = 0;
	high = leaf->size-1;

	*found_flag = 0;
	*index = -1;

	if(leaf->size < 1){
		*index = 0; /* 必須!! */
		return(NORMAL_RC);
	}

	while(low <= high){
		middle = (low + high)/2;
		comp = compar_row_col(leaf->row[middle], leaf->col[middle], row, col);
		if(comp == 0){
			*found_flag = 1;
			*index = middle;
			break;
		}else if(comp > 0){
			high = middle - 1;
			if(high <= low){
				if(high < 0) high = 0;
				comp = compar_row_col(leaf->row[high],
				                      leaf->col[high], row, col);
				if(comp == 0){
					*found_flag = 1;
					*index = high;
					break;
				}else if(comp > 0){
					*index = high;
				}else{
					if(leaf->size < high+1){
						return(UNKNOWN_ERROR_RC);
						/* *index = high; */
					}else{
						*index = high+1;
					}
				}
				break;
			}
		}else{ /* (comp < 0) */
			low = middle + 1;
			if(high <= low){
				comp = compar_row_col(leaf->row[low], leaf->col[low],
				                      row, col);
				if(comp == 0){
					*found_flag = 1;
					*index = high;
					break;
				}else if(comp > 0){
					*index = low;
				}else{ /* (comp < 0) */
					if(  (leaf->size == MATRIX_DB_LEAF_SIZE)
					  && (leaf->size <= low+1) ){
						*index = low;
						break;
					}else if(leaf->size < low+1){
						*index = low;
					}else{
						*index = low+1;
					}
				}
				break;
			}
		}
	}

	return(NORMAL_RC);
}


static int
is_root_node (MATRIX_DB_NODE *node)
{
	return(  (node->type == MATRIX_DB_ROOT_EDGE)
	       ||(node->type == MATRIX_DB_ROOT) );
}


static int
is_edge_node (MATRIX_DB_NODE *node)
{
	return(  (node->type == MATRIX_DB_ROOT_EDGE)
	       ||(node->type == MATRIX_DB_EDGE) );
}

static int
is_parent_edge (MATRIX_DB_LEAF *leaf, MATRIX_DB_NODE *edge_node)
{
	int ii1;

	if(leaf == NULL) return(0);
	if(edge_node == NULL) return(0);
	for(ii1=0; ii1<edge_node->size+1; ii1++){
		if( (leaf->file_self   == edge_node->file[ii1])
		 && (leaf->offset_self == edge_node->offset[ii1]) ){
			return(1);
		}
	}
	return(0);
}


static int
is_parent_node (MATRIX_DB_NODE *node, MATRIX_DB_NODE *parent)
{
	int ii1;

	if(node == NULL) return(0);
	if(parent == NULL) return(0);
	for(ii1=0; ii1<parent->size+1; ii1++){
		if( (node->file_self   == parent->file[ii1])
		 && (node->offset_self == parent->offset[ii1]) ){
			return(1);
		}
	}
	return(0);
}


#if 0
static int
is_range_node (MATRIX_DB_NODE *node, int index, int row, int col_min,
               int col_max)
{
	if(node == NULL) return(0);
	if( (index < 0) || (index >= node->size) ) return(0);
	if( (row < 0) || (col_min < 0) ) return(0);
	if(col_min > col_max) return(0);

	if( (node->row[index] == row) && (node->col[index] >= col_min)
	                              && (node->col[index] <  col_max) ){
		return(1);
	}

	return(0);
}
#endif


static int
is_range_leaf (MATRIX_DB_LEAF *leaf, int index, int row, int col_min,
               int col_max)
{
	if(leaf == NULL) return(0);
	if( (index < 0) || (index >= leaf->size) ) return(0);
	if( (row < 0) || (col_min < 0) ) return(0);
	if(col_min > col_max) return(0);

	if( (leaf->row[index] == row) && (leaf->col[index] >= col_min)
	                              && (leaf->col[index] <  col_max) ){
		return(1);
	}

	return(0);
}


/*
 * rowでsort して colでsort
 *  -1: key1 < key2
 *   0: key1 = key2
 *   1: key1 > key2
 */
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


/* 出力したuse_countは1足されることに注意 */
RC
print_matrix_db_tree (FILE *fp, MATRIX_DB_INFO *mdbi)
{
	int ii1, ii2, ii3;
	int edge_flag;;
	int size;
	int *file;
	WRAP64_INT *offset;
	int sub_size; /* subsrratum */
	int sub_alloc_size;
	int *sub_file = NULL;
	WRAP64_INT *sub_offset = NULL;

	RC_NULL_CHK(mdbi);

	/* RootNode を階層1としてセット */
	size = 1;
	RC_NULL_CHK( file = (int *)mm_alloc(sizeof(int)*(size)) );
	RC_NULL_CHK( offset = (WRAP64_INT *)mm_alloc(sizeof(WRAP64_INT)*(size)) );
	file[0]   = mdbi->root->file_self;
	offset[0] = mdbi->root->offset_self;
	sub_size = mdbi->root->size+1;
	sub_alloc_size = 0;
	if( is_edge_node(mdbi->root) ){
		edge_flag = 1;
	}else{
		edge_flag = 0;
	}

	for(ii1=1;; ii1++){
		sub_alloc_size += (sub_size*(MATRIX_DB_NODE_SIZE/2+1));
		sub_file = (int *)mm_alloc(sizeof(int)*sub_alloc_size);
		RC_NULL_CHK( sub_file );
		sub_offset = (WRAP64_INT *)mm_alloc(sizeof(WRAP64_INT)*sub_alloc_size);
		RC_NULL_CHK( sub_offset );
		for(ii2=0; ii2<sub_alloc_size; ii2++){
			sub_file[ii2] = 0;
			sub_offset[ii2] = 0;
		}

		fprintf(fp, "----- stratum %d (node) -----\n", ii1);
		for(ii2=0, sub_size=0; ii2<size; ii2++){
			MATRIX_DB_NODE *node;
			if((file[ii2] < 1) || (offset[ii2] < 0)) return(UNKNOWN_ERROR_RC);
			RC_TRY( read_node_start(mdbi, file[ii2], offset[ii2], &node) );
			RC_TRY( print_matrix_db_node(fp, node) );
			for(ii3=0; ii3<node->size+1; ii3++){
				sub_file[sub_size]   = node->file[ii3];
				sub_offset[sub_size] = node->offset[ii3];
				sub_size++;
			}
			if((sub_size+MATRIX_DB_NODE_SIZE+1) >= sub_alloc_size){
				sub_alloc_size += (MATRIX_DB_NODE_SIZE*2);
				RC_NULL_CHK( sub_file );
				sub_file = (int *)mm_realloc((void *)sub_file,
				            (size_t)(sizeof(int)*sub_alloc_size));
				RC_NULL_CHK( sub_file );
				sub_offset = (WRAP64_INT *)mm_realloc((void *)sub_offset,
				              (size_t)sizeof(WRAP64_INT)*sub_alloc_size);
				RC_NULL_CHK( sub_offset );
				for(ii3=sub_alloc_size-(MATRIX_DB_NODE_SIZE*2);
				    ii3<sub_alloc_size; ii3++){
					sub_file[ii3] = 0;
					sub_offset[ii3] = 0;
				}
			}
			if( is_edge_node(node) ){
				edge_flag = 1;
			}
			RC_TRY( read_node_end(node) );
		}

		/* swap */
		size = sub_size;
		sub_size = 0;
		RC_TRY( mm_free(file) );
		RC_TRY( mm_free(offset) );
		if(edge_flag){
			break;
		}else{
			file   = sub_file;
			offset = sub_offset;
			sub_file = NULL;
			sub_offset = NULL;
		}
	}

	log_printf(5, "[%d] size = %d\n", __LINE__, size);
	fprintf(fp, "----- stratum %d (leaf) -----\n", ++ii1);
	for(ii1=0; ii1<size; ii1++){
		MATRIX_DB_LEAF *leaf;
		log_printf(6, "sub_file = %d\n", sub_file[ii1]);
		log_printf(6, "sub_offset = %lld\n", sub_offset[ii1]);
		RC_TRY( read_leaf_start(mdbi, sub_file[ii1], sub_offset[ii1], &leaf) );
		RC_TRY( print_matrix_db_leaf(fp, leaf) );
		RC_TRY( read_leaf_end(leaf) );
	}
	RC_TRY( mm_free(sub_file) );
	RC_TRY( mm_free(sub_offset) );

	return(NORMAL_RC);
}


RC
print_matrix_db_cache_node_info (FILE *fp, MATRIX_DB_CACHE_NODE *cache)
{
	int ii1, ii2;

	RC_NULL_CHK(cache);

	for(ii1=0; ii1<cache->divisor; ii1++){
		fprintf(fp, "divisor : %2d\n", ii1);
		fprintf(fp, "asso dirty use     access file    offset type\n");
		for(ii2=0; ii2<cache->associate; ii2++){
			MATRIX_DB_NODE *node = &cache->cache[ii1][ii2];
			if(node->use_count == -1) continue;
#ifdef WIN32
			fprintf(fp, "[%2d]     %d  %2d %10d    %d %9I64d",
			             ii2,
			             node->dirty_flag, node->use_count,
			             cache->access_count[ii1][ii2],
			             node->file_self, node->offset_self);
#else /* WIN32 */
  #ifdef _LARGEFILE_SOURCE
			fprintf(fp, "[%2d]     %d  %2d %10d    %d %9lld",
			             ii2,
			             node->dirty_flag, node->use_count,
			             cache->access_count[ii1][ii2],
			             node->file_self, node->offset_self);
  #else
			fprintf(fp, "[%2d]     %d  %2d %10d    %d %9ld",
			             ii2,
			             node->dirty_flag, node->use_count,
			             cache->access_count[ii1][ii2],
			             node->file_self, node->offset_self);
  #endif
#endif  /* WIN32 */
			if(node->type == MATRIX_DB_INIT){
				fprintf(fp, " init\n");
			}else if(node->type == MATRIX_DB_ROOT){
				fprintf(fp, " root\n");
			}else if(node->type == MATRIX_DB_ROOT_EDGE){
				fprintf(fp, " root(edge)\n");
			}else if(node->type == MATRIX_DB_GRID){
				fprintf(fp, " grid\n");
			}else if(node->type == MATRIX_DB_EDGE){
				fprintf(fp, " edge\n");
			}else{
				DEBUG_PRINT;
				return(UNKNOWN_ERROR_RC);
			}
		}
	}

	return(NORMAL_RC);
}


RC
print_matrix_db_cache_leaf_info (FILE *fp, MATRIX_DB_CACHE_LEAF *cache)
{
	int ii1, ii2;

	RC_NULL_CHK(cache);

	for(ii1=0; ii1<cache->divisor; ii1++){
		fprintf(fp, "divisor : %2d\n", ii1);
		fprintf(fp, "asso dirty use     access file    offset\n");
		for(ii2=0; ii2<cache->associate; ii2++){
			MATRIX_DB_LEAF *leaf = &cache->cache[ii1][ii2];
			if(leaf->use_count == -1) continue;
#ifdef WIN32
			fprintf(fp, "[%2d]     %d  %2d %10d    %d %9I64d\n",
			             ii2,
			             leaf->dirty_flag, leaf->use_count,
			             cache->access_count[ii1][ii2],
			             leaf->file_self, leaf->offset_self);
#else /* WIN32 */
  #ifdef _LARGEFILE_SOURCE
			fprintf(fp, "[%2d]     %d  %2d %10d    %d %9lld\n",
			             ii2,
			             leaf->dirty_flag, leaf->use_count,
			             cache->access_count[ii1][ii2],
			             leaf->file_self, leaf->offset_self);
  #else
			fprintf(fp, "[%2d]     %d  %2d %10d    %d %9ld\n",
			             ii2,
			             leaf->dirty_flag, leaf->use_count,
			             cache->access_count[ii1][ii2],
			             leaf->file_self, leaf->offset_self);
  #endif
#endif /* WIN32 */
		}
	}

	return(NORMAL_RC);
}


RC
print_matrix_db_node (FILE *fp, MATRIX_DB_NODE *node)
{
	int ii1, ii2=0;
	int node_size;
	int size;

	RC_NULL_CHK(node);

	RC_TRY( print_matrix_db_node_info(fp, node) );
	
#if 1
	node_size = node->size+1;
	if(node_size%2){
		size = node_size/2+1;
	}else{
		size = node_size/2;
	}
	if(node->size == MATRIX_DB_NODE_SIZE) size = MATRIX_DB_NODE_SIZE/2+1;
	for(ii1=0; ii1<size; ii1++){
		ii2 = ii1+size;
#ifdef WIN32
		fprintf(fp, "file[%3d] = %d, offset[%3d] = %7I64d  "
		            "file[%3d] = %d, offset[%3d] = %7I64d\n",
		             ii1, node->file[ii1], ii1, node->offset[ii1],
		             ii2, node->file[ii2], ii2, node->offset[ii2]);
#else /* WIN32 */
  #ifdef _LARGEFILE_SOURCE
		fprintf(fp, "file[%3d] = %d, offset[%3d] = %7lld  "
		            "file[%3d] = %d, offset[%3d] = %7lld\n",
		             ii1, node->file[ii1], ii1, node->offset[ii1],
		             ii2, node->file[ii2], ii2, node->offset[ii2]);
  #else
		fprintf(fp, "file[%3d] = %d, offset[%3d] = %7ld  "
		            "file[%3d] = %d, offset[%3d] = %7ld\n",
		             ii1, node->file[ii1], ii1, node->offset[ii1],
		             ii2, node->file[ii2], ii2, node->offset[ii2]);
  #endif
#endif  /* WIN32 */
	}
#endif

	if(node->size == 0) return(NORMAL_RC);

	if(node->size%2){
		size = node->size/2+1;
	}else{
		size = node->size/2;
	}
	if(node->size == MATRIX_DB_NODE_SIZE) size = MATRIX_DB_NODE_SIZE/2;
	for(ii1=0; ii1<size; ii1++){
		ii2 = ii1+size;
		fprintf(fp, "row[%3d] = %4d, col[%3d] = %4d   "
		            "row[%3d] = %4d, col[%3d] = %4d\n",
		        ii1, node->row[ii1], ii1, node->col[ii1],
		        ii2, node->row[ii2], ii2, node->col[ii2]);
	}
	if(node->size == MATRIX_DB_NODE_SIZE){
		int ii2 = MATRIX_DB_NODE_SIZE-1;
		fprintf(fp, "                                   "
		            "row[%3d] = %4d, col[%3d] = %4d\n",
		        ii2, node->row[ii2], ii2, node->col[ii2]);
	}

	return(NORMAL_RC);
}


RC
print_matrix_db_node_info (FILE *fp, MATRIX_DB_NODE *node)
{
	RC_NULL_CHK(node);

	fprintf(fp, "size = %d\n", node->size);
	fprintf(fp, "dirty_flag = %d\n", node->dirty_flag);
	fprintf(fp, "delete_flag = %d\n", node->delete_flag);
	fprintf(fp, "use_count = %d\n", node->use_count);
	fprintf(fp, "key_row = %d\n", node->key_row);
	fprintf(fp, "key_col = %d\n", node->key_col);
	if(node->type == MATRIX_DB_INIT){
		fprintf(fp, "type = init\n");
	}else if(node->type == MATRIX_DB_ROOT){
		fprintf(fp, "type = root\n");
	}else if(node->type == MATRIX_DB_ROOT_EDGE){
		fprintf(fp, "type = root(edge)\n");
	}else if(node->type == MATRIX_DB_GRID){
		fprintf(fp, "type = grid\n");
	}else if(node->type == MATRIX_DB_EDGE){
		fprintf(fp, "type = edge\n");
	}else{
		fprintf(fp, "type = %d\n", node->type);
		return(UNKNOWN_ERROR_RC);
	}

#ifdef WIN32
	fprintf(fp, "file_self = %d\n",      node->file_self);
	fprintf(fp, "offset_self = %I64d\n", node->offset_self);
#else /* WIN32 */
  #ifdef _LARGEFILE_SOURCE
	fprintf(fp, "file_self = %d\n",     node->file_self);
	fprintf(fp, "offset_self = %lld\n", node->offset_self);
  #else /* _LARGEFILE_SOURCE */
	fprintf(fp, "file_self = %d\n",    node->file_self);
	fprintf(fp, "offset_self = %ld\n", node->offset_self);
  #endif /* _LARGEFILE_SOURCE */
#endif  /* WIN32 */

	return(NORMAL_RC);
}


RC
print_matrix_db_leaf (FILE *fp, MATRIX_DB_LEAF *leaf)
{
	int ii1, ii2=0;
	int size;

	RC_NULL_CHK(leaf);

	RC_TRY( print_matrix_db_leaf_info(fp, leaf) );
	if(leaf->size == 0) return(NORMAL_RC);
	
	if(leaf->size%2){
		size = leaf->size/2+1;
	}else{
		size = leaf->size/2;
	}
	if(leaf->size == MATRIX_DB_LEAF_SIZE) size = MATRIX_DB_LEAF_SIZE/2;
	for(ii1=0; ii1<size; ii1++){
		ii2 = ii1+size;
		fprintf(fp, "row[%3d] = %4d, col[%3d] = %4d   "
		            "row[%3d] = %4d, col[%3d] = %4d\n",
		        ii1, leaf->row[ii1], ii1, leaf->col[ii1],
		        ii2, leaf->row[ii2], ii2, leaf->col[ii2]);
	}
	if((leaf->size == MATRIX_DB_LEAF_SIZE)&&(MATRIX_DB_LEAF_SIZE%2)){
		int ii2 = MATRIX_DB_LEAF_SIZE-1;
		fprintf(fp, "                                   "
		            "row[%3d] = %4d, col[%3d] = %4d\n",
		        ii2, leaf->row[ii2], ii2, leaf->col[ii2]);
	}

#if 1
	for(ii1=0; ii1<leaf->size; ii1++){
		fprintf(fp, "v[%2d]\n", ii1);
		print_matrix33 (fp, leaf->v[ii1]);
	}
#endif

	return(NORMAL_RC);
}


RC
print_matrix_db_leaf_info (FILE *fp, MATRIX_DB_LEAF *leaf)
{
	RC_NULL_CHK(leaf);

	fprintf(fp, "size = %d\n", leaf->size);
	fprintf(fp, "dirty_flag = %d\n", leaf->dirty_flag);
	fprintf(fp, "delete_flag = %d\n", leaf->delete_flag);
	fprintf(fp, "use_count = %d\n", leaf->use_count);
	fprintf(fp, "key_row = %d\n", leaf->key_row);
	fprintf(fp, "key_col = %d\n", leaf->key_col);

#ifdef WIN32
	fprintf(fp, "file_self = %d\n",      leaf->file_self);
	fprintf(fp, "offset_self = %I64d\n", leaf->offset_self);
	fprintf(fp, "deleted_offset = %I64d\n", leaf->deleted_offset);
#else /* WIN32 */
  #ifdef _LARGEFILE_SOURCE
	fprintf(fp, "file_self = %d\n",     leaf->file_self);
	fprintf(fp, "offset_self = %lld\n", leaf->offset_self);
	fprintf(fp, "deleted_offset = %lld\n", leaf->deleted_offset);
  #else /* _LARGEFILE_SOURCE */
	fprintf(fp, "file_self = %d\n",    leaf->file_self);
	fprintf(fp, "offset_self = %ld\n", leaf->offset_self);
	fprintf(fp, "deleted_offset = %ld\n", leaf->deleted_offset);
  #endif /* _LARGEFILE_SOURCE */
#endif  /* WIN32 */

	return(NORMAL_RC);
}


void
print_matrix_db_struct_size (FILE *fp)
{
	fprintf(fp, " MATRIX_DB_NODE_SIZE    = %d\n", MATRIX_DB_NODE_SIZE);
	fprintf(fp, " MATRIX_DB_NODE_HEADER  = %d\n", MATRIX_DB_NODE_HEADER);
#ifdef WIN32
	fprintf(fp, " MATRIX_DB_NODE         = %d [Byte]\n",
	                                       sizeof(MATRIX_DB_NODE) );
#else /* WIN32 */
  #ifdef _LARGEFILE_SOURCE
	fprintf(fp, " MATRIX_DB_NODE         = %d [Byte]\n",
	                                       sizeof(MATRIX_DB_NODE) );
  #else /* _LARGEFILE_SOURCE */
	fprintf(fp, " MATRIX_DB_NODE         = %ld [Byte]\n",
	                                       sizeof(MATRIX_DB_NODE) );
  #endif /* _LARGEFILE_SOURCE */
#endif  /* WIN32 */
	fprintf(fp, " MATRIX_DB_LEAF_SIZE    = %d\n", MATRIX_DB_LEAF_SIZE);
	fprintf(fp, " MATRIX_DB_LEAF_HEADER  = %d\n", MATRIX_DB_LEAF_HEADER);
#ifdef WIN32
	fprintf(fp, " MATRIX_DB_LEAF         = %d [Byte]\n",
	                                       sizeof(MATRIX_DB_LEAF) );
#else /* WIN32 */
  #ifdef _LARGEFILE_SOURCE
	fprintf(fp, " MATRIX_DB_LEAF         = %d [Byte]\n",
	                                       sizeof(MATRIX_DB_LEAF) );
  #else /* _LARGEFILE_SOURCE */
	fprintf(fp, " MATRIX_DB_LEAF         = %ld [Byte]\n",
	                                       sizeof(MATRIX_DB_LEAF) );
  #endif /* _LARGEFILE_SOURCE */
#endif  /* WIN32 */
}


