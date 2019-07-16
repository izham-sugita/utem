/*********************************************************************
 * mtrix_db_plus.h
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

#ifndef MATRIX_DB_PLUS_H
#define MATRIX_DB_PLUS_H

#include <stdio.h>
#include "wrapper_64.h"
#include "scratch_io.h"
#include "rc.h"

typedef enum {
	MATRIX_DB_INIT,
	MATRIX_DB_ROOT,
	MATRIX_DB_ROOT_EDGE, /* 高さ1の木 */
	MATRIX_DB_GRID,
	MATRIX_DB_EDGE, /* index部の底辺ノード */
} MATRIX_DB_NODE_TYPE;

/* B+-Treeのindex(node)部 */
/* MATRIX_DB_NODE_SIZE は奇数であること */
#define MATRIX_DB_NODE_SIZE (5)
#define MATRIX_DB_NODE_HEADER (0)
typedef struct{
	int size;
	int dirty_flag; /* スクラッチファイル上と{一致:0, 不一致:1} */
	int delete_flag; /* 消去{可: 1, 不可:0} */
	int use_count;  /* <0:未使用, =0:参照なし, >0:使用中 */
	int key_row;
	int key_col;
	MATRIX_DB_NODE_TYPE type; /* root, node, edge */
	int header[MATRIX_DB_NODE_HEADER];
	int file_self;          /* 自分の書き込みファイル */
	int file[MATRIX_DB_NODE_SIZE+1];
	WRAP64_INT offset_self; /* 自分の書き込み位置 */
	WRAP64_INT offset[MATRIX_DB_NODE_SIZE+1];
	int row[MATRIX_DB_NODE_SIZE];
	int col[MATRIX_DB_NODE_SIZE];
} MATRIX_DB_NODE; 

/* B+-Treeのデータ(Leaf)部 */
#define MATRIX_DB_LEAF_SIZE (6)
#define MATRIX_DB_LEAF_HEADER (0)
typedef struct{
	int size;
	int dirty_flag; /* スクラッチファイル上と{一致:0, 不一致:1} */
	int delete_flag; /* 消去{可: 1, 不可:0} */
	int use_count;  /* <0:未使用, =0:参照なし, >0:使用中 */
	int key_row;
	int key_col;
	int header[MATRIX_DB_LEAF_HEADER];
	int file_self;            /* 自分の書き込みファイル */
	WRAP64_INT offset_self;   /* 自分の書き込み位置 */
	WRAP64_INT deleted_offset;/* 一つ前に削除されたleafの書き込み位置 */
	int row[MATRIX_DB_LEAF_SIZE];
	int col[MATRIX_DB_LEAF_SIZE];
	double v[MATRIX_DB_LEAF_SIZE][3][3];
} MATRIX_DB_LEAF;

#define MATRIX_DB_SF_SIZE (1) /* スクラッチファイル数 */
typedef struct{
	SCRATCH_FILE *scratch_array[MATRIX_DB_SF_SIZE];
	WRAP64_INT last_offset[MATRIX_DB_SF_SIZE]; /* scratch fileの末尾offset */
	WRAP64_INT blank_offset[MATRIX_DB_SF_SIZE];/* 削除された場所offset */
	WRAP64_INT current_offset[MATRIX_DB_SF_SIZE];/* 現在の読み書き位置 */
} MATRIX_DB_SCRATCH;

typedef struct{
	int associate;      /* キャッシュの連結数 */
	int divisor;        /* = (cache_size/associate) */
	MATRIX_DB_NODE **cache; /* cache配列 */
	int **access_count; /* cacheのアクセスカウンタ */
} MATRIX_DB_CACHE_NODE;

typedef struct{
	int associate;      /* キャッシュの連結数 */
	int divisor;        /* = (cache_size/associate) */
	MATRIX_DB_LEAF **cache; /* cache配列 */
	int **access_count; /* cacheのアクセスカウンタ */
} MATRIX_DB_CACHE_LEAF;


#define MATRIX_DB_ROUTE_SIZE (16) /* 経路記憶配列のサイズ */
typedef struct{
	MATRIX_DB_SCRATCH sf_node;
	MATRIX_DB_SCRATCH sf_leaf;
	MATRIX_DB_CACHE_NODE cache_node;
	MATRIX_DB_CACHE_LEAF cache_leaf;
	MATRIX_DB_NODE *root; /* RootNode へのポインタ */
	int route_size;  /* 通過経路数 */
	int route_count; /* 無限ループ回避用 */
	int *route_file; /* 経路情報 */
	WRAP64_INT *route_offset; /* 経路情報 */
} MATRIX_DB_INFO;

RC open_matrix_db(MATRIX_DB_INFO *mdbi, const char scratch_dir[],
                  WRAP64_INT index_cache_size, WRAP64_INT data_cache_size,
                  int index_associate, int data_assciate);
RC close_matrix_db(MATRIX_DB_INFO *mdbi);
RC add_v_matrix_db(int row, int col, double v[][3], MATRIX_DB_INFO *mdbi);
RC get_v_matrix_db(int row, int col, double v[][3], int *found_flag,
                   MATRIX_DB_INFO *mdbi);
RC delete_v_matrix_db(int row, int col, MATRIX_DB_INFO *mdbi);

RC min_col_matrix_db(MATRIX_DB_INFO *mdbi, int row, int col, int *min_col);
RC get_row_matrix_db(MATRIX_DB_INFO *mdbi, int row, int col_min, int col_max,
                     int *size, double array[][3][3], int col_array[]);
RC add_row_matrix_db(MATRIX_DB_INFO *mdbi, int row, int size,
                     double array[][3][3], int col_array[]);
#if 0
RC delete_row_matrix_db(MATRIX_DB_INFO *mdbi, int row, int col_min,int col_max);
RC resize_cache_matrix_db(MATRIX_DB_INFO *mdbi, WRAP64_INT size);
#endif

/* debug用 */
RC print_matrix_db_tree(FILE *fp, MATRIX_DB_INFO *mdbi);
RC print_matrix_db_cache_node_info(FILE *fp, MATRIX_DB_CACHE_NODE *cache);
RC print_matrix_db_cache_leaf_info(FILE *fp, MATRIX_DB_CACHE_LEAF *cache);
RC print_matrix_db_node(FILE *fp, MATRIX_DB_NODE *node);
RC print_matrix_db_node_info(FILE *fp, MATRIX_DB_NODE *node);
RC print_matrix_db_leaf(FILE *fp, MATRIX_DB_LEAF *leaf);
RC print_matrix_db_leaf_info(FILE *fp, MATRIX_DB_LEAF *leaf);
void print_matrix_db_struct_size(FILE *fp);

#endif /* MATRIX_DB_PLUS_H */

