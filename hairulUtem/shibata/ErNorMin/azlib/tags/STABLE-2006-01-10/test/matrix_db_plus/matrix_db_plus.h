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
	MATRIX_DB_NODE,
	MATRIX_DB_EDGE
} MATRIX_DB_NODE_TYPE;

/* B+-Treeのインデックス部 */
#define MATRIX_DB_INDEX_SIZE (5)
#define MATRIX_DB_INDEX_HEADER (0)
typedef struct{
	int size;
	int dirty_flag; /* スクラッチファイル上と{一致:0, 不一致:1} */
	int use_count;  /* <0:未使用, =0:参照なし, >0:使用中 */
	int key_row;
	int key_col;
	MATRIX_DB_NODE_TYPE type; /* root, node, edge */
	int header[MATRIX_DB_INDEX_HEADER];
	int file_self;          /* 自分の書き込みファイル */
	int file[MATRIX_DB_INDEX_SIZE+1];
	WRAP64_INT offset_self; /* 自分の書き込み位置 */
	WRAP64_INT offset[MATRIX_DB_INDEX_SIZE+1];
	int row[MATRIX_DB_INDEX_SIZE];
	int col[MATRIX_DB_INDEX_SIZE];
} MATRIX_DB_INDEX; 

/* B+-Treeのデータ(Leaf)部 */
#define MATRIX_DB_DATA_SIZE (10)
#define MATRIX_DB_DATA_HEADER (0)
typedef struct{
	int size;
	int dirty_flag; /* スクラッチファイル上と{一致:0, 不一致:1} */
	int use_count;  /* <0:未使用, =0:参照なし, >0:使用中 */
	int header[MATRIX_DB_DATA_HEADER];
	int file_self;            /* 自分の書き込みファイル */
	WRAP64_INT offset_self;   /* 自分の書き込み位置 */
	int row[MATRIX_DB_DATA_SIZE];
	int col[MATRIX_DB_DATA_SIZE];
	double v[MATRIX_DB_DATA_SIZE][3][3];
} MATRIX_DB_DATA;

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
	void **cache;       /* cache配列 */
	int **access_count; /* cacheのアクセスカウンタ */
} MATRIX_DB_CACHE;

#define MATRIX_DB_ROUTE_SIZE (16) /* 経路記憶配列のサイズ */
typedef struct{
	MATRIX_DB_SCRATCH sf_index;
	MATRIX_DB_SCRATCH sf_data;
	MATRIX_DB_CACHE ca_index;
	MATRIX_DB_CACHE ca_data;
	MATRIX_DB_INDEX *root; /* RootNode へのポインタ */
	int route_size;  /* 通過経路数 */
	int route_count; /* 無限ループ回避用 */
	int *route_file; /* 経路情報 */
	WRAP64_INT *route_offset; /* 経路情報 */
} MATRIX_DB_INFO;

void print_matrix_db_size(FILE *fp);

#endif /* MATRIX_DB_PLUS_H */

