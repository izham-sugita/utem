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

/* $Id: matrix_db.h 721 2006-01-05 07:38:13Z sasaoka $ */

#ifndef MATRIX_DB_H
#define MATRIX_DB_H

#include <stdio.h>
#include "wrapper_64.h"
#include "scratch_io.h"
#include "rc.h"
#include "fem_struct.h"


#define MATRIX_DB_SIZE (177)
#define MATRIX_DB_HEADER_SIZE (11)
/*
 [16KB]
#define MATRIX_DB_SIZE (177)
#define MATRIX_DB_HEADER_SIZE (11)
 [8KB]
#define MATRIX_DB_SIZE (87)
#define MATRIX_DB_HEADER_SIZE (33)
 [テスト用]
#define MATRIX_DB_SIZE (9)
#define MATRIX_DB_HEADER_SIZE (1827)
*/

typedef struct{
	int size;       /* row, col, v のサイズ */
	int dirty_flag;   /* スクラッチファイル上と{一致:0, 不一致:1} */
	int use_count;    /* <0:未使用, =0:参照なし, >0:使用中 */
	int key_row;
	int key_col;
	int header[MATRIX_DB_HEADER_SIZE];
	int p_file_parent;          /* 親の書き込みファイル */
	int p_file_self;            /* 自分の書き込みファイル */
	int p_file[MATRIX_DB_SIZE+1];
	WRAP64_INT p_offset_parent; /* 親の書き込み位置 */
	WRAP64_INT p_offset_self;   /* 自分の書き込み位置 */
	WRAP64_INT p_offset[MATRIX_DB_SIZE+1];
	int row[MATRIX_DB_SIZE];
	int col[MATRIX_DB_SIZE];
	double v[MATRIX_DB_SIZE][3][3];
} MATRIX_DB;

#define MATRIX_DB_ROUTE_SIZE (16) /* 経路記憶配列のサイズ */
#define MATRIX_DB_ASSOCIATE (64) /* キャッシュの深さ */
#define MATRIX_DB_SF_SIZE  (1)  /* スクラッチファイル数 */
typedef struct{
	SCRATCH_FILE *scratch_array[MATRIX_DB_SF_SIZE];
	WRAP64_INT last_p_offset[MATRIX_DB_SF_SIZE]; /* scratch fileの末尾offset */
	WRAP64_INT blank_p_offset[MATRIX_DB_SF_SIZE];/* 削除された場所offset */
	WRAP64_INT current_p_offset[MATRIX_DB_SF_SIZE];/* 現在の読み書き位置 */
	int cache_divisor; /* (cache size) = MDB_ASSOCIATE_SIZE * cache_divisor */
	MATRIX_DB (*cache)[MATRIX_DB_ASSOCIATE];  /* cache配列 */
	int (*access_count)[MATRIX_DB_ASSOCIATE]; /* cacheのアクセスカウンタ */
	MATRIX_DB *root; /* RootNode へのポインタ */
	int route_size; /* 通過経路数 */
	int route_count; /* 無限ループ回避用 */
	int *route_file; /* 経路情報 */
	WRAP64_INT *route_offset; /* 経路情報 */
} MATRIX_DB_INFO;


RC open_matrix_db(const char scratch_dir[], WRAP64_INT size,
                  MATRIX_DB_INFO *mdbi);
RC close_matrix_db(MATRIX_DB_INFO *mdbi);
RC get_v_matrix_db(int row, int col, double v[][3], int *found_flag,
                   MATRIX_DB_INFO *mdbi);
RC put_v_matrix_db(int row, int col, double v[][3], MATRIX_DB_INFO *mdbi);
RC add_v_matrix_db(int row, int col, double v[][3], MATRIX_DB_INFO *mdbi);
RC mul_v_matrix_db(int row, int col, double v[][3], MATRIX_DB_INFO *mdbi);
RC delete_v_matrix_db(int row, int col, MATRIX_DB_INFO *mdbi);
int min_col_matrix_db(MATRIX_DB_INFO *mdbi, int row);
int min_col_sub_matrix_db(MATRIX_DB_INFO *mdbi, int row, int col);
RC add_row_matrix_db(MATRIX_DB_INFO *mdbi, int row, int size,
                     double array[][3][3], int col_array[]);
RC get_row_matrix_db(MATRIX_DB_INFO *mdbi, int row, int *size,
                     double array[][3][3], int col_array[]);
RC get_row_col_matrix_db(MATRIX_DB_INFO *mdbi, int row, int col, int *size,
                         double array[][3][3], int col_array[]);
RC get_row_range_matrix_db(MATRIX_DB_INFO *mdbi, int row, int col_min,
                           int col_max, int *size, double array[][3][3],
                           int col_array[]);
RC delete_row_matrix_db(MATRIX_DB_INFO *mdbi, int row);
RC delete_row_col_matrix_db(MATRIX_DB_INFO *mdbi, int row, int col_min);
RC delete_row_range_matrix_db(MATRIX_DB_INFO *mdbi, int row, int col_min,
                              int col_max);
RC resize_cache_matrix_db(MATRIX_DB_INFO *mdbi, WRAP64_INT size);


/* B-Treeに関する関数 */
RC search_v_btree_matrix_db(MATRIX_DB_INFO *mdbi, MATRIX_DB *mdb, int row,
                            int col, int *found_flag, int *index, int *p_file,
                            WRAP64_INT *p_offset);
RC split_btree_matrix_db(MATRIX_DB *child, MATRIX_DB *new_child,
                         MATRIX_DB_INFO *mdbi);
RC delete_v_btree_matrix_db(int row, int col, MATRIX_DB_INFO *mdbi);

/* node */
RC balance_node_matrix_db(MATRIX_DB *mdb, MATRIX_DB_INFO *mdbi);    
RC move_from_left_node_matrix_db(MATRIX_DB *left, MATRIX_DB *right, 
                                 MATRIX_DB *parent, MATRIX_DB_INFO *mdbi,
                                 int index);
RC move_from_right_node_matrix_db(MATRIX_DB *left, MATRIX_DB *right,
                                  MATRIX_DB *parent, MATRIX_DB_INFO *mdbi,
                                  int index);
RC marge_node_matrix_db(MATRIX_DB *left, MATRIX_DB *right, MATRIX_DB *parent,
                        MATRIX_DB_INFO *mdbi, int index);
/* leaf */
RC balance_leaf_matrix_db(MATRIX_DB *mdb, MATRIX_DB_INFO *mdbi);
RC move_from_left_leaf_matrix_db(MATRIX_DB *left, MATRIX_DB *right,
                                 MATRIX_DB *parent, int index);
RC move_from_right_leaf_matrix_db(MATRIX_DB *left, MATRIX_DB *right,
                                  MATRIX_DB *parent, int index);
RC marge_leaf_matrix_db(MATRIX_DB *left, MATRIX_DB *right, MATRIX_DB *parent,
                        MATRIX_DB_INFO *mdbi, int index);

/* matrix_dbに関する低レベル関数 */
RC insert_v_node_matrix_db(int row, int col, double v[][3], int p_file,
                           WRAP64_INT p_offset, int index, MATRIX_DB *mdb);
RC insert_v_leaf_matrix_db(int row, int col, double v[][3],
                           int index, MATRIX_DB *mdb);
RC delete_v_node_matrix_db(int index, MATRIX_DB *mdb, MATRIX_DB_INFO *mdbi);
RC delete_v_leaf_matrix_db(int index, MATRIX_DB *mdb, MATRIX_DB_INFO *mdbi);
RC search_matrix_db(int row, int col, MATRIX_DB *mdb, int *found_flag,
                    int *index, int *p_file, WRAP64_INT *p_offset);

RC flush_cache_matrix_db(MATRIX_DB_INFO *mdbi);
RC delete_matrix_db(MATRIX_DB_INFO *mdbi, int db_file, WRAP64_INT db_offset);
RC new_matrix_db_start(MATRIX_DB **mdb, MATRIX_DB_INFO *mdbi);
RC new_matrix_db_end(MATRIX_DB *mdb);
RC write_matrix_db_start(int p_file, WRAP64_INT p_offset, MATRIX_DB **mdb,
                         MATRIX_DB_INFO *mdbi);
RC write_matrix_db_end(MATRIX_DB *mdb);
RC read_matrix_db_start(int p_file, WRAP64_INT p_offset, MATRIX_DB **mdb,
                        MATRIX_DB_INFO *mdbi);
RC read_matrix_db_end(MATRIX_DB *mdb);
RC put_matrix_db(MATRIX_DB_INFO *mdbi, int p_file, WRAP64_INT p_offset,
                 MATRIX_DB *ptr);
RC get_matrix_db(MATRIX_DB_INFO *mdbi, int p_file, WRAP64_INT p_offset,
                 MATRIX_DB *ptr);

int is_init_matrix_db(MATRIX_DB *mdb);
int is_leaf_matrix_db(MATRIX_DB *mdb);
int is_root_matrix_db(MATRIX_DB *mdb);
int is_node_matrix_db(MATRIX_DB *mdb);
int is_parent_matrix_db(MATRIX_DB *mdb, MATRIX_DB *parent);
int is_range_matrix_db(MATRIX_DB *mdb, int index, int row, int col_min,
                       int col_max);

RC set_v_matrix_db(int row, int col, double v[][3], int index, MATRIX_DB *mdb);
RC set_root_matrix_db(MATRIX_DB *root, MATRIX_DB_INFO *mdbi);
RC clean_matrix_db(MATRIX_DB *mdb);
RC init_matrix_db(MATRIX_DB *mdb);

/* debug print */
RC print_matrix_db_cache_info(FILE *fp, MATRIX_DB_INFO *mdbi);
RC print_stratum_info_matrix_db(FILE *fp, MATRIX_DB_INFO *mdbi);
RC print_btree_matrix_db(FILE *fp, MATRIX_DB_INFO *mdbi);
RC print_matrix_db(FILE *fp, MATRIX_DB *mdb);
RC print_matrix_db_debug(FILE *fp, MATRIX_DB *mdb);
RC print_matrix_db_info(FILE *fp, MATRIX_DB *mdb);
void print_matrix_db_size(FILE *fp);


/* matrix_db_amls.c */
RC result_chk_amls(int block_num, const int block_size_array[], int row_size,
                   int evects_size_array[], MATRIX_DB_INFO *mdbi_K,
                   MATRIX_DB_INFO *mdbi_M, MATRIX_DB_INFO *mdbi_U,
                   const int max_col_U[]);
RC make_UtKU_matrix_db_amls(int *block_size_array, int block_num, int row_size,
                            MATRIX_DB_INFO *mdbi_K, MATRIX_DB_INFO *mdbi_U,
                            int max_col_U[]);
RC make_KM_matrix_db_amls(NODE_ARRAY node, ELEMENT_ARRAY element,
                          MATERIAL_PROP_ARRAY material,
                          PHYSICAL_PROP_ARRAY physical,
                          MATRIX_DB_INFO *mdbi_K, MATRIX_DB_INFO *mdbi_M);
RC cal_block_array_size_amls (int *block_size_array, int block_num,
                              int row_size, MATRIX_DB_INFO *mdbi_K);



#endif /* MATRIX_DB_H */

