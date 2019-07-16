/*********************************************************************
 * vect_storage.h
 *
 * Copyright (C) 2005 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Thanks to following contributors.
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: vect_storage.h 505 2005-10-07 08:06:43Z sasaoka $ */

#ifndef VECT_STORAGE_H
#define VECT_STORAGE_H

#include <stdio.h>
#include "rc.h"
#include "scratch_io.h"
#include "nonzero_cg.h"


#define VECT_STORAGE_TMP_SIZE   (5)
#define VECT_STORAGE_BUF_SIZE   (8)
typedef struct{
	int array_size;
	int vect_size;
	SCRATCH_FILE *sf;
	int sf_handle[SCRATCH_DATA_SIZE];
	double *tmp[VECT_STORAGE_TMP_SIZE];
	double *buf[VECT_STORAGE_BUF_SIZE];
	int current_sf_index;  /* -1: buf[][]をテンポラリとして使用中 */
	int dirty_flag;   /* buf[][] の内容が 0: スクラッチファイル上と一致   */
	                  /*                  1: スクラッチファイル上と不一致 */
} VECT_STORAGE;


VECT_STORAGE *open_vect_storage(int array_size, int vect_size,
                                const char scratch_dir[]);
RC close_vect_storage(VECT_STORAGE *vs);
double *tmp_ptr_vect_storage(VECT_STORAGE *vs, int tmp_index);
double **buf_ptr_vect_storage(VECT_STORAGE *vs);
RC direct_write_vect_storage(VECT_STORAGE *vs, const double vect[], int index);
RC direct_read_vect_storage(VECT_STORAGE *vs, double vect[], int index);
RC buffered_write_vect_storage(VECT_STORAGE *vs,
                               const double vect[], int index);
RC buffered_read_vect_storage(VECT_STORAGE *vs, double vect[], int index);
RC get_buf_vect_storage(VECT_STORAGE *vs, int index);
double **start_use_buf4tmp_vect_storage(VECT_STORAGE *vs);
RC stop_use_buf4tmp_vect_storage(VECT_STORAGE *vs, int w_index, int r_index);
RC compaction_vect_storage(VECT_STORAGE *vs);
RC mul_nonzero3s_AtBA_vect_storage(VECT_STORAGE *A,
                                   NONZERO_MATRIX3 B, double **C);
RC B_orthogonalize_vect_storage(VECT_STORAGE *A, NONZERO_MATRIX3 B);

#endif /* VECT_STORAGE_H */

