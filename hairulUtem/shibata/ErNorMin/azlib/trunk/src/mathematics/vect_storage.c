/*********************************************************************
 * vect_storage.c
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

/* $Id: vect_storage.c 842 2006-12-04 06:21:20Z sasaoka $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rc.h"
#include "scratch_io.h"
#include "math_utl.h"
#include "memory_manager.h"
#include "nonzero_cg.h"
#include "vect_storage.h"


#define VS_ARRAY_SIZE  (8)
static VECT_STORAGE VSArray[VS_ARRAY_SIZE];
static int InitFlag = 0;

static RC flush_vect_storage(VECT_STORAGE *vs);


VECT_STORAGE *
open_vect_storage (int array_size, int vect_size, const char scratch_dir[])
{
	int ii1, ii2;
	VECT_STORAGE *vs;
	int alloc_size;

	if( (array_size <= 0)||(vect_size <= 0) ) return(NULL);
	if(array_size > SCRATCH_DATA_SIZE*VECT_STORAGE_BUF_SIZE) return(NULL);

	if(InitFlag == 0){
		for(ii1=0; ii1<VS_ARRAY_SIZE; ii1++){
			VSArray[ii1].array_size = 0;
			VSArray[ii1].vect_size = 0;
			VSArray[ii1].sf = NULL;
			for(ii2=0; ii2<SCRATCH_DATA_SIZE; ii2++){
				VSArray[ii1].sf_handle[ii2] = -1;
			}
			for(ii2=0; ii2<VECT_STORAGE_TMP_SIZE; ii2++){
				VSArray[ii1].tmp[ii2] = NULL;
			}
			for(ii2=0; ii2<VECT_STORAGE_BUF_SIZE; ii2++){
				VSArray[ii1].buf[ii2] = NULL;
			}
			VSArray[ii1].current_sf_index = -1;
			VSArray[ii1].dirty_flag = 0;
		}
		InitFlag = 1;
	}

	vs = NULL;
	for(ii1=0; ii1<VS_ARRAY_SIZE; ii1++){
		if(VSArray[ii1].array_size <= 0){
			vs = &(VSArray[ii1]);
			break;
		}
	}
	if(vs == NULL) return(NULL);

	vs->array_size = array_size;
	vs->vect_size = vect_size;

	vs->sf = open_scratch_file(scratch_dir);
	if(vs->sf == NULL) return(NULL);
	for(ii1=0; ii1<SCRATCH_DATA_SIZE; ii1++){
		vs->sf_handle[ii1] = -1;
	}

	for(ii1=0; ii1<VECT_STORAGE_TMP_SIZE; ii1++){
		if( allocate1D(vect_size, &(vs->tmp[ii1])) != NORMAL_RC) return(NULL);
	}

	alloc_size = VECT_STORAGE_BUF_SIZE;
	if( allocate1D(alloc_size*vect_size, &(vs->buf[0])) != NORMAL_RC){
		return(NULL);
	}
	for(ii1=1; ii1<alloc_size; ii1++){
		vs->buf[ii1] = &(vs->buf[0][ii1*vect_size]);
	}
	vs->current_sf_index = 0;
	vs->dirty_flag = 0;

	return(vs);
}


RC
close_vect_storage (VECT_STORAGE *vs)
{
	int ii1;
	int alloc_size;

	if(vs->sf != NULL){
		RC_TRY( close_scratch_file(vs->sf) );
	}
	vs->sf = NULL;

	for(ii1=0; ii1<SCRATCH_DATA_SIZE; ii1++){
		vs->sf_handle[ii1] = -1;
	}

	for(ii1=0; ii1<VECT_STORAGE_TMP_SIZE; ii1++){
		if(vs->tmp[ii1] != NULL){
			RC_TRY( free1D(vs->vect_size, &(vs->tmp[ii1])) );
		}
		vs->tmp[ii1] = NULL;
	}

	alloc_size = VECT_STORAGE_BUF_SIZE;
	RC_TRY( free1D(alloc_size*(vs->vect_size), &(vs->buf[0])) );
	for(ii1=0; ii1<VECT_STORAGE_BUF_SIZE; ii1++){
		vs->buf[ii1] = NULL;
	}
	vs->current_sf_index = -1;
	vs->dirty_flag = 0;

	vs->array_size = 0;
	vs->vect_size = 0;

	return(NORMAL_RC);
}


double *
tmp_ptr_vect_storage (VECT_STORAGE *vs, int tmp_index)
{
	if(vs->array_size <= 0) return(NULL);
	if( (tmp_index < 0)||(tmp_index >= VECT_STORAGE_TMP_SIZE) ){
		return(NULL);
	}

	return(vs->tmp[tmp_index]);
}


double **
buf_ptr_vect_storage (VECT_STORAGE *vs)
{
	if(vs->array_size <= 0) return(NULL);

	return(vs->buf);
}


/* vect[vect_size] を index に書き込み */
RC
direct_write_vect_storage (VECT_STORAGE *vs, const double vect[], int index)
{
	int ii1;
	int sf_index = index/VECT_STORAGE_BUF_SIZE;
	int buf_index = index%VECT_STORAGE_BUF_SIZE;
	FILE *fp;

	if(vs->array_size <= 0) return(ARG_ERROR_RC);
	if(vect == NULL) return(ARG_ERROR_RC);
	if( (index < 0)||(index >= vs->array_size) ) return(ARG_ERROR_RC);

	/* バッファ上にデータがあればそれも更新 */
	if(sf_index == vs->current_sf_index){
		for(ii1=0; ii1<(vs->vect_size); ii1++){
			vs->buf[buf_index][ii1] = vect[ii1];
		}
		/* 更新内容はスクラッチファイル上にも反映されるので */
		/* dirty_flag は変更しない */
	}

	if(vs->sf_handle[sf_index] >= 0){
		/* スクラッチファイル上の既存のデータを上書き */
		RC_NULL_CHK( fp = write_scratch_file(vs->sf, vs->sf_handle[sf_index]) );
		RC_TRY( wrap64_fseek(fp, buf_index*vs->vect_size*sizeof(double),
		                                                       SEEK_CUR) );
		if( fwrite(vect, sizeof(double), vs->vect_size, fp)
		                                 != (size_t)(vs->vect_size) ){
			return(WRITE_ERROR_RC);
		}
	}else{
		/* スクラッチファイル上に新規に書き込む */
		RC_NULL_CHK( fp = add_scratch_file_start(vs->sf) );
		if(sf_index == vs->current_sf_index){
			/* バッファ上にデータがあればその内容を書き込む */
			if( fwrite(vs->buf[0], sizeof(double),
			           VECT_STORAGE_BUF_SIZE*(vs->vect_size), fp)
			    != (size_t)(VECT_STORAGE_BUF_SIZE*(vs->vect_size)) ){
				return(WRITE_ERROR_RC);
			}
			vs->dirty_flag = 0;
		}else{
			/* バッファ上にデータが無ければ vect を直接書き込む */
			/* index 以外はダミーのデータを書き込む */
			for(ii1=0; ii1<VECT_STORAGE_BUF_SIZE; ii1++){
				void *ptr = (void *)(vs->tmp[VECT_STORAGE_TMP_SIZE - 1]);
				                                                 /* ダミー */

				if(ii1 == buf_index) ptr = (void *)vect;
				if( fwrite(ptr, sizeof(double), vs->vect_size, fp)
				                                  != (size_t)(vs->vect_size) ){
					return(WRITE_ERROR_RC);
				}
			}
		}
		vs->sf_handle[sf_index] = add_scratch_file_end(vs->sf);
		RC_NEG_CHK( vs->sf_handle[sf_index] );
	}

	return(NORMAL_RC);
}


/* vect[vect_size] に index のデータを読み込み */
RC
direct_read_vect_storage (VECT_STORAGE *vs, double vect[], int index)
{
	int ii1;
	int sf_index = index/VECT_STORAGE_BUF_SIZE;
	int buf_index = index%VECT_STORAGE_BUF_SIZE;
	FILE *fp;

	if(vs->array_size <= 0) return(ARG_ERROR_RC);
	if(vect == NULL) return(ARG_ERROR_RC);
	if( (index < 0)||(index >= vs->array_size) ) return(ARG_ERROR_RC);

	/* バッファ上にデータがあればそこからコピー */
	if(sf_index == vs->current_sf_index){
		for(ii1=0; ii1<(vs->vect_size); ii1++){
			vect[ii1] = vs->buf[buf_index][ii1];
		}
		return(NORMAL_RC);
	}

	if(vs->sf == NULL) return(READ_ERROR_RC);
	if(vs->sf_handle[sf_index] < 0) return(READ_ERROR_RC);

	/* スクラッチファイル上のデータを読み込み */
	RC_NULL_CHK( fp = read_scratch_file(vs->sf, vs->sf_handle[sf_index]) );
	RC_TRY(wrap64_fseek(fp, buf_index*vs->vect_size*sizeof(double),
	                                                      SEEK_CUR) );
	if( fread(vect, sizeof(double), vs->vect_size, fp)
	                        != (size_t)(vs->vect_size) ){
		return(READ_ERROR_RC);
	}

	return(NORMAL_RC);
}


/* vect[vect_size] を index に書き込み */
RC
buffered_write_vect_storage (VECT_STORAGE *vs, const double vect[], int index)
{
	int ii1;
	int buf_index = index%VECT_STORAGE_BUF_SIZE;

	if(vs->array_size <= 0) return(ARG_ERROR_RC);
	if(vect == NULL) return(ARG_ERROR_RC);
	if( (index < 0)||(index >= vs->array_size) ) return(ARG_ERROR_RC);
	if(vs->current_sf_index < 0) return(ARG_ERROR_RC);

	/* バッファ上に index を含むデータを読み込む */
	RC_TRY( get_buf_vect_storage(vs, index) );

	/* バッファ上のデータを更新 */
	for(ii1=0; ii1<(vs->vect_size); ii1++){
		vs->buf[buf_index][ii1] = vect[ii1];
	}
	vs->dirty_flag = 1;

	return(NORMAL_RC);
}


/* vect[vect_size] に index のデータを読み込み */
RC
buffered_read_vect_storage (VECT_STORAGE *vs, double vect[], int index)
{
	int ii1;
	int buf_index = index%VECT_STORAGE_BUF_SIZE;

	if(vs->array_size <= 0) return(ARG_ERROR_RC);
	if(vect == NULL) return(ARG_ERROR_RC);
	if( (index < 0)||(index >= vs->array_size) ) return(ARG_ERROR_RC);
	if(vs->current_sf_index < 0) return(ARG_ERROR_RC);

	/* バッファ上に index を含むデータを読み込む */
	RC_TRY( get_buf_vect_storage(vs, index) );

	/* バッファ上のデータからコピー */
	for(ii1=0; ii1<(vs->vect_size); ii1++){
		vect[ii1] = vs->buf[buf_index][ii1];
	}

	return(NORMAL_RC);
}


/* バッファ上のデータをスクラッチファイルに書き込み */
RC
flush_vect_storage (VECT_STORAGE *vs)
{
	FILE *fp;

	if(vs->dirty_flag == 0) return(NORMAL_RC);
	if(vs->current_sf_index < 0) return(NORMAL_RC);

	if(vs->sf_handle[vs->current_sf_index] < 0){
		RC_NULL_CHK( fp = add_scratch_file_start(vs->sf) );
	}else{
		RC_NULL_CHK( fp = write_scratch_file(vs->sf,
		                        vs->sf_handle[vs->current_sf_index]) );
	}
	if( fwrite(vs->buf[0], sizeof(double),
	           VECT_STORAGE_BUF_SIZE*(vs->vect_size), fp)
	        != (size_t)(VECT_STORAGE_BUF_SIZE*(vs->vect_size)) ){
		return(WRITE_ERROR_RC);
	}
	if(vs->sf_handle[vs->current_sf_index] < 0){
		RC_NEG_CHK( vs->sf_handle[vs->current_sf_index]
			               = add_scratch_file_end(vs->sf) );
	}
	vs->dirty_flag = 0;

	return(NORMAL_RC);
}


/* スクラッチファイルからバッファに index を含むデータを読み込む */
/* 必要に応じて自動的に既存のデータをスクラッチファイルに退避する */
RC
get_buf_vect_storage (VECT_STORAGE *vs, int index)
{
	FILE *fp;
	int sf_index = index/VECT_STORAGE_BUF_SIZE;
	int ii1;

	if(vs->current_sf_index < 0) return(ARG_ERROR_RC);

	/* バッファ上に該当データが既にある */
	if(sf_index == vs->current_sf_index){
		return(NORMAL_RC);
	}

	/* 既存のデータを書き込み */
	if(vs->dirty_flag == 1) RC_TRY( flush_vect_storage(vs) );

	if(vs->sf == NULL) return(READ_ERROR_RC);
	if(vs->sf_handle[sf_index] < 0){
		/* スクラッチファイル上にデータが存在しない場合は，*/
		/* ダミーデータでバッファを埋める */
		for(ii1=0; ii1<VECT_STORAGE_BUF_SIZE*(vs->vect_size); ii1++){
			vs->buf[0][ii1] = 0.0;
		}
		vs->dirty_flag = 0;
		vs->current_sf_index = sf_index;

		return(NORMAL_RC);
	}

	RC_NULL_CHK( fp = read_scratch_file(vs->sf, vs->sf_handle[sf_index]) );
	if( fread(vs->buf[0], sizeof(double),
	          VECT_STORAGE_BUF_SIZE*(vs->vect_size), fp)
	       != (size_t)(VECT_STORAGE_BUF_SIZE*(vs->vect_size)) ){
		return(READ_ERROR_RC);
	}
	vs->dirty_flag = 0;
	vs->current_sf_index = sf_index;

	return(NORMAL_RC);
}


/* buf をテンポラリとして使用開始 */
/* buf の内容は維持される */
double **
start_use_buf4tmp_vect_storage (VECT_STORAGE *vs)
{
	if(vs->array_size <= 0) return(NULL);
	if(vs->current_sf_index < 0) return(NULL); /* 既に使用中 */

	/* 既存のデータを書き込み */
	if(vs->dirty_flag == 1){
		if(flush_vect_storage(vs) != NORMAL_RC) return(NULL);
	}

	vs->current_sf_index = -1;

	return(vs->buf);
}


/* buf をテンポラリとして使用終了 */
/* buf の内容を w_index に書き込み，r_index の内容を読み込む */
/* w_index < 0 なら書き込まない */
RC
stop_use_buf4tmp_vect_storage (VECT_STORAGE *vs, int w_index, int r_index)
{
	if(vs->array_size <= 0) return(ARG_ERROR_RC);
	if(vs->current_sf_index >= 0) return(ARG_ERROR_RC);
	if(w_index >= vs->array_size) return(ARG_ERROR_RC);

	if(w_index >= 0){
		vs->current_sf_index = w_index/VECT_STORAGE_BUF_SIZE;
		vs->dirty_flag = 1;
		RC_TRY( flush_vect_storage(vs) );
	}else{
		/* get_buf_vect_storage() を呼び出すためのダミー処理 */
		vs->current_sf_index = (r_index/VECT_STORAGE_BUF_SIZE) + 1;
		vs->dirty_flag = 0;
	}

	RC_TRY( get_buf_vect_storage(vs, r_index) );

	return(NORMAL_RC);
}


/* A^T * B * A => C */
RC
mul_nonzero3s_AtBA_vect_storage (VECT_STORAGE *A, NONZERO_MATRIX3 B,
                                 double **C)
{
	int ii1, ii2, ii3;
	double *tmpv;
	double **buf;
	int dof = 3*B.size;
	int buf_size;
	int next_buf;

	if(A->vect_size != 3*B.size) return(ARG_ERROR_RC);

	RC_NULL_CHK( tmpv = tmp_ptr_vect_storage(A, 0) );

	for(ii1=0; ii1<A->array_size; ii1+=VECT_STORAGE_BUF_SIZE){
		/* ii1 : C[][] の列方向ループ */
		/*       VECT_STORAGE_BUF_SIZE でストリップマイニング */

		RC_TRY( get_buf_vect_storage(A, ii1) );
		RC_NULL_CHK( buf = start_use_buf4tmp_vect_storage(A) );

		/* buf * B => buf */
		buf_size = MIN2(VECT_STORAGE_BUF_SIZE, A->array_size - ii1);
		for(ii2=0; ii2<buf_size; ii2++){
			RC_TRY( mul_nonzero_matrix_vect3s(B, buf[ii2], tmpv) );
			for(ii3=0; ii3<dof; ii3++){
				buf[ii2][ii3] = tmpv[ii3];
			}
		}

		for(ii2=0; ii2<A->array_size; ii2++){
			/* C[][] の行方向ループ */
			/* ただし，OS のディスクキャッシュを期待して */
			/* 昇順，降順を交互に繰り返す */
			int ii2_sub = ii2;
			if((ii1/VECT_STORAGE_BUF_SIZE)%2 == 1){
				ii2_sub = A->array_size - 1 - ii2;
			}

			RC_TRY( direct_read_vect_storage(A, tmpv, ii2_sub) );
			for(ii3=0; ii3<buf_size; ii3++){
				C[ii1+ii3][ii2_sub] = inner_product(dof, buf[ii3], tmpv);
			}
		}

		next_buf = ii1 + VECT_STORAGE_BUF_SIZE;
		if(next_buf >= A->array_size) next_buf = 0;
		RC_TRY( stop_use_buf4tmp_vect_storage(A, -1, next_buf) );
	}

	return(NORMAL_RC);
}


RC
B_orthogonalize_vect_storage (VECT_STORAGE *A, NONZERO_MATRIX3 B)
{
	int ii1, ii2, ii3, ii4, ii5;
	double *tmpv;
	double vtBv, isq_vtBv;
	double v2tBv;
	double **buf;
	int dof = 3*B.size;
	int buf_size;
	double *Bv[VECT_STORAGE_TMP_SIZE - 1];
	const int Bv_size = sizeof(Bv)/sizeof(Bv[0]);
	int Bv_index;
	double *Bv_active;
	int next_buf;

	if(A->vect_size != 3*B.size) return(ARG_ERROR_RC);

	for(ii1=0; ii1<Bv_size; ii1++){
		RC_NULL_CHK( Bv[ii1] = tmp_ptr_vect_storage(A, ii1) );
	}
	RC_NULL_CHK( tmpv = tmp_ptr_vect_storage(A, Bv_size) );

	for(ii1=0; ii1<A->array_size; ii1+=VECT_STORAGE_BUF_SIZE){
		RC_TRY( get_buf_vect_storage(A, ii1) );
		RC_NULL_CHK( buf = start_use_buf4tmp_vect_storage(A) );

		buf_size = MIN2(VECT_STORAGE_BUF_SIZE, A->array_size - ii1);
		Bv_index = 0;
		for(ii2=0; ii2<buf_size; ii2++){
			Bv_active = Bv[Bv_index];

			/* buf 内のベクトルを直行化 */
			RC_TRY( mul_nonzero_matrix_vect3s(B, buf[ii2], Bv_active) );
			vtBv = inner_product(dof, buf[ii2], Bv_active);
			if(vtBv < ABS_TOL){
				RC_TRY( fill_random_array(dof, buf[ii2], ii1+ii2) );
				RC_TRY( mul_nonzero_matrix_vect3s(B, buf[ii2], Bv_active) );
				vtBv = inner_product(dof, buf[ii2], Bv_active);
				if(vtBv < ABS_TOL) return(CAL_ERROR_RC);
			}
			isq_vtBv = 1.0/(sqrt(vtBv));
			for(ii3=0; ii3<dof; ii3++){
				buf[ii2][ii3] *= isq_vtBv;
				Bv_active[ii3] *= isq_vtBv;
			}
			for(ii3=ii2+1; ii3<buf_size; ii3++){
				v2tBv = inner_product(dof, buf[ii3], Bv_active);
				for(ii4=0; ii4<dof; ii4++){
					buf[ii3][ii4] -= v2tBv * buf[ii2][ii4];
				}
			}

			Bv_index++;
			if((Bv_index >= Bv_size)||(ii2 == buf_size - 1) ){
				/* ディスク上のベクトルを直行化 */
				for(ii3=ii1+VECT_STORAGE_BUF_SIZE; ii3<A->array_size; ii3++){
					RC_TRY( direct_read_vect_storage(A, tmpv, ii3) );

					for(ii4=0; ii4<Bv_index; ii4++){
						v2tBv = inner_product(dof, tmpv, Bv[ii4]);
						for(ii5=0; ii5<dof; ii5++){
							tmpv[ii5]
							    -= v2tBv*buf[ii2 - Bv_index + ii4 + 1][ii5];
						}
					}
					RC_TRY( direct_write_vect_storage(A, tmpv, ii3) );
				}
				Bv_index = 0;
			}
		}

		next_buf = ii1 + VECT_STORAGE_BUF_SIZE;
		if(next_buf >= A->array_size) next_buf = 0;
		RC_TRY( stop_use_buf4tmp_vect_storage(A, ii1, next_buf) );
	}

	return(NORMAL_RC);
}


RC
compaction_vect_storage (VECT_STORAGE *vs)
{
	int ii1;
	void *ptr[VECT_STORAGE_TMP_SIZE + 1];

	for(ii1=0; ii1<VECT_STORAGE_TMP_SIZE; ii1++){
		ptr[ii1] = vs->tmp[ii1];
	}
	ptr[VECT_STORAGE_TMP_SIZE] = vs->buf[0];

	RC_TRY( mm_compaction_arr(VECT_STORAGE_TMP_SIZE + 1, ptr) );

	for(ii1=0; ii1<VECT_STORAGE_TMP_SIZE; ii1++){
		vs->tmp[ii1] = ptr[ii1];
	}

	vs->buf[0] = ptr[VECT_STORAGE_TMP_SIZE];
	for(ii1=1; ii1<VECT_STORAGE_BUF_SIZE; ii1++){
		/* buf[1?] も書き換えが必要 */
		vs->buf[ii1] = &(vs->buf[0][ii1*(vs->vect_size)]);
	}

	return(NORMAL_RC);
}

