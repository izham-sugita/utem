/*********************************************************************
 * memory_manager.c
 *
 * Copyright (C) 2004 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Thanks to following contributors.
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: memory_manager.c 1065 2016-05-04 09:39:57Z hayashi $ */

/* 単一のメモリー領域内でのメモリー管理
 * 
 * ヒープ領域の断片化を防ぎながら，malloc(), realloc(), free() のような，
 * メモリーの動的な確保・解放を行うために，単一のメモリー領域内で独自に
 * メモリー管理を行う．ただし，メモリー領域内が断片化する事は避けられない
 * ので適宜コンパクションを行う必要がある．
 *
 * (メモリーの管理方法)
 * 単一のメモリー領域は初期化時に確保され，それ自体は MM 構造体で管理する．
 * メモリー領域内において上位プログラムが使用中領域，未使用領域を
 * MM_LIST 構造体で管理する．メモリー領域内は次のような順序でデータが
 * 配置される．
 *
 * -------
 * MM_LIST 構造体(1)
 * MM_LIST 構造体(1) が管理する使用中/未使用領域
 * MM_LIST 構造体(2)
 * MM_LIST 構造体(2) が管理する使用中/未使用領域
 * .....
 * MM_LIST 構造体(N)
 * MM_LIST 構造体(N) が管理する使用中/未使用領域
 * MM_LIST 構造体(ダミー)
 * -------
 *
 *  ....でも，mm_init() の引数がゼロの場合は，
 *  malloc(), realloc(), free() を呼び出すだけ．
 */

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdarg.h>
#include "rc.h"
#include "memory_manager.h"

typedef struct mm_list{
	void *this;      /* 構造体の先頭アドレスまたは，
	                    コンパクション処理時のアドレス書き換え先 */
	struct mm_list *prev;   /* 一つ前方の構造体へのポインタ */
	size_t size;     /* 管理領域の大きさ + MM_LIST_SIZE */
	int alloc_flag;  /* 0: 管理領域は未使用，1: 使用中，2: 移動可能 */
	size_t chk_sum;  /* 構造体の破壊を検出するチェックサム */
} MM_LIST;

typedef struct{
	void *ptr;        /* メモリー領域先頭のポインタ */
	MM_LIST *list;    /* 最後のダミー MM_LIST 構造体のポインタ */
	size_t size;      /* メモリー領域全体の大きさ */
	size_t max_usage; /* メモリー領域の使用量の最大値 */
} MM;

static MM mm = {NULL, NULL, 0, 0};

#ifdef WIN32
#define ALIGNMENT   (8)    /* アライメントサイズ */
#else   /* WIN32 */
#define ALIGNMENT   (16)   /* アライメントサイズ */
#endif  /* WIN32 */

/* アライメントを考慮した MM_LIST 構造体の大きさ */
#define MM_LIST_SIZE (ALIGNMENT*((sizeof(MM_LIST)+(ALIGNMENT-1))/ALIGNMENT))

/*
 * [アライメントに関して]
 * 64bit OS環境では，メモリアクセスに対して，8の倍数の大きさでアクセスした
 * 方が高速になる(少なくとも遅くならない)場合が多い．32bit OSの場合，アラ
 * イメントの大きさは 4 で良いらしい．無意味に大きくしたり，中途半端な値を
 * 入れるとメモリの使用効率が極端に悪くなる．8 or 4 が Best.
 */


#define MAX_SEEK_COUNT   (10)   /* 適切な確保位置を探査する回数 */

static void *mm_alloc_sub(size_t size, size_t max_seek_count);
static MM_LIST *cal_next_ptr(const MM_LIST *list_ptr);
static MM_LIST *cal_list_ptr(const void *array_ptr);
static void *cal_array_ptr(const MM_LIST *list_ptr);
static int is_error_list(MM_LIST *list_ptr);
static size_t cal_chk_sum(const MM_LIST *list_ptr);
static MM_LIST *merge_free_list(MM_LIST *list_ptr);
static MM_LIST *compaction(MM_LIST *list_ptr);
static RC mm_compaction(void);
static MM_LIST *allocation(size_t size, MM_LIST *list_ptr);
static void update_max_usage(void);


/* データ構造を初期化して，メモリー領域を確保する．
 * 適切にコンパクションが行われて断片化されていなければ，少なくとも次式で
 * 計算される容量が使用できる．
 *
 * 使用領域(ALLIGN_SIZE 単位で切り上げ)の合計値
 *   = size - 2*MM_LIST_SIZE*(確保する個数 - 1)
 */
RC
mm_init(size_t size)
{
	MM_LIST *list_ptr;
	MM_LIST *next_ptr;

	if(mm.ptr != NULL) return(UNKNOWN_ERROR_RC);

	if(size <= 0){
		/* MM 構造体を初期化して素通り */
		mm.ptr = NULL;
		mm.list = NULL;
		mm.size = 0;
		mm.max_usage = 0;

		return(NORMAL_RC);
	}

	/* size を，アライメントを考慮した大きさに加えて，MM_LIST構造体2つ分
	 * 大きく確保する． 
	 * これは，先頭のMM_LIST構造体および最後のダミーMM_LIST構造体の分で
	 * ある． */
	size = ALIGNMENT*((size + (ALIGNMENT-1))/ALIGNMENT) + 2*MM_LIST_SIZE;

	/* メモリー領域を確保する． */
	mm.ptr = malloc(size);
	if(mm.ptr == NULL){
		RC_TRY( error_printf(90001, (int)(size/(1024*1024))) );
		return(ALLOC_ERROR_RC);
	}
	mm.size = size;
	mm.max_usage = 0;

	/* ダミー構造体以外の領域を未使用領域として格納する． 
	 * ? これを格納する理由 */
	list_ptr = (MM_LIST *)(mm.ptr);
	list_ptr->this = list_ptr;
	list_ptr->prev = NULL;
	list_ptr->size = size - MM_LIST_SIZE;
	list_ptr->alloc_flag = 0;
	list_ptr->chk_sum = cal_chk_sum(mm.ptr);

	/* 最後のダミー構造体を格納する． */
	next_ptr = cal_next_ptr(mm.ptr);
	if(next_ptr == NULL) return(UNKNOWN_ERROR_RC);
	next_ptr->this = next_ptr;
	next_ptr->prev = mm.ptr;
	next_ptr->size = 0;
	next_ptr->alloc_flag = 0;
	next_ptr->chk_sum = cal_chk_sum(next_ptr);

	mm.list = next_ptr;

	return(mm_chk_list(NULL));
}


/* すでに確保している領域を解放して，size[Bytes] で再確保する */
RC
mm_resize(size_t size)
{
	RC_TRY( mm_terminate() );

	return(mm_init(size));
}


/* すでに確保している領域を解放する．*/
RC
mm_terminate(void)
{
	if(mm.size == 0){
		return(NORMAL_RC);    /* 確保されていない */
	}

	if(mm.ptr == NULL) return(UNKNOWN_ERROR_RC);

	/* MM構造体ごと解放する． */
	free(mm.ptr);
	mm.ptr = NULL;
	mm.list = NULL;
	mm.size = 0;
	mm.max_usage = 0;

	return(NORMAL_RC);
}


/* size[Bytes] 確保して，そのポインタを返却する．
 * 確保に失敗すれば NULL を返却する．
 */
void *
mm_alloc(size_t size)
{
	return(mm_alloc_sub(size, MAX_SEEK_COUNT));
}


/* mm_alloc() の 64bit 引数版．
 * size_t が 32bit の環境で，オーバーフローの可能性がある場合に使用する．
 */
void *
mm_alloc64(WRAP64_INT size)
{
	if(sizeof(size_t) < 8){
		const size_t max_size_t = ~((size_t)0);

		if(size > (WRAP64_INT)max_size_t){
			error_printf(90003);
			return(NULL);
		}
	}

	return(mm_alloc_sub((size_t)size, MAX_SEEK_COUNT));
}


/* mm_alloc() と同じだが，極力後方の領域を使用する．
 */
void *
mm_alloc_tmp(size_t size)
{
	return(mm_alloc_sub(size, 0));
}


/* 後方の空き領域の大きさを返す．
 */
size_t
mm_last_allocatable(void)
{
	MM_LIST *list_ptr = mm.list;

	if(mm.size == 0) return(ULONG_MAX);

	if(is_error_list(list_ptr)) return(0);
	list_ptr = list_ptr->prev;
	if(is_error_list(list_ptr)) return(0);
	if(list_ptr->alloc_flag == 1) return(0);
	if(list_ptr->size <= MM_LIST_SIZE) return(0);

	return(list_ptr->size - MM_LIST_SIZE);
}


/* MM構造体の max_usage を更新する．
 */
static void
update_max_usage(void)
{
	MM_LIST *list_ptr = mm.list;
	size_t usage;

	if(mm.size == 0) return;

	if(is_error_list(list_ptr)) return;
	list_ptr = list_ptr->prev;
	if(is_error_list(list_ptr)) return;
	if(list_ptr->alloc_flag == 0){
		usage = mm.size - (list_ptr->size - MM_LIST_SIZE);
	}else{
		usage = mm.size;
	}
	if(usage > mm.max_usage) mm.max_usage = usage;

	return;
}

/* 最大 max_seek_count 回の未使用領域の探索を行い，メモリー領域を
 * 動的確保する． */
static void *
mm_alloc_sub(size_t size, size_t max_seek_count)
{
	MM_LIST *list_ptr;
	MM_LIST *best_ptr = NULL;
	size_t best_size = 0;
	size_t count = 0;
	void *array_ptr;

	/* MM構造体の大きさが 0 ならば，通常の malloc() を用いる． */
	if(mm.size == 0){
		return(malloc(size));
	}

	/* size を，アライメントと MM_LIST_SIZE を考慮した値に修正する． */
	if(size <= 0) return(NULL);
	size = ALIGNMENT*((size + (ALIGNMENT-1))/ALIGNMENT) + MM_LIST_SIZE;

	/* 後方から未使用領域を検索する． */
	list_ptr = mm.list;
	while(list_ptr != NULL){
		if(is_error_list(list_ptr)) return(NULL);

		if( (list_ptr->alloc_flag == 0)&&(list_ptr->size >= size) ){
			/* 可能な限り小さな未使用領域を選ぶ */
			if( (best_size == 0)||(list_ptr->size <= best_size) ){
				best_ptr = list_ptr;
				best_size = list_ptr->size;
			}
		}

		/* 未使用領域が見つかれば，一定回数以上検索を続けない．
		 * 未使用領域が見つかっても，少なくとも max_seek_count 回は検索を
		 * 繰り返す． */
		count++;
		if( (best_size > 0)&&(count >= max_seek_count) ) break;

		list_ptr = list_ptr->prev;
	}
	if(best_ptr == NULL){
		error_printf(90002);
		return(NULL);
	}

	list_ptr = allocation(size, best_ptr);
	if(list_ptr == NULL) return(NULL);
	array_ptr = cal_array_ptr(list_ptr);
	memset(array_ptr, 0, size - MM_LIST_SIZE);
	update_max_usage();

	return(array_ptr);
}


/* mm_alloc() で確保された領域 ptr を size[Bytes] に拡大または縮小し，
 * そのポインタを返却する．
 * 失敗すれば NULL を返却する．
 *
 * 再確保の手続きは以下の通り．
 *
 * 1. size が以前より小さければ，以下 1-1. ~ を実行する．
 *    1-1. 後方の領域が未使用であれば，領域を縮小した後で後方の領域を前方に
 *         拡大する．
 *
 *    (以下は，size が以前より大きい場合となる）
 *
 * 2. 後方の領域が未使用かつ拡大に十分な大きさがあれば，以下 2-1 ~ を実行する．
 *    2-1. 後方の未使用領域に拡大したうえでMM_LIST構造体を確保する大きさが
 *         あれば，余った領域を未使用領域として管理する．
 *    2-2. MM_LIST構造体を確保する大きさがなければ，未使用領域のすべてを管理
 *         領域とする．
 *
 * 3. 以上の処理で拡大ができなければ，以下 3-1. ~ を実行する．
 *    3-1. 前方の未使用領域から，拡大可能かつなるべく小さい未使用領域を検索する．
 *    3-2. 後方の未使用領域から，拡大可能かつなるべく小さい未使用領域を検索する．
 *    3-3. 移動先が見つかれば，確保・移動・解放を行う．
 *    3-4. 見つからなければ，mm_alloc() で確保し，移動・解放を行う．
 */
void *
mm_realloc(void *ptr, size_t size)
{
	MM_LIST *list_ptr;
	MM_LIST *next_ptr;
	MM_LIST *prev_ptr;
	void *ret;
	size_t old_size;
	size_t next_size;
	MM_LIST *best_ptr;
	size_t best_size;
	size_t count;

	/* 通常の realloc() */
	if(mm.size == 0){
		return(realloc(ptr, size)); 
	}
	/* 新しい size が 0 なら，解放する． */
	if(size <= 0){
		mm_free(ptr);
		return(NULL);
	}
	/* 確保されてなければ，mm_alloc() で確保する． */
	if(ptr == NULL){
		return(mm_alloc(size));
	}

	/* size を修正する． */
	size = ALIGNMENT*((size + (ALIGNMENT-1))/ALIGNMENT) + MM_LIST_SIZE;
	/* MM_LIST構造体のポインタを取得する． */
	list_ptr = cal_list_ptr(ptr);
	if(is_error_list(list_ptr)) return(NULL);
	/* 前方に詰める */
	list_ptr = compaction(list_ptr);
	if(list_ptr == NULL) return(NULL);

	/* 元々のsizeより小さければ，領域を縮小する． */
	if(list_ptr->size >= size){    
		old_size = list_ptr->size;
		MM_LIST *new_next_ptr;

		/* 縮小される領域のメモリを初期化する． */
		memset((char *)list_ptr + size, 0, list_ptr->size - size);

		/* 後方の領域が未使用であれば，前方に拡大する． */
		next_ptr = cal_next_ptr(list_ptr);
		if(next_ptr == NULL) return(NULL);
		if((next_ptr->alloc_flag == 0)&&(next_ptr->size > 0)){
			next_size = next_ptr->size + (list_ptr->size - size);
			MM_LIST *next_next_ptr = cal_next_ptr(next_ptr);

			/* 領域を縮小する． */
			list_ptr->size = size;
			list_ptr->chk_sum = cal_chk_sum(list_ptr);

			/* 後方の領域を前方に拡大する． */
			next_ptr = cal_next_ptr(list_ptr);
			if(next_ptr == NULL) return(NULL);
			next_ptr->this = next_ptr;
			next_ptr->prev = list_ptr;
			next_ptr->size = next_size;
			next_ptr->alloc_flag = 0;
			next_ptr->chk_sum = cal_chk_sum(next_ptr);

			/* 2つ後方のMM_LIST構造体を更新する． */
			if(next_next_ptr == NULL) return(NULL);
			next_next_ptr->prev = next_ptr;
			next_next_ptr->chk_sum = cal_chk_sum(next_next_ptr);

			return(cal_array_ptr(list_ptr));
		}

		/* 後方の領域が使用中かつ縮小後の空いた領域にMM_LIST構造体が収まらない
		 * 場合は，サイズを変更しない． */
		if(list_ptr->size <= size + MM_LIST_SIZE){
			return(cal_array_ptr(list_ptr));
		}

		/* 収まる場合は，後方に未使用領域を作成する． */
		list_ptr->size = size;
		list_ptr->alloc_flag = 1;
		list_ptr->chk_sum = cal_chk_sum(list_ptr);

		new_next_ptr = cal_next_ptr(list_ptr);
		if(new_next_ptr == NULL) return(NULL);
		new_next_ptr->this = new_next_ptr;
		new_next_ptr->prev = list_ptr;
		new_next_ptr->size = old_size - size;
		new_next_ptr->alloc_flag = 0;
		new_next_ptr->chk_sum = cal_chk_sum(new_next_ptr);

		next_ptr->prev = new_next_ptr;
		next_ptr->chk_sum = cal_chk_sum(next_ptr);

		return(cal_array_ptr(list_ptr));
	}

	/* 元々のsizeより大きければ，領域を拡大する． */
	old_size = list_ptr->size;
	next_ptr = cal_next_ptr(list_ptr);
	if(next_ptr == NULL) return(NULL);

	/* 後方の領域が未使用かつ，十分に拡大できる大きさがあれば，それを利用
	 * する． */
	if( (next_ptr->alloc_flag == 0)
	  &&(next_ptr->size > size - list_ptr->size) ){
		/* 後方の未使用領域に拡大したうえで，MM_LIST構造体を確保できる大きさが
		 * ある場合 */
		if(next_ptr->size > MM_LIST_SIZE + size - list_ptr->size){
			/* 後方の未使用領域を縮小する． */
			next_size = next_ptr->size - (size - list_ptr->size);
			MM_LIST *next_next_ptr = cal_next_ptr(next_ptr);

			list_ptr->size = size;
			list_ptr->chk_sum = cal_chk_sum(list_ptr);
			memset((char *)list_ptr + old_size, 0, size - old_size);

			/* 縮小した領域の後方の領域は新たにMM_LIST構造体を作って管理
			 * する． */
			next_ptr = cal_next_ptr(list_ptr);
			if(next_ptr == NULL) return(NULL);
			next_ptr->this = next_ptr;
			next_ptr->prev = list_ptr;
			next_ptr->size = next_size;
			next_ptr->alloc_flag = 0;
			next_ptr->chk_sum = cal_chk_sum(next_ptr);

			if(next_next_ptr == NULL) return(NULL);
			next_next_ptr->prev = next_ptr;
			next_next_ptr->chk_sum = cal_chk_sum(next_next_ptr);

			return(cal_array_ptr(list_ptr));
		}
		/* MM_LIST構造体を確保できない大きさであれば，未使用領域のすべてを拡大
		 * 対象にする． */
		list_ptr->size += next_ptr->size;
		list_ptr->chk_sum = cal_chk_sum(list_ptr);
		memset((char *)list_ptr + old_size, 0, size - old_size);

		next_ptr = cal_next_ptr(list_ptr);
		if(next_ptr == NULL) return(NULL);
		next_ptr->prev = list_ptr;
		next_ptr->chk_sum = cal_chk_sum(next_ptr);

		return(cal_array_ptr(list_ptr));
	}

	/* 前方の未使用領域を検索 */
	best_ptr = NULL;
	best_size = 0;
	count = 0;
	prev_ptr = list_ptr->prev;
	while(prev_ptr != NULL){
		if(is_error_list(prev_ptr)) return(NULL);

		if( (prev_ptr->alloc_flag == 0)&&(prev_ptr->size >= size) ){
			if( (best_size == 0)||(prev_ptr->size <= best_size) ){
				best_ptr = prev_ptr;
				best_size = prev_ptr->size;
			}
		}
		count++;
		if(count >= MAX_SEEK_COUNT) break;

		prev_ptr = prev_ptr->prev;
	}

	/* 後方の未使用領域を検索 */
	count = 0;
	next_ptr = cal_next_ptr(list_ptr);
	while(next_ptr != NULL){
		if(is_error_list(next_ptr)) return(NULL);

		if( (next_ptr->alloc_flag == 0)&&(next_ptr->size >= size) ){
			if( (best_size == 0)||(next_ptr->size < best_size) ){
				best_ptr = next_ptr;
				best_size = next_ptr->size;
			}
		}
		count++;
		if(count >= MAX_SEEK_COUNT) break;

		next_ptr = cal_next_ptr(next_ptr);
	}

	/* 移動先が見つかれば，新たに確保して移動し，元々の領域を解放する． */
	if(best_ptr != NULL){
		best_ptr = allocation(size, best_ptr);
		if(best_ptr == NULL) return(NULL);
		ret = cal_array_ptr(best_ptr);
		memmove(ret, cal_array_ptr(list_ptr), list_ptr->size - MM_LIST_SIZE);
		memset((char *)best_ptr + list_ptr->size, 0, size - list_ptr->size);
		mm_free(cal_array_ptr(list_ptr));
	}
	/* 見つからなければ，mm_alloc() で確保して移動し，元々の領域を解放する． */
	else{
		ret = mm_alloc(size-MM_LIST_SIZE);
		if(ret == NULL) return(NULL);
		memmove(ret, cal_array_ptr(list_ptr), list_ptr->size - MM_LIST_SIZE);
		mm_free(cal_array_ptr(list_ptr));
	}

	update_max_usage();

	return(ret);
}


/* mm_realloc() の 64bit 引数版．
 * size_t が 32bit の環境で，オーバーフローの可能性がある場合に使用する．
 */
void *
mm_realloc64(void *ptr, WRAP64_INT size)
{
	if(sizeof(size_t) < 8){
		const size_t max_size_t = ~((size_t)0);

		if(size > (WRAP64_INT)max_size_t){
			error_printf(90003);
			return(NULL);
		}
	}

	return( mm_realloc(ptr, (size_t)size) );
}


/* mm_alloc(), mm_realloc() で確保された領域 ptr を解放する．
 * 解放に成功すれば NORMAL_RC を，失敗すればその他を返却する．
 */
RC
mm_free(void *ptr)
{
	MM_LIST *list_ptr;

	if(mm.size == 0){
		free(ptr); /* 通常の free() */
		return(NORMAL_RC);
	}

	if(ptr == NULL) return(ARG_ERROR_RC);
	list_ptr = cal_list_ptr(ptr);
	if(is_error_list(list_ptr)) return(UNKNOWN_ERROR_RC);
	if(list_ptr->alloc_flag != 1) return(ARG_ERROR_RC);

	/* 領域を解放 */
	list_ptr->alloc_flag = 0;
	list_ptr->chk_sum = cal_chk_sum(list_ptr);
	memset(cal_array_ptr(list_ptr), 0, list_ptr->size - MM_LIST_SIZE);

	/* 解放された領域前後の未使用領域をマージ */
	list_ptr = merge_free_list(list_ptr);
	if(list_ptr == NULL) return(UNKNOWN_ERROR_RC);

	return(NORMAL_RC);
}


/* mm_alloc(), mm_realloc() で確保された領域を前方に移動し，
 * 移動後のポインタを返却する．
 * ptr とそれに続く引数は，各領域のポインタを格納している変数のポインタで，
 * 引数の最後は NULL とする．
 * 移動に成功すれば NORMAL_RC を，失敗すればその他を返却する．
 */
RC
mm_compaction_va(void **ptr, ...)
{
	va_list ap;
	void **array_ptr_ptr;
	MM_LIST *list_ptr;

	if(mm.size == 0) return(NORMAL_RC);

	/* 移動対象にフラグを立てる */
	va_start(ap, ptr);
	array_ptr_ptr = ptr;
	while(array_ptr_ptr != NULL){
		if(*array_ptr_ptr != NULL){
			list_ptr = cal_list_ptr(*array_ptr_ptr);
			if(is_error_list(list_ptr)) return(UNKNOWN_ERROR_RC);
			if(list_ptr->alloc_flag == 1){
				list_ptr->alloc_flag = 2;
				list_ptr->this = (void *)(array_ptr_ptr);
				list_ptr->chk_sum = cal_chk_sum(list_ptr);
			}else{
				return(ARG_ERROR_RC);
			}
		}
		array_ptr_ptr = va_arg(ap, void **);
	}
	va_end(ap);

	return(mm_compaction());
}


/* mm_alloc(), mm_realloc() で確保された領域 ptr[size] を前方に移動し，
 * 移動後のポインタで ptr[] を上書きする．
 * 移動に成功すれば NORMAL_RC を，失敗すればその他を返却する．
 */
RC
mm_compaction_arr(int ptr_num, void *ptr[])
{
	int ii1;
	MM_LIST *list_ptr;

	if(mm.size == 0) return(NORMAL_RC);

	/* 移動対象にフラグを立てる */
	for(ii1=0; ii1<ptr_num; ii1++){
		if(ptr[ii1] == NULL) continue;
		list_ptr = cal_list_ptr(ptr[ii1]);
		if(is_error_list(list_ptr)) return(UNKNOWN_ERROR_RC);
		if(list_ptr->alloc_flag == 1){
			list_ptr->alloc_flag = 2;
			list_ptr->this = (void *)(&(ptr[ii1]));
			list_ptr->chk_sum = cal_chk_sum(list_ptr);
		}else{
			return(ARG_ERROR_RC);
		}
	}

	return(mm_compaction());
}


/* mm_compaction_arr() の２配列版 */
RC
mm_compaction_arr2(int ptr1_num, void *ptr1[], int ptr2_num, void *ptr2[])
{
	int ii1;
	MM_LIST *list_ptr;

	if(mm.size == 0) return(NORMAL_RC);

	/* 移動対象にフラグを立てる */
	for(ii1=0; ii1<ptr1_num; ii1++){
		if(ptr1[ii1] == NULL) continue;
		list_ptr = cal_list_ptr(ptr1[ii1]);
		if(is_error_list(list_ptr)) return(UNKNOWN_ERROR_RC);
		if(list_ptr->alloc_flag == 1){
			list_ptr->alloc_flag = 2;
			list_ptr->this = (void *)(&(ptr1[ii1]));
			list_ptr->chk_sum = cal_chk_sum(list_ptr);
		}else{
			return(ARG_ERROR_RC);
		}
	}

	/* 移動対象にフラグを立てる */
	for(ii1=0; ii1<ptr2_num; ii1++){
		if(ptr2[ii1] == NULL) continue;
		list_ptr = cal_list_ptr(ptr2[ii1]);
		if(is_error_list(list_ptr)) return(UNKNOWN_ERROR_RC);
		if(list_ptr->alloc_flag == 1){
			list_ptr->alloc_flag = 2;
			list_ptr->this = (void *)(&(ptr2[ii1]));
			list_ptr->chk_sum = cal_chk_sum(list_ptr);
		}else{
			return(ARG_ERROR_RC);
		}
	}

	return(mm_compaction());
}


/* 移動対象を前方に詰める． 
 */
static RC
mm_compaction(void)
{
	/* ポインタをメモリー領域の先頭に合わせる． */
	MM_LIST *list_ptr = mm.ptr;

	if(mm.size == 0) return(NORMAL_RC);

	while(list_ptr != NULL){
		if(is_error_list(list_ptr)) return(UNKNOWN_ERROR_RC);
		/* 移動可能状態の場合のみ，前方に詰める． */
		if(list_ptr->alloc_flag == 2){
			/* 移動先（初期値は移動前） 
			 * list_ptr->this は void * 型なのに，なんで void ** 型に
			 * 変換している？*/
			void **dest_ptr = (void **)(list_ptr->this);

			list_ptr->alloc_flag = 1;
			list_ptr->this = list_ptr;
			list_ptr->chk_sum = cal_chk_sum(list_ptr);

			/* list_ptr はコンパクション後のアドレスになる． */
			list_ptr = compaction(list_ptr);
			if(list_ptr == NULL) return(UNKNOWN_ERROR_RC);
			*dest_ptr = cal_array_ptr(list_ptr);
		}
		list_ptr = cal_next_ptr(list_ptr);
	}

	return(NORMAL_RC);
}


/* MM, MM_LIST 構造体が正常かどうかを確認する．
 * fp != NULL なら，構造体の内容をファイルに出力する．
 * 構造体が正常なら NORMAL_RC を, 異常ならその他を返却する．
 */
RC
mm_chk_list(FILE *fp)
{
	MM_LIST *list_ptr;
	MM_LIST *prev_ptr = NULL;
	RC ret = NORMAL_RC;

	if(mm.size == 0) return(NORMAL_RC);

	if(fp != NULL){
		fprintf(fp, "ptr = %p, list = %p, size = %lu, MM_LIST_SIZE = %lu\n",
		        mm.ptr, mm.list, mm.size, MM_LIST_SIZE);
	}
	if( mm.list != (void *)((char *)(mm.ptr) + (mm.size - MM_LIST_SIZE)) ){
		return(UNKNOWN_ERROR_RC);
	}

	list_ptr = mm.ptr;
	if(fp != NULL){
		fprintf(fp, "   address   previous        size    checksum\n");
	}
	while(list_ptr != NULL){
		if(list_ptr->alloc_flag){
			if(fp != NULL) fprintf(fp, "A:");
		}else{
			if(fp != NULL) fprintf(fp, "F:");
		}
		if(fp != NULL){
			fprintf(fp, " %8p", list_ptr);
			if(list_ptr != list_ptr->this){
				fprintf(fp, "*");
			}else{
				fprintf(fp, " ");
			}
			fprintf(fp, " %8p ", list_ptr->prev);
			fprintf(fp, " %10lu ", list_ptr->size);
			fprintf(fp, " %10lu ", list_ptr->chk_sum);
			if(list_ptr->chk_sum != cal_chk_sum(list_ptr->this)){
				fprintf(fp, "*\n");
			}else{
				fprintf(fp, "\n");
			}
		}
		if( (list_ptr != list_ptr->this)
		  ||(list_ptr->prev != prev_ptr)
		  ||(list_ptr->chk_sum != cal_chk_sum(list_ptr)) ){
			ret = UNKNOWN_ERROR_RC;
		}
		prev_ptr = list_ptr;
		list_ptr = cal_next_ptr(list_ptr);
	}

	list_ptr = mm.list;
	while(list_ptr != NULL){
		if(fp != NULL) fprintf(fp, "%p\n", list_ptr);
		if(is_error_list(list_ptr)){
			ret = UNKNOWN_ERROR_RC;
			break;
		}
		list_ptr = list_ptr->prev;
	}

	return(ret);
}


/* 未使用領域の最大値を返却する．
 * 出力先 fp を指定すれば，管理の状況を出力する．
 */
size_t
mm_allocatable(FILE *fp)
{
	MM_LIST *list_ptr = mm.list;
	size_t ret = 0;

	if(mm.size == 0) return(ULONG_MAX);

	if(fp != NULL){
		fprintf(fp, "   address   previous        size    checksum\n");
	}

	while(list_ptr != NULL){
		if(is_error_list(list_ptr)){
			ret = 0;
			break;
		}
		if( (list_ptr->alloc_flag == 0)&&(list_ptr->size >= MM_LIST_SIZE) ){
			if((list_ptr->size - MM_LIST_SIZE) > ret){
				ret = list_ptr->size - MM_LIST_SIZE;
			}
			if(fp != NULL){
				fprintf(fp, "F: %8p ", list_ptr);
				fprintf(fp, " %8p ", list_ptr->prev);
				fprintf(fp, " %10lu ", list_ptr->size);
				fprintf(fp, " %10lu ", list_ptr->chk_sum);
				if(list_ptr->chk_sum != cal_chk_sum(list_ptr->this)){
					fprintf(fp, "*\n");
				}else{
					fprintf(fp, "\n");
				}
			}
		}
		list_ptr = list_ptr->prev;
	}

	return(ret);
}


/* list_ptr の一つ後方の MM_LIST 構造体のポインタを返却する．
 */
static MM_LIST *
cal_next_ptr (const MM_LIST *list_ptr)
{
	MM_LIST *ret;

	if(list_ptr == NULL) return(NULL);
	if(list_ptr->size < MM_LIST_SIZE) return(NULL);
	if( ((char *)(list_ptr) - (char *)(mm.ptr) + list_ptr->size)
	   > (mm.size - MM_LIST_SIZE) ){
		return(NULL);
	}
	ret = (void *)(((char *)list_ptr) + (list_ptr->size));

	return(ret);
}


/* array_ptr を管理する MM_LIST 構造体のポインタを返却する．
 */
static MM_LIST * 
cal_list_ptr (const void *array_ptr)
{
	MM_LIST *ret;

	if(array_ptr == NULL) return(NULL);
	ret = (MM_LIST *)(((char *)array_ptr) - MM_LIST_SIZE);

	return(ret);
}


/* list_ptr が管理する領域のポインタを返却する．
 */
static void *
cal_array_ptr(const MM_LIST *list_ptr)
{
	MM_LIST *ret;

	if(list_ptr == NULL) return(NULL);
	ret = (void *)(((char *)list_ptr) + MM_LIST_SIZE);

	return(ret);
}


/* list_ptr が正常かどうかを確認する．
 * 正常なら 0 を，異常なら 1 を返却する．
 */
static int
is_error_list(MM_LIST *list_ptr)
{
	/* ポインタが NULL でないかを確認する． */
	if(list_ptr == NULL){
		return(1);
	}
	/* チェックサムが一致するかを確認する． */
	if(list_ptr->chk_sum != cal_chk_sum(list_ptr)){
		return(1);
	}
	/* 確保状態が正常値かを確認する． */
	if( (list_ptr->alloc_flag != 0)
	  &&(list_ptr->alloc_flag != 1)
	  &&(list_ptr->alloc_flag != 2) ){
		return(1);
	}

	return(0);
}


/* MM_LIST構造体のチェックサムを計算する．
 */
static size_t
cal_chk_sum(const MM_LIST *list_ptr)
{
	size_t ret = 12345;

	ret ^= (size_t)(list_ptr->this);
	ret ^= (size_t)(list_ptr->prev);
	ret ^= (size_t)(list_ptr->size);
	ret ^= (size_t)(list_ptr->alloc_flag);

	return(ret);
}


/* list_ptr に隣接する未使用領域を併合する．
 * 併合後の list_ptr のアドレスを返却する．
 */
static MM_LIST *
merge_free_list(MM_LIST *list_ptr)
{
	/* 後方が未使用領域ならば併合する． */
	MM_LIST *next_ptr = cal_next_ptr(list_ptr);
	if(is_error_list(next_ptr)) return(NULL);
	if( (next_ptr->alloc_flag == 0)&&(next_ptr->size > 0) ){
		list_ptr->size += next_ptr->size;
		list_ptr->chk_sum = cal_chk_sum(list_ptr);
		memset(next_ptr, 0, MM_LIST_SIZE);
		next_ptr = cal_next_ptr(list_ptr);
		next_ptr->prev = list_ptr;
		next_ptr->chk_sum = cal_chk_sum(next_ptr);
	}

	/* 前方が未使用領域ならば併合する． */
	if(list_ptr->prev != NULL){
		if(is_error_list(list_ptr->prev)) return(NULL);
		if((list_ptr->prev)->alloc_flag == 0){
			MM_LIST *prev_ptr = list_ptr->prev;

			prev_ptr->size += list_ptr->size;
			prev_ptr->chk_sum = cal_chk_sum(list_ptr->prev);
			memset(list_ptr, 0, MM_LIST_SIZE);
			next_ptr->prev = prev_ptr;
			next_ptr->chk_sum = cal_chk_sum(next_ptr);
			list_ptr = prev_ptr;
		}
	}

	return(list_ptr);
}


/* list_ptr の前方が未使用領域ならば，前方に詰める．
 */
static MM_LIST *
compaction(MM_LIST *list_ptr)
{
	/* if(is_error_list(list_ptr)) return(NULL); 信頼できるものとして省略 */

	if(list_ptr->prev != NULL){
		if(is_error_list(list_ptr->prev)) return(NULL);
		if((list_ptr->prev)->alloc_flag == 0){
			MM_LIST *prev_ptr = list_ptr->prev;
			size_t prev_size = (list_ptr->prev)->size;
			size_t list_size = list_ptr->size;
			MM_LIST *next_ptr = cal_next_ptr(list_ptr);
			size_t st_size;

			prev_ptr->size = list_ptr->size;
			prev_ptr->alloc_flag = list_ptr->alloc_flag;
			prev_ptr->chk_sum = cal_chk_sum(prev_ptr);

			memmove(cal_array_ptr(prev_ptr),
			        cal_array_ptr(list_ptr),
			        prev_ptr->size - MM_LIST_SIZE);

			list_ptr = cal_next_ptr(prev_ptr);
			list_ptr->this = list_ptr;
			list_ptr->prev = prev_ptr;
			list_ptr->size = prev_size;
			list_ptr->alloc_flag = 0;
			list_ptr->chk_sum = cal_chk_sum(list_ptr);

			st_size = prev_size;
			if(list_size + MM_LIST_SIZE > st_size){
				st_size = list_size + MM_LIST_SIZE;
			}
			memset((char *)prev_ptr + st_size,
			       0, list_size + prev_size - st_size);

			next_ptr->prev = list_ptr;
			next_ptr->chk_sum = cal_chk_sum(next_ptr);

			list_ptr = merge_free_list(list_ptr);
			if(list_ptr == NULL) return(NULL);
			return(prev_ptr);
		}
	}

	return(list_ptr);
}


/* 未使用領域 list_ptr から size[Bytes] 確保する． 
 * 未使用領域の大きさが十分にあれば，残りを新たな未使用領域とする．
 *
 * size には MM_LIST_SIZE が考慮された値が入っていることを前提とする．
 */
static MM_LIST *
allocation(size_t size, MM_LIST *list_ptr)
{
//	if(is_error_list(list_ptr)) return(NULL); /* 信頼できるものとして省略 */
	if(list_ptr->size < size) return(NULL);

	/* 未使用領域が十分大きい場合は，残りを新たな未使用領域とする． */
	if(list_ptr->size > (size + MM_LIST_SIZE)){
		size_t old_size = list_ptr->size;
		MM_LIST *next_next_ptr = cal_next_ptr(list_ptr);
		MM_LIST *next_ptr;

		list_ptr->size = size;
		list_ptr->alloc_flag = 1;
		list_ptr->chk_sum = cal_chk_sum(list_ptr);

		/* next_ptr */
		next_ptr = cal_next_ptr(list_ptr);
		if(next_ptr == NULL) return(NULL);
		next_ptr->this = next_ptr;
		next_ptr->prev = list_ptr;
		next_ptr->size = old_size - size;
		next_ptr->alloc_flag = 0;
		next_ptr->chk_sum = cal_chk_sum(next_ptr);

		/* next_next_ptr */
		if(next_next_ptr == NULL) return(NULL);
		next_next_ptr->prev = next_ptr;
		next_next_ptr->chk_sum = cal_chk_sum(next_next_ptr);
	}
	/* 未使用領域の大きさが十分でなければ，領域全てを使用する． */
	else{
		list_ptr->alloc_flag = 1;
		list_ptr->chk_sum = cal_chk_sum(list_ptr);
	}

	return(list_ptr);
}


/* メモリ領域の最大使用量を返す．
 */
size_t
mm_max_usage(void)
{
	return(mm.max_usage);
}

/* メモリ領域全体の大きさを返す．
 */
size_t
mm_size(void)
{
	return(mm.size);
}

/* MM_LIST_SIZE の値を返す．
 */
size_t
mm_list_size(void)
{
	size_t ret = MM_LIST_SIZE;

	return(ret);
}



