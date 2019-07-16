/*
 * extract_geom_contact_surface() のボツアイディアバージョン
 * 区分木/区間木の考え方を元に配列をソーティングしながら絞り込みを
 * したが，バウンディングボックスで絞りこみをした方が効率的であった
 * ため，ボツにした．
 */
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include "base.h"
#include "mathematics.h"
#include "fem_struct.h"

/* 軸並行境界ボックス(axis-aligned bounding box) */
typedef struct{
	int index;  /* index of element array */
	VECT3D min;
	VECT3D max;
} GEOM_AABB;

static RC extract_geom_contact_surface_x(int *size1, GEOM_AABB **bb1,
                                         int *size2, GEOM_AABB **bb2);
static RC extract_geom_contact_surface_y(int *size1, GEOM_AABB **bb1,
                                         int *size2, GEOM_AABB **bb2);
static RC extract_geom_contact_surface_z(int *size1, GEOM_AABB **bb1,
                                         int *size2, GEOM_AABB **bb2);
static void sort_geom_aabb_array_x(int size, GEOM_AABB *bb);
static void sort_geom_aabb_array_y(int size, GEOM_AABB *bb);
static void sort_geom_aabb_array_z(int size, GEOM_AABB *bb);
static int compar_geom_aabb_min_x(const void *p1, const void *p2);
static int compar_geom_aabb_min_y(const void *p1, const void *p2);
static int compar_geom_aabb_min_z(const void *p1, const void *p2);
static RC cover_geom_aabb_array_x(int size, GEOM_AABB bb[]);
static RC cover_geom_aabb_array_y(int size, GEOM_AABB bb[]);
static RC cover_geom_aabb_array_z(int size, GEOM_AABB bb[]);
static int search_geom_aabb_min_x(int size, GEOM_AABB bb[], double min);
static int search_geom_aabb_min_y(int size, GEOM_AABB bb[], double min);
static int search_geom_aabb_min_z(int size, GEOM_AABB bb[], double min);
static int search_geom_aabb_max_x(int size, GEOM_AABB bb[], double max);
static int search_geom_aabb_max_y(int size, GEOM_AABB bb[], double max);
static int search_geom_aabb_max_z(int size, GEOM_AABB bb[], double max);

/*
static RC print_geom_aabb_array(int size, GEOM_AABB bb[]);
*/


/* 表面モデル surf1, surf2 の境界領域(Bounding Box)が交差している
 * ELEMENT_ARRAYのインデックス配列を返す関数(注意:labelではない)
 * 概して接触，内包，交差している可能性のある双方の表面要素群が抽出される．
 * モデル形状によるが，貫通している場合には正確に抽出できないこともある．
 * index_array{1,2} は予め表面要素数で確保しておくこと．
 */
RC
extract_geom_contact_surface(ELEMENT_ARRAY surf1, NODE_ARRAY node1,
                             ELEMENT_ARRAY surf2, NODE_ARRAY node2,
                             int *size1, int index_array1[],
                             int *size2, int index_array2[])
{
	int ii1;
	GEOM_AABB model1, model2; /* 各モデルのAABB */
	GEOM_AABB *bb1, *bb2; /* 要素ごとのAABB */
	GEOM_AABB *bb1_ptr, *bb2_ptr;
	double width[3];
	int bb1_size, bb2_size;
	int index[3] = {0, 1, 2};
	RC (*extract_func[3])(int *, GEOM_AABB **, int *, GEOM_AABB **)
	 = {extract_geom_contact_surface_x,
		extract_geom_contact_surface_y,
		extract_geom_contact_surface_z};

	if( (index_array1 == NULL) || (index_array2 == NULL)
	 || (size1 == NULL) || (size2 == NULL)){
		 return(ARG_ERROR_RC);
	}
	if( (surf1.size < 1) || (surf2.size < 1) ){
		*size1 = *size2 = 0;
		 return(NORMAL_RC);
	}

	RC_NULL_CHK( bb1 = mm_alloc(surf1.size*sizeof(GEOM_AABB)) );
	RC_NULL_CHK( bb2 = mm_alloc(surf2.size*sizeof(GEOM_AABB)) );
	/* 先頭アドレスを保存 */
	bb1_ptr = bb1;
	bb2_ptr = bb2;

	/* AABBを初期化 */
	model1.min.x = model1.min.y = model1.min.z =  DBL_MAX; 
	model1.max.x = model1.max.y = model1.max.z = -DBL_MAX; 
	model2.min.x = model2.min.y = model2.min.z =  DBL_MAX; 
	model2.max.x = model2.max.y = model2.max.z = -DBL_MAX; 

	/* bb1 配列に代入 */
	for(ii1=0, bb1_size=0; ii1<surf1.size; ii1++){
		if(surf1.array[ii1].label < 0) continue;
		bb1[ii1].index = ii1;
		RC_TRY( element_bbox(surf1.array[ii1], node1,
		                     &bb1[ii1].min, &bb1[ii1].max) );
		bb1_size++;
		model1.min.x = MIN2(model1.min.x, bb1[ii1].min.x);
		model1.min.y = MIN2(model1.min.y, bb1[ii1].min.y);
		model1.min.z = MIN2(model1.min.z, bb1[ii1].min.z);
		model1.max.x = MAX2(model1.max.x, bb1[ii1].max.x);
		model1.max.y = MAX2(model1.max.y, bb1[ii1].max.y);
		model1.max.z = MAX2(model1.max.z, bb1[ii1].max.z);
	}

	/* bb2 配列に代入 */
	for(ii1=0, bb2_size=0; ii1<surf2.size; ii1++){
		if(surf2.array[ii1].label < 0) continue;
		bb2[ii1].index = ii1;
		RC_TRY( element_bbox(surf2.array[ii1], node2,
		                     &bb2[ii1].min, &bb2[ii1].max) );
		bb2_size++;
		model2.min.x = MIN2(model2.min.x, bb2[ii1].min.x);
		model2.min.y = MIN2(model2.min.y, bb2[ii1].min.y);
		model2.min.z = MIN2(model2.min.z, bb2[ii1].min.z);
		model2.max.x = MAX2(model2.max.x, bb2[ii1].max.x);
		model2.max.y = MAX2(model2.max.y, bb2[ii1].max.y);
		model2.max.z = MAX2(model2.max.z, bb2[ii1].max.z);
	}

	/* モデル全域での交差領域の判定 */
	width[0] = MIN2(model1.max.x, model2.max.x)
	         - MAX2(model1.min.x, model2.min.x);
	width[1] = MIN2(model1.max.y, model2.max.y)
	         - MAX2(model1.min.y, model2.min.y);
	width[2] = MIN2(model1.max.z, model2.max.z)
	         - MAX2(model1.min.z, model2.min.z);
	/* width が0.0fでも接触している可能性はある */
	if( (width[0]+ABS_TOL < 0.0f) || (width[1]+ABS_TOL < 0.0f)
	 || (width[2]+ABS_TOL < 0.0f) ){
		/* 接触領域はない. 終了 */
		*size1 = *size2 = 0;
		RC_TRY( mm_free((void *)bb1_ptr) );
		RC_TRY( mm_free((void *)bb2_ptr) );
		return(NORMAL_RC);
	}

	/* 交差領域の少ない方向から順に計算する */
	if(width[index[0]] > width[index[1]]){
		RC_TRY( swap_int(&index[0], &index[1]) );
	}
	if(width[index[1]] > width[index[2]]){
		RC_TRY( swap_int(&index[1], &index[2]) );
	}
	if(width[index[0]] > width[index[1]]){
		RC_TRY( swap_int(&index[0], &index[1]) );
	}

	/* 選別効率の良い方向順に計算 */
	for(ii1=0; ii1<3; ii1++){
		RC_TRY( extract_func[index[ii1]](&bb1_size, &bb1,
		                                 &bb2_size, &bb2) );
	}

	*size1 = bb1_size;
	*size2 = bb2_size;

	for(ii1=0; ii1<bb1_size; ii1++){
		index_array1[ii1] = bb1[ii1].index;
	}
	RC_TRY( sort_int_array(*size1, index_array1) );

	for(ii1=0; ii1<bb2_size; ii1++){
		index_array2[ii1] = bb2[ii1].index;
	}
	RC_TRY( sort_int_array(*size2, index_array2) );

	RC_TRY( mm_free((void *)bb1_ptr) );
	RC_TRY( mm_free((void *)bb2_ptr) );

	return(NORMAL_RC);
}


#define FUNC_EXTRACT_GEOM_CONTACT_SURFACE(name, axis) \
static RC name(int *size1, GEOM_AABB **bb1, int *size2, GEOM_AABB **bb2) {\
	int index_min, index_max;\
	double section_min, section_max;\
	double min[2], max[2];\
\
	/* bbを min.axisの小さい順にソート */\
	sort_geom_aabb_array_##axis (*size1, *bb1);\
	/* max.axisを昇順にそれまでの最大値で塗りつぶす */\
	RC_TRY( cover_geom_aabb_array_##axis(*size1, *bb1) );\
\
	/* bbを min.axisの小さい順にソート */\
	sort_geom_aabb_array_##axis (*size2, *bb2);\
	/* max.axisを昇順にそれまでの最大値で塗りつぶす */\
	RC_TRY( cover_geom_aabb_array_##axis(*size2, *bb2) );\
\
	min[0] = (*bb1)[0].min.axis;\
	min[1] = (*bb2)[0].min.axis;\
	max[0] = (*bb1)[*size1-1].max.axis;\
	max[1] = (*bb2)[*size2-1].max.axis;\
	/* モデルの軸方向のBB幅が同じならこれ以降意味なし */ \
	if( (max[0] - REL_TOL*fabs(max[0]) - ABS_TOL < max[1]) \
	 && (max[0] + REL_TOL*fabs(max[0]) + ABS_TOL > max[1]) \
	 && (min[0] - REL_TOL*fabs(min[0]) - ABS_TOL < min[1]) \
	 && (min[0] + REL_TOL*fabs(min[0]) + ABS_TOL > min[1])) return(NORMAL_RC);\
\
	section_min = MAX2(min[0], min[1]);\
	section_max = MIN2(max[0], max[1]);\
\
	/* bb1[]の範囲を縮小(配列のメモリ確保領域は変わらない) */\
	RC_NEG_CHK( index_min = search_geom_aabb_max_##axis(*size1, *bb1,\
	                                                    section_min) );\
	RC_NEG_CHK( index_max = search_geom_aabb_min_##axis(*size1, *bb1,\
	                                                    section_max) );\
	RC_NEG_CHK( *size1 = index_max - index_min + 1 );\
	(*bb1) += index_min;\
\
	/* bb2[]の範囲を縮小(配列のメモリ確保領域は変わらない) */\
	RC_NEG_CHK( index_min = search_geom_aabb_max_##axis(*size2, *bb2,\
	                                                    section_min) );\
	RC_NEG_CHK( index_max = search_geom_aabb_min_##axis(*size2, *bb2,\
	                                                    section_max) );\
	RC_NEG_CHK( *size2 = index_max - index_min + 1 );\
	(*bb2) += index_min;\
\
	return(NORMAL_RC);\
}
FUNC_EXTRACT_GEOM_CONTACT_SURFACE(extract_geom_contact_surface_x, x)
FUNC_EXTRACT_GEOM_CONTACT_SURFACE(extract_geom_contact_surface_y, y)
FUNC_EXTRACT_GEOM_CONTACT_SURFACE(extract_geom_contact_surface_z, z)

static void
sort_geom_aabb_array_x(int size, GEOM_AABB bb[])
{
	qsort(bb, size, sizeof(GEOM_AABB), compar_geom_aabb_min_x);
}
static void
sort_geom_aabb_array_y(int size, GEOM_AABB bb[])
{
	qsort(bb, size, sizeof(GEOM_AABB), compar_geom_aabb_min_y);
}
static void
sort_geom_aabb_array_z(int size, GEOM_AABB bb[])
{
	qsort(bb, size, sizeof(GEOM_AABB), compar_geom_aabb_min_z);
}


static int
compar_geom_aabb_min_x(const void *p1, const void *p2)
{
	if( (((GEOM_AABB *)p1)->min.x) > (((GEOM_AABB *)p2)->min.x) ) return(1);
	if( (((GEOM_AABB *)p1)->min.x) < (((GEOM_AABB *)p2)->min.x) ) return(-1);
	return(0);
}
static int
compar_geom_aabb_min_y(const void *p1, const void *p2)
{
	if( (((GEOM_AABB *)p1)->min.y) > (((GEOM_AABB *)p2)->min.y) ) return(1);
	if( (((GEOM_AABB *)p1)->min.y) < (((GEOM_AABB *)p2)->min.y) ) return(-1);
	return(0);
}
static int
compar_geom_aabb_min_z(const void *p1, const void *p2)
{
	if( (((GEOM_AABB *)p1)->min.z) > (((GEOM_AABB *)p2)->min.z) ) return(1);
	if( (((GEOM_AABB *)p1)->min.z) < (((GEOM_AABB *)p2)->min.z) ) return(-1);
	return(0);
}


#define FUNC_COVER_GEOM_AABB_ARRAY(name, axis) \
static RC name(int size, GEOM_AABB bb[]) {\
	int ii1;\
	if(size <= 0) return(ARG_ERROR_RC);\
	if(bb == NULL) return(ARG_ERROR_RC);\
	for(ii1=1; ii1<size; ii1++){\
		if(bb[ii1-1].max.axis > bb[ii1].max.axis){\
			bb[ii1].max.axis = bb[ii1-1].max.axis;\
		}\
	}\
	return(NORMAL_RC);\
}
FUNC_COVER_GEOM_AABB_ARRAY(cover_geom_aabb_array_x, x)
FUNC_COVER_GEOM_AABB_ARRAY(cover_geom_aabb_array_y, y)
FUNC_COVER_GEOM_AABB_ARRAY(cover_geom_aabb_array_z, z)


/* bb[]から最小値が min 以下の要素を見つけ，そのindexを返す */
#define FUNC_SEARCH_GEOM_AABB_MIN(name, axis) \
static int name(int size, GEOM_AABB bb[], double min){ \
	int found_index = -1;\
	int low, high, middle;\
	if(size <= 0) return(-1);\
	if(bb == NULL) return(-1);\
	low = 0;\
	high = size -1;\
\
	while(low <= high){\
		middle = (low + high)/2;\
		if(bb[middle].min.axis <= min){\
			if(middle+1 >= size){\
				found_index = middle;\
				break;\
			}else if(bb[middle+1].min.axis > min){\
				found_index = middle;\
				break;\
			}else{\
				low = middle + 1;\
			}\
		}else{\
			high = middle - 1;\
		}\
	}\
	return(found_index);\
}
FUNC_SEARCH_GEOM_AABB_MIN(search_geom_aabb_min_x, x)
FUNC_SEARCH_GEOM_AABB_MIN(search_geom_aabb_min_y, y)
FUNC_SEARCH_GEOM_AABB_MIN(search_geom_aabb_min_z, z)

/* bb[]から最大値が max 以上の要素を見つけ，そのindexを返す */
#define FUNC_SEARCH_GEOM_AABB_MAX(name, axis) \
static int name(int size, GEOM_AABB bb[], double max) {\
	int found_index = -1;\
	int low, high, middle;\
	if(size <= 0) return(-1);\
	if(bb == NULL) return(-1);\
	low = 0;\
	high = size -1;\
\
	while(low <= high){\
		middle = (low + high)/2;\
		if(bb[middle].max.axis >= max){\
			if(middle-1 < 0){\
				found_index = middle;\
				break;\
			}else if(bb[middle-1].max.axis < max){\
				found_index = middle;\
				break;\
			}else{\
				high = middle - 1;\
			}\
		}else{\
			low = middle + 1;\
		}\
	}\
	return(found_index);\
}
FUNC_SEARCH_GEOM_AABB_MAX(search_geom_aabb_max_x, x)
FUNC_SEARCH_GEOM_AABB_MAX(search_geom_aabb_max_y, y)
FUNC_SEARCH_GEOM_AABB_MAX(search_geom_aabb_max_z, z)


