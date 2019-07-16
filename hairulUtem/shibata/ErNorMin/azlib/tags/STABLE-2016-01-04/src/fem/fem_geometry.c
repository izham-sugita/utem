/*********************************************************************
 * fem_geometry.c
 *
 * Copyright (C) 2004 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: fem_geometry.c 920 2009-03-05 05:31:26Z iwai $ */

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

static RC clip_geometry(int *size, GEOM_AABB bb[], GEOM_AABB clip,
                        GEOM_AABB *new_bb);
/*
static RC print_geom_aabb_array(int size, GEOM_AABB bb[]);
*/

static RC sort_segment_array(int size, GEOM_SEGMENT segment[]);
static int compar_segment_min(const void *min, const void *max);
static RC cover_segment_array(int size, GEOM_SEGMENT segment[]);
static int search_segment_min(int size, GEOM_SEGMENT segment[], double min);


RC
element_bbox (ELEMENT element, NODE_ARRAY node, VECT3D *min, VECT3D *max)
{
	int ii1;
	int index;

	if(element.label < 0) return(ARG_ERROR_RC);
	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	min->x = min->y = min->z = DBL_MAX;
	max->x = max->y = max->z = -DBL_MAX;

	for(ii1=0; ii1<element.node_num; ii1++){
		index = node_renum_index(node, element.node[ii1]);
		if(index < 0) return(SEEK_ERROR_RC);
		max->x = MAX2(node.array[index].p.x, max->x);
		max->y = MAX2(node.array[index].p.y, max->y);
		max->z = MAX2(node.array[index].p.z, max->z);
		min->x = MIN2(node.array[index].p.x, min->x);
		min->y = MIN2(node.array[index].p.y, min->y);
		min->z = MIN2(node.array[index].p.z, min->z);
	}

	return(NORMAL_RC);
}


/* doc/Memo/Closest_Point_Search/Contact_Geometry_Search.ppt 参照 */
RC
extract_geom_contact_surface(ELEMENT_ARRAY surf1, NODE_ARRAY node1,
                             ELEMENT_ARRAY surf2, NODE_ARRAY node2,
                             int *size1, int index_array1[],
                             int *size2, int index_array2[])
{
	int ii1, itmp;
	GEOM_AABB model1, model2; /* 各モデルのAABB */
	GEOM_AABB *bb1, *bb2; /* 要素ごとのAABB */
	GEOM_AABB clip;
	double width[3];
	int bb1_size, bb2_size;
	int index[] = {0, 1, 2};
	int cliped_flag[] = {0, 0, 0};

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

	for(ii1=0; ii1<3; ii1++){ /* XYZの3軸 */
		/* モデルの交差領域 == 切り抜き */
		if(cliped_flag[0] == 0){ /* 同じ方向を選別しない */
			clip.min.x = MAX2(model1.min.x, model2.min.x);
			clip.max.x = MIN2(model1.max.x, model2.max.x);
			width[0] = clip.max.x - clip.min.x;
		}
		if(cliped_flag[1] == 0){
			clip.min.y = MAX2(model1.min.y, model2.min.y);
			clip.max.y = MIN2(model1.max.y, model2.max.y);
			width[1] = clip.max.y - clip.min.y;
		}
		if(cliped_flag[2] == 0){
			clip.min.z = MAX2(model1.min.z, model2.min.z);
			clip.max.z = MIN2(model1.max.z, model2.max.z);
			width[2] = clip.max.z - clip.min.z;
		}

		/* width が0.0fでも接触している可能性はある(<=ではない) */
		if( (width[0]+ABS_TOL < 0.0f) || (width[1]+ABS_TOL < 0.0f)
		 || (width[2]+ABS_TOL < 0.0f) ){
			/* 接触領域はない. 終了 */
			*size1 = *size2 = 0;
			RC_TRY( mm_free((void *)bb1) );
			RC_TRY( mm_free((void *)bb2) );
			return(NORMAL_RC);
		}

		/* int index[] = {0, 1, 2};  0: x  1: y  2: z */
		/* 交差領域の少ない方向から順に計算する       */
		if(ii1 == 0){
			if(width[index[0]] > width[index[1]]){
				SWAP(index[0], index[1], itmp);
			}
			if(width[index[1]] > width[index[2]]){
				SWAP(index[1], index[2], itmp);
			}
			if(width[index[0]] > width[index[1]]){
				SWAP(index[0], index[1], itmp);
			}
		}else if(ii1 == 1){
			if(width[index[1]] > width[index[2]]){
				SWAP(index[1], index[2], itmp);
			}
		}
		cliped_flag[index[ii1]] = 1; /* 交差量が最小の所は範囲を変えない */

		/* clip範囲内の要素に絞り込み，配列を上方に詰める */
		RC_TRY( clip_geometry(&bb1_size, bb1, clip, &model1) );
		RC_TRY( clip_geometry(&bb2_size, bb2, clip, &model2) );
		if( (bb1_size == 0) || (bb2_size == 0) ){
			*size1 = *size2 = 0;
			RC_TRY( mm_free((void *)bb1) );
			RC_TRY( mm_free((void *)bb2) );
			return(NORMAL_RC);
		}
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

	RC_TRY( mm_free((void *)bb1) );
	RC_TRY( mm_free((void *)bb2) );

	return(NORMAL_RC);
}


static RC
clip_geometry(int *size, GEOM_AABB bb[], GEOM_AABB clip, GEOM_AABB *new_bb)
{
	int ii1;
	int pos;

	if( (size == NULL) || (bb == NULL) || (new_bb == NULL) )
		return(ARG_ERROR_RC);

	new_bb->min.x = new_bb->min.y = new_bb->min.z =  DBL_MAX; 
	new_bb->max.x = new_bb->max.y = new_bb->max.z = -DBL_MAX; 

	for(ii1=0, pos=0; ii1<*size; ii1++){
		if( (bb[ii1].min.x <= clip.max.x+REL_TOL*fabs(clip.min.x)+ABS_TOL)
		 && (bb[ii1].max.x >= clip.min.x-REL_TOL*fabs(clip.min.x)-ABS_TOL)
		 && (bb[ii1].min.y <= clip.max.y+REL_TOL*fabs(clip.min.y)+ABS_TOL)
		 && (bb[ii1].max.y >= clip.min.y-REL_TOL*fabs(clip.min.y)-ABS_TOL)
		 && (bb[ii1].min.z <= clip.max.z+REL_TOL*fabs(clip.min.z)+ABS_TOL)
		 && (bb[ii1].max.z >= clip.min.z-REL_TOL*fabs(clip.min.z)-ABS_TOL)){
			/* clip範囲に在るので配列の上方に詰める */
			if(pos != ii1) bb[pos] = bb[ii1];
			/* ついでに clip 範囲内の bb[] のAABBを求める */
			new_bb->min.x = MIN2(new_bb->min.x, bb[ii1].min.x);
			new_bb->min.y = MIN2(new_bb->min.y, bb[ii1].min.y);
			new_bb->min.z = MIN2(new_bb->min.z, bb[ii1].min.z);
			new_bb->max.x = MAX2(new_bb->max.x, bb[ii1].max.x);
			new_bb->max.y = MAX2(new_bb->max.y, bb[ii1].max.y);
			new_bb->max.z = MAX2(new_bb->max.z, bb[ii1].max.z);
			pos++;
		}
	}
	*size = pos;

	return(NORMAL_RC);
}


#if 0
static RC
print_geom_aabb_array(int size, GEOM_AABB bb[])
{
	int ii1;

	RC_TRY( log_printf(5, "num     index      min.x    max.x     min.y    "
	                      "max.y     min.z    max.z\n") );
	for(ii1=0; ii1<size; ii1++){
		RC_TRY( log_printf(5, "[%5d] [%5d] % 8.3f % 8.3f, % 8.3f % 8.3f,"
		                      " % 8.3f % 8.3f\n", ii1, bb[ii1].index,
		                   bb[ii1].min.x, bb[ii1].max.x,
		                   bb[ii1].min.y, bb[ii1].max.y,
		                   bb[ii1].min.z, bb[ii1].max.z) );
	}

	return(NORMAL_RC);
}
#endif


RC
make_geom_search_array(ELEMENT_ARRAY elem, NODE_ARRAY node,
                        int *size, GEOM_SEGMENT *segment[])
{
	int ii1, ii2;
	VECT3D min, max;

	if(!(*segment)) return(ARG_ERROR_RC);

	*size= count_valid_element(elem);
	if(*size < 1) return(NORMAL_RC);
	RC_NULL_CHK( segment[0] = mm_alloc(*size*sizeof(GEOM_SEGMENT)) );
	RC_NULL_CHK( segment[1] = mm_alloc(*size*sizeof(GEOM_SEGMENT)) );
	RC_NULL_CHK( segment[2] = mm_alloc(*size*sizeof(GEOM_SEGMENT)) );

	for(ii1=0, ii2=0; ii1<elem.size; ii1++){
		if(elem.array[ii1].label < 0) continue;
		/* Bounding Box */
		RC_TRY( element_bbox(elem.array[ii1], node, &min, &max) );
		segment[0][ii2].index = ii1;
		segment[0][ii2].min = min.x;
		segment[0][ii2].max = max.x;

		segment[1][ii2].index = ii1;
		segment[1][ii2].min = min.y;
		segment[1][ii2].max = max.y;

		segment[2][ii2].index = ii1;
		segment[2][ii2].min = min.z;
		segment[2][ii2].max = max.z;

		ii2++;
	}

	/* 各elementの大きさの配列を最小値で sort */
	sort_segment_array(*size, segment[0]);
	sort_segment_array(*size, segment[1]);
	sort_segment_array(*size, segment[2]);

	/* 最小値で sortした配列の最大値も昇順になるように塗りつぶし */
	RC_TRY( cover_segment_array(*size, segment[0]) );
	RC_TRY( cover_segment_array(*size, segment[1]) );
	RC_TRY( cover_segment_array(*size, segment[2]) );

	return(NORMAL_RC);
}


/* make_geom_search_array()の引数が違うバージョン        */
/* elem全体ではなく index_array[?]番目をターゲットにする */
RC
make_geom_search_array_2(ELEMENT_ARRAY elem, NODE_ARRAY node,
                         int index_size, int *index_array,
                         int *size, GEOM_SEGMENT *segment[])
{
	int ii1, ii2;
	VECT3D min, max;

	if(!(index_array)) return(ARG_ERROR_RC);

	*size= index_size;
	if(*size < 1) return(NORMAL_RC);
	RC_NULL_CHK( segment[0] = mm_alloc(*size*sizeof(GEOM_SEGMENT)) );
	RC_NULL_CHK( segment[1] = mm_alloc(*size*sizeof(GEOM_SEGMENT)) );
	RC_NULL_CHK( segment[2] = mm_alloc(*size*sizeof(GEOM_SEGMENT)) );

	for(ii1=0, ii2=0; ii1<index_size; ii1++){
		if(elem.array[index_array[ii1]].label < 0) continue;
		/* Bounding Box */
		RC_TRY( element_bbox(elem.array[index_array[ii1]], node, &min, &max) );
		segment[0][ii2].index = index_array[ii1];
		segment[0][ii2].min = min.x;
		segment[0][ii2].max = max.x;

		segment[1][ii2].index = index_array[ii1];
		segment[1][ii2].min = min.y;
		segment[1][ii2].max = max.y;

		segment[2][ii2].index = index_array[ii1];
		segment[2][ii2].min = min.z;
		segment[2][ii2].max = max.z;

		ii2++;
	}

	/* 各elementの大きさの配列を最小値で sort */
	sort_segment_array(*size, segment[0]);
	sort_segment_array(*size, segment[1]);
	sort_segment_array(*size, segment[2]);

	/* 最小値で sortした配列の最大値も昇順になるように塗りつぶし */
	RC_TRY( cover_segment_array(*size, segment[0]) );
	RC_TRY( cover_segment_array(*size, segment[1]) );
	RC_TRY( cover_segment_array(*size, segment[2]) );

	return(NORMAL_RC);
}


static RC
sort_segment_array(int size, GEOM_SEGMENT segment[])
{
	qsort(segment, size, sizeof(GEOM_SEGMENT), compar_segment_min);
	return(NORMAL_RC);
}


static int
compar_segment_min(const void *p1, const void *p2)
{
	if( (((GEOM_SEGMENT *)p1)->min) > (((GEOM_SEGMENT *)p2)->min) ) return(1);
	if( (((GEOM_SEGMENT *)p1)->min) < (((GEOM_SEGMENT *)p2)->min) ) return(-1);
	return(0);
}


static RC
cover_segment_array(int size, GEOM_SEGMENT segment[])
{
	int ii1;

	if(size <= 0) return(ARG_ERROR_RC);
	if(segment == NULL) return(ARG_ERROR_RC);

	for(ii1=1; ii1<size; ii1++){
		if(segment[ii1-1].max > segment[ii1].max){
			segment[ii1].max = segment[ii1-1].max;
		}
	}

	return(NORMAL_RC);
}


/* make_geom_search_array()で作られた配列から探査候補配列を作る */
/* segmentの中から min 以上 max 以下の elementのindex配列を作る */
RC
geom_search_candidate(int size, GEOM_SEGMENT *segment[],
                      VECT3D p, double range, int *cand_size, int cand[])

{
	int ii1;
	double tmin[3], tmax[3];
	int index[3];
	int break_num;

	if(range < 0.0) return(ARG_ERROR_RC);

	tmin[0] = p.x - range; tmin[1] = p.y - range; tmin[2] = p.z - range;
	tmax[0] = p.x + range; tmax[1] = p.y + range; tmax[2] = p.z + range;

	for(ii1=0; ii1<3; ii1++){
		index[ii1] = search_segment_min(size, segment[ii1], tmax[ii1]);
	}

	/* segment[]から最大値を降順に探査しtmin以上である要素を探す */
	/* x,y,zの内，先に終わったものを選ぶ                         */
	break_num = -1;
	*cand_size = 0;
	while(1){
		for(ii1=0; ii1<3; ii1++){
			if(index[ii1] < 0){
				break_num = ii1;
				break;
			}
			if(segment[ii1][index[ii1]].max < tmin[ii1]){
				break_num = ii1;
				break;
			}
			(index[ii1])--;
		}
		if(break_num >= 0) break;
		(*cand_size)++;
	}

	index[break_num] = search_segment_min(size, segment[break_num],
			tmax[break_num]);
	for(ii1=0; ii1<(*cand_size); ii1++){
		cand[ii1] = segment[break_num][index[break_num]].index;
		(index[break_num])--;
	}

	return(NORMAL_RC);
}


/* segment[]から最小値が min 以下の要素を見つけ，そのindexを返す */
static int
search_segment_min(int size, GEOM_SEGMENT segment[], double min)
{
	int found_index = -1;
	int low, high, middle;

	if(size <= 0) return(-1);
	if(segment == NULL) return(-1);

	low = 0;
	high = size -1;

	while(low <= high){
		middle = (low + high)/2;

		if(segment[middle].min <= min){
			if(middle+1 >= size){
				found_index = middle;
				break;
			}else if(segment[middle+1].min > min){
				found_index = middle;
				break;
			}else{
				low = middle + 1;
			}
		}else{
			high = middle - 1;
		}
	}

	return(found_index);
}


RC
print_geom_segment_array(FILE *fp, int size, GEOM_SEGMENT segment[])
{
	int ii1;

	if(segment == NULL) return(ARG_ERROR_RC);
	if(size < 0) return(ARG_ERROR_RC);

	for(ii1=0; ii1<size; ii1++){
		fprintf(fp, "[%5d] index = %5d  min = %11.3e  max = %11.3e\n",
		        ii1, segment[ii1].index, segment[ii1].min, segment[ii1].max);
	}

	return(NORMAL_RC);
}


