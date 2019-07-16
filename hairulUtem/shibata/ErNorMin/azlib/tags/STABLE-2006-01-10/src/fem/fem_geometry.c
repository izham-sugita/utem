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

/* $Id: fem_geometry.c 415 2005-07-15 06:04:24Z sasaoka $ */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include "rc.h"
#include "memory_manager.h"
#include "math_utl.h"
#include "fem_struct.h"

static void sort_segment_array(GEOM_SEGMENT_ARRAY *segment);
static int compar_segment_min(const void *min, const void *max);
static RC cover_segment_array(GEOM_SEGMENT_ARRAY *segment);
static int search_segment_min(GEOM_SEGMENT_ARRAY segment, double min);

RC
element_bbox (ELEMENT element, NODE_ARRAY node, VECT3D *min, VECT3D *max)
{
	int ii1;
	int index;

	if(element.label < 0) return(ARG_ERROR_RC);

	min->x = min->y = min->z = DBL_MAX;
	max->x = max->y = max->z = -DBL_MAX;

	for(ii1=0; ii1<element.node_num; ii1++){
		index = search_node_label(node, element.node[ii1]);
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


RC
make_geom_search_array (ELEMENT_ARRAY elem, NODE_ARRAY node,
                        GEOM_SEGMENT_ARRAY (*segment)[3])
{
	int ii1, ii2;
	int num;
	VECT3D min, max;

	if(!(*segment)) return(ARG_ERROR_RC);

	num = count_valid_element(elem);
	RC_TRY( allocate_geom_segment_array(num, &(*segment)[0]) );
	RC_TRY( allocate_geom_segment_array(num, &(*segment)[1]) );
	RC_TRY( allocate_geom_segment_array(num, &(*segment)[2]) );

	for(ii1=0, ii2=0; ii1<elem.size; ii1++){
		if(elem.array[ii1].label < 0) continue;
		/* Bounding Box */
		RC_TRY( element_bbox(elem.array[ii1], node, &min, &max) );

		(*segment)[0].array[ii2].min = min.x;
		(*segment)[0].array[ii2].max = max.x;
		(*segment)[0].array[ii2].label = elem.array[ii1].label;
		(*segment)[0].size++;

		(*segment)[1].array[ii2].min = min.y;
		(*segment)[1].array[ii2].max = max.y;
		(*segment)[1].array[ii2].label = elem.array[ii1].label;
		(*segment)[1].size++;

		(*segment)[2].array[ii2].min = min.z;
		(*segment)[2].array[ii2].max = max.z;
		(*segment)[2].array[ii2].label = elem.array[ii1].label;
		(*segment)[2].size++;

		ii2++;
	}

	/* 各elementの大きさの配列を最小値で sort */
	sort_segment_array(&(*segment)[0]);
	sort_segment_array(&(*segment)[1]);
	sort_segment_array(&(*segment)[2]);

	/* 最小値で sortした配列の最大値も昇順になるように塗りつぶし */
	RC_TRY( cover_segment_array(&(*segment)[0]) );
	RC_TRY( cover_segment_array(&(*segment)[1]) );
	RC_TRY( cover_segment_array(&(*segment)[2]) );

	return(NORMAL_RC);
}


static void
sort_segment_array (GEOM_SEGMENT_ARRAY *segment)
{
	qsort(segment->array, segment->size,
	      sizeof(GEOM_SEGMENT), compar_segment_min);
	segment->sort_flag = ARRAY_SORTED;
}


static int
compar_segment_min (const void *p1, const void *p2)
{
	if( nearly_eq( ((GEOM_SEGMENT *)p1)->min,
	               ((GEOM_SEGMENT *)p2)->min) ) return(0);

	return( (((GEOM_SEGMENT *)p1)->min > ((GEOM_SEGMENT *)p2)->min) ? 1 : -1 );
}


static RC
cover_segment_array (GEOM_SEGMENT_ARRAY *segment)
{
	int ii1;

	if(segment->size <= 0) return(ARG_ERROR_RC);
	if(segment->array == NULL) return(ARG_ERROR_RC);

	for(ii1=1; ii1<segment->size; ii1++){
		if(segment->array[ii1].label < 0) continue;
		if(segment->array[ii1-1].max > segment->array[ii1].max){
			segment->array[ii1].max = segment->array[ii1-1].max;
		}
	}

	return(NORMAL_RC);
}


/* make_geom_search_array()で作られた配列から探査候補配列を作る */
/* segmentの中から min 以上 max 以下の elementのlabel配列を作る */
RC
geom_search_candidate (const GEOM_SEGMENT_ARRAY segment[],
                       VECT3D min, VECT3D max, int *size, int **array)
{
	int ii1, ii2;
	double tmin[3], tmax[3];
	int index[3];
	int tsize[3];
	int break_num;

	if(min.x > max.x) return(ARG_ERROR_RC);
	if(min.y > max.y) return(ARG_ERROR_RC);
	if(min.z > max.z) return(ARG_ERROR_RC);

	tmin[0] = min.x; tmin[1] = min.y; tmin[2] = min.z;
	tmax[0] = max.x; tmax[1] = max.y; tmax[2] = max.z;

	for(ii1=0; ii1<3; ii1++){
		index[ii1] = search_segment_min(segment[ii1], tmax[ii1]);
		 /* 一つでも indexがなければ，探査領域に候補はない */
		if(index[ii1] < 0){
			*size = 0;
			*array = NULL;
			return(NORMAL_RC);
		}
		tsize[ii1] = 0;
	}

	/* segment->arrayから最大値を降順に探査しtmin以上である要素を探す */
	/* x,y,zの内，先に終わったものを選ぶ                              */
	break_num = -1;
	ii2 = 0;
	while(1){
		for(ii1=0; ii1<3; ii1++){
			int tindex = index[ii1] - ii2;
			if(tindex < 0){
				break_num = ii1;
				break;
			}
			if(segment[ii1].array[tindex].label < 0) continue;
			if(segment[ii1].array[tindex].max > tmin[ii1]){
				tsize[ii1]++;
			}else{
				break_num = ii1;
				break;
			}
		}
		ii2++;
		if(break_num != -1) break;
	}

	if(break_num == -1) return(UNKNOWN_ERROR_RC);

	*size = tsize[break_num];
	if(*size > 0){
		*array = (int *)mm_alloc(sizeof(int)*(*size));
		if(*array == NULL) return(ALLOC_ERROR_RC);
	}else{
		return(NORMAL_RC);
	}

	for(ii1=0, ii2=0; ii2<(*size); ii1++){
		int tindex = index[break_num] - ii1;
		int label = segment[break_num].array[tindex].label;
		if(label < 0) continue;
		(*array)[ii2++] = label;
	}

	return(NORMAL_RC);
}


/* segment->arrayから最小値が min 以下の要素を見つけ，そのindexを返す */
static int
search_segment_min (GEOM_SEGMENT_ARRAY segment, double min)
{
	int found_index = -1;
	int low, high, middle;

	if(segment.size <= 0) return(-1);
	if(segment.array == NULL) return(-1);

	low = 0;
	high = segment.size -1;

	while(low <= high){
		middle = (low + high)/2;

		if(segment.array[middle].min <= min){
			if(middle+1 >= segment.size){
				found_index = middle;
				break;
			}else if(segment.array[middle+1].min > min){
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


