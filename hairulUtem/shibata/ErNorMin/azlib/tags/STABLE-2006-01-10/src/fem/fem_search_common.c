/*********************************************************************
 * fem_search_common.c
 *
 * Copyright (C) 2005 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: fem_search_common.c 415 2005-07-15 06:04:24Z sasaoka $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rc.h"
#include "math_utl.h"
#include "fem_struct.h"
#include "memory_manager.h"


/*
 * labelの情報を手がかりに，モデル間の共通の節点を捜す．
 * 領域分割した2モデル間の共通節点(通常，境界にある)を探すのが主用途．
 * node1, node2 に不要な節点が含まれないことが前提．
 * 節点の幾何情報まで確認しない
 */
RC
search_common_node_label (NODE_ARRAY node1, NODE_ARRAY node2,
                          NODE_ARRAY *common)
{
	int ii1, ii2, ii3;
	int *label_array1, *label_array2;

	if( !((node1.size > 0)&&(node2.size > 0)) ) return(ARG_ERROR_RC);

	label_array1 = (int *)mm_alloc(node1.size*sizeof(int));
	if(label_array1 == NULL) return(ALLOC_ERROR_RC);

	label_array2 = (int *)mm_alloc(node2.size*sizeof(int));
	if(label_array2 == NULL) return(ALLOC_ERROR_RC);

	for(ii1=0; ii1<node1.size; ii1++){
		label_array1[ii1] = node1.array[ii1].label;
	}
	RC_TRY( sort_int_array(node1.size, label_array1) );

	for(ii1=0; ii1<node2.size; ii1++){
		label_array2[ii1] = node2.array[ii1].label;
	}
	RC_TRY( sort_int_array(node2.size, label_array2) );

	for(ii1=0, ii2=0, ii3=0; ii1<node1.size; ii1++){
		for(; ii2<node2.size; ii2++){
			if(label_array1[ii1] == label_array2[ii2]){
				/* 共通節点の情報をlabel_array1に上書き */
				label_array1[ii3] = label_array1[ii1];
				ii3++;
				break;
			}
			if(label_array1[ii1] < label_array2[ii2]){
				break;
			}
		}
		if(ii2 >= node2.size-1) break;
	}

	RC_TRY( allocate_node_array(ii3, common) );
	for(ii1=0; ii1<ii3; ii1++){
		int index = search_node_label(node1, label_array1[ii1]);
		if(index < 0){
			RC_TRY( mm_free(common) );
			return(NEGATIVE_ERROR_RC);
		}
		common->array[ii1] = node1.array[index];
	}
	common->size = ii3;
	common->source = node1.source;
	RC_TRY( clean_node_array(common) );

	RC_TRY( mm_free(label_array1) );
	RC_TRY( mm_free(label_array2) );

	return(NORMAL_RC);
}


/*
 * search_common_node_label()と同じ機能の関数．
 * 領域分割した2モデル間の共通節点を探すのが目的．
 * 要素表面の節点のみを検索の対象とするため，場合によっては上記関数より高速
 * になる．elem1, elem2はそれぞれnode1, node2を持つ要素．表面要素ではない．
 */
RC
search_common_node_label_element (NODE_ARRAY node1, NODE_ARRAY node2,
                                  ELEMENT_ARRAY elem1, ELEMENT_ARRAY elem2,
                                  NODE_ARRAY *common)
{
	int ii1, ii2, ii3;
	int *label_array1, *label_array2;
	int label_size1, label_size2;
	int index;
    ELEMENT_ARRAY surf1, surf2;


	if( !((node1.size > 0)&&(node2.size > 0)) ) return(ARG_ERROR_RC);
	if( !((elem1.size > 0)&&(elem2.size > 0)) ) return(ARG_ERROR_RC);

	RC_TRY( extract_surface(elem1, &surf1) );
	RC_TRY( extract_surface(elem2, &surf2) );

	/* quad2 : 8 node */
	label_array1 = (int *)mm_alloc(surf1.size*8*sizeof(int));
	if(label_array1 == NULL) return(ALLOC_ERROR_RC);

	label_array2 = (int *)mm_alloc(surf2.size*8*sizeof(int));
	if(label_array2 == NULL) return(ALLOC_ERROR_RC);

	/* 表面の全節点の重複しないリストを作る */
	for(ii1=0, ii3=0; ii1<surf1.size; ii1++){
		if(surf1.array[ii1].label < 0) continue;
		for(ii2=0; ii2<surf1.array[ii1].node_num; ii2++){
			RC_NEG_CHK(index = search_node_label(node1,
			                                     surf1.array[ii1].node[ii2]) );
			label_array1[ii3++] = node1.array[index].label;
		}
	}
	label_size1 = ii3;
	RC_TRY( sort_int_array(label_size1, label_array1) );
	RC_TRY( uniq_int_array(&label_size1, &label_array1) );
	RC_TRY( free_element_array(&surf1) );

	for(ii1=0, ii3=0; ii1<surf2.size; ii1++){
		if(surf2.array[ii1].label < 0) continue;
		for(ii2=0; ii2<surf2.array[ii1].node_num; ii2++){
			RC_NEG_CHK(index = search_node_label(node2,
			                                     surf2.array[ii1].node[ii2]) );
			label_array2[ii3++] = node2.array[index].label;
		}
	}
	label_size2 = ii3;
	RC_TRY( sort_int_array(label_size2, label_array2) );
	RC_TRY( uniq_int_array(&label_size2, &label_array2) );
	RC_TRY( free_element_array(&surf2) );

	/* label_array1とlabel_array2の同じ節点を見付け，label_array1に上書き */
	for(ii1=0, ii2=0, ii3=0; ii1<label_size1; ii1++){
		for(; ii2<label_size2; ii2++){
			if(label_array1[ii1] == label_array2[ii2]){
				/* 共通節点の情報をlabel_array1に上書き */
				label_array1[ii3] = label_array1[ii1];
				ii3++;
				break;
			}
			if(label_array1[ii1] < label_array2[ii2]){
				break;
			}
		}
		if(ii2 >= label_size2-1) break;
	}

	RC_TRY( allocate_node_array(ii3, common) );
	for(ii1=0; ii1<ii3; ii1++){
		int index = search_node_label(node1, label_array1[ii1]);
		if(index < 0){
			RC_TRY( mm_free(common) );
			return(NEGATIVE_ERROR_RC);
		}
		common->array[ii1] = node1.array[index];
	}
	common->size = ii3;
	common->source = node1.source;
	RC_TRY( clean_node_array(common) );

	RC_TRY( mm_free(label_array1) );
	RC_TRY( mm_free(label_array2) );

	return(NORMAL_RC);
}


