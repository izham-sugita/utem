/*********************************************************************
 * fem_renumber.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Thanks to following contributors.
 *  <Tomoyuki TSUBATA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: fem_renumber.c 842 2006-12-04 06:21:20Z sasaoka $ */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include "rc.h"
#include "memory_manager.h"
#include "math_utl.h"
#include "fem_struct.h"


/* 節点の「枝」情報 */
typedef struct {
	int size;
	int alloc_size;
	int *br;
	int weight;
} NODAL_BRANCH;

static int count_queue(int index, const NODAL_BRANCH *branch,
                       int queue_size, const int *queue);
static RC remove_queue(int rm_index, int *queue_size, int *queue,
                       NODAL_BRANCH *branch, int *next_index, int type_flag);
static RC add_queue(int new_index, int *queue_size, int *queue);
static RC init_branch(NODE_ARRAY node, ELEMENT_ARRAY element,
                      RIGID_ELEMENT_ARRAY rigid, NODAL_BRANCH *branch);
static RC remove_branch(NODAL_BRANCH *branch_ptr, int rm_br);
static RC add_branch(NODAL_BRANCH *branch_ptr, int new_br);


/* 節点番号の付けかえ（ダミー版） */
RC
dummy_renumber (NODE_ARRAY *node)
{
	int ii1;
	int count = 0;

	for(ii1=0; ii1<node->size; ii1++){
		if(node->array[ii1].label < 0){
			node->array[ii1].renum = -1;
		}else{
			node->array[ii1].renum = count;
			count++;
		}
	}
	node->renum_flag = RENUMBERED;

	return(NORMAL_RC);
}


/* type_flag: 0 高速，低性能版 */
/*            1 低速，高性能版 */
RC
fem_renumber0 (NODE_ARRAY *node, ELEMENT_ARRAY element,
               RIGID_ELEMENT_ARRAY rigid, int type_flag)
{
	NODAL_BRANCH *branch;
	int *queue;
	int queue_size;
	int ii1;
	int min_size;
	int next_index;
	int next_label;

	/* branch, queue の確保 */
	branch = (NODAL_BRANCH *)mm_alloc(node->size * sizeof(NODAL_BRANCH));
	if(branch == NULL) return(ALLOC_ERROR_RC);

	queue = (int *)mm_alloc(node->size * sizeof(int));
	if(queue == NULL) return(ALLOC_ERROR_RC);

	/* renum をリセット */
	for(ii1=0; ii1<node->size; ii1++){
		node->array[ii1].renum = -1;
	}
	node->renum_flag = UNRENUMBERED;

	/* branch を初期化し、最も size の小さい要素を next_index に */
	RC_TRY( init_branch(*node, element, rigid, branch) );
	min_size = INT_MAX;
	next_index = -1;
	for(ii1=0; ii1<(node->size); ii1++){
		if(node->array[ii1].label < 0)continue;
		if(branch[ii1].size <= 0) continue;
		if(branch[ii1].size < min_size){
			min_size = branch[ii1].size;
			next_index = ii1;
		}
	}
	if(next_index < 0) return(ARG_ERROR_RC);
	
	/* next_index を queue に登録 */
	queue_size = 0;
	RC_TRY( add_queue(next_index, &queue_size, queue) );

	for(ii1=0; ii1<(node->size); ii1++){
		RC_TRY( progress_printf(5, " [%6d/%6d]", ii1, node->size) );
		node->array[next_index].renum = ii1;
		RC_TRY( remove_queue(next_index, &queue_size, queue,
		                     branch, &next_index, type_flag) );
		if(next_index < 0) break;
	}
	next_label = ii1;
	RC_TRY( progress_printf(5, "                "));

	/* チェック */
	for(ii1=0; ii1<node->size; ii1++){
		if(node->array[ii1].label < 0)continue;
		if(node->array[ii1].renum < 0){
			RC_TRY( warning_printf(3002, node->array[ii1].label) );
			node->array[ii1].renum = next_label;
			next_label++;
		}
	}

	for(ii1=0; ii1<node->size; ii1++){
		RC_TRY( mm_free(branch[ii1].br) );
	}
	RC_TRY( mm_free(branch) );
	RC_TRY( mm_free(queue) );
	node->renum_flag = RENUMBERED;
	return(NORMAL_RC);
}


static RC
remove_queue (int rm_index, int *queue_size, int *queue,
              NODAL_BRANCH *branch, int *next_index, int type_flag)
{
	int ii1;
	int max_index;
	double max_score;

	/* queue から rm_index を削除 */
	for(ii1=0; ii1<(*queue_size); ii1++){
		if(queue[ii1] == rm_index) break;
	}
	if(ii1 >= (*queue_size)) return(SEEK_ERROR_RC);
	queue[ii1] = queue[(*queue_size) - 1];
	(*queue_size)--;

	/* rm_index の結合を削除                          */
	/* rm_index に結合する branch 要素を queue に追加 */
	for(ii1=0; ii1<branch[rm_index].size; ii1++){
		RC_TRY( remove_branch(&(branch[branch[rm_index].br[ii1]]), rm_index) );
		RC_TRY( add_queue(branch[rm_index].br[ii1], queue_size, queue) );
	}

	/* weight を増加させ、weight が高く、size が小さい要素を next_index に */
	max_index = -1;
	max_score = -(DBL_MAX/2.0);
	for(ii1=0; ii1<(*queue_size); ii1++){
		double  score = 0.0;

		branch[queue[ii1]].weight++;
		if(type_flag == 0){
			score = branch[queue[ii1]].weight - branch[queue[ii1]].size;
		}else{
			double tmp1 = (double)branch[queue[ii1]].weight;
			double tmp2 = (double)count_queue(queue[ii1], branch,
			                                  *queue_size, queue);
			score = tmp1*tmp1 - tmp2*tmp2;
		}

		if((score > max_score)||(max_index < 0)){
			max_score = score;
			max_index = queue[ii1];
		}
	}
	(*next_index) = max_index;

	return(NORMAL_RC);
}


/* branch[index].size の中でまだ queue に登録されていない節点をカウント */
static int
count_queue (int index, const NODAL_BRANCH *branch,
             int queue_size, const int *queue)
{
	int ii1, ii2;
	int ret = 0;

	for(ii1=0; ii1<(branch[index].size); ii1++){
		for(ii2=0; ii2<queue_size; ii2++){
			if(queue[ii2] == branch[index].br[ii1]){
				break;
			}
		}
		if(ii2 >= queue_size) ret++;
	}

	return(ret);
}


/* queue に new_index を登録 */
static RC
add_queue (int new_index, int *queue_size, int *queue)
{
	int ii1;

	for(ii1=0; ii1<(*queue_size); ii1++){
		if(queue[ii1] == new_index){   /* 登録済み */
			return(NORMAL_RC);
		}
	}
	queue[*queue_size] = new_index;
	(*queue_size)++;

	return(NORMAL_RC);
}


static RC
init_branch (NODE_ARRAY node, ELEMENT_ARRAY element,
             RIGID_ELEMENT_ARRAY rigid, NODAL_BRANCH *branch)
{
	int ii1, ii2, ii3;
	int index[MAX_NODE];
	int index2, index3;

	for(ii1=0; ii1<node.size; ii1++){
		branch[ii1].size = 0;
		branch[ii1].weight = 0;
		branch[ii1].alloc_size = MAX_NODE;
		branch[ii1].br = mm_alloc(sizeof(int)*branch[ii1].alloc_size);
		if(branch[ii1].br == NULL) return(ALLOC_ERROR_RC);
	}

	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;
		if(element.array[ii1].node_num > MAX_NODE) return(ARG_ERROR_RC);

		for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
			index[ii2] = search_node_label(node, element.array[ii1].node[ii2]);
			RC_NEG_CHK( index[ii2] );
		}

		for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
			for(ii3=0; ii3<element.array[ii1].node_num; ii3++){
				if(ii3 == ii2) continue;
				RC_TRY( add_branch(&(branch[index[ii2]]), index[ii3]) );
			}
		}
	}

	for(ii1=0; ii1<rigid.size; ii1++){
		if(rigid.array[ii1].label < 0) continue;
		if(rigid.array[ii1].node_num <= 0) return(ARG_ERROR_RC);

		for(ii2=0; ii2<rigid.array[ii1].node_num; ii2++){
			index2 = search_node_label(node, rigid.array[ii1].node[ii2]);
			RC_NEG_CHK( index2 );

			for(ii3=0; ii3<rigid.array[ii1].node_num; ii3++){
				if(ii3 == ii2) continue;
				index3 = search_node_label(node, rigid.array[ii1].node[ii3]);
				RC_NEG_CHK( index3 );
				RC_TRY( add_branch(&(branch[index2]), index3) );
			}
		}
	}

	return(NORMAL_RC);
}


static RC
remove_branch (NODAL_BRANCH *branch_ptr, int rm_br)
{
	int ii1;

	for(ii1=0; ii1<(branch_ptr->size); ii1++){
		if(branch_ptr->br[ii1] == rm_br){
			branch_ptr->br[ii1] = branch_ptr->br[branch_ptr->size - 1];
			(branch_ptr->size)--;
			return(NORMAL_RC);
		}
	}

	return(SEEK_ERROR_RC);
}


static RC
add_branch (NODAL_BRANCH *branch_ptr, int new_br)
{
	int ii1;

	for(ii1=0; ii1<(branch_ptr->size); ii1++){
		if(branch_ptr->br[ii1] == new_br){
			return(NORMAL_RC);
		}
	}

	if((branch_ptr->size) >= branch_ptr->alloc_size){
		branch_ptr->alloc_size += MAX_NODE;
		branch_ptr->br = mm_realloc(branch_ptr->br,
		                         sizeof(int)*branch_ptr->alloc_size);
		if(branch_ptr->br == NULL) return(ALLOC_ERROR_RC);
	}
	branch_ptr->br[branch_ptr->size] = new_br;
	(branch_ptr->size)++;

	return(NORMAL_RC);
}


/* node_label に対応する添字(renum を考慮) */
int
node_renum_index (NODE_ARRAY node, int node_label)
{
	int index;

	index = search_node_label(node, node_label);
	if(index < 0) return(-1);

	index = node.array[index].renum;
	if((index < 0)||(index >= node.size)) return(-2);

	return(index);
}


/* node_renum_index() のキャッシュ使用版                                  */
/* キャッシュ領域 cache[cache_size][2] は、事前に負値で初期化しておくこと */
int
node_renum_index_cache (NODE_ARRAY node, int node_label,
                        int cache_size, int cache[][2])
{
	int divisor = cache_size/4;
	int cache_index;
	int ret;

	if(node_label < 0) return(node_label);
	if(cache_size < 0) return(cache_size);

	/* 4 way set associative cache */
	cache_index = 4*(node_label%divisor);
	if(cache[cache_index+0][0] == node_label) return(cache[cache_index+0][1]);
	if(cache[cache_index+1][0] == node_label) return(cache[cache_index+1][1]);
	if(cache[cache_index+2][0] == node_label) return(cache[cache_index+2][1]);
	if(cache[cache_index+3][0] == node_label) return(cache[cache_index+3][1]);
	
	ret = node_renum_index(node, node_label);

	cache[cache_index+3][0] = cache[cache_index+2][0];
	cache[cache_index+3][1] = cache[cache_index+2][1];
	cache[cache_index+2][0] = cache[cache_index+1][0];
	cache[cache_index+2][1] = cache[cache_index+1][1];
	cache[cache_index+1][0] = cache[cache_index+0][0];
	cache[cache_index+1][1] = cache[cache_index+0][1];
	cache[cache_index+0][0] = node_label;
	cache[cache_index+0][1] = ret;

	return(ret);
}


#if 0

/* パラメータ */
#define DIST_WEIGHT (0.9)
#define MULTIPLIER  (0.9)

/* center からの距離に応じて節点番号を付けかえる */
/* 使用不可（性能悪すぎ！！） */
RC
sphere_renumber (NODE_ARRAY *node, TRANSLATION3D center)
{
	int ii1, ii2;
	TRANSLATION3D near = center;
	double score, min_score;
	int min_index;

	/* retnum をリセット */
	for(ii1=0; ii1<node->size; ii1++){
		node->array[ii1].renum = -1;
	}

	for(ii1=0; ii1<node->size; ii1++){
		min_score = DBL_MAX;
		min_index = -1;

		for(ii2=0; ii2<node->size; ii2++){
			if(node->array[ii2].label < 0) continue;
			if(node->array[ii2].renum >= 0) continue;

			score = DIST_WEIGHT
			        * dist_point3d(node->array[ii2].p, center)
	    	      + (1.0 - DIST_WEIGHT)
			        * dist_point3d(node->array[ii2].p, near);
			if(score < min_score){
				min_score = score;
				min_index = ii2;
			}
		}
		if(min_index < 0) break;
		node->array[min_index].renum = ii1;

		near = add_translation3d(
		       mul_scalar_translation3d(MULTIPLIER, node->array[min_index].p),
		       mul_scalar_translation3d(1.0-MULTIPLIER, near) );
	}
	node->renum_flag = RENUMBERED;

	return(NORMAL_RC);
}


#endif /* 0 */



