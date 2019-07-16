/*********************************************************************
 * b_tree.c
 *
 * Copyright (C) 2004 AzLib Developers Group
 *
 * Written by
 *  <Taiki AOYAMA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: b_tree.c 842 2006-12-04 06:21:20Z sasaoka $ */

/* B木の試作プログラム
 *
 * データ構造理解のため，B木のTree部分を作成．
 * 検索(search)，追加(insert)，削除(delete)関数を作成．
 *
 * ただし，リーフノードに格納するデータは簡単にint型の整数としている．
 * B木の階数はデフォルトで5である．
 */

#include <stdio.h>
#include <stdlib.h>
#include "b_tree.h"
#include "memory_manager.h"
#include "rc.h"


static NODE *insert_seek(NODE **node, int key, NODE **grown_node,
                                                 int *grown_key);
static int delete_seek(NODE *node, int key, int *result);
static int joining_node(NODE *node, int index);
static int locate_subtree(NODE *node, int key);


void
allocate_node (NODE **node)
{
	(*node) = (NODE *)mm_alloc(1*sizeof(NODE));
	if(*node == NULL){
		error_printf(9001, sizeof(NODE));
		exit(EXIT_FAILURE);
	}
}


NODE *
insert (NODE **root, int key)
{
	/* 未挿入 */
	if(*root == NULL){
		allocate_node(root);
		(*root)->kind = LEAF;
		(*root)->content.leaf.leaf_key = key;
		return(*root);
	}else{
		NODE *ret, *new_node, *grown_node;
		int grown_key;

		ret = insert_seek(root, key, &grown_node, &grown_key);

		/* 木が分割され新しく作られた部分木がある */
		/* 高さが１段上がる                       */
		if(grown_node != NULL){
			allocate_node(&new_node);
			new_node->kind = INTERNAL;
			new_node->content.internal.child_num = 2;
			new_node->content.internal.child[0] = *root;
			new_node->content.internal.child[1] = grown_node;
			new_node->content.internal.child_key[1] = grown_key;
			*root = new_node;
		}
		return(ret);
	}
}


/* 木をたどりデータを挿入すべきリーフノードへのポインタを得る   */
/* grown_nodeは挿入の仮定で新しく作られた木，リーフへのポインタ */
static NODE *
insert_seek (NODE **node, int key, NODE **grown_node, int *grown_key)
{
	NODE *node_temp;
	*grown_node = NULL;

	node_temp = *node;

	/* nodeがリーフノード */
	if(node_temp->kind == LEAF){

		/* 既存のリーフノードのKey値と一致 */
		if(node_temp->content.leaf.leaf_key == key){
			return(*node);
		/* Key値は一致せず新しくリーフノード作成 */
		}else{
			NODE *new_leaf;
			 
			allocate_node(&new_leaf);
			new_leaf->kind = LEAF;
			new_leaf->content.leaf.leaf_key = key;

			/* 以下は常に右隣り(grown_node)のKey値が常に */
			/* 大きくなるようにするための処理である      */

			/* nodeとnew_leafを入れ替え            */
			/* grown_nodeにnode_tempが来るよう変更 */
			/* node_tempのKey値 > new_leafのKey値  */
			if(key < node_temp->content.leaf.leaf_key){
				*node = new_leaf;
				*grown_key = node_temp->content.leaf.leaf_key;
				*grown_node = node_temp;
			/* grown_nodeにnew_leafを代入         */
			/* 入れ替えなし                       */
			/* node_tempのKey値 < new_leafのKey値 */
			}else{
				*grown_key = key;
				*grown_node = new_leaf;
			}

			return(new_leaf);
		}
	/* nodeが内部ノード */
	}else{
		int ii1, ii2;
		int pos;
		int new_key;
		NODE *new_node;
		NODE *ret;

		pos = locate_subtree(node_temp, key);

		/* 更に木をたどる(再帰呼び出し) */
		ret = insert_seek(&(node_temp->content.internal.child[pos]), key,
		                                             &new_node, &new_key);
		if(new_node == NULL) return(ret);

		/* 木が分割され新しく作られた部分木がある        */
		/* node_tempの子の数が新しく作られた子(new_node) */
		/* を加えてもMAX_CHILDを越えない                 */
		if(node_temp->content.internal.child_num < MAX_CHILD){
			for(ii1=node_temp->content.internal.child_num-1; ii1>pos; ii1--){
				node_temp->content.internal.child[ii1+1] =
				node_temp->content.internal.child[ii1];
				node_temp->content.internal.child_key[ii1+1] =
				node_temp->content.internal.child_key[ii1];
			}
			node_temp->content.internal.child[pos+1] = new_node;
			node_temp->content.internal.child_key[pos+1] = new_key;
			node_temp->content.internal.child_num++;
			return(ret);
		/* node_tempの子の数に新しく作られた子(new_node)を加えると        */
		/* オーバーする．node_temp の子をB木の条件を満足するよう2つに分割 */
		}else{
			NODE *divide_node;

			allocate_node(&divide_node);
			divide_node->kind = INTERNAL;

			/* new_nodeをnode_tempへ入れる */
			if( pos < (HALF_CHILD-1) ){
				ii2 = 0;
				/* node_tempの子(HALF_CHILDからMAX_CHILDをdivideへ移す */
				for(ii1=(HALF_CHILD-1); ii1<MAX_CHILD; ii1++){
					divide_node->content.internal.child[ii2] = 
					node_temp->content.internal.child[ii1];
					divide_node->content.internal.child_key[ii2] = 
					node_temp->content.internal.child_key[ii1];
					ii2++;
				}
				/* pos+1までを１つ右へ詰める */
				for(ii1=(HALF_CHILD-2); ii1>pos; ii1--){
					node_temp->content.internal.child[ii1+1] = 
					node_temp->content.internal.child[ii1];
					node_temp->content.internal.child_key[ii1+1] = 
					node_temp->content.internal.child_key[ii1];
				}
				/* new_nodeを代入 */
				node_temp->content.internal.child[pos+1] = new_node;
				node_temp->content.internal.child_key[pos+1] = new_key;
			/* new_nodeをdivide_nodeへ入れる */
			}else{
				/* node_tempの子(MAX_CHILDからHALF_CHILD+1)を  */
				/* divide_nodeへ代入適する位置にnew_nodeを入力 */
				ii2 = MAX_CHILD - HALF_CHILD;
				for(ii1=(MAX_CHILD-1); ii1>=HALF_CHILD; ii1--){
					if(ii1 == pos){
						divide_node->content.internal.child[ii2] = new_node;
						divide_node->content.internal.child_key[ii2] = new_key;
						ii2--;
					}
					divide_node->content.internal.child[ii2] = 
					node_temp->content.internal.child[ii1];
					divide_node->content.internal.child_key[ii2] = 
					node_temp->content.internal.child_key[ii1];
					ii2--;
				}
				if(pos < HALF_CHILD){
					divide_node->content.internal.child[0] = new_node;
					divide_node->content.internal.child_key[0] = new_key;
				}
			}

			/* node_tempとdivide_nodeの子の数を入力 */
			node_temp->content.internal.child_num = HALF_CHILD;
			divide_node->content.internal.child_num = (MAX_CHILD + 1)
			                                        - HALF_CHILD;

			/* 新しく作られた部分木としてdivide_nodeを返す */
			*grown_node = divide_node;
			*grown_key = divide_node->content.internal.child_key[0];
			return(ret);
		}
	}
}


NODE *
search (NODE *root, int key)
{
	int ii1;
	NODE *node;

	if(root == NULL) return(NULL);

	/* リーフノードまで内部ノードをたどる */
	node = root;
	while(node->kind == INTERNAL){
		ii1 = locate_subtree(node, key);
		node = node->content.internal.child[ii1];
	}

	/* 該当するリーフノード発見 */
	if(key == node->content.leaf.leaf_key){
		return(node);
	/* 発見できず */
	}else{
		return(NULL);
	}
}


int
delete (NODE **root, int key)
{
	int d_ret, result;
	NODE *temp;

	/* 未挿入 */
	if(root == NULL) return(0);

	d_ret = delete_seek(*root, key, &result);
	/* rootの子がすべて消去された */
	if(result == REMOVED){
		*root = NULL;
	/* rootの子がHALF_CHILDより少なく，1つのみであるから，高さを1段低くする */
	}else if(result == REFORM && (*root)->content.internal.child_num == 1){
		temp = *root;
		*root = (*root)->content.internal.child[0];
		if(mm_free(temp) != NORMAL_RC){
			error_printf(99999);
			exit(EXIT_FAILURE);
		}
	}

	return(d_ret);
}


/*
 * 木をたどり該当するKey値のデータを消去
 * 成功なら1，該当データ無しなら0を返す
 * resultには
 *  変化無し      : NOT_CHANGE
 *  結合可能性あり:REFORM
 *  消去          :REMOVED
 * が帰る
 */ 
static int
delete_seek (NODE *node, int key, int *result)
{
	/* 消去結果resultの初期化 */
	*result = NOT_CHANGE;

	/* リーフノードである */
	if(node->kind == LEAF){
		/* リーフノードが消去される */
		if(node->content.leaf.leaf_key == key){
			*result = REMOVED;
			if(mm_free(node) != NORMAL_RC){
				error_printf(99999);
				exit(EXIT_FAILURE);
			}
			return(1);
		/* 消去すべきKey値なし */
		}else{
			return(0);
		}
	/* 内部ノードである */
	}else{
		int ii1;
		int pos;
		int d_condition;
		int d_ret;
		int sub_pos;
		int docking_flag = 0;

		pos = locate_subtree(node, key);
		/* 更に木をたどる(再帰呼び出し) */
		d_ret = delete_seek(node->content.internal.child[pos],key,&d_condition);

		/* データに変化なし */
		if(d_condition == NOT_CHANGE){
			return(d_ret);
		}

		/* node->child[pos]とnode->child[pos-1]の結合または子の再分配 */
		/* pos=0の場合，node->child[pos] と node->child[pos+1]を      */
		/* 結合または再分配                                           */
		if(d_condition == REFORM){
			sub_pos = (pos == 0) ? 0 : pos - 1;
			docking_flag = joining_node(node, sub_pos);
			if(docking_flag){
				pos = sub_pos + 1;
			}
		}

		/* データが消去された */
		/* 子の結合があった   */
		if(d_condition == REMOVED || docking_flag){
			int n = node->content.internal.child_num;
			/* 空きがあるので１つ左に詰める */
			for(ii1=pos; ii1<(n-1); ii1++){
				node->content.internal.child[ii1] =
				node->content.internal.child[ii1+1];
				node->content.internal.child_key[ii1] =
				node->content.internal.child_key[ii1+1];
			}

			/* nodeの子と他のnode(左隣or右隣)の子との結合可能性の判定 */
			if(--node->content.internal.child_num < HALF_CHILD){
				*result = REFORM;
			}
		}

		return(d_ret);
	}
}


/* B木の条件より                                        */
/* node->child[index]とnode->child[index+1]の結合を行う */
/* 結合不能なら子の再分配を行う                         */
static int
joining_node (NODE *node, int index)
{
	int ii1;
	int n1;
	int n2;
	NODE *node1;
	NODE *node2;

	node1 = node->content.internal.child[index];
	node2 = node->content.internal.child[index+1];

	node2->content.internal.child_key[0] =
	node->content.internal.child_key[index+1];

	n1 = node1->content.internal.child_num;
	n2 = node2->content.internal.child_num;

	/* 結合可能 */
	if( (n1 + n2) <= MAX_CHILD ){
		for(ii1=0; ii1<n2; ii1++){
			node1->content.internal.child[n1+ii1] =
			node2->content.internal.child[ii1];
			node1->content.internal.child_key[n1+ii1] =
			node2->content.internal.child_key[ii1];
		}
		node1->content.internal.child_num += n2;
		if(mm_free(node2) != NORMAL_RC){
			error_printf(99999);
			exit(EXIT_FAILURE);
		}

		return(1);

	/* 結合不能 */
	/* node1とnode2の間で子を再分配 */
	}else{
		int n;           /* node1に分配すべき子の数 */
		int move_num;    /* 移動する子の数 */

		n = (n1 + n2) / 2;

		/* 部分木 node1から node2へ子を移動 */
		if(n1 > n){
			move_num = n1 - n;
			/*node1の子が入るので node2の子を右へずらす*/
			for(ii1=n2-1; ii1>=0; ii1--){
				node2->content.internal.child[move_num+ii1] =
				node2->content.internal.child[ii1];
				node2->content.internal.child_key[move_num+ii1] =
				node2->content.internal.child_key[ii1];
			}
			/* 部分木 node1 から node2 へ move_num 個だけ子を移動 */
			for(ii1=0; ii1<move_num; ii1++){
				node2->content.internal.child[ii1] =
				node1->content.internal.child[n+ii1];
				node2->content.internal.child_key[ii1] =
				node1->content.internal.child_key[n+ii1];
			}

		/* 部分木 node2 から node1 へ子を移動 */
		}else{
			move_num = n - n1;
			/* 部分木 node2 から node1へ move_num 個だけ子を移動 */
			for(ii1=0; ii1<move_num; ii1++){
				node1->content.internal.child[n1+ii1] =
				node2->content.internal.child[ii1];
				node1->content.internal.child_key[n1+ii1] =
				node2->content.internal.child_key[ii1];
			}
			/* 空きができたので node2の子を左へ詰める */
			for(ii1=0; ii1<(n2-move_num); ii1++){
				node2->content.internal.child[ii1] =
				node2->content.internal.child[move_num+ii1];
				node2->content.internal.child_key[ii1] =
				node2->content.internal.child_key[move_num+ii1];
			}
		}

		/* 子の数を更新する */
		node1->content.internal.child_num = n;
		node2->content.internal.child_num = n1 + n2 - n;
		node->content.internal.child_key[index+1] =
		node2->content.internal.child_key[0];

		return(0);
	}
}


/* 挿入または変更する位置探し(返り値はchild_keyのindex),  線形探索*/
static int
locate_subtree (NODE *node, int key)
{
	int ii1;

	for(ii1=node->content.internal.child_num-1; ii1>0; ii1--){
		if(key >= node->content.internal.child_key[ii1]) return(ii1);
	}

	return(0);
}
