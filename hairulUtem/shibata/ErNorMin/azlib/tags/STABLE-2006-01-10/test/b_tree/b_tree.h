/*********************************************************************
 * b_tree.h
 *
 * Copyright (C) 2004 AzLib Developers Group
 *
 * Written by
 *  <Taiki AOYAMA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: b_tree.h 455 2005-09-16 06:36:00Z sasaoka $ */

#ifndef B_TREE_H
#define B_TREE_H

#include <stdio.h>
#include <stdlib.h>

#define MAX_CHILD (5)                    /*B木の階数*/
#define HALF_CHILD ((MAX_CHILD + 1)/2)

typedef enum {
	INTERNAL,     /*内部ノード*/
	LEAF          /*リーフノード*/
} NODE_KIND;      /*ノードの種類*/


typedef enum {
	NOT_CHANGE,     /*データに変化無し*/
	REMOVED,        /*データ消去*/
	REFORM          /*隣りのノードとの子の結合or再分配*/
} DELETE_RESULT;    /*データの消去結果*/


typedef struct node {
	NODE_KIND kind;                             /*ノードの種類*/
	union {
		struct {
			int child_num;                      /*子の数*/
			int child_key[MAX_CHILD];           /*Key値の配列*/
			struct node *child[MAX_CHILD];      /*子のノードポインタの配列*/
		} internal;                             /*内部ノード*/
		struct {
			int leaf_key;                       /*リーフ(葉)のKey値*/
			int data;                           /*リーフの持つデータ(int型)*/
		} leaf;                                 /*リーフノード*/
	} content;
} NODE;


/* b_tree.c */
/*Key値の挿入(返り値はリーフノードへのポインタ*/
NODE *insert(NODE **root, int key);

/*該当するKey値のデータを消去．成功なら1，該当データ無しなら0を返す*/
int delete(NODE **root, int key);

/*Key値に対するリーフノードをサーチ(返り値はリーフノードへのポインタ)*/
NODE *search(NODE *root, int key);

/*NODE型の確保*/
void allocate_node(NODE **node);


#endif /* B_TREE_H */
