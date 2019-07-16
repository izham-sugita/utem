/*********************************************************************
 * extract_surface.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Tomoyuki TSUBATA>
 *
 * Thanks to following contributors.
 *  <Ryu SASAOKA> 
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: extract_surface.c 940 2010-01-20 06:15:24Z matsushita $ */

#include <stdio.h>
#include <stdlib.h>
#include "rc.h"
#include "memory_manager.h"
#include "math_utl.h"
#include "fem_solver.h"
#include "fem_struct.h"
#include "avs_utl.h"

#define I_RANDOM  16384
#define MAX_NODE_LIGHT  8
#define EDGE_MAX_FROM  16
#define CUT_MAX_NODE   81
static int i_random[I_RANDOM];

typedef struct {
	ELEM_TYPE type;
	int node_num;
	int node[MAX_NODE_LIGHT];
	int from_num;
	int from[2];   /* 元の要素label */
	int from_f[2]; /* cut_element 専用パラメータ */
	int tmp;       /* cut_element 専用パラメータ */
} ELEMENT_FACE;

typedef struct {
	ELEM_TYPE type;  /* ELEM_LINE1 or ELEM_LINE2 */
	int node_num;    /*     2      or     3      */
	int node[3];
	int from_num;
	int from[EDGE_MAX_FROM];   /* 元の要素label */
	int from_e[EDGE_MAX_FROM]; /* cut_element 専用パラメータ */
} ELEMENT_EDGE;

typedef struct {
	int label;
	ELEM_TYPE type;
	int material;
	int physical;
	int node_num;
	int node[CUT_MAX_NODE];
} ELEMENT_NODE_BOX;


/* local function */
static void init_i_random(int hash_size);
static int elem_hash_index(int hash_size, ELEMENT_FACE *hash,
                           int key_size, const int *keys);
static int equal_keys(int key_size, const int *keys1, const int *keys2);
static void sort_keys(int key_size, const int *source, int *dest);
static int first_index(int hash_size, int key_size, const int *keys);
static int cal_hash_size(ELEMENT_ARRAY elem);
static RC extract_face(ELEMENT_ARRAY elem, ELEMENT_ARRAY *face, int n);
static RC extract_face_light(ELEMENT_ARRAY elem, ELEMENT_FACE **hash,
                             int *hash_size);
static RC input_hash(ELEMENT_ARRAY elem, ELEMENT_FACE *hash, int hash_size);
static RC hash2elem(int hash_size, const ELEMENT_FACE *hash,
                    int n, ELEMENT_ARRAY *face);
static RC clean_hash(ELEMENT_FACE *hash, int *hash_size);
static RC set_hash_face(int hash_size, ELEMENT_FACE *hash,
                        int key_size, const int *keys, int label,
                        ELEM_TYPE type, int number);
static RC set_hash_tri1(int hash_size, ELEMENT_FACE *hash, ELEMENT elem);
static RC set_hash_tri2(int hash_size, ELEMENT_FACE *hash, ELEMENT elem);
static RC set_hash_quad1(int hash_size, ELEMENT_FACE *hash, ELEMENT elem);
static RC set_hash_quad2(int hash_size, ELEMENT_FACE *hash, ELEMENT elem);
static RC set_hash_tetra1(int hash_size, ELEMENT_FACE *hash, ELEMENT elem);
static RC set_hash_tetra2(int hash_size, ELEMENT_FACE *hash, ELEMENT elem);
static RC set_hash_penta1(int hash_size, ELEMENT_FACE *hash, ELEMENT elem);
static RC set_hash_penta2(int hash_size, ELEMENT_FACE *hash, ELEMENT elem);
static RC set_hash_hexa1(int hash_size, ELEMENT_FACE *hash, ELEMENT elem);
static RC set_hash_hexa2(int hash_size, ELEMENT_FACE *hash, ELEMENT elem);
static RC set_cut_elem(ELEMENT_ARRAY *elem, ELEMENT_NODE_BOX box,
                       int index_size, const int *index);
static RC input_cut_elem_line1(ELEMENT_ARRAY *elem, ELEMENT_NODE_BOX box);
static RC input_cut_elem_line2(ELEMENT_ARRAY *elem, ELEMENT_NODE_BOX box);
static RC input_cut_elem_tri1(ELEMENT_ARRAY *elem, ELEMENT_NODE_BOX box);
static RC input_cut_elem_tri2(ELEMENT_ARRAY *elem, ELEMENT_NODE_BOX box);
static RC input_cut_elem_quad1(ELEMENT_ARRAY *elem, ELEMENT_NODE_BOX box);
static RC input_cut_elem_quad2(ELEMENT_ARRAY *elem, ELEMENT_NODE_BOX box);
static RC input_cut_elem_tetra1(ELEMENT_ARRAY *elem, ELEMENT_NODE_BOX box);
static RC input_cut_elem_tetra2(ELEMENT_ARRAY *elem, ELEMENT_NODE_BOX box);
static RC input_cut_elem_penta1(ELEMENT_ARRAY *elem, ELEMENT_NODE_BOX box);
static RC input_cut_elem_penta2(ELEMENT_ARRAY *elem, ELEMENT_NODE_BOX box);
static RC input_cut_elem_hexa1(ELEMENT_ARRAY *elem, ELEMENT_NODE_BOX box);
static RC input_cut_elem_hexa2(ELEMENT_ARRAY *elem, ELEMENT_NODE_BOX box);
static RC extract_edge(ELEMENT_ARRAY elem, int *edge_size, ELEMENT_EDGE **edge);
static int count_node_box_size(ELEM_TYPE type);
static int cal_edge_size(ELEMENT_ARRAY elem);
static RC interpolate_coord_edge(NODE_ARRAY *node, ELEMENT_ARRAY elem,
                                 ELEMENT_EDGE edge,
                                 ELEMENT_NODE_BOX *box, int *new_label);
static RC interpolate_coord_face(NODE_ARRAY *node, ELEMENT_ARRAY elem,
                                 ELEMENT_FACE face,
                                 ELEMENT_NODE_BOX *box, int *new_label);
static RC interpolate_coord_elem(NODE_ARRAY *node, ELEMENT_ARRAY elem,
                                 ELEMENT element,
                                 ELEMENT_NODE_BOX *box, int *new_label);
static RC local_inside_coord(ELEM_TYPE type, double **coord, int *number);
static RC cal_interpolate_node(NODE_ARRAY node, const double *N, 
                               int node_num, const int *node_label,
                               VECT3D *new_coord);
static RC add_node_array(NODE_ARRAY *node, VECT3D coord, int new_label);
static int check_rotation_edge(ELEM_TYPE type, int *index);
static RC set_edge(int hash_size, ELEMENT_EDGE *hash, int key_size,
                   const int *keys, int label, ELEM_TYPE type, int number);
static int elem_edge_index(int hash_size, ELEMENT_EDGE *hash,
                           int key_size, const int *keys);
static RC set_edge_tri1(int hash_size, ELEMENT_EDGE *hash, ELEMENT elem);
static RC set_edge_tri2(int hash_size, ELEMENT_EDGE *hash, ELEMENT elem);
static RC set_edge_quad1(int hash_size, ELEMENT_EDGE *hash, ELEMENT elem);
static RC set_edge_quad2(int hash_size, ELEMENT_EDGE *hash, ELEMENT elem);
static RC set_edge_tetra1(int hash_size, ELEMENT_EDGE *hash, ELEMENT elem);
static RC set_edge_tetra2(int hash_size, ELEMENT_EDGE *hash, ELEMENT elem);
static RC set_edge_penta1(int hash_size, ELEMENT_EDGE *hash, ELEMENT elem);
static RC set_edge_penta2(int hash_size, ELEMENT_EDGE *hash, ELEMENT elem);
static RC set_edge_hexa1(int hash_size, ELEMENT_EDGE *hash, ELEMENT elem);
static RC set_edge_hexa2(int hash_size, ELEMENT_EDGE *hash, ELEMENT elem);


RC
extract_surface (ELEMENT_ARRAY elem, ELEMENT_ARRAY *surf)
{
	RC_TRY( extract_face(elem, surf, 1) );
	return(NORMAL_RC);
}


RC
extract_divider (ELEMENT_ARRAY elem, ELEMENT_ARRAY *div)
{
	RC_TRY( extract_face(elem, div, 2) );
	return(NORMAL_RC);
}


/* 実行結果は、実行前の節点、要素データを上書きする */
RC
cut_element (NODE_ARRAY *node, ELEMENT_ARRAY *elem)
{
	ELEMENT_EDGE *edge = NULL;
	ELEMENT_FACE *face = NULL;
	ELEMENT_NODE_BOX *box = NULL;
	int edge_size;
	int face_size;
	int box_size;
	int number;
	int ii1, ii2;

	RC_TRY( clean_element_array(elem) );
	
	/* 分割後の要素と節点の関係を記憶するための、構造体の配列を確保する */
	box_size = elem->size;
	box = mm_alloc( box_size * sizeof(ELEMENT_NODE_BOX) );
	if(box == NULL) return(ALLOC_ERROR_RC);

	for(ii1=0; ii1<box_size; ii1++){
		box[ii1].label = elem->array[ii1].label;
		box[ii1].type = elem->array[ii1].type;
		box[ii1].material = elem->array[ii1].material;
		box[ii1].physical = elem->array[ii1].physical;
		RC_NEG_CHK( box[ii1].node_num
		          = count_node_box_size(elem->array[ii1].type) );
		for(ii2=0; ii2<elem->array[ii1].node_num; ii2++){
			box[ii1].node[ii2] = elem->array[ii1].node[ii2];
		}
		for(ii2=elem->array[ii1].node_num; ii2<CUT_MAX_NODE; ii2++){
			box[ii1].node[ii2] = -1;
		}
	}

	/* 節点の最大label を求める */
	number = 0;
	for(ii1=0; ii1<node->size; ii1++){
		if(number < node->array[ii1].label) number = node->array[ii1].label;
	}

	/* extract_edge(),  extract_face_light() の必要性を調べる */
	/* edge, face内の点を、新たな節点として追加する           */
	switch( analysis_dim(*elem) ){
	case 1:
		break;
	case 2:
		RC_TRY( extract_edge(*elem, &edge_size, &edge) );
		for(ii2=0; ii2<edge_size; ii2++){
			RC_TRY( interpolate_coord_edge(node, *elem, edge[ii2], box,
			                               &number) );
		}
		RC_TRY( mm_free((void *)edge) );
		break;
	case 3:
		RC_TRY( extract_edge(*elem, &edge_size, &edge) );
		for(ii2=0; ii2<edge_size; ii2++){
			RC_TRY( interpolate_coord_edge(node, *elem, edge[ii2], box,
			                               &number) );
		}
		RC_TRY( mm_free((void *)edge) );

		RC_TRY( extract_face_light(*elem, &face, &face_size) );
		for(ii1=0; ii1<face_size; ii1++){
			RC_TRY( interpolate_coord_face(node, *elem, face[ii1], box,
			                               &number) );
		}
		RC_TRY( mm_free((void *)face) );
		break;
	default:
		return(ARG_ERROR_RC);
	}

	/* elem 内の点を、新たな節点として追加する */
	for(ii1=0; ii1<elem->size; ii1++){
		RC_TRY( interpolate_coord_elem(node, *elem, elem->array[ii1], box,
		                               &number) );
	}
	RC_TRY( clean_node_array(node) );
	RC_TRY( free_element_array(elem) );
	RC_TRY( allocate_element_array(0, elem) );

	/* 新たな要素情報を記憶する */
	for(ii1=0; ii1<box_size; ii1++){
		switch(box[ii1].type){
		case ELEM_LINE1:
			RC_TRY( input_cut_elem_line1(elem, box[ii1]) );
			break;
		case ELEM_LINE2:
			RC_TRY( input_cut_elem_line2(elem, box[ii1]) );
			break;
		case ELEM_TRI1:
			RC_TRY( input_cut_elem_tri1(elem, box[ii1]) );
			break;
		case ELEM_TRI2:
			RC_TRY( input_cut_elem_tri2(elem, box[ii1]) );
			break;
		case ELEM_QUAD1:
			RC_TRY( input_cut_elem_quad1(elem, box[ii1]) );
			break;
		case ELEM_QUAD2:
			RC_TRY( input_cut_elem_quad2(elem, box[ii1]) );
			break;
		case ELEM_TETRA1:
			RC_TRY( input_cut_elem_tetra1(elem, box[ii1]) );
			break;
		case ELEM_TETRA2:
			RC_TRY( input_cut_elem_tetra2(elem, box[ii1]) );
			break;
		case ELEM_PENTA1:
			RC_TRY( input_cut_elem_penta1(elem, box[ii1]) );
			break;
		case ELEM_PENTA2:
			RC_TRY( input_cut_elem_penta2(elem, box[ii1]) );
			break;
		case ELEM_HEXA1:
			RC_TRY( input_cut_elem_hexa1(elem, box[ii1]) );
			break;
		case ELEM_HEXA2:
			RC_TRY( input_cut_elem_hexa2(elem, box[ii1]) );
			break;
		default:
			return(IMPLEMENT_ERROR_RC);
		}
	}
	RC_TRY( mm_free((void *)box) );

	return(NORMAL_RC);
}


static void
init_i_random (int hash_size)
{
	int ii1;

	init_d_random(I_RANDOM);
	for(ii1=0; ii1<I_RANDOM; ii1++){
		i_random[ii1] = (int)((double)hash_size * d_random());
		if(i_random[ii1] < 0) i_random[ii1] = 0;
		if(i_random[ii1] >= hash_size) i_random[ii1] = hash_size;
	}
}


/* return index of hash */
static int
elem_hash_index (int hash_size, ELEMENT_FACE *hash,
                 int key_size, const int *keys)
{
	int skeys1[MAX_NODE_LIGHT];
	int skeys2[MAX_NODE_LIGHT];
	int ii1, ii2;
	int index;
	int inc;

	sort_keys(key_size, keys, skeys1);
	index = first_index(hash_size, key_size, skeys1);
	inc = 1 << (skeys1[0]%6);   /* 1, 2, 4, 8, 16 or 32 */

	for(ii1=0; ii1<hash_size; ii1++, index+=inc){
		index %= hash_size;
		if(hash[index].from_num < 0) break;
		if(hash[index].node_num != key_size) continue;
		sort_keys(key_size, hash[index].node, skeys2);
		if(equal_keys(key_size, skeys1, skeys2)) break;
	}
	if(ii1 >= hash_size){
		return(-1);     /* hash over flow */
	}

	/* set data */
	if(hash[index].from_num < 0){
		hash[index].from_num = 0;
		hash[index].node_num = key_size;
		for(ii2=0;ii2<key_size;ii2++){
			hash[index].node[ii2] = keys[ii2];
		}
	}
	hash[index].tmp = keys[0];

	return(index);
}


/* if keys1[] == keys2[] return 1 */
/* else return 0                  */
/* key?[] must be sorted          */
static int
equal_keys (int key_size, const int *keys1, const int *keys2)
{
	int ii1;

	for(ii1=0; ii1<key_size; ii1++){
		if(keys1[ii1] != keys2[ii1]) return(0);
	}
	return(1);
}


/* copy from source to dest and sort */
static void
sort_keys (int key_size, const int *source, int *dest)
{
	int ii1, ii2;
	int tmp;

	/* copy like a bubble sort algorithm */
	for(ii1=0; ii1<key_size; ii1++){
		dest[ii1] = source[ii1];
		for(ii2=ii1-1; ii2>=0; ii2--){
			if(dest[ii2] > dest[ii2+1]){
				tmp = dest[ii2];
				dest[ii2] = dest[ii2+1];
				dest[ii2+1] = tmp;
			}else{
				break;
			}
		}
	}
}


/* first hash index from keys */
static int
first_index (int hash_size, int key_size, const int *keys)
{
	int ii1;
	int key_sum = 0;
	int ret = 0;

	for(ii1=0;ii1<key_size;ii1++){
		key_sum = (key_sum + keys[ii1]) % I_RANDOM;
		ret += (i_random[key_sum] >> ii1);
	}
	ret %= hash_size;

	return(ret);
}


#define CONST_TERM   (128.0)
static int
cal_hash_size (ELEMENT_ARRAY elem)
{
	int ii1;
	double ret = CONST_TERM;

	for(ii1=0; ii1<elem.size; ii1++){
		switch(elem.array[ii1].type){
		case ELEM_TRI1:
		case ELEM_TRI2:
			ret += 3.0;
			break;
		case ELEM_QUAD1:
		case ELEM_QUAD2:
		case ELEM_QUAD2L:
			ret += 4.0;
			break;
		case ELEM_TETRA1:
		case ELEM_TETRA2:
			ret += 4.0;
			break;
		case ELEM_PENTA1:
		case ELEM_PENTA2:
			ret += 5.0;
			break;
		case ELEM_HEXA1:
		case ELEM_HEXA2:
			ret += 6.0;
			break;
		default :
			ret += 6.0;
		}
	}

	/* 32 の倍数 + 1 */
	return( (32*((int)(ret/32.0) + 1) + 1) );
}


/* n=1 では表面を、n=2 では表面以外を、n=3 では両方を求める  */ 
static RC
extract_face (ELEMENT_ARRAY elem, ELEMENT_ARRAY *face, int n)
{
	ELEMENT_FACE *hash = NULL;
	int hash_size;

	/* calculate size of hash */
	if(elem.size <= 0) return(ARG_ERROR_RC);

	hash_size = cal_hash_size(elem);

	hash = mm_alloc( hash_size * sizeof(ELEMENT_FACE) );
	if(hash == NULL) return(ALLOC_ERROR_RC);

	RC_TRY( input_hash(elem, hash, hash_size) );
	RC_TRY( hash2elem(hash_size, hash, n, face) );

	/* free hash[] */
	RC_TRY( mm_free(hash) );

	return(NORMAL_RC);
}


static RC
extract_face_light (ELEMENT_ARRAY elem, ELEMENT_FACE **hash, int *hash_size)
{
	/* calculate size of hash */
	if(elem.size <= 0) return(ARG_ERROR_RC);

	*hash_size = cal_hash_size(elem);

	(*hash) = mm_alloc( (*hash_size) * sizeof(ELEMENT_FACE) );
	if(*hash == NULL) return(ALLOC_ERROR_RC);

	RC_TRY( input_hash(elem, (*hash), *hash_size) );
	RC_TRY( clean_hash( (*hash), hash_size ) );

	return(NORMAL_RC);

}


static RC
input_hash (ELEMENT_ARRAY elem, ELEMENT_FACE *hash, int hash_size)
{
	int ii1;

	/* initialize hash */
	for(ii1=0; ii1<hash_size; ii1++){
		hash[ii1].from_num = -1;
	}

	/* initialize i_random */
	init_i_random(hash_size);

	for(ii1=0; ii1<elem.size; ii1++){
		if(elem.array[ii1].label < 0) continue;
		switch(elem.array[ii1].type){
		case ELEM_TRI1:
			RC_TRY( set_hash_tri1(hash_size, hash, elem.array[ii1]) );
			break;
		case ELEM_TRI2:
			RC_TRY( set_hash_tri2(hash_size, hash, elem.array[ii1]) );
			break;
		case ELEM_QUAD1:
			RC_TRY( set_hash_quad1(hash_size, hash, elem.array[ii1]) );
			break;
		case ELEM_QUAD2:
		case ELEM_QUAD2L:
			RC_TRY( set_hash_quad2(hash_size, hash, elem.array[ii1]) );
			break;
		case ELEM_TETRA1:
			RC_TRY( set_hash_tetra1(hash_size, hash, elem.array[ii1]) );
			break;
		case ELEM_TETRA2:
			RC_TRY( set_hash_tetra2(hash_size, hash, elem.array[ii1]) );
			break;
		case ELEM_PENTA1:
			RC_TRY( set_hash_penta1(hash_size, hash, elem.array[ii1]) );
			break;
		case ELEM_PENTA2:
			RC_TRY( set_hash_penta2(hash_size, hash, elem.array[ii1]) );
			break;
		case ELEM_HEXA1:
			RC_TRY( set_hash_hexa1(hash_size, hash, elem.array[ii1]) );
			break;
		case ELEM_HEXA2:
			RC_TRY( set_hash_hexa2(hash_size, hash, elem.array[ii1]) );
			break;
		case ELEM_BEAM1:
		case ELEM_BEAM2:
			break;
		default :
			return(IMPLEMENT_ERROR_RC);
		}
	}

	return(NORMAL_RC);

}


static RC
hash2elem (int hash_size, const ELEMENT_FACE *hash, int n, ELEMENT_ARRAY *face)
{
	int counter;
	int ii1, ii2;

	/* count num of face */
	counter = 0;
	for(ii1=0; ii1<hash_size; ii1++){
		if(hash[ii1].from_num == n) counter ++;
	}

	RC_TRY( allocate_element_array(counter, face) );

	face->source = FEM_GENERATED_SURFACE;
	face->size = 0;

	for(ii1=0; ii1<hash_size; ii1++){
		if(hash[ii1].from_num == n){
			init_element( &(face->array[face->size]) );
			face->array[face->size].label = face->size + 1;
			face->array[face->size].type = hash[ii1].type;
			face->array[face->size].node_num = hash[ii1].node_num;
			for(ii2=0; ii2<hash[ii1].node_num; ii2++){
				face->array[face->size].node[ii2] = hash[ii1].node[ii2];
			}
			for(ii2=0; ii2<n; ii2++){
				face->array[face->size].i_info[ii2] = hash[ii1].from[ii2];
				face->array[face->size].i_info[ii2+2] = hash[ii1].from_f[ii2];
			}
			(face->size)++;
		}
	}

	return(NORMAL_RC);
}


static RC
clean_hash (ELEMENT_FACE *hash, int *hash_size)
{
	int number;
	int ii1, ii2;

	number = *hash_size;
	*hash_size = 0;
	for(ii1=0; ii1<number; ii1++){
		if( hash[ii1].from_num <= 0 ) continue;
		hash[*hash_size].from_num = hash[ii1].from_num;
		hash[*hash_size].type = hash[ii1].type;
		hash[*hash_size].node_num = hash[ii1].node_num;
		hash[*hash_size].tmp = hash[ii1].tmp;
		for(ii2=0; ii2<hash[ii1].node_num; ii2++){
			hash[*hash_size].node[ii2] = hash[ii1].node[ii2];
		}
		for(ii2=0; ii2<hash[*hash_size].from_num; ii2++){
			if(hash[ii1].from[ii2] < 0) break;
			hash[*hash_size].from[ii2] = hash[ii1].from[ii2];
			hash[*hash_size].from_f[ii2] = hash[ii1].from_f[ii2];
		}
		(*hash_size) ++;
	}

	hash = mm_realloc(hash, (*hash_size) * sizeof(ELEMENT_FACE) );
	if(hash == NULL) return(ALLOC_ERROR_RC);

	return(NORMAL_RC);
}


static RC
set_hash_face (int hash_size, ELEMENT_FACE *hash, int key_size,
               const int *keys, int label, ELEM_TYPE type, int number)
{
	int index;

	RC_NEG_CHK( index = elem_hash_index(hash_size, hash, key_size, keys) );

	hash[index].from_num ++;
	if( (hash[index].from_num > 2)||(hash[index].from_num < 0) ){
		return(OVERFLOW_ERROR_RC);
	}	
	hash[index].from[hash[index].from_num - 1] = label;
	hash[index].from_f[hash[index].from_num - 1] = number;
	hash[index].type = type;

	return(NORMAL_RC);
}


static RC
set_hash_tri1 (int hash_size, ELEMENT_FACE *hash, ELEMENT elem)
{
	int keys[MAX_NODE_LIGHT];

	keys[0] = elem.node[0];
	keys[1] = elem.node[1];
	RC_TRY( set_hash_face(hash_size, hash, 2, keys, elem.label, ELEM_LINE1,
	                                                            -1) );

	keys[0] = elem.node[1];
	keys[1] = elem.node[2];
	RC_TRY( set_hash_face(hash_size, hash, 2, keys, elem.label, ELEM_LINE1,
	                                                            -1) );

	keys[0] = elem.node[2];
	keys[1] = elem.node[0];
	RC_TRY( set_hash_face(hash_size, hash, 2, keys, elem.label, ELEM_LINE1,
	                                                            -1) );

	return(NORMAL_RC);
}


static RC
set_hash_tri2 (int hash_size, ELEMENT_FACE *hash, ELEMENT elem)
{
	int keys[MAX_NODE_LIGHT];

	keys[0] = elem.node[0];
	keys[1] = elem.node[1];
	keys[2] = elem.node[5];
	RC_TRY( set_hash_face(hash_size, hash, 3, keys, elem.label, ELEM_LINE2,
	                                                            -1) );

	keys[0] = elem.node[1];
	keys[1] = elem.node[2];
	keys[2] = elem.node[3];
	RC_TRY( set_hash_face(hash_size, hash, 3, keys, elem.label, ELEM_LINE2,
	                                                            -1) );

	keys[0] = elem.node[2];
	keys[1] = elem.node[0];
	keys[2] = elem.node[4];
	RC_TRY( set_hash_face(hash_size, hash, 3, keys, elem.label, ELEM_LINE2,
	                                                            -1) );

	return(NORMAL_RC);
}


static RC
set_hash_quad1 (int hash_size, ELEMENT_FACE *hash, ELEMENT elem)
{
	int keys[MAX_NODE_LIGHT];

	keys[0] = elem.node[0];
	keys[1] = elem.node[1];
	RC_TRY( set_hash_face(hash_size, hash, 2, keys, elem.label, ELEM_LINE1,
	                                                            -1) );

	keys[0] = elem.node[1];
	keys[1] = elem.node[2];
	RC_TRY( set_hash_face(hash_size, hash, 2, keys, elem.label, ELEM_LINE1,
	                                                            -1) );

	keys[0] = elem.node[2];
	keys[1] = elem.node[3];
	RC_TRY( set_hash_face(hash_size, hash, 2, keys, elem.label, ELEM_LINE1,
	                                                            -1) );

	keys[0] = elem.node[3];
	keys[1] = elem.node[0];
	RC_TRY( set_hash_face(hash_size, hash, 2, keys, elem.label, ELEM_LINE1,
	                                                            -1) );

	return(NORMAL_RC);
}


static RC
set_hash_quad2 (int hash_size, ELEMENT_FACE *hash, ELEMENT elem)
{
	int keys[MAX_NODE_LIGHT];

	keys[0] = elem.node[0];
	keys[1] = elem.node[1];
	keys[2] = elem.node[4];
	RC_TRY( set_hash_face(hash_size, hash, 3, keys, elem.label, ELEM_LINE2,
	                                                            -1) );

	keys[0] = elem.node[1];
	keys[1] = elem.node[2];
	keys[2] = elem.node[5];
	RC_TRY( set_hash_face(hash_size, hash, 3, keys, elem.label, ELEM_LINE2,
	                                                            -1) );

	keys[0] = elem.node[2];
	keys[1] = elem.node[3];
	keys[2] = elem.node[6];
	RC_TRY( set_hash_face(hash_size, hash, 3, keys, elem.label, ELEM_LINE2,
	                                                            -1) );

	keys[0] = elem.node[3];
	keys[1] = elem.node[0];
	keys[2] = elem.node[7];
	RC_TRY( set_hash_face(hash_size, hash, 3, keys, elem.label, ELEM_LINE2,
	                                                            -1) );

	return(NORMAL_RC);
}


static RC
set_hash_tetra1 (int hash_size, ELEMENT_FACE *hash, ELEMENT elem)
{
	int keys[MAX_NODE_LIGHT];

	keys[0] = elem.node[0];
	keys[1] = elem.node[3];
	keys[2] = elem.node[1];
	RC_TRY( set_hash_face(hash_size, hash, 3, keys, elem.label, ELEM_TRI1,
	                                                            -1) );

	keys[0] = elem.node[1];
	keys[1] = elem.node[3];
	keys[2] = elem.node[2];
	RC_TRY( set_hash_face(hash_size, hash, 3, keys, elem.label, ELEM_TRI1,
	                                                            -1) );

	keys[0] = elem.node[0];
	keys[1] = elem.node[1];
	keys[2] = elem.node[2];
	RC_TRY( set_hash_face(hash_size, hash, 3, keys, elem.label, ELEM_TRI1,
	                                                            -1) );

	keys[0] = elem.node[2];
	keys[1] = elem.node[3];
	keys[2] = elem.node[0];
	RC_TRY( set_hash_face(hash_size, hash, 3, keys, elem.label, ELEM_TRI1,
	                                                            -1) );
	
	return(NORMAL_RC);
}


static RC
set_hash_tetra2 (int hash_size, ELEMENT_FACE *hash, ELEMENT elem)
{
	int keys[MAX_NODE_LIGHT];

	keys[0] = elem.node[0];
	keys[1] = elem.node[3];
	keys[2] = elem.node[1];
	keys[3] = elem.node[8];
	keys[4] = elem.node[6];
	keys[5] = elem.node[7];
	RC_TRY( set_hash_face(hash_size, hash, 6, keys, elem.label, ELEM_TRI2,
	                                                            22) );

	keys[0] = elem.node[1];
	keys[1] = elem.node[3];
	keys[2] = elem.node[2];
	keys[3] = elem.node[9];
	keys[4] = elem.node[4];
	keys[5] = elem.node[8];
	RC_TRY( set_hash_face(hash_size, hash, 6, keys, elem.label, ELEM_TRI2,
	                                                            25) );

	keys[0] = elem.node[0];
	keys[1] = elem.node[1];
	keys[2] = elem.node[2];
	keys[3] = elem.node[4];
	keys[4] = elem.node[5];
	keys[5] = elem.node[6];
	RC_TRY( set_hash_face(hash_size, hash, 6, keys, elem.label, ELEM_TRI2,
	                                                            28) );

	keys[0] = elem.node[2];
	keys[1] = elem.node[3];
	keys[2] = elem.node[0];
	keys[3] = elem.node[7];
	keys[4] = elem.node[5];
	keys[5] = elem.node[9];
	RC_TRY( set_hash_face(hash_size, hash, 6, keys, elem.label, ELEM_TRI2,
	                                                            31) );

	return(NORMAL_RC);
}


static RC
set_hash_penta1 (int hash_size, ELEMENT_FACE *hash, ELEMENT elem)
{
	int keys[MAX_NODE_LIGHT];

	keys[0] = elem.node[1];
	keys[1] = elem.node[0];
	keys[2] = elem.node[2];
	RC_TRY( set_hash_face(hash_size, hash, 3, keys, elem.label, ELEM_TRI1,
	                                                            -1) );

	keys[0] = elem.node[3];
	keys[1] = elem.node[4];
	keys[2] = elem.node[5];
	RC_TRY( set_hash_face(hash_size, hash, 3, keys, elem.label, ELEM_TRI1,
	                                                            -1) );

	keys[0] = elem.node[1];
	keys[1] = elem.node[2];
	keys[2] = elem.node[5];
	keys[3] = elem.node[4];
	RC_TRY( set_hash_face(hash_size, hash, 4, keys, elem.label, ELEM_QUAD1,
	                                                            15) );

	keys[0] = elem.node[0];
	keys[1] = elem.node[3];
	keys[2] = elem.node[5];
	keys[3] = elem.node[2];
	RC_TRY( set_hash_face(hash_size, hash, 4, keys, elem.label, ELEM_QUAD1,
	                                                            16) );

	keys[0] = elem.node[0];
	keys[1] = elem.node[1];
	keys[2] = elem.node[4];
	keys[3] = elem.node[3];
	RC_TRY( set_hash_face(hash_size, hash, 4, keys, elem.label, ELEM_QUAD1,
	                                                            17) );

	return(NORMAL_RC);
}


static RC
set_hash_penta2 (int hash_size, ELEMENT_FACE *hash, ELEMENT elem)
{
	int keys[MAX_NODE_LIGHT];

	keys[0] = elem.node[1];
	keys[1] = elem.node[0];
	keys[2] = elem.node[2];
	keys[3] = elem.node[7];
	keys[4] = elem.node[6];
	keys[5] = elem.node[8];
	RC_TRY( set_hash_face(hash_size, hash, 6, keys, elem.label, ELEM_TRI2,
	                                                            33) );

	keys[0] = elem.node[3];
	keys[1] = elem.node[4];
	keys[2] = elem.node[5];
	keys[3] = elem.node[9];
	keys[4] = elem.node[10];
	keys[5] = elem.node[11];
	RC_TRY( set_hash_face(hash_size, hash, 6, keys, elem.label, ELEM_TRI2,
	                                                            36) );

	keys[0] = elem.node[0];
	keys[1] = elem.node[1];
	keys[2] = elem.node[4];
	keys[3] = elem.node[3];
	keys[4] = elem.node[8];
	keys[5] = elem.node[13];
	keys[6] = elem.node[11];
	keys[7] = elem.node[12];
	RC_TRY( set_hash_face(hash_size, hash, 8, keys, elem.label, ELEM_QUAD2,
	                                                            39) );

	keys[0] = elem.node[1];
	keys[1] = elem.node[2];
	keys[2] = elem.node[5];
	keys[3] = elem.node[4];
	keys[4] = elem.node[6];
	keys[5] = elem.node[14];
	keys[6] = elem.node[9];
	keys[7] = elem.node[13];
	RC_TRY( set_hash_face(hash_size, hash, 8, keys, elem.label, ELEM_QUAD2,
	                                                            44) );

	keys[0] = elem.node[2];
	keys[1] = elem.node[0];
	keys[2] = elem.node[3];
	keys[3] = elem.node[5];
	keys[4] = elem.node[7];
	keys[5] = elem.node[12];
	keys[6] = elem.node[10];
	keys[7] = elem.node[14];
	RC_TRY( set_hash_face(hash_size, hash, 8, keys, elem.label, ELEM_QUAD2,
	                                                            49) );

	return(NORMAL_RC);
}


static RC
set_hash_hexa1 (int hash_size, ELEMENT_FACE *hash, ELEMENT elem)
{
	int keys[MAX_NODE_LIGHT];

	keys[0] = elem.node[0];
	keys[1] = elem.node[3];
	keys[2] = elem.node[2];
	keys[3] = elem.node[1];
	RC_TRY( set_hash_face(hash_size, hash, 4, keys, elem.label, ELEM_QUAD1,
	                                                            20) );

	keys[0] = elem.node[0];
	keys[1] = elem.node[1];
	keys[2] = elem.node[5];
	keys[3] = elem.node[4];
	RC_TRY( set_hash_face(hash_size, hash, 4, keys, elem.label, ELEM_QUAD1,
	                                                            21) );

	keys[0] = elem.node[1];
	keys[1] = elem.node[2];
	keys[2] = elem.node[6];
	keys[3] = elem.node[5];
	RC_TRY( set_hash_face(hash_size, hash, 4, keys, elem.label, ELEM_QUAD1,
	                                                            22) );

	keys[0] = elem.node[7];
	keys[1] = elem.node[4];
	keys[2] = elem.node[5];
	keys[3] = elem.node[6];
	RC_TRY( set_hash_face(hash_size, hash, 4, keys, elem.label, ELEM_QUAD1,
	                                                            23) );

	keys[0] = elem.node[3];
	keys[1] = elem.node[0];
	keys[2] = elem.node[4];
	keys[3] = elem.node[7];
	RC_TRY( set_hash_face(hash_size, hash, 4, keys, elem.label, ELEM_QUAD1,
	                                                            24) );

	keys[0] = elem.node[2];
	keys[1] = elem.node[3];
	keys[2] = elem.node[7];
	keys[3] = elem.node[6];
	RC_TRY( set_hash_face(hash_size, hash, 4, keys, elem.label, ELEM_QUAD1,
	                                                            25) );

	return(NORMAL_RC);
}


static RC
set_hash_hexa2 (int hash_size, ELEMENT_FACE *hash, ELEMENT elem)
{
	int keys[MAX_NODE_LIGHT];

	keys[0] = elem.node[0];
	keys[1] = elem.node[3];
	keys[2] = elem.node[2];
	keys[3] = elem.node[1];
	keys[4] = elem.node[11];
	keys[5] = elem.node[10];
	keys[6] = elem.node[9];
	keys[7] = elem.node[8];
	RC_TRY( set_hash_face(hash_size, hash, 8, keys, elem.label, ELEM_QUAD2,
	                                                            44) );

	keys[0] = elem.node[0];
	keys[1] = elem.node[1];
	keys[2] = elem.node[5];
	keys[3] = elem.node[4];
	keys[4] = elem.node[8];
	keys[5] = elem.node[17];
	keys[6] = elem.node[12];
	keys[7] = elem.node[16];
	RC_TRY( set_hash_face(hash_size, hash, 8, keys, elem.label, ELEM_QUAD2,
	                                                            49) );

	keys[0] = elem.node[1];
	keys[1] = elem.node[2];
	keys[2] = elem.node[6];
	keys[3] = elem.node[5];
	keys[4] = elem.node[9];
	keys[5] = elem.node[18];
	keys[6] = elem.node[13];
	keys[7] = elem.node[17];
	RC_TRY( set_hash_face(hash_size, hash, 8, keys, elem.label, ELEM_QUAD2,
	                                                            54) );

	keys[0] = elem.node[7];
	keys[1] = elem.node[4];
	keys[2] = elem.node[5];
	keys[3] = elem.node[6];
	keys[4] = elem.node[15];
	keys[5] = elem.node[12];
	keys[6] = elem.node[13];
	keys[7] = elem.node[14];
	RC_TRY( set_hash_face(hash_size, hash, 8, keys, elem.label, ELEM_QUAD2,
	                                                            59) );

	keys[0] = elem.node[3];
	keys[1] = elem.node[0];
	keys[2] = elem.node[4];
	keys[3] = elem.node[7];
	keys[4] = elem.node[11];
	keys[5] = elem.node[16];
	keys[6] = elem.node[15];
	keys[7] = elem.node[19];
	RC_TRY( set_hash_face(hash_size, hash, 8, keys, elem.label, ELEM_QUAD2,
	                                                            64) );

	keys[0] = elem.node[2];
	keys[1] = elem.node[3];
	keys[2] = elem.node[7];
	keys[3] = elem.node[6];
	keys[4] = elem.node[10];
	keys[5] = elem.node[19];
	keys[6] = elem.node[14];
	keys[7] = elem.node[18];
	RC_TRY( set_hash_face(hash_size, hash, 8, keys, elem.label, ELEM_QUAD2,
	                                                            69) );

	return(NORMAL_RC);
}


static RC
set_cut_elem (ELEMENT_ARRAY *elem, ELEMENT_NODE_BOX box,
              int index_size, const int *index)
{
	int ii1;

	RC_TRY( realloc_element_array(elem) );
	elem->array[elem->size].label = elem->size + 1;
	elem->array[elem->size].type = box.type;
	elem->array[elem->size].material = box.material;
	elem->array[elem->size].physical = box.physical;
	elem->array[elem->size].node_num = index_size;
	for(ii1=0; ii1<elem->array[elem->size].node_num; ii1++){
		elem->array[elem->size].node[ii1] = box.node[ index[ii1] ];
	}
	elem->size ++;

	return(NORMAL_RC);
}


static RC
input_cut_elem_line1 (ELEMENT_ARRAY *elem, ELEMENT_NODE_BOX box)
{
	int index[2];

	index[0] = 0;
	index[1] = 2;
	RC_TRY( set_cut_elem(elem, box, 2, index) );

	index[0] = 2;
	index[1] = 1;
	RC_TRY( set_cut_elem(elem, box, 2, index) );

	return(NORMAL_RC);
}


static RC
input_cut_elem_line2 (ELEMENT_ARRAY *elem, ELEMENT_NODE_BOX box)
{
	int index[3];

	index[0] = 0;
	index[1] = 2;
	index[2] = 3;
	RC_TRY( set_cut_elem(elem, box, 3, index) );

	index[0] = 2;
	index[1] = 1;
	index[2] = 4;
	RC_TRY( set_cut_elem(elem, box, 3, index) );

	return(NORMAL_RC);
}


static RC
input_cut_elem_tri1 (ELEMENT_ARRAY *elem, ELEMENT_NODE_BOX box)
{
	int index[3];

	index[0] = 0;
	index[1] = 3;
	index[2] = 5;
	RC_TRY( set_cut_elem(elem, box, 3, index) );

	index[0] = 1;
	index[1] = 4;
	index[2] = 3;
	RC_TRY( set_cut_elem(elem, box, 3, index) );

	index[0] = 2;
	index[1] = 5;
	index[2] = 4;
	RC_TRY( set_cut_elem(elem, box, 3, index) );

	index[0] = 3;
	index[1] = 4;
	index[2] = 5;
	RC_TRY( set_cut_elem(elem, box, 3, index) );

	return(NORMAL_RC);
}


static RC
input_cut_elem_tri2 (ELEMENT_ARRAY *elem, ELEMENT_NODE_BOX box)
{
	int index[6];

	index[0] = 0;
	index[1] = 5;
	index[2] = 4;
	index[3] = 12;
	index[4] = 11;
	index[5] = 6;
	RC_TRY( set_cut_elem(elem, box, 6, index) );

	index[0] = 1;
	index[1] = 3;
	index[2] = 5;
	index[3] = 13;
	index[4] = 7;
	index[5] = 8;
	RC_TRY( set_cut_elem(elem, box, 6, index) );

	index[0] = 2;
	index[1] = 4;
	index[2] = 3;
	index[3] = 14;
	index[4] = 9;
	index[5] = 10;
	RC_TRY( set_cut_elem(elem, box, 6, index) );

	index[0] = 3;
	index[1] = 4;
	index[2] = 5;
	index[3] = 12;
	index[4] = 13; 
	index[5] = 14;
	RC_TRY( set_cut_elem(elem, box, 6, index) );

	return(NORMAL_RC);
}


static RC
input_cut_elem_quad1 (ELEMENT_ARRAY *elem, ELEMENT_NODE_BOX box)
{
	int index[4];

	index[0] = 0;
	index[1] = 4;
	index[2] = 8;
	index[3] = 7;
	RC_TRY( set_cut_elem(elem, box, 4, index) );

	index[0] = 1;
	index[1] = 5;
	index[2] = 8;
	index[3] = 4;
	RC_TRY( set_cut_elem(elem, box, 4, index) );

	index[0] = 2;
	index[1] = 6;
	index[2] = 8;
	index[3] = 5;
	RC_TRY( set_cut_elem(elem, box, 4, index) );

	index[0] = 3;
	index[1] = 7;
	index[2] = 8;
	index[3] = 6;
	RC_TRY( set_cut_elem(elem, box, 4, index) );

	return(NORMAL_RC);
}


static RC
input_cut_elem_quad2 (ELEMENT_ARRAY *elem, ELEMENT_NODE_BOX box)
{
	int index[8];

	index[0] = 0;
	index[1] = 4;
	index[2] = 16;
	index[3] = 7;
	index[4] = 8;
	index[5] = 19;
	index[6] = 17;
	index[7] = 16;
	RC_TRY( set_cut_elem(elem, box, 8, index) );

	index[0] = 1;
	index[1] = 5;
	index[2] = 16;
	index[3] = 4;
	index[4] = 9;
	index[5] = 10;
	index[6] = 18;
	index[7] = 19;
	RC_TRY( set_cut_elem(elem, box, 8, index) );

	index[0] = 2;
	index[1] = 6;
	index[2] = 16;
	index[3] = 5;
	index[4] = 12;
	index[5] = 20;
	index[6] = 18;
	index[7] = 11;
	RC_TRY( set_cut_elem(elem, box, 8, index) );

	index[0] = 3;
	index[1] = 7;
	index[2] = 16;
	index[3] = 6;
	index[4] = 14;
	index[5] = 17;
	index[6] = 20;
	index[7] = 13;
	RC_TRY( set_cut_elem(elem, box, 8, index) );

	return(NORMAL_RC);
}


static RC
input_cut_elem_tetra1 (ELEMENT_ARRAY *elem, ELEMENT_NODE_BOX box)
{
	int index[4];

	index[0] = 0;
	index[1] = 6;
	index[2] = 5;
	index[3] = 7;
	RC_TRY( set_cut_elem(elem, box, 4, index) );

	index[0] = 1;
	index[1] = 4;
	index[2] = 6;
	index[3] = 8;
	RC_TRY( set_cut_elem(elem, box, 4, index) );

	index[0] = 2;
	index[1] = 5;
	index[2] = 4;
	index[3] = 9;
	RC_TRY( set_cut_elem(elem, box, 4, index) );

	index[0] = 7;
	index[1] = 8;
	index[2] = 9;
	index[3] = 3;
	RC_TRY( set_cut_elem(elem, box, 4, index) );

	index[0] = 7;
	index[1] = 4;
	index[2] = 9;
	index[3] = 8;
	RC_TRY( set_cut_elem(elem, box, 4, index) );

	index[0] = 7;
	index[1] = 9;
	index[2] = 4;
	index[3] = 5;
	RC_TRY( set_cut_elem(elem, box, 4, index) );

	index[0] = 7;
	index[1] = 6;
	index[2] = 4;
	index[3] = 8;
	RC_TRY( set_cut_elem(elem, box, 4, index) );

	index[0] = 7;
	index[1] = 4;
	index[2] = 6;
	index[3] = 5;
	RC_TRY( set_cut_elem(elem, box, 4, index) );

	return(NORMAL_RC);
}


static RC
input_cut_elem_tetra2 (ELEMENT_ARRAY *elem, ELEMENT_NODE_BOX box)
{
	int index[10];

	index[0] = 0;
	index[1] = 6;
	index[2] = 5;
	index[3] = 7;
	index[4] = 28;
	index[5] = 15;
	index[6] = 10;
	index[7] = 16;
	index[8] = 22;
	index[9] = 33;
	RC_TRY( set_cut_elem(elem, box, 10, index) );

	index[0] = 1;
	index[1] = 4;
	index[2] = 6;
	index[3] = 8;
	index[4] = 29;
	index[5] = 11;
	index[6] = 12;
	index[7] = 18;
	index[8] = 25;
	index[9] = 24;
	RC_TRY( set_cut_elem(elem, box, 10, index) );

	index[0] = 2;
	index[1] = 5;
	index[2] = 4;
	index[3] = 9;
	index[4] = 30;
	index[5] = 13;
	index[6] = 14;
	index[7] = 20;
	index[8] = 31;
	index[9] = 27;
	RC_TRY( set_cut_elem(elem, box, 10, index) );

	index[0] = 7;
	index[1] = 8;
	index[2] = 9;
	index[3] = 3;
	index[4] = 26;
	index[5] = 32;
	index[6] = 23;
	index[7] = 17;
	index[8] = 19;
	index[9] = 21;
	RC_TRY( set_cut_elem(elem, box, 10, index) );

	index[0] = 7;
	index[1] = 4;
	index[2] = 9;
	index[3] = 8;
	index[4] = 27;
	index[5] = 32;
	index[6] = 34;
	index[7] = 23;
	index[8] = 25;
	index[9] = 26;
	RC_TRY( set_cut_elem(elem, box, 10, index) );

	index[0] = 7;
	index[1] = 9;
	index[2] = 4;
	index[3] = 5;
	index[4] = 27;
	index[5] = 34;
	index[6] = 32;
	index[7] = 33;
	index[8] = 31;
	index[9] = 30;
	RC_TRY( set_cut_elem(elem, box, 10, index) );

	index[0] = 7;
	index[1] = 6;
	index[2] = 4;
	index[3] = 8;
	index[4] = 29;
	index[5] = 34;
	index[6] = 22;
	index[7] = 23;
	index[8] = 24;
	index[9] = 25;
	RC_TRY( set_cut_elem(elem, box, 10, index) );

	index[0] = 7;
	index[1] = 4;
	index[2] = 6;
	index[3] = 5;
	index[4] = 29;
	index[5] = 22;
	index[6] = 34;
	index[7] = 33;
	index[8] = 30;
	index[9] = 28;
	RC_TRY( set_cut_elem(elem, box, 10, index) );

	return(NORMAL_RC);
}


static RC
input_cut_elem_penta1 (ELEMENT_ARRAY *elem, ELEMENT_NODE_BOX box)
{
	int index[6];

	index[0] = 0;
	index[1] = 8;
	index[2] = 7;
	index[3] = 12;
	index[4] = 17;
	index[5] = 16;
	RC_TRY( set_cut_elem(elem, box, 6, index) );

	index[0] = 1;
	index[1] = 6;
	index[2] = 8;
	index[3] = 13;
	index[4] = 15;
	index[5] = 17;
	RC_TRY( set_cut_elem(elem, box, 6, index) );

	index[0] = 2;
	index[1] = 7;
	index[2] = 6;
	index[3] = 14;
	index[4] = 16;
	index[5] = 15;
	RC_TRY( set_cut_elem(elem, box, 6, index) );

	index[0] = 6;
	index[1] = 7;
	index[2] = 8;
	index[3] = 15;
	index[4] = 16;
	index[5] = 17;
	RC_TRY( set_cut_elem(elem, box, 6, index) );

	index[0] = 12;
	index[1] = 17;
	index[2] = 16;
	index[3] = 3;
	index[4] = 11;
	index[5] = 10;
	RC_TRY( set_cut_elem(elem, box, 6, index) );

	index[0] = 13;
	index[1] = 15;
	index[2] = 17;
	index[3] = 4;
	index[4] = 9;
	index[5] = 11;
	RC_TRY( set_cut_elem(elem, box, 6, index) );

	index[0] = 14;
	index[1] = 16;
	index[2] = 15;
	index[3] = 5;
	index[4] = 10;
	index[5] = 9;
	RC_TRY( set_cut_elem(elem, box, 6, index) );

	index[0] = 15;
	index[1] = 16;
	index[2] = 17;
	index[3] = 9;
	index[4] = 10;
	index[5] = 11;
	RC_TRY( set_cut_elem(elem, box, 6, index) );

	return(NORMAL_RC);
}


static RC
input_cut_elem_penta2 (ELEMENT_ARRAY *elem, ELEMENT_NODE_BOX box)
{
	int index[15];

	index[0] = 0;
	index[1] = 8;
	index[2] = 7;
	index[3] = 12;
	index[4] = 39;
	index[5] = 49;
	index[6] = 34;
	index[7] = 20;
	index[8] = 15;
	index[9] = 54;
	index[10] = 51;
	index[11] = 40;
	index[12] = 27;
	index[13] = 42;
	index[14] = 52;
	RC_TRY( set_cut_elem(elem, box, 15, index) );

	index[0] = 1;
	index[1] = 6;
	index[2] = 8;
	index[3] = 13;
	index[4] = 44;
	index[5] = 39;
	index[6] = 33;
	index[7] = 16;
	index[8] = 17;
	index[9] = 55;
	index[10] = 41;
	index[11] = 45;
	index[12] = 29;
	index[13] = 47;
	index[14] = 42;
	RC_TRY( set_cut_elem(elem, box, 15, index) );

	index[0] = 2;
	index[1] = 7;
	index[2] = 6;
	index[3] = 14;
	index[4] = 49;
	index[5] = 44;
	index[6] = 35;
	index[7] = 18;
	index[8] = 19;
	index[9] = 56;
	index[10] = 46;
	index[11] = 50;
	index[12] = 31;
	index[13] = 52;
	index[14] = 47;
	RC_TRY( set_cut_elem(elem, box, 15, index) );

	index[0] = 6;
	index[1] = 7;
	index[2] = 8;
	index[3] = 44;
	index[4] = 49;
	index[5] = 39;
	index[6] = 34;
	index[7] = 33;
	index[8] = 35;
	index[9] = 54;
	index[10] = 55;
	index[11] = 56;
	index[12] = 47;
	index[13] = 52;
	index[14] = 42;
	RC_TRY( set_cut_elem(elem, box, 15, index) );

	index[0] = 12;
	index[1] = 39;
	index[2] = 49;
	index[3] = 3;
	index[4] = 11;
	index[5] = 10;
	index[6] = 54;
	index[7] = 51;
	index[8] = 40;
	index[9] = 36;
	index[10] = 26;
	index[11] = 21;
	index[12] = 28;
	index[13] = 43;
	index[14] = 53;
	RC_TRY( set_cut_elem(elem, box, 15, index) );

	index[0] = 13;
	index[1] = 44;
	index[2] = 39;
	index[3] = 4;
	index[4] = 9;
	index[5] = 11;
	index[6] = 55;
	index[7] = 41;
	index[8] = 45;
	index[9] = 37;
	index[10] = 22;
	index[11] = 23;
	index[12] = 30;
	index[13] = 48;
	index[14] = 43;
	RC_TRY( set_cut_elem(elem, box, 15, index) );

	index[0] = 14;
	index[1] = 49;
	index[2] = 44;
	index[3] = 5;
	index[4] = 10;
	index[5] = 9;
	index[6] = 56;
	index[7] = 46;
	index[8] = 50;
	index[9] = 38;
	index[10] = 24;
	index[11] = 25;
	index[12] = 32;
	index[13] = 53;
	index[14] = 48;
	RC_TRY( set_cut_elem(elem, box, 15, index) );

	index[0] = 44;
	index[1] = 49;
	index[2] = 39;
	index[3] = 9;
	index[4] = 10;
	index[5] = 11;
	index[6] = 54;
	index[7] = 55;
	index[8] = 56;
	index[9] = 36;
	index[10] = 37;
	index[11] = 38;
	index[12] = 48;
	index[13] = 53;
	index[14] = 43;
	RC_TRY( set_cut_elem(elem, box, 15, index) );

	return(NORMAL_RC);
}


static RC
input_cut_elem_hexa1 (ELEMENT_ARRAY *elem, ELEMENT_NODE_BOX box)
{
	int index[8];

	index[0] = 0;
	index[1] = 8;
	index[2] = 20;
	index[3] = 11;
	index[4] = 16;
	index[5] = 21;
	index[6] = 26;
	index[7] = 24;
	RC_TRY( set_cut_elem(elem, box, 8, index) );

	index[0] = 8;
	index[1] = 1;
	index[2] = 9;
	index[3] = 20;
	index[4] = 21;
	index[5] = 17;
	index[6] = 22;
	index[7] = 26;
	RC_TRY( set_cut_elem(elem, box, 8, index) );

	index[0] = 20;
	index[1] = 9;
	index[2] = 2;
	index[3] = 10;
	index[4] = 26;
	index[5] = 22;
	index[6] = 18;
	index[7] = 25;
	RC_TRY( set_cut_elem(elem, box, 8, index) );

	index[0] = 11;
	index[1] = 20;
	index[2] = 10;
	index[3] = 3;
	index[4] = 24;
	index[5] = 26;
	index[6] = 25;
	index[7] = 19;
	RC_TRY( set_cut_elem(elem, box, 8, index) );

	index[0] = 16;
	index[1] = 21;
	index[2] = 26;
	index[3] = 24;
	index[4] = 4;
	index[5] = 12;
	index[6] = 23;
	index[7] = 15;
	RC_TRY( set_cut_elem(elem, box, 8, index) );

	index[0] = 21;
	index[1] = 17;
	index[2] = 22;
	index[3] = 26;
	index[4] = 12;
	index[5] = 5;
	index[6] = 13;
	index[7] = 23;
	RC_TRY( set_cut_elem(elem, box, 8, index) );

	index[0] = 26;
	index[1] = 22;
	index[2] = 18;
	index[3] = 25;
	index[4] = 23;
	index[5] = 13;
	index[6] = 6;
	index[7] = 14;
	RC_TRY( set_cut_elem(elem, box, 8, index) );

	index[0] = 24;
	index[1] = 26;
	index[2] = 25;
	index[3] = 19;
	index[4] = 15;
	index[5] = 23;
	index[6] = 14;
	index[7] = 7;
	RC_TRY( set_cut_elem(elem, box, 8, index) );

	return(NORMAL_RC);
}


static RC
input_cut_elem_hexa2 (ELEMENT_ARRAY *elem, ELEMENT_NODE_BOX box)
{
	int index[20];
	
	index[0] = 0;
	index[1] = 8;
	index[2] = 44;
	index[3] = 11;
	index[4] = 16;
	index[5] = 49;
	index[6] = 74;
	index[7] = 64;
	index[8] = 20;
	index[9] = 45;
	index[10] = 47;
	index[11] = 27;
	index[12] = 50;
	index[13] = 77;
	index[14] = 75;
	index[15] = 66;
	index[16] = 36;
	index[17] = 52;
	index[18] = 79;
	index[19] = 67;
	RC_TRY( set_cut_elem(elem, box, 20, index) );
	
	index[0] = 8;
	index[1] = 1;
	index[2] = 9;
	index[3] = 44;
	index[4] = 49;
	index[5] = 17;
	index[6] = 54;
	index[7] = 74;
	index[8] = 21;
	index[9] = 22;
	index[10] = 48;
	index[11] = 45;
	index[12] = 51;
	index[13] = 55;
	index[14] = 76;
	index[15] = 77;
	index[16] = 52;
	index[17] = 38;
	index[18] = 57;
	index[19] = 79;
	RC_TRY( set_cut_elem(elem, box, 20, index) );
	
	index[0] = 44;
	index[1] = 9;
	index[2] = 2;
	index[3] = 10;
	index[4] = 74;
	index[5] = 54;
	index[6] = 18;
	index[7] = 69;
	index[8] = 48;
	index[9] = 23;
	index[10] = 24;
	index[11] = 46;
	index[12] = 76;
	index[13] = 56;
	index[14] = 70;
	index[15] = 78;
	index[16] = 79;
	index[17] = 57;
	index[18] = 40;
	index[19] = 72;
	RC_TRY( set_cut_elem(elem, box, 20, index) );
	
	index[0] = 11;
	index[1] = 44;
	index[2] = 10;
	index[3] = 3;
	index[4] = 64;
	index[5] = 74;
	index[6] = 69;
	index[7] = 19;
	index[8] = 47;
	index[9] = 46;
	index[10] = 25;
	index[11] = 26;
	index[12] = 75;
	index[13] = 78;
	index[14] = 71;
	index[15] = 65;
	index[16] = 67;
	index[17] = 79;
	index[18] = 72;
	index[19] = 42;
	RC_TRY( set_cut_elem(elem, box, 20, index) );
	
	index[0] = 16;
	index[1] = 49;
	index[2] = 74;
	index[3] = 64;
	index[4] = 4;
	index[5] = 12;
	index[6] = 59;
	index[7] = 15;
	index[8] = 50;
	index[9] = 77;
	index[10] = 75;
	index[11] = 66;
	index[12] = 28;
	index[13] = 61;
	index[14] = 62;
	index[15] = 35;
	index[16] = 37;
	index[17] = 53;
	index[18] = 80;
	index[19] = 68;
	RC_TRY( set_cut_elem(elem, box, 20, index) );
	
	index[0] = 49;
	index[1] = 17;
	index[2] = 54;
	index[3] = 74;
	index[4] = 12;
	index[5] = 5;
	index[6] = 13;
	index[7] = 59;
	index[8] = 51;
	index[9] = 55;
	index[10] = 76;
	index[11] = 77;
	index[12] = 29;
	index[13] = 30;
	index[14] = 63;
	index[15] = 61;
	index[16] = 53;
	index[17] = 39;
	index[18] = 58;
	index[19] = 80;
	RC_TRY( set_cut_elem(elem, box, 20, index) );
	
	index[0] = 74;
	index[1] = 54;
	index[2] = 18;
	index[3] = 69;
	index[4] = 59;
	index[5] = 13;
	index[6] = 6;
	index[7] = 14;
	index[8] = 76;
	index[9] = 56;
	index[10] = 70;
	index[11] = 78;
	index[12] = 63;
	index[13] = 31;
	index[14] = 32;
	index[15] = 60;
	index[16] = 80;
	index[17] = 58;
	index[18] = 41;
	index[19] = 73;
	RC_TRY( set_cut_elem(elem, box, 20, index) );
	
	index[0] = 64;
	index[1] = 74;
	index[2] = 69;
	index[3] = 19;
	index[4] = 15;
	index[5] = 59;
	index[6] = 14;
	index[7] = 7;
	index[8] = 75;
	index[9] = 78;
	index[10] = 71;
	index[11] = 65;
	index[12] = 62;
	index[13] = 60;
	index[14] = 33;
	index[15] = 34;
	index[16] = 68;
	index[17] = 80;
	index[18] = 73;
	index[19] = 43;
	RC_TRY( set_cut_elem(elem, box, 20, index) );
	
	return(NORMAL_RC);
}


static RC
extract_edge (ELEMENT_ARRAY elem, int *edge_size, ELEMENT_EDGE **hash)
{
	int hash_size;
	int ii1, ii2;

	/* calculate size of edge */
	if(elem.size <= 0) return(ARG_ERROR_RC);
	hash_size = cal_edge_size(elem);

	*hash = mm_alloc( hash_size * sizeof(ELEMENT_EDGE) );
	if(*hash == NULL) return(ALLOC_ERROR_RC);

	/* initialize hash */
	for(ii1=0; ii1<hash_size; ii1++){
		(*hash)[ii1].from_num = -1;
		for(ii2=0; ii2<EDGE_MAX_FROM; ii2++){
			(*hash)[ii1].from[ii2] = -1;
			(*hash)[ii1].from_e[ii2] = -1;
		}
	}

	/* initialize i_random */
	init_i_random(hash_size);

	for(ii1=0; ii1<elem.size; ii1++){
		if(elem.array[ii1].label < 0) continue;
		switch(elem.array[ii1].type){
		case ELEM_TRI1:
			RC_TRY( set_edge_tri1(hash_size, *hash, elem.array[ii1]) );
			break;
		case ELEM_TRI2:
			RC_TRY( set_edge_tri2(hash_size, *hash, elem.array[ii1]) );
			break;
		case ELEM_QUAD1:
			RC_TRY( set_edge_quad1(hash_size, *hash, elem.array[ii1]) );
			break;
		case ELEM_QUAD2:
			RC_TRY( set_edge_quad2(hash_size, *hash, elem.array[ii1]) );
			break;
		case ELEM_TETRA1:
			RC_TRY( set_edge_tetra1(hash_size, *hash, elem.array[ii1]) );
			break;
		case ELEM_TETRA2:
			RC_TRY( set_edge_tetra2(hash_size, *hash, elem.array[ii1]) );
			break;
		case ELEM_PENTA1:
			RC_TRY( set_edge_penta1(hash_size, *hash, elem.array[ii1]) );
			break;
		case ELEM_PENTA2:
			RC_TRY( set_edge_penta2(hash_size, *hash, elem.array[ii1]) );
			break;
		case ELEM_HEXA1:
			RC_TRY( set_edge_hexa1(hash_size, *hash, elem.array[ii1]) );
			break;
		case ELEM_HEXA2:
			RC_TRY( set_edge_hexa2(hash_size, *hash, elem.array[ii1]) );
			break;
		default :
			return(IMPLEMENT_ERROR_RC);
		}
	}

	*edge_size = 0;
	for(ii1=0; ii1<hash_size; ii1++){
		if( (*hash)[ii1].from_num <= 0 ) continue;
		(*hash)[*edge_size].type = (*hash)[ii1].type;
		(*hash)[*edge_size].from_num = (*hash)[ii1].from_num;
		(*hash)[*edge_size].node_num = (*hash)[ii1].node_num;
		for(ii2=0; ii2<(*hash)[ii1].node_num; ii2++){
			(*hash)[*edge_size].node[ii2] = (*hash)[ii1].node[ii2];
		}
		for(ii2=0; ii2<(*hash)[ii1].from_num; ii2++){
			(*hash)[*edge_size].from[ii2] = (*hash)[ii1].from[ii2];
			(*hash)[*edge_size].from_e[ii2] = (*hash)[ii1].from_e[ii2];
		}
		(*edge_size) ++;
	}

	*hash = mm_realloc((*hash), (*edge_size) * sizeof(ELEMENT_EDGE) );
	if(*hash == NULL) return(ALLOC_ERROR_RC);

	return(NORMAL_RC);
}


static int
count_node_box_size (ELEM_TYPE type)
{
	switch(type){
	case ELEM_LINE1:
		return(3);
	case ELEM_LINE2:
		return(5);
	case ELEM_TRI1:
		return(6);
	case ELEM_TRI2:
		return(15);
	case ELEM_QUAD1:
		return(9);
	case ELEM_QUAD2:
		return(21);
	case ELEM_TETRA1:
		return(10);
	case ELEM_TETRA2:
		return(35);
	case ELEM_PENTA1:
		return(18);
	case ELEM_PENTA2:
		return(57);
	case ELEM_HEXA1:
		return(27);
	case ELEM_HEXA2:
		return(81);
	default:
		break;
	}

	return(-1);
}


static int
cal_edge_size (ELEMENT_ARRAY elem)
{
	int ii1;
	double ret = CONST_TERM;

	for(ii1=0; ii1<elem.size; ii1++){
		if(elem.array[ii1].label < 0) continue;
		switch(elem.array[ii1].type){
		case ELEM_TRI1:
		case ELEM_TRI2:
			ret += 3.0;
			break;
		case ELEM_QUAD1:
		case ELEM_QUAD2:
			ret += 4.0;
			break;
		case ELEM_TETRA1:
		case ELEM_TETRA2:
			ret += 6.0;
			break;
		case ELEM_PENTA1:
		case ELEM_PENTA2:
			ret += 9.0;
			break;
		case ELEM_HEXA1:
		case ELEM_HEXA2:
			ret += 12.0;
			break;
		default :
			ret += 12.0;
		}
	}

	/* 32 の倍数 + 1 */
	return( (32*((int)(ret/32.0) + 1) + 1) );

}


static RC
interpolate_coord_edge (NODE_ARRAY *node, ELEMENT_ARRAY elem,
                        ELEMENT_EDGE edge, ELEMENT_NODE_BOX *box,
                        int *new_label)
{
	double N[3];
	double **coord;
	VECT3D new_coord;
	int add_node_number;
	int index;
	int elem_index[2];
	int label_box[2];
	int ii1, ii2, ii3;

	RC_TRY( allocate2D(2, 1, &coord) );
	for(ii1=0; ii1<2; ii1++){
		elem_index[ii1] = -1;
		label_box[ii1] = -1;
	}

	switch(edge.type){
	case ELEM_LINE1:
	case ELEM_LINE2:
		RC_TRY( local_inside_coord(edge.type, coord, &add_node_number) );
		break;
	default:
		return(ARG_ERROR_RC);
	}

	for(ii1=0; ii1<add_node_number; ii1++){
		RC_TRY( set_N_vector(edge.type, coord[ii1], N) );
		RC_TRY( cal_interpolate_node(*node, N, edge.node_num, edge.node,
		                             &new_coord) );
		(*new_label) ++;
		RC_TRY( add_node_array(node, new_coord, *new_label) );
		label_box[ii1] = *new_label;
	}

	for(ii1=0; ii1<edge.from_num; ii1++){
		RC_NEG_CHK( index = search_element_label(elem, edge.from[ii1]) );
		if(edge.type == ELEM_LINE1){
			box[index].node[edge.from_e[ii1]] = label_box[0];
		}else{
			for(ii2=0; ii2<2; ii2++){
				for(ii3=0; ii3<elem.array[index].node_num; ii3++){
					if(edge.node[ii2] == elem.array[index].node[ii3]){
						elem_index[ii2] = ii3;
						break;
					}
				}
				RC_NEG_CHK(elem_index[ii2]);
			}
			switch( check_rotation_edge(elem.array[index].type, elem_index) ){
			case 0:
				box[index].node[edge.from_e[ii1] + 0] = label_box[0];
				box[index].node[edge.from_e[ii1] + 1] = label_box[1];
				break;
			case 1:
				box[index].node[edge.from_e[ii1] + 0] = label_box[1];
				box[index].node[edge.from_e[ii1] + 1] = label_box[0];
				break;
			default:
				return(ARG_ERROR_RC);
			}
		}
	}

	RC_TRY( free2D(2, 1, &coord) );

	return(NORMAL_RC);
}


static RC
interpolate_coord_face (NODE_ARRAY *node, ELEMENT_ARRAY elem,
                        ELEMENT_FACE face, ELEMENT_NODE_BOX *box,
                        int *new_label)
{
	double N[8];
	double **coord;
	VECT3D new_coord;
	int add_node_number;
	int index;
	int label_box[5];
	int ii1;

	RC_TRY( allocate2D(5, 3, &coord) );

	switch(face.type){
	case ELEM_TRI2:
	case ELEM_QUAD1:
	case ELEM_QUAD2:
		RC_TRY( local_inside_coord(face.type, coord, &add_node_number) );
		break;
	case ELEM_TRI1:
		return(NORMAL_RC);
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	for(ii1=0; ii1<add_node_number; ii1++){
		RC_TRY( set_N_vector(face.type, coord[ii1], N) );
		RC_TRY( cal_interpolate_node(*node, N, face.node_num, face.node,
		                             &new_coord) );
		(*new_label) ++;
		RC_TRY( add_node_array(node, new_coord, *new_label) );
		label_box[ii1] = *new_label;
	}
	
	for(ii1=0; ii1<face.from_num; ii1++){
		RC_NEG_CHK( index = search_element_label(elem, face.from[ii1]) );
		switch(face.type){
		case ELEM_TRI2:
			if(ii1 == 0){
				box[index].node[face.from_f[ii1] + 0] = label_box[0];
				box[index].node[face.from_f[ii1] + 1] = label_box[1];
				box[index].node[face.from_f[ii1] + 2] = label_box[2];
			}else{
				if(face.tmp == face.node[0]){
					box[index].node[face.from_f[ii1] + 0] = label_box[0];
					box[index].node[face.from_f[ii1] + 1] = label_box[2];
					box[index].node[face.from_f[ii1] + 2] = label_box[1];
				}else if(face.tmp == face.node[1]){
					box[index].node[face.from_f[ii1] + 0] = label_box[1];
					box[index].node[face.from_f[ii1] + 1] = label_box[0];
					box[index].node[face.from_f[ii1] + 2] = label_box[2];
				}else if(face.tmp == face.node[2]){
					box[index].node[face.from_f[ii1] + 0] = label_box[2];
					box[index].node[face.from_f[ii1] + 1] = label_box[1];
					box[index].node[face.from_f[ii1] + 2] = label_box[0];
				}else{
					return(ARG_ERROR_RC);
				}
			}
			break;
		case ELEM_QUAD1:
			box[index].node[ face.from_f[ii1] ] = label_box[0];
			break;
		case ELEM_QUAD2:
			if(ii1 == 0){
				box[index].node[face.from_f[ii1] + 1] = label_box[1];
				box[index].node[face.from_f[ii1] + 2] = label_box[2];
				box[index].node[face.from_f[ii1] + 3] = label_box[3];
				box[index].node[face.from_f[ii1] + 4] = label_box[4];
			}else{
				if(face.tmp == face.node[0]){
					box[index].node[face.from_f[ii1] + 1] = label_box[3];
					box[index].node[face.from_f[ii1] + 2] = label_box[4];
					box[index].node[face.from_f[ii1] + 3] = label_box[1];
					box[index].node[face.from_f[ii1] + 4] = label_box[2];
				}else if(face.tmp == face.node[1]){
					box[index].node[face.from_f[ii1] + 1] = label_box[2];
					box[index].node[face.from_f[ii1] + 2] = label_box[1];
					box[index].node[face.from_f[ii1] + 3] = label_box[3];
					box[index].node[face.from_f[ii1] + 4] = label_box[4];
				}else if(face.tmp == face.node[2]){
					box[index].node[face.from_f[ii1] + 1] = label_box[4];
					box[index].node[face.from_f[ii1] + 2] = label_box[3];
					box[index].node[face.from_f[ii1] + 3] = label_box[2];
					box[index].node[face.from_f[ii1] + 4] = label_box[1];
				}else if(face.tmp == face.node[3]){
					box[index].node[face.from_f[ii1] + 1] = label_box[1];
					box[index].node[face.from_f[ii1] + 2] = label_box[2];
					box[index].node[face.from_f[ii1] + 3] = label_box[4];
					box[index].node[face.from_f[ii1] + 4] = label_box[3];
				}else{
					return(ARG_ERROR_RC);
				}
			}
			box[index].node[ face.from_f[ii1] ] = label_box[0];
			break;
		default:
			return(ARG_ERROR_RC);
		}
	}

	RC_TRY( free2D(5, 3, &coord) );

	return(NORMAL_RC);
}


static RC
interpolate_coord_elem (NODE_ARRAY *node, ELEMENT_ARRAY elem,
                        ELEMENT element, ELEMENT_NODE_BOX *box,
                        int *new_label )
{
	double N[20];
	double **coord;
	VECT3D new_coord;
	int add_node_number;
	int number;
	int index;
	int ii1;

	RC_TRY( allocate2D(7, 4, &coord) );

	RC_TRY( local_inside_coord(element.type, coord, &add_node_number) );
	switch(element.type){
	case ELEM_LINE1:
		number = 2;
		break;
	case ELEM_LINE2:
		number = 3;
		break;
	case ELEM_TRI2:
		number = 12;
		break;
	case ELEM_QUAD1:
		number = 8;
		break;
	case ELEM_QUAD2:
		number = 16;
		break;
	case ELEM_TETRA2:
		number = 34;
		break;
	case ELEM_PENTA2:
		number = 54;
		break;
	case ELEM_HEXA1:
		number = 26;
		break;
	case ELEM_HEXA2:
		number = 74;
		break;
	case ELEM_TRI1:
	case ELEM_TETRA1:
	case ELEM_PENTA1:
		return(NORMAL_RC);
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	for(ii1=0; ii1<add_node_number; ii1++){
		RC_TRY( set_N_vector(element.type, coord[ii1], N) );

		RC_TRY( cal_interpolate_node(*node, N, element.node_num, element.node,
	                                                             &new_coord) );
		(*new_label) ++;
		RC_TRY( add_node_array(node, new_coord, *new_label) );
	
		RC_NEG_CHK( index = search_element_label(elem, element.label) );
		box[index].node[number + ii1]= *new_label;
	}

	RC_TRY( free2D(7, 4, &coord) );

	return(NORMAL_RC);
}


static RC
local_inside_coord (ELEM_TYPE type, double **coord, int *number)
{
	switch(type){
	case ELEM_LINE1:
		*number = 1;
		coord[0][0] = 0.0;
		break;
	case ELEM_LINE2:
		*number = 2;
		coord[0][0] = -0.5;
		coord[1][0] = 0.5;
		break;
	case ELEM_TRI2:
		*number = 3;
		coord[0][0] = 0.5;
		coord[0][1] = 0.25;
		coord[0][2] = 0.25;
		coord[1][0] = 0.25;
		coord[1][1] = 0.5;
		coord[1][2] = 0.25;
		coord[2][0] = 0.25;
		coord[2][1] = 0.25;
		coord[2][2] = 0.5;
		break;
	case ELEM_QUAD1:
		*number = 1;
		coord[0][0] = 0.0;
		coord[0][1] = 0.0;
		break;
	case ELEM_QUAD2:
		*number = 5;
		coord[0][0] = 0.0;
		coord[0][1] = 0.0;
		coord[1][0] = -0.5;
		coord[1][1] = 0.0;
		coord[2][0] = 0.5;
		coord[2][1] = 0.0;
		coord[3][0] = 0.0;
		coord[3][1] = -0.5;
		coord[4][0] = 0.0;
		coord[4][1] = 0.5;
		break;
	case ELEM_TETRA2:
		*number = 1;
		coord[0][0] = 0.25;
		coord[0][1] = 0.25;
		coord[0][2] = 0.25;
		coord[0][3] = 0.25;
		break;
	case ELEM_PENTA2:
		*number = 3;
		coord[0][0] = 0.5;
		coord[0][1] = 0.25;
		coord[0][2] = 0.25;
		coord[0][3] = 0.0;
		coord[1][0] = 0.25;
		coord[1][1] = 0.5;
		coord[1][2] = 0.25;
		coord[1][3] = 0.0;
		coord[2][0] = 0.25;
		coord[2][1] = 0.25;
		coord[2][2] = 0.5;
		coord[2][3] = 0.0;
		break;
	case ELEM_HEXA1:
		*number = 1;
		coord[0][0] = 0.0;
		coord[0][1] = 0.0;
		coord[0][2] = 0.0;
		break;
	case ELEM_HEXA2:
		*number = 7;
		coord[0][0] = 0.0;
		coord[0][1] = 0.0;
		coord[0][2] = 0.0;
		coord[1][0] = -0.5;
		coord[1][1] = 0.0;
		coord[1][2] = 0.0;
		coord[2][0] = 0.5;
		coord[2][1] = 0.0;
		coord[2][2] = 0.0;
		coord[3][0] = 0.0;
		coord[3][1] = -0.5;
		coord[3][2] = 0.0;
		coord[4][0] = 0.0;
		coord[4][1] = 0.5;
		coord[4][2] = 0.0;
		coord[5][0] = 0.0;
		coord[5][1] = 0.0;
		coord[5][2] = -0.5;
		coord[6][0] = 0.0;
		coord[6][1] = 0.0;
		coord[6][2] = 0.5;
		break;
	case ELEM_TRI1:
	case ELEM_TETRA1:
	case ELEM_PENTA1:
		return(NORMAL_RC);
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
cal_interpolate_node (NODE_ARRAY node, const double *N, 
                      int node_num, const int *node_label, VECT3D *new_coord)
{
	int ii1;
	int index;

	new_coord->x = 0.0;
	new_coord->y = 0.0;
	new_coord->z = 0.0;
	for(ii1=0; ii1<node_num; ii1++){
		RC_NEG_CHK( index = search_node_label(node, node_label[ii1]) ); 
		new_coord->x += N[ii1] * node.array[index].p.x;
		new_coord->y += N[ii1] * node.array[index].p.y;
		new_coord->z += N[ii1] * node.array[index].p.z;
	}

	return(NORMAL_RC);
}


static RC
add_node_array (NODE_ARRAY *node, VECT3D coord, int new_label)
{
	node->size ++;
	RC_TRY( realloc_node_array(node) );

	node->array[node->size - 1].label = new_label;
	node->array[node->size - 1].p.x = coord.x;
	node->array[node->size - 1].p.y = coord.y;
	node->array[node->size - 1].p.z = coord.z;

	return(NORMAL_RC);

}


static int
check_rotation_edge (ELEM_TYPE type, int *index)
{
	switch(type){
	case ELEM_TRI2:
		switch(index[0]){
		case 0:
			if(index[1] == 1) return(0);
			else if(index[1] == 2) return(1);
			else return(-1);
		case 1:
			if(index[1] == 2) return(0);
			else if(index[1] == 0) return(1);
			else return(-1);
		case 2:
			if(index[1] == 0) return(0);
			else if(index[1] == 1) return(1);
			else return(-1);
		default:
			return(-1);
		}
	case ELEM_QUAD2:
		switch(index[0]){
		case 0:
			if(index[1] == 1) return(0);
			else if(index[1] == 3) return(1);
			else return(-1);
		case 1:
			if(index[1] == 2) return(0);
			else if(index[1] == 0) return(1);
			else return(-1);
		case 2:
			if(index[1] == 3) return(0);
			else if(index[1] == 1) return(1);
			else return(-1);
		case 3:
			if(index[1] == 0) return(0);
			else if(index[1] == 2) return(1);
			else return(-1);
		default:
			return(-1);
		}
	case ELEM_TETRA2:
		switch(index[0]){
		case 0:
			if( (index[1] == 1) || (index[1] == 3) ) return(0);
			else if(index[1] == 2) return(1);
			else return(-1);
		case 1:
			if( (index[1] == 2) || (index[1] == 3) ) return(0);
			else if(index[1] == 0) return(1);
			else return(-1);
		case 2:
			if( (index[1] == 0) || (index[1] == 3) ) return(0);
			else if(index[1] == 1) return(1);
			else return(-1);
		case 3:
			return(1);
		default:
			return(-1);
		}
	case ELEM_PENTA2:
		switch(index[0]){
		case 0:
			if( (index[1] == 1) || (index[1] == 3) ) return(0);
			else if(index[1] == 2) return(1);
			else return(-1);
		case 1:
			if( (index[1] == 2) || (index[1] == 4) ) return(0);
			else if(index[1] == 0) return(1);
			else return(-1);
		case 2:
			if( (index[1] == 0) || (index[1] == 5) ) return(0);
			else if(index[1] == 1) return(1);
			else return(-1);
		case 3:
			if(index[1] == 4) return(0);
			else if( (index[1] == 0) || (index[1] == 5) ) return(1);
			else return(-1);
		case 4:
			if(index[1] == 5) return(0);
			else if( (index[1] == 1) || (index[1] == 3) ) return(1);
			else return(-1);
		case 5:
			if(index[1] == 3) return(0);
			else if( (index[1] == 2) || (index[1] == 4) ) return(1);
			else return(-1);
		default:
			return(-1);
		}
	case ELEM_HEXA2:
		switch(index[0]){
		case 0:
			if( (index[1] == 1) || (index[1] == 4) ) return(0);
			else if(index[1] == 3) return(1);
			else return(-1);
		case 1:
			if( (index[1] == 2) || (index[1] == 5) ) return(0);
			else if(index[1] == 0) return(1);
			else return(-1);
		case 2:
			if( (index[1] == 3) || (index[1] == 6) ) return(0);
			else if(index[1] == 1) return(1);
			else return(-1);
		case 3:
			if( (index[1] == 0) || (index[1] == 7) ) return(0);
			else if(index[1] == 2) return(1);
			else return(-1);
		case 4:
			if(index[1] == 5) return(0);
			else if( (index[1] == 0) || (index[1] == 7) ) return(1);
			else return(-1);
		case 5:
			if(index[1] == 6) return(0);
			else if( (index[1] == 1) || (index[1] == 4) ) return(1);
			else return(-1);
		case 6:
			if(index[1] == 7) return(0);
			else if( (index[1] == 2) || (index[1] == 5) ) return(1);
			else return(-1);
		case 7:
			if(index[1] == 4) return(0);
			else if( (index[1] == 3) || (index[1] == 6) ) return(1);
			else return(ARG_ERROR_RC);
		default:
			return(-1);
		}	
	default:
		return(-1);
	}
}


static RC
set_edge (int hash_size, ELEMENT_EDGE *hash, int key_size,
          const int *keys, int label, ELEM_TYPE type, int number)
{
	int index;

	RC_NEG_CHK( index = elem_edge_index(hash_size, hash, key_size, keys) );

	hash[index].from_num ++;
	if( (hash[index].from_num >= EDGE_MAX_FROM)
	 || (hash[index].from_num < 0) ) return(OVERFLOW_ERROR_RC);
	hash[index].from[ hash[index].from_num - 1 ] = label;
	hash[index].from_e[ hash[index].from_num - 1 ] = number;
	hash[index].type = type;

	return(NORMAL_RC);
}


/* return index of hash */
static int
elem_edge_index (int hash_size, ELEMENT_EDGE *hash,
                 int key_size, const int *keys)
{
	int skeys1[3];
	int skeys2[3];
	int ii1, ii2;
	int index;
	int inc;

	sort_keys(key_size, keys, skeys1);
	index = first_index(hash_size, key_size, skeys1);
	inc = 1 << (skeys1[0]%6);   /* 1, 2, 4, 8, 16 or 32 */

	for(ii1=0; ii1<hash_size; ii1++, index+=inc){
		index %= hash_size;
		if(hash[index].from_num < 0)break;
		if(hash[index].node_num != key_size)continue;
		sort_keys(key_size, hash[index].node, skeys2);
		if(equal_keys(key_size, skeys1, skeys2))break;
	}
	if(ii1 >= hash_size) return(-1); /* hash over flow */

	/* set data */
	if(hash[index].from_num < 0){
		hash[index].from_num = 0;
		hash[index].node_num = key_size;
		for(ii2=0;ii2<key_size;ii2++){
			hash[index].node[ii2] = keys[ii2];
		}
	}

	return(index);
}


static RC
set_edge_tri1 (int hash_size, ELEMENT_EDGE *hash, ELEMENT elem)
{
	int keys[2];

	keys[0] = elem.node[0];
	keys[1] = elem.node[1];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 3) );

	keys[0] = elem.node[1];
	keys[1] = elem.node[2];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 4) );

	keys[0] = elem.node[2];
	keys[1] = elem.node[0];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 5) );

	return(NORMAL_RC);
}


static RC
set_edge_tri2 (int hash_size, ELEMENT_EDGE *hash, ELEMENT elem)
{
	int keys[3];

	keys[0] = elem.node[0];
	keys[1] = elem.node[1];
	keys[2] = elem.node[5];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 6) );

	keys[0] = elem.node[1];
	keys[1] = elem.node[2];
	keys[2] = elem.node[3];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 8) );

	keys[0] = elem.node[2];
	keys[1] = elem.node[0];
	keys[2] = elem.node[4];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 10) );

	return(NORMAL_RC);
}


static RC
set_edge_quad1 (int hash_size, ELEMENT_EDGE *hash, ELEMENT elem)
{
	int keys[2];

	keys[0] = elem.node[0];
	keys[1] = elem.node[1];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 4) );

	keys[0] = elem.node[1];
	keys[1] = elem.node[2];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 5) );

	keys[0] = elem.node[2];
	keys[1] = elem.node[3];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 6) );

	keys[0] = elem.node[3];
	keys[1] = elem.node[0];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 7) );

	return(NORMAL_RC);
}


static RC
set_edge_quad2 (int hash_size, ELEMENT_EDGE *hash, ELEMENT elem)
{
	int keys[3];

	keys[0] = elem.node[0];
	keys[1] = elem.node[1];
	keys[2] = elem.node[4];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 9) );

	keys[0] = elem.node[1];
	keys[1] = elem.node[2];
	keys[2] = elem.node[5];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 11) );

	keys[0] = elem.node[2];
	keys[1] = elem.node[3];
	keys[2] = elem.node[6];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 13) );

	keys[0] = elem.node[3];
	keys[1] = elem.node[0];
	keys[2] = elem.node[7];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 15) );

	return(NORMAL_RC);
}


static RC
set_edge_tetra1 (int hash_size, ELEMENT_EDGE *hash, ELEMENT elem)
{
	int keys[2];

	keys[0] = elem.node[1];
	keys[1] = elem.node[2];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 4) );

	keys[0] = elem.node[2];
	keys[1] = elem.node[0];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 5) );

	keys[0] = elem.node[0];
	keys[1] = elem.node[1];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 6) );

	keys[0] = elem.node[0];
	keys[1] = elem.node[3];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 7) );
	
	keys[0] = elem.node[1];
	keys[1] = elem.node[3];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 8) );
	
	keys[0] = elem.node[2];
	keys[1] = elem.node[3];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 9) );
	
	return(NORMAL_RC);
}


static RC
set_edge_tetra2 (int hash_size, ELEMENT_EDGE *hash, ELEMENT elem)
{
	int keys[3];

	keys[0] = elem.node[0];
	keys[1] = elem.node[1];
	keys[2] = elem.node[6];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 10) );

	keys[0] = elem.node[1];
	keys[1] = elem.node[2];
	keys[2] = elem.node[4];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 12) );

	keys[0] = elem.node[2];
	keys[1] = elem.node[0];
	keys[2] = elem.node[5];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 14) );
	
	keys[0] = elem.node[0];
	keys[1] = elem.node[3];
	keys[2] = elem.node[7];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 16) );

	keys[0] = elem.node[1];
	keys[1] = elem.node[3];
	keys[2] = elem.node[8];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 18) );

	keys[0] = elem.node[2];
	keys[1] = elem.node[3];
	keys[2] = elem.node[9];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 20) );

	return(NORMAL_RC);
}


static RC
set_edge_penta1 (int hash_size, ELEMENT_EDGE *hash, ELEMENT elem)
{
	int keys[2];

	keys[0] = elem.node[1];
	keys[1] = elem.node[2];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 6) );

	keys[0] = elem.node[2];
	keys[1] = elem.node[0];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 7) );

	keys[0] = elem.node[0];
	keys[1] = elem.node[1];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 8) );

	keys[0] = elem.node[4];
	keys[1] = elem.node[5];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 9) );

	keys[0] = elem.node[5];
	keys[1] = elem.node[3];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 10) );

	keys[0] = elem.node[3];
	keys[1] = elem.node[4];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 11) );

	keys[0] = elem.node[0];
	keys[1] = elem.node[3];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 12) );

	keys[0] = elem.node[1];
	keys[1] = elem.node[4];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 13) );

	keys[0] = elem.node[2];
	keys[1] = elem.node[5];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 14) );

	return(NORMAL_RC);
}


static RC
set_edge_penta2 (int hash_size, ELEMENT_EDGE *hash, ELEMENT elem)
{
	int keys[3];

	keys[0] = elem.node[0];
	keys[1] = elem.node[1];
	keys[2] = elem.node[8];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 15) );

	keys[0] = elem.node[1];
	keys[1] = elem.node[2];
	keys[2] = elem.node[6];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 17) );

	keys[0] = elem.node[2];
	keys[1] = elem.node[0];
	keys[2] = elem.node[7];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 19) );

	keys[0] = elem.node[3];
	keys[1] = elem.node[4];
	keys[2] = elem.node[11];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 21) );

	keys[0] = elem.node[4];
	keys[1] = elem.node[5];
	keys[2] = elem.node[9];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 23) );

	keys[0] = elem.node[5];
	keys[1] = elem.node[3];
	keys[2] = elem.node[10];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 25) );

	keys[0] = elem.node[0];
	keys[1] = elem.node[3];
	keys[2] = elem.node[12];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 27) );

	keys[0] = elem.node[1];
	keys[1] = elem.node[4];
	keys[2] = elem.node[13];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 29) );

	keys[0] = elem.node[2];
	keys[1] = elem.node[5];
	keys[2] = elem.node[14];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 31) );

	return(NORMAL_RC);
}


static RC
set_edge_hexa1 (int hash_size, ELEMENT_EDGE *hash, ELEMENT elem)
{
	int keys[2];

	keys[0] = elem.node[0];
	keys[1] = elem.node[1];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 8) );

	keys[0] = elem.node[1];
	keys[1] = elem.node[2];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 9) );

	keys[0] = elem.node[2];
	keys[1] = elem.node[3];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 10) );

	keys[0] = elem.node[3];
	keys[1] = elem.node[0];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 11) );

	keys[0] = elem.node[4];
	keys[1] = elem.node[5];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 12) );

	keys[0] = elem.node[5];
	keys[1] = elem.node[6];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 13) );

	keys[0] = elem.node[6];
	keys[1] = elem.node[7];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 14) );

	keys[0] = elem.node[7];
	keys[1] = elem.node[4];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 15) );

	keys[0] = elem.node[0];
	keys[1] = elem.node[4];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 16) );

	keys[0] = elem.node[1];
	keys[1] = elem.node[5];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 17) );

	keys[0] = elem.node[2];
	keys[1] = elem.node[6];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 18) );

	keys[0] = elem.node[3];
	keys[1] = elem.node[7];
	RC_TRY( set_edge(hash_size, hash, 2, keys, elem.label, ELEM_LINE1, 19) );

	return(NORMAL_RC);
}


static RC
set_edge_hexa2 (int hash_size, ELEMENT_EDGE *hash, ELEMENT elem)
{
	int keys[3];

	keys[0] = elem.node[0];
	keys[1] = elem.node[1];
	keys[2] = elem.node[8];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 20) );

	keys[0] = elem.node[1];
	keys[1] = elem.node[2];
	keys[2] = elem.node[9];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 22) );

	keys[0] = elem.node[2];
	keys[1] = elem.node[3];
	keys[2] = elem.node[10];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 24) );

	keys[0] = elem.node[3];
	keys[1] = elem.node[0];
	keys[2] = elem.node[11];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 26) );

	keys[0] = elem.node[4];
	keys[1] = elem.node[5];
	keys[2] = elem.node[12];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 28) );

	keys[0] = elem.node[5];
	keys[1] = elem.node[6];
	keys[2] = elem.node[13];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 30) );

	keys[0] = elem.node[6];
	keys[1] = elem.node[7];
	keys[2] = elem.node[14];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 32) );

	keys[0] = elem.node[7];
	keys[1] = elem.node[4];
	keys[2] = elem.node[15];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 34) );

	keys[0] = elem.node[0];
	keys[1] = elem.node[4];
	keys[2] = elem.node[16];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 36) );

	keys[0] = elem.node[1];
	keys[1] = elem.node[5];
	keys[2] = elem.node[17];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 38) );

	keys[0] = elem.node[2];
	keys[1] = elem.node[6];
	keys[2] = elem.node[18];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 40) );

	keys[0] = elem.node[3];
	keys[1] = elem.node[7];
	keys[2] = elem.node[19];
	RC_TRY( set_edge(hash_size, hash, 3, keys, elem.label, ELEM_LINE2, 42) );

	return(NORMAL_RC);
}

