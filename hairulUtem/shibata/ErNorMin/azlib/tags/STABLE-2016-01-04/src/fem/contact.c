/*********************************************************************
 * contact.c
 *
 * Copyright (C) 2004 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: contact.c 849 2006-12-06 08:34:19Z sasaoka $ */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include "base.h"
#include "mathematics.h"
#include "fem_struct.h"

#define CONTACT_TERM_SIZE (10)

static RC contact_terminus_elem(ELEMENT element, NODE_ARRAY node, VECT3D point,
                                CONTACT_TERM *term);
static RC cal_center_point(ELEMENT element, VECT3D p[], VECT3D *center);
static RC contact_terminus_tri1(VECT3D p[], VECT3D point, CONTACT_TERM *term);
static RC contact_terminus_tri2(VECT3D p[], VECT3D point, CONTACT_TERM *term);
static RC contact_terminus_quad1(ELEMENT element, VECT3D p[],
                                 VECT3D point, CONTACT_TERM *term);
static RC contact_terminus_quad2(ELEMENT element, VECT3D p[],
                                 VECT3D point, CONTACT_TERM *term);
static RC allocate_contact_term_array(int num, CONTACT *contact);
static RC realloc_contact_term_array(int num, CONTACT *contact);
static void sort_contact_term_array(CONTACT *contact);
static int compar_contact_term_node(const void *p1, const void *p2);
static RC uniq_contact_term_array(CONTACT *contact);


/*
 * target_node の element に対する接触状況を計算する
 * (この関数自体は実用性は低いと思われる．サンプルとして残す．)
 * 探査配列を毎回作らないように，上位の関数でmake_geom_search_array()
 * を行い，それを使いまわすこと．
 */
#if 0
RC
make_contact_array (ELEMENT_ARRAY element, NODE_ARRAY node,
                    NODE_ARRAY target_node, double range,
                    GEOM_SEGMENT *segment[], CONTACT_ARRAY *contact)
{

	int ii1;
	int size;

	if(!contact) return(ARG_ERROR_RC);
	if(!segment) return(ARG_ERROR_RC);

	RC_TRY( make_geom_search_array(element, node, &size, segment) );

	if(contact->size == 0){
		int num_target = count_valid_node(target_node);
		RC_TRY( allocate_contact_array(num_target, contact) );
	}

	for(ii1=0; ii1<target_node.size; ii1++){
		int index;
		if(target_node.array[ii1].label < 0) continue;
		index = search_contact_node(*contact, target_node.array[ii1].label);
		if(index < 0){ /* 前回の計算結果がない target */
			RC_TRY( realloc_contact_array(contact) );
			RC_TRY( make_contact_terminus(element, node,
			                              target_node.array[ii1], range,
			                              size, segment,
			                              &(contact->array[contact->size])) );
			contact->size++;
		}else{
			RC_TRY( make_contact_terminus(element, node,
			                              target_node.array[ii1], range,
			                              size, segment,
			                              &(contact->array[index])) );
		}
	}

	RC_TRY( sort_contact_array(contact) );

	return(NORMAL_RC);
}
#endif


/* point+-range範囲内にある全 elementにおいて point と各要素との */
/* 最短距離にある交点を計算し，CONTACT_TERM_SIZE 個検出する */
RC
make_contact_terminus (ELEMENT_ARRAY element, NODE_ARRAY node, NODE point,
                       double range, int seg_size,
                       GEOM_SEGMENT *segment[], CONTACT *contact)
{
	int ii1, ii2;
	int size;
	int *cand;

	if(!contact) return(ARG_ERROR_RC);
	if(!segment) return(ARG_ERROR_RC);

	cand = (int *)mm_alloc(seg_size*sizeof(int));
	if(cand == NULL) return(ALLOC_ERROR_RC);

	/* 探査配列の中から，min < ? < max である elementのindex配列を抜き出す */
	RC_TRY( geom_search_candidate(seg_size, segment, point.p, range,
	                              &size, cand) );
	if(size <= 0){
		init_contact(contact);
		contact->node = point.label;
		RC_TRY( mm_free(cand) );
		return(NORMAL_RC);
	}

	RC_TRY( allocate_contact_term_array(size, contact) );

	for(ii1=0, ii2=0; ii1<size; ii1++){
		RC_TRY( contact_terminus_elem(element.array[cand[ii1]], node, point.p,
		                              &(contact->term[ii2])) );
		contact->size++;
		ii2++;
	}

	RC_TRY( mm_free(cand) );

	/* distance の昇順でsort */
	sort_contact_term_array(contact);

	/* 同一距離にあるもの(同一NODEが選ばれることが多い)を削除 */
	RC_TRY( uniq_contact_term_array(contact) );

	if(contact->size > CONTACT_TERM_SIZE){
		RC_TRY( realloc_contact_term_array(CONTACT_TERM_SIZE, contact) );
		contact->sort_flag = ARRAY_SORTED;
	}

	contact->node = point.label;

	return(NORMAL_RC);
}


#if 0
/*
 * 前回計算した最短交点の接触(予測)状態が CONTACT_PLANE であり，且つ
 * その同じ要素について距離を計算し，前回より短くなっている場合は point
 * が一つの要素に限りなく近い場合，もしくは その要素の垂線(法線)方向に
 * pointが==ほぼ==まっすぐ移動していることを示している．この場合，
 * ELEMENT_ARRAY 全体を再計算しても意味はない．重箱の隅のような部分に向
 * かって進んでいる場合，contact->term[]配列の上位に CONTACT_PLANEが幾つ
 * かあり，その後に CONTACT_POINT などが続く状態が考えられる．この場合，
 * 上位の CONTACT_PLANE 全てが前回の距離よりも短くなっている場合以外は
 * 再計算した方が無難．
 */
static RC
chk_recalculate (ELEMENT_ARRAY element, NODE_ARRAY node,
                 NODE point, CONTACT *contact, int *is_recal)
{
	int ii1;

	if(!contact) return(ARG_ERROR_RC);
	if(!contact->term) return(ARG_ERROR_RC);
	if(contact->node != point.label) return(ARG_ERROR_RC);

	for(ii1=0, *is_recal=1; ii1<contact->size; ii1++){
		if(contact->term[ii1].type == CONTACT_PLANE){
			CONTACT_TERM term;
			int index;
			index = search_element_label(element, contact->term[ii1].element);
			if(index < 0) return(SEEK_ERROR_RC);

			/* 前回最短距離にあった elementとの距離を計算 */
			RC_TRY( contact_terminus_elem(element.array[index],
			                              node, point.p, &term) );
			if(term.type == CONTACT_PLANE){
				if( fabs(term.distance) < fabs(contact->term[ii1].distance) ){
					*is_recal = 0; /* 再計算不要 */
					contact->term[ii1] = term;
				}else{
					*is_recal = 1; /* 再計算必要 */
					break;
				}
			}
		}else{
			break; /* CONTACT_PLANE以外が表れたので検査終了 */
		}
	}

	return(NORMAL_RC);
}
#endif


/* element に対する point の最短距離点を求める */
static RC
contact_terminus_elem (ELEMENT element, NODE_ARRAY node,
                       VECT3D point, CONTACT_TERM *term)
{
	int ii1;
	int index[6];
	VECT3D p[8];

	if(element.label < 0) return(ARG_ERROR_RC);
	if(element_dim(element.type) != 2) return(ARG_ERROR_RC);

	for(ii1=0; ii1<element.node_num; ii1++){
		index[ii1] = search_node_label(node, element.node[ii1]);
		if(index[ii1] < 0) return(SEEK_ERROR_RC);
		p[ii1] = node.array[index[ii1]].p;
	}

	if(element.type == ELEM_TRI1){
		RC_TRY( contact_terminus_tri1(p, point, term) );
	}else if(element.type == ELEM_TRI2){
		RC_TRY( contact_terminus_tri2(p, point, term) );
	}else if(element.type == ELEM_QUAD1){
		RC_TRY( contact_terminus_quad1(element, p, point, term) );
	}else if(element.type == ELEM_QUAD2){
		RC_TRY( contact_terminus_quad2(element, p, point, term) );
	}else{
		return(ARG_ERROR_RC);
	}

	term->element = element.label;

	return(NORMAL_RC);
}


/* p[3]で表される面 に対する point の最短点を求める */
static RC
contact_terminus_tri1 (VECT3D p[], VECT3D point, CONTACT_TERM *term)
{
	int ii1;
	int is_inside;
	int is_mins_distance = 0;

	plane_intersect_point(p, point, &term->p, &term->distance, &is_inside);
	if(is_inside){
		term->type = CONTACT_PLANE;
		return(NORMAL_RC);
	}

	/* 面が裏向きである．以下で求まる距離をマイナスにする */
	if(term->distance < 0){
		is_mins_distance = 1;
	}

	term->distance = DBL_MAX;
	term->p = point;
	for(ii1=0; ii1<3; ii1++){
		double tmp_dist;
		VECT3D tmp_term;

		/* 各辺への垂線距離を求める */
		line_intersect_point(p[ii1], p[(ii1==2?0:ii1+1)], point,
		                     &tmp_term, &tmp_dist, &is_inside);
		if(is_inside && (tmp_dist < term->distance)){
			term->distance = tmp_dist;
			term->p = tmp_term;
			term->type = CONTACT_LINE;
		}

		/* 各点への距離を求める */
		tmp_dist = abs_vect3d( sub_vect3d(p[ii1], point) );
		if(tmp_dist < term->distance){
			term->distance = tmp_dist;
			term->p = p[ii1];
			term->type = CONTACT_POINT;
		}
	}

	if(is_mins_distance) term->distance *= -1.0;

	return(NORMAL_RC);
}


static RC
contact_terminus_tri2 (VECT3D p[], VECT3D point, CONTACT_TERM *term)
{
	int ii1;
	VECT3D lp[3];
	CONTACT_TERM lterm[4];
	double dist = DBL_MAX;
	int near_term = 0;

	/* 二次要素を四つのTRI1に分解して計算 */
	/* 0-5-4 */
	lp[0] = p[0];  lp[1] = p[5];  lp[2] = p[4];
	RC_TRY( contact_terminus_tri1(lp, point, &lterm[0]) );
	/* 1-3-5 */
	lp[0] = p[1];  lp[1] = p[3];  lp[2] = p[5];
	RC_TRY( contact_terminus_tri1(lp, point, &lterm[1]) );
	/* 2-4-3 */
	lp[0] = p[2];  lp[1] = p[4];  lp[2] = p[3];
	RC_TRY( contact_terminus_tri1(lp, point, &lterm[2]) );
	/* 3-4-5 */
	lp[0] = p[3];  lp[1] = p[4];  lp[2] = p[5];
	RC_TRY( contact_terminus_tri1(lp, point, &lterm[3]) );
	for(ii1=0; ii1<4;ii1++){
		if( fabs(lterm[ii1].distance) < fabs(dist) ){
			dist = lterm[ii1].distance;
			near_term = ii1;
		}
	}
	*term = lterm[near_term];
	/* 3-4-5 平面を成す直線に近い時，交点は面内にある */
	if(term->type == CONTACT_LINE){
		if(near_term == 3){
			term->type = CONTACT_PLANE;
		}else{
			if( nearly_eq(term->distance, lterm[3].distance) )
			term->type = CONTACT_PLANE;
		}
	}

	return(NORMAL_RC);
}


static RC
contact_terminus_quad1 (ELEMENT element, VECT3D p[],
                        VECT3D point, CONTACT_TERM *term)
{
	int ii1;
	VECT3D lp[3];
	CONTACT_TERM lterm[4];
	double dist = DBL_MAX;
	int near_term = 0;
	VECT3D center;


	RC_TRY( cal_center_point(element, p, &center) );

	/* 0-1-c */
	lp[0] = p[0];  lp[1] = p[1];  lp[2] = center;
	RC_TRY( contact_terminus_tri1(lp, point, &lterm[0]) );
	/* 1-2-c */
	lp[0] = p[1];  lp[1] = p[2];  lp[2] = center;
	RC_TRY( contact_terminus_tri1(lp, point, &lterm[1]) );
	/* 2-3-c */
	lp[0] = p[2];  lp[1] = p[3];  lp[2] = center;
	RC_TRY( contact_terminus_tri1(lp, point, &lterm[2]) );
	/* 3-0-c */
	lp[0] = p[3];  lp[1] = p[0];  lp[2] = center;
	RC_TRY( contact_terminus_tri1(lp, point, &lterm[3]) );

	for(ii1=0; ii1<4;ii1++){
		if( fabs(lterm[ii1].distance) < fabs(dist) ){
			dist = lterm[ii1].distance;
			near_term = ii1;
		}
	}

	*term = lterm[near_term];
	/* 0-1-2-3 平面を成す直線に近い時，交点は面内にある */
	if(term->type == CONTACT_LINE){
		if(near_term == 0){
			if( !is_inside_line(p[0], p[1], term->p) )
				term->type = CONTACT_PLANE;
		}else if(near_term == 1){
			if( !is_inside_line(p[1], p[2], term->p) )
				term->type = CONTACT_PLANE;
		}else if(near_term == 2){
			if( !is_inside_line(p[2], p[3], term->p) )
				term->type = CONTACT_PLANE;
		}else if(near_term == 3){
			if( !is_inside_line(p[3], p[0], term->p) )
				term->type = CONTACT_PLANE;
		}
	}else if(term->type == CONTACT_POINT){
		if( dist_point3d(center, term->p) < ABS_ERROR )
			term->type = CONTACT_PLANE;
	}

	return(NORMAL_RC);
}


static RC
contact_terminus_quad2 (ELEMENT element, VECT3D p[],
                        VECT3D point, CONTACT_TERM *term)
{
	int ii1;
	VECT3D lp[3];
	CONTACT_TERM lterm[8];
	double dist = DBL_MAX;
	int near_term = 0;
	VECT3D center;


	RC_TRY( cal_center_point(element, p, &center) );

	/* 0-4-7 */
	lp[0] = p[0];  lp[1] = p[4];  lp[2] = p[7];
	RC_TRY( contact_terminus_tri1(lp, point, &lterm[0]) );
	/* 1-5-4 */
	lp[0] = p[1];  lp[1] = p[5];  lp[2] = p[4];
	RC_TRY( contact_terminus_tri1(lp, point, &lterm[1]) );
	/* 2-6-5 */
	lp[0] = p[2];  lp[1] = p[6];  lp[2] = p[5];
	RC_TRY( contact_terminus_tri1(lp, point, &lterm[2]) );
	/* 3-7-6 */
	lp[0] = p[3];  lp[1] = p[7];  lp[2] = p[6];
	RC_TRY( contact_terminus_tri1(lp, point, &lterm[3]) );
	/* 4-5-c */
	lp[0] = p[4];  lp[1] = p[5];  lp[2] = center;
	RC_TRY( contact_terminus_tri1(lp, point, &lterm[4]) );
	/* 5-6-c */
	lp[0] = p[5];  lp[1] = p[6];  lp[2] = center;
	RC_TRY( contact_terminus_tri1(lp, point, &lterm[5]) );
	/* 6-7-c */
	lp[0] = p[6];  lp[1] = p[7];  lp[2] = center;
	RC_TRY( contact_terminus_tri1(lp, point, &lterm[6]) );
	/* 7-4-c */
	lp[0] = p[7];  lp[1] = p[4];  lp[2] = center;
	RC_TRY( contact_terminus_tri1(lp, point, &lterm[7]) );

	for(ii1=0; ii1<8;ii1++){
		if( fabs(lterm[ii1].distance) < fabs(dist) ){
			dist = lterm[ii1].distance;
			near_term = ii1;
		}
	}

	*term = lterm[near_term];
	/* 4-5-6-7 を成す直線に近い時，交点は面内にある */
	if(term->type == CONTACT_LINE){
		if(  (near_term == 4) || (near_term == 5) 
		  || (near_term == 6) || (near_term == 7)) { 
			term->type = CONTACT_PLANE;
		}else if(near_term == 0){
			if( nearly_eq(term->distance, lterm[7].distance) )
			term->type = CONTACT_PLANE;
		}else if(near_term == 1){
			if( nearly_eq(term->distance, lterm[4].distance) )
			term->type = CONTACT_PLANE;
		}else if(near_term == 2){
			if( nearly_eq(term->distance, lterm[5].distance) )
			term->type = CONTACT_PLANE;
		}else if(near_term == 3){
			if( nearly_eq(term->distance, lterm[6].distance) )
			term->type = CONTACT_PLANE;
		}
	}else if(term->type == CONTACT_POINT){
		if( dist_point3d(center, term->p) < ABS_ERROR )
			term->type = CONTACT_PLANE;
	}

	return(NORMAL_RC);
}


static RC
cal_center_point (ELEMENT element, VECT3D p[], VECT3D *center)
{
	int ii1;
	double coord[4];
	double N[MAX_NODE];

	if(element.label < -1) return(ARG_ERROR_RC);
	if(element_dim(element.type) != 2) return(ARG_ERROR_RC);

    RC_TRY( set_center_point(element.type, coord) );
    RC_TRY( set_N_vector(element.type, coord, N) );
	RC_TRY( init_vect3d(center) );
    for(ii1=0; ii1<element.node_num; ii1++){
		center->x += N[ii1] * p[ii1].x;
		center->y += N[ii1] * p[ii1].y;
		center->z += N[ii1] * p[ii1].z;
    }

	return(NORMAL_RC);
}


static RC
allocate_contact_term_array (int num, CONTACT *contact)
{
	int ii1;

	if(num > 0){
		contact->term = mm_alloc(num * sizeof(CONTACT_TERM));
		if(contact->term == NULL) return(ALLOC_ERROR_RC);
	}else{
		contact->term = NULL;
	}

	contact->size = 0;
	contact->alloc_size = num;
	contact->sort_flag = ARRAY_UNSORTED;

	for(ii1=0; ii1<num; ii1++){
		init_contact_term( &((contact->term)[ii1]) );
	}

	return(NORMAL_RC);
}


static RC
realloc_contact_term_array (int num, CONTACT *contact)
{
	int ii1;

	if(num > 0){
		contact->term = mm_realloc(contact->term, num*sizeof(CONTACT_TERM));
		if(contact->term == NULL) return(ALLOC_ERROR_RC);
    }else{
        contact->term = NULL;
    }

	if(num > contact->size){
		for(ii1=contact->size-1; ii1<num; ii1++)
			init_contact_term( &((contact->term)[ii1]) );
	}

	contact->size = num;
	contact->alloc_size = num;
	contact->sort_flag = ARRAY_UNSORTED;

    return(NORMAL_RC);
}


static void
sort_contact_term_array (CONTACT *contact)
{
	qsort(contact->term, contact->size, sizeof(CONTACT_TERM),
	      compar_contact_term_node);
	contact->sort_flag = ARRAY_SORTED;
}

static int
compar_contact_term_node (const void *p1, const void *p2)
{
	if( nearly_eq( ((CONTACT_TERM *)p1)->distance,
	               ((CONTACT_TERM *)p2)->distance) ) return(0);

	return(  fabs(((CONTACT_TERM *)p1)->distance)
	       > fabs(((CONTACT_TERM *)p2)->distance) ? 1 : -1 );
}


/* contact->termは既にsortしてあるものとする． */
/* ほぼ同一の距離にあるterm[?]を削除する       */
static RC
uniq_contact_term_array (CONTACT *contact)
{
	int ii1, ii2, ii3;
	int  term_num;
	CONTACT_TERM *term;

	if(!contact) return(ARG_ERROR_RC);
	if(!contact->term) return(ARG_ERROR_RC);

	term = contact->term;

	for(ii1=0,ii2=1,term_num=0; ii2<contact->size; ii1++){
		int is_exist_same_node = 0;
		if( dist_point3d(term[ii1].p, term[ii2].p) < REL_ERROR ){
			ii2++;
			while(ii2 < contact->size){
				if( dist_point3d(term[ii1].p, term[ii2].p) > REL_ERROR ){
					term[ii1+1] = term[ii2];
					term_num++;
					break;
				}
				ii2++;
			}
		}else{
			term_num++;
			if(ii1+1 != ii2)
				term[ii1+1] = term[ii2];
		}
		/* 配列を遡って同じnodeがないか確認 */
		for(ii3=ii1-1; ii3>=0; ii3--){
			if(nearly_eq(term[ii3].distance, term[ii2].distance)){
				if(dist_point3d(term[ii3].p, term[ii2].p) < REL_ERROR){
					is_exist_same_node = 1;
					break;
				}
			}else{
				break;
			}
		}
		if(is_exist_same_node) ii1--;
		ii2++;
	}

	RC_TRY( realloc_contact_term_array(term_num+1, contact) );
	contact->sort_flag = ARRAY_SORTED;

	return(NORMAL_RC);
}

