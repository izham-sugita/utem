/*********************************************************************
 * stl_utl.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Thanks to following contributors.
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: stl_utl.c 842 2006-12-04 06:21:20Z sasaoka $ */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "rc.h"
#include "memory_manager.h"
#include "fem_struct.h"
#include "math_utl.h"
#include "stl_utl.h"
#include "string_utl.h"

#define BDIST 1E-9 /* not define yet */
#define BUFFER_SIZE  (512)
#define ALLOC_SIZE   (1024)
#define INIT_SIZE    (1024)


/*static char *stl_solid  = "solid";*/
static char *stl_eof  = "endsolid";

/* local function */
static RC cal_unit_normal_vect(VECT3D xyz[], VECT3D *un_vect);
static void write_stl_3dface(FILE *fp, const VECT3D xyz[], VECT3D un_vect);
static RC read_stl_3dface(FILE *fp, double *xyz, double *unvect);
static int compar_node(const void *p1, const void *p2);


/* output to STL file */
RC
output_stl_3dface (FILE *fp, char *layer,
                   ELEMENT_ARRAY elem, NODE_ARRAY node)
{
	int ii1,ii2;
	int index;
	VECT3D xyz[3];
	VECT3D un_vect;

	/* print headder of STL */
	fprintf(fp, "solid %s\n", layer);
	
	for(ii1=0;ii1<elem.size;ii1++){
		if(elem.array[ii1].label < 0)continue;
		switch(elem.array[ii1].type){
		case ELEM_TRI1:
		case ELEM_TRI2:
			for(ii2=0;ii2<3;ii2++){
				index = search_node_label(node, elem.array[ii1].node[ii2]);
				if(index < 0) return(SEEK_ERROR_RC);
				xyz[ii2] = node.array[index].p;
			}

			RC_TRY( cal_unit_normal_vect(xyz, &un_vect) );
			write_stl_3dface(fp, xyz, un_vect);
			break;
		case ELEM_QUAD1:
		case ELEM_QUAD2:
			/* 1つめの3角形 */
			index = search_node_label(node, elem.array[ii1].node[0]);
			if(index < 0) return(SEEK_ERROR_RC);
			xyz[0] = node.array[index].p;
			index = search_node_label(node, elem.array[ii1].node[1]);
			if(index < 0) return(SEEK_ERROR_RC);
			xyz[1] = node.array[index].p;
			index = search_node_label(node, elem.array[ii1].node[2]);
			if(index < 0) return(SEEK_ERROR_RC);
			xyz[2] = node.array[index].p;

			RC_TRY( cal_unit_normal_vect(xyz, &un_vect) );
			write_stl_3dface(fp, xyz, un_vect);

			/* 2つめの3角形 */
			index = search_node_label(node, elem.array[ii1].node[2]);
			if(index < 0) return(SEEK_ERROR_RC);
			xyz[0] = node.array[index].p;
			index = search_node_label(node, elem.array[ii1].node[3]);
			if(index < 0) return(SEEK_ERROR_RC);
			xyz[1] = node.array[index].p;
			index = search_node_label(node, elem.array[ii1].node[0]);
			if(index < 0) return(SEEK_ERROR_RC);
			xyz[2] = node.array[index].p;

			RC_TRY( cal_unit_normal_vect(xyz, &un_vect) );
			write_stl_3dface(fp, xyz, un_vect);

			break;
		default:
			return(IMPLEMENT_ERROR_RC);
		}
	}

	/* print footer of STL */
	fprintf(fp, "endsolid %s\n", layer);

	return(NORMAL_RC);
}


static void
write_stl_3dface (FILE *fp, const VECT3D xyz[], VECT3D un_vect)
{
	int ii1;

	fprintf(fp, "  facet normal %e %e %e\n    outer loop\n",
          	un_vect.x, un_vect.y, un_vect.z);
	
	for(ii1=0;ii1<3;ii1++){
		fprintf(fp, "      vertex %e %e %e\n",
		            xyz[ii1].x, xyz[ii1].y, xyz[ii1].z);
	}

	fprintf(fp, "    endloop\n  endfacet\n");
}


static RC
cal_unit_normal_vect (VECT3D xyz[], VECT3D *un_vect)
{

	VECT3D v1, v2, n_vect;
	double abs_vect;

	v1.x = xyz[1].x - xyz[0].x;
	v1.y = xyz[1].y - xyz[0].y;
	v1.z = xyz[1].z - xyz[0].z;
	v2.x = xyz[2].x - xyz[0].x;
	v2.y = xyz[2].y - xyz[0].y;
	v2.z = xyz[2].z - xyz[0].z;

	n_vect = outer_product3d(v1, v2);
	abs_vect = abs_vect3d(n_vect);
	if(nearly_eq(abs_vect, 0.0)){
		un_vect->x = 1.0;
		un_vect->y = 0.0;
		un_vect->z = 0.0;
	}else{
		un_vect->x = n_vect.x / abs_vect;
		un_vect->y = n_vect.y / abs_vect;
		un_vect->z = n_vect.z / abs_vect;
	}

	return(NORMAL_RC);
}


RC
input_stl_3dface (FILE *fp, ELEMENT_ARRAY *element, NODE_ARRAY *node)
{
	char *title = NULL ;
	return( input_stl_3dface_lay(fp, &title, element, node) );
}


RC
input_stl_3dface_lay (FILE *fp, char **title, 
                      ELEMENT_ARRAY *element, NODE_ARRAY *node)
{
	int ii1, ii2;
	RC rc;
	double xyz[12];
	double unvect[3];   /* 法線ベクトル用(未使用) */
	int counter;
	IDC idc[2];

	RC_TRY( allocate_element_array(INIT_SIZE, element) );
	RC_TRY( allocate_node_array(INIT_SIZE, node) );
	element->source = FEM_STL;
	node->source = FEM_STL;

	/* read title */
	rc = fgetidc(fp, "ss", idc);
	if(rc == END_RC) return(END_RC);
	if(rc != NORMAL_RC) return(rc);
	if(strncmp(idc[0].s, stl_eof, strlen(stl_eof)) == 0) return(END_RC);

	(*title) = (char *)mm_alloc(sizeof(strlen(idc[1].s)+1));
	if((*title) == NULL) return(ALLOC_ERROR_RC);
	strcpy((*title), idc[1].s);

	while(1){
		/* read 1 face */
		rc = read_stl_3dface(fp, xyz, unvect);
		if(rc == END_RC) break;
		if(rc != NORMAL_RC) return(rc);

		RC_TRY( realloc_element_array(element) );

		/* put element to array */
		element->array[element->size].label = (element->size) + 1;
		element->array[element->size].type = ELEM_TRI1;
		element->array[element->size].node_num = 3;

		/* put node to array */
		for(ii1=0;ii1<3;ii1++){
			RC_TRY( realloc_node_array(node) );
			node->array[node->size].label = (node->size) + 1;
			node->array[node->size].p.x= xyz[3*ii1];
			node->array[node->size].p.y= xyz[3*ii1+1];
			node->array[node->size].p.z= xyz[3*ii1+2];
			node->array[node->size].i_info[0] = element->size;  /* 逆引き用 */
			node->array[node->size].i_info[1] = ii1;            /* 逆引き用 */
			element->array[element->size].node[ii1] = (node->size) + 1;
			(node->size)++;
		}
		(element->size)++;
	}

	qsort((void *)(node->array), node->size, sizeof(NODE), compar_node);

	counter = 0;
	for(ii1=0; ii1<node->size; ii1++){
		int n_index = ii1;

		for(ii2=ii1-1; ii2>=0; ii2--){
			double dx, dy, dz;

			if(node->array[ii2].label < 0) continue;
			dx = node->array[ii2].p.x - node->array[ii1].p.x;
			if(fabs(dx) > BDIST){
				break;
			}
			dy = node->array[ii2].p.y - node->array[ii1].p.y;
			dz = node->array[ii2].p.z - node->array[ii1].p.z;
			if(sqrt(dx*dx + dy*dy + dz*dz) <= BDIST){
				n_index = ii2;
				break;
			}
		}
		if(n_index == ii1){
			counter++;
			node->array[ii1].label = counter;
			element->array[ node->array[ii1].i_info[0] ]
			         .node[ node->array[ii1].i_info[1] ] = counter;
			node->array[ii1].i_info[0] = 0;
			node->array[ii1].i_info[1] = 0;
		}else{
			element->array[ node->array[ii1].i_info[0] ]
			         .node[ node->array[ii1].i_info[1] ]
			= node->array[n_index].label;
			node->array[ii1].label = -1;
		}
	}
	RC_TRY( clean_element_array(element) );
	RC_TRY( clean_node_array(node) );

	return(NORMAL_RC);
}


static int
compar_node (const void *p1, const void *p2)
{
	if( ((NODE *)p1)->p.x < ((NODE *)p2)->p.x ) return(-1);
	if( ((NODE *)p1)->p.x > ((NODE *)p2)->p.x ) return(1);

	if( ((NODE *)p1)->p.y < ((NODE *)p2)->p.y ) return(-1);
	if( ((NODE *)p1)->p.y > ((NODE *)p2)->p.y ) return(1);

	if( ((NODE *)p1)->p.z < ((NODE *)p2)->p.z ) return(-1);
	if( ((NODE *)p1)->p.z > ((NODE *)p2)->p.z ) return(1);
	
	return(0);
}


static RC
read_stl_3dface (FILE *fp, double *xyz, double *unvect)
{
	int ii1,ii2;
	RC rc;
	IDC idc[5];
	char buf[BUFFER_SIZE];

	for(ii1=0;ii1<(3*3);ii1++){
		xyz[ii1] = 0.0;
	}
	for(ii1=0;ii1<3;ii1++){
		unvect[ii1] = 0.0;
	}

	/* read "facet normal [double] [double] [double]" */
	if(fgets(buf,sizeof(buf), fp) == NULL) return(READ_ERROR_RC);
	rc = sgetidc(buf, "ssddd", idc);
	if(rc == END_RC) return(END_RC);
	if(rc != NORMAL_RC){
		if(strncmp(idc[0].s, stl_eof, strlen(stl_eof)) == 0)
			return(END_RC);
		else
			return(rc);
	}
	for(ii1=0;ii1<3;ii1++){
		unvect[ii1] = idc[ii1+2].d;
	}

	/* read "outer loop" */
	rc = fgetidc(fp, "s", idc);
	if(rc == END_RC) return(END_RC);
	if(rc != NORMAL_RC) return(rc);

	/* read "vertex [double] [double] [double]" */
	for(ii1=0;ii1<3;ii1++){
		rc = fgetidc(fp, "sddd", idc);
		if(rc == END_RC) return(END_RC);
		if(rc != NORMAL_RC) return(rc);

		for(ii2=0;ii2<3;ii2++){
			xyz[3*ii1+ii2] = idc[ii2+1].d;
		}
	}

	/* read "outer loop" & "endfacet" */
	for(ii1=0;ii1<2;ii1++){
		rc = fgetidc(fp, "s", idc);
		if(rc == END_RC) return(END_RC);
		if(rc != NORMAL_RC) return(rc);
	}

	return(NORMAL_RC);
}

