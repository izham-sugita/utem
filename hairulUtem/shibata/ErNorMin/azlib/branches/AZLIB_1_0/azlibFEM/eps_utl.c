/*********************************************************************
 * eps_utl.c
 *
 * Copyright (C) 2001 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Ryu SASAOKA>
 *
 * Thanks to following contributors.
 *  <Tomoyuki TBUBATA> <Yasuo ASAGA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: eps_utl.c,v 1.8 2003/10/29 08:06:05 sasaoka Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "rc.h"
#include "math_utl.h"
#include "fem_struct.h"
#include "eps_utl.h"


static RC remove_back_face(FEM_ELEMENT_ARRAY *surf, FEM_ELEMENT_ARRAY elem,
                           FEM_NODE_ARRAY node, double **op_matrix);
static RC write_element_eps(FILE *fp, EPS_PROP prop,
                            FEM_ELEMENT elem, FEM_NODE_ARRAY node,
                            double **op_matrix);
static RC sort_element_depth(int *elem_layer, FEM_ELEMENT_ARRAY elem,
                            FEM_NODE_ARRAY node, double **op_matrix);
static RC elem_depth_quick_sort(int first,  int last,
                                double *depth, int *elem_layer);
static RC chk_size_obj(FEM_NODE_ARRAY node, double **op_matrix,
                       double *max_x, double *max_y,
                       double *min_x, double *min_y);


RC output_eps(FILE *fp, EPS_PROP prop, FEM_ELEMENT_ARRAY elem,
              FEM_NODE_ARRAY node, double **op_matrix)
{
	int ii1;
	int *elem_layer;
	double max_x, max_y, min_x, min_y;
	FEM_ELEMENT_ARRAY surf, welem;
	int surf_flag = 0;

	for(ii1=0; ii1<elem.size; ii1++){
		if(elem.array[ii1].label < 0) continue;

		if( (elem.array[ii1].type == ELEM_TETRA1)
		  ||(elem.array[ii1].type == ELEM_TETRA2)
		  ||(elem.array[ii1].type == ELEM_PENTA1)
		  ||(elem.array[ii1].type == ELEM_PENTA2)
		  ||(elem.array[ii1].type == ELEM_HEXA1)
		  ||(elem.array[ii1].type == ELEM_HEXA2) ){
			surf_flag = 1;
		}
	}
	if(surf_flag){
		RC_TRY( extract_surface(elem, &surf) );
		RC_TRY( remove_back_face(&surf, elem, node, op_matrix) );
		welem = surf;
	}else{
		welem = elem;
	}

	RC_NULL_CHK( elem_layer = (int *)malloc(sizeof(int) * welem.size) );
	for(ii1=0;ii1<welem.size;ii1++){
		elem_layer[ii1] = -1;
	}

	RC_TRY( sort_element_depth(elem_layer, welem, node, op_matrix) );

	RC_TRY( chk_size_obj(node, op_matrix, &max_x, &max_y, &min_x, &min_y) );

	fprintf(fp, "%%!PS-Adobe-3.1 EPSF-3.0\n"
	            "%%%%BoundingBox: %d %d %d %d\n"
	            "%3.1f setlinewidth\n"
	            "1 setlinejoin\n",
	            (int)(min_x - prop.line_width),
				(int)(min_y - prop.line_width),
				(int)(max_x + prop.line_width),
				(int)(max_y + prop.line_width), prop.line_width);

	for(ii1=0;ii1<welem.size;ii1++){
		if(elem_layer[ii1] < 0) continue;

		RC_TRY( write_element_eps(fp, prop, welem.array[elem_layer[ii1]],
		                          node, op_matrix) );
	}

	fprintf(fp, "showpage\n");

	if(surf_flag){
		RC_TRY( free_fem_element_array(&surf) );
	}

	return(NORMAL_RC);
}


static RC remove_back_face(FEM_ELEMENT_ARRAY *surf, FEM_ELEMENT_ARRAY elem,
                           FEM_NODE_ARRAY node, double **op_matrix)
{
	int ii1;
	int index;
	TRANSLATION3D vect;

	if(surf->source != FEM_GENERATED_SURFACE) return(ARG_ERROR_RC);

	for(ii1=0; ii1<surf->size; ii1++){
		if(surf->array[ii1].label < 0) continue;

		index = search_fem_element_label(elem, surf->array[ii1].i_info[0]);
		RC_NEG_CHK(index);

		RC_TRY( unit_normal_vect_center(surf->array[ii1], elem.array[index],
		                                node, &vect) );
		vect = operation3d(op_matrix, vect);
		if(vect.z <= 0.0){
			surf->array[ii1].label = -1;
		}
	}
	RC_TRY( clean_fem_element_array(surf) );

	return(NORMAL_RC);
}


static RC write_element_eps(FILE *fp, EPS_PROP prop,
                            FEM_ELEMENT elem, FEM_NODE_ARRAY node,
                            double **op_matrix)
{
	int ii1;
	int index;
	TRANSLATION3D point;

	for(ii1=0;ii1<elem.node_num;ii1++){
		RC_NEG_CHK( index = search_fem_node_label(node, elem.node[ii1]) );
		point = operation3d(op_matrix, node.array[index].p);

		if(ii1 == 0){
			fprintf(fp, "newpath\n"
			            "%7.3f %7.3f moveto\n", point.x, point.y);
		}else{
			fprintf(fp, "%7.3f %7.3f lineto\n", point.x, point.y);
		}
	}

	fprintf(fp, "closepath\n"
	            "gsave\n"
	            "%6.4f %6.4f %6.4f setrgbcolor fill\n"
	            "grestore\n"
	            "%6.4f %6.4f %6.4f setrgbcolor stroke\n",
	            prop.fill_rgb_color[0],
 	            prop.fill_rgb_color[1],
	            prop.fill_rgb_color[2],
	            prop.line_rgb_color[0],
	            prop.line_rgb_color[1],
	            prop.line_rgb_color[2]);

	return(NORMAL_RC);

}


static RC sort_element_depth(int *elem_layer, FEM_ELEMENT_ARRAY elem,
                            FEM_NODE_ARRAY node, double **op_matrix)
{


	int ii1,ii2;
	double *depth;
	TRANSLATION3D point;

	RC_NULL_CHK( depth =(double *)malloc(sizeof(double) * elem.size) );

	for(ii1=0;ii1<elem.size;ii1++){
		depth[ii1] = 0.0;

		if(elem.array[ii1].label < 0){
			elem_layer[ii1] = -1;
			continue;
		}

		elem_layer[ii1] = ii1;
		for(ii2=0; ii2<elem.array[ii1].node_num; ii2++){
			int index = search_fem_node_label(node, elem.array[ii1].node[ii2]);
			RC_NEG_CHK(index);

			point = operation3d(op_matrix, node.array[index].p);
			depth[ii1] += point.z;
		}
		depth[ii1] /= (double)(elem.array[ii1].node_num);
	}

	/* ・ッ・、・テ・ッ・ス。シ・ネ */
	RC_TRY( elem_depth_quick_sort(0, elem.size-1, depth, elem_layer)); 

	free((void *)depth);

	return(NORMAL_RC);

}


static RC elem_depth_quick_sort(int first,  int last,
                                double *depth, int *elem_layer)
{
	int top, tail, i_tmp;
	double key, d_tmp;

	key = depth[(first+last)/2];
	top = first;
	tail = last;

	for(;;top++,tail--){
		while(depth[top] < key) top++;
		while(depth[tail] > key) tail--;
		if(top >= tail) break;

		d_tmp = depth[tail];
		depth[tail] = depth[top];
		depth[top] = d_tmp;
		i_tmp = elem_layer[tail];
		elem_layer[tail] = elem_layer[top];
		elem_layer[top] = i_tmp;
	}

	if(first < top-1)
		elem_depth_quick_sort(first, top-1, depth, elem_layer);
	if(tail+1 < last)
		elem_depth_quick_sort(tail+1, last, depth, elem_layer);

	return(NORMAL_RC);
}


static RC chk_size_obj(FEM_NODE_ARRAY node, double **op_matrix,
                       double *max_x, double *max_y,
                       double *min_x, double *min_y)
{
	int ii1;
	TRANSLATION3D point;

	if(node.size <= 0) return(ARG_ERROR_RC);

	point = operation3d(op_matrix, node.array[0].p);
	*max_x = point.x;
	*max_y = point.y;
	*min_x = point.x;
	*min_y = point.y;

	for(ii1=1; ii1<node.size; ii1++){
		point = operation3d(op_matrix, node.array[ii1].p);
		if(*max_x < point.x) *max_x = point.x;
		if(*max_y < point.y) *max_y = point.y;
		if(*min_x > point.x) *min_x = point.x;
		if(*min_y > point.y) *min_y = point.y;
	}

	return(NORMAL_RC);
}


