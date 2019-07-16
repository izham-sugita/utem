/*********************************************************************
 * gtkgl_fem.c
 *
 * Copyright (C) 2002 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: gtkgl_fem.c 415 2005-07-15 06:04:24Z sasaoka $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <glib.h>
#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>
#include <gtk/gtkgl.h>
#include "gtk_utl.h"
#include "fem_struct.h"
#include "math_utl.h"
#include "rc.h"

/* Global variable */
extern GTKU_GL *gtku_gl;

static RC draw_fem_object_1D(GTKU_OBJECT *obj, GTKU_DRAW_TYPE draw_type);
static RC draw_fem_object_2D(GTKU_OBJECT *obj, GTKU_DRAW_TYPE draw_type);
static RC draw_fem_object_3D(GTKU_OBJECT *obj, GTKU_DRAW_TYPE draw_type);
static RC draw_edge(GTKU_OBJECT *obj, ELEMENT elem,
                    EDGE_TABLE e_table, GLdouble color[][4]);
static RC draw_face(GTKU_OBJECT *obj, ELEMENT elem,
                    FACE_TABLE f_table, GLdouble color[][4],
                    VECT3D nvect[]);
static RC draw_face_single(GTKU_OBJECT *obj, ELEMENT elem,
                           FACE_TABLE f_table, GLdouble color[][4]);
static RC cal_normal_vector(NODE_ARRAY node, ELEMENT elem,
                            VECT3D *nvect);
static RC drow_fem_polygon_2D(GTKU_OBJECT *obj);
static RC drow_fem_polygon_3D(GTKU_OBJECT *obj);
static RC drow_fem_solid_3D(GTKU_OBJECT *obj);
static RC drow_fem_wire_frame(GTKU_OBJECT *obj);
static RC drow_fem_grid_frame(GTKU_OBJECT *obj);
static RC drow_fem_surf_point(GTKU_OBJECT *obj);
static RC drow_fem_grid_point(GTKU_OBJECT *obj);
static RC drow_fem_hiddenline(GTKU_OBJECT *obj);
static RC draw_fem_number(GTKU_OBJECT *obj, GTKU_DRAW_TYPE draw_type);
static RC draw_fem_number_surf_point(GTKU_OBJECT *obj);
static RC draw_fem_number_grid_point(GTKU_OBJECT *obj);
static RC draw_fem_number_solid(GTKU_OBJECT *obj);
static RC draw_fem_number_polygon(GTKU_OBJECT *obj);


RC
gtku_gl_draw_fem_object (GTKU_OBJECT *obj)
{
	GTKU_DRAW_TYPE draw_type;

	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	RC_NULL_CHK(obj);
	draw_type = obj->d_prop->draw_type;
	if(draw_type & GTKU_DRAW_HIDE) return(NORMAL_RC);
	if(draw_type & GTKU_DRAW_DEFAULT)
		draw_type |= gtku_gl->draw_type;

	if(draw_type & GTKU_DRAW_BBOX)
		RC_TRY( gtku_gl_draw_bbox(obj->d_prop) );

	/* objectがnodeだけを持つとき */
	if(  (obj->type & GTKU_OBJECT_FEM_NODE) &&
	    !(obj->type & GTKU_OBJECT_FEM_ELEM) ){
		draw_type = GTKU_DRAW_GRID_POINT; /* cancel other draw_type */
		obj->dim = GTKU_OBJECT_3D;
	}

	/* Display node and element number */
	if( (obj->d_prop->draw_node_num_sw) ||
	    (obj->d_prop->draw_elem_num_sw) ||
	    (obj->d_prop->draw_surf_num_sw) )
		RC_TRY( draw_fem_number(obj, draw_type) );

	if(obj->dim == GTKU_OBJECT_1D){
		RC_TRY( draw_fem_object_1D(obj, draw_type) );
	}else if(obj->dim == GTKU_OBJECT_2D){
		RC_TRY( draw_fem_object_2D(obj, draw_type) );
	}else if(obj->dim == GTKU_OBJECT_3D){
		RC_TRY( draw_fem_object_3D(obj, draw_type) );
	}else{
		GTKU_WARNING("draw_fem_object()", "Unknown object dimension");
	}

	return(NORMAL_RC);
}


/* 1D(BAR, LINE, BEAM) */
static RC
draw_fem_object_1D (GTKU_OBJECT *obj, GTKU_DRAW_TYPE draw_type)
{
	if(draw_type & (GTKU_DRAW_GRID_POINT | GTKU_DRAW_SURF_POINT)){
		if(gtku_object_exist_chk(obj, 1, "node")){
			glPushName(GTKU_DRAW_GRID_POINT);
			RC_TRY(drow_fem_grid_point(obj));
			glPopName();
		}else GTKU_WARNING("draw_fem_object_1D()", "none object (Point)");
	}

	if(draw_type & (GTKU_DRAW_POLYGON    | GTKU_DRAW_SOLID      |
	                GTKU_DRAW_GRID_FRAME | GTKU_DRAW_WIRE_FRAME |
	                GTKU_DRAW_HIDDENLINE) ){
		if(gtku_object_exist_chk(obj, 2, "node", "element")){
			glPushName(GTKU_DRAW_WIRE_FRAME);
			RC_TRY(drow_fem_wire_frame(obj))
			glPopName();
		}else
			GTKU_WARNING("draw_fem_object_1D()", "none object (other POINT)");
	}

	return(NORMAL_RC);
}


/* 2D(SHELL, TRI, QUAD) */
static RC
draw_fem_object_2D (GTKU_OBJECT *obj, GTKU_DRAW_TYPE draw_type)
{
	if(draw_type & GTKU_DRAW_GRID_POINT){
		if( gtku_object_exist_chk(obj, 1, "node") ){
			glPushName(GTKU_DRAW_GRID_POINT);
			RC_TRY(drow_fem_grid_point(obj));
			glPopName();
		}else GTKU_WARNING("draw_fem_object_2D()", "none object (Grid Point)");
	}

	if(draw_type & GTKU_DRAW_SURF_POINT){
		if( gtku_object_exist_chk(obj, 2, "node", "surface") ){
			glPushName(GTKU_DRAW_SURF_POINT);
			RC_TRY(drow_fem_surf_point(obj));
			glPopName();
		}else GTKU_WARNING("draw_fem_object_2D()","none object (SurfacePoint)");
	}

	if(draw_type & GTKU_DRAW_GRID_FRAME){
		if( gtku_object_exist_chk(obj, 2, "node", "element") ){
			glPushName(GTKU_DRAW_GRID_FRAME);
			RC_TRY(drow_fem_grid_frame(obj));
			glPopName();
		}else GTKU_WARNING("draw_fem_object_2D()", "none object (Grid Frame)");
	}

	if(draw_type & GTKU_DRAW_WIRE_FRAME){
		if( gtku_object_exist_chk(obj, 2, "node", "surface") ){
			glPushName(GTKU_DRAW_WIRE_FRAME);
			RC_TRY(drow_fem_wire_frame(obj));
			glPopName();
		}else GTKU_WARNING("draw_fem_object_2D()", "none object (Wire Frame)");
	}

	if(draw_type & GTKU_DRAW_HIDDENLINE){
		if( gtku_object_exist_chk(obj, 2, "node", "element") ){
			glPushName(GTKU_DRAW_HIDDENLINE);
			RC_TRY(drow_fem_hiddenline(obj));
			glPopName();
		}else GTKU_WARNING("draw_fem_object_2D()", "none object (HiddenLine)");
	}

	if(draw_type & (GTKU_DRAW_POLYGON|GTKU_DRAW_SOLID)){
		if( gtku_object_exist_chk(obj, 2, "node", "element") ){
			glPushName(GTKU_DRAW_POLYGON);
			RC_TRY(drow_fem_polygon_2D(obj));
			glPopName();
		}else GTKU_WARNING("draw_fem_object_2D()","none object (Polygon)");
	}

	return(NORMAL_RC);
}


/* 3D(TETRA, PENTA, HEXA) */
static RC
draw_fem_object_3D (GTKU_OBJECT *obj, GTKU_DRAW_TYPE draw_type)
{
	if(draw_type & GTKU_DRAW_GRID_POINT){
		if( gtku_object_exist_chk(obj, 1, "node") ){
			glPushName(GTKU_DRAW_GRID_POINT);
			RC_TRY(drow_fem_grid_point(obj));
			glPopName();
		}else GTKU_WARNING("draw_fem_object_3D()", "none object (Grid Point)");
	}

	if(draw_type & GTKU_DRAW_SURF_POINT){
		if( gtku_object_exist_chk(obj, 2, "node", "surface") ){
			glPushName(GTKU_DRAW_SURF_POINT);
			RC_TRY(drow_fem_surf_point(obj));
			glPopName();
		}else GTKU_WARNING("draw_fem_object_3D()","none object (SurfacePoint)");
	}

	if(draw_type & GTKU_DRAW_GRID_FRAME){
		if( gtku_object_exist_chk(obj, 2, "node", "element") ){
			glPushName(GTKU_DRAW_GRID_FRAME);
			RC_TRY(drow_fem_grid_frame(obj));
			glPopName();
		}else GTKU_WARNING("draw_fem_object_3D()", "none object (Grid Frame)");
	}

	if(draw_type & GTKU_DRAW_WIRE_FRAME){
		if( gtku_object_exist_chk(obj, 2, "node", "surface") ){
			glPushName(GTKU_DRAW_WIRE_FRAME);
			RC_TRY(drow_fem_wire_frame(obj));
			glPopName();
		}else GTKU_WARNING("draw_fem_object_3D()", "none object(Wire Frame)");
	}

	if(draw_type & GTKU_DRAW_HIDDENLINE){
		if( gtku_object_exist_chk(obj, 3, "node", "element", "surface") ){
			glPushName(GTKU_DRAW_HIDDENLINE);
			RC_TRY(drow_fem_hiddenline(obj));
			glPopName();
		}else GTKU_WARNING("draw_fem_object_3D()", "none object (HiddenLine)");
	}

	if(draw_type & GTKU_DRAW_SOLID){
		if( gtku_object_exist_chk(obj, 2, "node", "element") ){
			glPushName(GTKU_DRAW_SOLID);
			RC_TRY(drow_fem_solid_3D(obj));
			glPopName();
		}else GTKU_WARNING("draw_fem_object_3D()", "none object (Solid)");
	}

	if(draw_type & GTKU_DRAW_POLYGON){
		if( gtku_object_exist_chk(obj, 3, "node", "element", "surface") ){
			glPushName(GTKU_DRAW_POLYGON);
			RC_TRY(drow_fem_polygon_3D(obj));
			glPopName();
		}else GTKU_WARNING("draw_fem_object_3D()","none object (Polygon)");
	}

	return(NORMAL_RC);
}


static RC
draw_edge (GTKU_OBJECT *obj, ELEMENT elem,
           EDGE_TABLE e_table, GLdouble color[][4])
{
	int ii1, ii2;

	if(element_order(elem.type) == 1){
		for(ii1=0; ii1<e_table.edge_num; ii1++){
			glBegin(GL_LINES);
			for(ii2=0; ii2<e_table.node_num[ii1]; ii2++){
				int label = elem.node[(e_table.table[ii1][ii2])];
				int node_index = search_node_label(*(obj->node), label);
				RC_NEG_CHK(node_index);
				if(color == NULL) glColor4dv(obj->d_prop->line_color);
				else glColor4dv(color[e_table.table[ii1][ii2]]);
				glVertex3d( obj->node->array[node_index].p.x,
				            obj->node->array[node_index].p.y,
				            obj->node->array[node_index].p.z );
			}
			glEnd();
		}
	}else if(element_order(elem.type) == 2){
		for(ii1=0; ii1<e_table.edge_num; ii1++){
			/* Fix!!                                               */
			/* 内挿関数使って 一線を8分割ぐらいして，通過点を求め，*/
			/* それをGL_LINEで結んだ方が色表現的にもいい．         */
			double ctlp[3][3];
			int label0 = elem.node[ e_table.table[ii1][0] ];
			int label1 = elem.node[ e_table.table[ii1][1] ];
			int label2 = elem.node[ e_table.table[ii1][2] ];
			int node_index0 = search_node_label(*(obj->node), label0);
			int node_index1 = search_node_label(*(obj->node), label1);
			int node_index2 = search_node_label(*(obj->node), label2);
			RC_NEG_CHK(node_index0);
			RC_NEG_CHK(node_index1);
			RC_NEG_CHK(node_index2);
			/* e_table->table[0][0] */
			ctlp[0][0] = obj->node->array[node_index0].p.x;
			ctlp[0][1] = obj->node->array[node_index0].p.y;
			ctlp[0][2] = obj->node->array[node_index0].p.z;
			/* e_table->table[0][1] */
			ctlp[2][0] = obj->node->array[node_index1].p.x;
			ctlp[2][1] = obj->node->array[node_index1].p.y;
			ctlp[2][2] = obj->node->array[node_index1].p.z;
			{ /* e_table->table[0][2] 2次要素中点 */
			 VECT3D mid = mid_point3d(obj->node->array[node_index0].p,
			                          obj->node->array[node_index1].p);
			 VECT3D sub = sub_vect3d(obj->node->array[node_index2].p, mid);
			 VECT3D add = add_vect3d(obj->node->array[node_index2].p, sub);
			 ctlp[1][0] = add.x;
			 ctlp[1][1] = add.y;
			 ctlp[1][2] = add.z;
			}
			
			glColor4dv(obj->d_prop->line_color);
			glMap1d(GL_MAP1_VERTEX_3, 0.0, 1.0, 3, 3, &ctlp[0][0]);
			glEnable(GL_MAP1_VERTEX_3);

			glMapGrid1d(10, 0.0, 1.0);
			glEvalMesh1(GL_LINE, 0, 10); /* GL_LINESではない */
			/* 上の2行は以下の方法と同じ
			glBegin(GL_LINE_STRIP);
			{
				int ii3;
				for(ii3=0; ii3<= 10; ii3++)
					glEvalCoord1f((GLfloat)ii3/10);
			}
			glEnd();
			*/
		}
	}else{
		GTKU_WARNING("draw_edge()", "Illegal element order");
	}

	return(NORMAL_RC);
}


static RC
draw_face (GTKU_OBJECT *obj, ELEMENT elem,
           FACE_TABLE f_table, GLdouble color[][4], VECT3D nvect[])
{
	if(element_dim(elem.type) == 2){
		glBegin(GL_POLYGON);
		if(nvect != NULL)
			glNormal3d(nvect->x, nvect->y, nvect->z);
		RC_TRY( draw_face_single(obj, elem, f_table, color) );
		glEnd();
	}else if(element_dim(elem.type) == 3){
		int ii1, ii2;
		ELEMENT local_elem;
		FACE_TABLE local_f_table;
		for(ii1=0; ii1<f_table.face_num; ii1++){
			local_elem.type =  f_table.face_type[ii1];
			RC_TRY( make_face_table(local_elem.type, &local_f_table) );
			for(ii2=0; ii2<f_table.node_num[ii1]; ii2++)
				local_elem.node[ii2] = elem.node[f_table.table[ii1][ii2]];
			glBegin(GL_POLYGON);
			if(nvect == NULL){
				VECT3D l_nvect;
				RC_TRY( cal_normal_vector(*(obj)->node, local_elem, &l_nvect) );
				glNormal3d(l_nvect.x, l_nvect.y, l_nvect.z);
			}else glNormal3d(nvect->x, nvect->y, nvect->z);
			RC_TRY( draw_face_single(obj, local_elem, local_f_table, color) );
			glEnd();
		}
	}else{
		GTKU_WARNING("draw_edge()", "Illegal element dimension");
	}

	return(NORMAL_RC);
}


static RC
draw_face_single (GTKU_OBJECT *obj, ELEMENT elem,
                  FACE_TABLE f_table, GLdouble color[][4])
{
	int ii1, ii2;

	/* Fix!! 2次要素はこれでは表示できない */
	for(ii1=0; ii1<f_table.face_num; ii1++){
		/* for(ii2=0; ii2<f_table.node_num[ii1]; ii2++){ */
		for(ii2=0; ii2<2; ii2++){
			int label = elem.node[ f_table.table[ii1][ii2] ];
			int node_index = search_node_label(*(obj->node), label);
			RC_NEG_CHK(node_index);
			if(color == NULL) glColor4dv(obj->d_prop->polygon_color);
			else glColor4dv(color[f_table.table[ii1][ii2]]);
			glVertex3d( obj->node->array[node_index].p.x,
			            obj->node->array[node_index].p.y,
			            obj->node->array[node_index].p.z );
		}
	}
	return(NORMAL_RC);
}


/* 2次元要素の法線を簡易的に計算 */
static RC
cal_normal_vector (NODE_ARRAY node, ELEMENT elem, VECT3D *nvect)
{
	int node_index0, node_index1, node_index2;
	VECT3D v1, v2;

	if( !((elem.type == ELEM_TRI1)  ||
	      (elem.type == ELEM_TRI2)  ||
	      (elem.type == ELEM_QUAD1) ||
	      (elem.type == ELEM_QUAD2)) ) return(ARG_ERROR_RC);

	RC_NEG_CHK( node_index0 = search_node_label(node, elem.node[0]) );
	RC_NEG_CHK( node_index1 = search_node_label(node, elem.node[1]) );
	RC_NEG_CHK( node_index2 = search_node_label(node, elem.node[2]) );
	v1 = sub_vect3d( node.array[node_index1].p,
	                 node.array[node_index0].p );
	v2 = sub_vect3d( node.array[node_index2].p,
	                 node.array[node_index0].p );
	*nvect = outer_product3d(v1, v2);
	return(NORMAL_RC);
}


static RC
drow_fem_polygon_2D (GTKU_OBJECT *obj)
{
	int ii1;
	FACE_TABLE f_table;

	RC_NULL_CHK(obj);
	RC_NULL_CHK(obj->node);
	RC_NULL_CHK(obj->elem);
	RC_NULL_CHK(obj->d_prop);

	/* translucent start */
	if(obj->d_prop->polygon_trans_sw) set_gl_translucent();

	glEnable(GL_LIGHTING);

	init_face_table(&f_table);
	for(ii1=0; ii1<obj->elem->size; ii1++){
		VECT3D nvect[1];
		ELEMENT elem = obj->elem->array[ii1];
		if(elem.label < 0) continue;
		/* Cal and Set Normal Vector */
		RC_TRY( cal_normal_vector(*(obj)->node, elem, &nvect[0]) );
		RC_TRY( make_face_table(elem.type, &f_table) );
		/* Drawing polygon */
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glPushName(elem.label);
		RC_TRY( draw_face(obj, elem, f_table, NULL, nvect) );
		glPopName();
	}

	/* translucent end */
	if(obj->d_prop->polygon_trans_sw) unset_gl_translucent();

	return(NORMAL_RC);
}


static RC
drow_fem_polygon_3D (GTKU_OBJECT *obj)
{
	int ii1;
	int elem_index;
	FACE_TABLE f_table;
	VECT3D nvect[1];

	RC_NULL_CHK(obj);
	RC_NULL_CHK(obj->node);
	RC_NULL_CHK(obj->elem);
	RC_NULL_CHK(obj->surf);
	RC_NULL_CHK(obj->d_prop);

	/* translucent start */
	if(obj->d_prop->polygon_trans_sw) set_gl_translucent();

	glEnable(GL_LIGHTING);

	init_face_table(&f_table);
	for(ii1=0; ii1<obj->surf->size; ii1++){
		ELEMENT surf = obj->surf->array[ii1];
		if(surf.label < 0) continue;
		/* Cal and Set Normal Vector */
		elem_index = search_element_label(*(obj)->elem, surf.i_info[0]);
		RC_NEG_CHK(elem_index);
		RC_TRY( unit_normal_vect_center(surf, obj->elem->array[elem_index],
		                                *(obj)->node, &nvect[0]) );
		RC_TRY( make_face_table(surf.type, &f_table) );
		/* Drawing polygon */
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glPushName(surf.label);
		RC_TRY( draw_face(obj, surf, f_table, NULL, nvect) );
		glPopName();
	}

	/* translucent end */
	if(obj->d_prop->polygon_trans_sw) unset_gl_translucent();

	return(NORMAL_RC);
}


static RC
drow_fem_solid_3D (GTKU_OBJECT *obj)
{
	int ii1;
	FACE_TABLE f_table;

	RC_NULL_CHK(obj);
	RC_NULL_CHK(obj->node);
	RC_NULL_CHK(obj->elem);
	RC_NULL_CHK(obj->d_prop);

	/* translucent start */
	if(obj->d_prop->polygon_trans_sw) set_gl_translucent();

	glEnable(GL_LIGHTING);

	init_face_table(&f_table);
	for(ii1=0; ii1<obj->elem->size; ii1++){
		ELEMENT elem = obj->elem->array[ii1];;
		if(elem.label < 0) continue;
		RC_TRY( make_face_table(elem.type, &f_table) );
		/* Drawing polygon */
		/* 法線に関してはとりあえず簡易的な方法で対処 */
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glPushName(elem.label);
		RC_TRY( draw_face(obj, elem, f_table, NULL, NULL) );
		glPopName();
	}

	/* translucent end */
	if(obj->d_prop->polygon_trans_sw) unset_gl_translucent();

	return(NORMAL_RC);
}


static RC
drow_fem_wire_frame (GTKU_OBJECT *obj)
{
	int ii1;
	EDGE_TABLE e_table;

	RC_NULL_CHK(obj);
	RC_NULL_CHK(obj->surf);
	RC_NULL_CHK(obj->d_prop);

	/* translucent start */
	if(obj->d_prop->line_trans_sw) set_gl_translucent();

	glDisable(GL_LIGHTING);

	init_edge_table(&e_table);
	for(ii1=0; ii1<obj->surf->size; ii1++){
		if(obj->surf->array[ii1].label < 0) continue;
		glLineWidth((float)obj->d_prop->line_size);
		if(element_dim(obj->surf->array[ii1].type) == 1){
			switch(element_order(obj->surf->array[ii1].type) ){
			case 1:
				e_table.edge_num = 1;
				e_table.node_num[0] = 2;
				e_table.edge_type[0] = ELEM_LINE1;
				e_table.table[0][0] = 0;
				e_table.table[0][1] = 1;
				break;
			case 2:
				e_table.edge_num = 1;
				e_table.node_num[0] = 2;
				e_table.edge_type[0] = ELEM_LINE2;
				e_table.table[0][0] = 0;
				e_table.table[0][1] = 1;
				e_table.table[0][2] = 2;
				break;
			default:
				break;
			}
		}else{
			RC_TRY( make_edge_table(obj->surf->array[ii1].type, &e_table) );
		}
		glPushName(obj->surf->array[ii1].label);
		RC_TRY( draw_edge(obj, obj->surf->array[ii1], e_table, NULL) );
		glPopName();
	}

	/* translucent end */
	if(obj->d_prop->line_trans_sw) unset_gl_translucent();

	return(NORMAL_RC);
}


static RC
drow_fem_grid_frame (GTKU_OBJECT *obj)
{
	int ii1;
	EDGE_TABLE e_table;

	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	RC_NULL_CHK(obj);
	RC_NULL_CHK(obj->elem);
	RC_NULL_CHK(obj->d_prop);

	/* translucent start */
	if(obj->d_prop->line_trans_sw) set_gl_translucent();

	glDisable(GL_LIGHTING);

	init_edge_table(&e_table);
	for(ii1=0; ii1<obj->elem->size; ii1++){
		if(obj->elem->array[ii1].label < 0) continue;
		glLineWidth((float)obj->d_prop->line_size);
		RC_TRY( make_edge_table(obj->elem->array[ii1].type, &e_table) );
		glPushName(obj->elem->array[ii1].label);
		RC_TRY( draw_edge(obj, obj->elem->array[ii1], e_table, NULL) );
		glPopName();
	}

	/* translucent end */
	if(obj->d_prop->line_trans_sw) unset_gl_translucent();

	return(NORMAL_RC);
}


static RC
drow_fem_surf_point (GTKU_OBJECT *obj)
{
	int ii1, ii2;

	RC_NULL_CHK(obj);
	RC_NULL_CHK(obj->d_prop);
	RC_NULL_CHK(obj->node);
	RC_NULL_CHK(obj->surf);

	glDisable(GL_LIGHTING);

	glColor4dv(obj->d_prop->point_color);
	glPointSize((float)obj->d_prop->point_size);
	if(obj->d_prop->point_trans_sw) set_gl_translucent();
	for(ii1=0; ii1<obj->surf->size; ii1++){
		ELEMENT surf = obj->surf->array[ii1];
		if(surf.label < 0) continue;
		for(ii2=0; ii2<surf.node_num; ii2++){
			int node_index = search_node_label(*(obj)->node,surf.node[ii2]);
			RC_NEG_CHK(node_index);
			glPushName(obj->node->array[node_index].label);
			glBegin(GL_POINTS);
			glVertex3d( obj->node->array[node_index].p.x,
			            obj->node->array[node_index].p.y,
			            obj->node->array[node_index].p.z );
			glEnd();
			glPopName();
		}
	}
	if(obj->d_prop->point_trans_sw) unset_gl_translucent();

	return(NORMAL_RC);
}


static RC
drow_fem_grid_point (GTKU_OBJECT *obj)
{
	int ii1;

	RC_NULL_CHK(obj);
	RC_NULL_CHK(obj->node);
	RC_NULL_CHK(obj->d_prop);

	glDisable(GL_LIGHTING);

	glColor4dv(obj->d_prop->point_color);
	glPointSize((float)obj->d_prop->point_size);
	if(obj->d_prop->point_trans_sw) set_gl_translucent();
	for(ii1=0; ii1<obj->node->size; ii1++){
		NODE node = obj->node->array[ii1];
		if(node.label < 0) continue;
		glPushName(node.label);
		glBegin(GL_POINTS);
		glVertex3d( node.p.x, node.p.y, node.p.z );
		glEnd();
		glPopName();
	}
	if(obj->d_prop->point_trans_sw) unset_gl_translucent();

	return(NORMAL_RC);
}


static RC
drow_fem_hiddenline (GTKU_OBJECT *obj)
{
	int ii1, ii2;
	GLdouble color[MAX_NODE][4];
	FACE_TABLE f_table;
	ELEMENT_ARRAY *surf = NULL;

	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	RC_NULL_CHK(obj);
	RC_NULL_CHK(obj->d_prop);
	if(obj->dim == GTKU_OBJECT_3D){
		RC_NULL_CHK(surf = obj->surf);
	}else if(obj->dim == GTKU_OBJECT_2D){
		RC_NULL_CHK(surf = obj->elem);
	}else{
		GTKU_WARNING("drow_fem_hiddenline()", "Illegal Dimension");
	}

	glDisable(GL_LIGHTING);

	init_face_table(&f_table);
	for(ii1=0; ii1<surf->size; ii1++){
		ELEMENT surface;
		if(surf->array[ii1].label < 0) continue;
		surface = surf->array[ii1];
		glPolygonMode(GL_FRONT, GL_LINE);
		RC_TRY( make_face_table(surface.type, &f_table) );
		RC_TRY( draw_face(obj, surface, f_table, NULL, NULL) );

		/* translucent start */
		if(obj->d_prop->polygon_trans_sw) set_gl_translucent();

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
		RC_TRY( make_face_table(surface.type, &f_table) );
		for(ii2=0; ii2<MAX_NODE; ii2++){
			color[ii2][0] = gtku_gl->bg_color0[0];
			color[ii2][1] = gtku_gl->bg_color0[1];
			color[ii2][2] = gtku_gl->bg_color0[2];
			color[ii2][3] = obj->d_prop->polygon_color[3];
		}
		glPushName(surface.label);
		RC_TRY( draw_face(obj, surface, f_table, color, NULL) );
		glPopName();
		glDisable(GL_POLYGON_OFFSET_FILL);

		/* translucent end */
		if(obj->d_prop->polygon_trans_sw) unset_gl_translucent();
	}

	return(NORMAL_RC);
}


/* Drawing node number, element number, surface element number */
static RC
draw_fem_number (GTKU_OBJECT *obj, GTKU_DRAW_TYPE draw_type)
{
	if(obj->d_prop->draw_node_num_sw){
		if(draw_type & (GTKU_DRAW_POLYGON    | GTKU_DRAW_WIRE_FRAME |
		                GTKU_DRAW_HIDDENLINE | GTKU_DRAW_SURF_POINT)){
			if(gtku_object_exist_chk(obj, 3, "node", "element", "surface"))
				RC_TRY( draw_fem_number_surf_point(obj) );
		}
		if(draw_type & (GTKU_DRAW_SOLID | GTKU_DRAW_GRID_FRAME |
		                GTKU_DRAW_GRID_POINT) ){
			if(gtku_object_exist_chk(obj, 2, "node", "element"))
				RC_TRY( draw_fem_number_grid_point(obj) );
		}
	}

	if(obj->d_prop->draw_elem_num_sw){
		if(gtku_object_exist_chk(obj, 2, "node", "element"))
		RC_TRY( draw_fem_number_solid(obj) );
	}

	if(obj->d_prop->draw_surf_num_sw){
		if(gtku_object_exist_chk(obj, 3, "node", "element", "surface"))
		RC_TRY( draw_fem_number_polygon(obj) );
	}

	return(NORMAL_RC);
}


static RC
draw_fem_number_surf_point (GTKU_OBJECT *obj)
{
	int ii1, ii2;
	char buf[32];

	RC_NULL_CHK(obj);
	RC_NULL_CHK(obj->node);
	RC_NULL_CHK(obj->surf);
	RC_NULL_CHK(obj->d_prop);
	glDisable(GL_LIGHTING);
	glColor4dv(obj->d_prop->text_color);
	if(obj->d_prop->text_trans_sw) set_gl_translucent();
	for(ii1=0; ii1<obj->surf->size; ii1++){
		ELEMENT surf = obj->surf->array[ii1];
		if(surf.label < 0) continue;
		for(ii2=0; ii2<surf.node_num; ii2++){
			int node_index;
			node_index=search_node_label(*(obj)->node, surf.node[ii2]);
			RC_NEG_CHK(node_index);
			g_snprintf(buf, sizeof(buf),"%d",
			           obj->node->array[node_index].label);
			buf[31] = '\0'; /* 念のため */
			gtku_gl_draw_string(obj->node->array[node_index].p, buf);
		}
	}
	if(obj->d_prop->text_trans_sw) unset_gl_translucent();
	return(NORMAL_RC);
}


static RC
draw_fem_number_grid_point (GTKU_OBJECT *obj)
{
	int ii1;
	char buf[32];

	RC_NULL_CHK(obj);
	RC_NULL_CHK(obj->node);
	RC_NULL_CHK(obj->d_prop);
	glDisable(GL_LIGHTING);
	glColor4dv(obj->d_prop->text_color);
	if(obj->d_prop->text_trans_sw) set_gl_translucent();
	for(ii1=0; ii1<obj->node->size; ii1++){
		NODE node = obj->node->array[ii1];
		if(node.label < 0) continue;
		g_snprintf(buf, sizeof(buf), "%d", node.label);
		buf[31] = '\0'; /* 念のため */
		gtku_gl_draw_string(node.p, buf);
	}
	if(obj->d_prop->text_trans_sw) unset_gl_translucent();
	return(NORMAL_RC);
}


static RC
draw_fem_number_solid (GTKU_OBJECT *obj)
{
	int ii1, ii2;
	char buf[32];

	RC_NULL_CHK(obj);
	RC_NULL_CHK(obj->node);
	RC_NULL_CHK(obj->elem);
	RC_NULL_CHK(obj->d_prop);

	glColor4dv(obj->d_prop->text_color);
	if(obj->d_prop->text_trans_sw) set_gl_translucent();
	for(ii1=0; ii1<obj->elem->size; ii1++){
		ELEMENT elem = obj->elem->array[ii1];
		VECT3D center;
		RC_TRY( init_vect3d(&center) );
		if(elem.label < 0) continue;
		for(ii2=0; ii2<elem.node_num; ii2++){
			int node_index;
			node_index=search_node_label(*(obj)->node, elem.node[ii2]);
			RC_NEG_CHK(node_index);
			center.x += obj->node->array[node_index].p.x; 
			center.y += obj->node->array[node_index].p.y; 
			center.z += obj->node->array[node_index].p.z; 
		}
		center.x /= (double)elem.node_num;
		center.y /= (double)elem.node_num;
		center.z /= (double)elem.node_num;
		g_snprintf(buf, sizeof(buf), "%d", elem.label);
		buf[31] = '\0'; /* 念のため */
		gtku_gl_draw_string(center, buf);
	}
	if(obj->d_prop->text_trans_sw) unset_gl_translucent();

	return(NORMAL_RC);
}


static RC
draw_fem_number_polygon (GTKU_OBJECT *obj)
{
	int ii1, ii2;
	char buf[32];

	RC_NULL_CHK(obj);
	RC_NULL_CHK(obj->node);
	RC_NULL_CHK(obj->surf);
	RC_NULL_CHK(obj->d_prop);

	glColor4dv(obj->d_prop->text_color);
	if(obj->d_prop->text_trans_sw) set_gl_translucent();
	for(ii1=0; ii1<obj->surf->size; ii1++){
		ELEMENT surf = obj->surf->array[ii1];
		VECT3D center;
		RC_TRY( init_vect3d(&center) );
		if(surf.label < 0) continue;
		for(ii2=0; ii2<surf.node_num; ii2++){
			int node_index;
			node_index=search_node_label(*(obj)->node, surf.node[ii2]);
			RC_NEG_CHK(node_index);
			center.x += obj->node->array[node_index].p.x; 
			center.y += obj->node->array[node_index].p.y; 
			center.z += obj->node->array[node_index].p.z; 
		}
		center.x /= (double)surf.node_num;
		center.y /= (double)surf.node_num;
		center.z /= (double)surf.node_num;
		g_snprintf(buf, sizeof(buf), "%d", surf.label);
		buf[31] = '\0'; /* 念のため */
		gtku_gl_draw_string(center, buf);
	}
	if(obj->d_prop->text_trans_sw) set_gl_translucent();

	return(NORMAL_RC);
}

