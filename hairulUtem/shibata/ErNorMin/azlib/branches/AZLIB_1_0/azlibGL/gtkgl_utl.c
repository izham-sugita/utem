/*********************************************************************
 * gtkgl_utl.c
 *
 * Copyright (C) 2002 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA>
 *
 * Thanks to following contributors.
 *  <Kenzen TAKEUCHI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: gtkgl_utl.c,v 1.12 2003/12/18 04:26:19 sasaoka Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <glib.h>
#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>
#include <gtk/gtkgl.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include "gtk_utl.h"
#include "fem_struct.h"
#include "math_utl.h"
#include "rc.h"

#define TIMEOUT_INTERVAL (10)
#define ANIMATE_THRESHOLD (100.0)
#define ANIMETION_LIMIT_THRESHOLD (5.0)
#define MOVE_THRESHOLD (3.0)

static gchar font_name[] = "courier 12";  /* Fix!! 決め打ち */

/* Global variable */
extern GTKU_GL *gtku_gl;

/* TOP Level OpenGL function */
static void set_gl_projection_matrix(void);
static void set_gl_model_view_matirx(GTKU_OBJECT *obj);
static void copy_form_matrix(double **src, double dest[4][4]);
static void set_gl_light_env(void);
static void set_gl_translucent(void);
static void unset_gl_translucent(void);
static void draw_coodinate_axis(void);
static void draw_center_marker(TRANSLATION3D center, double length);
static RC draw_world_plane(void);
static RC draw_object(GTKU_OBJECT *obj);

/* FEM Object */
static RC draw_fem_object(GTKU_OBJECT *obj);
static RC draw_fem_object_1D(GTKU_OBJECT *obj, GTKU_DRAW_TYPE draw_type);
static RC draw_fem_object_2D(GTKU_OBJECT *obj, GTKU_DRAW_TYPE draw_type);
static RC draw_fem_object_3D(GTKU_OBJECT *obj, GTKU_DRAW_TYPE draw_type);
static RC draw_edge(GTKU_OBJECT *obj, FEM_ELEMENT elem,
                    FEM_EDGE_TABLE e_table, GLdouble color[][4]);
static RC draw_face(GTKU_OBJECT *obj, FEM_ELEMENT elem,
                    FEM_FACE_TABLE f_table, GLdouble color[][4],
                    TRANSLATION3D nvect[]);
static RC draw_face_single(GTKU_OBJECT *obj, FEM_ELEMENT elem,
                           FEM_FACE_TABLE f_table, GLdouble color[][4]);
static RC cal_normal_vector(FEM_NODE_ARRAY node, FEM_ELEMENT elem,
                            TRANSLATION3D *nvect);
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

/* Image Object */
static RC draw_image_object(GTKU_OBJECT *obj);
static RC draw_pixel_texture(GRAPHIC *gr, GTKU_DRAW_PROP *d_prop);
static RC init_texture(GRAPHIC *gr, GTKU_DRAW_PROP *d_prop);
static void trans_rgb24_GLimage(GRAPHIC *gr,  GTKU_DRAW_PROP *d_prop,
                                GLubyte **image);
static RC init_image_direction(GRAPHIC *gr, GTKU_DRAW_PROP *d_prop);
static RC draw_pixel_cell(GRAPHIC *gr, GTKU_DRAW_PROP *d_prop);
static void draw_pixel_cell_single(int height, int x, int y);
static RC draw_voxel_cell(GRAPHIC *gr, GTKU_DRAW_PROP *d_prop);
static void draw_voxel_cell_single(int x, int y, int z);

/* GLAREA, Interface, Event ... */
static RC init_objects_drawing_attribute(void);
static RC init_objects_operation_matrix(void);
static RC init_view_volume(void);
static RC cal_objects_size(GTKU_OBJECT_ARRAY *object, TRANSLATION3D *center,
                           TRANSLATION3D *max, TRANSLATION3D *min);
static RC cal_object_size_fem(FEM_NODE_ARRAY *node, TRANSLATION3D *center,
                              TRANSLATION3D *max, TRANSLATION3D *min);
static RC cal_object_size_graphic(GRAPHIC *gr, TRANSLATION3D *center,
                                  TRANSLATION3D *max, TRANSLATION3D *min);
static RC cal_objects_view_scale(double *scale);

static RC glarea_mouse_rotation(GtkWidget *widget, GdkRectangle area);
static RC mouse_rotation_local(GTKU_DRAW_PROP *d_prop, double **rot);
static RC glarea_mouse_shift(GtkWidget *widget, GdkRectangle area);
static RC mouse_shift_local(GTKU_DRAW_PROP *d_prop, TRANSLATION3D disp);
static RC glarea_mouse_scale(GtkWidget *widget, GdkRectangle area,
                             GdkModifierType state);
static RC mouse_scale_local(GTKU_DRAW_PROP *d_prop, TRANSLATION3D rate);
static gboolean add_rotation_matrix_for_animation(GtkWidget *widget);
static gboolean glarea_key_scale(GtkWidget *widget, GdkEventKey *event);
static gboolean timeout(GtkWidget *widget);
static gboolean timeout_add(GtkWidget *widget);
static gboolean timeout_remove(GtkWidget *widget);


/* The "expose_event" signal handler. This is repeatedly called as the */
/* painting routine every time the "expose/draw" event is signalled.   */
gboolean gtku_gl_draw(GtkWidget *widget, GdkEventExpose *event)
{
	int ii1;
	GdkGLContext *glcontext = gtk_widget_get_gl_context(widget);
	GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable(widget);

	/* Draw only last expose. */
	if(event->count > 0) return(TRUE); 
	
	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

#ifdef GTKU_VERBOSE_CHECK
	/* OpenGL functions can be called only if make_current returns true */
	if(!gdk_gl_drawable_make_current(gldrawable, glcontext)){
		g_warning("glarea not drawable\n");
		return(FALSE);
	}
#endif

	/* OpenGL BEGIN */
	if(gdk_gl_drawable_gl_begin(gldrawable, glcontext)){
		/* lighting initialize */
		if(!gtku_gl->init_light){
			set_gl_light_env();
			gtku_gl->init_light = TRUE;
		}

		/* initialize color and depth buffer. clean screen */
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		/* Animation */
		if(gtku_gl->animation && !gtku_gl->anim_lock_sw)
			add_rotation_matrix_for_animation(widget);

		/* PROJECTION MATRIX */
		set_gl_projection_matrix();

		/* Drawing Object  */
		for(ii1=0;ii1<gtku_gl->object->size;ii1++){
			double length;
			GTKU_OBJECT *obj = gtku_gl->object->obj[ii1];
			set_gl_model_view_matirx(obj);
			GTKU_BOOLEAN_RC_TRY( draw_object(obj) );
			length = (  abs_translation3d(obj->d_prop->max)
					  - abs_translation3d(obj->d_prop->min) ) / 15.0;
			draw_center_marker(obj->d_prop->center, length);
		}

		/* Coordinate System */
		if(gtku_gl->draw_coord_sw) draw_coodinate_axis(); 

		/* Worl Plane */
		if(gtku_gl->draw_plane_sw){
			set_gl_model_view_matirx(NULL);
			RC_TRY( draw_world_plane() );
		}

		/* Swap backbuffer to front */
		if(gdk_gl_drawable_is_double_buffered(gldrawable)){
			gdk_gl_drawable_swap_buffers(gldrawable);
		}else{
			glFlush();
		}

		gdk_gl_drawable_gl_end(gldrawable);
	}/* OpenGL END */

	return TRUE;
}


/* projection matrix setting */
static void set_gl_projection_matrix(void)
{
	/* GTKU_GL existence check */
	g_return_if_fail(gtku_gl);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-(gtku_gl->view_size.x), (gtku_gl->view_size.x),
	        -(gtku_gl->view_size.y), (gtku_gl->view_size.y),
	        -(gtku_gl->view_size.z), (gtku_gl->view_size.z));
}


/* view matrix transformation */
static void set_gl_model_view_matirx(GTKU_OBJECT *obj)
{
	double m[4][4];
	double **mat;

	/* GTKU_GL existence check */
	g_return_if_fail(gtku_gl);

	GTKU_VOID_RC_TRY( allocate_operation_matrix(&mat) );

	/* local matrix operation */
	if(obj != NULL){
		mul_operation_matrix(obj->d_prop->rot_mat, mat);
		mul_operation_matrix(obj->d_prop->shift_mat, mat);
		mul_operation_matrix(obj->d_prop->scale_mat, mat);
	}

	/* global matrix operation */
	mul_operation_matrix(gtku_gl->initial_mat, mat);
	mul_operation_matrix(gtku_gl->rot_mat,     mat);
	mul_operation_matrix(gtku_gl->shift_mat,   mat);
	mul_operation_matrix(gtku_gl->scale_mat,   mat);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	copy_form_matrix(mat, m);
	glMultMatrixd(*m);

	free_operation_matrix(mat);
}


static void copy_form_matrix(double **src, double dest[4][4])
{
	/* 行列の並びに注意!!               */
	/* OpenGL Programing Guide P104参照 */
	int ii1, ii2;
	for(ii1=0;ii1<4;ii1++){
		for(ii2=0;ii2<4;ii2++){
			dest[ii2][ii1] = src[ii1][ii2];
		}
	}
}


/* object lighting and background color */
static void set_gl_light_env()
{
	/* GTKU_GL existence check */
	g_return_if_fail(gtku_gl);

	/* first light settting. lighting 0 */
	glLightfv(GL_LIGHT0, GL_POSITION, gtku_gl->light0_pos);
	glLightfv(GL_LIGHT0, GL_DIFFUSE,  gtku_gl->light0_color); 
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);

	/* remove back faces */
	glFrontFace(GL_CCW); /* 反時計回りのポリゴンを表と認識 */
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);

	/* speedups */
	glEnable(GL_DITHER);
	glShadeModel(GL_SMOOTH);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_FASTEST);
	glHint(GL_POINT_SMOOTH_HINT,           GL_FASTEST);
	glHint(GL_LINE_SMOOTH_HINT,            GL_FASTEST);
	glHint(GL_POLYGON_SMOOTH_HINT,         GL_FASTEST);
	/*
	 * sample for Anti-Aliasing. It is "TOO SLOW".
	 * glHint(GL_POINT_SMOOTH_HINT,   GL_NICEST);
	 * glHint(GL_LINE_SMOOTH_HINT,    GL_NICEST);
	 * glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
	 * glEnable(GL_POINT_SMOOTH);
	 * glEnable(GL_LINE_SMOOTH);
	 * glEnable(GL_POLYGON_SMOOTH);
	 */

	glDepthFunc(GL_LEQUAL);
	glEnable(GL_DEPTH_TEST);

	glMaterialfv(GL_FRONT, GL_SPECULAR,  gtku_gl->material_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, gtku_gl->material_shininess);
	glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);  

	/* background color */
	glClearColor(gtku_gl->bg_color[0], gtku_gl->bg_color[1],
                 gtku_gl->bg_color[2], gtku_gl->bg_color[3]);
	 
	/* 拡大縮小などの Matrix 操作を行うとき、OpenGL 側に法線ベクトルの */
	/* 単位化をさせる必要がある。glEnable(GL_NORMALIZE) は自動正規化を */
	/* ONにしている                                                    */
	glEnable(GL_NORMALIZE);

	/* GLUのTeapotなどで使われるNUBAS局面のの法線を自動計算する */
	glEnable(GL_AUTO_NORMAL);
}


static void set_gl_translucent(void)
{
	glEnable(GL_BLEND);
	glDepthMask(FALSE);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}


static void unset_gl_translucent(void)
{
	glDepthMask(TRUE);
	glDisable(GL_BLEND);
}


static void draw_coodinate_axis(void)
{
	double m[4][4];
	double **mat;
	TRANSLATION3D p; /* Position of the coordinate system Name */

	p.x = -4.0;
	p.y = -4.0;
	p.z = 60.0;
	
	/* GTKU_GL existence check */
	g_return_if_fail(gtku_gl);

	GTKU_VOID_RC_TRY( allocate_operation_matrix(&mat) );
	mul_operation_matrix(gtku_gl->rot_mat, mat);
	copy_form_matrix(mat, m);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(-(gtku_gl->view_size.x - 67),
	             -(gtku_gl->view_size.y - 67), 0.0);

	/* Origin of the coordinate system. Sphere */
	glPushMatrix();
	glMultMatrixd(*m);
	glEnable(GL_LIGHTING);
	glColor4d(0.7, 0.7, 0.7, 1.0);  /* Gray */
	gdk_gl_draw_sphere(TRUE, 12.0, 10, 10);
	glPopMatrix ();

	/* X-axis CONE */
	glPushMatrix();
	glMultMatrixd(*m);
	glTranslatef(5.0, 0.0, 0.0);
	glRotated(90.0, 0.0, 1.0, 0.0);
	glColor4d(0.8, 0.0, 0.0, 1.0);  /* Red */
	glEnable(GL_LIGHTING);
	gdk_gl_draw_cone(TRUE, 8.0, 50, 10, 10);
	glColor4d(1.0, 1.0, 1.0, 1.0);  /* White */
	gtku_gl_draw_string(p, "X");
	glPopMatrix();

	/* Y-axis CONE */
	glPushMatrix();
	glMultMatrixd(*m);
	glTranslatef(0.0, 5.0, 0.0);
	glRotated(-90.0, 1.0, 0.0, 0.0);
	glColor4d(0.0, 0.8, 0.0, 1.0);  /* Green */
	glEnable(GL_LIGHTING);
	gdk_gl_draw_cone(TRUE, 8.0, 50, 10, 10);
	glColor4f(1.0, 1.0, 1.0, 1.0);  /* White */
	gtku_gl_draw_string(p, "Y");
	glPopMatrix();

	/* Z-axis CONE */
	glPushMatrix();
	glMultMatrixd(*m);
	glTranslatef(0.0, 0.0, 5.0);
	glColor4d(0.0, 0.0, 0.8, 1.0);  /* Blue */
	glEnable(GL_LIGHTING);
	gdk_gl_draw_cone(TRUE, 8.0, 50, 10, 10);
	glColor4d(1.0, 1.0, 1.0, 1.0);  /* White */
	gtku_gl_draw_string(p, "Z");
	glPopMatrix();

	free_operation_matrix(mat);
}


static void draw_center_marker(TRANSLATION3D center, double length)
{
	glColor4d(1.0, 1.0, 0.0, 1.0);
	glDisable(GL_LIGHTING);
	glBegin(GL_LINES);
	glVertex3d(center.x-length, center.y,        center.z       );
	glVertex3d(center.x+length, center.y,        center.z       );
	glVertex3d(center.x,        center.y-length, center.z       );
	glVertex3d(center.x,        center.y+length, center.z       );
	glVertex3d(center.x,        center.y,        center.z-length);
	glVertex3d(center.x,        center.y,        center.z+length);
	glEnd();
}


#define WORLD_PLANE_DIVIDE (4) /* 格子数 */
#define WORLD_PLANE_FACTER (0.5)
static RC draw_world_plane(void)
{
	double ii1, ii2, ii3, max_dist, d;
	TRANSLATION3D p, min, max;

	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);
	
	min = gtku_gl->object->min;
	max = gtku_gl->object->max;
	p = (abs_translation3d(max) > abs_translation3d(min) ? max : min);
	max_dist = MAX2(MAX2(fabs(p.x), fabs(p.y)), fabs(p.z));
	d = max_dist/WORLD_PLANE_DIVIDE;

	glDisable(GL_LIGHTING);
	for(ii1=0; ii1<2; ii1++){
		if(ii1==0){ /* Plane */
			set_gl_translucent();
			glColor4d(0.7, 0.7, 0.7, 0.3);  /* Gray */
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glDisable(GL_CULL_FACE);
		}else{ /* Lince */
			glColor4d(0.8, 0.8, 0.8, 1.0);  /* Gray */
			glLineWidth(1.0);
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		}
		for(ii2=-max_dist; ii2<max_dist-d/2.0; ii2+=d){
			for(ii3=-max_dist; ii3<max_dist-d/2.0; ii3+=d){
				if(gtku_gl->draw_plane_x_sw){
					glBegin(GL_POLYGON);
						glVertex4d(0.0, ii2,   ii3,   WORLD_PLANE_FACTER);
						glVertex4d(0.0, ii2+d, ii3,   WORLD_PLANE_FACTER);
						glVertex4d(0.0, ii2+d, ii3+d, WORLD_PLANE_FACTER);
						glVertex4d(0.0, ii2,   ii3+d, WORLD_PLANE_FACTER);
					glEnd();
				}
				if(gtku_gl->draw_plane_y_sw){
					glBegin(GL_POLYGON);
						glVertex4d(ii2,   0.0, ii3,   WORLD_PLANE_FACTER);
						glVertex4d(ii2,   0.0, ii3+d, WORLD_PLANE_FACTER);
						glVertex4d(ii2+d, 0.0, ii3+d, WORLD_PLANE_FACTER);
						glVertex4d(ii2+d, 0.0, ii3,   WORLD_PLANE_FACTER);
					glEnd();
				}
				if(gtku_gl->draw_plane_z_sw){
					glBegin(GL_POLYGON);
						glVertex4d(ii2,   ii3,   0.0, WORLD_PLANE_FACTER);
						glVertex4d(ii2+d, ii3,   0.0, WORLD_PLANE_FACTER);
						glVertex4d(ii2+d, ii3+d, 0.0, WORLD_PLANE_FACTER);
						glVertex4d(ii2,   ii3+d, 0.0, WORLD_PLANE_FACTER);
					glEnd();
				}
			}
		}
		if(ii1==0){
			unset_gl_translucent();
			glCullFace(GL_BACK);
			glEnable(GL_CULL_FACE);
		}else{
			glColor4d(1.0, 1.0, 1.0, 1.0);  /* White */
			glLineWidth(2.0);
			if(gtku_gl->draw_plane_x_sw){
				glBegin(GL_LINES);
					glVertex4d(0.0, -max_dist, 0.0, WORLD_PLANE_FACTER);
					glVertex4d(0.0,  max_dist, 0.0, WORLD_PLANE_FACTER);
					glVertex4d(0.0, 0.0, -max_dist, WORLD_PLANE_FACTER);
					glVertex4d(0.0, 0.0,  max_dist, WORLD_PLANE_FACTER);
				glEnd();
			}
			if(gtku_gl->draw_plane_y_sw){
				glBegin(GL_LINES);
					glVertex4d(-max_dist, 0.0, 0.0, WORLD_PLANE_FACTER);
					glVertex4d( max_dist, 0.0, 0.0, WORLD_PLANE_FACTER);
					glVertex4d(0.0, 0.0, -max_dist, WORLD_PLANE_FACTER);
					glVertex4d(0.0, 0.0,  max_dist, WORLD_PLANE_FACTER);
				glEnd();
			}
			if(gtku_gl->draw_plane_z_sw){
				glBegin(GL_LINES);
					glVertex4d(-max_dist, 0.0, 0.0, WORLD_PLANE_FACTER);
					glVertex4d( max_dist, 0.0, 0.0, WORLD_PLANE_FACTER);
					glVertex4d(0.0, -max_dist, 0.0, WORLD_PLANE_FACTER);
					glVertex4d(0.0,  max_dist, 0.0, WORLD_PLANE_FACTER);
				glEnd();
			}
		}
	}

	return(NORMAL_RC);
}


static RC draw_object(GTKU_OBJECT *obj)
{
	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	if( (obj->type & GTKU_OBJECT_NONE) ){
		return(NORMAL_RC);
	}else if(obj->type & GTKU_OBJECT_FEM_NODE){ /* element does not nessarily */
		RC_TRY( draw_fem_object(obj) );
	}else if(obj->type & (GTKU_OBJECT_IMAGE_PIXEL|GTKU_OBJECT_IMAGE_VOXEL)){
		RC_TRY( draw_image_object(obj) );
	}else{
		GTKU_WARNING("draw_object()", "Unknown Object");
	}

	return(NORMAL_RC);
}


static RC draw_fem_object(GTKU_OBJECT *obj)
{
	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	RC_NULL_CHK(obj);
	GTKU_DRAW_TYPE draw_type = obj->d_prop->draw_type;
	if(draw_type & GTKU_DRAW_HIDE) return(NORMAL_RC);
	if(draw_type & GTKU_DRAW_DEFAULT)
		draw_type |= gtku_gl->draw_type;

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
static RC draw_fem_object_1D(GTKU_OBJECT *obj, GTKU_DRAW_TYPE draw_type)
{
	if(draw_type & (GTKU_DRAW_GRID_POINT | GTKU_DRAW_SURF_POINT)){
		if(gtku_object_exist_chk(obj, 1, "node")){
			RC_TRY(drow_fem_grid_point(obj));
		}else GTKU_WARNING("draw_fem_object_1D()", "none object (Point)");
	}

	if(draw_type & (GTKU_DRAW_POLYGON    | GTKU_DRAW_SOLID      |
	                GTKU_DRAW_GRID_FRAME | GTKU_DRAW_WIRE_FRAME |
	                GTKU_DRAW_HIDDENLINE) ){
		if(gtku_object_exist_chk(obj, 2, "node", "element")){
			RC_TRY(drow_fem_wire_frame(obj))
		}else
			GTKU_WARNING("draw_fem_object_1D()", "none object (other POINT)");
	}

	return(NORMAL_RC);
}


/* 2D(SHELL, TRI, QUAD) */
static RC draw_fem_object_2D(GTKU_OBJECT *obj, GTKU_DRAW_TYPE draw_type)
{
	if(draw_type & GTKU_DRAW_GRID_POINT){
		if( gtku_object_exist_chk(obj, 1, "node") ){
			RC_TRY(drow_fem_grid_point(obj));
		}else GTKU_WARNING("draw_fem_object_2D()", "none object (Grid Point)");
	}

	if(draw_type & GTKU_DRAW_SURF_POINT){
		if( gtku_object_exist_chk(obj, 2, "node", "surface") ){
			RC_TRY(drow_fem_surf_point(obj));
		}else GTKU_WARNING("draw_fem_object_2D()","none object (SurfacePoint)");
	}

	if(draw_type & GTKU_DRAW_GRID_FRAME){
		if( gtku_object_exist_chk(obj, 2, "node", "element") ){
			RC_TRY(drow_fem_grid_frame(obj));
		}else GTKU_WARNING("draw_fem_object_2D()", "none object (Grid Frame)");
	}

	if(draw_type & GTKU_DRAW_WIRE_FRAME){
		if( gtku_object_exist_chk(obj, 2, "node", "surface") ){
			RC_TRY(drow_fem_wire_frame(obj));
		}else GTKU_WARNING("draw_fem_object_2D()", "none object (Wire Frame)");
	}

	if(draw_type & GTKU_DRAW_HIDDENLINE){
		if( gtku_object_exist_chk(obj, 2, "node", "element") ){
			RC_TRY(drow_fem_hiddenline(obj));
		}else GTKU_WARNING("draw_fem_object_2D()", "none object (HiddenLine)");
	}

	if(draw_type & (GTKU_DRAW_POLYGON|GTKU_DRAW_SOLID)){
		if( gtku_object_exist_chk(obj, 2, "node", "element") ){
			RC_TRY(drow_fem_polygon_2D(obj));
		}else GTKU_WARNING("draw_fem_object_2D()","none object (Polygon)");
	}

	return(NORMAL_RC);
}


/* 3D(TETRA, PENTA, HEXA) */
static RC draw_fem_object_3D(GTKU_OBJECT *obj, GTKU_DRAW_TYPE draw_type)
{
	if(draw_type & GTKU_DRAW_GRID_POINT){
		if( gtku_object_exist_chk(obj, 1, "node") ){
			RC_TRY(drow_fem_grid_point(obj));
		}else GTKU_WARNING("draw_fem_object_3D()", "none object (Grid Point)");
	}

	if(draw_type & GTKU_DRAW_SURF_POINT){
		if( gtku_object_exist_chk(obj, 2, "node", "surface") ){
			RC_TRY(drow_fem_surf_point(obj));
		}else GTKU_WARNING("draw_fem_object_3D()","none object (SurfacePoint)");
	}

	if(draw_type & GTKU_DRAW_GRID_FRAME){
		if( gtku_object_exist_chk(obj, 2, "node", "element") ){
			RC_TRY(drow_fem_grid_frame(obj));
		}else GTKU_WARNING("draw_fem_object_3D()", "none object (Grid Frame)");
	}

	if(draw_type & GTKU_DRAW_WIRE_FRAME){
		if( gtku_object_exist_chk(obj, 2, "node", "surface") ){
			RC_TRY(drow_fem_wire_frame(obj));
		}else GTKU_WARNING("draw_fem_object_3D()", "none object(Wire Frame)");
	}

	if(draw_type & GTKU_DRAW_HIDDENLINE){
		if( gtku_object_exist_chk(obj, 3, "node", "element", "surface") ){
			RC_TRY(drow_fem_hiddenline(obj));
		}else GTKU_WARNING("draw_fem_object_3D()", "none object (HiddenLine)");
	}

	if(draw_type & GTKU_DRAW_SOLID){
		if( gtku_object_exist_chk(obj, 2, "node", "element") ){
			RC_TRY(drow_fem_solid_3D(obj));
		}else GTKU_WARNING("draw_fem_object_3D()", "none object (Solid)");
	}

	if(draw_type & GTKU_DRAW_POLYGON){
		if( gtku_object_exist_chk(obj, 3, "node", "element", "surface") ){
			RC_TRY(drow_fem_polygon_3D(obj));
		}else GTKU_WARNING("draw_fem_object_3D()","none object (Polygon)");
	}

	return(NORMAL_RC);
}


static RC draw_edge(GTKU_OBJECT *obj, FEM_ELEMENT elem,
                    FEM_EDGE_TABLE e_table, GLdouble color[][4])
{
	int ii1, ii2;

	if(element_order(elem.type) == 1){
		for(ii1=0; ii1<e_table.edge_num; ii1++){
			glBegin(GL_LINES);
			for(ii2=0; ii2<e_table.node_num[ii1]; ii2++){
				int label = elem.node[(e_table.table[ii1][ii2])];
				int node_index = search_fem_node_label(*(obj->node), label);
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
			int node_index0 = search_fem_node_label(*(obj->node), label0);
			int node_index1 = search_fem_node_label(*(obj->node), label1);
			int node_index2 = search_fem_node_label(*(obj->node), label2);
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
			 TRANSLATION3D mid = mid_point3d(obj->node->array[node_index0].p,
			                                 obj->node->array[node_index1].p);
			 TRANSLATION3D sub =
			     sub_translation3d(obj->node->array[node_index2].p, mid);
			 TRANSLATION3D add =
			     add_translation3d(obj->node->array[node_index2].p, sub);
			 ctlp[1][0] = add.x;
			 ctlp[1][1] = add.y;
			 ctlp[1][2] = add.z;
			}
			
			glColor4dv(obj->d_prop->line_color);
			glMap1d(GL_MAP1_VERTEX_3, 0.0, 1.0, 3, 3, &ctlp[0][0]);
			glEnable(GL_MAP1_VERTEX_3);
			/* なぜか この方法はうまくいかない
			glMapGrid1d(10, 0.0, 1.0);
			glEvalMesh1(GL_LINES, 0, 10);
			*/
			glBegin(GL_LINE_STRIP);
			{
				int ii3;
				for(ii3=0; ii3<= 10; ii3++)
					glEvalCoord1f((GLfloat)ii3/10);
			}
			glEnd();
		}
	}else{
		GTKU_WARNING("draw_edge()", "Illegal element order");
	}

	return(NORMAL_RC);
}


static RC draw_face(GTKU_OBJECT *obj, FEM_ELEMENT elem,
                    FEM_FACE_TABLE f_table, GLdouble color[][4],
                    TRANSLATION3D nvect[])
{
	if(element_dim(elem.type) == 2){
		glBegin(GL_POLYGON);
		if(nvect != NULL)
			glNormal3d(nvect->x, nvect->y, nvect->z);
		RC_TRY( draw_face_single(obj, elem, f_table, color) );
		glEnd();
	}else if(element_dim(elem.type) == 3){
		int ii1, ii2;
		FEM_ELEMENT local_elem;
		FEM_FACE_TABLE local_f_table;
		for(ii1=0; ii1<f_table.face_num; ii1++){
			local_elem.type =  f_table.face_type[ii1];
			RC_TRY( make_fem_face_table(local_elem.type, &local_f_table) );
			for(ii2=0; ii2<f_table.node_num[ii1]; ii2++)
				local_elem.node[ii2] = elem.node[f_table.table[ii1][ii2]];
			glBegin(GL_POLYGON);
			if(nvect == NULL){
				TRANSLATION3D l_nvect;
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


static RC draw_face_single(GTKU_OBJECT *obj, FEM_ELEMENT elem,
                           FEM_FACE_TABLE f_table, GLdouble color[][4])
{
	int ii1, ii2;

	/* Fix!! 2次要素はこれでは表示できない */
	for(ii1=0; ii1<f_table.face_num; ii1++){
		/* for(ii2=0; ii2<f_table.node_num[ii1]; ii2++){ */
		for(ii2=0; ii2<2; ii2++){
			int label = elem.node[ f_table.table[ii1][ii2] ];
			int node_index = search_fem_node_label(*(obj->node), label);
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
static RC cal_normal_vector(FEM_NODE_ARRAY node, FEM_ELEMENT elem,
                            TRANSLATION3D *nvect)
{
	int node_index0, node_index1, node_index2;
	TRANSLATION3D v1, v2;

	if( !((elem.type == ELEM_TRI1)  ||
	      (elem.type == ELEM_TRI2)  ||
	      (elem.type == ELEM_QUAD1) ||
	      (elem.type == ELEM_QUAD2)) ) return(ARG_ERROR_RC);

	RC_NEG_CHK( node_index0 = search_fem_node_label(node, elem.node[0]) );
	RC_NEG_CHK( node_index1 = search_fem_node_label(node, elem.node[1]) );
	RC_NEG_CHK( node_index2 = search_fem_node_label(node, elem.node[2]) );
	v1 = sub_translation3d(node.array[node_index1].p,node.array[node_index0].p);
	v2 = sub_translation3d(node.array[node_index2].p,node.array[node_index0].p);
	*nvect = outer_product3d(v1, v2);
	return(NORMAL_RC);
}


static RC drow_fem_polygon_2D(GTKU_OBJECT *obj)
{
	int ii1;
	FEM_FACE_TABLE f_table;

	RC_NULL_CHK(obj);
	RC_NULL_CHK(obj->node);
	RC_NULL_CHK(obj->elem);
	RC_NULL_CHK(obj->d_prop);

	/* translucent start */
	if(obj->d_prop->polygon_trans_sw) set_gl_translucent();

	glEnable(GL_LIGHTING);

	init_fem_face_table(&f_table);
	for(ii1=0; ii1<obj->elem->size; ii1++){
		TRANSLATION3D nvect[1];
		FEM_ELEMENT elem = obj->elem->array[ii1];
		if(elem.label < 0) continue;
		/* Cal and Set Normal Vector */
		RC_TRY( cal_normal_vector(*(obj)->node, elem, &nvect[0]) );
		RC_TRY( make_fem_face_table(elem.type, &f_table) );
		/* Drawing polygon */
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		RC_TRY( draw_face(obj, elem, f_table, NULL, nvect) );
	}

	/* translucent end */
	if(obj->d_prop->polygon_trans_sw) unset_gl_translucent();

	return(NORMAL_RC);
}


static RC drow_fem_polygon_3D(GTKU_OBJECT *obj)
{
	int ii1;
	int elem_index;
	FEM_FACE_TABLE f_table;
	TRANSLATION3D nvect[1];

	RC_NULL_CHK(obj);
	RC_NULL_CHK(obj->node);
	RC_NULL_CHK(obj->elem);
	RC_NULL_CHK(obj->surf);
	RC_NULL_CHK(obj->d_prop);

	/* translucent start */
	if(obj->d_prop->polygon_trans_sw) set_gl_translucent();

	glEnable(GL_LIGHTING);

	init_fem_face_table(&f_table);
	for(ii1=0; ii1<obj->surf->size; ii1++){
		FEM_ELEMENT surf = obj->surf->array[ii1];
		if(surf.label < 0) continue;
		/* Cal and Set Normal Vector */
		elem_index = search_fem_element_label(*(obj)->elem, surf.i_info[0]);
		RC_NEG_CHK(elem_index);
		RC_TRY( unit_normal_vect_center(surf, obj->elem->array[elem_index],
		                                *(obj)->node, &nvect[0]) );
		RC_TRY( make_fem_face_table(surf.type, &f_table) );
		/* Drawing polygon */
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		RC_TRY( draw_face(obj, surf, f_table, NULL, nvect) );
	}

	/* translucent end */
	if(obj->d_prop->polygon_trans_sw) unset_gl_translucent();

	return(NORMAL_RC);
}


static RC drow_fem_solid_3D(GTKU_OBJECT *obj)
{
	int ii1;
	FEM_FACE_TABLE f_table;

	RC_NULL_CHK(obj);
	RC_NULL_CHK(obj->node);
	RC_NULL_CHK(obj->elem);
	RC_NULL_CHK(obj->d_prop);

	/* translucent start */
	if(obj->d_prop->polygon_trans_sw) set_gl_translucent();

	glEnable(GL_LIGHTING);

	init_fem_face_table(&f_table);
	for(ii1=0; ii1<obj->elem->size; ii1++){
		FEM_ELEMENT elem = obj->elem->array[ii1];;
		if(elem.label < 0) continue;
		RC_TRY( make_fem_face_table(elem.type, &f_table) );
		/* Drawing polygon */
		/* 法線に関してはとりあえず簡易的な方法で対処 */
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		RC_TRY( draw_face(obj, elem, f_table, NULL, NULL) );
	}

	/* translucent end */
	if(obj->d_prop->polygon_trans_sw) unset_gl_translucent();

	return(NORMAL_RC);
}


static RC drow_fem_wire_frame(GTKU_OBJECT *obj)
{
	int ii1;
	FEM_EDGE_TABLE e_table;

	RC_NULL_CHK(obj);
	RC_NULL_CHK(obj->surf);
	RC_NULL_CHK(obj->d_prop);

	/* translucent start */
	if(obj->d_prop->line_trans_sw) set_gl_translucent();

	glDisable(GL_LIGHTING);

	init_fem_edge_table(&e_table);
	for(ii1=0; ii1<obj->surf->size; ii1++){
		if(obj->surf->array[ii1].label < 0) continue;
		glLineWidth(obj->d_prop->line_size);
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
			RC_TRY( make_fem_edge_table(obj->surf->array[ii1].type, &e_table) );
		}
		RC_TRY( draw_edge(obj, obj->surf->array[ii1], e_table, NULL) );
	}

	/* translucent end */
	if(obj->d_prop->line_trans_sw) unset_gl_translucent();

	return(NORMAL_RC);
}


static RC drow_fem_grid_frame(GTKU_OBJECT *obj)
{
	int ii1;
	FEM_EDGE_TABLE e_table;

	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	RC_NULL_CHK(obj);
	RC_NULL_CHK(obj->elem);
	RC_NULL_CHK(obj->d_prop);

	/* translucent start */
	if(obj->d_prop->line_trans_sw) set_gl_translucent();

	glDisable(GL_LIGHTING);

	init_fem_edge_table(&e_table);
	for(ii1=0; ii1<obj->elem->size; ii1++){
		if(obj->elem->array[ii1].label < 0) continue;
		glLineWidth(obj->d_prop->line_size);
		RC_TRY( make_fem_edge_table(obj->elem->array[ii1].type, &e_table) );
		RC_TRY( draw_edge(obj, obj->elem->array[ii1], e_table, NULL) );
	}

	/* translucent end */
	if(obj->d_prop->line_trans_sw) unset_gl_translucent();

	return(NORMAL_RC);
}


static RC drow_fem_surf_point(GTKU_OBJECT *obj)
{
	int ii1, ii2;

	RC_NULL_CHK(obj);
	RC_NULL_CHK(obj->d_prop);
	RC_NULL_CHK(obj->node);
	RC_NULL_CHK(obj->surf);

	glDisable(GL_LIGHTING);

	glColor4dv(obj->d_prop->point_color);
	glPointSize(obj->d_prop->point_size);
	if(obj->d_prop->point_trans_sw) set_gl_translucent();
	for(ii1=0; ii1<obj->surf->size; ii1++){
		FEM_ELEMENT surf = obj->surf->array[ii1];
		if(surf.label < 0) continue;
		for(ii2=0; ii2<surf.node_num; ii2++){
			int node_index = search_fem_node_label(*(obj)->node,surf.node[ii2]);
			RC_NEG_CHK(node_index);
			glBegin(GL_POINTS);
			glVertex3d( obj->node->array[node_index].p.x,
			            obj->node->array[node_index].p.y,
			            obj->node->array[node_index].p.z );
			glEnd();
		}
	}
	if(obj->d_prop->point_trans_sw) unset_gl_translucent();

	return(NORMAL_RC);
}


static RC drow_fem_grid_point(GTKU_OBJECT *obj)
{
	int ii1;

	RC_NULL_CHK(obj);
	RC_NULL_CHK(obj->node);
	RC_NULL_CHK(obj->d_prop);

	glDisable(GL_LIGHTING);

	glColor4dv(obj->d_prop->point_color);
	glPointSize(obj->d_prop->point_size);
	if(obj->d_prop->point_trans_sw) set_gl_translucent();
	for(ii1=0; ii1<obj->node->size; ii1++){
		FEM_NODE node = obj->node->array[ii1];
		if(node.label < 0) continue;
		glBegin(GL_POINTS);
		glVertex3d( node.p.x, node.p.y, node.p.z );
		glEnd();
	}
	if(obj->d_prop->point_trans_sw) unset_gl_translucent();

	return(NORMAL_RC);
}


static RC drow_fem_hiddenline(GTKU_OBJECT *obj)
{
	int ii1, ii2;
	GLdouble color[FEM_MAX_NODE][4];
	FEM_FACE_TABLE f_table;
	FEM_ELEMENT_ARRAY *surf = NULL;

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

	init_fem_face_table(&f_table);
	for(ii1=0; ii1<surf->size; ii1++){
		if(surf->array[ii1].label < 0) continue;
		FEM_ELEMENT surface = surf->array[ii1];
		glPolygonMode(GL_FRONT, GL_LINE);
		RC_TRY( make_fem_face_table(surface.type, &f_table) );
		RC_TRY( draw_face(obj, surface, f_table, NULL, NULL) );

		/* translucent start */
		if(obj->d_prop->polygon_trans_sw) set_gl_translucent();

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
		RC_TRY( make_fem_face_table(surface.type, &f_table) );
		for(ii2=0; ii2<FEM_MAX_NODE; ii2++){
			color[ii2][0] = gtku_gl->bg_color[0];
			color[ii2][1] = gtku_gl->bg_color[1];
			color[ii2][2] = gtku_gl->bg_color[2];
			color[ii2][3] = obj->d_prop->polygon_color[3];
		}
		RC_TRY( draw_face(obj, surface, f_table, color, NULL) );
		glDisable(GL_POLYGON_OFFSET_FILL);

		/* translucent end */
		if(obj->d_prop->polygon_trans_sw) unset_gl_translucent();
	}

	return(NORMAL_RC);
}


/* Drawing node number, element number, surface element number */
static RC draw_fem_number(GTKU_OBJECT *obj, GTKU_DRAW_TYPE draw_type)
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


static RC draw_fem_number_surf_point(GTKU_OBJECT *obj)
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
		FEM_ELEMENT surf = obj->surf->array[ii1];
		if(surf.label < 0) continue;
		for(ii2=0; ii2<surf.node_num; ii2++){
			int node_index;
			node_index=search_fem_node_label(*(obj)->node, surf.node[ii2]);
			RC_NEG_CHK(node_index);
			snprintf(buf, sizeof(buf),"%d", obj->node->array[node_index].label);
			buf[31] = '\0'; /* 念のため */
			gtku_gl_draw_string(obj->node->array[node_index].p, buf);
		}
	}
	if(obj->d_prop->text_trans_sw) unset_gl_translucent();
	return(NORMAL_RC);
}


static RC draw_fem_number_grid_point(GTKU_OBJECT *obj)
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
		FEM_NODE node = obj->node->array[ii1];
		if(node.label < 0) continue;
		snprintf(buf, sizeof(buf), "%d", node.label);
		buf[31] = '\0'; /* 念のため */
		gtku_gl_draw_string(node.p, buf);
	}
	if(obj->d_prop->text_trans_sw) unset_gl_translucent();
	return(NORMAL_RC);
}


static RC draw_fem_number_solid(GTKU_OBJECT *obj)
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
		FEM_ELEMENT elem = obj->elem->array[ii1];
		TRANSLATION3D center = init_translation3d();
		if(elem.label < 0) continue;
		for(ii2=0; ii2<elem.node_num; ii2++){
			int node_index;
			node_index=search_fem_node_label(*(obj)->node, elem.node[ii2]);
			RC_NEG_CHK(node_index);
			center.x += obj->node->array[node_index].p.x; 
			center.y += obj->node->array[node_index].p.y; 
			center.z += obj->node->array[node_index].p.z; 
		}
		center.x /= (double)elem.node_num;
		center.y /= (double)elem.node_num;
		center.z /= (double)elem.node_num;
		snprintf(buf, sizeof(buf), "%d", elem.label);
		buf[31] = '\0'; /* 念のため */
		gtku_gl_draw_string(center, buf);
	}
	if(obj->d_prop->text_trans_sw) unset_gl_translucent();

	return(NORMAL_RC);
}


static RC draw_fem_number_polygon(GTKU_OBJECT *obj)
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
		FEM_ELEMENT surf = obj->surf->array[ii1];
		TRANSLATION3D center = init_translation3d();
		if(surf.label < 0) continue;
		for(ii2=0; ii2<surf.node_num; ii2++){
			int node_index;
			node_index=search_fem_node_label(*(obj)->node, surf.node[ii2]);
			RC_NEG_CHK(node_index);
			center.x += obj->node->array[node_index].p.x; 
			center.y += obj->node->array[node_index].p.y; 
			center.z += obj->node->array[node_index].p.z; 
		}
		center.x /= (double)surf.node_num;
		center.y /= (double)surf.node_num;
		center.z /= (double)surf.node_num;
		snprintf(buf, sizeof(buf), "%d", surf.label);
		buf[31] = '\0'; /* 念のため */
		gtku_gl_draw_string(center, buf);
	}
	if(obj->d_prop->text_trans_sw) set_gl_translucent();

	return(NORMAL_RC);
}


static RC draw_image_object(GTKU_OBJECT *obj)
{
	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	RC_NULL_CHK(obj);
	GTKU_DRAW_TYPE draw_type = obj->d_prop->draw_type;

	if(draw_type & GTKU_DRAW_HIDE) return(NORMAL_RC);

	if(obj->type & GTKU_OBJECT_IMAGE_PIXEL){
		if(draw_type & GTKU_DRAW_DEFAULT) draw_type |= GTKU_DRAW_TEXTURE;
		if(draw_type & GTKU_DRAW_TEXTURE)
			RC_TRY( draw_pixel_texture(obj->gr, obj->d_prop) );
		if(draw_type & GTKU_DRAW_CELL)
			RC_TRY( draw_pixel_cell(obj->gr, obj->d_prop) );
	}else if(obj->type & GTKU_OBJECT_IMAGE_VOXEL){
		if( draw_type & (GTKU_DRAW_DEFAULT | GTKU_DRAW_CELL ) )
			RC_TRY( draw_voxel_cell(obj->gr, obj->d_prop) );
	}else{
		GTKU_WARNING("draw_image_object()", "Unknown object type");
	}

	return(NORMAL_RC);
}


/*
 * 2Dのピクセルイメージに関しては，z=0.0のX-Y平面に表示し，必要とあらば
 * OpenGL側で回転・移動をする．つまり実際の表示とは異なる場所にイメージが
 * あることになる．もし，積み重ねた輪切り写真を表示し，FEMモデルと重ね合わ
 * せた解析がしたいのならば専用の表示関数を作ること
 */
static RC draw_pixel_texture(GRAPHIC *gr, GTKU_DRAW_PROP *d_prop)
{
	double tcx, tcy; /* Texture Coord */
	TRANSLATION3D p;

	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	/* Texture Initialize */
	if(!(d_prop->init_texture)){
		RC_TRY( init_texture(gr, d_prop) );
		d_prop->init_texture = TRUE;
	}

	/* Initialize  Image Direction */
	if(!(d_prop->init_image_dir)){
		RC_TRY( init_image_direction(gr, d_prop) );
		d_prop->init_image_dir = TRUE;
		gtk_widget_queue_draw(gtku_gl->glarea);
		return(NORMAL_RC);
	}

	if(glIsTexture(d_prop->texname) == GL_FALSE){
		GTKU_WARNING("draw_pixel_texture()", "Texture has not setting");
		return(NORMAL_RC);
	}

	glEnable(GL_TEXTURE_2D);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
	glBindTexture(GL_TEXTURE_2D, d_prop->texname);
	glDisable(GL_CULL_FACE);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glBegin(GL_QUADS);

	tcx = gr->size.x/(double)d_prop->pack_size.x;
	tcy = gr->size.y/(double)d_prop->pack_size.y;
	p.x = (double)gr->size.x;
	p.y = (double)gr->size.y;

	glTexCoord2d(0.0, 0.0); glVertex2d(0.0, p.y);
	glTexCoord2d(0.0, tcy); glVertex2d(0.0, 0.0);
	glTexCoord2d(tcx, tcy); glVertex2d(p.x, 0.0);
	glTexCoord2d(tcx, 0.0); glVertex2d(p.x, p.y);

	/*
	if(d_prop->tex_dir == GTKU_TEXTURE_XY){ 
		glTexCoord2d(0.0, 0.0); glVertex3d(0.0, p.y, p.z);
		glTexCoord2d(0.0, tcy); glVertex3d(0.0, 0.0, p.z);
		glTexCoord2d(tcx, tcy); glVertex3d(p.x, 0.0, p.z);
		glTexCoord2d(tcx, 0.0); glVertex3d(p.x, p.y, p.z);
	}else if(d_prop->tex_dir == GTKU_TEXTURE_YZ){
		glTexCoord2d(0.0, 0.0); glVertex3d(p.z, p.y, p.x);
		glTexCoord2d(0.0, tcy); glVertex3d(p.z, 0.0, p.x);
		glTexCoord2d(tcx, tcy); glVertex3d(p.z, 0.0, 0.0);
		glTexCoord2d(tcx, 0.0); glVertex3d(p.z, p.y, 0.0);
	}else if(d_prop->tex_dir == GTKU_TEXTURE_ZX){
		glTexCoord2d(0.0, 0.0); glVertex3d(0.0, p.z, 0.0);
		glTexCoord2d(0.0, tcy); glVertex3d(0.0, p.z, p.y);
		glTexCoord2d(tcx, tcy); glVertex3d(p.x, p.z, p.y);
		glTexCoord2d(tcx, 0.0); glVertex3d(p.x, p.z, 0.0);
	}
	*/

	glEnd();
	glDisable(GL_TEXTURE_2D);

	return(NORMAL_RC);
}


/*
 * ここで行っている(ビット演算など)処理の意味は，
 * 「OpenGLプログラミングガイド」(9章)デクスチャマッピングの章を参照．
 * テクスチャーは必ず2^nの画像サイズでなければいけない．
 */
static RC init_texture(GRAPHIC *gr, GTKU_DRAW_PROP *d_prop)
{
	int ii1;
	GLubyte *image = NULL;
	unsigned long image_size;
	int x_size_bit = 0;
	int y_size_bit = 0;
	GLint max_tex_size;

	glGetIntegerv(GL_MAX_TEXTURE_SIZE, &max_tex_size);
	if( (gr->size.x > (unsigned long)max_tex_size) ||
	    (gr->size.y > (unsigned long)max_tex_size) ){
		GTKU_WARNING("init_texture()", "Too large Texture!!");
		g_warning("max texture size is %d on your environment\n", max_tex_size);
		return(NORMAL_RC);
	}

	for(ii1=0; ii1<sizeof(unsigned long)*8; ii1++){ /* 8 bit */ 
		if(gr->size.x >> ii1) x_size_bit = ii1; /* X size = 2^? */
		if(gr->size.y >> ii1) y_size_bit = ii1; /* Y size = 2^? */
	}

	if( (1<<x_size_bit) == (gr->size.x) ) d_prop->pack_size.x = 1<<x_size_bit;
	else d_prop->pack_size.x = 1<<(x_size_bit+1);

	if( (1<<y_size_bit) == (gr->size.y) ) d_prop->pack_size.y = 1<<y_size_bit;
	else d_prop->pack_size.y = 1<<(y_size_bit+1);

	d_prop->pack_size.z = 1;

	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glGenTextures(1, &d_prop->texname);
	glBindTexture(GL_TEXTURE_2D, d_prop->texname);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

	/* Fix!! mono16も読めるようにしよう */
	switch(gr->type){
	case GRAPHIC_RGB24:
		image_size = (d_prop->pack_size.x)*(d_prop->pack_size.y);
		RC_TRY( rc_malloc(sizeof(GLubyte)*image_size*3, (void *)(&image)) );
		trans_rgb24_GLimage(gr, d_prop, &image);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB,
		             d_prop->pack_size.x, d_prop->pack_size.y, 0,
		             GL_RGB, GL_UNSIGNED_BYTE, image); 
		break;
	default:
		 break;
	}

	return(NORMAL_RC);
}


static void trans_rgb24_GLimage(GRAPHIC *gr,  GTKU_DRAW_PROP *d_prop,
                                GLubyte **image)
{
	int ii1, ii2, ii3;

	ii3 = 0;
	for(ii1=0; ii1<(gr->size.y); ii1++){
		for(ii2=0; ii2<(gr->size.x); ii2++){
			(*image)[ii3++] = (GLubyte)gr->data.rgb24[0][ii1][ii2].r;
			(*image)[ii3++] = (GLubyte)gr->data.rgb24[0][ii1][ii2].g;
			(*image)[ii3++] = (GLubyte)gr->data.rgb24[0][ii1][ii2].b;
		}
		for(ii2=0; ii2<(d_prop->pack_size.x - gr->size.x); ii2++){
			(*image)[ii3++] = (GLubyte)0; /* dummy  padding */
			(*image)[ii3++] = (GLubyte)0; /* dummy  padding */
			(*image)[ii3++] = (GLubyte)0; /* dummy  padding */
		}
	}
	for(ii1=0; ii1<(d_prop->pack_size.y - gr->size.y); ii1++){
		for(ii2=0; ii2<(d_prop->pack_size.x); ii2++){
			(*image)[ii3++] = (GLubyte)0; /* dummy  padding */
			(*image)[ii3++] = (GLubyte)0; /* dummy  padding */
			(*image)[ii3++] = (GLubyte)0; /* dummy  padding */
		}
	}
}


static RC init_image_direction(GRAPHIC *gr, GTKU_DRAW_PROP *d_prop)
{
	RC_NULL_CHK(gr);
	RC_NULL_CHK(d_prop);

	init_operation_matrix(d_prop->shift_mat);
	init_operation_matrix(d_prop->rot_mat);

	/* Image Direction Initialize */
	if(d_prop->image_dir == GTKU_IMAGE_XY){ 
		shift_operation_matrix(0.0, 0.0, d_prop->image_depth,d_prop->shift_mat);
	}else if(d_prop->image_dir == GTKU_IMAGE_YZ){
		rotation_operation_matrix_y(PI/2, d_prop->rot_mat);
		shift_operation_matrix(d_prop->image_depth, 0.0, gr->size.x,
		                       d_prop->shift_mat);
	/* CT-imageだと逆かもしれない */
	}else if(d_prop->image_dir == GTKU_IMAGE_ZX){
		rotation_operation_matrix_y(-PI, d_prop->rot_mat);
		rotation_operation_matrix_x(PI/2, d_prop->rot_mat);
		shift_operation_matrix(gr->size.x, d_prop->image_depth, 0.0, 
		                       d_prop->shift_mat);
	}else{
		GTKU_WARNING("init_image_direction()", "Unknown Image Direction");
	}

	return(NORMAL_RC);
}


static RC draw_pixel_cell(GRAPHIC *gr, GTKU_DRAW_PROP *d_prop)
{
	/* Initialize  Image Direction */
	if(!(d_prop->init_image_dir)){
		RC_TRY( init_image_direction(gr, d_prop) );
		d_prop->init_image_dir = TRUE;
		gtk_widget_queue_draw(gtku_gl->glarea);
		return(NORMAL_RC);
	}

	/* Fix!! mono16なども読めるようにしよう */
	GRAPHIC_SIZE_LOOP_2D(gr->size.y, gr->size.x){
		switch(gr->type){
		case GRAPHIC_RGB24:
			glColor3d( gr->data.rgb24[0][ii1][ii2].r/(GLdouble)0xff,
			           gr->data.rgb24[0][ii1][ii2].g/(GLdouble)0xff,
			           gr->data.rgb24[0][ii1][ii2].b/(GLdouble)0xff );
			break;
		default:
			GTKU_WARNING("draw_pixel_cell()", "Unknown image type.");
		}
		draw_pixel_cell_single(gr->size.y, ii2, ii1);
	}GRAPHIC_SIZE_LOOP_2D_END;

	return(NORMAL_RC);
}


static void draw_pixel_cell_single(int height, int x, int y)
{
	int Y = height - y - 1;
	glBegin(GL_QUADS);
	glVertex2i( (x    ), (Y    ) );
	glVertex2i( (x + 1), (Y    ) );
	glVertex2i( (x + 1), (Y + 1) );
	glVertex2i( (x    ), (Y + 1) );
	glEnd();
}


static RC draw_voxel_cell(GRAPHIC *gr, GTKU_DRAW_PROP *d_prop)
{
	double value, min, max, avg;
	GRAPHIC_RGB color;

	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	/* Initialize  Image Direction */
	if(!(d_prop->init_image_dir)){
		RC_TRY( init_image_direction(gr, d_prop) );
		d_prop->init_image_dir = TRUE;
		gtk_widget_queue_draw(gtku_gl->glarea);
		return(NORMAL_RC);
	}

	if(gtku_gl->contour == NULL){
		gtku_gl->contour = (GRAPHIC_CONTOUR *)malloc(sizeof(GRAPHIC_CONTOUR));
		RC_NULL_CHK(gtku_gl->contour);
		RC_TRY( make_graphic_contour(gtku_gl->contour) );
	}

	if(gr->type == GRAPHIC_SCALAR)
		RC_TRY( chk_graphic_range_scalar(*gr, &min, &max, &avg) );

	GRAPHIC_SIZE_LOOP_3D(gr->size.z, gr->size.y, gr->size.x){
		switch(gr->type){
		case GRAPHIC_RGB24:
			color = gr->data.rgb24[ii1][ii2][ii3];
			glColor3d( color.r/(GLdouble)0xff,
			           color.g/(GLdouble)0xff,
			           color.b/(GLdouble)0xff );
			break;
		case GRAPHIC_SCALAR:
			value = gr->data.scalar[ii1][ii2][ii3];
			if(value < 128) continue; /* Fix!! GUIで操作するように */
			RC_TRY( get_graphic_contour_color24(value, min, max,
			               *(gtku_gl->contour), d_prop->cont_type, &color) );
			glColor3d( color.r/(GLdouble)0xff,
			           color.g/(GLdouble)0xff,
			           color.b/(GLdouble)0xff );
			break;
		default:
			GTKU_WARNING("draw_voxel_cell()", "Unknown image type.");
		}
		draw_voxel_cell_single(ii3, ii2, ii1);
	}GRAPHIC_SIZE_LOOP_3D_END;

	return(NORMAL_RC);
}


/* See element.eps, set_fem_face_table_hexa1() */
static void draw_voxel_cell_single(int x, int y, int z)
{
	glBegin(GL_QUADS); /* 0-3-2-1 */
	glVertex3i( (x    ), (y    ), (z    ) ); 
	glVertex3i( (x    ), (y + 1), (z    ) );
	glVertex3i( (x + 1), (y + 1), (z    ) );
	glVertex3i( (x + 1), (y    ), (z    ) );
	glEnd();
	glBegin(GL_QUADS); /* 0-1-5-4 */
	glVertex3i( (x    ), (y    ), (z    ) ); 
	glVertex3i( (x + 1), (y    ), (z    ) );
	glVertex3i( (x + 1), (y    ), (z + 1) );
	glVertex3i( (x    ), (y    ), (z + 1) );
	glEnd();
	glBegin(GL_QUADS); /* 1-2-6-5 */
	glVertex3i( (x + 1), (y    ), (z    ) ); 
	glVertex3i( (x + 1), (y + 1), (z    ) );
	glVertex3i( (x + 1), (y + 1), (z + 1) );
	glVertex3i( (x + 1), (y    ), (z + 1) );
	glEnd();
	glBegin(GL_QUADS); /* 7-4-5-6 */
	glVertex3i( (x    ), (y + 1), (z + 1) ); 
	glVertex3i( (x    ), (y    ), (z + 1) );
	glVertex3i( (x + 1), (y    ), (z + 1) );
	glVertex3i( (x + 1), (y + 1), (z + 1) );
	glEnd();
	glBegin(GL_QUADS); /* 3-0-4-7 */
	glVertex3i( (x    ), (y + 1), (z    ) ); 
	glVertex3i( (x    ), (y    ), (z    ) );
	glVertex3i( (x    ), (y    ), (z + 1) );
	glVertex3i( (x    ), (y + 1), (z + 1) );
	glEnd();
	glBegin(GL_QUADS); /* 2-3-7-6 */
	glVertex3i( (x + 1), (y + 1), (z    ) ); 
	glVertex3i( (x    ), (y + 1), (z    ) );
	glVertex3i( (x    ), (y + 1), (z + 1) );
	glVertex3i( (x + 1), (y + 1), (z + 1) );
	glEnd();
}


/* The "configure_event" signal handler. Almost always it will be used */
/* to resize the OpenGL viewport when the window is resized.           */
gboolean gtku_gl_reshape(GtkWidget *widget, GdkEventConfigure *event)
{
	GdkGLContext *glcontext = gtk_widget_get_gl_context(widget);
	GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable(widget);

	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

#ifdef GTKU_VERBOSE_CHECK
	/* OpenGL functions can be called only if make_current returns true */
	if(!gdk_gl_drawable_make_current(gldrawable, glcontext)){
		g_warning("glarea not drawable\n");
		return(FALSE);
	}
#endif

	/*** OpenGL BEGIN ***/
	if(gdk_gl_drawable_gl_begin(gldrawable, glcontext)){
		glViewport(0,0, widget->allocation.width,
		                widget->allocation.height);
		gtku_gl->view_size.x = widget->allocation.width/2;
		gtku_gl->view_size.y = widget->allocation.height/2;
		gdk_gl_drawable_gl_end(gldrawable);
	}/*** OpenGL END ***/

	return(TRUE);
}


/*
 * The "realize" signal handler.  Usually, it will be called one time.
 * All the OpenGL initialization should be performed here normally. But,
 * this function only calculate model(s) size and model center position
 * and so on.
 */
gboolean gtku_gl_realize(GtkWidget *widget)
{
	PangoFontDescription *font_desc;
	PangoFont *font;

	GdkGLContext  *glcontext  = gtk_widget_get_gl_context(widget);
	GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable(widget);

	/* OpenGL functions can be called only if make_current returns true */
	if(!gdk_gl_drawable_make_current(gldrawable, glcontext)){
		GTKU_WARNING("gtku_gl_realize()", "glarea not drawable");
		return FALSE; 
	}

	/* calculate object view size and so on.*/
	GTKU_BOOLEAN_RC_TRY( init_objects_drawing_attribute() );
	GTKU_BOOLEAN_RC_TRY( init_objects_operation_matrix() );
	GTKU_BOOLEAN_RC_TRY( init_view_volume() );

	/* Allocate ASCII character bitmap list */
	if(gtku_gl->char_list == 0){
		gtku_gl->char_list = glGenLists(128);
		g_return_val_if_fail(gtku_gl->char_list, FALSE);

		font_desc = pango_font_description_from_string(font_name);
		font = gdk_gl_font_use_pango_font(font_desc,0,128,gtku_gl->char_list);
		if(!font){
			GTKU_WARNING("gtku_gl_realize()", "Can't load font");
			return FALSE; 
		}
		pango_font_description_free(font_desc);
	}

	return(TRUE);
}


static RC init_objects_drawing_attribute(void)
{
	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	int ii1;
	for(ii1=0; ii1<gtku_gl->object->size; ii1++){
		gtku_gl->object->obj[ii1]->d_prop->init_image_dir = FALSE;
	}

	return(NORMAL_RC);
}


static RC init_objects_operation_matrix(void)
{
	int ii1;

	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	/* initialize global operation matrix */
	init_operation_matrix(gtku_gl->shift_mat);
	init_operation_matrix(gtku_gl->scale_mat);
	init_operation_matrix(gtku_gl->rot_mat);
	init_operation_matrix(gtku_gl->initial_mat);

	/* initialize local operation matrix */
	for(ii1=0; ii1<gtku_gl->object->size; ii1++){
		init_operation_matrix(gtku_gl->object->obj[ii1]->d_prop->shift_mat);
		init_operation_matrix(gtku_gl->object->obj[ii1]->d_prop->scale_mat);
		init_operation_matrix(gtku_gl->object->obj[ii1]->d_prop->rot_mat);
	}

	return(NORMAL_RC);
}


static RC init_view_volume(void)
{
	double rate;

	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	/* Objectの中心を求め、描画範囲に移動させる */
	cal_objects_size(gtku_gl->object, &(gtku_gl->object->center),
	                 &(gtku_gl->object->max), &(gtku_gl->object->min));
	shift_operation_matrix(-(gtku_gl->object->center.x),
	                       -(gtku_gl->object->center.y),
	                       -(gtku_gl->object->center.z), gtku_gl->initial_mat);

	/* Objectを描画範囲に合わせ、拡大縮小させる */
	RC_TRY( cal_objects_view_scale(&rate) );
	scale_operation_matrix(rate, rate, rate, gtku_gl->initial_mat);

	return(NORMAL_RC);
}


static RC cal_objects_size(GTKU_OBJECT_ARRAY *object, TRANSLATION3D *center,
                           TRANSLATION3D *max, TRANSLATION3D *min)
{
	int ii1;

	for(ii1=0; ii1<object->size; ii1++){
		GTKU_OBJECT *obj = object->obj[ii1];
		RC_NULL_CHK(obj);
		if(obj->type == GTKU_OBJECT_NONE) continue;
		else if(obj->type & GTKU_OBJECT_FEM_NODE){
			RC_TRY( cal_object_size_fem(obj->node, &obj->d_prop->center,
			                            &obj->d_prop->max, &obj->d_prop->min) );
		}else if( (obj->type & GTKU_OBJECT_IMAGE_PIXEL) ||
		          (obj->type & GTKU_OBJECT_IMAGE_VOXEL) ){
			RC_TRY( cal_object_size_graphic(obj->gr, &obj->d_prop->center,
			                                &obj->d_prop->max,
			                                &obj->d_prop->min) );
		}else{
			GTKU_WARNING("cal_objects_size()", "Unknown Object");
		}
	}

	/* object  check. if object is NONE? what do you doing!! fix */
	*max = object->obj[0]->d_prop->max;
	*min = object->obj[0]->d_prop->min;
	for(ii1=0; ii1<object->size; ii1++){
		max->x = MAX2(object->obj[ii1]->d_prop->max.x, max->x);
		max->y = MAX2(object->obj[ii1]->d_prop->max.y, max->y);
		max->z = MAX2(object->obj[ii1]->d_prop->max.z, max->z);
		min->x = MIN2(object->obj[ii1]->d_prop->min.x, min->x);
		min->y = MIN2(object->obj[ii1]->d_prop->min.y, min->y);
		min->z = MIN2(object->obj[ii1]->d_prop->min.z, min->z);
	}
	*center = mid_point3d(*max, *min);

	return(NORMAL_RC);
}


static RC cal_object_size_fem(FEM_NODE_ARRAY *node, TRANSLATION3D *center,
                              TRANSLATION3D *max, TRANSLATION3D *min)
{
	int ii1;
	*max = node->array[0].p;
	*min = node->array[0].p;
	for(ii1=1; ii1<node->size; ii1++){
		max->x = MAX2(node->array[ii1].p.x, max->x);
		max->y = MAX2(node->array[ii1].p.y, max->y);
		max->z = MAX2(node->array[ii1].p.z, max->z);
		min->x = MIN2(node->array[ii1].p.x, min->x);
		min->y = MIN2(node->array[ii1].p.y, min->y);
		min->z = MIN2(node->array[ii1].p.z, min->z);
	}
	*center = mid_point3d(*max, *min);

	return(NORMAL_RC);
}


static RC cal_object_size_graphic(GRAPHIC *gr, TRANSLATION3D *center,
                                  TRANSLATION3D *max, TRANSLATION3D *min)
{
	*min = init_translation3d();
	max->x = (double)gr->size.x;
	max->y = (double)gr->size.y;
	max->z = (double)gr->size.z;
	*center = mid_point3d(*max, *min);

	return(NORMAL_RC);
}


static RC cal_objects_view_scale(double *scale)
{
	double length;
	TRANSLATION3D view_box;

	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	view_box.x = gtku_gl->object->max.x - gtku_gl->object->min.x;
	view_box.y = gtku_gl->object->max.y - gtku_gl->object->min.y;
	view_box.z = gtku_gl->object->max.z - gtku_gl->object->min.z;
	length = abs_translation3d(view_box);
	if(gtku_gl->glarea_w_size < gtku_gl->glarea_h_size){
		*scale = gtku_gl->glarea_w_size/length;
	}else{
		*scale = gtku_gl->glarea_h_size/length;
	}

	return(NORMAL_RC);
}


/* The "motion_notify_event" signal handler. It will be called while  */
/* mouse cursor moving motion in widget area.                         */
gboolean gtku_gl_motion_notify(GtkWidget *widget, GdkEventMotion *event)
{
	GdkRectangle area;
	GdkModifierType state;
	gboolean redraw = FALSE;
	int x, y;

	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	if(gtku_gl->animation)
		return(TRUE);

	if(event->is_hint){
		gdk_window_get_pointer(event->window, &x, &y, &state);
	}else{
		x = event->x;
		y = event->y;
		state = event->state;
	}

	gtku_gl->d.x = (double)x - gtku_gl->pre.x;
	if(fabs(gtku_gl->d.x) < MOVE_THRESHOLD) gtku_gl->d.x = 0.0;

	gtku_gl->d.y = gtku_gl->pre.y - (double)y; /* GDKでは画面左上が(0,0) */
	if(fabs(gtku_gl->d.y) < MOVE_THRESHOLD) gtku_gl->d.y = 0.0;

	area.width  = widget->allocation.width;
	area.height = widget->allocation.height;

	if(state & GDK_BUTTON2_MASK){
		/* Middle Button [Rotation] */
		GTKU_BOOLEAN_RC_TRY(glarea_mouse_rotation(widget, area));
		redraw = TRUE;
	}else if(state & GDK_BUTTON3_MASK){
		/* Left Button [Translation] */
		GTKU_BOOLEAN_RC_TRY(glarea_mouse_shift(widget, area));
		redraw = TRUE;
	}else if(state & (GDK_BUTTON4_MASK|GDK_BUTTON5_MASK)){
		/* Scroll UP and Down [Scale] */
		GTKU_BOOLEAN_RC_TRY(glarea_mouse_scale(widget, area, state ));
		redraw = TRUE;
	}

	if(redraw && !gtku_gl->animation) gtk_widget_queue_draw(widget);

	gtku_gl->pre.x = (double)x;
	gtku_gl->pre.y = (double)y;

	return(TRUE);
}


#define ROT_ANG (PI/6)
/* motion_notify event sub call back */
static RC glarea_mouse_rotation(GtkWidget *widget, GdkRectangle area)
{
	int ii1;
	double **rot;
	double r, angle, angle1, angle2;


	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	if( (gtku_gl->d.x == 0.0) && (gtku_gl->d.y == 0.0) )
		return(NORMAL_RC);

	RC_TRY( allocate_operation_matrix(&rot) );

	/* Key Press Effect */
	r = abs_translation2d(gtku_gl->d);
	if(gtku_gl->key.shift == TRUE) angle = ROT_ANG;
	else                           angle = PI*r/area.height;

	if(gtku_gl->key.x == TRUE){ /* X */
		if( fabs(gtku_gl->d.x) > fabs(gtku_gl->d.y) ){
			if(gtku_gl->d.x > 0) rotation_operation_matrix_x( angle, rot);
			else                 rotation_operation_matrix_x(-angle, rot);
		}else{
			if(gtku_gl->d.y > 0) rotation_operation_matrix_x( angle, rot);
			else                 rotation_operation_matrix_x(-angle, rot);
		}
	}else if(gtku_gl->key.y == TRUE){ /* Y */
		if( fabs(gtku_gl->d.x) > fabs(gtku_gl->d.y) ){
			if(gtku_gl->d.x > 0) rotation_operation_matrix_y( angle, rot);
			else                 rotation_operation_matrix_y(-angle, rot);
		}else{
			if(gtku_gl->d.y > 0) rotation_operation_matrix_y( angle, rot);
			else                 rotation_operation_matrix_y(-angle, rot);
		}
	}else if(gtku_gl->key.z == TRUE){ /* Z */
		if( fabs(gtku_gl->d.x) > fabs(gtku_gl->d.y) ){
			if(gtku_gl->d.x > 0) rotation_operation_matrix_z( angle, rot);
			else                 rotation_operation_matrix_z(-angle, rot);
		}else{
			if(gtku_gl->d.y > 0) rotation_operation_matrix_z( angle, rot);
			else                 rotation_operation_matrix_z(-angle, rot);
		}
	}else if(gtku_gl->key.shift == TRUE){ /* SHIFT */
		if( fabs(gtku_gl->d.x) > 2*fabs(gtku_gl->d.y) ){
			if(gtku_gl->d.x > 0) rotation_operation_matrix_y( angle, rot);
			else                 rotation_operation_matrix_y(-angle, rot);
		}else{ 
			if(gtku_gl->d.y > 0) rotation_operation_matrix_x( angle, rot);
			else                 rotation_operation_matrix_x(-angle, rot);
		}
	}else{
		angle1 = acos(-gtku_gl->d.y/r);
		angle2 = r/area.height*5;
		if(gtku_gl->d.x < 0) angle1 *= (-1);
		rotation_operation_matrix_z(-angle1, rot);
		rotation_operation_matrix_x( angle2, rot);
		rotation_operation_matrix_z( angle1, rot);
	}

	if(gtku_gl->rot_lock_count == 0){ /* global rotation */
		RC_TRY( mul_operation_matrix(rot, gtku_gl->rot_mat) );
		/* グローバル座標で回転する
		if( (gtku_gl->key.x == TRUE) || (gtku_gl->key.y == TRUE) ||
		    (gtku_gl->key.z == TRUE) ){
			RC_TRY( mul_operation_matrix(gtku_gl->rot_mat, rot) );
			RC_TRY( copy_matrix(4, 4, rot, gtku_gl->rot_mat) );
		}
		*/
	}else{ /* local rotation */
		for(ii1=0; ii1<gtku_gl->object->size; ii1++){
			GTKU_DRAW_PROP *d_prop = gtku_gl->object->obj[ii1]->d_prop;
			if(d_prop->rot_lock_sw == TRUE) continue;
			RC_TRY( mouse_rotation_local(d_prop, rot) );
		}
	}

	free_operation_matrix(rot);

	return(NORMAL_RC);
}


static RC mouse_rotation_local(GTKU_DRAW_PROP *d_prop, double **rot)
{
	double **center, **center_inv;

	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	RC_TRY( allocate_operation_matrix(&center) );
	RC_TRY( allocate_operation_matrix(&center_inv) );
	shift_operation_matrix( -d_prop->center.x, -d_prop->center.y,
	                        -d_prop->center.z, center);
	RC_TRY( inverse_matrix(4, center, center_inv) );
	RC_TRY( mul_operation_matrix(center,           d_prop->rot_mat) );

	if(d_prop->use_local_axis_sw){
		RC_TRY( mul_operation_matrix(rot, d_prop->rot_mat) );
	}else{
		double **rot_inv;
		RC_TRY( allocate_operation_matrix(&rot_inv) );
		RC_TRY( inverse_matrix(4, gtku_gl->rot_mat, rot_inv) );
		RC_TRY( mul_operation_matrix(gtku_gl->rot_mat, d_prop->rot_mat) );
		RC_TRY( mul_operation_matrix(rot,              d_prop->rot_mat) );
		RC_TRY( mul_operation_matrix(rot_inv,          d_prop->rot_mat) );
		free_operation_matrix(rot_inv);
	}

	RC_TRY( mul_operation_matrix(center_inv, d_prop->rot_mat) );
	free_operation_matrix(center);
	free_operation_matrix(center_inv);

	return(NORMAL_RC);
}


/* motion_notify event sub call back    */
/* for 5 button mouse. anyone need this */
static RC glarea_mouse_scale(GtkWidget *widget, GdkRectangle area,
                             GdkModifierType state)
{
	int ii1;
	TRANSLATION3D rate;
	double d_rate;

	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	if( fabs(gtku_gl->d.x) > fabs(gtku_gl->d.y) ){
		d_rate = abs_translation2d(gtku_gl->d)/area.width;
		if(gtku_gl->d.x < 0) d_rate *= -1;
	}else{
		d_rate = abs_translation2d(gtku_gl->d)/area.height;
		if(gtku_gl->d.y < 0) d_rate *= -1;
	}

	/* Initialize */
	rate.x = rate.y = rate.z = 1.0;

	/* Key Press Effect */
	     if(gtku_gl->key.x == TRUE) rate.x += d_rate;
	else if(gtku_gl->key.y == TRUE) rate.y += d_rate;
	else if(gtku_gl->key.z == TRUE) rate.z += d_rate;
	else{ 
		rate.x += d_rate;
		rate.y += d_rate;
		rate.z += d_rate;
	}

	if(gtku_gl->scale_lock_count == 0){ /* global scale */
		/* グローバル座標で伸縮する
		if( (gtku_gl->key.x == TRUE) || (gtku_gl->key.y == TRUE) ||
		    (gtku_gl->key.z == TRUE) ){
			double **rot_inv;
			RC_TRY( allocate_operation_matrix(&rot_inv) );
			RC_TRY( inverse_matrix(4, gtku_gl->rot_mat, rot_inv));
			RC_TRY( mul_operation_matrix(rot_inv, gtku_gl->scale_mat) );
			scale_operation_matrix(rate.x, rate.y, rate.z, gtku_gl->scale_mat);
			RC_TRY( mul_operation_matrix(gtku_gl->rot_mat, gtku_gl->scale_mat));
			free_operation_matrix(rot_inv);
		}
		*/
		scale_operation_matrix(rate.x, rate.y, rate.z, gtku_gl->scale_mat);
	}else{ /* local scale */
		for(ii1=0; ii1<gtku_gl->object->size; ii1++){
			GTKU_DRAW_PROP *d_prop = gtku_gl->object->obj[ii1]->d_prop;
			if(d_prop->scale_lock_sw == TRUE) continue;
			RC_TRY( mouse_scale_local(d_prop, rate) );
		}
	}

	return(NORMAL_RC);
}


static RC mouse_scale_local(GTKU_DRAW_PROP *d_prop, TRANSLATION3D rate)
{
	double **ini_inv;

	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	RC_TRY( allocate_operation_matrix(&ini_inv) );
	RC_TRY( inverse_matrix(4, gtku_gl->initial_mat, ini_inv));
	RC_TRY( mul_operation_matrix(gtku_gl->initial_mat, d_prop->scale_mat) );

	if(d_prop->use_local_axis_sw){
		scale_operation_matrix(rate.x, rate.y, rate.z, d_prop->scale_mat);
	}else{
		double **rot_inv;
		RC_TRY( allocate_operation_matrix(&rot_inv) );
		RC_TRY( inverse_matrix(4, gtku_gl->rot_mat, rot_inv));
		RC_TRY( mul_operation_matrix(gtku_gl->rot_mat, d_prop->scale_mat) );
		scale_operation_matrix(rate.x, rate.y, rate.z, d_prop->scale_mat);
		RC_TRY( mul_operation_matrix(rot_inv, d_prop->scale_mat) );
		free_operation_matrix(rot_inv);
	}

	RC_TRY( mul_operation_matrix(ini_inv, d_prop->scale_mat) );
	free_operation_matrix(ini_inv);

	return(NORMAL_RC);
}




/* motion_notify event sub call back */
static RC glarea_mouse_shift(GtkWidget *widget, GdkRectangle area)
{
	TRANSLATION3D d;

	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	/* Initialize */
	d = init_translation3d(); 

	/* Key Press Effect */
	     if(gtku_gl->key.x == TRUE) d.x = gtku_gl->d.x;
	else if(gtku_gl->key.y == TRUE) d.y = gtku_gl->d.y;
	else if(gtku_gl->key.z == TRUE){
		if( fabs(gtku_gl->d.x) > fabs(gtku_gl->d.y) ){
			d.z = abs_translation2d(gtku_gl->d);
			if(gtku_gl->d.x < 0) d.z *= -1;
		}else{
			d.z = abs_translation2d(gtku_gl->d);
			if(gtku_gl->d.y < 0) d.z *= -1;
		}
	}else if(gtku_gl->key.shift == TRUE){
		if( fabs(gtku_gl->d.x) > fabs(gtku_gl->d.y) ) d.x = gtku_gl->d.x;
		else d.y = gtku_gl->d.y;
	}else{ /* NO Key */
		d.x = gtku_gl->d.x;
		d.y = gtku_gl->d.y;
	}

	{ /* 拡大率によって移動量を変更する */
		double **scale_inv;
		RC_TRY( allocate_operation_matrix(&scale_inv) );
		RC_TRY( inverse_matrix(4, gtku_gl->scale_mat, scale_inv));
		d =  operation3d(scale_inv, d);
		free_operation_matrix(scale_inv);
	}

	if(gtku_gl->shift_lock_count == 0){ /* global shift */
		/* グローバル座標で移動する
		if( (gtku_gl->key.x == TRUE) || (gtku_gl->key.y == TRUE) ||
		    (gtku_gl->key.z == TRUE) ){
			double **rot_inv;
			RC_TRY( allocate_operation_matrix(&rot_inv) );
			RC_TRY( inverse_matrix(4, gtku_gl->rot_mat, rot_inv));
			RC_TRY( mul_operation_matrix(rot_inv, gtku_gl->shift_mat) );
			shift_operation_matrix(d.x, d.y, d.z, gtku_gl->shift_mat);
			RC_TRY( mul_operation_matrix(gtku_gl->rot_mat, gtku_gl->shift_mat));
			free_operation_matrix(rot_inv);
		}
		*/
		shift_operation_matrix(d.x, d.y, d.z, gtku_gl->shift_mat);
	}else{ /* local shift */      
		int ii1;
		for(ii1=0; ii1<gtku_gl->object->size; ii1++){
			GTKU_DRAW_PROP *d_prop = gtku_gl->object->obj[ii1]->d_prop;
			if(d_prop->shift_lock_sw == TRUE) continue;
			RC_TRY( mouse_shift_local(d_prop, d) );
		}
	}

	return(NORMAL_RC);
}


static RC mouse_shift_local(GTKU_DRAW_PROP *d_prop, TRANSLATION3D disp)
{
	double **ini_inv;

	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	RC_TRY( allocate_operation_matrix(&ini_inv) );
	RC_TRY( inverse_matrix(4, gtku_gl->initial_mat, ini_inv));
	RC_TRY( mul_operation_matrix(gtku_gl->initial_mat, d_prop->shift_mat) );

	if(d_prop->use_local_axis_sw){
		shift_operation_matrix(disp.x, disp.y, disp.z, d_prop->shift_mat);
	}else{
		double **rot_inv;
		RC_TRY( allocate_operation_matrix(&rot_inv) );
		RC_TRY( inverse_matrix(4, gtku_gl->rot_mat, rot_inv));
		RC_TRY( mul_operation_matrix(gtku_gl->rot_mat, d_prop->shift_mat) );
		shift_operation_matrix(disp.x, disp.y, disp.z, d_prop->shift_mat);
		RC_TRY( mul_operation_matrix(rot_inv, d_prop->shift_mat) );
		free_operation_matrix(rot_inv);
	}

	RC_TRY( mul_operation_matrix(ini_inv, d_prop->shift_mat) );
	free_operation_matrix(ini_inv);

	return(NORMAL_RC);
}


static gboolean add_rotation_matrix_for_animation(GtkWidget *widget)
{
	GdkRectangle area;

	area.width  = widget->allocation.width;
	area.height = widget->allocation.height;
	glarea_mouse_rotation(widget, area);

	return(TRUE);
}


/* The "button_press_event" signal handler. */
gboolean gtku_gl_button_press(GtkWidget *widget, GdkEventButton *event)
{
	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	if(gtku_gl->animation){ 
		if( (event->button == 1)||(event->button == 2)||(event->button == 3))
		 gtku_gl_toggle_animation(widget);
	}else{
		gtku_gl->pre.x = event->x;
		gtku_gl->pre.y = event->y;
	}

	return(TRUE);
}


/* The "button_press_release" signal handler. */
gboolean gtku_gl_button_release(GtkWidget *widget, GdkEventButton *event)
{
	double n;

	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	if( (event->button == 1)||(event->button == 2)||(event->button == 3) ){
		gtku_gl->cur.x = event->x;
		gtku_gl->cur.y = event->y;
	}

	if(gtku_gl->animation){ 
		return(TRUE);
	}else if(!gtku_gl->anim_lock_sw){
		n = gtku_gl->d.x*gtku_gl->d.x + gtku_gl->d.y*gtku_gl->d.y;
		if( (event->button == 2) && (n > ANIMATE_THRESHOLD) )
			gtku_gl_toggle_animation(widget);
	}

	return(TRUE);
}


/* The "scroll" signal handler. */
gboolean gtku_gl_scroll(GtkWidget *widget, GdkEventScroll *event)
{
	TRANSLATION3D rate;
	double d_rate = 0.0;

	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	if(event->direction == GDK_SCROLL_UP) d_rate = 0.005;
	if(event->direction == GDK_SCROLL_DOWN) d_rate = -0.005;
	
	/* Initialize */
	rate.x = rate.y = rate.z = 1.0;

	/* Key Press Effect */
	     if(gtku_gl->key.x == TRUE) rate.x += d_rate;
	else if(gtku_gl->key.y == TRUE) rate.y += d_rate;
	else if(gtku_gl->key.z == TRUE) rate.z += d_rate;
	else{ 
		rate.x += d_rate;
		rate.y += d_rate;
		rate.z += d_rate;
	}

	if(gtku_gl->scale_lock_count == 0){ /* global scale */
		scale_operation_matrix(rate.x, rate.y, rate.z, gtku_gl->scale_mat);
	}else{ /* local scale */
		int ii1;
		for(ii1=0; ii1<gtku_gl->object->size; ii1++){
			GTKU_DRAW_PROP *d_prop = gtku_gl->object->obj[ii1]->d_prop;
			if(d_prop->scale_lock_sw == FALSE)
				scale_operation_matrix(rate.x,rate.y,rate.z,d_prop->scale_mat);
		}
	}

	if(!gtku_gl->animation) gtk_widget_queue_draw(widget);

	return(TRUE);
}


/* The "map_event" signal handler. */
gboolean gtku_gl_map(GtkWidget *widget, GdkEvent *event)
{
	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);
	if(gtku_gl->animation) timeout_add(widget);
	return(TRUE);
}


/* The "unmap_event" signal handler. */
gboolean gtku_gl_unmap(GtkWidget *widget, GdkEvent *event)
{
	timeout_remove(widget);
	return(TRUE);
}


/* The "visibility_notify_event" signal handler. */
gboolean gtku_gl_visibility_notify(GtkWidget *widget, GdkEventVisibility *event)
{
	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);
	if(gtku_gl->animation){
		if(event->state == GDK_VISIBILITY_FULLY_OBSCURED)
			timeout_remove(widget);
	}else timeout_add(widget);
	return(TRUE);
}


gboolean gtku_gl_destroy(GtkWidget *widget)
{
	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);
	if(gtku_gl->animation){
		timeout_remove(widget);
		free(gtku_gl);
	}
	return(TRUE);
}


/* The "key_press_event" signal handler. */
gboolean gtku_gl_key_press(GtkWidget *widget, GdkEventKey *event)
{
	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	switch (event->keyval){
	case GDK_Shift_L:
	case GDK_Shift_R:
		gtku_gl->key.shift = TRUE; break;
	case GDK_Control_L:
	case GDK_Control_R:
		gtku_gl->key.control = TRUE; break;
	case GDK_Alt_L:
	case GDK_Alt_R:
		gtku_gl->key.alt = TRUE; break;
	case GDK_x:
	case GDK_X:
		gtku_gl->key.x = TRUE; break;
	case GDK_y:
	case GDK_Y:
		gtku_gl->key.y = TRUE; break;
	case GDK_z:
	case GDK_Z:
		gtku_gl->key.z = TRUE; break;
	case GDK_plus:
	case GDK_minus:
		glarea_key_scale(widget, event); break;
	default:
		return TRUE;
	}

	return TRUE;
}


/* The "key_release_event" signal handler. */
gboolean gtku_gl_key_release(GtkWidget *widget,GdkEventKey *event)
{
	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	switch (event->keyval){
	case GDK_Shift_L:
	case GDK_Shift_R:
		gtku_gl->key.shift = FALSE; break;
	case GDK_Control_L:
	case GDK_Control_R:
		gtku_gl->key.control = FALSE; break;
	case GDK_Alt_L:
	case GDK_Alt_R:
		gtku_gl->key.alt = FALSE; break;
	case GDK_x:
	case GDK_X:
		gtku_gl->key.x = FALSE; break;
	case GDK_y:
	case GDK_Y:
		gtku_gl->key.y = FALSE; break;
	case GDK_z:
	case GDK_Z:
		gtku_gl->key.z = FALSE; break;
	default:
		return TRUE;
	}

	return TRUE;
}


/* The "key_release_event" sub call back */
static gboolean glarea_key_scale(GtkWidget *widget, GdkEventKey *event)
{
	TRANSLATION3D rate;
	double d_rate = 0.0;

	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	if(event->hardware_keycode & GDK_plus) d_rate = 0.005;
	else if(event->hardware_keycode & GDK_minus) d_rate = -0.005;
	
	/* Initialize */
	rate.x = rate.y = rate.z = 1.0;

	/* Key Press Effect */
	     if(gtku_gl->key.x == TRUE) rate.x += d_rate;
	else if(gtku_gl->key.y == TRUE) rate.y += d_rate;
	else if(gtku_gl->key.z == TRUE) rate.z += d_rate;
	else{ 
		rate.x += d_rate;
		rate.y += d_rate;
		rate.z += d_rate;
	}

	if(gtku_gl->scale_lock_count == 0){ /* global scale */
		scale_operation_matrix(rate.x, rate.y, rate.z, gtku_gl->scale_mat);
	}else{ /* local scale */
		int ii1;
		for(ii1=0; ii1<gtku_gl->object->size; ii1++){
			GTKU_DRAW_PROP *d_prop = gtku_gl->object->obj[ii1]->d_prop;
			if(d_prop->scale_lock_sw == TRUE) continue;
			RC_TRY( mouse_scale_local(d_prop, rate) );
		}
	}

	if(!gtku_gl->animation) gtk_widget_queue_draw(widget);

	return(TRUE);
}


static gboolean timeout(GtkWidget *widget)
{
	gtk_widget_queue_draw(widget);
	return(TRUE);
}


static gboolean timeout_add(GtkWidget *widget)
{
	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);
	if(gtku_gl->timeout_id == 0){
		gtku_gl->timeout_id = gtk_timeout_add(TIMEOUT_INTERVAL,
		                                     (GtkFunction)timeout, widget);
	}
  	return(TRUE);
}


static gboolean timeout_remove(GtkWidget *widget)
{
	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);
	if(gtku_gl->timeout_id != 0){
		gtk_timeout_remove(gtku_gl->timeout_id);
		gtku_gl->timeout_id = 0;
    }
  	return(TRUE);
}


gboolean gtku_gl_toggle_animation(GtkWidget *widget)
{
	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);
	gtku_gl->animation = !(gtku_gl->animation);

	if(gtku_gl->animation){
		double r = abs_translation2d(gtku_gl->d)/ANIMETION_LIMIT_THRESHOLD;
		if( (abs(gtku_gl->d.x) > ANIMETION_LIMIT_THRESHOLD) ||
		    (abs(gtku_gl->d.y) > ANIMETION_LIMIT_THRESHOLD) ){
			gtku_gl->d.x /= r;
			gtku_gl->d.y /= r;
		}
		timeout_add(widget);
	}else{
		gtku_gl->d.x = 0;
		gtku_gl->d.y = 0;
		timeout_remove(widget);
		gtk_widget_queue_draw(widget);
	}
	return(TRUE);
}


void gtku_gl_draw_string(TRANSLATION3D p, const char *str)
{
	/* GTKU_GL existence check */
	g_return_if_fail(gtku_gl);

    glDisable(GL_LIGHTING);
	glRasterPos3d(p.x, p.y, p.z);

	glPushAttrib(GL_LIST_BIT);
	glListBase(gtku_gl->char_list);
	glCallLists(strlen(str), GL_UNSIGNED_BYTE, str);
	glPopAttrib();
}


/* "..." = "node", "element", "surface", "pixcel", "voxel" */
gboolean gtku_object_exist_chk(GTKU_OBJECT *obj, int size, ...)
{
	int ii1;
	char *name;
	va_list ap;
	gboolean truth = TRUE;

	if(obj == NULL){
		GTKU_WARNING("gtku_object_exist_chk()", "object is NULL!!");
		g_return_val_if_fail(obj, FALSE);
	}

	va_start(ap, size);
	for(ii1=0; ii1<size; ii1++){
		name = va_arg(ap, char *);
		if(strncmp("node", name, 4) == 0){
			if(!(obj->node)) truth = FALSE;
			if(!(obj->type & GTKU_OBJECT_FEM_NODE)) truth = FALSE;
		}else if(strncmp("element", name, 7) == 0){
			if(!(obj->elem)) truth = FALSE;
			if(!(obj->type & GTKU_OBJECT_FEM_ELEM)) truth = FALSE;
		}else if(strncmp("surface", name, 7) == 0){
			if(!(obj->surf)) truth = FALSE;
			if(!(obj->type & GTKU_OBJECT_FEM_SURF)) truth = FALSE;
		}else if(strncmp("pixel", name, 5) == 0){
			if(!(obj->gr)) truth = FALSE;
			if(!(obj->type & GTKU_OBJECT_IMAGE_PIXEL)) truth = FALSE;
		}else if(strncmp("voxel", name, 5) == 0){
			if(!(obj->gr)) truth = FALSE;
			if(!(obj->type & GTKU_OBJECT_IMAGE_VOXEL)) truth = FALSE;
		}else{
			GTKU_WARNING("gtku_object_exist_chk()", "Unknown Object name.");
		}

	}
	va_end(ap);

	return(truth);
}

