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

/* $Id: gtkgl_utl.c 864 2007-01-17 10:02:04Z sasaoka $ */

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
#include "memory_manager.h"


#define TIMEOUT_INTERVAL (10)
#define ANIMATE_THRESHOLD (100.0)
#define ANIMETION_LIMIT_THRESHOLD (5.0)
#define MOVE_THRESHOLD (1.0)

static gchar font_name[] = "courier 12";  /* Fix!! 決め打ち */

/* Global variable */
extern GTKU_GL *gtku_gl;

/* TOP Level OpenGL function */
static gboolean is_render_mode(GTKU_RENDER_MODE mode);
static gboolean select_object(GLuint select_buf[]);
static void draw_gradation_backgroud(void);
static void set_gl_projection_matrix(void);
static void set_gl_model_view_matirx(GTKU_OBJECT *obj);
static void set_gl_light_env(void);
static RC draw_object(GTKU_OBJECT *obj);

/* GLAREA, Interface, Event ... */
static RC init_objects_drawing_attribute(void);
static RC init_objects_operation_matrix(void);
static RC init_view_volume(void);
static RC cal_objects_size(GTKU_OBJECT_ARRAY *object, VECT3D *center,
                           VECT3D *max, VECT3D *min);
static RC cal_object_size_fem(NODE_ARRAY *node, VECT3D *center,
                              VECT3D *max, VECT3D *min);
static RC cal_object_size_image(GRAPHIC *gr, VECT3D *center,
                                VECT3D *max, VECT3D *min);
static RC cal_objects_view_scale(double *scale);

static RC glarea_mouse_selection(GtkWidget *widget);
static RC glarea_mouse_rotation(GtkWidget *widget, GdkRectangle area);
static RC mouse_rotation_local(GTKU_DRAW_PROP *d_prop, double **rot);
static RC glarea_mouse_shift(GtkWidget *widget, GdkRectangle area);
static RC mouse_shift_local(GTKU_DRAW_PROP *d_prop, VECT3D disp);
static RC glarea_mouse_scale(GtkWidget *widget, GdkRectangle area,
                             GdkModifierType state);
static RC mouse_scale_local(GTKU_DRAW_PROP *d_prop, VECT3D rate);
static gboolean add_rotation_matrix_for_animation(GtkWidget *widget);
static gboolean glarea_key_scale(GtkWidget *widget, GdkEventKey *event);
static gboolean timeout(GtkWidget *widget);
static gboolean timeout_add(GtkWidget *widget);
static gboolean timeout_remove(GtkWidget *widget);
static GdkRectangle init_gdk_rectangle(void);


/* The "expose_event" signal handler. This is repeatedly called as the */
/* painting routine every time the "expose/draw" event is signalled.   */
gboolean
gtku_gl_draw (GtkWidget *widget, GdkEventExpose *event)
{
	int ii1;
	GLuint select_buf[512];
	GdkGLContext *glcontext;
	GdkGLDrawable *gldrawable;

	/* Draw only last expose. */
	if(event->count > 0) return(TRUE); 

	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	switch(gtku_gl->render_mode){
	case GTKU_RENDER_NONE:
		return(TRUE);
	case GTKU_RENDER_SELECT:
	case GTKU_RENDER_GL:
		glcontext  = gtk_widget_get_gl_context(widget);
		gldrawable = gtk_widget_get_gl_drawable(widget);
		break;
	case GTKU_RENDER_AREA:
		gdk_draw_drawable(GDK_DRAWABLE(widget->window),
		                  widget->style->mid_gc[GTK_WIDGET_STATE(widget)],
		                  GDK_DRAWABLE(gtku_gl->gl_pixmap), 0,0, 0,0,
		                  widget->allocation.width,
		                  widget->allocation.height);

		gdk_draw_rectangle(GDK_DRAWABLE(widget->window),
		                   widget->style->mid_gc[GTK_WIDGET_STATE(widget)],
		                   FALSE,
		                   gtku_gl->select_area.x,
		                   gtku_gl->select_area.y,
		                   gtku_gl->select_area.width,
		                   gtku_gl->select_area.height);
		return(TRUE);
	case GTKU_RENDER_PIXMAP:
		{
		  /* 描画エリアを GdkPixmapに変更 */
		  GdkGLConfig *glconfig;
		  glconfig = gdk_gl_config_new_by_mode(GDK_GL_MODE_RGB |
		                                       GDK_GL_MODE_DEPTH |
		                                       GDK_GL_MODE_SINGLE);
		  gldrawable = GDK_GL_DRAWABLE( gdk_pixmap_set_gl_capability(
		                                 gtku_gl->gl_pixmap, glconfig, NULL) );
		  glcontext = GDK_GL_CONTEXT( gdk_gl_context_new(gldrawable, NULL,
		                              TRUE, GDK_GL_RGBA_TYPE) );
		  g_object_unref( G_OBJECT(glconfig) );
		  /* 描画エリアが違うので光源設定やりなおし */
		  gtku_gl->init_light = FALSE;
		}
		break;
	default:
		GTKU_WARNING("gtku_gl_draw()", "Unknown Render Type");
		return(FALSE);
	}

	g_return_val_if_fail(gldrawable, FALSE);
	g_return_val_if_fail(glcontext, FALSE);

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

		if( is_render_mode(GTKU_RENDER_SELECT) ){
			glSelectBuffer(512, select_buf);
			glRenderMode(GL_SELECT);
			glInitNames();
			glPushName(0);
		}

		/* initialize color and depth buffer. clean screen */
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		/* Animation */
		if(gtku_gl->animation && !gtku_gl->anim_lock_sw)
			add_rotation_matrix_for_animation(widget);

		/* Draw Gradation Backgroud */
		draw_gradation_backgroud();

		/* PROJECTION MATRIX */
		set_gl_projection_matrix();

		/* Worl Plane */
		if(gtku_gl->draw_plane_sw && !is_render_mode(GTKU_RENDER_SELECT) ){
			set_gl_model_view_matirx(NULL);
			GTKU_BOOLEAN_RC_TRY( gtku_gl_draw_world_plane() );
		}

		/* Drawing Object  */
		for(ii1=0; ii1<gtku_gl->object->size; ii1++){
			GTKU_OBJECT *obj = gtku_gl->object->obj[ii1];
			set_gl_model_view_matirx(obj);
			if( is_render_mode(GTKU_RENDER_SELECT) ){
				glLoadName(ii1);
				GTKU_BOOLEAN_RC_TRY( draw_object(obj) );
			}else{
				if(gtku_gl->draw_center_mark_sw){
					double length;
					length = (  abs_vect3d(obj->d_prop->max)
					          - abs_vect3d(obj->d_prop->min) ) / 15.0;
					gtku_gl_draw_center_marker(obj->d_prop->center, length);
				}
				if((gtku_gl->motion)&&(obj->d_prop->motion_bbox_sw)){
					RC_TRY( gtku_gl_draw_bbox(obj->d_prop) );
				}else{
					GTKU_BOOLEAN_RC_TRY( draw_object(obj) );
				}
			}
			gdk_gl_drawable_wait_gl(gldrawable);
		}

		/* Coordinate System */
		if(gtku_gl->draw_coord_sw && !is_render_mode(GTKU_RENDER_SELECT) )
			gtku_gl_draw_coodinate_axis(); 

		if( is_render_mode(GTKU_RENDER_SELECT) ){
			glMatrixMode(GL_PROJECTION);
			glPopMatrix();
		}


		/* Pick up selected object */
		if( is_render_mode(GTKU_RENDER_SELECT) ){
			g_return_val_if_fail(select_object(select_buf), FALSE);
			gtku_gl->render_mode = GTKU_RENDER_GL;
		}

		/* Sync GL object rendering */
		gdk_gl_drawable_wait_gl(gldrawable);

		/* Swap backbuffer to front */
		if( gdk_gl_drawable_is_double_buffered(gldrawable) )
			gdk_gl_drawable_swap_buffers(gldrawable);
		else
			glFlush();

		gdk_gl_drawable_gl_end(gldrawable);
	}/* OpenGL END */


	if( is_render_mode(GTKU_RENDER_PIXMAP) ){
		g_object_unref( G_OBJECT(glcontext) );
		gtku_gl->render_mode = GTKU_RENDER_AREA;
	}

	/* Sync GDK rendering */
	gdk_gl_drawable_wait_gdk(gldrawable);

#if 0
	{ /* name stack size */
	 GLint value;
	 glGetIntegerv(GL_NAME_STACK_DEPTH, &value);
	 fprintf(stderr, "GL_NAME_STACK_DEPTH = %d\n",value);
	}
#endif

	return(TRUE);
}


static gboolean
is_render_mode (GTKU_RENDER_MODE mode)
{
	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);
	return(gtku_gl->render_mode == mode);
}


static gboolean
select_object (GLuint select_buf[])
{
	int ii1;
	unsigned int ii2;
	GLuint *ptr, names;
	GLint hits = glRenderMode(GL_RENDER);

	ptr = (GLuint *)select_buf;
	printf("hits = %d\n", hits);
	if(hits < 0) return(TRUE);
	names = *ptr;
	for(ii1=0; ii1<hits; ii1++){
		names = *ptr;
		printf("number of names =  %d\n", names); ptr++;
		printf("z1 =  %g ", (float)(*ptr)/0x7fffffff); ptr++;
		printf("z2 =  %g\n", (float)(*ptr)/0x7fffffff); ptr++;
		for(ii2=0; ii2<names; ii2++){
			printf("%d ", *ptr); ptr++;
		}
		printf("\n");
	}
	printf("\n");

	return(TRUE);
}

static void
draw_gradation_backgroud (void)
{
	int ii1, ii2;
	unsigned char bg_buf[4*64*64];
	GLfloat color_grad[3];

	color_grad[0] = (gtku_gl->bg_color1[0] - gtku_gl->bg_color0[0])/64;
	color_grad[1] = (gtku_gl->bg_color1[1] - gtku_gl->bg_color0[1])/64;
	color_grad[2] = (gtku_gl->bg_color1[2] - gtku_gl->bg_color0[2])/64;

	for(ii1=0; ii1<64; ii1++){
		unsigned char color[3];
		color[0] = (unsigned char)
		           (255 * (gtku_gl->bg_color0[0] + ii1 * color_grad[0]) );
		color[1] = (unsigned char)
		           (255 * (gtku_gl->bg_color0[1] + ii1 * color_grad[1]) );
		color[2] = (unsigned char)
		           (255 * (gtku_gl->bg_color0[2] + ii1 * color_grad[2]) );

		for(ii2=ii1*64*4; ii2<(ii1+1)*64*4;){
			bg_buf[ii2++] = color[0];
			bg_buf[ii2++] = color[1];
			bg_buf[ii2++] = color[2];
			bg_buf[ii2++] = 1;
		}
	}

	glDisable(GL_DEPTH_TEST);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0, gtku_gl->view_size.x*2,
	           0, gtku_gl->view_size.y*2);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glRasterPos2i(0, 0);
	glPixelZoom((GLfloat)(gtku_gl->view_size.x * 2 / 64),
	            (GLfloat)(gtku_gl->view_size.y * 2 / 64)); 
	glDrawPixels(64, 64, GL_RGBA, GL_UNSIGNED_BYTE, bg_buf);
	glClearDepth(1);
	glEnable(GL_DEPTH_TEST);
}


/* projection matrix setting */
static void
set_gl_projection_matrix (void)
{
	/* GTKU_GL existence check */
	g_return_if_fail(gtku_gl);

	glMatrixMode(GL_PROJECTION);

	if( is_render_mode(GTKU_RENDER_SELECT) ) glPushMatrix();

	glLoadIdentity();

	if( is_render_mode(GTKU_RENDER_SELECT) ){
		GLint viewport[4];
		glGetIntegerv(GL_VIEWPORT, viewport);
		gluPickMatrix( (GLdouble)gtku_gl->select_area.x,
		               (GLdouble)(viewport[3] - gtku_gl->select_area.y),
		               (GLdouble)gtku_gl->select_area.width,
		               (GLdouble)gtku_gl->select_area.height, viewport );
	}

	glOrtho(-(gtku_gl->view_size.x), (gtku_gl->view_size.x),
	        -(gtku_gl->view_size.y), (gtku_gl->view_size.y),
	        -(gtku_gl->view_size.z), (gtku_gl->view_size.z));

}


/* view matrix transformation */
static void
set_gl_model_view_matirx (GTKU_OBJECT *obj)
{
	double m[4][4];
	double **mat;

	/* GTKU_GL existence check */
	g_return_if_fail(gtku_gl);

	GTKU_VOID_RC_TRY( allocate_operation_matrix(&mat) );

	/* local matrix operation */
	if(obj != NULL){
		mul_operation_matrix(obj->d_prop->rot_mat,   mat);
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
	gtku_gl_conv_form_matrix(mat, m);
	glMultMatrixd(*m);

	free_operation_matrix(&mat);
}


/* object lighting and background color */
static void
set_gl_light_env (void)
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
	glClearColor(gtku_gl->bg_color0[0], gtku_gl->bg_color0[1],
	             gtku_gl->bg_color0[2], gtku_gl->bg_color0[3]);

	/* 拡大縮小などの Matrix 操作を行うとき、OpenGL 側に法線ベクトルの */
	/* 単位化をさせる必要がある。glEnable(GL_NORMALIZE) は自動正規化を */
	/* ONにしている                                                    */
	glEnable(GL_NORMALIZE);

	/* GLUのTeapotなどで使われるNUBAS局面のの法線を自動計算する */
	glEnable(GL_AUTO_NORMAL);
}


void
set_gl_translucent (void)
{
	glEnable(GL_BLEND);
	glDepthMask(FALSE);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}


void
unset_gl_translucent (void)
{
	glDepthMask(TRUE);
	glDisable(GL_BLEND);
}


static RC
draw_object (GTKU_OBJECT *obj)
{
	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	if( (obj->type & GTKU_OBJECT_NONE) ){
		return(NORMAL_RC);
	}else if(obj->type & (GTKU_OBJECT_IMAGE_PIXEL|GTKU_OBJECT_IMAGE_VOXEL)){
		RC_TRY( gtku_gl_draw_image_object(obj) ); /* see -> gtkgl_image.c */
	}else if(obj->type & GTKU_OBJECT_FEM_NODE){/* element does not nessarily */
		RC_TRY( gtku_gl_draw_fem_object(obj) ); /* see -> gtkgl_fem.c */
	}else{
		GTKU_WARNING("draw_object()", "Unknown Object");
	}

	return(NORMAL_RC);
}


/* The "configure_event" signal handler. Almost always it will be used */
/* to resize the OpenGL viewport when the window is resized.           */
gboolean
gtku_gl_reshape (GtkWidget *widget, GdkEventConfigure *event)
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

	if(gtku_gl->gl_pixmap != NULL)
		g_object_unref( G_OBJECT(gtku_gl->gl_pixmap) );

	gtku_gl->gl_pixmap = gdk_pixmap_new(widget->window,
	                                    widget->allocation.width,
	                                    widget->allocation.height, -1);

	return(TRUE);
}


/*
 * The "realize" signal handler.  Usually, it will be called one time.
 * All the OpenGL initialization should be performed here normally. But,
 * this function only calculate model(s) size and model center position
 * and so on.
 */
gboolean
gtku_gl_realize (GtkWidget *widget)
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

	/* calculate object view size and so on. */
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

	/* set default rendering mode */ 
	gtku_gl->render_mode = GTKU_RENDER_GL;

	return TRUE;
}


static RC
init_objects_drawing_attribute (void)
{
	int ii1;
	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	for(ii1=0; ii1<gtku_gl->object->size; ii1++){
		gtku_gl->object->obj[ii1]->d_prop->init_image_dir = FALSE;
	}

	return(NORMAL_RC);
}


static RC
init_objects_operation_matrix (void)
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
		GTKU_DRAW_PROP *d_prop = gtku_gl->object->obj[ii1]->d_prop;
		init_operation_matrix(d_prop->shift_mat);
		init_operation_matrix(d_prop->scale_mat);
		init_operation_matrix(d_prop->rot_mat);
	}

	return(NORMAL_RC);
}


static RC
init_view_volume (void)
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


static RC
cal_objects_size (GTKU_OBJECT_ARRAY *object,
                  VECT3D *center, VECT3D *max, VECT3D *min)
{
	int ii1;

	for(ii1=0; ii1<object->size; ii1++){
		GTKU_OBJECT *obj = object->obj[ii1];
		RC_NULL_CHK(obj);
		if(obj->type == GTKU_OBJECT_NONE) continue;
		else if(obj->type & GTKU_OBJECT_FEM_NODE){
			RC_TRY( cal_object_size_fem(obj->node, &obj->d_prop->center,
			                            &obj->d_prop->max, &obj->d_prop->min) );
		}else if( obj->type & (GTKU_OBJECT_IMAGE_PIXEL |
			                   GTKU_OBJECT_IMAGE_VOXEL) ){
			RC_TRY( cal_object_size_image(obj->gr, &obj->d_prop->center,
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
		GTKU_DRAW_PROP *d_prop = object->obj[ii1]->d_prop;
		max->x = MAX2(d_prop->max.x, max->x);
		max->y = MAX2(d_prop->max.y, max->y);
		max->z = MAX2(d_prop->max.z, max->z);
		min->x = MIN2(d_prop->min.x, min->x);
		min->y = MIN2(d_prop->min.y, min->y);
		min->z = MIN2(d_prop->min.z, min->z);
	}
	*center = mid_point3d(*max, *min);

	return(NORMAL_RC);
}


static RC
cal_object_size_fem (NODE_ARRAY *node, VECT3D *center, VECT3D *max, VECT3D *min)
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


static RC
cal_object_size_image (GRAPHIC *gr, VECT3D *center, VECT3D *max, VECT3D *min)
{
	RC_TRY( init_vect3d(min) );
	max->x = (double)gr->size.x;
	max->y = (double)gr->size.y;
	max->z = (double)gr->size.z;
	*center = mid_point3d(*max, *min);

	return(NORMAL_RC);
}


static RC
cal_objects_view_scale (double *scale)
{
	double length;
	VECT3D view_box;

	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	view_box.x = gtku_gl->object->max.x - gtku_gl->object->min.x;
	view_box.y = gtku_gl->object->max.y - gtku_gl->object->min.y;
	view_box.z = gtku_gl->object->max.z - gtku_gl->object->min.z;
	length = abs_vect3d(view_box);
	if(gtku_gl->glarea_w_size < gtku_gl->glarea_h_size){
		*scale = gtku_gl->glarea_w_size/length;
	}else{
		*scale = gtku_gl->glarea_h_size/length;
	}

	return(NORMAL_RC);
}


/* The "motion_notify_event" signal handler. It will be called while  */
/* mouse cursor moving motion in widget area.                         */
gboolean
gtku_gl_motion_notify (GtkWidget *widget, GdkEventMotion *event)
{
	GdkRectangle area;
	GdkModifierType state;
	gboolean redraw = FALSE;
	gint x, y;

	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	if(gtku_gl->animation) return(TRUE);

	if(event->is_hint)
		gdk_window_get_pointer(event->window, &x, &y, &state);
	else{
		x = (gint)event->x;
		y = (gint)event->y;
		state = event->state;
	}

	gtku_gl->cur.x = (double)x;
	gtku_gl->cur.y = (double)y;

	gtku_gl->d.x = gtku_gl->cur.x - gtku_gl->pre.x;
	if(fabs(gtku_gl->d.x) < MOVE_THRESHOLD)
		gtku_gl->d.x = 0.0;

	gtku_gl->d.y = gtku_gl->pre.y - gtku_gl->cur.y; /* GDKでは画面左上が(0,0) */
	if(fabs(gtku_gl->d.y) < MOVE_THRESHOLD)
		gtku_gl->d.y = 0.0;

	area.width  = widget->allocation.width;
	area.height = widget->allocation.height;

	if(state & GDK_BUTTON1_MASK){
		/* Left Button [Selection] */
		GTKU_BOOLEAN_RC_TRY(glarea_mouse_selection(widget));
		redraw = TRUE;
	}else if(state & GDK_BUTTON2_MASK){
		/* Middle Button [Rotation] */
		GTKU_BOOLEAN_RC_TRY(glarea_mouse_rotation(widget, area));
		redraw = TRUE;
	}else if(state & GDK_BUTTON3_MASK){
		/* Right Button [Translation] */
		GTKU_BOOLEAN_RC_TRY(glarea_mouse_shift(widget, area));
		redraw = TRUE;
	}else if(state & (GDK_BUTTON4_MASK|GDK_BUTTON5_MASK)){
		/* Scroll UP and Down [Scale] */
		GTKU_BOOLEAN_RC_TRY(glarea_mouse_scale(widget, area, state ));
		redraw = TRUE;
	}

	if(redraw && !gtku_gl->animation) gtk_widget_queue_draw(widget);

	gtku_gl->pre.x = gtku_gl->cur.x ;
	gtku_gl->pre.y = gtku_gl->cur.y ;

	return(TRUE);
}


/* motion_notify event sub call back */
static RC
glarea_mouse_selection (GtkWidget *widget)
{
	GdkRectangle *select;
	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	select = &(gtku_gl)->select_area;

	select->x = (gint)MIN2(gtku_gl->cur.x, gtku_gl->select_start.x);
	select->y = (gint)MIN2(gtku_gl->cur.y, gtku_gl->select_start.y);
	select->width  = (gint)fabs(gtku_gl->cur.x - gtku_gl->select_start.x);
	select->height = (gint)fabs(gtku_gl->cur.y - gtku_gl->select_start.y);

	return(NORMAL_RC);
}


#define ROT_ANG (PI/6)
/* motion_notify event sub call back */
static RC
glarea_mouse_rotation (GtkWidget *widget, GdkRectangle area)
{
	double **rot;
	double r, angle, angle1, angle2;


	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	if( (gtku_gl->d.x == 0.0) && (gtku_gl->d.y == 0.0) )
		return(NORMAL_RC);

	RC_TRY( allocate_operation_matrix(&rot) );

	/* Key Press Effect */
	r = abs_vect3d(gtku_gl->d);
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
		/* グローバル座標で回転させる場合
		if( (gtku_gl->key.x == TRUE) || (gtku_gl->key.y == TRUE) ||
		    (gtku_gl->key.z == TRUE) ){
			RC_TRY( mul_operation_matrix(gtku_gl->rot_mat, rot) );
			RC_TRY( gtku_gl_conv_form_matrix(4, 4, rot, gtku_gl->rot_mat) );
		}
		*/
	}else{ /* local rotation */
		int ii1;
		for(ii1=0; ii1<gtku_gl->object->size; ii1++){
			GTKU_DRAW_PROP *d_prop = gtku_gl->object->obj[ii1]->d_prop;
			if(d_prop->rot_lock_sw == TRUE) continue;
			RC_TRY( mouse_rotation_local(d_prop, rot) );
		}
	}

	free_operation_matrix(&rot);

	return(NORMAL_RC);
}


static RC
mouse_rotation_local (GTKU_DRAW_PROP *d_prop, double **rot)
{
	double **center, **center_inv;

	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	RC_TRY( allocate_operation_matrix(&center) );
	RC_TRY( allocate_operation_matrix(&center_inv) );
	shift_operation_matrix( -d_prop->center.x, -d_prop->center.y,
	                        -d_prop->center.z, center);
	shift_operation_matrix(  d_prop->center.x, d_prop->center.y,
	                         d_prop->center.z, center_inv);
	RC_TRY( mul_operation_matrix(center, d_prop->rot_mat) );

	if(d_prop->use_local_axis_sw){
		RC_TRY( mul_operation_matrix(rot, d_prop->rot_mat) );
	}else{
		double **rot_inv;
		RC_TRY( allocate_operation_matrix(&rot_inv) );
		RC_TRY( inverse_matrix(4, gtku_gl->rot_mat, rot_inv) );
		RC_TRY( mul_operation_matrix(gtku_gl->rot_mat, d_prop->rot_mat) );
		RC_TRY( mul_operation_matrix(rot,              d_prop->rot_mat) );
		RC_TRY( mul_operation_matrix(rot_inv,          d_prop->rot_mat) );
		free_operation_matrix(&rot_inv);
	}

	RC_TRY( mul_operation_matrix(center_inv, d_prop->rot_mat) );
	free_operation_matrix(&center);
	free_operation_matrix(&center_inv);

	return(NORMAL_RC);
}


/* motion_notify event sub call back    */
/* for 5 button mouse. anyone need this */
static RC
glarea_mouse_scale (GtkWidget *widget, GdkRectangle area,
                    GdkModifierType state)
{
	int ii1;
	VECT3D rate;
	double d_rate;

	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	if( fabs(gtku_gl->d.x) > fabs(gtku_gl->d.y) ){
		d_rate = abs_vect3d(gtku_gl->d)/area.width;
		if(gtku_gl->d.x < 0) d_rate *= -1;
	}else{
		d_rate = abs_vect3d(gtku_gl->d)/area.height;
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
		/* グローバル座標で伸縮させる場合
		if( (gtku_gl->key.x == TRUE) || (gtku_gl->key.y == TRUE) ||
		    (gtku_gl->key.z == TRUE) ){
			double **rot_inv;
			RC_TRY( allocate_operation_matrix(&rot_inv) );
			RC_TRY( inverse_matrix(4, gtku_gl->rot_mat, rot_inv));
			RC_TRY( mul_operation_matrix(rot_inv, gtku_gl->scale_mat) );
			scale_operation_matrix(rate.x, rate.y, rate.z, gtku_gl->scale_mat);
			RC_TRY( mul_operation_matrix(gtku_gl->rot_mat, gtku_gl->scale_mat));
			free_operation_matrix(&rot_inv);
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


static RC
mouse_scale_local (GTKU_DRAW_PROP *d_prop, VECT3D rate)
{

	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

#if 0
	double **center, **center_inv, **shift_inv, **ini_inv;
	RC_TRY( allocate_operation_matrix(&ini_inv) );
	RC_TRY( allocate_operation_matrix(&center) );
	RC_TRY( allocate_operation_matrix(&center_inv) );
	RC_TRY( allocate_operation_matrix(&shift_inv) );
	RC_TRY( inverse_matrix(4, gtku_gl->initial_mat, ini_inv) );
	RC_TRY( mul_operation_matrix(gtku_gl->initial_mat, d_prop->scale_mat) );
	RC_TRY( inverse_matrix(4, d_prop->shift_mat, shift_inv) );
	RC_TRY( mul_operation_matrix(shift_inv, center) );
	shift_operation_matrix( -d_prop->center.x, -d_prop->center.y,
	                        -d_prop->center.z, center);
	RC_TRY( mul_operation_matrix(center, d_prop->scale_mat) );
	RC_TRY( inverse_matrix(4, center, center_inv) );
#endif

	if(d_prop->use_local_axis_sw){
		scale_operation_matrix(rate.x, rate.y, rate.z, d_prop->scale_mat);
	}else{
		double **rot_inv;
		RC_TRY( allocate_operation_matrix(&rot_inv) );
		RC_TRY( inverse_matrix(4, gtku_gl->rot_mat, rot_inv));
		RC_TRY( mul_operation_matrix(gtku_gl->rot_mat, d_prop->scale_mat) );
		scale_operation_matrix(rate.x, rate.y, rate.z, d_prop->scale_mat);
		RC_TRY( mul_operation_matrix(rot_inv,          d_prop->scale_mat) );
		free_operation_matrix(&rot_inv);
	}

#if 0
	RC_TRY( mul_operation_matrix(center_inv, d_prop->scale_mat) );
	RC_TRY( mul_operation_matrix(ini_inv, d_prop->scale_mat) );
	free_operation_matrix(&center);
	free_operation_matrix(&center_inv);
	free_operation_matrix(&shift_inv);
	free_operation_matrix(&ini_inv);
#endif

	return(NORMAL_RC);
}


/* motion_notify event sub call back */
static RC
glarea_mouse_shift (GtkWidget *widget, GdkRectangle area)
{
	VECT3D d;

	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	/* Initialize */
	RC_TRY( init_vect3d(&d) ); 

	/* Key Press Effect */
	     if(gtku_gl->key.x == TRUE) d.x = gtku_gl->d.x;
	else if(gtku_gl->key.y == TRUE) d.y = gtku_gl->d.y;
	else if(gtku_gl->key.z == TRUE){
		d.z = abs_vect3d(gtku_gl->d);
		if( fabs(gtku_gl->d.x) > fabs(gtku_gl->d.y) ){
			if(gtku_gl->d.x < 0.0)
				d.z *= -1.0;
		}else{
			if(gtku_gl->d.y < 0.0)
				d.z *= -1.0;
		}
	}else if(gtku_gl->key.shift == TRUE){
		if( fabs(gtku_gl->d.x) > fabs(gtku_gl->d.y) )
			d.x = gtku_gl->d.x;
		else
			d.y = gtku_gl->d.y;
	}else{ /* NO Key */
		d.x = gtku_gl->d.x;
		d.y = gtku_gl->d.y;
	}
	
	{ /* 拡大率によって移動量を変更する */
		double **scale_inv;
		RC_TRY( allocate_operation_matrix(&scale_inv) );
		RC_TRY( inverse_matrix(4, gtku_gl->scale_mat, scale_inv));
		d =  operation3d(scale_inv, d);
		free_operation_matrix(&scale_inv);
	}

	if(gtku_gl->shift_lock_count == 0){ /* global shift */
		/* グローバル座標で移動させる場合
		if( (gtku_gl->key.x == TRUE) || (gtku_gl->key.y == TRUE) ||
		    (gtku_gl->key.z == TRUE) ){
			double **rot_inv;
			RC_TRY( allocate_operation_matrix(&rot_inv) );
			RC_TRY( inverse_matrix(4, gtku_gl->rot_mat, rot_inv));
			RC_TRY( mul_operation_matrix(rot_inv, gtku_gl->shift_mat) );
			shift_operation_matrix(d.x, d.y, d.z, gtku_gl->shift_mat);
			RC_TRY( mul_operation_matrix(gtku_gl->rot_mat, gtku_gl->shift_mat));
			free_operation_matrix(&rot_inv);
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


static RC
mouse_shift_local (GTKU_DRAW_PROP *d_prop, VECT3D disp)
{
	double **ini_inv;

	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	RC_TRY( allocate_operation_matrix(&ini_inv) );
	RC_TRY( inverse_matrix(4, gtku_gl->initial_mat, ini_inv));
	RC_TRY( mul_operation_matrix(gtku_gl->initial_mat, d_prop->shift_mat) );

	{ /* 拡大率によって移動量を変更する */
	  double **scale_inv;
	  RC_TRY( allocate_operation_matrix(&scale_inv) );
	  RC_TRY( inverse_matrix(4, d_prop->scale_mat, scale_inv) );
	  disp = operation3d(scale_inv, disp);
	  free_operation_matrix(&scale_inv);
	}

	if(d_prop->use_local_axis_sw){
		shift_operation_matrix(disp.x, disp.y, disp.z, d_prop->shift_mat);
	}else{
		double **rot_inv;
		RC_TRY( allocate_operation_matrix(&rot_inv) );
		RC_TRY( inverse_matrix(4, gtku_gl->rot_mat, rot_inv));
		RC_TRY( mul_operation_matrix(gtku_gl->rot_mat, d_prop->shift_mat) );
		shift_operation_matrix(disp.x, disp.y, disp.z, d_prop->shift_mat);
		RC_TRY( mul_operation_matrix(rot_inv, d_prop->shift_mat) );
		free_operation_matrix(&rot_inv);
	}

	RC_TRY( mul_operation_matrix(ini_inv, d_prop->shift_mat) );
	free_operation_matrix(&ini_inv);

	return(NORMAL_RC);
}


static gboolean
add_rotation_matrix_for_animation (GtkWidget *widget)
{
	GdkRectangle area;

	area.width  = widget->allocation.width;
	area.height = widget->allocation.height;
	glarea_mouse_rotation(widget, area);

	return(TRUE);
}


/* The "button_press_event" signal handler. */
gboolean
gtku_gl_button_press (GtkWidget *widget, GdkEventButton *event)
{
	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	if(gtku_gl->animation){ 
		if( (event->button == 1)||(event->button == 2)||(event->button == 3) )
			gtku_gl_toggle_animation(widget);
	}else{

		if(event->button == 1){
			gtku_gl->render_mode = GTKU_RENDER_PIXMAP;
			gtku_gl->select_start.x = (gint)event->x;
			gtku_gl->select_start.y = (gint)event->y;
			gtku_gl->select_area = init_gdk_rectangle();
		}else if( (event->button == 2)||(event->button == 3) ){
			gtku_gl->motion = TRUE;
		}
		gtku_gl->cur.x = event->x;
		gtku_gl->cur.y = event->y;
	}

#if 0
	GdkCursor* cursor;  
	/* unref NEED */
	if(event->button == 2){
		cursor = gdk_cursor_new(GDK_EXCHANGE);
		gdk_window_set_cursor(widget->window, cursor);
	}else{
		cursor = gdk_cursor_new(GDK_X_CURSOR);
		gdk_window_set_cursor(widget->window, cursor);
	}
#endif

	if(!gtku_gl->animation) gtk_widget_queue_draw(widget);

	return(TRUE);
}


/* The "button_press_release" signal handler. */
gboolean
gtku_gl_button_release (GtkWidget *widget, GdkEventButton *event)
{
	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	if( (event->button == 1)||(event->button == 2)||(event->button == 3) ){
		gtku_gl->cur.x = event->x;
		gtku_gl->cur.y = event->y;
		if(event->button == 1)
			gtku_gl->render_mode = GTKU_RENDER_SELECT;
		else if( (event->button == 2)||(event->button == 3) )
			gtku_gl->motion = FALSE;
	}

	if(gtku_gl->animation){ 
		return(TRUE);
	}else if(!gtku_gl->anim_lock_sw){
		VECT3D d = gtku_gl->d;
		double n = d.x*d.x + d.y*d.y;
		if( (event->button == 2) && (n > ANIMATE_THRESHOLD) )
			gtku_gl_toggle_animation(widget);
	}

	if(!gtku_gl->animation) gtk_widget_queue_draw(widget);

	return(TRUE);
}


/* The "scroll" signal handler. */
gboolean
gtku_gl_scroll (GtkWidget *widget, GdkEventScroll *event)
{
	VECT3D rate;
	double d_rate = 0.0;

	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	if(event->direction == GDK_SCROLL_UP)   d_rate =  0.005;
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
gboolean
gtku_gl_map (GtkWidget *widget, GdkEvent *event)
{
	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);
	if(gtku_gl->animation) timeout_add(widget);
	return(TRUE);
}


/* The "unmap_event" signal handler. */
gboolean
gtku_gl_unmap (GtkWidget *widget, GdkEvent *event)
{
	timeout_remove(widget);
	return(TRUE);
}


/* The "visibility_notify_event" signal handler. */
gboolean
gtku_gl_visibility_notify (GtkWidget *widget, GdkEventVisibility *event)
{
	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);
	if(gtku_gl->animation){
		if(event->state == GDK_VISIBILITY_FULLY_OBSCURED)
			timeout_remove(widget);
	}else timeout_add(widget);
	return(TRUE);
}


gboolean
gtku_gl_destroy (GtkWidget *widget)
{
	RC rc;
	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);
	if(gtku_gl->animation){
		timeout_remove(widget);
		rc = mm_free(gtku_gl);
		if(rc != NORMAL_RC) return(FALSE);
	}
	return(TRUE);
}


/* The "key_press_event" signal handler. */
gboolean
gtku_gl_key_press (GtkWidget *widget, GdkEventKey *event)
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
gboolean
gtku_gl_key_release (GtkWidget *widget,GdkEventKey *event)
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
static gboolean
glarea_key_scale (GtkWidget *widget, GdkEventKey *event)
{
	VECT3D rate;
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


static gboolean
timeout (GtkWidget *widget)
{
	gtk_widget_queue_draw(widget);
	return(TRUE);
}


static gboolean
timeout_add (GtkWidget *widget)
{
	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);
	if(gtku_gl->timeout_id == 0){
		gtku_gl->timeout_id = gtk_timeout_add(TIMEOUT_INTERVAL,
		                                     (GtkFunction)timeout, widget);
	}
	return(TRUE);
}


static gboolean
timeout_remove (GtkWidget *widget)
{
	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);
	if(gtku_gl->timeout_id != 0){
		gtk_timeout_remove(gtku_gl->timeout_id);
		gtku_gl->timeout_id = 0;
    }
  	return(TRUE);
}


gboolean
gtku_gl_toggle_animation (GtkWidget *widget)
{
	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);
	gtku_gl->animation = !(gtku_gl->animation);

	if(gtku_gl->animation){
		if( (fabs(gtku_gl->d.x) > ANIMETION_LIMIT_THRESHOLD) ||
		    (fabs(gtku_gl->d.y) > ANIMETION_LIMIT_THRESHOLD) ){
		double r = abs_vect3d(gtku_gl->d)/ANIMETION_LIMIT_THRESHOLD;
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


void
gtku_gl_draw_string (VECT3D p, const char *str)
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


void
gtku_gl_conv_form_matrix (double **src, double dest[4][4])
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


/* "..." = "node", "element", "surface", "pixcel", "voxel" */
gboolean
gtku_object_exist_chk (GTKU_OBJECT *obj, int size, ...)
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


static GdkRectangle
init_gdk_rectangle (void)
{
	GdkRectangle rect;

	rect.x = 0;
	rect.y = 0;
	rect.width = 0;
	rect.height = 0;

	return(rect);
}


