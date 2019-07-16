/*********************************************************************
 * gtkgl_image.c
 *
 * Copyright (C) 2002 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: gtkgl_image.c 415 2005-07-15 06:04:24Z sasaoka $ */

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

/* Image Object */
RC draw_image_object(GTKU_OBJECT *obj);
static RC draw_pixel_texture(GRAPHIC *gr, GTKU_DRAW_PROP *d_prop);
static RC init_texture(GRAPHIC *gr, GTKU_DRAW_PROP *d_prop);
static void trans_rgb24_GLimage(GRAPHIC *gr,  GTKU_DRAW_PROP *d_prop,
                                GLubyte **image);
static RC init_image_direction(GRAPHIC *gr, GTKU_DRAW_PROP *d_prop);
static RC draw_pixel_cell(GRAPHIC *gr, GTKU_DRAW_PROP *d_prop);
static void draw_pixel_cell_single(int height, int x, int y);
static RC draw_voxel_cell(GRAPHIC *gr, GTKU_DRAW_PROP *d_prop);
/*
static void draw_voxel_cell_single(int x, int y, int z);
*/


RC
gtku_gl_draw_image_object (GTKU_OBJECT *obj)
{
	GTKU_DRAW_TYPE draw_type;

	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	RC_NULL_CHK(obj);
	draw_type = obj->d_prop->draw_type;

	if(draw_type & GTKU_DRAW_HIDE) return(NORMAL_RC);

	if(draw_type & GTKU_DRAW_BBOX)
		RC_TRY( gtku_gl_draw_bbox(obj->d_prop) );

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
static RC
draw_pixel_texture (GRAPHIC *gr, GTKU_DRAW_PROP *d_prop)
{
	double tcx, tcy; /* Texture Coord */
	VECT3D p;

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
	glEnable(GL_CULL_FACE);

	return(NORMAL_RC);
}


/*
 * ここで行っている(ビット演算など)処理の意味は，
 * 「OpenGLプログラミングガイド」(9章)デクスチャマッピングの章を参照．
 * テクスチャーは必ず2^nの画像サイズでなければいけない．
 */
static RC
init_texture (GRAPHIC *gr, GTKU_DRAW_PROP *d_prop)
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

	if( (unsigned long)(1<<x_size_bit) == (gr->size.x) )
		d_prop->pack_size.x = 1<<x_size_bit;
	else
		d_prop->pack_size.x = 1<<(x_size_bit+1);

	if( (unsigned long)(1<<y_size_bit) == (gr->size.y) )
		d_prop->pack_size.y = 1<<y_size_bit;
	else
		d_prop->pack_size.y = 1<<(y_size_bit+1);

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
		image = (GLubyte *)mm_alloc(sizeof(GLubyte)*image_size*3);
		if(!image) return(ALLOC_ERROR_RC);
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


static void
trans_rgb24_GLimage (GRAPHIC *gr,  GTKU_DRAW_PROP *d_prop, GLubyte **image)
{
	unsigned long ii1, ii2, ii3;

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


static RC
init_image_direction (GRAPHIC *gr, GTKU_DRAW_PROP *d_prop)
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


static RC
draw_pixel_cell (GRAPHIC *gr, GTKU_DRAW_PROP *d_prop)
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


static void
draw_pixel_cell_single (int height, int x, int y)
{
	int Y = height - y - 1;
	glBegin(GL_QUADS);
	glVertex2i( (x    ), (Y    ) );
	glVertex2i( (x + 1), (Y    ) );
	glVertex2i( (x + 1), (Y + 1) );
	glVertex2i( (x    ), (Y + 1) );
	glEnd();
}


static RC
draw_voxel_cell (GRAPHIC *gr, GTKU_DRAW_PROP *d_prop)
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

	if(!gtku_gl->contour){
		gtku_gl->contour
			= (GRAPHIC_CONTOUR *)mm_alloc(sizeof(GRAPHIC_CONTOUR));
		if(!gtku_gl->contour) return(ALLOC_ERROR_RC);;
		RC_TRY( make_graphic_contour(gtku_gl->contour) );
	}

	if(gr->type == GRAPHIC_SCALAR){
		RC_TRY( chk_graphic_range_scalar(*gr, &min, &max, &avg) );
	}else if(gr->type == GRAPHIC_MONO16){
		long lmin, lmax, lavg;
		RC_TRY( chk_graphic_range_mono16(*gr, &lmin, &lmax, &lavg) );
		min = (double)lmin;
		max = (double)lmax;
	}

	/* translucent start */
	if(d_prop->img_trans_sw) set_gl_translucent();

	glDisable(GL_LIGHTING);
	glPointSize((float)d_prop->point_size);

	srand(3);
	GRAPHIC_SIZE_LOOP_3D(gr->size.z, gr->size.y, gr->size.x){
		if(d_prop->quantity == 0.0) return(NORMAL_RC);
		if(d_prop->quantity < 1.0){
			double random = (double)rand()/RAND_MAX;
			if( (d_prop->quantity - random) < 0.0)
				continue;
		}

		switch(gr->type){
		case GRAPHIC_RGB24:
			color = gr->data.rgb24[ii1][ii2][ii3];
			glColor4d( color.r/(GLdouble)0xff,
			           color.g/(GLdouble)0xff,
			           color.b/(GLdouble)0xff,
			           (GLdouble)d_prop->img_alpha );
			break;
		case GRAPHIC_MONO16:
			{
			  value = (double)gr->data.mono16[ii1][ii2][ii3];
			  if( ((value-min)/(max-min) < d_prop->img_threshold_min) ||
			      ((value-min)/(max-min) > d_prop->img_threshold_max) )
			  	continue;
			  RC_TRY( get_graphic_contour_rgb24( value, min, max,
			              *(gtku_gl->contour), d_prop->cont_type, &color) );
			  glColor4d( color.r/(GLdouble)0xff,
			             color.g/(GLdouble)0xff,
			             color.b/(GLdouble)0xff,
			             (GLdouble)d_prop->img_alpha );
			}
			break;
		case GRAPHIC_SCALAR:
			{
			  value = gr->data.scalar[ii1][ii2][ii3];
			  if( ((value-min)/(max-min) < d_prop->img_threshold_min) ||
			      ((value-min)/(max-min) > d_prop->img_threshold_max) )
			  	continue;
			  RC_TRY( get_graphic_contour_rgb24(value, min, max,
			                 *(gtku_gl->contour), d_prop->cont_type, &color) );
			  glColor4d( color.r/(GLdouble)0xff,
			             color.g/(GLdouble)0xff,
			             color.b/(GLdouble)0xff,
			             (GLdouble)d_prop->img_alpha );
			}
			break;
		default:
			GTKU_WARNING("draw_voxel_cell()", "Unknown image type.");
		}
		/*
		draw_voxel_cell_single(ii3, gr->size.y - ii2 - 1, ii1);
		*/
		glBegin(GL_POINTS);
		glVertex3d( ii3+0.5, (gr->size.y - ii2 - 0.5), (ii1 + 0.5) );
		glEnd();

	}GRAPHIC_SIZE_LOOP_3D_END;

	/* translucent end */
	if(d_prop->img_trans_sw) unset_gl_translucent();

	return(NORMAL_RC);
}


/* See element.eps, set_fem_face_table_hexa1() */
#if 0
static void
draw_voxel_cell_single(int x, int y, int z)
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
#endif
