/*********************************************************************
 * gtkgl_article.c
 *
 * Copyright (C) 2002 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: gtkgl_article.c 842 2006-12-04 06:21:20Z sasaoka $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <glib.h>
#include <gtk/gtk.h>
#include <gtk/gtkgl.h>
#include "gtk_utl.h"
#include "math_utl.h"
#include "rc.h"

/* Global variable */
extern GTKU_GL *gtku_gl;

void
gtku_gl_draw_coodinate_axis (void)
{
	double m[4][4];
	double **mat;
	VECT3D p; /* Position of the coordinate system Name */

	p.x = -4.0;
	p.y = -4.0;
	p.z = 60.0;
	
	/* GTKU_GL existence check */
	g_return_if_fail(gtku_gl);

	GTKU_VOID_RC_TRY( allocate_operation_matrix(&mat) );
	mul_operation_matrix(gtku_gl->rot_mat, mat);
	gtku_gl_conv_form_matrix(mat, m);

	/*
	glDepthFunc(GL_ALWAYS);
	*/

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslated(-(gtku_gl->view_size.x - 67.0),
	             -(gtku_gl->view_size.y - 67.0), 0.0);

	/* Origin of the coordinate system. Sphere */
	glPushMatrix();
	glMultMatrixd(*m);
	glEnable(GL_LIGHTING);
	glPolygonMode(GL_FRONT, GL_FILL);
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
	glPolygonMode(GL_FRONT, GL_FILL);
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
	glPolygonMode(GL_FRONT, GL_FILL);
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
	glPolygonMode(GL_FRONT, GL_FILL);
	gdk_gl_draw_cone(TRUE, 8.0, 50, 10, 10);
	glColor4d(1.0, 1.0, 1.0, 1.0);  /* White */
	gtku_gl_draw_string(p, "Z");
	glPopMatrix();

	free_operation_matrix(&mat);

	/*
	glDepthFunc(GL_LEQUAL);
	*/
}


void
gtku_gl_draw_center_marker (VECT3D center, double length)
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
RC
gtku_gl_draw_world_plane (void)
{
	double ii1, ii2, ii3, max_dist, d;
	VECT3D p, min, max;

	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);
	
	min = gtku_gl->object->min;
	max = gtku_gl->object->max;
	p = (abs_vect3d(max) > abs_vect3d(min) ? max : min);
	max_dist = MAX2(MAX2(fabs(p.x), fabs(p.y)), fabs(p.z));
	d = max_dist/WORLD_PLANE_DIVIDE;


	glDisable(GL_LIGHTING);
	for(ii1=0; ii1<2; ii1++){
		if(ii1==0){ /* Plane */
			set_gl_translucent();
			glColor4d(0.7, 0.7, 0.7, 0.3);  /* Gray */
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glDisable(GL_CULL_FACE);
		}else{ /* Line */
			glColor4d(0.8, 0.8, 0.8, 1.0);  /* Gray */
			glLineWidth(1.0);
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		}
		for(ii2=-max_dist; ii2<max_dist-d/2.0; ii2+=d){
			for(ii3=-max_dist; ii3<max_dist-d/2.0; ii3+=d){
				if(gtku_gl->draw_plane_yz_sw){
					glBegin(GL_POLYGON);
						glVertex4d(0.0, ii2,   ii3,   WORLD_PLANE_FACTER);
						glVertex4d(0.0, ii2+d, ii3,   WORLD_PLANE_FACTER);
						glVertex4d(0.0, ii2+d, ii3+d, WORLD_PLANE_FACTER);
						glVertex4d(0.0, ii2,   ii3+d, WORLD_PLANE_FACTER);
					glEnd();
				}
				if(gtku_gl->draw_plane_zx_sw){
					glBegin(GL_POLYGON);
						glVertex4d(ii2,   0.0, ii3,   WORLD_PLANE_FACTER);
						glVertex4d(ii2,   0.0, ii3+d, WORLD_PLANE_FACTER);
						glVertex4d(ii2+d, 0.0, ii3+d, WORLD_PLANE_FACTER);
						glVertex4d(ii2+d, 0.0, ii3,   WORLD_PLANE_FACTER);
					glEnd();
				}
				if(gtku_gl->draw_plane_xy_sw){
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
			if(gtku_gl->draw_plane_yz_sw){
				glBegin(GL_LINES);
					glVertex4d(0.0, -max_dist, 0.0, WORLD_PLANE_FACTER);
					glVertex4d(0.0,  max_dist, 0.0, WORLD_PLANE_FACTER);
					glVertex4d(0.0, 0.0, -max_dist, WORLD_PLANE_FACTER);
					glVertex4d(0.0, 0.0,  max_dist, WORLD_PLANE_FACTER);
				glEnd();
			}
			if(gtku_gl->draw_plane_zx_sw){
				glBegin(GL_LINES);
					glVertex4d(-max_dist, 0.0, 0.0, WORLD_PLANE_FACTER);
					glVertex4d( max_dist, 0.0, 0.0, WORLD_PLANE_FACTER);
					glVertex4d(0.0, 0.0, -max_dist, WORLD_PLANE_FACTER);
					glVertex4d(0.0, 0.0,  max_dist, WORLD_PLANE_FACTER);
				glEnd();
			}
			if(gtku_gl->draw_plane_xy_sw){
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


RC
gtku_gl_draw_bbox (GTKU_DRAW_PROP *prop)
{
	VECT3D min, max;

	RC_NULL_CHK(prop);

	min = prop->min;
	max = prop->max;

	glDisable(GL_LIGHTING);
	glLineWidth((float)prop->line_size);
	glColor4dv(prop->line_color);
	glBegin(GL_LINE_LOOP);
	 glVertex3d(min.x, min.y, min.z);
	 glVertex3d(max.x, min.y, min.z);
	 glVertex3d(max.x, max.y, min.z);
	 glVertex3d(min.x, max.y, min.z);
	glEnd();
	glBegin(GL_LINE_LOOP);
	 glVertex3d(min.x, min.y, max.z);
	 glVertex3d(max.x, min.y, max.z);
	 glVertex3d(max.x, max.y, max.z);
	 glVertex3d(min.x, max.y, max.z);
	glEnd();
	glBegin(GL_LINES);
	 glVertex3d(min.x, min.y, min.z);
	 glVertex3d(min.x, min.y, max.z);
	 glVertex3d(max.x, min.y, min.z);
	 glVertex3d(max.x, min.y, max.z);
	 glVertex3d(max.x, max.y, min.z);
	 glVertex3d(max.x, max.y, max.z);
	 glVertex3d(min.x, max.y, min.z);
	 glVertex3d(min.x, max.y, max.z);
	glEnd();

	return(NORMAL_RC);
}



