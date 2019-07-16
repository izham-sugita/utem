/*********************************************************************
 * gtk_toolbar.c
 *
 * Copyright (C) 2002 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: gtk_toolbar.c 1066 2016-05-19 09:53:33Z hayashi $ */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <glib.h>
#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>
#include <gtk/gtkgl.h>
#include "base.h"
#include "mathematics.h"
#include "fem.h"
#include "gtk_utl.h"
#include "gtk_icons.h"

/* Global variable */
extern GTKU_GL *gtku_gl;

/* local global variable */
static GtkWidget *draw_type_radio[7];

/* Toolbar Rotation */
gboolean gtku_toolbar_rotation(GtkWidget *window, GtkWidget *toolbar);
static gboolean set_model_view_direction_xy(GtkWidget *widget,
                                            gpointer cb_data);
static gboolean set_model_view_direction_yz(GtkWidget *widget,
                                            gpointer cb_data);
static gboolean set_model_view_direction_zx(GtkWidget *widget,
                                            gpointer cb_data);
static gboolean set_model_view_direction_xyz(GtkWidget *widget,
                                             gpointer cb_data);
/* Toolbar Drawwing Type */
static gboolean gtku_toolbar_drawtype(GtkWidget *window, GtkWidget *toolbar);
static gboolean set_model_draw_type_polygon(GtkWidget *widget,
                                            gpointer cb_data);
static gboolean set_model_draw_type_solid(GtkWidget *widget,
                                          gpointer cb_data);
static gboolean set_model_draw_type_wire_frame(GtkWidget *widget,
                                               gpointer cb_data);
static gboolean set_model_draw_type_grid_frame(GtkWidget *widget,
                                               gpointer cb_data);
static gboolean set_model_draw_type_surf_point(GtkWidget *widget,
                                               gpointer cb_data);
static gboolean set_model_draw_type_grid_point(GtkWidget *widget,
                                               gpointer cb_data);
static gboolean set_model_draw_type_hiddenline(GtkWidget *widget,
                                               gpointer cb_data);

/* gtk_model_view.c : Toolbar Clear */
gboolean gtku_toolbar_clear(GtkWidget *widget, GtkWidget *toolbar);
gboolean gtku_set_model_clear(GtkWidget *widget, gpointer cb_data);

/* gtk_model_view.c : Toolbar Scale */
static gboolean gtku_toolbar_scale(GtkWidget *window, GtkWidget *toolbar);
static gboolean set_model_scale_init(GtkWidget *widget, gpointer cb_data);
static gboolean set_model_scale_fit(GtkWidget *widget, gpointer cb_data);
static gboolean set_model_scale_zoomin(GtkWidget *widget, gpointer cb_data);
static gboolean set_model_scale_zoomout(GtkWidget *widget, gpointer cb_data);

/* gtk_model_view.c : Toolbar Animation */
GtkWidget *gtku_toolbar_animation(GtkWidget *widget, GtkWidget *toolbar);
gboolean gtku_set_model_animation(GtkWidget *widget, gpointer cb_data);


/* ツールバーウィジェットを作成する．
 */
GtkWidget *
gtku_create_toolbar(GtkWidget *window)
{
	GtkWidget *handlebox;
	GtkWidget *toolbar;

	/* Create Handlebox */
	handlebox = gtk_handle_box_new();
	g_return_val_if_fail(handlebox, NULL);
	gtk_widget_show(handlebox);

	/* ツールバーを作成する． */
	toolbar = gtk_toolbar_new();
	g_return_val_if_fail(toolbar, NULL);
	/** アイコンを水平方向に並べる． **/
	gtk_toolbar_set_orientation(GTK_TOOLBAR(toolbar),
	                            GTK_ORIENTATION_HORIZONTAL);
	/** アイコンとラベルの両方を表示する． **/
	gtk_toolbar_set_style(GTK_TOOLBAR(toolbar), GTK_TOOLBAR_BOTH);
	gtk_container_add(GTK_CONTAINER(handlebox), toolbar);
	gtk_widget_show(toolbar);

	/* 以下，ボタンのタイプごとにツールバーへの登録を行う． */
	/** 操作を戻すボタン **/
	g_return_val_if_fail(gtku_toolbar_clear(window, toolbar), NULL);
	gtk_toolbar_append_space(GTK_TOOLBAR(toolbar));

	/** 拡大・縮小系 **/
	g_return_val_if_fail(gtku_toolbar_scale(window, toolbar), NULL);
	gtk_toolbar_append_space(GTK_TOOLBAR(toolbar));

	/** 描画タイプ系 **/
	g_return_val_if_fail(gtku_toolbar_drawtype(window, toolbar), NULL);
	gtk_toolbar_append_space(GTK_TOOLBAR(toolbar));

	/** アニメーション系 **/
	g_return_val_if_fail(gtku_toolbar_animation(window, toolbar), NULL);
	gtk_toolbar_append_space(GTK_TOOLBAR(toolbar));

	/** 座標系回転系 **/
	g_return_val_if_fail(gtku_toolbar_rotation(window, toolbar), NULL);
	/*
	gtk_toolbar_append_space(GTK_TOOLBAR(toolbar));
	*/

	return(handlebox);
}


/* 座標系回転系のツールバーのボタンを作り，ツールバーに登録する．
 */
gboolean
gtku_toolbar_rotation (GtkWidget *widget, GtkWidget *toolbar)
{
	/* FRONT view */
	g_return_val_if_fail( gtk_toolbar_append_element( GTK_TOOLBAR(toolbar),
			GTK_TOOLBAR_CHILD_BUTTON, NULL, NULL, "Front View", NULL,
			gtku_create_gtkimage_from_xpm(widget, (gchar **)xy_xpm),
			GTK_SIGNAL_FUNC(set_model_view_direction_xy), NULL), FALSE);

	/* RIGHT view */
	g_return_val_if_fail( gtk_toolbar_append_element( GTK_TOOLBAR(toolbar),
			GTK_TOOLBAR_CHILD_BUTTON, NULL, NULL, "Side View", NULL,
			gtku_create_gtkimage_from_xpm(widget, (gchar **)yz_xpm),
			GTK_SIGNAL_FUNC(set_model_view_direction_yz), NULL), FALSE);

	/* TOP view */
	g_return_val_if_fail( gtk_toolbar_append_element( GTK_TOOLBAR(toolbar),
			GTK_TOOLBAR_CHILD_BUTTON, NULL, NULL, "Top View", NULL,
			gtku_create_gtkimage_from_xpm(widget, (gchar **)zx_xpm),
			GTK_SIGNAL_FUNC(set_model_view_direction_zx), NULL), FALSE);

	/* XYZ view */
	g_return_val_if_fail( gtk_toolbar_append_element( GTK_TOOLBAR(toolbar),
			GTK_TOOLBAR_CHILD_BUTTON, NULL, NULL, "XYZ View", NULL,
			gtku_create_gtkimage_from_xpm(widget, (gchar **)xyz_xpm),
			GTK_SIGNAL_FUNC(set_model_view_direction_xyz), NULL), FALSE);

	return(TRUE);
}

/* FRONT view */
static gboolean
set_model_view_direction_xy(GtkWidget *widget, gpointer cb_data)
{
	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	init_operation_matrix(gtku_gl->rot_mat);
	rotation_operation_matrix_y(0.0, gtku_gl->rot_mat);
	gtk_widget_queue_draw(gtku_gl->glarea);

	return(TRUE);
}

/* RIGHT view */
static gboolean
set_model_view_direction_yz(GtkWidget *widget, gpointer cb_data)
{
	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	init_operation_matrix(gtku_gl->rot_mat);
	rotation_operation_matrix_y(-PI/2, gtku_gl->rot_mat);
	gtk_widget_queue_draw(gtku_gl->glarea);

	return(TRUE);
}

/* TOP view */
static gboolean
set_model_view_direction_zx(GtkWidget *widget, gpointer cb_data)
{
	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);
	init_operation_matrix(gtku_gl->rot_mat);
	rotation_operation_matrix_x(PI/2, gtku_gl->rot_mat);
	gtk_widget_queue_draw(gtku_gl->glarea);

	return(TRUE);
}

/* XYZ view */
static gboolean
set_model_view_direction_xyz(GtkWidget *widget, gpointer cb_data)
{
	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	init_operation_matrix(gtku_gl->rot_mat);
	rotation_operation_matrix_y(-PI/4, gtku_gl->rot_mat);
	rotation_operation_matrix_x( PI/4, gtku_gl->rot_mat);
	gtk_widget_queue_draw(gtku_gl->glarea);

	return(TRUE);
}


/* 描画タイプ系のツールバーのボタンを作り，ツールバーに登録する．
 */
static gboolean
gtku_toolbar_drawtype(GtkWidget *widget, GtkWidget *toolbar)
{
	draw_type_radio[0] = gtk_toolbar_append_element(GTK_TOOLBAR(toolbar),
		GTK_TOOLBAR_CHILD_RADIOBUTTON, NULL, NULL, "Polygon", NULL,
		gtku_create_gtkimage_from_xpm(widget, (gchar **)polygon_xpm),
		(GtkSignalFunc)set_model_draw_type_polygon, NULL);
	g_return_val_if_fail(draw_type_radio[0], FALSE);

	draw_type_radio[1] = gtk_toolbar_append_element(GTK_TOOLBAR(toolbar),
		GTK_TOOLBAR_CHILD_RADIOBUTTON, draw_type_radio[0], NULL, "Solid",
		NULL, gtku_create_gtkimage_from_xpm(widget, (gchar **)solid_xpm),
		(GtkSignalFunc)set_model_draw_type_solid, NULL);
	g_return_val_if_fail(draw_type_radio[1], FALSE);

	draw_type_radio[2] = gtk_toolbar_append_element(GTK_TOOLBAR(toolbar),
		GTK_TOOLBAR_CHILD_RADIOBUTTON, draw_type_radio[1], NULL, "Wire Frame",
		NULL, gtku_create_gtkimage_from_xpm(widget, (gchar **)wire_frame_xpm),
		(GtkSignalFunc)set_model_draw_type_wire_frame, NULL);
	g_return_val_if_fail(draw_type_radio[2], FALSE);

	draw_type_radio[3] = gtk_toolbar_append_element(GTK_TOOLBAR(toolbar),
		GTK_TOOLBAR_CHILD_RADIOBUTTON, draw_type_radio[2], NULL, "Grid Frame",
		NULL, gtku_create_gtkimage_from_xpm(widget, (gchar **)grid_frame_xpm),
		(GtkSignalFunc)set_model_draw_type_grid_frame, NULL);
	g_return_val_if_fail(draw_type_radio[3], FALSE);

	draw_type_radio[4] = gtk_toolbar_append_element(GTK_TOOLBAR(toolbar),
		GTK_TOOLBAR_CHILD_RADIOBUTTON, draw_type_radio[3], NULL,
		"Surface Point", NULL,
		gtku_create_gtkimage_from_xpm(widget, (gchar **)surf_point_xpm),
		(GtkSignalFunc)set_model_draw_type_surf_point, NULL);
	g_return_val_if_fail(draw_type_radio[4], FALSE);

	draw_type_radio[5] = gtk_toolbar_append_element(GTK_TOOLBAR(toolbar),
		GTK_TOOLBAR_CHILD_RADIOBUTTON, draw_type_radio[4], NULL, "Grid Point",
		NULL, gtku_create_gtkimage_from_xpm(widget, (gchar **)grid_point_xpm),
		(GtkSignalFunc)set_model_draw_type_grid_point, NULL);
	g_return_val_if_fail(draw_type_radio[5], FALSE);

	draw_type_radio[6] = gtk_toolbar_append_element(GTK_TOOLBAR(toolbar),
		GTK_TOOLBAR_CHILD_RADIOBUTTON, draw_type_radio[5], NULL, "Hidden Line",
		NULL, gtku_create_gtkimage_from_xpm(widget, (gchar **)hiddenline_xpm),
		(GtkSignalFunc)set_model_draw_type_hiddenline, NULL);
	g_return_val_if_fail(draw_type_radio[6], FALSE);

	/* default on polygon */
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(draw_type_radio[0]), TRUE);

	return(TRUE);
}


static gboolean
set_model_draw_type_polygon(GtkWidget *widget, gpointer cb_data)
{
	g_return_val_if_fail(gtku_gl, FALSE);
	if(gtku_gl->draw_type != GTKU_DRAW_POLYGON){
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(draw_type_radio[0]),
		                             TRUE);
		gtku_gl->draw_type = GTKU_DRAW_POLYGON;
	}
	gtk_widget_queue_draw(gtku_gl->glarea);
	return(TRUE);
}


static gboolean
set_model_draw_type_solid(GtkWidget *widget, gpointer cb_data)
{                   
	g_return_val_if_fail(gtku_gl, FALSE);
	if(gtku_gl->draw_type != GTKU_DRAW_SOLID){
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(draw_type_radio[1]),
		                             TRUE);
		gtku_gl->draw_type = GTKU_DRAW_SOLID;
	}
	gtk_widget_queue_draw(gtku_gl->glarea);
	return(TRUE);
}


static gboolean
set_model_draw_type_wire_frame(GtkWidget *widget, gpointer cb_data)
{                   
	g_return_val_if_fail(gtku_gl, FALSE);
	if(gtku_gl->draw_type != GTKU_DRAW_WIRE_FRAME){
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(draw_type_radio[2]),
		                             TRUE);
		gtku_gl->draw_type = GTKU_DRAW_WIRE_FRAME;
	}
	gtk_widget_queue_draw(gtku_gl->glarea);
	return(TRUE);
}


static gboolean
set_model_draw_type_grid_frame(GtkWidget *widget, gpointer cb_data)
{                   
	g_return_val_if_fail(gtku_gl, FALSE);
	if(gtku_gl->draw_type != GTKU_DRAW_GRID_FRAME){
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(draw_type_radio[3]),
		                             TRUE);
		gtku_gl->draw_type = GTKU_DRAW_GRID_FRAME;
	}
	gtk_widget_queue_draw(gtku_gl->glarea);
	return(TRUE);
}


static gboolean
set_model_draw_type_surf_point(GtkWidget *widget, gpointer cb_data)
{                   
	g_return_val_if_fail(gtku_gl, FALSE);
	if(gtku_gl->draw_type != GTKU_DRAW_SURF_POINT){
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(draw_type_radio[4]),
		                             TRUE);
		gtku_gl->draw_type = GTKU_DRAW_SURF_POINT;
	}
	gtk_widget_queue_draw(gtku_gl->glarea);
	return(TRUE);
}


static gboolean
set_model_draw_type_grid_point(GtkWidget *widget, gpointer cb_data)
{                   
	g_return_val_if_fail(gtku_gl, FALSE);
	if(gtku_gl->draw_type != GTKU_DRAW_GRID_POINT){
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(draw_type_radio[5]),
		                             TRUE);
		gtku_gl->draw_type = GTKU_DRAW_GRID_POINT;
	}
	gtk_widget_queue_draw(gtku_gl->glarea);
	return(TRUE);
}


static gboolean
set_model_draw_type_hiddenline(GtkWidget *widget, gpointer cb_data)
{                   
	g_return_val_if_fail(gtku_gl, FALSE);
	if(gtku_gl->draw_type != GTKU_DRAW_HIDDENLINE){
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(draw_type_radio[6]),
		                             TRUE);
		gtku_gl->draw_type = GTKU_DRAW_HIDDENLINE;
	}
	gtk_widget_queue_draw(gtku_gl->glarea);
	return(TRUE);
}


/* 操作を戻すボタンを作り，ツールバーに登録する．
 */
gboolean
gtku_toolbar_clear(GtkWidget *widget, GtkWidget *toolbar)
{
	GtkWidget *elem;
	elem = gtk_toolbar_append_element(
	           GTK_TOOLBAR(toolbar), GTK_TOOLBAR_CHILD_BUTTON, NULL,
	           NULL, "Clear", NULL, gtk_image_new_from_stock(GTK_STOCK_CLEAR,
	           GTK_ICON_SIZE_LARGE_TOOLBAR),
	           GTK_SIGNAL_FUNC(gtku_set_model_clear), NULL );
	g_return_val_if_fail(elem, FALSE);
	return(TRUE);
}


gboolean
gtku_set_model_clear (GtkWidget *widget, gpointer cb_data)
{
	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);
	g_return_val_if_fail(gtku_gl_realize(gtku_gl->glarea), FALSE);
	gtk_widget_queue_draw(gtku_gl->glarea);
	return(TRUE);
}


gboolean
gtku_toolbar_scale (GtkWidget *widget, GtkWidget *toolbar)
{
	GtkWidget *elem;

	elem = gtk_toolbar_append_element(GTK_TOOLBAR(toolbar),
		GTK_TOOLBAR_CHILD_BUTTON, NULL, NULL, "Initialize", NULL,
		gtk_image_new_from_stock(GTK_STOCK_ZOOM_100,
		GTK_ICON_SIZE_LARGE_TOOLBAR),
		GTK_SIGNAL_FUNC(set_model_scale_init), NULL);
	g_return_val_if_fail(elem, FALSE);

	elem = gtk_toolbar_append_element(GTK_TOOLBAR(toolbar),
		GTK_TOOLBAR_CHILD_BUTTON, NULL, NULL, "Fit", NULL,
		gtk_image_new_from_stock(GTK_STOCK_ZOOM_FIT,
		GTK_ICON_SIZE_LARGE_TOOLBAR),
		GTK_SIGNAL_FUNC(set_model_scale_fit), NULL);
	g_return_val_if_fail(elem, FALSE);

	elem = gtk_toolbar_append_element(GTK_TOOLBAR(toolbar),
		GTK_TOOLBAR_CHILD_BUTTON, NULL, NULL, "Zoom IN", NULL,
		gtk_image_new_from_stock(GTK_STOCK_ZOOM_IN,
		GTK_ICON_SIZE_LARGE_TOOLBAR),
		GTK_SIGNAL_FUNC(set_model_scale_zoomin), NULL);
	g_return_val_if_fail(elem, FALSE);

	elem = gtk_toolbar_append_element(GTK_TOOLBAR(toolbar),
		GTK_TOOLBAR_CHILD_BUTTON, NULL, NULL, "Zoom OUT", NULL,
		gtk_image_new_from_stock(GTK_STOCK_ZOOM_OUT,
		GTK_ICON_SIZE_LARGE_TOOLBAR),
		GTK_SIGNAL_FUNC(set_model_scale_zoomout), NULL);
	g_return_val_if_fail(elem, FALSE);

	return(TRUE);
}


static gboolean
set_model_scale_init(GtkWidget *widget, gpointer cb_data)
{
	g_return_val_if_fail(gtku_gl, FALSE);

	init_operation_matrix(gtku_gl->scale_mat);
	init_operation_matrix(gtku_gl->shift_mat);
	gtk_widget_queue_draw(gtku_gl->glarea);

	return(TRUE);
}

static gboolean
set_model_scale_fit(GtkWidget *widget, gpointer cb_data)
{
	double fit_size = 0.0;

	g_return_val_if_fail(gtku_gl, FALSE);
	if(gtku_gl->glarea->allocation.width
		< gtku_gl->glarea->allocation.height){
		fit_size = (double)gtku_gl->glarea->allocation.width
		           / gtku_gl->glarea_w_size;
	}else{
		fit_size = (double)gtku_gl->glarea->allocation.height
		           / gtku_gl->glarea_h_size;
	}
	init_operation_matrix(gtku_gl->scale_mat);
	scale_operation_matrix(fit_size,fit_size,fit_size,gtku_gl->scale_mat);
	gtk_widget_queue_draw(gtku_gl->glarea);

	return(TRUE);
}


static gboolean
set_model_scale_zoomin(GtkWidget *widget, gpointer cb_data)
{
	double zoom_p = 1.1;

	g_return_val_if_fail(gtku_gl, FALSE);
	scale_operation_matrix(zoom_p, zoom_p, zoom_p, gtku_gl->scale_mat);
	gtk_widget_queue_draw(gtku_gl->glarea);

	return(TRUE);
}


static gboolean
set_model_scale_zoomout(GtkWidget *widget, gpointer cb_data)
{
	double zoom_m = 0.9;

	g_return_val_if_fail(gtku_gl, FALSE);

	scale_operation_matrix(zoom_m, zoom_m, zoom_m, gtku_gl->scale_mat);
	gtk_widget_queue_draw(gtku_gl->glarea);

	return(TRUE);
}


GtkWidget*
gtku_toolbar_animation (GtkWidget *widget, GtkWidget *toolbar)
{
	GtkWidget *elem;
	elem = gtk_toolbar_append_element(
	           GTK_TOOLBAR(toolbar), GTK_TOOLBAR_CHILD_TOGGLEBUTTON,
	           NULL, NULL, "Animation", NULL, gtku_create_gtkimage_from_xpm(
	           widget, (gchar **)anim_xpm),
	           GTK_SIGNAL_FUNC(gtku_set_model_animation), NULL);
	g_return_val_if_fail(elem, NULL);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(elem), FALSE);
	return(elem);
}


gboolean
gtku_set_model_animation (GtkWidget *widget, gpointer cb_data)
{
	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	gtku_gl->anim_lock_sw = !(gtku_gl->anim_lock_sw);
	g_return_val_if_fail(gtku_gl_toggle_animation(gtku_gl->glarea), FALSE);

	gtk_widget_queue_draw(gtku_gl->glarea);

	return(TRUE);
}


