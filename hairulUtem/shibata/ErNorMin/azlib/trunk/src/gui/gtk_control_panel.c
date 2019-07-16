/*********************************************************************
 * gtk_control_panel.c
 *
 * Copyright (C) 2006 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: gtk_control_panel.c 842 2006-12-04 06:21:20Z sasaoka $ */

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

/* Global variable */
extern GTKU_GL *gtku_gl;

/* Matrix Lock */
static GtkWidget *object_ctl_panel_matrix(GTKU_OBJECT *obj);
static gboolean set_local_rot_matrix_lock(GtkWidget *widget,
                                          gpointer cb_data);
static gboolean set_local_scale_matrix_lock(GtkWidget *widget,
                                            gpointer cb_data);
static gboolean set_local_shift_matrix_lock(GtkWidget *widget,
                                            gpointer cb_data);

/* Drawing Type */
static GtkWidget *object_ctl_panel_draw_type(GTKU_OBJECT *obj);
static gboolean set_local_object_draw_type(GtkWidget *widget, gpointer cb_data);

/* Quantity of Display */
static GtkWidget *object_ctl_panel_quantity(GTKU_OBJECT *obj);
static gboolean set_local_quantity(GtkWidget *widget, gpointer cb_data);

/* Image Direction */
static GtkWidget *object_ctl_panel_image_dir(GTKU_OBJECT *obj);
static GtkWidget *local_image_direction_option(GTKU_OBJECT *obj);
static gboolean set_local_image_direction(GtkWidget *widget, gpointer cb_data);

/* Image Contour */
static GtkWidget *object_ctl_panel_image_contour(GTKU_OBJECT *obj);
static GtkWidget *local_image_contour_option(GTKU_OBJECT *obj);
static gboolean set_local_image_contour(GtkWidget *widget, gpointer cb_data);

/* Image Display Thrshold */
static GtkWidget *object_ctl_panel_image_range(GTKU_OBJECT *obj);
static gboolean set_local_image_range(GtkWidget *widget, gpointer cb_data);

/* TransLucent of Image */
static GtkWidget *object_ctl_panel_image_trans(GTKU_OBJECT *obj);
static gboolean set_local_image_trans(GtkWidget *widget, gpointer cb_data);

/* Color&Alpha Selection */
static GtkWidget *object_ctl_panel_color(GTKU_OBJECT *obj);
static GtkWidget *object_ctl_panel_trans_check(GTKU_OBJECT *obj);
static GtkWidget *object_ctl_panel_color_button(GTKU_OBJECT *obj);
static gboolean local_object_color_selection_dialog(GtkWidget *widget,
                                                  gpointer cb_data);
static gboolean set_local_object_color(GtkWidget *widget, gpointer cb_data);
static gboolean unset_local_object_color(GtkWidget *widget, gpointer cb_data);

/* Drawing Number */
static GtkWidget *object_ctl_panel_draw_num(GTKU_OBJECT *obj);

/* Size (Point, Line) */
static GtkWidget *object_ctl_panel_size(GTKU_OBJECT *obj);
static gboolean set_local_object_size(GtkWidget *widget, gpointer cb_data);

static gboolean set_boolean_switch(GtkWidget *widget, gpointer cb_data);

/* etc */
static gboolean init_op_mat(GtkWidget *widget, gpointer cb_data);
static gboolean file_selection_save_op_mat(GtkWidget *widget,
                                           gpointer cb_data);
static gboolean filew_ok_save_op_mat(GtkWidget *widget, gpointer cb_data);
static gboolean save_op_mat(GtkWidget *widget, gpointer cb_data);
static RC print_object_prop_operation_mat (FILE *fp, GTKU_DRAW_PROP *d_prop);
static gboolean file_selection_load_op_mat(GtkWidget *widget,
                                           gpointer cb_data);
static gboolean filew_ok_load_op_mat(GtkWidget *widget, gpointer cb_data);
static gboolean load_op_mat(GtkWidget *widget, gpointer cb_data);
static RC read_object_prop_operation_mat(FILE *fp, GTKU_DRAW_PROP *dprop);

static GtkWidget *local_image_range_with_label(GTKU_OBJECT *obj,
                                               double *range, char *label);
static GtkWidget *local_object_size_spin_with_label(GTKU_OBJECT *obj,
                                                    GLdouble *size,
                                                    char *label);


GtkWidget *
gtku_object_ctl_panel (GTKU_OBJECT *obj)
{
	GtkWidget *vbox;

	/* Create main Vbox */
	vbox = gtk_vbox_new(FALSE, 2);
	gtk_widget_show(vbox);

	{ /* Create operation matrix lock frame */
	 GtkWidget *frame = object_ctl_panel_matrix(obj);
	 g_return_val_if_fail(frame, NULL);
	 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(frame), FALSE, FALSE, 0);
	 gtk_widget_show(frame);
	}

	{ /* Create local Drawing type */
	 GtkWidget *frame = object_ctl_panel_draw_type(obj);
	 g_return_val_if_fail(frame, NULL);
	 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(frame), FALSE, FALSE, 0);
	 gtk_widget_show(frame);
	}

	{ /* Create local Quantity of Display */
	 GtkWidget *frame = object_ctl_panel_quantity(obj);
	 g_return_val_if_fail(frame, NULL);
	 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(frame), FALSE, FALSE, 0);
	 gtk_widget_show(frame);
	}

	if(obj->type & (GTKU_OBJECT_IMAGE_PIXEL | GTKU_OBJECT_IMAGE_VOXEL)){
		{ /* Create local Image Dirction */
		 GtkWidget *frame = object_ctl_panel_image_dir(obj);
		 g_return_val_if_fail(frame, NULL);
		 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(frame), FALSE, FALSE, 0);
		 gtk_widget_show(frame);
		}

		{ /* Create local Use Contour Map */
		 GtkWidget *frame = object_ctl_panel_image_contour(obj);
		 g_return_val_if_fail(frame, NULL);
		 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(frame), FALSE, FALSE, 0);
		 gtk_widget_show(frame);
		}

		{ /* Create local Image Range */
		 GtkWidget *frame = object_ctl_panel_image_range(obj);
		 g_return_val_if_fail(frame, NULL);
		 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(frame), FALSE, FALSE, 0);
		 gtk_widget_show(frame);
		}

		{ /* Create local Image Trans */
		 GtkWidget *frame = object_ctl_panel_image_trans(obj);
		 g_return_val_if_fail(frame, NULL);
		 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(frame), FALSE, FALSE, 0);
		 gtk_widget_show(frame);
		}
	}else{
		{ /* Create local Drawing Color & Trans */
		 GtkWidget *frame = object_ctl_panel_color(obj);
		 g_return_val_if_fail(frame, NULL);
		 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(frame), FALSE, FALSE, 0);
		 gtk_widget_show(frame);
		}

		{ /* Create local Drawing Number */
		 GtkWidget *frame = object_ctl_panel_draw_num(obj);
		 g_return_val_if_fail(frame, NULL);
		 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(frame), FALSE, FALSE, 0);
		 gtk_widget_show(frame);
		}

		{ /* Create local Drawing Size */
		 GtkWidget *frame = object_ctl_panel_size(obj);
		 g_return_val_if_fail(frame, NULL);
		 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(frame), FALSE, FALSE, 0);
		 gtk_widget_show(frame);
		}
	}

	return(vbox);
}


GtkWidget *
object_ctl_panel_matrix (GTKU_OBJECT *obj)
{
	GtkWidget *vbox;
	GtkWidget *frame;

	/* Create main Vbox */
	vbox = gtk_vbox_new(FALSE, 0);
	gtk_widget_show(vbox);

	 /* Create frame */
	frame = gtk_frame_new("Operation Matrix Lock");
	g_return_val_if_fail(frame, NULL);
	gtk_frame_set_label_align(GTK_FRAME(frame), 0, 0.5);
	gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_ETCHED_IN);
	gtk_container_add(GTK_CONTAINER(frame), GTK_WIDGET(vbox));
	gtk_container_set_border_width(GTK_CONTAINER(frame), 10);
	gtk_widget_show(frame);

	{ /* Create check button rotation lock */
	 GtkWidget *check = gtk_check_button_new_with_label("Rotation");
	 g_return_val_if_fail(check, NULL);
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(check),
	                              obj->d_prop->rot_lock_sw);
	 g_signal_connect(G_OBJECT(check), "toggled",
	                  G_CALLBACK(set_local_rot_matrix_lock),
	                  (gpointer)&(obj->d_prop->rot_lock_sw));
	 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(check), FALSE, FALSE, 1);
	 gtk_widget_show(check);
	}

	{ /* Create check button scale lock */
	 GtkWidget *check = gtk_check_button_new_with_label("Scale");
	 g_return_val_if_fail(check, NULL);
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(check),
	                              obj->d_prop->scale_lock_sw);
	 g_signal_connect(G_OBJECT(check), "toggled",
	                  G_CALLBACK(set_local_scale_matrix_lock),
	                  (gpointer)&(obj->d_prop->scale_lock_sw));
	 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(check), FALSE, FALSE, 1);
	 gtk_widget_show(check);
	}

	{ /* Create check button shift lock */
	 GtkWidget *check = gtk_check_button_new_with_label("Shift");
	 g_return_val_if_fail(check, NULL);
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(check),
	                              obj->d_prop->shift_lock_sw);
	 g_signal_connect(G_OBJECT(check), "toggled",
	                  G_CALLBACK(set_local_shift_matrix_lock),
	                  (gpointer)&(obj->d_prop->shift_lock_sw));
	 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(check), FALSE, FALSE, 1);
	 gtk_widget_show(check);
	}

	{ /* Create check button  use local axis */
	 GtkWidget *check = gtk_check_button_new_with_label("Use Local Axis");
	 g_return_val_if_fail(check, NULL);
	 gtk_toggle_button_set_state(GTK_TOGGLE_BUTTON(check),
	                             obj->d_prop->use_local_axis_sw);
	 g_signal_connect(G_OBJECT(check), "toggled",
	                  G_CALLBACK(set_boolean_switch),
	                  (gpointer)&(obj->d_prop->use_local_axis_sw));
	 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(check), FALSE, FALSE, 1);
	 gtk_widget_show(check);
	}

	{ /* Initialize Operation Matrix */
	  GtkWidget *button = gtk_button_new_with_label("Init Operation Matrix");
	  g_return_val_if_fail(button, NULL);
	  g_signal_connect(G_OBJECT(button), "clicked",
	                   G_CALLBACK(init_op_mat), (gpointer)obj);
	  gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(button),
	                     FALSE, FALSE, 5);
	  gtk_widget_show(button);
	}

	{ /* Save Operation Matrix */
	  GtkWidget *button = gtk_button_new_with_label("Save Operation Matrix");
	  g_return_val_if_fail(button, NULL);
	  g_signal_connect(G_OBJECT(button), "clicked",
	                   G_CALLBACK(file_selection_save_op_mat), (gpointer)obj);
	  gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(button),
	                     FALSE, FALSE, 5);
	  gtk_widget_show(button);
	}

	{ /* Load Operation Matrix */
	  GtkWidget *button = gtk_button_new_with_label("Load Operation Matrix");
	  g_return_val_if_fail(button, NULL);
	  g_signal_connect(G_OBJECT(button), "clicked",
	                   G_CALLBACK(file_selection_load_op_mat), (gpointer)obj);
	  gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(button),
	                     FALSE, FALSE, 5);
	  gtk_widget_show(button);
	}

	return(frame);
}


#define GTKU_SET_LOCAL_MATRIX_LOCK(func_name, count)\
gboolean func_name(GtkWidget *widget, gpointer cb_data)\
{\
    gboolean *sw;\
    g_return_val_if_fail(gtku_gl, FALSE);\
    sw = (gboolean *)cb_data;\
    if(*sw) count++;\
    else count--;\
    *sw = !(*sw);\
    return(TRUE);\
}

GTKU_SET_LOCAL_MATRIX_LOCK(set_local_rot_matrix_lock,
                           gtku_gl->rot_lock_count)
GTKU_SET_LOCAL_MATRIX_LOCK(set_local_scale_matrix_lock,
                           gtku_gl->scale_lock_count)
GTKU_SET_LOCAL_MATRIX_LOCK(set_local_shift_matrix_lock,
                           gtku_gl->shift_lock_count)


static gboolean
init_op_mat (GtkWidget *widget, gpointer cb_data)
{
	GTKU_OBJECT *obj = (GTKU_OBJECT *)cb_data;

	g_return_val_if_fail(obj, FALSE);

	/* initialize local operation matrix */
	log_printf(5, "=== Init Operation Matix Data .... ");
	init_operation_matrix(obj->d_prop->shift_mat);
	init_operation_matrix(obj->d_prop->scale_mat);
	init_operation_matrix(obj->d_prop->rot_mat);
	log_printf(5, "Done\n");

	gtk_widget_queue_draw(gtku_gl->glarea);

	return(TRUE);
}


static gboolean
file_selection_save_op_mat (GtkWidget *widget, gpointer cb_data)
{
	GtkWidget *filew;

	filew =  gtk_file_selection_new("Save Operation Matrix File");
	g_return_val_if_fail(filew, FALSE);
	g_object_set_data(G_OBJECT(filew), "gtku_object", cb_data);
	g_signal_connect(G_OBJECT(GTK_FILE_SELECTION(filew)->ok_button),
	                 "clicked", G_CALLBACK(filew_ok_save_op_mat),
	                 (gpointer)filew);
	g_signal_connect_swapped(G_OBJECT(GTK_FILE_SELECTION(filew)->cancel_button),
	                         "clicked", G_CALLBACK(gtk_widget_destroy),
	                         G_OBJECT(filew));
	gtk_widget_show(filew);

	return(TRUE);
}


static gboolean
filew_ok_save_op_mat (GtkWidget *widget, gpointer cb_data)
{
	if(save_op_mat(widget, cb_data)){
		gtk_widget_destroy(GTK_WIDGET(cb_data));
		gtk_widget_queue_draw(gtku_gl->glarea);
	}

	return TRUE;
}


static gboolean
save_op_mat (GtkWidget *widget, gpointer cb_data)
{
	FILE *fp;
	char *filename;
	GTKU_OBJECT *obj;

	/* get Object */
	obj = (GTKU_OBJECT *)g_object_get_data(G_OBJECT(cb_data), "gtku_object");
	g_return_val_if_fail(obj, FALSE);

	/* get File Name */
	filename = (char *)gtk_file_selection_get_filename(
	                   GTK_FILE_SELECTION(cb_data));

	log_printf(5, "Save Operation Matrix File Name : \"%s\"\n", filename);
	log_printf(5, "=== Write Operation Matix Data .... ");
	GTKU_BOOLEAN_RC_TRY( rc_fopen(filename, "w", &fp) );
	GTKU_BOOLEAN_RC_TRY( print_object_prop_operation_mat(fp, obj->d_prop) );
	GTKU_BOOLEAN_RC_TRY( rc_fclose(fp) );
	log_printf(5, "Done\n");

	return(TRUE);
}

static RC
print_object_prop_operation_mat (FILE *fp, GTKU_DRAW_PROP *d_prop)
{
	if(d_prop == NULL) return(ARG_ERROR_RC);
	if( (d_prop->rot_mat   == NULL) ||
			(d_prop->scale_mat == NULL) ||
			(d_prop->shift_mat == NULL) ) return(ARG_ERROR_RC);

	fprintf(fp, "Rotation Matrix: \n");
	RC_TRY( print_operation_matrix(fp, d_prop->rot_mat));
	fprintf(fp, "Shift Matrix: \n");
	RC_TRY( print_operation_matrix(fp, d_prop->shift_mat));
	fprintf(fp, "Scale Matrix: \n");
	RC_TRY( print_operation_matrix(fp, d_prop->scale_mat));

	return(NORMAL_RC);

}


static gboolean
file_selection_load_op_mat (GtkWidget *widget, gpointer cb_data)
{
	GtkWidget *filew;

	filew =  gtk_file_selection_new("Load Operation Matrix File");
	g_return_val_if_fail(filew, FALSE);
	g_object_set_data(G_OBJECT(filew), "gtku_object", cb_data);
	g_signal_connect(G_OBJECT(GTK_FILE_SELECTION(filew)->ok_button),
	                 "clicked", G_CALLBACK(filew_ok_load_op_mat),
	                 (gpointer)filew);
	g_signal_connect_swapped(G_OBJECT(GTK_FILE_SELECTION(filew)->cancel_button),
	                         "clicked", G_CALLBACK(gtk_widget_destroy),
	                         G_OBJECT(filew));
	gtk_widget_show(filew);

	return(TRUE);
}


static gboolean
filew_ok_load_op_mat (GtkWidget *widget, gpointer cb_data)
{
	if(load_op_mat(widget, cb_data)){
		gtk_widget_destroy(GTK_WIDGET(cb_data));
		gtk_widget_queue_draw(gtku_gl->glarea);
	}

	return TRUE;
}


static gboolean
load_op_mat (GtkWidget *widget, gpointer cb_data)
{
	FILE *fp;
	char *filename;
	GTKU_OBJECT *obj;

	/* get Object */
	obj = (GTKU_OBJECT *)g_object_get_data(G_OBJECT(cb_data), "gtku_object");
	g_return_val_if_fail(obj, FALSE);

	/* get File Name */
	filename = (char *)gtk_file_selection_get_filename(
	                   GTK_FILE_SELECTION(cb_data));

	log_printf(5, "Load Operation Matrix File Name : \"%s\"\n", filename);
	log_printf(5, "=== Read Operation Matix Data .... \n");
	GTKU_BOOLEAN_RC_TRY( rc_fopen(filename, "r", &fp) );
	GTKU_BOOLEAN_RC_TRY( read_object_prop_operation_mat(fp, obj->d_prop) );
	GTKU_BOOLEAN_RC_TRY( rc_fclose(fp) );
	log_printf(5, "Done\n");

	return(TRUE);
}


static RC
read_object_prop_operation_mat (FILE *fp, GTKU_DRAW_PROP *dprop)
{
	int ii1;
	IDC idc[4];

	if(dprop == NULL) return(ARG_ERROR_RC);
	if( (dprop->rot_mat   == NULL) ||
	    (dprop->scale_mat == NULL) ||
	    (dprop->shift_mat == NULL) ) return(ARG_ERROR_RC);

	/* Rotation */
	RC_TRY( fgetidc(fp, "c", idc) );
	log_printf(5, "%s\n", idc[0].c);
	for(ii1=0; ii1<4; ii1++){
		RC_TRY( fgetidc(fp, "dddd", idc) );
		log_printf(5, " %15.7e %15.7e %15.7e %15.7e\n",
		           idc[0].d, idc[1].d, idc[2].d, idc[3].d);
		dprop->rot_mat[ii1][0] = idc[0].d;
		dprop->rot_mat[ii1][1] = idc[1].d;
		dprop->rot_mat[ii1][2] = idc[2].d;
		dprop->rot_mat[ii1][3] = idc[3].d;
	}

	/* Shift Matrix */
	RC_TRY( fgetidc(fp, "c", idc) );
	log_printf(5, "%s\n", idc[0].c);
	for(ii1=0; ii1<4; ii1++){
		RC_TRY( fgetidc(fp, "dddd", idc) );
		log_printf(5, " %15.7e %15.7e %15.7e %15.7e\n",
		           idc[0].d, idc[1].d, idc[2].d, idc[3].d);
		dprop->shift_mat[ii1][0] = idc[0].d;
		dprop->shift_mat[ii1][1] = idc[1].d;
		dprop->shift_mat[ii1][2] = idc[2].d;
		dprop->shift_mat[ii1][3] = idc[3].d;
	}

	/* Scale Matrix */
	RC_TRY( fgetidc(fp, "c", idc) );
	log_printf(5, "%s\n", idc[0].c);
	for(ii1=0; ii1<4; ii1++){
		RC_TRY( fgetidc(fp, "dddd", idc) );
		log_printf(5, " %15.7e %15.7e %15.7e %15.7e\n",
		           idc[0].d, idc[1].d, idc[2].d, idc[3].d);
		dprop->scale_mat[ii1][0] = idc[0].d;
		dprop->scale_mat[ii1][1] = idc[1].d;
		dprop->scale_mat[ii1][2] = idc[2].d;
		dprop->scale_mat[ii1][3] = idc[3].d;
	}

	return(NORMAL_RC);
}


GtkWidget *
object_ctl_panel_draw_type (GTKU_OBJECT *obj)
{
	int ii1;
	int type_num = 0;
	GtkWidget *frame;
	char **type_name = NULL;
	char *type_name1[] = {"Default", "Hide", "Bounding Box",
	                      "Polygon", "Solid", "Wire Frame",
	                      "Grid Frame", "Surface Point", "Grid Point",
	                      "Hiddenlin"};
	char *type_name2[] = {"Default", "Hide", "Bounding Box", "Grid Point"};
	char *type_name3[] = {"Default", "Hide", "Bounding Box",
	                      "Texture", "Cell"};
	char *type_name4[] = {"Default", "Hide", "Bounding Box", "Cell"};

	/* Create main Vbox */
	GtkWidget *vbox;
	vbox = gtk_vbox_new(FALSE, 0);
	gtk_widget_show(vbox);

	/* Create frame */
	frame = gtk_frame_new("Drawing Type");
	g_return_val_if_fail(frame, NULL);
	gtk_frame_set_label_align(GTK_FRAME(frame), 0, 0.5);
	gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_ETCHED_IN);
	gtk_container_add(GTK_CONTAINER(frame), GTK_WIDGET(vbox));
	gtk_container_set_border_width(GTK_CONTAINER(frame), 10);
	gtk_widget_show(frame);

	/* Create draw-type check box */
	if(obj->type & GTKU_OBJECT_FEM_ELEM){
		type_name = type_name1;
		type_num = 10;
	}else if( (obj->type & GTKU_OBJECT_FEM_NODE) &&
			 !(obj->type & GTKU_OBJECT_FEM_ELEM) ){ /* node only */
		type_name = type_name2;
		type_num = 4;
	}else if(obj->type & GTKU_OBJECT_IMAGE_PIXEL){
		type_name = type_name3;
		type_num = 5;
	}else if(obj->type & GTKU_OBJECT_IMAGE_VOXEL){
		type_name = type_name4;
		type_num = 4;
	}else{
		GTKU_WARNING("object_ctl_panel_draw_type()",
		             "Unknown object type");
	}

	for(ii1=0; ii1<type_num; ii1++){
		GtkWidget *check = gtk_check_button_new_with_label((type_name)[ii1]);
		g_return_val_if_fail(check, NULL);
		if(strncmp("Default", (type_name)[ii1], 7) == 0){
			gtk_toggle_button_set_state(GTK_TOGGLE_BUTTON(check), TRUE);
		}else{
			gtk_toggle_button_set_state(GTK_TOGGLE_BUTTON(check), FALSE);
		}
		g_object_set_data( G_OBJECT(check), "object_draw_type",
		                   &(obj->d_prop->draw_type) );
		g_signal_connect( G_OBJECT(check), "toggled",
		                  G_CALLBACK(set_local_object_draw_type),
		                  (gpointer)(type_name)[ii1] );
		gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(check), FALSE, FALSE, 1);
		gtk_widget_show(check);
	}

	return(frame);
}


gboolean
set_local_object_draw_type (GtkWidget *widget, gpointer cb_data)
{
	GTKU_DRAW_TYPE  *draw_type;

	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	draw_type = (GTKU_DRAW_TYPE *)g_object_get_data(G_OBJECT(widget),
	                                                "object_draw_type");
	g_return_val_if_fail(draw_type, FALSE);

	if(strncmp("Default", (char *)cb_data, 7) == 0){
		*draw_type ^= GTKU_DRAW_DEFAULT;
	}else if(strncmp("Hide", (char *)cb_data, 4) == 0){
		*draw_type ^= GTKU_DRAW_HIDE;
	}else if(strncmp("Bounding Box", (char *)cb_data, 12) == 0){
		*draw_type ^= GTKU_DRAW_BBOX;
	}else if(strncmp("Polygon", (char *)cb_data, 7) == 0){
		*draw_type ^= GTKU_DRAW_POLYGON;
	}else if(strncmp("Solid", (char *)cb_data, 5) == 0){
		*draw_type ^= GTKU_DRAW_SOLID;
	}else if(strncmp("Wire Frame", (char *)cb_data, 10) == 0){
		*draw_type ^= GTKU_DRAW_WIRE_FRAME;
	}else if(strncmp("Grid Frame", (char *)cb_data, 10) == 0){
		*draw_type ^= GTKU_DRAW_GRID_FRAME;
	}else if(strncmp("Surface Point", (char *)cb_data, 10) == 0){
		*draw_type ^= GTKU_DRAW_SURF_POINT;
	}else if(strncmp("Grid Point", (char *)cb_data, 10) == 0){
		*draw_type ^= GTKU_DRAW_GRID_POINT;
	}else if(strncmp("Hiddenlin", (char *)cb_data, 9) == 0){
		*draw_type ^= GTKU_DRAW_HIDDENLINE;
	}else if(strncmp("Texture", (char *)cb_data, 7) == 0){
		*draw_type ^= GTKU_DRAW_TEXTURE;
	}else if(strncmp("Cell", (char *)cb_data, 4) == 0){
		*draw_type ^= GTKU_DRAW_CELL;
	}else{
		GTKU_WARNING("set_local_object_draw_type()",
		             "Unknown Object Drawing Type!!");
	}

	gtk_widget_queue_draw(gtku_gl->glarea);

	return(TRUE);
}


GtkWidget *
object_ctl_panel_quantity (GTKU_OBJECT *obj)
{
	GtkWidget *vbox;
	GtkWidget *frame;
	GtkObject *adj;
	GtkWidget *scale;

	/* Create main Vbox */
	vbox = gtk_vbox_new(FALSE, 0);
	gtk_widget_show(vbox);

	/* Create frame */
	frame = gtk_frame_new("Quantity of Display");
	g_return_val_if_fail(frame, NULL);
	gtk_frame_set_label_align(GTK_FRAME(frame), 0, 0.5);
	gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_ETCHED_IN);
	gtk_container_add(GTK_CONTAINER(frame), GTK_WIDGET(vbox));
	gtk_container_set_border_width(GTK_CONTAINER(frame), 10);
	gtk_widget_show(frame);

	adj = gtk_adjustment_new(obj->d_prop->quantity, 0.0, 1.0, 0.05, 0.10, 0.0);
	g_return_val_if_fail(adj, NULL);
	scale = gtk_hscale_new( GTK_ADJUSTMENT(adj) );
	g_return_val_if_fail(scale, NULL);
	g_signal_connect(G_OBJECT(adj), "value_changed",
	                 G_CALLBACK(set_local_quantity),
	                 (gpointer)&obj->d_prop->quantity);

	gtk_widget_set_size_request(GTK_WIDGET(scale), 100, -1);
	gtk_range_set_update_policy(GTK_RANGE(scale), GTK_UPDATE_DISCONTINUOUS);
	gtk_scale_set_digits(GTK_SCALE(scale), 2);
	gtk_scale_set_value_pos(GTK_SCALE(scale), GTK_POS_TOP);
	gtk_scale_set_draw_value(GTK_SCALE(scale), TRUE);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(scale), FALSE, FALSE, 0);
	gtk_widget_show(scale);

	return(frame);
}

gboolean
set_local_quantity (GtkWidget *widget, gpointer cb_data)
{
	/* GTKU_GL existence check */
	double *size;
	g_return_val_if_fail(gtku_gl, FALSE);
	size = (double *)cb_data;
	*size = (double)(GTK_ADJUSTMENT(widget))->value;
	gtk_widget_queue_draw(gtku_gl->glarea);
	return(TRUE);
}


GtkWidget *
object_ctl_panel_image_dir (GTKU_OBJECT *obj)
{
	GtkWidget *vbox;
	GtkWidget *frame;

	/* Create main Vbox */
	vbox = gtk_vbox_new(FALSE, 0);
	gtk_widget_show(vbox);

	/* Create frame */
	frame = gtk_frame_new("Image Direction");
	g_return_val_if_fail(frame, NULL);
	gtk_frame_set_label_align(GTK_FRAME(frame), 0, 0.5);
	gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_ETCHED_IN);
	gtk_container_add(GTK_CONTAINER(frame), GTK_WIDGET(vbox));
	gtk_container_set_border_width(GTK_CONTAINER(frame), 10);
	gtk_widget_show(frame);

	/* Create image direction option */
	if( obj->type & (GTKU_OBJECT_IMAGE_PIXEL|GTKU_OBJECT_IMAGE_VOXEL) ){
		GtkWidget *option;
		option = local_image_direction_option(obj);
		g_return_val_if_fail(option, NULL);
		gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(option), FALSE, FALSE, 0);
		gtk_widget_show(option);
	}

	return(frame);
}


GtkWidget *
local_image_direction_option (GTKU_OBJECT *obj)
{
	int ii1;
	GSList *group = NULL;
	GtkWidget *option_menu, *menu;
	char *image_dir[] = {"XY-Plane","YZ-Plane", "ZX-Plane"};

	/* Create opetion menu */
	option_menu = gtk_option_menu_new();
	g_return_val_if_fail(option_menu, NULL);
	menu = gtk_menu_new();
	g_return_val_if_fail(menu, NULL);

	for(ii1=0; ii1<3; ii1++){
		GtkWidget *menuitem;
		menuitem = gtk_radio_menu_item_new_with_label(group, (image_dir)[ii1]);
		group = gtk_radio_menu_item_group(GTK_RADIO_MENU_ITEM(menuitem));
		gtk_menu_append(GTK_MENU(menu), menuitem);
		g_object_set_data(G_OBJECT(menuitem), "gtku_object", obj);
		g_signal_connect( G_OBJECT(menuitem), "activate",
		                  G_CALLBACK(set_local_image_direction),
		                  (gpointer)(image_dir)[ii1] );
		gtk_widget_show(menuitem);
	}

	gtk_option_menu_set_menu(GTK_OPTION_MENU(option_menu), menu);

	return(option_menu);
}


gboolean
set_local_image_direction (GtkWidget *widget, gpointer cb_data)
{
	GTKU_OBJECT *obj;

	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	obj = (GTKU_OBJECT *)g_object_get_data(G_OBJECT(widget), "gtku_object");
	g_return_val_if_fail(obj, FALSE);

	if(strncmp("XY-Plane", (char *)cb_data, 8) == 0){
		obj->d_prop->image_dir = GTKU_IMAGE_XY;
	}else if(strncmp("YZ-Plane", (char *)cb_data, 8) == 0){
		obj->d_prop->image_dir = GTKU_IMAGE_YZ;
	}else if(strncmp("ZX-Plane", (char *)cb_data, 8) == 0){
		obj->d_prop->image_dir = GTKU_IMAGE_ZX;
	}else{
		GTKU_WARNING("set_image_direction()",
		             "Unknown Image Object Direction!!");
	}
	obj->d_prop->init_image_dir = FALSE;
	gtk_widget_queue_draw(gtku_gl->glarea);

	return(TRUE);
}


GtkWidget *
object_ctl_panel_image_contour (GTKU_OBJECT *obj)
{
	GtkWidget *vbox;
	GtkWidget *frame;

	/* Create main Vbox */
	vbox = gtk_vbox_new(FALSE, 0);
	gtk_widget_show(vbox);

	/* Create frame */
	frame = gtk_frame_new("Contour Map");
	g_return_val_if_fail(frame, NULL);
	gtk_frame_set_label_align(GTK_FRAME(frame), 0, 0.5);
	gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_ETCHED_IN);
	gtk_container_add(GTK_CONTAINER(frame), GTK_WIDGET(vbox));
	gtk_container_set_border_width(GTK_CONTAINER(frame), 10);
	gtk_widget_show(frame);

	/* Create image direction option */
	if( obj->type & (GTKU_OBJECT_IMAGE_PIXEL|GTKU_OBJECT_IMAGE_VOXEL) ){
		GtkWidget *option;
		option = local_image_contour_option(obj);
		g_return_val_if_fail(option, NULL);
		gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(option), FALSE, FALSE, 0);
		gtk_widget_show(option);
	}

	return(frame);
}


GtkWidget *
local_image_contour_option (GTKU_OBJECT *obj)
{
	int ii1;
	GSList *group = NULL;
	GtkWidget *option_menu, *menu;
	char *contour_map[] = {"BlueGreenRed","BlueRed", "Hot",
	                       "Cold", "GrayScale"};

	/* Create opetion menu */
	option_menu = gtk_option_menu_new();
	g_return_val_if_fail(option_menu, NULL);
	menu = gtk_menu_new();
	g_return_val_if_fail(menu, NULL);

	for(ii1=0; ii1<5; ii1++){
		GtkWidget *menuitem;
		menuitem = gtk_radio_menu_item_new_with_label(group,
		                                              (contour_map)[ii1]);
		group = gtk_radio_menu_item_group(GTK_RADIO_MENU_ITEM(menuitem));
		gtk_menu_append(GTK_MENU(menu), menuitem);
		g_object_set_data(G_OBJECT(menuitem), "gtku_object", obj);
		g_signal_connect( G_OBJECT(menuitem), "activate",
		                  G_CALLBACK(set_local_image_contour),
		                  (gpointer)(contour_map)[ii1] );
		gtk_widget_show(menuitem);
	}

	gtk_option_menu_set_menu(GTK_OPTION_MENU(option_menu), menu);

	return(option_menu);
}


gboolean
set_local_image_contour (GtkWidget *widget, gpointer cb_data)
{
	GTKU_OBJECT *obj;

	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	obj = (GTKU_OBJECT *)g_object_get_data(G_OBJECT(widget), "gtku_object");
	g_return_val_if_fail(obj, FALSE);
	
	if(strncmp("BlueGreenRed", (char *)cb_data, 12) == 0){
		obj->d_prop->cont_type = GRAPHIC_CONT_BGR24;
	}else if(strncmp("BlueRed", (char *)cb_data, 7) == 0){
		obj->d_prop->cont_type = GRAPHIC_CONT_BR24;
	}else if(strncmp("Hot", (char *)cb_data, 3) == 0){
		obj->d_prop->cont_type = GRAPHIC_CONT_HOT24;
	}else if(strncmp("Cold", (char *)cb_data, 4) == 0){
		obj->d_prop->cont_type = GRAPHIC_CONT_COLD24;
	}else if(strncmp("GrayScale", (char *)cb_data, 9) == 0){
		obj->d_prop->cont_type = GRAPHIC_CONT_GRAYSCALE24;
	}else{
		GTKU_WARNING("set_image_direction()",
		             "Unknown Contour Type!!");
	}
	gtk_widget_queue_draw(gtku_gl->glarea);

	return(TRUE);
}


GtkWidget *
object_ctl_panel_image_range (GTKU_OBJECT *obj)
{
	GtkWidget *vbox;
	GtkWidget *frame;

	/* Create main Vbox */
	vbox = gtk_vbox_new(FALSE, 0);
	gtk_widget_show(vbox);

	/* Create frame */
	frame = gtk_frame_new("Display Range");
	g_return_val_if_fail(frame, NULL);
	gtk_frame_set_label_align(GTK_FRAME(frame), 0, 0.5);
	gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_ETCHED_IN);
	gtk_container_add(GTK_CONTAINER(frame), GTK_WIDGET(vbox));
	gtk_container_set_border_width(GTK_CONTAINER(frame), 10);
	gtk_widget_show(frame);

	{ /* Create max range  */
	 GtkWidget *scale = local_image_range_with_label(
	                       obj, &(obj->d_prop->img_threshold_max), "MAX");
	 g_return_val_if_fail(scale, NULL);
	 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(scale), FALSE, FALSE, 5);
	 gtk_widget_show(scale);
	}

	{ /* Create min range  */
	 GtkWidget *scale = local_image_range_with_label(
	                       obj, &(obj->d_prop->img_threshold_min), "MIN");
	 g_return_val_if_fail(scale, NULL);
	 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(scale), FALSE, FALSE, 5);
	 gtk_widget_show(scale);
	}

	return(frame);
}


static GtkWidget *
local_image_range_with_label (GTKU_OBJECT *obj, double *range, char *label)
{
	GtkObject *adj = gtk_adjustment_new(*range, 0.0, 1.0, 0.01, 0.05, 0.0);
	GtkWidget *scale = gtk_hscale_new( GTK_ADJUSTMENT(adj) );
	GtkWidget *hbox = gtk_hbox_new(FALSE, 0);
	GtkWidget *name = gtk_label_new(label);
	g_return_val_if_fail(adj, NULL);
	g_return_val_if_fail(scale, NULL);
	g_return_val_if_fail(hbox, NULL);
	g_return_val_if_fail(name, NULL);

	g_signal_connect(G_OBJECT(adj), "value_changed",
	                 G_CALLBACK(set_local_image_range),(gpointer)(range));

	gtk_widget_set_size_request(GTK_WIDGET(scale), 100, -1);
	gtk_range_set_update_policy(GTK_RANGE(scale), GTK_UPDATE_DISCONTINUOUS);
	gtk_scale_set_digits(GTK_SCALE(scale), 2);
	gtk_scale_set_value_pos(GTK_SCALE(scale), GTK_POS_TOP);
	gtk_scale_set_draw_value(GTK_SCALE(scale), TRUE);

	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(name),  TRUE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(scale), TRUE, TRUE, 0);
	gtk_widget_show(name);
	gtk_widget_show(scale);

	return(hbox);
}


gboolean
set_local_image_range (GtkWidget *widget, gpointer cb_data)
{
	/* GTKU_GL existence check */
	double *size;
	g_return_val_if_fail(gtku_gl, FALSE);
	size = (double *)cb_data;
	*size = (double)(GTK_ADJUSTMENT(widget))->value;
	gtk_widget_queue_draw(gtku_gl->glarea);
	return(TRUE);
}

GtkWidget *
object_ctl_panel_image_trans (GTKU_OBJECT *obj)
{
	GtkWidget *vbox;
	GtkWidget *frame;
	GtkObject *adj;
	GtkWidget *scale;

	/* Create main Vbox */
	vbox = gtk_vbox_new(FALSE, 0);
	gtk_widget_show(vbox);

	/* Create frame */
	frame = gtk_frame_new("TransLucent of Image");
	g_return_val_if_fail(frame, NULL);
	gtk_frame_set_label_align(GTK_FRAME(frame), 0, 0.5);
	gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_ETCHED_IN);
	gtk_container_add(GTK_CONTAINER(frame), GTK_WIDGET(vbox));
	gtk_container_set_border_width(GTK_CONTAINER(frame), 10);
	gtk_widget_show(frame);

	adj = gtk_adjustment_new(obj->d_prop->img_alpha, 0.0, 1.0, 0.05, 0.10, 0.0);
	g_return_val_if_fail(adj, NULL);
	scale = gtk_hscale_new( GTK_ADJUSTMENT(adj) );
	g_return_val_if_fail(scale, NULL);
	g_signal_connect(G_OBJECT(adj), "value_changed",
	                 G_CALLBACK(set_local_image_trans),
	                 (gpointer)&obj->d_prop->img_alpha);

	gtk_widget_set_size_request(GTK_WIDGET(scale), 100, -1);
	gtk_range_set_update_policy(GTK_RANGE(scale), GTK_UPDATE_DISCONTINUOUS);
	gtk_scale_set_digits(GTK_SCALE(scale), 2);
	gtk_scale_set_value_pos(GTK_SCALE(scale), GTK_POS_TOP);
	gtk_scale_set_draw_value(GTK_SCALE(scale), TRUE);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(scale), FALSE, FALSE, 0);
	gtk_widget_show(scale);

	return(frame);
}


gboolean
set_local_image_trans (GtkWidget *widget, gpointer cb_data)
{
	double *size;

	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);
	size = (double *)cb_data;
	*size = (double)(GTK_ADJUSTMENT(widget))->value;
	gtk_widget_queue_draw(gtku_gl->glarea);
	return(TRUE);
}


GtkWidget *
object_ctl_panel_color (GTKU_OBJECT *obj)
{
	GtkWidget *vbox;
	GtkWidget *frame;

	/* Create main Vbox */
	vbox = gtk_vbox_new(FALSE, 0);
	gtk_widget_show(vbox);

	/* Create frame */
	frame = gtk_frame_new("Color & Alpha");
	g_return_val_if_fail(frame, NULL);
	gtk_frame_set_label_align(GTK_FRAME(frame), 0, 0.5);
	gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_ETCHED_IN);
	gtk_container_add(GTK_CONTAINER(frame), GTK_WIDGET(vbox));
	gtk_container_set_border_width(GTK_CONTAINER(frame), 10);
	gtk_widget_show(frame);

	{ /* Create transration check button  */
	 GtkWidget *trans_frame = object_ctl_panel_trans_check(obj);
	 g_return_val_if_fail(trans_frame, NULL);
	 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(trans_frame), FALSE,FALSE, 1);
	 gtk_widget_show(trans_frame);
	}

	{ /* Create color selection button  */
	 GtkWidget *colorsel_frame = object_ctl_panel_color_button(obj);
	 g_return_val_if_fail(colorsel_frame, NULL);
	 gtk_box_pack_start(GTK_BOX(vbox),GTK_WIDGET(colorsel_frame),FALSE,FALSE,1);
	 gtk_widget_show(colorsel_frame);
	}

	return(frame);
}


GtkWidget *
object_ctl_panel_trans_check (GTKU_OBJECT *obj)
{
	int ii1;
	char *name[] = {"Point", "Line", "Polygon", "Text"};
	gboolean *sw[] = {&(obj->d_prop->point_trans_sw),
	                  &(obj->d_prop->line_trans_sw),
	                  &(obj->d_prop->polygon_trans_sw),
	                  &(obj->d_prop->text_trans_sw)};
	GtkWidget *vbox;
	GtkWidget *frame;

	/* Create main Vbox */
	vbox = gtk_vbox_new(FALSE, 0);
	gtk_widget_show(vbox);

	/* Create frame */
	frame = gtk_frame_new("TransLucent");
	g_return_val_if_fail(frame, NULL);
	gtk_frame_set_label_align(GTK_FRAME(frame), 0, 0.5);
	gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_ETCHED_IN);
	gtk_container_add(GTK_CONTAINER(frame), GTK_WIDGET(vbox));
	gtk_container_set_border_width(GTK_CONTAINER(frame), 10);

	/* Create transration check button */
	for(ii1=0; ii1<4; ii1++){
		GtkWidget *check = gtk_check_button_new_with_label(name[ii1]);
		g_return_val_if_fail(check, NULL);
		gtk_toggle_button_set_state(GTK_TOGGLE_BUTTON(check), FALSE);
		g_signal_connect(G_OBJECT(check), "toggled",
		                 G_CALLBACK(set_boolean_switch),
		                 (gpointer)sw[ii1]);
		gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(check), FALSE, FALSE, 1);
		gtk_widget_show(check);
	}

	return(frame);
}


GtkWidget *
object_ctl_panel_color_button (GTKU_OBJECT *obj)
{
	int ii1;
	char *name[] = {"Point", "Line", "Polygon", "Text"};
	GLdouble (*(color[]))[4] = {&(obj->d_prop->point_color),
	                           &(obj->d_prop->line_color),
	                           &(obj->d_prop->polygon_color),
	                           &(obj->d_prop->text_color)};
	GtkWidget *vbox;
	GtkWidget *frame;

	/* Create main Vbox */
	vbox = gtk_vbox_new(FALSE, 0);
	gtk_widget_show(vbox);

	/* Create frame */
	frame = gtk_frame_new("Color & Transparency");
	g_return_val_if_fail(frame, NULL);
	gtk_frame_set_label_align(GTK_FRAME(frame), 0, 0.5);
	gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_ETCHED_IN);
	gtk_container_add(GTK_CONTAINER(frame), GTK_WIDGET(vbox));
	gtk_container_set_border_width(GTK_CONTAINER(frame), 10);

	/* Create color selection dialog button */
	for(ii1=0; ii1<4; ii1++){
		GtkWidget *button = gtk_button_new_with_label(name[ii1]);
		gtk_widget_set_size_request(GTK_WIDGET(button), 20, -1);
		g_object_set_data(G_OBJECT(button), "object_d_prop_color",
		                  color[ii1]);
		g_signal_connect(G_OBJECT(button), "clicked",
		                 G_CALLBACK(local_object_color_selection_dialog),
		                 (gpointer)name[ii1]);
		gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(button), FALSE, FALSE, 0);
		gtk_widget_show(button);
	}

	return(frame);
}


gboolean
local_object_color_selection_dialog (GtkWidget *widget, gpointer cb_data)
{
	GLdouble (*color)[4];
	GtkWidget *cs_widget;
	GtkColorSelectionDialog *dialog;
	GtkColorSelection *colorsel;
	GdkColor current_color;
	guint16 alpha;

	color = (GLdouble (*)[4])g_object_get_data(G_OBJECT(widget),
	                                          "object_d_prop_color");

	cs_widget = gtk_color_selection_dialog_new((char *)cb_data);
	g_return_val_if_fail(cs_widget, FALSE);

	dialog = GTK_COLOR_SELECTION_DIALOG(cs_widget);
	g_return_val_if_fail(dialog, FALSE);

	colorsel = GTK_COLOR_SELECTION(dialog->colorsel);
	g_return_val_if_fail(colorsel, FALSE);

	current_color.red   = (guint16)((*color)[0] * 0xffff);
	current_color.green = (guint16)((*color)[1] * 0xffff);
	current_color.blue  = (guint16)((*color)[2] * 0xffff);
	              alpha = (guint16)((*color)[3] * 0xffff);

	gtk_color_selection_set_has_opacity_control(colorsel, TRUE);
	gtk_color_selection_set_update_policy(colorsel, GTK_UPDATE_DELAYED);
	gtk_color_selection_set_previous_color(colorsel, &current_color);
	gtk_color_selection_set_current_color(colorsel, &current_color);
	gtk_color_selection_set_previous_alpha(colorsel, alpha);
	gtk_color_selection_set_current_alpha(colorsel, alpha);
	gtk_color_selection_set_has_palette(colorsel, TRUE);

	/* Select */
	g_object_set_data(G_OBJECT(dialog->colorsel),
	                  "object_d_prop_color", color);
	g_signal_connect(G_OBJECT(dialog->colorsel), "color_changed",
	                 G_CALLBACK(set_local_object_color),
	                 (gpointer)colorsel);
	/* OK */
	g_object_set_data(G_OBJECT(dialog->ok_button),
	                  "object_d_prop_color", color);
	g_signal_connect(G_OBJECT(dialog->ok_button), "clicked",
	                 G_CALLBACK(set_local_object_color),
	                 (gpointer)colorsel);
	g_signal_connect_swapped(G_OBJECT(dialog->ok_button), "clicked",
	                         G_CALLBACK(gtk_widget_destroy), (gpointer)dialog);
	/* Cansel */
	g_object_set_data(G_OBJECT(dialog->cancel_button),
	                  "object_d_prop_color", color);
	g_object_set_data(G_OBJECT(dialog->cancel_button),
	                  "dialog_window", cs_widget);
	g_signal_connect(G_OBJECT(dialog->cancel_button), "clicked",
	                 G_CALLBACK(unset_local_object_color),
	                 G_OBJECT(colorsel));

	gtk_widget_show(cs_widget);
	gtk_main();

	return(TRUE);
}


gboolean
set_local_object_color (GtkWidget *widget, gpointer cb_data)
{
	GLdouble (*glcolor)[4];
	GtkColorSelection *colorsel;
	GdkColor gdkcolor;

	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	glcolor = (GLdouble (*)[4])g_object_get_data(G_OBJECT(widget),
	                                            "object_d_prop_color");
	g_return_val_if_fail(glcolor, FALSE);
	colorsel = (GtkColorSelection *)cb_data;
	gtk_color_selection_get_current_color(colorsel, &gdkcolor);
	(*glcolor)[0] = (gdkcolor.red)  /(GLdouble)0xffff;
	(*glcolor)[1] = (gdkcolor.green)/(GLdouble)0xffff;
	(*glcolor)[2] = (gdkcolor.blue) /(GLdouble)0xffff;
	(*glcolor)[3]
	    = gtk_color_selection_get_current_alpha(colorsel) / (GLdouble)0xffff;

	gtk_widget_queue_draw(gtku_gl->glarea);

	return(TRUE);
}


gboolean
unset_local_object_color (GtkWidget *widget, gpointer cb_data)
{
	GLdouble (*glcolor)[4];
	GtkColorSelection *colorsel;
	GdkColor gdkcolor;
	GtkWidget *dialog_window;

	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	glcolor = (GLdouble (*)[4])g_object_get_data(G_OBJECT(widget),
	                                     "object_d_prop_color");
	dialog_window = (GtkWidget *)g_object_get_data(G_OBJECT(widget),
	                                               "dialog_window");
	g_return_val_if_fail(glcolor, FALSE);
	g_return_val_if_fail(dialog_window, FALSE);
	colorsel = (GtkColorSelection *)cb_data;
	gtk_color_selection_get_previous_color(colorsel, &gdkcolor);
	(*glcolor)[0] = (gdkcolor.red)  /(GLdouble)0xffff;
	(*glcolor)[1] = (gdkcolor.green)/(GLdouble)0xffff;
	(*glcolor)[2] = (gdkcolor.blue) /(GLdouble)0xffff;
	(*glcolor)[3]
	      = gtk_color_selection_get_previous_alpha(colorsel) / (GLdouble)0xffff;

	gtk_widget_queue_draw(gtku_gl->glarea);
	gtk_widget_destroy(dialog_window);

	return(TRUE);
}


GtkWidget *
object_ctl_panel_draw_num (GTKU_OBJECT *obj)
{
	GtkWidget *vbox;
	GtkWidget *frame;

	/* Create main Vbox */
	vbox = gtk_vbox_new(FALSE, 0);
	gtk_widget_show(vbox);

	/* Create frame */
	frame = gtk_frame_new("Drawing Number");
	g_return_val_if_fail(frame, NULL);
	gtk_frame_set_label_align(GTK_FRAME(frame), 0, 0.5);
	gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_ETCHED_IN);
	gtk_container_add(GTK_CONTAINER(frame), GTK_WIDGET(vbox));
	gtk_container_set_border_width(GTK_CONTAINER(frame), 10);
	gtk_widget_show(frame);

	{ /* Create check button Node Number */
	 GtkWidget *check = gtk_check_button_new_with_label("Node");
	 g_return_val_if_fail(check, NULL);
	 gtk_toggle_button_set_state(GTK_TOGGLE_BUTTON(check), FALSE);
	 g_signal_connect(G_OBJECT(check), "toggled",
	                  G_CALLBACK(set_boolean_switch),
	                  (gpointer)&(obj->d_prop->draw_node_num_sw));
	 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(check), FALSE, FALSE, 1);
	 gtk_widget_show(check);
	}

	{ /* Create check button Element Number */
	 GtkWidget *check = gtk_check_button_new_with_label("Element");
	 g_return_val_if_fail(check, NULL);
	 gtk_toggle_button_set_state(GTK_TOGGLE_BUTTON(check), FALSE);
	 g_signal_connect(G_OBJECT(check), "toggled",
	                  G_CALLBACK(set_boolean_switch),
	                  (gpointer)&(obj->d_prop->draw_elem_num_sw));
	 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(check), FALSE, FALSE, 1);
	 gtk_widget_show(check);
	}

	{ /* Create check button Surface Number */
	 GtkWidget *check = gtk_check_button_new_with_label("Surface");
	 g_return_val_if_fail(check, NULL);
	 gtk_toggle_button_set_state(GTK_TOGGLE_BUTTON(check), FALSE);
	 g_signal_connect(G_OBJECT(check), "toggled",
	                  G_CALLBACK(set_boolean_switch),
	                  (gpointer)&(obj->d_prop->draw_surf_num_sw));
	 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(check), FALSE, FALSE, 1);
	 gtk_widget_show(check);
	}

	return(frame);
}


GtkWidget *
object_ctl_panel_size (GTKU_OBJECT *obj)
{
	GtkWidget *vbox;
	GtkWidget *frame;
	GtkWidget *hbox;

	/* Create main Vbox */
	vbox = gtk_vbox_new(FALSE, 0);
	gtk_widget_show(vbox);

	/* Create frame */
	frame = gtk_frame_new("Line & Point Size");
	g_return_val_if_fail(frame, NULL);
	gtk_frame_set_label_align(GTK_FRAME(frame), 0, 0.5);
	gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_ETCHED_IN);
	gtk_container_add(GTK_CONTAINER(frame), GTK_WIDGET(vbox));
	gtk_container_set_border_width(GTK_CONTAINER(frame), 10);
	gtk_widget_show(frame);

	/* Create main Hbox */
	hbox = gtk_hbox_new(FALSE, 0);
	gtk_widget_show(hbox);

	{ /* Create Point Size Spin Button */
	 GtkWidget *spin = local_object_size_spin_with_label(
	                       obj,&(obj->d_prop->point_size),"Point:");
	 g_return_val_if_fail(spin, NULL);
	 gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(spin), FALSE, FALSE, 5);
	 gtk_widget_show(spin);
	}

	{ /* Create Line Size Spin Button */
	 GtkWidget *spin = local_object_size_spin_with_label(
	                       obj,&(obj->d_prop->line_size),"Line:");
	 g_return_val_if_fail(spin, NULL);
	 gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(spin), FALSE, FALSE, 5);
	 gtk_widget_show(spin);
	}

	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 5);

	return(frame);
}


static GtkWidget *
local_object_size_spin_with_label (GTKU_OBJECT *obj, GLdouble *size,
                                        char *label)
{
	GtkObject *adjustment = gtk_adjustment_new(*size, 0.1, 20.0, 0.1, 0.5, 0.0);
	GtkWidget *spin = gtk_spin_button_new( GTK_ADJUSTMENT(adjustment), 0.1, 1);
	GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
	GtkWidget *name = gtk_label_new(label);
	g_return_val_if_fail(adjustment, NULL);
	g_return_val_if_fail(spin, NULL);
	g_return_val_if_fail(vbox, NULL);
	g_return_val_if_fail(name, NULL);

	g_signal_connect(G_OBJECT(adjustment), "value_changed",
	                 G_CALLBACK(set_local_object_size),(gpointer)(size));

	gtk_misc_set_alignment(GTK_MISC(name), 0, 0.5);
	gtk_spin_button_set_wrap(GTK_SPIN_BUTTON(spin), FALSE);
	gtk_spin_button_set_update_policy(GTK_SPIN_BUTTON(spin),
	                                  GTK_UPDATE_DISCONTINUOUS);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(name), FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(spin), FALSE, FALSE, 0);
	gtk_widget_show(name);
	gtk_widget_show(spin);

	return(vbox);
}


gboolean
set_local_object_size (GtkWidget *widget, gpointer cb_data)
{
	/* GTKU_GL existence check */
	GLdouble *size;
	g_return_val_if_fail(gtku_gl, FALSE);
	size = (GLdouble *)cb_data;
	*size = (GLdouble)(GTK_ADJUSTMENT(widget))->value;
	gtk_widget_queue_draw(gtku_gl->glarea);
	return(TRUE);
}


gboolean
set_boolean_switch (GtkWidget *widget, gpointer cb_data)
{
	/* GTKU_GL existence check */
	gboolean *sw;
	g_return_val_if_fail(gtku_gl, FALSE);
	sw = (gboolean *)cb_data;
	*sw = !(*sw);
	gtk_widget_queue_draw(gtku_gl->glarea);
	return(TRUE);
}


