/*********************************************************************
 * gtk_model_view.c
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

/* $Id: gtk_model_view.c 685 2005-12-16 09:59:51Z sasaoka $ */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <glib.h>
#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>
#include <gtk/gtkgl.h>
#include "gtk_utl.h"
#include "gtk_icons.h"
#include "rc.h"
#include "log_printf.h"
#include "memory_manager.h"
#include "math_utl.h"
#include "fem_struct.h"
#include "nst_component.h"
#include "eps_utl.h"
#include "macros.h"


#define INIT_WINDOW_SIZE_X (400)
#define INiT_WINDOW_SIZE_Y (400)

static GLdouble color[][4] = {
	{1.0, 1.0, 1.0, 0.5}, /* White */
	{1.0, 0.0, 0.0, 0.5}, /* Red */
	{0.0, 1.0, 0.0, 0.5}, /* Green */
	{0.0, 0.0, 1.0, 0.5}, /* Blue */
	{1.0, 1.0, 0.0, 0.5}, /* Yellow */
	{0.0, 1.0, 1.0, 0.5}, /* Cyan */
	{1.0, 0.0, 1.0, 0.5}, /* Magenta */
	{0.0, 0.0, 0.0, 0.5}  /* Black */
};

/* extern global variable                                           */
/*  There's no method except using gloval variable for signal-call, */
/*  toggle switch, and so on.                                       */
GTKU_GL *gtku_gl = NULL;

/* local global variable */
static GtkWidget *draw_type_radio[7];
static GtkItemFactory *item_factory = NULL;

/* Model view Window */
static GtkWidget *create_glarea(GtkWidget *window);
static GtkWidget *create_toolbar1(GtkWidget *window);
static GtkWidget *create_toolbar2(GtkWidget *window);
static RC create_ctl_tab(GtkWidget *notebook, GTKU_OBJECT *obj);
static gboolean file_selection_bdf(GtkWidget *widget, gpointer cb_data);
static gboolean filew_ok_bdf(GtkWidget *widget, gpointer cb_data);
static gboolean file_selection_eps(GtkWidget *widget, gpointer cb_data);
static gboolean filew_ok_eps(GtkWidget *widget, gpointer cb_data);
static gboolean save_bdf_file(GtkWidget *widget, gpointer cb_data);
static gboolean save_eps_file(GtkWidget *widget, gpointer cb_data);

/* CTRL Panel */
static gboolean init_op_mat(GtkWidget *widget, gpointer cb_data);
static gboolean file_selection_save_op_mat(GtkWidget *widget,
                                           gpointer cb_data);
static gboolean filew_ok_save_op_mat(GtkWidget *widget, gpointer cb_data);
static gboolean save_op_mat(GtkWidget *widget, gpointer cb_data);
static gboolean file_selection_load_op_mat(GtkWidget *widget,
                                           gpointer cb_data);
static gboolean filew_ok_load_op_mat(GtkWidget *widget, gpointer cb_data);
static gboolean load_op_mat(GtkWidget *widget, gpointer cb_data);

/* etc */
static gboolean format_chk(const char *format);
static RC set_element_dim(ELEMENT_ARRAY *elem, GTKU_OBJECT_DIM *dim);
static int drawtype2number(GTKU_DRAW_TYPE type);
static void gtku_set_model_scale_for_IF(
	gpointer cb_data, guint cb_action, GtkWidget *widget);
static void gtku_set_model_draw_type_for_IF(
	gpointer cb_data, guint cb_action, GtkWidget *widget);
static void gtku_set_model_view_direction_for_IF(
	gpointer cb_data, guint cb_action, GtkWidget *widget);
static void gtku_set_world_plane_for_IF(
	gpointer cb_data, guint cb_action, GtkWidget *widget);
static GtkWidget * gtku_local_image_range_with_label(
	GTKU_OBJECT *obj, double *range, char *label);
static GtkWidget *gtku_local_object_size_spin_with_label(
	GTKU_OBJECT *obj, GLdouble *size, char *label);


RC
gtk_model_view (void)
{
	int ii1;
	GtkWidget *window;
	GtkWidget *vbox;
	GtkWidget *tool_hbox;
	GtkWidget *hpaned;
	GtkWidget *glarea;
	GtkWidget *notebook;

	/* Create new top level Window. */
	window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
	gtk_window_set_title(GTK_WINDOW(window), "gtk_fem");
	gtk_container_set_border_width(GTK_CONTAINER(window), 0);
	gtk_window_set_default_size(GTK_WINDOW(window), 700, 600);
	gtk_quit_add_destroy(1, GTK_OBJECT(window));

	g_signal_connect(G_OBJECT(window), "delete_event",
	                 G_CALLBACK(gtku_window_destory_ask), window);
	g_signal_connect(G_OBJECT(window), "destroy",
	                 G_CALLBACK(gtk_main_quit), NULL);
	gtk_widget_show(GTK_WIDGET(window));


	/* Create main Vbox */
	vbox = gtk_vbox_new(FALSE, 2);
	gtk_container_add(GTK_CONTAINER(window), vbox);
	gtk_widget_show(vbox);

	{ /* Create Menubar */
	  GtkWidget *menubar;
	  RC_NULL_CHK( menubar = gtku_create_model_view_menu(window) );
	  gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(menubar), FALSE, FALSE, 0);
	  gtk_widget_show(menubar);
	  /*メーニューの追加の仕方 
	  GtkItemFactoryEntry item = {"/Option/hoge", NULL, NULL, 0, NULL};
	  gtk_item_factory_create_item(gtk_item_factory_from_widget(menubar),
	                               &item, NULL, 1);
	  */
	}
	
	/* Create Toolbar HBox */
	tool_hbox = gtk_hbox_new(FALSE, 0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(tool_hbox), FALSE, FALSE, 0);
	gtk_widget_show(tool_hbox);

	{/* Create Toolbar : Zoom */
	  GtkWidget *toolbar;
	  RC_NULL_CHK( toolbar = create_toolbar1(window) );
	  gtk_box_pack_start(GTK_BOX(tool_hbox),GTK_WIDGET(toolbar),FALSE,FALSE,0);
	  gtk_widget_show(toolbar);
	}

	{/* Create Toolbar 2 */
	  GtkWidget *toolbar;
	  RC_NULL_CHK( toolbar = create_toolbar2(window) );
	  gtk_box_pack_start(GTK_BOX(tool_hbox),GTK_WIDGET(toolbar),FALSE,FALSE,0);
	  gtk_widget_show(toolbar);
	}

	/* Create H-Paned Widget */
	RC_NULL_CHK( hpaned = gtk_hpaned_new() );
	gtk_container_add(GTK_CONTAINER(vbox), hpaned);
	gtk_widget_show(hpaned);

	/* Add H-Paned left widget */
	RC_NULL_CHK( glarea = create_glarea(window) );
	gtk_paned_pack1(GTK_PANED(hpaned), GTK_WIDGET(glarea), FALSE, FALSE);
	gtk_widget_show(glarea);


	/* Add H-Paned right widget */
	RC_NULL_CHK( notebook = gtk_notebook_new() );
	gtk_notebook_set_tab_pos(GTK_NOTEBOOK(notebook), GTK_POS_TOP);
	gtk_widget_set_size_request(notebook, 170, -1);
	/*
	gtk_widget_set_size_request(notebook, 0, -1);
	*/
	gtk_paned_pack2(GTK_PANED(hpaned), GTK_WIDGET(notebook), FALSE, TRUE);
	gtk_widget_show(notebook);

	/* Add NoteBook Tab */
	for(ii1=0; ii1<(gtku_gl->object->size); ii1++){
		RC_TRY( create_ctl_tab(notebook, gtku_gl->object->obj[ii1]) );
	}

	return(NORMAL_RC);
}


static GtkWidget *
create_toolbar1 (GtkWidget *window)
{
	/* Create Handlebox */
	GtkWidget *handlebox = gtk_handle_box_new();
	g_return_val_if_fail(handlebox, NULL);
	gtk_widget_show(handlebox);

	/* Create Toolbar */
	GtkWidget *toolbar = gtk_toolbar_new();
	g_return_val_if_fail(toolbar, NULL);
	gtk_toolbar_set_orientation(GTK_TOOLBAR(toolbar),
	                            GTK_ORIENTATION_HORIZONTAL);
	gtk_toolbar_set_style(GTK_TOOLBAR(toolbar), GTK_TOOLBAR_BOTH);
	gtk_container_add(GTK_CONTAINER(handlebox), toolbar);
	gtk_widget_show(toolbar);

	/* Clear */
	g_return_val_if_fail(gtku_toolbar_clear(window, toolbar), NULL);
	gtk_toolbar_append_space(GTK_TOOLBAR(toolbar));

	/* Scale */
	g_return_val_if_fail(gtku_toolbar_scale(window, toolbar), NULL);
	gtk_toolbar_append_space(GTK_TOOLBAR(toolbar));

	/* DrawType */
	g_return_val_if_fail(gtku_toolbar_drawtype(window, toolbar), NULL);
	gtk_toolbar_append_space(GTK_TOOLBAR(toolbar));

	/* Animation */
	g_return_val_if_fail(gtku_toolbar_animation(window, toolbar), NULL);
	/*
	gtk_toolbar_append_space(GTK_TOOLBAR(toolbar));
	*/

	return(handlebox);
}


static GtkWidget *
create_toolbar2 (GtkWidget *window)
{
	/* Create Handlebox */
	GtkWidget *handlebox = gtk_handle_box_new();
	g_return_val_if_fail(handlebox, NULL);
	gtk_widget_show(handlebox);

	/* Create Toolbar */
	GtkWidget *toolbar = gtk_toolbar_new();
	g_return_val_if_fail(toolbar, NULL);
	gtk_toolbar_set_orientation(GTK_TOOLBAR(toolbar),
	                            GTK_ORIENTATION_HORIZONTAL);
	gtk_toolbar_set_style(GTK_TOOLBAR(toolbar), GTK_TOOLBAR_BOTH);
	gtk_container_add(GTK_CONTAINER(handlebox), toolbar);
	gtk_widget_show(toolbar);

	/* Rotation */
	g_return_val_if_fail(gtku_toolbar_rotation(window, toolbar), NULL);
	/*
	gtk_toolbar_append_space(GTK_TOOLBAR(toolbar));
	*/

	return(handlebox);
}


static GtkWidget *
create_glarea (GtkWidget *window)
{	
	GtkWidget *frame;
	GtkWidget *gldrawarea;

	/* Create new FRAME */
	frame = gtk_frame_new("VIEW AREA");
	g_return_val_if_fail(frame, NULL);
	gtk_frame_set_label_align(GTK_FRAME(frame), 0.5, 0.5);
	gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_ETCHED_OUT);
	gtk_widget_show(frame);

	/* Create new GLAREA */
	gldrawarea = gtku_glarea_new(window);
	g_return_val_if_fail(gldrawarea, NULL);
	gtk_container_add(GTK_CONTAINER(frame), gldrawarea);
	gtk_widget_show(GTK_WIDGET(gldrawarea));

	return(frame);
}


static RC
create_ctl_tab (GtkWidget *notebook, GTKU_OBJECT *obj)
{
	GtkWidget *vbox;

	RC_NULL_CHK(obj);
	RC_NULL_CHK(notebook);


	/* Create Scrolled Window */
	GtkWidget *scr_window = gtk_scrolled_window_new(NULL, NULL);
	RC_NULL_CHK(scr_window);
	gtk_container_set_border_width(GTK_CONTAINER(scr_window), 0);
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scr_window),
	                               GTK_POLICY_AUTOMATIC,
	                               GTK_POLICY_AUTOMATIC);
	gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW(scr_window),
	                                    GTK_SHADOW_NONE);
	gtk_widget_show(scr_window);

	/* Create main Vbox */
	vbox = gtk_vbox_new(FALSE, 0);
	gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scr_window),
	                                      vbox);
	gtk_widget_show(vbox);

	{ /* Create Control Panel Widget */
	  GtkWidget *ctl_panel = gtku_object_ctl_panel(obj);
	  RC_NULL_CHK(ctl_panel);
	  gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(ctl_panel),
	                     FALSE, FALSE, 0);
	  gtk_widget_show(ctl_panel);
	}

	/* Create bdf Save Button */
	if(obj->type & GTKU_OBJECT_FEM_NODE){
	  GtkWidget *button = gtk_button_new_with_label("Save bdf");
	  RC_NULL_CHK(button);
	  g_signal_connect(G_OBJECT(button), "clicked",
	                   G_CALLBACK(file_selection_bdf), (gpointer)obj);
	  gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(button),
	                     FALSE, FALSE, 10);
	  gtk_widget_show(button);
	}

	/* Create eps Save Button */
	if( (obj->type & GTKU_OBJECT_FEM_NODE) &&
	    (obj->type & GTKU_OBJECT_FEM_ELEM) ){ 
	  GtkWidget *button = gtk_button_new_with_label("Save eps");
	  RC_NULL_CHK(button);
	  g_signal_connect(G_OBJECT(button), "clicked",
	                   G_CALLBACK(file_selection_eps), (gpointer)obj);
	  gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(button),
	                     FALSE, FALSE, 10);
	  gtk_widget_show(button);
	}

	{ /* Create Tab Label and Notebook Tab */
	  GtkWidget *tab_label = gtk_label_new(obj->title);
	  RC_NULL_CHK(tab_label);
	  gtk_widget_show(tab_label);
	  gtk_notebook_append_page(GTK_NOTEBOOK(notebook),
	                           GTK_WIDGET(scr_window), tab_label);
	}

	return(NORMAL_RC);
}


static gboolean
file_selection_bdf(GtkWidget *widget, gpointer cb_data)
{
	GtkWidget *filew;

	filew =  gtk_file_selection_new("Save BDF(NAS) File");
	g_return_val_if_fail(filew, FALSE);
	g_object_set_data(G_OBJECT(filew), "gtku_object", cb_data);
	g_signal_connect(G_OBJECT(GTK_FILE_SELECTION(filew)->ok_button), "clicked",
	                 G_CALLBACK(filew_ok_bdf), (gpointer)filew);
	g_signal_connect_swapped(G_OBJECT(GTK_FILE_SELECTION(filew)->cancel_button),
	                         "clicked", G_CALLBACK(gtk_widget_destroy),
	                         G_OBJECT(filew));
	gtk_widget_show(filew);

	return(TRUE);
}

static gboolean
filew_ok_bdf(GtkWidget *widget, gpointer cb_data)
{
	if(save_bdf_file(widget, cb_data))
		gtk_widget_destroy(GTK_WIDGET(cb_data));
	return TRUE;
}

static gboolean
file_selection_eps(GtkWidget *widget, gpointer cb_data)
{
	GtkWidget *filew;

	filew =  gtk_file_selection_new("Save  EPS File");
	g_return_val_if_fail(filew, FALSE);
	g_object_set_data(G_OBJECT(filew), "gtku_object", cb_data);
	g_signal_connect(G_OBJECT(GTK_FILE_SELECTION(filew)->ok_button), "clicked",
	                 G_CALLBACK(filew_ok_eps), (gpointer)filew);
	g_signal_connect_swapped(G_OBJECT(GTK_FILE_SELECTION(filew)->cancel_button),
	                         "clicked", G_CALLBACK(gtk_widget_destroy),
	                         G_OBJECT(filew));
	gtk_widget_show(filew);

	return(TRUE);
}

static gboolean
filew_ok_eps(GtkWidget *widget, gpointer cb_data)
{
	if(save_eps_file(widget, cb_data))
		gtk_widget_destroy(GTK_WIDGET(cb_data));
	return TRUE;
}



static gboolean
save_bdf_file(GtkWidget *widget, gpointer cb_data)
{
	int ii1;
	FILE *fp;
	char *filename;
	GTKU_OBJECT *obj;
	double **mat;
	NODE_ARRAY new_node;

	/* get Object */
	obj = (GTKU_OBJECT *)g_object_get_data(G_OBJECT(cb_data), "gtku_object");
	g_return_val_if_fail(obj, FALSE);

	/* get File Name */
	filename = (char *)gtk_file_selection_get_filename(
	                   GTK_FILE_SELECTION(cb_data));
	g_print("file name = \"%s\"\n", filename);

	GTKU_BOOLEAN_RC_TRY( allocate_operation_matrix(&mat) );


	g_return_val_if_fail(obj->node, FALSE);
	RC_TRY_MAIN( allocate_node_array(obj->node->size, &new_node) );
	RC_TRY_MAIN( copy_node_array(*(obj)->node, &new_node) );

	mul_operation_matrix(obj->d_prop->rot_mat,   mat);
	mul_operation_matrix(obj->d_prop->shift_mat, mat);
	mul_operation_matrix(obj->d_prop->scale_mat, mat);

	for(ii1=0; ii1<obj->node->size; ii1++)
		new_node.array[ii1].p = operation3d(mat, obj->node->array[ii1].p);

	RC_TRY_MAIN( rc_fopen(filename, "w", &fp) );
	fprintf(fp, "BEGIN BULK\n");
	log_printf(5, "Write Node Data (Nastran Bulk Data Format)\n");
	RC_TRY( nst_output_node(fp, new_node) );

	if(obj->elem){
		log_printf(5, "Write Element Data (Nastran Bulk Data Format)\n");
		RC_TRY( nst_output_element(fp, *(obj)->elem) );
	}
	fprintf(fp, "ENDDATA\n");
	RC_TRY_MAIN( rc_fclose(fp) );

	return(TRUE);
}


static gboolean
save_eps_file(GtkWidget *widget, gpointer cb_data)
{
	FILE *fp;
	char *filename;
	GTKU_OBJECT *obj;
	double **mat;
	NODE_ARRAY *node;
	ELEMENT_ARRAY *elem, *surf;
	EPS_PROP prop;

	/* get Object */
	obj = (GTKU_OBJECT *)g_object_get_data(G_OBJECT(cb_data), "gtku_object");
	g_return_val_if_fail(obj, FALSE);

	/* get File Name */
	filename = (char *)gtk_file_selection_get_filename(
	                   GTK_FILE_SELECTION(cb_data));
	g_print("file name = \"%s\"\n", filename);


	GTKU_BOOLEAN_RC_TRY( allocate_operation_matrix(&mat) );

	/* 線は黒，白で塗りつぶし */
	prop.line_rgb_color[0] = 0.0;
	prop.line_rgb_color[1] = 0.0;
	prop.line_rgb_color[2] = 0.0;
	prop.fill_rgb_color[0] = 1.0;
	prop.fill_rgb_color[1] = 1.0;
	prop.fill_rgb_color[2] = 1.0;
	prop.line_width = 0.5;

	/* GTK_GLAREA_INFO existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	node = obj->node;
	elem = obj->elem;
	surf = obj->surf;

	mul_operation_matrix(gtku_gl->rot_mat, mat);
	mul_operation_matrix(gtku_gl->shift_mat, mat);
	if(gtku_gl->glarea->allocation.width <
	   gtku_gl->glarea->allocation.height){
		shift_operation_matrix(gtku_gl->glarea->allocation.width/2,
		                       gtku_gl->glarea->allocation.width/2,
		                       0, mat);
	}else{
		shift_operation_matrix(gtku_gl->glarea->allocation.height/2,
		                       gtku_gl->glarea->allocation.height/2,
		                       0, mat);
	}

	mul_operation_matrix(gtku_gl->scale_mat, mat);

	if(gtku_gl->object->obj[0]->dim == 2){
		GTKU_BOOLEAN_RC_TRY( rc_fopen(filename, "w", &fp) );
		GTKU_BOOLEAN_RC_TRY( output_eps(fp, prop, *elem, *node, mat) );
		GTKU_BOOLEAN_RC_TRY( rc_fclose(fp) );
	}else if(gtku_gl->object->obj[0]->dim == 3){
		GTKU_BOOLEAN_RC_TRY( rc_fopen(filename, "w", &fp) );
		GTKU_BOOLEAN_RC_TRY( reduced_element_array_order(surf) );
		GTKU_BOOLEAN_RC_TRY( output_eps(fp, prop, *surf, *node, mat) );
		GTKU_BOOLEAN_RC_TRY( rc_fclose(fp) );
	}else{
		g_warning("illegal model dimension\n");
	}

	return(TRUE);
}


/* Fix!! 応力，感度などを追加*/
/*
 * t : title char *
 * n : NODE_ARRAY
 * e : ELEMENT_ARRAY
 * s : ELEMENT_ARRAY (surface)
 * p : 2D PIXEL
 * v : 3D VOXEL
 */
static const char *accept_gtku_obj_type = "tnespv";
RC
gtku_add_object (GTKU_OBJECT *obj, const char *format, ...)
{
	int ii1;
	va_list ap;

	/* check format */
	g_return_val_if_fail(format_chk(format), ARG_ERROR_RC);

	/* create GL area infomation and initialize */
 	if(gtku_gl == NULL){
		gtku_gl = (GTKU_GL *)mm_alloc(sizeof(GTKU_GL));
		if(!gtku_gl) return(ALLOC_ERROR_RC);
		RC_TRY( init_gtku_gl() );
	}

	/* create mother of all object(nodem, element, pixel) and intialize */
 	if(gtku_gl->object == NULL){
		gtku_gl->object
			= (GTKU_OBJECT_ARRAY *)mm_alloc(sizeof(GTKU_OBJECT_ARRAY));
		if(!gtku_gl->object) return(ALLOC_ERROR_RC);
		init_gtku_object_array(gtku_gl->object);
	}

	/* Add new object */
	if(obj == NULL){
		/* realloc object array */
		gtku_gl->object->size++;
		gtku_gl->object->obj = (GTKU_OBJECT **)mm_realloc(gtku_gl->object->obj,
		                       (gtku_gl->object->size)*sizeof(GTKU_OBJECT *));
		if(!gtku_gl->object->obj) return(ALLOC_ERROR_RC);

		/* realloc object array */
		obj = gtku_gl->object->obj[gtku_gl->object->size - 1]
		    = (GTKU_OBJECT *)mm_alloc(sizeof(GTKU_OBJECT));
		if(!obj) return(ALLOC_ERROR_RC);
		init_gtku_object(obj);

		obj->d_prop = (GTKU_DRAW_PROP *)mm_alloc(sizeof(GTKU_DRAW_PROP));
		if(!obj->d_prop) return(ALLOC_ERROR_RC);
		init_gtku_object_prop(obj->d_prop, gtku_gl->object->size - 1);
	}

	va_start(ap, format);
	for(ii1=0; ii1<(int)strlen(format); ii1++){
		switch(format[ii1]){
		case 't':
			obj->title = obj->file_name = va_arg(ap, char *);
			break;
		case 'n':
			obj->type |= GTKU_OBJECT_FEM_NODE;
			obj->node = va_arg(ap, NODE_ARRAY *);
			break;
		case 'e':
			obj->type |= GTKU_OBJECT_FEM_ELEM;
			obj->elem = va_arg(ap, ELEMENT_ARRAY *);
			/* set object(element) dimension */
			RC_TRY( set_element_dim(obj->elem, &(obj->dim)) );
			break;
		case 'p':
			obj->type |= GTKU_OBJECT_IMAGE_PIXEL;
			obj->gr = va_arg(ap, GRAPHIC *);
			obj->dim = GTKU_OBJECT_2D;
			break;
		case 'v':
			obj->type |= GTKU_OBJECT_IMAGE_VOXEL;
			obj->gr = va_arg(ap, GRAPHIC *);
			obj->dim = GTKU_OBJECT_3D;
			break;
		default:
			break;
		}
	}
	va_end(ap);

	/* is object title and file name NULL  */
	if(obj->title == NULL){
		int length = sizeof(char)*32;
		obj->title = (char *)mm_alloc(length*sizeof(char));
		if(!obj->title) return(ALLOC_ERROR_RC);
		RC_NEG_ZERO_CHK( g_snprintf(obj->title, length-1, "object%d",
		                            gtku_gl->object->size - 1) );
		obj->title[length-1] = (char)'\0';  /* 念の為 */
		obj->file_name = (char *)mm_alloc(length*sizeof(char));
		if(!obj->file_name) return(ALLOC_ERROR_RC);
		RC_NEG_ZERO_CHK( g_snprintf(obj->file_name, length-1, "save_file%d",
		                            gtku_gl->object->size - 1) );
		obj->file_name[length-1] = (char)'\0';  /* 念の為 */
	}

	/* surface need ? */
	if( ((obj->dim == GTKU_OBJECT_2D) || (obj->dim == GTKU_OBJECT_3D)) &&
	     (obj->type & GTKU_OBJECT_FEM_ELEM) &&
	    !(obj->type & GTKU_OBJECT_FEM_SURF) ){
		obj->surf = (ELEMENT_ARRAY *)mm_alloc(sizeof(ELEMENT_ARRAY));
		if(!obj->surf) return(ALLOC_ERROR_RC);
		RC_TRY(extract_surface(*(obj->elem), obj->surf));
		obj->type |= GTKU_OBJECT_FEM_SURF;
	}
#ifdef GTKU_VERBOSE_MESSAGE
	fprintf(stderr, " title = %s\n", obj->title);
	print_gtku_object_type(stderr, obj->type);
#endif /* GTKU_VERBOSE_MESSAGE */

	return(NORMAL_RC);
}


static gboolean
format_chk (const char *format)
{
	int ii1, ii2;
	gboolean accept_sw;
	unsigned int mask = 0;

	if(format == NULL) return(FALSE);
	if(strlen(format) == 0) return(FALSE);

	for(ii1=0; format[ii1]; ii1++){
		accept_sw =  FALSE;
		for(ii2=0; accept_gtku_obj_type[ii2]; ii2++){
			if(format[ii1] == accept_gtku_obj_type[ii2]){
				if((mask ^ (1<<ii2)) == 0){ /* same format type  ERROR!! */
					g_warning("Invalid format strings. Same tyep exist.\n");
					return(FALSE);
				}else{
					accept_sw = TRUE;
					mask += 1 << ii2;
				}
			}
		}
		if(accept_sw != TRUE){
			g_warning("Invalid format strings. Unknown type.\n");
			return(FALSE);
		}
	}

	return(TRUE);
}


static RC
set_element_dim (ELEMENT_ARRAY *elem, GTKU_OBJECT_DIM *dim)
{
	int dimension;

	RC_NULL_CHK(elem);
	RC_NEG_CHK( dimension = analysis_dim(*elem) ); 

	     if(dimension == 1) *dim = GTKU_OBJECT_1D;
	else if(dimension == 2) *dim = GTKU_OBJECT_2D;
	else if(dimension == 3) *dim = GTKU_OBJECT_3D;
	else return(UNKNOWN_ERROR_RC); 

	return(NORMAL_RC);
}


GtkWidget *
gtku_glarea_new (GtkWidget *window)
{
	GdkGLConfig *glconfig;
	
	/* get information about a OpenGL frame buffer configuration */
	glconfig = gtku_get_gdkglconfig();
	g_return_val_if_fail(glconfig, NULL);

#ifdef GTKU_VERBOSE_MESSAGE
	/* print OpenGL attibute */
	print_gl_config_attrib(stderr, glconfig);

	{ /* Query OpenGL extension version. */
	 gint major, minor;
	 gdk_gl_query_version(&major, &minor);
	 fprintf(stderr, "\nOpenGL extension version : <%d.%d>\n", major, minor);
	}
#endif /* GTKU_VERBOSE_MESSAGE */

	/* create drawing area and width and height setting */
	gtku_gl->glarea = GTK_WIDGET(gtk_drawing_area_new());
	g_return_val_if_fail(gtku_gl->glarea, NULL);
	gtk_widget_set_size_request(GTK_WIDGET(gtku_gl->glarea),
	                            gtku_gl->glarea_w_size,
	                            gtku_gl->glarea_h_size);

	/* set OpenGL-capability to the widget. */
	gtk_widget_set_gl_capability(GTK_WIDGET(gtku_gl->glarea), glconfig,
	                             NULL, TRUE, GDK_GL_RGBA_TYPE);
	g_object_unref( G_OBJECT(glconfig) );

	gtk_widget_set_events(GTK_WIDGET(gtku_gl->glarea),
	                      GDK_EXPOSURE_MASK |
	                      GDK_BUTTON_PRESS_MASK |
	                      GDK_BUTTON_RELEASE_MASK |
	                      GDK_POINTER_MOTION_MASK |
	                      GDK_POINTER_MOTION_HINT_MASK);

	/* シグナルはここに書かれた順番に試される．特にbutton_press */
	/* motion_notifyの順番には注意                              */
	g_signal_connect(G_OBJECT(gtku_gl->glarea), "expose_event",
	                 G_CALLBACK(gtku_gl_draw), NULL);
	g_signal_connect(G_OBJECT(gtku_gl->glarea), "configure_event",
	                 G_CALLBACK(gtku_gl_reshape), NULL);
	g_signal_connect(G_OBJECT(gtku_gl->glarea), "realize",
	                 G_CALLBACK(gtku_gl_realize), NULL);

	g_signal_connect(G_OBJECT(gtku_gl->glarea), "button_press_event",
	                 G_CALLBACK(gtku_gl_button_press), NULL);
	g_signal_connect(G_OBJECT(gtku_gl->glarea), "button_release_event",
	                 G_CALLBACK(gtku_gl_button_release), NULL);
	g_signal_connect(G_OBJECT(gtku_gl->glarea), "scroll_event",
	                 G_CALLBACK(gtku_gl_scroll), NULL);
	g_signal_connect(G_OBJECT(gtku_gl->glarea), "motion_notify_event",
	                 G_CALLBACK(gtku_gl_motion_notify), NULL);

	g_signal_connect(G_OBJECT(gtku_gl->glarea), "map_event",
	                 G_CALLBACK(gtku_gl_map), NULL);
	g_signal_connect(G_OBJECT(gtku_gl->glarea), "unmap_event",
	                 G_CALLBACK(gtku_gl_unmap), NULL);

	g_signal_connect(G_OBJECT(gtku_gl->glarea), "visibility_notify_event",
	                 G_CALLBACK(gtku_gl_visibility_notify), NULL);
	g_signal_connect(G_OBJECT(gtku_gl->glarea), "destroy",
	                 G_CALLBACK(gtku_gl_destroy), NULL);

	/* ショートカットメニューを優先するために"_after"は必須 */
	g_signal_connect_after(G_OBJECT(window), "key_press_event",
	                       G_CALLBACK(gtku_gl_key_press), NULL);
	g_signal_connect_after(G_OBJECT(window), "key_release_event",
	                       G_CALLBACK(gtku_gl_key_release), NULL);

	return(gtku_gl->glarea);
}


GdkGLConfig *
gtku_get_gdkglconfig (void)
{
	GdkGLConfig *glconfig;
	/* attribute listのTRUE/FALSEは関係ない．listに有るものがONとなる    */
	/* 例えば，GDK_GL_DOUBLEBUFFERを無効にしたければリストから外すこと． */ 
	const int attr_list[] = {GDK_GL_RGBA,         TRUE,
	                         GDK_GL_DOUBLEBUFFER, TRUE,
	                         GDK_GL_DEPTH_SIZE,   16,
	                         GDK_GL_ATTRIB_LIST_NONE};

	/* Check if OpenGL is supported */
	if( !gdk_gl_query_extension() )
		g_error("OpenGL not supported\n");

	/* Check double buffered mode */
	glconfig = gdk_gl_config_new(attr_list);
	if(glconfig == NULL){
		g_warning("*** Cannot find the double-buffered visual.\n");
		g_warning("*** Trying single-buffered visual.\n");

		/* Try single-buffered visual */
		glconfig = gdk_gl_config_new_by_mode(GDK_GL_MODE_RGB |
		                                     GDK_GL_MODE_DEPTH |
		                                     GDK_GL_MODE_SINGLE);
		if(glconfig == NULL)
			g_error("*** No appropriate OpenGL-capable visual found.\n");
	}

	return(glconfig);
}


/* ツールバーのボタンを作りツールバーに登録する関数 */
gboolean
gtku_toolbar_rotation (GtkWidget *widget, GtkWidget *toolbar)
{
	int ii1;
	GtkWidget *(*func[])(GtkWidget *, GtkWidget *, GtkSignalFunc, gpointer)
	   = { gtku_toolbar_rotation_xy,    /* FRONT view */
	       gtku_toolbar_rotation_yz,    /* RIGHT view */
	       gtku_toolbar_rotation_zx,    /* TOP   view */
	       gtku_toolbar_rotation_xyz }; /* XYZ   view */
     
	for(ii1=0; ii1<(sizeof(func)/sizeof(func[0])); ii1++){
		g_return_val_if_fail(func[ii1](widget, toolbar,
		                     GTK_SIGNAL_FUNC(gtku_set_model_view_direction),
		                     (gpointer)ii1), FALSE);
	}

	return(TRUE);
}

#define GTKU_TOOLBAR_ROTATION_COD(func_name, tip, xpm)\
GtkWidget *func_name(GtkWidget *widget, GtkWidget *toolbar,\
                     GtkSignalFunc func, gpointer data)\
{\
    return( gtk_toolbar_append_element( GTK_TOOLBAR(toolbar),\
            GTK_TOOLBAR_CHILD_BUTTON, NULL, NULL, tip, NULL,\
            gtku_create_gtkimage_from_xpm(widget, (gchar **)xpm),\
            (GtkSignalFunc)func, (gpointer)data) );\
}

GTKU_TOOLBAR_ROTATION_COD(gtku_toolbar_rotation_xy,  "Front View", xy_xpm)
GTKU_TOOLBAR_ROTATION_COD(gtku_toolbar_rotation_yz,  "Side View",  yz_xpm)
GTKU_TOOLBAR_ROTATION_COD(gtku_toolbar_rotation_zx,  "Top View",   zx_xpm)
GTKU_TOOLBAR_ROTATION_COD(gtku_toolbar_rotation_xyz, "XYZ View",  xyz_xpm)


gboolean
gtku_set_model_view_direction (GtkWidget *widget, gpointer cb_data)
{
	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	init_operation_matrix(gtku_gl->rot_mat);
	switch((gint)cb_data){
	case 0:
		rotation_operation_matrix_y(0.0, gtku_gl->rot_mat);
		break;
	case 1:
		rotation_operation_matrix_y(-PI/2, gtku_gl->rot_mat);
		break;
	case 2:
		rotation_operation_matrix_x(PI/2, gtku_gl->rot_mat);
		break;
	case 3:
		rotation_operation_matrix_y(-PI/4, gtku_gl->rot_mat);
		rotation_operation_matrix_x( PI/4, gtku_gl->rot_mat);
		break;
	default:
		break;
	}

	gtk_widget_queue_draw(gtku_gl->glarea);

	return(TRUE);
}


gboolean
gtku_toolbar_drawtype (GtkWidget *widget, GtkWidget *toolbar)
{
	int ii1;
	GtkWidget *group = NULL;
	GtkWidget *(*func[])(GtkWidget *, GtkWidget *, GtkWidget *,
	                     GtkSignalFunc, gpointer)
	   = { gtku_toolbar_drawtype_polygon,    gtku_toolbar_drawtype_solid,
		   gtku_toolbar_drawtype_wireframe,  gtku_toolbar_drawtype_gridframe,
		   gtku_toolbar_drawtype_surfpoint,  gtku_toolbar_drawtype_gridpoint,
		   gtku_toolbar_drawtype_hiddenline};

	GTKU_DRAW_TYPE type[]
	    = { GTKU_DRAW_POLYGON,    GTKU_DRAW_SOLID,      GTKU_DRAW_WIRE_FRAME,
	        GTKU_DRAW_GRID_FRAME, GTKU_DRAW_SURF_POINT, GTKU_DRAW_GRID_POINT,
	        GTKU_DRAW_HIDDENLINE};

	/* draw_type_radio[7] is global variable */
	for(ii1=0; ii1<(sizeof(func)/sizeof(func[0])); ii1++){
		if(ii1 != 0) group = draw_type_radio[ii1-1];
		draw_type_radio[ii1] = func[ii1](widget, toolbar, group,
		                       GTK_SIGNAL_FUNC(gtku_set_model_draw_type),
		                       (gpointer)(int)type[ii1]);
		g_return_val_if_fail(draw_type_radio[ii1], FALSE);
	}
	/* default on polygon */
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(draw_type_radio[0]), TRUE);

	return(TRUE);
}

#define GTKU_TOOLBAR_DRAWTYPE_TYPE(func_name, tip, xpm)\
GtkWidget *func_name(GtkWidget *widget, GtkWidget *toolbar, GtkWidget *group,\
                     GtkSignalFunc func, gpointer data)\
{\
    return( gtk_toolbar_append_element(GTK_TOOLBAR(toolbar),\
                GTK_TOOLBAR_CHILD_RADIOBUTTON, group, NULL, tip, NULL,\
                gtku_create_gtkimage_from_xpm(widget, (gchar **)xpm),\
                (GtkSignalFunc)func, (gpointer)data) );\
}

GTKU_TOOLBAR_DRAWTYPE_TYPE(
	gtku_toolbar_drawtype_polygon,    "Polygon",       polygon_xpm)
GTKU_TOOLBAR_DRAWTYPE_TYPE(
	gtku_toolbar_drawtype_solid,      "Solid",         solid_xpm)
GTKU_TOOLBAR_DRAWTYPE_TYPE(
	gtku_toolbar_drawtype_wireframe,  "Wire Frame",    wire_frame_xpm)
GTKU_TOOLBAR_DRAWTYPE_TYPE(
	gtku_toolbar_drawtype_gridframe,  "Grid Frame",    grid_frame_xpm)
GTKU_TOOLBAR_DRAWTYPE_TYPE(
	gtku_toolbar_drawtype_surfpoint,  "Surface Point", surf_point_xpm)
GTKU_TOOLBAR_DRAWTYPE_TYPE(
	gtku_toolbar_drawtype_gridpoint,  "Grid Point",    grid_point_xpm)
GTKU_TOOLBAR_DRAWTYPE_TYPE(
	gtku_toolbar_drawtype_hiddenline, "Hidden Line",   hiddenline_xpm)



/* この方法はあまりスマートとは言えない．メニューとボタン以外の第三の */
/* セレクションを追加したい場合，以下を書き換える必要がある．         */
/* タイマーを使った方が楽かもしれない                                 */
gboolean
gtku_set_model_draw_type (GtkWidget *widget, gpointer cb_data)
{
	static gboolean draw_button_sw = FALSE;
	static gboolean draw_menu_sw = FALSE;
	char *path[] ={"/View/DrawType/Polygon",
	               "/View/DrawType/Solid",
	               "/View/DrawType/Wire Frame",
	               "/View/DrawType/Grid Frame",
	               "/View/DrawType/Surface Point",
	               "/View/DrawType/Grid Point",
	               "/View/DrawType/HiddenLine"};

	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	/* radio switch はそのボタンなどが押された時，2回 toggle signalが */
	/* 送られる． 1回目のcallback dataは前回データが送られる．        */
	if( GTK_IS_RADIO_BUTTON(widget) ){ /* toolbarのボタンが押された時 */
		if(draw_menu_sw == TRUE){
			int num;
			GtkWidget *item;
			GtkCheckMenuItem *menu;
			gtku_gl->draw_type = (gint)cb_data;
			num = drawtype2number(gtku_gl->draw_type);
			item = gtk_item_factory_get_item(item_factory,path[num]);
			menu = GTK_CHECK_MENU_ITEM(item);
			if(!gtk_check_menu_item_get_active(menu)){
				gtk_check_menu_item_set_active(menu, TRUE);
				draw_menu_sw = !draw_menu_sw;
			}
		}else{
			draw_menu_sw = !draw_menu_sw;
		}
	}else if( GTK_IS_MENU_ITEM(widget) ){ /* メニューが選択された時 */
	    if(draw_button_sw == TRUE){
			int num;
			GtkToggleButton *button;
			gtku_gl->draw_type = (gint)cb_data;
			num = drawtype2number(gtku_gl->draw_type);
			button = GTK_TOGGLE_BUTTON(draw_type_radio[num]);
			if(!gtk_toggle_button_get_active(button)){
				gtk_toggle_button_set_active(button, TRUE);
				draw_button_sw = !draw_button_sw;
			}
		}else{
			draw_button_sw = !draw_button_sw;
		}
	}else{
		GTKU_WARNING("set_model_view_type()", "Unknown Widget");
	}

	if(draw_button_sw && draw_menu_sw){
		draw_button_sw = !draw_button_sw;
		draw_menu_sw   = !draw_menu_sw;
	}

	gtk_widget_queue_draw(gtku_gl->glarea);

	return(TRUE);
}


static int
drawtype2number (GTKU_DRAW_TYPE type)
{
	switch((int)type){
	case GTKU_DRAW_POLYGON:      return(0);   
	case GTKU_DRAW_SOLID:        return(1);
	case GTKU_DRAW_WIRE_FRAME:   return(2);
	case GTKU_DRAW_GRID_FRAME:   return(3);
	case GTKU_DRAW_SURF_POINT:   return(4);
	case GTKU_DRAW_GRID_POINT:   return(5);
	case GTKU_DRAW_HIDDENLINE:   return(6);
	default:
		GTKU_WARNING("draw_type2number()", "Illegal drawing type");
		return(0); /* for safe, don't return except 1-6 */
	}
}

gboolean
gtku_toolbar_clear (GtkWidget *widget, GtkWidget *toolbar)
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
	int ii1;
	char *tip[] = {"Initialize", "Fit", "Zoom IN", "Zoom OUT"};
	char *stock_id[] = { GTK_STOCK_ZOOM_100, GTK_STOCK_ZOOM_FIT,
	                     GTK_STOCK_ZOOM_IN,  GTK_STOCK_ZOOM_OUT };

	for(ii1=0; ii1<(sizeof(stock_id)/sizeof(stock_id[0])); ii1++){
		GtkWidget *elem;
		elem = gtk_toolbar_append_element(
		           GTK_TOOLBAR(toolbar), GTK_TOOLBAR_CHILD_BUTTON, NULL,
		           NULL, tip[ii1], NULL, gtk_image_new_from_stock(stock_id[ii1],
		           GTK_ICON_SIZE_LARGE_TOOLBAR),
		           GTK_SIGNAL_FUNC(gtku_set_model_scale), (gpointer)ii1 );
		g_return_val_if_fail(elem, FALSE);
	}

	return(TRUE);
}


gboolean
gtku_set_model_scale (GtkWidget *widget, gpointer cb_data)
{
	double zoom_p = 1.1;
	double zoom_m = 0.9;
	double fit_size = 0.0;

	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	/* ref  fit zomm re-zoom */
	switch((int)cb_data){
	case 0: /* 移動とスケールを元に戻す */
		init_operation_matrix(gtku_gl->scale_mat);
		init_operation_matrix(gtku_gl->shift_mat);
		break;
	case 1: /* 画面サイズに合わせる */
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
		break;
	case 2: /* Zoom IN  */
		scale_operation_matrix(zoom_p, zoom_p, zoom_p, gtku_gl->scale_mat);
		break;
	case 3: /* Zoom OUT */
		scale_operation_matrix(zoom_m, zoom_m, zoom_m, gtku_gl->scale_mat);
		break;
	default:
		break;
	}

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


gboolean
gtku_set_world_plane (GtkWidget *widget, gpointer cb_data)
{
	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	/* ref  fit zomm re-zoom */
	switch((int)cb_data){
	case 0: /* draw world plane  */
		gtku_gl->draw_plane_sw = !(gtku_gl->draw_plane_sw);
		break;
	case 1: /* yz world plane */
		gtku_gl->draw_plane_yz_sw = !(gtku_gl->draw_plane_yz_sw);
		break;
	case 2: /* zx world plane */
		gtku_gl->draw_plane_zx_sw = !(gtku_gl->draw_plane_zx_sw);
		break;
	case 3: /* xy world plane */
		gtku_gl->draw_plane_xy_sw = !(gtku_gl->draw_plane_xy_sw);
		break;
	default:
		break;
	}

	if(gtku_gl->glarea)
		gtk_widget_queue_draw(gtku_gl->glarea);

	return(TRUE);
}


/* ItemFactoryのコールバック関数(GtkItemFactoryCallback1)の引数が */
/* このファイル内のコールバック関数と違うのでそれを変換する       */
#define CONV_CB_FUNC_FOR_IF(func) \
static void func##_for_IF(gpointer cb_data, guint cb_action, GtkWidget *widget)\
{func(widget, (gpointer)cb_action);}

CONV_CB_FUNC_FOR_IF(gtku_set_model_scale)
CONV_CB_FUNC_FOR_IF(gtku_set_model_draw_type)
CONV_CB_FUNC_FOR_IF(gtku_set_model_view_direction)
CONV_CB_FUNC_FOR_IF(gtku_set_world_plane)


static GtkItemFactoryEntry menu_items[] = {
 /* File ->  */
 {("/_File"), NULL, NULL, 0, "<Branch>", NULL},
 { "/File/tear1", NULL, NULL, 0, "<Tearoff>", NULL },
 {("/File/_Quit"),"<control>Q",gtk_main_quit, 0, "<StockItem>", GTK_STOCK_QUIT},
 /* View ->  */
 {("/View"),          NULL, NULL, 0, "<Branch>", NULL},
 { "/View/tear2", NULL, NULL, 0, "<Tearoff>", NULL },
 /* View -> Scale */
 {("/View/Scale"), NULL, NULL, 0, "<Branch>", NULL},
 {("/View/Scale/Initialize"), NULL,
  gtku_set_model_scale_for_IF, 0, "<Item>", NULL},
 {("/View/Scale/Fit"), NULL,
  gtku_set_model_scale_for_IF, 1, "<Item>", NULL},
 {("/View/Scale/Zoom-IN"), NULL,
  gtku_set_model_scale_for_IF, 2, "<Item>", NULL},
 {("/View/Scale/Zoom-OUT"), NULL,
  gtku_set_model_scale_for_IF, 3, "<Item>", NULL},
 /* View -> Direction */
 {("/View/Direction/ XY-Plane"), NULL,
  gtku_set_model_view_direction_for_IF, 0, NULL, NULL},
 {("/View/Direction/ YZ-Plane"), NULL,
  gtku_set_model_view_direction_for_IF, 1, NULL, NULL},
 {("/View/Direction/ ZX-Plane"), NULL,
  gtku_set_model_view_direction_for_IF, 2, NULL, NULL},
 {("/View/Direction/XYZ-Plane"), NULL,
  gtku_set_model_view_direction_for_IF, 3, NULL, NULL},
 /* View -> DrawType */
 {("/View/DrawType"), NULL, NULL, 0, "<Branch>", NULL},
 {("/View/DrawType/Polygon"), NULL, gtku_set_model_draw_type_for_IF,
  GTKU_DRAW_POLYGON, "<RadioItem>", NULL},
 {("/View/DrawType/Solid"), NULL, gtku_set_model_draw_type_for_IF,
  GTKU_DRAW_SOLID, "/View/DrawType/Polygon", NULL},
 {("/View/DrawType/Wire Frame"), NULL, gtku_set_model_draw_type_for_IF,
  GTKU_DRAW_WIRE_FRAME, "/View/DrawType/Solid", NULL},
 {("/View/DrawType/Grid Frame"), NULL, gtku_set_model_draw_type_for_IF,
  GTKU_DRAW_GRID_FRAME, "/View/DrawType/Wire Frame", NULL},
 {("/View/DrawType/Surface Point"), NULL, gtku_set_model_draw_type_for_IF,
  GTKU_DRAW_SURF_POINT, "/View/DrawType/Grid Frame", NULL},
 {("/View/DrawType/Grid Point"), NULL, gtku_set_model_draw_type_for_IF,
  GTKU_DRAW_GRID_POINT, "/View/DrawType/Surface Point", NULL},
 {("/View/DrawType/HiddenLine"), NULL, gtku_set_model_draw_type_for_IF,
  GTKU_DRAW_HIDDENLINE, "/View/DrawType/Grid Point", NULL},
 /* View -> World Plane */
 {("/View/World Plane/Draw Plane"), NULL,
  gtku_set_world_plane_for_IF, 0, "<CheckItem>", NULL},
 {"/View/World Plane/separator", NULL, NULL, 0, "<Separator>", NULL},
 {("/View/World Plane/YZ-Plane"), NULL,
  gtku_set_world_plane_for_IF, 1, "<CheckItem>", NULL},
 {("/View/World Plane/ZX-Plane"), NULL,
  gtku_set_world_plane_for_IF, 2, "<CheckItem>", NULL},
 {("/View/World Plane/XY-Plane"), NULL,
  gtku_set_world_plane_for_IF, 3, "<CheckItem>", NULL},
 /* LastBranch */
 {"/_Help",       NULL, NULL, 0, "<LastBranch>"},
 {"/_Help/About", NULL, NULL, 0, NULL}
};


GtkWidget *
gtku_create_model_view_menu (GtkWidget *widget)
{
	GtkAccelGroup *accel_group;
	gint nmenu_items = sizeof(menu_items)/sizeof(menu_items[0]);

	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	accel_group = gtk_accel_group_new();
	item_factory = gtk_item_factory_new(GTK_TYPE_MENU_BAR,"<main>",accel_group);
	gtk_item_factory_create_items(item_factory, nmenu_items, menu_items, NULL);

	/* Attach the new accelerator group to the window. */
	gtk_window_add_accel_group(GTK_WINDOW(widget), accel_group);
 
	/* World Plane default setting */
	if(gtku_gl->draw_plane_sw){
		GtkWidget *item = gtk_item_factory_get_item(item_factory,
	                      "/View/World Plane/Draw Plane");
		gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(item), TRUE);
		gtku_gl->draw_plane_sw = !gtku_gl->draw_plane_sw;
	}
	if(gtku_gl->draw_plane_yz_sw){
		GtkWidget *item = gtk_item_factory_get_item(item_factory,
	                      "/View/World Plane/YZ-Plane");
		gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(item), TRUE);
		gtku_gl->draw_plane_yz_sw = !gtku_gl->draw_plane_yz_sw; 
	}
	if(gtku_gl->draw_plane_zx_sw){
		GtkWidget *item = gtk_item_factory_get_item(item_factory,
	                      "/View/World Plane/ZX-Plane");
		gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(item), TRUE);
		gtku_gl->draw_plane_zx_sw = !gtku_gl->draw_plane_zx_sw; 
	}
	if(gtku_gl->draw_plane_xy_sw){
		GtkWidget *item = gtk_item_factory_get_item(item_factory,
	                      "/View/World Plane/XY-Plane");
		gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(item), TRUE);
		gtku_gl->draw_plane_xy_sw = !gtku_gl->draw_plane_xy_sw; 
	}

	return(gtk_item_factory_get_widget(item_factory, "<main>"));
}


GtkWidget *
gtku_get_model_view_menu (void)
{
	GtkItemFactory *item_factory;
	gint nmenu_items = sizeof(menu_items)/sizeof(menu_items[0]);
	GtkWidget *menu;
     
	/* Same again, not bothering with the accelerators */
	item_factory = gtk_item_factory_new(GTK_TYPE_OPTION_MENU, "<main>", NULL);
	gtk_item_factory_create_items(item_factory, nmenu_items, menu_items, NULL);
	menu = gtk_item_factory_get_widget(item_factory, "<main>");

	return(menu);
}


GtkWidget *
gtku_object_ctl_panel (GTKU_OBJECT *obj)
{
	GtkWidget *vbox;

	/* Create main Vbox */
	vbox = gtk_vbox_new(FALSE, 2);
	gtk_widget_show(vbox);

	{ /* Create operation matrix lock frame */
	 GtkWidget *frame = gtku_object_ctl_panel_matrix(obj);
	 g_return_val_if_fail(frame, NULL);
	 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(frame), FALSE, FALSE, 0);
	 gtk_widget_show(frame);
	}

	{ /* Create local Drawing type */
	 GtkWidget *frame = gtku_object_ctl_panel_draw_type(obj);
	 g_return_val_if_fail(frame, NULL);
	 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(frame), FALSE, FALSE, 0);
	 gtk_widget_show(frame);
	}

	{ /* Create local Quantity of Display */
	 GtkWidget *frame = gtku_object_ctl_panel_quantity(obj);
	 g_return_val_if_fail(frame, NULL);
	 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(frame), FALSE, FALSE, 0);
	 gtk_widget_show(frame);
	}

	if(obj->type & (GTKU_OBJECT_IMAGE_PIXEL | GTKU_OBJECT_IMAGE_VOXEL)){
		{ /* Create local Image Dirction */
		 GtkWidget *frame = gtku_object_ctl_panel_image_dir(obj);
		 g_return_val_if_fail(frame, NULL);
		 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(frame), FALSE, FALSE, 0);
		 gtk_widget_show(frame);
		}

		{ /* Create local Use Contour Map */
		 GtkWidget *frame = gtku_object_ctl_panel_image_contour(obj);
		 g_return_val_if_fail(frame, NULL);
		 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(frame), FALSE, FALSE, 0);
		 gtk_widget_show(frame);
		}

		{ /* Create local Image Range */
		 GtkWidget *frame = gtku_object_ctl_panel_image_range(obj);
		 g_return_val_if_fail(frame, NULL);
		 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(frame), FALSE, FALSE, 0);
		 gtk_widget_show(frame);
		}

		{ /* Create local Image Trans */
		 GtkWidget *frame = gtku_object_ctl_panel_image_trans(obj);
		 g_return_val_if_fail(frame, NULL);
		 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(frame), FALSE, FALSE, 0);
		 gtk_widget_show(frame);
		}
	}else{
		{ /* Create local Drawing Color & Trans */
		 GtkWidget *frame = gtku_object_ctl_panel_color(obj);
		 g_return_val_if_fail(frame, NULL);
		 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(frame), FALSE, FALSE, 0);
		 gtk_widget_show(frame);
		}

		{ /* Create local Drawing Number */
		 GtkWidget *frame = gtku_object_ctl_panel_draw_num(obj);
		 g_return_val_if_fail(frame, NULL);
		 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(frame), FALSE, FALSE, 0);
		 gtk_widget_show(frame);
		}

		{ /* Create local Drawing Size */
		 GtkWidget *frame = gtku_object_ctl_panel_size(obj);
		 g_return_val_if_fail(frame, NULL);
		 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(frame), FALSE, FALSE, 0);
		 gtk_widget_show(frame);
		}
	}

	return(vbox);
}


GtkWidget *
gtku_object_ctl_panel_matrix (GTKU_OBJECT *obj)
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
	                  G_CALLBACK(gtku_set_local_rot_matrix_lock),
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
	                  G_CALLBACK(gtku_set_local_scale_matrix_lock),
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
	                  G_CALLBACK(gtku_set_local_shift_matrix_lock),
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
	                  G_CALLBACK(gtku_set_boolean_switch),
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

GTKU_SET_LOCAL_MATRIX_LOCK(gtku_set_local_rot_matrix_lock,
                           gtku_gl->rot_lock_count)
GTKU_SET_LOCAL_MATRIX_LOCK(gtku_set_local_scale_matrix_lock,
                           gtku_gl->scale_lock_count)
GTKU_SET_LOCAL_MATRIX_LOCK(gtku_set_local_shift_matrix_lock,
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
	GTKU_BOOLEAN_RC_TRY( print_gtku_object_prop_operation_mat(fp,
	                                                          obj->d_prop) );
	GTKU_BOOLEAN_RC_TRY( rc_fclose(fp) );
	log_printf(5, "Done\n");

	return(TRUE);
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
	GTKU_BOOLEAN_RC_TRY( read_gtku_object_prop_operation_mat(fp,
	                                                         obj->d_prop) );
	GTKU_BOOLEAN_RC_TRY( rc_fclose(fp) );
	log_printf(5, "Done\n");

	return(TRUE);
}



GtkWidget *
gtku_object_ctl_panel_draw_type (GTKU_OBJECT *obj)
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
		GTKU_WARNING("gtku_object_ctl_panel_draw_type()",
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
		g_object_set_data( G_OBJECT(check), "gtku_object_draw_type",
		                   &(obj->d_prop->draw_type) );
		g_signal_connect( G_OBJECT(check), "toggled",
		                  G_CALLBACK(gtku_set_local_object_draw_type),
		                  (gpointer)(type_name)[ii1] );
		gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(check), FALSE, FALSE, 1);
		gtk_widget_show(check);
	}

	return(frame);
}


gboolean
gtku_set_local_object_draw_type (GtkWidget *widget, gpointer cb_data)
{
	GTKU_DRAW_TYPE  *draw_type;

	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	draw_type = (GTKU_DRAW_TYPE *)g_object_get_data(G_OBJECT(widget),
	                                                "gtku_object_draw_type");
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
gtku_object_ctl_panel_quantity (GTKU_OBJECT *obj)
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
	                 G_CALLBACK(gtku_set_local_quantity),
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
gtku_set_local_quantity (GtkWidget *widget, gpointer cb_data)
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
gtku_object_ctl_panel_image_dir (GTKU_OBJECT *obj)
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
		option = gtku_local_image_direction_option(obj);
		g_return_val_if_fail(option, NULL);
		gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(option), FALSE, FALSE, 0);
		gtk_widget_show(option);
	}

	return(frame);
}


GtkWidget *
gtku_local_image_direction_option (GTKU_OBJECT *obj)
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
		                  G_CALLBACK(gtku_set_local_image_direction),
		                  (gpointer)(image_dir)[ii1] );
		gtk_widget_show(menuitem);
	}

	gtk_option_menu_set_menu(GTK_OPTION_MENU(option_menu), menu);

	return(option_menu);
}


gboolean
gtku_set_local_image_direction (GtkWidget *widget, gpointer cb_data)
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
		GTKU_WARNING("gtku_set_image_direction()",
		             "Unknown Image Object Direction!!");
	}
	obj->d_prop->init_image_dir = FALSE;
	gtk_widget_queue_draw(gtku_gl->glarea);

	return(TRUE);
}


GtkWidget *
gtku_object_ctl_panel_image_contour (GTKU_OBJECT *obj)
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
		option = gtku_local_image_contour_option(obj);
		g_return_val_if_fail(option, NULL);
		gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(option), FALSE, FALSE, 0);
		gtk_widget_show(option);
	}

	return(frame);
}


GtkWidget *
gtku_local_image_contour_option (GTKU_OBJECT *obj)
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
		                  G_CALLBACK(gtku_set_local_image_contour),
		                  (gpointer)(contour_map)[ii1] );
		gtk_widget_show(menuitem);
	}

	gtk_option_menu_set_menu(GTK_OPTION_MENU(option_menu), menu);

	return(option_menu);
}


gboolean
gtku_set_local_image_contour (GtkWidget *widget, gpointer cb_data)
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
		GTKU_WARNING("gtku_set_image_direction()",
		             "Unknown Contour Type!!");
	}
	gtk_widget_queue_draw(gtku_gl->glarea);

	return(TRUE);
}


GtkWidget *
gtku_object_ctl_panel_image_range (GTKU_OBJECT *obj)
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
	 GtkWidget *scale = gtku_local_image_range_with_label(
	                       obj, &(obj->d_prop->img_threshold_max), "MAX");
	 g_return_val_if_fail(scale, NULL);
	 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(scale), FALSE, FALSE, 5);
	 gtk_widget_show(scale);
	}

	{ /* Create min range  */
	 GtkWidget *scale = gtku_local_image_range_with_label(
	                       obj, &(obj->d_prop->img_threshold_min), "MIN");
	 g_return_val_if_fail(scale, NULL);
	 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(scale), FALSE, FALSE, 5);
	 gtk_widget_show(scale);
	}

	return(frame);
}


static GtkWidget *
gtku_local_image_range_with_label (GTKU_OBJECT *obj, double *range, char *label)
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
	                 G_CALLBACK(gtku_set_local_image_range),(gpointer)(range));

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
gtku_set_local_image_range (GtkWidget *widget, gpointer cb_data)
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
gtku_object_ctl_panel_image_trans (GTKU_OBJECT *obj)
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
	                 G_CALLBACK(gtku_set_local_image_trans),
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
gtku_set_local_image_trans (GtkWidget *widget, gpointer cb_data)
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
gtku_object_ctl_panel_color (GTKU_OBJECT *obj)
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
	 GtkWidget *trans_frame = gtku_object_ctl_panel_trans_check(obj);
	 g_return_val_if_fail(trans_frame, NULL);
	 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(trans_frame), FALSE,FALSE, 1);
	 gtk_widget_show(trans_frame);
	}

	{ /* Create color selection button  */
	 GtkWidget *colorsel_frame = gtku_object_ctl_panel_color_button(obj);
	 g_return_val_if_fail(colorsel_frame, NULL);
	 gtk_box_pack_start(GTK_BOX(vbox),GTK_WIDGET(colorsel_frame),FALSE,FALSE,1);
	 gtk_widget_show(colorsel_frame);
	}

	return(frame);
}


GtkWidget *
gtku_object_ctl_panel_trans_check (GTKU_OBJECT *obj)
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
		                 G_CALLBACK(gtku_set_boolean_switch),
		                 (gpointer)sw[ii1]);
		gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(check), FALSE, FALSE, 1);
		gtk_widget_show(check);
	}

	return(frame);
}


GtkWidget *
gtku_object_ctl_panel_color_button (GTKU_OBJECT *obj)
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
		g_object_set_data(G_OBJECT(button), "gtku_object_d_prop_color",
		                  color[ii1]);
		g_signal_connect(G_OBJECT(button), "clicked",
		                 G_CALLBACK(gtku_local_object_color_selection_dialog),
		                 (gpointer)name[ii1]);
		gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(button), FALSE, FALSE, 0);
		gtk_widget_show(button);
	}

	return(frame);
}


gboolean
gtku_local_object_color_selection_dialog (GtkWidget *widget, gpointer cb_data)
{
	GLdouble (*color)[4];
	GtkWidget *cs_widget;
	GtkColorSelectionDialog *dialog;
	GtkColorSelection *colorsel;
	GdkColor current_color;
	guint16 alpha;

	color = (GLdouble (*)[4])g_object_get_data(G_OBJECT(widget),
	                                          "gtku_object_d_prop_color");

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
	                  "gtku_object_d_prop_color", color);
	g_signal_connect(G_OBJECT(dialog->colorsel), "color_changed",
	                 G_CALLBACK(gtku_set_local_object_color),
	                 (gpointer)colorsel);
	/* OK */
	g_object_set_data(G_OBJECT(dialog->ok_button),
	                  "gtku_object_d_prop_color", color);
	g_signal_connect(G_OBJECT(dialog->ok_button), "clicked",
	                 G_CALLBACK(gtku_set_local_object_color),
	                 (gpointer)colorsel);
	g_signal_connect_swapped(G_OBJECT(dialog->ok_button), "clicked",
	                         G_CALLBACK(gtk_widget_destroy), (gpointer)dialog);
	/* Cansel */
	g_object_set_data(G_OBJECT(dialog->cancel_button),
	                  "gtku_object_d_prop_color", color);
	g_object_set_data(G_OBJECT(dialog->cancel_button),
	                  "dialog_window", cs_widget);
	g_signal_connect(G_OBJECT(dialog->cancel_button), "clicked",
	                 G_CALLBACK(gtku_unset_local_object_color),
	                 G_OBJECT(colorsel));

	gtk_widget_show(cs_widget);
	gtk_main();

	return(TRUE);
}


gboolean
gtku_set_local_object_color (GtkWidget *widget, gpointer cb_data)
{
	GLdouble (*glcolor)[4];
	GtkColorSelection *colorsel;
	GdkColor gdkcolor;

	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	glcolor = (GLdouble (*)[4])g_object_get_data(G_OBJECT(widget),
	                                            "gtku_object_d_prop_color");
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
gtku_unset_local_object_color (GtkWidget *widget, gpointer cb_data)
{
	GLdouble (*glcolor)[4];
	GtkColorSelection *colorsel;
	GdkColor gdkcolor;
	GtkWidget *dialog_window;

	/* GTKU_GL existence check */
	g_return_val_if_fail(gtku_gl, FALSE);

	glcolor = (GLdouble (*)[4])g_object_get_data(G_OBJECT(widget),
	                                     "gtku_object_d_prop_color");
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
gtku_object_ctl_panel_draw_num (GTKU_OBJECT *obj)
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
	                  G_CALLBACK(gtku_set_boolean_switch),
	                  (gpointer)&(obj->d_prop->draw_node_num_sw));
	 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(check), FALSE, FALSE, 1);
	 gtk_widget_show(check);
	}

	{ /* Create check button Element Number */
	 GtkWidget *check = gtk_check_button_new_with_label("Element");
	 g_return_val_if_fail(check, NULL);
	 gtk_toggle_button_set_state(GTK_TOGGLE_BUTTON(check), FALSE);
	 g_signal_connect(G_OBJECT(check), "toggled",
	                  G_CALLBACK(gtku_set_boolean_switch),
	                  (gpointer)&(obj->d_prop->draw_elem_num_sw));
	 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(check), FALSE, FALSE, 1);
	 gtk_widget_show(check);
	}

	{ /* Create check button Surface Number */
	 GtkWidget *check = gtk_check_button_new_with_label("Surface");
	 g_return_val_if_fail(check, NULL);
	 gtk_toggle_button_set_state(GTK_TOGGLE_BUTTON(check), FALSE);
	 g_signal_connect(G_OBJECT(check), "toggled",
	                  G_CALLBACK(gtku_set_boolean_switch),
	                  (gpointer)&(obj->d_prop->draw_surf_num_sw));
	 gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(check), FALSE, FALSE, 1);
	 gtk_widget_show(check);
	}

	return(frame);
}


GtkWidget *
gtku_object_ctl_panel_size (GTKU_OBJECT *obj)
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
	 GtkWidget *spin = gtku_local_object_size_spin_with_label(
	                       obj,&(obj->d_prop->point_size),"Point:");
	 g_return_val_if_fail(spin, NULL);
	 gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(spin), FALSE, FALSE, 5);
	 gtk_widget_show(spin);
	}

	{ /* Create Line Size Spin Button */
	 GtkWidget *spin = gtku_local_object_size_spin_with_label(
	                       obj,&(obj->d_prop->line_size),"Line:");
	 g_return_val_if_fail(spin, NULL);
	 gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(spin), FALSE, FALSE, 5);
	 gtk_widget_show(spin);
	}

	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 5);

	return(frame);
}


static GtkWidget *
gtku_local_object_size_spin_with_label (GTKU_OBJECT *obj, GLdouble *size,
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
	                 G_CALLBACK(gtku_set_local_object_size),(gpointer)(size));

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
gtku_set_local_object_size (GtkWidget *widget, gpointer cb_data)
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
gtku_set_boolean_switch (GtkWidget *widget, gpointer cb_data)
{
	/* GTKU_GL existence check */
	gboolean *sw;
	g_return_val_if_fail(gtku_gl, FALSE);
	sw = (gboolean *)cb_data;
	*sw = !(*sw);
	gtk_widget_queue_draw(gtku_gl->glarea);
	return(TRUE);
}


#define PRINT_GTKGLEXT_VISUAL_CONFIG(func) \
{ fprintf(fp, "  "#func" = %s\n", func ? "TRUE" : "FALSE"); }

#define PRINT_GTKGLEXT_CONFIG_ATTRIB(attr, is_boolean) \
{int value;\
    fprintf(fp, "  "#attr" = ");\
    if(gdk_gl_config_get_attrib (glconfig, attr, &value)){\
        if(is_boolean)\
            fprintf(fp, "%s\n", value == TRUE ? "TRUE" : "FALSE");\
        else\
            fprintf(fp, "%d\n", value);\
    } else\
        fprintf(fp, "*** Cannot get "#attr" attribute value\n");\
}

/* for check OpenGL configuration */
void
print_gl_config_attrib (FILE *fp, GdkGLConfig *glconfig)
{
	fprintf(fp, "\nGtkGLExt visual configurations :\n");
	PRINT_GTKGLEXT_VISUAL_CONFIG( gdk_gl_config_is_rgba(glconfig) );
	PRINT_GTKGLEXT_VISUAL_CONFIG( gdk_gl_config_is_double_buffered(glconfig) );
	PRINT_GTKGLEXT_VISUAL_CONFIG( gdk_gl_config_is_stereo(glconfig) );
	PRINT_GTKGLEXT_VISUAL_CONFIG( gdk_gl_config_has_alpha(glconfig) );
	PRINT_GTKGLEXT_VISUAL_CONFIG( gdk_gl_config_has_depth_buffer(glconfig) );
	PRINT_GTKGLEXT_VISUAL_CONFIG( gdk_gl_config_has_stencil_buffer(glconfig) );
	PRINT_GTKGLEXT_VISUAL_CONFIG( gdk_gl_config_has_accum_buffer(glconfig) );
	
	fprintf(fp, "\nGtkGLExt(OpenGL) attribute configurations :\n");
	PRINT_GTKGLEXT_CONFIG_ATTRIB(GDK_GL_USE_GL,           TRUE);
	PRINT_GTKGLEXT_CONFIG_ATTRIB(GDK_GL_BUFFER_SIZE,      FALSE);
	PRINT_GTKGLEXT_CONFIG_ATTRIB(GDK_GL_LEVEL,            FALSE);
	PRINT_GTKGLEXT_CONFIG_ATTRIB(GDK_GL_RGBA,             TRUE);
	PRINT_GTKGLEXT_CONFIG_ATTRIB(GDK_GL_DOUBLEBUFFER,     TRUE);
	PRINT_GTKGLEXT_CONFIG_ATTRIB(GDK_GL_STEREO,           TRUE);
	PRINT_GTKGLEXT_CONFIG_ATTRIB(GDK_GL_AUX_BUFFERS,      FALSE);
	PRINT_GTKGLEXT_CONFIG_ATTRIB(GDK_GL_RED_SIZE,         FALSE);
	PRINT_GTKGLEXT_CONFIG_ATTRIB(GDK_GL_GREEN_SIZE,       FALSE);
	PRINT_GTKGLEXT_CONFIG_ATTRIB(GDK_GL_BLUE_SIZE,        FALSE);
	PRINT_GTKGLEXT_CONFIG_ATTRIB(GDK_GL_ALPHA_SIZE,       FALSE);
	PRINT_GTKGLEXT_CONFIG_ATTRIB(GDK_GL_DEPTH_SIZE,       FALSE);
	PRINT_GTKGLEXT_CONFIG_ATTRIB(GDK_GL_STENCIL_SIZE,     FALSE);
	PRINT_GTKGLEXT_CONFIG_ATTRIB(GDK_GL_ACCUM_RED_SIZE,   FALSE);
	PRINT_GTKGLEXT_CONFIG_ATTRIB(GDK_GL_ACCUM_GREEN_SIZE, FALSE);
	PRINT_GTKGLEXT_CONFIG_ATTRIB(GDK_GL_ACCUM_BLUE_SIZE,  FALSE);
	PRINT_GTKGLEXT_CONFIG_ATTRIB(GDK_GL_ACCUM_ALPHA_SIZE, FALSE);
}


#define PRINT_GTKU_OBJECT_TYPE(type, obj_type) \
{fprintf(fp," [%c] "#obj_type"\n", ((type & obj_type) == obj_type) ? '*': ' ');}
void
print_gtku_object_type (FILE *fp, GTKU_OBJECT_TYPE type)
{
	PRINT_GTKU_OBJECT_TYPE(type, GTKU_OBJECT_FEM_NODE);
	PRINT_GTKU_OBJECT_TYPE(type, GTKU_OBJECT_FEM_ELEM);
	PRINT_GTKU_OBJECT_TYPE(type, GTKU_OBJECT_FEM_SURF);
	PRINT_GTKU_OBJECT_TYPE(type, GTKU_OBJECT_IMAGE_PIXEL);
	PRINT_GTKU_OBJECT_TYPE(type, GTKU_OBJECT_IMAGE_VOXEL);
}


RC
init_gtku_gl ()
{
	/* GTKU_GL existence check */
	RC_NULL_CHK(gtku_gl);

	gtku_gl->object = NULL;
	RC_TRY( allocate_operation_matrix( &(gtku_gl)->rot_mat) );
	RC_TRY( allocate_operation_matrix( &(gtku_gl)->scale_mat) );
	RC_TRY( allocate_operation_matrix( &(gtku_gl)->shift_mat) );
	RC_TRY( allocate_operation_matrix( &(gtku_gl)->initial_mat) );
	gtku_gl->rot_lock_count = 0;
	gtku_gl->scale_lock_count = 0;
	gtku_gl->shift_lock_count = 0;
	gtku_gl->anim_lock_sw = TRUE;
	gtku_gl->animation = FALSE;
	gtku_gl->motion = FALSE;
	gtku_gl->glarea = NULL;
	gtku_gl->glarea_w_size = INIT_WINDOW_SIZE_X;
	gtku_gl->glarea_h_size = INIT_WINDOW_SIZE_X;
	gtku_gl->timeout_id = 0;
	RC_TRY( init_vect3d(&gtku_gl->cur) );
	RC_TRY( init_vect3d(&gtku_gl->pre) );
	RC_TRY( init_vect3d(&gtku_gl->d) );
	gtku_gl->view_size.x = gtku_gl->glarea_w_size;
	gtku_gl->view_size.y = gtku_gl->glarea_h_size;
	gtku_gl->view_size.z = gtku_gl->glarea_w_size*32;
	gtku_gl->draw_type = GTKU_DRAW_POLYGON;
	gtku_gl->draw_coord_sw = TRUE;
	gtku_gl->draw_plane_sw = FALSE;
	gtku_gl->draw_plane_yz_sw = FALSE;
	gtku_gl->draw_plane_zx_sw = FALSE;
	gtku_gl->draw_plane_xy_sw = FALSE;
	gtku_gl->gl_pixmap = NULL;
	gtku_gl->select_start.x = 0;
	gtku_gl->select_start.y = 0;
	gtku_gl->select_area.x = 0;
	gtku_gl->select_area.y = 0;
	gtku_gl->select_area.width = 0;
	gtku_gl->select_area.height = 0;
	gtku_gl->render_mode = GTKU_RENDER_NONE;
	gtku_gl->key.shift = FALSE;
	gtku_gl->key.control = FALSE;
	gtku_gl->key.alt = FALSE;
	gtku_gl->key.x = FALSE;
	gtku_gl->key.y = FALSE;
	gtku_gl->key.z = FALSE;
	gtku_gl->init_light = FALSE;
	gtku_gl->light0_pos[0] = 0.0f;
	gtku_gl->light0_pos[1] = 0.0f;
	gtku_gl->light0_pos[2] = (float)gtku_gl->glarea_w_size*32;
	gtku_gl->light0_pos[3] = 0.0f;
	gtku_gl->light0_color[0] = 0.7f;
	gtku_gl->light0_color[1] = 0.7f;
	gtku_gl->light0_color[2] = 0.7f;
	gtku_gl->light0_color[3] = 1.0f;
	gtku_gl->material_specular[0] = 0.2f;
	gtku_gl->material_specular[1] = 0.2f;
	gtku_gl->material_specular[2] = 0.2f;
	gtku_gl->material_specular[3] = 1.0f;
	gtku_gl->material_shininess[0] = 1.0f;
	gtku_gl->bg_color0[0] = 0.37f;
	gtku_gl->bg_color0[1] = 0.37f;
	gtku_gl->bg_color0[2] = 0.37f;
	gtku_gl->bg_color0[3] = 1.00f;
	gtku_gl->bg_color1[0] = 0.82f;
	gtku_gl->bg_color1[1] = 0.82f;
	gtku_gl->bg_color1[2] = 0.82f;
	gtku_gl->bg_color1[3] = 1.00f;

	gtku_gl->char_list = 0;
	gtku_gl->contour = NULL;
	return(NORMAL_RC);
}
	

void
init_gtku_object_array (GTKU_OBJECT_ARRAY *object)
{
	object->size = 0;
	object->obj = NULL;
	init_vect3d(&object->center);
	init_vect3d(&object->min);
	init_vect3d(&object->max);
}


void
init_gtku_object (GTKU_OBJECT *obj)
{
	obj->title      = NULL;
	obj->file_name  = NULL;
	obj->type       = GTKU_OBJECT_NONE;
	obj->dim        = GTKU_OBJECT_0D;
	obj->d_prop     = NULL;
	obj->node       = NULL;
	obj->elem       = NULL;
	obj->surf       = NULL;
	obj->gr         = NULL;
}


RC
init_gtku_object_prop (GTKU_DRAW_PROP *prop, int obj_num)
{
	int ii1;
	RC_TRY( allocate_operation_matrix( &(prop)->rot_mat) );
	RC_TRY( allocate_operation_matrix( &(prop)->scale_mat) );
	RC_TRY( allocate_operation_matrix( &(prop)->shift_mat) );
	RC_TRY( init_vect3d(&prop->center) );
	RC_TRY( init_vect3d(&prop->min) );
	RC_TRY( init_vect3d(&prop->max) );
	prop->rot_lock_sw = TRUE;
	prop->scale_lock_sw = TRUE;
	prop->shift_lock_sw = TRUE;
	prop->use_local_axis_sw = FALSE;
	prop->cont_type = GRAPHIC_CONT_BGR24;
	prop->quantity = 0.02;
	prop->motion_bbox_sw = TRUE;
	/* FEM Object */
	prop->draw_type = GTKU_DRAW_DEFAULT;
	prop->draw_node_num_sw = FALSE;
	prop->draw_elem_num_sw = FALSE;
	prop->draw_surf_num_sw = FALSE;
	prop->point_trans_sw = FALSE;
	prop->line_trans_sw = FALSE;
	prop->polygon_trans_sw = FALSE;
	prop->text_trans_sw = FALSE;
	if(obj_num < 7){
		for(ii1=0;ii1<4;ii1++){
			prop->point_color[ii1]   = color[obj_num][ii1];
			prop->line_color[ii1]    = color[obj_num][ii1];
			prop->polygon_color[ii1] = color[obj_num][ii1];
			prop->text_color[ii1]    = (color[obj_num][ii1])/1.5;
		}
	}else{
		for(ii1=0;ii1<4;ii1++){
			prop->point_color[ii1]   = color[0][ii1];
			prop->line_color[ii1]    = color[0][ii1];
			prop->polygon_color[ii1] = color[0][ii1];
			prop->text_color[ii1]    = (color[0][ii1])/1.5;
		}
	}
	prop->point_size = 3.0;
	prop->line_size = 1.5;
	prop->texname = 0;
	/* Image Object */
	prop->init_texture = FALSE;
	prop->init_image_dir = FALSE;
	prop->img_trans_sw = TRUE;
	prop->image_dir = GTKU_IMAGE_XY;
	prop->pack_size.x = 0;
	prop->pack_size.y = 0;
	prop->pack_size.z = 0;
	prop->image_depth = 0.0;
	prop->img_threshold_min = 0.8;
	prop->img_threshold_max = 1.0;
	prop->img_alpha = 0.1;

	return(NORMAL_RC);
}

RC
print_gtku_object_prop_operation_mat (FILE *fp, GTKU_DRAW_PROP *d_prop)
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


RC
read_gtku_object_prop_operation_mat (FILE *fp, GTKU_DRAW_PROP *dprop)
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

