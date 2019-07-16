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

/* $Id: gtk_model_view.c 1066 2016-05-19 09:53:33Z hayashi $ */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <glib.h>
#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>
#include <gtk/gtkgl.h>
#include "gtk_utl.h"
#include "base.h"
#include "mathematics.h"
#include "fem.h"


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
static GtkItemFactory *item_factory = NULL;

/* Model view Window */
static GtkWidget *create_glarea(GtkWidget *window);
GdkGLConfig *gtku_get_gdkglconfig(void);
static RC create_ctl_tab(GtkWidget *notebook, GTKU_OBJECT *obj);
static gboolean file_selection_bdf(GtkWidget *widget, gpointer cb_data);
static gboolean filew_ok_bdf(GtkWidget *widget, gpointer cb_data);
static gboolean file_selection_eps(GtkWidget *widget, gpointer cb_data);
static gboolean filew_ok_eps(GtkWidget *widget, gpointer cb_data);
static gboolean save_bdf_file(GtkWidget *widget, gpointer cb_data);
static gboolean save_eps_file(GtkWidget *widget, gpointer cb_data);

/* Menu World Plane */
static void set_world_plane_draw(gpointer cb_data, guint cb_action,
                                 GtkWidget *widget);
static void set_world_plane_yz(gpointer cb_data, guint cb_action,
                               GtkWidget *widget);
static void set_world_plane_zx(gpointer cb_data, guint cb_action,
                               GtkWidget *widget);
static void set_world_plane_xy(gpointer cb_data, guint cb_action,
                               GtkWidget *widget);


/* Menu */
GtkWidget *gtku_create_model_view_menu(GtkWidget *widget);
GtkWidget *gtku_get_model_view_menu(void);

/* Infomation Print */
void print_gl_config_attrib(FILE *fp, GdkGLConfig *glconfig);
void print_gtku_object_type(FILE *fp, GTKU_OBJECT_TYPE type);

/* Initialize */
RC init_gtku_gl(void);
void init_gtku_object_array(GTKU_OBJECT_ARRAY *object);
void init_gtku_object(GTKU_OBJECT *obj);
RC init_gtku_object_prop(GTKU_DRAW_PROP *prop, int obj_num);

/* etc */
static gboolean format_chk(const char *format);
static RC set_element_dim(ELEMENT_ARRAY *elem, GTKU_OBJECT_DIM *dim);


/* 有限要素モデルビューワーを開く．
 * この関数は，main関数内で直接呼び出すことになる．
 * gtk_init() および gtk_gl_init() による初期化を行ったのち，
 * gtku_add_object() で表示するオブジェクトを追加してから呼び出す．
 */
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

	/* トップレベルウィンドウを作成する． */
	window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
	gtk_window_set_title(GTK_WINDOW(window), "gtk_fem");
	gtk_container_set_border_width(GTK_CONTAINER(window), 0);
	gtk_window_set_default_size(GTK_WINDOW(window), 700, 600);
	gtk_quit_add_destroy(1, GTK_OBJECT(window));

	/* ビューワーウィンドウを閉じるときのコールバック関数 */
	g_signal_connect(G_OBJECT(window), "delete_event",
	                 G_CALLBACK(gtku_window_destroy_ask), window);
	g_signal_connect(G_OBJECT(window), "destroy",
	                 G_CALLBACK(gtk_main_quit), NULL);
	gtk_widget_show(GTK_WIDGET(window));


	/* 垂直パッキングボックス（全体） */
	vbox = gtk_vbox_new(FALSE, 2);
	gtk_container_add(GTK_CONTAINER(window), vbox);
	gtk_widget_show(vbox);

	/* メニューバー */
	{ 
	  GtkWidget *menubar;
	  RC_NULL_CHK( menubar = gtku_create_model_view_menu(window) );
	  gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(menubar), FALSE, FALSE, 0);
	  gtk_widget_show(menubar);
	  /*メニューの追加の仕方 
	  GtkItemFactoryEntry item = {"/Option/hoge", NULL, NULL, 0, NULL};
	  gtk_item_factory_create_item(gtk_item_factory_from_widget(menubar),
	                               &item, NULL, 1);
	  */
	}
	
	/* ツールバー */
	/** 水平パッキングボックス **/
	tool_hbox = gtk_hbox_new(FALSE, 0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(tool_hbox), FALSE, FALSE, 0);
	gtk_widget_show(tool_hbox);
	{
	  GtkWidget *toolbar;
	  RC_NULL_CHK( toolbar = gtku_create_toolbar(window) );
	  gtk_box_pack_start(GTK_BOX(tool_hbox), GTK_WIDGET(toolbar),
	                      FALSE, FALSE, 0);
	  gtk_widget_show(toolbar);
	}

	/* 水平ペインウィジェット */
	RC_NULL_CHK( hpaned = gtk_hpaned_new() );
	gtk_container_add(GTK_CONTAINER(vbox), hpaned);
	gtk_widget_show(hpaned);
	/** 左側（OpenGL描画領域） **/
	RC_NULL_CHK( glarea = create_glarea(window) );
	gtk_paned_pack1(GTK_PANED(hpaned), GTK_WIDGET(glarea), FALSE, FALSE);
	gtk_widget_show(glarea);
	/** 右側（コントロールパネル領域） **/
	RC_NULL_CHK( notebook = gtk_notebook_new() );
	gtk_notebook_set_tab_pos(GTK_NOTEBOOK(notebook), GTK_POS_TOP);
	gtk_widget_set_size_request(notebook, 170, -1);
	gtk_paned_pack2(GTK_PANED(hpaned), GTK_WIDGET(notebook), FALSE, TRUE);
	/** Add NoteBook Tab **/
	for(ii1=0; ii1<(gtku_gl->object->size); ii1++){
		RC_TRY( create_ctl_tab(notebook, gtku_gl->object->obj[ii1]) );
	}
	gtk_widget_show(notebook);

	return(NORMAL_RC);
}


/* OpenGLの描画領域を作成し，そのポインタを返す．
 */
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


/* コントロールパネル領域を作成し，そのポインタを返す．
 */
static RC
create_ctl_tab (GtkWidget *notebook, GTKU_OBJECT *obj)
{
	GtkWidget *vbox;
	GtkWidget *scr_window;

	RC_NULL_CHK( notebook );
	RC_NULL_CHK( obj );

	/* スクロールバー付きウィンドウを作成する． */
	scr_window = gtk_scrolled_window_new(NULL, NULL);
	RC_NULL_CHK(scr_window);
	gtk_container_set_border_width(GTK_CONTAINER(scr_window), 0);
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scr_window),
	                               GTK_POLICY_AUTOMATIC,
	                               GTK_POLICY_AUTOMATIC);
	gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW(scr_window),
	                                    GTK_SHADOW_NONE);
	gtk_widget_show(scr_window);

	/* ノートブックを作り，タブにラベルを付ける． */
	{ 
		GtkWidget *tab_label = gtk_label_new(obj->title);
		RC_NULL_CHK( tab_label );
		gtk_widget_show(tab_label);
		gtk_notebook_append_page(GTK_NOTEBOOK(notebook),
			GTK_WIDGET(scr_window), tab_label);
	}

	/* 垂直パッキングボックス（全体） */
	vbox = gtk_vbox_new(FALSE, 0);
	gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scr_window),
	                                      vbox);
	gtk_widget_show(vbox);
	
	/* コントロールパネルを作成する． */
	{ 
		GtkWidget *ctl_panel = gtku_object_ctl_panel(obj);
		RC_NULL_CHK( ctl_panel );
		gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(ctl_panel),
				FALSE, FALSE, 0);
		gtk_widget_show(ctl_panel);
	}

	/* BDFファイルを保存するボタンを作成する． */
	if(obj->type & GTKU_OBJECT_FEM_NODE){
		GtkWidget *button = gtk_button_new_with_label("Save bdf");
		RC_NULL_CHK( button );
		g_signal_connect(G_OBJECT(button), "clicked",
				G_CALLBACK(file_selection_bdf), (gpointer)obj);
		gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(button),
				FALSE, FALSE, 10);
		gtk_widget_show(button);
	}

	/* EPSファイルを保存するボタンを作成する． */
	if( (obj->type & GTKU_OBJECT_FEM_NODE) 
			&& (obj->type & GTKU_OBJECT_FEM_ELEM) ){ 
		GtkWidget *button = gtk_button_new_with_label("Save eps");
		RC_NULL_CHK( button );
		g_signal_connect(G_OBJECT(button), "clicked",
				G_CALLBACK(file_selection_eps), (gpointer)obj);
		gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(button),
				FALSE, FALSE, 10);
		gtk_widget_show(button);
	}

	return(NORMAL_RC);
}


/* BDFファイルを保存するボタンのコールバック関数．
 * ファイル選択ウィンドウを開く．
 */
static gboolean
file_selection_bdf(GtkWidget *widget, gpointer cb_data)
{
	GtkWidget *filew;

	filew = gtk_file_selection_new("Save BDF(NAS) File");
	g_return_val_if_fail(filew, FALSE);
	g_object_set_data(G_OBJECT(filew), "gtku_object", cb_data);

	/* OKボタン */
	g_signal_connect(G_OBJECT(GTK_FILE_SELECTION(filew)->ok_button), "clicked",
	                 G_CALLBACK(filew_ok_bdf), (gpointer)filew);
	/* キャンセルボタン*/
	g_signal_connect_swapped(G_OBJECT(GTK_FILE_SELECTION(filew)->cancel_button),
	                         "clicked", G_CALLBACK(gtk_widget_destroy),
	                         G_OBJECT(filew));
	gtk_widget_show(filew);

	return(TRUE);
}

/* BDFファイル選択ウィンドウのOKボタンのコールバック関数．
 */
static gboolean
filew_ok_bdf(GtkWidget *widget, gpointer cb_data)
{
	if(save_bdf_file(widget, cb_data)) gtk_widget_destroy(GTK_WIDGET(cb_data));
	return(TRUE);
}


/* EPSファイルを保存するボタンのコールバック関数．
 * ファイル選択ウィンドウを開く．
 */
static gboolean
file_selection_eps(GtkWidget *widget, gpointer cb_data)
{
	GtkWidget *filew;

	filew =  gtk_file_selection_new("Save EPS File");
	g_return_val_if_fail(filew, FALSE);
	g_object_set_data(G_OBJECT(filew), "gtku_object", cb_data);
	/* OKボタン */
	g_signal_connect(G_OBJECT(GTK_FILE_SELECTION(filew)->ok_button), "clicked",
	                 G_CALLBACK(filew_ok_eps), (gpointer)filew);
	/* キャンセルボタン */
	g_signal_connect_swapped(G_OBJECT(GTK_FILE_SELECTION(filew)->cancel_button),
	                         "clicked", G_CALLBACK(gtk_widget_destroy),
	                         G_OBJECT(filew));
	gtk_widget_show(filew);

	return(TRUE);
}

/* EPSファイル選択ウィンドウのOKボタンのコールバック関数．
 */
static gboolean
filew_ok_eps(GtkWidget *widget, gpointer cb_data)
{
	if(save_eps_file(widget, cb_data)) gtk_widget_destroy(GTK_WIDGET(cb_data));
	return(TRUE);
}


/* BDFファイルを保存する．
 */
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


/* EPSファイルを保存する．
 */
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


/* obj に可変引数で記述される有限要素情報を格納する．
 * obj に NULL を指定すると，新しい情報として obj を確保する．
 *
 * format
 *   t : title char *
 *   n : NODE_ARRAY
 *   e : ELEMENT_ARRAY
 *   s : ELEMENT_ARRAY (surface)
 *   p : 2D PIXEL
 *   v : 3D VOXEL
 *   a : stress
 *
 * Fix!! 応力，感度などを追加 
 */
static const char *accept_gtku_obj_type = "tnespv";
RC
gtku_add_object (GTKU_OBJECT *obj, const char *format, ...)
{
	int ii1;
	va_list ap;

	/* 引数のフォーマットを確認する． */
	g_return_val_if_fail(format_chk(format), ARG_ERROR_RC);

	/* 描画情報を初期化する． */
 	if(gtku_gl == NULL){
		gtku_gl = (GTKU_GL *)mm_alloc(sizeof(GTKU_GL));
		if(!gtku_gl) return(ALLOC_ERROR_RC);
		RC_TRY( init_gtku_gl() );
	}

	/* create mother of all object(nodem, element, pixel) and initialize */
 	if(gtku_gl->object == NULL){
		gtku_gl->object
			= (GTKU_OBJECT_ARRAY *)mm_alloc(sizeof(GTKU_OBJECT_ARRAY));
		if(!gtku_gl->object) return(ALLOC_ERROR_RC);
		init_gtku_object_array(gtku_gl->object);
	}

	/* 新規登録 */
	if(obj == NULL){
		/* object array を拡大する． */
		gtku_gl->object->size++;
		gtku_gl->object->obj = (GTKU_OBJECT **)mm_realloc(gtku_gl->object->obj,
		                       (gtku_gl->object->size)*sizeof(GTKU_OBJECT *));
		if(gtku_gl->object->obj == NULL) return(ALLOC_ERROR_RC);

		/* gtku_object を確保する． */
		obj = gtku_gl->object->obj[gtku_gl->object->size - 1]
		    = (GTKU_OBJECT *)mm_alloc(sizeof(GTKU_OBJECT));
		if(obj == NULL) return(ALLOC_ERROR_RC);
		init_gtku_object(obj);

		/* gtku_object の d_prop を確保する． */
		obj->d_prop = (GTKU_DRAW_PROP *)mm_alloc(sizeof(GTKU_DRAW_PROP));
		if(obj->d_prop == NULL) return(ALLOC_ERROR_RC);
		init_gtku_object_prop(obj->d_prop, gtku_gl->object->size - 1);
	}

	/* 引数で指定された有限要素データを格納する． */
	va_start(ap, format);
	for(ii1=0; ii1<(int)strlen(format); ii1++){
		switch(format[ii1]){
		/* Title */
		case 't':
			obj->title = obj->file_name = va_arg(ap, char *);
			break;
		/* Node */
		case 'n':
			obj->type |= GTKU_OBJECT_FEM_NODE;
			obj->node = va_arg(ap, NODE_ARRAY *);
			break;
		/* Element */
		case 'e':
			obj->type |= GTKU_OBJECT_FEM_ELEM;
			obj->elem = va_arg(ap, ELEMENT_ARRAY *);
			/* set object(element) dimension */
			RC_TRY( set_element_dim(obj->elem, &(obj->dim)) );
			break;
		/* Element(surf) */
		case 's':
			obj->type |= GTKU_OBJECT_FEM_SURF;
			obj->surf = va_arg(ap, ELEMENT_ARRAY *);
			break;
		/* Stress */
		case 'a':
			obj->type |= GTKU_OBJECT_FEM_STRESS;
			obj->stress = va_arg(ap, STRESS_ARRAY *);
			break;
		/* Pixel */
		case 'p':
			obj->type |= GTKU_OBJECT_IMAGE_PIXEL;
			obj->gr = va_arg(ap, GRAPHIC *);
			obj->dim = GTKU_OBJECT_2D;
			break;
		/* Voxel */
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

	/* obj の title と filename がなければ，適当に決める． */
	if(obj->title == NULL){
		int length = sizeof(char)*32;

		/* obj->title */
		obj->title = (char *)mm_alloc(length*sizeof(char));
		if(obj->title == NULL) return(ALLOC_ERROR_RC);
		RC_NEG_ZERO_CHK( g_snprintf(obj->title, length-1, "object%d",
		                            gtku_gl->object->size - 1) );
		obj->title[length-1] = (char)'\0';  /* 念の為 */

		/* obj->filename */
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


/* gtku_add_object() に渡す書式が適切かどうかを調べる．
 */
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
				if((mask ^ (1<<ii2)) == 0){ /* same format type ERROR!! */
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


/* dim に要素の次元を格納する．
 */
static RC
set_element_dim (ELEMENT_ARRAY *elem, GTKU_OBJECT_DIM *dim)
{
	RC_NULL_CHK( elem );
	RC_NULL_CHK( dim );

	switch(analysis_dim(*elem)){
		case 1:
			*dim = GTKU_OBJECT_1D;
			break;
		case 2:
			*dim = GTKU_OBJECT_2D;
			break;
		case 3:
			*dim = GTKU_OBJECT_3D;
			break;
		default:
			return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


/* OpenGLの描画領域ウィジェットを作成し，そのポインタを返す．
 */
GtkWidget *
gtku_glarea_new (GtkWidget *window)
{
	GdkGLConfig *glconfig;
	
	/* get information about a OpenGL frame buffer configuration */
	/* 初期化する． */
	glconfig = gtku_get_gdkglconfig();
	g_return_val_if_fail(glconfig, NULL);

#ifdef GTKU_VERBOSE_MESSAGE
	/* OpenGLの属性値を表示する． */
	print_gl_config_attrib(stderr, glconfig);

	/* OpenGLのバージョンを表示する． */
	{ 
		gint major, minor;
		gdk_gl_query_version(&major, &minor);
		fprintf(stderr, "\nOpenGL extension version : <%d.%d>\n", major, minor);
	}
#endif /* GTKU_VERBOSE_MESSAGE */

	/* 描画領域ウィジェット */
	gtku_gl->glarea = GTK_WIDGET(gtk_drawing_area_new());
	g_return_val_if_fail(gtku_gl->glarea, NULL);
	gtk_widget_set_size_request(GTK_WIDGET(gtku_gl->glarea),
	                            gtku_gl->glarea_w_size,
	                            gtku_gl->glarea_h_size);

	/* 描画領域ウィジェットにGtkGLExtの設定を行う． */
	gtk_widget_set_gl_capability(GTK_WIDGET(gtku_gl->glarea), glconfig,
	                             NULL, TRUE, GDK_GL_RGBA_TYPE);
	g_object_unref(G_OBJECT(glconfig));

	gtk_widget_set_events(GTK_WIDGET(gtku_gl->glarea),
	                      GDK_EXPOSURE_MASK |
	                      GDK_BUTTON_PRESS_MASK |
	                      GDK_BUTTON_RELEASE_MASK |
	                      GDK_POINTER_MOTION_MASK |
	                      GDK_POINTER_MOTION_HINT_MASK);

	/* シグナルはここに書かれた順番に試される．特にbutton_press */
	/* motion_notifyの順番には注意                              */
	/** 描画タイミング **/
	g_signal_connect(G_OBJECT(gtku_gl->glarea), "expose_event",
	                 G_CALLBACK(gtku_gl_draw), NULL);
	g_signal_connect(G_OBJECT(gtku_gl->glarea), "configure_event",
	                 G_CALLBACK(gtku_gl_reshape), NULL);
	/** Xウィンドウが作成されたとき（一度だけ） **/
	g_signal_connect(G_OBJECT(gtku_gl->glarea), "realize",
	                 G_CALLBACK(gtku_gl_realize), NULL);

	/** マウスボタンが押されたとき **/
	g_signal_connect(G_OBJECT(gtku_gl->glarea), "button_press_event",
	                 G_CALLBACK(gtku_gl_button_press), NULL);
	/** マウスボタンが離されたとき **/
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


/* glconfig を初期化する．
 * 描画はダブルバッファで行う．
 */
GdkGLConfig *
gtku_get_gdkglconfig (void)
{
	GdkGLConfig *glconfig;
	/* attribute listのTRUE/FALSEは関係ない．listに有るものがONとなる    */
	/* 例えば，GDK_GL_DOUBLEBUFFERを無効にしたければリストから外すこと． */ 
	const int attr_list[] = {
		GDK_GL_RGBA,         TRUE,
		GDK_GL_DOUBLEBUFFER, TRUE,
		GDK_GL_DEPTH_SIZE,   16,
		GDK_GL_ATTRIB_LIST_NONE
	};

	/* OpenGLがサポートされているかを確認する． */
	if( !gdk_gl_query_extension() ) g_error("*** OpenGL is not supported.\n");

	/* ダブルバッファが利用できるかを確認する． */
	glconfig = gdk_gl_config_new(attr_list);
	if(glconfig == NULL){
		g_warning("*** Cannot find the double-buffered visual.\n");
		g_warning("*** Trying single-buffered visual.\n");

		/* シングルバッファが利用できるかを確認する． */
		glconfig = gdk_gl_config_new_by_mode(GDK_GL_MODE_RGB |
		                                     GDK_GL_MODE_DEPTH |
		                                     GDK_GL_MODE_SINGLE);
		if(glconfig == NULL)
			g_error("*** No appropriate OpenGL-capable visual found.\n");
	}

	return(glconfig);
}


/* 座標平面を描画する．
 */
static void
set_world_plane_draw(gpointer cb_data, guint cb_action, GtkWidget *widget)
{
	gtku_gl->draw_plane_sw = !(gtku_gl->draw_plane_sw);
	if(gtku_gl->glarea) gtk_widget_queue_draw(gtku_gl->glarea);
}


/* YZ平面を描画する．
 */
static void
set_world_plane_yz(gpointer cb_data, guint cb_action, GtkWidget *widget)
{
	gtku_gl->draw_plane_yz_sw = !(gtku_gl->draw_plane_yz_sw);
	if(gtku_gl->glarea) gtk_widget_queue_draw(gtku_gl->glarea);
}


/* ZX平面を描画する．
 */
static void
set_world_plane_zx(gpointer cb_data, guint cb_action, GtkWidget *widget)
{
	gtku_gl->draw_plane_xy_sw = !(gtku_gl->draw_plane_xy_sw);
	if(gtku_gl->glarea) gtk_widget_queue_draw(gtku_gl->glarea);
}


/* XY平面を描画する．
 */
static void
set_world_plane_xy(gpointer cb_data, guint cb_action, GtkWidget *widget)
{
	gtku_gl->draw_plane_xy_sw = !(gtku_gl->draw_plane_xy_sw);
	if(gtku_gl->glarea) gtk_widget_queue_draw(gtku_gl->glarea);
}


static GtkItemFactoryEntry menu_items[] = {
 /* File ->  */
 {("/_File"), NULL, NULL, 0, "<Branch>", NULL},
 { "/File/tear1", NULL, NULL, 0, "<Tearoff>", NULL },
 {("/File/_Quit"),"<control>Q",gtk_main_quit, 0, "<StockItem>", GTK_STOCK_QUIT},
 /* View ->  */
 {("/View"),      NULL, NULL, 0, "<Branch>", NULL},
 { "/View/tear2", NULL, NULL, 0, "<Tearoff>", NULL },
 /* View -> World Plane */
 {("/View/World Plane/Draw Plane"), NULL, set_world_plane_draw, 0, "<CheckItem>", NULL},
 {"/View/World Plane/separator", NULL, NULL, 0, "<Separator>", NULL},
 {("/View/World Plane/YZ-Plane"), NULL, set_world_plane_yz, 0, "<CheckItem>", NULL},
 {("/View/World Plane/ZX-Plane"), NULL, set_world_plane_zx, 0, "<CheckItem>", NULL},
 {("/View/World Plane/XY-Plane"), NULL, set_world_plane_xy, 0, "<CheckItem>", NULL}
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
	PRINT_GTKU_OBJECT_TYPE(type, GTKU_OBJECT_FEM_STRESS);
	PRINT_GTKU_OBJECT_TYPE(type, GTKU_OBJECT_IMAGE_PIXEL);
	PRINT_GTKU_OBJECT_TYPE(type, GTKU_OBJECT_IMAGE_VOXEL);
}


/* gtku_gl を初期化する．
 */
RC
init_gtku_gl (void)
{
	/* 確保済みかどうかを調べる． */
	RC_NULL_CHK( gtku_gl );

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


