/*********************************************************************
 * gtk_utl.h
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

/* $Id: gtk_utl.h 450 2005-09-12 05:07:22Z sasaoka $ */

#ifndef GTK_UTL_H
#define GTK_UTL_H

#if WIN32
#  include <windows.h> 
#endif

#include <GL/gl.h>
#include <GL/glu.h>
#include <glib.h>
#include <gtk/gtk.h>
#include <gtk/gtkgl.h>
#include "fem_struct.h"
#include "graphic_utl.h"

#define GTKU_VERBOSE_MESSAGE
#define GTKU_VERBOSE_CHECK

typedef enum{
	GTKU_OBJECT_NONE         = 0 << 0,
	GTKU_OBJECT_FEM_NODE     = 1 << 0,
	GTKU_OBJECT_FEM_ELEM     = 1 << 1,
	GTKU_OBJECT_FEM_SURF     = 1 << 2,
	GTKU_OBJECT_IMAGE_PIXEL  = 1 << 3,  
	GTKU_OBJECT_IMAGE_VOXEL  = 1 << 4   
} GTKU_OBJECT_TYPE;


typedef enum{
	GTKU_DRAW_DEFAULT    = 1 << 0, /* same global draw type  */
	GTKU_DRAW_HIDE       = 1 << 1, /* don't draw this object */
	GTKU_DRAW_BBOX       = 1 << 2,
	GTKU_DRAW_POLYGON    = 1 << 3,
	GTKU_DRAW_SOLID      = 1 << 4,
	GTKU_DRAW_WIRE_FRAME = 1 << 5,
	GTKU_DRAW_GRID_FRAME = 1 << 6,
	GTKU_DRAW_SURF_POINT = 1 << 7,
	GTKU_DRAW_GRID_POINT = 1 << 8,
	GTKU_DRAW_HIDDENLINE = 1 << 9,
	GTKU_DRAW_TEXTURE    = 1 << 10,
	GTKU_DRAW_CELL       = 1 << 11
} GTKU_DRAW_TYPE;


typedef enum{
	GTKU_RENDER_NONE,
	GTKU_RENDER_GL,
	GTKU_RENDER_AREA,
	GTKU_RENDER_SELECT,
	GTKU_RENDER_PIXMAP
} GTKU_RENDER_MODE;


typedef enum{
	GTKU_OBJECT_0D,
	GTKU_OBJECT_1D,
	GTKU_OBJECT_2D,
	GTKU_OBJECT_3D
} GTKU_OBJECT_DIM;

typedef enum{
	GTKU_IMAGE_XY,
	GTKU_IMAGE_YZ,
	GTKU_IMAGE_ZX
} GTKU_IMAGE_DIR;

typedef struct{
	double **rot_mat;                /* rotation matrix */ 
	double **scale_mat;              /* scale matrix */ 
	double **shift_mat;              /* shift matrix */ 
	VECT3D center;                   /* object center point */
	VECT3D min, max;                 /* object maximum, minimum point */
	gboolean rot_lock_sw;            /* rotation lock switch */
	gboolean scale_lock_sw;          /* scale lock switch */
	gboolean shift_lock_sw;          /* shift lock switch */
	gboolean use_local_axis_sw;      /* use local axis switch */
	GRAPHIC_CONTOUR_TYPE cont_type;  /* color contour type */
	double quantity;                 /* display quantity */
	gboolean motion_bbox_sw;         /* hide object instead of draw bbox) */
	/* FEM Object */
	GTKU_DRAW_TYPE draw_type;   /* drawing type of model */
	gboolean draw_node_num_sw;  /* draw node number */
	gboolean draw_elem_num_sw;  /* draw element number */
	gboolean draw_surf_num_sw;  /* draw surface number */
	gboolean point_trans_sw;    /* translucent point */
	gboolean line_trans_sw;     /* translucent line */
	gboolean polygon_trans_sw;  /* translucent polygon */
	gboolean text_trans_sw;     /* translucent text */
	GLdouble point_color[4];    /* point color */
	GLdouble line_color[4];     /* line color */
	GLdouble polygon_color[4];  /* polygon color */
	GLdouble text_color[4];     /* text color */
	GLdouble point_size;        /* point size */
	GLdouble line_size;         /* line size */
	GLuint texname;             /* Texture Name for GRAPHIC */
	/* Image Object */
	gboolean init_texture;      /* initialize texture setting */
	gboolean init_image_dir;    /* initialize image(pix,vox) direction */
	gboolean img_trans_sw;      /* translucent image */
	GTKU_IMAGE_DIR image_dir;   /* image view direction */
	GRAPHIC_SIZE pack_size;     /* texture packing size */
	double image_depth;         /* image view depth value */
	double img_threshold_min;   /* display range threshold */
	double img_threshold_max;   /* display range threshold */
	double img_alpha;           /* image translucent value */
} GTKU_DRAW_PROP;

typedef struct{
	char *title;              /* object name (default: NULL) */
	char *file_name;          /* file name   (default: NULL) */
	GTKU_OBJECT_TYPE type;    /* object type */
	GTKU_OBJECT_DIM dim;      /* object dimesion */
	GTKU_DRAW_PROP *d_prop;   /* object drawing property */
	NODE_ARRAY *node;     /* fem node array object not array */
	ELEMENT_ARRAY *elem;  /* fem element array object not array */
	ELEMENT_ARRAY *surf;  /* fem surface element array object not array */
	GRAPHIC *gr;              /* 2D-Pixel 3D-VOXEL object not array */
} GTKU_OBJECT;

typedef struct{
	int size;           /* number of  objects */
	GTKU_OBJECT **obj;  /* object array */
	VECT3D center;      /* all object center */
	VECT3D min, max;    /* all object maximum, minimum point */
} GTKU_OBJECT_ARRAY;

typedef struct{
	gboolean shift;     /* If press any key, value become TRUE, */
	gboolean control;   /* and  release key became FALSE        */
	gboolean alt;
	gboolean x;
	gboolean y;
	gboolean z;
} GTKU_KEY_MAP;

#if 0
typedef struct{
	GtkWidget *glarea;         /* OpenGL Rendering Area Widget */
	int glarea_w_size;         /* glarea window width size */
	int glarea_h_size;         /* glarea window height size */
	VECT3D cur, pre, d;        /* curent & pre-mouse potition and diff */
	VECT3D view_size;          /* view size for OpenGL func "glOrtho" */
	gboolean init_light;           /* true if initgl not yet called */
	GLfloat light0_pos[4];         /* lighting position */
	GLfloat light0_color[4];       /* light color */
	GLfloat material_specular[4];  /* material specular */
	GLfloat material_shininess[1]; /* material shininess */
	GLdouble bg_color[4];   
} GTKU_GLAREA;

typedef struct{
	int size;              /* number of  objects */
	GTKU_GLAREA **glarea;  /* object array */
} GTKU_GLAREA_ARRAY;
#endif

typedef struct{
	GTKU_OBJECT_ARRAY *object; /* object array */
	double **rot_mat;          /* global rotation matrix */ 
	double **scale_mat;        /* global scale matrix */ 
	double **shift_mat;        /* global shift matrix */ 
	double **initial_mat;      /* global initial matrix */ 
	int rot_lock_count;        /* global rotation lock count */
	int scale_lock_count;      /* global scale lock count */
	int shift_lock_count;      /* global scale lock count */
	gboolean anim_lock_sw;     /* global animation lock switch */
	gboolean animation;        /* animetion toggle switch */
	gboolean motion;           /* motion state switch */
	GtkWidget *glarea;         /* OpenGL Rendering Area Widget */
	int glarea_w_size;         /* glarea window width size */
	int glarea_h_size;         /* glarea window height size */
	unsigned int timeout_id;   /* animation timeout ID */
	VECT3D cur, pre, d;        /* curent & pre-mouse potition and diff */
	VECT3D view_size;          /* view size for OpenGL func "glOrtho" */
	GTKU_DRAW_TYPE draw_type;  /* global drawing type of model */
	gboolean draw_coord_sw;    /* draw coordinate system switch */
	gboolean draw_plane_sw;    /* draw world plane switch */
	gboolean draw_plane_yz_sw;  /* draw xy world plane switch */
	gboolean draw_plane_zx_sw;  /* draw yz world plane switch */
	gboolean draw_plane_xy_sw;  /* draw yz world plane switch */
	gboolean draw_center_mark_sw;  /* draw object center marker switch */
	GdkPixmap *gl_pixmap;      /* glarea off-screen rendering pixmap */
	GdkPoint select_start;     /* selection area */
	GdkRectangle select_area;  /* selection area */
	GTKU_RENDER_MODE render_mode; /* GL rendering mode */
	GTKU_KEY_MAP key;          /* event key map */
	gboolean init_light;           /* initialize lighting */
	GLfloat light0_pos[4];         /* lighting position */
	GLfloat light0_color[4];       /* light color */
	GLfloat material_specular[4];  /* material specular */
	GLfloat material_shininess[1]; /* material shininess */
	GLfloat bg_color0[4];          /* backgroud color (top)*/ 
	GLfloat bg_color1[4];          /* backgroud color (bottom)*/ 
	GLuint char_list;          /* disply-list of ASCII characters */
	GRAPHIC_CONTOUR *contour;  /* color contour */
} GTKU_GL;

/* gtk_model_view.c */
RC gtk_model_view(void);
RC gtku_add_object(GTKU_OBJECT *obj, const char *format, ...);
GtkWidget *gtku_glarea_new(GtkWidget *window);
GdkGLConfig *gtku_get_gdkglconfig(void);

/* gtk_model_view.c : Toolbar Rotation */
gboolean gtku_toolbar_rotation(GtkWidget *window, GtkWidget *toolbar);
GtkWidget *gtku_toolbar_rotation_xy(GtkWidget *widget, GtkWidget *toolbar,
                                    GtkSignalFunc func, gpointer data);
GtkWidget *gtku_toolbar_rotation_yz(GtkWidget *widget, GtkWidget *toolbar,
                                    GtkSignalFunc func, gpointer data);
GtkWidget *gtku_toolbar_rotation_zx(GtkWidget *widget, GtkWidget *toolbar,
                                    GtkSignalFunc func, gpointer data);
GtkWidget *gtku_toolbar_rotation_xyz(GtkWidget *widget, GtkWidget *toolbar,
                                     GtkSignalFunc func, gpointer data);
gboolean gtku_set_model_view_direction(GtkWidget *menu,gpointer cb_data);

/* gtk_model_view.c : Toolbar Drawwing Type */
gboolean gtku_toolbar_drawtype(GtkWidget *window, GtkWidget *toolbar);
GtkWidget *gtku_toolbar_drawtype_polygon(GtkWidget *widget, GtkWidget *toolbar,
                                         GtkWidget *group, GtkSignalFunc func,
                                         gpointer data);
GtkWidget *gtku_toolbar_drawtype_solid(GtkWidget *widget, GtkWidget *toolbar,
                                         GtkWidget *group, GtkSignalFunc func,
                                         gpointer data);
GtkWidget *gtku_toolbar_drawtype_gridframe(GtkWidget *widget,
                                           GtkWidget *group, GtkWidget *toolbar,
                                           GtkSignalFunc func, gpointer data);
GtkWidget *gtku_toolbar_drawtype_wireframe(GtkWidget *widget,
                                           GtkWidget *group, GtkWidget *toolbar,
                                           GtkSignalFunc func, gpointer data);
GtkWidget *gtku_toolbar_drawtype_hiddenline(GtkWidget *widget,
                                            GtkWidget *toolbar,
                                            GtkWidget *group,
                                            GtkSignalFunc func, gpointer data);
GtkWidget *gtku_toolbar_drawtype_gridpoint(GtkWidget *widget,GtkWidget *toolbar,
                                           GtkWidget *group, GtkSignalFunc func,
                                           gpointer data);
GtkWidget *gtku_toolbar_drawtype_surfpoint(GtkWidget *widget,GtkWidget *toolbar,
                                           GtkWidget *group, GtkSignalFunc func,
                                           gpointer data);
gboolean gtku_set_model_draw_type(GtkWidget *widget, gpointer cb_data);

/* gtk_model_view.c : Toolbar Clear */
gboolean gtku_toolbar_clear(GtkWidget *widget, GtkWidget *toolbar);
gboolean gtku_set_model_clear(GtkWidget *widget, gpointer cb_data);

/* gtk_model_view.c : Toolbar Scale */
gboolean gtku_toolbar_scale(GtkWidget *window, GtkWidget *toolbar);
gboolean gtku_set_model_scale(GtkWidget *widget, gpointer cb_data);

/* gtk_model_view.c : Toolbar Animation */
GtkWidget *gtku_toolbar_animation(GtkWidget *widget, GtkWidget *toolbar);
gboolean gtku_set_model_animation(GtkWidget *widget, gpointer cb_data);

/* gtk_model_view.c : Menu World Plane */
gboolean gtku_set_world_plane(GtkWidget *widget, gpointer cb_data);

/* gtk_model_view.c : Menu */
GtkWidget *gtku_create_model_view_menu(GtkWidget *widget);
GtkWidget *gtku_get_model_view_menu(void);

/* gtk_model_view.c : Control Panel */
GtkWidget *gtku_object_ctl_panel(GTKU_OBJECT *obj);

/* gtk_model_view.c : Control Panel Matrix Lock */
GtkWidget *gtku_object_ctl_panel_matrix(GTKU_OBJECT *obj);
gboolean gtku_set_local_rot_matrix_lock(GtkWidget *widget, gpointer cb_data);
gboolean gtku_set_local_scale_matrix_lock(GtkWidget *widget, gpointer cb_data);
gboolean gtku_set_local_shift_matrix_lock(GtkWidget *widget, gpointer cb_data);

/* gtk_model_view.c : Control Panel Drawing Type */
GtkWidget *gtku_object_ctl_panel_draw_type(GTKU_OBJECT *obj);
gboolean gtku_set_local_object_draw_type(GtkWidget *widget, gpointer cb_data);

/* gtk_model_view.c : Control Panel Quantity of Display */
GtkWidget *gtku_object_ctl_panel_quantity(GTKU_OBJECT *obj);
gboolean gtku_set_local_quantity(GtkWidget *widget, gpointer cb_data);

/* gtk_model_view.c : Control Panel Image Direction */
GtkWidget *gtku_object_ctl_panel_image_dir(GTKU_OBJECT *obj);
GtkWidget *gtku_local_image_direction_option(GTKU_OBJECT *obj);
gboolean gtku_set_local_image_direction(GtkWidget *widget, gpointer cb_data);

/* gtk_model_view.c : Control Panel Image Contour */
GtkWidget *gtku_object_ctl_panel_image_contour(GTKU_OBJECT *obj);
GtkWidget *gtku_local_image_contour_option(GTKU_OBJECT *obj);
gboolean gtku_set_local_image_contour(GtkWidget *widget, gpointer cb_data);


/* gtk_model_view.c : Control Panel Image Display Thrshold */
GtkWidget *gtku_object_ctl_panel_image_range(GTKU_OBJECT *obj);
gboolean gtku_set_local_image_range(GtkWidget *widget, gpointer cb_data);

/* gtk_model_view.c : Control Panel TransLucent of Image */
GtkWidget *gtku_object_ctl_panel_image_trans(GTKU_OBJECT *obj);
gboolean gtku_set_local_image_trans(GtkWidget *widget, gpointer cb_data);

/* gtk_model_view.c : Control Panel Color&Alpha Selection */
GtkWidget *gtku_object_ctl_panel_color(GTKU_OBJECT *obj);
GtkWidget *gtku_object_ctl_panel_trans_check(GTKU_OBJECT *obj);
GtkWidget *gtku_object_ctl_panel_color_button(GTKU_OBJECT *obj);
gboolean gtku_local_object_color_selection_dialog(GtkWidget *widget,
                                                  gpointer cb_data);
gboolean gtku_set_local_object_color(GtkWidget *widget, gpointer cb_data);
gboolean gtku_unset_local_object_color(GtkWidget *widget, gpointer cb_data);

/* gtk_model_view.c : Control Panel Drawing Number */
GtkWidget *gtku_object_ctl_panel_draw_num(GTKU_OBJECT *obj);

/* gtk_model_view.c : Control Panel Size (Point, Line) */
GtkWidget *gtku_object_ctl_panel_size(GTKU_OBJECT *obj);
gboolean gtku_set_local_object_size(GtkWidget *widget, gpointer cb_data);

gboolean gtku_set_boolean_switch(GtkWidget *widget, gpointer cb_data);

/* gtk_model_view.c : Infomation Print */
void print_gl_config_attrib(FILE *fp, GdkGLConfig *glconfig);
void print_gtku_object_type(FILE *fp, GTKU_OBJECT_TYPE type);
RC print_gtku_object_prop_operation_mat(FILE *fp, GTKU_DRAW_PROP *d_prop);

/* gtk_model_view.c : Initialize */
RC init_gtku_gl(void);
void init_gtku_object_array(GTKU_OBJECT_ARRAY *object);
void init_gtku_object(GTKU_OBJECT *obj);
RC init_gtku_object_prop(GTKU_DRAW_PROP *prop, int obj_num);

/* gtk_model_view.c : Infomation Read */
RC read_gtku_object_prop_operation_mat(FILE *fp, GTKU_DRAW_PROP *dprop);


/* gtk_graph.c */
GtkWidget *gtku_drow_graph(void);


/* gtk_utl.c */
gboolean gtku_window_destory(GtkWidget *widget);
gboolean gtku_window_destory_ask(GtkWidget *widget, GdkEvent *event,
                                 GtkWidget *parent);
GtkWidget *gtku_create_gtkimage_from_xpm(GtkWidget *widget, gchar **xpm);
void create_gtk_optionmenu_item(GtkWidget *menu, GSList *group,
                                GtkSignalFunc func, const gchar *menu_text,
                                gpointer cb_data);


/* gtkgl_utl.c */
gboolean gtku_gl_draw(GtkWidget *widget, GdkEventExpose *event);
gboolean gtku_gl_reshape(GtkWidget *widget, GdkEventConfigure *event);
gboolean gtku_gl_realize(GtkWidget *widget);
gboolean gtku_gl_motion_notify(GtkWidget *widget, GdkEventMotion *event);
gboolean gtku_gl_button_press(GtkWidget *widget, GdkEventButton *event);
gboolean gtku_gl_button_release(GtkWidget *widget, GdkEventButton *event);
gboolean gtku_gl_map(GtkWidget *widget, GdkEvent *event);
gboolean gtku_gl_unmap(GtkWidget *widget, GdkEvent *event);
gboolean gtku_gl_visibility_notify(GtkWidget *widget,GdkEventVisibility *event);
gboolean gtku_gl_destroy(GtkWidget *widget);
gboolean gtku_gl_key_press(GtkWidget *widget, GdkEventKey *event);
gboolean gtku_gl_key_release(GtkWidget *widget, GdkEventKey *event);
gboolean gtku_gl_scroll(GtkWidget *widget, GdkEventScroll *event);
gboolean gtku_gl_toggle_animation(GtkWidget *widget);
void set_gl_translucent(void);
void unset_gl_translucent(void);
void gtku_gl_draw_string(VECT3D p, const char *str);
void gtku_gl_conv_form_matrix(double **src, double dest[4][4]);
gboolean gtku_object_exist_chk(GTKU_OBJECT *obj, int size, ...);

/* gtkgl_article.c */
void gtku_gl_draw_coodinate_axis(void);
void gtku_gl_draw_center_marker(VECT3D center, double length);
RC gtku_gl_draw_world_plane(void);
RC gtku_gl_draw_bbox(GTKU_DRAW_PROP *prop);

/* gtkgl_draw_fem.c */
RC gtku_gl_draw_fem_object(GTKU_OBJECT *obj);

/* gtkgl_draw_image.c */
RC gtku_gl_draw_image_object(GTKU_OBJECT *obj);


#define GTKU_BOOLEAN_RC_TRY(func) \
{RC lrc = func;\
 if(lrc != NORMAL_RC){\
     rc_error_print(__FILE__, __LINE__, #func, lrc);\
     return(FALSE);\
}}

#define GTKU_VOID_RC_TRY(func) \
{RC lrc = func;\
 if(lrc != NORMAL_RC){\
     rc_error_print(__FILE__, __LINE__, #func, lrc);\
     g_critical("*** CRITICAL ERROR !! ***\n"); \
}}

#define GTKU_WARNING(func, message) \
{ g_warning("\n[%s : %d] %s\n    => %s\n", __FILE__, __LINE__,func, message); }

#endif /* GTK_UTL_H */

