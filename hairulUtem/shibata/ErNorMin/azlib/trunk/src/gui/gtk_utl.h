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

/* $Id: gtk_utl.h 1066 2016-05-19 09:53:33Z hayashi $ */

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
	GTKU_OBJECT_NONE         = 0 << 0, /* 未定 */
	GTKU_OBJECT_FEM_NODE     = 1 << 0, /* 節点 */
	GTKU_OBJECT_FEM_ELEM     = 1 << 1, /* 要素 */
	GTKU_OBJECT_FEM_SURF     = 1 << 2, /* 要素（表面） */
	GTKU_OBJECT_FEM_STRESS   = 1 << 3, /* 応力 */
	GTKU_OBJECT_IMAGE_PIXEL  = 1 << 4, /* ピクセル */
	GTKU_OBJECT_IMAGE_VOXEL  = 1 << 5  /* ボクセル */
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


/* 描画モード */
typedef enum{
	GTKU_RENDER_NONE,   /* 描画しない */
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
	char *title;            /* object name (default: NULL) */
	char *file_name;        /* file name   (default: NULL) */
	GTKU_OBJECT_TYPE type;  /* object type */
	GTKU_OBJECT_DIM dim;    /* object dimesion */
	GTKU_DRAW_PROP *d_prop; /* object drawing property */
	NODE_ARRAY *node;       /* fem node array object not array */
	ELEMENT_ARRAY *elem;    /* fem element array object not array */
	ELEMENT_ARRAY *surf;    /* fem surface element array object not array */
	STRESS_ARRAY *stress;   /* stress array object not array */
	GRAPHIC *gr;            /* 2D-Pixel 3D-VOXEL object not array */
} GTKU_OBJECT;

typedef struct{
	int size;          /* number of  objects */
	GTKU_OBJECT **obj; /* object array */
	VECT3D center;     /* all object center */
	VECT3D min, max;   /* all object maximum, minimum point */
} GTKU_OBJECT_ARRAY;

/* キーマッピング
 * 各キーが押されている間のみ，TRUE になる． */
typedef struct{
	gboolean shift;
	gboolean control;
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
	GTKU_OBJECT_ARRAY *object;     /* object array */
	double **rot_mat;              /* global rotation matrix */ 
	double **scale_mat;            /* global scale matrix */ 
	double **shift_mat;            /* global shift matrix */ 
	double **initial_mat;          /* global initial matrix */ 
	int rot_lock_count;            /* global rotation lock count */
	int scale_lock_count;          /* global scale lock count */
	int shift_lock_count;          /* global scale lock count */
	gboolean anim_lock_sw;         /* global animation lock switch */
	gboolean animation;            /* animation toggle switch */
	gboolean motion;               /* motion state switch */
	GtkWidget *glarea;             /* OpenGL Rendering Area Widget */
	int glarea_w_size;             /* glarea window width size */
	int glarea_h_size;             /* glarea window height size */
	unsigned int timeout_id;       /* animation timeout ID */
	VECT3D cur, pre, d;            /* curent & pre-mouse potition and diff */
	VECT3D view_size;              /* view size for OpenGL func "glOrtho" */
	GTKU_DRAW_TYPE draw_type;      /* global drawing type of model */
	gboolean draw_coord_sw;        /* draw coordinate system switch */
	gboolean draw_plane_sw;        /* draw world plane switch */
	gboolean draw_plane_yz_sw;     /* draw xy world plane switch */
	gboolean draw_plane_zx_sw;     /* draw yz world plane switch */
	gboolean draw_plane_xy_sw;     /* draw yz world plane switch */
	gboolean draw_center_mark_sw;  /* draw object center marker switch */
	GdkPixmap *gl_pixmap;          /* glarea off-screen rendering pixmap */
	GdkPoint select_start;         /* selection area */
	GdkRectangle select_area;      /* selection area */
	GTKU_RENDER_MODE render_mode;  /* GL rendering mode */
	GTKU_KEY_MAP key;              /* event key map */
	gboolean init_light;           /* initialize lighting */
	GLfloat light0_pos[4];         /* lighting position */
	GLfloat light0_color[4];       /* light color */
	GLfloat material_specular[4];  /* material specular */
	GLfloat material_shininess[1]; /* material shininess */
	GLfloat bg_color0[4];          /* backgroud color (top)*/ 
	GLfloat bg_color1[4];          /* backgroud color (bottom)*/ 
	GLuint char_list;              /* disply-list of ASCII characters */
	GRAPHIC_CONTOUR *contour;      /* color contour */
} GTKU_GL;

/* gtk_model_view.c */
RC gtk_model_view(void);
RC gtku_add_object(GTKU_OBJECT *obj, const char *format, ...);
GtkWidget *gtku_glarea_new(GtkWidget *window);

/* gtk_toolbar.c */
GtkWidget *gtku_create_toolbar(GtkWidget *window);

/* gtk_control_panel.c */
GtkWidget *gtku_object_ctl_panel(GTKU_OBJECT *obj);

/* gtk_utl.c */
gboolean gtku_window_destroy(GtkWidget *widget);
gboolean gtku_window_destroy_ask(GtkWidget *widget, GdkEvent *event,
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


