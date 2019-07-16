/*********************************************************************
 * graphic_utl.h
 *
 * Copyright (C) 2001 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA> <Yoshihito KURIZUKA>
 *
 * Thanks to following contributors.
 *  <Kenzen TAKEUCHI> <Tomoyuki TSUBATA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: graphic_utl.h,v 1.13 2004/02/08 08:21:48 sasaoka Exp $ */

#ifndef GRAPHIC_UTL_H
#define GRAPHIC_UTL_H

#include "rc.h"
#include "bin_io.h"
#include "math_utl.h"
#include "fem_struct.h"

typedef enum{
	GRAPHIC_RGB24,
	GRAPHIC_RGBA64,
	GRAPHIC_MONO8,
	GRAPHIC_MONO16,
	GRAPHIC_SCALAR,
	GRAPHIC_VECT2,
	GRAPHIC_VECT3,
	GRAPHIC_VECT6,
	GRAPHIC_NOTYPE
} GRAPHIC_TYPE;

typedef enum{
	GRAPHIC_BMP,
	GRAPHIC_TIFF,
	GRAPHIC_DICOM, 
	GRAPHIC_AVS_FIELD,
	GRAPHIC_NOSRC
} GRAPHIC_SOURCE;

typedef struct{
	unsigned long x;
	unsigned long y;
	unsigned long z;
} GRAPHIC_SIZE;

typedef struct{
	unsigned char r;  /* 8bit 256-depth each color */
	unsigned char g;
	unsigned char b;
} GRAPHIC_RGB;

typedef struct{
	unsigned short r;   /* 16bit 65536-depth each color */
	unsigned short g;
	unsigned short b;
	unsigned short a;  /* Alpha */
} GRAPHIC_RGBA;

typedef union{
	GRAPHIC_RGB ***rgb24;              
	GRAPHIC_RGBA ***rgba64;              
	unsigned char ***mono8;
	long ***mono16;
	double ***scalar;
	TRANSLATION2D ***vect2;
	TRANSLATION3D ***vect3;
	TRANS_ROTATION3D ***vect6;
} GRAPHIC_DATA;

typedef enum{
	GRAPHIC_UNIFORM,
	GRAPHIC_RECTILINEAR,
	GRAPHIC_IRREGULAR,
	GRAPHIC_NOCOORD
} GRAPHIC_COORD_TYPE;

#define GRAPHIC_INFO (10) 
typedef struct {
	GRAPHIC_SIZE size;                 /* $B2hA|$N%5%$%:(B [pixel] */
	TRANSLATION3D resolution;          /* $B2hA|2rA|EY(B [pixel/mm] */
	GRAPHIC_TYPE type;                 /* $B2hA|%?%$%W(B */
	GRAPHIC_SOURCE source;             /* $B2hA|%U%)!<%^%C%H(B */
	GRAPHIC_DATA data;                 /* $B%G!<%?(B */
	GRAPHIC_COORD_TYPE coord_type;     /* $B:BI8%?%$%W(B */
	char **c_info;                     /* the other informations(char) */
	int i_info[GRAPHIC_INFO];          /* the other informations(int) */
	double d_info[GRAPHIC_INFO];       /* the other informations(double) */
} GRAPHIC;

typedef struct {
	GRAPHIC_RGB *bgr24;
	GRAPHIC_RGB *br24;
	GRAPHIC_RGB *hot24;
	GRAPHIC_RGB *cold24;
	GRAPHIC_RGB *grayscale24;
} GRAPHIC_CONTOUR;

typedef enum{
	GRAPHIC_CONT_BGR24,       /* Blue->Cyan->Green->Yellow->Red (size: 1021) */
	GRAPHIC_CONT_BR24,        /* Blue -> Red (size: 511)*/
	GRAPHIC_CONT_HOT24,       /* Red -> Yellow (size: 256) */
	GRAPHIC_CONT_COLD24,      /* White -> Blue (size: 256) */
	GRAPHIC_CONT_GRAYSCALE24, /* White -> Black (size 256) */
} GRAPHIC_CONTOUR_TYPE;


/* graphic_utl.c */
RC allocate_graphic_data(GRAPHIC *gr); 

RC allocate_graphic_rgb24_2D(GRAPHIC_SIZE size, GRAPHIC_RGB ***rgb24);
RC allocate_graphic_rgba64_2D(GRAPHIC_SIZE size, GRAPHIC_RGBA ***rgba64);
RC allocate_graphic_mono8_2D(GRAPHIC_SIZE size, unsigned char ***mono8);
RC allocate_graphic_mono16_2D(GRAPHIC_SIZE size, long ***mono16);
RC allocate_graphic_scalar_2D(GRAPHIC_SIZE size, double ***scalar);
RC allocate_graphic_vect2_2D(GRAPHIC_SIZE size, TRANSLATION2D ***vect2);
RC allocate_graphic_vect3_2D(GRAPHIC_SIZE size, TRANSLATION3D ***vect3);
RC allocate_graphic_vect7_2D(GRAPHIC_SIZE size, TRANS_ROTATION3D ***vect6);

RC allocate_graphic_rgb24_3D(GRAPHIC_SIZE size, GRAPHIC_RGB ****rgb24);
RC allocate_graphic_rgba64_3D(GRAPHIC_SIZE size, GRAPHIC_RGBA ****rgba64);
RC allocate_graphic_mono8_3D(GRAPHIC_SIZE size, unsigned char ****mono8);
RC allocate_graphic_mono16_3D(GRAPHIC_SIZE size, long ****mono16);
RC allocate_graphic_scalar_3D(GRAPHIC_SIZE size, double ****scalar);
RC allocate_graphic_vect2_3D(GRAPHIC_SIZE size, TRANSLATION2D ****vect2);
RC allocate_graphic_vect3_3D(GRAPHIC_SIZE size, TRANSLATION3D ****vect3);
RC allocate_graphic_vect6_3D(GRAPHIC_SIZE size, TRANS_ROTATION3D ****vect6);

RC free_graphic_data(GRAPHIC gr); 

RC free_graphic_rgb24_2D(GRAPHIC_SIZE size, GRAPHIC_RGB **rgb24);
RC free_graphic_rgba64_2D(GRAPHIC_SIZE size, GRAPHIC_RGBA **rgba64);
RC free_graphic_mono8_2D(GRAPHIC_SIZE size, unsigned char **mono8);
RC free_graphic_mono16_2D(GRAPHIC_SIZE size, long **mono16);
RC free_graphic_scalar_2D(GRAPHIC_SIZE size, double **scalar);
RC free_graphic_vect2_2D(GRAPHIC_SIZE size, TRANSLATION2D **vect2);
RC free_graphic_vect3_2D(GRAPHIC_SIZE size, TRANSLATION3D **vect3);
RC free_graphic_vect6_2D(GRAPHIC_SIZE size, TRANS_ROTATION3D **vect6);

RC free_graphic_rgb24_3D(GRAPHIC_SIZE size, GRAPHIC_RGB ***rgb24);
RC free_graphic_rgba64_3D(GRAPHIC_SIZE size, GRAPHIC_RGBA ***rgba64);
RC free_graphic_mono8_3D(GRAPHIC_SIZE size, unsigned char ***mono8);
RC free_graphic_mono16_3D(GRAPHIC_SIZE size, long ***mono16);
RC free_graphic_scalar_3D(GRAPHIC_SIZE size, double ***scalar);
RC free_graphic_vect2_3D(GRAPHIC_SIZE size, TRANSLATION2D ***vect2);
RC free_graphic_vect3_3D(GRAPHIC_SIZE size, TRANSLATION3D ***vect3);
RC free_graphic_vect6_3D(GRAPHIC_SIZE size, TRANS_ROTATION3D ***vect6);


RC copy_graphic_info(GRAPHIC src, GRAPHIC *dest);
void init_graphic(GRAPHIC *gr);
GRAPHIC_SIZE init_graphic_size(void);
GRAPHIC_RGB init_graphic_rgb(void);
GRAPHIC_RGBA init_graphic_rgba(void);
void print_graphic_info(FILE *fp, GRAPHIC gr);

RC figure_rotation_2D(GRAPHIC src, GRAPHIC *dest, double angle);
RC matrix2grapic_normal(double **mat, int recood, int column, GRAPHIC *gr);
RC matrix2grapic_normal_area(double **mat, int row, int row1, int row2,
                             int col, int col1, int col2, GRAPHIC *gr);

RC make_graphic_contour(GRAPHIC_CONTOUR *cnt);
void free_graphic_contour(GRAPHIC_CONTOUR *cnt);
RC get_graphic_contour_rgb24(const double value, const double min,
                             const double max, GRAPHIC_CONTOUR cnt,
                             GRAPHIC_CONTOUR_TYPE type, GRAPHIC_RGB *color);
RC get_graphic_color_num_rgb24(const GRAPHIC_RGB color, 
                               GRAPHIC_CONTOUR_TYPE type, int *num);
RC chk_graphic_range_rgb24(GRAPHIC src, GRAPHIC_RGB *min, GRAPHIC_RGB*max,
                           GRAPHIC_RGB *avg);
RC chk_graphic_range_mono16(GRAPHIC src, long *min, long *max, long *avg);
RC chk_graphic_range_mono8(GRAPHIC src, unsigned char *min, unsigned char *max,
                           unsigned char *avg);
RC chk_graphic_range_scalar(GRAPHIC src, double *min, double *max, double *avg);
RC conv_mono16_to_rgb24_1(GRAPHIC src, GRAPHIC *dest);
RC conv_mono16_to_rgb24_2(GRAPHIC src, GRAPHIC *dest, double ratio);
RC extract_graphic_outline(GRAPHIC src, GRAPHIC *dest, const double *operator);
RC graphic_exist_chk(GRAPHIC gr);
RC inverse_graphic(GRAPHIC src, GRAPHIC *dest);
RC make_graphic_gradient_color_map(GRAPHIC src, GRAPHIC *dest,
                                   GRAPHIC_CONTOUR cnt,
                                   GRAPHIC_CONTOUR_TYPE type);
RC make_graphic_gradient_scalar(GRAPHIC src, GRAPHIC *dest);
RC conv_graphic_scalar2color_map(GRAPHIC src, GRAPHIC *dest,
                                 GRAPHIC_CONTOUR cnt,
                                 GRAPHIC_CONTOUR_TYPE type);
RC get_graphic_value_fem_node(const GRAPHIC src, const FEM_NODE node,
                              void *value);


/* avs_field_utl.c            */
/* KGT$B$N(BAVS_FIELD$B$N$_%5%]!<%H(B */
RC read_avs_field(FILE *fp, GRAPHIC *gr);
RC read_avs_field_heder(FILE *fp, GRAPHIC *gr);
RC read_avs_field_data(FILE *fp, GRAPHIC *gr);
RC read_avs_field_data_value(FILE *fp, GRAPHIC *gr);
RC read_avs_field_data_coord(FILE *fp, GRAPHIC *gr);
RC write_avs_field_header(FILE *fp, GRAPHIC gr);
RC write_avs_field_data(FILE *fp, GRAPHIC gr);

/*
 * bmp_utl.c
 * $B05=L$7$?(BBMP$B%$%a!<%8$O%5%]!<%H$7$J$$(B
 * Write $B$O(B24bit$B$N$_BP1~(B
 * Read $B$O(B24bit,8bit[palette]$B$N$_BP1~(B
 * (4bit,2bit$B$N2hA|$O%Q%C%G%#%s%0$N=hM}$,ITL@$N$?$aL$BP1~(B)
 */
RC write_bmp(FILE *fp, GRAPHIC gr);
RC write_bmp_header(FILE *fp, GRAPHIC gr);
RC write_bmp_image_rgb24(FILE *fp, GRAPHIC gr);
RC read_bmp(FILE *fp, GRAPHIC *gr);
RC read_bmp_header(FILE *fp, GRAPHIC *gr);

/*
 * tiff_utl.c
 * 1$B$D$N%U%!%$%kFb$KJ#?t$N2hA|$r;}$D(BTIFF$B%$%a!<%8(B
 * $B5Z$S!"05=L$7$?(BTIFF$B%$%a!<%8$O%5%]!<%H$7$J$$(B
 * $B2rA|EY$NC10L$K4X$9$k>pJs$,ITL@$N$?$a!"5-=R$7$F$$$kCM$r$=$N$^$^F~NO$7$F$$$k(B
 * $B8=:_BP1~$7$F$$$k$N$O(B "24bit","8bit$B%Q%l%C%H%+%i!<(B","8Bit$B%b%N%/%m(B","2bit"
 */ 
RC read_tiff(FILE *fp, GRAPHIC *gr);
RC read_tiff_header(FILE *fp, GRAPHIC *gr, ENDIAN *endian);
RC read_tiff_image(FILE *fp, GRAPHIC *gr, ENDIAN endian);

/* Macro */
#define GRAPHIC_SIZE_LOOP_2D(Y, X)\
{ int ii1, ii2;\
   for(ii1=0; (unsigned long)ii1<Y; ii1++){\
    for(ii2=0; (unsigned long)ii2<X; ii2++){

#define GRAPHIC_SIZE_LOOP_2D_END  } } }


#define GRAPHIC_SIZE_LOOP_3D(Z, Y, X)\
{ int ii1, ii2, ii3;\
   for(ii1=0; (unsigned long)ii1<Z; ii1++){\
    for(ii2=0; (unsigned long)ii2<Y; ii2++){\
     for(ii3=0; (unsigned long)ii3<X; ii3++){

#define GRAPHIC_SIZE_LOOP_3D_END  } } } }


#define GRAPHIC_SIZE_CHK_2D(SIZE)\
{\
    if( (SIZE.x < 1) || (SIZE.y < 1) ){\
        fprintf(stderr, "Bad graphic size X or Y.");\
        return(ARG_ERROR_RC);\
    }\
}

#define GRAPHIC_SIZE_CHK_3D(SIZE)\
{\
    if( (SIZE.x < 1) || (SIZE.y < 1) || (SIZE.z < 1) ){\
        fprintf(stderr, "Bad graphic size X or Y or Z.");\
        return(ARG_ERROR_RC);\
    }\
}


#endif /* GRAPHIC_UTL_H */



