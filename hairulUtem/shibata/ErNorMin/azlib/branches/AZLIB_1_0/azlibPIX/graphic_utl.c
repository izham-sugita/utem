/*********************************************************************
 * graphic_utl.c
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

/* $Id: graphic_utl.c,v 1.15 2004/02/08 08:21:33 sasaoka Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "rc.h"
#include "bin_io.h"
#include "fem_struct.h"
#include "graphic_utl.h"

static RC make_contour_bgr24(GRAPHIC_RGB **cnt);
static RC make_contour_br24(GRAPHIC_RGB **cnt);
static RC make_contour_hot24(GRAPHIC_RGB **cnt);
static RC make_contour_cold24(GRAPHIC_RGB **cnt);
static RC make_contour_grayscale24(GRAPHIC_RGB **cnt);

static RC conv_graphic_binary_scalar(GRAPHIC src, GRAPHIC *dest);

static int operator_mat[27][3] = { {-1, -1, -1}, {0, -1, -1}, {1, -1, -1},
                                   {-1,  0, -1}, {0,  0, -1}, {1,  0, -1},
                                   {-1,  1, -1}, {0,  1, -1}, {1,  1, -1},
                                   {-1, -1,  0}, {0, -1,  0}, {1, -1,  0},
                                   {-1,  0,  0}, {0,  0,  0}, {1,  0,  0},
                                   {-1,  1,  0}, {0,  1,  0}, {1,  1,  0},
                                   {-1, -1,  1}, {0, -1,  1}, {1, -1,  1},
                                   {-1,  0,  1}, {0,  0,  1}, {1,  0,  1},
                                   {-1,  1,  1}, {0,  1,  1}, {1,  1,  1} }; 

static double laplacian_op[9] = { 0.0, -1.0,  0.0,
                                 -1.0,  4.0, -1.0,
                                  0.0, -1.0,  0.0};

static double laplacian_op_3d[27] = { 0.0,  0.0,  0.0,
                                      0.0, -1.0,  0.0,
                                      0.0,  0.0,  0.0,
                                      0.0, -1.0,  0.0,
                                     -1.0,  6.0, -1.0,
                                      0.0, -1.0,  0.0,
                                      0.0,  0.0,  0.0,
                                      0.0, -1.0,  0.0,
                                      0.0,  0.0,  0.0 };

/* $B%5%$%:$J$I$N>pJs$O@h$K%;%C%H$5$l$F$$$k$b$N$H$9$k(B */
RC allocate_graphic_data(GRAPHIC *gr)
{
	switch(gr->type){
	case GRAPHIC_RGB24:
		RC_TRY( allocate_graphic_rgb24_3D(gr->size, &(gr->data.rgb24)) );
		break;
	case GRAPHIC_RGBA64:
		RC_TRY( allocate_graphic_rgba64_3D(gr->size, &(gr->data.rgba64)) );
		break;
	case GRAPHIC_MONO8:
		RC_TRY( allocate_graphic_mono8_3D(gr->size, &(gr->data.mono8)) );
		break;
	case GRAPHIC_MONO16:
		RC_TRY( allocate_graphic_mono16_3D(gr->size, &(gr->data.mono16)) );
		break;
	case GRAPHIC_SCALAR:
		RC_TRY( allocate_graphic_scalar_3D(gr->size, &(gr->data.scalar)) );
		break;
	case GRAPHIC_VECT2:
		RC_TRY( allocate_graphic_vect2_3D(gr->size, &(gr->data.vect2)) );
		break;
	case GRAPHIC_VECT3:
		RC_TRY( allocate_graphic_vect3_3D(gr->size, &(gr->data.vect3)) );
		break;
	case GRAPHIC_VECT6:
		RC_TRY( allocate_graphic_vect6_3D(gr->size, &(gr->data.vect6)) );
		break;
	default:
		fprintf(stderr,"unknown graphic type.\n");\
		return(UNKNOWN_ERROR_RC);
	}

	return(NORMAL_RC);
}


#define FUNC_ALLOC_GRAPHIC_DATA_2D(func, data_type, type, init_object)  \
RC func(GRAPHIC_SIZE size, data_type ***type)\
{\
	int ii1,ii2;\
	GRAPHIC_SIZE_CHK_2D(size);\
	(*type) = (data_type **)malloc(size.y * sizeof(data_type *));\
	RC_NULL_CHK(*type);\
	for(ii1=0; (unsigned long)ii1<size.y; ii1++){\
		(*type)[ii1] = (data_type *)malloc(size.x * sizeof(data_type));\
		RC_NULL_CHK((*type)[ii1]); \
		for(ii2=0; (unsigned long)ii2<size.x; ii2++){\
			(*type)[ii1][ii2] = init_object;\
		}\
	}\
	return(NORMAL_RC);\
}

FUNC_ALLOC_GRAPHIC_DATA_2D(allocate_graphic_rgb24_2D, GRAPHIC_RGB, rgb24,
                           init_graphic_rgb())
FUNC_ALLOC_GRAPHIC_DATA_2D(allocate_graphic_rgba64_2D, GRAPHIC_RGBA, rgb64,
                           init_graphic_rgba())
FUNC_ALLOC_GRAPHIC_DATA_2D(allocate_graphic_mono8_2D, unsigned char , mono8, 0)
FUNC_ALLOC_GRAPHIC_DATA_2D(allocate_graphic_mono16_2D, long, mono16, 0)
FUNC_ALLOC_GRAPHIC_DATA_2D(allocate_graphic_scalar_2D, double, scalar, 0.0)
FUNC_ALLOC_GRAPHIC_DATA_2D(allocate_graphic_vect2_2D,TRANSLATION2D,vect2,
                           init_translation2d())
FUNC_ALLOC_GRAPHIC_DATA_2D(allocate_graphic_vect3_2D,TRANSLATION3D,vect3,
                           init_translation3d())
FUNC_ALLOC_GRAPHIC_DATA_2D(allocate_graphic_vect6_2D, TRANS_ROTATION3D, vect6,
                           init_trans_rotation3d())


#define FUNC_ALLOC_GRAPHIC_DATA_3D(func, func_2d, data_type, type)  \
RC func(GRAPHIC_SIZE size, data_type ****type)\
{\
	int ii1;\
	GRAPHIC_SIZE_CHK_3D(size);\
	(*type) = (data_type ***)malloc(size.z * sizeof(data_type **));\
	RC_NULL_CHK(*type);\
	for(ii1=0; (unsigned long)ii1<size.z; ii1++){\
		RC_TRY( func_2d(size, &((*type)[ii1])) ); \
	}\
	return(NORMAL_RC);\
}

FUNC_ALLOC_GRAPHIC_DATA_3D(allocate_graphic_rgb24_3D,
                           allocate_graphic_rgb24_2D, GRAPHIC_RGB, rgb24)
FUNC_ALLOC_GRAPHIC_DATA_3D(allocate_graphic_rgba64_3D,
                           allocate_graphic_rgba64_2D, GRAPHIC_RGBA,rgba64)
FUNC_ALLOC_GRAPHIC_DATA_3D(allocate_graphic_mono8_3D,
                           allocate_graphic_mono8_2D, unsigned char, mono8)
FUNC_ALLOC_GRAPHIC_DATA_3D(allocate_graphic_mono16_3D,
                           allocate_graphic_mono16_2D, long, mono16)
FUNC_ALLOC_GRAPHIC_DATA_3D(allocate_graphic_scalar_3D,
                           allocate_graphic_scalar_2D, double, scalar)
FUNC_ALLOC_GRAPHIC_DATA_3D(allocate_graphic_vect2_3D,
                           allocate_graphic_vect2_2D, TRANSLATION2D, vect2)
FUNC_ALLOC_GRAPHIC_DATA_3D(allocate_graphic_vect3_3D,
                           allocate_graphic_vect3_2D, TRANSLATION3D, vect3)
FUNC_ALLOC_GRAPHIC_DATA_3D(allocate_graphic_vect6_3D,
                           allocate_graphic_vect6_2D, TRANS_ROTATION3D,vect6)


RC free_graphic_data(GRAPHIC gr)
{
	switch(gr.type){
	case GRAPHIC_RGB24:
		RC_TRY( free_graphic_rgb24_3D(gr.size, gr.data.rgb24) );
		break;
	case GRAPHIC_RGBA64:
		RC_TRY( free_graphic_rgba64_3D(gr.size, gr.data.rgba64) );
		break;
	case GRAPHIC_MONO8:
		RC_TRY( free_graphic_mono8_3D(gr.size, gr.data.mono8) );
		break;
	case GRAPHIC_MONO16:
		RC_TRY( free_graphic_mono16_3D(gr.size, gr.data.mono16) );
		break;
	case GRAPHIC_SCALAR:
		RC_TRY( free_graphic_scalar_3D(gr.size, gr.data.scalar) );
		break;
	case GRAPHIC_VECT2:
		RC_TRY( free_graphic_vect2_3D(gr.size, gr.data.vect2) );
		break;
	case GRAPHIC_VECT3:
		RC_TRY( free_graphic_vect3_3D(gr.size, gr.data.vect3) );
		break;
	case GRAPHIC_VECT6:
		RC_TRY( free_graphic_vect6_3D(gr.size, gr.data.vect6) );
		break;
	default:
		fprintf(stderr, "unknown graphic type.\n");\
		return(UNKNOWN_ERROR_RC);
	}

	return(NORMAL_RC);

}

#define FUNC_FREE_GRAPHIC_DATA_2D(func, data_type, type)  \
RC func(GRAPHIC_SIZE size, data_type **type)\
{\
	if(size.y > 0){\
		int ii1;\
		for(ii1=0; (unsigned long)ii1<size.y; ii1++){\
			free((void *)type[ii1]);\
			type[ii1] = (data_type *)NULL;\
		}\
	}\
	free((void *)type);\
	return(NORMAL_RC);\
}

FUNC_FREE_GRAPHIC_DATA_2D(free_graphic_rgb24_2D,  GRAPHIC_RGB, rgb24)
FUNC_FREE_GRAPHIC_DATA_2D(free_graphic_rgba64_2D, GRAPHIC_RGBA, rgb64)
FUNC_FREE_GRAPHIC_DATA_2D(free_graphic_mono8_2D,  unsigned char, mono8)
FUNC_FREE_GRAPHIC_DATA_2D(free_graphic_mono16_2D, long, mono16)
FUNC_FREE_GRAPHIC_DATA_2D(free_graphic_scalar_2D, double, scalar)
FUNC_FREE_GRAPHIC_DATA_2D(free_graphic_vect2_2D,  TRANSLATION2D, vect2)
FUNC_FREE_GRAPHIC_DATA_2D(free_graphic_vect3_2D,  TRANSLATION3D, vect3)
FUNC_FREE_GRAPHIC_DATA_2D(free_graphic_vect6_2D,  TRANS_ROTATION3D, vect6)


#define FUNC_FREE_GRAPHIC_DATA_3D(func, func_2d, data_type, type)  \
RC func(GRAPHIC_SIZE size, data_type ***type)\
{\
	if(size.z > 0){\
		int ii1;\
		for(ii1=0; (unsigned long)ii1<size.z; ii1++){\
			RC_TRY( func_2d(size, type[ii1]) );\
			type[ii1] = (data_type **)NULL;\
		}\
	}\
	free((void *)type);\
	return(NORMAL_RC);\
}

FUNC_FREE_GRAPHIC_DATA_3D(free_graphic_rgb24_3D,
                          free_graphic_rgb24_2D,  GRAPHIC_RGB, rgb24)
FUNC_FREE_GRAPHIC_DATA_3D(free_graphic_rgba64_3D,
                          free_graphic_rgba64_2D, GRAPHIC_RGBA, rgba64)
FUNC_FREE_GRAPHIC_DATA_3D(free_graphic_mono8_3D,
                          free_graphic_mono8_2D,  unsigned char, mono8)
FUNC_FREE_GRAPHIC_DATA_3D(free_graphic_mono16_3D,
                          free_graphic_mono16_2D, long, mono16)
FUNC_FREE_GRAPHIC_DATA_3D(free_graphic_scalar_3D,
                          free_graphic_scalar_2D, double, scalar)
FUNC_FREE_GRAPHIC_DATA_3D(free_graphic_vect2_3D,
                          free_graphic_vect2_2D,  TRANSLATION2D, vect2)
FUNC_FREE_GRAPHIC_DATA_3D(free_graphic_vect3_3D,
                          free_graphic_vect3_2D,  TRANSLATION3D, vect3)
FUNC_FREE_GRAPHIC_DATA_3D(free_graphic_vect6_3D,
                          free_graphic_vect6_2D,  TRANS_ROTATION3D, vect6)


RC copy_graphic_info(GRAPHIC src, GRAPHIC *dest)
{
	int ii1;
	RC_NULL_CHK(dest);
	dest->size = src.size;
	dest->resolution = src.resolution;
	dest->type = src.type;
	dest->source = src.source;
	dest->coord_type = src.coord_type;
	for(ii1=0; ii1<GRAPHIC_INFO; ii1++){
		dest->i_info[ii1] = src.i_info[ii1];
		dest->d_info[ii1] = src.d_info[ii1];
	}

	return(NORMAL_RC);
}


void init_graphic(GRAPHIC *gr)
{
	int ii1;
	gr->size = init_graphic_size();
	gr->resolution = init_translation3d();
	gr->type = GRAPHIC_NOTYPE;
	gr->source = GRAPHIC_NOSRC;
	gr->data.rgb24 = NULL;
	gr->coord_type = GRAPHIC_NOCOORD;
	gr->c_info = NULL;
	for(ii1=0; ii1<GRAPHIC_INFO; ii1++){
		gr->i_info[ii1] = 0;
		gr->d_info[ii1] = 0.0;
	}
}


GRAPHIC_SIZE init_graphic_size(void)
{
	GRAPHIC_SIZE ret = {0, 0, 0};
	return(ret);
}


GRAPHIC_RGB init_graphic_rgb(void)
{
	GRAPHIC_RGB ret = {0, 0, 0};
	return(ret);
}

GRAPHIC_RGBA init_graphic_rgba(void)
{
	GRAPHIC_RGBA ret = {0, 0, 0, 1};
	return(ret);
}

void print_graphic_info(FILE *fp, GRAPHIC gr)
{
	fprintf(stderr, "graphic size [x] = %ld\n", gr.size.x);
	fprintf(stderr, "graphic size [y] = %ld\n", gr.size.y);
	fprintf(stderr, "graphic size [z] = %ld\n", gr.size.z);

	fprintf(stderr, "graphic resolution[x] = %e\n", gr.resolution.x);
	fprintf(stderr, "graphic resolution[y] = %e\n", gr.resolution.y);
	fprintf(stderr, "graphic resolution[z] = %e\n", gr.resolution.z);

	if(gr.coord_type == GRAPHIC_UNIFORM){
		fprintf(stderr, "COORD TYPE = uniform\n");
	}else if(gr.coord_type == GRAPHIC_RECTILINEAR){
		fprintf(stderr, "COORD TYPE = rectilinear\n");
	}else if(gr.coord_type == GRAPHIC_IRREGULAR){
		fprintf(stderr, "COORD TYPE = irregular\n");
	}else{
		fprintf(stderr, "UNKNOWN COORD TYPE\n");
	}

	if(gr.source == GRAPHIC_BMP){
		fprintf(stderr,"GRAPHIC SOURCE = BMP\n");
	}else if(gr.source == GRAPHIC_TIFF){
		fprintf(stderr,"GRAPHIC SOURCE = TIFF\n");
	}else if(gr.source == GRAPHIC_DICOM){
		fprintf(stderr,"GRAPHIC SOURCE = DICOM\n");
	}else if(gr.source == GRAPHIC_AVS_FIELD){
		fprintf(stderr,"GRAPHIC SOURCE = AVS_FIELD\n");
	}else{
		fprintf(stderr,"UNKNOWN GRAPHIC SOURCE \n");
	}

	if(gr.type == GRAPHIC_RGB24){
		fprintf(stderr,"GRAPHIC TYPE = RGB24\n");
	}else if(gr.type == GRAPHIC_RGBA64){
		fprintf(stderr,"GRAPHIC TYPE = RGBA64\n");
	}else if(gr.type == GRAPHIC_MONO8){
		fprintf(stderr,"GRAPHIC TYPE = MONO8\n");
	}else if(gr.type == GRAPHIC_MONO16){
		fprintf(stderr,"GRAPHIC TYPE = MONO16\n");
	}else if(gr.type == GRAPHIC_SCALAR){
		fprintf(stderr,"GRAPHIC TYPE = SCALAR\n");
	}else if(gr.type == GRAPHIC_VECT2){
		fprintf(stderr,"GRAPHIC TYPE = VECT2\n");
	}else if(gr.type == GRAPHIC_VECT3){
		fprintf(stderr,"GRAPHIC TYPE = VECT3\n");
	}else{
		fprintf(stderr, "UNKNOWN GRAPHIC TYPE\n");
	}
}


RC figure_rotation_2D(GRAPHIC src, GRAPHIC *dest, double angle)
{
	if( (dest->type == GRAPHIC_RGB24)  || (dest->type == GRAPHIC_RGBA64) ||
	    (dest->type == GRAPHIC_MONO8)  || (dest->type == GRAPHIC_MONO16) ||
	    (dest->type == GRAPHIC_SCALAR) ){
		if(dest == NULL) RC_TRY( allocate_graphic_data(dest) );
		GRAPHIC_SIZE_CHK_2D(src.size);
		RC_TRY( copy_graphic_info(src, dest) );
	}else{
		fprintf(stderr, "Sorry. This type graphics haven't implement yet.\n");
		return(IMPLEMENT_ERROR_RC);
	}

	GRAPHIC_SIZE_LOOP_3D(dest->size.z, dest->size.y, dest->size.x)
	{
		unsigned int  kk1=0, kk2=0, kk3 = 0;
		int itmp2, itmp3;
		double jj2, jj3;

		kk1 = ii1;
		itmp2 = ii2 - (dest->size.y / 2);
		itmp3 = ii3 - (dest->size.x / 2);
		jj2 = (itmp3*sin(angle)) + (itmp2*cos(angle));
		jj3 = (itmp3*cos(angle)) - (itmp2*sin(angle));
				 
		jj2 += (dest->size.y / 2);
		jj3 += (dest->size.x / 2);
		if( jj2 > ((int)jj2 + 0.5) ){
			kk2 = (int)jj2 +1;
		}else{
			kk2 = (int)jj2;
		}
		if( jj3 > ((int)jj3 + 0.5) ){
			kk3 = (int)jj3 + 1;
		}else{
			kk3 = (int)jj3;
		}

		if( (kk2 >= 0) && (kk2 < dest->size.y) ){ 
			if( (kk3 >=0) && (kk3 < dest->size.x) ){
				switch(dest->type){
				case GRAPHIC_MONO8:
					dest->data.mono8[ii1][ii2][ii3]
					    = src.data.mono8[kk1][kk2][kk3];
					break;
				case GRAPHIC_MONO16:
					dest->data.mono16[ii1][ii2][ii3]
					    = src.data.mono16[kk1][kk2][kk3];
					break;
				case GRAPHIC_RGB24:
					dest->data.rgb24[ii1][ii2][ii3].r
					    = src.data.rgb24[kk1][kk2][kk3].r;
					dest->data.rgb24[ii1][ii2][ii3].g
					    = src.data.rgb24[kk1][kk2][kk3].g;
					dest->data.rgb24[ii1][ii2][ii3].b
					    = src.data.rgb24[kk1][kk2][kk3].b;
					break;
				case GRAPHIC_RGBA64:
					dest->data.rgba64[ii1][ii2][ii3].r
					    = src.data.rgba64[kk1][kk2][kk3].r;
					dest->data.rgba64[ii1][ii2][ii3].g
					    = src.data.rgba64[kk1][kk2][kk3].g;
					dest->data.rgba64[ii1][ii2][ii3].b
					    = src.data.rgba64[kk1][kk2][kk3].b;
					dest->data.rgba64[ii1][ii2][ii3].a
					    = src.data.rgba64[kk1][kk2][kk3].a;
					break;
				case GRAPHIC_SCALAR:
					dest->data.scalar[kk1][kk2][kk3]
					    = src.data.scalar[ii1][ii2][ii3];
					break;
				default:
					return(UNKNOWN_ERROR_RC);
					break;
				}
			}
		}
	}
	GRAPHIC_SIZE_LOOP_3D_END;

	return(NORMAL_RC);
}


RC matrix2grapic_normal(double **mat, int row, int col, GRAPHIC *gr)
{
	RC_TRY( matrix2grapic_normal_area(mat, row, 0, row-1, col, 0, col-1, gr) );

	return(NORMAL_RC);
}


RC matrix2grapic_normal_area(double **mat, int row, int row1, int row2,
                             int col, int col1, int col2, GRAPHIC *gr)
{
	int ii1, ii2;
	double max, min;
	GRAPHIC_CONTOUR contour;
	
	if(row2 < row1) return(ARG_ERROR_RC);
	if(col2 < col1) return(ARG_ERROR_RC);

	init_graphic(gr);
	gr->size.x = col2 - col1 + 1;
	gr->size.y = row2 - row1 + 1;
	gr->size.z = 1;
	gr->type = GRAPHIC_RGB24;

	RC_TRY( allocate_graphic_data(gr) );

	max = min = mat[0][0];
	for(ii1=0;ii1<row;ii1++){
		for(ii2=0;ii2<col;ii2++){
			if(max < mat[ii1][ii2]) max = mat[ii1][ii2];
			if(min > mat[ii1][ii2]) min = mat[ii1][ii2];
		}
	}

	RC_TRY( make_graphic_contour(&contour) );

	GRAPHIC_SIZE_LOOP_2D(gr->size.y, gr->size.x)
	{
		RC_TRY( get_graphic_contour_rgb24(mat[row1+ii1][col1+ii2], min, max,
		                                    contour, GRAPHIC_CONT_BGR24,
		                                    &(gr->data.rgb24[0][ii1][ii2])) );
	}
	GRAPHIC_SIZE_LOOP_2D_END;

	free_graphic_contour(&contour);

	return(NORMAL_RC);
}


RC make_graphic_contour(GRAPHIC_CONTOUR *cnt)
{
	RC_TRY( make_contour_bgr24(&(cnt)->bgr24) );
	RC_TRY( make_contour_br24(&(cnt)->br24));
	RC_TRY( make_contour_hot24(&(cnt)->hot24) );
	RC_TRY( make_contour_cold24(&(cnt)->cold24) );
	RC_TRY( make_contour_grayscale24(&(cnt)->grayscale24) );
	return(NORMAL_RC);
}


void free_graphic_contour(GRAPHIC_CONTOUR *cnt)
{
	free(cnt->bgr24);
	free(cnt->br24);
	free(cnt->hot24);
	free(cnt->cold24);
	free(cnt->grayscale24);
	cnt->bgr24 = NULL;
	cnt->br24 = NULL;
	cnt->hot24 = NULL;
	cnt->cold24 = NULL;
	cnt->grayscale24 = NULL;
}


/* Blue Green Red $B$N=g$GCM$,9b$/$J$k(B 1021$B3,D4(B */
static RC make_contour_bgr24(GRAPHIC_RGB **cnt)
{
	int ii1;
	*cnt = (GRAPHIC_RGB *)malloc(sizeof(GRAPHIC_RGB) * 1021);
	RC_NULL_CHK( *cnt );

	for(ii1=0;ii1<1021;ii1++){
		if(ii1 < 256){
			(*cnt)[ii1].r = 0;
			(*cnt)[ii1].g = ii1;
			(*cnt)[ii1].b = 255;
		}else if(ii1 < 511){
			(*cnt)[ii1].r = 0;
			(*cnt)[ii1].g = 255;
			(*cnt)[ii1].b = 510 - ii1;
		}else if(ii1 < 766){
			(*cnt)[ii1].r = ii1 - 510;
			(*cnt)[ii1].g = 255;
			(*cnt)[ii1].b = 0;
		}else{
			(*cnt)[ii1].r = 255;
			(*cnt)[ii1].g = 1020 - ii1;
			(*cnt)[ii1].b = 0;
		}
	}

	return(NORMAL_RC);
}


/* Blue Red $B$N=g$GCM$,9b$/$J$k(B 511$B3,D4(B */
static RC make_contour_br24(GRAPHIC_RGB **cnt)
{
	int ii1;
	RC_NULL_CHK( *cnt = (GRAPHIC_RGB *)malloc(sizeof(GRAPHIC_RGB) * 511) );
	for(ii1=0;ii1<511;ii1++){
		if(ii1 < 256){
			(*cnt)[ii1].r = ii1;
			(*cnt)[ii1].g = 0;
			(*cnt)[ii1].b = 255;
		}else{
			(*cnt)[ii1].r = 255;
			(*cnt)[ii1].g = 0;
			(*cnt)[ii1].b = 510 - ii1;
		}
	}
	return(NORMAL_RC);
}


/* Red Yeloow $B$N=g$GCM$,9b$/$J$k(B 256$B3,D4(B */
static RC make_contour_hot24(GRAPHIC_RGB **cnt)
{
	int ii1;
	RC_NULL_CHK( *cnt = (GRAPHIC_RGB *)malloc(sizeof(GRAPHIC_RGB) * 256) );
	for(ii1=0;ii1<256;ii1++){
		(*cnt)[ii1].r = 255;
		(*cnt)[ii1].g = ii1;
		(*cnt)[ii1].b = 0;
	}
	return(NORMAL_RC);
}


/* White Blue $B$N=g$GCM$,9b$/$J$k(B 256$B3,D4(B */
static RC make_contour_cold24(GRAPHIC_RGB **cnt)
{
	int ii1;
	RC_NULL_CHK( *cnt = (GRAPHIC_RGB *)malloc(sizeof(GRAPHIC_RGB) * 256) );
	for(ii1=0;ii1<256;ii1++){
		(*cnt)[ii1].r = 255 - ii1;
		(*cnt)[ii1].g = 255 - ii1;
		(*cnt)[ii1].b = 255;
	}
	return(NORMAL_RC);
}


/* Black White $B$N=g$GCM$,9b$/$J$k(B 256$B3,D4(B */
static RC make_contour_grayscale24(GRAPHIC_RGB **cnt)
{
	int ii1;
	RC_NULL_CHK( *cnt = (GRAPHIC_RGB *)malloc(sizeof(GRAPHIC_RGB) * 256) );
	for(ii1=0;ii1<256;ii1++){
		(*cnt)[ii1].r = ii1;
		(*cnt)[ii1].g = ii1;
		(*cnt)[ii1].b = ii1;
	}
	return(NORMAL_RC);
}


#define GET_GRAPHIC_CONTOUR_COLOR(contour, size) \
{\
    RC_NULL_CHK(contour);\
    color_num = (int)( ((size-1)*(value - min)) / (max - min) );\
    if((color_num < 0) || (color_num > (size-1))) return(CAL_ERROR_RC);\
    else *color = contour[color_num];\
}
/* mix < value < max$B$G$"$k$H$-$N!$(Bvalue$B$N%+%i!<%3%s%??'(B(color)$B$rJV$9!%(B*/
RC get_graphic_contour_rgb24(const double value, const double min,
                             const double max, GRAPHIC_CONTOUR cnt,
                             GRAPHIC_CONTOUR_TYPE type, GRAPHIC_RGB *color)
{
	int color_num;
	if(max < min)   return(ARG_ERROR_RC);
	if(value < min) return(ARG_ERROR_RC);
	if(value > max) return(ARG_ERROR_RC);
	if(nearly_eq(max, min)){
		color->r = color->g = color->b = 0;
		return(NORMAL_RC);
	}
	switch(type){
	case GRAPHIC_CONT_BGR24:
		GET_GRAPHIC_CONTOUR_COLOR(cnt.bgr24, 1021)
		break;
	case GRAPHIC_CONT_BR24:
		GET_GRAPHIC_CONTOUR_COLOR(cnt.br24, 511)
		break;
	case GRAPHIC_CONT_HOT24:
		GET_GRAPHIC_CONTOUR_COLOR(cnt.hot24, 256)
		break;
	case GRAPHIC_CONT_COLD24:
		GET_GRAPHIC_CONTOUR_COLOR(cnt.cold24, 256)
		break;
	case GRAPHIC_CONT_GRAYSCALE24:
		GET_GRAPHIC_CONTOUR_COLOR(cnt.grayscale24, 256)
		break;
	default:
		fprintf(stderr,"unknown color-contour type.\n");
		return(ARG_ERROR_RC);
		break;
	}
	return(NORMAL_RC);
}


/* color$B$G$"$?$($i$l$??'$N%+%i!<%3%s%?>e$NHV9f(B($BBg$-$5(B)$B$rJV$9!%(B*/
RC get_graphic_color_num_rgb24(const GRAPHIC_RGB color,
                               GRAPHIC_CONTOUR_TYPE type, int *num)
{
	switch(type){
	case GRAPHIC_CONT_BGR24:
		if(color.b == 255){
			*num = color.g;
		}else if(color.g == 255){
			if(color.r == 0){
				*num = 510 - color.b;
			}else{
				*num = 510 + color.r;
			}
		}else{
			*num = 1020 -  color.g;
		}
		break;
	case GRAPHIC_CONT_BR24:
		if(color.b == 255){
			*num = color.r;
		}else{
			*num = 510 - color.b;
		}
		break;
	case GRAPHIC_CONT_HOT24:
		*num = color.g;
		break;
	case GRAPHIC_CONT_COLD24:
		*num = 255 - color.r;
		break;
	case GRAPHIC_CONT_GRAYSCALE24:
		*num = color.r;
		break;
	default:
		fprintf(stderr,"unknown color-contour type.\n");
		return(ARG_ERROR_RC);
		break;
	}
	return(NORMAL_RC);
}


/*
 * $B51EY$r85$K(BRGB$B$N:GBgCM$r9M$($k!%J?6Q$O3F?'$NJ?6Q$G$"$k!%(B
 * CIE($B9q:]>HL@0Q0w2q(B)XYZ$BI=?'7O(B Y$B$,51EY$K$"$?$k(B
 * X = 0.49000R + 0.31000G + 0.20000B
 * Y = 0.17697R + 0.81240G + 0.01063B
 * Z = 0.00000R + 0.01000G + 0.99000B
 */
RC chk_graphic_range_rgb24(GRAPHIC src, GRAPHIC_RGB *min, GRAPHIC_RGB *max,
                           GRAPHIC_RGB *avg)
{
	double Y_min, Y_max, t_avg_r, t_avg_g, t_avg_b;
	double size = (src.size.z)*(src.size.y)*(src.size.x);
	if(size == 0) return(ARG_ERROR_RC);
	RC_NULL_CHK(src.data.rgb24);
	Y_min = 255;
	Y_max = 0;
	t_avg_r = t_avg_g = t_avg_b = 0.0;
	GRAPHIC_SIZE_LOOP_3D(src.size.z, src.size.y, src.size.x)
	{
		GRAPHIC_RGB rgb = src.data.rgb24[ii1][ii2][ii3];
		double Y = 0.17697*rgb.r + 0.81240*rgb.g + 0.01063*rgb.b;
		if(Y < Y_min){
			Y_min = Y;
			*min = rgb;
		}
		if(Y > Y_max){
			Y_max = Y;
			*max = rgb;
		}
		t_avg_r += rgb.r/size;
		t_avg_g += rgb.g/size;
		t_avg_b += rgb.b/size;
	}GRAPHIC_SIZE_LOOP_3D_END;
	avg->r = (unsigned char)t_avg_r;
	avg->g = (unsigned char)t_avg_g;
	avg->b = (unsigned char)t_avg_b;

	return(NORMAL_RC);
}


RC chk_graphic_range_mono16(GRAPHIC src, long *min, long *max, long *avg)
{
	double t_avg = 0.0;
	double size = (src.size.z)*(src.size.y)*(src.size.x);
	if(size == 0) return(ARG_ERROR_RC);
	RC_NULL_CHK(src.data.mono16);
	*min = (1<<16);
	*max = -((1<<15)-1);
	GRAPHIC_SIZE_LOOP_3D(src.size.z, src.size.y, src.size.x)
	{
		long cl = src.data.mono16[ii1][ii2][ii3];
		if(cl < *min) *min = cl;
		if(cl > *max) *max = cl;
		t_avg += cl/size;
	}GRAPHIC_SIZE_LOOP_3D_END;
	*avg = (long)t_avg;
	return(NORMAL_RC);
}


RC chk_graphic_range_mono8(GRAPHIC src, unsigned char *min, unsigned char *max,
                           unsigned char *avg)
{
	double t_avg = 0.0;
	double size = (src.size.z)*(src.size.y)*(src.size.x);
	if(size == 0) return(ARG_ERROR_RC);
	RC_NULL_CHK(src.data.mono8);
	*min = 255;
	*max = 0;
	GRAPHIC_SIZE_LOOP_3D(src.size.z, src.size.y, src.size.x)
	{
		unsigned char uc = src.data.mono8[ii1][ii2][ii3];
		if(uc < *min) *min = uc;
		if(uc > *max) *max = uc;
		t_avg += uc/size;
	}GRAPHIC_SIZE_LOOP_3D_END;
	*avg = (unsigned char)t_avg;
	return(NORMAL_RC);
}


RC chk_graphic_range_scalar(GRAPHIC src, double *min, double *max, double *avg)
{
	double t_avg = 0.0;
	double size = (src.size.z)*(src.size.y)*(src.size.x);
	if(size == 0) return(ARG_ERROR_RC);
	RC_NULL_CHK(src.data.scalar);
	*min = *max = src.data.scalar[0][0][0];
	GRAPHIC_SIZE_LOOP_3D(src.size.z, src.size.y, src.size.x)
	{
		double value = src.data.scalar[ii1][ii2][ii3];
		if(value < *min) *min = value;
		if(value > *max) *max = value;
		t_avg += value;
	}GRAPHIC_SIZE_LOOP_3D_END;
	*avg = t_avg/size;
	return(NORMAL_RC);
}


/*
 * $B852hA|(B(src)$B$O(B mono16 $BJQ498e2hA|(B(dest)$B$O(B rgb24
 * $B%b%N%/%m(B(16bit)$B$N:GBgCM:G>.CM$r4p$K!$$=$NI}$r%j%K%"$K(BRGB(8bit)$B$KJQ49$9$k(B
 * $B2hA|G;EY$NJ,I[$rL5;k$9$k$N$G!$(BCT$B$J$I$G$O$*A&$a$7$J$$(B
 */ 
RC conv_mono16_to_rgb24_1(GRAPHIC src, GRAPHIC *dest)
{
	long min, max, avg;

	if(src.type != GRAPHIC_MONO16){
		fprintf(stderr,"wrong graphic type.\n");
		return(ARG_ERROR_RC);
	}
	RC_NULL_CHK(src.data.mono16);
	RC_TRY( graphic_exist_chk(src) );
	RC_TRY( copy_graphic_info(src, dest) );
	dest->type = GRAPHIC_RGB24;
	RC_TRY( allocate_graphic_data(dest) );
	RC_TRY( chk_graphic_range_mono16(src, &min, &max, &avg) );

	GRAPHIC_SIZE_LOOP_3D(src.size.z, src.size.y, src.size.x){
		long cl = (255*(src.data.mono16[ii1][ii2][ii3]-min)) / (max-min);
		if((cl < 0) || (cl > 255)) return(CAL_ERROR_RC);
		dest->data.rgb24[ii1][ii2][ii3].r = (unsigned char)cl;
		dest->data.rgb24[ii1][ii2][ii3].g = (unsigned char)cl;
		dest->data.rgb24[ii1][ii2][ii3].b = (unsigned char)cl;
	}GRAPHIC_SIZE_LOOP_3D_END;

	return(NORMAL_RC);
}


/*
 * $B852hA|(B(src)$B$O(B mono16 $BJQ498e2hA|(B(dest)$B$O(B rgb24
 * $B%b%N%/%m(B(16bit)$B$N:GBgCM!$:G>.CM!$J?6QCM$r4p$K!$(B"$BJ?6QCM(B - $B:GBgCM(B"$B$b$7$/$O(B
 * "$BJ?6QCM(B - $B:G>.CM(B"$B$N:9$,>/$J$$J}$K(Bratio$B$r3]$1$?HO0O$rI=8=HO0O$H$9$k!%(B
 * $BJ?6Q$+$i!^(Bband$BHO0OFb$K$"$kCM$N$_$rBP>]$H$7$F$=$NCM$r%j%K%"$K(BRGB(8bit)$B$K(B
 * $BJQ49$9$k!%@55,J,I[$J$I$r;H$&$h$j$O9bB.!%(B
 */ 
RC conv_mono16_to_rgb24_2(GRAPHIC src, GRAPHIC *dest, double ratio)
{
	long min, max, avg;
	long band;

	if(src.type != GRAPHIC_MONO16){
		fprintf(stderr,"wrong graphic type.\n");
		return(ARG_ERROR_RC);
	}
	RC_NULL_CHK(src.data.mono16);
	RC_TRY( graphic_exist_chk(src) );
	RC_TRY( copy_graphic_info(src, dest) );
	dest->type = GRAPHIC_RGB24;	
	RC_TRY( allocate_graphic_data(dest) );
	RC_TRY( chk_graphic_range_mono16(src, &min, &max, &avg) );

	if(abs(max-avg) > abs(avg-min)) band = (long)ratio*(avg-min);
	else band = (long)ratio*(max-avg);

	min = avg - band;
	max = avg + band;
	band *= 2;
	GRAPHIC_SIZE_LOOP_3D(src.size.z, src.size.y, src.size.x){
		long cl = src.data.mono16[ii1][ii2][ii3];
		if(cl > max) cl = 255;
		else if(cl < min) cl = 0;
		else{
			cl = (255*(cl-min))/band;
			if((cl < 0) || (cl > 255)) return(CAL_ERROR_RC);
		}	
		dest->data.rgb24[ii1][ii2][ii3].r = (unsigned char)cl;
		dest->data.rgb24[ii1][ii2][ii3].g = (unsigned char)cl;
		dest->data.rgb24[ii1][ii2][ii3].b = (unsigned char)cl;
	}GRAPHIC_SIZE_LOOP_3D_END;

	return(NORMAL_RC);
}


/*
 * 3$B<!852hA|$NNX3TCj=P$r9T$&!%(Boperator $B$,(B NULL$B$N$H$-$O%G%U%)%k%H$N(B
 * Laplacian Operator$B$rMQ$$$k!%(Boperator$B$O(B3$B<!852hA|$J$i(B27$B!$(B2$B<!852hA|$J$i(B
 * 9$B8D$NCM$r;}$DG[Ns$G$J$1$l$P!$$=$NF0:n$OJ]>Z$7$J$$!%(B
 * $B8|$_$,(B3$B0J2<(B(size.z<3)$B$N;~!$(B2$B<!85$N2hA|$NNX3TCj=P$r7+$jJV$9!%(B
 * $B2hA|6-3&$NItJ,$O2hAG$N:GDcCM(B($BDL>o9u(B)$B$GEI$j$D$V$9!%(B
 * $B%G%U%)%k%H$N%*%Z%l!<%?$G$O2hA|$rFsCM2=$7$F$J$$$H$&$^$/$$$+$J$$!%(B
 */
RC extract_graphic_outline(GRAPHIC src, GRAPHIC *dest, const double *operator)
{
	double *op = NULL;
	int op_size;
	long lmin, lmax, lavg;
	double dmin, dmax, davg;
	double r, g, b;
	double cel;

	RC_TRY( graphic_exist_chk(src) );
	RC_TRY( copy_graphic_info(src, dest) );
	RC_TRY( allocate_graphic_data(dest) );

	if(operator == NULL){
		if(src.size.z < 3){
			op = laplacian_op;
			op_size = 9;
		}else{
			op = laplacian_op_3d;
			op_size = 27;
		}
	}else{
		if(src.size.z < 3) op_size = 9;
		else op_size = 27;
		*op = *operator;
	}

	if(src.type == GRAPHIC_MONO16) 
		RC_TRY( chk_graphic_range_mono16(src, &lmin, &lmax, &lavg) );

	if(src.type == GRAPHIC_SCALAR) 
		RC_TRY( chk_graphic_range_scalar(src, &dmin, &dmax, &davg) );

	GRAPHIC_SIZE_LOOP_3D(src.size.z, src.size.y, src.size.x){
		unsigned int *x, *y, *z;
		int ii4;
		int boundary = 0 ; /* 0: FLASE, 1: TRUE */

		RC_NULL_CHK( x = (unsigned int *)malloc(sizeof(unsigned int)*op_size) );
		RC_NULL_CHK( y = (unsigned int *)malloc(sizeof(unsigned int)*op_size) );
		RC_NULL_CHK( z = (unsigned int *)malloc(sizeof(unsigned int)*op_size) );

		if(op_size == 9){
			for(ii4=0; ii4<9; ii4++){
				if( ((operator_mat[ii4+9][0]+ii3) < 0) ||
				    ((operator_mat[ii4+9][1]+ii2) < 0) ){
					boundary = 1;
					break;
				}
				x[ii4] = operator_mat[ii4+9][0] + ii3;
				y[ii4] = operator_mat[ii4+9][1] + ii2;
				z[ii4] = ii1;
				if( (x[ii4]>src.size.x-1) || (y[ii4]>src.size.y-1) ){
					boundary = 1;
					break;
				}
			}
		}else{
			for(ii4=0; ii4<27; ii4++){
				if( ((operator_mat[ii4][0]+ii3) < 0) ||
				    ((operator_mat[ii4][1]+ii2) < 0) ||
				    ((operator_mat[ii4][2]+ii1) < 0) ){
					boundary = 1;
					break;
				}
				x[ii4] = operator_mat[ii4][0] + ii3;
				y[ii4] = operator_mat[ii4][1] + ii2;
				z[ii4] = operator_mat[ii4][2] + ii1;
				if( (x[ii4]>src.size.x-1) || (y[ii4]>src.size.y-1) ||
				    (z[ii4]>src.size.z-1) ){
					boundary = 1;
					break;
				}
			}
		}

		if(boundary){ /* $B6-3&$O:GDcCM$GEI$j$D$V$7(B */
			switch(src.type){
			case GRAPHIC_MONO8:
				dest->data.mono8[ii1][ii2][ii3] = 0;
				break;
			case GRAPHIC_MONO16:
				dest->data.mono16[ii1][ii2][ii3] = 0; 
				break;
			case GRAPHIC_RGB24:
				dest->data.rgb24[ii1][ii2][ii3] = init_graphic_rgb();
				break;
			case GRAPHIC_RGBA64:
				dest->data.rgba64[ii1][ii2][ii3] = init_graphic_rgba();
				break;
			case GRAPHIC_SCALAR:
				dest->data.scalar[ii1][ii2][ii3] = 0.0;
				break;
			default:
				fprintf(stderr, "unknown graphic type\n");
				break;
			}
			continue;
		}

		switch(src.type){
		case GRAPHIC_MONO8:
			cel = 0.0;
			for(ii4=0; ii4<op_size; ii4++)
				cel += op[ii4] * src.data.mono8[z[ii4]][y[ii4]][x[ii4]];
			if(cel > 255.0) cel = 255.0;
			if(cel < 0.0)   cel = 0.0;
			dest->data.mono8[ii1][ii2][ii3] = (unsigned char)cel;
			break;
		case GRAPHIC_MONO16:
			cel = 0.0;
			for(ii4=0; ii4<op_size; ii4++)
				cel += op[ii4] * src.data.mono16[z[ii4]][y[ii4]][x[ii4]];
			if(cel > lmax) cel = lmax;
			if(cel < lmin) cel = lmin;
			dest->data.mono16[ii1][ii2][ii3] = (long)cel;
			break;
		case GRAPHIC_SCALAR:
			cel = 0.0;
			for(ii4=0; ii4<op_size; ii4++)
				cel += op[ii4] * src.data.scalar[z[ii4]][y[ii4]][x[ii4]];
			if(cel > dmax) cel = dmax;
			if(cel < dmin) cel = dmin;
			dest->data.scalar[ii1][ii2][ii3] = cel;
			break;
		case GRAPHIC_RGB24:
			r = g = b = 0.0;
			for(ii4=0; ii4<op_size; ii4++){
				GRAPHIC_RGB rgb24 = src.data.rgb24[z[ii4]][y[ii4]][x[ii4]];
				r += op[ii4] * rgb24.r;
				g += op[ii4] * rgb24.g;
				b += op[ii4] * rgb24.b;
			}
			if(r > 255.0) r = 255.0;
			if(g > 255.0) g = 255.0;
			if(b > 255.0) b = 255.0;
			if(r < 0.0) r = 0.0;
			if(g < 0.0) g = 0.0;
			if(b < 0.0) b = 0.0;
			dest->data.rgb24[ii1][ii2][ii3].r = (unsigned char)r;
			dest->data.rgb24[ii1][ii2][ii3].g = (unsigned char)g;
			dest->data.rgb24[ii1][ii2][ii3].b = (unsigned char)b;
			break;
		case GRAPHIC_RGBA64:
			r = g = b = 0.0;
			for(ii4=0; ii4<op_size; ii4++){
				GRAPHIC_RGBA rgba64 = src.data.rgba64[z[ii4]][y[ii4]][x[ii4]];
				r += op[ii4] * rgba64.r;
				g += op[ii4] * rgba64.g;
				b += op[ii4] * rgba64.b;
			}
			if(r > 65535.0) r = 65535.0;
			if(g > 65535.0) g = 65535.0;
			if(b > 65535.0) b = 65535.0;
			if(r < 0.0) r = 0.0;
			if(g < 0.0) g = 0.0;
			if(b < 0.0) b = 0.0;
			dest->data.rgba64[ii1][ii2][ii3].r = (unsigned short)r;
			dest->data.rgba64[ii1][ii2][ii3].g = (unsigned short)g;
			dest->data.rgba64[ii1][ii2][ii3].b = (unsigned short)b;
			dest->data.rgba64[ii1][ii2][ii3].a
				= src.data.rgba64[ii1][ii2][ii3].a;
			break;
		case GRAPHIC_VECT2:
		case GRAPHIC_VECT3:
		case GRAPHIC_VECT6:
			fprintf(stderr, "non-support graphic type\n");
			return(ARG_ERROR_RC);
		default:
			fprintf(stderr, "unknown graphic type\n");
			return(ARG_ERROR_RC);
			break;
		}

		free(x);
		free(y);
		free(z);
	}GRAPHIC_SIZE_LOOP_3D_END;

	return(NORMAL_RC);
}


/* type $B$r85$K%5%$%:$N3NG'$H%G!<%?$,(B NULL$B$G$J$$$+%A%'%C%/$9$k!%(B*/
RC graphic_exist_chk(GRAPHIC gr)
{
	GRAPHIC_SIZE_CHK_3D(gr.size);

	switch(gr.type){
	case GRAPHIC_MONO8:
		RC_NULL_CHK(gr.data.mono8);
		break;
	case GRAPHIC_MONO16:
		RC_NULL_CHK(gr.data.mono16);
		break;
	case GRAPHIC_RGB24:
		RC_NULL_CHK(gr.data.rgb24);
		break;
	case GRAPHIC_RGBA64:
		RC_NULL_CHK(gr.data.rgba64);
		break;
	case GRAPHIC_SCALAR:
		RC_NULL_CHK(gr.data.scalar);
		break;
	case GRAPHIC_VECT2:
		RC_NULL_CHK(gr.data.vect2);
		break;
	case GRAPHIC_VECT3:
		RC_NULL_CHK(gr.data.vect3);
		break;
	case GRAPHIC_VECT6:
		RC_NULL_CHK(gr.data.vect6);
		break;
	case GRAPHIC_NOTYPE:
		fprintf(stderr, "graphic type does not set...\n");
		break;
	default:
		fprintf(stderr, "unknown graphic type\n");
		break;
	}

	return(NORMAL_RC);
}


/* $B2hAGCM$NH?E>$r9T$&!%(B*/
RC inverse_graphic(GRAPHIC src, GRAPHIC *dest)
{
	long lmin, lmax, lavg;
	double dmin, dmax, davg;

	RC_TRY( graphic_exist_chk(src) );
	RC_TRY( copy_graphic_info(src, dest) );
	RC_TRY( allocate_graphic_data(dest) );

	if(src.type == GRAPHIC_MONO16) 
		RC_TRY( chk_graphic_range_mono16(src, &lmin, &lmax, &lavg) );

	if(src.type == GRAPHIC_SCALAR) 
		RC_TRY( chk_graphic_range_scalar(src, &dmin, &dmax, &davg) );

	GRAPHIC_SIZE_LOOP_3D(src.size.z, src.size.y, src.size.x){
		switch(src.type){
		case GRAPHIC_MONO8:
			dest->data.mono8[ii1][ii2][ii3]
				= 255 - src.data.mono8[ii1][ii2][ii3];
			break;
		case GRAPHIC_MONO16:
			dest->data.mono16[ii1][ii2][ii3] 
				= lmax - src.data.mono16[ii1][ii2][ii3];
			break;
		case GRAPHIC_SCALAR:
			dest->data.scalar[ii1][ii2][ii3]
				= dmax - src.data.scalar[ii1][ii2][ii3];
			break;
		case GRAPHIC_RGB24:
			{
				GRAPHIC_RGB s_rgb = src.data.rgb24[ii1][ii2][ii3];
				GRAPHIC_RGB *d_rgb = &(dest)->data.rgb24[ii1][ii2][ii3];
				d_rgb->r = 255 - s_rgb.r;
				d_rgb->g = 255 - s_rgb.g;
				d_rgb->b = 255 - s_rgb.b;
			}
			break;
		case GRAPHIC_RGBA64:
			{
				GRAPHIC_RGBA s_rgba = src.data.rgba64[ii1][ii2][ii3];
				GRAPHIC_RGBA *d_rgba = &(dest)->data.rgba64[ii1][ii2][ii3];
				d_rgba->r = 65535 - s_rgba.r;
				d_rgba->g = 65535 - s_rgba.g;
				d_rgba->b = 65535 - s_rgba.b;
			}
			break;
		case GRAPHIC_VECT2:
		case GRAPHIC_VECT3:
		case GRAPHIC_VECT6:
			fprintf(stderr, "non-support graphic type\n");
			return(ARG_ERROR_RC);
		default:
			fprintf(stderr, "unknown graphic type\n");
			break;
		}
	}GRAPHIC_SIZE_LOOP_3D_END;

	return(NORMAL_RC);
}


/*
 * 3$B<!852hA|$NNX3T@~$+$iEy9b@~$KEv$?$kLOMM$rIA$/!%(B
 * $BNX3T@~$O2hAG$N:GDcCM$r;}$A!$B>$NItJ,$O:G9bCM$r;}$D$h$&$KFsCM2=$7$F$*$/(B
 * $B$3$H$,K>$^$7$$!%A*Br$5$l$?%3%s%?$NI=8=8B3&(B($B?'?t(B)$B$,>/$J$$>l9g!$<JLOMM$K(B
 * $B$J$j$d$9$$!%(Bsrc$B$O(B GRAPHIC_VECT?$B0J30$N2hA|!$(Bdest$B$O(B GRAPHIC_RGB24
 */
RC make_graphic_gradient_color_map(GRAPHIC src, GRAPHIC *dest,
                                   GRAPHIC_CONTOUR cnt,
                                   GRAPHIC_CONTOUR_TYPE type)
{
	GRAPHIC tmp;
	RC_TRY( make_graphic_gradient_scalar(src, &tmp) );
	RC_TRY( conv_graphic_scalar2color_map(tmp, dest, cnt, type) );
	free_graphic_data(tmp);

	return(NORMAL_RC);
}


/* $B2hA|$NNX3T@~$+$iEy9b@~$KEv$?$k%9%+%i!<CM(B(double$BCM(B)$B$N2hA|$rIA$/!%(B */
/* src$B$O(BGRAPHIC_VECT?$B0J30$N2hA|!$(Bdest$B$O(BGRAPHIC_SCALAR$B!%(B             */
RC make_graphic_gradient_scalar(GRAPHIC src, GRAPHIC *dest)
{
	GRAPHIC tmp;
	int exist_draw; 

	RC_TRY( graphic_exist_chk(src) );

	/* src($BG$0U(B)$B$N2hA|$r:GDcCM(B(0.0)$B$H$=$l0J30(B(-1.0)$B$K?6$jJ,$1$k(B */
	RC_TRY( conv_graphic_binary_scalar(src, dest) );

	/* $B:GDcCM6aJU$+$i1o<h$j$rKd$a$F$$$/(B */
	RC_TRY( copy_graphic_info(*dest, &tmp) );
	RC_TRY( allocate_graphic_data(&tmp) );
	for(exist_draw = 1; exist_draw; ){
		exist_draw = 0;

		GRAPHIC_SIZE_LOOP_3D(src.size.z, src.size.y, src.size.x){
			int ii4;
			unsigned int x, y, z;
			double min = DBL_MAX;
			if(dest->data.scalar[ii1][ii2][ii3] >= 0.0) continue;
			for(ii4=0; ii4<27; ii4++){
				if(ii4 == 13) continue; /* $BCf?42hAG(B */
				if( ((operator_mat[ii4][0]+ii3) < 0) ||
				    ((operator_mat[ii4][1]+ii2) < 0) ||
				    ((operator_mat[ii4][2]+ii1) < 0) ) continue; /* O.B. */
				x = operator_mat[ii4][0] + ii3;
				y = operator_mat[ii4][1] + ii2;
				z = operator_mat[ii4][2] + ii1;
				if((x > src.size.x-1)||(y > src.size.y-1)||(z > src.size.z-1)){
					 continue;/* O.B. */
				}
				if(dest->data.scalar[z][y][x] < 0.0){
					exist_draw = 1;
					continue;
				}
				min = MIN2(min, dest->data.scalar[z][y][x]);
			}
			if(min != DBL_MAX) tmp.data.scalar[ii1][ii2][ii3] = min + 2.0;
		}GRAPHIC_SIZE_LOOP_3D_END;

		GRAPHIC_SIZE_LOOP_3D(src.size.z, src.size.y, src.size.x){
			dest->data.scalar[ii1][ii2][ii3] += tmp.data.scalar[ii1][ii2][ii3];
			tmp.data.scalar[ii1][ii2][ii3] = 0.0;
		}GRAPHIC_SIZE_LOOP_3D_END;
	}
	free_graphic_data(tmp);

	return(NORMAL_RC);
}


/* src(GRAPHIC_SCALAR)$B$r(B $B;XDj$5$l$?%3%s%?$G(Bdest(GRAPHIC_RGB24)$B$KJQ49(B */
RC conv_graphic_scalar2color_map(GRAPHIC src, GRAPHIC *dest,
                                 GRAPHIC_CONTOUR cnt, GRAPHIC_CONTOUR_TYPE type)
{
	double min, max, avg;

	RC_TRY( graphic_exist_chk(src) );
	if(src.type != GRAPHIC_SCALAR){
		fprintf(stderr, "Illegal graphic type.\n");
		return(ARG_ERROR_RC);
	}
	RC_TRY( copy_graphic_info(src, dest) );
	dest->type = GRAPHIC_RGB24;
	RC_TRY( allocate_graphic_data(dest) );
	RC_TRY( chk_graphic_range_scalar(src, &min, &max, &avg) );

	GRAPHIC_SIZE_LOOP_3D(src.size.z, src.size.y, src.size.x){
		RC_TRY( get_graphic_contour_rgb24(src.data.scalar[ii1][ii2][ii3],
		                                  min, max, cnt, type,
		                                  &(dest)->data.rgb24[ii1][ii2][ii3]) );
	}GRAPHIC_SIZE_LOOP_3D_END;

	return(NORMAL_RC);
}


/* $B2hAG$NJ?6QCM$r4p$KDc$$J}$r(B0.0 $B$=$l0J30$NItJ,$O(B-1.0$B$N%9%+%i!<2hA|$K$9$k(B */
static RC conv_graphic_binary_scalar(GRAPHIC src, GRAPHIC *dest)
{
	unsigned char cmin, cmax, cavg;
	double dmin, dmax, davg;
	long lmin, lmax, lavg;
	GRAPHIC_RGB min3, max3, avg3;

	RC_TRY( graphic_exist_chk(src) );
	RC_TRY( copy_graphic_info(src, dest) );
	dest->type = GRAPHIC_SCALAR;
	RC_TRY( allocate_graphic_data(dest) );

	switch(src.type){
	case GRAPHIC_MONO8:
		RC_TRY( chk_graphic_range_mono8(src, &cmin, &cmax, &cavg) );
		GRAPHIC_SIZE_LOOP_3D(src.size.z, src.size.y, src.size.x){
			if(src.data.mono8[ii1][ii2][ii3] < cavg)
				dest->data.scalar[ii1][ii2][ii3] = 0.0;
			else
				dest->data.scalar[ii1][ii2][ii3] = -1.0;
		}GRAPHIC_SIZE_LOOP_3D_END;
		break;
	case GRAPHIC_MONO16:
		RC_TRY( chk_graphic_range_mono16(src, &lmin, &lmax, &lavg) );
		GRAPHIC_SIZE_LOOP_3D(src.size.z, src.size.y, src.size.x){
			if(src.data.mono16[ii1][ii2][ii3] < lavg)
				dest->data.scalar[ii1][ii2][ii3] = 0.0;
			else
				dest->data.scalar[ii1][ii2][ii3] = -1.0;
		}GRAPHIC_SIZE_LOOP_3D_END;
		break;
	case GRAPHIC_SCALAR:
		RC_TRY( chk_graphic_range_scalar(src, &dmin, &dmax, &davg) );
		GRAPHIC_SIZE_LOOP_3D(src.size.z, src.size.y, src.size.x){
			if(src.data.scalar[ii1][ii2][ii3] < davg)
				dest->data.scalar[ii1][ii2][ii3] = 0.0;
			else
				dest->data.scalar[ii1][ii2][ii3] = -1.0;
		}GRAPHIC_SIZE_LOOP_3D_END;
		break;
	case GRAPHIC_RGB24:
		RC_TRY( chk_graphic_range_rgb24(src, &min3, &max3, &avg3) );
		GRAPHIC_SIZE_LOOP_3D(src.size.z, src.size.y, src.size.x){
			if( (src.data.rgb24[ii1][ii2][ii3].r < (min3.r+max3.r)/2) &&
			    (src.data.rgb24[ii1][ii2][ii3].g < (min3.g+max3.g)/2) &&
			    (src.data.rgb24[ii1][ii2][ii3].b < (min3.b+max3.b)/2) )
				dest->data.scalar[ii1][ii2][ii3] = 0.0;
			else
				dest->data.scalar[ii1][ii2][ii3] = -1.0;
		}GRAPHIC_SIZE_LOOP_3D_END;
		break;
	case GRAPHIC_RGBA64: /* $B:GDcCM$r5a$a$k4X?t$,$J$$(B */
		GRAPHIC_SIZE_LOOP_3D(src.size.z, src.size.y, src.size.x){
			if( (src.data.rgba64[ii1][ii2][ii3].r == 0) &&
			    (src.data.rgba64[ii1][ii2][ii3].g == 0) &&
			    (src.data.rgba64[ii1][ii2][ii3].b == 0) )
				dest->data.scalar[ii1][ii2][ii3] = 0.0;
			else
				dest->data.scalar[ii1][ii2][ii3] = -1.0;
		}GRAPHIC_SIZE_LOOP_3D_END;
		break;
	case GRAPHIC_VECT2:
	case GRAPHIC_VECT3:
	case GRAPHIC_VECT6:
		fprintf(stderr, "non-support graphic type\n");
		return(ARG_ERROR_RC);
	default:
		fprintf(stderr, "unknown graphic type\n");
		return(ARG_ERROR_RC);
		break;
	}

	return(NORMAL_RC);
}


/*
 * node$B:BI8$K0lCW$9$k(Bsrc$B2hA|$N2hAGCM$rJV$94X?t!%(B
 * $B2hA|$H(BFEM$B$N:BI886E@$,0c$&$3$H$rG'<1$7$F$*$/$3$H!%(B0.0 =< node.p.? < 1.0
 * $B$G$"$k$H$-(B(0, Y_SIZE, 0)$B$N2hAGCM$,F@$i$l$k!%2hAG$+$i$O$_$@$7$F$$$k(Bnode
 * $B$G$O(B UNDERFLOW_ERROR_RC,$B$^$?$O!$(BOVERFLOW_ERROR_RC$B$,JV$5$l$k!%(B
 */
RC get_graphic_value_fem_node(const GRAPHIC src, const FEM_NODE node,
                              void *value)
{
	unsigned int x, y, z;

	RC_TRY( graphic_exist_chk(src) );

	if( (node.p.x < 0.0) || (node.p.y < 0.0) || (node.p.z < 0.0) )
		return(UNDERFLOW_ERROR_RC);

	x = node.p.x;
	y = (src.size.y - 1) - (unsigned int)node.p.y;
	z = node.p.z;

	if(src.size.z != 1){ /* src is voxel(3D) */
		if( (x > src.size.x-1) || (y > src.size.y-1) || (z > src.size.z-1) )
			return(OVERFLOW_ERROR_RC);
	}else{ /* src is pixel(2D) */
		if( (x > src.size.x-1) || (y > src.size.y-1) )
			return(OVERFLOW_ERROR_RC);
		z = 0;
	}

	switch(src.type){
	case GRAPHIC_MONO8:
		{
			unsigned char *c = (unsigned char *)value;
			*c = src.data.mono8[z][y][x];
		}
		break;
	case GRAPHIC_MONO16:
		{
			long *l = (long *)value;
			*l = src.data.mono16[z][y][x];
		}
		break;
	case GRAPHIC_SCALAR:
		{
			double *d = (double *)value;
			*d = (src).data.scalar[z][y][x];
		}
		break;
	case GRAPHIC_RGB24:
		{
			GRAPHIC_RGB *rgb = (GRAPHIC_RGB *)value;
			*rgb = src.data.rgb24[z][y][x];
		}
		break;
	case GRAPHIC_RGBA64:
		{
			GRAPHIC_RGBA *rgba = (GRAPHIC_RGBA *)value;
			*rgba = src.data.rgba64[z][y][x];
		}
		break;
	case GRAPHIC_VECT2:
		{
			TRANSLATION2D *p = (TRANSLATION2D *)value;
			*p = src.data.vect2[z][y][x];
		}
		break;
	case GRAPHIC_VECT3:
		{
			TRANSLATION3D *p = (TRANSLATION3D *)value;
			*p = src.data.vect3[z][y][x];
		}
		break;
	case GRAPHIC_VECT6:
		{
			TRANS_ROTATION3D *p = (TRANS_ROTATION3D *)value;
			*p = src.data.vect6[z][y][x];
		}
		break;
	default:
		fprintf(stderr, "unknown graphic type\n");
		return(ARG_ERROR_RC);
		break;
	}

	return(NORMAL_RC);
}

