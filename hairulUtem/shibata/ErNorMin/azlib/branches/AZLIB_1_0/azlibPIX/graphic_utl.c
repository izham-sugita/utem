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

/* サイズなどの情報は先にセットされているものとする */
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


/* Blue Green Red の順で値が高くなる 1021階調 */
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


/* Blue Red の順で値が高くなる 511階調 */
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


/* Red Yeloow の順で値が高くなる 256階調 */
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


/* White Blue の順で値が高くなる 256階調 */
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


/* Black White の順で値が高くなる 256階調 */
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
/* mix < value < maxであるときの，valueのカラーコンタ色(color)を返す．*/
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


/* colorであたえられた色のカラーコンタ上の番号(大きさ)を返す．*/
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
 * 輝度を元にRGBの最大値を考える．平均は各色の平均である．
 * CIE(国際照明委員会)XYZ表色系 Yが輝度にあたる
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
 * 元画像(src)は mono16 変換後画像(dest)は rgb24
 * モノクロ(16bit)の最大値最小値を基に，その幅をリニアにRGB(8bit)に変換する
 * 画像濃度の分布を無視するので，CTなどではお薦めしない
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
 * 元画像(src)は mono16 変換後画像(dest)は rgb24
 * モノクロ(16bit)の最大値，最小値，平均値を基に，"平均値 - 最大値"もしくは
 * "平均値 - 最小値"の差が少ない方にratioを掛けた範囲を表現範囲とする．
 * 平均から±band範囲内にある値のみを対象としてその値をリニアにRGB(8bit)に
 * 変換する．正規分布などを使うよりは高速．
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
 * 3次元画像の輪郭抽出を行う．operator が NULLのときはデフォルトの
 * Laplacian Operatorを用いる．operatorは3次元画像なら27，2次元画像なら
 * 9個の値を持つ配列でなければ，その動作は保証しない．
 * 厚みが3以下(size.z<3)の時，2次元の画像の輪郭抽出を繰り返す．
 * 画像境界の部分は画素の最低値(通常黒)で塗りつぶす．
 * デフォルトのオペレータでは画像を二値化してないとうまくいかない．
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

		if(boundary){ /* 境界は最低値で塗りつぶし */
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


/* type を元にサイズの確認とデータが NULLでないかチェックする．*/
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


/* 画素値の反転を行う．*/
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
 * 3次元画像の輪郭線から等高線に当たる模様を描く．
 * 輪郭線は画素の最低値を持ち，他の部分は最高値を持つように二値化しておく
 * ことが望ましい．選択されたコンタの表現限界(色数)が少ない場合，縞模様に
 * なりやすい．srcは GRAPHIC_VECT?以外の画像，destは GRAPHIC_RGB24
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


/* 画像の輪郭線から等高線に当たるスカラー値(double値)の画像を描く． */
/* srcはGRAPHIC_VECT?以外の画像，destはGRAPHIC_SCALAR．             */
RC make_graphic_gradient_scalar(GRAPHIC src, GRAPHIC *dest)
{
	GRAPHIC tmp;
	int exist_draw; 

	RC_TRY( graphic_exist_chk(src) );

	/* src(任意)の画像を最低値(0.0)とそれ以外(-1.0)に振り分ける */
	RC_TRY( conv_graphic_binary_scalar(src, dest) );

	/* 最低値近辺から縁取りを埋めていく */
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
				if(ii4 == 13) continue; /* 中心画素 */
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


/* src(GRAPHIC_SCALAR)を 指定されたコンタでdest(GRAPHIC_RGB24)に変換 */
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


/* 画素の平均値を基に低い方を0.0 それ以外の部分は-1.0のスカラー画像にする */
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
	case GRAPHIC_RGBA64: /* 最低値を求める関数がない */
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
 * node座標に一致するsrc画像の画素値を返す関数．
 * 画像とFEMの座標原点が違うことを認識しておくこと．0.0 =< node.p.? < 1.0
 * であるとき(0, Y_SIZE, 0)の画素値が得られる．画素からはみだしているnode
 * では UNDERFLOW_ERROR_RC,または，OVERFLOW_ERROR_RCが返される．
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

