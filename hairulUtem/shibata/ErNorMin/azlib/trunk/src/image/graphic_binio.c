/*********************************************************************
 * graphic_binio.c
 *
 * Copyright (C) 2004 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: graphic_binio.c 864 2007-01-17 10:02:04Z sasaoka $ */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "rc.h"
#include "memory_manager.h"
#include "graphic_utl.h"


#define FUNC_WRITE_GRAPHIC(func, src_type, data_type, img)\
RC func(FILE *fp, GRAPHIC gr)\
{\
    int ii1, ii2;\
    GRAPHIC_SIZE_CHK_3D(gr.size);\
    if(gr.type != src_type) return(ARG_ERROR_RC); \
    if( fwrite(&(gr), sizeof(GRAPHIC), 1, fp) != (size_t)1 )\
         return(WRITE_ERROR_RC);\
    if( fwrite(gr.data.img, sizeof(data_type **),\
               gr.size.z, fp) != (size_t)gr.size.z )\
        return(WRITE_ERROR_RC);\
    for(ii1=0; (unsigned long)ii1<gr.size.z; ii1++){\
        if( fwrite(gr.data.img[ii1], sizeof(data_type *),\
                   gr.size.y, fp) != (size_t)gr.size.y )\
            return(WRITE_ERROR_RC);\
        for(ii2=0; (unsigned long)ii2<gr.size.y; ii2++){\
             if( fwrite(gr.data.img[ii1][ii2], sizeof(data_type),\
                        gr.size.x, fp) != (size_t)gr.size.x )\
                 return(WRITE_ERROR_RC);\
        }\
    }\
    return(NORMAL_RC);\
}
FUNC_WRITE_GRAPHIC(write_graphic_rgb24,  GRAPHIC_RGB24,  GRAPHIC_RGB,   rgb24)
FUNC_WRITE_GRAPHIC(write_graphic_rgba64, GRAPHIC_RGBA64, GRAPHIC_RGBA,  rgba64)
FUNC_WRITE_GRAPHIC(write_graphic_mono8,  GRAPHIC_MONO8,  unsigned char, mono8)
FUNC_WRITE_GRAPHIC(write_graphic_mono16, GRAPHIC_MONO16, long,          mono16)
FUNC_WRITE_GRAPHIC(write_graphic_scalar, GRAPHIC_SCALAR, double,        scalar)



#define FUNC_READ_GRAPHIC(func, src_type, data_type, img)\
RC func(FILE *fp, GRAPHIC *gr)\
{\
    int ii1, ii2;\
    if( fread(gr, sizeof(GRAPHIC), 1, fp) != (size_t)1 )\
        return(READ_ERROR_RC);\
    if(gr->type != src_type) return(ARG_ERROR_RC);\
    gr->data.img = (data_type ***)mm_alloc(gr->size.z*sizeof(data_type **));\
    if(gr->data.img == NULL){\
        RC_TRY( error_printf(9001, gr->size.z*sizeof(data_type **)) );\
        return(ALLOC_ERROR_RC);\
    }\
    if( fread(gr->data.img, sizeof(data_type **),\
              gr->size.z, fp) != (size_t)gr->size.z ){\
        return(READ_ERROR_RC);\
    }\
    for(ii1=0; (unsigned long)ii1<gr->size.z; ii1++){\
        (gr->data.img)[ii1]\
            = (data_type **)mm_alloc(gr->size.y*sizeof(data_type *));\
        if((gr->data.img)[ii1] == NULL){\
            RC_TRY( error_printf(9001, gr->size.y*sizeof(data_type *)) );\
            return(ALLOC_ERROR_RC);\
        }\
        if( fread((gr->data.img)[ii1], sizeof(data_type *),\
                   gr->size.y, fp) != (size_t)gr->size.y ){\
            return(READ_ERROR_RC);\
        }\
        for(ii2=0; (unsigned long)ii2<gr->size.y; ii2++){\
            (gr->data.img)[ii1][ii2]\
                 = (data_type *)mm_alloc(gr->size.x*sizeof(data_type));\
            if((gr->data.img)[ii1][ii2] == NULL){\
                RC_TRY( error_printf(9001, gr->size.x*sizeof(data_type)) );\
                return(ALLOC_ERROR_RC);\
            }\
            if( fread((gr->data.img)[ii1][ii2], sizeof(data_type),\
                       gr->size.x, fp) != (size_t)gr->size.x ){\
                return(READ_ERROR_RC);\
            }\
        }\
    }\
    return(NORMAL_RC);\
}
FUNC_READ_GRAPHIC(read_graphic_rgb24,  GRAPHIC_RGB24,  GRAPHIC_RGB,   rgb24)
FUNC_READ_GRAPHIC(read_graphic_rgba64, GRAPHIC_RGBA64, GRAPHIC_RGBA,  rgba64)
FUNC_READ_GRAPHIC(read_graphic_mono8,  GRAPHIC_MONO8,  unsigned char, mono8)
FUNC_READ_GRAPHIC(read_graphic_mono16, GRAPHIC_MONO16, long,          mono16)
FUNC_READ_GRAPHIC(read_graphic_scalar, GRAPHIC_SCALAR, double,        scalar)


