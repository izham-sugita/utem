/*********************************************************************
 * voxel_utl.c
 *
 * Copyright (C) 2004 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: voxel_utl.c 432 2005-08-24 03:00:31Z sasaoka $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rc.h"
#include "log_printf.h"
#include "memory_manager.h"
#include "dicom_utl.h"
#include "graphic_utl.h"

static RC read_one_file(GRAPHIC *gr, char *file, GRAPHIC_SOURCE source);

RC
make_voxel_from_file (GRAPHIC *gr, GRAPHIC_SOURCE source, char *files[], 
                     int file_num, double thickness, double resolution,
                     int fill)
{
	int ii1, z_pos;
	GRAPHIC test;
	void **data;

	/* files read test */
	RC_TRY( read_one_file(&test, files[0], source) );

	/* GRAPHIC parameter setting */
	gr->size.x = test.size.x;
	gr->size.y = test.size.y;
	if(test.resolution.x != test.resolution.y)
		RC_TRY( error_printf(7004) );
	gr->size.z = (unsigned long)(file_num * thickness / resolution);
	gr->type = test.type;
	gr->resolution.x = test.resolution.x;
	gr->resolution.y = test.resolution.y;
	gr->resolution.z = resolution;
	gr->source = test.source;

	/* free open test data */
	free_graphic_data(test);

	/* import all files */
	data = (void **)mm_alloc(gr->size.z * sizeof(void *));
	if(data == NULL) return(ALLOC_ERROR_RC);

	for(z_pos=0,ii1=0; ii1<file_num; ii1++){
		GRAPHIC target;
		RC_TRY( read_one_file(&target, files[ii1], source) );
		RC_TRY( log_printf(4, "file[ii1] = %d\n", ii1) );
		GRAPHIC_SIZE_CHK_2D(target.size);
		RC_TRY( graphic_exist_chk(target) );

		while(z_pos*resolution/thickness < (double)(ii1+1)){
			GRAPHIC *copy = NULL;
			if(fill){
				copy = (GRAPHIC *)mm_alloc(sizeof(GRAPHIC));
				if(data == NULL) return(ALLOC_ERROR_RC);
				RC_TRY( copy_graphic(target, copy) );
			}else{
				copy = &target;
			}
			switch(target.type){
			case GRAPHIC_RGB24:
				data[z_pos] = (void *)copy->data.rgb24[0];
				break;
			case GRAPHIC_RGBA64:
				data[z_pos] = (void *)copy->data.rgba64[0];
				break;
			case GRAPHIC_MONO8:
				data[z_pos] = (void *)copy->data.mono8[0];
				break;
			case GRAPHIC_MONO16:
				data[z_pos] = (void *)copy->data.mono16[0];
				break;
			case GRAPHIC_SCALAR:
				data[z_pos] = (void *)copy->data.scalar[0];
				break;
			default:
				return(IMPLEMENT_ERROR_RC);
			}
			z_pos++;
			if((unsigned long)z_pos > (gr->size.z-1)) break;
		}
		if(fill) RC_TRY( free_graphic_data(target) );
	}

	switch(gr->type){
	case GRAPHIC_RGB24:
		gr->data.rgb24 = (GRAPHIC_RGB ***)data;
		break;
	case GRAPHIC_RGBA64:
		gr->data.rgba64 = (GRAPHIC_RGBA ***)data;
		break;
	case GRAPHIC_MONO8:
		gr->data.mono8 = (unsigned char ***)data;
		break;
	case GRAPHIC_MONO16:
		gr->data.mono16 = (long ***)data;
		break;
	case GRAPHIC_SCALAR:
		gr->data.scalar = (double ***)data;
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
read_one_file (GRAPHIC *gr, char *file, GRAPHIC_SOURCE source)
{
	FILE *fp;

	RC_TRY( rc_fopen(file, "r", &fp) );
	switch(source){
	case GRAPHIC_BMP:
		RC_TRY( log_printf(4, "Read BMP file : %s\n", file) );
		RC_TRY( read_bmp(fp, gr) );
		break;
	case GRAPHIC_TIFF:
		RC_TRY( log_printf(4, "Read TIFF file : %s\n", file) );
		RC_TRY( read_tiff(fp, gr) );
		break;
	case GRAPHIC_DICOM:
		{
			DICOM_HEADER dicom;
			RC_TRY( log_printf(4, "Read DICOM file : %s\n", file) );
			RC_TRY( read_dicom(fp, &dicom, gr) );
			/* 2の補数表現をしている */
			if( (gr->type == GRAPHIC_MONO16) && (dicom.pixel_representation) ){
				GRAPHIC conv;
				RC_TRY(log_printf(4, "This DICOM file is complement"
				                     "representation\n") );
				RC_TRY( conv_2complement_to_mono16(*gr, &conv, dicom) );
				RC_TRY( free_graphic_data(*gr) );
				gr->data.mono16 = conv.data.mono16;
			}
		}
		break;
	case GRAPHIC_AVS_FIELD:
		RC_TRY( log_printf(4, "Read AVS(FIELD DATA) file : %s\n", file) );
		RC_TRY( read_avs_field(fp, gr) );
		break;
	case GRAPHIC_NOSRC:
		return(ARG_ERROR_RC);
	default:
		return(IMPLEMENT_ERROR_RC);
	}
	RC_TRY( rc_fclose(fp) );

	return(NORMAL_RC);
}

