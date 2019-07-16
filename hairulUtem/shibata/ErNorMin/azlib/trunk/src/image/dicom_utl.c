/*********************************************************************
 * dicom_utl.c
 *
 * Copyright (C) 2001 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA> <Yoshihito KURIZUKA> <Kazuma ADACHI>
 *
 * Thanks to following contributors.
 *  <Kenzen TAKEUCHI> <Tomoyuki TSUBATA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: dicom_utl.c 864 2007-01-17 10:02:04Z sasaoka $ */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "rc.h"
#include "memory_manager.h"
#include "bin_io.h"
#include "dicom_utl.h"
#include "dicom_element.h"
#include "dicom_tag.h"
#include "graphic_utl.h"

#define BUFFER_SIZE (512)

static RC read_tag(FILE *fp, DICOM_HEADER *dicom, unsigned long *ul);
static RC read_tag_elem(FILE *fp, char tmp[], unsigned long tag,
                        int vr_flag);
static RC skip_tag_elemOB(FILE *fp);
static RC skip_tag_elemAE(FILE *fp, int vr_flag);
static RC skip_tag_elemNONE(FILE *fp);
static int search_dicom_tag(const int size, const unsigned long target);
static RC check_vr(FILE *fp, int *flag);

static int dicom_elem_size = (sizeof(Dicom_elem)/sizeof(DICOM_ELEM));

RC
read_dicom (FILE *fp, DICOM_HEADER *dicom, GRAPHIC *gr)
{
	int pimage;

	init_dicom(dicom);
	init_graphic(gr);
	RC_TRY( read_dicom_header(fp, dicom, &pimage) );
	RC_TRY( set_dicom_header(*dicom, gr) );
	RC_TRY( read_dicom_image(fp, *dicom, gr, pimage) );
	/*
	print_dicom_header(stderr, *dicom);
	RC_TRY( print_dicom_tag(stderr, *dicom) );
	*/
	return(NORMAL_RC);
}

RC
read_dicom_header (FILE *fp, DICOM_HEADER *dicom, int *pimage)
{
	unsigned long ul;
	unsigned int ui;
	unsigned int vr;
	int index;
	char ctmp[DICOM_STR_SIZE];
	char c;
	float ftmp, ftmp1, ftmp2;
	int vr_flag ;
	RC rc;

	rc = read_tag(fp, dicom, &ul);
	if(rc == END_RC){
		return(NORMAL_RC);
	}else if(rc != NORMAL_RC){
		return(rc);
	}

	if(ul == Group_0000_Length){
		if(fseek(fp, 128L, SEEK_CUR) != 0) return(SEEK_ERROR_RC);
	}else{
		rewind(fp);
	}

	while(1){
		rc = read_tag(fp, dicom, &ul);
		if(rc == END_RC){
			return(NORMAL_RC);
		}else if(rc != NORMAL_RC){
			return(rc);
		}
		RC_TRY( check_vr(fp, &vr_flag) );

		switch(ul){
		case Group_0008_Length_to_End:
			RC_TRY( lread_4bytes(fp, &ul) );
			RC_TRY( lread_4bytes(fp, &ul) );
			dicom->length_to_end = ul;
			break;

		case Study_Date: 
			RC_TRY( read_tag_elem(fp, dicom->study_date, ul, vr_flag) );
			break;

		case Study_Time:
			RC_TRY( read_tag_elem(fp, dicom->study_time, ul, vr_flag) );
			break;

		case Modality:
			RC_TRY( read_tag_elem(fp, dicom->modality, ul, vr_flag) );
			break;

		case Manufacturer:
			RC_TRY( read_tag_elem(fp, dicom->manufacturer, ul, vr_flag) );
			break;

		case Station_Name:
			RC_TRY( read_tag_elem(fp, dicom->station_name, ul, vr_flag) );
			break;

		case Patients_Sex:
			RC_TRY( read_tag_elem(fp, dicom->patients_sex, ul, vr_flag) );
			break;

		case Patients_Age:
			RC_TRY( read_tag_elem(fp, dicom->patients_age, ul, vr_flag) );
			break;
			
		case Slice_Thickness:
			RC_TRY( read_tag_elem(fp, ctmp, ul, vr_flag) );
			if(sscanf(ctmp, "%f", &ftmp) != 1) return(READ_ERROR_RC);
			dicom->slice_thickness = (double)ftmp;
			break;

		case Gantry_Detector_Tilt:
			RC_TRY( read_tag_elem(fp, ctmp, ul, vr_flag) );
			if(sscanf(ctmp, "%f", &ftmp) != 1) return(READ_ERROR_RC);
		    dicom->gantry_detector = (double)ftmp;
			break;

		case Table_Height:
			RC_TRY( read_tag_elem(fp, ctmp, ul, vr_flag) );
			if(sscanf(ctmp, "%f", &ftmp) != 1) return(READ_ERROR_RC);
		    dicom->table_height = (double)ftmp;
			break;

		case Rotation_Direction:
			RC_TRY( read_tag_elem(fp, dicom->rotation_direction,
			                      ul, vr_flag) );
			break;

		case Patient_Position:
			RC_TRY( read_tag_elem(fp, dicom->pation_position, ul, vr_flag) );
			break;

		case Patient_Orientation:
			RC_TRY( read_tag_elem(fp, dicom->patient_orientation,
			                      ul,vr_flag) );
			break;

		case Smallest_Image_Pixel_Value:
			RC_TRY( lread_2bytes(fp, &ui) );
			RC_TRY( lread_2bytes(fp, &ui) );
			RC_TRY( lread_2bytes(fp, &ui) );
			dicom->smallest_image_pixel_value = (unsigned long)ui;
			break;

		case Largest_Image_Pixel_Value:
			RC_TRY( lread_2bytes(fp, &ui) );
			RC_TRY( lread_2bytes(fp, &ui) );
			RC_TRY( lread_2bytes(fp, &ui) );
			dicom->largest_image_pixel_value = (unsigned long)ui;
			break;

		case Slice_Location:
			RC_TRY( read_tag_elem(fp, ctmp, ul, vr_flag) );
			if(sscanf(ctmp, "%f", &ftmp) != 1) return(READ_ERROR_RC);
			dicom->slice_location = (double)ftmp;
			break;

		case Samples_per_Pixel:
			if(vr_flag == 1){
				RC_TRY( lread_2bytes(fp, &vr) );
				RC_TRY( lread_2bytes(fp, &ui) );
				ul = (unsigned long)ui;
			}else if(vr_flag == 0){
				RC_TRY( lread_4bytes(fp, &ul) );
			}else{
				return(UNKNOWN_ERROR_RC);
			}
			if(ul == 2){
				RC_TRY( lread_2bytes(fp, &ui) );
				dicom->samples_per_pixel = (unsigned long)ui;
			}else if(ul == 4){
				RC_TRY( lread_4bytes(fp, &ul) );
				dicom->samples_per_pixel = ul;
			}else{
				return(READ_ERROR_RC);
			}
			break;

		case Photometric_Interpretation:
			RC_TRY( read_tag_elem(fp, dicom->photometric_interpretation,
			                      ul, vr_flag) );
			break;

		case Planar_Configuration:
			if(vr_flag == 1){
				RC_TRY( lread_2bytes(fp, &vr) );
				RC_TRY( lread_2bytes(fp, &ui) );
				ul = (unsigned long)ui;
			}else if(vr_flag == 0){
				RC_TRY( lread_4bytes(fp, &ul) );
			}else{
				return(UNKNOWN_ERROR_RC);
			}
			if(ul == 2){
				RC_TRY( lread_2bytes(fp, &ui) );
				dicom->planar_configuration = ui;
			}else if(ul == 4){
				RC_TRY( lread_4bytes(fp, &ul) );
				dicom->planar_configuration = (unsigned int)ul;
			}else{
				return(READ_ERROR_RC);
			}
			break;

		case Rows:
			if(vr_flag == 1){
				RC_TRY( lread_2bytes(fp, &vr) );
				RC_TRY( lread_2bytes(fp, &ui) );
				ul = (unsigned long)ui;
			}else{
				RC_TRY( lread_4bytes(fp, &ul) );
			}

			if(ul == 2){
				RC_TRY( lread_2bytes(fp, &ui) );
				dicom->rows = (unsigned long)ui;
			}else if(ul == 4){
				RC_TRY( lread_4bytes(fp, &ul) );
				dicom->rows = ul;
			}else{
				return(READ_ERROR_RC);
			}
			break;

		case Columns:
			if(vr_flag == 1){
				RC_TRY( lread_2bytes(fp, &vr) );
				RC_TRY( lread_2bytes(fp, &ui) );
				ul = (unsigned long)ui;
			}else{
				RC_TRY( lread_4bytes(fp, &ul) );
			}
			
			if(ul == 2){
				RC_TRY( lread_2bytes(fp, &ui) );
				dicom->columns = (unsigned long)ui;
			}else if(ul == 4){
				RC_TRY( lread_4bytes(fp, &ul) );
				dicom->columns = ul;
			}else{
				return(READ_ERROR_RC);
			}
			break;

		case Pixel_Spacing:
			RC_TRY( read_tag_elem(fp, ctmp, ul, vr_flag) );
			if(sscanf(ctmp, "%f%c%f", &ftmp1, &c, &ftmp2) != 3){
				return(READ_ERROR_RC);
			}
			dicom->pixcel_spacing_x = (double)ftmp1;
			dicom->pixcel_spacing_y = (double)ftmp2;
			break;

		case Pixel_Aspect_Ratio:
			RC_TRY( read_tag_elem(fp, ctmp, ul, vr_flag) );
			break;

		case Bits_Allocated:
			if(vr_flag == 1){
			RC_TRY( lread_2bytes(fp, &vr) );
				RC_TRY( lread_2bytes(fp, &ui) );
				ul = (unsigned long)ui;
			}else{
				RC_TRY( lread_4bytes(fp, &ul) );
			}

			if(ul == 2){
				RC_TRY( lread_2bytes(fp, &ui) );
				dicom->bits_allocated = (unsigned long)ui;
			}else if(ul == 4){
				RC_TRY( lread_4bytes(fp, &ul) );
				dicom->bits_allocated = ul;
			}else{
				return(READ_ERROR_RC);
			}
			break;

		case Bits_Stored:
			if(vr_flag == 1){
				RC_TRY( lread_2bytes(fp, &vr) );
				RC_TRY( lread_2bytes(fp, &ui) );
				ul = (unsigned long)ui;
			}else{
				RC_TRY( lread_4bytes(fp, &ul) );
			}

			if(ul == 2){
				RC_TRY( lread_2bytes(fp, &ui) );
				dicom->bits_stored = (unsigned long)ui;
			}else if(ul == 4){
				RC_TRY( lread_4bytes(fp, &ul) );
				dicom->bits_stored = ul;
			}else{
				return(READ_ERROR_RC);
			}
			break;

		case Pixel_Representation:
			if(vr_flag == 1){
				RC_TRY( lread_2bytes(fp, &vr) );
				RC_TRY( lread_2bytes(fp, &ui) );
				ul = (unsigned long)ui;
			}else{
				RC_TRY( lread_4bytes(fp, &ul) );
			}

			if(ul == 2){
				RC_TRY( lread_2bytes(fp, &ui) );
				dicom->pixel_representation = (unsigned long)ui;
			}else if(ul == 4){
				RC_TRY( lread_4bytes(fp, &ul) );
				dicom->pixel_representation = ul;
			}else{
				return(READ_ERROR_RC);
			}
			break;

		case Pixel_Data:
			if(vr_flag == 1){
				RC_TRY( lread_2bytes(fp, &vr) );
				RC_TRY( lread_2bytes(fp, &ui) );
				RC_TRY( lread_4bytes(fp, &ul) );
			}else{
				RC_TRY( lread_4bytes(fp, &ul) );
			}
			
			dicom->pixel_data = ul;
			*pimage = ftell(fp);
			if(fseek(fp, ul, SEEK_CUR) != 0){
				return(SEEK_ERROR_RC);
			}
			break;

		default:
			index = search_dicom_tag(dicom_elem_size, ul);
			if(index >= 0){
				if( (Dicom_elem[index].vr == OB) ||
					(Dicom_elem[index].vr == OW) ||
					(Dicom_elem[index].vr == SQ) ||
					(Dicom_elem[index].vr == UN)){
					RC_TRY( skip_tag_elemOB(fp) );
				}else if(Dicom_elem[index].vr == NONE){
					RC_TRY( skip_tag_elemNONE(fp) );
				}else{
					RC_TRY( skip_tag_elemAE(fp, vr_flag) );
				}
			}else{
				/* 奇数番号を持つ私的データ要素に対する暫定措置 */
				if(((ul>>16) % 2) != 0){
					RC_TRY( skip_tag_elemAE(fp, vr_flag) );
				/* 不明番号タグを表示して終了 */
				}else{
					RC_TRY( error_printf(7301, ul) );
					return(IMPLEMENT_ERROR_RC);
				}
			}
			break;
		}
	}

	return(NORMAL_RC);
}


RC
set_dicom_header (DICOM_HEADER dicom, GRAPHIC *gr)
{
	if((dicom.rows < 1) || (dicom.columns < 1)) return(READ_ERROR_RC);
	
	gr->source = GRAPHIC_DICOM;
	gr->size.x = dicom.columns;
	gr->size.y = dicom.rows;
	gr->size.z = 1;
	gr->resolution.x = dicom.pixcel_spacing_x;
	gr->resolution.y = dicom.pixcel_spacing_y;
	gr->resolution.z = 0.0;
	
	if(strncmp(dicom.photometric_interpretation, "MONOCHROME1", 11) == 0){
		if(dicom.bits_stored == 16){
			gr->type = GRAPHIC_MONO16;
			gr->i_info[0] = 0;
		}else{
			RC_TRY( error_printf(7001) );
			return(IMPLEMENT_ERROR_RC);
		}
	}else if(strncmp(dicom.photometric_interpretation ,"MONOCHROME2", 11) == 0){
		if(dicom.bits_stored == 16){
			gr->type = GRAPHIC_MONO16;
			gr->i_info[0] = 1;
		}else{
			RC_TRY( error_printf(7001) );
			return(IMPLEMENT_ERROR_RC);
		}
	}else if(strncmp(dicom.photometric_interpretation,"PALETTE COLOR",13) == 0){
		if(dicom.bits_stored == 8){
			gr->type = GRAPHIC_RGB24;
			gr->i_info[0] = 2;
		}else{
			RC_TRY( error_printf(7001) );
			return(IMPLEMENT_ERROR_RC);
		}
	}else if(strncmp(dicom.photometric_interpretation ,"RGB", 3) == 0){
		if(dicom.bits_stored == 8){
			gr->type = GRAPHIC_RGB24;
			gr->i_info[0] = 3;
		}else{
			RC_TRY( error_printf(7001) );
			return(IMPLEMENT_ERROR_RC);
		}
	}else if(strncmp(dicom.photometric_interpretation ,"HSV", 3) == 0){
		if(dicom.bits_stored == 8){
			gr->type = GRAPHIC_RGB24;
			gr->i_info[0] = 4;
		}else{
			RC_TRY( error_printf(7001) );
			return(IMPLEMENT_ERROR_RC);
		}
	}else if(strncmp(dicom.photometric_interpretation ,"ARGB", 4) == 0){
		if(dicom.bits_stored == 8){
			gr->type = GRAPHIC_RGB24;
			gr->i_info[0] = 5;
		}else{
			RC_TRY( error_printf(7001) );
			return(IMPLEMENT_ERROR_RC);
		}
	}else if(strncmp(dicom.photometric_interpretation ,"CMYK", 4) == 0){
		if(dicom.bits_stored == 8){
			gr->type = GRAPHIC_RGB24;
			gr->i_info[0] = 6;
		}else{
			RC_TRY( error_printf(7001) );
			return(IMPLEMENT_ERROR_RC);
		}
	}else if(strncmp(dicom.photometric_interpretation ,"YBR_FULL", 8) == 0){
		if(dicom.bits_stored == 8){
			gr->type = GRAPHIC_RGB24;
			gr->i_info[0] = 7;
		}else{
			RC_TRY( error_printf(7001) );
			return(IMPLEMENT_ERROR_RC);
		}
	}else if(strncmp(dicom.photometric_interpretation ,"YBR_FULL_422", 12)== 0){
		if(dicom.bits_stored == 8){
			gr->type = GRAPHIC_RGB24;
			gr->i_info[0] = 8;
		}else{
			RC_TRY( error_printf(7001) );
			return(IMPLEMENT_ERROR_RC);
		}
	}else if(strncmp(dicom.photometric_interpretation , "\0", 1)== 0){
		if(dicom.bits_stored == 16){
			gr->type = GRAPHIC_MONO16;
			gr->i_info[0] = 9;
		}else{
			RC_TRY( error_printf(7001) );
			return(IMPLEMENT_ERROR_RC);
		}
	}else{
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


RC
read_dicom_image (FILE *fp, DICOM_HEADER dicom, GRAPHIC *gr, int pimage)
{
	unsigned int ui;
	int cl;
	
	if(pimage <= 0) return(READ_ERROR_RC);
	if(fseek(fp, pimage, SEEK_SET) != 0) return(SEEK_ERROR_RC);

	if((gr->type == GRAPHIC_RGB24) || (gr->type == GRAPHIC_MONO16)){
		if(gr != NULL) RC_TRY( allocate_graphic_data(gr) );
	}else{
		RC_TRY( error_printf(7001) );
		return(IMPLEMENT_ERROR_RC);
	}
		
	switch(gr->i_info[0]){
	case 1:
	case 9:
		if((dicom.bits_stored == 8) && (dicom.bits_allocated == 8)){
			GRAPHIC_SIZE_LOOP_2D(gr->size.y, gr->size.x){
				if((cl = fgetc(fp)) == EOF) return(READ_ERROR_RC);
				gr->data.mono16[0][ii1][ii2] = (long)cl;
			}GRAPHIC_SIZE_LOOP_2D_END;
		}else if((dicom.bits_stored == 16) && (dicom.bits_allocated == 16)){
			GRAPHIC_SIZE_LOOP_2D(gr->size.y, gr->size.x){
				RC_TRY( lread_2bytes(fp, &ui) );
				gr->data.mono16[0][ii1][ii2] = (long)ui;
			}GRAPHIC_SIZE_LOOP_2D_END;
		}else{
			RC_TRY( error_printf(7001) );
			return(READ_ERROR_RC);
		}
		break;
	case 3:
		if((dicom.bits_stored == 8) && (dicom.bits_allocated == 8)){
			if(dicom.planar_configuration == 0){ 
				GRAPHIC_SIZE_LOOP_2D(gr->size.y, gr->size.x){
				  if( (cl=fgetc(fp)) == EOF ) return(READ_ERROR_RC);
				  gr->data.rgb24[0][ii1][ii2].r = (unsigned char)cl;
				  if( (cl=fgetc(fp)) == EOF ) return(READ_ERROR_RC);
				  gr->data.rgb24[0][ii1][ii2].g = (unsigned char)cl;
				  if( (cl=fgetc(fp)) == EOF ) return(READ_ERROR_RC);
				  gr->data.rgb24[0][ii1][ii2].b = (unsigned char)cl;
				}GRAPHIC_SIZE_LOOP_2D_END;
			}else if(dicom.planar_configuration > 0){
				int ii1, ii2;
				for(ii1=0; (unsigned long)ii1<gr->size.y; ii1++){
					for(ii2=0; (unsigned long)ii2<gr->size.x; ii2++){
						if( (cl=fgetc(fp)) == EOF ) return(READ_ERROR_RC);
						gr->data.rgb24[0][ii1][ii2].r = (unsigned char)cl;
					}
					for(ii2=0; (unsigned long)ii2<gr->size.x; ii2++){
						if( (cl=fgetc(fp)) == EOF) return(READ_ERROR_RC);
						gr->data.rgb24[0][ii1][ii2].g = (unsigned char)cl;
					}
					for(ii2=0; (unsigned long)ii2<gr->size.x; ii2++){
						if( (cl=fgetc(fp)) == EOF) return(READ_ERROR_RC);
						gr->data.rgb24[0][ii1][ii2].b = (unsigned char)cl;
					}
				}
			}else{
				RC_TRY( error_printf(7001) );
				return(READ_ERROR_RC);
			}
		}else{
			RC_TRY( error_printf(7001) );
			return(READ_ERROR_RC);
		}
		break;
	default: 
		return(UNKNOWN_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
read_tag (FILE *fp, DICOM_HEADER *dicom, unsigned long *ul)
{
	unsigned int ui1, ui2;
	int current_size;

	if(fgetc(fp) == EOF) return(END_RC);
	if(fseek(fp, -1, SEEK_CUR) != 0) return(SEEK_ERROR_RC);

	RC_TRY( lread_2bytes(fp, &ui1) );
	RC_TRY( lread_2bytes(fp, &ui2) );

	*ul = (unsigned long)ui1;
	*ul <<= 16;
	*ul += ui2;
	
	current_size = dicom->realloc_tag_size;
	if(dicom->tag_size == current_size){
		dicom->realloc_tag_size += DICOM_TAG_SIZE;
		dicom->tag = (unsigned long *)mm_realloc((void *)(dicom->tag),
		              sizeof(unsigned long) * (dicom->realloc_tag_size));
		if((dicom->tag) == (unsigned long *)NULL){
			RC_TRY( error_printf(9001,
			        sizeof(unsigned long)*(dicom->realloc_tag_size)) );
			return(ALLOC_ERROR_RC);
	    }
		current_size = dicom->realloc_tag_size;
	}

	dicom->tag[dicom->tag_size] = *ul;
	dicom->tag_size++;

	return(NORMAL_RC);
}


/* OB,OW,SQ,またはUNの明示的VRをもつデータ要素 */
static RC
skip_tag_elemOB (FILE *fp)
{
	unsigned long ul;
	unsigned int ui;
	unsigned int vr;
	
	RC_TRY( lread_2bytes(fp, &ui) );
	RC_TRY( lread_2bytes(fp, &vr) );
	RC_TRY( lread_4bytes(fp, &ul) );
	if (ul == 0) return(NORMAL_RC);
	if(fseek(fp, ul, SEEK_CUR) != 0) return(SEEK_ERROR_RC);
	return(NORMAL_RC);
}


/* OB,OW,SQ,またはUN以外の明示的VRをもつデータ要素 */
static RC
skip_tag_elemAE (FILE *fp, int vr_flag)
{
	unsigned int ui;
	unsigned int vr;
	unsigned long ul;

	if(vr_flag == 1){
		RC_TRY( lread_2bytes(fp, &vr) );
		RC_TRY( lread_2bytes(fp, &ui) );
		if(fseek(fp, ui, SEEK_CUR) != 0) return(SEEK_ERROR_RC);
	}else if(vr_flag == 0){
		RC_TRY( lread_4bytes(fp, &ul) );
		if(fseek(fp, ul, SEEK_CUR) != 0) return(SEEK_ERROR_RC);
	}else{
		return(UNKNOWN_ERROR_RC);
	}
	return(NORMAL_RC);
}


/* 暗黙的VRをもつデータ要素 */
static RC
skip_tag_elemNONE (FILE * fp)
{
	unsigned long ul;

	RC_TRY( lread_4bytes(fp, &ul) );
	if (ul == 0) return(NORMAL_RC);
	if(fseek(fp, ul, SEEK_CUR) != 0) return(SEEK_ERROR_RC);
	return(NORMAL_RC);
}


static RC
read_tag_elem (FILE *fp, char tmp[], unsigned long tag, int vr_flag)
{
	unsigned long ul;
	unsigned int  ui;
	unsigned int  vr;
	int index;

	if (vr_flag == 1){
		index = search_dicom_tag(dicom_elem_size, tag);
		if(index >= 0){
			if( (Dicom_elem[index].vr == OB) ||
				(Dicom_elem[index].vr == OW) ||
				(Dicom_elem[index].vr == SQ) ||
				(Dicom_elem[index].vr == UN)){
				RC_TRY( lread_2bytes(fp, &vr) );
				RC_TRY( lread_2bytes(fp, &ui) );
				RC_TRY( lread_4bytes(fp, &ul) );
			}else if(Dicom_elem[index].vr == NONE){
				RC_TRY( lread_2bytes(fp, &vr) );
				RC_TRY( lread_2bytes(fp, &ui) );
				ul = (unsigned long)ui;
			}else{
				RC_TRY( lread_2bytes(fp, &vr) );
				RC_TRY( lread_2bytes(fp, &ui) );
				ul = (unsigned long)ui;
			}
		}
	}else if(vr_flag == 0){
		RC_TRY( lread_4bytes(fp, &ul) );
	}else{
		return(UNKNOWN_ERROR_RC);
	}

	if(ul == 0){
		tmp[0] = (char)'\0';
		return(NORMAL_RC);
	}else if(ul >= DICOM_STR_SIZE){
		if( fread((void *)tmp, sizeof(char), DICOM_STR_SIZE-1, fp)
		    != DICOM_STR_SIZE-1){
			return(READ_ERROR_RC);
		}	
		tmp[DICOM_STR_SIZE-1] = (char)'\0';
		if(fseek(fp, (long)(ul - (DICOM_STR_SIZE-1)), SEEK_CUR) != 0){
			return(SEEK_ERROR_RC);
		}
	}else{
		if(fread((void *)tmp,sizeof(char),ul,fp) != ul) return(READ_ERROR_RC);
		tmp[ul] = (char)'\0';
	}

	return(NORMAL_RC);
}


static int
search_dicom_tag (const int size, const unsigned long target)
{
	int found_index = -1;
	int low = 0;
	int high = size;
	int middle;

	while(low <= high){
		middle = (low + high) / 2;
		if(target == Dicom_elem[middle].tag){
			found_index = middle;
			break;
		}else if(target < Dicom_elem[middle].tag){
			high = middle - 1;
		}else{
			low = middle + 1;
		}
	}
	return(found_index);
}


void
init_dicom (DICOM_HEADER *dicom)
{
	dicom->tag = NULL;
	dicom->tag_size = 0;
	dicom->realloc_tag_size = 0;
	dicom->length_to_end = 0;
	dicom->study_date[0] = (char)'\0';
	dicom->study_time[0] = (char)'\0';
	dicom->modality[0] = (char)'\0';
	dicom->manufacturer[0] = (char)'\0';
	dicom->station_name[0] = (char)'\0';
	dicom->patients_sex[0] = (char)'\0';
	dicom->patients_age[0] = (char)'\0';
	dicom->slice_thickness = 0.0; 
	dicom->gantry_detector = 0.0;
	dicom->table_height = 0.0;
	dicom->rotation_direction[0] = (char)'\0';
	dicom->pation_position[0] = (char)'\0';
	dicom->patient_orientation[0] = (char)'\0';
	dicom->slice_location = 0.0;
	dicom->samples_per_pixel = 0;
	dicom->photometric_interpretation[0] = (char)'\0';
	dicom->planar_configuration = 0;
	dicom->rows = 0;
	dicom->columns = 0;
	dicom->pixcel_spacing_x = 0.0;
	dicom->pixcel_spacing_y = 0.0;
	dicom->bits_allocated = 0;
	dicom->bits_stored = 0;
	dicom->pixel_representation = 0;
	dicom->smallest_image_pixel_value = 0;
	dicom->largest_image_pixel_value = 0;
	dicom->pixel_data = 0;
}

void
print_dicom_header (FILE *fp, DICOM_HEADER dicom)
{
/*
	fprintf(fp, "Length to End = %ld\n", dicom.length_to_end);
	fprintf(fp, "検査日付 = %s\n", dicom.study_date );
	fprintf(fp, "検査時間 = %s\n", dicom.study_time );
	fprintf(fp, "モダリティ = %s\n", dicom.modality);
	fprintf(fp, "製造者名 = %s\n", dicom.manufacturer);
	fprintf(fp, "ステーション名 = %s\n", dicom.station_name);
	fprintf(fp, "患者性別 = %s\n", dicom.patients_sex);
	fprintf(fp, "患者年令 = %s\n", dicom.patients_age);
	fprintf(fp, "スライス厚さ = %8.3f [mm]\n", dicom.slice_thickness);
	fprintf(fp, "架台/検出器傾き = %8.3f [deg]\n", dicom.gantry_detector);
	fprintf(fp, "テーブル高さ = %8.3f [mm]\n", dicom.table_height);
	fprintf(fp, "回転方向 = %s: CW:時計回り  CC:反時計回り＼n",
	             dicom.rotation_direction);
	fprintf(fp, "患者位置 = %s: HFP:Head First-Prone(うつむき)\n"
	            "                 HFS:Head First-Supinei(仰向け)\n"
	            "                 HFDR:Head First-Decubitus Right(臥床姿勢)\n"
	            "                 HFDL:Head First-Decubitus Left\n"
	            "                 FFDL:Feet First_Decubitus Left\n"
	            "                 FFDL:Feet First_Decubitus Right\n"
	            "                 FFP:Feet First_Prone\n"
	            "                 FFS:Feet First_Supine\n",
	             dicom.pation_position);
	fprintf(fp, "患者方向 = %s: A(anterior:前面), P(posterior:背面)\n"
	            "                 R(right:右), L(left:左)\n"
	            "                 H(head:頭), F(foot:足)\n",
	             dicom.patient_orientation);
	fprintf(fp, "スライス位置 = %8.3f\n", dicom.slice_location);
	fprintf(fp, "画素あたりのサンプル = %ld\n", dicom.samples_per_pixel);
	fprintf(fp, "光度測定解釈 = %s\n", dicom.photometric_interpretation);
	fprintf(fp, "面構成 = %d\n", dicom.planar_configuration);
	fprintf(fp, "横行X size = %ld\n", dicom.rows);
	fprintf(fp, "縦行Y size = %ld\n", dicom.columns);
	fprintf(fp, "画素間隔X = %6.4g [mm]\n", dicom.pixcel_spacing_x);
	fprintf(fp, "画素間隔Y = %6.4g [mm]\n", dicom.pixcel_spacing_y);
	fprintf(fp, "割り当てビット = %ld [bit]\n", dicom.bits_allocated);
	fprintf(fp, "格納ビット = %ld [bit]\n", dicom.bits_stored);
	fprintf(fp, "画素表現 = %ld\n", dicom.pixel_representation);
	fprintf(fp, "最小画素数 = %ld\n", dicom.smallest_image_pixel_value);
	fprintf(fp, "最大画素数 = %ld\n", dicom.largest_image_pixel_value);
	fprintf(fp, "pixel size[X*Y*ColorBit] = %ld\n", dicom.pixel_data);
*/
	fprintf(fp, "Length to End = %ld\n", dicom.length_to_end);
	fprintf(fp, "Study Date = %s\n", dicom.study_date );
	fprintf(fp, "Study Time = %s\n", dicom.study_time );
	fprintf(fp, "Modality = %s\n", dicom.modality);
	fprintf(fp, "Manufacturer = %s\n", dicom.manufacturer);
	fprintf(fp, "Station Name = %s\n", dicom.station_name);
	fprintf(fp, "Patients Sex = %s\n", dicom.patients_sex);
	fprintf(fp, "Patients Age = %s\n", dicom.patients_age);
	fprintf(fp, "Slice Thickness = %8.3f [mm]\n", dicom.slice_thickness);
	fprintf(fp, "Gantry/Detector Tilt = %8.3f [deg]\n", dicom.gantry_detector);
	fprintf(fp, "Table Height = %8.3f [mm]\n", dicom.table_height);
	fprintf(fp, "Rotation Direction = %s: CW  CC\n",
	             dicom.rotation_direction);
	fprintf(fp, "Patient Position = %s: HFP:Head First-Prone\n"
	            "                 HFS:Head First-Supinei\n"
	            "                 HFDR:Head First-Decubitus Right\n"
	            "                 HFDL:Head First-Decubitus Left\n"
	            "                 FFDL:Feet First_Decubitus Left\n"
	            "                 FFDL:Feet First_Decubitus Right\n"
	            "                 FFP:Feet First_Prone\n"
	            "                 FFS:Feet First_Supine\n",
	             dicom.pation_position);
	fprintf(fp, "Patient Orientation = %s: A(anterior), P(posterior)\n"
	            "                 R(right), L(left)\n"
	            "                 H(head), F(foot)\n",
	             dicom.patient_orientation);
	fprintf(fp, "Slice Location = %8.3f\n", dicom.slice_location);
	fprintf(fp, "Samples Per Pixel = %ld\n", dicom.samples_per_pixel);
	fprintf(fp, "Photometric Interpretation = %s\n",
	             dicom.photometric_interpretation);
	fprintf(fp, "Planar Configuration = %d\n", dicom.planar_configuration);
	fprintf(fp, "RowsX Size = %ld\n", dicom.rows);
	fprintf(fp, "ColumnsY Size = %ld\n", dicom.columns);
	fprintf(fp, "Pixel SpacingX = %6.4g [mm]\n", dicom.pixcel_spacing_x);
	fprintf(fp, "Pixel SpacingY = %6.4g [mm]\n", dicom.pixcel_spacing_y);
	fprintf(fp, "Bits Allocated = %ld [bit]\n", dicom.bits_allocated);
	fprintf(fp, "Bits Stored = %ld [bit]\n", dicom.bits_stored);
	fprintf(fp, "Pixel Representation = %ld\n", dicom.pixel_representation);
	fprintf(fp, "Smallest Image Pixel Value = %ld\n",
	             dicom.smallest_image_pixel_value);
	fprintf(fp, "Largest Image Pixel Value = %ld\n",
	             dicom.largest_image_pixel_value);
	fprintf(fp, "pixel size[X*Y*ColorBit] = %ld\n", dicom.pixel_data);
}


RC
print_dicom_tag (FILE *fp, DICOM_HEADER dicom)
{
	int ii1, index;
	unsigned long current_num;
	int dicom_elem_size = (sizeof(Dicom_elem)/sizeof(DICOM_ELEM));

	if(dicom.tag == NULL) return(UNKNOWN_ERROR_RC);
	
	current_num = dicom.tag[0]>>16;
	for(ii1=0;ii1<dicom.tag_size;ii1++){
		index = search_dicom_tag(dicom_elem_size, dicom.tag[ii1]);
		if(index >= 0){
			if(current_num != (dicom.tag[ii1]>>16)){
				fprintf(fp,"\n");
				current_num = dicom.tag[ii1]>>16;
			}
			fprintf(fp, "%0#10lx [%s]\n",dicom.tag[ii1],Dicom_elem[index].name);
		}else{
			if((dicom.tag[ii1]>>16 % 2) != 0){
				if(current_num != (dicom.tag[ii1]>>16)){
					fprintf(fp,"\n");
					current_num = dicom.tag[ii1]>>16;
				}
				/*
				fprintf(fp, "%0#10lx [私的要素]\n", dicom.tag[ii1]);
				*/
				fprintf(fp, "%0#10lx [Pricate Element]\n", dicom.tag[ii1]);
			}else{
				RC_TRY( error_printf(7301, dicom.tag[ii1]) );
				return(IMPLEMENT_ERROR_RC);
			}
		}
	}
	return(NORMAL_RC);
}


static RC
check_vr (FILE *fp, int *vr_flag)
{
	int ii1;
	char tmp[3];
	char *vr[] = {"OB","OW","SQ","UN","AE","AS",
	              "AT","CS","DA","DS","DT","FL",
	              "IS","LO","LT","PN","SH","SL",
	              "SS","ST","TM","UI","UL","US","UT"};

	if((tmp[0] = (char)fgetc(fp)) == EOF) return(READ_ERROR_RC);
	if((tmp[1] = (char)fgetc(fp)) == EOF) return(READ_ERROR_RC);

	tmp[2] = (char)'\0';
	if(fseek(fp, -2L, SEEK_CUR) != 0) return(SEEK_ERROR_RC);

	for(ii1=0;ii1<25;ii1++){
		if(strncmp(vr[ii1], tmp, 2) == 0){
			*vr_flag = 1;
			return(NORMAL_RC);
		}
	}
	*vr_flag = 0;
	return(NORMAL_RC);
}


/*
 * Pixel_Representation(画素表現)がセットされていると，そのモノクロ
 * 画像は 2の補数表現をしている．それをazlibPIX向けのmono16に変換する
 * Pixel_Representationがセットされていても，2の補数でないピクセルが
 * 含まれていることがあるので注意．各ピクセルごとに調べるしかない．
 */ 
RC
conv_2complement_to_mono16 (GRAPHIC src, GRAPHIC *dest, DICOM_HEADER dicom)
{
	RC_NULL_CHK( dest );
	RC_TRY( graphic_exist_chk(src) );
	RC_TRY( copy_graphic_info(src, dest) );
	RC_TRY( allocate_graphic_data(dest) );

	GRAPHIC_SIZE_LOOP_3D(src.size.z, src.size.y, src.size.x)
	{
		/* src 2の補数表現かどうか2重のチェック */
		if( (dicom.pixel_representation != 0) &&
			((src.data.mono16[ii1][ii2][ii3] & (1 << 15)) != 0) ){
			dest->data.mono16[ii1][ii2][ii3] =
			  -( (~(src.data.mono16[ii1][ii2][ii3]-1))&(0xFFFF) );
		}else
			dest->data.mono16[ii1][ii2][ii3] = src.data.mono16[ii1][ii2][ii3];
	}GRAPHIC_SIZE_LOOP_3D_END;
	return(NORMAL_RC);
}


