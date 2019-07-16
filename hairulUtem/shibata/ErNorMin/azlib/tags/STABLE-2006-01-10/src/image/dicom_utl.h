/*********************************************************************
 * dicom_utl.h
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

/* $Id: dicom_utl.h 360 2004-07-27 08:09:03Z sasaoka $ */

#ifndef DICOM_UTL_H
#define DICOM_UTL_H

#include <stdio.h>
#include "rc.h"
#include "graphic_utl.h"

#define DICOM_STR_SIZE (32)
#define DICOM_TAG_SIZE (20)

typedef struct {
	unsigned long *tag;
	int tag_size;
	int realloc_tag_size;
	unsigned long length_to_end;
	char study_date[DICOM_STR_SIZE];
	char study_time[DICOM_STR_SIZE];
	char modality[DICOM_STR_SIZE];
	char manufacturer[DICOM_STR_SIZE];
	char station_name[DICOM_STR_SIZE];
	char patients_sex[DICOM_STR_SIZE];
	char patients_age[DICOM_STR_SIZE];
	double slice_thickness;
	double gantry_detector;
	double table_height;
	char rotation_direction[DICOM_STR_SIZE];
	char pation_position[DICOM_STR_SIZE];
	char patient_orientation[DICOM_STR_SIZE];
	double slice_location;
	unsigned long samples_per_pixel;
	char photometric_interpretation[DICOM_STR_SIZE];
	unsigned int planar_configuration;
	unsigned long rows;
	unsigned long columns;
	double pixcel_spacing_x;
	double pixcel_spacing_y;
	unsigned long bits_allocated;
	unsigned long bits_stored;
	unsigned long pixel_representation;
	unsigned long smallest_image_pixel_value;
	unsigned long largest_image_pixel_value;
	unsigned long pixel_data;
} DICOM_HEADER;

RC read_dicom(FILE *fp, DICOM_HEADER *dicom, GRAPHIC *gr);
RC read_dicom_header(FILE *fp, DICOM_HEADER *dicom, int *pimage);
RC read_dicom_image(FILE *fp, DICOM_HEADER dicom, GRAPHIC *gr, int pimage);
void init_dicom(DICOM_HEADER *dicom);
void print_dicom_header(FILE *fp, DICOM_HEADER dicom);
RC print_dicom_tag(FILE *fp, DICOM_HEADER dicom);
RC set_dicom_header(DICOM_HEADER dicom, GRAPHIC *gr);
RC conv_2complement_to_mono16(GRAPHIC src, GRAPHIC *dst, DICOM_HEADER dicom);

#endif /*  DICOM_UTL_H  */
 
