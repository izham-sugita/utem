/*********************************************************************
 * nst_pch_component.c
 *
 * Copyright (C) 2001 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: nst_pch_component.c 842 2006-12-04 06:21:20Z sasaoka $ */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "nst_component.h"

#define BUFFER_SIZE     (512)
#define INIT_SIZE      (1024)

static RC seek_pch(FILE *fp, int *pch_data, int *subcase, ELEM_TYPE *elem_type);
static void idc2ss_3d(const IDC idc1[], const IDC idc2[], const IDC idc3[],
                      const IDC idc4[], const IDC idc5[], const IDC idc6[],
                      const IDC idc7[], int element_id,
                      int node_id, STRESS_STRAIN *ss);
static RC read_pch_ss(FILE *fp, int subcase, int pch_sw, int pch_data,
                      STRAIN_ARRAY *ss);


RC
read_pch_disp (FILE *fp, int subcase, int pch_sw, DISP_ARRAY *disp)
{
	char buf[BUFFER_SIZE];
	int tmp_pch_data, tmp_subcase;
	ELEM_TYPE tmp_elem_type;
	IDC idc1[6];
	IDC idc2[5];
	RC rc;
	int eof_flag;
	long file_position;

	RC_TRY( allocate_disp_array(INIT_SIZE, disp) );
	disp->source = FEM_NASTRAN;

	rewind(fp);
	file_position = ftell(fp);
	eof_flag = 0;
	while(eof_flag == 0){
		if(fseek(fp, file_position, SEEK_SET) != 0) return(SEEK_ERROR_RC);
		while(1){
			rc = seek_pch(fp, &tmp_pch_data, &tmp_subcase, &tmp_elem_type);
			if(rc == END_RC){
				eof_flag = 1;
				break;
			}else if(rc != NORMAL_RC){
				return(rc);
			}
			if( (tmp_pch_data == PCH_DATA_DISP)
			  &&( (!(pch_sw & PCH_SW_SUBCASE))||(tmp_subcase == subcase) ) ){
				break;
			}
		}
		if(eof_flag) break;

		while(1){
			file_position = ftell(fp);
			if(fgets(buf, sizeof(buf), fp) == NULL){
				if(feof(fp)){
					eof_flag = 1;
					break;
				}
				return(READ_ERROR_RC);
			}
			if(buf[0] == (char)'$') break;
			if(sgetidc(buf, "iadddi", idc1) != NORMAL_RC) break;
			RC_TRY( fgetidc(fp, "sdddi", idc2) );

			RC_TRY( realloc_disp_array(disp) );
			disp->array[disp->size].node = idc1[0].i;
			disp->array[disp->size].v.x = idc1[2].d;
			disp->array[disp->size].v.y = idc1[3].d;
			disp->array[disp->size].v.z = idc1[4].d;
			disp->array[disp->size].v.xy = idc2[1].d;
			disp->array[disp->size].v.yz = idc2[2].d;
			disp->array[disp->size].v.zx = idc2[3].d;
			disp->array[disp->size].i_info[0] = tmp_subcase;
			(disp->size)++;
		}
	}

	if(disp->size == 0) return(END_RC);

	RC_TRY( clean_disp_array(disp) );
	return(NORMAL_RC);
}


RC
read_pch_stress (FILE *fp, int subcase, int pch_sw, STRESS_ARRAY *stress)
{
	STRAIN_ARRAY false_strain;

	false_strain = stress2strain_array(*stress);
	RC_TRY( read_pch_ss(fp, subcase, pch_sw, PCH_DATA_STRESS, &false_strain) );
	*stress = strain2stress_array(false_strain);

	return(NORMAL_RC);
}


RC
read_pch_strain (FILE *fp, int subcase, int pch_sw, STRAIN_ARRAY *strain)
{
	RC_TRY( read_pch_ss(fp, subcase, pch_sw, PCH_DATA_STRAIN, strain) );

	return(NORMAL_RC);
}


static RC
read_pch_ss (FILE *fp, int subcase, int pch_sw, int pch_data, STRAIN_ARRAY *ss)
{
	char buf[BUFFER_SIZE];
	int tmp_pch_data, tmp_subcase;
	ELEM_TYPE tmp_elem_type;
	IDC idc0[5], idc1[5], idc2[5], idc3[5];
	IDC idc4[5], idc5[5], idc6[5], idc7[5];
	RC rc;
	int eof_flag;
	int ii1;
	int elem_id = -1;
	long file_position;
	int node_num;

	RC_TRY( allocate_strain_array(INIT_SIZE, ss) );
	ss->source = FEM_NASTRAN;

	rewind(fp);
	file_position = ftell(fp);
	eof_flag = 0;
	while(eof_flag == 0){
		if(fseek(fp, file_position, SEEK_SET) != 0) return(SEEK_ERROR_RC);
		while(1){
			rc = seek_pch(fp, &tmp_pch_data, &tmp_subcase, &tmp_elem_type);
			if(rc == END_RC){
				eof_flag = 1;
				break;
			}else if(rc != NORMAL_RC){
				return(rc);
			}
			if( (tmp_pch_data == pch_data)
			  &&( (!(pch_sw & PCH_SW_SUBCASE))||(tmp_subcase == subcase) ) ){
				break;
			}
		}
		if(eof_flag) break;

		if( (tmp_elem_type == ELEM_TETRA1)
		  ||(tmp_elem_type == ELEM_PENTA1)
		  ||(tmp_elem_type == ELEM_HEXA1) ){
			while(1){
				file_position = ftell(fp);
				if(fgets(buf, sizeof(buf), fp) == NULL){
					if(feof(fp)){
						eof_flag = 1;
						break;
					}
					return(READ_ERROR_RC);
				}
				if(buf[0] == (char)'$'){
					break;
				}
				if(sgetidc(buf, "idddi", idc0) != NORMAL_RC) break;
				elem_id = idc0[0].i;

				if(tmp_elem_type == ELEM_TETRA1){
					node_num = 4;
				}else if(tmp_elem_type == ELEM_PENTA1){
					node_num = 6;
				}else{   /* ELEM_HEXA1 */
					node_num = 8;
				}

				for(ii1=-1; ii1<node_num; ii1++){
					RC_TRY( fgetidc(fp,  "sdddi", idc1) );
					RC_TRY( fgetidc(fp,  "sdddi", idc2) );
					RC_TRY( fgetidc(fp,  "sdddi", idc3) );
					RC_TRY( fgetidc(fp,  "sdddi", idc4) );
					RC_TRY( fgetidc(fp,  "sdddi", idc5) );
					RC_TRY( fgetidc(fp,  "sdddi", idc6) );
					RC_TRY( fgetidc(fp,  "sdddi", idc7) );

					RC_TRY( realloc_strain_array(ss) );
					idc2ss_3d(idc1, idc2, idc3, idc4, idc5, idc6, idc7,
					          elem_id, ii1, &(ss->array[ss->size]));
					(ss->size)++;
				}
			}
		}else{
			return(IMPLEMENT_ERROR_RC);
		}
	}

	if(ss->size == 0) return(END_RC);

	RC_TRY( clean_strain_array(ss) );
	return(NORMAL_RC);
}


static void
idc2ss_3d (const IDC idc1[], const IDC idc2[], const IDC idc3[],
           const IDC idc4[], const IDC idc5[], const IDC idc6[],
           const IDC idc7[], int element_id, int node_id, STRESS_STRAIN *ss)
{
	ss->element = element_id;
	ss->node = node_id;
	ss->v.x = idc1[2].d;
	ss->v.y = idc4[1].d;
	ss->v.z = idc6[1].d;
	ss->v.xy = idc1[3].d;
	ss->v.yz = idc4[2].d;
	ss->v.zx = idc6[2].d;

	ss->d_info[0] = idc3[2].d;     /* MEAN PRESSURE */
	ss->i_info[0] = 2;
}


static RC
seek_pch (FILE *fp, int *pch_data, int *subcase, ELEM_TYPE *elem_type)
{
	char buf[BUFFER_SIZE];
	long file_position;
	IDC idc;

	/* seek "$********************" line */
	while(1){
		file_position = ftell(fp);

		if(fgets(buf, sizeof(buf), fp) == NULL){
			if(feof(fp)) return(END_RC);
			return(READ_ERROR_RC);
		}
		if(buf[0] == (char)'$') break;
	}
	if(fseek(fp, file_position, SEEK_SET) != 0) return(SEEK_ERROR_RC);

	*pch_data = PCH_DATA_UNKNOWN;
	*subcase = -1;
	*elem_type = ELEM_RBAR;
	while(1){
		file_position = ftell(fp);

		if(fgets(buf, sizeof(buf), fp) == NULL){
			if(feof(fp)) return(END_RC);
			return(READ_ERROR_RC);
		}
		if(buf[0] != (char)'$') break;

		      if(strncmp(&(buf[1]), "DISPLACEMENTS", 13) == 0){
			*pch_data = PCH_DATA_DISP;
		}else if(strncmp(&(buf[1]), "ELEMENT STRESSES", 16) == 0){
			*pch_data = PCH_DATA_STRESS;
		}else if(strncmp(&(buf[1]), "ELEMENT STRAINS", 15) == 0){
			*pch_data = PCH_DATA_STRAIN;
		}else if(strncmp(&(buf[1]), "ELEMENT TYPE", 12) == 0){
			RC_TRY( sgetidc(&(buf[13]), "i", &idc) );
			if(idc.i == 39){
				*elem_type = ELEM_TETRA1;
			}else if(idc.i == 68){
				*elem_type = ELEM_PENTA1;
			}else if(idc.i == 67){
				*elem_type = ELEM_HEXA1;
			}else{
				return(IMPLEMENT_ERROR_RC);
			}
		}else if(strncmp(&(buf[1]), "SUBCASE ID", 10) == 0){
			RC_TRY( sgetidc(&(buf[11]), "i", &idc) );
			*subcase = idc.i;
		}
	}
	if(fseek(fp, file_position, SEEK_SET) != 0){
		return(SEEK_ERROR_RC);
	}

	return(NORMAL_RC);
}


