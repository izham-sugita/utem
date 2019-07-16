/*********************************************************************
 * nst_f06_component.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Thanks to following contributors.
 *  <Yasuo ASAGA> <Tomoyuki TSUBATA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: nst_f06_component.c,v 1.8 2003/07/23 04:06:31 sasaoka Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rc.h"
#include "string_utl.h"
#include "nst_component.h"

#define BUFFER_SIZE  (512)
#define INIT_SIZE    (1024)

typedef enum {
	F06_3D,
	F06_2D_1O,
	F06_2D_2O,
	F06_CBAR
} F06STAT;


static RC seek_f06(FILE *fp, const char *key, char *buf, int buf_size,
                   int *subcase, double *load_step, int *mode);
static RC check_f06key(const char *buf, const char *key);
static int check_f06subcase_load(int f06_sw, int sub1, int sub2,
                                 double load1, double load2,
                                 int mode1, int mode2);
static RC approach_f06eigen_value(FILE *fp);
static RC approach_f06disp(FILE *fp);
static RC check_f06eigen_value(char *buf);
static RC check_f06disp(char *buf);
static void idc2strain_2d(const IDC *idc, int elem_id, int node_id,
                          FEM_STRESS_STRAIN *strain);
static void idc2strain_3d(const IDC *idc, const IDC *idc2, const IDC *idc3,
                          int elem_id, int node_id, FEM_STRESS_STRAIN *strain);
static RC approach_f06strain(FILE *fp, FEM_ELEM_TYPE elem_type);
static RC check_f06_2D(char *buf);
static RC check_f06_3D(char *buf);
static RC skip_file(FILE *fp, int record);
static RC read_f06ss(FILE *fp, int subcase, double load_step, int mode,
                     int f06_sw, FEM_STRAIN_ARRAY *strain,
                     const char *key, int f06_ss);



RC nst_input_disp(FILE *fp, FEM_DISP_ARRAY *disp)
{
	RC_TRY( read_f06disp(fp, 1, 1.0, -1, 0, disp, F06_KEY_DISP) );

	return(NORMAL_RC);
}


RC nst_input_eigenvector_value(FILE *fp, FEM_DISP_ARRAY *evect, double value)
{
	RC_TRY( read_f06disp(fp, 1, value, -1, F06_SW_LOADSTEP,
	                     evect, F06_KEY_EIGENVECTOR) );

	return(NORMAL_RC);
}


RC nst_input_eigenvector_mode(FILE *fp, FEM_DISP_ARRAY *evect, int mode)
{
	RC_TRY( read_f06disp(fp, 1, 1.0, mode, F06_SW_MODE,
	                     evect, F06_KEY_EIGENVECTOR) );

	return(NORMAL_RC);
}


RC nst_input_eigen_value(FILE *fp, int *array_size, double **value_array)
{
	RC_TRY( read_f06eigen_value(fp, 1, F06_SW_NULL,
	                            array_size, value_array) );
	return(NORMAL_RC);
}


RC read_f06stress(FILE *fp, int subcase, double load_step, int mode,
                  int f06_sw, FEM_STRESS_ARRAY *stress, int f06_ss)
{
	FEM_STRAIN_ARRAY false_strain;

	false_strain = fem_stress2strain_array(*stress);
	RC_TRY( read_f06ss(fp, subcase, load_step, mode, f06_sw,
	                   &false_strain, F06_KEY_STRESS, f06_ss) );
	*stress = fem_strain2stress_array(false_strain);

	return(NORMAL_RC);
}


RC read_f06strain(FILE *fp, int subcase, double load_step, int mode,
                  int f06_sw, FEM_STRAIN_ARRAY *strain, int f06_ss)
{
	RC_TRY( read_f06ss(fp, subcase, load_step, mode, f06_sw,
	                   strain, F06_KEY_STRAIN, f06_ss) );

	return(NORMAL_RC);
}


RC read_f06eigen_value(FILE *fp, int subcase, int f06_sw, int *array_size,
                       double **value_array)
{
	int subc0;
	double loads0;
	int mode0;
	IDC idc[7];
	RC rc;
	double *tmp;
	char buf[BUFFER_SIZE];

	rewind(fp);

	while(1){
		RC_TRY( seek_f06(fp, F06_KEY_EIGENVALUE, buf, sizeof(buf),
		                 &subc0, &loads0, &mode0) );
		if(f06_sw & F06_SW_SUBCASE){
			if(subc0 != subcase) continue;
		}
		break;
	}
	RC_TRY( approach_f06eigen_value(fp) );

	*value_array = NULL;
	*array_size = 0;
	while(1){
		if(fgets(buf, sizeof(buf), fp) == NULL){
			if(feof(fp)) break;
			fprintf(stderr, "fgets < ");
			return(READ_ERROR_RC);
		}
		if((int)buf[0] == '1') break;

		rc = sgetidc(buf, "iiddddd", idc);
		if(rc != NORMAL_RC)break;

		tmp = (double *)realloc(*value_array,
		                        ((*array_size)+1)*sizeof(double));
		if(tmp == NULL){
			fprintf(stderr, "realloc < ");
			return(ALLOC_ERROR_RC);
		}
		*value_array = tmp;

		(*value_array)[*array_size] = idc[2].d;
		(*array_size)++;
	}

	return(NORMAL_RC);
}


RC read_f06disp(FILE *fp, int subcase, double load_step, int mode, int f06_sw,
                FEM_DISP_ARRAY *disp, const char *key)
{
	int subc0;
	double loads0;
	int mode0;
	IDC idc[8];
	char buf[BUFFER_SIZE];
	long file_position;
	int end_flag = 0;
	RC rc;

	rewind(fp);

	RC_TRY(allocate_fem_disp_array(INIT_SIZE, disp) );
	disp->source = FEM_NASTRAN;
	
	while(1){
		/* seek target displacement */
		while(1){
			rc = seek_f06(fp, key, buf, sizeof(buf), &subc0, &loads0, &mode0);
			if(rc == END_RC){
				end_flag = 1;
				break;
			}else if(rc == NORMAL_RC){
				if(check_f06subcase_load(f06_sw, subc0, subcase,
		                     loads0, load_step, mode0, mode))break;
			}else{
				fprintf(stderr,"seek_f06 < ");
				return(rc);
			}
		}
		if(end_flag) break;
		RC_TRY( approach_f06disp(fp) );

		/* read target displacement */
		while(1){
			/* store current file position */
			file_position = ftell(fp);

			if(fgets(buf, sizeof(buf), fp) == NULL){
				fprintf(stderr,"fgets[2] < ");
				return(READ_ERROR_RC);
			}

			if((buf[0] == '0')||(buf[0] == '1')){
				/* restore position */
				if(fseek(fp, file_position, SEEK_SET) != 0){
					fprintf(stderr,"fseek < ");
					return(SEEK_ERROR_RC);
				}
				break;
			}

			RC_TRY( realloc_fem_disp_array(disp) );
			RC_TRY( sgetidc(buf, "iadddddd", idc) );
			disp->array[disp->size].node = idc[0].i;
			disp->array[disp->size].v.x = idc[2].d;
			disp->array[disp->size].v.y = idc[3].d;
			disp->array[disp->size].v.z = idc[4].d;
			disp->array[disp->size].v.yz = idc[5].d;
			disp->array[disp->size].v.zx = idc[6].d;
			disp->array[disp->size].v.xy = idc[7].d;
			disp->array[disp->size].i_info[0] = subc0;
			disp->array[disp->size].i_info[1] = mode0;
			disp->array[disp->size].d_info[0] = loads0;
			(disp->size)++;
		}
	}

	RC_TRY( clean_fem_disp_array(disp) );

	if(disp->size == 0) return(ZERO_ERROR_RC);
	return(NORMAL_RC);
}


/* ページ先頭（1から始まる行）から５行目に文字列 key を含む行を探し、 */
/* 該当行を buf[buf_size] に格納する */
static RC seek_f06(FILE *fp, const char *key, char *buf, int buf_size,
                   int *subcase, double *load_step, int *mode)
{
	char buf_1[BUFFER_SIZE];
	char buf_dum[BUFFER_SIZE];
	char buf_sub[BUFFER_SIZE];
	char buf_step[BUFFER_SIZE];
	IDC idc[2];
	long file_position;

	while(1){
		/* seek "1***************" line */
		while(1){
			if(fgets(buf_1, sizeof(buf_1), fp) == NULL){
				if(feof(fp))return(END_RC);
				fprintf(stderr,"fgets[0] < ");
				return(READ_ERROR_RC);
			}
			if(buf_1[0] == '1')break;
		}

		/* store current positon */
		file_position = ftell(fp);

		if(fgets(buf_dum, sizeof(buf_dum), fp) == NULL){
			if(feof(fp))return(END_RC);
			fprintf(stderr,"fgets[1] < ");
			return(READ_ERROR_RC);
		}

		if(fgets(buf_sub, sizeof(buf_sub), fp) == NULL){
			if(feof(fp))return(END_RC);
			fprintf(stderr,"fgets[2] < ");
			return(READ_ERROR_RC);
		}
		if(check_f06key(buf_sub, "SUBCASE") != NORMAL_RC){
			/* restore position */
			if(fseek(fp, file_position, SEEK_SET) != 0){
				fprintf(stderr,"fseek < ");
				return(SEEK_ERROR_RC);
			}
			continue;
		}
		RC_TRY( sgetidc(buf_sub, "ii", idc) );
		(*subcase) = idc[1].i;

		if(fgets(buf_step, sizeof(buf_step), fp) == NULL){
			if(feof(fp))return(END_RC);
			fprintf(stderr,"fgets[3] < ");
			return(READ_ERROR_RC);
		}
		if( sgetidc(buf_step, "d", idc) == NORMAL_RC ){
			(*load_step) = idc[0].d;
		}else{
			(*load_step) = 1.0;
		}

		if(fgets(buf, buf_size, fp) == NULL){
			if(feof(fp))return(END_RC);
			fprintf(stderr,"fgets[4] < ");
			return(READ_ERROR_RC);
		}
		if(check_f06key(buf, key) != NORMAL_RC){
			/* restore position */
			if(fseek(fp, file_position, SEEK_SET) != 0){
				fprintf(stderr,"fseek < ");
				return(SEEK_ERROR_RC);
			}
			continue;
		}
		if(check_f06key(buf, "CYCLES") == NORMAL_RC){
			if( sgetidc(buf, "di", idc) == NORMAL_RC){
				(*mode) = idc[1].i;
			}else{
				(*mode) = -1;
			}
		}else{
			if( sgetidc(buf, "i", idc) == NORMAL_RC){
				(*mode) = idc[0].i;
			}else{
				(*mode) = -1;
			}
		}

		break;
	}

	return(NORMAL_RC);
}


static RC check_f06key(const char *buf, const char *key)
{
	int len;
	int ii = 0;

	len = strlen(key);
	while(buf[ii] != '\0'){
		if(buf[ii] == key[0]){   /* for reduction calculating time */
			if(strncmp(&(buf[ii]), key, len) == 0){
				return(NORMAL_RC);
			}
		}
		ii++;
	}

	return(END_RC);
}


/* check subcase, loadset, mode */
static int check_f06subcase_load(int f06_sw, int sub1, int sub2,
                                 double load1, double load2,
                                 int mode1, int mode2)
{
	if(f06_sw & F06_SW_SUBCASE){
		if(sub1 != sub2)return(0);
	}

	if(f06_sw & F06_SW_LOADSTEP){
		if( nearly_eq(load1, load2) == 0){
			return(0);
		}
	}

	if(f06_sw & F06_SW_MODE){
		if(mode1 != mode2)return(0);
	}

	return(1);
}


static RC approach_f06eigen_value(FILE *fp)
{
	char buf[BUFFER_SIZE];

	while(1){
		if(fgets(buf, sizeof(buf), fp) == NULL){
			if(feof(fp)){
				return(END_RC);
			}else{
				fprintf(stderr,"fgets[1] < ");
				return(READ_ERROR_RC);
			}
		}
		if(check_f06eigen_value(buf) == NORMAL_RC) break;
	}

	if(fgets(buf, sizeof(buf), fp) == NULL){
		if(feof(fp)){
			return(END_RC);
		}else{
			fprintf(stderr,"fgets[1] < ");
			return(READ_ERROR_RC);
		}
	}

	return(NORMAL_RC);
}


static RC approach_f06disp(FILE *fp)
{
	char buf[BUFFER_SIZE];

	while(1){
		if(fgets(buf, sizeof(buf), fp) == NULL){
			if(feof(fp)){
				return(END_RC);
			}else{
				fprintf(stderr,"fgets[1] < ");
				return(READ_ERROR_RC);
			}
		}
		if(check_f06disp(buf) == NORMAL_RC) break;
	}

	return(NORMAL_RC);
}


static RC check_f06eigen_value(char *buf)
{
	if(check_f06key(buf,"MODE") == NORMAL_RC){
		if((check_f06key(buf,"EXTRACTION") == NORMAL_RC)
				 &&(check_f06key(buf,"EIGENVALUE") == NORMAL_RC)
				 &&(check_f06key(buf,"RADIANS") == NORMAL_RC)
				 &&(check_f06key(buf,"CYCLES") == NORMAL_RC)
				 &&(check_f06key(buf,"GENERALIZED") == NORMAL_RC) ){
			return(NORMAL_RC);
		}else{
			return(SPECIAL_RC);
		}
	}else{
		return(SPECIAL_RC);
	}
}


static RC check_f06disp(char *buf)
{
	if(check_f06key(buf,"POINT ID.") == NORMAL_RC){
		if((check_f06key(buf,"TYPE") == NORMAL_RC)
				 &&(check_f06key(buf,"T1") == NORMAL_RC)
				 &&(check_f06key(buf,"T2") == NORMAL_RC)
				 &&(check_f06key(buf,"T3") == NORMAL_RC)
				 &&(check_f06key(buf,"R1") == NORMAL_RC)
				 &&(check_f06key(buf,"R2") == NORMAL_RC)
				 &&(check_f06key(buf,"R3") == NORMAL_RC)){
			return(NORMAL_RC);
		}else{
			return(SPECIAL_RC);
		}
	}else{
		return(SPECIAL_RC);
	}
}


/* 応力またはひずみの読み込み */
/* f06_ss : F06_SS_SIMPLE : 要素中心の情報のみ */
/*          F06_SS_CORNER : 要素の各頂点の情報も読み込む */
/*          F06_SS_FIBER  : 全 fiber distance, strain curvature を読み込む */
static RC read_f06ss(FILE *fp, int subcase, double load_step, int mode,
                     int f06_sw, FEM_STRAIN_ARRAY *strain,
                     const char *key, int f06_ss)
{
	int subc0;
	double loads0;
	int mode0;
	IDC idc[11];
	char buf[BUFFER_SIZE];
	char buf_key[BUFFER_SIZE];
	long file_position;
	int end_flag = 0;
	RC rc;
	FEM_ELEM_TYPE elem_type;
	int elem_id = -1;
	int node_id = -1;

	rewind(fp);

	RC_TRY( allocate_fem_strain_array(INIT_SIZE, strain) );
	strain->source = FEM_NASTRAN;
	
	while(1){
		/* seek target strain */
		while(1){
			rc = seek_f06(fp, key, buf_key, sizeof(buf_key),
			              &subc0, &loads0, &mode0);
			if(rc == END_RC){
				end_flag = 1;
				break;
			}else if(rc == NORMAL_RC){
				if(check_f06subcase_load(f06_sw, subc0, subcase,
		                     loads0, load_step, mode0, mode))break;
			}else{
				fprintf(stderr,"seek_f06 < ");
				return(rc);
			}
		}
		if(end_flag) break;

		if(check_f06key(buf_key, "T E T R A") == NORMAL_RC){
			elem_type = ELEM_TETRA1;
		}else if(check_f06key(buf_key, "P E N T A") == NORMAL_RC){
			elem_type = ELEM_PENTA1;
		}else if(check_f06key(buf_key, "H E X A") == NORMAL_RC){
			elem_type = ELEM_HEXA1;
		}else if(check_f06key(buf_key, "T R I A 3") == NORMAL_RC){
			elem_type = ELEM_TRI1;
		}else if(check_f06key(buf_key, "T R I A 6") == NORMAL_RC){
			elem_type = ELEM_TRI2;
		}else if(check_f06key(buf_key, "Q U A D 4") == NORMAL_RC){
			elem_type = ELEM_QUAD1;
		}else if(check_f06key(buf_key, "Q U A D 8") == NORMAL_RC){
			elem_type = ELEM_QUAD2;
		}else{
			fprintf(stderr, "[%s]\n", buf_key);
			return(IMPLEMENT_ERROR_RC);
		}

		RC_TRY( approach_f06strain(fp, elem_type) );

		/* read target strain */
		while(1){
			/* store current file position */
			file_position = ftell(fp);

			if(fgets(buf, sizeof(buf), fp) == NULL){
				fprintf(stderr,"fgets[1] < ");
				return(READ_ERROR_RC);
			}

			if(buf[0] == '1'){
				/* restore position */
				if(fseek(fp, file_position, SEEK_SET) != 0){
					fprintf(stderr,"fseek < ");
					return(SEEK_ERROR_RC);
				}
				break;
			}

			if(elem_type == ELEM_TRI1){
				RC_TRY( realloc_fem_strain_array(strain) );
				RC_TRY( sgetidc(buf, "iidddddddd", idc) );
				elem_id = idc[1].i;
				node_id = -1;
				idc2strain_2d(&(idc[2]), elem_id, node_id,
				              &(strain->array[strain->size]));
				(strain->size)++;

				RC_TRY( fgetidc(fp, "dddddddd", idc) );
				if(f06_ss & F06_SS_FIBRE){
					RC_TRY( realloc_fem_strain_array(strain) );
					idc2strain_2d(idc, elem_id, node_id,
					              &(strain->array[strain->size]));
					(strain->size)++;
				}
			}else if( (elem_type == ELEM_TRI2)
			        ||(elem_type == ELEM_QUAD1)
			        ||(elem_type == ELEM_QUAD2) ){
				if(buf[0] == '0'){
					RC_TRY( sgetidc(buf, "iiidddddddd", idc) );
					elem_id = idc[1].i;
					node_id = -1;
					RC_TRY( realloc_fem_strain_array(strain) );
					idc2strain_2d(&(idc[3]), elem_id, node_id,
					              &(strain->array[strain->size]));
					(strain->size)++;
				}else{
					RC_TRY( sgetidc(buf, "idddddddd", idc) );
					node_id = idc[0].i;
					if(f06_ss & F06_SS_CORNER){
						RC_TRY( realloc_fem_strain_array(strain) );
						idc2strain_2d(&(idc[1]), elem_id, node_id,
						              &(strain->array[strain->size]));
						(strain->size)++;
					}
				}

				RC_TRY( fgetidc(fp, "dddddddd", idc) );
				if(f06_ss & F06_SS_FIBRE){
					RC_TRY( realloc_fem_strain_array(strain) );
					idc2strain_2d(idc, elem_id, node_id,
					              &(strain->array[strain->size]));
					(strain->size)++;
				}

				RC_TRY( skip_file(fp, 1) );
			}else if( (elem_type == ELEM_TETRA1)
			        ||(elem_type == ELEM_PENTA1)
			        ||(elem_type == ELEM_HEXA1) ){
				IDC idc2[6];
				IDC idc3[6];

				if(check_f06key(buf, "GRID CS") == NORMAL_RC){
					RC_TRY( sgetidc(buf, "ii", idc) );
					elem_id = idc[1].i;
					continue;
				}

				if(check_f06key(buf, "CENTER") == NORMAL_RC){
					RC_TRY( realloc_fem_strain_array(strain) );
					node_id = -1;
					RC_TRY( sgetidc(buf, "idddddddd", idc) );
					RC_TRY( fgetidc(fp, "dddddd", idc2) );
					RC_TRY( fgetidc(fp, "dddddd", idc3) );
					idc2strain_3d(&(idc[1]), idc2, idc3, elem_id, node_id,
					              &(strain->array[strain->size]));
					(strain->size)++;
				}else{
					RC_TRY( sgetidc(buf, "iidddddddd", idc) );
					RC_TRY( fgetidc(fp, "dddddd", idc2) );
					RC_TRY( fgetidc(fp, "dddddd", idc3) );
					if(f06_ss & F06_SS_CORNER){
						RC_TRY( realloc_fem_strain_array(strain) );
						node_id = idc[1].i;
						idc2strain_3d(&(idc[2]), idc2, idc3, elem_id, node_id,
						              &(strain->array[strain->size]));
						(strain->size)++;
					}
				}
			}else{
				return(IMPLEMENT_ERROR_RC);
			}
		}
	}

	RC_TRY( clean_fem_strain_array(strain) );

	if(strain->size == 0) return(ZERO_ERROR_RC);
	return(NORMAL_RC);
}


static void idc2strain_2d(const IDC *idc, int elem_id, int node_id,
                          FEM_STRESS_STRAIN *strain)
{
	strain->element = elem_id;
	strain->node = node_id;
	strain->v.x = idc[1].d;
	strain->v.y = idc[2].d;
	strain->v.xy = idc[3].d;

	strain->i_info[0] = 1;
	strain->d_info[0] = idc[0].d; /* FIBRE DISTANCE or STRAIN CURVATURE */
}


static void idc2strain_3d(const IDC *idc, const IDC *idc2, const IDC *idc3,
                          int elem_id, int node_id, FEM_STRESS_STRAIN *strain)
{
	strain->element = elem_id;
	strain->node = node_id;
	strain->v.x = idc[0].d;
	strain->v.y = idc2[0].d;
	strain->v.z = idc3[0].d;
	strain->v.xy = idc[1].d;
	strain->v.yz = idc2[1].d;
	strain->v.zx = idc3[1].d;

	strain->i_info[0] = 2;
	strain->d_info[0] = idc[6].d; /* MEAN PRESSURE */
}


static RC approach_f06strain(FILE *fp, FEM_ELEM_TYPE elem_type)
{
	char buf[BUFFER_SIZE];

	while(1){
		if(fgets(buf, sizeof(buf), fp) == NULL){
			if(feof(fp)){
				return(END_RC);
			}else{
				fprintf(stderr,"fgets[1] < ");
				return(READ_ERROR_RC);
			}
		}
		if( (elem_type == ELEM_TETRA1)
		  ||(elem_type == ELEM_PENTA1)
		  ||(elem_type == ELEM_HEXA1) ){
			if(check_f06_3D(buf) == NORMAL_RC) break;
		}else if( (elem_type == ELEM_TRI1)
		        ||(elem_type == ELEM_TRI2)
		        ||(elem_type == ELEM_QUAD1)
		        ||(elem_type == ELEM_QUAD2) ){
			if(check_f06_2D(buf) == NORMAL_RC) break;
		}else{
			fprintf(stderr, "elem_type < ");
			return(IMPLEMENT_ERROR_RC);
		}
	}

	return(NORMAL_RC);
}


static RC check_f06_2D(char *buf)
{
	if( (check_f06key(buf,"NORMAL-X") == NORMAL_RC)
	  &&(check_f06key(buf,"NORMAL-Y") == NORMAL_RC)
	  &&(check_f06key(buf,"SHEAR-XY") == NORMAL_RC)
	  &&(check_f06key(buf,"ANGLE") == NORMAL_RC)
	  &&(check_f06key(buf,"MAJOR") == NORMAL_RC)
	  &&(check_f06key(buf,"MINOR") == NORMAL_RC)
	  &&(check_f06key(buf,"VON MISES") == NORMAL_RC) ){
		return(NORMAL_RC);
	}

	return(SPECIAL_RC);
}


static RC check_f06_3D(char *buf)
{
	if( (check_f06key(buf,"ELEMENT-ID") == NORMAL_RC)
	  &&(check_f06key(buf,"GRID-ID") == NORMAL_RC)
	  &&(check_f06key(buf,"NORMAL") == NORMAL_RC)
	  &&(check_f06key(buf,"SHEAR") == NORMAL_RC)
	  &&(check_f06key(buf,"PRINCIPAL") == NORMAL_RC)
	  &&(check_f06key(buf,"-A-") == NORMAL_RC)
	  &&(check_f06key(buf,"-B-") == NORMAL_RC)
	  &&(check_f06key(buf,"-C-") == NORMAL_RC)
	  &&(check_f06key(buf,"PRESSURE") == NORMAL_RC)
	  &&(check_f06key(buf,"VON MISES") == NORMAL_RC) ){
		return(NORMAL_RC);
	}

	return(SPECIAL_RC);
}


static RC skip_file(FILE *fp, int record)
{
	char buf[BUFFER_SIZE];
	int ii1;

	for(ii1=0; ii1<record; ii1++){
		if(fgets(buf, sizeof(buf), fp) == NULL){
			if(feof(fp)) return(END_RC);
			return(READ_ERROR_RC);
		}
	}

	return(NORMAL_RC);
}


#if 0
static void set_tensor_stress(FEM_STRESS_ARRAY *stress)
{
	int ii1;

	for(ii1=0; ii1<(stress->size); ii1++){
		stress->array[ii1].tensor[0][0] = stress->array[ii1].v.x;
		stress->array[ii1].tensor[0][1] = stress->array[ii1].tensor[1][0]
		                                = stress->array[ii1].v.xy;
		stress->array[ii1].tensor[0][2] = stress->array[ii1].tensor[2][0]
		                                = stress->array[ii1].v.zx;
		stress->array[ii1].tensor[1][1] = stress->array[ii1].v.y;
		stress->array[ii1].tensor[1][2] = stress->array[ii1].tensor[2][1]
		                                = stress->array[ii1].v.yz;
		stress->array[ii1].tensor[2][2] = stress->array[ii1].v.z;
	}
}
#endif /* 0 */


#if 0
static void set_tensor_strain(FEM_STRAIN_ARRAY *strain)
{
	int ii1;

	for(ii1=0; ii1<(strain->size); ii1++){
		strain->array[ii1].tensor[0][0] = strain->array[ii1].v.x;
		strain->array[ii1].tensor[0][1] = strain->array[ii1].tensor[1][0]
		                                = (strain->array[ii1].v.xy)/2.0;
		strain->array[ii1].tensor[0][2] = strain->array[ii1].tensor[2][0]
		                                = (strain->array[ii1].v.zx)/2.0;
		strain->array[ii1].tensor[1][1] = strain->array[ii1].v.y;
		strain->array[ii1].tensor[1][2] = strain->array[ii1].tensor[2][1]
		                                = (strain->array[ii1].v.yz)/2.0;
		strain->array[ii1].tensor[2][2] = strain->array[ii1].v.z;
	}
}
#endif /* 0 */


