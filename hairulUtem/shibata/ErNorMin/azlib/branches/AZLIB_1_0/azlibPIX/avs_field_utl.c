/*********************************************************************
 * avs_field_utl.c
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

/* $Id: avs_field_utl.c,v 1.10 2003/11/26 08:10:15 sasaoka Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rc.h"
#include "bin_io.h"
#include "string_utl.h"
#include "graphic_utl.h"

#define BUFFER_SIZE (512)
#define STRING_SIZE (512)

/***********************************************************************
 *    AVS_FIELDでのi_infoの内容
 *    
 *    i_info[0] = 0   BYTE 
 *              = 1   SHORT
 *              = 2   INTEGER
 *              = 3   FLOAT
 *              = 4   DOUBLE
 *
 *    i_info[1] = ndim
 *    i_info[2] = nspace
 *    i_info[3] = veclen
 *
 ***********************************************************************/
 
static RC cast_value(GRAPHIC *gr, const char *str1, const char *str2);
#if DEBUG
static void print_value(GRAPHIC *gr);
#endif /* DEBUG */

RC read_avs_field(FILE *fp, GRAPHIC *gr)
{
	RC_TRY(read_avs_field_heder(fp, gr));
	RC_TRY(read_avs_field_data(fp, gr));
	/*print_value(gr);*/
	return(NORMAL_RC);
}


RC read_avs_field_heder(FILE *fp, GRAPHIC *gr)
{
	long int current;
	char tmp[BUFFER_SIZE];
	
	if(fgets(tmp, BUFFER_SIZE, fp) == NULL) return(READ_ERROR_RC);
	if(strncmp(tmp, "# AVS", 5) != 0){
		fprintf(stderr, "not AVS file\n");
		return(READ_ERROR_RC);
	}
	gr->source = GRAPHIC_AVS_FIELD;

	while(1){
		current = ftell(fp);
		if(fgets(tmp, BUFFER_SIZE, fp) == NULL) return(READ_ERROR_RC);
		if(tmp[0] == '#'){
			continue;
		}else if((tmp[0] == 0x000c) && (tmp[1] == 0x000c)){
			if(fseek(fp, current, SEEK_SET) != 0) return(SEEK_ERROR_RC);
			if(fgets(tmp, 3, fp) == NULL) return(READ_ERROR_RC);
			return(NORMAL_RC);
		}else{ /* hoge = fuga を切り分け，値を格納する */ 
			STRING_ARRAY s_array;
			init_string_array(&s_array);
			/* これをやってはいけない labelがくっついてしまう
			RC_TRY( strip_string_char(tmp, " ") );
			*/
			RC_TRY( split_string_delim(&s_array, tmp, "=") );
			RC_TRY( cast_value(gr, s_array.str[0], s_array.str[1]));
			RC_TRY( free_string_array(&s_array) );
		}
	}

	/* dim?の値がそのまま格納されているので，nspaceに合わせて修正 */
	     if(gr->i_info[2] == 1) gr->size.y = gr->size.z = 1;
	else if(gr->i_info[2] == 2) gr->size.z = 1;
	else if(gr->i_info[2] == 3) ;
	else{ 
		fprintf(stderr, "Unknown dimension.\n");
		return(ARG_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC cast_value(GRAPHIC *gr, const char *str1, const char *str2)
{
	IDC idc[10];

	     if(strncmp(str1, "ndim", 4) == 0)   gr->i_info[1] = atoi(str2);
	else if(strncmp(str1, "dim1", 4) == 0)   gr->size.x = atoi(str2);
	else if(strncmp(str1, "dim2", 4) == 0)   gr->size.y = atoi(str2);
	else if(strncmp(str1, "dim3", 4) == 0)   gr->size.z = atoi(str2);
	else if(strncmp(str1, "nspace", 6) == 0) gr->i_info[2] = atoi(str2);
	else if(strncmp(str1, "veclen", 6) == 0){
		gr->i_info[3] = atoi(str2);
		     if(gr->i_info[3] == 1) gr->type = GRAPHIC_SCALAR;
		else if(gr->i_info[3] == 2) gr->type = GRAPHIC_VECT2;
		else if(gr->i_info[3] == 3) gr->type = GRAPHIC_VECT3;
		else if(gr->i_info[3] == 6) gr->type = GRAPHIC_VECT6;
		else{
			fprintf(stderr, "Don't read\n");
			return(READ_ERROR_RC);
		}
	}else if(strncmp(str1, "data", 4) == 0){
		RC_TRY( sgetidc(str2, "s", idc) );
		     if(strncmp(idc[0].s, "byte",    4) == 0) gr->i_info[0] = 0;
		else if(strncmp(idc[0].s, "short",   5) == 0) gr->i_info[0] = 1;
		else if(strncmp(idc[0].s, "integer", 7) == 0) gr->i_info[0] = 2;
		else if(strncmp(idc[0].s, "float",   5) == 0) gr->i_info[0] = 3;
		else if(strncmp(idc[0].s, "double",  6) == 0) gr->i_info[0] = 4;
		else{
			fprintf(stderr,"UNKNOWN DATA TYPE\n");
			return(UNKNOWN_ERROR_RC);
		}
	}else if(strncmp(str1, "field", 5) == 0){
		RC_TRY( sgetidc(str2, "s", idc) );
		if(strncmp(idc[0].s, "uniform", 7) == 0){
			gr->coord_type = GRAPHIC_UNIFORM;
		}else if(strncmp(idc[0].s, "rectilinear", 10) == 0){
			gr->coord_type = GRAPHIC_RECTILINEAR;
		}else if(strncmp(idc[0].s, "irregular", 9) == 0){
			gr->coord_type = GRAPHIC_IRREGULAR;
		}else{
			fprintf(stderr,"UNKNOWN DATA TYPE\n");
			return(UNKNOWN_ERROR_RC);
		}
/* サンプルデータが無い，もしくは必要性が無いために保留
	}else if(strncmp(str1, "min_ext", 7) == 0){
		gr->i_info[4] = gr->i_info[2];
		     if(gr->i_info[2] == 1) RC_TRY( sgetidc(str2, "d",   idc) )
		else if(gr->i_info[2] == 2) RC_TRY( sgetidc(str2, "dd",  idc) )
		else if(gr->i_info[2] == 3) RC_TRY( sgetidc(str2, "ddd", idc) )
		else{
			fprintf(stderr,"UNKNOWN FIELD TYPE");
			return(UNKNOWN_ERROR_RC);
		}
		int max_ext[100];
		for(ii1=0;ii1<gr->i_info[2];ii1++){
			min_ext[ii1] = idc[ii1].d;
		}
	}else if(strncmp(str1, "max_ext", 7) == 0){
		     if(gr->i_info[2] == 1) RC_TRY( sgetidc(str2, "d",   idc) )
		else if(gr->i_info[2] == 2) RC_TRY( sgetidc(str2, "dd",  idc) )
		else if(gr->i_info[2] == 3) RC_TRY( sgetidc(str2, "ddd", idc) )
		else{
			fprintf(stderr,"UNKNOWN FIELD TYPE");
			return(UNKNOWN_ERROR_RC);
		}
		int min_ext[100];
		for(ii1=0;ii1<gr->i_info[2];ii1++){
			max_ext[ii1] = idc[ii1].d;
		}
	}else if(strncmp(str1, "label", 5) == 0){
		int itmp = strlen(str2);
		char *label = NULL;
		if(label == (char *)NULL){
			RC_NULL_CHK( label = (char *)malloc(itmp+1) );
			strncpy(label, str2, itmp);
			label[itmp+1] = '\0';
		}else{
			RC_NULL_CHK( label = (char *)realloc(label, strlen(label)+itmp+1) );
			label[strlen(label)-1] = ' ';
			strncat(label, str2, itmp);
		}
	}else if(strncmp(str1, "unit", 4) == 0){
		char *unit ;
		RC_NULL_CHK( unit = (char *)malloc(STRING_SIZE) );
		strncpy(unit, str2, STRING_SIZE-1);
		unit[STRING_SIZE] = '\0';
	}else if(strncmp(str1, "min_val", 7) == 0){
		char *min_val;
		RC_NULL_CHK( min_val = (char *)malloc(STRING_SIZE) );
		strncpy(min_val, str2, STRING_SIZE-1);
		min_val[STRING_SIZE] = '\0';
	}else if(strncmp(str1, "max_val", 7) == 0){
		char *max_val;
		RC_NULL_CHK( max_val = (char *)malloc(STRING_SIZE) );
		strncpy(max_val, str2, STRING_SIZE-1);
		max_val[STRING_SIZE] = '\0';
	}else if(strncmp(str1, "nstep", 5) == 0){
		int nstep;
		nstep = atoi(str2);
*/
	}else{
		fprintf(stderr, "Wrong file format!!\n");
		return(UNKNOWN_ERROR_RC);
	}
	return(NORMAL_RC);
}


#if DEBUG
static void print_value(GRAPHIC *gr)
{
	fprintf(stderr,"ndim = %d\n", gr->i_info[1]);
	fprintf(stderr,"dim1 = %ld\n", gr->size.x);
	fprintf(stderr,"dim2 = %ld\n", gr->size.y);
	fprintf(stderr,"dim3 = %ld\n", gr->size.z);
	fprintf(stderr,"nspace = %d\n", gr->i_info[2]);
	fprintf(stderr,"veclen = %d\n", gr->i_info[3]);

	     if(gr->i_info[0] == 0) fprintf(stderr, "data = byte\n");
	else if(gr->i_info[0] == 1) fprintf(stderr, "data = short\n");
	else if(gr->i_info[0] == 2) fprintf(stderr, "data = integer\n");
	else if(gr->i_info[0] == 3) fprintf(stderr, "data = float\n");
	else if(gr->i_info[0] == 4) fprintf(stderr, "data = double\n");
	else fprintf(stderr,"UNKNOWN DATA TYPE\n");
	
	switch(gr->coord_type){
	case GRAPHIC_UNIFORM:     fprintf(stderr, "field = uniform\n");     break;
	case GRAPHIC_RECTILINEAR: fprintf(stderr, "field = rectilinear\n"); break;
	case GRAPHIC_IRREGULAR:   fprintf(stderr, "field = irregular\n");   break;
	default: fprintf(stderr, "UNKNOWN FIELD TYPE\n");
	}
}
#endif /* DEBUG */


RC read_avs_field_data(FILE *fp, GRAPHIC *gr)
{
	RC_TRY( read_avs_field_data_value(fp, gr) );
	/* サンプルデータがないために保留
	if(gr->coord_type != UNIFORM){
		RC_TRY( read_avs_field_data_coord(fp, gr) );
	}else{
		fprintf(stderr, "Don't read,yet!\n");
		return(READ_ERROR_RC);
	}
	*/
	return(NORMAL_RC);
}


RC read_avs_field_data_value(FILE *fp, GRAPHIC *gr)
{
	unsigned int ui;
	unsigned long ul;

	switch(gr->type){
	case GRAPHIC_SCALAR:
		RC_TRY( allocate_graphic_data(gr) );
		if(gr->i_info[0] == 0){
			GRAPHIC_SIZE_LOOP_3D(gr->size.z, gr->size.y, gr->size.x){
				if((ui = (unsigned int)fgetc(fp)) == EOF) return(READ_ERROR_RC);
				gr->data.scalar[ii1][ii2][ii3] = (double)ui;
			}GRAPHIC_SIZE_LOOP_3D_END;
		}else if( (gr->i_info[0] == 1) || (gr->i_info[0] == 2) ){
			GRAPHIC_SIZE_LOOP_3D(gr->size.z, gr->size.y, gr->size.x){
				RC_TRY(hread_2bytes(fp, &(ui)) );
				gr->data.scalar[ii1][ii2][ii3] = (double)ui;
			}GRAPHIC_SIZE_LOOP_3D_END;
		}else if( (gr->i_info[0] == 3) || (gr->i_info[0] == 4) ){
			GRAPHIC_SIZE_LOOP_3D(gr->size.z, gr->size.y, gr->size.x){
				RC_TRY(hread_4bytes(fp, &(ul)) );
				gr->data.scalar[ii1][ii2][ii3] = (double)ul;
			}GRAPHIC_SIZE_LOOP_3D_END;
		}else{
			fprintf(stderr, "Unknown Field Data Type\n");
			return(READ_ERROR_RC);
		}
		break;
	default:
		break;
	}

	return(NORMAL_RC);
}


RC read_avs_field_data_coord(FILE *fp, GRAPHIC *gr)
{
	if(gr->i_info[2] <= 0){
		fprintf(stderr,"NO DIMENSION");
		return(UNKNOWN_ERROR_RC);
	}

	RC_TRY( allocate_graphic_data(gr) );

	GRAPHIC_SIZE_LOOP_3D( (unsigned int)(gr->i_info[2]),
	                      (unsigned int)(gr->i_info[1]), gr->size.x){
		gr->data.rgb24[ii1][ii2][ii3].r = 0;
		gr->data.rgb24[ii1][ii2][ii3].g = 0;
		gr->data.rgb24[ii1][ii2][ii3].b = 0;
	}GRAPHIC_SIZE_LOOP_3D_END;

	return(NORMAL_RC);
}


RC write_avs_field_header(FILE *fp, GRAPHIC gr)
{
	fprintf(fp, "# AVS field file \n");
	fprintf(fp, "ndim = %d\n", gr.i_info[1]);

	if(gr.i_info[1] == 0){
		fprintf(stderr, "NO DIMENSION!\n");
		return(READ_ERROR_RC);
	}
	if(gr.i_info[1]>0) fprintf(fp, "dim1 = %ld\n", gr.size.x);
	if(gr.i_info[1]>1) fprintf(fp, "dim2 = %ld\n", gr.size.y);
	if(gr.i_info[1]>2) fprintf(fp, "dim3 = %ld\n", gr.size.z);

	fprintf(fp, "nspace = %d\n", gr.i_info[2]);
	fprintf(fp, "veclen = %d\n", gr.i_info[3]);

	     if(gr.i_info[0] == 0) fprintf(fp, "data = byte\n");
	else if(gr.i_info[0] == 1) fprintf(fp, "data = short\n");
	else if(gr.i_info[0] == 2) fprintf(fp, "data = integer\n");
	else if(gr.i_info[0] == 3) fprintf(fp, "data = float\n");
	else if(gr.i_info[0] == 4) fprintf(fp, "data = double\n");
	else fprintf(fp, "UNKNOWN DATA TYPE\n");

	switch(gr.coord_type){
	case GRAPHIC_UNIFORM:     fprintf(fp, "field = uniform\n");     break;
	case GRAPHIC_RECTILINEAR: fprintf(fp, "field = rectilinear\n"); break;
	case GRAPHIC_IRREGULAR:   fprintf(fp, "field = irregular\n");   break;
	default: fprintf(stderr, "UNKNOWN FIELD TYPE\n");
	}
	/* オプションのサポート サンプルデータが無いので保留
	 * この部分はコメントアウトを外せば使えるものではない
	if(field->min_ext_size != 0){
	fprintf(stderr, "min_ext =");
	for(ii1=0;ii1<field->min_ext_size;ii1++){
		fprintf(stderr, " %f", field->min_ext[ii1]); 
	}
	fprintf(stderr, "\n");
	}
	if(field->max_ext_size != 0){
	fprintf(stderr, "max_ext =");
	for(ii1=0;ii1<field->max_ext_size;ii1++){
		fprintf(stderr, " %f", field->max_ext[ii1]); 
	}
	fprintf(stderr, "\n");
	}
		
	if(field->label != (char *)NULL){
		fprintf(stderr, "label = %s\n", field->label);
	}
	if(field->unit != (char *)NULL){
		fprintf(stderr, "unit = %s\n", field->unit);
	}
	if(field->min_val != (char *)NULL){
		fprintf(stderr, "min_val = %s\n", field->min_val);
	}
	if(field->max_val != (char *)NULL){
		fprintf(stderr, "max_val = %s\n", field->max_val);
	}
	if(field->nstep != 0){
		fprintf(stderr, "nstep = %d\n", field->nstep);
	}
	*/
	fprintf(fp, "\f\f");/*----------セパレータ----------*/
	return(NORMAL_RC);
}


RC write_avs_field_data(FILE *fp, GRAPHIC gr)
{
	unsigned int ui;
	unsigned long ul;

	switch(gr.type){
	case GRAPHIC_SCALAR:
		if(gr.i_info[0] == 0){
			GRAPHIC_SIZE_LOOP_3D(gr.size.z, gr.size.y, gr.size.x){
				fputc((char)gr.data.scalar[ii1][ii2][ii3], fp);
			}GRAPHIC_SIZE_LOOP_3D_END;
		}else if((gr.i_info[0] == 1) || (gr.i_info[0] == 2)){
			GRAPHIC_SIZE_LOOP_3D(gr.size.z, gr.size.y, gr.size.x){
				ui = (unsigned int)gr.data.scalar[ii1][ii2][ii3];
				RC_TRY(lwrite_2bytes(fp, ui) );
			}GRAPHIC_SIZE_LOOP_3D_END;
		}else if((gr.i_info[0] == 3) || (gr.i_info[0] == 4)){
			GRAPHIC_SIZE_LOOP_3D(gr.size.z, gr.size.y, gr.size.x){
				ul = (unsigned long)gr.data.scalar[ii1][ii2][ii3];
				RC_TRY(lwrite_4bytes(fp, ul) );
			}GRAPHIC_SIZE_LOOP_3D_END;
		}else{
			fprintf(stderr, "Unknown Field Data Type\n");
			return(READ_ERROR_RC);
		}
		break;
	case GRAPHIC_RGB24:
		if(gr.i_info[0] == 0){
			GRAPHIC_SIZE_LOOP_3D(gr.size.z, gr.size.y, gr.size.x){
				fputc((char)gr.data.rgb24[ii1][ii2][ii3].g, fp);
			}GRAPHIC_SIZE_LOOP_3D_END;
		}else{
			fprintf(stderr, "Unknown Field Data Type\n");
			return(READ_ERROR_RC);
		}
		break;
	case GRAPHIC_VECT2:
	case GRAPHIC_VECT3:
	case GRAPHIC_VECT6:
		fprintf(stderr, "Don't read,yet!\n");
		return(READ_ERROR_RC);
		break;
	default:
		break;
	}

	return(NORMAL_RC);
}

