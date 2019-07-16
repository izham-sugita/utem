/*********************************************************************
 * sysnoise_component.c
 *
 * Copyright (C) 2001 AzLib Developers Group
 *
 * Written by
 *  <Yasuhiro MATSUURA> <Masanobu SUYAMA>
 *
 * Thanks to following contributors.
 *  <Kenzen TAKEUCHI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: sysnoise_component.c 415 2005-07-15 06:04:24Z sasaoka $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rc.h"
#include "sysnoise_component.h"
#include "fem_struct.h"

#define BUFFER_SIZE (512)

static RC seek_str( FILE *fp, char *str, int num );
static RC sysnoise_output_node( FILE *fp, NODE_ARRAY node );
static RC elem_type( int type_num, ELEM_TYPE *type );
static RC set_line( ELEMENT *element, int *node_num );
static RC set_tri( ELEMENT *element, int *node_num );
static RC set_quad( ELEMENT *element, int *node_num );
static RC set_tetra( ELEMENT *element, int *node_num );
static RC set_penta( ELEMENT *element, int *node_num );
static RC set_hexa( ELEMENT *element, int *node_num );
static RC sysnoise_output_element( FILE *fp, ELEMENT_ARRAY element );
static RC put_line( FILE *fp, ELEMENT element );
static RC put_tri( FILE *fp, ELEMENT element );
static RC put_quad( FILE *fp, ELEMENT element );
static RC put_tetra( FILE *fp, ELEMENT element );
static RC put_penta( FILE *fp, ELEMENT element );
static RC put_hexa( FILE *fp, ELEMENT element );
static int elem_type_number( ELEM_TYPE type );


RC
sysnoise_input_pressure ( FILE *fp, SOUND_PRESSURE_ARRAY *sound )
{
	int num, label;
	double real, imag, magnitude, phase, db;
	char buf[BUFFER_SIZE];

	RC_TRY( allocate_sound_pressure_array(0, sound) );

	rewind(fp);
	while(1){
		if(fgets( buf, sizeof(buf), fp ) == NULL){
			if(feof(fp)){
				break;
			} else {
				return(READ_ERROR_RC);
			}
		}

		RC_TRY( realloc_sound_pressure_array(sound) );

		if(sscanf( buf, " %d %d %le %le %le %le %le ", 
		           &num, &label, &real, &imag, &magnitude,
		           &phase, &db ) == 7){
			sound->array[sound->size].node = label;
			sound->array[sound->size].re = real;  /* real part */
			sound->array[sound->size].im = imag;  /* imaginary part */
			(sound->size)++;
		}
	}

	RC_TRY( clean_sound_pressure_array(sound) );

	return(NORMAL_RC);
}


static RC
seek_str ( FILE *fp, char *str, int num )
{
	char buf[BUFFER_SIZE];

	rewind(fp);
	while(1){
		if(fgets( buf, sizeof(buf), fp ) == NULL){
			if(feof(fp)){
				return(END_RC);
			} else {
				return(READ_ERROR_RC);
			}
		}
		if(strncmp( buf, str, num ) == 0){
			break;
		}
	}

	return(NORMAL_RC);
}


RC
sysnoise_input_node ( FILE *fp, NODE_ARRAY *node )
{
	int ii1, num, node_num, elem_num, elem_type;
	char buf[BUFFER_SIZE];

	rewind(fp);
	for(ii1=0; ii1<5; ii1++){
		RC_NULL_CHK( fgets(buf, sizeof(buf), fp) ); 
	}

	if(sscanf( buf, "%d %d %d ", &node_num, &elem_num, &elem_type ) != 3){
		fprintf(stderr, "Error : sscanf 1\n");
		return(READ_ERROR_RC);
	}

	RC_TRY( allocate_node_array(node_num, node) );
	node->size = node_num;
	node->source = FEM_SYSNOISE;

	RC_TRY( seek_str(fp, "NODES", 5) );

	for(ii1=0; ii1<node->size; ii1++){
		RC_NULL_CHK( fgets(buf, sizeof(buf), fp) ); 

		if(sscanf( buf, "%d %d %le %le %le ",
		           &num, &node->array[ii1].label, &node->array[ii1].p.x,
		           &node->array[ii1].p.y, &node->array[ii1].p.z ) != 5){
			fprintf(stderr, "Error : sscanf 2\n");
			return(READ_ERROR_RC);
		}
	}

	RC_TRY( clean_node_array(node) );

	return(NORMAL_RC);
}


static RC
sysnoise_output_node ( FILE *fp, NODE_ARRAY node )
{
	int ii1;

	fprintf(fp, "NODES\n");
	for(ii1=0; ii1<node.size; ii1++){
		if(node.array[ii1].label < 0) continue;

		fprintf(fp, " %9d %9d %19.8E %19.8E %19.8E\n",
		        ii1+1, node.array[ii1].label,
		        node.array[ii1].p.x, node.array[ii1].p.y,
		        node.array[ii1].p.z);
	}

	return(NORMAL_RC);
}


RC
sysnoise_input_element ( FILE *fp, ELEMENT_ARRAY *element )
{
	int ii1, num, typenum;
	int node[20], node_num, elem_num, type;
	char buf[BUFFER_SIZE];

	rewind(fp);
	for(ii1=0; ii1<5; ii1++){
		RC_NULL_CHK( fgets(buf, sizeof(buf), fp) ); 
	}
	
	if(sscanf( buf, "%d %d %d ", &node_num, &elem_num, &type ) != 3){
		fprintf(stderr, "Error : sscanf 1\n");
		return(READ_ERROR_RC);
	}

	RC_TRY( allocate_element_array(elem_num, element) );
	element->size = elem_num;
	element->source = FEM_SYSNOISE;

	RC_TRY( seek_str(fp, "ELEMENTS", 8) );

	for(ii1=0; ii1<element->size; ii1++){
		element->array[ii1].material = 1; 
		element->array[ii1].physical = 1; 
		RC_NULL_CHK( fgets(buf, sizeof(buf), fp) ); 

		if(sscanf( buf, "%d %d %d %d %d %d %d %d",
		           &num, &element->array[ii1].label, &typenum,
		           &element->array[ii1].node_num,
		           &node[0], &node[1], &node[2], &node[3] ) < 6){
			fprintf(stderr, "Error : sscanf 2\n");
			return(READ_ERROR_RC);
		}
   
		RC_TRY( elem_type(typenum, &(element->array[ii1].type)) );

		if(element->array[ii1].node_num > 4){
			RC_NULL_CHK( fgets(buf, sizeof(buf), fp) ); 

			if(sscanf( buf, "%d %d %d %d %d %d %d %d",
			           &node[4], &node[5], &node[6],
			           &node[7], &node[8], &node[9],
			           &node[10], &node[11] ) == 0){
				fprintf(stderr, "Error : sscanf 3\n");
				return(READ_ERROR_RC);
			}
			if(element->array[ii1].node_num > 12){
				RC_NULL_CHK( fgets(buf, sizeof(buf), fp) ); 

				if(sscanf( buf, "%d %d %d %d %d %d %d %d",
				           &node[12], &node[13], &node[14],
			    	       &node[15], &node[16], &node[17],
			        	   &node[18], &node[19] ) == 0){
					fprintf(stderr, "Error : sscanf 4\n");
					return(READ_ERROR_RC);
				}
			}
		}

		switch(element->array[ii1].type){
		case ELEM_LINE1:
		case ELEM_LINE2:
			RC_TRY( set_line(&(element->array[ii1]), node) );
			break;
		case ELEM_TRI1:
		case ELEM_TRI2:
			RC_TRY( set_tri(&(element->array[ii1]), node) );
			break;
		case ELEM_QUAD1:
		case ELEM_QUAD2:
			RC_TRY( set_quad(&(element->array[ii1]), node) );
			break;
		case ELEM_TETRA1:
		case ELEM_TETRA2:
			RC_TRY( set_tetra(&(element->array[ii1]), node) );
			break;
		case ELEM_PENTA1:
		case ELEM_PENTA2:
			RC_TRY( set_penta(&(element->array[ii1]), node) );
			break;
		case ELEM_HEXA1:
		case ELEM_HEXA2:
			RC_TRY( set_hexa(&(element->array[ii1]), node) );
			break;
		default:
			return(IMPLEMENT_ERROR_RC);
        }
	}

	RC_TRY( clean_element_array(element) );

	return(NORMAL_RC);
}


static RC
elem_type ( int type_num, ELEM_TYPE *type )
{
	switch(type_num){
	case 2:
		*type = ELEM_LINE1;
		break;
	case 3:
		*type = ELEM_LINE2;
		break;
	case 4:
		*type = ELEM_TRI1;
		break;
	case 5:
		*type = ELEM_TRI2;
		break;
	case 6:
		*type = ELEM_QUAD1;
		break;
	case 7:
		*type = ELEM_QUAD2;
		break;
	case 8:
		*type = ELEM_TETRA1;
		break;
	case 9:
		*type = ELEM_TETRA2;
		break;
	case 10:
		*type = ELEM_HEXA1;
		break;
	case 11:
		*type = ELEM_HEXA2;
		break;
	case 12:
		*type = ELEM_PENTA1;
		break;
	case 13:
		*type = ELEM_PENTA2;
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
set_line ( ELEMENT *element, int *node_num )
{

	if(element->type == ELEM_LINE1){
		element->node[0] = node_num[0];
		element->node[1] = node_num[1];
	}else if(element->type == ELEM_LINE2){
		element->node[0] = node_num[0];
		element->node[2] = node_num[1];
		element->node[1] = node_num[2];
	}else{
		return(ARG_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
set_tri ( ELEMENT *element, int *node_num )
{

	if(element->type == ELEM_TRI1){
		element->node[0] = node_num[0];
		element->node[1] = node_num[1];
		element->node[2] = node_num[2];
	}else if(element->type == ELEM_TRI2){
		element->node[0] = node_num[0];
		element->node[5] = node_num[1];
		element->node[1] = node_num[2];
		element->node[3] = node_num[3];
		element->node[2] = node_num[4];
		element->node[4] = node_num[5];
	}else{
		return(ARG_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
set_quad ( ELEMENT *element, int *node_num )
{

	if(element->type == ELEM_QUAD1){
		element->node[0] = node_num[0];
		element->node[1] = node_num[1];
		element->node[2] = node_num[2];
		element->node[3] = node_num[3];
	}else if(element->type == ELEM_QUAD2){
		element->node[0] = node_num[0];
		element->node[4] = node_num[1];
		element->node[1] = node_num[2];
		element->node[5] = node_num[3];
		element->node[2] = node_num[4];
		element->node[6] = node_num[5];
		element->node[3] = node_num[6];
		element->node[7] = node_num[7];
	}else{
		return(ARG_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
set_tetra ( ELEMENT *element, int *node_num )
{

	if(element->type == ELEM_TETRA1){
		element->node[0] = node_num[0];
		element->node[1] = node_num[1];
		element->node[2] = node_num[2];
		element->node[3] = node_num[3];
	}else if(element->type == ELEM_TETRA2){
		element->node[0] = node_num[0];
		element->node[6] = node_num[1];
		element->node[1] = node_num[2];
		element->node[4] = node_num[3];
		element->node[2] = node_num[4];
		element->node[5] = node_num[5];
		element->node[7] = node_num[6];
		element->node[8] = node_num[7];
		element->node[9] = node_num[8];
		element->node[3] = node_num[9];
	}else{
		return(ARG_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
set_penta ( ELEMENT *element, int *node_num )
{

	if(element->type == ELEM_PENTA1){
		element->node[0] = node_num[0];
		element->node[1] = node_num[1];
		element->node[2] = node_num[2];
		element->node[3] = node_num[3];
		element->node[4] = node_num[4];
		element->node[5] = node_num[5];
	}else if(element->type == ELEM_PENTA2){
		element->node[0] = node_num[0];
		element->node[8] = node_num[1];
		element->node[1] = node_num[2];
		element->node[6] = node_num[3];
		element->node[2] = node_num[4];
		element->node[7] = node_num[5];
		element->node[12] = node_num[6];
		element->node[13] = node_num[7];
		element->node[14] = node_num[8];
		element->node[3] = node_num[9];
		element->node[11] = node_num[10];
		element->node[4] = node_num[11];
		element->node[9] = node_num[12];
		element->node[5] = node_num[13];
		element->node[10] = node_num[14];
	}else{
		return(ARG_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
set_hexa ( ELEMENT *element, int *node_num )
{

	if(element->type == ELEM_HEXA1){
		element->node[0] = node_num[0];
		element->node[1] = node_num[1];
		element->node[2] = node_num[2];
		element->node[3] = node_num[3];
		element->node[4] = node_num[4];
		element->node[5] = node_num[5];
		element->node[6] = node_num[6];
		element->node[7] = node_num[7];
	}else if(element->type == ELEM_HEXA2){
		element->node[0] = node_num[0];
		element->node[8] = node_num[1];
		element->node[1] = node_num[2];
		element->node[9] = node_num[3];
		element->node[2] = node_num[4];
		element->node[10] = node_num[5];
		element->node[3] = node_num[6];
		element->node[11] = node_num[7];
		element->node[16] = node_num[8];
		element->node[17] = node_num[9];
		element->node[18] = node_num[10];
		element->node[19] = node_num[11];
		element->node[4] = node_num[12];
		element->node[12] = node_num[13];
		element->node[5] = node_num[14];
		element->node[13] = node_num[15];
		element->node[6] = node_num[16];
		element->node[14] = node_num[17];
		element->node[7] = node_num[18];
		element->node[15] = node_num[19];
	}else{
		return(ARG_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
sysnoise_output_element ( FILE *fp, ELEMENT_ARRAY element )
{
	int ii1, typenum;

	fprintf(fp, "ELEMENTS\n"); 
	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		typenum = elem_type_number(element.array[ii1].type); 

		fprintf(fp, " %9d %9d %9d %9d",
		        ii1+1, element.array[ii1].label,
		        typenum, element.array[ii1].node_num);

		switch(element.array[ii1].type){
		case ELEM_LINE1:
		case ELEM_LINE2:
			RC_TRY( put_line(fp, element.array[ii1]) );
			break;
		case ELEM_TRI1:
		case ELEM_TRI2:
			RC_TRY( put_tri(fp, element.array[ii1]) );
			break;
		case ELEM_QUAD1:
		case ELEM_QUAD2:
			RC_TRY( put_quad(fp, element.array[ii1]) );
			break;
		case ELEM_TETRA1:
		case ELEM_TETRA2:
			RC_TRY( put_tetra(fp, element.array[ii1]) );
			break;
		case ELEM_PENTA1:
		case ELEM_PENTA2:
			RC_TRY( put_penta(fp, element.array[ii1]) );
			break;
		case ELEM_HEXA1:
		case ELEM_HEXA2:
			RC_TRY( put_hexa(fp, element.array[ii1]) );
			break;
		default:
			return(IMPLEMENT_ERROR_RC);
        }

	}

	return(NORMAL_RC);
}


static RC
put_line ( FILE *fp, ELEMENT element )
{

	if(element.type == ELEM_LINE1){
		fprintf(fp, " %9d", element.node[0]);
		fprintf(fp, " %9d", element.node[1]);
		fprintf(fp, "\n");
	}else if(element.type == ELEM_LINE2){
		fprintf(fp, " %9d", element.node[0]);
		fprintf(fp, " %9d", element.node[2]);
		fprintf(fp, " %9d", element.node[1]);
		fprintf(fp, "\n");
	}else{
		return(ARG_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
put_tri ( FILE *fp, ELEMENT element )
{

	if(element.type == ELEM_TRI1){
		fprintf(fp, " %9d", element.node[0]);
		fprintf(fp, " %9d", element.node[1]);
		fprintf(fp, " %9d", element.node[2]);
		fprintf(fp, "\n");
	}else if(element.type == ELEM_TRI2){
		fprintf(fp, " %9d", element.node[0]);
		fprintf(fp, " %9d", element.node[5]);
		fprintf(fp, " %9d", element.node[1]);
		fprintf(fp, " %9d", element.node[3]);
		fprintf(fp, "\n");
		fprintf(fp, " %9d", element.node[2]);
		fprintf(fp, " %9d", element.node[4]);
		fprintf(fp, "\n");
	}else{
		return(ARG_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
put_quad ( FILE *fp, ELEMENT element )
{

	if(element.type == ELEM_QUAD1){
		fprintf(fp, " %9d", element.node[0]);
		fprintf(fp, " %9d", element.node[1]);
		fprintf(fp, " %9d", element.node[2]);
		fprintf(fp, " %9d", element.node[3]);
		fprintf(fp, "\n");
	}else if(element.type == ELEM_QUAD2){
		fprintf(fp, " %9d", element.node[0]);
		fprintf(fp, " %9d", element.node[4]);
		fprintf(fp, " %9d", element.node[1]);
		fprintf(fp, " %9d", element.node[5]);
		fprintf(fp, "\n");
		fprintf(fp, " %9d", element.node[2]);
		fprintf(fp, " %9d", element.node[6]);
		fprintf(fp, " %9d", element.node[3]);
		fprintf(fp, " %9d", element.node[7]);
		fprintf(fp, "\n");
	}else{
		return(ARG_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
put_tetra ( FILE *fp, ELEMENT element )
{

	if(element.type == ELEM_TETRA1){
		fprintf(fp, " %9d", element.node[0]);
		fprintf(fp, " %9d", element.node[1]);
		fprintf(fp, " %9d", element.node[2]);
		fprintf(fp, " %9d", element.node[3]);
		fprintf(fp, "\n");
	}else if(element.type == ELEM_TETRA2){
		fprintf(fp, " %9d", element.node[0]);
		fprintf(fp, " %9d", element.node[6]);
		fprintf(fp, " %9d", element.node[1]);
		fprintf(fp, " %9d", element.node[4]);
		fprintf(fp, "\n");
		fprintf(fp, " %9d", element.node[2]);
		fprintf(fp, " %9d", element.node[5]);
		fprintf(fp, " %9d", element.node[7]);
		fprintf(fp, " %9d", element.node[8]);
		fprintf(fp, " %9d", element.node[9]);
		fprintf(fp, " %9d", element.node[3]);
		fprintf(fp, "\n");
	}else{
		return(ARG_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
put_penta ( FILE *fp, ELEMENT element )
{

	if(element.type == ELEM_PENTA1){
		fprintf(fp, " %9d", element.node[0]);
		fprintf(fp, " %9d", element.node[1]);
		fprintf(fp, " %9d", element.node[2]);
		fprintf(fp, " %9d", element.node[3]);
		fprintf(fp, "\n");
		fprintf(fp, " %9d", element.node[4]);
		fprintf(fp, " %9d", element.node[5]);
		fprintf(fp, "\n");
	}else if(element.type == ELEM_PENTA2){
		fprintf(fp, " %9d", element.node[0]);
		fprintf(fp, " %9d", element.node[8]);
		fprintf(fp, " %9d", element.node[1]);
		fprintf(fp, " %9d", element.node[6]);
		fprintf(fp, "\n");
		fprintf(fp, " %9d", element.node[2]);
		fprintf(fp, " %9d", element.node[7]);
		fprintf(fp, " %9d", element.node[12]);
		fprintf(fp, " %9d", element.node[13]);
		fprintf(fp, " %9d", element.node[14]);
		fprintf(fp, " %9d", element.node[3]);
		fprintf(fp, " %9d", element.node[11]);
		fprintf(fp, " %9d", element.node[4]);
		fprintf(fp, "\n");
		fprintf(fp, " %9d", element.node[9]);
		fprintf(fp, " %9d", element.node[5]);
		fprintf(fp, " %9d", element.node[10]);
		fprintf(fp, "\n");
	}else{
		return(ARG_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
put_hexa ( FILE *fp, ELEMENT element )
{

	if(element.type == ELEM_HEXA1){
		fprintf(fp, " %9d", element.node[0]);
		fprintf(fp, " %9d", element.node[1]);
		fprintf(fp, " %9d", element.node[2]);
		fprintf(fp, " %9d", element.node[3]);
		fprintf(fp, "\n");
		fprintf(fp, " %9d", element.node[4]);
		fprintf(fp, " %9d", element.node[5]);
		fprintf(fp, " %9d", element.node[6]);
		fprintf(fp, " %9d", element.node[7]);
		fprintf(fp, "\n");
	}else if(element.type == ELEM_HEXA2){
		fprintf(fp, " %9d", element.node[0]);
		fprintf(fp, " %9d", element.node[8]);
		fprintf(fp, " %9d", element.node[1]);
		fprintf(fp, " %9d", element.node[9]);
		fprintf(fp, "\n");
		fprintf(fp, " %9d", element.node[2]);
		fprintf(fp, " %9d", element.node[10]);
		fprintf(fp, " %9d", element.node[3]);
		fprintf(fp, " %9d", element.node[11]);
		fprintf(fp, " %9d", element.node[16]);
		fprintf(fp, " %9d", element.node[17]);
		fprintf(fp, " %9d", element.node[18]);
		fprintf(fp, " %9d", element.node[19]);
		fprintf(fp, "\n");
		fprintf(fp, " %9d", element.node[4]);
		fprintf(fp, " %9d", element.node[12]);
		fprintf(fp, " %9d", element.node[5]);
		fprintf(fp, " %9d", element.node[13]);
		fprintf(fp, " %9d", element.node[6]);
		fprintf(fp, " %9d", element.node[14]);
		fprintf(fp, " %9d", element.node[7]);
		fprintf(fp, " %9d", element.node[15]);
		fprintf(fp, "\n");
	}else{
		return(ARG_ERROR_RC);
	}

	return(NORMAL_RC);
}


static int
elem_type_number ( ELEM_TYPE type )
{

	switch(type){
	case ELEM_LINE1:
		return(2);
	case ELEM_LINE2:
		return(3);
	case ELEM_TRI1:
		return(4);
	case ELEM_TRI2:
		return(5);
	case ELEM_QUAD1:
		return(6);
	case ELEM_QUAD2:
		return(7);
	case ELEM_TETRA1:
		return(8);
	case ELEM_TETRA2:
		return(9);
	case ELEM_HEXA1:
		return(10);
	case ELEM_HEXA2:
		return(11);
	case ELEM_PENTA1:
		return(12);
	case ELEM_PENTA2:
		return(13);
	default:
		break;
	}

	return(-1);
}


RC
sysnoise_output_model_file ( FILE *fp, NODE_ARRAY node,
                             ELEMENT_ARRAY element )
{
	int ii1, nodenum;

	fprintf(fp, "SYSNOISE MESH FILE\n");
	fprintf(fp, "Rev 5.4\n");
	fprintf(fp, "SYSNOISE Default Model\n");
	fprintf(fp, "\n");

	nodenum = 0;
	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		if(nodenum < element.array[ii1].node_num){
			nodenum = element.array[ii1].node_num;
		}
	}
	fprintf(fp, " %9d %9d %9d\n", count_valid_node(node),
	                              count_valid_element(element),
	                              nodenum );

	RC_TRY( sysnoise_output_node(fp, node) );
	RC_TRY( sysnoise_output_element(fp, element) );

	return(NORMAL_RC);
}


/* SYSNOISE コマンドファイル(.cmd)作成用 */
RC
sysnoise_output_option ( FILE *fp, char *name, int model_num, int dim )
{ 
	fprintf(fp, "New Name '%s' Model %d File \"%s.sdb\" Return \n",
	                                         name, model_num, name);
	fprintf(fp, "Option FEM Frequency Fluid Return \n");
	fprintf(fp, "Import Mesh Format Free File \"femmodel.fre\" Return \n");

	if(dim == 2) 
		fprintf(fp, "TwoDimensional Return \n");

	fprintf(fp, "Renumber \n");
	fprintf(fp, "    Node All \n");
	fprintf(fp, "    Return \n");
	fprintf(fp, "Check Mesh Return \n");
	fprintf(fp, "Material Fluid \n");
	fprintf(fp, "    Name 'Air' \n"); 
/*	fprintf(fp, "    Sound 3.4000e+02  Imaginary 3.400e+00  Rho 1.2250e+00 \n");
*/	fprintf(fp, "    Sound 3.4000e+02  Rho 1.2250e+00 \n");
	fprintf(fp, "    Elements All \n");
	fprintf(fp, "    Return \n");

	return(NORMAL_RC);
}


RC
sysnoise_output_point ( FILE *fp )
{ 

	fprintf(fp, "Point Mesh Format Free File \"receive.fre\" Return \n");

	return(NORMAL_RC);
}


RC
sysnoise_output_parameter ( FILE *fp, int a_num )
{

	fprintf(fp,  "Parameter Model 1 \n");
	fprintf(fp,  "    Vector 0 \n");
	fprintf(fp,  "    Save Potentials Step 1 \n");
	if(a_num == 1) 
		fprintf(fp, "    Save Results Step 1 \n");
	else 
		fprintf(fp, "    Save Results none \n");
	fprintf(fp, "    Store Potentials none \n");
	fprintf(fp, "    Return \n");
	fprintf(fp, "    Store Results none \n");
	fprintf(fp, "    Return \n");
	fprintf(fp, "    NoFlow \n");
	fprintf(fp, "    Return \n");
		
	return(NORMAL_RC);
}


RC
sysnoise_output_solve_pres ( FILE *fp, char *name, double minfreq,
                             double maxfreq, double linestep, int a_num )
{
	int ii1;
	int f_count;
	double freq;

	if(minfreq == maxfreq){
		fprintf(fp, "Solve Frequency %f Return \n", minfreq); 
	} else{
		fprintf(fp, "Solve Frequency %f To %f LinStep %f Return \n", 
		                                  minfreq, maxfreq, linestep); 
	}

	if(nearly_eq(minfreq, maxfreq)){
		f_count = 1;
	} else{
		freq = (maxfreq - minfreq) / linestep;
		f_count = (int)freq + 1;
	}

	freq = minfreq;
	for(ii1=0; ii1<f_count; ii1++){
		fprintf(fp, "Search Potentials Frequency %f Return\n", freq);
		fprintf(fp, "LogFile File %s%d.dat Return \n", name, ii1);
		fprintf(fp, "Extract Pressure Nodes All Return \n");
/*		if(a_num == 1 && strcmp(name, "pres") == 0){
*/		if(a_num == 1){
			fprintf(fp, "Search Results Frequency %f Return\n", freq);
			fprintf(fp, "LogFile File rec_%s%d.dat Return \n", name, ii1);
			fprintf(fp, "Extract Pressure Point All Return \n");
		}
		freq += linestep;
	}
		
	return(NORMAL_RC);
}


RC
sysnoise_output_solve_freq_data ( FILE *fp, char *name, double freq,
                                  int counter, int a_num )
{

	fprintf(fp, "Solve Frequency %f Return \n", freq); 
	fprintf(fp, "Search Potentials Frequency %f Return\n", freq);
	fprintf(fp, "LogFile File %s%d.dat Return \n", name, counter);
	fprintf(fp, "Extract Pressure Nodes All Return \n");
	if(a_num == 1){
		fprintf(fp, "Search Results Frequency %f Return\n", freq);
		fprintf(fp, "LogFile File rec_%s%d.dat Return \n", name, counter);
		fprintf(fp, "Extract Pressure Point All Return \n");
	}
		
	return(NORMAL_RC);
}


RC
sysnoise_output_solve_mode ( FILE *fp, char *name, double min, double max )
{

	fprintf(fp, "Modes \n"); 
	fprintf(fp, "    From %f To %f Tolerance 1.0000e-6 Iteration 100 \n",
	                                                            min, max);
	fprintf(fp, "    BLanczos \n");
	fprintf(fp, "    Return \n");
	fprintf(fp, "LogFile File  %s.dat Return \n", name);
	fprintf(fp, "Extract Frequency Return \n");
		
	return(NORMAL_RC);
}


RC
sysnoise_output_structure_freq ( FILE *fp, double minfreq, double maxfreq,
                                 double linestep )
{
	if(minfreq == maxfreq){
		fprintf(fp,  "    Frequency %f \n", minfreq); 
	} else{
		fprintf(fp,  "    Frequency %f To %f LinStep %f \n",
		                         minfreq, maxfreq, linestep); 
	}

    return(NORMAL_RC);
}


RC
sysnoise_output_structure_mesh ( FILE *fp, char *name )
{
	fprintf(fp,  "    Mesh File  %s.fre Format Free \n", name);
	fprintf(fp,  "    Algorithm  1  Tolerance 0.100000 Average 1 \n");
	fprintf(fp,  "    Return \n");

    return(NORMAL_RC);
}


RC
sysnoise_output_structure_data ( FILE *fp, char *name )
{
	fprintf(fp,  "Generate \n");
	fprintf(fp,  "    Face Set 1 \n");
	fprintf(fp,  "    From Displacements File  %s.pch Format Nastran \n", name);

    return(NORMAL_RC);
}


RC
sysnoise_output_set_face ( FILE *fp, char *name, int set_num,
                           SOUND_PRESSURE_ARRAY sound )
{
	int ii1;
	int counter = 0;

	fprintf(fp, "Set  %d Name  \" %s\" \n", set_num, name);
	fprintf(fp, "    Faces");
	for(ii1=0; ii1<sound.size; ii1++){
		if(sound.array[ii1].node < 0) continue;

		if(counter == 20){
			fprintf(fp, " \n");
			fprintf(fp, "    Faces");
			counter = 0;
		}
		fprintf(fp, " %d", sound.array[ii1].node);
		counter++;
	}
	fprintf(fp, " \n");
	fprintf(fp, "    Return \n");

    return(NORMAL_RC);
}


RC
sysnoise_output_structure_option ( FILE *fp, char *name, int model_num,
                                   double thickness )
{ 
	fprintf(fp, "New Name '%s' Model %d File \"%s.sdb\" Return \n",
	                                         name, model_num, name);
	fprintf(fp, "Option FEM Frequency Structure Return \n");
	fprintf(fp, "Import Mesh Format Free File \"panel.fre\" Return \n");
	fprintf(fp, "Renumber \n");
	fprintf(fp, "    Node All \n");
	fprintf(fp, "    Return \n");
	fprintf(fp, "Check Mesh Return \n");
	fprintf(fp, "Material Shell \n");
	fprintf(fp, "    Name 'shell' \n"); 
	fprintf(fp, "    Young 2.1e11  Poisson 0.3  Rho 7800 \n");
	fprintf(fp, "    Elements All \n");
	fprintf(fp, "    Return \n");
	fprintf(fp, "Geometry Thickness %f \n", thickness);
	fprintf(fp, "    Elements All \n");
	fprintf(fp, "    Return \n");

	return(NORMAL_RC);
}


RC
sysnoise_output_bc_option ( FILE *fp, int set_num )
{ 
	fprintf(fp, "Set %d Name  \"Envelope\" Envelope Return \n", set_num);
	fprintf(fp, "Boundary \n");
	fprintf(fp, "     UX Real 0 Imag 0 UY Real 0 Imag 0 UZ Real 0 Imag 0 \n");
	fprintf(fp, "    Nodes Set %d \n", set_num);
	fprintf(fp, "    Return \n");

	return(NORMAL_RC);
}


RC
sysnoise_output_vector_mode ( FILE *fp, char *name, int mode_num,
                              int damping_flag )
{
	int ii1;

	fprintf(fp, "Modes \n"); 
	fprintf(fp, "    Vector %d Tolerance 1.0000e-6 Iteration 100 \n", mode_num);
	fprintf(fp, "    Lanczos \n");
	fprintf(fp, "    Return \n");
	fprintf(fp, "LogFile File  %s.dat Return \n", name);
	fprintf(fp, "Extract Frequency Return \n");

	if(damping_flag == 1){
		for(ii1=0; ii1<mode_num; ii1++){
			fprintf(fp, "Damping Mode %d KSI: 0.020000 Return \n", ii1+1);
		}
	}
		
	return(NORMAL_RC);
}


RC
sysnoise_output_load_option ( FILE *fp, BC_ARRAY force )
{ 
	int ii1;

	fprintf(fp, "Boundary \n");
	for(ii1=0; ii1<force.size; ii1++){
		if(force.array[ii1].node < 0) continue;

		if(force.array[ii1].v_type.x == BC_FORCE){
			fprintf(fp, "     FX Real %f Imag 0 \n", force.array[ii1].v.x);
		}
		if(force.array[ii1].v_type.y == BC_FORCE){
			fprintf(fp, "     FY Real %f Imag 0 \n", force.array[ii1].v.y);
		}
		if(force.array[ii1].v_type.z == BC_FORCE){
			fprintf(fp, "     FZ Real %f Imag 0 \n", force.array[ii1].v.z);
		}
		fprintf(fp, "    Nodes %d \n", force.array[ii1].node);
	}
	fprintf(fp, "    Return \n");

	return(NORMAL_RC);
}


RC
sysnoise_output_save_option ( FILE *fp, char *name )
{ 
	fprintf(fp, "Save File  %s.sdb Return \n", name);

	return(NORMAL_RC);
}


RC
sysnoise_output_end_option ( FILE *fp, char *name )
{ 
	fprintf(fp, "Exit Save Journal \"%s.jnl\" Save Models \n", name);

	return(NORMAL_RC);
}


RC
sysnoise_output_link_option ( FILE *fp, int struct_mode_num,
                              int fluid_mode_num, int set_num, int a_num )
{ 
	fprintf(fp, "Link \n");
	fprintf(fp, "    Model 1 \n");
	fprintf(fp, "    Elements All \n");
	fprintf(fp, "    To \n");
	fprintf(fp, "    Model 2 \n");
	fprintf(fp, "    Faces Set %d \n", set_num);
	fprintf(fp, "    Behavior FLUID-STRUCTURE \n");
	fprintf(fp, "    Return \n");
	fprintf(fp, "Parameter Model 1 \n");
	fprintf(fp, "    Vector %d \n", struct_mode_num);
	fprintf(fp, "    Save Displacements Step 1 \n");
	fprintf(fp, "    Store Displacements none \n");
	fprintf(fp, "    Return \n");
	fprintf(fp, "    Return \n");
	fprintf(fp, "Parameter Model 2 \n");
	fprintf(fp, "    Vector %d \n", fluid_mode_num);
	fprintf(fp, "    Save Potentials Step 1 \n");
	if(a_num == 1) 
		fprintf(fp, "    Save Results Step 1 \n");
	else 
		fprintf(fp, "    Save Results none \n");
	fprintf(fp, "    Store Potentials none \n");
	fprintf(fp, "    Return \n");
	fprintf(fp, "    Store Results none \n");
	fprintf(fp, "    Return \n");
	fprintf(fp, "    NoFlow \n");
	fprintf(fp, "    Return \n");

	return(NORMAL_RC);
}


RC
sysnoise_output_file_end ( FILE *fp, char *name )
{

	fprintf(fp, "Save Return \n");
	fprintf(fp, "Close Return \n");
	fprintf(fp, "Exit Save Journal \"%s.jnl\" NoSave Models \n", name);
	fprintf(fp, "\n");

	return(NORMAL_RC);
}


RC
sysnoise_input_source_point ( FILE *fp, NODE_ARRAY *node,
                              SOUND_PRESSURE_ARRAY *sound )
{
	int label;
	double x, y, z, real, imag;
	char buf[BUFFER_SIZE];

	RC_TRY( allocate_node_array(0, node) );
	RC_TRY( allocate_sound_pressure_array(0, sound) );

	rewind(fp);
	RC_TRY( seek_str(fp, " SOURCE", 7) );

	while(1){
		if(fgets( buf, sizeof(buf), fp ) == NULL){
			if(feof(fp)){
				break;
			} else{
				return(READ_ERROR_RC);
			}
		}

		RC_TRY( realloc_node_array(node) );
		RC_TRY( realloc_sound_pressure_array(sound) );

		if(sscanf( buf, " Source Id Number   %d ",&label ) == 1){
			sound->array[sound->size].node = label;
			node->array[node->size].label = label;
			RC_NULL_CHK( fgets(buf, sizeof(buf), fp) ); 

			if(sscanf( buf, " Source amplitude real : %le   imag : %le ",
			            &real, &imag ) == 2){
				sound->array[sound->size].re = real;
				sound->array[sound->size].im = imag;
			}
			while(1){
				RC_NULL_CHK( fgets(buf, sizeof(buf), fp) ); 

				if(sscanf( buf, " Point Source Position :"
					        "  X = %le  Y = %le  Z = %le ",
				                         &x, &y, &z ) == 3){
					node->array[node->size].p.x = x;
					node->array[node->size].p.y = y;
					node->array[node->size].p.z = z;
					break;
				}
			}
		(sound->size)++;
		(node->size)++;
		}
	}

	RC_TRY( clean_node_array(node) );
	RC_TRY( clean_sound_pressure_array(sound) );

	return(NORMAL_RC);
}


RC
sysnoise_output_source_point ( FILE *fp, NODE_ARRAY node,
                               SOUND_PRESSURE_ARRAY sound, int dim )
{
	int ii1, index;

	for(ii1=0; ii1<node.size; ii1++){
		if(node.array[ii1].label < 0) continue;

		RC_NEG_CHK( index = search_sound_pressure_node(sound,
		                                 node.array[ii1].label) );

		if(dim == 2){
			fprintf(fp, "Source Cylindrical \n");
		} else if(dim == 3){
			fprintf(fp, "Source Spherical \n");
		} else{
			return(READ_ERROR_RC);
		}
		fprintf(fp, "    Amplitude Real %15.7E Imag %15.7E \n",
		        sound.array[index].re, sound.array[index].im );
		fprintf(fp, "    Position %15.7E %15.7E %15.7E \n",
		        node.array[ii1].p.x, node.array[ii1].p.y,
		        node.array[ii1].p.z);
		if(dim == 2){
			fprintf(fp, "    Vector 0 0 1 \n");
		}
		fprintf(fp, "Return \n");
	}

	return(NORMAL_RC);
}


RC
sysnoise_input_boundary_pressure ( FILE *fp, SOUND_PRESSURE_ARRAY *sound )
{
	int num, inter_num, exter_num;
	double real, imag;
	char buf[BUFFER_SIZE];

	RC_TRY( allocate_sound_pressure_array(0, sound) );

	rewind(fp);
	RC_TRY( seek_str(fp, " PRESCRIBED NODAL PRESSURE", 26) );

	while(1){
		if(fgets( buf, sizeof(buf), fp ) == NULL){
			if(feof(fp)){
				break;
			} else{
				return(READ_ERROR_RC);
			}
		}

		RC_TRY( realloc_sound_pressure_array(sound) );

		if(sscanf( buf, " %d %d %d %le %le ",
		           &num, &inter_num, &exter_num, &real, &imag ) == 5){
			sound->array[sound->size].node = exter_num;
			sound->array[sound->size].re = real;
			sound->array[sound->size].im = imag;
			(sound->size)++;
		} else{
			if(strncmp( buf, " PRESCRIBED", 11 ) == 0) break;
		}
	}

	RC_TRY( clean_sound_pressure_array(sound) );

	return(NORMAL_RC);
}


RC
sysnoise_output_boundary_pressure ( FILE *fp, SOUND_PRESSURE_ARRAY sound )
{
	int ii1;

	for(ii1=0; ii1<sound.size; ii1++){
		if(sound.array[ii1].node < 0) continue;

		fprintf(fp, "Boundary Pressure Real %E Imag %E \n",
		         sound.array[ii1].re, sound.array[ii1].im);
		fprintf(fp, "    Nodes %4d \n", sound.array[ii1].node);
		fprintf(fp, "    Return \n");
	}

	return(NORMAL_RC);
}


RC
sysnoise_input_boundary_velocity ( FILE *fp, SOUND_PRESSURE_ARRAY *sound )
{
	int num, facenum, elemnum, nodenum;
	double real, imag;
	char buf[BUFFER_SIZE];

	RC_TRY( allocate_sound_pressure_array(0, sound) );

	rewind(fp);
	RC_TRY( seek_str(fp, " PRESCRIBED FACE VELOCITY", 25) );

	while(1){
		if(fgets( buf, sizeof(buf), fp ) == NULL){
			if(feof(fp)){
				break;
			} else{
				return(READ_ERROR_RC);
			}
		}

		RC_TRY( realloc_sound_pressure_array(sound) );

		if(sscanf( buf," %d %d %d %d %le %le ",
		           &num, &facenum, &elemnum, &nodenum,
		           &real, &imag ) == 6){
			sound->array[sound->size].node = facenum;
			sound->array[sound->size].re = real;
			sound->array[sound->size].im = imag;
			(sound->size)++;
		} else{
			if(strncmp( buf, " PRESCRIBED", 11 ) == 0) break;
		}
	}

	RC_TRY( clean_sound_pressure_array(sound) );

	return(NORMAL_RC);
}


RC
sysnoise_output_boundary_velocity ( FILE *fp, SOUND_PRESSURE_ARRAY sound )
{
	int ii1;

	for(ii1=0; ii1<sound.size; ii1++){
		if(sound.array[ii1].node < 0) continue;

		fprintf(fp, "Boundary Velocity Real %E Imag %E \n",
		         sound.array[ii1].re, sound.array[ii1].im);
		fprintf(fp, "    Faces %4d \n", sound.array[ii1].node);
		fprintf(fp, "    Return \n");
	}

	return(NORMAL_RC);
}


RC
sysnoise_input_absorbent_panels ( FILE *fp, SOUND_PRESSURE_ARRAY *sound )
{
	int num, facenum, elemnum;
	double real, imag;
	char buf[BUFFER_SIZE];

	RC_TRY( allocate_sound_pressure_array(0, sound) );

	rewind(fp);
	RC_TRY( seek_str(fp, " PRESCRIBED FACE ADMITTANCE", 27) );

	while(1){
		if(fgets( buf, sizeof(buf), fp ) == NULL){
			if(feof(fp)) {
				break;
			} else {
				return(READ_ERROR_RC);
			}
		}

		RC_TRY( realloc_sound_pressure_array(sound) );

		if(sscanf( buf, " %d %d %d %le %le ",
		           &num, &facenum, &elemnum, &real, &imag ) == 5){
			sound->array[sound->size].node = facenum;
			sound->array[sound->size].re = real;
			sound->array[sound->size].im = imag;
			(sound->size)++;
		} else{
			if(strncmp( buf, " PRESCRIBED", 11 ) == 0) break;
		}
	}

	RC_TRY( clean_sound_pressure_array(sound) );

	return(NORMAL_RC);
}


RC
sysnoise_output_absorbent_panels ( FILE *fp, SOUND_PRESSURE_ARRAY sound )
{
	int ii1;

	for(ii1=0; ii1<sound.size; ii1++){
		if(sound.array[ii1].node < 0) continue;

		fprintf(fp, "Boundary Admittance Real %E Imag %E \n",
		        sound.array[ii1].re, sound.array[ii1].im);
		fprintf(fp, "    Faces %4d \n", sound.array[ii1].node);
		fprintf(fp, "    Return \n");
	}

	return(NORMAL_RC);
}


RC
sysnoise_input_symmetry_plane ( FILE *fp, VECT3D *p, int *counter )
{
	double x, y, z;
	char buf[BUFFER_SIZE];

	p->x = p->y = p->z = FSYSNOISEDEFAULT;
	*counter = 0;

	RC_NULL_CHK( fgets(buf, sizeof(buf), fp) ); 
	
	if(sscanf( buf," x = %le ", &x  ) == 1){
		p->x = x;
		*counter += 1;
	}
	RC_NULL_CHK( fgets(buf, sizeof(buf), fp) ); 

	if(sscanf( buf," y = %le ", &y  ) == 1){ 
		p->y = y;
		*counter += 2;
	}
	RC_NULL_CHK( fgets(buf, sizeof(buf), fp) ); 

	if(sscanf( buf," z = %le ", &z  ) == 1){ 
		p->z = z;
		*counter += 4;
	}

	return(NORMAL_RC);
}


RC
sysnoise_output_symmetry_plane ( FILE *fp, VECT3D p, int counter )
{
	switch(counter){
	case 0:
		break;
	case 1:
		fprintf(fp, "Symmetry No Return Symmetry Plane X = %f Return\n", p.x);
		break;
	case 2:
		fprintf(fp, "Symmetry No Return Symmetry Plane Y = %f Return\n", p.y);
		break;
	case 3:
		fprintf(fp, "Symmetry No Return Symmetry Plane X = %f Return "
		                   "Symmetry Plane Y = %f Return\n", p.x, p.y);
		break;
	case 4:
		fprintf(fp, "Symmetry No Return Symmetry Plane Z = %f Return\n", p.z);
		break;
	case 5:
		fprintf(fp, "Symmetry No Return Symmetry Plane X = %f Return "
		                   "Symmetry Plane Z = %f Return\n", p.x, p.z);
		break;
	case 6:
		fprintf(fp, "Symmetry No Return Symmetry Plane Y = %f Return "
		                   "Symmetry Plane Z = %f Return\n", p.y, p.z);
		break;
	case 7:
		fprintf(fp, "Symmetry No Return Symmetry Plane X = %f Return "
		       "Symmetry Plane Y = %f Return Symmetry Plane Z = %f Return \n",
		        p.x, p.y, p.z);
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


RC
sysnoise_input_eigen_frequency ( FILE *fp, double *eigen_freq )
{
	int label;
	int counter;
	double freq, phase;
	char buf[BUFFER_SIZE];

	rewind(fp);
	RC_TRY( seek_str(fp, " ACOUSTIC EIGEN FREQUENCIES", 27) );

	counter = 0;
	while(1){
		if(fgets( buf, sizeof(buf), fp ) == NULL){
			if(feof(fp)){
				break;
			} else{
				return(READ_ERROR_RC);
			}
		}

		if(sscanf(buf, " MODE  %d     FREQUENCY (Hz)  %le      PHASE   %le",
		          &label, &freq, &phase) == 3){
			eigen_freq[counter] = freq - 0.2;
			counter++;
			eigen_freq[counter] = freq;
			counter++;
			eigen_freq[counter] = freq + 0.2;
			counter++;
		}
	}

	return(NORMAL_RC);
}


int
count_eigen_frequency ( FILE *fp )
{
	int label;
	int count = 0;
	double freq, phase;
	char buf[BUFFER_SIZE];

	rewind(fp);
	while(1){
		if(fgets( buf, sizeof(buf), fp ) == NULL){
			if(feof(fp)){
				fprintf(stderr, "EIGEN FREQUENCY ZERO\n");
				return(count);
			}
		}
		if(strncmp( buf, " ACOUSTIC EIGEN FREQUENCIES", 27 ) == 0){
			break;
		}
	}

	while(1){
		if(fgets( buf, sizeof(buf), fp ) == NULL){
			if(feof(fp)) break;
		}

		if(sscanf(buf, " MODE  %d     FREQUENCY (Hz)  %le      PHASE   %le",
		          &label, &freq, &phase) == 3){
			count++;
		}
	}

/*	return(count);
*/	return(count*3);
}


RC
sysnoise_input_couple_eigen_frequency ( FILE *fp, double *eigen_freq,
                                        double min, double max )
{
	int label;
	int counter;
	double freq, phase;
	char buf[BUFFER_SIZE];

	rewind(fp);
	RC_TRY( seek_str(fp, " COUPLED EIGEN FREQUENCIES", 26) );

	counter = 0;
	while(1){
		if(fgets( buf, sizeof(buf), fp ) == NULL){
			if(feof(fp)){
				break;
			} else{
				return(READ_ERROR_RC);
			}
		}

		if(sscanf(buf, " MODE  %d     FREQUENCY (Hz)  %le    EIGEN VALUE  %le",
		          &label, &freq, &phase) == 3){
			if(min < freq && freq < max){
				eigen_freq[counter] = freq;
				counter++;
			}
		}
	}

	return(NORMAL_RC);
}


int
count_couple_eigen_frequency ( FILE *fp, double min, double max )
{
	int label;
	int count = 0;
	double freq, phase;
	char buf[BUFFER_SIZE];

	rewind(fp);
	while(1){
		if(fgets( buf, sizeof(buf), fp ) == NULL){
			if(feof(fp)){
				fprintf(stderr, "EIGEN FREQUENCY ZERO\n");
				return(count);
			}
		}
		if(strncmp( buf, " COUPLED EIGEN FREQUENCIES", 26 ) == 0){
			break;
		}
	}

	while(1){
		if(fgets( buf, sizeof(buf), fp ) == NULL){
			if(feof(fp)) break;
		}

		if(sscanf(buf, " MODE  %d     FREQUENCY (Hz)  %le    EIGEN VALUE  %le",
		          &label, &freq, &phase) == 3){
			if(min < freq && freq < max){
				count++;
			}
		}
	}

	return(count);
}


