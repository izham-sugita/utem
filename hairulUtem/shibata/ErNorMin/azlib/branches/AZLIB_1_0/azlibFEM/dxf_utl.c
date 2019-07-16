/*********************************************************************
 * dxf_utl.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Ryu SASAOKA> <Takaaki NAGATANI>
 *
 * Thanks to following contributors.
 *  <Tomoyuki TSUBATA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: dxf_utl.c,v 1.11 2003/11/04 06:55:35 nagatani Exp $ */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "rc.h"
#include "string_utl.h"
#include "fem_struct.h"
#include "dxf_utl.h"

static char *dxf_section = "SECTION";
static char *dxf_3dface = "3DFACE";
static char *dxf_vertex = "VERTEX";
static char *dxf_polyline = "POLYLINE";
static char *dxf_eof = "EOF";
static char *default_layer_name = "UNKNOWN LAYER";

#define BDIST        (1.0E-9)
#define BUFFER_SIZE  (512)
#define INIT_SIZE	 (1024)


/* local function */
static RC output_dxf_3dface_base(FILE *fp, char *def_layer, int def_color,
                                 char **layer, int num_layer,
                                 FEM_ELEMENT_ARRAY elem, FEM_NODE_ARRAY node);
static void write_dxf_section(FILE *fp);
static void write_dxf_vertex(FILE *fp);
static void write_dxf_end(FILE *fp);
static void write_dxf_3dface(FILE *fp, double *nodes, char *layer, int color);
static RC read_dxf_3dface(FILE *fp, char **layer, int *num_node,
                                         double *xyz, int *color);
static void clean_str(char *str);
static int seek_0(FILE *fp);
static int layer_index(int *numlayer, char ***layer, char *str);
static int compar_node(const void *p1, const void *p2);



/* output to DXF(3DFACE) file */
RC output_dxf_3dface_lay(FILE *fp, char **layer, int num_layer,
                         FEM_ELEMENT_ARRAY elem, FEM_NODE_ARRAY node)
{
	return( output_dxf_3dface_base(fp, default_layer_name,
	                               DXF_HUE_RED+DXF_BRIGHT_NORMAL,
	                               layer, num_layer, elem, node) );
}

RC output_dxf_3dface(FILE *fp, char *def_layer, int def_color,
                     FEM_ELEMENT_ARRAY elem, FEM_NODE_ARRAY node)
{
	return( output_dxf_3dface_base(fp, def_layer, def_color,
	                               NULL, 0, elem, node) );
}

static RC output_dxf_3dface_base(FILE *fp, char *def_layer, int def_color,
                                 char **layer, int num_layer,
                                 FEM_ELEMENT_ARRAY elem, FEM_NODE_ARRAY node)
{
	int ii1,ii2;
	int index;
	double xyz[12];
	char *lay = NULL;
	int col = 0;
	
	if((def_layer == NULL)||(def_color < 0)||(def_color > 259)){
		return(ARG_ERROR_RC);
	}

	/* print headder of DXF */
	write_dxf_section(fp);

	for(ii1=0;ii1<elem.size;ii1++){
		if(elem.array[ii1].label < 0)continue;
		switch(elem.array[ii1].type){
		case ELEM_TRI1:
		case ELEM_TRI2:
			for(ii2=0;ii2<3;ii2++){
				index = search_fem_node_label(node, elem.array[ii1].node[ii2]);
				if(index < 0){
					fprintf(stderr,
					        "search_fem_node_label[%d|%d] < ", ii1, ii2);
					return(SEEK_ERROR_RC);
				}
				xyz[3*ii2  ] = node.array[index].p.x;
				xyz[3*ii2+1] = node.array[index].p.y;
				xyz[3*ii2+2] = node.array[index].p.z;
			}

			/* repeat last one */
			xyz[3*3  ] = node.array[index].p.x;
			xyz[3*3+1] = node.array[index].p.y;
			xyz[3*3+2] = node.array[index].p.z;
			break;
		case ELEM_QUAD1:
		case ELEM_QUAD2:
			for(ii2=0;ii2<4;ii2++){
				index = search_fem_node_label(node, elem.array[ii1].node[ii2]);
				if(index < 0){
					fprintf(stderr,
					        "search_fem_node_label[%d|%d] < ", ii1, ii2);
					return(SEEK_ERROR_RC);
				}
				xyz[3*ii2  ] = node.array[index].p.x;
				xyz[3*ii2+1] = node.array[index].p.y;
				xyz[3*ii2+2] = node.array[index].p.z;
			}
			break;
		default:
			return(IMPLEMENT_ERROR_RC);
		}
		if(elem.source == FEM_DXF){
			if( (elem.array[ii1].i_info[0] < 0)
			  ||(elem.array[ii1].i_info[0] > 259) ){
				fprintf(stderr, "elem.array[%d].i_info[0]: %d < \n",
				        ii1, elem.array[ii1].i_info[0]);
				return(ARG_ERROR_RC);
			}
			col = elem.array[ii1].i_info[0];

			if(layer != NULL){
				if( (elem.array[ii1].i_info[1] < 0)
			      ||(elem.array[ii1].i_info[1] >= num_layer) ){
					fprintf(stderr, "elem.array[%d].i_info[1]: %d/%d < \n",
					        ii1, elem.array[ii1].i_info[1], num_layer);
					return(ARG_ERROR_RC);
				}
				lay = layer[elem.array[ii1].i_info[1]];
			}else{
				lay = def_layer;
			}
		}else{
			col = def_color;
			lay = def_layer;
		}

		write_dxf_3dface(fp, xyz, lay, col);
	}

	/* print footer of DXF */
	write_dxf_end(fp);

	return(NORMAL_RC);
}


/* output to DXF(POLYLINE) file */
RC output_dxf_polyline(FILE *fp, FEM_ELEMENT_ARRAY element, FEM_NODE_ARRAY node)
{
	return( output_dxf_polyline_lay(fp, "0", element, node) );
}

RC output_dxf_polyline_lay(FILE *fp, char *layer,
                            FEM_ELEMENT_ARRAY elem, FEM_NODE_ARRAY node)
{
	int ii1;

	/* print headder of DXF(POLYLINE) */
	write_dxf_section(fp);
	fprintf(fp, "0\n" "POLYLINE\n" "  8\n" "%s\n", layer);
	fprintf(fp, "66\n" "1\n" "70\n" "64\n");
	fprintf(fp, "71\n" "%d\n" "72\n" "%d\n", node.size, elem.size);
	
	/* print node_data */
	for(ii1 = 0; ii1 < node.size; ii1++){
		write_dxf_vertex(fp);
		fprintf(fp, "%s\n", layer);
		fprintf(fp, "10\n" "%14.7e\n", node.array[ii1].p.x);
		fprintf(fp, "20\n" "%14.7e\n", node.array[ii1].p.y);
		fprintf(fp, "30\n" "%14.7e\n", node.array[ii1].p.z);
		fprintf(fp, "70\n" "192\n");
	}

	/* print element_data */
	for(ii1 = 0; ii1 < elem.size; ii1++){
		switch(elem.array[ii1].type){
		case ELEM_TRI1:
		case ELEM_TRI2:
		case ELEM_QUAD1:
		case ELEM_QUAD2:
			write_dxf_vertex(fp);
			fprintf(fp, "%s\n", layer);
			fprintf(fp, "10\n" "0\n" "20\n" "0\n" "30\n" "0\n");
			fprintf(fp, "70\n" "128\n");

			fprintf(fp, "71\n" " %d\n", elem.array[ii1].node[0]);
			fprintf(fp, "72\n" " %d\n", elem.array[ii1].node[1]);
			fprintf(fp, "73\n" " %d\n", elem.array[ii1].node[2]);
			if(elem.array[ii1].type == ELEM_QUAD1 ||
			   elem.array[ii1].type == ELEM_QUAD2){ /* QUAD の場合1つ追加 */
				write_dxf_vertex(fp);
				fprintf(fp, "%s\n", layer);
				fprintf(fp, "10\n" "0\n" "20\n" "0\n" "30\n" "0\n");
				fprintf(fp, "70\n" "128\n");

				fprintf(fp, "71\n" " %d\n", elem.array[ii1].node[0]);
				fprintf(fp, "72\n" " %d\n", elem.array[ii1].node[2]);
				fprintf(fp, "73\n" " %d\n", elem.array[ii1].node[3]);
			}
			break;
		default:
			return(IMPLEMENT_ERROR_RC);
		}
	}

	/* print footer of DXF */
	fprintf(fp, "  0\n" "SEQEND\n");
	write_dxf_end(fp);

	return(NORMAL_RC);
}

static void write_dxf_section(FILE *fp)
{
	fprintf(fp, "  0\n"
	            "%s\n"
	            "  2\n"
	            "ENTITIES\n", dxf_section);
}

static void write_dxf_vertex(FILE *fp)
{
	fprintf(fp, "  0\n"
    	        "%s\n"
        	    "  8\n", dxf_vertex);
}

static void write_dxf_end(FILE *fp)
{
	fprintf(fp, "  0\n"
	            "ENDSEC\n"
	            "  0\n"
	            "%s\n", dxf_eof);
}


static void write_dxf_3dface(FILE *fp, double *xyz, char *layer, int color)
{
	int ii1, ii2;
	int index;

	fprintf(fp, "  0\n"
	            "%s\n"
	            "  8\n"
	            "%s\n"
	            " 62\n"
	            "%3d\n", dxf_3dface, layer, color);
	
	index = 0;
	for(ii1=0;ii1<4;ii1++){
		for(ii2=0;ii2<3;ii2++){
			fprintf(fp,"%3d\n"
			           "%14.7e\n", 10 * (ii2 + 1) + ii1, xyz[index]);
			index++;
		}
	}
}

RC input_dxf(FILE *fp, FEM_ELEMENT_ARRAY *element, FEM_NODE_ARRAY *node)
{
	int ret_zero, dummy;
	int face_flag = 0, poly_flag = 0;

	while(1){
		if((ret_zero = seek_0(fp)) < 0){
			break;
		} else if(ret_zero == 2){	 /* 3DFACE */
			face_flag = 1;
		} else if(ret_zero == 4){	 /* POLYLINE */
			poly_flag = 1; 
		}
	}

	rewind(fp);
	if(face_flag == 1 && poly_flag == 0){
		return( input_dxf_3dface_lay(fp,  NULL, &dummy, element, node) );
	} else if(face_flag == 0 && poly_flag == 1){
		return( input_dxf_polyline_lay(fp,  NULL, &dummy, element, node) );
	} else {
		fprintf(stderr, "NO_SURPPORT_FILE_FORMAT < \n");
		return(IMPLEMENT_ERROR_RC);
	}
	return(NORMAL_RC);
}


RC input_dxf_3dface(FILE *fp, FEM_ELEMENT_ARRAY *element, FEM_NODE_ARRAY *node)
{
	int dummy;

	return( input_dxf_3dface_lay(fp,  NULL, &dummy, element, node) );
}


/* DXF ファイルより、要素、節点データを読み込む (3DFACE版)*/
RC input_dxf_3dface_lay(FILE *fp, char ***layer, int *num_layer,
                        FEM_ELEMENT_ARRAY *element, FEM_NODE_ARRAY *node)
{
	int ii1, ii2;
	RC rc;
	int nn, col;
	double xyz[12];
	int l_index;
	int counter;
	char *lay;

	RC_TRY( allocate_fem_element_array(INIT_SIZE, element) );
	RC_TRY( allocate_fem_node_array(INIT_SIZE, node) );
	element->source = FEM_DXF;
	node->source = FEM_DXF;
	if(layer != NULL){
		*num_layer = 0;
		*layer = NULL;
	}

	while(1){
		rc = read_dxf_3dface(fp, &lay, &nn, xyz, &col);
		if(rc == END_RC) break;
		if(rc != NORMAL_RC){
			fprintf(stderr, "read_dxf_3dface() < ");
			return(rc);
		}

		RC_TRY( realloc_fem_element_array(element) );

		/* put element */
		element->array[element->size].label = (element->size) + 1;
		element->array[element->size].i_info[0] = col;
		if(nn == 3){
			element->array[element->size].type = ELEM_TRI1;
			element->array[element->size].node_num = 3;
		}else if(nn == 4){
			element->array[element->size].type = ELEM_QUAD1;
			element->array[element->size].node_num = 4;
		}else{
			fprintf(stderr, "???");
			return(UNKNOWN_ERROR_RC);
		}

		/* layer が NULL の場合はレイヤ名を記憶しない */
		if(layer != NULL){
			l_index = layer_index(num_layer, layer, lay);
			if(l_index < 0){
				fprintf(stderr, "layer_index() < ");
				return(ALLOC_ERROR_RC);
			}
			element->array[element->size].i_info[1] = l_index;
		}else{
			element->array[element->size].i_info[1] = -1;
		}

		/* put node */
		for(ii1=0; ii1<nn; ii1++){
			RC_TRY( realloc_fem_node_array(node) );
			node->array[node->size].label = (node->size) + 1;
			node->array[node->size].p.x= xyz[3*ii1];
			node->array[node->size].p.y= xyz[3*ii1+1];
			node->array[node->size].p.z= xyz[3*ii1+2];
			node->array[node->size].i_info[0] = element->size;  /* 逆引き用 */
			node->array[node->size].i_info[1] = ii1;            /* 逆引き用 */
			element->array[element->size].node[ii1] = (node->size) + 1;
			(node->size)++;
		}
		(element->size)++;
	}
	qsort((void *)(node->array), node->size, sizeof(FEM_NODE), compar_node);

	counter = 0;
	for(ii1=0; ii1<node->size; ii1++){
		int n_index = ii1;

		for(ii2=ii1-1; ii2>=0; ii2--){
			double dx, dy, dz;

			if(node->array[ii2].label < 0) continue;
			dx = node->array[ii2].p.x - node->array[ii1].p.x;
			if(fabs(dx) > BDIST){
				break;
			}
			dy = node->array[ii2].p.y - node->array[ii1].p.y;
			dz = node->array[ii2].p.z - node->array[ii1].p.z;
			if(sqrt(dx*dx + dy*dy + dz*dz) <= BDIST){
				n_index = ii2;
				break;
			}
		}
		if(n_index == ii1){
			counter++;
			node->array[ii1].label = counter;
			element->array[ node->array[ii1].i_info[0] ]
			         .node[ node->array[ii1].i_info[1] ] = counter;
			node->array[ii1].i_info[0] = 0;
			node->array[ii1].i_info[1] = 0;
		}else{
			element->array[ node->array[ii1].i_info[0] ]
			         .node[ node->array[ii1].i_info[1] ]
			= node->array[n_index].label;
			node->array[ii1].label = -1;
		}
	}
	RC_TRY( clean_fem_element_array(element) );
	RC_TRY( clean_fem_node_array(node) );

	return(NORMAL_RC);
}


static int compar_node(const void *p1, const void *p2)
{
	if( ((FEM_NODE *)p1)->p.x < ((FEM_NODE *)p2)->p.x ) return(-1);
	if( ((FEM_NODE *)p1)->p.x > ((FEM_NODE *)p2)->p.x ) return(1);

	if( ((FEM_NODE *)p1)->p.y < ((FEM_NODE *)p2)->p.y ) return(-1);
	if( ((FEM_NODE *)p1)->p.y > ((FEM_NODE *)p2)->p.y ) return(1);

	if( ((FEM_NODE *)p1)->p.z < ((FEM_NODE *)p2)->p.z ) return(-1);
	if( ((FEM_NODE *)p1)->p.z > ((FEM_NODE *)p2)->p.z ) return(1);
	
	return(0);
}


static int layer_index(int *numlayer, char ***layer, char *str)
{
	int ii1;
	char **lay_tmp;

	if(*numlayer < 0) return(-1);
	if(layer == NULL) return(-1);

	for(ii1=0;ii1<(*numlayer);ii1++){
		if(strcmp((*layer)[ii1], str) == 0) return(ii1);
	}

	(*numlayer)++;
	lay_tmp = (char **)realloc((*layer), (*numlayer)*sizeof(char *));
	if(lay_tmp == NULL){
		fprintf(stderr, "realloc() < ");
		return(-1);
	}
	(*layer) = lay_tmp;

	(*layer)[(*numlayer) - 1] = (char *)malloc(strlen(str)+1);
	if((*layer)[(*numlayer) - 1] == NULL){
		fprintf(stderr, "malloc() < ");
		return(-1);
	}
	strcpy((*layer)[(*numlayer) - 1], str);

	return((*numlayer) - 1);
}


/* 3DFACE を1個読み込む */
/* layer には、静的配列のポインタが代入される */
static RC read_dxf_3dface(FILE *fp, char **layer, int *num_node,
                                         double *xyz, int *color)
{
	int ret_zero;
	int ii1;
	RC rc;
	IDC idc;
	int num;
	long file_position;
	double dx, dy, dz;
	double dist;
	static char lay_buf[128];

	while(1){
		if((ret_zero = seek_0(fp)) < 0) return(END_RC);
		if(ret_zero == 2) break;  /* 3DFACE */
	}

	for(ii1=0;ii1<(4*3);ii1++){
		xyz[ii1] = 0.0;
	}
	*color = 256;
	*layer = default_layer_name;
	*num_node = 0;

	while(1){
		file_position = ftell(fp);
		rc = fgetidc(fp, "i", &idc);
		if(rc == END_RC) return(END_RC);
		if(rc != NORMAL_RC){
			fprintf(stderr, "fgetidc1 < ");
			return(rc);
		}

		num = idc.i;
		switch(num){
		case 8:
			if(fgets(lay_buf, sizeof(lay_buf), fp) == NULL){
				fprintf(stderr, "fgets() < ");
				return(READ_ERROR_RC);
			}
			clean_str(lay_buf);
			(*layer) = lay_buf;
			break;
		case 62:
			rc = fgetidc(fp, "i", &idc);
			if(rc != NORMAL_RC){
				fprintf(stderr, "fgetidc3 < ");
				return(rc);
			}
			(*color) = idc.i;
			break;
		case 10:
		case 20:
		case 30:
		case 11:
		case 21:
		case 31:
		case 12:
		case 22:
		case 32:
		case 13:
		case 23:
		case 33:
			rc = fgetidc(fp, "d", &idc);
			if(rc != NORMAL_RC){
				fprintf(stderr, "fgetidc4 < ");
				return(rc);
			}
			xyz[(3 * (num % 10)) + (num / 10 - 1)] = idc.d;
			break;
		case 0:
			if(fseek(fp, file_position, SEEK_SET) != 0){
				fprintf(stderr, "fseek < ");
				return(SEEK_ERROR_RC);
			}
			dx = xyz[3*3+0] - xyz[3*2+0];
			dy = xyz[3*3+1] - xyz[3*2+1];
			dz = xyz[3*3+2] - xyz[3*2+2];
			dist = sqrt(dx*dx + dy*dy + dz*dz);
			if(dist < BDIST){
				(*num_node) = 3;
			}else{
				(*num_node) = 4;
			}
			return(NORMAL_RC);
		default:
			break;
		}
	}
}

RC input_dxf_polyline(FILE *fp, FEM_ELEMENT_ARRAY *element, FEM_NODE_ARRAY *node)
{
	int dummy;

	return( input_dxf_polyline_lay(fp,  NULL, &dummy, element, node) );
}


/* DXF ファイルより、要素、節点データを読み込む (POLYLINE版)*/
RC input_dxf_polyline_lay(FILE *fp, char ***layer, int *num_layer,
                        FEM_ELEMENT_ARRAY *element, FEM_NODE_ARRAY *node)
{
	int ret_zero;
	int ii1;
	int num;
	int flag;
	IDC idc;
	long file_position;
	double xyz[3];
	int elem_no[3];
	static char lay_buf[128];

	RC_TRY( allocate_fem_element_array(INIT_SIZE, element) );
	RC_TRY( allocate_fem_node_array(INIT_SIZE, node) );
	element->source = FEM_DXF;
	node->source = FEM_DXF;
	if(layer != NULL){
		*num_layer = 0;
		*layer = NULL;
	}
	for(ii1=0; ii1 < 3; ii1++){
		xyz[ii1] = 0.0;
		elem_no[ii1] = 0;
	}

	/* POLYLINEを読み込む */
	while(1){
		flag = 0;
		while(1){
			if((ret_zero = seek_0(fp)) < 0) return(NORMAL_RC);
			if(ret_zero == 3) break;  /* VERTEX */
		}
		while(1){
			file_position = ftell(fp);
			RC_TRY( fgetidc(fp, "i", &idc) );

			num = idc.i;
			switch(num){
			case 8:
				if(fgets(lay_buf, sizeof(lay_buf), fp) == NULL){
					fprintf(stderr, "fgets() < ");
					return(READ_ERROR_RC);
				}
				/* layer が NULL の場合はレイヤ名を記憶しない */
				if(layer != NULL){ 
					clean_str(lay_buf);
					*(layer[(*num_layer)]) = lay_buf;
					(*num_layer)++;
				}
				break;
			case 70:	/* node,elementを判定 */
				RC_TRY( fgetidc(fp, "i", &idc) );

				if(idc.i == 192){		/* put node */
					(node->size)++;
					node->array[node->size - 1].label = (node->size);
					RC_TRY( realloc_fem_node_array(node) );
					node->array[node->size - 1].p.x= xyz[0];
					node->array[node->size - 1].p.y= xyz[1];
					node->array[node->size - 1].p.z= xyz[2];
				}else if(idc.i == 128){	 /* put element */
					(element->size)++;
					RC_TRY( realloc_fem_element_array(element) );
					element->array[element->size - 1].label = (element->size);
					element->array[element->size - 1].type = ELEM_TRI1;
					element->array[element->size - 1].node_num = 3;
				} else {
					fprintf(stderr, "FORMAT_ERROR < \n");
					return(IMPLEMENT_ERROR_RC);
				}
				break;
			case 10:	/* x */
			case 20:	/* y */
			case 30:	/* z */
				RC_TRY( fgetidc(fp, "d", &idc) );
				xyz[num / 10 - 1] = idc.d;
				break;
			case 71:	/* element_No... */
			case 72:
			case 73:
				RC_TRY( fgetidc(fp, "i", &idc) );
				element->array[element->size - 1].node[num % 70 - 1] = idc.i;
				break;
			case 0:
				if(fseek(fp, file_position, SEEK_SET) != 0){
					fprintf(stderr, "fseek < ");
					return(SEEK_ERROR_RC);
				}
				flag = 1;
				break;
			default:
				break;
			}
			if (flag == 1) break;
		}
	}
	RC_TRY( clean_fem_element_array(element) );
	RC_TRY( clean_fem_node_array(node) );

	return(NORMAL_RC);
}

static void clean_str(char *str)
{
	int rp = 0;
	int wp = 0;

	while( (str[rp] != (char)'\0')&&(str[rp] != (char)'\n') ){
		if( (wp != 0)||(str[rp] != (char)' ') ){
			str[wp] = str[rp];
			wp++;
		}
		rp++;
	}
	str[wp] = (char)'\0';
}


static int seek_0(FILE *fp)
{
	char buf[BUFFER_SIZE];
	RC rc;
	IDC idc;

	idc.i = 1;
	while(1){
		if(fgets(buf, sizeof(buf), fp) == NULL){
			return(-1);
		}
		rc = sgetidc(buf, "i", &idc);
		if((rc == NORMAL_RC)&&(idc.i == 0)) break;
	}

	if(fgets(buf, sizeof(buf), fp) == NULL){
		return(-1);
	}
	rc = sgetidc(buf, "s", &idc);
	if(rc != NORMAL_RC) return(0);
	
	if(strncmp(buf,dxf_section,strlen(dxf_section)) == 0)return(1);
	if(strncmp(buf,dxf_3dface,strlen(dxf_3dface)) == 0)return(2);
	if(strncmp(buf,dxf_vertex,strlen(dxf_vertex)) == 0)return(3);
	if(strncmp(buf,dxf_polyline,strlen(dxf_polyline)) == 0)return(4);
	if(strncmp(buf,dxf_eof,strlen(dxf_eof)) == 0)return(-2);

	return(0);
}

