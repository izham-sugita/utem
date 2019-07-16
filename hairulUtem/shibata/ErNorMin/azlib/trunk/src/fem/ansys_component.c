/*********************************************************************
 * ansys_component.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Tomoyuki TSUBATA> <Syoji ITO>
 *
 * Thanks to following contributors.
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: ansys_component.c 842 2006-12-04 06:21:20Z sasaoka $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ansys_component.h"
#include "fem_struct.h"
#include "rc.h"
#include "string_utl.h"

static RC input_element_node(int *node, int num, ELEMENT *elem);
static void ansys_num_bc(int *numnod, int labnod, int num_bc, BC *bc);
static void ansys_output_cdb_file_end( FILE *fp);


RC
ansys_input_node (FILE *fp, NODE_ARRAY *node)
{
	IDC idc[9];

	RC_TRY(allocate_node_array(0, node));
	node->source = FEM_ANSYS;

	rewind(fp);
	while(1){
		if( (fgetidc(fp,"sssiiiddd", idc) ) == NORMAL_RC){
			if(strncmp(idc[0].s, "N", 1) == 0){
				RC_TRY( realloc_node_array(node) );
				node->array[node->size].label     = idc[3].i; /* node label */
				node->array[node->size].i_info[0] = idc[4].i;
				node->array[node->size].i_info[1] = idc[5].i;
				node->array[node->size].p.x       = idc[6].d; /* x */
				node->array[node->size].p.y       = idc[7].d; /* y */
				node->array[node->size].p.z       = idc[8].d; /* z */
				(node->size)++;
			}
		} else if( feof(fp) ){
			break;
		}
	}
	RC_TRY( clean_node_array(node) );

	return(NORMAL_RC);
}


RC
ansys_output_node (FILE *fp, NODE_ARRAY node)
{
	int ii1;
	
	RC_NULL_CHK(node.array);
	
	for(ii1=0; ii1<node.size; ii1++) {
		if(node.array[ii1].label < 0) continue;
		if(node.source != FEM_ANSYS){
			fprintf( fp, "N,R5.3,LOC, %7d,%7d,%7d, %9.8E , %9.8E , %9.8E\n", 
			              node.array[ii1].label, node.array[ii1].i_info[0], 
			              node.array[ii1].i_info[1], node.array[ii1].p.x,
		                  node.array[ii1].p.y, node.array[ii1].p.z);
		}else{
			fprintf( fp, "N,R5.3,LOC, %7d,%7d,%7d, %9.8E , %9.8E , %9.8E\n", 
			              node.array[ii1].label, 0, 0, node.array[ii1].p.x,
		                  node.array[ii1].p.y, node.array[ii1].p.z);
		}
	}

	return(NORMAL_RC);
}

	
static RC
input_element_node (int *node, int num, ELEMENT *elem)
{
	int	ii1, et_num, type_num;

	et_num   = elem[num].i_info[0];
	type_num = elem[num].i_info[1];

	if( (type_num == 42) || (type_num == 141) ) {
		if(node[2] == node[3]) {
			elem[num].type = ELEM_TRI1;
			elem[num].node_num = 3;
			for(ii1=0; ii1<elem[num].node_num; ii1++) {
				elem[num].node[ii1] = node[ii1];
			}
		}else{
			elem[num].type = ELEM_QUAD1;
			for(ii1=0; ii1<elem[num].node_num; ii1++) {
				elem[num].node[ii1] = node[ii1];
			}
		}
	}else if(type_num == 82) {
		if(node[2] == node[3]) {
			elem[num].type = ELEM_TRI2;
			elem[num].node_num = 6;
			for(ii1=0; ii1<3; ii1++) {
				elem[num].node[ii1] = node[ii1];
			}
			elem[num].node[3] = node[5];
			elem[num].node[4] = node[7];
			elem[num].node[5] = node[4];
		}else{
			elem[num].type = ELEM_QUAD2;
			for(ii1=0; ii1<elem[num].node_num; ii1++) {
				elem[num].node[ii1] = node[ii1];
			}
		}
	}else if( (type_num == 45) || (type_num == 142) ) {
		if(node[2] == node[3]) {
			if(node[4] == node[5]) {
				elem[num].type = ELEM_TETRA1;
				elem[num].node_num = 4;
				elem[num].node[0] = node[0];
				elem[num].node[1] = node[2];
				elem[num].node[2] = node[1];
				elem[num].node[3] = node[4];
			}else{
				elem[num].type = ELEM_PENTA1;
				elem[num].node_num = 5;
				for(ii1=0; ii1<3; ii1++) {
					elem[num].node[ii1] = node[ii1];
				}
				for(ii1 = 4; ii1<7; ii1++) {
					elem[num].node[ii1-1] = node[ii1];
				}
			}
		} else {
			elem[num].type = ELEM_HEXA1;
			for(ii1=0; ii1<elem[num].node_num; ii1++) {
				elem[num].node[ii1] = node[ii1];
			}
		}
	}else if(type_num == 92) {
		elem[num].type = ELEM_TETRA2;
		elem[num].node_num = 10;
		elem[num].node[0] = node[0];
		elem[num].node[1] = node[2];
		elem[num].node[2] = node[1];
		elem[num].node[3] = node[3];
		elem[num].node[4] = node[5];
		elem[num].node[5] = node[4];
		elem[num].node[6] = node[6];
		elem[num].node[7] = node[7];
		elem[num].node[8] = node[9];
		elem[num].node[9] = node[8];
	} else if( type_num == 95 ) {
		elem[num].type = ELEM_HEXA2;
		elem[num].node_num = 20;
		for(ii1=0; ii1<20; ii1++) {
			elem[num].node[ii1] = node[ii1];
		}
	} else {
		fputs("input_element_node : ",stderr);
		return(READ_ERROR_RC);
	}
	return(NORMAL_RC);

}


RC
ansys_input_element (FILE *fp, ELEMENT_ARRAY *elem)
{
	int node, et_num, tmp[20];
	IDC idc[13];

	node = 0;

	RC_TRY( allocate_element_array(0, elem));

	rewind(fp);
	/* ANSYS用要素タイプ識別番号の入力 */
	while(1){
		if( ( (fgetidc(fp,"sii", idc ) ) == NORMAL_RC ) &&
		      (strncmp(idc[0].s, "ET", 2) == 0) ){
			et_num = idc[2].i;
			break;
		}else if(feof(fp)){
			fprintf(stderr,"ansys_input_element ( ET_NUM ) <- \n");
			return(READ_ERROR_RC);
		}
	}

	rewind(fp);
	elem->source = FEM_ANSYS;
	while(1){
		if(  ((fgetidc(fp,"sssiiiiiiiii", idc )) == NORMAL_RC)
		  && (strncmp( idc[2].s, "ATTR", 4 ) == 0)){
			RC_TRY(realloc_element_array(elem));
			/* 要素ラベル */
			elem->array[elem->size].label = idc[9].i;
			/* 材料特性参照番号 */
			elem->array[elem->size].material = idc[4].i;
			/* real constant number */
			elem->array[elem->size].physical = idc[6].i;
			/* 要素内の節点数*/
			elem->array[elem->size].node_num = idc[3].i;
			/* 要素タイプの参照番号 */
			elem->array[elem->size].i_info[0] = idc[5].i;
			/* ANSYS用要素タイプ識別番号 */
			elem->array[elem->size].i_info[1] = et_num;
			node = idc[3].i;

			switch(node){
			case 8:
				if(  (fgetidc(fp,"sssiiiiiiii", idc) == NORMAL_RC)
				  && (strncmp( idc[2].s, "NODE", 4 ) == 0)){
					RC_TRY( realloc_element_array(elem) );
					tmp[0] = idc[3].i;
					tmp[1] = idc[4].i;
					tmp[2] = idc[5].i;
					tmp[3] = idc[6].i;
					tmp[4] = idc[7].i;
					tmp[5] = idc[8].i;
					tmp[6] = idc[9].i;
					tmp[7] = idc[10].i;
					RC_TRY(input_element_node(tmp, elem->size, elem->array));
					(elem->size)++;
				}else{
					return(READ_ERROR_RC);
				}
				break;
			case 4:
				if( ( fgetidc(fp,"sssiiii", idc) ) == NORMAL_RC){
					RC_TRY(realloc_element_array(elem));
					tmp[0] = idc[3].i;
					tmp[1] = idc[4].i;
					tmp[2] = idc[5].i;
					tmp[3] = idc[6].i;
					RC_TRY(input_element_node(tmp, elem->size, elem->array));
					(elem->size)++;
				}else{
					return(READ_ERROR_RC);
				}
				break;
			case 6:
				if( ( fgetidc(fp,"sssiiiiii", idc) ) == NORMAL_RC){
					RC_TRY(realloc_element_array(elem));
					tmp[0] = idc[3].i;
					tmp[1] = idc[4].i;
					tmp[2] = idc[5].i;
					tmp[3] = idc[6].i;
					tmp[4] = idc[7].i;
					tmp[5] = idc[8].i;
					RC_TRY(input_element_node(tmp, elem->size, elem->array));
					(elem->size)++;
				}else{
					return(READ_ERROR_RC);
				}
				break;
			case 10:
				if( ( fgetidc(fp,"sssiiiiiiii", idc) ) == NORMAL_RC){
					RC_TRY( realloc_element_array(elem));
					tmp[0] = idc[3].i;
					tmp[1] = idc[4].i;
					tmp[2] = idc[5].i;
					tmp[3] = idc[6].i;
					tmp[4] = idc[7].i;
					tmp[5] = idc[8].i;
					tmp[6] = idc[9].i;
					tmp[7] = idc[10].i;
					if((fgetidc(fp,"sssii", idc )) == NORMAL_RC){
						RC_TRY(realloc_element_array(elem));
						tmp[8] = idc[3].i;
						tmp[9] = idc[4].i;
					}
					RC_TRY(input_element_node(tmp, elem->size, elem->array));
					(elem->size)++;
				}else{
					return(READ_ERROR_RC);
				}
				break;
			case 20:
				if( ( fgetidc(fp,"sssiiiiiiii", idc) ) == NORMAL_RC){
					RC_TRY(realloc_element_array(elem));
					tmp[0] = idc[3].i;
					tmp[1] = idc[4].i;
					tmp[2] = idc[5].i;
					tmp[3] = idc[6].i;
					tmp[4] = idc[7].i;
					tmp[5] = idc[8].i;
					tmp[6] = idc[9].i;
					tmp[7] = idc[10].i;
					if((fgetidc(fp,"sssiiiiiiii", idc )) == NORMAL_RC){
						RC_TRY(realloc_element_array(elem));
						tmp[8] = idc[3].i;
						tmp[9] = idc[4].i;
						tmp[10] = idc[5].i;
						tmp[11] = idc[6].i;
						tmp[12] = idc[7].i;
						tmp[13] = idc[8].i;
						tmp[14] = idc[9].i;
						tmp[15] = idc[10].i;
						if((fgetidc(fp,"sssiiii", idc )) == NORMAL_RC){
							RC_TRY(realloc_element_array(elem));
							tmp[16] = idc[3].i;
							tmp[17] = idc[4].i;
							tmp[18] = idc[5].i;
							tmp[19] = idc[6].i;
						}
					}
					RC_TRY(input_element_node(tmp, elem->size, elem->array));
					(elem->size)++;
				}else{
					return(READ_ERROR_RC);
				}
				break;
			default:
				fprintf(stderr, "ansys_input_element( node = %4d ):", node);
				return(READ_ERROR_RC);
			}
		}else if( feof(fp) ){
			break;
		}
	}

	return(NORMAL_RC);
}


RC
ansys_output_element (FILE *fp, ELEMENT_ARRAY elem)
{
	int ii1, ii2;
	int node, type;

	node = 0;
	if( elem.array == NULL )
		return(ARG_ERROR_RC);

	for( ii1 = 0 ; ii1 < elem.size ; ii1++ ) {
		if( elem.array[ii1].label < 0 ) continue;

		type = elem.array[ii1].type;
		node = elem.array[ii1].node_num;

		if(type == ELEM_TRI1){
			fprintf(fp, "EN,R5.5,ATTR,      4,%7d,%7d,%7d,%7d,"
			            "      0,%7d,      0, 0, 0\n",
			            elem.array[ii1].material,
			            elem.array[ii1].i_info[0],
			            elem.array[ii1].physical,
			            elem.array[ii1].label,elem.array[ii1].label); 

			fprintf(fp,"EN,R5.5,NODE,");
			for(ii2 = 0 ; ii2<node ; ii2++){
				fprintf(fp, "%7d,",elem.array[ii1].node[ii2]);
			}
			fprintf(fp,"%7d,",elem.array[ii1].node[2]);
		}
		if(type == ELEM_TRI2){
			fprintf(fp, "EN,R5.5,ATTR,      8,%7d,%7d,%7d,%7d,"
			            "      0,%7d,      0, 0, 0\n",
			            elem.array[ii1].material,
			            elem.array[ii1].i_info[0],
			            elem.array[ii1].physical,
			            elem.array[ii1].label,elem.array[ii1].label); 

			fprintf(fp, "EN,R5.5,NODE,");
			fprintf(fp, "%7d,%7d,%7d,%7d,%7d,%7d,%7d,%7d,",
			            elem.array[ii1].node[0],
			            elem.array[ii1].node[1], elem.array[ii1].node[2],
			            elem.array[ii1].node[2], elem.array[ii1].node[5],
			            elem.array[ii1].node[3], elem.array[ii1].node[2],
			            elem.array[ii1].node[4] );
		}
		if((type == ELEM_QUAD1)||(type == ELEM_QUAD2)||(type == ELEM_HEXA1)){
			fprintf(fp, "EN,R5.5,ATTR,%7d,%7d,%7d,%7d,%7d,"
			            "      0,%7d,      0, 0, 0\n",
			            elem.array[ii1].node_num,
			            elem.array[ii1].material,
			            elem.array[ii1].i_info[0],
			            elem.array[ii1].physical,
			            elem.array[ii1].label,
			            elem.array[ii1].label); 
			fprintf(fp, "EN,R5.5,NODE,");
			for( ii2 = 0 ; ii2 < node ; ii2++ ){
				fprintf(fp, "%7d,", elem.array[ii1].node[ii2]);
			}
		}
		if(type == ELEM_TETRA1){
			fprintf(fp, "EN,R5.5,ATTR,      8,%7d,%7d,%7d,%7d,"
			            "      0,%7d,      0, 0, 0\n",
			            elem.array[ii1].material,
			            elem.array[ii1].i_info[0],
			            elem.array[ii1].physical,
			            elem.array[ii1].label,elem.array[ii1].label); 
			fprintf(fp, "EN,R5.5,NODE,");
			fprintf(fp, "%7d,%7d,%7d,%7d,%7d,%7d,%7d,%7d,",
			            elem.array[ii1].node[0],
			            elem.array[ii1].node[2],elem.array[ii1].node[1],
			            elem.array[ii1].node[1],elem.array[ii1].node[3],
			            elem.array[ii1].node[3],elem.array[ii1].node[3],
			            elem.array[ii1].node[3] );
		}
		if(type == ELEM_TETRA2){
			fprintf(fp, "EN,R5.5,ATTR,     20,%7d,%7d,%7d,%7d,"
			            "      0,%7d,      0, 0, 0\n",
			            elem.array[ii1].material,
			            elem.array[ii1].i_info[0],
			            elem.array[ii1].physical,
			            elem.array[ii1].label,elem.array[ii1].label); 

			fprintf(fp, "EN,R5.5,NODE,");
			for( ii2 = 0 ; ii2 < 8 ; ii2++ ){
				fprintf(fp, "%7d,", elem.array[ii1].node[ii2]);
			}
			fprintf(fp, "\n");
			fprintf(fp, "EN,R5.5,NODE,");
			for(ii2=8; ii2<10; ii2++){
				fprintf(fp, "%7d,", elem.array[ii1].node[ii2]);
			}
		}
		if(type == ELEM_HEXA2){
			fprintf(fp, "EN,R5.5,ATTR,%7d,%7d,%7d,%7d,%7d,"
			            "      0,%7d,      0, 0, 0\n",
			            elem.array[ii1].node_num,
			            elem.array[ii1].material,
			            elem.array[ii1].i_info[0],
			            elem.array[ii1].physical,
			            elem.array[ii1].label,
			            elem.array[ii1].label); 

			fprintf(fp, "EN,R5.5,NODE,");
			for( ii2 = 0 ; ii2 < 8 ; ii2++ ){
				fprintf(fp, "%7d,", elem.array[ii1].node[ii2]);
			}
			fprintf(fp, "\n");
			fprintf(fp, "EN,R5.5,NODE,");
			for( ii2 = 8 ; ii2 < 16 ; ii2++ ){
				fprintf(fp, "%7d,", elem.array[ii1].node[ii2]);
			}
			fprintf(fp, "\n");
			fprintf(fp, "EN,R5.5,NODE,");
			for(ii2=16; ii2<20; ii2++){
				fprintf(fp, "%7d,", elem.array[ii1].node[ii2]);
			}
		}
		fprintf(fp, "\n");
	}
	return(NORMAL_RC);
}


static void
ansys_num_bc (int *numnod, int labnod, int num_bc, BC *bc)
{
	int ii1;

	(*numnod) = -1;
	for(ii1 = 0 ; ii1 < num_bc ; ii1++ ) {
		if(labnod == bc[ii1].node) {
			(*numnod) = ii1;
			break;
		}
	}
}


RC
ansys_input_force (FILE *fp, BC_ARRAY *bc)
{
	int numnod;
	IDC tmp[5];
	
	RC_TRY( allocate_bc_array(0, bc) );
	bc->source = FEM_ANSYS;

	rewind(fp);
	while(1) {
		if((fgetidc(fp,"sisdd", tmp )) == NORMAL_RC){
			if(strncmp( tmp[0].s, "F", 1 ) == 0){
				RC_TRY(realloc_bc_array(bc));
				ansys_num_bc( &(numnod), tmp[1].i, bc->size, bc->array );
				if( numnod < 0 ) {
					bc->array[bc->size].node = tmp[1].i;
					numnod = bc->size;
					(bc->size)++;
				}	
				if( strncmp( tmp[2].s, "FX", 2 ) == 0 )  {
					bc->array[numnod].v_type.x = BC_FORCE;
					bc->array[numnod].v.x = tmp[3].d;
				} else if( strncmp( tmp[2].s, "FY", 2 ) == 0 ) {
					bc->array[numnod].v_type.y = BC_FORCE;
					bc->array[numnod].v.y = tmp[3].d;
				} else if( strncmp( tmp[2].s, "FZ", 2 ) == 0 ) {
					bc->array[numnod].v_type.z = BC_FORCE;
					bc->array[numnod].v.z = tmp[3].d;
				} else {
					return(READ_ERROR_RC);
				}
			}
		} else if(feof(fp)){
			break;
		}
	}
	RC_TRY(clean_bc_array(bc));
	return(NORMAL_RC);
}


RC
ansys_output_force (FILE *fp, BC_ARRAY bc)
{
	int ii1;

	for( ii1 = 0 ; ii1 < bc.size ; ii1++ ) {
		if( bc.array[ii1].v_type.x == BC_FORCE ) {
			fprintf(fp,"F, %6d,FX  , %10.9E    , %10.9E \n",
			        bc.array[ii1].node, bc.array[ii1].v.x, 0.0);
		} 
		if( bc.array[ii1].v_type.y == BC_FORCE ) {
			fprintf(fp,"F, %6d,FY  , %10.9E    , %10.9E \n",
			        bc.array[ii1].node, bc.array[ii1].v.y, 0.0);
		}
		if( bc.array[ii1].v_type.z == BC_FORCE ) {
			fprintf(fp,"F, %6d,FZ  , %10.9E    , %10.9E \n",
			        bc.array[ii1].node, bc.array[ii1].v.z, 0.0);
		}
	}
	return(NORMAL_RC);
}


RC
ansys_input_bc_displacement (FILE *fp, BC_ARRAY *bc)
{
	int     numnod;
	IDC     tmp[5];
	
	RC_TRY( allocate_bc_array(0, bc) );
	bc->source = FEM_ANSYS;

	rewind(fp);
	while(1) {
		if((fgetidc(fp,"sisd", tmp )) == NORMAL_RC){
			if((strncmp( tmp[0].s, "D", 1 ) == 0) && (tmp[1].i != 0) ){
				RC_TRY( realloc_bc_array(bc));
				ansys_num_bc( &(numnod), tmp[1].i, bc->size, bc->array );
				if( numnod < 0 ) {
					bc->array[bc->size].node = tmp[1].i;
					numnod = bc->size;
					(bc->size)++;
				}
				if(strncmp( tmp[2].s, "UX", 2 ) == 0){
					bc->array[numnod].v_type.x  = BC_FIX;
					bc->array[numnod].v.x       = tmp[3].d;
				} else if(strncmp( tmp[2].s, "UY", 2 ) == 0){
					bc->array[numnod].v_type.y  = BC_FIX;
					bc->array[numnod].v.y       = tmp[3].d;
				} else if(strncmp( tmp[2].s, "UZ", 2 ) == 0){
					bc->array[numnod].v_type.z  = BC_FIX;
					bc->array[numnod].v.z       = tmp[3].d;
				} else if(strncmp( tmp[2].s, "ROTX", 4 ) == 0){
					bc->array[numnod].v_type.yz = BC_FIX;
					bc->array[numnod].v.yz      = tmp[3].d;
				} else if(strncmp( tmp[2].s, "ROTY", 4 ) == 0){
					bc->array[numnod].v_type.zx = BC_FIX;
					bc->array[numnod].v.zx      = tmp[3].d;
				} else if(strncmp( tmp[2].s, "ROTZ", 4 ) == 0){
					bc->array[numnod].v_type.xy = BC_FIX;
					bc->array[numnod].v.xy      = tmp[3].d;
				} else if(strncmp( tmp[2].s, "VX", 2 ) == 0){
					bc->array[numnod].v_type.x  = BC_VELO;
					bc->array[numnod].v.x       = tmp[3].d;
				} else if(strncmp( tmp[2].s, "VY", 2 ) == 0){
					bc->array[numnod].v_type.y  = BC_VELO;
					bc->array[numnod].v.y       = tmp[3].d;
				} else if(strncmp( tmp[2].s, "VZ", 2 ) == 0){
					bc->array[numnod].v_type.z  = BC_VELO;
					bc->array[numnod].v.z       = tmp[3].d;
				} else if(strncmp( tmp[2].s, "PRES", 4 ) == 0){
					bc->array[numnod].s_type[0] = BC_PRESS;
					bc->array[numnod].s[0]      = tmp[3].d;
				} else {
				}
			}
		} else if(feof(fp)){
			break;
		} 
	}
	RC_TRY( clean_bc_array(bc));

	return(NORMAL_RC);
}


RC
ansys_input_bc_temperature (FILE *fp, BC_ARRAY *bc)
{
	int     numnod;
	IDC     tmp[5];
	
	RC_TRY( allocate_bc_array(0, bc) );
	bc->source = FEM_ANSYS;

	rewind(fp);
	while(1) {
		if((fgetidc(fp,"sisd", tmp )) == NORMAL_RC){
			if((strncmp( tmp[0].s, "BF", 2 ) == 0) && (tmp[1].i != 0)){
				RC_TRY( realloc_bc_array(bc));
				ansys_num_bc( &(numnod), tmp[1].i, bc->size, bc->array );
				if( numnod < 0 ) {
					bc->array[bc->size].node = tmp[1].i;
					numnod = bc->size;
					(bc->size)++;
				}
				if(strncmp( tmp[2].s, "TEMP", 4 ) == 0){
					bc->array[numnod].s_type[0] = BC_TEMP;
					bc->array[numnod].s[0] = tmp[3].d;
				}
			}
		} else if(feof(fp)){
			break;
		} 
	}
	RC_TRY( clean_bc_array(bc));

	return(NORMAL_RC);
}


RC
ansys_output_bc (FILE *fp, BC_ARRAY bc)
{
	int ii1;

	for( ii1 = 0 ; ii1 < bc.size ; ii1++ ) {
		if( bc.array[ii1].v_type.x == BC_FIX ) {
			fprintf( fp,"D, %6d,UX  , %10.9E    , %10.9E \n",
					bc.array[ii1].node, bc.array[ii1].v.x, 0.0 );
		} 
		if( bc.array[ii1].v_type.y == BC_FIX ) {
			fprintf( fp,"D, %6d,UY  , %10.9E    , %10.9E \n",
					bc.array[ii1].node, bc.array[ii1].v.y, 0.0 );
		}
		if( bc.array[ii1].v_type.z == BC_FIX ) {
			fprintf( fp,"D, %6d,UZ  , %10.9E    , %10.9E \n",
					bc.array[ii1].node, bc.array[ii1].v.z, 0.0 );
		}
		if( bc.array[ii1].v_type.x == BC_VELO ) {
			fprintf( fp,"D, %6d,VX  , %10.9E    , %10.9E \n",
					bc.array[ii1].node, bc.array[ii1].v.x, 0.0 );
		}
		if( bc.array[ii1].v_type.y == BC_VELO ) {
			fprintf( fp,"D, %6d,VY  , %10.9E    , %10.9E \n",
					bc.array[ii1].node, bc.array[ii1].v.y, 0.0 );
		}
		if( bc.array[ii1].v_type.z == BC_VELO ) {
			fprintf( fp,"D, %6d,VZ  , %10.9E    , %10.9E \n",
					bc.array[ii1].node, bc.array[ii1].v.z, 0.0 );
		}
		if( bc.array[ii1].s_type[0] == BC_PRESS ) {
			fprintf( fp,"D, %6d,PRES, %10.9E    , %10.9E \n",
					bc.array[ii1].node, bc.array[ii1].s[0], 0.0 );
		}
	}
	return(NORMAL_RC);
}

RC
ansys_input_restraint (FILE *fp, BC_ARRAY *rest)
{
	int     numnod;
	IDC     tmp[5];
	
	RC_TRY( allocate_bc_array(0, rest) );
	rest->source = FEM_ANSYS;

	rewind(fp);
	while(1) {
		if((fgetidc(fp,"sisd", tmp )) == NORMAL_RC){
			if((strncmp( tmp[0].s, "D", 1 ) == 0) && (tmp[1].i != 0) ){
				RC_TRY( realloc_bc_array(rest));
				ansys_num_bc( &(numnod), tmp[1].i, rest->size, rest->array );
				if( numnod < 0 ) {
					rest->array[rest->size].node = tmp[1].i;
					numnod = rest->size;
					(rest->size)++;
				}
				if(strncmp( tmp[2].s, "UX", 2 ) == 0){
					rest->array[numnod].v_type.x  = BC_FIX;
					rest->array[numnod].v.x       = tmp[3].d;
				} else if(strncmp( tmp[2].s, "UY", 2 ) == 0){
					rest->array[numnod].v_type.y  = BC_FIX;
					rest->array[numnod].v.y       = tmp[3].d;
				} else if(strncmp( tmp[2].s, "UZ", 2 ) == 0){
					rest->array[numnod].v_type.z  = BC_FIX;
					rest->array[numnod].v.z       = tmp[3].d;
				} else if(strncmp( tmp[2].s, "ROTX", 4 ) == 0){
					rest->array[numnod].v_type.yz = BC_FIX;
					rest->array[numnod].v.yz      = tmp[3].d;
				} else if(strncmp( tmp[2].s, "ROTY", 4 ) == 0){
					rest->array[numnod].v_type.zx = BC_FIX;
					rest->array[numnod].v.zx      = tmp[3].d;
				} else if(strncmp( tmp[2].s, "ROTZ", 4 ) == 0){
					rest->array[numnod].v_type.xy = BC_FIX;
					rest->array[numnod].v.xy      = tmp[3].d;
				} else {
				}
			}
		} else if(feof(fp)){
			break;
		} 
	}
	RC_TRY( clean_bc_array(rest));
	return(NORMAL_RC);
}


RC
ansys_input_material (FILE *fp, MATERIAL_PROP_ARRAY *material)
{
	int counter, flag;
	IDC tmp[10];

	counter = 0;
	flag = 0;

	RC_TRY( allocate_material_prop_array(0, material) );
	material->source = FEM_ANSYS;

	rewind(fp);
	while(1) {
		if(((fgetidc(fp,"ssisiid",tmp)) == NORMAL_RC) && 
			(strncmp(tmp[0].s, "MPDATA", 6) == 0)){
			RC_TRY(realloc_material_prop_array(material));
			if( flag == 0 ) {
				init_material_prop( &(material->array[counter]) );
				material->array[counter].label = tmp[4].i;
				flag++;
			}else if(   (flag != 0)
			         && (material->array[counter].label != tmp[4].i) ) {
				counter++;
				init_material_prop(&(material->array[counter]));
				material->array[counter].label = tmp[4].i;
			}
			if( strncmp( tmp[3].s, "EX", 2 ) == 0 ) {
				material->array[counter].E = tmp[6].d;
				material->array[counter].E_vect.x = tmp[6].d;
			} else if( strncmp( tmp[3].s, "EY", 2 ) == 0 ) {
				material->array[counter].E_vect.y = tmp[6].d;
			} else if( strncmp( tmp[3].s, "EZ", 2 ) == 0 ) {
				material->array[counter].E_vect.z = tmp[6].d;
			} else if( strncmp( tmp[3].s, "NUXY", 4 ) == 0 ) {
				material->array[counter].nu = tmp[6].d;
				material->array[counter].nu_vect.xy = tmp[6].d;
			} else if( strncmp( tmp[3].s, "NUYZ", 4 ) == 0 ) {
				material->array[counter].nu_vect.yz = tmp[6].d;
			} else if( strncmp( tmp[3].s, "NUXZ", 4 ) == 0 ) {
				material->array[counter].nu_vect.zx = tmp[6].d;
			} else if( strncmp( tmp[3].s, "DENS", 4) == 0 ) {
				material->array[counter].rho = tmp[6].d;
			} else if( strncmp( tmp[3].s, "GXY", 3 ) == 0 ) {
				material->array[counter].G_vect.xy = tmp[6].d;
			} else {
				return(READ_ERROR_RC);
			}
			material->size = counter + 1; 
		} else if(feof(fp)){
			break;
		}
	}
	RC_TRY(clean_material_prop_array(material));
	return(NORMAL_RC);
}


RC
ansys_output_material (FILE *fp, MATERIAL_PROP_ARRAY material)
{
	int	ii1;

	for(ii1=0; ii1<material.size; ii1++) {
		if(material.array[ii1].E != FANSYSDEFAULT) {
			fprintf(fp, "MPTEMP,R5.0, 1, 1, 0.000000000E+00, \n");
			fprintf(fp, "MPDATA,R5.0, 1,EX  , %5d, 1,  %9.6E    , \n", 
			             ii1+1, material.array[ii1].E);
		}
		if( material.array[ii1].nu != FANSYSDEFAULT ) {
			fprintf(fp, "MPTEMP,R5.0, 1, 1, 0.000000000E+00, \n");
			fprintf(fp, "MPDATA,R5.0, 1,NUXY, %5d, 1,  %9.6E    , \n", 
			             ii1+1, material.array[ii1].nu);
		}
		if( material.array[ii1].rho != FANSYSDEFAULT ) {
			fprintf(fp, "MPTEMP,R5.0, 1, 1, 0.000000000E+00, \n");
			fprintf(fp, "MPDATA,R5.0, 1,DENS, %5d, 1,  %9.6E    , \n", 
			            ii1+1, material.array[ii1].rho);
		}
	}
	return(NORMAL_RC);
}


RC
ansys_output_cdb_file (FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY element,
                       MATERIAL_PROP_ARRAY material,
                       BC_ARRAY bc, BC_ARRAY force)
{
	fprintf(fp,
	    "/COM,ANSYS RELEASE  5.6 \n"
	    "/PREP7 \n"
	    "/NOPR \n"
	    "/TITLE, \n"
	    "*IF,_CDRDOFF,EQ,1,THEN     !if solid model was read in \n"
	    "_CDRDOFF=             !reset flag, numoffs already performed \n"
	    "*ELSE              !offset database for the following FE model \n");


	if(node.size > 0) {
		fprintf(fp, "NUMOFF,NODE, %5d \n", node.size);    
	}
	if(element.size > 0) {
		fprintf(fp, "NUMOFF,ELEM, %5d \n", element.size);
	} 
	if(material.size > 0) {   
		fprintf(fp, "NUMOFF,MAT , %5d \n", material.size);    
	}
	if(element.array[0].i_info[0] > 0) {
		fprintf(fp, "NUMOFF,TYPE, %5d \n", element.array[0].i_info[0]);
	}    
	fprintf(fp, "*ENDIF \n");
	fprintf(fp, "ET, %5d,%3d \n", element.array[0].i_info[0],
	                              element.array[0].i_info[1]);
	RC_TRY(ansys_output_node(fp,node));
	RC_TRY(ansys_output_element(fp,element));
	RC_TRY(ansys_output_material(fp,material));
	RC_TRY(ansys_output_bc(fp,bc));
	RC_TRY(ansys_output_force(fp,force));
	ansys_output_cdb_file_end( fp );

	return(NORMAL_RC);
}


RC
ansys_output_model_cdb (FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY element,
                        BC_ARRAY bc)
{
	fprintf(fp,
	    "/COM,ANSYS RELEASE  5.6 \n"
	    "/PREP7 \n"
	    "/NOPR \n"
	    "/TITLE, \n"
	    "*IF,_CDRDOFF,EQ,1,THEN     !if solid model was read in \n"
	    "_CDRDOFF=             !reset flag, numoffs already performed \n"
	    "*ELSE              !offset database for the following FE model \n");

	if(node.size > 0) {
		fprintf(fp, "NUMOFF,NODE, %5d \n", node.size);    
	}
	if(element.size > 0) {
		fprintf(fp, "NUMOFF,ELEM, %5d \n", element.size );
	} 
	if(element.array[0].i_info[0] > 0) {
		fprintf(fp, "NUMOFF,TYPE, %5d \n", element.array[0].i_info[0]);
	}    
	fprintf(fp, "*ENDIF \n");
	fprintf(fp, "ET, %5d,%3d \n", element.array[0].i_info[0],
	                              element.array[0].i_info[1]);
	RC_TRY(ansys_output_node(fp,node));
	RC_TRY(ansys_output_element(fp,element));
	RC_TRY(ansys_output_bc(fp,bc));
	ansys_output_cdb_file_end( fp );

	return(NORMAL_RC);
}


static void
ansys_output_cdb_file_end (FILE *fp)
{
	fprintf( fp, "/GO \n" );
	fprintf( fp, "FINISH \n" );
}


/* 解析結果ファイルの読み込み */
RC
ansys_read_vp_result (FILE *fp_velo, FILE *fp_press, VELO_ARRAY *vp)
{
	int ii1 = 0;
	int vn, pn;
	double vx, vy, vz, vsum, p;
	char buffer[512];

	RC_TRY(allocate_velo_array(0, vp));

	/* 流速結果ファイルの読み込み */
	while(1){
		if( fgets( buffer, sizeof(buffer) ,fp_velo ) == NULL ) break;
		if( sscanf(buffer, "%d %le %le %le %le",
			               &vn, &vx, &vy, &vz, &vsum) == 5 ) {
			RC_TRY( realloc_velo_array(vp));
			vp->array[vp->size].node = vn;
			vp->array[vp->size].v.x  = vx;
			vp->array[vp->size].v.y  = vy;
			vp->array[vp->size].v.z  = vz;
			(vp->size)++;
		}else if(feof(fp_velo)){
			break;
		}
	}

	/* 圧力結果ファイルの読み込み */
	while(1){
		if( fgets( buffer, sizeof(buffer) ,fp_press ) == NULL )
			break;
		if( sscanf(buffer,"%d %le ", &pn, &p) == 2 ) {
			vp->array[ii1].s = p;
			ii1++;
		}else if(feof(fp_press)){
			break;
		}
	}	
/*
	for(ii1=0;ii1<vp->size;ii1++){
		if(((fgetidc(fp_press,"id", press )) == NORMAL_RC)&&
				(vp->array[ii1].node = press[0].i)){
			vp->array[ii1].s = press[1].d;
			fprintf(stderr, "%f \n", vp->array[ii1].s);
		}else if(feof(fp_press)){
			break;
		}
	}
*/
	RC_TRY(clean_velo_array(vp));
 
	return(NORMAL_RC);
}


RC
ansys_input_bc (FILE *fp, BC_ARRAY *bc)
{
	int     numnod;
	IDC     tmp[5];
	
	RC_TRY( allocate_bc_array(0, bc) );
	bc->source = FEM_ANSYS;

	rewind(fp);
	while(1) {
		if((fgetidc(fp,"sisd", tmp )) == NORMAL_RC){
			if((strncmp( tmp[0].s, "D", 1 ) == 0) && (tmp[1].i != 0) ){
				RC_TRY( realloc_bc_array(bc));
				ansys_num_bc( &(numnod), tmp[1].i, bc->size, bc->array );
				if( numnod < 0 ) {
					bc->array[bc->size].node = tmp[1].i;
					numnod = bc->size;
					(bc->size)++;
				}
				if(strncmp( tmp[2].s, "UX", 2 ) == 0){
					bc->array[numnod].v_type.x  = BC_FIX;
					bc->array[numnod].v.x       = tmp[3].d;
				} else if(strncmp( tmp[2].s, "UY", 2 ) == 0){
					bc->array[numnod].v_type.y  = BC_FIX;
					bc->array[numnod].v.y       = tmp[3].d;
				} else if(strncmp( tmp[2].s, "UZ", 2 ) == 0){
					bc->array[numnod].v_type.z  = BC_FIX;
					bc->array[numnod].v.z       = tmp[3].d;
				} else if(strncmp( tmp[2].s, "ROTX", 4 ) == 0){
					bc->array[numnod].v_type.yz = BC_FIX;
					bc->array[numnod].v.yz      = tmp[3].d;
				} else if(strncmp( tmp[2].s, "ROTY", 4 ) == 0){
					bc->array[numnod].v_type.zx = BC_FIX;
					bc->array[numnod].v.zx      = tmp[3].d;
				} else if(strncmp( tmp[2].s, "ROTZ", 4 ) == 0){
					bc->array[numnod].v_type.xy = BC_FIX;
					bc->array[numnod].v.xy      = tmp[3].d;
				} else if(strncmp( tmp[2].s, "VX", 2 ) == 0){
					bc->array[numnod].v_type.x  = BC_VELO;
					bc->array[numnod].v.x       = tmp[3].d;
				} else if(strncmp( tmp[2].s, "VY", 2 ) == 0){
					bc->array[numnod].v_type.y  = BC_VELO;
					bc->array[numnod].v.y       = tmp[3].d;
				} else if(strncmp( tmp[2].s, "VZ", 2 ) == 0){
					bc->array[numnod].v_type.z  = BC_VELO;
					bc->array[numnod].v.z       = tmp[3].d;
				} else if(strncmp( tmp[2].s, "PRES", 4 ) == 0){
					bc->array[numnod].s_type[0] = BC_PRESS;
					bc->array[numnod].s[0]      = tmp[3].d;
				} else {
				}
			}else if((strncmp( tmp[0].s, "BF", 2 ) == 0) && (tmp[1].i != 0)){
				RC_TRY( realloc_bc_array(bc));
				ansys_num_bc( &(numnod), tmp[1].i, bc->size, bc->array );
				if( numnod < 0 ) {
					bc->array[bc->size].node = tmp[1].i;
					numnod = bc->size;
					(bc->size)++;
				}
				if(strncmp( tmp[2].s, "TEMP", 4 ) == 0){
					bc->array[numnod].s_type[0] = BC_TEMP;
					bc->array[numnod].s[0] = tmp[3].d;
				}
			}
		} else if(feof(fp)){
			break;
		} 
	}
	RC_TRY( clean_bc_array(bc));
	return(NORMAL_RC);
}

