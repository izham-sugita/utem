/*********************************************************************
 * nst_blk_component.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> 
 *
 * Thanks to following contributors.
 *  <Takaaki NAGATANI> <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: nst_blk_component.c 1093 2016-12-02 08:24:58Z okada $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nst_component.h"
#include "fem_struct.h"
#include "base.h"
#include "mathematics.h"


#define TOKS_SIZE       (512)
#define TOKS_SIZE_LARGE (4096)

#define INIT_SIZE       (1024)
#define INIT_SIZE_SMALL (4)

#define FORCE_FLAG   (1)
#define MOMENT_FLAG  (2)

/* local functions */
static RC set_grid(NODE_ARRAY *node, char **toks, int num_toks);
static RC set_ctria3(ELEMENT_ARRAY *element, char **toks, int num_toks);
static RC set_ctria6(ELEMENT_ARRAY *element, char **toks, int num_toks);
static RC set_cquad4(ELEMENT_ARRAY *element, char **toks, int num_toks);
static RC set_cquad8(ELEMENT_ARRAY *element, char **toks, int num_toks);
static RC set_ctetra(ELEMENT_ARRAY *element, char **toks, int num_toks);
static RC set_cpenta(ELEMENT_ARRAY *element, char **toks, int num_toks);
static RC set_chexa(ELEMENT_ARRAY *element, char **toks, int num_toks);
static RC set_cbar(ELEMENT_ARRAY *element, char **toks, int num_toks);
static RC set_cbeam(ELEMENT_ARRAY *element, char **toks, int num_toks);
static RC set_mat1(MATERIAL_PROP_ARRAY *material, char **toks, int num_toks);
static RC set_mat2(MATERIAL_PROP_ARRAY *material, char **toks, int num_toks);
static RC set_mat4(MATERIAL_PROP_ARRAY *material, char **toks, int num_toks);
static RC set_mat5(MATERIAL_PROP_ARRAY *material, char **toks, int num_toks);
static RC set_mat8(MATERIAL_PROP_ARRAY *material, char **toks, int num_toks);
static RC set_mat9(MATERIAL_PROP_ARRAY *material, char **toks, int num_toks);
static RC set_mat10(MATERIAL_PROP_ARRAY *material, char **toks, int num_toks);
static RC set_pshell(PHYSICAL_PROP_ARRAY *phys, char **toks, int num_toks);
static RC set_psolid(PHYSICAL_PROP_ARRAY *phys, char **toks, int num_toks);
static RC set_pbar(PHYSICAL_PROP_ARRAY *phys, char **toks, int num_toks);
static RC set_pbarl(PHYSICAL_PROP_ARRAY *phys, char **toks, int num_toks);
static RC set_pelas(PHYSICAL_PROP_ARRAY *phys, char **toks, int num_toks);
static RC set_pbeam (PHYSICAL_PROP_ARRAY *phys, char **toks, int num_toks);
static RC set_force(BC_ARRAY *force, char **toks, int num_toks);
static RC set_moment(BC_ARRAY *force, char **toks, int num_toks);
static RC set_grav(DEFAULT_BC_ARRAY *b_force, char **toks, int num_toks);
static RC set_sload(BC_ARRAY *force, char **toks, int num_toks);
static RC set_temp(BC_ARRAY *temp, char **toks, int num_toks);
static RC set_tempd(DEFAULT_BC_ARRAY *def_temp, char **toks, int num_toks);
static RC set_qload(BC_ARRAY *flux, char **toks, int num_toks);
static RC set_spc(BC_ARRAY *rest, char **toks, int num_toks);
static RC set_spc1(BC_ARRAY *rest, char **toks, int num_toks);
/*static RC set_tempp1(ELEM_BC_ARRAY *elem_temp, char **toks, int num_toks);*/
static RC set_spcadd(BC_SET_ARRAY *bc_set, char **toks, int num_toks);
static RC set_mpcadd(BC_SET_ARRAY *bc_set, char **toks, int num_toks);
static RC set_load(BC_SET_ARRAY *bc_set, char **toks, int num_toks);
static RC set_mpc(RIGID_ELEMENT_ARRAY *element, char **toks, int num_toks);
static RC set_rbar(RIGID_ELEMENT_ARRAY *element, char **toks, int num_toks);
static RC set_rbe2(RIGID_ELEMENT_ARRAY *element, char **toks, int num_toks);
static RC set_celas2(RIGID_ELEMENT_ARRAY *element,char **toks,int num_toks);
static RC set_celas1(RIGID_ELEMENT_ARRAY *element,char **toks,int num_toks);
static RC set_cord2r(LOCAL_COORD_ARRAY *coord, char **toks, int num_toks);
static RC set_pload4(ELEM_BC_ARRAY *press, char **toks, int num_toks,
                     ELEMENT_ARRAY element);
static RC set_conm2(RIGID_ELEMENT_ARRAY *element, char **toks,int num_toks);
#if 0
static RC set_rbe3(RIGID_ELEMENT_ARRAY *element,
                   char **toks, int num_toks);
#endif /* 0 */
static RC set_param(NST_PARAM_ARRAY *param, char **toks, int num_toks);

static RC put_ctria3(FILE *fp, ELEMENT element, DATA_TYPE source);
static RC put_ctria6(FILE *fp, ELEMENT element, DATA_TYPE source);
static RC put_cquad4(FILE *fp, ELEMENT element, DATA_TYPE source);
static RC put_cquad8(FILE *fp, ELEMENT element, DATA_TYPE source);
static RC put_ctetra(FILE *fp, ELEMENT element, DATA_TYPE source);
static RC put_cpenta(FILE *fp, ELEMENT element, DATA_TYPE source);
static RC put_chexa(FILE *fp, ELEMENT element, DATA_TYPE source);
static RC put_cbar(FILE *fp, ELEMENT element, DATA_TYPE source);
static RC put_cbeam(FILE *fp, ELEMENT element, DATA_TYPE source);
static RC put_mat1(FILE *fp, MATERIAL_PROP mat, DATA_TYPE source);
static RC put_mat2(FILE *fp, MATERIAL_PROP mat, DATA_TYPE source);
static RC put_mat8(FILE *fp, MATERIAL_PROP mat, DATA_TYPE source);
static RC put_mat9(FILE *fp, MATERIAL_PROP mat, DATA_TYPE source);
static RC put_mat10(FILE *fp, MATERIAL_PROP mat, DATA_TYPE source);
static RC put_pshell(FILE *fp, PHYSICAL_PROP phys, DATA_TYPE source);
static RC put_psolid(FILE *fp, PHYSICAL_PROP phys, DATA_TYPE source);
static RC put_pbar(FILE *fp, PHYSICAL_PROP phys, DATA_TYPE source);
static RC put_pbarl(FILE *fp, PHYSICAL_PROP phys, DATA_TYPE source);
static RC put_pelas(FILE *fp, PHYSICAL_PROP phys, DATA_TYPE source);
static RC put_pbeam(FILE *fp, PHYSICAL_PROP phys, DATA_TYPE source);
static RC put_force(FILE *fp, BC force, DATA_TYPE source);
static RC put_moment(FILE *fp, BC force, DATA_TYPE source);
static RC put_grav(FILE *fp, DEFAULT_BC b_force, DATA_TYPE source);
static RC put_sload (FILE *fp, BC_ARRAY force);
static RC put_spc1(FILE *fp, BC rest);
static RC put_spc_all(FILE *fp, BC rest, double di);
static RC put_spc(FILE *fp, BC rest);
static RC put_temp(FILE *fp, BC temp);
static RC put_tempd(FILE *fp, DEFAULT_BC def_temp);
/*static RC put_tempp1(FILE *fp, ELEM_BC elem_temp);*/
static RC put_spcadd(FILE *fp, BC_SET bc_set);
static RC put_mpcadd(FILE *fp, BC_SET bc_set);
static RC put_load(FILE *fp, BC_SET bc_set);
static RC put_mpc(FILE *fp, RIGID_ELEMENT element, DATA_TYPE source);
static RC put_rbar(FILE *fp, RIGID_ELEMENT element, DATA_TYPE source);
static RC put_celas(FILE *fp, RIGID_ELEMENT element, DATA_TYPE source);
static RC put_cord2r(FILE *fp, LOCAL_COORD coord);
static RC put_pload4(FILE *fp, ELEM_BC press, ELEMENT_ARRAY element,
                     DATA_TYPE source);
static RC put_conm2(FILE *fp, RIGID_ELEMENT element, DATA_TYPE source);
#if 0
static RC put_rbe3(FILE *fp, const RIGID_ELEMENT *element,
                   DATA_TYPE source);
#endif /* 0 */

static BC_TYPE3D nst_bc_type3d(int num);
static int bc_type3d_number(BC_TYPE3D v_type);
static VECT3DR nst_rest_value(BC_TYPE3D type3d, double value);
static SECTION_TYPE section_type(int toks_size, int tok_num,
                                 char **toks, int *dim_num);
static RC nst_read_case_control(FILE *fp, BC_SET_ARRAY *bc_set);
#if 0
static int exist_period(const char *str);
#endif /* 0 */
static int offset_strncmp(const char *buf, const char *key,
                          size_t n, int *offset);
static RC get_case_ctrl(NST_FILE *nst_file, int num_keys,
                        const NST_BULK_KEY *keys,
                        int *hit_key, char **buf_ptr, STRING_ARRAY *stra);
static int chk_spc1(BC rest);
static int chk_spc_all(BC rest, double *di);
static RC ident_face(ELEMENT elem, int g1, int g34,
                     int *face_id, int *tri_flag);
static RC merge_moment(BC_ARRAY *force, BC_ARRAY moment);


RC
nst_input_node (FILE *fp, NODE_ARRAY *node)
{
	const int num_keys = 1;
	NST_BULK_KEY keys[] = {{"GRID", 4}};
	int hit_key;
	char *toks[TOKS_SIZE];
	int num_toks;
	NST_FILE nst_file;

	RC_TRY( nst_allocate_tok(toks, TOKS_SIZE) );
	RC_TRY( init_nst_file(fp, &nst_file) );
	RC_TRY( nst_seek_begin_bulk(&nst_file) );
	RC_TRY( allocate_node_array(INIT_SIZE, node) );
	node->source = FEM_NASTRAN;

	while(1){
		RC rc = nst_blk_get_tok(&nst_file, num_keys, keys, &hit_key,
		                        TOKS_SIZE, toks, &num_toks);
		if(rc != NORMAL_RC){
			if(rc == END_RC) break;
			return(rc);
		}

		RC_TRY( realloc_node_array(node) );

		switch(hit_key){
		case 0:
			RC_TRY( set_grid(node, toks, num_toks) );
			break;
		default:
			return(UNKNOWN_ERROR_RC);
		}
		(node->size)++;
	}

	RC_TRY( clean_node_array(node) );
	RC_TRY( term_nst_file(&nst_file) );
	RC_TRY( nst_free_tok(toks, TOKS_SIZE) );

	return(NORMAL_RC);
}


static RC
set_grid (NODE_ARRAY *node, char **toks, int num_toks)
{
	int itmp;
	double dtmp;

	itmp = nst_tok2i(num_toks, 1, toks);
	if(itmp < 0){
		fprintf(stderr, "tok2i <");
		return(CONVERT_ERROR_RC);
	}
	node->array[node->size].label = itmp; /* ID */

	node->array[node->size].i_info[0]
	        = nst_tok2i(num_toks, 2, toks);  /* CP */

	dtmp = nst_tok2d(num_toks, 3, toks);  /* X1 */
	if(dtmp > CR_FNSDEFAULT){
		node->array[node->size].p.x = dtmp;
	}else{
		node->array[node->size].p.x = 0.0;
	}

	dtmp = nst_tok2d(num_toks, 4, toks);  /* X2 */
	if(dtmp > CR_FNSDEFAULT){
		node->array[node->size].p.y = dtmp;
	}else{
		node->array[node->size].p.y = 0.0;
	}

	dtmp = nst_tok2d(num_toks, 5, toks);  /* X3 */
	if(dtmp > CR_FNSDEFAULT){
		node->array[node->size].p.z = dtmp;
	}else{
		node->array[node->size].p.z = 0.0;
	}

	node->array[node->size].i_info[1]
			= nst_tok2i(num_toks, 6, toks);  /* CD */	

	node->array[node->size].i_info[2]
			= nst_tok2i(num_toks, 7, toks);  /* PS */	

	node->array[node->size].i_info[3]
			= nst_tok2i(num_toks, 8, toks);  /* SEID */	

	return(NORMAL_RC);
}


RC
nst_output_node (FILE *fp, NODE_ARRAY node)
{
	int ii1;
	IDC tmp[10];

	for(ii1=0;ii1<node.size;ii1++){
		if(node.array[ii1].label < 0) continue;

		tmp[0].c = "GRID*";
		tmp[1].i = node.array[ii1].label;     /* ID */
		if(node.source == FEM_NASTRAN){
			tmp[2].i = node.array[ii1].i_info[0];  /* CP */
		}else{
			tmp[2].i = NSDEFAULT;
		}
		tmp[3].d = node.array[ii1].p.x;       /* X */
		tmp[4].d = node.array[ii1].p.y;       /* Y */
		tmp[5].c = "*";
		tmp[6].d = node.array[ii1].p.z;       /* Z */
		if(node.source == FEM_NASTRAN){
			tmp[7].i = node.array[ii1].i_info[1];  /* CD */
			tmp[8].i = node.array[ii1].i_info[2];  /* PS */
			tmp[9].i = node.array[ii1].i_info[3];  /* SEID */
		}else{
			tmp[7].i = NSDEFAULT;
			tmp[8].i = NSDEFAULT;
			tmp[9].i = NSDEFAULT;
		}

		RC_TRY( nst_print(fp, "cIIDDrcDIIIr", tmp) );
	}

	return(NORMAL_RC);
}


RC
nst_input_element (FILE *fp, ELEMENT_ARRAY *element)
{
	const int num_keys = 9;
	NST_BULK_KEY keys[] = {{"CTRIA3", 6},{"CTRIA6", 6},{"CQUAD4", 6},
	                       {"CQUAD8", 6},{"CTETRA", 6},{"CPENTA", 6},
	                       {"CHEXA", 5}, {"CBAR", 4}, {"CBEAM", 5} };
	int hit_key;
	char *toks[TOKS_SIZE];
	int num_toks;
	NST_FILE nst_file;

	RC_TRY( nst_allocate_tok(toks, TOKS_SIZE) );
	RC_TRY( init_nst_file(fp, &nst_file) );
	RC_TRY( nst_seek_begin_bulk(&nst_file) );
	RC_TRY( allocate_element_array(INIT_SIZE, element) );
	element->source = FEM_NASTRAN;

	while(1){
		RC rc = nst_blk_get_tok(&nst_file, num_keys, keys, &hit_key,
		                        TOKS_SIZE, toks, &num_toks);
		if(rc != NORMAL_RC){
			if(rc == END_RC) break;
			return(rc);
		}

		RC_TRY( realloc_element_array(element) );

		switch(hit_key){
		case 0:
			RC_TRY( set_ctria3(element, toks, num_toks) );
			break;
		case 1:
			RC_TRY( set_ctria6(element, toks, num_toks) );
			break;
		case 2:
			RC_TRY( set_cquad4(element, toks, num_toks) );
			break;
		case 3:
			RC_TRY( set_cquad8(element, toks, num_toks) );
			break;
		case 4:
			RC_TRY( set_ctetra(element, toks, num_toks) );
			break;
		case 5:
			RC_TRY( set_cpenta(element, toks, num_toks) );
			break;
		case 6:
			RC_TRY( set_chexa(element, toks, num_toks) );
			break;
		case 7:
			RC_TRY( set_cbar(element, toks, num_toks) );
			break;
		case 8:
			RC_TRY( set_cbeam(element, toks, num_toks) );
			break;
		default:
			return(UNKNOWN_ERROR_RC);
		}
		(element->size)++;
	}

	RC_TRY( clean_element_array(element) );
	RC_TRY( term_nst_file(&nst_file) );
	RC_TRY( nst_free_tok(toks, TOKS_SIZE) );

	return(NORMAL_RC);
}


RC
nst_output_element (FILE *fp, ELEMENT_ARRAY element)
{
	int ii1;

	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		switch(element.array[ii1].type){
		case ELEM_TRI1:
			RC_TRY( put_ctria3(fp, element.array[ii1], element.source) );
			break;
		case ELEM_TRI2:
			RC_TRY( put_ctria6(fp, element.array[ii1], element.source) );
			break;
		case ELEM_QUAD1:
			RC_TRY( put_cquad4(fp, element.array[ii1], element.source) );
			break;
		case ELEM_QUAD2:
			RC_TRY( put_cquad8(fp, element.array[ii1], element.source) );
			break;
		case ELEM_TETRA1:
		case ELEM_TETRA2:
			RC_TRY( put_ctetra(fp, element.array[ii1], element.source) );
			break;
		case ELEM_PENTA1:
		case ELEM_PENTA2:
			RC_TRY( put_cpenta(fp, element.array[ii1], element.source) );
			break;
		case ELEM_HEXA1:
		case ELEM_HEXA2:
			RC_TRY( put_chexa(fp, element.array[ii1], element.source) );
			break;
		case ELEM_BEAM1:
			if(element.array[ii1].i_info[4] == 1){ /* CBEAM */
				RC_TRY( put_cbeam(fp, element.array[ii1], element.source) );
			}else if(element.array[ii1].i_info[4] == -1){ /* CBAR */
				RC_TRY( put_cbar(fp, element.array[ii1], element.source) );
			}else{
				return(ARG_ERROR_RC);
			}
			break;
		default:
			return(IMPLEMENT_ERROR_RC);
		}
	}

	return(NORMAL_RC);
}


static RC
set_ctria3 (ELEMENT_ARRAY *element, char **toks, int num_toks)
{
	int ii1;
	int itmp;

	itmp = nst_tok2i(num_toks, 1, toks);    /* EID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	element->array[element->size].label = itmp;

	element->array[element->size].type = ELEM_TRI1;

	element->array[element->size].material = -1;

	itmp = nst_tok2i(num_toks, 2, toks);  /* PID */	
	if(itmp != NSDEFAULT){
		element->array[element->size].physical = itmp;
	}else{
		element->array[element->size].physical
		       = element->array[element->size].label;
	}

	element->array[element->size].node_num = 3;

	element->array[element->size].node[0]
	        = nst_tok2i(num_toks, 3, toks);  /* G1 */

	element->array[element->size].node[1]
	        = nst_tok2i(num_toks, 4, toks);  /* G2 */

	element->array[element->size].node[2]
	        = nst_tok2i(num_toks, 5, toks);  /* G3 */

	for(ii1=0; ii1<3; ii1++){
		if(element->array[element->size].node[ii1] <= 0){
			return(CONVERT_ERROR_RC);
		}
	}

	element->array[element->size].d_info[0]
			= nst_tok2d(num_toks, 6, toks);  /* THETA */

	element->array[element->size].d_info[1]
			= nst_tok2d(num_toks, 6, toks);  /* ZOFFS */

	element->array[element->size].d_info[2]
			= nst_tok2d(num_toks, 11, toks);  /* T1 */

	element->array[element->size].d_info[3]
			= nst_tok2d(num_toks, 12, toks);  /* T2 */

	element->array[element->size].d_info[4]
			= nst_tok2d(num_toks, 13, toks);  /* T3 */

	return(NORMAL_RC);
}


static RC
put_ctria3 (FILE *fp, ELEMENT element, DATA_TYPE source)
{
	RC rc;
	IDC tmp[8];

	tmp[0].c = "CTRIA3";
	tmp[1].i = element.label;     /* EID */
	tmp[2].i = element.physical;  /* PID */
	if(tmp[2].i <= 0) tmp[2].i = 1;
	tmp[3].i = element.node[0];   /* G1 */
	tmp[4].i = element.node[1];   /* G2 */
	tmp[5].i = element.node[2];   /* G3 */
	if(source == FEM_NASTRAN){
		tmp[6].d = element.d_info[0]; /* THETA */
		tmp[7].d = element.d_info[1]; /* ZOFFS */
	}else{
		tmp[6].d = FNSDEFAULT;
		tmp[7].d = FNSDEFAULT;
	}
	rc = nst_print(fp, "ciiiiiddr", tmp);
	if(rc != NORMAL_RC){
		fprintf(stderr,"nsout < ");
		return(rc);
	}

	if((source != FEM_NASTRAN)||( (element.d_info[2] < CR_FNSDEFAULT)
	                            &&(element.d_info[3] < CR_FNSDEFAULT)
	                            &&(element.d_info[4] < CR_FNSDEFAULT) )){
		return(NORMAL_RC);
	}
	                            
	tmp[0].d = element.d_info[2];
	tmp[1].d = element.d_info[3];
	tmp[2].d = element.d_info[4];
	
	rc = nst_print(fp, "---dddr", tmp);
	if(rc != NORMAL_RC) return(rc);

	return(NORMAL_RC);
}


static RC
set_ctria6 (ELEMENT_ARRAY *element, char **toks, int num_toks)
{
	int ii1;
	int itmp;

	itmp = nst_tok2i(num_toks, 1, toks); /* EID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	element->array[element->size].label = itmp;

	element->array[element->size].type = ELEM_TRI2;

	element->array[element->size].material = -1;

	itmp = nst_tok2i(num_toks, 2, toks);  /* PID */	
	if(itmp < 0) return(CONVERT_ERROR_RC);
	element->array[element->size].physical = itmp;

	element->array[element->size].node_num = 6;

	element->array[element->size].node[0]
	        = nst_tok2i(num_toks, 3, toks);  /* G1 */
	
	element->array[element->size].node[1]
	        = nst_tok2i(num_toks, 4, toks);  /* G2 */
	
	element->array[element->size].node[2]
	        = nst_tok2i(num_toks, 5, toks);  /* G3 */

	element->array[element->size].node[5]
	        = nst_tok2i(num_toks, 6, toks);  /* G4 */
	
	element->array[element->size].node[3]
	        = nst_tok2i(num_toks, 7, toks);  /* G5 */
	
	element->array[element->size].node[4]
	        = nst_tok2i(num_toks, 8, toks);  /* G6 */
	
	for(ii1=0; ii1<6; ii1++){
		if(element->array[element->size].node[ii1] <= 0){
			return(CONVERT_ERROR_RC);
		}
	}

	element->array[element->size].d_info[0]
	        = nst_tok2d(num_toks, 9, toks);  /* THETA */

	element->array[element->size].d_info[1]
	        = nst_tok2d(num_toks, 10, toks);  /* ZOFFS */

	element->array[element->size].d_info[2]
	        = nst_tok2d(num_toks, 11, toks);  /* T1 */

	element->array[element->size].d_info[3]
	        = nst_tok2d(num_toks, 12, toks);  /* T2 */

	element->array[element->size].d_info[4]
	        = nst_tok2d(num_toks, 13, toks);  /* T3 */

	return(NORMAL_RC);
}


static RC
put_ctria6 (FILE *fp, ELEMENT element, DATA_TYPE source)
{
	IDC tmp[9];

	tmp[0].c = "CTRIA6";
	tmp[1].i = element.label;     /* EID */
	tmp[2].i = element.physical;  /* PID */
	if(tmp[2].i <= 0) tmp[2].i = 1;
	tmp[3].i = element.node[0];   /* G1 */
	tmp[4].i = element.node[1];   /* G2 */
	tmp[5].i = element.node[2];   /* G3 */
	tmp[6].i = element.node[5];   /* G4 */
	tmp[7].i = element.node[3];   /* G5 */
	tmp[8].i = element.node[4];   /* G6 */
	RC_TRY( nst_print(fp, "ciiiiiiiir", tmp) );

	if((source != FEM_NASTRAN)||( (element.d_info[0] < CR_FNSDEFAULT)
	                            &&(element.d_info[1] < CR_FNSDEFAULT)
	                            &&(element.d_info[2] < CR_FNSDEFAULT)
	                            &&(element.d_info[3] < CR_FNSDEFAULT)
	                            &&(element.d_info[4] < CR_FNSDEFAULT) )){
		return(NORMAL_RC);
	}

	tmp[0].d = element.d_info[0];   /* THETA */
	tmp[1].d = element.d_info[1];   /* ZOFFS */
	tmp[2].d = element.d_info[2];   /* T1 */
	tmp[3].d = element.d_info[3];   /* T2 */
	tmp[4].d = element.d_info[4];   /* T3 */
	RC_TRY( nst_print(fp, "-dddddr", tmp) );

	return(NORMAL_RC);
}


static RC
set_cquad4 (ELEMENT_ARRAY *element, char **toks, int num_toks)
{
	int ii1;
	int itmp;

	itmp = nst_tok2i(num_toks, 1, toks);     /* EID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	element->array[element->size].label = itmp;
	element->array[element->size].type = ELEM_QUAD1;
	element->array[element->size].material = -1;

	itmp = nst_tok2i(num_toks, 2, toks);     /* PID */
	if(itmp != NSDEFAULT){
		element->array[element->size].physical = itmp;
	}else{
		element->array[element->size].physical
		       = element->array[element->size].label;
	}

	element->array[element->size].node_num = 4;

	element->array[element->size].node[0]
	        = nst_tok2i(num_toks, 3, toks);  /* G1 */
	
	element->array[element->size].node[1]
	        = nst_tok2i(num_toks, 4, toks);  /* G2 */
	
	element->array[element->size].node[2]
	        = nst_tok2i(num_toks, 5, toks);  /* G3 */

	element->array[element->size].node[3]
	        = nst_tok2i(num_toks, 6, toks);  /* G4 */

	for(ii1=0; ii1<4; ii1++){
		if(element->array[element->size].node[ii1] <= 0){
			return(CONVERT_ERROR_RC);
		}
	}

	element->array[element->size].d_info[0]
	        = nst_tok2d(num_toks, 7, toks);  /* THETA */

	element->array[element->size].d_info[1]
	        = nst_tok2d(num_toks, 8, toks);  /* ZOFFS */

	element->array[element->size].d_info[2]
	        = nst_tok2d(num_toks, 11, toks);  /* T1 */

	element->array[element->size].d_info[3]
	        = nst_tok2d(num_toks, 12, toks);  /* T2 */

	element->array[element->size].d_info[4]
	        = nst_tok2d(num_toks, 13, toks);  /* T3 */

	element->array[element->size].d_info[5]
	        = nst_tok2d(num_toks, 14, toks);  /* T4 */

	return(NORMAL_RC);
}


static RC
put_cquad4 (FILE *fp, ELEMENT element, DATA_TYPE source)
{
	IDC tmp[9];

	tmp[0].c = "CQUAD4";
	tmp[1].i = element.label;     /* EID */
	tmp[2].i = element.physical;  /* PID */
	if(tmp[2].i <= 0) tmp[2].i = 1;
	tmp[3].i = element.node[0];   /* G1 */
	tmp[4].i = element.node[1];   /* G2 */
	tmp[5].i = element.node[2];   /* G3 */
	tmp[6].i = element.node[3];   /* G4 */
	if(source == FEM_NASTRAN){
		tmp[7].d = element.d_info[0];   /* THETA */
		tmp[8].d = element.d_info[1];   /* ZOFFS */
	}else{
		tmp[7].d = FNSDEFAULT;
		tmp[8].d = FNSDEFAULT;
	}
	RC_TRY( nst_print(fp, "ciiiiiiddr", tmp) );

	if((source != FEM_NASTRAN)||( (element.d_info[2] < CR_FNSDEFAULT)
	                            &&(element.d_info[3] < CR_FNSDEFAULT)
	                            &&(element.d_info[4] < CR_FNSDEFAULT)
	                            &&(element.d_info[5] < CR_FNSDEFAULT) )){
		return(NORMAL_RC);
	}

	tmp[0].d = element.d_info[2];   /* T1 */
	tmp[1].d = element.d_info[3];   /* T2 */
	tmp[2].d = element.d_info[4];   /* T3 */
	tmp[3].d = element.d_info[5];   /* T4 */
	RC_TRY( nst_print(fp, "---ddddr", tmp) );

	return(NORMAL_RC);
}


static RC
set_cquad8 (ELEMENT_ARRAY *element, char **toks, int num_toks)
{
	int ii1;
	int itmp;

	itmp = nst_tok2i(num_toks, 1, toks);     /* EID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	element->array[element->size].label = itmp;
	element->array[element->size].type = ELEM_QUAD2;
	element->array[element->size].material = -1;

	itmp = nst_tok2i(num_toks, 2, toks);  /* PID */	
	if(itmp < 0) return(CONVERT_ERROR_RC);
	element->array[element->size].physical = itmp;

	element->array[element->size].node_num = 8;

	element->array[element->size].node[0]
	        = nst_tok2i(num_toks, 3, toks);   /* G1 */
	
	element->array[element->size].node[1]
	        = nst_tok2i(num_toks, 4, toks);   /* G2 */
	
	element->array[element->size].node[2]
	        = nst_tok2i(num_toks, 5, toks);   /* G3 */

	element->array[element->size].node[3]
	        = nst_tok2i(num_toks, 6, toks);   /* G4 */
	
	element->array[element->size].node[4]
	        = nst_tok2i(num_toks, 7, toks);   /* G5 */
	
	element->array[element->size].node[5]
	        = nst_tok2i(num_toks, 8, toks);   /* G6 */

	element->array[element->size].node[6]
	        = nst_tok2i(num_toks, 9, toks);   /* G7 */

	element->array[element->size].node[7]
	        = nst_tok2i(num_toks, 10, toks);  /* G8 */

	for(ii1=0; ii1<8; ii1++){
		if(element->array[element->size].node[ii1] <= 0){
			return(CONVERT_ERROR_RC);
		}
	}

	element->array[element->size].d_info[0]
	        = nst_tok2d(num_toks, 15, toks);  /* THETA */

	element->array[element->size].d_info[1]
	        = nst_tok2d(num_toks, 16, toks);  /* ZOFFS */

	element->array[element->size].d_info[2]
	        = nst_tok2d(num_toks, 11, toks);  /* T1 */

	element->array[element->size].d_info[3]
	        = nst_tok2d(num_toks, 12, toks);  /* T2 */

	element->array[element->size].d_info[4]
	        = nst_tok2d(num_toks, 13, toks);  /* T3 */

	element->array[element->size].d_info[5]
	        = nst_tok2d(num_toks, 14, toks);  /* T4 */

	return(NORMAL_RC);
}


static RC
put_cquad8 (FILE *fp, ELEMENT element, DATA_TYPE source)
{
	IDC tmp[9];

	tmp[0].c = "CQUAD8";
	tmp[1].i = element.label;     /* EID */
	tmp[2].i = element.physical;  /* PID */
	if(tmp[2].i <= 0) tmp[2].i = 1;
	tmp[3].i = element.node[0];   /* G1 */
	tmp[4].i = element.node[1];   /* G2 */
	tmp[5].i = element.node[2];   /* G3 */
	tmp[6].i = element.node[3];   /* G4 */
	tmp[7].i = element.node[4];   /* G5 */
	tmp[8].i = element.node[5];   /* G6 */
	RC_TRY( nst_print(fp, "ciiiiiiiir", tmp) );

	tmp[0].i = element.node[6];   /* G7 */
	tmp[1].i = element.node[7];   /* G8 */
	if(source == FEM_NASTRAN){
		tmp[2].d = element.d_info[2];   /* T1 */
		tmp[3].d = element.d_info[3];   /* T2 */
		tmp[4].d = element.d_info[4];   /* T3 */
		tmp[5].d = element.d_info[5];   /* T4 */
		tmp[6].d = element.d_info[0];   /* THETA */
		tmp[7].d = element.d_info[1];   /* ZOFFS */
	}else{
		tmp[2].d = FNSDEFAULT;
		tmp[3].d = FNSDEFAULT;
		tmp[4].d = FNSDEFAULT;
		tmp[5].d = FNSDEFAULT;
		tmp[6].d = FNSDEFAULT;
		tmp[7].d = FNSDEFAULT;
	}

	RC_TRY( nst_print(fp, "-iiddddddr", tmp) );

	return(NORMAL_RC);
}


static RC
set_ctetra (ELEMENT_ARRAY *element, char **toks, int num_toks)
{
	int ii1;
	int itmp;

	itmp = nst_tok2i(num_toks, 1, toks);     /* EID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	element->array[element->size].label = itmp;
	element->array[element->size].type = ELEM_TETRA2;
	element->array[element->size].material = -1;

	itmp = nst_tok2i(num_toks, 2, toks);     /* PID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	element->array[element->size].physical = itmp;
	element->array[element->size].node_num = 10;
	element->array[element->size].node[0]
	        = nst_tok2i(num_toks, 3,  toks);    /* G1 */
	
	element->array[element->size].node[2]
	        = nst_tok2i(num_toks, 4,  toks);    /* G2 */
	
	element->array[element->size].node[1]
	        = nst_tok2i(num_toks, 5,  toks);    /* G3 */

	element->array[element->size].node[3]
	        = nst_tok2i(num_toks, 6,  toks);    /* G4 */

	element->array[element->size].node[5]
	        = nst_tok2i(num_toks, 7,  toks);    /* G5 */
	
	element->array[element->size].node[4]
	        = nst_tok2i(num_toks, 8,  toks);    /* G6 */
	
	element->array[element->size].node[6]
	        = nst_tok2i(num_toks, 9,  toks);    /* G7 */

	element->array[element->size].node[7]
	        = nst_tok2i(num_toks,10,  toks);    /* G8 */

	element->array[element->size].node[9]
	        = nst_tok2i(num_toks,11,  toks);    /* G9 */

	element->array[element->size].node[8]
	        = nst_tok2i(num_toks,12,  toks);    /* G10 */

	for(ii1=0; ii1<4; ii1++){
		if(element->array[element->size].node[ii1] <= 0){
			return(CONVERT_ERROR_RC);
		}
	}

	for(ii1=4; ii1<10; ii1++){
		if(element->array[element->size].node[ii1] <= 0){
			element->array[element->size].node[ii1] = -1;
			element->array[element->size].type = ELEM_TETRA1;
			element->array[element->size].node_num = 4;
		}
	}

	return(NORMAL_RC);
}


static RC
put_ctetra (FILE *fp, ELEMENT element, DATA_TYPE source)
{
	IDC tmp[13];

	tmp[0].c = "CTETRA";
	tmp[1].i = element.label;     /* EID */
	tmp[2].i = element.physical;  /* PID */
	if(tmp[2].i <= 0) tmp[2].i = 1;
	tmp[3].i = element.node[0];   /* G1 */
	tmp[4].i = element.node[2];   /* G2 */
	tmp[5].i = element.node[1];   /* G3 */
	tmp[6].i = element.node[3];   /* G4 */

	if(element.type == ELEM_TETRA1){
		RC_TRY( nst_print(fp, "ciiiiiir", tmp) );
	}else if(element.type == ELEM_TETRA2){
		tmp[7].i = element.node[5];   /* G5 */
		tmp[8].i = element.node[4];   /* G6 */
		tmp[9].i = element.node[6];   /* G7 */
		tmp[10].i = element.node[7];   /* G8 */
		tmp[11].i = element.node[9];   /* G9 */
		tmp[12].i = element.node[8];   /* G10 */
		RC_TRY( nst_print(fp, "ciiiiiiiir-iiiir", tmp) );
	}else{
		return(ARG_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
set_cpenta (ELEMENT_ARRAY *element, char **toks, int num_toks)
{
	int ii1;
	int itmp;

	itmp = nst_tok2i(num_toks, 1, toks);     /* EID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	element->array[element->size].label = itmp;
	element->array[element->size].type = ELEM_PENTA2;
	element->array[element->size].material = -1;

	itmp = nst_tok2i(num_toks, 2, toks);     /* PID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	element->array[element->size].physical = itmp;

	element->array[element->size].node_num = 15;

	element->array[element->size].node[0]
	        = nst_tok2i(num_toks, 3,  toks);    /* G1 */
	
	element->array[element->size].node[1]
	        = nst_tok2i(num_toks, 4,  toks);    /* G2 */
	
	element->array[element->size].node[2]
	        = nst_tok2i(num_toks, 5,  toks);    /* G3 */

	element->array[element->size].node[3]
	        = nst_tok2i(num_toks, 6,  toks);    /* G4 */

	element->array[element->size].node[4]
	        = nst_tok2i(num_toks, 7,  toks);    /* G5 */

	element->array[element->size].node[5]
	        = nst_tok2i(num_toks, 8,  toks);    /* G6 */

	element->array[element->size].node[8]
	        = nst_tok2i(num_toks, 9,  toks);    /* G7 */
	
	element->array[element->size].node[6]
	        = nst_tok2i(num_toks,10,  toks);    /* G8 */
	
	element->array[element->size].node[7]
	        = nst_tok2i(num_toks,11,  toks);    /* G9 */

	element->array[element->size].node[12]
	        = nst_tok2i(num_toks,12,  toks);    /* G10 */

	element->array[element->size].node[13]
	        = nst_tok2i(num_toks,13,  toks);    /* G11 */

	element->array[element->size].node[14]
	        = nst_tok2i(num_toks,14,  toks);    /* G12 */

	element->array[element->size].node[11]
	        = nst_tok2i(num_toks,15,  toks);    /* G13 */

	element->array[element->size].node[9]
	        = nst_tok2i(num_toks,16,  toks);    /* G14 */

	element->array[element->size].node[10]
	        = nst_tok2i(num_toks,17,  toks);    /* G15 */

	for(ii1=0; ii1<6; ii1++){
		if(element->array[element->size].node[ii1] <= 0){
			return(CONVERT_ERROR_RC);
		}
	}

	for(ii1=6; ii1<15; ii1++){
		if(element->array[element->size].node[ii1] <= 0){
			element->array[element->size].node[ii1] = -1;
			element->array[element->size].type = ELEM_PENTA1;
			element->array[element->size].node_num = 6;
		}
	}

	return(NORMAL_RC);
}


static RC
put_cpenta (FILE *fp, ELEMENT element, DATA_TYPE source)
{
	IDC tmp[18];

	tmp[0].c = "CPENTA";
	tmp[1].i = element.label;     /* EID */
	tmp[2].i = element.physical;  /* PID */
	if(tmp[2].i <= 0) tmp[2].i = 1;
	tmp[3].i = element.node[0];    /* G1 */
	tmp[4].i = element.node[1];    /* G2 */
	tmp[5].i = element.node[2];    /* G3 */
	tmp[6].i = element.node[3];    /* G4 */
	tmp[7].i = element.node[4];    /* G5 */
	tmp[8].i = element.node[5];    /* G6 */

	if(element.type == ELEM_PENTA1){
		RC_TRY( nst_print(fp, "ciiiiiiiir", tmp) );
	}else if(element.type == ELEM_PENTA2){
		tmp[9].i = element.node[8];    /* G7 */
		tmp[10].i = element.node[6];   /* G8 */
		tmp[11].i = element.node[7];   /* G9 */
		tmp[12].i = element.node[12];  /* G10 */
		tmp[13].i = element.node[13];  /* G11 */
		tmp[14].i = element.node[14];  /* G12 */
		tmp[15].i = element.node[11];  /* G13 */
		tmp[16].i = element.node[9];   /* G14 */
		tmp[17].i = element.node[10];  /* G15 */
		RC_TRY( nst_print(fp, "ciiiiiiiir-iiiiiiiir-ir", tmp) );
	}else{
		return(ARG_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
set_chexa (ELEMENT_ARRAY *element, char **toks, int num_toks)
{
	int ii1;
	int itmp;

	itmp = nst_tok2i(num_toks, 1, toks);     /* EID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	element->array[element->size].label = itmp;
	element->array[element->size].type = ELEM_HEXA2;
	element->array[element->size].material = -1;

	itmp = nst_tok2i(num_toks, 2, toks);     /* PID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	element->array[element->size].physical = itmp;

	element->array[element->size].node_num = 20;

	element->array[element->size].node[0]
	        = nst_tok2i(num_toks, 3,  toks);    /* G1 */
	
	element->array[element->size].node[1]
	        = nst_tok2i(num_toks, 4,  toks);    /* G2 */
	
	element->array[element->size].node[2]
	        = nst_tok2i(num_toks, 5,  toks);    /* G3 */

	element->array[element->size].node[3]
	        = nst_tok2i(num_toks, 6,  toks);    /* G4 */

	element->array[element->size].node[4]
	        = nst_tok2i(num_toks, 7,  toks);    /* G5 */

	element->array[element->size].node[5]
	        = nst_tok2i(num_toks, 8,  toks);    /* G6 */

	element->array[element->size].node[6]
	        = nst_tok2i(num_toks, 9,  toks);    /* G7 */
	
	element->array[element->size].node[7]
	        = nst_tok2i(num_toks,10,  toks);    /* G8 */
	
	element->array[element->size].node[8]
	        = nst_tok2i(num_toks,11,  toks);    /* G9 */

	element->array[element->size].node[9]
	        = nst_tok2i(num_toks,12,  toks);    /* G10 */

	element->array[element->size].node[10]
	        = nst_tok2i(num_toks,13,  toks);    /* G11 */

	element->array[element->size].node[11]
	        = nst_tok2i(num_toks,14,  toks);    /* G12 */

	element->array[element->size].node[16]
	        = nst_tok2i(num_toks,15,  toks);    /* G13 */

	element->array[element->size].node[17]
	        = nst_tok2i(num_toks,16,  toks);    /* G14 */

	element->array[element->size].node[18]
	        = nst_tok2i(num_toks,17,  toks);    /* G15 */

	element->array[element->size].node[19]
	        = nst_tok2i(num_toks,18,  toks);    /* G16 */

	element->array[element->size].node[12]
	        = nst_tok2i(num_toks,19,  toks);    /* G17 */

	element->array[element->size].node[13]
	        = nst_tok2i(num_toks,20,  toks);    /* G18 */

	element->array[element->size].node[14]
	        = nst_tok2i(num_toks,21,  toks);    /* G19 */

	element->array[element->size].node[15]
	        = nst_tok2i(num_toks,22,  toks);    /* G20 */

	for(ii1=0; ii1<8; ii1++){
		if(element->array[element->size].node[ii1] <= 0){
			return(CONVERT_ERROR_RC);
		}
	}

	for(ii1=8; ii1<20; ii1++){
		if(element->array[element->size].node[ii1] <= 0){
			element->array[element->size].node[ii1] = -1;
			element->array[element->size].type = ELEM_HEXA1;
			element->array[element->size].node_num = 8;
		}
	}

	return(NORMAL_RC);
}


static RC
put_chexa (FILE *fp, ELEMENT element, DATA_TYPE source)
{
	IDC tmp[23];

	tmp[0].c = "CHEXA";
	tmp[1].i = element.label;     /* EID */
	tmp[2].i = element.physical;  /* PID */
	if(tmp[2].i <= 0) tmp[2].i = 1;
	tmp[3].i = element.node[0];    /* G1 */
	tmp[4].i = element.node[1];    /* G2 */
	tmp[5].i = element.node[2];    /* G3 */
	tmp[6].i = element.node[3];    /* G4 */
	tmp[7].i = element.node[4];    /* G5 */
	tmp[8].i = element.node[5];    /* G6 */
	tmp[9].i = element.node[6];    /* G7 */
	tmp[10].i = element.node[7];   /* G8 */

	if(element.type == ELEM_HEXA1){
		RC_TRY( nst_print(fp, "ciiiiiiiir-iir", tmp) );
	}else if(element.type == ELEM_HEXA2){
		tmp[11].i = element.node[8];   /* G9 */
		tmp[12].i = element.node[9];   /* G10 */
		tmp[13].i = element.node[10];  /* G11 */
		tmp[14].i = element.node[11];  /* G12 */
		tmp[15].i = element.node[16];  /* G13 */
		tmp[16].i = element.node[17];  /* G14 */
		tmp[17].i = element.node[18];  /* G15 */
		tmp[18].i = element.node[19];  /* G16 */
		tmp[19].i = element.node[12];  /* G17 */
		tmp[20].i = element.node[13];  /* G18 */
		tmp[21].i = element.node[14];  /* G19 */
		tmp[22].i = element.node[15];  /* G20 */
		RC_TRY( nst_print(fp, "ciiiiiiiir-iiiiiiiir-iiiiiir", tmp) );
	}else{
		return(ARG_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
set_cbar (ELEMENT_ARRAY *element, char **toks, int num_toks)
{
	int ii1;
	int itmp;

	itmp = nst_tok2i(num_toks, 1, toks);    /* EID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	element->array[element->size].label = itmp;
	element->array[element->size].type = ELEM_BEAM1;
	element->array[element->size].material = -1;

	itmp = nst_tok2i(num_toks, 2, toks);     /* PID */
	if(itmp < 0) itmp = element->array[element->size].label;
	element->array[element->size].physical = itmp;

	element->array[element->size].node_num = 2;

	element->array[element->size].node[0]
	        = nst_tok2i(num_toks, 3, toks);  /* GA */

	element->array[element->size].node[1]
	        = nst_tok2i(num_toks, 4, toks);  /* GB */

	for(ii1=0; ii1<2; ii1++){
		if(element->array[element->size].node[ii1] <= 0){
			return(CONVERT_ERROR_RC);
		}
	}

	element->array[element->size].d_info[0]
	        = nst_tok2d(num_toks, 5, toks);  /* X1 */

	element->array[element->size].d_info[1]
	        = nst_tok2d(num_toks, 6, toks);  /* X2 */

	element->array[element->size].d_info[2]
	        = nst_tok2d(num_toks, 7, toks);  /* X3 */

	element->array[element->size].i_info[4] = -1;   /* AZLIB CBAR flag */

	return(NORMAL_RC);
}


static RC
put_cbar (FILE *fp, ELEMENT element, DATA_TYPE source)
{

	IDC tmp[8];

	tmp[0].c = "CBAR";
	tmp[1].i = element.label;        /* EID */
	tmp[2].i = element.physical;     /* PID */
	if(tmp[2].i <= 0) tmp[2].i = 1;
	tmp[3].i = element.node[0];      /* GA */
	tmp[4].i = element.node[1];      /* GB */
	tmp[5].d = element.d_info[0];    /* X1 */
	tmp[6].d = element.d_info[1];    /* X2 */
	tmp[7].d = element.d_info[2];    /* X3 */

	RC_TRY( nst_print(fp, "ciiiidddr", tmp) );

	return(NORMAL_RC);
}	


static RC
set_cbeam (ELEMENT_ARRAY *element, char **toks, int num_toks)
{
	int ii1;
	int itmp;

	itmp = nst_tok2i(num_toks, 1, toks);    /* EID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	element->array[element->size].label = itmp;
	element->array[element->size].type = ELEM_BEAM1;
	element->array[element->size].material = -1;

	itmp = nst_tok2i(num_toks, 2, toks);     /* PID */
	if(itmp < 0) itmp = element->array[element->size].label;
	element->array[element->size].physical = itmp;

	element->array[element->size].node_num = 2;

	element->array[element->size].node[0]
	        = nst_tok2i(num_toks, 3, toks);  /* GA */

	element->array[element->size].node[1]
	        = nst_tok2i(num_toks, 4, toks);  /* GB */

	for(ii1=0; ii1<2; ii1++){
		if(element->array[element->size].node[ii1] <= 0){
			return(CONVERT_ERROR_RC);
		}
	}

	element->array[element->size].d_info[0]
	        = nst_tok2d(num_toks, 5, toks);  /* X1 */

	element->array[element->size].d_info[1]
	        = nst_tok2d(num_toks, 6, toks);  /* X2 */

	element->array[element->size].d_info[2]
	        = nst_tok2d(num_toks, 7, toks);  /* X3 */

	element->array[element->size].i_info[0]
	        = nst_tok2i(num_toks, 8, toks);  /* PA */

	element->array[element->size].i_info[1]
	        = nst_tok2i(num_toks, 9, toks);  /* PB */

	element->array[element->size].d_info[3]
	        = nst_tok2d(num_toks, 10, toks);  /* W1A */

	element->array[element->size].d_info[4]
	        = nst_tok2d(num_toks, 11, toks);  /* W2A */

	element->array[element->size].d_info[5]
	        = nst_tok2d(num_toks, 12, toks);  /* W3A */

	element->array[element->size].d_info[6]
	        = nst_tok2d(num_toks, 13, toks);  /* W1B */

	element->array[element->size].d_info[7]
	        = nst_tok2d(num_toks, 14, toks);  /* W2B */

	element->array[element->size].d_info[8]
	        = nst_tok2d(num_toks, 15, toks);  /* W3B */

	element->array[element->size].i_info[2]
	        = nst_tok2i(num_toks, 16, toks);  /* SA */

	element->array[element->size].i_info[3]
	        = nst_tok2i(num_toks, 17, toks);  /* SB */

	element->array[element->size].i_info[4] = 1; /* AZLIB CBEAM flag */

	return(NORMAL_RC);
}


static RC
put_cbeam (FILE *fp, ELEMENT element, DATA_TYPE source)
{

	IDC tmp[18];

	tmp[0].c = "CBEAM";
	tmp[1].i = element.label;        /* EID */
	tmp[2].i = element.physical;     /* PID */
	if(tmp[2].i <= 0) tmp[2].i = 1;
	tmp[3].i = element.node[0];      /* GA */
	tmp[4].i = element.node[1];      /* GB */
	tmp[5].d = element.d_info[0];    /* X1 */
	tmp[6].d = element.d_info[1];    /* X2 */
	tmp[7].d = element.d_info[2];    /* X3 */
	tmp[8].i = element.i_info[0];    /* PA */
	tmp[9].i = element.i_info[1];    /* PB */
	tmp[10].d = element.d_info[3];   /* W1A */
	tmp[11].d = element.d_info[4];   /* W2A */
	tmp[12].d = element.d_info[5];   /* W3A */
	tmp[13].d = element.d_info[6];   /* W1B */
	tmp[14].d = element.d_info[7];   /* W2B */
	tmp[15].d = element.d_info[8];   /* W3B */
	tmp[16].i = element.i_info[2];   /* SA */
	tmp[17].i = element.i_info[3];   /* SB */

	RC_TRY( nst_print(fp, "ciiiidddr-iiddddddr-iir", tmp) );

	return(NORMAL_RC);
}	


RC
nst_input_material_prop (FILE *fp, MATERIAL_PROP_ARRAY *material)
{
	const int num_keys = 7;
	NST_BULK_KEY keys[] = { {"MAT1", 4}, {"MAT2", 4}, {"MAT4", 4}, 
                            {"MAT8", 4}, {"MAT9", 4}, {"MAT5", 4}, 
                            {"MAT10", 5} };
	int hit_key,ii1;
	char *toks[TOKS_SIZE];
	int num_toks;
	NST_FILE nst_file;

	RC_TRY( nst_allocate_tok(toks, TOKS_SIZE) );
	RC_TRY( init_nst_file(fp, &nst_file) );
	RC_TRY( nst_seek_begin_bulk(&nst_file) );
	RC_TRY( allocate_material_prop_array(INIT_SIZE_SMALL, material) );
	material->source = FEM_NASTRAN;

	while(1){
		RC rc = nst_blk_get_tok(&nst_file, num_keys, keys, &hit_key,
		                        TOKS_SIZE, toks, &num_toks);
		if(rc != NORMAL_RC){
			if(rc == END_RC) break;
			return(rc);
		}

		RC_TRY( realloc_material_prop_array(material) );

		switch(hit_key){
		case 0:
			RC_TRY( set_mat1(material, toks, num_toks) );
			break;
		case 1:
			RC_TRY( set_mat2(material, toks, num_toks) );
			break;
		case 2:
			RC_TRY( set_mat4(material, toks, num_toks) );
			break;
		case 3:
			RC_TRY( set_mat8(material, toks, num_toks) );
			break;
		case 4:
			RC_TRY( set_mat9(material, toks, num_toks) );
			break;
		case 5:
			RC_TRY( set_mat5(material, toks, num_toks) );
			break;
        case 6:
			RC_TRY( set_mat10(material, toks, num_toks) );
            break;
		default:
			return(UNKNOWN_ERROR_RC);
		}
		material->size++;
	}

	for(ii1=0; ii1<material->size;ii1++){
		RC_TRY( fill_material(&(material->array[ii1])) );
	}
	
	RC_TRY( clean_material_prop_array(material) );
	RC_TRY( term_nst_file(&nst_file) );
	RC_TRY( nst_free_tok(toks, TOKS_SIZE) );

	return(NORMAL_RC);
}


RC
nst_output_material_prop (FILE *fp, MATERIAL_PROP_ARRAY material)
{
	int ii1;

	for(ii1=0; ii1<material.size; ii1++){
		if(material.array[ii1].label < 0) continue;

        if(material.array[ii1].mat_state == STATE_SOLID){
    		if(material.array[ii1].mat_type == MAT_ISOTROPIC){
    			RC_TRY( put_mat1(fp, material.array[ii1], material.source) );
    		}else if(material.array[ii1].mat_type == MAT_ORTHOTROPIC){
    			if( (material.array[ii1].ana_type == ANA_PLANE_STRESS)
    			  ||(material.array[ii1].ana_type == ANA_PLANE_STRAIN) ){
    				RC_TRY( put_mat8(fp, material.array[ii1], material.source) );
    			}else{
    				return(IMPLEMENT_ERROR_RC);
    			}
	    	}else if(material.array[ii1].mat_type == MAT_ANISOTROPIC){
	    		if( (material.array[ii1].ana_type == ANA_PLANE_STRESS)
	    		  ||(material.array[ii1].ana_type == ANA_PLANE_STRAIN) ){
	    			RC_TRY( put_mat2(fp, material.array[ii1], material.source) );
	    		}else if(material.array[ii1].ana_type == ANA_3D){
	    			RC_TRY( put_mat9(fp, material.array[ii1], material.source) );
	    		}else{
	    			return(IMPLEMENT_ERROR_RC);
	    		}
	    	}else{
	    		return(IMPLEMENT_ERROR_RC);
	    	}
        }else if(material.array[ii1].mat_state == STATE_FLUID){
            RC_TRY( put_mat10(fp, material.array[ii1], material.source) );
        }
	}

	return(NORMAL_RC);
}


static RC
set_mat1 (MATERIAL_PROP_ARRAY *material, char **toks, int num_toks)
{
	int itmp, index, tmp = 0;

	itmp = nst_tok2i(num_toks, 1, toks);    /* MID */
	if(itmp < 0) return(CONVERT_ERROR_RC);

	index = search_material_prop_label(*material, itmp);
	if(index < 0){
		tmp = material->size;
	}else{
		tmp = index;
		(material->size)--;
	}

	material->array[tmp].label = itmp;

	material->array[tmp].mat_type = MAT_ISOTROPIC;

	material->array[tmp].mat_state = STATE_SOLID;

	material->array[tmp].ana_type = ANA_3D;

	material->array[tmp].E
	        = nst_tok2d(num_toks, 2, toks);  /* E */

	material->array[tmp].G
	        = nst_tok2d(num_toks, 3, toks);  /* G */

	material->array[tmp].nu
	        = nst_tok2d(num_toks, 4, toks);  /* NU */

	material->array[tmp].rho
	        = nst_tok2d(num_toks, 5, toks);  /* RHO */
	if(material->array[tmp].rho <= CR_FNSDEFAULT){
		material->array[tmp].rho = 0.0;
	}

	material->array[tmp].alpha
	        = nst_tok2d(num_toks, 6, toks);  /* ALPHA */
	if(material->array[tmp].alpha <= CR_FNSDEFAULT){
		material->array[tmp].alpha = 0.0;
	}

	material->array[tmp].tref
	        = nst_tok2d(num_toks, 7, toks);  /* TREF */

	material->array[tmp].damp
	        = nst_tok2d(num_toks, 8, toks);  /* GE */

	material->array[tmp].d_info[1]
	        = nst_tok2d(num_toks, 9, toks);  /* ST */

	material->array[tmp].d_info[2]
	        = nst_tok2d(num_toks,10, toks);  /* SC */

	material->array[tmp].d_info[3]
	        = nst_tok2d(num_toks,11, toks);  /* SS */

	material->array[tmp].i_info[0]
	        = nst_tok2i(num_toks,12, toks);  /* MCSID */

	return(NORMAL_RC);
}


static RC
put_mat1 (FILE *fp, MATERIAL_PROP mat, DATA_TYPE source)
{
	IDC tmp[9];

	tmp[0].c = "MAT1";
	tmp[1].i = mat.label;     /* MID */
	tmp[2].d = mat.E;         /* E */
	tmp[3].d = mat.G;         /* G */
	tmp[4].d = mat.nu;        /* NU */
	tmp[5].d = mat.rho;       /* RHO */
	tmp[6].d = mat.alpha;     /* A */
	tmp[7].d = mat.tref;      /* TREF */
	if(source == FEM_NASTRAN){
		tmp[8].d = mat.damp;  /* GE */
		if(nearly_eq(tmp[5].d, 0.0)) tmp[5].d = FNSDEFAULT; /* RHO */
		if(nearly_eq(tmp[6].d, 0.0)) tmp[6].d = FNSDEFAULT; /* A */
	}else{
		tmp[8].d = FNSDEFAULT;
	}
	RC_TRY( nst_print(fp, "cidddddddr", tmp) );

	if((source != FEM_NASTRAN)||((mat.d_info[1] <= CR_FNSDEFAULT)
	                           &&(mat.d_info[2] <= CR_FNSDEFAULT)
	                           &&(mat.d_info[3] <= CR_FNSDEFAULT)
	                           &&(mat.i_info[0] == NSDEFAULT))){
		return(NORMAL_RC);
	}

	tmp[0].d = mat.d_info[1];  /* ST */
	tmp[1].d = mat.d_info[2];  /* SC */
	tmp[2].d = mat.d_info[3];  /* SS */
	tmp[3].i = mat.i_info[0];  /* MCSID */
	RC_TRY( nst_print(fp, "-dddir", tmp) );

	return(NORMAL_RC);
}


static RC
set_mat2 (MATERIAL_PROP_ARRAY *material, char **toks, int num_toks)
{
	int itmp;

	itmp = nst_tok2i(num_toks, 1, toks);    /* MID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	material->array[material->size].label = itmp;

	material->array[material->size].mat_type = MAT_ANISOTROPIC;

	material->array[material->size].mat_state = STATE_SOLID;

	material->array[material->size].ana_type = ANA_PLANE_STRESS;

	material->array[material->size].D_matrix[0][0]
	        = nst_tok2d(num_toks, 2, toks);  /* G11 */

	material->array[material->size].D_matrix[0][1]
	        = material->array[material->size].D_matrix[1][0]
	        = nst_tok2d(num_toks, 3, toks);  /* G12 */

	material->array[material->size].D_matrix[0][2]
	        = material->array[material->size].D_matrix[2][0]
	        = nst_tok2d(num_toks, 4, toks);  /* G13 */

	material->array[material->size].D_matrix[1][1]
	        = nst_tok2d(num_toks, 5, toks);  /* G22 */

	material->array[material->size].D_matrix[1][2]
	        = material->array[material->size].D_matrix[2][1]
	        = nst_tok2d(num_toks, 6, toks);  /* G23 */

	material->array[material->size].D_matrix[2][2]
	        = nst_tok2d(num_toks, 7, toks);  /* G33 */

	material->array[material->size].rho
	        = nst_tok2d(num_toks, 8, toks);  /* RHO */
	if(material->array[material->size].rho <= CR_FNSDEFAULT){
		material->array[material->size].rho = 0.0;
	}

	material->array[material->size].alpha_vect.x
	        = nst_tok2d(num_toks, 9, toks);  /* A1 */
	if(material->array[material->size].alpha_vect.x <= CR_FNSDEFAULT){
		material->array[material->size].alpha_vect.x = 0.0;
	}

	material->array[material->size].alpha_vect.y
	        = nst_tok2d(num_toks,10, toks);  /* A2 */
	if(material->array[material->size].alpha_vect.y <= CR_FNSDEFAULT){
		material->array[material->size].alpha_vect.y = 0.0;
	}

	material->array[material->size].alpha_vect.xy
	        = nst_tok2d(num_toks,11, toks);  /* A3 */
	if(material->array[material->size].alpha_vect.xy <= CR_FNSDEFAULT){
		material->array[material->size].alpha_vect.xy = 0.0;
	}

	material->array[material->size].tref
	        = nst_tok2d(num_toks,12, toks);  /* TREF */

	material->array[material->size].d_info[0]
	        = nst_tok2d(num_toks,13, toks);  /* GE */

	material->array[material->size].d_info[1]
	        = nst_tok2d(num_toks,14, toks);  /* ST */

	material->array[material->size].d_info[2]
	        = nst_tok2d(num_toks,15, toks);  /* SC */

	material->array[material->size].d_info[3]
	        = nst_tok2d(num_toks,16, toks);  /* SS */

	material->array[material->size].i_info[0]
	        = nst_tok2i(num_toks,17, toks);  /* MCSID */

	return(NORMAL_RC);
}


static RC
put_mat2 (FILE *fp, MATERIAL_PROP mat, DATA_TYPE source)
{
	IDC tmp[17];

	tmp[0].c = "MAT2";
	tmp[1].i = mat.label;             /* MID */
	tmp[2].d = mat.D_matrix[0][0];    /* G11 */
	tmp[3].d = mat.D_matrix[0][1];    /* G12 */
	tmp[4].d = mat.D_matrix[0][2];    /* G13 */
	tmp[5].d = mat.D_matrix[1][1];    /* G22 */
	tmp[6].d = mat.D_matrix[1][2];    /* G23 */
	tmp[7].d = mat.D_matrix[2][2];    /* G33 */
	tmp[8].d = mat.rho;               /* RHO */
	tmp[9].d = mat.alpha_vect.x;      /* A1 */
	tmp[10].d = mat.alpha_vect.y;     /* A2 */
	tmp[11].d = mat.alpha_vect.xy;    /* A3 */
	tmp[12].d = mat.tref;             /* TREF */
	if(source == FEM_NASTRAN){
		tmp[13].d = mat.d_info[0];     /* GE */
		tmp[14].d = mat.d_info[1];     /* ST */
		tmp[15].d = mat.d_info[2];    /* SC */
		tmp[16].d = mat.d_info[3];    /* SS */
	}else{
		tmp[13].d = FNSDEFAULT;
		tmp[14].d = FNSDEFAULT;
		tmp[15].d = FNSDEFAULT;
		tmp[16].d = FNSDEFAULT;
	}
	RC_TRY( nst_print(fp, "cidddddddr-ddddddddr", tmp) );

	if((source != FEM_NASTRAN)||(mat.i_info[0] == NSDEFAULT)){
		return(NORMAL_RC);
	}

	tmp[0].i = mat.i_info[0];    /* MCSID */
	RC_TRY( nst_print(fp, "-ir", tmp) );

	return(NORMAL_RC);
}

static RC
set_mat4 (MATERIAL_PROP_ARRAY *material, char **toks, int num_toks)
{
	int itmp, index, tmp = 0;

	itmp = nst_tok2i(num_toks, 1, toks);    /* MID */
	if(itmp < 0) return(CONVERT_ERROR_RC);

	index = search_material_prop_label(*material, itmp);
	if(index < 0){
		tmp = material->size;
	}else{
		tmp = index;
		(material->size)--;
	}
	
	material->array[tmp].label = itmp;

	material->array[tmp].mat_type = MAT_ISOTROPIC;

	material->array[tmp].mat_state = STATE_SOLID;

	material->array[tmp].ana_type = ANA_3D;

	material->array[tmp].K
	        = nst_tok2d(num_toks, 2, toks);  /* K */
	
	material->array[tmp].sheat
	        = nst_tok2d(num_toks, 3, toks);  /* sheat */

	material->array[tmp].rho
	        = nst_tok2d(num_toks, 4, toks);  /* RHO */

	return(NORMAL_RC);
}


static RC
set_mat5 (MATERIAL_PROP_ARRAY *material, char **toks, int num_toks)
{
	int itmp, index, tmp = 0;

	itmp = nst_tok2i(num_toks, 1, toks);    /* MID */
	if(itmp < 0) return(CONVERT_ERROR_RC);

	index = search_material_prop_label(*material, itmp);
	if(index < 0){
		tmp = material->size;
	}else{
		tmp = index;
		(material->size)--;
	}
	
	material->array[tmp].label = itmp;
	
	material->array[tmp].mat_type = MAT_ORTHOTROPIC;
	
	material->array[tmp].ana_type = ANA_PLANE_STRESS;
	
	material->array[tmp].K_vect.x  = nst_tok2d(num_toks, 2, toks);
	material->array[tmp].K_vect.xy = nst_tok2d(num_toks, 3, toks);
	material->array[tmp].K_vect.zx = nst_tok2d(num_toks, 4, toks);
	material->array[tmp].K_vect.y  = nst_tok2d(num_toks, 5, toks);
	material->array[tmp].K_vect.yz = nst_tok2d(num_toks, 6, toks);
	material->array[tmp].K_vect.z  = nst_tok2d(num_toks, 7, toks);
	
	return(NORMAL_RC);
}


static RC
set_mat8 (MATERIAL_PROP_ARRAY *material, char **toks, int num_toks)
{
	int itmp, index, tmp = 0;

	itmp = nst_tok2i(num_toks, 1, toks);    /* MID */
	if(itmp < 0) return(CONVERT_ERROR_RC);

	index = search_material_prop_label(*material, itmp);
	if(index < 0){
		tmp = material->size;
	}else{
		tmp = index;
		(material->size)--;
	}

	material->array[tmp].label = itmp;

	material->array[tmp].mat_type = MAT_ORTHOTROPIC;

	material->array[tmp].mat_state = STATE_SOLID;

	material->array[tmp].ana_type = ANA_PLANE_STRESS;

	material->array[tmp].E_vect.x
	        = nst_tok2d(num_toks, 2, toks);  /* E1 */

	material->array[tmp].E_vect.y
	        = nst_tok2d(num_toks, 3, toks);  /* E2 */

	material->array[tmp].nu_vect.xy
	        = nst_tok2d(num_toks, 4, toks);  /* NU12 */

	material->array[tmp].G_vect.xy
	        = nst_tok2d(num_toks, 5, toks);  /* G12 */

	material->array[tmp].G_vect.zx
	        = nst_tok2d(num_toks, 6, toks);  /* G1Z */

	material->array[tmp].G_vect.yz
	        = nst_tok2d(num_toks, 7, toks);  /* G2Z */

	material->array[tmp].rho
	        = nst_tok2d(num_toks, 8, toks);  /* RHO */
	if(material->array[tmp].rho <= CR_FNSDEFAULT){
		material->array[tmp].rho = 0.0;
	}

	material->array[tmp].alpha_vect.x
	        = nst_tok2d(num_toks, 9, toks);  /* A1 */
	if(material->array[tmp].alpha_vect.x <= CR_FNSDEFAULT){
		material->array[tmp].alpha_vect.x = 0.0;
	}

	material->array[tmp].alpha_vect.y
	        = nst_tok2d(num_toks,10, toks);  /* A2 */
	if(material->array[tmp].alpha_vect.y <= CR_FNSDEFAULT){
		material->array[tmp].alpha_vect.y = 0.0;
	}

	material->array[tmp].tref
	        = nst_tok2d(num_toks,11, toks);  /* TREF */

	material->array[tmp].d_info[1]
	        = nst_tok2d(num_toks,12, toks);  /* Xt */

	material->array[tmp].d_info[2]
	        = nst_tok2d(num_toks,13, toks);  /* Xc */

	material->array[tmp].d_info[4]
	        = nst_tok2d(num_toks,14, toks);  /* Yt */

	material->array[tmp].d_info[5]
	        = nst_tok2d(num_toks,15, toks);  /* Yc */

	material->array[tmp].d_info[3]
	        = nst_tok2d(num_toks,16, toks);  /* S */

	material->array[tmp].d_info[0]
	        = nst_tok2d(num_toks,17, toks);  /* GE */

	material->array[tmp].d_info[6]
	        = nst_tok2d(num_toks,18, toks);  /* F12 */

	material->array[tmp].d_info[7]
	        = nst_tok2d(num_toks,19, toks);  /* STRN */

	return(NORMAL_RC);
}


static RC
put_mat8 (FILE *fp, MATERIAL_PROP mat, DATA_TYPE source)
{
	IDC tmp[17];

	tmp[0].c = "MAT8";
	tmp[1].i = mat.label;             /* MID */
	tmp[2].d = mat.E_vect.x;          /* E1 */
	tmp[3].d = mat.E_vect.y;          /* E2 */
	tmp[4].d = mat.nu_vect.xy;        /* NU12 */
	tmp[5].d = mat.G_vect.xy;         /* G12 */
	tmp[6].d = mat.G_vect.zx;         /* G1Z */
	tmp[7].d = mat.G_vect.yz;         /* G2Z */
	tmp[8].d = mat.rho;               /* RHO */
	tmp[9].d = mat.alpha_vect.x;      /* A1 */
	tmp[10].d = mat.alpha_vect.y;     /* A2 */
	tmp[11].d = mat.tref;             /* TREF */
	if(source == FEM_NASTRAN){
		tmp[12].d = mat.d_info[1];    /* Xt */
		tmp[13].d = mat.d_info[2];    /* Xc */
		tmp[14].d = mat.d_info[4];    /* Yt */
		tmp[15].d = mat.d_info[5];    /* Yc */
		tmp[16].d = mat.d_info[3];    /* S */
	}else{
		tmp[12].d = FNSDEFAULT;
		tmp[13].d = FNSDEFAULT;
		tmp[14].d = FNSDEFAULT;
		tmp[15].d = FNSDEFAULT;
		tmp[16].d = FNSDEFAULT;
	}
	RC_TRY( nst_print(fp, "cidddddddr-ddddddddr", tmp) );

	if( (source != FEM_NASTRAN)||((mat.d_info[0] <= CR_FNSDEFAULT)
	                            &&(mat.d_info[6] <= CR_FNSDEFAULT)
	                            &&(mat.d_info[7] <= CR_FNSDEFAULT)) ){
		return(NORMAL_RC);
	}

	tmp[0].d = mat.d_info[0];    /* GE */
	tmp[1].d = mat.d_info[6];    /* F12 */
	tmp[2].d = mat.d_info[7];    /* STRN */
	RC_TRY( nst_print(fp, "-dddr", tmp) );

	return(NORMAL_RC);
}


static RC
set_mat9 (MATERIAL_PROP_ARRAY *material, char **toks, int num_toks)
{
	int itmp;

	itmp = nst_tok2i(num_toks, 1, toks);    /* MID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	material->array[material->size].label = itmp;

	material->array[material->size].mat_type = MAT_ANISOTROPIC;

	material->array[material->size].mat_state = STATE_SOLID;

	material->array[material->size].ana_type = ANA_3D;

	material->array[material->size].D_matrix[0][0]
	        = nst_tok2d(num_toks, 2, toks);  /* G11 */

	material->array[material->size].D_matrix[0][1]
	        = material->array[material->size].D_matrix[1][0]
	        = nst_tok2d(num_toks, 3, toks);  /* G12 */

	material->array[material->size].D_matrix[0][2]
	        = material->array[material->size].D_matrix[2][0]
	        = nst_tok2d(num_toks, 4, toks);  /* G13 */

	material->array[material->size].D_matrix[0][5]
	        = material->array[material->size].D_matrix[5][0]
	        = nst_tok2d(num_toks, 5, toks);  /* G14 */

	material->array[material->size].D_matrix[0][3]
	        = material->array[material->size].D_matrix[3][0]
	        = nst_tok2d(num_toks, 6, toks);  /* G15 */

	material->array[material->size].D_matrix[0][4]
	        = material->array[material->size].D_matrix[4][0]
	        = nst_tok2d(num_toks, 7, toks);  /* G16 */

	material->array[material->size].D_matrix[1][1]
	        = nst_tok2d(num_toks, 8, toks);  /* G22 */

	material->array[material->size].D_matrix[1][2]
	        = material->array[material->size].D_matrix[2][1]
	        = nst_tok2d(num_toks, 9, toks);  /* G23 */

	material->array[material->size].D_matrix[1][5]
	        = material->array[material->size].D_matrix[5][1]
	        = nst_tok2d(num_toks,10, toks);  /* G24 */

	material->array[material->size].D_matrix[1][3]
	        = material->array[material->size].D_matrix[3][1]
	        = nst_tok2d(num_toks,11, toks);  /* G25 */

	material->array[material->size].D_matrix[1][4]
	        = material->array[material->size].D_matrix[4][1]
	        = nst_tok2d(num_toks,12, toks);  /* G26 */

	material->array[material->size].D_matrix[2][2]
	        = nst_tok2d(num_toks,13, toks);  /* G33 */

	material->array[material->size].D_matrix[2][5]
	        = material->array[material->size].D_matrix[5][2]
	        = nst_tok2d(num_toks,14, toks);  /* G34 */

	material->array[material->size].D_matrix[2][3]
	        = material->array[material->size].D_matrix[3][2]
	        = nst_tok2d(num_toks,15, toks);  /* G35 */

	material->array[material->size].D_matrix[2][4]
	        = material->array[material->size].D_matrix[4][2]
	        = nst_tok2d(num_toks,16, toks);  /* G36 */

	material->array[material->size].D_matrix[5][5]
	        = nst_tok2d(num_toks,17, toks);  /* G44 */

	material->array[material->size].D_matrix[5][3]
	        = material->array[material->size].D_matrix[3][5]
	        = nst_tok2d(num_toks,18, toks);  /* G45 */

	material->array[material->size].D_matrix[5][4]
	        = material->array[material->size].D_matrix[4][5]
	        = nst_tok2d(num_toks,19, toks);  /* G46 */

	material->array[material->size].D_matrix[3][3]
	        = nst_tok2d(num_toks,20, toks);  /* G55 */

	material->array[material->size].D_matrix[3][4]
	        = material->array[material->size].D_matrix[4][3]
	        = nst_tok2d(num_toks,21, toks);  /* G56 */

	material->array[material->size].D_matrix[4][4]
	        = nst_tok2d(num_toks,22, toks);  /* G66 */

	material->array[material->size].rho
	        = nst_tok2d(num_toks,23, toks);  /* RHO */
	if(material->array[material->size].rho <= CR_FNSDEFAULT){
		material->array[material->size].rho = 0.0;
	}

	material->array[material->size].alpha_vect.x
	        = nst_tok2d(num_toks,24, toks);  /* A1 */
	if(material->array[material->size].alpha_vect.x <= CR_FNSDEFAULT){
		material->array[material->size].alpha_vect.x = 0.0;
	}

	material->array[material->size].alpha_vect.y
	        = nst_tok2d(num_toks,25, toks);  /* A2 */
	if(material->array[material->size].alpha_vect.y <= CR_FNSDEFAULT){
		material->array[material->size].alpha_vect.y = 0.0;
	}

	material->array[material->size].alpha_vect.z
	        = nst_tok2d(num_toks,26, toks);  /* A3 */
	if(material->array[material->size].alpha_vect.z <= CR_FNSDEFAULT){
		material->array[material->size].alpha_vect.z = 0.0;
	}

	material->array[material->size].alpha_vect.xy
	        = nst_tok2d(num_toks,27, toks);  /* A4 */
	if(material->array[material->size].alpha_vect.xy <= CR_FNSDEFAULT){
		material->array[material->size].alpha_vect.xy = 0.0;
	}

	material->array[material->size].alpha_vect.yz
	        = nst_tok2d(num_toks,28, toks);  /* A5 */
	if(material->array[material->size].alpha_vect.yz <= CR_FNSDEFAULT){
		material->array[material->size].alpha_vect.yz = 0.0;
	}

	material->array[material->size].alpha_vect.zx
	        = nst_tok2d(num_toks,29, toks);  /* A6 */
	if(material->array[material->size].alpha_vect.zx <= CR_FNSDEFAULT){
		material->array[material->size].alpha_vect.zx = 0.0;
	}

	material->array[material->size].tref
	        = nst_tok2d(num_toks,30, toks);  /* TREF */

	material->array[material->size].d_info[0]
	        = nst_tok2d(num_toks,31, toks);  /* GE */

	return(NORMAL_RC);
}


static RC
put_mat9 (FILE *fp, MATERIAL_PROP mat, DATA_TYPE source)
{
	IDC tmp[32];

	tmp[0].c = "MAT9";
	tmp[1].i = mat.label;            /* MID */
	tmp[2].d = mat.D_matrix[0][0];   /* G11 */
	tmp[3].d = mat.D_matrix[0][1];   /* G12 */
	tmp[4].d = mat.D_matrix[0][2];   /* G13 */
	tmp[5].d = mat.D_matrix[0][5];   /* G14 */
	tmp[6].d = mat.D_matrix[0][3];   /* G15 */
	tmp[7].d = mat.D_matrix[0][4];   /* G16 */
	tmp[8].d = mat.D_matrix[1][1];   /* G22 */
	tmp[9].d = mat.D_matrix[1][2];   /* G23 */
	tmp[10].d = mat.D_matrix[1][5];  /* G24 */
	tmp[11].d = mat.D_matrix[1][3];  /* G25 */
	tmp[12].d = mat.D_matrix[1][4];  /* G26 */
	tmp[13].d = mat.D_matrix[2][2];  /* G33 */
	tmp[14].d = mat.D_matrix[2][5];  /* G34 */
	tmp[15].d = mat.D_matrix[2][3];  /* G35 */
	tmp[16].d = mat.D_matrix[2][4];  /* G36 */
	tmp[17].d = mat.D_matrix[5][5];  /* G44 */
	tmp[18].d = mat.D_matrix[3][5];  /* G45 */
	tmp[19].d = mat.D_matrix[4][5];  /* G46 */
	tmp[20].d = mat.D_matrix[3][3];  /* G55 */
	tmp[21].d = mat.D_matrix[3][4];  /* G56 */
	tmp[22].d = mat.D_matrix[4][4];  /* G66 */
	tmp[23].d = mat.rho;             /* RHO */
	tmp[24].d = mat.alpha_vect.x;    /* A1 */
	tmp[25].d = mat.alpha_vect.y;    /* A2 */
	tmp[26].d = mat.alpha_vect.z;    /* A3 */
	tmp[27].d = mat.alpha_vect.xy;   /* A4 */
	tmp[28].d = mat.alpha_vect.yz;   /* A5 */
	tmp[29].d = mat.alpha_vect.zx;   /* A6 */
	tmp[30].d = mat.tref;            /* TREF */
	if(source == FEM_NASTRAN){
		tmp[31].d = mat.d_info[0];   /* GE */
	}else{
		tmp[31].d = FNSDEFAULT;
	}
	RC_TRY( nst_print(fp, "cidddddddr-ddddddddr-ddddddddr-dddddddr", tmp) );

	return(NORMAL_RC);
}

static RC
set_mat10 (MATERIAL_PROP_ARRAY *material, char **toks, int num_toks)
{
	int itmp, index, tmp = 0;

	itmp = nst_tok2i(num_toks, 1, toks);    /* MID */
	if(itmp < 0) return(CONVERT_ERROR_RC);

	index = search_material_prop_label(*material, itmp);
	if(index < 0){
		tmp = material->size;
	}else{
		tmp = index;
		(material->size)--;
	}

	material->array[tmp].label = itmp;

	material->array[tmp].mat_type = MAT_ISOTROPIC;

	material->array[tmp].mat_state = STATE_FLUID;

	material->array[tmp].ana_type = ANA_3D;

	material->array[tmp].bulk
	        = nst_tok2d(num_toks, 2, toks);  /* BULK */

	material->array[tmp].rho
	        = nst_tok2d(num_toks, 3, toks);  /* RHO */

	material->array[tmp].speed
	        = nst_tok2d(num_toks, 4, toks);  /* C */

	material->array[tmp].damp
	        = nst_tok2d(num_toks, 5, toks);  /* GE */

	material->array[tmp].d_info[0]
	        = nst_tok2d(num_toks, 6, toks);  /* ALPHA */

	return(NORMAL_RC);
}

static RC
put_mat10 (FILE *fp, MATERIAL_PROP mat, DATA_TYPE source)
{
	IDC tmp[8];

	tmp[0].c = "MAT10";
	tmp[1].i = mat.label;      /* MID */
	tmp[2].d = mat.bulk;       /* BULK */
	tmp[3].d = mat.rho;        /* RHO */
	tmp[4].d = mat.speed;      /* C */
	tmp[5].d = mat.damp;       /* GE */
	tmp[6].d = mat.d_info[0];  /* ALPHA */
	RC_TRY( nst_print(fp, "ciddddr", tmp) );

	return(NORMAL_RC);
}

RC
nst_input_physical_prop (FILE *fp, PHYSICAL_PROP_ARRAY *phys)
{
	const int num_keys = 6;
	NST_BULK_KEY keys[] = { {"PSHELL", 6}, {"PSOLID", 6},
	                        {"PBAR", 4},   {"PBARL",  5},
	                        {"PELAS",  5}, {"PBEAM", 5} };
	int hit_key;
	char *toks[TOKS_SIZE];
	int num_toks;
	NST_FILE nst_file;

	RC_TRY( nst_allocate_tok(toks, TOKS_SIZE) );
	RC_TRY( init_nst_file(fp, &nst_file) );
	RC_TRY( nst_seek_begin_bulk(&nst_file) );
	RC_TRY( allocate_physical_prop_array(INIT_SIZE_SMALL, phys) );
	phys->source = FEM_NASTRAN;

	while(1){
		RC rc = nst_blk_get_tok(&nst_file, num_keys, keys, &hit_key,
		                        TOKS_SIZE, toks, &num_toks);
		if(rc != NORMAL_RC){
			if(rc == END_RC) break;
			return(rc);
		}

		RC_TRY( realloc_physical_prop_array(phys) );

		switch(hit_key){
		case 0:
			RC_TRY( set_pshell(phys, toks, num_toks) );
			break;
		case 1:
			RC_TRY( set_psolid(phys, toks, num_toks) );
			break;
		case 2:
			RC_TRY( set_pbar(phys, toks, num_toks) );
			break;
		case 3:
			RC_TRY( set_pbarl(phys, toks, num_toks) );
			break;
		case 4:
			RC_TRY( set_pelas(phys, toks, num_toks) );
			break;
		case 5:
			RC_TRY( set_pbeam(phys, toks, num_toks) );
			break;
		default:
			return(UNKNOWN_ERROR_RC);
		}
		(phys->size)++;
	}

	RC_TRY( clean_physical_prop_array(phys) );
	RC_TRY( term_nst_file(&nst_file) );
	RC_TRY( nst_free_tok(toks, TOKS_SIZE) );

	return(NORMAL_RC);
}


RC
nst_output_physical_prop (FILE *fp, PHYSICAL_PROP_ARRAY physical)
{
	int ii1;

	for(ii1=0; ii1<physical.size; ii1++){
		if(physical.array[ii1].label < 0) continue;

		if(physical.source == FEM_NASTRAN){
			switch(physical.array[ii1].i_info[1]){
			case NST_PROP_PSHELL:
				RC_TRY( put_pshell(fp, physical.array[ii1], physical.source) );
				break;
			case NST_PROP_PSOLID:
				RC_TRY( put_psolid(fp, physical.array[ii1], physical.source) );
				break;
			case NST_PROP_PBAR:
				RC_TRY( put_pbar(fp, physical.array[ii1], physical.source) );
				break;
			case NST_PROP_PBARL:
				RC_TRY( put_pbarl(fp, physical.array[ii1], physical.source) );
				break;
			case NST_PROP_PELAS:
				RC_TRY( put_pelas(fp, physical.array[ii1], physical.source) );
				break;
			case NST_PROP_PBEAM:
				RC_TRY( put_pbeam(fp, physical.array[ii1], physical.source) );
				break;
			default:
				return(IMPLEMENT_ERROR_RC);
			}
		}else{
			return(IMPLEMENT_ERROR_RC);
		}
	}

	return(NORMAL_RC);
}


static RC
set_pshell (PHYSICAL_PROP_ARRAY *phys, char **toks, int num_toks)
{
	int itmp;

	itmp = nst_tok2i(num_toks, 1, toks);    /* PID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	phys->array[phys->size].label = itmp;

	phys->array[phys->size].i_info[1] = NST_PROP_PSHELL;

	phys->array[phys->size].i_info[0]
	    = nst_tok2i(num_toks, 2, toks);     /* MID1 */

	phys->array[phys->size].d_info[0]
	    = nst_tok2d(num_toks, 3, toks);     /* T */

	phys->array[phys->size].i_info[2]
	    = nst_tok2i(num_toks, 4, toks);     /* MID2 */

	phys->array[phys->size].d_info[1]
	    = nst_tok2d(num_toks, 5, toks);     /* 12I/T**3 */

	phys->array[phys->size].i_info[3]
	    = nst_tok2i(num_toks, 6, toks);     /* MID3 */

	phys->array[phys->size].d_info[2]
	    = nst_tok2d(num_toks, 7, toks);     /* TS/T */

	phys->array[phys->size].d_info[3]
	    = nst_tok2d(num_toks, 8, toks);     /* NSM */

	phys->array[phys->size].d_info[4]
	    = nst_tok2d(num_toks, 0, toks);     /* Z1 */

	phys->array[phys->size].d_info[5]
	    = nst_tok2d(num_toks,10, toks);     /* Z2 */

	phys->array[phys->size].i_info[4]
	    = nst_tok2i(num_toks,11, toks);     /* MID4 */

	return(NORMAL_RC);
}


static RC
put_pshell (FILE *fp, PHYSICAL_PROP phys, DATA_TYPE source)
{
	IDC tmp[9];

	tmp[0].c = "PSHELL";
	tmp[1].i = phys.label;           /* PID */
	if(source == FEM_NASTRAN){
		tmp[2].i = phys.i_info[0];   /* MID1 */
		tmp[3].d = phys.d_info[0];   /* T */
		tmp[4].i = phys.i_info[2];   /* MID2 */
		tmp[5].d = phys.d_info[1];   /* 12I/T**3 */
		tmp[6].i = phys.i_info[3];   /* MID3 */
		tmp[7].d = phys.d_info[2];   /* TS/T */
		tmp[8].d = phys.d_info[3];   /* NSM */
	}else{
		tmp[2].i = NSDEFAULT;
		tmp[3].d = FNSDEFAULT;
		tmp[4].i = NSDEFAULT;
		tmp[5].d = FNSDEFAULT;
		tmp[6].i = NSDEFAULT;
		tmp[7].d = FNSDEFAULT;
		tmp[8].d = FNSDEFAULT;
	}
	RC_TRY( nst_print(fp, "ciidididdr", tmp) );

	if( (source != FEM_NASTRAN)||( (phys.d_info[4] <= CR_FNSDEFAULT)
	                             &&(phys.d_info[5] <= CR_FNSDEFAULT)
	                             &&(phys.i_info[4] == NSDEFAULT) ) ){
		return(NORMAL_RC);
	}

	tmp[0].d = phys.d_info[4];   /* Z1 */
	tmp[1].d = phys.d_info[5];   /* Z2 */
	tmp[2].i = phys.i_info[4];   /* MID4 */
	RC_TRY( nst_print(fp, "-ddir", tmp) );

	return(NORMAL_RC);
}


static RC
set_psolid (PHYSICAL_PROP_ARRAY *phys, char **toks, int num_toks)
{
	int itmp;

	itmp = nst_tok2i(num_toks, 1, toks);    /* PID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	phys->array[phys->size].label = itmp;

	phys->array[phys->size].i_info[1] = NST_PROP_PSOLID;

	phys->array[phys->size].i_info[0]
	    = nst_tok2i(num_toks, 2, toks);     /* MID1 */

	phys->array[phys->size].i_info[2]
	    = nst_tok2i(num_toks, 3, toks);     /* CORDM */

	phys->array[phys->size].i_info[3]
	    = nst_tok2i(num_toks, 4, toks);     /* IN */

	phys->array[phys->size].i_info[4]
	    = nst_tok2i(num_toks, 5, toks);     /* STRESS */

	phys->array[phys->size].i_info[5]
	    = nst_tok2i(num_toks, 6, toks);     /* ISOP */

    /* FCTN */
    if(!nst_tok_strncmp(num_toks, 7, toks, "FFLUID", 6)){
    	phys->array[phys->size].i_info[6] = 1;
    }else if(!nst_tok_strncmp(num_toks, 7, toks, "PFLUID", 6)){
    	phys->array[phys->size].i_info[6] = 2;
    }else if(!nst_tok_strncmp(num_toks, 7, toks, "POLO", 4)){
    	phys->array[phys->size].i_info[6] = 3;
    }else{
    	phys->array[phys->size].i_info[6] = 0;
    }

	return(NORMAL_RC);
}


static RC
put_psolid (FILE *fp, PHYSICAL_PROP phys, DATA_TYPE source)
{
	IDC tmp[8];

	tmp[0].c = "PSOLID";
	tmp[1].i = phys.label;           /* PID */
	if(source == FEM_NASTRAN){
		tmp[2].i = phys.i_info[0];   /* MID1 */
		tmp[3].i = phys.i_info[2];   /* CORDM */
		tmp[4].i = phys.i_info[3];   /* IN */
		tmp[5].i = phys.i_info[4];   /* STRESS */
		tmp[6].i = phys.i_info[5];   /* ISOP */
		tmp[7].i = phys.i_info[6];   /* FCTN */
	}else{
		tmp[2].i = NSDEFAULT;
		tmp[3].i = NSDEFAULT;
		tmp[4].i = NSDEFAULT;
		tmp[5].i = NSDEFAULT;
		tmp[6].i = NSDEFAULT;
		tmp[7].i = NSDEFAULT;
	}
	RC_TRY( nst_print(fp, "ciiiiiiir", tmp) );

	return(NORMAL_RC);
}


static RC
set_pbar (PHYSICAL_PROP_ARRAY *phys, char **toks, int num_toks)
{
	int itmp;
	double dtmp;

	itmp = nst_tok2i(num_toks, 1, toks);    /* PID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	phys->array[phys->size].label = itmp;

	phys->array[phys->size].i_info[1] = NST_PROP_PBAR;

	phys->array[phys->size].i_info[0]
	    = nst_tok2i(num_toks, 2, toks);     /* MID */

	dtmp = nst_tok2d(num_toks, 3, toks);    /* A */
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	phys->array[phys->size].d_info[0] = dtmp;

	dtmp = nst_tok2d(num_toks, 4, toks);    /* I1 */
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	phys->array[phys->size].d_info[1] = dtmp;

	dtmp = nst_tok2d(num_toks, 5, toks);    /* I2 */
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	phys->array[phys->size].d_info[2] = dtmp;

	dtmp = nst_tok2d(num_toks, 6, toks);    /* J */
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	phys->array[phys->size].d_info[3] = dtmp;

	dtmp = nst_tok2d(num_toks, 7, toks);    /* NSM */
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	phys->array[phys->size].d_info[4] = dtmp;

	return(NORMAL_RC);
}


static RC
put_pbar (FILE *fp, PHYSICAL_PROP phys, DATA_TYPE source)
{
	IDC tmp[8];

	tmp[0].c = "PBAR";
	tmp[1].i = phys.label;       /* PID */
	tmp[2].i = phys.i_info[0];   /* MID1 */
	tmp[3].d = phys.d_info[0];   /* A */
	tmp[4].d = phys.d_info[1];   /* I1 */
	tmp[5].d = phys.d_info[2];   /* I2 */
	tmp[6].d = phys.d_info[3];   /* J */
	tmp[7].d = phys.d_info[4];   /* NSM */
	
	RC_TRY( nst_print(fp, "ciidddddr", tmp) );

	return(NORMAL_RC);
}


static RC
set_pbarl (PHYSICAL_PROP_ARRAY *phys, char **toks, int num_toks)
{
	int itmp;
	int dim_num;
	int ii1;
	int info_index;

	itmp = nst_tok2i(num_toks, 1, toks);    /* PID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	phys->array[phys->size].label = itmp;

	phys->array[phys->size].i_info[1] = NST_PROP_PBARL;

	phys->array[phys->size].i_info[0]
	    = nst_tok2i(num_toks, 2, toks);     /* MID */

	phys->array[phys->size].i_info[2]
	    = nst_tok2i(num_toks, 3, toks);     /* GROUP */
	
	phys->array[phys->size].section
	    = section_type(num_toks, 4,toks, &dim_num);  /* TYPE */
	if(phys->array[phys->size].section == SECTION_UNKNOWN){
		return(CONVERT_ERROR_RC);
	}

	phys->array[phys->size].i_info[3] = dim_num;
	info_index = 1;
	for(ii1=9; ii1<9+dim_num; ii1++){
		phys->array[phys->size].d_info[info_index]
		    = nst_tok2d(num_toks, ii1, toks);  /* DIMi */
		info_index++;
	}

	phys->array[phys->size].d_info[0]
	    = nst_tok2d(num_toks, 9+dim_num, toks);  /* NSM */

	return(NORMAL_RC);
}


static RC
put_pbarl (FILE *fp, PHYSICAL_PROP phys, DATA_TYPE source)
{
	IDC tmp[13];

	tmp[0].c = "PBARL";
	tmp[1].i = phys.label;           /* PID */
	if(source == FEM_NASTRAN){
		tmp[2].i = phys.i_info[0];   /* MID1 */
		tmp[3].i = phys.i_info[2];   /* GROUP */
		switch(phys.section){
		case SECTION_ROD:
			tmp[4].c = "ROD";
			tmp[5].c = "+";
			tmp[6].c = "+";
			tmp[7].d = phys.d_info[1];
			RC_TRY( nst_print(fp, "ciiic----crcdr", tmp) );
			break;

		case SECTION_TUBE:
			tmp[4].c = "TUBE";
			tmp[5].c = "+";
			tmp[6].c = "+";
			tmp[7].d = phys.d_info[1];
			tmp[8].d = phys.d_info[2];
			RC_TRY( nst_print(fp, "ciiic----crcddr", tmp) );
			break;

		case SECTION_I:
			tmp[4].c = "I";
			tmp[5].c = "+";
			tmp[6].c = "+";
			tmp[7].d = phys.d_info[1];
			tmp[8].d = phys.d_info[2];
			tmp[9].d = phys.d_info[3];
			tmp[10].d = phys.d_info[4];
			tmp[11].d = phys.d_info[5];
			tmp[12].d = phys.d_info[6];
			RC_TRY( nst_print(fp, "ciiic----crcddddddr", tmp) );
			break;

		case SECTION_CHAN:
			tmp[4].c = "CHAN";
			tmp[5].c = "+";
			tmp[6].c = "+";
			tmp[7].d = phys.d_info[1];
			tmp[8].d = phys.d_info[2];
			tmp[9].d = phys.d_info[3];
			tmp[10].d = phys.d_info[4];
			RC_TRY( nst_print(fp, "ciiic----crcddddr", tmp) );
			break;

		case SECTION_T:
			tmp[4].c = "T";
			tmp[5].c = "+";
			tmp[6].c = "+";
			tmp[7].d = phys.d_info[1];
			tmp[8].d = phys.d_info[2];
			tmp[9].d = phys.d_info[3];
			tmp[10].d = phys.d_info[4];
			RC_TRY( nst_print(fp, "ciiic----crcddddr", tmp) );
			break;

		case SECTION_BOX:
			tmp[4].c = "BOX";
			tmp[5].c = "+";
			tmp[6].c = "+";
			tmp[7].d = phys.d_info[1];
			tmp[8].d = phys.d_info[2];
			tmp[9].d = phys.d_info[3];
			tmp[10].d = phys.d_info[4];
			RC_TRY( nst_print(fp, "ciiic----crcddddr", tmp) );
			break;

		case SECTION_BAR:
			tmp[4].c = "BAR";
			tmp[5].c = "+";
			tmp[6].c = "+";
			tmp[7].d = phys.d_info[1];
			tmp[8].d = phys.d_info[2];
			RC_TRY( nst_print(fp, "ciiic----crcddr", tmp) );
			break;

		case SECTION_CROSS:
			tmp[4].c = "CROSS";
			tmp[5].c = "+";
			tmp[6].c = "+";
			tmp[7].d = phys.d_info[1];
			tmp[8].d = phys.d_info[2];
			tmp[9].d = phys.d_info[3];
			tmp[10].d = phys.d_info[4];
			RC_TRY( nst_print(fp, "ciiic----crcddddr", tmp) );
			break;

		case SECTION_H:
			tmp[4].c = "H";
			tmp[5].c = "+";
			tmp[6].c = "+";
			tmp[7].d = phys.d_info[1];
			tmp[8].d = phys.d_info[2];
			tmp[9].d = phys.d_info[3];
			tmp[10].d = phys.d_info[4];
			RC_TRY( nst_print(fp, "ciiic----crcddddr", tmp) );
			break;

		case SECTION_T1:
			tmp[4].c = "T1";
			tmp[5].c = "+";
			tmp[6].c = "+";
			tmp[7].d = phys.d_info[1];
			tmp[8].d = phys.d_info[2];
			tmp[9].d = phys.d_info[3];
			tmp[10].d = phys.d_info[4];
			RC_TRY( nst_print(fp, "ciiic----crcddddr", tmp) );
			break;

		case SECTION_I1:
			tmp[4].c = "I1";
			tmp[5].c = "+";
			tmp[6].c = "+";
			tmp[7].d = phys.d_info[1];
			tmp[8].d = phys.d_info[2];
			tmp[9].d = phys.d_info[3];
			tmp[10].d = phys.d_info[4];
			RC_TRY( nst_print(fp, "ciiic----crcddddr", tmp) );
			break;

		case SECTION_CHAN1:
			tmp[4].c = "CHAN1";
			tmp[5].c = "+";
			tmp[6].c = "+";
			tmp[7].d = phys.d_info[1];
			tmp[8].d = phys.d_info[2];
			tmp[9].d = phys.d_info[3];
			tmp[10].d = phys.d_info[4];
			RC_TRY( nst_print(fp, "ciiic----crcddddr", tmp) );
			break;

		case SECTION_Z:
			tmp[4].c = "Z";
			tmp[5].c = "+";
			tmp[6].c = "+";
			tmp[7].d = phys.d_info[1];
			tmp[8].d = phys.d_info[2];
			tmp[9].d = phys.d_info[3];
			tmp[10].d = phys.d_info[4];
			RC_TRY( nst_print(fp, "ciiic----crcddddr", tmp) );
			break;

		case SECTION_CHAN2:
			tmp[4].c = "CHAN2";
			tmp[5].c = "+";
			tmp[6].c = "+";
			tmp[7].d = phys.d_info[1];
			tmp[8].d = phys.d_info[2];
			tmp[9].d = phys.d_info[3];
			tmp[10].d = phys.d_info[4];
			RC_TRY( nst_print(fp, "ciiic----crcddddr", tmp) );
			break;

		case SECTION_T2:
			tmp[4].c = "T2";
			tmp[5].c = "+";
			tmp[6].c = "+";
			tmp[7].d = phys.d_info[1];
			tmp[8].d = phys.d_info[2];
			tmp[9].d = phys.d_info[3];
			tmp[10].d = phys.d_info[4];
			RC_TRY( nst_print(fp, "ciiic----crcddddr", tmp) );
			break;

		case SECTION_BOX1:
			tmp[4].c = "BOX1";
			tmp[5].c = "+";
			tmp[6].c = "+";
			tmp[7].d = phys.d_info[1];
			tmp[8].d = phys.d_info[2];
			tmp[9].d = phys.d_info[3];
			tmp[10].d = phys.d_info[4];
			tmp[11].d = phys.d_info[5];
			tmp[12].d = phys.d_info[6];
			RC_TRY( nst_print(fp, "ciiic----crcddddddr", tmp) );
			break;

		case SECTION_HEXA:
			tmp[4].c = "HEXA";
			tmp[5].c = "+";
			tmp[6].c = "+";
			tmp[7].d = phys.d_info[1];
			tmp[8].d = phys.d_info[2];
			tmp[9].d = phys.d_info[3];
			RC_TRY( nst_print(fp, "ciiic----crcdddr", tmp) );
			break;

		case SECTION_HAT:
			tmp[4].c = "HAT";
			tmp[5].c = "+";
			tmp[6].c = "+";
			tmp[7].d = phys.d_info[1];
			tmp[8].d = phys.d_info[2];
			tmp[9].d = phys.d_info[3];
			tmp[10].d = phys.d_info[4];
			RC_TRY( nst_print(fp, "ciiic----crcddddr", tmp) );
			break;

		default:
			return(CONVERT_ERROR_RC);
		
		}
	}else{
		tmp[2].i = NSDEFAULT;
		tmp[3].i = NSDEFAULT;
		tmp[4].i = NSDEFAULT;
		tmp[5].i = NSDEFAULT;
		tmp[6].i = NSDEFAULT;
		tmp[7].i = NSDEFAULT;
	}

	return(NORMAL_RC);
}


static SECTION_TYPE
section_type (int toks_size, int tok_num, char **toks, int *dim_num)
{
	if(tok_num >= toks_size){
		*dim_num = 0;
		return(SECTION_UNKNOWN);
	}else if(strncmp(toks[tok_num], "ROD ", 4) == 0){
		*dim_num = 1;
		return(SECTION_ROD);
	}else if(strncmp(toks[tok_num], "TUBE ", 5) == 0){
		*dim_num = 2;
		return(SECTION_TUBE);
	}else if(strncmp(toks[tok_num], "I ", 2) == 0){
		*dim_num = 6;
		return(SECTION_I);
	}else if(strncmp(toks[tok_num], "CHAN ", 5) == 0){
		*dim_num = 4;
		return(SECTION_CHAN);
	}else if(strncmp(toks[tok_num], "T ", 2) == 0){
		*dim_num = 4;
		return(SECTION_T);
	}else if(strncmp(toks[tok_num], "BOX ", 4) == 0){
		*dim_num = 4;
		return(SECTION_BOX);
	}else if(strncmp(toks[tok_num], "BAR ", 4) == 0){
		*dim_num = 2;
		return(SECTION_BAR);
	}else if(strncmp(toks[tok_num], "CROSS ", 6) == 0){
		*dim_num = 4;
		return(SECTION_CROSS);
	}else if(strncmp(toks[tok_num], "H ", 2) == 0){
		*dim_num = 4;
		return(SECTION_H);
	}else if(strncmp(toks[tok_num], "T1 ", 3) == 0){
		*dim_num = 3;
		return(SECTION_T1);
	}else if(strncmp(toks[tok_num], "I1 ", 3) == 0){
		*dim_num = 4;
		return(SECTION_I1);
	}else if(strncmp(toks[tok_num], "CHAN1 ", 6) == 0){
		*dim_num = 4;
		return(SECTION_CHAN1);
	}else if(strncmp(toks[tok_num], "Z ", 2) == 0){
		*dim_num = 4;
		return(SECTION_Z);
	}else if(strncmp(toks[tok_num], "CHAN2 ", 6) == 0){
		*dim_num = 4;
		return(SECTION_CHAN2);
	}else if(strncmp(toks[tok_num], "T2 ", 3) == 0){
		*dim_num = 4;
		return(SECTION_T2);
	}else if(strncmp(toks[tok_num], "BOX1 ", 5) == 0){
		*dim_num = 6;
		return(SECTION_BOX1);
	}else if(strncmp(toks[tok_num], "HEXA ", 5) == 0){
		*dim_num = 3;
		return(SECTION_HEXA);
	}else if(strncmp(toks[tok_num], "HAT ", 4) == 0){
		*dim_num = 4;
		return(SECTION_HAT);
	}

	*dim_num = 0;
	return(SECTION_UNKNOWN);
}


static RC
set_pelas (PHYSICAL_PROP_ARRAY *phys, char **toks, int num_toks)
{
	int itmp;

	itmp = nst_tok2i(num_toks, 1, toks);    /* PID1 */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	phys->array[phys->size].label = itmp;

	phys->array[phys->size].i_info[0] = -1;

	phys->array[phys->size].i_info[1] = NST_PROP_PELAS;

	phys->array[phys->size].d_info[0]
	    = nst_tok2d(num_toks, 2, toks);     /* K1 */

	phys->array[phys->size].d_info[1]
	    = nst_tok2d(num_toks, 3, toks);     /* GE1 */

	phys->array[phys->size].d_info[2]
	    = nst_tok2d(num_toks, 4, toks);     /* S1 */


	itmp = nst_tok2i(num_toks, 5, toks);    /* PID2 */
	if(itmp < 0) return(NORMAL_RC);
	(phys->size)++;
	RC_TRY( realloc_physical_prop_array(phys) );
	phys->array[phys->size].label = itmp;

	phys->array[phys->size].i_info[0] = -1;

	phys->array[phys->size].i_info[1] = NST_PROP_PELAS;

	phys->array[phys->size].d_info[0]
	    = nst_tok2d(num_toks, 6, toks);     /* K2 */

	phys->array[phys->size].d_info[1]
	    = nst_tok2d(num_toks, 7, toks);     /* GE2 */

	phys->array[phys->size].d_info[2]
	    = nst_tok2d(num_toks, 8, toks);     /* S2 */

	return(NORMAL_RC);
}


static RC
put_pelas (FILE *fp, PHYSICAL_PROP phys, DATA_TYPE source)
{
	IDC tmp[5];

	tmp[0].c = "PELAS";
	tmp[1].i = phys.label;           /* PID */
	if(source == FEM_NASTRAN){
		tmp[2].d = phys.d_info[0];   /* K1 */
		tmp[3].d = phys.d_info[1];   /* GE1 */
		tmp[4].d = phys.d_info[2];   /* S1 */
	}else{
		tmp[2].d = FNSDEFAULT;
		tmp[3].d = FNSDEFAULT;
		tmp[4].d = FNSDEFAULT;
	}
	RC_TRY( nst_print(fp, "cidddr", tmp) );

	return(NORMAL_RC);
}


static RC
set_pbeam (PHYSICAL_PROP_ARRAY *phys, char **toks, int num_toks)
{
	int itmp;
	double dtmp;

	itmp = nst_tok2i(num_toks, 1, toks);    /* PID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	phys->array[phys->size].label = itmp;

	phys->array[phys->size].i_info[1] = NST_PROP_PBEAM;

	phys->array[phys->size].i_info[0]
	    = nst_tok2i(num_toks, 2, toks);     /* MID */

	dtmp = nst_tok2d(num_toks, 3, toks);    /* A */
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	phys->array[phys->size].d_info[0] = dtmp;

	dtmp = nst_tok2d(num_toks, 4, toks);    /* I1 */
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	phys->array[phys->size].d_info[1] = dtmp;

	dtmp = nst_tok2d(num_toks, 5, toks);    /* I2 */
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	phys->array[phys->size].d_info[2] = dtmp;

	dtmp = nst_tok2d(num_toks, 6, toks);    /* I12 */
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	phys->array[phys->size].d_info[3] = dtmp;

	dtmp = nst_tok2d(num_toks, 7, toks);    /* J */
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	phys->array[phys->size].d_info[4] = dtmp;

	dtmp = nst_tok2d(num_toks, 8, toks);    /* NSM */
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	phys->array[phys->size].d_info[5] = dtmp;

	dtmp = nst_tok2d(num_toks, 9, toks);    /* C1 */
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	phys->array[phys->size].d_info[6] = dtmp;

	dtmp = nst_tok2d(num_toks, 10, toks);    /* C2 */
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	phys->array[phys->size].d_info[7] = dtmp;

	dtmp = nst_tok2d(num_toks, 11, toks);    /* D1 */
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	phys->array[phys->size].d_info[8] = dtmp;

	dtmp = nst_tok2d(num_toks, 12, toks);    /* D2 */
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	phys->array[phys->size].d_info[9] = dtmp;

	dtmp = nst_tok2d(num_toks, 13, toks);    /* E1 */
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	phys->array[phys->size].d_info[10] = dtmp;

	dtmp = nst_tok2d(num_toks, 14, toks);    /* E2 */
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	phys->array[phys->size].d_info[11] = dtmp;

	dtmp = nst_tok2d(num_toks, 15, toks);    /* F1 */
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	phys->array[phys->size].d_info[12] = dtmp;

	dtmp = nst_tok2d(num_toks, 16, toks);    /* F2 */
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	phys->array[phys->size].d_info[13] = dtmp;

	return(NORMAL_RC);
}


static RC
put_pbeam (FILE *fp, PHYSICAL_PROP phys, DATA_TYPE source)
{
	IDC tmp[9];

	tmp[0].c = "PBEAM";
	tmp[1].i = phys.label;         /* PID */
	tmp[2].i = phys.i_info[0];     /* MID1 */
	tmp[3].d = phys.d_info[0];     /* A */
	tmp[4].d = phys.d_info[1];     /* I1 */
	tmp[5].d = phys.d_info[2];     /* I2 */
	tmp[6].d = phys.d_info[3];     /* I12 */
	tmp[7].d = phys.d_info[4];     /* J */
	tmp[8].d = phys.d_info[5];     /* NSM */
	
	RC_TRY( nst_print(fp, "ciiddddddr", tmp) );

	tmp[0].d = phys.d_info[6];     /* C1 */
	tmp[1].d = phys.d_info[7];     /* C2 */
	tmp[2].d = phys.d_info[8];     /* D1 */
	tmp[3].d = phys.d_info[9];     /* D2 */
	tmp[4].d = phys.d_info[10];    /* E1 */
	tmp[5].d = phys.d_info[11];    /* E2 */
	tmp[6].d = phys.d_info[12];    /* F1 */
	tmp[7].d = phys.d_info[13];    /* F2 */

	RC_TRY( nst_print(fp, "-ddddddddr", tmp) );

	return(NORMAL_RC);
}


RC
nst_input_force (FILE *fp, BC_ARRAY *force, DEFAULT_BC_ARRAY *b_force)
{
	const int num_keys = 4;
	NST_BULK_KEY keys[] = { {"FORCE", 5}, {"MOMENT", 6}, {"GRAV", 4}, 
                            {"SLOAD", 5} };
	int hit_key;
	char *toks[TOKS_SIZE];
	int num_toks;
	NST_FILE nst_file;
	BC_ARRAY moment;

	RC_TRY( nst_allocate_tok(toks, TOKS_SIZE) );
	RC_TRY( init_nst_file(fp, &nst_file) );
	RC_TRY( nst_seek_begin_bulk(&nst_file) );
	RC_TRY( allocate_bc_array(INIT_SIZE, force) );
	RC_TRY( allocate_bc_array(INIT_SIZE, &moment) );
	RC_TRY( allocate_default_bc_array(INIT_SIZE_SMALL, b_force) );
	force->source = FEM_NASTRAN;
	b_force->source = FEM_NASTRAN;

	while(1){
		RC rc = nst_blk_get_tok(&nst_file, num_keys, keys, &hit_key,
		                        TOKS_SIZE, toks, &num_toks);
		if(rc != NORMAL_RC){
			if(rc == END_RC) break;
			return(rc);
		}

		switch(hit_key){
		case 0:
			RC_TRY( realloc_bc_array(force) );
			RC_TRY( set_force(force, toks, num_toks) );
			(force->size)++;
			break;
		case 1:
			RC_TRY( realloc_bc_array(&moment) );
			RC_TRY( set_moment(&moment, toks, num_toks) );
			(moment.size)++;
			break;
		case 2:
			RC_TRY( realloc_default_bc_array(b_force) );
			RC_TRY( set_grav(b_force, toks, num_toks) );
			(b_force->size)++;
			break;
		case 3:
			RC_TRY( realloc_bc_array(force) );
			RC_TRY( set_sload(force, toks, num_toks) );
			(force->size)++;
			break;
		default:
			return(UNKNOWN_ERROR_RC);
		}
	}

	RC_TRY( clean_bc_array(force) );
	RC_TRY( clean_bc_array(&moment) );
	RC_TRY( merge_moment(force, moment) );
	RC_TRY( clean_bc_array(&moment) );
	RC_TRY( clean_default_bc_array(b_force) );
	RC_TRY( term_nst_file(&nst_file) );
	RC_TRY( nst_free_tok(toks, TOKS_SIZE) );

	return(NORMAL_RC);
}


static RC
merge_moment (BC_ARRAY *force, BC_ARRAY moment)
{
	int ii1, ii2;

	for(ii1=0; ii1<moment.size; ii1++){
		RC_TRY( realloc_bc_array(force) );
		force->array[force->size] = moment.array[ii1];
		force->size++;
	}
	RC_TRY( clean_bc_array(force) );

	for(ii1=0; ii1<force->size; ii1++){
		if(force->array[ii1].node < 0) continue;
		for(ii2=ii1+1; ii2<force->size; ii2++){
			if(force->array[ii2].node < 0) continue;
			if(force->array[ii2].node != force->array[ii1].node) break;
			if( (force->array[ii1].node == force->array[ii2].node)
			  &&(force->array[ii1].set_id == force->array[ii2].set_id) ){
				if(force->array[ii1].i_info[2] == MOMENT_FLAG){
					BC tmp = force->array[ii1];
					force->array[ii1] = force->array[ii2];
					force->array[ii2] = tmp;
				}
				if( (force->array[ii1].i_info[2] == FORCE_FLAG)
				  &&(force->array[ii2].i_info[2] == MOMENT_FLAG) ){
					force->array[ii1].i_info[1] = force->array[ii2].i_info[1];
					force->array[ii1].v_type.yz = force->array[ii2].v_type.yz;
					force->array[ii1].v_type.zx = force->array[ii2].v_type.zx;
					force->array[ii1].v_type.xy = force->array[ii2].v_type.xy;
					force->array[ii1].v.yz = force->array[ii2].v.yz;
					force->array[ii1].v.zx = force->array[ii2].v.zx;
					force->array[ii1].v.xy = force->array[ii2].v.xy;
					force->array[ii2].node = -1;
				}
			}
		}
	}
	RC_TRY( clean_bc_array(force) );

	return(NORMAL_RC);
}


static RC
set_force (BC_ARRAY *force, char **toks, int num_toks)
{
	int itmp;
	double dtmp;
	double fact;

	force->array[force->size].i_info[2] = FORCE_FLAG;

	itmp = nst_tok2i(num_toks, 1, toks);    /* SID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	force->array[force->size].set_id = itmp;

	itmp = nst_tok2i(num_toks, 2, toks);    /* G */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	force->array[force->size].node = itmp;

	force->array[force->size].i_info[0]
	     = nst_tok2i(num_toks, 3, toks);    /* CID */

	fact = nst_tok2d(num_toks, 4, toks);    /* F */
	if(fact < CR_FNSDEFAULT) return(CONVERT_ERROR_RC);

	dtmp = nst_tok2d(num_toks, 5, toks);    /* N1 */
	if(dtmp < CR_FNSDEFAULT) return(CONVERT_ERROR_RC);
	if(nearly_eq(dtmp, 0.0)){
		force->array[force->size].v_type.x = BC_FREE;
	}else{
		force->array[force->size].v_type.x = BC_FORCE;
		force->array[force->size].v.x = dtmp * fact;
	}

	dtmp = nst_tok2d(num_toks, 6, toks);    /* N2 */
	if(dtmp < CR_FNSDEFAULT) return(CONVERT_ERROR_RC);
	if(nearly_eq(dtmp, 0.0)){
		force->array[force->size].v_type.y = BC_FREE;
	}else{
		force->array[force->size].v_type.y = BC_FORCE;
		force->array[force->size].v.y = dtmp * fact;
	}

	dtmp = nst_tok2d(num_toks, 7, toks);    /* N3 */
	if(dtmp < CR_FNSDEFAULT) return(CONVERT_ERROR_RC);
	if(nearly_eq(dtmp, 0.0)){
		force->array[force->size].v_type.z = BC_FREE;
	}else{
		force->array[force->size].v_type.z = BC_FORCE;
		force->array[force->size].v.z = dtmp * fact;
	}

	return(NORMAL_RC);
}


static RC
set_moment (BC_ARRAY *force, char **toks, int num_toks)
{
	int itmp;
	double dtmp;
	double fact;

	force->array[force->size].i_info[2] = MOMENT_FLAG;

	itmp = nst_tok2i(num_toks, 1, toks);    /* SID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	force->array[force->size].set_id = itmp;

	itmp = nst_tok2i(num_toks, 2, toks);    /* G */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	force->array[force->size].node = itmp;

	force->array[force->size].i_info[1]
	     = nst_tok2i(num_toks, 3, toks);    /* CID */

	fact = nst_tok2d(num_toks, 4, toks);    /* M */
	if(fact < CR_FNSDEFAULT) return(CONVERT_ERROR_RC);

	dtmp = nst_tok2d(num_toks, 5, toks);    /* N1 */
	if(dtmp < CR_FNSDEFAULT) return(CONVERT_ERROR_RC);
	if(nearly_eq(dtmp, 0.0)){
		force->array[force->size].v_type.yz = BC_FREE;
	}else{
		force->array[force->size].v_type.yz = BC_FORCE;
		force->array[force->size].v.yz = dtmp * fact;
	}

	dtmp = nst_tok2d(num_toks, 6, toks);    /* N2 */
	if(dtmp < CR_FNSDEFAULT) return(CONVERT_ERROR_RC);
	if(nearly_eq(dtmp, 0.0)){
		force->array[force->size].v_type.zx = BC_FREE;
	}else{
		force->array[force->size].v_type.zx = BC_FORCE;
		force->array[force->size].v.zx = dtmp * fact;
	}

	dtmp = nst_tok2d(num_toks, 7, toks);    /* N3 */
	if(dtmp < CR_FNSDEFAULT) return(CONVERT_ERROR_RC);
	if(nearly_eq(dtmp, 0.0)){
		force->array[force->size].v_type.xy = BC_FREE;
	}else{
		force->array[force->size].v_type.xy = BC_FORCE;
		force->array[force->size].v.xy = dtmp * fact;
	}

	return(NORMAL_RC);
}


static RC
set_grav (DEFAULT_BC_ARRAY *b_force, char **toks, int num_toks)
{
	int itmp;
	double dtmp;

	itmp = nst_tok2i(num_toks, 1, toks);    /* SID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	b_force->array[b_force->size].set_id = itmp;

	b_force->array[b_force->size].i_info[0]
	     = nst_tok2i(num_toks, 2, toks);    /* SID */

	dtmp = nst_tok2d(num_toks, 3, toks);    /* A */
	if(dtmp < CR_FNSDEFAULT) return(CONVERT_ERROR_RC);
	b_force->array[b_force->size].d_info[0] = dtmp;

	dtmp = nst_tok2d(num_toks, 4, toks);    /* N1 */
	if(dtmp < CR_FNSDEFAULT) return(CONVERT_ERROR_RC);
	b_force->array[b_force->size].grav.x
	     = dtmp * b_force->array[b_force->size].d_info[0];
	b_force->array[b_force->size].d_info[1] = dtmp;

	dtmp = nst_tok2d(num_toks, 5, toks);    /* N2 */
	if(dtmp < CR_FNSDEFAULT) return(CONVERT_ERROR_RC);
	b_force->array[b_force->size].grav.y
	     = dtmp * b_force->array[b_force->size].d_info[0];
	b_force->array[b_force->size].d_info[2] = dtmp;

	dtmp = nst_tok2d(num_toks, 6, toks);    /* N3 */
	if(dtmp < CR_FNSDEFAULT) return(CONVERT_ERROR_RC);
	b_force->array[b_force->size].grav.z
	     = dtmp * b_force->array[b_force->size].d_info[0];
	b_force->array[b_force->size].d_info[3] = dtmp;
	
	b_force->array[b_force->size].i_info[1]
	     = nst_tok2i(num_toks, 7, toks);    /* MB */

	return(NORMAL_RC);
}


static RC
set_sload (BC_ARRAY *force, char **toks, int num_toks)
{
	int itmp;
    int set_id;
	double dtmp;

	force->array[force->size].i_info[2] = FORCE_FLAG;

	itmp = nst_tok2i(num_toks, 1, toks);    /* SID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
    set_id = itmp;
	force->array[force->size].set_id = set_id;

	itmp = nst_tok2i(num_toks, 2, toks);    /* S1 */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	force->array[force->size].node = itmp;

    dtmp = nst_tok2d(num_toks, 3, toks);    /* F1 */
    if(dtmp < NSDEFAULT) return(CONVERT_ERROR_RC);
    if(nearly_eq(dtmp, 0.0)){
        force->array[force->size].s_type[0] = BC_FREE;
    }else{
        force->array[force->size].s_type[0] = BC_FORCE;
        force->array[force->size].s[0] = dtmp;
    }

    force->size++;

	itmp = nst_tok2i(num_toks, 4, toks);    /* S2 */
	if(itmp < 0) return(NORMAL_RC);
    RC_TRY( realloc_bc_array(force) );
	force->array[force->size].set_id = set_id;
	force->array[force->size].node = itmp;

    dtmp = nst_tok2d(num_toks, 5, toks);    /* F2 */
    if(dtmp < CR_FNSDEFAULT) return(CONVERT_ERROR_RC);
    if(nearly_eq(dtmp, 0.0)){
        force->array[force->size].s_type[0] = BC_FREE;
    }else{
        force->array[force->size].s_type[0] = BC_FORCE;
        force->array[force->size].s[0] = dtmp;
    }

    force->size++;

	itmp = nst_tok2i(num_toks, 6, toks);    /* S3 */
	if(itmp < 0) return(NORMAL_RC);
    RC_TRY( realloc_bc_array(force) );
	force->array[force->size].set_id = set_id;
	force->array[force->size].node = itmp;

    dtmp = nst_tok2d(num_toks, 7, toks);    /* F3 */
    if(dtmp < CR_FNSDEFAULT) return(CONVERT_ERROR_RC);
    if(nearly_eq(dtmp, 0.0)){
        force->array[force->size].s_type[0] = BC_FREE;
    }else{
        force->array[force->size].s_type[0] = BC_FORCE;
        force->array[force->size].s[0] = dtmp;
    }

	return(NORMAL_RC);
}


RC
nst_output_force (FILE *fp, BC_ARRAY force, DEFAULT_BC_ARRAY b_force)
{
	int ii1;
	
	for(ii1=0;ii1<force.size;ii1++){
		if(force.array[ii1].node < 0) continue;

		RC_TRY( put_force(fp, force.array[ii1], force.source) );
		RC_TRY( put_moment(fp, force.array[ii1], force.source) );
	}
	for(ii1=0;ii1<b_force.size;ii1++){
		if(b_force.array[ii1].set_id < 0) continue;

		RC_TRY( put_grav(fp, b_force.array[ii1], b_force.source) );
	}
	RC_TRY( put_sload(fp, force) );

	return(NORMAL_RC);
}


static RC
put_force (FILE *fp, BC force, DATA_TYPE source)
{
	double v_abs;
	IDC tmp[9];
	VECT3D v_tmp;

	tmp[0].c = "FORCE*";
	tmp[1].i = force.set_id;          /* SID */
	tmp[2].i = force.node;            /* G */

	v_tmp.x = force.v.x;
	v_tmp.y = force.v.y;
	v_tmp.z = force.v.z;
	if(force.v_type.x != BC_FORCE) v_tmp.x = 0.0;
	if(force.v_type.y != BC_FORCE) v_tmp.y = 0.0;
	if(force.v_type.z != BC_FORCE) v_tmp.z = 0.0;
	v_abs = abs_vect3d(v_tmp);
	if( nearly_eq(v_abs, 0.0) ) return(NORMAL_RC);

	if(source == FEM_NASTRAN){
		tmp[3].i = force.i_info[0];   /* CID */
	}else{
		tmp[3].i = NSDEFAULT;         /* CID */
	}
	tmp[4].d = v_abs;                 /* F */
	tmp[5].c = "*";
	tmp[6].d = v_tmp.x/v_abs;         /* N1 */
	tmp[7].d = v_tmp.y/v_abs;         /* N2 */
	tmp[8].d = v_tmp.z/v_abs;         /* N3 */

	RC_TRY( nst_print(fp, "cIIIDrcDDDr", tmp) );

	return(NORMAL_RC);
}


static RC
put_moment (FILE *fp, BC force, DATA_TYPE source)
{
	double v_abs;
	IDC tmp[9];
	VECT3D v_tmp;

	tmp[0].c = "MOMENT*";
	tmp[1].i = force.set_id;          /* SID */
	tmp[2].i = force.node;            /* G */

	if(force.v_type.yz == BC_FORCE){
		v_tmp.x = force.v.yz;
	}else{
		v_tmp.x = 0.0;
	}
	if(force.v_type.zx == BC_FORCE){
		v_tmp.y = force.v.zx;
	}else{
		v_tmp.y = 0.0;
	}
	if(force.v_type.xy == BC_FORCE){
		v_tmp.z = force.v.xy;
	}else{
		v_tmp.z = 0.0;
	}
	v_abs = abs_vect3d(v_tmp);
	if( nearly_eq(v_abs, 0.0) ) return(NORMAL_RC);

	if(source == FEM_NASTRAN){
		tmp[3].i = force.i_info[1];   /* CID */
	}else{
		tmp[3].i = NSDEFAULT;         /* CID */
	}
	tmp[4].d = v_abs;                 /* F */
	tmp[5].c = "*";
	tmp[6].d = v_tmp.x/v_abs;         /* N1 */
	tmp[7].d = v_tmp.y/v_abs;         /* N2 */
	tmp[8].d = v_tmp.z/v_abs;         /* N3 */

	RC_TRY( nst_print(fp, "cIIIDrcDDDr", tmp) );

	return(NORMAL_RC);
}


static RC
put_grav (FILE *fp, DEFAULT_BC b_force, DATA_TYPE source)
{
	IDC tmp[9];

	tmp[0].c = "GRAV*";
	tmp[1].i = b_force.set_id;            /* SID */
	if(source == FEM_NASTRAN){
		tmp[2].i = b_force.i_info[0];     /* CID */
		tmp[3].d = b_force.d_info[0];     /* A */
		tmp[4].d = b_force.d_info[1];     /* N1 */
		tmp[5].c = "*";
		tmp[6].d = b_force.d_info[2];     /* N2 */
		tmp[7].d = b_force.d_info[3];     /* N3 */
		tmp[8].i = b_force.i_info[1];     /* MB */
	}else{
		double v_abs = abs_vect3d(b_force.grav);

		if( nearly_eq(v_abs, 0.0) ) return(NORMAL_RC);

		tmp[2].i = NSDEFAULT;             /* CID */
		tmp[3].d = 1.0;                   /* A */
		tmp[4].d = b_force.grav.x/v_abs;  /* N1 */
		tmp[5].c = "*";
		tmp[6].d = b_force.grav.y/v_abs;  /* N2 */
		tmp[7].d = b_force.grav.z/v_abs;  /* N3 */
		tmp[8].i = NSDEFAULT;             /* MB */
	}

	RC_TRY( nst_print(fp, "cIIDDrcDDIr", tmp) );

	return(NORMAL_RC);
}


static RC
put_sload (FILE *fp, BC_ARRAY force)
{
    int index = 0;
    int entry_count;
    IDC tmp[9];

    if(force.size == 0) return (NORMAL_RC);

    while(1){
	    tmp[0].c = "SLOAD";
	    tmp[1].i = force.array[index].set_id;          /* SID */

        /* 1行分の情報をtmpに格納する． */
        for(entry_count = 0; (entry_count < 3) && (index < force.size); 
                index++){
            if(force.array[index].s_type[0] == BC_FORCE){
                if(force.array[index].set_id != tmp[1].i){
                    break;
                }
                tmp[2+2*entry_count].i = force.array[index].node;
                                                       /* S1, S2, S3 */
                tmp[3+2*entry_count].d = force.array[index].s[0];
                                                       /* F1, F2, F3 */
                entry_count++;
            }
        }

        switch(entry_count){
            case 0:
                break;
            case 1:
                RC_TRY( nst_print(fp, "ciidr", tmp) );
                break;
            case 2:
                RC_TRY( nst_print(fp, "ciididr", tmp) );
                break;
            case 3:
                RC_TRY( nst_print(fp, "ciidididr", tmp) );
                break;
            default:
                return(IMPLEMENT_ERROR_RC);
        }

        if(index >= force.size) break;
    }

	return(NORMAL_RC);
}


RC
nst_input_restraint (FILE *fp, BC_ARRAY *rest)
{
	const int num_keys = 3;
	NST_BULK_KEY keys[] = { {"SPC", 3}, {"SPCD", 4}, {"SPC1", 4} };
	int hit_key;
	char *toks[TOKS_SIZE_LARGE];
	int num_toks;
	NST_FILE nst_file;

	RC_TRY( nst_allocate_tok(toks, TOKS_SIZE_LARGE) );
	RC_TRY( init_nst_file(fp, &nst_file) );
	RC_TRY( nst_seek_begin_bulk(&nst_file) );
	RC_TRY( allocate_bc_array(INIT_SIZE, rest) );
	rest->source = FEM_NASTRAN;

	while(1){
		RC rc;
		
		rc = nst_blk_get_tok(&nst_file, num_keys, keys, &hit_key,
		                     TOKS_SIZE_LARGE, toks, &num_toks);
		if(rc != NORMAL_RC){
			if(rc == END_RC) break;
			return(rc);
		}

		switch(hit_key){
		case 0:
			RC_TRY( set_spc(rest, toks, num_toks) );
			break;
		case 1:
			RC_TRY( set_spc(rest, toks, num_toks) );
			break;
		case 2:
			RC_TRY( set_spc1(rest, toks, num_toks) );
			break;
		default:
			return(UNKNOWN_ERROR_RC);
		}
	}

	RC_TRY( clean_bc_array(rest) );
	RC_TRY( term_nst_file(&nst_file) );
	RC_TRY( nst_free_tok(toks, TOKS_SIZE_LARGE) );

	return(NORMAL_RC);
}


#if 0
RC
nst_input_restraint (FILE *fp, BC_ARRAY *rest)
{
	const int num_keys = 3;
	NST_BULK_KEY keys[] = { {"SPC", 3}, {"SPCD", 4}, {"SPC1", 4} };
	int hit_key;
	char *toks[TOKS_SIZE_LARGE];
	int num_toks;
	NST_FILE nst_file;
	BC_ARRAY rest_d;
	int ii1;

	RC_TRY( allocate_bc_array(INIT_SIZE, &rest_d) );

	RC_TRY( nst_allocate_tok(toks, TOKS_SIZE_LARGE) );
	RC_TRY( init_nst_file(fp, &nst_file) );
	RC_TRY( nst_seek_begin_bulk(&nst_file) );
	RC_TRY( allocate_bc_array(INIT_SIZE, rest) );
	rest->source = FEM_NASTRAN;

	while(1){
		RC rc;
		
		rc = nst_blk_get_tok(&nst_file, num_keys, keys, &hit_key,
		                     TOKS_SIZE_LARGE, toks, &num_toks);
		if(rc != NORMAL_RC){
			if(rc == END_RC) break;
			return(rc);
		}

		switch(hit_key){
		case 0:
			RC_TRY( set_spc(rest, toks, num_toks) );
			break;
		case 1:
			RC_TRY( set_spc(&rest_d, toks, num_toks) );
			break;
		case 2:
			RC_TRY( set_spc1(rest, toks, num_toks) );
			break;
		default:
			return(UNKNOWN_ERROR_RC);
		}
	}

	RC_TRY( clean_bc_array(rest) );
	RC_TRY( clean_bc_array(&rest_d) );
	RC_TRY( term_nst_file(&nst_file) );
	RC_TRY( nst_free_tok(toks, TOKS_SIZE_LARGE) );

	for(ii1=0; ii1<rest_d.size; ii1++){
		int index = search_bc_node(*rest, rest_d.array[ii1].node);

		RC_NEG_CHK(index);
		while(rest->array[index].set_id != rest_d.array[ii1].set_id){
			if(rest->array[index].node != rest_d.array[ii1].node){
				return(SEEK_ERROR_RC);
			}
			index++;
			if(index >= rest->size) return(SEEK_ERROR_RC);
		}

		if(rest_d.array[ii1].v_type.x == BC_FIX){
			rest->array[index].v_type.x = BC_FIX;
			rest->array[index].v.x = rest_d.array[ii1].v.x;
		}
		if(rest_d.array[ii1].v_type.y == BC_FIX){
			rest->array[index].v_type.y = BC_FIX;
			rest->array[index].v.y = rest_d.array[ii1].v.y;
		}
		if(rest_d.array[ii1].v_type.z == BC_FIX){
			rest->array[index].v_type.z = BC_FIX;
			rest->array[index].v.z = rest_d.array[ii1].v.z;
		}
		if(rest_d.array[ii1].v_type.xy == BC_FIX){
			rest->array[index].v_type.xy = BC_FIX;
			rest->array[index].v.xy = rest_d.array[ii1].v.xy;
		}
		if(rest_d.array[ii1].v_type.yz == BC_FIX){
			rest->array[index].v_type.yz = BC_FIX;
			rest->array[index].v.yz = rest_d.array[ii1].v.yz;
		}
		if(rest_d.array[ii1].v_type.zx == BC_FIX){
			rest->array[index].v_type.zx = BC_FIX;
			rest->array[index].v.zx = rest_d.array[ii1].v.zx;
		}
	}
	RC_TRY( free_bc_array(&rest_d) );

	return(NORMAL_RC);
}
#endif /* 0 */


RC
nst_output_restraint (FILE *fp, BC_ARRAY rest)
{
	int ii1;
	double di;

	for(ii1=0; ii1<rest.size; ii1++){
		if(rest.array[ii1].node < 0) continue;

		if(chk_spc1(rest.array[ii1])){
			RC_TRY( put_spc1(fp, rest.array[ii1]) );
		}else if(chk_spc_all(rest.array[ii1], &di)){
			RC_TRY( put_spc_all(fp, rest.array[ii1], di) );
		}else{
			RC_TRY( put_spc(fp, rest.array[ii1]) );
		}
	}

	return(NORMAL_RC);
}


static RC
set_spc (BC_ARRAY *rest, char **toks, int num_toks)
{
	int itmp;
	double dtmp;

	RC_TRY( realloc_bc_array(rest) );

	itmp = nst_tok2i(num_toks, 1, toks);    /* SID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	rest->array[rest->size].set_id = itmp;

	itmp = nst_tok2i(num_toks, 2, toks);    /* G1 */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	if(itmp > 100000){
		rest->array[rest->size].node = -1;
		return(NORMAL_RC);
	}
	rest->array[rest->size].node = itmp;

	itmp = nst_tok2i(num_toks, 3, toks);    /* C1 */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	rest->array[rest->size].v_type = nst_bc_type3d(itmp);

	dtmp = nst_tok2d(num_toks, 4, toks);    /* D1 */
	rest->array[rest->size].v
	    = nst_rest_value(rest->array[rest->size].v_type, dtmp);

	(rest->size)++;

	itmp = nst_tok2i(num_toks, 5, toks);    /* G2 */
	if(itmp < 0) return(NORMAL_RC);

	RC_TRY( realloc_bc_array(rest) );
	rest->array[rest->size].set_id = rest->array[rest->size - 1].set_id;
	rest->array[rest->size].node = itmp;

	itmp = nst_tok2i(num_toks, 6, toks);    /* C2 */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	rest->array[rest->size].v_type = nst_bc_type3d(itmp);

	dtmp = nst_tok2d(num_toks, 7, toks);    /* D2 */
	rest->array[rest->size].v
	    = nst_rest_value(rest->array[rest->size].v_type, dtmp);

	(rest->size)++;

	return(NORMAL_RC);
}


static RC
set_spc1 (BC_ARRAY *rest, char **toks, int num_toks)
{
	int sid, c;
	int g = -1;
	int ii1;
	BC_TYPE3D c_type3d;

	sid = nst_tok2i(num_toks, 1, toks);	/* SID */
	if(sid < 0) return(CONVERT_ERROR_RC);

	c = nst_tok2i(num_toks, 2, toks);   /* C */
	if(c < 0) return(CONVERT_ERROR_RC);
	
	c_type3d = nst_bc_type3d(c);

	ii1 = 0;
	while(1){
		if(nst_tok_strncmp(num_toks, 3+ii1, toks, "THRU", 4) == 0){
			int g2;
			int ii2;

			g2 = nst_tok2i(num_toks, 4+ii1, toks);
			if(g < 0) return(CONVERT_ERROR_RC);
			if(g2 < 0) return(CONVERT_ERROR_RC);
			if(g >= g2) return(CONVERT_ERROR_RC);

			for(ii2=g+1; ii2<=g2; ii2++){
				RC_TRY( realloc_bc_array(rest) );

				rest->array[rest->size].set_id = sid;
				rest->array[rest->size].node = ii2;
				rest->array[rest->size].v_type = c_type3d;
				rest->array[rest->size].v
				    = nst_rest_value(rest->array[rest->size].v_type, 0.0);
				(rest->size)++;
			}
			g = -1;
			ii1 += 2;
		}else{
			g = nst_tok2i(num_toks, 3+ii1, toks);
			if(g < 0) break;

			RC_TRY( realloc_bc_array(rest) );

			rest->array[rest->size].set_id = sid;
			rest->array[rest->size].node = g;
			rest->array[rest->size].v_type = c_type3d;
			rest->array[rest->size].v
			    = nst_rest_value(rest->array[rest->size].v_type, 0.0);
			(rest->size)++;
			ii1++;
		}
	}

	return(NORMAL_RC);
}


static int
chk_spc1 (BC rest)
{
	if(rest.v_type.x == BC_FIX){
		if( !nearly_eq(rest.v.x, 0.0) ) return(0);
	}
	if(rest.v_type.y == BC_FIX){
		if( !nearly_eq(rest.v.y, 0.0) ) return(0);
	}
	if(rest.v_type.z == BC_FIX){
		if( !nearly_eq(rest.v.z, 0.0) ) return(0);
	}
	if(rest.v_type.xy == BC_FIX){
		if( !nearly_eq(rest.v.xy, 0.0) ) return(0);
	}
	if(rest.v_type.yz == BC_FIX){
		if( !nearly_eq(rest.v.yz, 0.0) ) return(0);
	}
	if(rest.v_type.zx == BC_FIX){
		if( !nearly_eq(rest.v.zx, 0.0) ) return(0);
	}

	return(1);
}


static int
chk_spc_all (BC rest, double *di)
{
	int hit_flag = 0;
	double val = 0.0;

	if(rest.v_type.x == BC_FIX){
		hit_flag = 1;
		val = rest.v.x;
	}

	if(rest.v_type.y == BC_FIX){
		if(hit_flag){
			if( !nearly_eq(rest.v.y, val) ) return(0);
		}else{
			hit_flag = 1;
			val = rest.v.y;
		}
	}

	if(rest.v_type.z == BC_FIX){
		if(hit_flag){
			if( !nearly_eq(rest.v.z, val) ) return(0);
		}else{
			hit_flag = 1;
			val = rest.v.z;
		}
	}

	if(rest.v_type.yz == BC_FIX){
		if(hit_flag){
			if( !nearly_eq(rest.v.yz, val) ) return(0);
		}else{
			hit_flag = 1;
			val = rest.v.yz;
		}
	}

	if(rest.v_type.zx == BC_FIX){
		if(hit_flag){
			if( !nearly_eq(rest.v.zx, val) ) return(0);
		}else{
			hit_flag = 1;
			val = rest.v.zx;
		}
	}

	if(rest.v_type.xy == BC_FIX){
		if(hit_flag){
			if( !nearly_eq(rest.v.xy, val) ) return(0);
		}else{
			hit_flag = 1;
			val = rest.v.xy;
		}
	}

	*di = val;

	return(1);
}


static RC
put_spc1 (FILE *fp, BC rest)
{
	IDC tmp[4];

	tmp[0].c = "SPC1";
	tmp[1].i = rest.set_id;
	tmp[2].i = bc_type3d_number(rest.v_type);
	if(tmp[2].i <= 0) return(NORMAL_RC);
	tmp[3].i = rest.node;

	RC_TRY( nst_print(fp, "ciiir", tmp) );
	
	return(NORMAL_RC);
}


static RC
put_spc_all (FILE *fp, BC rest, double di)
{
	IDC tmp[5];

	tmp[0].c = "SPC";
	tmp[1].i = rest.set_id;
	tmp[2].i = rest.node;
	tmp[3].i = bc_type3d_number(rest.v_type);
	if(tmp[3].i <= 0) return(NORMAL_RC);
	tmp[4].d = di;

	RC_TRY( nst_print(fp, "ciiidr", tmp) );

	return(NORMAL_RC);
}


static RC
put_spc (FILE *fp, BC rest)
{
	IDC tmp[5];

	tmp[0].c = "SPC";
	tmp[1].i = rest.set_id;
	tmp[2].i = rest.node;

	if(rest.v_type.x == BC_FIX){
		tmp[3].i = 1;
		tmp[4].d = rest.v.x;
		RC_TRY( nst_print(fp, "ciiidr", tmp) );
	}
	if(rest.v_type.y == BC_FIX){
		tmp[3].i = 2;
		tmp[4].d = rest.v.y;
		RC_TRY( nst_print(fp, "ciiidr", tmp) );
	}
	if(rest.v_type.z == BC_FIX){
		tmp[3].i = 3;
		tmp[4].d = rest.v.z;
		RC_TRY( nst_print(fp, "ciiidr", tmp) );
	}
	if(rest.v_type.yz == BC_FIX){
		tmp[3].i = 4;
		tmp[4].d = rest.v.yz;
		RC_TRY( nst_print(fp, "ciiidr", tmp) );
	}
	if(rest.v_type.zx == BC_FIX){
		tmp[3].i = 5;
		tmp[4].d = rest.v.zx;
		RC_TRY( nst_print(fp, "ciiidr", tmp) );
	}
	if(rest.v_type.xy == BC_FIX){
		tmp[3].i = 6;
		tmp[4].d = rest.v.xy;
		RC_TRY( nst_print(fp, "ciiidr", tmp) );
	}

	return(NORMAL_RC);
}


static int
bc_type3d_number (BC_TYPE3D v_type)
{
	int ii1=0;
	
	if(v_type.x == BC_FIX) ii1 = 1;
	if(v_type.y == BC_FIX) ii1 = ii1*10+2;
	if(v_type.z == BC_FIX) ii1 = ii1*10+3;
	if(v_type.yz == BC_FIX) ii1 = ii1*10+4;
	if(v_type.zx == BC_FIX) ii1 = ii1*10+5;
	if(v_type.xy == BC_FIX) ii1 = ii1*10+6;

	return(ii1);
}


static BC_TYPE3D
nst_bc_type3d (int num)
{
	BC_TYPE3D ret;

	ret.x = ret.y = ret.z = ret.xy = ret.yz = ret.zx = BC_FREE;

	while(num != 0){
		switch(num % 10){
		case 1:
			ret.x = BC_FIX;
			break;
		case 2:
			ret.y = BC_FIX;
			break;
		case 3:
			ret.z = BC_FIX;
			break;
		case 4:
			ret.yz = BC_FIX;
			break;
		case 5:
			ret.zx = BC_FIX;
			break;
		case 6:
			ret.xy = BC_FIX;
			break;
		default:
			break;
		}
		num /= 10;
	}

	return(ret);
}


static VECT3DR
nst_rest_value (BC_TYPE3D type3d, double value)
{
	VECT3DR v={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

	if(type3d.x == BC_FIX) v.x = value;
	if(type3d.y == BC_FIX) v.y = value;
	if(type3d.z == BC_FIX) v.z = value;
	if(type3d.yz == BC_FIX) v.yz = value;
	if(type3d.zx == BC_FIX) v.zx = value;
	if(type3d.xy == BC_FIX) v.xy = value;
	
	return(v);
}


RC
nst_input_temperature (FILE *fp, BC_ARRAY *temp, ELEM_BC_ARRAY *elem_temp,
                       DEFAULT_BC_ARRAY *def_temp)
{
	const int num_keys = 3;
	NST_BULK_KEY keys[] = { {"TEMP", 4}, {"TEMPD", 5}, {"TEMPP1", 6} };
	int hit_key;
	char *toks[TOKS_SIZE];
	int num_toks;
	NST_FILE nst_file;

	RC_TRY( nst_allocate_tok(toks, TOKS_SIZE) );
	RC_TRY( init_nst_file(fp, &nst_file) );
	RC_TRY( nst_seek_begin_bulk(&nst_file) );
	RC_TRY( allocate_bc_array(INIT_SIZE, temp) );
	RC_TRY( allocate_elem_bc_array(INIT_SIZE, elem_temp) );
	RC_TRY( allocate_default_bc_array(INIT_SIZE_SMALL, def_temp) );
	temp->source = FEM_NASTRAN;
	elem_temp->source = FEM_NASTRAN;
	def_temp->source = FEM_NASTRAN;

	while(1){
		RC rc = nst_blk_get_tok(&nst_file, num_keys, keys, &hit_key,
		                        TOKS_SIZE, toks, &num_toks);
		if(rc != NORMAL_RC){
			if(rc == END_RC) break;
			return(rc);
		}

		switch(hit_key){
		case 0:
			RC_TRY( set_temp(temp, toks, num_toks) );
			break;
		case 1:
			RC_TRY( set_tempd(def_temp, toks, num_toks) );
			break;
		case 2:
			return(IMPLEMENT_ERROR_RC);
			/* RC_TRY( set_tempp1(elem_temp, toks, num_toks) );*/
			/*break;*/
		default:
			return(UNKNOWN_ERROR_RC);
		}
	}

	RC_TRY( clean_bc_array(temp) );
	RC_TRY( clean_elem_bc_array(elem_temp) );
	RC_TRY( clean_default_bc_array(def_temp) );
	RC_TRY( term_nst_file(&nst_file) );
	RC_TRY( nst_free_tok(toks, TOKS_SIZE) );

	return(NORMAL_RC);
}


static RC
set_temp (BC_ARRAY *temp, char **toks, int num_toks)
{
	int sid;
	int itmp;
	double dtmp;
	int ii1;

	sid = nst_tok2i(num_toks, 1, toks);	           /* SID */
	if(sid < 0) return(CONVERT_ERROR_RC);

	for(ii1=2; ii1<(num_toks-1); ii1+=2){
		RC_TRY( realloc_bc_array(temp) );

		itmp = nst_tok2i(num_toks, ii1, toks);     /* G? */
		if(itmp < 0 || itmp > 100000) continue;

		dtmp = nst_tok2d(num_toks, ii1+1, toks);   /* T? */
		if(dtmp < CR_FNSDEFAULT) continue;

		temp->array[temp->size].node = itmp;
		temp->array[temp->size].set_id = sid;
		temp->array[temp->size].s_type[0] = BC_TEMP;
		temp->array[temp->size].s[0] = dtmp;
		(temp->size)++;
	}

	return(NORMAL_RC);
}


static RC
set_tempd (DEFAULT_BC_ARRAY *def_temp, char **toks, int num_toks)
{
	int itmp;
	double dtmp;
	int ii1;

	for(ii1=1; ii1<(num_toks-1); ii1+=2){
		RC_TRY( realloc_default_bc_array(def_temp) );

		itmp = nst_tok2i(num_toks, ii1, toks);     /* SID? */
		if(itmp < 0) continue;

		dtmp = nst_tok2d(num_toks, ii1+1, toks);   /* T? */
		if(dtmp < CR_FNSDEFAULT) continue;

		def_temp->array[def_temp->size].set_id = itmp;
		def_temp->array[def_temp->size].temp = dtmp;
		(def_temp->size)++;
	}

	return(NORMAL_RC);
}


#if 0
static RC
set_tempp1 (ELEM_BC_ARRAY *elem_temp, char **toks, int num_toks)
{
	int sid;
	int eid;
	double tbar,tprime,T1,T2;
	int ii1;

	sid = nst_tok2i(num_toks, 1, toks);               /* SID */
	if(sid < 0) return(CONVERT_ERROR_RC);

	eid = nst_tok2i(num_toks, 2, toks);              /* EID1 */
	if(eid < 0) return(CONVERT_ERROR_RC);

	tbar = nst_tok2d(num_toks, 3, toks);             /* TBAR */
	if(tbar < CR_FNSDEFAULT) return(CONVERT_ERROR_RC);

	tprime = nst_tok2d(num_toks, 4, toks);           /* TPRIME */
	if(tprime < CR_FNSDEFAULT) return(CONVERT_ERROR_RC);

	T1 = nst_tok2d(num_toks, 5, toks);               /* T1 */

	T2 = nst_tok2d(num_toks, 6, toks);               /* T2 */

	RC_TRY( realloc_elem_bc_array(elem_temp) );
	elem_temp->array[elem_temp->size].set_id = sid;
	elem_temp->array[elem_temp->size].element = eid;
	elem_temp->array[elem_temp->size].s_type[0] = BC_TEMP;
	elem_temp->array[elem_temp->size].s[0] = tbar;
	elem_temp->array[elem_temp->size].s_type[1] = BC_TEMP_GR;
	elem_temp->array[elem_temp->size].s[1] = tprime;
	elem_temp->array[elem_temp->size].d_info[0] = T1;
	elem_temp->array[elem_temp->size].d_info[1] = T2;
	(elem_temp->size)++;

	
	ii1 = 9;
	eid = -1;
	while(1){
		if(nst_tok_strncmp(num_toks, ii1, toks, "THRU", 4) == 0){
			int eid2;
			int ii2;

			eid2 = nst_tok2i(num_toks, ii1+1, toks);
			if(eid < 0) return(CONVERT_ERROR_RC);
			if(eid2 < 0) return(CONVERT_ERROR_RC);
			if(eid >= eid2) return(CONVERT_ERROR_RC);

			for(ii2=eid+1; ii2<=eid2; ii2++){
				RC_TRY( realloc_elem_bc_array(elem_temp) );
				elem_temp->array[elem_temp->size].set_id = sid;
				elem_temp->array[elem_temp->size].element = ii2;
				elem_temp->array[elem_temp->size].s_type[0] = BC_TEMP;
				elem_temp->array[elem_temp->size].s[0] = tbar;
				elem_temp->array[elem_temp->size].s_type[1] = BC_TEMP_GR;
				elem_temp->array[elem_temp->size].s[1] = tprime;
				elem_temp->array[elem_temp->size].d_info[0] = T1;
				elem_temp->array[elem_temp->size].d_info[1] = T2;
				(elem_temp->size)++;
			}
			eid = -1;
			ii1 += 2;
		}else{
			eid = nst_tok2i(num_toks, ii1, toks);
			if(eid < 0) break;

			RC_TRY( realloc_elem_bc_array(elem_temp) );
			elem_temp->array[elem_temp->size].set_id = sid;
			elem_temp->array[elem_temp->size].element = eid;
			elem_temp->array[elem_temp->size].s_type[0] = BC_TEMP;
			elem_temp->array[elem_temp->size].s[0] = tbar;
			elem_temp->array[elem_temp->size].s_type[1] = BC_TEMP_GR;
			elem_temp->array[elem_temp->size].s[1] = tprime;
			elem_temp->array[elem_temp->size].d_info[0] = T1;
			elem_temp->array[elem_temp->size].d_info[1] = T2;
			(elem_temp->size)++;

			ii1++;
		}
	}

	return(NORMAL_RC);
}
#endif


RC
nst_output_temperature (FILE *fp, BC_ARRAY temp, ELEM_BC_ARRAY elem_temp,
                        DEFAULT_BC_ARRAY def_temp)
{
	int ii1;

	for(ii1=0; ii1<temp.size; ii1++){
		if(temp.source == FEM_NASTRAN){
			RC_TRY( put_temp(fp, temp.array[ii1]) );
		}else{
			return(IMPLEMENT_ERROR_RC);
		}
	}
/*	
	for(ii1=0; ii1<elem_temp.size; ii1++){
		if(elem_temp.source == FEM_NASTRAN){
			RC_TRY( put_tempp1(fp, elem_temp.array[ii1]) );
		}else{
			return(IMPLEMENT_ERROR_RC);
		}
	}	
*/
	
	for(ii1=0; ii1<def_temp.size; ii1++){
		if(def_temp.source == FEM_NASTRAN){
			RC_TRY( put_tempd(fp, def_temp.array[ii1]) );
		}else{
			return(IMPLEMENT_ERROR_RC);
		}
	}

	return(NORMAL_RC);
	
}


static RC
put_temp (FILE *fp, BC temp)
{
	IDC tmp[4];

	tmp[0].c = "TEMP";
	tmp[1].i = temp.set_id;         /* SID */
	tmp[2].i = temp.node;           /* G? */
	tmp[3].d = temp.s[0];           /* T? */

	RC_TRY( nst_print(fp, "ciidr", tmp) );

	return(NORMAL_RC);

}	


static RC
put_tempd (FILE *fp, DEFAULT_BC def_temp)
{
	IDC tmp[3];

	tmp[0].c = "TEMPD";
	tmp[1].i = def_temp.set_id;         /* SID? */
	tmp[2].d = def_temp.temp;           /* T? */

	RC_TRY( nst_print(fp, "cidr", tmp) );

	return(NORMAL_RC);
}

#if 0
static RC
put_tempp1 (FILE *fp, ELEM_BC elem_temp)
{
	IDC tmp[7];

	tmp[0].c = "TEMPP1";
	tmp[1].i = elem_temp.set_id;         /* SID */
	tmp[2].i = elem_temp.element;        /* EID1 */
	tmp[3].d = elem_temp.s[0];           /* TBAR */
	tmp[4].d = elem_temp.s[1];           /* TPRIME */
	tmp[5].d = elem_temp.d_info[0];      /* T1 */
	tmp[6].d = elem_temp.d_info[1];      /* T2 */

	RC_TRY( nst_print(fp, "ciiddddr", tmp) );

	return(NORMAL_RC);
}
#endif


RC
nst_input_bc_set (FILE *fp, BC_SET_ARRAY *bc_set)
{
	RC rc;
	const int num_keys = 3;
	NST_BULK_KEY keys[] = {{"SPCADD", 6}, {"LOAD", 4}, {"MPCADD", 6}};
	int hit_key;
	char *toks[TOKS_SIZE];
	int num_toks;
	NST_FILE nst_file;

	RC_TRY( nst_read_case_control(fp, bc_set) );

	RC_TRY( nst_allocate_tok(toks, TOKS_SIZE) );
	RC_TRY( init_nst_file(fp, &nst_file) );
	RC_TRY( nst_seek_begin_bulk(&nst_file) );

	while(1){
		rc = nst_blk_get_tok(&nst_file, num_keys, keys, &hit_key,
		                     TOKS_SIZE, toks, &num_toks);
		if(rc != NORMAL_RC){
			if(rc == END_RC) break;
			return(rc);
		}

		switch(hit_key){
		case 0:
			RC_TRY( set_spcadd(bc_set, toks, num_toks) );
			break;
		case 1:
			RC_TRY( set_load(bc_set, toks, num_toks) );
			break;
		case 2:
			RC_TRY( set_mpcadd(bc_set, toks, num_toks) );
			break;
		default:
			return(UNKNOWN_ERROR_RC);
		}
	}

	RC_TRY( term_nst_file(&nst_file) );
	RC_TRY( nst_free_tok(toks, TOKS_SIZE) );

	return(NORMAL_RC);
}


static RC
set_spcadd (BC_SET_ARRAY *bc_set, char **toks, int num_toks)
{
	int itmp;
	int sid;
	int ii1, ii2;
	int index;

	sid = nst_tok2i(num_toks, 1, toks);    /* SID */
	if(sid < 0) return(CONVERT_ERROR_RC);
	
	for(ii1=0; ii1<(bc_set->size); ii1++){
		if(bc_set->array[ii1].label < 0) continue;

		if(bc_set->array[ii1].id_fix[0] == sid){
			index = 1;
			for(ii2=2; ii2<num_toks; ii2++){
				itmp = nst_tok2i(num_toks, ii2, toks); /* S? */
				if(itmp < 0) break;

				if(index >= MAX_BC_SET) return(OVERFLOW_ERROR_RC);
				bc_set->array[ii1].id_fix[index] = itmp;
				bc_set->array[ii1].scale_fix[index] = 1.0;
				index++;
			}
			bc_set->array[ii1].num_fix = index;
		}
	}

	return(NORMAL_RC);
}


static RC
set_mpcadd (BC_SET_ARRAY *bc_set, char **toks, int num_toks)
{
	int itmp;
	int sid;
	int ii1, ii2;
	int index;

	sid = nst_tok2i(num_toks, 1, toks);    /* SID */
	if(sid < 0) return(CONVERT_ERROR_RC);
	
	for(ii1=0; ii1<(bc_set->size); ii1++){
		if(bc_set->array[ii1].label < 0) continue;

		if(bc_set->array[ii1].id_mpc[0] == sid){
			index = 1;
			for(ii2=2; ii2<num_toks; ii2++){
				itmp = nst_tok2i(num_toks, ii2, toks); /* S? */
				if(itmp < 0) break;

				if(index >= MAX_BC_SET) return(OVERFLOW_ERROR_RC);
				bc_set->array[ii1].id_mpc[index] = itmp;
				bc_set->array[ii1].scale_mpc[index] = 1.0;
				index++;
			}
			bc_set->array[ii1].num_mpc = index;
		}
	}

	return(NORMAL_RC);
}



static RC
set_load (BC_SET_ARRAY *bc_set, char **toks, int num_toks)
{
	int itmp;
	int sid;
	double factor;
	double factor_i;
	int ii1, ii2;
	int index;

	sid = nst_tok2i(num_toks, 1, toks);       /* SID */
	if(sid < 0) return(CONVERT_ERROR_RC);

	factor = nst_tok2d(num_toks, 2, toks);    /* S */
	if(factor < CR_FNSDEFAULT) return(CONVERT_ERROR_RC);

	for(ii1=0; ii1<(bc_set->size); ii1++){
		if(bc_set->array[ii1].label < 0) continue;

		if(bc_set->array[ii1].id_force[0] == sid){
			index = 1;
			for(ii2=3; ii2<(num_toks-1); ii2+=2){
				factor_i = nst_tok2d(num_toks, ii2, toks); /* S? */
				if(factor_i < CR_FNSDEFAULT) break;

				itmp = nst_tok2i(num_toks, ii2+1, toks);   /* L? */
				if(itmp < 0) break;

				if(index >= MAX_BC_SET) return(OVERFLOW_ERROR_RC);
				bc_set->array[ii1].id_force[index] = itmp;
				bc_set->array[ii1].scale_force[index] = factor * factor_i;
				index++;
			}
			bc_set->array[ii1].num_force = index;
		}
	}

	return(NORMAL_RC);
}


RC
nst_write_case_control (FILE *fp, BC_SET_ARRAY bc_set)
{
	int ii1, ii2;
	int com_index;

	if(bc_set.source == FEM_NASTRAN){
		for(ii1=0; ii1<bc_set.s_info.size; ii1++){
			fprintf(fp, "%s\n", bc_set.s_info.str[ii1]);
		}
	}

	com_index = -1;
	for(ii1=0; ii1<bc_set.size; ii1++){
		if(bc_set.array[ii1].label >= 0){
			com_index = ii1;
			break;
		}
	}
	if(com_index < 0) return(ARG_ERROR_RC);
	
	if(bc_set.com_fix == BC_COMMON){
		if(bc_set.array[com_index].num_fix < 1) return(ARG_ERROR_RC);
		fprintf(fp, "SPC = %d\n", bc_set.array[com_index].id_fix[0]);
	}
	if(bc_set.com_mpc == BC_COMMON){
		if(bc_set.array[com_index].num_mpc < 1) return(ARG_ERROR_RC);
		fprintf(fp, "MPC = %d\n", bc_set.array[com_index].id_mpc[0]);
	}
	if(bc_set.com_force == BC_COMMON){
		if(bc_set.array[com_index].num_force < 1) return(ARG_ERROR_RC);
		fprintf(fp, "LOAD = %d\n", bc_set.array[com_index].id_force[0]);
	}
	if(bc_set.com_temp == BC_COMMON){
		if(bc_set.array[com_index].num_temp < 1) return(ARG_ERROR_RC);
		fprintf(fp, "TEMP = %d\n", bc_set.array[com_index].id_temp[0]);
	}

	for(ii1=0; ii1<bc_set.size; ii1++){
		if(bc_set.array[ii1].label < 0) continue;

		fprintf(fp, "SUBCASE %d\n", bc_set.array[ii1].label);
		if(bc_set.source == FEM_NASTRAN){
			for(ii2=0; ii2<bc_set.array[ii1].s_info.size; ii2++){
				fprintf(fp, "    %s\n", bc_set.array[ii1].s_info.str[ii2]);
			}
		}

		if( (bc_set.com_fix != BC_COMMON)
		  &&(bc_set.array[ii1].num_fix >= 1) ){
			fprintf(fp, "    SPC = %d\n", bc_set.array[ii1].id_fix[0]);
		}
		if( (bc_set.com_mpc != BC_COMMON)
		  &&(bc_set.array[ii1].num_mpc >= 1) ){
			fprintf(fp, "    MPC = %d\n", bc_set.array[ii1].id_mpc[0]);
		}
		if( (bc_set.com_force != BC_COMMON)
		  &&(bc_set.array[ii1].num_force >= 1) ){
			fprintf(fp, "    LOAD = %d\n", bc_set.array[ii1].id_force[0]);
		}
		if( (bc_set.com_temp != BC_COMMON)
		  &&(bc_set.array[ii1].num_temp >= 1) ){
			fprintf(fp, "    TEMP = %d\n", bc_set.array[ii1].id_temp[0]);
		}
	}

	fprintf(fp, "BEGIN BULK\n");

	return(NORMAL_RC);
}


RC
nst_output_bc_set (FILE *fp, BC_SET_ARRAY bc_set)
{
	int ii1;
	int com_index;

	com_index = -1;
	for(ii1=0; ii1<bc_set.size; ii1++){
		if(bc_set.array[ii1].label >= 0){
			com_index = ii1;
			break;
		}
	}
	if(com_index < 0) return(ARG_ERROR_RC);
	
	if(bc_set.com_fix == BC_COMMON){
		RC_TRY( put_spcadd(fp, bc_set.array[com_index]) );
	}
	if(bc_set.com_force == BC_COMMON){
		RC_TRY( put_load(fp, bc_set.array[com_index]) );
	}
	if(bc_set.com_mpc == BC_COMMON){
		RC_TRY( put_mpcadd(fp, bc_set.array[com_index]) );
	}

	for(ii1=0; ii1<bc_set.size; ii1++){
		if(bc_set.array[ii1].label < 0) continue;

		if(bc_set.com_fix != BC_COMMON){
			RC_TRY( put_spcadd(fp, bc_set.array[ii1]) );
		}
		if(bc_set.com_force != BC_COMMON){
			RC_TRY( put_load(fp, bc_set.array[ii1]) );
		}
		if(bc_set.com_mpc != BC_COMMON){
			RC_TRY( put_mpcadd(fp, bc_set.array[ii1]) );
		}
	}

	return(NORMAL_RC);
}


static RC
put_spcadd (FILE *fp, BC_SET bc_set)
{
	IDC tmp[2];
	int ii1;
	int r_count;

	if(bc_set.num_fix <= 1) return(NORMAL_RC);

	tmp[0].c = "SPCADD";
	tmp[1].i = bc_set.id_fix[0];  /* SID */
	RC_TRY( nst_print(fp, "ci", tmp) );
	
	r_count = 2;
	for(ii1=1; ii1<bc_set.num_fix; ii1++){
		if(r_count >= 9){
			RC_TRY( nst_print(fp, "r-", tmp) );
			r_count = 1;
		}

		tmp[0].i = bc_set.id_fix[ii1];   /* S? */
		RC_TRY( nst_print(fp, "i", tmp) );
		r_count++;
	}
	RC_TRY( nst_print(fp, "r", tmp) );

	return(NORMAL_RC);
}


static RC
put_mpcadd (FILE *fp, BC_SET bc_set)
{
	IDC tmp[2];
	int ii1;
	int r_count;

	if(bc_set.num_mpc <= 1) return(NORMAL_RC);

	tmp[0].c = "MPCADD";
	tmp[1].i = bc_set.id_mpc[0];  /* SID */
	RC_TRY( nst_print(fp, "ci", tmp) );
	
	r_count = 2;
	for(ii1=1; ii1<bc_set.num_mpc; ii1++){
		if(r_count >= 9){
			RC_TRY( nst_print(fp, "r-", tmp) );
			r_count = 1;
		}

		tmp[0].i = bc_set.id_mpc[ii1];   /* S? */
		RC_TRY( nst_print(fp, "i", tmp) );
		r_count++;
	}
	RC_TRY( nst_print(fp, "r", tmp) );

	return(NORMAL_RC);
}


static RC
put_load (FILE *fp, BC_SET bc_set)
{
	IDC tmp[3];
	int ii1;
	int r_count;

	if(bc_set.num_force <= 1) return(NORMAL_RC);

	tmp[0].c = "LOAD";
	tmp[1].i = bc_set.id_force[0];  /* SID */
	tmp[2].d = 1.0;                 /* S */
	RC_TRY( nst_print(fp, "cid", tmp) );
	
	r_count = 3;
	for(ii1=1; ii1<bc_set.num_force; ii1++){
		if(r_count >= 9){
			RC_TRY( nst_print(fp, "r-", tmp) );
			r_count = 1;
		}

		tmp[0].d = bc_set.scale_force[ii1];   /* S? */
		tmp[1].i = bc_set.id_force[ii1];      /* L? */
		RC_TRY( nst_print(fp, "di", tmp) );
		r_count += 2;
	}
	RC_TRY( nst_print(fp, "r", tmp) );

	return(NORMAL_RC);
}


static RC
nst_read_case_control (FILE *fp, BC_SET_ARRAY *bc_set)
{
	RC rc;
	int itmp;
	char *ptr;
	BC_SET com_bc_set;
	NST_FILE nst_file;
	int hit_key;
	int array_index = -1;
	const int num_keys = 5;
	NST_BULK_KEY keys[] = {{"SUBCASE", 7}, {"SPC", 3},
	                       {"LOAD", 4},    {"TEMP", 4}, {"MPC", 3}};
	STRING_ARRAY *stra_ptr;

	RC_TRY( init_nst_file(fp, &nst_file) );
	RC_TRY( nst_seek_cend(&nst_file) );
	RC_TRY( allocate_bc_set_array(INIT_SIZE_SMALL, bc_set) );
	init_bc_set(&com_bc_set);

	bc_set->source = FEM_NASTRAN;

	stra_ptr = &(bc_set->s_info);
	while(1){
		rc = get_case_ctrl(&nst_file, num_keys, keys,
		                   &hit_key, &ptr, stra_ptr);
		if(rc != NORMAL_RC){
			if(rc == END_RC) break;
			return(rc);
		}

		switch(hit_key){
		case 0:
			if(sscanf(ptr, "%d", &itmp) != 1) return(CONVERT_ERROR_RC);
			if(itmp < 0) return(CONVERT_ERROR_RC);
			RC_TRY( realloc_bc_set_array(bc_set) );
			array_index++;
			(bc_set->size)++;

			bc_set->array[array_index] = com_bc_set;
			bc_set->array[array_index].label = itmp;
			stra_ptr = &(bc_set->array[array_index].s_info);

			break;
		case 1:
			if(sscanf(ptr, "%d", &itmp) != 1) return(CONVERT_ERROR_RC);
			if(itmp < 0) return(CONVERT_ERROR_RC);

			if(array_index < 0){
				bc_set->com_fix = BC_COMMON;
				com_bc_set.num_fix = 1;
				com_bc_set.scale_fix[0] = 1.0;
				com_bc_set.id_fix[0] = itmp;
			}else{
				bc_set->array[array_index].num_fix = 1;
				bc_set->array[array_index].scale_fix[0] = 1.0;
				bc_set->array[array_index].id_fix[0] = itmp;
			}

			break;
		case 2:
			if(sscanf(ptr, "%d", &itmp) != 1) return(CONVERT_ERROR_RC);
			if(itmp < 0) return(CONVERT_ERROR_RC);

			if(array_index < 0){
				bc_set->com_force = BC_COMMON;
				com_bc_set.num_force = 1;
				com_bc_set.scale_force[0] = 1.0;
				com_bc_set.id_force[0] = itmp;
			}else{
				bc_set->array[array_index].num_force = 1;
				bc_set->array[array_index].scale_force[0] = 1.0;
				bc_set->array[array_index].id_force[0] = itmp;
			}

			break;
		case 3:
			if(sscanf(ptr, "%d", &itmp) != 1) return(CONVERT_ERROR_RC);
			if(itmp < 0) return(CONVERT_ERROR_RC);

			if(array_index < 0){
				bc_set->com_temp = BC_COMMON;
				com_bc_set.num_temp = 1;
				com_bc_set.scale_temp[0] = 1.0;
				com_bc_set.id_temp[0] = itmp;
			}else{
				bc_set->array[array_index].num_temp = 1;
				bc_set->array[array_index].scale_temp[0] = 1.0;
				bc_set->array[array_index].id_temp[0] = itmp;
			}

			break;
		case 4:
			if(sscanf(ptr, "%d", &itmp) != 1) return(CONVERT_ERROR_RC);
			if(itmp < 0) return(CONVERT_ERROR_RC);

			if(array_index < 0){
				bc_set->com_mpc = BC_COMMON;
				com_bc_set.num_mpc = 1;
				com_bc_set.scale_mpc[0] = 1.0;
				com_bc_set.id_mpc[0] = itmp;
			}else{
				bc_set->array[array_index].num_mpc = 1;
				bc_set->array[array_index].scale_mpc[0] = 1.0;
				bc_set->array[array_index].id_mpc[0] = itmp;
			}

			break;
		default:
			return(UNKNOWN_ERROR_RC);
		}
	}

	if(bc_set->size <= 0){
		RC_TRY( realloc_bc_set_array(bc_set) );

		bc_set->size = 1;
		bc_set->array[0] = com_bc_set;
		bc_set->array[0].label = 1;
	}

	RC_TRY( clean_bc_set_array(bc_set) );
	RC_TRY( term_nst_file(&nst_file) );

	return(NORMAL_RC);
}


/* keys と一致する行が現れるまで読み込んで、該当行を   */
/* 格納するバッファのデータ部のポインタ buf_ptr に格納 */
static RC
get_case_ctrl (NST_FILE *nst_file, int num_keys, const NST_BULK_KEY *keys,
               int *hit_key, char **buf_ptr, STRING_ARRAY *stra)
{
	int ii1;
	int offset;

	while(1){
		RC_TRY( nst_gets(buf_ptr, nst_file) );
		if(nst_chk_begin_bulk(*buf_ptr)) break;

		for(ii1=0; ii1<num_keys; ii1++){
			if(offset_strncmp(*buf_ptr, keys[ii1].key,
			                  keys[ii1].key_size, &offset) == 0){
				*hit_key = ii1;
				*buf_ptr += offset;
				return(NORMAL_RC);
			}
		}
		if( ((*buf_ptr)[0] != '$')&&((*buf_ptr)[0] != '\0') ){
			char *str_ptr = *buf_ptr;
			while(*str_ptr == ' ') str_ptr++;
			RC_TRY( add_string_array(stra, str_ptr) );
		}
	}

	return(END_RC);
}


/* buf の先端の空白を無視して key と比較する              */
/* 一致していれば、buf の先頭の空白と'='と key の文字数を */
/* offset に代入して 0 を返却                             */
static int
offset_strncmp (const char *buf, const char *key, size_t n, int *offset)
{
	int ret = 1;

	*offset = 0;
	while((int)(buf[*offset]) == ' '){
		(*offset)++;
	}

	ret = strncmp(&(buf[*offset]), key, n);
	if(ret != 0) return(ret);
	
	(*offset) += n;

	/* buf 内の文字列が長ければ不一致とする */
	if( ( ((char)'A' <= buf[*offset])&&(buf[*offset] <= (char)'Z') )
	  ||( ((char)'a' <= buf[*offset])&&(buf[*offset] <= (char)'z') ) ){
		return(1);
	}

	while( ((buf[*offset]) == (char)' ')
	     ||((buf[*offset]) == (char)'=') ){
		(*offset)++;
	}

	return(0);
}


RC
nst_input_rigid_element (FILE *fp, RIGID_ELEMENT_ARRAY *element)
{
	RC rc;
	const int num_keys = 6;
	NST_BULK_KEY keys[] = {{"MPC", 3},   {"RBAR", 4},  {"RBE2", 4},
	                       {"CELAS1", 6},{"CELAS2", 6},{"CONM2", 5} };
	int hit_key;
	char *toks[TOKS_SIZE_LARGE];
	int num_toks;
	NST_FILE nst_file;

	RC_TRY( nst_allocate_tok(toks, TOKS_SIZE_LARGE) );
	RC_TRY( init_nst_file(fp, &nst_file) );
	RC_TRY( nst_seek_begin_bulk(&nst_file) );
	RC_TRY( allocate_rigid_element_array(INIT_SIZE_SMALL, element) );
	element->source = FEM_NASTRAN;

	while(1){
		rc = nst_blk_get_tok(&nst_file, num_keys, keys, &hit_key,
		                     TOKS_SIZE_LARGE, toks, &num_toks);
		if(rc != NORMAL_RC){
			if(rc == END_RC) break;
			return(rc);
		}

		RC_TRY( realloc_rigid_element_array(element) );

		switch(hit_key){
		case 0:
			RC_TRY( set_mpc(element, toks, num_toks) );
			break;
		case 1:
			RC_TRY( set_rbar(element, toks, num_toks) );
			break;
		case 2:
			RC_TRY( set_rbe2(element, toks, num_toks) );
			break;
		case 3:
			RC_TRY( set_celas1(element, toks, num_toks) );
			break;
		case 4:
			RC_TRY( set_celas2(element, toks, num_toks) );
			break;
		case 5:
			RC_TRY( set_conm2(element, toks, num_toks) );
			break;
		default:
			return(UNKNOWN_ERROR_RC);
		}
		(element->size)++;
	}

	RC_TRY( clean_rigid_element_array(element) );
	RC_TRY( term_nst_file(&nst_file) );
	RC_TRY( nst_free_tok(toks, TOKS_SIZE_LARGE) );

	return(NORMAL_RC);
}


static RC
set_mpc (RIGID_ELEMENT_ARRAY *element, char **toks, int num_toks)
{
	int ii1;
	int count;
	int itmp;
	double dtmp;

	/* 節点数をカウント */
	count = 1;
	ii1 = 5;
	while(1){
		itmp = nst_tok2i(num_toks, ii1, toks);  /* G? */
		if(itmp <= 0) break;

		count++;
		if(count%2 == 0){
			ii1 += 5;
		}else{
			ii1 += 3;
		}
	}

	/* 配列確保 */
	element->array[element->size].label = (element->size) + 1;
	element->array[element->size].node_num = count;
	element->array[element->size].node = mm_alloc(count * sizeof(int));
	if(element->array[element->size].node == NULL) return(ALLOC_ERROR_RC);
	element->array[element->size].weight_num = count;
	element->array[element->size].weight = mm_alloc(count * sizeof(double));
	if(element->array[element->size].weight == NULL) return(ALLOC_ERROR_RC);
	element->array[element->size].dof_num = count;
	element->array[element->size].dof = mm_alloc(count * sizeof(int));
	if(element->array[element->size].dof == NULL) return(ALLOC_ERROR_RC);

	itmp = nst_tok2i(num_toks, 1, toks);    /* SID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	element->array[element->size].set_id = itmp;
	element->array[element->size].type = ELEM_MPC;

	itmp = nst_tok2i(num_toks, 2, toks);   /* GM */
	if(itmp <= 0) return(CONVERT_ERROR_RC);
	element->array[element->size].node[0] = itmp;

	itmp = nst_tok2i(num_toks, 3, toks);    /* CM */
	if(itmp <= 0) return(CONVERT_ERROR_RC);
	element->array[element->size].dof[0] = itmp - 1;

	dtmp = nst_tok2d(num_toks, 4, toks);    /* Am */
	if(dtmp < CR_FNSDEFAULT) return(CONVERT_ERROR_RC);
	element->array[element->size].weight[0] = dtmp;

	count = 1;
	ii1 = 5;
	while(1){
		itmp = nst_tok2i(num_toks, ii1, toks);  /* G? */
		if(itmp <= 0) break;
		if(count >= element->array[element->size].node_num){
			return(UNKNOWN_ERROR_RC);
		}
		element->array[element->size].node[count] = itmp;

		itmp = nst_tok2i(num_toks, ii1+1, toks);    /* C? */
		if(itmp <= 0) return(CONVERT_ERROR_RC);
		element->array[element->size].dof[count] = itmp - 1;

		dtmp = nst_tok2d(num_toks, ii1+2, toks);    /* Am */
		if(dtmp < CR_FNSDEFAULT) return(CONVERT_ERROR_RC);
		element->array[element->size].weight[count] = dtmp;

		count++;
		if(count%2 == 0){
			ii1 += 5;
		}else{
			ii1 += 3;
		}
	}
	element->array[element->size].node_num = count;

	return(NORMAL_RC);
}


static RC
set_rbar (RIGID_ELEMENT_ARRAY *element, char **toks, int num_toks)
{
	int cna, cnb;
	int cma, cmb;
	int itmp;

	/* 配列確保 */
	element->array[element->size].node_num = 2;
	element->array[element->size].node = mm_alloc(2 * sizeof(int));
	if(element->array[element->size].node == NULL) return(ALLOC_ERROR_RC);
	element->array[element->size].weight_num = 1;
	element->array[element->size].weight = mm_alloc(1 * sizeof(double));
	if(element->array[element->size].weight == NULL) return(ALLOC_ERROR_RC);
	element->array[element->size].dof_num = 1;
	element->array[element->size].dof = mm_alloc(1 * sizeof(int));
	if(element->array[element->size].dof == NULL) return(ALLOC_ERROR_RC);

	element->array[element->size].set_id = -1;
	element->array[element->size].type = ELEM_RBAR;
	element->array[element->size].i_info[0] = NST_RBAR_RBAR;

	itmp = nst_tok2i(num_toks, 1, toks);   /* EID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	element->array[element->size].label = itmp;

	itmp = nst_tok2i(num_toks, 2, toks);   /* GA */
	if(itmp <= 0) return(CONVERT_ERROR_RC);
	element->array[element->size].node[0] = itmp;

	itmp = nst_tok2i(num_toks, 3, toks);   /* GB */
	if(itmp <= 0) return(CONVERT_ERROR_RC);
	element->array[element->size].node[1] = itmp;

	cna = nst_tok2i(num_toks, 4, toks);   /* CNA */
	cnb = nst_tok2i(num_toks, 5, toks);   /* CNB */
	cma = nst_tok2i(num_toks, 6, toks);   /* CMA */
	cmb = nst_tok2i(num_toks, 7, toks);   /* CMB */

	if( (cna > 0)&&(cnb > 0) ) return(IMPLEMENT_ERROR_RC);
	if( (cna <= 0)&&(cnb <= 0) ) return(CONVERT_ERROR_RC);
	if(cna > 0){
		element->array[element->size].dof[0] = cna;
		element->array[element->size].i_info[1] = cnb;
		element->array[element->size].i_info[2] = cma;
		element->array[element->size].i_info[3] = cmb;
	}else{
		itmp = element->array[element->size].node[0];
		element->array[element->size].node[0]
		       = element->array[element->size].node[1];
		element->array[element->size].node[1] = itmp;
		element->array[element->size].dof[0] = cnb;
		element->array[element->size].i_info[1] = cna;
		element->array[element->size].i_info[2] = cmb;
		element->array[element->size].i_info[3] = cma;
	}

	return(NORMAL_RC);
}


static RC
set_rbe2 (RIGID_ELEMENT_ARRAY *element, char **toks, int num_toks)
{
	int ii1, count;
	int itmp;

	/* 節点数をカウント */
	count = 1;
	while(1){
		itmp = nst_tok2i(num_toks, 3+count, toks);  /* GM? */
		if(itmp <= 0) break;
		count++;
	}

	/* 配列確保 */
	element->array[element->size].node_num = count;
	element->array[element->size].node = mm_alloc(count * sizeof(int));
	if(element->array[element->size].node == NULL) return(ALLOC_ERROR_RC);
	element->array[element->size].weight_num = 1;
	element->array[element->size].weight = mm_alloc(1 * sizeof(double));
	if(element->array[element->size].weight == NULL) return(ALLOC_ERROR_RC);
	element->array[element->size].dof_num = 1;
	element->array[element->size].dof = mm_alloc(1 * sizeof(int));
	if(element->array[element->size].dof == NULL) return(ALLOC_ERROR_RC);

	element->array[element->size].set_id = -1;
	element->array[element->size].type = ELEM_RBAR;
	element->array[element->size].i_info[0] = NST_RBAR_RBE2;

	itmp = nst_tok2i(num_toks, 1, toks);    /* EID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	element->array[element->size].label = itmp;

	itmp = nst_tok2i(num_toks, 2, toks);   /* GN */
	if(itmp <= 0) return(CONVERT_ERROR_RC);
	element->array[element->size].node[0] = itmp;

	itmp = nst_tok2i(num_toks, 3, toks);    /* CM */
	if(itmp <= 0) return(CONVERT_ERROR_RC);
	element->array[element->size].dof[0] = itmp;
	element->array[element->size].weight[0] = 1.0;

	ii1 = 1;
	while(1){
		itmp = nst_tok2i(num_toks, 3+ii1, toks);  /* GM? */
		if(itmp <= 0) break;
		if(ii1 >= element->array[element->size].node_num){
			return(UNKNOWN_ERROR_RC);
		}
		element->array[element->size].node[ii1] = itmp;
		ii1++;
	}
	element->array[element->size].node_num = ii1;

	return(NORMAL_RC);
}


static RC
set_celas1 (RIGID_ELEMENT_ARRAY *element, char **toks, int num_toks)
{
	int count = 0;
	int itmp;
	int chk_count;
	int ii1;

	if(nst_tok2i(num_toks, 3, toks) > 0) count++;
	if(nst_tok2i(num_toks, 5, toks) > 0) count++;
	if(count <= 0) return(CONVERT_ERROR_RC);
	element->array[element->size].node_num = count;
	element->array[element->size].weight_num = count;
	element->array[element->size].dof_num = count;

	element->array[element->size].node = mm_alloc(sizeof(int)*count);
	if(element->array[element->size].node == NULL) return(ALLOC_ERROR_RC);
	element->array[element->size].weight = mm_alloc(sizeof(double)*count);
	if(element->array[element->size].weight == NULL) return(ALLOC_ERROR_RC);
	element->array[element->size].dof = mm_alloc(sizeof(int)*count);
	if(element->array[element->size].dof == NULL) return(ALLOC_ERROR_RC);

	element->array[element->size].set_id = -1;
	element->array[element->size].type = ELEM_SPRING;
	element->array[element->size].i_info[0] = NST_SPRING_CELAS1;

	itmp = nst_tok2i(num_toks, 1, toks);    /* EID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	element->array[element->size].label = itmp;

	itmp = nst_tok2i(num_toks, 2, toks);    /* PID */
	if(itmp < 0) itmp = element->array[element->size].label;
	element->array[element->size].physical = itmp;

	for(ii1=0; ii1<count; ii1++){
		element->array[element->size].weight[ii1] = 0.0;
	}

	chk_count = 0;
	itmp = nst_tok2i(num_toks, 3, toks);    /* G1 */
	if(itmp > 0){
		element->array[element->size].node[chk_count] = itmp;
		itmp = nst_tok2i(num_toks, 4, toks);   /* C1 */
		if(itmp <= 0) return(IMPLEMENT_ERROR_RC);
		element->array[element->size].dof[chk_count] = itmp;
		chk_count++;
	}

	itmp = nst_tok2i(num_toks, 5, toks);   /* G2 */
	if(itmp > 0){
		if(chk_count >= count) return(UNKNOWN_ERROR_RC);
		element->array[element->size].node[chk_count] = itmp;
		itmp = nst_tok2i(num_toks, 6, toks);   /* C2 */
		if(itmp <= 0) return(IMPLEMENT_ERROR_RC);
		element->array[element->size].dof[chk_count] = itmp;
		chk_count++;
	}
	if(chk_count != count) return(UNKNOWN_ERROR_RC);

	return(NORMAL_RC);
}


static RC
set_celas2 (RIGID_ELEMENT_ARRAY *element, char **toks, int num_toks)
{
	int count = 0;
	int itmp;
	int chk_count;
	int ii1;
	double dtmp;

	if(nst_tok2i(num_toks, 3, toks) > 0) count++;
	if(nst_tok2i(num_toks, 5, toks) > 0) count++;
	if(count <= 0) return(CONVERT_ERROR_RC);
	element->array[element->size].node_num = count;
	element->array[element->size].weight_num = count;
	element->array[element->size].dof_num = count;

	element->array[element->size].node = mm_alloc(sizeof(int)*count);
	if(element->array[element->size].node == NULL) return(ALLOC_ERROR_RC);
	element->array[element->size].weight = mm_alloc(sizeof(double)*count);
	if(element->array[element->size].weight == NULL) return(ALLOC_ERROR_RC);
	element->array[element->size].dof = mm_alloc(sizeof(int)*count);
	if(element->array[element->size].dof == NULL) return(ALLOC_ERROR_RC);

	element->array[element->size].set_id = -1;
	element->array[element->size].type = ELEM_SPRING;
	element->array[element->size].i_info[0] = NST_SPRING_CELAS2;

	itmp = nst_tok2i(num_toks, 1, toks);    /* EID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	element->array[element->size].label = itmp;

	dtmp = nst_tok2d(num_toks, 2, toks);    /* K */
	if(dtmp < CR_FNSDEFAULT) return(CONVERT_ERROR_RC);
	for(ii1=0; ii1<count; ii1++){
		element->array[element->size].weight[ii1] = dtmp;
	}

	chk_count = 0;
	itmp = nst_tok2i(num_toks, 3, toks);    /* G1 */
	if(itmp > 0){
		element->array[element->size].node[chk_count] = itmp;
		itmp = nst_tok2i(num_toks, 4, toks);   /* C1 */
		if(itmp <= 0) return(IMPLEMENT_ERROR_RC);
		element->array[element->size].dof[chk_count] = itmp;
		chk_count++;
	}

	itmp = nst_tok2i(num_toks, 5, toks);   /* G2 */
	if(itmp > 0){
		if(chk_count >= count) return(UNKNOWN_ERROR_RC);
		element->array[element->size].node[chk_count] = itmp;
		itmp = nst_tok2i(num_toks, 6, toks);   /* C2 */
		if(itmp <= 0) return(IMPLEMENT_ERROR_RC);
		element->array[element->size].dof[chk_count] = itmp;
		chk_count++;
	}
	if(chk_count != count) return(UNKNOWN_ERROR_RC);
	
	element->array[element->size].d_info[0]
	       = nst_tok2d(num_toks, 7, toks);   /* GE */

	element->array[element->size].d_info[1]
	       = nst_tok2d(num_toks, 8, toks);   /* S */

	return(NORMAL_RC);
}


static RC
put_celas (FILE *fp, RIGID_ELEMENT element, DATA_TYPE source)
{
	IDC tmp[9];

	if(element.node_num <= 0) return(ARG_ERROR_RC);
	if(source != FEM_NASTRAN) return(IMPLEMENT_ERROR_RC);

	if(element.i_info[0] == NST_SPRING_CELAS1){
		tmp[0].c = "CELAS1";
		tmp[1].i = element.label;         /* EID */
		tmp[2].i = element.physical;      /* PID */
		tmp[3].i = element.node[0];       /* G1 */
		tmp[4].i = element.dof[0];        /* C1 */
		if(element.node_num >= 2){
			tmp[5].i = element.node[1];   /* G2 */
			tmp[6].i = element.dof[1];    /* C2 */
		}else{
			tmp[5].i = NSDEFAULT;         /* G2 */
			tmp[6].i = NSDEFAULT;         /* C2 */
		}
		RC_TRY( nst_print(fp, "ciiiiiir", tmp) );
	}else if(element.i_info[0] == NST_SPRING_CELAS2){
		tmp[0].c = "CELAS2";
		tmp[1].i = element.label;         /* EID */
		tmp[2].d = element.weight[0];     /* K */
		tmp[3].i = element.node[0];       /* G1 */
		tmp[4].i = element.dof[0];        /* C1 */
		if(element.node_num >= 2){
			tmp[5].i = element.node[1];   /* G2 */
			tmp[6].i = element.dof[1];    /* C2 */
		}else{
			tmp[5].i = NSDEFAULT;         /* G2 */
			tmp[6].i = NSDEFAULT;         /* C2 */
		}
		tmp[7].d = element.d_info[0];     /* GE */
		tmp[8].d = element.d_info[1];     /* S */
		RC_TRY( nst_print(fp, "cidiiiiddr", tmp) );
	}else{
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
set_conm2 (RIGID_ELEMENT_ARRAY *element, char **toks, int num_toks)
{
	int itmp;
	double dtmp;
	double *mass;
	double m, x1, x2, x3;
	double i11, i21, i22, i31, i32, i33;
	int ii1;

	element->array[element->size].set_id = -1;
	element->array[element->size].type = ELEM_MASS;
	element->array[element->size].i_info[0] = NST_MASS_CONM2;

	element->array[element->size].node_num = 1;
	element->array[element->size].node = mm_alloc(sizeof(int)*1);
	if( element->array[element->size].node ) return(ALLOC_ERROR_RC);

	element->array[element->size].weight_num = 21;
	mass = mm_alloc(sizeof(double)*21);
	if( mass ) return(ALLOC_ERROR_RC);
	element->array[element->size].weight = mass;
	for(ii1=0; ii1<21; ii1++){
		mass[ii1] = 0.0;
	}

	element->array[element->size].dof_num = 0;
	element->array[element->size].dof = NULL;

	itmp = nst_tok2i(num_toks, 1, toks);    /* EID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	element->array[element->size].label = itmp;

	itmp = nst_tok2i(num_toks, 2, toks);    /* G */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	element->array[element->size].node[0] = itmp;

	itmp = nst_tok2i(num_toks, 3, toks);	/* CID */
	if((itmp == -1)||(itmp > 0)) return(IMPLEMENT_ERROR_RC);
	element->array[element->size].i_info[1] = itmp;

	dtmp = nst_tok2d(num_toks, 4, toks);    /* M */
	if(dtmp < CR_FNSDEFAULT) return(CONVERT_ERROR_RC);
	element->array[element->size].d_info[0] = m = dtmp;

	dtmp = nst_tok2d(num_toks, 5, toks);    /* X1 */
	element->array[element->size].d_info[1] = dtmp;
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	x1 = dtmp;

	dtmp = nst_tok2d(num_toks, 6, toks);    /* X2 */
	element->array[element->size].d_info[2] = dtmp;
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	x2 = dtmp;

	dtmp = nst_tok2d(num_toks, 7, toks);    /* X3 */
	element->array[element->size].d_info[3] = dtmp;
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	x3 = dtmp;

	dtmp = nst_tok2d(num_toks, 9, toks);    /* I11 */
	element->array[element->size].d_info[4] = dtmp;
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	i11 = dtmp;

	dtmp = nst_tok2d(num_toks, 10, toks);   /* I21 */
	element->array[element->size].d_info[5] = dtmp;
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	i21 = dtmp;

	dtmp = nst_tok2d(num_toks, 11, toks);   /* I22 */
	element->array[element->size].d_info[6] = dtmp;
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	i22 = dtmp;

	dtmp = nst_tok2d(num_toks, 12, toks);   /* I31 */
	element->array[element->size].d_info[7] = dtmp;
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	i31 = dtmp;

	dtmp = nst_tok2d(num_toks, 13, toks);   /* I32 */
	element->array[element->size].d_info[8] = dtmp;
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	i32 = dtmp;

	dtmp = nst_tok2d(num_toks, 14, toks);   /* I33 */
	element->array[element->size].d_info[9] = dtmp;
	if(dtmp < CR_FNSDEFAULT) dtmp = 0.0;
	i33 = dtmp;

	mass[0] = mass[2] = mass[5] = m;
	mass[9]  = m*(x2*x2 + x3*x3) + i11;
	mass[14] = m*(x1*x1 + x3*x3) + i22;
	mass[20] = m*(x1*x1 + x2*x2) + i33;
	mass[13] = -(m*(x1*x2) + i21);
	mass[18] = -(m*(x1*x3) + i31);
	mass[19] = -(m*(x2*x3) + i32);

	return(NORMAL_RC);
}


static RC
put_conm2 (FILE *fp, RIGID_ELEMENT element, DATA_TYPE source)
{
	IDC tmp[9];

	if(element.node_num != 1) return(ARG_ERROR_RC);
	if( (source != FEM_NASTRAN) ||(element.i_info[0] != NST_MASS_CONM2) ){
		return(IMPLEMENT_ERROR_RC);
	}

	tmp[0].c = "CONM2";
	tmp[1].i = element.label;         /* EID */
	tmp[2].i = element.node[0];       /* G */
	tmp[3].i = element.i_info[1];     /* CID */
	tmp[4].d = element.d_info[0];     /* M */
	tmp[5].d = element.d_info[1];     /* X1 */
	tmp[6].d = element.d_info[2];     /* X2 */
	tmp[7].d = element.d_info[3];     /* X3 */

	RC_TRY( nst_print(fp, "ciiiddddr", tmp) );

	if( (element.d_info[4] < CR_FNSDEFAULT)
	  &&(element.d_info[5] < CR_FNSDEFAULT)
	  &&(element.d_info[6] < CR_FNSDEFAULT) ){
		return(NORMAL_RC);
	}
	tmp[0].d = element.d_info[4];     /* I11 */
	tmp[1].d = element.d_info[5];     /* I21 */
	tmp[2].d = element.d_info[6];     /* I22 */
	tmp[3].d = element.d_info[7];     /* I31 */
	tmp[4].d = element.d_info[8];     /* I32 */
	tmp[5].d = element.d_info[9];     /* I33 */

	RC_TRY( nst_print(fp, "-ddddddr", tmp) );

	return(NORMAL_RC);
}


#if 0
static RC
set_rbe3 (RIGID_ELEMENT_ARRAY *element, char **toks, int num_toks)
{
	int ii1;
	int itmp;
	double weight;
	int compo;
	int node_count;
	int um_mode;

	itmp = nst_tok2i(num_toks, 1, toks);	/* EID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	element->array[element->size].label = itmp;

	element->array[element->size].type = ELEM_RFORCE;

	itmp = nst_tok2i(num_toks, 3, toks);   /* REFGRID */
	if(itmp <= 0) return(CONVERT_ERROR_RC);
	element->array[element->size].node[0] = itmp;

	element->array[element->size].i_info[0]
	    = nst_tok2i(num_toks, 4, toks);    /* CM */

	weight = 0.0;
	compo = 0;
	node_count = 1;
	um_mode = 0;
	for(ii1=5; ii1<num_toks; ii1++){
		if(exist_period(toks[ii1])){
			weight = nst_tok2d(num_toks, ii1, toks);   /* WT? */
			compo = nst_tok2i(num_toks, ii1+1, toks);  /* C? */
			ii1++;
			continue;
		}

		if(nst_tok_strncmp(num_toks, ii1, toks, "UM", 2) == 0){
			um_mode = 1;
			element->array[element->size].i_info[RIGID_MAX_NODE]
			                                           = node_count;
			continue;
		}

		itmp = nst_tok2i(num_toks, ii1, toks);    /* G?,? or GM? */
		if(itmp < 0) continue;
		if(node_count >= RIGID_MAX_NODE){
			fprintf(stderr, "RIGID_MAX_NODE < \n");
			return(OVERFLOW_ERROR_RC);
		}
		element->array[element->size].node[node_count] = itmp;

		if(um_mode){
			compo = nst_tok2i(num_toks, ii1+1, toks);    /* CM? */
			ii1++;
			element->array[element->size].d_info[node_count] = 1.0;
			element->array[element->size].i_info[node_count] = compo;
		}else{
			element->array[element->size].d_info[node_count] = weight;
			element->array[element->size].i_info[node_count] = compo;
		}
		node_count++;
	}
	element->array[element->size].node_num = node_count;
	if(um_mode == 0){
		element->array[element->size].i_info[RIGID_MAX_NODE]
		                                           = node_count;
	}

	return(NORMAL_RC);
}
#endif /* 0 */

#if 0
static int exist_period(const char *str)
{
	int ii1 = 0;

	while(str[ii1] != (char)'\0'){
		if(str[ii1] == (char)'.'){
			return(1);
		}
		ii1++;
	}

	return(0);
}
#endif /* 0 */

RC
nst_output_rigid_element (FILE *fp, RIGID_ELEMENT_ARRAY element)
{
	int ii1;

	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		switch(element.array[ii1].type){
		case ELEM_MPC:
			RC_TRY( put_mpc(fp, element.array[ii1], element.source) );
			break;
		case ELEM_RBAR:
			RC_TRY( put_rbar(fp, element.array[ii1], element.source) );
			break;
		case ELEM_SPRING:
			RC_TRY( put_celas(fp, element.array[ii1], element.source) );
			break;
		case ELEM_MASS:
			RC_TRY( put_conm2(fp, element.array[ii1], element.source) );
			break;
		default:
			return(IMPLEMENT_ERROR_RC);
		}
	}

	return(NORMAL_RC);
}


static RC
put_mpc (FILE *fp, RIGID_ELEMENT element, DATA_TYPE source)
{
	IDC tmp[5];
	int ii1;

	if(element.node_num <= 1) return(ARG_ERROR_RC);

	tmp[0].c = "MPC";
	tmp[1].i = element.set_id;     /* SID */
	tmp[2].i = element.node[0];    /* GM */
	tmp[3].i = element.dof[0] + 1; /* CM */
	tmp[4].d = element.weight[0];  /* Am */
	RC_TRY( nst_print(fp, "ciiid", tmp) );
	
	for(ii1=1; ii1<element.node_num; ii1++){
		if(ii1%2 == 0){
			RC_TRY( nst_print(fp, "-r--", tmp) );
		}

		tmp[0].i = element.node[ii1];    /* GM */
		tmp[1].i = element.dof[ii1] + 1; /* CM */
		tmp[2].d = element.weight[ii1];  /* Am */
		RC_TRY( nst_print(fp, "iid", tmp) );
	}
	RC_TRY( nst_print(fp, "r", tmp) );

	return(NORMAL_RC);
}


static RC
put_rbar (FILE *fp, RIGID_ELEMENT element, DATA_TYPE source)
{
	IDC tmp[8];
	int ii1;
	int r_count;

	if(source != FEM_NASTRAN) return(IMPLEMENT_ERROR_RC);
	if(element.i_info[0] == NST_RBAR_RBAR){
		tmp[0].c = "RBAR";
		tmp[1].i = element.label;      /* EID */
		tmp[2].i = element.node[0];    /* GA */
		tmp[3].i = element.node[1];    /* GB */
		tmp[4].i = element.dof[0];     /* CNA */
		tmp[5].i = element.i_info[1];  /* CNB */
		tmp[6].i = element.i_info[2];  /* CMA */
		tmp[7].i = element.i_info[3];  /* CMB */
		RC_TRY( nst_print(fp, "ciiiiiiir", tmp) );
	}else if(element.i_info[0] == NST_RBAR_RBE2){
		tmp[0].c = "RBE2";
		tmp[1].i = element.label;      /* EID */
		tmp[2].i = element.node[0];    /* GN */
		tmp[3].i = element.dof[0];     /* CM */
		RC_TRY( nst_print(fp, "ciii", tmp) );
	
		r_count = 4;
		for(ii1=1; ii1<element.node_num; ii1++){
			if(r_count >= 9){
				RC_TRY( nst_print(fp, "r-", tmp) );
				r_count = 1;
			}
	
			tmp[0].i = element.node[ii1];
			RC_TRY( nst_print(fp, "i", tmp) );
			r_count++;
		}
		RC_TRY( nst_print(fp, "r", tmp) );
	}else{
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


#if 0
static RC
put_rbe3(FILE *fp, const RIGID_ELEMENT *element, DATA_TYPE source)
{
	IDC tmp[4];
	int ii1, ii_max;
	int r_count;
	int c, old_c;
	double w, old_w;

	tmp[0].c = "RBE3";
	tmp[1].i = element->label;     /* EID */
	tmp[2].i = element->node[0];   /* REFGRID */
	if(source == FEM_NASTRAN){
		tmp[3].i = element->i_info[0]; /* REFC */
	}else{
		tmp[3].i = NSDEFAULT;
	}
	RC_TRY( nst_print(fp, "ci-ii", tmp) );
	
	r_count = 5;
	ii_max = element->node_num;
	if(source == FEM_NASTRAN){
		ii_max = element->i_info[RIGID_MAX_NODE];
	}
	c = old_c = 0;
	w = old_w = 0.0;
	for(ii1=1; ii1<ii_max; ii1++){
		if(source == FEM_NASTRAN){
			c = element->i_info[ii1];
			w = element->d_info[ii1];
		}else{
			c = NSDEFAULT;
			w = 1.0;
		}
		if( (ii1 == 1) || (c != old_c) || !(nearly_eq(w, old_w)) ){
			if(r_count >= 9){
				RC_TRY( nst_print(fp, "r-", tmp) );
				r_count = 1;
			}
			tmp[0].d = w;
			RC_TRY( nst_print(fp, "d", tmp) );
			r_count++;

			if(r_count >= 9){
				RC_TRY( nst_print(fp, "r-", tmp) );
				r_count = 1;
			}
			tmp[0].i = c;
			RC_TRY( nst_print(fp, "i", tmp) );
			r_count++;
		}

		if(r_count >= 9){
			RC_TRY( nst_print(fp, "r-", tmp) );
			r_count = 1;
		}
		tmp[0].i = element->node[ii1];
		RC_TRY( nst_print(fp, "i", tmp) );
		r_count++;

		old_c = c;
		old_w = w;
	}
	RC_TRY( nst_print(fp, "r", tmp) );

	if(ii_max == element->node_num) return(NORMAL_RC);
	

	tmp[0].c = "UM";
	RC_TRY( nst_print(fp, "-c", tmp) );
	
	r_count = 2;
	for(ii1=ii_max; ii1<element->node_num; ii1++){
		if(source == FEM_NASTRAN){
			c = element->i_info[ii1];
		}else{
			c = NSDEFAULT;
		}
		if(r_count >= 8){
			RC_TRY( nst_print(fp, "r--", tmp) );
			r_count = 2;
		}
		tmp[0].i = element->node[ii1];
		tmp[1].i = element->i_info[ii1];
		RC_TRY( nst_print(fp, "ii", tmp) );
		r_count += 2;
	}
	RC_TRY( nst_print(fp, "r", tmp) );

	return(NORMAL_RC);
}
#endif /* 0 */


RC
nst_write_exective_control (FILE *fp, NST_PARAM_ARRAY param)
{
	int ii1;

	for(ii1=0; ii1<param.exective.size; ii1++){
		fprintf(fp, "%s\n", param.exective.str[ii1]);
	}
	fprintf(fp, "CEND\n");

	return(NORMAL_RC);
}


RC
nst_input_param (FILE *fp, NST_PARAM_ARRAY *param)
{
	const int num_keys = 1;
	NST_BULK_KEY keys[] = {{"PARAM", 5}};
	int hit_key;
	char *toks[TOKS_SIZE];
	int num_toks;
	NST_FILE nst_file;
	int alloc_size;
	char *buf;

	RC_TRY( nst_allocate_tok(toks, TOKS_SIZE) );
	RC_TRY( init_nst_file(fp, &nst_file) );

	/* Exective Control Section */
	RC_TRY( init_string_array(&(param->exective)) );
	while(1){
		RC_TRY( nst_gets(&buf, &nst_file) );
		if( (buf[0] == '$')||(buf[0] == '\0') ) continue;
		if( nst_chk_cend(buf) ) break;
		RC_TRY( add_string_array(&(param->exective), buf) );
	}

	RC_TRY( nst_seek_begin_bulk(&nst_file) );

	param->array = (NST_PARAM *)mm_alloc(5*sizeof(NST_PARAM));
	if((param->array) == NULL) return(ALLOC_ERROR_RC);
	alloc_size = 5;
	param->size = 0;

	/* PARAM in Bulk Data Section */
	while(1){
		RC rc = nst_blk_get_tok(&nst_file, num_keys, keys, &hit_key,
		                        TOKS_SIZE, toks, &num_toks);
		if(rc != NORMAL_RC){
			if(rc == END_RC) break;
			return(rc);
		}

		if( alloc_size <= (param->size) ){
			param->array = (NST_PARAM *)mm_realloc( (void *)(param->array),
			               (alloc_size + 5)*sizeof(NST_PARAM) );
			if((param->array) == NULL) return(ALLOC_ERROR_RC);
			alloc_size += 5;
		}

		if(hit_key != 0) return(UNKNOWN_ERROR_RC);
		RC_TRY( set_param(param, toks, num_toks) );
		(param->size)++;
	}

	RC_TRY( term_nst_file(&nst_file) );
	RC_TRY( nst_free_tok(toks, TOKS_SIZE_LARGE) );

	return(NORMAL_RC);
}


static RC
set_param (NST_PARAM_ARRAY *param, char **toks, int num_toks)
{
	if(num_toks < 3) return(ARG_ERROR_RC);

	strncpy(param->array[param->size].name, toks[1], MAX_IDC_S_SIZE);
	param->array[param->size].name[MAX_IDC_S_SIZE] = '\0';

	strncpy(param->array[param->size].val_1.s, toks[2], MAX_IDC_S_SIZE);
	param->array[param->size].val_1.s[MAX_IDC_S_SIZE] = '\0';

	if(num_toks >= 4){
		strncpy(param->array[param->size].val_2.s, toks[3], MAX_IDC_S_SIZE);
		param->array[param->size].val_2.s[MAX_IDC_S_SIZE] = '\0';
	}else{
		param->array[param->size].val_2.s[0] = '\0';
	}

	return(NORMAL_RC);
}


RC
nst_output_param (FILE *fp, NST_PARAM_ARRAY param)
{
	int ii1;
	IDC tmp[4];

	for(ii1=0; ii1<param.size; ii1++){
		tmp[0].c = "PARAM";
		tmp[1].c = param.array[ii1].name;
		tmp[2] = param.array[ii1].val_1;
		tmp[3] = param.array[ii1].val_2;
		nst_print(fp, "ccssr", tmp);
	}

	return(NORMAL_RC);
}


RC
nst_input_NL_param (FILE *fp, NST_NL_PARAM *param)
{
	const int num_keys = 1;
	NST_BULK_KEY keys[] = {{"NLPARM", 6}};
	int hit_key;
	char *toks[TOKS_SIZE];
	int num_toks;
	NST_FILE nst_file;
	int itmp;
	double dtmp;

	RC_TRY( nst_allocate_tok(toks, TOKS_SIZE) );
	RC_TRY( init_nst_file(fp, &nst_file) );
	RC_TRY( nst_seek_begin_bulk(&nst_file) );

	while(1){
		RC rc = nst_blk_get_tok(&nst_file, num_keys, keys, &hit_key,
		                        TOKS_SIZE, toks, &num_toks);
		if(rc != NORMAL_RC){
			if(rc == END_RC) break;
			return(rc);
		}

		/* NLPARAM 読み込み */
		if(hit_key == 0){
			itmp = nst_tok2i(num_toks, 1, toks);
			if(itmp < 0) return(CONVERT_ERROR_RC);

			param->ninc = nst_tok2i(num_toks, 2, toks); /* NINC */

			/* convergence conclusion parameters */
			strncpy(param->conv, toks[7], 3);
			param->conv[3] = '\0';

			/* input of convergence parameters */
			dtmp = nst_tok2d(num_toks, 9, toks);    /* EPSU */
			if(dtmp > CR_FNSDEFAULT){
				param->eps[0] = dtmp;
			}else{
				param->eps[0] = 1.0E-3;
			}

			dtmp = nst_tok2d(num_toks, 10, toks);   /* EPSP */
			if(dtmp > CR_FNSDEFAULT){
				param->eps[1] = dtmp;
			}else{
				param->eps[1] = 1.0E-3;
			}

			dtmp = nst_tok2d(num_toks, 11, toks);   /* EPSW */
			if(dtmp > CR_FNSDEFAULT){
				param->eps[2] = dtmp;
			}else{
				param->eps[2] = 1.0E-5;
			}
		}else{
			return(UNKNOWN_ERROR_RC);
		}
	}

	RC_TRY( term_nst_file(&nst_file) );
	RC_TRY( nst_free_tok(toks, TOKS_SIZE_LARGE) );

	return(NORMAL_RC);
}


RC
nst_input_local_coord (FILE *fp, LOCAL_COORD_ARRAY *coord)
{
	const int num_keys = 1;
	NST_BULK_KEY keys[] = {{"CORD2R", 6}};
	int hit_key;
	char *toks[TOKS_SIZE];
	int num_toks;
	NST_FILE nst_file;

	RC_TRY( nst_allocate_tok(toks, TOKS_SIZE) );
	RC_TRY( init_nst_file(fp, &nst_file) );
	RC_TRY( nst_seek_begin_bulk(&nst_file) );
	RC_TRY( allocate_local_coord_array(INIT_SIZE_SMALL, coord) );
	coord->source = FEM_NASTRAN;

	while(1){
		RC rc = nst_blk_get_tok(&nst_file, num_keys, keys, &hit_key,
		                        TOKS_SIZE, toks, &num_toks);
		if(rc != NORMAL_RC){
			if(rc == END_RC) break;
			return(rc);
		}

		RC_TRY( realloc_local_coord_array(coord) );

		switch(hit_key){
		case 0:
			RC_TRY( set_cord2r(coord, toks, num_toks) );
			break;
		default:
			return(UNKNOWN_ERROR_RC);
		}
		(coord->size)++;
	}

	RC_TRY( clean_local_coord_array(coord) );
	RC_TRY( term_nst_file(&nst_file) );
	RC_TRY( nst_free_tok(toks, TOKS_SIZE) );

	return(NORMAL_RC);
}


static RC
set_cord2r (LOCAL_COORD_ARRAY *coord, char **toks, int num_toks)
{
	int ii1;
	int itmp;
	double dtmp;
	VECT3D pa, pb, pc;
	VECT3D vtmp;

	itmp = nst_tok2i(num_toks, 1, toks);    /* CID */
	if(itmp < 0) return(CONVERT_ERROR_RC);
	coord->array[coord->size].label = itmp;

	itmp = nst_tok2i(num_toks, 2, toks);    /* RID */
	coord->array[coord->size].i_info[0] = itmp;
	coord->array[coord->size].i_info[1] = NST_LC_2R;

	for(ii1=0; ii1<9; ii1++){
		dtmp = nst_tok2d(num_toks, ii1+3, toks);  /* A?, B?, C? */
		coord->array[coord->size].d_info[ii1] = dtmp;
	}

	if(coord->array[coord->size].i_info[0] > 0){
		/* 参照座標系が別のローカル座標系であってはならない． */
		return(IMPLEMENT_ERROR_RC);
	}
	pa = array2vect3d(&(coord->array[coord->size].d_info[0]));
	pb = array2vect3d(&(coord->array[coord->size].d_info[3]));
	pc = array2vect3d(&(coord->array[coord->size].d_info[6]));

	coord->array[coord->size].origin = pa;

	vtmp = sub_vect3d(pb, pa);
	vtmp = mul_scalar_vect3d(1.0/(ABS_ERROR+abs_vect3d(vtmp)), vtmp);
	coord->array[coord->size].d_cos[2] = vtmp;

	vtmp = outer_product3d( sub_vect3d(pb, pa), sub_vect3d(pc, pa) );
	vtmp = mul_scalar_vect3d(1.0/(ABS_ERROR+abs_vect3d(vtmp)), vtmp);
	coord->array[coord->size].d_cos[1] = vtmp;

	vtmp = outer_product3d( coord->array[coord->size].d_cos[1],
	                        coord->array[coord->size].d_cos[2] );
	vtmp = mul_scalar_vect3d(1.0/(ABS_ERROR+abs_vect3d(vtmp)), vtmp);
	coord->array[coord->size].d_cos[0] = vtmp;

	return(NORMAL_RC);
}


RC
nst_output_local_coord (FILE *fp, LOCAL_COORD_ARRAY coord)
{
	int ii1;

	if(coord.source != FEM_NASTRAN) return(IMPLEMENT_ERROR_RC);
	for(ii1=0; ii1<coord.size; ii1++){
		if(coord.array[ii1].label < 0) continue;

		switch(coord.array[ii1].i_info[1]){
		case NST_LC_2R:
			RC_TRY( put_cord2r(fp, coord.array[ii1]) );
			break;
		default:
			return(IMPLEMENT_ERROR_RC);
		}
	}

	return(NORMAL_RC);
}


static RC
put_cord2r (FILE *fp, LOCAL_COORD coord)
{
	IDC tmp[12];

	tmp[0].c  = "CORD2R";
	tmp[1].i  = coord.label;       /* CID */
	tmp[2].i  = coord.i_info[0];   /* RID */
	tmp[3].d  = coord.d_info[0];   /* A1 */
	tmp[4].d  = coord.d_info[1];   /* A2 */
	tmp[5].d  = coord.d_info[2];   /* A3 */
	tmp[6].d  = coord.d_info[3];   /* B1 */
	tmp[7].d  = coord.d_info[4];   /* B2 */
	tmp[8].d  = coord.d_info[5];   /* B3 */
	tmp[9].d  = coord.d_info[6];   /* C1 */
	tmp[10].d = coord.d_info[7];   /* C2 */
	tmp[11].d = coord.d_info[8];   /* C3 */
	nst_print(fp, "ciiddddddr-dddr", tmp);

	return(NORMAL_RC);
}


RC
nst_input_pressure (FILE *fp, ELEM_BC_ARRAY *press, ELEMENT_ARRAY element)
{
	const int num_keys = 1;
	NST_BULK_KEY keys[] = { {"PLOAD4", 6} };
	int hit_key;
	char *toks[TOKS_SIZE_LARGE];
	int num_toks;
	NST_FILE nst_file;

	RC_TRY( nst_allocate_tok(toks, TOKS_SIZE_LARGE) );
	RC_TRY( init_nst_file(fp, &nst_file) );
	RC_TRY( nst_seek_begin_bulk(&nst_file) );
	RC_TRY( allocate_elem_bc_array(INIT_SIZE, press) );
	press->source = FEM_NASTRAN;

	while(1){
		RC rc = nst_blk_get_tok(&nst_file, num_keys, keys, &hit_key,
		                        TOKS_SIZE_LARGE, toks, &num_toks);
		if(rc != NORMAL_RC){
			if(rc == END_RC) break;
			return(rc);
		}

		switch(hit_key){
		case 0:
			RC_TRY( set_pload4(press, toks, num_toks, element) );
			break;
		default:
			return(UNKNOWN_ERROR_RC);
		}
	}

	RC_TRY( clean_elem_bc_array(press) );
	RC_TRY( term_nst_file(&nst_file) );
	RC_TRY( nst_free_tok(toks, TOKS_SIZE_LARGE) );

	return(NORMAL_RC);
}


static RC
set_pload4 (ELEM_BC_ARRAY *press, char **toks, int num_toks,
            ELEMENT_ARRAY element)
{
	int itmp;
	int index;
	int g1, g34;
	int tri_flag = 0;

	RC_TRY( realloc_elem_bc_array(press) );

	itmp = nst_tok2i(num_toks, 1, toks);   /* SID */
	if(itmp <= 0) return(CONVERT_ERROR_RC);
	press->array[press->size].set_id = itmp;

	itmp = nst_tok2i(num_toks, 2, toks);   /* EID */
	if(itmp <= 0) return(CONVERT_ERROR_RC);
	press->array[press->size].element = itmp;

	g1 = nst_tok2i(num_toks, 7, toks);   /* G1 */
	if(g1 <= 0) return(CONVERT_ERROR_RC);

	g34 = nst_tok2i(num_toks, 8, toks);   /* G3 or G4 */
	if(g34 <= 0) return(IMPLEMENT_ERROR_RC);

	index = search_element_label(element, press->array[press->size].element);
	if(index < 0) return(SEEK_ERROR_RC);
	switch(element.array[index].type){
	case ELEM_TETRA1:
	case ELEM_TETRA2:
		tri_flag = 1;
		if(      element.array[index].node[2] == g34){
			press->array[press->size].face_id = 0;
		}else if(element.array[index].node[0] == g34){
			press->array[press->size].face_id = 1;
		}else if(element.array[index].node[3] == g34){
			press->array[press->size].face_id = 2;
		}else if(element.array[index].node[1] == g34){
			press->array[press->size].face_id = 3;
		}else{
			return(UNKNOWN_ERROR_RC);
		}
		break;
	case ELEM_PENTA1:
	case ELEM_PENTA2:
	case ELEM_HEXA1:
	case ELEM_HEXA2:
		RC_TRY( ident_face(element.array[index], g1, g34,
		                   &(press->array[press->size].face_id), &tri_flag) );
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	if(nst_tok2d(num_toks, 3, toks) <= CR_FNSDEFAULT){
		return(CONVERT_ERROR_RC);
	}

	press->array[press->size].type = BC_PRESS;
	if( (nst_tok2d(num_toks, 4, toks) > CR_FNSDEFAULT)
	  &&(nst_tok2d(num_toks, 5, toks) > CR_FNSDEFAULT)
	  &&(tri_flag||(nst_tok2d(num_toks, 6, toks) > CR_FNSDEFAULT)) ){
		press->array[press->size].type = BC_PRESS_D1;
	}
	if(press->array[press->size].type == BC_PRESS){
		press->array[press->size].val.s
		     = nst_tok2d(num_toks, 3, toks);     /* P1 */
	}else{
		press->array[press->size].val.ns[0]
		     = nst_tok2d(num_toks, 3, toks);     /* P1 */
		press->array[press->size].val.ns[1]
		     = nst_tok2d(num_toks, 4, toks);     /* P2 */
		press->array[press->size].val.ns[2]
		     = nst_tok2d(num_toks, 5, toks);     /* P3 */
		if(tri_flag == 0){
			press->array[press->size].val.ns[3]
			     = nst_tok2d(num_toks, 6, toks); /* P4 */
		}
	}

	if( (nst_tok2d(num_toks, 10, toks) > CR_FNSDEFAULT)
	  &&(nst_tok2d(num_toks, 11, toks) > CR_FNSDEFAULT)
	  &&(nst_tok2d(num_toks, 12, toks) > CR_FNSDEFAULT) ){
		if(press->array[press->size].type == BC_PRESS){
			press->array[press->size].type = BC_TRACTION;
		}else{
			press->array[press->size].type = BC_TRACTION_D1;
		}
		press->array[press->size].dir.x
		     = nst_tok2d(num_toks, 10, toks);  /* N1 */
		press->array[press->size].dir.y
		     = nst_tok2d(num_toks, 11, toks);  /* N2 */
		press->array[press->size].dir.z
		     = nst_tok2d(num_toks, 12, toks);  /* N3 */
		press->array[press->size].i_info[0]
		     = nst_tok2i(num_toks, 9, toks);   /* CID */
	}

	(press->size)++;

	return(NORMAL_RC);
}


static RC
ident_face (ELEMENT elem, int g1, int g34, int *face_id, int *tri_flag)
{
	int ii1;
	FACE_TABLE f_table;

	RC_TRY( make_face_table(elem.type, &f_table) );
	for(ii1=0; ii1<f_table.face_num; ii1++){
		switch(f_table.face_type[ii1]){
		case ELEM_QUAD1:
		case ELEM_QUAD2:
			if( ( (f_table.table[ii1][0] == g1)
			    &&(f_table.table[ii1][2] == g34) )
			  ||( (f_table.table[ii1][2] == g1)
			    &&(f_table.table[ii1][0] == g34) )
			  ||( (f_table.table[ii1][1] == g1)
			    &&(f_table.table[ii1][3] == g34) )
			  ||( (f_table.table[ii1][3] == g1)
			    &&(f_table.table[ii1][1] == g34) ) ){
				*face_id = ii1;
				*tri_flag = 0;
				return(NORMAL_RC);
			}
			break;
		case ELEM_TETRA1:
		case ELEM_TETRA2:
			if( (g34 <= 0) &&( (f_table.table[ii1][0] == g1)
			                 ||(f_table.table[ii1][1] == g1)
			                 ||(f_table.table[ii1][2] == g1) ) ){
				*face_id = ii1;
				*tri_flag = 1;
				return(NORMAL_RC);
			}
			if( ( (f_table.table[ii1][0] == g1)
			    &&(f_table.table[ii1][1] == g34) )
			  ||( (f_table.table[ii1][1] == g1)
			    &&(f_table.table[ii1][0] == g34) )
			  ||( (f_table.table[ii1][1] == g1)
			    &&(f_table.table[ii1][2] == g34) )
			  ||( (f_table.table[ii1][2] == g1)
			    &&(f_table.table[ii1][1] == g34) )
			  ||( (f_table.table[ii1][2] == g1)
			    &&(f_table.table[ii1][0] == g34) )
			  ||( (f_table.table[ii1][0] == g1)
			    &&(f_table.table[ii1][2] == g34) ) ){
				*face_id = ii1;
				*tri_flag = 1;
				return(NORMAL_RC);
			}
			break;
		default:
			return(IMPLEMENT_ERROR_RC);
		}
	}

	return(UNKNOWN_ERROR_RC);
}


RC
nst_output_pressure (FILE *fp, ELEM_BC_ARRAY press, ELEMENT_ARRAY element)
{
	int ii1;

	for(ii1=0; ii1<press.size; ii1++){
		if(press.array[ii1].element < 0) continue;

		if( (press.array[ii1].type == BC_PRESS)
		  ||(press.array[ii1].type == BC_PRESS_D1)
		  ||(press.array[ii1].type == BC_TRACTION)
		  ||(press.array[ii1].type == BC_TRACTION_D1) ){
			RC_TRY( put_pload4(fp, press.array[ii1], element, press.source) );
		}else{
			return(IMPLEMENT_ERROR_RC);
		}
	}

	return(NORMAL_RC);
}


static RC
put_pload4 (FILE *fp, ELEM_BC press, ELEMENT_ARRAY element, DATA_TYPE source)
{
	IDC tmp[9];
	FACE_TABLE f_table;
	int index = search_element_label(element, press.element);

	if( (press.type != BC_PRESS)
	  &&(press.type != BC_PRESS_D1)
	  &&(press.type != BC_TRACTION)
	  &&(press.type != BC_TRACTION_D1) ){
		return(IMPLEMENT_ERROR_RC);
	}

	RC_TRY( make_face_table(element.array[index].type, &f_table) );
	if(f_table.face_num <= press.face_id) return(ARG_ERROR_RC);

	tmp[0].c  = "PLOAD4";
	tmp[1].i  = press.set_id;    /* SID */
	tmp[2].i  = press.element;   /* EID */
	if( (press.type == BC_PRESS)||(press.type == BC_TRACTION) ){
		tmp[3].d  = press.val.s;     /* P1 */
		tmp[4].d  = FNSDEFAULT;      /* P2 */
		tmp[5].d  = FNSDEFAULT;      /* P3 */
		tmp[6].d  = FNSDEFAULT;      /* P4 */
	}else{
		tmp[3].d  = press.val.ns[0]; /* P1 */
		tmp[4].d  = press.val.ns[1]; /* P2 */
		tmp[5].d  = press.val.ns[2]; /* P3 */
		if( (f_table.face_type[press.face_id] == ELEM_QUAD1)
		  ||(f_table.face_type[press.face_id] == ELEM_QUAD2) ){
			tmp[6].d  = press.val.ns[3]; /* P4 */
		}else{
			tmp[6].d  = FNSDEFAULT;      /* P4 */
		}
	}

	RC_NEG_CHK( index );
	switch(element.array[index].type){
	case ELEM_TETRA1:
	case ELEM_TETRA2:
		if(      press.face_id == 0){
			tmp[7].i = element.array[index].node[0];    /* G1 */
			tmp[8].i = element.array[index].node[2];    /* G3 or G4 */
		}else if(press.face_id == 1){
			tmp[7].i = element.array[index].node[1];    /* G1 */
			tmp[8].i = element.array[index].node[0];    /* G3 or G4 */
		}else if(press.face_id == 2){
			tmp[7].i = element.array[index].node[0];    /* G1 */
			tmp[8].i = element.array[index].node[3];    /* G3 or G4 */
		}else if(press.face_id == 3){
			tmp[7].i = element.array[index].node[0];    /* G1 */
			tmp[8].i = element.array[index].node[1];    /* G3 or G4 */
		}else{
			return(UNKNOWN_ERROR_RC);
		}
		break;
	case ELEM_PENTA1:
	case ELEM_PENTA2:
		if(      press.face_id == 0){
			tmp[7].i = element.array[index].node[0];    /* G1 */
			tmp[8].i = NSDEFAULT;                       /* G3 or G4 */
		}else if(press.face_id == 1){
			tmp[7].i = element.array[index].node[3];    /* G1 */
			tmp[8].i = NSDEFAULT;                       /* G3 or G4 */
		}else if(press.face_id == 2){
			tmp[7].i = element.array[index].node[1];    /* G1 */
			tmp[8].i = element.array[index].node[5];    /* G3 or G4 */
		}else if(press.face_id == 3){
			tmp[7].i = element.array[index].node[0];    /* G1 */
			tmp[8].i = element.array[index].node[5];    /* G3 or G4 */
		}else if(press.face_id == 4){
			tmp[7].i = element.array[index].node[0];    /* G1 */
			tmp[8].i = element.array[index].node[4];    /* G3 or G4 */
		}else{
			return(UNKNOWN_ERROR_RC);
		}
		break;
	case ELEM_HEXA1:
	case ELEM_HEXA2:
		if(      press.face_id == 0){
			tmp[7].i = element.array[index].node[0];    /* G1 */
			tmp[8].i = element.array[index].node[2];    /* G3 or G4 */
		}else if(press.face_id == 1){
			tmp[7].i = element.array[index].node[0];    /* G1 */
			tmp[8].i = element.array[index].node[5];    /* G3 or G4 */
		}else if(press.face_id == 2){
			tmp[7].i = element.array[index].node[1];    /* G1 */
			tmp[8].i = element.array[index].node[6];    /* G3 or G4 */
		}else if(press.face_id == 3){
			tmp[7].i = element.array[index].node[7];    /* G1 */
			tmp[8].i = element.array[index].node[5];    /* G3 or G4 */
		}else if(press.face_id == 4){
			tmp[7].i = element.array[index].node[3];    /* G1 */
			tmp[8].i = element.array[index].node[4];    /* G3 or G4 */
		}else if(press.face_id == 5){
			tmp[7].i = element.array[index].node[2];    /* G1 */
			tmp[8].i = element.array[index].node[7];    /* G3 or G4 */
		}else{
			return(UNKNOWN_ERROR_RC);
		}
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}
			
	nst_print(fp, "ciiddddiir", tmp);

	if( (press.type == BC_TRACTION)||(press.type == BC_TRACTION_D1) ){
		if(source == FEM_NASTRAN){
			tmp[0].i = press.i_info[0];    /* CID */
		}else{
			tmp[0].i = NSDEFAULT;          /* CID */
		}
		tmp[1].d = press.dir.x;            /* N1 */
		tmp[2].d = press.dir.y;            /* N2 */
		tmp[3].d = press.dir.z;            /* N3 */
		nst_print(fp, "-idddr", tmp);
	}

	return(NORMAL_RC);
}


RC
nst_input_node_matrix (FILE *fp, char *name, NODE_MATRIX_ARRAY *nmat)
{
	const int num_keys = 1;
	NST_BULK_KEY keys[] = { {"DMIG", 4} };
	int hit_key;
	char *toks[TOKS_SIZE_LARGE];
	int num_toks;
	NST_FILE nst_file;
	int name_len;
	int gj, cj;
	int g1, c1;
	double a1, b1;
	int index, ii1, ii2;
	int sub_row, sub_col;
	int shift_row, shift_col;
	int col_index;

	if(name == NULL) return(ARG_ERROR_RC);
	name_len = strlen(name);
	if(name_len > 8) name_len = 8;

	RC_TRY( nst_allocate_tok(toks, TOKS_SIZE_LARGE) );
	RC_TRY( init_nst_file(fp, &nst_file) );
	RC_TRY( nst_seek_begin_bulk(&nst_file) );
	RC_TRY( allocate_node_matrix_array(INIT_SIZE, nmat) );
	nmat->source = FEM_NASTRAN;

	while(1){
		RC rc = nst_blk_get_tok(&nst_file, num_keys, keys, &hit_key,
		                        TOKS_SIZE_LARGE, toks, &num_toks);
		if(rc != NORMAL_RC){
			if(rc == END_RC) break;
			return(rc);
		}

		switch(hit_key){
		case 0:
			if(nst_tok_strncmp(num_toks, 1, toks, name, name_len)) continue;
			if(nst_tok2i(num_toks, 2, toks) <= 0) continue;

			gj = nst_tok2i(num_toks, 2, toks);   /* GJ */
			cj = nst_tok2i(num_toks, 3, toks);   /* CJ */
			if( (cj < 0)||(cj > 6) ) return(READ_ERROR_RC);
			shift_row = 0;
			sub_row = cj - 1;
			if(cj >= 4){
				shift_row = 1;
				sub_row -= 3;
			}

			RC_TRY( realloc_node_matrix_array(nmat) );
			index = nmat->size;
			for(ii1=nmat->size-1; ii1>=0; ii1--){
				if( (nmat->array[ii1].node == gj)
				  &&(nmat->array[ii1].shift == shift_row) ){
					index = ii1;
					break;
				}
			}
			if(index == nmat->size){
				nmat->array[index].node = gj;
				nmat->array[index].shift = shift_row;
				nmat->size++;
			}

			for(ii1=5; ii1<num_toks; ii1+=4){
				g1 = nst_tok2i(num_toks, ii1, toks);     /* G? */
				c1 = nst_tok2i(num_toks, ii1+1, toks);   /* C? */
				a1 = nst_tok2d(num_toks, ii1+2, toks);   /* A? */
				b1 = nst_tok2d(num_toks, ii1+3, toks);   /* B? */
				if( (g1 <= 0)||(c1 <= 0) ) break;
				if(a1 < CR_FNSDEFAULT) return(READ_ERROR_RC);
				if(c1 > 6) return(READ_ERROR_RC);
				shift_col = 0;
				sub_col = c1 - 1;
				if(c1 >= 4){
					shift_col = 1;
					sub_col -= 3;
				}

				RC_TRY( realloc_node_matrix(&(nmat->array[index])) );
				col_index = nmat->array[index].size;
				for(ii2=nmat->array[index].size-1; ii2>=0; ii2--){
					if( (nmat->array[index].col_node[ii2] == g1)
					  &&(nmat->array[index].col_shift[ii2] == shift_col) ){
						col_index = ii2;
						break;
					}
				}
				if(col_index == nmat->array[index].size){
					nmat->array[index].col_node[col_index] = g1;
					nmat->array[index].col_shift[col_index] = shift_col;
					nmat->array[index].size++;
				}
				nmat->array[index].matrix[col_index][sub_row][sub_col] += a1;
				if( (gj == g1)&&(shift_col == shift_row)
				                  &&(sub_col != sub_row) ){
					nmat->array[index].matrix[col_index][sub_col][sub_row]+=a1;
				}
			}
			if(nmat->array[index].alloc_size > nmat->array[index].size){
				nmat->array[index].col_node
				    = mm_realloc(nmat->array[index].col_node,
				         nmat->array[index].size * sizeof(int));
				if(nmat->array[index].col_node == NULL)
					return(ALLOC_ERROR_RC);

				nmat->array[index].col_shift
				    = mm_realloc(nmat->array[index].col_shift,
				                 nmat->array[index].size*sizeof(int));
				if(nmat->array[index].col_shift == NULL) return(ALLOC_ERROR_RC);
				nmat->array[index].matrix
				    = mm_realloc(nmat->array[index].matrix,
				                 nmat->array[index].size*sizeof(double[3][3]));
				if(nmat->array[index].matrix == NULL) return(ALLOC_ERROR_RC);

				nmat->array[index].alloc_size = nmat->array[index].size;
			}
			break;
		default:
			return(UNKNOWN_ERROR_RC);
		}
	}

	RC_TRY( clean_node_matrix_array(nmat) );
	RC_TRY( term_nst_file(&nst_file) );
	RC_TRY( nst_free_tok(toks, TOKS_SIZE_LARGE) );

	return(NORMAL_RC);
}


RC
nst_output_node_matrix (FILE *fp, char *name, NODE_MATRIX_ARRAY nmat)
{
	int ii1, ii2, ii3;
	IDC tmp[6];

	if(nmat.size <= 0) return(NORMAL_RC);

	tmp[0].c = "DMIG";
	tmp[1].c = name;
	tmp[2].i = 0;
	tmp[3].i = 6;
	tmp[4].i = 2;
	tmp[5].i = 0;
	RC_TRY( nst_print(fp, "cciiiir", tmp) );

	for(ii1=0; ii1<nmat.size; ii1++){
		if(nmat.array[ii1].node < 0) continue;
		for(ii2=0; ii2<3; ii2++){
			tmp[0].c = "DMIG*";
			tmp[1].c = name;
			tmp[2].i = nmat.array[ii1].node;
			tmp[3].i = 1 + ii2 + 3*nmat.array[ii1].shift;
			RC_TRY( nst_print(fp, "cc-II=r", tmp) );
			for(ii3=0; ii3<nmat.array[ii1].size; ii3++){
				int dia_flag = 0;
				if((nmat.array[ii1].node == nmat.array[ii1].col_node[ii3])
				 &&(nmat.array[ii1].shift == nmat.array[ii1].col_shift[ii3])){
					dia_flag = 1;
				}

				tmp[0].c = "*";
				tmp[1].i = nmat.array[ii1].col_node[ii3];
				tmp[2].i = 1 + 3*nmat.array[ii1].col_shift[ii3];
				tmp[3].d = nmat.array[ii1].matrix[ii3][ii2][0];
				RC_TRY( nst_print(fp, "cIIG=r", tmp) );
				if(dia_flag && (ii2 == 0)) continue;

				tmp[0].c = "*";
				tmp[1].i = nmat.array[ii1].col_node[ii3];
				tmp[2].i = 2 + 3*nmat.array[ii1].col_shift[ii3];
				tmp[3].d = nmat.array[ii1].matrix[ii3][ii2][1];
				RC_TRY( nst_print(fp, "cIIG=r", tmp) );
				if(dia_flag && (ii2 == 1)) continue;

				tmp[0].c = "*";
				tmp[1].i = nmat.array[ii1].col_node[ii3];
				tmp[2].i = 3 + 3*nmat.array[ii1].col_shift[ii3];
				tmp[3].d = nmat.array[ii1].matrix[ii3][ii2][2];
				RC_TRY( nst_print(fp, "cIIG=r", tmp) );
				if(dia_flag && (ii2 == 2)) continue;
			}
		}
	}

	return(NORMAL_RC);
}


RC
nst_input_heatflux (FILE *fp, BC_ARRAY *flux)
{
	const int num_keys = 1;
	NST_BULK_KEY keys[] = { {"QHBDY",5} };
	int hit_key;
	char *toks[TOKS_SIZE];
	int num_toks;
	NST_FILE nst_file;
	
	RC_TRY( nst_allocate_tok(toks, TOKS_SIZE) );
	RC_TRY( init_nst_file(fp, &nst_file) );
	RC_TRY( nst_seek_begin_bulk(&nst_file) );
	RC_TRY( allocate_bc_array(INIT_SIZE, flux) );

	flux->source = FEM_NASTRAN;
	
	while(1){
		RC rc = nst_blk_get_tok(&nst_file, num_keys, keys, &hit_key,
		                        TOKS_SIZE, toks, &num_toks);
		if(rc != NORMAL_RC){
			if(rc == END_RC) break;
			return(rc);
		}
		
		switch(hit_key){
		case 0:
			RC_TRY( set_qload(flux, toks, num_toks) );
			break;
		default:
			return(UNKNOWN_ERROR_RC);
		}
	}
	
	RC_TRY( clean_bc_array(flux) );
	RC_TRY( term_nst_file(&nst_file) );
	RC_TRY( nst_free_tok(toks, TOKS_SIZE) );
	
	return(NORMAL_RC);
}


static RC
set_qload (BC_ARRAY *flux, char **toks, int num_toks)
{
	int sid;       /* data type */
	int itmp;      /* node label */
	double dtmp;   /* heat flux */
	
	sid = nst_tok2i(num_toks, 1, toks);
	if(sid < 0) return(CONVERT_ERROR_RC);
	
	itmp = nst_tok2i(num_toks, 5, toks);
	dtmp = nst_tok2d(num_toks, 3, toks);
	
	flux->array[flux->size].node = itmp;
	flux->array[flux->size].set_id = sid;
	flux->array[flux->size].s[0] = dtmp;
	(flux->size)++;
	return(NORMAL_RC);
}


RC
nst_input_heattransfer (FILE *fp, ELEMENT_ARRAY surf, ELEMENT_ARRAY *transfer)
{
	const int num_keys = 4;
	NST_BULK_KEY keys[] = { {"MAT4", 4}, {"SPC", 3},
	                        {"CONV", 4}, {"CHBDYP", 6}};
	int hit_key;
	char *toks[TOKS_SIZE];
	int num_toks;
	int sid = -1;
	NST_FILE nst_file;
	double coefficient = 0.0;
	double airtmp = 0.0;
	int bc_label, elem_label = -1;
	int node_label[4];
	int ii1, ii2;
	
	RC_TRY( nst_allocate_tok(toks, TOKS_SIZE) );
	RC_TRY( init_nst_file(fp, &nst_file) );
	RC_TRY( nst_seek_begin_bulk(&nst_file) );
	RC_TRY( allocate_element_array(INIT_SIZE, transfer) );
	
	transfer->source = FEM_NASTRAN;
	
	while(1){
		RC rc = nst_blk_get_tok(&nst_file, num_keys, keys, &hit_key,
		                        TOKS_SIZE, toks, &num_toks);
		if(rc != NORMAL_RC){
			if(rc == END_RC) break;
			return(rc);
		}
		
		switch(hit_key){
		case 0:
			sid = nst_tok2i(num_toks, 1, toks);
			if(sid < 0) return(CONVERT_ERROR_RC);
			coefficient = nst_tok2d(num_toks, 5, toks);
			break;
		case 1:
			airtmp = nst_tok2d(num_toks, 4, toks);
			break;
		case 2:
			transfer->array[transfer->size].i_info[0] = sid;
			transfer->array[transfer->size].d_info[0] = coefficient;
			transfer->array[transfer->size].d_info[1] = airtmp;
			transfer->array[transfer->size].i_info[1] = nst_tok2i(num_toks, 1,
			                                                      toks);
			(transfer->size)++;
			break;
		case 3:
			bc_label = nst_tok2i(num_toks, 1, toks);
			node_label[0] = nst_tok2i(num_toks, 6, toks);
			node_label[1] = nst_tok2i(num_toks, 7, toks);
			node_label[2] = nst_tok2i(num_toks, 8, toks);
			
			for(ii1=0; ii1<surf.size; ii1++){
				for(ii2=0; ii2<surf.array[ii1].node_num; ii2++){
					if(surf.array[ii1].node[ii2] == node_label[0]){
						for(ii2=0; ii2<surf.array[ii1].node_num; ii2++){
							if(surf.array[ii1].node[ii2] == node_label[1]){
								elem_label = surf.array[ii1].label;
							}
						}
					}
				}
			}
			/* ごちゃごちゃしすぎ、もっと簡潔にする */
			for(ii1=0; ii1<(transfer->size); ii1++){
				if( transfer->array[ii1].i_info[1] == bc_label){
					transfer->array[ii1].label = elem_label;
					transfer->array[ii1].node[0] = node_label[0];
					transfer->array[ii1].node[1] = node_label[1];
					transfer->array[ii1].node[2] = node_label[2];
				}
			}
			/* もしLINEだったら1とかにして場合わけするか？ */
			break;
		default:
			return(UNKNOWN_ERROR_RC);
		}
	}

	RC_TRY( clean_element_array(transfer) );
	RC_TRY( term_nst_file(&nst_file) );
	RC_TRY( nst_free_tok(toks, TOKS_SIZE) );
	
	return(NORMAL_RC);
}


