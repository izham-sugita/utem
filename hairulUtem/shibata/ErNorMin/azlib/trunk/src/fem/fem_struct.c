/*********************************************************************
 * fem_struct.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Tomoyuki TSUBATA>
 *
 * Thanks to following contributors.
 *  <Yasuo ASAGA> <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: fem_struct.c 1089 2016-11-30 09:36:47Z hayashi $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rc.h"
#include "math_utl.h"
#include "fem_struct.h"
#include "memory_manager.h"

#define SHEAR_FACT (5.0 / 6.0)     /* 修正せん断係数 */

/* local functions */
static void clean_D_matrix(MATERIAL_PROP *prop);
static void clean_C(MATERIAL_PROP *prop);
static RC put_iso_D_mat(MATERIAL_PROP *prop);
static RC put_ortho_D_mat(MATERIAL_PROP *prop);


/* id1[num1], scale1[num1] と id2[num2], scale2[num2] 比較 */
int
same_bc_set (int num1, const int id1[], double scale1[],
             int num2, const int id2[], double scale2[])
{
	int ii1;

	if(num1 != num2) return(0);

	for(ii1=0; ii1<num1; ii1++){
		if(id1[ii1] != id2[ii1]) return(0);
		if(!nearly_eq(scale1[ii1], scale2[ii1])) return(0);
	}

	return(1);
}


/* 要素 s_elem の全ての節点が、v_elem にも存在する時    */
/* s_elem 上の局所座標 s_coords[coord_size][4] を       */
/* v_elem 上の局所座標 v_coords[coord_size][4] に変換   */
RC
trans_coord_s2v(ELEMENT s_elem, ELEMENT v_elem, int coord_size,
                double **s_coords, double **v_coords)
{
	int ii1, ii2, ii3;
	double **v_node_coord;
	double s_N[MAX_NODE];
	int v_index;

	RC_TRY( allocate2D(MAX_NODE, 4, &v_node_coord) );
	RC_TRY( set_local_node_points(v_elem.type, v_node_coord) );

	for(ii1=0; ii1<coord_size; ii1++){
		for(ii2=0; ii2<4; ii2++){
			v_coords[ii1][ii2] = 0.0;
		}
		RC_TRY( set_N_vector(s_elem.type, s_coords[ii1], s_N) );

		/* v_elem の各節点での局所座標値を s_elem 内で内挿 */
		for(ii2=0; ii2<s_elem.node_num; ii2++){
			v_index = -1;
			for(ii3=0; ii3<v_elem.node_num; ii3++){
				if(v_elem.node[ii3] == s_elem.node[ii2]){
					v_index = ii3;
					break;
				}
			}
			RC_NEG_CHK( v_index );

			for(ii3=0; ii3<4; ii3++){
				v_coords[ii1][ii3] += s_N[ii2] * v_node_coord[v_index][ii3];
			}
		}
	}

	RC_TRY( free2D(MAX_NODE, 4, &v_node_coord) );

	return(NORMAL_RC);
}


RC
allocate_work_set (WORK_SET *ws)
{
	int ii1;

	for(ii1=0; ii1<WS_SIZE; ii1++){
		RC_TRY( allocate1D(MAX_NODE, &(ws->vec_N[ii1])) );
		RC_TRY( allocate1D(MAX_GAUSS_POINT, &(ws->vec_G[ii1])) );

		RC_TRY( allocate2D(3, 3, &(ws->mat_3_3[ii1])) );
		RC_TRY( allocate2D(6, 6, &(ws->mat_6_6[ii1])) );
		RC_TRY( allocate2D(9, 9, &(ws->mat_9_9[ii1])) );

		RC_TRY( allocate2D(3, MAX_NODE, &(ws->mat_3_N[ii1])) );
		RC_TRY( allocate2D(3, 3*MAX_NODE, &(ws->mat_3_3N[ii1])) );
		RC_TRY( allocate2D(6, 3*MAX_NODE, &(ws->mat_6_3N[ii1])) );
		RC_TRY( allocate2D(9, 3*MAX_NODE, &(ws->mat_9_3N[ii1])) );
		RC_TRY( allocate2D(3*MAX_NODE, 3*MAX_NODE,
		                                      &(ws->mat_3N_3N[ii1])) );

		RC_TRY( allocate2D(MAX_NODE, 3, &(ws->mat_N_3[ii1])) );
		RC_TRY( allocate2D(MAX_NODE, 4, &(ws->mat_N_4[ii1])) );
		RC_TRY( allocate2D(MAX_NODE, 6, &(ws->mat_N_6[ii1])) );
		RC_TRY( allocate2D(MAX_GAUSS_POINT, 4, &(ws->mat_G_4[ii1])) );
	}

	return(NORMAL_RC);
}


RC
free_work_set (WORK_SET *ws)
{
	int ii1;

	for(ii1=0; ii1<WS_SIZE; ii1++){
		RC_TRY( mm_free(ws->vec_N[ii1]) );
		RC_TRY( mm_free(ws->vec_G[ii1]) );
		ws->vec_N[ii1] = NULL;
		ws->vec_G[ii1] = NULL;

		RC_TRY( free2D(3, 3, &(ws->mat_3_3[ii1])) );
		RC_TRY( free2D(6, 6, &(ws->mat_6_6[ii1])) );
		RC_TRY( free2D(9, 9, &(ws->mat_9_9[ii1])) );
		ws->mat_3_3[ii1] = NULL;
		ws->mat_6_6[ii1] = NULL;
		ws->mat_9_9[ii1] = NULL;

		RC_TRY( free2D(3, MAX_NODE,   &(ws->mat_3_N[ii1])) );
		RC_TRY( free2D(3, 3*MAX_NODE, &(ws->mat_3_3N[ii1])) );
		RC_TRY( free2D(6, 3*MAX_NODE, &(ws->mat_6_3N[ii1])) );
		RC_TRY( free2D(9, 3*MAX_NODE, &(ws->mat_9_3N[ii1])) );
		RC_TRY( free2D(3*MAX_NODE, 3*MAX_NODE,
		                                  &(ws->mat_3N_3N[ii1])) );
		ws->mat_3_N[ii1] = NULL;
		ws->mat_3_3N[ii1] = NULL;
		ws->mat_6_3N[ii1] = NULL;
		ws->mat_9_3N[ii1] = NULL;
		ws->mat_3N_3N[ii1] = NULL;

		RC_TRY( free2D(MAX_NODE, 3, &(ws->mat_N_3[ii1])) );
		RC_TRY( free2D(MAX_NODE, 4, &(ws->mat_N_4[ii1])) );
		RC_TRY( free2D(MAX_NODE, 6, &(ws->mat_N_6[ii1])) );
		RC_TRY( free2D(MAX_GAUSS_POINT, 4, &(ws->mat_G_4[ii1])) );
		ws->mat_N_3[ii1] = NULL;
		ws->mat_N_4[ii1] = NULL;
		ws->mat_N_6[ii1] = NULL;
		ws->mat_G_4[ii1] = NULL;
	}

	return(NORMAL_RC);
}


/* initialization of NODE */
void
init_node (NODE *node)
{
	int ii1;
	
	node->label = -1;
	init_vect3d(&node->p);
	node->renum = -1;
	node->dummy_label = -1;
	for(ii1=0;ii1<NODE_INFO;ii1++){
		node->i_info[ii1] = 0;
		node->d_info[ii1] = 0.0;
	}
}


/* initialization of ELEMENT */
void
init_element (ELEMENT *elem)
{
	int ii;

	elem->label = -1;
	elem->type = ELEM_TRI1;
	elem->material = -1;
	elem->physical = -1;
	elem->node_num = 0;
	for(ii=0;ii<MAX_NODE;ii++){
		elem->node[ii] = -1;
	}
	elem->volume = 0.0;
	for(ii=0; ii<ELEMENT_INFO; ii++){
		elem->i_info[ii] = 0;
		elem->d_info[ii] = 0.0;
	}
}


/* initialization of RIGID_ELEMENT */
void
init_rigid_element (RIGID_ELEMENT *elem)
{
	int ii;

	elem->label = -1;
	elem->set_id = -1;
	elem->type = ELEM_RBAR;
	elem->physical = -1;
	elem->node_num = 0;
	elem->node = NULL;
	elem->weight_num = 0;
	elem->weight = NULL;
	elem->dof_num = 0;
	elem->dof = NULL;
	for(ii=0; ii<RIGID_ELEMENT_INFO; ii++){
		elem->i_info[ii] = 0;
		elem->d_info[ii] = 0.0;
	}
}


RC
free_rigid_element (RIGID_ELEMENT *elem)
{
	if( (elem->label >= 0)&&(elem->node_num > 0) ){
		if(elem->node != NULL)   RC_TRY( mm_free(elem->node) );
		if(elem->weight != NULL) RC_TRY( mm_free(elem->weight) );
		if(elem->dof != NULL)    RC_TRY( mm_free(elem->dof) );
	}
	init_rigid_element(elem);

	return(NORMAL_RC);
}


RC
copy_rigid_element (RIGID_ELEMENT src, RIGID_ELEMENT *dest)
{
	int ii1;

	*dest = src;
	if(src.node_num > 0){
		dest->node = mm_alloc(src.node_num * sizeof(int));
		if(dest->node == NULL) return(ALLOC_ERROR_RC); 
		for(ii1=0; ii1<src.node_num; ii1++){
			dest->node[ii1] = src.node[ii1];
		}
	}else{
		dest->node = NULL;
	}

	if(src.weight_num > 0){
		dest->weight = mm_alloc(src.weight_num * sizeof(double));
		if(dest->weight == NULL) return(ALLOC_ERROR_RC); 
		for(ii1=0; ii1<src.weight_num; ii1++){
			dest->weight[ii1] = src.weight[ii1];
		}
	}else{
		dest->weight = NULL;
	}

	if(src.dof_num > 0){
		dest->dof = mm_alloc(src.dof_num * sizeof(int));
		if(dest->dof == NULL) return(ALLOC_ERROR_RC); 
		for(ii1=0; ii1<src.dof_num; ii1++){
			dest->dof[ii1] = src.dof[ii1];
		}
	}else{
		dest->dof = NULL;
	}

	return(NORMAL_RC);
}


/* initialization of MATERIAL_PROP */
void
init_material_prop (MATERIAL_PROP *prop)
{
	int ii1;

	prop->label = -1;
	prop->mat_type = MAT_ISOTROPIC;
	prop->ana_type = ANA_3D;

	prop->E = 0.0;
	prop->G = 0.0;
	prop->nu = 0.0;
	prop->alpha = 0.0;
	prop->K = 0.0;

	init_vect3d(&prop->E_vect);
	init_vect3dr(&prop->G_vect);
	init_vect3dr(&prop->nu_vect);
	init_vect3dr(&prop->alpha_vect);
	init_vect3dr(&prop->K_vect);

	clean_D_matrix(prop);

	clean_C(prop);
	prop->rho = 0.0;
	prop->tref = 0.0;
	for(ii1=0;ii1<MATERIAL_INFO;ii1++){
		prop->i_info[ii1] = 0;
		prop->d_info[ii1] = 0.0;
	}
}


/* initialization of D_matrix */
static void
clean_D_matrix (MATERIAL_PROP *prop)
{
	int ii1,ii2;

	for(ii1=0;ii1<6;ii1++){
		for(ii2=0;ii2<6;ii2++){
			prop->D_matrix[ii1][ii2] = 0.0;
		}
	}
}


/* initialization of C of material */
static void
clean_C (MATERIAL_PROP *prop)
{
	int ii1,ii2,ii3,ii4;

	for(ii1=0;ii1<3;ii1++){
		for(ii2=0;ii2<3;ii2++){
			for(ii3=0;ii3<3;ii3++){
				for(ii4=0;ii4<3;ii4++){
					prop->C[ii1][ii2][ii3][ii4] = 0.0;
				}
			}
		}
	}
}


/* initialization of PHYSICAL_PROP */
void
init_physical_prop (PHYSICAL_PROP *prop)
{
	int ii1;

	prop->label = -1;
	prop->section = SECTION_ROD;
	for(ii1=0;ii1<PHYSICAL_INFO;ii1++){
		prop->i_info[ii1] = 0;
		prop->d_info[ii1] = 0.0;
	}
}


/* initialization of BC */
void
init_bc (BC *bc)
{
	int ii1;

	bc->node = -1;
	bc->set_id = -1;
	bc->v_type.x = BC_FREE;
	bc->v_type.y = BC_FREE;
	bc->v_type.z = BC_FREE;
	bc->v_type.yz = BC_FREE;
	bc->v_type.zx = BC_FREE;
	bc->v_type.xy = BC_FREE;
	init_vect3dr(&bc->v);
	for(ii1=0;ii1<BC_SCALAR;ii1++){
		bc->s_type[ii1] = BC_FREE;
		bc->s[ii1] = 0.0;
	}
	for(ii1=0;ii1<BC_INFO;ii1++){
		bc->i_info[ii1] = 0;
		bc->d_info[ii1] = 0.0;
	}
}


/* initialization of ELEM_BC */
void
init_elem_bc (ELEM_BC *bc)
{
	int ii1;

	bc->element = -1;
	bc->set_id = -1;
	bc->type = BC_FREE;
	bc->face_id = -1;
	for(ii1=0; ii1<MAX_NODE; ii1++){
		bc->val.ns[ii1] = 0.0;
	}
	bc->dir.x = bc->dir.y = bc->dir.z = 0.0;
	for(ii1=0;ii1<ELEM_BC_INFO;ii1++){
		bc->i_info[ii1] = 0;
		bc->d_info[ii1] = 0.0;
	}
}


/* initialization of DISP_VELO */
void
init_disp_velo (DISP_VELO *disp)
{
	int ii1;

	disp->node = -1;
	init_vect3dr(&disp->v);
	disp->s = 0.0;
	for(ii1=0;ii1<DISP_VELO_INFO;ii1++){
		disp->i_info[ii1] = 0;
		disp->d_info[ii1] = 0.0;
	}
}


/* initialization of STRESS_STRAIN */
void
init_stress_strain (STRESS_STRAIN *ss)
{
	int ii1;

	ss->element = -1;
	ss->node = -1;
	init_vect3dr(&ss->v);
	for(ii1=0;ii1<STRESS_STRAIN_INFO;ii1++){
		ss->i_info[ii1] = 0;
		ss->d_info[ii1] = 0.0;
	}
}


void
init_sensitivity (SENSITIVITY *sens)
{
	sens->element = -1;
	sens->node = -1;
	sens->s = 0.0;
	init_vect3dr(&sens->v);
}


void
init_default_bc (DEFAULT_BC *def_bc)
{
	int ii1;

	def_bc->set_id = -1;
	def_bc->temp = 0.0;
	init_vect3d(&def_bc->grav);

	for(ii1=0; ii1<DEFAULT_BC_INFO; ii1++){
		def_bc->i_info[ii1] = 0;
		def_bc->d_info[ii1] = 0.0;
	}
}


void
init_bc_set (BC_SET *bc_set)
{
	int ii1;

	bc_set->label = -1;
	bc_set->num_fix = 0;
	bc_set->num_mpc = 0;
	bc_set->num_force = 0;
	bc_set->num_temp = 0;

	for(ii1=0; ii1<MAX_BC_SET; ii1++){
		bc_set->id_fix[ii1] = -1;
		bc_set->id_mpc[ii1] = -1;
		bc_set->id_force[ii1] = -1;
		bc_set->id_temp[ii1] = -1;

		bc_set->scale_fix[ii1] = 0.0;
		bc_set->scale_mpc[ii1] = 0.0;
		bc_set->scale_force[ii1] = 0.0;
		bc_set->scale_temp[ii1] = 0.0;
	}

	for(ii1=0; ii1<BC_SET_INFO; ii1++){
		bc_set->i_info[ii1] = 0;
		bc_set->d_info[ii1] = 0.0;
	}
	init_string_array(&(bc_set->s_info));
}


RC
free_bc_set (BC_SET *bc_set)
{
	RC_TRY( free_string_array(&(bc_set->s_info)) );
	init_bc_set(bc_set);

	return(NORMAL_RC);
}


void
init_sound_pressure (SOUND_PRESSURE *sound)
{
	int ii1;

	sound->node = -1;
	sound->re = 0.0;
	sound->im = 0.0;

	for(ii1=0; ii1<SOUND_PRESSURE_INFO; ii1++){
		sound->i_info[ii1] = 0;
		sound->d_info[ii1] = 0.0;
	}
}


void
init_load_curve (LOAD_CURVE *curve)
{
	int ii1;

	curve->label = -1;
	curve->max_time = 0.0;
	curve->max_value = 0.0;
	curve->total_step = 0;
	for(ii1=0; ii1<LOAD_CURVE_MAX_INC; ii1++){
		curve->inc_value[ii1] = 0.0;
	}
	curve->approx_size = 0;
	for(ii1=0; ii1<LOAD_CURVE_MAX_APPROX; ii1++){
		curve->approx_time[ii1] = 0.0;
		curve->approx_value[ii1] = 0.0;
	}
}


void
init_local_coord (LOCAL_COORD *local_coord)
{
	int ii1;

	local_coord->label = -1;
	init_vect3d(&local_coord->origin);
	init_vect3d(&local_coord->d_cos[0]);
	init_vect3d(&local_coord->d_cos[1]);
	init_vect3d(&local_coord->d_cos[2]);
	for(ii1=0; ii1<LOCAL_COORD_INFO; ii1++){
		local_coord->i_info[ii1] = 0;
		local_coord->d_info[ii1] = 0.0;
	}
}


void
init_elem_matrix (ELEM_MATRIX *elem_matrix)
{
	elem_matrix->element = -1;
	elem_matrix->size = 0;
	elem_matrix->dim = 0;
	elem_matrix->index = NULL;
	elem_matrix->matrix = NULL;
}


RC
free_elem_matrix (ELEM_MATRIX *elem_matrix)
{
	if((elem_matrix->element >= 0)&&(elem_matrix->size > 0)){
		if(elem_matrix->index == NULL) return(ARG_ERROR_RC);

		RC_TRY( mm_free(elem_matrix->index) );
		RC_TRY( free2D(elem_matrix->size,
		               elem_matrix->size, &elem_matrix->matrix) );
	}
	init_elem_matrix(elem_matrix);

	return(NORMAL_RC);
}


void
init_elem_lookup (ELEM_LOOKUP *elem_lookup)
{
	elem_lookup->node = -1;
	elem_lookup->size = 0;
	elem_lookup->alloc_size = 0;
	elem_lookup->element = NULL;
	elem_lookup->r_size = 0;
	elem_lookup->r_alloc_size = 0;
	elem_lookup->r_element = NULL;
}


void
init_node_matrix (NODE_MATRIX *nmat)
{
	nmat->node = -1;
	nmat->shift = 0;
	nmat->size = 0;
	nmat->alloc_size = 0;
	nmat->col_node = NULL;
	nmat->col_shift = NULL;
	nmat->matrix = NULL;
}


RC
free_elem_lookup (ELEM_LOOKUP *elem_lookup)
{
	if((elem_lookup->node >= 0)&&(elem_lookup->alloc_size > 0)){
		if(elem_lookup->element == NULL) return(ARG_ERROR_RC);
		RC_TRY( mm_free(elem_lookup->element) );
	}
	if((elem_lookup->node >= 0)&&(elem_lookup->r_alloc_size > 0)){
		if(elem_lookup->r_element == NULL) return(ARG_ERROR_RC);
		RC_TRY( mm_free(elem_lookup->r_element) );
	}
	init_elem_lookup(elem_lookup);

	return(NORMAL_RC);
}


RC
free_node_matrix (NODE_MATRIX *nmat)
{
	if(nmat->alloc_size > 0){
		if(nmat->col_node == NULL) return(ARG_ERROR_RC);
		if(nmat->col_shift == NULL) return(ARG_ERROR_RC);
		if(nmat->matrix == NULL) return(ARG_ERROR_RC);
		RC_TRY( mm_free(nmat->col_node) );
		RC_TRY( mm_free(nmat->col_shift) );
		RC_TRY( mm_free(nmat->matrix) );
	}
	init_node_matrix(nmat);

	return(NORMAL_RC);
}


void
init_contact_term (CONTACT_TERM *term)
{
	term->element = -1;
	init_vect3d(&term->p);
	term->distance = -1.0*DBL_MAX;
	term->type = CONTACT_POINT;
}


void
init_contact (CONTACT *contact)
{
	contact->node = -1;
	contact->size = 0;
	contact->alloc_size = 0;
	contact->term = NULL;
	contact->sort_flag = ARRAY_UNSORTED;
}


RC
free_contact (CONTACT *contact)
{
	if(contact->alloc_size > 0){
		if(contact->term == NULL) return(ARG_ERROR_RC);
		RC_TRY( mm_free(contact->term) );
		init_contact(contact);
	}
	return(NORMAL_RC);
}


void
init_geom_segment (GEOM_SEGMENT *segment)
{
	segment->index = -1;
	segment->min =  DBL_MAX;
	segment->max = -DBL_MAX;
}


STRESS_ARRAY
strain2stress_array (STRAIN_ARRAY strain)
{
	STRESS_ARRAY ret;

	ret.size = strain.size;
	ret.alloc_size = strain.alloc_size;
	ret.source = strain.source;
	ret.array = strain.array;
	ret.sort_flag = strain.sort_flag;

	return(ret);
}


STRAIN_ARRAY
stress2strain_array (STRESS_ARRAY stress)
{
	STRAIN_ARRAY ret;

	ret.size = stress.size;
	ret.alloc_size = stress.alloc_size;
	ret.source = stress.source;
	ret.array = stress.array;
	ret.sort_flag = stress.sort_flag;

	return(ret);
}


/* generate suitable stiffness tensor and stress-strain matrix */
RC
fill_material (MATERIAL_PROP *prop)
{
	switch(prop->mat_type){
	case MAT_ISOTROPIC: /* 等方性材料 */
		prop->E_vect.x = prop->E;
		prop->E_vect.y = prop->E;
		prop->E_vect.z = prop->E;
		prop->K_vect.x = prop->K;
		prop->K_vect.y = prop->K;
		prop->K_vect.z = prop->K;
		prop->nu_vect.yz = prop->nu;
		prop->nu_vect.zx = prop->nu;
		prop->nu_vect.xy = prop->nu;
		prop->alpha_vect.x = prop->alpha;
		prop->alpha_vect.y = prop->alpha;
		prop->alpha_vect.z = prop->alpha;
		prop->alpha_vect.yz = 0.0;
		prop->alpha_vect.zx = 0.0;
		prop->alpha_vect.xy = 0.0;
		RC_TRY( put_iso_D_mat(prop) );
		break;
	case MAT_ORTHOTROPIC: /* 直交異方性材料 */
		RC_TRY( put_ortho_D_mat(prop) );
		break;
	case MAT_ANISOTROPIC: /* 異方性材料 */
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	clean_C(prop);
	switch(prop->ana_type){
	case ANA_3D:
		prop->C[0][0][0][0] = prop->D_matrix[0][0];
		prop->C[1][1][1][1] = prop->D_matrix[1][1];
		prop->C[2][2][2][2] = prop->D_matrix[2][2];
		prop->C[0][0][1][1] = prop->D_matrix[0][1];
		prop->C[1][1][0][0] = prop->D_matrix[1][0];
		prop->C[1][1][2][2] = prop->D_matrix[1][2];
		prop->C[2][2][1][1] = prop->D_matrix[2][1];
		prop->C[2][2][0][0] = prop->D_matrix[2][0];
		prop->C[0][0][2][2] = prop->D_matrix[0][2];
		prop->C[0][1][0][1] = prop->C[0][1][1][0] = prop->D_matrix[3][3];
		prop->C[1][0][0][1] = prop->C[1][0][1][0] = prop->D_matrix[3][3];
		prop->C[0][1][1][2] = prop->C[0][1][2][1] = prop->D_matrix[3][4];
		prop->C[1][0][1][2] = prop->C[1][0][2][1] = prop->D_matrix[3][4];
		prop->C[0][1][2][0] = prop->C[0][1][0][2] = prop->D_matrix[3][5];
		prop->C[0][1][2][0] = prop->C[0][1][0][2] = prop->D_matrix[3][5];
		prop->C[1][2][0][1] = prop->C[1][2][1][0] = prop->D_matrix[4][3];
		prop->C[2][1][0][1] = prop->C[2][1][1][0] = prop->D_matrix[4][3];
		prop->C[1][2][1][2] = prop->C[1][2][2][1] = prop->D_matrix[4][4];
		prop->C[2][1][1][2] = prop->C[2][1][2][1] = prop->D_matrix[4][4];
		prop->C[2][1][2][0] = prop->C[2][1][0][2] = prop->D_matrix[4][5];
		prop->C[1][2][2][0] = prop->C[1][2][0][2] = prop->D_matrix[4][5];
		prop->C[2][0][0][1] = prop->C[2][0][1][0] = prop->D_matrix[5][3];
		prop->C[0][2][0][1] = prop->C[0][2][1][0] = prop->D_matrix[5][3];
		prop->C[2][0][1][2] = prop->C[2][0][2][1] = prop->D_matrix[5][4];
		prop->C[0][2][1][2] = prop->C[0][2][2][1] = prop->D_matrix[5][4];
		prop->C[2][0][2][0] = prop->C[2][0][0][2] = prop->D_matrix[5][5];
		prop->C[0][2][2][0] = prop->C[0][2][0][2] = prop->D_matrix[5][5];
		prop->C[0][0][0][1] = prop->C[0][0][1][0] = prop->D_matrix[0][3];
		prop->C[0][1][0][0] = prop->C[1][0][0][0] = prop->D_matrix[3][0];
		prop->C[0][0][1][2] = prop->C[0][0][2][1] = prop->D_matrix[0][4];
		prop->C[1][2][0][0] = prop->C[2][1][0][0] = prop->D_matrix[4][0];
		prop->C[0][0][2][0] = prop->C[0][0][0][2] = prop->D_matrix[0][5];
		prop->C[2][0][0][0] = prop->C[0][2][0][0] = prop->D_matrix[5][0];
		prop->C[1][1][0][1] = prop->C[1][1][1][0] = prop->D_matrix[1][3];
		prop->C[0][1][1][1] = prop->C[1][0][1][1] = prop->D_matrix[3][1];
		prop->C[1][1][1][2] = prop->C[1][1][2][1] = prop->D_matrix[1][4];
		prop->C[1][2][1][1] = prop->C[2][1][1][1] = prop->D_matrix[4][1];
		prop->C[1][1][2][0] = prop->C[1][1][0][2] = prop->D_matrix[1][5];
		prop->C[2][0][1][1] = prop->C[0][2][1][1] = prop->D_matrix[5][1];
		prop->C[2][2][0][1] = prop->C[2][2][1][0] = prop->D_matrix[2][3];
		prop->C[0][1][2][2] = prop->C[1][0][2][2] = prop->D_matrix[3][2];
		prop->C[2][2][1][2] = prop->C[2][2][2][1] = prop->D_matrix[2][4];
		prop->C[1][2][2][2] = prop->C[2][1][2][2] = prop->D_matrix[4][2];
		prop->C[2][2][2][0] = prop->C[2][2][0][2] = prop->D_matrix[2][5];
		prop->C[2][0][2][2] = prop->C[0][2][2][2] = prop->D_matrix[5][2];
		break;
	case ANA_PLANE_STRESS:
	case ANA_PLANE_STRAIN:
	case ANA_PLATE_BEND:
		prop->C[0][0][0][0] = prop->D_matrix[0][0];
		prop->C[0][0][1][1] = prop->D_matrix[0][1];
		prop->C[1][1][0][0] = prop->D_matrix[1][0];
		prop->C[1][1][1][1] = prop->D_matrix[1][1];
		prop->C[0][1][0][1] = prop->C[0][1][1][0] = prop->D_matrix[2][2];
		prop->C[1][0][0][1] = prop->C[1][0][1][0] = prop->D_matrix[2][2];
		prop->C[0][0][0][1] = prop->C[0][0][1][0] = prop->D_matrix[0][2];
		prop->C[0][1][0][0] = prop->C[1][0][0][0] = prop->D_matrix[2][0];
		prop->C[1][1][0][1] = prop->C[1][1][1][0] = prop->D_matrix[1][2];
		prop->C[0][1][1][1] = prop->C[1][0][1][1] = prop->D_matrix[2][1];
		break;
	case ANA_AXISYMMETRY:
		prop->C[0][0][0][0] = prop->D_matrix[0][0];
		prop->C[0][0][1][1] = prop->D_matrix[0][1];
		prop->C[0][0][2][2] = prop->D_matrix[0][2];
		prop->C[1][1][0][0] = prop->D_matrix[1][0];
		prop->C[1][1][1][1] = prop->D_matrix[1][1];
		prop->C[1][1][2][2] = prop->D_matrix[1][2];
		prop->C[2][2][0][0] = prop->D_matrix[2][0];
		prop->C[2][2][1][1] = prop->D_matrix[2][1];
		prop->C[2][2][2][2] = prop->D_matrix[2][2];
		prop->C[0][1][0][1] = prop->C[0][1][1][0] = prop->D_matrix[3][3];
		prop->C[1][0][0][1] = prop->C[1][0][1][0] = prop->D_matrix[3][3];
		prop->C[0][0][0][1] = prop->C[0][0][1][0] = prop->D_matrix[0][3];
		prop->C[0][1][0][0] = prop->C[1][0][0][0] = prop->D_matrix[3][0];
		prop->C[1][1][0][1] = prop->C[1][1][1][0] = prop->D_matrix[1][3];
		prop->C[0][1][1][1] = prop->C[1][0][1][1] = prop->D_matrix[3][1];
		prop->C[2][2][0][1] = prop->C[2][2][1][0] = prop->D_matrix[2][3];
		prop->C[0][1][2][2] = prop->C[1][0][2][2] = prop->D_matrix[3][2];
		break;
	case ANA_PLATE_SHEAR:
		prop->C[0][0][0][0] = prop->D_matrix[0][0];
		prop->C[0][0][1][1] = prop->C[1][1][0][0] = prop->D_matrix[0][1];
		prop->C[1][1][1][1] = prop->D_matrix[1][1];
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}



/* set D_matrix for isotropic material */
static RC
put_iso_D_mat (MATERIAL_PROP *prop)
{
	double fact;
	double div_tmp;

	clean_D_matrix(prop);
	switch(prop->ana_type){
	case ANA_PLANE_STRESS:
		div_tmp = 1.0 - ((prop->nu) * (prop->nu));
		if(nearly_eq(div_tmp, 0.0)) return(CAL_ERROR_RC);
		fact = (prop->E)/div_tmp;
		prop->D_matrix[0][0] = prop->D_matrix[1][1] = fact;
		prop->D_matrix[0][1] = prop->D_matrix[1][0] = fact * (prop->nu);
		prop->D_matrix[2][2] = 0.5 * fact * (1.0 - prop->nu);
		break;
	case ANA_PLANE_STRAIN:
		div_tmp = (1.0 + prop->nu) * (1.0 - 2.0 * prop->nu);
		if(nearly_eq(div_tmp, 0.0)) return(CAL_ERROR_RC);
		fact = (prop->E) * (1.0 - (prop->nu)) / div_tmp;
		div_tmp = 1.0 - prop->nu;
		if(nearly_eq(div_tmp, 0.0)) return(CAL_ERROR_RC);
		prop->D_matrix[0][0] = prop->D_matrix[1][1] = fact;
		prop->D_matrix[0][1] = prop->D_matrix[1][0]
		                     = fact * (prop->nu) / div_tmp;
		prop->D_matrix[2][2] = 0.5 * fact * (1.0 - 2.0*prop->nu) / div_tmp;
		break;
	case ANA_AXISYMMETRY:
		/* z,r,theta,rz */
		div_tmp = (1.0 + prop->nu) * (1.0 - 2.0 * prop->nu);
		if(nearly_eq(div_tmp, 0.0)) return(CAL_ERROR_RC);
		fact = (prop->E) * (1.0 - (prop->nu)) / div_tmp;
		div_tmp = 1.0 - prop->nu;
		if(nearly_eq(div_tmp, 0.0)) return(CAL_ERROR_RC);
		prop->D_matrix[0][0] = prop->D_matrix[1][1]
		                     = prop->D_matrix[2][2] = fact;
		prop->D_matrix[0][1] = prop->D_matrix[1][0] 
		                     = prop->D_matrix[0][2] = prop->D_matrix[2][0] 
		                     = prop->D_matrix[1][2] = prop->D_matrix[2][1] 
		                     = fact * (prop->nu) / div_tmp;
		prop->D_matrix[3][3] = 0.5 * fact * (1.0 - 2.0*prop->nu) / div_tmp;
		break;
	case ANA_3D:
		div_tmp = (1.0 + prop->nu) * (1.0 - 2.0 * prop->nu);
		if(nearly_eq(div_tmp, 0.0)) return(CAL_ERROR_RC);
		fact = (prop->E) * (1.0 - (prop->nu)) / div_tmp;
		div_tmp = 1.0 - prop->nu;
		if(nearly_eq(div_tmp, 0.0)) return(CAL_ERROR_RC);
		prop->D_matrix[0][0] = prop->D_matrix[1][1]
		                     = prop->D_matrix[2][2] = fact;
		prop->D_matrix[0][1] = prop->D_matrix[1][0] 
		                     = prop->D_matrix[0][2] = prop->D_matrix[2][0] 
		                     = prop->D_matrix[1][2] = prop->D_matrix[2][1] 
		                     = fact * (prop->nu) / div_tmp;
		prop->D_matrix[3][3] = prop->D_matrix[4][4] = prop->D_matrix[5][5]
		                     = 0.5 * fact * (1.0 - 2.0*prop->nu) / div_tmp;
		break;
	case ANA_PLATE_BEND:
		div_tmp = 1.0 - ((prop->nu) * (prop->nu));
		if(nearly_eq(div_tmp, 0.0)) return(CAL_ERROR_RC);
		fact = (prop->E) / div_tmp;
		prop->D_matrix[0][0] = prop->D_matrix[1][1] = fact;
		prop->D_matrix[0][1] = prop->D_matrix[1][0] = fact * (prop->nu);
		prop->D_matrix[2][2] = 0.5 * fact * (1.0 - prop->nu);
		break;
	case ANA_PLATE_SHEAR:
		div_tmp = 2.0 * ( 1.0 + (prop->nu) );
		if(nearly_eq(div_tmp, 0.0)) return(CAL_ERROR_RC);
		fact = (prop->E) / div_tmp;
		prop->D_matrix[0][0] = prop->D_matrix[1][1] = fact * SHEAR_FACT;
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


/* set D_matrix for orthotropic material */
static RC
put_ortho_D_mat (MATERIAL_PROP *prop)
{
	double div_tmp;
	double nu_yx;
	double n, m;
	double nu1;
	double fact_tmp;

	clean_D_matrix(prop);
	switch(prop->ana_type){
	case ANA_PLANE_STRESS:
		if(nearly_eq(prop->E_vect.x, 0.0)) return(CAL_ERROR_RC);
		nu_yx = (prop->E_vect.y)*(prop->nu_vect.xy)/(prop->E_vect.x);
		div_tmp = 1.0 - nu_yx * (prop->nu_vect.xy);
		if(nearly_eq(div_tmp, 0.0)) return(CAL_ERROR_RC);
		prop->D_matrix[0][0] = (prop->E_vect.x) / div_tmp;
		prop->D_matrix[1][1] = (prop->E_vect.y) / div_tmp;
		prop->D_matrix[0][1] = prop->D_matrix[1][0]
		                     = (prop->E_vect.x) * nu_yx / div_tmp;
		prop->D_matrix[2][2] = prop->G_vect.xy;
		break;
	case ANA_PLANE_STRAIN:
		if(nearly_eq(prop->E_vect.y, 0.0)) return(CAL_ERROR_RC);
		n = (prop->E_vect.x)/(prop->E_vect.y);
		m = (prop->G_vect.xy)/(prop->E_vect.y);
		if(nearly_eq(prop->G_vect.zx, 0.0)) return(CAL_ERROR_RC);
		nu1 = 0.5 * (prop->E_vect.x)/(prop->G_vect.zx) - 1.0;
		div_tmp = (1 + nu1)
		         *(1 - nu1
		         - 2.0 * n * (prop->nu_vect.xy) * (prop->nu_vect.xy));
		if(nearly_eq(div_tmp, 0.0)) return(CAL_ERROR_RC);
		fact_tmp = (prop->E_vect.y)/div_tmp;
		prop->D_matrix[0][0] = n * (1.0 - n * (prop->nu_vect.xy)
		                                    * (prop->nu_vect.xy)) * fact_tmp;
		prop->D_matrix[0][1] = prop->D_matrix[1][0]
		                     = n * (prop->nu_vect.xy) * (1.0 + nu1) * fact_tmp;
		prop->D_matrix[1][1] = (1.0 - nu1 * nu1) * fact_tmp;
		prop->D_matrix[2][2] = prop->G_vect.xy;
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


RC
Gauss_const_max (ELEM_TYPE type, double **points,
                 double *weight, int *num_point)
{
	switch(type){
	case ELEM_RBAR:
	case ELEM_BEAM1:
	case ELEM_BEAM2:
	case ELEM_LINE1:
	case ELEM_LINE2:
		*num_point = 4;
		RC_TRY( Gauss_const_line(*num_point, points, weight) );
		break;
	case ELEM_TRI1:
	case ELEM_TRI2:
		*num_point = 12;
		RC_TRY( Gauss_const_tri(*num_point, points, weight) );
		break;
	case ELEM_QUAD1:
	case ELEM_QUAD2:
		*num_point = 9;
		RC_TRY( Gauss_const_quad(*num_point, points, weight) );
		break;
	case ELEM_QUAD2L:
		*num_point = 16;
		RC_TRY( Gauss_const_quad(*num_point, points, weight) );
		break;
	case ELEM_TETRA1:
	case ELEM_TETRA2:
		*num_point = 15;
		RC_TRY( Gauss_const_tetra(*num_point, points, weight) );
		break;
	case ELEM_PENTA1:
	case ELEM_PENTA2:
		*num_point = 21;
		RC_TRY( Gauss_const_penta(*num_point, points, weight) );
		break;
	case ELEM_HEXA1:
	case ELEM_HEXA2:
		*num_point = 14;
		RC_TRY( Gauss_const_hexa(*num_point, points, weight) );
		break;
	default :
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


RC
Gauss_const_default (ELEM_TYPE type, double **points,
                     double *weight, int *num_point)
{
	switch(type){
	case ELEM_RBAR:
	case ELEM_BEAM1:
	case ELEM_BEAM2:
	case ELEM_LINE1:
	case ELEM_LINE2:
		*num_point = 3;
		RC_TRY( Gauss_const_line(*num_point, points, weight) );
		break;
	case ELEM_TRI1:
	case ELEM_TRI2:
		*num_point = 3;
		RC_TRY( Gauss_const_tri(*num_point, points, weight) );
		break;
	case ELEM_QUAD1:
		*num_point = 4;
		RC_TRY( Gauss_const_quad(*num_point, points, weight) );
		break;
	case ELEM_QUAD2:
	case ELEM_QUAD2L:
		*num_point = 9;
		RC_TRY( Gauss_const_quad(*num_point, points, weight) );
		break;
	case ELEM_TETRA1:
	case ELEM_TETRA2:
		*num_point = 4;
		RC_TRY( Gauss_const_tetra(*num_point, points, weight) );
		break;
	case ELEM_PENTA1:
		*num_point = 6;
		RC_TRY( Gauss_const_penta(*num_point, points, weight) );
		break;
	case ELEM_PENTA2:
		*num_point = 9;
		RC_TRY( Gauss_const_penta(*num_point, points, weight) );
		break;
	case ELEM_HEXA1:
		*num_point = 8;
		RC_TRY( Gauss_const_hexa(*num_point, points, weight) );
		break;
	case ELEM_HEXA2:
		*num_point = 14;
		RC_TRY( Gauss_const_hexa(*num_point, points, weight) );
		break;
	default :
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


RC
Gauss_const_RI (ELEM_TYPE type, double **points,
                double *weight, int *num_point)
{
	switch(type){
	case ELEM_BEAM1:
		*num_point = 1;
		RC_TRY( Gauss_const_line(*num_point, points, weight) );
		break;
	case ELEM_BEAM2:
		*num_point = 2;
		RC_TRY( Gauss_const_line(*num_point, points, weight) );
		break;
	case ELEM_TRI1:
		*num_point = 1;
		RC_TRY( Gauss_const_tri(*num_point, points, weight) );
		break;
	case ELEM_TRI2:
		*num_point = 3;
		RC_TRY( Gauss_const_tri(*num_point, points, weight) );
		break;
	case ELEM_QUAD1:
		*num_point = 1;
		RC_TRY( Gauss_const_quad(*num_point, points, weight) );
		break;
	case ELEM_QUAD2:
		*num_point = 4;
		RC_TRY( Gauss_const_quad(*num_point, points, weight) );
		break;
	default :
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


RC
Gauss_const_1 (ELEM_TYPE type, double **points, double *weight)
{
	switch(type){
	case ELEM_RBAR:
	case ELEM_BEAM1:
	case ELEM_BEAM2:
	case ELEM_LINE1:
	case ELEM_LINE2:
		RC_TRY( Gauss_const_line(1, points, weight) );
		break;
	case ELEM_TRI1:
	case ELEM_TRI2:
		RC_TRY( Gauss_const_tri(1, points, weight) );
		break;
	case ELEM_QUAD1:
	case ELEM_QUAD2:
	case ELEM_QUAD2L:
		RC_TRY( Gauss_const_quad(1, points, weight) );
		break;
	case ELEM_TETRA1:
	case ELEM_TETRA2:
		RC_TRY( Gauss_const_tetra(1, points, weight) );
		break;
	case ELEM_PENTA1:
	case ELEM_PENTA2:
		RC_TRY( Gauss_const_penta(1, points, weight) );
		break;
	case ELEM_HEXA1:
	case ELEM_HEXA2:
		RC_TRY( Gauss_const_hexa(1, points, weight) );
		break;
	default :
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


/* 問題のひずみ、応力の数を返却 */
int
element_ssnum (ELEM_TYPE type)
{
	switch(type){
	case ELEM_TRI1:
	case ELEM_TRI2:
	case ELEM_QUAD1:
	case ELEM_QUAD2:
	case ELEM_QUAD2L:
		return(3);   /* xx, yy, xy */
	case ELEM_TETRA1:
	case ELEM_TETRA2:
	case ELEM_PENTA1:
	case ELEM_PENTA2:
	case ELEM_HEXA1:
	case ELEM_HEXA2:
		return(6);   /* xx, yy, zz, xy, yz, zx */
	default:
		break;
	}

	return(-1);
}


/* 要素の局所座標系の次元を返却する．
 */
int
element_dim (ELEM_TYPE type)
{
	switch(type){
	case ELEM_LINE1:
	case ELEM_LINE2:
		return(1);
	case ELEM_TRI1:
	case ELEM_TRI2:
	case ELEM_QUAD1:
	case ELEM_QUAD2:
	case ELEM_QUAD2L:
		return(2);
	case ELEM_TETRA1:
	case ELEM_TETRA2:
	case ELEM_PENTA1:
	case ELEM_PENTA2:
	case ELEM_HEXA1:
	case ELEM_HEXA2:
		return(3);
	case ELEM_BEAM1:
	case ELEM_BEAM2:
	case ELEM_SPRING:
		return(6);
	default:
		break;
	}

	return(-1);
}


/* 要素の次数（1次要素か2次要素か）を返却する．
 */
int
element_order (ELEM_TYPE type)
{
	switch(type){
	case ELEM_SHELL1:
	case ELEM_LINE1:
	case ELEM_BEAM1:
	case ELEM_TRI1:
	case ELEM_QUAD1:
	case ELEM_TETRA1:
	case ELEM_PENTA1:
	case ELEM_HEXA1:
		return(1);
	case ELEM_SHELL2:
	case ELEM_LINE2:
	case ELEM_BEAM2:
	case ELEM_TRI2:
	case ELEM_QUAD2:
	case ELEM_QUAD2L:
	case ELEM_TETRA2:
	case ELEM_PENTA2:
	case ELEM_HEXA2:
		return(2);
	default:
		break;
	}

	return(-1);
}


/* 局所座標 local_coord における，要素 element のJacobi行列式を計算する．
 * その過程で dN およびJacobi行列を計算するので，それらも渡す．
 */
RC
cal_det_J (ELEMENT element, const double local_coord[], double **node_location, 
		double *det_J, double **dN, double **Jacobi)
{
	RC_NULL_CHK( node_location );
	RC_NULL_CHK( det_J );

	int dim = element_dim(element.type);
    static double **local_dN = NULL;
    static double **local_Jacobi = NULL;
    static int init_flag = 0;
	VECT3D d_A;

    if(init_flag == 0){
        RC_TRY( allocate2D(3, MAX_NODE, &local_dN) );
        RC_TRY( allocate2D(3, 3, &local_Jacobi) );
        init_flag = 1;
    }

    init_matrix(3, MAX_NODE, local_dN);
    init_matrix(3, 3, local_Jacobi);

	switch(dim){
        case 1:
			RC_TRY( set_dN_matrix(element.type, local_coord, local_dN) );
			mul_matrix(3, element.node_num, 3, local_dN, node_location, 
                    local_Jacobi);
            *det_J = sqrt(inner_product(3, local_Jacobi[0], local_Jacobi[0]));
            break;
		case 2:
			RC_TRY( cal_d_A(element, local_coord, node_location, &d_A, 
                        local_dN, local_Jacobi) );
			*det_J = sqrt(inner_product3d(d_A, d_A));
			break;
		case 3:
			RC_TRY( set_dN_matrix(element.type, local_coord, local_dN) );
			mul_matrix(dim, element.node_num, dim, local_dN, node_location, 
                    local_Jacobi);
			*det_J = determinant(3, local_Jacobi);
			break;
		default:
			return(IMPLEMENT_ERROR_RC);
	}

    if(dN != NULL){
        RC_TRY( copy_matrix(3, MAX_NODE, local_dN, dN) );
    }
    if(Jacobi != NULL){
        RC_TRY( copy_matrix(3, 3, local_Jacobi, Jacobi) );
    }

	return(NORMAL_RC);
}


/* 境界の要素 surf 上の局所座標 local_coord における法線ベクトルを計算する．
 * その過程で dN およびJacobi行列を計算するので，それぞれ
 *   ws_mat_3_N,
 *   ws_mat_3_3
 * で渡す．
 * surf は線要素もしくは平面要素でなければならない．
 *
 * この関数で計算する法線ベクトルの大きさは，surfのJacobi行列式に等しい．
 * したがって，Gauss積分でこの関数で計算した法線ベクトルを用いるときは，
 * Jacobi行列式は必要ない．
 */
RC
cal_d_A (ELEMENT surf, const double local_coord[], double **node_location, 
        VECT3D *d_A, double **dN, double **Jacobi)
{
	int dim;
	VECT3D d_xi, d_eta;

    static double **local_dN = NULL;
    static double **local_Jacobi = NULL;
    static int init_flag = 0;

    if(init_flag == 0){
        RC_TRY( allocate2D(3, MAX_NODE, &local_dN) );
        RC_TRY( allocate2D(3, 3, &local_Jacobi) );
        init_flag = 1;
    }

	/* 前回格納した分を消去しきれない場合があるため，初期化が必要． */
	init_matrix(3, MAX_NODE, local_dN);  
	init_matrix(3, 3, local_Jacobi);  
	RC_TRY( set_dN_matrix(surf.type, local_coord, local_dN) );
	mul_matrix(3, surf.node_num, 3, local_dN, node_location, local_Jacobi);
	dim = element_dim(surf.type);
	if(dim == 2){
		d_xi.x = local_Jacobi[0][0];
		d_xi.y = local_Jacobi[0][1];
		d_xi.z = local_Jacobi[0][2];
		d_eta.x = local_Jacobi[1][0];
		d_eta.y = local_Jacobi[1][1];
		d_eta.z = local_Jacobi[1][2];
		*d_A = outer_product3d(d_xi, d_eta);
	}else if(dim == 1){
		d_A->x = local_Jacobi[0][1];
		d_A->y = -local_Jacobi[0][0];
		d_A->z = 0.0;
	}else{
		return(IMPLEMENT_ERROR_RC);
	}

    if(dN != NULL){
        RC_TRY( copy_matrix(3, MAX_NODE, local_dN, dN) );
    }
    if(Jacobi != NULL){
        RC_TRY( copy_matrix(3, 3, local_Jacobi, Jacobi) );
    }

	return(NORMAL_RC);
}

/* 形状関数の全体座標系における微分 dNxyz を計算する．
 * その過程で計算する dN，Jacobi，Jacobi_inv は，必要であればポインタを
 * 渡すことで格納する．
 */
RC
cal_dNxyz (ELEMENT element, const double local_coord[], double **node_location, 
        double **dNxyz, double **dN, double **Jacobi, double **Jacobi_inv)
{
    RC_NULL_CHK( dNxyz );

    static int init_flag = 0;
    static double **local_dN = NULL;
    static double **local_Jacobi = NULL;
    static double **local_Jacobi_inv = NULL;

    if(init_flag == 0){
        RC_TRY( allocate2D(3, MAX_NODE, &local_dN) );
        RC_TRY( allocate2D(3, 3, &local_Jacobi) );
        RC_TRY( allocate2D(3, 3, &local_Jacobi_inv) );
        init_flag = 1;
    }

    init_matrix(3, MAX_NODE, local_dN); 
    RC_TRY( set_dN_matrix(element.type, local_coord, local_dN) );
    mul_matrix(3, element.node_num, 3, local_dN, node_location, local_Jacobi);
    if(element_dim(element.type) == 2){
        local_Jacobi[2][0] = local_Jacobi[0][1]*local_Jacobi[1][2] 
            - local_Jacobi[0][2]*local_Jacobi[1][1];
        local_Jacobi[2][1] = local_Jacobi[0][2]*local_Jacobi[1][0] 
            - local_Jacobi[0][0]*local_Jacobi[1][2];
        local_Jacobi[2][2] = local_Jacobi[0][0]*local_Jacobi[1][1] 
            - local_Jacobi[0][1]*local_Jacobi[1][0];
    }
    RC_TRY( inverse_matrix(3, local_Jacobi, local_Jacobi_inv) );
    mul_matrix(3, 3, element.node_num, local_Jacobi_inv, local_dN, dNxyz);

    if(dN != NULL){
        RC_TRY( copy_matrix(3, element.node_num, local_dN, dN) );
    }
    if(Jacobi != NULL){
        RC_TRY( copy_matrix(3, 3, local_Jacobi, Jacobi) );
    }
    if(Jacobi_inv != NULL){
        RC_TRY( copy_matrix(3, 3, local_Jacobi_inv, Jacobi_inv) );
    }

    return(NORMAL_RC);
}


/* 要素 element 上の局所座標 local_coord における変位 disp を内挿する．
 * disp_matrix は fill_disp_matrix() で計算する．
 */
RC
interpolate_disp (ELEMENT element, const double local_coord[],
                  double **disp_matrix, VECT3DR *disp)
{
	int ii1;
	double N[MAX_NODE];

	disp->x = 0.0;
	disp->y = 0.0;
	disp->z = 0.0;
	disp->yz = 0.0;
	disp->zx = 0.0;
	disp->xy = 0.0;

	RC_TRY( set_N_vector(element.type, local_coord, N) );
	for(ii1=0; ii1<element.node_num; ii1++){
		disp->x += N[ii1] * disp_matrix[ii1][0];
		disp->y += N[ii1] * disp_matrix[ii1][1];
		disp->z += N[ii1] * disp_matrix[ii1][2];
		disp->yz += N[ii1] * disp_matrix[ii1][3];
		disp->zx += N[ii1] * disp_matrix[ii1][4];
		disp->xy += N[ii1] * disp_matrix[ii1][5];
	}

	return(NORMAL_RC);
}


RC
realloc_node_matrix (NODE_MATRIX *nmat)
{
	int ii1;
	int add_size;

	if(nmat->alloc_size <= 0){
		nmat->alloc_size = 0;
		nmat->size = 0;
	}

	if( (nmat->alloc_size > 0)&&(nmat->alloc_size > nmat->size) ){
		return(NORMAL_RC);
	}

	add_size = nmat->alloc_size;
	if(add_size <= 0) add_size = 1;
	if(add_size > 256) add_size = 256;

	nmat->alloc_size += add_size;

	nmat->col_node = mm_realloc(nmat->col_node,
	                            (nmat->alloc_size)*sizeof(int));
	if(nmat->col_node == NULL) return(ALLOC_ERROR_RC);

	nmat->col_shift = mm_realloc(nmat->col_shift,
	                             (nmat->alloc_size)*sizeof(int));
	if(nmat->col_shift == NULL) return(ALLOC_ERROR_RC);

	nmat->matrix = mm_realloc(nmat->matrix,
	                          (nmat->alloc_size)*sizeof(double[3][3]));
	if(nmat->matrix == NULL) return(ALLOC_ERROR_RC);

	for(ii1=nmat->size; ii1<(nmat->alloc_size); ii1++){
		nmat->col_node[ii1] = -1;
		nmat->col_shift[ii1] = -1;
		nmat->matrix[ii1][0][0] = 0.0;
		nmat->matrix[ii1][0][1] = 0.0;
		nmat->matrix[ii1][0][2] = 0.0;
		nmat->matrix[ii1][1][0] = 0.0;
		nmat->matrix[ii1][1][1] = 0.0;
		nmat->matrix[ii1][1][2] = 0.0;
		nmat->matrix[ii1][2][0] = 0.0;
		nmat->matrix[ii1][2][1] = 0.0;
		nmat->matrix[ii1][2][2] = 0.0;
	}

	return(NORMAL_RC);
}


