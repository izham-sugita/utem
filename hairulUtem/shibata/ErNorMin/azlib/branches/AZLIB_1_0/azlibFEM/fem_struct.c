/*********************************************************************
 * fem_struct.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Tomoyuki TSUBATA>
 *
 * Thanks to following contributors.
 *  <Yasuo ASAGA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: fem_struct.c,v 1.26 2003/07/22 11:44:01 adachi Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rc.h"
#include "math_utl.h"
#include "fem_struct.h"

#define SHEAR_FACT (1.2)     /* $BEy2A$;$sCG78?t(B($B$^$?$O=$@5$;$sCG78?t(B) */

/* local functions */
static void clean_D_matrix(FEM_MATERIAL_PROP *prop);
static void clean_C(FEM_MATERIAL_PROP *prop);
static RC put_iso_D_mat(FEM_MATERIAL_PROP *prop);
static RC put_ortho_D_mat(FEM_MATERIAL_PROP *prop);
static RC set_fem_face_table_tri1(FEM_FACE_TABLE *f_table);
static RC set_fem_face_table_tri2(FEM_FACE_TABLE *f_table);
static RC set_fem_face_table_quad1(FEM_FACE_TABLE *f_table);
static RC set_fem_face_table_quad2(FEM_FACE_TABLE *f_table);
static RC set_fem_face_table_tetra1(FEM_FACE_TABLE *f_table);
static RC set_fem_face_table_tetra2(FEM_FACE_TABLE *f_table);
static RC set_fem_face_table_penta1(FEM_FACE_TABLE *f_table);
static RC set_fem_face_table_penta2(FEM_FACE_TABLE *f_table);
static RC set_fem_face_table_hexa1(FEM_FACE_TABLE *f_table);
static RC set_fem_face_table_hexa2(FEM_FACE_TABLE *f_table);
static RC set_fem_edge_table_tri1(FEM_EDGE_TABLE *e_table);
static RC set_fem_edge_table_tri2(FEM_EDGE_TABLE *e_table);
static RC set_fem_edge_table_quad1(FEM_EDGE_TABLE *e_table);
static RC set_fem_edge_table_quad2(FEM_EDGE_TABLE *e_table);
static RC set_fem_edge_table_tetra1(FEM_EDGE_TABLE *e_table);
static RC set_fem_edge_table_tetra2(FEM_EDGE_TABLE *e_table);
static RC set_fem_edge_table_penta1(FEM_EDGE_TABLE *e_table);
static RC set_fem_edge_table_penta2(FEM_EDGE_TABLE *e_table);
static RC set_fem_edge_table_hexa1(FEM_EDGE_TABLE *e_table);
static RC set_fem_edge_table_hexa2(FEM_EDGE_TABLE *e_table);

/* local variables */
static TRANSLATION3D ZeroTranslation3d = {0.0, 0.0, 0.0};
static ROTATION3D ZeroRotation3d = {0.0, 0.0, 0.0};
static TRANS_ROTATION3D ZeroTransRotation3d = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
static I_TRANS_ROTATION3D ZeroITransRotation3d = {0, 0, 0, 0, 0, 0};


/* id1[num1], scale1[num1] $B$H(B id2[num2], scale2[num2] $BHf3S(B */
int same_bc_set(int num1, const int id1[], double scale1[],
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


/* $BMWAG(B s_elem $B$NA4$F$N@aE@$,!"(Bv_elem $B$K$bB8:_$9$k;~(B */
/* s_elem $B>e$N6I=j:BI8(B s_coords[coord_size][4] $B$r(B       */
/* v_elem $B>e$N6I=j:BI8(B v_coords[coord_size][4] $B$KJQ49(B   */
RC trans_coord_s2v(FEM_ELEMENT s_elem, FEM_ELEMENT v_elem,
                   int coord_size, double **s_coords, double **v_coords)
{
	int ii1, ii2, ii3;
	double **v_node_coord;
	double s_N[FEM_MAX_NODE];
	int v_index;

	RC_TRY( allocate2D(FEM_MAX_NODE, 4, &v_node_coord) );
	RC_TRY( set_local_node_points(v_elem.type, v_node_coord) );

	for(ii1=0; ii1<coord_size; ii1++){
		for(ii2=0; ii2<4; ii2++){
			v_coords[ii1][ii2] = 0.0;
		}
		RC_TRY( set_N_vector(s_elem.type, s_coords[ii1], s_N) );

		/* v_elem $B$N3F@aE@$G$N6I=j:BI8CM$r(B s_elem $BFb$GFbA^(B */
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

	RC_TRY( free2D(FEM_MAX_NODE, 4, v_node_coord) );

	return(NORMAL_RC);
}


RC allocate_fem_work_set(FEM_WORK_SET *ws)
{
	int ii1;

	for(ii1=0; ii1<FEM_WS_SIZE; ii1++){
		RC_TRY( allocate1D(FEM_MAX_NODE, &(ws->vec_N[ii1])) );
		RC_TRY( allocate1D(MAX_GAUSS_POINT, &(ws->vec_G[ii1])) );

		RC_TRY( allocate2D(3, 3, &(ws->mat_3_3[ii1])) );
		RC_TRY( allocate2D(6, 6, &(ws->mat_6_6[ii1])) );
		RC_TRY( allocate2D(9, 9, &(ws->mat_9_9[ii1])) );

		RC_TRY( allocate2D(3, FEM_MAX_NODE, &(ws->mat_3_N[ii1])) );
		RC_TRY( allocate2D(3, 3*FEM_MAX_NODE, &(ws->mat_3_3N[ii1])) );
		RC_TRY( allocate2D(6, 3*FEM_MAX_NODE, &(ws->mat_6_3N[ii1])) );
		RC_TRY( allocate2D(9, 3*FEM_MAX_NODE, &(ws->mat_9_3N[ii1])) );
		RC_TRY( allocate2D(3*FEM_MAX_NODE, 3*FEM_MAX_NODE,
		                                      &(ws->mat_3N_3N[ii1])) );

		RC_TRY( allocate2D(FEM_MAX_NODE, 3, &(ws->mat_N_3[ii1])) );
		RC_TRY( allocate2D(FEM_MAX_NODE, 4, &(ws->mat_N_4[ii1])) );
		RC_TRY( allocate2D(FEM_MAX_NODE, 6, &(ws->mat_N_6[ii1])) );
		RC_TRY( allocate2D(MAX_GAUSS_POINT, 4, &(ws->mat_G_4[ii1])) );
	}

	return(NORMAL_RC);
}


RC free_fem_work_set(FEM_WORK_SET *ws)
{
	int ii1;

	for(ii1=0; ii1<FEM_WS_SIZE; ii1++){
		free(ws->vec_N[ii1]);
		free(ws->vec_G[ii1]);
		ws->vec_N[ii1] = NULL;
		ws->vec_G[ii1] = NULL;

		RC_TRY( free2D(3, 3, (ws->mat_3_3[ii1])) );
		RC_TRY( free2D(6, 6, (ws->mat_6_6[ii1])) );
		RC_TRY( free2D(9, 9, (ws->mat_9_9[ii1])) );
		ws->mat_3_3[ii1] = NULL;
		ws->mat_6_6[ii1] = NULL;
		ws->mat_9_9[ii1] = NULL;

		RC_TRY( free2D(3, FEM_MAX_NODE, (ws->mat_3_N[ii1])) );
		RC_TRY( free2D(3, 3*FEM_MAX_NODE, (ws->mat_3_3N[ii1])) );
		RC_TRY( free2D(6, 3*FEM_MAX_NODE, (ws->mat_6_3N[ii1])) );
		RC_TRY( free2D(9, 3*FEM_MAX_NODE, (ws->mat_9_3N[ii1])) );
		RC_TRY( free2D(3*FEM_MAX_NODE, 3*FEM_MAX_NODE,
		                                  (ws->mat_3N_3N[ii1])) );
		ws->mat_3_N[ii1] = NULL;
		ws->mat_3_3N[ii1] = NULL;
		ws->mat_6_3N[ii1] = NULL;
		ws->mat_9_3N[ii1] = NULL;
		ws->mat_3N_3N[ii1] = NULL;

		RC_TRY( free2D(FEM_MAX_NODE, 3, (ws->mat_N_3[ii1])) );
		RC_TRY( free2D(FEM_MAX_NODE, 4, (ws->mat_N_4[ii1])) );
		RC_TRY( free2D(FEM_MAX_NODE, 6, (ws->mat_N_6[ii1])) );
		RC_TRY( free2D(MAX_GAUSS_POINT, 4, (ws->mat_G_4[ii1])) );
		ws->mat_N_3[ii1] = NULL;
		ws->mat_N_4[ii1] = NULL;
		ws->mat_N_6[ii1] = NULL;
		ws->mat_G_4[ii1] = NULL;
	}

	return(NORMAL_RC);
}


/* initialization of FEM_NODE */
void init_fem_node(FEM_NODE *node)
{
	int ii1;
	
	node->label = -1;
	node->p = ZeroTranslation3d;
	node->renum = -1;
	for(ii1=0;ii1<FEM_NODE_INFO;ii1++){
		node->i_info[ii1] = 0;
		node->d_info[ii1] = 0.0;
	}
}


/* initialization of FEM_ELEMENT */
void init_fem_element(FEM_ELEMENT *elem)
{
	int ii;

	elem->label = -1;
	elem->type = ELEM_TRI1;
	elem->material = -1;
	elem->physical = -1;
	elem->node_num = 0;
	for(ii=0;ii<FEM_MAX_NODE;ii++){
		elem->node[ii] = -1;
	}
	elem->volume = 0.0;
	for(ii=0; ii<FEM_ELEMENT_INFO; ii++){
		elem->i_info[ii] = 0;
		elem->d_info[ii] = 0.0;
	}
}


/* initialization of FEM_RIGID_ELEMENT */
void init_fem_rigid_element(FEM_RIGID_ELEMENT *elem)
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
	for(ii=0; ii<FEM_RIGID_ELEMENT_INFO; ii++){
		elem->i_info[ii] = 0;
		elem->d_info[ii] = 0.0;
	}
}


RC free_fem_rigid_element(FEM_RIGID_ELEMENT *elem)
{
	if( (elem->label >= 0)&&(elem->node_num > 0) ){
		if(elem->node != NULL) free(elem->node);
		if(elem->weight != NULL) free(elem->weight);
		if(elem->dof != NULL) free(elem->dof);
	}
	init_fem_rigid_element(elem);

	return(NORMAL_RC);
}


RC copy_fem_rigid_element(FEM_RIGID_ELEMENT src, FEM_RIGID_ELEMENT *dest)
{
	int ii1;

	*dest = src;
	if(src.node_num > 0){
		RC_NULL_CHK( dest->node = malloc(src.node_num * sizeof(int)) );
		for(ii1=0; ii1<src.node_num; ii1++){
			dest->node[ii1] = src.node[ii1];
		}
	}else{
		dest->node = NULL;
	}

	if(src.weight_num > 0){
		RC_NULL_CHK( dest->weight = malloc(src.weight_num * sizeof(double)) );
		for(ii1=0; ii1<src.weight_num; ii1++){
			dest->weight[ii1] = src.weight[ii1];
		}
	}else{
		dest->weight = NULL;
	}

	if(src.dof_num > 0){
		RC_NULL_CHK( dest->dof = malloc(src.dof_num * sizeof(int)) );
		for(ii1=0; ii1<src.dof_num; ii1++){
			dest->dof[ii1] = src.dof[ii1];
		}
	}else{
		dest->dof = NULL;
	}

	return(NORMAL_RC);
}


/* initialization of FEM_MATERIAL_PROP */
void init_fem_material_prop(FEM_MATERIAL_PROP *prop)
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

	prop->E_vect = ZeroTranslation3d;
	prop->G_vect = ZeroRotation3d;
	prop->nu_vect = ZeroRotation3d;
	prop->alpha_vect = ZeroTransRotation3d;

	clean_D_matrix(prop);

	clean_C(prop);
	prop->rho = 0.0;
	prop->tref = 0.0;
	for(ii1=0;ii1<FEM_MATERIAL_INFO;ii1++){
		prop->i_info[ii1] = 0;
		prop->d_info[ii1] = 0.0;
	}
}


/* initialization of D_matrix */
static void clean_D_matrix(FEM_MATERIAL_PROP *prop)
{
	int ii1,ii2;

	for(ii1=0;ii1<6;ii1++){
		for(ii2=0;ii2<6;ii2++){
			prop->D_matrix[ii1][ii2] = 0.0;
		}
	}
}


/* initialization of C of material */
static void clean_C(FEM_MATERIAL_PROP *prop)
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


/* initialization of FEM_PHYSICAL_PROP */
void init_fem_physical_prop(FEM_PHYSICAL_PROP *prop)
{
	int ii1;

	prop->label = -1;
	prop->section = SECTION_ROD;
	for(ii1=0;ii1<FEM_PHYSICAL_INFO;ii1++){
		prop->i_info[ii1] = 0;
		prop->d_info[ii1] = 0.0;
	}
}


/* initialization of FEM_BC */
void init_fem_bc(FEM_BC *bc)
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
	bc->v = ZeroTransRotation3d;
	bc->v_nl = ZeroITransRotation3d;
	for(ii1=0;ii1<FEM_BC_SCALAR;ii1++){
		bc->s_type[ii1] = BC_FREE;
		bc->s[ii1] = 0.0;
		bc->s_nl[ii1] = 0;
	}
	for(ii1=0;ii1<FEM_BC_INFO;ii1++){
		bc->i_info[ii1] = 0;
		bc->d_info[ii1] = 0.0;
	}
}


/* initialization of FEM_ELEM_BC */
void init_fem_elem_bc(FEM_ELEM_BC *bc)
{
	int ii1;

	bc->element = -1;
	bc->set_id = -1;
	bc->type = BC_FREE;
	bc->face_id = -1;
	for(ii1=0; ii1<FEM_MAX_NODE; ii1++){
		bc->val.ns[ii1] = 0.0;
	}
	bc->dir.x = bc->dir.y = bc->dir.z = 0.0;
	for(ii1=0;ii1<FEM_ELEM_BC_INFO;ii1++){
		bc->i_info[ii1] = 0;
		bc->d_info[ii1] = 0.0;
	}
}


/* initialization of FEM_DISP_VELO */
void init_fem_disp_velo(FEM_DISP_VELO *disp)
{
	int ii1;

	disp->node = -1;
	disp->v = ZeroTransRotation3d;
	disp->s = 0.0;
	for(ii1=0;ii1<FEM_DISP_VELO_INFO;ii1++){
		disp->i_info[ii1] = 0;
		disp->d_info[ii1] = 0.0;
	}
}


/* initialization of FEM_STRESS_STRAIN */
void init_fem_stress_strain(FEM_STRESS_STRAIN *ss)
{
	int ii1;

	ss->element = -1;
	ss->node = -1;
	ss->v = ZeroTransRotation3d;
	for(ii1=0;ii1<FEM_STRESS_STRAIN_INFO;ii1++){
		ss->i_info[ii1] = 0;
		ss->d_info[ii1] = 0.0;
	}
}


void init_fem_sensitivity(FEM_SENSITIVITY *sens)
{
	sens->element = -1;
	sens->node = -1;
	sens->s = 0.0;
	sens->v = ZeroTransRotation3d;
}


void init_fem_default_bc(FEM_DEFAULT_BC *def_bc)
{
	int ii1;

	def_bc->set_id = -1;
	def_bc->temp = 0.0;
	def_bc->grav = ZeroTranslation3d;

	for(ii1=0; ii1<FEM_DEFAULT_BC_INFO; ii1++){
		def_bc->i_info[ii1] = 0;
		def_bc->d_info[ii1] = 0.0;
	}
}


void init_fem_bc_set(FEM_BC_SET *bc_set)
{
	int ii1;

	bc_set->label = -1;
	bc_set->num_fix = 0;
	bc_set->num_mpc = 0;
	bc_set->num_force = 0;
	bc_set->num_temp = 0;

	for(ii1=0; ii1<FEM_MAX_BC_SET; ii1++){
		bc_set->id_fix[ii1] = -1;
		bc_set->id_mpc[ii1] = -1;
		bc_set->id_force[ii1] = -1;
		bc_set->id_temp[ii1] = -1;

		bc_set->scale_fix[ii1] = 0.0;
		bc_set->scale_mpc[ii1] = 0.0;
		bc_set->scale_force[ii1] = 0.0;
		bc_set->scale_temp[ii1] = 0.0;
	}

	for(ii1=0; ii1<FEM_BC_SET_INFO; ii1++){
		bc_set->i_info[ii1] = 0;
		bc_set->d_info[ii1] = 0.0;
	}
	init_string_array(&(bc_set->s_info));
}


RC free_fem_bc_set(FEM_BC_SET *bc_set)
{
	RC_TRY( free_string_array(&(bc_set->s_info)) );
	init_fem_bc_set(bc_set);

	return(NORMAL_RC);
}


void init_fem_sound_pressure(FEM_SOUND_PRESSURE *sound)
{
	int ii1;

	sound->node = -1;
	sound->re = 0.0;
	sound->im = 0.0;

	for(ii1=0; ii1<FEM_SOUND_PRESSURE_INFO; ii1++){
		sound->i_info[ii1] = 0;
		sound->d_info[ii1] = 0.0;
	}
}


void init_fem_load_curve(FEM_LOAD_CURVE *curve)
{
	int ii1;

	curve->label = -1;
	curve->max_time = 0.0;
	curve->max_value = 0.0;
	curve->total_step = 0;
	for(ii1=0; ii1<FEM_LOAD_CURVE_MAX_INC; ii1++){
		curve->inc_value[ii1] = 0.0;
	}
	curve->approx_size = 0;
	for(ii1=0; ii1<FEM_LOAD_CURVE_MAX_APPROX; ii1++){
		curve->approx_time[ii1] = 0.0;
		curve->approx_value[ii1] = 0.0;
	}
}


void init_fem_local_coord(FEM_LOCAL_COORD *local_coord)
{
	int ii1;

	local_coord->label = -1;
	local_coord->origin = ZeroTranslation3d;
	local_coord->d_cos[0] = ZeroTranslation3d;
	local_coord->d_cos[1] = ZeroTranslation3d;
	local_coord->d_cos[2] = ZeroTranslation3d;
	for(ii1=0; ii1<FEM_LOCAL_COORD_INFO; ii1++){
		local_coord->i_info[ii1] = 0;
		local_coord->d_info[ii1] = 0.0;
	}
}


void init_fem_elem_matrix(FEM_ELEM_MATRIX *elem_matrix)
{
	elem_matrix->element = -1;
	elem_matrix->size = 0;
	elem_matrix->dim = 0;
	elem_matrix->index = NULL;
	elem_matrix->matrix = NULL;
}


void init_fem_face_table(FEM_FACE_TABLE *f_table)
{
	int ii1, ii2;

	f_table->face_num = 0;
	for(ii1=0; ii1<FEM_MAX_FACE; ii1++){
		f_table->node_num[ii1] = 0;
		f_table->face_type[ii1] = ELEM_LINE1;
		for(ii2=0; ii2<FEM_MAX_FACE_NODE; ii2++){
			f_table->table[ii1][ii2] = 0;
		}
	}
}


RC make_fem_face_table(FEM_ELEM_TYPE type, FEM_FACE_TABLE *f_table)
{
	init_fem_face_table(f_table);

	switch(type){
	case ELEM_TRI1:
		RC_TRY( set_fem_face_table_tri1(f_table) );
		break;
	case ELEM_TRI2:
		RC_TRY( set_fem_face_table_tri2(f_table) );
		break;
	case ELEM_QUAD1:
		RC_TRY( set_fem_face_table_quad1(f_table) );
		break;
	case ELEM_QUAD2:
		RC_TRY( set_fem_face_table_quad2(f_table) );
		break;
	case ELEM_TETRA1:
		RC_TRY( set_fem_face_table_tetra1(f_table) );
		break;
	case ELEM_TETRA2:
		RC_TRY( set_fem_face_table_tetra2(f_table) );
		break;
	case ELEM_PENTA1:
		RC_TRY( set_fem_face_table_penta1(f_table) );
		break;
	case ELEM_PENTA2:
		RC_TRY( set_fem_face_table_penta2(f_table) );
		break;
	case ELEM_HEXA1:
		RC_TRY( set_fem_face_table_hexa1(f_table) );
		break;
	case ELEM_HEXA2:
		RC_TRY( set_fem_face_table_hexa2(f_table) );
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC set_fem_face_table_tri1(FEM_FACE_TABLE *f_table)
{
	f_table->face_num = 3;

	f_table->node_num[0] = 2;
	f_table->face_type[0] = ELEM_LINE1;
	f_table->table[0][0] = 0;
	f_table->table[0][1] = 1;

	f_table->node_num[1] = 2;
	f_table->face_type[1] = ELEM_LINE1;
	f_table->table[1][0] = 1;
	f_table->table[1][1] = 2;

	f_table->node_num[2] = 2;
	f_table->face_type[2] = ELEM_LINE1;
	f_table->table[2][0] = 2;
	f_table->table[2][1] = 0;

	return(NORMAL_RC);
}


static RC set_fem_face_table_tri2(FEM_FACE_TABLE *f_table)
{
	f_table->face_num = 3;

	f_table->node_num[0] = 3;
	f_table->face_type[0] = ELEM_LINE2;
	f_table->table[0][0] = 0;
	f_table->table[0][1] = 1;
	f_table->table[0][2] = 5;

	f_table->node_num[1] = 3;
	f_table->face_type[1] = ELEM_LINE2;
	f_table->table[1][0] = 1;
	f_table->table[1][1] = 2;
	f_table->table[1][2] = 3;

	f_table->node_num[2] = 3;
	f_table->face_type[2] = ELEM_LINE2;
	f_table->table[2][0] = 2;
	f_table->table[2][1] = 0;
	f_table->table[2][2] = 4;

	return(NORMAL_RC);
}


static RC set_fem_face_table_quad1(FEM_FACE_TABLE *f_table)
{
	f_table->face_num = 4;

	f_table->node_num[0] = 2;
	f_table->face_type[0] = ELEM_LINE1;
	f_table->table[0][0] = 0;
	f_table->table[0][1] = 1;

	f_table->node_num[1] = 2;
	f_table->face_type[1] = ELEM_LINE1;
	f_table->table[1][0] = 1;
	f_table->table[1][1] = 2;

	f_table->node_num[2] = 2;
	f_table->face_type[2] = ELEM_LINE1;
	f_table->table[2][0] = 2;
	f_table->table[2][1] = 3;

	f_table->node_num[3] = 2;
	f_table->face_type[3] = ELEM_LINE1;
	f_table->table[3][0] = 3;
	f_table->table[3][1] = 0;

	return(NORMAL_RC);
}


static RC set_fem_face_table_quad2(FEM_FACE_TABLE *f_table)
{
	f_table->face_num = 4;

	f_table->node_num[0] = 3;
	f_table->face_type[0] = ELEM_LINE2;
	f_table->table[0][0] = 0;
	f_table->table[0][1] = 1;
	f_table->table[0][2] = 4;

	f_table->node_num[1] = 3;
	f_table->face_type[1] = ELEM_LINE2;
	f_table->table[1][0] = 1;
	f_table->table[1][1] = 2;
	f_table->table[1][2] = 5;

	f_table->node_num[2] = 3;
	f_table->face_type[2] = ELEM_LINE2;
	f_table->table[2][0] = 2;
	f_table->table[2][1] = 3;
	f_table->table[2][2] = 6;

	f_table->node_num[3] = 3;
	f_table->face_type[3] = ELEM_LINE2;
	f_table->table[3][0] = 3;
	f_table->table[3][1] = 0;
	f_table->table[3][2] = 7;

	return(NORMAL_RC);
}


static RC set_fem_face_table_tetra1(FEM_FACE_TABLE *f_table)
{
	f_table->face_num = 4;

	f_table->node_num[0] = 3;
	f_table->face_type[0] = ELEM_TRI1;
	f_table->table[0][0] = 0;
	f_table->table[0][1] = 3;
	f_table->table[0][2] = 1;

	f_table->node_num[1] = 3;
	f_table->face_type[1] = ELEM_TRI1;
	f_table->table[1][0] = 1;
	f_table->table[1][1] = 3;
	f_table->table[1][2] = 2;

	f_table->node_num[2] = 3;
	f_table->face_type[2] = ELEM_TRI1;
	f_table->table[2][0] = 0;
	f_table->table[2][1] = 1;
	f_table->table[2][2] = 2;

	f_table->node_num[3] = 3;
	f_table->face_type[3] = ELEM_TRI1;
	f_table->table[3][0] = 2;
	f_table->table[3][1] = 3;
	f_table->table[3][2] = 0;

	return(NORMAL_RC);
}


static RC set_fem_face_table_tetra2(FEM_FACE_TABLE *f_table)
{
	f_table->face_num = 4;

	f_table->node_num[0] = 6;
	f_table->face_type[0] = ELEM_TRI2;
	f_table->table[0][0] = 0;
	f_table->table[0][1] = 3;
	f_table->table[0][2] = 1;
	f_table->table[0][3] = 8;
	f_table->table[0][4] = 6;
	f_table->table[0][5] = 7;

	f_table->node_num[1] = 6;
	f_table->face_type[1] = ELEM_TRI2;
	f_table->table[1][0] = 1;
	f_table->table[1][1] = 3;
	f_table->table[1][2] = 2;
	f_table->table[1][3] = 9;
	f_table->table[1][4] = 4;
	f_table->table[1][5] = 8;

	f_table->node_num[2] = 6;
	f_table->face_type[2] = ELEM_TRI2;
	f_table->table[2][0] = 0;
	f_table->table[2][1] = 1;
	f_table->table[2][2] = 2;
	f_table->table[2][3] = 4;
	f_table->table[2][4] = 5;
	f_table->table[2][5] = 6;

	f_table->node_num[3] = 6;
	f_table->face_type[3] = ELEM_TRI2;
	f_table->table[3][0] = 2;
	f_table->table[3][1] = 3;
	f_table->table[3][2] = 0;
	f_table->table[3][3] = 7;
	f_table->table[3][4] = 5;
	f_table->table[3][5] = 9;

	return(NORMAL_RC);
}


static RC set_fem_face_table_penta1(FEM_FACE_TABLE *f_table)
{
	f_table->face_num = 5;

	f_table->node_num[0] = 3;
	f_table->face_type[0] = ELEM_TRI1;
	f_table->table[0][0] = 1;
	f_table->table[0][1] = 0;
	f_table->table[0][2] = 2;

	f_table->node_num[1] = 3;
	f_table->face_type[1] = ELEM_TRI1;
	f_table->table[1][0] = 3;
	f_table->table[1][1] = 4;
	f_table->table[1][2] = 5;

	f_table->node_num[2] = 4;
	f_table->face_type[2] = ELEM_QUAD1;
	f_table->table[2][0] = 1;
	f_table->table[2][1] = 2;
	f_table->table[2][2] = 5;
	f_table->table[2][3] = 4;

	f_table->node_num[3] = 4;
	f_table->face_type[3] = ELEM_QUAD1;
	f_table->table[3][0] = 0;
	f_table->table[3][1] = 3;
	f_table->table[3][2] = 5;
	f_table->table[3][3] = 2;

	f_table->node_num[4] = 4;
	f_table->face_type[4] = ELEM_QUAD1;
	f_table->table[4][0] = 0;
	f_table->table[4][1] = 1;
	f_table->table[4][2] = 4;
	f_table->table[4][3] = 3;

	return(NORMAL_RC);
}


static RC set_fem_face_table_penta2(FEM_FACE_TABLE *f_table)
{
	f_table->face_num = 5;

	f_table->node_num[0] = 6;
	f_table->face_type[0] = ELEM_TRI2;
	f_table->table[0][0] = 1;
	f_table->table[0][1] = 0;
	f_table->table[0][2] = 2;
	f_table->table[0][3] = 7;
	f_table->table[0][4] = 6;
	f_table->table[0][5] = 8;

	f_table->node_num[1] = 6;
	f_table->face_type[1] = ELEM_TRI2;
	f_table->table[1][0] = 3;
	f_table->table[1][1] = 4;
	f_table->table[1][2] = 5;
	f_table->table[1][3] = 9;
	f_table->table[1][4] = 10;
	f_table->table[1][5] = 11;

	f_table->node_num[2] = 8;
	f_table->face_type[2] = ELEM_QUAD2;
	f_table->table[2][0] = 1;
	f_table->table[2][1] = 2;
	f_table->table[2][2] = 5;
	f_table->table[2][3] = 4;
	f_table->table[2][4] = 6;
	f_table->table[2][5] = 14;
	f_table->table[2][6] = 9;
	f_table->table[2][7] = 13;

	f_table->node_num[3] = 8;
	f_table->face_type[3] = ELEM_QUAD2;
	f_table->table[3][0] = 0;
	f_table->table[3][1] = 3;
	f_table->table[3][2] = 5;
	f_table->table[3][3] = 2;
	f_table->table[3][4] = 12;
	f_table->table[3][5] = 10;
	f_table->table[3][6] = 14;
	f_table->table[3][7] = 7;

	f_table->node_num[4] = 8;
	f_table->face_type[4] = ELEM_QUAD2;
	f_table->table[4][0] = 0;
	f_table->table[4][1] = 1;
	f_table->table[4][2] = 4;
	f_table->table[4][3] = 3;
	f_table->table[4][4] = 8;
	f_table->table[4][5] = 13;
	f_table->table[4][6] = 11;
	f_table->table[4][7] = 12;

	return(NORMAL_RC);
}


static RC set_fem_face_table_hexa1(FEM_FACE_TABLE *f_table)
{
	f_table->face_num = 6;

	f_table->node_num[0] = 4;
	f_table->face_type[0] = ELEM_QUAD1;
	f_table->table[0][0] = 0;
	f_table->table[0][1] = 3;
	f_table->table[0][2] = 2;
	f_table->table[0][3] = 1;

	f_table->node_num[1] = 4;
	f_table->face_type[1] = ELEM_QUAD1;
	f_table->table[1][0] = 0;
	f_table->table[1][1] = 1;
	f_table->table[1][2] = 5;
	f_table->table[1][3] = 4;

	f_table->node_num[2] = 4;
	f_table->face_type[2] = ELEM_QUAD1;
	f_table->table[2][0] = 1;
	f_table->table[2][1] = 2;
	f_table->table[2][2] = 6;
	f_table->table[2][3] = 5;

	f_table->node_num[3] = 4;
	f_table->face_type[3] = ELEM_QUAD1;
	f_table->table[3][0] = 7;
	f_table->table[3][1] = 4;
	f_table->table[3][2] = 5;
	f_table->table[3][3] = 6;

	f_table->node_num[4] = 4;
	f_table->face_type[4] = ELEM_QUAD1;
	f_table->table[4][0] = 3;
	f_table->table[4][1] = 0;
	f_table->table[4][2] = 4;
	f_table->table[4][3] = 7;

	f_table->node_num[5] = 4;
	f_table->face_type[5] = ELEM_QUAD1;
	f_table->table[5][0] = 2;
	f_table->table[5][1] = 3;
	f_table->table[5][2] = 7;
	f_table->table[5][3] = 6;

	return(NORMAL_RC);
}


static RC set_fem_face_table_hexa2(FEM_FACE_TABLE *f_table)
{
	f_table->face_num = 6;

	f_table->node_num[0] = 8;
	f_table->face_type[0] = ELEM_QUAD2;
	f_table->table[0][0] = 0;
	f_table->table[0][1] = 3;
	f_table->table[0][2] = 2;
	f_table->table[0][3] = 1;
	f_table->table[0][4] = 11;
	f_table->table[0][5] = 10;
	f_table->table[0][6] = 9;
	f_table->table[0][7] = 8;

	f_table->node_num[1] = 8;
	f_table->face_type[1] = ELEM_QUAD2;
	f_table->table[1][0] = 0;
	f_table->table[1][1] = 1;
	f_table->table[1][2] = 5;
	f_table->table[1][3] = 4;
	f_table->table[1][4] = 8;
	f_table->table[1][5] = 17;
	f_table->table[1][6] = 12;
	f_table->table[1][7] = 16;

	f_table->node_num[2] = 8;
	f_table->face_type[2] = ELEM_QUAD2;
	f_table->table[2][0] = 1;
	f_table->table[2][1] = 2;
	f_table->table[2][2] = 6;
	f_table->table[2][3] = 5;
	f_table->table[2][4] = 9;
	f_table->table[2][5] = 18;
	f_table->table[2][6] = 13;
	f_table->table[2][7] = 17;

	f_table->node_num[3] = 8;
	f_table->face_type[3] = ELEM_QUAD2;
	f_table->table[3][0] = 7;
	f_table->table[3][1] = 4;
	f_table->table[3][2] = 5;
	f_table->table[3][3] = 6;
	f_table->table[3][4] = 15;
	f_table->table[3][5] = 12;
	f_table->table[3][6] = 13;
	f_table->table[3][7] = 14;

	f_table->node_num[4] = 8;
	f_table->face_type[4] = ELEM_QUAD2;
	f_table->table[4][0] = 3;
	f_table->table[4][1] = 0;
	f_table->table[4][2] = 4;
	f_table->table[4][3] = 7;
	f_table->table[4][4] = 11;
	f_table->table[4][5] = 16;
	f_table->table[4][6] = 15;
	f_table->table[4][7] = 19;

	f_table->node_num[5] = 8;
	f_table->face_type[5] = ELEM_QUAD2;
	f_table->table[5][0] = 2;
	f_table->table[5][1] = 3;
	f_table->table[5][2] = 7;
	f_table->table[5][3] = 6;
	f_table->table[5][4] = 10;
	f_table->table[5][5] = 19;
	f_table->table[5][6] = 14;
	f_table->table[5][7] = 18;

	return(NORMAL_RC);
}


void init_fem_edge_table(FEM_EDGE_TABLE *e_table)
{
	int ii1, ii2;

	e_table->edge_num = 0;
	for(ii1=0; ii1<FEM_MAX_EDGE; ii1++){
		e_table->node_num[ii1] = 0;
		e_table->edge_type[ii1] = ELEM_LINE1;
		for(ii2=0; ii2<FEM_MAX_EDGE_NODE; ii2++){
			e_table->table[ii1][ii2] = 0;
		}
	}
}


RC make_fem_edge_table(FEM_ELEM_TYPE type, FEM_EDGE_TABLE *e_table)
{
	init_fem_edge_table(e_table);

	switch(type){
	case ELEM_TRI1:
		RC_TRY( set_fem_edge_table_tri1(e_table) );
		break;
	case ELEM_TRI2:
		RC_TRY( set_fem_edge_table_tri2(e_table) );
		break;
	case ELEM_QUAD1:
		RC_TRY( set_fem_edge_table_quad1(e_table) );
		break;
	case ELEM_QUAD2:
		RC_TRY( set_fem_edge_table_quad2(e_table) );
		break;
	case ELEM_TETRA1:
		RC_TRY( set_fem_edge_table_tetra1(e_table) );
		break;
	case ELEM_TETRA2:
		RC_TRY( set_fem_edge_table_tetra2(e_table) );
		break;
	case ELEM_PENTA1:
		RC_TRY( set_fem_edge_table_penta1(e_table) );
		break;
	case ELEM_PENTA2:
		RC_TRY( set_fem_edge_table_penta2(e_table) );
		break;
	case ELEM_HEXA1:
		RC_TRY( set_fem_edge_table_hexa1(e_table) );
		break;
	case ELEM_HEXA2:
		RC_TRY( set_fem_edge_table_hexa2(e_table) );
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC set_fem_edge_table_tri1(FEM_EDGE_TABLE *e_table)
{
	e_table->edge_num = 3;

	e_table->node_num[0] = 2;
	e_table->edge_type[0] = ELEM_LINE1;
	e_table->table[0][0] = 1;
	e_table->table[0][1] = 2;

	e_table->node_num[1] = 2;
	e_table->edge_type[1] = ELEM_LINE1;
	e_table->table[1][0] = 2;
	e_table->table[1][1] = 0;

	e_table->node_num[2] = 2;
	e_table->edge_type[2] = ELEM_LINE1;
	e_table->table[2][0] = 0;
	e_table->table[2][1] = 1;

	return(NORMAL_RC);
}


static RC set_fem_edge_table_tri2(FEM_EDGE_TABLE *e_table)
{
	e_table->edge_num = 3;

	e_table->node_num[0] = 3;
	e_table->edge_type[0] = ELEM_LINE2;
	e_table->table[0][0] = 1;
	e_table->table[0][1] = 2;
	e_table->table[0][2] = 3;

	e_table->node_num[1] = 3;
	e_table->edge_type[1] = ELEM_LINE2;
	e_table->table[1][0] = 2;
	e_table->table[1][1] = 0;
	e_table->table[1][2] = 4;

	e_table->node_num[2] = 3;
	e_table->edge_type[2] = ELEM_LINE2;
	e_table->table[2][0] = 0;
	e_table->table[2][1] = 1;
	e_table->table[2][2] = 5;

	return(NORMAL_RC);
}


static RC set_fem_edge_table_quad1(FEM_EDGE_TABLE *e_table)
{
	e_table->edge_num = 4;

	e_table->node_num[0] = 2;
	e_table->edge_type[0] = ELEM_LINE1;
	e_table->table[0][0] = 0;
	e_table->table[0][1] = 1;

	e_table->node_num[1] = 2;
	e_table->edge_type[1] = ELEM_LINE1;
	e_table->table[1][0] = 1;
	e_table->table[1][1] = 2;

	e_table->node_num[2] = 2;
	e_table->edge_type[2] = ELEM_LINE1;
	e_table->table[2][0] = 2;
	e_table->table[2][1] = 3;

	e_table->node_num[3] = 2;
	e_table->edge_type[3] = ELEM_LINE1;
	e_table->table[3][0] = 3;
	e_table->table[3][1] = 0;

	return(NORMAL_RC);
}


static RC set_fem_edge_table_quad2(FEM_EDGE_TABLE *e_table)
{
	e_table->edge_num = 4;

	e_table->node_num[0] = 3;
	e_table->edge_type[0] = ELEM_LINE2;
	e_table->table[0][0] = 0;
	e_table->table[0][1] = 1;
	e_table->table[0][2] = 4;

	e_table->node_num[1] = 3;
	e_table->edge_type[1] = ELEM_LINE2;
	e_table->table[1][0] = 1;
	e_table->table[1][1] = 2;
	e_table->table[1][2] = 5;

	e_table->node_num[2] = 3;
	e_table->edge_type[2] = ELEM_LINE2;
	e_table->table[2][0] = 2;
	e_table->table[2][1] = 3;
	e_table->table[2][2] = 6;

	e_table->node_num[3] = 3;
	e_table->edge_type[3] = ELEM_LINE2;
	e_table->table[3][0] = 3;
	e_table->table[3][1] = 0;
	e_table->table[3][2] = 7;

	return(NORMAL_RC);
}


static RC set_fem_edge_table_tetra1(FEM_EDGE_TABLE *e_table)
{
	e_table->edge_num = 6;

	e_table->node_num[0] = 2;
	e_table->edge_type[0] = ELEM_LINE1;
	e_table->table[0][0] = 1;
	e_table->table[0][1] = 2;

	e_table->node_num[1] = 2;
	e_table->edge_type[1] = ELEM_LINE1;
	e_table->table[1][0] = 2;
	e_table->table[1][1] = 0;

	e_table->node_num[2] = 2;
	e_table->edge_type[2] = ELEM_LINE1;
	e_table->table[2][0] = 0;
	e_table->table[2][1] = 1;

	e_table->node_num[3] = 2;
	e_table->edge_type[3] = ELEM_LINE1;
	e_table->table[3][0] = 3;
	e_table->table[3][1] = 0;

	e_table->node_num[4] = 2;
	e_table->edge_type[4] = ELEM_LINE1;
	e_table->table[4][0] = 1;
	e_table->table[4][1] = 3;

	e_table->node_num[5] = 2;
	e_table->edge_type[5] = ELEM_LINE1;
	e_table->table[5][0] = 2;
	e_table->table[5][1] = 3;

	return(NORMAL_RC);
}


static RC set_fem_edge_table_tetra2(FEM_EDGE_TABLE *e_table)
{
	e_table->edge_num = 6;

	e_table->node_num[0] = 3;
	e_table->edge_type[0] = ELEM_LINE2;
	e_table->table[0][0] = 1;
	e_table->table[0][1] = 2;
	e_table->table[0][2] = 4;

	e_table->node_num[1] = 3;
	e_table->edge_type[1] = ELEM_LINE2;
	e_table->table[1][0] = 2;
	e_table->table[1][1] = 0;
	e_table->table[1][2] = 5;

	e_table->node_num[2] = 3;
	e_table->edge_type[2] = ELEM_LINE2;
	e_table->table[2][0] = 0;
	e_table->table[2][1] = 1;
	e_table->table[2][2] = 6;

	e_table->node_num[3] = 3;
	e_table->edge_type[3] = ELEM_LINE2;
	e_table->table[3][0] = 3;
	e_table->table[3][1] = 0;
	e_table->table[3][2] = 7;

	e_table->node_num[4] = 3;
	e_table->edge_type[4] = ELEM_LINE2;
	e_table->table[4][0] = 1;
	e_table->table[4][1] = 3;
	e_table->table[4][2] = 8;

	e_table->node_num[5] = 3;
	e_table->edge_type[5] = ELEM_LINE2;
	e_table->table[5][0] = 2;
	e_table->table[5][1] = 3;
	e_table->table[5][2] = 9;

	return(NORMAL_RC);
}


static RC set_fem_edge_table_penta1(FEM_EDGE_TABLE *e_table)
{
	e_table->edge_num = 9;

	e_table->node_num[0] = 2;
	e_table->edge_type[0] = ELEM_LINE1;
	e_table->table[0][0] = 1;
	e_table->table[0][1] = 2;

	e_table->node_num[1] = 2;
	e_table->edge_type[1] = ELEM_LINE1;
	e_table->table[1][0] = 2;
	e_table->table[1][1] = 0;

	e_table->node_num[2] = 2;
	e_table->edge_type[2] = ELEM_LINE1;
	e_table->table[2][0] = 0;
	e_table->table[2][1] = 1;

	e_table->node_num[3] = 2;
	e_table->edge_type[3] = ELEM_LINE1;
	e_table->table[3][0] = 4;
	e_table->table[3][1] = 5;

	e_table->node_num[4] = 2;
	e_table->edge_type[4] = ELEM_LINE1;
	e_table->table[4][0] = 5;
	e_table->table[4][1] = 3;

	e_table->node_num[5] = 2;
	e_table->edge_type[5] = ELEM_LINE1;
	e_table->table[5][0] = 3;
	e_table->table[5][1] = 4;

	e_table->node_num[6] = 2;
	e_table->edge_type[6] = ELEM_LINE1;
	e_table->table[6][0] = 0;
	e_table->table[6][1] = 3;

	e_table->node_num[7] = 2;
	e_table->edge_type[7] = ELEM_LINE1;
	e_table->table[7][0] = 1;
	e_table->table[7][1] = 4;

	e_table->node_num[8] = 2;
	e_table->edge_type[8] = ELEM_LINE1;
	e_table->table[8][0] = 2;
	e_table->table[8][1] = 5;

	return(NORMAL_RC);
}


static RC set_fem_edge_table_penta2(FEM_EDGE_TABLE *e_table)
{
	e_table->edge_num = 9;

	e_table->node_num[0] = 3;
	e_table->edge_type[0] = ELEM_LINE2;
	e_table->table[0][0] = 1;
	e_table->table[0][1] = 2;
	e_table->table[0][2] = 6;

	e_table->node_num[1] = 3;
	e_table->edge_type[1] = ELEM_LINE2;
	e_table->table[1][0] = 2;
	e_table->table[1][1] = 0;
	e_table->table[1][2] = 7;

	e_table->node_num[2] = 3;
	e_table->edge_type[2] = ELEM_LINE2;
	e_table->table[2][0] = 0;
	e_table->table[2][1] = 1;
	e_table->table[2][2] = 8;

	e_table->node_num[3] = 3;
	e_table->edge_type[3] = ELEM_LINE2;
	e_table->table[3][0] = 4;
	e_table->table[3][1] = 5;
	e_table->table[3][2] = 9;

	e_table->node_num[4] = 3;
	e_table->edge_type[4] = ELEM_LINE2;
	e_table->table[4][0] = 5;
	e_table->table[4][1] = 3;
	e_table->table[4][2] = 10;

	e_table->node_num[5] = 3;
	e_table->edge_type[5] = ELEM_LINE2;
	e_table->table[5][0] = 3;
	e_table->table[5][1] = 4;
	e_table->table[5][2] = 11;

	e_table->node_num[6] = 3;
	e_table->edge_type[6] = ELEM_LINE2;
	e_table->table[6][0] = 0;
	e_table->table[6][1] = 3;
	e_table->table[6][2] = 12;

	e_table->node_num[7] = 3;
	e_table->edge_type[7] = ELEM_LINE2;
	e_table->table[7][0] = 1;
	e_table->table[7][1] = 4;
	e_table->table[7][2] = 13;

	e_table->node_num[8] = 3;
	e_table->edge_type[8] = ELEM_LINE2;
	e_table->table[8][0] = 2;
	e_table->table[8][1] = 5;
	e_table->table[8][2] = 14;

	return(NORMAL_RC);
}


static RC set_fem_edge_table_hexa1(FEM_EDGE_TABLE *e_table)
{
	e_table->edge_num = 12;

	e_table->node_num[0] = 2;
	e_table->edge_type[0] = ELEM_LINE1;
	e_table->table[0][0] = 0;
	e_table->table[0][1] = 1;

	e_table->node_num[1] = 2;
	e_table->edge_type[1] = ELEM_LINE1;
	e_table->table[1][0] = 1;
	e_table->table[1][1] = 2;

	e_table->node_num[2] = 2;
	e_table->edge_type[2] = ELEM_LINE1;
	e_table->table[2][0] = 2;
	e_table->table[2][1] = 3;

	e_table->node_num[3] = 2;
	e_table->edge_type[3] = ELEM_LINE1;
	e_table->table[3][0] = 3;
	e_table->table[3][1] = 0;

	e_table->node_num[4] = 2;
	e_table->edge_type[4] = ELEM_LINE1;
	e_table->table[4][0] = 4;
	e_table->table[4][1] = 5;

	e_table->node_num[5] = 2;
	e_table->edge_type[5] = ELEM_LINE1;
	e_table->table[5][0] = 5;
	e_table->table[5][1] = 6;

	e_table->node_num[6] = 2;
	e_table->edge_type[6] = ELEM_LINE1;
	e_table->table[6][0] = 6;
	e_table->table[6][1] = 7;

	e_table->node_num[7] = 2;
	e_table->edge_type[7] = ELEM_LINE1;
	e_table->table[7][0] = 7;
	e_table->table[7][1] = 4;

	e_table->node_num[8] = 2;
	e_table->edge_type[8] = ELEM_LINE1;
	e_table->table[8][0] = 0;
	e_table->table[8][1] = 4;

	e_table->node_num[9] = 2;
	e_table->edge_type[9] = ELEM_LINE1;
	e_table->table[9][0] = 1;
	e_table->table[9][1] = 5;

	e_table->node_num[10] = 2;
	e_table->edge_type[10] = ELEM_LINE1;
	e_table->table[10][0] = 2;
	e_table->table[10][1] = 6;

	e_table->node_num[11] = 2;
	e_table->edge_type[11] = ELEM_LINE1;
	e_table->table[11][0] = 3;
	e_table->table[11][1] = 7;

	return(NORMAL_RC);
}


static RC set_fem_edge_table_hexa2(FEM_EDGE_TABLE *e_table)
{
	e_table->edge_num = 12;

	e_table->node_num[0] = 3;
	e_table->edge_type[0] = ELEM_LINE2;
	e_table->table[0][0] = 0;
	e_table->table[0][1] = 1;
	e_table->table[0][2] = 8;

	e_table->node_num[1] = 3;
	e_table->edge_type[1] = ELEM_LINE2;
	e_table->table[1][0] = 1;
	e_table->table[1][1] = 2;
	e_table->table[1][2] = 9;

	e_table->node_num[2] = 3;
	e_table->edge_type[2] = ELEM_LINE2;
	e_table->table[2][0] = 2;
	e_table->table[2][1] = 3;
	e_table->table[2][2] = 10;

	e_table->node_num[3] = 3;
	e_table->edge_type[3] = ELEM_LINE2;
	e_table->table[3][0] = 3;
	e_table->table[3][1] = 0;
	e_table->table[3][2] = 11;

	e_table->node_num[4] = 3;
	e_table->edge_type[4] = ELEM_LINE2;
	e_table->table[4][0] = 4;
	e_table->table[4][1] = 5;
	e_table->table[4][2] = 12;

	e_table->node_num[5] = 3;
	e_table->edge_type[5] = ELEM_LINE2;
	e_table->table[5][0] = 5;
	e_table->table[5][1] = 6;
	e_table->table[5][2] = 13;

	e_table->node_num[6] = 3;
	e_table->edge_type[6] = ELEM_LINE2;
	e_table->table[6][0] = 6;
	e_table->table[6][1] = 7;
	e_table->table[6][2] = 14;

	e_table->node_num[7] = 3;
	e_table->edge_type[7] = ELEM_LINE2;
	e_table->table[7][0] = 7;
	e_table->table[7][1] = 4;
	e_table->table[7][2] = 15;

	e_table->node_num[8] = 3;
	e_table->edge_type[8] = ELEM_LINE2;
	e_table->table[8][0] = 0;
	e_table->table[8][1] = 4;
	e_table->table[8][2] = 16;

	e_table->node_num[9] = 3;
	e_table->edge_type[9] = ELEM_LINE2;
	e_table->table[9][0] = 1;
	e_table->table[9][1] = 5;
	e_table->table[9][2] = 17;

	e_table->node_num[10] = 3;
	e_table->edge_type[10] = ELEM_LINE2;
	e_table->table[10][0] = 2;
	e_table->table[10][1] = 6;
	e_table->table[10][2] = 18;

	e_table->node_num[11] = 3;
	e_table->edge_type[11] = ELEM_LINE2;
	e_table->table[11][0] = 3;
	e_table->table[11][1] = 7;
	e_table->table[11][2] = 19;

	return(NORMAL_RC);
}
RC free_fem_elem_matrix(FEM_ELEM_MATRIX *elem_matrix)
{
	if((elem_matrix->element >= 0)&&(elem_matrix->size > 0)){
		if(elem_matrix->index == NULL) return(ARG_ERROR_RC);

		free(elem_matrix->index);
		RC_TRY( free2D(elem_matrix->size,
		               elem_matrix->size, elem_matrix->matrix) );
	}
	init_fem_elem_matrix(elem_matrix);

	return(NORMAL_RC);
}


void init_fem_elem_lookup(FEM_ELEM_LOOKUP *elem_lookup)
{
	elem_lookup->node = -1;
	elem_lookup->size = 0;
	elem_lookup->alloc_size = 0;
	elem_lookup->element = NULL;
	elem_lookup->r_size = 0;
	elem_lookup->r_alloc_size = 0;
	elem_lookup->r_element = NULL;
}


void init_fem_node_matrix(FEM_NODE_MATRIX *nmat)
{
	nmat->node = -1;
	nmat->shift = 0;
	nmat->size = 0;
	nmat->alloc_size = 0;
	nmat->col_node = NULL;
	nmat->col_shift = NULL;
	nmat->matrix = NULL;
}



RC free_fem_elem_lookup(FEM_ELEM_LOOKUP *elem_lookup)
{
	if((elem_lookup->node >= 0)&&(elem_lookup->alloc_size > 0)){
		if(elem_lookup->element == NULL) return(ARG_ERROR_RC);

		free(elem_lookup->element);
	}
	if((elem_lookup->node >= 0)&&(elem_lookup->r_alloc_size > 0)){
		if(elem_lookup->r_element == NULL) return(ARG_ERROR_RC);

		free(elem_lookup->r_element);
	}
	init_fem_elem_lookup(elem_lookup);

	return(NORMAL_RC);
}


RC free_fem_node_matrix(FEM_NODE_MATRIX *nmat)
{
	if(nmat->alloc_size > 0){
		if(nmat->col_node == NULL) return(ARG_ERROR_RC);
		if(nmat->col_shift == NULL) return(ARG_ERROR_RC);
		if(nmat->matrix == NULL) return(ARG_ERROR_RC);
		free(nmat->col_node);
		free(nmat->col_shift);
		free(nmat->matrix);
	}
	init_fem_node_matrix(nmat);

	return(NORMAL_RC);
}


FEM_STRESS_ARRAY fem_strain2stress_array(FEM_STRAIN_ARRAY strain)
{
	FEM_STRESS_ARRAY ret;

	ret.size = strain.size;
	ret.alloc_size = strain.alloc_size;
	ret.source = strain.source;
	ret.array = strain.array;
	ret.sort_flag = strain.sort_flag;

	return(ret);
}


FEM_STRAIN_ARRAY fem_stress2strain_array(FEM_STRESS_ARRAY stress)
{
	FEM_STRAIN_ARRAY ret;

	ret.size = stress.size;
	ret.alloc_size = stress.alloc_size;
	ret.source = stress.source;
	ret.array = stress.array;
	ret.sort_flag = stress.sort_flag;

	return(ret);
}


/* generate suitable stiffness tensor and stress-strain matrix */
RC fem_fill_material(FEM_MATERIAL_PROP *prop)
{
	RC rc;

	switch(prop->mat_type){
	case MAT_ISOTROPIC:
		prop->E_vect.x = prop->E;
		prop->E_vect.y = prop->E;
		prop->E_vect.z = prop->E;
		prop->nu_vect.yz = prop->nu;
		prop->nu_vect.zx = prop->nu;
		prop->nu_vect.xy = prop->nu;
		prop->alpha_vect.x = prop->alpha;
		prop->alpha_vect.y = prop->alpha;
		prop->alpha_vect.z = prop->alpha;
		prop->alpha_vect.yz = 0.0;
		prop->alpha_vect.zx = 0.0;
		prop->alpha_vect.xy = 0.0;

		if((rc = put_iso_D_mat(prop)) != NORMAL_RC){
			fprintf(stderr,"put_iso_D_mat < ");
			return(rc);
		}
		break;
	case MAT_ORTHOTROPIC:
		if((rc = put_ortho_D_mat(prop)) != NORMAL_RC){
			fprintf(stderr,"put_ortho_D_mat < ");
			return(rc);
		}
		break;
	case MAT_ANISOTROPIC:
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
	case ANA_PLANE_BEND:
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
	case ANA_PLANE_SHEAR:
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
static RC put_iso_D_mat(FEM_MATERIAL_PROP *prop)
{
	double fact;
	double div_tmp;

	clean_D_matrix(prop);
	switch(prop->ana_type){
	case ANA_PLANE_STRESS:
		div_tmp = 1.0 - ((prop->nu) * (prop->nu));
		if(nearly_eq(div_tmp, 0.0)){
			fprintf(stderr,"nearly_eq:1 < ");
			return(CAL_ERROR_RC);
		}
		fact = (prop->E)/div_tmp;
		prop->D_matrix[0][0] = prop->D_matrix[1][1] = fact;
		prop->D_matrix[0][1] = prop->D_matrix[1][0] = fact * (prop->nu);
		prop->D_matrix[2][2] = 0.5 * fact * (1.0 - prop->nu);
		break;
	case ANA_PLANE_STRAIN:
		div_tmp = (1.0 + prop->nu) * (1.0 - 2.0 * prop->nu);
		if(nearly_eq(div_tmp, 0.0)){
			fprintf(stderr,"nearly_eq:2 < ");
			return(CAL_ERROR_RC);
		}
		fact = (prop->E) * (1.0 - (prop->nu)) / div_tmp;
		div_tmp = 1.0 - prop->nu;
		if(nearly_eq(div_tmp, 0.0)){
			fprintf(stderr,"nearly_eq:3 < ");
			return(CAL_ERROR_RC);
		}
		prop->D_matrix[0][0] = prop->D_matrix[1][1] = fact;
		prop->D_matrix[0][1] = prop->D_matrix[1][0] 
		                     = fact * (prop->nu) / div_tmp;
		prop->D_matrix[2][2] = 0.5 * fact * (1.0 - 2.0*prop->nu) / div_tmp;
		break;
	case ANA_AXISYMMETRY:
		/* z,r,theta,rz */
		div_tmp = (1.0 + prop->nu) * (1.0 - 2.0 * prop->nu);
		if(nearly_eq(div_tmp, 0.0)){
			fprintf(stderr,"nearly_eq:4 < ");
			return(CAL_ERROR_RC);
		}
		fact = (prop->E) * (1.0 - (prop->nu)) / div_tmp;
		div_tmp = 1.0 - prop->nu;
		if(nearly_eq(div_tmp, 0.0)){
			fprintf(stderr,"nearly_eq:5 < ");
			return(CAL_ERROR_RC);
		}
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
		if(nearly_eq(div_tmp, 0.0)){
			fprintf(stderr,"nearly_eq:6 < ");
			return(CAL_ERROR_RC);
		}
		fact = (prop->E) * (1.0 - (prop->nu)) / div_tmp;
		div_tmp = 1.0 - prop->nu;
		if(nearly_eq(div_tmp, 0.0)){
			fprintf(stderr,"nearly_eq:7 < ");
			return(CAL_ERROR_RC);
		}
		prop->D_matrix[0][0] = prop->D_matrix[1][1]
		                     = prop->D_matrix[2][2] = fact;
		prop->D_matrix[0][1] = prop->D_matrix[1][0] 
		                     = prop->D_matrix[0][2] = prop->D_matrix[2][0] 
		                     = prop->D_matrix[1][2] = prop->D_matrix[2][1] 
							 = fact * (prop->nu) / div_tmp;
		prop->D_matrix[3][3] = prop->D_matrix[4][4] = prop->D_matrix[5][5]
		                     = 0.5 * fact * (1.0 - 2.0*prop->nu) / div_tmp;
		break;
	case ANA_PLANE_BEND:
		div_tmp = 12.0 * ( 1.0 - ((prop->nu) * (prop->nu)) );
		if(nearly_eq(div_tmp, 0.0)){
			fprintf(stderr,"nearly_eq:8 < ");
			return(CAL_ERROR_RC);
		}
		fact = (prop->E) / div_tmp;
		prop->D_matrix[0][0] = prop->D_matrix[1][1] = fact;
		prop->D_matrix[0][1] = prop->D_matrix[1][0] = fact * (prop->nu);
		prop->D_matrix[2][2] = 0.5 * fact * (1.0 - prop->nu);
		break;
	case ANA_PLANE_SHEAR:
		div_tmp = 2.0 * SHEAR_FACT * ( 1.0 + (prop->nu) );
		if(nearly_eq(div_tmp, 0.0)){
			fprintf(stderr,"nearly_eq:9 < ");
			return(CAL_ERROR_RC);
		}
		fact = (prop->E) / div_tmp;
		prop->D_matrix[0][0] = prop->D_matrix[1][1] = fact;
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


/* set D_matrix for orthotropic material */
static RC put_ortho_D_mat(FEM_MATERIAL_PROP *prop)
{
	double div_tmp;
	double nu_yx;
	double n, m;
	double nu1;
	double fact_tmp;

	clean_D_matrix(prop);
	switch(prop->ana_type){
	case ANA_PLANE_STRESS:
		if(nearly_eq(prop->E_vect.x, 0.0)){
			fprintf(stderr,"nearly_eq:1 < ");
			return(CAL_ERROR_RC);
		}
		nu_yx = (prop->E_vect.y)*(prop->nu_vect.xy)/(prop->E_vect.x);
		div_tmp = 1.0 - nu_yx * (prop->nu_vect.xy);
		if(nearly_eq(div_tmp, 0.0)){
			fprintf(stderr,"nearly_eq:2 < ");
			return(CAL_ERROR_RC);
		}
		prop->D_matrix[0][0] = (prop->E_vect.x) / div_tmp;
		prop->D_matrix[1][1] = (prop->E_vect.y) / div_tmp;
		prop->D_matrix[0][1] = prop->D_matrix[1][0]
		                     = (prop->E_vect.x) * nu_yx / div_tmp;
		prop->D_matrix[2][2] = prop->G_vect.xy;
		break;
	case ANA_PLANE_STRAIN:
		if(nearly_eq(prop->E_vect.y, 0.0)){
			fprintf(stderr,"nearly_eq:3 < ");
			return(CAL_ERROR_RC);
		}
		n = (prop->E_vect.x)/(prop->E_vect.y);
		m = (prop->G_vect.xy)/(prop->E_vect.y);
		if(nearly_eq(prop->G_vect.zx, 0.0)){
			fprintf(stderr,"nearly_eq:4 < ");
			return(CAL_ERROR_RC);
		}
		nu1 = 0.5 * (prop->E_vect.x)/(prop->G_vect.zx) - 1.0;
		div_tmp = (1 + nu1)
		         *(1 - nu1
				 - 2.0 * n * (prop->nu_vect.xy) * (prop->nu_vect.xy));
		if(nearly_eq(div_tmp, 0.0)){
			fprintf(stderr,"nearly_eq:5 < ");
			return(CAL_ERROR_RC);
		}
		fact_tmp = (prop->E_vect.y)/div_tmp;
		prop->D_matrix[0][0] = n * (1.0 - n * (prop->nu_vect.xy)
		                                    * (prop->nu_vect.xy))
											* fact_tmp;
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


RC Gauss_const_max(FEM_ELEM_TYPE type, double **points,
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


RC Gauss_const_default(FEM_ELEM_TYPE type, double **points,
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


RC Gauss_const_1(FEM_ELEM_TYPE type, double **points, double *weight)
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


/* $BLdBj$N$R$:$_!"1~NO$N?t$rJV5Q(B */
int element_ssnum(FEM_ELEM_TYPE type)
{
	switch(type){
	case ELEM_TRI1:
	case ELEM_TRI2:
	case ELEM_QUAD1:
	case ELEM_QUAD2:
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


int element_dim(FEM_ELEM_TYPE type)
{
	switch(type){
	case ELEM_LINE1:
	case ELEM_LINE2:
		return(1);
	case ELEM_TRI1:
	case ELEM_TRI2:
	case ELEM_QUAD1:
	case ELEM_QUAD2:
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


int element_order(FEM_ELEM_TYPE type)
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
	case ELEM_TETRA2:
	case ELEM_PENTA2:
	case ELEM_HEXA2:
		return(2);
	default:
		break;
	}

	return(-1);
}


RC cal_det_J(FEM_ELEMENT element, const double local_coord[],
             double **node_location, double *det_J,
             double **ws_mat_3_N, double **ws_mat_3_3)
{
	double **dN = ws_mat_3_N;
	double **jacobi = ws_mat_3_3;
	int dim;

	RC_TRY( set_dN_matrix(element.type, local_coord, dN) );
	dim = element_dim(element.type);
	if( (dim < 2)||(dim > 3) ) return(ARG_ERROR_RC);
	mul_matrix(dim, element.node_num, dim, dN, node_location, jacobi);
	*det_J = determinant(dim, jacobi);

	return(NORMAL_RC);
}


RC cal_d_A(FEM_ELEMENT surf, const double local_coord[],
           double **node_location, TRANSLATION3D *d_A,
           double **ws_mat_3_N, double **ws_mat_3_3)
{
	double **dN = ws_mat_3_N;
	double **jacobi = ws_mat_3_3;
	int dim;
	TRANSLATION3D d_xi, d_eta;

	init_matrix(3, FEM_MAX_NODE, dN);  /* $B=i4|2=$,I,MW(B */
	RC_TRY( set_dN_matrix(surf.type, local_coord, dN) );
	mul_matrix(3, surf.node_num, 3, dN, node_location, jacobi);
	dim = element_dim(surf.type);
	if(dim == 2){
		d_xi.x = jacobi[0][0];
		d_xi.y = jacobi[0][1];
		d_xi.z = jacobi[0][2];
		d_eta.x = jacobi[1][0];
		d_eta.y = jacobi[1][1];
		d_eta.z = jacobi[1][2];
		*d_A = outer_product3d(d_xi, d_eta);
	}else if(dim == 1){
		d_A->x = jacobi[0][1];
		d_A->y = -jacobi[0][0];
		d_A->z = 0.0;
	}else{
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


RC interpolate_disp(FEM_ELEMENT element, const double local_coord[],
                    double **disp_matrix, TRANS_ROTATION3D *disp)
{
	int ii1;
	double N[FEM_MAX_NODE];

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


RC realloc_fem_node_matrix(FEM_NODE_MATRIX *nmat)
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

	nmat->col_node = realloc(nmat->col_node, (nmat->alloc_size)*sizeof(int));
	RC_NULL_CHK( nmat->col_node );

	nmat->col_shift = realloc(nmat->col_shift, (nmat->alloc_size)*sizeof(int));
	RC_NULL_CHK( nmat->col_shift );

	nmat->matrix = realloc(nmat->matrix,
	                       (nmat->alloc_size)*sizeof(double[3][3]));
	RC_NULL_CHK( nmat->matrix );

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

