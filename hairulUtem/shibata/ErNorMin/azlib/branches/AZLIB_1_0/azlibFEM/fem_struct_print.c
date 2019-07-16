/*********************************************************************
 * fem_struct_print.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Thanks to following contributors.
 *  <Tomoyuki TSUBATA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: fem_struct_print.c,v 1.15 2003/07/22 11:44:01 adachi Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "rc.h"
#include "math_utl.h"
#include "fem_struct.h"

static void print_i_info(FILE *fp, int size, const int i_info[]);
static void print_d_info(FILE *fp, int size, const double d_info[]);


static void print_i_info(FILE *fp, int size, const int *i_info)
{
	int ii1;

	fprintf(fp, "i_info: ");
	for(ii1=0; ii1<size; ii1++){
		fprintf(fp, "%10d ", i_info[ii1]);
	}
	fprintf(fp, "\n");
}


static void print_d_info(FILE *fp, int size, const double *d_info)
{
	int ii1;

	fprintf(fp, "d_info: ");
	for(ii1=0; ii1<size; ii1++){
		fprintf(fp, "%10.3e ", d_info[ii1]);
	}
	fprintf(fp, "\n");
}


void print_fem_node(FILE *fp, FEM_NODE node)
{
	fprintf(fp, "label: %5d\n", node.label);
	fprintf(fp, "p: ");
	print_translation3d(fp, node.p);
	fprintf(fp, "renum: %5d\n", node.renum);
	print_i_info(fp, FEM_NODE_INFO, node.i_info);
	print_d_info(fp, FEM_NODE_INFO, node.d_info);
}


void print_fem_elem_type(FILE *fp, FEM_ELEM_TYPE type)
{
	switch(type){
	case ELEM_RBAR:
		fprintf(fp, "ELEM_RBAR\n");
		break;
	case ELEM_RBODY:
		fprintf(fp, "ELEM_RBODY\n");
		break;
	case ELEM_RFORCE:
		fprintf(fp, "ELEM_RFORCE\n");
		break;
	case ELEM_MPC:
		fprintf(fp, "ELEM_MPC\n");
		break;
	case ELEM_SPRING:
		fprintf(fp, "ELEM_SPRING\n");
		break;
	case ELEM_MASS:
		fprintf(fp, "ELEM_MASS\n");
		break;
	case ELEM_BEAM1:
		fprintf(fp, "ELEM_BEAM1\n");
		break;
	case ELEM_BEAM2:
		fprintf(fp, "ELEM_BEAM2\n");
		break;
	case ELEM_SHELL1:
		fprintf(fp, "ELEM_SHELL1\n");
		break;
	case ELEM_SHELL2:
		fprintf(fp, "ELEM_SHELL2\n");
		break;
	case ELEM_LINE1:
		fprintf(fp, "ELEM_LINE1\n");
		break;
	case ELEM_LINE2:
		fprintf(fp, "ELEM_LINE2\n");
		break;
	case ELEM_TRI1:
		fprintf(fp, "ELEM_TRI1\n");
		break;
	case ELEM_TRI2:
		fprintf(fp, "ELEM_TRI2\n");
		break;
	case ELEM_QUAD1:
		fprintf(fp, "ELEM_QUAD1\n");
		break;
	case ELEM_QUAD2:
		fprintf(fp, "ELEM_QUAD2\n");
		break;
	case ELEM_TETRA1:
		fprintf(fp, "ELEM_TETRA1\n");
		break;
	case ELEM_TETRA2:
		fprintf(fp, "ELEM_TETRA2\n");
		break;
	case ELEM_PENTA1:
		fprintf(fp, "ELEM_PENTA1\n");
		break;
	case ELEM_PENTA2:
		fprintf(fp, "ELEM_PENTA2\n");
		break;
	case ELEM_HEXA1:
		fprintf(fp, "ELEM_HEXA1\n");
		break;
	case ELEM_HEXA2:
		fprintf(fp, "ELEM_HEXA2\n");
		break;
	default:
		fprintf(fp, "unknown\n");
	}
}


void print_fem_element(FILE *fp, FEM_ELEMENT elem)
{
	int ii1;

	fprintf(fp, "label: %5d\n", elem.label);
	fprintf(fp, "type: ");
	print_fem_elem_type(fp, elem.type);
	fprintf(fp, "material: %5d\n", elem.material);
	fprintf(fp, "physical: %5d\n", elem.physical);
	fprintf(fp, "node_num: %5d\n", elem.node_num);
	fprintf(fp, "node:       ");
	for(ii1=0;ii1<FEM_MAX_NODE;ii1++){
		fprintf(fp, "%5d ", elem.node[ii1]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "volume: %15.7e\n", elem.volume);
	print_i_info(fp, FEM_ELEMENT_INFO, elem.i_info);
	print_d_info(fp, FEM_ELEMENT_INFO, elem.d_info);
}


void print_fem_rigid_element(FILE *fp, const FEM_RIGID_ELEMENT *elem)
{
	int ii1;

	fprintf(fp, "label: %5d\n", elem->label);
	fprintf(fp, "set_id: %5d\n", elem->set_id);
	fprintf(fp, "type: ");
	print_fem_elem_type(fp, elem->type);
	fprintf(fp, "physical: %5d\n", elem->physical);
	fprintf(fp, "node_num: %5d\n", elem->node_num);
	fprintf(fp, "node:");
	for(ii1=0;ii1<elem->node_num;ii1++){
		fprintf(fp, " %5d", elem->node[ii1]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "weight:");
	for(ii1=0;ii1<elem->weight_num;ii1++){
		fprintf(fp, " %10.2e", elem->weight[ii1]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "dof:");
	for(ii1=0;ii1<elem->dof_num;ii1++){
		fprintf(fp, " %5d", elem->dof[ii1]);
	}
	fprintf(fp, "\n");
	print_i_info(fp, FEM_RIGID_ELEMENT_INFO, elem->i_info);
	print_d_info(fp, FEM_RIGID_ELEMENT_INFO, elem->d_info);
}


void print_fem_mat_type(FILE *fp, FEM_MAT_TYPE type)
{
	switch(type){
	case MAT_ISOTROPIC:
		fprintf(fp, "MAT_ISOTROPIC\n");
		break;
	case MAT_ORTHOTROPIC:
		fprintf(fp, "MAT_ORTHOTROPIC\n");
		break;
	case MAT_ANISOTROPIC:
		fprintf(fp, "MAT_ANISOTROPIC\n");
		break;
	default:
		fprintf(fp, "unknown\n");
	}
}


void print_fem_analysis_type(FILE *fp, FEM_ANALYSIS_TYPE type)
{
	switch(type){
	case ANA_3D:
		fprintf(fp, "ANA_3D\n");
		break;
	case ANA_PLANE_STRESS:
		fprintf(fp, "ANA_PLANE_STRESS\n");
		break;
	case ANA_PLANE_STRAIN:
		fprintf(fp, "ANA_PLANE_STRAIN\n");
		break;
	case ANA_AXISYMMETRY:
		fprintf(fp, "NA_AXISYMMETRY\n");
		break;
	default:
		fprintf(fp, "unknown\n");
	}
}

void print_fem_material_prop(FILE *fp, FEM_MATERIAL_PROP m_prop)
{
	int ii1,ii2,ii3,ii4;

	fprintf(fp, "label: %5d\n", m_prop.label);
	fprintf(fp, "mat_type: ");
	print_fem_mat_type(fp, m_prop.mat_type);
	fprintf(fp, "ana_type: ");
	print_fem_analysis_type(fp, m_prop.ana_type);
	fprintf(fp, "E: %15.7e\n",m_prop.E);
	fprintf(fp, "G: %15.7e\n",m_prop.G);
	fprintf(fp, "nu: %15.7e\n",m_prop.nu);
	fprintf(fp, "alpha: %15.7e\n",m_prop.alpha);
	fprintf(fp, "K: %15.7e\n",m_prop.K);
	fprintf(fp, "E_vect: ");
	print_translation3d(fp, m_prop.E_vect);
	fprintf(fp, "G_vect: ");
	print_rotation3d(fp, m_prop.G_vect);
	fprintf(fp, "nu_vect: ");
	print_rotation3d(fp, m_prop.nu_vect);
	fprintf(fp, "alpha_vect: ");
	print_trans_rotation3d(fp, m_prop.alpha_vect);
	fprintf(fp, "D_matrix: ");
	for(ii1=0;ii1<3;ii1++){
		for(ii2=0;ii2<3;ii2++){
			fprintf(fp, "%10.3e ", m_prop.D_matrix[ii1][ii2]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
	fprintf(fp, "C: ");
	for(ii1=0;ii1<3;ii1++){
		for(ii3=0;ii3<3;ii3++){
			for(ii2=0;ii2<3;ii2++){
				for(ii4=0;ii4<3;ii4++){
					fprintf(fp, "%10.3e ", m_prop.C[ii1][ii2][ii3][ii4]);
				}
			}
			fprintf(fp, "\n");
		}
	}
	fprintf(fp, "\n");
	fprintf(fp, "rho: %15.7e\n",m_prop.rho);
	fprintf(fp, "tref: %15.7e\n",m_prop.tref);
	print_i_info(fp, FEM_MATERIAL_INFO, m_prop.i_info);
	print_d_info(fp, FEM_MATERIAL_INFO, m_prop.d_info);

}

void print_fem_section_type(FILE *fp, FEM_SECTION_TYPE type)
{
        switch(type){
	case SECTION_ROD:
		fprintf(fp, "SECTION_ROD\n");
		break;
	case SECTION_TUBE:
		fprintf(fp, "SECTION_TUBE\n");
		break;
	case SECTION_BAR:
		fprintf(fp, "SECTION_BAR\n");
		break;
	case SECTION_BOX:
		fprintf(fp, "SECTION_BOX\n");
		break;
	case SECTION_BOX1:
		fprintf(fp, "SECTION_BOX1\n");
		break;
	case SECTION_CHAN:
		fprintf(fp, "SECTION_CHAN\n");
		break;
	case SECTION_CHAN1:
		fprintf(fp, "SECTION_CHAN1\n");
		break;
	case SECTION_CHAN2:
		fprintf(fp, "SECTION_CHAN2\n");
		break;
	case SECTION_I:
		fprintf(fp, "SECTION_I\n");
		break;
	case SECTION_I1:
		fprintf(fp, "SECTION_I1\n");
		break;
	case SECTION_T:
		fprintf(fp, "SECTION_T\n");
		break;
	case SECTION_T1:
		fprintf(fp, "SECTION_T1\n");
		break;
	case SECTION_T2:
		fprintf(fp, "SECTION_T2\n");
		break;
	case SECTION_H:
		fprintf(fp, "SECTION_H\n");
		break;
	case SECTION_Z:
		fprintf(fp, "SECTION_Z\n");
		break;
	case SECTION_CROSS:
		fprintf(fp, "SECTION_CROSS\n");
		break;
	case SECTION_HEXA:
		fprintf(fp, "SECTION_HEXA\n");
		break;
	case SECTION_HAT:
		fprintf(fp, "SECTION_HAT\n");
		break;
	default:
		fprintf(fp, "unknown\n");
	}

}

void print_fem_physical_prop(FILE *fp, FEM_PHYSICAL_PROP p_prop)
{
	fprintf(fp, "label: %5d\n", p_prop.label);
	fprintf(fp, "section: ");
	print_fem_section_type(fp,p_prop.section);
	print_i_info(fp, FEM_PHYSICAL_INFO, p_prop.i_info);
	print_d_info(fp, FEM_PHYSICAL_INFO, p_prop.d_info);
}

void print_fem_bc_type(FILE *fp, FEM_BC_TYPE type)
{
	switch(type){
	case BC_FREE:
		fprintf(fp, "BC_FREE");
		break;
	case BC_FIX:
		fprintf(fp, "BC_FIX");
		break;
	case BC_FIX_NL:
		fprintf(fp, "BC_FIX_NL");
		break;
	case BC_MPC:
		fprintf(fp, "BC_MPC");
		break;
	case BC_FORCE:
		fprintf(fp, "BC_FORCE");
		break;
	case BC_FORCE_NL:
		fprintf(fp, "BC_FORCE_NL");
		break;
	case BC_VELO:
		fprintf(fp, "BC_VELO");
		break;
	case BC_PRESS:
		fprintf(fp, "BC_PRESS");
		break;
	case BC_PRESS_D1:
		fprintf(fp, "BC_PRESS_D1");
		break;
	case BC_PRESS_D2:
		fprintf(fp, "BC_PRESS_D2");
		break;
	case BC_TRACTION:
		fprintf(fp, "BC_TRACTION");
		break;
	case BC_TRACTION_D1:
		fprintf(fp, "BC_TRACTION_D1");
		break;
	case BC_TRACTION_D2:
		fprintf(fp, "BC_TRACTION_D2");
		break;
	case BC_TEMP:
		fprintf(fp, "BC_TEMP");
		break;
	case BC_TEMP_NL:
		fprintf(fp, "BC_TEMP_NL");
		break;
	case BC_K_DIAGONAL:
		fprintf(fp, "BC_K_DIAGONAL");
		break;
	default:
		fprintf(fp, "unknown");
	}
}

void print_fem_bc_type_3d(FILE *fp, FEM_BC_TYPE3D type)
{
	fprintf(fp, "  x:");
	print_fem_bc_type(fp, type.x);
	fprintf(fp, "\n  y:");
	print_fem_bc_type(fp, type.y);
	fprintf(fp, "\n  z:");
	print_fem_bc_type(fp, type.z);
	fprintf(fp, "\n  xy:");
	print_fem_bc_type(fp, type.xy);
	fprintf(fp, "\n  yz:");
	print_fem_bc_type(fp, type.yz);
	fprintf(fp, "\n  zx:");
	print_fem_bc_type(fp, type.zx);
	fprintf(fp, "\n");
}

void print_fem_bc(FILE *fp, FEM_BC bc)
{
	int ii1;

	fprintf(fp, "node:   %d\n", bc.node);
	fprintf(fp, "set_id: %d\n", bc.set_id);
	fprintf(fp, "v_type:\n");
	print_fem_bc_type_3d(fp, bc.v_type);
	fprintf(fp, "v:\n");
	print_trans_rotation3d(fp, bc.v);
	fprintf(fp, "v_nl:\n");
	print_i_trans_rotation3d(fp, bc.v_nl);
	fprintf(fp, "s_type: ");
	for(ii1=0; ii1<FEM_BC_SCALAR; ii1++){
		print_fem_bc_type(fp, bc.s_type[ii1]);
		fprintf(fp, "  ");
	}
	fprintf(fp, "\n");
	fprintf(fp, "s: ");
	for(ii1=0; ii1<FEM_BC_SCALAR; ii1++){
		fprintf(fp, " %10.3e ", bc.s[ii1]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "s_nl: ");
	for(ii1=0; ii1<FEM_BC_SCALAR; ii1++){
		fprintf(fp, " %5d ", bc.s_nl[ii1]);
	}
	fprintf(fp, "\n");
	print_i_info(fp, FEM_BC_INFO, bc.i_info);
	print_d_info(fp, FEM_BC_INFO, bc.d_info);
	fprintf(fp, "\n");
}

void print_fem_elem_bc(FILE *fp, FEM_ELEM_BC bc)
{
	int ii1;

	fprintf(fp, "element: %5d\n", bc.element);
	fprintf(fp, "set_id: %5d\n", bc.set_id);
	fprintf(fp, "type: ");
	print_fem_bc_type(fp, bc.type);
	fprintf(fp, "face_id: %5d\n", bc.face_id);
	fprintf(fp, "val.s %15.7e\n", bc.val.s);
	fprintf(fp, "val.v: ");
	print_trans_rotation3d(fp, bc.val.v);
	fprintf(fp, "val.ns: ");
	for(ii1=0; ii1<FEM_MAX_NODE; ii1++){
		fprintf(fp, "%10.2e ", bc.val.ns[ii1]);
	}
	fprintf(fp, "\n");
	print_translation3d(fp, bc.dir);
	print_i_info(fp, FEM_ELEM_BC_INFO, bc.i_info);
	print_d_info(fp, FEM_ELEM_BC_INFO, bc.d_info);
}

void print_fem_disp_velo(FILE *fp, FEM_DISP_VELO d_velo)
{
	fprintf(fp, "node: %5d\n", d_velo.node);
	fprintf(fp, "v: ");
	print_trans_rotation3d(fp, d_velo.v);
	fprintf(fp, "s: %15.7e\n", d_velo.s);
	print_i_info(fp, FEM_DISP_VELO_INFO, d_velo.i_info);
	print_d_info(fp, FEM_DISP_VELO_INFO, d_velo.d_info);
}

void print_fem_stress_strain(FILE *fp, FEM_STRESS_STRAIN ss)
{
	fprintf(fp, "element: %5d\n", ss.element);
	fprintf(fp, "node: %5d\n", ss.node);
	fprintf(fp, "v: ");
	print_trans_rotation3d(fp, ss.v);
	print_i_info(fp, FEM_STRESS_STRAIN_INFO, ss.i_info);
	print_d_info(fp, FEM_STRESS_STRAIN_INFO, ss.d_info);
}


void print_fem_sensitivity(FILE *fp, FEM_SENSITIVITY sensi)
{
	fprintf(fp, "element: %5d\n", sensi.element);
	fprintf(fp, "node: %5d\n", sensi.node);
	fprintf(fp, "v: ");
	print_trans_rotation3d(fp, sensi.v);
	fprintf(fp, "s: %15.7e\n", sensi.s);
}


void print_fem_sound_pressure(FILE *fp, FEM_SOUND_PRESSURE sound)
{
	fprintf(fp, "node: %5d\n", sound.node);
	fprintf(fp, "re: %15.7e\n", sound.re);
	fprintf(fp, "im: %15.7e\n", sound.im);

	print_i_info(fp, FEM_SOUND_PRESSURE_INFO, sound.i_info);
	print_d_info(fp, FEM_SOUND_PRESSURE_INFO, sound.d_info);
}


void print_fem_load_curve(FILE *fp, const FEM_LOAD_CURVE *curve)
{
	int ii1;

	fprintf(fp, "label: %5d\n", curve->label);
	fprintf(fp, "max_time:  %15.7e\n", curve->max_time);
	fprintf(fp, "max_value: %15.7e\n", curve->max_value);
	fprintf(fp, "total_step: %5d\n", curve->total_step);

	fprintf(fp, "inc_value:");
	for(ii1=0; ii1<(curve->total_step); ii1++){
		fprintf(fp, " %10.2e", curve->inc_value[ii1]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "approx_size: %5d\n", curve->approx_size);

	fprintf(fp, "approx_time:");
	for(ii1=0; ii1<(curve->approx_size); ii1++){
		fprintf(fp, " %10.2e\n", curve->approx_time[ii1]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "approx_value:");
	for(ii1=0; ii1<(curve->approx_size); ii1++){
		fprintf(fp, " %10.2e\n", curve->approx_value[ii1]);
	}
	fprintf(fp, "\n");
}


void print_fem_local_coord(FILE *fp, FEM_LOCAL_COORD local_coord)
{
	fprintf(fp, "label: %5d\n", local_coord.label);
	fprintf(fp, "origin: ");
	print_translation3d(fp, local_coord.origin);
	fprintf(fp, "d_cos[0]: ");
	print_translation3d(fp, local_coord.d_cos[0]);
	fprintf(fp, "d_cos[1]: ");
	print_translation3d(fp, local_coord.d_cos[1]);
	fprintf(fp, "d_cos[2]: ");
	print_translation3d(fp, local_coord.d_cos[2]);

	print_i_info(fp, FEM_LOCAL_COORD_INFO, local_coord.i_info);
	print_d_info(fp, FEM_LOCAL_COORD_INFO, local_coord.d_info);
}


void print_fem_elem_matrix(FILE *fp, FEM_ELEM_MATRIX elem_matrix)
{
	int ii1;

	fprintf(fp, "element: %5d\n", elem_matrix.element);
	fprintf(fp, "size: %5d\n", elem_matrix.size);
	fprintf(fp, "dim:  %5d\n", elem_matrix.dim);

	fprintf(fp, "index: ");
	if(elem_matrix.index != NULL){
		for(ii1=0; ii1<elem_matrix.size; ii1++){
			fprintf(fp, "%5d ", elem_matrix.index[ii1]);
		}
		fprintf(fp, "\n");
	}else{
		fprintf(fp, "NULL\n");
	}

	fprintf(fp, "matrix: \n");
	if(elem_matrix.matrix != NULL){
		print_matrix(fp, elem_matrix.size,
		                 elem_matrix.size, elem_matrix.matrix);
	}else{
		fprintf(fp, "NULL\n");
	}
}


void print_fem_elem_lookup(FILE *fp, FEM_ELEM_LOOKUP elem_lookup)
{
	int ii1;

	fprintf(fp, "node: %5d\n", elem_lookup.node);

	fprintf(fp, "size: %5d\n", elem_lookup.size);
	fprintf(fp, "alloc_size: %5d\n", elem_lookup.alloc_size);
	fprintf(fp, "element: ");
	if(elem_lookup.element != NULL){
		for(ii1=0; ii1<elem_lookup.size; ii1++){
			fprintf(fp, "%5d ", elem_lookup.element[ii1]);
		}
		fprintf(fp, "\n");
	}else{
		fprintf(fp, "NULL\n");
	}

	fprintf(fp, "r_size: %5d\n", elem_lookup.r_size);
	fprintf(fp, "r_alloc_size: %5d\n", elem_lookup.r_alloc_size);
	fprintf(fp, "r_element: ");
	if(elem_lookup.r_element != NULL){
		for(ii1=0; ii1<elem_lookup.r_size; ii1++){
			fprintf(fp, "%5d ", elem_lookup.r_element[ii1]);
		}
		fprintf(fp, "\n");
	}else{
		fprintf(fp, "NULL\n");
	}
}


void print_fem_data_type(FILE *fp, FEM_DATA_TYPE type)
{
	switch(type){
	case FEM_ANSYS:
		fprintf(fp, "FEM_ANSYS\n");
		break;
	case FEM_IDEAS:
		fprintf(fp, "FEM_IDEAS\n");
		break;
	case FEM_NASTRAN:
		fprintf(fp, "FEM_NASTRAN\n");
		break;
	case FEM_SYSNOISE:
		fprintf(fp, "FEM_SYSNOISE\n");
		break;
	case FEM_DXF:
		fprintf(fp, "FEM_DXF\n");
		break;
	case FEM_QUINT:
		fprintf(fp, "FEM_QUINT\n");
		break;
	case FEM_GENERATED_SURFACE:
		fprintf(fp, "FEM_GENERATED_SURFACE\n");
		break;
	case FEM_K_SOL:
		fprintf(fp, "FEM_K_SOL\n");
		break;
	case FEM_ADV:
		fprintf(fp, "FEM_ADV\n");
		break;
	case FEM_NIKE:
		fprintf(fp, "FEM_NIKE\n");
		break;
	case FEM_UNKNOWN:
		fprintf(fp, "UNKNOWN\n");
		break;
	default:
		fprintf(fp, "unknown\n");
	}
}

void print_fem_array_sort_flag(FILE *fp, FEM_ARRAY_SORT_FLAG sort_flag)
{
	fprintf(fp, "sort_flag : ");

	switch(sort_flag){
	case ARRAY_SORTED:
		fprintf(fp, "ARRAY_SORTED\n");
		break;
	case ARRAY_UNSORTED:
		fprintf(fp, "ARRAY_UNSORTED\n");
		break;
	default:
		fprintf(fp, "unknown\n");
	}
}

void print_fem_renumber_flag(FILE *fp, FEM_RENUMBER_FLAG renum_flag)
{
	fprintf(fp, "renum_flag : ");

	switch(renum_flag){
	case RENUMBERED:
		fprintf(fp, "RENUMBERED\n");
		break;
	case UNRENUMBERED:
		fprintf(fp, "UNRENUMBERED\n");
		break;
	default:
		fprintf(fp, "unknown\n");
	}
}

void print_fem_node_array(FILE *fp, FEM_NODE_ARRAY node)
{
	int ii1;

	fprintf(fp, "size: %10d\n", node.size);
	fprintf(fp, "alloc_size: %10d\n", node.alloc_size);
	print_fem_data_type(fp, node.source);

	if(node.array == NULL){
		fprintf(fp,"node_array: NULL\n");
	}else{
		for(ii1=0; ii1<node.size; ii1++){
			print_fem_node(fp, node.array[ii1]);
		}
	}

	print_fem_renumber_flag(fp, node.renum_flag);
	print_fem_array_sort_flag(fp, node.sort_flag);
}

void print_fem_element_array(FILE *fp, FEM_ELEMENT_ARRAY element)
{
	int ii1;

	fprintf(fp, "size: %10d\n", element.size);
	fprintf(fp, "alloc_size: %10d\n", element.alloc_size);
	print_fem_data_type(fp, element.source);

	if(element.array == NULL){
		fprintf(fp,"array: NULL\n");
	}else{
		for(ii1=0; ii1<element.size; ii1++){
			print_fem_element(fp, element.array[ii1]);
		}
	}

	print_fem_array_sort_flag(fp, element.sort_flag);
}	

void print_fem_rigid_element_array(FILE *fp, FEM_RIGID_ELEMENT_ARRAY element)
{
	int ii1;

	fprintf(fp, "size: %10d\n", element.size);
	fprintf(fp, "alloc_size: %10d\n", element.alloc_size);
	print_fem_data_type(fp, element.source);

	if(element.array == NULL){
		fprintf(fp,"array: NULL\n");
	}else{
		for(ii1=0; ii1<element.size; ii1++){
			print_fem_rigid_element(fp, &(element.array[ii1]));
		}
	}

	print_fem_array_sort_flag(fp, element.sort_flag);
}

void print_fem_material_prop_array(FILE *fp, FEM_MATERIAL_PROP_ARRAY m_prop)
{
	int ii1;

	fprintf(fp, "size: %10d\n", m_prop.size);
	fprintf(fp, "alloc_size: %10d\n", m_prop.alloc_size);
	print_fem_data_type(fp, m_prop.source);

	if(m_prop.array == NULL){
		fprintf(fp,"array: NULL\n");
	}else{
		for(ii1=0; ii1<m_prop.size; ii1++){
			print_fem_material_prop(fp, m_prop.array[ii1]);
		}
	}

	print_fem_array_sort_flag(fp, m_prop.sort_flag);
}

void print_fem_physical_prop_array(FILE *fp, FEM_PHYSICAL_PROP_ARRAY p_prop)
{
	int ii1;

	fprintf(fp, "size: %10d\n", p_prop.size);
	fprintf(fp, "alloc_size: %10d\n", p_prop.alloc_size);
	print_fem_data_type(fp, p_prop.source);

	if(p_prop.array == NULL){
		fprintf(fp,"array: NULL\n");
	}else{
		for(ii1=0; ii1<p_prop.size; ii1++){
			print_fem_physical_prop(fp, p_prop.array[ii1]);
		}
	}

	print_fem_array_sort_flag(fp, p_prop.sort_flag);
}

void print_fem_bc_array(FILE *fp, FEM_BC_ARRAY bc)
{
	int ii1;

	fprintf(fp, "size: %10d\n", bc.size);
	fprintf(fp, "alloc_size: %10d\n", bc.alloc_size);
	print_fem_data_type(fp, bc.source);

	if(bc.array == NULL){
		fprintf(fp,"bc_array: NULL\n");
	}else{
		for(ii1=0; ii1<bc.size; ii1++){
			print_fem_bc(fp, bc.array[ii1]);
		}
	}

	print_fem_array_sort_flag(fp, bc.sort_flag);
}

void print_fem_elem_bc_array(FILE *fp, FEM_ELEM_BC_ARRAY bc)
{
	int ii1;

	fprintf(fp, "size: %10d\n", bc.size);
	fprintf(fp, "alloc_size: %10d\n", bc.alloc_size);
	print_fem_data_type(fp, bc.source);

	if(bc.array == NULL){
		fprintf(fp,"bc_array: NULL\n");
	}else{
		for(ii1=0; ii1<bc.size; ii1++){
			print_fem_elem_bc(fp, bc.array[ii1]);
		}
	}

	print_fem_array_sort_flag(fp, bc.sort_flag);
}

void print_fem_disp_array(FILE *fp, FEM_DISP_ARRAY disp)
{
	int ii1;

	fprintf(fp, "size: %10d\n", disp.size);
	fprintf(fp, "alloc_size: %10d\n", disp.alloc_size);
	print_fem_data_type(fp, disp.source);

	if(disp.array == NULL){
		fprintf(fp,"disp_array: NULL\n");
	}else{
		for(ii1=0; ii1<disp.size; ii1++){
			print_fem_disp_velo(fp, disp.array[ii1]);
		}
	}

	print_fem_array_sort_flag(fp, disp.sort_flag);
}

void print_fem_velo_array(FILE *fp, FEM_VELO_ARRAY velo)
{
	int ii1;

	fprintf(fp, "size: %10d\n", velo.size);
	fprintf(fp, "alloc_size: %10d\n", velo.alloc_size);
	print_fem_data_type(fp, velo.source);

	if(velo.array == NULL){
		fprintf(fp,"velo_array: NULL\n");
	}else{
		for(ii1=0; ii1<velo.size; ii1++){
			print_fem_disp_velo(fp, velo.array[ii1]);
		}
	}

	print_fem_array_sort_flag(fp, velo.sort_flag);
}

void print_fem_react_array(FILE *fp, FEM_REACT_ARRAY react)
{
	int ii1;

	fprintf(fp, "size: %10d\n", react.size);
	fprintf(fp, "alloc_size: %10d\n", react.alloc_size);
	print_fem_data_type(fp, react.source);

	if(react.array == NULL){
		fprintf(fp,"react_array: NULL\n");
	}else{
		for(ii1=0; ii1<react.size; ii1++){
			print_fem_disp_velo(fp, react.array[ii1]);
		}
	}

	print_fem_array_sort_flag(fp, react.sort_flag);
}

void print_fem_stress_array(FILE *fp, FEM_STRESS_ARRAY stress)
{
	int ii1;

	fprintf(fp, "size: %10d\n", stress.size);
	fprintf(fp, "alloc_size: %10d\n", stress.alloc_size);
	print_fem_data_type(fp, stress.source);

	if(stress.array == NULL){
		fprintf(fp,"stress_array: NULL\n");
	}else{
		for(ii1=0; ii1<stress.size; ii1++){
			print_fem_stress_strain(fp, stress.array[ii1]);
		}
	}

	print_fem_array_sort_flag(fp, stress.sort_flag);
}

void print_fem_strain_array(FILE *fp, FEM_STRAIN_ARRAY strain)
{
	int ii1;

	fprintf(fp, "size: %10d\n", strain.size);
	fprintf(fp, "alloc_size: %10d\n", strain.alloc_size);
	print_fem_data_type(fp, strain.source);

	if(strain.array == NULL){
 		fprintf(fp,"strain_array: NULL\n");
	}else{
		for(ii1=0; ii1<strain.size; ii1++){
			print_fem_stress_strain(fp, strain.array[ii1]);
		}
	}

	print_fem_array_sort_flag(fp, strain.sort_flag);
}

void print_fem_sensitivity_type(FILE *fp, FEM_SENSITIVITY_TYPE type)
{
	switch(type){
	case SENS_SHAPE_GRAD_DENSITY:
		fprintf(fp, "SENS_SHAPE_GRAD_DENSITY\n");
		break;
	case SENS_SHAPE_GRAD:
		fprintf(fp, "SENS_SHAPE_GRAD\n");
		break;
	case SENS_THICKNESS:
		fprintf(fp, "SENS_THICKNESS\n");
		break;
	case SENS_BEAM_XY:
		fprintf(fp, "SENS_BEAM_XY\n");
		break;
	case SENS_UNKNOWN:
		fprintf(fp, "SENS_UNKNOWN\n");
		break;
	default:
		fprintf(fp, "unknown\n");
	}
}

void print_fem_sensitivity_array(FILE *fp, FEM_SENSITIVITY_ARRAY sensi)
{
	int ii1;

	fprintf(fp, "size: %10d\n", sensi.size);
	fprintf(fp, "alloc_size: %10d\n", sensi.alloc_size);
	print_fem_sensitivity_type(fp, sensi.type);

	if(sensi.array == NULL){
		fprintf(fp,"array: NULL");
	}else{
		for(ii1=0; ii1<sensi.size; ii1++){
			print_fem_sensitivity(fp, sensi.array[ii1]);
		}
	}

	print_fem_array_sort_flag(fp, sensi.sort_flag);
}


void print_fem_default_bc(FILE *fp, FEM_DEFAULT_BC def)
{
	fprintf(fp, "set_id: %5d\n", def.set_id);
	fprintf(fp, "temp: %15.7e\n", def.temp);
	fprintf(fp, "grav: ");
	print_translation3d(fp, def.grav);
	print_i_info(fp, FEM_DEFAULT_BC_INFO, def.i_info);
	print_d_info(fp, FEM_DEFAULT_BC_INFO, def.d_info);	
}


void print_fem_default_bc_array(FILE *fp, FEM_DEFAULT_BC_ARRAY def)
{
	int ii1;

	fprintf(fp, "size: %10d\n", def.size);
	fprintf(fp, "alloc_size: %10d\n", def.alloc_size);
	print_fem_data_type(fp, def.source);

	if(def.array == NULL){
		fprintf(fp,"array: NULL\n");
	}else{
		for(ii1=0; ii1<def.size; ii1++){
			print_fem_default_bc(fp, def.array[ii1]);
		}
	}

	print_fem_array_sort_flag(fp, def.sort_flag);
}


void print_fem_sound_pressure_array(FILE *fp, FEM_SOUND_PRESSURE_ARRAY sound)
{
	int ii1;

	fprintf(fp, "size: %10d\n", sound.size);
	fprintf(fp, "alloc_size: %10d\n", sound.alloc_size);
	print_fem_data_type(fp, sound.source);

	if(sound.array == NULL){
		fprintf(fp,"array: NULL\n");
	}else{
		for(ii1=0; ii1<sound.size; ii1++){
			print_fem_sound_pressure(fp, sound.array[ii1]);
		}
	}

	print_fem_array_sort_flag(fp, sound.sort_flag);
}


void print_fem_load_curve_array(FILE *fp, FEM_LOAD_CURVE_ARRAY curve)
{
	int ii1;

	fprintf(fp, "size: %10d\n", curve.size);
	fprintf(fp, "alloc_size: %10d\n", curve.alloc_size);
	print_fem_data_type(fp, curve.source);

	if(curve.array == NULL){
		fprintf(fp,"array: NULL\n");
	}else{
		for(ii1=0; ii1<curve.size; ii1++){
			print_fem_load_curve(fp, &(curve.array[ii1]));
		}
	}

	print_fem_array_sort_flag(fp, curve.sort_flag);
}


void print_fem_elem_matrix_array(FILE *fp, FEM_ELEM_MATRIX_ARRAY elem_matrix)
{
	int ii1;

	fprintf(fp, "size: %10d\n", elem_matrix.size);
	fprintf(fp, "alloc_size: %10d\n", elem_matrix.alloc_size);
	print_fem_data_type(fp, elem_matrix.source);

	if(elem_matrix.array == NULL){
		fprintf(fp,"array: NULL\n");
	}else{
		for(ii1=0; ii1<elem_matrix.size; ii1++){
			print_fem_elem_matrix(fp, elem_matrix.array[ii1]);
		}
	}

	print_fem_array_sort_flag(fp, elem_matrix.sort_flag);
}


void print_fem_elem_lookup_array(FILE *fp, FEM_ELEM_LOOKUP_ARRAY elem_lookup)
{
	int ii1;

	fprintf(fp, "size: %10d\n", elem_lookup.size);
	fprintf(fp, "alloc_size: %10d\n", elem_lookup.alloc_size);
	print_fem_data_type(fp, elem_lookup.source);

	if(elem_lookup.array == NULL){
		fprintf(fp,"array: NULL\n");
	}else{
		for(ii1=0; ii1<elem_lookup.size; ii1++){
			print_fem_elem_lookup(fp, elem_lookup.array[ii1]);
		}
	}

	print_fem_array_sort_flag(fp, elem_lookup.sort_flag);
}


void print_fem_local_coord_array(FILE *fp, FEM_LOCAL_COORD_ARRAY local_coord)
{
	int ii1;

	fprintf(fp, "size: %10d\n", local_coord.size);
	fprintf(fp, "alloc_size: %10d\n", local_coord.alloc_size);
	print_fem_data_type(fp, local_coord.source);

	if(local_coord.array == NULL){
		fprintf(fp,"array: NULL\n");
	}else{
		for(ii1=0; ii1<local_coord.size; ii1++){
			print_fem_local_coord(fp, local_coord.array[ii1]);
		}
	}

	print_fem_array_sort_flag(fp, local_coord.sort_flag);
}


void print_fem_bc_set(FILE *fp, FEM_BC_SET bc_set)
{
	int ii1;

	fprintf(fp, "label: %5d\n", bc_set.label);

	fprintf(fp, "num_fix: %5d\n", bc_set.num_fix);
	fprintf(fp, "scale_fix: ");
	for(ii1=0; ii1<FEM_MAX_BC_SET; ii1++){
		fprintf(fp, "%10.3e ", bc_set.scale_fix[ii1]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "id_fix: ");
	for(ii1=0; ii1<FEM_MAX_BC_SET; ii1++){
		fprintf(fp, "%5d ", bc_set.id_fix[ii1]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "num_mpc: %5d\n", bc_set.num_mpc);
	fprintf(fp, "scale_mpc: ");
	for(ii1=0; ii1<FEM_MAX_BC_SET; ii1++){
		fprintf(fp, "%10.3e ", bc_set.scale_mpc[ii1]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "id_mpc: ");
	for(ii1=0; ii1<FEM_MAX_BC_SET; ii1++){
		fprintf(fp, "%5d ", bc_set.id_mpc[ii1]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "num_force: %5d\n", bc_set.num_force);
	fprintf(fp, "scale_force: ");
	for(ii1=0; ii1<FEM_MAX_BC_SET; ii1++){
		fprintf(fp, "%10.3e ", bc_set.scale_force[ii1]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "id_force: ");
	for(ii1=0; ii1<FEM_MAX_BC_SET; ii1++){
		fprintf(fp, "%5d ", bc_set.id_force[ii1]);
	}
	fprintf(fp, "\n");	

	fprintf(fp, "num_temp: %5d\n", bc_set.num_temp);
	fprintf(fp, "scale_temp: ");
	for(ii1=0; ii1<FEM_MAX_BC_SET; ii1++){
		fprintf(fp, "%10.3e ", bc_set.scale_temp[ii1]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "id_temp: ");
	for(ii1=0; ii1<FEM_MAX_BC_SET; ii1++){
		fprintf(fp, "%5d ", bc_set.id_temp[ii1]);
	}
	fprintf(fp, "\n");

	print_i_info(fp, FEM_BC_SET_INFO, bc_set.i_info);
	print_d_info(fp, FEM_BC_SET_INFO, bc_set.d_info);	
}		


void print_fem_bc_common_flag(FILE *fp, FEM_BC_COMMON_FLAG flag)
{
	if(flag == FEM_BC_NOT_COMMON){
		fprintf(fp, "FEM_BC_NOT_COMMON\n");
	}else if(flag == FEM_BC_COMMON){
		fprintf(fp, "FEM_BC_COMMON\n");
	}
}


void print_fem_bc_set_array(FILE *fp, FEM_BC_SET_ARRAY bc_set)
{
	int ii1;

	fprintf(fp, "size: %10d\n", bc_set.size);
	fprintf(fp, "alloc_size: %10d\n", bc_set.alloc_size);
	print_fem_data_type(fp, bc_set.source);
	if(bc_set.array == NULL){
		fprintf(fp,"array: NULL\n");
		return;
	}

	fprintf(fp, "com_fix: ");
	print_fem_bc_common_flag(fp, bc_set.com_fix);
	fprintf(fp, "com_mpc: ");
	print_fem_bc_common_flag(fp, bc_set.com_mpc);
	fprintf(fp, "com_force: ");
	print_fem_bc_common_flag(fp, bc_set.com_force);
	fprintf(fp, "com_temp: ");
	print_fem_bc_common_flag(fp, bc_set.com_temp);

	for(ii1=0; ii1<bc_set.size; ii1++){
		print_fem_bc_set(fp, bc_set.array[ii1]);
	}
	print_fem_array_sort_flag(fp, bc_set.sort_flag);
}


void print_system_info(FILE *fp)
{
	fprintf(fp, "-- Integers --\n");
	fprintf(fp, "char[%d]\n", (int)(sizeof(char)));
	fprintf(fp, "short[%d]\n", (int)(sizeof(short)));
	fprintf(fp, "int[%d]\n", (int)(sizeof(int)));
	fprintf(fp, "long[%d]\n", (int)(sizeof(long)));
	/* fprintf(fp, "long long[%d]\n", (int)(sizeof(long long))); */
	fprintf(fp, "-- Floating points --\n");
	fprintf(fp, "float[%d]: %15.7e : %15.7e\n",
	        (int)(sizeof(float)), FLT_EPSILON, FLT_MAX);
	fprintf(fp, "double[%d]: %15.7e : %15.7e\n",
	        (int)(sizeof(double)), DBL_EPSILON, DBL_MAX);
	fprintf(fp, "long double[%d]: %15.7Le : %15.7Le\n",
	        (int)(sizeof(long double)),
	        (long double)LDBL_EPSILON, (long double)LDBL_MAX);
	fprintf(fp, "-- FEM structs --\n");
	fprintf(fp, "FEM_NODE[%d]\n", (int)(sizeof(FEM_NODE)));
	fprintf(fp, "FEM_ELEMENT[%d]\n", (int)(sizeof(FEM_ELEMENT)));
	fprintf(fp, "FEM_RIGID_ELEMENT[%d]\n", (int)(sizeof(FEM_RIGID_ELEMENT)));
	fprintf(fp, "FEM_MATERIAL_PROP[%d]\n", (int)(sizeof(FEM_MATERIAL_PROP)));
	fprintf(fp, "FEM_PHYSICAL_PROP[%d]\n", (int)(sizeof(FEM_PHYSICAL_PROP)));
	fprintf(fp, "FEM_BC[%d]\n", (int)(sizeof(FEM_BC)));
	fprintf(fp, "FEM_ELEM_BC[%d]\n", (int)(sizeof(FEM_ELEM_BC)));
	fprintf(fp, "FEM_DISP_VELO[%d]\n", (int)(sizeof(FEM_DISP_VELO)));
	fprintf(fp, "FEM_STRESS_STRAIN[%d]\n", (int)(sizeof(FEM_STRESS_STRAIN)));
	fprintf(fp, "FEM_SENSITIVITY[%d]\n", (int)(sizeof(FEM_SENSITIVITY)));
	fprintf(fp, "FEM_DEFAULT_BC[%d]\n", (int)(sizeof(FEM_DEFAULT_BC)));
	fprintf(fp, "FEM_BC_SET[%d]\n", (int)(sizeof(FEM_BC_SET)));
	fprintf(fp, "FEM_SOUND_PRESSURE[%d]\n", (int)(sizeof(FEM_SOUND_PRESSURE)));
}


void print_fem_node_matrix(FILE *fp, FEM_NODE_MATRIX nmat)
{
	int ii1, ii2;

	fprintf(fp, "node : %d\n", nmat.node);
	fprintf(fp, "shift : %d\n", nmat.shift);
	fprintf(fp, "size : %d\n", nmat.size);
	fprintf(fp, "alloc_size : %d\n", nmat.alloc_size);

	fprintf(fp, "col_node/col_shift/matrix :\n");
	for(ii1=0; ii1<nmat.size; ii1++){
		fprintf(fp, "%5d[%d]\n", nmat.col_node[ii1], nmat.col_shift[ii1]);
		for(ii2=0; ii2<3; ii2++){
			fprintf(fp, "%15.7e %15.7e %15.7e\n", nmat.matrix[ii1][ii2][0],
			        nmat.matrix[ii1][ii2][1], nmat.matrix[ii1][ii2][2]);
		}
	}
}


void print_fem_node_matrix_array(FILE *fp, FEM_NODE_MATRIX_ARRAY nmat)
{
	int ii1;

	fprintf(fp, "size: %10d\n", nmat.size);
	fprintf(fp, "alloc_size: %10d\n", nmat.alloc_size);
	print_fem_data_type(fp, nmat.source);

	if(nmat.array == NULL){
		fprintf(fp,"array: NULL\n");
	}else{
		for(ii1=0; ii1<nmat.size; ii1++){
			print_fem_node_matrix(fp, nmat.array[ii1]);
		}
	}

	print_fem_array_sort_flag(fp, nmat.sort_flag);
}



