/*********************************************************************
 * fem_struct_binio.c
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

/* $Id: fem_struct_binio.c,v 1.10 2003/07/22 11:44:01 adachi Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include "rc.h"
#include "fem_struct.h"

/* local functions */

/* local variables */

#define FUNC_WRITE_ARRAY(name, array_type, type)  \
RC name(FILE *fp, array_type arg){\
	if( (arg.size > arg.alloc_size)||(arg.size < 0) ) return(ARG_ERROR_RC);\
	if( fwrite(&arg, sizeof(array_type), 1, fp) != (size_t)1 ){\
		return(WRITE_ERROR_RC);\
	}\
	if(arg.size > 0){\
		if( fwrite(arg.array, sizeof(type), arg.size, fp) != (size_t)arg.size){\
			return(WRITE_ERROR_RC);\
		}\
	}\
	return(NORMAL_RC);\
}

FUNC_WRITE_ARRAY(write_fem_node_array, FEM_NODE_ARRAY, FEM_NODE)

FUNC_WRITE_ARRAY(write_fem_element_array, FEM_ELEMENT_ARRAY, FEM_ELEMENT)

FUNC_WRITE_ARRAY(write_fem_material_prop_array,
                 FEM_MATERIAL_PROP_ARRAY, FEM_MATERIAL_PROP)

FUNC_WRITE_ARRAY(write_fem_physical_prop_array,
                 FEM_PHYSICAL_PROP_ARRAY, FEM_PHYSICAL_PROP)

FUNC_WRITE_ARRAY(write_fem_bc_array, FEM_BC_ARRAY, FEM_BC)

FUNC_WRITE_ARRAY(write_fem_elem_bc_array, FEM_ELEM_BC_ARRAY, FEM_ELEM_BC)

FUNC_WRITE_ARRAY(write_fem_disp_array, FEM_DISP_ARRAY, FEM_DISP_VELO)

FUNC_WRITE_ARRAY(write_fem_velo_array, FEM_VELO_ARRAY, FEM_DISP_VELO)

FUNC_WRITE_ARRAY(write_fem_react_array, FEM_REACT_ARRAY, FEM_DISP_VELO)

FUNC_WRITE_ARRAY(write_fem_stress_array, FEM_STRESS_ARRAY, FEM_STRESS_STRAIN)

FUNC_WRITE_ARRAY(write_fem_strain_array, FEM_STRAIN_ARRAY, FEM_STRESS_STRAIN)

FUNC_WRITE_ARRAY(write_fem_sensitivity_array,
                 FEM_SENSITIVITY_ARRAY, FEM_SENSITIVITY)

FUNC_WRITE_ARRAY(write_fem_default_bc_array,
                 FEM_DEFAULT_BC_ARRAY, FEM_DEFAULT_BC)

FUNC_WRITE_ARRAY(write_fem_sound_pressure_array,
                 FEM_SOUND_PRESSURE_ARRAY, FEM_SOUND_PRESSURE)

FUNC_WRITE_ARRAY(write_fem_bc_set_array, FEM_BC_SET_ARRAY, FEM_BC_SET)

FUNC_WRITE_ARRAY(write_fem_load_curve_array, FEM_LOAD_CURVE_ARRAY,
                 FEM_LOAD_CURVE)

FUNC_WRITE_ARRAY(write_fem_local_coord_array, FEM_LOCAL_COORD_ARRAY,
                 FEM_LOCAL_COORD)



#define FUNC_READ_ARRAY(name, array_type, type)  \
RC name(FILE *fp, array_type *arg){\
	if( fread(arg, sizeof(array_type), 1, fp) != (size_t)1){\
		return(READ_ERROR_RC);\
	}\
	if(arg->alloc_size > 0){\
		arg->array = malloc(arg->alloc_size * sizeof(type));\
		if(arg->array == NULL){\
			fprintf(stderr, "malloc() < ");\
			return(ALLOC_ERROR_RC);\
		}\
	}else{\
		arg->alloc_size = 0;\
		arg->array = NULL;\
	}\
	if(arg->size > arg->alloc_size) arg->size = arg->alloc_size;\
	if(arg->size < 0) arg->size = 0;\
	if(arg->size > 0){\
		if( fread(arg->array, sizeof(type), arg->size, fp)\
			                        != (size_t)(arg->size) ){\
			return(WRITE_ERROR_RC);\
		}\
	}\
	return(NORMAL_RC);\
}


FUNC_READ_ARRAY(read_fem_node_array, FEM_NODE_ARRAY, FEM_NODE)

FUNC_READ_ARRAY(read_fem_element_array, FEM_ELEMENT_ARRAY, FEM_ELEMENT)

FUNC_READ_ARRAY(read_fem_material_prop_array,
                FEM_MATERIAL_PROP_ARRAY, FEM_MATERIAL_PROP)

FUNC_READ_ARRAY(read_fem_physical_prop_array,
                FEM_PHYSICAL_PROP_ARRAY, FEM_PHYSICAL_PROP)

FUNC_READ_ARRAY(read_fem_bc_array, FEM_BC_ARRAY, FEM_BC)

FUNC_READ_ARRAY(read_fem_elem_bc_array, FEM_ELEM_BC_ARRAY, FEM_ELEM_BC)

FUNC_READ_ARRAY(read_fem_disp_array, FEM_DISP_ARRAY, FEM_DISP_VELO)

FUNC_READ_ARRAY(read_fem_velo_array, FEM_VELO_ARRAY, FEM_DISP_VELO)

FUNC_READ_ARRAY(read_fem_react_array, FEM_REACT_ARRAY, FEM_DISP_VELO)

FUNC_READ_ARRAY(read_fem_stress_array, FEM_STRESS_ARRAY, FEM_STRESS_STRAIN)

FUNC_READ_ARRAY(read_fem_strain_array, FEM_STRAIN_ARRAY, FEM_STRESS_STRAIN)

FUNC_READ_ARRAY(read_fem_sensitivity_array,
                FEM_SENSITIVITY_ARRAY, FEM_SENSITIVITY)

FUNC_READ_ARRAY(read_fem_default_bc_array,
                FEM_DEFAULT_BC_ARRAY, FEM_DEFAULT_BC)

FUNC_READ_ARRAY(read_fem_sound_pressure_array,
                FEM_SOUND_PRESSURE_ARRAY, FEM_SOUND_PRESSURE)

FUNC_READ_ARRAY(read_fem_bc_set_array, FEM_BC_SET_ARRAY, FEM_BC_SET)

FUNC_READ_ARRAY(read_fem_load_curve_array, FEM_LOAD_CURVE_ARRAY,
                FEM_LOAD_CURVE)

FUNC_READ_ARRAY(read_fem_local_coord_array, FEM_LOCAL_COORD_ARRAY,
                FEM_LOCAL_COORD)

