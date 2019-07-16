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

/* $Id: fem_struct_binio.c 316 2004-06-21 03:54:58Z sasaoka $ */

#include <stdio.h>
#include <stdlib.h>
#include "rc.h"
#include "memory_manager.h"
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

FUNC_WRITE_ARRAY(write_node_array, NODE_ARRAY, NODE)

FUNC_WRITE_ARRAY(write_element_array, ELEMENT_ARRAY, ELEMENT)

FUNC_WRITE_ARRAY(write_material_prop_array,
                 MATERIAL_PROP_ARRAY, MATERIAL_PROP)

FUNC_WRITE_ARRAY(write_physical_prop_array,
                 PHYSICAL_PROP_ARRAY, PHYSICAL_PROP)

FUNC_WRITE_ARRAY(write_bc_array, BC_ARRAY, BC)

FUNC_WRITE_ARRAY(write_elem_bc_array, ELEM_BC_ARRAY, ELEM_BC)

FUNC_WRITE_ARRAY(write_disp_array, DISP_ARRAY, DISP_VELO)

FUNC_WRITE_ARRAY(write_velo_array, VELO_ARRAY, DISP_VELO)

FUNC_WRITE_ARRAY(write_react_array, REACT_ARRAY, DISP_VELO)

FUNC_WRITE_ARRAY(write_stress_array, STRESS_ARRAY, STRESS_STRAIN)

FUNC_WRITE_ARRAY(write_strain_array, STRAIN_ARRAY, STRESS_STRAIN)

FUNC_WRITE_ARRAY(write_sensitivity_array,
                 SENSITIVITY_ARRAY, SENSITIVITY)

FUNC_WRITE_ARRAY(write_default_bc_array,
                 DEFAULT_BC_ARRAY, DEFAULT_BC)

FUNC_WRITE_ARRAY(write_sound_pressure_array,
                 SOUND_PRESSURE_ARRAY, SOUND_PRESSURE)

FUNC_WRITE_ARRAY(write_bc_set_array, BC_SET_ARRAY, BC_SET)

FUNC_WRITE_ARRAY(write_load_curve_array, LOAD_CURVE_ARRAY,
                 LOAD_CURVE)

FUNC_WRITE_ARRAY(write_local_coord_array, LOCAL_COORD_ARRAY,
                 LOCAL_COORD)



#define FUNC_READ_ARRAY(name, array_type, type)  \
RC name(FILE *fp, array_type *arg){\
	if( fread(arg, sizeof(array_type), 1, fp) != (size_t)1){\
		return(READ_ERROR_RC);\
	}\
	if(arg->alloc_size > 0){\
		arg->array = mm_alloc(arg->alloc_size * sizeof(type));\
		if(arg->array == NULL) return(ALLOC_ERROR_RC);\
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


FUNC_READ_ARRAY(read_node_array, NODE_ARRAY, NODE)

FUNC_READ_ARRAY(read_element_array, ELEMENT_ARRAY, ELEMENT)

FUNC_READ_ARRAY(read_material_prop_array,
                MATERIAL_PROP_ARRAY, MATERIAL_PROP)

FUNC_READ_ARRAY(read_physical_prop_array,
                PHYSICAL_PROP_ARRAY, PHYSICAL_PROP)

FUNC_READ_ARRAY(read_bc_array, BC_ARRAY, BC)

FUNC_READ_ARRAY(read_elem_bc_array, ELEM_BC_ARRAY, ELEM_BC)

FUNC_READ_ARRAY(read_disp_array, DISP_ARRAY, DISP_VELO)

FUNC_READ_ARRAY(read_velo_array, VELO_ARRAY, DISP_VELO)

FUNC_READ_ARRAY(read_react_array, REACT_ARRAY, DISP_VELO)

FUNC_READ_ARRAY(read_stress_array, STRESS_ARRAY, STRESS_STRAIN)

FUNC_READ_ARRAY(read_strain_array, STRAIN_ARRAY, STRESS_STRAIN)

FUNC_READ_ARRAY(read_sensitivity_array,
                SENSITIVITY_ARRAY, SENSITIVITY)

FUNC_READ_ARRAY(read_default_bc_array,
                DEFAULT_BC_ARRAY, DEFAULT_BC)

FUNC_READ_ARRAY(read_sound_pressure_array,
                SOUND_PRESSURE_ARRAY, SOUND_PRESSURE)

FUNC_READ_ARRAY(read_bc_set_array, BC_SET_ARRAY, BC_SET)

FUNC_READ_ARRAY(read_load_curve_array, LOAD_CURVE_ARRAY,
                LOAD_CURVE)

FUNC_READ_ARRAY(read_local_coord_array, LOCAL_COORD_ARRAY,
                LOCAL_COORD)

