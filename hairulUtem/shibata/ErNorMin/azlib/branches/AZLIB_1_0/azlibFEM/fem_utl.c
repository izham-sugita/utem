/*********************************************************************
 * fem_utl.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Yasuo ASAGA> <Tomoyuki TSUBATA> <Masanobu SUYAMA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: fem_utl.c,v 1.24 2003/10/02 12:32:05 adachi Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "rc.h"
#include "math_utl.h"
#include "fem_struct.h"

#define REL_ERROR (DBL_EPSILON * 1.0e+7)
#define ABS_ERROR (DBL_MIN/DBL_EPSILON)

#define BB_NUM 8  /* 双方向バブルソートの処理回数 */
/*#define REALLOC_SIZE 1024 */
#define MAX_INC_SIZE  (4096)

/* local function */
static int compar_label(int label1, int label2);
static RC bbqsort(void *base, size_t nmemb, size_t size,
                  int (*compar)(const void *, const void *));
static int compar_node_label(const void *p1, const void *p2);
static int compar_node_renum(const void *p1, const void *p2);
static int compar_element_label(const void *p1, const void *p2);
static int compar_rigid_element_label(const void *p1, const void *p2);
static int compar_material_prop_label(const void *p1, const void *p2);
static int compar_physical_prop_label(const void *p1, const void *p2);
static int compar_bc_node(const void *p1, const void *p2);
static int compar_elem_bc_element(const void *p1, const void *p2);
static int compar_disp_velo_node(const void *p1, const void *p2);
static int compar_stress_strain_element(const void *p1, const void *p2);
static int compar_sensitivity_element(const void *p1, const void *p2);
static int compar_default_bc_set_id(const void *p1, const void *p2);
static int compar_sound_pressure_node(const void *p1, const void *p2);
static int compar_bc_set_label(const void *p1, const void *p2);
static int compar_elem_matrix_element(const void *p1, const void *p2);
static int compar_elem_lookup_node(const void *p1, const void *p2);
static int compar_node_matrix_node(const void *p1, const void *p2);
static int compar_load_curve_label(const void *p1, const void *p2);
static int compar_local_coord_label(const void *p1, const void *p2);
static int search(const void *base, size_t nmemb, size_t size,
                  int (*compar)(const void *, const void *),
                  const void *target, FEM_ARRAY_SORT_FLAG sort_flag);
/* static int zero2max(int val, int max); */
static RC overwrite_bc(FEM_BC_TYPE *type1, double *v1,
                       FEM_BC_TYPE type2,  double v2);


/* 配列の解放 */
#define FUNC_FREE_ARRAY(name, type)   \
RC name(type *arg){\
	if(arg->alloc_size > 0){\
		free(arg->array);\
	}\
	arg->array = NULL;\
	arg->alloc_size = 0;\
	arg->size = 0;\
	return(NORMAL_RC);\
}

FUNC_FREE_ARRAY(free_fem_node_array, FEM_NODE_ARRAY)

FUNC_FREE_ARRAY(free_fem_element_array, FEM_ELEMENT_ARRAY)

FUNC_FREE_ARRAY(free_fem_material_prop_array, FEM_MATERIAL_PROP_ARRAY)

FUNC_FREE_ARRAY(free_fem_physical_prop_array, FEM_PHYSICAL_PROP_ARRAY)

FUNC_FREE_ARRAY(free_fem_bc_array, FEM_BC_ARRAY)

FUNC_FREE_ARRAY(free_fem_elem_bc_array, FEM_ELEM_BC_ARRAY)

FUNC_FREE_ARRAY(free_fem_disp_array, FEM_DISP_ARRAY)

FUNC_FREE_ARRAY(free_fem_velo_array, FEM_VELO_ARRAY)

FUNC_FREE_ARRAY(free_fem_react_array, FEM_REACT_ARRAY)

FUNC_FREE_ARRAY(free_fem_stress_array, FEM_STRESS_ARRAY)

FUNC_FREE_ARRAY(free_fem_strain_array, FEM_STRAIN_ARRAY)

FUNC_FREE_ARRAY(free_fem_sensitivity_array, FEM_SENSITIVITY_ARRAY)

FUNC_FREE_ARRAY(free_fem_default_bc_array, FEM_DEFAULT_BC_ARRAY)

FUNC_FREE_ARRAY(free_fem_sound_pressure_array, FEM_SOUND_PRESSURE_ARRAY)

FUNC_FREE_ARRAY(free_fem_load_curve_array, FEM_LOAD_CURVE_ARRAY)

FUNC_FREE_ARRAY(free_fem_local_coord_array, FEM_LOCAL_COORD_ARRAY)


/* 配列の解放(各要素の解放が必要な場合用) */
#define FUNC_FREE_ARRAY_F(name, type, free_func)   \
RC name(type *arg){\
	int ii1;\
	if(arg->alloc_size > 0){\
		for(ii1=0; ii1<(arg->size); ii1++){\
			RC_TRY( free_func(&(arg->array[ii1])) );\
		}\
		free(arg->array);\
	}\
	arg->array = NULL;\
	arg->alloc_size = 0;\
	arg->size = 0;\
	return(NORMAL_RC);\
}

FUNC_FREE_ARRAY_F(free_fem_elem_matrix_array, FEM_ELEM_MATRIX_ARRAY, 
                  free_fem_elem_matrix)

FUNC_FREE_ARRAY_F(free_fem_elem_lookup_array, FEM_ELEM_LOOKUP_ARRAY,
                  free_fem_elem_lookup)

FUNC_FREE_ARRAY_F(free_fem_node_matrix_array, FEM_NODE_MATRIX_ARRAY,
                  free_fem_node_matrix)

FUNC_FREE_ARRAY_F(free_fem_rigid_element_array, FEM_RIGID_ELEMENT_ARRAY,
                  free_fem_rigid_element)

FUNC_FREE_ARRAY_F(free_fem_bc_set_array, FEM_BC_SET_ARRAY,
                  free_fem_bc_set)


/* 配列のコピー(コピー先の配列は既に確保されているものとする) */
#define FUNC_COPY_ARRAY(name, array_type, type)  \
RC name(array_type src, array_type *dest){\
	int ii1;\
	type *array_tmp;\
	int alloc_size_tmp;\
\
	if(src.size > dest->alloc_size) return(OVERFLOW_ERROR_RC);\
	if(dest->array == NULL) return(NULL_ERROR_RC);\
\
	array_tmp = dest->array;\
	alloc_size_tmp = dest->alloc_size;\
	*dest = src;\
	dest->array = array_tmp;\
	dest->alloc_size = alloc_size_tmp;\
\
	for(ii1=0; ii1<(dest->size); ii1++){\
		dest->array[ii1] = src.array[ii1];\
	}\
\
	return(NORMAL_RC);\
}

FUNC_COPY_ARRAY(copy_fem_node_array, FEM_NODE_ARRAY, FEM_NODE)

FUNC_COPY_ARRAY(copy_fem_element_array, FEM_ELEMENT_ARRAY, FEM_ELEMENT)

FUNC_COPY_ARRAY(copy_fem_material_prop_array,
                FEM_MATERIAL_PROP_ARRAY, FEM_MATERIAL_PROP)

FUNC_COPY_ARRAY(copy_fem_physical_prop_array,
                FEM_PHYSICAL_PROP_ARRAY, FEM_PHYSICAL_PROP)

FUNC_COPY_ARRAY(copy_fem_bc_array, FEM_BC_ARRAY, FEM_BC)

FUNC_COPY_ARRAY(copy_fem_elem_bc_array, FEM_ELEM_BC_ARRAY, FEM_ELEM_BC)

FUNC_COPY_ARRAY(copy_fem_disp_array, FEM_DISP_ARRAY, FEM_DISP_VELO)

FUNC_COPY_ARRAY(copy_fem_velo_array, FEM_VELO_ARRAY, FEM_DISP_VELO)

FUNC_COPY_ARRAY(copy_fem_react_array, FEM_REACT_ARRAY, FEM_DISP_VELO)

FUNC_COPY_ARRAY(copy_fem_stress_array, FEM_STRESS_ARRAY, FEM_STRESS_STRAIN)

FUNC_COPY_ARRAY(copy_fem_strain_array, FEM_STRAIN_ARRAY, FEM_STRESS_STRAIN)

FUNC_COPY_ARRAY(copy_fem_sensitivity_array,
                FEM_SENSITIVITY_ARRAY, FEM_SENSITIVITY)

FUNC_COPY_ARRAY(copy_fem_default_bc_array,
                FEM_DEFAULT_BC_ARRAY, FEM_DEFAULT_BC)

FUNC_COPY_ARRAY(copy_fem_sound_pressure_array,
                FEM_SOUND_PRESSURE_ARRAY, FEM_SOUND_PRESSURE)

FUNC_COPY_ARRAY(copy_fem_load_curve_array,
                FEM_LOAD_CURVE_ARRAY, FEM_LOAD_CURVE)

FUNC_COPY_ARRAY(copy_fem_local_coord_array,
                FEM_LOCAL_COORD_ARRAY, FEM_LOCAL_COORD)


#define FUNC_CLEAN_ARRAY(name, array_type, type, sort_func, elem)  \
RC name(array_type *ptr){\
	int ii1;\
	RC rc;\
\
	if(ptr->size > ptr->alloc_size) return(ARG_ERROR_RC);\
	if(ptr->size < 0) return(ARG_ERROR_RC);\
\
	rc = sort_func(ptr);\
	if(rc != NORMAL_RC){\
		fprintf(stderr, #sort_func": ");\
		return(rc);\
	}\
	for(ii1=(ptr->size)-1; ii1>=0; ii1--){\
		if(ptr->array[ii1].elem < 0){\
			ptr->size = ii1;\
		}else{\
			break;\
		}\
	}\
\
	if(ptr->size < ptr->alloc_size){\
		if(ptr->size > 0){\
			void *tmp;\
\
			tmp = realloc(ptr->array, ptr->size*sizeof(type));\
			if(tmp == NULL){\
				fprintf(stderr, "realloc < ");\
				return(ALLOC_ERROR_RC);\
			}\
			ptr->array = (type *)tmp;\
			ptr->alloc_size = ptr->size;\
		}else{\
			free(ptr->array);\
			ptr->array = NULL;\
			ptr->size = 0;\
			ptr->alloc_size = 0;\
		}\
	}\
\
	return(NORMAL_RC);\
}

FUNC_CLEAN_ARRAY(clean_fem_node_array, FEM_NODE_ARRAY, FEM_NODE,
                 sort_fem_node_array, label)

FUNC_CLEAN_ARRAY(clean_fem_element_array, FEM_ELEMENT_ARRAY, FEM_ELEMENT,
                 sort_fem_element_array, label)

FUNC_CLEAN_ARRAY(clean_fem_rigid_element_array, FEM_RIGID_ELEMENT_ARRAY,
                 FEM_RIGID_ELEMENT, sort_fem_rigid_element_array, label)

FUNC_CLEAN_ARRAY(clean_fem_material_prop_array, FEM_MATERIAL_PROP_ARRAY,
                 FEM_MATERIAL_PROP, sort_fem_material_prop_array, label)

FUNC_CLEAN_ARRAY(clean_fem_physical_prop_array, FEM_PHYSICAL_PROP_ARRAY,
                 FEM_PHYSICAL_PROP, sort_fem_physical_prop_array, label)

FUNC_CLEAN_ARRAY(clean_fem_bc_array, FEM_BC_ARRAY,
                 FEM_BC, sort_fem_bc_array, node)

FUNC_CLEAN_ARRAY(clean_fem_elem_bc_array, FEM_ELEM_BC_ARRAY,
                 FEM_ELEM_BC, sort_fem_elem_bc_array, element)

FUNC_CLEAN_ARRAY(clean_fem_disp_array, FEM_DISP_ARRAY,
                 FEM_DISP_VELO, sort_fem_disp_array, node)

FUNC_CLEAN_ARRAY(clean_fem_velo_array, FEM_VELO_ARRAY,
                 FEM_DISP_VELO, sort_fem_velo_array, node)

FUNC_CLEAN_ARRAY(clean_fem_react_array, FEM_REACT_ARRAY,
                 FEM_DISP_VELO, sort_fem_react_array, node)

FUNC_CLEAN_ARRAY(clean_fem_stress_array, FEM_STRESS_ARRAY,
                 FEM_STRESS_STRAIN, sort_fem_stress_array, element)

FUNC_CLEAN_ARRAY(clean_fem_strain_array, FEM_STRAIN_ARRAY,
                 FEM_STRESS_STRAIN, sort_fem_strain_array, element)

FUNC_CLEAN_ARRAY(clean_fem_sensitivity_array, FEM_SENSITIVITY_ARRAY,
                 FEM_SENSITIVITY, sort_fem_sensitivity_array, element)

FUNC_CLEAN_ARRAY(clean_fem_default_bc_array, FEM_DEFAULT_BC_ARRAY,
                 FEM_DEFAULT_BC, sort_fem_default_bc_array, set_id)

FUNC_CLEAN_ARRAY(clean_fem_sound_pressure_array, FEM_SOUND_PRESSURE_ARRAY,
                 FEM_SOUND_PRESSURE, sort_fem_sound_pressure_array, node)

FUNC_CLEAN_ARRAY(clean_fem_bc_set_array, FEM_BC_SET_ARRAY,
                 FEM_BC_SET, sort_fem_bc_set_array, label)

FUNC_CLEAN_ARRAY(clean_fem_elem_matrix_array, FEM_ELEM_MATRIX_ARRAY,
                 FEM_ELEM_MATRIX, sort_fem_elem_matrix_array, element)

FUNC_CLEAN_ARRAY(clean_fem_elem_lookup_array, FEM_ELEM_LOOKUP_ARRAY,
                 FEM_ELEM_LOOKUP, sort_fem_elem_lookup_array, node)

FUNC_CLEAN_ARRAY(clean_fem_node_matrix_array, FEM_NODE_MATRIX_ARRAY,
                 FEM_NODE_MATRIX, sort_fem_node_matrix_array, node)

FUNC_CLEAN_ARRAY(clean_fem_load_curve_array, FEM_LOAD_CURVE_ARRAY,
                 FEM_LOAD_CURVE, sort_fem_load_curve_array, label)

FUNC_CLEAN_ARRAY(clean_fem_local_coord_array, FEM_LOCAL_COORD_ARRAY,
                 FEM_LOCAL_COORD, sort_fem_local_coord_array, label)


/* realloc_*() は、最低でも要素１つ追加できるよう array を再確保する */
#define FUNC_REALLOC_ARRAY(name, array_type, type, init_func)  \
RC name(array_type *arg){\
	int ii1;\
	type *tmp;\
	int realloc_size;\
	int inc_size;\
\
	if((arg->size < 0)||(arg->alloc_size < 0)){\
		arg->size = 0;\
		arg->alloc_size = 0;\
		arg->array = NULL;\
	}\
\
	if(arg->size >= arg->alloc_size){\
		inc_size = arg->alloc_size;\
		if(inc_size <= 0) inc_size = 1;\
		if(inc_size > MAX_INC_SIZE) inc_size = MAX_INC_SIZE;\
		realloc_size = arg->alloc_size + inc_size;\
\
		tmp = (type *)realloc((void *)arg->array,\
		                      realloc_size*sizeof(type));\
		if(tmp == NULL){\
			fprintf(stderr, "realloc < ");\
			return(ALLOC_ERROR_RC);\
		}\
		arg->array = tmp;\
\
		for(ii1=arg->alloc_size; ii1<realloc_size; ii1++){\
			init_func( &((arg->array)[ii1]) );\
		}\
		arg->alloc_size = realloc_size;\
	}\
\
	return(NORMAL_RC);\
}

FUNC_REALLOC_ARRAY(realloc_fem_node_array, FEM_NODE_ARRAY,
                   FEM_NODE, init_fem_node)

FUNC_REALLOC_ARRAY(realloc_fem_element_array, FEM_ELEMENT_ARRAY,
                   FEM_ELEMENT, init_fem_element)

FUNC_REALLOC_ARRAY(realloc_fem_rigid_element_array, FEM_RIGID_ELEMENT_ARRAY,
                   FEM_RIGID_ELEMENT, init_fem_rigid_element)

FUNC_REALLOC_ARRAY(realloc_fem_material_prop_array, FEM_MATERIAL_PROP_ARRAY,
                   FEM_MATERIAL_PROP, init_fem_material_prop)

FUNC_REALLOC_ARRAY(realloc_fem_physical_prop_array, FEM_PHYSICAL_PROP_ARRAY,
                   FEM_PHYSICAL_PROP, init_fem_physical_prop)

FUNC_REALLOC_ARRAY(realloc_fem_bc_array, FEM_BC_ARRAY,
                   FEM_BC, init_fem_bc)

FUNC_REALLOC_ARRAY(realloc_fem_elem_bc_array, FEM_ELEM_BC_ARRAY,
                   FEM_ELEM_BC, init_fem_elem_bc)

FUNC_REALLOC_ARRAY(realloc_fem_disp_array, FEM_DISP_ARRAY,
                   FEM_DISP_VELO, init_fem_disp_velo)

FUNC_REALLOC_ARRAY(realloc_fem_velo_array, FEM_VELO_ARRAY,
                   FEM_DISP_VELO, init_fem_disp_velo)

FUNC_REALLOC_ARRAY(realloc_fem_react_array, FEM_REACT_ARRAY,
                   FEM_DISP_VELO, init_fem_disp_velo)

FUNC_REALLOC_ARRAY(realloc_fem_stress_array, FEM_STRESS_ARRAY,
                   FEM_STRESS_STRAIN, init_fem_stress_strain)

FUNC_REALLOC_ARRAY(realloc_fem_strain_array, FEM_STRAIN_ARRAY,
                   FEM_STRESS_STRAIN, init_fem_stress_strain)

FUNC_REALLOC_ARRAY(realloc_fem_sensitivity_array, FEM_SENSITIVITY_ARRAY,
                   FEM_SENSITIVITY, init_fem_sensitivity)

FUNC_REALLOC_ARRAY(realloc_fem_default_bc_array, FEM_DEFAULT_BC_ARRAY,
                   FEM_DEFAULT_BC, init_fem_default_bc)

FUNC_REALLOC_ARRAY(realloc_fem_sound_pressure_array, FEM_SOUND_PRESSURE_ARRAY,
                   FEM_SOUND_PRESSURE, init_fem_sound_pressure)

FUNC_REALLOC_ARRAY(realloc_fem_bc_set_array, FEM_BC_SET_ARRAY,
                   FEM_BC_SET, init_fem_bc_set)

FUNC_REALLOC_ARRAY(realloc_fem_elem_matrix_array, FEM_ELEM_MATRIX_ARRAY,
                   FEM_ELEM_MATRIX, init_fem_elem_matrix)

FUNC_REALLOC_ARRAY(realloc_fem_elem_lookup_array, FEM_ELEM_LOOKUP_ARRAY,
                   FEM_ELEM_LOOKUP, init_fem_elem_lookup)

FUNC_REALLOC_ARRAY(realloc_fem_node_matrix_array, FEM_NODE_MATRIX_ARRAY,
                   FEM_NODE_MATRIX, init_fem_node_matrix)

FUNC_REALLOC_ARRAY(realloc_fem_load_curve_array, FEM_LOAD_CURVE_ARRAY,
                   FEM_LOAD_CURVE, init_fem_load_curve)

FUNC_REALLOC_ARRAY(realloc_fem_local_coord_array, FEM_LOCAL_COORD_ARRAY,
                   FEM_LOCAL_COORD, init_fem_local_coord)


RC allocate_fem_node_array(int num, FEM_NODE_ARRAY *node)
{
	int ii1;

	if(num > 0){
		node->array = malloc(num * sizeof(FEM_NODE));
		if(node->array == NULL){
			fprintf(stderr,"malloc < ");
			return(ALLOC_ERROR_RC);
		}
	}else{
		node->array = NULL;
	}

	node->size = 0;
	node->alloc_size = num;
	node->source = FEM_UNKNOWN;
	node->renum_flag = UNRENUMBERED;
	node->sort_flag = ARRAY_UNSORTED;

	for(ii1=0;ii1<num;ii1++){
		init_fem_node( &((node->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC allocate_fem_element_array(int num, FEM_ELEMENT_ARRAY *element)
{
	int ii1;

	if(num > 0){
		element->array = malloc(num * sizeof(FEM_ELEMENT));
		if(element->array == NULL){
			fprintf(stderr,"malloc < ");
			return(ALLOC_ERROR_RC);
		}
	}else{
		element->array = NULL;
	}

	element->size = 0;
	element->alloc_size = num;
	element->source = FEM_UNKNOWN;
	element->sort_flag = ARRAY_UNSORTED;

	for(ii1=0;ii1<num;ii1++){
		init_fem_element( &((element->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC allocate_fem_rigid_element_array(int num, FEM_RIGID_ELEMENT_ARRAY *element)
{
	int ii1;

	if(num > 0){
		element->array = malloc(num * sizeof(FEM_RIGID_ELEMENT));
		if(element->array == NULL){
			fprintf(stderr,"malloc < ");
			return(ALLOC_ERROR_RC);
		}
	}else{
		element->array = NULL;
	}

	element->size = 0;
	element->alloc_size = num;
	element->source = FEM_UNKNOWN;
	element->sort_flag = ARRAY_UNSORTED;

	for(ii1=0;ii1<num;ii1++){
		init_fem_rigid_element( &((element->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC allocate_fem_material_prop_array(int num, FEM_MATERIAL_PROP_ARRAY *mat)
{
	int ii1;

	if(num > 0){
		mat->array = malloc(num * sizeof(FEM_MATERIAL_PROP));
		if(mat->array == NULL){
			fprintf(stderr,"malloc < ");
			return(ALLOC_ERROR_RC);
		}
	}else{
		mat->array = NULL;
	}

	mat->size = 0;
	mat->alloc_size = num;
	mat->source = FEM_UNKNOWN;
	mat->sort_flag = ARRAY_UNSORTED;

	for(ii1=0;ii1<num;ii1++){
		init_fem_material_prop( &((mat->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC allocate_fem_physical_prop_array(int num, FEM_PHYSICAL_PROP_ARRAY *phys)
{
	int ii1;

	if(num > 0){
		phys->array = malloc(num * sizeof(FEM_PHYSICAL_PROP));
		if(phys->array == NULL){
			fprintf(stderr,"malloc < ");
			return(ALLOC_ERROR_RC);
		}
	}else{
		phys->array = NULL;
	}

	phys->size = 0;
	phys->alloc_size = num;
	phys->source = FEM_UNKNOWN;
	phys->sort_flag = ARRAY_UNSORTED;

	for(ii1=0;ii1<num;ii1++){
		init_fem_physical_prop( &((phys->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC allocate_fem_bc_array(int num, FEM_BC_ARRAY *bc)
{
	int ii1;

	if(num > 0){
		bc->array = malloc(num * sizeof(FEM_BC));
		if(bc->array == NULL){
			fprintf(stderr,"malloc < ");
			return(ALLOC_ERROR_RC);
		}
	}else{
		bc->array = NULL;
	}

	bc->size = 0;
	bc->alloc_size = num;
	bc->source = FEM_UNKNOWN;
	bc->sort_flag = ARRAY_UNSORTED;

	for(ii1=0;ii1<num;ii1++){
		init_fem_bc( &((bc->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC allocate_fem_elem_bc_array(int num, FEM_ELEM_BC_ARRAY *bc)
{
	int ii1;

	if(num > 0){
		bc->array = malloc(num * sizeof(FEM_ELEM_BC));
		if(bc->array == NULL){
			fprintf(stderr,"malloc < ");
			return(ALLOC_ERROR_RC);
		}
	}else{
		bc->array = NULL;
	}

	bc->size = 0;
	bc->alloc_size = num;
	bc->source = FEM_UNKNOWN;
	bc->sort_flag = ARRAY_UNSORTED;

	for(ii1=0;ii1<num;ii1++){
		init_fem_elem_bc( &((bc->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC allocate_fem_disp_array(int num, FEM_DISP_ARRAY *disp)
{
	int ii1;

	if(num > 0){
		disp->array = malloc(num * sizeof(FEM_DISP_VELO));
		if(disp->array == NULL){
			fprintf(stderr,"malloc < ");
			return(ALLOC_ERROR_RC);
		}
	}else{
		disp->array = NULL;
	}

	disp->size = 0;
	disp->alloc_size = num;
	disp->source = FEM_UNKNOWN;
	disp->sort_flag = ARRAY_UNSORTED;

	for(ii1=0;ii1<num;ii1++){
		init_fem_disp_velo( &((disp->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC allocate_fem_velo_array(int num, FEM_VELO_ARRAY *velo)
{
	int ii1;

	if(num > 0){
		velo->array = malloc(num * sizeof(FEM_DISP_VELO));
		if(velo->array == NULL){
			fprintf(stderr,"malloc < ");
			return(ALLOC_ERROR_RC);
		}
	}else{
		velo->array = NULL;
	}
	
	velo->size = 0;
	velo->alloc_size = num;
	velo->source = FEM_UNKNOWN;
	velo->sort_flag = ARRAY_UNSORTED;

	for(ii1=0;ii1<num;ii1++){
		init_fem_disp_velo( &((velo->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC allocate_fem_react_array(int num, FEM_REACT_ARRAY *react)
{
	int ii1;

	if(num > 0){
		react->array = malloc(num * sizeof(FEM_DISP_VELO));
		if(react->array == NULL){
			fprintf(stderr,"malloc < ");
			return(ALLOC_ERROR_RC);
		}
	}else{
		react->array = NULL;
	}
	
	react->size = 0;
	react->alloc_size = num;
	react->source = FEM_UNKNOWN;
	react->sort_flag = ARRAY_UNSORTED;

	for(ii1=0;ii1<num;ii1++){
		init_fem_disp_velo( &((react->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC allocate_fem_stress_array(int num, FEM_STRESS_ARRAY *stress)
{
	int ii1;

	if(num > 0){
		stress->array = malloc(num * sizeof(FEM_STRESS_STRAIN));
		if(stress->array == NULL){
			fprintf(stderr,"malloc < ");
			return(ALLOC_ERROR_RC);
		}
	}else{
		stress->array = NULL;
	}

	stress->size = 0;
	stress->alloc_size = num;
	stress->source = FEM_UNKNOWN;
	stress->sort_flag = ARRAY_UNSORTED;

	for(ii1=0;ii1<num;ii1++){
		init_fem_stress_strain( &((stress->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC allocate_fem_strain_array(int num, FEM_STRAIN_ARRAY *strain)
{
	int ii1;

	if(num > 0){
		strain->array = malloc(num * sizeof(FEM_STRESS_STRAIN));
		if(strain->array == NULL){
			fprintf(stderr,"malloc < ");
			return(ALLOC_ERROR_RC);
		}
	}else{
		strain->array = NULL;
	}

	strain->size = 0;
	strain->alloc_size = num;
	strain->source = FEM_UNKNOWN;
	strain->sort_flag = ARRAY_UNSORTED;

	for(ii1=0;ii1<num;ii1++){
		init_fem_stress_strain( &((strain->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC allocate_fem_sensitivity_array(int num, FEM_SENSITIVITY_ARRAY *sens)
{
	int ii1;

	if(num > 0){
		sens->array = malloc(num * sizeof(FEM_SENSITIVITY));
		if(sens->array == NULL){
			fprintf(stderr,"malloc < ");
			return(ALLOC_ERROR_RC);
		}
	}else{
		sens->array = NULL;
	}

	sens->size = 0;
	sens->alloc_size = num;
	sens->type = SENS_UNKNOWN;
	sens->sort_flag = ARRAY_UNSORTED;

	for(ii1=0;ii1<num;ii1++){
		init_fem_sensitivity( &((sens->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC allocate_fem_default_bc_array(int num, FEM_DEFAULT_BC_ARRAY *def)
{
	int ii1;

	if(num > 0){
		def->array = malloc(num * sizeof(FEM_DEFAULT_BC));
		if(def->array == NULL){
			fprintf(stderr,"malloc < ");
			return(ALLOC_ERROR_RC);
		}
	}else{
		def->array = NULL;
	}

	def->size = 0;
	def->alloc_size = num;
	def->source = FEM_UNKNOWN;
	def->sort_flag = ARRAY_UNSORTED;

	for(ii1=0; ii1<num; ii1++){
		init_fem_default_bc( &((def->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC allocate_fem_sound_pressure_array(int num, FEM_SOUND_PRESSURE_ARRAY *sound)
{
	int ii1;

	if(num > 0){
		sound->array = malloc(num * sizeof(FEM_SOUND_PRESSURE));
		if(sound->array == NULL){
			fprintf(stderr,"malloc < ");
			return(ALLOC_ERROR_RC);
		}
	}else{
		sound->array = NULL;
	}

	sound->size = 0;
	sound->alloc_size = num;
	sound->source = FEM_UNKNOWN;
	sound->sort_flag = ARRAY_UNSORTED;

	for(ii1=0; ii1<num; ii1++){
		init_fem_sound_pressure( &((sound->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC allocate_fem_bc_set_array(int num, FEM_BC_SET_ARRAY *bc_set)
{
	int ii1;

	if(num > 0){
		bc_set->array = malloc(num * sizeof(FEM_BC_SET));
		if(bc_set->array == NULL){
			fprintf(stderr,"malloc < ");
			return(ALLOC_ERROR_RC);
		}
	}else{
		bc_set->array = NULL;
	}

	bc_set->size = 0;
	bc_set->alloc_size = num;
	bc_set->source = FEM_UNKNOWN;
	bc_set->sort_flag = ARRAY_UNSORTED;
	bc_set->com_fix = FEM_BC_NOT_COMMON;
	bc_set->com_mpc = FEM_BC_NOT_COMMON;
	bc_set->com_force = FEM_BC_NOT_COMMON;
	bc_set->com_temp = FEM_BC_NOT_COMMON;

	for(ii1=0; ii1<num; ii1++){
		init_fem_bc_set( &((bc_set->array)[ii1]) );
	}
	RC_TRY( init_string_array(&(bc_set->s_info)) );

	return(NORMAL_RC);
}


RC allocate_fem_elem_matrix_array(int num, FEM_ELEM_MATRIX_ARRAY *elem_matrix)
{
	int ii1;

	if(num > 0){
		elem_matrix->array = malloc(num * sizeof(FEM_ELEM_MATRIX));
		if(elem_matrix->array == NULL){
			fprintf(stderr,"malloc() < ");
			return(ALLOC_ERROR_RC);
		}
	}else{
		elem_matrix->array = NULL;
	}

	elem_matrix->size = 0;
	elem_matrix->alloc_size = num;
	elem_matrix->source = FEM_UNKNOWN;
	elem_matrix->sort_flag = ARRAY_UNSORTED;

	for(ii1=0; ii1<num; ii1++){
		init_fem_elem_matrix( &((elem_matrix->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC allocate_fem_elem_lookup_array(int num, FEM_ELEM_LOOKUP_ARRAY *elem_lookup)
{
	int ii1;

	if(num > 0){
		elem_lookup->array = malloc(num * sizeof(FEM_ELEM_LOOKUP));
		if(elem_lookup->array == NULL){
			fprintf(stderr,"malloc() < ");
			return(ALLOC_ERROR_RC);
		}
	}else{
		elem_lookup->array = NULL;
	}

	elem_lookup->size = 0;
	elem_lookup->alloc_size = num;
	elem_lookup->source = FEM_UNKNOWN;
	elem_lookup->sort_flag = ARRAY_UNSORTED;

	for(ii1=0; ii1<num; ii1++){
		init_fem_elem_lookup( &((elem_lookup->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC allocate_fem_node_matrix_array(int num, FEM_NODE_MATRIX_ARRAY *nmat)
{
	int ii1;

	if(num > 0){
		nmat->array = malloc(num * sizeof(FEM_NODE_MATRIX));
		RC_NULL_CHK( nmat->array );
	}else{
		nmat->array = NULL;
	}

	nmat->size = 0;
	nmat->alloc_size = num;
	nmat->source = FEM_UNKNOWN;
	nmat->sort_flag = ARRAY_UNSORTED;

	for(ii1=0; ii1<num; ii1++){
		init_fem_node_matrix( &((nmat->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC allocate_fem_load_curve_array(int num, FEM_LOAD_CURVE_ARRAY *curve)
{
	int ii1;

	if(num > 0){
		curve->array = malloc(num * sizeof(FEM_LOAD_CURVE));
		if(curve->array == NULL){
			fprintf(stderr,"malloc() < ");
			return(ALLOC_ERROR_RC);
		}
	}else{
		curve->array = NULL;
	}

	curve->size = 0;
	curve->alloc_size = num;
	curve->source = FEM_UNKNOWN;
	curve->sort_flag = ARRAY_UNSORTED;

	for(ii1=0; ii1<num; ii1++){
		init_fem_load_curve( &((curve->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC allocate_fem_local_coord_array(int num, FEM_LOCAL_COORD_ARRAY *local)
{
	int ii1;

	if(num > 0){
		local->array = malloc(num * sizeof(FEM_LOCAL_COORD));
		if(local->array == NULL){
			fprintf(stderr,"malloc() < ");
			return(ALLOC_ERROR_RC);
		}
	}else{
		local->array = NULL;
	}

	local->size = 0;
	local->alloc_size = num;
	local->source = FEM_UNKNOWN;
	local->sort_flag = ARRAY_UNSORTED;

	for(ii1=0; ii1<num; ii1++){
		init_fem_local_coord( &((local->array)[ii1]) );
	}

	return(NORMAL_RC);
}


/* count_valid_*() は，配列中で実際に有効な要素数を返却 */
#define FUNC_COUNT_VALID(name, array_type, elem)  \
int name(array_type arg){\
	int ii1;\
	int ret = 0;\
	for(ii1=0; ii1<arg.size; ii1++){\
		if(arg.array[ii1].elem >= 0){\
			ret++;\
		}\
	}\
	return(ret);\
}

FUNC_COUNT_VALID(count_valid_node, FEM_NODE_ARRAY, label)

FUNC_COUNT_VALID(count_valid_element, FEM_ELEMENT_ARRAY, label)

FUNC_COUNT_VALID(count_valid_rigid_element, FEM_RIGID_ELEMENT_ARRAY, label)

FUNC_COUNT_VALID(count_valid_material_prop, FEM_MATERIAL_PROP_ARRAY, label)

FUNC_COUNT_VALID(count_valid_physical_prop, FEM_PHYSICAL_PROP_ARRAY, label)

FUNC_COUNT_VALID(count_valid_bc, FEM_BC_ARRAY, node)

FUNC_COUNT_VALID(count_valid_elem_bc, FEM_ELEM_BC_ARRAY, element)

FUNC_COUNT_VALID(count_valid_disp, FEM_DISP_ARRAY, node)

FUNC_COUNT_VALID(count_valid_velo, FEM_VELO_ARRAY, node)

FUNC_COUNT_VALID(count_valid_react, FEM_REACT_ARRAY, node)

FUNC_COUNT_VALID(count_valid_load_curve, FEM_LOAD_CURVE_ARRAY, label)

FUNC_COUNT_VALID(count_valid_local_coord, FEM_LOCAL_COORD_ARRAY, label)


/* special comparison of label */
static int compar_label(int label1, int label2)
{
	if(label1 < 0){
		if(label2 < 0){
			return(0);
		}else{
			return(1);
		}
	}else{
		if(label2 < 0){
			return(-1);
		}
	}

	return(label1 - label2);
}


/* qsort() 互換のソートプログラム */
/* 大部分がソートされたデータを想定して、双方向バブルソートを数回行う */
/* それでも完全にソートできない場合にのみ qsort() を使う */
static RC bbqsort(void *base, size_t nmemb, size_t size,
                  int (*compar)(const void *, const void *))
{
	int ii1, ii2;
	int swap_flag;
	void *swap_area;

	if((swap_area = malloc(size)) == NULL){
		fprintf(stderr, "malloc < ");
		return(ALLOC_ERROR_RC);
	}
	
	swap_flag = 0;
	for(ii1=0;ii1<BB_NUM;ii1++){
		/* 正方向 */
		ii2 = swap_flag + 1;
		swap_flag = 0;
		for(;ii2<((int)nmemb-ii1);ii2++){
			void *elem1 = (void *)((char *)base + (ii2 - 1)*size);
			void *elem2 = (void *)((char *)base + ii2*size);

			if((*compar)(elem1, elem2) > 0){
				memcpy(swap_area, elem1, size);
				memcpy(elem1, elem2, size);
				memcpy(elem2, swap_area, size);
				swap_flag = ii2;
			}
		}
		if(swap_flag == 0) break;

		/* 逆方向 */
		ii2 = swap_flag - 1;
		swap_flag = 0;
		for(;ii2>ii1;ii2--){
			void *elem1 = (void *)((char *)base + (ii2 - 1)*size);
			void *elem2 = (void *)((char *)base + ii2*size);

			if((*compar)(elem1, elem2) > 0){
				memcpy(swap_area, elem1, size);
				memcpy(elem1, elem2, size);
				memcpy(elem2, swap_area, size);
				swap_flag = ii2;
			}
		}
		if(swap_flag == 0) break;
	}
	free(swap_area);
	if(swap_flag == 0) return(NORMAL_RC);
	
	nmemb -= 2*BB_NUM;
	if(nmemb <= 1) return(NORMAL_RC);
	qsort((void *)((char *)base + BB_NUM*size), nmemb, size, compar);

	return(NORMAL_RC);
}


RC sort_fem_node_array(FEM_NODE_ARRAY *node)
{
	RC_TRY( bbqsort(node->array, node->size, sizeof(FEM_NODE),
	                                    compar_node_label) );
	node->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC sort_fem_node_renum(FEM_NODE_ARRAY *node)
{
	RC_TRY( bbqsort(node->array, node->size, sizeof(FEM_NODE),
	                                    compar_node_renum) );
	node->sort_flag = ARRAY_UNSORTED;
	
	return(NORMAL_RC);
}


RC sort_fem_element_array(FEM_ELEMENT_ARRAY *elem)
{
	RC_TRY( bbqsort(elem->array, elem->size, sizeof(FEM_ELEMENT),
	                                    compar_element_label) );
	elem->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC sort_fem_rigid_element_array(FEM_RIGID_ELEMENT_ARRAY *elem)
{
	RC_TRY( bbqsort(elem->array, elem->size, sizeof(FEM_RIGID_ELEMENT),
	                                    compar_rigid_element_label) );
	elem->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC sort_fem_material_prop_array(FEM_MATERIAL_PROP_ARRAY *mat)
{
	RC_TRY( bbqsort(mat->array, mat->size, sizeof(FEM_MATERIAL_PROP),
	                                       compar_material_prop_label) );
	mat->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC sort_fem_physical_prop_array(FEM_PHYSICAL_PROP_ARRAY *phys)
{
	RC_TRY( bbqsort(phys->array, phys->size, sizeof(FEM_PHYSICAL_PROP),
	                                         compar_physical_prop_label) );
	phys->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC sort_fem_bc_array(FEM_BC_ARRAY *bc)
{
	RC_TRY( bbqsort(bc->array, bc->size, sizeof(FEM_BC),
	                                 compar_bc_node) );
	bc->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC sort_fem_elem_bc_array(FEM_ELEM_BC_ARRAY *bc)
{
	RC_TRY( bbqsort(bc->array, bc->size, sizeof(FEM_ELEM_BC),
	                                 compar_elem_bc_element) );
	bc->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC sort_fem_disp_array(FEM_DISP_ARRAY *disp)
{
	RC_TRY( bbqsort(disp->array, disp->size, sizeof(FEM_DISP_VELO),
	                                     compar_disp_velo_node) );
	disp->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC sort_fem_velo_array(FEM_VELO_ARRAY *velo)
{
	RC_TRY( bbqsort(velo->array, velo->size, sizeof(FEM_DISP_VELO),
	                                     compar_disp_velo_node) );
	velo->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC sort_fem_react_array(FEM_REACT_ARRAY *react)
{
	RC_TRY( bbqsort(react->array, react->size, sizeof(FEM_DISP_VELO),
	                                     compar_disp_velo_node) );
	react->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC sort_fem_stress_array(FEM_STRESS_ARRAY *stress)
{
	RC_TRY( bbqsort(stress->array, stress->size, sizeof(FEM_STRESS_STRAIN),
	                                      compar_stress_strain_element) );
	stress->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC sort_fem_strain_array(FEM_STRAIN_ARRAY *strain)
{
	RC_TRY( bbqsort(strain->array, strain->size, sizeof(FEM_STRESS_STRAIN),
	                                      compar_stress_strain_element) );
	strain->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC sort_fem_sensitivity_array(FEM_SENSITIVITY_ARRAY *sens)
{
	RC_TRY( bbqsort(sens->array, sens->size, sizeof(FEM_SENSITIVITY),
	                                  compar_sensitivity_element) );
	sens->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC sort_fem_default_bc_array(FEM_DEFAULT_BC_ARRAY *def)
{
	RC_TRY( bbqsort(def->array, def->size, sizeof(FEM_DEFAULT_BC),
	                                  compar_default_bc_set_id) );
	def->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC sort_fem_sound_pressure_array(FEM_SOUND_PRESSURE_ARRAY *sound)
{
	RC_TRY( bbqsort(sound->array, sound->size, sizeof(FEM_SOUND_PRESSURE),
	                                  compar_sound_pressure_node) );
	sound->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC sort_fem_bc_set_array(FEM_BC_SET_ARRAY *bc_set)
{
	RC_TRY( bbqsort(bc_set->array, bc_set->size, sizeof(FEM_BC_SET),
	                                  compar_bc_set_label) );
	bc_set->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC sort_fem_elem_matrix_array(FEM_ELEM_MATRIX_ARRAY *elem_matrix)
{
	RC_TRY( bbqsort(elem_matrix->array, elem_matrix->size,
	                sizeof(FEM_ELEM_MATRIX),
	                compar_elem_matrix_element) );
	elem_matrix->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC sort_fem_elem_lookup_array(FEM_ELEM_LOOKUP_ARRAY *elem_lookup)
{
	RC_TRY( bbqsort(elem_lookup->array, elem_lookup->size,
	                sizeof(FEM_ELEM_LOOKUP), compar_elem_lookup_node) );
	elem_lookup->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC sort_fem_node_matrix_array(FEM_NODE_MATRIX_ARRAY *nmat)
{
	RC_TRY( bbqsort(nmat->array, nmat->size, sizeof(FEM_NODE_MATRIX),
	                compar_node_matrix_node) );
	nmat->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC sort_fem_load_curve_array(FEM_LOAD_CURVE_ARRAY *curve)
{
	RC_TRY( bbqsort(curve->array, curve->size,
	                sizeof(FEM_LOAD_CURVE),
	                compar_load_curve_label) );
	curve->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC sort_fem_local_coord_array(FEM_LOCAL_COORD_ARRAY *coord)
{
	RC_TRY( bbqsort(coord->array, coord->size,
	                sizeof(FEM_LOCAL_COORD),
	                compar_local_coord_label) );
	coord->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


/* for bbqsort() */
static int compar_node_label(const void *p1, const void *p2)
{
	return(compar_label( ((FEM_NODE *)p1)->label,
	                     ((FEM_NODE *)p2)->label ));
}

static int compar_node_renum(const void *p1, const void *p2)
{
	return(compar_label( ((FEM_NODE *)p1)->renum,
	                     ((FEM_NODE *)p2)->renum ));
}

static int compar_element_label(const void *p1, const void *p2)
{
	return(compar_label( ((FEM_ELEMENT *)p1)->label,
	                     ((FEM_ELEMENT *)p2)->label ));
}

static int compar_rigid_element_label(const void *p1, const void *p2)
{
	return(compar_label( ((FEM_RIGID_ELEMENT *)p1)->label,
	                     ((FEM_RIGID_ELEMENT *)p2)->label ));
}

static int compar_material_prop_label(const void *p1, const void *p2)
{
	return(compar_label( ((FEM_MATERIAL_PROP *)p1)->label,
	                     ((FEM_MATERIAL_PROP *)p2)->label ));
}

static int compar_physical_prop_label(const void *p1, const void *p2)
{
	return(compar_label( ((FEM_PHYSICAL_PROP *)p1)->label,
	                     ((FEM_PHYSICAL_PROP *)p2)->label ));
}

static int compar_bc_node(const void *p1, const void *p2)
{
	return(compar_label( ((FEM_BC *)p1)->node,
	                     ((FEM_BC *)p2)->node ));
}

static int compar_elem_bc_element(const void *p1, const void *p2)
{
	return(compar_label( ((FEM_ELEM_BC *)p1)->element,
	                     ((FEM_ELEM_BC *)p2)->element ));
}

static int compar_disp_velo_node(const void *p1, const void *p2)
{
	return(compar_label( ((FEM_DISP_VELO *)p1)->node,
	                     ((FEM_DISP_VELO *)p2)->node ));
}

static int compar_stress_strain_element(const void *p1, const void *p2)
{
	return(compar_label( ((FEM_STRESS_STRAIN *)p1)->element,
	                     ((FEM_STRESS_STRAIN *)p2)->element ));
}

static int compar_sensitivity_element(const void *p1, const void *p2)
{
	return(compar_label( ((FEM_SENSITIVITY *)p1)->element,
	                     ((FEM_SENSITIVITY *)p2)->element ));
}

static int compar_default_bc_set_id(const void *p1, const void *p2)
{
	return(compar_label( ((FEM_DEFAULT_BC *)p1)->set_id,
	                     ((FEM_DEFAULT_BC *)p2)->set_id ));
}

static int compar_sound_pressure_node(const void *p1, const void *p2)
{
	return(compar_label( ((FEM_SOUND_PRESSURE *)p1)->node,
	                     ((FEM_SOUND_PRESSURE *)p2)->node ));
}

static int compar_bc_set_label(const void *p1, const void *p2)
{
	return(compar_label( ((FEM_BC_SET *)p1)->label,
	                     ((FEM_BC_SET *)p2)->label ));
}

static int compar_elem_matrix_element(const void *p1, const void *p2)
{
	return(compar_label( ((FEM_ELEM_MATRIX *)p1)->element,
	                     ((FEM_ELEM_MATRIX *)p2)->element ));
}

static int compar_elem_lookup_node(const void *p1, const void *p2)
{
	return(compar_label( ((FEM_ELEM_LOOKUP *)p1)->node,
	                     ((FEM_ELEM_LOOKUP *)p2)->node ));
}

static int compar_node_matrix_node(const void *p1, const void *p2)
{
	return(compar_label( ((FEM_NODE_MATRIX *)p1)->node,
	                     ((FEM_NODE_MATRIX *)p2)->node ));
}

static int compar_load_curve_label(const void *p1, const void *p2)
{
	return(compar_label( ((FEM_LOAD_CURVE *)p1)->label,
	                     ((FEM_LOAD_CURVE *)p2)->label ));
}

static int compar_local_coord_label(const void *p1, const void *p2)
{
	return(compar_label( ((FEM_LOCAL_COORD *)p1)->label,
	                     ((FEM_LOCAL_COORD *)p2)->label ));
}


/* 探査プログラム */
/* ポインタ base から要素サイズ size、要素数 nmemb の配列から、 */
/* target と同じ属性の要素を探し、そのインデックスを返却 */
/* 要素同士の比較は compar()を使用 */
/* 探査に失敗すれば、-1 を返却 */
static int search(const void *base, size_t nmemb, size_t size,
                  int (*compar)(const void *, const void *),
                  const void *target, FEM_ARRAY_SORT_FLAG sort_flag)
{
	int found_index = -1;

	if((base == NULL)||(target == NULL)||(nmemb <= 0)) return(-2);

	if(sort_flag == ARRAY_SORTED){
		/* ソートされている場合には二分探査を行う */
		int low = 0;
		int high = nmemb;
		int middle, comp;

		while(low <= high){
			middle = (low + high)/2;
			comp = compar((void *)((char *)base + middle*size), target);

			if(comp == 0){
				found_index = middle;
				break;
			}else if(comp > 0){
				high = middle - 1;
			}else{
				low = middle + 1;
			}
		}
	}else{
		/* 線形探査 */
		int ii1;

		for(ii1=0;ii1<(int)nmemb;ii1++){
			if(compar((void *)((char *)base + ii1*size), target) == 0){
				found_index = ii1;
				break;
			}
		}
	}

	return(found_index);
}


/* 探査プログラムは label と index が一致している場合を想定して、 */
/* まず該当する index を検査し、一致していない場合に限って、search() を */
/* 使用する。 */
/* 同じ label が連続している場合はできる限り小さい index を返却する */

#define FUNC_SEARCH_ARRAY(name, array_type, type, elem, init_f, compar_f)  \
int name(array_type arg, int label)\
{\
	int index0, index1, index2, index3, itmp;\
	int found_index = -1;\
	type dum;\
\
	if(arg.size <= 0) return(-2);\
	index0 = label - MAX2(0, arg.array[0].elem);\
	index1 = MIN_MAX(0, index0, arg.size-1);\
	itmp = arg.array[arg.size-1].elem - arg.size + 1;\
	index2 = label - MAX2(0, itmp);\
	index3 = MIN_MAX(0, index2, arg.size-1);\
	if(compar_label(arg.array[index1].elem, label) == 0){\
		found_index = index1;\
	}else if(compar_label(arg.array[index3].elem, label) == 0){\
		found_index = index3;\
	}else{\
		init_f(&dum);\
		dum.elem = label;\
		found_index = search(arg.array, arg.size, sizeof(type),\
		                     compar_f, &dum, arg.sort_flag);\
	}\
	for(; found_index>0; found_index--){\
		if(compar_label(arg.array[found_index - 1].elem, label) != 0) break;\
	}\
	return(found_index);\
}


FUNC_SEARCH_ARRAY(search_fem_node_label, FEM_NODE_ARRAY, FEM_NODE, label,
                  init_fem_node, compar_node_label)

FUNC_SEARCH_ARRAY(search_fem_element_label, FEM_ELEMENT_ARRAY, FEM_ELEMENT,
                  label, init_fem_element, compar_element_label)

FUNC_SEARCH_ARRAY(search_fem_rigid_element_label, FEM_RIGID_ELEMENT_ARRAY,
                  FEM_RIGID_ELEMENT, label, init_fem_rigid_element,
                  compar_rigid_element_label)

FUNC_SEARCH_ARRAY(search_fem_material_prop_label, FEM_MATERIAL_PROP_ARRAY,
                  FEM_MATERIAL_PROP, label, init_fem_material_prop,
                  compar_material_prop_label)

FUNC_SEARCH_ARRAY(search_fem_physical_prop_label, FEM_PHYSICAL_PROP_ARRAY,
                  FEM_PHYSICAL_PROP, label, init_fem_physical_prop,
                  compar_physical_prop_label)

FUNC_SEARCH_ARRAY(search_fem_bc_node, FEM_BC_ARRAY, FEM_BC, node,
                  init_fem_bc, compar_bc_node)

FUNC_SEARCH_ARRAY(search_fem_elem_bc_element, FEM_ELEM_BC_ARRAY, FEM_ELEM_BC,
                  element, init_fem_elem_bc, compar_elem_bc_element)

FUNC_SEARCH_ARRAY(search_fem_disp_node, FEM_DISP_ARRAY, FEM_DISP_VELO,
                  node, init_fem_disp_velo, compar_disp_velo_node)

FUNC_SEARCH_ARRAY(search_fem_velo_node, FEM_VELO_ARRAY, FEM_DISP_VELO,
                  node, init_fem_disp_velo, compar_disp_velo_node)

FUNC_SEARCH_ARRAY(search_fem_react_node, FEM_REACT_ARRAY, FEM_DISP_VELO,
                  node, init_fem_disp_velo, compar_disp_velo_node)

FUNC_SEARCH_ARRAY(search_fem_stress_element, FEM_STRESS_ARRAY,
                  FEM_STRESS_STRAIN, element, init_fem_stress_strain,
                  compar_stress_strain_element)

FUNC_SEARCH_ARRAY(search_fem_strain_element, FEM_STRAIN_ARRAY,
                  FEM_STRESS_STRAIN, element, init_fem_stress_strain,
                  compar_stress_strain_element)

FUNC_SEARCH_ARRAY(search_fem_sensitivity_element, FEM_SENSITIVITY_ARRAY,
                  FEM_SENSITIVITY, element, init_fem_sensitivity,
                  compar_sensitivity_element)

FUNC_SEARCH_ARRAY(search_fem_default_bc_set_id, FEM_DEFAULT_BC_ARRAY,
                  FEM_DEFAULT_BC, set_id, init_fem_default_bc,
                  compar_default_bc_set_id)

FUNC_SEARCH_ARRAY(search_fem_sound_pressure_node, FEM_SOUND_PRESSURE_ARRAY,
                  FEM_SOUND_PRESSURE, node, init_fem_sound_pressure,
                  compar_sound_pressure_node)

FUNC_SEARCH_ARRAY(search_fem_bc_set_label, FEM_BC_SET_ARRAY,
                  FEM_BC_SET, label, init_fem_bc_set,
                  compar_bc_set_label)

FUNC_SEARCH_ARRAY(search_fem_elem_matrix_element, FEM_ELEM_MATRIX_ARRAY,
                  FEM_ELEM_MATRIX, element, init_fem_elem_matrix,
                  compar_elem_matrix_element)

FUNC_SEARCH_ARRAY(search_fem_elem_lookup_node, FEM_ELEM_LOOKUP_ARRAY,
                  FEM_ELEM_LOOKUP, node, init_fem_elem_lookup,
                  compar_elem_lookup_node)

FUNC_SEARCH_ARRAY(search_fem_node_matrix_node, FEM_NODE_MATRIX_ARRAY,
                  FEM_NODE_MATRIX, node, init_fem_node_matrix,
                  compar_node_matrix_node)

FUNC_SEARCH_ARRAY(search_fem_load_curve_label, FEM_LOAD_CURVE_ARRAY,
                  FEM_LOAD_CURVE, label, init_fem_load_curve,
                  compar_load_curve_label)

FUNC_SEARCH_ARRAY(search_fem_local_coord_label, FEM_LOCAL_COORD_ARRAY,
                  FEM_LOCAL_COORD, label, init_fem_local_coord,
                  compar_local_coord_label)


int search_fem_stress_element_node(FEM_STRESS_ARRAY stress, int element,
                                   int node)
{
	int found_index;

	found_index = search_fem_stress_element(stress, element);
	if(found_index < 0) return(found_index);
	
	for(; found_index<stress.size; found_index++){
		if(compar_label(stress.array[found_index].element, element) == 0){
			if(compar_label(stress.array[found_index].node, node) == 0){
				return(found_index);
			}
		}else{
			if(stress.sort_flag == ARRAY_SORTED) return(-3);
		}
	}

	return(-4);
}


int search_fem_strain_element_node(FEM_STRAIN_ARRAY strain, int element,
                                   int node)
{
	int found_index;

	found_index = search_fem_strain_element(strain, element);
	if(found_index < 0) return(found_index);
	
	for(; found_index<strain.size; found_index++){
		if(compar_label(strain.array[found_index].element, element) == 0){
			if(compar_label(strain.array[found_index].node, node) == 0){
				return(found_index);
			}
		}else{
			if(strain.sort_flag == ARRAY_SORTED) return(-3);
		}
	}

	return(-4);
}


int search_fem_sensitivity_element_node(FEM_SENSITIVITY_ARRAY sens, int element,
                                        int node)
{
	int found_index;

	found_index = search_fem_sensitivity_element(sens, element);
	if(found_index < 0) return(found_index);
	
	for(; found_index<sens.size; found_index++){
		if(compar_label(sens.array[found_index].element, element) == 0){
			if(compar_label(sens.array[found_index].node, node) == 0){
				return(found_index);
			}
		}else{
			if(sens.sort_flag == ARRAY_SORTED) return(-3);
		}
	}

	return(-4);
}


/* disp の weight 倍を disp_total に付加 */
/* disp_total に存在しない要素は新規に追加する */
RC add_fem_disp_array(double weight, FEM_DISP_ARRAY disp,
                      FEM_DISP_ARRAY *disp_total)
{
	int ii1, index, num;
	TRANS_ROTATION3D tmp;

	for(ii1=0; ii1<disp.size; ii1++){
		if(disp.array[ii1].node < 0) continue;

		index = search_fem_disp_node(*disp_total, disp.array[ii1].node);
		if(index >= 0){
			tmp = mul_scalar_trans_rotation3d(weight, disp.array[ii1].v);
			disp_total->array[index].v = add_trans_rotation3d(tmp,
			                               disp_total->array[index].v);
			disp_total->array[index].s += weight * disp.array[ii1].s;
		}else{
			RC_TRY( realloc_fem_disp_array(disp_total) );
			num = disp_total->size;
			disp_total->array[num] = disp.array[ii1];
			disp_total->array[num].v = mul_scalar_trans_rotation3d(weight,
			                                            disp.array[ii1].v);
			disp_total->array[num].s = weight * disp.array[ii1].s;
			(disp_total->size)++;
			RC_TRY( clean_fem_disp_array(disp_total) );
		}
	}

	return(NORMAL_RC);
}


/* sens の weight 倍を sens_total に付加 */
/* sens_total に存在しない要素は新規に追加する */
RC add_fem_sensitivity_array(double weight, FEM_SENSITIVITY_ARRAY sens,
                             FEM_SENSITIVITY_ARRAY *sens_total)
{
	int ii1, index, num;
	TRANS_ROTATION3D tmp;

	for(ii1=0; ii1<sens.size; ii1++){
		if(sens.array[ii1].element < 0) continue;

		index = search_fem_sensitivity_element_node(*sens_total,
		                  sens.array[ii1].element, sens.array[ii1].node);
		if(index >= 0){
			tmp = mul_scalar_trans_rotation3d(weight, sens.array[ii1].v);
			sens_total->array[index].v = add_trans_rotation3d(tmp,
			                               sens_total->array[index].v);
			sens_total->array[index].s += weight * sens.array[ii1].s;
		}else{
			RC_TRY( realloc_fem_sensitivity_array(sens_total) );
			num = sens_total->size;
			sens_total->array[num] = sens.array[ii1];
			sens_total->array[num].v = mul_scalar_trans_rotation3d(weight,
			                                            sens.array[ii1].v);
			sens_total->array[num].s = weight * sens.array[ii1].s;
			(sens_total->size)++;
			RC_TRY( clean_fem_sensitivity_array(sens_total) );
		}
	}

	return(NORMAL_RC);
}


/* 要素 element に対応する材料 mat のインデックスを返却 */
int search_fem_material_prop_element(FEM_MATERIAL_PROP_ARRAY mat,
                                     FEM_PHYSICAL_PROP_ARRAY phys,
                                     FEM_ELEMENT *element)
{
	int phys_index;

	if((phys.source == FEM_NASTRAN)&&(element->material < 0)){
		/* FEM_NASTRAN の場合は、一度 phys を検索する必要有り */
		phys_index = search_fem_physical_prop_label(phys, element->physical);
		if(phys_index < 0){
			return(phys_index);
		}

		/* 次回の検索を高速化する */
		element->material = phys.array[phys_index].i_info[0];
	}

	if(element->material < 0) return(-2);

	return(search_fem_material_prop_label(mat, element->material));
}


/* 要素 element を構成する節点座標値のマトリックスを作成 */
RC node_location_matrix(FEM_ELEMENT *element,
                        FEM_NODE_ARRAY node, double **matrix)
{
	int ii1;

	for(ii1=0; ii1<element->node_num; ii1++){
		int index = search_fem_node_label(node, element->node[ii1]);
		if(index < 0){
			fprintf(stderr, "search_fem_node_label < \n");
			return(SEEK_ERROR_RC);
		}

		matrix[ii1][0] = node.array[index].p.x;
		matrix[ii1][1] = node.array[index].p.y;
		matrix[ii1][2] = node.array[index].p.z;
	}

	return(NORMAL_RC);
}


RC fill_disp_matrix(FEM_ELEMENT element,
                    FEM_DISP_ARRAY disp, double **matrix)
{
	int ii1;

	for(ii1=0; ii1<element.node_num; ii1++){
		int index;

		RC_NEG_CHK( index = search_fem_disp_node(disp, element.node[ii1]) );

		matrix[ii1][0] = disp.array[index].v.x;
		matrix[ii1][1] = disp.array[index].v.y;
		matrix[ii1][2] = disp.array[index].v.z;
		matrix[ii1][3] = disp.array[index].v.yz;
		matrix[ii1][4] = disp.array[index].v.zx;
		matrix[ii1][5] = disp.array[index].v.xy;
	}

	return(NORMAL_RC);
}


RC fill_sound_pres_matrix(FEM_ELEMENT element,
                          FEM_SOUND_PRESSURE_ARRAY pres, double **matrix)
{
	int ii1;

	for(ii1=0; ii1<element.node_num; ii1++){
		int index;

		RC_NEG_CHK( index = search_fem_sound_pressure_node(pres,
		                                    element.node[ii1]) );
		matrix[ii1][0] = pres.array[index].re;
		matrix[ii1][1] = pres.array[index].im;
	}

	return(NORMAL_RC);
}


RC unit_normal_vect_center(FEM_ELEMENT surf, FEM_ELEMENT inner_elem,
                           FEM_NODE_ARRAY node, TRANSLATION3D *vect)
{
	double center_coord[4];

	RC_TRY( set_center_point(surf.type, center_coord) );
	RC_TRY( unit_normal_vect(center_coord, surf, inner_elem, node, vect) );

	return(NORMAL_RC);
}


/* 要素 surf 上の座標 local_coord における外向き単位法線ベクトル */
/* inner_elem の反対向きを外向きとする */
RC unit_normal_vect(double *local_coord, FEM_ELEMENT surf,
                    FEM_ELEMENT inner_elem, FEM_NODE_ARRAY node,
                    TRANSLATION3D *vect)
{
	int ii1;
	int index;
	double s_coord[4];
	double i_coord[4];
	double s_N[FEM_MAX_NODE];
	double i_N[FEM_MAX_NODE];
	TRANSLATION3D s_center;
	TRANSLATION3D i_center;
	TRANSLATION3D d_xi, d_eta;
	TRANSLATION3D s_to_i;
	static int init_flag = 0;
	static double **dN = NULL;
	static double **node_location = NULL;
	static double **jacobi = NULL;
	double abs_vect;

	if(init_flag == 0){
		RC_TRY( allocate2D(3, FEM_MAX_NODE, &dN) );
		RC_TRY( allocate2D(FEM_MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(3, 3, &jacobi) );

		init_flag = 1;
	}

	/* 法線ベクトルの計算 */
	init_matrix(3, FEM_MAX_NODE, dN);   /* 初期化が必要 */
	RC_TRY( node_location_matrix(&surf, node, node_location) );
	RC_TRY( set_dN_matrix(surf.type, local_coord, dN) );
	mul_matrix(3, surf.node_num, 3, dN, node_location, jacobi);
	switch(surf.type){
	case ELEM_TRI1:
	case ELEM_TRI2:
	case ELEM_QUAD1:
	case ELEM_QUAD2:
		d_xi.x = jacobi[0][0];
		d_xi.y = jacobi[0][1];
		d_xi.z = jacobi[0][2];
		d_eta.x = jacobi[1][0];
		d_eta.y = jacobi[1][1];
		d_eta.z = jacobi[1][2];
		*vect = outer_product3d(d_xi, d_eta);
		break;
	case ELEM_RBAR:
	case ELEM_BEAM1:
	case ELEM_LINE1:
	case ELEM_BEAM2:
	case ELEM_LINE2:
		vect->x = -jacobi[0][1];
		vect->y = jacobi[0][0];
		vect->z = 0.0;
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	/* 法線ベクトルを単位化 */
	abs_vect = abs_translation3d(*vect);
	if(nearly_eq(abs_vect, 0.0)){
		fprintf(stderr, "abs_translation3d < ");
		return(CAL_ERROR_RC);
	}
	vect->x /= abs_vect;
	vect->y /= abs_vect;
	vect->z /= abs_vect;

	RC_TRY( set_center_point(surf.type, s_coord) );
	RC_TRY( set_N_vector(surf.type, s_coord, s_N) );
	s_center.x = 0.0;
	s_center.y = 0.0;
	s_center.z = 0.0;
	for(ii1=0; ii1<surf.node_num; ii1++){
		index = search_fem_node_label(node, surf.node[ii1]);
		if(index < 0){
			fprintf(stderr, "search_fem_node_label < ");
			return(SEEK_ERROR_RC);
		}
		s_center.x += s_N[ii1] * node.array[index].p.x;
		s_center.y += s_N[ii1] * node.array[index].p.y;
		s_center.z += s_N[ii1] * node.array[index].p.z;
	}

	RC_TRY( set_center_point(inner_elem.type, i_coord) );
	RC_TRY( set_N_vector(inner_elem.type, i_coord, i_N) );
	i_center.x = 0.0;
	i_center.y = 0.0;
	i_center.z = 0.0;
	for(ii1=0; ii1<inner_elem.node_num; ii1++){
		index = search_fem_node_label(node, inner_elem.node[ii1]);
		if(index < 0){
			fprintf(stderr, "search_fem_node_label < ");
			return(SEEK_ERROR_RC);
		}
		i_center.x += i_N[ii1] * node.array[index].p.x;
		i_center.y += i_N[ii1] * node.array[index].p.y;
		i_center.z += i_N[ii1] * node.array[index].p.z;
	}

	/* 法線ベクトルの向きをチェック */
	s_to_i = sub_translation3d(i_center, s_center);
	if(nearly_eq(abs_translation3d(s_to_i), 0.0)){
		fprintf(stderr, "abs_translation3d < ");
		return(CAL_ERROR_RC);
	}
	if(inner_product3d(*vect, s_to_i) > 0.0){
		vect->x *= -1.0;
		vect->y *= -1.0;
		vect->z *= -1.0;
	}

	return(NORMAL_RC);
}


RC extract_rigid_set(int num, const int id[], double scale[],
                     FEM_RIGID_ELEMENT_ARRAY src, FEM_RIGID_ELEMENT_ARRAY *dest)
{
	int ii1, ii2;
	int set_flag;

	RC_TRY( allocate_fem_rigid_element_array(0, dest) );
	dest->source = src.source;
	
	for(ii1=0; ii1<src.size; ii1++){
		if(src.array[ii1].label < 0) continue;

		set_flag = 0;
		if(src.array[ii1].set_id < 0) set_flag = 1;
		for(ii2=0; ii2<num; ii2++){
			if(src.array[ii1].set_id == id[ii2]) set_flag = 1;
		}

		if(set_flag){
			RC_TRY( realloc_fem_rigid_element_array(dest) );
			RC_TRY( copy_fem_rigid_element(src.array[ii1],
			                               &(dest->array[dest->size])) );
			(dest->size)++;
		}
	}
	RC_TRY( clean_fem_rigid_element_array(dest) );

	return(NORMAL_RC);
}


/* 境界条件 src 内より、id[num] に
   該当する要素を scale[] 倍して dest にコピー */
RC extract_bc_set(int num, const int id[], double scale[],
                  FEM_BC_ARRAY src, FEM_BC_ARRAY *dest)
{
	int ii1, ii2, ii3;

	RC_TRY( allocate_fem_bc_array(0, dest) );
	dest->source = src.source;
	
	for(ii1=0; ii1<src.size; ii1++){
		if(src.array[ii1].node < 0) continue;

		for(ii2=0; ii2<num; ii2++){
			if(src.array[ii1].set_id == id[ii2]){
				RC_TRY( realloc_fem_bc_array(dest) );
				dest->array[dest->size] = src.array[ii1];
				dest->array[dest->size].v
					= mul_scalar_trans_rotation3d(scale[ii2],
					               dest->array[dest->size].v);
				for(ii3=0; ii3<FEM_BC_SCALAR; ii3++){
					dest->array[dest->size].s[ii3] *= scale[ii2];
				}
				(dest->size)++;
				break;
			}
		}
	}
	RC_TRY( compress_fem_bc_array(dest) );
	RC_TRY( clean_fem_bc_array(dest) );

	return(NORMAL_RC);
}


RC compress_fem_bc_array(FEM_BC_ARRAY *dest)
{
	int ii1, ii2, ii3;

	/* 重複する節点をマージ */
	for(ii1=0; ii1<dest->size;){
		if(dest->array[ii1].node < 0) continue;
		for(ii2=ii1+1; ii2<dest->size; ii2++){
			if(dest->array[ii1].node != dest->array[ii2].node) break;
			RC_TRY( overwrite_bc(
			           &(dest->array[ii1].v_type.x), &(dest->array[ii1].v.x),
			                dest->array[ii2].v_type.x, dest->array[ii2].v.x) );
			RC_TRY( overwrite_bc(
			           &(dest->array[ii1].v_type.y), &(dest->array[ii1].v.y),
			                dest->array[ii2].v_type.y, dest->array[ii2].v.y) );
			RC_TRY( overwrite_bc(
			           &(dest->array[ii1].v_type.z), &(dest->array[ii1].v.z),
			                dest->array[ii2].v_type.z, dest->array[ii2].v.z) );
			RC_TRY( overwrite_bc(
			         &(dest->array[ii1].v_type.yz), &(dest->array[ii1].v.yz),
			              dest->array[ii2].v_type.yz, dest->array[ii2].v.yz) );
			RC_TRY( overwrite_bc(
			         &(dest->array[ii1].v_type.zx), &(dest->array[ii1].v.zx),
			              dest->array[ii2].v_type.zx, dest->array[ii2].v.zx) );
			RC_TRY( overwrite_bc(
			         &(dest->array[ii1].v_type.xy), &(dest->array[ii1].v.xy),
			              dest->array[ii2].v_type.xy, dest->array[ii2].v.xy) );
			for(ii3=0; ii3<FEM_BC_SCALAR; ii3++){
				RC_TRY( overwrite_bc(&(dest->array[ii1].s_type[ii3]),
				                     &(dest->array[ii1].s[ii3]),
				                     dest->array[ii2].s_type[ii3],
				                     dest->array[ii2].s[ii3]) );
			}
			dest->array[ii2].node = -1;
		}
		ii1 = ii2;
	}

	/* 完全に Free の節点を削除 */
	for(ii1=0; ii1<dest->size; ii1++){
		int flag = 0;

		if(dest->array[ii1].node < 0) continue;
		if( (dest->array[ii1].v_type.x != BC_FREE)
		  ||(dest->array[ii1].v_type.y != BC_FREE)
		  ||(dest->array[ii1].v_type.z != BC_FREE)
		  ||(dest->array[ii1].v_type.yz != BC_FREE)
		  ||(dest->array[ii1].v_type.zx != BC_FREE)
		  ||(dest->array[ii1].v_type.xy != BC_FREE) ){
			flag = 1;
		}
		for(ii2=0; ii2<FEM_BC_SCALAR; ii2++){
			if(dest->array[ii1].s_type[ii2] != BC_FREE) flag = 1;
		}
		if(flag == 0){
			dest->array[ii1].node = -1;
		}
	}

	RC_TRY( clean_fem_bc_array(dest) );

	return(NORMAL_RC);
}


static RC overwrite_bc(FEM_BC_TYPE *type1, double *v1,
                       FEM_BC_TYPE type2,  double v2)
{
	/* 同じ境界条件なら境界値を加算 */
	if(*type1 == type2){
		if(fabs(*v1) < fabs(v2)) *v1 += v2;
		return(NORMAL_RC);
	}

	/* BC_FREE なら上書き */
	if( (*type1 == BC_FREE)&&(type2 != BC_FREE) ){
		*type1 = type2;
		*v1 = v2;
		return(NORMAL_RC);
	}

	/* BC_FORCE(_NL) より BC_FIX(_NL) が優先 */
	if( ((*type1 == BC_FORCE)&&(type2 == BC_FIX))
	  ||((*type1 == BC_FORCE_NL)&&(type2 == BC_FIX_NL)) ){
		*type1 = type2;
		*v1 = v2;
		return(NORMAL_RC);
	}

	return(NORMAL_RC);
}


RC extract_default_bc_set(int num, const int id[], double scale[],
                          FEM_DEFAULT_BC_ARRAY src,
                          FEM_DEFAULT_BC_ARRAY *dest)
{
	int ii1, ii2;

	RC_TRY( allocate_fem_default_bc_array(0, dest) );
	dest->source = src.source;
	
	for(ii1=0; ii1<src.size; ii1++){
		if(src.array[ii1].set_id < 0) continue;

		for(ii2=0; ii2<num; ii2++){
			if(src.array[ii1].set_id == id[ii2]){
				RC_TRY( realloc_fem_default_bc_array(dest) );
				dest->array[dest->size] = src.array[ii1];
				dest->array[dest->size].temp *= scale[ii2];
				dest->array[dest->size].grav
					= mul_scalar_translation3d(scale[ii2],
					               dest->array[dest->size].grav);
				(dest->size)++;
				break;
			}
		}
	}

	RC_TRY( clean_fem_default_bc_array(dest) );

	return(NORMAL_RC);
}


RC extract_elem_bc_set(int num, const int id[], double scale[],
                       FEM_ELEM_BC_ARRAY src, FEM_ELEM_BC_ARRAY *dest)
{
	int ii1, ii2, ii3;

	RC_TRY( allocate_fem_elem_bc_array(0, dest) );
	dest->source = src.source;
	
	for(ii1=0; ii1<src.size; ii1++){
		if(src.array[ii1].element < 0) continue;

		for(ii2=0; ii2<num; ii2++){
			if(src.array[ii1].set_id == id[ii2]){
				RC_TRY( realloc_fem_elem_bc_array(dest) );
				dest->array[dest->size] = src.array[ii1];
				if( (dest->array[dest->size].type == BC_PRESS)
				  ||(dest->array[dest->size].type == BC_TRACTION) ){
					dest->array[dest->size].val.s *= scale[ii2];
				}else if( (dest->array[dest->size].type == BC_PRESS_D1)
				        ||(dest->array[dest->size].type == BC_TRACTION_D1)
				        ||(dest->array[dest->size].type == BC_PRESS_D2)
				        ||(dest->array[dest->size].type == BC_TRACTION_D2) ){
					for(ii3=0; ii3<FEM_MAX_NODE; ii3++){
						dest->array[dest->size].val.ns[ii3] *= scale[ii2];
					}
				}else{
					return(IMPLEMENT_ERROR_RC);
				}
				(dest->size)++;
				break;
			}
		}
	}

	RC_TRY( clean_fem_elem_bc_array(dest) );

	return(NORMAL_RC);
}


/* elem_matrixより作成される全体マトリックスと vector[total_dof] の掛け算 */
/* 結果は answer[total_dof] に格納される */
/* 可能な限りエラーチェックを入れた低速版 */
RC mul_elem_matrices_vector(int total_dof, FEM_ELEM_MATRIX_ARRAY elem_matrix,
                            const double *vector, double *answer)
{
	int ii1, ii2, ii3;

	if(total_dof <= 0) return(ARG_ERROR_RC);

	for(ii1=0; ii1<total_dof; ii1++){
		answer[ii1] = 0.0;
	}

	for(ii1=0; ii1<elem_matrix.size; ii1++){
		if(elem_matrix.array[ii1].element < 0) continue;
		if( (elem_matrix.array[ii1].index == NULL)
		  ||(elem_matrix.array[ii1].matrix == NULL) ){
			return(ARG_ERROR_RC);
		}

		for(ii2=0; ii2<elem_matrix.array[ii1].size; ii2++){
			if( (elem_matrix.array[ii1].index[ii2] < 0)
			  ||(elem_matrix.array[ii1].index[ii2] >= total_dof) ){
				return(SEEK_ERROR_RC);
			}
		}
		for(ii2=0; ii2<elem_matrix.array[ii1].size; ii2++){
			for(ii3=0; ii3<elem_matrix.array[ii1].size; ii3++){
				answer[elem_matrix.array[ii1].index[ii2]]
				      += vector[elem_matrix.array[ii1].index[ii3]]
				       * elem_matrix.array[ii1].matrix[ii2][ii3];
			}
		}
	}

	return(NORMAL_RC);
}


/* elem_matrix の m 番目の解が v となるように vectors[vec_size][] を修正 */
/* vectors[][m] への代入は 上位関数が行う必要あり */
RC modify_elem_matrix_n(FEM_ELEM_MATRIX elem_matrix, int m, double v,
                        double **vectors, int vec_size)
{
	int ii1, ii2;

	if((m < 0)||(elem_matrix.size < m)) return(ARG_ERROR_RC);

	for(ii1=0; ii1<elem_matrix.size; ii1++){
		if(ii1 == m){
			elem_matrix.matrix[m][ii1] = 1.0;
		}else{
			elem_matrix.matrix[m][ii1] = 0.0;
			for(ii2=0; ii2<vec_size; ii2++){
				vectors[ii2][elem_matrix.index[ii1]]
				          -= elem_matrix.matrix[ii1][m] * v;
			}
			elem_matrix.matrix[ii1][m] = 0.0;
		}
	}

	return(NORMAL_RC);
}


RC make_elem_lookup(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                    FEM_ELEM_LOOKUP_ARRAY *elem_lookup)
{
	int ii1, ii2;

	/* elem_lookup の確保と初期化 */
	RC_TRY( allocate_fem_elem_lookup_array(node.size, elem_lookup) );
	for(ii1=0; ii1<node.size; ii1++){
		elem_lookup->array[ii1].node = node.array[ii1].label;
		elem_lookup->array[ii1].size = 0;
		elem_lookup->array[ii1].alloc_size = 0;
		elem_lookup->array[ii1].r_size = 0;
		elem_lookup->array[ii1].r_alloc_size = 0;
	}
	elem_lookup->size = node.size;
	RC_TRY( clean_fem_elem_lookup_array(elem_lookup) );

	/* 各テーブルのサイズを決定 */
	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
			int index = search_fem_elem_lookup_node(*elem_lookup,
			                          element.array[ii1].node[ii2]);
			RC_NEG_CHK(index);
			elem_lookup->array[index].alloc_size++;
		}
	}

	/* 各テーブルを確保 */
	for(ii1=0; ii1<node.size; ii1++){
		if(elem_lookup->array[ii1].node < 0) continue;
		if(elem_lookup->array[ii1].alloc_size > 0){
			elem_lookup->array[ii1].element
			     = malloc(sizeof(int)*(elem_lookup->array[ii1].alloc_size));
			RC_NULL_CHK( elem_lookup->array[ii1].element );
		}
	}

	/* 各テーブルに代入 */
	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
			int size;
			int index = search_fem_elem_lookup_node(*elem_lookup,
			                          element.array[ii1].node[ii2]);
			RC_NEG_CHK(index);

			size = elem_lookup->array[index].size;
			if((size < 0)||(size >= elem_lookup->array[index].alloc_size)){
				return(UNKNOWN_ERROR_RC);
			}
			elem_lookup->array[index].element[size] = element.array[ii1].label;
			elem_lookup->array[index].size++;
		}
	}

	return(NORMAL_RC);
}


/* rigid_elem も考慮して elem_lookup を構築 */
RC make_elem_lookup_rigid(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                          FEM_RIGID_ELEMENT_ARRAY rigid_elem,
                          FEM_ELEM_LOOKUP_ARRAY *elem_lookup)
{
	int ii1, ii2;

	/* elem_lookup の確保と初期化，element分の処理 */
	RC_TRY( make_elem_lookup(node, element, elem_lookup) );

	/* 各テーブルのサイズを決定 */
	for(ii1=0; ii1<rigid_elem.size; ii1++){
		if(rigid_elem.array[ii1].label < 0) continue;

		for(ii2=0; ii2<rigid_elem.array[ii1].node_num; ii2++){
			int index;
			index = search_fem_elem_lookup_node(*elem_lookup,
			                          rigid_elem.array[ii1].node[ii2]);
			if(index < 0){
				RC_TRY( realloc_fem_elem_lookup_array(elem_lookup) );
				elem_lookup->array[elem_lookup->size].node
				           = rigid_elem.array[ii1].node[ii2];
				(elem_lookup->size)++;
				RC_TRY( clean_fem_elem_lookup_array(elem_lookup) );
				index = search_fem_elem_lookup_node(*elem_lookup,
				                          rigid_elem.array[ii1].node[ii2]);
			}
			RC_NEG_CHK( index );
			elem_lookup->array[index].r_alloc_size++;
		}
	}

	/* 各テーブルを確保 */
	for(ii1=0; ii1<elem_lookup->size; ii1++){
		if(elem_lookup->array[ii1].node < 0) continue;
		if(elem_lookup->array[ii1].r_alloc_size > 0){
			elem_lookup->array[ii1].r_element
			     = malloc(sizeof(int)*(elem_lookup->array[ii1].r_alloc_size));
			RC_NULL_CHK( elem_lookup->array[ii1].r_element );
		}
	}

	/* 各テーブルに代入 */
	for(ii1=0; ii1<rigid_elem.size; ii1++){
		if(rigid_elem.array[ii1].label < 0) continue;

		for(ii2=0; ii2<rigid_elem.array[ii1].node_num; ii2++){
			int size;
			int index = search_fem_elem_lookup_node(*elem_lookup,
			                          rigid_elem.array[ii1].node[ii2]);
			RC_NEG_CHK(index);

			size = elem_lookup->array[index].r_size;
			if((size < 0)||(size >= elem_lookup->array[index].r_alloc_size)){
				return(UNKNOWN_ERROR_RC);
			}
			elem_lookup->array[index].r_element[size]
			                     = rigid_elem.array[ii1].label;
			elem_lookup->array[index].r_size++;
		}
	}

	return(NORMAL_RC);
}

/* 問題の次元 */
/* 2...２次元問題、3...３次元問題 */
int fem_analysis_dim(FEM_ELEMENT_ARRAY element)
{
	int dim, tmp_dim;
	int ii1;

	dim = 0;
	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		tmp_dim = element_dim(element.array[ii1].type);
		if(dim == 0){
			dim = tmp_dim;
		}else if(dim != tmp_dim){
			return(-1);
		}
	}
	if(dim <= 0){
		return(-1);
	}

	return(dim);
}


/* 要素の次数 */
/* 1...１次要素、2...２次要素 */
/* 混在している場合は大きい次数を返却 */
int fem_analysis_order(FEM_ELEMENT_ARRAY element)
{
	int ii1;
	int order = -1;
	int tmp_order;

	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;
		tmp_order = element_order(element.array[ii1].type);
		if(tmp_order > order) order = tmp_order;
	}
#if 0
	int order, tmp_order;
	int ii1;
	order = 0;
	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		tmp_order = element_order(element.array[ii1].type);
		if(order == 0){
			order = tmp_order;
		}else if(order != tmp_order){
			return(-1);
		}
	}
	if(order <= 0){
		return(-1);
	}
#endif

	return(order);
}


/* 全ての要素の次数を二次から一次へ変換 */
RC reduced_element_array_order(FEM_ELEMENT_ARRAY *elem)
{
	FEM_ELEM_TYPE new_type;
	int number;
	int ii1;

	for(ii1=0; ii1<elem->size; ii1++){
		if(elem->array[ii1].label < 0) continue;

		RC_TRY( reduced_element_type(elem->array[ii1].type, &new_type,
		                             &number) );
		elem->array[ii1].node_num = number;
		elem->array[ii1].type = new_type;
	}

	return(NORMAL_RC);
}


/* 1つの二次要素の type から、一次要素の type と node_num を受け取る */
RC reduced_element_type(FEM_ELEM_TYPE old_type, FEM_ELEM_TYPE *new_type,
                        int *new_node_num)
{
	switch(old_type){
	case ELEM_BEAM1:
	case ELEM_BEAM2:
		*new_type = ELEM_BEAM1;
		*new_node_num = 2;
		break;
	case ELEM_LINE1:
	case ELEM_LINE2:
		*new_type = ELEM_LINE1;
		*new_node_num = 2;
		break;
	case ELEM_TRI1:
	case ELEM_TRI2:
		*new_type = ELEM_TRI1;
		*new_node_num = 3;
		break;
	case ELEM_QUAD1:
	case ELEM_QUAD2:
		*new_type = ELEM_QUAD1;
		*new_node_num = 4;
		break;
	case ELEM_TETRA1:
	case ELEM_TETRA2:
		*new_type = ELEM_TETRA1;
		*new_node_num = 4;
		break;
	case ELEM_PENTA1:
	case ELEM_PENTA2:
		*new_type = ELEM_PENTA1;
		*new_node_num = 6;
		break;
	case ELEM_HEXA1:
	case ELEM_HEXA2:
		*new_type = ELEM_HEXA1;
		*new_node_num = 8;
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


/* 距離が cr_dist 以下の節点を重複しているとみなして削除 */
RC delete_duplicate_node(FEM_ELEMENT_ARRAY element, FEM_NODE_ARRAY *node,
                         double cr_dist)
{
	int ii1, ii2, ii3, ii4;
	FEM_ELEM_LOOKUP_ARRAY elem_lookup;
	double dist;
	int index, elem_index;

	RC_TRY( make_elem_lookup(*node, element, &elem_lookup) );

	for(ii1=0; ii1<node->size; ii1++){
		if(node->array[ii1].label < 0) continue;

		for(ii2=ii1+1; ii2<node->size; ii2++){
			if(node->array[ii2].label < 0) continue;

			dist = dist_point3d(node->array[ii1].p, node->array[ii2].p);
			if(dist >= cr_dist) continue;

			index = search_fem_elem_lookup_node(elem_lookup,
			                                    node->array[ii2].label);
			for(ii3=0; ii3<elem_lookup.array[index].size; ii3++){
				elem_index = search_fem_element_label(element,
				                  elem_lookup.array[index].element[ii3]);
				RC_NEG_CHK(elem_index);
				
				for(ii4=0; ii4<element.array[elem_index].node_num; ii4++){
					if(element.array[elem_index].node[ii4]
					             == node->array[ii2].label){
						element.array[elem_index].node[ii4]
					 	                = node->array[ii1].label;
					}
				}
			}
			node->array[ii2].label = -1;
		}
	}

	clean_fem_node_array(node);
	RC_TRY( free_fem_elem_lookup_array(&elem_lookup) );

	return(NORMAL_RC);
}


/* 節点間距離の最小値 */
RC minimum_node_distance(FEM_ELEMENT_ARRAY element, FEM_NODE_ARRAY node,
                         double *dist)
{
	int ii1, ii2, ii3;
	double dist_temp;
	int index0, index1;

	*dist = DBL_MAX;
	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;
		if(element.array[ii1].node_num <= 1) continue;

		for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
			index0 = search_fem_node_label(node,
			                element.array[ii1].node[ii2]);
			RC_NEG_CHK(index0);

			for(ii3=0; ii3<element.array[ii1].node_num; ii3++){
				if(ii2 == ii3) continue;
				index1 = search_fem_node_label(node,
				                element.array[ii1].node[ii3]);
				RC_NEG_CHK(index1);

				dist_temp = dist_point3d(node.array[index0].p,
				                         node.array[index1].p);
				if(dist_temp < *dist){
					*dist = dist_temp;
				}
			}
		}
	}

	return(NORMAL_RC);
}


/* elem.array[].node[] に使用されていない節点を消去する */
/* 消去後の node は clean_fem_node_array() を実行される */ 
RC delete_unused_node(FEM_ELEMENT_ARRAY elem, FEM_NODE_ARRAY *node)
{
	int *flags;
	int index;
	int ii1, ii2;

	RC_NULL_CHK( flags = (int *)malloc( node->size * sizeof(int) ) );
	for(ii1=0; ii1<node->size; ii1++){
		flags[ii1] = 0;
	}

	for(ii1=0; ii1<elem.size; ii1++){
		if(elem.array[ii1].label < 0) continue;
		for(ii2=0; ii2<elem.array[ii1].node_num; ii2++){
			RC_NEG_CHK( index = search_fem_node_label(*node,
			                     elem.array[ii1].node[ii2]) );
			flags[index] = 1;
		}
	}
	for(ii1=0; ii1<node->size; ii1++){
		if(flags[ii1] == 0) node->array[ii1].label = -1;
	}

	RC_TRY( clean_fem_node_array(node) );
	node->renum_flag = UNRENUMBERED;
	free( (void *)flags );

	return(NORMAL_RC);
}


RC delete_unused_node_r(FEM_ELEMENT_ARRAY elem, FEM_RIGID_ELEMENT_ARRAY rigid,
                        FEM_NODE_ARRAY *node)
{
	int *flags;
	int index;
	int ii1, ii2;

	RC_NULL_CHK( flags = (int *)malloc( node->size * sizeof(int) ) );
	for(ii1=0; ii1<node->size; ii1++){
		flags[ii1] = 0;
	}

	for(ii1=0; ii1<elem.size; ii1++){
		if(elem.array[ii1].label < 0) continue;
		for(ii2=0; ii2<elem.array[ii1].node_num; ii2++){
			RC_NEG_CHK( index = search_fem_node_label(*node,
			                     elem.array[ii1].node[ii2]) );
			flags[index] = 1;
		}
	}

	for(ii1=0; ii1<rigid.size; ii1++){
		if(rigid.array[ii1].label < 0) continue;
		for(ii2=0; ii2<rigid.array[ii1].node_num; ii2++){
			RC_NEG_CHK( index = search_fem_node_label(*node,
			                     rigid.array[ii1].node[ii2]) );
			flags[index] = 1;
		}
	}

	for(ii1=0; ii1<node->size; ii1++){
		if(flags[ii1] == 0) node->array[ii1].label = -1;
	}

	RC_TRY( clean_fem_node_array(node) );
	node->renum_flag = UNRENUMBERED;
	free( (void *)flags );

	return(NORMAL_RC);
}


RC delete_unused_bc(FEM_NODE_ARRAY node, FEM_BC_ARRAY *bc)
{
	int index;
	int ii1;

	for(ii1=0; ii1<bc->size; ii1++){
		if(bc->array[ii1].node < 0) continue;
		index = search_fem_node_label(node, bc->array[ii1].node);
		if(index < 0){
			bc->array[ii1].node = -1;
		}
	}

	RC_TRY( clean_fem_bc_array(bc) );

	return(NORMAL_RC);
}


/* 中間節点を両端の中点に移動させる */
RC locate_middle_node(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY elem)
{
	int ii1, ii2;
	int table[FEM_MAX_NODE][2];
	int n_index, index0, index1;

	for(ii1=0; ii1<elem.size; ii1++){
		if(elem.array[ii1].label < 0) continue;

		/* ２次要素以外は処理しない */
		if( element_order(elem.array[ii1].type) != 2) continue;

		RC_TRY( interpolate_table(elem.array[ii1].type, table) );
		for(ii2=0; ii2<elem.array[ii1].node_num; ii2++){
			RC_NEG_CHK( n_index = search_fem_node_label(node,
			                       elem.array[ii1].node[ii2]) );
			RC_NEG_CHK( index0 = search_fem_node_label(node,
			                       elem.array[ii1].node[table[ii2][0]]) );
			RC_NEG_CHK( index1 = search_fem_node_label(node,
			                       elem.array[ii1].node[table[ii2][1]]) );
			node.array[n_index].p = mid_point3d(node.array[index0].p,
			                                    node.array[index1].p);
		}
	}

	return(NORMAL_RC);
}


RC interpolate_table(FEM_ELEM_TYPE type, int table[][2])
{
	switch(type){
	case ELEM_BEAM2:
	case ELEM_LINE2:
		table[0][0] = 0;   table[0][1] = 0;
		table[1][0] = 1;   table[1][1] = 1;

		table[2][0] = 0;   table[2][1] = 1;
		break;
	case ELEM_TRI2:
		table[0][0] = 0;   table[0][1] = 0;
		table[1][0] = 1;   table[1][1] = 1;
		table[2][0] = 2;   table[2][1] = 2;

		table[3][0] = 1;   table[3][1] = 2;
		table[4][0] = 0;   table[4][1] = 2;
		table[5][0] = 0;   table[5][1] = 1;
		break;
	case ELEM_QUAD2:
		table[0][0] = 0;   table[0][1] = 0;
		table[1][0] = 1;   table[1][1] = 1;
		table[2][0] = 2;   table[2][1] = 2;
		table[3][0] = 3;   table[3][1] = 3;

		table[4][0] = 0;   table[4][1] = 1;
		table[5][0] = 1;   table[5][1] = 2;
		table[6][0] = 2;   table[6][1] = 3;
		table[7][0] = 0;   table[7][1] = 3;
		break;
	case ELEM_TETRA2:
		table[0][0] = 0;   table[0][1] = 0;
		table[1][0] = 1;   table[1][1] = 1;
		table[2][0] = 2;   table[2][1] = 2;
		table[3][0] = 3;   table[3][1] = 3;

		table[4][0] = 1;   table[4][1] = 2;
		table[5][0] = 0;   table[5][1] = 2;
		table[6][0] = 0;   table[6][1] = 1;
		table[7][0] = 0;   table[7][1] = 3;
		table[8][0] = 1;   table[8][1] = 3;
		table[9][0] = 2;   table[9][1] = 3;
		break;
	case ELEM_PENTA2:
		table[0][0] = 0;   table[0][1] = 0;
		table[1][0] = 1;   table[1][1] = 1;
		table[2][0] = 2;   table[2][1] = 2;
		table[3][0] = 3;   table[3][1] = 3;
		table[4][0] = 4;   table[4][1] = 4;
		table[5][0] = 5;   table[5][1] = 5;

		table[6][0] = 1;   table[6][1] = 2;
		table[7][0] = 0;   table[7][1] = 2;
		table[8][0] = 0;   table[8][1] = 1;
		table[9][0] = 4;   table[9][1] = 5;
		table[10][0] = 3;   table[10][1] = 5;
		table[11][0] = 3;   table[11][1] = 4;
		table[12][0] = 0;   table[12][1] = 3;
		table[13][0] = 1;   table[13][1] = 4;
		table[14][0] = 2;   table[14][1] = 5;
		break;
	case ELEM_HEXA2:
		table[0][0] = 0;   table[0][1] = 0;
		table[1][0] = 1;   table[1][1] = 1;
		table[2][0] = 2;   table[2][1] = 2;
		table[3][0] = 3;   table[3][1] = 3;
		table[4][0] = 4;   table[4][1] = 4;
		table[5][0] = 5;   table[5][1] = 5;
		table[6][0] = 6;   table[6][1] = 6;
		table[7][0] = 7;   table[7][1] = 7;

		table[8][0] = 0;   table[8][1] = 1;
		table[9][0] = 1;   table[9][1] = 2;
		table[10][0] = 2;   table[10][1] = 3;
		table[11][0] = 0;   table[11][1] = 3;
		table[12][0] = 4;   table[12][1] = 5;
		table[13][0] = 5;   table[13][1] = 6;
		table[14][0] = 6;   table[14][1] = 7;
		table[15][0] = 4;   table[15][1] = 7;
		table[16][0] = 0;   table[16][1] = 4;
		table[17][0] = 1;   table[17][1] = 5;
		table[18][0] = 2;   table[18][1] = 6;
		table[19][0] = 3;   table[19][1] = 7;
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


/* １要素の平均長さ(代表長さ）を計算 */
/* 具体的には edge の長さ            */
RC element_average_length(FEM_ELEMENT element, FEM_NODE_ARRAY node,
		                  double *length)
{
	int table[12][2]; /* 12->HEXA2最大edge数 */
	int *index;
	int edge_num;
	int ii1;

	RC_TRY( allocate1D_i(element.node_num, &index) );
	for(ii1=0; ii1<element.node_num; ii1++){
		index[ii1] = search_fem_node_label(node, element.node[ii1]);
		RC_NEG_CHK(index[ii1]);
	}

	RC_TRY( element_edge_table(element.type, &edge_num, table) );

	*length = 0.0;
	for(ii1=0; ii1<edge_num; ii1++){
		*length += dist_point3d(node.array[ index[ table[ii1][0] ] ].p,
		                        node.array[ index[ table[ii1][1] ] ].p);
	}

	*length /= (double)edge_num;

	free( (void *)index );

	return(NORMAL_RC);
}


/* ２次要素の中間節点は考慮していない(直線扱い）*/
/* 中間節点の扱いは検討中                       */
RC element_edge_table(FEM_ELEM_TYPE type, int *edge_num, int table[][2])
{
	switch(type){
	case ELEM_BEAM1:
	case ELEM_LINE1:
		*edge_num = 1;
		table[0][0] = 0;    table[0][1] = 1;
		break;
	case ELEM_BEAM2:
	case ELEM_LINE2:
		*edge_num = 1;
		table[0][0] = 0;    table[0][1] = 1; /* table[0][2] = 2; */
		break;
	case ELEM_TRI1:
		*edge_num = 3;
		table[0][0] = 1;    table[0][1] = 2;
		table[1][0] = 2;    table[1][1] = 0;
		table[2][0] = 0;    table[2][1] = 1;
		break;
	case ELEM_TRI2:
		*edge_num = 3;
		table[0][0] = 1;    table[0][1] = 2; /* table[0][2] = 3; */
		table[1][0] = 2;    table[1][1] = 0; /* table[1][2] = 4; */
		table[2][0] = 0;    table[2][1] = 1; /* table[2][2] = 5; */
		break;
	case ELEM_QUAD1:
		*edge_num = 4;
		table[0][0] = 0;    table[0][1] = 1;
		table[1][0] = 1;    table[1][1] = 2;
		table[2][0] = 2;    table[2][1] = 3;
		table[3][0] = 3;    table[3][1] = 0;
		break;
	case ELEM_QUAD2:
		*edge_num = 4;
		table[0][0] = 0;    table[0][1] = 1; /* table[0][2] = 4; */
		table[1][0] = 1;    table[1][1] = 2; /* table[1][2] = 5; */
		table[2][0] = 2;    table[2][1] = 3; /* table[2][2] = 6; */
		table[3][0] = 3;    table[3][1] = 0; /* table[3][2] = 7; */
		break;
	case ELEM_TETRA1:
		*edge_num = 6;
		table[0][0] = 1;    table[0][1] = 2;
		table[1][0] = 2;    table[1][1] = 0;
		table[2][0] = 0;    table[2][1] = 1;
		table[3][0] = 0;    table[3][1] = 3;
		table[4][0] = 1;    table[4][1] = 3;
		table[5][0] = 2;    table[5][1] = 3;
		break;
	case ELEM_TETRA2:
		*edge_num = 6;
		table[0][0] = 1;    table[0][1] = 2; /* table[0][2] = 4; */
		table[1][0] = 2;    table[1][1] = 0; /* table[1][2] = 5; */
		table[2][0] = 0;    table[2][1] = 1; /* table[2][2] = 6; */
		table[3][0] = 0;    table[3][1] = 3; /* table[3][2] = 7; */
		table[4][0] = 1;    table[4][1] = 3; /* table[4][2] = 8; */
		table[5][0] = 2;    table[5][1] = 3; /* table[5][2] = 9; */
		break;
	case ELEM_PENTA1:
		*edge_num = 9;
		table[0][0] = 1;    table[0][1] = 2;
		table[1][0] = 2;    table[1][1] = 0;
		table[2][0] = 0;    table[2][1] = 1;
		table[3][0] = 4;    table[3][1] = 5;
		table[4][0] = 5;    table[4][1] = 3;
		table[5][0] = 3;    table[5][1] = 4;
		table[6][0] = 0;    table[6][1] = 3;
		table[7][0] = 1;    table[7][1] = 4;
		table[8][0] = 2;    table[8][1] = 5;
		break;
	case ELEM_PENTA2:
		*edge_num = 9;
		table[0][0] = 1;    table[0][1] = 2; /* table[0][2] = 6; */
		table[1][0] = 2;    table[1][1] = 0; /* table[1][2] = 7; */
		table[2][0] = 0;    table[2][1] = 1; /* table[2][2] = 8; */
		table[3][0] = 4;    table[3][1] = 5; /* table[3][2] = 9; */
		table[4][0] = 5;    table[4][1] = 3; /* table[4][2] = 10;*/
		table[5][0] = 3;    table[5][1] = 4; /* table[5][2] = 11;*/
		table[6][0] = 0;    table[6][1] = 3; /* table[6][2] = 12;*/
		table[7][0] = 1;    table[7][1] = 4; /* table[7][2] = 13;*/
		table[8][0] = 2;    table[8][1] = 5; /* table[8][2] = 14;*/
		break;
	case ELEM_HEXA1:
		*edge_num = 12;
		table[0][0] = 0;    table[0][1] = 1;
		table[1][0] = 1;    table[1][1] = 2;
		table[2][0] = 2;    table[2][1] = 3;
		table[3][0] = 3;    table[3][1] = 0;
		table[4][0] = 4;    table[4][1] = 5;
		table[5][0] = 5;    table[5][1] = 6;
		table[6][0] = 6;    table[6][1] = 7;
		table[7][0] = 7;    table[7][1] = 4;
		table[8][0] = 0;    table[8][1] = 4;
		table[9][0] = 1;    table[9][1] = 5;
		table[10][0] = 2;    table[10][1] = 6;
		table[11][0] = 3;    table[11][1] = 7;
		break;
	case ELEM_HEXA2:
		*edge_num = 12;
		table[0][0] = 0;    table[0][1] = 1; /* table[0][2] = 8; */
		table[1][0] = 1;    table[1][1] = 2; /* table[1][2] = 9; */
		table[2][0] = 2;    table[2][1] = 3; /* table[2][2] = 10; */
		table[3][0] = 3;    table[3][1] = 0; /* table[3][2] = 11; */
		table[4][0] = 4;    table[4][1] = 5; /* table[4][2] = 12; */
		table[5][0] = 5;    table[5][1] = 6; /* table[5][2] = 13; */
		table[6][0] = 6;    table[6][1] = 7; /* table[6][2] = 14; */
		table[7][0] = 7;    table[7][1] = 4; /* table[7][2] = 15; */
		table[8][0] = 0;    table[8][1] = 4; /* table[8][2] = 16; */
		table[9][0] = 1;    table[9][1] = 5; /* table[9][2] = 17; */
		table[10][0] = 2;    table[10][1] = 6; /* table[10][2] = 18;*/
		table[11][0] = 3;    table[11][1] = 7; /* table[11][2] = 19;*/
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


/* 連立方程式の解を、要素マトリックスの情報で求める */
/* solution[] に近似解(初期値)を代入すると、処理後は厳密解が代入されている */
/* 初期値と厳密解の差や、全体マトリックスの値や大きさによって、 */
/* ループ回数は変化する */ 
RC elem_matrix_cg(int total_dof, FEM_ELEM_MATRIX_ARRAY elem_matrix,
                  const double *vector, double *solution)
{
	double *r;             /* 残差ベクトル */
	double *p;             /* 最急勾配方向ベクトル */
	double *ar;            /* ar = a * r   a : 全体マトリックス */
	double *ap;            /* ap = a * p */ 
	double alpha;
	double beta;
	long double r_r;
	long double r1_r1;
	long double p_ar;
	double vector_value;
	double relative_error; /* 相対誤差 */
	int counter;           /* ループ回数 */
	int ii1;

	RC_TRY( allocate1D(total_dof, &r) );
	RC_TRY( allocate1D(total_dof, &p) );
	RC_TRY( allocate1D(total_dof, &ar) );
	RC_TRY( allocate1D(total_dof, &ap) );

	vector_value = 0.0;
	for(ii1=0; ii1<total_dof; ii1++){
		vector_value += vector[ii1] * vector[ii1];
	}
	
	counter = 0;

	/* 残差ベクトルの計算 (p[] = r[] = vector[] - a[][] * solution[]) */
	RC_TRY(mul_elem_matrices_vector(total_dof, elem_matrix, solution, ar));
	r_r = 0.0;
	for(ii1=0; ii1<total_dof; ii1++){
		p[ii1] = r[ii1] = vector[ii1] - ar[ii1];
		r_r += r[ii1] * r[ii1];
	}	
		
	while(1){
		RC_TRY( mul_elem_matrices_vector(total_dof, elem_matrix, r, ar) );
		p_ar = 0.0;
		for(ii1=0; ii1<total_dof; ii1++){
			p_ar += p[ii1] * ar[ii1];
		}

		/* 係数ベクトルのオーダー、ベクトルの次元数を考慮した */
		/* 相対誤差の計算 */
		relative_error = sqrt( (r_r / vector_value) / total_dof );
fprintf(stderr, "%d:%15.7e:%15.7e",
        counter, relative_error, REL_ERROR);
		if( relative_error <= REL_ERROR ) break;

		if( fabs(p_ar) <= ABS_ERROR ){
			fprintf(stderr, "p_ar < ");
			return (CAL_ERROR_RC);
		}
		alpha = r_r / p_ar;
fprintf(stderr, "%15.7e\n", alpha);

		/* solution[] の修正 */
		/* r[] の更新 */
		RC_TRY( mul_elem_matrices_vector(total_dof, elem_matrix, p, ap) );
		for(ii1=0; ii1<total_dof; ii1++){
			solution[ii1] += alpha * p[ii1];
			r[ii1] -= alpha * ap[ii1];
		}

		/* beta の計算 */
		r1_r1 = 0.0;
		for(ii1=0; ii1<total_dof; ii1++){
			r1_r1 += r[ii1] * r[ii1];
		}
		if( fabs(r_r) <= ABS_ERROR ){
			fprintf(stderr, "r_r < ");
			return (CAL_ERROR_RC);
		}
		beta = r1_r1 / r_r;

		/* r_r を更新 */
		r_r = r1_r1;
		
		/* 次回の修正方向の計算 */
		for(ii1=0; ii1<total_dof; ii1++){
			p[ii1] = r[ii1] + beta * p[ii1];
		}

		counter ++;
	}

	free( (void*)r );
	free( (void*)p );
	free( (void*)ar );
	free( (void*)ap );

	return(NORMAL_RC);
}
