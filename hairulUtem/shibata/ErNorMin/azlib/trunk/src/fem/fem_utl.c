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

/* $Id: fem_utl.c 1089 2016-11-30 09:36:47Z hayashi $ */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "base.h"
#include "mathematics.h"
#include "fem_struct.h"

#define BB_NUM 8  /* 双方向バブルソートの処理回数 */
/*#define REALLOC_SIZE 1024 */
#define MAX_INC_SIZE  (4096)

/* local function */
static int compar_label(int label1, int label2);
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
static int compar_contact_node(const void *p1, const void *p2);
static int search(const void *base, size_t nmemb, size_t size,
                  int (*compar)(const void *, const void *),
                  const void *target, ARRAY_SORT_FLAG sort_flag);
/* static int zero2max(int val, int max); */
static RC overwrite_bc(BC_TYPE *type1, double *v1,
                       BC_TYPE type2,  double v2);


/* 配列の解放 */
#define FUNC_FREE_ARRAY(name, type)   \
RC name(type *arg){\
	if(arg->alloc_size > 0) RC_TRY( mm_free(arg->array) );\
	arg->array = NULL;\
	arg->alloc_size = 0;\
	arg->size = 0;\
	return(NORMAL_RC);\
}

FUNC_FREE_ARRAY(free_node_array, NODE_ARRAY)

FUNC_FREE_ARRAY(free_element_array, ELEMENT_ARRAY)

FUNC_FREE_ARRAY(free_material_prop_array, MATERIAL_PROP_ARRAY)

FUNC_FREE_ARRAY(free_physical_prop_array, PHYSICAL_PROP_ARRAY)

FUNC_FREE_ARRAY(free_bc_array, BC_ARRAY)

FUNC_FREE_ARRAY(free_elem_bc_array, ELEM_BC_ARRAY)

FUNC_FREE_ARRAY(free_disp_array, DISP_ARRAY)

FUNC_FREE_ARRAY(free_velo_array, VELO_ARRAY)

FUNC_FREE_ARRAY(free_react_array, REACT_ARRAY)

FUNC_FREE_ARRAY(free_stress_array, STRESS_ARRAY)

FUNC_FREE_ARRAY(free_strain_array, STRAIN_ARRAY)

FUNC_FREE_ARRAY(free_sensitivity_array, SENSITIVITY_ARRAY)

FUNC_FREE_ARRAY(free_default_bc_array, DEFAULT_BC_ARRAY)

FUNC_FREE_ARRAY(free_sound_pressure_array, SOUND_PRESSURE_ARRAY)

FUNC_FREE_ARRAY(free_load_curve_array, LOAD_CURVE_ARRAY)

FUNC_FREE_ARRAY(free_local_coord_array, LOCAL_COORD_ARRAY)


/* 配列の解放(各要素の解放が必要な場合用) */
#define FUNC_FREE_ARRAY_F(name, type, free_func)   \
RC name(type *arg){\
	int ii1;\
	if(arg->alloc_size > 0){\
		for(ii1=0; ii1<(arg->size); ii1++){\
			RC_TRY( free_func(&(arg->array[ii1])) );\
		}\
		RC_TRY( mm_free(arg->array) );\
	}\
	arg->array = NULL;\
	arg->alloc_size = 0;\
	arg->size = 0;\
	return(NORMAL_RC);\
}

FUNC_FREE_ARRAY_F(free_elem_matrix_array, ELEM_MATRIX_ARRAY, free_elem_matrix)

FUNC_FREE_ARRAY_F(free_elem_lookup_array, ELEM_LOOKUP_ARRAY, free_elem_lookup)

FUNC_FREE_ARRAY_F(free_node_matrix_array, NODE_MATRIX_ARRAY, free_node_matrix)

FUNC_FREE_ARRAY_F(free_rigid_element_array, RIGID_ELEMENT_ARRAY,
                  free_rigid_element)

FUNC_FREE_ARRAY_F(free_bc_set_array, BC_SET_ARRAY, free_bc_set)

FUNC_FREE_ARRAY_F(free_contact_array, CONTACT_ARRAY, free_contact)


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

FUNC_COPY_ARRAY(copy_node_array, NODE_ARRAY, NODE)

FUNC_COPY_ARRAY(copy_element_array, ELEMENT_ARRAY, ELEMENT)

FUNC_COPY_ARRAY(copy_material_prop_array, MATERIAL_PROP_ARRAY, MATERIAL_PROP)

FUNC_COPY_ARRAY(copy_physical_prop_array, PHYSICAL_PROP_ARRAY, PHYSICAL_PROP)

FUNC_COPY_ARRAY(copy_bc_array, BC_ARRAY, BC)

FUNC_COPY_ARRAY(copy_elem_bc_array, ELEM_BC_ARRAY, ELEM_BC)

FUNC_COPY_ARRAY(copy_disp_array, DISP_ARRAY, DISP_VELO)

FUNC_COPY_ARRAY(copy_velo_array, VELO_ARRAY, DISP_VELO)

FUNC_COPY_ARRAY(copy_react_array, REACT_ARRAY, DISP_VELO)

FUNC_COPY_ARRAY(copy_stress_array, STRESS_ARRAY, STRESS_STRAIN)

FUNC_COPY_ARRAY(copy_strain_array, STRAIN_ARRAY, STRESS_STRAIN)

FUNC_COPY_ARRAY(copy_sensitivity_array, SENSITIVITY_ARRAY, SENSITIVITY)

FUNC_COPY_ARRAY(copy_default_bc_array, DEFAULT_BC_ARRAY, DEFAULT_BC)

FUNC_COPY_ARRAY(copy_sound_pressure_array, SOUND_PRESSURE_ARRAY, SOUND_PRESSURE)

FUNC_COPY_ARRAY(copy_load_curve_array, LOAD_CURVE_ARRAY, LOAD_CURVE)

FUNC_COPY_ARRAY(copy_local_coord_array, LOCAL_COORD_ARRAY, LOCAL_COORD)


/* clean_*() は，構造体の中の配列で使用していない要素のメモリを解放する．
 */
#define FUNC_CLEAN_ARRAY(name, array_type, type, sort_func, elem)  \
RC name(array_type *ptr){\
	int ii1;\
\
	if(ptr->size > ptr->alloc_size) return(ARG_ERROR_RC);\
	if(ptr->size < 0) return(ARG_ERROR_RC);\
\
	RC_TRY( sort_func(ptr) );\
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
			tmp = mm_realloc(ptr->array, ptr->size*sizeof(type));\
			if(tmp == NULL) return(ALLOC_ERROR_RC);\
			ptr->array = (type *)tmp;\
			ptr->alloc_size = ptr->size;\
		}else{\
			RC_TRY( mm_free(ptr->array) );\
			ptr->array = NULL;\
			ptr->size = 0;\
			ptr->alloc_size = 0;\
		}\
	}\
\
	return(NORMAL_RC);\
}

FUNC_CLEAN_ARRAY(clean_node_array, NODE_ARRAY, NODE,
                 sort_node_array, label)

FUNC_CLEAN_ARRAY(clean_element_array, ELEMENT_ARRAY, ELEMENT,
                 sort_element_array, label)

FUNC_CLEAN_ARRAY(clean_rigid_element_array, RIGID_ELEMENT_ARRAY,
                 RIGID_ELEMENT, sort_rigid_element_array, label)

FUNC_CLEAN_ARRAY(clean_material_prop_array, MATERIAL_PROP_ARRAY,
                 MATERIAL_PROP, sort_material_prop_array, label)

FUNC_CLEAN_ARRAY(clean_physical_prop_array, PHYSICAL_PROP_ARRAY,
                 PHYSICAL_PROP, sort_physical_prop_array, label)

FUNC_CLEAN_ARRAY(clean_bc_array, BC_ARRAY,
                 BC, sort_bc_array, node)

FUNC_CLEAN_ARRAY(clean_elem_bc_array, ELEM_BC_ARRAY,
                 ELEM_BC, sort_elem_bc_array, element)

FUNC_CLEAN_ARRAY(clean_disp_array, DISP_ARRAY,
                 DISP_VELO, sort_disp_array, node)

FUNC_CLEAN_ARRAY(clean_velo_array, VELO_ARRAY,
                 DISP_VELO, sort_velo_array, node)

FUNC_CLEAN_ARRAY(clean_react_array, REACT_ARRAY,
                 DISP_VELO, sort_react_array, node)

FUNC_CLEAN_ARRAY(clean_stress_array, STRESS_ARRAY,
                 STRESS_STRAIN, sort_stress_array, element)

FUNC_CLEAN_ARRAY(clean_strain_array, STRAIN_ARRAY,
                 STRESS_STRAIN, sort_strain_array, element)

FUNC_CLEAN_ARRAY(clean_sensitivity_array, SENSITIVITY_ARRAY,
                 SENSITIVITY, sort_sensitivity_array, element)

FUNC_CLEAN_ARRAY(clean_default_bc_array, DEFAULT_BC_ARRAY,
                 DEFAULT_BC, sort_default_bc_array, set_id)

FUNC_CLEAN_ARRAY(clean_sound_pressure_array, SOUND_PRESSURE_ARRAY,
                 SOUND_PRESSURE, sort_sound_pressure_array, node)

FUNC_CLEAN_ARRAY(clean_bc_set_array, BC_SET_ARRAY,
                 BC_SET, sort_bc_set_array, label)

FUNC_CLEAN_ARRAY(clean_elem_matrix_array, ELEM_MATRIX_ARRAY,
                 ELEM_MATRIX, sort_elem_matrix_array, element)

FUNC_CLEAN_ARRAY(clean_elem_lookup_array, ELEM_LOOKUP_ARRAY,
                 ELEM_LOOKUP, sort_elem_lookup_array, node)

FUNC_CLEAN_ARRAY(clean_node_matrix_array, NODE_MATRIX_ARRAY,
                 NODE_MATRIX, sort_node_matrix_array, node)

FUNC_CLEAN_ARRAY(clean_load_curve_array, LOAD_CURVE_ARRAY,
                 LOAD_CURVE, sort_load_curve_array, label)

FUNC_CLEAN_ARRAY(clean_local_coord_array, LOCAL_COORD_ARRAY,
                 LOCAL_COORD, sort_local_coord_array, label)

FUNC_CLEAN_ARRAY(clean_contact_array, CONTACT_ARRAY,
                 CONTACT, sort_contact_array, node)


/* realloc_*() は、最低でも要素１つ追加できるよう array を再確保する．
 * 通常は，確保する要素の数を倍にする．ただし，最大値は MAX_INC_SIZE である．
 */
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
		tmp = (type *)mm_realloc((void *)arg->array,\
		                         realloc_size*sizeof(type));\
		if(tmp == NULL) return(ALLOC_ERROR_RC);\
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

FUNC_REALLOC_ARRAY(realloc_node_array, NODE_ARRAY,
                   NODE, init_node)

FUNC_REALLOC_ARRAY(realloc_element_array, ELEMENT_ARRAY,
                   ELEMENT, init_element)

FUNC_REALLOC_ARRAY(realloc_rigid_element_array, RIGID_ELEMENT_ARRAY,
                   RIGID_ELEMENT, init_rigid_element)

FUNC_REALLOC_ARRAY(realloc_material_prop_array, MATERIAL_PROP_ARRAY,
                   MATERIAL_PROP, init_material_prop)

FUNC_REALLOC_ARRAY(realloc_physical_prop_array, PHYSICAL_PROP_ARRAY,
                   PHYSICAL_PROP, init_physical_prop)

FUNC_REALLOC_ARRAY(realloc_bc_array, BC_ARRAY,
                   BC, init_bc)

FUNC_REALLOC_ARRAY(realloc_elem_bc_array, ELEM_BC_ARRAY,
                   ELEM_BC, init_elem_bc)

FUNC_REALLOC_ARRAY(realloc_disp_array, DISP_ARRAY,
                   DISP_VELO, init_disp_velo)

FUNC_REALLOC_ARRAY(realloc_velo_array, VELO_ARRAY,
                   DISP_VELO, init_disp_velo)

FUNC_REALLOC_ARRAY(realloc_react_array, REACT_ARRAY,
                   DISP_VELO, init_disp_velo)

FUNC_REALLOC_ARRAY(realloc_stress_array, STRESS_ARRAY,
                   STRESS_STRAIN, init_stress_strain)

FUNC_REALLOC_ARRAY(realloc_strain_array, STRAIN_ARRAY,
                   STRESS_STRAIN, init_stress_strain)

FUNC_REALLOC_ARRAY(realloc_sensitivity_array, SENSITIVITY_ARRAY,
                   SENSITIVITY, init_sensitivity)

FUNC_REALLOC_ARRAY(realloc_default_bc_array, DEFAULT_BC_ARRAY,
                   DEFAULT_BC, init_default_bc)

FUNC_REALLOC_ARRAY(realloc_sound_pressure_array, SOUND_PRESSURE_ARRAY,
                   SOUND_PRESSURE, init_sound_pressure)

FUNC_REALLOC_ARRAY(realloc_bc_set_array, BC_SET_ARRAY,
                   BC_SET, init_bc_set)

FUNC_REALLOC_ARRAY(realloc_elem_matrix_array, ELEM_MATRIX_ARRAY,
                   ELEM_MATRIX, init_elem_matrix)

FUNC_REALLOC_ARRAY(realloc_elem_lookup_array, ELEM_LOOKUP_ARRAY,
                   ELEM_LOOKUP, init_elem_lookup)

FUNC_REALLOC_ARRAY(realloc_node_matrix_array, NODE_MATRIX_ARRAY,
                   NODE_MATRIX, init_node_matrix)

FUNC_REALLOC_ARRAY(realloc_load_curve_array, LOAD_CURVE_ARRAY,
                   LOAD_CURVE, init_load_curve)

FUNC_REALLOC_ARRAY(realloc_local_coord_array, LOCAL_COORD_ARRAY,
                   LOCAL_COORD, init_local_coord)

FUNC_REALLOC_ARRAY(realloc_contact_array, CONTACT_ARRAY,
                   CONTACT, init_contact)


RC
allocate_node_array (int num, NODE_ARRAY *node)
{
	int ii1;

	if(num > 0){
		node->array = mm_alloc(num * sizeof(NODE));
		if(node->array == NULL) return(ALLOC_ERROR_RC);
	}else{
		node->array = NULL;
	}

	node->size = 0;
	node->alloc_size = num;
	node->source = FEM_UNKNOWN;
	node->renum_flag = UNRENUMBERED;
	node->sort_flag = ARRAY_UNSORTED;

	for(ii1=0;ii1<num;ii1++){
		init_node( &((node->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC
allocate_element_array (int num, ELEMENT_ARRAY *element)
{
	int ii1;

	if(num > 0){
		element->array = mm_alloc(num * sizeof(ELEMENT));
		if(element->array == NULL) return(ALLOC_ERROR_RC);
	}else{
		element->array = NULL;
	}

	element->size = 0;
	element->alloc_size = num;
	element->source = FEM_UNKNOWN;
	element->sort_flag = ARRAY_UNSORTED;

	for(ii1=0;ii1<num;ii1++){
		init_element( &((element->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC
allocate_rigid_element_array (int num, RIGID_ELEMENT_ARRAY *element)
{
	int ii1;

	if(num > 0){
		element->array = mm_alloc(num * sizeof(RIGID_ELEMENT));
		if(element->array == NULL) return(ALLOC_ERROR_RC);
	}else{
		element->array = NULL;
	}

	element->size = 0;
	element->alloc_size = num;
	element->source = FEM_UNKNOWN;
	element->sort_flag = ARRAY_UNSORTED;

	for(ii1=0;ii1<num;ii1++){
		init_rigid_element( &((element->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC
allocate_material_prop_array (int num, MATERIAL_PROP_ARRAY *mat)
{
	int ii1;

	if(num > 0){
		mat->array = mm_alloc(num * sizeof(MATERIAL_PROP));
		if(mat->array == NULL) return(ALLOC_ERROR_RC);
	}else{
		mat->array = NULL;
	}

	mat->size = 0;
	mat->alloc_size = num;
	mat->source = FEM_UNKNOWN;
	mat->sort_flag = ARRAY_UNSORTED;

	for(ii1=0;ii1<num;ii1++){
		init_material_prop( &((mat->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC
allocate_physical_prop_array (int num, PHYSICAL_PROP_ARRAY *phys)
{
	int ii1;

	if(num > 0){
		phys->array = mm_alloc(num * sizeof(PHYSICAL_PROP));
		if(phys->array == NULL) return(ALLOC_ERROR_RC);
	}else{
		phys->array = NULL;
	}

	phys->size = 0;
	phys->alloc_size = num;
	phys->source = FEM_UNKNOWN;
	phys->sort_flag = ARRAY_UNSORTED;

	for(ii1=0;ii1<num;ii1++){
		init_physical_prop( &((phys->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC
allocate_bc_array (int num, BC_ARRAY *bc)
{
	int ii1;

	if(num > 0){
		bc->array = mm_alloc(num * sizeof(BC));
		if(bc->array == NULL) return(ALLOC_ERROR_RC);
	}else{
		bc->array = NULL;
	}

	bc->size = 0;
	bc->alloc_size = num;
	bc->source = FEM_UNKNOWN;
	bc->sort_flag = ARRAY_UNSORTED;

	for(ii1=0;ii1<num;ii1++){
		init_bc( &((bc->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC
allocate_elem_bc_array (int num, ELEM_BC_ARRAY *bc)
{
	int ii1;

	if(num > 0){
		bc->array = mm_alloc(num * sizeof(ELEM_BC));
		if(bc->array == NULL) return(ALLOC_ERROR_RC);
	}else{
		bc->array = NULL;
	}

	bc->size = 0;
	bc->alloc_size = num;
	bc->source = FEM_UNKNOWN;
	bc->sort_flag = ARRAY_UNSORTED;

	for(ii1=0;ii1<num;ii1++){
		init_elem_bc( &((bc->array)[ii1]) );
	}

	return(NORMAL_RC);
}


/* num個の要素を持つ disp の配列を動的確保する．
 */
RC
allocate_disp_array (int num, DISP_ARRAY *disp)
{
	int ii1;

	if(num > 0){
		disp->array = mm_alloc(num * sizeof(DISP_VELO));
		if(disp->array == NULL) return(ALLOC_ERROR_RC);
	}else{
		disp->array = NULL;
	}

	disp->size = 0;
	disp->alloc_size = num;
	disp->source = FEM_UNKNOWN;
	disp->sort_flag = ARRAY_UNSORTED;

	for(ii1=0;ii1<num;ii1++){
		init_disp_velo( &((disp->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC
allocate_velo_array (int num, VELO_ARRAY *velo)
{
	int ii1;

	if(num > 0){
		velo->array = mm_alloc(num * sizeof(DISP_VELO));
		if(velo->array == NULL) return(ALLOC_ERROR_RC);
	}else{
		velo->array = NULL;
	}
	
	velo->size = 0;
	velo->alloc_size = num;
	velo->source = FEM_UNKNOWN;
	velo->sort_flag = ARRAY_UNSORTED;

	for(ii1=0;ii1<num;ii1++){
		init_disp_velo( &((velo->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC
allocate_react_array (int num, REACT_ARRAY *react)
{
	int ii1;

	if(num > 0){
		react->array = mm_alloc(num * sizeof(DISP_VELO));
		if(react->array == NULL) return(ALLOC_ERROR_RC);
	}else{
		react->array = NULL;
	}
	
	react->size = 0;
	react->alloc_size = num;
	react->source = FEM_UNKNOWN;
	react->sort_flag = ARRAY_UNSORTED;

	for(ii1=0;ii1<num;ii1++){
		init_disp_velo( &((react->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC
allocate_stress_array (int num, STRESS_ARRAY *stress)
{
	int ii1;

	if(num > 0){
		stress->array = mm_alloc(num * sizeof(STRESS_STRAIN));
		if(stress->array == NULL) return(ALLOC_ERROR_RC);
	}else{
		stress->array = NULL;
	}

	stress->size = 0;
	stress->alloc_size = num;
	stress->source = FEM_UNKNOWN;
	stress->sort_flag = ARRAY_UNSORTED;

	for(ii1=0;ii1<num;ii1++){
		init_stress_strain( &((stress->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC
allocate_strain_array (int num, STRAIN_ARRAY *strain)
{
	int ii1;

	if(num > 0){
		strain->array = mm_alloc(num * sizeof(STRESS_STRAIN));
		if(strain->array == NULL) return(ALLOC_ERROR_RC);
	}else{
		strain->array = NULL;
	}

	strain->size = 0;
	strain->alloc_size = num;
	strain->source = FEM_UNKNOWN;
	strain->sort_flag = ARRAY_UNSORTED;

	for(ii1=0;ii1<num;ii1++){
		init_stress_strain( &((strain->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC
allocate_sensitivity_array (int num, SENSITIVITY_ARRAY *sens)
{
	int ii1;

	if(num > 0){
		sens->array = mm_alloc(num * sizeof(SENSITIVITY));
		if(sens->array == NULL) return(ALLOC_ERROR_RC);
	}else{
		sens->array = NULL;
	}

	sens->size = 0;
	sens->alloc_size = num;
	sens->type = SENS_UNKNOWN;
	sens->sort_flag = ARRAY_UNSORTED;

	for(ii1=0;ii1<num;ii1++){
		init_sensitivity( &((sens->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC
allocate_default_bc_array (int num, DEFAULT_BC_ARRAY *def)
{
	int ii1;

	if(num > 0){
		def->array = mm_alloc(num * sizeof(DEFAULT_BC));
		if(def->array == NULL) return(ALLOC_ERROR_RC);
	}else{
		def->array = NULL;
	}

	def->size = 0;
	def->alloc_size = num;
	def->source = FEM_UNKNOWN;
	def->sort_flag = ARRAY_UNSORTED;

	for(ii1=0; ii1<num; ii1++){
		init_default_bc( &((def->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC
allocate_sound_pressure_array (int num, SOUND_PRESSURE_ARRAY *sound)
{
	int ii1;

	if(num > 0){
		sound->array = mm_alloc(num * sizeof(SOUND_PRESSURE));
		if(sound->array == NULL) return(ALLOC_ERROR_RC);
	}else{
		sound->array = NULL;
	}

	sound->size = 0;
	sound->alloc_size = num;
	sound->source = FEM_UNKNOWN;
	sound->sort_flag = ARRAY_UNSORTED;

	for(ii1=0; ii1<num; ii1++){
		init_sound_pressure( &((sound->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC
allocate_bc_set_array (int num, BC_SET_ARRAY *bc_set)
{
	int ii1;

	if(num > 0){
		bc_set->array = mm_alloc(num * sizeof(BC_SET));
		if(bc_set->array == NULL) return(ALLOC_ERROR_RC);
	}else{
		bc_set->array = NULL;
	}

	bc_set->size = 0;
	bc_set->alloc_size = num;
	bc_set->source = FEM_UNKNOWN;
	bc_set->sort_flag = ARRAY_UNSORTED;
	bc_set->com_fix = BC_NOT_COMMON;
	bc_set->com_mpc = BC_NOT_COMMON;
	bc_set->com_force = BC_NOT_COMMON;
	bc_set->com_temp = BC_NOT_COMMON;

	for(ii1=0; ii1<num; ii1++){
		init_bc_set( &((bc_set->array)[ii1]) );
	}
	RC_TRY( init_string_array(&(bc_set->s_info)) );

	return(NORMAL_RC);
}


RC
allocate_elem_matrix_array (int num, ELEM_MATRIX_ARRAY *elem_matrix)
{
	int ii1;

	if(num > 0){
		elem_matrix->array = mm_alloc(num * sizeof(ELEM_MATRIX));
		if(elem_matrix->array == NULL) return(ALLOC_ERROR_RC);
	}else{
		elem_matrix->array = NULL;
	}

	elem_matrix->size = 0;
	elem_matrix->alloc_size = num;
	elem_matrix->source = FEM_UNKNOWN;
	elem_matrix->sort_flag = ARRAY_UNSORTED;

	for(ii1=0; ii1<num; ii1++){
		init_elem_matrix( &((elem_matrix->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC
allocate_elem_lookup_array (int num, ELEM_LOOKUP_ARRAY *elem_lookup)
{
	int ii1;

	if(num > 0){
		elem_lookup->array = mm_alloc(num * sizeof(ELEM_LOOKUP));
		if(elem_lookup->array == NULL) return(ALLOC_ERROR_RC);
	}else{
		elem_lookup->array = NULL;
	}

	elem_lookup->size = 0;
	elem_lookup->alloc_size = num;
	elem_lookup->source = FEM_UNKNOWN;
	elem_lookup->sort_flag = ARRAY_UNSORTED;

	for(ii1=0; ii1<num; ii1++){
		init_elem_lookup( &((elem_lookup->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC
allocate_node_matrix_array (int num, NODE_MATRIX_ARRAY *nmat)
{
	int ii1;

	if(num > 0){
		nmat->array = mm_alloc(num * sizeof(NODE_MATRIX));
		if(nmat->array == NULL) return(ALLOC_ERROR_RC);
	}else{
		nmat->array = NULL;
	}

	nmat->size = 0;
	nmat->alloc_size = num;
	nmat->source = FEM_UNKNOWN;
	nmat->sort_flag = ARRAY_UNSORTED;

	for(ii1=0; ii1<num; ii1++){
		init_node_matrix( &((nmat->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC
allocate_load_curve_array (int num, LOAD_CURVE_ARRAY *curve)
{
	int ii1;

	if(num > 0){
		curve->array = mm_alloc(num * sizeof(LOAD_CURVE));
		if(curve->array == NULL) return(ALLOC_ERROR_RC);
	}else{
		curve->array = NULL;
	}

	curve->size = 0;
	curve->alloc_size = num;
	curve->source = FEM_UNKNOWN;
	curve->sort_flag = ARRAY_UNSORTED;

	for(ii1=0; ii1<num; ii1++){
		init_load_curve( &((curve->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC
allocate_local_coord_array (int num, LOCAL_COORD_ARRAY *local)
{
	int ii1;

	if(num > 0){
		local->array = mm_alloc(num * sizeof(LOCAL_COORD));
		if(local->array == NULL) return(ALLOC_ERROR_RC);
	}else{
		local->array = NULL;
	}

	local->size = 0;
	local->alloc_size = num;
	local->source = FEM_UNKNOWN;
	local->sort_flag = ARRAY_UNSORTED;

	for(ii1=0; ii1<num; ii1++){
		init_local_coord( &((local->array)[ii1]) );
	}

	return(NORMAL_RC);
}


RC
allocate_contact_array (int num, CONTACT_ARRAY *contact)
{
	int ii1;

	if(num > 0){
		contact->array = mm_alloc(num * sizeof(CONTACT));
		if(contact->array == NULL) return(ALLOC_ERROR_RC);
	}else{
		contact->array = NULL;
	}

	contact->size = 0;
	contact->alloc_size = num;
	contact->sort_flag = ARRAY_UNSORTED;

	for(ii1=0; ii1<num; ii1++){
		init_contact( &((contact->array)[ii1]) );
	}

	return(NORMAL_RC);
}


/* count_valid_*() は，配列中で使用中の要素（labelなどの値が非負であるもの）の
 * 数を返す．
 */
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

FUNC_COUNT_VALID(count_valid_node, NODE_ARRAY, label)

FUNC_COUNT_VALID(count_valid_element, ELEMENT_ARRAY, label)

FUNC_COUNT_VALID(count_valid_rigid_element, RIGID_ELEMENT_ARRAY, label)

FUNC_COUNT_VALID(count_valid_material_prop, MATERIAL_PROP_ARRAY, label)

FUNC_COUNT_VALID(count_valid_physical_prop, PHYSICAL_PROP_ARRAY, label)

FUNC_COUNT_VALID(count_valid_bc, BC_ARRAY, node)

FUNC_COUNT_VALID(count_valid_elem_bc, ELEM_BC_ARRAY, element)

FUNC_COUNT_VALID(count_valid_disp, DISP_ARRAY, node)

FUNC_COUNT_VALID(count_valid_velo, VELO_ARRAY, node)

FUNC_COUNT_VALID(count_valid_react, REACT_ARRAY, node)

FUNC_COUNT_VALID(count_valid_load_curve, LOAD_CURVE_ARRAY, label)

FUNC_COUNT_VALID(count_valid_local_coord, LOCAL_COORD_ARRAY, label)

FUNC_COUNT_VALID(count_valid_contact, CONTACT_ARRAY, node)


/* special comparison of label */
static int
compar_label (int label1, int label2)
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


RC
sort_node_array (NODE_ARRAY *node)
{
	RC_TRY( bbqsort(node->array, node->size, sizeof(NODE),
	                                    compar_node_label) );
	node->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC
sort_node_renum (NODE_ARRAY *node)
{
	RC_TRY( bbqsort(node->array, node->size, sizeof(NODE),
	                                    compar_node_renum) );
	node->sort_flag = ARRAY_UNSORTED;
	
	return(NORMAL_RC);
}


RC
sort_element_array (ELEMENT_ARRAY *elem)
{
	RC_TRY( bbqsort(elem->array, elem->size, sizeof(ELEMENT),
	                                    compar_element_label) );
	elem->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC
sort_rigid_element_array (RIGID_ELEMENT_ARRAY *elem)
{
	RC_TRY( bbqsort(elem->array, elem->size, sizeof(RIGID_ELEMENT),
	                                    compar_rigid_element_label) );
	elem->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC
sort_material_prop_array (MATERIAL_PROP_ARRAY *mat)
{
	RC_TRY( bbqsort(mat->array, mat->size, sizeof(MATERIAL_PROP),
	                                       compar_material_prop_label) );
	mat->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC
sort_physical_prop_array (PHYSICAL_PROP_ARRAY *phys)
{
	RC_TRY( bbqsort(phys->array, phys->size, sizeof(PHYSICAL_PROP),
	                                         compar_physical_prop_label) );
	phys->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC
sort_bc_array (BC_ARRAY *bc)
{
	RC_TRY( bbqsort(bc->array, bc->size, sizeof(BC),
	                                 compar_bc_node) );
	bc->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC
sort_elem_bc_array (ELEM_BC_ARRAY *bc)
{
	RC_TRY( bbqsort(bc->array, bc->size, sizeof(ELEM_BC),
	                                 compar_elem_bc_element) );
	bc->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC
sort_disp_array (DISP_ARRAY *disp)
{
	RC_TRY( bbqsort(disp->array, disp->size, sizeof(DISP_VELO),
	                                     compar_disp_velo_node) );
	disp->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC
sort_velo_array (VELO_ARRAY *velo)
{
	RC_TRY( bbqsort(velo->array, velo->size, sizeof(DISP_VELO),
	                                     compar_disp_velo_node) );
	velo->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC
sort_react_array (REACT_ARRAY *react)
{
	RC_TRY( bbqsort(react->array, react->size, sizeof(DISP_VELO),
	                                     compar_disp_velo_node) );
	react->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC
sort_stress_array (STRESS_ARRAY *stress)
{
	RC_TRY( bbqsort(stress->array, stress->size, sizeof(STRESS_STRAIN),
	                                      compar_stress_strain_element) );
	stress->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC
sort_strain_array (STRAIN_ARRAY *strain)
{
	RC_TRY( bbqsort(strain->array, strain->size, sizeof(STRESS_STRAIN),
	                                      compar_stress_strain_element) );
	strain->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC
sort_sensitivity_array (SENSITIVITY_ARRAY *sens)
{
	RC_TRY( bbqsort(sens->array, sens->size, sizeof(SENSITIVITY),
	                                  compar_sensitivity_element) );
	sens->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC
sort_default_bc_array (DEFAULT_BC_ARRAY *def)
{
	RC_TRY( bbqsort(def->array, def->size, sizeof(DEFAULT_BC),
	                                  compar_default_bc_set_id) );
	def->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC
sort_sound_pressure_array (SOUND_PRESSURE_ARRAY *sound)
{
	RC_TRY( bbqsort(sound->array, sound->size,
	                sizeof(SOUND_PRESSURE), compar_sound_pressure_node) );
	sound->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC
sort_bc_set_array (BC_SET_ARRAY *bc_set)
{
	RC_TRY( bbqsort(bc_set->array, bc_set->size, sizeof(BC_SET),
	                                        compar_bc_set_label) );
	bc_set->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC
sort_elem_matrix_array (ELEM_MATRIX_ARRAY *elem_matrix)
{
	RC_TRY( bbqsort(elem_matrix->array, elem_matrix->size,
	                sizeof(ELEM_MATRIX), compar_elem_matrix_element) );
	elem_matrix->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC
sort_elem_lookup_array (ELEM_LOOKUP_ARRAY *elem_lookup)
{
	RC_TRY( bbqsort(elem_lookup->array, elem_lookup->size,
	                sizeof(ELEM_LOOKUP), compar_elem_lookup_node) );
	elem_lookup->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC
sort_node_matrix_array (NODE_MATRIX_ARRAY *nmat)
{
	RC_TRY( bbqsort(nmat->array, nmat->size, sizeof(NODE_MATRIX),
	                                     compar_node_matrix_node) );
	nmat->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC
sort_load_curve_array (LOAD_CURVE_ARRAY *curve)
{
	RC_TRY( bbqsort(curve->array, curve->size, sizeof(LOAD_CURVE),
	                                      compar_load_curve_label) );
	curve->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC
sort_local_coord_array (LOCAL_COORD_ARRAY *coord)
{
	RC_TRY( bbqsort(coord->array, coord->size,
	                sizeof(LOCAL_COORD), compar_local_coord_label) );
	coord->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


RC
sort_contact_array (CONTACT_ARRAY *contact)
{
	RC_TRY( bbqsort(contact->array, contact->size,
	                sizeof(CONTACT), compar_contact_node) );
	contact->sort_flag = ARRAY_SORTED;
	
	return(NORMAL_RC);
}


/* for bbqsort() */
static int
compar_node_label (const void *p1, const void *p2)
{
	return(compar_label( ((NODE *)p1)->label,
	                     ((NODE *)p2)->label ));
}

static int
compar_node_renum (const void *p1, const void *p2)
{
	return(compar_label( ((NODE *)p1)->renum,
	                     ((NODE *)p2)->renum ));
}

static int
compar_element_label (const void *p1, const void *p2)
{
	return(compar_label( ((ELEMENT *)p1)->label,
	                     ((ELEMENT *)p2)->label ));
}

static int
compar_rigid_element_label (const void *p1, const void *p2)
{
	return(compar_label( ((RIGID_ELEMENT *)p1)->label,
	                     ((RIGID_ELEMENT *)p2)->label ));
}

static int
compar_material_prop_label (const void *p1, const void *p2)
{
	return(compar_label( ((MATERIAL_PROP *)p1)->label,
	                     ((MATERIAL_PROP *)p2)->label ));
}

static int
compar_physical_prop_label (const void *p1, const void *p2)
{
	return(compar_label( ((PHYSICAL_PROP *)p1)->label,
	                     ((PHYSICAL_PROP *)p2)->label ));
}

static int
compar_bc_node (const void *p1, const void *p2)
{
	return(compar_label( ((BC *)p1)->node,
	                     ((BC *)p2)->node ));
}

static int
compar_elem_bc_element (const void *p1, const void *p2)
{
	return(compar_label( ((ELEM_BC *)p1)->element,
	                     ((ELEM_BC *)p2)->element ));
}

static int
compar_disp_velo_node (const void *p1, const void *p2)
{
	return(compar_label( ((DISP_VELO *)p1)->node,
	                     ((DISP_VELO *)p2)->node ));
}

static int
compar_stress_strain_element (const void *p1, const void *p2)
{
	return(compar_label( ((STRESS_STRAIN *)p1)->element,
	                     ((STRESS_STRAIN *)p2)->element ));
}

static int
compar_sensitivity_element (const void *p1, const void *p2)
{
	return(compar_label( ((SENSITIVITY *)p1)->element,
	                     ((SENSITIVITY *)p2)->element ));
}

static int
compar_default_bc_set_id (const void *p1, const void *p2)
{
	return(compar_label( ((DEFAULT_BC *)p1)->set_id,
	                     ((DEFAULT_BC *)p2)->set_id ));
}

static int
compar_sound_pressure_node (const void *p1, const void *p2)
{
	return(compar_label( ((SOUND_PRESSURE *)p1)->node,
	                     ((SOUND_PRESSURE *)p2)->node ));
}

static int
compar_bc_set_label (const void *p1, const void *p2)
{
	return(compar_label( ((BC_SET *)p1)->label,
	                     ((BC_SET *)p2)->label ));
}

static int
compar_elem_matrix_element (const void *p1, const void *p2)
{
	return(compar_label( ((ELEM_MATRIX *)p1)->element,
	                     ((ELEM_MATRIX *)p2)->element ));
}

static int
compar_elem_lookup_node (const void *p1, const void *p2)
{
	return(compar_label( ((ELEM_LOOKUP *)p1)->node,
	                     ((ELEM_LOOKUP *)p2)->node ));
}

static int
compar_node_matrix_node (const void *p1, const void *p2)
{
	return(compar_label( ((NODE_MATRIX *)p1)->node,
	                     ((NODE_MATRIX *)p2)->node ));
}

static int
compar_load_curve_label (const void *p1, const void *p2)
{
	return(compar_label( ((LOAD_CURVE *)p1)->label,
	                     ((LOAD_CURVE *)p2)->label ));
}

static int
compar_local_coord_label (const void *p1, const void *p2)
{
	return(compar_label( ((LOCAL_COORD *)p1)->label,
	                     ((LOCAL_COORD *)p2)->label ));
}

static int
compar_contact_node (const void *p1, const void *p2)
{
	return(compar_label( ((CONTACT *)p1)->node,
	                     ((CONTACT *)p2)->node ));
}


/* --探査プログラム--
 * ポインタ base から要素サイズ(size)，要素数(nmemb)の配列から targetと
 * 同じ属性の要素を探し，そのインデックスを返却する．要素同士の比較は
 * compar()を使用する．探査に失敗すれば、-1 を返却する．
 */ 
static int
search (const void *base, size_t nmemb, size_t size,
        int (*compar)(const void *, const void *),
        const void *target, ARRAY_SORT_FLAG sort_flag)
{
	int found_index = -1;

	if((base == NULL)||(target == NULL)||(nmemb <= 0)) return(-2);

	if(sort_flag == ARRAY_SORTED){
		/* ソートされている場合には二分探査を行う */
		int low = 0;
		int high = nmemb - 1;
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


/* 探査プログラムは labelと indexが一致している場合を想定して，まず
 * 該当する index を検査し，一致していない場合に限って，search()を
 * 使用する．同じ label が連続している場合はできる限り小さい index を
 * 返却する．
 */
#define FUNC_SEARCH_ARRAY(name, array_type, type, elem, init_f, compar_f)  \
int name(array_type arg, int label)\
{\
	int index0, index1, index2, index3, itmp;\
	int found_index = -1;\
	type dum;\
\
	if(arg.array == NULL) return(-2);\
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


FUNC_SEARCH_ARRAY(search_node_label, NODE_ARRAY, NODE, label,
                  init_node, compar_node_label)

FUNC_SEARCH_ARRAY(search_element_label, ELEMENT_ARRAY, ELEMENT,
                  label, init_element, compar_element_label)

FUNC_SEARCH_ARRAY(search_rigid_element_label, RIGID_ELEMENT_ARRAY,
                  RIGID_ELEMENT, label, init_rigid_element,
                  compar_rigid_element_label)

FUNC_SEARCH_ARRAY(search_material_prop_label, MATERIAL_PROP_ARRAY,
                  MATERIAL_PROP, label, init_material_prop,
                  compar_material_prop_label)

FUNC_SEARCH_ARRAY(search_physical_prop_label, PHYSICAL_PROP_ARRAY,
                  PHYSICAL_PROP, label, init_physical_prop,
                  compar_physical_prop_label)

FUNC_SEARCH_ARRAY(search_bc_node, BC_ARRAY, BC, node,
                  init_bc, compar_bc_node)

FUNC_SEARCH_ARRAY(search_elem_bc_element, ELEM_BC_ARRAY, ELEM_BC,
                  element, init_elem_bc, compar_elem_bc_element)

FUNC_SEARCH_ARRAY(search_disp_node, DISP_ARRAY, DISP_VELO,
                  node, init_disp_velo, compar_disp_velo_node)

FUNC_SEARCH_ARRAY(search_velo_node, VELO_ARRAY, DISP_VELO,
                  node, init_disp_velo, compar_disp_velo_node)

FUNC_SEARCH_ARRAY(search_react_node, REACT_ARRAY, DISP_VELO,
                  node, init_disp_velo, compar_disp_velo_node)

FUNC_SEARCH_ARRAY(search_stress_element, STRESS_ARRAY,
                  STRESS_STRAIN, element, init_stress_strain,
                  compar_stress_strain_element)

FUNC_SEARCH_ARRAY(search_strain_element, STRAIN_ARRAY,
                  STRESS_STRAIN, element, init_stress_strain,
                  compar_stress_strain_element)

FUNC_SEARCH_ARRAY(search_sensitivity_element, SENSITIVITY_ARRAY,
                  SENSITIVITY, element, init_sensitivity,
                  compar_sensitivity_element)

FUNC_SEARCH_ARRAY(search_default_bc_set_id, DEFAULT_BC_ARRAY,
                  DEFAULT_BC, set_id, init_default_bc,
                  compar_default_bc_set_id)

FUNC_SEARCH_ARRAY(search_sound_pressure_node, SOUND_PRESSURE_ARRAY,
                  SOUND_PRESSURE, node, init_sound_pressure,
                  compar_sound_pressure_node)

FUNC_SEARCH_ARRAY(search_bc_set_label, BC_SET_ARRAY,
                  BC_SET, label, init_bc_set,
                  compar_bc_set_label)

FUNC_SEARCH_ARRAY(search_elem_matrix_element, ELEM_MATRIX_ARRAY,
                  ELEM_MATRIX, element, init_elem_matrix,
                  compar_elem_matrix_element)

FUNC_SEARCH_ARRAY(search_elem_lookup_node, ELEM_LOOKUP_ARRAY,
                  ELEM_LOOKUP, node, init_elem_lookup,
                  compar_elem_lookup_node)

FUNC_SEARCH_ARRAY(search_node_matrix_node, NODE_MATRIX_ARRAY,
                  NODE_MATRIX, node, init_node_matrix,
                  compar_node_matrix_node)

FUNC_SEARCH_ARRAY(search_load_curve_label, LOAD_CURVE_ARRAY,
                  LOAD_CURVE, label, init_load_curve,
                  compar_load_curve_label)

FUNC_SEARCH_ARRAY(search_local_coord_label, LOCAL_COORD_ARRAY,
                  LOCAL_COORD, label, init_local_coord,
                  compar_local_coord_label)

FUNC_SEARCH_ARRAY(search_contact_node, CONTACT_ARRAY, CONTACT,
                  node, init_contact, compar_contact_node)


/*
 */
int
search_stress_element_node (STRESS_ARRAY stress, int element, int node)
{
	int found_index;

	found_index = search_stress_element(stress, element);
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


int
search_strain_element_node (STRAIN_ARRAY strain, int element, int node)
{
	int found_index;

	found_index = search_strain_element(strain, element);
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


int
search_sensitivity_element_node (SENSITIVITY_ARRAY sens, int element, int node)
{
	int found_index;

	found_index = search_sensitivity_element(sens, element);
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
RC
add_disp_array (double weight, DISP_ARRAY disp, DISP_ARRAY *disp_total)
{
	int ii1, index, num;
	VECT3DR tmp;

	for(ii1=0; ii1<disp.size; ii1++){
		if(disp.array[ii1].node < 0) continue;

		index = search_disp_node(*disp_total, disp.array[ii1].node);
		if(index >= 0){
			tmp = mul_scalar_vect3dr(weight, disp.array[ii1].v);
			disp_total->array[index].v = add_vect3dr(tmp,
			                               disp_total->array[index].v);
			disp_total->array[index].s += weight * disp.array[ii1].s;
		}else{
			RC_TRY( realloc_disp_array(disp_total) );
			num = disp_total->size;
			disp_total->array[num] = disp.array[ii1];
			disp_total->array[num].v = mul_scalar_vect3dr(weight,
			                                            disp.array[ii1].v);
			disp_total->array[num].s = weight * disp.array[ii1].s;
			(disp_total->size)++;
			RC_TRY( clean_disp_array(disp_total) );
		}
	}

	return(NORMAL_RC);
}


/* sens の weight 倍を sens_total に付加 */
/* sens_total に存在しない要素は新規に追加する */
RC
add_sensitivity_array (double weight, SENSITIVITY_ARRAY sens,
                       SENSITIVITY_ARRAY *sens_total)
{
	int ii1, index, num;
	VECT3DR tmp;

	for(ii1=0; ii1<sens.size; ii1++){
		if(sens.array[ii1].element < 0) continue;

		index = search_sensitivity_element_node(*sens_total,
		                  sens.array[ii1].element, sens.array[ii1].node);
		if(index >= 0){
			tmp = mul_scalar_vect3dr(weight, sens.array[ii1].v);
			sens_total->array[index].v = add_vect3dr(tmp,
			                                 sens_total->array[index].v);
			sens_total->array[index].s += weight * sens.array[ii1].s;
		}else{
			RC_TRY( realloc_sensitivity_array(sens_total) );
			num = sens_total->size;
			sens_total->array[num] = sens.array[ii1];
			sens_total->array[num].v = mul_scalar_vect3dr(weight,
			                                              sens.array[ii1].v);
			sens_total->array[num].s = weight * sens.array[ii1].s;
			(sens_total->size)++;
			RC_TRY( clean_sensitivity_array(sens_total) );
		}
	}

	return(NORMAL_RC);
}


/* 要素 element に対応する材料 mat のインデックスを返却 */
int
search_material_prop_element (MATERIAL_PROP_ARRAY mat,
                              PHYSICAL_PROP_ARRAY phys,
                              ELEMENT *element)
{
	int phys_index;

	if((phys.source == FEM_NASTRAN)&&(element->material < 0)){
		/* NASTRAN の場合は、一度 phys を検索する必要有り */
		phys_index = search_physical_prop_label(phys, element->physical);
		if(phys_index < 0){
			return(phys_index);
		}

		/* 次回の検索を高速化する */
		element->material = phys.array[phys_index].i_info[0];
	}

	if(element->material < 0) return(-2);

	return(search_material_prop_label(mat, element->material));
}


/* 要素 element を構成する節点座標値のマトリックスを作成する． 
 * 
 * matrix[n][3] = [
 *   [p1.x, p1.y, p1.z],
 *   [p2.x, p2.y, p2.z],
 *   [p3.x, p3.y, p3.z],
 *   ...
 *   [pn.x, pn.y, pn.z],
 * ]
 *
 * なぜ element はポインタ？
 */
RC
node_location_matrix (ELEMENT *element, NODE_ARRAY node, double **matrix)
{
	int ii1;

	for(ii1=0; ii1<element->node_num; ii1++){
		int index = search_node_label(node, element->node[ii1]);
		if(index < 0){
			fprintf(stderr, "search_node_label < \n");
			return(SEEK_ERROR_RC);
		}

		matrix[ii1][0] = node.array[index].p.x;
		matrix[ii1][1] = node.array[index].p.y;
		matrix[ii1][2] = node.array[index].p.z;
	}

	return(NORMAL_RC);
}


/* 全体座標系節点座標値マトリックスを局所座標系節点座標値マトリックスへ */
RC
local_node_location (double L[][3], int node_num, double **node_location)
{
	int ii1;

	for(ii1=0; ii1<node_num; ii1++){
		mul_matrix33_vect(L, node_location[ii1]);
	}

	return(NORMAL_RC);
}


/* 局所座標系節点座標値マトリックスを全体座標系節点座標値マトリックスへ */
RC
global_node_location (double L[][3], int node_num, double **node_location)
{
	int ii1;

	for(ii1=0; ii1<node_num; ii1++){
		mul_matrix33t_vect(L, node_location[ii1]);
	}

	return(NORMAL_RC);
}

/* 要素 element と変位配列 disp から，disp_matrix を作成する．
 */
RC
fill_disp_matrix (ELEMENT element, DISP_ARRAY disp, double **matrix)
{
	int ii1;

	for(ii1=0; ii1<element.node_num; ii1++){
		int index;

		RC_NEG_CHK( index = search_disp_node(disp, element.node[ii1]) );

		matrix[ii1][0] = disp.array[index].v.x;
		matrix[ii1][1] = disp.array[index].v.y;
		matrix[ii1][2] = disp.array[index].v.z;
		matrix[ii1][3] = disp.array[index].v.yz;
		matrix[ii1][4] = disp.array[index].v.zx;
		matrix[ii1][5] = disp.array[index].v.xy;
	}

	return(NORMAL_RC);
}


RC
fill_sound_pres_matrix (ELEMENT element, SOUND_PRESSURE_ARRAY pres,
                        double **matrix)
{
	int ii1;

	for(ii1=0; ii1<element.node_num; ii1++){
		int index;

		RC_NEG_CHK( index = search_sound_pressure_node(pres,
		                                    element.node[ii1]) );
		matrix[ii1][0] = pres.array[index].re;
		matrix[ii1][1] = pres.array[index].im;
	}

	return(NORMAL_RC);
}


RC
unit_normal_vect_center (ELEMENT surf, ELEMENT inner_elem,
                         NODE_ARRAY node, VECT3D *vect)
{
	double center_coord[4];

	RC_TRY( set_center_point(surf.type, center_coord) );
	RC_TRY( unit_normal_vect(center_coord, surf, inner_elem, node, vect) );

	return(NORMAL_RC);
}


/* 要素 surf 上の座標 local_coord における外向き単位法線ベクトルを計算する．
 * inner_elem の反対向きを外向きとする．
 */
RC
unit_normal_vect (double *local_coord, ELEMENT surf,
                  ELEMENT inner_elem, NODE_ARRAY node, VECT3D *vect)
{
	int ii1;
	int index;
	double s_coord[4];
	double i_coord[4];
	double s_N[MAX_NODE];
	double i_N[MAX_NODE];
	VECT3D s_center;
	VECT3D i_center;
	VECT3D d_xi, d_eta;
	VECT3D s_to_i;
	static int init_flag = 0;
	static double **dN = NULL;
	static double **node_location = NULL;
	static double **jacobi = NULL;
    static double **s_local_coord;
    static double **v_local_coord;
	double abs_vect;

	if(init_flag == 0){
		RC_TRY( allocate2D(3, MAX_NODE, &dN) );
		RC_TRY( allocate2D(MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(3, 3, &jacobi) );
        RC_TRY( allocate2D(1, 4, &s_local_coord) );
        RC_TRY( allocate2D(1, 4, &v_local_coord) );

		init_flag = 1;
	}

	/* surf の local_coord における法線ベクトルを計算する． */
	init_matrix(3, MAX_NODE, dN);   /* 初期化が必要 */
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
        /* 線要素の接線ベクトル */
		d_xi.x = jacobi[0][0];
		d_xi.y = jacobi[0][1];
		d_xi.z = jacobi[0][2];

        /* 平面要素の法線ベクトル */
        /** surf上の局所座標値をinner_elem上の局所座標値に変換する． **/
	    RC_TRY( node_location_matrix(&inner_elem, node, node_location) );
        RC_TRY( copy_vect(4, local_coord, s_local_coord[0]) );
        RC_TRY( trans_coord_s2v(surf, inner_elem, 1, s_local_coord, 
                    v_local_coord) );
        RC_TRY( cal_d_A(inner_elem, v_local_coord[0], node_location, &d_eta,
                   NULL, NULL) ); 

		*vect = outer_product3d(d_xi, d_eta);
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	/* 法線ベクトルを単位化 */
	abs_vect = abs_vect3d(*vect);
	if(nearly_eq(abs_vect, 0.0)) return(CAL_ERROR_RC);
	vect->x /= abs_vect;
	vect->y /= abs_vect;
	vect->z /= abs_vect;

	RC_TRY( set_center_point(surf.type, s_coord) );
	RC_TRY( set_N_vector(surf.type, s_coord, s_N) );
	s_center.x = 0.0;
	s_center.y = 0.0;
	s_center.z = 0.0;
	for(ii1=0; ii1<surf.node_num; ii1++){
		index = search_node_label(node, surf.node[ii1]);
		if(index < 0) return(SEEK_ERROR_RC);
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
		index = search_node_label(node, inner_elem.node[ii1]);
		if(index < 0) return(SEEK_ERROR_RC);
		i_center.x += i_N[ii1] * node.array[index].p.x;
		i_center.y += i_N[ii1] * node.array[index].p.y;
		i_center.z += i_N[ii1] * node.array[index].p.z;
	}

	/* 法線ベクトルの向きをチェック */
	s_to_i = sub_vect3d(i_center, s_center);
	if(nearly_eq(abs_vect3d(s_to_i), 0.0)) return(CAL_ERROR_RC);
	if(inner_product3d(*vect, s_to_i) > 0.0){
		vect->x *= -1.0;
		vect->y *= -1.0;
		vect->z *= -1.0;
	}

	return(NORMAL_RC);
}


RC
extract_rigid_set (int num, const int id[], double scale[],
                   RIGID_ELEMENT_ARRAY src, RIGID_ELEMENT_ARRAY *dest)
{
	int ii1, ii2;
	int set_flag;

	RC_TRY( allocate_rigid_element_array(0, dest) );
	dest->source = src.source;
	
	for(ii1=0; ii1<src.size; ii1++){
		if(src.array[ii1].label < 0) continue;

		set_flag = 0;
		if(src.array[ii1].set_id < 0) set_flag = 1;
		for(ii2=0; ii2<num; ii2++){
			if(src.array[ii1].set_id == id[ii2]) set_flag = 1;
		}

		if(set_flag){
			RC_TRY( realloc_rigid_element_array(dest) );
			RC_TRY( copy_rigid_element(src.array[ii1],
			                               &(dest->array[dest->size])) );
			(dest->size)++;
		}
	}
	RC_TRY( clean_rigid_element_array(dest) );

	return(NORMAL_RC);
}


/* 境界条件 src 内より、id[num] に該当する要素を */
/* scale[] 倍して dest にコピー                  */
RC
extract_bc_set (int num, const int id[], double scale[],
                BC_ARRAY src, BC_ARRAY *dest)
{
	int ii1, ii2, ii3;

	RC_TRY( allocate_bc_array(0, dest) );
	dest->source = src.source;
	
	for(ii1=0; ii1<src.size; ii1++){
		if(src.array[ii1].node < 0) continue;

		for(ii2=0; ii2<num; ii2++){
			if(src.array[ii1].set_id == id[ii2]){
				RC_TRY( realloc_bc_array(dest) );
				dest->array[dest->size] = src.array[ii1];
				dest->array[dest->size].v
					= mul_scalar_vect3dr(scale[ii2], dest->array[dest->size].v);
				for(ii3=0; ii3<BC_SCALAR; ii3++){
					dest->array[dest->size].s[ii3] *= scale[ii2];
				}
				(dest->size)++;
				break;
			}
		}
	}
	RC_TRY( compress_bc_array(dest) );
	RC_TRY( clean_bc_array(dest) );

	return(NORMAL_RC);
}


RC
compress_bc_array (BC_ARRAY *dest)
{
	int ii1, ii2, ii3;

	/* 重複する節点をマージ */
	for(ii1=0; ii1<dest->size;){
		if(dest->array[ii1].node < 0) continue;
		for(ii2=ii1+1; ii2<dest->size; ii2++){
			if(dest->array[ii1].node != dest->array[ii2].node) break;
			RC_TRY( overwrite_bc(
			           &(dest->array[ii1].v_type.x), &(dest->array[ii1].v.x),
			             dest->array[ii2].v_type.x,    dest->array[ii2].v.x) );
			RC_TRY( overwrite_bc(
			           &(dest->array[ii1].v_type.y), &(dest->array[ii1].v.y),
			             dest->array[ii2].v_type.y,    dest->array[ii2].v.y) );
			RC_TRY( overwrite_bc(
			           &(dest->array[ii1].v_type.z), &(dest->array[ii1].v.z),
			             dest->array[ii2].v_type.z,    dest->array[ii2].v.z) );
			RC_TRY( overwrite_bc(
			         &(dest->array[ii1].v_type.yz), &(dest->array[ii1].v.yz),
			           dest->array[ii2].v_type.yz,    dest->array[ii2].v.yz) );
			RC_TRY( overwrite_bc(
			         &(dest->array[ii1].v_type.zx), &(dest->array[ii1].v.zx),
			           dest->array[ii2].v_type.zx,    dest->array[ii2].v.zx) );
			RC_TRY( overwrite_bc(
			         &(dest->array[ii1].v_type.xy), &(dest->array[ii1].v.xy),
			           dest->array[ii2].v_type.xy,    dest->array[ii2].v.xy) );
			for(ii3=0; ii3<BC_SCALAR; ii3++){
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
		for(ii2=0; ii2<BC_SCALAR; ii2++){
			if(dest->array[ii1].s_type[ii2] != BC_FREE) flag = 1;
		}
		if(flag == 0){
			dest->array[ii1].node = -1;
		}
	}

	RC_TRY( clean_bc_array(dest) );

	return(NORMAL_RC);
}


static RC
overwrite_bc (BC_TYPE *type1, double *v1, BC_TYPE type2,  double v2)
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


RC
extract_default_bc_set (int num, const int id[], double scale[],
                        DEFAULT_BC_ARRAY src, DEFAULT_BC_ARRAY *dest)
{
	int ii1, ii2;

	RC_TRY( allocate_default_bc_array(0, dest) );
	dest->source = src.source;
	
	for(ii1=0; ii1<src.size; ii1++){
		if(src.array[ii1].set_id < 0) continue;

		for(ii2=0; ii2<num; ii2++){
			if(src.array[ii1].set_id == id[ii2]){
				RC_TRY( realloc_default_bc_array(dest) );
				dest->array[dest->size] = src.array[ii1];
				dest->array[dest->size].temp *= scale[ii2];
				dest->array[dest->size].grav
				  = mul_scalar_vect3d(scale[ii2], dest->array[dest->size].grav);
				(dest->size)++;
				break;
			}
		}
	}

	RC_TRY( clean_default_bc_array(dest) );

	return(NORMAL_RC);
}


RC
extract_elem_bc_set (int num, const int id[], double scale[],
                     ELEM_BC_ARRAY src, ELEM_BC_ARRAY *dest)
{
	int ii1, ii2, ii3;

	RC_TRY( allocate_elem_bc_array(0, dest) );
	dest->source = src.source;
	
	for(ii1=0; ii1<src.size; ii1++){
		if(src.array[ii1].element < 0) continue;

		for(ii2=0; ii2<num; ii2++){
			if(src.array[ii1].set_id == id[ii2]){
				RC_TRY( realloc_elem_bc_array(dest) );
				dest->array[dest->size] = src.array[ii1];
				if( (dest->array[dest->size].type == BC_PRESS)
				  ||(dest->array[dest->size].type == BC_TRACTION) ){
					dest->array[dest->size].val.s *= scale[ii2];
				}else if( (dest->array[dest->size].type == BC_PRESS_D1)
				        ||(dest->array[dest->size].type == BC_TRACTION_D1)
				        ||(dest->array[dest->size].type == BC_PRESS_D2)
				        ||(dest->array[dest->size].type == BC_TRACTION_D2) ){
					for(ii3=0; ii3<MAX_NODE; ii3++){
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

	RC_TRY( clean_elem_bc_array(dest) );

	return(NORMAL_RC);
}


/* elem_matrixより作成される全体マトリックスと vector[total_dof] の掛け算 */
/* 結果は answer[total_dof] に格納される                                  */
/* 可能な限りエラーチェックを入れた低速版                                 */
RC
mul_elem_matrices_vector (int total_dof, ELEM_MATRIX_ARRAY elem_matrix,
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
/* vectors[][m] への代入は 上位関数が行う必要あり                        */
RC
modify_elem_matrix_n (ELEM_MATRIX elem_matrix, int m, double v,
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


RC
make_elem_lookup (NODE_ARRAY node, ELEMENT_ARRAY element,
                  ELEM_LOOKUP_ARRAY *elem_lookup)
{
    int ii1, ii2;

	/* elem_lookup の確保と初期化 */
	RC_TRY( allocate_elem_lookup_array(node.size, elem_lookup) );
	for(ii1=0; ii1<node.size; ii1++){
		elem_lookup->array[ii1].node = node.array[ii1].label;
		elem_lookup->array[ii1].size = 0;
		elem_lookup->array[ii1].alloc_size = 0;
		elem_lookup->array[ii1].element = NULL;
		elem_lookup->array[ii1].r_size = 0;
		elem_lookup->array[ii1].r_alloc_size = 0;
		elem_lookup->array[ii1].r_element = NULL;
	}
	elem_lookup->size = node.size;
	RC_TRY( clean_elem_lookup_array(elem_lookup) );

	/* 各テーブルのサイズを決定 */
	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
			int index = search_elem_lookup_node(*elem_lookup,
			                          element.array[ii1].node[ii2]);
			RC_NEG_CHK(index);
			elem_lookup->array[index].alloc_size++;
		}
	}

	/* 各テーブルを確保 */
	for(ii1=0; ii1<(elem_lookup->size); ii1++){
		if(elem_lookup->array[ii1].node < 0){
            continue;
        }
		if(elem_lookup->array[ii1].alloc_size > 0){
            RC_TRY( allocate1I(elem_lookup->array[ii1].alloc_size, 
                        &(elem_lookup->array[ii1].element)) );
		}
	}

	/* 各テーブルに代入 */
	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
			int size;
			int index = search_elem_lookup_node(*elem_lookup,
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
RC
make_elem_lookup_rigid (NODE_ARRAY node, ELEMENT_ARRAY element,
                        RIGID_ELEMENT_ARRAY rigid_elem,
                        ELEM_LOOKUP_ARRAY *elem_lookup)
{
	int ii1, ii2;

	/* elem_lookup の確保と初期化，element分の処理 */
	RC_TRY( make_elem_lookup(node, element, elem_lookup) );

	/* 各テーブルのサイズを決定 */
	for(ii1=0; ii1<rigid_elem.size; ii1++){
		if(rigid_elem.array[ii1].label < 0) continue;

		for(ii2=0; ii2<rigid_elem.array[ii1].node_num; ii2++){
			int index;
			index = search_elem_lookup_node(*elem_lookup,
			                          rigid_elem.array[ii1].node[ii2]);
			if(index < 0){
				RC_TRY( realloc_elem_lookup_array(elem_lookup) );
				elem_lookup->array[elem_lookup->size].node
				           = rigid_elem.array[ii1].node[ii2];
				(elem_lookup->size)++;
				RC_TRY( clean_elem_lookup_array(elem_lookup) );
				index = search_elem_lookup_node(*elem_lookup,
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
			    = mm_alloc(sizeof(int)*(elem_lookup->array[ii1].r_alloc_size));
			if(elem_lookup->array[ii1].r_element == NULL)
				return(ALLOC_ERROR_RC);
		}
	}

	/* 各テーブルに代入 */
	for(ii1=0; ii1<rigid_elem.size; ii1++){
		if(rigid_elem.array[ii1].label < 0) continue;

		for(ii2=0; ii2<rigid_elem.array[ii1].node_num; ii2++){
			int size;
			int index = search_elem_lookup_node(*elem_lookup,
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


/* 問題の次元                     */
/* 2...２次元問題、3...３次元問題 */
int
analysis_dim (ELEMENT_ARRAY element)
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


/* 要素の次数                         */
/* 1...１次要素、2...２次要素         */
/* 混在している場合は大きい次数を返却 */
int
analysis_order (ELEMENT_ARRAY element)
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
RC
reduced_element_array_order (ELEMENT_ARRAY *elem)
{
	ELEM_TYPE new_type;
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
RC
reduced_element_type (ELEM_TYPE old_type, ELEM_TYPE *new_type,
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
RC
delete_duplicate_node (ELEMENT_ARRAY element, NODE_ARRAY *node, double cr_dist)
{
	int ii1, ii2, ii3, ii4;
	ELEM_LOOKUP_ARRAY elem_lookup;
	double dist;
	int index, elem_index;

	RC_TRY( make_elem_lookup(*node, element, &elem_lookup) );

	for(ii1=0; ii1<node->size; ii1++){
		if(node->array[ii1].label < 0) continue;

		for(ii2=ii1+1; ii2<node->size; ii2++){
			if(node->array[ii2].label < 0) continue;

			dist = dist_point3d(node->array[ii1].p, node->array[ii2].p);
			if(dist >= cr_dist) continue;

			index = search_elem_lookup_node(elem_lookup,
			                                    node->array[ii2].label);
			for(ii3=0; ii3<elem_lookup.array[index].size; ii3++){
				elem_index = search_element_label(element,
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

	clean_node_array(node);
	RC_TRY( free_elem_lookup_array(&elem_lookup) );

	return(NORMAL_RC);
}


/* 節点間距離の最小値 */
RC
minimum_node_distance (ELEMENT_ARRAY element, NODE_ARRAY node, double *dist)
{
	int ii1, ii2, ii3;
	double dist_temp;
	int index0, index1;

	*dist = DBL_MAX;
	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;
		if(element.array[ii1].node_num <= 1) continue;

		for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
			index0 = search_node_label(node,
			                element.array[ii1].node[ii2]);
			RC_NEG_CHK(index0);

			for(ii3=0; ii3<element.array[ii1].node_num; ii3++){
				if(ii2 == ii3) continue;
				index1 = search_node_label(node, element.array[ii1].node[ii3]);
				RC_NEG_CHK(index1);
				dist_temp = dist_point3d(node.array[index0].p,
				                         node.array[index1].p);
				if(dist_temp < *dist) *dist = dist_temp;
			}
		}
	}

	return(NORMAL_RC);
}


/* elem.array[].node[] に使用されていない節点を消去する */
/* 消去後の node は clean_node_array() を実行される     */ 
RC
delete_unused_node (ELEMENT_ARRAY elem, NODE_ARRAY *node)
{
	int *flags;
	int index;
	int ii1, ii2;

	flags = (int *)mm_alloc( node->size * sizeof(int) );
	if(flags == NULL) return(ALLOC_ERROR_RC);

	for(ii1=0; ii1<node->size; ii1++){
		flags[ii1] = 0;
	}

	for(ii1=0; ii1<elem.size; ii1++){
		if(elem.array[ii1].label < 0) continue;
		for(ii2=0; ii2<elem.array[ii1].node_num; ii2++){
			RC_NEG_CHK( index = search_node_label(*node,
			                     elem.array[ii1].node[ii2]) );
			flags[index] = 1;
		}
	}
	for(ii1=0; ii1<node->size; ii1++){
		if(flags[ii1] == 0) node->array[ii1].label = -1;
	}

	RC_TRY( clean_node_array(node) );
	node->renum_flag = UNRENUMBERED;
	RC_TRY( mm_free((void *)flags) );

	return(NORMAL_RC);
}


RC
delete_unused_node_r (ELEMENT_ARRAY elem, RIGID_ELEMENT_ARRAY rigid,
                      NODE_ARRAY *node)
{
	int *flags;
	int index;
	int ii1, ii2;

	flags = (int *)mm_alloc( node->size * sizeof(int) );
	if(flags == NULL) return(ALLOC_ERROR_RC);
	for(ii1=0; ii1<node->size; ii1++){
		flags[ii1] = 0;
	}

	for(ii1=0; ii1<elem.size; ii1++){
		if(elem.array[ii1].label < 0) continue;
		for(ii2=0; ii2<elem.array[ii1].node_num; ii2++){
			RC_NEG_CHK( index = search_node_label(*node,
			                     elem.array[ii1].node[ii2]) );
			flags[index] = 1;
		}
	}

	for(ii1=0; ii1<rigid.size; ii1++){
		if(rigid.array[ii1].label < 0) continue;
		for(ii2=0; ii2<rigid.array[ii1].node_num; ii2++){
			RC_NEG_CHK( index = search_node_label(*node,
			                     rigid.array[ii1].node[ii2]) );
			flags[index] = 1;
		}
	}

	for(ii1=0; ii1<node->size; ii1++){
		if(flags[ii1] == 0) node->array[ii1].label = -1;
	}

	RC_TRY( clean_node_array(node) );
	node->renum_flag = UNRENUMBERED;
	RC_TRY( mm_free((void *)flags) );

	return(NORMAL_RC);
}


RC
delete_unused_bc (NODE_ARRAY node, BC_ARRAY *bc)
{
	int index;
	int ii1;

	for(ii1=0; ii1<bc->size; ii1++){
		if(bc->array[ii1].node < 0) continue;
		index = search_node_label(node, bc->array[ii1].node);
		if(index < 0){
			bc->array[ii1].node = -1;
		}
	}

	RC_TRY( clean_bc_array(bc) );

	return(NORMAL_RC);
}


/* 中間節点を両端の中点に移動させる */
RC
locate_middle_node (NODE_ARRAY node, ELEMENT_ARRAY elem)
{
	int ii1, ii2;
	int table[MAX_NODE][2];
	int n_index, index0, index1;

	for(ii1=0; ii1<elem.size; ii1++){
		if(elem.array[ii1].label < 0) continue;

		/* ２次要素以外は処理しない */
		if( element_order(elem.array[ii1].type) != 2) continue;

		RC_TRY( interpolate_table(elem.array[ii1].type, table) );
		for(ii2=0; ii2<elem.array[ii1].node_num; ii2++){
			RC_NEG_CHK( n_index = search_node_label(node,
			                       elem.array[ii1].node[ii2]) );
			RC_NEG_CHK( index0 = search_node_label(node,
			                       elem.array[ii1].node[table[ii2][0]]) );
			RC_NEG_CHK( index1 = search_node_label(node,
			                       elem.array[ii1].node[table[ii2][1]]) );
			node.array[n_index].p = mid_point3d(node.array[index0].p,
			                                    node.array[index1].p);
		}
	}

	return(NORMAL_RC);
}


RC
interpolate_table (ELEM_TYPE type, int table[][2])
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
RC
element_average_length(ELEMENT element, NODE_ARRAY node, double *length)
{
	int table[12][2]; /* 12->HEXA2最大edge数 */
	int *index;
	int edge_num;
	int ii1;

	RC_TRY( allocate1I(element.node_num, &index) );
	for(ii1=0; ii1<element.node_num; ii1++){
		index[ii1] = search_node_label(node, element.node[ii1]);
		RC_NEG_CHK(index[ii1]);
	}

	RC_TRY( element_edge_table(element.type, &edge_num, table) );

	*length = 0.0;
	for(ii1=0; ii1<edge_num; ii1++){
		*length += dist_point3d(node.array[ index[ table[ii1][0] ] ].p,
		                        node.array[ index[ table[ii1][1] ] ].p);
	}

	*length /= (double)edge_num;

	RC_TRY( mm_free((void *)index) );

	return(NORMAL_RC);
}


/* １要素の最大長さを計算 */
/* 具体的には edge の長さ */
double
element_length(int node_num, double **node_location)
{
	int ii1, ii2;
	double ret = ABS_TOL;

	for(ii1=0; ii1<node_num; ii1++){
		for(ii2=ii1+1; ii2<node_num; ii2++){
			double len = dist_point3d( array2vect3d(node_location[ii1]),
			                           array2vect3d(node_location[ii2]) );
			if(len > ret) ret = len;
		}
	}

	return(ret);
}


/* ２次要素の中間節点は考慮していない(直線扱い）*/
/* 中間節点の扱いは検討中                       */
RC
element_edge_table (ELEM_TYPE type, int *edge_num, int table[][2])
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


/*
 * 連立方程式の解を、要素マトリックスの情報で求める．
 * solution[] に近似解(初期値)を代入すると、処理後は厳密解が代入され
 * ている初期値と厳密解の差や、全体マトリックスの値や大きさによって、
 * ループ回数は変化する
 */
RC
elem_matrix_cg (int total_dof, ELEM_MATRIX_ARRAY elem_matrix,
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
		/* 相対誤差の計算                                     */
		relative_error = sqrt( (r_r / vector_value) / total_dof );
		RC_TRY( log_printf(5, "%d:%15.7e:%15.7e", counter,
		                   relative_error, REL_ERROR) );
		if( relative_error <= REL_ERROR ) break;

		if( fabs(p_ar) <= ABS_ERROR ) return (CAL_ERROR_RC);
		alpha = r_r / p_ar;
		RC_TRY( log_printf(5, "%15.7e\n", alpha) );

		/* solution[] の修正 */
		/* r[] の更新        */
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
		if( fabs(r_r) <= ABS_ERROR ) return (CAL_ERROR_RC);
		beta = r1_r1 / r_r;

		/* r_r を更新 */
		r_r = r1_r1;
		
		/* 次回の修正方向の計算 */
		for(ii1=0; ii1<total_dof; ii1++){
			p[ii1] = r[ii1] + beta * p[ii1];
		}

		counter ++;
	}

	RC_TRY( mm_free((void*)r) );
	RC_TRY( mm_free((void*)p) );
	RC_TRY( mm_free((void*)ar) );
	RC_TRY( mm_free((void*)ap) );

	return(NORMAL_RC);
}


/* elemを構成している節点ラベルをソートして抽出 */
RC
get_element_construct_node(int *num, int **node_label, ELEMENT_ARRAY elem)
{
	int ii1, ii2;
	int size;

	if( (num == NULL) || (node_label == NULL) ) return(ARG_ERROR_RC);
	if(elem.size <= 0){
		*num = 0;
		return(NORMAL_RC);
	}

	size = elem.size * 2;
	*node_label = (int *)mm_alloc(size * sizeof(int) );

	for(ii1=0, *num=0; ii1<elem.size; ii1++){
		if(elem.array[ii1].label < 0) continue;
		for(ii2=0; ii2<elem.array[ii1].node_num; ii2++){
			(*node_label)[(*num)++] = elem.array[ii1].node[ii2];
			if(*num == size){
				size += MAX_NODE*10;
				RC_NULL_CHK( *node_label = (int *)mm_realloc(*node_label,
				                                        size * sizeof(int) ) );
			}
		}
	}
	RC_NULL_CHK( *node_label = (int *)mm_realloc(*node_label,
	                                             *num * sizeof(int) ) );
	RC_TRY( sort_int_array(*num, *node_label) );
	RC_TRY( uniq_int_array(num, node_label) );

	return(NORMAL_RC);
}

