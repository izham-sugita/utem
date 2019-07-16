#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fem.h"
#include "base.h"
#include "file_io.h"


RC
input_bdf (FILE *fp, NODE_ARRAY *node, ELEMENT_ARRAY *element,
           MATERIAL_PROP_ARRAY *material, PHYSICAL_PROP_ARRAY *physical,
           NST_PARAM_ARRAY *bdf_param, BC_SET_ARRAY *bc_set, BC_ARRAY *rest,
           BC_ARRAY *force, DEFAULT_BC_ARRAY *b_force, BC_ARRAY *shape_rest,
           BC_ARRAY *temp, DEFAULT_BC_ARRAY *def_temp)

{
	BC_ARRAY all_rest;
	ELEM_BC_ARRAY elem_temp;

	RC_TRY( log_printf(4, "---------- model data ----------\n") );

	RC_TRY( nst_input_node(fp, node) );
	RC_TRY( log_printf(4, "Number of nodes : %d\n", node->size) );

	RC_TRY( nst_input_element(fp, element) );
	RC_TRY( log_printf(4, "Number of elements : %d\n", element->size) );

	RC_TRY( nst_input_material_prop(fp, material) );
	RC_TRY( log_printf(4, "Number of materials : %d\n", material->size) );

	RC_TRY( nst_input_physical_prop(fp, physical) );
	RC_TRY( log_printf(4, "Number of physical properties : %d\n",
	                   physical->size) );

	RC_TRY( nst_input_param( fp, bdf_param) );
	RC_TRY( log_printf(4, "Number of bdf_param: %d\n", bdf_param->size) );

	RC_TRY( nst_input_bc_set(fp, bc_set) );
	RC_TRY( log_printf(4, "Number of bc_set : %d\n", bc_set->size) );

	RC_TRY( nst_input_force(fp, force, b_force) );
	RC_TRY( log_printf(4, "Number of forces : %d\n", force->size) );
	RC_TRY( log_printf(4, "Number of b_force : %d\n", b_force->size) );

	RC_TRY( nst_input_restraint(fp, &all_rest) );
	RC_TRY( extract_bc_set(bc_set->array[0].num_fix ,
	                       bc_set->array[0].id_fix ,
	                       bc_set->array[0].scale_fix ,
	                       all_rest, rest ) );
	RC_TRY( log_printf(4, "Number of restraints : %d\n", rest->size) );

	RC_TRY( extract_bc_set(bc_set->array[bc_set->size-1].num_fix ,
	                       bc_set->array[bc_set->size-1].id_fix ,
	                       bc_set->array[bc_set->size-1].scale_fix ,
	                       all_rest, shape_rest ) );
	RC_TRY( log_printf(4, "Number of shape restraints : %d\n",
	                   shape_rest->size) );

	RC_TRY( nst_input_temperature(fp, temp, &elem_temp, def_temp) );

	RC_TRY( free_bc_array( &all_rest ) );
	RC_TRY( free_elem_bc_array( &elem_temp ) );


	return(NORMAL_RC);
}


RC
output_bdf (FILE *fpw, NODE_ARRAY node, ELEMENT_ARRAY element,
            MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
            NST_PARAM_ARRAY bdf_param, BC_SET_ARRAY bc_set, BC_ARRAY rest,
            BC_ARRAY force, DEFAULT_BC_ARRAY b_force, BC_ARRAY shape_rest)
{
	fprintf(fpw, "SOL 106\n");
	fprintf(fpw, "CEND\n");
	RC_TRY( nst_write_case_control( fpw , bc_set ) );
	RC_TRY( nst_output_param( fpw, bdf_param) );
	fprintf(fpw, "NLPARM   1       2               ITER    1       25      "
	             "        NO      \n");
	RC_TRY( nst_output_physical_prop(fpw, physical) );
	RC_TRY( nst_output_node(fpw, node) );
	RC_TRY( nst_output_element(fpw, element) );
	/*
	RC_TRY( nst_output_material_prop(fpw, material) );
	*/
	RC_TRY( nst_output_bc_set( fpw , bc_set ) );
	RC_TRY( nst_output_restraint(fpw, rest) );
	RC_TRY( nst_output_force(fpw, force, b_force) );
	fprintf(fpw, "$---------- shape_rest ----------\n");
	RC_TRY( nst_output_restraint(fpw, shape_rest) );
	fprintf(fpw, "ENDDATA\n");

	return(NORMAL_RC);
}


RC
output_dxf (FILE *fpw, NODE_ARRAY node, ELEMENT_ARRAY element)
{
	ELEMENT_ARRAY surf;

	switch( analysis_dim(element) ){
	case 3:
		RC_TRY( extract_surface(element, &surf) );
		RC_TRY( output_dxf_polyline(fpw, surf, node) );
		/* RC_TRY( output_dxf_3dface(fpw, "result", 0, surf, node) ); */

		RC_TRY( free_element_array( &surf ) );
		break;
	case 2:
		RC_TRY( output_dxf_polyline(fpw, element, node) );
		/* RC_TRY( output_dxf_3dface(fpw, "result", 0, element, node) ); */
		break;
	default:
		return(WRITE_ERROR_RC);
	}

	return(NORMAL_RC);
}

