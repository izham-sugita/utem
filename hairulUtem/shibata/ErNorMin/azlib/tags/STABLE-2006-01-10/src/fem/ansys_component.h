/*********************************************************************
 * ansys_component.h
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

/* $Id: ansys_component.h 300 2004-06-09 11:00:54Z sasaoka $ */

#ifndef ANSYS_COMPONENT_H
#define	ANSYS_COMPONENT_H

#include <stdio.h>
#include <limits.h>
#include <float.h>
#include "fem_struct.h"
#include "rc.h"
#include "string_utl.h"


#define ANSYSDEFAULT  (INT_MIN)
#define FANSYSDEFAULT (-(DBL_MAX/2.0))


RC ansys_input_node( FILE *fp, NODE_ARRAY *node );
RC ansys_output_node( FILE *fp, NODE_ARRAY node );
RC ansys_input_element( FILE *fp, ELEMENT_ARRAY *elem );
RC ansys_output_element( FILE *fp, ELEMENT_ARRAY elem );
RC ansys_input_force( FILE *fp, BC_ARRAY *bc );
RC ansys_output_force( FILE *fp, BC_ARRAY bc );
RC ansys_input_bc_displacement( FILE *fp, BC_ARRAY *bc );
RC ansys_input_bc_temperature( FILE *fp, BC_ARRAY *bc );
RC ansys_output_bc( FILE *fp, BC_ARRAY bc );
RC ansys_input_restraint(FILE *fp, BC_ARRAY *rest);
RC ansys_input_material( FILE *fp, MATERIAL_PROP_ARRAY *material );
RC ansys_output_material( FILE *fp, MATERIAL_PROP_ARRAY material );
RC ansys_output_cdb_file( FILE *fp, 
                          NODE_ARRAY node,
                          ELEMENT_ARRAY element,
                          MATERIAL_PROP_ARRAY material,
                          BC_ARRAY bc,
                          BC_ARRAY force );
RC ansys_output_model_cdb( FILE *fp,
                          NODE_ARRAY node,
                          ELEMENT_ARRAY element,
                          BC_ARRAY bc);

RC ansys_read_vp_result(FILE *fp_velo, FILE *fp_press, VELO_ARRAY *vp);

RC ansys_input_bc( FILE *fp, BC_ARRAY *bc );


#endif /* ANSYS_COMPONENT_H */

/* 
memo
*/






