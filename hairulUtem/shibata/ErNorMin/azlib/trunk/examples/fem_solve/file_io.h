#ifndef FILE_IO_H
#define FILE_IO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fem_struct.h"
#include "nst_component.h"
#include "rc.h"

RC input_bdf(FILE *fp, NODE_ARRAY *node, ELEMENT_ARRAY *element,
             MATERIAL_PROP_ARRAY *material, PHYSICAL_PROP_ARRAY *physical,
             NST_PARAM_ARRAY *bdf_param, BC_SET_ARRAY *bc_set, BC_ARRAY *rest,
             BC_ARRAY *force, DEFAULT_BC_ARRAY *b_force, BC_ARRAY *shape_rest,
             BC_ARRAY *temp, DEFAULT_BC_ARRAY *def_temp);

RC output_bdf(FILE *fpw, NODE_ARRAY node, ELEMENT_ARRAY element,
              MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
              NST_PARAM_ARRAY bdf_param, BC_SET_ARRAY bc_set, BC_ARRAY rest,
              BC_ARRAY force, DEFAULT_BC_ARRAY b_force, BC_ARRAY shape_rest);

RC output_dxf(FILE *fpw, NODE_ARRAY node, ELEMENT_ARRAY element);

#endif /* FILE_IO_H */



