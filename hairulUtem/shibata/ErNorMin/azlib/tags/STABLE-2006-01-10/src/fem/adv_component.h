/*********************************************************************
 * adv_component.h
 *
 * Copyright (C) 2001 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Yasuo ASAGA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: adv_component.h 410 2005-06-27 07:16:06Z sasaoka $ */

#ifndef ADV_COMPONENT_H
#define ADV_COMPONENT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rc.h"
#include "fem_struct.h"
#include "Adv/AdvDocument.h"


RC adv_input_model(const char *fname, FEM_NODE_ARRAY *node,
                   FEM_ELEMENT_ARRAY *elem, FEM_BC_ARRAY *rest,
                   FEM_BC_ARRAY *force, FEM_MATERIAL_PROP_ARRAY *material,
                   FEM_PHYSICAL_PROP_ARRAY *physical);
RC adv_input_disp(const char *adv_model, const char *adv_result,
                  FEM_DISP_ARRAY *disp);
RC adv_input_stress(const char *adv_model, const char *adv_result,
                    FEM_STRESS_ARRAY *stress);
RC adv_input_strain(const char *adv_model, const char *adv_result,
                    FEM_STRAIN_ARRAY *strain);
RC adv_input_dens_ratio(AdvDatabox *dbox, FEM_MATERIAL_PROP_ARRAY *dens);
RC adv_input_element(AdvDatabox *dbox, FEM_ELEMENT_ARRAY *elem);
RC adv_input_node(AdvDatabox *dbox, FEM_NODE_ARRAY *node);
RC adv_input_restraint(AdvDatabox *dbox, FEM_BC_ARRAY *rest);
RC adv_input_force(AdvDatabox *dbox, FEM_BC_ARRAY *force);
RC adv_input_material(AdvDatabox *dbox, FEM_MATERIAL_PROP_ARRAY *material);
RC adv_input_physical(AdvDatabox *dbox, FEM_PHYSICAL_PROP_ARRAY *physical);

RC adv_output_model(char *fname, FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY elem,
                    FEM_BC_ARRAY rest, FEM_BC_ARRAY force,
                    FEM_MATERIAL_PROP_ARRAY material);
RC adv_output_element(AdvDocFile *dfile, FEM_ELEMENT_ARRAY elem,
                      FEM_NODE_ARRAY node);
RC adv_output_node(AdvDocFile *dfile, FEM_NODE_ARRAY node);
RC adv_output_restraint(AdvDocFile *dfile, FEM_BC_ARRAY rest,
                        FEM_NODE_ARRAY node);
RC adv_output_force(AdvDocFile *dfile, FEM_BC_ARRAY force, FEM_NODE_ARRAY node);
RC adv_output_material(AdvDocFile *dfile, FEM_MATERIAL_PROP_ARRAY material,
                       FEM_ELEMENT_ARRAY elem);
RC adv_output_dens_ratio(AdvDocFile *dfile, FEM_MATERIAL_PROP_ARRAY dens);
RC adv_output_dens_ratio_to_every_part(char *adv_model, char *adv_result, 
                                       FEM_ELEMENT_ARRAY elem,
                                       FEM_MATERIAL_PROP_ARRAY dens);

RC adv_fem(const char *adv_model, const char *adv_result,
           const char *adv_metis, const char *adv_solid,
           FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY elem,
           FEM_MATERIAL_PROP_ARRAY material, FEM_BC_ARRAY rest,
           FEM_BC_ARRAY force, FEM_DISP_ARRAY *disp, FEM_STRESS_ARRAY *stress,
           FEM_STRAIN_ARRAY *strain);

RC adv_speed_analysis(const char *adv_model, const char *adv_result,
                      const char *adv_metis, const char *adv_solid,
                      int num_case, FEM_NODE_ARRAY node,
                      FEM_ELEMENT_ARRAY elem, FEM_MATERIAL_PROP_ARRAY material,
                      FEM_BC_ARRAY rest, const FEM_BC_ARRAY *trac,
                      FEM_DISP_ARRAY *disp, FEM_STRAIN_ARRAY *strain);

#endif /* ADV_COMPONENT_H */
