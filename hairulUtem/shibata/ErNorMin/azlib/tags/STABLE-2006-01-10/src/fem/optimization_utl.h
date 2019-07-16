/*********************************************************************
 * optimization_utl.h
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Tomoyuki TSUBATA> <Takaaki NAGATANI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: optimization_utl.h 343 2004-07-16 01:54:46Z sasaoka $ */

#ifndef OPTIMIZATION_UTL_H
#define OPTIMIZATION_UTL_H

#include <stdio.h>
#include "rc.h"
#include "fem_struct.h"

#define IGNORE_MAT_LABEL    (900000)
#define NONDESIGN_MAT_LABEL (800000)

RC velocity_analysis(NODE_ARRAY node, ELEMENT_ARRAY element,
                     MATERIAL_PROP_ARRAY material,
                     PHYSICAL_PROP_ARRAY physical,
                     BC_ARRAY shape_rest, int num_case,
                     const BC_ARRAY tractions[], DISP_ARRAY disps[],
                     int solver_flag);

RC functional_variation(BC_ARRAY sens_force,
                        DISP_ARRAY disp, double *delta_f);

RC volume_functional(NODE_ARRAY node, ELEMENT_ARRAY element,
                     double *volume);
RC compliance_functional(BC_ARRAY force, DISP_ARRAY disp,
                         double *compliance);
RC potential_functional(BC_ARRAY force, DISP_ARRAY disp,
                        BC_ARRAY rest, REACT_ARRAY react,
                        double *potential);

RC dissipation_functional(double fluid_param, NODE_ARRAY node,
                          ELEMENT_ARRAY element,
                          STRAIN_ARRAY strain_velo, double *dissipation);

RC dens_sens_comp_disp(const double power, DISP_ARRAY disp,
                       NODE_ARRAY node, ELEMENT_ARRAY element,
                       MATERIAL_PROP_ARRAY material,
                       PHYSICAL_PROP_ARRAY physical,
                       SENSITIVITY_ARRAY *sens,
                       MATERIAL_PROP_ARRAY dens);
RC shape_sens_compliance_disp(ELEMENT_ARRAY surf, DISP_ARRAY disp,
                              NODE_ARRAY node, ELEMENT_ARRAY element,
                              MATERIAL_PROP_ARRAY material,
                              PHYSICAL_PROP_ARRAY physical,
                              SENSITIVITY_ARRAY *sens);
RC shape_sens_compliance_ss(ELEMENT_ARRAY surf,
                            STRESS_ARRAY stress, STRAIN_ARRAY strain,
                            SENSITIVITY_ARRAY *sens);
RC shape_sens_dissipation(double viscosity, ELEMENT_ARRAY surf,
                          VELO_ARRAY flow, VELO_ARRAY adjflow,
                          STRAIN_ARRAY strain, STRAIN_ARRAY adjstrain,
                          SENSITIVITY_ARRAY *sens);
RC shape_sens_volume(ELEMENT_ARRAY surf, SENSITIVITY_ARRAY *sens);

RC sens2traction(NODE_ARRAY node, ELEMENT_ARRAY surf,
                 BC_ARRAY shape_rest, SENSITIVITY_ARRAY sens,
                 BC_ARRAY *trac, int reduce_flag);

RC add_disp_node(NODE_ARRAY node, int disp_num, const double *factors,
                 const DISP_ARRAY *disps);

RC max_pr_strain(int str_num, const double *factors,
                 const STRAIN_ARRAY *strain, double *max_str);

RC max_pr_strain_disp(int str_num, const double factors[],
                      NODE_ARRAY node, ELEMENT_ARRAY element,
                      const DISP_ARRAY disp[], double *max_str,
                      int sol_ss_flag);

RC make_surface_stabilizer(NODE_ARRAY node, ELEMENT_ARRAY surf,
                           double ****matrix);

RC make_trac_material_const(MATERIAL_PROP_ARRAY material,
                            MATERIAL_PROP_ARRAY *material_trac);

RC make_trac_material(MATERIAL_PROP_ARRAY material,
                      MATERIAL_PROP_ARRAY *material_trac);

RC max_mises_stress(NODE_ARRAY node, ELEMENT_ARRAY element,
                    MATERIAL_PROP_ARRAY material,
                    PHYSICAL_PROP_ARRAY physical,
                    DISP_ARRAY disp, double *max_val);
RC max_mises_stress_therm(NODE_ARRAY node, ELEMENT_ARRAY element,
                          MATERIAL_PROP_ARRAY material,
                          PHYSICAL_PROP_ARRAY physical,
                          BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                          DISP_ARRAY disp, double *max_val);

RC mises_integration(NODE_ARRAY node, ELEMENT_ARRAY element,
                     MATERIAL_PROP_ARRAY material,
                     PHYSICAL_PROP_ARRAY physical,
                     DISP_ARRAY disp, double rho, double *integ_val,
                     double *volume);

RC max_mises_adj_force(NODE_ARRAY node, ELEMENT_ARRAY element,
                       MATERIAL_PROP_ARRAY material,
                       PHYSICAL_PROP_ARRAY physical, DISP_ARRAY disp,
                       double rho, BC_ARRAY *adj_force);

RC shape_sens_max_mises(ELEMENT_ARRAY surf, DISP_ARRAY disp,
                        DISP_ARRAY adj_disp,
                        NODE_ARRAY node, ELEMENT_ARRAY element,
                        MATERIAL_PROP_ARRAY material,
                        PHYSICAL_PROP_ARRAY physical,
                        double rho,
                        SENSITIVITY_ARRAY *sens);
RC shape_trac_max_mises(ELEMENT_ARRAY surf, DISP_ARRAY disp,
                        DISP_ARRAY adj_disp,
                        NODE_ARRAY node, ELEMENT_ARRAY element,
                        MATERIAL_PROP_ARRAY material,
                        PHYSICAL_PROP_ARRAY physical,
                        double rho, BC_ARRAY *sens_trac, int reduce_flag);
RC shape_trac_compliance(ELEMENT_ARRAY surf, DISP_ARRAY disp,
                         NODE_ARRAY node, ELEMENT_ARRAY element,
                         MATERIAL_PROP_ARRAY material,
                         PHYSICAL_PROP_ARRAY physical,
                         BC_ARRAY *sens_trac, int reduce_flag);
RC shape_trac_volume(ELEMENT_ARRAY surf, NODE_ARRAY node,
                     BC_ARRAY *sens_trac, int reduce_flag);

RC shape_sens_therm(ELEMENT_ARRAY surf, DISP_ARRAY disp,
                    NODE_ARRAY node, ELEMENT_ARRAY element,
                    MATERIAL_PROP_ARRAY material,
                    PHYSICAL_PROP_ARRAY physical,
                    BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                    SENSITIVITY_ARRAY *sens);

RC add_shape_rest(ELEMENT_ARRAY element, MATERIAL_PROP_ARRAY material,
                  PHYSICAL_PROP_ARRAY physical, BC_ARRAY *shape_rest);


#endif /* OPTIMIZATION_UTL_H */
