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

/* $Id: optimization_utl.h,v 1.18 2004/01/08 11:31:23 nagatani Exp $ */

#ifndef OPTIMIZATION_UTL_H
#define OPTIMIZATION_UTL_H

#include <stdio.h>
#include "rc.h"
#include "fem_struct.h"

#define IGNORE_MAT_LABEL    (900000)
#define NONDESIGN_MAT_LABEL (800000)

RC speed_analysis(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                  FEM_MATERIAL_PROP_ARRAY material,
                  FEM_PHYSICAL_PROP_ARRAY physical,
                  FEM_BC_ARRAY shape_rest, int num_case,
                  const FEM_BC_ARRAY tractions[], FEM_DISP_ARRAY disps[],
                  int solver_flag);

RC functional_variation(FEM_BC_ARRAY sens_force,
                        FEM_DISP_ARRAY disp, double *delta_f);

RC volume_functional(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                     double *volume);
RC compliance_functional(FEM_BC_ARRAY force, FEM_DISP_ARRAY disp,
                         double *compliance);
RC potential_functional(FEM_BC_ARRAY force, FEM_DISP_ARRAY disp,
                        FEM_BC_ARRAY rest, FEM_REACT_ARRAY react,
                        double *potential);

RC dissipation_functional(double fluid_param, FEM_NODE_ARRAY node,
                          FEM_ELEMENT_ARRAY element,
                          FEM_STRAIN_ARRAY strain_velo, double *dissipation);

RC dens_sens_comp_disp(const double power, FEM_DISP_ARRAY disp,
                       FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                       FEM_MATERIAL_PROP_ARRAY material,
                       FEM_PHYSICAL_PROP_ARRAY physical,
                       FEM_SENSITIVITY_ARRAY *sens,
                       FEM_MATERIAL_PROP_ARRAY dens);
RC fem_local_stress_strain_therm(FEM_ELEMENT element,
                             const double *local_point,
                             FEM_NODE_ARRAY node, FEM_DISP_ARRAY disp,
                             FEM_MATERIAL_PROP_ARRAY material,
                             FEM_PHYSICAL_PROP_ARRAY physical,
                             FEM_BC_ARRAY temp, FEM_DEFAULT_BC_ARRAY def_temp,
                             double fact,
                             FEM_STRESS_STRAIN *stress,
                             FEM_STRESS_STRAIN *strain);
RC shape_sens_compliance_disp(FEM_ELEMENT_ARRAY surf, FEM_DISP_ARRAY disp,
                              FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                              FEM_MATERIAL_PROP_ARRAY material,
                              FEM_PHYSICAL_PROP_ARRAY physical,
                              FEM_SENSITIVITY_ARRAY *sens);
RC shape_sens_compliance_ss(FEM_ELEMENT_ARRAY surf,
                            FEM_STRESS_ARRAY stress, FEM_STRAIN_ARRAY strain,
                            FEM_SENSITIVITY_ARRAY *sens);
RC shape_sens_dissipation(double viscosity, FEM_ELEMENT_ARRAY surf,
                          FEM_VELO_ARRAY flow, FEM_VELO_ARRAY adjflow,
                          FEM_STRAIN_ARRAY strain, FEM_STRAIN_ARRAY adjstrain,
                          FEM_SENSITIVITY_ARRAY *sens);
RC shape_sens_volume(FEM_ELEMENT_ARRAY surf, FEM_SENSITIVITY_ARRAY *sens);

RC sens2traction(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY surf,
                 FEM_BC_ARRAY shape_rest, FEM_SENSITIVITY_ARRAY sens,
                 FEM_BC_ARRAY *trac, int reduce_flag);

RC add_disp_node(FEM_NODE_ARRAY node, int disp_num, const double *factors,
                 const FEM_DISP_ARRAY *disps);

RC max_pr_strain(int str_num, const double *factors,
                 const FEM_STRAIN_ARRAY *strain, double *max_str);

RC max_pr_strain_disp(int str_num, const double factors[],
                      FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                      const FEM_DISP_ARRAY disp[], double *max_str,
                      int fem_sol_ss_flag);

RC make_surface_stabilizer(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY surf,
                           double ****matrix);

RC make_trac_material_const(FEM_MATERIAL_PROP_ARRAY material,
                            FEM_MATERIAL_PROP_ARRAY *material_trac);

RC make_trac_material(FEM_MATERIAL_PROP_ARRAY material,
                      FEM_MATERIAL_PROP_ARRAY *material_trac);

RC max_mises_stress(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                    FEM_MATERIAL_PROP_ARRAY material,
                    FEM_PHYSICAL_PROP_ARRAY physical,
                    FEM_DISP_ARRAY disp, double *max_val);
RC max_mises_stress_therm(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                          FEM_MATERIAL_PROP_ARRAY material,
                          FEM_PHYSICAL_PROP_ARRAY physical,
                          FEM_BC_ARRAY temp, FEM_DEFAULT_BC_ARRAY def_temp,
                          FEM_DISP_ARRAY disp, double *max_val);

RC mises_integration(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                     FEM_MATERIAL_PROP_ARRAY material,
                     FEM_PHYSICAL_PROP_ARRAY physical,
                     FEM_DISP_ARRAY disp, double rho, double *integ_val,
                     double *volume);

RC max_mises_adj_force(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                       FEM_MATERIAL_PROP_ARRAY material,
                       FEM_PHYSICAL_PROP_ARRAY physical, FEM_DISP_ARRAY disp,
                       double rho, FEM_BC_ARRAY *adj_force);

RC shape_sens_max_mises(FEM_ELEMENT_ARRAY surf, FEM_DISP_ARRAY disp,
                        FEM_DISP_ARRAY adj_disp,
                        FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                        FEM_MATERIAL_PROP_ARRAY material,
                        FEM_PHYSICAL_PROP_ARRAY physical,
                        double rho,
                        FEM_SENSITIVITY_ARRAY *sens);
RC shape_trac_max_mises(FEM_ELEMENT_ARRAY surf, FEM_DISP_ARRAY disp,
                        FEM_DISP_ARRAY adj_disp,
                        FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                        FEM_MATERIAL_PROP_ARRAY material,
                        FEM_PHYSICAL_PROP_ARRAY physical,
                        double rho, FEM_BC_ARRAY *sens_trac, int reduce_flag);
RC shape_trac_compliance(FEM_ELEMENT_ARRAY surf, FEM_DISP_ARRAY disp,
                         FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                         FEM_MATERIAL_PROP_ARRAY material,
                         FEM_PHYSICAL_PROP_ARRAY physical,
                         FEM_BC_ARRAY *sens_trac, int reduce_flag);
RC shape_trac_volume(FEM_ELEMENT_ARRAY surf, FEM_NODE_ARRAY node,
                     FEM_BC_ARRAY *sens_trac, int reduce_flag);

RC shape_sens_therm(FEM_ELEMENT_ARRAY surf, FEM_DISP_ARRAY disp,
                    FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                    FEM_MATERIAL_PROP_ARRAY material,
                    FEM_PHYSICAL_PROP_ARRAY physical,
                    FEM_BC_ARRAY temp, FEM_DEFAULT_BC_ARRAY def_temp,
                    FEM_SENSITIVITY_ARRAY *sens);

RC add_shape_rest(FEM_ELEMENT_ARRAY element, FEM_MATERIAL_PROP_ARRAY material,
                  FEM_PHYSICAL_PROP_ARRAY physical, FEM_BC_ARRAY *shape_rest);


#endif /* OPTIMIZATION_UTL_H */
