/*********************************************************************
 * fem_solver.h
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Yasuo ASAGA> <Takaaki NAGATANI>
 *
 * Thanks to following contributors.
 *  <Tomoyuki TSUBATA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: fem_solver.h,v 1.33 2004/02/27 06:48:03 nagatani Exp $ */

#ifndef FEM_SOLVER_H
#define FEM_SOLVER_H

#include <stdio.h>
#include "fem_struct.h"
#include "rc.h"
#include "nonzero_cg.h"

#define FEM_SOL_SS_CENTER  (1)   /* 要素の中心での値 */
#define FEM_SOL_SS_NODE    (2)   /* 要素の各節点での値 */
#define FEM_SOL_SS_NODE_AV (4)   /* 要素の各節点での値を平均 */

/* 全体剛性マトリックス */
typedef struct {
	int dim;
	long dof;
	long array_size;
	long *index1;
	long *index2;
	double *array;
} FEM_SKYLINE_MATRIX;

/* fem_solver.c */
RC fem_react_force(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                   FEM_MATERIAL_PROP_ARRAY material,
                   FEM_PHYSICAL_PROP_ARRAY physical,
                   FEM_DISP_ARRAY disp, FEM_REACT_ARRAY *react);

RC fem_solve_iccg3(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                   FEM_MATERIAL_PROP_ARRAY material,
                   FEM_PHYSICAL_PROP_ARRAY physical,
                   FEM_BC_ARRAY rest, FEM_BC_ARRAY force,
                   FEM_DEFAULT_BC_ARRAY b_force,
                   FEM_BC_ARRAY temp, FEM_DEFAULT_BC_ARRAY def_temp,
                   FEM_DISP_ARRAY *disp, FEM_STRESS_ARRAY *stress,
                   FEM_STRAIN_ARRAY *strain, int fem_sol_ss_flag);
RC fem_solve(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
             FEM_MATERIAL_PROP_ARRAY material,
             FEM_PHYSICAL_PROP_ARRAY physical,
             FEM_BC_ARRAY rest, FEM_BC_ARRAY force,
             FEM_BC_ARRAY temp, FEM_DEFAULT_BC_ARRAY def_temp,
             FEM_DISP_ARRAY *disp, FEM_STRESS_ARRAY *stress,
             FEM_STRAIN_ARRAY *strain, int fem_sol_ss_flag);
RC fem_solve_multi_force(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                         FEM_MATERIAL_PROP_ARRAY material,
                         FEM_PHYSICAL_PROP_ARRAY physical,
                         FEM_BC_ARRAY rest,
                         int num_force, const FEM_BC_ARRAY *forces,
                         FEM_BC_ARRAY temp, FEM_DEFAULT_BC_ARRAY def_temp,
                         FEM_DISP_ARRAY *disps, FEM_STRESS_ARRAY *stresses,
                         FEM_STRAIN_ARRAY *strains, int fem_sol_ss_flag);

RC fem_allocate_vector(int dim, FEM_NODE_ARRAY node, double **f_vect);
RC fem_zero_vector(int dim, FEM_NODE_ARRAY node, double *f_vect);
RC fem_free_vector(double **vect);
RC fem_force_vector2bc(const double *f_vect, FEM_ELEMENT_ARRAY element,
                       FEM_NODE_ARRAY node, FEM_BC_ARRAY *bc);
RC fem_add_nodal_force(int dim, FEM_NODE_ARRAY node,
                       FEM_BC_ARRAY force, double *f_vect);
RC fem_add_temp_force(int dim, FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                      FEM_MATERIAL_PROP_ARRAY material,
                      FEM_PHYSICAL_PROP_ARRAY physical, FEM_BC_ARRAY temp,
                      FEM_DEFAULT_BC_ARRAY def_temp, double *f_vect);
RC fem_add_body_force(int dim, FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                      FEM_MATERIAL_PROP_ARRAY material,
                      FEM_PHYSICAL_PROP_ARRAY physical,
                      FEM_DEFAULT_BC_ARRAY def_b_force, double *f_vect);
RC fem_add_init_strain_force(int dim, FEM_NODE_ARRAY node,
                             FEM_ELEMENT_ARRAY element,
                             FEM_MATERIAL_PROP_ARRAY material,
                             FEM_PHYSICAL_PROP_ARRAY physical,
                             FEM_STRAIN_ARRAY init_strain, double *f_vect);

RC fem_stress_strain(int fem_sol_ss_flag, FEM_NODE_ARRAY node,
                     FEM_ELEMENT_ARRAY element, FEM_DISP_ARRAY disp,
                     FEM_MATERIAL_PROP_ARRAY material,
                     FEM_PHYSICAL_PROP_ARRAY physical,
                     FEM_STRESS_ARRAY *stress,
                     FEM_STRAIN_ARRAY *strain);
RC fem_local_stress_strain(FEM_ELEMENT element, const double *local_point,
                           FEM_NODE_ARRAY node, FEM_DISP_ARRAY disp,
                           FEM_MATERIAL_PROP_ARRAY material,
                           FEM_PHYSICAL_PROP_ARRAY physical,
                           FEM_STRESS_STRAIN *stress,
                           FEM_STRESS_STRAIN *strain); 
RC fem_local_stress_strain_therm2(FEM_ELEMENT element,
                             const double *local_point,
                             FEM_NODE_ARRAY node, FEM_DISP_ARRAY disp,
                             FEM_MATERIAL_PROP_ARRAY material,
                             FEM_PHYSICAL_PROP_ARRAY physical,
                             FEM_BC_ARRAY temp, FEM_DEFAULT_BC_ARRAY def_temp,
                             FEM_STRESS_STRAIN *stress,
                             FEM_STRESS_STRAIN *strain);
RC fem_local_strain(FEM_ELEMENT element, const double *local_point,
                    FEM_NODE_ARRAY node, FEM_DISP_ARRAY disp,
                    FEM_STRESS_STRAIN *strain);
RC fem_local_strain_sub(FEM_ELEMENT element, const double local_point[],
                        double **node_location, const double disp_v[],
                        double **ws_mat_6_3N, FEM_STRESS_STRAIN *strain);
RC fem_local_strain2local_stress(FEM_ELEMENT element,
                                 FEM_MATERIAL_PROP_ARRAY material,
                                 FEM_PHYSICAL_PROP_ARRAY physical,
                                 FEM_STRESS_STRAIN strain,
                                 FEM_STRESS_STRAIN *stress);
TRANS_ROTATION3D d_mises_stress(FEM_STRESS_STRAIN stress);
double mises_stress(FEM_STRESS_STRAIN stress);


RC fem_make_g_matrix(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                     FEM_MATERIAL_PROP_ARRAY material,
                     FEM_PHYSICAL_PROP_ARRAY physical,
                     FEM_SKYLINE_MATRIX *g_matrix);
RC fem_allocate_g_matrix(int dim, FEM_NODE_ARRAY node,
                         FEM_ELEMENT_ARRAY element,
                         FEM_SKYLINE_MATRIX *g_matrix);
RC fem_free_g_matrix(FEM_SKYLINE_MATRIX *g_matrix);
RC fem_modify_g_matrix(FEM_NODE_ARRAY node, FEM_BC_ARRAY rest,
                       FEM_SKYLINE_MATRIX g_matrix, double *f_vect);
RC fem_modify_g_matrix_n(FEM_NODE_ARRAY node, FEM_BC_ARRAY rest,
                         FEM_SKYLINE_MATRIX g_matrix, double **f_vect,
                         int f_vec_size);
RC fem_decomp_g_matrix(FEM_SKYLINE_MATRIX g_matrix);
RC fem_subst_g_matrix(FEM_NODE_ARRAY node, FEM_SKYLINE_MATRIX g_matrix,
                      double *f_vect, FEM_DISP_ARRAY *disp);

RC fem_impose_elem_matrix(FEM_SKYLINE_MATRIX g_matrix,
                          FEM_ELEM_MATRIX elem_matrix, double fact);

RC fem_make_elem_matrix(FEM_ELEMENT element, FEM_NODE_ARRAY node,
                        FEM_MATERIAL_PROP_ARRAY material,
                        FEM_PHYSICAL_PROP_ARRAY physical,
                        FEM_ELEM_MATRIX *elem_matrix);


RC fem_total_dof(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                 int *total_dof);

RC fem_allocate_nonzero3(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                         NONZERO_MATRIX3 *matrix);

RC fem_impose_nonzero3(NONZERO_MATRIX3 matrix,
                       FEM_ELEM_MATRIX elem_matrix, double fact);

RC fem_modify_nonzero3(FEM_NODE_ARRAY node, FEM_BC_ARRAY rest,
                       NONZERO_MATRIX3 matrix, double *f_vect);

RC fem_modify_nonzero3_n(FEM_NODE_ARRAY node, FEM_BC_ARRAY rest,
                         NONZERO_MATRIX3 matrix, double **f_vect,
                         int vect_size);

RC fem_rest2vect(int dim, FEM_NODE_ARRAY node, FEM_BC_ARRAY rest,
                 double **d_vect);
RC fem_disp2vect(int dim, FEM_NODE_ARRAY node, FEM_DISP_ARRAY disp,
                 double **d_vect);
RC fem_vect2disp(int dim, FEM_NODE_ARRAY node, const double *d_vect,
                 FEM_DISP_ARRAY *disp);
RC fem_vect2react(int dim, FEM_NODE_ARRAY node, const double *r_vect,
                  FEM_REACT_ARRAY *react);

RC fem_iccg_nonzero3(NONZERO_MATRIX3 matrix, const double *f_vect,
                     double **d_vect, const double *init_vect);
RC fem_scale_iccg_nonzero3(NONZERO_MATRIX3 matrix, double *f_vect,
                           double **d_vect, const double *init_vect);

RC fem_make_mass_matrix(FEM_ELEMENT element, FEM_NODE_ARRAY node,
                        FEM_MATERIAL_PROP_ARRAY material,
                        FEM_PHYSICAL_PROP_ARRAY physical,
                        FEM_ELEM_MATRIX *elem_matrix);

RC vect2stress_strain(int vect_size, const double *vect,
                      FEM_STRESS_STRAIN *ss);
RC stress_strain2vect(int vect_size, FEM_STRESS_STRAIN ss, double *vect);
RC fill_disp_vector(FEM_ELEMENT element, FEM_DISP_ARRAY disp, double *disp_v);

RC make_B_matrix(FEM_ELEMENT element, const double *g_point,
                 double **node_location, double **B_matrix, double *det_J);
RC get_spring_const(FEM_PHYSICAL_PROP_ARRAY physical, int p_label,
                    double *spring_const);
RC get_thickness(FEM_ELEMENT element,
                 FEM_PHYSICAL_PROP_ARRAY physical, double *thickness);
RC get_D_matrix(FEM_ELEMENT element, FEM_MATERIAL_PROP_ARRAY material,
                FEM_PHYSICAL_PROP_ARRAY physical, double **D_matrix);
RC get_D_matrix_therm(FEM_ELEMENT element, FEM_MATERIAL_PROP_ARRAY material,
                      FEM_PHYSICAL_PROP_ARRAY physical, double **D_matrix);

RC fem_make_elem_matrix_therm(FEM_ELEMENT element, FEM_NODE_ARRAY node,
                              FEM_MATERIAL_PROP_ARRAY material,
                              FEM_PHYSICAL_PROP_ARRAY physical,
                              FEM_ELEM_MATRIX *elem_matrix);
RC fem_modify_g_matrix_therm(FEM_NODE_ARRAY node, FEM_BC_ARRAY temp,
                             FEM_SKYLINE_MATRIX g_matrix, double *s_vect);
RC fem_subst_g_matrix_therm(FEM_NODE_ARRAY node, FEM_SKYLINE_MATRIX g_matrix,
                            double *s_vect, FEM_BC_ARRAY *temp);
RC fem_solve_therm(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                   FEM_MATERIAL_PROP_ARRAY material,
                   FEM_PHYSICAL_PROP_ARRAY physical,
                   FEM_BC_ARRAY bc_temp, FEM_BC_ARRAY *result_temp);


/* fem_solver_6dof */
RC fem_solve_iccg6(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                   FEM_RIGID_ELEMENT_ARRAY rigid,
                   FEM_MATERIAL_PROP_ARRAY material,
                   FEM_PHYSICAL_PROP_ARRAY physical,
                   FEM_LOCAL_COORD_ARRAY coord,
                   FEM_BC_ARRAY rest, FEM_BC_ARRAY force,
                   FEM_BC_ARRAY temp, FEM_DEFAULT_BC_ARRAY def_temp,
                   FEM_DEFAULT_BC_ARRAY def_body_force,
                   FEM_DISP_ARRAY *disp, FEM_STRESS_ARRAY *stress,
                   FEM_STRAIN_ARRAY *strain, int fem_sol_ss_flag);
RC fem_allocate_nonzero6(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                         FEM_RIGID_ELEMENT_ARRAY rigid_elem,
                         NONZERO_MATRIX3 *matrix);
RC fem_impose_nonzero6(NONZERO_MATRIX3 matrix,
                       FEM_ELEM_MATRIX elem_matrix, double fact);
RC fem_modify_node_matrix6(FEM_NODE_ARRAY node, FEM_NODE_MATRIX_ARRAY nmat,
                           double fact, NONZERO_MATRIX3 matrix);
RC fem_modify_spring6(FEM_NODE_ARRAY node, FEM_RIGID_ELEMENT_ARRAY src_rigid,
                      FEM_PHYSICAL_PROP_ARRAY physical,
                      NONZERO_MATRIX3 matrix);
RC fem_modify_mass6(FEM_NODE_ARRAY node, FEM_RIGID_ELEMENT_ARRAY rigid,
                    double fact, NONZERO_MATRIX3 matrix);
RC fem_modify_nonzero6(FEM_NODE_ARRAY node, FEM_LOCAL_COORD_ARRAY coord,
                       FEM_BC_ARRAY rest, FEM_RIGID_ELEMENT_ARRAY src_rigid,
                       NONZERO_MATRIX3 matrix,
                       double **f_vect, int vect_size, int active_dof_flag[]);
RC fem_rigid_recover6(FEM_NODE_ARRAY node,
                      FEM_RIGID_ELEMENT_ARRAY src_rigid, double d_vect[]);
RC fem_add_nodal_force6(FEM_NODE_ARRAY node, FEM_LOCAL_COORD_ARRAY coord,
                        FEM_BC_ARRAY force, double *f_vect);
RC fem_add_pressure_force6(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                           FEM_LOCAL_COORD_ARRAY coord,
                           FEM_ELEM_BC_ARRAY press, double *f_vect);
RC fem_force_local_pre6(FEM_NODE_ARRAY node, FEM_LOCAL_COORD_ARRAY coord,
                        double f_vect[]);
RC fem_disp_local_recover6(FEM_NODE_ARRAY node, FEM_LOCAL_COORD_ARRAY coord,
                           double d_vect[]);
RC fem_nonzero2skyline6(NONZERO_MATRIX3 nonzero, FEM_SKYLINE_MATRIX *skyline,
                        long allocatable_size);
RC fem_make_beam_matrix(FEM_ELEMENT element, FEM_NODE_ARRAY node,
                        FEM_MATERIAL_PROP_ARRAY material,
                        FEM_PHYSICAL_PROP_ARRAY physical,
                        FEM_ELEM_MATRIX *elem_matrix);
RC fem_make_beam_mass_matrix(FEM_ELEMENT element, FEM_NODE_ARRAY node,
                             FEM_MATERIAL_PROP_ARRAY material,
                             FEM_PHYSICAL_PROP_ARRAY physical,
                             FEM_ELEM_MATRIX *elem_matrix);
RC fem_solve_shell_iccg6(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                         FEM_RIGID_ELEMENT_ARRAY rigid,
                         FEM_MATERIAL_PROP_ARRAY material,
                         FEM_PHYSICAL_PROP_ARRAY physical,
                         FEM_LOCAL_COORD_ARRAY coord,
                         FEM_BC_ARRAY rest, FEM_BC_ARRAY force,
                         FEM_BC_ARRAY temp, FEM_DEFAULT_BC_ARRAY def_temp,
                         FEM_DEFAULT_BC_ARRAY def_body_force,
                         FEM_ELEM_BC_ARRAY press, double distance[],
                         FEM_DISP_ARRAY *disp, FEM_STRESS_ARRAY stress[],
                         FEM_STRAIN_ARRAY strain[], FEM_STRESS_ARRAY *iforce,
                         FEM_STRESS_ARRAY *moment,
                         FEM_STRESS_ARRAY *shear_force,
                         FEM_STRAIN_ARRAY *curv, int fem_sol_ss_flag);
RC fem_make_shell_matrix(FEM_ELEMENT element, FEM_NODE_ARRAY node,
                         FEM_MATERIAL_PROP_ARRAY material,
                         FEM_PHYSICAL_PROP_ARRAY physical,
                         FEM_ELEM_MATRIX *elem_matrix);
RC fem_stress_strain_shell(int fem_sol_ss_flag, double *distance,
                           FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                           FEM_DISP_ARRAY disp,
                           FEM_MATERIAL_PROP_ARRAY material,
                           FEM_PHYSICAL_PROP_ARRAY physical,
                           FEM_STRESS_ARRAY stress[],
                           FEM_STRAIN_ARRAY strain[]);
RC fem_resultant_stress_strain_shell(int fem_sol_ss_flag, FEM_NODE_ARRAY node,
                                     FEM_ELEMENT_ARRAY element,
                                     FEM_DISP_ARRAY disp,
                                     FEM_MATERIAL_PROP_ARRAY material,
                                     FEM_PHYSICAL_PROP_ARRAY physical,
                                     FEM_STRESS_ARRAY *force,
                                     FEM_STRESS_ARRAY *moment,
                                     FEM_STRESS_ARRAY *shear_force,
                                     FEM_STRAIN_ARRAY *curv);
RC fem_local_stress_strain_shell(double distance[], FEM_ELEMENT element,
                                 const double *local_point,
                                 FEM_NODE_ARRAY node, FEM_DISP_ARRAY disp,
                                 FEM_MATERIAL_PROP_ARRAY material,
                                 FEM_PHYSICAL_PROP_ARRAY physical,
                                 FEM_WORK_SET ws, FEM_STRESS_STRAIN stress[],
                                 FEM_STRESS_STRAIN strain[]);
RC fem_local_resultant_stress_strain(FEM_ELEMENT element,
                                     const double *local_point,
                                     FEM_NODE_ARRAY node, FEM_DISP_ARRAY disp,
                                     FEM_MATERIAL_PROP_ARRAY material,
                                     FEM_PHYSICAL_PROP_ARRAY physical,
                                     FEM_WORK_SET ws, FEM_STRESS_STRAIN *force,
                                     FEM_STRESS_STRAIN *moment,
                                     FEM_STRESS_STRAIN *shear_force,
                                     FEM_STRESS_STRAIN *curv);

/* nl_solver.c */
RC fem_NLsolve_iccg3(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                     FEM_MATERIAL_PROP_ARRAY material,
                     FEM_PHYSICAL_PROP_ARRAY physical,
                     FEM_BC_ARRAY rest, int inc_num,
                     FEM_BC_ARRAY force, FEM_DEFAULT_BC_ARRAY b_force,
                     FEM_DISP_ARRAY *disp, FEM_STRESS_ARRAY *stress,
                     FEM_STRAIN_ARRAY *strain, int fem_sol_ss_flag);
RC fem_add_NLreact_force(FEM_ELEMENT_ARRAY element, FEM_NODE_ARRAY node,
                         FEM_MATERIAL_PROP_ARRAY material,
                         FEM_PHYSICAL_PROP_ARRAY physical, FEM_DISP_ARRAY disp,
                         double *f_vect);
RC fem_stress_Green_strain(int fem_sol_ss_flag, FEM_NODE_ARRAY node,
                           FEM_ELEMENT_ARRAY element, FEM_DISP_ARRAY disp,
                           FEM_MATERIAL_PROP_ARRAY material,
                           FEM_PHYSICAL_PROP_ARRAY physical,
                           FEM_STRESS_ARRAY *stress, FEM_STRAIN_ARRAY *strain);
RC fem_local_Green_strain(FEM_ELEMENT element, const double *local_point,
                          FEM_NODE_ARRAY node, FEM_DISP_ARRAY disp,
                          FEM_STRESS_STRAIN *strain);
RC make_AG_matrix(FEM_ELEMENT element, const double *g_point,
                  double **node_location, FEM_DISP_ARRAY disp,
                  double **A_matrix, double **G_matrix);
RC fem_make_tangent_matrix(FEM_ELEMENT element, FEM_NODE_ARRAY node,
                           FEM_MATERIAL_PROP_ARRAY material,
                           FEM_PHYSICAL_PROP_ARRAY physical,
                           FEM_DISP_ARRAY disp, FEM_ELEM_MATRIX *elem_matrix);
RC set_M_matrix(int dim, double *stress, double **M_matrix);


/* eigen_solver.c */
RC fem_eigen_solve_iccg3(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                         FEM_MATERIAL_PROP_ARRAY material,
                         FEM_PHYSICAL_PROP_ARRAY physical,
                         double *evals, FEM_DISP_ARRAY *disps,
                         double sigma, int max_e_num,
                         int *e_num, int max_step);

/* fem_solver_shell.c */
RC get_D_matrix_shell(FEM_ELEMENT element, FEM_MATERIAL_PROP_ARRAY material,
                      FEM_PHYSICAL_PROP_ARRAY physical, double **D_matrix);
RC fem_impose_shell_elem_matrix(int plane_dim, FEM_ELEMENT element,
                                FEM_MATERIAL_PROP_ARRAY material,
                                FEM_PHYSICAL_PROP_ARRAY physical,
                                FEM_ELEM_MATRIX elem_matrix,
                                FEM_ELEM_MATRIX *total_elem_matrix);
RC fem_make_shell_elem_matrix(FEM_ELEMENT element, FEM_NODE_ARRAY node,
                              FEM_MATERIAL_PROP_ARRAY material,
                              FEM_PHYSICAL_PROP_ARRAY physical,
                              FEM_ELEM_MATRIX *elem_matrix,
                              FEM_WORK_SET ws);
RC make_B_matrix_shell(FEM_ELEMENT element, FEM_MATERIAL_PROP_ARRAY material,
                       FEM_PHYSICAL_PROP_ARRAY physical, FEM_WORK_SET ws,
                       const double *g_point, double **node_location,
                       double **B_matrix, double *det_J);
RC fem_make_shell_matrix(FEM_ELEMENT element, FEM_NODE_ARRAY node,
                         FEM_MATERIAL_PROP_ARRAY material,
                         FEM_PHYSICAL_PROP_ARRAY physical,
                         FEM_ELEM_MATRIX *elem_matrix);
RC fem_stress_strain_shell(int fem_sol_ss_flag, double *distance,
                           FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                           FEM_DISP_ARRAY disp,
                           FEM_MATERIAL_PROP_ARRAY material,
                           FEM_PHYSICAL_PROP_ARRAY physical,
                           FEM_STRESS_ARRAY stress[],
                           FEM_STRAIN_ARRAY strain[]);
RC fem_resultant_stress_strain_shell(int fem_sol_ss_flag, FEM_NODE_ARRAY node,
                                     FEM_ELEMENT_ARRAY element,
                                     FEM_DISP_ARRAY disp,
                                     FEM_MATERIAL_PROP_ARRAY material,
                                     FEM_PHYSICAL_PROP_ARRAY physical,
                                     FEM_STRESS_ARRAY *force,
                                     FEM_STRESS_ARRAY *moment,
                                     FEM_STRESS_ARRAY *shear_force,
                                     FEM_STRAIN_ARRAY *curv);
RC fem_local_stress_strain_shell(double distance[], FEM_ELEMENT element,
                                 const double *local_point,
                                 FEM_NODE_ARRAY node, FEM_DISP_ARRAY disp,
                                 FEM_MATERIAL_PROP_ARRAY material,
                                 FEM_PHYSICAL_PROP_ARRAY physical,
                                 FEM_WORK_SET ws, FEM_STRESS_STRAIN stress[],
                                 FEM_STRESS_STRAIN strain[]);
RC fem_local_resultant_stress_strain(FEM_ELEMENT element,
                                     const double *local_point,
                                     FEM_NODE_ARRAY node, FEM_DISP_ARRAY disp,
                                     FEM_MATERIAL_PROP_ARRAY material,
                                     FEM_PHYSICAL_PROP_ARRAY physical,
                                     FEM_WORK_SET ws, FEM_STRESS_STRAIN *force,
                                     FEM_STRESS_STRAIN *moment,
                                     FEM_STRESS_STRAIN *shear_force,
                                     FEM_STRESS_STRAIN *curv);

/* flow_solver.c */
RC steady_flow_analysis(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY elem,
                        FEM_BC_ARRAY bc_velo,FEM_BC_ARRAY bc_press, 
                        int velo_order, int press_order,
                        double density, double viscosity,
                        FEM_VELO_ARRAY *flow);

RC steady_adjflow_analysis(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY elem,
                           FEM_BC_ARRAY bc_velo, FEM_BC_ARRAY bc_press,
                           int velo_order, int press_order,
                           double density, double viscosity,
                           FEM_VELO_ARRAY flow, FEM_VELO_ARRAY *adjoint_flow);

RC flow_solve_for_optimization(double delta_t, int iteration, double Re,
                               FEM_NODE_ARRAY node,
                               FEM_ELEMENT_ARRAY element,
                               FEM_BC_ARRAY bc_velo, FEM_BC_ARRAY bc_press,
                               FEM_VELO_ARRAY *vp);
RC flow_solve_non_steady_flow(int plot_step, double delta_t, int iteration,
                              double Re, FEM_NODE_ARRAY node,
                              FEM_ELEMENT_ARRAY element, FEM_BC_ARRAY bc_velo,
                              FEM_BC_ARRAY bc_press, FEM_VELO_ARRAY *vp);

RC flow_adjoint_analysis(double Re, FEM_NODE_ARRAY node,
                         FEM_ELEMENT_ARRAY element, FEM_BC_ARRAY bc,
                         FEM_VELO_ARRAY vp, FEM_VELO_ARRAY *adjoint_vp);
RC flow_adjoint_analysis2(double Re, FEM_NODE_ARRAY node,
                          FEM_ELEMENT_ARRAY element, FEM_BC_ARRAY bc,
                          FEM_VELO_ARRAY vp, FEM_VELO_ARRAY *adjoint_vp);

RC strain_velocity(int fem_sol_ss_flag, FEM_NODE_ARRAY node,
                   FEM_ELEMENT_ARRAY element, FEM_VELO_ARRAY velo,
                   FEM_STRAIN_ARRAY *velo_strain);
RC pick_out_wall_bc_node(FEM_BC_ARRAY bc, FEM_NODE_ARRAY *wall_node);
RC change_bc_type4flow(FEM_BC_ARRAY *bc);


#endif /* FEM_SOLVER_H */
