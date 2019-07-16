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

/* $Id: fem_solver.h 476 2005-09-22 07:16:38Z aoyama $ */

#ifndef FEM_SOLVER_H
#define FEM_SOLVER_H

#include <stdio.h>
#include "fem_struct.h"
#include "rc.h"
#include "nonzero_cg.h"

#define SOL_SS_CENTER  (1)   /* 要素の中心での値 */
#define SOL_SS_NODE    (2)   /* 要素の各節点での値 */
#define SOL_SS_NODE_AV (4)   /* 要素の各節点での値を平均 */

/* 全体剛性マトリックス */
typedef struct {
	int dim;
	long dof;
	long array_size;
	long *index1;
	long *index2;
	double *array;
} SKYLINE_MATRIX;

/* fem_solver.c */
RC react_force(NODE_ARRAY node, ELEMENT_ARRAY element,
               MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
               DISP_ARRAY disp, REACT_ARRAY *react);

RC fem_solve_iccg3(NODE_ARRAY node, ELEMENT_ARRAY element,
                   MATERIAL_PROP_ARRAY material,
                   PHYSICAL_PROP_ARRAY physical,
                   BC_ARRAY rest, BC_ARRAY force,
                   DEFAULT_BC_ARRAY b_force,
                   BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                   DISP_ARRAY *disp, STRESS_ARRAY *stress,
                   STRAIN_ARRAY *strain, int sol_ss_flag);
RC fem_solve(NODE_ARRAY node, ELEMENT_ARRAY element,
             MATERIAL_PROP_ARRAY material,
             PHYSICAL_PROP_ARRAY physical,
             BC_ARRAY rest, BC_ARRAY force,
             BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
             DISP_ARRAY *disp, STRESS_ARRAY *stress,
             STRAIN_ARRAY *strain, int sol_ss_flag);
RC fem_solve_multi_force(NODE_ARRAY node, ELEMENT_ARRAY element,
                         MATERIAL_PROP_ARRAY material,
                         PHYSICAL_PROP_ARRAY physical,
                         BC_ARRAY rest, int num_force, const BC_ARRAY *forces,
                         BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                         DISP_ARRAY *disps, STRESS_ARRAY *stresses,
                         STRAIN_ARRAY *strains, int sol_ss_flag);

RC fem_allocate_vector(int dim, NODE_ARRAY node, double **f_vect);
RC fem_zero_vector(int dim, NODE_ARRAY node, double *f_vect);
RC fem_free_vector(double **vect);
RC fem_force_vector2bc(const double *f_vect, ELEMENT_ARRAY element,
                       NODE_ARRAY node, BC_ARRAY *bc);
RC add_nodal_force(int dim, NODE_ARRAY node, BC_ARRAY force, double *f_vect);
RC add_temp_force(int dim, NODE_ARRAY node, ELEMENT_ARRAY element,
                  MATERIAL_PROP_ARRAY material,
                  PHYSICAL_PROP_ARRAY physical, BC_ARRAY temp,
                  DEFAULT_BC_ARRAY def_temp, double *f_vect);
RC add_body_force(int dim, NODE_ARRAY node, ELEMENT_ARRAY element,
                  MATERIAL_PROP_ARRAY material,
                  PHYSICAL_PROP_ARRAY physical,
                  DEFAULT_BC_ARRAY def_b_force, double *f_vect);
RC add_init_strain_force(int dim, NODE_ARRAY node, ELEMENT_ARRAY element,
                         MATERIAL_PROP_ARRAY material,
                         PHYSICAL_PROP_ARRAY physical,
                         STRAIN_ARRAY init_strain, double *f_vect);

RC stress_strain(int sol_ss_flag, NODE_ARRAY node,
                 ELEMENT_ARRAY element, DISP_ARRAY disp,
                 MATERIAL_PROP_ARRAY material,
                 PHYSICAL_PROP_ARRAY physical,
                 STRESS_ARRAY *stress, STRAIN_ARRAY *strain);
RC local_points_stress_strain(ELEMENT element, const double *local_point,
                              NODE_ARRAY node, DISP_ARRAY disp,
                              MATERIAL_PROP_ARRAY material,
                              PHYSICAL_PROP_ARRAY physical,
                              STRESS_STRAIN *stress, STRESS_STRAIN *strain); 
RC local_points_stress_strain_therm(ELEMENT element, const double *local_point,
                                    NODE_ARRAY node, DISP_ARRAY disp,
                                    MATERIAL_PROP_ARRAY material,
                                    PHYSICAL_PROP_ARRAY physical,
                                    BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                                    double fact,
                                    STRESS_STRAIN *stress,
                                    STRESS_STRAIN *strain);
RC local_points_strain(ELEMENT element, const double *local_point,
                       NODE_ARRAY node, DISP_ARRAY disp, STRESS_STRAIN *strain);
RC local_points_strain_sub(ELEMENT element, const double local_point[],
                           double **node_location, const double disp_v[],
                           double **ws_mat_6_3N, STRESS_STRAIN *strain);
RC local_points_strain_therm(ELEMENT element, const double *local_point,
                             NODE_ARRAY node, DISP_ARRAY disp,
                             MATERIAL_PROP_ARRAY material,
                             PHYSICAL_PROP_ARRAY physical,
                             BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                             double fact,
                             STRESS_STRAIN *strain);
RC local_strain2local_stress(ELEMENT element,
                             MATERIAL_PROP_ARRAY material,
                             PHYSICAL_PROP_ARRAY physical,
                             STRESS_STRAIN strain, STRESS_STRAIN *stress);
VECT3DR d_mises_stress(STRESS_STRAIN stress);
double mises_stress(STRESS_STRAIN stress);


RC make_g_matrix(NODE_ARRAY node, ELEMENT_ARRAY element,
                 MATERIAL_PROP_ARRAY material,
                 PHYSICAL_PROP_ARRAY physical,
                 SKYLINE_MATRIX *g_matrix);
RC allocate_g_matrix(int dim, NODE_ARRAY node, ELEMENT_ARRAY element,
                     SKYLINE_MATRIX *g_matrix);
RC free_g_matrix(SKYLINE_MATRIX *g_matrix);
RC modify_g_matrix(NODE_ARRAY node, BC_ARRAY rest,
                   SKYLINE_MATRIX g_matrix, double *f_vect);
RC modify_g_matrix_n(NODE_ARRAY node, BC_ARRAY rest,
                     SKYLINE_MATRIX g_matrix, double **f_vect, int f_vec_size);
RC decomp_g_matrix(SKYLINE_MATRIX g_matrix);
RC subst_g_matrix(NODE_ARRAY node, SKYLINE_MATRIX g_matrix,
                  double *f_vect, DISP_ARRAY *disp);

RC impose_elem_matrix(SKYLINE_MATRIX g_matrix,
                      ELEM_MATRIX elem_matrix, double fact);
RC make_elem_matrix(ELEMENT element, NODE_ARRAY node,
                    MATERIAL_PROP_ARRAY material,
                    PHYSICAL_PROP_ARRAY physical,
                    ELEM_MATRIX *elem_matrix);

RC fem_total_dof(NODE_ARRAY node, ELEMENT_ARRAY element, int *total_dof);

RC allocate_nonzero3(NODE_ARRAY node, ELEMENT_ARRAY element,
                     NONZERO_MATRIX3 *matrix);
RC impose_nonzero3(NONZERO_MATRIX3 matrix,
                   ELEM_MATRIX elem_matrix, double fact);
RC modify_nonzero3(NODE_ARRAY node, BC_ARRAY rest,
                   NONZERO_MATRIX3 matrix, double *f_vect);
RC modify_nonzero3_n(NODE_ARRAY node, BC_ARRAY rest,
                     NONZERO_MATRIX3 matrix, double **f_vect, int vect_size);

RC rest2vect(int dim, NODE_ARRAY node, BC_ARRAY rest, double **d_vect);
RC disp2vect(int dim, NODE_ARRAY node, DISP_ARRAY disp, double **d_vect);
RC vect2disp(int dim, NODE_ARRAY node, const double *d_vect, DISP_ARRAY *disp);
RC vect2react(int dim, NODE_ARRAY node, const double *r_vect,
              REACT_ARRAY *react);

RC iccg_nonzero3(NONZERO_MATRIX3 matrix, const double *f_vect,
                 double **d_vect, const double *init_vect);
RC scale_iccg_nonzero3(NONZERO_MATRIX3 matrix, double *f_vect,
                       double **d_vect, const double *init_vect);

RC make_mass_matrix(ELEMENT element, NODE_ARRAY node,
                    MATERIAL_PROP_ARRAY material,
                    PHYSICAL_PROP_ARRAY physical,
                    ELEM_MATRIX *elem_matrix);

RC vect2stress_strain(int vect_size, const double *vect,
                      STRESS_STRAIN *ss);
RC stress_strain2vect(int vect_size, STRESS_STRAIN ss, double *vect);
RC fill_disp_vector(ELEMENT element, DISP_ARRAY disp, double *disp_v);

RC make_B_matrix(ELEMENT element, const double *g_point,
                 double **node_location, double **B_matrix, double *det_J);
RC get_spring_const(PHYSICAL_PROP_ARRAY physical, int p_label,
                    double *spring_const);
RC get_thickness(ELEMENT element,
                 PHYSICAL_PROP_ARRAY physical, double *thickness);
RC get_D_matrix(ELEMENT element, MATERIAL_PROP_ARRAY material,
                PHYSICAL_PROP_ARRAY physical, double **D_matrix);
RC get_D_matrix_therm(ELEMENT element, MATERIAL_PROP_ARRAY material,
                      PHYSICAL_PROP_ARRAY physical, double **D_matrix);

RC make_elem_matrix_therm(ELEMENT element, NODE_ARRAY node,
                          MATERIAL_PROP_ARRAY material,
                          PHYSICAL_PROP_ARRAY physical,
                          ELEM_MATRIX *elem_matrix);
RC modify_g_matrix_therm(NODE_ARRAY node, BC_ARRAY temp,
                         SKYLINE_MATRIX g_matrix, double *s_vect);
RC subst_g_matrix_therm(NODE_ARRAY node, SKYLINE_MATRIX g_matrix,
                        double *s_vect, BC_ARRAY *temp);
RC fem_solve_therm(NODE_ARRAY node, ELEMENT_ARRAY element,
                   MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
                   BC_ARRAY bc_temp, BC_ARRAY *result_temp);
double node_temp(BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp, int node_label);

/* fem_solver_6dof */
RC fem_solve_iccg6(NODE_ARRAY node, ELEMENT_ARRAY element,
                   RIGID_ELEMENT_ARRAY rigid,
                   MATERIAL_PROP_ARRAY material,
                   PHYSICAL_PROP_ARRAY physical,
                   LOCAL_COORD_ARRAY coord,
                   BC_ARRAY rest, BC_ARRAY force,
                   BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                   DEFAULT_BC_ARRAY def_body_force,
                   DISP_ARRAY *disp, STRESS_ARRAY *stress,
                   STRAIN_ARRAY *strain, int sol_ss_flag);
RC allocate_nonzero6(NODE_ARRAY node, ELEMENT_ARRAY element,
                     RIGID_ELEMENT_ARRAY rigid_elem, NONZERO_MATRIX3 *matrix);
RC impose_nonzero6(NONZERO_MATRIX3 matrix, ELEM_MATRIX elem_matrix,double fact);
RC modify_node_matrix6(NODE_ARRAY node, NODE_MATRIX_ARRAY nmat,
                       double fact, NONZERO_MATRIX3 matrix);
RC modify_spring6(NODE_ARRAY node, RIGID_ELEMENT_ARRAY src_rigid,
                  PHYSICAL_PROP_ARRAY physical, NONZERO_MATRIX3 matrix);
RC modify_mass6(NODE_ARRAY node, RIGID_ELEMENT_ARRAY rigid,
                double fact, NONZERO_MATRIX3 matrix);
RC modify_nonzero6(NODE_ARRAY node, LOCAL_COORD_ARRAY coord,
                   BC_ARRAY rest, RIGID_ELEMENT_ARRAY src_rigid,
                   NONZERO_MATRIX3 matrix,
                   double **f_vect, int vect_size, int active_dof_flag[]);
RC rigid_recover6(NODE_ARRAY node,
                  RIGID_ELEMENT_ARRAY src_rigid, double d_vect[]);
RC add_nodal_force6(NODE_ARRAY node, LOCAL_COORD_ARRAY coord,
                    BC_ARRAY force, double *f_vect);
RC add_pressure_force6(NODE_ARRAY node, ELEMENT_ARRAY element,
                       LOCAL_COORD_ARRAY coord,
                       ELEM_BC_ARRAY press, double *f_vect);
RC force_local_pre6(NODE_ARRAY node, LOCAL_COORD_ARRAY coord,
                   double f_vect[]);
RC disp_local_recover6(NODE_ARRAY node, LOCAL_COORD_ARRAY coord,
                      double d_vect[]);
RC nonzero2skyline6(NONZERO_MATRIX3 nonzero, SKYLINE_MATRIX *skyline,
                    long allocatable_size);
RC make_beam_matrix(ELEMENT element, NODE_ARRAY node,
                    MATERIAL_PROP_ARRAY material,
                    PHYSICAL_PROP_ARRAY physical,
                    ELEM_MATRIX *elem_matrix);
RC make_beam_mass_matrix(ELEMENT element, NODE_ARRAY node,
                         MATERIAL_PROP_ARRAY material,
                         PHYSICAL_PROP_ARRAY physical,
                         ELEM_MATRIX *elem_matrix);
RC fem_solve_shell_iccg6(NODE_ARRAY node, ELEMENT_ARRAY element,
                         RIGID_ELEMENT_ARRAY rigid,
                         MATERIAL_PROP_ARRAY material,
                         PHYSICAL_PROP_ARRAY physical,
                         LOCAL_COORD_ARRAY coord,
                         BC_ARRAY rest, BC_ARRAY force,
                         BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                         DEFAULT_BC_ARRAY def_body_force,
                         ELEM_BC_ARRAY press, double distance[],
                         DISP_ARRAY *disp, STRESS_ARRAY stress[],
                         STRAIN_ARRAY strain[], STRESS_ARRAY *iforce,
                         STRESS_ARRAY *moment,
                         STRESS_ARRAY *shear_force,
                         STRAIN_ARRAY *curv, int sol_ss_flag);
RC make_shell_matrix(ELEMENT element, NODE_ARRAY node,
                     MATERIAL_PROP_ARRAY material,
                     PHYSICAL_PROP_ARRAY physical,
                     ELEM_MATRIX *elem_matrix);
RC stress_strain_shell(int sol_ss_flag, double *distance,
                       NODE_ARRAY node, ELEMENT_ARRAY element,
                       DISP_ARRAY disp,
                       MATERIAL_PROP_ARRAY material,
                       PHYSICAL_PROP_ARRAY physical,
                       STRESS_ARRAY stress[], STRAIN_ARRAY strain[]);
RC resultant_stress_strain_shell(int sol_ss_flag, NODE_ARRAY node,
                                 ELEMENT_ARRAY element, DISP_ARRAY disp,
                                 MATERIAL_PROP_ARRAY material,
                                 PHYSICAL_PROP_ARRAY physical,
                                 STRESS_ARRAY *force, STRESS_ARRAY *moment,
                                 STRESS_ARRAY *shear_force, STRAIN_ARRAY *curv);
RC local_stress_strain_shell(double distance[], ELEMENT element,
                             const double *local_point,
                             NODE_ARRAY node, DISP_ARRAY disp,
                             MATERIAL_PROP_ARRAY material,
                             PHYSICAL_PROP_ARRAY physical,
                             WORK_SET ws, STRESS_STRAIN stress[],
                             STRESS_STRAIN strain[]);
RC local_resultant_stress_strain(ELEMENT element,
                                 const double *local_point,
                                 NODE_ARRAY node, DISP_ARRAY disp,
                                 MATERIAL_PROP_ARRAY material,
                                 PHYSICAL_PROP_ARRAY physical,
                                 WORK_SET ws, STRESS_STRAIN *force,
                                 STRESS_STRAIN *moment,
                                 STRESS_STRAIN *shear_force,
                                 STRESS_STRAIN *curv);

/* nl_solver.c */
RC fem_cylindrical_AL_solve(NODE_ARRAY node, ELEMENT_ARRAY element,
                            MATERIAL_PROP_ARRAY material,
                            PHYSICAL_PROP_ARRAY physical,
                            BC_ARRAY rest, BC_ARRAY force,
                            DEFAULT_BC_ARRAY b_force,
                            BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                            double load_factor,
                            DISP_ARRAY *disp, STRESS_ARRAY *stress,
                            STRAIN_ARRAY *strain, int sol_ss_flag);
RC fem_NL_solve_iccg3(NODE_ARRAY node, ELEMENT_ARRAY element,
                      MATERIAL_PROP_ARRAY material,
                      PHYSICAL_PROP_ARRAY physical,
                      BC_ARRAY rest, BC_ARRAY force,
                      DEFAULT_BC_ARRAY b_force,
                      BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                      DISP_ARRAY *disp, STRESS_ARRAY *stress,
                      STRAIN_ARRAY *strain, int sol_ss_flag);
RC make_NL_g_matrix(NODE_ARRAY node, ELEMENT_ARRAY element,
                    MATERIAL_PROP_ARRAY material,
                    PHYSICAL_PROP_ARRAY physical,
                    DISP_ARRAY disp,
                    SKYLINE_MATRIX *K);
RC make_NL_nonzero_matrix3(NODE_ARRAY node, ELEMENT_ARRAY element,
                           MATERIAL_PROP_ARRAY material,
                           PHYSICAL_PROP_ARRAY physical,
                           DISP_ARRAY disp,
                           NONZERO_MATRIX3 *K);
RC make_NL_react_force(ELEMENT_ARRAY element, NODE_ARRAY node,
                       MATERIAL_PROP_ARRAY material,
                       PHYSICAL_PROP_ARRAY physical,
                       DISP_ARRAY disp, double f_vect[]);
RC make_tangent_elem_matrix(ELEMENT element, NODE_ARRAY node,
                            MATERIAL_PROP_ARRAY material,
                            PHYSICAL_PROP_ARRAY physical,
                            WORK_SET ws, DISP_ARRAY disp,
                            ELEM_MATRIX *elem_matrix);
RC Cauchy_stress_Almansi_strain(int sol_ss_flag, NODE_ARRAY node,
                                ELEMENT_ARRAY element, DISP_ARRAY disp,
                                MATERIAL_PROP_ARRAY material,
                                PHYSICAL_PROP_ARRAY physical,
                                STRESS_ARRAY *stress,
                                STRAIN_ARRAY *strain);
RC PK2_stress_Green_strain(int sol_ss_flag, NODE_ARRAY node,
                           ELEMENT_ARRAY element, DISP_ARRAY disp,
                           MATERIAL_PROP_ARRAY material,
                           PHYSICAL_PROP_ARRAY physical,
                           STRESS_ARRAY *stress, STRAIN_ARRAY *strain);
RC local_Cauchy_stress_Almansi_strain(ELEMENT element,
                                      const double *local_points,
                                      NODE_ARRAY node, DISP_ARRAY disp,
                                      MATERIAL_PROP_ARRAY material,
                                      PHYSICAL_PROP_ARRAY physical,
                                      WORK_SET ws,
                                      STRESS_STRAIN *stress,
                                      STRESS_STRAIN *strain);
RC local_PK2_stress_Green_strain(ELEMENT element,
                                 const double *local_points,
                                 NODE_ARRAY node, DISP_ARRAY disp,
                                 MATERIAL_PROP_ARRAY material,
                                 PHYSICAL_PROP_ARRAY physical,
                                 WORK_SET ws,
                                 STRESS_STRAIN *stress,
                                 STRESS_STRAIN *strain);
RC local_Almansi_strain(ELEMENT element, const double *local_loint,
                        NODE_ARRAY node, WORK_SET ws, DISP_ARRAY disp,
                        STRESS_STRAIN *Almansi_strain);
RC local_Green_strain(ELEMENT element, const double *local_point,
                      NODE_ARRAY node, WORK_SET ws, DISP_ARRAY disp,
                      STRESS_STRAIN *Green_strain);
RC make_G_matrix(ELEMENT element, const double *g_point,
                 double **node_location, double **G_matrix);
RC set_M_matrix(int dim, double *stress, double **M_matrix);
RC copy_disp2vect(int dim, NODE_ARRAY node, DISP_ARRAY disp, double d_vect[],
                  double factor);
RC copy_vect2disp(int dim, NODE_ARRAY node, const double d_vect[],
                  DISP_ARRAY disp, double factor);
RC add_disp2vect(int dim, NODE_ARRAY node, DISP_ARRAY disp, double d_vect[],
                 double factor);
RC add_vect2disp(int dim, NODE_ARRAY node, const double d_vect[],
                 DISP_ARRAY disp, double factor);

/* eigen_solver.c */
RC fem_eigen_solve_iccg3(NODE_ARRAY node, ELEMENT_ARRAY element,
                         MATERIAL_PROP_ARRAY material,
                         PHYSICAL_PROP_ARRAY physical,
                         double shift, int max_step,
                         int max_e_num, int *e_num,
                         double evals[], DISP_ARRAY mode_disps[]);
RC fem_eigen_block_solve_sky3(NODE_ARRAY node, ELEMENT_ARRAY element,
                              MATERIAL_PROP_ARRAY material,
                              PHYSICAL_PROP_ARRAY physical,
                              BC_ARRAY rest,
                              double shift, int b_size, int max_step,
                              int max_e_num, int *e_num,
                              double evals[], DISP_ARRAY mode_disps[]);

/* fem_solver_shell.c */
RC get_D_matrix_shell(ELEMENT element, MATERIAL_PROP_ARRAY material,
                      PHYSICAL_PROP_ARRAY physical, double **D_matrix);
RC impose_shell_elem_matrix(int plane_dim, ELEMENT element,
                            MATERIAL_PROP_ARRAY material,
                            PHYSICAL_PROP_ARRAY physical,
                            ELEM_MATRIX elem_matrix,
                            ELEM_MATRIX *total_elem_matrix);
RC make_shell_elem_matrix(ELEMENT element, NODE_ARRAY node,
                          MATERIAL_PROP_ARRAY material,
                          PHYSICAL_PROP_ARRAY physical,
                          ELEM_MATRIX *elem_matrix, WORK_SET ws);
RC make_B_matrix_shell(ELEMENT element, MATERIAL_PROP_ARRAY material,
                       PHYSICAL_PROP_ARRAY physical, WORK_SET ws,
                       const double *g_point, double **node_location,
                       double **B_matrix, double *det_J);
RC make_shell_matrix(ELEMENT element, NODE_ARRAY node,
                     MATERIAL_PROP_ARRAY material,
                     PHYSICAL_PROP_ARRAY physical,
                     ELEM_MATRIX *elem_matrix);
RC stress_strain_shell(int sol_ss_flag, double *distance,
                       NODE_ARRAY node, ELEMENT_ARRAY element,
                       DISP_ARRAY disp,
                       MATERIAL_PROP_ARRAY material,
                       PHYSICAL_PROP_ARRAY physical,
                       STRESS_ARRAY stress[], STRAIN_ARRAY strain[]);
RC resultant_stress_strain_shell(int sol_ss_flag, NODE_ARRAY node,
                                 ELEMENT_ARRAY element,
                                 DISP_ARRAY disp,
                                 MATERIAL_PROP_ARRAY material,
                                 PHYSICAL_PROP_ARRAY physical,
                                 STRESS_ARRAY *force,
                                 STRESS_ARRAY *moment,
                                 STRESS_ARRAY *shear_force,
                                 STRAIN_ARRAY *curv);
RC local_stress_strain_shell(double distance[], ELEMENT element,
                             const double *local_point,
                             NODE_ARRAY node, DISP_ARRAY disp,
                             MATERIAL_PROP_ARRAY material,
                             PHYSICAL_PROP_ARRAY physical,
                             WORK_SET ws, STRESS_STRAIN stress[],
                             STRESS_STRAIN strain[]);
RC local_resultant_stress_strain(ELEMENT element,
                                 const double *local_point,
                                 NODE_ARRAY node, DISP_ARRAY disp,
                                 MATERIAL_PROP_ARRAY material,
                                 PHYSICAL_PROP_ARRAY physical,
                                 WORK_SET ws, STRESS_STRAIN *force,
                                 STRESS_STRAIN *moment,
                                 STRESS_STRAIN *shear_force,
                                 STRESS_STRAIN *curv);

/* flow_solver.c */
RC steady_flow_analysis(NODE_ARRAY node, ELEMENT_ARRAY elem,
                        BC_ARRAY bc_velo,BC_ARRAY bc_press, 
                        int velo_order, int press_order,
                        double density, double viscosity,
                        VELO_ARRAY *flow);

RC steady_adjflow_analysis(NODE_ARRAY node, ELEMENT_ARRAY elem,
                           BC_ARRAY bc_velo, BC_ARRAY bc_press,
                           int velo_order, int press_order,
                           double density, double viscosity,
                           VELO_ARRAY flow, VELO_ARRAY *adjoint_flow);

RC flow_solve_for_optimization(double delta_t, int iteration, double Re,
                               NODE_ARRAY node,
                               ELEMENT_ARRAY element,
                               BC_ARRAY bc_velo, BC_ARRAY bc_press,
                               VELO_ARRAY *vp);
RC flow_solve_non_steady_flow(int plot_step, double delta_t, int iteration,
                              double Re, NODE_ARRAY node,
                              ELEMENT_ARRAY element, BC_ARRAY bc_velo,
                              BC_ARRAY bc_press, VELO_ARRAY *vp);

RC flow_adjoint_analysis(double Re, NODE_ARRAY node,
                         ELEMENT_ARRAY element, BC_ARRAY bc,
                         VELO_ARRAY vp, VELO_ARRAY *adjoint_vp);
RC flow_adjoint_analysis2(double Re, NODE_ARRAY node,
                          ELEMENT_ARRAY element, BC_ARRAY bc,
                          VELO_ARRAY vp, VELO_ARRAY *adjoint_vp);

RC strain_velocity(int sol_ss_flag, NODE_ARRAY node,
                   ELEMENT_ARRAY element, VELO_ARRAY velo,
                   STRAIN_ARRAY *velo_strain);
RC pick_out_wall_bc_node(BC_ARRAY bc, NODE_ARRAY *wall_node);
RC change_bc_type4flow(BC_ARRAY *bc);


#endif /* FEM_SOLVER_H */
