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

/* $Id: fem_solver.h 1022 2012-04-30 05:00:38Z aoyama $ */

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

/* fem_solver_shell.c */
#define MITC_SIZE        (9)
#define MAX_SHELL_NODE   (9)
#define TYING_MAX        (8)

#define SOL_SS_BOTTOM  (1)   /* シェル要素の下面での値 */
#define SOL_SS_MID     (2)   /* シェル要素の中立面での値 */
#define SOL_SS_TOP     (3)   /* シェル要素の上面での値 */

typedef struct {
	int label;   /* label of director vector(if < 0 then unused) */
	VECT3D v1;
	VECT3D v2;
	VECT3D vn;
} DIRECTOR_VECT;

typedef struct {
	int size;
	int alloc_size;
	DIRECTOR_VECT *array;
	ARRAY_SORT_FLAG sort_flag;
} DIRECTOR_VECT_ARRAY;

typedef struct {
	double *vec_N[MITC_SIZE];
	double *vec_T[MITC_SIZE];
	double **mat_T_3[MITC_SIZE];
	double **mat_3_N[MITC_SIZE];
	double **mat_3_5N[MITC_SIZE];
	double **mat_5N_5N[MITC_SIZE];
} MITC_TYING_SET;


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

RC fem_solve_add_density(NODE_ARRAY node, ELEMENT_ARRAY element,
             MATERIAL_PROP_ARRAY material,
             PHYSICAL_PROP_ARRAY physical,
             BC_ARRAY rest, BC_ARRAY force,
             BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
             DISP_ARRAY *disp, STRESS_ARRAY *stress,
             STRAIN_ARRAY *strain, int sol_ss_flag
			 ,double *density);

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
                       NODE_ARRAY node, DISP_ARRAY disp,
                       STRESS_STRAIN *strain);
RC local_points_strain_sub(ELEMENT element, const double local_point[],
                           double **node_location, const double disp_v[],
                           double **ws_mat_6_3N, STRESS_STRAIN *strain);
RC stress_strain_therm(int sol_ss_flag, NODE_ARRAY node,
                       ELEMENT_ARRAY element, DISP_ARRAY disp,
                       MATERIAL_PROP_ARRAY material,
                       PHYSICAL_PROP_ARRAY physical,
                       BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                       double fact,
                       STRESS_ARRAY *stress, STRAIN_ARRAY *strain);
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

RC make_g_matrix_add_density(NODE_ARRAY node, ELEMENT_ARRAY element,
                 MATERIAL_PROP_ARRAY material,
                 PHYSICAL_PROP_ARRAY physical,
                 SKYLINE_MATRIX *g_matrix,
				 double *density);

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

RC make_elem_matrix_add_density(ELEMENT element, NODE_ARRAY node,
                  MATERIAL_PROP_ARRAY material, 
				  PHYSICAL_PROP_ARRAY physical,
                  ELEM_MATRIX *elem_matrix,double *density);

double N_func_ex(int i,double x,double y);

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
RC make_mass_matrix_add_density (ELEMENT element, NODE_ARRAY node,
                  MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
                  ELEM_MATRIX *elem_matrix,double *density);

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
RC set_L_matrix(ELEMENT element, double **node_location, double L[][3]);
RC set_L_matrix_shell_solid(ELEMENT type, double **node_location,
                            double L[][3]);
RC set_L_matrix_beam(ELEMENT element, double **node_location, double L[][3]);
RC add_temp_force1(int dim, NODE_ARRAY node, ELEMENT_ARRAY element,
                    MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
                    BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                    double basis_temp, double *f_vect);
RC fem_solve1(NODE_ARRAY node, ELEMENT_ARRAY element,
               MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
               BC_ARRAY rest, BC_ARRAY force, BC_ARRAY temp,
               DEFAULT_BC_ARRAY def_temp, double basis_temp, DISP_ARRAY *disp,
               STRESS_ARRAY *stress, STRAIN_ARRAY *strain, int sol_ss_flag);
RC fem_solve_thermal(NODE_ARRAY node, ELEMENT_ARRAY element,
                     ELEMENT_ARRAY surf, DEFAULT_BC_ARRAY def_temp,
                     MATERIAL_PROP_ARRAY material,
                     PHYSICAL_PROP_ARRAY physical,
                     BC_ARRAY bc_temp, BC_ARRAY flux,
                     ELEMENT_ARRAY transfer, BC_ARRAY *result_temp);
RC make_K_matrix_thermal(SKYLINE_MATRIX g_matrix1,
                         NODE_ARRAY node, ELEMENT_ARRAY element,
                         ELEMENT_ARRAY surf, ELEMENT_ARRAY transfer,
                         MATERIAL_PROP_ARRAY material,
                         PHYSICAL_PROP_ARRAY physical);
RC make_elem_matrix_transfer(ELEMENT surf, NODE_ARRAY node,
                             ELEMENT_ARRAY transfer, ELEM_MATRIX *elem_matrix);
RC add_therm_vector(NODE_ARRAY node, BC_ARRAY flux, double *f_vect);
RC add_transfer_vector(ELEMENT_ARRAY surf, NODE_ARRAY node,
                       ELEMENT_ARRAY transfer, double *s_vect);
RC fem_solve_thermal_transient(NODE_ARRAY node, ELEMENT_ARRAY element,
                               ELEMENT_ARRAY surf, DEFAULT_BC_ARRAY def_temp,
                               MATERIAL_PROP_ARRAY material,
                               PHYSICAL_PROP_ARRAY physical,
                               BC_ARRAY bc_temp, BC_ARRAY flux,
                               ELEMENT_ARRAY transfer, BC_ARRAY **result_temp,
                               int step_num, double delta_t, double theta);
RC make_K_matrix_transient(SKYLINE_MATRIX g_matrix1, SKYLINE_MATRIX g_matrix2,
                           NODE_ARRAY node, ELEMENT_ARRAY element,
                           ELEMENT_ARRAY surf, ELEMENT_ARRAY transfer,
                           MATERIAL_PROP_ARRAY material, double delta_t,
                           PHYSICAL_PROP_ARRAY physical, double theta);
RC make_elem_matrix_heatcapacity(ELEMENT element, NODE_ARRAY node,
                                 MATERIAL_PROP_ARRAY material,
                                 PHYSICAL_PROP_ARRAY physical,
                                 ELEM_MATRIX *elem_matrix);
RC init_cond_temp(BC_ARRAY *init_temp, DEFAULT_BC_ARRAY def_temp,
                  BC_ARRAY bc_temp, NODE_ARRAY node);
RC add_transient_therm_vector(double *s_vect, SKYLINE_MATRIX g_matrix,
                              NODE_ARRAY node, BC_ARRAY temp);

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
RC make_NL_update_react_force(ELEMENT_ARRAY element, NODE_ARRAY node,
                              MATERIAL_PROP_ARRAY material,
                              PHYSICAL_PROP_ARRAY physical,
                              DISP_ARRAY disp, double f_vect[]);
RC make_tangent_elem_matrix(ELEMENT element, NODE_ARRAY node,
                            MATERIAL_PROP_ARRAY material,
                            PHYSICAL_PROP_ARRAY physical,
                            DISP_ARRAY disp, ELEM_MATRIX *elem_matrix,
                            WORK_SET ws);
RC make_update_tangent_elem_matrix(ELEMENT element, NODE_ARRAY node,
                                   MATERIAL_PROP_ARRAY material,
                                   PHYSICAL_PROP_ARRAY physical,
                                   DISP_ARRAY disp, ELEM_MATRIX *elem_matrix,
                                   WORK_SET ws);
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
                                      STRESS_STRAIN *stress,
                                      STRESS_STRAIN *strain,
                                      WORK_SET ws);
RC local_PK2_stress_Green_strain(ELEMENT element,
                                 const double *local_points,
                                 NODE_ARRAY node, DISP_ARRAY disp,
                                 MATERIAL_PROP_ARRAY material,
                                 PHYSICAL_PROP_ARRAY physical,
                                 STRESS_STRAIN *stress,
                                 STRESS_STRAIN *strain,
                                 WORK_SET ws);
RC local_Almansi_strain(ELEMENT element, const double *local_loint,
                        NODE_ARRAY node, DISP_ARRAY disp,
                        STRESS_STRAIN *Almansi_strain, WORK_SET ws);
RC local_Green_strain(ELEMENT element, const double *local_point,
                      NODE_ARRAY node, DISP_ARRAY disp,
                      STRESS_STRAIN *Green_strain, WORK_SET ws);
RC make_AG_matrix(ELEMENT element, const double *g_point, WORK_SET ws,
                  double **node_location, DISP_ARRAY disp,
                  double **A_matrix, double **G_matrix);
RC make_G_matrix(ELEMENT element, const double *g_point,
                 double **node_location, double **G_matrix);
RC set_A_matrix(int dim, double **duvw, double **A_matrix);
RC set_G_matrix(int dim, int node_num, double **dNxyz, double **G_matrix);
RC set_M_matrix(int dim, double *stress, double **M_matrix);
RC copy_disp2vect(int dim, NODE_ARRAY node, DISP_ARRAY disp, double d_vect[],
                  double factor);
RC copy_vect2disp(int dim, NODE_ARRAY node, const double d_vect[],
                  DISP_ARRAY disp, double factor);
RC add_disp2vect(int dim, NODE_ARRAY node, DISP_ARRAY disp, double d_vect[],
                 double factor);
RC add_vect2disp(int dim, NODE_ARRAY node, const double d_vect[],
                 DISP_ARRAY disp, double factor);
RC update_node_location_matrix(ELEMENT *element, NODE_ARRAY node,
                               DISP_ARRAY disp, double **matrix);
RC fem_hyperelastic_solve(NODE_ARRAY node, ELEMENT_ARRAY element,
                          MATERIAL_PROP_ARRAY material,
                          PHYSICAL_PROP_ARRAY physical,
                          BC_ARRAY rest, BC_ARRAY force,
                          DEFAULT_BC_ARRAY b_force,
                          BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                          DISP_ARRAY *disp, STRESS_ARRAY *stress,
                          STRAIN_ARRAY *strain, int sol_ss_flag);
RC make_hyperelastic_tangent_elem_matrix(ELEMENT element, NODE_ARRAY node,
                                         MATERIAL_PROP_ARRAY material,
                                         PHYSICAL_PROP_ARRAY physical,
                                         DISP_ARRAY disp,
                                         ELEM_MATRIX *elem_matrix,
                                         WORK_SET ws);
RC make_incompressible_elem_matrix(ELEMENT element, NODE_ARRAY node,
                                   MATERIAL_PROP_ARRAY material,
                                   PHYSICAL_PROP_ARRAY physical,
                                   DISP_ARRAY disp,
                                   ELEM_MATRIX *elem_matrix, int offset,
                                   WORK_SET ws);
RC make_hyperelestic_react_force(ELEMENT_ARRAY element, NODE_ARRAY node,
                                 MATERIAL_PROP_ARRAY material,
                                 PHYSICAL_PROP_ARRAY physical,
                                 DISP_ARRAY disp, double f_vect[]);
RC make_incompressible_force(ELEMENT_ARRAY element, NODE_ARRAY node,
                             MATERIAL_PROP_ARRAY material,
                             PHYSICAL_PROP_ARRAY physical,
                             DISP_ARRAY disp, double f_vect[]);
RC hyperelastic_PK2_stress_Green_strain(int sol_ss_flag, NODE_ARRAY node,
                                        ELEMENT_ARRAY element, DISP_ARRAY disp,
                                        MATERIAL_PROP_ARRAY material,
                                        PHYSICAL_PROP_ARRAY physical,
                                        STRESS_ARRAY *stress,
                                        STRAIN_ARRAY *strain);
RC local_hyperelastic_PK2_stress_Green_strain(ELEMENT element,
                                              const double *local_point,
                                              NODE_ARRAY node, DISP_ARRAY disp,
                                              MATERIAL_PROP_ARRAY material,
                                              PHYSICAL_PROP_ARRAY physical,
                                              STRESS_STRAIN *stress,
                                              STRESS_STRAIN *strain,
                                              WORK_SET ws);
RC local_hyperelastic_PK2_elastic_stress_Green_strain(ELEMENT element,
                                                   const double *local_point,
                                                   NODE_ARRAY node,
                                                   DISP_ARRAY disp,
                                                   MATERIAL_PROP_ARRAY material,
                                                   PHYSICAL_PROP_ARRAY physical,
                                                   STRESS_STRAIN *stress,
                                                   STRESS_STRAIN *strain,
                                                   WORK_SET ws);
RC
local_hyperelastic_PK2_hydrostatic_stress_Green_strain(ELEMENT element,
                                                   const double *local_point,
                                                   NODE_ARRAY node,
                                                   DISP_ARRAY disp,
                                                   MATERIAL_PROP_ARRAY material,
                                                   PHYSICAL_PROP_ARRAY physical,
                                                   STRESS_STRAIN *stress,
                                                   STRESS_STRAIN *strain,
                                                   WORK_SET ws);
double first_invariant_3x3tensor(double **A);
double second_invariant_3x3tensor(double **A);
double third_invariant_3x3tensor(double **A);

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
RC fem_solve_shell_iccg6(NODE_ARRAY node, ELEMENT_ARRAY element,
                         RIGID_ELEMENT_ARRAY rigid,
                         MATERIAL_PROP_ARRAY material,
                         PHYSICAL_PROP_ARRAY physical,
                         BC_ARRAY rest, BC_ARRAY force,
                         BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                         DEFAULT_BC_ARRAY def_body_force,
                         DISP_ARRAY *disp,
                         STRESS_ARRAY *stress, STRAIN_ARRAY *strain,
                         int sol_ss_flag, int sol_ss_face);
RC fem_solve_shell(NODE_ARRAY node, ELEMENT_ARRAY element,
                   RIGID_ELEMENT_ARRAY rigid,
                   MATERIAL_PROP_ARRAY material,
                   PHYSICAL_PROP_ARRAY physical,
                   BC_ARRAY rest, BC_ARRAY force,
                   BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                   DEFAULT_BC_ARRAY def_body_force,
                   DISP_ARRAY *disp,
                   STRESS_ARRAY *stress, STRAIN_ARRAY *strain,
                   int sol_ss_flag, int sol_ss_face);
RC modify_g_matrix6(NODE_ARRAY node, BC_ARRAY rest, SKYLINE_MATRIX g_matrix,
                    double f_vect[]);
RC modify_g_matrix6_n(NODE_ARRAY node, BC_ARRAY rest,
                      SKYLINE_MATRIX g_matrix, double **f_vect,
                      int f_vect_size);
RC make_shell_elem_matrix(ELEMENT element, NODE_ARRAY node,
                          MATERIAL_PROP_ARRAY material,
                          PHYSICAL_PROP_ARRAY physical,
                          DIRECTOR_VECT_ARRAY d_vect,
                          ELEM_MATRIX *elem_matrix,
                          MITC_TYING_SET mitc,
                          WORK_SET ws1, WORK_SET ws2);
RC make_shell_mass_matrix(ELEMENT element, NODE_ARRAY node,
                          MATERIAL_PROP_ARRAY material,
                          PHYSICAL_PROP_ARRAY physical,
                          DIRECTOR_VECT_ARRAY d_vect,
                          ELEM_MATRIX *elem_matrix,
                          WORK_SET ws1, WORK_SET ws2);
RC add_body_force_shell(NODE_ARRAY node, ELEMENT_ARRAY element,
                        MATERIAL_PROP_ARRAY material,
                        PHYSICAL_PROP_ARRAY physical,
                        DEFAULT_BC_ARRAY def_b_force,
                        DIRECTOR_VECT_ARRAY d_vect, double f_vect[]);
RC add_temp_force_shell(NODE_ARRAY node, ELEMENT_ARRAY element,
                        MATERIAL_PROP_ARRAY material,
                        PHYSICAL_PROP_ARRAY physical,
                        BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                        DIRECTOR_VECT_ARRAY d_vect, double f_vect[]);
int analysis_dim_shell(ELEMENT_ARRAY element);
int element_dim_shell(ELEM_TYPE type);
int element_ssnum_shell(ELEM_TYPE type);
RC stress_strain_shell(int sol_ss_flag, int sol_ss_face,
                       NODE_ARRAY node, ELEMENT_ARRAY element,
                       MATERIAL_PROP_ARRAY material,
                       PHYSICAL_PROP_ARRAY physical,
                       DIRECTOR_VECT_ARRAY d_vect, DISP_ARRAY disp,
                       STRESS_ARRAY *stress, STRAIN_ARRAY *strain);
RC local_stress_strain_shell(ELEMENT element, const double *local_point,
                             NODE_ARRAY node,
                             MATERIAL_PROP_ARRAY material,
                             PHYSICAL_PROP_ARRAY physical,
                             DIRECTOR_VECT_ARRAY d_vect,
                             DISP_ARRAY disp,
                             STRESS_STRAIN *stress,
                             STRESS_STRAIN *strain,
                             MITC_TYING_SET mitc,
                             WORK_SET ws1, WORK_SET ws2);
RC make_nodal_director_vect(NODE_ARRAY node, ELEMENT_ARRAY element,
                            DIRECTOR_VECT_ARRAY *vect);
RC set_element_direcor_vect(ELEMENT element, DIRECTOR_VECT_ARRAY vect,
                            DIRECTOR_VECT elem_vect[]);
int search_director_vect_label(DIRECTOR_VECT_ARRAY vect, int label);
RC allocate_director_vect_array(int num, DIRECTOR_VECT_ARRAY *vect);
RC free_director_vect_array(DIRECTOR_VECT_ARRAY *vect);
RC copy_director_vect_array(DIRECTOR_VECT_ARRAY src,
                            DIRECTOR_VECT_ARRAY *dest);
RC allocate_mitc_tying_set(MITC_TYING_SET *mitc);
RC free_mitc_tying_set(MITC_TYING_SET *mitc);
RC make_D_matrix_shell(ELEMENT element, MATERIAL_PROP_ARRAY material,
                       PHYSICAL_PROP_ARRAY physical, VECT3D g[],
                       double **D_matrix);
RC set_center_point_shell(int sol_ss_face, ELEM_TYPE type,
                          double *lcoord);
RC set_local_node_points_shell(int sol_ss_face, ELEM_TYPE type,
                               double **lcoord);
RC make_G1G2G3_basis(ELEM_TYPE elem_type, int elem_node_num,
                     DIRECTOR_VECT elem_vect[], double thickness,
                     double **node_location, const double *local_point,
                     VECT3D G[], WORK_SET ws);
RC make_G1G2G3g1g2g3_basis(ELEM_TYPE elem_type, int elem_node_num,
                           DIRECTOR_VECT elem_vect[], double thickness,
                           double **node_location, const double *local_point,
                           VECT3D G[], VECT3D g[], WORK_SET ws);
RC make_g1g2g3_basis(ELEM_TYPE elem_type, int elem_node_num,
                     DIRECTOR_VECT elem_vect[], double thickness,
                     double **node_location, const double *local_point,
                     VECT3D g[], WORK_SET ws);
RC make_B_matrix_shell(ELEM_TYPE elem_type, int elem_node_num,
                       DIRECTOR_VECT elem_vect[], double thickness,
                       double **node_location,
                       const double *local_point,
                       double **B_matrix,
                       MITC_TYING_SET mitc, WORK_SET ws);
RC set_tying_data_mitc3(const double *local_point, int *ErtEst_tying_num,
                        double **Ert_rt_tying_point,
                        double **Ert_st_tying_point,
                        double **Est_rt_tying_point,
                        double **Est_st_tying_point,
                        double *Ert_rt_h, double *Ert_st_h,
                        double *Est_rt_h, double *Est_st_h);
RC set_tying_data_mitc6(const double *local_point,
                        int *ErrEss_tying_num,
                        int *Ers_tying_num, int *ErtEst_tying_num,
                        double **Err_tying_point,
                        double **Ess_tying_point,
                        double **Ers_rr_tying_point,
                        double **Ers_ss_tying_point,
                        double **Ers_rs_tying_point,
                        double **Ert_rt_tying_point,
                        double **Ert_st_tying_point,
                        double **Est_rt_tying_point,
                        double **Est_st_tying_point,
                        double *Err_h, double *Ess_h, double *Ers_rr_h,
                        double *Ers_ss_h, double *Ers_rs_h,
                        double *Ert_rt_h, double *Ert_st_h,
                        double *Est_rt_h, double *Est_st_h);
RC set_tying_data_mitc4(const double *local_point,
                        int *ErtEst_tying_num,
                        double **Ert_tying_point,
                        double **Est_tying_point,
                        double *Ert_h, double *Est_h);
RC set_tying_data_mitc8(const double *local_point, int *Ers_tying_num,
                        int *ErrEssErtEst_tying_num,
                        double **Err_tying_point,
                        double **Ess_tying_point,
                        double **Ers_tying_point,
                        double **Ert_tying_point,
                        double **Est_tying_point,
                        double *Err_h, double *Ess_h, double *Ers_h,
                        double *Ert_h, double *Est_h);
RC set_tying_data_mitc9(const double *local_point, int *Ers_tying_num,
                        int *ErrEssErtEst_tying_num,
                        double **Err_tying_point,
                        double **Ess_tying_point,
                        double **Ers_tying_point,
                        double **Ert_tying_point,
                        double **Est_tying_point,
                        double *Err_h, double *Ess_h, double *Ers_h,
                        double *Ert_h, double *Est_h);
RC Gauss_const_default_shell(ELEM_TYPE type, double **points, double *weight,
                             int *num_point);
RC set_N_vector_shell(ELEM_TYPE type, const double *coord, double *N);
RC set_dN_matrix_shell(ELEM_TYPE type, const double *coord, double **dN);

/* fem_solver_plate_shell.c */
RC fem_solve_plate_shell_iccg6(NODE_ARRAY node, ELEMENT_ARRAY element,
                               RIGID_ELEMENT_ARRAY rigid,
                               MATERIAL_PROP_ARRAY material,
                               PHYSICAL_PROP_ARRAY physical,
                               BC_ARRAY rest, BC_ARRAY force,
                               BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                               DEFAULT_BC_ARRAY b_force,
                               DISP_ARRAY *disp, STRESS_ARRAY *stress,
                               STRAIN_ARRAY *strain,
                               int sol_ss_flag, int sol_ss_face);
RC fem_solve_plate_shell(NODE_ARRAY node, ELEMENT_ARRAY element,
                         RIGID_ELEMENT_ARRAY rigid,
                         MATERIAL_PROP_ARRAY material,
                         PHYSICAL_PROP_ARRAY physical,
                         BC_ARRAY rest, BC_ARRAY force,
                         BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                         DEFAULT_BC_ARRAY b_force,
                         DISP_ARRAY *disp, STRESS_ARRAY *stress,
                         STRAIN_ARRAY *strain,
                         int sol_ss_flag, int sol_ss_face);
RC make_plate_shell_elem_matrix(ELEMENT element, NODE_ARRAY node,
                                MATERIAL_PROP_ARRAY material,
                                PHYSICAL_PROP_ARRAY physical,
                                ELEM_MATRIX *elem_matrix,
                                WORK_SET ws);
RC fill_disp_vect_shell(ELEMENT element, DISP_ARRAY disp, double *disp_v);
RC stress_strain_plate_shell(int sol_ss_flag, int sol_ss_face,
                             NODE_ARRAY node, ELEMENT_ARRAY element,
                             MATERIAL_PROP_ARRAY material,
                             PHYSICAL_PROP_ARRAY physical,
                             DISP_ARRAY disp,
                             STRESS_ARRAY *stress, STRAIN_ARRAY *strain);
RC local_stress_strain_plate_shell(ELEMENT element, const double *local_point,
                                   NODE_ARRAY node,
                                   MATERIAL_PROP_ARRAY material,
                                   PHYSICAL_PROP_ARRAY physical,
                                   DISP_ARRAY disp,
                                   STRESS_STRAIN *stress,
                                   STRESS_STRAIN *strain,
                                   int sol_ss_face, WORK_SET ws);
RC set_trans_matrix4stress_strain (double L[][3], double **Q_matrix);
RC make_B_matrix_plate_shell(ELEMENT element, const double *local_point,
                             NODE_ARRAY node,
                             double **Bm_matrix, double **Bb_matrix,
                             double **Bs_matrix, double *det_J);
RC make_plate_shell_mass_matrix(ELEMENT element, NODE_ARRAY node,
                                MATERIAL_PROP_ARRAY material,
                                PHYSICAL_PROP_ARRAY physical,
                                ELEM_MATRIX *elem_matrix,
                                WORK_SET ws);

/* nl_solver_shell.c */
RC fem_NL_solve_shell_iccg6LCM(NODE_ARRAY node, ELEMENT_ARRAY element,
                               RIGID_ELEMENT_ARRAY rigid,
                               MATERIAL_PROP_ARRAY material,
                               PHYSICAL_PROP_ARRAY physical,
                               BC_ARRAY rest, BC_ARRAY force,
                               BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                               DEFAULT_BC_ARRAY def_body_force,
                               DIRECTOR_VECT_ARRAY *initial_vect,
                               DIRECTOR_VECT_ARRAY *current_vect,
                               DISP_ARRAY *disp,
                               STRESS_ARRAY *stress, STRAIN_ARRAY *strain,
                               int sol_ss_flag, int sol_ss_face, int lcm_num);
RC fem_NL_solve_shell_LCM(NODE_ARRAY node, ELEMENT_ARRAY element,
                          RIGID_ELEMENT_ARRAY rigid,
                          MATERIAL_PROP_ARRAY material,
                          PHYSICAL_PROP_ARRAY physical,
                          BC_ARRAY rest, BC_ARRAY force,
                          BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                          DEFAULT_BC_ARRAY def_body_force,
                          DIRECTOR_VECT_ARRAY *initial_vect,
                          DIRECTOR_VECT_ARRAY *current_vect,
                          DISP_ARRAY *disp,
                          STRESS_ARRAY *stress, STRAIN_ARRAY *strain,
                          int sol_ss_flag, int sol_ss_face, int lcm_num);
RC make_NL_shell_react_force(ELEMENT_ARRAY element, NODE_ARRAY node,
                             MATERIAL_PROP_ARRAY material,
                             PHYSICAL_PROP_ARRAY physical,
                             DIRECTOR_VECT_ARRAY initial_vect,
                             DIRECTOR_VECT_ARRAY current_vect,
                             DISP_ARRAY disp, double react_vect[]);
RC make_NL_shell_elem_react_force(ELEMENT element, NODE_ARRAY node,
                                  MATERIAL_PROP_ARRAY material,
                                  PHYSICAL_PROP_ARRAY physical,
                                  DIRECTOR_VECT_ARRAY initial_vect,
                                  DIRECTOR_VECT_ARRAY current_vect,
                                  DISP_ARRAY disp,
                                  int elem_force_index[], double elem_force[],
                                  MITC_TYING_SET mitc,
                                  WORK_SET ws1, WORK_SET ws2);
RC make_NL_shell_tangent_elem_matrix(ELEMENT element, NODE_ARRAY node,
                                     MATERIAL_PROP_ARRAY material,
                                     PHYSICAL_PROP_ARRAY physical,
                                     DIRECTOR_VECT_ARRAY initial_vect,
                                     DIRECTOR_VECT_ARRAY current_vect,
                                     DISP_ARRAY disp,
                                     ELEM_MATRIX *elem_matrix,
                                     MITC_TYING_SET mitc,
                                     WORK_SET ws1, WORK_SET ws2);
RC update_current_director_vect(DISP_ARRAY delta_disp,
                                DIRECTOR_VECT_ARRAY *current_vect);
RC make_large_rotation_matrix(ELEM_TYPE elem_type, int elem_node_num,
                              DIRECTOR_VECT elem_curr_vect[],
                              double thickness,
                              double **update_location,
                              const double *local_point,
                              double stress_vect[],
                              double *Kex,
                              MITC_TYING_SET mitc, WORK_SET ws);
RC make_SZtZ_matrix(ELEM_TYPE elem_type, int elem_node_num,
                    DIRECTOR_VECT elem_curr_vect[],
                    double thickness, const double *local_point,
                    double stress_vect[], double **SZtZ_matrix,
                    MITC_TYING_SET mitc, WORK_SET ws);
RC local_covariant_Green_strain_vect(ELEM_TYPE elem_type,
                                     int elem_node_num,
                                     DIRECTOR_VECT elem_init_vect[],
                                     DIRECTOR_VECT elem_curr_vect[],
                                     double thickness,
                                     double **node_location,
                                     double **update_location,
                                     const double *local_point,
                                     double strain_vect[],
                                     MITC_TYING_SET mitc, WORK_SET ws);
RC local_PK2_stress_Green_strain_shell(ELEMENT element,
                                       const double *local_point,
                                       NODE_ARRAY node,
                                       MATERIAL_PROP_ARRAY material,
                                       PHYSICAL_PROP_ARRAY physical,
                                       DIRECTOR_VECT_ARRAY initial_vect,
                                       DIRECTOR_VECT_ARRAY current_vect,
                                       DISP_ARRAY disp,
                                       STRESS_STRAIN *stress,
                                       STRESS_STRAIN *strain,
                                       MITC_TYING_SET mitc,
                                       WORK_SET ws1, WORK_SET ws2);
RC PK2_stress_Green_strain_shell(int sol_ss_flag, int sol_ss_face,
                                 NODE_ARRAY node, ELEMENT_ARRAY element,
                                 MATERIAL_PROP_ARRAY material,
                                 PHYSICAL_PROP_ARRAY physical,
                                 DIRECTOR_VECT_ARRAY initial_vect,
                                 DIRECTOR_VECT_ARRAY current_vect,
                                 DISP_ARRAY disp,
                                 STRESS_ARRAY *stress, STRAIN_ARRAY *strain);

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
