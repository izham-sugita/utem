#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fem.h"
#include "mathematics.h"
#include "base.h"

#define JOB_INIT (-1) 
#define JOB_END (-2) 
#define DIM_DENSITY (1) /* 密度変動は1自由度 */ 
#define DENS_UPPER_LIMIT (1.0)     /* phi^pの上限値（phi_1) */
#define UPPER_DENS_PHI (1.0)
//#define INIT_DENS (0.0005)     /* 初期のphiの値（0<=phi<=phi_1) */
#define INIT_DENS (0.9999999)     /* 初期のphiの値（0<=phi<=phi_1) */

typedef struct {
	int    it_max;
//	int    nr_max;
	int    sol_ss_flag;
	double dens_pow;
	double c0;
	double c1;
	double max_dens;
//	double subj_ratio;
//	double subj_rel_error;
//	double max_divisor;
	double obj_rel_error;
} PARAM;

/* main.c */
double node_weight(BC_ARRAY weight, int node_label);
RC nl_strain_analysis(NODE_ARRAY node, ELEMENT_ARRAY element,
                      MATERIAL_PROP_ARRAY material,
                      PHYSICAL_PROP_ARRAY physical,
                      BC_ARRAY rest, BC_ARRAY dens_rest, BC_ARRAY force,
                      BC_ARRAY density, BC_ARRAY alpha,
                      NODE_ARRAY target_node, ELEMENT_ARRAY target_surf,
                      DISP_ARRAY *disp, STRAIN_ARRAY *strain, PARAM param);
RC cal_obj_func(NODE_ARRAY node, ELEMENT_ARRAY element,
                MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
                BC_ARRAY rest, BC_ARRAY force, BC_ARRAY density, BC_ARRAY alpha,
                NODE_ARRAY target_node, ELEMENT_ARRAY target_surf,
                DISP_ARRAY *disp, STRAIN_ARRAY *strain,
                double *obj, BC_ARRAY *sens_trac, PARAM param);
RC integral_error_norm(NODE_ARRAY node, NODE_ARRAY target_node,
                       ELEMENT_ARRAY target_surf, DISP_ARRAY disp,
                       BC_ARRAY alpha, double *obj);
RC adj_NL_solve_topology(NODE_ARRAY node, ELEMENT_ARRAY element,
                         MATERIAL_PROP_ARRAY material,
                         PHYSICAL_PROP_ARRAY physical,
                         DISP_ARRAY disp, BC_ARRAY rest, BC_ARRAY alpha,
                         BC_ARRAY density, double dens_pow,
                         NODE_ARRAY target_node, ELEMENT_ARRAY target_surf,
                         DISP_ARRAY *adj_disp);
RC make_adj_force(int dim, NODE_ARRAY node, NODE_ARRAY target_node,
                  ELEMENT_ARRAY target_surf, DISP_ARRAY disp,
                  BC_ARRAY alpha, double *f_vect);
RC make_adj_nonzero1a_ccs_sym (NODE_ARRAY node, ELEMENT_ARRAY element,
                               MATERIAL_PROP_ARRAY material,
                               PHYSICAL_PROP_ARRAY physical,
                               DISP_ARRAY disp, BC_ARRAY alpha,
                               BC_ARRAY density, double dens_pow,
                               ELEMENT_ARRAY target_surf,
                               NONZERO_MATRIX_CCS *K);
RC make_NtN_elem_matrix(int dim, ELEMENT surf_elem, NODE_ARRAY node,
                        PHYSICAL_PROP_ARRAY physical, BC_ARRAY alpha,
                        WORK_SET ws, ELEM_MATRIX *elem_matrix);
RC make_robin_react_force(ELEMENT_ARRAY surf, NODE_ARRAY node,
                          MATERIAL_PROP_ARRAY material,
                          PHYSICAL_PROP_ARRAY physical,
                          BC_ARRAY weight, DISP_ARRAY disp,
                          double robin_force[]);
RC sens_trac_obj(NODE_ARRAY node, ELEMENT_ARRAY element,
                 DISP_ARRAY disp, DISP_ARRAY adj_disp,
                 MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
                 BC_ARRAY density, double dens_pow,
                 BC_ARRAY *sens_trac);
RC local_delta_Green_strain(ELEMENT element, const double *local_points,
                            NODE_ARRAY node, DISP_ARRAY disp,
                            DISP_ARRAY adj_disp,
                            WORK_SET ws,
                            STRESS_STRAIN *Green_strain);
RC velocity_analysis_topology (NODE_ARRAY node, ELEMENT_ARRAY element,
                               PHYSICAL_PROP_ARRAY physical,
                               BC_ARRAY dens_rest, int num_case,
                               const BC_ARRAY tractions[], BC_ARRAY density[],
                               PARAM param);
RC make_ccs_matrix_topology (NODE_ARRAY node, ELEMENT_ARRAY element,
                             PHYSICAL_PROP_ARRAY physical,
                             NONZERO_MATRIX_CCS *matrix, PARAM param);
RC make_H1_elem_matrix_topology(ELEMENT element, NODE_ARRAY node,
                                PHYSICAL_PROP_ARRAY physical,
                                ELEM_MATRIX *elem_matrix,
                                WORK_SET ws, double c0, double c1);
RC modify_ccs_format_n (int dim, NODE_ARRAY node, BC_ARRAY rest, CCS *ccs,
                        double **f_vect, int vect_size);
RC modify_ccs_format_dof_n (int index_xyz, double rest_value, CCS *ccs,
                            double **f_vect, int vect_size);
RC vect2bc_density (NODE_ARRAY node, const double *d_vect, BC_ARRAY *density);
RC functional_variation_topology(BC_ARRAY sens_force, BC_ARRAY sens_dens,
                                 double *delta_f);
RC max_density(int dens_num, const double factors[],
               NODE_ARRAY node, ELEMENT_ARRAY element,
               const BC_ARRAY density[], double *max_dens);
RC add_bc_density(BC_ARRAY density_total, int dens_num,
                  const double *factors, const BC_ARRAY *densities);
int chk_armijo (double old_obj, double new_obj, double G_rho, double delta_s,
                double *divisor, PARAM param);


/* nl_solver_topology.c */
RC fem_NL_solve_mumps_topology (NODE_ARRAY node, ELEMENT_ARRAY element,
                                MATERIAL_PROP_ARRAY material,
                                PHYSICAL_PROP_ARRAY physical,
                                BC_ARRAY rest, BC_ARRAY force, BC_ARRAY density,
                                DEFAULT_BC_ARRAY b_force,
                                BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                                DISP_ARRAY *disp, STRESS_ARRAY *stress,
                                STRAIN_ARRAY *strain, double dens_pow,
                                int sol_ss_flag);
RC make_nonzero1a_ccs_sym_topology (NODE_ARRAY node, ELEMENT_ARRAY element,
                                    MATERIAL_PROP_ARRAY material,
                                    PHYSICAL_PROP_ARRAY physical,
                                    BC_ARRAY density, double dens_pow,
                                    NONZERO_MATRIX_CCS *matrix);
RC make_NL_nonzero1a_ccs_sym_topology (NODE_ARRAY node, ELEMENT_ARRAY element,
                                       MATERIAL_PROP_ARRAY material,
                                       PHYSICAL_PROP_ARRAY physical,
                                       DISP_ARRAY disp,
                                       BC_ARRAY density, double dens_pow,
                                       NONZERO_MATRIX_CCS *K);
RC fem_NL_solve_mumps_topology_robin_LCM (NODE_ARRAY node,
                                      ELEMENT_ARRAY element,
                                      ELEMENT_ARRAY surf,
                                      MATERIAL_PROP_ARRAY material,
                                      PHYSICAL_PROP_ARRAY physical,
                                      BC_ARRAY rest, BC_ARRAY force,
                                      BC_ARRAY density, BC_ARRAY weight,
                                      DEFAULT_BC_ARRAY b_force,
                                      BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                                      DISP_ARRAY *disp, STRESS_ARRAY *stress,
                                      STRAIN_ARRAY *strain, double dens_pow,
                                      int sol_ss_flag, int lcm_num);
RC fem_NL_solve_mumps_topology_robin (NODE_ARRAY node, ELEMENT_ARRAY element,
                                      ELEMENT_ARRAY surf,
                                      MATERIAL_PROP_ARRAY material,
                                      PHYSICAL_PROP_ARRAY physical,
                                      BC_ARRAY rest, BC_ARRAY force,
                                      BC_ARRAY density, BC_ARRAY weight,
                                      DEFAULT_BC_ARRAY b_force,
                                      BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                                      DISP_ARRAY *disp, STRESS_ARRAY *stress,
                                      STRAIN_ARRAY *strain, double dens_pow,
                                      int sol_ss_flag);
RC make_nonzero1a_ccs_sym_topology_robin (NODE_ARRAY node,
                                          ELEMENT_ARRAY element,
                                          ELEMENT_ARRAY surf,
                                          MATERIAL_PROP_ARRAY material,
                                          PHYSICAL_PROP_ARRAY physical,
                                          BC_ARRAY density, BC_ARRAY weight,
                                          double dens_pow,
                                          NONZERO_MATRIX_CCS *matrix);
RC make_NL_nonzero1a_ccs_sym_topology_robin (NODE_ARRAY node,
                                             ELEMENT_ARRAY element,
                                             ELEMENT_ARRAY surf,
                                             MATERIAL_PROP_ARRAY material,
                                             PHYSICAL_PROP_ARRAY physical,
                                             DISP_ARRAY disp, BC_ARRAY density,
                                             BC_ARRAY weight, double dens_pow,
                                             NONZERO_MATRIX_CCS *K);
RC make_elem_matrix_topology (ELEMENT element, NODE_ARRAY node,
                              MATERIAL_PROP_ARRAY material,
                              PHYSICAL_PROP_ARRAY physical,
                              BC_ARRAY density, double dens_pow,
                              ELEM_MATRIX *elem_matrix);
RC make_tangent_elem_matrix_topology (ELEMENT element, NODE_ARRAY node,
                                      MATERIAL_PROP_ARRAY material,
                                      PHYSICAL_PROP_ARRAY physical,
                                      WORK_SET ws, DISP_ARRAY disp,
                                      BC_ARRAY density, double dens_pow,
                                      ELEM_MATRIX *elem_matrix);
RC make_NL_react_force_topology (ELEMENT_ARRAY element, NODE_ARRAY node,
                                 MATERIAL_PROP_ARRAY material,
                                 PHYSICAL_PROP_ARRAY physical,
                                 BC_ARRAY density, double dens_pow,
                                 DISP_ARRAY disp, double f_vect[]);
double node_density (BC_ARRAY density, int node_label);
