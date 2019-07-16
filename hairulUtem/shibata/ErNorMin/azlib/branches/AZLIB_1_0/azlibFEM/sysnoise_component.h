/*********************************************************************
 * sysnoise_component.h
 *
 * Copyright (C) 2001 AzLib Developers Group
 *
 * Written by
 *  <Yasuhiro MATSUURA> <Masanobu SUYAMA>
 *
 * Thanks to following contributors.
 *  <Kenzen TAKEUCHI> <Tomoyuki TSUBATA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: sysnoise_component.h,v 1.8 2003/07/23 04:57:45 sasaoka Exp $ */

#ifndef SYSNOISE_COMPONENT_H
#define SYSNOISE_COMPONENT_H

#include <float.h>
#include <limits.h>
#include "rc.h"
#include "fem_struct.h"

#define FSYSNOISEDEFAULT (-(DBL_MAX/2.0))

RC sysnoise_input_pressure( FILE *fp, FEM_SOUND_PRESSURE_ARRAY *sound );
RC sysnoise_input_node( FILE *fp, FEM_NODE_ARRAY *node );
RC sysnoise_input_element( FILE *fp, FEM_ELEMENT_ARRAY *element );
RC sysnoise_output_model_file( FILE *fp, FEM_NODE_ARRAY node,
                               FEM_ELEMENT_ARRAY element );
RC sysnoise_output_option( FILE *fp, char *name, int model_num, int dim );
RC sysnoise_output_point( FILE *fp );
RC sysnoise_output_parameter( FILE *fp, int a_num );
RC sysnoise_output_solve_pres( FILE *fp, char *name, double minfreq,
                               double maxfreq, double linstep, int a_num );
RC sysnoise_output_solve_freq_data( FILE *fp, char *name, double freq,
                                    int counter, int a_num );
RC sysnoise_output_solve_mode( FILE *fp, char *name, double min, double max );
RC sysnoise_output_structure_freq( FILE *fp, double minfreq, double maxfreq,
                                   double linestep );
RC sysnoise_output_structure_mesh( FILE *fp, char *name );
RC sysnoise_output_structure_data( FILE *fp, char *name );
RC sysnoise_output_set_face( FILE *fp, char *name, int set_num,
                             FEM_SOUND_PRESSURE_ARRAY sound );
RC sysnoise_output_structure_option( FILE *fp, char *name, int model_num,
                                     double thickness );
RC sysnoise_output_bc_option( FILE *fp, int set_num );
RC sysnoise_output_vector_mode( FILE *fp, char *name, int mode_num,
                                int damping_flag );
RC sysnoise_output_load_option( FILE *fp, FEM_BC_ARRAY force );
RC sysnoise_output_save_option( FILE *fp, char *name );
RC sysnoise_output_end_option( FILE *fp, char *name );
RC sysnoise_output_link_option( FILE *fp, int struct_mode_num,
                                int fluid_mode_num, int set_num, int a_num );
RC sysnoise_output_file_end( FILE *fp, char *name );
RC sysnoise_input_source_point( FILE *fp, FEM_NODE_ARRAY *node,
                                FEM_SOUND_PRESSURE_ARRAY *sound );
RC sysnoise_output_source_point( FILE *fp, FEM_NODE_ARRAY node,
                                 FEM_SOUND_PRESSURE_ARRAY sound, int dim );
RC sysnoise_input_boundary_pressure( FILE *fp,
                                     FEM_SOUND_PRESSURE_ARRAY *sound );
RC sysnoise_output_boundary_pressure( FILE *fp,
                                      FEM_SOUND_PRESSURE_ARRAY sound );
RC sysnoise_input_boundary_velocity( FILE *fp,
                                     FEM_SOUND_PRESSURE_ARRAY *sound );
RC sysnoise_output_boundary_velocity( FILE *fp,
                                      FEM_SOUND_PRESSURE_ARRAY sound );
RC sysnoise_input_absorbent_panels( FILE *fp, FEM_SOUND_PRESSURE_ARRAY *sound );
RC sysnoise_output_absorbent_panels( FILE *fp, FEM_SOUND_PRESSURE_ARRAY sound );
RC sysnoise_input_symmetry_plane( FILE *fp, TRANSLATION3D *p, int *counter );
RC sysnoise_output_symmetry_plane( FILE *fp, TRANSLATION3D p, int counter );
RC sysnoise_input_eigen_frequency( FILE *fp, double *eigen_freq );
int count_eigen_frequency( FILE *fp );
RC sysnoise_input_couple_eigen_frequency( FILE *fp, double *eigen_freq,
                                          double min, double max );
int count_couple_eigen_frequency( FILE *fp, double min, double max );

#endif /* SYSNOISE_COMPNENT_H */


