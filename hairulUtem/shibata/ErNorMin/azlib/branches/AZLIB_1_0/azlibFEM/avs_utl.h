/*********************************************************************
 * avs_utl.h
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Syoji ITO> <Tomoyuki TSUBATA> <Yasuo ASAGA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: avs_utl.h,v 1.9 2003/10/07 11:41:32 sasaoka Exp $ */

#ifndef AVS_UTL_H
#define AVS_UTL_H

#include "rc.h"
#include "fem_struct.h"

/* 1step 用 */
RC output_avs_disp(FILE *fp, FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY elem,
                   FEM_DISP_ARRAY disp);
RC output_avs_velo(FILE *fp, FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY elem,
                   FEM_VELO_ARRAY velo);
RC output_avs_elem_density(FILE *fp, FEM_NODE_ARRAY node, 
                           FEM_ELEMENT_ARRAY elem, FEM_VELO_ARRAY density);
RC output_avs_bc_force(FILE *fp, FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY elem,
                       FEM_BC_ARRAY bc);
RC output_avs_bc_temp(FILE *fp, FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY elem,
                      FEM_BC_ARRAY temp, double def_temp);
RC output_avs_sens(FILE *fp, FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY elem,
                   FEM_SENSITIVITY_ARRAY sens);
RC output_avs_stress(FILE *fp, FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY elem,
                     FEM_STRESS_ARRAY stress);
RC output_avs_material_rho(FILE *fp, FEM_NODE_ARRAY node,
                           FEM_ELEMENT_ARRAY elem, FEM_MATERIAL_PROP_ARRAY mat,
                           FEM_PHYSICAL_PROP_ARRAY phys);

/* 複数ステップ用 */

/* 形状固定、データ変化 */
FILE *fopen_avs_step_data(const char *FileName, FEM_NODE_ARRAY fix_node,
                          FEM_ELEMENT_ARRAY fix_elem);
RC output_avs_step_disp(FILE *fp_avs, FEM_DISP_ARRAY disp);
RC output_avs_step_velo(FILE *fp_avs, FEM_VELO_ARRAY velo);
RC output_avs_step_elem_density(FILE *fp_avs, FEM_VELO_ARRAY density);
RC output_avs_step_bc_force(FILE *fp_avs, FEM_NODE_ARRAY node, FEM_BC_ARRAY bc);


/* 形状、データ共に変化 */
FILE *fopen_avs_step_datageom(const char *FileName);
RC output_avs_step_disp_model(FILE *fp_avs, FEM_NODE_ARRAY node,
                              FEM_ELEMENT_ARRAY elem, FEM_DISP_ARRAY disp);
RC output_avs_step_velo_model(FILE *fp_avs, FEM_NODE_ARRAY node,
                              FEM_ELEMENT_ARRAY elem, FEM_VELO_ARRAY velo);
RC output_avs_step_sens_model(FILE *fp_avs, FEM_NODE_ARRAY node,
                              FEM_ELEMENT_ARRAY elem,
                              FEM_SENSITIVITY_ARRAY sens);
RC output_avs_step_stress_model(FILE *fp_avs, FEM_NODE_ARRAY node,
                                FEM_ELEMENT_ARRAY elem,
                                FEM_STRESS_ARRAY stress);
RC output_avs_step_material_rho_model(FILE *fp_avs, FEM_NODE_ARRAY node,
                                      FEM_ELEMENT_ARRAY elem,
                                      FEM_MATERIAL_PROP_ARRAY mat,
                                      FEM_PHYSICAL_PROP_ARRAY phys);
RC output_avs_step_bc_force_model(FILE *fp_avs, FEM_NODE_ARRAY node,
                                  FEM_ELEMENT_ARRAY elem, FEM_BC_ARRAY bc);

/* mgf format */
RC output_avs_fe_model(FILE *fp_mgf, FEM_NODE_ARRAY node,
                       FEM_ELEMENT_ARRAY elem);

#endif /* AVS_UTL_H */

/*
memo
*/
