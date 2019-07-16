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

/* $Id: avs_utl.h 842 2006-12-04 06:21:20Z sasaoka $ */

#ifndef AVS_UTL_H
#define AVS_UTL_H

#include "rc.h"
#include "fem_struct.h"

/* 1step 用 */
RC output_avs_disp(FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY elem,
                   DISP_ARRAY disp);
RC output_avs_velo(FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY elem,
                   VELO_ARRAY velo);
RC output_avs_elem_density(FILE *fp, NODE_ARRAY node, 
                           ELEMENT_ARRAY elem, VELO_ARRAY density);
RC output_avs_bc_force(FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY elem,
                       BC_ARRAY bc);
RC output_avs_bc_temp(FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY elem,
                      BC_ARRAY temp, double def_temp);
RC output_avs_sens(FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY elem,
                   SENSITIVITY_ARRAY sens);
RC output_avs_stress(FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY elem,
                     STRESS_ARRAY stress);
RC output_avs_material_rho(FILE *fp, NODE_ARRAY node,
                           ELEMENT_ARRAY elem, MATERIAL_PROP_ARRAY mat,
                           PHYSICAL_PROP_ARRAY phys);

/* 複数ステップ用 */

/* 形状固定、データ変化 */
FILE *fopen_avs_step_data(const char *FileName, NODE_ARRAY fix_node,
                          ELEMENT_ARRAY fix_elem);
RC output_avs_step_disp(FILE *fp_avs, DISP_ARRAY disp);
RC output_avs_step_velo(FILE *fp_avs, VELO_ARRAY velo);
RC output_avs_step_elem_density(FILE *fp_avs, VELO_ARRAY density);
RC output_avs_step_bc_force(FILE *fp_avs, NODE_ARRAY node, BC_ARRAY bc);
RC output_avs_step_bc_temp(FILE *fp_avs, NODE_ARRAY node, BC_ARRAY temp,
                           double def_temp);


/* 形状、データ共に変化 */
FILE *fopen_avs_step_datageom(const char *FileName);
RC output_avs_step_disp_model(FILE *fp_avs, NODE_ARRAY node,
                              ELEMENT_ARRAY elem, DISP_ARRAY disp);
RC output_avs_step_velo_model(FILE *fp_avs, NODE_ARRAY node,
                              ELEMENT_ARRAY elem, VELO_ARRAY velo);
RC output_avs_step_sens_model(FILE *fp_avs, NODE_ARRAY node,
                              ELEMENT_ARRAY elem,
                              SENSITIVITY_ARRAY sens);
RC output_avs_step_stress_model(FILE *fp_avs, NODE_ARRAY node,
                                ELEMENT_ARRAY elem,
                                STRESS_ARRAY stress);
RC output_avs_step_material_rho_model(FILE *fp_avs, NODE_ARRAY node,
                                      ELEMENT_ARRAY elem,
                                      MATERIAL_PROP_ARRAY mat,
                                      PHYSICAL_PROP_ARRAY phys);
RC output_avs_step_bc_force_model(FILE *fp_avs, NODE_ARRAY node,
                                  ELEMENT_ARRAY elem, BC_ARRAY bc);

/* mgf format */
RC output_avs_fe_model(FILE *fp_mgf, NODE_ARRAY node,
                       ELEMENT_ARRAY elem);

#endif /* AVS_UTL_H */

/*
memo
*/
