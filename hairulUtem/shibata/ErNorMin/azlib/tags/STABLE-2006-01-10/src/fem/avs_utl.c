/*********************************************************************
  iontact.obj \
 * avs_utl.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Syoji ITO> <Tomoyuki TSUBATA> <Yasuo ASAGA>
 *
 * Thanks to following contributors.
 *  <Takaaki NAGATANI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: avs_utl.c 415 2005-07-15 06:04:24Z sasaoka $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "rc.h"
#include "avs_utl.h"
#include "fem_solver.h"
#include "fem_struct.h"

/* ULONG_MAX の桁数 */
#define INT_PLACE_NUMBER ( (int)log10( (double)ULONG_MAX ) + 1 )

/* local functions */
static RC avs_printelem_type(FILE *fp, ELEM_TYPE type);
static RC avs_printelement(FILE *fp, int label, ELEM_TYPE type,
                        int node_num, const int *node);
static RC avs_printfe_model(FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY elem);
static RC avs_printfe_model4discont_data(FILE *fp, NODE_ARRAY node,
                                       ELEMENT_ARRAY elem);
static RC avs_printvelo(FILE *fp, VELO_ARRAY velo);
static RC avs_printdisp(FILE *fp, DISP_ARRAY disp);
static RC avs_printelem_density(FILE *fp, VELO_ARRAY density);
static RC avs_printbc_force(FILE *fp, BC_ARRAY bc_force, NODE_ARRAY node);
static RC avs_printbc_temp(FILE *fp, BC_ARRAY temp, double def_temp,
                        NODE_ARRAY node);
static RC avs_printmaterial_rho(FILE *fp, ELEMENT_ARRAY elem,
                             MATERIAL_PROP_ARRAY mat,
                             PHYSICAL_PROP_ARRAY phys);
static RC avs_printsens(FILE *fp, ELEMENT_ARRAY elem,
                     SENSITIVITY_ARRAY sens);
static RC avs_printstress(FILE *fp, ELEMENT_ARRAY elem,
                       STRESS_ARRAY stress);
static int check_step(FILE *fp);
static RC output_avs_fe_model_face(FILE *fp, NODE_ARRAY node,
                                   ELEMENT_ARRAY elem);
static RC output_avs_fe_model_line(FILE *fp, NODE_ARRAY node,
                                   ELEMENT_ARRAY elem);

static RC
avs_printelem_type (FILE *fp, ELEM_TYPE type)
{
	switch(type){
	case ELEM_LINE1:
		fprintf(fp, "line   ");
		break;
	case ELEM_LINE2:
		fprintf(fp, "line2  ");
		break;
	case ELEM_TRI1:
		fprintf(fp, "tri    ");
		break;
	case ELEM_TRI2:
		fprintf(fp, "tri2   ");
		break;
	case ELEM_QUAD1:
		fprintf(fp, "quad   ");
		break;
	case ELEM_QUAD2:
		fprintf(fp, "quad2  ");
		break;
	case ELEM_TETRA1:
		fprintf(fp, "tet    ");
		break;
	case ELEM_TETRA2:
		fprintf(fp, "tet2   ");
		break;
	case ELEM_PENTA1:
		fprintf(fp, "prism  ");
		break;
	case ELEM_PENTA2:
		fprintf(fp, "prism2 ");
		break;
	case ELEM_HEXA1:
		fprintf(fp, "hex    ");
		break;
	case ELEM_HEXA2:
		fprintf(fp, "hex2   ");
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
avs_printelement (FILE *fp, int label, ELEM_TYPE type,
                  int node_num, const int *node)
{
	int ii1;

	fprintf(fp, "%10d %d ", label, 1);

	RC_TRY( avs_printelem_type(fp, type) ); 

	switch(type){
	case ELEM_TRI2:
		fprintf(fp, "%10d %10d %10d %10d %10d %10d\n",
		        node[0], node[1], node[2], node[5], node[3], node[4]);
		break;
	case ELEM_TETRA2:
		fprintf(fp, "%10d %10d %10d %10d %10d %10d %10d %10d %10d %10d\n",
		        node[0], node[1], node[2], node[3], node[6],
		        node[5], node[7], node[4], node[9], node[8]);
		break;
	case ELEM_PENTA2:
		fprintf(fp, "%10d %10d %10d %10d %10d %10d %10d %10d %10d %10d"
		            "%10d %10d %10d %10d %10d\n",
		        node[0], node[1],  node[2],  node[3],  node[4],
		        node[5], node[8],  node[6],  node[7],  node[11],
		        node[9], node[10], node[12], node[13], node[14]);
		break;
	default:
		for(ii1=0; ii1<node_num; ii1++){
			fprintf(fp, "%10d ", node[ii1]);
		}
		fprintf(fp,"\n");
	}

	return(NORMAL_RC);
}


static RC
avs_printfe_model (FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY elem)
{
	int ii1;

	fprintf( fp, "%d %d\n", count_valid_node(node), count_valid_element(elem) );

	if( node.size <= 0 ) return(WRITE_ERROR_RC);
	for(ii1=0; ii1<node.size; ii1++){
		if(node.array[ii1].label < 0) continue;

		fprintf(fp, "%10d %15.7E %15.7E %15.7E\n", node.array[ii1].label,
		        node.array[ii1].p.x, node.array[ii1].p.y,
		        node.array[ii1].p.z);
	}

	if( elem.size <= 0 ) return(WRITE_ERROR_RC);
	for(ii1=0; ii1<elem.size; ii1++){
		if(elem.array[ii1].label < 0) continue;

		RC_TRY( avs_printelement(fp, elem.array[ii1].label,
		                          elem.array[ii1].type,
		                          elem.array[ii1].node_num,
		                          elem.array[ii1].node) );
	}

	return(NORMAL_RC);
}


static RC
avs_printfe_model4discont_data (FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY elem)
{
	int ii1, ii2;
	int index;
	int counter = 0;
	int new_node_size = 0;
	int elem_node[MAX_NODE];

	for(ii1=0; ii1<elem.size; ii1++){
		if(elem.array[ii1].label < 0) continue;
		new_node_size += elem.array[ii1].node_num;
	}
	fprintf(fp,"%d %d\n", new_node_size, count_valid_element(elem));

	for(ii1=0; ii1<elem.size; ii1++){
		if(elem.array[ii1].label < 0) continue;
		for(ii2=0; ii2<elem.array[ii1].node_num; ii2++){
			counter ++;
			index = search_node_label(node, elem.array[ii1].node[ii2]);
			RC_NEG_CHK(index);
			fprintf(fp, "%10d %15.7e %15.7e %15.7e\n", counter,
			                                           node.array[index].p.x,
			                                           node.array[index].p.y,
			                                           node.array[index].p.z);
		}
	}

	counter = 0;
	for(ii1=0; ii1<elem.size; ii1++){
		if(elem.array[ii1].label < 0) continue;

		for(ii2=0; ii2<elem.array[ii1].node_num; ii2++){
			counter ++;
			elem_node[ii2] = counter;
		}
		RC_TRY( avs_printelement(fp, elem.array[ii1].label,
		                          elem.array[ii1].type,
		                          elem.array[ii1].node_num,
		                          elem_node) );
	}

	return(NORMAL_RC);
}


static RC
avs_printvelo (FILE *fp, VELO_ARRAY velo)
{
	int ii1;
	
	fprintf(fp, "4 0\n");
	
	fprintf(fp, "4 1 1 1 1\n");
	fprintf(fp, "velocity_x,\n");
	fprintf(fp, "velocity_y,\n");
	fprintf(fp, "velocity_z,\n");
	fprintf(fp, "pressure,\n");
		
	for(ii1=0; ii1<velo.size; ii1++){
		if(velo.array[ii1].node < 0) continue;
		fprintf(fp, "%10d %15.7E %15.7E %15.7E %15.7E\n", velo.array[ii1].node,
		                                                  velo.array[ii1].v.x,
		                                                  velo.array[ii1].v.y,
		                                                  velo.array[ii1].v.z,
		                                                  velo.array[ii1].s);
	}

	return(NORMAL_RC);
}


static RC
avs_printdisp (FILE *fp, DISP_ARRAY disp)
{
	int ii1;
	
	fprintf(fp, "3 0\n");
	fprintf(fp, "3 1 1 1\n");
	fprintf(fp, "displacement_x,\n");
	fprintf(fp, "displacement_y,\n");
	fprintf(fp, "displacement_z,\n");

	for(ii1=0; ii1<disp.size; ii1++){
		if(disp.array[ii1].node < 0) continue;
		fprintf(fp, "%10d %15.7E %15.7E %15.7E\n", disp.array[ii1].node,
		                                           disp.array[ii1].v.x,
		                                           disp.array[ii1].v.y,
		                                           disp.array[ii1].v.z);
	}

	return(NORMAL_RC);
}


static RC
avs_printelem_density (FILE *fp, VELO_ARRAY density)
{
	int ii1;

	fprintf(fp, "0 1\n");
	fprintf(fp, "1 1\n");
	fprintf(fp, "density(element data),\n");

	for(ii1=0; ii1<density.size; ii1++){
		if(density.array[ii1].node < 0) continue;
		fprintf(fp, "%10d %15.7E\n", density.array[ii1].node,
		                             density.array[ii1].s);
	}

	return(NORMAL_RC);
}


static RC
avs_printmaterial_rho (FILE *fp, ELEMENT_ARRAY elem,
                       MATERIAL_PROP_ARRAY mat, PHYSICAL_PROP_ARRAY phys)
{
	int ii1, ii2;
	int counter;
	int index;

	fprintf(fp, "1 0\n");
	fprintf(fp, "1 1\n");
	fprintf(fp, "material_rho,\n");

	counter = 0;
	for(ii1=0; ii1<elem.size; ii1++){
		if(elem.array[ii1].label < 0) continue;
		index = search_material_prop_element(mat, phys, &(elem.array[ii1]));
		RC_NEG_CHK(index);
		for(ii2=0; ii2<elem.array[ii1].node_num; ii2++){
			counter ++;
			fprintf(fp, "%10d %15.7e\n", counter, mat.array[index].rho);
		}
	}

	return(NORMAL_RC);
}


static RC
avs_printbc_force (FILE *fp, BC_ARRAY bc_force, NODE_ARRAY node)
{
	int index;
	int ii1;

	fprintf(fp, "3 0\n");

	fprintf(fp, "3 1 1 1\n");
	fprintf(fp, "bc_force_x,\n");
	fprintf(fp, "bc_force_y,\n");
	fprintf(fp, "bc_force_z,\n");

	for(ii1=0; ii1<node.size; ii1++){
		if(node.array[ii1].label < 0) continue;
		index = search_bc_node(bc_force, node.array[ii1].label);
		if(index < 0){
			fprintf(fp, "%10d %15.7e %15.7e %15.7e\n", node.array[ii1].label,
			                                         0.0, 0.0, 0.0);
			continue;
		}
		fprintf(fp, "%10d ", bc_force.array[index].node);
		if(bc_force.array[index].v_type.x == BC_FORCE){
			fprintf(fp, "%15.7e ", bc_force.array[index].v.x);
		}else{
			fprintf(fp, "%15.7e ", 0.0);
		}
		if(bc_force.array[index].v_type.y == BC_FORCE){
			fprintf(fp, "%15.7e ", bc_force.array[index].v.y);
		}else{
			fprintf(fp, "%15.7e ", 0.0);
		}
		if(bc_force.array[index].v_type.z == BC_FORCE){
			fprintf(fp, "%15.7e", bc_force.array[index].v.z);
		}else{
			fprintf(fp, "%15.7e", 0.0);
		}
		fprintf(fp, "\n");
	}

	return(NORMAL_RC);
}


static RC
avs_printbc_temp (FILE *fp, BC_ARRAY temp, double def_temp, NODE_ARRAY node)
{
	int ii1, ii2;
	int index;
	double temp_val;

	fprintf(fp, "1 0\n");

	fprintf(fp, "1 1\n");
	fprintf(fp, "bc_temp,\n");

	for(ii1=0; ii1<node.size; ii1++){
		if(node.array[ii1].label < 0) continue;

		temp_val = def_temp;

		index = search_bc_node(temp, node.array[ii1].label);
		if(index >= 0){
			for(ii2=0; ii2<BC_SCALAR; ii2++){
				if(temp.array[index].s_type[ii2] == BC_TEMP){
					temp_val = temp.array[index].s[ii2];
					break;
				}
			}
		}

		fprintf(fp, "%10d %15.7e\n", node.array[ii1].label, temp_val);
	}

	return(NORMAL_RC);
}


static RC
avs_printsens (FILE *fp, ELEMENT_ARRAY elem, SENSITIVITY_ARRAY sens)
{
	int ii1, ii2;
	int counter;
	int index;

	fprintf(fp, "1 0\n");
	fprintf(fp, "1 1\n");
	fprintf(fp, "sensitivity,\n");

	counter = 0;
	for(ii1=0; ii1<elem.size; ii1++){
		if(elem.array[ii1].label < 0) continue;
		for(ii2=0; ii2<elem.array[ii1].node_num; ii2++){
			counter ++;
			RC_NEG_CHK( index = search_sensitivity_element_node( sens,
			                    elem.array[ii1].label,
			                    elem.array[ii1].node[ii2]) );
			fprintf(fp, "%10d %15.7e\n", counter, sens.array[index].s);
		}
	}
	return(NORMAL_RC);
}


static RC
avs_printstress (FILE *fp, ELEMENT_ARRAY elem, STRESS_ARRAY stress)
{
	int ii1, ii2;
	int counter;
	int index;
	double pr_val[3];

	fprintf(fp, "9 0\n");
	fprintf(fp, "9 1 1 1 1 1 1 1 1 1\n");
	fprintf(fp, "mises_stress,\n");
	fprintf(fp, "max_pr_stress,\n");
	fprintf(fp, "min_pr_stress,\n");
	fprintf(fp, "stress_x,\n");
	fprintf(fp, "stress_y,\n");
	fprintf(fp, "stress_z,\n");
	fprintf(fp, "stress_yz,\n");
	fprintf(fp, "stress_zx,\n");
	fprintf(fp, "stress_xy,\n");

	counter = 0;
	for(ii1=0; ii1<elem.size; ii1++){
		if(elem.array[ii1].label < 0) continue;
		for(ii2=0; ii2<elem.array[ii1].node_num; ii2++){
			counter ++;
			index = search_stress_element_node(stress,
			                                       elem.array[ii1].label,
			                                       elem.array[ii1].node[ii2]);
			RC_NEG_CHK(index);
			RC_TRY( principal3d(stress.array[index].v, pr_val, NULL, NULL,
			                    NULL) );

			fprintf(fp, "%10d %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e"
			            "%15.7e %15.7e\n",
			             counter,
			             mises_stress(stress.array[index]),
			             MAX2(MAX2(pr_val[0], pr_val[1]), pr_val[2]),
			             MIN2(MIN2(pr_val[0], pr_val[1]), pr_val[2]),
			             stress.array[index].v.x,
			             stress.array[index].v.y,
			             stress.array[index].v.z,
			             stress.array[index].v.yz,
			             stress.array[index].v.zx,
			             stress.array[index].v.xy);
		}
	}
	return(NORMAL_RC);
}


static int
check_step (FILE *fp_avs)
{
	IDC idc[1];

	fflush(fp_avs);

	if( fseek(fp_avs, 0L, SEEK_SET) < 0) return(-1);

	/* 現在までのステップ数を調べる */
	if(fgetidc(fp_avs, "i", idc ) != NORMAL_RC) return(-1);

	if( fseek(fp_avs, 0L, SEEK_SET) < 0) return(-1);

	/* 現在のステップ数 = 過去のステップ数 + 1 */
	(idc[0].i)++;
	fprintf(fp_avs, "%d", idc[0].i);
	fflush(fp_avs);
	if( fseek(fp_avs, 0L, SEEK_END) < 0) return(-1);

	return(idc[0].i);
}



RC
output_avs_disp (FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY elem, DISP_ARRAY disp)
{
	fprintf(fp, "1\n"); 
	fprintf(fp, "data_geom\n");
	fprintf(fp, "step1\n");
	RC_TRY( avs_printfe_model(fp, node, elem) );
	RC_TRY( avs_printdisp(fp, disp) );

	fflush(fp);
	
	return(NORMAL_RC);
}


RC
output_avs_velo (FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY elem,
                 VELO_ARRAY velo)
{
	fprintf(fp, "1\n"); 
	fprintf(fp, "data_geom\n");
	fprintf(fp, "step1\n");
	RC_TRY( avs_printfe_model(fp, node, elem) );
	RC_TRY( avs_printvelo(fp, velo) );

	fflush(fp);
	
	return(NORMAL_RC);
}


RC
output_avs_elem_density (FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY elem,
                         VELO_ARRAY density)
{
	fprintf(fp, "1\n"); 
	fprintf(fp, "data_geom\n");
	fprintf(fp, "step1\n");
	RC_TRY( avs_printfe_model(fp, node, elem) );
	RC_TRY( avs_printelem_density(fp, density) );

	fflush(fp);
	
	return(NORMAL_RC);
}


RC
output_avs_material_rho (FILE *fp, NODE_ARRAY node,
                         ELEMENT_ARRAY elem, MATERIAL_PROP_ARRAY mat,
                         PHYSICAL_PROP_ARRAY phys)
{
	if(elem.size <= 0) return(WRITE_ERROR_RC);

	fprintf(fp,"1\n");
	fprintf(fp,"data\n");
	fprintf(fp,"step1\n");

	RC_TRY( avs_printfe_model4discont_data(fp, node, elem) );
	RC_TRY( avs_printmaterial_rho(fp, elem, mat, phys) );

	fflush(fp);

	return(NORMAL_RC);
}


/* 節点力境界条件表示 */
RC
output_avs_bc_force (FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY elem,
                     BC_ARRAY bc)
{
	fprintf(fp, "1\n"); 
	fprintf(fp, "data_geom\n");
	fprintf(fp, "step1\n");
	RC_TRY( avs_printfe_model(fp, node, elem) );
	RC_TRY( avs_printbc_force(fp, bc, node) );

	fflush(fp);

	return(NORMAL_RC);
}


/* 節点温度境界条件表示 */
RC
output_avs_bc_temp (FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY elem,
                    BC_ARRAY temp, double def_temp)
{
	fprintf(fp, "1\n"); 
	fprintf(fp, "data_geom\n");
	fprintf(fp, "step1\n");
	RC_TRY( avs_printfe_model(fp, node, elem) );
	RC_TRY( avs_printbc_temp(fp, temp, def_temp, node) );

	fflush(fp);

	return(NORMAL_RC);
}


RC
output_avs_sens (FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY elem,
                 SENSITIVITY_ARRAY sens)
{
	if(elem.size <= 0) return(WRITE_ERROR_RC);

	fprintf(fp,"1\n");
	fprintf(fp,"data\n");
	fprintf(fp,"step1\n");

	RC_TRY( avs_printfe_model4discont_data(fp, node, elem) );
	RC_TRY( avs_printsens(fp, elem, sens) );

	return(NORMAL_RC);
}


RC
output_avs_stress (FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY elem,
                   STRESS_ARRAY stress)
{
	if(elem.size <= 0) return(WRITE_ERROR_RC);

	fprintf(fp,"1\n");
	fprintf(fp,"data\n");
	fprintf(fp,"step1\n");

	RC_TRY( avs_printfe_model4discont_data(fp, node, elem) );
	RC_TRY( avs_printstress(fp, elem, stress) );

	fflush(fp);

	return(NORMAL_RC);
}



/* 各ステップごとにデータのみ変化する場合の fopen */
/* FileName = "***.inp"                           */
FILE *
fopen_avs_step_data (const char *FileName, NODE_ARRAY fix_node,
                     ELEMENT_ARRAY fix_elem)
{
	FILE *fp;
	RC rc;
	int ii1;

	fp = fopen(FileName, "w+");
	if(fp == NULL) return(NULL);

	fprintf(fp, "0");
	for(ii1=0; ii1<INT_PLACE_NUMBER; ii1++){
		fprintf(fp, " ");
	}
	fprintf(fp, "\n");

	fprintf(fp, "data\n");
	fprintf(fp, "step1\n");

	rc = avs_printfe_model(fp, fix_node, fix_elem);
	if(rc != NORMAL_RC) return(NULL);
	
	return(fp);
}


/* FILE *fopen_avs_step_data() 専用 */
RC
output_avs_step_disp (FILE *fp_avs, DISP_ARRAY disp)
{
	int step_num;

	RC_NEG_ZERO_CHK( step_num = check_step(fp_avs) );
	if(step_num != 1) fprintf(fp_avs, "step%d\n", step_num);

	RC_TRY( avs_printdisp(fp_avs, disp) );

	fflush(fp_avs);

	return(NORMAL_RC);
}


/* FILE *fopen_avs_step_data() 専用 */
RC
output_avs_step_velo (FILE *fp_avs, VELO_ARRAY velo)
{
	int step_num;

	RC_NEG_ZERO_CHK( step_num = check_step(fp_avs) );
	if(step_num != 1) fprintf(fp_avs, "step%d\n", step_num);

	RC_TRY( avs_printvelo(fp_avs, velo) );

	fflush(fp_avs);

	return(NORMAL_RC);
}


/* FILE *fopen_avs_step_data() 専用 */
RC
output_avs_step_elem_density (FILE *fp_avs, VELO_ARRAY density)
{
	int step_num;

	RC_NEG_ZERO_CHK( step_num = check_step(fp_avs) );
	if(step_num != 1) fprintf(fp_avs, "step%d\n", step_num);

	RC_TRY( avs_printelem_density(fp_avs, density) );

	fflush(fp_avs);

	return(NORMAL_RC);
}


/* FILE *fopen_avs_step_data() 専用 */
RC
output_avs_step_bc_force (FILE *fp_avs, NODE_ARRAY node, BC_ARRAY bc)
{
	int step_num;

	RC_NEG_ZERO_CHK( step_num = check_step(fp_avs) );
	if(step_num != 1) fprintf(fp_avs, "step%d\n", step_num);

	RC_TRY( avs_printbc_force(fp_avs, bc, node) );

	fflush(fp_avs);

	return(NORMAL_RC);
}


/* 各ステップごとに形状、データ共に変化する場合の fopen */
/* FileName = "***.inp"                                 */
FILE *
fopen_avs_step_datageom (const char *FileName)
{
	FILE *fp;
	int ii1;

	fp = fopen(FileName, "w+");
	if(fp == NULL) return(NULL);

	fprintf(fp, "0");
	for(ii1=0; ii1<INT_PLACE_NUMBER; ii1++){
		fprintf(fp, " ");
	}
	fprintf(fp, "\n");

	fprintf(fp, "data_geom\n");

	return(fp);
}


/* FILE *fopen_avs_step_datageom() 用 */
RC
output_avs_step_disp_model (FILE *fp_avs, NODE_ARRAY node,
                            ELEMENT_ARRAY elem, DISP_ARRAY disp)
{
	int step_num;

	RC_NEG_ZERO_CHK( step_num = check_step(fp_avs) );

	fprintf(fp_avs, "step%d\n", step_num);
	RC_TRY( avs_printfe_model(fp_avs, node, elem) );
	RC_TRY( avs_printdisp(fp_avs, disp) );

	fflush(fp_avs);

	return(NORMAL_RC);
}


/* FILE *fopen_avs_step_datageom() 用 */
RC
output_avs_step_velo_model (FILE *fp_avs, NODE_ARRAY node,
                            ELEMENT_ARRAY elem, VELO_ARRAY velo)
{
	int step_num;

	RC_NEG_ZERO_CHK( step_num = check_step(fp_avs) );

	fprintf(fp_avs, "step%d\n", step_num);
	RC_TRY( avs_printfe_model(fp_avs, node, elem) );
	RC_TRY( avs_printvelo(fp_avs, velo) );

	fflush(fp_avs);

	return(NORMAL_RC);
}


/* FILE *fopen_avs_step_datageom() 用 */
RC
output_avs_step_sens_model (FILE *fp_avs, NODE_ARRAY node,
                            ELEMENT_ARRAY elem,
                            SENSITIVITY_ARRAY sens)
{
	int step_num;

	RC_NEG_ZERO_CHK( step_num = check_step(fp_avs) );

	fprintf(fp_avs, "step%d\n", step_num);
	RC_TRY( avs_printfe_model4discont_data(fp_avs, node, elem) );
	RC_TRY( avs_printsens(fp_avs, elem, sens) );

	fflush(fp_avs);

	return(NORMAL_RC);
}


/* FILE *fopen_avs_step_datageom() 用 */
RC
output_avs_step_stress_model (FILE *fp_avs, NODE_ARRAY node,
                              ELEMENT_ARRAY elem,
                              STRESS_ARRAY stress)
{
	int step_num;

	RC_NEG_ZERO_CHK( step_num = check_step(fp_avs) );

	fprintf(fp_avs, "step%d\n", step_num);
	RC_TRY( avs_printfe_model4discont_data(fp_avs, node, elem) );
	RC_TRY( avs_printstress(fp_avs, elem, stress) );

	fflush(fp_avs);

	return(NORMAL_RC);
}


/* FILE *fopen_avs_step_datageom() 用 */
RC
output_avs_step_material_rho_model (FILE *fp_avs, NODE_ARRAY node,
                                    ELEMENT_ARRAY elem,
                                    MATERIAL_PROP_ARRAY mat,
                                    PHYSICAL_PROP_ARRAY phys)
{
	int step_num;

	RC_NEG_ZERO_CHK( step_num = check_step(fp_avs) );

	fprintf(fp_avs, "step%d\n", step_num);
	RC_TRY( avs_printfe_model4discont_data(fp_avs, node, elem) );
	RC_TRY( avs_printmaterial_rho(fp_avs, elem, mat, phys) );

	fflush(fp_avs);

	return(NORMAL_RC);
}


/* FILE *fopen_avs_step_datageom() 用 */
RC
output_avs_step_bc_force_model (FILE *fp_avs, NODE_ARRAY node,
                                ELEMENT_ARRAY elem, BC_ARRAY bc)
{
	int step_num;

	RC_NEG_ZERO_CHK( step_num = check_step(fp_avs) );

	fprintf(fp_avs, "step%d\n", step_num);
	RC_TRY( avs_printfe_model(fp_avs, node, elem) );
	RC_TRY( avs_printbc_force(fp_avs, bc, node) );

	fflush(fp_avs);

	return(NORMAL_RC);
}


/* output_avs_fe_model_****() の上位関数 */
/* elem の次元数を見て処理 */
RC
output_avs_fe_model (FILE *fp_mgf, NODE_ARRAY node, ELEMENT_ARRAY elem)
{
	ELEMENT_ARRAY surf;

	switch( analysis_dim(elem) ){
	case 3:
		RC_TRY( extract_surface(elem, &surf) );
		RC_TRY( output_avs_fe_model_face(fp_mgf, node, surf) );
		RC_TRY( free_element_array(&surf) );
		break;
	case 2:
		RC_TRY( output_avs_fe_model_face(fp_mgf, node, elem) );
		break;
	case 1:
		RC_TRY( output_avs_fe_model_line(fp_mgf, node, elem) );
		break;
	default:
		return(WRITE_ERROR_RC);
	}

	return(NORMAL_RC);
}


/*
 * 2次元要素表示専用
 * 作成するファイルの拡張子は .mgf
 * 3次元要素の場合、elem は事前に extract_surface() などに通す必要有り
 * 二次要素による曲線は、直線として表示される
 */
static RC
output_avs_fe_model_face (FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY elem)
{
	int node_label[8];
	int index;
	int ii1, ii2;

	if(analysis_dim(elem) != 2) return(WRITE_ERROR_RC);

	fprintf(fp, "# Micro AVS Geom:2.00\n");
	fprintf(fp, "polyhedron\n");
	fprintf(fp, "fe_model_face\n");
	fprintf(fp, "facet\n");
	fprintf(fp, "vertex\n");
	
	fprintf(fp, "%d\n", node.size);
	
	for(ii1=0; ii1<node.size; ii1++){
		if(node.array[ii1].label < 0){
			fprintf(fp, "%15.7E %15.7E %15.7E\n", 0.0, 0.0, 0.0);
		}else{
			fprintf(fp, "%15.7E %15.7E %15.7E\n", node.array[ii1].p.x,
			                                      node.array[ii1].p.y,
			                                      node.array[ii1].p.z);
		}
	}

	fprintf( fp, "%d\n", count_valid_element(elem) );

	for(ii1=0; ii1<elem.size; ii1++){
		if(elem.array[ii1].label < 0) continue;
		switch(elem.array[ii1].type){
		case ELEM_TRI2:
			node_label[0] = elem.array[ii1].node[0];
			node_label[1] = elem.array[ii1].node[5];
			node_label[2] = elem.array[ii1].node[1];
			node_label[3] = elem.array[ii1].node[3];
			node_label[4] = elem.array[ii1].node[2];
			node_label[5] = elem.array[ii1].node[4];
			break;
		case ELEM_QUAD2:
			node_label[0] = elem.array[ii1].node[0];
			node_label[1] = elem.array[ii1].node[4];
			node_label[2] = elem.array[ii1].node[1];
			node_label[3] = elem.array[ii1].node[5];
			node_label[4] = elem.array[ii1].node[2];
			node_label[5] = elem.array[ii1].node[6];
			node_label[6] = elem.array[ii1].node[3];
			node_label[7] = elem.array[ii1].node[7];
			break;
		case ELEM_TRI1:
		case ELEM_QUAD1:
			for(ii2=0; ii2<elem.array[ii1].node_num; ii2++){
				node_label[ii2] = elem.array[ii1].node[ii2];
			}
			break;
		default:
			return(ARG_ERROR_RC);
		}
		fprintf(fp, "%d\n", elem.array[ii1].node_num);
		for(ii2=0; ii2<elem.array[ii1].node_num; ii2++){
			RC_NEG_CHK( index = search_node_label(node, node_label[ii2]) );
			fprintf(fp, "%10d ", index + 1);
		}
		fprintf(fp,"\n");
	}

	return(NORMAL_RC);
}


/* １次元要素表示専用                                                   */
/* 作成するファイルの拡張子は .mgf                                      */
/* ２次元要素の場合、elem は事前に extract_surface() などに通す必要有り */
static RC
output_avs_fe_model_line (FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY elem)
{
	int node_index[3];
	int ii1, ii2;

	if(analysis_dim(elem) != 1) return(WRITE_ERROR_RC);
	
	fprintf(fp, "# Micro AVS Geom:2.00\n");
	fprintf(fp, "disjoint line\n");
	fprintf(fp, "fe_model_line\n");
	fprintf(fp, "vertex\n");

	switch( analysis_order(elem) ){
	case 2:
		fprintf(fp, "%d\n", count_valid_element(elem) * 4);
		break;
	case 1:
		fprintf(fp, "%d\n", count_valid_element(elem) * 2);
		break;
	default:
		return(WRITE_ERROR_RC);
	}

	for(ii1=0; ii1<elem.size; ii1++){
		if(elem.array[ii1].label < 0) continue;
		for(ii2=0; ii2<elem.array[ii1].node_num; ii2++){
			node_index[ii2] = search_node_label(node,
			                                        elem.array[ii1].node[ii2]);
			RC_NEG_CHK(node_index[ii2]);
		}

		switch(elem.array[ii1].type){
		case ELEM_BEAM2:
		case ELEM_LINE2:
			fprintf(fp, "%15.7e %15.7e %15.7e\n%15.7e %15.7e %15.7e\n\n",
			            node.array[node_index[0]].p.x,
			            node.array[node_index[0]].p.y,
			            node.array[node_index[0]].p.z,
			            node.array[node_index[2]].p.x,
			            node.array[node_index[2]].p.y,
			            node.array[node_index[2]].p.z);
			fprintf(fp, "%15.7e %15.7e %15.7e\n%15.7e %15.7e %15.7e\n\n",
			            node.array[node_index[2]].p.x,
			            node.array[node_index[2]].p.y,
			            node.array[node_index[2]].p.z,
			            node.array[node_index[1]].p.x,
			            node.array[node_index[1]].p.y,
			            node.array[node_index[1]].p.z);
			break;
		case ELEM_BEAM1:
		case ELEM_LINE1:
			for(ii2=0; ii2<2; ii2++){
				fprintf(fp, "%15.7e %15.7e %15.7e\n",
				            node.array[node_index[ii2]].p.x,
				            node.array[node_index[ii2]].p.y,
				            node.array[node_index[ii2]].p.z);
			}
			fprintf(fp,"\n");
			break;
		default:
			return(WRITE_ERROR_RC);
		}

	}

	return(NORMAL_RC);
}


