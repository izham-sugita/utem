/*********************************************************************
 * optimization_utl.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Thanks to following contributors.
 *  <Yasuo ASAGA> <Tomoyuki TSUBATA> <Takaaki NAGATANI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: optimization_utl.c,v 1.27 2004/01/08 11:31:07 nagatani Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "rc.h"
#include "math_utl.h"
#include "fem_struct.h"
#include "optimization_utl.h"
#include "fem_solver.h"
#include "sky_cholesky.h"

#define EXP_LIMIT        (0.99*log(DBL_MAX/2.0))
#define MIN_RHO          (10.0*DBL_MIN)


static RC make_matrix(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                      FEM_MATERIAL_PROP_ARRAY material,
                      FEM_PHYSICAL_PROP_ARRAY physical,
                      NONZERO_MATRIX3 *nonzero, FEM_SKYLINE_MATRIX *skyline);
static RC restdisp_work(FEM_BC_ARRAY rest, FEM_REACT_ARRAY react,
                        double *compliance);


/* solver_flag = 0: スカイライン             */
/*             = 1: ICCG                     */
RC speed_analysis(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                  FEM_MATERIAL_PROP_ARRAY material,
                  FEM_PHYSICAL_PROP_ARRAY physical,
                  FEM_BC_ARRAY shape_rest, int num_case,
                  const FEM_BC_ARRAY tractions[], FEM_DISP_ARRAY disps[],
                  int solver_flag)
{
	NONZERO_MATRIX3 nonzero;
	FEM_SKYLINE_MATRIX skyline;
	double **f_vect;
	int ii1;
	int dim;

	RC_NULL_CHK( f_vect = (double **)malloc(num_case * sizeof(double *)) );

	/* 剛性マトリックス作成 */
	if(solver_flag == 0){
		RC_TRY( make_matrix(node, element, material, physical,
		                    NULL, &skyline) );
	}else if(solver_flag == 1){
		RC_TRY( make_matrix(node, element, material, physical,
		                    &nonzero, NULL) );
	}else{
		return(IMPLEMENT_ERROR_RC);
	}

	/* 荷重ベクトル作成 */
	dim = fem_analysis_dim(element);
	for(ii1=0; ii1<num_case; ii1++){
		RC_TRY( fem_allocate_vector(dim, node, &(f_vect[ii1])) );
		RC_TRY( fem_add_nodal_force(dim, node, tractions[ii1], f_vect[ii1]) );
	}

	/* 拘束条件組み込み */
	if(solver_flag == 0){
		RC_TRY( fem_modify_g_matrix_n(node, shape_rest,
		                              skyline, f_vect, num_case) );
	}else{         /* solver_flag == 1 */
		RC_TRY( fem_modify_nonzero3_n(node, shape_rest,
		                              nonzero, f_vect, num_case) );
	}

	/* スカイラインの場合，LU 分解 */
	if(solver_flag == 0) RC_TRY( fem_decomp_g_matrix(skyline) );

	/* 変位の計算 */
	for(ii1=0; ii1<num_case; ii1++){
		if(solver_flag == 0){
			RC_TRY( fem_subst_g_matrix(node, skyline,
			                           f_vect[ii1], &(disps[ii1])) );
		}else{         /* solver_flag == 1 */
			double *d_vect;

			RC_TRY( fem_iccg_nonzero3(nonzero, f_vect[ii1], &d_vect, NULL) );
			RC_TRY( fem_vect2disp(3, node, d_vect, &(disps[ii1])) );
			RC_TRY( fem_free_vector(&d_vect) );
		}
		RC_TRY( fem_free_vector(&(f_vect[ii1])) );
	}

	/* マトリックスの開放 */
	if(solver_flag == 0){
		RC_TRY( fem_free_g_matrix(&skyline) );
	}else{         /* solver_flag == 1 */
		RC_TRY( free_nonzero_matrix3(&nonzero) );
	}
	free((void *)f_vect);

	return(NORMAL_RC);
}


static RC make_matrix(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                      FEM_MATERIAL_PROP_ARRAY material,
                      FEM_PHYSICAL_PROP_ARRAY physical,
                      NONZERO_MATRIX3 *nonzero, FEM_SKYLINE_MATRIX *skyline)
{
	int ii1;
	FEM_ELEM_MATRIX elem_matrix;

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	/* 剛性マトリックスの作成と剛性マトリックス対角項合計値の計算 */
	if(nonzero != NULL){
		RC_TRY( fem_allocate_nonzero3(node, element, nonzero) );
		for(ii1=0; ii1<element.size; ii1++){
			if(element.array[ii1].label < 0) continue;

			RC_TRY( fem_make_elem_matrix(element.array[ii1], node,
			                             material, physical, &elem_matrix) );
			RC_TRY( fem_impose_nonzero3(*nonzero, elem_matrix, 1.0) );
			RC_TRY( free_fem_elem_matrix(&elem_matrix) );
		}
	}
	if(skyline != NULL){
		RC_TRY( fem_make_g_matrix(node, element,
		                          material, physical, skyline) );
	}

	return(NORMAL_RC);
}


/* 材料番号 NONDESIGN_MAT_LABEL 以上の要素に属する節点に対して */
/* 完全拘束条件を生成し，shape_rest に加える */
RC add_shape_rest(FEM_ELEMENT_ARRAY element, FEM_MATERIAL_PROP_ARRAY material,
                  FEM_PHYSICAL_PROP_ARRAY physical, FEM_BC_ARRAY *shape_rest)
{
	int mat_index;
	int ii1, ii2;

	for(ii1=0; ii1<element.size; ii1++){
		RC_NEG_CHK( mat_index = search_fem_material_prop_element(material,
		                           physical, &(element.array[ii1])) );
		if(material.array[mat_index].label < NONDESIGN_MAT_LABEL) continue;
		for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
			RC_TRY( realloc_fem_bc_array(shape_rest) );
			shape_rest->array[shape_rest->size].node
			                 = element.array[ii1].node[ii2];
			shape_rest->array[shape_rest->size].v_type.x = BC_FIX;
			shape_rest->array[shape_rest->size].v_type.y = BC_FIX;
			shape_rest->array[shape_rest->size].v_type.z = BC_FIX;
			shape_rest->array[shape_rest->size].v_type.yz = BC_FIX;
			shape_rest->array[shape_rest->size].v_type.zx = BC_FIX;
			shape_rest->array[shape_rest->size].v_type.xy = BC_FIX;
			shape_rest->array[shape_rest->size].v.x = 0.0;
			shape_rest->array[shape_rest->size].v.y = 0.0;
			shape_rest->array[shape_rest->size].v.z = 0.0;
			shape_rest->array[shape_rest->size].v.yz = 0.0;
			shape_rest->array[shape_rest->size].v.zx = 0.0;
			shape_rest->array[shape_rest->size].v.xy = 0.0;
			(shape_rest->size)++;
		}
	}
	RC_TRY( clean_fem_bc_array(shape_rest) );
	RC_TRY( compress_fem_bc_array(shape_rest) );

	return(NORMAL_RC);
}


/* 固有値に関する感度 */
/* mode は，質量マトリックスに対して正規化しておく */
RC shape_trac_eigen_value(FEM_ELEMENT_ARRAY surf, FEM_DISP_ARRAY mode,
                          FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                          FEM_MATERIAL_PROP_ARRAY material,
                          FEM_PHYSICAL_PROP_ARRAY physical,
                          double lambda, FEM_BC_ARRAY *sens_trac)
{
	int ii1, ii2, ii3;
	FEM_WORK_SET ws;

	RC_TRY( allocate_fem_work_set(&ws) );

	RC_TRY( allocate_fem_bc_array(node.size, sens_trac) );
	for(ii1=0; ii1<node.size; ii1++){
		sens_trac->array[ii1].node = node.array[ii1].label;
	}
	sens_trac->size = node.size;
	RC_TRY( clean_fem_bc_array(sens_trac) );

	for(ii1=0; ii1<surf.size; ii1++){
		double *weights = ws.vec_G[0];
		double **gauss_points = ws.mat_G_4[0];
		double **v_gauss_points = ws.mat_G_4[1];
		double **node_location = ws.mat_N_3[0];
		double **disp_matrix = ws.mat_N_6[0];
		int num_points;
		TRANSLATION3D f[FEM_MAX_NODE];
		int elem_index, mat_index;

		if(surf.array[ii1].label < 0) continue;

		RC_NEG_CHK( elem_index = search_fem_element_label(element,
		                                      surf.array[ii1].i_info[0]) );
		RC_NEG_CHK( mat_index = search_fem_material_prop_element(material,
		                           physical, &(element.array[elem_index])) );
		if(material.array[mat_index].label >= IGNORE_MAT_LABEL) continue;

		RC_TRY( Gauss_const_max(surf.array[ii1].type,
		                        gauss_points, weights, &num_points) );
		RC_TRY( trans_coord_s2v(surf.array[ii1], element.array[elem_index],
		                        num_points, gauss_points, v_gauss_points) );

		RC_TRY( node_location_matrix(&(surf.array[ii1]), node, node_location) );
		RC_TRY( fill_disp_matrix(element.array[elem_index],
		                         mode, disp_matrix) );

		for(ii2=0; ii2<FEM_MAX_NODE; ii2++){
			f[ii2].x = f[ii2].y = f[ii2].z = 0.0;
		}

		for(ii2=0; ii2<num_points; ii2++){
			double *N = ws.vec_N[0];
			TRANSLATION3D d_A;
			FEM_STRESS_STRAIN local_stress, local_strain;
			double g1, g2, G;
			TRANS_ROTATION3D disp_v;

			/* 積分の準備 */
			RC_TRY( set_N_vector(surf.array[ii1].type, gauss_points[ii2], N) );
			RC_TRY( interpolate_disp(element.array[elem_index],
			                    v_gauss_points[ii2], disp_matrix, &disp_v) );
			RC_TRY( cal_d_A(surf.array[ii1], gauss_points[ii2], node_location,
			                &d_A, ws.mat_3_N[0], ws.mat_3_3[0]) );

			/* 感度 G の計算 */
			RC_TRY( fem_local_stress_strain(element.array[elem_index],
			           v_gauss_points[ii2], node, mode, material, physical,
			           &local_stress, &local_strain) );
			g1 = (local_stress.v.x * local_strain.v.x
			    + local_stress.v.y * local_strain.v.y
			    + local_stress.v.z * local_strain.v.z
			    + local_stress.v.yz * local_strain.v.yz
			    + local_stress.v.zx * local_strain.v.zx
			    + local_stress.v.xy * local_strain.v.xy);
			g2 = - lambda * material.array[mat_index].rho
			     *( (disp_v.x*disp_v.x) + (disp_v.y*disp_v.y)
			                            + (disp_v.z*disp_v.z) );
			G = g1 + g2;

			for(ii3=0; ii3<surf.array[ii1].node_num; ii3++){
				f[ii3].x += weights[ii2] * G * N[ii3] * d_A.x;
				f[ii3].y += weights[ii2] * G * N[ii3] * d_A.y;
				f[ii3].z += weights[ii2] * G * N[ii3] * d_A.z;
			}
		}

		/* f -> sens_trac */
		for(ii2=0; ii2<surf.array[ii1].node_num; ii2++){
			int index = search_fem_bc_node(*sens_trac,
			                               surf.array[ii1].node[ii2]);
			RC_NEG_CHK(index);
			sens_trac->array[index].v_type.x = BC_FORCE;
			sens_trac->array[index].v_type.y = BC_FORCE;
			sens_trac->array[index].v_type.z = BC_FORCE;
			sens_trac->array[index].v.x += f[ii2].x;
			sens_trac->array[index].v.y += f[ii2].y;
			sens_trac->array[index].v.z += f[ii2].z;
		}
	}

	RC_TRY( compress_fem_bc_array(sens_trac) );

	RC_TRY( free_fem_work_set(&ws) );

	return(NORMAL_RC);
}


RC shape_trac_max_mises(FEM_ELEMENT_ARRAY surf, FEM_DISP_ARRAY disp,
                        FEM_DISP_ARRAY adj_disp,
                        FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                        FEM_MATERIAL_PROP_ARRAY material,
                        FEM_PHYSICAL_PROP_ARRAY physical,
                        double rho, FEM_BC_ARRAY *sens_trac, int reduce_flag)
{
	int ii1, ii2, ii3;
	FEM_WORK_SET ws;
	double m_integ, volume;
	double max_val;

	RC_TRY( allocate_fem_work_set(&ws) );

	RC_TRY( allocate_fem_bc_array(node.size, sens_trac) );
	for(ii1=0; ii1<node.size; ii1++){
		sens_trac->array[ii1].node = node.array[ii1].label;
	}
	sens_trac->size = node.size;
	RC_TRY( clean_fem_bc_array(sens_trac) );

	RC_TRY( mises_integration(node, element, material, physical,
	                                        disp, rho, &m_integ, &volume) );

	RC_TRY( max_mises_stress(node, element, material,
	                         physical, disp, &max_val) );

	for(ii1=0; ii1<surf.size; ii1++){
		double *weights = ws.vec_G[0];
		double **gauss_points = ws.mat_G_4[0];
		double **v_gauss_points = ws.mat_G_4[1];
		double **node_location = ws.mat_N_3[0];
		int num_points;
		TRANSLATION3D f[FEM_MAX_NODE];
		int elem_index, mat_index;
		FEM_ELEMENT surf_elem;

		if(surf.array[ii1].label < 0) continue;
		surf_elem = surf.array[ii1];
		if( reduce_flag && (element_order(surf.array[ii1].type) == 2) ){
			/* 中間節点に荷重を振り分けない */
			RC_TRY( reduced_element_type(surf_elem.type, &(surf_elem.type),
			                             &(surf_elem.node_num)) );
		}

		RC_NEG_CHK( elem_index = search_fem_element_label(element,
		                                      surf_elem.i_info[0]) );
		RC_NEG_CHK( mat_index = search_fem_material_prop_element(material,
		                           physical, &(element.array[elem_index])) );
		if(material.array[mat_index].label >= IGNORE_MAT_LABEL) continue;

		RC_TRY( Gauss_const_max(surf_elem.type,
		                        gauss_points, weights, &num_points) );
		RC_TRY( trans_coord_s2v(surf_elem, element.array[elem_index],
		                        num_points, gauss_points, v_gauss_points) );

		RC_TRY( node_location_matrix(&(surf_elem), node, node_location) );

		for(ii2=0; ii2<FEM_MAX_NODE; ii2++){
			f[ii2].x = f[ii2].y = f[ii2].z = 0.0;
		}

		for(ii2=0; ii2<num_points; ii2++){
			double *N = ws.vec_N[0];
			TRANSLATION3D d_A;
			FEM_STRESS_STRAIN local_stress, local_strain;
			FEM_STRESS_STRAIN adj_local_stress, adj_local_strain;
			double mises_str, g1, g2, g3, G;

			/* 積分の準備 */
			RC_TRY( set_N_vector(surf_elem.type, gauss_points[ii2], N) );
			RC_TRY( cal_d_A(surf_elem, gauss_points[ii2], node_location,
			                &d_A, ws.mat_3_N[0], ws.mat_3_3[0]) );

			/* 感度 G の計算 */
			RC_TRY( fem_local_stress_strain(element.array[elem_index],
			           v_gauss_points[ii2], node, disp, material, physical,
			           &local_stress, &local_strain) );
/*fprintf(stdout, "%3d:%3d:%15.7e %15.7e %15.7e\n"
                "%15.7e %15.7e %15.7e\n"
                "%15.7e\n", element.array[elem_index].label,
                ii2, local_stress.v.x, local_stress.v.y, local_stress.v.z,
                local_stress.v.yz, local_stress.v.zx, local_stress.v.xy,
                mises_stress(local_stress));*/
			RC_TRY( fem_local_stress_strain(element.array[elem_index],
			           v_gauss_points[ii2], node, adj_disp, material, physical,
			           &adj_local_stress, &adj_local_strain) );
			mises_str = mises_stress(local_stress);
			if(mises_str > 0.99*max_val) mises_str = 0.99*max_val;
			mises_str *= 0.99;
			g1 = - (local_stress.v.x * adj_local_strain.v.x
			      + local_stress.v.y * adj_local_strain.v.y
			      + local_stress.v.z * adj_local_strain.v.z
			      + local_stress.v.yz * adj_local_strain.v.yz
			      + local_stress.v.zx * adj_local_strain.v.zx
			      + local_stress.v.xy * adj_local_strain.v.xy)*(1.0+REL_TOL);
			g2 = - 1.0/(rho*volume);
			g3 = exp(rho*mises_str)/(rho*m_integ);
			G = g1 + g2 + g3;

			for(ii3=0; ii3<surf_elem.node_num; ii3++){
				f[ii3].x += weights[ii2] * G * N[ii3] * d_A.x;
				f[ii3].y += weights[ii2] * G * N[ii3] * d_A.y;
				f[ii3].z += weights[ii2] * G * N[ii3] * d_A.z;
			}
		}

		/* f -> sens_trac */
		for(ii2=0; ii2<surf_elem.node_num; ii2++){
			int index = search_fem_bc_node(*sens_trac, surf_elem.node[ii2]);
			RC_NEG_CHK(index);
			sens_trac->array[index].v_type.x = BC_FORCE;
			sens_trac->array[index].v_type.y = BC_FORCE;
			sens_trac->array[index].v_type.z = BC_FORCE;
			sens_trac->array[index].v.x += f[ii2].x;
			sens_trac->array[index].v.y += f[ii2].y;
			sens_trac->array[index].v.z += f[ii2].z;
		}
	}

	RC_TRY( compress_fem_bc_array(sens_trac) );

	RC_TRY( free_fem_work_set(&ws) );

	return(NORMAL_RC);
}


RC shape_trac_compliance(FEM_ELEMENT_ARRAY surf, FEM_DISP_ARRAY disp,
                         FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                         FEM_MATERIAL_PROP_ARRAY material,
                         FEM_PHYSICAL_PROP_ARRAY physical,
                         FEM_BC_ARRAY *sens_trac, int reduce_flag)
{
	int ii1, ii2, ii3;
	FEM_WORK_SET ws;

	RC_TRY( allocate_fem_work_set(&ws) );

	RC_TRY( allocate_fem_bc_array(node.size, sens_trac) );
	for(ii1=0; ii1<node.size; ii1++){
		sens_trac->array[ii1].node = node.array[ii1].label;
	}
	sens_trac->size = node.size;
	RC_TRY( clean_fem_bc_array(sens_trac) );

	for(ii1=0; ii1<surf.size; ii1++){
		double *weights = ws.vec_G[0];
		double **gauss_points = ws.mat_G_4[0];
		double **v_gauss_points = ws.mat_G_4[1];
		double **node_location = ws.mat_N_3[0];
		int num_points;
		TRANSLATION3D f[FEM_MAX_NODE];
		int elem_index;
		FEM_ELEMENT surf_elem;

		if(surf.array[ii1].label < 0) continue;
		surf_elem = surf.array[ii1];
		if( reduce_flag && (element_order(surf.array[ii1].type) == 2) ){
			/* 中間節点に荷重を振り分けない */
			RC_TRY( reduced_element_type(surf_elem.type, &(surf_elem.type),
			                             &(surf_elem.node_num)) );
		}

		RC_NEG_CHK( elem_index = search_fem_element_label(element,
		                                      surf_elem.i_info[0]) );

		RC_TRY( Gauss_const_max(surf_elem.type,
		                        gauss_points, weights, &num_points) );
		RC_TRY( trans_coord_s2v(surf_elem, element.array[elem_index],
		                        num_points, gauss_points, v_gauss_points) );

		RC_TRY( node_location_matrix(&(surf_elem), node, node_location) );

		for(ii2=0; ii2<FEM_MAX_NODE; ii2++){
			f[ii2].x = f[ii2].y = f[ii2].z = 0.0;
		}

		for(ii2=0; ii2<num_points; ii2++){
			double *N = ws.vec_N[0];
			TRANSLATION3D d_A;
			FEM_STRESS_STRAIN local_stress, local_strain;
			double G;

			/* 積分の準備 */
			RC_TRY( set_N_vector(surf_elem.type, gauss_points[ii2], N) );
			RC_TRY( cal_d_A(surf_elem, gauss_points[ii2], node_location,
			                &d_A, ws.mat_3_N[0], ws.mat_3_3[0]) );

			/* 感度 G の計算 */
			RC_TRY( fem_local_stress_strain(element.array[elem_index],
			           v_gauss_points[ii2], node, disp, material, physical,
			           &local_stress, &local_strain) );
			G = - (local_stress.v.x * local_strain.v.x
			     + local_stress.v.y * local_strain.v.y
			     + local_stress.v.z * local_strain.v.z
			     + local_stress.v.yz * local_strain.v.yz
			     + local_stress.v.zx * local_strain.v.zx
			     + local_stress.v.xy * local_strain.v.xy);

			for(ii3=0; ii3<surf_elem.node_num; ii3++){
				f[ii3].x += weights[ii2] * G * N[ii3] * d_A.x;
				f[ii3].y += weights[ii2] * G * N[ii3] * d_A.y;
				f[ii3].z += weights[ii2] * G * N[ii3] * d_A.z;
			}
		}

		/* f -> sens_trac */
		for(ii2=0; ii2<surf_elem.node_num; ii2++){
			int index = search_fem_bc_node(*sens_trac, surf_elem.node[ii2]);
			RC_NEG_CHK(index);
			sens_trac->array[index].v_type.x = BC_FORCE;
			sens_trac->array[index].v_type.y = BC_FORCE;
			sens_trac->array[index].v_type.z = BC_FORCE;
			sens_trac->array[index].v.x += f[ii2].x;
			sens_trac->array[index].v.y += f[ii2].y;
			sens_trac->array[index].v.z += f[ii2].z;
		}
	}

	RC_TRY( compress_fem_bc_array(sens_trac) );

	RC_TRY( free_fem_work_set(&ws) );

	return(NORMAL_RC);
}


RC shape_trac_volume(FEM_ELEMENT_ARRAY surf, FEM_NODE_ARRAY node,
                     FEM_BC_ARRAY *sens_trac, int reduce_flag)
{
	int ii1, ii2, ii3;
	FEM_WORK_SET ws;

	RC_TRY( allocate_fem_work_set(&ws) );

	RC_TRY( allocate_fem_bc_array(node.size, sens_trac) );
	for(ii1=0; ii1<node.size; ii1++){
		sens_trac->array[ii1].node = node.array[ii1].label;
	}
	sens_trac->size = node.size;
	RC_TRY( clean_fem_bc_array(sens_trac) );

	for(ii1=0; ii1<surf.size; ii1++){
		double *weights = ws.vec_G[0];
		double **gauss_points = ws.mat_G_4[0];
		double **node_location = ws.mat_N_3[0];
		int num_points;
		TRANSLATION3D f[FEM_MAX_NODE];
		FEM_ELEMENT surf_elem;

		if(surf.array[ii1].label < 0) continue;
		surf_elem = surf.array[ii1];
		if( reduce_flag && (element_order(surf.array[ii1].type) == 2) ){
			/* 中間節点に荷重を振り分けない */
			RC_TRY( reduced_element_type(surf_elem.type, &(surf_elem.type),
			                             &(surf_elem.node_num)) );
		}

		RC_TRY( Gauss_const_default(surf_elem.type,
		                            gauss_points, weights, &num_points) );
		RC_TRY( node_location_matrix(&(surf_elem), node, node_location) );

		for(ii2=0; ii2<FEM_MAX_NODE; ii2++){
			f[ii2].x = f[ii2].y = f[ii2].z = 0.0;
		}

		for(ii2=0; ii2<num_points; ii2++){
			double *N = ws.vec_N[0];
			TRANSLATION3D d_A;
			double G;

			/* 積分の準備 */
			RC_TRY( set_N_vector(surf_elem.type, gauss_points[ii2], N) );
			RC_TRY( cal_d_A(surf_elem, gauss_points[ii2], node_location,
			                &d_A, ws.mat_3_N[0], ws.mat_3_3[0]) );

			G = 1.0;
			for(ii3=0; ii3<surf_elem.node_num; ii3++){
				f[ii3].x += weights[ii2] * G * N[ii3] * d_A.x;
				f[ii3].y += weights[ii2] * G * N[ii3] * d_A.y;
				f[ii3].z += weights[ii2] * G * N[ii3] * d_A.z;
			}
		}

		/* f -> sens_trac */
		for(ii2=0; ii2<surf_elem.node_num; ii2++){
			int index = search_fem_bc_node(*sens_trac, surf_elem.node[ii2]);
			RC_NEG_CHK(index);
			sens_trac->array[index].v_type.x = BC_FORCE;
			sens_trac->array[index].v_type.y = BC_FORCE;
			sens_trac->array[index].v_type.z = BC_FORCE;
			sens_trac->array[index].v.x += f[ii2].x;
			sens_trac->array[index].v.y += f[ii2].y;
			sens_trac->array[index].v.z += f[ii2].z;
		}
	}

	RC_TRY( compress_fem_bc_array(sens_trac) );

	RC_TRY( free_fem_work_set(&ws) );

	return(NORMAL_RC);
}


RC shape_sens_max_mises(FEM_ELEMENT_ARRAY surf, FEM_DISP_ARRAY disp,
                        FEM_DISP_ARRAY adj_disp,
                        FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                        FEM_MATERIAL_PROP_ARRAY material,
                        FEM_PHYSICAL_PROP_ARRAY physical,
                        double rho, FEM_SENSITIVITY_ARRAY *sens)
{
	int ii1, ii2, ii3;
	static int init_flag = 0;
	static double **local_points;
	FEM_STRESS_STRAIN local_stress, local_strain;
	FEM_STRESS_STRAIN adj_local_stress, adj_local_strain;
	int elem_index, mat_index;
	double m_integ, volume, max_val;
	double mises_str, exp_limit = EXP_LIMIT;
	double g1, g2, g3;

	if(init_flag == 0){
		RC_TRY( allocate2D(FEM_MAX_NODE, 4, &local_points) );
		init_flag = 1;
	}

	RC_TRY( allocate_fem_sensitivity_array(0, sens) );
	sens->size = 0;
	sens->type = SENS_SHAPE_GRAD_DENSITY;

	RC_TRY( mises_integration(node, element, material, physical, disp, rho,
	                          &m_integ, &volume) );
	/*RC_TRY( volume_functional(node, element, &volume) );*/
	RC_TRY( max_mises_stress(node, element, material, physical, disp,
	                         &max_val) );

	for(ii1=0; ii1<surf.size; ii1++){
		if(surf.array[ii1].label < 0) continue;

		elem_index = search_fem_element_label(element,
		                                      surf.array[ii1].i_info[0]);
		if(elem_index < 0) return(SEEK_ERROR_RC);
		mat_index = search_fem_material_prop_element(material, physical,
		                                       &(element.array[elem_index]));
		RC_NEG_CHK( mat_index );

		RC_TRY( set_local_node_points(element.array[elem_index].type,
		                                                local_points) );
		for(ii2=0; ii2<surf.array[ii1].node_num; ii2++){
			RC_TRY( realloc_fem_sensitivity_array(sens) );
			sens->array[sens->size].element = surf.array[ii1].label;
			sens->array[sens->size].node = surf.array[ii1].node[ii2];
			if( material.array[mat_index].label >= IGNORE_MAT_LABEL ){
				sens->array[sens->size].s = 0.0;
				(sens->size)++;
				continue;
			}

			for(ii3=0; ii3<element.array[elem_index].node_num; ii3++){
				if(surf.array[ii1].node[ii2]
				   == element.array[elem_index].node[ii3]) break;
			}
			if(ii3 >= element.array[elem_index].node_num){
				return(SEEK_ERROR_RC);
			}

			RC_TRY( fem_local_stress_strain(element.array[elem_index],
			           local_points[ii3], node, disp, material, physical,
			           &local_stress, &local_strain) );
			RC_TRY( fem_local_stress_strain(element.array[elem_index],
			           local_points[ii3], node, adj_disp, material, physical,
			           &adj_local_stress, &adj_local_strain) );
			mises_str = rho * mises_stress(local_stress);
			if(mises_str > exp_limit) mises_str = exp_limit;
			g1 = - (local_stress.v.x * adj_local_strain.v.x
			      + local_stress.v.y * adj_local_strain.v.y
			      + local_stress.v.z * adj_local_strain.v.z
			      + local_stress.v.yz * adj_local_strain.v.yz
			      + local_stress.v.zx * adj_local_strain.v.zx
			      + local_stress.v.xy * adj_local_strain.v.xy);
			g2 = - 1.0/(rho*volume);
			g3 = exp(mises_str)/(rho*m_integ);
			sens->array[sens->size].s = g1 + g2 + g3;
			sens->size++;
		}
	}

	RC_TRY( clean_fem_sensitivity_array(sens) );

	return(NORMAL_RC);
}


RC max_mises_adj_force(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                       FEM_MATERIAL_PROP_ARRAY material,
                       FEM_PHYSICAL_PROP_ARRAY physical, FEM_DISP_ARRAY disp,
                       double rho, FEM_BC_ARRAY *adj_force)
{
	int ii1, ii2, ii3;
	double **node_location = NULL;
	double **gauss_point = NULL;
	double weight[MAX_GAUSS_POINT];
	double **D_matrix = NULL;
	double **B_matrix = NULL;
	double **trans_B = NULL;
	double **BtD = NULL;
	int num_point;
	int dim, ssnum;
	FEM_STRESS_STRAIN local_stress, local_strain;
	double det_J;
	double mises_str;
	double exp_limit = EXP_LIMIT;
	TRANS_ROTATION3D adj_strain;
	double thickness;
	double adj_strain_v[6];
	double force_v[3*FEM_MAX_NODE];
	double force_v_sum[3*FEM_MAX_NODE];
	double m_integ, volume;
	int mat_index;

	RC_TRY( allocate2D(MAX_GAUSS_POINT, 4, &gauss_point) );
	RC_TRY( allocate2D(FEM_MAX_NODE, 3, &node_location) );
	RC_TRY( allocate2D(6, 6, &D_matrix) );
	RC_TRY( allocate2D(6, 3*FEM_MAX_NODE, &B_matrix) );
	RC_TRY( allocate2D(3*FEM_MAX_NODE, 6, &trans_B) );
	RC_TRY( allocate2D(3*FEM_MAX_NODE, 6, &BtD) );

	RC_TRY( allocate_fem_bc_array(node.size, adj_force) );
	for(ii1=0; ii1<node.size; ii1++){
		if(node.array[ii1].label < 0) continue;

		adj_force->array[ii1].node = node.array[ii1].label;
	}
	adj_force->size = node.size;
	RC_TRY( clean_fem_bc_array(adj_force) );

	RC_TRY( mises_integration(node, element, material, physical, disp, rho,
	                          &m_integ, &volume) );

	for(ii1=0; ii1<element.size; ii1++){
		if( element.array[ii1].label < 0 ) continue;

		RC_TRY( Gauss_const_max(element.array[ii1].type,
		                            gauss_point, weight, &num_point) );
		RC_TRY( node_location_matrix(&(element.array[ii1]), node,
		                               node_location) ); 
		dim = element_dim(element.array[ii1].type);
		ssnum = element_ssnum(element.array[ii1].type);
		RC_TRY( get_D_matrix(element.array[ii1],
		                     material, physical, D_matrix) );
		RC_TRY( get_thickness(element.array[ii1], physical, &thickness) );
		for(ii2=0; ii2<dim*element.array[ii1].node_num; ii2++){
			force_v_sum[ii2] = 0.0;
		}

		mat_index = search_fem_material_prop_element(material, physical,
		                                             &(element.array[ii1]));
		RC_NEG_CHK( mat_index );
		if( material.array[mat_index].label >= IGNORE_MAT_LABEL ) continue;

		for(ii2=0; ii2<num_point; ii2++){
			RC_TRY( make_B_matrix(element.array[ii1], gauss_point[ii2],
			                      node_location, B_matrix, &det_J) );

			RC_TRY( fem_local_stress_strain(element.array[ii1],
			                 gauss_point[ii2], node, disp, material, physical,
			                 &local_stress, &local_strain) );

			mises_str = rho * mises_stress(local_stress);
			if(mises_str > exp_limit) mises_str = exp_limit;
			adj_strain = mul_scalar_trans_rotation3d(exp(mises_str)/m_integ,
			                                     d_mises_stress(local_stress));
			if(dim == 2){
				adj_strain_v[0] = adj_strain.x;
				adj_strain_v[1] = adj_strain.y;
				adj_strain_v[2] = adj_strain.xy;
			}else{   /* dim == 3 */
				adj_strain_v[0] = adj_strain.x;
				adj_strain_v[1] = adj_strain.y;
				adj_strain_v[2] = adj_strain.z;
				adj_strain_v[3] = adj_strain.yz;
				adj_strain_v[4] = adj_strain.zx;
				adj_strain_v[5] = adj_strain.xy;
			}
			transpose_matrix(ssnum, dim*element.array[ii1].node_num,
			                 B_matrix, trans_B);
			mul_matrix(dim*element.array[ii1].node_num, ssnum, ssnum,
			           trans_B, D_matrix, BtD);
			mul_matrix_vect(dim*element.array[ii1].node_num, ssnum,
			                BtD, adj_strain_v, force_v);
			for(ii3=0; ii3<dim*element.array[ii1].node_num; ii3++){
				force_v_sum[ii3] += force_v[ii3]*det_J*weight[ii2]*thickness;
			}
		}
		for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
			int index = search_fem_bc_node(*adj_force,
			                               element.array[ii1].node[ii2]);
			RC_NEG_CHK(index);
			adj_force->array[index].v.x += force_v_sum[dim*ii2];
			adj_force->array[index].v.y += force_v_sum[dim*ii2+1];
			if(dim == 3) adj_force->array[index].v.z += force_v_sum[dim*ii2+2];
			adj_force->array[index].v_type.x = BC_FORCE;
			adj_force->array[index].v_type.y = BC_FORCE;
			if(dim == 3) adj_force->array[index].v_type.z = BC_FORCE;
		}
	}

	RC_TRY( free2D(MAX_GAUSS_POINT, 4, gauss_point) );
	RC_TRY( free2D(FEM_MAX_NODE, 3, node_location) );
	RC_TRY( free2D(6, 6, D_matrix) );
	RC_TRY( free2D(6, 3*FEM_MAX_NODE, B_matrix) );
	RC_TRY( free2D(3*FEM_MAX_NODE, 6, trans_B) );
	RC_TRY( free2D(3*FEM_MAX_NODE, 6, BtD) );

	return(NORMAL_RC);
}


RC max_mises_stress(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                    FEM_MATERIAL_PROP_ARRAY material,
                    FEM_PHYSICAL_PROP_ARRAY physical,
                    FEM_DISP_ARRAY disp, double *max_val)
{
	int ii1, ii2;
	double **gauss_point = NULL;
	double weight[MAX_GAUSS_POINT];
	int num_point;
	FEM_STRESS_STRAIN stress, strain;
	double mises_str;
	int mat_index;
	double max_val_all;

	RC_TRY( allocate2D(MAX_GAUSS_POINT, 4, &gauss_point) );

	*max_val = ABS_TOL;
	max_val_all = ABS_TOL;
	for(ii1=0; ii1<element.size; ii1++){
		if( element.array[ii1].label < 0 ) continue;

		mat_index = search_fem_material_prop_element(material, physical,
		                                             &(element.array[ii1]));
		RC_NEG_CHK( mat_index );

		RC_TRY( Gauss_const_max(element.array[ii1].type,
		                            gauss_point, weight, &num_point) );
		for(ii2=0; ii2<num_point; ii2++){
			RC_TRY( fem_local_stress_strain(element.array[ii1],
			                 gauss_point[ii2], node, disp, material, physical,
			                 &stress, &strain) );
			mises_str = mises_stress(stress);
			if(mises_str > max_val_all) max_val_all = mises_str;
			if( material.array[mat_index].label >= IGNORE_MAT_LABEL ) continue;
			if(mises_str > *max_val) *max_val = mises_str;
		}
	}
/*	fprintf(stderr, "%15.7e -> %15.7e\n", max_val_all, *max_val); */

	RC_TRY( free2D(MAX_GAUSS_POINT, 4, gauss_point) );

	return(NORMAL_RC);
}


/* 熱応力を考慮したミーゼス応力の最大値 */
RC max_mises_stress_therm(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                          FEM_MATERIAL_PROP_ARRAY material,
                          FEM_PHYSICAL_PROP_ARRAY physical,
                          FEM_BC_ARRAY temp, FEM_DEFAULT_BC_ARRAY def_temp,
                          FEM_DISP_ARRAY disp, double *max_val)
{
	int ii1, ii2;
	double **node_location = NULL;
	double **gauss_point = NULL;
	double **dN = NULL;
	double **jacobi = NULL;
	double weight[MAX_GAUSS_POINT];
	int num_point;
	int dim;
	FEM_STRESS_STRAIN stress, strain;
	double det_J;
	double mises_str;

	RC_TRY( allocate2D(MAX_GAUSS_POINT, 4, &gauss_point) );
	RC_TRY( allocate2D(FEM_MAX_NODE, 3, &node_location) );
	RC_TRY( allocate2D(3, FEM_MAX_NODE, &dN) );
	RC_TRY( allocate2D(3, 3, &jacobi) );

	*max_val = 0.0;
	for(ii1=0; ii1<element.size; ii1++){
		if( element.array[ii1].label < 0 ) continue;

		RC_TRY( Gauss_const_max(element.array[ii1].type,
		                            gauss_point, weight, &num_point) );
		RC_TRY( node_location_matrix(&(element.array[ii1]), node,
		                               node_location) ); 
		dim = element_dim(element.array[ii1].type);
		for(ii2=0; ii2<num_point; ii2++){
			RC_TRY( set_dN_matrix(element.array[ii1].type,
			                      gauss_point[ii2], dN) );
			mul_matrix(dim, element.array[ii1].node_num,
			           dim, dN, node_location, jacobi);
			det_J = determinant(dim, jacobi);

			RC_TRY( fem_local_stress_strain_therm(element.array[ii1],
			                 gauss_point[ii2], node, disp, material, physical,
			                 temp, def_temp, 1.0, &stress, &strain) );
			mises_str = mises_stress(stress);
			if(mises_str > *max_val) *max_val = mises_str;
		}
	}

	RC_TRY( free2D(MAX_GAUSS_POINT, 4, gauss_point) );
	RC_TRY( free2D(FEM_MAX_NODE, 3, node_location) );
	RC_TRY( free2D(3, FEM_MAX_NODE, dN) );
	RC_TRY( free2D(3, 3, jacobi) );

	return(NORMAL_RC);
}


RC mises_integration(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                     FEM_MATERIAL_PROP_ARRAY material,
                     FEM_PHYSICAL_PROP_ARRAY physical,
                     FEM_DISP_ARRAY disp, double rho, double *integ_val,
                     double *vol)
{
	int ii1, ii2;
	double **node_location = NULL;
	double **gauss_point = NULL;
	double **dN = NULL;
	double **jacobi = NULL;
	double weight[MAX_GAUSS_POINT];
	int num_point;
	long double sum, sum_v;
	long double integ_val_long;
	int dim;
	FEM_STRESS_STRAIN stress, strain;
	double det_J;
	double mises_str;
	int mat_index;
	double integ_val_all, vol_all;

	RC_TRY( allocate2D(MAX_GAUSS_POINT, 4, &gauss_point) );
	RC_TRY( allocate2D(FEM_MAX_NODE, 3, &node_location) );
	RC_TRY( allocate2D(3, FEM_MAX_NODE, &dN) );
	RC_TRY( allocate2D(3, 3, &jacobi) );

	*vol = ABS_TOL;
	vol_all = ABS_TOL;
	integ_val_long = ABS_TOL;
	integ_val_all = ABS_TOL;
	for(ii1=0; ii1<element.size; ii1++){
		if( element.array[ii1].label < 0 ) continue;
		mat_index = search_fem_material_prop_element(material, physical,
		                                             &(element.array[ii1]));
		RC_NEG_CHK( mat_index );

		RC_TRY( Gauss_const_max(element.array[ii1].type,
		                            gauss_point, weight, &num_point) );
		RC_TRY( node_location_matrix(&(element.array[ii1]), node,
		                               node_location) ); 
		sum = 0.0;
		sum_v = 0.0;
		dim = element_dim(element.array[ii1].type);
		for(ii2=0; ii2<num_point; ii2++){
			RC_TRY( set_dN_matrix(element.array[ii1].type,
			                      gauss_point[ii2], dN) );
			mul_matrix(dim, element.array[ii1].node_num,
			           dim, dN, node_location, jacobi);
			det_J = determinant(dim, jacobi);

			RC_TRY( fem_local_stress_strain(element.array[ii1],
			                 gauss_point[ii2], node, disp, material, physical,
			                 &stress, &strain) );
			mises_str = rho * mises_stress(stress);
			sum += weight[ii2] * det_J * exp(mises_str);
			sum_v += weight[ii2] * det_J;
		}
		integ_val_all += sum;
		vol_all += sum_v;
		if( material.array[mat_index].label >= IGNORE_MAT_LABEL ) continue;
		integ_val_long += sum;
		*vol += sum_v;
	}
	*integ_val = (double)integ_val_long;
/*	fprintf(stderr, "%15.7e -> %15.7e (integ.)\n", integ_val_all, *integ_val);
	fprintf(stderr, "%15.7e -> %15.7e (vol.)\n", vol_all, *vol);*/

	RC_TRY( free2D(MAX_GAUSS_POINT, 4, gauss_point) );
	RC_TRY( free2D(FEM_MAX_NODE, 3, node_location) );
	RC_TRY( free2D(3, FEM_MAX_NODE, dN) );
	RC_TRY( free2D(3, 3, jacobi) );

	return(NORMAL_RC);
}


/* 節点が disp だけ変動した場合の汎関数の変動量を感度に対応する
   荷重 sens_force より予測 */
RC functional_variation(FEM_BC_ARRAY sens_force,
                        FEM_DISP_ARRAY disp, double *delta_f)
{
	/* 実際はコンプライアンスの計算と同様 */
	return( compliance_functional(sens_force, disp, delta_f) );
}


RC volume_functional(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                     double *volume)
{
	return( total_element_volume(element, node, volume) );
}


/* 外力による仕事 - 強制変位による仕事 */
RC potential_functional(FEM_BC_ARRAY force, FEM_DISP_ARRAY disp,
                        FEM_BC_ARRAY rest, FEM_REACT_ARRAY react,
                        double *potential)
{
	double comp1, comp2;

	RC_TRY( compliance_functional(force, disp, &comp1) );
	RC_TRY( restdisp_work(rest, react, &comp2) );

	*potential = comp1 - comp2;

	return(NORMAL_RC);
}


/* 強制変位による仕事 */
static RC restdisp_work(FEM_BC_ARRAY rest, FEM_REACT_ARRAY react,
                        double *compliance)
{
	int ii1;
	int index;
	long double sum = 0.0;

	for(ii1=0; ii1<rest.size; ii1++){
		if(rest.array[ii1].node < 0)continue;

		index = search_fem_react_node(react, rest.array[ii1].node);
		RC_NEG_CHK(index);

		if(rest.array[ii1].v_type.x == BC_FIX){
			sum += rest.array[ii1].v.x * react.array[index].v.x;
		}
		if(rest.array[ii1].v_type.y == BC_FIX){
			sum += rest.array[ii1].v.y * react.array[index].v.y;
		}
		if(rest.array[ii1].v_type.z == BC_FIX){
			sum += rest.array[ii1].v.z * react.array[index].v.z;
		}
	}
	*compliance = (double)sum;

	return(NORMAL_RC);
}


/* 外力仕事 */
RC compliance_functional(FEM_BC_ARRAY force, FEM_DISP_ARRAY disp,
                         double *compliance)
{
	int ii1;
	int index;
	long double sum = 0.0;

	for(ii1=0; ii1<force.size; ii1++){
		if(force.array[ii1].node < 0)continue;

		index = search_fem_disp_node(disp, force.array[ii1].node);
		RC_NEG_CHK(index);

		if(force.array[ii1].v_type.x == BC_FORCE){
			sum += force.array[ii1].v.x * disp.array[index].v.x;
		}
		if(force.array[ii1].v_type.y == BC_FORCE){
			sum += force.array[ii1].v.y * disp.array[index].v.y;
		}
		if(force.array[ii1].v_type.z == BC_FORCE){
			sum += force.array[ii1].v.z * disp.array[index].v.z;
		}
	}
	*compliance = (double)sum;

	return(NORMAL_RC);
}


/* 散逸エネルギー */
/* fluid_param : 無次元化された問題ではレイノルズ数の逆数、*/
/*               それ以外は粘性係数を代入 */
RC dissipation_functional(double fluid_param, FEM_NODE_ARRAY node,
                          FEM_ELEMENT_ARRAY element,
                          FEM_STRAIN_ARRAY strain_velo, double *dissipation)
{
	int ii1, ii2, ii3;
	int index;
	int g_point_num;
	int dim;
	double sum;
	double local_dissipation;
	double det_J;

	double **node_location;
	double **dN;
	double **Jacobi;
	double **g_points;
	double *N;
	double *weight;

	RC_TRY( allocate2D(FEM_MAX_NODE, 3, &node_location) );
	RC_TRY( allocate2D(3, FEM_MAX_NODE, &dN) );
	RC_TRY( allocate2D(3, 3, &Jacobi) );
	RC_TRY( allocate2D(MAX_GAUSS_POINT, 4, &g_points) );

	RC_TRY( allocate1D(FEM_MAX_NODE, &N) );
	RC_TRY( allocate1D(MAX_GAUSS_POINT, &weight) );

	*dissipation = 0.0;

	dim = fem_analysis_dim(element);

	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		RC_TRY( Gauss_const_max(element.array[ii1].type, g_points, weight,
		                        &g_point_num) );
		RC_TRY( node_location_matrix(&(element.array[ii1]), node,
		                                                    node_location) );

		local_dissipation = 0.0;
		for(ii2=0; ii2<g_point_num; ii2++){
			RC_TRY( set_N_vector(element.array[ii1].type, g_points[ii2], N) );
			RC_TRY( set_dN_matrix(element.array[ii1].type, g_points[ii2], dN) );
			mul_matrix(dim, element.array[ii1].node_num, dim,
			           dN, node_location, Jacobi);
			det_J = determinant(dim, Jacobi);
			RC_NEG_CHK(det_J);

			sum = 0.0;
			for(ii3=0; ii3<element.array[ii1].node_num; ii3++){
				index = search_fem_strain_element_node(strain_velo,
				                                  element.array[ii1].label,
				                                  element.array[ii1].node[ii3]);
				RC_NEG_CHK(index);
				sum += (2.0 * strain_velo.array[index].v.x
				            * strain_velo.array[index].v.x
				      + 2.0 * strain_velo.array[index].v.y
				            * strain_velo.array[index].v.y
				      + 2.0 * strain_velo.array[index].v.z
				            * strain_velo.array[index].v.z
				      +       strain_velo.array[index].v.xy
				            * strain_velo.array[index].v.xy
				      +       strain_velo.array[index].v.yz
				            * strain_velo.array[index].v.yz
				      +       strain_velo.array[index].v.zx
				            * strain_velo.array[index].v.zx)*fluid_param*N[ii3];
			}
			local_dissipation += sum * weight[ii2] * det_J;
		}
		*dissipation += local_dissipation;
	}

	RC_TRY( free2D(FEM_MAX_NODE, 3, node_location) );
	RC_TRY( free2D(3, FEM_MAX_NODE, dN) );
	RC_TRY( free2D(3, 3, Jacobi) );
	RC_TRY( free2D(MAX_GAUSS_POINT, 4, g_points) );
	free( (void *)N );
	free( (void *)weight );

	return(NORMAL_RC);
}


RC dens_sens_comp_disp(const double power, FEM_DISP_ARRAY disp,
                       FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                       FEM_MATERIAL_PROP_ARRAY material,
                       FEM_PHYSICAL_PROP_ARRAY physical,
                       FEM_SENSITIVITY_ARRAY *sens,
                       FEM_MATERIAL_PROP_ARRAY dens)
{
	int ii1, ii2;
	int dim;
	static int num_point;
	double sum;
	double sens_of_point;
	double det_J;
	double elem_volume;
	FEM_WORK_SET ws;
	FEM_STRESS_STRAIN local_stress, local_strain;

	RC_TRY( allocate_fem_work_set(&ws) );

	RC_TRY( allocate_fem_sensitivity_array(0, sens) );
	sens->size = 0;
	sens->type = SENS_SHAPE_GRAD_DENSITY;

	RC_NEG_CHK( dim = fem_analysis_dim(element) );
	for(ii1=0; ii1<element.size; ii1++){
		double *weights = ws.vec_G[0];
		double **gauss_point = ws.mat_G_4[0];
		double **node_location = ws.mat_N_3[0];
		double **dN = ws.mat_3_N[0];
		double **jacobi = ws.mat_3_3[0];

		if( element.array[ii1].label < 0 ) continue;

		RC_TRY( Gauss_const_default(element.array[ii1].type, gauss_point,
		                            weights, &num_point) );
		RC_TRY( node_location_matrix(&(element.array[ii1]), node,
		                             node_location) );
		RC_TRY( realloc_fem_sensitivity_array(sens) );
		sum = 0.0;
		elem_volume = 0.0;
		for(ii2=0; ii2<num_point; ii2++){
			init_matrix(3, FEM_MAX_NODE, dN);
			RC_TRY( set_dN_matrix(element.array[ii1].type,
			                      gauss_point[ii2], dN) );
			mul_matrix(3, element.array[ii1].node_num, 3, dN,
			           node_location, jacobi);
			RC_NEG_CHK( dim );
			det_J = determinant(dim, jacobi);
			if( det_J < 0.0){
				fprintf(stderr, "warning: negative determinant Jacobian "
				                "[%d-%d]\n", element.array[ii1].label, ii1);
			}
			RC_TRY( fem_local_stress_strain(element.array[ii1],
			                                gauss_point[ii2], node, disp,
			                                material, physical, &local_stress,
			                                &local_strain) );
			sens_of_point = - (local_stress.v.x * local_strain.v.x
			                 + local_stress.v.y * local_strain.v.y
			                 + local_stress.v.z * local_strain.v.z
			                + local_stress.v.yz * local_strain.v.yz
			                + local_stress.v.zx * local_strain.v.zx
			                + local_stress.v.xy * local_strain.v.xy);
			sens_of_point *= ( power * pow(dens.array[ii1].rho, power - 1.0) );
			sum += sens_of_point * det_J * weights[ii2];
			elem_volume += 1.0 * det_J * weights[ii2];
		}
		sens->array[sens->size].element = element.array[ii1].label;
		sens->array[sens->size].s = sum / elem_volume;
		sens->size++;
	}
	RC_TRY( clean_fem_sensitivity_array(sens) );

	RC_TRY( free_fem_work_set(&ws) );

	return(NORMAL_RC);
}


/* 熱歪による仕事の感度（簡易テスト版） */
RC shape_sens_therm(FEM_ELEMENT_ARRAY surf, FEM_DISP_ARRAY disp,
                    FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                    FEM_MATERIAL_PROP_ARRAY material,
                    FEM_PHYSICAL_PROP_ARRAY physical,
                    FEM_BC_ARRAY temp, FEM_DEFAULT_BC_ARRAY def_temp,
                    FEM_SENSITIVITY_ARRAY *sens)
{
	int ii1, ii2, ii3;
	static int init_flag = 0;
	static double **local_points;
	FEM_STRESS_STRAIN local_stress, local_strain;
	FEM_STRESS_STRAIN local_stress2, local_strain2;
	int elem_index;

	if(init_flag == 0){
		RC_TRY( allocate2D(FEM_MAX_NODE, 4, &local_points) );
		init_flag = 1;
	}

	RC_TRY( allocate_fem_sensitivity_array(0, sens) );
	sens->size = 0;
	sens->type = SENS_SHAPE_GRAD_DENSITY;

	for(ii1=0; ii1<surf.size; ii1++){
		if(surf.array[ii1].label < 0) continue;

		elem_index = search_fem_element_label(element,
		                                      surf.array[ii1].i_info[0]);
		if(elem_index < 0) return(SEEK_ERROR_RC);

		RC_TRY( set_local_node_points(element.array[elem_index].type,
		                                      local_points) );
		for(ii2=0; ii2<surf.array[ii1].node_num; ii2++){
			RC_TRY( realloc_fem_sensitivity_array(sens) );
			sens->array[sens->size].element = surf.array[ii1].label;
			sens->array[sens->size].node = surf.array[ii1].node[ii2];

			for(ii3=0; ii3<element.array[elem_index].node_num; ii3++){
				if(surf.array[ii1].node[ii2]
				   == element.array[elem_index].node[ii3]) break;
			}
			if(ii3 >= element.array[elem_index].node_num){
				return(SEEK_ERROR_RC);
			}
			RC_TRY( fem_local_stress_strain(element.array[elem_index],
			           local_points[ii3], node, disp, material, physical,
			           &local_stress, &local_strain) );
			RC_TRY( fem_local_stress_strain_therm(element.array[elem_index],
			           local_points[ii3], node, disp, material, physical,
			           temp, def_temp, 2.0, &local_stress2, &local_strain2) );
			sens->array[sens->size].s =  (local_stress.v.x * local_strain.v.x
			                          + local_stress.v.y * local_strain.v.y
			                          + local_stress.v.z * local_strain.v.z
			                          + local_stress.v.yz * local_strain.v.yz
			                          + local_stress.v.zx * local_strain.v.zx
			                          + local_stress.v.xy * local_strain.v.xy);
			sens->size++;
		}
	}

	RC_TRY( clean_fem_sensitivity_array(sens) );

	return(NORMAL_RC);
}


/* 変位 -> 応力、歪 -> 平均コンプライアンスの感度 */
RC shape_sens_compliance_disp(FEM_ELEMENT_ARRAY surf, FEM_DISP_ARRAY disp,
                              FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                              FEM_MATERIAL_PROP_ARRAY material,
                              FEM_PHYSICAL_PROP_ARRAY physical,
                              FEM_SENSITIVITY_ARRAY *sens)
{
	int ii1, ii2, ii3;
	static int init_flag = 0;
	static double **local_points;
	FEM_STRESS_STRAIN local_stress, local_strain;
	int elem_index;

	if(init_flag == 0){
		RC_TRY( allocate2D(FEM_MAX_NODE, 4, &local_points) );
		init_flag = 1;
	}

	RC_TRY( allocate_fem_sensitivity_array(0, sens) );
	sens->size = 0;
	sens->type = SENS_SHAPE_GRAD_DENSITY;

	for(ii1=0; ii1<surf.size; ii1++){
		if(surf.array[ii1].label < 0) continue;

		elem_index = search_fem_element_label(element,
		                                      surf.array[ii1].i_info[0]);
		if(elem_index < 0) return(SEEK_ERROR_RC);

		RC_TRY( set_local_node_points(element.array[elem_index].type,
		                                      local_points) );
		for(ii2=0; ii2<surf.array[ii1].node_num; ii2++){
			RC_TRY( realloc_fem_sensitivity_array(sens) );
			sens->array[sens->size].element = surf.array[ii1].label;
			sens->array[sens->size].node = surf.array[ii1].node[ii2];

			for(ii3=0; ii3<element.array[elem_index].node_num; ii3++){
				if(surf.array[ii1].node[ii2]
				   == element.array[elem_index].node[ii3]) break;
			}
			if(ii3 >= element.array[elem_index].node_num){
				return(SEEK_ERROR_RC);
			}
			RC_TRY( fem_local_stress_strain(element.array[elem_index],
			           local_points[ii3], node, disp, material, physical,
			           &local_stress, &local_strain) );
			sens->array[sens->size].s = - (local_stress.v.x * local_strain.v.x
			                           + local_stress.v.y * local_strain.v.y
			                           + local_stress.v.z * local_strain.v.z
			                           + local_stress.v.yz * local_strain.v.yz
			                           + local_stress.v.zx * local_strain.v.zx
			                           + local_stress.v.xy * local_strain.v.xy);
			sens->size++;
		}
	}

	RC_TRY( clean_fem_sensitivity_array(sens) );

	return(NORMAL_RC);
}


/* 応力、歪 -> 平均コンプライアンスの感度 */
RC shape_sens_compliance_ss(FEM_ELEMENT_ARRAY surf,
                            FEM_STRESS_ARRAY stress, FEM_STRAIN_ARRAY strain,
                            FEM_SENSITIVITY_ARRAY *sens)
{
	int ii1, ii2;
	int index_stress;
	int index_strain;

	if(surf.source != FEM_GENERATED_SURFACE){
		fprintf(stderr, "surf.source < \n");
		return(ARG_ERROR_RC);
	}
	
	RC_TRY( allocate_fem_sensitivity_array(0, sens) );
	sens->size = 0;
	sens->type = SENS_SHAPE_GRAD_DENSITY;

	for(ii1=0; ii1<surf.size; ii1++){
		if(surf.array[ii1].label < 0) continue;

		for(ii2=0; ii2<surf.array[ii1].node_num; ii2++){
			RC_TRY( realloc_fem_sensitivity_array(sens) );

			sens->array[sens->size].element = surf.array[ii1].label;
			sens->array[sens->size].node = surf.array[ii1].node[ii2];

			index_stress = search_fem_stress_element_node(stress,
			                           surf.array[ii1].i_info[0],
			                           surf.array[ii1].node[ii2]);
			index_strain = search_fem_strain_element_node(strain,
			                           surf.array[ii1].i_info[0],
			                           surf.array[ii1].node[ii2]);
			if((index_stress < 0)||(index_strain < 0)){
				fprintf(stderr,
				        "search_fem_[stress|strain]_element_node() < \n");
				return(SEEK_ERROR_RC);
			}

			sens->array[sens->size].s = - (stress.array[index_stress].v.x
			                              *strain.array[index_strain].v.x
			                             + stress.array[index_stress].v.y
			                              *strain.array[index_strain].v.y
			                             + stress.array[index_stress].v.z
			                              *strain.array[index_strain].v.z
			                             + stress.array[index_stress].v.yz
			                              *strain.array[index_strain].v.yz
			                             + stress.array[index_stress].v.zx
			                              *strain.array[index_strain].v.zx
			                             + stress.array[index_stress].v.xy
			                              *strain.array[index_strain].v.xy);
			(sens->size)++;
		}
	}

	RC_TRY( clean_fem_sensitivity_array(sens) );

	return(NORMAL_RC);
}


/* strain : ひずみ流速 */
RC shape_sens_dissipation(double fluid_param, FEM_ELEMENT_ARRAY surf,
                          FEM_VELO_ARRAY flow, FEM_VELO_ARRAY adjflow,
                          FEM_STRAIN_ARRAY strain,
                          FEM_STRAIN_ARRAY adjstrain,
                          FEM_SENSITIVITY_ARRAY *sens)
{
	int ii1, ii2;
	int index_v;
	int index_w;
	int index_p;
	int index_q;
	double tmp;

	if(surf.source != FEM_GENERATED_SURFACE){
		fprintf(stderr, "surf.source < \n");
		return(ARG_ERROR_RC);
	}
	RC_TRY( allocate_fem_sensitivity_array(0, sens) );
	sens->size = 0;
	sens->type = SENS_SHAPE_GRAD_DENSITY;
	for(ii1=0; ii1<surf.size; ii1++){
		if(surf.array[ii1].label < 0) continue;
		for(ii2=0; ii2<surf.array[ii1].node_num; ii2++){
			RC_TRY( realloc_fem_sensitivity_array(sens) );
			sens->array[sens->size].element = surf.array[ii1].label;
			sens->array[sens->size].node = surf.array[ii1].node[ii2];
			index_v = search_fem_strain_element_node(strain,
			                                         surf.array[ii1].i_info[0],
			                                         surf.array[ii1].node[ii2]);
			RC_NEG_CHK(index_v);
			index_w = search_fem_strain_element_node(adjstrain,
			                                         surf.array[ii1].i_info[0],
			                                         surf.array[ii1].node[ii2]);
			RC_NEG_CHK(index_w);
			index_p = search_fem_velo_node(flow, surf.array[ii1].node[ii2]);
			RC_NEG_CHK(index_p);
			index_q = search_fem_velo_node(adjflow, surf.array[ii1].node[ii2]);
			RC_NEG_CHK(index_q);

			tmp = ( 2.0 * strain.array[index_v].v.x
			            * strain.array[index_v].v.x
			      + 2.0 * strain.array[index_v].v.y
			            * strain.array[index_v].v.y
			      + 2.0 * strain.array[index_v].v.z
			            * strain.array[index_v].v.z
			            + strain.array[index_v].v.xy
			            * strain.array[index_v].v.xy
			            + strain.array[index_v].v.yz
			            * strain.array[index_v].v.yz
			            + strain.array[index_v].v.zx
			            * strain.array[index_v].v.zx ) * fluid_param
			    - ( 2.0 * strain.array[index_v].v.x
			            * adjstrain.array[index_w].v.x
			      + 2.0 * strain.array[index_v].v.y
			            * adjstrain.array[index_w].v.y
			      + 2.0 * strain.array[index_v].v.z
			            * adjstrain.array[index_w].v.z
			            + strain.array[index_v].v.xy
			            * adjstrain.array[index_w].v.xy
			            + strain.array[index_v].v.yz
			            * adjstrain.array[index_w].v.yz
			            + strain.array[index_v].v.zx
			            * adjstrain.array[index_w].v.zx ) * fluid_param
			     + ( adjstrain.array[index_w].v.x
			       + adjstrain.array[index_w].v.y
			       + adjstrain.array[index_w].v.z ) * flow.array[index_p].s
			     + ( strain.array[index_v].v.x
			       + strain.array[index_v].v.y
			       + strain.array[index_v].v.z ) * adjflow.array[index_q].s;
			sens->array[sens->size].s = - tmp;
			(sens->size)++;
		}
	}

	RC_TRY( clean_fem_sensitivity_array(sens) );

	return(NORMAL_RC);
}


RC shape_sens_volume(FEM_ELEMENT_ARRAY surf, FEM_SENSITIVITY_ARRAY *sens)
{
	int ii1, ii2;

	RC_TRY( allocate_fem_sensitivity_array(0, sens) );
	sens->size = 0;
	sens->type = SENS_SHAPE_GRAD_DENSITY;

	for(ii1=0; ii1<surf.size; ii1++){
		if(surf.array[ii1].label < 0) continue;

		for(ii2=0; ii2<surf.array[ii1].node_num; ii2++){
			RC_TRY( realloc_fem_sensitivity_array(sens) );

			sens->array[sens->size].element = surf.array[ii1].label;
			sens->array[sens->size].node = surf.array[ii1].node[ii2];
			sens->array[sens->size].s = 1.0;

			(sens->size)++;
		}
	}

	RC_TRY( clean_fem_sensitivity_array(sens) );

	return(NORMAL_RC);
}


/* 感度（形状勾配密度関数）に対応する節点荷重 */
/* sens は surf の各節点で与える */
/* trac は自動的に確保される */
RC sens2traction(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY surf,
                 FEM_BC_ARRAY shape_rest, FEM_SENSITIVITY_ARRAY sens,
                 FEM_BC_ARRAY *trac, int reduce_flag)
{
	int ii1, ii2;
	double p[FEM_MAX_NODE];
	int p_size;
	TRANSLATION3D f[FEM_MAX_NODE];
	int index;
	FEM_ELEMENT tmp_surf;
	FEM_ELEM_TYPE tmp_type;
	int tmp_node_num;

	if(sens.type != SENS_SHAPE_GRAD_DENSITY) return(ARG_ERROR_RC);
	
	RC_TRY( allocate_fem_bc_array(node.size, trac) );
	trac->size = node.size;
	for(ii1=0; ii1<(trac->size); ii1++){
		trac->array[ii1].node = node.array[ii1].label;
	}
	RC_TRY( sort_fem_bc_array(trac) );

	for(ii1=0; ii1<(surf.size); ii1++){
		if(surf.array[ii1].label < 0) continue;

		/* sens -> p */
		p_size = surf.array[ii1].node_num;
		for(ii2=0; ii2<p_size; ii2++){
			index = search_fem_sensitivity_element_node(sens,
			                           surf.array[ii1].label,
			                           surf.array[ii1].node[ii2]);
			if(index < 0){
				/* サーチ失敗、２次要素の場合で
				   １次要素分までサーチできていれば p_size を適当に変更 */
				if(surf.array[ii1].type == ELEM_TRI2){
					if(ii2 >= 3){
						p_size = 3;
						break;
					}
				}else if(surf.array[ii1].type == ELEM_QUAD2){
					if(ii2 >= 4){
						p_size = 4;
						break;
					}
				}else if( (surf.array[ii1].type == ELEM_BEAM2)
				        ||(surf.array[ii1].type == ELEM_LINE2) ){
					if(ii2 >= 2){
						p_size = 2;
						break;
					}
				}
				p_size = 1;   /* 要素中心で再度サーチ */
				break;
			}
			p[ii2] = sens.array[index].s;
		}
		if(p_size == 1){
			index = search_fem_sensitivity_element_node(sens,
			                           surf.array[ii1].label, -1);
			if(index < 0){
				fprintf(stderr, "search_fem_sensitivity_element_node() < \n");
				return(SEEK_ERROR_RC);
			}
		}

		tmp_surf = surf.array[ii1];
		if(reduce_flag){
			if(element_order(tmp_surf.type) == 2){
				RC_TRY( reduced_element_type(tmp_surf.type, &tmp_type,
				                             &tmp_node_num) );
				tmp_surf.type = tmp_type;
				tmp_surf.node_num = tmp_node_num;
			}
		}

		/* p -> f */
		RC_TRY( pressure_force(p_size, p, &tmp_surf, node, f) );

		/* f -> trac */
		for(ii2=0; ii2<tmp_surf.node_num; ii2++){
			index = search_fem_bc_node(*trac, tmp_surf.node[ii2]);
			if(index < 0){
				fprintf(stderr, "search_fem_bc_node < ");
				return(SEEK_ERROR_RC);
			}
			trac->array[index].v_type.x = BC_FORCE;
			trac->array[index].v_type.y = BC_FORCE;
			trac->array[index].v_type.z = BC_FORCE;
			trac->array[index].v.x += f[ii2].x;
			trac->array[index].v.y += f[ii2].y;
			trac->array[index].v.z += f[ii2].z;
		}
	}

	/* 形状拘束を受ける自由度の荷重をゼロにする */
	for(ii1=0; ii1<(shape_rest.size); ii1++){
		if(shape_rest.array[ii1].node < 0) continue;

		index = search_fem_bc_node(*trac, shape_rest.array[ii1].node);
		if(index < 0){
			fprintf(stderr, "search_fem_bc_node < ");
			return(SEEK_ERROR_RC);
		}

		if(shape_rest.array[ii1].v_type.x == BC_FIX){
			trac->array[index].v_type.x = BC_FREE;
			trac->array[index].v.x = 0.0;
		}
		if(shape_rest.array[ii1].v_type.y == BC_FIX){
			trac->array[index].v_type.y = BC_FREE;
			trac->array[index].v.y = 0.0;
		}
		if(shape_rest.array[ii1].v_type.z == BC_FIX){
			trac->array[index].v_type.z = BC_FREE;
			trac->array[index].v.z = 0.0;
		}
	}

	/* 全自由度に荷重を受けない節点を荷重条件より除外 */
	for(ii1=0; ii1<(trac->size); ii1++){
		if(trac->array[ii1].node < 0) continue;

		if( (trac->array[ii1].v_type.x == BC_FREE)
		  &&(trac->array[ii1].v_type.y == BC_FREE)
		  &&(trac->array[ii1].v_type.z == BC_FREE) ){
			trac->array[ii1].node = -1; 
		}
	}
	RC_TRY( clean_fem_bc_array(trac) );

	return(NORMAL_RC);
}


/* 変位 disps[disp_num] の factors[disp_num] 倍を節点 node の座標に加算 */
RC add_disp_node(FEM_NODE_ARRAY node, int disp_num, const double *factors,
                 const FEM_DISP_ARRAY *disps)
{
	int ii1, ii2;
	int index;
	TRANSLATION3D *sum;

	sum = (TRANSLATION3D *)malloc((node.size)*sizeof(TRANSLATION3D));
	if(sum == NULL){
		fprintf(stderr, "malloc < ");
		return(ALLOC_ERROR_RC);
	}
	for(ii1=0; ii1<node.size; ii1++){
		sum[ii1].x = 0.0;
		sum[ii1].y = 0.0;
		sum[ii1].z = 0.0;
	}

	/* 桁落ちを避けるために、加算する変位を一時的に sum に格納 */
	for(ii1=0; ii1<disp_num; ii1++){
		for(ii2=0; ii2<disps[ii1].size; ii2++){
			if(disps[ii1].array[ii2].node < 0) continue;

			index = search_fem_node_label(node, disps[ii1].array[ii2].node);
			if(index < 0){
				fprintf(stderr, "search_fem_node_label < ");
				return(SEEK_ERROR_RC);
			}
			sum[index].x += factors[ii1] * disps[ii1].array[ii2].v.x;
			sum[index].y += factors[ii1] * disps[ii1].array[ii2].v.y;
			sum[index].z += factors[ii1] * disps[ii1].array[ii2].v.z;
		}
	}

	/* sum を node に加算 */
	for(ii1=0; ii1<node.size; ii1++){
		node.array[ii1].p.x += sum[ii1].x;
		node.array[ii1].p.y += sum[ii1].y;
		node.array[ii1].p.z += sum[ii1].z;
	}

	free((void *)sum);

	return(NORMAL_RC);
}


RC max_pr_strain_disp(int str_num, const double factors[],
                      FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                      const FEM_DISP_ARRAY disp[], double *max_str,
                      int fem_sol_ss_flag)
{
	int ii1, ii2, ii3;
	double **local_points;
	FEM_STRESS_STRAIN local_strain;
	TRANS_ROTATION3D ss_v;
	double pr_val[3];
	double abs_v;

	RC_TRY( allocate2D(FEM_MAX_NODE, 4, &local_points) );

	*max_str = 0.0;
	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		if(fem_sol_ss_flag & FEM_SOL_SS_CENTER){
			RC_TRY( set_center_point(element.array[ii1].type,
			                         local_points[0]) );
			ss_v.x = ss_v.y = ss_v.z = ss_v.xy = ss_v.yz = ss_v.zx = 0.0;
			for(ii2=0; ii2<str_num; ii2++){
				RC_TRY( fem_local_strain(element.array[ii1],
				            local_points[0], node, disp[ii2], 
				            &local_strain) );
				ss_v.x += factors[ii2] * local_strain.v.x;
				ss_v.y += factors[ii2] * local_strain.v.y;
				ss_v.z += factors[ii2] * local_strain.v.z;
				ss_v.yz += factors[ii2] * local_strain.v.yz;
				ss_v.zx += factors[ii2] * local_strain.v.zx;
				ss_v.xy += factors[ii2] * local_strain.v.xy;
			}
			RC_TRY( principal3d(ss_v, pr_val, NULL, NULL, NULL) );
			for(ii2=0; ii2<3; ii2++){
				abs_v = fabs(pr_val[ii2]);
				if(abs_v > (*max_str)) (*max_str) = abs_v;
			}
		}
		if(fem_sol_ss_flag & FEM_SOL_SS_NODE){
			RC_TRY( set_local_node_points(element.array[ii1].type,
			                              local_points) );
			for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
				ss_v.x = ss_v.y = ss_v.z = ss_v.xy = ss_v.yz = ss_v.zx = 0.0;
				for(ii3=0; ii3<str_num; ii3++){
					RC_TRY( fem_local_strain(element.array[ii1],
					            local_points[ii2], node, disp[ii3],
					            &local_strain) );
					ss_v.x += factors[ii3] * local_strain.v.x;
					ss_v.y += factors[ii3] * local_strain.v.y;
					ss_v.z += factors[ii3] * local_strain.v.z;
					ss_v.yz += factors[ii3] * local_strain.v.yz;
					ss_v.zx += factors[ii3] * local_strain.v.zx;
					ss_v.xy += factors[ii3] * local_strain.v.xy;
				}
				RC_TRY( principal3d(ss_v, pr_val, NULL, NULL, NULL) );
				for(ii3=0; ii3<3; ii3++){
					abs_v = fabs(pr_val[ii3]);
					if(abs_v > (*max_str)) (*max_str) = abs_v;
				}
			}
		}
	}

	RC_TRY( free2D(FEM_MAX_NODE, 4, local_points) );

	return(NORMAL_RC);
}


/* ひずみ strain[str_num] の factors[str_num] 倍を重ね合わせた際の */
/* 最大主歪みの絶対値 */
/* strain[x].array[*].element と strain[x].array[*].node の値は x によらず */
/* 一定でなければならない */
RC max_pr_strain(int str_num, const double *factors,
                 const FEM_STRAIN_ARRAY *strain, double *max_str)
{
	int ii1, ii2;
	FEM_STRESS_STRAIN ss;
	double abs_v;
	double pr_val[3];

	if(str_num <= 0) return(ARG_ERROR_RC);
	for(ii1=1; ii1<str_num; ii1++){
		if(strain[ii1].size != strain[0].size) return(ARG_ERROR_RC);
	}

	*max_str = 0.0;
	for(ii1=0; ii1<(strain[0].size); ii1++){
		if(strain[0].array[ii1].element < 0) continue;

		ss = strain[0].array[ii1];
		ss.v.x = 0.0;
		ss.v.y = 0.0;
		ss.v.z = 0.0;
		ss.v.yz = 0.0;
		ss.v.zx = 0.0;
		ss.v.xy = 0.0;

		for(ii2=0; ii2<str_num; ii2++){
			if(strain[ii2].array[ii1].element != ss.element){
				return(ARG_ERROR_RC);
			}
			if(strain[ii2].array[ii1].node != ss.node){
				return(ARG_ERROR_RC);
			}

			ss.v.x += factors[ii2] * strain[ii2].array[ii1].v.x;
			ss.v.y += factors[ii2] * strain[ii2].array[ii1].v.y;
			ss.v.z += factors[ii2] * strain[ii2].array[ii1].v.z;
			ss.v.yz += factors[ii2] * strain[ii2].array[ii1].v.yz;
			ss.v.zx += factors[ii2] * strain[ii2].array[ii1].v.zx;
			ss.v.xy += factors[ii2] * strain[ii2].array[ii1].v.xy;
		}
		RC_TRY( principal3d(ss.v, pr_val, NULL, NULL, NULL) );

		for(ii2=0; ii2<3; ii2++){
			abs_v = fabs(pr_val[ii2]);
			if(abs_v > (*max_str)) (*max_str) = abs_v;
		}
	}

	return(NORMAL_RC);
}


RC make_trac_material_const(FEM_MATERIAL_PROP_ARRAY material,
                            FEM_MATERIAL_PROP_ARRAY *material_trac)
{
	int ii1;

	RC_TRY( allocate_fem_material_prop_array(material.size, material_trac) );
	RC_TRY( copy_fem_material_prop_array(material, material_trac) );

	for(ii1=0; ii1<material_trac->size; ii1++){
		if(material_trac->array[ii1].label < 0) continue;

		material_trac->array[ii1].mat_type = MAT_ISOTROPIC;
		material_trac->array[ii1].E = 1000.0;
		material_trac->array[ii1].nu = 0.0;
		material_trac->array[ii1].G = material_trac->array[ii1].E/2.0;
		material_trac->array[ii1].alpha = 0.0;
	}

	return(NORMAL_RC);
}


RC make_trac_material(FEM_MATERIAL_PROP_ARRAY material,
                      FEM_MATERIAL_PROP_ARRAY *material_trac)
{
	int ii1;

	RC_TRY( allocate_fem_material_prop_array(material.size, material_trac) );
	RC_TRY( copy_fem_material_prop_array(material, material_trac) );

	for(ii1=0; ii1<material_trac->size; ii1++){
		if(material_trac->array[ii1].label < 0) continue;
		switch(material_trac->array[ii1].mat_type){
		case MAT_ORTHOTROPIC:
			material_trac->array[ii1].E
			  = (material_trac->array[ii1].E_vect.x
			   + material_trac->array[ii1].E_vect.y
			   + material_trac->array[ii1].E_vect.z)/3.0;
			break;
		case MAT_ANISOTROPIC:
			material_trac->array[ii1].E
			  = (material_trac->array[ii1].D_matrix[0][0]
			   + material_trac->array[ii1].D_matrix[1][1]
			   + material_trac->array[ii1].D_matrix[2][2])/3.0;
			break;
		default:
			break;
		}
		material_trac->array[ii1].mat_type = MAT_ISOTROPIC;
		material_trac->array[ii1].nu = 0.0;
		material_trac->array[ii1].G = material_trac->array[ii1].E/2.0;
		material_trac->array[ii1].alpha = 0.0;
	}

	return(NORMAL_RC);
}
