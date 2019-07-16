/*********************************************************************
 * fem_solver_shell.c
 *
 * Copyright (C) 2003 AzLib Developers Group
 *
 * Written by
 *  <Yasuo ASAGA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: fem_solver_shell.c 415 2005-07-15 06:04:24Z sasaoka $*/

#include <stdio.h>
#include <stdlib.h>
#include "fem_solver.h"
#include "fem_struct.h"
#include "math_utl.h"
#include "rc.h"
#include "nonzero_cg.h"

#define DRILLING         (0.03)  /* ドリリング回転剛性のパラメータ */
#define REDUCED_INTEGRAL (0)     /* 0:完全積分 */
                                 /* 1:次数低減積分 */
                                 /* 2:選択型低減積分 */


static RC set_B_matrix_shell(int dim, int node_num,
                             ANALYSIS_TYPE ana_type, double *N,
                             double **dNxyz, double **B_matrix);
static RC trans_coord_node_location(NODE_ARRAY node, ELEMENT element,
                                    WORK_SET ws, double **coord);
static RC set_Gauss_const4RI(ELEM_TYPE type, double **points,
                             double *weight, int *num_point);
static RC set_drilling(int dim, ELEMENT element,
                       MATERIAL_PROP_ARRAY material,
                       PHYSICAL_PROP_ARRAY physical,
                       ELEM_MATRIX *elem_matrix);
static RC trans_coord_elem_matrix(NODE_ARRAY node, ELEMENT element,
                                  WORK_SET ws,
                                  ELEM_MATRIX *elem_matrix);
static RC fill_disp_vector4shell(int dim, ANALYSIS_TYPE ana_type,
                                 WORK_SET ws, NODE_ARRAY node,
                                 ELEMENT element, DISP_ARRAY disp,
                                 double disp_v[]);
static RC vect2stress_strain_shell(ANALYSIS_TYPE ana_type,
                                   const double *vect, STRESS_STRAIN *ss);
static RC vect2resultant_stress_strain(ANALYSIS_TYPE ana_type,
                                       const double *vect,
                                       STRESS_STRAIN *ss);
static RC fill_physical(MATERIAL_PROP *prop, ELEMENT element,
                        PHYSICAL_PROP_ARRAY physical);
static int plane_element_ssnum(ELEMENT element,
                               MATERIAL_PROP_ARRAY material,
                               PHYSICAL_PROP_ARRAY physical);
static int plane_element_dim(ELEMENT element,
                             MATERIAL_PROP_ARRAY material,
                             PHYSICAL_PROP_ARRAY physical);
static RC set_mid_point4quad(double **coords);
static RC set_elem_dcosine(NODE_ARRAY node, ELEMENT element,
                           TRANSLATION3D d_cos[]);
static void vect2dcosine_array(TRANSLATION3D v[], double **array);


/* 種々の剛性（膜、曲げ、せん断）マトリクスを */
/* 対応する箇所に組み込み、要素剛性マトリクスを作成（シェル解析用） */
RC fem_impose_shell_elem_matrix(int plane_dim, ELEMENT element,
                                MATERIAL_PROP_ARRAY material,
                                PHYSICAL_PROP_ARRAY physical,
                                ELEM_MATRIX elem_matrix,
                                ELEM_MATRIX *total_elem_matrix)
{
	int fact;
	int offset;
	int pos1, pos2;
	int ii1, ii2, ii3, ii4;

	switch(plane_dim){
	case 2:
		fact = 3;
		offset = 0;
		break;
	case 3:
		fact = 2;
		offset = 2;
		break;
	default:
		return(ARG_ERROR_RC);
	}

	for(ii1=0; ii1<elem_matrix.size; ii1+=plane_dim){
		pos1 = fact*elem_matrix.index[ii1] + offset;
		for(ii2=0; ii2<elem_matrix.size; ii2+=plane_dim){
			pos2 = fact*elem_matrix.index[ii2] + offset;
			for(ii3=0; ii3<plane_dim; ii3++){
				for(ii4=0; ii4<plane_dim; ii4++){
					total_elem_matrix->matrix[ii3+pos1][ii4+pos2]
					    += elem_matrix.matrix[ii1+ii3][ii2+ii4];
				}
			}
		}
	}

	return(NORMAL_RC);
}


/* 要素剛性マトリックスの作成（シェル解析用） */
RC fem_make_shell_matrix(ELEMENT element, NODE_ARRAY node,
                         MATERIAL_PROP_ARRAY material,
                         PHYSICAL_PROP_ARRAY physical,
                         ELEM_MATRIX *elem_matrix)
{
	int ii1, ii2;
	int plane_dim;
	int m_index, phys_index, renum_index;
	WORK_SET ws;
	ELEM_MATRIX tmp_elem_matrix;   /* each elem_matrix */

	RC_TRY( allocate_work_set(&ws) );

	/* elem_matrix->matrix[][] の確保 */
	/* preparation of total elem_matrix */
	elem_matrix->dim = 6;
	elem_matrix->element = element.label;
	elem_matrix->size = elem_matrix->dim * element.node_num;
	RC_TRY( allocate2D(elem_matrix->size, elem_matrix->size,
	                   &(elem_matrix->matrix)) );
	/* elem_matirx->index[] の確保と代入 */
	RC_TRY( allocate1D_i((elem_matrix->size), &(elem_matrix->index)) );
	for(ii1=0; ii1<element.node_num; ii1++){
		renum_index = node_renum_index(node, element.node[ii1]);
		RC_NEG_CHK(renum_index);
		for(ii2=0; ii2<(elem_matrix->dim); ii2++){
			elem_matrix->index[ii1*(elem_matrix->dim)+ii2]
			           = renum_index*(elem_matrix->dim) + ii2;
		}
	}

	phys_index = search_physical_prop_label(physical, element.physical);
	RC_NEG_CHK(phys_index);

	/* 膜剛性 */
	m_index = search_material_prop_element(material, physical, &(element) );
	RC_NEG_CHK(m_index);
	material.array[m_index].ana_type = ANA_PLANE_STRESS;
	plane_dim = 2;
	RC_TRY( fem_make_shell_elem_matrix(element, node, material,
	                                    physical, &tmp_elem_matrix, ws) );
	RC_TRY( fem_impose_shell_elem_matrix(plane_dim, element, material,
	                                     physical, tmp_elem_matrix,
	                                     elem_matrix) );

	RC_TRY( free_elem_matrix(&tmp_elem_matrix) );

	/* 曲げ剛性 */
	m_index = search_material_prop_element(material, physical, &(element) );
	RC_NEG_CHK(m_index);
	material.array[m_index].ana_type = ANA_PLANE_BEND;
	plane_dim = 3;
	RC_TRY( fem_make_shell_elem_matrix(element, node, material,
	                                    physical, &tmp_elem_matrix, ws) );
	RC_TRY( fem_impose_shell_elem_matrix(plane_dim, element, material,
	                                     physical, tmp_elem_matrix,
	                                     elem_matrix) );
	RC_TRY( free_elem_matrix(&tmp_elem_matrix) );

	/* せん断剛性 */
	m_index = search_material_prop_element(material, physical, &(element) );
	RC_NEG_CHK(m_index);
	material.array[m_index].ana_type = ANA_PLANE_SHEAR;
	plane_dim = 3; 
	RC_TRY( fem_make_shell_elem_matrix(element, node, material,
	                                    physical, &tmp_elem_matrix, ws) );
	RC_TRY( fem_impose_shell_elem_matrix(plane_dim, element, material,
	                                     physical, tmp_elem_matrix,
	                                     elem_matrix) );
	RC_TRY( free_elem_matrix(&tmp_elem_matrix) );

	/* ドリリング回転剛性 */
	RC_TRY( element_volume(&element, node) );
	RC_TRY( set_drilling(elem_matrix->dim, element, material, physical,
	                     elem_matrix) );

	/* 要素剛性マトリクスの座標変換 */
	RC_TRY( trans_coord_elem_matrix(node, element, ws, elem_matrix) );

	RC_TRY( free_work_set(&ws) );

	return(NORMAL_RC);
}


/* 種々の要素剛性マトリックスの作成（シェル解析用） */
RC fem_make_shell_elem_matrix(ELEMENT element, NODE_ARRAY node,
                              MATERIAL_PROP_ARRAY material,
                              PHYSICAL_PROP_ARRAY physical,
                              ELEM_MATRIX *elem_matrix,
                              WORK_SET ws)
{
	int ii1, ii2, ii3;
	int num_g_point;
	int ssnum;
	int m_index;
	double thickness;
	double det_J;
	double *weights = ws.vec_G[0];
	double **node_location = ws.mat_N_3[0];
	double **BtDB = ws.mat_3N_3N[0];
	double **D_matrix = ws.mat_6_6[0];
	double **g_points = ws.mat_G_4[0];
	double **B_matrix = ws.mat_6_3N[0];
	WORK_SET ws_tmp;

	if(element.label < 0) return(ARG_ERROR_RC);

	RC_TRY( allocate_work_set(&ws_tmp) );

	RC_NEG_CHK( ssnum = plane_element_ssnum(element, material, physical) );

	elem_matrix->dim = plane_element_dim(element, material, physical);
	RC_NEG_CHK(elem_matrix->dim);
	elem_matrix->element = element.label;
	elem_matrix->size = (element.node_num) * (elem_matrix->dim);

	/* elem_matrix->array[][] の確保 */
	RC_TRY( allocate2D(elem_matrix->size,
	                   elem_matrix->size, &(elem_matrix->matrix)) );
	/* elem_matirx->index[] の確保と代入 */
	RC_TRY( allocate1D_i((elem_matrix->size), &(elem_matrix->index)) );

	/* total_elem_matrixに代入するためのインデックス */
	for(ii1=0; ii1<element.node_num; ii1++){
		for(ii2=0; ii2<(elem_matrix->dim); ii2++){
			elem_matrix->index[ii1*(elem_matrix->dim)+ii2]
			             = ii1*(elem_matrix->dim) + ii2;
		}
	}

	/* D マトリックス */
	RC_TRY( get_D_matrix_shell(element, material, physical, D_matrix) );
	RC_TRY( get_thickness(element, physical, &thickness) );

	/* （局所座標系での）節点座標値のマトリックスを作成 */
	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( trans_coord_node_location(node, element, ws_tmp, node_location));

	/* ガウスポイントの設定 */
	if(REDUCED_INTEGRAL == 0){   /* 完全積分 */
		RC_TRY( Gauss_const_default(element.type, g_points, weights,
		                            &num_g_point));
	}else if(REDUCED_INTEGRAL == 1){   /* 次数低減積分 */
		RC_TRY( set_Gauss_const4RI(element.type, g_points, weights,
		                           &num_g_point) );
	}else{   /* 選択型次数低減積分(せん断項のみを次数低減積分) */
		RC_NEG_CHK( m_index = search_material_prop_element(material,
		                            physical, &(element)) );
		if(material.array[m_index].ana_type == ANA_PLANE_SHEAR){
			RC_TRY( set_Gauss_const4RI(element.type, g_points, weights,
		                               &num_g_point) );
		}else{
			RC_TRY( Gauss_const_default(element.type, g_points, weights,
		                                &num_g_point));
		}
	}

	for(ii1=0; ii1<num_g_point; ii1++){
		/* B マトリクス */
		RC_TRY( make_B_matrix_shell(element, material, physical, ws_tmp,
		                            g_points[ii1], node_location, B_matrix,
		                            &det_J) );
		/* B^T * D * B => BtDB */
		mul_matrix_AtBA(ssnum, (elem_matrix->size), B_matrix,
		                D_matrix, BtDB);
		/* BDB * detJ * weight * thickness */
		for(ii2=0; ii2<(elem_matrix->size); ii2++){
			for(ii3=0; ii3<(elem_matrix->size); ii3++){
				elem_matrix->matrix[ii2][ii3]
			         += BtDB[ii2][ii3] * det_J * weights[ii1] * thickness;
			}
		}
	}
	RC_TRY( free_work_set(&ws_tmp) );

	return(NORMAL_RC);
}


/* D マトリックス（シェル解析用） */
RC get_D_matrix_shell(ELEMENT element, MATERIAL_PROP_ARRAY material,
                      PHYSICAL_PROP_ARRAY physical, double **D_matrix)
{
	int ii1, ii2;
	int dim;
	int m_index;
	int ssnum;
	MATERIAL_PROP tmp_mat;

	m_index = search_material_prop_element(material, physical, &element);
	RC_NEG_CHK(m_index);

	if((dim = plane_element_dim(element, material, physical)) <= 0){
		fprintf(stderr, "plane_element_dim < ");
		return(ARG_ERROR_RC);
	}

	if((ssnum = plane_element_ssnum(element, material, physical)) <= 0){
		fprintf(stderr, "plane_element_ssnum < ");
		return(ARG_ERROR_RC);
	}

	tmp_mat = material.array[m_index];
	RC_TRY( fem_fill_material(&tmp_mat) );
	RC_TRY( fem_fill_physical(&tmp_mat, element, physical) )

	/* tmp_mat.D_matrix -> D_matrix */
	for(ii1=0;ii1<ssnum;ii1++){
		for(ii2=0;ii2<ssnum;ii2++){
			D_matrix[ii1][ii2] = tmp_mat.D_matrix[ii1][ii2];
		}
	}

	return(NORMAL_RC);
}


/* B マトリクスと |J|（シェル解析用） */
RC make_B_matrix_shell(ELEMENT element, MATERIAL_PROP_ARRAY material,
                       PHYSICAL_PROP_ARRAY physical, WORK_SET ws,
                       const double *g_point, double **node_location,
                       double **B_matrix, double *det_J)
{
	int elem_dim, plane_dim;
	int m_index;
	double *N = ws.vec_N[0];
	double **dN = ws.mat_3_N[0];
	double **dNxyz = ws.mat_3_N[1];
	double **Jacobi = ws.mat_3_3[0];
	double **inverse_J = ws.mat_3_3[1];

	elem_dim = element_dim(element.type);
	if(elem_dim != 2) return(ARG_ERROR_RC);
	m_index = search_material_prop_element(material, physical, &element);
	RC_NEG_CHK(m_index);

	/* g_point ->  N */
	/* g_point -> dN */
	/* dN * node_location => Jacobi */
	/* Jacobi^-1 => inverse_J */
	/* inverse_J * dN => dNxy */
	RC_TRY( set_N_vector(element.type, g_point, N) );
	RC_TRY( set_dN_matrix(element.type, g_point, dN) );
	mul_matrix(elem_dim, element.node_num, elem_dim, dN, node_location, Jacobi);
	RC_TRY( inverse_matrix(elem_dim, Jacobi, inverse_J) );
	mul_matrix(elem_dim, elem_dim, element.node_num, inverse_J, dN, dNxyz);

	/* B マトリクス */
	plane_dim = plane_element_dim(element, material, physical);
	if( (plane_dim != 2)&&(plane_dim != 3) ) return(ARG_ERROR_RC);
	RC_TRY( set_B_matrix_shell(plane_dim, element.node_num,
	                           material.array[m_index].ana_type,  N, dNxyz,
	                           B_matrix) );

	if(det_J != NULL){
		/* 注意 : このヤコビアンは局所座標によるもの */
		*det_J = determinant(elem_dim, Jacobi);
	}

	return(NORMAL_RC);
}


/* N, dNxyz を用いて B マトリクスを作成(シェル解析用) */
/* dNxyz[0][i]...Ni の x 方向微分 */
/* dNxyz[1][i]...Ni の y 方向微分 */
/* dNxyz[2][i]...Ni の z 方向微分 */
/* dim=2 ...2次元, dim=3 ...3次元 */
static RC set_B_matrix_shell(int dim, int node_num,
                             ANALYSIS_TYPE ana_type, double *N,
                             double **dNxyz, double **B_matrix)
{
	int ii1;

	switch(ana_type){
	case ANA_PLANE_STRESS:
	case ANA_PLANE_STRAIN:
		for(ii1=0;ii1<node_num;ii1++){
			B_matrix[0][dim*ii1]   = dNxyz[0][ii1];
			B_matrix[0][dim*ii1+1] = 0.0;
			B_matrix[1][dim*ii1]   = 0.0;
			B_matrix[1][dim*ii1+1] = dNxyz[1][ii1];
			B_matrix[2][dim*ii1]   = dNxyz[1][ii1];
			B_matrix[2][dim*ii1+1] = dNxyz[0][ii1];
		}
		break;
	case ANA_PLANE_BEND:
		for(ii1=0;ii1<node_num;ii1++){
			B_matrix[0][dim*ii1]   = 0.0;
			B_matrix[0][dim*ii1+1] = 0.0;
			B_matrix[0][dim*ii1+2] = -dNxyz[0][ii1];
			B_matrix[1][dim*ii1]   = 0.0;
			B_matrix[1][dim*ii1+1] = dNxyz[1][ii1];
			B_matrix[1][dim*ii1+2] = 0.0;
			B_matrix[2][dim*ii1]   = 0.0;
			B_matrix[2][dim*ii1+1] = dNxyz[0][ii1];
			B_matrix[2][dim*ii1+2] = -dNxyz[1][ii1];
		}
		break;
	case ANA_PLANE_SHEAR:
		for(ii1=0;ii1<node_num;ii1++){
			B_matrix[0][dim*ii1]   = dNxyz[0][ii1];
			B_matrix[0][dim*ii1+1] = 0.0;
			B_matrix[0][dim*ii1+2] = N[ii1];
			B_matrix[1][dim*ii1]   = dNxyz[1][ii1];
			B_matrix[1][dim*ii1+1] = -N[ii1];
			B_matrix[1][dim*ii1+2] = 0.0;
		}
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


/* 要素剛性マトリクスの座標変換 */
/* [K_GLOBAL] = [T]^T * [K_LOCAL] * [T] */
/* (スパース性を考慮した行列積を利用) */
static RC trans_coord_elem_matrix(NODE_ARRAY node, ELEMENT element,
                                  WORK_SET ws, ELEM_MATRIX *elem_matrix)
{
	int ix1, ix2;
	int ii1, ii2, ii3;
	double **dcosine_matrix = ws.mat_3_3[0];
	double **sub_matrix = ws.mat_3_3[1];
	TRANSLATION3D v_dcosine[3];

	RC_TRY( set_elem_dcosine(node, element, v_dcosine) );
	vect2dcosine_array(v_dcosine, dcosine_matrix);

	for(ix1=0; ix1<elem_matrix->size; ix1+=3){
		for(ix2=0; ix2<elem_matrix->size; ix2+=3){
			for(ii1=0; ii1<3; ii1++){
				for(ii2=0; ii2<3; ii2++){
					sub_matrix[ii1][ii2]
					   = elem_matrix->matrix[ix1+ii1][ix2+ii2];
					elem_matrix->matrix[ix1+ii1][ix2+ii2] = 0.0;
				}
			}
			for(ii1=0; ii1<3; ii1++){
				for(ii2=0; ii2<3; ii2++){
					double sum = 0.0;
					double *A_ii2 = dcosine_matrix[ii2];

					for(ii3=0; ii3<3; ii3++){
						sum += dcosine_matrix[ii3][ii1] * sub_matrix[ii3][ii2];
					}
					for(ii3=0; ii3<3; ii3++){
						elem_matrix->matrix[ix1+ii1][ix2+ii3]
						    += sum * A_ii2[ii3];
					}
				}
			}
		}
	}

	return(NORMAL_RC);
}


static void vect2dcosine_array(TRANSLATION3D v[], double **array)
{
	array[0][0] = v[0].x;
	array[0][1] = v[0].y;
	array[0][2] = v[0].z;
	array[1][0] = v[1].x;
	array[1][1] = v[1].y;
	array[1][2] = v[1].z;
	array[2][0] = v[2].x;
	array[2][1] = v[2].y;
	array[2][2] = v[2].z;
}


/* （局所座標系での）節点座標値のマトリックスを作成 */
static RC trans_coord_node_location(NODE_ARRAY node, ELEMENT element,
                                    WORK_SET ws, double **coord)
{
	int ii1;
	double *N = ws.vec_N[0];
	double **local_points = ws.mat_N_4[1];
	double **dcosine_matrix = ws.mat_3_3[1];
	double **vtmp = ws.mat_N_3[1];
	TRANSLATION3D center_point;
	TRANSLATION3D v_dcosine[3];

	switch(element.type){
		case ELEM_TRI1:
		case ELEM_TRI2:
			vtmp[0][0] = 0.0;
			vtmp[1][0] = 0.0;
			vtmp[2][0] = 0.0;
			for(ii1=1; ii1<element.node_num; ii1++){
				vtmp[ii1][0] = coord[ii1][0] - coord[0][0];
				vtmp[ii1][1] = coord[ii1][1] - coord[0][1];
				vtmp[ii1][2] = coord[ii1][2] - coord[0][2];
			}
			break;
		case ELEM_QUAD1:
		case ELEM_QUAD2:
			RC_TRY( set_center_point(element.type, local_points[0]) );
			RC_TRY( set_N_vector(element.type, local_points[0], N) );
			center_point.x = center_point.y = center_point.z = 0.0;
			for(ii1=0; ii1<element.node_num; ii1++){
				center_point.x += ( coord[ii1][0] * N[ii1] );
				center_point.y += ( coord[ii1][1] * N[ii1] );
				center_point.z += ( coord[ii1][2] * N[ii1] );
			}
			for(ii1=0; ii1<element.node_num; ii1++){
				vtmp[ii1][0] = coord[ii1][0] - center_point.x;
				vtmp[ii1][1] = coord[ii1][1] - center_point.y;
				vtmp[ii1][2] = coord[ii1][2] - center_point.z;
			}
			break;
		default:
			return(IMPLEMENT_ERROR_RC);
	}

	/* 方向余弦マトリクスの設定 */
	RC_TRY( set_elem_dcosine(node, element, v_dcosine) );
	vect2dcosine_array(v_dcosine, dcosine_matrix);

	transpose_overwrite_matrix(3, dcosine_matrix);
	/* 座標変換 */
	mul_matrix(element.node_num, 3, 3, vtmp, dcosine_matrix, coord);

	return(NORMAL_RC);
}


/* ドリリング回転剛性 (面の法線回りの回転剛性) */
static RC set_drilling(int dim, ELEMENT element,
                       MATERIAL_PROP_ARRAY material,
                       PHYSICAL_PROP_ARRAY physical,
                       ELEM_MATRIX *elem_matrix)
{
	int ii1, ii2;
	int ix1, ix2;
	int m_index;
	double coeff1, coeff2; /* 係数 */
	double thickness;

	/* 係数の計算 */
	m_index = search_material_prop_element(material, physical, &element);
	RC_NEG_CHK(m_index);
	RC_TRY( get_thickness(element, physical, &thickness) );
	coeff1 = DRILLING * material.array[m_index].E * thickness * element.volume;
	coeff2 = - coeff1 * ( 1.0/(double)(element.node_num-1) );

	/* 対応する行と列に組み込む */
	for(ii1=0; ii1<element.node_num; ii1++){
		ix1 = dim*(ii1+1) - 1;
		for(ii2=0; ii2<element.node_num; ii2++){
			ix2 = dim*(ii2+1) - 1;
			elem_matrix->matrix[ix1][ix2] = coeff2;
		}
	}
	/* 対角要素に代入 */
	for(ii1=0; ii1<element.node_num; ii1++){
		ix1 =  dim*(ii1+1) - 1;
		elem_matrix->matrix[ix1][ix1] = coeff1;
	}

	return(NORMAL_RC);
}


/* 次数低減積分用に積分点を設定 */
static RC set_Gauss_const4RI(ELEM_TYPE type, double **points,
                             double *weight, int *num_point)
{
	switch(type){
	case ELEM_TRI1:
		*num_point = 1;
		RC_TRY( Gauss_const_tri(*num_point, points, weight) );
		break;
	case ELEM_TRI2:
		*num_point = 3;
		RC_TRY( Gauss_const_tri(*num_point, points, weight) );
		break;
	case ELEM_QUAD1:
		*num_point = 1;
		RC_TRY( Gauss_const_quad(*num_point, points, weight) );
		break;
	case ELEM_QUAD2:
		*num_point = 4;
		RC_TRY( Gauss_const_quad(*num_point, points, weight) );
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);

}


/* 応力、歪の計算 */
/* stress, strain が NULL なら代入しない */
/* 応力、歪を計算する位置（層）をdistance[] で指定 */
/* distanceがNULLなら、-thickness/2とthickness/2 (裏と表)に設定 */
RC fem_stress_strain_shell(int fem_sol_ss_flag, double *distance,
                           NODE_ARRAY node, ELEMENT_ARRAY element,
                           DISP_ARRAY disp,
                           MATERIAL_PROP_ARRAY material,
                           PHYSICAL_PROP_ARRAY physical,
                           STRESS_ARRAY stress[], STRAIN_ARRAY strain[])
{
	int ii1, ii2;
	int node_index;
	int *count = NULL;
	double thickness;
	double default_distance[2];
	double **local_points = NULL;
	TRANS_ROTATION3D *stress_sum[2];
	TRANS_ROTATION3D *strain_sum[2];
	TRANS_ROTATION3D init_v = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	STRESS_STRAIN local_stress[2], local_strain[2];
	WORK_SET ws;

	if( ((stress == NULL)&&(strain == NULL)) || (fem_sol_ss_flag == 0) ){
		return(NORMAL_RC);
	}

	if(fem_sol_ss_flag & SOL_SS_NODE_AV){
		RC_NULL_CHK( stress_sum[0] = (TRANS_ROTATION3D *)
                              malloc((node.size)*sizeof(TRANS_ROTATION3D)) );
		RC_NULL_CHK( stress_sum[1] = (TRANS_ROTATION3D *)
                              malloc((node.size)*sizeof(TRANS_ROTATION3D)) );
		RC_NULL_CHK( strain_sum[0] = (TRANS_ROTATION3D *)
                              malloc((node.size)*sizeof(TRANS_ROTATION3D)) );
		RC_NULL_CHK( strain_sum[1] = (TRANS_ROTATION3D *)
                              malloc((node.size)*sizeof(TRANS_ROTATION3D)) );
		RC_NULL_CHK( count = (int *) malloc((node.size)*sizeof(int)) );
		for(ii1=0; ii1<node.size; ii1++){
			count[ii1] = 0;
			stress_sum[0][ii1] = stress_sum[1][ii1] = init_v;
			strain_sum[0][ii1] = strain_sum[1][ii1] = init_v;
		}
	}

	if(stress != NULL){
		RC_TRY( allocate_stress_array(0, &(stress[0])) );
		RC_TRY( allocate_stress_array(0, &(stress[1])) );
		stress[0].source = stress[1].source = K_SOL;
	}
	if(strain != NULL){
		RC_TRY( allocate_strain_array(0, &(strain[0])) );
		RC_TRY( allocate_strain_array(0, &(strain[1])) );
		strain[0].source = strain[1].source = K_SOL;
	}

	RC_TRY( allocate_work_set(&ws) );
	local_points = ws.mat_N_4[0];

	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		if(distance == NULL){ /* 応力、歪を計算する位置を板の裏と表に設定 */
			RC_TRY( get_thickness(element.array[ii1], physical, &thickness) );
			default_distance[0] = -0.5*thickness;
			default_distance[1] =  0.5*thickness;
		}else{ /* 応力、歪を計算する位置(層)をユーザ指定値に設定 */
			default_distance[0] = distance[0];
			default_distance[1] = distance[1];
		}

		/* 要素中心 */
		if(fem_sol_ss_flag & SOL_SS_CENTER){
			RC_TRY( set_center_point(element.array[ii1].type,
                                     local_points[0]) );
			RC_TRY( fem_local_stress_strain_shell(default_distance,
			           element.array[ii1], local_points[0], node, disp,
			           material, physical, ws, local_stress, local_strain) );
			if(stress != NULL){
				RC_TRY( realloc_stress_array(&(stress[0])) );
				RC_TRY( realloc_stress_array(&(stress[1])) );
				stress[0].array[stress[0].size] = local_stress[0];
				stress[1].array[stress[1].size] = local_stress[1];
				stress[0].array[stress[0].size].node = -1;  /* center */
				stress[1].array[stress[1].size].node = -1;  /* center */
				(stress[0].size)++;
				(stress[1].size)++;
			}
			if(strain != NULL){
				RC_TRY( realloc_strain_array(&(strain[0])) );
				RC_TRY( realloc_strain_array(&(strain[1])) );
				strain[0].array[strain[0].size] = local_strain[0];
				strain[1].array[strain[1].size] = local_strain[1];
				strain[0].array[strain[0].size].node = -1;  /* center */
				strain[1].array[strain[1].size].node = -1;  /* center */
				(strain[0].size)++;
				(strain[1].size)++;
			}
		}

		/* 各節点 */
		if( (fem_sol_ss_flag & SOL_SS_NODE)
          ||(fem_sol_ss_flag & SOL_SS_NODE_AV) ){
			RC_TRY( set_local_node_points(element.array[ii1].type,
			                              local_points) );

			for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
				RC_TRY( fem_local_stress_strain_shell(default_distance,
				            element.array[ii1], local_points[ii2], node, disp,
				            material, physical, ws, local_stress,
				            local_strain) );

				if(fem_sol_ss_flag & SOL_SS_NODE_AV){
					node_index = search_node_label(node,
					                   element.array[ii1].node[ii2]);
					RC_NEG_CHK( node_index );
					count[node_index]++;
					stress_sum[0][node_index]
					    = add_trans_rotation3d(stress_sum[0][node_index],
					                           local_stress[0].v);
					stress_sum[1][node_index]
					    = add_trans_rotation3d(stress_sum[1][node_index],
					                           local_stress[1].v);
					strain_sum[0][node_index]
					    = add_trans_rotation3d(strain_sum[0][node_index],
					                           local_strain[0].v);
					strain_sum[1][node_index]
					    = add_trans_rotation3d(strain_sum[1][node_index],
					                           local_strain[1].v);
				}

				if(stress != NULL){
					RC_TRY( realloc_stress_array(&(stress[0])) );
					RC_TRY( realloc_stress_array(&(stress[1])) );
					stress[0].array[stress[0].size] = local_stress[0];
					stress[1].array[stress[1].size] = local_stress[1];
					stress[0].array[stress[0].size].node
					      = element.array[ii1].node[ii2];
					stress[1].array[stress[1].size].node
					      = element.array[ii1].node[ii2];
					(stress[0].size)++;
					(stress[1].size)++;
				}

				if(strain != NULL){
					RC_TRY( realloc_strain_array(&(strain[0])) );
					RC_TRY( realloc_strain_array(&(strain[1])) );
					strain[0].array[strain[0].size] = local_strain[0];
					strain[1].array[strain[1].size] = local_strain[1];
					strain[0].array[strain[0].size].node
					      = element.array[ii1].node[ii2];
					strain[1].array[strain[1].size].node
					      = element.array[ii1].node[ii2];
					(strain[0].size)++;
					(strain[1].size)++;
				}
			}
		}
	}

	if(stress != NULL){
		RC_TRY( clean_stress_array(&(stress[0])) );
		RC_TRY( clean_stress_array(&(stress[1])) );
	}
	if(strain != NULL){
		RC_TRY( clean_strain_array(&(strain[0])) );
		RC_TRY( clean_strain_array(&(strain[1])) );
	}

	if(fem_sol_ss_flag & SOL_SS_NODE_AV){
		if(stress != NULL){
			for(ii1=0; ii1<stress[0].size; ii1++){
				if(stress[0].array[ii1].element < 0) continue;
				node_index = search_node_label(node,
				                                   stress[0].array[ii1].node);
				if(count[node_index] <= 0) return(UNKNOWN_ERROR_RC);
				stress[0].array[ii1].v
				 = mul_scalar_trans_rotation3d(1.0/(double)(count[node_index]),
				                               stress_sum[0][node_index]);
				stress[1].array[ii1].v
				 = mul_scalar_trans_rotation3d(1.0/(double)(count[node_index]),
				                               stress_sum[1][node_index]);
			}
		}
		if(strain != NULL){
			for(ii1=0; ii1<strain->size; ii1++){
				if(strain->array[ii1].element < 0) continue;
				node_index = search_node_label(node,
				                                   strain->array[ii1].node);
				if(count[node_index] <= 0) return(UNKNOWN_ERROR_RC);
				strain[0].array[ii1].v
				 = mul_scalar_trans_rotation3d(1.0/(double)(count[node_index]),
				                               strain_sum[0][node_index]);
				strain[1].array[ii1].v
				 = mul_scalar_trans_rotation3d(1.0/(double)(count[node_index]),
				                               strain_sum[1][node_index]);
			}
		}
		free(stress_sum[0]);
		free(stress_sum[1]);
		free(strain_sum[0]);
		free(strain_sum[1]);
		free(count);
	}

	RC_TRY( free_work_set(&ws) );

	return(NORMAL_RC);
}


/* シェル要素の局所応力、ひずみ */
RC fem_local_stress_strain_shell(double distance[], ELEMENT element,
                                 const double *local_point,
                                 NODE_ARRAY node, DISP_ARRAY disp,
                                 MATERIAL_PROP_ARRAY material,
                                 PHYSICAL_PROP_ARRAY physical,
                                 WORK_SET ws, STRESS_STRAIN stress[],
                                 STRESS_STRAIN strain[])
{
	int ii1;
	int dim, m_index, phys_index;
	int ssnum;
	double thickness, default_distance[2];
	double disp_v[3*MAX_NODE];
	double strain_v[2][6];
	double stress_v[2][6];
	double **node_location = ws.mat_N_3[0];
	double **B_matrix = ws.mat_6_3N[0];
	double **D_matrix = ws.mat_6_6[0];
	WORK_SET ws_tmp;

	if(element.label < 0) return(ARG_ERROR_RC);

	m_index = search_material_prop_element(material, physical, &(element));
	RC_NEG_CHK(m_index);

	phys_index = search_physical_prop_label(physical, element.physical);
	RC_NEG_CHK(phys_index)

	RC_TRY( allocate_work_set(&ws_tmp) );

	for(ii1=0; ii1<2; ii1++){
		init_stress_strain(&(stress[ii1]));
		init_stress_strain(&(strain[ii1]));
		stress[ii1].element = strain[ii1].element = element.label;
	}

	if(distance == NULL){ /* 応力、歪を計算する位置を板の裏と表に設定 */
		RC_TRY( get_thickness(element, physical, &thickness) );
		default_distance[0] = -0.5*thickness;
		default_distance[1] =  0.5*thickness;
	}else{ /* 応力、歪を計算する位置をユーザ指定値に設定 */
		default_distance[0] = distance[0];
		default_distance[1] = distance[1];
	}

	/* membrane */
	dim = 2;
	ssnum = 3;
	material.array[m_index].ana_type = ANA_PLANE_STRESS;

	/* B_matrix, disp_v */
	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( trans_coord_node_location(node, element, ws_tmp, node_location));
	RC_TRY( make_B_matrix_shell(element, material, physical, ws_tmp,
	                            local_point, node_location, B_matrix, NULL) );
	RC_TRY( fill_disp_vector4shell(dim, ANA_PLANE_STRESS, ws_tmp, node, element,
	                               disp, disp_v) );
	/* B_matrix, disp_v -> strain_v */
	mul_matrix_vect(ssnum, dim*element.node_num, B_matrix, disp_v, strain_v[0]);
	/* D_matrix */
	RC_TRY( get_D_matrix_shell(element, material, physical, D_matrix) );

	/* D_matrix, strain_v -> stress_v */
	mul_matrix_vect(ssnum, ssnum, D_matrix, strain_v[0], stress_v[0]);

	/* strain_v, stress_v -> strain, stress */
	for(ii1=0; ii1<2; ii1++){
		RC_TRY( vect2stress_strain_shell(ANA_PLANE_STRESS, strain_v[0],
		                                 &(strain[ii1])) );
		RC_TRY( vect2stress_strain_shell(ANA_PLANE_STRESS, stress_v[0],
		                                 &(stress[ii1])) );
	}

	/* bending */
	dim = 3;
	ssnum = 3;
	material.array[m_index].ana_type = ANA_PLANE_BEND;

	/* B_matrix, disp_v */
	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( trans_coord_node_location(node, element, ws_tmp, node_location));
	RC_TRY( make_B_matrix_shell(element, material, physical, ws_tmp,
	                            local_point, node_location, B_matrix, NULL) );
	RC_TRY( fill_disp_vector4shell(dim, ANA_PLANE_BEND, ws_tmp, node, element,
	                               disp, disp_v) );
	/* B_matrix, disp_v -> strain_v */
	mul_matrix_vect(ssnum, dim*element.node_num, B_matrix, disp_v, strain_v[0]);
	for(ii1=0; ii1<ssnum; ii1++){
		strain_v[1][ii1] = strain_v[0][ii1]*default_distance[1];
		strain_v[0][ii1] *= default_distance[0];
	}

	/* D_matrix */
	material.array[m_index].ana_type = ANA_PLANE_STRESS;
	RC_TRY( get_D_matrix_shell(element, material, physical, D_matrix) );

	/* D_matrix, strain_v -> stress_v */
	mul_matrix_vect(ssnum, ssnum, D_matrix, strain_v[0], stress_v[0]);
	mul_matrix_vect(ssnum, ssnum, D_matrix, strain_v[1], stress_v[1]);

	/* strain_v, stress_v -> strain, stress */
	for(ii1=0; ii1<2; ii1++){
		RC_TRY( vect2stress_strain_shell(ANA_PLANE_BEND, strain_v[ii1],
		                                 &(strain[ii1])) );
		RC_TRY( vect2stress_strain_shell(ANA_PLANE_BEND, stress_v[ii1],
		                                  &(stress[ii1])) );
	}

	/* shear */
	dim = 3;
	ssnum = 2;
	material.array[m_index].ana_type = ANA_PLANE_SHEAR;

	/* B_matrix, disp_v */
	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( trans_coord_node_location(node, element, ws_tmp, node_location));
	RC_TRY( make_B_matrix_shell(element, material, physical, ws_tmp,
	                            local_point, node_location, B_matrix, NULL) );
	RC_TRY( fill_disp_vector4shell(dim, ANA_PLANE_SHEAR, ws_tmp, node, element,
	                               disp, disp_v) );

	/* B_matrix, disp_v -> strain_v */
	mul_matrix_vect(ssnum, dim*element.node_num, B_matrix, disp_v,
	                strain_v[0]);
	/* D_matrix */
	RC_TRY( get_D_matrix_shell(element, material, physical, D_matrix) );

	/* D_matrix, strain_v -> stress_v */
	mul_matrix_vect(ssnum, ssnum, D_matrix, strain_v[0], stress_v[0]);

	/* strain_v, stress_v -> strain, stress */
	for(ii1=0; ii1<2; ii1++){
		RC_TRY( vect2stress_strain_shell(ANA_PLANE_SHEAR, strain_v[0],
		                                 &(strain[ii1])) );
		RC_TRY( vect2stress_strain_shell(ANA_PLANE_SHEAR, stress_v[0],
		                                 &(stress[ii1])) );
	}

	RC_TRY( free_work_set(&ws_tmp) );

	return(NORMAL_RC);
}


/* 合応力、曲率の計算 */
/*  force 〜 curv が NULL なら代入しない */
RC fem_resultant_stress_strain_shell(int fem_sol_ss_flag, NODE_ARRAY node,
                                     ELEMENT_ARRAY element,
                                     DISP_ARRAY disp,
                                     MATERIAL_PROP_ARRAY material,
                                     PHYSICAL_PROP_ARRAY physical,
                                     STRESS_ARRAY *force,
                                     STRESS_ARRAY *moment,
                                     STRESS_ARRAY *shear_force,
                                     STRAIN_ARRAY *curv)
{
	int ii1, ii2;
	int node_index;
	int *count = NULL;
	double **local_points = NULL;
	TRANS_ROTATION3D *force_sum = NULL, *moment_sum = NULL;
	TRANS_ROTATION3D *shear_force_sum = NULL, *curv_sum = NULL;
	TRANS_ROTATION3D init_v = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	STRESS_STRAIN local_force, local_moment, local_shear_force;
	STRESS_STRAIN local_curv;
	WORK_SET ws;

	if( ((force == NULL)&&(moment == NULL)&&(shear_force == NULL)
	     &&(curv == NULL)) || (fem_sol_ss_flag == 0) ){
		return(NORMAL_RC);
	}

	if(fem_sol_ss_flag & SOL_SS_NODE_AV){
		RC_NULL_CHK( force_sum = (TRANS_ROTATION3D *)
                              malloc((node.size)*sizeof(TRANS_ROTATION3D)) );
		RC_NULL_CHK( moment_sum = (TRANS_ROTATION3D *)
                              malloc((node.size)*sizeof(TRANS_ROTATION3D)) );
		RC_NULL_CHK( shear_force_sum = (TRANS_ROTATION3D *)
                              malloc((node.size)*sizeof(TRANS_ROTATION3D)) );
		RC_NULL_CHK( curv_sum = (TRANS_ROTATION3D *)
                              malloc((node.size)*sizeof(TRANS_ROTATION3D)) );
		RC_NULL_CHK( count = (int *) malloc((node.size)*sizeof(int)) );
		for(ii1=0; ii1<node.size; ii1++){
			count[ii1] = 0;
			force_sum[ii1] = moment_sum[ii1] = shear_force_sum[ii1]
			               = curv_sum[ii1] = init_v;
		}
	}

	if(force != NULL){
		RC_TRY( allocate_stress_array(0, force) );
		force->source = K_SOL;
	}
	if(moment != NULL){
		RC_TRY( allocate_stress_array(0, moment) );
		moment->source = K_SOL;
	}
	if(shear_force != NULL){
		RC_TRY( allocate_stress_array(0, shear_force) );
		shear_force->source = K_SOL;
	}
	if(curv != NULL){
		RC_TRY( allocate_strain_array(0, curv) );
		curv->source = K_SOL;
	}

	RC_TRY( allocate_work_set(&ws) );
	local_points = ws.mat_N_4[0];

	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		/* 要素中心 */
		if(fem_sol_ss_flag & SOL_SS_CENTER){
			RC_TRY( set_center_point(element.array[ii1].type,
                                     local_points[0]) );
			RC_TRY( fem_local_resultant_stress_strain(element.array[ii1],
			           local_points[0], node, disp, material, physical, ws,
			           &local_force, &local_moment, &local_shear_force,
			           &local_curv) );
			if(force != NULL){
				RC_TRY( realloc_stress_array(force) );
				force->array[force->size] = local_force;
				force->array[force->size].node = -1;  /* center */
				(force->size)++;
			}
			if(moment != NULL){
				RC_TRY( realloc_stress_array(moment) );
				moment->array[moment->size] = local_moment;
				moment->array[moment->size].node = -1;  /* center */
				(moment->size)++;
			}
			if(shear_force != NULL){
				RC_TRY( realloc_stress_array(shear_force) );
				shear_force->array[shear_force->size] = local_shear_force;
				shear_force->array[shear_force->size].node = -1;  /* center */
				(shear_force->size)++;
			}
			if(curv != NULL){
				RC_TRY( realloc_strain_array(curv) );
				curv->array[curv->size] = local_curv;
				curv->array[curv->size].node = -1;  /* center */
				(curv->size)++;
			}
		}

		/* 各節点 */
		if( (fem_sol_ss_flag & SOL_SS_NODE)
          ||(fem_sol_ss_flag & SOL_SS_NODE_AV) ){
			RC_TRY( set_local_node_points(element.array[ii1].type,
			                              local_points) );

			for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
				RC_TRY( fem_local_resultant_stress_strain(element.array[ii1],
				            local_points[ii2], node, disp, material, physical,
				            ws, &local_force, &local_moment,
				            &local_shear_force, &local_curv) );

				if(fem_sol_ss_flag & SOL_SS_NODE_AV){
					node_index = search_node_label(node,
					                   element.array[ii1].node[ii2]);
					RC_NEG_CHK( node_index );
					count[node_index]++;
					force_sum[node_index]
					    = add_trans_rotation3d(force_sum[node_index],
					                           local_force.v);
					moment_sum[node_index]
					    = add_trans_rotation3d(moment_sum[node_index],
					                           local_moment.v);
					shear_force_sum[node_index]
					    = add_trans_rotation3d(shear_force_sum[node_index],
					                           local_shear_force.v);
					curv_sum[node_index]
					    = add_trans_rotation3d(curv_sum[node_index],
					                           local_curv.v);
				}

				if(force != NULL){
					RC_TRY( realloc_stress_array(force) );
					force->array[force->size] = local_force;
					force->array[force->size].node 
					     = element.array[ii1].node[ii2];
					(force->size)++;
				}
				if(moment != NULL){
					RC_TRY( realloc_stress_array(moment) );
					moment->array[moment->size] = local_moment;
					moment->array[moment->size].node
					      = element.array[ii1].node[ii2];
					(moment->size)++;
				}
				if(shear_force != NULL){
					RC_TRY( realloc_stress_array(shear_force) );
					shear_force->array[shear_force->size] = local_shear_force;
					shear_force->array[shear_force->size].node
					      = element.array[ii1].node[ii2];
					(shear_force->size)++;
				}
				if(curv != NULL){
					RC_TRY( realloc_strain_array(curv) );
					curv->array[curv->size] = local_curv;
					curv->array[curv->size].node
					      = element.array[ii1].node[ii2];
					(curv->size)++;
				}
			}
		}
	}

	if(force != NULL) RC_TRY( clean_stress_array(force) );
	if(moment != NULL) RC_TRY( clean_stress_array(moment) );
	if(shear_force != NULL) RC_TRY( clean_stress_array(shear_force) );
	if(curv != NULL) RC_TRY( clean_strain_array(curv) );

	if(fem_sol_ss_flag & SOL_SS_NODE_AV){
		if(force != NULL){
			for(ii1=0; ii1<force->size; ii1++){
				if(force->array[ii1].element < 0) continue;
				node_index = search_node_label(node,
				                                   force->array[ii1].node);
				if(count[node_index] <= 0) return(UNKNOWN_ERROR_RC);
				force->array[ii1].v
				 = mul_scalar_trans_rotation3d(1.0/(double)(count[node_index]),
				                               force_sum[node_index]);
			}
		}
		if(moment != NULL){
			for(ii1=0; ii1<moment->size; ii1++){
				if(moment->array[ii1].element < 0) continue;
				node_index = search_node_label(node,
				                                   moment->array[ii1].node);
				if(count[node_index] <= 0) return(UNKNOWN_ERROR_RC);
				moment->array[ii1].v
				 = mul_scalar_trans_rotation3d(1.0/(double)(count[node_index]),
				                               moment_sum[node_index]);
			}
		}
		if(shear_force != NULL){
			for(ii1=0; ii1<shear_force->size; ii1++){
				if(shear_force->array[ii1].element < 0) continue;
				node_index = search_node_label(node,
				                                  shear_force->array[ii1].node);
				if(count[node_index] <= 0) return(UNKNOWN_ERROR_RC);
				shear_force->array[ii1].v
				 = mul_scalar_trans_rotation3d(1.0/(double)(count[node_index]),
				                               shear_force_sum[node_index]);
			}
		}
		if(curv != NULL){
			for(ii1=0; ii1<curv->size; ii1++){
				if(curv->array[ii1].element < 0) continue;
				node_index = search_node_label(node,
				                                   curv->array[ii1].node);
				if(count[node_index] <= 0) return(UNKNOWN_ERROR_RC);
				curv->array[ii1].v
				 = mul_scalar_trans_rotation3d(1.0/(double)(count[node_index]),
				                               curv_sum[node_index]);
			}
		}
		free(force_sum);
		free(moment_sum);
		free(shear_force_sum);
		free(curv_sum);
		free(count);
	}

	RC_TRY( free_work_set(&ws) );

	return(NORMAL_RC);
}


/* シェル要素の局所合応力、曲率 */
RC fem_local_resultant_stress_strain(ELEMENT element,
                                     const double *local_point,
                                     NODE_ARRAY node, DISP_ARRAY disp,
                                     MATERIAL_PROP_ARRAY material,
                                     PHYSICAL_PROP_ARRAY physical,
                                     WORK_SET ws, STRESS_STRAIN *force,
                                     STRESS_STRAIN *moment,
                                     STRESS_STRAIN *shear_force,
                                     STRESS_STRAIN *curv)
{
	int dim, m_index, phys_index,  ssnum;
	double disp_v[3*MAX_NODE];
	double strain_v[6];
	double stress_v[6];
	double **node_location = ws.mat_N_3[0];
	double **B_matrix = ws.mat_6_3N[0];
	double **D_matrix = ws.mat_6_6[0];
	WORK_SET ws_tmp;

	if(element.label < 0) return(ARG_ERROR_RC);

	m_index = search_material_prop_element(material, physical, &(element));
	RC_NEG_CHK(m_index);

	phys_index = search_physical_prop_label(physical, element.physical);
	RC_NEG_CHK(phys_index)

	RC_TRY( allocate_work_set(&ws_tmp) );

	init_stress_strain(force);
	init_stress_strain(moment);
	init_stress_strain(shear_force);
	force->element = moment->element = shear_force->element = element.label;

	init_stress_strain(curv);
	curv->element = element.label;

	/* membrane */
	dim = 2;
	ssnum = 3;
	material.array[m_index].ana_type = ANA_PLANE_STRESS;

	/* B_matrix, disp_v */
	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( trans_coord_node_location(node, element, ws, node_location));
	RC_TRY( make_B_matrix_shell(element, material, physical, ws_tmp,
	                            local_point, node_location, B_matrix, NULL) );
	RC_TRY( fill_disp_vector4shell(dim, ANA_PLANE_STRESS, ws, node, element,
	                               disp, disp_v) );
	/* B_matrix, disp_v -> strain_v */
	mul_matrix_vect(ssnum, dim*element.node_num, B_matrix, disp_v,
	                strain_v);
	/* D_matrix */
	RC_TRY( get_D_matrix_shell(element, material, physical, D_matrix) );

	/* D_matrix, strain_v -> stress_v */
	mul_matrix_vect(ssnum, ssnum, D_matrix, strain_v, stress_v);

	/* strain_v, stress_v -> strain, stress */
	RC_TRY( vect2resultant_stress_strain(ANA_PLANE_STRESS, stress_v, force) );

	/* bending */
	dim = 3;
	ssnum = 3;
	material.array[m_index].ana_type = ANA_PLANE_BEND;

	/* B_matrix, disp_v */
	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( trans_coord_node_location(node, element, ws, node_location));
	RC_TRY( make_B_matrix_shell(element, material, physical, ws_tmp, 
	                            local_point, node_location, B_matrix, NULL) );
	RC_TRY( fill_disp_vector4shell(dim, ANA_PLANE_BEND, ws, node, element,
	                               disp, disp_v) );
	/* B_matrix, disp_v -> strain_v */
	mul_matrix_vect(ssnum, dim*element.node_num, B_matrix, disp_v, strain_v);

	/* D_matrix */
	RC_TRY( get_D_matrix_shell(element, material, physical, D_matrix) );

	/* D_matrix, strain_v -> stress_v */
	mul_matrix_vect(ssnum, ssnum, D_matrix, strain_v, stress_v);

	/* strain_v, stress_v -> strain, stress */
	RC_TRY( vect2resultant_stress_strain(ANA_PLANE_BEND, stress_v, moment) );
	RC_TRY( vect2resultant_stress_strain(ANA_PLANE_BEND, strain_v, curv) );

	/* shear */
	dim = 3;
	ssnum = 2;
	material.array[m_index].ana_type = ANA_PLANE_SHEAR;

	/* B_matrix, disp_v */
	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( trans_coord_node_location(node, element, ws, node_location));
	RC_TRY( make_B_matrix_shell(element, material, physical, ws_tmp,
	                            local_point, node_location, B_matrix, NULL) );
	RC_TRY( fill_disp_vector4shell(dim, ANA_PLANE_SHEAR, ws, node, element,
	                               disp, disp_v) );

	/* B_matrix, disp_v -> strain_v */
	mul_matrix_vect(ssnum, dim*element.node_num, B_matrix, disp_v,
	                strain_v);
	/* D_matrix */
	RC_TRY( get_D_matrix_shell(element, material, physical, D_matrix) );

	/* D_matrix, strain_v -> stress_v */
	mul_matrix_vect(ssnum, ssnum, D_matrix, strain_v, stress_v);

	/* strain_v, stress_v -> strain, stress */
	RC_TRY( vect2resultant_stress_strain(ANA_PLANE_SHEAR, stress_v,
	                                     shear_force) );

	RC_TRY( free_work_set(&ws_tmp) );

	return(NORMAL_RC);
}


static RC fill_disp_vector4shell(int dim, ANALYSIS_TYPE ana_type,
                                 WORK_SET ws, NODE_ARRAY node,
                                 ELEMENT element, DISP_ARRAY disp,
                                 double disp_v[])
{
	int ii1;
	int node_index;
	double **dcosine_matrix = ws.mat_3_3[1];
	double *local_disp_v = ws.vec_N[1];
	TRANSLATION3D v_dcosine[3];

	/* 方向余弦マトリクスの設定 */
	RC_TRY( set_elem_dcosine(node, element, v_dcosine) );
	vect2dcosine_array(v_dcosine, dcosine_matrix);

	switch(ana_type){
	case ANA_PLANE_STRESS:
		for(ii1=0; ii1<element.node_num; ii1++){
			node_index = search_disp_node(disp, element.node[ii1]);
			RC_NEG_CHK( node_index );
			local_disp_v[0] = disp.array[node_index].v.x;
			local_disp_v[1] = disp.array[node_index].v.y;
			local_disp_v[2] = disp.array[node_index].v.z;
			/* 座標変換 */
			disp_v[dim*ii1] = local_disp_v[0]*dcosine_matrix[0][0]
			                + local_disp_v[1]*dcosine_matrix[0][1]
			                + local_disp_v[2]*dcosine_matrix[0][2];
			disp_v[dim*ii1+1] = local_disp_v[0]*dcosine_matrix[1][0]
			                  + local_disp_v[1]*dcosine_matrix[1][1]
			                  + local_disp_v[2]*dcosine_matrix[1][2];
		}
		break;
	case ANA_PLANE_BEND:
		for(ii1=0; ii1<element.node_num; ii1++){
			node_index = search_disp_node(disp, element.node[ii1]);
			RC_NEG_CHK( node_index );
			local_disp_v[0] = disp.array[node_index].v.x;
			local_disp_v[1] = disp.array[node_index].v.y;
			local_disp_v[2] = disp.array[node_index].v.z;
			local_disp_v[3] = disp.array[node_index].v.yz;
			local_disp_v[4] = disp.array[node_index].v.zx;
			local_disp_v[5] = disp.array[node_index].v.xy;
			/* 座標変換 */
			disp_v[dim*ii1] = local_disp_v[0]*dcosine_matrix[2][0]
			                + local_disp_v[1]*dcosine_matrix[2][1]
			                + local_disp_v[2]*dcosine_matrix[2][2];
			disp_v[dim*ii1+1] = local_disp_v[3]*dcosine_matrix[0][0]
			                  + local_disp_v[4]*dcosine_matrix[0][1]
			                  + local_disp_v[5]*dcosine_matrix[0][2];
			disp_v[dim*ii1+2] = local_disp_v[3]*dcosine_matrix[1][0]
			                  + local_disp_v[4]*dcosine_matrix[1][1]
			                  + local_disp_v[5]*dcosine_matrix[1][2];
		}

		break;
	case ANA_PLANE_SHEAR:
		for(ii1=0; ii1<element.node_num; ii1++){
			node_index = search_disp_node(disp, element.node[ii1]);
			RC_NEG_CHK( node_index );
			local_disp_v[0] = disp.array[node_index].v.x;
			local_disp_v[1] = disp.array[node_index].v.y;
			local_disp_v[2] = disp.array[node_index].v.z;
			local_disp_v[3] = disp.array[node_index].v.yz;
			local_disp_v[4] = disp.array[node_index].v.zx;
			local_disp_v[5] = disp.array[node_index].v.xy;
			/* 座標変換 */
			disp_v[dim*ii1] = local_disp_v[0]*dcosine_matrix[2][0]
			                + local_disp_v[1]*dcosine_matrix[2][1]
			                + local_disp_v[2]*dcosine_matrix[2][2];
			disp_v[dim*ii1+1] = local_disp_v[3]*dcosine_matrix[0][0]
			                  + local_disp_v[4]*dcosine_matrix[0][1]
			                  + local_disp_v[5]*dcosine_matrix[0][2];
			disp_v[dim*ii1+2] = local_disp_v[3]*dcosine_matrix[1][0]
			                  + local_disp_v[4]*dcosine_matrix[1][1]
			                  + local_disp_v[5]*dcosine_matrix[1][2];
		}
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC vect2stress_strain_shell(ANALYSIS_TYPE ana_type,
                                   const double *vect, STRESS_STRAIN *ss)
{
	switch(ana_type){
	case ANA_PLANE_STRESS:
		ss->v.x  = vect[0];
		ss->v.y  = vect[1];
		ss->v.xy = vect[2];
		break;
	case ANA_PLANE_BEND:
		ss->v.x  -= vect[0];
		ss->v.y  -= vect[1];
		ss->v.xy -= vect[2];
		break;
	case ANA_PLANE_SHEAR:
		ss->v.yz = vect[0];
		ss->v.zx = vect[1];
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC vect2resultant_stress_strain(ANALYSIS_TYPE ana_type,
                                       const double *vect,
                                       STRESS_STRAIN *ss)
{
	switch(ana_type){
	case ANA_PLANE_STRESS:
		ss->v.x  = vect[0];
		ss->v.y  = vect[1];
		ss->v.xy = vect[2];
		break;
	case ANA_PLANE_BEND:
		ss->v.x = vect[0];
		ss->v.y = vect[1];
		ss->v.xy = vect[2];
		break;
	case ANA_PLANE_SHEAR:
		ss->v.yz = vect[0];
		ss->v.zx = vect[1];
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC fem_fill_physical(MATERIAL_PROP *prop, ELEMENT element,
                            PHYSICAL_PROP_ARRAY physical)
{
	double thickness;

	RC_TRY( get_thickness(element, physical, &thickness) );

	switch(prop->ana_type){
	case ANA_PLANE_BEND:
		prop->D_matrix[0][0] *= (thickness * thickness * thickness);
		prop->D_matrix[0][1] *= (thickness * thickness * thickness);
		prop->D_matrix[1][0] *= (thickness * thickness * thickness);
		prop->D_matrix[1][1] *= (thickness * thickness * thickness);
		prop->D_matrix[2][2] *= (thickness * thickness * thickness);
		break;
	case ANA_PLANE_SHEAR:
		prop->D_matrix[0][0] *= thickness;
		prop->D_matrix[1][1] *= thickness;
		break;
	default:
		return(NORMAL_RC);
	}

	return(NORMAL_RC);
}


/* 板問題のひずみ、応力の数を返却 */
static int plane_element_ssnum(ELEMENT element,
                               MATERIAL_PROP_ARRAY material,
                               PHYSICAL_PROP_ARRAY physical)
{
	int ssnum;
	int m_index;

	m_index = search_material_prop_element(material, physical, &(element));
	switch( material.array[m_index].ana_type ){
	case ANA_PLANE_STRESS:
	case ANA_PLANE_STRAIN:
	case ANA_PLANE_BEND:
		ssnum = 3;  /* membrane(xx,yy,xy), bending(xx,yy,xy) */
		break;
	case ANA_PLANE_SHEAR:
		ssnum = 2;  /* shear(yz,zx) */
		break;
	default:
		ssnum = -1;
		break;
	}

	return(ssnum);
}


static int plane_element_dim(ELEMENT element,
                             MATERIAL_PROP_ARRAY material,
                             PHYSICAL_PROP_ARRAY physical)
{
	int dim;
	int m_index;

	m_index = search_material_prop_element(material, physical, &(element));
	switch( material.array[m_index].ana_type ){
	case ANA_PLANE_STRESS:
	case ANA_PLANE_STRAIN:
		dim = 2;
		break;
	case ANA_PLANE_BEND:
	case ANA_PLANE_SHEAR:
		dim = 3;
		break;
	default:
		dim = -1;
		break;
	}

	return(dim);
}


/* 要素の方向余弦マトリクスの設定 */
static RC set_elem_dcosine(NODE_ARRAY node, ELEMENT element,
                           TRANSLATION3D d_cos[])
{
	int ii1, ii2;
	int node_index1;
	int node_index2;
	double *N = NULL;
	double **node_location = NULL;
	double **local_points = NULL;
	TRANSLATION3D vtmp[3];       /* 節点間の距離(ベクトル) */
	TRANSLATION3D mid_point[4];  /* 各辺の中点 */

	switch(element.type){
	case ELEM_TRI1:
	case ELEM_TRI2:
		node_index1 = search_node_label(node, element.node[0]);
		RC_NEG_CHK( node_index1 );
		node_index2 = search_node_label(node, element.node[1]);
		RC_NEG_CHK( node_index2 );
		vtmp[0] = sub_translation3d(node.array[node_index2].p,
		                            node.array[node_index1].p);

		node_index2 = search_node_label(node, element.node[2]);
		RC_NEG_CHK( node_index2 );
		vtmp[1] = sub_translation3d(node.array[node_index2].p,
		                            node.array[node_index1].p);
		break;
	case ELEM_QUAD1:
	case ELEM_QUAD2:
		RC_TRY( allocate1D(MAX_NODE, &N) );
		RC_TRY( allocate2D(MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(MAX_NODE, 4, &local_points) );
		RC_TRY( node_location_matrix(&element, node, node_location) );
		RC_TRY( set_mid_point4quad(local_points) );
		for(ii1=0; ii1<4; ii1++){
			RC_TRY( set_N_vector(element.type, local_points[ii1], N) );
			mid_point[ii1].x = mid_point[ii1].y = mid_point[ii1].z = 0.0;
			for(ii2=0; ii2<element.node_num; ii2++){
				mid_point[ii1].x += node_location[ii2][0] * N[ii2];
				mid_point[ii1].y += node_location[ii2][1] * N[ii2];
				mid_point[ii1].z += node_location[ii2][2] * N[ii2];
			}
		}
		free(N);
		RC_TRY( free2D(MAX_NODE, 3, node_location) );
		RC_TRY( free2D(MAX_NODE, 4, local_points) );

		vtmp[0] = sub_translation3d(mid_point[1], mid_point[3]);
		vtmp[1] = sub_translation3d(mid_point[2], mid_point[0]);
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	d_cos[0] = mul_scalar_translation3d(
	           1.0/(abs_translation3d(vtmp[0])), vtmp[0]);

	vtmp[2] = outer_product3d(vtmp[0], vtmp[1]);
	d_cos[2] = mul_scalar_translation3d(
	           1.0/(abs_translation3d(vtmp[2])), vtmp[2]);

	d_cos[1] = outer_product3d(d_cos[2], d_cos[0]);

	return(NORMAL_RC);
}


/* ローカル座標を四角形要素の各辺の中点に設定する */
static RC set_mid_point4quad(double **coords)
{
	coords[0][0] =  0.0;
	coords[0][1] = -1.0;
	coords[1][0] =  1.0;
	coords[1][1] =  0.0;
	coords[2][0] =  0.0;
	coords[2][1] =  1.0;
	coords[3][0] = -1.0;
	coords[3][1] =  0.0;

	return(NORMAL_RC);
}






