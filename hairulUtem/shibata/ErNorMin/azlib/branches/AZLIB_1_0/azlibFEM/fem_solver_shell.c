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

/* $Id: fem_solver_shell.c,v 1.2 2003/07/22 10:03:00 aoyama Exp $*/

#include <stdio.h>
#include <stdlib.h>
#include "fem_solver.h"
#include "fem_struct.h"
#include "math_utl.h"
#include "rc.h"
#include "nonzero_cg.h"

#define DRILLING         (0.03)  /* $B%I%j%j%s%02sE>9d@-$N%Q%i%a!<%?(B */
#define REDUCED_INTEGRAL (0)     /* 0:$B40A4@QJ,(B */
                                 /* 1:$B<!?tDc8:@QJ,(B */
                                 /* 2:$BA*Br7?Dc8:@QJ,(B */


static RC set_B_matrix_shell(int dim, int node_num,
                             FEM_ANALYSIS_TYPE ana_type, double *N,
                             double **dNxyz, double **B_matrix);
static RC trans_coord_node_location(FEM_NODE_ARRAY node, FEM_ELEMENT element,
                                    FEM_WORK_SET ws, double **coord);
static RC set_Gauss_const4RI(FEM_ELEM_TYPE type, double **points,
                             double *weight, int *num_point);
static RC set_drilling(int dim, FEM_ELEMENT element,
                       FEM_MATERIAL_PROP_ARRAY material,
                       FEM_PHYSICAL_PROP_ARRAY physical,
                       FEM_ELEM_MATRIX *elem_matrix);
static RC trans_coord_elem_matrix(FEM_NODE_ARRAY node, FEM_ELEMENT element,
                                  FEM_WORK_SET ws,
                                  FEM_ELEM_MATRIX *elem_matrix);
static RC fill_disp_vector4shell(int dim, FEM_ANALYSIS_TYPE ana_type,
                                 FEM_WORK_SET ws, FEM_NODE_ARRAY node,
                                 FEM_ELEMENT element, FEM_DISP_ARRAY disp,
                                 double disp_v[]);
static RC vect2stress_strain_shell(FEM_ANALYSIS_TYPE ana_type,
                                   const double *vect, FEM_STRESS_STRAIN *ss);
static RC vect2resultant_stress_strain(FEM_ANALYSIS_TYPE ana_type,
                                       const double *vect,
                                       FEM_STRESS_STRAIN *ss);
static RC fem_fill_physical(FEM_MATERIAL_PROP *prop, FEM_ELEMENT element,
                            FEM_PHYSICAL_PROP_ARRAY physical);
static int plane_element_ssnum(FEM_ELEMENT element,
                               FEM_MATERIAL_PROP_ARRAY material,
                               FEM_PHYSICAL_PROP_ARRAY physical);
static int plane_element_dim(FEM_ELEMENT element,
                             FEM_MATERIAL_PROP_ARRAY material,
                             FEM_PHYSICAL_PROP_ARRAY physical);
static RC set_mid_point4quad(double **coords);
static RC set_elem_dcosine(FEM_NODE_ARRAY node, FEM_ELEMENT element,
                           TRANSLATION3D d_cos[]);
static void vect2dcosine_array(TRANSLATION3D v[], double **array);


/* $B<o!9$N9d@-!JKl!"6J$2!"$;$sCG!K%^%H%j%/%9$r(B */
/* $BBP1~$9$k2U=j$KAH$_9~$_!"MWAG9d@-%^%H%j%/%9$r:n@.!J%7%'%k2r@OMQ!K(B */
RC fem_impose_shell_elem_matrix(int plane_dim, FEM_ELEMENT element,
                                FEM_MATERIAL_PROP_ARRAY material,
                                FEM_PHYSICAL_PROP_ARRAY physical,
                                FEM_ELEM_MATRIX elem_matrix,
                                FEM_ELEM_MATRIX *total_elem_matrix)
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


/* $BMWAG9d@-%^%H%j%C%/%9$N:n@.!J%7%'%k2r@OMQ!K(B */
RC fem_make_shell_matrix(FEM_ELEMENT element, FEM_NODE_ARRAY node,
                         FEM_MATERIAL_PROP_ARRAY material,
                         FEM_PHYSICAL_PROP_ARRAY physical,
                         FEM_ELEM_MATRIX *elem_matrix)
{
	int ii1, ii2;
	int plane_dim;
	int m_index, phys_index, renum_index;
	FEM_WORK_SET ws;
	FEM_ELEM_MATRIX tmp_elem_matrix;   /* each elem_matrix */

	RC_TRY( allocate_fem_work_set(&ws) );

	/* elem_matrix->matrix[][] $B$N3NJ](B */
	/* preparation of total elem_matrix */
	elem_matrix->dim = 6;
	elem_matrix->element = element.label;
	elem_matrix->size = elem_matrix->dim * element.node_num;
	RC_TRY( allocate2D(elem_matrix->size, elem_matrix->size,
	                   &(elem_matrix->matrix)) );
	/* elem_matirx->index[] $B$N3NJ]$HBeF~(B */
	RC_TRY( allocate1D_i((elem_matrix->size), &(elem_matrix->index)) );
	for(ii1=0; ii1<element.node_num; ii1++){
		renum_index = node_renum_index(node, element.node[ii1]);
		RC_NEG_CHK(renum_index);
		for(ii2=0; ii2<(elem_matrix->dim); ii2++){
			elem_matrix->index[ii1*(elem_matrix->dim)+ii2]
			           = renum_index*(elem_matrix->dim) + ii2;
		}
	}

	phys_index = search_fem_physical_prop_label(physical, element.physical);
	RC_NEG_CHK(phys_index);

	/* $BKl9d@-(B */
	m_index = search_fem_material_prop_element(material, physical, &(element) );
	RC_NEG_CHK(m_index);
	material.array[m_index].ana_type = ANA_PLANE_STRESS;
	plane_dim = 2;
	RC_TRY( fem_make_shell_elem_matrix(element, node, material,
	                                    physical, &tmp_elem_matrix, ws) );
	RC_TRY( fem_impose_shell_elem_matrix(plane_dim, element, material,
	                                     physical, tmp_elem_matrix,
	                                     elem_matrix) );

	RC_TRY( free_fem_elem_matrix(&tmp_elem_matrix) );

	/* $B6J$29d@-(B */
	m_index = search_fem_material_prop_element(material, physical, &(element) );
	RC_NEG_CHK(m_index);
	material.array[m_index].ana_type = ANA_PLANE_BEND;
	plane_dim = 3;
	RC_TRY( fem_make_shell_elem_matrix(element, node, material,
	                                    physical, &tmp_elem_matrix, ws) );
	RC_TRY( fem_impose_shell_elem_matrix(plane_dim, element, material,
	                                     physical, tmp_elem_matrix,
	                                     elem_matrix) );
	RC_TRY( free_fem_elem_matrix(&tmp_elem_matrix) );

	/* $B$;$sCG9d@-(B */
	m_index = search_fem_material_prop_element(material, physical, &(element) );
	RC_NEG_CHK(m_index);
	material.array[m_index].ana_type = ANA_PLANE_SHEAR;
	plane_dim = 3; 
	RC_TRY( fem_make_shell_elem_matrix(element, node, material,
	                                    physical, &tmp_elem_matrix, ws) );
	RC_TRY( fem_impose_shell_elem_matrix(plane_dim, element, material,
	                                     physical, tmp_elem_matrix,
	                                     elem_matrix) );
	RC_TRY( free_fem_elem_matrix(&tmp_elem_matrix) );

	/* $B%I%j%j%s%02sE>9d@-(B */
	RC_TRY( element_volume(&element, node) );
	RC_TRY( set_drilling(elem_matrix->dim, element, material, physical,
	                     elem_matrix) );

	/* $BMWAG9d@-%^%H%j%/%9$N:BI8JQ49(B */
	RC_TRY( trans_coord_elem_matrix(node, element, ws, elem_matrix) );

	RC_TRY( free_fem_work_set(&ws) );

	return(NORMAL_RC);
}


/* $B<o!9$NMWAG9d@-%^%H%j%C%/%9$N:n@.!J%7%'%k2r@OMQ!K(B */
RC fem_make_shell_elem_matrix(FEM_ELEMENT element, FEM_NODE_ARRAY node,
                              FEM_MATERIAL_PROP_ARRAY material,
                              FEM_PHYSICAL_PROP_ARRAY physical,
                              FEM_ELEM_MATRIX *elem_matrix,
                              FEM_WORK_SET ws)
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
	FEM_WORK_SET ws_tmp;

	if(element.label < 0) return(ARG_ERROR_RC);

	RC_TRY( allocate_fem_work_set(&ws_tmp) );

	RC_NEG_CHK( ssnum = plane_element_ssnum(element, material, physical) );

	elem_matrix->dim = plane_element_dim(element, material, physical);
	RC_NEG_CHK(elem_matrix->dim);
	elem_matrix->element = element.label;
	elem_matrix->size = (element.node_num) * (elem_matrix->dim);

	/* elem_matrix->array[][] $B$N3NJ](B */
	RC_TRY( allocate2D(elem_matrix->size,
	                   elem_matrix->size, &(elem_matrix->matrix)) );
	/* elem_matirx->index[] $B$N3NJ]$HBeF~(B */
	RC_TRY( allocate1D_i((elem_matrix->size), &(elem_matrix->index)) );

	/* total_elem_matrix$B$KBeF~$9$k$?$a$N%$%s%G%C%/%9(B */
	for(ii1=0; ii1<element.node_num; ii1++){
		for(ii2=0; ii2<(elem_matrix->dim); ii2++){
			elem_matrix->index[ii1*(elem_matrix->dim)+ii2]
			             = ii1*(elem_matrix->dim) + ii2;
		}
	}

	/* D $B%^%H%j%C%/%9(B */
	RC_TRY( get_D_matrix_shell(element, material, physical, D_matrix) );
	RC_TRY( get_thickness(element, physical, &thickness) );

	/* $B!J6I=j:BI87O$G$N!K@aE@:BI8CM$N%^%H%j%C%/%9$r:n@.(B */
	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( trans_coord_node_location(node, element, ws_tmp, node_location));

	/* $B%,%&%9%]%$%s%H$N@_Dj(B */
	if(REDUCED_INTEGRAL == 0){   /* $B40A4@QJ,(B */
		RC_TRY( Gauss_const_default(element.type, g_points, weights,
		                            &num_g_point));
	}else if(REDUCED_INTEGRAL == 1){   /* $B<!?tDc8:@QJ,(B */
		RC_TRY( set_Gauss_const4RI(element.type, g_points, weights,
		                           &num_g_point) );
	}else{   /* $BA*Br7?<!?tDc8:@QJ,(B($B$;$sCG9`$N$_$r<!?tDc8:@QJ,(B) */
		RC_NEG_CHK( m_index = search_fem_material_prop_element(material,
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
		/* B $B%^%H%j%/%9(B */
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
	RC_TRY( free_fem_work_set(&ws_tmp) );

	return(NORMAL_RC);
}


/* D $B%^%H%j%C%/%9!J%7%'%k2r@OMQ!K(B */
RC get_D_matrix_shell(FEM_ELEMENT element, FEM_MATERIAL_PROP_ARRAY material,
                      FEM_PHYSICAL_PROP_ARRAY physical, double **D_matrix)
{
	int ii1, ii2;
	int dim;
	int m_index;
	int ssnum;
	FEM_MATERIAL_PROP tmp_mat;

	m_index = search_fem_material_prop_element(material, physical, &element);
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


/* B $B%^%H%j%/%9$H(B |J|$B!J%7%'%k2r@OMQ!K(B */
RC make_B_matrix_shell(FEM_ELEMENT element, FEM_MATERIAL_PROP_ARRAY material,
                       FEM_PHYSICAL_PROP_ARRAY physical, FEM_WORK_SET ws,
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
	m_index = search_fem_material_prop_element(material, physical, &element);
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

	/* B $B%^%H%j%/%9(B */
	plane_dim = plane_element_dim(element, material, physical);
	if( (plane_dim != 2)&&(plane_dim != 3) ) return(ARG_ERROR_RC);
	RC_TRY( set_B_matrix_shell(plane_dim, element.node_num,
	                           material.array[m_index].ana_type,  N, dNxyz,
	                           B_matrix) );

	if(det_J != NULL){
		/* $BCm0U(B : $B$3$N%d%3%S%"%s$O6I=j:BI8$K$h$k$b$N(B */
		*det_J = determinant(elem_dim, Jacobi);
	}

	return(NORMAL_RC);
}


/* N, dNxyz $B$rMQ$$$F(B B $B%^%H%j%/%9$r:n@.(B($B%7%'%k2r@OMQ(B) */
/* dNxyz[0][i]...Ni $B$N(B x $BJ}8~HyJ,(B */
/* dNxyz[1][i]...Ni $B$N(B y $BJ}8~HyJ,(B */
/* dNxyz[2][i]...Ni $B$N(B z $BJ}8~HyJ,(B */
/* dim=2 ...2$B<!85(B, dim=3 ...3$B<!85(B */
static RC set_B_matrix_shell(int dim, int node_num,
                             FEM_ANALYSIS_TYPE ana_type, double *N,
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


/* $BMWAG9d@-%^%H%j%/%9$N:BI8JQ49(B */
/* [K_GLOBAL] = [T]^T * [K_LOCAL] * [T] */
/* ($B%9%Q!<%9@-$r9MN8$7$?9TNs@Q$rMxMQ(B) */
static RC trans_coord_elem_matrix(FEM_NODE_ARRAY node, FEM_ELEMENT element,
                                  FEM_WORK_SET ws, FEM_ELEM_MATRIX *elem_matrix)
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


/* $B!J6I=j:BI87O$G$N!K@aE@:BI8CM$N%^%H%j%C%/%9$r:n@.(B */
static RC trans_coord_node_location(FEM_NODE_ARRAY node, FEM_ELEMENT element,
                                    FEM_WORK_SET ws, double **coord)
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

	/* $BJ}8~M>89%^%H%j%/%9$N@_Dj(B */
	RC_TRY( set_elem_dcosine(node, element, v_dcosine) );
	vect2dcosine_array(v_dcosine, dcosine_matrix);

	transpose_overwrite_matrix(3, dcosine_matrix);
	/* $B:BI8JQ49(B */
	mul_matrix(element.node_num, 3, 3, vtmp, dcosine_matrix, coord);

	return(NORMAL_RC);
}


/* $B%I%j%j%s%02sE>9d@-(B ($BLL$NK!@~2s$j$N2sE>9d@-(B) */
static RC set_drilling(int dim, FEM_ELEMENT element,
                       FEM_MATERIAL_PROP_ARRAY material,
                       FEM_PHYSICAL_PROP_ARRAY physical,
                       FEM_ELEM_MATRIX *elem_matrix)
{
	int ii1, ii2;
	int ix1, ix2;
	int m_index;
	double coeff1, coeff2; /* $B78?t(B */
	double thickness;

	/* $B78?t$N7W;;(B */
	m_index = search_fem_material_prop_element(material, physical, &element);
	RC_NEG_CHK(m_index);
	RC_TRY( get_thickness(element, physical, &thickness) );
	coeff1 = DRILLING * material.array[m_index].E * thickness * element.volume;
	coeff2 = - coeff1 * ( 1.0/(double)(element.node_num-1) );

	/* $BBP1~$9$k9T$HNs$KAH$_9~$`(B */
	for(ii1=0; ii1<element.node_num; ii1++){
		ix1 = dim*(ii1+1) - 1;
		for(ii2=0; ii2<element.node_num; ii2++){
			ix2 = dim*(ii2+1) - 1;
			elem_matrix->matrix[ix1][ix2] = coeff2;
		}
	}
	/* $BBP3QMWAG$KBeF~(B */
	for(ii1=0; ii1<element.node_num; ii1++){
		ix1 =  dim*(ii1+1) - 1;
		elem_matrix->matrix[ix1][ix1] = coeff1;
	}

	return(NORMAL_RC);
}


/* $B<!?tDc8:@QJ,MQ$K@QJ,E@$r@_Dj(B */
static RC set_Gauss_const4RI(FEM_ELEM_TYPE type, double **points,
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


/* $B1~NO!"OD$N7W;;(B */
/* stress, strain $B$,(B NULL $B$J$iBeF~$7$J$$(B */
/* $B1~NO!"OD$r7W;;$9$k0LCV!JAX!K$r(Bdistance[] $B$G;XDj(B */
/* distance$B$,(BNULL$B$J$i!"(B-thickness/2$B$H(Bthickness/2 ($BN"$HI=(B)$B$K@_Dj(B */
RC fem_stress_strain_shell(int fem_sol_ss_flag, double *distance,
                           FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                           FEM_DISP_ARRAY disp,
                           FEM_MATERIAL_PROP_ARRAY material,
                           FEM_PHYSICAL_PROP_ARRAY physical,
                           FEM_STRESS_ARRAY stress[], FEM_STRAIN_ARRAY strain[])
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
	FEM_STRESS_STRAIN local_stress[2], local_strain[2];
	FEM_WORK_SET ws;

	if( ((stress == NULL)&&(strain == NULL)) || (fem_sol_ss_flag == 0) ){
		return(NORMAL_RC);
	}

	if(fem_sol_ss_flag & FEM_SOL_SS_NODE_AV){
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
		RC_TRY( allocate_fem_stress_array(0, &(stress[0])) );
		RC_TRY( allocate_fem_stress_array(0, &(stress[1])) );
		stress[0].source = stress[1].source = FEM_K_SOL;
	}
	if(strain != NULL){
		RC_TRY( allocate_fem_strain_array(0, &(strain[0])) );
		RC_TRY( allocate_fem_strain_array(0, &(strain[1])) );
		strain[0].source = strain[1].source = FEM_K_SOL;
	}

	RC_TRY( allocate_fem_work_set(&ws) );
	local_points = ws.mat_N_4[0];

	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		if(distance == NULL){ /* $B1~NO!"OD$r7W;;$9$k0LCV$rHD$NN"$HI=$K@_Dj(B */
			RC_TRY( get_thickness(element.array[ii1], physical, &thickness) );
			default_distance[0] = -0.5*thickness;
			default_distance[1] =  0.5*thickness;
		}else{ /* $B1~NO!"OD$r7W;;$9$k0LCV(B($BAX(B)$B$r%f!<%6;XDjCM$K@_Dj(B */
			default_distance[0] = distance[0];
			default_distance[1] = distance[1];
		}

		/* $BMWAGCf?4(B */
		if(fem_sol_ss_flag & FEM_SOL_SS_CENTER){
			RC_TRY( set_center_point(element.array[ii1].type,
                                     local_points[0]) );
			RC_TRY( fem_local_stress_strain_shell(default_distance,
			           element.array[ii1], local_points[0], node, disp,
			           material, physical, ws, local_stress, local_strain) );
			if(stress != NULL){
				RC_TRY( realloc_fem_stress_array(&(stress[0])) );
				RC_TRY( realloc_fem_stress_array(&(stress[1])) );
				stress[0].array[stress[0].size] = local_stress[0];
				stress[1].array[stress[1].size] = local_stress[1];
				stress[0].array[stress[0].size].node = -1;  /* center */
				stress[1].array[stress[1].size].node = -1;  /* center */
				(stress[0].size)++;
				(stress[1].size)++;
			}
			if(strain != NULL){
				RC_TRY( realloc_fem_strain_array(&(strain[0])) );
				RC_TRY( realloc_fem_strain_array(&(strain[1])) );
				strain[0].array[strain[0].size] = local_strain[0];
				strain[1].array[strain[1].size] = local_strain[1];
				strain[0].array[strain[0].size].node = -1;  /* center */
				strain[1].array[strain[1].size].node = -1;  /* center */
				(strain[0].size)++;
				(strain[1].size)++;
			}
		}

		/* $B3F@aE@(B */
		if( (fem_sol_ss_flag & FEM_SOL_SS_NODE)
          ||(fem_sol_ss_flag & FEM_SOL_SS_NODE_AV) ){
			RC_TRY( set_local_node_points(element.array[ii1].type,
			                              local_points) );

			for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
				RC_TRY( fem_local_stress_strain_shell(default_distance,
				            element.array[ii1], local_points[ii2], node, disp,
				            material, physical, ws, local_stress,
				            local_strain) );

				if(fem_sol_ss_flag & FEM_SOL_SS_NODE_AV){
					node_index = search_fem_node_label(node,
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
					RC_TRY( realloc_fem_stress_array(&(stress[0])) );
					RC_TRY( realloc_fem_stress_array(&(stress[1])) );
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
					RC_TRY( realloc_fem_strain_array(&(strain[0])) );
					RC_TRY( realloc_fem_strain_array(&(strain[1])) );
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
		RC_TRY( clean_fem_stress_array(&(stress[0])) );
		RC_TRY( clean_fem_stress_array(&(stress[1])) );
	}
	if(strain != NULL){
		RC_TRY( clean_fem_strain_array(&(strain[0])) );
		RC_TRY( clean_fem_strain_array(&(strain[1])) );
	}

	if(fem_sol_ss_flag & FEM_SOL_SS_NODE_AV){
		if(stress != NULL){
			for(ii1=0; ii1<stress[0].size; ii1++){
				if(stress[0].array[ii1].element < 0) continue;
				node_index = search_fem_node_label(node,
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
				node_index = search_fem_node_label(node,
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

	RC_TRY( free_fem_work_set(&ws) );

	return(NORMAL_RC);
}


/* $B%7%'%kMWAG$N6I=j1~NO!"$R$:$_(B */
RC fem_local_stress_strain_shell(double distance[], FEM_ELEMENT element,
                                 const double *local_point,
                                 FEM_NODE_ARRAY node, FEM_DISP_ARRAY disp,
                                 FEM_MATERIAL_PROP_ARRAY material,
                                 FEM_PHYSICAL_PROP_ARRAY physical,
                                 FEM_WORK_SET ws, FEM_STRESS_STRAIN stress[],
                                 FEM_STRESS_STRAIN strain[])
{
	int ii1;
	int dim, m_index, phys_index;
	int ssnum;
	double thickness, default_distance[2];
	double disp_v[3*FEM_MAX_NODE];
	double strain_v[2][6];
	double stress_v[2][6];
	double **node_location = ws.mat_N_3[0];
	double **B_matrix = ws.mat_6_3N[0];
	double **D_matrix = ws.mat_6_6[0];
	FEM_WORK_SET ws_tmp;

	if(element.label < 0) return(ARG_ERROR_RC);

	m_index = search_fem_material_prop_element(material, physical, &(element));
	RC_NEG_CHK(m_index);

	phys_index = search_fem_physical_prop_label(physical, element.physical);
	RC_NEG_CHK(phys_index)

	RC_TRY( allocate_fem_work_set(&ws_tmp) );

	for(ii1=0; ii1<2; ii1++){
		init_fem_stress_strain(&(stress[ii1]));
		init_fem_stress_strain(&(strain[ii1]));
		stress[ii1].element = strain[ii1].element = element.label;
	}

	if(distance == NULL){ /* $B1~NO!"OD$r7W;;$9$k0LCV$rHD$NN"$HI=$K@_Dj(B */
		RC_TRY( get_thickness(element, physical, &thickness) );
		default_distance[0] = -0.5*thickness;
		default_distance[1] =  0.5*thickness;
	}else{ /* $B1~NO!"OD$r7W;;$9$k0LCV$r%f!<%6;XDjCM$K@_Dj(B */
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

	RC_TRY( free_fem_work_set(&ws_tmp) );

	return(NORMAL_RC);
}


/* $B9g1~NO!"6JN($N7W;;(B */
/*  force $B!A(B curv $B$,(B NULL $B$J$iBeF~$7$J$$(B */
RC fem_resultant_stress_strain_shell(int fem_sol_ss_flag, FEM_NODE_ARRAY node,
                                     FEM_ELEMENT_ARRAY element,
                                     FEM_DISP_ARRAY disp,
                                     FEM_MATERIAL_PROP_ARRAY material,
                                     FEM_PHYSICAL_PROP_ARRAY physical,
                                     FEM_STRESS_ARRAY *force,
                                     FEM_STRESS_ARRAY *moment,
                                     FEM_STRESS_ARRAY *shear_force,
                                     FEM_STRAIN_ARRAY *curv)
{
	int ii1, ii2;
	int node_index;
	int *count = NULL;
	double **local_points = NULL;
	TRANS_ROTATION3D *force_sum = NULL, *moment_sum = NULL;
	TRANS_ROTATION3D *shear_force_sum = NULL, *curv_sum = NULL;
	TRANS_ROTATION3D init_v = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	FEM_STRESS_STRAIN local_force, local_moment, local_shear_force;
	FEM_STRESS_STRAIN local_curv;
	FEM_WORK_SET ws;

	if( ((force == NULL)&&(moment == NULL)&&(shear_force == NULL)
	     &&(curv == NULL)) || (fem_sol_ss_flag == 0) ){
		return(NORMAL_RC);
	}

	if(fem_sol_ss_flag & FEM_SOL_SS_NODE_AV){
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
		RC_TRY( allocate_fem_stress_array(0, force) );
		force->source = FEM_K_SOL;
	}
	if(moment != NULL){
		RC_TRY( allocate_fem_stress_array(0, moment) );
		moment->source = FEM_K_SOL;
	}
	if(shear_force != NULL){
		RC_TRY( allocate_fem_stress_array(0, shear_force) );
		shear_force->source = FEM_K_SOL;
	}
	if(curv != NULL){
		RC_TRY( allocate_fem_strain_array(0, curv) );
		curv->source = FEM_K_SOL;
	}

	RC_TRY( allocate_fem_work_set(&ws) );
	local_points = ws.mat_N_4[0];

	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		/* $BMWAGCf?4(B */
		if(fem_sol_ss_flag & FEM_SOL_SS_CENTER){
			RC_TRY( set_center_point(element.array[ii1].type,
                                     local_points[0]) );
			RC_TRY( fem_local_resultant_stress_strain(element.array[ii1],
			           local_points[0], node, disp, material, physical, ws,
			           &local_force, &local_moment, &local_shear_force,
			           &local_curv) );
			if(force != NULL){
				RC_TRY( realloc_fem_stress_array(force) );
				force->array[force->size] = local_force;
				force->array[force->size].node = -1;  /* center */
				(force->size)++;
			}
			if(moment != NULL){
				RC_TRY( realloc_fem_stress_array(moment) );
				moment->array[moment->size] = local_moment;
				moment->array[moment->size].node = -1;  /* center */
				(moment->size)++;
			}
			if(shear_force != NULL){
				RC_TRY( realloc_fem_stress_array(shear_force) );
				shear_force->array[shear_force->size] = local_shear_force;
				shear_force->array[shear_force->size].node = -1;  /* center */
				(shear_force->size)++;
			}
			if(curv != NULL){
				RC_TRY( realloc_fem_strain_array(curv) );
				curv->array[curv->size] = local_curv;
				curv->array[curv->size].node = -1;  /* center */
				(curv->size)++;
			}
		}

		/* $B3F@aE@(B */
		if( (fem_sol_ss_flag & FEM_SOL_SS_NODE)
          ||(fem_sol_ss_flag & FEM_SOL_SS_NODE_AV) ){
			RC_TRY( set_local_node_points(element.array[ii1].type,
			                              local_points) );

			for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
				RC_TRY( fem_local_resultant_stress_strain(element.array[ii1],
				            local_points[ii2], node, disp, material, physical,
				            ws, &local_force, &local_moment,
				            &local_shear_force, &local_curv) );

				if(fem_sol_ss_flag & FEM_SOL_SS_NODE_AV){
					node_index = search_fem_node_label(node,
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
					RC_TRY( realloc_fem_stress_array(force) );
					force->array[force->size] = local_force;
					force->array[force->size].node 
					     = element.array[ii1].node[ii2];
					(force->size)++;
				}
				if(moment != NULL){
					RC_TRY( realloc_fem_stress_array(moment) );
					moment->array[moment->size] = local_moment;
					moment->array[moment->size].node
					      = element.array[ii1].node[ii2];
					(moment->size)++;
				}
				if(shear_force != NULL){
					RC_TRY( realloc_fem_stress_array(shear_force) );
					shear_force->array[shear_force->size] = local_shear_force;
					shear_force->array[shear_force->size].node
					      = element.array[ii1].node[ii2];
					(shear_force->size)++;
				}
				if(curv != NULL){
					RC_TRY( realloc_fem_strain_array(curv) );
					curv->array[curv->size] = local_curv;
					curv->array[curv->size].node
					      = element.array[ii1].node[ii2];
					(curv->size)++;
				}
			}
		}
	}

	if(force != NULL) RC_TRY( clean_fem_stress_array(force) );
	if(moment != NULL) RC_TRY( clean_fem_stress_array(moment) );
	if(shear_force != NULL) RC_TRY( clean_fem_stress_array(shear_force) );
	if(curv != NULL) RC_TRY( clean_fem_strain_array(curv) );

	if(fem_sol_ss_flag & FEM_SOL_SS_NODE_AV){
		if(force != NULL){
			for(ii1=0; ii1<force->size; ii1++){
				if(force->array[ii1].element < 0) continue;
				node_index = search_fem_node_label(node,
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
				node_index = search_fem_node_label(node,
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
				node_index = search_fem_node_label(node,
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
				node_index = search_fem_node_label(node,
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

	RC_TRY( free_fem_work_set(&ws) );

	return(NORMAL_RC);
}


/* $B%7%'%kMWAG$N6I=j9g1~NO!"6JN((B */
RC fem_local_resultant_stress_strain(FEM_ELEMENT element,
                                     const double *local_point,
                                     FEM_NODE_ARRAY node, FEM_DISP_ARRAY disp,
                                     FEM_MATERIAL_PROP_ARRAY material,
                                     FEM_PHYSICAL_PROP_ARRAY physical,
                                     FEM_WORK_SET ws, FEM_STRESS_STRAIN *force,
                                     FEM_STRESS_STRAIN *moment,
                                     FEM_STRESS_STRAIN *shear_force,
                                     FEM_STRESS_STRAIN *curv)
{
	int dim, m_index, phys_index,  ssnum;
	double disp_v[3*FEM_MAX_NODE];
	double strain_v[6];
	double stress_v[6];
	double **node_location = ws.mat_N_3[0];
	double **B_matrix = ws.mat_6_3N[0];
	double **D_matrix = ws.mat_6_6[0];
	FEM_WORK_SET ws_tmp;

	if(element.label < 0) return(ARG_ERROR_RC);

	m_index = search_fem_material_prop_element(material, physical, &(element));
	RC_NEG_CHK(m_index);

	phys_index = search_fem_physical_prop_label(physical, element.physical);
	RC_NEG_CHK(phys_index)

	RC_TRY( allocate_fem_work_set(&ws_tmp) );

	init_fem_stress_strain(force);
	init_fem_stress_strain(moment);
	init_fem_stress_strain(shear_force);
	force->element = moment->element = shear_force->element = element.label;

	init_fem_stress_strain(curv);
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

	RC_TRY( free_fem_work_set(&ws_tmp) );

	return(NORMAL_RC);
}


static RC fill_disp_vector4shell(int dim, FEM_ANALYSIS_TYPE ana_type,
                                 FEM_WORK_SET ws, FEM_NODE_ARRAY node,
                                 FEM_ELEMENT element, FEM_DISP_ARRAY disp,
                                 double disp_v[])
{
	int ii1;
	int node_index;
	double **dcosine_matrix = ws.mat_3_3[1];
	double *local_disp_v = ws.vec_N[1];
	TRANSLATION3D v_dcosine[3];

	/* $BJ}8~M>89%^%H%j%/%9$N@_Dj(B */
	RC_TRY( set_elem_dcosine(node, element, v_dcosine) );
	vect2dcosine_array(v_dcosine, dcosine_matrix);

	switch(ana_type){
	case ANA_PLANE_STRESS:
		for(ii1=0; ii1<element.node_num; ii1++){
			node_index = search_fem_disp_node(disp, element.node[ii1]);
			RC_NEG_CHK( node_index );
			local_disp_v[0] = disp.array[node_index].v.x;
			local_disp_v[1] = disp.array[node_index].v.y;
			local_disp_v[2] = disp.array[node_index].v.z;
			/* $B:BI8JQ49(B */
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
			node_index = search_fem_disp_node(disp, element.node[ii1]);
			RC_NEG_CHK( node_index );
			local_disp_v[0] = disp.array[node_index].v.x;
			local_disp_v[1] = disp.array[node_index].v.y;
			local_disp_v[2] = disp.array[node_index].v.z;
			local_disp_v[3] = disp.array[node_index].v.yz;
			local_disp_v[4] = disp.array[node_index].v.zx;
			local_disp_v[5] = disp.array[node_index].v.xy;
			/* $B:BI8JQ49(B */
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
			node_index = search_fem_disp_node(disp, element.node[ii1]);
			RC_NEG_CHK( node_index );
			local_disp_v[0] = disp.array[node_index].v.x;
			local_disp_v[1] = disp.array[node_index].v.y;
			local_disp_v[2] = disp.array[node_index].v.z;
			local_disp_v[3] = disp.array[node_index].v.yz;
			local_disp_v[4] = disp.array[node_index].v.zx;
			local_disp_v[5] = disp.array[node_index].v.xy;
			/* $B:BI8JQ49(B */
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


static RC vect2stress_strain_shell(FEM_ANALYSIS_TYPE ana_type,
                                   const double *vect, FEM_STRESS_STRAIN *ss)
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


static RC vect2resultant_stress_strain(FEM_ANALYSIS_TYPE ana_type,
                                       const double *vect,
                                       FEM_STRESS_STRAIN *ss)
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


static RC fem_fill_physical(FEM_MATERIAL_PROP *prop, FEM_ELEMENT element,
                            FEM_PHYSICAL_PROP_ARRAY physical)
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


/* $BHDLdBj$N$R$:$_!"1~NO$N?t$rJV5Q(B */
static int plane_element_ssnum(FEM_ELEMENT element,
                               FEM_MATERIAL_PROP_ARRAY material,
                               FEM_PHYSICAL_PROP_ARRAY physical)
{
	int ssnum;
	int m_index;

	m_index = search_fem_material_prop_element(material, physical, &(element));
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


static int plane_element_dim(FEM_ELEMENT element,
                             FEM_MATERIAL_PROP_ARRAY material,
                             FEM_PHYSICAL_PROP_ARRAY physical)
{
	int dim;
	int m_index;

	m_index = search_fem_material_prop_element(material, physical, &(element));
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


/* $BMWAG$NJ}8~M>89%^%H%j%/%9$N@_Dj(B */
static RC set_elem_dcosine(FEM_NODE_ARRAY node, FEM_ELEMENT element,
                           TRANSLATION3D d_cos[])
{
	int ii1, ii2;
	int node_index1;
	int node_index2;
	double *N = NULL;
	double **node_location = NULL;
	double **local_points = NULL;
	TRANSLATION3D vtmp[3];       /* $B@aE@4V$N5wN%(B($B%Y%/%H%k(B) */
	TRANSLATION3D mid_point[4];  /* $B3FJU$NCfE@(B */

	switch(element.type){
	case ELEM_TRI1:
	case ELEM_TRI2:
		node_index1 = search_fem_node_label(node, element.node[0]);
		RC_NEG_CHK( node_index1 );
		node_index2 = search_fem_node_label(node, element.node[1]);
		RC_NEG_CHK( node_index2 );
		vtmp[0] = sub_translation3d(node.array[node_index2].p,
		                            node.array[node_index1].p);

		node_index2 = search_fem_node_label(node, element.node[2]);
		RC_NEG_CHK( node_index2 );
		vtmp[1] = sub_translation3d(node.array[node_index2].p,
		                            node.array[node_index1].p);
		break;
	case ELEM_QUAD1:
	case ELEM_QUAD2:
		RC_TRY( allocate1D(FEM_MAX_NODE, &N) );
		RC_TRY( allocate2D(FEM_MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(FEM_MAX_NODE, 4, &local_points) );
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
		RC_TRY( free2D(FEM_MAX_NODE, 3, node_location) );
		RC_TRY( free2D(FEM_MAX_NODE, 4, local_points) );

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


/* $B%m!<%+%k:BI8$r;M3Q7AMWAG$N3FJU$NCfE@$K@_Dj$9$k(B */
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






