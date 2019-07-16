#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "mathematics.h"
#include "fem.h"
#include "base.h"


RC make_elem_matrix4shell(ELEMENT element, NODE_ARRAY node,
                          MATERIAL_PROP_ARRAY material,
                          PHYSICAL_PROP_ARRAY physical,
                          WORK_SET ws,
                          ELEM_MATRIX *elem_matrix);
static RC make_drilling_matrix(ELEMENT element, NODE_ARRAY node,
                               MATERIAL_PROP_ARRAY material,
                               PHYSICAL_PROP_ARRAY physical,
                               double thickness,
                               double **drilling_matrix);
static RC set_B_matrix_bend(int node_num, double **dNxy, double **B_matrix);
static RC set_B_matrix_shear(ELEM_TYPE type, int node_num, double **dNw_xy,
                             double *N, double **B_matrix);
static RC set_B_matrix_membrane(int node_num, double **dNxy, double **B_matrix);
static RC get_D_matrix4shell(ANALYSIS_TYPE type, ELEMENT element,
                             MATERIAL_PROP_ARRAY material,
                             PHYSICAL_PROP_ARRAY physical,
                             double **D_matrix);
static RC set_N_vector4shell(ELEM_TYPE type, const double *coord,
                             double **node_location, double *Nw, double *N);
static RC N_tri1_shell(const double *l1_l2_l3, double **xy, double *Nw,
                       double *N);
static RC N_tri2_shell(const double *l1_l2_l3, double *Nw, double *N);
#if 0
static RC N_tri2_shell_ver2(const double *l1_l2_l3, double *Nw, double *N);
#endif
static RC N_quad1_shell(const double *xi_eta, double **xy, double *Nw,
                        double *N);
static RC N_quad2_shell(const double *xi_eta, double **xy, double *Nw,
                        double *N);
static RC set_dN_matrix4shell(ELEM_TYPE type, const double *coord,
                              double **node_location, double **dNw,
                              double **dN);
static RC dN_tri1_shell(const double *l1_l2_l3, double **xy, double **dNw,
                        double **dN);
static RC dN_tri2_shell(const double *l1_l2_l3, double **dNw, double **dN);
#if 0
static RC dN_tri2_shell_ver2(const double *l1_l2_l3, double **dNw, double **dN);
#endif
static RC dN_quad1_shell(const double *xi_eta, double **xy, double **dNw,
                         double **dN);
static RC dN_quad1_shell_sub(const double P[], double **xy, double *dNw,
                             double *dN);
static RC dN_quad2_shell(const double *xi_eta, double **xy, double **dNw,
                         double **dN);
static RC dN_quad2_shell_sub(const double *xi_eta, double **xy,
                             double M[][9], double D[]);
static RC set_L_matrix_shell(ELEM_TYPE type, double **node_location,
                             double L[][3]);
static RC local_node_location(int node_num, double L[][3],
                              double **node_location);
static double element_length(int node_num, double **node_location);
int main(int argc, char **argv);



int main(int argc, char **argv)
{
#if 0
	if(argc < 2){
		fprintf(stderr, "Usage : %s [r:model(.bdf)]\n", argv[0]);
		return(EXIT_FAILURE);
	}
#endif
	return(EXIT_SUCCESS);
}


RC make_elem_matrix4shell (ELEMENT element, NODE_ARRAY node,
                           MATERIAL_PROP_ARRAY material,
                           PHYSICAL_PROP_ARRAY physical,
                           WORK_SET ws,
                           ELEM_MATRIX *elem_matrix)
{
	int ii1, ii2, ii3, ii4, ii5;
	int dim;
	int g_num_point;
	double thickness;
	double length;
	double k;
	double beta;
	double L[3][3];
	double *g_weights = ws.vec_G[0];
	double **g_points = ws.mat_G_4[0];
	double **node_location = ws.mat_N_3[0];
	double **drilling_matrix = ws.mat_3N_3N[0];


	if(element.label < 0) return(ARG_ERROR_RC);
	if(element_dim(element.type) != 2) return(ARG_ERROR_RC);

	dim = elem_matrix->dim = 6;
	elem_matrix->element = element.label;
	elem_matrix->size = (element.node_num) * dim;
	RC_TRY( allocate2D(elem_matrix->size, elem_matrix->size,
	                   &(elem_matrix->matrix)) );
	RC_TRY( allocate1I(elem_matrix->size, &(elem_matrix->index)) );
	for(ii1=0; ii1<element.node_num; ii1++){
		int renum_index = node_renum_index(node, element.node[ii1]);
		for(ii2=0; ii2<dim; ii2++){
			elem_matrix->index[ii1*dim+ii2] = renum_index * dim + ii2;
		}
	}

	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( set_L_matrix_shell(element.type, node_location, L) );
	RC_TRY( local_node_location (element.node_num, L, node_location) );
	RC_TRY( get_thickness(element, physical, &thickness) );
	RC_TRY( Gauss_const_default(element.type, g_points, g_weights,
	                            &g_num_point) );
	/* 曲げ */
	for(ii1=0; ii1<g_num_point; ii1++){
		double det_J;
		double factor;
		double **dN_temp = ws.mat_3_N[0];
		double **Jacobi = ws.mat_3_3[0];
		double **inverse_J = ws.mat_3_3[1];
		double **dNw = ws.mat_3_3N[0];
		double **dN = ws.mat_3_N[1];
		double **dNxy = ws.mat_3_N[2];
		double **B_matrix = ws.mat_6_3N[0];
		double **D_matrix = ws.mat_6_6[0];
		double **BtDB = ws.mat_3N_3N[1];

		RC_TRY( set_dN_matrix(element.type, g_points[ii1], dN_temp) );
		mul_matrix(2, element.node_num, 2, dN_temp, node_location, Jacobi);
		det_J = determinant(2, Jacobi);
		if(det_J < 0.0) return(CAL_ERROR_RC);
		RC_TRY( inverse_matrix(2, Jacobi, inverse_J) );
		RC_TRY( set_dN_matrix4shell(element.type, g_points[ii1], node_location,
		                            dNw, dN) );
		mul_matrix(2, 2, 3*element.node_num, inverse_J, dN, dNxy);
		RC_TRY( set_B_matrix_bend(element.node_num, dNxy, B_matrix) );
		RC_TRY( get_D_matrix4shell(ANA_PLANE_BEND, element, material, physical,
		                           D_matrix) );
		mul_matrix_AtBA(3, 3*element.node_num, B_matrix, D_matrix, BtDB);

		factor = det_J * g_weights[ii1]
		       * ((thickness*thickness*thickness) / 12.0);
		for(ii2=0; ii2<element.node_num; ii2++){
			for(ii3=0; ii3<element.node_num; ii3++){
				for(ii4=0; ii4<3; ii4++){
					for(ii5=0; ii5<3; ii5++){
						elem_matrix->matrix[6*ii2+2+ii4][6*ii3+2+ii5]
						+= factor*BtDB[3*ii2+ii4][3*ii3+ii5];
					}
				}
			}
		}
	}

	/* 剪断 */
	k = 5.0 / 6.0;
	beta = 0.2;
	length = element_length(element.node_num, node_location);
	if( (thickness / length) < beta ){
		k *= (thickness/(beta*length)) * (thickness/(beta*length));
	}
	for(ii1=0; ii1<g_num_point; ii1++){
		double det_J;
		double factor;
		double Nw[3*MAX_NODE];
		double N[3*MAX_NODE];
		double **dN_temp = ws.mat_3_N[0];
		double **Jacobi = ws.mat_3_3[0];
		double **inverse_J = ws.mat_3_3[1];
		double **dNw = ws.mat_3_3N[0];
		double **dNw_xy = ws.mat_3_3N[1];
		double **dN = ws.mat_3_N[1];
		double **B_matrix = ws.mat_6_3N[0];
		double **D_matrix = ws.mat_6_6[0];
		double **BtDB = ws.mat_3N_3N[1];

		RC_TRY( set_dN_matrix(element.type, g_points[ii1], dN_temp) );
		mul_matrix(2, element.node_num, 2, dN_temp, node_location, Jacobi);
		det_J = determinant(2, Jacobi);
		if(det_J < 0.0) return(CAL_ERROR_RC);
		RC_TRY( inverse_matrix(2, Jacobi, inverse_J) );
		RC_TRY( set_N_vector4shell(element.type, g_points[ii1], node_location,
		                           Nw, N) );
		RC_TRY( set_dN_matrix4shell(element.type, g_points[ii1], node_location,
		                            dNw, dN) );
		mul_matrix(2, 2, 3*element.node_num, inverse_J, dNw, dNw_xy);
		RC_TRY( set_B_matrix_shear(element.type, element.node_num, dNw_xy, N,
		                           B_matrix) );
		RC_TRY( get_D_matrix4shell(ANA_PLANE_SHEAR, element, material, physical,
		                           D_matrix) );
		mul_matrix_AtBA(2, 3*element.node_num, B_matrix, D_matrix, BtDB);


		factor = det_J * g_weights[ii1] * k * thickness;
		for(ii2=0; ii2<element.node_num; ii2++){
			for(ii3=0; ii3<element.node_num; ii3++){
				for(ii4=0; ii4<3; ii4++){
					for(ii5=0; ii5<3; ii5++){
						elem_matrix->matrix[6*ii2+2+ii4][6*ii3+2+ii5]
						+= factor*BtDB[3*ii2+ii4][3*ii3+ii5];
					}
				}
			}
		}
	}

	/* 膜 */
	for(ii1=0; ii1<g_num_point; ii1++){
		double det_J;
		double factor;
		double **dN_temp = ws.mat_3_N[0];
		double **dNxy_temp = ws.mat_3_N[1];
		double **Jacobi = ws.mat_3_3[0];
		double **inverse_J = ws.mat_3_3[1];
		double **B_matrix = ws.mat_6_3N[0];
		double **D_matrix = ws.mat_6_6[0];
		double **BtDB = ws.mat_3N_3N[1];

		RC_TRY( set_dN_matrix(element.type, g_points[ii1], dN_temp) );
		mul_matrix(2, element.node_num, 2, dN_temp, node_location, Jacobi);
		det_J = determinant(2, Jacobi);
		if(det_J < 0.0) return(CAL_ERROR_RC);
		RC_TRY( inverse_matrix(2, Jacobi, inverse_J) );
		mul_matrix(2, 2, 3*element.node_num, inverse_J, dN_temp, dNxy_temp);
		RC_TRY( set_B_matrix_membrane(element.node_num, dNxy_temp, B_matrix) );
		RC_TRY( get_D_matrix4shell(ANA_PLANE_STRESS, element, material,
		                           physical, D_matrix) );
		mul_matrix_AtBA(3, 3*element.node_num, B_matrix, D_matrix, BtDB);

		factor = det_J * g_weights[ii1] * thickness;
		for(ii2=0; ii2<element.node_num; ii2++){
			for(ii3=0; ii3<element.node_num; ii3++){
				for(ii4=0; ii4<2; ii4++){
					for(ii5=0; ii5<2; ii5++){
						elem_matrix->matrix[6*ii2+ii4][6*ii3+ii5]
						+= factor*BtDB[2*ii2+ii4][2*ii3+ii5];
					}
				}
			}
		}
	}

	/* ドリリング仮想剛性 */
	RC_TRY( make_drilling_matrix(element, node, material, physical, thickness,
	                             drilling_matrix) );
	for(ii1=0; ii1<element.node_num; ii1++){
		for(ii2=0; ii2<element.node_num; ii2++){
			elem_matrix->matrix[6*ii1+5][6*ii2+5] += drilling_matrix[ii1][ii2];
		}
	}

	/* 局所座標系 -> 全体座標系 */
	RC_TRY( mul_matrix33_LtAL(elem_matrix->size, elem_matrix->matrix, L) );

	return(NORMAL_RC);
}


static RC
make_drilling_matrix (ELEMENT element, NODE_ARRAY node,
                      MATERIAL_PROP_ARRAY material,
                      PHYSICAL_PROP_ARRAY physical,
                      double thickness,
                      double **drilling_matrix)
{
	int ii1, ii2;
	int m_index;
	double drilling_factor;
	double drilling_value;
	double alpha = 0.03;
	MATERIAL_PROP temp_mat;


	m_index = search_material_prop_element(material, physical, &element);
	if(m_index < 0) return(SEEK_ERROR_RC);
	temp_mat = material.array[m_index];
	temp_mat.ana_type = ANA_PLANE_STRESS;
	RC_TRY( fill_material(&temp_mat) );

	RC_TRY( element_volume(&element, node) );

	drilling_factor = alpha * temp_mat.E * thickness * element.volume;
	drilling_value = (-1.0 / (double)(element.node_num - 1)) * drilling_factor;
	for(ii1=0; ii1<element.node_num; ii1++){
		for(ii2=ii1+1; ii2<element.node_num; ii2++){
			drilling_matrix[ii1][ii2] = drilling_matrix[ii2][ii1]
			                          = drilling_value;
		}
		drilling_matrix[ii1][ii1] = drilling_factor;
	}

	return(NORMAL_RC);
}


static RC
set_B_matrix_bend (int node_num, double **dNxy, double **B_matrix)
{
	int ii1;
	int dim = 3;


	for(ii1=0;ii1<node_num;ii1++){
		B_matrix[0][dim*ii1]   = 0.0;
		B_matrix[0][dim*ii1+1] = 0.0;
		B_matrix[0][dim*ii1+2] = -dNxy[0][ii1];
		B_matrix[1][dim*ii1]   = 0.0;
		B_matrix[1][dim*ii1+1] = dNxy[1][ii1];
		B_matrix[1][dim*ii1+2] = 0.0;
		B_matrix[2][dim*ii1]   = 0.0;
		B_matrix[2][dim*ii1+1] = dNxy[0][ii1];
		B_matrix[2][dim*ii1+2] = -dNxy[1][ii1];
	}

	return(NORMAL_RC);
}


static RC
set_B_matrix_shear (ELEM_TYPE type, int node_num, double **dNw_xy, double *N,
                    double **B_matrix)
{
	int ii1;
	int dim = 3;


	switch(type){
	case ELEM_TRI1:
	case ELEM_QUAD1:
		for(ii1=0; ii1<node_num; ii1++){
			B_matrix[0][dim*ii1]   = dNw_xy[0][ii1];
			B_matrix[0][dim*ii1+1] = dNw_xy[0][dim*ii1+1];
			B_matrix[0][dim*ii1+2] = dNw_xy[0][dim*ii1+2] + N[ii1];
			B_matrix[1][dim*ii1]   = dNw_xy[1][ii1];
			B_matrix[1][dim*ii1+1] = dNw_xy[1][dim*ii1+1] - N[ii1];
			B_matrix[1][dim*ii1+2] = dNw_xy[1][dim*ii1+2];
		}
		break;
	case ELEM_TRI2:
	case ELEM_QUAD2:
		for(ii1=0; ii1<node_num; ii1++){
			B_matrix[0][dim*ii1]   = dNw_xy[0][ii1];
			B_matrix[0][dim*ii1+1] = 0.0;
			B_matrix[0][dim*ii1+2] = N[ii1];
			B_matrix[1][dim*ii1]   = dNw_xy[1][ii1];
			B_matrix[1][dim*ii1+1] = -N[ii1];
			B_matrix[1][dim*ii1+2] = 0.0;
		}
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
set_B_matrix_membrane (int node_num, double **dNxy, double **B_matrix)
{
	int ii1;
	int dim = 3;


	for(ii1=0;ii1<node_num;ii1++){
		B_matrix[0][dim*ii1]   = dNxy[0][ii1];
		B_matrix[0][dim*ii1+1] = 0.0;
		B_matrix[1][dim*ii1]   = 0.0;
		B_matrix[1][dim*ii1+1] = dNxy[1][ii1];
		B_matrix[2][dim*ii1]   = dNxy[1][ii1];
		B_matrix[2][dim*ii1+1] = dNxy[0][ii1];
	}

	return(NORMAL_RC);
}


static RC
get_D_matrix4shell (ANALYSIS_TYPE type, ELEMENT element,
                    MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
                    double **D_matrix)
{
	int ii1, ii2;
	int dim;
	int m_index;
	MATERIAL_PROP temp_mat;

	m_index = search_material_prop_element(material, physical, &element);
	if(m_index < 0) return(SEEK_ERROR_RC);

	dim = element_dim(element.type);
	if(dim != 2) return(ARG_ERROR_RC);

	temp_mat = material.array[m_index];

	switch(type){
	case ANA_PLANE_STRESS:
	case ANA_PLANE_BEND:
		temp_mat.ana_type = ANA_PLANE_STRESS;
		RC_TRY( fill_material(&temp_mat) );
		for(ii1=0; ii1<3; ii1++){
			for(ii2=0; ii2<3; ii2++){
				D_matrix[ii1][ii2] = temp_mat.D_matrix[ii1][ii2];
			}
		}
		break;
	case ANA_PLANE_SHEAR:
		temp_mat.ana_type = ANA_PLANE_SHEAR;
		RC_TRY( fill_material(&temp_mat) );
		for(ii1=0; ii1<2; ii1++){
			for(ii2=0; ii2<2; ii2++){
				D_matrix[ii1][ii2] = 0.0;
			}
			D_matrix[ii1][ii1] = (temp_mat.E) / (2.0*(1.0 + temp_mat.nu));
		}
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


/* 平面シェル要素用内挿関数 */
static RC
set_N_vector4shell (ELEM_TYPE type, const double *coord,
                    double **node_location, double *Nw, double *N)
{
	if( (coord == NULL) || (Nw == NULL) || (N == NULL) ) return(ARG_ERROR_RC);

	switch(type){
	case ELEM_TRI1:
		RC_TRY( N_tri1_shell(coord, node_location, Nw, N) );
		break;
	case ELEM_TRI2:
		RC_TRY( log_printf(5, "Warning Locking!!\n") );
		RC_TRY( N_tri2_shell(coord, Nw, N) );
		break;
	case ELEM_QUAD1:
		RC_TRY( N_quad1_shell(coord, node_location, Nw, N) );
		break;
	case ELEM_QUAD2:
		RC_TRY( N_quad2_shell(coord, node_location, Nw, N) );
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
N_tri1_shell (const double *l1_l2_l3, double **xy, double *Nw, double *N)
{
	int ii1;
	int i, j, k;
	const double *L = l1_l2_l3;


	for(ii1=0; ii1<3; ii1++){
		i = ii1;
		j = ii1 + 1;
		if(j >= 3) j -= 3;
		k = ii1 + 2;
		if(k >= 3) k -= 3;

		Nw[3*ii1] = L[ii1];
		Nw[3*ii1+1] = -(L[i]*L[j]*(xy[i][1] - xy[j][1])/2.0
		            + L[i]*L[k]*(xy[i][1] - xy[k][1])/2.0);
		Nw[3*ii1+2] = L[i]*L[j]*(xy[i][0] - xy[j][0])/2.0
		            + L[i]*L[k]*(xy[i][0] - xy[k][0])/2.0;
		N[ii1] = L[ii1];
	}

	return(NORMAL_RC);
}


static RC
N_tri2_shell (const double *l1_l2_l3, double *Nw, double *N)
{
	int ii1;
	const double *L = l1_l2_l3;
	double P[6];


	P[0] = L[0] * (2.0 * L[0] - 1.0);
	P[1] = L[1] * (2.0 * L[1] - 1.0);
	P[2] = L[2] * (2.0 * L[2] - 1.0);
	P[3] = 4.0 * L[0] * L[1];
	P[4] = 4.0 * L[1] * L[2];
	P[5] = 4.0 * L[2] * L[0];

	for(ii1=0; ii1<6; ii1++){
		Nw[ii1] = N[ii1] = P[ii1];
	}

	return(NORMAL_RC);
}


#if 0
static RC
N_tri2_shell_ver2 (const double *l1_l2_l3, double *Nw, double *N)
{
	int ii1;
	double xi = l1_l2_l3[0];
	double eta = l1_l2_l3[1];
	double m_xi = 1.0 - xi;
	double m_eta = 1.0 - eta;
	double p_xi = 1.0 + xi;
	double p_eta = 1.0 + eta;
	double m_xi2 = 1.0 - xi * xi;
	double m_eta2 = 1.0 - eta * eta;
	double xi_m = xi - 1.0;
	double eta_m = eta - 1.0;
	double M[9];
	double P[6];


	M[0] = 0.25 * m_xi * m_eta * xi * eta;
	M[1] = 0.25 * p_xi * eta_m * xi * eta;
	M[2] = 0.25 * p_xi * p_eta * xi * eta;
	M[3] = 0.25 * xi_m * p_eta * xi * eta;
	M[4] = 0.5 * eta * eta_m * m_xi2;
	M[5] = 0.5 * xi * p_xi * m_eta2;
	M[6] = 0.5 * eta * p_eta * m_xi2;
	M[7] = 0.5 * xi * xi_m * m_eta2;
	M[8] = m_xi2 * m_eta2;

	P[0] = M[0] - 0.125*M[8];
	P[1] = M[1] - 0.125*M[8];
	P[2] = M[2] + M[3] + M[6];
	P[3] = M[4] + 0.25*M[8];
	P[4] = M[5] + 0.5*M[8];
	P[5] = M[7] + 0.5*M[8];

	for(ii1=0; ii1<6; ii1++){
		Nw[ii1] = N[ii1] = P[ii1];
	}

	return(NORMAL_RC);
}
#endif


static RC
N_quad1_shell (const double *xi_eta, double **xy, double *Nw, double *N)
{
	int ii1;
	int i, j, k, m;
	int sign;
	double xi = xi_eta[0];
	double eta = xi_eta[1];
	double m_xi = 1.0 - xi;
	double m_eta = 1.0 - eta;
	double p_xi = 1.0 + xi;
	double p_eta = 1.0 + eta;
	double m_xi2 = 1.0 - xi * xi;
	double m_eta2 = 1.0 - eta * eta;
	double P[9];


	P[0] = 0.25 * m_xi * m_eta;
	P[1] = 0.25 * p_xi * m_eta;
	P[2] = 0.25 * p_xi * p_eta;
	P[3] = 0.25 * m_xi * p_eta;
	P[4] = m_xi2 * m_eta;
	P[5] = p_xi * m_eta2;
	P[6] = m_xi2 * p_eta;
	P[7] = m_xi * m_eta2;
	P[8] = m_xi2 * m_eta2;

	for(ii1=0; ii1<4; ii1++){
		i = ii1;
		j = ii1 + 1;
		if(j >= 4) j -= 4;
		k = ii1 + 2;
		if(k >= 4) k -= 4;
		m = ii1 + 3;
		if(m >= 4) m -= 4;
		sign = 1;
		if(ii1%2 == 1) sign = -1;

		Nw[3*ii1] = P[ii1];
		Nw[3*ii1+1] = -(P[i+4]*(xy[i][1] - xy[j][1])/16.0
		            + P[m+4]*(xy[i][1] - xy[m][1])/16.0 + sign
		            * P[8]*(xy[1][1] + xy[3][1] - xy[0][1] - xy[2][1])/32.0);
		Nw[3*ii1+2] = P[i+4]*(xy[i][0] - xy[j][0])/16.0
		            + P[m+4]*(xy[i][0] - xy[m][0])/16.0 + sign
		            * P[8]*(xy[1][0] + xy[3][0] - xy[0][0] - xy[2][0])/32.0;
		N[ii1] = P[ii1];
	}

	return(NORMAL_RC);
}


static RC
N_quad2_shell (const double *xi_eta, double **xy, double *Nw, double *N)
{
	int ii1;
	int i, j, k, m;
	double xi = xi_eta[0];
	double eta = xi_eta[1];
	double m_xi = 1.0 - xi;
	double m_eta = 1.0 - eta;
	double p_xi = 1.0 + xi;
	double p_eta = 1.0 + eta;
	double m_xi2 = 1.0 - xi * xi;
	double m_eta2 = 1.0 - eta * eta;
	double xi_m = xi - 1.0;
	double eta_m = eta - 1.0;
	double M[9];
	double D[4];
	double Di_Dk;


	M[0] = 0.25 * m_xi * m_eta * xi * eta;
	M[1] = 0.25 * p_xi * eta_m * xi * eta;
	M[2] = 0.25 * p_xi * p_eta * xi * eta;
	M[3] = 0.25 * xi_m * p_eta * xi * eta;
	M[4] = 0.5 * eta * eta_m * m_xi2;
	M[5] = 0.5 * xi * p_xi * m_eta2;
	M[6] = 0.5 * eta * p_eta * m_xi2;
	M[7] = 0.5 * xi * xi_m * m_eta2;
	M[8] = m_xi2 * m_eta2;

	for(ii1=0; ii1<4; ii1++){
		i = ii1;
		j = ii1 + 1;
		if(j >= 4) j -= 4;
		m = ii1 + 3;
		if(m >= 4) m -= 4;

		D[i] = (xy[j][0] - xy[i][0])*(xy[m][1] - xy[i][1])
		     - (xy[m][0] - xy[i][0])*(xy[j][1] - xy[i][1]);
	}

	for(ii1=0; ii1<4; ii1++){
		i = ii1;
		j = ii1 + 1;
		if(j >= 4) j -= 4;
		k = ii1 + 2;
		if(k >= 4) k -= 4;
		m = ii1 + 3;
		if(m >= 4) m -= 4;

		Di_Dk = D[i] + D[k];
		Nw[ii1] = N[ii1] = M[ii1] - (0.25 + (D[k] - D[i])/(8.0*Di_Dk))*M[8];
		Nw[ii1+4] = N[ii1+4] = M[ii1+4]
		                     + (0.5 + (D[m] - D[i])/(4.0*Di_Dk))*M[8];
	}

	return(NORMAL_RC);
}


/* 平面シェル要素用内挿関数の偏導関数 */
static RC
set_dN_matrix4shell (ELEM_TYPE type, const double *coord,
                     double **node_location, double **dNw, double **dN)
{
	if(coord == NULL) return(ARG_ERROR_RC);
	if(node_location == NULL) return(ARG_ERROR_RC);

	switch(type){
	case ELEM_TRI1:
		RC_TRY( dN_tri1_shell(coord, node_location, dNw, dN) );
		break;
	case ELEM_TRI2:
		RC_TRY( log_printf(5, "Warning Locking!!\n") );
		RC_TRY( dN_tri2_shell(coord, dNw, dN) );
		break;
	case ELEM_QUAD1:
		RC_TRY( dN_quad1_shell(coord, node_location, dNw, dN) );
		break;
	case ELEM_QUAD2:
		RC_TRY( dN_quad2_shell(coord, node_location, dNw, dN) );
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}


	return(NORMAL_RC);
}


static RC
dN_tri1_shell (const double *l1_l2_l3, double **xy, double **dNw, double **dN)
{
	int ii1, ii2;
	const double *L = l1_l2_l3;
	double dP[3][9];
	double const_matrix[2][3] = {{ 1.0,  0.0, -1.0},
	                             { 0.0,  1.0, -1.0}};
	double temp_matrix[3][3] = {{ 1.0,  0.0, 0.0},
	                            { 0.0,  1.0, 0.0},
	                            { 0.0,  0.0, 1.0}};


	dP[0][0] = 1.0;
	dP[0][1] = -(L[1]*(xy[0][1] - xy[1][1])/2.0
	         + L[2]*(xy[0][1] - xy[2][1])/2.0);
	dP[0][2] = L[1]*(xy[0][0] - xy[1][0])/2.0
	         + L[2]*(xy[0][0] - xy[2][0])/2.0;
	dP[0][3] = 0.0;
	dP[0][4] = -(L[1]*(xy[1][1] - xy[0][1])/2.0);
	dP[0][5] = L[1]*(xy[1][0] - xy[0][0])/2.0;
	dP[0][6] = 0.0;
	dP[0][7] = -(L[2]*(xy[2][1] - xy[0][1])/2.0);
	dP[0][8] = L[2]*(xy[2][0] - xy[0][0])/2.0;

	dP[1][0] = 0.0;
	dP[1][1] = -(L[0]*(xy[0][1] - xy[1][1])/2.0);
	dP[1][2] = L[0]*(xy[0][0] - xy[1][0])/2.0;
	dP[1][3] = 1.0;
	dP[1][4] = -(L[2]*(xy[1][1] - xy[2][1])/2.0
	         + L[0]*(xy[1][1] - xy[0][1])/2.0);
	dP[1][5] = L[2]*(xy[1][0] - xy[2][0])/2.0
	         + L[0]*(xy[1][0] - xy[0][0])/2.0;
	dP[1][6] = 0.0;
	dP[1][7] = -(L[2]*(xy[2][1] - xy[1][1])/2.0);
	dP[1][8] = L[2]*(xy[2][0] - xy[1][0])/2.0;

	dP[2][0] = 0.0;
	dP[2][1] = -(L[0]*(xy[0][1] - xy[2][1])/2.0);
	dP[2][2] = L[0]*(xy[0][0] - xy[2][0])/2.0;
	dP[2][3] = 0.0;
	dP[2][4] = -(L[1]*(xy[1][1] - xy[2][1])/2.0);
	dP[2][5] = L[1]*(xy[1][0] - xy[2][0])/2.0;
	dP[2][6] = 1.0;
	dP[2][7] = -(L[0]*(xy[2][1] - xy[0][1])/2.0
	         + L[1]*(xy[2][1] - xy[1][1])/2.0);
	dP[2][8] = L[0]*(xy[2][0] - xy[0][0])/2.0
	         + L[1]*(xy[2][0] - xy[1][0])/2.0;

	for(ii1=0; ii1<2; ii1++){
		for(ii2=0; ii2<9; ii2++){
			dNw[ii1][ii2] = const_matrix[ii1][0]*dN[0][ii2]
			              + const_matrix[ii1][1]*dN[1][ii2]
			              + const_matrix[ii1][2]*dN[2][ii2];
		}
	}

	for(ii1=0; ii1<2; ii1++){
		for(ii2=0; ii2<3; ii2++){
			dN[ii1][ii2] = const_matrix[ii1][0] * temp_matrix[0][ii2]
			             + const_matrix[ii1][1] * temp_matrix[1][ii2]
			             + const_matrix[ii1][2] * temp_matrix[2][ii2];
		}
	}

	return(NORMAL_RC);
}


static RC
dN_tri2_shell (const double *l1_l2_l3, double **dNw, double **dN)
{
	int ii1, ii2;
	double l1 = l1_l2_l3[0];
	double l2 = l1_l2_l3[1];
	double l3 = l1_l2_l3[2];
	double const_matrix[2][3] = {{ 1.0,  0.0, -1.0},
	                             { 0.0,  1.0, -1.0}};
	double tmp_matrix[3][6];


	tmp_matrix[0][0] = 4.0 * l1 - 1.0;
	tmp_matrix[0][1] = 0.0;
	tmp_matrix[0][2] = 0.0;
	tmp_matrix[0][3] = 4.0 * l2;
	tmp_matrix[0][4] = 0.0;
	tmp_matrix[0][5] = 4.0 * l3;

	tmp_matrix[1][0] = 0.0;
	tmp_matrix[1][1] = 4.0 * l2 - 1.0;
	tmp_matrix[1][2] = 0.0;
	tmp_matrix[1][3] = 4.0 * l1;
	tmp_matrix[1][4] = 4.0 * l3;
	tmp_matrix[1][5] = 0.0;

	tmp_matrix[2][0] = 0.0;
	tmp_matrix[2][1] = 0.0;
	tmp_matrix[2][2] = 4.0 * l3 - 1.0;
	tmp_matrix[2][3] = 0.0;
	tmp_matrix[2][4] = 4.0 * l2;
	tmp_matrix[2][5] = 4.0 * l1;

	for(ii1=0; ii1<2; ii1++){
		for(ii2=0; ii2<6; ii2++){
			dNw[ii1][ii2] = dN[ii1][ii2]
			              = const_matrix[ii1][0] * tmp_matrix[0][ii2]
			              + const_matrix[ii1][1] * tmp_matrix[1][ii2]
			              + const_matrix[ii1][2] * tmp_matrix[2][ii2];
		}
	}

	return(NORMAL_RC);
}


#if 0
static RC
dN_tri2_shell_ver2 (const double *l1_l2_l3, double **dNw, double **dN)
{
	int ii1, ii2;
	double xi = l1_l2_l3[0];
	double eta = l1_l2_l3[1];
	double m_xi = 1.0 - xi;
	double p_xi = 1.0 + xi;
	double m_eta = 1.0 - eta;
	double p_eta = 1.0 + eta;
	double m_xi2 = 1.0 - xi * xi;
	double m_eta2 = 1.0 - eta * eta;
	double xi_m = xi - 1.0;
	double eta_m = eta - 1.0;
	double M[2][9];
	double dP[2][6];


	M[0][0] = 0.25 * (-xi + m_xi) * m_eta * eta;
	M[0][1] = 0.25 * (xi + p_xi) * eta_m * eta;
	M[0][2] = 0.25 * (xi + p_xi) * p_eta * eta;
	M[0][3] = 0.25 * (xi + xi_m) * p_eta * eta;
	M[0][4] = 0.5 * eta * eta_m * (-2.0*xi);
	M[0][5] = 0.5 * (p_xi + xi) * m_eta2;
	M[0][6] = 0.5 * eta * p_eta * (-2.0*xi);
	M[0][7] = 0.5 * (xi_m + xi) * m_eta2;
	M[0][8] = (-2.0*xi) * m_eta2;

	M[1][0] = 0.25 * m_xi * xi * (-eta + m_eta);
	M[1][1] = 0.25 * p_xi * xi * (eta + eta_m);
	M[1][2] = 0.25 * p_xi * xi * (eta + p_eta);
	M[1][3] = 0.25 * xi_m * xi * (eta + p_eta);
	M[1][4] = 0.5 * (eta_m + eta) * m_xi2;
	M[1][5] = 0.5 * xi * p_xi * (-2.0*eta);
	M[1][6] = 0.5 * (p_eta + eta) * m_xi2;
	M[1][7] = 0.5 * xi * xi_m * (-2.0*eta);
	M[1][8] = m_xi2 * (-2.0*eta);

	dP[0][0] = M[0][0] - 0.125*M[0][8];
	dP[0][1] = M[0][1] - 0.125*M[0][8];
	dP[0][2] = M[0][2] + M[0][3] + M[0][6];
	dP[0][3] = M[0][4] + 0.25*M[0][8];
	dP[0][4] = M[0][5] + 0.5*M[0][8];
	dP[0][5] = M[0][7] + 0.5*M[0][8];

	dP[1][0] = M[1][0] - 0.125*M[1][8];
	dP[1][1] = M[1][1] - 0.125*M[1][8];
	dP[1][2] = M[1][2] + M[1][3] + M[1][6];
	dP[1][3] = M[1][4] + 0.25*M[1][8];
	dP[1][4] = M[1][5] + 0.5*M[1][8];
	dP[1][5] = M[1][7] + 0.5*M[1][8];

	for(ii1=0; ii1<2; ii1++){
		for(ii2=0; ii2<6; ii2++){
			dNw[ii1][ii2] = dN[ii1][ii2] = dP[ii1][ii2];
		}
	}

	return(NORMAL_RC);
}
#endif


static RC
dN_quad1_shell (const double *xi_eta, double **xy, double **dNw, double **dN)
{
	int ii1;
	double xi = xi_eta[0];
	double eta = xi_eta[1];
	double P[2][9];

	P[0][0] = -0.25*(1.0 - eta);
	P[0][1] =  0.25*(1.0 - eta);
	P[0][2] =  0.25*(1.0 + eta);
	P[0][3] = -0.25*(1.0 + eta);
	P[0][4] = -2.0*xi*(1.0 - eta);
	P[0][5] =  (1.0 - eta*eta);
	P[0][6] = -2.0*xi*(1.0 + eta);
	P[0][7] = -(1.0 - eta*eta);
	P[0][8] = -2.0*xi*(1.0 - eta*eta);

	P[1][0] = -0.25*(1.0 - xi);
	P[1][1] = -0.25*(1.0 + xi);
	P[1][2] =  0.25*(1.0 + xi);
	P[1][3] =  0.25*(1.0 - xi);
	P[1][4] = -(1.0 - xi*xi);
	P[1][5] = -2.0*eta*(1.0 + xi);
	P[1][6] =  (1.0 - xi*xi);
	P[1][7] = -2.0*eta*(1.0 - xi);
	P[1][8] = -2.0*eta*(1.0 - xi*xi);

	for(ii1=0; ii1<2; ii1++){
		RC_TRY( dN_quad1_shell_sub(P[ii1], xy, dNw[ii1], dN[ii1]) );
	}

	return(NORMAL_RC);
}


static RC
dN_quad1_shell_sub (const double P[], double **xy, double *dNw, double *dN)
{
	int ii1;
	int i, j, k, m;
	int sign;


	for(ii1=0; ii1<4; ii1++){
		i = ii1;
		j = ii1 + 1;
		if(j >= 4) j -= 4;
		k = ii1 + 2;
		if(k >= 4) k -= 4;
		m = ii1 + 3;
		if(m >= 4) m -= 4;
		sign = 1;
		if(ii1%2 == 1) sign = -1;

		dNw[3*ii1] = P[ii1];
		dNw[3*ii1+1] = -(P[i+4]*(xy[i][1] - xy[j][1])/16.0
		            + P[m+4]*(xy[i][1] - xy[m][1])/16.0
		            + sign*P[8]*(xy[1][1] + xy[3][1]
		            - xy[0][1] - xy[2][1])/32.0);
		dNw[3*ii1+2] = P[i+4]*(xy[i][0] - xy[j][0])/16.0
		            + P[m+4]*(xy[i][0] - xy[m][0])/16.0
		            + sign*P[8]*(xy[1][0] + xy[3][0]
		            - xy[0][0] - xy[2][0])/32.0;
		dN[ii1] = P[ii1];
	}

	return(NORMAL_RC);
}


static RC
dN_quad2_shell (const double *xi_eta, double **xy, double **dNw, double **dN)
{
	int ii1, ii2;
	int i, j, k, m;
	double Di_Dk;
	double D[4];
	double dP[2][8];
	double M[2][9];


	RC_TRY( dN_quad2_shell_sub(xi_eta, xy, M, D) );
	for(ii1=0; ii1<4; ii1++){
		i = ii1;
		j = ii1 + 1;
		if(j >= 4) j -= 4;
		k = ii1 + 2;
		if(k >= 4) k -= 4;
		m = ii1 + 3;
		if(m >= 4) m -= 4;

		Di_Dk = D[i] + D[k];
		Di_Dk += ABS_TOL * SIGN(Di_Dk);

		dP[0][ii1] = M[0][ii1] - (0.25 + (D[k] - D[i])/(8.0*Di_Dk))*M[0][8];
		dP[0][ii1+4] = M[0][ii1+4] + (0.5 + (D[m] - D[i])/(4.0*Di_Dk))*M[0][8];

		dP[1][ii1] = M[1][ii1] - (0.25 + (D[k] - D[i])/(8.0*Di_Dk))*M[1][8];
		dP[1][ii1+4] = M[1][ii1+4] + (0.5 + (D[m] - D[i])/(4.0*Di_Dk))*M[1][8];
	}

	for(ii1=0; ii1<2; ii1++){
		for(ii2=0; ii2<8; ii2++){
			dNw[ii1][ii2] = dN[ii1][ii2] = dP[ii1][ii2];
		}
	}

	return(NORMAL_RC);
}


static RC
dN_quad2_shell_sub (const double *xi_eta, double **xy, double M[][9],
                    double D[])
{
	int ii1;
	int i, j, m;
	double xi = xi_eta[0];
	double eta = xi_eta[1];
	double m_xi = 1.0 - xi;
	double p_xi = 1.0 + xi;
	double m_eta = 1.0 - eta;
	double p_eta = 1.0 + eta;
	double m_xi2 = 1.0 - xi * xi;
	double m_eta2 = 1.0 - eta * eta;
	double xi_m = xi - 1.0;
	double eta_m = eta - 1.0;


	M[0][0] = 0.25 * (-xi + m_xi) * m_eta * eta;
	M[0][1] = 0.25 * (xi + p_xi) * eta_m * eta;
	M[0][2] = 0.25 * (xi + p_xi) * p_eta * eta;
	M[0][3] = 0.25 * (xi + xi_m) * p_eta * eta;
	M[0][4] = 0.5 * eta * eta_m * (-2.0*xi);
	M[0][5] = 0.5 * (p_xi + xi) * m_eta2;
	M[0][6] = 0.5 * eta * p_eta * (-2.0*xi);
	M[0][7] = 0.5 * (xi_m + xi) * m_eta2;
	M[0][8] = (-2.0*xi) * m_eta2;

	M[1][0] = 0.25 * m_xi * xi * (-eta + m_eta);
	M[1][1] = 0.25 * p_xi * xi * (eta + eta_m);
	M[1][2] = 0.25 * p_xi * xi * (eta + p_eta);
	M[1][3] = 0.25 * xi_m * xi * (eta + p_eta);
	M[1][4] = 0.5 * (eta_m + eta) * m_xi2;
	M[1][5] = 0.5 * xi * p_xi * (-2.0*eta);
	M[1][6] = 0.5 * (p_eta + eta) * m_xi2;
	M[1][7] = 0.5 * xi * xi_m * (-2.0*eta);
	M[1][8] = m_xi2 * (-2.0*eta);

	for(ii1=0; ii1<4; ii1++){
		i = ii1;
		j = ii1 + 1;
		if(j >= 4) j -= 4;
		m = ii1 + 3;
		if(m >= 4) m -= 4;

		D[i] = (xy[j][0] - xy[i][0])*(xy[m][1] - xy[i][1])
		     - (xy[m][0] - xy[i][0])*(xy[j][1] - xy[i][1]);
	}

	return(NORMAL_RC);
}


/* 座標変換マトリクス L[3][3] */
static RC
set_L_matrix_shell (ELEM_TYPE type, double **node_location, double L[][3])
{
	VECT3D temp_v1, temp_v2;
	VECT3D ex, ey, ez;

	switch(type){
		case ELEM_TRI1:
		case ELEM_TRI2:
			temp_v1.x = node_location[1][0] - node_location[0][0];
			temp_v1.y = node_location[1][1] - node_location[0][1];
			temp_v1.z = node_location[1][2] - node_location[0][2];
			if(nearly_eq(abs_vect3d(temp_v1), 0.0)) return(CAL_ERROR_RC);
			ex = unit_vect3d(temp_v1);

			temp_v2.x = node_location[2][0] - node_location[0][0];
			temp_v2.y = node_location[2][1] - node_location[0][1];
			temp_v2.z = node_location[2][2] - node_location[0][2];
			temp_v2 = gs_ortho3d(temp_v2, ex);
			if(nearly_eq(abs_vect3d(temp_v2), 0.0)) return(CAL_ERROR_RC);
			ey = unit_vect3d(temp_v2);

			ez = unit_vect3d( outer_product3d(ex, ey) );
			break;
		case ELEM_QUAD1:
		case ELEM_QUAD2:
			temp_v1.x = node_location[1][0] - node_location[3][0];
			temp_v1.y = node_location[1][1] - node_location[3][1];
			temp_v1.z = node_location[1][2] - node_location[3][2];

			temp_v2.x = node_location[2][0] - node_location[0][0];
			temp_v2.y = node_location[2][1] - node_location[0][1];
			temp_v2.z = node_location[2][2] - node_location[0][2];

			if(nearly_eq(abs_vect3d(temp_v1), 0.0)) return(CAL_ERROR_RC);
			if(nearly_eq(abs_vect3d(temp_v2), 0.0)) return(CAL_ERROR_RC);

			ex = unit_vect3d( add_vect3d(unit_vect3d(temp_v2),
			                             unit_vect3d(temp_v1)) );
			ey = unit_vect3d( sub_vect3d(unit_vect3d(temp_v2),
			                             unit_vect3d(temp_v1)) );
			ey = unit_vect3d( gs_ortho3d(ey, ex) );
			ez = unit_vect3d( outer_product3d(ex, ey) );
			break;
		default:
			return(IMPLEMENT_ERROR_RC);
	}

	L[0][0] = ex.x;
	L[0][1] = ex.y;
	L[0][2] = ex.z;
	L[1][0] = ey.x;
	L[1][1] = ey.y;
	L[1][2] = ey.z;
	L[2][0] = ez.x;
	L[2][1] = ez.y;
	L[2][2] = ez.z;

	return(NORMAL_RC);
}


static RC
local_node_location (int node_num, double L[][3], double **node_location)
{
	int ii1;

	
	for(ii1=0; ii1<node_num; ii1++){
		mul_matrix33_vect(L, node_location[ii1]);
	}

	return(NORMAL_RC);
}


static double
element_length (int node_num, double **node_location)
{
	int ii1, ii2;
	double ret;


	ret = 0.0;
	for(ii1=0; ii1<node_num; ii1++){
		for(ii2=ii1+1; ii2<node_num; ii2++){
			double length;

			length = dist_point3d( array2vect3d(node_location[ii1]),
			                       array2vect3d(node_location[ii2]) );
			if(length > ret) ret = length;
		}
	}

	return(ret);
}




