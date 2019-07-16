/*********************************************************************
 * interpolation.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Thanks to following contributors.
 *  <Tomoyuki TSUBATA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: interpolation.c,v 1.7 2003/07/22 11:44:01 adachi Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include "rc.h"
#include "math_utl.h"
#include "fem_struct.h"


static RC dN_line1(const double *xi, double **dN);
static RC dN_line2(const double *xi, double **dN);
static RC dN_tri1(const double *l1_l2_l3, double **dN);
static RC dN_tri2(const double *l1_l2_l3, double **dN);
static RC dN_quad1(const double *xi_eta, double **dN); 
static RC dN_quad2(const double *xi_eta, double **dN);
static RC dN_tetra1(const double *l1_l2_l3_l4, double **dN);
static RC dN_tetra2(const double *l1_l2_l3_l4, double **dN);
static RC dN_penta1(const double *l1_l2_l3_zeta, double **dN);
static RC dN_penta2(const double *l1_l2_l3_zeta, double **dN);
static RC dN_hexa1(const double *xi_eta_zeta, double **dN);
static RC dN_hexa2(const double *xi_eta_zeta, double **dN);
static RC N_line1(const double *xi, double *N);
static RC N_line2(const double *xi, double *N);
static RC N_tri1(const double *l1_l2_l3, double *N);
static RC N_tri2(const double *l1_l2_l3, double *N);
static RC N_quad1(const double *xi_eta, double *N);
static RC N_quad2(const double *xi_eta, double *N);
static RC N_tetra1(const double *l1_l2_l3_l4, double *N);
static RC N_tetra2(const double *l1_l2_l3_l4, double *N);
static RC N_penta1(const double *l1_l2_l3_zeta, double *N);
static RC N_penta2(const double *l1_l2_l3_zeta, double *N);
static RC N_hexa1(const double *xi_eta_zeta, double *N);
static RC N_hexa2(const double *xi_eta_zeta, double *N);


RC set_local_node_points(FEM_ELEM_TYPE type, double **coords)
{
	switch(type){
	case ELEM_RBAR:
	case ELEM_BEAM1:
	case ELEM_LINE1:
		coords[0][0] = -1.0;
		coords[1][0] = 1.0;
		break;
	case ELEM_BEAM2:
	case ELEM_LINE2:
		coords[0][0] = -1.0;
		coords[1][0] = 1.0;
		coords[2][0] = 0.0;
		break;
	case ELEM_TRI1:
		coords[0][0] = 1.0;
		coords[0][1] = 0.0;
		coords[0][2] = 0.0;
		coords[1][0] = 0.0;
		coords[1][1] = 1.0;
		coords[1][2] = 0.0;
		coords[2][0] = 0.0;
		coords[2][1] = 0.0;
		coords[2][2] = 1.0;
		break;
	case ELEM_TRI2:
		coords[0][0] = 1.0;
		coords[0][1] = 0.0;
		coords[0][2] = 0.0;
		coords[1][0] = 0.0;
		coords[1][1] = 1.0;
		coords[1][2] = 0.0;
		coords[2][0] = 0.0;
		coords[2][1] = 0.0;
		coords[2][2] = 1.0;
		coords[3][0] = 0.0;
		coords[3][1] = 0.5;
		coords[3][2] = 0.5;
		coords[4][0] = 0.5;
		coords[4][1] = 0.0;
		coords[4][2] = 0.5;
		coords[5][0] = 0.5;
		coords[5][1] = 0.5;
		coords[5][2] = 0.0;
		break;
	case ELEM_QUAD1:
		coords[0][0] = -1.0;
		coords[0][1] = -1.0;
		coords[1][0] = 1.0;
		coords[1][1] = -1.0;
		coords[2][0] = 1.0;
		coords[2][1] = 1.0;
		coords[3][0] = -1.0;
		coords[3][1] = 1.0;
		break;
	case ELEM_QUAD2:
		coords[0][0] = -1.0;
		coords[0][1] = -1.0;
		coords[1][0] = 1.0;
		coords[1][1] = -1.0;
		coords[2][0] = 1.0;
		coords[2][1] = 1.0;
		coords[3][0] = -1.0;
		coords[3][1] = 1.0;
		coords[4][0] = 0.0;
		coords[4][1] = -1.0;
		coords[5][0] = 1.0;
		coords[5][1] = 0.0;
		coords[6][0] = 0.0;
		coords[6][1] = 1.0;
		coords[7][0] = -1.0;
		coords[7][1] = 0.0;
		break;
	case ELEM_TETRA1:
		coords[0][0] = 1.0;
		coords[0][1] = 0.0;
		coords[0][2] = 0.0;
		coords[0][3] = 0.0;
		coords[1][0] = 0.0;
		coords[1][1] = 1.0;
		coords[1][2] = 0.0;
		coords[1][3] = 0.0;
		coords[2][0] = 0.0;
		coords[2][1] = 0.0;
		coords[2][2] = 1.0;
		coords[2][3] = 0.0;
		coords[3][0] = 0.0;
		coords[3][1] = 0.0;
		coords[3][2] = 0.0;
		coords[3][3] = 1.0;
		break;
	case ELEM_TETRA2:
		coords[0][0] = 1.0;
		coords[0][1] = 0.0;
		coords[0][2] = 0.0;
		coords[0][3] = 0.0;
		coords[1][0] = 0.0;
		coords[1][1] = 1.0;
		coords[1][2] = 0.0;
		coords[1][3] = 0.0;
		coords[2][0] = 0.0;
		coords[2][1] = 0.0;
		coords[2][2] = 1.0;
		coords[2][3] = 0.0;
		coords[3][0] = 0.0;
		coords[3][1] = 0.0;
		coords[3][2] = 0.0;
		coords[3][3] = 1.0;
		coords[4][0] = 0.0;
		coords[4][1] = 0.5;
		coords[4][2] = 0.5;
		coords[4][3] = 0.0;
		coords[5][0] = 0.5;
		coords[5][1] = 0.0;
		coords[5][2] = 0.5;
		coords[5][3] = 0.0;
		coords[6][0] = 0.5;
		coords[6][1] = 0.5;
		coords[6][2] = 0.0;
		coords[6][3] = 0.0;
		coords[7][0] = 0.5;
		coords[7][1] = 0.0;
		coords[7][2] = 0.0;
		coords[7][3] = 0.5;
		coords[8][0] = 0.0;
		coords[8][1] = 0.5;
		coords[8][2] = 0.0;
		coords[8][3] = 0.5;
		coords[9][0] = 0.0;
		coords[9][1] = 0.0;
		coords[9][2] = 0.5;
		coords[9][3] = 0.5;
		break;
	case ELEM_PENTA1:
		coords[0][0] = 1.0;
		coords[0][1] = 0.0;
		coords[0][2] = 0.0;
		coords[0][3] = -1.0;
		coords[1][0] = 0.0;
		coords[1][1] = 1.0;
		coords[1][2] = 0.0;
		coords[1][3] = -1.0;
		coords[2][0] = 0.0;
		coords[2][1] = 0.0;
		coords[2][2] = 1.0;
		coords[2][3] = -1.0;
		coords[3][0] = 1.0;
		coords[3][1] = 0.0;
		coords[3][2] = 0.0;
		coords[3][3] = 1.0;
		coords[4][0] = 0.0;
		coords[4][1] = 1.0;
		coords[4][2] = 0.0;
		coords[4][3] = 1.0;
		coords[5][0] = 0.0;
		coords[5][1] = 0.0;
		coords[5][2] = 1.0;
		coords[5][3] = 1.0;
		break;
	case ELEM_PENTA2:
		coords[0][0] = 1.0;
		coords[0][1] = 0.0;
		coords[0][2] = 0.0;
		coords[0][3] = -1.0;
		coords[1][0] = 0.0;
		coords[1][1] = 1.0;
		coords[1][2] = 0.0;
		coords[1][3] = -1.0;
		coords[2][0] = 0.0;
		coords[2][1] = 0.0;
		coords[2][2] = 1.0;
		coords[2][3] = -1.0;
		coords[3][0] = 1.0;
		coords[3][1] = 0.0;
		coords[3][2] = 0.0;
		coords[3][3] = 1.0;
		coords[4][0] = 0.0;
		coords[4][1] = 1.0;
		coords[4][2] = 0.0;
		coords[4][3] = 1.0;
		coords[5][0] = 0.0;
		coords[5][1] = 0.0;
		coords[5][2] = 1.0;
		coords[5][3] = 1.0;
		coords[6][0] = 0.0;
		coords[6][1] = 0.5;
		coords[6][2] = 0.5;
		coords[6][3] = -1.0;
		coords[7][0] = 0.5;
		coords[7][1] = 0.0;
		coords[7][2] = 0.5;
		coords[7][3] = -1.0;
		coords[8][0] = 0.5;
		coords[8][1] = 0.5;
		coords[8][2] = 0.0;
		coords[8][3] = -1.0;
		coords[9][0] = 0.0;
		coords[9][1] = 0.5;
		coords[9][2] = 0.5;
		coords[9][3] = 1.0;
		coords[10][0] = 0.5;
		coords[10][1] = 0.0;
		coords[10][2] = 0.5;
		coords[10][3] = 1.0;
		coords[11][0] = 0.5;
		coords[11][1] = 0.5;
		coords[11][2] = 0.0;
		coords[11][3] = 1.0;
		coords[12][0] = 1.0;
		coords[12][1] = 0.0;
		coords[12][2] = 0.0;
		coords[12][3] = 0.0;
		coords[13][0] = 0.0;
		coords[13][1] = 1.0;
		coords[13][2] = 0.0;
		coords[13][3] = 0.0;
		coords[14][0] = 0.0;
		coords[14][1] = 0.0;
		coords[14][2] = 1.0;
		coords[14][3] = 0.0;
		break;
	case ELEM_HEXA1:
		coords[0][0] = -1.0;
		coords[0][1] = -1.0;
		coords[0][2] = -1.0;
		coords[1][0] = 1.0;
		coords[1][1] = -1.0;
		coords[1][2] = -1.0;
		coords[2][0] = 1.0;
		coords[2][1] = 1.0;
		coords[2][2] = -1.0;
		coords[3][0] = -1.0;
		coords[3][1] = 1.0;
		coords[3][2] = -1.0;
		coords[4][0] = -1.0;
		coords[4][1] = -1.0;
		coords[4][2] = 1.0;
		coords[5][0] = 1.0;
		coords[5][1] = -1.0;
		coords[5][2] = 1.0;
		coords[6][0] = 1.0;
		coords[6][1] = 1.0;
		coords[6][2] = 1.0;
		coords[7][0] = -1.0;
		coords[7][1] = 1.0;
		coords[7][2] = 1.0;
		break;
	case ELEM_HEXA2:
		coords[0][0] = -1.0;
		coords[0][1] = -1.0;
		coords[0][2] = -1.0;
		coords[1][0] = 1.0;
		coords[1][1] = -1.0;
		coords[1][2] = -1.0;
		coords[2][0] = 1.0;
		coords[2][1] = 1.0;
		coords[2][2] = -1.0;
		coords[3][0] = -1.0;
		coords[3][1] = 1.0;
		coords[3][2] = -1.0;
		coords[4][0] = -1.0;
		coords[4][1] = -1.0;
		coords[4][2] = 1.0;
		coords[5][0] = 1.0;
		coords[5][1] = -1.0;
		coords[5][2] = 1.0;
		coords[6][0] = 1.0;
		coords[6][1] = 1.0;
		coords[6][2] = 1.0;
		coords[7][0] = -1.0;
		coords[7][1] = 1.0;
		coords[7][2] = 1.0;
		coords[8][0] = 0.0;
		coords[8][1] = -1.0;
		coords[8][2] = -1.0;
		coords[9][0] = 1.0;
		coords[9][1] = 0.0;
		coords[9][2] = -1.0;
		coords[10][0] = 0.0;
		coords[10][1] = 1.0;
		coords[10][2] = -1.0;
		coords[11][0] = -1.0;
		coords[11][1] = 0.0;
		coords[11][2] = -1.0;
		coords[12][0] = 0.0;
		coords[12][1] = -1.0;
		coords[12][2] = 1.0;
		coords[13][0] = 1.0;
		coords[13][1] = 0.0;
		coords[13][2] = 1.0;
		coords[14][0] = 0.0;
		coords[14][1] = 1.0;
		coords[14][2] = 1.0;
		coords[15][0] = -1.0;
		coords[15][1] = 0.0;
		coords[15][2] = 1.0;
		coords[16][0] = -1.0;
		coords[16][1] = -1.0;
		coords[16][2] = 0.0;
		coords[17][0] = 1.0;
		coords[17][1] = -1.0;
		coords[17][2] = 0.0;
		coords[18][0] = 1.0;
		coords[18][1] = 1.0;
		coords[18][2] = 0.0;
		coords[19][0] = -1.0;
		coords[19][1] = 1.0;
		coords[19][2] = 0.0;
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


RC set_center_point(FEM_ELEM_TYPE type, double *coord)
{
	double w;

	switch(type){
	case ELEM_RBAR:
	case ELEM_BEAM1:
	case ELEM_BEAM2:
	case ELEM_LINE1:
	case ELEM_LINE2:
		RC_TRY( Gauss_const_line(1, &coord, &w) );
		break;
	case ELEM_TRI1:
	case ELEM_TRI2:
		RC_TRY( Gauss_const_tri(1, &coord, &w) );
		break;
	case ELEM_QUAD1:
	case ELEM_QUAD2:
		RC_TRY( Gauss_const_quad(1, &coord, &w) );
		break;
	case ELEM_TETRA1:
	case ELEM_TETRA2:
		RC_TRY( Gauss_const_tetra(1, &coord, &w) );
		break;
	case ELEM_PENTA1:
	case ELEM_PENTA2:
		RC_TRY( Gauss_const_penta(1, &coord, &w) );
		break;
	case ELEM_HEXA1:
	case ELEM_HEXA2:
		RC_TRY( Gauss_const_hexa(1, &coord, &w) );
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


RC set_dN_matrix(FEM_ELEM_TYPE type, const double *coord, double **dN)
{
	if( (coord == NULL)||(dN == NULL) ) return(ARG_ERROR_RC);

	switch(type){
	case ELEM_RBAR:
	case ELEM_BEAM1:
	case ELEM_LINE1:
		RC_TRY( dN_line1(coord, dN) );
		break;
	case ELEM_BEAM2:
	case ELEM_LINE2:
		RC_TRY( dN_line2(coord, dN) );
		break;
	case ELEM_TRI1:
		RC_TRY( dN_tri1(coord, dN) );
		break;
	case ELEM_TRI2:
		RC_TRY( dN_tri2(coord, dN) );
		break;
	case ELEM_QUAD1:
		RC_TRY( dN_quad1(coord, dN) );
		break;
	case ELEM_QUAD2:
		RC_TRY( dN_quad2(coord, dN) );
		break;
	case ELEM_TETRA1:
		RC_TRY( dN_tetra1(coord, dN) );
		break;
	case ELEM_TETRA2:
		RC_TRY( dN_tetra2(coord, dN) );
		break;
	case ELEM_PENTA1:
		RC_TRY( dN_penta1(coord, dN) );
		break;
	case ELEM_PENTA2:
		RC_TRY( dN_penta2(coord, dN) );
		break;
	case ELEM_HEXA1:
		RC_TRY( dN_hexa1(coord, dN) );
		break;
	case ELEM_HEXA2:
		RC_TRY( dN_hexa2(coord, dN) );
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


RC set_N_vector(FEM_ELEM_TYPE type, const double *coord, double *N)
{
	if( (coord == NULL)||(N == NULL) ) return(ARG_ERROR_RC);

	switch(type){
	case ELEM_RBAR:
	case ELEM_BEAM1:
	case ELEM_LINE1:
		RC_TRY( N_line1(coord, N) );
		break;
	case ELEM_BEAM2:
	case ELEM_LINE2:
		RC_TRY( N_line2(coord, N) );
		break;
	case ELEM_TRI1:
		RC_TRY( N_tri1(coord, N) );
		break;
	case ELEM_TRI2:
		RC_TRY( N_tri2(coord, N) );
		break;
	case ELEM_QUAD1:
		RC_TRY( N_quad1(coord, N) );
		break;
	case ELEM_QUAD2:
		RC_TRY( N_quad2(coord, N) );
		break;
	case ELEM_TETRA1:
		RC_TRY( N_tetra1(coord, N) );
		break;
	case ELEM_TETRA2:
		RC_TRY( N_tetra2(coord, N) );
		break;
	case ELEM_PENTA1:
		RC_TRY( N_penta1(coord, N) );
		break;
	case ELEM_PENTA2:
		RC_TRY( N_penta2(coord, N) );
		break;
	case ELEM_HEXA1:
		RC_TRY( N_hexa1(coord, N) );
		break;
	case ELEM_HEXA2:
		RC_TRY( N_hexa2(coord, N) );
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC dN_hexa2(const double *xi_eta_zeta, double **dN)
{
	double xi = xi_eta_zeta[0];
	double eta = xi_eta_zeta[1];
	double zeta = xi_eta_zeta[2];
	double m_xi = 1.0 - xi;
	double p_xi = 1.0 + xi;
	double m_eta = 1.0 - eta;
	double p_eta = 1.0 + eta;
	double m_zeta = 1.0 - zeta;
	double p_zeta = 1.0 + zeta;
	double xi2 = 1.0 - xi * xi;
	double eta2 = 1.0 - eta * eta;
	double zeta2 = 1.0 - zeta * zeta;

	dN[0][0] = 0.125 * m_eta * m_zeta * (1 + 2.0 * xi + eta + zeta);
	dN[0][1] = 0.125 * m_eta * m_zeta * (-1 + 2.0 * xi - eta - zeta);
	dN[0][2] = 0.125 * p_eta * m_zeta * (-1 + 2.0 * xi + eta - zeta);
	dN[0][3] = 0.125 * p_eta * m_zeta * (1 + 2.0 * xi - eta + zeta);
	dN[0][4] = 0.125 * m_eta * p_zeta * (1 + 2.0 * xi + eta - zeta);
	dN[0][5] = 0.125 * m_eta * p_zeta * (-1 + 2.0 * xi - eta + zeta);
	dN[0][6] = 0.125 * p_eta * p_zeta * (-1 + 2.0 * xi + eta + zeta);
	dN[0][7] = 0.125 * p_eta * p_zeta * (1 + 2.0 * xi - eta - zeta);
	dN[0][8] = - 0.5 * xi * m_eta * m_zeta;
	dN[0][9] = 0.25 * eta2 * m_zeta;
	dN[0][10] = - 0.5 * xi * p_eta * m_zeta;
	dN[0][11] = - 0.25 * eta2 * m_zeta;
	dN[0][12] = - 0.5 * xi * m_eta * p_zeta;
	dN[0][13] = 0.25 * eta2 * p_zeta;
	dN[0][14] = - 0.5 * xi * p_eta * p_zeta;   
	dN[0][15] = - 0.25 * eta2 * p_zeta;
	/*dN[0][16] = - 0.25 * p_eta * zeta2;*/
	dN[0][16] = - 0.25 * m_eta * zeta2;
	dN[0][17] = 0.25 * m_eta * zeta2;
	dN[0][18] = 0.25 * p_eta * zeta2;
	dN[0][19] = - 0.25 * p_eta * zeta2;

	dN[1][0] = 0.125 * m_zeta * m_xi * (1 + xi + 2.0 * eta + zeta);
	dN[1][1] = 0.125 * m_zeta * p_xi * (1 - xi + 2.0 * eta + zeta);
	dN[1][2] = 0.125 * m_zeta * p_xi * (-1 + xi + 2.0 * eta - zeta);
	dN[1][3] = 0.125 * m_zeta * m_xi * (-1 - xi + 2.0 * eta - zeta);
	dN[1][4] = 0.125 * p_zeta * m_xi * (1 + xi + 2.0 * eta - zeta);
	dN[1][5] = 0.125 * p_zeta * p_xi * (1 - xi + 2.0 * eta - zeta);
	dN[1][6] = 0.125 * p_zeta * p_xi * (-1 + xi + 2.0 * eta + zeta);
	dN[1][7] = 0.125 * p_zeta * m_xi * (-1 - xi + 2.0 * eta + zeta);
	dN[1][8] = - 0.25 * xi2 * m_zeta;
	dN[1][9] = - 0.5 * eta * p_xi * m_zeta;
	dN[1][10] = 0.25 * xi2 * m_zeta;
	dN[1][11] = - 0.5 * eta * m_xi * m_zeta;
	dN[1][12] = - 0.25 * xi2 * p_zeta;
	dN[1][13] = - 0.5 * eta * p_xi * p_zeta;
	dN[1][14] = 0.25 * xi2 * p_zeta;  
	dN[1][15] = - 0.5 * eta * m_xi * p_zeta; 
	dN[1][16] = - 0.25 * m_xi * zeta2;
	dN[1][17] = - 0.25 * p_xi * zeta2;
	dN[1][18] = 0.25 * p_xi * zeta2;
	dN[1][19] = 0.25 * m_xi * zeta2;

	dN[2][0] = 0.125 * m_xi * m_eta * (1 + xi + eta + 2.0 * zeta);
	dN[2][1] = 0.125 * p_xi * m_eta * (1 - xi + eta + 2.0 * zeta);
	dN[2][2] = 0.125 * p_xi * p_eta * (1 - xi - eta + 2.0 * zeta);
	dN[2][3] = 0.125 * m_xi * p_eta * (1 + xi - eta + 2.0 * zeta);
	dN[2][4] = 0.125 * m_xi * m_eta * (-1 - xi - eta + 2.0 * zeta);
	dN[2][5] = 0.125 * p_xi * m_eta * (-1 + xi - eta + 2.0 * zeta);
	dN[2][6] = 0.125 * p_xi * p_eta * (-1 + xi + eta + 2.0 * zeta);
	dN[2][7] = 0.125 * m_xi * p_eta * (-1 - xi + eta + 2.0 * zeta);
	dN[2][8] = - 0.25 * xi2 * m_eta;
	dN[2][9] = - 0.25 * p_xi * eta2;
	dN[2][10] = - 0.25 * xi2 * p_eta;
	dN[2][11] = - 0.25 * m_xi * eta2;
	dN[2][12] = 0.25 * xi2 * m_eta;   
	dN[2][13] = 0.25 * p_xi * eta2;  
	dN[2][14] = 0.25 * xi2 * p_eta;  
	dN[2][15] = 0.25 * m_xi * eta2;  
	dN[2][16] = - 0.5 * zeta * m_xi * m_eta;
	dN[2][17] = - 0.5 * zeta * p_xi * m_eta;
	dN[2][18] = - 0.5 * zeta * p_xi * p_eta;
	dN[2][19] = - 0.5 * zeta * m_xi * p_eta;

	return(NORMAL_RC);
}


static RC dN_hexa1(const double *xi_eta_zeta, double **dN)
{
	double xi = xi_eta_zeta[0];
	double eta = xi_eta_zeta[1];
	double zeta = xi_eta_zeta[2];
	double m_xi = 1.0 - xi;
	double p_xi = 1.0 + xi;
	double m_eta = 1.0 - eta;
	double p_eta = 1.0 + eta;
	double m_zeta = 1.0 - zeta;
	double p_zeta = 1.0 + zeta;

	dN[0][0] = - 0.125 * m_eta * m_zeta;
	dN[0][1] =   0.125 * m_eta * m_zeta;
	dN[0][2] =   0.125 * p_eta * m_zeta;
	dN[0][3] = - 0.125 * p_eta * m_zeta;
	dN[0][4] = - 0.125 * m_eta * p_zeta;
	dN[0][5] =   0.125 * m_eta * p_zeta;
	dN[0][6] =   0.125 * p_eta * p_zeta;
	dN[0][7] = - 0.125 * p_eta * p_zeta;

	dN[1][0] = - 0.125 * m_zeta * m_xi;
	dN[1][1] = - 0.125 * m_zeta * p_xi;
	dN[1][2] =   0.125 * m_zeta * p_xi;
	dN[1][3] =   0.125 * m_zeta * m_xi;
	dN[1][4] = - 0.125 * p_zeta * m_xi;
	dN[1][5] = - 0.125 * p_zeta * p_xi;
	dN[1][6] =   0.125 * p_zeta * p_xi;
	dN[1][7] =   0.125 * p_zeta * m_xi;

	dN[2][0] = - 0.125 * m_xi * m_eta;
	dN[2][1] = - 0.125 * p_xi * m_eta;
	dN[2][2] = - 0.125 * p_xi * p_eta;
	dN[2][3] = - 0.125 * m_xi * p_eta;
	dN[2][4] =   0.125 * m_xi * m_eta;
	dN[2][5] =   0.125 * p_xi * m_eta;
	dN[2][6] =   0.125 * p_xi * p_eta;
	dN[2][7] =   0.125 * m_xi * p_eta;

	return(NORMAL_RC);
}


static RC dN_quad2(const double *xi_eta, double **dN)
{
	double xi = xi_eta[0];
	double eta = xi_eta[1];
	double m_xi = 1.0 - xi;
	double p_xi = 1.0 + xi;
	double m_eta = 1.0 - eta;
	double p_eta = 1.0 + eta;
	double eta2 = 1.0 - eta * eta;
	double xi2 = 1.0 - xi * xi;

	dN[0][0] = 0.25 * m_eta * (2.0 * xi + eta);
	dN[0][1] = 0.25 * m_eta * (2.0 * xi - eta);
	dN[0][2] = 0.25 * p_eta * (2.0 * xi + eta);
	dN[0][3] = 0.25 * p_eta * (2.0 * xi - eta);
	dN[0][4] = - xi * m_eta;
	dN[0][5] = 0.5 * eta2;
	dN[0][6] = - xi * p_eta;
	dN[0][7] = - 0.5 * eta2;
	
	dN[1][0] = 0.25 * m_xi * (2.0 * eta + xi);
	dN[1][1] = 0.25 * p_xi * (2.0 * eta - xi);
	dN[1][2] = 0.25 * p_xi * (2.0 * eta + xi);
	dN[1][3] = 0.25 * m_xi * (2.0 * eta - xi);
	dN[1][4] = - 0.5 * xi2;
	dN[1][5] = - eta * p_xi;
	dN[1][6] = 0.5 * xi2;
	dN[1][7] = - eta * m_xi;

	return(NORMAL_RC);
}


static RC dN_quad1(const double *xi_eta, double **dN)
{
	double xi = xi_eta[0];
	double eta = xi_eta[1];
	double m_xi = 1.0 - xi;
	double p_xi = 1.0 + xi;
	double m_eta = 1.0 - eta;
	double p_eta = 1.0 + eta;

	dN[0][0] = - 0.25 * m_eta;
	dN[0][1] = 0.25 * m_eta;
	dN[0][2] = 0.25 * p_eta;
	dN[0][3] = -0.25 * p_eta;
	
	dN[1][0] = - 0.25 * m_xi;
	dN[1][1] = - 0.25 * p_xi;
	dN[1][2] = 0.25 * p_xi;
	dN[1][3] = 0.25 * m_xi;

	return(NORMAL_RC);
}


static RC dN_penta2(const double *l1_l2_l3_zeta, double **dN)
{
	int ii1, ii2;
	double l1 = l1_l2_l3_zeta[0];
	double l2 = l1_l2_l3_zeta[1];
	double l3 = l1_l2_l3_zeta[2];
	double zeta = l1_l2_l3_zeta[3];
	double m_zeta = 1.0 - zeta;
	double p_zeta = 1.0 + zeta;
	double zeta2 = 1.0 - zeta * zeta;
	double const_matrix[3][4] = {{ 1.0,  0.0, -1.0,  0.0},
	                             { 0.0,  1.0, -1.0,  0.0},
	                             { 0.0,  0.0,  0.0,  1.0}};
	double tmp_matrix[4][15];

	tmp_matrix[0][0] = l1*m_zeta + 0.5*m_zeta*(-2.0 + 2.0*l1 - zeta);
	tmp_matrix[0][1] = 0.0;
	tmp_matrix[0][2] = 0.0;
	tmp_matrix[0][3] = l1*p_zeta + 0.5*p_zeta*(-2.0 + 2.0*l1 + zeta);
	tmp_matrix[0][4] = 0.0;
	tmp_matrix[0][5] = 0.0;
	tmp_matrix[0][6] = 0.0;
	tmp_matrix[0][7] = 2.0 * l3 * m_zeta;
	tmp_matrix[0][8] = 2.0 * l2 * m_zeta;
	tmp_matrix[0][9] = 0.0;
	tmp_matrix[0][10] = 2.0 * l3 * p_zeta;
	tmp_matrix[0][11] = 2.0 * l2 * p_zeta;
	tmp_matrix[0][12] = zeta2;
	tmp_matrix[0][13] = 0.0;
	tmp_matrix[0][14] = 0.0;

	tmp_matrix[1][0] = 0.0;
	tmp_matrix[1][1] = l2*m_zeta + 0.5*m_zeta*(-2.0 + 2.0*l2 - zeta);
	tmp_matrix[1][2] = 0.0;
	tmp_matrix[1][3] = 0.0;
	tmp_matrix[1][4] = l2*p_zeta + 0.5*p_zeta*(-2.0 + 2.0*l2 + zeta);
	tmp_matrix[1][5] = 0.0;
	tmp_matrix[1][6] = 2.0 * l3 * m_zeta;
	tmp_matrix[1][7] = 0.0;
	tmp_matrix[1][8] = 2.0 * l1 * m_zeta;
	tmp_matrix[1][9] = 2.0 * l3 * p_zeta;
	tmp_matrix[1][10] = 0.0;
	tmp_matrix[1][11] = 2.0 * l1 * p_zeta;
	tmp_matrix[1][12] = 0.0;
	tmp_matrix[1][13] = zeta2;
	tmp_matrix[1][14] = 0.0;

	tmp_matrix[2][0] = 0.0;
	tmp_matrix[2][1] = 0.0;
	tmp_matrix[2][2] = l3*m_zeta + 0.5*m_zeta*(-2.0 + 2.0*l3 - zeta);
	tmp_matrix[2][3] = 0.0;
	tmp_matrix[2][4] = 0.0;
	tmp_matrix[2][5] = l3*p_zeta + 0.5*p_zeta*(-2.0 + 2.0*l3 + zeta);
	tmp_matrix[2][6] = 2.0 * l2 * m_zeta;
	tmp_matrix[2][7] = 2.0 * l1 * m_zeta;
	tmp_matrix[2][8] = 0.0;
	tmp_matrix[2][9] = 2.0 * l2 * p_zeta;
	tmp_matrix[2][10] = 2.0 * l1 * p_zeta;
	tmp_matrix[2][11] = 0.0;
	tmp_matrix[2][12] = 0.0;
	tmp_matrix[2][13] = 0.0;
	tmp_matrix[2][14] = zeta2;

	tmp_matrix[3][0] = - 0.5*l1*m_zeta - 0.5*l1*(-2.0 + 2.0*l1 - zeta);
	tmp_matrix[3][1] = - 0.5*l2*m_zeta - 0.5*l2*(-2.0 + 2.0*l2 - zeta);
	tmp_matrix[3][2] = - 0.5*l3*m_zeta - 0.5*l3*(-2.0 + 2.0*l3 - zeta);
	tmp_matrix[3][3] = 0.5*l1*p_zeta + 0.5*l1*(-2.0 + 2.0*l1 + zeta);
	tmp_matrix[3][4] = 0.5*l2*p_zeta + 0.5*l2*(-2.0 + 2.0*l2 + zeta);
	tmp_matrix[3][5] = 0.5*l3*p_zeta + 0.5*l3*(-2.0 + 2.0*l3 + zeta);
	tmp_matrix[3][6] = - 2.0 * l2 * l3;
	tmp_matrix[3][7] = - 2.0 * l3 * l1;
	tmp_matrix[3][8] = - 2.0 * l1 * l2;
	tmp_matrix[3][9] = 2.0 * l2 * l3;
	tmp_matrix[3][10] = 2.0 * l3 * l1;
	tmp_matrix[3][11] = 2.0 * l1 * l2;
	tmp_matrix[3][12] = - 2.0 * l1 * zeta;
	tmp_matrix[3][13] = - 2.0 * l2 * zeta;
	tmp_matrix[3][14] = - 2.0 * l3 * zeta;

	/* const_matrix * tmp_matrix => dN */
	for(ii1=0; ii1<3; ii1++){
		for(ii2=0; ii2<15; ii2++){
			dN[ii1][ii2] = const_matrix[ii1][0] * tmp_matrix[0][ii2]
                         + const_matrix[ii1][1] * tmp_matrix[1][ii2]
                         + const_matrix[ii1][2] * tmp_matrix[2][ii2]
                         + const_matrix[ii1][3] * tmp_matrix[3][ii2];
		}
	}

	return(NORMAL_RC);
}


static RC dN_penta1(const double *l1_l2_l3_zeta, double **dN)
{
	int ii1, ii2;
	double l1 = l1_l2_l3_zeta[0];
	double l2 = l1_l2_l3_zeta[1];
	double l3 = l1_l2_l3_zeta[2];
	double zeta = l1_l2_l3_zeta[3];
	double m_zeta = 1.0 - zeta;
	double p_zeta = 1.0 + zeta;
	double const_matrix[3][4] = {{ 1.0,  0.0, -1.0,  0.0},
	                             { 0.0,  1.0, -1.0,  0.0},
	                             { 0.0,  0.0,  0.0,  1.0}};
	double tmp_matrix[4][6];

	tmp_matrix[0][0] = 0.5 * m_zeta;
	tmp_matrix[0][1] = 0.0;
	tmp_matrix[0][2] = 0.0;
	tmp_matrix[0][3] = 0.5 * p_zeta;
	tmp_matrix[0][4] = 0.0;
	tmp_matrix[0][5] = 0.0;

	tmp_matrix[1][0] = 0.0;
	tmp_matrix[1][1] = 0.5 * m_zeta;
	tmp_matrix[1][2] = 0.0;
	tmp_matrix[1][3] = 0.0;
	tmp_matrix[1][4] = 0.5 * p_zeta;
	tmp_matrix[1][5] = 0.0;

	tmp_matrix[2][0] = 0.0;
	tmp_matrix[2][1] = 0.0;
	tmp_matrix[2][2] = 0.5 * m_zeta;
	tmp_matrix[2][3] = 0.0;
	tmp_matrix[2][4] = 0.0;
	tmp_matrix[2][5] = 0.5 * p_zeta;

	tmp_matrix[3][0] = - 0.5 * l1;
	tmp_matrix[3][1] = - 0.5 * l2;
	tmp_matrix[3][2] = - 0.5 * l3;
	tmp_matrix[3][3] = 0.5 * l1;
	tmp_matrix[3][4] = 0.5 * l2;
	tmp_matrix[3][5] = 0.5 * l3;

	/* const_matrix * tmp_matrix => dN */
	for(ii1=0; ii1<3; ii1++){
		for(ii2=0; ii2<6; ii2++){
			dN[ii1][ii2] = const_matrix[ii1][0] * tmp_matrix[0][ii2]
                         + const_matrix[ii1][1] * tmp_matrix[1][ii2]
                         + const_matrix[ii1][2] * tmp_matrix[2][ii2]
                         + const_matrix[ii1][3] * tmp_matrix[3][ii2];
		}
	}

	return(NORMAL_RC);
}


static RC dN_tetra2(const double *l1_l2_l3_l4, double **dN)
{
	int ii1, ii2;
	double l1 = l1_l2_l3_l4[0];
	double l2 = l1_l2_l3_l4[1];
	double l3 = l1_l2_l3_l4[2];
	double l4 = l1_l2_l3_l4[3];
	double const_matrix[3][4] = {{ 1.0,  0.0,  0.0, -1.0},
	                             { 0.0,  1.0,  0.0, -1.0},
	                             { 0.0,  0.0,  1.0, -1.0}};
	double tmp_matrix[4][10];

	tmp_matrix[0][0] = 4.0 * l1 - 1.0;
	tmp_matrix[0][1] = 0.0;
	tmp_matrix[0][2] = 0.0;
	tmp_matrix[0][3] = 0.0;
	tmp_matrix[0][4] = 0.0;
	tmp_matrix[0][5] = 4.0 * l3;
	tmp_matrix[0][6] = 4.0 * l2;
	tmp_matrix[0][7] = 4.0 * l4;
	tmp_matrix[0][8] = 0.0;
	tmp_matrix[0][9] = 0.0;

	tmp_matrix[1][0] = 0.0;
	tmp_matrix[1][1] = 4.0 * l2 - 1.0;
	tmp_matrix[1][2] = 0.0;
	tmp_matrix[1][3] = 0.0;
	tmp_matrix[1][4] = 4.0 * l3;
	tmp_matrix[1][5] = 0.0;
	tmp_matrix[1][6] = 4.0 * l1;
	tmp_matrix[1][7] = 0.0;
	tmp_matrix[1][8] = 4.0 * l4;
	tmp_matrix[1][9] = 0.0;

	tmp_matrix[2][0] = 0.0;
	tmp_matrix[2][1] = 0.0;
	tmp_matrix[2][2] = 4.0 * l3 - 1.0;
	tmp_matrix[2][3] = 0.0;
	tmp_matrix[2][4] = 4.0 * l2;
	tmp_matrix[2][5] = 4.0 * l1;
	tmp_matrix[2][6] = 0.0;
	tmp_matrix[2][7] = 0.0;
	tmp_matrix[2][8] = 0.0;
	tmp_matrix[2][9] = 4.0 * l4;

	tmp_matrix[3][0] = 0.0;
	tmp_matrix[3][1] = 0.0;
	tmp_matrix[3][2] = 0.0;
	tmp_matrix[3][3] = 4.0 * l4 - 1.0;
	tmp_matrix[3][4] = 0.0;
	tmp_matrix[3][5] = 0.0;
	tmp_matrix[3][6] = 0.0;
	tmp_matrix[3][7] = 4.0 * l1;
	tmp_matrix[3][8] = 4.0 * l2;
	tmp_matrix[3][9] = 4.0 * l3;

	/* const_matrix * tmp_matrix => dN */
	for(ii1=0; ii1<3; ii1++){
		for(ii2=0; ii2<10; ii2++){
			dN[ii1][ii2] = const_matrix[ii1][0] * tmp_matrix[0][ii2]
                         + const_matrix[ii1][1] * tmp_matrix[1][ii2]
                         + const_matrix[ii1][2] * tmp_matrix[2][ii2]
                         + const_matrix[ii1][3] * tmp_matrix[3][ii2];
		}
	}

	return(NORMAL_RC);
}


static RC dN_tetra1(const double *l1_l2_l3_l4, double **dN)
{
	int ii1, ii2;
	double const_matrix[3][4] = {{ 1.0,  0.0,  0.0, -1.0},
	                             { 0.0,  1.0,  0.0, -1.0},
	                             { 0.0,  0.0,  1.0, -1.0}};
	double tmp_matrix[4][4]   = {{ 1.0,  0.0,  0.0,  0.0},
	                             { 0.0,  1.0,  0.0,  0.0},
	                             { 0.0,  0.0,  1.0,  0.0},
	                             { 0.0,  0.0,  0.0,  1.0}};

	/* const_matrix * tmp_matrix => dN */
	for(ii1=0; ii1<3; ii1++){
		for(ii2=0; ii2<4; ii2++){
			dN[ii1][ii2] = const_matrix[ii1][0] * tmp_matrix[0][ii2]
                         + const_matrix[ii1][1] * tmp_matrix[1][ii2]
                         + const_matrix[ii1][2] * tmp_matrix[2][ii2]
                         + const_matrix[ii1][3] * tmp_matrix[3][ii2];
		}
	}

	return(NORMAL_RC);
}


static RC dN_tri2(const double *l1_l2_l3, double **dN)
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
	tmp_matrix[0][3] = 0.0;
	tmp_matrix[0][4] = 4.0 * l3;
	tmp_matrix[0][5] = 4.0 * l2;

	tmp_matrix[1][0] = 0.0;
	tmp_matrix[1][1] = 4.0 * l2 - 1.0;
	tmp_matrix[1][2] = 0.0;
	tmp_matrix[1][3] = 4.0 * l3;
	tmp_matrix[1][4] = 0.0;
	tmp_matrix[1][5] = 4.0 * l1;

	tmp_matrix[2][0] = 0.0;
	tmp_matrix[2][1] = 0.0;
	tmp_matrix[2][2] = 4.0 * l3 - 1.0;
	tmp_matrix[2][3] = 4.0 * l2;
	tmp_matrix[2][4] = 4.0 * l1;
	tmp_matrix[2][5] = 0.0;

	/* const_matrix * tmp_matrix => dN */
	for(ii1=0; ii1<2; ii1++){
		for(ii2=0; ii2<6; ii2++){
			dN[ii1][ii2] = const_matrix[ii1][0] * tmp_matrix[0][ii2]
                         + const_matrix[ii1][1] * tmp_matrix[1][ii2]
                         + const_matrix[ii1][2] * tmp_matrix[2][ii2];
		}
	}

	return(NORMAL_RC);
}


static RC dN_tri1(const double *l1_l2_l3, double **dN)
{
	int ii1, ii2;
	double const_matrix[2][3] = {{ 1.0,  0.0, -1.0},
	                             { 0.0,  1.0, -1.0}};
	double tmp_matrix[3][3] =   {{ 1.0,  0.0,  0.0},
	                             { 0.0,  1.0,  0.0},
	                             { 0.0,  0.0,  1.0}};

	/* const_matrix * tmp_matrix => dN */
	for(ii1=0; ii1<2; ii1++){
		for(ii2=0; ii2<3; ii2++){
			dN[ii1][ii2] = const_matrix[ii1][0] * tmp_matrix[0][ii2]
                         + const_matrix[ii1][1] * tmp_matrix[1][ii2]
                         + const_matrix[ii1][2] * tmp_matrix[2][ii2];
		}
	}

	return(NORMAL_RC);
}


static RC dN_line1(const double *xi, double **dN)
{
	dN[0][0] = -0.5;
	dN[0][1] = 0.5;

	return(NORMAL_RC);
}


static RC dN_line2(const double *xi, double **dN)
{
	dN[0][0] = -0.5 + xi[0];
	dN[0][1] = 0.5 + xi[0];
	dN[0][2] = -2.0 * xi[0];

	return(NORMAL_RC);
}


static RC N_line1(const double *xi, double *N)
{
	N[0] = 0.5 * (1.0 - xi[0]);
	N[1] = 0.5 * (1.0 + xi[0]);

	return(NORMAL_RC);
}


static RC N_line2(const double *xi, double *N)
{
	N[0] = -0.5 * xi[0] * (1.0 - xi[0]);
	N[1] = 0.5 * xi[0] * (1.0 + xi[0]);
	N[2] = 1.0 - (xi[0] * xi[0]);

	return(NORMAL_RC);
}


static RC N_tri1(const double *l1_l2_l3, double *N)
{
	double l1 = l1_l2_l3[0];
	double l2 = l1_l2_l3[1];
	double l3 = l1_l2_l3[2];

	N[0] = l1;
	N[1] = l2;
	N[2] = l3;

	return(NORMAL_RC);
}


static RC N_tri2(const double *l1_l2_l3, double *N)
{
	double l1 = l1_l2_l3[0];
	double l2 = l1_l2_l3[1];
	double l3 = l1_l2_l3[2];

	N[0] = l1 * (2.0 * l1 - 1.0);
	N[1] = l2 * (2.0 * l2 - 1.0);
	N[2] = l3 * (2.0 * l3 - 1.0);
	N[3] = 4.0 * l2 * l3;
	N[4] = 4.0 * l3 * l1;
	N[5] = 4.0 * l1 * l2;

	return(NORMAL_RC);
}


static RC N_quad1(const double *xi_eta, double *N)
{
	double xi = xi_eta[0];
	double eta = xi_eta[1];
	double m_xi = 1.0 - xi;
	double p_xi = 1.0 + xi;
	double m_eta = 1.0 - eta;
	double p_eta = 1.0 + eta;

	N[0] = 0.25 * m_xi * m_eta;
	N[1] = 0.25 * p_xi * m_eta;
	N[2] = 0.25 * p_xi * p_eta;
	N[3] = 0.25 * m_xi * p_eta;

	return(NORMAL_RC);
}


static RC N_quad2(const double *xi_eta, double *N)
{
	double xi = xi_eta[0];
	double eta = xi_eta[1];
	double m_xi = 1.0 - xi;
	double p_xi = 1.0 + xi;
	double m_eta = 1.0 - eta;
	double p_eta = 1.0 + eta;

	N[0] = 0.25 * m_xi * m_eta * (-1.0 - xi - eta);
	N[1] = 0.25 * p_xi * m_eta * (-1.0 + xi - eta);
	N[2] = 0.25 * p_xi * p_eta * (-1.0 + xi + eta);
	N[3] = 0.25 * m_xi * p_eta * (-1.0 - xi + eta);
	N[4] = 0.5 * m_eta * (1 - xi * xi);
	N[5] = 0.5 * p_xi * (1 - eta * eta);
	N[6] = 0.5 * p_eta * (1 - xi * xi);
	N[7] = 0.5 * m_xi * (1 - eta * eta);

	return(NORMAL_RC);
}


static RC N_tetra1(const double *l1_l2_l3_l4, double *N)
{
	double l1 = l1_l2_l3_l4[0];
	double l2 = l1_l2_l3_l4[1];
	double l3 = l1_l2_l3_l4[2];
	double l4 = l1_l2_l3_l4[3];

	N[0] = l1;
	N[1] = l2;
	N[2] = l3;
	N[3] = l4;

	return(NORMAL_RC);
}


static RC N_tetra2(const double *l1_l2_l3_l4, double *N)
{
	double l1 = l1_l2_l3_l4[0];
	double l2 = l1_l2_l3_l4[1];
	double l3 = l1_l2_l3_l4[2];
	double l4 = l1_l2_l3_l4[3];

	N[0] = l1 * (2.0 * l1 - 1.0);
	N[1] = l2 * (2.0 * l2 - 1.0);
	N[2] = l3 * (2.0 * l3 - 1.0);
	N[3] = l4 * (2.0 * l4 - 1.0);
	N[4] = 4.0 * l2 * l3;
	N[5] = 4.0 * l3 * l1;
	N[6] = 4.0 * l1 * l2;
	N[7] = 4.0 * l1 * l4;
	N[8] = 4.0 * l2 * l4;
	N[9] = 4.0 * l3 * l4;

	return(NORMAL_RC);
}


static RC N_penta1(const double *l1_l2_l3_zeta, double *N)
{
	double l1 = l1_l2_l3_zeta[0];
	double l2 = l1_l2_l3_zeta[1];
	double l3 = l1_l2_l3_zeta[2];
	double zeta = l1_l2_l3_zeta[3];
	double m_zeta = 1.0 - zeta;
	double p_zeta = 1.0 + zeta;

	N[0] = 0.5 * l1 * m_zeta;
	N[1] = 0.5 * l2 * m_zeta;
	N[2] = 0.5 * l3 * m_zeta;
	N[3] = 0.5 * l1 * p_zeta;
	N[4] = 0.5 * l2 * p_zeta;
	N[5] = 0.5 * l3 * p_zeta;

	return(NORMAL_RC);
}


static RC N_penta2(const double *l1_l2_l3_zeta, double *N)
{
	double l1 = l1_l2_l3_zeta[0];
	double l2 = l1_l2_l3_zeta[1];
	double l3 = l1_l2_l3_zeta[2];
	double zeta = l1_l2_l3_zeta[3];
	double m_zeta = 1.0 - zeta;
	double p_zeta = 1.0 + zeta;

	N[0] = 0.5 * l1 * (2.0 * l1 - 2.0 - zeta) * m_zeta;
	N[1] = 0.5 * l2 * (2.0 * l2 - 2.0 - zeta) * m_zeta;
	N[2] = 0.5 * l3 * (2.0 * l3 - 2.0 - zeta) * m_zeta;
	N[3] = 0.5 * l1 * (2.0 * l1 - 2.0 + zeta) * p_zeta;
	N[4] = 0.5 * l2 * (2.0 * l2 - 2.0 + zeta) * p_zeta;
	N[5] = 0.5 * l3 * (2.0 * l3 - 2.0 + zeta) * p_zeta;
	N[6] = 2.0 * l2 * l3 * m_zeta;
	N[7] = 2.0 * l3 * l1 * m_zeta;
	N[8] = 2.0 * l1 * l2 * m_zeta;
	N[9] = 2.0 * l2 * l3 * p_zeta;
	N[10] = 2.0 * l3 * l1 * p_zeta;
	N[11] = 2.0 * l1 * l2 * p_zeta;
	N[12] = l1 * (1.0 - zeta * zeta);
	N[13] = l2 * (1.0 - zeta * zeta);
	N[14] = l3 * (1.0 - zeta * zeta);

	return(NORMAL_RC);
}


static RC N_hexa1(const double *xi_eta_zeta, double *N)
{
	double xi = xi_eta_zeta[0];
	double eta = xi_eta_zeta[1];
	double zeta = xi_eta_zeta[2];
	double m_xi = 1.0 - xi;
	double p_xi = 1.0 + xi;
	double m_eta = 1.0 - eta;
	double p_eta = 1.0 + eta;
	double m_zeta = 1.0 - zeta;
	double p_zeta = 1.0 + zeta;

	N[0] = 0.125 * m_xi * m_eta * m_zeta;
	N[1] = 0.125 * p_xi * m_eta * m_zeta;
	N[2] = 0.125 * p_xi * p_eta * m_zeta;
	N[3] = 0.125 * m_xi * p_eta * m_zeta;
	N[4] = 0.125 * m_xi * m_eta * p_zeta;
	N[5] = 0.125 * p_xi * m_eta * p_zeta;
	N[6] = 0.125 * p_xi * p_eta * p_zeta;
	N[7] = 0.125 * m_xi * p_eta * p_zeta;

	return(NORMAL_RC);
}


static RC N_hexa2(const double *xi_eta_zeta, double *N)
{
	double xi = xi_eta_zeta[0];
	double eta = xi_eta_zeta[1];
	double zeta = xi_eta_zeta[2];
	double m_xi = 1.0 - xi;
	double p_xi = 1.0 + xi;
	double m_eta = 1.0 - eta;
	double p_eta = 1.0 + eta;
	double m_zeta = 1.0 - zeta;
	double p_zeta = 1.0 + zeta;
	double xi2 = 1.0 - xi * xi;
	double eta2 = 1.0 - eta * eta;
	double zeta2 = 1.0 - zeta * zeta;

	N[0] = -0.125 * m_xi * m_eta * m_zeta * (2.0 + xi + eta + zeta);
	N[1] = -0.125 * p_xi * m_eta * m_zeta * (2.0 - xi + eta + zeta);
	N[2] = -0.125 * p_xi * p_eta * m_zeta * (2.0 - xi - eta + zeta);
	N[3] = -0.125 * m_xi * p_eta * m_zeta * (2.0 + xi - eta + zeta);
	N[4] = -0.125 * m_xi * m_eta * p_zeta * (2.0 + xi + eta - zeta);
	N[5] = -0.125 * p_xi * m_eta * p_zeta * (2.0 - xi + eta - zeta);
	N[6] = -0.125 * p_xi * p_eta * p_zeta * (2.0 - xi - eta - zeta);
	N[7] = -0.125 * m_xi * p_eta * p_zeta * (2.0 + xi - eta - zeta);
	N[8] = 0.25 * xi2 * m_eta * m_zeta;
	N[9] = 0.25 * p_xi * eta2 * m_zeta;
	N[10] = 0.25 * xi2 * p_eta * m_zeta;
	N[11] = 0.25 * m_xi * eta2 * m_zeta;
	N[12] = 0.25 * xi2 * m_eta * p_zeta;
	N[13] = 0.25 * p_xi * eta2 * p_zeta;
	N[14] = 0.25 * xi2 * p_eta * p_zeta;
	N[15] = 0.25 * m_xi * eta2 * p_zeta;
	N[16] = 0.25 * m_xi * m_eta * zeta2;
	N[17] = 0.25 * p_xi * m_eta * zeta2;
	N[18] = 0.25 * p_xi * p_eta * zeta2;
	N[19] = 0.25 * m_xi * p_eta * zeta2;

	return(NORMAL_RC);
}


