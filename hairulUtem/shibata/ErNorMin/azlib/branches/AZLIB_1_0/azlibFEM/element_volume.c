/*********************************************************************
 * element_volume.c
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

/* $Id: element_volume.c,v 1.6 2003/07/22 09:16:32 nagatani Exp $ */

#include <stdio.h>
#include <math.h>
#include "rc.h"
#include "math_utl.h"
#include "fem_struct.h"


static RC Gauss_const4volume(FEM_ELEM_TYPE type, double **points,
                                 double *weight, int *num_point);


RC total_element_volume(FEM_ELEMENT_ARRAY element,
                        FEM_NODE_ARRAY node, double *total_volume)
{
	int ii1;
	long double sum = 0.0;

	for(ii1=0; ii1<(element.size); ii1++){
		RC_TRY( element_volume(&(element.array[ii1]), node) );
		sum += element.array[ii1].volume;
	}
	*total_volume = (double)sum;

	return(NORMAL_RC);
}


RC element_volume(FEM_ELEMENT *element, FEM_NODE_ARRAY node)
{
	static int allocate_flag = 0;
	static FEM_ELEM_TYPE latest_type;
	static int num_point;
	static double **node_location = NULL;
	static double **gauss_point = NULL;
	static double weight[MAX_GAUSS_POINT];
	static double **dN = NULL;
	static double **jacobi = NULL;
	double sum;
	int ii1;
	double f;
	TRANSLATION3D d_xi, d_eta, d_A;

	if(allocate_flag == 0){
		RC_TRY( allocate2D(FEM_MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(MAX_GAUSS_POINT, 4, &gauss_point) );
		RC_TRY( allocate2D(3, FEM_MAX_NODE, &dN) );
		RC_TRY( allocate2D(3, 3, &jacobi) );

		RC_TRY( Gauss_const4volume(element->type,
		                           gauss_point, weight, &num_point) );
		latest_type = element->type;
		allocate_flag = 1;
	}

	if(element->type != latest_type){
		RC_TRY( Gauss_const4volume(element->type,
		                           gauss_point, weight, &num_point) );
		latest_type = element->type;
	}

	RC_TRY( node_location_matrix(element, node, node_location) );

	sum = 0.0;
	for(ii1=0;ii1<num_point;ii1++){
		init_matrix(3, FEM_MAX_NODE, dN);
		RC_TRY( set_dN_matrix(element->type, gauss_point[ii1], dN) );

		mul_matrix(3, element->node_num, 3, dN, node_location, jacobi);
		switch(element->type){
		case ELEM_TETRA1:
		case ELEM_TETRA2:
		case ELEM_PENTA1:
		case ELEM_PENTA2:
		case ELEM_HEXA1:
		case ELEM_HEXA2:
			f = determinant3(jacobi);
			break;
		case ELEM_TRI1:
		case ELEM_TRI2:
		case ELEM_QUAD1:
		case ELEM_QUAD2:
			d_xi.x = jacobi[0][0];
			d_xi.y = jacobi[0][1];
			d_xi.z = jacobi[0][2];
			d_eta.x = jacobi[1][0];
			d_eta.y = jacobi[1][1];
			d_eta.z = jacobi[1][2];
			d_A = outer_product3d(d_xi, d_eta);
			f = sqrt(d_A.x*d_A.x + d_A.y*d_A.y + d_A.z*d_A.z);
			break;
		case ELEM_RBAR:
		case ELEM_BEAM1:
		case ELEM_LINE1:
		case ELEM_BEAM2:
		case ELEM_LINE2:
			d_xi.x = jacobi[0][0];
			d_xi.y = jacobi[0][1];
			d_xi.z = jacobi[0][2];
			f = sqrt(d_xi.x*d_xi.x + d_xi.y*d_xi.y + d_xi.z*d_xi.z);
			break;
		default:
			return(IMPLEMENT_ERROR_RC);
		}

		if(f < 0.0){
			fprintf(stderr, "warning: negative volume[%d-%d]\n",
			        element->label, ii1);
		}

		sum += 1.0 * f * weight[ii1];
	}
	element->volume = sum;

	return(NORMAL_RC);
}


static RC Gauss_const4volume(FEM_ELEM_TYPE type, double **points,
                                 double *weight, int *num_point)
{
	switch(type){
	case ELEM_RBAR:
	case ELEM_BEAM1:
	case ELEM_BEAM2:
	case ELEM_LINE1:
	case ELEM_LINE2:
		*num_point = 2;
		RC_TRY( Gauss_const_line(*num_point, points, weight) );
		break;
	case ELEM_TRI1:
	case ELEM_TRI2:
		*num_point = 3;
		RC_TRY( Gauss_const_tri(*num_point, points, weight) );
		break;
	case ELEM_QUAD1:
	case ELEM_QUAD2:
		*num_point = 4;
		RC_TRY( Gauss_const_quad(*num_point, points, weight) );
		break;
	case ELEM_TETRA1:
	case ELEM_TETRA2:
		*num_point = 4;
		RC_TRY( Gauss_const_tetra(*num_point, points, weight) );
		break;
	case ELEM_PENTA1:
	case ELEM_PENTA2:
		*num_point = 6;
		RC_TRY( Gauss_const_penta(*num_point, points, weight) );
		break;
	case ELEM_HEXA1:
	case ELEM_HEXA2:
		*num_point = 8;
		RC_TRY( Gauss_const_hexa(*num_point, points, weight) );
		break;
	default : 
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}

