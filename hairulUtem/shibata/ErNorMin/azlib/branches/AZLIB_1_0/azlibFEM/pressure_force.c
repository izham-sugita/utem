/*********************************************************************
 * pressure_force.c
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

/* $Id: pressure_force.c,v 1.9 2003/07/23 04:06:32 sasaoka Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rc.h"
#include "math_utl.h"
#include "fem_struct.h"


static RC Gauss_const4pressure(FEM_ELEM_TYPE type, double **points,
                                    double *weight, int *num_point);


static RC Gauss_const4pressure(FEM_ELEM_TYPE type, double **points,
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
	default : 
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}



/* $B#2<!85$^$?$O#1<!85MWAG(B element $B$K05NO(B p[p_size] $B$,:nMQ$9$k>l9g$N(B */
/* $B@aE@2Y=E%Y%/%H%k$r5a$a$k!#$?$@$7!"0zD%$jB&$r@5$H$9$k(B */
/* p_size $B$,#1$N>l9g$O05NO$,MWAGFb$G0lDj!"(B*/
/* P_size > 1 $B$N>l9g$O!"3F@aE@$G$NCM$H$7$FFbA^$9$k(B */
RC pressure_force(int p_size, const double *p, FEM_ELEMENT *element,
                  FEM_NODE_ARRAY node, TRANSLATION3D *f)
{
	static int allocate_flag = 0;
	static FEM_ELEM_TYPE latest_type;
	static int num_point;
	static double **node_location = NULL;
	static double **gauss_point = NULL;
	static double weight[MAX_GAUSS_POINT];
	static double **dN = NULL;
	static double **jacobi = NULL;
	double Np[FEM_MAX_NODE];
	double Nf[FEM_MAX_NODE];
	int ii1, ii2;
	TRANSLATION3D d_xi, d_eta, d_A;
	FEM_ELEM_TYPE p_type = ELEM_TRI1;

	if(allocate_flag == 0){
		RC_TRY( allocate2D(FEM_MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(MAX_GAUSS_POINT, 4, &gauss_point) );
		RC_TRY( allocate2D(3, FEM_MAX_NODE, &dN) );
		RC_TRY( allocate2D(3, 3, &jacobi) );

		RC_TRY( Gauss_const4pressure(element->type,
		                           gauss_point, weight, &num_point) );
		latest_type = element->type;
		allocate_flag = 1;
	}

	if(element->type != latest_type){
		RC_TRY( Gauss_const4pressure(element->type,
		                           gauss_point, weight, &num_point) );
		latest_type = element->type;
	}

	RC_TRY( node_location_matrix(element, node, node_location) );

	if(p_size != 1){
		switch(element->type){
		case ELEM_TRI1:
		case ELEM_TRI2:
			if(p_size == 3){
				p_type = ELEM_TRI1;
			}else if(p_size == 6){
				p_type = ELEM_TRI2;
			}else{
				return(IMPLEMENT_ERROR_RC);
			}
			break;
		case ELEM_QUAD1:
		case ELEM_QUAD2:
			if(p_size == 4){
				p_type = ELEM_QUAD1;
			}else if(p_size == 8){
				p_type = ELEM_QUAD2;
			}else{
				return(IMPLEMENT_ERROR_RC);
			}
			break;
		case ELEM_RBAR:
		case ELEM_BEAM1:
		case ELEM_BEAM2:
		case ELEM_LINE1:
		case ELEM_LINE2:
			if(p_size == 2){
				p_type = ELEM_LINE1;
			}else if(p_size == 3){
				p_type = ELEM_LINE2;
			}else{
				return(IMPLEMENT_ERROR_RC);
			}
			break;
		default:
			return(IMPLEMENT_ERROR_RC);
		}
	}

	for(ii1=0; ii1<element->node_num; ii1++){
		f[ii1].x = f[ii1].y = f[ii1].z = 0.0;
	}

	for(ii1=0; ii1<num_point; ii1++){
		double lp;

		init_matrix(3, FEM_MAX_NODE, dN);  /* $B=i4|2=$,I,MW(B */
		RC_TRY( set_dN_matrix(element->type, gauss_point[ii1], dN) );
		RC_TRY( set_N_vector(element->type, gauss_point[ii1], Nf) );
		mul_matrix(3, element->node_num, 3, dN, node_location, jacobi);

		/* p $B$rFbA^(B */
		if(p_size == 1){
			lp = p[0];
		}else{
			RC_TRY( set_N_vector(p_type, gauss_point[ii1], Np) );
			lp = 0.0;
			for(ii2=0; ii2 < p_size; ii2++){
				lp += p[ii2]* Np[ii2];
			}
		}

		switch(element->type){
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
			break;
		case ELEM_RBAR:
		case ELEM_BEAM1:
		case ELEM_LINE1:
		case ELEM_BEAM2:
		case ELEM_LINE2:
/*			d_A.x = -jacobi[0][1];*/
			d_A.x = jacobi[0][1];
/*			d_A.y = jacobi[0][0];*/
			d_A.y = -jacobi[0][0];
			d_A.z = 0.0;
			break;
		default:
			return(IMPLEMENT_ERROR_RC);
		}

		for(ii2=0; ii2<element->node_num; ii2++){
			f[ii2].x += weight[ii1] * lp * Nf[ii2] * d_A.x;
			f[ii2].y += weight[ii1] * lp * Nf[ii2] * d_A.y;
			f[ii2].z += weight[ii1] * lp * Nf[ii2] * d_A.z;
		}
	}

	return(NORMAL_RC);
}



/* $B#2<!85$^$?$O#1<!85MWAG(B element $B$KI=LLNO$,:nMQ$9$k>l9g$N(B */
/* $B@aE@2Y=E%Y%/%H%k$r5a$a$k!#(B*/
/* $B$?$@$7!"(Bp[p_size] $B$OI=LLNO$N@dBPCM!$(Bdir $B$OI=LLNO$NJ}8~(B */
/* p_size $B$,#1$N>l9g$OI=LLNO$,MWAGFb$G0lDj!"(B*/
/* P_size > 1 $B$N>l9g$O!"3F@aE@$G$NCM$H$7$FFbA^$9$k(B */
RC traction_force(int p_size, const double *p, TRANSLATION3D dir,
            FEM_ELEMENT *element, FEM_NODE_ARRAY node, TRANSLATION3D *f)
{
	static int allocate_flag = 0;
	static FEM_ELEM_TYPE latest_type;
	static int num_point;
	static double **node_location = NULL;
	static double **gauss_point = NULL;
	static double weight[MAX_GAUSS_POINT];
	static double **dN = NULL;
	static double **jacobi = NULL;
	double Np[FEM_MAX_NODE];
	double Nf[FEM_MAX_NODE];
	int ii1, ii2;
	FEM_ELEM_TYPE p_type = ELEM_TRI1;

	if(allocate_flag == 0){
		RC_TRY( allocate2D(FEM_MAX_NODE, 3, &node_location) );
		RC_TRY( allocate2D(MAX_GAUSS_POINT, 4, &gauss_point) );
		RC_TRY( allocate2D(3, FEM_MAX_NODE, &dN) );
		RC_TRY( allocate2D(3, 3, &jacobi) );

		RC_TRY( Gauss_const4pressure(element->type,
		                           gauss_point, weight, &num_point) );
		latest_type = element->type;
		allocate_flag = 1;
	}

	if(element->type != latest_type){
		RC_TRY( Gauss_const4pressure(element->type,
		                           gauss_point, weight, &num_point) );
		latest_type = element->type;
	}

	RC_TRY( node_location_matrix(element, node, node_location) );

	if(p_size != 1){
		switch(element->type){
		case ELEM_TRI1:
		case ELEM_TRI2:
			if(p_size == 3){
				p_type = ELEM_TRI1;
			}else if(p_size == 6){
				p_type = ELEM_TRI2;
			}else{
				return(IMPLEMENT_ERROR_RC);
			}
			break;
		case ELEM_QUAD1:
		case ELEM_QUAD2:
			if(p_size == 4){
				p_type = ELEM_QUAD1;
			}else if(p_size == 8){
				p_type = ELEM_QUAD2;
			}else{
				return(IMPLEMENT_ERROR_RC);
			}
			break;
		case ELEM_RBAR:
		case ELEM_BEAM1:
		case ELEM_BEAM2:
		case ELEM_LINE1:
		case ELEM_LINE2:
			if(p_size == 2){
				p_type = ELEM_LINE1;
			}else if(p_size == 3){
				p_type = ELEM_LINE2;
			}else{
				return(IMPLEMENT_ERROR_RC);
			}
			break;
		default:
			return(IMPLEMENT_ERROR_RC);
		}
	}

	for(ii1=0; ii1<element->node_num; ii1++){
		f[ii1].x = f[ii1].y = f[ii1].z = 0.0;
	}

	for(ii1=0; ii1<num_point; ii1++){
		double lp;

		init_matrix(3, FEM_MAX_NODE, dN);  /* $B=i4|2=$,I,MW(B */
		RC_TRY( set_dN_matrix(element->type, gauss_point[ii1], dN) );
		RC_TRY( set_N_vector(element->type, gauss_point[ii1], Nf) );
		mul_matrix(3, element->node_num, 3, dN, node_location, jacobi);

		/* p $B$rFbA^(B */
		if(p_size == 1){
			lp = p[0];
		}else{
			RC_TRY( set_N_vector(p_type, gauss_point[ii1], Np) );
			lp = 0.0;
			for(ii2=0; ii2 < p_size; ii2++){
				lp += p[ii2]* Np[ii2];
			}
		}

		for(ii2=0; ii2<element->node_num; ii2++){
			f[ii2].x += weight[ii1] * lp * Nf[ii2] * dir.x;
			f[ii2].y += weight[ii1] * lp * Nf[ii2] * dir.y;
			f[ii2].z += weight[ii1] * lp * Nf[ii2] * dir.z;
		}
	}

	return(NORMAL_RC);
}

