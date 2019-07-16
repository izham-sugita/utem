/*********************************************************************
 * operation.c.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: operation.c 442 2005-09-06 07:28:09Z sasaoka $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "rc.h"
#include "math_utl.h"


RC
allocate_operation_matrix (double ***matrix)
{
	RC_TRY( allocate2D(4, 4, matrix) );
	RC_TRY( init_operation_matrix(*matrix) );

	return(NORMAL_RC);
}


void
free_operation_matrix (double ***matrix)
{
	free2D(4, 4, matrix);
}


RC
init_operation_matrix (double **matrix)
{
	int ii1, ii2;

	RC_NULL_CHK(matrix);

	for(ii1=0; ii1<4; ii1++){
		for(ii2=0; ii2<4; ii2++){
			matrix[ii1][ii2] = 0.0;
		}
		matrix[ii1][ii1] = 1.0;
	}

	return(NORMAL_RC);
}

RC
print_operation_matrix (FILE *fp, double **matrix)
{
	int ii1;

	RC_NULL_CHK(matrix);
	for(ii1=0; ii1<4; ii1++){
		fprintf(fp, " %15.7e  %15.7e  %15.7e  %15.7e\n",
		             matrix[ii1][0], matrix[ii1][1],
		             matrix[ii1][2], matrix[ii1][3]);
	}

	return(NORMAL_RC);
}


RC
set_operation_matrix_pxpy (double fact, double **mat)
{
	mat[0][0] = fact; mat[0][1] = 0.0;  mat[0][2] = 0.0;  mat[0][3] = 0.0;
	mat[1][0] = 0.0;  mat[1][1] = fact; mat[1][2] = 0.0;  mat[1][3] = 0.0;
	mat[2][0] = 0.0;  mat[2][1] = 0.0;  mat[2][2] = fact; mat[2][3] = 0.0;
	mat[3][0] = 0.0;  mat[3][1] = 0.0;  mat[3][2] = 0.0;  mat[3][3] = 1.0;

	return(NORMAL_RC);
}


RC
set_operation_matrix_mzpy (double fact, double **mat)
{
	mat[0][0] = 0.0;  mat[0][1] = 0.0;  mat[0][2] =-fact; mat[0][3] = 0.0;
	mat[1][0] = 0.0;  mat[1][1] = fact; mat[1][2] = 0.0;  mat[1][3] = 0.0;
	mat[2][0] = fact; mat[2][1] = 0.0;  mat[2][2] = 0.0;  mat[2][3] = 0.0;
	mat[3][0] = 0.0;  mat[3][1] = 0.0;  mat[3][2] = 0.0;  mat[3][3] = 1.0;

	return(NORMAL_RC);
}


RC
set_operation_matrix_pxmz (double fact, double **mat)
{
	mat[0][0] = fact; mat[0][1] = 0.0;  mat[0][2] = 0.0;  mat[0][3] = 0.0;
	mat[1][0] = 0.0;  mat[1][1] = 0.0;  mat[1][2] =-fact; mat[1][3] = 0.0;
	mat[2][0] = 0.0;  mat[2][1] = fact; mat[2][2] = 0.0;  mat[2][3] = 0.0;
	mat[3][0] = 0.0;  mat[3][1] = 0.0;  mat[3][2] = 0.0;  mat[3][3] = 1.0;

	return(NORMAL_RC);
}


RC
shift_operation_matrix (double dx, double dy, double dz, double **mat)
{
	mat[0][0] += dx*mat[3][0];
	mat[0][1] += dx*mat[3][1];
	mat[0][2] += dx*mat[3][2];
	mat[0][3] += dx*mat[3][3];
	mat[1][0] += dy*mat[3][0];
	mat[1][1] += dy*mat[3][1];
	mat[1][2] += dy*mat[3][2];
	mat[1][3] += dy*mat[3][3];
	mat[2][0] += dz*mat[3][0];
	mat[2][1] += dz*mat[3][1];
	mat[2][2] += dz*mat[3][2];
	mat[2][3] += dz*mat[3][3];

	return(NORMAL_RC);
}


RC
scale_operation_matrix (double x, double y, double z, double **mat)
{
	int ii1;
		 
	for(ii1=0;ii1<4;ii1++){
		mat[0][ii1] *= x; 
		mat[1][ii1] *= y;
		mat[2][ii1] *= z;
	}

	return(NORMAL_RC);
}


RC
rotation_operation_matrix_x (double rad, double **mat)
{
	double m[2][4];
	double c,s;
	int ii1;

	for(ii1=0;ii1<4;ii1++){
		m[0][ii1] = mat[1][ii1];
		m[1][ii1] = mat[2][ii1];
	}
	c = cos(rad);
	s = sin(rad);

	mat[1][0] = c*m[0][0] - s*m[1][0];
	mat[1][1] = c*m[0][1] - s*m[1][1];
	mat[1][2] = c*m[0][2] - s*m[1][2];
	mat[1][3] = c*m[0][3] - s*m[1][3];
	mat[2][0] = s*m[0][0] + c*m[1][0];
	mat[2][1] = s*m[0][1] + c*m[1][1];
	mat[2][2] = s*m[0][2] + c*m[1][2];
	mat[2][3] = s*m[0][3] + c*m[1][3];

	return(NORMAL_RC);
}


RC
rotation_operation_matrix_y (double rad, double **mat)
{
	double m[2][4];
	double c,s;
	int ii1;

	for(ii1=0;ii1<4;ii1++){
		m[0][ii1] = mat[0][ii1];
		m[1][ii1] = mat[2][ii1];
	}
	c = cos(rad);
	s = sin(rad);

	mat[0][0] = c*m[0][0] + s*m[1][0];
	mat[0][1] = c*m[0][1] + s*m[1][1];
	mat[0][2] = c*m[0][2] + s*m[1][2];
	mat[0][3] = c*m[0][3] + s*m[1][3];
	mat[2][0] = -s*m[0][0] + c*m[1][0];
	mat[2][1] = -s*m[0][1] + c*m[1][1];
	mat[2][2] = -s*m[0][2] + c*m[1][2];
	mat[2][3] = -s*m[0][3] + c*m[1][3];

	return(NORMAL_RC);
}


RC
rotation_operation_matrix_z (double rad, double **mat)
{
	double m[2][4];
	double c,s;
	int ii1;

	for(ii1=0;ii1<4;ii1++){
		m[0][ii1] = mat[0][ii1];
		m[1][ii1] = mat[1][ii1];
	}
	c = cos(rad);
	s = sin(rad);

	mat[0][0] = c*m[0][0] - s*m[1][0];
	mat[0][1] = c*m[0][1] - s*m[1][1];
	mat[0][2] = c*m[0][2] - s*m[1][2];
	mat[0][3] = c*m[0][3] - s*m[1][3];
	mat[1][0] = s*m[0][0] + c*m[1][0];
	mat[1][1] = s*m[0][1] + c*m[1][1];
	mat[1][2] = s*m[0][2] + c*m[1][2];
	mat[1][3] = s*m[0][3] + c*m[1][3];

	return(NORMAL_RC);
}


/* 変換マトリックス B[4][4] に A[4][4] を作用させる */
RC
mul_operation_matrix (double **A, double **B)
{
	double m[4][4];
	int ii1, ii2;

	for(ii1=0;ii1<4;ii1++){
		for(ii2=0;ii2<4;ii2++){
			m[ii1][ii2] = B[ii1][ii2];
		}
	}

	for(ii1=0;ii1<4;ii1++){
		for(ii2=0;ii2<4;ii2++){
			B[ii1][ii2] = A[ii1][0] * m[0][ii2] + A[ii1][1] * m[1][ii2]
			            + A[ii1][2] * m[2][ii2] + A[ii1][3] * m[3][ii2];
		}
	}

	return(NORMAL_RC);
}


/* ベクトル v に変換マトリックス matrix[4][4] を作用させる */
VECT3D
operation3d (double **matrix, VECT3D v)
{
	VECT3D ret;

	ret.x = matrix[0][0]*v.x + matrix[0][1]*v.y
	      + matrix[0][2]*v.z + matrix[0][3];
	ret.y = matrix[1][0]*v.x + matrix[1][1]*v.y
	      + matrix[1][2]*v.z + matrix[1][3];
	ret.z = matrix[2][0]*v.x + matrix[2][1]*v.y
	      + matrix[2][2]*v.z + matrix[2][3];

	return(ret);
}

