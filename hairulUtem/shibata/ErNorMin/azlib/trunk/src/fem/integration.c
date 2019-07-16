/*********************************************************************
 * integration.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Mitsunori HIROSE> <Syoji ITO>
 *  <Tomoyuki TSUBATA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: integration.c 970 2010-08-25 03:10:08Z aoyama $ */

#include <stdio.h>
#include <math.h>
#include "math_utl.h"
#include "rc.h"


RC
Gauss_const_hexa (int number_of_point, double **xi_eta_zeta, double *weight)
{
	int inc;
	int ii1, ii2, ii3;
	double p[4];
	double w[4];

	switch(number_of_point){
	case 1:  /* m = 1 */
		xi_eta_zeta[0][0] = 0.000000000000000;
		xi_eta_zeta[0][1] = 0.000000000000000;
		xi_eta_zeta[0][2] = 0.000000000000000;
		weight[0] = 8.000000000000000;
		break;
	case 8:  /* m = 3 */
		xi_eta_zeta[0][0] = -0.577350269189626;
		xi_eta_zeta[0][1] = -0.577350269189626;
		xi_eta_zeta[0][2] = -0.577350269189626;
		weight[0] = 1.000000000000000;
		xi_eta_zeta[1][0] = +0.577350269189626;
		xi_eta_zeta[1][1] = -0.577350269189626;
		xi_eta_zeta[1][2] = -0.577350269189626;
		weight[1] = 1.000000000000000;
		xi_eta_zeta[2][0] = +0.577350269189626;
		xi_eta_zeta[2][1] = +0.577350269189626;
		xi_eta_zeta[2][2] = -0.577350269189626;
		weight[2] = 1.000000000000000;
		xi_eta_zeta[3][0] = -0.577350269189626;
		xi_eta_zeta[3][1] = +0.577350269189626;
		xi_eta_zeta[3][2] = -0.577350269189626;
		weight[3] = 1.000000000000000;
		xi_eta_zeta[4][0] = -0.577350269189626;
		xi_eta_zeta[4][1] = -0.577350269189626;
		xi_eta_zeta[4][2] = +0.577350269189626;
		weight[4] = 1.000000000000000;
		xi_eta_zeta[5][0] = +0.577350269189626;
		xi_eta_zeta[5][1] = -0.577350269189626;
		xi_eta_zeta[5][2] = +0.577350269189626;
		weight[5] = 1.000000000000000;
		xi_eta_zeta[6][0] = +0.577350269189626;
		xi_eta_zeta[6][1] = +0.577350269189626;
		xi_eta_zeta[6][2] = +0.577350269189626;
		weight[6] = 1.000000000000000;
		xi_eta_zeta[7][0] = -0.577350269189626;
		xi_eta_zeta[7][1] = +0.577350269189626;
		xi_eta_zeta[7][2] = +0.577350269189626;
		weight[7] = 1.000000000000000;
		break;
	case 14:  /* m = 5 */
		xi_eta_zeta[0][0] = -0.758686910639328;
		xi_eta_zeta[0][1] = -0.758686910639328;
		xi_eta_zeta[0][2] = -0.758686910639328;
		weight[0] = 0.335180055401662;
		xi_eta_zeta[1][0] = +0.758686910639328;
		xi_eta_zeta[1][1] = -0.758686910639328;
		xi_eta_zeta[1][2] = -0.758686910639328;
		weight[1] = 0.335180055401662;
		xi_eta_zeta[2][0] = +0.758686910639328;
		xi_eta_zeta[2][1] = +0.758686910639328;
		xi_eta_zeta[2][2] = -0.758686910639328;
		weight[2] = 0.335180055401662;
		xi_eta_zeta[3][0] = -0.758686910639328;
		xi_eta_zeta[3][1] = +0.758686910639328;
		xi_eta_zeta[3][2] = -0.758686910639328;
		weight[3] = 0.335180055401662;
		xi_eta_zeta[4][0] = -0.758686910639328;
		xi_eta_zeta[4][1] = -0.758686910639328;
		xi_eta_zeta[4][2] = +0.758686910639328;
		weight[4] = 0.335180055401662;
		xi_eta_zeta[5][0] = +0.758686910639328;
		xi_eta_zeta[5][1] = -0.758686910639328;
		xi_eta_zeta[5][2] = +0.758686910639328;
		weight[5] = 0.335180055401662;
		xi_eta_zeta[6][0] = +0.758686910639328;
		xi_eta_zeta[6][1] = +0.758686910639328;
		xi_eta_zeta[6][2] = +0.758686910639328;
		weight[6] = 0.335180055401662;
		xi_eta_zeta[7][0] = -0.758686910639328;
		xi_eta_zeta[7][1] = +0.758686910639328;
		xi_eta_zeta[7][2] = +0.758686910639328;
		weight[7] = 0.335180055401662;
		xi_eta_zeta[8][0] = 0.000000000000000;
		xi_eta_zeta[8][1] = 0.000000000000000;
		xi_eta_zeta[8][2] = -0.795822425754222;
		weight[8] = 0.886426592797784;
		xi_eta_zeta[9][0] = 0.000000000000000;
		xi_eta_zeta[9][1] = -0.795822425754222;
		xi_eta_zeta[9][2] = 0.000000000000000;
		weight[9] = 0.886426592797784;
		xi_eta_zeta[10][0] = +0.795822425754222;
		xi_eta_zeta[10][1] = 0.000000000000000;
		xi_eta_zeta[10][2] = 0.000000000000000;
		weight[10] = 0.886426592797784;
		xi_eta_zeta[11][0] = 0.000000000000000;
		xi_eta_zeta[11][1] = +0.795822425754222;
		xi_eta_zeta[11][2] = 0.000000000000000;
		weight[11] = 0.886426592797784;
		xi_eta_zeta[12][0] = -0.795822425754222;
		xi_eta_zeta[12][1] = 0.000000000000000;
		xi_eta_zeta[12][2] = 0.000000000000000;
		weight[12] = 0.886426592797784;
		xi_eta_zeta[13][0] = 0.000000000000000;
		xi_eta_zeta[13][1] = 0.000000000000000;
		xi_eta_zeta[13][2] = +0.795822425754222;
		weight[13] = 0.886426592797784;
		break;
	case 27:
		p[0] = -0.774596669241483;
		p[1] = +0.774596669241483;
		p[2] =  0.000000000000000;
		w[0] = 0.555555555555556;
		w[1] = 0.555555555555556;
		w[2] = 2.0 - (w[0] + w[1]);
		inc = 0;
		for(ii1=0; ii1<3; ii1++){
			for(ii2=0; ii2<3; ii2++){
				for(ii3=0; ii3<3; ii3++){
					xi_eta_zeta[inc][0] = p[ii1];
					xi_eta_zeta[inc][1] = p[ii2];
					xi_eta_zeta[inc][2] = p[ii3];
					weight[inc] = w[ii1] * w[ii2] * w[ii3];
					inc++;
				}
			}
		}
		break;
	case 64:
		p[0] = -0.861136311594053;
		p[1] = +0.861136311594053;
		p[2] = -0.339981043584856;
		p[3] = +0.339981043584856;
		w[0] = 0.347854845137454;
		w[1] = 0.347854845137454;
		w[2] = 0.652145154862546;
		w[3] = 0.652145154862546;
		inc = 0;
		for(ii1=0; ii1<4; ii1++){
			for(ii2=0; ii2<4; ii2++){
				for(ii3=0; ii3<4; ii3++){
					xi_eta_zeta[inc][0] = p[ii1];
					xi_eta_zeta[inc][1] = p[ii2];
					xi_eta_zeta[inc][2] = p[ii3];
					weight[inc] = w[ii1] * w[ii2] * w[ii3];
					inc++;
				}
			}
		}
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


RC
Gauss_const_penta(int number_of_point, double **L1_L2_L3_zeta, double *weight)
{
	double z1, z2, zw1, zw2;
	double a, b, w1, w2, w3;

	switch(number_of_point){
	case 1:  /* m = 1 */
		L1_L2_L3_zeta[0][0] = 0.333333333333333;
		L1_L2_L3_zeta[0][1] = 0.333333333333333;
		L1_L2_L3_zeta[0][2] = 0.333333333333333;
		L1_L2_L3_zeta[0][3] = 0.000000000000000;
		weight[0] = 2.000000000000000/2.0;
		break;
	case 6:  /* m = 2, 3 */
		L1_L2_L3_zeta[0][0] = 0.666666666666667;
		L1_L2_L3_zeta[0][1] = 0.166666666666667;
		L1_L2_L3_zeta[0][2] = 0.166666666666667;
		L1_L2_L3_zeta[0][3] = -0.577350269189626;
		weight[0] = 0.333333333333333/2.0;
		L1_L2_L3_zeta[1][0] = 0.166666666666667;
		L1_L2_L3_zeta[1][1] = 0.666666666666667;
		L1_L2_L3_zeta[1][2] = 0.166666666666667;
		L1_L2_L3_zeta[1][3] = -0.577350269189626;
		weight[1] = 0.333333333333333/2.0;
		L1_L2_L3_zeta[2][0] = 0.166666666666667;
		L1_L2_L3_zeta[2][1] = 0.166666666666667;
		L1_L2_L3_zeta[2][2] = 0.666666666666667;
		L1_L2_L3_zeta[2][3] = -0.577350269189626;
		weight[2] = 0.333333333333333/2.0;
		L1_L2_L3_zeta[3][0] = 0.166666666666667;
		L1_L2_L3_zeta[3][1] = 0.166666666666667;
		L1_L2_L3_zeta[3][2] = 0.666666666666667;
		L1_L2_L3_zeta[3][3] = 0.577350269189626;
		weight[3] = 0.333333333333333/2.0;
		L1_L2_L3_zeta[4][0] = 0.166666666666667;
		L1_L2_L3_zeta[4][1] = 0.166666666666667;
		L1_L2_L3_zeta[4][2] = 0.666666666666667;
		L1_L2_L3_zeta[4][3] = 0.577350269189626;
		weight[4] = 0.333333333333333/2.0;
		L1_L2_L3_zeta[5][0] = 0.166666666666667;
		L1_L2_L3_zeta[5][1] = 0.166666666666667;
		L1_L2_L3_zeta[5][2] = 0.666666666666667;
		L1_L2_L3_zeta[5][3] = 0.577350269189626;
		weight[5] = 0.333333333333333/2.0;
		break;
	case 9:  /* m = 2, 5 */
		L1_L2_L3_zeta[0][0] = 0.666666666666667;
		L1_L2_L3_zeta[0][1] = 0.166666666666667;
		L1_L2_L3_zeta[0][2] = 0.166666666666667;
		L1_L2_L3_zeta[0][3] = -0.774596669241483;
		weight[0] = 0.555555555555556/6.0;
		L1_L2_L3_zeta[1][0] = 0.166666666666667;
		L1_L2_L3_zeta[1][1] = 0.666666666666667;
		L1_L2_L3_zeta[1][2] = 0.166666666666667;
		L1_L2_L3_zeta[1][3] = -0.774596669241483;
		weight[1] = 0.555555555555556/6.0;
		L1_L2_L3_zeta[2][0] = 0.166666666666667;
		L1_L2_L3_zeta[2][1] = 0.166666666666667;
		L1_L2_L3_zeta[2][2] = 0.666666666666667;
		L1_L2_L3_zeta[2][3] = -0.774596669241483;
		weight[2] = 0.555555555555556/6.0;
		L1_L2_L3_zeta[3][0] = 0.166666666666667;
		L1_L2_L3_zeta[3][1] = 0.166666666666667;
		L1_L2_L3_zeta[3][2] = 0.666666666666667;
		L1_L2_L3_zeta[3][3] = 0.000000000000000;
		weight[3] = 0.888888888888889/6.0;
		L1_L2_L3_zeta[4][0] = 0.166666666666667;
		L1_L2_L3_zeta[4][1] = 0.166666666666667;
		L1_L2_L3_zeta[4][2] = 0.666666666666667;
		L1_L2_L3_zeta[4][3] = 0.000000000000000;
		weight[4] = 0.888888888888889/6.0;
		L1_L2_L3_zeta[5][0] = 0.166666666666667;
		L1_L2_L3_zeta[5][1] = 0.166666666666667;
		L1_L2_L3_zeta[5][2] = 0.666666666666667;
		L1_L2_L3_zeta[5][3] = 0.000000000000000;
		weight[5] = 0.888888888888889/6.0;
		L1_L2_L3_zeta[6][0] = 0.666666666666667;
		L1_L2_L3_zeta[6][1] = 0.166666666666667;
		L1_L2_L3_zeta[6][2] = 0.166666666666667;
		L1_L2_L3_zeta[6][3] = +0.774596669241483;
		weight[6] = 0.555555555555556/6.0;
		L1_L2_L3_zeta[7][0] = 0.166666666666667;
		L1_L2_L3_zeta[7][1] = 0.666666666666667;
		L1_L2_L3_zeta[7][2] = 0.166666666666667;
		L1_L2_L3_zeta[7][3] = +0.774596669241483;
		weight[7] = 0.555555555555556/6.0;
		L1_L2_L3_zeta[8][0] = 0.166666666666667;
		L1_L2_L3_zeta[8][1] = 0.166666666666667;
		L1_L2_L3_zeta[8][2] = 0.666666666666667;
		L1_L2_L3_zeta[8][3] = +0.774596669241483;
		weight[8] = 0.555555555555556/6.0;
	case 21:  /* m = 5, 5 */
		z1 = 0.774596669241483;
		z2 = 0.000000000000000;
		zw1 = 0.555555555555556;
		zw2 = 0.888888888888889;
		a = (6.0+sqrt(15.0))/21.0;
		b = 4.0/7.0 - a;
		w1 = 9.0/80.0;
		w2 = (155.0 + sqrt(15.0))/2400.0;
		w3 = 31.0/240.0 - w2;
		L1_L2_L3_zeta[0][0] = 1.0/3.0;
		L1_L2_L3_zeta[0][1] = 1.0/3.0;
		L1_L2_L3_zeta[0][2] = 1.0/3.0;
		L1_L2_L3_zeta[0][3] = z1;
		weight[0] = w1*zw1;
		L1_L2_L3_zeta[1][0] = a;
		L1_L2_L3_zeta[1][1] = a;
		L1_L2_L3_zeta[1][2] = 1.0 - 2.0*a;
		L1_L2_L3_zeta[1][3] = z1;
		weight[1] = w2*zw1;
		L1_L2_L3_zeta[2][0] = 1.0 - 2.0*a;
		L1_L2_L3_zeta[2][1] = a;
		L1_L2_L3_zeta[2][2] = a;
		L1_L2_L3_zeta[2][3] = z1;
		weight[2] = w2*zw1;
		L1_L2_L3_zeta[3][0] = a;
		L1_L2_L3_zeta[3][1] = 1.0 - 2.0*a;
		L1_L2_L3_zeta[3][2] = a;
		L1_L2_L3_zeta[3][3] = z1;
		weight[3] = w2*zw1;
		L1_L2_L3_zeta[4][0] = b;
		L1_L2_L3_zeta[4][1] = b;
		L1_L2_L3_zeta[4][2] = 1.0 - 2.0*b;
		L1_L2_L3_zeta[4][3] = z1;
		weight[4] = w3*zw1;
		L1_L2_L3_zeta[5][0] = 1.0 - 2.0*b;
		L1_L2_L3_zeta[5][1] = b;
		L1_L2_L3_zeta[5][2] = b;
		L1_L2_L3_zeta[5][3] = z1;
		weight[5] = w3*zw1;
		L1_L2_L3_zeta[6][0] = b;
		L1_L2_L3_zeta[6][1] = 1.0 - 2.0*b;
		L1_L2_L3_zeta[6][2] = b;
		L1_L2_L3_zeta[6][3] = z1;
		weight[6] = w3*zw1;
		L1_L2_L3_zeta[7][0] = 1.0/3.0;
		L1_L2_L3_zeta[7][1] = 1.0/3.0;
		L1_L2_L3_zeta[7][2] = 1.0/3.0;
		L1_L2_L3_zeta[7][3] = z2;
		weight[7] = w1*zw2;
		L1_L2_L3_zeta[8][0] = a;
		L1_L2_L3_zeta[8][1] = a;
		L1_L2_L3_zeta[8][2] = 1.0 - 2.0*a;
		L1_L2_L3_zeta[8][3] = z2;
		weight[8] = w2*zw2;
		L1_L2_L3_zeta[9][0] = 1.0 - 2.0*a;
		L1_L2_L3_zeta[9][1] = a;
		L1_L2_L3_zeta[9][2] = a;
		L1_L2_L3_zeta[9][3] = z2;
		weight[9] = w2*zw2;
		L1_L2_L3_zeta[10][0] = a;
		L1_L2_L3_zeta[10][1] = 1.0 - 2.0*a;
		L1_L2_L3_zeta[10][2] = a;
		L1_L2_L3_zeta[10][3] = z2;
		weight[10] = w2*zw2;
		L1_L2_L3_zeta[11][0] = b;
		L1_L2_L3_zeta[11][1] = b;
		L1_L2_L3_zeta[11][2] = 1.0 - 2.0*b;
		L1_L2_L3_zeta[11][3] = z2;
		weight[11] = w3*zw2;
		L1_L2_L3_zeta[12][0] = 1.0 - 2.0*b;
		L1_L2_L3_zeta[12][1] = b;
		L1_L2_L3_zeta[12][2] = b;
		L1_L2_L3_zeta[12][3] = z2;
		weight[12] = w3*zw2;
		L1_L2_L3_zeta[13][0] = b;
		L1_L2_L3_zeta[13][1] = 1.0 - 2.0*b;
		L1_L2_L3_zeta[13][2] = b;
		L1_L2_L3_zeta[13][3] = z2;
		weight[13] = w3*zw2;
		L1_L2_L3_zeta[14][0] = 1.0/3.0;
		L1_L2_L3_zeta[14][1] = 1.0/3.0;
		L1_L2_L3_zeta[14][2] = 1.0/3.0;
		L1_L2_L3_zeta[14][3] = -z1;
		weight[14] = w1*zw1;
		L1_L2_L3_zeta[15][0] = a;
		L1_L2_L3_zeta[15][1] = a;
		L1_L2_L3_zeta[15][2] = 1.0 - 2.0*a;
		L1_L2_L3_zeta[15][3] = -z1;
		weight[15] = w2*zw1;
		L1_L2_L3_zeta[16][0] = 1.0 - 2.0*a;
		L1_L2_L3_zeta[16][1] = a;
		L1_L2_L3_zeta[16][2] = a;
		L1_L2_L3_zeta[16][3] = -z1;
		weight[16] = w2*zw1;
		L1_L2_L3_zeta[17][0] = a;
		L1_L2_L3_zeta[17][1] = 1.0 - 2.0*a;
		L1_L2_L3_zeta[17][2] = a;
		L1_L2_L3_zeta[17][3] = -z1;
		weight[17] = w2*zw1;
		L1_L2_L3_zeta[18][0] = b;
		L1_L2_L3_zeta[18][1] = b;
		L1_L2_L3_zeta[18][2] = 1.0 - 2.0*b;
		L1_L2_L3_zeta[18][3] = -z1;
		weight[18] = w3*zw1;
		L1_L2_L3_zeta[19][0] = 1.0 - 2.0*b;
		L1_L2_L3_zeta[19][1] = b;
		L1_L2_L3_zeta[19][2] = b;
		L1_L2_L3_zeta[19][3] = -z1;
		weight[19] = w3*zw1;
		L1_L2_L3_zeta[20][0] = b;
		L1_L2_L3_zeta[20][1] = 1.0 - 2.0*b;
		L1_L2_L3_zeta[20][2] = b;
		L1_L2_L3_zeta[20][3] = -z1;
		weight[20] = w3*zw1;
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


RC
Gauss_const_tetra (int number_of_point, double **L1_L2_L3_L4, double *weight)
{
	double a, b1, b2, c1, c2, d, e;
	double w0, w1, w2, w3;

	switch(number_of_point){
	case 1:  /* m = 1 */
		L1_L2_L3_L4[0][0] = 0.250000000000000;
		L1_L2_L3_L4[0][1] = 0.250000000000000;
		L1_L2_L3_L4[0][2] = 0.250000000000000;
		L1_L2_L3_L4[0][3] = 0.250000000000000;
		weight[0] = 1.000000000000000/6.0;
		break;
	case 4:  /* m = 2 */
		L1_L2_L3_L4[0][0] = 0.585410196624968;
		L1_L2_L3_L4[0][1] = 0.138196601125010;
		L1_L2_L3_L4[0][2] = 0.138196601125010;
		L1_L2_L3_L4[0][3] = 0.138196601125010;
		weight[0] = 0.250000000000000/6.0;
		L1_L2_L3_L4[1][0] = 0.138196601125010;
		L1_L2_L3_L4[1][1] = 0.585410196624968;
		L1_L2_L3_L4[1][2] = 0.138196601125010;
		L1_L2_L3_L4[1][3] = 0.138196601125010;
		weight[1] = 0.250000000000000/6.0;
		L1_L2_L3_L4[2][0] = 0.138196601125010;
		L1_L2_L3_L4[2][1] = 0.138196601125010;
		L1_L2_L3_L4[2][2] = 0.585410196624968;
		L1_L2_L3_L4[2][3] = 0.138196601125010;
		weight[2] = 0.250000000000000/6.0;
		L1_L2_L3_L4[3][0] = 0.138196601125010;
		L1_L2_L3_L4[3][1] = 0.138196601125010;
		L1_L2_L3_L4[3][2] = 0.138196601125010;
		L1_L2_L3_L4[3][3] = 0.585410196624968;
		weight[3] = 0.250000000000000/6.0;
		break;
	case 5:  /* m = 3 */
		L1_L2_L3_L4[0][0] = 0.500000000000000;
		L1_L2_L3_L4[0][1] = 0.166666666666666;
		L1_L2_L3_L4[0][2] = 0.166666666666666;
		L1_L2_L3_L4[0][3] = 0.166666666666666;
		weight[0] = 0.450000000000000/6.0;
		L1_L2_L3_L4[1][0] = 0.166666666666666;
		L1_L2_L3_L4[1][1] = 0.500000000000000;
		L1_L2_L3_L4[1][2] = 0.166666666666666;
		L1_L2_L3_L4[1][3] = 0.166666666666666;
		weight[1] = 0.450000000000000/6.0;
		L1_L2_L3_L4[2][0] = 0.166666666666666;
		L1_L2_L3_L4[2][1] = 0.166666666666666;
		L1_L2_L3_L4[2][2] = 0.500000000000000;
		L1_L2_L3_L4[2][3] = 0.166666666666666;
		weight[2] = 0.450000000000000/6.0;
		L1_L2_L3_L4[3][0] = 0.166666666666666;
		L1_L2_L3_L4[3][1] = 0.166666666666666;
		L1_L2_L3_L4[3][2] = 0.166666666666666;
		L1_L2_L3_L4[3][3] = 0.500000000000000;
		weight[3] = 0.450000000000000/6.0;
		L1_L2_L3_L4[4][0] = 0.250000000000000;
		L1_L2_L3_L4[4][1] = 0.250000000000000;
		L1_L2_L3_L4[4][2] = 0.250000000000000;
		L1_L2_L3_L4[4][3] = 0.250000000000000;
		weight[4] = -0.800000000000000/6.0;
		break;
	case 15:  /* m = 5 */
		a = 1.0/4.0;
		w0 = 112.0/5670.0;
		b1 = (7.0 - sqrt(15.0))/34.0;
		b2 = (7.0 + sqrt(15.0))/34.0;
		c1 = (13.0 + 3.0*sqrt(15.0))/34.0;
		c2 = (13.0 - 3.0*sqrt(15.0))/34.0;
		w1 = (2665.0+14.0*sqrt(15.0))/226800.0;
		w2 = (2665.0-14.0*sqrt(15.0))/226800.0;
		d = (5.0 - sqrt(15.0))/20.0;
		e = (5.0 + sqrt(15.0))/20.0;
		w3 = 5.0/567.0;
		L1_L2_L3_L4[0][0] = a;
		L1_L2_L3_L4[0][1] = a;
		L1_L2_L3_L4[0][2] = a;
		L1_L2_L3_L4[0][3] = a;
		weight[0] = w0;
		L1_L2_L3_L4[1][0] = b1;
		L1_L2_L3_L4[1][1] = b1;
		L1_L2_L3_L4[1][2] = b1;
		L1_L2_L3_L4[1][3] = c1;
		weight[1] = w1;
		L1_L2_L3_L4[2][0] = b1;
		L1_L2_L3_L4[2][1] = b1;
		L1_L2_L3_L4[2][2] = c1;
		L1_L2_L3_L4[2][3] = b1;
		weight[2] = w1;
		L1_L2_L3_L4[3][0] = b1;
		L1_L2_L3_L4[3][1] = c1;
		L1_L2_L3_L4[3][2] = b1;
		L1_L2_L3_L4[3][3] = b1;
		weight[3] = w1;
		L1_L2_L3_L4[4][0] = c1;
		L1_L2_L3_L4[4][1] = b1;
		L1_L2_L3_L4[4][2] = b1;
		L1_L2_L3_L4[4][3] = b1;
		weight[4] = w1;
		L1_L2_L3_L4[5][0] = b2;
		L1_L2_L3_L4[5][1] = b2;
		L1_L2_L3_L4[5][2] = b2;
		L1_L2_L3_L4[5][3] = c2;
		weight[5] = w2;
		L1_L2_L3_L4[6][0] = b2;
		L1_L2_L3_L4[6][1] = b2;
		L1_L2_L3_L4[6][2] = c2;
		L1_L2_L3_L4[6][3] = b2;
		weight[6] = w2;
		L1_L2_L3_L4[7][0] = b2;
		L1_L2_L3_L4[7][1] = c2;
		L1_L2_L3_L4[7][2] = b2;
		L1_L2_L3_L4[7][3] = b2;
		weight[7] = w2;
		L1_L2_L3_L4[8][0] = c2;
		L1_L2_L3_L4[8][1] = b2;
		L1_L2_L3_L4[8][2] = b2;
		L1_L2_L3_L4[8][3] = b2;
		weight[8] = w2;
		L1_L2_L3_L4[9][0] = d;
		L1_L2_L3_L4[9][1] = d;
		L1_L2_L3_L4[9][2] = e;
		L1_L2_L3_L4[9][3] = e;
		weight[9] = w3;
		L1_L2_L3_L4[10][0] = d;
		L1_L2_L3_L4[10][1] = e;
		L1_L2_L3_L4[10][2] = d;
		L1_L2_L3_L4[10][3] = e;
		weight[10] = w3;
		L1_L2_L3_L4[11][0] = e;
		L1_L2_L3_L4[11][1] = d;
		L1_L2_L3_L4[11][2] = d;
		L1_L2_L3_L4[11][3] = e;
		weight[11] = w3;
		L1_L2_L3_L4[12][0] = d;
		L1_L2_L3_L4[12][1] = e;
		L1_L2_L3_L4[12][2] = e;
		L1_L2_L3_L4[12][3] = d;
		weight[12] = w3;
		L1_L2_L3_L4[13][0] = e;
		L1_L2_L3_L4[13][1] = d;
		L1_L2_L3_L4[13][2] = e;
		L1_L2_L3_L4[13][3] = d;
		weight[13] = w3;
		L1_L2_L3_L4[14][0] = e;
		L1_L2_L3_L4[14][1] = e;
		L1_L2_L3_L4[14][2] = d;
		L1_L2_L3_L4[14][3] = d;
		weight[14] = w3;
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


RC
Gauss_const_quad (int number_of_point, double **xi_eta, double *weight)
{	
	int inc;
	int ii1, ii2;
	double p[6];
	double w[6];

	switch(number_of_point){
	case 1:  /* m = 1 */
		xi_eta[0][0] = 0.000000000000000;
		xi_eta[0][1] = 0.000000000000000;
		weight[0] = 4.000000000000000;
		break;
	case 4:  /* m = 3 */
		xi_eta[0][0] = -0.577350269189626;
		xi_eta[0][1] = -0.577350269189626;
		weight[0] = 1.000000000000000;
		xi_eta[1][0] = +0.577350269189626;
		xi_eta[1][1] = -0.577350269189626;
		weight[1] = 1.000000000000000;
		xi_eta[2][0] = +0.577350269189626;
		xi_eta[2][1] = +0.577350269189626;
		weight[2] = 1.000000000000000;
		xi_eta[3][0] = -0.577350269189626;
		xi_eta[3][1] = +0.577350269189626;
		weight[3] = 1.000000000000000;
		break;
	case 9:  /* m = 5 */
		xi_eta[0][0] = -0.774596669241483;
		xi_eta[0][1] = -0.774596669241483;
		weight[0] = 0.308641975308642;
		xi_eta[1][0] = +0.774596669241483;
		xi_eta[1][1] = -0.774596669241483;
		weight[1] = 0.308641975308642;
		xi_eta[2][0] = +0.774596669241483;
		xi_eta[2][1] = +0.774596669241483;
		weight[2] = 0.308641975308642;
		xi_eta[3][0] = -0.774596669241483;
		xi_eta[3][1] = +0.774596669241483;
		weight[3] = 0.308641975308642;
		xi_eta[4][0] = 0.000000000000000;
		xi_eta[4][1] = -0.774596669241483;
		weight[4] = 0.493827160493827;
		xi_eta[5][0] = +0.774596669241483;
		xi_eta[5][1] = 0.000000000000000;
		weight[5] = 0.493827160493827;
		xi_eta[6][0] = 0.000000000000000;
		xi_eta[6][1] = +0.774596669241483;
		weight[6] = 0.493827160493827;
		xi_eta[7][0] = -0.774596669241483;
		xi_eta[7][1] = 0.000000000000000;
		weight[7] = 0.493827160493827;
		xi_eta[8][0] = 0.000000000000000;
		xi_eta[8][1] = 0.000000000000000;
		weight[8] = 0.790123456790123;
		break;
	case 16:
		p[0] = -0.861136311594053;
		p[1] = +0.861136311594053;
		p[2] = -0.339981043584856;
		p[3] = +0.339981043584856;
		w[0] = 0.347854845137454;
		w[1] = 0.347854845137454;
		w[2] = 0.652145154862546;
		w[3] = 0.652145154862546;
		inc = 0;
		for(ii1=0; ii1<4; ii1++){
			for(ii2=0; ii2<4; ii2++){
				xi_eta[inc][0] = p[ii1];
				xi_eta[inc][1] = p[ii2];
				weight[inc] = w[ii1] * w[ii2];
				inc++;
			}
		}
		break;
	case 25:
		p[0] = -0.906179845938664;
		p[1] = +0.906179845938664;
		p[2] = -0.538469310105683;
		p[3] = +0.538469310105683;
		p[4] = 0.000000000000000;
		w[0] = 0.236926885056189;
		w[1] = 0.236926885056189;
		w[2] = 0.478628670499366;
		w[3] = 0.478628670499366;
		w[4] = 0.568888888888889;
		inc = 0;
		for(ii1=0; ii1<5; ii1++){
			for(ii2=0; ii2<5; ii2++){
				xi_eta[inc][0] = p[ii1];
				xi_eta[inc][1] = p[ii2];
				weight[inc] = w[ii1] * w[ii2];
				inc++;
			}
		}
		break;
	case 36:
		p[0] = -0.932469514203152;
		p[1] = +0.932469514203152;
		p[2] = -0.661209386466265;
		p[3] = +0.661209386466265;
		p[4] = -0.238619186083197;
		p[5] = +0.238619186083197;
		w[0] = 0.171324492379170;
		w[1] = 0.171324492379170;
		w[2] = 0.360761573048139;
		w[3] = 0.360761573048139;
		w[4] = 0.467913934572691;
		w[5] = 0.467913934572691;
		inc = 0;
		for(ii1=0; ii1<6; ii1++){
			for(ii2=0; ii2<6; ii2++){
				xi_eta[inc][0] = p[ii1];
				xi_eta[inc][1] = p[ii2];
				weight[inc] = w[ii1] * w[ii2];
				inc++;
			}
		}
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


RC
Gauss_const_tri (int number_of_point, double **L1_L2_L3, double *weight)
{
	double a, b, c, d, w1, w2, w3;

	switch(number_of_point){
	case 1:   /* m = 1 */
		L1_L2_L3[0][0] = 0.333333333333333;
		L1_L2_L3[0][1] = 0.333333333333333;
		L1_L2_L3[0][2] = 0.333333333333333;
		weight[0] = 1.000000000000000/2.0;
		break;
	case 3:   /* m = 2 */
		L1_L2_L3[0][0] = 0.666666666666667;
		L1_L2_L3[0][1] = 0.166666666666667;
		L1_L2_L3[0][2] = 0.166666666666667;
		weight[0] = 0.333333333333333/2.0;
		L1_L2_L3[1][0] = 0.166666666666667;
		L1_L2_L3[1][1] = 0.666666666666667;
		L1_L2_L3[1][2] = 0.166666666666667;
		weight[1] = 0.333333333333333/2.0;
		L1_L2_L3[2][0] = 0.166666666666667;
		L1_L2_L3[2][1] = 0.166666666666667;
		L1_L2_L3[2][2] = 0.666666666666667;
		weight[2] = 0.333333333333333/2.0;
		break;
	case 4:   /* m = 3 */
		L1_L2_L3[0][0] = 1.0/3.0;
		L1_L2_L3[0][1] = 1.0/3.0;
		L1_L2_L3[0][2] = 1.0/3.0;
		weight[0] = -27.0/96.0;
		L1_L2_L3[1][0] = 1.0/5.0;
		L1_L2_L3[1][1] = 1.0/5.0;
		L1_L2_L3[1][2] = 3.0/5.0;
		weight[1] = 25.0/96.0;
		L1_L2_L3[2][0] = 3.0/5.0;
		L1_L2_L3[2][1] = 1.0/5.0;
		L1_L2_L3[2][2] = 1.0/5.0;
		weight[2] = 25.0/96.0;
		L1_L2_L3[3][0] = 1.0/5.0;
		L1_L2_L3[3][1] = 3.0/5.0;
		L1_L2_L3[3][2] = 1.0/5.0;
		weight[3] = 25.0/96.0;
		break;
	case 6:   /* m = 4 */
		a = 0.445948490915965;
		b = 0.091576213509771;
		w1 = 0.111690794839005;
		w2 = 0.054975871827661;
		L1_L2_L3[0][0] = a;
		L1_L2_L3[0][1] = a;
		L1_L2_L3[0][2] = 1.0 - 2.0*a;
		weight[0] = w1;
		L1_L2_L3[1][0] = 1.0 - 2.0*a;
		L1_L2_L3[1][1] = a;
		L1_L2_L3[1][2] = a;
		weight[1] = w1;
		L1_L2_L3[2][0] = a;
		L1_L2_L3[2][1] = 1.0 - 2.0*a;
		L1_L2_L3[2][2] = a;
		weight[2] = w1;
		L1_L2_L3[3][0] = b;
		L1_L2_L3[3][1] = b;
		L1_L2_L3[3][2] = 1.0 - 2.0*b;
		weight[3] = w2;
		L1_L2_L3[4][0] = 1.0 - 2.0*b;
		L1_L2_L3[4][1] = b;
		L1_L2_L3[4][2] = b;
		weight[4] = w2;
		L1_L2_L3[5][0] = b;
		L1_L2_L3[5][1] = 1.0 - 2.0*b;
		L1_L2_L3[5][2] = b;
		weight[5] = w2;
		break;
	case 7:   /* m = 5 */
		a = (6.0+sqrt(15.0))/21.0;
		b = 4.0/7.0 - a;
		w1 = 9.0/80.0;
		w2 = (155.0 + sqrt(15.0))/2400.0;
		w3 = 31.0/240.0 - w2;
		L1_L2_L3[0][0] = 1.0/3.0;
		L1_L2_L3[0][1] = 1.0/3.0;
		L1_L2_L3[0][2] = 1.0/3.0;
		weight[0] = w1;
		L1_L2_L3[1][0] = a;
		L1_L2_L3[1][1] = a;
		L1_L2_L3[1][2] = 1.0 - 2.0*a;
		weight[1] = w2;
		L1_L2_L3[2][0] = 1.0 - 2.0*a;
		L1_L2_L3[2][1] = a;
		L1_L2_L3[2][2] = a;
		weight[2] = w2;
		L1_L2_L3[3][0] = a;
		L1_L2_L3[3][1] = 1.0 - 2.0*a;
		L1_L2_L3[3][2] = a;
		weight[3] = w2;
		L1_L2_L3[4][0] = b;
		L1_L2_L3[4][1] = b;
		L1_L2_L3[4][2] = 1.0 - 2.0*b;
		weight[4] = w3;
		L1_L2_L3[5][0] = 1.0 - 2.0*b;
		L1_L2_L3[5][1] = b;
		L1_L2_L3[5][2] = b;
		weight[5] = w3;
		L1_L2_L3[6][0] = b;
		L1_L2_L3[6][1] = 1.0 - 2.0*b;
		L1_L2_L3[6][2] = b;
		weight[6] = w3;
		break;
	case 12:   /* m = 6 */
		a = 0.063089014491502;
		b = 0.249286745170910;
		c = 0.310352451033785;
		d = 0.053145049844816;
		w1 = 0.025422453185103;
		w2 = 0.058393137863189;
		w3 = 0.041425537809187;
		L1_L2_L3[0][0] = a;
		L1_L2_L3[0][1] = a;
		L1_L2_L3[0][2] = 1.0 - 2.0*a;
		weight[0] = w1;
		L1_L2_L3[1][0] = 1.0 - 2.0*a;
		L1_L2_L3[1][1] = a;
		L1_L2_L3[1][2] = a;
		weight[1] = w1;
		L1_L2_L3[2][0] = a;
		L1_L2_L3[2][1] = 1.0 - 2.0*a;
		L1_L2_L3[2][2] = a;
		weight[2] = w1;
		L1_L2_L3[3][0] = b;
		L1_L2_L3[3][1] = b;
		L1_L2_L3[3][2] = 1.0 - 2.0*b;
		weight[3] = w2;
		L1_L2_L3[4][0] = 1.0 - 2.0*b;
		L1_L2_L3[4][1] = b;
		L1_L2_L3[4][2] = b;
		weight[4] = w2;
		L1_L2_L3[5][0] = b;
		L1_L2_L3[5][1] = 1.0 - 2.0*b;
		L1_L2_L3[5][2] = b;
		weight[5] = w2;
		L1_L2_L3[6][0] = c;
		L1_L2_L3[6][1] = d;
		L1_L2_L3[6][2] = 1.0 - (c + d);
		weight[6] = w3;
		L1_L2_L3[7][0] = d;
		L1_L2_L3[7][1] = c;
		L1_L2_L3[7][2] = 1.0 - (c + d);
		weight[7] = w3;
		L1_L2_L3[8][0] = 1.0 - (c + d);
		L1_L2_L3[8][1] = c;
		L1_L2_L3[8][2] = d;
		weight[8] = w3;
		L1_L2_L3[9][0] = 1.0 - (c + d);
		L1_L2_L3[9][1] = d;
		L1_L2_L3[9][2] = c;
		weight[9] = w3;
		L1_L2_L3[10][0] = c;
		L1_L2_L3[10][1] = 1.0 - (c + d);
		L1_L2_L3[10][2] = d;
		weight[10] = w3;
		L1_L2_L3[11][0] = d;
		L1_L2_L3[11][1] = 1.0 - (c + d);
		L1_L2_L3[11][2] = c;
		weight[11] = w3;
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


RC
Gauss_const_line (int number_of_point, double **xi, double *weight)
 {
	switch(number_of_point){
	case 1:  /* m = 1 */
		xi[0][0] = 0.000000000000000;
		weight[0] = 2.000000000000000;
		break;
	case 2:  /* m = 3 */
		xi[0][0] = -0.577350269189626;
		weight[0] = 1.000000000000000;
		xi[1][0] = +0.577350269189626;
		weight[1] = 1.000000000000000;
		break;
	case 3:  /* m = 5 */
		xi[0][0] = -0.774596669241483;
		weight[0] = 0.555555555555556;
		xi[1][0] = +0.774596669241483;
		weight[1] = 0.555555555555556;
		xi[2][0] = 0.000000000000000;
		weight[2] = 0.888888888888889;
		break;
	case 4:  /* m = 7 */
		xi[0][0] = -0.861136311594053;
		weight[0] = 0.347854845137454;
		xi[1][0] =  0.861136311594053;
		weight[1] = 0.347854845137454;
		xi[2][0] = -0.339981043584856;
		weight[2] = 0.652145154862546;
		xi[3][0] =  0.339981043584856;
		weight[3] = 0.652145154862546;
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}

