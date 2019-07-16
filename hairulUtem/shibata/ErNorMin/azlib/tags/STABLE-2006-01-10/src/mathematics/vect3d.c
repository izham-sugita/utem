/*********************************************************************
 * vect3d.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Syoji ITO> <Masanobu SUYAMA> <Ryu SASAOKA>
 *  <Takaaki NAGATANI> <Yasuo ASAGA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: vect3d.c 719 2006-01-04 09:22:04Z sasaoka $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "rc.h"
#include "math_utl.h"
#include "memory_manager.h"
#include "log_printf.h"


RC
init_vect3d (VECT3D *v)
{
	v->x = v->y = v->z = 0.0;
	return(NORMAL_RC);
}


RC
init_vect3dr (VECT3DR *v)
{
	v->x = v->y = v->z = 0.0;
	v->yz = v->zx = v->xy = 0.0;
	return(NORMAL_RC);
}


/* a に直行するベクトルの一つを求める */
VECT3D
arbitrary_ortho_vect3d (VECT3D a)
{
	int ii1;
	int max;
	double a_v, max_a_v;
	VECT3D v[3] = { {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0} };

	a_v = abs_vect3d( outer_product3d(a, v[0]) );
	max = 0;
	max_a_v = a_v;
	for(ii1=1; ii1<3; ii1++){
		a_v = abs_vect3d( outer_product3d(a, v[ii1]) );
		if(a_v > max_a_v){
			max = ii1;
			max_a_v = a_v;
		}
	}

	return(unit_vect3d(outer_product3d(a, v[max])));
}


/* ベクトル a を b に直行するするように修正 */
VECT3D
gs_ortho3d (VECT3D a, VECT3D b)
{
	VECT3D unit_b = unit_vect3d(b);

	return( sub_vect3d(a,
	            mul_scalar_vect3d(inner_product3d(a, unit_b), unit_b)) );
}


VECT3D
center3d_p4array (const double p1[], const double p2[],
                  const double p3[], const double p4[])
{
	VECT3D ret;

	ret.x = (p1[0] + p2[0] + p3[0] + p4[0])/4.0;
	ret.y = (p1[1] + p2[1] + p3[1] + p4[1])/4.0;
	ret.z = (p1[2] + p2[2] + p3[2] + p4[2])/4.0;

	return(ret);
}


VECT3D
center3d_p3array (const double p1[], const double p2[], const double p3[])
{
	VECT3D ret;

	ret.x = (p1[0] + p2[0] + p3[0])/3.0;
	ret.y = (p1[1] + p2[1] + p3[1])/3.0;
	ret.z = (p1[2] + p2[2] + p3[2])/3.0;

	return(ret);
}


/* p1, p2 の中央に位置する点 */
VECT3D
mid_point3d (VECT3D p1, VECT3D p2)
{
	VECT3D ret;

	ret.x = (p1.x + p2.x)/2.0;
	ret.y = (p1.y + p2.y)/2.0;
	ret.z = (p1.z + p2.z)/2.0;

	return(ret);
}


/* ベクトルの絶対値 〜VECT3D ver.〜 */
double
abs_vect3d (VECT3D v)
{
	return(sqrt(v.x*v.x + v.y*v.y + v.z*v.z));
}


/* 単位ベクトル */
VECT3D
unit_vect3d (VECT3D v)
{
	double abs = abs_vect3d(v);
	VECT3D ret;

	if(abs > ABS_TOL){
		ret.x = v.x/abs;
		ret.y = v.y/abs;
		ret.z = v.z/abs;
	}else{
		ret.x = 1.0;
		ret.y = 0.0;
		ret.z = 0.0;
	}

	return(ret);
}


/* p1, p2 間の距離 */
double
dist_point3d (VECT3D p1, VECT3D p2)
{
	double ret;

	double dx = p1.x - p2.x;
	double dy = p1.y - p2.y;
	double dz = p1.z - p2.z;
	ret = sqrt(dx*dx + dy*dy + dz*dz);

	return(ret);
}

/* VECT3DR -> VECT3D */
VECT3D
extract_vect3d (VECT3DR v)
{
	VECT3D ret;

	ret.x = v.x;
	ret.y = v.y;
	ret.z = v.z;

	return(ret);
}


/* ベクトルの足し算 (v1 + v2) */
VECT3D
add_vect3d (VECT3D v1, VECT3D v2)
{
	VECT3D v;

	v.x = v2.x + v1.x;
	v.y = v2.y + v1.y;
	v.z = v2.z + v1.z;

	return(v);
}


/* ベクトルの引き算 (v1 - v2) */
VECT3D
sub_vect3d (VECT3D v1, VECT3D v2)
{
	VECT3D v;

	v.x = v1.x - v2.x;
	v.y = v1.y - v2.y;
	v.z = v1.z - v2.z;

	return(v);
}


/* ベクトルの足し算 (v1 + v2) */
VECT3DR
add_vect3dr (VECT3DR v1, VECT3DR v2)
{
	VECT3DR v;

	v.x = v2.x + v1.x;
	v.y = v2.y + v1.y;
	v.z = v2.z + v1.z;
	v.yz = v2.yz + v1.yz;
	v.zx = v2.zx + v1.zx;
	v.xy = v2.xy + v1.xy;

	return(v);
}


/* ベクトルの引き算 (v1 - v2) */
VECT3DR
sub_vect3dr (VECT3DR v1, VECT3DR v2)
{
	VECT3DR v;

	v.x = v1.x - v2.x;
	v.y = v1.y - v2.y;
	v.z = v1.z - v2.z;
	v.yz = v1.yz - v2.yz;
	v.zx = v1.zx - v2.zx;
	v.xy = v1.xy - v2.xy;

	return(v);
}


VECT3DR
ave_vect3dr (VECT3DR v1, VECT3DR v2)
{
	VECT3DR v;

	v.x = (v2.x + v1.x)/2.0;
	v.y = (v2.y + v1.y)/2.0;
	v.z = (v2.z + v1.z)/2.0;
	v.yz = (v2.yz + v1.yz)/2.0;
	v.zx = (v2.zx + v1.zx)/2.0;
	v.xy = (v2.xy + v1.xy)/2.0;

	return(v);
}


/* スカラー三重積 v1 (dot) v2 (cross) v3 */
double
scalar_triple_product3d (VECT3D v1, VECT3D v2, VECT3D v3)
{
	double s;
	VECT3D v;

	v = outer_product3d(v2, v3);
	s = inner_product3d(v1, v);

	return(s);
}


/* ベクトルの内積 v1 (dot) v2 */
double
inner_product3d (VECT3D v1, VECT3D v2)
{
	return(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}


/* ベクトルの内積 v1 (dot) v2 */
double
inner_product3dr (VECT3DR v1, VECT3DR v2)
{
	return(v1.x * v2.x   + v1.y * v2.y   + v1.z * v2.z
	     + v1.yz * v2.yz + v1.zx * v2.zx + v1.xy * v2.xy);
}


/* ベクトルの外積 v1 (cross) v2 */
VECT3D
outer_product3d (VECT3D v1, VECT3D v2)
{
	VECT3D v;

	v.x = (v1.y * v2.z) - (v1.z * v2.y);
	v.y = (v1.z * v2.x) - (v1.x * v2.z);
	v.z = (v1.x * v2.y) - (v1.y * v2.x);

	return(v);
}


/* cosθ = (a・b)/(|a||b|) */
double
cos_vect3d (VECT3D v1, VECT3D v2)
{
	return( inner_product3d(v1, v2)
	        / keep_away_zero( abs_vect3d(v1)*abs_vect3d(v2) ) );
}


/* ベクトルのテンソル積 v1 * v2 => ans[][] */
double **
tensor_product3d (VECT3D v1, VECT3D v2, double **ans)
{
	ans[0][0] = v1.x * v2.x;
	ans[0][1] = v1.x * v2.y;
	ans[0][2] = v1.x * v2.z;
	ans[1][0] = v1.y * v2.x;
	ans[1][1] = v1.y * v2.y;
	ans[1][2] = v1.y * v2.z;
	ans[2][0] = v1.z * v2.x;
	ans[2][1] = v1.z * v2.y;
	ans[2][2] = v1.z * v2.z;

	return(ans);
}


/* VECT3D の出力 */
RC
print_vect3d (FILE *fp, VECT3D v)
{
	fprintf(fp, "  x: %11.3e, y: %11.3e, z: %11.3e\n", v.x, v.y, v.z);

	return(NORMAL_RC);
}


/* VECT3DR の出力 */
RC
print_vect3dr (FILE *fp, VECT3DR v)
{
	fprintf(fp, "   x:%11.3e  y:%11.3e  z:%11.3e\n"
	            "  yz:%11.3e zx:%11.3e xy:%11.3e\n",
	             v.x, v.y, v.z, v.yz, v.zx, v.xy);

	return(NORMAL_RC);
}


/* 3x3のマトリックスの出力 */
RC
print_matrix33 (FILE *fp, double mat[][3])
{
	fprintf(fp, "%15.7e %15.7e %15.7e\n", mat[0][0], mat[0][1], mat[0][2]);
	fprintf(fp, "%15.7e %15.7e %15.7e\n", mat[1][0], mat[1][1], mat[1][2]);
	fprintf(fp, "%15.7e %15.7e %15.7e\n", mat[2][0], mat[2][1], mat[2][2]);

	return(NORMAL_RC);
}


/* fa * A + fb * B */
VECT3D
wsum_vect3d (double fa, VECT3D A, double fb, VECT3D B)
{
	VECT3D ret;

	ret.x = fa * A.x + fb * B.x;
	ret.y = fa * A.y + fb * B.y;
	ret.z = fa * A.z + fb * B.z;

	return(ret);
}


/* (fa * A + fb * B) */
VECT3DR
wsum_vect3dr (double fa, VECT3DR A, double fb, VECT3DR B)
{
	VECT3DR ret;

	ret.x = fa * A.x + fb * B.x;
	ret.y = fa * A.y + fb * B.y;
	ret.z = fa * A.z + fb * B.z;
	ret.yz = fa * A.yz + fb * B.yz;
	ret.zx = fa * A.zx + fb * B.zx;
	ret.xy = fa * A.xy + fb * B.xy;

	return(ret);
}


/* fa * A */
VECT3D
mul_scalar_vect3d (double fa, VECT3D A)
{
	VECT3D ret;

	ret.x = fa * A.x;
	ret.y = fa * A.y;
	ret.z = fa * A.z;

	return(ret);
}


/* fa * A */
VECT3DR
mul_scalar_vect3dr (double fa, VECT3DR A)
{
	VECT3DR ret;

	ret.x = fa * A.x;
	ret.y = fa * A.y;
	ret.z = fa * A.z;
	ret.yz = fa * A.yz;
	ret.zx = fa * A.zx;
	ret.xy = fa * A.xy;

	return(ret);
}


void
init_matrix33 (double matrix[][3])
{
	matrix[0][0] = matrix[0][1] = matrix[0][2] = 0.0;
	matrix[1][0] = matrix[1][1] = matrix[1][2] = 0.0;
	matrix[2][0] = matrix[2][1] = matrix[2][2] = 0.0;
}


void
unit_matrix33 (double matrix[][3])
{
	matrix[0][0] = 1.0;
	matrix[0][1] = 0.0;
	matrix[0][2] = 0.0;
	matrix[1][0] = 0.0;
	matrix[1][1] = 1.0;
	matrix[1][2] = 0.0;
	matrix[2][0] = 0.0;
	matrix[2][1] = 0.0;
	matrix[2][2] = 1.0;
}


/* matrix[3][3] と v の積 */
VECT3D
mul_matrix33_vect3d (double matrix[][3], VECT3D v)
{
	VECT3D ret;

	ret.x = matrix[0][0]*v.x + matrix[0][1]*v.y + matrix[0][2]*v.z;
	ret.y = matrix[1][0]*v.x + matrix[1][1]*v.y + matrix[1][2]*v.z;
	ret.z = matrix[2][0]*v.x + matrix[2][1]*v.y + matrix[2][2]*v.z;

	return(ret);
}


/* matrix[3][3] と v の積 */
VECT3DR
mul_matrix33_vect3dr (double matrix[][3], VECT3DR v)
{
	VECT3DR ret;

	ret.x = matrix[0][0]*v.x + matrix[0][1]*v.y + matrix[0][2]*v.z;
	ret.y = matrix[1][0]*v.x + matrix[1][1]*v.y + matrix[1][2]*v.z;
	ret.z = matrix[2][0]*v.x + matrix[2][1]*v.y + matrix[2][2]*v.z;

	ret.yz = matrix[0][0]*v.yz + matrix[0][1]*v.zx + matrix[0][2]*v.xy;
	ret.zx = matrix[1][0]*v.yz + matrix[1][1]*v.zx + matrix[1][2]*v.xy;
	ret.xy = matrix[2][0]*v.yz + matrix[2][1]*v.zx + matrix[2][2]*v.xy;

	return(ret);
}


/* matrix[3][3] と v の積 -> v */
void
mul_matrix33_vect (double matrix[][3], double v[3])
{
	VECT3D ret;

	ret.x = matrix[0][0]*v[0] + matrix[0][1]*v[1] + matrix[0][2]*v[2];
	ret.y = matrix[1][0]*v[0] + matrix[1][1]*v[1] + matrix[1][2]*v[2];
	ret.z = matrix[2][0]*v[0] + matrix[2][1]*v[1] + matrix[2][2]*v[2];

	v[0] = ret.x;
	v[1] = ret.y;
	v[2] = ret.z;
}


/* matrix[3][3]^T と v の積 */
VECT3D
mul_matrix33t_vect3d (double matrix[][3], VECT3D v)
{
	VECT3D ret;

	ret.x = matrix[0][0]*v.x + matrix[1][0]*v.y + matrix[2][0]*v.z;
	ret.y = matrix[0][1]*v.x + matrix[1][1]*v.y + matrix[2][1]*v.z;
	ret.z = matrix[0][2]*v.x + matrix[1][2]*v.y + matrix[2][2]*v.z;

	return(ret);
}


/* matrix[3][3]^T と v の積 */
VECT3DR
mul_matrix33t_vect3dr (double matrix[][3], VECT3DR v)
{
	VECT3DR ret;

	ret.x = matrix[0][0]*v.x + matrix[1][0]*v.y + matrix[2][0]*v.z;
	ret.y = matrix[0][1]*v.x + matrix[1][1]*v.y + matrix[2][1]*v.z;
	ret.z = matrix[0][2]*v.x + matrix[1][2]*v.y + matrix[2][2]*v.z;

	ret.yz = matrix[0][0]*v.yz + matrix[1][0]*v.zx + matrix[2][0]*v.xy;
	ret.zx = matrix[0][1]*v.yz + matrix[1][1]*v.zx + matrix[2][1]*v.xy;
	ret.xy = matrix[0][2]*v.yz + matrix[1][2]*v.zx + matrix[2][2]*v.xy;

	return(ret);
}


/* matrix[3][3]^T と v の積 -> v */
void
mul_matrix33t_vect (double matrix[][3], double v[3])
{
	VECT3D ret;

	ret.x = matrix[0][0]*v[0] + matrix[1][0]*v[1] + matrix[2][0]*v[2];
	ret.y = matrix[0][1]*v[0] + matrix[1][1]*v[1] + matrix[2][1]*v[2];
	ret.z = matrix[0][2]*v[0] + matrix[1][2]*v[1] + matrix[2][2]*v[2];

	v[0] = ret.x;
	v[1] = ret.y;
	v[2] = ret.z;
}


/* 配列を VECT3D に変換 */
VECT3D
array2vect3d (double *array)
{
	VECT3D ret;

	ret.x = array[0];
	ret.y = array[1];
	ret.z = array[2];

	return(ret);
}


/* 配列を VECT3DR に変換 */
VECT3DR
array2vect3dr (double *array)
{
	VECT3DR ret;

	ret.x = array[0];
	ret.y = array[1];
	ret.z = array[2];
	ret.yz = array[3];
	ret.zx = array[4];
	ret.xy = array[5];

	return(ret);
}


RC
principal3d (VECT3DR ss_v, double pr_ss[],
             double pr_dir0[], double pr_dir1[], double pr_dir2[])
{
	double *v_matrix[3];
	double *dir_matrix[3];
	double v_matrix_arr[9];
	double dir_matrix_arr[9];
	int ii1;

	v_matrix[0] = &(v_matrix_arr[0]);
	v_matrix[1] = &(v_matrix_arr[3]);
	v_matrix[2] = &(v_matrix_arr[6]);
	dir_matrix[0] = &(dir_matrix_arr[0]);
	dir_matrix[1] = &(dir_matrix_arr[3]);
	dir_matrix[2] = &(dir_matrix_arr[6]);

	v_matrix[0][0] = ss_v.x;
	v_matrix[1][1] = ss_v.y;
	v_matrix[2][2] = ss_v.z;
	v_matrix[0][1] = v_matrix[1][0] = ss_v.xy;
	v_matrix[1][2] = v_matrix[2][1] = ss_v.yz;
	v_matrix[2][0] = v_matrix[0][2] = ss_v.zx;

	/* 固有値 -> pr_ss                   */
	/* 固有ベクトル -> dir_matrix の各列 */
	RC_TRY( eigen_jacobi_trans(3, v_matrix, pr_ss, dir_matrix) );

	for(ii1=0; ii1<3; ii1++){
		if(pr_dir0 != NULL) pr_dir0[ii1] = dir_matrix[ii1][0];
		if(pr_dir1 != NULL) pr_dir1[ii1] = dir_matrix[ii1][1];
		if(pr_dir2 != NULL) pr_dir2[ii1] = dir_matrix[ii1][2];
	}

	return(NORMAL_RC);

}


RC
principal3d_val (VECT3DR ss_v, double pr_ss[])
{
	double xx = ss_v.x;
	double yy = ss_v.y;
	double zz = ss_v.z;
	double xy = ss_v.xy;
	double yz = ss_v.yz;
	double zx = ss_v.zx;
	double a1 = -(xx + yy + zz);
	double a2 = -(-xx*yy - xx*zz - yy*zz + yz*yz + zx*zx + xy*xy);
	double a3 = -(-xx*yz*yz - yy*zx*zx - zz*xy*xy + 2.0*xy*yz*zx + xx*yy*zz);

	RC_TRY( root_poly3(a1, a2, a3, &(pr_ss[0]), &(pr_ss[1]), &(pr_ss[2])) );

	return(NORMAL_RC);
}


/* L1[3][3] L2[3][3]を対角ブロックとする行列 L1m, L2m と A[n][n] の積 */
/* L1m*A*L2m^T => A                                                   */
RC
mul_matrix33_L1AL2t (int n, double **A, double L1[][3], double L2[][3])
{
	int ii1, ii2, ii3, ii4;
	double tmp[3][3];

	if(n%3 != 0) return(ARG_ERROR_RC);

	n /= 3;
	for(ii1=0; ii1<n; ii1++){
		for(ii2=0; ii2<n; ii2++){
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					tmp[ii3][ii4] = L1[ii3][0] * A[3*ii1  ][3*ii2+ii4]
					              + L1[ii3][1] * A[3*ii1+1][3*ii2+ii4]
					              + L1[ii3][2] * A[3*ii1+2][3*ii2+ii4];
				}
			}
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					A[3*ii1+ii3][3*ii2+ii4] = tmp[ii3][ii4];
				}
			}
		}
	}

	for(ii1=0; ii1<n; ii1++){
		for(ii2=0; ii2<n; ii2++){
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					tmp[ii3][ii4] = A[3*ii1+ii3][3*ii2  ] * L2[ii4][0]
					              + A[3*ii1+ii3][3*ii2+1] * L2[ii4][1]
					              + A[3*ii1+ii3][3*ii2+2] * L2[ii4][2];
				}
			}
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					A[3*ii1+ii3][3*ii2+ii4] = tmp[ii3][ii4];
				}
			}
		}
	}

	return(NORMAL_RC);
}


/* L[3][3] を対角ブロックとする行列 Lm と A[n][n] の積 */
/* Lm*A*Lm^T => A                                      */
RC
mul_matrix33_LALt (int n, double **A, double L[][3])
{
	int ii1, ii2, ii3, ii4;
	double tmp[3][3];

	if(n%3 != 0) return(ARG_ERROR_RC);

	n /= 3;
	for(ii1=0; ii1<n; ii1++){
		for(ii2=0; ii2<n; ii2++){
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					tmp[ii3][ii4] = L[ii3][0] * A[3*ii1  ][3*ii2+ii4]
					              + L[ii3][1] * A[3*ii1+1][3*ii2+ii4]
					              + L[ii3][2] * A[3*ii1+2][3*ii2+ii4];
				}
			}
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					A[3*ii1+ii3][3*ii2+ii4] = tmp[ii3][ii4];
				}
			}
		}
	}

	for(ii1=0; ii1<n; ii1++){
		for(ii2=0; ii2<n; ii2++){
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					tmp[ii3][ii4] = A[3*ii1+ii3][3*ii2  ] * L[ii4][0]
					              + A[3*ii1+ii3][3*ii2+1] * L[ii4][1]
					              + A[3*ii1+ii3][3*ii2+2] * L[ii4][2];
				}
			}
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					A[3*ii1+ii3][3*ii2+ii4] = tmp[ii3][ii4];
				}
			}
		}
	}

	return(NORMAL_RC);
}


/* L[3][3] を対角ブロックとする行列 Lm と A[n][n] の積 */
/* Lm^T*A*Lm => A                                      */
RC
mul_matrix33_LtAL (int n, double **A, double L[][3])
{
	int ii1, ii2, ii3, ii4;
	double tmp[3][3];

	if(n%3 != 0) return(ARG_ERROR_RC);

	n /= 3;
	for(ii1=0; ii1<n; ii1++){
		for(ii2=0; ii2<n; ii2++){
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					tmp[ii3][ii4] = L[0][ii3] * A[3*ii1  ][3*ii2+ii4]
					              + L[1][ii3] * A[3*ii1+1][3*ii2+ii4]
					              + L[2][ii3] * A[3*ii1+2][3*ii2+ii4];
				}
			}
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					A[3*ii1+ii3][3*ii2+ii4] = tmp[ii3][ii4];
				}
			}
		}
	}

	for(ii1=0; ii1<n; ii1++){
		for(ii2=0; ii2<n; ii2++){
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					tmp[ii3][ii4] = A[3*ii1+ii3][3*ii2  ] * L[0][ii4]
					              + A[3*ii1+ii3][3*ii2+1] * L[1][ii4]
					              + A[3*ii1+ii3][3*ii2+2] * L[2][ii4];
				}
			}
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					A[3*ii1+ii3][3*ii2+ii4] = tmp[ii3][ii4];
				}
			}
		}
	}

	return(NORMAL_RC);
}


/* L[3][3] を対角ブロックとする行列 L と D[n] を対角行列とする Dm の積 */
/* Lm^T*Dm*Lm => A                                                     */
RC
mul_matrix33_LtDL (int n, const double D[], double L[][3], double **A)
{
	int ii1, ii2, ii3;
	double tmp[3][3];

	if(n%3 != 0) return(ARG_ERROR_RC);

	n /= 3;
	for(ii1=0; ii1<n; ii1++){
		for(ii2=0; ii2<3; ii2++){
			for(ii3=0; ii3<3; ii3++){
				tmp[ii2][ii3] = L[ii3][ii2] * D[3*ii1 + ii3];
			}
		}

		for(ii2=0; ii2<3*n; ii2++){
			A[3*ii1  ][ii2] = 0.0;
			A[3*ii1+1][ii2] = 0.0;
			A[3*ii1+2][ii2] = 0.0;
		}

		for(ii2=0; ii2<3; ii2++){
			for(ii3=0; ii3<3; ii3++){
				A[3*ii1+ii2][3*ii1+ii3] = tmp[ii2][0] * L[0][ii3]
				                        + tmp[ii2][1] * L[1][ii3]
				                        + tmp[ii2][2] * L[2][ii3];
			}
		}
	}

	return(NORMAL_RC);
}


/* A[3][3] * B[3][3] => C[3][3] */
RC
mul_matrix33 (double A[][3], double B[][3], double C[][3])
{
	int ii1, ii2;

	for(ii1=0; ii1<3; ii1++){
		for(ii2=0; ii2<3; ii2++){
			C[ii1][ii2] = A[ii1][0]*B[0][ii2]
			            + A[ii1][1]*B[1][ii2]
			            + A[ii1][2]*B[2][ii2];
		}
	}

	return(NORMAL_RC);
}


/* A[3][3] * B^T[3][3] => C[3][3] */
RC
mul_matrix33_ABt (double A[][3], double B[][3], double C[][3])
{
	int ii1, ii2;

	for(ii1=0; ii1<3; ii1++){
		for(ii2=0; ii2<3; ii2++){
			C[ii1][ii2] = A[ii1][0]*B[ii2][0]
			            + A[ii1][1]*B[ii2][1]
			            + A[ii1][2]*B[ii2][2];
		}
	}

	return(NORMAL_RC);
}


/* A^T[3][3] * B[3][3] => C[3][3] */
RC
mul_matrix33_AtB (double A[][3], double B[][3], double C[][3])
{
	int ii1, ii2;

	for(ii1=0; ii1<3; ii1++){
		for(ii2=0; ii2<3; ii2++){
			C[ii1][ii2] = A[0][ii1]*B[0][ii2]
			            + A[1][ii1]*B[1][ii2]
			            + A[2][ii1]*B[2][ii2];
		}
	}

	return(NORMAL_RC);
}


/* A[3][3] * B[3][3] + C[3][3] => C[3][3] */
RC
mul_add_matrix33 (double A[][3], double B[][3], double C[][3])
{
	int ii1, ii2;

	for(ii1=0; ii1<3; ii1++){
		for(ii2=0; ii2<3; ii2++){
			C[ii1][ii2] += A[ii1][0]*B[0][ii2]
			             + A[ii1][1]*B[1][ii2]
			             + A[ii1][2]*B[2][ii2];
		}
	}

	return(NORMAL_RC);
}


/* A[3][3] * B^T[3][3] + C[3][3] => C[3][3] */
RC
mul_add_matrix33_ABt (double A[][3], double B[][3], double C[][3])
{
	int ii1, ii2;

	for(ii1=0; ii1<3; ii1++){
		for(ii2=0; ii2<3; ii2++){
			C[ii1][ii2] += A[ii1][0]*B[ii2][0]
			             + A[ii1][1]*B[ii2][1]
			             + A[ii1][2]*B[ii2][2];
		}
	}

	return(NORMAL_RC);
}


/* A^T[3][3] * B[3][3] => C[3][3] */
RC
mul_add_matrix33_AtB (double A[][3], double B[][3], double C[][3])
{
	int ii1, ii2;

	for(ii1=0; ii1<3; ii1++){
		for(ii2=0; ii2<3; ii2++){
			C[ii1][ii2] += A[0][ii1]*B[0][ii2]
			             + A[1][ii1]*B[1][ii2]
			             + A[2][ii1]*B[2][ii2];
		}
	}

	return(NORMAL_RC);
}


/* A = A + B */
void add_matrix33(double A[][3], double B[][3])
{
	A[0][0] += B[0][0];
	A[0][1] += B[0][1];
	A[0][2] += B[0][2];
	A[1][0] += B[1][0];
	A[1][1] += B[1][1];
	A[1][2] += B[1][2];
	A[2][0] += B[2][0];
	A[2][1] += B[2][1];
	A[2][2] += B[2][2];
}


/* A = A + B^t */
void add_matrix33t(double A[][3], double B[][3])
{
	A[0][0] += B[0][0];
	A[0][1] += B[1][0];
	A[0][2] += B[2][0];
	A[1][0] += B[0][1];
	A[1][1] += B[1][1];
	A[1][2] += B[2][1];
	A[2][0] += B[0][2];
	A[2][1] += B[1][2];
	A[2][2] += B[2][2];
}


/* A[3][3] + B[3][3] => C[3][3] */
RC
add_matrix33a (double A[][3], double B[][3], double C[][3])
{
	C[0][0] = A[0][0] + B[0][0];
	C[0][1] = A[0][1] + B[0][1];
	C[0][2] = A[0][2] + B[0][2];
	C[1][0] = A[1][0] + B[1][0];
	C[1][1] = A[1][1] + B[1][1];
	C[1][2] = A[1][2] + B[1][2];
	C[2][0] = A[2][0] + B[2][0];
	C[2][1] = A[2][1] + B[2][1];
	C[2][2] = A[2][2] + B[2][2];

	return(NORMAL_RC);
}


/* factor * A[3][3] => B[3][3] */
RC
mul_factor33 (double factor, double A[][3], double B[][3])
{
	B[0][0] = factor * A[0][0];
	B[0][1] = factor * A[0][1];
	B[0][2] = factor * A[0][2];
	B[1][0] = factor * A[1][0];
	B[1][1] = factor * A[1][1];
	B[1][2] = factor * A[1][2];
	B[2][0] = factor * A[2][0];
	B[2][1] = factor * A[2][1];
	B[2][2] = factor * A[2][2];

	return(NORMAL_RC);
}


/* A[3][3] => B[3][3] */
RC
copy_matrix33 (double A[][3], double B[][3])
{
	B[0][0] = A[0][0];
	B[0][1] = A[0][1];
	B[0][2] = A[0][2];
	B[1][0] = A[1][0];
	B[1][1] = A[1][1];
	B[1][2] = A[1][2];
	B[2][0] = A[2][0];
	B[2][1] = A[2][1];
	B[2][2] = A[2][2];

	return(NORMAL_RC);
}


/* A^t[3][3] => At[3][3] */
RC
transpose_matrix33 (double A[][3], double At[][3])
{
	At[0][0] = A[0][0];
	At[1][0] = A[0][1];
	At[2][0] = A[0][2];
	At[0][1] = A[1][0];
	At[1][1] = A[1][1];
	At[2][1] = A[1][2];
	At[0][2] = A[2][0];
	At[1][2] = A[2][1];
	At[2][2] = A[2][2];

    return(NORMAL_RC);
}


RC
allocate1D33 (int record, double (**array)[3][3])
{
	int ii1, ii2, ii3;

	if(record <= 0) return(ARG_ERROR_RC);

	*array = (double (*)[3][3])mm_alloc(record * sizeof(double [3][3]));
	RC_NULL_CHK( *array );

	for(ii1=0; ii1<record; ii1++){
		for(ii2=0; ii2<3; ii2++){
			for(ii3=0; ii3<3; ii3++){
				(*array)[ii1][ii2][ii3] = 0.0;
			}
		}
	}

	return(NORMAL_RC);
}


RC 
allocate2D33 (int record, int column, double (***array)[3][3])
{
	int ii1, ii2, ii3, ii4;
	double (*ptr)[3][3];

	if((record <= 0)||(column <= 0)) return(ARG_ERROR_RC);

	*array = (double (**)[3][3])mm_alloc(record * sizeof(double (*)[3][3]));
	RC_NULL_CHK( *array );

	ptr = (double (*)[3][3])mm_alloc(record*column * sizeof(double [3][3]));
	RC_NULL_CHK( ptr );

	for(ii1=0; ii1<record; ii1++){
		(*array)[ii1] = &(ptr[ii1*column]);
	}

	for(ii1=0; ii1<record; ii1++){
		for(ii2=0; ii2<column; ii2++){
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					(*array)[ii1][ii2][ii3][ii4] = 0.0;
				}
			}
		}
	}

	return(NORMAL_RC);
}


RC
free1D33 (int record, double (**array)[3][3])
{
	if(record <= 0) return(ARG_ERROR_RC);

	RC_TRY( mm_free(*array) );
	*array = NULL;

	return(NORMAL_RC);
}


RC
free2D33 (int record, int column, double (***array)[3][3])
{
	if((record <= 0)||(column <= 0)) return(ARG_ERROR_RC);

	RC_TRY( mm_free((*array)[0]) );
	RC_TRY( mm_free(*array) );
	*array = NULL;

	return(NORMAL_RC);
}


int
nonzero_chk33 (double array[][3])
{
	int ii1, ii2;

	for(ii1=0; ii1<3; ii1++){
		for(ii2=0; ii2<3; ii2++){
			if(nearly_eq(array[ii1][ii2], 0.0) == 0) return(1);
		}
	}

	return(0);
}


RC
print_0133 (FILE *fp, int record, int column, double (**array)[3][3])
{
	int ii1, ii2;

	for(ii1=0; ii1<record; ii1++){
		for(ii2=0; ii2<column; ii2++){
			if(nonzero_chk33(array[ii1][ii2])){
				if(ii1 == ii2){
					fprintf(fp, "@ ");
				}else{
					fprintf(fp, "1 ");
				}
			}else{
				fprintf(fp, "0 ");
			}
		}
		fprintf(fp, "\n");
	}

	return(NORMAL_RC);
}

