/*********************************************************************
 * math_utl.h
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Syoji ITO> <Ryu SASAOKA>
 *  <Takaaki NAGATANI> <Yasuo ASAGA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: math_utl.h,v 1.19 2003/12/02 10:20:47 sasaoka Exp $ */

#ifndef MATH_UTL_H
#define MATH_UTL_H

#include <stdio.h>
#include <float.h>
#include "rc.h"

#ifndef PI
#define PI    (3.14159265358979323846)
#endif /* PI */

#define ABS_TOL  (DBL_MIN / DBL_EPSILON)
#define REL_TOL  (DBL_EPSILON * 1.0e+3)

#define MAX_GAUSS_POINT (21)

#define SIGN(a)    ((a) >= 0.0 ? 1.0 : -1.0)
#define MAX2(a, b) ((a) > (b) ? (a) : (b))
#define MIN2(a, b) ((a) < (b) ? (a) : (b))
#define MIN_MAX(min, a, max) (MIN2(MAX2((min), (a)), (max)))


/* 乱数用パラメータ */
#define RANDOM_IM1 2147483563
#define RANDOM_IM2 2147483399
#define RANDOM_AM (1.0/RANDOM_IM1)       /* 乱数列の最小値 */
#define RANDOM_IMM1 (RANDOM_IM1-1)
#define RANDOM_IA1 40014
#define RANDOM_IA2 40692
#define RANDOM_IQ1 53668
#define RANDOM_IQ2 52774
#define RANDOM_IR1 12211
#define RANDOM_IR2 3791
#define RANDOM_NTAB 32
#define RANDOM_NDIV (1+RANDOM_IMM1/RANDOM_NTAB)
#define RANDOM_EPS 1.2e-7
#define RANDOM_RNMX (1.0-RANDOM_EPS)     /* 乱数列の最大値 */

typedef struct {
	double x;
	double y;
} TRANSLATION2D;


typedef struct {
	double x;
	double y;
	double z;
} TRANSLATION3D;


typedef struct {
	double yz;
	double zx;
	double xy;
} ROTATION3D;


typedef struct {
	double x;
	double y;
	double z;
	double yz;
	double zx;
	double xy;
} TRANS_ROTATION3D;


typedef struct {
	int ix;
	int iy;
	int iz;
	int iyz;
	int izx;
	int ixy;
} I_TRANS_ROTATION3D;


/* math_utl.h */
void my_random_init(long idum);
double my_random(void);
int nearly_eq(double v1, double v2);

RC free2D(int record, int column, double **array);
RC free3D(int record1, int record2, int record3, double ***array);
RC free4D(int record1, int record2,
          int record3, int record4, double ****array);
RC allocate1D(int size, double **array);
RC allocate2D(int record, int column, double ***array);
RC allocate3D(int record1, int record2, int record3, double ****array);
RC allocate4D(int record1, int record2,
              int record3, int record4, double *****array);
RC free2D_i(int record, int column, int **array);
RC free3D_i(int record1, int record2, int record3, int ***array);
RC free4D_i(int record1, int record2,
            int record3, int record4, int ****array);
RC allocate1D_i(int size, int **array);
RC allocate2D_i(int record, int column, int ***array);
RC allocate3D_i(int record1, int record2, int record3, int ****array);
RC allocate4D_i(int record1, int record2,
                int record3, int record4, int *****array);

void mul_matrix(int l, int m, int n, double **A, double **B, double **C);
void mul_matrix_vect(int m, int n, double **A, const double *X, double *Y);
void mul_trans_matrix_vect(int m, int n, double **A, const double *X,
                           double *Y);
void mul_vect_matrix(int m, int n, const double *X, double **A, double *Y);
double determinant(int dim, double **matrix);
double determinant2(double **matrix);
double determinant3(double **matrix);
RC inverse_matrix(int dim, double **A, double **inverseA);
RC add_matrix(int m, int n, double wa, double **A, double wb, double **B,
              double **C);
RC copy_matrix(int m, int n, double **A, double **B);
RC lu_solve(double **A, int n, const double *b, double *x);
RC lu_decomp(double **A, int n, int *index, double *d);
RC lu_subst(double **A, int n, int *index, double *B);
void transpose_matrix(int row, int col, double **A, double **B);
void transpose_overwrite_matrix(int n, double **A);
void init_matrix(int row, int col, double **A);
void unit_matrix(int n, double **A);

TRANSLATION2D init_translation2d(void);
TRANSLATION3D init_translation3d(void);
ROTATION3D init_rotation3d(void);
TRANS_ROTATION3D init_trans_rotation3d(void);
I_TRANS_ROTATION3D init_i_trans_rotation3d(void);

TRANSLATION2D mid_point2d(TRANSLATION2D p1, TRANSLATION2D p2);
TRANSLATION3D mid_point3d(TRANSLATION3D p1, TRANSLATION3D p2);
double abs_translation2d(TRANSLATION2D v);
double abs_translation3d(TRANSLATION3D v);
double abs_rotation3d(ROTATION3D v);
double dist_point2d(TRANSLATION2D p1, TRANSLATION2D p2);
double dist_point3d(TRANSLATION3D p1, TRANSLATION3D p2);
TRANSLATION2D sub_translation2d(TRANSLATION2D v1, TRANSLATION2D v2);
TRANSLATION3D sub_translation3d(TRANSLATION3D v1, TRANSLATION3D v2);
TRANSLATION2D add_translation2d(TRANSLATION2D v1, TRANSLATION2D v2);
TRANSLATION3D add_translation3d(TRANSLATION3D v1, TRANSLATION3D v2);
TRANS_ROTATION3D add_trans_rotation3d(TRANS_ROTATION3D v1, TRANS_ROTATION3D v2);
TRANS_ROTATION3D ave_trans_rotation3d(TRANS_ROTATION3D v1, TRANS_ROTATION3D v2);
double scalar_triple_product3d(TRANSLATION3D v1,
                               TRANSLATION3D v2, TRANSLATION3D v3);
double inner_product2d(TRANSLATION2D v1, TRANSLATION2D v2);
double inner_product3d(TRANSLATION3D v1, TRANSLATION3D v2);
double outer_product2d(TRANSLATION2D v1, TRANSLATION2D v2);
TRANSLATION3D outer_product3d(TRANSLATION3D v1, TRANSLATION3D v2);
double **tensor_product3d(TRANSLATION3D v1, TRANSLATION3D v2, double **ans);

RC allocate_operation_matrix(double ***matrix);
void free_operation_matrix(double **matrix);
RC init_operation_matrix(double **matrix);
RC set_operation_matrix_pxpy(double fact, double **mat);
RC set_operation_matrix_mzpy(double fact, double **mat);
RC set_operation_matrix_pxmz(double fact, double **mat);
RC shift_operation_matrix(double dx, double dy, double dz, double **mat);
RC scale_operation_matrix(double x, double y, double z, double **mat);
RC rotation_operation_matrix_x(double rad, double **mat);
RC rotation_operation_matrix_y(double rad, double **mat);
RC rotation_operation_matrix_z(double rad, double **mat);
RC mul_operation_matrix(double **A, double **B);
TRANSLATION3D operation3d(double **matrix, TRANSLATION3D v);

RC print_translation2d(FILE *fp, TRANSLATION2D v);
RC print_translation3d(FILE *fp, TRANSLATION3D v);
RC print_rotation3d(FILE *fp, ROTATION3D v);
RC print_trans_rotation3d(FILE *fp, TRANS_ROTATION3D v);
RC print_i_trans_rotation3d(FILE *fp, I_TRANS_ROTATION3D v);

RC tensor_product(int m, int n, const double *v1,
                                const double *v2, double **A);
double inner_product(int n, const double *v1, const double *v2);
void mul_matrix_AtBA(int m, int n, double **A, double **B, double **C);
void print_matrix(FILE *fp, int m, int n, double **A);
TRANSLATION2D mul_scalar_translation2d(double fa, TRANSLATION2D A);
TRANSLATION3D mul_scalar_translation3d(double fa, TRANSLATION3D A);
TRANS_ROTATION3D mul_scalar_trans_rotation3d(double fa, TRANS_ROTATION3D A);
TRANSLATION3D transrot2trans(TRANS_ROTATION3D A);
ROTATION3D transrot2rot(TRANS_ROTATION3D A);
TRANSLATION3D array2translation3d(double *array);

RC principal3d(TRANS_ROTATION3D ss_v, double pr_ss[],
               double pr_dir0[], double pr_dir1[], double pr_dir2[]);
RC eigen_jacobi_trans(int n, double **a, double *e_vals, double **e_vects);

RC double2int_pos(double dv, int *iv);
long double factorial(int i);
long double permutation(int n, int r);
long double combination(int n, int r);

RC qr_factor(int n, double **A, double **Q, double **R);
RC qr_subst(int n, double **Q, double **R, const double *b, double *x);
RC qr_solve(double **A, int n, const double *b, double *x);

void my_sort_int(int n, int array[]);


/* integration.c */
RC Gauss_const_line(int number_of_point,
                    double **xi, double *weight);
RC Gauss_const_tri(int number_of_point,
                   double **L1_L2_L3, double *weight);
RC Gauss_const_quad(int number_of_point,
                    double **xi_eta, double *weight);
RC Gauss_const_tetra(int number_of_point,
                     double **L1_L2_L3_L4, double *weight);
RC Gauss_const_penta(int number_of_point,
                     double **L1_L2_L3_zeta, double *weight);
RC Gauss_const_hexa(int number_of_point,
                    double **xi_eta_zeta, double *weight);


#endif /* MATH_UTL_H */

