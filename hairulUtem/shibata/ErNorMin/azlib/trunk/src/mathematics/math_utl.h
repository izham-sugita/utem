/*********************************************************************
 * math_utl.h
 *
 * Copyright (C) 2016 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Syoji ITO> <Ryu SASAOKA>
 *  <Takaaki NAGATANI> <Yasuo ASAGA>
 *  <Takuya HAYASHI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: math_utl.h 1094 2016-12-05 07:00:41Z hayashi $ */

#ifndef MATH_UTL_H
#define MATH_UTL_H

#include <stdio.h>
#include <float.h>
#include "rc.h"
#include "nonzero_cg.h"

#ifndef PI
#define PI    (3.14159265358979323846)
#endif /* PI */

#define ABS_TOL   (DBL_MIN / DBL_EPSILON) /* 絶対許容誤差 */
#define REL_TOL   (DBL_EPSILON * 1.0e+3)  /* 相対許容誤差 */
#define ABS_ERROR (DBL_MIN/DBL_EPSILON)   /* 絶対許容誤差(互換性のため) */
#define REL_ERROR (DBL_EPSILON * 1.0e+3)  /* 相対許容誤差(互換性のため) */
/* メモ: 絶対誤差 = 近似値 - 真値 ,  相対誤差 = 絶対誤差 / 真値 */
#define EXP_LIMIT    (0.99*log(DBL_MAX/2.0))

#define SIGN(a)    ((a) >= 0.0 ? 1.0 : -1.0)
#define ISIGN(a)   ((a) >= 0.0 ? 0 : 1)
#define MAX2(a, b) ((a) > (b) ? (a) : (b))
#define MIN2(a, b) ((a) < (b) ? (a) : (b))
#define MIN_MAX(min, a, max) (MIN2(MAX2((min), (a)), (max)))

#define KEEP_AWAY_ZERO(a) ((a) >= 0.0 ? ((a) + ABS_TOL) : ((a) - ABS_TOL))

#define SWAP(a, b, swp) ((swp) = (a), (a) = (b), (b) = (swp))

#define ERROR_SUM(sum, x, r, tmp)  \
        (r)+=(x), (tmp)=(sum)+(r), (r)-=(tmp)-(sum), (sum)=(tmp);

typedef struct {
	double x;
	double y;
	double z;
} VECT3D;

typedef struct {
	double x;
	double y;
	double z;
	double yz;
	double zx;
	double xy;
} VECT3DR;


/* math_utl.c */
double keep_away_zero(double v);
int nearly_eq(double v1, double v2);
RC sort_int_array(int n, int array[]);
RC sort_double_array(int n, double array[]);
RC uniq_int_array(int *n, int **array);
int search_int_array_bin(int n, const int array[], int i);
int search_int_array_lin(int n, const int array[], int i);
long double error_sum(long double sum, long double x, long double *r);
long double factorial(int i);
long double permutation(int n, int r);
long double combination(int n, int r);
RC bbqsort(void *base, size_t nmemb, size_t size,
           int (*compar)(const void *, const void *));
RC root_poly2(double a, double b, double c, double *x1, double *x2);
RC root_poly3(double a1, double a2, double a3,
              double *x1, double *x2, double *x3);
RC allocate_bit_flags(int n, unsigned int **flags);
RC free_bit_flags(int n, unsigned int **flags);
int chk_bit_flags(const unsigned int flags[], int i);
RC set_bit_flags(unsigned int flags[], int i);
RC unset_bit_flags(unsigned int flags[], int i);
RC allset_bit_flags(int n, unsigned int flags[]);
RC double2int_pos(double dv, int *iv);
RC chk_range_double(int size, double array[],
                    double *min, double *max, double *avg);
RC chk_range_int(int size, int array[], int *min, int *max, int *avg);
RC print_histogram_double(FILE *stream, int size, double array[]);
RC print_histogram_int(FILE *stream, int size, int array[]);
RC print_histogram_bar(FILE *stream, double hist_percent,
                       double percent, unsigned long num);


/* random.c */
void init_d_random(long idum);
double d_random(void);
long random_seed_time(void);
RC fill_random_array(int n, double array[], int seed);


/* matrix.c */
RC allocate1D(int size, double **array);
RC allocate2D(int record, int column, double ***array);
RC allocate3D(int record1, int record2, int record3, double ****array);
RC allocate1I(int size, int **array);
RC allocate2I(int record, int column, int ***array);
RC allocate3I(int record1, int record2, int record3, int ****array);
RC free1D(int size, double **array);
RC free2D(int record, int column, double ***array);
RC free3D(int record1, int record2, int record3, double ****array);
RC free1I(int size, int **array);
RC free2I(int record, int column, int ***array);
RC free3I(int record1, int record2, int record3, int ****array);
RC array2matrix(int record, int column, double array[], double *matrix[]);
void mul_scalar_matrix(int m, int n, double wa, double **A, double **B);
void mul_matrix(int l, int m, int n, double **A, double **B, double **C);
void mul_matrix_vect(int m, int n, double **A, const double *X, double *Y);
void mul_trans_matrix_vect(int m, int n, double **A, const double *X,double *Y);
void mul_vect_matrix(int m, int n, const double *X, double **A, double *Y);
double determinant(int dim, double **matrix);
RC inverse_matrix(int dim, double **A, double **inverseA);
RC add_matrix(int m, int n, double wa, double **A, double wb, double **B,
              double **C);
RC copy_vect(int size, double *vect1, double *vect2);
RC copy_matrix(int m, int n, double **A, double **B);
void transpose_matrix(int row, int col, double **A, double **B);
void transpose_overwrite_matrix(int n, double **A);
void init_vect(int n, double v[]);
void init_matrix(int row, int col, double **A);
RC unit_vect(int n, double *v);
void unit_matrix(int n, double **A);
RC tensor_product(int m, int n, const double *v1, const double *v2, double **A);
double inner_product(int n, const double *v1, const double *v2);
RC outer_product(const double *v1, const double *v2, double *ret);
long double inner_product_ac(int n, const long double *v1,
                                    const long double *v2);
double norm_vect(int size, double *vect);
void mul_matrix_AtB(int l, int m, int n, double **A, double **B, double **C);
void mul_matrix_AtBA(int m, int n, double **A, double **B, double **C);

void print_vect(FILE *fp, int n, const double v[]);
void print_matrix(FILE *fp, int m, int n, double **A);
RC static_condens_matrix(int n, double **A, double y[], int i);
RC row_col_add_zero3_matrix(int n, double **A, double y[], int i1, double w1,
                            int i2, double w2, int i3, double w3, int j);
RC row_col_add_zero_matrix(int n, double **A, double y[], int i_size,
                           const int i[], const double w[], int j);
RC symmetrize_matrix(int n, double **A);
RC flag_convert_matrix(int n, double **A, int flag_a[], int flag_b[],
                       double ***tmp_matrix);
RC eigen_jacobi_trans(int n, double **a, double *e_vals, double **e_vects);

/* lu.c */
RC lu_solve(double **A, int n, const double *b, double *x);
RC lu_decomp(double **A, int n, int *index, double *d);
RC lu_subst(double **A, int n, int *index, double *B);

/* qr.c */
RC qr_solve(double **A, int n, const double *b, double *x);
RC qr_factor(int n, double **A, double **Q, double **R);
RC qr_subst(int n, double **Q, double **R, const double *b, double *x);

/* lq.c */
RC lq_solve(double **A, int n, const double *b, double *x, int active[]);
RC lq_factor(int n, double **A, double **L, double **Q);
RC lq_subst(int n, double **L, double **Q, const double *b, double *x);
RC inverse_matrix_lq(int n, double **A, double **inverseA);


/* vect3d.c */
RC init_vect3d(VECT3D *v);
RC init_vect3dr(VECT3DR *v);
VECT3D arbitrary_ortho_vect3d(VECT3D a);
VECT3D gs_ortho3d(VECT3D a, VECT3D b);
VECT3D center3d_p4array(const double p1[], const double p2[],
                        const double p3[], const double p4[]);
VECT3D center3d_p3array(const double p1[], const double p2[],const double p3[]);
VECT3D mid_point3d(VECT3D p1, VECT3D p2);
double abs_vect3d(VECT3D v);
VECT3D unit_vect3d(VECT3D v);
double dist_point3d(VECT3D p1, VECT3D p2);
VECT3D extract_vect3d(VECT3DR v);
VECT3D add_vect3d(VECT3D v1, VECT3D v2);
VECT3D sub_vect3d(VECT3D v1, VECT3D v2);
VECT3DR add_vect3dr(VECT3DR v1, VECT3DR v2);
VECT3DR sub_vect3dr(VECT3DR v1, VECT3DR v2);
VECT3DR ave_vect3dr(VECT3DR v1, VECT3DR v2);
double scalar_triple_product3d(VECT3D v1, VECT3D v2, VECT3D v3);
double inner_product3d(VECT3D v1, VECT3D v2);
double inner_product3dr(VECT3DR v1, VECT3DR v2);
VECT3D outer_product3d(VECT3D v1, VECT3D v2);
double cos_vect3d(VECT3D v1, VECT3D v2);
double **tensor_product3d(VECT3D v1, VECT3D v2, double **ans);
RC print_vect3d(FILE *fp, VECT3D v);
RC print_vect3dr(FILE *fp, VECT3DR v);
RC print_matrix33(FILE *fp, double mat[][3]);
VECT3D wsum_vect3d(double fa, VECT3D A, double fb, VECT3D B);
VECT3DR wsum_vect3dr(double fa, VECT3DR A, double fb, VECT3DR B);
VECT3D mul_scalar_vect3d(double fa, VECT3D A);
VECT3DR mul_scalar_vect3dr(double fa, VECT3DR A);
VECT3D array2vect3d(double *array);

void init_matrix33(double matrix[][3]);
void unit_matrix33(double matrix[][3]);
VECT3D mul_matrix33_vect3d(double matrix[][3], VECT3D v);
VECT3D mul_matrix33t_vect3d(double matrix[][3], VECT3D v);
VECT3DR mul_matrix33_vect3dr(double matrix[][3], VECT3DR v);
VECT3DR mul_matrix33t_vect3dr(double matrix[][3], VECT3DR v);
void mul_matrix33_vect(double matrix[][3], double v[3]);
void mul_matrix33t_vect(double matrix[][3], double v[3]);
VECT3D array2vect3d(double *array);
VECT3DR array2vect3dr(double *array);
RC principal3d(VECT3DR ss_v, double pr_ss[],
               double pr_dir0[], double pr_dir1[], double pr_dir2[]);
RC principal3d_val(VECT3DR ss_v, double pr_ss[]);
RC mul_matrix33_L1AL2t(int n, double **A, double L1[][3], double L2[][3]);
RC mul_matrix33_LALt(int n, double **A, double L[][3]);
RC mul_matrix33_LtAL(int n, double **A, double L[][3]);
RC mul_matrix33_LtDL(int n, const double D[], double L[][3], double **A);
RC mul_matrix33(double A[][3], double B[][3], double C[][3]);
RC mul_matrix33_ABt(double A[][3], double B[][3], double C[][3]);
RC mul_matrix33_AtB(double A[][3], double B[][3], double C[][3]);
RC mul_add_matrix33(double A[][3], double B[][3], double C[][3]);
RC mul_add_matrix33_ABt(double A[][3], double B[][3], double C[][3]);
RC mul_add_matrix33_AtB(double A[][3], double B[][3], double C[][3]);
void add_matrix33(double A[][3], double B[][3]);
void add_matrix33t(double A[][3], double B[][3]);
RC add_matrix33a(double A[][3], double B[][3], double C[][3]);
RC mul_factor33(double factor, double A[][3], double B[][3]);
RC copy_matrix33(double A[][3], double B[][3]);
RC transpose_matrix33(double A[][3], double At[][3]);
RC allocate1D33(int record, double (**array)[3][3]);
RC allocate2D33(int record, int column, double (***array)[3][3]);
RC free1D33(int record, double (**array)[3][3]);
RC free2D33(int record, int column, double (***array)[3][3]);
RC print_0133(FILE *fp, int record, int column, double (**array)[3][3]);
int nonzero_chk33(double array[][3]);

/* geometry.c */
void plane_intersect_point(const VECT3D p[], VECT3D point,
                           VECT3D *term, double *distance, int *inside);
int is_inside_plane(const VECT3D p[], VECT3D point);
int is_inside_tetra(const VECT3D p[], VECT3D point);
double tri_area(VECT3D p1, VECT3D p2, VECT3D p3);        
double tetra_volume(VECT3D p1, VECT3D p2, VECT3D p3, VECT3D p4);
void line_intersect_point(VECT3D p0, VECT3D p1, VECT3D point,          
                          VECT3D *term, double *distance, int *inside);
int is_inside_line(VECT3D p0, VECT3D p1, VECT3D point);
int is_cross_aabb(VECT3D min1, VECT3D max1, VECT3D min2, VECT3D max2);

/* operation.c */
RC allocate_operation_matrix(double ***matrix);
void free_operation_matrix(double ***matrix);
RC init_operation_matrix(double **matrix);
RC print_operation_matrix (FILE *fp, double **matrix);
RC set_operation_matrix_pxpy(double fact, double **mat);
RC set_operation_matrix_mzpy(double fact, double **mat);
RC set_operation_matrix_pxmz(double fact, double **mat);
RC shift_operation_matrix(double dx, double dy, double dz, double **mat);
RC scale_operation_matrix(double x, double y, double z, double **mat);
RC rotation_operation_matrix_x(double rad, double **mat);
RC rotation_operation_matrix_y(double rad, double **mat);
RC rotation_operation_matrix_z(double rad, double **mat);
RC mul_operation_matrix(double **A, double **B);
VECT3D operation3d(double **matrix, VECT3D v);


#endif /* MATH_UTL_H */
