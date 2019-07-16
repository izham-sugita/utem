/*********************************************************************
 * sky_cholesky.h
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: sky_cholesky.h 415 2005-07-15 06:04:24Z sasaoka $ */

#ifndef SKY_CHOLESKY_H
#define SKY_CHOLESKY_H

#include "rc.h"

/* 配列の確保                                 */
/* index1[] を参照し、index2[],array[] を確保 */
/* 確保した配列は適当に初期化される           */
RC s_cholesky_alloc(long n, const long *index1, long **index2,
                    double **array, long *array_size);

/* 配列の解放 */
RC s_cholesky_free(long *index2, double *array);

/* m 番目の解が常に v となるよう、array[],vector[] を修正 */
RC s_cholesky_mod(long n, long *index1, const long *index2,
                  double *array, double *vector, long m, double v);
RC s_cholesky_mod_n(long n, long *index1, const long *index2,
                    double *array, double **vectors, long m, double v,
                    int vec_size);

/* 改訂コレスキー法による分解 */
RC s_cholesky_decomp(long n, const long *index1, const long *index2,
                     double *array, int p_flag);
RC s_cholesky_decomp3(long n, const long *index1, const long *index2,
                      double *array, int p_flag);

/* 改訂コレスキー法により得られた下三角行列、対角行列を用いて */
/* 前進代入、後退代入を行なう                                 */
RC s_cholesky_subst(long n, const long *index1, const long *index2,
                    const double *array, double *vector, int p_flag);

/* 配列の内容を表示 (debug用) */
void s_cholesky_print(FILE *stream, long n, const long *index1,
                      const long *index2, const double *array);

double s_cholesky_get(long n, const long *index1, const long *index2, 
                      const double *array, int row, int col, RC *rc);

double *s_cholesky_ptr(long n, const long *index1, const long *index2, 
                       double *array, int row, int col, RC *rc);

RC s_cholesky_mul_vect(long n, const long *index1, const long *index2,
                       const double *array, const double *vector,
                       double *answer);


#endif /* SKY_CHOLESKY_H */

