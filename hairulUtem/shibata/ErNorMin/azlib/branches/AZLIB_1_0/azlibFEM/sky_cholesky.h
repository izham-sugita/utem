/*********************************************************************
 * sky_cholesky.h
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: sky_cholesky.h,v 1.6 2003/07/23 04:06:32 sasaoka Exp $ */

#ifndef SKY_CHOLESKY_H
#define SKY_CHOLESKY_H

#include "rc.h"

/* $BG[Ns$N3NJ](B */
/* index1[] $B$r;2>H$7!"(Bindex2[],array[] $B$r3NJ](B */
/* $B3NJ]$7$?G[Ns$OE,Ev$K=i4|2=$5$l$k(B */
RC s_cholesky_alloc(long n, const long *index1, long **index2,
                    double **array, long *array_size);

/* $BG[Ns$N2rJ|(B */
RC s_cholesky_free(long *index2, double *array);

/* m $BHVL\$N2r$,>o$K(B v $B$H$J$k$h$&!"(Barray[],vector[] $B$r=$@5(B */
RC s_cholesky_mod(long n, long *index1, const long *index2,
                  double *array, double *vector, long m, double v);
RC s_cholesky_mod_n(long n, long *index1, const long *index2,
                    double *array, double **vectors, long m, double v,
                    int vec_size);

/* $B2~D{%3%l%9%-!<K!$K$h$kJ,2r(B */
RC s_cholesky_decomp(long n, const long *index1, const long *index2,
                                         double *array, int p_flag);
RC s_cholesky_decomp3(long n, const long *index1, const long *index2,
                                               double *array, int p_flag);

/* $B2~D{%3%l%9%-!<K!$K$h$jF@$i$l$?2<;03Q9TNs!"BP3Q9TNs$rMQ$$$F(B */
/* $BA0?JBeF~!"8eB`BeF~$r9T$J$&(B */
RC s_cholesky_subst(long n, const long *index1, const long *index2,
                            const double *array, double *vector, int p_flag);

/* $BG[Ns$NFbMF$rI=<((B (debug$BMQ(B) */
void s_cholesky_print(long n, const long *index1,
                              const long *index2, const double *array);

double s_cholesky_get(long n, const long *index1, const long *index2, 
                      const double *array, int row, int col, RC *rc);

double *s_cholesky_ptr(long n, const long *index1, const long *index2, 
                       double *array, int row, int col, RC *rc);
                       
RC s_cholesky_mul_vect(long n, const long *index1, const long *index2,
                       const double *array, const double *vector,
                       double *answer);


#endif /* SKY_CHOLESKY_H */


