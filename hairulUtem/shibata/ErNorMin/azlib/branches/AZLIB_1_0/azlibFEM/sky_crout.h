/*********************************************************************
 * sky_crout.h
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: sky_crout.h,v 1.6 2003/07/23 04:06:32 sasaoka Exp $ */

#ifndef SKY_CROUT_H
#define SKY_CROUT_H

#include "rc.h"

RC s_crout_alloc(long n, const long *index1, long **index2,
                 double **array, long *array_size);
RC s_crout_free(long *index2, double *array);
double s_crout_get(long n, const long *index1,
                           const long *index2, const double *array,
                           int row, int col, RC *rc);
double *s_crout_ptr(long n, const long *index1,
                            const long *index2, const double *array,
                            int row, int col, RC *rc);
void s_crout_print(long n, const long *index1,
                           const long *index2, const double *array);
RC s_crout_mod(long n, long *index1, const long *index2,
                          double *array, double *vector, long m, double v);
RC s_crout_decomp(long n, const long *index1, const long *index2,
                                              double *array, int p_flag);
RC s_crout_subst(long n, const long *index1, const long *index2,
                         const double *array, double *vector, int p_flag);

#endif /* SKY_CROUT_H */


