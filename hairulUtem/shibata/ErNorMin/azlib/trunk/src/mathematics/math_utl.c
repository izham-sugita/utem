/*********************************************************************
 * math_utl.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Syoji ITO> <Masanobu SUYAMA> <Ryu SASAOKA>
 *  <Takaaki NAGATANI> <Yasuo ASAGA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: math_utl.c 1075 2016-10-31 09:48:59Z hayashi $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "rc.h"
#include "math_utl.h"
#include "memory_manager.h"

#define ADM_REL_ERROR (2.0*DBL_EPSILON)
#define ADM_ABS_ERROR (10.0*DBL_MIN)
#define HISTOGRAM_NUM (20)
#define HISTOGRAM_HEIGHT (40)

/* ゼロ割回避 */
double
keep_away_zero (double v)
{
	if(v >= 0.0) return(v + ABS_TOL);

	return(v - ABS_TOL);
}


int
nearly_eq (double v1, double v2)
{
	double abs_v1_v2;
	double mean_v1_v2;
	double adm_error;

	abs_v1_v2 = fabs(v1 - v2);
	mean_v1_v2 = fabs(v1) + fabs(v2);
	adm_error = ABS_TOL + (0.5*REL_TOL*mean_v1_v2);

	if(abs_v1_v2 < adm_error) return(1);

	return(0);
}

RC
sort_int_array (int n, int array[])
{
	int ii1, ii2;
	int bb_num = 1;     /* bb_num >= 1 */
	int swap_flag = 0;
	int mval;
	int left, right;
	int swap;

	if(n == 0) return(NORMAL_RC);
	if(n < 0) return(ARG_ERROR_RC);
	if(array == NULL) return(ARG_ERROR_RC);

	for(ii1=0; ii1<bb_num; ii1++){
		/* 正方向 */
		ii2 = swap_flag + 1;
		swap_flag = 0;
		for(; ii2<(n-ii1); ii2++){
			if(array[ii2 - 1] > array[ii2]){
				swap = array[ii2 - 1];
				array[ii2 - 1] = array[ii2];
				array[ii2] = swap;
				swap_flag = ii2;
			}
		}
		if(swap_flag == 0) break;

		/* 逆方向 */
		ii2 = swap_flag - 1;
		swap_flag = 0;
		for(; ii2>ii1; ii2--){
			if(array[ii2 - 1] > array[ii2]){
				swap = array[ii2 - 1];
				array[ii2 - 1] = array[ii2];
				array[ii2] = swap;
				swap_flag = ii2;
			}
		}
		if(swap_flag == 0) break;
	}
	if(swap_flag == 0) return(NORMAL_RC);

	mval = (array[swap_flag - 1] + array[n - bb_num])/2;
	n -= (swap_flag + bb_num);
	array += swap_flag;

	left = 0;
	right = n - 1;
	while(1){
		while(array[left] < mval) left++;
		while(array[right] > mval) right--;
		if(left >= right) break;
		swap = array[left];
		array[left] = array[right];
		array[right] = swap;
		left++;
		right--;
	}
	if(left > 1) RC_TRY( sort_int_array(left, array) );
	if(n - (right + 1) > 1) RC_TRY( sort_int_array(n - (right + 1),
	                                               &(array[right + 1])) );

	return(NORMAL_RC);
}

RC
sort_double_array (int n, double array[])
{
	int ii1, ii2;
	int bb_num = 1;     /* bb_num >= 1 */
	int swap_flag = 0;
	double mval;
	int left, right;
	double swap;

	if(n == 0) return(NORMAL_RC);
	if(n < 0) return(ARG_ERROR_RC);
	if(array == NULL) return(ARG_ERROR_RC);

	for(ii1=0; ii1<bb_num; ii1++){
		/* 正方向 */
		ii2 = swap_flag + 1;
		swap_flag = 0;
		for(; ii2<(n-ii1); ii2++){
			if(array[ii2 - 1] > array[ii2]){
				swap = array[ii2 - 1];
				array[ii2 - 1] = array[ii2];
				array[ii2] = swap;
				swap_flag = ii2;
			}
		}
		if(swap_flag == 0) break;

		/* 逆方向 */
		ii2 = swap_flag - 1;
		swap_flag = 0;
		for(; ii2>ii1; ii2--){
			if(array[ii2 - 1] > array[ii2]){
				swap = array[ii2 - 1];
				array[ii2 - 1] = array[ii2];
				array[ii2] = swap;
				swap_flag = ii2;
			}
		}
		if(swap_flag == 0) break;
	}
	if(swap_flag == 0) return(NORMAL_RC);

	mval = (array[swap_flag - 1] + array[n - bb_num])/2;
	n -= (swap_flag + bb_num);
	array += swap_flag;

	left = 0;
	right = n - 1;
	while(1){
		while(array[left] < mval) left++;
		while(array[right] > mval) right--;
		if(left >= right) break;
		swap = array[left];
		array[left] = array[right];
		array[right] = swap;
		left++;
		right--;
	}
	if(left > 1) RC_TRY( sort_double_array(left, array) );
	if(n - (right + 1) > 1) RC_TRY( sort_double_array(n - (right + 1),
	                                               &(array[right + 1])) );

	return(NORMAL_RC);
}


/* 配列は既にソートされているものとする． */
/* 配列のメモリ領域の大きさは変更されない． */
RC
uniq_int_array (int *n, int **array)
{
	int ii1, ii2, ii3;
	int array_size = 0;

	if(*n < 0) return(ARG_ERROR_RC);
	if(*array == NULL) return(ARG_ERROR_RC);
	if(*n <= 1) return(NORMAL_RC);

	for(ii1=0,ii2=1; ii2<*n; ii1++){
		if( (*array)[ii1] == (*array)[ii2] ){
			while((*array)[ii1] == (*array)[ii2]){
				ii2++;
				if(ii2 >= *n){
					array_size = ii1+1; /* ii1より後ろは全て同じ値 */
					break;
				}
			}
			if(ii2 >= *n) break;
			/* ii2の値で ii1とii2の間を塗りつぶす */
			for(ii3=ii1+1; ii3<ii2; ii3++)
				(*array)[ii3] = (*array)[ii2];
		}else{
			array_size = ii1+2; /* 少なくともii1の隣までは値が異なる */
			ii2++;
		}
	}

	*n = array_size;

	return(NORMAL_RC);
}


/* バイナリーサーチ */
int
search_int_array_bin (int n, const int array[], int i)
{
	int found_index = -1;
	int low = 0;
	int high = n - 1;
	int middle;

	if((array == NULL)||(n <= 0)) return(-2);

	while(low <= high){
		middle = (low + high)/2;

		if(array[middle] == i){
			found_index = middle;
			break;
		}else if(array[middle] > i){
			high = middle - 1;
		}else{
			low = middle + 1;
		}
	}

	return(found_index);
}


/* 線形サーチ */
int
search_int_array_lin (int n, const int array[], int i)
{
	int ii1;
	int found_index = -1;

	for(ii1=0; ii1<n; ii1++){
		if(array[ii1] == i){
			found_index = ii1;
			break;
		}
	}

	return(found_index);
}


/* x + sum + r を返却  */
/* 桁落ち分を r に代入 */
long double
error_sum (long double sum, long double x, long double *r)
{
	long double ret;

	(*r) += x;
	ret = sum + (*r);
	(*r) -= (ret - sum);

	return(ret);
}


long double 
factorial (int i)
{
	int ii1;
	long double ret = 1.0L;

	if(i < 0) return(0.0L);
	for(ii1=1; ii1<=i; ii1++){
		ret *= (long double)ii1;
	}

	return(ret);
}


/* nPr */
long double
permutation (int n, int r)
{
	int ii1;
	long double ret = 1.0L;

	if( (n < 0) || (r < 0) || (r > n) ) return(0.0L);
	for(ii1=n-r+1; ii1<=n; ii1++){
		ret *= (long double)ii1;
	}

	return(ret);
}


/* nCr */
long double
combination (int n, int r)
{
	long double ret;

	if( (n < 0) || (r < 0) || (r > n) ) return(0.0L);
	if((n - r) < r) r = n - r;
	ret = permutation(n, r)/factorial(r);

	return(ret);
}


/* qsort()互換のソートプログラム                                      */
/* 大部分がソートされたデータを想定して，双方向バブルソートを数回行う */
/* それでも完全にソートできない場合にのみ qsort()を使う               */
RC
bbqsort (void *base, size_t nmemb, size_t size,
         int (*compar)(const void *, const void *))
{
	int ii1, ii2;
	int swap_flag;
	void *swap_area;
	int bb_num = 4 + (int)(log10((double)nmemb));

	if(nmemb <= 0) return(NORMAL_RC);
	RC_NULL_CHK( swap_area = mm_alloc(size) );

	swap_flag = 0;
	for(ii1=0; ii1<bb_num; ii1++){
		/* 正方向 */
		ii2 = swap_flag + 1;
		swap_flag = 0;
		for(; ii2<((int)nmemb-ii1); ii2++){
			void *elem1 = (void *)((char *)base + (ii2 - 1)*size);
			void *elem2 = (void *)((char *)base + ii2*size);

			if((*compar)(elem1, elem2) > 0){
				memcpy(swap_area, elem1, size);
				memcpy(elem1, elem2, size);
				memcpy(elem2, swap_area, size);
				swap_flag = ii2;
			}
		}
		if(swap_flag == 0) break;

		/* 逆方向 */
		ii2 = swap_flag - 1;
		swap_flag = 0;
		for(; ii2>ii1; ii2--){
			void *elem1 = (void *)((char *)base + (ii2 - 1)*size);
			void *elem2 = (void *)((char *)base + ii2*size);

			if((*compar)(elem1, elem2) > 0){
				memcpy(swap_area, elem1, size);
				memcpy(elem1, elem2, size);
				memcpy(elem2, swap_area, size);
				swap_flag = ii2;
			}
		}
		if(swap_flag == 0) break;
	}
	RC_TRY( mm_free(swap_area) );
	if(swap_flag == 0) return(NORMAL_RC);

	if(nmemb <= (size_t)(2*bb_num+1)) return(NORMAL_RC);
	nmemb -= 2*bb_num;
	qsort((void *)((char *)base + bb_num*size), nmemb, size, compar);

	return(NORMAL_RC);
}


/* ax^2 + bx + c = 0 の解 */
RC
root_poly2 (double a, double b, double c, double *x1, double *x2)
{
	double q = -0.5*(b + SIGN(b)*sqrt(b*b - 4*a*c));

	*x1 = q/a;
	*x2 = c/q;

	return(NORMAL_RC);
}


/* x^3 + a1 x^2 + a2 x + a3 = 0 の解 */
RC
root_poly3 (double a1, double a2, double a3,
            double *x1, double *x2, double *x3)
{
	double q = (a1*a1 - 3.0*a2)/9.0;
	double r = (2.0*a1*a1*a1 - 9.0*a1*a2 + 27.0*a3)/54.0;
	double q3 = q*q*q;

	if(q3 - r*r >= 0.0){
		double alpha = r/(ABS_TOL+sqrt(fabs(q3)));
		double m2sqrt_q = -2.0*sqrt(q);
		double theta;

		if(alpha >= (1.0 - 2.0*DBL_EPSILON)) alpha = 1.0 - 2.0*DBL_EPSILON;
		if(alpha <= -(1.0 - 2.0*DBL_EPSILON)) alpha = -(1.0 - 2.0*DBL_EPSILON);
		theta = acos(alpha);

		*x1 = m2sqrt_q*cos(theta/3.0) - a1/3.0;
		*x2 = m2sqrt_q*cos((theta + 2.0*PI)/3.0) - a1/3.0;
		*x3 = m2sqrt_q*cos((theta + 4.0*PI)/3.0) - a1/3.0;
	}else{
		double pow_v = pow(sqrt(fabs(r*r - q3)) + fabs(r), 1.0/3.0);

		*x1 = *x2 = *x3 = -SIGN(r)*(pow_v + q/(pow_v + ABS_TOL)) - a1/3.0;
	}
/*	fprintf(stderr, "%15.7e %15.7e %15.7e (%15.7e)\n",
	        (*x1)*(*x1)*(*x1) + a1*(*x1)*(*x1) + a2*(*x1) + a3,
	        (*x2)*(*x2)*(*x2) + a1*(*x2)*(*x2) + a2*(*x2) + a3,
	        (*x3)*(*x3)*(*x3) + a1*(*x3)*(*x3) + a2*(*x3) + a3,
	        fabs((*x1)*(*x1)*(*x1)) + fabs(a1*(*x1)*(*x1))
	        + fabs(a2*(*x1)) + fabs(a3));*/

	return(NORMAL_RC);
}


RC
allocate_bit_flags (int n, unsigned int **flags)
{
	int size;

	if(flags == NULL) return(ARG_ERROR_RC);
	if(n <= 0){
		*flags = NULL;
		return(NORMAL_RC);
	}

	size = (n + sizeof(unsigned int) - 1)/sizeof(unsigned int);
	RC_NULL_CHK( *flags = mm_alloc(size * sizeof(unsigned int)) );

	return(NORMAL_RC);
}


RC
free_bit_flags (int n, unsigned int **flags)
{
	if(flags == NULL) return(ARG_ERROR_RC);

	RC_TRY( mm_free(*flags) );
	*flags = NULL;

	return(NORMAL_RC);
}


int
chk_bit_flags (const unsigned int flags[], int i)
{
	int index = i/sizeof(unsigned int);
	unsigned int mask = ((unsigned int)1)<<(i%sizeof(unsigned int));

	if(i < 0) return(0);
	if(flags == NULL) return(0);

	return(flags[index] & mask);
}


RC
set_bit_flags (unsigned int flags[], int i)
{
	int index = i/sizeof(unsigned int);
	unsigned int mask = ((unsigned int)1)<<(i%sizeof(unsigned int));

	if(i < 0) return(ARG_ERROR_RC);
	if(flags == NULL) return(ARG_ERROR_RC);

	flags[index] |= mask;

	return(NORMAL_RC);
}


RC
unset_bit_flags (unsigned int flags[], int i)
{
	int index = i/sizeof(unsigned int);
	unsigned int mask = ((unsigned int)1)<<(i%sizeof(unsigned int));

	if(i < 0) return(ARG_ERROR_RC);
	if(flags == NULL) return(ARG_ERROR_RC);

	flags[index] &= ~mask;

	return(NORMAL_RC);
}


RC
allset_bit_flags (int n, unsigned int flags[])
{
	int ii1;
	int size;
	unsigned int bit_set = ~((unsigned int)0);

	if(flags == NULL) return(ARG_ERROR_RC);
	if(n <= 0) return(NORMAL_RC);

	size = (n + sizeof(unsigned int) - 1)/sizeof(unsigned int);
	for(ii1=0; ii1<size; ii1++){
		flags[ii1] = bit_set;
	}

	return(NORMAL_RC);
}


RC
double2int_pos (double dv, int *iv)
{
	if(dv < -ADM_ABS_ERROR) return(CAL_ERROR_RC);
	if(dv < ADM_ABS_ERROR){
		*iv = 0;
		return(NORMAL_RC);
	}

	*iv = (int)(dv + 0.4999);
	if( !nearly_eq((double)(*iv), dv) ) return(CAL_ERROR_RC);

	return(NORMAL_RC);
}


RC
chk_range_double (int size, double array[],
                  double *min, double *max, double *avg)
{
	int ii1;
	long double t_avg = 0.0;

	if(size < 1) return(ARG_ERROR_RC);
	*min = *max = *avg = array[0];
	for(ii1=0; ii1<size; ii1++){
		if(array[ii1] < *min) *min = array[ii1];
		if(array[ii1] > *max) *max = array[ii1];
		t_avg += array[ii1];
	}
	*avg = (double)(t_avg/size);
	return(NORMAL_RC);
}

RC
chk_range_int (int size, int array[], int *min, int *max, int *avg)
{
	int ii1;
	long int t_avg = 0L;

	if(size < 1) return(ARG_ERROR_RC);
	*min = *max = *avg = array[0];
	for(ii1=0; ii1<size; ii1++){
		if(array[ii1] < *min) *min = array[ii1];
		if(array[ii1] > *max) *max = array[ii1];
		t_avg += array[ii1];
	}
	*avg = (int)(t_avg/size);
	return(NORMAL_RC);
}


RC
print_histogram_double (FILE *stream, int size, double array[])
{
	int ii1;
	unsigned long hist[HISTOGRAM_NUM], hist_max = 0L;
	double min, max, avg;
	double band;

	if(size < 1) return(ARG_ERROR_RC); 

	for(ii1=0; ii1<HISTOGRAM_NUM; ii1++) hist[ii1] = 0;

	RC_TRY( chk_range_double(size, array, &min, &max, &avg) );
	band = (max - min) / HISTOGRAM_NUM;

	for(ii1=0; ii1<size; ii1++){
		int num = (int)( (array[ii1]-min) / band );
		if(num < 0) num = 0;
		if(num >= HISTOGRAM_NUM) num = HISTOGRAM_NUM-1;
		hist[num]++;
	}

	fprintf(stream, "\n * size: [%10d]\n"
	                " * min: [%11.4e]\n"
	                " * avg: [%11.4e]\n"
	                " * max: [%11.4e]\n\n"
	                " [Class(min-)]"
	                " [Relative Frequency (%%)]                      "
	                " [Frequency]\n", size, min, avg, max);
	for(ii1=0; ii1<HISTOGRAM_NUM; ii1++){
		hist_max = MAX2(hist[ii1], hist_max);
	}
	for(ii1=0; ii1<HISTOGRAM_NUM; ii1++){
		fprintf(stream, "  %11.4e", min + band * ii1);
		RC_TRY( print_histogram_bar(stream, (double)hist[ii1]/hist_max,
		                            (double)hist[ii1]/size, hist[ii1]) );
	}
	fprintf(stream, "\n");

	return(NORMAL_RC);
}


RC
print_histogram_int (FILE *stream, int size, int array[])
{
	int ii1;
	unsigned long hist[HISTOGRAM_NUM], hist_max = 0L;
	int min, max, avg;
	double band;

	if(size < 1) return(ARG_ERROR_RC); 

	for(ii1=0; ii1<HISTOGRAM_NUM; ii1++) hist[ii1] = 0;

	RC_TRY( chk_range_int(size, array, &min, &max, &avg) );
	band = (max - min) / (double)HISTOGRAM_NUM;

	for(ii1=0; ii1<size; ii1++){
		int num = (int)((array[ii1] - min) / band);
		if(num < 0) num = 0;
		if(num >= HISTOGRAM_NUM) num = HISTOGRAM_NUM-1;
		hist[num]++;
	}

	fprintf(stream, "\n * size: [%10d]\n"
	                " * min: [%10d]\n"
	                " * avg: [%10d]\n"
	                " * max: [%10d]\n\n"
	                " [Class(min-)]"
	                " [Relative Frequency (%%)]                      "
	                " [Frequency]\n", size, min, avg, max);
	for(ii1=0; ii1<HISTOGRAM_NUM; ii1++){
		hist_max = MAX2(hist[ii1], hist_max);
	}
	for(ii1=0; ii1<HISTOGRAM_NUM; ii1++){
		fprintf(stream, "   %10d", (int)(min + band * ii1));
		RC_TRY( print_histogram_bar(stream, (double)hist[ii1]/hist_max,
		                            (double)hist[ii1]/size, hist[ii1]) );
	}
	fprintf(stream, "\n");

	return(NORMAL_RC);
}


RC
print_histogram_bar (FILE *stream, double hist_percent,
                     double percent, unsigned long num)
{
	int ii1, str_num;
	char buf[80];

	if( hist_percent < 0.0) hist_percent = 0.0;
	if( hist_percent > 1.0) hist_percent = 1.0;
	if( percent < 0.0) percent = 0.0;
	if( percent > 1.0) percent = 1.0;

	str_num = sprintf(&buf[0], " |");
	if(str_num < 0) return(WRITE_ERROR_RC);

	for(ii1=str_num; ii1<(HISTOGRAM_HEIGHT+str_num); ii1++){
		if(ii1 < hist_percent*HISTOGRAM_HEIGHT){
			buf[ii1] = '=';
		}else{
			buf[ii1] = ' ';
		}
	}
	str_num = sprintf(&buf[ii1], "[%3d%%] %11ld\n", (int)(percent*100), num);
	if(str_num < 0) return(WRITE_ERROR_RC);

	buf[str_num+ii1+1] = '\0';

	if(fputs(buf, stream) == EOF) return(WRITE_ERROR_RC);
	
	return(NORMAL_RC);
}



