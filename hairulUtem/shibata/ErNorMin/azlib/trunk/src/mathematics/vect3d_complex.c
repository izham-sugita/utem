/*********************************************************************
 * cvect3d.c
 *
 * Copyright (C) 2016 AzLib Developers Group
 *
 * Written by
 *  <Takuya HAYASHI>
 *
 *  Refer to log documents for details.
 *********************************************************************/

/* $Id: vect3d_complex.c 1094 2016-12-05 07:00:41Z hayashi $ */

#include <stdio.h>
#include "rc.h"
#include "complex_utl.h"

/* 複素3Dベクトル vect を複素共役にする．
 */
RC conj_cvect3d (CVECT3D *vect)
{
    vect->x = conj(vect->x);
    vect->y = conj(vect->y);
    vect->z = conj(vect->z);

    return(NORMAL_RC);
}

/* 複素3Dベクトル vect の各方向成分に 0.0 を代入する．
 */
RC init_cvect3d (CVECT3D *vect)
{
	vect->x = vect->y = vect->z = 0.0;

	return(NORMAL_RC);
}

/* 複素3DRベクトル vect の各方向成分に 0.0 を代入する．
 */
RC init_cvect3dr (CVECT3DR *vect)
{
	vect->x = vect->y = vect->z = 0.0;
	vect->yz = vect->zx = vect->xy = 0.0;

	return(NORMAL_RC);
}

/* 複素3Dベクトル vect を表示する．
 */
RC print_cvect3d (FILE *fp, CVECT3D vect)
{
    if(fp == NULL) return(ARG_ERROR_RC);

    fprintf(fp, "  x: %11.3e + %11.3e i\n", creal(vect.x), cimag(vect.x));
    fprintf(fp, "  y: %11.3e + %11.3e i\n", creal(vect.y), cimag(vect.y));
    fprintf(fp, "  z: %11.3e + %11.3e i\n", creal(vect.z), cimag(vect.z));

    return(NORMAL_RC);
}

/* 複素3DRベクトル vect を表示する．
 */
RC print_cvect3dr (FILE *fp, CVECT3DR vect)
{
    if(fp == NULL) return(ARG_ERROR_RC);

    fprintf(fp, "   x: %11.3e + %11.3e i\n", creal(vect.x), cimag(vect.x));
    fprintf(fp, "   y: %11.3e + %11.3e i\n", creal(vect.y), cimag(vect.y));
    fprintf(fp, "   z: %11.3e + %11.3e i\n", creal(vect.z), cimag(vect.z));
    fprintf(fp, "  yz: %11.3e + %11.3e i\n", creal(vect.yz), cimag(vect.yz));
    fprintf(fp, "  zx: %11.3e + %11.3e i\n", creal(vect.zx), cimag(vect.zx));
    fprintf(fp, "  xy: %11.3e + %11.3e i\n", creal(vect.xy), cimag(vect.xy));

    return(NORMAL_RC);
}


/* 複素3Dベクトル vect の絶対値を計算する． 
 */
double abs_cvect3d (CVECT3D vect)
{
	return(sqrt(sq_abs_cvect3d(vect)));
}


/* 複素3Dベクトル vect の絶対値2乗を計算する． 
 */
double sq_abs_cvect3d (CVECT3D vect)
{
	return((double)(conj(vect.x)*vect.x + conj(vect.y)*vect.y 
					+ conj(vect.z)*vect.z));
}

/* ret = w1*vect1 + w2*vect2 を計算する．
 */
RC wadd_cvect3d (COMPLEX w1, CVECT3D vect1, COMPLEX w2, CVECT3D vect2, 
		CVECT3D *ret)
{
	RC_NULL_CHK( ret );

	ret->x = w1*vect1.x + w2*vect2.x;
	ret->y = w1*vect1.y + w2*vect2.y;
	ret->z = w1*vect1.z + w2*vect2.z;

	return(NORMAL_RC);
}

COMPLEX inner_product_cvect3d (CVECT3D vect1, CVECT3D vect2)
{
	return(vect1.x*vect2.x + vect1.y*vect2.y + vect1.z*vect2.z);
		
}

COMPLEX inner_product_cvect3dr (CVECT3DR vect1, CVECT3DR vect2)
{
	return(vect1.x*vect2.x + vect1.y*vect2.y + vect1.z*vect2.z
			+ vect1.yz*vect2.yz + vect1.zx*vect2.zx + vect1.xy*vect2.xy);
}

RC tensor_product_cvect3d (CVECT3D vect1, CVECT3D vect2, COMPLEX **tensor)
{
    RC_NULL_CHK( tensor );

    tensor[0][0] = vect1.x*vect2.x;
    tensor[0][1] = vect1.x*vect2.y;
    tensor[0][2] = vect1.x*vect2.z;
    tensor[1][0] = vect1.y*vect2.x;
    tensor[1][1] = vect1.y*vect2.y;
    tensor[1][2] = vect1.y*vect2.z;
    tensor[2][0] = vect1.z*vect2.x;
    tensor[2][1] = vect1.z*vect2.y;
    tensor[2][2] = vect1.z*vect2.z;

    return(NORMAL_RC);
}


