/*********************************************************************
 * geometry.c
 *
 * Copyright (C) 2004 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: geometry.c 842 2006-12-04 06:21:20Z sasaoka $ */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include "rc.h"
#include "math_utl.h"

static void plane_eq_factor(const VECT3D p[], VECT3D *nvect, double *d);
static void plane_distance(VECT3D point, VECT3D nvect, double d,
                           double *distance);
static void plane_intersect(VECT3D point, VECT3D nvect, double distance,
                            VECT3D *ipoint);

/* p[3]を通る平面に垂線を降ろした時，その交点について座標(term) */
/* 距離(distance)，交点が三角形の内部にあるか(inside)を調べる． */
/* 裏向きの面に対して交点を求めた場合，distanceは負になる       */
void
plane_intersect_point (const VECT3D p[], VECT3D point,
                       VECT3D *term, double *distance, int *inside)
{
	VECT3D nvect;
	double d, cos;

	plane_eq_factor(p, &nvect, &d);
	plane_distance(point, nvect, d, distance);

	cos = cos_vect3d(nvect, sub_vect3d(point, p[0]));
	if(cos < 0) *distance *= -1.0;

	plane_intersect(point, nvect, *distance, term);
	*inside = is_inside_plane(p, *term);
}


/* p[3]を通る平面方程式(ax + by + cz + d = 0)の係数を求める */
/* a,b,cは面の法線ベクトル(nvect)になる                     */
static void
plane_eq_factor (const VECT3D p[], VECT3D *nvect, double *d)
{
	VECT3D v1, v2;
	v1 = sub_vect3d( p[1], p[0] );
	v2 = sub_vect3d( p[2], p[0] );
	*nvect = outer_product3d(v1, v2);
	*d = p[0].x * (p[1].y * p[2].z - p[1].z * p[2].y)
	   + p[0].y * (p[1].z * p[2].x - p[1].x * p[2].z)
	   + p[0].z * (p[1].x * p[2].y - p[1].y * p[2].x);
}


/* point から平面への距離を計算                         */
/* nvect,dはplane_eq_factor()で得られた平面方程式の係数 */
/* distanceは pointから平面までの距離                   */
static void
plane_distance (VECT3D point, VECT3D nvect, double d, double *distance)
{
	*distance = fabs(-(  nvect.x * point.x
	                   + nvect.y * point.y
	                   + nvect.z * point.z) + d)
	          / keep_away_zero( abs_vect3d(nvect) );
}


/* point から平面への垂線をおろしたときの交点を計算    */
/* 得られた交点は目標とする面内にあるとは限らない．    */
/* nvectは面の法線，distanceは pointから平面までの距離 */
static void
plane_intersect (VECT3D point, VECT3D nvect, double distance, VECT3D *ipoint)
{
	VECT3D uvect = unit_vect3d(nvect);
	ipoint->x = point.x - distance * uvect.x;
	ipoint->y = point.y - distance * uvect.y;
	ipoint->z = point.z - distance * uvect.z;
}


/* p[3]で表される三角形の内に pointがあるか確認する  */
int
is_inside_plane (const VECT3D p[], VECT3D point)
{
	VECT3D center;
	double area1, area2;

	/* 三角形の中心から三点を含む円内にあるか確認する */
	center.x =  (p[0].x + p[1].x +  p[2].x)/3;
	center.y =  (p[0].y + p[1].y +  p[2].y)/3;
	center.z =  (p[0].z + p[1].z +  p[2].z)/3;
	if( (dist_point3d(center, p[0])+ABS_ERROR) < dist_point3d(center, point) ){
		return(0);
	}

	/* 点が重なっていたら平面上にある */
	if( dist_point3d(p[0], point) < ABS_ERROR) return(1);
	if( dist_point3d(p[1], point) < ABS_ERROR) return(1);
	if( dist_point3d(p[2], point) < ABS_ERROR) return(1);

	/* pointを頂点とする三角形を三つ計算し，その和が三角形の面積に */
	/* 等しいかを判定し，pointが面内にあるか判定する．ヘロンの公式 */
	/* を使うよりも外積計算の方が高速                              */
	area1 = tri_area(p[0], p[1], p[2]);
	area2 = tri_area(point, p[0], p[1])
	      + tri_area(point, p[1], p[2])
	      + tri_area(point, p[2], p[0]);

	if(nearly_eq(area1, area2)) return(1);

	return(0);
}


/* 三点を結ぶ三角形の面積計算 */
double
tri_area (VECT3D p1, VECT3D p2, VECT3D p3)
{
	return( 0.5 * abs_vect3d( outer_product3d(sub_vect3d(p1, p2),
	                                          sub_vect3d(p1, p3)) ) );
}


#if 0
/* ヘロンの公式 : 三点を結ぶ三角形の面積を計算 */
static double 
Heron_area (VECT3D p1, VECT3D p2, VECT3D p3)
{
	double a, b, c, s;
	a = abs_vect3d(sub_vect3d(p1, p2));
	b = abs_vect3d(sub_vect3d(p2, p3));
	c = abs_vect3d(sub_vect3d(p3, p1));
	s = 0.5*(a+b+c);
	return( sqrt(s*(s-a)*(s-b)*(s-c)) );
}
#endif


/* p0,p1を結ぶ直線とpointからの垂線の交点とその距離，得られた */
/* 交点が p0-p1間にあるかどうか判定する                       */
void 
line_intersect_point (VECT3D p0, VECT3D p1, VECT3D point,
                      VECT3D *term, double *distance, int *inside)
{
	VECT3D v_0_1 = sub_vect3d(p1, p0);
	VECT3D v_p_t = gs_ortho3d(sub_vect3d(p0, point), v_0_1);

	*distance = abs_vect3d(v_p_t);
	*term = add_vect3d(v_p_t, point);
	*inside = is_inside_line(p0, p1, *term);
}


/* p0, p1間に pointがあるか確認する  */
int
is_inside_line (VECT3D p0, VECT3D p1, VECT3D point)
{
	return( nearly_eq( dist_point3d(p0, point) + dist_point3d(p1, point),
	                   dist_point3d(p0, p1)) );
}

/* min-max(対角) で表されるBounding Boxの領域が交差しているかどうかの判定 */
/* aabb = 軸並行境界ボックス(axis-aligned bounding box) */
int
is_cross_aabb(VECT3D min1, VECT3D max1, VECT3D min2, VECT3D max2)
{
	/* 0 : 交差してない */
	if(MIN2(max1.x, max2.x) - MAX2(min1.x, min2.x) + ABS_TOL < 0.0f) return(0);
	if(MIN2(max1.y, max2.y) - MAX2(min1.y, min2.y) + ABS_TOL < 0.0f) return(0);
	if(MIN2(max1.z, max2.z) - MAX2(min1.z, min2.z) + ABS_TOL < 0.0f) return(0);
	return(1);
}


