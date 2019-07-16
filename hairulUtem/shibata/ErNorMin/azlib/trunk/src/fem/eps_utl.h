/*********************************************************************
 * eps_utl.h
 *
 * Copyright (C) 2001 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Ryu SASAOKA>
 *
 * Thanks to following contributors.
 *  <Tomoyuki TSUBATA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: eps_utl.h 842 2006-12-04 06:21:20Z sasaoka $ */

#ifndef EPS_UTL_H
#define EPS_UTL_H

#include <stdio.h>
#include "rc.h"
#include "fem_struct.h"

typedef struct {
	double line_rgb_color[3];
	double fill_rgb_color[3];
	double line_width;
} EPS_PROP;


RC output_eps(FILE *fp, EPS_PROP prop, ELEMENT_ARRAY elem,
              NODE_ARRAY node, double **op_matrix);


#endif /* EPS_UTL_H */

