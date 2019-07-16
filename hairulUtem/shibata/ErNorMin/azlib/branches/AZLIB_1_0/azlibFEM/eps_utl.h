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

/* $Id: eps_utl.h,v 1.5 2003/10/07 11:41:32 sasaoka Exp $ */

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


RC output_eps(FILE *fp, EPS_PROP prop, FEM_ELEMENT_ARRAY elem,
              FEM_NODE_ARRAY node, double **op_matrix);


#endif /* EPS_UTL_H */

