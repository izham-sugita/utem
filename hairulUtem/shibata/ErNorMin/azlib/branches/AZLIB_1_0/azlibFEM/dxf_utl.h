/*********************************************************************
 * dxf_utl.h
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Ryu SASAOKA> <Takaaki NAGATANI>
 *
 * Thanks to following contributors.
 *  <Tomoyuki TSUBATA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: dxf_utl.h,v 1.7 2003/07/22 09:20:31 nagatani Exp $ */

#ifndef DXF_UTL_H
#define DXF_UTL_H

#include <stdio.h>
#include "rc.h"
#include "fem_struct.h"

/* basic color of DXF */
#define DXF_RED     1
#define DXF_YELLOW  2
#define DXF_GREEN   3
#define DXF_CYAN    4
#define DXF_BLUE    5
#define DXF_MAGENTA 6
#define DXF_WHITE   7

/* On DXF 256 color is defined with  hue + brightness + saturation */
/* hue of 256 color */
#define DXF_HUE_RED              10
#define DXF_HUE_REDNESSORANGE    20
#define DXF_HUE_ORANGE           30
#define DXF_HUE_FLAVESCENTORANGE 40
#define DXF_HUE_YELLOW           50
#define DXF_HUE_VERDANTGREEN     70
#define DXF_HUE_GREEN            90
#define DXF_HUE_INDIGOBLUE      110
#define DXF_HUE_SKYBLUE         130
#define DXF_HUE_VERDURE         150
#define DXF_HUE_BLUE            170
#define DXF_HUE_PURPLE          190
#define DXF_HUE_DARKCERISE      210
#define DXF_HUE_RUSSET          230
#define DXF_HUE_GRAY            250

/* brightness of 256 color */
#define DXF_BRIGHT_MOSTBRIGHT     0
#define DXF_BRIGHT_BRIGHT         2
#define DXF_BRIGHT_NORMAL         4
#define DXF_BRIGHT_DARK           6
#define DXF_BRIGHT_MOSTDARK       8

/* saturation of 256 color */
#define DXF_SAT_WITHOUT_HITE      0
#define DXF_SAT_WITH_WHITE        1

/* declaration of functions */
RC output_dxf_3dface(FILE *fp, char *def_layer, int def_color,
                     FEM_ELEMENT_ARRAY elem, FEM_NODE_ARRAY node);
RC output_dxf_3dface_lay(FILE *fp, char **layer, int num_layer,
                         FEM_ELEMENT_ARRAY elem, FEM_NODE_ARRAY node);
RC output_dxf_polyline(FILE *fp, FEM_ELEMENT_ARRAY element,
                       FEM_NODE_ARRAY node);
RC output_dxf_polyline_lay(FILE *fp, char *layer,
                            FEM_ELEMENT_ARRAY elem, FEM_NODE_ARRAY node);

RC input_dxf(FILE *fp, FEM_ELEMENT_ARRAY *element, FEM_NODE_ARRAY *node);
RC input_dxf_3dface(FILE *fp, FEM_ELEMENT_ARRAY *element, FEM_NODE_ARRAY *node);
RC input_dxf_3dface_lay(FILE *fp, char ***layer, int *num_layer,
                        FEM_ELEMENT_ARRAY *element, FEM_NODE_ARRAY *node);

RC input_dxf_polyline(FILE *fp, FEM_ELEMENT_ARRAY *element, FEM_NODE_ARRAY *node);
RC input_dxf_polyline_lay(FILE *fp, char ***layer, int *num_layer,
                          FEM_ELEMENT_ARRAY *element, FEM_NODE_ARRAY *node);

#endif /* DXF_UTL_H */


