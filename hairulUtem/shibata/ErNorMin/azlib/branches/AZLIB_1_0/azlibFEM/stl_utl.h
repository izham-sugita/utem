/*********************************************************************
 * stl_utl.h 
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: stl_utl.h,v 1.5 2003/07/23 04:06:32 sasaoka Exp $ */

#ifndef STL_UTL_H
#define STL_UTL_H

/* declaration of functions */
RC output_stl_3dface(FILE *fp, char *layer,
                     FEM_ELEMENT_ARRAY elem, FEM_NODE_ARRAY node);

RC input_stl_3dface(FILE *fp, FEM_ELEMENT_ARRAY *element, FEM_NODE_ARRAY *node);
RC input_stl_3dface_lay(FILE *fp, char **title, 
                        FEM_ELEMENT_ARRAY *element, FEM_NODE_ARRAY *node);

#endif /* STL_UTL_H */

