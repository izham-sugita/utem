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

/* $Id: stl_utl.h 300 2004-06-09 11:00:54Z sasaoka $ */

#ifndef STL_UTL_H
#define STL_UTL_H

/* declaration of functions */
RC output_stl_3dface(FILE *fp, char *layer,
                     ELEMENT_ARRAY elem, NODE_ARRAY node);

RC input_stl_3dface(FILE *fp, ELEMENT_ARRAY *element, NODE_ARRAY *node);
RC input_stl_3dface_lay(FILE *fp, char **title, 
                        ELEMENT_ARRAY *element, NODE_ARRAY *node);

#endif /* STL_UTL_H */

