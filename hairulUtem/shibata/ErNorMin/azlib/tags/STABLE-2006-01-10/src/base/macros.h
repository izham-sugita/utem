/*********************************************************************
 * macros.h
 *
 * Copyright (C) 2005 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: macros.h 422 2005-08-19 09:52:15Z sasaoka $ */

#ifndef MACROS_H
#define MACROS_H

#define KILO_BYTE (1024)
#define MEGA_BYTE (1024*1024)
#define GIGA_BYTE (1024*1024*1024)

#define DEBUG_PRINT { fprintf(stderr, "[%s : %d]\n", __FILE__, __LINE__); }

#endif /* MACROS_H */
