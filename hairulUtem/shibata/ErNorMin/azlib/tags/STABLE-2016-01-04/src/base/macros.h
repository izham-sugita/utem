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

/* $Id: macros.h 877 2007-02-16 12:05:47Z sasaoka $ */

#ifndef MACROS_H
#define MACROS_H

#define KILO_BYTE (1024)
#define MEGA_BYTE (1024*1024)
#define GIGA_BYTE (1024*1024*1024)

#define DEBUG_PRINT { log_printf(5, "[%8.8s:%d]\n", __FILE__, __LINE__); }

#ifdef USE_ISOC99
	/* C99 â¬ïœà¯êîÉ}ÉNÉç */
	#define DEBUG_LINE(format, ...) \
	 { log_printf(5, "[%8.8s:%d] " ## format, \
	              __FILE__, __LINE__, __VA_ARGS__); }
	#define DEBUG_TIME(format, ...) \
	 { tlog_printf(5, "[%8.8s:%d] " ## format, \
	              __FILE__, __LINE__, __VA_ARGS__); }
#else /* USE_ISOC99 */
	#define DEBUG_LINE(format, arg) \
	 { log_printf(5, "[%8.8s:%d] " ## format, __FILE__, __LINE__, arg); }
	#define DEBUG_TIME(format, arg) \
	 { tlog_printf(5, "[%8.8s:%d] " ## format, __FILE__, __LINE__, arg); }
#endif /* USE_ISOC99 */

#endif /* MACROS_H */
