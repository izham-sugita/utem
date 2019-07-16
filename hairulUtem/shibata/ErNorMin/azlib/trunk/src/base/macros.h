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

/* $Id: macros.h 1043 2016-02-28 23:58:30Z hayashi $ */

#ifndef MACROS_H
#define MACROS_H

/* 容量 */
#define KILO_BYTE (1024)
#define MEGA_BYTE (1024*1024)
#define GIGA_BYTE (1024*1024*1024)

#ifdef USE_ISOC99
	/* C99 ‰Â•Ïˆø”ƒ}ƒNƒ */
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
