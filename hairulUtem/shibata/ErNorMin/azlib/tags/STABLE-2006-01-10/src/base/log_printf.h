/*********************************************************************
 * log_printf.h
 *
 * Copyright (C) 2005 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: log_printf.h 424 2005-08-21 12:14:11Z sasaoka $ */

#ifndef LOG_PRINTF_H
#define LOG_PRINTF_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "rc.h"


RC set_log_level(int level);
int get_log_level(void);
RC set_log_file(int i, FILE *fp);

RC log_printf(int level, const char *format, ...);
RC tlog_printf(int level, const char *format, ...);
RC progress_printf(int level, const char *format, ...);
RC progress_bar_printf(int level, double percent, const char *format, ...);
RC error_printf(int error_id, ...);
RC warning_printf(int warning_id, ...);

void set_log_printf_func(RC (*func)(FILE *, const char *, va_list));
void set_tlog_printf_func(RC (*func)(FILE *,const char *,va_list));
void set_progress_printf_func(RC (*func)(FILE *, const char *, va_list));
void set_progress_bar_printf_func(RC (*func)(FILE *, double,
                                             const char *, va_list));
void set_error_printf_func(RC (*func)(FILE *, int, va_list));
void set_warning_printf_func(RC (*func)(FILE *, int, va_list));

void reset_log_printf_func(void);
void reset_tlog_printf_func(void);
void reset_progress_printf_func(void);
void reset_progress_bar_printf_func(void);
void reset_error_printf_func(void);
void reset_warning_printf_func(void);

RC log_flush(void);

void rc_error_print_log(const char *file_name, int line,
                        const char *rc_func, RC rc);

#endif /* LOG_PRINTF_H */

