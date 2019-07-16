/*********************************************************************
 * rc.h
 *
 * Copyright (C) 2001 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Ryu SASAOKA> <Yasuo ASAGA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: rc.h,v 1.9 2003/12/05 06:26:25 sasaoka Exp $ */

#ifndef RC_H
#define RC_H

#include <stdio.h>
#include <stdlib.h>

typedef enum {
	NORMAL_RC,
	END_RC,
	IMPLEMENT_ERROR_RC,
	READ_ERROR_RC,
	WRITE_ERROR_RC, 
	CONVERT_ERROR_RC,
	NULL_ERROR_RC,
	ALLOC_ERROR_RC,
	OPEN_ERROR_RC,
	CLOSE_ERROR_RC,
	SEEK_ERROR_RC,
	FLUSH_ERROR_RC,
	ARG_ERROR_RC,
	CAL_ERROR_RC, 
	OVERFLOW_ERROR_RC,
	UNDERFLOW_ERROR_RC,
	ZERO_ERROR_RC,
	NEGATIVE_ERROR_RC,
	SUBMIT_ERROR_RC,
	UNKNOWN_ERROR_RC,
	SPECIAL_RC
} RC;


RC rc_malloc(size_t size, void **ptr);
RC rc_fopen(const char *path, const char *mode, FILE **fp);
RC rc_fclose(FILE *fp);

void rc_push_func(void (*func)(const char *, int, const char *, RC *));
void rc_pop_func(void);
void rc_error_print(const char *file_name, int line,
                    const char *rc_func, RC *rc);

void print_rc(FILE *stream, RC rc);

void rc_error_print_default(const char *file_name, int line,
                            const char *rc_func, RC *rc);
void rc_error_print_silent(const char *file_name, int line,
                           const char *rc_func, RC *rc);



#define RC_TRY_MAIN(func) {RC lrc = func;\
                           if(lrc != NORMAL_RC){\
                               rc_error_print(__FILE__, __LINE__, #func, &lrc);\
                               return(EXIT_FAILURE);\
                           }}

#define RC_TRY(func) {RC lrc = func;\
                      if(lrc != NORMAL_RC){\
                          rc_error_print(__FILE__, __LINE__, #func, &lrc);\
                          return(lrc);\
                      }}

#define RC_NEG_CHK(expr)  {if((expr) < 0){\
                              RC lrc = NEGATIVE_ERROR_RC;\
                              rc_error_print(__FILE__, __LINE__, #expr, &lrc);\
                               return(lrc);\
                           }}

#define RC_NEG_ZERO_CHK(expr)  {if((expr) <= 0){\
                                   RC lrc = NEGATIVE_ERROR_RC;\
                                   rc_error_print(__FILE__, __LINE__,\
                                                  #expr, &lrc);\
                                    return(lrc);\
                                }}

#define RC_NULL_CHK(expr)  {if((expr) == NULL){\
                               RC lrc = NULL_ERROR_RC;\
                               rc_error_print(__FILE__, __LINE__,\
                                              #expr, &lrc);\
                                return(lrc);\
                            }}

#define DEBUG_PRINT { fprintf(stderr, "[%s : %d]\n", __FILE__, __LINE__); }

#endif /* RC_H */
