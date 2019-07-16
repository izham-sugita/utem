/*********************************************************************
 * rc.h
 *
 * Copyright (C) 2001 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Ryu SASAOKA> <Yasuo ASAGA>
 *  <Takuya HAYASHI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: rc.h 1045 2016-03-09 07:00:14Z hayashi $ */

#ifndef RC_H
#define RC_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>


typedef enum {
	NORMAL_RC = 0,
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
	OVERRUN_ERROR_RC,
	IGNORE_ERROR_RC,
	UNKNOWN_ERROR_RC,
	SPECIAL_RC
} RC;


/* rc.c */
RC rc_fopen(const char *path, const char *mode, FILE **fp);
RC rc_fclose(FILE *fp);
RC rc_touch(const char *path);

void rc_push_func(void (*func)(const char *, int, const char *, RC));
void rc_pop_func(void);
void rc_error_print(const char *file_name, int line,
                    const char *rc_func, RC rc);

void print_rc(FILE *stream, RC rc);

void rc_error_print_default(const char *file_name, int line,
                            const char *rc_func, RC rc);
void rc_error_print_silent(const char *file_name, int line,
                           const char *rc_func, RC rc);

RC rc_fwrite(const void *ptr, size_t size, size_t nmemb, FILE *fp);
RC rc_fwrite_fn(const void *ptr, size_t size, size_t nmemb, FILE *fp,
                const char file_name[]);
RC rc_fread(void *ptr, size_t size, size_t nmemb, FILE *fp);
RC rc_fread_fn(void *ptr, size_t size, size_t nmemb, FILE *fp,
               const char file_name[]);

/* log_printf.c */
RC set_log_level(int level);
int get_log_level(void);
RC set_log_file(int i, FILE *fp);
FILE *get_log_file(int i);
RC tlog_printf(int level, const char *format, ...);
RC log_printf(int level, const char *format, ...);
RC log_printf_line(int level, char line_char, const char *format, ...);
RC error_printf(int error_id, ...);
RC warning_printf(int warning_id, ...);
RC log_flush(void);

RC progress_printf(int level, const char *format, ...);
RC progress_bar_printf(int level, double percent, const char *format, ...);

void rc_error_print_log(const char *file_name, int line,
                        const char *rc_func, RC rc);

#define AZLIB_VERSION \
{\
	RC_TRY( log_printf(1, "**************************************************\n") );\
	RC_TRY( log_printf(1, "*                                                *\n") );\
	RC_TRY( log_printf(1, "*                                                *\n") );\
	RC_TRY( log_printf(1, "*                                                *\n") );\
	RC_TRY( log_printf(1, "*                                                *\n") );\
	RC_TRY( log_printf(1, "*                                                *\n") );\
	RC_TRY( log_printf(1, "*                                                *\n") );\
	RC_TRY( log_printf(1, "*                                                *\n") );\
	RC_TRY( log_printf(1, "**************************************************\n") );\
}


#define RC_TRY_MAIN(func) \
{\
    RC lrc = func;\
    if(lrc != NORMAL_RC){\
        rc_error_print(__FILE__, __LINE__, #func, lrc);\
        return(EXIT_FAILURE);\
    }\
}

#define RC_TRY(func) \
{\
    RC lrc = func;\
    if(lrc != NORMAL_RC){\
        rc_error_print(__FILE__, __LINE__, #func, lrc);\
        return(lrc);\
    }\
}

#define RC_TRY1(func, rc1) \
{\
    RC lrc = func;\
    if((lrc != NORMAL_RC)&&(lrc != rc1)){\
        rc_error_print(__FILE__, __LINE__, #func, lrc);\
        return(lrc);\
    }\
}

#define RC_TRY2(func, rc1, rc2) \
{\
    RC lrc = func;\
    if((lrc != NORMAL_RC)&&(lrc != rc1)&&(lrc != rc2)){\
         rc_error_print(__FILE__, __LINE__, #func, lrc);\
         return(lrc);\
    }\
}

#define RC_NEG_CHK(expr) \
{\
    if((expr) < 0){\
        RC lrc = NEGATIVE_ERROR_RC;\
        rc_error_print(__FILE__, __LINE__, #expr, lrc);\
        return(lrc);\
    }\
}

#define RC_NEG_ZERO_CHK(expr) \
{\
    if((expr) <= 0){\
        RC lrc = NEGATIVE_ERROR_RC;\
        rc_error_print(__FILE__, __LINE__, #expr, lrc);\
        return(lrc);\
    }\
}

#define RC_ZERO_CHK(expr) \
{\
    if((expr) == 0){\
        RC lrc = ZERO_ERROR_RC;\
        rc_error_print(__FILE__, __LINE__, #expr, lrc);\
        return(lrc);\
    }\
}

#define RC_NULL_CHK(expr) \
{\
    if((expr) == NULL){\
        RC lrc = NULL_ERROR_RC;\
        rc_error_print(__FILE__, __LINE__, #expr, lrc);\
        return(lrc);\
    }\
}

#define RC_NAN_CHK(expr) \
{\
	if(isnan((double)expr)){\
		RC lrc = CAL_ERROR_RC;\
		rc_error_print(__FILE__, __LINE__, #expr, lrc);\
		return(lrc);\
	}\
}

#define RC_INDEX_CHK(ARRAY_SIZE, ARRAY_INDEX) \
{\
	if(ARRAY_SIZE < 0){\
		return(ARG_ERROR_RC);\
	}\
	if(ARRAY_INDEX >= ARRAY_SIZE){\
		fprintf(stderr, "[%s : %d] <<< OVERRUN IN THIS LOOP!!\n",\
		        __FILE__, __LINE__);\
		return(OVERRUN_ERROR_RC);\
	}\
	if(ARRAY_INDEX < 0){\
		fprintf(stderr, "[%s : %d] <<< NEGATIVE INDEX IN THIS LOOP!!\n",\
		        __FILE__, __LINE__);\
		return(NEGATIVE_ERROR_RC);\
	}\
}

#define DEBUG_PRINT { log_printf(1, "[%8.8s:%d]\n", __FILE__, __LINE__); }


#endif /* RC_H */
