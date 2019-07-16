/*********************************************************************
 * time_utl.c
 *
 * Copyright (C) 2004 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: time_utl.c 842 2006-12-04 06:21:20Z sasaoka $ */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "rc.h"

#define TIMER_STACK_SIZE 8
static time_t timer_start_stack[TIMER_STACK_SIZE];
static time_t timer_end_stack[TIMER_STACK_SIZE];


/* RFC 3339(Date and Time on the Internet: Timestamps) */
/* ex: 2004-05-10 02:49:05+09:00                       */
/* string size = 25                                    */
RC
time_string_stamp (char *str, size_t size)
{
	time_t tt;
	struct tm *tm_ptr;
	int ret;

	if(size < 26) return(ARG_ERROR_RC);

	tt = time(NULL);
	if(tt == (time_t)-1) return(UNKNOWN_ERROR_RC);

	tm_ptr = localtime(&tt);
	if(!tm_ptr) return(UNKNOWN_ERROR_RC);

	ret = strftime(str, size, "%Y-%m-%d %H:%M:%S+09:00", tm_ptr);
	if(ret != 25) RC_TRY( warning_printf(0, "Can't wrire full strings.\n") );
	str[25] = (char)'\0';

	return(NORMAL_RC);
}


/* RFC822:(Standard for the Format of ARPA Internet Text Messages) */
/* ex: Fri, 10 Apr 2004 09:49:39 +0900                             */
/* string size = 31                                                */
RC
time_string_mail (char *str, size_t size)
{
	time_t tt;
	struct tm *tm_ptr;
	int ret;

	if(size < 32) return(ARG_ERROR_RC);

	tt = time(NULL);
	if(tt == (time_t)-1) return(UNKNOWN_ERROR_RC);

	tm_ptr = localtime(&tt);
	if(!tm_ptr) return(UNKNOWN_ERROR_RC);

	ret = strftime(str, size, "%a, %d %b %Y %H:%M:%S +0900", tm_ptr);
	if(ret != 31) RC_TRY( warning_printf(0, "Can't wrire fullstrings.\n") );
	str[31] = (char)'\0';

	return(NORMAL_RC);
}


/* ISO 8601:1988 Date/Time Representations */
/* ex: 20040510T024905                     */
/* string size = 15                        */
RC
time_string_log (char *str, size_t size)
{
	time_t tt;
	struct tm *tm_ptr;
	int ret;

	if(size < 16) return(ARG_ERROR_RC);

	tt = time(NULL);
	if(tt == (time_t)-1) return(UNKNOWN_ERROR_RC);

	tm_ptr = localtime(&tt);
	if(!tm_ptr) return(UNKNOWN_ERROR_RC);

	ret = strftime(str, size, "%Y%m%dT%H%M%S", tm_ptr);
	if(ret != 15) RC_TRY( warning_printf(0, "Can't wrire full strings.\n") );
	str[15] = (char)'\0';

	return(NORMAL_RC);
}


/* ex: 20040510 */
RC
time_string_day (char *str, size_t size)
{
	time_t tt;
	struct tm *tm_ptr;
	int ret;

	if(size < 9) return(ARG_ERROR_RC);

	tt = time(NULL);
	if(tt == (time_t)-1) return(UNKNOWN_ERROR_RC);

	tm_ptr = localtime(&tt);
	if(!tm_ptr) return(UNKNOWN_ERROR_RC);

	ret = strftime(str, size, "%Y%m%d", tm_ptr);
	if(ret != 8) RC_TRY( warning_printf(0, "Can't wrire full strings.\n") );
	str[8] = (char)'\0';

	return(NORMAL_RC);
}


/* ex: 09:49:39 */
RC
time_string_times (char *str, size_t size)
{
	time_t tt;
	struct tm *tm_ptr;
	int ret;

	if(size < 9) return(ARG_ERROR_RC);

	tt = time(NULL);
	if(tt == (time_t)-1) return(UNKNOWN_ERROR_RC);

	tm_ptr = localtime(&tt);
	if(!tm_ptr) return(UNKNOWN_ERROR_RC);

	ret = strftime(str, size, "%H:%M:%S", tm_ptr);
	if(ret != 8) RC_TRY( warning_printf(0, "Can't wrire full strings.\n") );
	str[8] = (char)'\0';

	return(NORMAL_RC);
}


void
timer_start (int num)
{
	if(num < 0) num = 0;
	if(num > TIMER_STACK_SIZE-1) num = TIMER_STACK_SIZE-1;

	timer_start_stack[num] = time(NULL);
}


int
timer_end (int num)
{
	if(num < 0) num = 0;
	if(num > TIMER_STACK_SIZE-1) num = TIMER_STACK_SIZE-1;

	timer_end_stack[num] = time(NULL);
	if(timer_start_stack[num] == (time_t)-1) return(-1);
	if(timer_end_stack[num] == (time_t)-1) return(-1);
	
	return( (int)(timer_end_stack[num] - timer_start_stack[num]));
}


