/*********************************************************************
 * time_utl.h
 *
 * Copyright (C) 2004 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA>
 *  <Takuya HAYASHI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: time_utl.h 1034 2016-01-13 07:43:41Z hayashi $ */

#ifndef TIME_UTL_H
#define TIME_UTL_H

RC time_string_stamp(char *str, size_t size);
RC time_string_mail(char *str, size_t size);
RC time_string_log(char *str, size_t size);
RC time_string_day(char *str, size_t size);
RC time_string_times(char *str, size_t size);

void timer_start(int num);
int timer_end(int num);
void print_timer(FILE *fp, int time);

#endif /* TIME_UTL_H */
