/*********************************************************************
 * bin_io.h
 *
 * Copyright (C) 2001 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: bin_io.h 842 2006-12-04 06:21:20Z sasaoka $ */

#ifndef BIN_IO_H
#define BIN_IO_H

#include <stdio.h>
#include "rc.h"

typedef enum {
	LH,
	HL
} ENDIAN;


RC read_2bytes(FILE *fp, ENDIAN endian, unsigned int *ui);
RC read_4bytes(FILE *fp, ENDIAN endian, unsigned long *ul);

RC lread_2bytes(FILE *fp,unsigned int *value);
RC lread_4bytes(FILE *fp,unsigned long *value);
RC lwrite_2bytes(FILE *fp,unsigned int value);
RC lwrite_4bytes(FILE *fp,unsigned long value);
RC hread_2bytes(FILE *fp,unsigned int *value);
RC hread_4bytes(FILE *fp,unsigned long *value);
RC hwrite_2bytes(FILE *fp,unsigned int value);
RC hwrite_4bytes(FILE *fp,unsigned long value);

#endif /* BIN_IO_H */

