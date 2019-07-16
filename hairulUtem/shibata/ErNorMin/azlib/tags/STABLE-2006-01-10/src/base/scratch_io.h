/*********************************************************************
 * scratch_io.h
 *
 * Copyright (C) 2005 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Thanks to following contributors.
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: scratch_io.h 433 2005-08-29 07:33:36Z sasaoka $ */

#ifndef SCRATCH_IO_H
#define SCRATCH_IO_H

#include <stdio.h>
#include "rc.h"
#include "wrapper_64.h"

#define SCRATCH_DATA_SIZE (256)
typedef struct{
	FILE *fp;
	char *file_name;
	WRAP64_INT position[SCRATCH_DATA_SIZE + 1];
	int use_flag[SCRATCH_DATA_SIZE];  /* 0:未使用, 1:使用, 2:書き込み中 */
	int stack_ptr;
} SCRATCH_FILE;

SCRATCH_FILE *open_scratch_file(const char scratch_dir[]);
RC close_scratch_file(SCRATCH_FILE *sf);
FILE *add_scratch_file_start(SCRATCH_FILE *sf);
int add_scratch_file_end(SCRATCH_FILE *sf);
FILE *read_scratch_file(SCRATCH_FILE *sf, int handle);
FILE *write_scratch_file(SCRATCH_FILE *sf, int handle);
RC delete_scratch_file(SCRATCH_FILE *sf, int handle);
RC print_scratch_file_info(FILE *fp, SCRATCH_FILE *sf);

#endif /* SCRATCH_IO_H */

