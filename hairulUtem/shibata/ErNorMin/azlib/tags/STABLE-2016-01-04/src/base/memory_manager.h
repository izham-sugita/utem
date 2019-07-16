/*********************************************************************
 * memory_manager.h
 *
 * Copyright (C) 2004 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Thanks to following contributors.
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: memory_manager.h 842 2006-12-04 06:21:20Z sasaoka $ */

#ifndef MEMORY_MANAGER_H
#define MEMORY_MANAGER_H

#include <stdio.h>
#include "rc.h"
#include "wrapper_64.h"

RC mm_init(size_t size);
RC mm_resize(size_t size);
RC mm_terminate(void);
void *mm_alloc(size_t size);
void *mm_alloc64(WRAP64_INT size);
void *mm_alloc_tmp(size_t size);
void *mm_realloc(void *ptr, size_t size);
void *mm_realloc64(void *ptr, WRAP64_INT size);
RC mm_free(void *ptr);
RC mm_compaction_va(void **ptr, ...);
RC mm_compaction_arr(int ptr_num, void *ptr[]);
RC mm_compaction_arr2(int ptr1_num, void *ptr1[], int ptr2_num, void *ptr2[]);
RC mm_chk_list(FILE *fp);
size_t mm_allocatable(FILE *fp);
size_t mm_last_allocatable(void);
size_t mm_max_usage(void);
size_t mm_size(void);
size_t mm_list_size(void);

#endif /* MEMORY_MANAGER_H */

