/*********************************************************************
 * pokecom2.h
 *
 * Copyright (C) 2004 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: pokecom2.h 842 2006-12-04 06:21:20Z sasaoka $ */

#ifndef POKECOM2_H
#define POKECOM2_H

#include "rc.h"

typedef enum {
	PC_INT,
	PC_DOUBLE,
	PC_STRING,
	PC_ARRAY,
	PC_UNKNOWN
} PC_VARIABLE_TYPE;


typedef struct pc_variable{
	char *name;
	PC_VARIABLE_TYPE type;
	union {
		long int i;
		long double d;
		char *s;
		struct pc_variable *array;
	} v;
	int array_size;
	int alloc_size;
} PC_VARIABLE;


RC eval_add(FILE *fp, char *buf, int buf_size, int *ptr,
            PC_VARIABLE **varray, int *varray_size, PC_VARIABLE *val,
            int *line, int *subst_flag);
RC pc_print_variable_array(FILE *fp, PC_VARIABLE *varray, int varray_size);
RC pc_free_variable(PC_VARIABLE *var);
PC_VARIABLE pc_init_variable(void);
int pc_search_varray(const char vname[],
                     int varray_size, const PC_VARIABLE varray[]);
RC pc_varray2str(const char vname[],
                 int varray_size, const PC_VARIABLE varray[], char **v);
RC pc_varray2int(const char vname[],
                 int varray_size, const PC_VARIABLE varray[],
                 int *v, int v_min, int v_max, int v_default);
RC pc_variable2int(int index, int varray_size, const PC_VARIABLE varray[],
                   int *v, int v_min, int v_max, int v_default,
                   const char vname[]);
RC pc_varray2double(const char vname[],
                    int varray_size, const PC_VARIABLE varray[],
                    double *v, double v_min, double v_max, double v_default);
RC pc_variable2double(int index, int varray_size, const PC_VARIABLE varray[],
                      double *v, double v_min, double v_max, double v_default,
                      const char vname[]);
RC pc_variable2int_array(int index, int varray_size, const PC_VARIABLE varray[],
                         int v[], int *v_size, int v_size_max,
                         int v_min, int v_max, int v_default,
                         const char vname[]);
RC pc_variable2int_array_var(int index, int varray_size,
                             const PC_VARIABLE varray[], int **v, int *v_size,
                             int v_min, int v_max, int v_default,
                             const char vname[]);

#endif /* POKECOM2_H */


