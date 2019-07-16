/*********************************************************************
 * pokecom.h
 *
 * Copyright (C) 2002 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Thanks to following contributors.
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: pokecom.h,v 1.6 2003/07/22 10:21:27 adachi Exp $ */

#ifndef POKECOM_H
#define POKECOM_H

#include "string_utl.h"

typedef enum {
	POKE_INTEGER,
	POKE_FLOATING_POINT,
	POKE_FLOATING_POINT_ARRAY,
	POKE_CHARACTER_STRING,
	POKE_UNKNOWN
} POKE_VARIABLE_TYPE;

typedef struct poke_variable{
	char *name;
	POKE_VARIABLE_TYPE type;
	IDC v;
	struct poke_variable *next;
} POKE_VARIABLE;

RC print_poke_variable(FILE *fp, POKE_VARIABLE val);
RC print_poke_variables(FILE *fp, POKE_VARIABLE *vlist_ptr);
RC free_poke_variables(POKE_VARIABLE *vlist_ptr);
RC init_poke_variable(POKE_VARIABLE *val);
double poke_var2double(const char *vname, POKE_VARIABLE *vlist_ptr,
                       double def_val);
int poke_var2int(const char *vname, POKE_VARIABLE *vlist_ptr,
                 int def_val);
POKE_VARIABLE *poke_search_variable(const char *vname, 
                                    POKE_VARIABLE *vlist_ptr);
POKE_VARIABLE *poke_search_add_variable(const char *vname,
                                        POKE_VARIABLE *vlist_ptr);
RC eval_add(const char *str, int *ptr, POKE_VARIABLE **vlist_ptr,
            POKE_VARIABLE *val);
RC poke_val2float(POKE_VARIABLE *val);

RC eval_expr(const char *str, long double *v);


#endif /* POKECOM_H */

