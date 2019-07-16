/*********************************************************************
 * pokecom.c
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

/* $Id: pokecom.c,v 1.10 2003/10/08 08:27:04 sasaoka Exp $ */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <limits.h>
#include "rc.h"
#include "math_utl.h"
#include "pokecom.h"
#include "string_utl.h"


static RC eval_mul(const char *str, int *ptr, POKE_VARIABLE **vlist_ptr,
                   POKE_VARIABLE *val);
static RC eval_expo(const char *str, int *ptr, POKE_VARIABLE **vlist_ptr,
                    POKE_VARIABLE *val);
static RC call_func(const char *str, int *ptr, POKE_VARIABLE **vlist_ptr,
                    const char *func_str,
                    RC (*func)(long double, long double *),
                    POKE_VARIABLE *val);
static RC eval_func(const char *str, int *ptr, POKE_VARIABLE **vlist_ptr,
                    POKE_VARIABLE *val);
static RC eval_fact(const char *str, int *ptr, POKE_VARIABLE **vlist_ptr,
                    POKE_VARIABLE *val);
static RC eval_int(const char *str, int *ptr, int *val);
static char *get_name(const char *str, int *ptr);
static int case_strncmp(const char *s1, const char *s2, size_t n);
static int case_strcmp(const char *s1, const char *s2);
static RC val2float(POKE_VARIABLE *val);
static RC val2int(POKE_VARIABLE *val);
static RC call_sin(long double v, long double *ret);
static RC call_cos(long double v, long double *ret);
static RC call_tan(long double v, long double *ret);




static RC read_int(const char *str, int *ptr, int *val);
static RC read_fact(const char *str, int *ptr, long double *val);
static RC expr_add(const char *str, int *ptr, long double *val);
static RC expr_mul(const char *str, int *ptr, long double *val);
static RC expr_expo(const char *str, int *ptr, long double *val);
static RC expr_func(const char *str, int *ptr, long double *val);


double poke_var2double(const char *vname, POKE_VARIABLE *vlist_ptr,
                       double def_val)
{
	double ret = def_val;
	POKE_VARIABLE *ptr = poke_search_variable(vname, vlist_ptr);

	if(ptr != NULL){
		if(ptr->type == POKE_FLOATING_POINT) ret = (double)(ptr->v.ld);
	}

	return(ret);
}


int poke_var2int(const char *vname, POKE_VARIABLE *vlist_ptr,
                 int def_val)
{
	int ret = def_val;
	POKE_VARIABLE *ptr = poke_search_variable(vname, vlist_ptr);

	if(ptr != NULL){
		if(ptr->type == POKE_INTEGER) ret = (int)(ptr->v.li);
	}

	return(ret);
}


POKE_VARIABLE *poke_search_add_variable(const char *vname,
                                        POKE_VARIABLE *vlist_ptr)
{
	POKE_VARIABLE *pre_ptr = NULL;
	POKE_VARIABLE *alloc_ptr = NULL;

	while(vlist_ptr != NULL){
		if(case_strcmp(vname, vlist_ptr->name) == 0) return(vlist_ptr);

		pre_ptr = vlist_ptr;
		vlist_ptr = vlist_ptr->next;
	}

	alloc_ptr = (POKE_VARIABLE *)malloc(sizeof(POKE_VARIABLE));
	if(alloc_ptr == NULL) return(NULL);
	alloc_ptr->name = (char *)malloc(strlen(vname)+1);
	if(alloc_ptr->name == NULL) return(NULL);
	strcpy(alloc_ptr->name, vname);
	alloc_ptr->type = POKE_INTEGER;
	alloc_ptr->v.li = 0L;
	alloc_ptr->next = NULL;
	if(pre_ptr != NULL) pre_ptr->next = alloc_ptr;

	return(alloc_ptr);
}


POKE_VARIABLE *poke_search_variable(const char *vname,
                                    POKE_VARIABLE *vlist_ptr)
{
	while(vlist_ptr != NULL){
		if(case_strcmp(vname, vlist_ptr->name) == 0) break;

		vlist_ptr = vlist_ptr->next;
	}

	return(vlist_ptr);
}


RC print_poke_variables(FILE *fp, POKE_VARIABLE *vlist_ptr)
{
	while(vlist_ptr != NULL){
		print_poke_variable(fp, *vlist_ptr);
		vlist_ptr = vlist_ptr->next;
	}
	
	return(NORMAL_RC);
}


RC print_poke_variable(FILE *fp, POKE_VARIABLE val)
{
	int ii1;

	if(val.name){
		fprintf(fp, "[%s:", val.name);
	}else{
		fprintf(fp, "[-:");
	}

	if(val.type == POKE_INTEGER){
		fprintf(fp, " I]: %ld\n", val.v.li);
	}else if(val.type == POKE_FLOATING_POINT){
		fprintf(fp, " F]: %25.17Le\n", val.v.ld);
	}else if(val.type == POKE_FLOATING_POINT_ARRAY){
		fprintf(fp, " F%d]:", val.v.ldarr.size);
		for(ii1=0; ii1<val.v.ldarr.size; ii1++){
			fprintf(fp, " %15.7Le", val.v.ldarr.v[ii1]);
		}
		fprintf(fp, "\n");
	}else if(val.type == POKE_CHARACTER_STRING){
		fprintf(fp, " S]: %s\n", val.v.c);
	}else if(val.type == POKE_UNKNOWN){
		fprintf(fp, " UNKNOWN]\n");
	}else{
		return(IMPLEMENT_ERROR_RC);;
	}

	return(NORMAL_RC);
}


RC free_poke_variables(POKE_VARIABLE *vlist_ptr)
{
	POKE_VARIABLE *old_ptr;

	while(vlist_ptr != NULL){
		if(vlist_ptr->name) free(vlist_ptr->name);
		vlist_ptr->name = NULL;
		if(vlist_ptr->type == POKE_FLOATING_POINT_ARRAY){
			if(vlist_ptr->v.ldarr.v != NULL) free(vlist_ptr->v.ldarr.v);
			vlist_ptr->v.ldarr.v = NULL;
		}
		old_ptr = vlist_ptr;
		vlist_ptr = vlist_ptr->next;
		free(old_ptr);
	}
	
	return(NORMAL_RC);
}


RC init_poke_variable(POKE_VARIABLE *val)
{
	val->name = NULL;
	val->type = POKE_UNKNOWN;
	val->next = NULL;

	return(NORMAL_RC);
}


RC eval_expr(const char *str, long double *v)
{
	int ptr = 0;
	POKE_VARIABLE *vlist_ptr = NULL;
	POKE_VARIABLE *old_ptr = NULL;
	POKE_VARIABLE val;

	RC_TRY( eval_add(str, &ptr, &vlist_ptr, &val) );
	RC_TRY( val2float(&val) );
	*v = val.v.ld;

	while(vlist_ptr != NULL){
		fprintf(stdout, "%s:", vlist_ptr->name);
		free(vlist_ptr->name);
		if(vlist_ptr->type == POKE_INTEGER){
			fprintf(stdout, "I:%d\n", (int)(vlist_ptr->v.li));
		}else if(vlist_ptr->type == POKE_FLOATING_POINT){
			fprintf(stdout, "F:%25.17Le\n", (vlist_ptr->v.ld));
		}else{
			fprintf(stdout, "???\n");
		}
		old_ptr = vlist_ptr;
		vlist_ptr = vlist_ptr->next;
		free(old_ptr);
	}

	return(NORMAL_RC);
}


RC eval_add(const char *str, int *ptr, POKE_VARIABLE **vlist_ptr,
            POKE_VARIABLE *val)
{
	POKE_VARIABLE v1, v2;

	RC_TRY( init_poke_variable(val) );
	RC_TRY( eval_mul(str, ptr, vlist_ptr, &v1) );
	while(1){
		while(str[*ptr] == (char)' ') (*ptr)++;

		if( (str[*ptr] == (char)'+')
		  ||(str[*ptr] == (char)'-') ){
			char cal_symbol = str[*ptr];

			(*ptr)++;
			RC_TRY( eval_mul(str, ptr, vlist_ptr, &v2) );

			if( (v1.type == POKE_FLOATING_POINT)
			  ||(v2.type == POKE_FLOATING_POINT) ){
				RC_TRY( val2float(&v1) );
				RC_TRY( val2float(&v2) );
				if(cal_symbol == '+'){
					v1.v.ld += v2.v.ld;
				}else{
					v1.v.ld -= v2.v.ld;
				}
			}else if( (v1.type == POKE_INTEGER)
			        &&(v2.type == POKE_INTEGER) ){
				if(cal_symbol == '+'){
					v1.v.li += v2.v.li;
				}else{
					v1.v.li -= v2.v.li;
				}
			}else{
				return(CAL_ERROR_RC);
			}
		}else{
			break;
		}
	}
	*val = v1;

	return(NORMAL_RC);
}


static RC eval_mul(const char *str, int *ptr, POKE_VARIABLE **vlist_ptr,
                   POKE_VARIABLE *val)
{
	POKE_VARIABLE v1, v2;

	RC_TRY( init_poke_variable(val) );
	RC_TRY( eval_expo(str, ptr, vlist_ptr, &v1) );
	while(1){
		while(str[*ptr] == (char)' ') (*ptr)++;

		if( (str[*ptr] == (char)'*')
		  ||(str[*ptr] == (char)'/') ){
			char cal_symbol = str[*ptr];

			(*ptr)++;
			RC_TRY( eval_expo(str, ptr, vlist_ptr, &v2) );

			if( (v1.type == POKE_FLOATING_POINT)
			  ||(v2.type == POKE_FLOATING_POINT) ){
				RC_TRY( val2float(&v1) );
				RC_TRY( val2float(&v2) );
				if(cal_symbol == '*'){
					v1.v.ld *= v2.v.ld;
				}else{
					if( nearly_eq((double)(v2.v.ld), 0.0) ){
						return(CAL_ERROR_RC);
					}
					v1.v.ld /= v2.v.ld;
				}
			}else if( (v1.type == POKE_INTEGER)
			        &&(v2.type == POKE_INTEGER) ){
				if(cal_symbol == '*'){
					v1.v.li *= v2.v.li;
				}else{
					if(v2.v.li == 0L) return(CAL_ERROR_RC);
					v1.v.li /= v2.v.li;
				}
			}else{
				return(CAL_ERROR_RC);
			}
		}else{
			break;
		}
	}
	*val = v1;

	return(NORMAL_RC);
}


static RC eval_expo(const char *str, int *ptr, POKE_VARIABLE **vlist_ptr,
                    POKE_VARIABLE *val)
{
	POKE_VARIABLE v1;

	RC_TRY( init_poke_variable(val) );
	RC_TRY( eval_func(str, ptr, vlist_ptr, val) );

	while((int)(str[*ptr]) == ' ') (*ptr)++;

	if((int)(str[*ptr]) == '^'){
		(*ptr)++;
	}else if( ((int)(str[*ptr]) == '*')&&((int)(str[(*ptr)+1]) == '*') ){
		(*ptr) += 2;
	}else{
		return(NORMAL_RC);
	}

	RC_TRY( eval_expo(str, ptr, vlist_ptr, &v1) );
	RC_TRY( val2float(val) );
	RC_TRY( val2float(&v1) );
	val->v.ld = (long double)pow((double)val->v.ld, (double)v1.v.ld);

	return(NORMAL_RC);

}


static RC call_func(const char *str, int *ptr, POKE_VARIABLE **vlist_ptr,
                    const char *func_str,
                    RC (*func)(long double, long double *),
                    POKE_VARIABLE *val)
{
	int func_len = strlen(func_str);
	int next_ptr = (*ptr) + func_len;
	POKE_VARIABLE v1;

	RC_TRY( init_poke_variable(val) );
	if(func_len <= 0) return(ARG_ERROR_RC);
	if(case_strncmp(func_str, &(str[*ptr]), func_len) != 0) return(SPECIAL_RC);
	if( (((char)'A' <= str[next_ptr])&&(str[next_ptr] <= (char)'Z'))
	  ||(((char)'a' <= str[next_ptr])&&(str[next_ptr] <= (char)'z'))
	  ||(str[next_ptr] == (char)'_') ){
		return(SPECIAL_RC);
	}

	*ptr = next_ptr;
	RC_TRY( eval_func(str, ptr, vlist_ptr, &v1) );
	RC_TRY( val2float(&v1) );
	val->type = POKE_FLOATING_POINT;
	RC_TRY( func(v1.v.ld, &(val->v.ld)) );

	return(NORMAL_RC);
}


static RC call_sin(long double v, long double *ret)
{
#ifndef C99
	*ret = (long double)sin(v);
#else   /* C99 */
	*ret = sinl(v);
#endif  /* C99 */

	return(NORMAL_RC);
}

static RC call_cos(long double v, long double *ret)
{
#ifndef C99
	*ret = (long double)cos(v);
#else   /* C99 */
	*ret = cosl(v);
#endif  /* C99 */

	return(NORMAL_RC);
}

static RC call_tan(long double v, long double *ret)
{
#ifndef C99
	*ret = (long double)tan(v);
#else   /* C99 */
	*ret = tanl(v);
#endif  /* C99 */

	return(NORMAL_RC);
}

#define RC_TRY_SP(func)  {RC lrc; lrc = func;\
                          if(lrc == NORMAL_RC) return(NORMAL_RC);\
                          if(lrc != SPECIAL_RC){\
                              fprintf(stderr, #func" <\n"); return(lrc);\
                         }}


static RC eval_func(const char *str, int *ptr, POKE_VARIABLE **vlist_ptr,
                    POKE_VARIABLE *val)
{
	POKE_VARIABLE v1;
	int iv1, iv2;
	char open_ch;

	RC_TRY( init_poke_variable(val) );
	while((int)(str[*ptr]) == ' ') (*ptr)++;
	
	RC_TRY_SP( call_func(str, ptr, vlist_ptr,
	                     "SIN", call_sin, val) );
	
	RC_TRY_SP( call_func(str, ptr, vlist_ptr,
	                     "COS", call_cos, val) );
	
	RC_TRY_SP( call_func(str, ptr, vlist_ptr,
	                     "TAN", call_tan, val) );
	

	if(case_strncmp("NPR", &(str[*ptr]), 3) == 0){
		(*ptr) += 3;
		while((int)(str[*ptr]) == ' ') (*ptr)++;
		if( (str[*ptr] != (char)'(')
		  &&(str[*ptr] != (char)'{')
		  &&(str[*ptr] != (char)'[')
		  &&(str[*ptr] != (char)'<') ) return(READ_ERROR_RC);
		open_ch = str[*ptr];
		(*ptr)++;
		RC_TRY( eval_add(str, ptr, vlist_ptr, &v1) );
		RC_TRY( val2int(&v1) );
		iv1 = v1.v.li;
		while((int)(str[*ptr]) == ' ') (*ptr)++;
		if(str[*ptr] != (char)',') return(READ_ERROR_RC);
		(*ptr)++;
		RC_TRY( eval_add(str, ptr, vlist_ptr, &v1) );
		RC_TRY( val2int(&v1) );
		iv2 = v1.v.li;
		while((int)(str[*ptr]) == ' ') (*ptr)++;
		if( ((open_ch == (char)'(')&&(str[*ptr] != (char)')'))
		  ||((open_ch == (char)'{')&&(str[*ptr] != (char)'}'))
		  ||((open_ch == (char)'[')&&(str[*ptr] != (char)']'))
		  ||((open_ch == (char)'<')&&(str[*ptr] != (char)'>')) ){
			return(READ_ERROR_RC);
		}
		(*ptr)++;
		val->v.ld = permutation(iv1, iv2);
		val->type = POKE_FLOATING_POINT;
		return(NORMAL_RC);
	}

	RC_TRY( eval_fact(str, ptr, vlist_ptr, &v1) );
	*val = v1;

	return(NORMAL_RC);
}


static RC eval_fact(const char *str, int *ptr, POKE_VARIABLE **vlist_ptr,
                    POKE_VARIABLE *val)
{
	int sign = 1;
	int start_flag, stop_flag, point_flag;
	long double point_val;
	int iexp;
	int ii1;

	RC_TRY( init_poke_variable(val) );

	/* 符号の処理 */
	while(1){
		if(str[*ptr] == (char)'-'){
			sign *= -1;
		}else if( (str[*ptr] != (char)'+')
		        &&(str[*ptr] != (char)' ') ){
			break;
		}
		(*ptr)++;
	}

	/* 配列の処理 */
	if(str[*ptr] == (char)'{'){
		(*ptr)++;
		val->type = POKE_FLOATING_POINT_ARRAY;
		val->v.ldarr.size = 0;
		val->v.ldarr.v = NULL;
		while(1){
			POKE_VARIABLE tmp_val;

			while(str[*ptr] == (char)' ') (*ptr)++;
			if(str[*ptr] == (char)'}'){
				(*ptr)++;
				break;
			}
			RC_TRY( eval_add(str, ptr, vlist_ptr, &tmp_val) );
			(val->v.ldarr.size)++;
			val->v.ldarr.v = realloc(val->v.ldarr.v,
			                         (val->v.ldarr.size)*sizeof(long double));
			RC_NULL_CHK( val->v.ldarr.v );
			RC_TRY( val2float(&tmp_val) );
			val->v.ldarr.v[val->v.ldarr.size - 1] = sign*(tmp_val.v.ld);
			while(str[*ptr] == (char)' ') (*ptr)++;
			if(str[*ptr] == (char)',') (*ptr)++;
		}

		return(NORMAL_RC);
	}

	/* 括弧の処理 */
	if( (str[*ptr] == (char)'(')
	  ||(str[*ptr] == (char)'[')
	  ||(str[*ptr] == (char)'<') ){
		char open = str[*ptr];

		(*ptr)++;
		RC_TRY( eval_add(str, ptr, vlist_ptr, val) );

		while(str[*ptr] == (char)' ') (*ptr)++;

		if( ((open == (char)'(')&&(str[*ptr] != (char)')'))
		  ||((open == (char)'[')&&(str[*ptr] != (char)']'))
		  ||((open == (char)'<')&&(str[*ptr] != (char)'>')) ){
			return(READ_ERROR_RC);
		}
		(*ptr)++;

		if((val->type) == POKE_INTEGER){
			val->v.li *= (long int)sign;
		}else if((val->type) == POKE_FLOATING_POINT){
			val->v.ld *= (long double)sign;
		}else{
			return(READ_ERROR_RC);
		}

		return(NORMAL_RC);
	}

	/* 変数の処理 */
	if( (((char)'A' <= str[*ptr])&&(str[*ptr] <= (char)'Z'))
	  ||(((char)'a' <= str[*ptr])&&(str[*ptr] <= (char)'z'))
	  ||(str[*ptr] == (char)'_') ){
		char *vname;
		POKE_VARIABLE *v_ptr;

		RC_NULL_CHK( vname = get_name(str, ptr) );
		while(str[*ptr] == (char)' ') (*ptr)++;

		if(str[*ptr] == (char)'='){
			/* 代入演算子の処理 */
			(*ptr)++;
			RC_TRY( eval_add(str, ptr, vlist_ptr, val) );
			v_ptr = poke_search_add_variable(vname, *vlist_ptr);
			if((*vlist_ptr) == NULL) (*vlist_ptr) = v_ptr;
			free(vname);
			v_ptr->type = val->type;
			v_ptr->v = val->v;
			if(val->type == POKE_FLOATING_POINT_ARRAY){
				v_ptr->v.ldarr.v = malloc( (val->v.ldarr.size)
				                           *sizeof(long double) );
				RC_NULL_CHK(v_ptr->v.ldarr.v);
				for(ii1=0; ii1<(val->v.ldarr.size); ii1++){
					v_ptr->v.ldarr.v[ii1] = val->v.ldarr.v[ii1];
				}
			}
		}else{
			v_ptr = poke_search_variable(vname, *vlist_ptr);
			if(v_ptr == NULL){
				RC lrc = SEEK_ERROR_RC;
				rc_error_print(__FILE__, __LINE__, vname, &lrc);
				free(vname); 
				return(READ_ERROR_RC);
			}
			free(vname);
			val->type = v_ptr->type;
			val->v = v_ptr->v;
		}

		if((val->type) == POKE_INTEGER){
			val->v.li *= (long int)sign;
		}else if((val->type) == POKE_FLOATING_POINT){
			val->v.ld *= (long double)sign;
		}else if((val->type) == POKE_FLOATING_POINT_ARRAY){
			for(ii1=0; ii1<val->v.ldarr.size; ii1++){
				val->v.ldarr.v[ii1] *= (long double)sign;
			}
		}else{
			return(READ_ERROR_RC);
		}

		return(NORMAL_RC);
	}

	(val->type) = POKE_INTEGER;
	(val->v.li) = 0;
	start_flag = 0;
	stop_flag = 0;
	point_val = 1.0L;
	point_flag = 0;
	while(1){
		switch((int)(str[*ptr])){
		case '0':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7':
		case '8':
		case '9':
			start_flag = 1;
			if(point_flag){
				point_val *= 0.1L;
				val->v.ld += point_val * (long double)((int)(str[*ptr]) - '0');
			}else{
				if( ((val->type) == POKE_INTEGER)
				  &&( ((val->v.li) >= (LONG_MAX/10L))
				    ||((val->v.li) <= (LONG_MIN/10L)) ) ){
					RC_TRY( val2float(val) );
				}
				if((val->type) == POKE_INTEGER){
					(val->v.li) *= 10L;
					(val->v.li) += (long int)((int)(str[*ptr]) - '0');
				}else{
					(val->v.ld) *= 10.0L;
					(val->v.ld) += (long double)((int)(str[*ptr]) - '0');
				}
			}
			break;
		case '.':
			if(point_flag){
				stop_flag = 1;
				break;
			}
			RC_TRY( val2float(val) );
			start_flag = 1;
			point_flag = 1;
			break;
		case 'E':
		case 'e':
			RC_TRY( val2float(val) );
			(*ptr)++;
			RC_TRY( eval_int(str, ptr, &iexp) );
			(val->v.ld) *= (long double)(pow(10.0, (double)iexp));
			start_flag = 1;
			stop_flag = 1;
			break;
		default:
			stop_flag = 1;
		}
		if(stop_flag) break;
		(*ptr)++;
	}
	if(start_flag == 0){
		/* 変換できなかった時 */
		(val->type) = POKE_UNKNOWN;
		(val->v.li) = LONG_MIN;
		return(NORMAL_RC);
	}
	
	if((val->type) == POKE_INTEGER){
		val->v.li *= (long int)sign;
	}else if((val->type) == POKE_FLOATING_POINT){
		val->v.ld *= (long double)sign;
	}else{
		return(READ_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC eval_int(const char *str, int *ptr, int *val)
{
	int sign = 1;
	int start_flag = 0;
	int stop_flag = 0;

	while((int)(str[*ptr]) == ' ') (*ptr)++;

	*val = 0;
	while(1){
		switch((int)(str[*ptr])){
		case '+':
			if(start_flag) stop_flag = 1;
			break;
		case '-':
			if(start_flag){
				stop_flag = 1;
			}else{
				sign *= -1;
			}
			break;
		case '0':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7':
		case '8':
		case '9':
			start_flag = 1;
			*val *= 10;
			*val += (int)(str[*ptr]) - '0';
			break;
		case ' ':
		case '\t':
			break;
		default:
			stop_flag = 1;
		}
		if(stop_flag) break;
		(*ptr)++;
	}
	if(start_flag == 0) return(READ_ERROR_RC);

	*val *= sign;

	return(NORMAL_RC);
}


static char *get_name(const char *str, int *ptr)
{
	int name_size = 0;
	int start_ptr = *ptr;
	char *ret = NULL;
	int ii1;
		
	while( (((char)'A' <= str[*ptr])&&(str[*ptr] <= (char)'Z'))
	     ||(((char)'a' <= str[*ptr])&&(str[*ptr] <= (char)'z'))
	     ||(((char)'0' <= str[*ptr])&&(str[*ptr] <= (char)'9'))
	     ||(str[*ptr] == (char)'_') ){
		(*ptr)++;
		name_size++;
	}

	ret = (char *)malloc(sizeof(char)*(name_size + 1));
	if(ret == NULL) return(NULL);

	for(ii1=0; ii1<name_size; ii1++){
		ret[ii1] = str[start_ptr + ii1];
	}
	ret[name_size] = (char)'\0';

	return(ret);
}


static int case_strcmp(const char *s1, const char *s2)
{
	int ii1 = 0;

	while( (s1[ii1] != (char)'\0')
	     &&(s2[ii1] != (char)'\0') ){
		int up1 = toupper((int)(s1[ii1]));
		int up2 = toupper((int)(s2[ii1]));

		if(up1 != up2) return( up1 - up2 );
		ii1++;
	}

	return(0);
}


static int case_strncmp(const char *s1, const char *s2, size_t n)
{
	int ii1;

	for(ii1=0; ii1<(int)n; ii1++){
		int up1 = toupper((int)(s1[ii1]));
		int up2 = toupper((int)(s2[ii1]));

		if(up1 != up2) return( up1 - up2 );
		if(s1[ii1] == (char)'\0') return(0);
	}

	return(0);
}


RC poke_val2float(POKE_VARIABLE *val)
{
	if(val == NULL) return(NULL_ERROR_RC);
	return(val2float(val));
}


static RC val2float(POKE_VARIABLE *val)
{
	if((val->type) == POKE_FLOATING_POINT){
		return(NORMAL_RC);
	}else if((val->type) == POKE_FLOATING_POINT_ARRAY){
		long double sum = 0.0;
		int ii1;
		for(ii1=0; ii1<val->v.ldarr.size; ii1++){
			sum += val->v.ldarr.v[ii1];
		}
		if(val->v.ldarr.v != NULL) free(val->v.ldarr.v);
		(val->v.ld) = sum;
	}else if((val->type) == POKE_INTEGER){
		(val->v.ld) = (long double)(val->v.li);
	}else if((val->type) == POKE_UNKNOWN){
		(val->v.ld) = (long double)(DBL_MIN);
	}else{
		return(ARG_ERROR_RC);
	}
	(val->type) = POKE_FLOATING_POINT;

	return(NORMAL_RC);
}


static RC val2int(POKE_VARIABLE *val)
{
	if((val->type) == POKE_INTEGER){
		return(NORMAL_RC);
	}if((val->type) == POKE_FLOATING_POINT_ARRAY){
		RC_TRY( val2float(val) );
		(val->v.li) = (long int)(ceil((double)(val->v.ld)) + 0.0001);
	}else if((val->type) == POKE_FLOATING_POINT){
		(val->v.li) = (long int)(ceil((double)(val->v.ld)) + 0.0001);
	}else if((val->type) == POKE_UNKNOWN){
		(val->v.ld) = (long int)(INT_MIN);
	}else{
		return(ARG_ERROR_RC);
	}
	(val->type) = POKE_INTEGER;

	return(NORMAL_RC);
}




static RC expr_add(const char *str, int *ptr, long double *val)
{
	long double v1, v2;
	int stop_flag = 0;

	while((int)(str[*ptr]) == ' ') (*ptr)++;

	RC_TRY( expr_mul(str, ptr, &v1) );
	while(1){
		switch((int)(str[*ptr])){
		case '+':
			(*ptr)++;
			RC_TRY( expr_mul(str, ptr, &v2) );
			v1 += v2;
			break;
		case '-':
			(*ptr)++;
			RC_TRY( expr_mul(str, ptr, &v2) );
			v1 -= v2;
			break;
		case ' ':
			(*ptr)++;
			break;
		default:
			stop_flag = 1;
		}
		if(stop_flag) break;
	}

	*val = v1;

	return(NORMAL_RC);
}


static RC expr_mul(const char *str, int *ptr, long double *val)
{
	long double v1, v2;
	int stop_flag = 0;

	while((int)(str[*ptr]) == ' ') (*ptr)++;

	RC_TRY( expr_expo(str, ptr, &v1) );
	while(1){
		switch((int)(str[*ptr])){
		case '*':
			(*ptr)++;
			RC_TRY( expr_expo(str, ptr, &v2) );
			v1 *= v2;
			break;
		case '/':
			(*ptr)++;
			RC_TRY( expr_expo(str, ptr, &v2) );
			if( nearly_eq(v2, 0.0) ) return(CAL_ERROR_RC);
			v1 /= v2;
			break;
		case ' ':
			(*ptr)++;
			break;
		default:
			stop_flag = 1;
		}
		if(stop_flag) break;
	}

	*val = v1;

	return(NORMAL_RC);
}


static RC expr_expo(const char *str, int *ptr, long double *val)
{
	long double v1, v2;
	int stop_flag = 0;

	while((int)(str[*ptr]) == ' ') (*ptr)++;

	RC_TRY( expr_func(str, ptr, &v1) );
	while(1){
		switch((int)(str[*ptr])){
		case '^':
			(*ptr)++;
			RC_TRY( expr_expo(str, ptr, &v2) );
			v1 = pow(v1,v2);
			break;
		case '*':
			if(str[(*ptr)+1] != (char)'*'){
				stop_flag = 1;
				break;
			}
			(*ptr) += 2;
			RC_TRY( expr_expo(str, ptr, &v2) );
			v1 = pow(v1,v2);
			break;
		case ' ':
			(*ptr)++;
			break;
		default:
			stop_flag = 1;
		}
		if(stop_flag) break;
	}

	*val = v1;

	return(NORMAL_RC);
}


static RC expr_func(const char *str, int *ptr, long double *val)
{
	long double v1;
	int iv1, iv2;
	char open_ch;

	while((int)(str[*ptr]) == ' ') (*ptr)++;

	if(case_strncmp("SINH", &(str[*ptr]), 4) == 0){
		(*ptr) += 4;
		RC_TRY( expr_func(str, ptr, &v1) );
		*val = sinh(v1);
		return(NORMAL_RC);
	}

	if(case_strncmp("HSN", &(str[*ptr]), 3) == 0){
		(*ptr) += 3;
		RC_TRY( expr_func(str, ptr, &v1) );
		*val = sinh(v1);
		return(NORMAL_RC);
	}

	if(case_strncmp("COSH", &(str[*ptr]), 4) == 0){
		(*ptr) += 4;
		RC_TRY( expr_func(str, ptr, &v1) );
		*val = cosh(v1);
		return(NORMAL_RC);
	}

	if(case_strncmp("HCS", &(str[*ptr]), 3) == 0){
		(*ptr) += 3;
		RC_TRY( expr_func(str, ptr, &v1) );
		*val = cosh(v1);
		return(NORMAL_RC);
	}

	if(case_strncmp("TANH", &(str[*ptr]), 4) == 0){
		(*ptr) += 4;
		RC_TRY( expr_func(str, ptr, &v1) );
		*val = tanh(v1);
		return(NORMAL_RC);
	}

	if(case_strncmp("HTN", &(str[*ptr]), 3) == 0){
		(*ptr) += 3;
		RC_TRY( expr_func(str, ptr, &v1) );
		*val = tanh(v1);
		return(NORMAL_RC);
	}

	if(case_strncmp("SIN", &(str[*ptr]), 3) == 0){
		(*ptr) += 3;
		RC_TRY( expr_func(str, ptr, &v1) );
		*val = sin(v1);
		return(NORMAL_RC);
	}

	if(case_strncmp("ASIN", &(str[*ptr]), 4) == 0){
		(*ptr) += 4;
		RC_TRY( expr_func(str, ptr, &v1) );
		if((v1 > 1.0L)||(v1 < -1.0L)) return(CAL_ERROR_RC);
		*val = asin(v1);
		return(NORMAL_RC);
	}

	if(case_strncmp("ASN", &(str[*ptr]), 3) == 0){
		(*ptr) += 3;
		RC_TRY( expr_func(str, ptr, &v1) );
		if((v1 > 1.0L)||(v1 < -1.0L)) return(CAL_ERROR_RC);
		*val = asin(v1);
		return(NORMAL_RC);
	}

	if(case_strncmp("COS", &(str[*ptr]), 3) == 0){
		(*ptr) += 3;
		RC_TRY( expr_func(str, ptr, &v1) );
		*val = cos(v1);
		return(NORMAL_RC);
	}

	if(case_strncmp("ACOS", &(str[*ptr]), 4) == 0){
		(*ptr) += 4;
		RC_TRY( expr_func(str, ptr, &v1) );
		if((v1 > 1.0L)||(v1 < -1.0L)) return(CAL_ERROR_RC);
		*val = acos(v1);
		return(NORMAL_RC);
	}

	if(case_strncmp("ACS", &(str[*ptr]), 3) == 0){
		(*ptr) += 3;
		RC_TRY( expr_func(str, ptr, &v1) );
		if((v1 > 1.0L)||(v1 < -1.0L)) return(CAL_ERROR_RC);
		*val = acos(v1);
		return(NORMAL_RC);
	}

	if(case_strncmp("TAN", &(str[*ptr]), 3) == 0){
		(*ptr) += 3;
		RC_TRY( expr_func(str, ptr, &v1) );
		*val = tan(v1);
		return(NORMAL_RC);
	}

	if(case_strncmp("ATAN", &(str[*ptr]), 4) == 0){
		(*ptr) += 4;
		RC_TRY( expr_func(str, ptr, &v1) );
		*val = atan(v1);
		return(NORMAL_RC);
	}

	if(case_strncmp("ATN", &(str[*ptr]), 3) == 0){
		(*ptr) += 3;
		RC_TRY( expr_func(str, ptr, &v1) );
		*val = atan(v1);
		return(NORMAL_RC);
	}

	if(case_strncmp("FACT", &(str[*ptr]), 4) == 0){
		(*ptr) += 4;
		RC_TRY( expr_func(str, ptr, &v1) );
		RC_TRY( double2int_pos((double)v1, &iv1) );
		*val = factorial(iv1);
		return(NORMAL_RC);
	}

	if(case_strncmp("RCP", &(str[*ptr]), 3) == 0){
		(*ptr) += 3;
		RC_TRY( expr_func(str, ptr, &v1) );
		if( nearly_eq((double)v1, 0.0) )return(CAL_ERROR_RC);
		*val = 1.0/v1;
		return(NORMAL_RC);
	}

	if(case_strncmp("EXP", &(str[*ptr]), 3) == 0){
		(*ptr) += 3;
		RC_TRY( expr_func(str, ptr, &v1) );
		*val = exp(v1);
		return(NORMAL_RC);
	}

	if(case_strncmp("LN", &(str[*ptr]), 2) == 0){
		(*ptr) += 2;
		RC_TRY( expr_func(str, ptr, &v1) );
		*val = log(v1);
		return(NORMAL_RC);
	}

	if(case_strncmp("TEN", &(str[*ptr]), 3) == 0){
		(*ptr) += 3;
		RC_TRY( expr_func(str, ptr, &v1) );
		*val = pow(10.0, v1);
		return(NORMAL_RC);
	}

	if(case_strncmp("LOG", &(str[*ptr]), 3) == 0){
		(*ptr) += 3;
		RC_TRY( expr_func(str, ptr, &v1) );
		*val = log10(v1);
		return(NORMAL_RC);
	}

	if(case_strncmp("SQRT", &(str[*ptr]), 4) == 0){
		(*ptr) += 4;
		RC_TRY( expr_func(str, ptr, &v1) );
		if(v1 < 0.0L) return(CAL_ERROR_RC);
		*val = sqrt(v1);
		return(NORMAL_RC);
	}

	if(case_strncmp("SQR", &(str[*ptr]), 3) == 0){
		(*ptr) += 3;
		RC_TRY( expr_func(str, ptr, &v1) );
		if(v1 < 0.0L) return(CAL_ERROR_RC);
		*val = sqrt(v1);
		return(NORMAL_RC);
	}

	if(case_strncmp("CBRT", &(str[*ptr]), 4) == 0){
		(*ptr) += 4;
		RC_TRY( expr_func(str, ptr, &v1) );
		/* ISO C99 READY?   */     
		/* *val = cbrt(v1); */
		*val = pow(v1, 1.0/3.0);

		return(NORMAL_RC);
	}

	if(case_strncmp("CUR", &(str[*ptr]), 3) == 0){
		(*ptr) += 3;
		RC_TRY( expr_func(str, ptr, &v1) );
		/* ISO C99 READY?   */     
		/* *val = cbrt(v1); */
		*val = pow(v1, 1.0/3.0);
		return(NORMAL_RC);
	}

	if(case_strncmp("SQU", &(str[*ptr]), 3) == 0){
		(*ptr) += 3;
		RC_TRY( expr_func(str, ptr, &v1) );
		*val = v1 * v1;
		return(NORMAL_RC);
	}

	if(case_strncmp("CUB", &(str[*ptr]), 3) == 0){
		(*ptr) += 3;
		RC_TRY( expr_func(str, ptr, &v1) );
		*val = v1 * v1 * v1;
		return(NORMAL_RC);
	}

	if(case_strncmp("NPR", &(str[*ptr]), 3) == 0){
		(*ptr) += 3;
		while((int)(str[*ptr]) == ' ') (*ptr)++;
		if( (str[*ptr] != (char)'(')
		  &&(str[*ptr] != (char)'{')
		  &&(str[*ptr] != (char)'[')
		  &&(str[*ptr] != (char)'<') ) return(READ_ERROR_RC);
		open_ch = str[*ptr];
		(*ptr)++;
		RC_TRY( expr_add(str, ptr, &v1) );
		RC_TRY( double2int_pos((double)v1, &iv1) );
		while((int)(str[*ptr]) == ' ') (*ptr)++;
		if(str[*ptr] != (char)',') return(READ_ERROR_RC);
		(*ptr)++;
		RC_TRY( expr_add(str, ptr, &v1) );
		RC_TRY( double2int_pos((double)v1, &iv2) );
		while((int)(str[*ptr]) == ' ') (*ptr)++;
		if( ((open_ch == (char)'(')&&(str[*ptr] != (char)')'))
		  ||((open_ch == (char)'{')&&(str[*ptr] != (char)'}'))
		  ||((open_ch == (char)'[')&&(str[*ptr] != (char)']'))
		  ||((open_ch == (char)'<')&&(str[*ptr] != (char)'>')) ){
			return(READ_ERROR_RC);
		}
		(*ptr)++;
		*val = permutation(iv1, iv2);
		return(NORMAL_RC);
	}

	if(case_strncmp("NCR", &(str[*ptr]), 3) == 0){
		(*ptr) += 3;
		while((int)(str[*ptr]) == ' ') (*ptr)++;
		if( (str[*ptr] != (char)'(')
		  &&(str[*ptr] != (char)'{')
		  &&(str[*ptr] != (char)'[')
		  &&(str[*ptr] != (char)'<') ) return(READ_ERROR_RC);
		open_ch = str[*ptr];
		(*ptr)++;
		while((int)(str[*ptr]) == ' ') (*ptr)++;
		RC_TRY( expr_add(str, ptr, &v1) );
		RC_TRY( double2int_pos((double)v1, &iv1) );
		while((int)(str[*ptr]) == ' ') (*ptr)++;
		if(str[*ptr] != (char)',') return(READ_ERROR_RC);
		(*ptr)++;
		RC_TRY( expr_add(str, ptr, &v1) );
		RC_TRY( double2int_pos((double)v1, &iv2) );
		while((int)(str[*ptr]) == ' ') (*ptr)++;
		if( ((open_ch == (char)'(')&&(str[*ptr] != (char)')'))
		  ||((open_ch == (char)'{')&&(str[*ptr] != (char)'}'))
		  ||((open_ch == (char)'[')&&(str[*ptr] != (char)']'))
		  ||((open_ch == (char)'<')&&(str[*ptr] != (char)'>')) ){
			return(READ_ERROR_RC);
		}
		(*ptr)++;
		*val = combination(iv1, iv2);
		return(NORMAL_RC);
	}

	RC_TRY( read_fact(str, ptr, &v1) );
	*val = v1;

	return(NORMAL_RC);
}


static RC read_fact(const char *str, int *ptr, long double *val)
{
	int exp;
	int sign = 1;
	int start_flag = 0;
	int point_flag = 0;
	int stop_flag = 0;
	long double point_val = 1.0L;

	while((int)(str[*ptr]) == ' ') (*ptr)++;

	if(case_strncmp("PI", &(str[*ptr]), 2) == 0){
		(*ptr) += 2;
		*val = 3.1415926535897932384626433832795029L;
		return(NORMAL_RC);
	}

	if(case_strncmp("E", &(str[*ptr]), 1) == 0){
		(*ptr) += 1;
		*val = 2.7182818284590452353602874713526625L;
		return(NORMAL_RC);
	}

	*val = 0.0;
	while(1){

		switch((int)(str[*ptr])){
		case '+':
			if(start_flag) stop_flag = 1;
			break;
		case '-':
			if(start_flag){
				stop_flag = 1;
			}else{
				sign *= -1;
			}
			break;
		case '(':
		case '{':
		case '[':
		case '<':
			if(!start_flag){
				char open = str[*ptr];

				(*ptr)++;
				RC_TRY( expr_add(str, ptr, val) );
				if( ((open == (char)'(')&&(str[*ptr] != (char)')'))
				  ||((open == (char)'{')&&(str[*ptr] != (char)'}'))
				  ||((open == (char)'[')&&(str[*ptr] != (char)']'))
				  ||((open == (char)'<')&&(str[*ptr] != (char)'>')) ){
					return(READ_ERROR_RC);
				}
				start_flag = 1;
				(*ptr)++;
			}
			stop_flag = 1;
			break;
		case '0':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7':
		case '8':
		case '9':
			start_flag = 1;
			if(point_flag){
				point_val *= 0.1L;
				*val += point_val * (long double)((int)(str[*ptr]) - '0');
			}else{
				*val *= 10.0L;
				*val += (long double)((int)(str[*ptr]) - '0');
			}
			break;
		case '.':
			if(point_flag){
				stop_flag = 1;
				break;
			}
			start_flag = 1;
			point_flag = 1;
			break;
		case 'E':
		case 'e':
			if(start_flag){
				(*ptr)++;
				RC_TRY( read_int(str, ptr, &exp) );
				*val *= (long double)(pow(10.0, (double)exp));
			}
			stop_flag = 1;
			break;
		case ' ':
		case '\t':
			break;
		default:
			stop_flag = 1;
		}
		if(stop_flag) break;
		(*ptr)++;
	}
	if(start_flag == 0) return(READ_ERROR_RC);
	
	*val *= (long double)sign;

	return(NORMAL_RC);
}


static RC read_int(const char *str, int *ptr, int *val)
{
	int sign = 1;
	int start_flag = 0;
	int stop_flag = 0;

	while((int)(str[*ptr]) == ' ') (*ptr)++;

	*val = 0;
	while(1){
		switch((int)(str[*ptr])){
		case '+':
			if(start_flag) stop_flag = 1;
			break;
		case '-':
			if(start_flag){
				stop_flag = 1;
			}else{
				sign *= -1;
			}
			break;
		case '0':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7':
		case '8':
		case '9':
			start_flag = 1;
			*val *= 10;
			*val += (int)(str[*ptr]) - '0';
			break;
		case ' ':
		case '\t':
			break;
		default:
			stop_flag = 1;
		}
		if(stop_flag) break;
		(*ptr)++;
	}
	if(start_flag == 0) return(READ_ERROR_RC);

	*val *= sign;

	return(NORMAL_RC);
}


