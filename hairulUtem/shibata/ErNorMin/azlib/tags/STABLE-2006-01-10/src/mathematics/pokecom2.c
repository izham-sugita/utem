/*********************************************************************
 * pokecom2.c
 *
 * Copyright (C) 2004 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: pokecom2.c 432 2005-08-24 03:00:31Z sasaoka $ */

#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "rc.h"
#include "log_printf.h"
#include "pokecom2.h"
#include "memory_manager.h"
#include "math_utl.h"
#include "string_utl.h"


typedef enum {
	PC_PLUS,
	PC_MINUS,
	PC_MUL,
	PC_DIV,
	PC_POW,
	PC_REMAIN
} PC_OPERATOR;


/* エラー処理用 */
static int FileLine = 0;    /* ファイルの行数 */
static int SubstFlag = 0;   /* 代入・未代入フラグ */

static RC mul_double_variable(long double d, PC_VARIABLE *var);
static RC mul_int_variable(long int i, PC_VARIABLE *var);
static RC copy_variable(PC_VARIABLE src, PC_VARIABLE *dest);
/*
static long int variable2int(PC_VARIABLE var);
*/
static long double variable2double(PC_VARIABLE var);
static RC pc_operation(PC_VARIABLE var1, PC_VARIABLE var2, PC_VARIABLE *ans,
                       PC_OPERATOR operator);
static char *get_name(const char *buf, int *ptr);
static RC eval_fact(FILE *fp, char *buf, int buf_size, int *ptr,
                    PC_VARIABLE **varray, int *varray_size, PC_VARIABLE *val);
static RC eval_int(const char *buf, int *ptr, int *val);
static RC eval_mul(FILE *fp, char *buf, int buf_size, int *ptr,
                   PC_VARIABLE **varray, int *varray_size, PC_VARIABLE *val);
static RC eval_expo(FILE *fp, char *buf, int buf_size, int *ptr,
                    PC_VARIABLE **varray, int *varray_size, PC_VARIABLE *val);
static RC print_variable(FILE *fp, PC_VARIABLE var);
static RC next_token(FILE *fp, char *buf, int buf_size, int *ptr);


static RC
mul_double_variable (long double d, PC_VARIABLE *var)
{
	int ii1;

	if(var->type == PC_INT){
		long double dtmp = (long double)var->v.i;
		var->type = PC_DOUBLE;
		var->v.d = d * dtmp;
	}else if(var->type == PC_DOUBLE){
		var->v.d *= d;
	}else if(var->type == PC_ARRAY){
		for(ii1=0; ii1<(var->array_size); ii1++){
			RC_TRY( mul_double_variable(d, &(var->v.array[ii1])) );
		}
	}

	return(NORMAL_RC);
}


static RC
mul_int_variable (long int i, PC_VARIABLE *var)
{
	int ii1;

	if(var->type == PC_INT){
		if(i != 0){
			if(labs(var->v.i) >= LONG_MAX/labs(i)){
				RC_TRY( warning_printf(11001, FileLine) );
				var->v.d = (long double)(var->v.i);
				var->v.d *= (long double)i;
				var->type = PC_DOUBLE;
			}
		}
		var->v.i *= i;
	}else if(var->type == PC_DOUBLE){
		var->v.d *= (long double)i;
	}else if(var->type == PC_ARRAY){
		for(ii1=0; ii1<(var->array_size); ii1++){
			RC_TRY( mul_int_variable(i, &(var->v.array[ii1])) );
		}
	}

	return(NORMAL_RC);
}


PC_VARIABLE
pc_init_variable (void)
{
	PC_VARIABLE ret;

	ret.name = NULL;
	ret.type = PC_UNKNOWN;
	ret.v.d = (long double)0.0;
	ret.array_size = 0;
	ret.alloc_size = 0;

	return(ret);
}


RC
pc_free_variable (PC_VARIABLE *var)
{
	int ii1;

	if(var == NULL) return(ARG_ERROR_RC);

	if(var->name != NULL){
		RC_TRY( mm_free(var->name) );
		var->name = NULL;
	}
	if(var->type == PC_ARRAY){
		if(var->array_size > var->alloc_size) return(ARG_ERROR_RC);

		for(ii1=0; ii1<var->array_size; ii1++){
			RC_TRY( pc_free_variable(&(var->v.array[ii1])) );
		}
		if(var->v.array != NULL) RC_TRY( mm_free(var->v.array) );
		var->v.array = NULL;
		var->array_size = 0;
		var->alloc_size = 0;
	}else if(var->type == PC_STRING){
		if(var->v.s != NULL){
			RC_TRY( mm_free(var->v.s) );
		}
	}
	*var = pc_init_variable();

	return(NORMAL_RC);
}


static RC
copy_variable (PC_VARIABLE src, PC_VARIABLE *dest)
{
	int ii1;

	*dest = src;
	dest->name = NULL;
	if(src.type == PC_ARRAY){
		if(src.alloc_size > 0){
			if(src.array_size > src.alloc_size) return(ARG_ERROR_RC);

			dest->v.array = mm_alloc(src.alloc_size*sizeof(PC_VARIABLE));
			RC_NULL_CHK( dest->v.array );
			for(ii1=0; ii1<src.array_size; ii1++){
				RC_TRY( copy_variable(src.v.array[ii1],
				                      &(dest->v.array[ii1])) );
			}
		}else{
			src.alloc_size = 0;
			src.array_size = 0;
			src.v.array = NULL;
		}
	}
	if(src.type == PC_STRING){
		RC_NULL_CHK( dest->v.s = mm_alloc(sizeof(char)*(strlen(src.v.s) + 1)) );
		strcpy(dest->v.s, src.v.s);
	}

	return(NORMAL_RC);
}


/*
static long int
variable2int (PC_VARIABLE var)
{
	if(var.type == PC_INT) return(var.v.i);
	if(var.type == PC_DOUBLE) return((long int)(var.v.d));

	return((long int)LONG_MAX);
}
*/


static long double
variable2double (PC_VARIABLE var)
{
	if(var.type == PC_INT) return((long double)(var.v.i));
	if(var.type == PC_DOUBLE) return(var.v.d);

	return((long double)LDBL_MAX);
}


static RC
pc_operation (PC_VARIABLE var1, PC_VARIABLE var2, PC_VARIABLE *ans,
              PC_OPERATOR operator)
{
	int ii1;
	*ans = pc_init_variable();

	if( (var1.type == PC_ARRAY)&&(var2.type == PC_ARRAY) ){
		if(var1.array_size != var2.array_size) return(ARG_ERROR_RC);
		ans->type = PC_ARRAY;
		ans->v.array = mm_alloc(var1.array_size * sizeof(PC_VARIABLE));
		RC_NULL_CHK( ans->v.array );
		ans->array_size = ans->alloc_size = var1.array_size;

		for(ii1=0; ii1<var1.array_size; ii1++){
			RC_TRY( pc_operation(var1.v.array[ii1], var2.v.array[ii1],
			                     &(ans->v.array[ii1]), operator) );
		}

		return(NORMAL_RC);
	}

	if( (var1.type == PC_INT)&&(var2.type == PC_ARRAY) ){
		if(operator == PC_MUL){
			RC_TRY( copy_variable(var2, ans) );
			RC_TRY( mul_int_variable(var1.v.i, ans) );
		}else{
			return(ARG_ERROR_RC);
		}
		return(NORMAL_RC);
	}

	if( (var1.type == PC_DOUBLE)&&(var2.type == PC_ARRAY) ){
		if(operator == PC_MUL){
			RC_TRY( copy_variable(var2, ans) );
			RC_TRY( mul_double_variable(var1.v.d, ans) );
		}else{
			return(ARG_ERROR_RC);
		}
		return(NORMAL_RC);
	}

	if( (var1.type == PC_ARRAY)||(var2.type == PC_ARRAY) ){
		return(ARG_ERROR_RC);
	}

	if( (var1.type == PC_DOUBLE)||(var2.type == PC_DOUBLE) ){
		ans->type = PC_DOUBLE;
		if(operator == PC_PLUS){
			ans->v.d = variable2double(var1) + variable2double(var2);
		}else if(operator == PC_MINUS){
			ans->v.d = variable2double(var1) - variable2double(var2);
		}else if(operator == PC_MUL){
			ans->v.d = variable2double(var1) * variable2double(var2);
		}else if(operator == PC_DIV){
			ans->v.d = variable2double(var1) / variable2double(var2);
		}else if(operator == PC_POW){
			ans->v.d = (long double)pow(variable2double(var1),
			                            variable2double(var2));
		}else{
			return(ARG_ERROR_RC);
		}
		return(NORMAL_RC);
	}

	if( (var1.type == PC_INT)&&(var2.type == PC_INT) ){
		ans->type = PC_INT;
		if(operator == PC_PLUS){
			ans->v.i = var1.v.i + var2.v.i;
			if(SIGN(var1.v.i)*SIGN(var2.v.i) == 1){
				if(labs(var1.v.i) >= LONG_MAX - labs(var2.v.i)){
					RC_TRY( warning_printf(11001, FileLine) );
					ans->type = PC_DOUBLE;
					ans->v.d = (long double)var1.v.i + (long double)var2.v.i;
				}
			}
		}else if(operator == PC_MINUS){
			ans->v.i = var1.v.i - var2.v.i;
			if(SIGN(var1.v.i)*SIGN(var2.v.i) == -1){
				if(labs(var1.v.i) >= LONG_MAX - labs(var2.v.i)){
					RC_TRY( warning_printf(11001, FileLine) );
					ans->type = PC_DOUBLE;
					ans->v.d = (long double)var1.v.i - (long double)var2.v.i;
				}
			}
		}else if(operator == PC_MUL){
			ans->v = var2.v;
			RC_TRY( mul_int_variable(var1.v.i, ans) );
		}else if(operator == PC_DIV){
			ans->v.i = var1.v.i / var2.v.i;
		}else if(operator == PC_POW){
			double d_result = pow(fabs((double)(var1.v.i)),
			                      (double)(var2.v.i));
			if( (var2.v.i < 0)||(d_result >= (double)(LONG_MAX/2)) ){
				RC_TRY( warning_printf(11001, FileLine) );
				ans->type = PC_DOUBLE;
				ans->v.d = pow((double)(var1.v.i), (double)(var2.v.i));
			}else{
				ans->v.i = 1;
				for(ii1=0; ii1<(var2.v.i); ii1++){
					ans->v.i *= var1.v.i;
				}
			}
		}else if(operator == PC_REMAIN){
			ans->v.i = var1.v.i % var2.v.i;
		}else{
			return(ARG_ERROR_RC);
		}
		return(NORMAL_RC);
	}

	if( (var1.type == PC_STRING)&&(var2.type == PC_STRING) ){
		ans->type = PC_STRING;
		if(operator == PC_PLUS){
			ans->v.s = mm_alloc(sizeof(char)*(strlen(var1.v.s)
			                              + strlen(var2.v.s) + 1) );
			RC_NULL_CHK( ans->v.s );
			sprintf(ans->v.s, "%s%s", var1.v.s, var2.v.s);
		}else{
			return(ARG_ERROR_RC);
		}
		return(NORMAL_RC);
	}

	return(ARG_ERROR_RC);
}


static char *
get_name (const char *buf, int *ptr)
{
	int name_size = 0;
	int start_ptr = *ptr;
	char *ret = NULL;
	int ii1;
		
	while( (((char)'A' <= buf[*ptr])&&(buf[*ptr] <= (char)'Z'))
	     ||(((char)'a' <= buf[*ptr])&&(buf[*ptr] <= (char)'z'))
	     ||(((char)'0' <= buf[*ptr])&&(buf[*ptr] <= (char)'9'))
	     ||(buf[*ptr] == (char)'_') ){
		(*ptr)++;
		name_size++;
	}

	ret = (char *)mm_alloc(sizeof(char)*(name_size + 1));
	if(ret == NULL) return(NULL);

	for(ii1=0; ii1<name_size; ii1++){
		ret[ii1] = (char)toupper((int)(buf[start_ptr + ii1]));
	}
	ret[name_size] = (char)'\0';

	return(ret);
}


static RC
eval_fact (FILE *fp, char *buf, int buf_size, int *ptr,
           PC_VARIABLE **varray, int *varray_size, PC_VARIABLE *val)
{
	int sign = 1;
	int start_flag, stop_flag, point_flag;
	long double point_val, dtmp;
	int iexp;
	int ii1;
	int sign_flag;

	SubstFlag = 0;
	*val = pc_init_variable();

	/* 符号の処理 */
	sign_flag = 0;
	while(1){
		RC_TRY( next_token(fp, buf, buf_size, ptr) );
		if(buf[*ptr] == (char)'-'){
			sign *= -1;
			(*ptr)++;
			sign_flag = 1;
		}else if(buf[*ptr] == (char)'+'){
			(*ptr)++;
			sign_flag = 1;
		}else{
			break;
		}
	}
	RC_TRY( next_token(fp, buf, buf_size, ptr) );

	/* 文字列の処理 */
	if(buf[*ptr] == (char)'\"'){
		int count = 0;

		(*ptr)++;
		while(buf[(*ptr) + count] != '\"'){
			if( (buf[(*ptr) + count] == '\0')||(buf[(*ptr) + count] == '\n') ){
				RC_TRY( error_printf(11002, FileLine) );
				return(END_RC);
			}
			count++;
		}
		val->type = PC_STRING;
		RC_NULL_CHK( val->v.s = mm_alloc(sizeof(char)*(count + 1)) );

		count = 0;
		while(buf[*ptr] != '\"'){
			val->v.s[count] = buf[*ptr];
			count++;
			(*ptr)++;
		}
		val->v.s[count] = '\0';
		if(buf[*ptr] == '\"') (*ptr)++;

		RC_TRY( mul_int_variable((long int)sign, val) );
		return(NORMAL_RC);
	}

	/* 配列の処理 */
	if(buf[*ptr] == (char)'{'){
		(*ptr)++;
		val->type = PC_ARRAY;
		val->array_size = 0;
		val->alloc_size = 0;
		val->v.array = NULL;
		while(1){
			RC_TRY( next_token(fp, buf, buf_size, ptr) );
			if(buf[*ptr] == (char)'}'){
				(*ptr)++;
				break;
			}
			if( (buf[*ptr] == (char)'\0')||(buf[*ptr] == (char)'\n') ){
				RC_TRY( error_printf(11001, FileLine) );
				return(END_RC);
			}
			if(val->array_size >= val->alloc_size){
				val->v.array = mm_realloc(val->v.array,
				                 sizeof(PC_VARIABLE)*(4 + val->alloc_size) );
				RC_NULL_CHK( val->v.array );
				val->alloc_size += 4;
				if(val->array_size >= val->alloc_size) return(UNKNOWN_ERROR_RC);
			}
			RC_TRY( eval_add(fp, buf, buf_size, ptr, varray, varray_size,
			                 &(val->v.array[val->array_size]), NULL, NULL) );
			val->array_size++;
			RC_TRY( next_token(fp, buf, buf_size, ptr) );
			if(buf[*ptr] == (char)',') (*ptr)++;
		}

		RC_TRY( mul_int_variable((long int)sign, val) );
		return(NORMAL_RC);
	}

	/* 括弧の処理 */
	if( (buf[*ptr] == (char)'(')
	  ||(buf[*ptr] == (char)'[')
	  ||(buf[*ptr] == (char)'<') ){
		char open = buf[*ptr];

		(*ptr)++;
		RC_TRY( eval_add(fp, buf, buf_size, ptr,
		                 varray, varray_size, val, NULL, NULL) );

		RC_TRY( next_token(fp, buf, buf_size, ptr) );

		if( ((open == (char)'(')&&(buf[*ptr] != (char)')'))
		  ||((open == (char)'[')&&(buf[*ptr] != (char)']'))
		  ||((open == (char)'<')&&(buf[*ptr] != (char)'>')) ){
			RC_TRY( error_printf(11011, FileLine) );
			return(READ_ERROR_RC);
		}
		(*ptr)++;

		RC_TRY( mul_int_variable((long int)sign, val) );
		return(NORMAL_RC);
	}

	/* 変数の処理 */
	if( (((char)'A' <= buf[*ptr])&&(buf[*ptr] <= (char)'Z'))
	  ||(((char)'a' <= buf[*ptr])&&(buf[*ptr] <= (char)'z'))
	  ||(buf[*ptr] == (char)'_') ){
		char *vname;
		int index;
		RC_NULL_CHK( vname = get_name(buf, ptr) );

		RC_TRY( next_token(fp, buf, buf_size, ptr) );
		if(buf[*ptr] == (char)'='){ /* 代入演算子の処理 */
			(*ptr)++;
			RC_TRY( eval_add(fp, buf, buf_size, ptr,
			                 varray, varray_size, val, NULL, NULL) );
/*			if(val->type == PC_UNKNOWN){
				RC_TRY( warning_printf(11002, FileLine, vname) );
			}*/
			index = -1;
			for(ii1=0; ii1<(*varray_size); ii1++){
				if(case_strcmp((*varray)[ii1].name, vname) == 0){
					index = ii1;
					break;
				}
			}
			if(index < 0){
				*varray = mm_realloc(*varray,
				                     (*varray_size + 1)*sizeof(PC_VARIABLE));
				RC_NULL_CHK( *varray );
				index = *varray_size;
				(*varray)[index] = pc_init_variable();
				*varray_size += 1;
			}else{
				RC_TRY( pc_free_variable(&((*varray)[index])) );
			}
			RC_TRY( copy_variable(*val, &((*varray)[index])) );
			RC_NULL_CHK( (*varray)[index].name
			              = mm_alloc(sizeof(char)*(strlen(vname) + 1)) );
			strcpy((*varray)[index].name, vname);
			if(sign_flag == 0) SubstFlag = 1;
		}else{
			index = -1;
			for(ii1=0; ii1<(*varray_size); ii1++){
				if(case_strcmp((*varray)[ii1].name, vname) == 0){
					index = ii1;
					break;
				}
			}
			if(index < 0){
				RC_TRY( error_printf(11003, FileLine, vname) );
				RC_TRY( mm_free(vname) ); 
				return(SEEK_ERROR_RC);
			}
			RC_TRY( copy_variable((*varray)[index], val) );
		}
		RC_TRY( mm_free(vname) ); 
		RC_TRY( mul_int_variable((long int)sign, val) );

		return(NORMAL_RC);
	}


	(val->type) = PC_INT;
	(val->v.i) = 0;
	start_flag = 0;
	stop_flag = 0;
	point_val = 1.0L;
	point_flag = 0;
	while(1){
		switch((int)(buf[*ptr])){
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
				val->v.d += point_val * (long double)((int)(buf[*ptr]) - '0');
			}else{
				if( ((val->type) == PC_INT)
				  &&( (val->v.i) > (LONG_MAX - 9L)/10L) ){
					RC_TRY( warning_printf(11001, FileLine) );
					dtmp = (long double)(val->v.i);
					val->type = PC_DOUBLE;
					val->v.d = dtmp;
				}
				if((val->type) == PC_INT){
					(val->v.i) *= 10L;
					(val->v.i) += (long int)((int)(buf[*ptr]) - '0');
				}else{
					(val->v.d) *= 10.0L;
					(val->v.d) += (long double)((int)(buf[*ptr]) - '0');
				}
			}
			break;
		case '.':
			if(point_flag){
				stop_flag = 1;
				break;
			}
			if((val->type) == PC_INT){
				dtmp = (long double)(val->v.i);
				val->type = PC_DOUBLE;
				val->v.d = dtmp;
			}
			start_flag = 1;
			point_flag = 1;
			break;
		case 'E':
		case 'e':
			if((val->type) == PC_INT){
				dtmp = (long double)(val->v.i);
				val->type = PC_DOUBLE;
				val->v.d = dtmp;
			}
			(*ptr)++;
			RC_TRY( eval_int(buf, ptr, &iexp) );
			(val->v.d) *= (long double)(pow(10.0, (double)iexp));
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
		(val->type) = PC_UNKNOWN;
		(val->v.i) = LONG_MIN;
		return(READ_ERROR_RC);
	}
	
	RC_TRY( mul_int_variable((long int)sign, val) );

	return(NORMAL_RC);
}


static RC
eval_int (const char *buf, int *ptr, int *val)
{
	int sign = 1;
	int start_flag = 0;
	int stop_flag = 0;

	*val = 0;
	while(1){
		switch((int)(buf[*ptr])){
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
			*val += (int)(buf[*ptr]) - '0';
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


RC
eval_add (FILE *fp, char *buf, int buf_size, int *ptr,
          PC_VARIABLE **varray, int *varray_size, PC_VARIABLE *val,
          int *line, int *subst_flag)
{
	PC_VARIABLE v1, v2;

	if(line != NULL) FileLine = *line;
	*val = pc_init_variable();

	RC_TRY( next_token(fp, buf, buf_size, ptr) );
	if( (buf[*ptr] == (char)'\n')||(buf[*ptr] == (char)'\0') ){
		if(line != NULL) *line = FileLine;
		return(NORMAL_RC);
	}

	RC_TRY( eval_mul(fp, buf, buf_size, ptr, varray, varray_size, &v1) );

	while(1){
		RC_TRY( next_token(fp, buf, buf_size, ptr) );

		if(buf[*ptr] == (char)'+'){
			SubstFlag = 0;
			(*ptr)++;
			RC_TRY( eval_mul(fp, buf, buf_size, ptr,
			                 varray, varray_size, &v2) );
			RC_TRY( pc_operation(v1, v2, val, PC_PLUS) );
			RC_TRY( pc_free_variable(&v1) );
			RC_TRY( pc_free_variable(&v2) );
			v1 = *val;
		}else if(buf[*ptr] == (char)'-'){
			SubstFlag = 0;
			(*ptr)++;
			RC_TRY( eval_mul(fp, buf, buf_size, ptr,
			                 varray, varray_size, &v2) );
			RC_TRY( pc_operation(v1, v2, val, PC_MINUS) );
			RC_TRY( pc_free_variable(&v1) );
			RC_TRY( pc_free_variable(&v2) );
			v1 = *val;
		}else{
			break;
		}
	}
	*val = v1;
	if(line != NULL) *line = FileLine;
	if(subst_flag != NULL) *subst_flag = SubstFlag;

	return(NORMAL_RC);
}


static RC
eval_mul (FILE *fp, char *buf, int buf_size, int *ptr,
          PC_VARIABLE **varray, int *varray_size, PC_VARIABLE *val)
{
	PC_VARIABLE v1, v2;

	*val = pc_init_variable();
	RC_TRY( eval_expo(fp, buf, buf_size, ptr, varray, varray_size, &v1) );
	while(1){
		RC_TRY( next_token(fp, buf, buf_size, ptr) );

		if(buf[*ptr] == (char)'*'){
			SubstFlag = 0;
			(*ptr)++;
			RC_TRY( eval_expo(fp, buf, buf_size, ptr,
			                  varray, varray_size, &v2) );
			RC_TRY( pc_operation(v1, v2, val, PC_MUL) );
			RC_TRY( pc_free_variable(&v1) );
			RC_TRY( pc_free_variable(&v2) );
			v1 = *val;
		}else if(buf[*ptr] == (char)'/'){
			SubstFlag = 0;
			(*ptr)++;
			RC_TRY( eval_expo(fp, buf, buf_size, ptr,
			                  varray, varray_size, &v2) );
			RC_TRY( pc_operation(v1, v2, val, PC_DIV) );
			RC_TRY( pc_free_variable(&v1) );
			RC_TRY( pc_free_variable(&v2) );
			v1 = *val;
		}else if(buf[*ptr] == (char)'%'){
			SubstFlag = 0;
			(*ptr)++;
			RC_TRY( eval_expo(fp, buf, buf_size, ptr,
			                  varray, varray_size, &v2) );
			RC_TRY( pc_operation(v1, v2, val, PC_REMAIN) );
			RC_TRY( pc_free_variable(&v1) );
			RC_TRY( pc_free_variable(&v2) );
			v1 = *val;
		}else{
			break;
		}
	}
	*val = v1;

	return(NORMAL_RC);
}


static RC
eval_expo (FILE *fp, char *buf, int buf_size, int *ptr,
           PC_VARIABLE **varray, int *varray_size, PC_VARIABLE *val)
{
	PC_VARIABLE v1, v2;

	*val = pc_init_variable();
	RC_TRY( eval_fact(fp, buf, buf_size, ptr, varray, varray_size, val) );

	RC_TRY( next_token(fp, buf, buf_size, ptr) );

	if((int)(buf[*ptr]) == '^'){
		(*ptr)++;
	}else if( ((int)(buf[*ptr]) == '*')&&((int)(buf[(*ptr)+1]) == '*') ){
		(*ptr) += 2;
	}else{
		return(NORMAL_RC);
	}

	SubstFlag = 0;
	RC_TRY( eval_expo(fp, buf, buf_size, ptr, varray, varray_size, &v1) );
	RC_TRY( pc_operation(*val, v1, &v2, PC_POW) );
	RC_TRY( pc_free_variable(val) );
	RC_TRY( pc_free_variable(&v1) );
	*val = v2;

	return(NORMAL_RC);

}


RC
pc_print_variable_array (FILE *fp, PC_VARIABLE *varray, int varray_size)
{
	int ii1;

	for(ii1=0; ii1<varray_size; ii1++){
		RC_TRY( print_variable(fp, varray[ii1]) );
		fprintf(fp, "\n");
	}

	return(NORMAL_RC);
}


static RC
print_variable (FILE *fp, PC_VARIABLE var)
{
	int ii1;

	if(var.name) fprintf(fp, "%s = ", var.name);
	if(var.type == PC_INT){
		fprintf(fp, "%ld", var.v.i);
	}else if(var.type == PC_DOUBLE){
		fprintf(fp, "%13.5e", (double)(var.v.d));
	}else if(var.type == PC_STRING){
		fprintf(fp, "\"%s\"", var.v.s);
	}else if(var.type == PC_ARRAY){
		fprintf(fp, "{");
		for(ii1=0; ii1<var.array_size; ii1++){
			RC_TRY( print_variable(fp, var.v.array[ii1]) );
			if(ii1 < var.array_size - 1) fprintf(fp, ", ");
		}
		fprintf(fp, "}");
	}else{
		fprintf(fp, "???");
	}

	return(NORMAL_RC);
}


static RC
next_token (FILE *fp, char *buf, int buf_size, int *ptr)
{
	while(1){
		while( (buf[*ptr] == (char)' ')
		     ||(buf[*ptr] == (char)'\t') ) (*ptr)++;

		if(buf[*ptr] == (char)'\\'){
			if(fp == NULL) return(END_RC);
			if(fgets(buf, buf_size, fp) == NULL){
				if(feof(fp)){
					return(END_RC);
				}else{
					return(READ_ERROR_RC);
				}
			}
			clean_cr_string(buf);
			FileLine++;
			*ptr = 0;
			continue;
		}

		break;
	}

	return(NORMAL_RC);
}


int
pc_search_varray (const char vname[],
                  int varray_size, const PC_VARIABLE varray[])
{
	int ii1;
	int index = -1;

	for(ii1=0; ii1<varray_size; ii1++){
		if(case_strcmp(varray[ii1].name, vname) == 0){
			index = ii1;
			break;
		}
	}

	return(index);
}


RC
pc_varray2str (const char vname[], int varray_size,
               const PC_VARIABLE varray[], char **v)
{
	int index = pc_search_varray(vname, varray_size, varray);
	int len;

	if((index < 0)||(index >= varray_size)){
		*v = NULL;
		return(NORMAL_RC);
	}
	if(varray[index].type == PC_UNKNOWN){
		*v = NULL;
		return(NORMAL_RC);
	}

	if(varray[index].type != PC_STRING){
		RC_TRY( error_printf(11010, vname) );
		return(CONVERT_ERROR_RC);
	}

	len = strlen(varray[index].v.s) + 1;
	if(len <= 1){
		*v = NULL;
		return(NORMAL_RC);
	}
	RC_NULL_CHK( *v = mm_alloc(len*sizeof(char)) );
	strcpy(*v, varray[index].v.s);

	return(NORMAL_RC);
}


RC
pc_varray2int (const char vname[], int varray_size, const PC_VARIABLE varray[],
               int *v, int v_min, int v_max, int v_default)
{
	int index = pc_search_varray(vname, varray_size, varray);

	RC_TRY( pc_variable2int(index, varray_size, varray,
	                        v, v_min, v_max, v_default, vname) );

	return(NORMAL_RC);
}


RC
pc_variable2int (int index, int varray_size, const PC_VARIABLE varray[],
                 int *v, int v_min, int v_max, int v_default,
                 const char vname[])
{
	if((index < 0)||(index >= varray_size)){
		*v = v_default;
		return(NORMAL_RC);
	}
	if(varray[index].type == PC_UNKNOWN){
		*v = v_default;
		return(NORMAL_RC);
	}

	if(varray[index].type != PC_INT){
		RC_TRY( error_printf(11004, vname) );
		return(CONVERT_ERROR_RC);
	}

	if(varray[index].v.i < v_min){
		*v = v_min;
		RC_TRY( warning_printf(11005, vname, v_min, v_max, *v) );
	}else if(varray[index].v.i > v_max){
		*v = v_max;
		RC_TRY( warning_printf(11005, vname, v_min, v_max, *v) );
	}else{
		*v = (int)(varray[index].v.i);
	}

	return(NORMAL_RC);
}


RC
pc_varray2double (const char vname[], int varray_size,
                  const PC_VARIABLE varray[],
                  double *v, double v_min, double v_max, double v_default)
{
	int index = pc_search_varray(vname, varray_size, varray);

	RC_TRY( pc_variable2double(index, varray_size, varray,
	                           v, v_min, v_max, v_default, vname) );

	return(NORMAL_RC);
}


RC
pc_variable2double (int index, int varray_size, const PC_VARIABLE varray[],
                    double *v, double v_min, double v_max, double v_default,
                    const char vname[])
{
	if((index < 0)||(index >= varray_size)){
		*v = v_default;
		return(NORMAL_RC);
	}
	if(varray[index].type == PC_UNKNOWN){
		*v = v_default;
		return(NORMAL_RC);
	}

	if(varray[index].type == PC_INT){
		*v = (double)(varray[index].v.i);
	}else if(varray[index].type == PC_DOUBLE){
		*v = (varray[index].v.d);
	}else{
		RC_TRY( error_printf(11005, vname) );
		return(CONVERT_ERROR_RC);
	}
	if(*v < v_min){
		*v = v_min;
		RC_TRY( warning_printf(11006, vname, v_min, v_max, *v) );
	}else if(*v > v_max){
		*v = v_max;
		RC_TRY( warning_printf(11006, vname, v_min, v_max, *v) );
	}

	return(NORMAL_RC);
}


RC
pc_variable2int_array (int index, int varray_size, const PC_VARIABLE varray[],
                       int v[], int *v_size, int v_size_max,
                       int v_min, int v_max, int v_default, const char vname[])
{
	int ii1;

	if((index < 0)||(index >= varray_size)){
		*v_size = 0;
		return(NORMAL_RC);
	}
	if(varray[index].type == PC_UNKNOWN){
		*v_size = 0;
		return(NORMAL_RC);
	}

	if(varray[index].type == PC_INT){
		if(v_size_max < 1) return(OVERRUN_ERROR_RC);
		RC_TRY( pc_variable2int(index, varray_size, varray, &(v[0]),
		                        v_min, v_max, v_default, vname) );
		*v_size = 1;
		return(NORMAL_RC);
	}

	if(varray[index].type != PC_ARRAY){
		RC_TRY( error_printf(11004, vname) );
		return(CONVERT_ERROR_RC);
	}
	if(varray[index].array_size > v_size_max){
		return(OVERRUN_ERROR_RC);
	}

	*v_size = varray[index].array_size;
	for(ii1=0; ii1<varray[index].array_size; ii1++){
		RC_TRY( pc_variable2int(ii1, varray[index].array_size,
		        varray[index].v.array, &(v[ii1]), v_min, v_max,
		        v_default, vname) );
	}

	return(NORMAL_RC);
}


RC
pc_variable2int_array_var (int index, int varray_size,
                           const PC_VARIABLE varray[], int **v, int *v_size,
                           int v_min, int v_max, int v_default,
                           const char vname[])
{
	int ii1;

	*v_size = 0;
	*v = NULL;
	if((index < 0)||(index >= varray_size)) return(NORMAL_RC);
	if(varray[index].type == PC_UNKNOWN) return(NORMAL_RC);

	if(varray[index].type == PC_INT){
		*v_size = 1;
		RC_NULL_CHK( *v = mm_alloc(1*sizeof(int)) );
		RC_TRY( pc_variable2int(index, varray_size, varray, &((*v)[0]),
		                        v_min, v_max, v_default, vname) );
		return(NORMAL_RC);
	}

	if(varray[index].type != PC_ARRAY){
		RC_TRY( error_printf(11004, vname) );
		return(CONVERT_ERROR_RC);
	}

	*v_size = varray[index].array_size;
	if(*v_size > 0) RC_NULL_CHK( *v = mm_alloc(*v_size * sizeof(int)) );
	for(ii1=0; ii1<(*v_size); ii1++){
		RC_TRY( pc_variable2int(ii1, varray[index].array_size,
		        varray[index].v.array, &((*v)[ii1]), v_min, v_max,
		        v_default, vname) );
	}

	return(NORMAL_RC);
}

