/*********************************************************************
 * nst_component.h
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Thanks to following contributors.
 *  <Takaaki NAGATANI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: nst_component.h 300 2004-06-09 11:00:54Z sasaoka $ */

#ifndef NST_COMPONENT_H
#define NST_COMPONENT_H

#include <limits.h>
#include <float.h>
#include "fem_struct.h"
#include "math_utl.h"
#include "rc.h"
#include "string_utl.h"

#define NSDEFAULT     (INT_MIN)
#define FNSDEFAULT    (-(DBL_MAX/2.0))
#define CR_FNSDEFAULT (FNSDEFAULT/2.0)

#define NST_PROP_PSHELL  (1)
#define NST_PROP_PSOLID  (2)
#define NST_PROP_PBAR    (3)
#define NST_PROP_PBARL   (4)
#define NST_PROP_PELAS   (5)

#define NST_LC_2R        (1)

#define NST_SPRING_CELAS1 (1)
#define NST_SPRING_CELAS2 (2)

#define NST_RBAR_RBAR     (1)
#define NST_RBAR_RBE2     (2)

#define NST_MASS_CONM2    (1)

#define F06_KEY_DISP          "D I S P L A C E M E N T"
#define F06_KEY_EIGENVECTOR   "E I G E N V E C T O R"
#define F06_KEY_EIGENVALUE    "R E A L   E I G E N V A L U E S"
#define F06_KEY_STRESS        "S T R E S S E S"
#define F06_KEY_STRAIN        "S T R A I N S"

#define F06_SW_NULL      (0)
#define F06_SW_SUBCASE   (1)
#define F06_SW_LOADSTEP  (2)
#define F06_SW_MODE      (4)

#define F06_SS_SIMPLE    (0)
#define F06_SS_CORNER    (1)
#define F06_SS_FIBRE     (2)

#define PCH_SW_NULL      (0)
#define PCH_SW_SUBCASE   (1)

#define PCH_DATA_UNKNOWN (0)
#define PCH_DATA_DISP    (1)
#define PCH_DATA_STRESS  (2)
#define PCH_DATA_STRAIN  (3)

#define NST_FILE_BUFS (4)
typedef struct {
	FILE *parent;
	FILE *child;
	char *bufs[NST_FILE_BUFS];
	int r_point;
	int w_point;
} NST_FILE;


typedef struct {
	char *key;
	int key_size;
} NST_BULK_KEY;


typedef struct {
	char name[MAX_IDC_S_SIZE+1];
	IDC val_1;
	IDC val_2;
} NST_PARAM;


typedef struct {
	STRING_ARRAY exective;
	int size;
	NST_PARAM *array;
} NST_PARAM_ARRAY;


typedef struct{
	int ninc;
	char conv[4];
	double eps[3];
} NST_NL_PARAM;


/* nst_blk_llio.c */
RC init_nst_file(FILE *fp, NST_FILE *nst_file);
RC term_nst_file(NST_FILE *nst_file);

RC nst_allocate_tok(char **toks, int size);
RC nst_free_tok(char **toks, int size);
void nst_print_tok(char **tosk, int size);

RC nst_seek_begin_bulk(NST_FILE *nst_file);
RC nst_seek_cend(NST_FILE *nst_file);
RC nst_gets(char **s, NST_FILE *nst_file);
RC nst_ungets(NST_FILE *nst_file);
int nst_chk_begin_bulk(const char *buf);
int nst_chk_cend(const char *buf);
RC nst_blk_get_tok(NST_FILE *nst_file,
                   int num_keys, NST_BULK_KEY *keys, int *hit_key,
                   int max_toks, char **toks, int *num_toks);

double nst_tok2d(int tok_size, int tok_num, char **toks);
int nst_tok2i(int tok_size, int tok_num, char **toks);
int nst_tok_strncmp(int tok_size, int tok_num, char **toks, 
                    const char *str, int size);
RC nst_get_param(const char name[], const char format[],
                 NST_PARAM_ARRAY param, void *ptr1, void *ptr2);

RC nst_print(FILE *stream, const char *format, const IDC *data);


/* nst_blk_component.c */
RC nst_input_node(FILE *fp, NODE_ARRAY *node);
RC nst_input_element(FILE *fp, ELEMENT_ARRAY *element);
RC nst_input_rigid_element(FILE *fp, RIGID_ELEMENT_ARRAY *element);
RC nst_input_material_prop(FILE *fp, MATERIAL_PROP_ARRAY *material);
RC nst_input_physical_prop(FILE *fp, PHYSICAL_PROP_ARRAY *phys);
RC nst_input_force(FILE *fp, BC_ARRAY *force,
                   DEFAULT_BC_ARRAY *b_force);
RC nst_input_restraint(FILE *fp, BC_ARRAY *rest);
RC nst_input_temperature(FILE *fp, BC_ARRAY *temp, 
                         ELEM_BC_ARRAY *elem_temp,
                         DEFAULT_BC_ARRAY *def_temp);
RC nst_input_param(FILE *fp, NST_PARAM_ARRAY *param);
RC nst_input_bc_set(FILE *fp, BC_SET_ARRAY *bc_set);
RC nst_input_local_coord(FILE *fp, LOCAL_COORD_ARRAY *coord);
RC nst_input_NL_param(FILE *fp, NST_NL_PARAM *param);
RC nst_input_pressure(FILE *fp, ELEM_BC_ARRAY *press,
                      ELEMENT_ARRAY element);
RC nst_input_node_matrix(FILE *fp, char *name, NODE_MATRIX_ARRAY *nmat);

RC nst_output_node(FILE *fp, NODE_ARRAY node);
RC nst_output_element(FILE *fp, ELEMENT_ARRAY element);
RC nst_output_rigid_element(FILE *fp, RIGID_ELEMENT_ARRAY element);
RC nst_output_material_prop(FILE *fp, MATERIAL_PROP_ARRAY material);
RC nst_output_physical_prop(FILE *fp, PHYSICAL_PROP_ARRAY physical);
RC nst_output_force(FILE *fp, BC_ARRAY force,
                    DEFAULT_BC_ARRAY b_force);
RC nst_output_restraint(FILE *fp, BC_ARRAY rest);
RC nst_output_temperature(FILE *fp, BC_ARRAY temp, 
                          ELEM_BC_ARRAY elem_temp,
                          DEFAULT_BC_ARRAY def_temp);
RC nst_output_param(FILE *fp, NST_PARAM_ARRAY param);
RC nst_output_bc_set(FILE *fp, BC_SET_ARRAY bc_set);
RC nst_output_local_coord(FILE *fp, LOCAL_COORD_ARRAY coord);
RC nst_output_pressure(FILE *fp, ELEM_BC_ARRAY press,
                       ELEMENT_ARRAY element);
RC nst_output_node_matrix(FILE *fp, char *name,
                          NODE_MATRIX_ARRAY nmat);

RC nst_write_exective_control(FILE *fp, NST_PARAM_ARRAY param);
RC nst_write_case_control(FILE *fp, BC_SET_ARRAY bc_set);


/* nst_f06_component.c */
RC nst_input_disp(FILE *fp, DISP_ARRAY *disp);
RC nst_input_eigenvector_value(FILE *fp,
                               DISP_ARRAY *evect, double value);
RC nst_input_eigenvector_mode(FILE *fp, DISP_ARRAY *evect, int mode);
RC nst_input_eigen_value(FILE *fp, int *array_size, double **value_array);

RC read_f06disp(FILE *fp, int subcase, double load_step, int mode, int f06_sw,
                DISP_ARRAY *disp, const char *key);
RC read_f06eigen_value(FILE *fp, int subcase, int f06_sw, int *array_size,
                       double **value_array);
RC read_f06stress(FILE *fp, int subcase, double load_step, int mode,
                  int f06_sw, STRESS_ARRAY *stress, int f06_ss);
RC read_f06strain(FILE *fp, int subcase, double load_step, int mode,
                  int f06_sw, STRAIN_ARRAY *strain, int f06_ss);


/* nst_pch_component.c */
RC read_pch_disp(FILE *fp, int subcase, int pch_sw, DISP_ARRAY *disp);
RC read_pch_stress(FILE *fp, int subcase, int pch_sw,
                   STRESS_ARRAY *stress);
RC read_pch_strain(FILE *fp, int subcase, int pch_sw,
                   STRAIN_ARRAY *strain);


#endif  /* NST_COMPONENT_H */

