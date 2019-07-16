/*********************************************************************
 * fem_struct.h
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Tmoyuki TSUBATA>
 *
 * Thanks to following contributors.
 *  <Yasuo ASAGA> <Masanobu SUYAMA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: fem_struct.h,v 1.30 2003/07/22 11:44:01 adachi Exp $ */

#ifndef FEM_STRUCT_H
#define FEM_STRUCT_H

#include <stdio.h>
#include "rc.h"
#include "string_utl.h"
#include "math_utl.h"


/* node */
#define FEM_NODE_INFO    4
typedef struct {
	int label;           /* label of node(if < 0 then unused) */
	TRANSLATION3D p;     /* x,y,z coordinate value */
	int renum;           /* for renumbering */
	int i_info[FEM_NODE_INFO];     /* the other informations(int) */
	double d_info[FEM_NODE_INFO];  /* the other informations(double) */
} FEM_NODE;


/* type of element */
typedef enum {
	ELEM_RBAR,
	ELEM_RBODY,
	ELEM_RFORCE,
	ELEM_MPC,
	ELEM_SPRING,
	ELEM_MASS,
	ELEM_BEAM1,
	ELEM_BEAM2,
	ELEM_SHELL1,
	ELEM_SHELL2,
	ELEM_LINE1,
	ELEM_LINE2,
	ELEM_TRI1,
	ELEM_TRI2,
	ELEM_QUAD1,
	ELEM_QUAD2,
	ELEM_TETRA1,
	ELEM_TETRA2,
	ELEM_PENTA1,
	ELEM_PENTA2,
	ELEM_HEXA1,
	ELEM_HEXA2
} FEM_ELEM_TYPE;


/* element */
#define FEM_MAX_NODE     20
#define FEM_ELEMENT_INFO  8
typedef struct {
	int label;                     /* label of element(if < 0 then unused) */
	FEM_ELEM_TYPE type;            /* element type */
	int material;                  /* material property label */
	int physical;                  /* physical property label */
	int node_num;                  /* number of node */
	int node[FEM_MAX_NODE];        /* node labels */
	double volume;                 /* volume or area of element */
	int i_info[FEM_ELEMENT_INFO];  /* another information (int) */
	double d_info[FEM_ELEMENT_INFO]; /* another information (double) */
} FEM_ELEMENT;


/* rigid or mpc element */
#define FEM_RIGID_ELEMENT_INFO    10
typedef struct {
	int label;             /* label of element(if < 0 then unused) */
	int set_id;            /* ID of bc set (if < 0 then permanent element) */
	FEM_ELEM_TYPE type;    /* element type */
	int physical;          /* physical property label */
	int node_num;          /* number of node */
	int *node;             /* node labels */
	int weight_num;        /* size of weight[] */
	double *weight;        /* weight */
	int dof_num;           /* size of dof[] */
	int *dof;              /* DOF flag (0:x, 1:y, 2:z, 3:yz, 4:zx, 5:xy) */
	int i_info[FEM_RIGID_ELEMENT_INFO];    /* information (int) */
	double d_info[FEM_RIGID_ELEMENT_INFO]; /* information (double) */
} FEM_RIGID_ELEMENT;


/* type of material */
typedef enum {
	MAT_ISOTROPIC,
	MAT_ORTHOTROPIC,
	MAT_ANISOTROPIC
} FEM_MAT_TYPE;


/* type of analysis */
typedef enum {
	ANA_3D,
	ANA_PLANE_STRESS,
	ANA_PLANE_STRAIN,
	ANA_AXISYMMETRY,
	ANA_PLANE_BEND,
	ANA_PLANE_SHEAR
} FEM_ANALYSIS_TYPE;


/* material property */
#define FEM_MATERIAL_INFO 8
typedef struct {
	int label;              /* label of material (if < 0 then unused) */
	FEM_MAT_TYPE mat_type;  /* material type */
	FEM_ANALYSIS_TYPE ana_type;  /* analysis type */

	/* for mat_type == MAT_ISOTROPIC */
	double E;               /* Youg's modulus */
	double G;               /* Shear modulus */
	double nu;              /* Poisson's ratio */ 
	double alpha;           /* coefficient of thermal expansion */
	double K;               /* 熱伝導率 */

	/* for mat_type == MAT_ORTHOTROPIC */
	TRANSLATION3D E_vect;   /* Young's modulus(x,y,z) */
	ROTATION3D G_vect;      /* Shear mdulus(yz,zx,xy) */
	ROTATION3D nu_vect;     /* Poisson's ratio(yz,zx,xy) */
	TRANS_ROTATION3D alpha_vect; /* coefficient of thermal expansion */

	/* for mat_type == MAT_UNISOTROPIC */
	double D_matrix[6][6];  /* stress - strain matrix */

	double C[3][3][3][3];   /* stiffness tensor */
	double rho;             /* density */
	double tref;            /* refference temperature */
	int i_info[FEM_MATERIAL_INFO];    /* another information (int) */
	double d_info[FEM_MATERIAL_INFO]; /* another information (double) */
} FEM_MATERIAL_PROP;


typedef enum{
	SECTION_ROD,
	SECTION_TUBE,
	SECTION_BAR,
	SECTION_BOX,
	SECTION_BOX1,
	SECTION_CHAN,
	SECTION_CHAN1,
	SECTION_CHAN2,
	SECTION_I,
	SECTION_I1,
	SECTION_T,
	SECTION_T1,
	SECTION_T2,
	SECTION_H,
	SECTION_Z,
	SECTION_CROSS,
	SECTION_HEXA,
	SECTION_HAT,
	SECTION_UNKNOWN
} FEM_SECTION_TYPE;


/* physical property */
#define FEM_PHYSICAL_INFO 24
typedef struct {
	int label;           /* label of physical property (if < 0 then unused) */
	FEM_SECTION_TYPE section;         /* cross section type of beam */
	int i_info[FEM_PHYSICAL_INFO];    /* another information (int) */
	double d_info[FEM_PHYSICAL_INFO]; /* another information (double) */
} FEM_PHYSICAL_PROP;


/* type of b.c. */
typedef enum {
	BC_FREE,
	BC_FIX,
	BC_FIX_NL,
	BC_MPC,
	BC_FORCE,
	BC_FORCE_NL,
	BC_VELO,
	BC_PRESS,        /* 一様圧力 */
	BC_PRESS_D1,     /* １次分布圧力 */
	BC_PRESS_D2,     /* ２次分布圧力 */
	BC_TRACTION,     /* 一様分布力 */
	BC_TRACTION_D1,  /* １次分布力 */
	BC_TRACTION_D2,  /* ２次分布力 */
	BC_TEMP,
	BC_TEMP_GR,
	BC_TEMP_NL,
	BC_K_DIAGONAL
} FEM_BC_TYPE;


/* type of b.c. on each DOF */
typedef struct {
	FEM_BC_TYPE x;
	FEM_BC_TYPE y;
	FEM_BC_TYPE z;
	FEM_BC_TYPE yz;
	FEM_BC_TYPE zx;
	FEM_BC_TYPE xy;
} FEM_BC_TYPE3D;


/* b.c. on node */
#define FEM_BC_INFO    4
#define FEM_BC_SCALAR  2
typedef struct {
	int node;                     /* label of node(if < 0 then unused) */
	int set_id;                   /* ID of bc set */
	FEM_BC_TYPE3D v_type;         /* type of bc on each DOF */
	TRANS_ROTATION3D v;           /* bc on each DOF */
	I_TRANS_ROTATION3D v_nl;      /* for nonlinear */
	FEM_BC_TYPE s_type[FEM_BC_SCALAR];  /* type of scalar bc */
	double s[FEM_BC_SCALAR];      /* scalar values */
	int s_nl[FEM_BC_SCALAR];      /* for nonlinear */
	int i_info[FEM_BC_INFO];      /* another information (int) */
	double d_info[FEM_BC_INFO];   /* another information (double) */
} FEM_BC;


/* b.c. on element */
#define FEM_ELEM_BC_INFO    4
#define FEM_ELEM_BC_SCALAR  2
typedef struct {
	int element;                  /* label of element(if < 0 then unused) */
	int set_id;                   /* ID of bc set */
	FEM_BC_TYPE type;             /* type of B.C. */
	int face_id;                  /* ID of face */
	union {
		double s;                 /* uniform scalar value */
		TRANS_ROTATION3D v;       /* uniform vector value */
		double ns[FEM_MAX_NODE];  /* scalar value defined on each node */
	} val;
	TRANSLATION3D dir;
	int i_info[FEM_ELEM_BC_INFO];      /* another information (int) */
	double d_info[FEM_ELEM_BC_INFO];   /* another information (double) */
} FEM_ELEM_BC;


/* displacement or velocity on node */
#define FEM_DISP_VELO_INFO    4
typedef struct {
	int node;                     /* label of node(if < 0 then unused) */
	TRANS_ROTATION3D v;           /* displacement or velocity */
	double s;                     /* pressure */
	int i_info[FEM_DISP_VELO_INFO];    /* another information (int) */
	double d_info[FEM_DISP_VELO_INFO]; /* another information (double) */
} FEM_DISP_VELO;


/* stress or strain on node of element */
#define FEM_STRESS_STRAIN_INFO 4
typedef struct {
	int element;     /* label of element(if < 0 then unused) */
	int node;        /* label of node (if < 0 then center of element) */
	TRANS_ROTATION3D v;   /* stress(or strain) vector x,y,z,yz,zx,xy */
	/*double tensor[3][3];*/  /* stress(or strain) tensor */
	/*double mises;*/         /* von mises */
	/*double principal[3];*/  /* principal stress(or strain) */
	/*double p_dir[3][3];*/   /* direction of principal stress(or strain) */
	int i_info[FEM_STRESS_STRAIN_INFO];    /* another information (int) */
	double d_info[FEM_STRESS_STRAIN_INFO]; /* another information (double) */
} FEM_STRESS_STRAIN;


/* sensitivity on node of element */
typedef struct {
	int element;     /* label of element(if < 0 then unused) */
	int node;        /* label of node (if < 0 then center of element) */
	TRANS_ROTATION3D v; /* sensitivity(vector) */
	double s;        /* sensitivity(scalar) */
} FEM_SENSITIVITY;


/* default temperature and gravity force */
#define FEM_DEFAULT_BC_INFO 4
typedef struct {
	int set_id;         /* ID of bc set (if < 0 then unused) */
	double temp;        /* default temperature */
	TRANSLATION3D grav;              /* gravity force vector */
	int i_info[FEM_DEFAULT_BC_INFO];    /* another information (int) */
	double d_info[FEM_DEFAULT_BC_INFO]; /* another information (double) */
} FEM_DEFAULT_BC;


/* b.c. case set */
#define FEM_MAX_BC_SET 32
#define FEM_BC_SET_INFO 4
typedef struct {
	int label;
	int num_fix;
	double scale_fix[FEM_MAX_BC_SET];
	int id_fix[FEM_MAX_BC_SET];
	int num_mpc;
	double scale_mpc[FEM_MAX_BC_SET];
	int id_mpc[FEM_MAX_BC_SET];
	int num_force;
	double scale_force[FEM_MAX_BC_SET];
	int id_force[FEM_MAX_BC_SET];
	int num_temp;
	double scale_temp[FEM_MAX_BC_SET];
	int id_temp[FEM_MAX_BC_SET];
	int i_info[FEM_BC_SET_INFO];
	double d_info[FEM_BC_SET_INFO];
	STRING_ARRAY s_info;
} FEM_BC_SET;


/* pressure of sound field */
#define FEM_SOUND_PRESSURE_INFO 4
typedef struct{
	int node;    /* label of node (if < 0 then unused) */
	double re;   /* real part */
	double im;   /* imaginary part */
	int i_info[FEM_SOUND_PRESSURE_INFO];
	double d_info[FEM_SOUND_PRESSURE_INFO];
} FEM_SOUND_PRESSURE;


#define FEM_LOAD_CURVE_MAX_INC    1024
#define FEM_LOAD_CURVE_MAX_APPROX   50
typedef struct {
	int label;
	double max_time;            /* load(disp) maximum time  */
	double max_value;           /* load(disp) maximum value */
	int total_step;             /* number of load(disp) step */
	double inc_value[FEM_LOAD_CURVE_MAX_INC];
	                            /* increment load(disp)value at every step */
	int approx_size;            /* number of load(disp) division point */
	double approx_time[FEM_LOAD_CURVE_MAX_APPROX];
	double approx_value[FEM_LOAD_CURVE_MAX_APPROX];
} FEM_LOAD_CURVE;


#define FEM_LOCAL_COORD_INFO 16
typedef struct {
	int label;
	TRANSLATION3D origin;   /* 原点 */
	TRANSLATION3D d_cos[3]; /* 方向余弦 */
	int i_info[FEM_LOCAL_COORD_INFO];
	double d_info[FEM_LOCAL_COORD_INFO];
} FEM_LOCAL_COORD;


/* matrix of element */
typedef struct{
	int element;      /* label of element (if < 0 then unused) */
	int size;         /* total size */
	int dim;          /* dof on each node */
	int *index;       /* index on global matrix */
	double **matrix;  /* element matrix */
} FEM_ELEM_MATRIX;


/* node -> element lookup table */
typedef struct{
	int node;         /* label of node (if < 0 then unused) */
	int size;         /* size of lookup table */
	int alloc_size;   /* allocated size of lookup table */
	int *element;     /* lookup table of element label */
	int r_size;       /* size of lookup table for rigid element */
	int r_alloc_size; /* allocated size of lookup table for rigid element */
	int *r_element;   /* lookup table of element label for rigid element */
} FEM_ELEM_LOOKUP;


typedef struct{
	int node;
	int shift;               /* 0: 並進自由度， 1: 回転自由度 */
	int size;
	int alloc_size;
	int *col_node;
	int *col_shift;          /* 0: 並進自由度， 1: 回転自由度 */
	double (*matrix)[3][3];
} FEM_NODE_MATRIX;


/* FEM data type */
typedef enum {
	FEM_ANSYS,
	FEM_IDEAS,
	FEM_NASTRAN,
	FEM_SYSNOISE,
	FEM_DXF,
	FEM_STL,
	FEM_QUINT,
	FEM_GENERATED_SURFACE,
	FEM_K_SOL,
	FEM_ADV,
	FEM_NIKE,
	FEM_UNKNOWN
} FEM_DATA_TYPE;


/* Sort flag */
typedef enum {
	ARRAY_SORTED,
	ARRAY_UNSORTED
} FEM_ARRAY_SORT_FLAG;


/* renumber flag */
typedef enum {
	RENUMBERED,
	UNRENUMBERED
} FEM_RENUMBER_FLAG;


/* array of node */
typedef struct {
	int size;
	int alloc_size;
	FEM_DATA_TYPE source;
	FEM_NODE *array;
	FEM_RENUMBER_FLAG renum_flag;
	FEM_ARRAY_SORT_FLAG sort_flag;
} FEM_NODE_ARRAY;


/* array of element */
typedef struct {
	int size;
	int alloc_size;
	FEM_DATA_TYPE source;
	FEM_ELEMENT *array;
	FEM_ARRAY_SORT_FLAG sort_flag;
} FEM_ELEMENT_ARRAY;


/* array of rigid element */
typedef struct {
	int size;
	int alloc_size;
	FEM_DATA_TYPE source;
	FEM_RIGID_ELEMENT *array;
	FEM_ARRAY_SORT_FLAG sort_flag;
} FEM_RIGID_ELEMENT_ARRAY;


/* array of material */
typedef struct {
	int size;
	int alloc_size;
	FEM_DATA_TYPE source;
	FEM_MATERIAL_PROP *array;
	FEM_ARRAY_SORT_FLAG sort_flag;
} FEM_MATERIAL_PROP_ARRAY;


/* array of physical property */
typedef struct {
	int size;
	int alloc_size;
	FEM_DATA_TYPE source;
	FEM_PHYSICAL_PROP *array;
	FEM_ARRAY_SORT_FLAG sort_flag;
} FEM_PHYSICAL_PROP_ARRAY;


/* array of b.c. */
typedef struct {
	int size;
	int alloc_size;
	FEM_DATA_TYPE source;
	FEM_BC *array;
	FEM_ARRAY_SORT_FLAG sort_flag;
} FEM_BC_ARRAY;


/* array of b.c. */
typedef struct {
	int size;
	int alloc_size;
	FEM_DATA_TYPE source;
	FEM_ELEM_BC *array;
	FEM_ARRAY_SORT_FLAG sort_flag;
} FEM_ELEM_BC_ARRAY;


/* array of displacement */
typedef struct{
	int size;
	int alloc_size;
	FEM_DATA_TYPE source;
	FEM_DISP_VELO *array;
	FEM_ARRAY_SORT_FLAG sort_flag;
} FEM_DISP_ARRAY;


/* array of velocity */
typedef struct{
	int size;
	int alloc_size;
	FEM_DATA_TYPE source;
	FEM_DISP_VELO *array;
	FEM_ARRAY_SORT_FLAG sort_flag;
} FEM_VELO_ARRAY;


/* array of reaction force */
typedef struct{
	int size;
	int alloc_size;
	FEM_DATA_TYPE source;
	FEM_DISP_VELO *array;
	FEM_ARRAY_SORT_FLAG sort_flag;
} FEM_REACT_ARRAY;


/* array of stress */
typedef struct {
	int size;
	int alloc_size;
	FEM_DATA_TYPE source;
	FEM_STRESS_STRAIN *array;
	FEM_ARRAY_SORT_FLAG sort_flag;
} FEM_STRESS_ARRAY;


/* array of strain */
typedef struct {
	int size;
	int alloc_size;
	FEM_DATA_TYPE source;
	FEM_STRESS_STRAIN *array;
	FEM_ARRAY_SORT_FLAG sort_flag;
} FEM_STRAIN_ARRAY;


typedef enum {
	SENS_SHAPE_GRAD_DENSITY,
	SENS_SHAPE_GRAD,
	SENS_DENSITY,
	SENS_THICKNESS,
	SENS_BEAM_XY,
	SENS_UNKNOWN
} FEM_SENSITIVITY_TYPE;


/* array of sensitivity */
typedef struct {
	int size;
	int alloc_size;
	FEM_SENSITIVITY_TYPE type;
	FEM_SENSITIVITY *array;
	FEM_ARRAY_SORT_FLAG sort_flag;
} FEM_SENSITIVITY_ARRAY;


/* array of default b.c. */
typedef struct {
	int size;
	int alloc_size;
	FEM_DATA_TYPE source;
	FEM_DEFAULT_BC *array;
	FEM_ARRAY_SORT_FLAG sort_flag;
} FEM_DEFAULT_BC_ARRAY;


/* array of sound pressure */
typedef struct {
	int size;
	int alloc_size;
	FEM_DATA_TYPE source;
	FEM_SOUND_PRESSURE *array;
	FEM_ARRAY_SORT_FLAG sort_flag;
} FEM_SOUND_PRESSURE_ARRAY;


/* array of element matrix */
typedef struct {
	int size;
	int alloc_size;
	FEM_DATA_TYPE source;
	FEM_ELEM_MATRIX *array;
	FEM_ARRAY_SORT_FLAG sort_flag;
} FEM_ELEM_MATRIX_ARRAY;


/* array of element lookup array */
typedef struct {
	int size;
	int alloc_size;
	FEM_DATA_TYPE source;
	FEM_ELEM_LOOKUP *array;
	FEM_ARRAY_SORT_FLAG sort_flag;
} FEM_ELEM_LOOKUP_ARRAY;


typedef struct {
	int size;
	int alloc_size;
	FEM_DATA_TYPE source;
	FEM_NODE_MATRIX *array;
	FEM_ARRAY_SORT_FLAG sort_flag;
} FEM_NODE_MATRIX_ARRAY;


/* array of load curve */
typedef struct {
	int size;
	int alloc_size;
	FEM_DATA_TYPE source;
	FEM_LOAD_CURVE *array;
	FEM_ARRAY_SORT_FLAG sort_flag;
} FEM_LOAD_CURVE_ARRAY;


typedef struct {
	int size;
	int alloc_size;
	FEM_DATA_TYPE source;
	FEM_LOCAL_COORD *array;
	FEM_ARRAY_SORT_FLAG sort_flag;
} FEM_LOCAL_COORD_ARRAY;


typedef enum {
	FEM_BC_NOT_COMMON,
	FEM_BC_COMMON
} FEM_BC_COMMON_FLAG;


/* array of b.c. set */
typedef struct {
	int size;
	int alloc_size;
	FEM_BC_COMMON_FLAG com_fix;
	FEM_BC_COMMON_FLAG com_mpc;
	FEM_BC_COMMON_FLAG com_force;
	FEM_BC_COMMON_FLAG com_temp;
	FEM_DATA_TYPE source;
	FEM_BC_SET *array;
	FEM_ARRAY_SORT_FLAG sort_flag;
	STRING_ARRAY s_info;
} FEM_BC_SET_ARRAY;


#define FEM_MAX_FACE       6
#define FEM_MAX_FACE_NODE  8
typedef struct {
	int face_num;
	int node_num[FEM_MAX_FACE];
	FEM_ELEM_TYPE face_type[FEM_MAX_FACE];
	int table[FEM_MAX_FACE][FEM_MAX_FACE_NODE];
} FEM_FACE_TABLE;


#define FEM_MAX_EDGE       12
#define FEM_MAX_EDGE_NODE  3
typedef struct {
	int edge_num;
	int node_num[FEM_MAX_EDGE];
	FEM_ELEM_TYPE edge_type[FEM_MAX_EDGE];
	int table[FEM_MAX_EDGE][FEM_MAX_EDGE_NODE];
} FEM_EDGE_TABLE;


/* 作業配列セット */
#define FEM_WS_SIZE    4
typedef struct {
	double *vec_N[FEM_WS_SIZE];
	double *vec_G[FEM_WS_SIZE];

	double **mat_3_3[FEM_WS_SIZE];
	double **mat_6_6[FEM_WS_SIZE];
	double **mat_9_9[FEM_WS_SIZE];

	double **mat_3_N[FEM_WS_SIZE];
	double **mat_3_3N[FEM_WS_SIZE];
	double **mat_6_3N[FEM_WS_SIZE];
	double **mat_9_3N[FEM_WS_SIZE];
	double **mat_3N_3N[FEM_WS_SIZE];

	double **mat_N_3[FEM_WS_SIZE];
	double **mat_N_4[FEM_WS_SIZE];
	double **mat_N_6[FEM_WS_SIZE];
	double **mat_G_4[FEM_WS_SIZE];
} FEM_WORK_SET;



/* fem_struct.c */
void init_fem_node(FEM_NODE *node);
void init_fem_element(FEM_ELEMENT *elem);
void init_fem_rigid_element(FEM_RIGID_ELEMENT *elem);
void init_fem_material_prop(FEM_MATERIAL_PROP *prop);
void init_fem_physical_prop(FEM_PHYSICAL_PROP *prop);
void init_fem_bc(FEM_BC *bc);
void init_fem_elem_bc(FEM_ELEM_BC *bc);
void init_fem_disp_velo(FEM_DISP_VELO *disp);
void init_fem_stress_strain(FEM_STRESS_STRAIN *ss);
void init_fem_sensitivity(FEM_SENSITIVITY *sens);
void init_fem_default_bc(FEM_DEFAULT_BC *def);
void init_fem_bc_set(FEM_BC_SET *bc_set);
void init_fem_sound_pressure(FEM_SOUND_PRESSURE *sound);
void init_fem_load_curve(FEM_LOAD_CURVE *curve);
void init_fem_elem_matrix(FEM_ELEM_MATRIX *elem_matrix);
void init_fem_elem_lookup(FEM_ELEM_LOOKUP *elem_lookup);
void init_fem_node_matrix(FEM_NODE_MATRIX *nmat);
void init_fem_local_coord(FEM_LOCAL_COORD *local_coord);
RC free_fem_elem_matrix(FEM_ELEM_MATRIX *elem_matrix);
RC free_fem_elem_lookup(FEM_ELEM_LOOKUP *elem_lookup);
RC free_fem_node_matrix(FEM_NODE_MATRIX *nmat);
RC free_fem_rigid_element(FEM_RIGID_ELEMENT *elem);
RC free_fem_bc_set(FEM_BC_SET *bc_set);
RC copy_fem_rigid_element(FEM_RIGID_ELEMENT src, FEM_RIGID_ELEMENT *dest);
FEM_STRESS_ARRAY fem_strain2stress_array(FEM_STRAIN_ARRAY strain);
FEM_STRAIN_ARRAY fem_stress2strain_array(FEM_STRESS_ARRAY stress);
RC fem_fill_material(FEM_MATERIAL_PROP *prop);
RC Gauss_const_default(FEM_ELEM_TYPE type, double **points,
                       double *weight, int *num_point);
RC Gauss_const_max(FEM_ELEM_TYPE type, double **points,
                   double *weight, int *num_point);
RC Gauss_const_1(FEM_ELEM_TYPE type, double **points, double *weight);
int element_ssnum(FEM_ELEM_TYPE type);
int element_dim(FEM_ELEM_TYPE type);
int element_order(FEM_ELEM_TYPE type);
void init_fem_face_table(FEM_FACE_TABLE *f_table);
RC make_fem_face_table(FEM_ELEM_TYPE type, FEM_FACE_TABLE *f_table);
void init_fem_edge_table(FEM_EDGE_TABLE *e_table);
RC make_fem_edge_table(FEM_ELEM_TYPE type, FEM_EDGE_TABLE *e_table);
RC trans_coord_s2v(FEM_ELEMENT s_elem, FEM_ELEMENT v_elem,
                   int coord_size, double **s_coords, double **v_coords);
RC allocate_fem_work_set(FEM_WORK_SET *ws);
RC free_fem_work_set(FEM_WORK_SET *ws);
int same_bc_set(int num1, const int id1[], double scale1[],
                int num2, const int id2[], double scale2[]);
RC cal_det_J(FEM_ELEMENT element, const double local_coord[],
             double **node_location, double *det_J,
             double **ws_mat_3_N, double **ws_mat_3_3);
RC cal_d_A(FEM_ELEMENT surf, const double local_coord[],
           double **node_location, TRANSLATION3D *d_A,
           double **ws_mat_3_N, double **ws_mat_3_3);
RC interpolate_disp(FEM_ELEMENT element, const double local_coord[],
                    double **disp_matrix, TRANS_ROTATION3D *disp);
RC realloc_fem_node_matrix(FEM_NODE_MATRIX *nmat);


/* fem_struct_print.c */
void print_fem_node(FILE *fp, FEM_NODE node);
void print_fem_elem_type(FILE *fp, FEM_ELEM_TYPE type);
void print_fem_element(FILE *fp, FEM_ELEMENT elem);
void print_fem_rigid_element(FILE *fp, const FEM_RIGID_ELEMENT *elem);
void print_fem_mat_type(FILE *fp, FEM_MAT_TYPE type);
void print_fem_analysis_type(FILE *fp, FEM_ANALYSIS_TYPE type);
void print_fem_material_prop(FILE *fp, FEM_MATERIAL_PROP m_prop);
void print_fem_section_type(FILE *fp, FEM_SECTION_TYPE type);
void print_fem_physical_prop(FILE *fp, FEM_PHYSICAL_PROP p_prop);
void print_fem_bc_type(FILE *fp, FEM_BC_TYPE type);
void print_fem_bc_type_3d(FILE *fp, FEM_BC_TYPE3D type);
void print_fem_bc(FILE *fp, FEM_BC bc);
void print_fem_elem_bc(FILE *fp, FEM_ELEM_BC bc);
void print_fem_disp_velo(FILE *fp, FEM_DISP_VELO d_velo);
void print_fem_stress_strain(FILE *fp, FEM_STRESS_STRAIN ss);
void print_fem_sensitivity(FILE *fp, FEM_SENSITIVITY sensi);
void print_fem_sound_pressure(FILE *fp, FEM_SOUND_PRESSURE sound);
void print_fem_load_curve(FILE *fp, const FEM_LOAD_CURVE *curve);
void print_fem_elem_matrix(FILE *fp, FEM_ELEM_MATRIX elem_matrix);
void print_fem_elem_lookup(FILE *fp, FEM_ELEM_LOOKUP elem_lookup);
void print_fem_local_coord(FILE *fp, FEM_LOCAL_COORD coord);
void print_fem_data_type(FILE *fp, FEM_DATA_TYPE type);
void print_fem_array_sort_flag(FILE *fp, FEM_ARRAY_SORT_FLAG sort_flag);
void print_fem_renumber_flag(FILE *fp, FEM_RENUMBER_FLAG renum_flag);
void print_fem_node_array(FILE *fp, FEM_NODE_ARRAY node);
void print_fem_element_array(FILE *fp, FEM_ELEMENT_ARRAY element);
void print_fem_rigid_element_array(FILE *fp, FEM_RIGID_ELEMENT_ARRAY element);
void print_fem_material_prop_array(FILE *fp, FEM_MATERIAL_PROP_ARRAY m_prop);
void print_fem_physical_prop_array(FILE *fp, FEM_PHYSICAL_PROP_ARRAY p_prop);
void print_fem_bc_array(FILE *fp, FEM_BC_ARRAY bc);
void print_fem_elem_bc_array(FILE *fp, FEM_ELEM_BC_ARRAY bc);
void print_fem_disp_array(FILE *fp, FEM_DISP_ARRAY disp);
void print_fem_velo_array(FILE *fp, FEM_VELO_ARRAY velo);
void print_fem_react_array(FILE *fp, FEM_REACT_ARRAY react);
void print_fem_stress_array(FILE *fp, FEM_STRESS_ARRAY stress);
void print_fem_strain_array(FILE *fp, FEM_STRAIN_ARRAY strain);
void print_fem_sensitivity_type(FILE *fp, FEM_SENSITIVITY_TYPE type);
void print_fem_sensitivity_array(FILE *fp, FEM_SENSITIVITY_ARRAY sensi);
void print_fem_default_bc(FILE *fp, FEM_DEFAULT_BC def);
void print_fem_default_bc_array(FILE *fp, FEM_DEFAULT_BC_ARRAY def);
void print_fem_bc_set(FILE *fp, FEM_BC_SET bc_set);
void print_fem_bc_set_array(FILE *fp, FEM_BC_SET_ARRAY bc_set);
void print_fem_bc_common_flag(FILE *fp, FEM_BC_COMMON_FLAG flag);
void print_fem_sound_pressure_array(FILE *fp, FEM_SOUND_PRESSURE_ARRAY sound);
void print_fem_load_curve_array(FILE *fp, FEM_LOAD_CURVE_ARRAY curve);
void print_fem_elem_matrix_array(FILE *fp, FEM_ELEM_MATRIX_ARRAY elem_matrix);
void print_fem_elem_lookup_array(FILE *fp, FEM_ELEM_LOOKUP_ARRAY elem_lookup);
void print_fem_local_coord_array(FILE *fp, FEM_LOCAL_COORD_ARRAY coord);
void print_fem_node_matrix(FILE *fp, FEM_NODE_MATRIX nmat);
void print_fem_node_matrix_array(FILE *fp, FEM_NODE_MATRIX_ARRAY nmat);
void print_system_info(FILE *fp);


/* fem_struct_binio.c */
RC write_fem_node_array(FILE *fp, FEM_NODE_ARRAY node);
RC write_fem_element_array(FILE *fp, FEM_ELEMENT_ARRAY element);
RC write_fem_material_prop_array(FILE *fp, FEM_MATERIAL_PROP_ARRAY m_prop);
RC write_fem_physical_prop_array(FILE *fp, FEM_PHYSICAL_PROP_ARRAY p_prop);
RC write_fem_bc_array(FILE *fp, FEM_BC_ARRAY bc);
RC write_fem_elem_bc_array(FILE *fp, FEM_ELEM_BC_ARRAY bc);
RC write_fem_disp_array(FILE *fp, FEM_DISP_ARRAY disp);
RC write_fem_velo_array(FILE *fp, FEM_VELO_ARRAY velo);
RC write_fem_react_array(FILE *fp, FEM_REACT_ARRAY react);
RC write_fem_stress_array(FILE *fp, FEM_STRESS_ARRAY stress);
RC write_fem_strain_array(FILE *fp, FEM_STRAIN_ARRAY strain);
RC write_fem_sensitivity_array(FILE *fp, FEM_SENSITIVITY_ARRAY sensi);
RC write_fem_default_bc_array(FILE *fp, FEM_DEFAULT_BC_ARRAY def);
RC write_fem_sound_pressure_array(FILE *fp, FEM_SOUND_PRESSURE_ARRAY sound);
RC write_fem_bc_set_array(FILE *fp, FEM_BC_SET_ARRAY bc_set);
RC write_fem_load_curve_array(FILE *fp, FEM_LOAD_CURVE_ARRAY curve);
RC write_fem_local_coord_array(FILE *fp, FEM_LOCAL_COORD_ARRAY coord);

RC read_fem_node_array(FILE *fp, FEM_NODE_ARRAY *node);
RC read_fem_element_array(FILE *fp, FEM_ELEMENT_ARRAY *element);
RC read_fem_material_prop_array(FILE *fp, FEM_MATERIAL_PROP_ARRAY *m_prop);
RC read_fem_physical_prop_array(FILE *fp, FEM_PHYSICAL_PROP_ARRAY *p_prop);
RC read_fem_bc_array(FILE *fp, FEM_BC_ARRAY *bc);
RC read_fem_elem_bc_array(FILE *fp, FEM_ELEM_BC_ARRAY *bc);
RC read_fem_disp_array(FILE *fp, FEM_DISP_ARRAY *disp);
RC read_fem_velo_array(FILE *fp, FEM_VELO_ARRAY *velo);
RC read_fem_react_array(FILE *fp, FEM_REACT_ARRAY *react);
RC read_fem_stress_array(FILE *fp, FEM_STRESS_ARRAY *stress);
RC read_fem_strain_array(FILE *fp, FEM_STRAIN_ARRAY *strain);
RC read_fem_sensitivity_array(FILE *fp, FEM_SENSITIVITY_ARRAY *sensi);
RC read_fem_default_bc_array(FILE *fp, FEM_DEFAULT_BC_ARRAY *def);
RC read_fem_sound_pressure_array(FILE *fp, FEM_SOUND_PRESSURE_ARRAY *sound);
RC read_fem_bc_set_array(FILE *fp, FEM_BC_SET_ARRAY *bc_set);
RC read_fem_load_curve_array(FILE *fp, FEM_LOAD_CURVE_ARRAY *curve);
RC read_fem_local_coord_array(FILE *fp, FEM_LOCAL_COORD_ARRAY *coord);


/* fem_utl.c */
RC free_fem_node_array(FEM_NODE_ARRAY *node);
RC free_fem_element_array(FEM_ELEMENT_ARRAY *element);
RC free_fem_rigid_element_array(FEM_RIGID_ELEMENT_ARRAY *element);
RC free_fem_material_prop_array(FEM_MATERIAL_PROP_ARRAY *mat);
RC free_fem_physical_prop_array(FEM_PHYSICAL_PROP_ARRAY *phys);
RC free_fem_bc_array(FEM_BC_ARRAY *bc);
RC free_fem_elem_bc_array(FEM_ELEM_BC_ARRAY *bc);
RC free_fem_disp_array(FEM_DISP_ARRAY *disp);
RC free_fem_velo_array(FEM_VELO_ARRAY *velo);
RC free_fem_react_array(FEM_REACT_ARRAY *react);
RC free_fem_stress_array(FEM_STRESS_ARRAY *stress);
RC free_fem_strain_array(FEM_STRAIN_ARRAY *strain);
RC free_fem_sensitivity_array(FEM_SENSITIVITY_ARRAY *sens);
RC free_fem_default_bc_array(FEM_DEFAULT_BC_ARRAY *def);
RC free_fem_sound_pressure_array(FEM_SOUND_PRESSURE_ARRAY *sound);
RC free_fem_bc_set_array(FEM_BC_SET_ARRAY *bc_set);
RC free_fem_load_curve_array(FEM_LOAD_CURVE_ARRAY *curve);
RC free_fem_elem_matrix_array(FEM_ELEM_MATRIX_ARRAY *elem_matrix);
RC free_fem_elem_lookup_array(FEM_ELEM_LOOKUP_ARRAY *elem_lookup);
RC free_fem_node_matrix_array(FEM_NODE_MATRIX_ARRAY *node_matrix);
RC free_fem_local_coord_array(FEM_LOCAL_COORD_ARRAY *coord);

RC copy_fem_node_array(FEM_NODE_ARRAY src, FEM_NODE_ARRAY *dest);
RC copy_fem_element_array(FEM_ELEMENT_ARRAY src, FEM_ELEMENT_ARRAY *dest);
RC copy_fem_material_prop_array(FEM_MATERIAL_PROP_ARRAY src,
                                FEM_MATERIAL_PROP_ARRAY *dest);
RC copy_fem_physical_prop_array(FEM_PHYSICAL_PROP_ARRAY src,
                                FEM_PHYSICAL_PROP_ARRAY *dest);
RC copy_fem_bc_array(FEM_BC_ARRAY src, FEM_BC_ARRAY *dest);
RC copy_fem_elem_bc_array(FEM_ELEM_BC_ARRAY src, FEM_ELEM_BC_ARRAY *dest);
RC copy_fem_disp_array(FEM_DISP_ARRAY src, FEM_DISP_ARRAY *dest);
RC copy_fem_velo_array(FEM_VELO_ARRAY src, FEM_VELO_ARRAY *dest);
RC copy_fem_react_array(FEM_REACT_ARRAY src, FEM_REACT_ARRAY *dest);
RC copy_fem_stress_array(FEM_STRESS_ARRAY src, FEM_STRESS_ARRAY *dest);
RC copy_fem_strain_array(FEM_STRAIN_ARRAY src, FEM_STRAIN_ARRAY *dest);
RC copy_fem_sensitivity_array(FEM_SENSITIVITY_ARRAY src,
                              FEM_SENSITIVITY_ARRAY *dest);
RC copy_fem_default_bc_array(FEM_DEFAULT_BC_ARRAY src,
                             FEM_DEFAULT_BC_ARRAY *dest);
RC copy_fem_sound_pressure_array(FEM_SOUND_PRESSURE_ARRAY src,
                                 FEM_SOUND_PRESSURE_ARRAY *dest);
RC copy_fem_load_curve_array(FEM_LOAD_CURVE_ARRAY src,
                             FEM_LOAD_CURVE_ARRAY *dest);
RC copy_fem_local_coord_array(FEM_LOCAL_COORD_ARRAY src,
                              FEM_LOCAL_COORD_ARRAY *dest);

RC clean_fem_node_array(FEM_NODE_ARRAY *node);
RC clean_fem_element_array(FEM_ELEMENT_ARRAY *elem);
RC clean_fem_rigid_element_array(FEM_RIGID_ELEMENT_ARRAY *elem);
RC clean_fem_material_prop_array(FEM_MATERIAL_PROP_ARRAY *mat);
RC clean_fem_physical_prop_array(FEM_PHYSICAL_PROP_ARRAY *phys);
RC clean_fem_bc_array(FEM_BC_ARRAY *bc);
RC clean_fem_elem_bc_array(FEM_ELEM_BC_ARRAY *bc);
RC clean_fem_disp_array(FEM_DISP_ARRAY *disp);
RC clean_fem_velo_array(FEM_VELO_ARRAY *velo);
RC clean_fem_react_array(FEM_REACT_ARRAY *react);
RC clean_fem_stress_array(FEM_STRESS_ARRAY *stress);
RC clean_fem_strain_array(FEM_STRAIN_ARRAY *strain);
RC clean_fem_sensitivity_array(FEM_SENSITIVITY_ARRAY *sens);
RC clean_fem_default_bc_array(FEM_DEFAULT_BC_ARRAY *def);
RC clean_fem_sound_pressure_array(FEM_SOUND_PRESSURE_ARRAY *sound);
RC clean_fem_bc_set_array(FEM_BC_SET_ARRAY *bc_set);
RC clean_fem_elem_matrix_array(FEM_ELEM_MATRIX_ARRAY *elem_matrix);
RC clean_fem_elem_lookup_array(FEM_ELEM_LOOKUP_ARRAY *elem_lookup);
RC clean_fem_node_matrix_array(FEM_NODE_MATRIX_ARRAY *node_matrix);
RC clean_fem_load_curve_array(FEM_LOAD_CURVE_ARRAY *curve);
RC clean_fem_local_coord_array(FEM_LOCAL_COORD_ARRAY *coord);

RC realloc_fem_node_array(FEM_NODE_ARRAY *node);
RC realloc_fem_element_array(FEM_ELEMENT_ARRAY *elem);
RC realloc_fem_rigid_element_array(FEM_RIGID_ELEMENT_ARRAY *elem);
RC realloc_fem_material_prop_array(FEM_MATERIAL_PROP_ARRAY *mat);
RC realloc_fem_physical_prop_array(FEM_PHYSICAL_PROP_ARRAY *phys);
RC realloc_fem_bc_array(FEM_BC_ARRAY *bc);
RC realloc_fem_elem_bc_array(FEM_ELEM_BC_ARRAY *bc);
RC realloc_fem_disp_array(FEM_DISP_ARRAY *disp);
RC realloc_fem_velo_array(FEM_VELO_ARRAY *velo);
RC realloc_fem_react_array(FEM_REACT_ARRAY *react);
RC realloc_fem_stress_array(FEM_STRESS_ARRAY *stress);
RC realloc_fem_strain_array(FEM_STRAIN_ARRAY *strain);
RC realloc_fem_sensitivity_array(FEM_SENSITIVITY_ARRAY *sens);
RC realloc_fem_default_bc_array(FEM_DEFAULT_BC_ARRAY *def);
RC realloc_fem_sound_pressure_array(FEM_SOUND_PRESSURE_ARRAY *sound);
RC realloc_fem_bc_set_array(FEM_BC_SET_ARRAY *bc_set);
RC realloc_fem_elem_matrix_array(FEM_ELEM_MATRIX_ARRAY *elem_matrix);
RC realloc_fem_elem_lookup_array(FEM_ELEM_LOOKUP_ARRAY *elem_lookup);
RC realloc_fem_node_matrix_array(FEM_NODE_MATRIX_ARRAY *node_matrix);
RC realloc_fem_load_curve_array(FEM_LOAD_CURVE_ARRAY *curve);
RC realloc_fem_local_coord_array(FEM_LOCAL_COORD_ARRAY *local);

RC allocate_fem_node_array(int num, FEM_NODE_ARRAY *node);
RC allocate_fem_element_array(int num, FEM_ELEMENT_ARRAY *element);
RC allocate_fem_rigid_element_array(int num, FEM_RIGID_ELEMENT_ARRAY *element);
RC allocate_fem_material_prop_array(int num, FEM_MATERIAL_PROP_ARRAY *mat);
RC allocate_fem_physical_prop_array(int num, FEM_PHYSICAL_PROP_ARRAY *phys);
RC allocate_fem_bc_array(int num, FEM_BC_ARRAY *bc);
RC allocate_fem_elem_bc_array(int num, FEM_ELEM_BC_ARRAY *bc);
RC allocate_fem_disp_array(int num, FEM_DISP_ARRAY *disp);
RC allocate_fem_velo_array(int num, FEM_VELO_ARRAY *velo);
RC allocate_fem_react_array(int num, FEM_REACT_ARRAY *react);
RC allocate_fem_stress_array(int num, FEM_STRESS_ARRAY *stress);
RC allocate_fem_strain_array(int num, FEM_STRAIN_ARRAY *strain);
RC allocate_fem_sensitivity_array(int num, FEM_SENSITIVITY_ARRAY *sens);
RC allocate_fem_default_bc_array(int num, FEM_DEFAULT_BC_ARRAY *def);
RC allocate_fem_sound_pressure_array(int num, FEM_SOUND_PRESSURE_ARRAY *sound);
RC allocate_fem_bc_set_array(int num, FEM_BC_SET_ARRAY *bc_set);
RC allocate_fem_elem_matrix_array(int num, FEM_ELEM_MATRIX_ARRAY *elem_matrix);
RC allocate_fem_elem_lookup_array(int num, FEM_ELEM_LOOKUP_ARRAY *elem_lookup);
RC allocate_fem_node_matrix_array(int num, FEM_NODE_MATRIX_ARRAY *node_matrix);
RC allocate_fem_load_curve_array(int num, FEM_LOAD_CURVE_ARRAY *curve);
RC allocate_fem_local_coord_array(int num, FEM_LOCAL_COORD_ARRAY *coord);

int count_valid_node(FEM_NODE_ARRAY node);
int count_valid_element(FEM_ELEMENT_ARRAY elem);
int count_valid_rigid_element(FEM_RIGID_ELEMENT_ARRAY elem);
int count_valid_material_prop(FEM_MATERIAL_PROP_ARRAY material);
int count_valid_physical_prop(FEM_PHYSICAL_PROP_ARRAY phys);
int count_valid_bc(FEM_BC_ARRAY bc);
int count_valid_elem_bc(FEM_ELEM_BC_ARRAY bc);
int count_valid_disp(FEM_DISP_ARRAY disp);
int count_valid_velo(FEM_VELO_ARRAY velo);
int count_valid_react(FEM_REACT_ARRAY react);
int count_valid_load_curve(FEM_LOAD_CURVE_ARRAY curve);
int count_valid_local_coord(FEM_LOCAL_COORD_ARRAY coord);

RC sort_fem_node_array(FEM_NODE_ARRAY *node);
RC sort_fem_node_renum(FEM_NODE_ARRAY *node);
RC sort_fem_element_array(FEM_ELEMENT_ARRAY *elem);
RC sort_fem_rigid_element_array(FEM_RIGID_ELEMENT_ARRAY *elem);
RC sort_fem_material_prop_array(FEM_MATERIAL_PROP_ARRAY *mat);
RC sort_fem_physical_prop_array(FEM_PHYSICAL_PROP_ARRAY *phys);
RC sort_fem_bc_array(FEM_BC_ARRAY *bc);
RC sort_fem_elem_bc_array(FEM_ELEM_BC_ARRAY *bc);
RC sort_fem_disp_array(FEM_DISP_ARRAY *disp);
RC sort_fem_velo_array(FEM_VELO_ARRAY *velo);
RC sort_fem_react_array(FEM_REACT_ARRAY *react);
RC sort_fem_stress_array(FEM_STRESS_ARRAY *stress);
RC sort_fem_strain_array(FEM_STRAIN_ARRAY *strain);
RC sort_fem_sensitivity_array(FEM_SENSITIVITY_ARRAY *sens);
RC sort_fem_default_bc_array(FEM_DEFAULT_BC_ARRAY *def);
RC sort_fem_sound_pressure_array(FEM_SOUND_PRESSURE_ARRAY *sound);
RC sort_fem_bc_set_array(FEM_BC_SET_ARRAY *bc_set);
RC sort_fem_elem_matrix_array(FEM_ELEM_MATRIX_ARRAY *elem_matrix);
RC sort_fem_elem_lookup_array(FEM_ELEM_LOOKUP_ARRAY *elem_lookup);
RC sort_fem_node_matrix_array(FEM_NODE_MATRIX_ARRAY *node_matrix);
RC sort_fem_load_curve_array(FEM_LOAD_CURVE_ARRAY *curve);
RC sort_fem_local_coord_array(FEM_LOCAL_COORD_ARRAY *coord);

int search_fem_node_label(FEM_NODE_ARRAY node, int label);
int search_fem_element_label(FEM_ELEMENT_ARRAY elem, int label);
int search_fem_rigid_element_label(FEM_RIGID_ELEMENT_ARRAY elem, int label);
int search_fem_material_prop_label(FEM_MATERIAL_PROP_ARRAY mat, int label);
int search_fem_physical_prop_label(FEM_PHYSICAL_PROP_ARRAY phys, int label);
int search_fem_bc_node(FEM_BC_ARRAY bc, int node);
int search_fem_elem_bc_element(FEM_ELEM_BC_ARRAY bc, int element);
int search_fem_disp_node(FEM_DISP_ARRAY disp, int node);
int search_fem_velo_node(FEM_VELO_ARRAY velo, int node);
int search_fem_react_node(FEM_REACT_ARRAY react, int node);
int search_fem_stress_element(FEM_STRESS_ARRAY stress, int element);
int search_fem_strain_element(FEM_STRAIN_ARRAY strain, int element);
int search_fem_sensitivity_element(FEM_SENSITIVITY_ARRAY sens, int element);
int search_fem_default_bc_set_id(FEM_DEFAULT_BC_ARRAY def, int set_id);
int search_fem_sound_pressure_node(FEM_SOUND_PRESSURE_ARRAY sound, int node);
int search_fem_bc_set_label(FEM_BC_SET_ARRAY bc_set, int label);
int search_fem_elem_matrix_element(FEM_ELEM_MATRIX_ARRAY elem_matrix,
                                   int element);
int search_fem_elem_lookup_node(FEM_ELEM_LOOKUP_ARRAY elem_lookup, int node);
int search_fem_node_matrix_node(FEM_NODE_MATRIX_ARRAY node_matrix, int node);
int search_fem_load_curve_label(FEM_LOAD_CURVE_ARRAY curve, int label);
int search_fem_local_coord_label(FEM_LOCAL_COORD_ARRAY coord, int label);
int search_fem_material_prop_element(FEM_MATERIAL_PROP_ARRAY mat,
                                     FEM_PHYSICAL_PROP_ARRAY phys,
                                     FEM_ELEMENT *element);
int search_fem_stress_element_node(FEM_STRESS_ARRAY stress,
                                   int element, int node);
int search_fem_strain_element_node(FEM_STRAIN_ARRAY stress,
                                   int element, int node);
int search_fem_sensitivity_element_node(FEM_SENSITIVITY_ARRAY sens,
                                        int element, int node);

RC add_fem_disp_array(double weight, FEM_DISP_ARRAY disp,
                      FEM_DISP_ARRAY *disp_total);
RC add_fem_sensitivity_array(double weight, FEM_SENSITIVITY_ARRAY sens,
                             FEM_SENSITIVITY_ARRAY *sens_total);
RC node_location_matrix(FEM_ELEMENT *element, 
                        FEM_NODE_ARRAY node, double **matrix);
RC fill_disp_matrix(FEM_ELEMENT element,
                    FEM_DISP_ARRAY disp, double **matrix);
RC fill_sound_pres_matrix(FEM_ELEMENT element,
                          FEM_SOUND_PRESSURE_ARRAY pres, double **matrix);
RC unit_normal_vect_center(FEM_ELEMENT surf, FEM_ELEMENT inner_elem,
                           FEM_NODE_ARRAY node, TRANSLATION3D *vect);
RC unit_normal_vect(double *local_coord, FEM_ELEMENT surf, 
                    FEM_ELEMENT inner_elem, FEM_NODE_ARRAY node,
                    TRANSLATION3D *vect);
RC extract_bc_set(int num, const int id[], double scale[],
                  FEM_BC_ARRAY src, FEM_BC_ARRAY *dest);
RC extract_default_bc_set(int num, const int id[], double scale[],
                          FEM_DEFAULT_BC_ARRAY src, 
                          FEM_DEFAULT_BC_ARRAY *dest);
RC extract_rigid_set(int num, const int id[], double scale[],
                  FEM_RIGID_ELEMENT_ARRAY src, FEM_RIGID_ELEMENT_ARRAY *dest);
RC extract_elem_bc_set(int num, const int id[], double scale[],
                       FEM_ELEM_BC_ARRAY src, FEM_ELEM_BC_ARRAY *dest);
RC compress_fem_bc_array(FEM_BC_ARRAY *dest);
RC mul_elem_matrices_vector(int total_dof, FEM_ELEM_MATRIX_ARRAY elem_matrix,
                            const double *vector, double *answer);
RC make_elem_lookup(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                    FEM_ELEM_LOOKUP_ARRAY *elem_lookup);
RC make_elem_lookup_rigid(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                          FEM_RIGID_ELEMENT_ARRAY rigid_elem,
                          FEM_ELEM_LOOKUP_ARRAY *elem_lookup);
int fem_analysis_dim(FEM_ELEMENT_ARRAY element);
int fem_analysis_order(FEM_ELEMENT_ARRAY element);
RC reduced_element_array_order(FEM_ELEMENT_ARRAY *elem);
RC reduced_element_type(FEM_ELEM_TYPE old_type, FEM_ELEM_TYPE *new_type,
                        int *new_node_num);
RC delete_duplicate_node(FEM_ELEMENT_ARRAY element, FEM_NODE_ARRAY *node,
                         double cr_dist);
RC minimum_node_distance(FEM_ELEMENT_ARRAY element, FEM_NODE_ARRAY node,
                         double *dist);
RC delete_unused_node(FEM_ELEMENT_ARRAY elem, FEM_NODE_ARRAY *node);
RC delete_unused_node_r(FEM_ELEMENT_ARRAY elem, FEM_RIGID_ELEMENT_ARRAY rigid,
                        FEM_NODE_ARRAY *node);
RC delete_unused_bc(FEM_NODE_ARRAY node, FEM_BC_ARRAY *bc);
RC locate_middle_node(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY elem);
RC interpolate_table(FEM_ELEM_TYPE type, int table[][2]);
RC element_average_length(FEM_ELEMENT element, FEM_NODE_ARRAY node,
                          double *length);
RC element_edge_table(FEM_ELEM_TYPE type, int *edge_num, int table[][2]);
RC elem_matrix_cg(int total_dof, FEM_ELEM_MATRIX_ARRAY elem_matrix,
                  const double *vector, double *solution);


/* fem_renumber.c */
RC dummy_renumber(FEM_NODE_ARRAY *node);
RC fem_renumber0(FEM_NODE_ARRAY *node, FEM_ELEMENT_ARRAY element,
                 FEM_RIGID_ELEMENT_ARRAY rigid,
                 int type_flag, int verbose_flag);
int node_renum_index(FEM_NODE_ARRAY node, int node_label);
int node_renum_index_cache(FEM_NODE_ARRAY node, int node_label,
                           int cache_size, int cache[][2]);


/* element_volume.c */
RC total_element_volume(FEM_ELEMENT_ARRAY element,
                        FEM_NODE_ARRAY node, double *total_volume);
RC element_volume(FEM_ELEMENT *element, FEM_NODE_ARRAY node);


/* extract_surface.c */
RC extract_surface(FEM_ELEMENT_ARRAY elem, FEM_ELEMENT_ARRAY *surf);
RC extract_divider(FEM_ELEMENT_ARRAY elem, FEM_ELEMENT_ARRAY *div);
RC cut_element(FEM_NODE_ARRAY *node, FEM_ELEMENT_ARRAY *elem);


/* interpolation.c */
RC set_local_node_points(FEM_ELEM_TYPE type, double **coords);
RC set_center_point(FEM_ELEM_TYPE type, double *coord);
RC set_N_vector(FEM_ELEM_TYPE type, const double *coord, double *N);
RC set_dN_matrix(FEM_ELEM_TYPE type, const double *coord, double **dN);


/* pressure_force.c */
RC pressure_force(int p_size, const double *p, FEM_ELEMENT *element,
                  FEM_NODE_ARRAY node, TRANSLATION3D *f);
RC traction_force(int p_size, const double *p, TRANSLATION3D dir,
            FEM_ELEMENT *element, FEM_NODE_ARRAY node, TRANSLATION3D *f);


#endif /* FEM_STRUCT_H */


