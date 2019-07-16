/*********************************************************************
 * fem_struct.h
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Tmoyuki TSUBATA>
 *
 * Thanks to following contributors.
 *  <Yasuo ASAGA> <Masanobu SUYAMA> <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: fem_struct.h 1082 2016-11-17 14:12:12Z hayashi $ */

#ifndef FEM_STRUCT_H
#define FEM_STRUCT_H

#include <stdio.h>
#include "rc.h"
#include "string_utl.h"
#include "math_utl.h"

#define MAX_GAUSS_POINT (64)

/* node */
#define NODE_INFO    4
typedef struct {
	int label;   /* label of node(if < 0 then unused) */
	VECT3D p;    /* x,y,z coordinate value */
	int renum;   /* for renumbering */
	int dummy_label;  /* reference label(if < 0 then not dummy) */
	int i_info[NODE_INFO];     /* the other informations(int) */
	double d_info[NODE_INFO];  /* the other informations(double) */
} NODE;


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
	ELEM_QUAD2L,
	ELEM_TETRA1,
	ELEM_TETRA2,
	ELEM_PENTA1,
	ELEM_PENTA2,
	ELEM_HEXA1,
	ELEM_HEXA2
} ELEM_TYPE;


/* element */
#define MAX_NODE     20
#define ELEMENT_INFO  9
typedef struct {
	int label;                   /* label of element(if < 0 then unused) */
	ELEM_TYPE type;              /* element type */
	int material;                /* material property label */
	int physical;                /* physical property label */
	int node_num;                /* number of node */
	int node[MAX_NODE];          /* node labels */
	double volume;               /* volume or area of element */
	int i_info[ELEMENT_INFO];    /* another information (int) */
	double d_info[ELEMENT_INFO]; /* another information (double) */
} ELEMENT;


/* rigid or mpc element */
#define RIGID_ELEMENT_INFO    10
typedef struct {
	int label;        /* label of element(if < 0 then unused) */
	int set_id;       /* ID of bc set (if < 0 then permanent element) */
	ELEM_TYPE type;   /* element type */
	int physical;     /* physical property label */
	int node_num;     /* number of node */
	int *node;        /* node labels */
	int weight_num;   /* size of weight[] */
	double *weight;   /* weight */
	int dof_num;      /* size of dof[] */
	int *dof;         /* DOF flag (0:x, 1:y, 2:z, 3:yz, 4:zx, 5:xy) */
	int i_info[RIGID_ELEMENT_INFO];    /* information (int) */
	double d_info[RIGID_ELEMENT_INFO]; /* information (double) */
} RIGID_ELEMENT;


/* type of material */
typedef enum {
	MAT_ISOTROPIC,
	MAT_ORTHOTROPIC,
	MAT_ANISOTROPIC,
	MAT_MOONEY_RIVLIN,        /* Mooney-Rivlin model */
	MAT_MOONEY_RIVLIN_CUBIC   /* Mooney-Rivlin model (cubic equation) */
} MAT_TYPE;

/* state of material */
typedef enum {
    STATE_SOLID,
    STATE_FLUID
} MAT_STATE;


/* type of analysis */
typedef enum {
	ANA_3D,
	ANA_PLANE_STRESS,
	ANA_PLANE_STRAIN,
	ANA_AXISYMMETRY,
	ANA_PLATE_BEND,
	ANA_PLATE_SHEAR
} ANALYSIS_TYPE;


/* material property */
#define MATERIAL_INFO 9
typedef struct {
	int label;            /* label of material (if < 0 then unused) */
	MAT_TYPE mat_type;    /* material type */
    MAT_STATE mat_state;  /* material state */
	ANALYSIS_TYPE ana_type;  /* analysis type */

	/* for mat_type == MAT_ISOTROPIC */
	double E;               /* Young's modulus */
	double G;               /* Shear modulus */
	double nu;              /* Poisson's ratio */ 
	double alpha;           /* coefficient of thermal expansion */
	double K;               /* 熱伝導率 */

	/* for mat_type == MAT_ORTHOTROPIC */
	VECT3D E_vect;       /* Young's modulus(x,y,z) */
	VECT3DR G_vect;      /* Shear modulus(yz,zx,xy) */
	VECT3DR nu_vect;     /* Poisson's ratio(yz,zx,xy) */
	VECT3DR alpha_vect;  /* coefficient of thermal expansion */
	VECT3DR K_vect;      /* (コメント文字化け) */

	/* for mat_type == MAT_UNISOTROPIC */
	double D_matrix[6][6];  /* stress - strain matrix */

	/* for hyper-elastic */
	double c[MATERIAL_INFO];   /* coefficient of strain energy function */

	double C[3][3][3][3];   /* stiffness tensor */
	double rho;             /* mass density */
	double sheat;           /* specific heat*/
	double tref;            /* reference temperature */

    /* for fluid */
    double bulk;            /* bulk modulus */
    double speed;           /* speed of sound */
    double damp;            /* element damping coefficient */

	int i_info[MATERIAL_INFO];    /* another information (int) */
	double d_info[MATERIAL_INFO]; /* another information (double) */
} MATERIAL_PROP;



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
} SECTION_TYPE;


/* physical property */
#define PHYSICAL_INFO 24
typedef struct {
	int label;           /* label of physical property (if < 0 then unused) */
	SECTION_TYPE section;         /* cross section type of beam */
	int i_info[PHYSICAL_INFO];    /* another information (int) */
	double d_info[PHYSICAL_INFO]; /* another information (double) */
} PHYSICAL_PROP;


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
	BC_K_DIAGONAL,
	BC_DENSITY,       /* 位相最適化の密度 */
	BC_WEIGHT         /* 重み */
} BC_TYPE;


/* type of b.c. on each DOF */
typedef struct {
	BC_TYPE x;
	BC_TYPE y;
	BC_TYPE z;
	BC_TYPE yz;
	BC_TYPE zx;
	BC_TYPE xy;
} BC_TYPE3D;


/* b.c. on node */
#define BC_INFO    4
#define BC_SCALAR  2
typedef struct {
	int node;                   /* label of node(if < 0 then unused) */
	int set_id;                 /* ID of bc set */
	BC_TYPE3D v_type;           /* type of bc on each DOF */
	VECT3DR v;                  /* bc on each DOF */
	BC_TYPE s_type[BC_SCALAR];  /* type of scalar bc */
	double s[BC_SCALAR];        /* scalar values */
	int i_info[BC_INFO];        /* another information (int) */
	double d_info[BC_INFO];     /* another information (double) */
} BC;


/* b.c. on element */
#define ELEM_BC_INFO    4
#define ELEM_BC_SCALAR  2
typedef struct {
	int element;              /* label of element(if < 0 then unused) */
	int set_id;               /* ID of bc set */
	BC_TYPE type;             /* type of B.C. */
	int face_id;              /* ID of face */
	union {
		double s;             /* uniform scalar value */
		VECT3DR v;            /* uniform vector value */
		double ns[MAX_NODE];  /* scalar value defined on each node */
	} val;
	VECT3D dir;
	int i_info[ELEM_BC_INFO];      /* another information (int) */
	double d_info[ELEM_BC_INFO];   /* another information (double) */
} ELEM_BC;


/* displacement or velocity on node */
#define DISP_VELO_INFO    4
typedef struct {
	int node;                      /* label of node(if < 0 then unused) */
	VECT3DR v;                     /* displacement or velocity */
	double s;                      /* pressure, scalar value */
	int i_info[DISP_VELO_INFO];    /* another information (int) */
	double d_info[DISP_VELO_INFO]; /* another information (double) */
} DISP_VELO;


/* stress or strain on node of element */
#define STRESS_STRAIN_INFO 4
typedef struct {
	int element;     /* label of element(if < 0 then unused) */
	int node;        /* label of node (if -1 then center of element) */
	VECT3DR v;       /* stress(or strain) vector (Voigt notation)*/
	/*double tensor[3][3];*/  /* stress(or strain) tensor */
	/*double mises;*/         /* von mises */
	/*double principal[3];*/  /* principal stress(or strain) */
	/*double p_dir[3][3];*/   /* direction of principal stress(or strain) */
	int i_info[STRESS_STRAIN_INFO];    /* another information (int) */
	double d_info[STRESS_STRAIN_INFO]; /* another information (double) */
} STRESS_STRAIN;


/* sensitivity on node of element */
typedef struct {
	int element;     /* label of element(if < 0 then unused) */
	int node;        /* label of node (if < 0 then center of element) */
	VECT3DR v;       /* sensitivity(vector) */
	double s;        /* sensitivity(scalar) */
} SENSITIVITY;


/* default temperature and gravity force */
#define DEFAULT_BC_INFO 4
typedef struct {
	int set_id;         /* ID of bc set (if < 0 then unused) */
	double temp;        /* default temperature */
	VECT3D grav;        /* gravity force vector */
	int i_info[DEFAULT_BC_INFO];    /* another information (int) */
	double d_info[DEFAULT_BC_INFO]; /* another information (double) */
} DEFAULT_BC;


/* b.c. case set */
#define MAX_BC_SET 32
#define BC_SET_INFO 4
typedef struct {
	int label;
	int num_fix;
	double scale_fix[MAX_BC_SET];
	int id_fix[MAX_BC_SET];
	int num_mpc;
	double scale_mpc[MAX_BC_SET];
	int id_mpc[MAX_BC_SET];
	int num_force;
	double scale_force[MAX_BC_SET];
	int id_force[MAX_BC_SET];
	int num_temp;
	double scale_temp[MAX_BC_SET];
	int id_temp[MAX_BC_SET];
	int i_info[BC_SET_INFO];
	double d_info[BC_SET_INFO];
	STRING_ARRAY s_info;
} BC_SET;


/* pressure of sound field */
#define SOUND_PRESSURE_INFO 4
typedef struct{
	int node;    /* label of node (if < 0 then unused) */
	double re;   /* real part */
	double im;   /* imaginary part */
	int i_info[SOUND_PRESSURE_INFO];
	double d_info[SOUND_PRESSURE_INFO];
} SOUND_PRESSURE;


#define LOAD_CURVE_MAX_INC    1024
#define LOAD_CURVE_MAX_APPROX   50
typedef struct {
	int label;
	double max_time;            /* load(disp) maximum time  */
	double max_value;           /* load(disp) maximum value */
	int total_step;             /* number of load(disp) step */
	double inc_value[LOAD_CURVE_MAX_INC];
	                            /* increment load(disp)value at every step */
	int approx_size;            /* number of load(disp) division point */
	double approx_time[LOAD_CURVE_MAX_APPROX];
	double approx_value[LOAD_CURVE_MAX_APPROX];
} LOAD_CURVE;


#define LOCAL_COORD_INFO 16
typedef struct {
	int label;
	VECT3D origin;   /* 原点 */
	VECT3D d_cos[3]; /* 方向余弦 */
	int i_info[LOCAL_COORD_INFO];
	double d_info[LOCAL_COORD_INFO];
} LOCAL_COORD;


/* matrix of element */
typedef struct{
	int element;      /* label of element (if < 0 then unused) */
	int size;         /* total size */
	int dim;          /* dof on each node */
	int *index;       /* index on global matrix */
	double **matrix;  /* element matrix */
} ELEM_MATRIX;


/* node -> element lookup table */
typedef struct{
	int node;         /* label of node (if < 0 then unused) */
	int size;         /* size of lookup table */
	int alloc_size;   /* allocated size of lookup table */
	int *element;     /* lookup table of element label */
	int r_size;       /* size of lookup table for rigid element */
	int r_alloc_size; /* allocated size of lookup table for rigid element */
	int *r_element;   /* lookup table of element label for rigid element */
} ELEM_LOOKUP;


typedef struct{
	int node;
	int shift;               /* 0: 並進自由度， 1: 回転自由度 */
	int size;
	int alloc_size;
	int *col_node;
	int *col_shift;          /* 0: 並進自由度， 1: 回転自由度 */
	double (*matrix)[3][3];
} NODE_MATRIX;


/* contact type */
typedef enum{
	CONTACT_PLANE,
	CONTACT_LINE,
	CONTACT_POINT
} CONTACT_TYPE;


/* contact terminus */
typedef struct{
	int element;       /* label of target element (if < 0 then unused) */
	VECT3D p;          /* terminus point of the target node */ 
	double distance;   /* distance of node to term */
	CONTACT_TYPE type; /* type of contact point */
} CONTACT_TERM;


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
} DATA_TYPE;


/* renumber flag */
typedef enum {
	RENUMBERED,
	UNRENUMBERED
} RENUMBER_FLAG;


/* Sort flag */
typedef enum {
	ARRAY_SORTED,
	ARRAY_UNSORTED
} ARRAY_SORT_FLAG;


/* array of node */
typedef struct {
	int size;
	int alloc_size;
	DATA_TYPE source;
	NODE *array;
	RENUMBER_FLAG renum_flag;
	ARRAY_SORT_FLAG sort_flag;
} NODE_ARRAY;


/* array of element */
typedef struct {
	int size;
	int alloc_size;
	DATA_TYPE source;
	ELEMENT *array;
	ARRAY_SORT_FLAG sort_flag;
} ELEMENT_ARRAY;


/* array of rigid element */
typedef struct {
	int size;
	int alloc_size;
	DATA_TYPE source;
	RIGID_ELEMENT *array;
	ARRAY_SORT_FLAG sort_flag;
} RIGID_ELEMENT_ARRAY;


/* array of material */
typedef struct {
	int size;
	int alloc_size;
	DATA_TYPE source;
	MATERIAL_PROP *array;
	ARRAY_SORT_FLAG sort_flag;
} MATERIAL_PROP_ARRAY;


/* array of physical property */
typedef struct {
	int size;
	int alloc_size;
	DATA_TYPE source;
	PHYSICAL_PROP *array;
	ARRAY_SORT_FLAG sort_flag;
} PHYSICAL_PROP_ARRAY;


/* array of b.c. */
typedef struct {
	int size;
	int alloc_size;
	DATA_TYPE source;
	BC *array;
	ARRAY_SORT_FLAG sort_flag;
} BC_ARRAY;


/* array of b.c. */
typedef struct {
	int size;
	int alloc_size;
	DATA_TYPE source;
	ELEM_BC *array;
	ARRAY_SORT_FLAG sort_flag;
} ELEM_BC_ARRAY;


/* array of displacement */
typedef struct{
	int size;
	int alloc_size;
	DATA_TYPE source;
	DISP_VELO *array;
	ARRAY_SORT_FLAG sort_flag;
} DISP_ARRAY;


/* array of velocity */
typedef struct{
	int size;
	int alloc_size;
	DATA_TYPE source;
	DISP_VELO *array;
	ARRAY_SORT_FLAG sort_flag;
} VELO_ARRAY;


/* array of reaction force */
typedef struct{
	int size;
	int alloc_size;
	DATA_TYPE source;
	DISP_VELO *array;
	ARRAY_SORT_FLAG sort_flag;
} REACT_ARRAY;


/* array of stress */
typedef struct {
	int size;
	int alloc_size;
	DATA_TYPE source;
	STRESS_STRAIN *array;
	ARRAY_SORT_FLAG sort_flag;
} STRESS_ARRAY;


/* array of strain */
typedef struct {
	int size;
	int alloc_size;
	DATA_TYPE source;
	STRESS_STRAIN *array;
	ARRAY_SORT_FLAG sort_flag;
} STRAIN_ARRAY;


typedef enum {
	SENS_SHAPE_GRAD_DENSITY,
	SENS_SHAPE_GRAD,
	SENS_DENSITY,
	SENS_THICKNESS,
	SENS_BEAM_XY,
	SENS_UNKNOWN
} SENSITIVITY_TYPE;


/* array of sensitivity */
typedef struct {
	int size;
	int alloc_size;
	SENSITIVITY_TYPE type;
	SENSITIVITY *array;
	ARRAY_SORT_FLAG sort_flag;
} SENSITIVITY_ARRAY;


/* array of default b.c. */
typedef struct {
	int size;
	int alloc_size;
	DATA_TYPE source;
	DEFAULT_BC *array;
	ARRAY_SORT_FLAG sort_flag;
} DEFAULT_BC_ARRAY;


/* array of sound pressure */
typedef struct {
	int size;
	int alloc_size;
	DATA_TYPE source;
	SOUND_PRESSURE *array;
	ARRAY_SORT_FLAG sort_flag;
} SOUND_PRESSURE_ARRAY;


/* array of element matrix */
typedef struct {
	int size;
	int alloc_size;
	DATA_TYPE source;
	ELEM_MATRIX *array;
	ARRAY_SORT_FLAG sort_flag;
} ELEM_MATRIX_ARRAY;


/* array of element lookup array */
typedef struct {
	int size;
	int alloc_size;
	DATA_TYPE source;
	ELEM_LOOKUP *array;
	ARRAY_SORT_FLAG sort_flag;
} ELEM_LOOKUP_ARRAY;


typedef struct {
	int size;
	int alloc_size;
	DATA_TYPE source;
	NODE_MATRIX *array;
	ARRAY_SORT_FLAG sort_flag;
} NODE_MATRIX_ARRAY;


/* array of load curve */
typedef struct {
	int size;
	int alloc_size;
	DATA_TYPE source;
	LOAD_CURVE *array;
	ARRAY_SORT_FLAG sort_flag;
} LOAD_CURVE_ARRAY;


typedef struct {
	int size;
	int alloc_size;
	DATA_TYPE source;
	LOCAL_COORD *array;
	ARRAY_SORT_FLAG sort_flag;
} LOCAL_COORD_ARRAY;


typedef enum {
	BC_NOT_COMMON,
	BC_COMMON
} BC_COMMON_FLAG;


/* array of b.c. set */
typedef struct {
	int size;
	int alloc_size;
	BC_COMMON_FLAG com_fix;
	BC_COMMON_FLAG com_mpc;
	BC_COMMON_FLAG com_force;
	BC_COMMON_FLAG com_temp;
	DATA_TYPE source;
	BC_SET *array;
	ARRAY_SORT_FLAG sort_flag;
	STRING_ARRAY s_info;
} BC_SET_ARRAY;


#define MAX_FACE       6
#define MAX_FACE_NODE  8
typedef struct {
	int face_num;
	int node_num[MAX_FACE];
	ELEM_TYPE face_type[MAX_FACE];
	int table[MAX_FACE][MAX_FACE_NODE];
} FACE_TABLE;


#define MAX_EDGE       12
#define MAX_EDGE_NODE  3
typedef struct {
	int edge_num;
	int node_num[MAX_EDGE];
	ELEM_TYPE edge_type[MAX_EDGE];
	int table[MAX_EDGE][MAX_EDGE_NODE];
} EDGE_TABLE;


/* contact condition of node */
typedef struct{
	int node;             /* label of target node (if < 0 then unused) */
	int size;             /* size of contact table */
	int alloc_size;       /* allocated size of table */
	CONTACT_TERM *term;   /* table of predict terminus point of node */
	ARRAY_SORT_FLAG sort_flag;
} CONTACT;

/* array of target nodes */
typedef struct{
	int size;
	int alloc_size;
	CONTACT *array;
	ARRAY_SORT_FLAG sort_flag;
} CONTACT_ARRAY;


typedef struct{
	int index;  /* index(element, node,..) */
	double min;
	double max;
} GEOM_SEGMENT;


/* 作業配列セット */
#define WS_SIZE    5
typedef struct {
	double *vec_N[WS_SIZE];
	double *vec_G[WS_SIZE];

	double **mat_3_3[WS_SIZE];
	double **mat_6_6[WS_SIZE];
	double **mat_9_9[WS_SIZE];

	double **mat_3_N[WS_SIZE];
	double **mat_3_3N[WS_SIZE];
	double **mat_6_3N[WS_SIZE];
	double **mat_9_3N[WS_SIZE];
	double **mat_3N_3N[WS_SIZE];

	double **mat_N_3[WS_SIZE];
	double **mat_N_4[WS_SIZE];
	double **mat_N_6[WS_SIZE];
	double **mat_G_4[WS_SIZE];
} WORK_SET;


/* fem_struct.c */
void init_node(NODE *node);
void init_element(ELEMENT *elem);
void init_rigid_element(RIGID_ELEMENT *elem);
void init_material_prop(MATERIAL_PROP *prop);
void init_physical_prop(PHYSICAL_PROP *prop);
void init_bc(BC *bc);
void init_elem_bc(ELEM_BC *bc);
void init_disp_velo(DISP_VELO *disp);
void init_stress_strain(STRESS_STRAIN *ss);
void init_sensitivity(SENSITIVITY *sens);
void init_default_bc(DEFAULT_BC *def);
void init_bc_set(BC_SET *bc_set);
void init_sound_pressure(SOUND_PRESSURE *sound);
void init_load_curve(LOAD_CURVE *curve);
void init_elem_matrix(ELEM_MATRIX *elem_matrix);
void init_elem_lookup(ELEM_LOOKUP *elem_lookup);
void init_node_matrix(NODE_MATRIX *nmat);
void init_local_coord(LOCAL_COORD *local_coord);
void init_contact_term(CONTACT_TERM *term);
void init_contact(CONTACT *contact);
void init_geom_segment(GEOM_SEGMENT *segment);
RC free_elem_matrix(ELEM_MATRIX *elem_matrix);
RC free_elem_lookup(ELEM_LOOKUP *elem_lookup);
RC free_node_matrix(NODE_MATRIX *nmat);
RC free_rigid_element(RIGID_ELEMENT *elem);
RC free_bc_set(BC_SET *bc_set);
RC free_contact(CONTACT *contact);
RC copy_rigid_element(RIGID_ELEMENT src, RIGID_ELEMENT *dest);
STRESS_ARRAY strain2stress_array(STRAIN_ARRAY strain);
STRAIN_ARRAY stress2strain_array(STRESS_ARRAY stress);
RC fill_material(MATERIAL_PROP *prop);
RC Gauss_const_default(ELEM_TYPE type, double **points,
                       double *weight, int *num_point);
RC Gauss_const_RI(ELEM_TYPE type, double **points,
                  double *weight, int *num_point);
RC Gauss_const_max(ELEM_TYPE type, double **points,
                   double *weight, int *num_point);
RC Gauss_const_1(ELEM_TYPE type, double **points, double *weight);
int element_ssnum(ELEM_TYPE type);
int element_dim(ELEM_TYPE type);
int element_order(ELEM_TYPE type);
RC trans_coord_s2v(ELEMENT s_elem, ELEMENT v_elem,
                   int coord_size, double **s_coords, double **v_coords);
RC allocate_work_set(WORK_SET *ws);
RC free_work_set(WORK_SET *ws);
int same_bc_set(int num1, const int id1[], double scale1[],
                int num2, const int id2[], double scale2[]);
RC cal_det_J(ELEMENT element, const double local_coord[],
             double **node_location, double *det_J,
             double **ws_mat_3_N, double **ws_mat_3_3);
RC cal_d_A(ELEMENT surf, const double local_coord[],
           double **node_location, VECT3D *d_A,
           double **ws_mat_3_N, double **ws_mat_3_3);
RC cal_dNxyz(ELEMENT element, const double local_coord[],
        double **node_location, double **dNxyz, double **dN, double **Jacobi,
        double **Jacobi_inv);
RC interpolate_disp(ELEMENT element, const double local_coord[],
                    double **disp_matrix, VECT3DR *disp);
RC realloc_node_matrix(NODE_MATRIX *nmat);


/* fem_struct_print.c */
void print_node(FILE *fp, NODE node);
void print_elem_type(FILE *fp, ELEM_TYPE type);
void print_element(FILE *fp, ELEMENT elem);
void print_rigid_element(FILE *fp, const RIGID_ELEMENT *elem);
void print_mat_type(FILE *fp, MAT_TYPE type);
void print_mat_state(FILE *fp, MAT_STATE state);
void print_analysis_type(FILE *fp, ANALYSIS_TYPE type);
void print_material_prop(FILE *fp, MATERIAL_PROP m_prop);
void print_section_type(FILE *fp, SECTION_TYPE type);
void print_physical_prop(FILE *fp, PHYSICAL_PROP p_prop);
void print_bc_type(FILE *fp, BC_TYPE type);
void print_bc_type_3d(FILE *fp, BC_TYPE3D type);
void print_bc(FILE *fp, BC bc);
void print_elem_bc(FILE *fp, ELEM_BC bc);
void print_disp_velo(FILE *fp, DISP_VELO d_velo);
void print_stress_strain(FILE *fp, STRESS_STRAIN ss);
void print_sensitivity(FILE *fp, SENSITIVITY sensi);
void print_sound_pressure(FILE *fp, SOUND_PRESSURE sound);
void print_load_curve(FILE *fp, const LOAD_CURVE *curve);
void print_elem_matrix(FILE *fp, ELEM_MATRIX elem_matrix);
void print_elem_lookup(FILE *fp, ELEM_LOOKUP elem_lookup);
void print_local_coord(FILE *fp, LOCAL_COORD coord);
void print_data_type(FILE *fp, DATA_TYPE type);
void print_array_sort_flag(FILE *fp, ARRAY_SORT_FLAG sort_flag);
void print_renumber_flag(FILE *fp, RENUMBER_FLAG renum_flag);
void print_contact_term(FILE *fp, CONTACT_TERM term);
void print_contact(FILE *fp, CONTACT contact);

void print_node_array(FILE *fp, NODE_ARRAY node);
void print_element_array(FILE *fp, ELEMENT_ARRAY element);
void print_rigid_element_array(FILE *fp, RIGID_ELEMENT_ARRAY element);
void print_material_prop_array(FILE *fp, MATERIAL_PROP_ARRAY m_prop);
void print_physical_prop_array(FILE *fp, PHYSICAL_PROP_ARRAY p_prop);
void print_bc_array(FILE *fp, BC_ARRAY bc);
void print_elem_bc_array(FILE *fp, ELEM_BC_ARRAY bc);
void print_disp_array(FILE *fp, DISP_ARRAY disp);
void print_velo_array(FILE *fp, VELO_ARRAY velo);
void print_react_array(FILE *fp, REACT_ARRAY react);
void print_stress_array(FILE *fp, STRESS_ARRAY stress);
void print_strain_array(FILE *fp, STRAIN_ARRAY strain);
void print_sensitivity_type(FILE *fp, SENSITIVITY_TYPE type);
void print_sensitivity_array(FILE *fp, SENSITIVITY_ARRAY sensi);
void print_default_bc(FILE *fp, DEFAULT_BC def);
void print_default_bc_array(FILE *fp, DEFAULT_BC_ARRAY def);
void print_bc_set(FILE *fp, BC_SET bc_set);
void print_bc_set_array(FILE *fp, BC_SET_ARRAY bc_set);
void print_bc_common_flag(FILE *fp, BC_COMMON_FLAG flag);
void print_sound_pressure_array(FILE *fp, SOUND_PRESSURE_ARRAY sound);
void print_load_curve_array(FILE *fp, LOAD_CURVE_ARRAY curve);
void print_elem_matrix_array(FILE *fp, ELEM_MATRIX_ARRAY elem_matrix);
void print_elem_lookup_array(FILE *fp, ELEM_LOOKUP_ARRAY elem_lookup);
void print_local_coord_array(FILE *fp, LOCAL_COORD_ARRAY coord);
void print_node_matrix(FILE *fp, NODE_MATRIX nmat);
void print_node_matrix_array(FILE *fp, NODE_MATRIX_ARRAY nmat);
void print_system_info(FILE *fp);
void print_contact_array(FILE *fp, CONTACT_ARRAY contact);


/* fem_struct_binio.c */
RC write_node_array(FILE *fp, NODE_ARRAY node);
RC write_element_array(FILE *fp, ELEMENT_ARRAY element);
RC write_material_prop_array(FILE *fp, MATERIAL_PROP_ARRAY m_prop);
RC write_physical_prop_array(FILE *fp, PHYSICAL_PROP_ARRAY p_prop);
RC write_bc_array(FILE *fp, BC_ARRAY bc);
RC write_elem_bc_array(FILE *fp, ELEM_BC_ARRAY bc);
RC write_disp_array(FILE *fp, DISP_ARRAY disp);
RC write_velo_array(FILE *fp, VELO_ARRAY velo);
RC write_react_array(FILE *fp, REACT_ARRAY react);
RC write_stress_array(FILE *fp, STRESS_ARRAY stress);
RC write_strain_array(FILE *fp, STRAIN_ARRAY strain);
RC write_sensitivity_array(FILE *fp, SENSITIVITY_ARRAY sensi);
RC write_default_bc_array(FILE *fp, DEFAULT_BC_ARRAY def);
RC write_sound_pressure_array(FILE *fp, SOUND_PRESSURE_ARRAY sound);
RC write_bc_set_array(FILE *fp, BC_SET_ARRAY bc_set);
RC write_load_curve_array(FILE *fp, LOAD_CURVE_ARRAY curve);
RC write_local_coord_array(FILE *fp, LOCAL_COORD_ARRAY coord);

RC read_node_array(FILE *fp, NODE_ARRAY *node);
RC read_element_array(FILE *fp, ELEMENT_ARRAY *element);
RC read_material_prop_array(FILE *fp, MATERIAL_PROP_ARRAY *m_prop);
RC read_physical_prop_array(FILE *fp, PHYSICAL_PROP_ARRAY *p_prop);
RC read_bc_array(FILE *fp, BC_ARRAY *bc);
RC read_elem_bc_array(FILE *fp, ELEM_BC_ARRAY *bc);
RC read_disp_array(FILE *fp, DISP_ARRAY *disp);
RC read_velo_array(FILE *fp, VELO_ARRAY *velo);
RC read_react_array(FILE *fp, REACT_ARRAY *react);
RC read_stress_array(FILE *fp, STRESS_ARRAY *stress);
RC read_strain_array(FILE *fp, STRAIN_ARRAY *strain);
RC read_sensitivity_array(FILE *fp, SENSITIVITY_ARRAY *sensi);
RC read_default_bc_array(FILE *fp, DEFAULT_BC_ARRAY *def);
RC read_sound_pressure_array(FILE *fp, SOUND_PRESSURE_ARRAY *sound);
RC read_bc_set_array(FILE *fp, BC_SET_ARRAY *bc_set);
RC read_load_curve_array(FILE *fp, LOAD_CURVE_ARRAY *curve);
RC read_local_coord_array(FILE *fp, LOCAL_COORD_ARRAY *coord);


/* fem_utl.c */
RC free_node_array(NODE_ARRAY *node);
RC free_element_array(ELEMENT_ARRAY *element);
RC free_rigid_element_array(RIGID_ELEMENT_ARRAY *element);
RC free_material_prop_array(MATERIAL_PROP_ARRAY *mat);
RC free_physical_prop_array(PHYSICAL_PROP_ARRAY *phys);
RC free_bc_array(BC_ARRAY *bc);
RC free_elem_bc_array(ELEM_BC_ARRAY *bc);
RC free_disp_array(DISP_ARRAY *disp);
RC free_velo_array(VELO_ARRAY *velo);
RC free_react_array(REACT_ARRAY *react);
RC free_stress_array(STRESS_ARRAY *stress);
RC free_strain_array(STRAIN_ARRAY *strain);
RC free_sensitivity_array(SENSITIVITY_ARRAY *sens);
RC free_default_bc_array(DEFAULT_BC_ARRAY *def);
RC free_sound_pressure_array(SOUND_PRESSURE_ARRAY *sound);
RC free_bc_set_array(BC_SET_ARRAY *bc_set);
RC free_load_curve_array(LOAD_CURVE_ARRAY *curve);
RC free_elem_matrix_array(ELEM_MATRIX_ARRAY *elem_matrix);
RC free_elem_lookup_array(ELEM_LOOKUP_ARRAY *elem_lookup);
RC free_node_matrix_array(NODE_MATRIX_ARRAY *node_matrix);
RC free_local_coord_array(LOCAL_COORD_ARRAY *coord);
RC free_contact_array(CONTACT_ARRAY *contact);

RC copy_node_array(NODE_ARRAY src, NODE_ARRAY *dest);
RC copy_element_array(ELEMENT_ARRAY src, ELEMENT_ARRAY *dest);
RC copy_material_prop_array(MATERIAL_PROP_ARRAY src,
                            MATERIAL_PROP_ARRAY *dest);
RC copy_physical_prop_array(PHYSICAL_PROP_ARRAY src,
                            PHYSICAL_PROP_ARRAY *dest);
RC copy_bc_array(BC_ARRAY src, BC_ARRAY *dest);
RC copy_elem_bc_array(ELEM_BC_ARRAY src, ELEM_BC_ARRAY *dest);
RC copy_disp_array(DISP_ARRAY src, DISP_ARRAY *dest);
RC copy_velo_array(VELO_ARRAY src, VELO_ARRAY *dest);
RC copy_react_array(REACT_ARRAY src, REACT_ARRAY *dest);
RC copy_stress_array(STRESS_ARRAY src, STRESS_ARRAY *dest);
RC copy_strain_array(STRAIN_ARRAY src, STRAIN_ARRAY *dest);
RC copy_sensitivity_array(SENSITIVITY_ARRAY src, SENSITIVITY_ARRAY *dest);
RC copy_default_bc_array(DEFAULT_BC_ARRAY src, DEFAULT_BC_ARRAY *dest);
RC copy_sound_pressure_array(SOUND_PRESSURE_ARRAY src,
                             SOUND_PRESSURE_ARRAY *dest);
RC copy_load_curve_array(LOAD_CURVE_ARRAY src, LOAD_CURVE_ARRAY *dest);
RC copy_local_coord_array(LOCAL_COORD_ARRAY src, LOCAL_COORD_ARRAY *dest);

RC clean_node_array(NODE_ARRAY *node);
RC clean_element_array(ELEMENT_ARRAY *elem);
RC clean_rigid_element_array(RIGID_ELEMENT_ARRAY *elem);
RC clean_material_prop_array(MATERIAL_PROP_ARRAY *mat);
RC clean_physical_prop_array(PHYSICAL_PROP_ARRAY *phys);
RC clean_bc_array(BC_ARRAY *bc);
RC clean_elem_bc_array(ELEM_BC_ARRAY *bc);
RC clean_disp_array(DISP_ARRAY *disp);
RC clean_velo_array(VELO_ARRAY *velo);
RC clean_react_array(REACT_ARRAY *react);
RC clean_stress_array(STRESS_ARRAY *stress);
RC clean_strain_array(STRAIN_ARRAY *strain);
RC clean_sensitivity_array(SENSITIVITY_ARRAY *sens);
RC clean_default_bc_array(DEFAULT_BC_ARRAY *def);
RC clean_sound_pressure_array(SOUND_PRESSURE_ARRAY *sound);
RC clean_bc_set_array(BC_SET_ARRAY *bc_set);
RC clean_elem_matrix_array(ELEM_MATRIX_ARRAY *elem_matrix);
RC clean_elem_lookup_array(ELEM_LOOKUP_ARRAY *elem_lookup);
RC clean_node_matrix_array(NODE_MATRIX_ARRAY *node_matrix);
RC clean_load_curve_array(LOAD_CURVE_ARRAY *curve);
RC clean_local_coord_array(LOCAL_COORD_ARRAY *coord);
RC clean_contact_array(CONTACT_ARRAY *contact);

RC realloc_node_array(NODE_ARRAY *node);
RC realloc_element_array(ELEMENT_ARRAY *elem);
RC realloc_rigid_element_array(RIGID_ELEMENT_ARRAY *elem);
RC realloc_material_prop_array(MATERIAL_PROP_ARRAY *mat);
RC realloc_physical_prop_array(PHYSICAL_PROP_ARRAY *phys);
RC realloc_bc_array(BC_ARRAY *bc);
RC realloc_elem_bc_array(ELEM_BC_ARRAY *bc);
RC realloc_disp_array(DISP_ARRAY *disp);
RC realloc_velo_array(VELO_ARRAY *velo);
RC realloc_react_array(REACT_ARRAY *react);
RC realloc_stress_array(STRESS_ARRAY *stress);
RC realloc_strain_array(STRAIN_ARRAY *strain);
RC realloc_sensitivity_array(SENSITIVITY_ARRAY *sens);
RC realloc_default_bc_array(DEFAULT_BC_ARRAY *def);
RC realloc_sound_pressure_array(SOUND_PRESSURE_ARRAY *sound);
RC realloc_bc_set_array(BC_SET_ARRAY *bc_set);
RC realloc_elem_matrix_array(ELEM_MATRIX_ARRAY *elem_matrix);
RC realloc_elem_lookup_array(ELEM_LOOKUP_ARRAY *elem_lookup);
RC realloc_node_matrix_array(NODE_MATRIX_ARRAY *node_matrix);
RC realloc_load_curve_array(LOAD_CURVE_ARRAY *curve);
RC realloc_local_coord_array(LOCAL_COORD_ARRAY *local);
RC realloc_contact_array(CONTACT_ARRAY *contact);

RC allocate_node_array(int num, NODE_ARRAY *node);
RC allocate_element_array(int num, ELEMENT_ARRAY *element);
RC allocate_rigid_element_array(int num, RIGID_ELEMENT_ARRAY *element);
RC allocate_material_prop_array(int num, MATERIAL_PROP_ARRAY *mat);
RC allocate_physical_prop_array(int num, PHYSICAL_PROP_ARRAY *phys);
RC allocate_bc_array(int num, BC_ARRAY *bc);
RC allocate_elem_bc_array(int num, ELEM_BC_ARRAY *bc);
RC allocate_disp_array(int num, DISP_ARRAY *disp);
RC allocate_velo_array(int num, VELO_ARRAY *velo);
RC allocate_react_array(int num, REACT_ARRAY *react);
RC allocate_stress_array(int num, STRESS_ARRAY *stress);
RC allocate_strain_array(int num, STRAIN_ARRAY *strain);
RC allocate_sensitivity_array(int num, SENSITIVITY_ARRAY *sens);
RC allocate_default_bc_array(int num, DEFAULT_BC_ARRAY *def);
RC allocate_sound_pressure_array(int num, SOUND_PRESSURE_ARRAY *sound);
RC allocate_bc_set_array(int num, BC_SET_ARRAY *bc_set);
RC allocate_elem_matrix_array(int num, ELEM_MATRIX_ARRAY *elem_matrix);
RC allocate_elem_lookup_array(int num, ELEM_LOOKUP_ARRAY *elem_lookup);
RC allocate_node_matrix_array(int num, NODE_MATRIX_ARRAY *node_matrix);
RC allocate_load_curve_array(int num, LOAD_CURVE_ARRAY *curve);
RC allocate_local_coord_array(int num, LOCAL_COORD_ARRAY *coord);
RC allocate_contact_array(int num, CONTACT_ARRAY *contact);

int count_valid_node(NODE_ARRAY node);
int count_valid_element(ELEMENT_ARRAY elem);
int count_valid_rigid_element(RIGID_ELEMENT_ARRAY elem);
int count_valid_material_prop(MATERIAL_PROP_ARRAY material);
int count_valid_physical_prop(PHYSICAL_PROP_ARRAY phys);
int count_valid_bc(BC_ARRAY bc);
int count_valid_elem_bc(ELEM_BC_ARRAY bc);
int count_valid_disp(DISP_ARRAY disp);
int count_valid_velo(VELO_ARRAY velo);
int count_valid_react(REACT_ARRAY react);
int count_valid_load_curve(LOAD_CURVE_ARRAY curve);
int count_valid_local_coord(LOCAL_COORD_ARRAY coord);
int count_valid_contact(CONTACT_ARRAY contact);

RC sort_node_array(NODE_ARRAY *node);
RC sort_node_renum(NODE_ARRAY *node);
RC sort_element_array(ELEMENT_ARRAY *elem);
RC sort_rigid_element_array(RIGID_ELEMENT_ARRAY *elem);
RC sort_material_prop_array(MATERIAL_PROP_ARRAY *mat);
RC sort_physical_prop_array(PHYSICAL_PROP_ARRAY *phys);
RC sort_bc_array(BC_ARRAY *bc);
RC sort_elem_bc_array(ELEM_BC_ARRAY *bc);
RC sort_disp_array(DISP_ARRAY *disp);
RC sort_velo_array(VELO_ARRAY *velo);
RC sort_react_array(REACT_ARRAY *react);
RC sort_stress_array(STRESS_ARRAY *stress);
RC sort_strain_array(STRAIN_ARRAY *strain);
RC sort_sensitivity_array(SENSITIVITY_ARRAY *sens);
RC sort_default_bc_array(DEFAULT_BC_ARRAY *def);
RC sort_sound_pressure_array(SOUND_PRESSURE_ARRAY *sound);
RC sort_bc_set_array(BC_SET_ARRAY *bc_set);
RC sort_elem_matrix_array(ELEM_MATRIX_ARRAY *elem_matrix);
RC sort_elem_lookup_array(ELEM_LOOKUP_ARRAY *elem_lookup);
RC sort_node_matrix_array(NODE_MATRIX_ARRAY *node_matrix);
RC sort_load_curve_array(LOAD_CURVE_ARRAY *curve);
RC sort_local_coord_array(LOCAL_COORD_ARRAY *coord);
RC sort_contact_array(CONTACT_ARRAY *contact);

int search_node_label(NODE_ARRAY node, int label);
int search_element_label(ELEMENT_ARRAY elem, int label);
int search_rigid_element_label(RIGID_ELEMENT_ARRAY elem, int label);
int search_material_prop_label(MATERIAL_PROP_ARRAY mat, int label);
int search_physical_prop_label(PHYSICAL_PROP_ARRAY phys, int label);
int search_bc_node(BC_ARRAY bc, int node);
int search_elem_bc_element(ELEM_BC_ARRAY bc, int element);
int search_disp_node(DISP_ARRAY disp, int node);
int search_velo_node(VELO_ARRAY velo, int node);
int search_react_node(REACT_ARRAY react, int node);
int search_stress_element(STRESS_ARRAY stress, int element);
int search_strain_element(STRAIN_ARRAY strain, int element);
int search_sensitivity_element(SENSITIVITY_ARRAY sens, int element);
int search_default_bc_set_id(DEFAULT_BC_ARRAY def, int set_id);
int search_sound_pressure_node(SOUND_PRESSURE_ARRAY sound, int node);
int search_bc_set_label(BC_SET_ARRAY bc_set, int label);
int search_elem_matrix_element(ELEM_MATRIX_ARRAY elem_matrix, int element);
int search_elem_lookup_node(ELEM_LOOKUP_ARRAY elem_lookup, int node);
int search_node_matrix_node(NODE_MATRIX_ARRAY node_matrix, int node);
int search_load_curve_label(LOAD_CURVE_ARRAY curve, int label);
int search_local_coord_label(LOCAL_COORD_ARRAY coord, int label);
int search_material_prop_element(MATERIAL_PROP_ARRAY mat,
                                 PHYSICAL_PROP_ARRAY phys, ELEMENT *element);
int search_stress_element_node(STRESS_ARRAY stress, int element, int node);
int search_strain_element_node(STRAIN_ARRAY stress, int element, int node);
int search_sensitivity_element_node(SENSITIVITY_ARRAY sens,
                                    int element, int node);
int search_contact_node(CONTACT_ARRAY contact, int node);

RC add_disp_array(double weight, DISP_ARRAY disp, DISP_ARRAY *disp_total);
RC add_sensitivity_array(double weight, SENSITIVITY_ARRAY sens,
                         SENSITIVITY_ARRAY *sens_total);
RC node_location_matrix(ELEMENT *element, NODE_ARRAY node, double **matrix);
RC local_node_location(double L[][3], int node_num, double **node_location);
RC global_node_location(double L[][3], int node_num, double **node_location);
RC fill_disp_matrix(ELEMENT element, DISP_ARRAY disp, double **matrix);
RC fill_sound_pres_matrix(ELEMENT element,
                          SOUND_PRESSURE_ARRAY pres, double **matrix);
RC unit_normal_vect_center(ELEMENT surf, ELEMENT inner_elem,
                           NODE_ARRAY node, VECT3D *vect);
RC unit_normal_vect(double *local_coord, ELEMENT surf, 
                    ELEMENT inner_elem, NODE_ARRAY node, VECT3D *vect);
RC extract_bc_set(int num, const int id[], double scale[],
                  BC_ARRAY src, BC_ARRAY *dest);
RC extract_default_bc_set(int num, const int id[], double scale[],
                          DEFAULT_BC_ARRAY src, DEFAULT_BC_ARRAY *dest);
RC extract_rigid_set(int num, const int id[], double scale[],
                     RIGID_ELEMENT_ARRAY src, RIGID_ELEMENT_ARRAY *dest);
RC extract_elem_bc_set(int num, const int id[], double scale[],
                       ELEM_BC_ARRAY src, ELEM_BC_ARRAY *dest);
RC compress_bc_array(BC_ARRAY *dest);
RC mul_elem_matrices_vector(int total_dof, ELEM_MATRIX_ARRAY elem_matrix,
                            const double *vector, double *answer);
RC make_elem_lookup(NODE_ARRAY node, ELEMENT_ARRAY element,
                    ELEM_LOOKUP_ARRAY *elem_lookup);
RC make_elem_lookup_rigid(NODE_ARRAY node, ELEMENT_ARRAY element,
                          RIGID_ELEMENT_ARRAY rigid_elem,
                          ELEM_LOOKUP_ARRAY *elem_lookup);
int analysis_dim(ELEMENT_ARRAY element);
int analysis_order(ELEMENT_ARRAY element);
RC reduced_element_array_order(ELEMENT_ARRAY *elem);
RC reduced_element_type(ELEM_TYPE old_type, ELEM_TYPE *new_type,
                        int *new_node_num);
RC delete_duplicate_node(ELEMENT_ARRAY element, NODE_ARRAY *node,
                         double cr_dist);
RC minimum_node_distance(ELEMENT_ARRAY element, NODE_ARRAY node,
                         double *dist);
RC delete_unused_node(ELEMENT_ARRAY elem, NODE_ARRAY *node);
RC delete_unused_node_r(ELEMENT_ARRAY elem, RIGID_ELEMENT_ARRAY rigid,
                        NODE_ARRAY *node);
RC delete_unused_bc(NODE_ARRAY node, BC_ARRAY *bc);
RC locate_middle_node(NODE_ARRAY node, ELEMENT_ARRAY elem);
RC interpolate_table(ELEM_TYPE type, int table[][2]);
RC element_average_length(ELEMENT element, NODE_ARRAY node, double *length);
double element_length(int node_num, double **node_location);
RC element_edge_table(ELEM_TYPE type, int *edge_num, int table[][2]);
RC elem_matrix_cg(int total_dof, ELEM_MATRIX_ARRAY elem_matrix,
                  const double *vector, double *solution);
RC get_element_construct_node(int *num, int **node_label, ELEMENT_ARRAY elem);


/* fem_search_common.c */
RC search_common_node_label(NODE_ARRAY node1, NODE_ARRAY node2,
                            NODE_ARRAY *common);

RC search_common_node_label_element(NODE_ARRAY node1, NODE_ARRAY node2,
                                    ELEMENT_ARRAY elem1, ELEMENT_ARRAY elem2,
                                    NODE_ARRAY *common);


/* fem_renumber.c */
RC dummy_renumber(NODE_ARRAY *node);
RC fem_renumber0(NODE_ARRAY *node, ELEMENT_ARRAY element,
                 RIGID_ELEMENT_ARRAY rigid, int type_flag);
int node_renum_index(NODE_ARRAY node, int node_label);
int node_renum_index_cache(NODE_ARRAY node, int node_label,
                           int cache_size, int cache[][2]);


/* element_volume.c */
RC total_element_volume(ELEMENT_ARRAY element,
                        NODE_ARRAY node, double *total_volume);
RC element_volume(ELEMENT *element, NODE_ARRAY node);


/* extract_surface.c */
RC extract_surface(ELEMENT_ARRAY elem, ELEMENT_ARRAY *surf);
RC extract_divider(ELEMENT_ARRAY elem, ELEMENT_ARRAY *div);
RC cut_element(NODE_ARRAY *node, ELEMENT_ARRAY *elem);


/* face_edge_table.c */
void init_face_table(FACE_TABLE *f_table);
RC make_face_table(ELEM_TYPE type, FACE_TABLE *f_table);
void init_edge_table(EDGE_TABLE *e_table);
RC make_edge_table(ELEM_TYPE type, EDGE_TABLE *e_table);


/* interpolation.c */
RC set_local_node_points(ELEM_TYPE type, double **coords);
RC set_center_point(ELEM_TYPE type, double *coord);
RC set_N_vector(ELEM_TYPE type, const double *coord, double *N);
RC set_dN_matrix(ELEM_TYPE type, const double *coord, double **dN);
RC set_d2N_matrix(ELEM_TYPE type, const double *coord, double **d2N);


/* integration.c */
RC Gauss_const_line(int number_of_point, double **xi, double *weight);
RC Gauss_const_tri(int number_of_point, double **L1_L2_L3, double *weight);
RC Gauss_const_quad(int number_of_point, double **xi_eta, double *weight);
RC Gauss_const_tetra(int number_of_point, double **L1_L2_L3_L4, double *weight);
RC Gauss_const_penta(int number_of_point,double **L1_L2_L3_zeta,double *weight);
RC Gauss_const_hexa(int number_of_point, double **xi_eta_zeta, double *weight);


/* pressure_force.c */
RC pressure_force(int p_size, const double *p, ELEMENT *element,
                  NODE_ARRAY node, VECT3D *f);
RC traction_force(int p_size, const double *p, VECT3D dir,
                  ELEMENT *element, NODE_ARRAY node, VECT3D *f);


/* contact.c */
RC make_contact_terminus(ELEMENT_ARRAY element, NODE_ARRAY node,
                         NODE point, double range, int seg_size,
                         GEOM_SEGMENT *segment[], CONTACT *contact);


/* fem_geometry.c */
RC element_bbox(ELEMENT element, NODE_ARRAY node, VECT3D *min, VECT3D *max);
RC make_geom_search_array(ELEMENT_ARRAY elem, NODE_ARRAY node,
                          int *size, GEOM_SEGMENT *segment[]);
RC make_geom_search_array_2(ELEMENT_ARRAY elem, NODE_ARRAY node,
                            int index_size, int *index_array,
                            int *size, GEOM_SEGMENT *segment[]);
RC geom_search_candidate(int size, GEOM_SEGMENT **segment,            
                         VECT3D p, double range, int *cand_size, int cand[]);
RC extract_geom_contact_surface(ELEMENT_ARRAY surf1, NODE_ARRAY node1,
                                ELEMENT_ARRAY surf2, NODE_ARRAY node2,
                                int *size1, int index_array1[],
                                int *size2, int index_array2[]);
RC print_geom_segment_array(FILE *fp, int size, GEOM_SEGMENT segment[]);


#endif /* FEM_STRUCT_H */

