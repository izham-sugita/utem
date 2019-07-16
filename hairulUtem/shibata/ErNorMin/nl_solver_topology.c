#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "base.h"
#include "fem.h"
#include "mathematics.h"
#include "ccs_format.h"
#include "dmumps_c.h"
#include "topology_opt.h"

#define NR_MAX (100000)
#define PRI_NODE (46943)      /* 荷重・変位曲線出力する節点番号 */

typedef struct{
	int num_g_point;
	double **g_points;
	double weight[MAX_GAUSS_POINT+1];
	int B_matrix_row;
	int B_matrix_col;
	double **B_matrix[MAX_GAUSS_POINT+1];
	double det_J[MAX_GAUSS_POINT+1];
} G_PROP;


static RC set_rest_zero(NODE_ARRAY node, BC_ARRAY rest, double resid_v[]);
static RC allocate_g_prop (G_PROP *g_prop);
static RC set_Gauss_const (ELEM_TYPE type, G_PROP *g_prop);
static RC set_g_prop (ELEMENT element, double **node_location, G_PROP *g_prop);

/* 位相最適化用微小歪み大変形解析 Total Lagrange法ソルバー */
RC
fem_NL_solve_mumps_topology (NODE_ARRAY node, ELEMENT_ARRAY element,
                             MATERIAL_PROP_ARRAY material,
                             PHYSICAL_PROP_ARRAY physical,
                             BC_ARRAY rest, BC_ARRAY force, BC_ARRAY density,
                             DEFAULT_BC_ARRAY b_force,
                             BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                             DISP_ARRAY *disp, STRESS_ARRAY *stress,
                             STRAIN_ARRAY *strain, double dens_pow,
                             int sol_ss_flag)
{
	CCS ccs;
	NONZERO_MATRIX_CCS tan_K;
	int ii1, ii2;
	int dim, total_dof;
	double beta1 = 1.0e-3;
	double beta2 = 1.0e-7;
	double *f_vect;
	double *resid_v;
	double *react_v;
	double *delta_d;

	double omp_begin;
	double omp_end;

	/* MUMPS */
	DMUMPS_STRUC_C id;


	dim = analysis_dim(element);
	if(dim != 3) return(ARG_ERROR_RC);

	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

	/* 全自由度 */
	RC_TRY( fem_total_dof(node, element, &total_dof) );

	RC_TRY( allocate1D(total_dof, &f_vect) );
	RC_TRY( allocate1D(total_dof, &resid_v) );
	RC_TRY( allocate1D(total_dof, &react_v) );
	RC_TRY( allocate1D(total_dof, &delta_d) );
	RC_TRY( vect2disp(dim, node, delta_d, disp) );

	/* 節点荷重，重力入力 */
	RC_TRY( add_nodal_force(dim, node, force, f_vect) );
	RC_TRY( add_body_force(dim, node, element, material, physical, b_force,
                           f_vect) );
	/* 熱ひずみ(初期ひずみ)による荷重入力 */
	RC_TRY( add_temp_force(dim, node, element, material, physical, temp,
                           def_temp, f_vect) );

	for(ii1=0; ii1<total_dof; ii1++){
		resid_v[ii1] = f_vect[ii1];
	}
	for(ii1=0; ii1<NR_MAX; ii1++){
		double du_norm, du_r, u_norm, u_f;

		RC_TRY( log_printf(5, " +-------------------------+\n") );
		RC_TRY( log_printf(5, " | N-R iteration NUM : %3d |\n", ii1+1) );
		RC_TRY( log_printf(5, " +-------------------------+\n") );
		if(ii1 == 0){
			RC_TRY( make_nonzero1a_ccs_sym_topology(node, element, material,
                                                    physical, density,
													dens_pow, &tan_K) );
		}else{
			RC_TRY( make_NL_nonzero1a_ccs_sym_topology(node, element, material,
                                                       physical, *disp,
													   density, dens_pow, 
													   &tan_K) );
		}
		RC_TRY( copy_nonzero_matrix1a_ccs2ccs_format(&tan_K, &ccs) );
		RC_TRY( free_nonzero_matrix1a_ccs(&tan_K) );
		RC_TRY( modify_ccs_format(dim, node, rest, &ccs, resid_v) );
		
		omp_begin = omp_get_wtime();
		/* MUMPS */
		/* Inititalization */
		id.job = JOB_INIT;
		id.par = 1;
		/* Symmetric or Non-symmetric parameter */
		id.sym = 1;
//  id.comm_fortran=USE_COMM_WORLD;
		dmumps_c(&id);
		/* Set Matrix data */
		id.n = total_dof;
		id.nz = ccs.array_size;
		id.irn = ccs.row_index;
		id.jcn = ccs.col_index;
		for(ii2=0; ii2<id.nz; ii2++){
			id.irn[ii2]++;
			id.jcn[ii2]++;
		}
		id.a = ccs.array;
		id.rhs = resid_v;
		/* No outputs */
		id.icntl[0]=-1; id.icntl[1]=-1; id.icntl[2]=-1; id.icntl[3]=0;
		id.icntl[21]=0;
		/* Symbolic + Factorize + Solve */
		id.job = 6;
		dmumps_c(&id);
		/* End Job */
		id.job = JOB_END;
		dmumps_c(&id);
		omp_end = omp_get_wtime();
		printf("MUMPS OMP Solve  OMPTIME %5.3lf sec\n", (omp_end-omp_begin));

		for(ii2=0; ii2<total_dof; ii2++){
			delta_d[ii2] = resid_v[ii2];
		}
		
		RC_TRY( free_nonzero_ccs(&ccs) );
		RC_TRY( add_vect2disp(dim, node, delta_d, *disp, 1.0) );

		/* 反力計算 */
		RC_TRY( make_NL_react_force_topology(element, node, material, physical,
                                             density, dens_pow,
		                                     *disp, react_v) );
		/* 残差荷重ベクトル */
		for(ii2=0; ii2<total_dof; ii2++){
			resid_v[ii2] = f_vect[ii2] - react_v[ii2];
		}
		RC_TRY( set_rest_zero(node, rest, resid_v) );

		/* 変位，エネルギー収束判定 */
		du_norm = sqrt( inner_product(total_dof, delta_d, delta_d) );
		du_r = 0.0;
		for(ii2=0; ii2<total_dof; ii2++){
			du_r += fabs(delta_d[ii2] * resid_v[ii2]);
		}
		RC_TRY( copy_disp2vect(dim, node, *disp, delta_d, 1.0) );
		u_norm = sqrt( inner_product(total_dof, delta_d, delta_d) );
		u_f = 0.0;
		for(ii2=0; ii2<total_dof; ii2++){
			u_f += fabs(delta_d[ii2] * f_vect[ii2]);
		}
		RC_TRY( log_printf(5, "    Disp   error :%15.7e\n", du_norm/u_norm) );
		RC_TRY( log_printf(5, "    Energy error :%15.7e\n", du_r/u_f) );
		/* 収束判定 */
		if(du_norm < beta1*u_norm + ABS_TOL
           || du_r < beta2*u_f + ABS_TOL){
			break;
		}
		RC_TRY( log_printf(5, "\n") );
	}
	if(ii1 == NR_MAX){
		log_printf(3, "N-R iteration NOT converged!!\n");
	}
	/* メモリ開放 */
	RC_TRY( mm_free(f_vect) );
	RC_TRY( mm_free(resid_v) );
	RC_TRY( mm_free(react_v) );
	RC_TRY( mm_free(delta_d) );

	/* 第2Piola-Kirchhoff応力，Green歪み(応力解析は位相最適化用にしていない) */
	RC_TRY( PK2_stress_Green_strain(sol_ss_flag, node, element, *disp,
                                    material, physical, stress, strain) );

	return(NORMAL_RC);
}

/* 密度を考慮した線形剛性マトリックスの計算 */
RC
make_nonzero1a_ccs_sym_topology (NODE_ARRAY node, ELEMENT_ARRAY element,
                                 MATERIAL_PROP_ARRAY material,
                                 PHYSICAL_PROP_ARRAY physical,
							     BC_ARRAY density, double dens_pow,
                                 NONZERO_MATRIX_CCS *matrix)
{
	int ii1, ii2, ii3;
	int dim;
	double *ptr;
	ELEM_MATRIX elem_matrix;


	dim = analysis_dim(element);

	RC_TRY( allocate_nonzero_matrix1a_ccs(node.size*dim, matrix) );
	for(ii1=0; ii1<element.size; ii1++){

		if(element.array[ii1].label < 0) continue;

		RC_TRY( make_elem_matrix_topology(element.array[ii1], node,
                                          material, physical,
		                                  density, dens_pow,
		                                  &elem_matrix) );

		for(ii2=0; ii2<elem_matrix.size; ii2++){
			for(ii3=0; ii3<=ii2; ii3++){
				ptr = nonzero_matrix1a_ccs_ptr(matrix, elem_matrix.index[ii2],
                                               elem_matrix.index[ii3]);
				(*ptr) += elem_matrix.matrix[ii2][ii3];
			}
		}
		RC_TRY( free_elem_matrix(&elem_matrix) );
	}

	return(NORMAL_RC);
}

/* 密度を考慮した非線形剛性マトリックスの計算 */
RC
make_NL_nonzero1a_ccs_sym_topology (NODE_ARRAY node, ELEMENT_ARRAY element,
                                    MATERIAL_PROP_ARRAY material,
                                    PHYSICAL_PROP_ARRAY physical,
                                    DISP_ARRAY disp,
								    BC_ARRAY density, double dens_pow,
                                    NONZERO_MATRIX_CCS *K)
{
	int ii1, ii2, ii3;
	int dim;
	double *ptr;
	WORK_SET ws;
	ELEM_MATRIX elem_matrix;

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	dim = analysis_dim(element);
	if(dim != 3) return(ARG_ERROR_RC);

	RC_TRY( allocate_nonzero_matrix1a_ccs(node.size*dim, K) );
	RC_TRY( allocate_work_set(&ws) );

	for(ii1=0; ii1<element.size; ii1++){

		if(element.array[ii1].label < 0) continue;

		RC_TRY( make_tangent_elem_matrix_topology(element.array[ii1], node,
                                                  material, physical, ws, disp,
		                                          density, dens_pow,
                                                  &elem_matrix) );
		for(ii2=0; ii2<elem_matrix.size; ii2++){
			for(ii3=0; ii3<=ii2; ii3++){
				ptr = nonzero_matrix1a_ccs_ptr(K, elem_matrix.index[ii2],
                                               elem_matrix.index[ii3]);
				(*ptr) += elem_matrix.matrix[ii2][ii3];
			}
		}

		RC_TRY( free_elem_matrix(&elem_matrix) );
	}

	RC_TRY( free_work_set(&ws) );

	return(NORMAL_RC);
}

#if 0
/* robin型ver. with Load Control Method (Modified Newton-Raphson) */
RC
fem_NL_solve_mumps_topology_robin_modLCM (NODE_ARRAY node,
                                       ELEMENT_ARRAY element,
                                       ELEMENT_ARRAY surf,
                                       MATERIAL_PROP_ARRAY material,
                                       PHYSICAL_PROP_ARRAY physical,
                                       BC_ARRAY rest, BC_ARRAY force,
								       BC_ARRAY density, BC_ARRAY weight,
                                       DEFAULT_BC_ARRAY b_force,
                                       BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                                       DISP_ARRAY *disp, STRESS_ARRAY *stress,
                                       STRAIN_ARRAY *strain, double dens_pow,
                                       int sol_ss_flag, int lcm_num)
{
	CCS ccs;
	NONZERO_MATRIX_CCS tan_K;
	int ii0, ii1, ii2;
	int dim, total_dof;
	double beta1 = 1.0e-3;
	double beta2 = 1.0e-6;
	double *f_vect;
	double *df_vect;
	double *resid_v;
	double *react_v;
	double *delta_d;
	double *d0_vect;

	double omp_begin;
	double omp_end;

	/* MUMPS */
	DMUMPS_STRUC_C id;


	dim = analysis_dim(element);
	if(dim != 3) return(ARG_ERROR_RC);

	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

	/* 全自由度 */
	RC_TRY( fem_total_dof(node, element, &total_dof) );

	RC_TRY( allocate1D(total_dof, &f_vect) );
	RC_TRY( allocate1D(total_dof, &df_vect) );
	RC_TRY( allocate1D(total_dof, &resid_v) );
	RC_TRY( allocate1D(total_dof, &react_v) );
	RC_TRY( allocate1D(total_dof, &delta_d) );
	RC_TRY( allocate1D(total_dof, &d0_vect) );
	RC_TRY( vect2disp(dim, node, delta_d, disp) );

	/* 節点荷重，重力入力 */
	RC_TRY( add_nodal_force(dim, node, force, f_vect) );
	RC_TRY( add_body_force(dim, node, element, material, physical, b_force,
                           f_vect) );
	/* 熱ひずみ(初期ひずみ)による荷重入力 */
	RC_TRY( add_temp_force(dim, node, element, material, physical, temp,
                           def_temp, f_vect) );

	for(ii1=0; ii1<total_dof; ii1++){
		df_vect[ii1] = f_vect[ii1] / (double)lcm_num;
		f_vect[ii1] = 0.0;
	}
	for(ii0=0; ii0<lcm_num; ii0++){
		for(ii1=0; ii1<total_dof; ii1++){
			f_vect[ii1] += df_vect[ii1];
			resid_v[ii1] = f_vect[ii1] - react_v[ii1];
		}

		RC_TRY( make_NL_nonzero1a_ccs_sym_topology_robin(node, element,
		                                            surf, material,
		                                            physical, *disp,
		                                            density, weight,
		                                            dens_pow, &tan_K));
		RC_TRY( copy_nonzero_matrix1a_ccs2ccs_format(&tan_K, &ccs) );
		RC_TRY( free_nonzero_matrix1a_ccs(&tan_K) );
		RC_TRY( modify_ccs_format(dim, node, rest, &ccs, resid_v) );

		omp_begin = omp_get_wtime();
		/* MUMPS */
		/* Inititalization */
		id.job = JOB_INIT;
		id.par = 1;
		/* Symmetric or Non-symmetric parameter */
		id.sym = 1;
		//  id.comm_fortran=USE_COMM_WORLD;
		dmumps_c(&id);
		/* Set Matrix data */
		id.n = total_dof;
		id.nz = ccs.array_size;
		id.irn = ccs.row_index;
		id.jcn = ccs.col_index;
		for(ii2=0; ii2<id.nz; ii2++){
			id.irn[ii2]++;
			id.jcn[ii2]++;
		}
		id.a = ccs.array;
		/* No outputs */
		id.icntl[0]=-1; id.icntl[1]=-1; id.icntl[2]=-1; id.icntl[3]=0;
		id.icntl[21]=0;
		/* Symbolic + Factorize */
		id.job = 4;
		dmumps_c(&id);
		omp_end = omp_get_wtime();
		printf("MUMPS OMP DECOMP  OMPTIME %5.3lf sec\n",
		       (omp_end-omp_begin));

		for(ii1=0; ii1<NR_MAX; ii1++){
			double du_norm, du_r, u_norm, u_f;
			double u0_norm;

			RC_TRY( log_printf(5, " +-------------------------+\n") );
			RC_TRY( log_printf(5, " |--------   %3d   --------|\n", ii0+1) );
			RC_TRY( log_printf(5, " +-------------------------+\n") );
			RC_TRY( log_printf(5, " | N-R iteration NUM : %3d |\n", ii1+1) );
			RC_TRY( log_printf(5, " +-------------------------+\n") );

			/* Solve */
			id.rhs = resid_v;
			id.job = 3;
			dmumps_c(&id);

			for(ii2=0; ii2<total_dof; ii2++){
				delta_d[ii2] = resid_v[ii2];
			}

			if(ii1 == 0){
				for(ii2=0; ii2<total_dof; ii2++){
					d0_vect[ii2] = delta_d[ii2];
				}
			}

			RC_TRY( add_vect2disp(dim, node, delta_d, *disp, 1.0) );

			/* 反力計算 */
			RC_TRY( make_NL_react_force_topology(element, node, material,
			                                     physical, density, dens_pow,
			                                     *disp, react_v) );
			/* 残差荷重ベクトル */
			for(ii2=0; ii2<total_dof; ii2++){
				resid_v[ii2] = f_vect[ii2] - react_v[ii2];
			}
			RC_TRY( set_rest_zero(node, rest, resid_v) );

			/* 変位，エネルギー収束判定 */
			du_norm = sqrt( inner_product(total_dof, delta_d, delta_d) );
			du_r = 0.0;
			for(ii2=0; ii2<total_dof; ii2++){
				du_r += fabs(delta_d[ii2] * resid_v[ii2]);
			}
			RC_TRY( copy_disp2vect(dim, node, *disp, delta_d, 1.0) );
			u_norm = sqrt( inner_product(total_dof, delta_d, delta_d) );
			u0_norm = sqrt( inner_product(total_dof, d0_vect, d0_vect) );
			u_f = 0.0;
			for(ii2=0; ii2<total_dof; ii2++){
				u_f += fabs(delta_d[ii2] * f_vect[ii2]);
			}
			RC_TRY( log_printf(5, " Disp   error :%15.7e\n", du_norm/u0_norm) );
			RC_TRY( log_printf(5, " Energy error :%15.7e\n", du_r/u_f) );
			/* 収束判定 */
			if(du_norm < beta1*u0_norm + ABS_TOL
			&& du_r < beta2*u_f + ABS_TOL){
				break;
			}
			RC_TRY( log_printf(5, "\n") );

#if 1
			/* Modify for CCS */
			for(ii2=0; ii2<ccs.array_size; ii2++){
				ccs.row_index[ii2]--;
				ccs.col_index[ii2]--;
			}
			RC_TRY( modify_ccs_format(dim, node, rest, &ccs, resid_v) );
			for(ii2=0; ii2<ccs.array_size; ii2++){
				ccs.row_index[ii2]++;
				ccs.col_index[ii2]++;
			}
#endif
		}
		/* End Job */
		id.job = JOB_END;
		dmumps_c(&id);

		RC_TRY( free_nonzero_ccs(&ccs) );

		if(ii1 == NR_MAX){
			log_printf(3, "N-R iteration NOT converged!!\n");
		}
		RC_TRY( log_printf(5, "---\n") );
		RC_TRY( log_printf(5, "\n") );
	}
	/* メモリ開放 */
	RC_TRY( mm_free(f_vect) );
	RC_TRY( mm_free(df_vect) );
	RC_TRY( mm_free(resid_v) );
	RC_TRY( mm_free(react_v) );
	RC_TRY( mm_free(delta_d) );
	RC_TRY( mm_free(d0_vect) );

	/* 第2Piola-Kirchhoff応力，Green歪み(応力解析は位相最適化用にしていない) */
	RC_TRY( PK2_stress_Green_strain(sol_ss_flag, node, element, *disp,
                                    material, physical, stress, strain) );

	return(NORMAL_RC);
}
#endif

#if 1
/* robin型ver. with Load Control Method*/
RC
fem_NL_solve_mumps_topology_robin_LCM (NODE_ARRAY node, ELEMENT_ARRAY element,
                                       ELEMENT_ARRAY surf,
                                       MATERIAL_PROP_ARRAY material,
                                       PHYSICAL_PROP_ARRAY physical,
                                       BC_ARRAY rest, BC_ARRAY force,
								       BC_ARRAY density, BC_ARRAY weight,
                                       DEFAULT_BC_ARRAY b_force,
                                       BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                                       DISP_ARRAY *disp, STRESS_ARRAY *stress,
                                       STRAIN_ARRAY *strain, double dens_pow,
                                       int sol_ss_flag, int lcm_num)
{
	FILE *fp, *fp1;
	CCS ccs;
	NONZERO_MATRIX_CCS tan_K;
	int ii0, ii1, ii2;
	int dim, total_dof;
	double beta1 = 1.0e-4;
	double beta2 = 1.0e-7;
	double *f_vect;
	double *df_vect;
	double *resid_v;
	double *react_v;
	double *robin_v;
	double *delta_d;
	double *d0_vect;
	double temp_f, temp_d;

	double omp_begin;
	double omp_end;

	/* MUMPS */
	DMUMPS_STRUC_C id;


	dim = analysis_dim(element);
	if(dim != 3) return(ARG_ERROR_RC);

	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );
	/* 荷重・変位曲線出力準備 */
	RC_TRY( rc_fopen("./Output_data/u-p.dat", "w", &fp1) );

	/* 全自由度 */
	RC_TRY( fem_total_dof(node, element, &total_dof) );

	RC_TRY( allocate1D(total_dof, &f_vect) );
	RC_TRY( allocate1D(total_dof, &df_vect) );
	RC_TRY( allocate1D(total_dof, &resid_v) );
	RC_TRY( allocate1D(total_dof, &react_v) );
	RC_TRY( allocate1D(total_dof, &robin_v) );
	RC_TRY( allocate1D(total_dof, &delta_d) );
	RC_TRY( allocate1D(total_dof, &d0_vect) );
	RC_TRY( vect2disp(dim, node, delta_d, disp) );

	/* 節点荷重，重力入力 */
	RC_TRY( add_nodal_force(dim, node, force, f_vect) );
#if 0
	RC_TRY( add_body_force(dim, node, element, material, physical, b_force,
                           f_vect) );
	/* 熱ひずみ(初期ひずみ)による荷重入力 */
	RC_TRY( add_temp_force(dim, node, element, material, physical, temp,
                           def_temp, f_vect) );
#endif

	for(ii1=0; ii1<total_dof; ii1++){
		df_vect[ii1] = f_vect[ii1] / (double)lcm_num;
		f_vect[ii1] = 0.0;
	}
	for(ii0=0; ii0<lcm_num; ii0++){
		for(ii1=0; ii1<total_dof; ii1++){
			f_vect[ii1] += df_vect[ii1];
			resid_v[ii1] = f_vect[ii1] - react_v[ii1] - robin_v[ii1];
		}
		for(ii1=0; ii1<NR_MAX; ii1++){
			double du_norm, du_r, u_norm, u_f;
			double u0_norm;

			RC_TRY( log_printf(5, " +-------------------------+\n") );
			RC_TRY( log_printf(5, " |--------   %3d   --------|\n", ii0+1) );
			RC_TRY( log_printf(5, " +-------------------------+\n") );
			RC_TRY( log_printf(5, " | N-R iteration NUM : %3d |\n", ii1+1) );
			RC_TRY( log_printf(5, " +-------------------------+\n") );
			if( (ii1 == 0) && (ii0 == 0) ){
				RC_TRY( make_nonzero1a_ccs_sym_topology_robin(node, element,
				                                             surf, material,
				                                             physical, density,
				                                             weight, dens_pow,
				                                             &tan_K) );
			}else{
				RC_TRY( make_NL_nonzero1a_ccs_sym_topology_robin(node, element,
				                                            surf, material,
				                                            physical, *disp,
				                                            density, weight,
				                                            dens_pow, &tan_K));
			}
			RC_TRY( copy_nonzero_matrix1a_ccs2ccs_format(&tan_K, &ccs) );
			RC_TRY( free_nonzero_matrix1a_ccs(&tan_K) );
			RC_TRY( modify_ccs_format(dim, node, rest, &ccs, resid_v) );

			omp_begin = omp_get_wtime();
			/* MUMPS */
			/* Inititalization */
			id.job = JOB_INIT;
			id.par = 1;
			/* Symmetric or Non-symmetric parameter */
			id.sym = 1;
			//  id.comm_fortran=USE_COMM_WORLD;
			dmumps_c(&id);
			/* Set Matrix data */
			id.n = total_dof;
			id.nz = ccs.array_size;
			id.irn = ccs.row_index;
			id.jcn = ccs.col_index;
			for(ii2=0; ii2<id.nz; ii2++){
				id.irn[ii2]++;
				id.jcn[ii2]++;
			}
			id.a = ccs.array;
			id.rhs = resid_v;
			/* No outputs */
			id.icntl[0]=-1; id.icntl[1]=-1; id.icntl[2]=-1; id.icntl[3]=0;
			id.icntl[21]=0;
			/* Symbolic + Factorize + Solve */
			id.job = 6;
			dmumps_c(&id);
			/* End Job */
			id.job = JOB_END;
			dmumps_c(&id);
			omp_end = omp_get_wtime();
			printf("MUMPS OMP Solve  OMPTIME %5.3lf sec\n",
			       (omp_end-omp_begin));

			for(ii2=0; ii2<total_dof; ii2++){
				delta_d[ii2] = resid_v[ii2];
			}

			if(ii1 == 0){
				for(ii2=0; ii2<total_dof; ii2++){
					d0_vect[ii2] = delta_d[ii2];
				}
			}

			RC_TRY( free_nonzero_ccs(&ccs) );
			RC_TRY( add_vect2disp(dim, node, delta_d, *disp, 1.0) );

			/* 反力計算 */
			RC_TRY( make_NL_react_force_topology(element, node, material,
			                                     physical, density, dens_pow,
			                                     *disp, react_v) );
			RC_TRY( make_robin_react_force(surf, node, material, physical,
                                           weight, *disp, robin_v) );
			/* 残差荷重ベクトル */
			for(ii2=0; ii2<total_dof; ii2++){
				resid_v[ii2] = f_vect[ii2] - react_v[ii2] - robin_v[ii2];
			}
			RC_TRY( set_rest_zero(node, rest, resid_v) );

			/* 変位，エネルギー収束判定 */
			du_norm = sqrt( inner_product(total_dof, delta_d, delta_d) );
			du_r = 0.0;
			for(ii2=0; ii2<total_dof; ii2++){
				du_r += fabs(delta_d[ii2] * resid_v[ii2]);
			}
			RC_TRY( copy_disp2vect(dim, node, *disp, delta_d, 1.0) );
			u_norm = sqrt( inner_product(total_dof, delta_d, delta_d) );
			u0_norm = sqrt( inner_product(total_dof, d0_vect, d0_vect) );
			u_f = 0.0;
			for(ii2=0; ii2<total_dof; ii2++){
				u_f += fabs(delta_d[ii2] * f_vect[ii2]);
			}
			RC_TRY( log_printf(5, " Disp   error :%15.7e\n", du_norm/u0_norm) );
			RC_TRY( log_printf(5, " Energy error :%15.7e\n", du_r/u_f) );
			/* 収束判定 */
			if(du_norm < beta1*u0_norm + ABS_TOL
			&& du_r < beta2*u_f + ABS_TOL){
#if 0
				/* 条件数計算用の変位情報出力 */
				RC_TRY( rc_fopen("./Output_data/interim_disp_binary.data", 
							     "w", &fp) );
				RC_TRY( write_disp_array(fp, *disp) );
				RC_TRY( rc_fclose(fp) );
#endif
				break;
			}
			RC_TRY( log_printf(5, "\n") );
		}
		if(ii1 == NR_MAX){
			log_printf(3, "N-R iteration NOT converged!!\n");
		}
		RC_TRY( log_printf(5, "---\n") );
		RC_TRY( log_printf(5, "\n") );
		/* 荷重・変位曲線出力 */
		temp_d = disp->array[PRI_NODE].v.x;
		temp_f = f_vect[dim*PRI_NODE];
		fprintf(fp1, "%lf %lf\n", temp_f, temp_d);

	}
	/* メモリ開放 */
	RC_TRY( mm_free(f_vect) );
	RC_TRY( mm_free(df_vect) );
	RC_TRY( mm_free(resid_v) );
	RC_TRY( mm_free(react_v) );
	RC_TRY( mm_free(robin_v) );
	RC_TRY( mm_free(delta_d) );
	RC_TRY( mm_free(d0_vect) );

	/* 第2Piola-Kirchhoff応力，Green歪み(応力解析は位相最適化用にしていない) */
	RC_TRY( PK2_stress_Green_strain(sol_ss_flag, node, element, *disp,
                                    material, physical, stress, strain) );
	RC_TRY( rc_fclose(fp1) );
	return(NORMAL_RC);
}
#endif

/* robin型ver. */
RC
fem_NL_solve_mumps_topology_robin (NODE_ARRAY node, ELEMENT_ARRAY element,
                                   ELEMENT_ARRAY surf,
                                   MATERIAL_PROP_ARRAY material,
                                   PHYSICAL_PROP_ARRAY physical,
                                   BC_ARRAY rest, BC_ARRAY force,
								   BC_ARRAY density, BC_ARRAY weight,
                                   DEFAULT_BC_ARRAY b_force,
                                   BC_ARRAY temp, DEFAULT_BC_ARRAY def_temp,
                                   DISP_ARRAY *disp, STRESS_ARRAY *stress,
                                   STRAIN_ARRAY *strain, double dens_pow,
                                   int sol_ss_flag)
{
	CCS ccs;
	NONZERO_MATRIX_CCS tan_K;
	int ii1, ii2;
	int dim, total_dof;
	double beta1 = 1.0e-3;
	double beta2 = 1.0e-6;
	double *f_vect;
	double *resid_v;
	double *react_v;
	double *robin_v;
	double *delta_d;

	double omp_begin;
	double omp_end;

	/* MUMPS */
	DMUMPS_STRUC_C id;


	dim = analysis_dim(element);
	if(dim != 3) return(ARG_ERROR_RC);

	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

	/* 全自由度 */
	RC_TRY( fem_total_dof(node, element, &total_dof) );

	RC_TRY( allocate1D(total_dof, &f_vect) );
	RC_TRY( allocate1D(total_dof, &resid_v) );
	RC_TRY( allocate1D(total_dof, &react_v) );
	RC_TRY( allocate1D(total_dof, &robin_v) );
	RC_TRY( allocate1D(total_dof, &delta_d) );
	RC_TRY( vect2disp(dim, node, delta_d, disp) );

	/* 節点荷重，重力入力 */
	RC_TRY( add_nodal_force(dim, node, force, f_vect) );
	RC_TRY( add_body_force(dim, node, element, material, physical, b_force,
                           f_vect) );
	/* 熱ひずみ(初期ひずみ)による荷重入力 */
	RC_TRY( add_temp_force(dim, node, element, material, physical, temp,
                           def_temp, f_vect) );

	for(ii1=0; ii1<total_dof; ii1++){
		resid_v[ii1] = f_vect[ii1];
	}
	for(ii1=0; ii1<NR_MAX; ii1++){
		double du_norm, du_r, u_norm, u_f;

		RC_TRY( log_printf(5, " +-------------------------+\n") );
		RC_TRY( log_printf(5, " | N-R iteration NUM : %3d |\n", ii1+1) );
		RC_TRY( log_printf(5, " +-------------------------+\n") );
		if(ii1 == 0){
			RC_TRY( make_nonzero1a_ccs_sym_topology_robin(node, element,
			                                              surf, material,
                                                          physical, density,
		                                                  weight, dens_pow,
			                                              &tan_K) );
		}else{
			RC_TRY( make_NL_nonzero1a_ccs_sym_topology_robin(node, element,
			                                                 surf, material,
                                                             physical, *disp,
			                                                 density, weight,
			                                                 dens_pow, &tan_K));
		}
		RC_TRY( copy_nonzero_matrix1a_ccs2ccs_format(&tan_K, &ccs) );
		RC_TRY( free_nonzero_matrix1a_ccs(&tan_K) );
		RC_TRY( modify_ccs_format(dim, node, rest, &ccs, resid_v) );
		
		omp_begin = omp_get_wtime();
		/* MUMPS */
		/* Inititalization */
		id.job = JOB_INIT;
		id.par = 1;
		/* Symmetric or Non-symmetric parameter */
		id.sym = 1;
//  id.comm_fortran=USE_COMM_WORLD;
		dmumps_c(&id);
		/* Set Matrix data */
		id.n = total_dof;
		id.nz = ccs.array_size;
		id.irn = ccs.row_index;
		id.jcn = ccs.col_index;
		for(ii2=0; ii2<id.nz; ii2++){
			id.irn[ii2]++;
			id.jcn[ii2]++;
		}
		id.a = ccs.array;
		id.rhs = resid_v;
		/* No outputs */
		id.icntl[0]=-1; id.icntl[1]=-1; id.icntl[2]=-1; id.icntl[3]=0;
		id.icntl[21]=0;
		/* Symbolic + Factorize + Solve */
		id.job = 6;
		dmumps_c(&id);
		/* End Job */
		id.job = JOB_END;
		dmumps_c(&id);
		omp_end = omp_get_wtime();
		printf("MUMPS OMP Solve  OMPTIME %5.3lf sec\n", (omp_end-omp_begin));

		for(ii2=0; ii2<total_dof; ii2++){
			delta_d[ii2] = resid_v[ii2];
		}
		
		RC_TRY( free_nonzero_ccs(&ccs) );
		RC_TRY( add_vect2disp(dim, node, delta_d, *disp, 1.0) );

		/* 反力計算 */
		RC_TRY( make_NL_react_force_topology(element, node, material, physical,
                                             density, dens_pow,
		                                     *disp, react_v) );
		RC_TRY( make_robin_react_force(surf, node, material, physical,
                                       weight, *disp, robin_v) );
		/* 残差荷重ベクトル */
		for(ii2=0; ii2<total_dof; ii2++){
			resid_v[ii2] = f_vect[ii2] - react_v[ii2] - robin_v[ii2];
		}
		RC_TRY( set_rest_zero(node, rest, resid_v) );

		/* 変位，エネルギー収束判定 */
		du_norm = sqrt( inner_product(total_dof, delta_d, delta_d) );
		du_r = 0.0;
		for(ii2=0; ii2<total_dof; ii2++){
			du_r += fabs(delta_d[ii2] * resid_v[ii2]);
		}
		RC_TRY( copy_disp2vect(dim, node, *disp, delta_d, 1.0) );
		u_norm = sqrt( inner_product(total_dof, delta_d, delta_d) );
		u_f = 0.0;
		for(ii2=0; ii2<total_dof; ii2++){
			u_f += fabs(delta_d[ii2] * f_vect[ii2]);
		}
		RC_TRY( log_printf(5, "    Disp   error :%15.7e\n", du_norm/u_norm) );
		RC_TRY( log_printf(5, "    Energy error :%15.7e\n", du_r/u_f) );
		/* 収束判定 */
		if(du_norm < beta1*u_norm + ABS_TOL
        && du_r < beta2*u_f + ABS_TOL){
			break;
		}
		RC_TRY( log_printf(5, "\n") );
	}
	if(ii1 == NR_MAX){
		log_printf(3, "N-R iteration NOT converged!!\n");
	}
	/* メモリ開放 */
	RC_TRY( mm_free(f_vect) );
	RC_TRY( mm_free(resid_v) );
	RC_TRY( mm_free(react_v) );
	RC_TRY( mm_free(robin_v) );
	RC_TRY( mm_free(delta_d) );

	/* 第2Piola-Kirchhoff応力，Green歪み(応力解析は位相最適化用にしていない) */
	RC_TRY( PK2_stress_Green_strain(sol_ss_flag, node, element, *disp,
                                    material, physical, stress, strain) );

	return(NORMAL_RC);
}


RC
make_robin_react_force (ELEMENT_ARRAY surf, NODE_ARRAY node,
                        MATERIAL_PROP_ARRAY material,
                        PHYSICAL_PROP_ARRAY physical,
                        BC_ARRAY weight, DISP_ARRAY disp,
                        double robin_force[])
{
	int ii1, ii2, ii3;
	int dim;
	int total_dof;
	double *d_vect;
	ELEM_MATRIX elem_matrix;
	WORK_SET ws;

	dim = analysis_dim(surf) + 1;
	if(dim != 3){
		fprintf(stderr, "surf_elem_dim() < ");
		return(ARG_ERROR_RC);
	}

	total_dof = dim * node.size;

	RC_TRY( disp2vect(dim, node, disp, &d_vect) );
	RC_TRY( allocate_work_set(&ws) );

	for(ii1=0; ii1<total_dof; ii1++) robin_force[ii1] = 0.0;

	for(ii1=0; ii1<surf.size; ii1++){
		if(surf.array[ii1].label < 0) continue;

		RC_TRY( make_NtN_elem_matrix(dim, surf.array[ii1], node, physical,
		                             weight, ws, &elem_matrix) );

		for(ii2=0; ii2<elem_matrix.size; ii2++){
			for(ii3=0; ii3<elem_matrix.size; ii3++){
				robin_force[elem_matrix.index[ii2]]
				+= elem_matrix.matrix[ii2][ii3]*d_vect[elem_matrix.index[ii3]];
			}
		}

		RC_TRY( free_elem_matrix(&elem_matrix) );
	}

	RC_TRY( free_work_set(&ws) );
	RC_TRY( free1D(total_dof, &d_vect) );

	return(NORMAL_RC);
}


RC
make_nonzero1a_ccs_sym_topology_robin (NODE_ARRAY node, ELEMENT_ARRAY element,
                                       ELEMENT_ARRAY surf,
                                       MATERIAL_PROP_ARRAY material,
                                       PHYSICAL_PROP_ARRAY physical,
								       BC_ARRAY density, BC_ARRAY weight,
                                       double dens_pow,
                                       NONZERO_MATRIX_CCS *matrix)
{
	int ii1, ii2, ii3;
	int dim;
	double *ptr;
	ELEM_MATRIX elem_matrix;
	WORK_SET ws;

	dim = analysis_dim(element);
    if(dim != 3) return(ARG_ERROR_RC);

	RC_TRY( allocate_work_set(&ws) );
	RC_TRY( allocate_nonzero_matrix1a_ccs(node.size*dim, matrix) );
	for(ii1=0; ii1<element.size; ii1++){

		if(element.array[ii1].label < 0) continue;

		RC_TRY( make_elem_matrix_topology(element.array[ii1], node,
                                          material, physical,
		                                  density, dens_pow,
		                                  &elem_matrix) );

		for(ii2=0; ii2<elem_matrix.size; ii2++){
			for(ii3=0; ii3<=ii2; ii3++){
				ptr = nonzero_matrix1a_ccs_ptr(matrix, elem_matrix.index[ii2],
                                               elem_matrix.index[ii3]);
				(*ptr) += elem_matrix.matrix[ii2][ii3];
			}
		}
		RC_TRY( free_elem_matrix(&elem_matrix) );
	}

	for(ii1=0; ii1<surf.size; ii1++){
		if(surf.array[ii1].label < 0) continue;

		RC_TRY( make_NtN_elem_matrix(dim, surf.array[ii1], node,
		                             physical, weight, ws, &elem_matrix) );

		for(ii2=0; ii2<elem_matrix.size; ii2++){
			for(ii3=0; ii3<=ii2; ii3++){
				ptr = nonzero_matrix1a_ccs_ptr(matrix, elem_matrix.index[ii2],
						elem_matrix.index[ii3]);
				(*ptr) += elem_matrix.matrix[ii2][ii3];
			}
		}

		RC_TRY( free_elem_matrix(&elem_matrix) );
	}
	RC_TRY( free_work_set(&ws) );

	return(NORMAL_RC);
}

RC
make_NL_nonzero1a_ccs_sym_topology_robin (NODE_ARRAY node,
                                          ELEMENT_ARRAY element,
                                          ELEMENT_ARRAY surf,
                                          MATERIAL_PROP_ARRAY material,
                                          PHYSICAL_PROP_ARRAY physical,
                                          DISP_ARRAY disp, BC_ARRAY density,
                                          BC_ARRAY weight, double dens_pow,
                                          NONZERO_MATRIX_CCS *K)
{
	int ii1, ii2, ii3;
	int dim;
	double *ptr;
	WORK_SET ws;
	ELEM_MATRIX elem_matrix;

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	dim = analysis_dim(element);
	if(dim != 3) return(ARG_ERROR_RC);

	RC_TRY( allocate_nonzero_matrix1a_ccs(node.size*dim, K) );
	RC_TRY( allocate_work_set(&ws) );

	for(ii1=0; ii1<element.size; ii1++){

		if(element.array[ii1].label < 0) continue;

		RC_TRY( make_tangent_elem_matrix_topology(element.array[ii1], node,
                                                  material, physical, ws, disp,
		                                          density, dens_pow,
                                                  &elem_matrix) );
		for(ii2=0; ii2<elem_matrix.size; ii2++){
			for(ii3=0; ii3<=ii2; ii3++){
				ptr = nonzero_matrix1a_ccs_ptr(K, elem_matrix.index[ii2],
                                               elem_matrix.index[ii3]);
				(*ptr) += elem_matrix.matrix[ii2][ii3];
			}
		}

		RC_TRY( free_elem_matrix(&elem_matrix) );
	}

	for(ii1=0; ii1<surf.size; ii1++){
		if(surf.array[ii1].label < 0) continue;

		RC_TRY( make_NtN_elem_matrix(dim, surf.array[ii1], node,
		                             physical, weight, ws, &elem_matrix) );

		for(ii2=0; ii2<elem_matrix.size; ii2++){
			for(ii3=0; ii3<=ii2; ii3++){
				ptr = nonzero_matrix1a_ccs_ptr(K, elem_matrix.index[ii2],
						elem_matrix.index[ii3]);
				(*ptr) += elem_matrix.matrix[ii2][ii3];
			}
		}

		RC_TRY( free_elem_matrix(&elem_matrix) );
	}
	RC_TRY( free_work_set(&ws) );

	return(NORMAL_RC);
}

/* 密度を考慮した線形の要素剛性マトリックスの作成 */
RC
make_elem_matrix_topology (ELEMENT element, NODE_ARRAY node,
                           MATERIAL_PROP_ARRAY material,
                           PHYSICAL_PROP_ARRAY physical,
                           BC_ARRAY density, double dens_pow,
						   ELEM_MATRIX *elem_matrix)
{
	int ii1, ii2, ii3;
	double thickness;
	double local_density;
	double local_phi;
	double node_density_vect[MAX_NODE];
	double N[MAX_NODE];
	static double **D_matrix = NULL;
	static double **BtDB = NULL;
	static int init_flag = 0;
	static double **node_location = NULL;
	static G_PROP g_prop;

	if(init_flag == 0){
		RC_TRY( allocate2D(6, 6, &D_matrix) );
		RC_TRY( allocate2D(3*MAX_NODE, 3*MAX_NODE, &BtDB) );
		RC_TRY( allocate2D(MAX_NODE, 3, &node_location) );
		RC_TRY( allocate_g_prop(&g_prop) );
		init_flag = 1;
	}

	if((elem_matrix->dim = element_dim(element.type)) <= 0){
		return(ARG_ERROR_RC);
	}

	if(element.label < 0) return(ARG_ERROR_RC);
	elem_matrix->element = element.label;
	elem_matrix->size = (element.node_num) * (elem_matrix->dim);

	/* elem_matrix->array[][] の確保 */
	RC_TRY( allocate2D(elem_matrix->size,
	                   elem_matrix->size, &(elem_matrix->matrix)) );

	/* elem_matirx->index[] の確保と代入 */

	elem_matrix->index = mm_alloc(sizeof(int)*(elem_matrix->size));
	if(elem_matrix->index == NULL) return(ALLOC_ERROR_RC);

	for(ii1=0; ii1<element.node_num; ii1++){
		int index = node_renum_index(node, element.node[ii1]);
		if(index < 0) return(SEEK_ERROR_RC);

		for(ii2=0; ii2<(elem_matrix->dim); ii2++){
			elem_matrix->index[ii1*(elem_matrix->dim)+ii2]
			                   = index*(elem_matrix->dim) + ii2;
		}
		/* 節点上の密度のベクトルを作成 */
		node_density_vect[ii1] = node_density(density, element.node[ii1]);
	}

	/* D マトリックス */
	RC_TRY( get_D_matrix(element, material, physical, D_matrix) );
	RC_TRY( get_thickness(element, physical, &thickness) );

	/* B マトリックス */
	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( set_Gauss_const(element.type, &g_prop) );
	RC_TRY( set_g_prop(element, node_location, &g_prop) );
//	RC_TRY( RI_B_matrix(element, g_prop) );

	for(ii1=0; ii1<g_prop.num_g_point; ii1++){
		/* B^T * D * B => BtDB */
		mul_matrix_AtBA(g_prop.B_matrix_row, g_prop.B_matrix_col,
		                g_prop.B_matrix[ii1], D_matrix, BtDB);

		/* ガウスポイント上の密度を求める */
		RC_TRY( set_N_vector(element.type, g_prop.g_points[ii1], N) );
		local_density = 0.0;
		for(ii2=0; ii2<element.node_num; ii2++){
			local_density += node_density_vect[ii2] * N[ii2];
		}
		local_phi
		= DENS_UPPER_LIMIT * pow( (tanh(local_density)+1.0)*0.5 , dens_pow );

		/* BDB * detJ * weight * thickness */
		for(ii2=0; ii2<(elem_matrix->size); ii2++){
			for(ii3=0; ii3<(elem_matrix->size); ii3++){
				elem_matrix->matrix[ii2][ii3]
				     += BtDB[ii2][ii3] * local_phi * g_prop.det_J[ii1]
				      * g_prop.weight[ii1] * thickness;
			}
		}
	}

	return(NORMAL_RC);
}

/* 密度を考慮した非線形の要素剛性マトリックスの作成 */
RC
make_tangent_elem_matrix_topology (ELEMENT element, NODE_ARRAY node,
                                   MATERIAL_PROP_ARRAY material,
                                   PHYSICAL_PROP_ARRAY physical,
                                   WORK_SET ws, DISP_ARRAY disp,
                                   BC_ARRAY density, double dens_pow,
                                   ELEM_MATRIX *elem_matrix)
{
	int dim, ssnum;
	int ii1, ii2, ii3;
	int g_num_point;
	double thickness;
	double det_J;
	double local_density; /* 密度を制御する変数 */
	double local_phi; /* 密度(定義域0〜1)を表すファイ */
	double node_density_vect[MAX_NODE];
	double N[MAX_NODE];
	double PK2_stress_v[6];
	double Green_strain_v[6];
	double disp_v[3*MAX_NODE];
	double *g_weights = ws.vec_G[0];
	double **node_location = ws.mat_N_3[0];
	double **D_matrix = ws.mat_6_6[0];
	double **B_matrix = ws.mat_6_3N[0];
	double **A_matrix = ws.mat_9_9[0];
	double **G_matrix = ws.mat_9_3N[0];
	double **AG = ws.mat_6_3N[1];
	double **B_strain = ws.mat_6_3N[2];
	double **M_matrix = ws.mat_9_9[1];
	double **BtDB = ws.mat_3N_3N[0];
	double **GtMG = ws.mat_3N_3N[1];
	double **g_points = ws.mat_G_4[0];

	if(element.label < 0) return(ARG_ERROR_RC);

	if( (elem_matrix->dim = element_dim(element.type)) <= 0 ||
	             (ssnum = element_ssnum(element.type)) <= 0 ){
		return(ARG_ERROR_RC);
	}

	dim = elem_matrix->dim;
	elem_matrix->element = element.label;
	elem_matrix->size = (element.node_num) * dim;

	/* elem_matrix->array[][]の確保 */
	RC_TRY( allocate2D(elem_matrix->size, elem_matrix->size,
	                   &(elem_matrix->matrix)) );
	RC_TRY( allocate1I(elem_matrix->size, &(elem_matrix->index)) );

	for(ii1=0; ii1<element.node_num; ii1++){
		int index = node_renum_index(node, element.node[ii1]);
		if(index < 0) return(SEEK_ERROR_RC);
		for(ii2=0; ii2<dim; ii2++){
			elem_matrix->index[ii1*dim+ii2] = index*dim + ii2;
		}
		/* 節点上の密度のベクトルを作成 */
		node_density_vect[ii1] = node_density(density, element.node[ii1]);
	}

	RC_TRY( Gauss_const_default(element.type, g_points, g_weights,
	                            &g_num_point) );
	RC_TRY( fill_disp_vector(element, disp, disp_v) );
	RC_TRY( node_location_matrix(&element, node, node_location) );
	RC_TRY( get_D_matrix(element, material, physical, D_matrix) );
	RC_TRY( get_thickness(element, physical, &thickness) );

	for(ii1=0; ii1<g_num_point; ii1++){
		STRESS_STRAIN Green_strain;

		RC_TRY( make_B_matrix(element, g_points[ii1], node_location, B_matrix,
		                     &det_J) );
		RC_TRY( make_AG_matrix(element, g_points[ii1], ws, node_location, disp,
		                       A_matrix, G_matrix) );
		mul_matrix(ssnum,dim*dim,dim*element.node_num,A_matrix,G_matrix,AG);
		RC_TRY( add_matrix(ssnum, dim*element.node_num, 1.0, B_matrix, 0.5, AG,
		                   B_strain) );
		mul_matrix_vect(ssnum, dim*element.node_num, B_strain, disp_v,
		                Green_strain_v);
		RC_TRY( vect2stress_strain(ssnum, Green_strain_v, &Green_strain) );
		RC_TRY( add_matrix(ssnum, dim*element.node_num, 1.0, B_matrix, 1.0, AG,
		                   B_matrix) );
		
		/* ガウスポイント上の密度を求める */
		RC_TRY( set_N_vector(element.type, g_points[ii1], N) );
		local_density = 0.0;
		for(ii2=0; ii2<element.node_num; ii2++){
			local_density += node_density_vect[ii2] * N[ii2];
		}
		local_phi
		= DENS_UPPER_LIMIT * pow( (tanh(local_density)+1.0)*0.5 , dens_pow );

		if(dim == 2){
			for(ii2=0; ii2<ssnum; ii2++){
				PK2_stress_v[ii2] = D_matrix[ii2][0]*Green_strain.v.x
				                  + D_matrix[ii2][1]*Green_strain.v.y
				                  + D_matrix[ii2][2]*Green_strain.v.xy;
			}
		}else{ /* dim == 3 */
			for(ii2=0; ii2<ssnum; ii2++){
				PK2_stress_v[ii2] = D_matrix[ii2][0]*Green_strain.v.x
				                  + D_matrix[ii2][1]*Green_strain.v.y
				                  + D_matrix[ii2][2]*Green_strain.v.z
				                  + D_matrix[ii2][3]*Green_strain.v.yz
				                  + D_matrix[ii2][4]*Green_strain.v.zx
				                  + D_matrix[ii2][5]*Green_strain.v.xy;
			}
		}

		RC_TRY( set_M_matrix(dim, PK2_stress_v, M_matrix) );

		mul_matrix_AtBA(ssnum, dim*element.node_num, B_matrix, D_matrix, BtDB);
		mul_matrix_AtBA(dim*dim,dim*element.node_num,G_matrix,M_matrix,GtMG);

		for(ii2=0; ii2<(elem_matrix->size); ii2++){
			for(ii3=0; ii3<(elem_matrix->size); ii3++){
				elem_matrix->matrix[ii2][ii3] += local_phi
				                               * (BtDB[ii2][ii3]+GtMG[ii2][ii3])
				                               * det_J*g_weights[ii1]*thickness;
			}
		}
	}

	return(NORMAL_RC);
}

/* 密度を考慮した内力マトリックスの作成 */
RC
make_NL_react_force_topology (ELEMENT_ARRAY element, NODE_ARRAY node,
                              MATERIAL_PROP_ARRAY material,
                              PHYSICAL_PROP_ARRAY physical,
                              BC_ARRAY density, double dens_pow,
                              DISP_ARRAY disp, double f_vect[])
{
	int dim, ssnum;
	int ii1, ii2, ii3;
	int total_dof;
	int g_num_point;
	double thickness;
	double det_J;
	double local_density; /* 密度を制御する変数 */
	double local_phi; /* 密度(定義域0〜1)を表すファイ */
	double node_density_vect[MAX_NODE];
	double N[MAX_NODE];
	double react_v[3*MAX_NODE];
	double PK2_stress_v[6];
	double Green_strain_v[6];
	double disp_v[3*MAX_NODE];
	double *g_weights = NULL;
	double **D_matrix = NULL;
	double **B_matrix = NULL;
	double **A_matrix = NULL;
	double **G_matrix = NULL;
	double **AG = NULL;
	double **B_strain = NULL;
	double **node_location = NULL;
	double **g_points = NULL;
	WORK_SET ws;

	if( (dim = analysis_dim(element)) <= 1 || dim > 3){
		return(ARG_ERROR_RC);
	}

	RC_TRY( allocate_work_set(&ws) );
	node_location = ws.mat_N_3[0];
	g_weights = ws.vec_G[0];
	D_matrix = ws.mat_6_6[0];
	B_matrix = ws.mat_6_3N[0];
	A_matrix = ws.mat_9_9[0];
	G_matrix = ws.mat_9_3N[0];
	AG = ws.mat_6_3N[1];
	B_strain = ws.mat_6_3N[2];
	g_points = ws.mat_G_4[0];

	RC_TRY( fem_total_dof(node, element, &total_dof) );

	for(ii1=0; ii1<total_dof; ii1++) f_vect[ii1] = 0.0;
	
	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		if( (ssnum = element_ssnum(element.array[ii1].type)) <= 0 ){
			return(ARG_ERROR_RC);
		}
		
		for(ii2=0; ii2<element.array[ii1].node_num; ii2++){
			node_density_vect[ii2] = node_density(density,
			                                      element.array[ii1].node[ii2]);
		}

		RC_TRY( Gauss_const_default(element.array[ii1].type, g_points,
					g_weights, &g_num_point) );
		RC_TRY( fill_disp_vector(element.array[ii1], disp, disp_v) );
		RC_TRY( node_location_matrix(&(element.array[ii1]), node,
					node_location) );
		RC_TRY( get_D_matrix(element.array[ii1],material,physical,D_matrix) );
		RC_TRY( get_thickness(element.array[ii1], physical, &thickness) );

		for(ii2=0; ii2<g_num_point; ii2++){
			STRESS_STRAIN Green_strain;

			RC_TRY( make_B_matrix(element.array[ii1], g_points[ii2],
						node_location, B_matrix, &det_J) );
			RC_TRY( make_AG_matrix(element.array[ii1], g_points[ii2], ws,
						node_location, disp, A_matrix, G_matrix) );
			mul_matrix(ssnum, dim*dim, dim*element.array[ii1].node_num,
					A_matrix, G_matrix, AG);
			RC_TRY( add_matrix(ssnum, dim*element.array[ii1].node_num, 1.0,
						B_matrix, 0.5, AG, B_strain) );
			mul_matrix_vect(ssnum, dim*element.array[ii1].node_num, B_strain,
					disp_v, Green_strain_v);
			RC_TRY( vect2stress_strain(ssnum, Green_strain_v, &Green_strain) );
			RC_TRY( add_matrix(ssnum, dim*element.array[ii1].node_num, 1.0,
						B_matrix, 1.0, AG, B_matrix) );
			
			/* ガウスポイント上の密度を求める */
			RC_TRY( set_N_vector(element.array[ii1].type, g_points[ii2], N) );
			local_density = 0.0;
			for(ii3=0; ii3<element.array[ii1].node_num; ii3++){
				local_density += node_density_vect[ii3] * N[ii3];
			}
			local_phi
			= DENS_UPPER_LIMIT * pow((tanh(local_density)+1.0)*0.5 ,dens_pow);

			if(dim==2){
				for(ii3=0; ii3<ssnum; ii3++){
					PK2_stress_v[ii3] = local_phi
					                  *(D_matrix[ii3][0]*Green_strain.v.x
					                  + D_matrix[ii3][1]*Green_strain.v.y
					                  + D_matrix[ii3][2]*Green_strain.v.xy);
				}
			}else{  /* dim == 3 */
				for(ii3=0; ii3<ssnum; ii3++){
					PK2_stress_v[ii3] = local_phi
					                  *(D_matrix[ii3][0]*Green_strain.v.x
					                  + D_matrix[ii3][1]*Green_strain.v.y
					                  + D_matrix[ii3][2]*Green_strain.v.z
					                  + D_matrix[ii3][3]*Green_strain.v.yz
					                  + D_matrix[ii3][4]*Green_strain.v.zx
					                  + D_matrix[ii3][5]*Green_strain.v.xy);
				}
			}

			mul_trans_matrix_vect(ssnum, dim*element.array[ii1].node_num,
					B_matrix, PK2_stress_v, react_v);

			for(ii3=0; ii3<element.array[ii1].node_num; ii3++){
				int index;

				index = node_renum_index(node, element.array[ii1].node[ii3]);
				if(index < 0) return(SEEK_ERROR_RC);

				f_vect[index*dim] += g_weights[ii2] * det_J * thickness
					* react_v[ii3*dim];
				f_vect[index*dim+1] += g_weights[ii2] * det_J * thickness
					* react_v[ii3*dim+1];
				if(dim == 3){
					f_vect[index*dim+2] += g_weights[ii2] * det_J * thickness
					* react_v[ii3*dim+2];
				}
			}
		}
	}

	RC_TRY( free_work_set(&ws) );

	return(NORMAL_RC);
}


/* 節点上の密度 */
double
node_density (BC_ARRAY density, int node_label)
{
	int index;
	int ii1, ii2;

	/* 密度が変数でない節点では密度をとりあえず0.0とする */
	index = search_bc_node(density, node_label);
	if(index < 0) return(0.0);

	for(ii1=0; ii1<BC_SCALAR; ii1++){
		if(density.array[index].s_type[ii1] == BC_DENSITY){
			return(density.array[index].s[ii1]);
		}
	}

	return(0.0);
}


/* 拘束点にゼロ代入 */
static RC
set_rest_zero (NODE_ARRAY node, BC_ARRAY rest, double resid_v[])
{
	int ii1;
	int index;

	for(ii1=0; ii1<rest.size; ii1++){

		if(rest.array[ii1].node < 0) continue;

		index = node_renum_index(node, rest.array[ii1].node);
		if(index < 0) return(SEEK_ERROR_RC);

		if(rest.array[ii1].v_type.x == BC_FIX) resid_v[3*index] = 0.0;
		if(rest.array[ii1].v_type.y == BC_FIX) resid_v[3*index+1] = 0.0;
		if(rest.array[ii1].v_type.z == BC_FIX) resid_v[3*index+2] = 0.0;
	}

	return(NORMAL_RC);
}

static RC
allocate_g_prop (G_PROP *g_prop)
{
	int ii1;

	RC_TRY( allocate2D(MAX_GAUSS_POINT+1, 4, &(g_prop->g_points)) );

	for(ii1=0; ii1<MAX_GAUSS_POINT+1; ii1++){
		RC_TRY( allocate2D(6, 3*MAX_NODE, &(g_prop->B_matrix[ii1])) );
	}

	return(NORMAL_RC);
}

static RC
set_Gauss_const (ELEM_TYPE type, G_PROP *g_prop)
{
	switch(type){
		case ELEM_TRI1:
		case ELEM_TRI2:
			g_prop->num_g_point = 3;
			RC_TRY( Gauss_const_tri(g_prop->num_g_point,
						g_prop->g_points, g_prop->weight) );
			RC_TRY( Gauss_const_tri(1, &(g_prop->g_points[MAX_GAUSS_POINT]),
						&(g_prop->weight[MAX_GAUSS_POINT])) );
			break;
		case ELEM_QUAD1:
			g_prop->num_g_point = 4;
			RC_TRY( Gauss_const_quad(g_prop->num_g_point,
						g_prop->g_points, g_prop->weight) );
			RC_TRY( Gauss_const_quad(1, &(g_prop->g_points[MAX_GAUSS_POINT]),
						&(g_prop->weight[MAX_GAUSS_POINT])) );
			break;
		case ELEM_QUAD2:
			g_prop->num_g_point = 9;
			RC_TRY( Gauss_const_quad(g_prop->num_g_point,
						g_prop->g_points, g_prop->weight) );
			RC_TRY( Gauss_const_quad(1, &(g_prop->g_points[MAX_GAUSS_POINT]),
						&(g_prop->weight[MAX_GAUSS_POINT])) );
			break;
		case ELEM_TETRA1:
			g_prop->num_g_point = 4;
			RC_TRY( Gauss_const_tetra(g_prop->num_g_point,
						g_prop->g_points, g_prop->weight) );
			RC_TRY( Gauss_const_tetra(1, &(g_prop->g_points[MAX_GAUSS_POINT]),
						&(g_prop->weight[MAX_GAUSS_POINT])) );
			break;
		case ELEM_TETRA2:
			g_prop->num_g_point = 5;
			RC_TRY( Gauss_const_tetra(g_prop->num_g_point,
						g_prop->g_points, g_prop->weight) );
			RC_TRY( Gauss_const_tetra(1, &(g_prop->g_points[MAX_GAUSS_POINT]),
						&(g_prop->weight[MAX_GAUSS_POINT])) );
			break;
		case ELEM_PENTA1:
			g_prop->num_g_point = 6;
			RC_TRY( Gauss_const_penta(g_prop->num_g_point,
						g_prop->g_points, g_prop->weight) );
			RC_TRY( Gauss_const_penta(1, &(g_prop->g_points[MAX_GAUSS_POINT]),
						&(g_prop->weight[MAX_GAUSS_POINT])) );
			break;
		case ELEM_HEXA1:
			g_prop->num_g_point = 8;
			RC_TRY( Gauss_const_hexa(g_prop->num_g_point,
						g_prop->g_points, g_prop->weight) );
			RC_TRY( Gauss_const_hexa(1, &(g_prop->g_points[MAX_GAUSS_POINT]),
						&(g_prop->weight[MAX_GAUSS_POINT])) );
			break;
		case ELEM_HEXA2:
			g_prop->num_g_point = 14;
			RC_TRY( Gauss_const_hexa(g_prop->num_g_point,
						g_prop->g_points, g_prop->weight) );
			RC_TRY( Gauss_const_hexa(1, &(g_prop->g_points[MAX_GAUSS_POINT]),
						&(g_prop->weight[MAX_GAUSS_POINT])) );
			break;
		default:
			return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}

static RC
set_g_prop (ELEMENT element, double **node_location, G_PROP *g_prop)
{
	int ii1;
	int dim;

	if( (g_prop->B_matrix_row = element_ssnum(element.type)) <= 0){
		return(ARG_ERROR_RC);
	}

	if((dim = element_dim(element.type)) <= 0) return(ARG_ERROR_RC);

	g_prop->B_matrix_col = dim * element.node_num;

	/* 各ガウスポイント */
	for(ii1=0; ii1<(g_prop->num_g_point); ii1++){
		RC_TRY( make_B_matrix(element, g_prop->g_points[ii1], node_location,
					g_prop->B_matrix[ii1], &(g_prop->det_J[ii1])) );
	}

	/* 要素中心 */
	RC_TRY( make_B_matrix(element, g_prop->g_points[MAX_GAUSS_POINT],
				node_location, g_prop->B_matrix[MAX_GAUSS_POINT],
				&(g_prop->det_J[MAX_GAUSS_POINT])) );

	return(NORMAL_RC);
}

