/*********************************************************************
 * mtrix_db_amls.c
 *
 * Copyright (C) 2005 AzLib Developers Group
 *
 * Written by
 *  <Taiki AOYAMA> <Kenzen TAKEUCHI> 
 *
 * Thanks to following contributors.
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: matrix_db_amls.c 663 2005-12-12 09:51:36Z sasaoka $ */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rc.h"
#include "macros.h"
#include "memory_manager.h"
#include "nonzero_cg.h"
#include "matrix_db.h"
#include "sky_cholesky.h"
#include "fem_struct.h"
#include "fem_solver.h"
#include "amls.h"


static RC cal_eval_evect(int size, double **subK, double **subM,
                         double *e_vals, double **e_vects);
static RC cholesky_decomp(int n, double **matrix, double **L_matrix);
static void eigen_sort(int n, double *e_vals, double **e_vects);
static RC matrix_dbi2skyline3(int row_size, int row0_index, int block_size,
                              MATRIX_DB_INFO *mdbi, SKYLINE_MATRIX *skyline);
static RC s_cholesky_subst_sp3(long n, const long *index1, const long *index2,
                      const double *array, double *vector0,
                      double *vector1, double *vector2, int p_flag);
#if 0
static RC mul_AtBA_matrix_db(int row_size, int At_row_size,
                              MATRIX_DB_INFO *mdbi_At, MATRIX_DB_INFO *mdbi_B,
                              MATRIX_DB_INFO *mdbi_AtBA);
#endif
static RC mul_matrix_vect_db_sp33(int row_size, MATRIX_DB_INFO *A,
                                  const int min_col_A[], int vect_num,
                                  const int vect_size[], int **vect_row,
                                  double (**vect)[3][3], double (**ans)[3][3]);
static RC mul_AtBA_matrix_db(int row_size, int At_row_size,
                              MATRIX_DB_INFO *mdbi_At, MATRIX_DB_INFO *mdbi_B,
                              double **AtBA);
#if 1
static RC mul_ZtUt_matrix_db(int row_size, int block_num,
              const int block_size_array[], const int evects_size_array[],
              MATRIX_DB_INFO *mdbi_Zt, MATRIX_DB_INFO *mdbi_U,
              const int max_col_U[]);
#endif
static RC mul_ZtUt_matrix_db_old(int row_size, int block_num,
              int block_size_array[], const int evects_size_array[],
              MATRIX_DB_INFO *mdbi_Zt, MATRIX_DB_INFO *mdbi_U);
static RC chk_elim_block33(int row_size, int current_pos, int next_pos,
                           int elim_pos, int next_elim_pos,
                           MATRIX_DB_INFO *mdbi_K, int *chk_flag);
static RC cal_UtKU(MATRIX_DB_INFO *mdbi_K, int row_size, int c_index,
                   int current_pos, int elim_pos, int next_elim_pos,
                   int *col_array, double array[][3][3], int block_U_size,
                   int *block_U_index, double (*block_U_array)[3][3],
                   int min_col_index[]);
static RC make_U_matrix_db(MATRIX_DB_INFO *mdbi_U, int c_index, int elim_pos,
                           int next_elim_pos, int block_U_size,
                           int *block_U_index, double (*block_U_array)[3][3],
                           int max_col_U[]);
static RC vect2array33spm(int v_size, const double v0[], const double v1[],
                          const double v2[], int offset, int index[],
                          double array[][3][3], int *size);
static RC inner_product33sp(int size1, double array1[][3][3],
                const int index1[], int size2, double array2[][3][3],
                const int index2[], double v[3][3], int *set_flag);
static RC block_gauss_diag_K(int block_size, int g_block_num, 
                             int g_block_size_array[],
                             int row_size,
                             int elim_pos, int next_elim_pos, int current_pos,
                             int next_pos, int min_col_index[],
                             int col_array[], double array[][3][3],
                             MATRIX_DB_INFO *mdbi_K, MATRIX_DB_INFO *mdbi_G,
                             MATRIX_DB_INFO *mdbi_b, SKYLINE_MATRIX *skyline[]);
static RC block_gauss_subst3(int block_size, int g_block_num,
                             int g_block_size_array[], int col_array[],
                             double array[][3][3], MATRIX_DB_INFO *mdbi_G,
                             SKYLINE_MATRIX skyline[],
                             double v0[], double v1[], double v2[]);
static RC
cal_UtL (MATRIX_DB_INFO *mdbi_G, int block_size, int g_current_pos,
         int g_next_pos, int g_elim_pos, int g_next_elim_pos,
         int *col_array, double array[][3][3], int g_block_U_size[],
         int **g_block_U_index, double (**g_block_U_array)[3][3],
         int g_min_col_index[]);
static RC
cal_Utb (MATRIX_DB_INFO *mdbi_b, int row_size, int current_pos, int elim_pos,
         int g_current_pos, int g_next_pos, int g_elim_pos,
         int g_next_elim_pos, int *col_array, double array[][3][3],
         int g_block_U_size[], int **g_block_U_index,
         double (**g_block_U_array)[3][3],
         int g_min_col_index[], int min_col_index[]);


RC
result_chk_amls_bg (int block_num, int block_size_array[], int row_size,
                    int evects_size_array[], MATRIX_DB_INFO *mdbi_K,
                    MATRIX_DB_INFO *mdbi_M, MATRIX_DB_INFO *mdbi_U,
                    const int max_col_U[])
{
	int ii1, ii2, ii3, ii4, ii5;
	int vect_size;
	int real_vect_size;
	double *evals;
	double **evects;
	double (*array)[3][3];
	int *col_array;
	int col_size;
	int row_pos;
	int col_pos;
	int col_pos3;
	int diag_pos3;
	int next_diag_pos3;
	double **block_K;
	double **block_M;
	double **ZtUtKUZ;
	double **ZtUtMUZ;
	MATRIX_DB_INFO mdbi_Zt;


	RC_TRY( open_matrix_db(".", 128*1024*1024, &mdbi_Zt) );

	/* Make Z, ZtUtKUZ */
	vect_size = 0;
	for(ii1=0; ii1<block_num; ii1++){
		vect_size += evects_size_array[ii1] * 3;
	}
	RC_TRY( allocate2D(vect_size, vect_size, &ZtUtKUZ) );
	diag_pos3 = 0;
	next_diag_pos3 = 0;
	row_pos = 0;
	col_pos = 0;
	col_pos3 = 0;
	RC_TRY( allocate1D33(row_size, &array) );
	RC_TRY( allocate1I(row_size, &col_array) );
	for(ii1=0; ii1<block_num; ii1++){
		double v[3][3];
		int found_flag;
		int block_size = block_size_array[ii1] * 3;
		int block_size3 = block_size_array[ii1];

		next_diag_pos3 += block_size_array[ii1];

		RC_TRY( allocate1D(block_size, &evals) );
		RC_TRY( allocate2D(block_size, block_size, &evects) );
		RC_TRY( allocate2D(block_size, block_size, &block_K) );
		RC_TRY( allocate2D(block_size, block_size, &block_M) );

		for(ii2=0; ii2<block_size3; ii2++){
			RC_TRY( get_row_range_matrix_db(mdbi_K,ii2+diag_pos3, diag_pos3,
			                                diag_pos3+block_size3, &col_size,
			                                array, col_array) );
			for(ii3=0; ii3<(col_size-1); ii3++){
				for(ii4=0; ii4<3; ii4++){
					for(ii5=0; ii5<3; ii5++){
						block_K[ii2*3+ii4][(col_array[ii3]-diag_pos3)*3+ii5]
						= array[ii3][ii4][ii5];
						block_K[(col_array[ii3]-diag_pos3)*3+ii5][ii2*3+ii4]
						= array[ii3][ii4][ii5];
					}
				}
			}
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					block_K[ii2*3+ii3][(col_array[col_size-1]-diag_pos3)*3+ii4]
					= array[col_size-1][ii3][ii4];
				}
			}
		}
		/* block_M に元のM(集中質量)の対角ブロック代入 */
		for(ii2=diag_pos3; ii2<next_diag_pos3; ii2++){
			double v[3][3];

			RC_TRY( get_v_matrix_db(ii2, ii2, v, &found_flag, mdbi_M) );
			if(!found_flag) return(UNKNOWN_ERROR_RC);
			for(ii3=0; ii3<3; ii3++){
				block_M[(ii2-diag_pos3)*3+ii3]
				       [(ii2-diag_pos3)*3+ii3] = v[ii3][ii3];
			}
		}

		RC_TRY( cal_eval_evect(block_size, block_K, block_M, evals, evects) );

		/* ZtUtKUZ */
		for(ii2=0; ii2<evects_size_array[ii1]*3; ii2++){
			ZtUtKUZ[col_pos+ii2][col_pos+ii2] = evals[ii2];
		}
		/* mdbi_Zt */
		for(ii2=0; ii2<evects_size_array[ii1]; ii2++){
			for(ii3=0; ii3<block_size_array[ii1]; ii3++){
				/* Zt[ii2][ii3] -> v */
				for(ii4=0; ii4<3; ii4++){
					for(ii5=0; ii5<3; ii5++){
						v[ii4][ii5] = evects[ii3*3+ii5][ii2*3+ii4];
					}
				}
				RC_TRY(put_v_matrix_db(ii2+col_pos3,ii3+diag_pos3,v,&mdbi_Zt));
			}
		}
		RC_TRY( free1D(block_size, &evals) );
		RC_TRY( free2D(block_size, block_size, &evects) );
		RC_TRY( free2D(block_size, block_size, &block_K) );
		RC_TRY( free2D(block_size, block_size, &block_M) );

		diag_pos3 += block_size_array[ii1];
		row_pos += block_size_array[ii1] * 3;
		col_pos += evects_size_array[ii1] * 3;
		col_pos3 += evects_size_array[ii1];
	}
	RC_TRY( free1D33(row_size, &array) );
	RC_TRY( free1I(row_size, &col_array) );

	/* ZtUtMUZ */
	real_vect_size = 0;
	for(ii1=0; ii1<block_num; ii1++){
		real_vect_size += 3*(evects_size_array[ii1]);
	}

#if 0
	RC_TRY( mul_ZtUt_matrix_db(row_size, block_num, block_size_array,
	                           evects_size_array, &mdbi_Zt, mdbi_U,
	                           max_col_U) );
#endif
#if 1
	RC_TRY( mul_ZtUt_matrix_db_old(row_size, block_num, block_size_array,
	                           evects_size_array, &mdbi_Zt, mdbi_U) );
#endif

#if 0
	RC_TRY( open_matrix_db(".", 64*1024*1024, &mdbi_m) );
	RC_TRY( mul_AtBA_matrix_db(row_size, real_vect_size/3, &mdbi_Zt, mdbi_M,
	                           &mdbi_m) );
	RC_TRY( close_matrix_db(&mdbi_m) );
#endif

	RC_TRY( allocate2D(real_vect_size, real_vect_size, &ZtUtMUZ) );
	RC_TRY( mul_AtBA_matrix_db(row_size, real_vect_size/3, &mdbi_Zt, mdbi_M,
	                           ZtUtMUZ) );
#if 0
	RC_TRY( allocate1I(real_vect_size/3, &col_array) );
	RC_TRY( allocate1D33(real_vect_size/3, &array) );
	for(ii1=0; ii1<real_vect_size/3; ii1++){
		RC_TRY( get_row_matrix_db(&mdbi_m, ii1, &col_size, array, col_array) );
		for(ii2=0; ii2<col_size; ii2++){
			if(ii1 >= col_array[ii2]){
				for(ii3=0; ii3<3; ii3++){
					for(ii4=0; ii4<3; ii4++){
						ZtUtMUZ[ii1*3+ii3][col_array[ii2]*3+ii4]
						    = array[ii2][ii3][ii4];
						ZtUtMUZ[col_array[ii2]*3+ii4][ii1*3+ii3]
						    = array[ii2][ii3][ii4];
					}
				}
			}
		}
	}
	RC_TRY( free1I(real_vect_size/3, &col_array) );
	RC_TRY( free1D33(real_vect_size/3, &array) );
#endif

	RC_TRY( close_matrix_db(&mdbi_Zt) );

	RC_TRY( allocate1D(real_vect_size, &evals) );
	RC_TRY( allocate2D(real_vect_size, real_vect_size, &evects) );

	tlog_printf(5, "[%d]\n", __LINE__);
	RC_TRY( cal_eval_evect(real_vect_size, ZtUtKUZ, ZtUtMUZ, evals, evects) );
	tlog_printf(5, "[%d]\n", __LINE__);
	for(ii1=0; ii1<10; ii1++){
		RC_TRY( log_printf(5, "evals[%d] = %15.7e %15.7e\n", ii1+1, evals[ii1],
		                   sqrt(fabs(evals[ii1]))/(2.0*PI)) );
	}
	RC_TRY( free1D(real_vect_size, &evals) );
	RC_TRY( free2D(real_vect_size, real_vect_size, &evects) );
	RC_TRY( free2D(vect_size, vect_size, &ZtUtKUZ) );
	RC_TRY( free2D(real_vect_size, real_vect_size, &ZtUtMUZ) );

	return(NORMAL_RC);
}


#if 1
/* ZtUt -> Zt */
#define SUB_U_SIZE  (30)
static RC
mul_ZtUt_matrix_db (int row_size, int block_num, const int block_size_array[],
                    const int evects_size_array[], MATRIX_DB_INFO *mdbi_Zt,
                    MATRIX_DB_INFO *mdbi_U, const int max_col_U[])
{
	int ii1, ii2, ii3, ii4, ii5, ii6;
	int *block_index;
	int *evects_index;
	double (**sub_U)[3][3];
	int **sub_U_col;
	int *sub_U_size;
	int *sub_U_row;
	double (*Z_array)[3][3];
	int *Z_col_array;
	int Z_size;


	if(block_num <= 0) return(ARG_ERROR_RC);

	RC_TRY( allocate1I(block_num + 1, &block_index) );
	RC_TRY( allocate1I(block_num + 1, &evects_index) );

	block_index[0] = 0;
	evects_index[0] = 0;
	for(ii1=0; ii1<block_num; ii1++){
		block_index[ii1+1]  = block_index[ii1]  + block_size_array[ii1];
		evects_index[ii1+1] = evects_index[ii1] + evects_size_array[ii1];
	}

	for(ii1=block_num-2; ii1>=0; ii1--){
		tlog_printf(4, "[%d] ii1 = %d\n", __LINE__, ii1);
		for(ii2=block_num-1; ii2>ii1; ii2--){
			int sub_U_count = 0;

			RC_TRY( allocate2D33(SUB_U_SIZE, block_size_array[ii2], &sub_U) );
			RC_TRY( allocate2I(SUB_U_SIZE, block_size_array[ii2], &sub_U_col) );
			RC_TRY( allocate1I(SUB_U_SIZE, &sub_U_size) );
			RC_TRY( allocate1I(SUB_U_SIZE, &sub_U_row) );
			RC_TRY( allocate1D33(block_size_array[ii2], &Z_array) );
			RC_TRY( allocate1I(block_size_array[ii2], &Z_col_array) );

			for(ii3=ii2; ii3<block_num; ii3++){
				if(evects_index[ii3+1] - evects_index[ii3] <= 0) continue;

				for(ii4=block_index[ii1]; ii4<block_index[ii1+1];
				                          ii4+=sub_U_count){
					sub_U_count = 0;
					for(ii5=ii4; ii5<block_index[ii1+1]; ii5++){
						if(max_col_U[ii5] < block_index[ii2]) continue;

						RC_TRY( get_row_range_matrix_db(mdbi_U, ii5,
						                block_index[ii2], block_index[ii2+1],
						                &(sub_U_size[sub_U_count]),
						                sub_U[sub_U_count],
						                sub_U_col[sub_U_count]) );
						if(sub_U_size[sub_U_count] <= 0) continue;
						sub_U_row[sub_U_count] = ii5;
						sub_U_count++;
						if(sub_U_count >= SUB_U_SIZE) break;
					}
					if(sub_U_count <= 0) break;

					for(ii5=evects_index[ii3]; ii5<evects_index[ii3+1]; ii5++){
						RC_TRY( get_row_range_matrix_db(mdbi_Zt, ii5,
						                block_index[ii2], block_index[ii2+1],
						                &Z_size, Z_array, Z_col_array) );
						for(ii6=0; ii6<sub_U_count; ii6++){
							double ZU[3][3];
							int set_flag;

							RC_TRY( inner_product33sp(Z_size, Z_array,
							        Z_col_array, sub_U_size[ii6], sub_U[ii6],
							        sub_U_col[ii6], ZU, &set_flag) );
							if(set_flag){
								RC_TRY( add_v_matrix_db(ii5, sub_U_row[ii6],
								                        ZU, mdbi_Zt) );
							}
						}
					}
				}
			}

			RC_TRY( free2D33(SUB_U_SIZE, block_size_array[ii2], &sub_U) );
			RC_TRY( free2I(SUB_U_SIZE, block_size_array[ii2], &sub_U_col) );
			RC_TRY( free1I(SUB_U_SIZE, &sub_U_size) );
			RC_TRY( free1I(SUB_U_SIZE, &sub_U_row) );
			RC_TRY( free1D33(block_size_array[ii2], &Z_array) );
			RC_TRY( free1I(block_size_array[ii2], &Z_col_array) );
		}


	}
	RC_TRY( free1I(block_num + 1, &block_index) );
	RC_TRY( free1I(block_num + 1, &evects_index) );

	return(NORMAL_RC);
}
#endif


static RC
mul_ZtUt_matrix_db_old (int row_size, int block_num, int block_size_array[],
                        const int evects_size_array[], MATRIX_DB_INFO *mdbi_Zt,
                        MATRIX_DB_INFO *mdbi_U)
{
	int ii1, ii2, ii3, ii4, ii5, ii6, ii7, ii8;
	int elim_pos;
	int next_elim_pos;
	int Zt_col_size;
	int U_col_size;
	int *Zt_col_array;
	int *U_col_array;
	int Zt_row_size;
	double (*temp_array)[3][3];
	double (*Zt_array)[3][3];
	double (*Zt_vect)[3][3];
	double (*U_array)[3][3];


	RC_TRY( allocate1D33(row_size, &Zt_array) );
	RC_TRY( allocate1I(row_size, &Zt_col_array) );
	RC_TRY( allocate1D33(row_size, &U_array) );
	RC_TRY( allocate1I(row_size, &U_col_array) );
	RC_TRY( allocate1D33(row_size, &Zt_vect) );

	Zt_row_size = 0;
	for(ii1=0; ii1<block_num; ii1++) Zt_row_size += evects_size_array[ii1];

	next_elim_pos = row_size - block_size_array[block_num-1];
	elim_pos = next_elim_pos;

	for(ii1=block_num-2; ii1>=0; ii1--){
		int Zt_current_pos;
		int current_pos;
		int next_pos;

		elim_pos -= block_size_array[ii1];
		Zt_current_pos = Zt_row_size;
		current_pos = row_size;
		next_pos = row_size;
		RC_TRY( allocate1D33(next_elim_pos-elim_pos, &temp_array) );
		for(ii2=block_num-1; ii2>ii1; ii2--){
			int chk_flag;

			current_pos -= block_size_array[ii2];
			Zt_current_pos -= evects_size_array[ii2];
			/* Block U Check */
			chk_flag = 0;
			RC_TRY(chk_elim_block33(row_size, elim_pos, next_elim_pos,
			                        current_pos, next_pos, mdbi_U, &chk_flag));
			if(chk_flag == 0){
				next_pos -= block_size_array[ii2];
				continue;
			}

			for(ii3=Zt_current_pos; ii3<Zt_row_size; ii3++){
				RC_TRY( get_row_col_matrix_db(mdbi_Zt, ii3, current_pos,
				                              &Zt_col_size, Zt_array,
				                              Zt_col_array) );
				for(ii4=0; ii4<Zt_col_size; ii4++){
					for(ii5=0; ii5<3; ii5++){
						for(ii6=0; ii6<3; ii6++){
							Zt_vect[Zt_col_array[ii4]][ii5][ii6]
							= Zt_array[ii4][ii5][ii6];
						}
					}
				}
				for(ii4=elim_pos; ii4<next_elim_pos; ii4++){
					int temp_ii4 = ii4 - elim_pos;
					RC_TRY( get_row_range_matrix_db(mdbi_U, ii4, current_pos,
					                                next_pos, &U_col_size,
					                                U_array, U_col_array) );
					/* col_sizeで繰り返し */
					for(ii5=0; ii5<U_col_size; ii5++){
						if(U_col_array[ii5] < Zt_col_array[0]) continue;
						if(U_col_array[ii5] > Zt_col_array[Zt_col_size-1]){
							break;
						}
						for(ii6=0; ii6<3; ii6++){
							for(ii7=0; ii7<3; ii7++){
								for(ii8=0; ii8<3; ii8++){
									temp_array[temp_ii4][ii6][ii7]
									+= Zt_vect[U_col_array[ii5]][ii6][ii8]
									 * U_array[ii5][ii7][ii8];
								}
							}
						}
					}
				}
				for(ii4=0; ii4<Zt_col_size; ii4++){
					for(ii5=0; ii5<3; ii5++){
						for(ii6=0; ii6<3; ii6++){
							Zt_vect[Zt_col_array[ii4]][ii5][ii6] = 0.0;
						}
					}
				}
				for(ii4=0; ii4<next_elim_pos-elim_pos; ii4++){
					if( nonzero_chk33(temp_array[ii4]) != 0 ){
					RC_TRY( add_v_matrix_db(ii3, ii4+elim_pos,
					                        temp_array[ii4], mdbi_Zt) );
					init_matrix33(temp_array[ii4]);
					}
				}
			}

		next_pos -= block_size_array[ii2];
		}
		RC_TRY( free1D33(next_elim_pos-elim_pos, &temp_array) );
		next_elim_pos -= block_size_array[ii1];
	}
	RC_TRY( free1D33(row_size, &Zt_array) );
	RC_TRY( free1I(row_size, &Zt_col_array) );
	RC_TRY( free1D33(row_size, &U_array) );
	RC_TRY( free1I(row_size, &U_col_array) );
	RC_TRY( free1D33(row_size, &Zt_vect) );


	return(NORMAL_RC);
}


static RC
chk_elim_block33 (int row_size, int current_pos, int next_pos, int elim_pos,
                  int next_elim_pos, MATRIX_DB_INFO *mdbi_K, int *chk_flag)
{
	int ii1;


	*chk_flag = 0;
	for(ii1=current_pos; ii1<next_pos; ii1++){
		int min_col = min_col_sub_matrix_db(mdbi_K, ii1, elim_pos);

		if(min_col >= next_elim_pos) continue;
		if(min_col < next_elim_pos){
			*chk_flag = 1;
			break;
		}
	}

	return(NORMAL_RC);
}


/* mdbi_B is symmetric and low tridiagonal */
/* mdbi_AtBA is symmetric and low tridiagonal */
#if 0
static RC
mul_AtBA_matrix_db (int row_size, int At_row_size, MATRIX_DB_INFO *mdbi_At,
                    MATRIX_DB_INFO *mdbi_B, MATRIX_DB_INFO *mdbi_AtBA)
{
	int ii1, ii2, ii3, ii4, ii5, ii6;
	int At_col_size;
	int B_col_size;
	int *At_col_array;
	int *B_col_array;
	double (*At_array)[3][3];
	double (*At_vect)[3][3];
	double (*B_array)[3][3];
	double (*AtB_array)[3][3];
	double (*AtBA_array)[3][3];


	RC_TRY( allocate1I(row_size, &At_col_array) );
	RC_TRY( allocate1D33(row_size, &At_array) );
	RC_TRY( allocate1D33(row_size, &At_vect) );
	RC_TRY( allocate1I(row_size, &B_col_array) );
	RC_TRY( allocate1D33(row_size, &B_array) );
	RC_TRY( allocate1D33(row_size, &AtB_array) );
	RC_TRY( allocate1D33(At_row_size, &AtBA_array) );

	for(ii1=0; ii1<At_row_size; ii1++){
		/* AtB */
		RC_TRY( get_row_matrix_db(mdbi_At,ii1, &At_col_size,
		                          At_array, At_col_array) );
		for(ii2=0; ii2<At_col_size; ii2++){
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					At_vect[At_col_array[ii2]][ii3][ii4]
					= At_array[ii2][ii3][ii4];
				}
			}
		}
		for(ii2=0; ii2<row_size; ii2++){
			RC_TRY( get_row_matrix_db(mdbi_B,ii2, &B_col_size, B_array,
			                          B_col_array) );

			for(ii3=0; ii3<B_col_size-1; ii3++){
				/* 左下3角部分 */
				if(B_col_array[ii3] >= At_col_array[0]
				&& B_col_array[ii3] <= At_col_array[At_col_size-1]){
					for(ii4=0; ii4<3; ii4++){
						for(ii5=0; ii5<3; ii5++){
							for(ii6=0; ii6<3; ii6++){
								AtB_array[ii2][ii5][ii4] +=
								B_array[ii3][ii4][ii6]
								* At_vect[B_col_array[ii3]][ii5][ii6];
							}
						}
					}
				}
				/* 右上3角部分 */
				if(ii2 < At_col_array[0]
				|| ii2 > At_col_array[At_col_size-1]) continue;
				for(ii4=0; ii4<3; ii4++){
					for(ii5=0; ii5<3; ii5++){
						for(ii6=0; ii6<3; ii6++){
							AtB_array[B_col_array[ii3]][ii5][ii4] +=
							B_array[ii3][ii6][ii4]
							* At_vect[ii2][ii5][ii6];
						}
					}
				}
			}
			/* 対角部分 */
			if(B_col_array[B_col_size-1] < At_col_array[0]
			|| B_col_array[B_col_size-1] > At_col_array[At_col_size-1]){
				continue;
			}
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					for(ii5=0; ii5<3; ii5++){
						AtB_array[ii2][ii4][ii3] +=
						B_array[B_col_size-1][ii3][ii5]
						* At_vect[B_col_array[B_col_size-1]][ii4][ii5];
					}
				}
			}
		}

		/* AtBA */
		for(ii2=0; ii2<=ii1; ii2++){
			RC_TRY( get_row_matrix_db(mdbi_At,ii2, &At_col_size, At_array,
			                          At_col_array) );
			for(ii3=0; ii3<At_col_size; ii3++){
				for(ii4=0; ii4<3; ii4++){
					for(ii5=0; ii5<3; ii5++){
						for(ii6=0; ii6<3; ii6++){
							AtBA_array[ii2][ii4][ii5] +=
							AtB_array[At_col_array[ii3]][ii4][ii6]
							* At_array[ii3][ii5][ii6];
						}
					}
				}
			}
		}
		for(ii2=0; ii2<row_size; ii2++){
			if( nonzero_chk33(AtB_array[ii2]) != 0 ){
				init_matrix33(AtB_array[ii2]);
			}
			if( nonzero_chk33(At_vect[ii2]) != 0 ){
				init_matrix33(At_vect[ii2]);
			}
		}
		for(ii2=0; ii2<=ii1; ii2++){
			if( nonzero_chk33(AtBA_array[ii2]) != 0 ){
				RC_TRY(put_v_matrix_db(ii1, ii2, AtBA_array[ii2],
				                       mdbi_AtBA));
				init_matrix33(AtBA_array[ii2]);
			}
		}
	}

	RC_TRY( free1I(row_size, &At_col_array) );
	RC_TRY( free1D33(row_size, &At_array) );
	RC_TRY( free1D33(row_size, &At_vect) );
	RC_TRY( free1I(row_size, &B_col_array) );
	RC_TRY( free1D33(row_size, &B_array) );
	RC_TRY( free1D33(row_size, &AtB_array) );
	RC_TRY( free1D33(At_row_size, &AtBA_array) );

	return(NORMAL_RC);
}
#endif


/* フルマトリクスの固有値解析 */
static RC
cal_eval_evect (int size, double **subK, double **subM,
                double *e_vals, double **e_vects)
{
	double **mat;
	double **inverse_mat;
	double **L_matrix;


	RC_TRY( allocate2D(size, size, &L_matrix) );
	RC_TRY( allocate2D(size, size, &mat) );
	RC_TRY( allocate2D(size, size, &inverse_mat) );

	RC_TRY( cholesky_decomp(size, subM, L_matrix) );
	RC_TRY( inverse_matrix(size, L_matrix, inverse_mat) );
	transpose_overwrite_matrix(size, inverse_mat);
	mul_matrix_AtBA(size, size, inverse_mat, subK, mat);
	RC_TRY( eigen_jacobi_trans(size, mat, e_vals, L_matrix) );
	mul_matrix(size, size, size, inverse_mat, L_matrix, e_vects);
	eigen_sort(size, e_vals, e_vects);

	RC_TRY( free2D(size, size, &L_matrix) );
	RC_TRY( free2D(size, size, &mat) );
	RC_TRY( free2D(size, size, &inverse_mat) );

	return(NORMAL_RC);
}


static RC
cholesky_decomp (int n, double **matrix, double **L_matrix)
{
	int ii1, ii2, ii3;


	if(n < 2) return(ARG_ERROR_RC);

	/* 正定値性，ゼロ割り判定 */
	if(matrix[0][0] < 0.0 || nearly_eq(matrix[0][0], 0.0)){
		L_matrix[0][0] = ABS_TOL;
	}else{
		L_matrix[0][0] = sqrt(matrix[0][0]);
	}

	for(ii1=1; ii1<n; ii1++){
		L_matrix[ii1][0] = matrix[ii1][0] / L_matrix[0][0];
	}
	for(ii1=1; ii1<n; ii1++){
		double root_val = matrix[ii1][ii1];
		for(ii2=0; ii2<ii1; ii2++){
			root_val -= L_matrix[ii1][ii2] * L_matrix[ii1][ii2];
		}

		/* 正定値性，ゼロ割り判定 */
		if( root_val < 0.0 || nearly_eq(root_val, 0.0) ){
			L_matrix[ii1][ii1] = ABS_TOL;
		}else{
			L_matrix[ii1][ii1] = sqrt(root_val);
		}

		for(ii2=(ii1+1); ii2<n; ii2++){
			double val = matrix[ii2][ii1];
			for(ii3=0; ii3<ii1; ii3++){
				val -= L_matrix[ii2][ii3] * L_matrix[ii1][ii3];
			}
			L_matrix[ii2][ii1] = val / L_matrix[ii1][ii1];
		}
	}

	return(NORMAL_RC);
}


static void
eigen_sort (int n, double *e_vals, double **e_vects)
{
	int ii1, ii2;
	int half = n / 2;
	double temp;

	for(ii1=0; ii1<half; ii1++){
		temp = e_vals[ii1];
		e_vals[ii1] = e_vals[n-1-ii1];
		e_vals[n-1-ii1] = temp;
		for(ii2=0; ii2<n; ii2++){
			temp = e_vects[ii2][ii1];
			e_vects[ii2][ii1] = e_vects[ii2][n-1-ii1];
			e_vects[ii2][n-1-ii1] = temp;
		}
	}
}


RC
make_UtKU_matrix_db_amls_bg (int *block_size_array, int block_num, int row_size,
                             MATRIX_DB_INFO *mdbi_K, MATRIX_DB_INFO *mdbi_U,
                             int max_col_U[])
{
	int ii1, ii2, ii3, ii4, ii5, ii6;
	/*double **block_U;
	int **zero_flag;*/
	double **temp_vect;
	int col_size;
	int elim_pos;
	int next_elim_pos;
	int *col_array;
	double (*array)[3][3];
	int nonzero_count;
	int *min_col_index;
	MATRIX_DB_INFO mdbi_G;
	MATRIX_DB_INFO mdbi_b;


	RC_TRY( allocate1D33(row_size, &array) );
	RC_TRY( allocate1I(row_size, &col_array) );
	RC_TRY( allocate1I(row_size, &min_col_index) );
	for(ii1=0; ii1<row_size; ii1++){
		RC_NEG_CHK( min_col_index[ii1] = min_col_matrix_db(mdbi_K, ii1) );
		max_col_U[ii1] = -1;
	}

	elim_pos = 0;
	next_elim_pos = 0;
	for(ii1=0; ii1<block_num-1; ii1++){
		SKYLINE_MATRIX *skyline;
		int decomp_flag = 1;
		int current_pos;
		int next_pos;
		int g_block_num;
		int *g_block_size_array;
		long total_sky_dof = 0L;
		int *block_U_index;
		int block_U_size;
		double (*block_U_array)[3][3];

		tlog_printf(4, "[%d] ii1 = %d\n", __LINE__, ii1);

		next_elim_pos += block_size_array[ii1];
		current_pos = next_elim_pos;
		next_pos = next_elim_pos;

		g_block_num = 3; /* 仮の分割数 */
		tlog_printf(4, "[%d] gb = %d\n", __LINE__, g_block_num);

		RC_TRY( allocate1I(g_block_num, &g_block_size_array) );
		RC_NULL_CHK( skyline = mm_alloc(sizeof(SKYLINE_MATRIX)* g_block_num) );
		RC_TRY( allocate1I(block_size_array[ii1], &block_U_index) );
		RC_TRY( allocate1D33(block_size_array[ii1], &block_U_array) );
		RC_TRY( open_matrix_db("./", 64*1024*1024, &mdbi_G) );
		RC_TRY( open_matrix_db("./", 64*1024*1024, &mdbi_b) );
		for(ii2=ii1+1; ii2<block_num; ii2++){
			int block_U_row_size = block_size_array[ii1] * 3;
			int chk_flag;

			next_pos += block_size_array[ii2];

			/* chk_flag = 1だと消去範囲に非ゼロ有り */
			/* chk_flag = 0だと消去範囲ゼロのみ */
			chk_flag = 0;
			for(ii3=current_pos; ii3<next_pos; ii3++){
				if(min_col_index[ii3] < next_elim_pos){
					chk_flag = 1;
					break;
				}
			}
			if(chk_flag == 0){
				current_pos += block_size_array[ii2];
				log_printf(6, "[%d] Skip [ii1:ii2]=[%d : %d]\n",
				           __LINE__, ii1, ii2);
				continue;
			}

			/* ブロックガウス消去 */
			if(decomp_flag){
#if 1
				RC_TRY( block_gauss_diag_K(block_size_array[ii1], g_block_num,
				                           g_block_size_array, row_size,
				                           elim_pos, next_elim_pos, current_pos,
				                           next_pos, min_col_index,
				                           col_array, array,
				                           mdbi_K, &mdbi_G, &mdbi_b,
				                           &skyline) );
#endif
#if 0
/* g_block_num = 1 */
RC_TRY( matrix_dbi2skyline3(row_size, elim_pos,
                            block_size_array[ii1],
                            mdbi_K, &(skyline[0])) );
RC_TRY( decomp_g_matrix(skyline[0]) );
#endif
				RC_TRY( allocate2D(3, block_size_array[ii1]*3, &temp_vect) );
				for(ii3=0; ii3<g_block_num; ii3++){
					total_sky_dof += skyline[ii3].dof;
				}
				decomp_flag = 0;
			}

			for(ii3=current_pos; ii3<next_pos; ii3++){
				if(min_col_index[ii3] >= next_elim_pos){
					continue;
				}

#if 1
				RC_TRY( get_row_range_matrix_db(&mdbi_b, ii3, elim_pos,
				                                next_elim_pos, &col_size,
				                                array, col_array) );
#endif
#if 0
/* g_block_num = 1 */
RC_TRY( get_row_range_matrix_db(mdbi_K, ii3, elim_pos,
                                next_elim_pos, &col_size,
                                array, col_array) );
#endif

				nonzero_count = 0;
				for(ii4=0; ii4<col_size; ii4++){
					int temp_ii4_3;

					nonzero_count++;
					temp_ii4_3 = (col_array[ii4] - elim_pos) * 3;
					for(ii5=0; ii5<3; ii5++){
						for(ii6=0; ii6<3; ii6++){
							temp_vect[ii5][temp_ii4_3+ii6]
							        = array[ii4][ii5][ii6];
						}
					}
				}
				if(nonzero_count == 0){
					continue;
				}

				/* ブロックガウス求解 */
#if 1
				RC_TRY( block_gauss_subst3(block_size_array[ii1], g_block_num,
				                           g_block_size_array, col_array, array,
				                           &mdbi_G, skyline, temp_vect[0],
				                           temp_vect[1], temp_vect[2]) );
#endif
#if 0
/* g_block_num = 1 */
RC_TRY( s_cholesky_subst_sp3(skyline[0].dof, skyline[0].index1,
                             skyline[0].index2, skyline[0].array,
                             temp_vect[0], temp_vect[1],
                             temp_vect[2], 1) );
#endif

				RC_TRY( vect2array33spm(total_sky_dof, temp_vect[0],
				                        temp_vect[1], temp_vect[2], elim_pos,
				                        block_U_index, block_U_array,
				                        &block_U_size) );
				for(ii4=0; ii4<block_U_row_size; ii4++){
					temp_vect[0][ii4] = 0.0;
					temp_vect[1][ii4] = 0.0;
					temp_vect[2][ii4] = 0.0;
				}

				/* ---------- UtKU ---------- */
				RC_TRY( cal_UtKU(mdbi_K, row_size, ii3, current_pos, elim_pos,
				                 next_elim_pos, col_array, array, block_U_size,
				                 block_U_index, block_U_array, min_col_index) );
				/* ----------- U ------------ */
				RC_TRY( make_U_matrix_db(mdbi_U, ii3, elim_pos, next_elim_pos,
				                         block_U_size, block_U_index,
				                         block_U_array, max_col_U) );
			}
			current_pos += block_size_array[ii2];
		}
		elim_pos += block_size_array[ii1];
		if(decomp_flag == 0){
			RC_TRY( free2D(3, block_size_array[ii1]*3, &temp_vect) );
			for(ii2=0; ii2<g_block_num; ii2++){
				RC_TRY( free_g_matrix(&(skyline[ii2])) );
			}
		}
		RC_TRY( free1I(g_block_num, &g_block_size_array) );
		RC_TRY( free1I(block_size_array[ii1], &block_U_index) );
		RC_TRY( free1D33(block_size_array[ii1], &block_U_array) );
		RC_TRY( close_matrix_db(&mdbi_G) );
		RC_TRY( close_matrix_db(&mdbi_b) );
	}
	RC_TRY( free1I(row_size, &min_col_index) );
	RC_TRY( free1I(row_size, &col_array) );
	RC_TRY( free1D33(row_size, &array) );

	return(NORMAL_RC);
}


static RC
cal_UtKU (MATRIX_DB_INFO *mdbi_K, int row_size, int c_index, int current_pos,
          int elim_pos, int next_elim_pos, int *col_array,
          double array[][3][3], int block_U_size, int *block_U_index,
          double (*block_U_array)[3][3], int min_col_index[])
{
	int ii1;
	int col_size;

	/* ------------------ KU ------------------ */
//	for(ii1=current_pos; ii1<row_size; ii1++){
	for(ii1=c_index; ii1<row_size; ii1++){
		int set_flag;
		double v[3][3];

		if(min_col_index[ii1] >= next_elim_pos) continue;

		RC_TRY( get_row_range_matrix_db(mdbi_K, ii1, elim_pos, next_elim_pos,
		                                &col_size, array, col_array) );

//		if(ii1 < c_index) continue;
		RC_TRY( inner_product33sp(col_size, array, col_array,
		                          block_U_size, block_U_array, block_U_index,
		                          v, &set_flag) );

		if(set_flag){
			RC_TRY( add_v_matrix_db(ii1, c_index, v, mdbi_K) );
		}
	}

	/* ------------------ UtK ------------------ */
	/* UtKは計算してもdeleteしか起きない */
	if(min_col_index[c_index] < next_elim_pos){
		for(ii1=min_col_index[c_index]; ii1<next_elim_pos; ii1++){
			RC_TRY( delete_v_matrix_db(c_index, ii1, mdbi_K) );
		}
	}
	RC_NEG_CHK( min_col_index[c_index] = min_col_matrix_db(mdbi_K, c_index) );

	return(NORMAL_RC);
}


static RC
make_U_matrix_db (MATRIX_DB_INFO *mdbi_U, int c_index, int elim_pos,
                  int next_elim_pos, int block_U_size,
                  int *block_U_index, double (*block_U_array)[3][3],
                  int max_col_U[])
{
	int ii1;
	int index;


	index = 0;
	for(ii1=elim_pos; ii1<next_elim_pos; ii1++){
		double vt[3][3];

//		if(index >= block_U_size) continue;
		if(index >= block_U_size) break;
		if(block_U_index[index] != ii1) continue;

		RC_TRY( transpose_matrix33( block_U_array[index], vt) );
		RC_TRY( put_v_matrix_db(ii1, c_index, vt, mdbi_U) );
		if(c_index > max_col_U[ii1]) max_col_U[ii1] = c_index;
		index++;
	}


	return(NORMAL_RC);
}


static RC
block_gauss_diag_K(int block_size, int g_block_num, 
                   int g_block_size_array[],
                   int row_size, int elim_pos, int next_elim_pos,
                   int current_pos,
                   int next_pos, int min_col_index[],
                   int col_array[], double array[][3][3],
                   MATRIX_DB_INFO *mdbi_K, MATRIX_DB_INFO *mdbi_G,
                   MATRIX_DB_INFO *mdbi_b, SKYLINE_MATRIX *skyline[])
{
	int ii1, ii2, ii3, ii4, ii5, ii6;
	int col_size;
	int g_elim_pos;
	int g_next_elim_pos;
	int nonzero_count;
	int *g_min_col_index;
	double **temp_vect;


	if(g_block_num < 1) return(ARG_ERROR_RC);

	/* mdbi_Gの初期化 */
	RC_TRY( allocate1I(block_size, &g_min_col_index) );
	for(ii1=0; ii1<block_size; ii1++){
		RC_TRY( get_row_range_matrix_db(mdbi_K, ii1+elim_pos,
		                                min_col_index[ii1+elim_pos],
		                                elim_pos+ii1+1, &col_size, array,
		                                col_array) );
#if 0
		RC_TRY( get_row_range_matrix_db(mdbi_K, ii1+elim_pos,
		                                elim_pos,
		                                next_elim_pos, &col_size, array,
		                                col_array) );
#endif
//		if(col_size <= 0) continue;
		g_min_col_index[ii1] = col_array[0] - elim_pos;
		for(ii2=0; ii2<col_size; ii2++){
			RC_TRY( put_v_matrix_db(ii1, col_array[ii2]-elim_pos, array[ii2],
			                        mdbi_G) );
		}
	}
	/* mdbi_Gの各ブロックサイズ */
	if(g_block_num == 1){
		g_block_size_array[0] = block_size;
	}else{
		RC_TRY( cal_block_array_size_amls_bg(g_block_size_array, g_block_num,
		                                     block_size, mdbi_G) );
	}
	/* 右辺ベクトル mdbi_bの初期化 */
	for(ii1=next_elim_pos; ii1<row_size; ii1++){
		if(min_col_index[ii1] >= next_elim_pos) continue;

		RC_TRY( get_row_range_matrix_db(mdbi_K, ii1, elim_pos, next_elim_pos,
		                                &col_size, array, col_array) );
		for(ii2=0; ii2<col_size; ii2++){
			RC_TRY( put_v_matrix_db(ii1, col_array[ii2], array[ii2], mdbi_b) );
		}
	}

	g_elim_pos = 0;
	g_next_elim_pos = 0;
	for(ii1=0; ii1<g_block_num-1; ii1++){
		int g_decomp_flag = 1;
		int g_current_pos;
		int g_next_pos;

		g_next_elim_pos += g_block_size_array[ii1];
		g_current_pos = g_next_elim_pos;
		g_next_pos = g_next_elim_pos;
		for(ii2=ii1+1; ii2<g_block_num; ii2++){
			int g_block_U_row_size = g_block_size_array[ii1] * 3;
			int chk_flag;
			int *g_block_U_size;
			int **g_block_U_index;
			double (**g_block_U_array)[3][3];

			g_next_pos += g_block_size_array[ii2];

			/* 対角ブロックLU分解(求解時必須) */
			if(g_decomp_flag){
				RC_TRY( matrix_dbi2skyline3(block_size, g_elim_pos, 
				                            g_block_size_array[ii1],
				                            mdbi_G, &((*skyline)[ii1])) );
				RC_TRY( decomp_g_matrix((*skyline)[ii1]) );
				RC_TRY( allocate2D(3, g_block_size_array[ii1]*3, &temp_vect) );
				g_decomp_flag = 0;
			}
			/* chk_flag = 1だと消去範囲に非ゼロ有り */
			/* chk_flag = 0だと消去範囲ゼロのみ */
			chk_flag = 0;
			for(ii3=g_current_pos; ii3<g_next_pos; ii3++){
				if(g_min_col_index[ii3] < g_next_elim_pos){
					chk_flag = 1;
					break;
				}
			}
			if(chk_flag == 0){
				g_current_pos += g_block_size_array[ii2];
				log_printf(6, "[%d] Skip [ii1:ii2]=[%d : %d]\n",
				           __LINE__, ii1, ii2);
				continue;
			}

			RC_NULL_CHK( g_block_U_array = mm_alloc(sizeof(double(*)[3][3])
			                                       *g_block_size_array[ii2]) );
			RC_NULL_CHK( g_block_U_index = mm_alloc(sizeof(int *)
			                                       *g_block_size_array[ii2]) );
			RC_NULL_CHK( g_block_U_size = mm_alloc(sizeof(int)
			                                      *g_block_size_array[ii2]) );

			for(ii3=g_current_pos; ii3<g_next_pos; ii3++){
				if(g_min_col_index[ii3] >= g_next_elim_pos){
					g_block_U_size[ii3-g_current_pos] = 0;
					g_block_U_index[ii3-g_current_pos] = NULL;
					g_block_U_array[ii3-g_current_pos] = NULL;
					continue;
				}

				RC_TRY( get_row_range_matrix_db(mdbi_G, ii3, g_elim_pos,
				                g_next_elim_pos, &col_size, array, col_array) );

				nonzero_count = 0;
				for(ii4=0; ii4<col_size; ii4++){
					int temp_ii4_3;

					nonzero_count++;
					temp_ii4_3 = (col_array[ii4] - g_elim_pos) * 3;
					for(ii5=0; ii5<3; ii5++){
						for(ii6=0; ii6<3; ii6++){
							temp_vect[ii5][temp_ii4_3+ii6]
							        = array[ii4][ii5][ii6];
						}
					}
				}
				if(nonzero_count == 0){
					g_block_U_size[ii3-g_current_pos] = 0;
					g_block_U_index[ii3-g_current_pos] = NULL;
					g_block_U_array[ii3-g_current_pos] = NULL;
					continue;
				}

				RC_TRY( s_cholesky_subst_sp3((*skyline)[ii1].dof,
				                             (*skyline)[ii1].index1,
				                             (*skyline)[ii1].index2,
				                             (*skyline)[ii1].array,
				                             temp_vect[0], temp_vect[1],
				                             temp_vect[2], 1) );

				RC_TRY( allocate1I(g_block_size_array[ii1],
				                   &g_block_U_index[ii3-g_current_pos]) );
				RC_TRY( allocate1D33(g_block_size_array[ii1],
				                     &g_block_U_array[ii3-g_current_pos]) );

				RC_TRY( vect2array33spm((*skyline)[ii1].dof, temp_vect[0],
				                        temp_vect[1], temp_vect[2], g_elim_pos,
				                        g_block_U_index[ii3-g_current_pos],
				                        g_block_U_array[ii3-g_current_pos],
				                        &g_block_U_size[ii3-g_current_pos]) );

				for(ii4=0; ii4<g_block_U_row_size; ii4++){
					temp_vect[0][ii4] = 0.0;
					temp_vect[1][ii4] = 0.0;
					temp_vect[2][ii4] = 0.0;
				}
				if(g_block_U_size[ii3-g_current_pos] <= 0){
					RC_TRY( mm_free(g_block_U_index[ii3-g_current_pos]) );
					g_block_U_index[ii3-g_current_pos] = NULL;
					RC_TRY( mm_free(g_block_U_array[ii3-g_current_pos]) );
					g_block_U_array[ii3-g_current_pos] = NULL;
				}else{
					RC_NULL_CHK( g_block_U_index[ii3-g_current_pos]
					= mm_realloc(g_block_U_index[ii3-g_current_pos],
					g_block_U_size[ii3-g_current_pos]*sizeof(int)) );
					RC_NULL_CHK(g_block_U_array[ii3-g_current_pos]
					= mm_realloc(g_block_U_array[ii3-g_current_pos],
					g_block_U_size[ii3-g_current_pos]*sizeof(double[3][3])));
				}
			}
			RC_TRY( cal_UtL(mdbi_G, block_size, g_current_pos, g_next_pos,
			                g_elim_pos, g_next_elim_pos, col_array, array,
			                g_block_U_size, g_block_U_index, g_block_U_array,
			                g_min_col_index));
			for(ii3=g_current_pos; ii3<g_next_pos; ii3++){
				for(ii4=0; ii4<g_block_U_size[ii3-g_current_pos];ii4++){
					g_block_U_index[ii3-g_current_pos][ii4] += elim_pos;
				}
			}
			RC_TRY( cal_Utb(mdbi_b, row_size, current_pos, elim_pos,
			                g_current_pos, g_next_pos, g_elim_pos,
			                g_next_elim_pos, col_array, array, g_block_U_size,
			                g_block_U_index, g_block_U_array,
			                g_min_col_index, min_col_index) );

			for(ii3=g_current_pos; ii3<g_next_pos; ii3++){
				if(g_block_U_size[ii3-g_current_pos] > 0){
					RC_TRY( mm_free(g_block_U_index[ii3-g_current_pos]) );
					RC_TRY( mm_free(g_block_U_array[ii3-g_current_pos]) );
				}
			}
			RC_TRY( mm_free(g_block_U_array) );
			RC_TRY( mm_free(g_block_U_index) );
			RC_TRY( mm_free(g_block_U_size) );

			g_current_pos += g_block_size_array[ii2];
		}

		g_elim_pos += g_block_size_array[ii1];
		if(g_decomp_flag == 0){
			RC_TRY( free2D(3, g_block_size_array[ii1]*3, &temp_vect) );
		}
	}
	RC_TRY( free1I(block_size, &g_min_col_index) );

	/* 最終対角ブロックLU分解 */
	RC_TRY( matrix_dbi2skyline3(block_size, g_elim_pos, 
	                            g_block_size_array[g_block_num-1],
	                            mdbi_G, &((*skyline)[g_block_num-1])) );
	RC_TRY( decomp_g_matrix((*skyline)[g_block_num-1]) );

	return(NORMAL_RC);
}


static RC
block_gauss_subst3(int block_size, int g_block_num,
                   int g_block_size_array[], int col_array[],
                   double array[][3][3], MATRIX_DB_INFO *mdbi_G,
                   SKYLINE_MATRIX skyline[],
                   double v0[], double v1[], double v2[])
{
	int ii1, ii2, ii3, ii4;
	int ii2x;
	int col_size;
	int sol_pos;
	int old_pos;
	int *min_index;
	double *sub_v0;
	double *sub_v1;
	double *sub_v2;


	RC_TRY( allocate1I(block_size-1, &min_index) );
	RC_TRY( allocate1D(block_size * 3, &sub_v0) );
	RC_TRY( allocate1D(block_size * 3, &sub_v1) );
	RC_TRY( allocate1D(block_size * 3, &sub_v2) );


	/* 後退代入 */
	sol_pos = block_size;
	old_pos = block_size;
	for(ii1=g_block_num-1; ii1>=0; ii1--){
		int col_pos;
		int old_col_pos;

		sol_pos -= g_block_size_array[ii1];
		for(ii2=sol_pos,ii2x=0; ii2<old_pos; ii2++,ii2x++){
			int ii2_3 = ii2 * 3;
			int ii2x_3 = ii2x * 3;

			sub_v0[ii2x_3]   = v0[ii2_3];
			sub_v0[ii2x_3+1] = v0[ii2_3+1];
			sub_v0[ii2x_3+2] = v0[ii2_3+2];
			sub_v1[ii2x_3]   = v1[ii2_3];
			sub_v1[ii2x_3+1] = v1[ii2_3+1];
			sub_v1[ii2x_3+2] = v1[ii2_3+2];
			sub_v2[ii2x_3]   = v2[ii2_3];
			sub_v2[ii2x_3+1] = v2[ii2_3+1];
			sub_v2[ii2x_3+2] = v2[ii2_3+2];
		}
#if 0
		for(ii2=old_pos; ii2<block_size; ii2++){
			RC_NEG_CHK( min_index[ii2-old_pos] = min_col_sub_matrix_db(mdbi_G,
			                                     ii2, sol_pos) );
		}
#endif

		col_pos = block_size;
		old_col_pos = block_size;
		for(ii2=g_block_num-1; ii2>=ii1+1; ii2--){
			col_pos -= g_block_size_array[ii2];
			for(ii3=col_pos; ii3<old_col_pos; ii3++){
				int ii3_3 = ii3 * 3;

//				if(min_index[ii3-old_pos] >= old_pos) continue;
				RC_TRY( get_row_range_matrix_db(mdbi_G, ii3, sol_pos, old_pos,
		                                        &col_size, array, col_array) );
				for(ii4=0; ii4<col_size; ii4++){
					int ii4_col3 = (col_array[ii4] - sol_pos) * 3;
					double (*array_v)[3] = array[ii4];

					sub_v0[ii4_col3]   -= array_v[0][0] * v0[ii3_3]
					                    + array_v[1][0] * v0[ii3_3+1]
					                    + array_v[2][0] * v0[ii3_3+2];
					sub_v0[ii4_col3+1] -= array_v[0][1] * v0[ii3_3]
					                    + array_v[1][1] * v0[ii3_3+1]
					                    + array_v[2][1] * v0[ii3_3+2];
					sub_v0[ii4_col3+2] -= array_v[0][2] * v0[ii3_3]
					                    + array_v[1][2] * v0[ii3_3+1]
					                    + array_v[2][2] * v0[ii3_3+2];

					sub_v1[ii4_col3]   -= array_v[0][0] * v1[ii3_3]
					                    + array_v[1][0] * v1[ii3_3+1]
					                    + array_v[2][0] * v1[ii3_3+2];
					sub_v1[ii4_col3+1] -= array_v[0][1] * v1[ii3_3]
					                    + array_v[1][1] * v1[ii3_3+1]
					                    + array_v[2][1] * v1[ii3_3+2];
					sub_v1[ii4_col3+2] -= array_v[0][2] * v1[ii3_3]
					                    + array_v[1][2] * v1[ii3_3+1]
					                    + array_v[2][2] * v1[ii3_3+2];

					sub_v2[ii4_col3]   -= array_v[0][0] * v2[ii3_3]
					                    + array_v[1][0] * v2[ii3_3+1]
					                    + array_v[2][0] * v2[ii3_3+2];
					sub_v2[ii4_col3+1] -= array_v[0][1] * v2[ii3_3]
					                    + array_v[1][1] * v2[ii3_3+1]
					                    + array_v[2][1] * v2[ii3_3+2];
					sub_v2[ii4_col3+2] -= array_v[0][2] * v2[ii3_3]
					                    + array_v[1][2] * v2[ii3_3+1]
					                    + array_v[2][2] * v2[ii3_3+2];
				}
			}
			old_col_pos = col_pos;
		}
		RC_TRY(s_cholesky_subst_sp3(skyline[ii1].dof, skyline[ii1].index1,
		                            skyline[ii1].index2, skyline[ii1].array,
		                            sub_v0, sub_v1, sub_v2, 1));

		for(ii2=sol_pos,ii2x=0; ii2<old_pos; ii2++,ii2x++){
			int ii2_3 = ii2 * 3;
			int ii2x_3 = ii2x * 3;

			v0[ii2_3]   = sub_v0[ii2x_3];
			v0[ii2_3+1] = sub_v0[ii2x_3+1];
			v0[ii2_3+2] = sub_v0[ii2x_3+2];
			v1[ii2_3]   = sub_v1[ii2x_3];
			v1[ii2_3+1] = sub_v1[ii2x_3+1];
			v1[ii2_3+2] = sub_v1[ii2x_3+2];
			v2[ii2_3]   = sub_v2[ii2x_3];
			v2[ii2_3+1] = sub_v2[ii2x_3+1];
			v2[ii2_3+2] = sub_v2[ii2x_3+2];

			sub_v0[ii2x_3]   = 0.0;
			sub_v0[ii2x_3+1] = 0.0;
			sub_v0[ii2x_3+2] = 0.0;
			sub_v1[ii2x_3]   = 0.0;
			sub_v1[ii2x_3+1] = 0.0;
			sub_v1[ii2x_3+2] = 0.0;
			sub_v2[ii2x_3]   = 0.0;
			sub_v2[ii2x_3+1] = 0.0;
			sub_v2[ii2x_3+2] = 0.0;
		}
		old_pos = sol_pos;
	}
	RC_TRY( free1I(block_size-1, &min_index) );
	RC_TRY( free1D(block_size * 3, &sub_v0) );
	RC_TRY( free1D(block_size * 3, &sub_v1) );
	RC_TRY( free1D(block_size * 3, &sub_v2) );

	return(NORMAL_RC);
}


static RC
cal_UtL (MATRIX_DB_INFO *mdbi_G, int block_size, int g_current_pos,
         int g_next_pos, int g_elim_pos, int g_next_elim_pos,
         int *col_array, double array[][3][3], int g_block_U_size[],
         int **g_block_U_index, double (**g_block_U_array)[3][3],
         int g_min_col_index[])
{
	int ii1, ii2;
	int col_size;

	/* ------------------ UtL ------------------ */
	for(ii1=g_current_pos; ii1<block_size; ii1++){
		if(g_min_col_index[ii1] >= g_next_elim_pos) continue;
		RC_TRY( get_row_range_matrix_db(mdbi_G, ii1, g_elim_pos,
		                                g_next_elim_pos, &col_size, array,
		                                col_array) );
		for(ii2=g_current_pos; ii2<g_next_pos; ii2++){
			int set_flag;
			double v[3][3];

			if(ii1 < ii2) continue;
			RC_TRY( inner_product33sp(col_size, array, col_array,
			                          g_block_U_size[ii2-g_current_pos],
			                          g_block_U_array[ii2-g_current_pos],
			                          g_block_U_index[ii2-g_current_pos], v,
			                          &set_flag) );

			if(set_flag){
				RC_TRY( add_v_matrix_db(ii1, ii2, v, mdbi_G) );
			}
		}
	}

	/* ダミー消去 */
	for(ii1=g_current_pos; ii1<g_next_pos; ii1++){
		if(g_min_col_index[ii1] >= g_next_elim_pos) continue;
		RC_NEG_CHK( g_min_col_index[ii1] = min_col_sub_matrix_db(mdbi_G, ii1,
		                                   g_next_elim_pos) );
	}

	return(NORMAL_RC);
}


static RC
cal_Utb (MATRIX_DB_INFO *mdbi_b, int row_size, int current_pos, int elim_pos,
         int g_current_pos, int g_next_pos, int g_elim_pos,
         int g_next_elim_pos, int *col_array, double array[][3][3],
         int g_block_U_size[], int **g_block_U_index,
         double (**g_block_U_array)[3][3],
         int g_min_col_index[], int min_col_index[])
{
	int ii1, ii2;
	int col_size;
	int col_pos0 = elim_pos + g_elim_pos;
	int col_pos_max = elim_pos + g_next_elim_pos;


	for(ii1=current_pos; ii1<row_size; ii1++){
		if(min_col_index[ii1] >= col_pos_max) continue;

		RC_TRY( get_row_range_matrix_db(mdbi_b, ii1, col_pos0,
		                                col_pos_max, &col_size, array,
		                                col_array) );
		for(ii2=g_current_pos; ii2<g_next_pos; ii2++){
			int set_flag;
			double v[3][3];
			double vt[3][3];

			RC_TRY( inner_product33sp(g_block_U_size[ii2-g_current_pos],
			                          g_block_U_array[ii2-g_current_pos],
			                          g_block_U_index[ii2-g_current_pos],
			                          col_size, array, col_array,
			                          v, &set_flag) );

			if(set_flag){
				RC_TRY( transpose_matrix33(v, vt) );
				RC_TRY( add_v_matrix_db(ii1, elim_pos+ii2, vt, mdbi_b) );
			}
		}
	}

	return(NORMAL_RC);
}


static RC
s_cholesky_subst_sp3 (long n, const long *index1, const long *index2,
                      const double *array, double *vector0,
                      double *vector1, double *vector2, int p_flag)
{
	long ii1, ii2;
	long nz_begin;
	long nz_end;

	if((n < 1L)||((void *)array == NULL)
			||((void *)index1 == NULL)
			||((void *)index2 == NULL)
			||((void *)vector0 == NULL)
			||((void *)vector1 == NULL)
			||((void *)vector2 == NULL)) return(ARG_ERROR_RC);

	/* Set nonzero begin index nz_begin */
	/* Set nonzero end index nz_end */
	nz_begin = n;
	for(ii1=0L; ii1<n; ii1++){
		if( (fabs(vector0[ii1]) > ABS_TOL) 
		  ||(fabs(vector1[ii1]) > ABS_TOL) 
		  ||(fabs(vector2[ii1]) > ABS_TOL) ){
			nz_begin = ii1;
			break;
		}
	}
	nz_end = 0L;
	for(ii1=n-1L; ii1>=nz_begin; ii1--){
		if( (fabs(vector0[ii1]) > ABS_TOL) 
		  ||(fabs(vector1[ii1]) > ABS_TOL) 
		  ||(fabs(vector2[ii1]) > ABS_TOL) ){
			nz_end = ii1;
			break;
		}
	}

	if(p_flag) RC_TRY( log_printf(6, "[B]") );
	for(ii1=nz_begin;ii1<n;ii1++){
		long ii2_start = index1[ii1];
		long ii2_finish = ii1;

		if(ii2_start < nz_begin) ii2_start = nz_begin;
		if(ii2_finish > nz_end + 1L) ii2_finish = nz_end + 1L;
		if(ii2_start < ii2_finish){
			double vector0_ii1 = vector0[ii1];
			double vector1_ii1 = vector1[ii1];
			double vector2_ii1 = vector2[ii1];

			for(ii2=ii2_start;ii2<ii2_finish;ii2++){
				double array_v = array[index2[ii1] - ii1 + ii2];

				vector0_ii1 -= vector0[ii2] * array_v;
				vector1_ii1 -= vector1[ii2] * array_v;
				vector2_ii1 -= vector2[ii2] * array_v;
			}
			vector0[ii1] = vector0_ii1;
			vector1[ii1] = vector1_ii1;
			vector2[ii1] = vector2_ii1;
			if(ii1 < nz_begin) nz_begin = ii1;
			if(ii1 > nz_end) nz_end = ii1;
		}
	}

	for(ii1=nz_begin;ii1<=nz_end;ii1++){
		if(fabs(array[index2[ii1]]) < ABS_TOL) return(CAL_ERROR_RC);
		vector0[ii1] /= array[index2[ii1]];
		vector1[ii1] /= array[index2[ii1]];
		vector2[ii1] /= array[index2[ii1]];
	}

	for(ii1=nz_end;ii1>=nz_begin;ii1--){
		double vector0_ii1 = vector0[ii1];
		double vector1_ii1 = vector1[ii1];
		double vector2_ii1 = vector2[ii1];

		for(ii2=index1[ii1];ii2<ii1;ii2++){
			double array_v = array[index2[ii1] - ii1 + ii2];

			vector0[ii2] -= vector0_ii1 * array_v;
			vector1[ii2] -= vector1_ii1 * array_v;
			vector2[ii2] -= vector2_ii1 * array_v;
		}
		if(index1[ii1] < nz_begin) nz_begin = index1[ii1];
	}
	if(p_flag) RC_TRY( log_printf(6, "...finished.\n") );

	return(NORMAL_RC);
}


/* v0[v_size], v1[v_size], v2[v_size] のマイナス値を */
/* ブロック圧縮（v_size は３の倍数）*/
/* offset : ブロックインデックスのオフセット */
static RC
vect2array33spm (int v_size, const double v0[], const double v1[],
                 const double v2[], int offset, int index[],
                 double array[][3][3], int *size)
{
	int ii1;

	*size = 0;
	for(ii1=0; ii1<v_size; ii1+=3){
		array[*size][0][0] = -v0[ii1    ];
		array[*size][0][1] = -v0[ii1 + 1];
		array[*size][0][2] = -v0[ii1 + 2];
		array[*size][1][0] = -v1[ii1    ];
		array[*size][1][1] = -v1[ii1 + 1];
		array[*size][1][2] = -v1[ii1 + 2];
		array[*size][2][0] = -v2[ii1    ];
		array[*size][2][1] = -v2[ii1 + 1];
		array[*size][2][2] = -v2[ii1 + 2];
		if(nonzero_chk33(array[*size])){
			index[*size] = offset + ii1/3;
			(*size)++;
		}
	}

	return(NORMAL_RC);
}


static RC
inner_product33sp (int size1, double array1[][3][3], const int index1[],
                   int size2, double array2[][3][3], const int index2[],
                   double v[3][3], int *set_flag)
{
	int ii1 = 0;
	int ii2 = 0;

	*set_flag = 0;
	init_matrix33(v);
	while((ii1 < size1)&&(ii2 < size2)){
		if(index1[ii1] == index2[ii2]){
			double (*v1)[3] = array1[ii1];
			double (*v2)[3] = array2[ii2];

			v[0][0] += v1[0][0]*v2[0][0]+v1[0][1]*v2[0][1]+v1[0][2]*v2[0][2];
			v[0][1] += v1[0][0]*v2[1][0]+v1[0][1]*v2[1][1]+v1[0][2]*v2[1][2];
			v[0][2] += v1[0][0]*v2[2][0]+v1[0][1]*v2[2][1]+v1[0][2]*v2[2][2];
			v[1][0] += v1[1][0]*v2[0][0]+v1[1][1]*v2[0][1]+v1[1][2]*v2[0][2];
			v[1][1] += v1[1][0]*v2[1][0]+v1[1][1]*v2[1][1]+v1[1][2]*v2[1][2];
			v[1][2] += v1[1][0]*v2[2][0]+v1[1][1]*v2[2][1]+v1[1][2]*v2[2][2];
			v[2][0] += v1[2][0]*v2[0][0]+v1[2][1]*v2[0][1]+v1[2][2]*v2[0][2];
			v[2][1] += v1[2][0]*v2[1][0]+v1[2][1]*v2[1][1]+v1[2][2]*v2[1][2];
			v[2][2] += v1[2][0]*v2[2][0]+v1[2][1]*v2[2][1]+v1[2][2]*v2[2][2];
			ii1++;
			ii2++;
			*set_flag = 1;
		}else if(index1[ii1] < index2[ii2]){
			ii1++;
		}else{   /* index1[ii1] > index2[ii2] */
			ii2++;
		}
	}

	return(NORMAL_RC);
}


static RC
matrix_dbi2skyline3 (int row_size, int row0_index, int block_size,
                     MATRIX_DB_INFO *mdbi, SKYLINE_MATRIX *skyline)
{
	int ii1, ii2, ii3, ii4;
	int *col_array;
	double (*array)[3][3];


	skyline->dim = 3;
	skyline->dof = block_size * 3; 
	RC_NULL_CHK( skyline->index1 = (long *)mm_alloc(skyline->dof*sizeof(long)));

	log_printf(6, "block_size = %d\n", block_size);
	/* fill index1[] */
	for(ii1=0; ii1<block_size; ii1++){
		int col0_index;
		int row_col0 = min_col_sub_matrix_db(mdbi, row0_index+ii1, row0_index);
		log_printf(6, "row_col0 = %d\n", row_col0);
		if(row_col0 < 0) return(ARG_ERROR_RC);

		col0_index = row_col0 - row0_index;
		skyline->index1[3*ii1 + 0] = 3 * col0_index;
		skyline->index1[3*ii1 + 1] = 3 * col0_index;
		skyline->index1[3*ii1 + 2] = 3 * col0_index;
	}

	RC_TRY( s_cholesky_alloc(skyline->dof, skyline->index1, &(skyline->index2),
	                         &(skyline->array), &(skyline->array_size)) );

	RC_TRY( allocate1I(row_size, &col_array) );
	RC_TRY( allocate1D33(row_size, &array) );
	for(ii1=0; ii1<block_size; ii1++){
		int size;

		RC_TRY( get_row_col_matrix_db(mdbi, row0_index+ii1, row0_index,
		                              &size, array, col_array) );

		for(ii2=0; ii2<size; ii2++){
			int row3 = 3*ii1;
			int col3;
			int dif_col = col_array[ii2] - row0_index;

			if(dif_col < 0) continue;
			if(dif_col >= block_size) break;
			col3 = 3*(col_array[ii2] - row0_index);

			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					int row = row3 + ii3;
					int col = col3 + ii4;

					if(row < col) continue;

					skyline->array[skyline->index2[row] - row + col]
					        = array[ii2][ii3][ii4];
				}
			}
		}
	}
	RC_TRY( free1I(row_size, &col_array) );
	RC_TRY( free1D33(row_size, &array) );

	return(NORMAL_RC);
}


/* 下三角で記憶 */
RC
make_KM_matrix_db_amls_bg (NODE_ARRAY node, ELEMENT_ARRAY element,
                           MATERIAL_PROP_ARRAY material,
                           PHYSICAL_PROP_ARRAY physical,
                           MATRIX_DB_INFO *mdbi_K, MATRIX_DB_INFO *mdbi_M)
{
	int ii1, ii2, ii3, ii4, ii5;

	for(ii1=0; ii1<element.size; ii1++){
		double v[3][3];
		ELEM_MATRIX elem_matrix;

		if(element.array[ii1].label < 0) continue;
		
		RC_TRY( make_elem_matrix(element.array[ii1], node,
		                         material, physical, &elem_matrix) );
		for(ii2=0; ii2<elem_matrix.size; ii2+=3){
			for(ii3=0; ii3<elem_matrix.size; ii3+=3){
				int ix2 = elem_matrix.index[ii2] / 3;
				int ix3 = elem_matrix.index[ii3] / 3;

				if(ix2 < ix3) continue;

				for(ii4=0; ii4<3; ii4++){
					for(ii5=0; ii5<3; ii5++){
						v[ii4][ii5] = elem_matrix.matrix[ii2+ii4][ii3+ii5];
					}
				}

				RC_TRY( add_v_matrix_db(ix2, ix3, v, mdbi_K) );
			}
		}
		RC_TRY( free_elem_matrix(&elem_matrix) );

		RC_TRY( make_mass_matrix(element.array[ii1], node,
		                         material, physical, &elem_matrix) );
		/* 要素集中質量 */
		for(ii2=0; ii2<elem_matrix.size; ii2++){
			double lump_val = 0.0;
			for(ii3=0; ii3<elem_matrix.size; ii3++){
				lump_val += elem_matrix.matrix[ii2][ii3];
				elem_matrix.matrix[ii2][ii3] = 0.0;
			}
			elem_matrix.matrix[ii2][ii2] = lump_val;
		}
		init_matrix33(v);
		for(ii2=0; ii2<elem_matrix.size; ii2+=3){
			int ix2 = elem_matrix.index[ii2] / 3;

			for(ii3=0; ii3<3; ii3++){
				v[ii3][ii3] = elem_matrix.matrix[ii2+ii3][ii2+ii3];
			}
			RC_TRY( add_v_matrix_db(ix2, ix2, v, mdbi_M) );
		}
		RC_TRY( free_elem_matrix(&elem_matrix) );
	}

	return(NORMAL_RC);
}


RC
cal_block_array_size_amls_bg (int *block_size_array, int block_num,
                              int row_size, MATRIX_DB_INFO *mdbi_K)
{
	int ii1, ii2;
	int col_size;
	int *col_array;
	int *nonzero_size;
	double (*array)[3][3];
	WRAP64_INT sum;
	int block_avg;
	int nonzero_avg;
	int nonzero_num;
	int diag_weight;

	if(row_size < 1) return(ARG_ERROR_RC);
	if(row_size < block_num) return(ARG_ERROR_RC);
	if(block_num < 2) return(ARG_ERROR_RC);
	/* 多すぎるブロックサイズは危険 */
	if(block_num > row_size/2) return(ARG_ERROR_RC);
	RC_NULL_CHK( block_size_array );
	RC_NULL_CHK( mdbi_K );

	RC_TRY( allocate1I(row_size, &nonzero_size) );

	RC_TRY( allocate1I(row_size, &col_array) );
	RC_TRY( allocate1D33(row_size, &array) );
	for(ii1=0, sum=0; ii1<row_size; ii1++){
		RC_TRY( get_row_matrix_db(mdbi_K, ii1, &col_size, array, col_array) );
		if(col_size <= 0) return(CAL_ERROR_RC);
		nonzero_size[ii1] += col_size;
		for(ii2=0; ii2<(col_size - 1); ii2++){
			nonzero_size[col_array[ii2]]++;
		}
		sum += col_size + col_size - 1;
	}
	RC_TRY( free1I(row_size, &col_array) );
	RC_TRY( free1D33(row_size, &array) );

	block_size_array[0] = 0;
	for(ii1=0; ii1<row_size; ii1++){
		if(nonzero_size[ii1] == 1){
			(block_size_array[0])++;
		}else{
			break;
		}
	}
	if(block_size_array[0] <= 2){
		block_size_array[0] = 0;
	}
	sum -= block_size_array[0];

	diag_weight = (int)(sum/(row_size - block_size_array[0]));
	nonzero_avg
	    = (int)((sum + (row_size - block_size_array[0])*diag_weight)/row_size);
	block_avg
	    = (int)((sum + (row_size - block_size_array[0])*diag_weight)
	             /(block_num - 1));

	ii2 = 1;
	if(block_size_array[0] == 0){
		ii2 = 0;
	}
	for(ii1=block_size_array[0], nonzero_num=0; ii1<row_size; ii1++){
		if(ii2+1 >= block_num){
			block_size_array[ii2]++;
			continue;
		}
		if((nonzero_num + nonzero_avg + diag_weight) < block_avg){
			nonzero_num += nonzero_size[ii1] + diag_weight;
			block_size_array[ii2]++;
			continue;
		}
		if(block_size_array[ii2] == 0){
			block_size_array[ii2]++;
			nonzero_num += nonzero_size[ii1] + diag_weight;
			continue;
		}
		ii2++;
		block_size_array[ii2]++;
		nonzero_num = nonzero_size[ii1] + diag_weight;
	}

	{
		int size_sum = 0;

		for(ii1=0; ii1<block_num; ii1++){
			RC_TRY( log_printf(5, "block_size_array[%4d] = %3d\n", ii1, 
			                   block_size_array[ii1]) );
			size_sum += block_size_array[ii1];
		}
		RC_TRY( log_printf(5, "block size sum = %3d\n", size_sum) );
	}

	RC_TRY( free1I(row_size, &nonzero_size) );

	return(NORMAL_RC);
}


/* A * vect[vect_num] => ans[vect_num] */
/* vect_num <= MAX_VECT_NUM */
#define MAX_VECT_NUM   (20)
static RC
mul_matrix_vect_db_sp33(int row_size, MATRIX_DB_INFO *A, const int min_col_A[],
                int vect_num, const int vect_size[], int **vect_row,
                double (**vect)[3][3], double (**ans)[3][3])
{
	int ii1, ii2, ii3;
	int A_size;
	double (*A_array)[3][3];
	int *A_col_array;
	int rp_vect[MAX_VECT_NUM];
	int min_vect_row = row_size - 1;
	int max_vect_row = 0;

	if(vect_num > MAX_VECT_NUM) return(ARG_ERROR_RC);

	for(ii1=0; ii1<vect_num; ii1++){
		for(ii2=0; ii2<row_size; ii2++){
			init_matrix33(ans[ii1][ii2]);
		}
		if(vect_row[ii1][0] < min_vect_row){
			min_vect_row = vect_row[ii1][0];
		}
		if(vect_row[ii1][vect_size[ii1]-1] > max_vect_row){
			max_vect_row = vect_row[ii1][vect_size[ii1]-1];
		}
		if(vect_size[ii1] <= 0) return(NORMAL_RC);
	}

	RC_TRY( allocate1D33(row_size, &A_array) );
	RC_TRY( allocate1I(row_size, &A_col_array) );

	for(ii1=0; ii1<vect_num; ii1++) rp_vect[ii1] = 0;

	for(ii1=min_vect_row; ii1<row_size; ii1++){
		double v[3][3];
		int set_flag;

		if(max_vect_row < min_col_A[ii1]) continue;

		RC_TRY( get_row_matrix_db(A, ii1, &A_size, A_array, A_col_array) );
		for(ii2=0; ii2<vect_num; ii2++){
			RC_TRY( inner_product33sp(A_size, A_array, A_col_array,
			                          vect_size[ii2], vect[ii2], vect_row[ii2],
			                          v, &set_flag) );
			if(set_flag) add_matrix33t(ans[ii2][ii1], v);

			if(rp_vect[ii2] >= vect_size[ii2]) continue;
			if(ii1 != vect_row[ii2][rp_vect[ii2]]) continue;

			for(ii3=0; ii3<(A_size-1); ii3++){
				RC_TRY( mul_add_matrix33(vect[ii2][rp_vect[ii2]], A_array[ii3],
				                         ans[ii2][A_col_array[ii3]]) );
			}
			(rp_vect[ii2])++;
		}
	}

	RC_TRY( free1D33(row_size, &A_array) );
	RC_TRY( free1I(row_size, &A_col_array) );

	return(NORMAL_RC);
}


/* mdbi_B is symmetric and low tridiagonal */
static RC
mul_AtBA_matrix_db(int row_size, int At_row_size, MATRIX_DB_INFO *mdbi_At,
                   MATRIX_DB_INFO *mdbi_B, double **AtBA)
{
	int ii1, ii2, ii3, ii4;
	int At_size[MAX_VECT_NUM];
	int At_size2;
	int *min_col_B;
	double (**BA)[3][3];
	double (**At_array)[3][3];
	int **At_col_array;
	double (*At_array2)[3][3];
	int *At_col_array2;
	double v[3][3];
	int At_num;

	RC_TRY( allocate1I(row_size, &min_col_B) );
	RC_TRY( allocate2D33(MAX_VECT_NUM, row_size, &BA) );
	RC_TRY( allocate2D33(MAX_VECT_NUM, row_size, &At_array) );
	RC_TRY( allocate2I(MAX_VECT_NUM, row_size, &At_col_array) );
	RC_TRY( allocate1D33(row_size, &At_array2) );
	RC_TRY( allocate1I(row_size, &At_col_array2) );

	for(ii1=0; ii1<row_size; ii1++){
		RC_NEG_CHK( min_col_B[ii1] = min_col_matrix_db(mdbi_B, ii1) );
	}

	init_matrix(3*At_row_size, 3*At_row_size, AtBA);
	for(ii1=0; ii1<At_row_size; ii1+=MAX_VECT_NUM){
		At_num = 0;

		for(ii2=0; ii2<MAX_VECT_NUM; ii2++){
			if(ii1+ii2 >= At_row_size) break;
			RC_TRY( get_row_matrix_db(mdbi_At, ii1 + ii2, &(At_size[ii2]),
			                          At_array[ii2], At_col_array[ii2]) );
			At_num++;
		}

		RC_TRY( mul_matrix_vect_db_sp33(row_size, mdbi_B, min_col_B, At_num,
		                           At_size, At_col_array, At_array, BA) );

		for(ii2=0; ii2<(ii1+At_num); ii2++){
			RC_TRY( get_row_matrix_db(mdbi_At,ii2, &At_size2,
			                          At_array2, At_col_array2) );
			for(ii3=0; ii3<At_num; ii3++){
				int index = ii1+ii3;
				if(index < ii2) continue;

				init_matrix33(v);
				for(ii4=0; ii4<At_size2; ii4++){
					RC_TRY( mul_add_matrix33_ABt(BA[ii3][At_col_array2[ii4]],
					                             At_array2[ii4], v) );
				}
				if(ii2 < index){
					for(ii4=0; ii4<3; ii4++){
						AtBA[3*index+ii4][3*ii2  ] = v[ii4][0];
						AtBA[3*index+ii4][3*ii2+1] = v[ii4][1];
						AtBA[3*index+ii4][3*ii2+2] = v[ii4][2];
					}
					for(ii4=0; ii4<3; ii4++){
						AtBA[3*ii2  ][3*index+ii4] = v[ii4][0];
						AtBA[3*ii2+1][3*index+ii4] = v[ii4][1];
						AtBA[3*ii2+2][3*index+ii4] = v[ii4][2];
					}
				}else{  /* ii2 == index */
					AtBA[3*index  ][3*index  ] = v[0][0];
					AtBA[3*index+1][3*index+1] = v[1][1];
					AtBA[3*index+2][3*index+2] = v[2][2];

					AtBA[3*index  ][3*index+1] = AtBA[3*index+1][3*index  ]
					                           = (v[0][1] + v[1][0])/2.0;
					AtBA[3*index  ][3*index+2] = AtBA[3*index+2][3*index  ]
					                           = (v[0][2] + v[2][0])/2.0;
					AtBA[3*index+1][3*index+2] = AtBA[3*index+2][3*index+1]
					                           = (v[1][2] + v[2][1])/2.0;
				}
			}
		}
	}

	RC_TRY( free1I(row_size, &min_col_B) );
	RC_TRY( free2D33(MAX_VECT_NUM, row_size, &BA) );
	RC_TRY( free2D33(MAX_VECT_NUM, row_size, &At_array) );
	RC_TRY( free2I(MAX_VECT_NUM, row_size, &At_col_array) );
	RC_TRY( free1D33(row_size, &At_array2) );
	RC_TRY( free1I(row_size, &At_col_array2) );

	return(NORMAL_RC);
}
