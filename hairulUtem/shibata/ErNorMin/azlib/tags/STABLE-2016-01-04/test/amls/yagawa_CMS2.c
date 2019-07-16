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
static RC mul_AtBA_matrix_db_old(int row_size, int At_row_size,
                                 MATRIX_DB_INFO *mdbi_At,
                                 MATRIX_DB_INFO *mdbi_B,
                                 MATRIX_DB_INFO *mdbi_AtBA);
static RC mul_matrix_vect_db_sp33(int row_size, MATRIX_DB_INFO *A,
                                  const int min_col_A[], int vect_num,
                                  const int vect_size[], int **vect_row,
                                  double (**vect)[3][3], double (**ans)[3][3]);
static RC mul_AtBA_matrix_db(int row_size, int At_row_size,
                              MATRIX_DB_INFO *mdbi_At, MATRIX_DB_INFO *mdbi_B,
                              double **AtBA);
static RC inner_product33sp(int size1, double array1[][3][3],
                const int index1[], int size2, double array2[][3][3],
                const int index2[], double v[3][3], int *set_flag);

static double INIT33[3][3];

static RC make_matrix4amls(NODE_ARRAY node, ELEMENT_ARRAY element,
                           MATERIAL_PROP_ARRAY material,
                           PHYSICAL_PROP_ARRAY physical,
                           int block_num,
                           MATRIX_DB_INFO *mdbi_K, MATRIX_DB_INFO *mdbi_M);
RC result_chk_cms (int block_num, int block_size_array[], int row_size,
                   int evects_size_array[], MATRIX_DB_INFO *mdbi_K,
                   MATRIX_DB_INFO *mdbi_M);
int main(int argc, char **argv);

int main(int argc, char **argv)
{
	FILE *fp;
	NODE_ARRAY node;
	ELEMENT_ARRAY element;
	MATERIAL_PROP_ARRAY material;
	PHYSICAL_PROP_ARRAY physical;
	int block_num = 7; /* row sizeの分割数(ブロック数) */
	int cache_M = 64;
	int cache_K = 64;
	WRAP64_INT cache_M_size;
	WRAP64_INT cache_K_size;
	MATRIX_DB_INFO mdbi_K;
	MATRIX_DB_INFO mdbi_M;

	RIGID_ELEMENT_ARRAY rigid;

	if(argc < 3){
		fprintf(stderr, "Usage : %s [r:model(.bdf)] [scratch dir]\n", argv[0]);
		return(EXIT_FAILURE);
	}

	/* MM, Log */
	RC_TRY_MAIN( set_log_level(5) );
	RC_TRY_MAIN( mm_init((size_t)700*MEGA_BYTE) );
	RC_TRY_MAIN( set_log_file(0, stderr) );
	
	/* グローバル変数初期化 */
	init_matrix33(INIT33);

	/* Model Input */
	RC_TRY_MAIN( rc_fopen(argv[1], "r", &fp) );
	RC_TRY_MAIN( nst_input_node(fp, &node) );
	RC_TRY_MAIN( nst_input_element(fp, &element) );
	RC_TRY_MAIN( nst_input_material_prop(fp, &material) );
	RC_TRY_MAIN( nst_input_physical_prop(fp, &physical) );
	RC_TRY_MAIN( nst_input_rigid_element(fp, &rigid) )
	RC_TRY_MAIN( rc_fclose(fp) );
#if 0
	RC_TRY_MAIN( cut_element(&node, &element) );
	RC_TRY_MAIN( cut_element(&node, &element) );
	RC_TRY_MAIN( cut_element(&node, &element) );
	RC_TRY_MAIN( fem_renumber0(&node, element, rigid, 1) );
#endif
	RC_TRY_MAIN( dummy_renumber(&node) );

	fprintf(stderr, "node = %d\n", node.size);

	print_matrix_db_size(stderr);
	cache_M_size = sizeof(MATRIX_DB)*MATRIX_DB_ASSOCIATE*cache_M;
	RC_TRY_MAIN( open_matrix_db(argv[2], cache_M_size, &mdbi_M) );

	cache_K_size = sizeof(MATRIX_DB)*MATRIX_DB_ASSOCIATE*cache_K;
	RC_TRY_MAIN( open_matrix_db(argv[2], cache_K_size, &mdbi_K) );


	RC_TRY_MAIN( make_matrix4amls(node, element, material, physical,
	                              block_num, &mdbi_K, &mdbi_M) );


	RC_TRY_MAIN( close_matrix_db(&mdbi_K) );
	RC_TRY_MAIN( close_matrix_db(&mdbi_M) );

	RC_TRY_MAIN( free_node_array(&node) );
	RC_TRY_MAIN( free_element_array(&element) );
	RC_TRY_MAIN( free_material_prop_array(&material) );
	RC_TRY_MAIN( free_physical_prop_array(&physical) );

	RC_TRY_MAIN( mm_terminate() );

	return(EXIT_SUCCESS);
}


static RC
make_matrix4amls (NODE_ARRAY node, ELEMENT_ARRAY element,
                  MATERIAL_PROP_ARRAY material , PHYSICAL_PROP_ARRAY physical,
                  int block_num, MATRIX_DB_INFO *mdbi_K, MATRIX_DB_INFO *mdbi_M)
{
	int ii1;
	int row_size;
	int *block_size_array;
	int *evects_size_array;

	
	/* -------------全体剛性，全体質量マトリクス------------------- */
	RC_TRY( make_KM_matrix_db_amls_bg(node, element, material, physical, mdbi_K,
	                                  mdbi_M) );


	/* ブロックサイズ */
	row_size = count_valid_node(node);
	RC_TRY( log_printf(5, "row_size = %d\n", row_size) );
	RC_TRY( allocate1I(block_num, &block_size_array) );
	RC_TRY( cal_block_array_size_amls_bg(block_size_array, block_num, row_size,
	                                     mdbi_K) );


	tlog_printf(4, "[%d]\n", __LINE__);

	RC_TRY( allocate1I(block_num, &evects_size_array) );
	/* 縮退時の固有ベクトル数 */
	for(ii1=0; ii1<block_num; ii1++){
		evects_size_array[ii1] = block_size_array[ii1];
	}
#if 0
	for(ii1=1; ii1<block_num; ii1+=2){
		evects_size_array[ii1] += 1;
	}
#endif
#if 1
	RC_TRY( result_chk_cms(block_num, block_size_array, row_size,
	                       evects_size_array, mdbi_K, mdbi_M) );
#endif

	RC_TRY( free1I(block_num, &block_size_array) );
	RC_TRY( free1I(block_num, &evects_size_array) );

	return(NORMAL_RC);

}


RC
result_chk_cms (int block_num, int block_size_array[], int row_size,
                int evects_size_array[], MATRIX_DB_INFO *mdbi_K,
                MATRIX_DB_INFO *mdbi_M)
{
	int ii1, ii2, ii3, ii4, ii5;
	int vect_size;
	int real_vect_size;
	double *evals;
	double **evects;
	double (*array)[3][3];
	int *col_array;
	int *add_b_size_array;
	int col_size;
	int col_pos3;
	int diag_pos3;
	int next_diag_pos3;
	double **block_K;
	double **block_M;
	double **ZtKZ;
	double **ZtMZ;
	MATRIX_DB_INFO mdbi_Zt;
	MATRIX_DB_INFO mdbi_m;
	MATRIX_DB_INFO mdbi_k;


	RC_TRY( open_matrix_db(".", 128*1024*1024, &mdbi_Zt) );
	RC_TRY( open_matrix_db(".", 32*1024*1024, &mdbi_m) );
	RC_TRY( open_matrix_db(".", 32*1024*1024, &mdbi_k) );

	RC_TRY( allocate1I(block_num, &add_b_size_array) );
	add_b_size_array[0] = 0;
	for(ii1=0; ii1<block_num; ii1++){
		add_b_size_array[ii1] = add_b_size_array[ii1] + block_size_array[ii1];
	}

	/* Make Z */
	vect_size = 0;
	for(ii1=0; ii1<block_num; ii1++){
		vect_size += evects_size_array[ii1] * 3;
	}
	diag_pos3 = 0;
	next_diag_pos3 = 0;
	col_pos3 = 0;
	RC_TRY( allocate1D33(row_size, &array) );
	RC_TRY( allocate1I(row_size, &col_array) );
	/* 対角ブロック */
	for(ii1=0; ii1<block_num; ii1++){
		double v[3][3];
		int found_flag;
		int block_size;
		int block_size3;

		next_diag_pos3 += block_size_array[ii1];

		if((ii1%2) == 0){
			block_size = block_size_array[ii1] * 3;
			block_size3 = block_size_array[ii1];

			RC_TRY( allocate1D(block_size, &evals) );
			RC_TRY( allocate2D(block_size, block_size, &evects) );
			RC_TRY( allocate2D(block_size, block_size, &block_K) );
			RC_TRY( allocate2D(block_size, block_size, &block_M) );

			for(ii2=0; ii2<block_size3; ii2++){
				RC_TRY( get_row_range_matrix_db(mdbi_K,ii2+diag_pos3, diag_pos3,
				                                diag_pos3+block_size3,
				                                &col_size, array, col_array) );
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
						block_K[ii2*3+ii3]
						[(col_array[col_size-1]-diag_pos3)*3+ii4]
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

			RC_TRY(cal_eval_evect(block_size, block_K, block_M, evals, evects));

			/* mdbi_Zt */
			for(ii2=0; ii2<evects_size_array[ii1]; ii2++){
				for(ii3=0; ii3<block_size3; ii3++){
					/* Zt[ii2][ii3] -> v */
					for(ii4=0; ii4<3; ii4++){
						for(ii5=0; ii5<3; ii5++){
							v[ii4][ii5] = evects[ii3*3+ii5][ii2*3+ii4];
						}
					}
					RC_TRY(add_v_matrix_db(ii2+col_pos3, ii3+diag_pos3, v,
					                       &mdbi_Zt));
				}
			}
			RC_TRY( free1D(block_size, &evals) );
			RC_TRY( free2D(block_size, block_size, &evects) );
			RC_TRY( free2D(block_size, block_size, &block_K) );
			RC_TRY( free2D(block_size, block_size, &block_M) );
			col_pos3 += evects_size_array[ii1];
		}
		diag_pos3 += block_size_array[ii1];
	}
	/* 対角ブロック and 隣接ブロック結合 */
#if 1
	diag_pos3 = 0;
	next_diag_pos3 = 0;
	for(ii1=0; ii1<block_num; ii1++){
		double v[3][3];
		int found_flag;
		int block_size;
		int block_size3;

		next_diag_pos3 += block_size_array[ii1];

		if(ii1%2==0 && ii1 > 0){ 
fprintf(stderr, "ii1 = %d\n", ii1);
			block_size3 = next_diag_pos3 - diag_pos3;
			block_size = block_size3 * 3;

			RC_TRY( allocate1D(block_size, &evals) );
			RC_TRY( allocate2D(block_size, block_size, &evects) );
			RC_TRY( allocate2D(block_size, block_size, &block_K) );
			RC_TRY( allocate2D(block_size, block_size, &block_M) );

			for(ii2=0; ii2<block_size3; ii2++){
				RC_TRY( get_row_range_matrix_db(mdbi_K,ii2+diag_pos3, diag_pos3,
				                                diag_pos3+block_size3,
				                                &col_size, array, col_array) );
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
						block_K[ii2*3+ii3]
						[(col_array[col_size-1]-diag_pos3)*3+ii4]
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

			RC_TRY(cal_eval_evect(block_size, block_K, block_M, evals, evects));

//			evects_size_array[ii1-1] += 1;
			/* mdbi_Zt */
			for(ii2=0; ii2<evects_size_array[ii1-1]; ii2++){
				int row0 = block_size_array[ii1-2];
				int row1 = row0 + block_size_array[ii1-1];
				for(ii3=0; ii3<block_size3; ii3++){
//				for(ii3=row0; ii3<row1; ii3++){
					/* Zt[ii2][ii3] -> v */
					for(ii4=0; ii4<3; ii4++){
						for(ii5=0; ii5<3; ii5++){
							v[ii4][ii5] = evects[ii3*3+ii5][ii2*3+ii4];
						}
					}
					RC_TRY(add_v_matrix_db(ii2+col_pos3, ii3+diag_pos3, v,
					                       &mdbi_Zt));
				}
			}
			RC_TRY( free1D(block_size, &evals) );
			RC_TRY( free2D(block_size, block_size, &evects) );
			RC_TRY( free2D(block_size, block_size, &block_K) );
			RC_TRY( free2D(block_size, block_size, &block_M) );
			col_pos3 += evects_size_array[ii1-1];
			diag_pos3 = next_diag_pos3 - block_size_array[ii1];
		}
	}
#endif
	RC_TRY( free1D33(row_size, &array) );
	RC_TRY( free1I(row_size, &col_array) );

	fprintf(stderr, "col_pos = %d\n", col_pos3);




	real_vect_size = 0;
	for(ii1=0; ii1<block_num; ii1++){
		real_vect_size += 3*(evects_size_array[ii1]);
	}
	fprintf(stderr, "real_pos = %d\n", real_vect_size/3);
	fprintf(stderr, "row_size = %d\n", row_size);

	/* ZtKZ, ZtMZ */
	RC_TRY( allocate2D(real_vect_size, real_vect_size, &ZtKZ) );
	RC_TRY( allocate2D(real_vect_size, real_vect_size, &ZtMZ) );
#if 1
	RC_TRY( mul_AtBA_matrix_db(row_size, real_vect_size/3, &mdbi_Zt, mdbi_K,
	                           ZtKZ) );
	RC_TRY( mul_AtBA_matrix_db(row_size, real_vect_size/3, &mdbi_Zt, mdbi_M,
	                           ZtMZ) );
	for(ii1=0; ii1<real_vect_size; ii1++){
		if( nearly_eq(ZtMZ[ii1][ii1], 0.0) || nearly_eq(ZtKZ[ii1][ii1], 0.0) ){
			return(CAL_ERROR_RC);
		}
	}
#if 0
	for(ii1=0; ii1<real_vect_size; ii1++){
		for(ii2=0; ii2<real_vect_size; ii2++){
			fprintf(stderr, "ZtMZ[%d][%d] = %e\n", ii1, ii2, ZtMZ[ii1][ii2]);
		}
	}
#endif
#endif
#if 0
	RC_TRY( mul_AtBA_matrix_db_old(row_size, real_vect_size/3, &mdbi_Zt, mdbi_M,
	                               &mdbi_m) );
	RC_TRY( mul_AtBA_matrix_db_old(row_size, real_vect_size/3, &mdbi_Zt, mdbi_K,
	                               &mdbi_k) );
DEBUG_PRINT;
	RC_TRY( allocate1D33(real_vect_size/3, &array) );
	RC_TRY( allocate1I(real_vect_size/3, &col_array) );
	for(ii1=0; ii1<real_vect_size/3; ii1++){
		RC_TRY( get_row_range_matrix_db(&mdbi_m, ii1, 0, ii1+1, &col_size,
		                                array, col_array) );
		for(ii2=0; ii2<col_size; ii2++){
			if(ii1 >= col_array[ii2]){
				for(ii3=0; ii3<3; ii3++){
					for(ii4=0; ii4<3; ii4++){
						ZtMZ[ii1*3+ii3][col_array[ii2]*3+ii4]
						= array[ii2][ii3][ii4];
						ZtMZ[col_array[ii2]*3+ii4][ii1*3+ii3]
						= array[ii2][ii3][ii4];
					}
				}
			}
		}
		RC_TRY( get_row_range_matrix_db(&mdbi_k, ii1, 0, ii1+1, &col_size,
		                                array, col_array) );
		for(ii2=0; ii2<col_size; ii2++){
			if(ii1 >= col_array[ii2]){
				for(ii3=0; ii3<3; ii3++){
					for(ii4=0; ii4<3; ii4++){
						ZtKZ[ii1*3+ii3][col_array[ii2]*3+ii4]
						= array[ii2][ii3][ii4];
						ZtKZ[col_array[ii2]*3+ii4][ii1*3+ii3]
						= array[ii2][ii3][ii4];
					}
				}
			}
		}
	}
	RC_TRY( free1D33(real_vect_size/3, &array) );
	RC_TRY( free1I(real_vect_size/3, &col_array) );
DEBUG_PRINT;
#endif

	RC_TRY( close_matrix_db(&mdbi_Zt) );

	RC_TRY( allocate1D(real_vect_size, &evals) );
	RC_TRY( allocate2D(real_vect_size, real_vect_size, &evects) );

	tlog_printf(5, "[%d]\n", __LINE__);
	RC_TRY( cal_eval_evect(real_vect_size, ZtKZ, ZtMZ, evals, evects) );
	tlog_printf(5, "[%d]\n", __LINE__);
	for(ii1=0; ii1<10; ii1++){
		RC_TRY( log_printf(5, "evals[%d] = %15.7e %15.7e\n", ii1+1, evals[ii1],
		                   sqrt(fabs(evals[ii1]))/(2.0*PI)) );
	}
	RC_TRY( free1D(real_vect_size, &evals) );
	RC_TRY( free2D(real_vect_size, real_vect_size, &evects) );
	RC_TRY( free2D(vect_size, vect_size, &ZtKZ) );
	RC_TRY( free2D(real_vect_size, real_vect_size, &ZtMZ) );
	RC_TRY( allocate1I(block_num, &add_b_size_array) );

	return(NORMAL_RC);
}


/* mdbi_B is symmetric and low tridiagonal */
/* mdbi_AtBA is symmetric and low tridiagonal */
static RC
mul_AtBA_matrix_db_old (int row_size, int At_row_size, MATRIX_DB_INFO *mdbi_At,
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
