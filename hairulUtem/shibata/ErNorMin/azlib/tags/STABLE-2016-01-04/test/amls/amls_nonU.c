#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "mathematics.h"
#include "fem.h"
#include "base.h"


static double INIT33[3][3];
/* 結果確認 */
static RC result_chk(int block_num, int *block_size_array, int row_size,
                     int * evects_size_array,
                     MATRIX_DB_INFO *mdbi_K, MATRIX_DB_INFO *mdbi_M,
                     MATRIX_DB_INFO *mdbi_U);
/* 変換後のK，Mによる固有値解析(結果確認用) */
static RC cal_eval_evect(int size, double **subK, double **subM,
                         double *e_vals, double **e_vects);
static RC cholesky_decomp(int n, double **matrix, double **L_matrix);
static void eigen_sort (int n, double *e_vals, double **e_vects);
/* K，Mの作成 */
static RC make_matrix4amls(NODE_ARRAY node, ELEMENT_ARRAY element,
                           MATERIAL_PROP_ARRAY material,
                           PHYSICAL_PROP_ARRAY physical,
                           int block_num,
                           MATRIX_DB_INFO *mdbi_K, MATRIX_DB_INFO *mdbi_M,
                           MATRIX_DB_INFO *mdbi_U);
static RC make_KM_matrix_db(NODE_ARRAY node, ELEMENT_ARRAY element,
                            MATERIAL_PROP_ARRAY material,
                            PHYSICAL_PROP_ARRAY physical,
                            MATRIX_DB_INFO *mdbi_K, MATRIX_DB_INFO *mdbi_M);
/* K，Mの変換 */
static RC make_UtKU_matrix_db(int *block_size_array, int block_num,
                              int row_size, MATRIX_DB_INFO *mdbi_K,
                              MATRIX_DB_INFO *mdbi_U);
static RC initial_U_matrix_db(int row_size, MATRIX_DB_INFO *mdbi_U);
static RC chk_elim_block33(int row_size, int current_pos, int next_pos,
                           int elim_pos, int next_elim_pos,
                           MATRIX_DB_INFO *mdbi_K, int *chk_flag);
static RC matrix_dbi2skyline3(int row_size, int row0_index, int block_size,
                              MATRIX_DB_INFO *mdbi, SKYLINE_MATRIX *skyline);
static RC cal_block_array_size (int *block_size_array, int block_num,
                                int row_size, MATRIX_DB_INFO *mdbi_K);
static RC lumped_matrix33(double matrix[][3]);
static RC make_U_matrix_db (MATRIX_DB_INFO *mdbi_U, int current_pos,
                            int next_pos, int elim_pos, int next_elim_pos,
                            int *col_array, double array[][3][3],
                            double **block_U, int **zero_flag);
static RC cal_UtKU(MATRIX_DB_INFO *mdbi_K, int row_size, int current_pos,
                   int next_pos, int elim_pos, int next_elim_pos,
                   int *col_array, double array[][3][3], double **block_U,
                   int **zero_flag);
static RC block_U_zero_flag(int current_pos, int next_pos, int elim_pos,
                            int next_elim_pos, double **block_U,
                            int **zero_flag);
static RC s_cholesky_subst_sp (long n, const long *index1, const long *index2,
                               const double *array, double *vector, int p_flag);
static RC s_cholesky_subst_sp3 (long n, const long *index1, const long *index2,
                                const double *array, double *vector0,
                                double *vector1, double *vector2, int p_flag);
static RC mul_ZtUt_matrix_db(int row_size, int block_num,
                             int block_size_array[], int evects_size_array[],
                             MATRIX_DB_INFO *mdbi_Zt, MATRIX_DB_INFO *mdbi_U);
static RC mul_AtBA_matrix_db (int row_size, int At_row_size,
                              MATRIX_DB_INFO *mdbi_At, MATRIX_DB_INFO *mdbi_B,
                              MATRIX_DB_INFO *mdbi_AtBA);

int main(int argc, char **argv);

int main(int argc, char **argv)
{
	FILE *fp;
	NODE_ARRAY node;
	ELEMENT_ARRAY element;
	MATERIAL_PROP_ARRAY material;
	PHYSICAL_PROP_ARRAY physical;
	int block_num = 5; /* row sizeの分割数(ブロック数) */
	int cache_M = 64;
	int cache_K = 128;
	int cache_U = 128;
	WRAP64_INT cache_M_size;
	WRAP64_INT cache_K_size;
	WRAP64_INT cache_U_size;
	MATRIX_DB_INFO mdbi_K;
	MATRIX_DB_INFO mdbi_M;
	MATRIX_DB_INFO mdbi_U;

	RIGID_ELEMENT_ARRAY rigid;

	if(argc < 3){
		fprintf(stderr, "Usage : %s [r:model(.bdf)] [scratch dir]\n", argv[0]);
		return(EXIT_FAILURE);
	}

	/* MM, Log */
	RC_TRY_MAIN( set_log_level(4) );
	RC_TRY_MAIN( mm_init((size_t)900*MEGA_BYTE) );
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

	cache_M_size = sizeof(MATRIX_DB)*MATRIX_DB_ASSOCIATE*cache_M;
	RC_TRY_MAIN( open_matrix_db(argv[2], cache_M_size, &mdbi_M) );

	cache_K_size = sizeof(MATRIX_DB)*MATRIX_DB_ASSOCIATE*cache_K;
	RC_TRY_MAIN( open_matrix_db(argv[2], cache_K_size, &mdbi_K) );

	cache_U_size = sizeof(MATRIX_DB)*MATRIX_DB_ASSOCIATE*cache_U;
	RC_TRY_MAIN( open_matrix_db(argv[2], cache_U_size, &mdbi_U) );

	RC_TRY_MAIN( make_matrix4amls(node, element, material, physical,
	                              block_num, &mdbi_K, &mdbi_M, &mdbi_U) );


	RC_TRY_MAIN( close_matrix_db(&mdbi_K) );
	RC_TRY_MAIN( close_matrix_db(&mdbi_M) );
	RC_TRY_MAIN( close_matrix_db(&mdbi_U) );

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
                  int block_num,
                  MATRIX_DB_INFO  *mdbi_K, MATRIX_DB_INFO  *mdbi_M,
                  MATRIX_DB_INFO  *mdbi_U)
{
	int ii1;
	int row_size;
	int *block_size_array;
	int *evects_size_array;

	
	/* -------------全体剛性，全体質量マトリクス------------------- */
	RC_TRY( make_KM_matrix_db(node, element, material, physical, mdbi_K,
	                          mdbi_M) );


	/* ブロックサイズ */
	row_size = count_valid_node(node);
	RC_TRY( log_printf(5, "row_size = %d\n", row_size) );
	RC_TRY( allocate1I(block_num, &block_size_array) );
	RC_TRY( cal_block_array_size(block_size_array, block_num, row_size,
	                             mdbi_K) );

#if 0
	RC_TRY( print_btree_matrix_db(stderr, mdbi_M) );
	RC_TRY( print_btree_matrix_db(stderr, mdbi_K) );
#endif
	RC_TRY( make_UtKU_matrix_db(block_size_array, block_num, row_size,
	                                 mdbi_K, mdbi_U) );

	tlog_printf(5, "[%d]\n", __LINE__);

	RC_TRY( allocate1I(block_num, &evects_size_array) );
	/* 縮退時の固有ベクトル数 */
#if 1
	for(ii1=0; ii1<block_num; ii1++){
		evects_size_array[ii1] = block_size_array[ii1];
	}
#endif
#if 0
	for(ii1=0; ii1<block_num; ii1++){
		evects_size_array[ii1] = 10;
	}
#endif
DEBUG_PRINT;
#if 1
	RC_TRY( result_chk(block_num, block_size_array, row_size, evects_size_array,
                       mdbi_K, mdbi_M, mdbi_U) );
#endif

	RC_TRY( free1I(block_num, &block_size_array) );
	RC_TRY( free1I(block_num, &evects_size_array) );

	return(NORMAL_RC);

}


static RC
result_chk (int block_num, int *block_size_array, int row_size,
            int *evects_size_array,
            MATRIX_DB_INFO *mdbi_K, MATRIX_DB_INFO *mdbi_M,
            MATRIX_DB_INFO *mdbi_U)
{
	int ii1, ii2, ii3, ii4;
	int ii5, ii6, ii7;
	int total_size = row_size * 3;
	int vect_size;
	int vect_size3;
	double **K;
	double **M;
	double **ZtUtMUZ;
	double **Z;
	double *evals;
	double **evects;
	double **block_K;
	double **block_M;
	double (*array)[3][3];
	int *col_array;
	int col_size;
	int row_pos;
	int col_pos;
	int row_pos3;
	int col_pos3;
	int diag_pos3;
	int next_diag_pos3;
	double **ZtKZ;
	MATRIX_DB_INFO mdbi_Zt;
	MATRIX_DB_INFO mdbi_m;
	int cache_Zt = 64;
	int cache_m = 64;
	WRAP64_INT cache_Zt_size;
	WRAP64_INT cache_m_size;


	/* MATRIX_DB_INFO */
	cache_Zt_size = sizeof(MATRIX_DB)*MATRIX_DB_ASSOCIATE*cache_Zt;
	RC_TRY_MAIN( open_matrix_db("./", cache_Zt_size, &mdbi_Zt) );
	cache_m_size = sizeof(MATRIX_DB)*MATRIX_DB_ASSOCIATE*cache_m;
	RC_TRY_MAIN( open_matrix_db("./", cache_m_size, &mdbi_m) );


	
DEBUG_PRINT;
	RC_TRY( allocate2D(total_size, total_size, &K) );
	RC_TRY( allocate2D(total_size, total_size, &M) );

	RC_TRY( allocate1D33(row_size, &array) );
	RC_TRY( allocate1I(row_size, &col_array) );
	for(ii1=0; ii1<row_size; ii1++){
		/* mdbi_K(UtKU) -> K */
		RC_TRY( get_row_matrix_db(mdbi_K, ii1, &col_size, array, col_array) );
		for(ii2=0; ii2<col_size; ii2++){
			if(ii1 >= col_array[ii2]){
				for(ii3=0; ii3<3; ii3++){
					for(ii4=0; ii4<3; ii4++){
						K[ii1*3+ii3][col_array[ii2]*3+ii4]
						    = array[ii2][ii3][ii4];
						K[col_array[ii2]*3+ii4][ii1*3+ii3]
						    = array[ii2][ii3][ii4];
					}
				}
			}
		}
		/* mdbi_M(UtMU) -> M */
		RC_TRY( get_row_matrix_db(mdbi_M, ii1, &col_size, array, col_array) );
		for(ii2=0; ii2<col_size; ii2++){
			if(ii1 >= col_array[ii2]){
				for(ii3=0; ii3<3; ii3++){
					for(ii4=0; ii4<3; ii4++){
						M[ii1*3+ii3][col_array[ii2]*3+ii4]
						    = array[ii2][ii3][ii4];
						M[col_array[ii2]*3+ii4][ii1*3+ii3]
						    = array[ii2][ii3][ii4];
					}
				}
			}
		}
	}
	RC_TRY( free1D33(row_size, &array) );
	RC_TRY( free1I(row_size, &col_array) );


DEBUG_PRINT;
	/* Make Z */
	vect_size = 0;
	vect_size3 = 0;
	for(ii1=0; ii1<block_num; ii1++){
		vect_size += evects_size_array[ii1] * 3;
		vect_size3 += evects_size_array[ii1];
	}
	RC_TRY( allocate2D(vect_size, vect_size, &ZtKZ) );
	diag_pos3 = 0;
	next_diag_pos3 = 0;
	row_pos = 0;
	col_pos = 0;
	row_pos3 = 0;
	col_pos3 = 0;
	RC_TRY( allocate2D(total_size, vect_size, &Z) );
	RC_TRY( allocate1D33(row_size, &array) );
	RC_TRY( allocate1I(row_size, &col_array) );
	for(ii1=0; ii1<block_num; ii1++){
		double v[3][3];
		int found_flag;
		int block_size = block_size_array[ii1] * 3;

		next_diag_pos3 += block_size_array[ii1];

		RC_TRY( allocate1D(block_size, &evals) );
		RC_TRY( allocate2D(block_size, block_size, &evects) );
		RC_TRY( allocate2D(block_size, block_size, &block_K) );
		RC_TRY( allocate2D(block_size, block_size, &block_M) );
		for(ii2=0; ii2<block_size; ii2++){
			for(ii3=0; ii3<block_size; ii3++){
				block_K[ii2][ii3] = K[row_pos+ii2][row_pos+ii3];
			}
		}
#if 1
		if(ii1==0){
#endif
#if 0
		if(ii1<block_num){
#endif
			for(ii2=diag_pos3; ii2<next_diag_pos3; ii2++){
				double v[3][3];

				RC_TRY( get_v_matrix_db(ii2, ii2, v, &found_flag, mdbi_M) );
				if(!found_flag) return(UNKNOWN_ERROR_RC);
				for(ii3=0; ii3<3; ii3++){
					block_M[(ii2-diag_pos3)*3+ii3]
					       [(ii2-diag_pos3)*3+ii3] = v[ii3][ii3];
				}
			}
		}else{
			for(ii2=0; ii2<row_size; ii2++){
				double v[3][3];

				RC_TRY( get_v_matrix_db(ii2, ii2, v, &found_flag, mdbi_M) );
				if(!found_flag) continue;
//				RC_TRY(get_row_matrix_db(mdbi_U,ii2,&col_size,array,col_array));
				RC_TRY(get_row_col_matrix_db(mdbi_U, ii2, diag_pos3,
				                             &col_size, array, col_array));

				for(ii3=0; ii3<col_size; ii3++){
//					if(col_array[ii3] < diag_pos3) continue;
					if(col_array[ii3] >= next_diag_pos3) break;;
					for(ii4=0; ii4<=ii3; ii4++){
//						if(col_array[ii4] < diag_pos3) continue;
//						if(col_array[ii4] >= next_diag_pos3) break;;
						for(ii5=0; ii5<3; ii5++){
							for(ii6=0; ii6<3; ii6++){
								for(ii7=0; ii7<3; ii7++){
									block_M
									[(col_array[ii3]-diag_pos3)*3+ii6]
									[(col_array[ii4]-diag_pos3)*3+ii7]
									+= array[ii3][ii5][ii6] * v[ii5][ii5]
									 * array[ii4][ii5][ii7];
								}
							}
						}
					}
				}
			}
			for(ii2=0; ii2<block_size_array[ii1]; ii2++){
				for(ii3=ii2+1; ii3<block_size_array[ii1]; ii3++){
					for(ii4=0; ii4<3; ii4++){
						for(ii5=0; ii5<3; ii5++){
							block_M[ii2*3+ii4][ii3*3+ii5]
							= block_M[ii3*3+ii5][ii2*3+ii4];
						}
					}
				}
			}
		}

		RC_TRY( cal_eval_evect(block_size, block_K, block_M, evals, evects) );
#if 0
		if(ii1==block_num-1){
			for(ii2=0; ii2<block_size; ii2++){
				fprintf(stderr, "e0[%d] = %e\n", ii2, evals[ii2]);
			}
			exit(0);
		}
#endif


		for(ii2=0; ii2<evects_size_array[ii1]*3; ii2++){
			ZtKZ[col_pos+ii2][col_pos+ii2] = evals[ii2];
		}
		for(ii2=0; ii2<block_size; ii2++){
			for(ii3=0; ii3<evects_size_array[ii1]*3; ii3++){
				Z[row_pos+ii2][col_pos+ii3] = evects[ii2][ii3];
			}
		}
		for(ii2=0; ii2<evects_size_array[ii1]; ii2++){
			for(ii3=0; ii3<block_size_array[ii1]; ii3++){
				/* Zt[ii2][ii3] -> v */
				for(ii4=0; ii4<3; ii4++){
					for(ii5=0; ii5<3; ii5++){
#if 1
						v[ii4][ii5] = evects[ii3*3+ii5][ii2*3+ii4];
#endif
#if 0
						v[ii4][ii5] = evects[ii2*3+ii4][ii3*3+ii5];
#endif
					}
				}
#if 1
				RC_TRY( put_v_matrix_db(ii2+col_pos3,ii3+row_pos3,v,&mdbi_Zt) );
#endif
			}
		}
		RC_TRY( free1D(block_size, &evals) );
		RC_TRY( free2D(block_size, block_size, &evects) );
		RC_TRY( free2D(block_size, block_size, &block_K) );
		RC_TRY( free2D(block_size, block_size, &block_M) );

		diag_pos3 += block_size_array[ii1];
		row_pos += block_size_array[ii1] * 3;
		col_pos += evects_size_array[ii1] * 3;
		row_pos3 += block_size_array[ii1];
		col_pos3 += evects_size_array[ii1];
	}
	RC_TRY( free1D33(row_size, &array) );
	RC_TRY( free1I(row_size, &col_array) );
	RC_TRY( free2D(total_size, total_size, &K) );
DEBUG_PRINT;


	/* ZtUtMUZ */
	RC_TRY( allocate1I(vect_size3, &col_array) );
	RC_TRY( allocate1D33(vect_size3, &array) );
	RC_TRY( allocate2D(vect_size, vect_size, &ZtUtMUZ) );
#if 0
	RC_TRY( mul_AtBA_matrix_db(row_size, vect_size3, &mdbi_Zt, mdbi_K,
	                           &mdbi_m) );
	RC_TRY( print_btree_matrix_db(stderr, &mdbi_m) );
exit(0);
#endif

	RC_TRY( mul_ZtUt_matrix_db(row_size, block_num, block_size_array,
	                           evects_size_array, &mdbi_Zt, mdbi_U) );
	RC_TRY( mul_AtBA_matrix_db(row_size, vect_size3, &mdbi_Zt, mdbi_M,
	                           &mdbi_m) );
#if 0
	RC_TRY( print_btree_matrix_db(stderr, &mdbi_m) );
	exit(0);
#endif
	for(ii1=0; ii1<vect_size3; ii1++){
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


	RC_TRY( free1I(vect_size3, &col_array) );
	RC_TRY( free1D33(vect_size3, &array) );
	RC_TRY( close_matrix_db(&mdbi_Zt) );
	RC_TRY( close_matrix_db(&mdbi_m) );


	/* 固有値のみ確認のためfree */
	RC_TRY( free2D(total_size, vect_size, &Z) );

	RC_TRY( allocate1D(vect_size, &evals) );
	RC_TRY( allocate2D(vect_size, vect_size, &evects) );
	RC_TRY( cal_eval_evect(vect_size, ZtKZ, ZtUtMUZ, evals, evects) );
	for(ii1=0; ii1<10; ii1++){
//		RC_TRY( log_printf(5, "evals[%d] = % e\n", ii1, evals[ii1]) );
		fprintf(stderr, "evals[%d] = % e\n", ii1, evals[ii1]);
	}

	RC_TRY( free1D(vect_size, &evals) );
	RC_TRY( free2D(vect_size, vect_size, &evects) );
	RC_TRY( free2D(total_size, total_size, &M) );
	RC_TRY( free2D(vect_size, vect_size, &ZtKZ) );
	RC_TRY( allocate2D(vect_size, vect_size, &ZtUtMUZ) );

	return(NORMAL_RC);
}


/* ZtUt -> Zt */
static RC
mul_ZtUt_matrix_db (int row_size, int block_num, int block_size_array[],
                    int evects_size_array[], MATRIX_DB_INFO *mdbi_Zt,
                    MATRIX_DB_INFO *mdbi_U)
{
	int ii1, ii2, ii3, ii4, ii5, ii6, ii7, ii8;
	int current_pos;
	int next_pos;
	int Zt_current_pos; 
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

	current_pos = row_size;
	next_pos = row_size;
	Zt_current_pos = Zt_row_size;
	for(ii1=block_num-1; ii1>=1; ii1--){
		int next_elim_pos;
		int elim_pos;

		Zt_current_pos -= evects_size_array[ii1];
		current_pos -= block_size_array[ii1];
		next_elim_pos = current_pos;
		elim_pos = next_elim_pos;
		for(ii2=ii1-1; ii2>=0; ii2--){
			int chk_flag;

			elim_pos -= block_size_array[ii2];
			/* Block U Check */
#if 1
			chk_flag = 0;
			RC_TRY(chk_elim_block33(row_size, elim_pos, next_elim_pos,
			                        current_pos, next_pos, mdbi_U, &chk_flag));
			if(chk_flag == 0){
				next_elim_pos -= block_size_array[ii2];
				continue;
			}
#endif

			RC_TRY( allocate1D33(next_elim_pos-elim_pos, &temp_array) );
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
					RC_TRY( get_row_col_matrix_db(mdbi_U, ii4, current_pos,
					                              &U_col_size, U_array,
					                              U_col_array) );
					/* col_sizeで繰り返し */
					for(ii5=0; ii5<U_col_size; ii5++){
#if 1
						if(U_col_array[ii5] < Zt_col_array[0]) continue;
						if(U_col_array[ii5] > Zt_col_array[Zt_col_size-1]){
							break;
						}
#endif
						if(U_col_array[ii5] >= next_pos) break;
						for(ii6=0; ii6<3; ii6++){
							for(ii7=0; ii7<3; ii7++){
								for(ii8=0; ii8<3; ii8++){
#if 0
									temp_array
									[(U_col_array[ii5]-current_pos)*3+ii8]
									[ii6][ii8]
									+= Zt_array[ii5][ii6][ii8]
									 * U_array[ii5][ii7][ii8];
#endif
									temp_array[(ii4-elim_pos)][ii6][ii7]
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
			RC_TRY( free1D33(next_elim_pos-elim_pos, &temp_array) );

			next_elim_pos -= block_size_array[ii2];
		}
//		next_pos = current_pos;
		next_pos -= block_size_array[ii1];
	}
	RC_TRY( free1D33(row_size, &Zt_array) );
	RC_TRY( free1I(row_size, &Zt_col_array) );
	RC_TRY( free1D33(row_size, &U_array) );
	RC_TRY( free1I(row_size, &U_col_array) );
	RC_TRY( free1D33(row_size, &Zt_vect) );

	return(NORMAL_RC);
}


/* mdbi_B is symmetric and low tridiagonal */
/* mdbi_AtBA is symmetric and low tridiagonal */
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
#if 0
				for(ii4=0; ii4<3; ii4++){
					for(ii5=0; ii5<3; ii5++){
						for(ii6=0; ii6<3; ii6++){
							AtB_array[ii2][ii5][ii4] +=
							B_array[ii3][ii4][ii6]
							* At_vect[B_col_array[ii3]][ii5][ii6];
							AtB_array[B_col_array[ii3]][ii5][ii4] +=
							B_array[ii3][ii6][ii4]
							* At_vect[ii2][ii5][ii6];
						}
					}
				}
#endif
#if 1
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
#endif
			}
#if 1
			/* 対角部分 */
			if(B_col_array[B_col_size-1] < At_col_array[0]
			|| B_col_array[B_col_size-1] > At_col_array[At_col_size-1]){
				continue;
			}
#endif
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
	if(matrix[0][0] < 0.0 || nearly_eq(matrix[0][0], 0.0)) return(CAL_ERROR_RC);
	L_matrix[0][0] = sqrt(matrix[0][0]);

	for(ii1=1; ii1<n; ii1++){
		L_matrix[ii1][0] = matrix[ii1][0] / L_matrix[0][0];
	}
	for(ii1=1; ii1<n; ii1++){
		double root_val = matrix[ii1][ii1];
		for(ii2=0; ii2<ii1; ii2++){
			root_val -= L_matrix[ii1][ii2] * L_matrix[ii1][ii2];
		}

		/* 正定値性，ゼロ割り判定 */
		if( root_val < 0.0 || nearly_eq(L_matrix[0][0], 0.0) ){
			return(CAL_ERROR_RC);
		}

		L_matrix[ii1][ii1] = sqrt(root_val);
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
make_UtKU_matrix_db (int *block_size_array, int block_num, int row_size,
                     MATRIX_DB_INFO *mdbi_K, MATRIX_DB_INFO *mdbi_U)
{
	int ii1, ii2, ii3, ii4, ii5, ii6;
	double **block_U;
	double **temp_vect;
	SKYLINE_MATRIX skyline;
	int col_size;
	int current_pos = 0;
	int next_pos;
	int *col_array;
	double (*array)[3][3];
	int **zero_flag;
	int nonzero_count;


	RC_TRY( allocate1D33(row_size, &array) );
	RC_TRY( allocate1I(row_size, &col_array) );

	/* Uの初期化(対角マトリクスへ) */
	RC_TRY( initial_U_matrix_db(row_size, mdbi_U) );

	for(ii1=1; ii1<block_num; ii1++){
		tlog_printf(5, "[%d] ii1 = %d\n", __LINE__, ii1);
		int elim_pos = 0;
		int next_elim_pos = 0;
		current_pos += block_size_array[ii1-1];
		next_pos = current_pos + block_size_array[ii1];
		log_printf(6, "[%d] current_pos = %d\n", __LINE__, current_pos);
		log_printf(6, "[%d] next_pos = %d\n", __LINE__, next_pos);

		for(ii2=0; ii2<ii1; ii2++){
			log_printf(6, "[%d] ii2 =%d\n", __LINE__, ii2);
			int chk_flag;

			next_elim_pos += block_size_array[ii2];

			/* chk_flag = 1だと消去範囲に非ゼロ有り */
			/* chk_flag = 0だと消去範囲ゼロのみ */
			chk_flag = 0;
			RC_TRY( chk_elim_block33(row_size, current_pos, next_pos, elim_pos,
			                         next_elim_pos, mdbi_K, &chk_flag) );
			if(chk_flag == 0){
				for(ii3=current_pos; ii3<next_pos; ii3++){
					for(ii4=elim_pos; ii4<next_elim_pos; ii4++){
						RC_TRY( delete_v_matrix_db(ii3, ii4, mdbi_K) );
					}
				}
				elim_pos += block_size_array[ii2];
				log_printf(6, "[%d] Skip [ii1:ii2]=[%d : %d]\n", __LINE__,
				           ii1, ii2);
				continue;
			}
			log_printf(6, "[%d] elim_pos = %d\n", __LINE__, elim_pos);


			RC_TRY( matrix_dbi2skyline3(row_size, elim_pos, 
			                            block_size_array[ii2],
			                            mdbi_K, &skyline) );
			RC_TRY( decomp_g_matrix(skyline) );
			RC_TRY( allocate2D(3, block_size_array[ii2]*3, &temp_vect) );
			RC_TRY( allocate2D(block_size_array[ii2]*3,
			                   block_size_array[ii1]*3, &block_U) );
			for(ii3=current_pos; ii3<next_pos; ii3++){
				RC_TRY( get_row_col_matrix_db(mdbi_K, ii3, elim_pos,
				                              &col_size, array, col_array) );
				nonzero_count = 0;
				for(ii4=0; ii4<col_size; ii4++){
					if(col_array[ii4] >= next_elim_pos) break;
					nonzero_count++;
					for(ii5=0; ii5<3; ii5++){
						for(ii6=0; ii6<3; ii6++){
							temp_vect[ii5][(col_array[ii4]-elim_pos)*3+ii6]
							= array[ii4][ii5][ii6];
						}
					}
				}
				if(nonzero_count == 0) continue;

#if 1
				RC_TRY( s_cholesky_subst_sp3(skyline.dof, skyline.index1,
				                             skyline.index2, skyline.array,
				                             temp_vect[0], temp_vect[1],
				                             temp_vect[2], 1) );
#endif
				for(ii4=0; ii4<3; ii4++){
#if 0
					RC_TRY( s_cholesky_subst(skyline.dof, skyline.index1,
				                             skyline.index2, skyline.array,
				                             temp_vect[ii4], 1) );
#endif
#if 0
					RC_TRY( s_cholesky_subst_sp(skyline.dof, skyline.index1,
				                                skyline.index2, skyline.array,
				                                temp_vect[ii4], 1) );
#endif
					for(ii5=0; ii5<(next_elim_pos-elim_pos)*3; ii5++){
						block_U[ii5][(ii3-current_pos)*3+ii4]
						= -temp_vect[ii4][ii5];
						temp_vect[ii4][ii5] = 0.0;
					}
				}
			}
			RC_TRY( free2D(3, block_size_array[ii2]*3, &temp_vect) );
			RC_TRY( free_g_matrix(&skyline) );


			/* Zero check for block_U[][] */
			RC_TRY( allocate2I(next_elim_pos-elim_pos, next_pos-current_pos,
			                   &zero_flag) );
			RC_TRY( block_U_zero_flag(current_pos, next_pos, elim_pos,
			                          next_elim_pos, block_U, zero_flag) );

			/* ------------------ UtKU ------------------ */
			RC_TRY( cal_UtKU(mdbi_K, row_size, current_pos, next_pos,
			                 elim_pos, next_elim_pos, col_array, array,
			                 block_U, zero_flag) );
			/* ------------------ block_U -> U ----------------- */
			RC_TRY( make_U_matrix_db(mdbi_U, current_pos, next_pos, elim_pos,
			                         next_elim_pos, col_array, array,
			                         block_U, zero_flag) );

			RC_TRY( free2I(next_elim_pos-elim_pos, next_pos-current_pos,
			               &zero_flag) );
			RC_TRY( free2D(block_size_array[ii2]*3,
			               block_size_array[ii1]*3, &block_U) );
			elim_pos += block_size_array[ii2];
		}
	}

	RC_TRY( free1I(row_size, &col_array) );
	RC_TRY( free1D33(row_size, &array) );

	return(NORMAL_RC);
}


static RC
s_cholesky_subst_sp (long n, const long *index1, const long *index2,
                     const double *array, double *vector, int p_flag)
{
	long ii1, ii2;
	long nz_begin;
	long nz_end;


	/* Set nonzero begin index nz_begin */
	/* Set nonzero end index nz_end */
	nz_begin = n;
	for(ii1=0L; ii1<n; ii1++){
		if( !nearly_eq(vector[ii1], 0.0) ){
			nz_begin = ii1;
			break;
		}
	}
	nz_end = 0L;
	for(ii1=n-1L; ii1>=nz_begin; ii1--){
		if( !nearly_eq(vector[ii1], 0.0) ){
			nz_end = ii1;
			break;
		}
	}

	if((n < 1L)||((void *)array == NULL)
			||((void *)index1 == NULL)
			||((void *)index2 == NULL)
			||((void *)vector == NULL)) return(ARG_ERROR_RC);

	/* 前進代入 */
	if(p_flag) RC_TRY( log_printf(6, "[B]") );
	for(ii1=0L;ii1<n;ii1++){
		long ii2_start = index1[ii1];
		long ii2_finish = ii1;

		if(ii2_start < nz_begin){
			ii2_start = nz_begin;
		}
		if(ii2_finish > nz_end + 1L){
			ii2_finish = nz_end + 1L;
		}
		for(ii2=ii2_start;ii2<ii2_finish;ii2++){
			vector[ii1] -= vector[ii2] * array[index2[ii1] - ii1 + ii2];
		}   
		if(ii2_start < ii2_finish){
			if(ii1 < nz_begin) nz_begin = ii1;
			if(ii1 > nz_end) nz_end = ii1;
		}
	}

	/* 対角行列 */
	for(ii1=nz_begin;ii1<=nz_end;ii1++){
		if(fabs(array[index2[ii1]]) < ABS_TOL) return(CAL_ERROR_RC);
		vector[ii1] /= array[index2[ii1]];
	}

	/* 後退代入 */
	for(ii1=nz_end;ii1>=0;ii1--){
		for(ii2=index1[ii1];ii2<ii1;ii2++){
			vector[ii2] -= vector[ii1] * array[index2[ii1] - ii1 + ii2];
		}   
	}
	if(p_flag) RC_TRY( log_printf(6, "...finished.\n") );

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


static RC
cal_UtKU (MATRIX_DB_INFO *mdbi_K, int row_size, int current_pos,
          int next_pos, int elim_pos, int next_elim_pos, int *col_array,
          double array[][3][3], double **block_U, int **zero_flag)
{

	int ii1, ii2, ii3, ii4, ii5, ii6;
	int col_size;
	double (*temp_array)[3][3];

	/* ------------------ KU ------------------ */
	RC_TRY( allocate1D33(next_pos-current_pos, &temp_array) );
	for(ii1=current_pos; ii1<row_size; ii1++){
		RC_TRY( get_row_col_matrix_db(mdbi_K, ii1, elim_pos,
		                              &col_size, array, col_array) );
		if(col_size == 0) continue;
		for(ii2=0; ii2<col_size; ii2++){
			if(col_array[ii2] >= next_elim_pos) break;
			for(ii3=current_pos; ii3<next_pos; ii3++){
				if(zero_flag[col_array[ii2]-elim_pos][ii3-current_pos]){
					continue;
				}
				for(ii4=0; ii4<3; ii4++){
					for(ii5=0; ii5<3; ii5++){
						for(ii6=0; ii6<3; ii6++){
							temp_array[ii3-current_pos][ii4][ii6]
							+= array[ii2][ii4][ii5]
							 * block_U[(col_array[ii2]-elim_pos)*3+ii5]
							          [(ii3-current_pos)*3+ii6];
						}
					}
				}
			}
		}
		/* 書き込み */
		for(ii2=0; ii2<next_pos-current_pos; ii2++){
			if(ii1 >= ii2+current_pos){
				if( nonzero_chk33(temp_array[ii2]) != 0 ){
					RC_TRY( add_v_matrix_db(ii1, ii2+current_pos,
					temp_array[ii2], mdbi_K) );
				}
			}
			init_matrix33(temp_array[ii2]);
		}
	}
	RC_TRY( free1D33(next_pos-current_pos, &temp_array) );

	/* ------------------ UtK ------------------ */
	/* UtKは計算してもdeleteしか起きない */
	for(ii1=current_pos; ii1<next_pos; ii1++){
		for(ii2=elim_pos; ii2<next_elim_pos; ii2++){
			RC_TRY( delete_v_matrix_db(ii1, ii2, mdbi_K) );
		}
	}

	return(NORMAL_RC);
}


static RC
make_U_matrix_db (MATRIX_DB_INFO *mdbi_U, int current_pos, int next_pos,
                  int elim_pos, int next_elim_pos, int *col_array,
                  double array[][3][3], double **block_U, int **zero_flag)
{
	int ii1, ii2, ii3, ii4;
	double (*temp_array)[3][3];

	/* Put block_U (Fill-in) */
	RC_TRY( allocate1D33(next_pos-current_pos, &temp_array) );
	for(ii1=elim_pos; ii1<next_elim_pos; ii1++){
		for(ii2=0; ii2<(next_pos-current_pos); ii2++){
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					temp_array[ii2][ii3][ii4]
					    = block_U[(ii1-elim_pos)*3+ii3]
					             [ii2*3+ii4];

				}
			}
		}
		/* Put (fill-in) */
		/* 書き込む部分がゼロ成分のみ */
		for(ii2=0; ii2<(next_pos-current_pos); ii2++){
			if( nonzero_chk33(temp_array[ii2]) != 0 ){
				RC_TRY(put_v_matrix_db(ii1, ii2+current_pos,
				                       temp_array[ii2], mdbi_U));
			}
		}
	}
	RC_TRY( free1D33(next_pos-current_pos, &temp_array) );

	return(NORMAL_RC);
}


static RC block_U_zero_flag(int current_pos, int next_pos, int elim_pos,
                            int next_elim_pos, double **block_U,
                            int **zero_flag)
{
	int ii1, ii2, ii3, ii4;


	for(ii1=0; ii1<(next_elim_pos-elim_pos); ii1++){
		for(ii2=0; ii2<(next_pos-current_pos); ii2++){
			int flag = 1;
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					if(nearly_eq(0.0, block_U[ii1*3+ii3][ii2*3+ii4]) != 1){
						flag = 0;
						break;
					}
				}
				if(flag == 0) break;
			}
			zero_flag[ii1][ii2] = flag;
		}
	}

	return(NORMAL_RC);
}


/* Initialize U */
static RC
initial_U_matrix_db (int row_size, MATRIX_DB_INFO *mdbi_U)
{
	int ii1, ii2;
	double v[3][3];


	for(ii1=0; ii1<3; ii1++){
		for(ii2=0; ii2<3; ii2++){
			v[ii1][ii2] = 0.0;
		}
		v[ii1][ii1] = 1.0;
	}
	for(ii1=0; ii1<row_size; ii1++){
		RC_TRY( put_v_matrix_db(ii1, ii1, v, mdbi_U) );
	}

	return(NORMAL_RC);
}


/* 非ゼロがあったら chk_flag = 1 */
static RC
chk_elim_block33 (int row_size, int current_pos, int next_pos, int elim_pos,
                  int next_elim_pos, MATRIX_DB_INFO *mdbi_K, int *chk_flag)
{
	int ii1, ii2;
	int col_size;
	int *col_array;
	double (*array)[3][3];


	RC_TRY( allocate1I(row_size, &col_array) );
	RC_TRY( allocate1D33(row_size, &array) );

	*chk_flag = 0;
	for(ii1=current_pos; ii1<next_pos; ii1++){
		RC_TRY( get_row_matrix_db(mdbi_K, ii1, &col_size, array, col_array) );
		if(col_array[0] >= next_elim_pos) continue;

		for(ii2=0; ii2<col_size; ii2++){
			if(col_array[ii2] < elim_pos) continue;
			if(col_array[ii2] >= next_elim_pos) break;
			if( nonzero_chk33(array[ii2]) != 0 ){
				*chk_flag = 1;
				break;
			}
		}
		if(*chk_flag == 1) break;
	}
	RC_TRY( free1I(row_size, &col_array) );
	RC_TRY( free1D33(row_size, &array) );

	return(NORMAL_RC);
}


/* 下三角で記憶 */
static RC
make_KM_matrix_db (NODE_ARRAY node, ELEMENT_ARRAY element,
                   MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
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
		/*
		log_printf(5, "row0_index+ii1 = %d+%d\n", row0_index, ii1);
		*/
		RC_TRY( get_row_col_matrix_db(mdbi, row0_index+ii1, row0_index,
		                              &size, array, col_array) );
		/*
		log_printf(5, "size = %d\n", size);
		for(ii2=0; ii2<size; ii2++){
			log_printf(5, "col_array[%d] = %d\n", ii2, col_array[ii2]);
			RC_TRY( print_matrix33(stderr, array[ii2]) );
		}
		*/

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


static RC
cal_block_array_size (int *block_size_array, int block_num, int row_size,
                      MATRIX_DB_INFO *mdbi_K)
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


	if(row_size < 1) return(ARG_ERROR_RC);
	if(row_size < block_num) return(ARG_ERROR_RC);
	if(block_num < 1) return(ARG_ERROR_RC);
	/* 多すぎるブロックサイズは危険 */
	if(block_num/row_size > 0.5) return(ARG_ERROR_RC);
	RC_NULL_CHK( block_size_array );
	RC_NULL_CHK( mdbi_K );

	RC_TRY( allocate1I(row_size, &nonzero_size) );

	RC_TRY( allocate1I(row_size, &col_array) );
	RC_TRY( allocate1D33(row_size, &array) );
	for(ii1=0, sum=0; ii1<row_size; ii1++){
		RC_TRY( get_row_matrix_db(mdbi_K, ii1, &col_size, array, col_array) );
		for(ii2=0; ii2<col_size; ii2++){
			if( nonzero_chk33(array[ii2]) != 0 ){
				nonzero_size[ii1]++;
			}
		}
		sum += nonzero_size[ii1];
	}
	RC_TRY( free1I(row_size, &col_array) );
	RC_TRY( free1D33(row_size, &array) );

#if 0
	for(ii1=0; ii1<row_size; ii1++){
		RC_TRY( log_printf(5, "row = %4d, nonzero_size = %d\n", ii1, 
		                   nonzero_size[ii1]) );
	}
#endif

	nonzero_avg = (int)sum/row_size;
	block_avg = (int)sum/block_num;

#if 0
	RC_TRY( log_printf(5, "sum = %lld\n", sum) );
	RC_TRY( log_printf(5, "nonzero_avg = %d\n", nonzero_avg) );
	RC_TRY( log_printf(5, "block_avg(ideal) = %d\n", block_avg) );
#endif

	for(ii1=0, ii2=0, nonzero_num=0; ii1<row_size; ii1++){
		/*
		int next_size = (ii1+1<row_size ? nonzero_size[ii1+1] : 0);
		if(block_avg > (nonzero_num + next_size)){
		if(block_avg > (nonzero_num + next_size/2)){
		*/
		if(block_avg > (nonzero_num + nonzero_avg)){
			nonzero_num += nonzero_size[ii1];
			block_size_array[ii2]++;
		}else if(ii2+1 == block_num){
			block_size_array[ii2]++;
		}else{
			if(block_size_array[ii2] == 0){
				block_size_array[ii2]++;
				nonzero_num = 0;
			}else{;
				if(ii2+1 >= block_num) break;
				block_size_array[ii2+1]++;
				nonzero_num = nonzero_size[ii1];
			}
			ii2++;
			if(ii2 >= block_num) break;
		}
	}

#if 1
	{
		int ii3, tmp;
		for(ii1=0, ii3=0; ii1<block_num; ii1++){
			int tmp=0;
			RC_TRY( log_printf(5, "block_size_array[%4d] = %3d", ii1, 
			                   block_size_array[ii1]) );
			for(ii2=0; ii2<block_size_array[ii1]; ii2++){
				tmp += nonzero_size[ii3];
				ii3++;
			}
			RC_TRY( log_printf(5, "  size = %d \n", tmp) );
		}

		for(ii1=0, ii3=0; ii1<block_num; ii1++){
			ii3 += block_size_array[ii1];
		}
		RC_TRY( log_printf(5, "[Check] row sum = %3d\n", ii3) );

		for(ii1=0, ii3=0, tmp=0; ii1<block_num; ii1++){
			for(ii2=0; ii2<block_size_array[ii1]; ii2++){
				tmp += nonzero_size[ii3];
				ii3++;
			}
		}
		RC_TRY( log_printf(5, "[Check] nonzero sum = %3d\n", tmp) );
	}
#endif

	RC_TRY( free1I(row_size, &nonzero_size) );

	return(NORMAL_RC);
}


static RC
lumped_matrix33(double matrix[][3])
{
	int ii1, ii2;


	for(ii1=0; ii1<3; ii1++){
		double lump_val = 0.0;

		for(ii2=0; ii2<3; ii2++){
			lump_val += matrix[ii1][ii2];
			matrix[ii1][ii2] = 0.0;
		}
		matrix[ii1][ii1] = lump_val;
	}

	return(NORMAL_RC);
}
