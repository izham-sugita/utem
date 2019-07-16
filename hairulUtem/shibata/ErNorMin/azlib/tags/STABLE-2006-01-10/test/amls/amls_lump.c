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
static RC make_KM_matrix4full_mat(NODE_ARRAY node, ELEMENT_ARRAY element,
                        MATERIAL_PROP_ARRAY material,
                        PHYSICAL_PROP_ARRAY physical,
                        double **K, double **M);
static RC result4full_chk (int block_num, int row_size, int *block_size,
                           int *evects_size_array, double **K,
                           double **M, MATRIX_DB_INFO *mdbi_U);
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
static RC make_UtKU_UtMU_matrix_db(int *block_size_array, int block_num,
                                   int row_size, MATRIX_DB_INFO *mdbi_K,
                                   MATRIX_DB_INFO *mdbi_M,
                                   MATRIX_DB_INFO *mdbi_U);
static RC lumped_matrix33(double matrix[][3]);
static RC initial_U_matrix_db(int row_size, MATRIX_DB_INFO *mdbi_U);
static RC chk_array33(int row_size, int current_pos, int next_pos,
                      int elim_pos, int next_elim_pos, MATRIX_DB_INFO *mdbi_K,
                      int *chk_flag);
static RC matrix_dbi2skyline3(int row_size, int row0_index, int block_size,
                              MATRIX_DB_INFO *mdbi, SKYLINE_MATRIX *skyline);
static RC cal_block_array_size (int *block_size_array, int block_num,
                                int row_size, MATRIX_DB_INFO *mdbi_K);

int main(int argc, char **argv);

int main(int argc, char **argv)
{
	FILE *fp;
	NODE_ARRAY node;
	ELEMENT_ARRAY element;
	MATERIAL_PROP_ARRAY material;
	PHYSICAL_PROP_ARRAY physical;
	int block_num = 10; /* row sizeの分割数(ブロック数) */
	int cache_M = 128;
	int cache_K = 128;
	int cache_U = 128;
	WRAP64_INT cache_M_size;
	WRAP64_INT cache_K_size;
	WRAP64_INT cache_U_size;
	MATRIX_DB_INFO mdbi_K;
	MATRIX_DB_INFO mdbi_M;
	MATRIX_DB_INFO mdbi_U;


	if(argc < 3){
		fprintf(stderr, "Usage : %s [r:model(.bdf)] [scratch dir]\n", argv[0]);
		return(EXIT_FAILURE);
	}

	/* MM, Log */
	RC_TRY_MAIN( set_log_level(5) );
	RC_TRY_MAIN( mm_init((size_t)600*MEGA_BYTE) );
	RC_TRY_MAIN( set_log_file(0, stderr) );
	
	/* グローバル変数初期化 */
	init_matrix33(INIT33);

	/* Model Input */
	RC_TRY_MAIN( rc_fopen(argv[1], "r", &fp) );
	RC_TRY_MAIN( nst_input_node(fp, &node) );
	RC_TRY_MAIN( nst_input_element(fp, &element) );
	RC_TRY_MAIN( nst_input_material_prop(fp, &material) );
	RC_TRY_MAIN( nst_input_physical_prop(fp, &physical) );
	RC_TRY_MAIN( rc_fclose(fp) );
#if 0
	RC_TRY_MAIN( cut_element(&node, &element) );
	RC_TRY_MAIN( cut_element(&node, &element) );
	RC_TRY_MAIN( cut_element(&node, &element) );
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
DEBUG_PRINT;
	RC_TRY( cal_block_array_size(block_size_array, block_num, row_size,
	                             mdbi_K) );
DEBUG_PRINT;

	for(ii1=0; ii1<block_num; ii1++){
		log_printf(5, "block_size[%d] = %d\n", ii1, block_size_array[ii1]);
	}
#if 0
	RC_TRY( print_btree_matrix_db(stderr, mdbi_M) );
	RC_TRY( print_btree_matrix_db(stderr, mdbi_K) );
	RC_TRY( before_chk4full_mat(node, element, material, physical) );
	RC_TRY( before_chk(row_size, mdbi_K, mdbi_M) );
#endif
	RC_TRY( make_UtKU_UtMU_matrix_db(block_size_array, block_num, row_size,
	                                 mdbi_K, mdbi_M, mdbi_U) );

//	RC_TRY( print_btree_matrix_db(stderr, mdbi_M) );
//	exit(0);
	RC_TRY( allocate1I(block_num, &evects_size_array) );
	/* 縮退時の固有ベクトル数 */
	for(ii1=0; ii1<block_num; ii1++){
		evects_size_array[ii1] = block_size_array[ii1] * 3;
	}
#if 1
{
	double **K;
	double **M;
	RC_TRY( allocate2D(row_size*3, row_size*3, &K) );
	RC_TRY( allocate2D(row_size*3, row_size*3, &M) );
	RC_TRY( make_KM_matrix4full_mat(node, element, material, physical, K, M) );
	RC_TRY( result4full_chk(block_num, row_size, block_size_array,
	                        evects_size_array,
	                        K, M, mdbi_U) );
	RC_TRY( free2D(row_size*3, row_size*3, &K) );
	RC_TRY( free2D(row_size*3, row_size*3, &M) );
}
#endif

	RC_TRY( allocate1I(block_num, &evects_size_array) );
	/* 縮退時の固有ベクトル数 */
	for(ii1=0; ii1<block_num; ii1++){
		evects_size_array[ii1] = block_size_array[ii1] * 3;
	}
#if 1
	RC_TRY( result_chk(block_num, block_size_array, row_size, evects_size_array,
                       mdbi_K, mdbi_M, mdbi_U) );
#endif

	RC_TRY( free1I(block_num, &block_size_array) );
	RC_TRY( free1I(block_num, &evects_size_array) );

	return(NORMAL_RC);

}


static RC
make_KM_matrix4full_mat (NODE_ARRAY node, ELEMENT_ARRAY element,
                         MATERIAL_PROP_ARRAY material,
                         PHYSICAL_PROP_ARRAY physical,
                         double **K, double **M)
{
	int ii1, ii2, ii3;

	for(ii1=0; ii1<element.size; ii1++){
		ELEM_MATRIX elem_matrix;
		if(element.array[ii1].label < 0) continue;

		RC_TRY( make_elem_matrix(element.array[ii1], node,
		                         material, physical, &elem_matrix) );
		for(ii2=0; ii2<elem_matrix.size; ii2++){
			for(ii3=0; ii3<elem_matrix.size; ii3++){
				K[elem_matrix.index[ii2]][elem_matrix.index[ii3]]
				+= elem_matrix.matrix[ii2][ii3];
			}
		}
		RC_TRY( free_elem_matrix(&elem_matrix) );

		RC_TRY( make_mass_matrix(element.array[ii1], node,
		                         material, physical, &elem_matrix) );
		for(ii2=0; ii2<elem_matrix.size; ii2++){
			for(ii3=0; ii3<elem_matrix.size; ii3++){
				M[elem_matrix.index[ii2]][elem_matrix.index[ii3]]
				+= elem_matrix.matrix[ii2][ii3];
			}
		}
		RC_TRY( free_elem_matrix(&elem_matrix) );
	}

	return(NORMAL_RC);
}


static RC
result4full_chk (int block_num, int row_size, int *block_size,
                 int *evects_size_array,
                 double **K, double **M, MATRIX_DB_INFO *mdbi_U)
{
	int ii1, ii2, ii3, ii4;
	int total_size = row_size * 3;
	double *evals;
	double **evects;
	double **U;
	double **Z;
	double **UtMU;
	double **UtKU;
	double (*array)[3][3];
	int *col_array;
	int col_size;
	int vect_size;
	double **block_K;
	double **block_M;
	double **sub_mat;
	int row_pos;
	int col_pos;


	RC_TRY( allocate2D(total_size, total_size, &U) );
	RC_TRY( allocate1D33(row_size, &array) );
	RC_TRY( allocate1I(row_size, &col_array) );
	for(ii1=0; ii1<row_size; ii1++){
		/* mdbi_U -> U */
		RC_TRY( get_row_matrix_db(mdbi_U, ii1, &col_size,
		                          array, col_array) );
		for(ii2=0; ii2<col_size; ii2++){
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					U[ii1*3+ii3][col_array[ii2]*3+ii4] = array[ii2][ii3][ii4];
				}
			}
		}
	}
#if 0
	double **UtM;
	RC_TRY( allocate2D(total_size, total_size, &UtM) );
	mul_matrix_AtB(total_size, total_size, total_size, U, M, UtM);
	for(ii1=block_size[0]*3; ii1<total_size; ii1++){
		for(ii2=0; ii2<block_size[0]*3; ii2++){
			RC_TRY( log_printf(5, "UtM[%3d][%3d] = % e\n",
			ii1-block_size[0]*3,ii2, UtM[ii1][ii2]) );
		}
	}
	exit(0);
#endif
#if 0
	double **KU;
	RC_TRY( allocate2D(total_size, total_size, &KU) );
	mul_matrix(total_size, total_size, total_size, K, U, KU);
	for(ii1=block_size[0]*3; ii1<total_size; ii1++){
		for(ii2=block_size[0]*3; ii2<total_size; ii2++){
			RC_TRY( log_printf(5, "KU[%3d][%3d] = % e\n",
			ii1-block_size[0]*3,ii2-block_size[0]*3,KU[ii1][ii2]) );
		}
	}
	exit(0);
#endif
#if 0
	double **MU;
	RC_TRY( allocate2D(total_size, total_size, &MU) );
	mul_matrix(total_size, total_size, total_size, M, U, MU);
	for(ii1=block_size[0]*3; ii1<total_size; ii1++){
		for(ii2=block_size[0]*3; ii2<total_size; ii2++){
			RC_TRY( log_printf(5, "MU[%3d][%3d] = % e\n",
			ii1-block_size[0]*3,ii2-block_size[0]*3,MU[ii1][ii2]) );
		}
	}
	exit(0);
#endif
#if 0
	for(ii1=0; ii1<total_size; ii1++){
		double val = 0.0;
		for(ii2=0; ii2<total_size; ii2++){
			val += M[ii1][ii2];
			M[ii1][ii2] = 0.0;
		}
		M[ii1][ii1] = val;
	}
#endif

#if 0
	RC_TRY( allocate2D(total_size, total_size, &UtKU) );
	mul_matrix_AtBA(total_size, total_size, U, K, UtKU); 
	RC_TRY( allocate2D(total_size, total_size, &UtMU) );
	mul_matrix_AtBA(total_size, total_size, U, M, UtMU); 
	RC_TRY( allocate1D(total_size, &evals) );
	RC_TRY( allocate2D(total_size, total_size, &evects) );
	RC_TRY( cal_eval_evect(total_size, UtKU, UtMU, evals, evects) );
	for(ii1=0; ii1<30; ii1++){
		RC_TRY( log_printf(5, "evals[%d] = % e\n",ii1,evals[ii1]) );
	}
	RC_TRY( free1D(total_size, &evals) );
	RC_TRY( free2D(total_size, total_size, &evects) );
	RC_TRY( free2D(total_size, total_size, &UtMU) );
	RC_TRY( free2D(total_size, total_size, &UtKU) );
#endif
#if 1
	RC_TRY( allocate2D(total_size, total_size, &UtKU) );
	mul_matrix_AtBA(total_size, total_size, U, K, UtKU); 
	RC_TRY( allocate2D(total_size, total_size, &UtMU) );
	mul_matrix_AtBA(total_size, total_size, U, M, UtMU); 
	vect_size = 0;
	for(ii1=0; ii1<block_num; ii1++) vect_size += evects_size_array[ii1];
	row_pos = 0;
	col_pos = 0;
	RC_TRY( allocate2D(total_size, vect_size, &Z) );
	for(ii1=0; ii1<block_num; ii1++){
		int b_size = block_size[ii1] * 3;

		RC_TRY( allocate1D(b_size, &evals) );
		RC_TRY( allocate2D(b_size, b_size, &evects) );
		RC_TRY( allocate2D(b_size, b_size, &block_K) );
		RC_TRY( allocate2D(b_size, b_size, &block_M) );
		for(ii2=0; ii2<b_size; ii2++){
			for(ii3=0; ii3<b_size; ii3++){
				block_K[ii2][ii3] = UtKU[row_pos+ii2][row_pos+ii3];
				block_M[ii2][ii3] = UtMU[row_pos+ii2][row_pos+ii3];
			}
		}

		RC_TRY( cal_eval_evect(b_size, block_K, block_M, evals, evects) );

		for(ii2=0; ii2<b_size; ii2++){
			for(ii3=0; ii3<evects_size_array[ii1]; ii3++){
				Z[row_pos+ii2][col_pos+ii3] = evects[ii2][ii3];
			}
		}
		RC_TRY( free1D(b_size, &evals) );
		RC_TRY( free2D(b_size, b_size, &evects) );
		RC_TRY( free2D(b_size, b_size, &block_K) );
		RC_TRY( free2D(b_size, b_size, &block_M) );

		row_pos += block_size[ii1] * 3;
		col_pos += evects_size_array[ii1];
	}

	/* ZtUtKUZ */
	RC_TRY( allocate2D(total_size, total_size, &sub_mat) );
	for(ii1=0; ii1<total_size; ii1++){
		for(ii2=0; ii2<total_size; ii2++){
			sub_mat[ii1][ii2] = UtKU[ii1][ii2];
		}
	}
	mul_matrix_AtBA(total_size, vect_size, Z, sub_mat, UtKU);

	/* ZtUtMUZ */
	for(ii1=0; ii1<total_size; ii1++){
		for(ii2=0; ii2<total_size; ii2++){
			sub_mat[ii1][ii2] = UtMU[ii1][ii2];
		}
	}
	mul_matrix_AtBA(total_size, vect_size, Z, sub_mat, UtMU);
	RC_TRY( free2D(total_size, total_size, &sub_mat) );


	/* 固有値のみ確認のためfree */
	RC_TRY( free2D(total_size, total_size, &U) );
	RC_TRY( free2D(total_size, vect_size, &Z) );


	RC_TRY( allocate1D(vect_size, &evals) );
	RC_TRY( allocate2D(vect_size, vect_size, &evects) );
	RC_TRY( cal_eval_evect(vect_size, UtKU, UtMU, evals, evects) );
	for(ii1=0; ii1<vect_size; ii1++){
		RC_TRY( log_printf(5, "evals[%d] = % e\n", ii1, evals[ii1]) );
	}
	exit(0);

#endif
#if 0
	for(ii1=0; ii1<total_size; ii1++){
		double val = 0.0;
		for(ii2=0; ii2<total_size; ii2++){
			val += UtMU[ii1][ii2];
			UtMU[ii1][ii2] = 0.0;
		}
		UtMU[ii1][ii1] = val;
		if(val < 0.0){
			log_printf(5, "UtMU[%d][%d] = %e\n", ii1, ii1, UtMU[ii1][ii1]);
		}
	}
#endif
#if 0
	RC_TRY( allocate1D33(row_size, &array) );
	RC_TRY( allocate1D(total_size, &evals) );
	RC_TRY( allocate2D(total_size, total_size, &evects) );
	RC_TRY( eigen_jacobi_trans(total_size, UtKU, evals, evects) );
	for(ii1=0; ii1<30; ii1++){
		RC_TRY( log_printf(5, "evals[%d] = % e\n",ii1,evals[ii1]) );
	}
	RC_TRY( free1D(total_size, &evals) );
	RC_TRY( free2D(total_size, total_size, &evects) );
#endif

	RC_TRY( free2D(total_size, total_size, &U) );
	RC_TRY( free1D33(row_size, &array) );
	RC_TRY( free1I(row_size, &col_array) );

	return(NORMAL_RC);
}


static RC
result_chk (int block_num, int *block_size_array, int row_size,
            int *evects_size_array,
            MATRIX_DB_INFO *mdbi_K, MATRIX_DB_INFO *mdbi_M,
            MATRIX_DB_INFO *mdbi_U)
{
	int ii1, ii2, ii3, ii4;
	int total_size = row_size * 3;
	int vect_size;
	double **K;
	double **M;
	double **U;
	double **Z;
	double **sub_mat;
	double *evals;
	double **evects;
	double **block_K;
	double **block_M;
	double (*array)[3][3];
	int *col_array;
	int col_size;
	int row_pos;
	int col_pos;
	
	RC_TRY( allocate2D(total_size, total_size, &K) );
	RC_TRY( allocate2D(total_size, total_size, &M) );
	RC_TRY( allocate2D(total_size, total_size, &U) );

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
		/* mdbi_U -> U */
		RC_TRY( get_row_matrix_db(mdbi_U, ii1, &col_size, array, col_array) );
		for(ii2=0; ii2<col_size; ii2++){
			if(ii1 >= col_array[ii2]){
				for(ii3=0; ii3<3; ii3++){
					for(ii4=0; ii4<3; ii4++){
						U[ii1*3+ii3][col_array[ii2]*3+ii4]
						     = array[ii2][ii3][ii4];
					}
				}
			}
		}
	}
	RC_TRY( free1I(row_size, &col_array) );
	RC_TRY( free1D33(row_size, &array) );

	for(ii1=0; ii1<total_size; ii1++){
		if(M[ii1][ii1] < 0.0){
			log_printf(5, "m[%d][%d] = %e\n", ii1, ii1, M[ii1][ii1]);
		}
	}

#if 0
	for(ii1=block_size_array[0]*3; ii1<total_size; ii1++){
		for(ii2=block_size_array[0]*3; ii2<total_size; ii2++){
			RC_TRY( log_printf(5, "MU[%3d][%3d] = % e\n",
			ii1-block_size_array[0]*3,ii2-block_size_array[0]*3, M[ii1][ii2]) );
		}
	}
	exit(0);
#endif
#if 0
	for(ii1=block_size_array[0]*3; ii1<total_size; ii1++){
		for(ii2=block_size_array[0]*3; ii2<total_size; ii2++){
			RC_TRY( log_printf(5, "KU[%3d][%3d] = % e\n",
			ii1-block_size_array[0]*3,ii2-block_size_array[0]*3, K[ii1][ii2]) );
		}
	}
	exit(0);
#endif
#if 0
	for(ii1=block_size_array[0]*3; ii1<total_size; ii1++){
		for(ii2=0; ii2<block_size_array[0]*3; ii2++){
			RC_TRY( log_printf(5, "UtM[%3d][%3d] = % e\n",
			ii1-block_size_array[0]*3,ii2, M[ii1][ii2]) );
		}
	}
	exit(0);
#endif
#if 0
	/* UtMUデバッグ */
	RC_TRY( allocate1D(total_size, &evals) );
	RC_TRY( allocate2D(total_size, total_size, &evects) );
	RC_TRY( eigen_jacobi_trans(total_size, M, evals, evects) );
	for(ii1=0; ii1<30; ii1++){
		RC_TRY( log_printf(5, "evals[%d] = % e\n",ii1,evals[ii1]) );
	}
	exit(0);
#endif
#if 1
	/* UtKUデバッグ */
	RC_TRY( allocate1D(total_size, &evals) );
	RC_TRY( allocate2D(total_size, total_size, &evects) );
	RC_TRY( eigen_jacobi_trans(total_size, K, evals, evects) );
	for(ii1=0; ii1<30; ii1++){
		RC_TRY( log_printf(5, "evals[%d] = % e\n",ii1,evals[ii1]) );
	}
	exit(0);
#endif


	/* Make Z */
	vect_size = 0;
	for(ii1=0; ii1<block_num; ii1++) vect_size += evects_size_array[ii1];
	row_pos = 0;
	col_pos = 0;
	RC_TRY( allocate2D(total_size, vect_size, &Z) );
	for(ii1=0; ii1<block_num; ii1++){
		int block_size = block_size_array[ii1] * 3;

		RC_TRY( allocate1D(block_size, &evals) );
		RC_TRY( allocate2D(block_size, block_size, &evects) );
		RC_TRY( allocate2D(block_size, block_size, &block_K) );
		RC_TRY( allocate2D(block_size, block_size, &block_M) );
		for(ii2=0; ii2<block_size; ii2++){
			for(ii3=0; ii3<block_size; ii3++){
				block_K[ii2][ii3] = K[row_pos+ii2][row_pos+ii3];
				block_M[ii2][ii3] = M[row_pos+ii2][row_pos+ii3];
			}
		}

		RC_TRY( cal_eval_evect(block_size, block_K, block_M, evals, evects) );

		for(ii2=0; ii2<block_size; ii2++){
			for(ii3=0; ii3<evects_size_array[ii1]; ii3++){
				Z[row_pos+ii2][col_pos+ii3] = evects[ii2][ii3];
			}
		}
		RC_TRY( free1D(block_size, &evals) );
		RC_TRY( free2D(block_size, block_size, &evects) );
		RC_TRY( free2D(block_size, block_size, &block_K) );
		RC_TRY( free2D(block_size, block_size, &block_M) );

		row_pos += block_size_array[ii1] * 3;
		col_pos += evects_size_array[ii1];
	}

	/* ZtUtKUZ */
	RC_TRY( allocate2D(total_size, total_size, &sub_mat) );
	for(ii1=0; ii1<total_size; ii1++){
		for(ii2=0; ii2<total_size; ii2++){
			sub_mat[ii1][ii2] = K[ii1][ii2];
		}
	}
	mul_matrix_AtBA(total_size, vect_size, Z, sub_mat, K);

	/* ZtUtMUZ */
	for(ii1=0; ii1<total_size; ii1++){
		for(ii2=0; ii2<total_size; ii2++){
			sub_mat[ii1][ii2] = M[ii1][ii2];
		}
	}
	mul_matrix_AtBA(total_size, vect_size, Z, sub_mat, M);
	RC_TRY( free2D(total_size, total_size, &sub_mat) );


	/* 固有値のみ確認のためfree */
	RC_TRY( free2D(total_size, total_size, &U) );
	RC_TRY( free2D(total_size, vect_size, &Z) );


	RC_TRY( allocate1D(vect_size, &evals) );
	RC_TRY( allocate2D(vect_size, vect_size, &evects) );
	RC_TRY( cal_eval_evect(vect_size, K, M, evals, evects) );
	for(ii1=0; ii1<vect_size; ii1++){
		RC_TRY( log_printf(5, "evals[%d] = % e\n", ii1, evals[ii1]) );
	}

	RC_TRY( free1D(vect_size, &evals) );
	RC_TRY( free2D(vect_size, vect_size, &evects) );
	RC_TRY( free2D(total_size, total_size, &K) );
	RC_TRY( free2D(total_size, total_size, &M) );

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
make_UtKU_UtMU_matrix_db (int *block_size_array, int block_num, int row_size,
                          MATRIX_DB_INFO *mdbi_K, MATRIX_DB_INFO *mdbi_M,
                          MATRIX_DB_INFO *mdbi_U)
{
	int ii1, ii2, ii3, ii4, ii5, ii6, ii7, ii8;
	double **d_block_inv;
	double **block_U;
	SKYLINE_MATRIX skyline;
	int col_size;
	int current_pos = 0;
	int next_pos;
	int *col_array;
	double (*array)[3][3];


	RC_TRY( allocate1D33(row_size, &array) );
	RC_TRY( allocate1I(row_size, &col_array) );

	/* Uの初期化(対角マトリクスへ) */
	RC_TRY( initial_U_matrix_db(row_size, mdbi_U) );
	/*
	print_btree_matrix_db(stderr, mdbi_U);
	*/
	for(ii1=1; ii1<block_num; ii1++){
		log_printf(5, "[%d] ii1 =%d\n", __LINE__, ii1);
		int elim_pos = 0;
		int next_elim_pos = 0;
		current_pos += block_size_array[ii1-1];
		next_pos = current_pos + block_size_array[ii1];
		log_printf(6, "[%d] current_pos = %d\n", __LINE__, current_pos);
		log_printf(6, "[%d] next_pos = %d\n", __LINE__, next_pos);

		for(ii2=0; ii2<ii1; ii2++){
			log_printf(6, "[%d] ii2 =%d\n", __LINE__, ii2);
			int chk_flag;
			double (*temp_array)[3][3];
			double (*temp_array2)[3][3];

			next_elim_pos += block_size_array[ii2];

			/* chk_flag = 1だと消去範囲に非ゼロ有り */
			/* chk_flag = 0だと消去範囲ゼロのみ */
			chk_flag = 0;
			RC_TRY( chk_array33(row_size, current_pos, next_pos, elim_pos,
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

			/* Kii^-1 */
			RC_TRY( allocate2D(block_size_array[ii2]*3,
			                   block_size_array[ii2]*3, &d_block_inv) );
			RC_TRY( matrix_dbi2skyline3(row_size, elim_pos, 
			                            block_size_array[ii2],
			                            mdbi_K, &skyline) );
			RC_TRY( decomp_g_matrix(skyline) );
			/* Initialize diagonal matrix */
			for(ii3=0; ii3<block_size_array[ii2]*3; ii3++){
				for(ii4=0; ii4<block_size_array[ii2]*3; ii4++){
					d_block_inv[ii3][ii4] = 0.0;
				}
				d_block_inv[ii3][ii3] = 1.0;
			}
			for(ii3=0; ii3<block_size_array[ii2]*3; ii3++){
				RC_TRY( s_cholesky_subst(skyline.dof, skyline.index1,
				                         skyline.index2, skyline.array,
				                         d_block_inv[ii3], 1) );
			}
			RC_TRY( free_g_matrix(&skyline) );

			/* block_U */
			RC_TRY( allocate2D(block_size_array[ii2]*3,
			                   block_size_array[ii1]*3, &block_U) );
			for(ii3=0; ii3<block_size_array[ii1]; ii3++){
				RC_TRY( get_row_col_matrix_db(mdbi_K, current_pos+ii3, elim_pos,
				                              &col_size, array, col_array) );
				for(ii4=0; ii4<block_size_array[ii2]; ii4++){
					for(ii5=0; ii5<col_size; ii5++){
						if(col_array[ii5] < elim_pos) continue;
						if(col_array[ii5] >= next_elim_pos) break;
						for(ii6=0; ii6<3; ii6++){
							for(ii7=0; ii7<3; ii7++){
								for(ii8=0; ii8<3; ii8++){
									block_U[ii4*3+ii7][ii3*3+ii6]
									  -= d_block_inv[ii4*3+ii7]
									     [(col_array[ii5]-elim_pos)*3+ii8]
									   * array[ii5][ii6][ii8];
								}
							}
						}
					}
				}
			}

			/* ------------------ KU ------------------ */
			RC_TRY( allocate1D33(next_pos-current_pos, &temp_array) );
			for(ii3=current_pos; ii3<row_size; ii3++){
				RC_TRY( get_row_matrix_db(mdbi_K, ii3, &col_size,
				                          array, col_array) );
				for(ii4=current_pos; ii4<next_pos; ii4++){
					for(ii5=0; ii5<col_size; ii5++){
						if(col_array[ii5] < elim_pos) continue;
						if(col_array[ii5] >= next_elim_pos) break;
						for(ii6=0; ii6<3; ii6++){
							for(ii7=0; ii7<3; ii7++){
								for(ii8=0; ii8<3; ii8++){
									temp_array[ii4-current_pos][ii6][ii7]
									+= array[ii5][ii6][ii8]
									* block_U[(col_array[ii5]-elim_pos)*3+ii8]
									         [(ii4-current_pos)*3+ii7];
								}
							}
						}
					}
				}
				/* 書き込む部分(arrayの後半)のデータをtemp_arrayに */
				/* マージしてからput */
				for(ii4=0; ii4<col_size; ii4++){
					if(col_array[ii4] < current_pos) continue;
					if(col_array[ii4] >= next_pos) break;
					for(ii5=0; ii5<3; ii5++){
						for(ii6=0; ii6<3; ii6++){
							temp_array[col_array[ii4]-current_pos][ii5][ii6]
							+= array[ii4][ii5][ii6];
						}
					}
				}
				for(ii4=0; ii4<next_pos-current_pos; ii4++){
					if(ii3 >= ii4+current_pos){
						if( nonzero_chk33(temp_array[ii4]) != 0 ){
							RC_TRY( put_v_matrix_db(ii3, ii4+current_pos,
							                        temp_array[ii4], mdbi_K) );
						}
					}
					init_matrix33(temp_array[ii4]);
				}
			}
			RC_TRY( free1D33(next_pos-current_pos, &temp_array) );

			/* ------------------ UtK ------------------ */
			/* UtKは計算してもdeleteしか起きない(?) */
			for(ii3=current_pos; ii3<next_pos; ii3++){
				for(ii4=elim_pos; ii4<next_elim_pos; ii4++){
					RC_TRY( delete_v_matrix_db(ii3, ii4, mdbi_K) );
				}
			}
			/* ------------------ U = U * block_U  ----------------- */
			/* Multiplication of U (Renewal) */
			RC_TRY( allocate1D33(next_pos-current_pos, &temp_array) );
			for(ii3=0; ii3<elim_pos; ii3++){
				RC_TRY( get_row_matrix_db(mdbi_U, ii3, &col_size,
				                          array, col_array) );
				for(ii4=current_pos; ii4<next_pos; ii4++){
					for(ii5=0; ii5<col_size; ii5++){
						if(col_array[ii5] < elim_pos) continue;
						if(col_array[ii5] >= next_elim_pos) break;
						for(ii6=0; ii6<3; ii6++){
							for(ii7=0; ii7<3; ii7++){
								for(ii8=0; ii8<3; ii8++){
									temp_array[ii4-current_pos][ii6][ii7]
									+= array[ii5][ii6][ii8]
									 * block_U[(col_array[ii5]-elim_pos)*3+ii8]
									          [(ii4-current_pos)*3+ii7];
								}
							}
						}
					}
				}
				/* 書き込む部分のデータをtemp_arrayに */
				/* マージしてからput */
				for(ii4=0; ii4<col_size; ii4++){
					if(col_array[ii4] < current_pos) continue;
					if(col_array[ii4] >= next_pos) break;
					for(ii5=0; ii5<3; ii5++){
						for(ii6=0; ii6<3; ii6++){
							temp_array[col_array[ii4]-current_pos][ii5][ii6]
							+= array[ii4][ii5][ii6];
						}
					}
				}
				for(ii4=0; ii4<next_pos-current_pos; ii4++){
					if( nonzero_chk33(temp_array[ii4]) != 0 ){
						RC_TRY( put_v_matrix_db(ii3, ii4+current_pos,
						                        temp_array[ii4], mdbi_U) );
					}
					init_matrix33(temp_array[ii4]);
				}
			}
			/* Multiplication of U (Fill-in) */
			for(ii3=elim_pos; ii3<next_elim_pos; ii3++){
				for(ii4=0; ii4<(next_pos-current_pos); ii4++){
					for(ii5=0; ii5<3; ii5++){
						for(ii6=0; ii6<3; ii6++){
							temp_array[ii4][ii5][ii6]
							    = block_U[(ii3-elim_pos)*3+ii5]
							             [ii4*3+ii6];

						}
					}
				}
				/* Put (fill-in) */
				/* 書き込む部分がゼロ成分のみ */
				for(ii4=0; ii4<(next_pos-current_pos); ii4++){
					if( nonzero_chk33(temp_array[ii4]) != 0 ){
						RC_TRY(put_v_matrix_db(ii3, ii4+current_pos,
						                       temp_array[ii4], mdbi_U));
					}
				}
			}
			RC_TRY( free1D33(next_pos-current_pos, &temp_array) );

			/* ------------------ UtMU ------------------ */
			/* --------フィルイン時対角項へ足し込み------ */
			RC_TRY( allocate1D33(next_pos-current_pos, &temp_array) );
			RC_TRY( allocate1D33(next_pos-current_pos, &temp_array2) );
			for(ii3=elim_pos; ii3<next_elim_pos; ii3++){
				double v[3][3];
				int found_flag = 0;

				RC_TRY( get_v_matrix_db(ii3, ii3, v, &found_flag, mdbi_M) );
if(ii3==elim_pos)log_printf(5, "ii3 = %d : v[0][0] = %e\n", ii3,v[0][0]);
if(ii3==elim_pos)log_printf(5, "ii3 = %d : v[1][1] = %e\n", ii3,v[1][1]);
if(ii3==elim_pos)log_printf(5, "ii3 = %d : v[2][2] = %e\n", ii3,v[2][2]);
				if(!found_flag) return(UNKNOWN_ERROR_RC);
				for(ii4=current_pos; ii4<next_pos; ii4++){
					for(ii6=0; ii6<3; ii6++){
						for(ii7=0; ii7<3; ii7++){
							for(ii8=0; ii8<3; ii8++){
								temp_array[ii4-current_pos][ii6][ii7]
								 += v[ii6][ii8]
								 * block_U[(ii3-elim_pos)*3+ii8]
								          [(ii4-current_pos)*3+ii7];
							}
						}
					}
				}
				for(ii4=current_pos; ii4<next_pos; ii4++){
					double tmp[3][3];

					RC_TRY( add_matrix33(v, temp_array[ii4-current_pos], v) );

					RC_TRY( get_v_matrix_db(ii4,ii4,tmp,&found_flag,mdbi_M) );
#if 0
if(ii3==elim_pos&&ii4==current_pos)
log_printf(5, "ii3 = %d : tmp[0][0] = %e\n", ii3,tmp[0][0]);
if(ii3==elim_pos&&ii4==current_pos)
log_printf(5, "ii3 = %d : tmp[1][1] = %e\n", ii3,tmp[1][1]);
if(ii3==elim_pos&&ii4==current_pos)
log_printf(5, "ii3 = %d : tmp[2][2] = %e\n", ii3,tmp[2][2]);
#endif
					if(found_flag){
						for(ii5=0; ii5<3; ii5++){
							for(ii6=0; ii6<3; ii6++){
								tmp[ii5][ii6] +=
								temp_array[ii4-current_pos][ii6][ii5];
							}
						}
						RC_TRY( lumped_matrix33(tmp) );
#if 0
if(tmp[0][0] < 0.0)log_printf(5, "ii3 = %d : tmp[0][0] = %e\n", ii3,tmp[0][0]);
if(tmp[1][1] < 0.0)log_printf(5, "ii3 = %d : tmp[1][1] = %e\n", ii3,tmp[1][1]);
if(tmp[2][2] < 0.0)log_printf(5, "ii3 = %d : tmp[2][2] = %e\n", ii3,tmp[2][2]);
#endif
						RC_TRY( put_v_matrix_db(ii4, ii4, tmp, mdbi_M) );
					}else{
						return(UNKNOWN_ERROR_RC);
					}
				}
				RC_TRY( lumped_matrix33(v) );
if(v[0][0] < 0.0)log_printf(5, "v[0][0] = %e\n", v[0][0]);
if(v[1][1] < 0.0)log_printf(5, "v[1][1] = %e\n", v[1][1]);
if(v[2][2] < 0.0)log_printf(5, "v[2][2] = %e\n", v[2][2]);
if(ii3 == elim_pos)log_printf(5, "v[0][0] = %e\n", v[0][0]);
if(ii3 == elim_pos)log_printf(5, "v[1][1] = %e\n", v[1][1]);
if(ii3 == elim_pos)log_printf(5, "v[2][2] = %e\n", v[2][2]);
				RC_TRY( put_v_matrix_db(ii3, ii3, v, mdbi_M) );
#if 1
				for(ii4=current_pos; ii4<next_pos; ii4++){
					double tmp[3][3];

					for(ii5=current_pos; ii5<next_pos; ii5++){
						for(ii6=0; ii6<3; ii6++){
							for(ii7=0; ii7<3; ii7++){
								for(ii8=0; ii8<3; ii8++){
									temp_array2[ii5-current_pos][ii7][ii8]
									+= block_U[(ii3-elim_pos)*3+ii6]
									          [(ii4-current_pos)*3+ii7]
									* temp_array[ii5-current_pos][ii6][ii8];
								}
							}
						}
					}
					RC_TRY( get_v_matrix_db(ii4,ii4,tmp,&found_flag,mdbi_M) );
					if(!found_flag) return(UNKNOWN_ERROR_RC);
					for(ii5=0; ii5<(next_pos-current_pos); ii5++){
						RC_TRY(add_matrix33(tmp, temp_array2[ii5], tmp));
					}
					RC_TRY( lumped_matrix33(tmp) );
#if 0
if(tmp[0][0] < 0.0)log_printf(5, "ii3 = %d : tmp[0][0] = %e\n", ii3,tmp[0][0]);
if(tmp[1][1] < 0.0)log_printf(5, "ii3 = %d : tmp[1][1] = %e\n", ii3,tmp[1][1]);
if(tmp[2][2] < 0.0)log_printf(5, "ii3 = %d : tmp[2][2] = %e\n", ii3,tmp[2][2]);
#endif
					RC_TRY( put_v_matrix_db(ii4, ii4, tmp, mdbi_M) );
					for(ii5=0; ii5<(next_pos-current_pos); ii5++){
						init_matrix33(temp_array2[ii5]);
					}
				}
#endif

				for(ii4=0; ii4<(next_pos-current_pos); ii4++){
					init_matrix33(temp_array[ii4]);
				}
			}
			RC_TRY( free1D33(next_pos-current_pos, &temp_array) );
			RC_TRY( free1D33(next_pos-current_pos, &temp_array2) );
		
			elim_pos += block_size_array[ii2];
			RC_TRY( free2D(block_size_array[ii2]*3,
			               block_size_array[ii1]*3, &block_U) );
			RC_TRY( free2D(block_size_array[ii2]*3,
			               block_size_array[ii2]*3, &d_block_inv) );
		}
	}


	RC_TRY( free1I(row_size, &col_array) );
	RC_TRY( free1D33(row_size, &array) );

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
chk_array33 (int row_size, int current_pos, int next_pos, int elim_pos,
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
		double tmp[3][3];

		if(element.array[ii1].label < 0) continue;
		
		RC_TRY( make_elem_matrix(element.array[ii1], node,
		                         material, physical, &elem_matrix) );
		for(ii2=0; ii2<elem_matrix.size; ii2+=3){
			for(ii3=0; ii3<elem_matrix.size; ii3+=3){
				int found_flag;
				int ix2 = elem_matrix.index[ii2] / 3;
				int ix3 = elem_matrix.index[ii3] / 3;

				if(ix2 < ix3) continue;

				for(ii4=0; ii4<3; ii4++){
					for(ii5=0; ii5<3; ii5++){
						v[ii4][ii5] = elem_matrix.matrix[ii2+ii4][ii3+ii5];
					}
				}
				RC_TRY( get_v_matrix_db(ix2, ix3, tmp, &found_flag, mdbi_K) );
				if(found_flag){
					for(ii4=0; ii4<3; ii4++){
						for(ii5=0; ii5<3; ii5++){
							v[ii4][ii5] += tmp[ii4][ii5];
						}
					}
					RC_TRY( put_v_matrix_db(ix2, ix3, v, mdbi_K) );
				}else{
					RC_TRY( put_v_matrix_db(ix2, ix3, v, mdbi_K) );
				}
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

		for(ii2=0; ii2<elem_matrix.size; ii2+=3){
			int found_flag;
			int ix2 = elem_matrix.index[ii2] / 3;

			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					v[ii3][ii4] = elem_matrix.matrix[ii2+ii3][ii2+ii4];
				}
			}
			RC_TRY( get_v_matrix_db(ix2, ix2, tmp, &found_flag, mdbi_M) );
			if(found_flag){
				for(ii3=0; ii3<3; ii3++){
					for(ii4=0; ii4<3; ii4++){
						v[ii3][ii4] += tmp[ii3][ii4];
					}
				}
				RC_TRY( put_v_matrix_db(ix2, ix2, v, mdbi_M) );
			}else{
				RC_TRY( put_v_matrix_db(ix2, ix2, v, mdbi_M) );
			}
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

#if 0
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

