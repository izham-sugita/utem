#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "mathematics.h"
#include "fem.h"
#include "base.h"
#include "amls.h"


static double INIT33[3][3];

static RC make_matrix4amls(NODE_ARRAY node, ELEMENT_ARRAY element,
                           MATERIAL_PROP_ARRAY material,
                           PHYSICAL_PROP_ARRAY physical,
                           int block_num,
                           MATRIX_DB_INFO  *mdbi_K, MATRIX_DB_INFO  *mdbi_M,
                           MATRIX_DB_INFO  *mdbi_U);
int main(int argc, char **argv);

int main(int argc, char **argv)
{
	FILE *fp;
	NODE_ARRAY node;
	ELEMENT_ARRAY element;
	MATERIAL_PROP_ARRAY material;
	PHYSICAL_PROP_ARRAY physical;
	int block_num = 10; /* row sizeの分割数(ブロック数) */
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
	RC_TRY_MAIN( set_log_level(5) );
	RC_TRY_MAIN( mm_init((size_t)800*MEGA_BYTE) );
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
	int *max_col_U;
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

#if 0
	RC_TRY( print_btree_matrix_db(stderr, mdbi_M) );
	RC_TRY( print_btree_matrix_db(stderr, mdbi_K) );
#endif
	RC_TRY( allocate1I(row_size, &max_col_U) );
	RC_TRY( make_UtKU_matrix_db_amls_bg(block_size_array, block_num, row_size,
	                                    mdbi_K, mdbi_U, max_col_U) );

	tlog_printf(4, "[%d]\n", __LINE__);

	RC_TRY( allocate1I(block_num, &evects_size_array) );
	/* 縮退時の固有ベクトル数 */
#if 1
	for(ii1=0; ii1<block_num; ii1++){
		evects_size_array[ii1] = block_size_array[ii1] / 1;
	}
#endif
#if 1
	RC_TRY( result_chk_amls_bg(block_num, block_size_array, row_size,
	                           evects_size_array, mdbi_K, mdbi_M, mdbi_U,
	                           max_col_U) );
#endif

	RC_TRY( free1I(block_num, &block_size_array) );
	RC_TRY( free1I(block_num, &evects_size_array) );

	return(NORMAL_RC);

}
