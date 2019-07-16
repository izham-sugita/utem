#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fem.h"
#include "base.h"
#include "mathematics.h"
#include "file_io.h"


int main(int argc, char **argv);

int main(int argc, char **argv)
{
	int ii1;
	FILE *fp, *fp_log, *fpw_bdf, *fpw_mgf, *fpw_inp, *fpw_dxf;
	NODE_ARRAY          node;
	ELEMENT_ARRAY       element;
	MATERIAL_PROP_ARRAY material;
	PHYSICAL_PROP_ARRAY physical;
	BC_ARRAY            rest, force, shape_rest, temp;
	BC_SET_ARRAY bc_set;
	DEFAULT_BC_ARRAY b_force, def_temp;
	DISP_ARRAY disp;
	STRESS_ARRAY stress;
	STRAIN_ARRAY strain;
	NST_PARAM_ARRAY bdf_param;
	char buf[128];

	if(argc < 3){
		fprintf(stderr, "Usage : %s [r:model_file(.bdf)]"
		                "[w:result_file_name(to bdr, dxf, inp)]\n", argv[0]);
		return(EXIT_FAILURE);
	}

	/* ログ出力レベル */
	RC_TRY_MAIN( set_log_level(5) );

	/* ログ出力先 */
	RC_TRY_MAIN( rc_fopen("analysis.log", "w", &fp_log) );
	RC_TRY_MAIN( set_log_file(0, stderr) );
	RC_TRY_MAIN( set_log_file(1, fp_log) );

	/* 使用メモリの上限値設定 */
	RC_TRY_MAIN( mm_init((size_t)256*MEGA_BYTE) );

	/* 解析モデルの読み込み */
	RC_TRY_MAIN( rc_fopen(argv[1], "r", &fp) );
	RC_TRY_MAIN( input_bdf(fp, &node, &element, &material, &physical,
	                       &bdf_param, &bc_set, &rest, &force,
	                       &b_force, &shape_rest, &temp, &def_temp));
	RC_TRY_MAIN( rc_fclose(fp) );

	/*
	 * 結果ファイル書き出しのための前処理
	 * 長時間計算の後に解析結果を出力しようとしてファイル処理でエラーに
	 * なると無駄になるため，予防的に解析前にファイルを開いておく．
	 */ 
	fprintf(stderr, "Result model file Open ...");
#ifdef WIN32
	RC_NEG_CHK( sprintf(&buf[0],  "%s_result.bdf", argv[2]));
	RC_TRY_MAIN( rc_fopen(buf, "w", &fpw_bdf) );
	RC_NEG_CHK( sprintf(&buf[0], "%s_result.mgf", argv[2]));
	RC_TRY_MAIN( rc_fopen(buf, "w", &fpw_mgf) );
	RC_NEG_CHK( sprintf(&buf[0], "%s_result.inp", argv[2]));
	RC_TRY_MAIN( rc_fopen(buf, "w", &fpw_inp) );
	RC_NEG_CHK( sprintf(&buf[0], "%s_result.dxf", argv[2]));
	RC_TRY_MAIN( rc_fopen(buf, "w", &fpw_dxf) );
#else
	RC_NEG_CHK( snprintf(&buf[0], sizeof(buf), "%s_result.bdf", argv[2]));
	RC_TRY_MAIN( rc_fopen(buf, "w", &fpw_bdf) );
	RC_NEG_CHK( snprintf(&buf[0], sizeof(buf), "%s_result.mgf", argv[2]));
	RC_TRY_MAIN( rc_fopen(buf, "w", &fpw_mgf) );
	RC_NEG_CHK( snprintf(&buf[0], sizeof(buf), "%s_result.inp", argv[2]));
	RC_TRY_MAIN( rc_fopen(buf, "w", &fpw_inp) );
	RC_NEG_CHK( snprintf(&buf[0], sizeof(buf), "%s_result.dxf", argv[2]));
	RC_TRY_MAIN( rc_fopen(buf, "w", &fpw_dxf) );
#endif
	fprintf(stdout, "OK.\n");

	fprintf(stderr, "FEM Solve ...");
#if 1
	/* スカイライン版ソルバーで計算 */
	RC_TRY_MAIN( fem_solve(node, element, material, physical,
	                       rest, force, temp, def_temp, 
	                       &disp, &stress, &strain, SOL_SS_NODE) );
#else
	/* ICCG版ソルバーで計算 */
	RC_TRY_MAIN( fem_solve_iccg3(node, element, material, physical,
	                             rest, force, b_force, temp, def_temp, 
	                             &disp, &stress, &strain, SOL_SS_NODE) );
#endif
	fprintf(stdout, "OK.\n");

	/* 変位が大き過ぎるまたは小さ過ぎる場合の処理 */
	for(ii1=0; ii1<node.size; ii1++){
		int index1, index2;

		RC_NEG_CHK(index1 = search_node_label(node, node.array[ii1].label));
		RC_NEG_CHK(index2 = search_disp_node(disp, node.array[ii1].label) );

		node.array[index1].p.x += disp.array[index2].v.x * 10.0;
		node.array[index1].p.y += disp.array[index2].v.y * 10.0;
		node.array[index1].p.z += disp.array[index2].v.z * 10.0;
	}

	/* 解析結果の出力 (BDF, MGF, INP, DXF) */
	fprintf(stderr, "Result model output ");
	RC_TRY_MAIN( output_bdf(fpw_bdf, node, element, material, physical,
                            bdf_param, bc_set, rest, force, b_force,
	                        shape_rest) );
	RC_TRY_MAIN( rc_fclose(fpw_bdf) );

	fprintf(stderr, ".");
	RC_TRY_MAIN( output_avs_fe_model(fpw_mgf, node, element) );
	RC_TRY_MAIN( rc_fclose(fpw_mgf) );

	fprintf(stderr, ".");
	RC_TRY_MAIN( output_avs_stress(fpw_inp, node, element, stress) );
	RC_TRY_MAIN( rc_fclose(fpw_inp) );

	fprintf(stderr, ".");
	RC_TRY_MAIN( output_dxf(fpw_dxf, node, element) );
	RC_TRY_MAIN( rc_fclose(fpw_dxf) );
	fprintf(stderr, " OK !!\n");


	RC_TRY_MAIN( mm_terminate() );

	RC_TRY_MAIN( rc_fclose(fp_log) );

	return(EXIT_SUCCESS);
}

