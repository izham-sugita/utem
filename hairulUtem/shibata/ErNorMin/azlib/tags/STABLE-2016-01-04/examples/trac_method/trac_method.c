#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "base.h"
#include "mathematics.h"
#include "fem.h"


#define SAFETY_STR_RATIO  (0.90)    /* SAFETY_STR_RATIO * MAX_STR で概算 */
#define MAX_DIVISOR       (65536.0)
#define DIV_DIVISOR       (1.0)
#define ARMIJO            (0.50)
#define DIA_FACT          (0.01)

/* PARAM構造体 */
typedef struct{
	int it_max;
	int nr_max;
	int sol_flag;
	int opt_flag;
	double max_str;
	double subj_ratio;
	double subj_rel_error;
	double obj_rel_error;
} PARAM;


int main(int argc, char **argv);
RC optimization(NODE_ARRAY node, ELEMENT_ARRAY element,
                MATERIAL_PROP_ARRAY material,
                PHYSICAL_PROP_ARRAY physical, BC_ARRAY rest,
                BC_ARRAY force, BC_ARRAY shape_rest, 
                BC_ARRAY temp, DEFAULT_BC_ARRAY b_force,
                DEFAULT_BC_ARRAY def_temp, PARAM param);
RC cal_compliance(NODE_ARRAY node, ELEMENT_ARRAY element,
                  MATERIAL_PROP_ARRAY material,
                  PHYSICAL_PROP_ARRAY physical, BC_ARRAY rest,
                  BC_ARRAY force, ELEMENT_ARRAY surf,
                  BC_ARRAY shape_rest, BC_ARRAY temp, 
                  DEFAULT_BC_ARRAY b_force, DEFAULT_BC_ARRAY def_temp,
                  double *comp, BC_ARRAY *sens_trac, int sol_flag);
RC cal_volume(NODE_ARRAY node, ELEMENT_ARRAY element,
              ELEMENT_ARRAY surf, BC_ARRAY shape_rest,
              double *volume, BC_ARRAY *sens_trac);
int chk_armijo(double old_obj, double new_obj, double sbj0, double old_sbj,
               double new_sbj, double G_nu[][2],
               double delta_s[], double *divisor, PARAM param);
int subj_compare(double comp, double comp0, PARAM param);
RC input_param(FILE *fp, PARAM *param);
RC input_bdf(FILE *fp, NODE_ARRAY *node, ELEMENT_ARRAY *element,
             MATERIAL_PROP_ARRAY *material,
             PHYSICAL_PROP_ARRAY *physical,
             BC_ARRAY *rest, BC_ARRAY *force,
             DEFAULT_BC_ARRAY *b_force,
             BC_ARRAY *temp, DEFAULT_BC_ARRAY *def_temp,
             RIGID_ELEMENT_ARRAY *rigid);
RC output_bdf(FILE *fpw, NODE_ARRAY node, ELEMENT_ARRAY element,
              MATERIAL_PROP_ARRAY material,
              PHYSICAL_PROP_ARRAY physical, BC_ARRAY rest,
              BC_ARRAY force, DEFAULT_BC_ARRAY b_force);
RC output_dxf(FILE *fpw, NODE_ARRAY node, ELEMENT_ARRAY element);
RC output_stl(FILE *fpw, NODE_ARRAY node, ELEMENT_ARRAY element);
RC output_log(FILE *stream, int iteration, double obj0, double obj, double sbj0,
              double sbj);
RC char_chain(char base_name[], char end_name[], char *result_fname);


int
main (int argc, char **argv)
{
	FILE *fp, *fp_log, *fpw_bdf, *fpw_dxf, *fpw_mgf, *fpw_stl;
	NODE_ARRAY node;
	ELEMENT_ARRAY element;
	MATERIAL_PROP_ARRAY material;
	PHYSICAL_PROP_ARRAY physical;
	BC_ARRAY force, rest, shape_rest, temp;
	DEFAULT_BC_ARRAY b_force, def_temp;
	RIGID_ELEMENT_ARRAY rigid;
	PARAM param;
	char bdf_result_fname[256];
	char dxf_result_fname[256];
	char mgf_result_fname[256];
	char stl_result_fname[256];

	if( argc < 5 ){
		fprintf(stderr, "Usage : %s [r:parameter_file] [r:fem_model(*.bdf)] "
		                "[r:rest_model(*.bdf)] [w:result_file_name"
		                "(to bdf, dxf, mgf and stl)]\n", argv[0]);
		return(EXIT_FAILURE);
	}

	/* ログ出力レベル */
	RC_TRY_MAIN( set_log_level(5) );

	/* ログ出力先 */
	RC_TRY_MAIN( rc_fopen("iteration.log", "w", &fp_log) );
	RC_TRY_MAIN( set_log_file(0, stderr) );
	RC_TRY_MAIN( set_log_file(1, fp_log) );

	/* 使用メモリの上限値設定 */
	RC_TRY_MAIN( mm_init((size_t)256*MEGA_BYTE) );

	/* パラメータ読み込み */
	RC_TRY_MAIN( rc_fopen(argv[1], "r", &fp) );
	RC_TRY_MAIN( input_param(fp, &param) );
	RC_TRY_MAIN( rc_fclose(fp) );

	/* モデルデータ読み込み */
	RC_TRY_MAIN( rc_fopen(argv[2], "r", &fp) );
	RC_TRY_MAIN( input_bdf(fp, &node, &element, &material, &physical,
	                       &rest, &force, &b_force, &temp, &def_temp, &rigid) );
	RC_TRY_MAIN( rc_fclose(fp) );

	/* 拘束条件データ読み込み */
	RC_TRY_MAIN( rc_fopen(argv[3], "r", &fp) );
	RC_TRY_MAIN( nst_input_restraint(fp, &shape_rest) );
	RC_TRY_MAIN( rc_fclose(fp) );
	
	fprintf(stdout, "Number of shape restraints : %d\n", shape_rest.size);

	/* renumber処理 */
	fprintf(stderr, "renumbering...");
	if(param.sol_flag == 0){ /* スカイライン用 */
		RC_TRY_MAIN( fem_renumber0(&node, element, rigid, 1) );
	}else if(param.sol_flag == 1){ /* ICCG用 */
		RC_TRY_MAIN( dummy_renumber(&node) );
	}else{
		fprintf(stderr, "NG !!\n");
		return(1);
	}
	fprintf(stderr, "OK\n");

	/* ファイル書き出しのための前処理 */
	fprintf(stderr, "Trac result model file Open ... ");
	RC_TRY_MAIN( char_chain(argv[4], "_trac_result.bdf", bdf_result_fname) );
	RC_TRY_MAIN( char_chain(argv[4], "_trac_result.dxf", dxf_result_fname) );
	RC_TRY_MAIN( char_chain(argv[4], "_trac_result.mgf", mgf_result_fname) );
	RC_TRY_MAIN( char_chain(argv[4], "_trac_result.stl", stl_result_fname) );
	RC_TRY_MAIN( rc_fopen(bdf_result_fname, "w", &fpw_bdf) );
	RC_TRY_MAIN( rc_fopen(dxf_result_fname, "w", &fpw_dxf) );
	RC_TRY_MAIN( rc_fopen(mgf_result_fname, "w", &fpw_mgf) );
	RC_TRY_MAIN( rc_fopen(stl_result_fname, "w", &fpw_stl) );
	fprintf(stdout, "OK\n\n");

	fprintf(stdout, "### optimization START ! ###\n");
	RC_TRY_MAIN( optimization(node, element, material, physical, rest, force,
	                          shape_rest, temp, b_force, def_temp, param) );
	fprintf(stdout, "### optimization FINISH !! ###\n");

	/* ファイル書き出し処理 */
	fprintf(stdout, "\nTrac result model Output ... ");
	RC_TRY_MAIN( output_bdf(fpw_bdf, node, element, material, physical, rest,
	                        force, b_force) );
	RC_TRY_MAIN( rc_fclose(fpw_bdf) );

	RC_TRY_MAIN( output_dxf(fpw_dxf, node, element) );
	RC_TRY_MAIN( rc_fclose(fpw_dxf) );

	RC_TRY_MAIN( output_avs_fe_model(fpw_mgf, node, element) );
	RC_TRY_MAIN( rc_fclose(fpw_mgf) );

	RC_TRY_MAIN( output_stl(fpw_stl, node, element) );
	RC_TRY_MAIN( rc_fclose(fpw_stl) );
	fprintf(stdout, "OK !!\n\n");

	fprintf(stdout, "++++++++++++ Analysis Finish !! ++++++++++++\n\n");

	RC_TRY_MAIN( mm_terminate() );

	RC_TRY_MAIN( rc_fclose(fp_log) );

	return(EXIT_SUCCESS);
}


RC
optimization (NODE_ARRAY node, ELEMENT_ARRAY element,
              MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
              BC_ARRAY rest, BC_ARRAY force, BC_ARRAY shape_rest, 
              BC_ARRAY temp, DEFAULT_BC_ARRAY b_force,
              DEFAULT_BC_ARRAY def_temp, PARAM param)
{
	int ii1, ii2, ii3;
	int active_flag = 0;
	int over_flag = 0;
	int armijo_result = 0;
	double delta_s[2] = {0.0, 0.0};
	double new_delta_s[2] = {0.0, 0.0};
	double G_nu[2][2] = {{0.0, 0.0},{0.0, 0.0}};
	double sbj0 = 0.0, sbj, sbj_new = 0.0;
	double obj0 = 0.0, obj, obj_new = 0.0;
	double violation;
	double max_str;
	double divisor = 1.0;
	FILE *fp_log;
	NODE_ARRAY node_new;
	ELEMENT_ARRAY surf;
	BC_ARRAY sens_trac[2];
	BC_ARRAY new_sens_trac[2];
	DISP_ARRAY sens_disp[2];
	DISP_ARRAY disp;

	/* extracting surface */
	RC_TRY_MAIN( extract_surface(element, &surf) );

	RC_TRY( allocate_node_array(node.size, &node_new) );
	for(ii1=0; ii1<2; ii1++){
		RC_TRY( allocate_bc_array(0, &(sens_trac[ii1])) );
		RC_TRY( allocate_disp_array(0, &(sens_disp[ii1])) );
		RC_TRY( allocate_bc_array(0, &(new_sens_trac[ii1])) );
	}

	/* 状態方程式を解いて、目的汎関数、制約汎関数を計算 */
	/* 感度に比例した外力(sens_trac[]) も計算           */
	if(param.opt_flag == 1){ /* 剛性最大化 */
		RC_TRY( cal_compliance(node, element, material, physical, rest, force,
		                       surf, shape_rest, temp, b_force, def_temp,
		                       &obj, &(sens_trac[0]), param.sol_flag) );
		RC_TRY( cal_volume(node, element, surf, shape_rest, &sbj,
		                   &(sens_trac[1])) );
	}else{ /* 体積最小化 */
		RC_TRY( cal_volume(node, element, surf, shape_rest, &obj,
		                   &(sens_trac[0])) );
		RC_TRY( cal_compliance(node, element, material, physical, rest, force,
		                       surf, shape_rest, temp, b_force, def_temp,
		                       &sbj, &(sens_trac[1]), param.sol_flag) );
	}
	obj0 = obj;
	sbj0 = sbj;
	if( (obj0 <= ABS_ERROR)||(sbj0 <= ABS_ERROR)){
		fprintf(stderr, "obj0 or sbj0 < \n");
		return(CAL_ERROR_RC);
	}

	/***** 速度場解析 *****/
	RC_TRY( velocity_analysis(node, element, material, physical, shape_rest,
	                          2, sens_trac, sens_disp, param.sol_flag) );

	RC_TRY( rc_fopen("iteration_history.dat", "w", &fp_log) );
	/* 各汎関数の各変位(sens_disp)に対する変動量 G_nu[][] を計算 */
	for(ii1=0; ii1<param.it_max; ii1++){
		fprintf(stderr, "\n");
		RC_TRY( output_log(fp_log, ii1+1, obj0, obj, sbj0, sbj) );
		for(ii2=0; ii2<2; ii2++){
			for(ii3=0; ii3<2; ii3++){
				RC_TRY( functional_variation(sens_trac[ii2], sens_disp[ii3],
				                             &(G_nu[ii2][ii3])) );
				fprintf(stderr, "G_nu[%d][%d]: %15.7e\n", ii2, ii3,
				        G_nu[ii2][ii3]);
			}
			if(G_nu[ii2][ii2] <= ABS_ERROR){
				fprintf(stderr, "G_nu[%d][%d] < \n", ii2, ii2);
				return(CAL_ERROR_RC);
			}
		}

		/* 仮の荷重係数を計算 */
		delta_s[0] = -1.0;
		delta_s[1] = -(delta_s[0]*G_nu[1][0])/G_nu[1][1];

		/* 仮の荷重係数における最大主ひずみ(max_str)を計算 */
		RC_TRY( max_pr_strain_disp(2, delta_s, node, element, sens_disp,
		                           &max_str, SOL_SS_CENTER) );
		if(max_str <= ABS_ERROR){
			fprintf(stderr, "max_str < \n");
			return(CAL_ERROR_RC);
		}

		/* 正式な荷重係数を求める */
		delta_s[0] = -(SAFETY_STR_RATIO * param.max_str)/(divisor * max_str);
		violation = sbj - param.subj_ratio*sbj0;
		delta_s[1] = -(delta_s[0]*G_nu[1][0] + violation)/G_nu[1][1];
		if(delta_s[1] > 0.0){
			delta_s[1] = 0.0;
			active_flag = 0;
		}else{
			active_flag = 1;
		}
		fprintf(stderr, "delta_s[0]:%15.7e  delta_s[1]:%15.7e\n",
		        delta_s[0], delta_s[1]);

		/* 最大ひずみが設定値を越えていれば、再度修正 */
		RC_TRY( max_pr_strain_disp(2, delta_s, node, element, sens_disp,
		                           &max_str, SOL_SS_CENTER) );
		fprintf(stderr, "max_str:%15.7e\n", max_str);
		if(max_str > (param.max_str/divisor)){
			delta_s[0] *= param.max_str / (divisor * max_str);
			delta_s[1] *= param.max_str / (divisor * max_str);
			fprintf(stderr, "delta_s[0]:%15.7e  delta_s[1]:%15.7e\n", delta_s[0],
			        delta_s[1]);
			RC_TRY( max_pr_strain_disp(2, delta_s, node, element, sens_disp,
			                           &max_str, SOL_SS_CENTER) );
			fprintf(stderr, "max_str:%15.7e\n", max_str);
			over_flag = 1;
		}else{
			over_flag = 0;
		}

		/* 変位(sens_disp) の delta_s 倍を節点座標(node_new)に付加 */
		RC_TRY( copy_node_array(node, &node_new) );
		RC_TRY( add_disp_node(node_new, 2, delta_s, sens_disp) );

		RC_TRY( free_bc_array(&(new_sens_trac[0])) );
		RC_TRY( free_bc_array(&(new_sens_trac[1])) );

		/* 状態方程式を解いて、目的汎関数、制約汎関数を計算 */
		if(param.opt_flag == 1){ /* 剛性最大化 */
			RC_TRY( cal_compliance(node_new, element, material, physical, rest,
			                       force, surf, shape_rest, temp, b_force,
			                       def_temp, &obj_new, NULL, param.sol_flag) );
			RC_TRY( cal_volume(node_new, element, surf, shape_rest, &sbj_new,
			                   &new_sens_trac[1]) );
		}else{ /* 体積最小化 */
			RC_TRY( cal_volume(node_new, element, surf, shape_rest,
			                   &obj_new, NULL) );
			RC_TRY( cal_compliance(node_new, element, material, physical, rest,
			                       force, surf, shape_rest, temp, b_force,
			                       def_temp, &sbj_new, &new_sens_trac[1],
			                       param.sol_flag) );
		}

		fprintf(stderr, "%d: new obj =%15.7e\n", ii1+1, obj_new);
		fprintf(stderr, "%d: new sbj =%15.7e\n", ii1+1, sbj_new);
		if((active_flag == 0)&&(sbj_new > param.subj_ratio*sbj0)){
			/* 制約が inactive と予想されたが、実際には active だった */
			active_flag = 1;
			over_flag = 0;
		}

		/***** Armijo チェック *****/
		armijo_result = chk_armijo(obj, obj_new, sbj0, sbj, sbj_new, G_nu,
		                           delta_s, &divisor, param);
		if(armijo_result == 0){
			ii1--;
			continue;
		}else if(armijo_result < 0){
			break;
		}

		new_delta_s[0] = delta_s[0];
		new_delta_s[1] = delta_s[1];
		if((over_flag == 0)&&(active_flag == 1)){
			int subj_flag = 0;
			double d_delta_s1 = 0.0;
			double old_delta_s1;
			double new_G_nu11;
			double sbj_old = sbj_new;
			double armijo;

			fprintf(stderr, "##### in newton #####\n");
			for(ii2=0; ii2<param.nr_max; ii2++){
				RC_TRY( functional_variation(new_sens_trac[1], sens_disp[1],
				                             &(new_G_nu11)) );
				/* 荷重係数、節点座標を修正 */
				/*d_delta_s1 -= (sbj_new - SUBJ_RATIO*sbj0) / G_nu[1][1];*/
				d_delta_s1 -= (sbj_new - param.subj_ratio*sbj0) / new_G_nu11;
				old_delta_s1 = new_delta_s[1];
				new_delta_s[1] = delta_s[1] + d_delta_s1;
				if(new_delta_s[1] > 0.0) new_delta_s[1] = 0.0;
				if( nearly_eq(new_delta_s[1], 0.0)
				 && nearly_eq(old_delta_s1, 0.0) ){
					subj_flag = 1;
					break;    /* 不等式制約を満足 */
				}

				RC_TRY( max_pr_strain_disp(2, new_delta_s, node, element,
				                           sens_disp, &max_str,
				                           SOL_SS_CENTER) );
				if(max_str > param.max_str){
					fprintf(stderr, "max_str:%15.7e\n", max_str);
					break;
				}
				RC_TRY( copy_node_array(node, &node_new) );
				RC_TRY( add_disp_node(node_new, 2, new_delta_s, sens_disp) );

				/* コンプライアンス制約のチェック */
				RC_TRY( free_bc_array(&(new_sens_trac[1])) );
				if(param.opt_flag == 1){
					RC_TRY( cal_volume(node_new, element, surf, shape_rest,
					                   &sbj_new, &(new_sens_trac[1])) );
				}else{
					RC_TRY( cal_compliance(node_new, element, material,
					                       physical, rest, force, surf,
					                       shape_rest, temp, b_force, def_temp,
					                       &sbj_new, &(new_sens_trac[1]),
					                       param.sol_flag) );
				}
				fprintf(stderr, "%d-%d: new sbj =%15.7e\n", ii1+1, ii2+1,
				        sbj_new);

				if( subj_compare(sbj_new, sbj0, param) ){
					subj_flag = 1;
					break;   /* 等式制約を満足 */
				}
				armijo = (sbj_new - sbj_old) / (d_delta_s1*G_nu[1][1]);
				if((armijo < ARMIJO)||(armijo > 1.0/ARMIJO)){
					fprintf(stderr, "armijo: %15.7e\n", armijo);
					break;
				}
			}
			if(subj_flag == 0){
				divisor *= 2.0;
				ii1--;
				fprintf(stderr, "divisor: %f\n", divisor);
				continue;
			}
		}

		RC_TRY( free_bc_array(&(new_sens_trac[0])) );
		if(param.opt_flag == 1){
			RC_TRY( cal_compliance(node_new, element, material, physical, rest,
			                       force, surf, shape_rest, temp, b_force,
			                       def_temp, &obj_new, &(new_sens_trac[0]),
			                       param.sol_flag) );
		}else{
			RC_TRY( cal_volume(node_new, element, surf, shape_rest, &obj_new,
			                   &(new_sens_trac[0])) );
		}
		fprintf(stderr, "%d: new obj =%15.7e\n", ii1+1, obj_new);

		/* Armijo チェック */
		armijo_result = chk_armijo(obj, obj_new, sbj0, sbj, sbj_new, G_nu,
		                           new_delta_s, &divisor, param);
		if(armijo_result == 0){
			ii1--;
			continue;
		}else if(armijo_result < 0){
			break;
		}

		/* 節点座標を更新 */
		RC_TRY( copy_node_array(node_new, &node) );
		for(ii2=0; ii2<2; ii2++){
			RC_TRY( free_bc_array(&(sens_trac[ii2])) );
			sens_trac[ii2] = new_sens_trac[ii2];
			RC_TRY( allocate_bc_array(0, &(new_sens_trac[ii2])) );
		}
		obj = obj_new;
		sbj = sbj_new;

		/* 速度場解析 */
		/* 速度場 sens_disp[], 速度場のひずみ sens_strain[] を更新 */
		RC_TRY( free_disp_array(&(sens_disp[0])) );
		RC_TRY( free_disp_array(&(sens_disp[1])) );
		RC_TRY( velocity_analysis(node, element, material, physical,
		        shape_rest, 2, sens_trac, sens_disp, param.sol_flag) );
		divisor /= DIV_DIVISOR; 
		if(divisor < 1.0) divisor = 1.0;
		fprintf(stderr, "divisor: %f\n", divisor);
	}
	RC_TRY( output_log(fp_log, ii1+1, obj0, obj, sbj0, sbj) );
	RC_TRY( rc_fclose(fp_log) );

	/***** fem analysis of optimum shape *****/
	if(param.sol_flag == 0){
		RC_TRY( fem_solve(node, element, material, physical, rest, force, temp,
		                  def_temp, &disp, NULL, NULL, SOL_SS_CENTER) );
	}else if(param.sol_flag == 1){
		RC_TRY( fem_solve_iccg3(node, element, material, physical, rest,
		                        force, b_force, temp, def_temp, &disp,
		                        NULL, NULL, SOL_SS_CENTER) );
	}else{
		return(IMPLEMENT_ERROR_RC);
	}

	if(armijo_result < 0){
		fprintf(stderr, ": converged!\n");
	}else{
		fprintf(stderr, ": Not converged!\n");
	}

	for(ii2=0; ii2<2; ii2++){
		RC_TRY( free_bc_array(&(sens_trac[ii2])) );
		RC_TRY( free_bc_array(&(new_sens_trac[ii2])) );
		RC_TRY( free_disp_array(&(sens_disp[ii2])) );
	}
	RC_TRY( free_element_array(&surf) );

	return(NORMAL_RC);
}


/* コンプライアンスを計算 */
/* sens_trac != NULL なら感度に比例する荷重も計算 */
RC
cal_compliance (NODE_ARRAY node, ELEMENT_ARRAY element,
                MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
                BC_ARRAY rest, BC_ARRAY force, ELEMENT_ARRAY surf,
                BC_ARRAY shape_rest, BC_ARRAY temp, DEFAULT_BC_ARRAY b_force,
                DEFAULT_BC_ARRAY def_temp, double *comp, BC_ARRAY *sens_trac,
                int sol_flag)
{
	DISP_ARRAY disp;

	fprintf(stderr, "cal_compliance...\n");
	if(sol_flag == 0){
		RC_TRY( fem_solve(node, element, material, physical, rest, force, temp,
                          def_temp, &disp, NULL, NULL, SOL_SS_CENTER) );
	}else if(sol_flag == 1){
		RC_TRY( fem_solve_iccg3(node, element, material, physical, rest,
		                        force, b_force, temp, def_temp, &disp,
		                        NULL, NULL, SOL_SS_CENTER) );
	}else{
		return(IMPLEMENT_ERROR_RC);
	}

	if(comp == NULL) return(ARG_ERROR_RC);
	RC_TRY( compliance_functional(force, disp, comp) );

	if(sens_trac == NULL){
		RC_TRY( free_disp_array(&disp) );
	}else{
		SENSITIVITY_ARRAY sens;
		RC_TRY( shape_sens_compliance_disp(surf, disp, node, element,
		                                   material, physical, &sens) );
		RC_TRY( free_disp_array(&disp) );
		RC_TRY( sens2traction(node, surf, shape_rest, sens, sens_trac, 0) );
		RC_TRY( free_sensitivity_array(&sens) );
	}

	return(NORMAL_RC);
}


/* 体積を計算 */
/* sens_trac != NULL なら感度に比例する荷重も計算 */
RC
cal_volume (NODE_ARRAY node, ELEMENT_ARRAY element, ELEMENT_ARRAY surf,
            BC_ARRAY shape_rest, double *volume, BC_ARRAY *sens_trac)
{
	fprintf(stderr, "cal_volume...\n");

	if(volume == NULL) return(ARG_ERROR_RC);
	RC_TRY( volume_functional(node, element, volume) );

	if(sens_trac != NULL){
		SENSITIVITY_ARRAY sens;
		RC_TRY( shape_sens_volume(surf, &sens) );
		RC_TRY( sens2traction(node, surf, shape_rest, sens, sens_trac, 0) );
		RC_TRY( free_sensitivity_array(&sens) );
	}

	return(NORMAL_RC);
}


/* Armijo の条件をチェック */
/* 1... 合格、0... 失格、-1... エラー、-2...収束 */
int
chk_armijo (double old_obj, double new_obj, double sbj0, double old_sbj,
            double new_sbj, double G_nu[][2], double delta_s[],
            double *divisor, PARAM param)
{
	double d_sbj = delta_s[0] * G_nu[1][0] + delta_s[1] * G_nu[1][1];
	double d_obj = delta_s[0] * G_nu[0][0] + delta_s[1] * G_nu[0][1];

	double lambda;
	double d_L, delta_L;

	lambda = delta_s[1]/delta_s[0];
	/* lambda = -G_nu[1][0]/G_nu[1][1]; */
	fprintf(stderr, "lambda = %15.7e\n", lambda);

	d_L = d_obj + lambda*d_sbj;
	delta_L = (new_obj - old_obj) + lambda*(new_sbj - old_sbj);

	fprintf(stderr, "d_L = %15.7e, delta_L = %15.7e\n", d_L, delta_L);
	if( (-delta_L < (param.obj_rel_error*fabs(new_obj) + ABS_ERROR))
	  &&( (subj_compare(new_sbj, sbj0, param))
	    ||(new_sbj < param.subj_ratio*sbj0) ) ){
		fprintf(stderr, "converged: objective function\n");
		return(-2);
	}

	if(nearly_eq(d_L, 0.0)){
		fprintf(stderr, "converged: d_L\n");
		return(-2);
	}
	fprintf(stderr, "delta_L/d_L = %15.7e\n", delta_L/d_L);
	if((delta_L/d_L < ARMIJO)||(delta_L/d_L > 1.0/ARMIJO)){
		(*divisor) *= 2.0;
		fprintf(stderr, "divisor: %f\n", (*divisor));
		if((*divisor) > MAX_DIVISOR){
			fprintf(stderr, "converged: divisor\n");
			return(-2);
		}
		return(0);
	}

	return(1);
}


/* 制約(sbj)が等式制約(SUBJ_RATIO*sbj0)を */
/* 満足していれば 1 を返却 */
int
subj_compare (double sbj, double sbj0, PARAM param)
{
	double subj_comp = param.subj_ratio * sbj0;

	if( fabs(sbj - subj_comp)
	     < (param.subj_rel_error*fabs(subj_comp) + ABS_ERROR) ){
		return(1);
	}

	return(0);
}


RC
input_param (FILE *fp, PARAM *param)
{
	char buf[1024];
	IDC idc[2];

	while(1){
		if(fgets(buf, sizeof(buf), fp) == NULL){
			if(feof(fp)) break;
			return(READ_ERROR_RC);
		}
		if( !strncmp("IT_MAX", buf, 6) ){
			RC_TRY( sgetidc(buf, "ci", idc) );
			param->it_max = idc[1].i;
			fprintf(stdout, "IT_MAX : %d\n", param->it_max);
		}else if( !strncmp("NR_MAX", buf, 6) ){
			RC_TRY( sgetidc(buf, "ci", idc) );
			param->nr_max = idc[1].i;
			fprintf(stdout, "NR_MAX : %d\n", param->nr_max);
		}else if( !strncmp("SOL_FLAG", buf, 8) ){
			RC_TRY( sgetidc(buf, "ci", idc) );
			param->sol_flag = idc[1].i;
			fprintf(stdout, "SOL_FLAG : %d\n", param->sol_flag);
		}else if( !strncmp("OPT_FLAG", buf, 8) ){
			RC_TRY( sgetidc(buf, "ci", idc) );
			param->opt_flag = idc[1].i;
			fprintf(stdout, "OPT_FLAG : %d\n", param->opt_flag);
		}else if( !strncmp("MAX_STR", buf, 7) ){
			RC_TRY( sgetidc(buf, "cd", idc) );
			param->max_str = idc[1].d;
			fprintf(stdout, "MAX_STR : %15.7e\n", param->max_str);
		}else if( !strncmp("SUBJ_RATIO", buf, 10) ){
			RC_TRY( sgetidc(buf, "cd", idc) );
			param->subj_ratio = idc[1].d;
			fprintf(stdout, "SUBJ_RATIO : %15.7e\n", param->subj_ratio);
		}else if( !strncmp("SUBJ_REL_ERROR", buf, 14) ){
			RC_TRY( sgetidc(buf, "cd", idc) );
			param->subj_rel_error = idc[1].d;
			fprintf(stdout, "SUBJ_REL_ERROR : %15.7e\n", param->subj_rel_error);
		}else if( !strncmp("OBJ_REL_ERROR", buf, 13) ){
			RC_TRY( sgetidc(buf, "cd", idc) );
			param->obj_rel_error = idc[1].d;
			fprintf(stdout, "OBJ_REL_ERROR : %15.7e\n\n", param->obj_rel_error);
		}
	}

	if(param->it_max <= 0){
		fprintf(stdout, "Invalid data > IT_MAX\n");
		return(ARG_ERROR_RC);
	}
	if(param->nr_max <= 0){
		fprintf(stdout, "Invalid data > NR_MAX\n");
		return(ARG_ERROR_RC);
	}
	if( (param->sol_flag < 0) || (param->sol_flag > 1) ){
		fprintf(stdout, "Invalid data > SOL_FLAG\n");
		return(ARG_ERROR_RC);
	}
	if( (param->opt_flag < 0) || (param->opt_flag > 1) ){
		fprintf(stdout, "Invalid data > OPT_FLAG\n");
		return(ARG_ERROR_RC);
	}
	if(param->max_str <= 0.0){
		fprintf(stdout, "Invalid data > MAX_STR\n");
		return(ARG_ERROR_RC);
	}
	if(param->subj_ratio <= 0.0){
		fprintf(stdout, "Invalid data > SUBJ_RATIO\n");
		return(ARG_ERROR_RC);
	}
	if(param->subj_rel_error <= 0.0){
		fprintf(stdout, "Invalid data > SUBJ_REL_ERROR\n");
		return(ARG_ERROR_RC);
	}
	if(param->obj_rel_error <= 0.0){
		fprintf(stdout, "Invalid data > OBJ_REL_ERROR\n");
		return(ARG_ERROR_RC);
	}

	return(NORMAL_RC);
}


RC
input_bdf (FILE *fp, NODE_ARRAY *node, ELEMENT_ARRAY *element,
           MATERIAL_PROP_ARRAY *material, PHYSICAL_PROP_ARRAY *physical,
           BC_ARRAY *rest, BC_ARRAY *force, DEFAULT_BC_ARRAY *b_force,
           BC_ARRAY *temp, DEFAULT_BC_ARRAY *def_temp,
           RIGID_ELEMENT_ARRAY *rigid)
{
	ELEM_BC_ARRAY elem_temp;

	RC_TRY( nst_input_node(fp, node) );
	fprintf(stdout, "Number of nodes : %d\n", node->size);

	RC_TRY( nst_input_element(fp, element) );
	fprintf(stdout, "Number of elements : %d\n", element->size);

	RC_TRY( nst_input_material_prop(fp, material) );
	fprintf(stdout, "Number of materials : %d\n", material->size);

	RC_TRY( nst_input_physical_prop(fp, physical) );
	fprintf(stdout, "Number of physical properties : %d\n", physical->size);

	RC_TRY( nst_input_restraint(fp, rest) );
	fprintf(stdout, "Number of restraints : %d\n", rest->size);

	RC_TRY( nst_input_force(fp, force, b_force) );
	fprintf(stdout, "Number of forces : %d\n", force->size);
	fprintf(stdout, "Number of b_force : %d\n", b_force->size);

	RC_TRY( nst_input_temperature(fp, temp, &elem_temp, def_temp) );
	fprintf(stdout, "Number of temperature on node : %d\n", temp->size);
	fprintf(stdout, "Number of default temperature : %d\n", def_temp->size);

	RC_TRY( nst_input_rigid_element(fp, rigid) );
	fprintf(stdout, "Number of rigid : %d\n\n", rigid->size);

	return(NORMAL_RC);
}


RC
output_bdf (FILE *fpw, NODE_ARRAY node, ELEMENT_ARRAY element,
            MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
            BC_ARRAY rest, BC_ARRAY force, DEFAULT_BC_ARRAY b_force)
{
	fprintf(fpw, "BEGIN BULK\n");
	RC_TRY( nst_output_node(fpw, node) );
	RC_TRY( nst_output_element(fpw, element) );
	RC_TRY( nst_output_material_prop(fpw, material) );
	RC_TRY( nst_output_physical_prop(fpw, physical) );
	RC_TRY( nst_output_restraint(fpw, rest) );
	RC_TRY( nst_output_force(fpw, force, b_force) );
	fprintf(fpw, "ENDDATA\n");

	return(NORMAL_RC);
}

RC
output_dxf (FILE *fpw, NODE_ARRAY node, ELEMENT_ARRAY element)
{
	ELEMENT_ARRAY surf;

	switch( analysis_dim(element) ){
	case 3:
		RC_TRY( extract_surface(element, &surf) );
		RC_TRY( output_dxf_polyline(fpw, surf, node) );
		/* RC_TRY( output_dxf_3dface(fpw, "result", 0, surf, node) ); */

		RC_TRY( free_element_array( &surf ) );
		break;
	case 2:
		RC_TRY( output_dxf_polyline(fpw, element, node) );
		/* RC_TRY( output_dxf_3dface(fpw, "result", 0, element, node) ); */
		break;
	default:
		return(WRITE_ERROR_RC);
	}

	return(NORMAL_RC);
}


RC
output_stl (FILE *fpw, NODE_ARRAY node, ELEMENT_ARRAY element)
{
	ELEMENT_ARRAY surf;

	switch( analysis_dim(element) ){
	case 3:
		RC_TRY( extract_surface(element, &surf) );
		RC_TRY( output_stl_3dface(fpw, "NOTITLE", surf, node) );
		RC_TRY( free_element_array( &surf ) );
		break;
	case 2:
		RC_TRY( output_stl_3dface(fpw, "NOTITLE", element, node) );
		break;
	default:
		return(WRITE_ERROR_RC);
	}

	return(NORMAL_RC);
}


RC
output_log (FILE *stream, int iteration, double obj0, double obj,
            double sbj0, double sbj)
{
	fprintf(stream, "%3d : obj %15.7e ( %6.2f %%) sbj %15.7e ( %6.2f %%)\n",
	        iteration, obj, 100.0*obj/obj0, sbj, 100.0*sbj/sbj0);
	fprintf(stdout, "%3d : obj %15.7e ( %6.2f %%) sbj %15.7e ( %6.2f %%)\n",
	        iteration, obj, 100.0*obj/obj0, sbj, 100.0*sbj/sbj0);
                
	return(NORMAL_RC);
}


RC
char_chain (char base_name[], char end_name[], char *result_fname)
{
	int ii1, ii2;

	for (ii1=0; base_name[ii1] != '\0'; ii1++){
		result_fname[ii1] = base_name[ii1];
	}
	for (ii2=0;  end_name[ii2] != '\0'; ii1++, ii2++){
		result_fname[ii1] =  end_name[ii2];
	}
	result_fname[ii1] = '\0';

	return(NORMAL_RC);
}



