#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "base.h"
#include "mathematics.h"
#include "fem.h"
#include "ccs_format.h"
#include "dmumps_c.h"
#include "topology_opt.h"

#include "omp.h"

#define SAFETY_STR_RATIO  (0.90)    /* SAFETY_STR_RATIO * MAX_DENS で概算 */
#if 1 /* 脊柱モデル用divisor(ステップ幅大) */
#define MAX_DIVISOR       (16.0)
//#define INIT_DIVISOR      (7.8125e-2)
#define INIT_DIVISOR      (7.8125e-3)
#endif
#if 0 /* 3次元柱モデル用divisor(ステップ幅小) */
#define MAX_DIVISOR       (256.0)
#define INIT_DIVISOR      (1.0)
#endif
#define DIV_DIVISOR       (1.0)
#define ARMIJO            (0.1)
#define INIT_SIZE         (1024)

RC set_density_array(NODE_ARRAY node, BC_ARRAY *density);
RC set_force_array(NODE_ARRAY node, ELEMENT_ARRAY element,
                   MATERIAL_PROP_ARRAY material,
                   PHYSICAL_PROP_ARRAY physical,
                   BC_ARRAY rest, BC_ARRAY *force,
                   BC_ARRAY density, BC_ARRAY *alpha,
                   NODE_ARRAY target_node,
                   ELEMENT_ARRAY target_surf,
                   double dens_pow);
RC set_density_rest(BC_ARRAY rest, BC_ARRAY *dens_rest);
RC get_alpha (int dim, int total_dof,
              NODE_ARRAY node, ELEMENT_ARRAY element,
              MATERIAL_PROP_ARRAY material,
              PHYSICAL_PROP_ARRAY physical,
              NODE_ARRAY target_node, ELEMENT_ARRAY target_surf,
              BC_ARRAY rest, BC_ARRAY force, BC_ARRAY density,
              BC_ARRAY alpha, double dens_pow);
RC get_target_vect_L2_norm(int dim, NODE_ARRAY node, NODE_ARRAY target_node,
                           ELEMENT_ARRAY target_surf, double *vect,
                           double *norm);
RC dist_u_inner_product_L2(int dim, NODE_ARRAY node, NODE_ARRAY target_node,
                           ELEMENT_ARRAY target_surf,
                           double *u_vect, double *L2_disp_dist);
RC cal_force (NODE_ARRAY node, NODE_ARRAY target_node,
              ELEMENT_ARRAY target_surf, BC_ARRAY *force,
              BC_ARRAY alpha);
RC extract_ref_point(NODE_ARRAY node, NODE_ARRAY *target_node,
		             ELEMENT_ARRAY elem, ELEMENT_ARRAY *target_surf);
RC input_model(FILE *fp, FILE *fp1, NODE_ARRAY *node, ELEMENT_ARRAY *element,
               MATERIAL_PROP_ARRAY *material,
               PHYSICAL_PROP_ARRAY *physical,
               BC_ARRAY *rest);
RC input_param (FILE *fp, PARAM *param);

static RC avs_printfe_model (FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY elem);
RC output_avs_bc_theta (FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY elem,
                        BC_ARRAY density);
static RC avs_printbc_theta (FILE *fp, BC_ARRAY density, NODE_ARRAY node);
RC output_avs_bc_density (FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY elem,
                          BC_ARRAY density);
RC output_avs_step_bc_density (FILE *fp_avs, NODE_ARRAY node, BC_ARRAY density);
static RC avs_printbc_density (FILE *fp, BC_ARRAY density, NODE_ARRAY node);
RC output_log (FILE *stream, int iteration, double obj0, double obj);
static RC output_avs_strain (FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY elem,
                             STRAIN_ARRAY strain);
static RC avs_printstrain (FILE *fp, ELEMENT_ARRAY elem, STRAIN_ARRAY strain);
static RC avs_printfe_model4discont_data (FILE *fp, NODE_ARRAY node,
                                          ELEMENT_ARRAY elem);
static RC avs_printelement (FILE *fp, int label, ELEM_TYPE type,
                            int node_num, const int *node);
static RC avs_printelem_type (FILE *fp, ELEM_TYPE type);
static double mises_strain (STRESS_STRAIN strain);
static double volume_strain (STRESS_STRAIN strain);
static int check_step (FILE *fp_avs);


int main(int argc, char **argv)
{

/* Setting runtime openmp thread counts */
//omp_set_num_threads(2); //mumps stills use all threads available

    FILE *fp, *fp1,  *fp_log;
	int ii1;
	
    NODE_ARRAY node;
    ELEMENT_ARRAY element;
    BC_ARRAY rest, dens_rest;
    BC_ARRAY force;
    BC_ARRAY density;
    BC_ARRAY alpha;
    MATERIAL_PROP_ARRAY material;
    PHYSICAL_PROP_ARRAY physical;
    DISP_ARRAY disp;
    STRAIN_ARRAY strain;

    NODE_ARRAY target_node_tmp;
    ELEMENT_ARRAY target_surf_tmp;

	NODE_ARRAY target_node;
	ELEMENT_ARRAY target_surf;


    PARAM param;

    /* Memory Manager */
    RC_TRY_MAIN( mm_init((size_t)800*1024*1024) );

    /* Log Print Level */
    RC_TRY_MAIN( set_log_level(5) );
    RC_TRY_MAIN( rc_fopen("iteration.log", "w", &fp_log) );
    RC_TRY_MAIN( set_log_file(0, stderr) );
    RC_TRY_MAIN( set_log_file(1, fp_log) );

    /* Input Model */
#if 0
	/* Same material */
	RC_TRY_MAIN( rc_fopen("./Input_data/model.bdf", "r", &fp) );
#endif
#if 1 
	/* Different material */
	RC_TRY_MAIN( rc_fopen("./Input_data/spine-20101215.bdf", "r", &fp) );
	RC_TRY_MAIN( rc_fopen("./Input_data/model_bone.bdf", "r", &fp1) );
#endif
    RC_TRY_MAIN( input_model(fp, fp1, &node, &element,
                             &material, &physical, &rest) );
    RC_TRY_MAIN( rc_fclose(fp) );
	RC_TRY_MAIN( rc_fclose(fp1) );
#if 0
    /* Input Target Node (binary file) */
    RC_TRY_MAIN( rc_fopen("./Input_data/target_node.bin", "r", &fp) );
    RC_TRY_MAIN( read_node_array(fp, &target_node) );
    RC_TRY_MAIN( rc_fclose(fp) );
    printf("Number of target nodes : %d\n", target_node.size);
    
    /* Input Target Surface (binary file) */
    RC_TRY_MAIN( rc_fopen("./Input_data/target_surf.bin", "r", &fp) );
    RC_TRY_MAIN( read_element_array(fp, &target_surf) );
    RC_TRY_MAIN( rc_fclose(fp) );
    printf("Number of target surfaces : %d\n", target_surf.size);
#endif

    disp.array = NULL;
    strain.array = NULL;
	density.array = NULL; 
#if 1
	/* Input Target Node and Surface (bdf file) */
	RC_TRY_MAIN( rc_fopen("./Input_data/renew_new_l1.bdf", "r", &fp) );
	RC_TRY_MAIN( nst_input_node(fp, &target_node_tmp) );
	RC_TRY_MAIN( nst_input_element(fp, &target_surf_tmp) );
	RC_TRY_MAIN( rc_fclose(fp) );
	RC_TRY_MAIN( extract_ref_point(target_node_tmp, &target_node, 
				                   target_surf_tmp, &target_surf) );
	printf("Number of target all nodes : %d\n", target_node_tmp.size);
	printf("Number of target all elemtns : %d\n", target_surf_tmp.size);
	printf("Number of target nodes : %d\n", target_node.size);
	printf("Number of target surfaces : %d\n", target_surf.size);
#endif

	/* Input PARAM file */
    RC_TRY_MAIN( rc_fopen("./Input_data/PARAM", "r", &fp) );
    RC_TRY_MAIN( input_param(fp, &param) );
    RC_TRY_MAIN( rc_fclose(fp) );
#if 1
    RC_TRY_MAIN( rc_fopen("./Output_data/reference_points.bdf", "w", &fp) );
    RC_TRY_MAIN( nst_output_node (fp, target_node) );
    RC_TRY_MAIN( nst_output_element (fp, target_surf) );
    RC_TRY_MAIN( rc_fclose(fp) );
#endif
#if 1
    /* Set Density (all default values are 0.0) */
    RC_TRY_MAIN( set_density_array(node, &density) );
#endif

	/* Set Force (force proportional to distance to reference points) */
    RC_TRY_MAIN( set_force_array(node, element, material, physical,
                                 rest, &force, density, &alpha,
                                 target_node, target_surf,
                                 param.dens_pow) );

    /* Set Density rest */
    RC_TRY_MAIN( set_density_rest(rest, &dens_rest) );


    /* Renumbering */
    RC_TRY_MAIN( dummy_renumber(&node) );


#if 0	
	RC_TRY_MAIN( rc_fopen("./Output_data/disp_binary.data", "r", &fp) );
	RC_TRY_MAIN( read_disp_array(fp, &disp) );
	RC_TRY_MAIN( rc_fclose(fp) );
	RC_TRY_MAIN( rc_fopen("./Output_data/density_binary.data", "r", &fp) );
	RC_TRY_MAIN( read_bc_array(fp, &density) );
	RC_TRY_MAIN( rc_fclose(fp) );
	
	for(ii1=0; ii1<node.size; ii1++){
		node.array[ii1].p.x += disp.array[ii1].v.x;
		node.array[ii1].p.y += disp.array[ii1].v.y;
		node.array[ii1].p.z += disp.array[ii1].v.z;
	}
#endif
#if 0
    RC_TRY_MAIN( rc_fopen("./Output_data/reference_points.bdf", "w", &fp) );
    RC_TRY_MAIN( nst_output_node (fp, target_node) );
    RC_TRY_MAIN( nst_output_element (fp, target_surf) );
    RC_TRY_MAIN( rc_fclose(fp) );
#endif	
    /* Strain Analysis with Geometrical Nonlinear */
    RC_TRY_MAIN( nl_strain_analysis(node, element, material, physical,
                                    rest, dens_rest, force, density, alpha,
                                    target_node, target_surf, &disp, &strain,
                                    param) );

    /* ----------Analysis Result------------- */
    /* Disp */
    RC_TRY_MAIN( rc_fopen("./Output_data/opt_disp.inp", "w", &fp) );
    RC_TRY_MAIN( output_avs_disp(fp, node, element, disp) );
    RC_TRY_MAIN( rc_fclose(fp) );
	
	RC_TRY_MAIN( rc_fopen("./Output_data/disp_binary.data", "w", &fp) );
	RC_TRY_MAIN( write_disp_array(fp, disp) );
	RC_TRY_MAIN( rc_fclose(fp) );

    /* Strain */
    RC_TRY_MAIN( rc_fopen("./Output_data/opt_strain.inp", "w", &fp) );
    RC_TRY_MAIN( output_avs_strain(fp, node, element, strain) );
    RC_TRY_MAIN( rc_fclose(fp) );

    /* Density */
    RC_TRY_MAIN( rc_fopen("./Output_data/opt_density.inp", "w", &fp) );
    RC_TRY_MAIN( output_avs_bc_density(fp, node, element, density) );
    RC_TRY_MAIN( rc_fclose(fp) );
	
	RC_TRY_MAIN( rc_fopen("./Output_data/density_binary.data", "w", &fp) );
	RC_TRY_MAIN( write_bc_array(fp, density) );
	RC_TRY_MAIN( rc_fclose(fp) );

    /* theta */
    RC_TRY_MAIN( rc_fopen("./Output_data/opt_theta.inp", "w", &fp) );
    RC_TRY_MAIN( output_avs_bc_theta(fp, node, element, density) );
    RC_TRY_MAIN( rc_fclose(fp) );

    /* Memory free */
    RC_TRY_MAIN( free_node_array(&node) );
    RC_TRY_MAIN( free_node_array(&target_node) );
    RC_TRY_MAIN( free_element_array(&element) );
    RC_TRY_MAIN( free_element_array(&target_surf) );
    RC_TRY_MAIN( free_material_prop_array(&material) );
    RC_TRY_MAIN( free_physical_prop_array(&physical) );
    RC_TRY_MAIN( free_bc_array(&rest) );
    RC_TRY_MAIN( free_bc_array(&dens_rest) );
    RC_TRY_MAIN( free_bc_array(&force) );
    RC_TRY_MAIN( free_bc_array(&alpha) );
    RC_TRY_MAIN( free_bc_array(&density) );
    RC_TRY_MAIN( free_disp_array(&disp) );
    RC_TRY_MAIN( free_strain_array(&strain) );

    /* Memory Manager */
    RC_TRY_MAIN( mm_terminate() );
    RC_TRY_MAIN( rc_fclose(fp_log) );

    return(0);
}

/*--------------------- Set Array -----------------------*/

/* density_arrayの確保
   とりあえず初期密度は一律の値としてある． */
RC
set_density_array(NODE_ARRAY node, BC_ARRAY *density)
{
    int ii1;

    RC_TRY( allocate_bc_array(node.size, density) );
    density->source = FEM_K_SOL;
    density->size = node.size;

    for(ii1=0; ii1<density->size; ii1++){
        density->array[ii1].node = node.array[ii1].label;
        density->array[ii1].s_type[0] = BC_DENSITY;
        density->array[ii1].s[0] = atanh(2.0*INIT_DENS - 1.0);
    }

    return(NORMAL_RC);
}

/* force_arrayの確保． 
   初めに重みalphaの値を1.0とし荷重で弾性変形問題の1ステップ分だけ計算し，
   求められた変位のノルムと初期の変位誤差ノルムの比を求める．
   その比がalphaの値となる． */
RC
set_force_array(NODE_ARRAY node, ELEMENT_ARRAY element,
                MATERIAL_PROP_ARRAY material,
                PHYSICAL_PROP_ARRAY physical,
                BC_ARRAY rest, BC_ARRAY *force,
                BC_ARRAY density, BC_ARRAY *alpha,
                NODE_ARRAY target_node, ELEMENT_ARRAY target_surf,
                double dens_pow)
{
    int ii1;
    int dim, total_dof;
    WORK_SET ws;

    /* 力forceと重みalphaの確保＆設定 */
    RC_TRY( allocate_work_set(&ws) );
    RC_TRY( allocate_bc_array(target_node.size, force) );
    RC_TRY( allocate_bc_array(target_node.size, alpha) );
    force->source = FEM_K_SOL;
    alpha->source = FEM_K_SOL;
#if 1    
	force->size = target_node.size;
    alpha->size = target_node.size;
#endif
	
    dim = analysis_dim(element);
    if(dim != 3) return(ARG_ERROR_RC);

    /* 全自由度 */
    RC_TRY( fem_total_dof(node, element, &total_dof) );

    for(ii1=0; ii1<force->size; ii1++){
        force->array[ii1].i_info[2] = 1;
        force->array[ii1].set_id = 1;
        force->array[ii1].node = target_node.array[ii1].label;
        force->array[ii1].i_info[0] = 0;
    }

    for(ii1=0; ii1<alpha->size; ii1++){
        alpha->array[ii1].node = target_node.array[ii1].label;
        alpha->array[ii1].s_type[0] = BC_WEIGHT;
        /* 仮の重みの値 */
        alpha->array[ii1].s[0] = 1.0;
    }

    /* 仮の力を求める　(alpha = 1.0) */
    RC_TRY( cal_force(node, target_node, target_surf, force, *alpha) );

    /* 仮の力で弾性変形問題を解き，求められた変位からalphaを求める． */
    RC_TRY( get_alpha(dim, total_dof, node, element, material, physical,
                      target_node, target_surf, rest, *force, density,
                      *alpha, dens_pow) );
    
    /* 求めたalphaを使ってもう一度力を求める */
    RC_TRY( cal_force(node, target_node, target_surf, force, *alpha) );

    RC_TRY( free_work_set(&ws) );

    return(NORMAL_RC);
}

/* H1勾配法用の拘束条件.弾性問題と同じ拘束場所にして値だけ変える */
RC
set_density_rest(BC_ARRAY rest, BC_ARRAY *dens_rest)
{
    int ii1;

    RC_TRY( allocate_bc_array(rest.size, dens_rest) );
    RC_TRY( copy_bc_array(rest, dens_rest) );

    for(ii1=0; ii1<dens_rest->size; ii1++){
        /* 拘束は1自由度に直す */
        dens_rest->array[ii1].v_type.x = BC_FIX;
        dens_rest->array[ii1].v_type.y = BC_FREE;
        dens_rest->array[ii1].v_type.z = BC_FREE;
        
        /* 拘束条件の値はとりあえず0.0にした */
        dens_rest->array[ii1].v.x = 0.0;
    }

    return(NORMAL_RC);
}

/* 荷重の係数を決定．
   弾性変形問題を1ステップ分だけ解き，
   得られた変位と初期変位誤差ノルムの比を計算し，
   その比をalphaの基準値とする */
RC
get_alpha (int dim, int total_dof,
           NODE_ARRAY node, ELEMENT_ARRAY element,
           MATERIAL_PROP_ARRAY material,
           PHYSICAL_PROP_ARRAY physical,
           NODE_ARRAY target_node, ELEMENT_ARRAY target_surf,
           BC_ARRAY rest, BC_ARRAY force, BC_ARRAY density,
           BC_ARRAY alpha, double dens_pow)
{
    CCS ccs;
    NONZERO_MATRIX_CCS tan_K;

    int ii1;
    double *f_vect;
    double *u_vect;
    double u_norm;
    double disp_dist;
    double temp_alpha;

    double omp_begin;
    double omp_end;

    /* MUMPS */
    DMUMPS_STRUC_C id;

    if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

    RC_TRY( allocate1D(total_dof, &f_vect) );
    RC_TRY( allocate1D(total_dof, &u_vect) );

    RC_TRY( add_nodal_force(dim, node, force, f_vect) );

    /* 弾性変形問題1ステップ分の計算 */
    RC_TRY( make_nonzero1a_ccs_sym_topology_robin (node, element, target_surf,
                                                   material, physical, density,
                                                   alpha, dens_pow, &tan_K) );
    RC_TRY( copy_nonzero_matrix1a_ccs2ccs_format(&tan_K, &ccs) );
    RC_TRY( free_nonzero_matrix1a_ccs(&tan_K) );
    RC_TRY( modify_ccs_format(dim, node, rest, &ccs, f_vect) );

printf("Clear!\n");

    omp_begin = omp_get_wtime();
    /* MUMPS */
    /* Inititalization */
    id.job = JOB_INIT;
    id.par = 1;
    /* Symmetric or Non-symmetric parameter */
    id.sym = 1;
    
 printf("Clear!\n");   
               //id.comm_fortran=USE_COMM_WORLD;
    dmumps_c(&id);
    /* Set Matrix data */
    id.n = total_dof;
    id.nz = ccs.array_size;
    id.irn = ccs.row_index;
    id.jcn = ccs.col_index;
 
 printf("Clear!\n");   
 
 
    for(ii1=0; ii1<id.nz; ii1++){
        id.irn[ii1]++;
        id.jcn[ii1]++;
    }
    id.a = ccs.array;
    id.rhs = f_vect;
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
 
 
    
    for(ii1=0; ii1<total_dof; ii1++){
        u_vect[ii1] = f_vect[ii1];
    }

    RC_TRY( free_nonzero_ccs(&ccs) );

    /* 弾性変形問題の解（変位）のL2ノルムを求める */
    RC_TRY( get_target_vect_L2_norm(dim, node, target_node, target_surf,
                                    u_vect, &u_norm) );
    RC_TRY( log_printf(5, " u_norm:%15.7e\n", u_norm) );

    /* 弾性変形問題の解（変位）と初期変位誤差との内積 */
    RC_TRY( dist_u_inner_product_L2(dim, node, target_node,
                                   target_surf, u_vect, &disp_dist) );
    RC_TRY( log_printf(5, " disp_dist:%15.7e\n", disp_dist) );

    /* 弾性変形問題の解（変位）と初期変位誤差からalphaを求める */
    temp_alpha = disp_dist / (u_norm * u_norm);
    RC_TRY( log_printf(5, " alpha:%15.7e\n", temp_alpha) );

    /* 元々決めてあった重み係数の値に掛け合わせる */
    for(ii1=0; ii1<alpha.size; ii1++){
        if(alpha.array[ii1].s_type[0] == BC_WEIGHT){
            alpha.array[ii1].s[0] *= temp_alpha;
        }
    }

    RC_TRY( mm_free(f_vect) );
    RC_TRY( mm_free(u_vect) );

    return(NORMAL_RC);
}

/* ベクトルのL2ノルムを求める */
RC
get_target_vect_L2_norm(int dim, NODE_ARRAY node, NODE_ARRAY target_node,
                        ELEMENT_ARRAY target_surf, double *vect, double *norm)
{
    int ii1, ii2, ii3;
    double sum = 0.0;
    WORK_SET ws;

    RC_TRY( allocate_work_set(&ws) );
    /* 既定変位までの距離を計算 */
    for(ii1=0; ii1<target_surf.size; ii1++){
        double *weights = ws.vec_G[0];
        double **gauss_points = ws.mat_G_4[0];
        double **node_location = ws.mat_N_3[0];
        int num_points;
        VECT3D f[MAX_NODE];
        ELEMENT surf_elem;

        if(target_surf.array[ii1].label < 0) continue;
        surf_elem = target_surf.array[ii1];
        if( (element_order(target_surf.array[ii1].type) == 2) ){
            /* 中間節点に荷重を振り分けない */
            RC_TRY( reduced_element_type(surf_elem.type, &(surf_elem.type),
                        &(surf_elem.node_num)) );
        }

        RC_TRY( Gauss_const_default(surf_elem.type,
                                    gauss_points, weights, &num_points) );
        RC_TRY( node_location_matrix(&(surf_elem), node, node_location) );

        for(ii2=0; ii2<num_points; ii2++){
            int index;
            double *N = ws.vec_N[0];
            VECT3D d_A;
            VECT3D disp_v;
            double disp_squa_norm;
            double det_J;

            /* 積分の準備 */
            RC_TRY( set_N_vector(surf_elem.type, gauss_points[ii2], N) );
            RC_TRY( cal_d_A(surf_elem, gauss_points[ii2], node_location,
                        &d_A, ws.mat_3_N[0], ws.mat_3_3[0]) );
            /* ヤコビアンdet_J */
            det_J = abs_vect3d(d_A);

            RC_TRY( init_vect3d(&disp_v) );
            /* 既定変位までの距離を求める */
            for(ii3=0; ii3<surf_elem.node_num; ii3++){
                index = node_renum_index(node, surf_elem.node[ii3]);
                disp_v.x += vect[index*dim] * N[ii3];
                disp_v.y += vect[index*dim+1] * N[ii3];
                disp_v.z += vect[index*dim+2] * N[ii3];
            }

            disp_squa_norm = inner_product3d(disp_v, disp_v);

            sum += weights[ii2] * disp_squa_norm * det_J;
        }
    }

    *norm = sqrt(sum);

    RC_TRY( free_work_set(&ws) );

    return(NORMAL_RC);
}

/* 変位誤差と変位ベクトルの内積を計算する */
RC
dist_u_inner_product_L2(int dim, NODE_ARRAY node, NODE_ARRAY target_node,
                        ELEMENT_ARRAY target_surf,
                        double *u_vect, double *L2_disp_dist)
{
    int ii1, ii2, ii3;
    double sum = 0.0;
    WORK_SET ws;

    RC_TRY( allocate_work_set(&ws) );
    /* 既定変位までの距離を計算 */
    for(ii1=0; ii1<target_surf.size; ii1++){
        double *weights = ws.vec_G[0];
        double **gauss_points = ws.mat_G_4[0];
        double **node_location = ws.mat_N_3[0];
        double **target_node_location = ws.mat_N_3[1];
        int num_points;
        VECT3D f[MAX_NODE];
        ELEMENT surf_elem;

        if(target_surf.array[ii1].label < 0) continue;
        surf_elem = target_surf.array[ii1];
        if( (element_order(target_surf.array[ii1].type) == 2) ){
            /* 中間節点に荷重を振り分けない */
            RC_TRY( reduced_element_type(surf_elem.type, &(surf_elem.type),
                        &(surf_elem.node_num)) );
        }

        RC_TRY( Gauss_const_default(surf_elem.type,
                                    gauss_points, weights, &num_points) );
        RC_TRY( node_location_matrix(&(surf_elem), target_node,
                                     target_node_location) );
        RC_TRY( node_location_matrix(&(surf_elem), node, node_location) );

        for(ii2=0; ii2<num_points; ii2++){
            int index;
            double *N = ws.vec_N[0];
            VECT3D d_A;
            VECT3D disp_v;
            VECT3D dist;
            double disp_dist;
            double det_J;

            /* 積分の準備 */
            RC_TRY( set_N_vector(surf_elem.type, gauss_points[ii2], N) );
            RC_TRY( cal_d_A(surf_elem, gauss_points[ii2], node_location,
                            &d_A, ws.mat_3_N[0], ws.mat_3_3[0]) );
            /* ヤコビアンdet_J */
            det_J = abs_vect3d(d_A);

            RC_TRY( init_vect3d(&disp_v) );
            RC_TRY( init_vect3d(&dist) );
            /* 既定変位までの距離を求める */
            for(ii3=0; ii3<surf_elem.node_num; ii3++){
                index = node_renum_index(node, surf_elem.node[ii3]);
                disp_v.x += u_vect[index*dim] * N[ii3];
                disp_v.y += u_vect[index*dim+1] * N[ii3];
                disp_v.z += u_vect[index*dim+2] * N[ii3];

                dist.x += (target_node_location[ii3][0] * N[ii3])
                        - (node_location[ii3][0] * N[ii3]);
                dist.y += (target_node_location[ii3][1] * N[ii3])
                        - (node_location[ii3][1] * N[ii3]);
                dist.z += (target_node_location[ii3][2] * N[ii3])
                        - (node_location[ii3][2] * N[ii3]);
            }

            disp_dist = inner_product3d(disp_v, dist);

            sum += weights[ii2] * disp_dist * det_J;
        }
    }

    *L2_disp_dist = sum;

    RC_TRY( free_work_set(&ws) );

    return(NORMAL_RC);
}

/* 参照点までの距離に比例した力を求める */
RC
cal_force (NODE_ARRAY node, NODE_ARRAY target_node,
           ELEMENT_ARRAY target_surf, BC_ARRAY *force,
           BC_ARRAY alpha)
{
    int ii1, ii2, ii3;
    WORK_SET ws;

    RC_TRY( allocate_work_set(&ws) );

    for(ii1=0; ii1<force->size; ii1++){
        force->array[ii1].v_type.x = BC_FORCE;
        force->array[ii1].v.x = 0.0;
        force->array[ii1].v_type.y = BC_FORCE;
        force->array[ii1].v.y = 0.0;
        force->array[ii1].v_type.z = BC_FORCE;
        force->array[ii1].v.z = 0.0;
    }

    for(ii1=0; ii1<target_surf.size; ii1++){
        double local_alpha;
        double *weights = ws.vec_G[0];
        double *node_alpha_v = ws.vec_N[0];
        double **gauss_points = ws.mat_G_4[0];
        double **node_location = ws.mat_N_3[0];
        double **target_node_location = ws.mat_N_3[1];
        int num_points;
        VECT3D f[MAX_NODE];
        ELEMENT surf_elem;

        if(target_surf.array[ii1].label < 0) continue;
        surf_elem = target_surf.array[ii1];
        if( (element_order(target_surf.array[ii1].type) == 2) ){
            /* 中間節点に荷重を振り分けない */
            RC_TRY( reduced_element_type(surf_elem.type, &(surf_elem.type),
                        &(surf_elem.node_num)) );
        }

        for(ii2=0; ii2<surf_elem.node_num; ii2++){
            node_alpha_v[ii2] = node_weight(alpha, surf_elem.node[ii2]);
        }

        RC_TRY( Gauss_const_default(surf_elem.type,
                                    gauss_points, weights, &num_points) );
        RC_TRY( node_location_matrix(&(surf_elem), node,
                                     node_location) );
        RC_TRY( node_location_matrix(&(surf_elem), target_node,
                                     target_node_location) );

        for(ii2=0; ii2<MAX_NODE; ii2++){
            f[ii2].x = f[ii2].y = f[ii2].z = 0.0;
        }

        for(ii2=0; ii2<num_points; ii2++){
            double *N = ws.vec_N[1];
            VECT3D d_A;
            VECT3D dist;
            double det_J;

            /* 積分の準備 */
            RC_TRY( set_N_vector(surf_elem.type, gauss_points[ii2], N) );
            RC_TRY( cal_d_A(surf_elem, gauss_points[ii2], node_location,
                            &d_A, ws.mat_3_N[0], ws.mat_3_3[0]) );
            /* ヤコビアンdet_J */
            det_J = abs_vect3d(d_A);

            RC_TRY( init_vect3d(&dist) );
            local_alpha = 0.0;
            /* 既定変位までの距離を求める */
            for(ii3=0; ii3<surf_elem.node_num; ii3++){
                dist.x += (target_node_location[ii3][0] * N[ii3])
                        - (node_location[ii3][0] * N[ii3]);
                dist.y += (target_node_location[ii3][1] * N[ii3])
                        - (node_location[ii3][1] * N[ii3]);
                dist.z += (target_node_location[ii3][2] * N[ii3])
                        - (node_location[ii3][2] * N[ii3]);

                local_alpha += node_alpha_v[ii3] * N[ii3];
            }

            for(ii3=0; ii3<surf_elem.node_num; ii3++){
                f[ii3].x += local_alpha * N[ii3] * dist.x
                          * weights[ii2] * det_J;
                f[ii3].y += local_alpha * N[ii3] * dist.y
                          * weights[ii2] * det_J;
                f[ii3].z += local_alpha * N[ii3] * dist.z
                          * weights[ii2] * det_J;
            }
        }

        for(ii2=0; ii2<surf_elem.node_num; ii2++){
            int index = search_bc_node(*force, surf_elem.node[ii2]);
            RC_NEG_CHK(index);
            force->array[index].v.x += f[ii2].x;
            force->array[index].v.y += f[ii2].y;
            force->array[index].v.z += f[ii2].z;
        }
    }

    RC_TRY( free_work_set(&ws) );

    return(NORMAL_RC);
}

/* 節点上の重み */
double
node_weight (BC_ARRAY weight, int node_label)
{
    int index;
    int ii1, ii2;

    /* 重みを定義していない節点では重みを1.0とする */
    index = search_bc_node(weight, node_label);
    if(index < 0) return(1.0);

    for(ii1=0; ii1<BC_SCALAR; ii1++){
        if(weight.array[index].s_type[ii1] == BC_WEIGHT){
            return(weight.array[index].s[ii1]);
        }
    }

    return(1.0);
}



/* -------------------------- Optimization ------------------------------- */

/* 密度（ヤング率）を設計変数とした変位誤差ノルム最小化問題 */
/* 最終的に求めたいものは、密度ではなく最適な密度分布のときの変位 */
RC nl_strain_analysis(NODE_ARRAY node, ELEMENT_ARRAY element,
                      MATERIAL_PROP_ARRAY material,
                      PHYSICAL_PROP_ARRAY physical,
                      BC_ARRAY rest, BC_ARRAY dens_rest, BC_ARRAY force,
                      BC_ARRAY density, BC_ARRAY alpha,
                      NODE_ARRAY target_node, ELEMENT_ARRAY target_surf,
                      DISP_ARRAY *disp, STRAIN_ARRAY *strain, PARAM param)
{
    FILE *fp, *fpd, *fp_log, *fp_logd, *fp_inte;
    int ii1;

    BC_ARRAY sens_trac, new_sens_trac;
    BC_ARRAY sens_dens;
    BC_ARRAY dens_new, dens_temp;
    double delta_s;
    double max_dens;
    double G_rho;
    double obj0 = 0.0, obj, obj_new = 0.0;
    double divisor = INIT_DIVISOR; 
    int set_delta_s_flag = 0;
    int armijo_result;

    if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

    RC_TRY( allocate_bc_array(density.size, &dens_new) );

    RC_TRY( allocate_bc_array(0, &new_sens_trac) );

    /* 各ステップごとに変化するファイルの準備 */
    fpd = fopen_avs_step_data("./Output_data/step_sens_density.inp",
                              node, element);

    RC_TRY( log_printf(5, "OPTIMIZATION ANALYSIS START !!\n") );
    /* 目的汎関数に対するsens_tracを計算 */
    RC_TRY( cal_obj_func(node, element, material, physical, rest,
                         force, density, alpha, target_node, target_surf,
                         disp, strain, &obj, &sens_trac, param) );
    
    /* 初期変位，初期ひずみを出力 */
    RC_TRY( log_printf(5, "---------OUTPUT INIT_DISP & INIT_STRAIN DATA--------\n") );
    RC_TRY( log_printf(5, "\n") );
    RC_TRY( rc_fopen("./Output_data/ini_disp.inp", "w", &fp) );
    RC_TRY( output_avs_disp(fp, node, element, *disp) );
    RC_TRY( rc_fclose(fp) );
    RC_TRY( rc_fopen("./Output_data/ini_strain.inp", "w", &fp) );
    RC_TRY( output_avs_strain(fp, node, element, *strain) );
    RC_TRY( rc_fclose(fp) );

    RC_TRY( free_disp_array(disp) );
    RC_TRY( free_strain_array(strain) );

    obj0 = obj;
    if(obj0 <= ABS_ERROR){
        RC_TRY( log_printf(5, "obj0 < \n") );
        return(CAL_ERROR_RC);
    }

    /* 密度変動を求める */
    RC_TRY( velocity_analysis_topology(node, element, physical, dens_rest, 1,
                                       &sens_trac, &sens_dens, param) );
    /* 各ステップごとの密度の出力の準備 */
    RC_TRY( output_avs_step_bc_density (fpd, node, sens_dens) );

    /* 収束履歴の出力 */
    RC_TRY( rc_fopen("./Output_data/iteration_history.dat", "w", &fp_log) );
    /* metapost用 */
    RC_TRY( rc_fopen("./Output_data/iteration_history.d", "w", &fp_logd) );

    /*----------反復計算開始------------*/
    for(ii1=0; ii1<param.it_max; ii1++){
    RC_TRY( log_printf(5, "-------------%3d -------------\n", ii1+1) );
    RC_TRY( output_log(fp_log, ii1+1, obj0, obj) );
    fprintf(fp_logd, "%d %.8lf \n", ii1, obj/obj0);

    /* 感度に比例した力と密度変動を掛け合わせる */
    RC_TRY( functional_variation_topology(sens_trac, sens_dens, &G_rho) );
    RC_TRY( log_printf(5, "G_rho: %15.7e\n", G_rho) );

    if(G_rho <= ABS_ERROR){
        RC_TRY( log_printf(5, "G_rho < \n") );
        return(CAL_ERROR_RC);
    }
    
    delta_s = -1.0;
    /* 仮の荷重係数における最大密度変動量を計算 */
    RC_TRY( max_density(1, &delta_s, node, element, &sens_dens, &max_dens) );
    if(max_dens <= ABS_ERROR){
        RC_TRY( log_printf(5, "max_dens < \n") );
        return(CAL_ERROR_RC);
    }

    /* 正式な荷重係数を求める */
    delta_s = -(SAFETY_STR_RATIO * param.max_dens)/(divisor * max_dens);
    RC_TRY( log_printf(5, "delta_s:%15.7e \n", delta_s) );

    /* 最大ひずみが設定値を越えていれば、再度修正 */
    RC_TRY( max_density(1, &delta_s, node, element, &sens_dens, &max_dens) );
    RC_TRY( log_printf(5, "max_dens:%15.7e\n", max_dens) );
    if(max_dens > (param.max_dens/divisor)){
        delta_s *= param.max_dens / (divisor * max_dens);
        RC_TRY( log_printf(5, "delta_s:%15.7e\n", delta_s) );
        RC_TRY( max_density(1, &delta_s, node, element, &sens_dens, &max_dens));
        RC_TRY( log_printf(5, "max_dens:%15.7e\n", max_dens) );
    }

    /* 密度変動量(sens_dens) の delta_s 倍を密度(dens_new)に付加 */
    RC_TRY( copy_bc_array(density, &dens_new) );
    RC_TRY( add_bc_density(dens_new, 1, &delta_s, &sens_dens) );

    RC_TRY( free_bc_array(&new_sens_trac) );

    /* 目的汎関数に対するsens_tracを計算 */
    RC_TRY( cal_obj_func(node, element, material, physical, rest,
                         force, dens_new, alpha, target_node, target_surf,
                         disp, strain, &obj_new, &new_sens_trac, param) );

    RC_TRY( log_printf(5, "%d: new obj =%15.7e (%8.1lf%%)\n",
            ii1+1, obj_new, obj_new/obj0*100) );

    RC_TRY( free_disp_array(disp) );
    RC_TRY( free_strain_array(strain) );
    
    /***** Armijoチェック *****/
    armijo_result = chk_armijo(obj, obj_new, G_rho, delta_s, &divisor, param);

    if(armijo_result == 0){
        RC_TRY( log_printf(5, "\n armijoで引っかかる \n") );
        ii1--;
        continue;
    }else if(armijo_result < 0){
        break;
    }

    /* 節点座標を更新 */
    RC_TRY( copy_bc_array(dens_new, &density) );
    RC_TRY( free_bc_array(&sens_trac) );
    sens_trac = new_sens_trac;
    RC_TRY( allocate_bc_array(0, &new_sens_trac) );
    
	obj = obj_new;

    /* 速度場解析 */
    RC_TRY( free_bc_array(&sens_dens) );
    RC_TRY( velocity_analysis_topology(node, element, physical, dens_rest, 1,
                                       &sens_trac, &sens_dens, param) );
    RC_TRY( output_avs_step_bc_density (fpd, node, sens_dens) );

    divisor /= DIV_DIVISOR;
    if(divisor < INIT_DIVISOR) divisor = INIT_DIVISOR;
    RC_TRY( log_printf(5, "divisor: %f\n", divisor) );
    }
    RC_TRY( output_log(fp_log, ii1+1, obj0, obj) );
    RC_TRY( rc_fclose(fp_log) );
    RC_TRY( rc_fclose(fp_logd) );
    /*-------------反復計算終了----------------*/

    DEFAULT_BC_ARRAY b_force;
    BC_ARRAY temp;
    DEFAULT_BC_ARRAY def_temp;
    
    temp.size = 0;
    temp.array = NULL;
    def_temp.size = 0;
    def_temp.array = NULL;
    b_force.size = 0;
    b_force.array = NULL;

    /* 最後に弾性変形問題を解いて最適分布における変位とひずみを求める */
    RC_TRY( fem_NL_solve_mumps_topology_robin(node, element, target_surf,
                                              material, physical, rest,
                                              force, density, alpha,
                                              b_force, temp, def_temp, disp,
                                              NULL, strain, param.dens_pow,
                                              param.sol_ss_flag) );

    if(armijo_result < 0){
        RC_TRY( log_printf(5, ": converged!\n") );
    }else{
        RC_TRY( log_printf(5, ": Not converged!\n") );
    }

    RC_TRY( free_bc_array(&sens_trac) );
    RC_TRY( free_bc_array(&new_sens_trac) );
    RC_TRY( free_bc_array(&sens_dens) );

    RC_TRY( rc_fclose(fpd) );

    RC_TRY( free_bc_array(&dens_new) );

    return(NORMAL_RC);
}

/* 順解析　→　目的汎関数の計算　→　随伴解析　→　感度の計算 */
RC
cal_obj_func(NODE_ARRAY node, ELEMENT_ARRAY element,
             MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
             BC_ARRAY rest, BC_ARRAY force, BC_ARRAY density, BC_ARRAY alpha,
             NODE_ARRAY target_node, ELEMENT_ARRAY target_surf,
             DISP_ARRAY *disp, STRAIN_ARRAY *strain,
             double *obj, BC_ARRAY *sens_trac, PARAM param)
{
    FILE *fp;
    int ii1;

    DEFAULT_BC_ARRAY b_force;
    BC_ARRAY temp;
    DEFAULT_BC_ARRAY def_temp;
    
    DISP_ARRAY adj_disp;

    temp.size = 0;
    temp.array = NULL;
    def_temp.size = 0;
    def_temp.array = NULL;
    b_force.size = 0;
    b_force.array = NULL;
    
    /* 順解析 */
    fprintf(stderr, "\n");
    fprintf(stderr, "---------FEM NL SOLVE--------\n");
    fprintf(stderr, "\n");

/* 荷重制御法使用 */
#if 1
    RC_TRY( fem_NL_solve_mumps_topology_robin_LCM(node, element, target_surf, 
                                                  material, physical, rest,
                                                  force, density, alpha,
                                                  b_force, temp, def_temp, disp,
                                                  NULL, strain, param.dens_pow,
                                                  param.sol_ss_flag, 40) );
#endif

#if 0
    RC_TRY( fem_NL_solve_mumps_topology_robin(node, element, target_surf, 
                                              material, physical, rest,
                                              force, density, alpha,
                                              b_force, temp, def_temp, disp,
                                              NULL, strain, param.dens_pow,
                                              param.sol_ss_flag) );
#endif

    /* 目的汎関数 */
    RC_TRY( integral_error_norm(node, target_node, target_surf,
                                *disp, alpha, obj) );
    RC_TRY( log_printf(5, "obj : %15.7e\n", *obj) );

    /* 随伴解析 */
    fprintf(stderr, "\n");
    fprintf(stderr, "---------ADJOINT ANALYSIS--------\n");
    fprintf(stderr, "\n");
    RC_TRY( adj_NL_solve_topology(node, element, material, physical,
                                  *disp, rest, alpha, density, param.dens_pow,
                                  target_node, target_surf, &adj_disp) );

    RC_TRY( allocate_bc_array(node.size, sens_trac) );
    sens_trac->size = node.size;
    for(ii1=0; ii1<sens_trac->size; ii1++){
        sens_trac->array[ii1].node = node.array[ii1].label;
    }
    RC_TRY( clean_bc_array(sens_trac) );

    /* G(sens_trac)の計算(H1勾配法の右辺項の計算) */
    RC_TRY( sens_trac_obj(node, element, *disp, adj_disp, material,
                          physical, density, param.dens_pow, sens_trac) );
    RC_TRY( compress_bc_array(sens_trac) );

    /* sens_tracの出力 */
    fprintf(stderr, "---------OUTPUT SENS_TRAC DATA--------\n");
    fprintf(stderr, "\n");
    RC_TRY( rc_fopen("./Output_data/sens_trac.inp", "w", &fp) );
    RC_TRY( output_avs_bc_force(fp, node, element, *sens_trac) );
    RC_TRY( rc_fclose(fp) );

    RC_TRY( free_disp_array(&adj_disp) );

    return(NORMAL_RC);
}

/* 変位誤差ノルムの積分値(目的汎関数) */
RC
integral_error_norm(NODE_ARRAY node, NODE_ARRAY target_node,
                    ELEMENT_ARRAY target_surf, DISP_ARRAY disp,
                    BC_ARRAY alpha, double *obj)
{
    int ii1, ii2, ii3;
    double sum = 0.0;
    WORK_SET ws;

    RC_TRY( allocate_work_set(&ws) );
    for(ii1=0; ii1<target_surf.size; ii1++){
        double local_alpha;
        double *node_alpha_v = ws.vec_N[0];
        double *weights = ws.vec_G[0];
        double **gauss_points = ws.mat_G_4[0];
        double **node_location = ws.mat_N_3[0];
        double **target_node_location = ws.mat_N_3[1];
        int num_points;
        ELEMENT surf_elem;

        /* gauss積分の準備 */
        if(target_surf.array[ii1].label < 0) continue;
        surf_elem = target_surf.array[ii1];
        if( (element_order(target_surf.array[ii1].type) == 2) ){
            /* 中間節点に荷重を振り分けない */
            RC_TRY( reduced_element_type(surf_elem.type, &(surf_elem.type),
                                         &(surf_elem.node_num)) );
        }

        for(ii2=0; ii2<surf_elem.node_num; ii2++){
            node_alpha_v[ii2] = node_weight(alpha, surf_elem.node[ii2]);
        }

        RC_TRY( Gauss_const_default(surf_elem.type,
                                    gauss_points, weights, &num_points) );
        RC_TRY( node_location_matrix(&(surf_elem), target_node,
                                     target_node_location) );
        RC_TRY( node_location_matrix(&(surf_elem), node, node_location) );

        for(ii2=0; ii2<num_points; ii2++){
            int index;
            double *N = ws.vec_N[1];
            VECT3D d_A;
            VECT3D disp_v;
            VECT3D dist;
            VECT3D disp_error;
            double disp_error_square;
            double det_J;

            /* 積分の準備 */
            RC_TRY( set_N_vector(surf_elem.type, gauss_points[ii2], N) );
            RC_TRY( cal_d_A(surf_elem, gauss_points[ii2], node_location,
                            &d_A, ws.mat_3_N[0], ws.mat_3_3[0]) );
            /* ヤコビアンdet_J */
            det_J = abs_vect3d(d_A);

            RC_TRY( init_vect3d(&disp_v) );
            RC_TRY( init_vect3d(&dist) );
            local_alpha = 0.0;
            /* 既定変位までの距離と順解析によって求められた変位の差を求める */
            for(ii3=0; ii3<surf_elem.node_num; ii3++){
                index = search_disp_node(disp, surf_elem.node[ii3]);
                disp_v.x += disp.array[index].v.x * N[ii3];
                disp_v.y += disp.array[index].v.y * N[ii3];
                disp_v.z += disp.array[index].v.z * N[ii3];

                dist.x += (target_node_location[ii3][0] * N[ii3])
                        - (node_location[ii3][0] * N[ii3]);
                dist.y += (target_node_location[ii3][1] * N[ii3])
                        - (node_location[ii3][1] * N[ii3]);
                dist.z += (target_node_location[ii3][2] * N[ii3])
                        - (node_location[ii3][2] * N[ii3]);
                
                local_alpha += node_alpha_v[ii3] * N[ii3];
            }

            disp_error.x = dist.x - disp_v.x;
            disp_error.y = dist.y - disp_v.y;
            disp_error.z = dist.z - disp_v.z;

            disp_error_square = inner_product3d(disp_error, disp_error);

            sum += local_alpha * disp_error_square * weights[ii2] * det_J;
        }
    }

    *obj = sum;

    RC_TRY( free_work_set(&ws) );

    return(NORMAL_RC);
}

/* 随伴解析(robin型であることに注意) */
RC
adj_NL_solve_topology(NODE_ARRAY node, ELEMENT_ARRAY element,
                      MATERIAL_PROP_ARRAY material,
                      PHYSICAL_PROP_ARRAY physical,
                      DISP_ARRAY disp, BC_ARRAY rest, BC_ARRAY alpha,
                      BC_ARRAY density, double dens_pow,
                      NODE_ARRAY target_node, ELEMENT_ARRAY target_surf,
                      DISP_ARRAY *adj_disp)
{
    CCS ccs;
    NONZERO_MATRIX_CCS tan_K;
    int ii1;
    int dim, total_dof;
    double *adj_f;
    double *adj_u;

    /* MUMPS */
    DMUMPS_STRUC_C id;

    dim = analysis_dim(element);
    if(dim != 3) return(ARG_ERROR_RC);

    if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

    /* 全自由度 */
    RC_TRY( fem_total_dof(node, element, &total_dof) );

    RC_TRY( allocate1D(total_dof, &adj_f) );
    RC_TRY( allocate1D(total_dof, &adj_u) );
    RC_TRY( vect2disp(dim, node, adj_u, adj_disp) );

    /* 随伴方程式における右辺項を計算 */
    RC_TRY( make_adj_force(dim, node, target_node, target_surf,
                           disp, alpha, adj_f) );

    /* 随伴方程式の剛性マトリックスを計算 */
    RC_TRY( make_adj_nonzero1a_ccs_sym(node, element, material,
                                       physical, disp, alpha, density,
                                       dens_pow, target_surf, &tan_K) );

    RC_TRY( copy_nonzero_matrix1a_ccs2ccs_format(&tan_K, &ccs) );
    RC_TRY( free_nonzero_matrix1a_ccs(&tan_K) );
    RC_TRY( modify_ccs_format(dim, node, rest, &ccs, adj_f) );

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
    for(ii1=0; ii1<id.nz; ii1++){
        id.irn[ii1]++;
        id.jcn[ii1]++;
    }
    id.a = ccs.array;
    id.rhs = adj_f;
    /* No outputs */
    id.icntl[0]=-1; id.icntl[1]=-1; id.icntl[2]=-1; id.icntl[3]=0;
    id.icntl[21]=0;
    /* Symbolic + Factorize + Solve */
    id.job = 6;
    dmumps_c(&id);
    /* End Job */
    id.job = JOB_END;
    dmumps_c(&id);

    for(ii1=0; ii1<total_dof; ii1++){
        adj_u[ii1] = adj_f[ii1];
    }

    RC_TRY( free_nonzero_ccs(&ccs) );
    RC_TRY( copy_vect2disp(dim, node, adj_u, *adj_disp, 1.0) );

    /* メモリ開放 */
    RC_TRY( mm_free(adj_f) );
    RC_TRY( mm_free(adj_u) );

    return(NORMAL_RC);
}

/* 随伴方程式における右辺項を計算 */
RC
make_adj_force(int dim, NODE_ARRAY node, NODE_ARRAY target_node,
               ELEMENT_ARRAY target_surf, DISP_ARRAY disp,
               BC_ARRAY alpha, double *f_vect)
{
    int ii1, ii2, ii3;
    WORK_SET ws;

    RC_TRY( allocate_work_set(&ws) );
    for(ii1=0; ii1<target_surf.size; ii1++){
        double local_alpha;
        double *node_alpha_v = ws.vec_N[0];
        double *weights = ws.vec_G[0];
        double **gauss_points = ws.mat_G_4[0];
        double **node_location = ws.mat_N_3[0];
        double **target_node_location = ws.mat_N_3[1];
        int num_points;
        VECT3D f[MAX_NODE];
        ELEMENT surf_elem;

        /* gauss積分の準備 */
        if(target_surf.array[ii1].label < 0) continue;
        surf_elem = target_surf.array[ii1];
        if( (element_order(target_surf.array[ii1].type) == 2) ){
            /* 中間節点に荷重を振り分けない */
            RC_TRY( reduced_element_type(surf_elem.type, &(surf_elem.type),
                        &(surf_elem.node_num)) );
        }

        RC_TRY( Gauss_const_default(surf_elem.type,
                                    gauss_points, weights, &num_points) );
        RC_TRY( node_location_matrix(&(surf_elem), node,
                                     node_location) );
        RC_TRY( node_location_matrix(&(surf_elem), target_node,
                                     target_node_location) );
        
        /* 節点上の荷重係数alphaを代入 */
        for(ii2=0; ii2<surf_elem.node_num; ii2++){
            node_alpha_v[ii2] = node_weight(alpha, surf_elem.node[ii2]);
        }

        for(ii2=0; ii2<MAX_NODE; ii2++){
            f[ii2].x = f[ii2].y = f[ii2].z = 0.0;
        }

        for(ii2=0; ii2<num_points; ii2++){
            double *N = ws.vec_N[1];
            VECT3D d_A;
            VECT3D dist;
            VECT3D disp_v;
            VECT3D disp_error;
            double det_J;

            /* 積分の準備 */
            RC_TRY( set_N_vector(surf_elem.type, gauss_points[ii2], N) );
            RC_TRY( cal_d_A(surf_elem, gauss_points[ii2], node_location,
                            &d_A, ws.mat_3_N[0], ws.mat_3_3[0]) );
            /* ヤコビアンdet_J */
            det_J = abs_vect3d(d_A);

            RC_TRY( init_vect3d(&dist) );
            RC_TRY( init_vect3d(&disp_v) );
            local_alpha = 0.0;
            /* 既定変位までの距離と順解析で求めた変位との差を求める */
            for(ii3=0; ii3<surf_elem.node_num; ii3++){
                int index = search_disp_node(disp, surf_elem.node[ii3]);
                disp_v.x += disp.array[index].v.x * N[ii3];
                disp_v.y += disp.array[index].v.y * N[ii3];
                disp_v.z += disp.array[index].v.z * N[ii3];

                dist.x += (target_node_location[ii3][0] * N[ii3])
                        - (node_location[ii3][0] * N[ii3]);
                dist.y += (target_node_location[ii3][1] * N[ii3])
                        - (node_location[ii3][1] * N[ii3]);
                dist.z += (target_node_location[ii3][2] * N[ii3])
                        - (node_location[ii3][2] * N[ii3]);
                
                local_alpha += node_alpha_v[ii3] * N[ii3];
            }

            disp_error.x = 2.0 * (disp_v.x - dist.x);
            disp_error.y = 2.0 * (disp_v.y - dist.y);
            disp_error.z = 2.0 * (disp_v.z - dist.z);

            /* gauss積分 */
            for(ii3=0; ii3<surf_elem.node_num; ii3++){
                f[ii3].x += local_alpha * N[ii3] * disp_error.x
                          * weights[ii2] * det_J;
                f[ii3].y += local_alpha * N[ii3] * disp_error.y
                          * weights[ii2] * det_J;
                f[ii3].z += local_alpha * N[ii3] * disp_error.z
                          * weights[ii2] * det_J;
            }
        }

        for(ii2=0; ii2<surf_elem.node_num; ii2++){
            int index = node_renum_index(node, surf_elem.node[ii2]);
            RC_NEG_CHK(index);
            f_vect[index*dim] += f[ii2].x;
            f_vect[index*dim+1] += f[ii2].y;
            f_vect[index*dim+2] += f[ii2].z;
        }
    }

    RC_TRY( free_work_set(&ws) );

    return(NORMAL_RC);
}

/* 随伴方程式における剛性マトリックスの作成 */
RC
make_adj_nonzero1a_ccs_sym (NODE_ARRAY node, ELEMENT_ARRAY element,
                            MATERIAL_PROP_ARRAY material,
                            PHYSICAL_PROP_ARRAY physical,
                            DISP_ARRAY disp, BC_ARRAY alpha,
                            BC_ARRAY density, double dens_pow,
                            ELEMENT_ARRAY target_surf,
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

    /* 随伴方程式の左辺項の計算 */
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
    
    /* 方程式の右辺第二項を移行した式の計算 */
    for(ii1=0; ii1<target_surf.size; ii1++){
        if(target_surf.array[ii1].label < 0) continue;

        RC_TRY( make_NtN_elem_matrix(dim, target_surf.array[ii1], node,
                                     physical, alpha, ws, &elem_matrix) );
        
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

/* NtNの項の計算 */
RC
make_NtN_elem_matrix(int dim, ELEMENT surf_elem, NODE_ARRAY node,
                     PHYSICAL_PROP_ARRAY physical, BC_ARRAY alpha,
                     WORK_SET ws, ELEM_MATRIX *elem_matrix)
{
    int ii1, ii2, ii3;
    int g_num_point;
    double thickness;
    double local_alpha;
    double *node_alpha_v = ws.vec_N[0];
    double *N_vec = ws.vec_N[1];
    double *g_weights = ws.vec_G[0];
    double **N_matrix = ws.mat_3_3N[0];
    double **node_location = ws.mat_N_3[0];
    double **g_points = ws.mat_G_4[0];
    double **NtN;


    if((elem_matrix->dim = dim) <= 0){
        return(ARG_ERROR_RC);
    }

    if(surf_elem.label < 0) return(ARG_ERROR_RC);
    elem_matrix->element = surf_elem.label;
    elem_matrix->size = (surf_elem.node_num) * (elem_matrix->dim);

    /* elem_matrix->array[][] の確保 */
    RC_TRY( allocate2D(elem_matrix->size, elem_matrix->size,
                       &(elem_matrix->matrix)) );
    RC_TRY( allocate2D(elem_matrix->size, elem_matrix->size, &NtN) );

    /* elem_matirx->index[] の確保と代入 */
    elem_matrix->index = mm_alloc(sizeof(int)*(elem_matrix->size));
    if(elem_matrix->index == NULL) return(ALLOC_ERROR_RC);

    for(ii1=0; ii1<surf_elem.node_num; ii1++){
        int index = node_renum_index(node, surf_elem.node[ii1]);
        if(index < 0) return(SEEK_ERROR_RC);

        for(ii2=0; ii2<(elem_matrix->dim); ii2++){
            elem_matrix->index[ii1*(elem_matrix->dim)+ii2]
                               = index*(elem_matrix->dim) + ii2;
        }
        /* 節点上の荷重係数(alpha)の計算 */
        node_alpha_v[ii1] = node_weight(alpha, surf_elem.node[ii1]);
    }

    RC_TRY( node_location_matrix(&surf_elem, node, node_location) );
    RC_TRY( Gauss_const_default(surf_elem.type, g_points, g_weights,
                                &g_num_point) );
//  RC_TRY( get_thickness(surf_elem, physical, &thickness) );
    thickness = 1.0;

    for(ii1=0; ii1<g_num_point; ii1++){
        VECT3D d_A;
        double det_J;

        RC_TRY( set_N_vector(surf_elem.type, g_points[ii1], N_vec) );
        RC_TRY( cal_d_A(surf_elem, g_points[ii2], node_location, &d_A,
                        ws.mat_3_N[0], ws.mat_3_3[0]) );
        /* ヤコビアンdet_J */
        det_J = abs_vect3d(d_A);

        local_alpha = 0.0;

        /* NtNの計算 */
        for(ii2=0; ii2<surf_elem.node_num; ii2++){
            N_matrix[0][ii2*3+0] = N_vec[ii2];
            N_matrix[0][ii2*3+1] = 0.0;
            N_matrix[0][ii2*3+2] = 0.0;
            N_matrix[1][ii2*3+0] = 0.0;
            N_matrix[1][ii2*3+1] = N_vec[ii2];
            N_matrix[1][ii2*3+2] = 0.0;
            N_matrix[2][ii2*3+0] = 0.0;
            N_matrix[2][ii2*3+1] = 0.0;
            N_matrix[2][ii2*3+2] = N_vec[ii2];
                
            local_alpha += node_alpha_v[ii2] * N_vec[ii2];
        }
        init_matrix(elem_matrix->size, elem_matrix->size, NtN);
        mul_matrix_AtB(dim, elem_matrix->size, elem_matrix->size,
                       N_matrix, N_matrix, NtN);

        for(ii2=0; ii2<(elem_matrix->size); ii2++){
            for(ii3=0; ii3<(elem_matrix->size); ii3++){
                elem_matrix->matrix[ii2][ii3]
                += local_alpha*NtN[ii2][ii3]*det_J*g_weights[ii1]*thickness;
            }
        }
    }

    RC_TRY( free2D(elem_matrix->size, elem_matrix->size, &NtN) );

    return(NORMAL_RC);
}

/* sens_trac(H1勾配法の右辺項)の計算 */
RC
sens_trac_obj(NODE_ARRAY node, ELEMENT_ARRAY element,
              DISP_ARRAY disp, DISP_ARRAY adj_disp,
              MATERIAL_PROP_ARRAY material, PHYSICAL_PROP_ARRAY physical,
              BC_ARRAY density, double dens_pow,
              BC_ARRAY *sens_trac)
{
    int ii1, ii2, ii3;
    WORK_SET ws1, ws2;

    RC_TRY( allocate_work_set(&ws1) );
    RC_TRY( allocate_work_set(&ws2) );

    for(ii1=0; ii1<element.size; ii1++){
        ELEMENT elem;
        double local_density;
        double local_phi;
        double f[MAX_NODE];
        double *node_density_vect = ws1.vec_N[0];
        double *weights = ws1.vec_G[0];
        double **gauss_points = ws1.mat_G_4[0];
        double **node_location = ws1.mat_N_3[0];
        int num_points;

        if(element.array[ii1].label < 0) continue;
        elem = element.array[ii1];

        RC_TRY( Gauss_const_max(elem.type, gauss_points,
                                weights, &num_points) );
        RC_TRY( node_location_matrix(&elem, node, node_location) );
        
        for(ii2=0; ii2<elem.node_num; ii2++){
            /* 節点上の密度のベクトルを作成 */
            node_density_vect[ii2] = node_density(density, elem.node[ii2]);
        }
        
        for(ii2=0; ii2<MAX_NODE; ii2++){
            f[ii2] = 0.0;
        }

        for(ii2=0; ii2<num_points; ii2++){
            double *N = ws1.vec_N[1];
            STRESS_STRAIN local_stress, local_strain;
            STRESS_STRAIN adj_local_stress, adj_local_strain;
            double det_J;
            double G;

            RC_TRY( set_N_vector(elem.type, gauss_points[ii2], N) );
            RC_TRY( cal_det_J(elem, gauss_points[ii2], node_location,
                              &det_J, ws1.mat_3_N[0], ws1.mat_3_3[0]) );

            /* 各節点上の密度の計算 */
            local_density = 0.0;
            for(ii3=0; ii3<elem.node_num; ii3++){
                local_density += node_density_vect[ii3] * N[ii3];
            }
            local_phi = DENS_UPPER_LIMIT
                      * pow( (tanh(local_density)+1.0)*0.5 , dens_pow - 1.0 );

            /* S(u)の計算*/
            RC_TRY( local_PK2_stress_Green_strain(elem, gauss_points[ii2],
                                                  node, disp, material,
                                                  physical, &local_strain,
                                                  &local_stress,
                                                  ws2) );
#if 1
            /* δE(u,v)の計算  */
            RC_TRY( local_delta_Green_strain(elem, gauss_points[ii2],
                                             node, disp, adj_disp,
                                             ws2, &adj_local_strain) );
#endif
            /* 感度G */
            G = - dens_pow * local_phi * 0.5 
              / (cosh(local_density) * cosh(local_density))
              * ( local_stress.v.x * adj_local_strain.v.x
              + local_stress.v.y * adj_local_strain.v.y
              + local_stress.v.z * adj_local_strain.v.z
              + local_stress.v.yz * adj_local_strain.v.yz
              + local_stress.v.zx * adj_local_strain.v.zx
              + local_stress.v.xy * adj_local_strain.v.xy);
            
            for(ii3=0; ii3<elem.node_num; ii3++){
                f[ii3] += weights[ii2] * G * N[ii3] * det_J;
            }
        }

        for(ii2=0; ii2<elem.node_num; ii2++){
            int index = search_bc_node(*sens_trac, elem.node[ii2]);
            if(index < 0) return(SEEK_ERROR_RC);
            sens_trac->array[index].v_type.x = BC_FORCE;
            sens_trac->array[index].v.x += f[ii2];
        }
    }

    RC_TRY( free_work_set(&ws1) );
    RC_TRY( free_work_set(&ws2) );

    return(NORMAL_RC);
}

/* δE(u,v)の計算  */
RC
local_delta_Green_strain(ELEMENT element, const double *local_points,
                         NODE_ARRAY node, DISP_ARRAY disp,
                         DISP_ARRAY adj_disp,
                         WORK_SET ws,
                         STRESS_STRAIN *Green_strain)
{
    int dim, ssnum;
    double det_J;
    double Green_strain_v[6];
    double adj_disp_v[3*MAX_NODE];
    double *g_weights = ws.vec_G[0];
    double **node_location = ws.mat_N_3[0];
    double **B_matrix = ws.mat_6_3N[0];
    double **A_matrix = ws.mat_9_9[0];
    double **G_matrix = ws.mat_9_3N[0];
    double **AG = ws.mat_6_3N[1];
    double **g_points = ws.mat_G_4[0];

    if(element.label < 0) return(ARG_ERROR_RC);

    dim = element_dim(element.type);
    ssnum = element_ssnum(element.type);
    init_stress_strain(Green_strain);

    RC_TRY( fill_disp_vector(element, adj_disp, adj_disp_v) );
    RC_TRY( node_location_matrix(&element, node, node_location) );

    RC_TRY( make_B_matrix(element, local_points, node_location, B_matrix,
                          &det_J) );
    RC_TRY( make_AG_matrix(element, local_points, ws, node_location, disp,
                           A_matrix, G_matrix) );
    mul_matrix(ssnum,dim*dim,dim*element.node_num,A_matrix,G_matrix,AG);
    RC_TRY( add_matrix(ssnum, dim*element.node_num, 1.0, B_matrix, 1.0, AG,
                       B_matrix) );
    mul_matrix_vect(ssnum, dim*element.node_num, B_matrix, adj_disp_v,
                    Green_strain_v);
    RC_TRY( vect2stress_strain(ssnum, Green_strain_v, Green_strain) );

    return(NORMAL_RC);
}
#if 1
/* mumps使用 */
/* 密度変動量の計算 */
RC
velocity_analysis_topology (NODE_ARRAY node, ELEMENT_ARRAY element,
                            PHYSICAL_PROP_ARRAY physical,
                            BC_ARRAY dens_rest, int num_case,
                            const BC_ARRAY tractions[], BC_ARRAY dens[],
                            PARAM param)
{
    NONZERO_MATRIX_CCS matrix;
    CCS ccs;

    double **f_vect;
    double *d_vect;
    int ii1, ii2;
    int dim, total_dof;

    /* MUMPS */
    DMUMPS_STRUC_C id;


    fprintf(stderr, "---------VELOCITY ANALYSIS--------\n");
    fprintf(stderr, "\n");
    f_vect = (double **)mm_alloc(num_case * sizeof(double *));
    if(f_vect == NULL ) return(ALLOC_ERROR_RC);

    /* 剛性マトリックス作成 */
    RC_TRY( make_ccs_matrix_topology(node, element, physical,
                                     &matrix, param) );
    
    /* 荷重ベクトル作成 */
    dim = DIM_DENSITY;
    total_dof = dim * node.size;
    for(ii1=0; ii1<num_case; ii1++){
        RC_TRY( fem_allocate_vector(dim, node, &(f_vect[ii1])) );
        RC_TRY( add_nodal_force(dim, node, tractions[ii1], f_vect[ii1]) );
    }
    RC_TRY( allocate1D(total_dof, &d_vect) );

    RC_TRY( copy_nonzero_matrix1a_ccs2ccs_format(&matrix, &ccs) );
    RC_TRY( free_nonzero_matrix1a_ccs(&matrix) );
    /* 拘束条件組み込み */
    RC_TRY( modify_ccs_format_n(dim, node, dens_rest, &ccs, f_vect, num_case) );
    /* 密度変動量の計算 */
    for(ii1=0; ii1<num_case; ii1++){
        /* MUMPS */
        /* Inititalization */
        id.job = JOB_INIT;
        id.par = 1;
        /* Symmetric or Non-symmetric parameter */
        id.sym = 1;
        //   id.comm_fortran=USE_COMM_WORLD;
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
        id.rhs = f_vect[ii1];

        /* No outputs */
        id.icntl[0]=-1; id.icntl[1]=-1; id.icntl[2]=-1; id.icntl[3]=0;

        /* OOC(Out-Of-Core) option */
        /* 0:In Core,  1:Out-Of-Core */
        id.icntl[21]=0;

        /* Symbolic + Factorize + Solve */
        id.job = 6;
        dmumps_c(&id);

        /* End Job */
        id.job = JOB_END;
        dmumps_c(&id);
        
        for(ii2=0; ii2<total_dof; ii2++){
            d_vect[ii2] = f_vect[ii1][ii2];
        }

        RC_TRY( vect2bc_density(node, d_vect, &(dens[ii1])) );
    }
    /* マトリックスの開放 */
    RC_TRY( free_nonzero_ccs(&ccs) );
    RC_TRY( free2D(num_case, total_dof, &f_vect) );
    RC_TRY( mm_free(d_vect) );

    return(NORMAL_RC);
}
#endif

/* 密度最適化用H1勾配法の左辺マトリックスの計算 */
RC
make_ccs_matrix_topology (NODE_ARRAY node, ELEMENT_ARRAY element,
                          PHYSICAL_PROP_ARRAY physical,
                          NONZERO_MATRIX_CCS *matrix, PARAM param)
{
    int ii1, ii2, ii3;
    int dim = DIM_DENSITY;
    ELEM_MATRIX elem_matrix;
    WORK_SET ws;

    RC_TRY( allocate_work_set(&ws) );
    RC_TRY( allocate_nonzero_matrix1a_ccs(node.size*dim, matrix) );
    for(ii1=0; ii1<element.size; ii1++){
        double *ptr;

        if(element.array[ii1].label < 0) continue;

        RC_TRY( make_H1_elem_matrix_topology(element.array[ii1], node,
                                             physical, &elem_matrix,
                                             ws, param.c0, param.c1) );
        for(ii2=0; ii2<elem_matrix.size; ii2++){
            for(ii3=0; ii3<elem_matrix.size; ii3++){
                if(elem_matrix.index[ii3] > elem_matrix.index[ii2]) continue;
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

/* 密度最適化用H1勾配法の要素マトリックスの計算 */
RC
make_H1_elem_matrix_topology(ELEMENT element, NODE_ARRAY node,
                             PHYSICAL_PROP_ARRAY physical,
                             ELEM_MATRIX *elem_matrix,
                             WORK_SET ws, double c0, double c1)
{
    int ii1, ii2, ii3;
    int num_points;
    int dim;
    double det_J;
    double thickness;
    double *N_vec = ws.vec_N[0];
    double *weights = ws.vec_G[0];
    double **dN = ws.mat_3_N[0];
    double **dNxyz = ws.mat_3_N[1];
    double **Jacobi = ws.mat_3_3[0];
    double **inverse_J = ws.mat_3_3[1];
    double **node_location = ws.mat_N_3[0];
    double **gauss_points = ws.mat_G_4[0];
    double **dNtdN;
    double **NtN;

    if((dim = element_dim(element.type)) <= 0) return(ARG_ERROR_RC);

    elem_matrix->dim = DIM_DENSITY; /* 密度変動は1自由度 */
    if(element.label < 0) return(ARG_ERROR_RC);
    elem_matrix->element = element.label;
    elem_matrix->size = (element.node_num) * (elem_matrix->dim);

    /* elem_matrix->array[][] の確保 */
    RC_TRY( allocate2D(elem_matrix->size,
                       elem_matrix->size, &(elem_matrix->matrix)) );
    RC_TRY( allocate2D(elem_matrix->size, elem_matrix->size, &dNtdN) );
    RC_TRY( allocate2D(elem_matrix->size, elem_matrix->size, &NtN) );

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
    }

    RC_TRY( Gauss_const_default(element.type, gauss_points,
                                weights, &num_points) );
    RC_TRY( get_thickness(element, physical, &thickness) );

    RC_TRY( node_location_matrix(&element, node, node_location) );

    /* c0 * dNtdN + c1 * NtNの計算 */
    for(ii1=0; ii1<num_points; ii1++){
        RC_TRY( set_dN_matrix(element.type, gauss_points[ii1], dN) );
        mul_matrix(dim, element.node_num, dim, dN, node_location, Jacobi);
        RC_TRY( inverse_matrix(dim, Jacobi, inverse_J) );
        mul_matrix(dim, dim, element.node_num, inverse_J, dN, dNxyz);

        mul_matrix_AtB(dim, element.node_num, element.node_num,
                       dNxyz, dNxyz, dNtdN);

        RC_TRY( set_N_vector(element.type, gauss_points[ii1], N_vec) );

        RC_TRY( tensor_product(element.node_num, element.node_num,
                               N_vec, N_vec, NtN) );

        det_J = determinant(dim, Jacobi);

        for(ii2=0; ii2<elem_matrix->size; ii2++){
            for(ii3=0; ii3<elem_matrix->size; ii3++){
                elem_matrix->matrix[ii2][ii3]
                    += (c0 * dNtdN[ii2][ii3] + c1 * NtN[ii2][ii3])
                     * weights[ii1] * det_J * thickness;
            }
        }
    }

    RC_TRY( free2D(elem_matrix->size, elem_matrix->size, &dNtdN) );
    RC_TRY( free2D(elem_matrix->size, elem_matrix->size, &NtN) );

    return(NORMAL_RC);
}

/* 拘束条件が設けられている節点の判別 */
RC
modify_ccs_format_n (int dim, NODE_ARRAY node, BC_ARRAY rest, CCS *ccs,
                     double **f_vect, int vect_size)
{
    int ii1;
    int index;

    for(ii1=0; ii1<rest.size; ii1++){
        if(rest.array[ii1].node < 0) continue;
        index = node_renum_index(node, rest.array[ii1].node);
        RC_NEG_CHK(index);

        if(rest.array[ii1].v_type.x == BC_FIX){
            RC_TRY( modify_ccs_format_dof_n(dim*index + 0, rest.array[ii1].v.x,
                                            ccs, f_vect, vect_size) );
        }
        if(rest.array[ii1].v_type.y == BC_FIX){
            RC_TRY( modify_ccs_format_dof_n(dim*index + 1, rest.array[ii1].v.y,
                                            ccs, f_vect, vect_size) );
        }
        if(dim == 3){
            if(rest.array[ii1].v_type.z == BC_FIX){
                RC_TRY( modify_ccs_format_dof_n(dim*index + 2,
                                                rest.array[ii1].v.z, ccs,
                                                f_vect, vect_size) );
            }
        }
    }

    return(NORMAL_RC);
}

/* 拘束条件による荷重ベクトル，剛性マトリックスの修正 */
RC
modify_ccs_format_dof_n (int index_xyz, double rest_value, CCS *ccs,
                         double **f_vect, int vect_size)
{
    int ii1, ii2;
    int index;


    for(ii1=ccs->col_offset[index_xyz];
            ii1<ccs->col_offset[index_xyz+1]; ii1++){
        for(ii2=0; ii2<vect_size; ii2++){
            f_vect[ii2][ccs->row_index[ii1]] -= ccs->array[ii1] * rest_value;
        }
        ccs->array[ii1] = 0.0;
    }

    for(ii1=0; ii1<vect_size; ii1++){
        f_vect[ii1][index_xyz] = rest_value;
    }

    for(ii1=0; ii1<ccs->size; ii1++){
        if(ccs->row_index[ccs->col_offset[ii1]] > index_xyz
           || ccs->row_index[ccs->col_offset[ii1+1]-1] < index_xyz) continue;

        index = search_nonzero_ccs_array_index(*ccs, index_xyz, ii1);
        if(index < 0) continue;

        ccs->array[index] = 0.0;
    }

    index = search_nonzero_ccs_array_index(*ccs, index_xyz, index_xyz);
    if(index < 0) return(SEEK_ERROR_RC);
    ccs->array[index] = 1.0;

    return(NORMAL_RC);
}

/* 1次元配列に格納してある密度をBC_ARRAY型に格納する */
RC
vect2bc_density (NODE_ARRAY node, const double *d_vect, BC_ARRAY *density)
{
    int ii1;
    int index;

    RC_TRY( allocate_bc_array(node.size, density) );
    density->source = FEM_K_SOL;
    for(ii1=0; ii1<node.size; ii1++){
        if(node.array[ii1].label < 0){
            density->array[ii1].node = -1;
            continue;
        }
        index = node.array[ii1].renum;
        if((index < 0)||(index >= node.size)) return(ARG_ERROR_RC);

        density->array[ii1].node = node.array[ii1].label;
        density->array[ii1].s_type[0] = BC_DENSITY;
        density->array[ii1].s[0] = d_vect[index];

    }
    density->size = node.size;

    RC_TRY( clean_bc_array(density) );

    return(NORMAL_RC);
}

/* sens_tracとsens_densを掛け合わせる */
RC
functional_variation_topology(BC_ARRAY sens_force, BC_ARRAY sens_dens,
                              double *delta_f)
{
    int ii1, ii2, ii3;
    int index;
    long double sum = 0.0;

    for(ii1=0; ii1<sens_force.size; ii1++){
        if(sens_force.array[ii1].node < 0) continue;

        index = search_bc_node(sens_dens, sens_force.array[ii1].node);
        if(index < 0) return(SEEK_ERROR_RC);

        if(sens_force.array[ii1].v_type.x == BC_FORCE){
            for(ii2=0; ii2<BC_SCALAR; ii2++){
                if(sens_dens.array[index].s_type[ii2] == BC_DENSITY){
                    sum += sens_force.array[ii1].v.x
                         * sens_dens.array[index].s[ii2];
                    break;
                }
            }
        }

    }

    *delta_f = (double)sum;

    return(NORMAL_RC);
}

/* 全節点上の密度の中で最大の密度の値を求める */
RC
max_density(int dens_num, const double factors[],
            NODE_ARRAY node, ELEMENT_ARRAY element,
            const BC_ARRAY density[], double *max_dens)
{
    int ii1, ii2, ii3;
    double temp_dens;

    *max_dens=0.0;

    for(ii1=0; ii1<node.size; ii1++){
        temp_dens = 0.0;

        for(ii2=0; ii2<dens_num; ii2++){
            int index = search_bc_node(density[ii2], node.array[ii1].label);
            if(index < 0) continue;

            for(ii3=0; ii3<BC_SCALAR; ii3++){
                if(density[ii2].array[index].s_type[ii3] == BC_DENSITY){
                    temp_dens += factors[ii2]*density[ii2].array[index].s[ii3];
                    break;
                }
            }
        }
        temp_dens = fabs(temp_dens);
        *max_dens = MAX2(temp_dens, *max_dens);
    }

    return(NORMAL_RC);
}

/* 密度変動量を現密度に足す */
RC
add_bc_density(BC_ARRAY density_total, int dens_num,
               const double *factors, const BC_ARRAY *densities)
{
    int ii1, ii2, ii3;
    int index;
    double *sum;

    sum = (double *)mm_alloc((density_total.size)*sizeof(double));
    if(sum == NULL) return(ALLOC_ERROR_RC);

    for(ii1=0; ii1<density_total.size; ii1++){
        sum[ii1] = 0.0;
    }

    /* 桁落ちを避けるために、加算する変位を一時的に sum に格納 */
    for(ii1=0; ii1<dens_num; ii1++){
        for(ii2=0; ii2<densities[ii1].size; ii2++){
            if(densities[ii1].array[ii2].node < 0) continue;

            for(ii3=0; ii3<BC_SCALAR; ii3++){
                if(densities[ii1].array[ii2].s_type[ii3] == BC_DENSITY){
                    index = search_bc_node(density_total,
                                           densities[ii1].array[ii2].node);
                    if(index < 0) return(SEEK_ERROR_RC);
                    sum[index] += factors[ii1]
                                * densities[ii1].array[ii2].s[ii3];
                    break;
                }
            }
        }
    }

    for(ii1=0; ii1<density_total.size; ii1++){
        for(ii2=0; ii2<BC_SCALAR; ii2++){
            if(density_total.array[ii1].s_type[ii2] == BC_DENSITY){
                density_total.array[ii1].s[ii2] += sum[ii1];
            }
        }
    }

    RC_TRY( mm_free((void *)sum) );

    return(NORMAL_RC);
}


/* ターゲットとなる患者モデルの表面を抽出 */
RC
extract_ref_point(NODE_ARRAY node, NODE_ARRAY *target_node,
		          ELEMENT_ARRAY elem, ELEMENT_ARRAY *target_surf)
{
	int ii1, ii2;
	int label;
	int index1, index2;
	
	if(elem.size < 0) return(ARG_ERROR_RC);

	RC_TRY( extract_surface(elem, target_surf) );
	RC_TRY( allocate_node_array(INIT_SIZE, target_node) );
	for(ii1=0; ii1<target_surf->size; ii1++){
		
		if(target_surf->array[ii1].label < 0) continue;

		for(ii2=0; ii2 <target_surf->array[ii1].node_num; ii2++){
			if(target_surf->array[ii1].node[ii2] < 0) continue;

			label = target_surf->array[ii1].node[ii2];
			index1 = search_node_label(node, label);

			if(search_node_label(*target_node, label) >= 0) continue;

			(target_node->size)++;

			RC_TRY( realloc_node_array(target_node) );

			target_node->array[(target_node->size) - 1].label = label;

			target_node->array[(target_node->size) - 1].p.x
				             = node.array[index1].p.x;
			target_node->array[(target_node->size) - 1].p.y
				             = node.array[index1].p.y;
			target_node->array[(target_node->size) - 1].p.z
				             = node.array[index1].p.z;

			target_node->source = node.source;

			target_node->array[(target_node->size) - 1].i_info[0]
				                             = node.array[index1].i_info[0];
			target_node->array[(target_node->size) - 1].i_info[1]
				                             = node.array[index1].i_info[1];
			target_node->array[(target_node->size) - 1].i_info[2]
				                             = node.array[index1].i_info[2];
			target_node->array[(target_node->size) - 1].i_info[3]
				                             = node.array[index1].i_info[3];
		}
	}

	return(NORMAL_RC);
}



/* armijoの規準 */
int
chk_armijo (double old_obj, double new_obj, double G_rho, double delta_s,
            double *divisor, PARAM param)
{
    double d_L = delta_s * G_rho;

    double delta_L = new_obj - old_obj;

    RC_TRY( log_printf(5, "d_L = %15.7e, delta_L = %15.7e\n", d_L, delta_L) );
    RC_TRY( log_printf(5, "delta_L/d_L = %15.7e\n", delta_L/d_L) );

    if( -delta_L < ABS_ERROR){
        if((*divisor) > MAX_DIVISOR){
            RC_TRY( log_printf(5, "converged: objective function\n") );
            return(-2);
        }
        (*divisor) *= 2.0;
        RC_TRY( log_printf(5, "divisor: %f\n", *divisor) );
        return(0);
    }else if( -delta_L < (param.obj_rel_error*fabs(new_obj) + ABS_ERROR)){
        RC_TRY( log_printf(5, "converged: objective function\n") );
        return(-2);
    }else if(nearly_eq(d_L, 0.0)){
        RC_TRY( log_printf(5, "converged: d_L\n") );
        return(-2);
    }else if(delta_L/d_L < ARMIJO){
        (*divisor) *= 2.0;
        RC_TRY( log_printf(5, "divisor: %f\n", (*divisor)) );
        if((*divisor) > MAX_DIVISOR){
            RC_TRY( log_printf(5, "converged: divisor\n") );
            return(-2);
        }
        return(0);
    }

    return(1);
}   



/*--------------- Input & Output ----------------*/
/* bdfファイルからの読み込み */
RC input_model(FILE *fp, FILE *fp1, NODE_ARRAY *node, ELEMENT_ARRAY *element,
               MATERIAL_PROP_ARRAY *material,
               PHYSICAL_PROP_ARRAY *physical,
               BC_ARRAY *rest)
{
    RC_TRY( nst_input_node(fp1, node) );
    printf("Number of nodes : %d\n", node->size);

    RC_TRY( nst_input_element(fp1, element) );
    printf("Number of elements : %d\n", element->size);

    RC_TRY( nst_input_material_prop(fp, material) );
    printf("Number of materials : %d\n", material->size);

    RC_TRY( nst_input_physical_prop(fp, physical) );
    printf("Number of physical properties : %d\n", physical->size);

    RC_TRY( nst_input_restraint(fp, rest) );
    printf("Number of restraints : %d\n", rest->size);

    return(NORMAL_RC);

}

/* PARAMファイルの読み込み */
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
            fprintf(stdout, "IT_MAX           : %3d\n", param->it_max);
        }else if( !strncmp("SOL_SS_FLAG", buf, 11) ){
            RC_TRY( sgetidc(buf, "ci", idc) );
            param->sol_ss_flag = idc[1].i;
            fprintf(stdout, "SOL_SS_FLAG       : %3d\n", param->sol_ss_flag);
        }else if( !strncmp("DENS_POW", buf, 7) ){
            RC_TRY( sgetidc(buf, "cd", idc) );
            param->dens_pow = idc[1].d;
            fprintf(stdout, "DENS_POW          : %15.7e\n", param->dens_pow);
        }else if( !strncmp("C0", buf, 2) ){
            RC_TRY( sgetidc(buf, "cd", idc) );
            param->c0 = idc[1].d;
            fprintf(stdout, "C0                : %15.7e\n", param->c0);
        }else if( !strncmp("C1", buf, 2) ){
            RC_TRY( sgetidc(buf, "cd", idc) );
            param->c1 = idc[1].d;
            fprintf(stdout, "C1                : %15.7e\n", param->c1);
        }else if( !strncmp("MAX_DENS", buf, 7) ){
            RC_TRY( sgetidc(buf, "cd", idc) );
            param->max_dens = idc[1].d;
            fprintf(stdout, "MAX_DENS          : %15.7e\n", param->max_dens);
        }else if( !strncmp("OBJ_REL_ERROR", buf, 13) ){
            RC_TRY( sgetidc(buf, "cd", idc) );
            param->obj_rel_error = idc[1].d;
            fprintf(stdout, "OBJ_REL_ERROR    : %15.7e\n\n",
                    param->obj_rel_error);
        }
    }

    if(param->it_max <= 0){
        fprintf(stdout, "Invalid data > IT_MAX\n");
        return(ARG_ERROR_RC);
    }
    if( (param->sol_ss_flag != 1) && (param->sol_ss_flag != 2)
         && (param->sol_ss_flag != 4) ){
        fprintf(stdout, "Invalid data > SOL_SS_FLAG\n");
        return(ARG_ERROR_RC);
    }
    if(param->dens_pow <= 0.0){
        fprintf(stdout, "Invalid data > DENS_POW\n");
        return(ARG_ERROR_RC);
    }
    if(param->c0 <= 0.0){
        fprintf(stdout, "Invalid data > C0\n");
        return(ARG_ERROR_RC);
    }
    if(param->c1 <= 0.0){
        fprintf(stdout, "Invalid data > c1\n");
        return(ARG_ERROR_RC);
    }
    if(param->max_dens <= 0.0){
        fprintf(stdout, "Invalid data > MAX_DENS\n");
        return(ARG_ERROR_RC);
    }
    if(param->obj_rel_error <= 0.0){
        fprintf(stdout, "Invalid data > OBJ_REL_ERROR\n");
        return(ARG_ERROR_RC);
    }

    return(NORMAL_RC);
}

static RC
avs_printfe_model (FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY elem)
{
    int ii1;

    fprintf( fp, "%d %d\n", count_valid_node(node), count_valid_element(elem) );

    if( node.size <= 0 ) return(WRITE_ERROR_RC);
    for(ii1=0; ii1<node.size; ii1++){
        if(node.array[ii1].label < 0) continue;

        fprintf(fp, "%10d %15.7E %15.7E %15.7E\n", node.array[ii1].label,
                node.array[ii1].p.x, node.array[ii1].p.y,
                node.array[ii1].p.z);
    }

    if( elem.size <= 0 ) return(WRITE_ERROR_RC);
    for(ii1=0; ii1<elem.size; ii1++){
        if(elem.array[ii1].label < 0) continue;

        RC_TRY( avs_printelement(fp, elem.array[ii1].label,
                    elem.array[ii1].type,
                    elem.array[ii1].node_num,
                    elem.array[ii1].node) );
    }

    return(NORMAL_RC);
}

/* 各ステップの密度の値の出力 */
RC
output_avs_step_bc_density (FILE *fp_avs, NODE_ARRAY node, BC_ARRAY density)
{
    int step_num;

    RC_NEG_ZERO_CHK( step_num = check_step(fp_avs) );
    if(step_num != 1) fprintf(fp_avs, "step%d\n", step_num);

    RC_TRY( avs_printbc_density(fp_avs, density, node) );

    fflush(fp_avs);

    return(NORMAL_RC);
}

/* thetaの値の出力 */
RC
output_avs_bc_theta (FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY elem,
                       BC_ARRAY density)
{
    fprintf(fp, "1\n");
    fprintf(fp, "data_geom\n");
    fprintf(fp, "step1\n");
    RC_TRY( avs_printfe_model(fp, node, elem) );
    RC_TRY( avs_printbc_theta(fp, density, node) );

    fflush(fp);

    return(NORMAL_RC);
}


static RC
avs_printbc_theta (FILE *fp, BC_ARRAY density, NODE_ARRAY node)
{
    int ii1, ii2;
    int index;
    double density_val = 0.0;

    fprintf(fp, "1 0\n");

    fprintf(fp, "1 1\n");
    fprintf(fp, "density,\n");

    for(ii1=0; ii1<node.size; ii1++){
        if(node.array[ii1].label < 0) continue;

        index = search_bc_node(density, node.array[ii1].label);
        if(index >= 0){
            for(ii2=0; ii2<BC_SCALAR; ii2++){
                if(density.array[index].s_type[ii2] == BC_DENSITY){
                    density_val = density.array[index].s[ii2];
                    break;
                }
            }
        }

        fprintf(fp, "%10d %15.7e\n", node.array[ii1].label, density_val);
    }

    return(NORMAL_RC);
}

/* 密度の値の出力 */
RC
output_avs_bc_density (FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY elem,
                       BC_ARRAY density)
{
    fprintf(fp, "1\n");
    fprintf(fp, "data_geom\n");
    fprintf(fp, "step1\n");
    RC_TRY( avs_printfe_model(fp, node, elem) );
    RC_TRY( avs_printbc_density(fp, density, node) );

    fflush(fp);

    return(NORMAL_RC);
}


static RC
avs_printbc_density (FILE *fp, BC_ARRAY density, NODE_ARRAY node)
{
    int ii1, ii2;
    int index;
    double density_val = 0.0;

    fprintf(fp, "1 0\n");

    fprintf(fp, "1 1\n");
    fprintf(fp, "density,\n");

    for(ii1=0; ii1<node.size; ii1++){
        if(node.array[ii1].label < 0) continue;

        index = search_bc_node(density, node.array[ii1].label);
        if(index >= 0){
            for(ii2=0; ii2<BC_SCALAR; ii2++){
                if(density.array[index].s_type[ii2] == BC_DENSITY){
                    density_val = (tanh(density.array[index].s[ii2]) + 1.0)/2.0;
                    break;
                }
            }
        }

        fprintf(fp, "%10d %15.7e\n", node.array[ii1].label, density_val);
    }

    return(NORMAL_RC);
}

/* 収束履歴の出力 */
RC
output_log (FILE *stream, int iteration, double obj0, double obj)
{
    fprintf(stream, "%3d : obj %15.7e ( %6.2f %%)\n",
            iteration, obj, 100.0*obj/obj0);
    fprintf(stdout, "%3d : obj %15.7e ( %6.2f %%)\n",
            iteration, obj, 100.0*obj/obj0);

    return(NORMAL_RC);
}

/* ひずみの出力 */
static RC
output_avs_strain (FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY elem,
                   STRAIN_ARRAY strain)
{
    if(elem.size <= 0) return(WRITE_ERROR_RC);

    fprintf(fp,"1\n");
    fprintf(fp,"data\n");
    fprintf(fp,"step1\n");

    RC_TRY( avs_printfe_model4discont_data(fp, node, elem) );
    RC_TRY( avs_printstrain(fp, elem, strain) );

    fflush(fp);

    return(NORMAL_RC);
}

static RC
avs_printstrain (FILE *fp, ELEMENT_ARRAY elem, STRAIN_ARRAY strain)
{
    int ii1, ii2;
    int counter;
    int index;
    double pr_val[3];

    fprintf(fp, "8 0\n");
    fprintf(fp, "8 1 1 1 1 1 1 1 1\n");
    fprintf(fp, "mises_strain,\n");
    fprintf(fp, "volume_strain,\n");
    fprintf(fp, "strain_x,\n");
    fprintf(fp, "strain_y,\n");
    fprintf(fp, "strain_z,\n");
    fprintf(fp, "strain_yz,\n");
    fprintf(fp, "strain_zx,\n");
    fprintf(fp, "strain_xy,\n");

    counter = 0;
    for(ii1=0; ii1<elem.size; ii1++){
        if(elem.array[ii1].label < 0) continue;
        for(ii2=0; ii2<elem.array[ii1].node_num; ii2++){
            counter ++;
            index = search_strain_element_node(strain,
                    elem.array[ii1].label,
                    elem.array[ii1].node[ii2]);
            RC_NEG_CHK(index);
            RC_TRY( principal3d(strain.array[index].v, pr_val, NULL, NULL,
                                NULL) );
            fprintf(fp, "%10d %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e"
                    "%15.7e %15.7e\n",
                    counter,
                    mises_strain(strain.array[index]),
                    volume_strain(strain.array[index]),
                    strain.array[index].v.x,
                    strain.array[index].v.y,
                    strain.array[index].v.z,
                    strain.array[index].v.yz,
                    strain.array[index].v.zx,
                    strain.array[index].v.xy);
        }
    }
    return(NORMAL_RC);
}

static RC
avs_printfe_model4discont_data (FILE *fp, NODE_ARRAY node, ELEMENT_ARRAY elem)
{
    int ii1, ii2;
    int index;
    int counter = 0;
    int new_node_size = 0;
    int elem_node[MAX_NODE];

    for(ii1=0; ii1<elem.size; ii1++){
        if(elem.array[ii1].label < 0) continue;
        new_node_size += elem.array[ii1].node_num;
    }
    fprintf(fp,"%d %d\n", new_node_size, count_valid_element(elem));

    for(ii1=0; ii1<elem.size; ii1++){
        if(elem.array[ii1].label < 0) continue;
        for(ii2=0; ii2<elem.array[ii1].node_num; ii2++){
            counter ++;
            index = search_node_label(node, elem.array[ii1].node[ii2]);
            RC_NEG_CHK(index);
            fprintf(fp, "%10d %15.7e %15.7e %15.7e\n", counter,
                    node.array[index].p.x,
                    node.array[index].p.y,
                    node.array[index].p.z);
        }
    }

    counter = 0;
    for(ii1=0; ii1<elem.size; ii1++){
        if(elem.array[ii1].label < 0) continue;

        for(ii2=0; ii2<elem.array[ii1].node_num; ii2++){
            counter ++;
            elem_node[ii2] = counter;
        }
        RC_TRY( avs_printelement(fp, elem.array[ii1].label,
                                 elem.array[ii1].type,
                                 elem.array[ii1].node_num,
                                 elem_node) );
    }

    return(NORMAL_RC);
}

static RC
avs_printelement (FILE *fp, int label, ELEM_TYPE type,
        int node_num, const int *node)
{
    int ii1;

    fprintf(fp, "%10d %d ", label, 1);

    RC_TRY( avs_printelem_type(fp, type) );

    switch(type){
        case ELEM_TRI2:
            fprintf(fp, "%10d %10d %10d %10d %10d %10d\n",
                    node[0], node[1], node[2], node[5], node[3], node[4]);
            break;
        case ELEM_TETRA2:
            fprintf(fp, "%10d %10d %10d %10d %10d %10d %10d %10d %10d %10d\n",
                    node[0], node[1], node[2], node[3], node[6],
                    node[5], node[7], node[4], node[9], node[8]);
            break;
        case ELEM_PENTA2:
            fprintf(fp, "%10d %10d %10d %10d %10d %10d %10d %10d %10d %10d"
                    "%10d %10d %10d %10d %10d\n",
                    node[0], node[1],  node[2],  node[3],  node[4],
                    node[5], node[8],  node[6],  node[7],  node[11],
                    node[9], node[10], node[12], node[13], node[14]);
            break;
        default:
            for(ii1=0; ii1<node_num; ii1++){
                fprintf(fp, "%10d ", node[ii1]);
            }
            fprintf(fp,"\n");
    }

    return(NORMAL_RC);
}

static RC
avs_printelem_type (FILE *fp, ELEM_TYPE type)
{
    switch(type){
        case ELEM_LINE1:
            fprintf(fp, "line   ");
            break;
        case ELEM_LINE2:
            fprintf(fp, "line2  ");
            break;
        case ELEM_TRI1:
            fprintf(fp, "tri    ");
            break;
        case ELEM_TRI2:
            fprintf(fp, "tri2   ");
            break;
        case ELEM_QUAD1:
            fprintf(fp, "quad   ");
            break;
        case ELEM_QUAD2:
            fprintf(fp, "quad2  ");
            break;
        case ELEM_TETRA1:
            fprintf(fp, "tet    ");
            break;
        case ELEM_TETRA2:
            fprintf(fp, "tet2   ");
            break;
        case ELEM_PENTA1:
            fprintf(fp, "prism  ");
            break;
        case ELEM_PENTA2:
            fprintf(fp, "prism2 ");
            break;
        case ELEM_HEXA1:
            fprintf(fp, "hex    ");
            break;
        case ELEM_HEXA2:
            fprintf(fp, "hex2   ");
            break;
        default:
            return(IMPLEMENT_ERROR_RC);
    }

    return(NORMAL_RC);
}

/* ミーゼスひずみの計算 */
static double
mises_strain (STRESS_STRAIN strain)
{
    double ret;

    ret = 0.5*((strain.v.x - strain.v.y)*(strain.v.x - strain.v.y)
        + (strain.v.y - strain.v.z)*(strain.v.y - strain.v.z)
        + (strain.v.z - strain.v.x)*(strain.v.z - strain.v.x))
        + 3.0*((strain.v.yz * strain.v.yz)
        +(strain.v.zx * strain.v.zx)
        +(strain.v.xy * strain.v.xy));
    ret = sqrt(1.0/9.0*ret+ABS_TOL);

    return(ret);
}

/* 体積ひずみの計算 */
static double
volume_strain (STRESS_STRAIN strain)
{
    double ret;

    ret = strain.v.x + strain.v.y + strain.v.z;

    return(ret);
}

static int
check_step (FILE *fp_avs)
{
    IDC idc[1];

    fflush(fp_avs);

    if( fseek(fp_avs, 0L, SEEK_SET) < 0) return(-1);

    /* 現在までのステップ数を調べる */
    if(fgetidc(fp_avs, "i", idc ) != NORMAL_RC) return(-1);

    if( fseek(fp_avs, 0L, SEEK_SET) < 0) return(-1);

    /* 現在のステップ数 = 過去のステップ数 + 1 */
    (idc[0].i)++;
    fprintf(fp_avs, "%d", idc[0].i);
    fflush(fp_avs);
    if( fseek(fp_avs, 0L, SEEK_END) < 0) return(-1);

    return(idc[0].i);
}

