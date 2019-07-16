/*********************************************************************
 * fem_solver_6dof.c
 *
 * Copyright (C) 2002 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Thanks to following contributors.
 *  <Yasuo ASAGA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: fem_solver_6dof.c,v 1.15 2003/07/22 11:44:01 adachi Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include "fem_solver.h"
#include "fem_struct.h"
#include "math_utl.h"
#include "rc.h"
#include "nst_component.h"
#include "nonzero_cg.h"
#include "sky_cholesky.h"


static RC fill_row_nodes6(FEM_ELEM_LOOKUP elem_lookup,
                          FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                          int row_nodes1[], int *row_size1,
                          int row_nodes2[], int *row_size2,
                          int cache_size, int cache[][2]);
static RC rigid_modify_nonzero6(NONZERO_MATRIX3 matrix,
               FEM_NODE_ARRAY node, FEM_RIGID_ELEMENT_ARRAY rigid,
               int slave_flag[], int master_flag[], double master_weight[],
               int *master_count, double **f_vect, int vect_size,
               double weight, int label0, int dof0, int label, int dof,
               const int slave_table[]);
static RC add_master_dof(NONZERO_MATRIX3 matrix, int master_count,
                        int master_flag[], double master_weight[],
                        int index0, int dof0);
static int dof_chk(int dof, int i);
static void init_flags(int flags[], int size, int i);
static RC expand_rigid(FEM_NODE_ARRAY node, FEM_RIGID_ELEMENT_ARRAY src,
                       FEM_RIGID_ELEMENT_ARRAY *dest, int slave_table[]);
static RC expand_rigid_sub(FEM_NODE_ARRAY node, FEM_RIGID_ELEMENT rigid,
                           FEM_RIGID_ELEMENT_ARRAY *dest, int slave_table[]);
static RC cal_index_dof6(FEM_NODE_ARRAY node, int label, int src_dof,
                         int *index, int *dof);
static TRANSLATION3D trans_local2global(FEM_LOCAL_COORD coord,
                                        TRANSLATION3D local);
static void local_coord2s_matrix(FEM_LOCAL_COORD coord, double s_matrix[][3]);
static void local_coord2s_matrix_t(FEM_LOCAL_COORD coord, double s_matrix[][3]);
static void mul_matrix33_vect(double s_matrix[][3], double vect[]);
static RC rigid_recover_sub(int gdim, FEM_NODE_ARRAY node, int slave_flag[],
                            FEM_RIGID_ELEMENT_ARRAY rigid, double weight,
                            int index0, int dof0, int index1, int dof1,
                            double d_vect[], const int slave_table[]);
static RC fem_rigid_recover(int gdim, FEM_NODE_ARRAY node,
                            FEM_RIGID_ELEMENT_ARRAY rigid, double d_vect[],
                            const int slave_table[]);



/* 各節点６自由度の線形弾性静解析ソルバー */
/* ただし、剛性マトリックスは 3*3/ブロック */
RC fem_solve_iccg6(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                   FEM_RIGID_ELEMENT_ARRAY rigid,
                   FEM_MATERIAL_PROP_ARRAY material,
                   FEM_PHYSICAL_PROP_ARRAY physical,
                   FEM_LOCAL_COORD_ARRAY coord,
                   FEM_BC_ARRAY rest, FEM_BC_ARRAY force,
                   FEM_BC_ARRAY temp, FEM_DEFAULT_BC_ARRAY def_temp,
                   FEM_DEFAULT_BC_ARRAY def_body_force,
                   FEM_DISP_ARRAY *disp, FEM_STRESS_ARRAY *stress,
                   FEM_STRAIN_ARRAY *strain, int fem_sol_ss_flag)
{
	int ii1;
	NONZERO_MATRIX3 matrix;
	double *f_vect;
	double *d_vect;
	FEM_ELEM_MATRIX elem_matrix;

	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

	/* 全体剛性マトリックスの作成 */
	RC_TRY( fem_allocate_nonzero6(node, element, rigid, &matrix) );
	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		RC_TRY( fem_make_elem_matrix(element.array[ii1], node,
		                             material, physical, &elem_matrix) );
		RC_TRY( fem_impose_nonzero6(matrix, elem_matrix, 1.0) );
		RC_TRY( free_fem_elem_matrix(&elem_matrix) );
	}

	/* 荷重ベクトルの作成 */
	RC_TRY( fem_allocate_vector(6, node, &f_vect) );
	RC_TRY( fem_add_nodal_force6(node, coord, force, f_vect) );
	RC_TRY( fem_add_temp_force(6, node, element, material, physical,
	                           temp, def_temp, f_vect) );
	RC_TRY( fem_add_body_force(6, node, element, material, physical,
	                           def_body_force, f_vect) );
	RC_TRY( fem_force_local_pre6(node, coord, f_vect) );

	/* 変位の計算 */
	RC_TRY( fem_modify_nonzero6(node, coord, rest, rigid, matrix,
	                            &f_vect, 1, NULL) );
	RC_TRY( fem_scale_iccg_nonzero3(matrix, f_vect, &d_vect, NULL) );
	RC_TRY( free_nonzero_matrix3(&matrix) );
	RC_TRY( fem_free_vector(&f_vect) );
	RC_TRY( fem_rigid_recover6(node, rigid, d_vect) );
	RC_TRY( fem_disp_local_recover6(node, coord, d_vect) );
	RC_TRY( fem_vect2disp(6, node, d_vect, disp) );
	RC_TRY( fem_free_vector(&d_vect) );

	/* ひずみ、応力の計算 */
	RC_TRY( fem_stress_strain(fem_sol_ss_flag, node, element, *disp,
	                          material, physical, stress, strain) );

	return(NORMAL_RC);
}


/* 各節点６自由度の線形弾性静解析シェルソルバー */
/* ただし、剛性マトリックスは 3*3/ブロック */
RC fem_solve_shell_iccg6(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                         FEM_RIGID_ELEMENT_ARRAY rigid,
                         FEM_MATERIAL_PROP_ARRAY material,
                         FEM_PHYSICAL_PROP_ARRAY physical,
                         FEM_LOCAL_COORD_ARRAY coord,
                         FEM_BC_ARRAY rest, FEM_BC_ARRAY force,
                         FEM_BC_ARRAY temp, FEM_DEFAULT_BC_ARRAY def_temp,
                         FEM_DEFAULT_BC_ARRAY def_body_force,
                         FEM_ELEM_BC_ARRAY press, double distance[],
                         FEM_DISP_ARRAY *disp, FEM_STRESS_ARRAY stress[],
                         FEM_STRAIN_ARRAY strain[], FEM_STRESS_ARRAY *iforce,
                         FEM_STRESS_ARRAY *moment,
                         FEM_STRESS_ARRAY *shear_force,
                         FEM_STRAIN_ARRAY *curv, int fem_sol_ss_flag)
{
	int ii1;
	double *f_vect;
	double *d_vect;
	NONZERO_MATRIX3 matrix;
	FEM_ELEM_MATRIX elem_matrix;

	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );

	/* 全体剛性マトリックスの作成 */
	RC_TRY( fem_allocate_nonzero6(node, element, rigid, &matrix) );

	for(ii1=0; ii1<element.size; ii1++){
		if(element.array[ii1].label < 0) continue;

		fprintf(stderr, " [%6d/%6d]\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",
		        element.array[ii1].label, element.size);
		RC_TRY( fem_make_shell_matrix(element.array[ii1], node, material,
		                              physical, &elem_matrix) );
		RC_TRY( fem_impose_nonzero6(matrix, elem_matrix, 1.0) );

		RC_TRY( free_fem_elem_matrix(&elem_matrix) );
	}
	fprintf(stderr, "\n");

	/* 荷重ベクトルの作成 */
	RC_TRY( fem_allocate_vector(6, node, &f_vect) );
	RC_TRY( fem_add_nodal_force6(node, coord, force, f_vect) );
/***未対応****
	RC_TRY( fem_add_temp_force(6, node, element, material, physical,
	                           temp, def_temp, f_vect) );
**************/
	RC_TRY( fem_add_body_force(6, node, element, material, physical,
	                           def_body_force, f_vect) );
	RC_TRY( fem_force_local_pre6(node, coord, f_vect) );

	/* 変位の計算 */
	RC_TRY( fem_modify_nonzero6(node, coord, rest, rigid, matrix, &f_vect,
	                           1, NULL) );
	RC_TRY( fem_scale_iccg_nonzero3(matrix, f_vect, &d_vect, NULL) );
	RC_TRY( free_nonzero_matrix3(&matrix) );
	RC_TRY( fem_free_vector(&f_vect) );
	RC_TRY( fem_rigid_recover6(node, rigid, d_vect) );
	RC_TRY( fem_disp_local_recover6(node, coord, d_vect) );
	RC_TRY( fem_vect2disp(6, node, d_vect, disp) );
	RC_TRY( fem_free_vector(&d_vect) );

	/* ひずみ、応力の計算 */
	RC_TRY( fem_stress_strain_shell(fem_sol_ss_flag, distance, node, element,
	                                *disp, material, physical, stress, strain));
	/* 曲率、合応力の計算 */
	RC_TRY( fem_resultant_stress_strain_shell(fem_sol_ss_flag, node, element,
	            *disp, material, physical, iforce, moment, shear_force, curv) );

	return(NORMAL_RC);
}


/* 全体マトリックスの確保(ただし，rigid_elem 分は無視する) */
RC fem_allocate_nonzero6(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                         FEM_RIGID_ELEMENT_ARRAY rigid_elem,
                         NONZERO_MATRIX3 *matrix)
{
	int ii1, ii2;
	FEM_ELEM_LOOKUP_ARRAY elem_lookup;
	int *row_nodes1;
	int *row_nodes2;
	int row_size1, row_size2;
	int syn_size1, syn_size2;
	int index1, index2;
	int node_size;
	const int cache_size = 2048;
	int (*cache)[2];

	if(node.renum_flag != RENUMBERED) RC_TRY( dummy_renumber(&node) );
	node_size = count_valid_node(node);

	/* 作業用配列 */
	RC_NULL_CHK( row_nodes1 = (int *)malloc(2*node_size*sizeof(int)) );
	RC_NULL_CHK( row_nodes2 = (int *)malloc(2*node_size*sizeof(int)) );
	RC_NULL_CHK( cache = (int(*)[2])malloc(2*cache_size*sizeof(int)) );
	for(ii1=0; ii1<cache_size; ii1++){
		cache[ii1][0] = -1;
		cache[ii1][1] = -1;
	}

	RC_TRY( make_elem_lookup_rigid(node, element, rigid_elem, &elem_lookup) );

	matrix->size = 2*node_size;
	matrix->row = (NONZERO_ROW3 *)malloc(matrix->size*sizeof(NONZERO_ROW3));
	RC_NULL_CHK(matrix->row);
	for(ii1=0; ii1<matrix->size; ii1++){
		matrix->row[ii1].size = 0;
		matrix->row[ii1].elem = NULL;
		matrix->row[ii1].syn_size = 0;
		matrix->row[ii1].syn = NULL;
	}

	for(ii1=0; ii1<elem_lookup.size; ii1++){
		if(elem_lookup.array[ii1].node < 0) continue;

		RC_NEG_CHK( index1 = 2 * node_renum_index(node,
			                 elem_lookup.array[ii1].node) );
		index2 = index1 + 1;

		RC_TRY( fill_row_nodes6(elem_lookup.array[ii1], node, element,
		                        row_nodes1, &row_size1,
		                        row_nodes2, &row_size2, cache_size, cache) );
		syn_size1 = 0;
		for(ii2=0; ii2<row_size1; ii2++){
			if(row_nodes1[ii2] > index1) syn_size1++;
		}
		syn_size2 = 0;
		for(ii2=0; ii2<row_size2; ii2++){
			if(row_nodes2[ii2] > index2) syn_size2++;
		}

		matrix->row[index1].syn_size = syn_size1;
		if(syn_size1 > 0){
			RC_NULL_CHK( matrix->row[index1].syn
			             = malloc(syn_size1*sizeof(int)) );
			for(ii2=0; ii2<syn_size1; ii2++){
				matrix->row[index1].syn[ii2]
				    = row_nodes1[row_size1 - syn_size1 + ii2];
			}
		}
		matrix->row[index2].syn_size = syn_size2;
		if(syn_size2 > 0){
			RC_NULL_CHK( matrix->row[index2].syn
			             = malloc(syn_size2*sizeof(int)) );
			for(ii2=0; ii2<syn_size2; ii2++){
				matrix->row[index2].syn[ii2]
				    = row_nodes2[row_size2 - syn_size2 + ii2];
			}
		}

		row_size1 -= syn_size1;
		if(row_size1 <= 0) return(UNKNOWN_ERROR_RC);
		if(row_nodes1[row_size1 - 1] != index1) return(UNKNOWN_ERROR_RC);
		matrix->row[index1].size = row_size1;
		RC_NULL_CHK( matrix->row[index1].elem
		             = malloc(row_size1*sizeof(NONZERO_ELEM3)) );
		for(ii2=0; ii2<row_size1; ii2++){
			matrix->row[index1].elem[ii2].col = row_nodes1[ii2];
			matrix->row[index1].elem[ii2].v[0][0] = 0.0;
			matrix->row[index1].elem[ii2].v[0][1] = 0.0;
			matrix->row[index1].elem[ii2].v[0][2] = 0.0;
			matrix->row[index1].elem[ii2].v[1][0] = 0.0;
			matrix->row[index1].elem[ii2].v[1][1] = 0.0;
			matrix->row[index1].elem[ii2].v[1][2] = 0.0;
			matrix->row[index1].elem[ii2].v[2][0] = 0.0;
			matrix->row[index1].elem[ii2].v[2][1] = 0.0;
			matrix->row[index1].elem[ii2].v[2][2] = 0.0;
		}
		row_size2 -= syn_size2;
		if(row_size2 <= 0) return(UNKNOWN_ERROR_RC);
		if(row_nodes2[row_size2 - 1] != index2) return(UNKNOWN_ERROR_RC);
		matrix->row[index2].size = row_size2;
		RC_NULL_CHK( matrix->row[index2].elem
		             = malloc(row_size2*sizeof(NONZERO_ELEM3)) );
		for(ii2=0; ii2<row_size2; ii2++){
			matrix->row[index2].elem[ii2].col = row_nodes2[ii2];
			matrix->row[index2].elem[ii2].v[0][0] = 0.0;
			matrix->row[index2].elem[ii2].v[0][1] = 0.0;
			matrix->row[index2].elem[ii2].v[0][2] = 0.0;
			matrix->row[index2].elem[ii2].v[1][0] = 0.0;
			matrix->row[index2].elem[ii2].v[1][1] = 0.0;
			matrix->row[index2].elem[ii2].v[1][2] = 0.0;
			matrix->row[index2].elem[ii2].v[2][0] = 0.0;
			matrix->row[index2].elem[ii2].v[2][1] = 0.0;
			matrix->row[index2].elem[ii2].v[2][2] = 0.0;
		}
	}

	free(row_nodes1);
	free(row_nodes2);
	free(cache);
	RC_TRY( free_fem_elem_lookup_array(&elem_lookup) );

	return(NORMAL_RC);
}


/* 該当行に関連する列のインデックスを row_nodes に登録 */
static RC fill_row_nodes6(FEM_ELEM_LOOKUP elem_lookup,
                          FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                          int row_nodes1[], int *row_size1,
                          int row_nodes2[], int *row_size2,
                          int cache_size, int cache[][2])
{
	int ii1, ii2, ii3;
	int flag, col_index;
	int row_index1;
	int elem_index;
	int elem_dim;

	/* 少なくとも対角要素だけは確保する */
	RC_NEG_CHK( row_index1 = 2 * node_renum_index_cache(node, elem_lookup.node,
	                                                    cache_size, cache) );
	*row_size1 = 1;
	*row_size2 = 1;
	row_nodes1[0] = row_index1;
	row_nodes2[0] = row_index1+1;

	for(ii1=0; ii1<elem_lookup.size; ii1++){
		elem_index = search_fem_element_label(element,
		                                      elem_lookup.element[ii1]);
		RC_NEG_CHK(elem_index);
		elem_dim = element_dim(element.array[elem_index].type);

		for(ii2=0; ii2<element.array[elem_index].node_num; ii2++){
			RC_NEG_CHK( col_index = 2 * node_renum_index_cache(node,
			                            element.array[elem_index].node[ii2],
			                            cache_size, cache) );
			flag = 0;
			for(ii3=0; ii3<(*row_size1); ii3++){
				if(row_nodes1[ii3] == col_index){
					flag = 1;
					break;
				}
			}
			if(flag == 0){
				row_nodes1[*row_size1] = col_index;   /* 登録 */
				(*row_size1)++;
			}

			if(elem_dim == 3) continue;

			/* 三次元ソリッド要素以外の処理 */
			flag = 0;
			for(ii3=0; ii3<(*row_size1); ii3++){
				if(row_nodes1[ii3] == col_index+1){
					flag = 1;
					break;
				}
			}
			if(flag == 0){
				row_nodes1[*row_size1] = col_index+1;   /* 登録 */
				(*row_size1)++;
			}

			flag = 0;
			for(ii3=0; ii3<(*row_size2); ii3++){
				if(row_nodes2[ii3] == col_index){
					flag = 1;
					break;
				}
			}
			if(flag == 0){
				row_nodes2[*row_size2] = col_index;   /* 登録 */
				(*row_size2)++;
			}

			flag = 0;
			for(ii3=0; ii3<(*row_size2); ii3++){
				if(row_nodes2[ii3] == col_index+1){
					flag = 1;
					break;
				}
			}
			if(flag == 0){
				row_nodes2[*row_size2] = col_index+1;   /* 登録 */
				(*row_size2)++;
			}
		}
	}

	/* ソートしておく */
	my_sort_int(*row_size1, row_nodes1);
	my_sort_int(*row_size2, row_nodes2);

	return(NORMAL_RC);
}


RC fem_impose_nonzero6(NONZERO_MATRIX3 matrix,
                       FEM_ELEM_MATRIX elem_matrix, double fact)
{
	int ii1, ii2, ii3, ii4;
	int ix1, ix2;
	int col_index;
	int dim = elem_matrix.dim;

	if((dim != 3)&&(dim != 6)) return(ARG_ERROR_RC);

	for(ii1=0; ii1<elem_matrix.size; ii1+=3){
		for(ii2=0; ii2<elem_matrix.size; ii2+=3){
			ix1 = elem_matrix.index[ii1]/3;
			ix2 = elem_matrix.index[ii2]/3;
			if(ix1 < ix2) continue;
			if(dim == 3){
				ix1 *= 2;
				ix2 *= 2;
			}

			if(ix1 >= matrix.size) return(ARG_ERROR_RC);
			col_index = -1;
			for(ii3=0; ii3<matrix.row[ix1].size; ii3++){
				if(matrix.row[ix1].elem[ii3].col == ix2){
					col_index = ii3;
					break;
				}
			}
			RC_NEG_CHK( col_index );

			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					matrix.row[ix1].elem[col_index].v[ii3][ii4]
					      += fact * elem_matrix.matrix[ii1+ii3][ii2+ii4];
				}
			}
		}
	}

	return(NORMAL_RC);
}


static RC rigid_modify_nonzero6(NONZERO_MATRIX3 matrix,
               FEM_NODE_ARRAY node, FEM_RIGID_ELEMENT_ARRAY rigid,
               int slave_flag[], int master_flag[], double master_weight[],
               int *master_count, double **f_vect, int vect_size,
               double weight, int label0, int dof0, int label, int dof,
               const int slave_table[])
{
	int ii1;
	int index, ndof;

	RC_TRY( cal_index_dof6(node, label, dof, &index, &ndof) );

	if(slave_table[3*index+ndof] >= 0){
		int r_index = slave_table[3*index+ndof];

		if(slave_flag[r_index]) return(UNKNOWN_ERROR_RC);
		slave_flag[r_index] = 1;
		for(ii1=1; ii1<rigid.array[r_index].node_num; ii1++){
			int label_new = rigid.array[r_index].node[ii1];
			int dof_new = rigid.array[r_index].dof[ii1];
			double weight_new = rigid.array[r_index].weight[ii1];

			RC_TRY( rigid_modify_nonzero6(matrix, node, rigid, slave_flag,
			                master_flag, master_weight, master_count,
			                f_vect, vect_size, weight*weight_new,
			                label0, dof0, label_new, dof_new, slave_table) );
		}
	}else{
		int index0, ndof0;

		RC_TRY( cal_index_dof6(node, label0, dof0, &index0, &ndof0) );
		master_flag[*master_count] = 3*index+ndof;
		master_weight[*master_count] = weight;
		(*master_count)++;
		RC_TRY( row_col_add_nonzero3s(matrix, index0, ndof0,
		                              weight, index, ndof) );
		for(ii1=0; ii1<vect_size; ii1++){
			f_vect[ii1][3*index+ndof] += weight * f_vect[ii1][3*index0+ndof0];
		}
	}

	return(NORMAL_RC);
}


static RC expand_rigid(FEM_NODE_ARRAY node, FEM_RIGID_ELEMENT_ARRAY src,
                       FEM_RIGID_ELEMENT_ARRAY *dest, int slave_table[])
{
	int ii1;
	int index, dof;

	for(ii1=0; ii1<6*node.size; ii1++){
		slave_table[ii1] = -1;
	}

	RC_TRY( allocate_fem_rigid_element_array(src.size, dest) );

	for(ii1=0; ii1<src.size; ii1++){
		if((src.array[ii1].type == ELEM_MPC)&&(src.array[ii1].node_num > 1)){
			RC_TRY( realloc_fem_rigid_element_array(dest) );
			RC_TRY( copy_fem_rigid_element(src.array[ii1],
			                               &(dest->array[dest->size])) );
			dest->array[dest->size].label = (dest->size) + 1;

			RC_TRY( cal_index_dof6(node, src.array[ii1].node[0],
			                       src.array[ii1].dof[0], &index, &dof) );
			slave_table[3*index+dof] = dest->size;

			(dest->size)++;
		}else if(src.array[ii1].type == ELEM_RBAR){
			RC_TRY( expand_rigid_sub(node, src.array[ii1], dest, slave_table) );
		}else if(src.array[ii1].type == ELEM_SPRING){
			continue;
		}else if(src.array[ii1].type == ELEM_MASS){
			continue;
		}else{
			return(IMPLEMENT_ERROR_RC);
		}
	}

	RC_TRY( clean_fem_rigid_element_array(dest) );

	return(NORMAL_RC);
}


static RC expand_rigid_sub(FEM_NODE_ARRAY node, FEM_RIGID_ELEMENT rigid,
                           FEM_RIGID_ELEMENT_ARRAY *dest, int slave_table[])
{
	int ii1, ii2;
	TRANSLATION3D p0, p;
	int n_index;
	int index, dof;

	if(rigid.type != ELEM_RBAR) return(ARG_ERROR_RC);

	for(ii1=1; ii1<rigid.node_num; ii1++){
		RC_NEG_CHK( n_index = search_fem_node_label(node, rigid.node[0]) );
		p0 = node.array[n_index].p;

		RC_NEG_CHK( n_index = search_fem_node_label(node, rigid.node[ii1]) );
		p = node.array[n_index].p;

		if(dof_chk(rigid.dof[0], 1)){
			RC_TRY( realloc_fem_rigid_element_array(dest) );
			dest->array[dest->size].label = (dest->size) + 1;
			dest->array[dest->size].set_id = rigid.set_id;
			dest->array[dest->size].type = ELEM_MPC;
			dest->array[dest->size].node_num = 4;
			dest->array[dest->size].weight_num = 4;
			dest->array[dest->size].dof_num = 4;
			RC_NULL_CHK( dest->array[dest->size].node
			                 = malloc(4 * sizeof(int)) );
			RC_NULL_CHK( dest->array[dest->size].weight
			                 = malloc(4 * sizeof(double)) );
			RC_NULL_CHK( dest->array[dest->size].dof
			                 = malloc(4 * sizeof(int)) );
			dest->array[dest->size].node[0] = rigid.node[ii1];
			dest->array[dest->size].node[1] = rigid.node[0];
			dest->array[dest->size].node[2] = rigid.node[0];
			dest->array[dest->size].node[3] = rigid.node[0];
			dest->array[dest->size].dof[0] = 0;
			dest->array[dest->size].dof[1] = 0;
			dest->array[dest->size].dof[2] = 5;
			dest->array[dest->size].dof[3] = 4;
			dest->array[dest->size].weight[0] = -1.0;
			dest->array[dest->size].weight[1] = 1.0;
			dest->array[dest->size].weight[2] = -(p.y - p0.y);
			dest->array[dest->size].weight[3] = (p.z - p0.z);
			RC_TRY( cal_index_dof6(node, dest->array[dest->size].node[0],
			                  dest->array[dest->size].dof[0], &index, &dof) );
			slave_table[3*index+dof] = dest->size;
			(dest->size)++;
		}
		if(dof_chk(rigid.dof[0], 2)){
			RC_TRY( realloc_fem_rigid_element_array(dest) );
			dest->array[dest->size].label = (dest->size) + 1;
			dest->array[dest->size].set_id = rigid.set_id;
			dest->array[dest->size].type = ELEM_MPC;
			dest->array[dest->size].node_num = 4;
			dest->array[dest->size].weight_num = 4;
			dest->array[dest->size].dof_num = 4;
			RC_NULL_CHK( dest->array[dest->size].node
			                 = malloc(4 * sizeof(int)) );
			RC_NULL_CHK( dest->array[dest->size].weight
			                 = malloc(4 * sizeof(double)) );
			RC_NULL_CHK( dest->array[dest->size].dof
			                 = malloc(4 * sizeof(int)) );
			dest->array[dest->size].node[0] = rigid.node[ii1];
			dest->array[dest->size].node[1] = rigid.node[0];
			dest->array[dest->size].node[2] = rigid.node[0];
			dest->array[dest->size].node[3] = rigid.node[0];
			dest->array[dest->size].dof[0] = 1;
			dest->array[dest->size].dof[1] = 1;
			dest->array[dest->size].dof[2] = 3;
			dest->array[dest->size].dof[3] = 5;
			dest->array[dest->size].weight[0] = -1.0;
			dest->array[dest->size].weight[1] = 1.0;
			dest->array[dest->size].weight[2] = -(p.z - p0.z);
			dest->array[dest->size].weight[3] = (p.x - p0.x);
			RC_TRY( cal_index_dof6(node, dest->array[dest->size].node[0],
			                  dest->array[dest->size].dof[0], &index, &dof) );
			slave_table[3*index+dof] = dest->size;
			(dest->size)++;
		}
		if(dof_chk(rigid.dof[0], 3)){
			RC_TRY( realloc_fem_rigid_element_array(dest) );
			dest->array[dest->size].label = (dest->size) + 1;
			dest->array[dest->size].set_id = rigid.set_id;
			dest->array[dest->size].type = ELEM_MPC;
			dest->array[dest->size].node_num = 4;
			dest->array[dest->size].weight_num = 4;
			dest->array[dest->size].dof_num = 4;
			RC_NULL_CHK( dest->array[dest->size].node
			                 = malloc(4 * sizeof(int)) );
			RC_NULL_CHK( dest->array[dest->size].weight
			                 = malloc(4 * sizeof(double)) );
			RC_NULL_CHK( dest->array[dest->size].dof
			                 = malloc(4 * sizeof(int)) );
			dest->array[dest->size].node[0] = rigid.node[ii1];
			dest->array[dest->size].node[1] = rigid.node[0];
			dest->array[dest->size].node[2] = rigid.node[0];
			dest->array[dest->size].node[3] = rigid.node[0];
			dest->array[dest->size].dof[0] = 2;
			dest->array[dest->size].dof[1] = 2;
			dest->array[dest->size].dof[2] = 4;
			dest->array[dest->size].dof[3] = 3;
			dest->array[dest->size].weight[0] = -1.0;
			dest->array[dest->size].weight[1] = 1.0;
			dest->array[dest->size].weight[2] = -(p.x - p0.x);
			dest->array[dest->size].weight[3] = (p.y - p0.y);
			RC_TRY( cal_index_dof6(node, dest->array[dest->size].node[0],
			                  dest->array[dest->size].dof[0], &index, &dof) );
			slave_table[3*index+dof] = dest->size;
			(dest->size)++;
		}
		for(ii2=3; ii2<6; ii2++){
			if(dof_chk(rigid.dof[0], 1 + ii2)){
				RC_TRY( realloc_fem_rigid_element_array(dest) );
				dest->array[dest->size].label = (dest->size) + 1;
				dest->array[dest->size].set_id = rigid.set_id;
				dest->array[dest->size].type = ELEM_MPC;
				dest->array[dest->size].node_num = 2;
				dest->array[dest->size].weight_num = 2;
				dest->array[dest->size].dof_num = 2;
				RC_NULL_CHK( dest->array[dest->size].node
				                 = malloc(2 * sizeof(int)) );
				RC_NULL_CHK( dest->array[dest->size].weight
				                 = malloc(2 * sizeof(double)) );
				RC_NULL_CHK( dest->array[dest->size].dof
				                 = malloc(2 * sizeof(int)) );
				dest->array[dest->size].node[0] = rigid.node[ii1];
				dest->array[dest->size].node[1] = rigid.node[0];
				dest->array[dest->size].dof[0] = ii2;
				dest->array[dest->size].dof[1] = ii2;
				dest->array[dest->size].weight[0] = -1.0;
				dest->array[dest->size].weight[1] = 1.0;
				RC_TRY( cal_index_dof6(node, dest->array[dest->size].node[0],
				            dest->array[dest->size].dof[0], &index, &dof) );
				slave_table[3*index+dof] = dest->size;
				(dest->size)++;
			}
		}
	}

	return(NORMAL_RC);
}


static TRANSLATION3D trans_local2global(FEM_LOCAL_COORD coord,
                                        TRANSLATION3D local)
{
	double s_matrix[3][3];
	TRANSLATION3D ret;

	local_coord2s_matrix(coord, s_matrix);
	ret.x = s_matrix[0][0]*local.x
	      + s_matrix[0][1]*local.y
	      + s_matrix[0][2]*local.z;
	ret.y = s_matrix[1][0]*local.x
	      + s_matrix[1][1]*local.y
	      + s_matrix[1][2]*local.z;
	ret.z = s_matrix[2][0]*local.x
	      + s_matrix[2][1]*local.y
	      + s_matrix[2][2]*local.z;

	return(ret);
}


static void local_coord2s_matrix(FEM_LOCAL_COORD coord, double s_matrix[][3])
{
	s_matrix[0][0] = coord.d_cos[0].x;
	s_matrix[0][1] = coord.d_cos[1].x;
	s_matrix[0][2] = coord.d_cos[2].x;
	s_matrix[1][0] = coord.d_cos[0].y;
	s_matrix[1][1] = coord.d_cos[1].y;
	s_matrix[1][2] = coord.d_cos[2].y;
	s_matrix[2][0] = coord.d_cos[0].z;
	s_matrix[2][1] = coord.d_cos[1].z;
	s_matrix[2][2] = coord.d_cos[2].z;
}


static void local_coord2s_matrix_t(FEM_LOCAL_COORD coord, double s_matrix[][3])
{
	s_matrix[0][0] = coord.d_cos[0].x;
	s_matrix[0][1] = coord.d_cos[0].y;
	s_matrix[0][2] = coord.d_cos[0].z;
	s_matrix[1][0] = coord.d_cos[1].x;
	s_matrix[1][1] = coord.d_cos[1].y;
	s_matrix[1][2] = coord.d_cos[1].z;
	s_matrix[2][0] = coord.d_cos[2].x;
	s_matrix[2][1] = coord.d_cos[2].y;
	s_matrix[2][2] = coord.d_cos[2].z;
}


#if 0
RC fem_make_beam_mass_matrix(FEM_ELEMENT element, FEM_NODE_ARRAY node,
                             FEM_MATERIAL_PROP_ARRAY material,
                             FEM_PHYSICAL_PROP_ARRAY physical,
                             FEM_ELEM_MATRIX *elem_matrix)
{
	int ii1, ii2, ii3, ii4;
	double rho;
	double A, Ixx, Iyy, J;
	double len, len2;
	int p_index, m_index, r_index, n_index;
	TRANSLATION3D p1, p2, pd;
	TRANSLATION3D v1, v2, vtmp;
	double T[3][3];
	double Tt[3][3];
	double Ttt[3][3];

	if(element.label < 0) return(ARG_ERROR_RC);
	if(element.type != ELEM_BEAM1) return(ARG_ERROR_RC);
	elem_matrix->element = element.label;
	elem_matrix->size = 12;
	elem_matrix->dim = 6;

	m_index = search_fem_material_prop_element(material, physical, &element);
	RC_NEG_CHK( m_index );
	if(material.array[m_index].mat_type != MAT_ISOTROPIC) return(ARG_ERROR_RC);
	rho = material.array[m_index].rho;

	p_index = search_fem_physical_prop_label(physical, element.physical);
	RC_NEG_CHK( p_index );
	A   = physical.array[p_index].d_info[0];
	Ixx = physical.array[p_index].d_info[1];
	Iyy = physical.array[p_index].d_info[2];
	J   = physical.array[p_index].d_info[3];
	v1.x= element.d_info[0];
	v1.y= element.d_info[1];
	v1.z= element.d_info[2];

	/* elem_matrix->array[][] の確保 */
	RC_TRY( allocate2D(elem_matrix->size,
	                   elem_matrix->size, &(elem_matrix->matrix)) );

	/* elem_matirx->index[] の確保と代入 */
	RC_TRY( rc_malloc(sizeof(int)*(elem_matrix->size),
	                  (void**)(&elem_matrix->index)) );
	for(ii1=0; ii1<element.node_num; ii1++){
		RC_NEG_CHK( r_index = node_renum_index(node, element.node[ii1]) );
		for(ii2=0; ii2<(elem_matrix->dim); ii2++){
			elem_matrix->index[ii1*(elem_matrix->dim)+ii2]
			             = r_index*(elem_matrix->dim) + ii2;
		}
	}

	RC_NEG_CHK( n_index = search_fem_node_label(node, element.node[0]) );
	p1 = node.array[n_index].p;

	RC_NEG_CHK( n_index = search_fem_node_label(node, element.node[1]) );
	p2 = node.array[n_index].p;

	len = abs_translation3d( pd = sub_translation3d(p2, p1) );
	len2 = len * len;
	pd = mul_scalar_translation3d(1.0/len, pd);
	vtmp = mul_scalar_translation3d(inner_product3d(v1, pd), pd);
	v1 = sub_translation3d(v1, vtmp);
	v1 = mul_scalar_translation3d(1.0/abs_translation3d(v1), v1);
	v2 = outer_product3d(pd, v1);

	elem_matrix->matrix[0][0] = 13.0/35.0 + 6*Iyy/(5.0*A*len2);
	elem_matrix->matrix[0][4] = -(-11.0*len/210.0 + Iyy/(10.0*A*len));
	elem_matrix->matrix[0][6] = 9.0/70.0 - 6*Iyy/(5.0*A*len2);
	elem_matrix->matrix[0][10]= -(13.0*len/420.0 - Iyy/(10.0*A*len));
	elem_matrix->matrix[1][1] = 13.0/35.0 + 6.0*Ixx/(5.0*A*len2);
	elem_matrix->matrix[1][3] = -(11.0*len/210.0 + Ixx/(10.0*A*len));
	elem_matrix->matrix[1][7] = 9.0/70.0 - 6.0*Ixx/(5.0*A*len2);
	elem_matrix->matrix[1][9] = -(-13.0*len/420.0 + Ixx/(10.0*A*len));
	elem_matrix->matrix[2][2] = 1.0/3.0;
	elem_matrix->matrix[2][8] = 1.0/6.0;
	elem_matrix->matrix[3][3] = len2/105.0 + 2.0*Ixx/(15.0*A);
	elem_matrix->matrix[3][7] = -(13.0*len/420.0 - Ixx/(10.0*A*len));
	elem_matrix->matrix[3][9] = -len2/140.0 - Ixx/(30.0*A);
	elem_matrix->matrix[4][4] = len2/105.0 + 2.0*Iyy/(15.0*A);
	elem_matrix->matrix[4][6] = -(-13.0*len/420.0 + Iyy/(10.0*A*len));
	elem_matrix->matrix[4][10]= -len2/140.0 - Iyy/(30.0*A);
	elem_matrix->matrix[5][5] = J/(3.0*A);
	elem_matrix->matrix[5][11]= J/(6.0*A);
	elem_matrix->matrix[6][6] = 13.0/35.0 + 6.0*Iyy/(5.0*A*len2);
	elem_matrix->matrix[6][10]= -(11.0*len/210.0 + Iyy/(10.0*A*len));
	elem_matrix->matrix[7][7] = 13.0/35.0 + 6.0*Ixx/(5.0*A*len2);
	elem_matrix->matrix[7][9] = -(-11.0*len/210.0 - Ixx/(10.0*A*len));
	elem_matrix->matrix[8][8] = 1.0/3.0;
	elem_matrix->matrix[9][9] = len2/105.0 + 2.0*Ixx/(15.0*A);
	elem_matrix->matrix[10][10]=len2/105.0 + 2.0*Iyy/(15.0*A);
	elem_matrix->matrix[11][11]=J/(3.0*A);

	for(ii1=0; ii1<12; ii1++){
		elem_matrix->matrix[ii1][ii1] *= rho*A*len;
		for(ii2=ii1+1; ii2<12; ii2++){
			elem_matrix->matrix[ii1][ii2] *= rho*A*len;
			elem_matrix->matrix[ii2][ii1] = elem_matrix->matrix[ii1][ii2];
		}
	}
	T[0][0] = v1.x;
	T[1][0] = v2.x;
	T[2][0] = pd.x;
	T[0][1] = v1.y;
	T[1][1] = v2.y;
	T[2][1] = pd.y;
	T[0][2] = v1.z;
	T[1][2] = v2.z;
	T[2][2] = pd.z;

	for(ii1=0; ii1<4; ii1++){
		for(ii2=0; ii2<4; ii2++){
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					Tt[ii3][ii4] = elem_matrix->matrix[3*ii1+ii3][3*ii2+ii4];
				}
			}
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					Ttt[ii3][ii4] = T[0][ii3]*Tt[0][ii4]
					              + T[1][ii3]*Tt[1][ii4]
					              + T[2][ii3]*Tt[2][ii4];
				}
			}
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					elem_matrix->matrix[3*ii1+ii3][3*ii2+ii4]
					                   = Ttt[ii3][0]*T[0][ii4]
					                   + Ttt[ii3][1]*T[1][ii4]
					                   + Ttt[ii3][2]*T[2][ii4];
				}
			}
		}
	}

	return(NORMAL_RC);
}
#endif


/* 梁要素の質量マトリックス */
RC fem_make_beam_mass_matrix(FEM_ELEMENT element, FEM_NODE_ARRAY node,
                             FEM_MATERIAL_PROP_ARRAY material,
                             FEM_PHYSICAL_PROP_ARRAY physical,
                             FEM_ELEM_MATRIX *elem_matrix)
{
	int ii1, ii2, ii3, ii4;
	double A, rho, len;
	int p_index, m_index, r_index, n_index;
	TRANSLATION3D p1, p2, pd;
	TRANSLATION3D v1, v2, vtmp;
	double T[3][3];
	double Tt[3][3];
	double Ttt[3][3];

	if(element.label < 0) return(ARG_ERROR_RC);
	if(element.type != ELEM_BEAM1) return(ARG_ERROR_RC);
	elem_matrix->element = element.label;
	elem_matrix->size = 12;
	elem_matrix->dim = 6;

	m_index = search_fem_material_prop_element(material, physical, &element);
	RC_NEG_CHK( m_index );
	if(material.array[m_index].mat_type != MAT_ISOTROPIC) return(ARG_ERROR_RC);
	rho = material.array[m_index].rho;

	p_index = search_fem_physical_prop_label(physical, element.physical);
	RC_NEG_CHK( p_index );
	A   = physical.array[p_index].d_info[0];
	v1.x= element.d_info[0];
	v1.y= element.d_info[1];
	v1.z= element.d_info[2];

	/* elem_matrix->array[][] の確保 */
	RC_TRY( allocate2D(elem_matrix->size,
	                   elem_matrix->size, &(elem_matrix->matrix)) );

	/* elem_matirx->index[] の確保と代入 */
	RC_TRY( rc_malloc(sizeof(int)*(elem_matrix->size),
	                  (void**)(&elem_matrix->index)) );
	for(ii1=0; ii1<element.node_num; ii1++){
		RC_NEG_CHK( r_index = node_renum_index(node, element.node[ii1]) );
		for(ii2=0; ii2<(elem_matrix->dim); ii2++){
			elem_matrix->index[ii1*(elem_matrix->dim)+ii2]
			             = r_index*(elem_matrix->dim) + ii2;
		}
	}

	RC_NEG_CHK( n_index = search_fem_node_label(node, element.node[0]) );
	p1 = node.array[n_index].p;

	RC_NEG_CHK( n_index = search_fem_node_label(node, element.node[1]) );
	p2 = node.array[n_index].p;

	len = abs_translation3d( pd = sub_translation3d(p2, p1) );
	pd = mul_scalar_translation3d(1.0/len, pd);
	vtmp = mul_scalar_translation3d(inner_product3d(v1, pd), pd);
	v1 = sub_translation3d(v1, vtmp);
	v1 = mul_scalar_translation3d(1.0/abs_translation3d(v1), v1);
	v2 = outer_product3d(pd, v1);

	elem_matrix->matrix[0][0]   = 0.5*A*rho*len;
	elem_matrix->matrix[1][1]   = 0.5*A*rho*len;
	elem_matrix->matrix[2][2]   = 0.5*A*rho*len;
	/*elem_matrix->matrix[3][3]   = (A*rho*len)*(0.5*len);
	elem_matrix->matrix[4][4]   = (A*rho*len)*(0.5*len);
	elem_matrix->matrix[5][5]   = (A*rho*len)*(0.5*len);*/
	elem_matrix->matrix[3][3]   = (A*rho*len)*(0.5*len)*(0.5*len);
	elem_matrix->matrix[4][4]   = (A*rho*len)*(0.5*len)*(0.5*len);
	elem_matrix->matrix[5][5]   = 0.0;
	elem_matrix->matrix[6][6]   = 0.5*A*rho*len;
	elem_matrix->matrix[7][7]   = 0.5*A*rho*len;
	elem_matrix->matrix[8][8]   = 0.5*A*rho*len;
	/*elem_matrix->matrix[9][9]   = (A*rho*len)*(0.5*len);
	elem_matrix->matrix[10][10] = (A*rho*len)*(0.5*len);
	elem_matrix->matrix[11][11] = (A*rho*len)*(0.5*len);*/
	elem_matrix->matrix[9][9]   = (A*rho*len)*(0.5*len)*(0.5*len);
	elem_matrix->matrix[10][10] = (A*rho*len)*(0.5*len)*(0.5*len);
	elem_matrix->matrix[11][11] = 0.0;

	for(ii1=0; ii1<12; ii1++){
		for(ii2=ii1+1; ii2<12; ii2++){
			elem_matrix->matrix[ii2][ii1] = elem_matrix->matrix[ii1][ii2];
		}
	}
	T[0][0] = v1.x;
	T[1][0] = v2.x;
	T[2][0] = pd.x;
	T[0][1] = v1.y;
	T[1][1] = v2.y;
	T[2][1] = pd.y;
	T[0][2] = v1.z;
	T[1][2] = v2.z;
	T[2][2] = pd.z;

	for(ii1=0; ii1<4; ii1++){
		for(ii2=0; ii2<4; ii2++){
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					Tt[ii3][ii4] = elem_matrix->matrix[3*ii1+ii3][3*ii2+ii4];
				}
			}
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					Ttt[ii3][ii4] = T[0][ii3]*Tt[0][ii4]
					              + T[1][ii3]*Tt[1][ii4]
					              + T[2][ii3]*Tt[2][ii4];
				}
			}
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					elem_matrix->matrix[3*ii1+ii3][3*ii2+ii4]
					                   = Ttt[ii3][0]*T[0][ii4]
					                   + Ttt[ii3][1]*T[1][ii4]
					                   + Ttt[ii3][2]*T[2][ii4];
				}
			}
		}
	}

	return(NORMAL_RC);
}


/* 要素剛性マトリックスの作成 */
RC fem_make_beam_matrix(FEM_ELEMENT element, FEM_NODE_ARRAY node,
                        FEM_MATERIAL_PROP_ARRAY material,
                        FEM_PHYSICAL_PROP_ARRAY physical,
                        FEM_ELEM_MATRIX *elem_matrix)
{
	int ii1, ii2, ii3, ii4;
	double E, G;
	double A, Ixx, Iyy, J;
	double len, len2, len3;
	int p_index, m_index, r_index, n_index;
	TRANSLATION3D p1, p2, pd;
	TRANSLATION3D v1, v2, vtmp;
	double T[3][3];
	double Tt[3][3];
	double Ttt[3][3];

	if(element.label < 0) return(ARG_ERROR_RC);
	if(element.type != ELEM_BEAM1) return(ARG_ERROR_RC);
	elem_matrix->element = element.label;
	elem_matrix->size = 12;
	elem_matrix->dim = 6;

	m_index = search_fem_material_prop_element(material, physical, &element);
	RC_NEG_CHK( m_index );
	if(material.array[m_index].mat_type != MAT_ISOTROPIC) return(ARG_ERROR_RC);
	E = material.array[m_index].E;
	G = E/(2.0*(1.0+material.array[m_index].nu));

	p_index = search_fem_physical_prop_label(physical, element.physical);
	RC_NEG_CHK( p_index );
	A   = physical.array[p_index].d_info[0];
	Ixx = physical.array[p_index].d_info[1];
	Iyy = physical.array[p_index].d_info[2];
	J   = physical.array[p_index].d_info[3];
	v1.x= element.d_info[0];
	v1.y= element.d_info[1];
	v1.z= element.d_info[2];

	/* elem_matrix->array[][] の確保 */
	RC_TRY( allocate2D(elem_matrix->size,
	                   elem_matrix->size, &(elem_matrix->matrix)) );

	/* elem_matirx->index[] の確保と代入 */
	RC_TRY( rc_malloc(sizeof(int)*(elem_matrix->size),
	                  (void**)(&elem_matrix->index)) );
	for(ii1=0; ii1<element.node_num; ii1++){
		RC_NEG_CHK( r_index = node_renum_index(node, element.node[ii1]) );
		for(ii2=0; ii2<(elem_matrix->dim); ii2++){
			elem_matrix->index[ii1*(elem_matrix->dim)+ii2]
			             = r_index*(elem_matrix->dim) + ii2;
		}
	}

	RC_NEG_CHK( n_index = search_fem_node_label(node, element.node[0]) );
	p1 = node.array[n_index].p;

	RC_NEG_CHK( n_index = search_fem_node_label(node, element.node[1]) );
	p2 = node.array[n_index].p;

	len = abs_translation3d( pd = sub_translation3d(p2, p1) );
	len2 = len * len;
	len3 = len * len * len;
	pd = mul_scalar_translation3d(1.0/len, pd);
	vtmp = mul_scalar_translation3d(inner_product3d(v1, pd), pd);
	v1 = sub_translation3d(v1, vtmp);
	v1 = mul_scalar_translation3d(1.0/abs_translation3d(v1), v1);
	v2 = outer_product3d(pd, v1);

	elem_matrix->matrix[0][0] = 12.0*E*Ixx/len3;
	elem_matrix->matrix[0][4] = 6.0*E*Ixx/len2;
	elem_matrix->matrix[0][6] = -12.0*E*Ixx/len3;
	elem_matrix->matrix[0][10] = 6.0*E*Ixx/len2;
	elem_matrix->matrix[1][1] = 12.0*E*Iyy/len3;
	elem_matrix->matrix[1][3] = -6.0*E*Iyy/len2;
	elem_matrix->matrix[1][7] = -12.0*E*Iyy/len3;
	elem_matrix->matrix[1][9] = -6.0*E*Iyy/len2;
	elem_matrix->matrix[2][2] = E*A/len;
	elem_matrix->matrix[2][8] = -E*A/len;
	elem_matrix->matrix[3][3] = 4.0*E*Iyy/len;
	elem_matrix->matrix[3][7] = 6.0*E*Iyy/len2;
	elem_matrix->matrix[3][9] = 2.0*E*Iyy/len;
	elem_matrix->matrix[4][4] = 4.0*E*Ixx/len;
	elem_matrix->matrix[4][6] = -6.0*E*Ixx/len2;
	elem_matrix->matrix[4][10] = 2.0*E*Ixx/len;
	elem_matrix->matrix[5][5] = G*J/len;
	elem_matrix->matrix[5][11] = -G*J/len;
	elem_matrix->matrix[6][6] = 12.0*E*Ixx/len3;
	elem_matrix->matrix[6][10] = -6.0*E*Ixx/len2;
	elem_matrix->matrix[7][7] = 12.0*E*Iyy/len3;
	elem_matrix->matrix[7][9] = 6.0*E*Iyy/len2;
	elem_matrix->matrix[8][8] = E*A/len;
	elem_matrix->matrix[9][9] = 4.0*E*Iyy/len;
	elem_matrix->matrix[10][10] = 4.0*E*Ixx/len;
	elem_matrix->matrix[11][11] = G*J/len;

	for(ii1=0; ii1<12; ii1++){
		for(ii2=ii1+1; ii2<12; ii2++){
			elem_matrix->matrix[ii2][ii1] = elem_matrix->matrix[ii1][ii2];
		}
	}
	T[0][0] = v1.x;
	T[1][0] = v2.x;
	T[2][0] = pd.x;
	T[0][1] = v1.y;
	T[1][1] = v2.y;
	T[2][1] = pd.y;
	T[0][2] = v1.z;
	T[1][2] = v2.z;
	T[2][2] = pd.z;

	for(ii1=0; ii1<4; ii1++){
		for(ii2=0; ii2<4; ii2++){
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					Tt[ii3][ii4] = elem_matrix->matrix[3*ii1+ii3][3*ii2+ii4];
				}
			}
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					Ttt[ii3][ii4] = T[0][ii3]*Tt[0][ii4]
					              + T[1][ii3]*Tt[1][ii4]
					              + T[2][ii3]*Tt[2][ii4];
				}
			}
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					elem_matrix->matrix[3*ii1+ii3][3*ii2+ii4]
					                   = Ttt[ii3][0]*T[0][ii4]
					                   + Ttt[ii3][1]*T[1][ii4]
					                   + Ttt[ii3][2]*T[2][ii4];
				}
			}
		}
	}

	return(NORMAL_RC);
}


/* 節点バネ要素の処理 */
RC fem_modify_spring6(FEM_NODE_ARRAY node, FEM_RIGID_ELEMENT_ARRAY src_rigid,
                      FEM_PHYSICAL_PROP_ARRAY physical,
                      NONZERO_MATRIX3 matrix)
{
	int ii1, ii2, ii3;
	int index1, index2;
	int dof1, dof2;
	NONZERO_ELEM3 *ptr;
	double spring_const = 0.0;

	for(ii1=0; ii1<src_rigid.size; ii1++){
		if(src_rigid.array[ii1].label < 0) continue;
		if(src_rigid.array[ii1].type != ELEM_SPRING) continue;
		if(src_rigid.array[ii1].physical < 0){
			spring_const = src_rigid.array[ii1].weight[0];
		}else{
			RC_TRY( get_spring_const(physical, src_rigid.array[ii1].physical,
			                         &spring_const) );
		}
		for(ii2=0; ii2<src_rigid.array[ii1].node_num; ii2++){
			RC_TRY( cal_index_dof6(node, src_rigid.array[ii1].node[ii2],
			                 src_rigid.array[ii1].dof[ii2]-1, &index1, &dof1) );
			ptr = nonzero_matrix3s_ptr(matrix, index1, index1);
			RC_NULL_CHK(ptr);
			(ptr->v)[dof1][dof1] += spring_const;

			for(ii3=0; ii3<src_rigid.array[ii1].node_num; ii3++){
				if(ii2 == ii3) continue;
				RC_TRY( cal_index_dof6(node, src_rigid.array[ii1].node[ii3],
				            src_rigid.array[ii1].dof[ii3]-1, &index2, &dof2) );
				if(index1 < index2) continue;
				ptr = nonzero_matrix3s_ptr(matrix, index1, index2);
				RC_NULL_CHK(ptr);
				(ptr->v)[dof1][dof2] -= spring_const;
			}
		}
	}

	return(NORMAL_RC);
}


RC fem_modify_node_matrix6(FEM_NODE_ARRAY node, FEM_NODE_MATRIX_ARRAY nmat,
                           double fact, NONZERO_MATRIX3 matrix)
{
	int ii1, ii2;
	int row, col;
	NONZERO_ELEM3 *ptr;

	for(ii1=0; ii1<nmat.size; ii1++){
		if(nmat.array[ii1].node < 0) continue;
		RC_NEG_CHK( row = 2*node_renum_index(node, nmat.array[ii1].node) );
		row += nmat.array[ii1].shift;
		for(ii2=0; ii2<nmat.array[ii1].size; ii2++){
			col = 2*node_renum_index(node, nmat.array[ii1].col_node[ii2]);
			RC_NEG_CHK( col );
			col += nmat.array[ii1].col_shift[ii2];
			if(col > row){
				RC_NULL_CHK( ptr = nonzero_matrix3s_ptr(matrix, col, row) );
				ptr->v[0][0] += fact * nmat.array[ii1].matrix[ii2][0][0];
				ptr->v[0][1] += fact * nmat.array[ii1].matrix[ii2][1][0];
				ptr->v[0][2] += fact * nmat.array[ii1].matrix[ii2][2][0];
				ptr->v[1][0] += fact * nmat.array[ii1].matrix[ii2][0][1];
				ptr->v[1][1] += fact * nmat.array[ii1].matrix[ii2][1][1];
				ptr->v[1][2] += fact * nmat.array[ii1].matrix[ii2][2][1];
				ptr->v[2][0] += fact * nmat.array[ii1].matrix[ii2][0][2];
				ptr->v[2][1] += fact * nmat.array[ii1].matrix[ii2][1][2];
				ptr->v[2][2] += fact * nmat.array[ii1].matrix[ii2][2][2];
			}else{
				RC_NULL_CHK( ptr = nonzero_matrix3s_ptr(matrix, row, col) );
				ptr->v[0][0] += fact * nmat.array[ii1].matrix[ii2][0][0];
				ptr->v[0][1] += fact * nmat.array[ii1].matrix[ii2][0][1];
				ptr->v[0][2] += fact * nmat.array[ii1].matrix[ii2][0][2];
				ptr->v[1][0] += fact * nmat.array[ii1].matrix[ii2][1][0];
				ptr->v[1][1] += fact * nmat.array[ii1].matrix[ii2][1][1];
				ptr->v[1][2] += fact * nmat.array[ii1].matrix[ii2][1][2];
				ptr->v[2][0] += fact * nmat.array[ii1].matrix[ii2][2][0];
				ptr->v[2][1] += fact * nmat.array[ii1].matrix[ii2][2][1];
				ptr->v[2][2] += fact * nmat.array[ii1].matrix[ii2][2][2];
			}
		}
	}

	return(NORMAL_RC);
}


/* 集中質量要素の処理 */
/* fact: 係数 */
RC fem_modify_mass6(FEM_NODE_ARRAY node, FEM_RIGID_ELEMENT_ARRAY rigid,
                    double fact, NONZERO_MATRIX3 matrix)
{
	int ii1;
	int index;
	NONZERO_ELEM3 *ptr;

	for(ii1=0; ii1<rigid.size; ii1++){
		if(rigid.array[ii1].label < 0) continue;
		if(rigid.array[ii1].type != ELEM_MASS) continue;
		if(rigid.array[ii1].node_num != 1) return(ARG_ERROR_RC);
		index = 2*node_renum_index(node, rigid.array[ii1].node[0]);
		RC_NEG_CHK( index );

		RC_NULL_CHK( ptr = nonzero_matrix3s_ptr(matrix, index, index) );
		ptr->v[0][0] += fact * rigid.array[ii1].weight[0];
		ptr->v[0][1] += ptr->v[1][0] = fact * rigid.array[ii1].weight[1];
		ptr->v[1][1] += fact * rigid.array[ii1].weight[2];
		ptr->v[0][2] += ptr->v[2][0] = fact * rigid.array[ii1].weight[3];
		ptr->v[1][2] += ptr->v[2][1] = fact * rigid.array[ii1].weight[4];
		ptr->v[2][2] += fact * rigid.array[ii1].weight[5];

		RC_NULL_CHK( ptr = nonzero_matrix3s_ptr(matrix, index+1, index) );
		ptr->v[0][0] += fact * rigid.array[ii1].weight[6];
		ptr->v[0][1] += fact * rigid.array[ii1].weight[7];
		ptr->v[0][2] += fact * rigid.array[ii1].weight[8];
		ptr->v[1][0] += fact * rigid.array[ii1].weight[10];
		ptr->v[1][1] += fact * rigid.array[ii1].weight[11];
		ptr->v[1][2] += fact * rigid.array[ii1].weight[12];
		ptr->v[2][0] += fact * rigid.array[ii1].weight[15];
		ptr->v[2][1] += fact * rigid.array[ii1].weight[16];
		ptr->v[2][2] += fact * rigid.array[ii1].weight[17];

		RC_NULL_CHK( ptr = nonzero_matrix3s_ptr(matrix, index+1, index+1) );
		ptr->v[0][0] += fact * rigid.array[ii1].weight[9];
		ptr->v[0][1] += ptr->v[1][0] = fact * rigid.array[ii1].weight[13];
		ptr->v[1][1] += fact * rigid.array[ii1].weight[14];
		ptr->v[0][2] += ptr->v[2][0] = fact * rigid.array[ii1].weight[18];
		ptr->v[1][2] += ptr->v[2][1] = fact * rigid.array[ii1].weight[19];
		ptr->v[2][2] += fact * rigid.array[ii1].weight[20];
	}

	return(NORMAL_RC);
}


RC fem_modify_nonzero6(FEM_NODE_ARRAY node, FEM_LOCAL_COORD_ARRAY coord,
                       FEM_BC_ARRAY rest, FEM_RIGID_ELEMENT_ARRAY src_rigid,
                       NONZERO_MATRIX3 matrix,
                       double **f_vect, int vect_size, int active_dof_flag[])
{
	int ii1, ii2;
	int label0, label;
	int index0, index;
	double weight;
	int dof0, dof;
	int *slave_table = NULL;
	int *slave_flag  = NULL;
	int *master_flag = NULL;
	double *master_weight = NULL;
	int master_count;
	FEM_RIGID_ELEMENT_ARRAY rigid;
	double s_matrix[3][3];

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	if(active_dof_flag != NULL){
		for(ii1=0; ii1<3*matrix.size; ii1++){
			active_dof_flag[ii1] = 1;
		}
	}

	/* ローカル座標系の処理 */
	for(ii1=0; ii1<node.size; ii1++){
		if(node.array[ii1].label < 0) continue;

		if( (node.source == FEM_NASTRAN)
		  &&(node.array[ii1].i_info[1] > 0) ){
			int coord_index
			     = search_fem_local_coord_label(coord,
			                                    node.array[ii1].i_info[1]);
			index = 2 * node_renum_index(node, node.array[ii1].label);

			RC_NEG_CHK(coord_index);
			RC_NEG_CHK(index);
			local_coord2s_matrix_t(coord.array[coord_index], s_matrix);
			RC_TRY( row_col_mul_nonzero3s(matrix, index, s_matrix) );
			RC_TRY( row_col_mul_nonzero3s(matrix, index+1, s_matrix) );
			for(ii2=0; ii2<vect_size; ii2++){
				mul_matrix33_vect(s_matrix, &(f_vect[ii2][3*index]));
				mul_matrix33_vect(s_matrix, &(f_vect[ii2][3*(index+1)]));
			}
		}
	}

	if(src_rigid.size > 0){
		RC_NULL_CHK( slave_table = (int *)malloc(6*node.size*sizeof(int)) );
		RC_NULL_CHK( master_flag = (int *)malloc(6*node.size*sizeof(int)) );
		RC_NULL_CHK( master_weight
		             = (double *)malloc(6*node.size*sizeof(double)) );
		RC_TRY( expand_rigid(node, src_rigid, &rigid, slave_table) );
		RC_NULL_CHK( slave_flag  = (int *)malloc(rigid.size*sizeof(int)) );
	}else{
		rigid.size = 0;
	}

	/* 剛体要素，多点拘束 */
	for(ii1=0; ii1<rigid.size; ii1++){
		if(rigid.array[ii1].label < 0) continue;
		if(rigid.array[ii1].type != ELEM_MPC) return(IMPLEMENT_ERROR_RC);
		if(rigid.array[ii1].node_num <= 1) return(ARG_ERROR_RC);

		label0 = rigid.array[ii1].node[0];
		dof0 = rigid.array[ii1].dof[0];
		master_count = 0;
		for(ii2=1; ii2<rigid.array[ii1].node_num; ii2++){
			label = rigid.array[ii1].node[ii2];
			dof = rigid.array[ii1].dof[ii2];
			weight = rigid.array[ii1].weight[ii2];
			init_flags(slave_flag, rigid.size, ii1);
			RC_TRY( rigid_modify_nonzero6(matrix, node, rigid,
			             slave_flag, master_flag, master_weight,
			             &master_count, f_vect, vect_size, weight,
			             label0, dof0, label, dof, slave_table) );
		}
		RC_TRY( cal_index_dof6(node, rigid.array[ii1].node[0],
		                       rigid.array[ii1].dof[0], &index0, &dof0) );
		RC_TRY( add_master_dof(matrix, master_count,
		                       master_flag, master_weight, index0, dof0) );

		/* slave 自由度を消去 */
		RC_TRY( row_col_zero_nonzero3s(matrix, index0, dof0) );
		for(ii2=0; ii2<vect_size; ii2++){
			f_vect[ii2][3*index0+dof0] = 0.0;
		}
		if(active_dof_flag != NULL) active_dof_flag[3*index0+dof0] = 0;
	}

	/* 単点拘束 */
	RC_TRY( delete_zero_matrix3s(matrix) );
	for(ii1=0; ii1<rest.size; ii1++){
		int index;
		if(rest.array[ii1].node < 0) continue;

		RC_NEG_CHK( index = 2 * node_renum_index(node, rest.array[ii1].node) );

		if(rest.array[ii1].v_type.x == BC_FIX){
			RC_TRY( modify_nonzero3s_n(matrix, index, 0, rest.array[ii1].v.x,
			                           f_vect, vect_size) );
			if(active_dof_flag != NULL) active_dof_flag[3*index + 0] = 0;
		}
		if(rest.array[ii1].v_type.y == BC_FIX){
			RC_TRY( modify_nonzero3s_n(matrix, index, 1, rest.array[ii1].v.y,
			                           f_vect, vect_size) );
			if(active_dof_flag != NULL) active_dof_flag[3*index + 1] = 0;
		}
		if(rest.array[ii1].v_type.z == BC_FIX){
			RC_TRY( modify_nonzero3s_n(matrix, index, 2, rest.array[ii1].v.z,
			                           f_vect, vect_size) );
			if(active_dof_flag != NULL) active_dof_flag[3*index + 2] = 0;
		}

		if(rest.array[ii1].v_type.yz == BC_FIX){
			RC_TRY( modify_nonzero3s_n(matrix, index+1, 0, rest.array[ii1].v.yz,
			                           f_vect, vect_size) );
			if(active_dof_flag != NULL) active_dof_flag[3*(index+1) + 0] = 0;
		}
		if(rest.array[ii1].v_type.zx == BC_FIX){
			RC_TRY( modify_nonzero3s_n(matrix, index+1, 1, rest.array[ii1].v.zx,
			                           f_vect, vect_size) );
			if(active_dof_flag != NULL) active_dof_flag[3*(index+1) + 1] = 0;
		}
		if(rest.array[ii1].v_type.xy == BC_FIX){
			RC_TRY( modify_nonzero3s_n(matrix, index+1, 2, rest.array[ii1].v.xy,
			                           f_vect, vect_size) );
			if(active_dof_flag != NULL) active_dof_flag[3*(index+1) + 2] = 0;
		}
	}

	if(rigid.size > 0){
		free(slave_table);
		free(slave_flag);
		free(master_flag);
		free(master_weight);
		RC_TRY( free_fem_rigid_element_array(&rigid) );
	}

	/* 対角要素に小さな正数を加算し、正定値性を確保 */
	RC_TRY( delete_zero_matrix3s(matrix) );
	for(ii1=0; ii1<matrix.size; ii1++){
		int dia_index = matrix.row[ii1].size - 1;

		matrix.row[ii1].elem[dia_index].v[0][0] += ABS_TOL;
		matrix.row[ii1].elem[dia_index].v[1][1] += ABS_TOL;
		matrix.row[ii1].elem[dia_index].v[2][2] += ABS_TOL;
		if(active_dof_flag != NULL){
			if(matrix.row[ii1].elem[dia_index].v[0][0] < 10.0*ABS_TOL){
				active_dof_flag[3*ii1 + 0] = 0;
			}
			if(matrix.row[ii1].elem[dia_index].v[1][1] < 10.0*ABS_TOL){
				active_dof_flag[3*ii1 + 1] = 0;
			}
			if(matrix.row[ii1].elem[dia_index].v[2][2] < 10.0*ABS_TOL){
				active_dof_flag[3*ii1 + 2] = 0;
			}
		}
	}

	return(NORMAL_RC);
}

static RC cal_index_dof6(FEM_NODE_ARRAY node, int label, int src_dof,
                         int *index, int *dof)
{
	if((src_dof < 0)||(src_dof >= 6)) return(ARG_ERROR_RC);
	*index = 2 * node_renum_index(node, label);
	RC_NEG_CHK(*index);
	*dof = src_dof;
	if(*dof >= 3){
		*index += 1;
		*dof -= 3;
	}

	return(NORMAL_RC);
}


/* master_flag[], master_weight[] を基にマトリックスを修正 */
static RC add_master_dof(NONZERO_MATRIX3 matrix, int master_count,
                        int master_flag[], double master_weight[],
                        int index0, int dof0)
{
	int count = 0;
	int ii1, ii2, ii3;
	double dia0 = matrix.row[index0]
	                    .elem[matrix.row[index0].size - 1].v[dof0][dof0];

	for(ii1=0; ii1<master_count; ii1++){
		for(ii2=0; ii2<master_count; ii2++){
			int row = master_flag[ii1];
			int col = master_flag[ii2];

			if(ii1 == ii2) continue;
			if(row/3 >= col/3){
				int index = -1;
				for(ii3=0; ii3<matrix.row[row/3].size; ii3++){
					if(matrix.row[row/3].elem[ii3].col == col/3){
						index = ii3;
						break;
					}
				}
				RC_NEG_CHK( index );
				matrix.row[row/3].elem[index].v[row%3][col%3]
				   += master_weight[ii1] * master_weight[ii2] * dia0;
				count++;
			}
		}
	}

	return(NORMAL_RC);
}


static int dof_chk(int dof, int i)
{
	while(dof != 0){
		if(dof%10 == i) return(1);
		dof /= 10;
	}
	
fprintf(stdout, "???\n");
	return(0);
}

static void init_flags(int flags[], int size, int i)
{
	int ii1;

	for(ii1=0; ii1<size; ii1++){
		flags[ii1] = 0;
	}
	if(i < size) flags[i] = 1;
}


/* スレーブ自由度にマスター自由度を反映させる */
RC fem_rigid_recover6(FEM_NODE_ARRAY node,
                      FEM_RIGID_ELEMENT_ARRAY src_rigid, double d_vect[])
{
	FEM_RIGID_ELEMENT_ARRAY rigid;
	int *slave_table;

	RC_NULL_CHK( slave_table = (int *)malloc(6*node.size*sizeof(int)) );
	RC_TRY( expand_rigid(node, src_rigid, &rigid, slave_table) );
	RC_TRY( fem_rigid_recover(6, node, rigid, d_vect, slave_table) );
	RC_TRY( free_fem_rigid_element_array(&rigid) );
	free(slave_table);

	return(NORMAL_RC);
}


/* スレーブ自由度にマスター自由度を反映させる */
static RC fem_rigid_recover(int gdim, FEM_NODE_ARRAY node,
                            FEM_RIGID_ELEMENT_ARRAY rigid, double d_vect[],
                            const int slave_table[])
{
	int ii1, ii2, ii3;
	int index0, index1;
	int dof0, dof1;
	int *slave_flag = NULL;

	if(rigid.size <= 0) return(NORMAL_RC);

	RC_NULL_CHK( slave_flag  = (int *)malloc(rigid.size*sizeof(int)) );

	for(ii1=0; ii1<rigid.size; ii1++){
		if(rigid.array[ii1].label < 0) continue;
		if(rigid.array[ii1].type != ELEM_MPC) continue;

		RC_NEG_CHK( index0 = node_renum_index(node,
					              rigid.array[ii1].node[0]) );
		dof0 = rigid.array[ii1].dof[0];
		if((dof0 < 0)||(dof0 > gdim)) return(ARG_ERROR_RC);

		for(ii2=1; ii2<rigid.array[ii1].node_num; ii2++){
			RC_NEG_CHK( index1 = node_renum_index(node,
						              rigid.array[ii1].node[ii2]) );
			dof1 = rigid.array[ii1].dof[ii2];
			if((dof1 < 0)||(dof1 > gdim)) return(ARG_ERROR_RC);

			for(ii3=0; ii3<rigid.size; ii3++){
				slave_flag[ii3] = 0;
			}
			RC_TRY( rigid_recover_sub(gdim, node, slave_flag, rigid,
                                  rigid.array[ii1].weight[ii2], index0, dof0,
			                      index1, dof1, d_vect, slave_table) );
		}
	}

	free(slave_flag);

	return(NORMAL_RC);
}


static RC rigid_recover_sub(int gdim, FEM_NODE_ARRAY node, int slave_flag[],
                            FEM_RIGID_ELEMENT_ARRAY rigid, double weight,
                            int index0, int dof0, int index1, int dof1,
                            double d_vect[], const int slave_table[])
{
	int ii1;

	if(slave_table[gdim*index1 + dof1] >= 0){
		int r_index = slave_table[gdim*index1 + dof1];

		if(slave_flag[r_index]) return(UNKNOWN_ERROR_RC);
		slave_flag[r_index] = 1;
		for(ii1=1; ii1<rigid.array[r_index].node_num; ii1++){
			int index1_new, dof1_new;
			double weight_new;

			RC_NEG_CHK( index1_new = node_renum_index(node,
			                         rigid.array[r_index].node[ii1]) );
			dof1_new = rigid.array[r_index].dof[ii1];
			weight_new = rigid.array[r_index].weight[ii1];
			if((dof1_new < 0)||(dof1_new > gdim)) return(ARG_ERROR_RC);
			RC_TRY( rigid_recover_sub(gdim, node, slave_flag, rigid,
			                     weight*weight_new, index0, dof0,
			                     index1_new, dof1_new, d_vect, slave_table) );
		}
	}else{
		d_vect[gdim*index0 + dof0] += weight * d_vect[gdim*index1 + dof1];
	}

	return(NORMAL_RC);
}


/* 荷重ベクトルに節点荷重を加算 */
RC fem_add_nodal_force6(FEM_NODE_ARRAY node, FEM_LOCAL_COORD_ARRAY coord,
                        FEM_BC_ARRAY force, double *f_vect)
{
	int index;
	int ii1;

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	for(ii1=0; ii1<force.size; ii1++){
		double tmp_v1[3] = {0.0, 0.0, 0.0};
		double tmp_v2[3] = {0.0, 0.0, 0.0};
		if(force.array[ii1].node < 0) continue;

		RC_NEG_CHK( index = 6*node_renum_index(node, force.array[ii1].node) );
		if(force.array[ii1].v_type.x == BC_FORCE){
			tmp_v1[0] = force.array[ii1].v.x;
		}
		if(force.array[ii1].v_type.y == BC_FORCE){
			tmp_v1[1] = force.array[ii1].v.y;
		}
		if(force.array[ii1].v_type.z == BC_FORCE){
			tmp_v1[2] = force.array[ii1].v.z;
		}
		if(force.array[ii1].v_type.yz == BC_FORCE){
			tmp_v2[0] = force.array[ii1].v.yz;
		}
		if(force.array[ii1].v_type.zx == BC_FORCE){
			tmp_v2[1] = force.array[ii1].v.zx;
		}
		if(force.array[ii1].v_type.xy == BC_FORCE){
			tmp_v2[2] = force.array[ii1].v.xy;
		}
		/* ローカル座標系の処理 */
		if( (force.source == FEM_NASTRAN)
		  &&(force.array[ii1].i_info[0] > 0) ){
			double s_matrix[3][3];
			int coord_index = search_fem_local_coord_label(coord,
			                          force.array[ii1].i_info[0]);
			RC_NEG_CHK(coord_index);
			local_coord2s_matrix(coord.array[coord_index], s_matrix);
			mul_matrix33_vect(s_matrix, tmp_v1);
			mul_matrix33_vect(s_matrix, tmp_v2);
		}
		f_vect[index]     += tmp_v1[0];
		f_vect[index + 1] += tmp_v1[1];
		f_vect[index + 2] += tmp_v1[2];
		f_vect[index + 3] += tmp_v2[0];
		f_vect[index + 4] += tmp_v2[1];
		f_vect[index + 5] += tmp_v2[2];
	}

	return(NORMAL_RC);
}


/* 荷重ベクトルに圧力/表面荷重を加算 */
RC fem_add_pressure_force6(FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY element,
                           FEM_LOCAL_COORD_ARRAY coord,
                           FEM_ELEM_BC_ARRAY press, double *f_vect)
{
	int ii1, ii2;
	int p_size = 1;
	int elem_index;
	FEM_FACE_TABLE f_table;
	FEM_ELEMENT face;
	double p[FEM_MAX_NODE];
	TRANSLATION3D f[FEM_MAX_NODE];

	if(node.renum_flag != RENUMBERED) return(ARG_ERROR_RC);

	for(ii1=0; ii1<press.size; ii1++){
		if(press.array[ii1].element < 0) continue;
		if( (press.array[ii1].type != BC_PRESS)
		  &&(press.array[ii1].type != BC_PRESS_D1)
		  &&(press.array[ii1].type != BC_PRESS_D2)
		  &&(press.array[ii1].type != BC_TRACTION)
		  &&(press.array[ii1].type != BC_TRACTION_D1)
		  &&(press.array[ii1].type != BC_TRACTION_D2) ) continue;

		elem_index = search_fem_element_label(element,
		                        press.array[ii1].element);
		RC_NEG_CHK( elem_index );
		RC_TRY( make_fem_face_table(element.array[elem_index].type, &f_table) );
		if(press.array[ii1].face_id >= f_table.face_num) return(ARG_ERROR_RC);
		face.label = 1;
		face.node_num = f_table.node_num[press.array[ii1].face_id];
		face.type = f_table.face_type[press.array[ii1].face_id];
		for(ii2=0; ii2<f_table.node_num[press.array[ii1].face_id]; ii2++){
			face.node[ii2] = element.array[elem_index]
			                .node[f_table.table[press.array[ii1].face_id][ii2]];
		}
		if( (press.array[ii1].type == BC_PRESS)
		  ||(press.array[ii1].type == BC_TRACTION) ){
			p_size = 1;
			p[0] = -(press.array[ii1].val.s);
		}else if( (press.array[ii1].type == BC_PRESS_D1)
		        ||(press.array[ii1].type == BC_TRACTION_D1) ){
			if( (face.type == ELEM_QUAD1)||(face.type == ELEM_QUAD2) ){
				p_size = 4;
			}else if( (face.type == ELEM_TRI1)||(face.type == ELEM_TRI2) ){
				p_size = 3;
			}else{
				return(IMPLEMENT_ERROR_RC);
			}
			for(ii2=0; ii2<p_size; ii2++){
				p[ii2] = -(press.array[ii1].val.ns[ii2]);
			}
		}else if( (press.array[ii1].type == BC_PRESS_D2)
		        ||(press.array[ii1].type == BC_TRACTION_D2) ){
			if( (face.type == ELEM_QUAD1)||(face.type == ELEM_QUAD2) ){
				p_size = 8;
			}else if( (face.type == ELEM_TRI1)||(face.type == ELEM_TRI2) ){
				p_size = 6;
			}else{
				return(IMPLEMENT_ERROR_RC);
			}
			for(ii2=0; ii2<p_size; ii2++){
				p[ii2] = -(press.array[ii1].val.ns[ii2]);
			}
		}

		if( (press.array[ii1].type == BC_PRESS)
		  ||(press.array[ii1].type == BC_PRESS_D1)
		  ||(press.array[ii1].type == BC_PRESS_D2) ){
			RC_TRY( pressure_force(p_size, p, &face, node, f) );
		}else{
			TRANSLATION3D gdir = press.array[ii1].dir;

			if( (press.source == FEM_NASTRAN)
			  &&(press.array[ii1].i_info[0] > 0) ){
				int index = search_fem_local_coord_label(coord,
				                       press.array[ii1].i_info[0]);
				RC_NEG_CHK(index);
				gdir = trans_local2global(coord.array[index],
				                          press.array[ii1].dir);
			}
			RC_TRY( traction_force(p_size, p, gdir, &face, node, f) );
		}

		for(ii2=0; ii2<face.node_num; ii2++){
			int index = node_renum_index(node, face.node[ii2]);
			RC_NEG_CHK( index );
			f_vect[6*index] += f[ii2].x;
			f_vect[6*index + 1] += f[ii2].y;
			f_vect[6*index + 2] += f[ii2].z;
		}
	}

	return(NORMAL_RC);
}


static void mul_matrix33_vect(double s_matrix[][3], double vect[])
{
	double tmp_vect[3];

	tmp_vect[0] = s_matrix[0][0] * vect[0]
	            + s_matrix[0][1] * vect[1]
	            + s_matrix[0][2] * vect[2];
	tmp_vect[1] = s_matrix[1][0] * vect[0]
	            + s_matrix[1][1] * vect[1]
	            + s_matrix[1][2] * vect[2];
	tmp_vect[2] = s_matrix[2][0] * vect[0]
	            + s_matrix[2][1] * vect[1]
	            + s_matrix[2][2] * vect[2];
	vect[0] = tmp_vect[0];
	vect[1] = tmp_vect[1];
	vect[2] = tmp_vect[2];
}


/* 荷重ベクトルを節点の座標系に合わせる */
RC fem_force_local_pre6(FEM_NODE_ARRAY node, FEM_LOCAL_COORD_ARRAY coord,
                        double f_vect[])
{
	int ii1;
	int coord_index;
	int index;
	double s_matrix[3][3];

	for(ii1=0; ii1<node.size; ii1++){
		if(node.array[ii1].label < 0) continue;

		if( (node.source == FEM_NASTRAN)&&(node.array[ii1].i_info[1] > 0) ){
			coord_index = search_fem_local_coord_label(coord,
			                                    node.array[ii1].i_info[1]);
			index = 6 * node_renum_index(node, node.array[ii1].label);

			RC_NEG_CHK(coord_index);
			RC_NEG_CHK(index);
			local_coord2s_matrix_t(coord.array[coord_index], s_matrix);
			mul_matrix33_vect(s_matrix, &(f_vect[index]));
			mul_matrix33_vect(s_matrix, &(f_vect[index+3]));
		}
	}

	return(NORMAL_RC);
}


/* 変位ベクトルを全体座標系に戻す */
RC fem_disp_local_recover6(FEM_NODE_ARRAY node, FEM_LOCAL_COORD_ARRAY coord,
                           double d_vect[])
{
	int ii1;
	int coord_index;
	int index;
	double s_matrix[3][3];

	for(ii1=0; ii1<node.size; ii1++){
		if(node.array[ii1].label < 0) continue;

		if( (node.source == FEM_NASTRAN)&&(node.array[ii1].i_info[1] > 0) ){
			coord_index = search_fem_local_coord_label(coord,
			                                    node.array[ii1].i_info[1]);
			index = 6 * node_renum_index(node, node.array[ii1].label);

			RC_NEG_CHK(coord_index);
			RC_NEG_CHK(index);
			local_coord2s_matrix(coord.array[coord_index], s_matrix);
			mul_matrix33_vect(s_matrix, &(d_vect[index]));
			mul_matrix33_vect(s_matrix, &(d_vect[index+3]));
		}
	}

	return(NORMAL_RC);
}


/* 非ゼロ形式からスカイライン形式に変換 */
/* allocatable_size : 確保する要素数の上限または -1 */
RC fem_nonzero2skyline6(NONZERO_MATRIX3 nonzero, FEM_SKYLINE_MATRIX *skyline,
                        long allocatable_size)
{
	int ii1, ii2, ii3, ii4;
	long alloc_size = 0L;

	skyline->dim = 6;
	skyline->dof = 3 * nonzero.size;
	RC_NULL_CHK( skyline->index1 = (long *)malloc(skyline->dof*sizeof(long)) );

	for(ii1=0; ii1<nonzero.size; ii1++){
		if(nonzero.row[ii1].size <= 0) return(ARG_ERROR_RC);

		skyline->index1[3*ii1 + 0] = 3*(nonzero.row[ii1].elem[0].col);
		alloc_size += (long)((3*ii1+0)-3*(nonzero.row[ii1].elem[0].col)+1);

		skyline->index1[3*ii1 + 1] = 3*(nonzero.row[ii1].elem[0].col);
		alloc_size += (long)((3*ii1+1)-3*(nonzero.row[ii1].elem[0].col)+1);

		skyline->index1[3*ii1 + 2] = 3*(nonzero.row[ii1].elem[0].col);
		alloc_size += (long)((3*ii1+2)-3*(nonzero.row[ii1].elem[0].col)+1);
		if(allocatable_size > 0){
			if(alloc_size > (allocatable_size/(long)sizeof(double))){
				return(SPECIAL_RC);
			}
		}
	}
	RC_TRY( s_cholesky_alloc(skyline->dof, skyline->index1, &(skyline->index2),
	                         &(skyline->array), &(skyline->array_size)) );

	for(ii1=0; ii1<nonzero.size; ii1++){
		for(ii2=0; ii2<nonzero.row[ii1].size; ii2++){
			int row0 = 3*ii1;
			int col0 = 3*(nonzero.row[ii1].elem[ii2].col);

			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					int row = row0 + ii3;
					int col = col0 + ii4;

					if(row < col) continue;

					skyline->array[skyline->index2[row] - row + col]
					        = nonzero.row[ii1].elem[ii2].v[ii3][ii4];
				}
			}
		}
	}

	return(NORMAL_RC);
}


