#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "base.h"
#include "mathematics.h"
#include "fem.h"
#include "matrix_db.h"

int
main(int argc, char **argv)
{
	int ii1;
	FILE *fp;
	ELEMENT_ARRAY elem1, elem2, surf1, surf2;
	ELEMENT_ARRAY new2;
	NODE_ARRAY node1, node2;
	int size1, size2;
	int *index_array1, *index_array2;
	GEOM_SEGMENT *segment[3];
	int seg_size;
	int *cand;
	int *search_elem_index;
	int elem_index_size;

	if(argc < 3){
		fprintf(stderr, "Usage : %s [r:model1(.NAS)] [r:model2(.NAS)]\n",
		        argv[0]);
		return(EXIT_FAILURE);
	}

	/* Log Print Level */
	RC_TRY_MAIN( set_log_level(5) );
	RC_TRY_MAIN( set_log_file(0, stderr) )

	/* Memory Manager */
	RC_TRY_MAIN( mm_init((size_t)1024*MEGA_BYTE) );

	/* target */
	RC_TRY_MAIN( rc_fopen(argv[1], "r", &fp) );
	RC_TRY_MAIN( nst_input_node(fp, &node1) );
	RC_TRY_MAIN( nst_input_element(fp, &elem1) );
	RC_TRY_MAIN( rc_fclose(fp) );
	RC_TRY_MAIN( dummy_renumber(&node1) );
	log_printf(5, "node1.size = %d\n", node1.size);
	log_printf(5, "elem1.size = %d\n", elem1.size);


	/* hit */
	RC_TRY_MAIN( rc_fopen(argv[2], "r", &fp) );
	RC_TRY_MAIN( nst_input_node(fp, &node2) );
	RC_TRY_MAIN( nst_input_element(fp, &elem2) );
	RC_TRY_MAIN( rc_fclose(fp) );
	RC_TRY_MAIN( dummy_renumber(&node2) );
	log_printf(5, "node2.size = %d\n", node2.size);
	log_printf(5, "elem2.size = %d\n", elem2.size);
#if 1
	{
		for(ii1=0; ii1<node2.size; ii1++){
			node2.array[ii1].p.y -= 0.001f; /* キューブ */
			//node2.array[ii1].p.x += 200.0f; /* Knuckle */
		}
	}
#endif

	DEBUG_PRINT;
	RC_TRY_MAIN( extract_surface(elem1, &surf1) );
	RC_TRY_MAIN( extract_surface(elem2, &surf2) );
	log_printf(5, "surf1.size = %d\n", surf1.size);
	log_printf(5, "surf2.size = %d\n", surf2.size);

	DEBUG_PRINT;
	RC_NULL_CHK( index_array1 = mm_alloc(surf1.size*sizeof(int)) );
	RC_NULL_CHK( index_array2 = mm_alloc(surf2.size*sizeof(int)) );
	DEBUG_PRINT;
	RC_TRY_MAIN( extract_geom_contact_surface(surf1, node1, surf2, node2,
	                                          &size1, index_array1,
	                                          &size2, index_array2) );
	RC_TRY_MAIN( log_printf(5, "size1 =  %d\n", size1) );
	RC_TRY_MAIN( log_printf(5, "size2 =  %d\n", size2) );
	DEBUG_PRINT;

	RC_TRY_MAIN( allocate_element_array(size2, &new2) );
	for(ii1=0; ii1<size2; ii1++){
		log_printf(5, "index_array2[%d] = %d\n", ii1, index_array2[ii1]);
		new2.array[ii1] = surf2.array[index_array2[ii1]];
	}
	new2.size = size2;
	RC_TRY_MAIN( clean_element_array(&new2) );
	print_element_array(stderr, new2);
	//print_node_array(stderr, node2);
	print_node(stderr, node2.array[10]);

#if 0
	RC_TRY_MAIN( make_geom_search_array(new2, node2, &seg_size, segment) );
#else
	RC_TRY_MAIN( make_geom_search_array_2(surf2, node2, size2, index_array2,
	                                      &seg_size, segment) );
#endif
	fprintf(stderr, "seg_size = %d\n", seg_size);
	if(seg_size < 1) return(EXIT_SUCCESS);

	for(ii1=0; ii1<seg_size; ii1++){
		fprintf(stderr, "segment[0][%d].index = %d\n",
		        ii1, segment[0][ii1].index);
		fprintf(stderr, "segment[1][%d].index = %d\n",
		        ii1, segment[1][ii1].index);
		fprintf(stderr, "segment[2][%d].index = %d\n\n",
		        ii1, segment[2][ii1].index);
	}

	/* 披接触側要素のindex配列の整頓 */
	DEBUG_PRINT;
	RC_NULL_CHK( search_elem_index = mm_alloc(size1*sizeof(int)) );
	for(ii1=0, elem_index_size=0; ii1<size1; ii1++){
		int label = surf1.array[index_array1[ii1]].i_info[0];
		int index = search_element_label(elem1, label);
		RC_NEG_CHK( search_elem_index[elem_index_size] = index );
		elem_index_size++;
		//fprintf(stderr, "label = %d\n", label);
	}
	DEBUG_PRINT;
	RC_TRY_MAIN( sort_int_array(elem_index_size, search_elem_index) );
	RC_TRY_MAIN( uniq_int_array(&elem_index_size, &search_elem_index) );
	DEBUG_PRINT;
	fprintf(stderr, "elem_index_size = %d\n", elem_index_size);
	for(ii1=0; ii1<elem_index_size; ii1++){
		fprintf(stderr, "search_elem_index[%d] = %d\n", ii1,
		                 search_elem_index[ii1]);
	}


	DEBUG_PRINT;
	/* 衝突可能性のある相手要素から整頓された節点のindex配列(cand[])を作る */
	RC_NULL_CHK( cand = mm_alloc(seg_size*sizeof(int)) );
	for(ii1=0; ii1<elem_index_size; ii1++){
		int ii2, ii3;
		VECT3D min, max, center;
		int cand_size;
		int node_size;
		int *contact_node_array;
		/* 披接触側要素のAABB */
		RC_TRY_MAIN( element_bbox(elem1.array[search_elem_index[ii1]],
		                          node1, &min, &max) );
		//fprintf(stderr, "[%d]\n", ii1);
		//RC_TRY_MAIN( print_vect3d(stderr, min) );
		//RC_TRY_MAIN( print_vect3d(stderr, max) );
		center = mid_point3d(min, max);
		/* 接触側要素の接触候補要素index */
		RC_TRY_MAIN( geom_search_candidate(seg_size, segment, center,
		                                   dist_point3d(center, min)+ABS_TOL,
		                                   &cand_size, cand) );
		fprintf(stderr, "cand_size = %d\n", cand_size);
		/* 表面要素の最大節点数 = 8 */
		RC_NULL_CHK( contact_node_array = mm_alloc(cand_size*8*sizeof(int)) );
		for(ii2=0, node_size=0; ii2<cand_size; ii2++){
#if 0
			for(ii3=0; ii3<new2.array[cand[ii2]].node_num;ii3++){
				int index = node_renum_index(node2,
				                             new2.array[cand[ii2]].node[ii3]);
				RC_NEG_CHK(index);
				contact_node_array[node_size++] = index;
			}
#else
			for(ii3=0; ii3<surf2.array[cand[ii2]].node_num;ii3++){
				int index = node_renum_index(node2,
				                             surf2.array[cand[ii2]].node[ii3]);
				RC_NEG_CHK(index);
				contact_node_array[node_size++] = index;
			}
#endif
		}
		RC_TRY_MAIN( sort_int_array(node_size, contact_node_array) );
		RC_TRY_MAIN( uniq_int_array(&node_size, &contact_node_array) );
		fprintf(stderr, "node_size = %d\n", node_size);
		for(ii2=0; ii2<node_size; ii2++){
			/* elem1.array[search_elem_index[ii1]]に */
			/* 接触する可能性が高い節点のindex配列   */
			fprintf(stderr, "contact_node_array[%d] = %d\n",
			        ii2, contact_node_array[ii2]);
		}
#if 0
		for(ii2=0; ii2<cand_size; ii2++){
			fprintf(stderr, "cand[%d] = %d\n", ii2, cand[ii2]);
		}
#endif
	}
	RC_TRY_MAIN( mm_free((void *)cand) );
	RC_TRY_MAIN( mm_free((void *)search_elem_index) );


	/*
	CONTACT contact;
	RC_TRY_MAIN( make_contact_terminus(surf1, node1, node2.array[10], 0.1f,
	                                   seg_size, segment, &contact) );
	print_contact(stderr, contact);
	 */

	/*
	RC_TRY_MAIN( get_element_construct_node(&num, &node_label, new2) );
	fprintf(stderr, "num = %d\n", num);
	for(ii1=0; ii1<num; ii1++){
		fprintf(stderr, "node_label[%d] = %d\n", ii1, node_label[ii1]);
	}
	*/

	RC_TRY_MAIN( mm_free((void *)index_array1) );
	RC_TRY_MAIN( mm_free((void *)index_array2) );
	RC_TRY_MAIN( free_node_array(&node1) );
	RC_TRY_MAIN( free_node_array(&node2) );
	RC_TRY_MAIN( free_element_array(&elem1) );
	RC_TRY_MAIN( free_element_array(&elem2) );
	RC_TRY_MAIN( free_element_array(&new2) );


	/* Memory Manager */
	RC_TRY_MAIN( mm_terminate() );

	return(EXIT_SUCCESS);
}

