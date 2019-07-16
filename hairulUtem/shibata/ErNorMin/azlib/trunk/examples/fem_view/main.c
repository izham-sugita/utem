#include <stdio.h>
#include <stdlib.h>
#include <GL/gl.h>
#include <GL/glut.h>
#include "base.h"
#include "mathematics.h"
#include "fem.h"
#include "glut_utl.h"

int main(int argc, char **argv);
RC input_bdf(FILE *fp, NODE_ARRAY *node, ELEMENT_ARRAY *element);

int main (int argc, char **argv)
{
	FILE *fp;
	NODE_ARRAY node;
	ELEMENT_ARRAY element, surf;

	if(argc < 2){
		fprintf(stderr, "Usage : %s [r:fem_model(*.bdf)]\n", argv[0]);
		return(EXIT_FAILURE);
	}


	/* 使用メモリの上限値設定 */
	RC_TRY_MAIN( mm_init((size_t)512*MEGA_BYTE) );

	/* モデルデータ読み込み */
	RC_TRY_MAIN( rc_fopen(argv[1], "r", &fp) );
	RC_TRY_MAIN( input_bdf(fp, &node, &element) );
	RC_TRY_MAIN( rc_fclose(fp) );

	/* モデル表示 */
	glutInit(&argc, argv);
	if(analysis_dim(element) == 3){
		RC_TRY_MAIN( extract_surface(element, &surf) );
		RC_TRY_MAIN( reduced_element_array_order(&surf) );
		RC_TRY_MAIN( show_fem_model(surf, node) );
	}else{
		RC_TRY_MAIN( reduced_element_array_order(&element) );
		RC_TRY_MAIN( show_fem_model(element, node) );
	}

	RC_TRY_MAIN( mm_terminate() );

	return(EXIT_SUCCESS);
}


RC input_bdf (FILE *fp, NODE_ARRAY *node, ELEMENT_ARRAY *element)
{
	/* 節点 */
	RC_TRY( nst_input_node(fp, node) );
	fprintf(stdout, "Number of nodes    : %d\n", node->size);

	/* 要素 */
	RC_TRY( nst_input_element(fp, element) );
	fprintf(stdout, "Number of elements : %d\n", element->size);

	return(NORMAL_RC);
}



