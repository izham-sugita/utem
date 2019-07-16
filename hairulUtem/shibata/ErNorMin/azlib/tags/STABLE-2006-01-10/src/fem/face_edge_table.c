/*********************************************************************
 * face_edge_table.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Tomoyuki TSUBATA>
 *
 * Thanks to following contributors.
 *  <Yasuo ASAGA> <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: face_edge_table.c 431 2005-08-24 02:59:22Z sasaoka $ */

#include <stdio.h>
#include <stdlib.h>
#include "rc.h"
#include "fem_struct.h"

/* local functions */
static RC set_face_table_tri1(FACE_TABLE *f_table);
static RC set_face_table_tri2(FACE_TABLE *f_table);
static RC set_face_table_quad1(FACE_TABLE *f_table);
static RC set_face_table_quad2(FACE_TABLE *f_table);
static RC set_face_table_tetra1(FACE_TABLE *f_table);
static RC set_face_table_tetra2(FACE_TABLE *f_table);
static RC set_face_table_penta1(FACE_TABLE *f_table);
static RC set_face_table_penta2(FACE_TABLE *f_table);
static RC set_face_table_hexa1(FACE_TABLE *f_table);
static RC set_face_table_hexa2(FACE_TABLE *f_table);
static RC set_edge_table_tri1(EDGE_TABLE *e_table);
static RC set_edge_table_tri2(EDGE_TABLE *e_table);
static RC set_edge_table_quad1(EDGE_TABLE *e_table);
static RC set_edge_table_quad2(EDGE_TABLE *e_table);
static RC set_edge_table_tetra1(EDGE_TABLE *e_table);
static RC set_edge_table_tetra2(EDGE_TABLE *e_table);
static RC set_edge_table_penta1(EDGE_TABLE *e_table);
static RC set_edge_table_penta2(EDGE_TABLE *e_table);
static RC set_edge_table_hexa1(EDGE_TABLE *e_table);
static RC set_edge_table_hexa2(EDGE_TABLE *e_table);


void
init_face_table (FACE_TABLE *f_table)
{
	int ii1, ii2;

	f_table->face_num = 0;
	for(ii1=0; ii1<MAX_FACE; ii1++){
		f_table->node_num[ii1] = 0;
		f_table->face_type[ii1] = ELEM_LINE1;
		for(ii2=0; ii2<MAX_FACE_NODE; ii2++){
			f_table->table[ii1][ii2] = 0;
		}
	}
}


RC
make_face_table (ELEM_TYPE type, FACE_TABLE *f_table)
{
	init_face_table(f_table);

	switch(type){
	case ELEM_TRI1:
		RC_TRY( set_face_table_tri1(f_table) );
		break;
	case ELEM_TRI2:
		RC_TRY( set_face_table_tri2(f_table) );
		break;
	case ELEM_QUAD1:
		RC_TRY( set_face_table_quad1(f_table) );
		break;
	case ELEM_QUAD2:
		RC_TRY( set_face_table_quad2(f_table) );
		break;
	case ELEM_TETRA1:
		RC_TRY( set_face_table_tetra1(f_table) );
		break;
	case ELEM_TETRA2:
		RC_TRY( set_face_table_tetra2(f_table) );
		break;
	case ELEM_PENTA1:
		RC_TRY( set_face_table_penta1(f_table) );
		break;
	case ELEM_PENTA2:
		RC_TRY( set_face_table_penta2(f_table) );
		break;
	case ELEM_HEXA1:
		RC_TRY( set_face_table_hexa1(f_table) );
		break;
	case ELEM_HEXA2:
		RC_TRY( set_face_table_hexa2(f_table) );
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
set_face_table_tri1 (FACE_TABLE *f_table)
{
	f_table->face_num = 3;

	f_table->node_num[0] = 2;
	f_table->face_type[0] = ELEM_LINE1;
	f_table->table[0][0] = 0;
	f_table->table[0][1] = 1;

	f_table->node_num[1] = 2;
	f_table->face_type[1] = ELEM_LINE1;
	f_table->table[1][0] = 1;
	f_table->table[1][1] = 2;

	f_table->node_num[2] = 2;
	f_table->face_type[2] = ELEM_LINE1;
	f_table->table[2][0] = 2;
	f_table->table[2][1] = 0;

	return(NORMAL_RC);
}


static RC
set_face_table_tri2 (FACE_TABLE *f_table)
{
	f_table->face_num = 3;

	f_table->node_num[0] = 3;
	f_table->face_type[0] = ELEM_LINE2;
	f_table->table[0][0] = 0;
	f_table->table[0][1] = 1;
	f_table->table[0][2] = 5;

	f_table->node_num[1] = 3;
	f_table->face_type[1] = ELEM_LINE2;
	f_table->table[1][0] = 1;
	f_table->table[1][1] = 2;
	f_table->table[1][2] = 3;

	f_table->node_num[2] = 3;
	f_table->face_type[2] = ELEM_LINE2;
	f_table->table[2][0] = 2;
	f_table->table[2][1] = 0;
	f_table->table[2][2] = 4;

	return(NORMAL_RC);
}


static RC
set_face_table_quad1 (FACE_TABLE *f_table)
{
	f_table->face_num = 4;

	f_table->node_num[0] = 2;
	f_table->face_type[0] = ELEM_LINE1;
	f_table->table[0][0] = 0;
	f_table->table[0][1] = 1;

	f_table->node_num[1] = 2;
	f_table->face_type[1] = ELEM_LINE1;
	f_table->table[1][0] = 1;
	f_table->table[1][1] = 2;

	f_table->node_num[2] = 2;
	f_table->face_type[2] = ELEM_LINE1;
	f_table->table[2][0] = 2;
	f_table->table[2][1] = 3;

	f_table->node_num[3] = 2;
	f_table->face_type[3] = ELEM_LINE1;
	f_table->table[3][0] = 3;
	f_table->table[3][1] = 0;

	return(NORMAL_RC);
}


static RC
set_face_table_quad2 (FACE_TABLE *f_table)
{
	f_table->face_num = 4;

	f_table->node_num[0] = 3;
	f_table->face_type[0] = ELEM_LINE2;
	f_table->table[0][0] = 0;
	f_table->table[0][1] = 1;
	f_table->table[0][2] = 4;

	f_table->node_num[1] = 3;
	f_table->face_type[1] = ELEM_LINE2;
	f_table->table[1][0] = 1;
	f_table->table[1][1] = 2;
	f_table->table[1][2] = 5;

	f_table->node_num[2] = 3;
	f_table->face_type[2] = ELEM_LINE2;
	f_table->table[2][0] = 2;
	f_table->table[2][1] = 3;
	f_table->table[2][2] = 6;

	f_table->node_num[3] = 3;
	f_table->face_type[3] = ELEM_LINE2;
	f_table->table[3][0] = 3;
	f_table->table[3][1] = 0;
	f_table->table[3][2] = 7;

	return(NORMAL_RC);
}


static RC
set_face_table_tetra1 (FACE_TABLE *f_table)
{
	f_table->face_num = 4;

	f_table->node_num[0] = 3;
	f_table->face_type[0] = ELEM_TRI1;
	f_table->table[0][0] = 0;
	f_table->table[0][1] = 3;
	f_table->table[0][2] = 1;

	f_table->node_num[1] = 3;
	f_table->face_type[1] = ELEM_TRI1;
	f_table->table[1][0] = 1;
	f_table->table[1][1] = 3;
	f_table->table[1][2] = 2;

	f_table->node_num[2] = 3;
	f_table->face_type[2] = ELEM_TRI1;
	f_table->table[2][0] = 0;
	f_table->table[2][1] = 1;
	f_table->table[2][2] = 2;

	f_table->node_num[3] = 3;
	f_table->face_type[3] = ELEM_TRI1;
	f_table->table[3][0] = 2;
	f_table->table[3][1] = 3;
	f_table->table[3][2] = 0;

	return(NORMAL_RC);
}


static RC
set_face_table_tetra2 (FACE_TABLE *f_table)
{
	f_table->face_num = 4;

	f_table->node_num[0] = 6;
	f_table->face_type[0] = ELEM_TRI2;
	f_table->table[0][0] = 0;
	f_table->table[0][1] = 3;
	f_table->table[0][2] = 1;
	f_table->table[0][3] = 8;
	f_table->table[0][4] = 6;
	f_table->table[0][5] = 7;

	f_table->node_num[1] = 6;
	f_table->face_type[1] = ELEM_TRI2;
	f_table->table[1][0] = 1;
	f_table->table[1][1] = 3;
	f_table->table[1][2] = 2;
	f_table->table[1][3] = 9;
	f_table->table[1][4] = 4;
	f_table->table[1][5] = 8;

	f_table->node_num[2] = 6;
	f_table->face_type[2] = ELEM_TRI2;
	f_table->table[2][0] = 0;
	f_table->table[2][1] = 1;
	f_table->table[2][2] = 2;
	f_table->table[2][3] = 4;
	f_table->table[2][4] = 5;
	f_table->table[2][5] = 6;

	f_table->node_num[3] = 6;
	f_table->face_type[3] = ELEM_TRI2;
	f_table->table[3][0] = 2;
	f_table->table[3][1] = 3;
	f_table->table[3][2] = 0;
	f_table->table[3][3] = 7;
	f_table->table[3][4] = 5;
	f_table->table[3][5] = 9;

	return(NORMAL_RC);
}


static RC
set_face_table_penta1 (FACE_TABLE *f_table)
{
	f_table->face_num = 5;

	f_table->node_num[0] = 3;
	f_table->face_type[0] = ELEM_TRI1;
	f_table->table[0][0] = 1;
	f_table->table[0][1] = 0;
	f_table->table[0][2] = 2;

	f_table->node_num[1] = 3;
	f_table->face_type[1] = ELEM_TRI1;
	f_table->table[1][0] = 3;
	f_table->table[1][1] = 4;
	f_table->table[1][2] = 5;

	f_table->node_num[2] = 4;
	f_table->face_type[2] = ELEM_QUAD1;
	f_table->table[2][0] = 1;
	f_table->table[2][1] = 2;
	f_table->table[2][2] = 5;
	f_table->table[2][3] = 4;

	f_table->node_num[3] = 4;
	f_table->face_type[3] = ELEM_QUAD1;
	f_table->table[3][0] = 0;
	f_table->table[3][1] = 3;
	f_table->table[3][2] = 5;
	f_table->table[3][3] = 2;

	f_table->node_num[4] = 4;
	f_table->face_type[4] = ELEM_QUAD1;
	f_table->table[4][0] = 0;
	f_table->table[4][1] = 1;
	f_table->table[4][2] = 4;
	f_table->table[4][3] = 3;

	return(NORMAL_RC);
}


static RC
set_face_table_penta2 (FACE_TABLE *f_table)
{
	f_table->face_num = 5;

	f_table->node_num[0] = 6;
	f_table->face_type[0] = ELEM_TRI2;
	f_table->table[0][0] = 1;
	f_table->table[0][1] = 0;
	f_table->table[0][2] = 2;
	f_table->table[0][3] = 7;
	f_table->table[0][4] = 6;
	f_table->table[0][5] = 8;

	f_table->node_num[1] = 6;
	f_table->face_type[1] = ELEM_TRI2;
	f_table->table[1][0] = 3;
	f_table->table[1][1] = 4;
	f_table->table[1][2] = 5;
	f_table->table[1][3] = 9;
	f_table->table[1][4] = 10;
	f_table->table[1][5] = 11;

	f_table->node_num[2] = 8;
	f_table->face_type[2] = ELEM_QUAD2;
	f_table->table[2][0] = 1;
	f_table->table[2][1] = 2;
	f_table->table[2][2] = 5;
	f_table->table[2][3] = 4;
	f_table->table[2][4] = 6;
	f_table->table[2][5] = 14;
	f_table->table[2][6] = 9;
	f_table->table[2][7] = 13;

	f_table->node_num[3] = 8;
	f_table->face_type[3] = ELEM_QUAD2;
	f_table->table[3][0] = 0;
	f_table->table[3][1] = 3;
	f_table->table[3][2] = 5;
	f_table->table[3][3] = 2;
	f_table->table[3][4] = 12;
	f_table->table[3][5] = 10;
	f_table->table[3][6] = 14;
	f_table->table[3][7] = 7;

	f_table->node_num[4] = 8;
	f_table->face_type[4] = ELEM_QUAD2;
	f_table->table[4][0] = 0;
	f_table->table[4][1] = 1;
	f_table->table[4][2] = 4;
	f_table->table[4][3] = 3;
	f_table->table[4][4] = 8;
	f_table->table[4][5] = 13;
	f_table->table[4][6] = 11;
	f_table->table[4][7] = 12;

	return(NORMAL_RC);
}


static RC
set_face_table_hexa1 (FACE_TABLE *f_table)
{
	f_table->face_num = 6;

	f_table->node_num[0] = 4;
	f_table->face_type[0] = ELEM_QUAD1;
	f_table->table[0][0] = 0;
	f_table->table[0][1] = 3;
	f_table->table[0][2] = 2;
	f_table->table[0][3] = 1;

	f_table->node_num[1] = 4;
	f_table->face_type[1] = ELEM_QUAD1;
	f_table->table[1][0] = 0;
	f_table->table[1][1] = 1;
	f_table->table[1][2] = 5;
	f_table->table[1][3] = 4;

	f_table->node_num[2] = 4;
	f_table->face_type[2] = ELEM_QUAD1;
	f_table->table[2][0] = 1;
	f_table->table[2][1] = 2;
	f_table->table[2][2] = 6;
	f_table->table[2][3] = 5;

	f_table->node_num[3] = 4;
	f_table->face_type[3] = ELEM_QUAD1;
	f_table->table[3][0] = 7;
	f_table->table[3][1] = 4;
	f_table->table[3][2] = 5;
	f_table->table[3][3] = 6;

	f_table->node_num[4] = 4;
	f_table->face_type[4] = ELEM_QUAD1;
	f_table->table[4][0] = 3;
	f_table->table[4][1] = 0;
	f_table->table[4][2] = 4;
	f_table->table[4][3] = 7;

	f_table->node_num[5] = 4;
	f_table->face_type[5] = ELEM_QUAD1;
	f_table->table[5][0] = 2;
	f_table->table[5][1] = 3;
	f_table->table[5][2] = 7;
	f_table->table[5][3] = 6;

	return(NORMAL_RC);
}


static RC
set_face_table_hexa2 (FACE_TABLE *f_table)
{
	f_table->face_num = 6;

	f_table->node_num[0] = 8;
	f_table->face_type[0] = ELEM_QUAD2;
	f_table->table[0][0] = 0;
	f_table->table[0][1] = 3;
	f_table->table[0][2] = 2;
	f_table->table[0][3] = 1;
	f_table->table[0][4] = 11;
	f_table->table[0][5] = 10;
	f_table->table[0][6] = 9;
	f_table->table[0][7] = 8;

	f_table->node_num[1] = 8;
	f_table->face_type[1] = ELEM_QUAD2;
	f_table->table[1][0] = 0;
	f_table->table[1][1] = 1;
	f_table->table[1][2] = 5;
	f_table->table[1][3] = 4;
	f_table->table[1][4] = 8;
	f_table->table[1][5] = 17;
	f_table->table[1][6] = 12;
	f_table->table[1][7] = 16;

	f_table->node_num[2] = 8;
	f_table->face_type[2] = ELEM_QUAD2;
	f_table->table[2][0] = 1;
	f_table->table[2][1] = 2;
	f_table->table[2][2] = 6;
	f_table->table[2][3] = 5;
	f_table->table[2][4] = 9;
	f_table->table[2][5] = 18;
	f_table->table[2][6] = 13;
	f_table->table[2][7] = 17;

	f_table->node_num[3] = 8;
	f_table->face_type[3] = ELEM_QUAD2;
	f_table->table[3][0] = 7;
	f_table->table[3][1] = 4;
	f_table->table[3][2] = 5;
	f_table->table[3][3] = 6;
	f_table->table[3][4] = 15;
	f_table->table[3][5] = 12;
	f_table->table[3][6] = 13;
	f_table->table[3][7] = 14;

	f_table->node_num[4] = 8;
	f_table->face_type[4] = ELEM_QUAD2;
	f_table->table[4][0] = 3;
	f_table->table[4][1] = 0;
	f_table->table[4][2] = 4;
	f_table->table[4][3] = 7;
	f_table->table[4][4] = 11;
	f_table->table[4][5] = 16;
	f_table->table[4][6] = 15;
	f_table->table[4][7] = 19;

	f_table->node_num[5] = 8;
	f_table->face_type[5] = ELEM_QUAD2;
	f_table->table[5][0] = 2;
	f_table->table[5][1] = 3;
	f_table->table[5][2] = 7;
	f_table->table[5][3] = 6;
	f_table->table[5][4] = 10;
	f_table->table[5][5] = 19;
	f_table->table[5][6] = 14;
	f_table->table[5][7] = 18;

	return(NORMAL_RC);
}


void
init_edge_table (EDGE_TABLE *e_table)
{
	int ii1, ii2;

	e_table->edge_num = 0;
	for(ii1=0; ii1<MAX_EDGE; ii1++){
		e_table->node_num[ii1] = 0;
		e_table->edge_type[ii1] = ELEM_LINE1;
		for(ii2=0; ii2<MAX_EDGE_NODE; ii2++){
			e_table->table[ii1][ii2] = 0;
		}
	}
}


RC
make_edge_table (ELEM_TYPE type, EDGE_TABLE *e_table)
{
	init_edge_table(e_table);

	switch(type){
	case ELEM_TRI1:
		RC_TRY( set_edge_table_tri1(e_table) );
		break;
	case ELEM_TRI2:
		RC_TRY( set_edge_table_tri2(e_table) );
		break;
	case ELEM_QUAD1:
		RC_TRY( set_edge_table_quad1(e_table) );
		break;
	case ELEM_QUAD2:
		RC_TRY( set_edge_table_quad2(e_table) );
		break;
	case ELEM_TETRA1:
		RC_TRY( set_edge_table_tetra1(e_table) );
		break;
	case ELEM_TETRA2:
		RC_TRY( set_edge_table_tetra2(e_table) );
		break;
	case ELEM_PENTA1:
		RC_TRY( set_edge_table_penta1(e_table) );
		break;
	case ELEM_PENTA2:
		RC_TRY( set_edge_table_penta2(e_table) );
		break;
	case ELEM_HEXA1:
		RC_TRY( set_edge_table_hexa1(e_table) );
		break;
	case ELEM_HEXA2:
		RC_TRY( set_edge_table_hexa2(e_table) );
		break;
	default:
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC
set_edge_table_tri1 (EDGE_TABLE *e_table)
{
	e_table->edge_num = 3;

	e_table->node_num[0] = 2;
	e_table->edge_type[0] = ELEM_LINE1;
	e_table->table[0][0] = 1;
	e_table->table[0][1] = 2;

	e_table->node_num[1] = 2;
	e_table->edge_type[1] = ELEM_LINE1;
	e_table->table[1][0] = 2;
	e_table->table[1][1] = 0;

	e_table->node_num[2] = 2;
	e_table->edge_type[2] = ELEM_LINE1;
	e_table->table[2][0] = 0;
	e_table->table[2][1] = 1;

	return(NORMAL_RC);
}


static RC
set_edge_table_tri2 (EDGE_TABLE *e_table)
{
	e_table->edge_num = 3;

	e_table->node_num[0] = 3;
	e_table->edge_type[0] = ELEM_LINE2;
	e_table->table[0][0] = 1;
	e_table->table[0][1] = 2;
	e_table->table[0][2] = 3;

	e_table->node_num[1] = 3;
	e_table->edge_type[1] = ELEM_LINE2;
	e_table->table[1][0] = 2;
	e_table->table[1][1] = 0;
	e_table->table[1][2] = 4;

	e_table->node_num[2] = 3;
	e_table->edge_type[2] = ELEM_LINE2;
	e_table->table[2][0] = 0;
	e_table->table[2][1] = 1;
	e_table->table[2][2] = 5;

	return(NORMAL_RC);
}


static RC
set_edge_table_quad1 (EDGE_TABLE *e_table)
{
	e_table->edge_num = 4;

	e_table->node_num[0] = 2;
	e_table->edge_type[0] = ELEM_LINE1;
	e_table->table[0][0] = 0;
	e_table->table[0][1] = 1;

	e_table->node_num[1] = 2;
	e_table->edge_type[1] = ELEM_LINE1;
	e_table->table[1][0] = 1;
	e_table->table[1][1] = 2;

	e_table->node_num[2] = 2;
	e_table->edge_type[2] = ELEM_LINE1;
	e_table->table[2][0] = 2;
	e_table->table[2][1] = 3;

	e_table->node_num[3] = 2;
	e_table->edge_type[3] = ELEM_LINE1;
	e_table->table[3][0] = 3;
	e_table->table[3][1] = 0;

	return(NORMAL_RC);
}


static RC
set_edge_table_quad2 (EDGE_TABLE *e_table)
{
	e_table->edge_num = 4;

	e_table->node_num[0] = 3;
	e_table->edge_type[0] = ELEM_LINE2;
	e_table->table[0][0] = 0;
	e_table->table[0][1] = 1;
	e_table->table[0][2] = 4;

	e_table->node_num[1] = 3;
	e_table->edge_type[1] = ELEM_LINE2;
	e_table->table[1][0] = 1;
	e_table->table[1][1] = 2;
	e_table->table[1][2] = 5;

	e_table->node_num[2] = 3;
	e_table->edge_type[2] = ELEM_LINE2;
	e_table->table[2][0] = 2;
	e_table->table[2][1] = 3;
	e_table->table[2][2] = 6;

	e_table->node_num[3] = 3;
	e_table->edge_type[3] = ELEM_LINE2;
	e_table->table[3][0] = 3;
	e_table->table[3][1] = 0;
	e_table->table[3][2] = 7;

	return(NORMAL_RC);
}


static RC
set_edge_table_tetra1 (EDGE_TABLE *e_table)
{
	e_table->edge_num = 6;

	e_table->node_num[0] = 2;
	e_table->edge_type[0] = ELEM_LINE1;
	e_table->table[0][0] = 1;
	e_table->table[0][1] = 2;

	e_table->node_num[1] = 2;
	e_table->edge_type[1] = ELEM_LINE1;
	e_table->table[1][0] = 2;
	e_table->table[1][1] = 0;

	e_table->node_num[2] = 2;
	e_table->edge_type[2] = ELEM_LINE1;
	e_table->table[2][0] = 0;
	e_table->table[2][1] = 1;

	e_table->node_num[3] = 2;
	e_table->edge_type[3] = ELEM_LINE1;
	e_table->table[3][0] = 3;
	e_table->table[3][1] = 0;

	e_table->node_num[4] = 2;
	e_table->edge_type[4] = ELEM_LINE1;
	e_table->table[4][0] = 1;
	e_table->table[4][1] = 3;

	e_table->node_num[5] = 2;
	e_table->edge_type[5] = ELEM_LINE1;
	e_table->table[5][0] = 2;
	e_table->table[5][1] = 3;

	return(NORMAL_RC);
}


static RC
set_edge_table_tetra2 (EDGE_TABLE *e_table)
{
	e_table->edge_num = 6;

	e_table->node_num[0] = 3;
	e_table->edge_type[0] = ELEM_LINE2;
	e_table->table[0][0] = 1;
	e_table->table[0][1] = 2;
	e_table->table[0][2] = 4;

	e_table->node_num[1] = 3;
	e_table->edge_type[1] = ELEM_LINE2;
	e_table->table[1][0] = 2;
	e_table->table[1][1] = 0;
	e_table->table[1][2] = 5;

	e_table->node_num[2] = 3;
	e_table->edge_type[2] = ELEM_LINE2;
	e_table->table[2][0] = 0;
	e_table->table[2][1] = 1;
	e_table->table[2][2] = 6;

	e_table->node_num[3] = 3;
	e_table->edge_type[3] = ELEM_LINE2;
	e_table->table[3][0] = 3;
	e_table->table[3][1] = 0;
	e_table->table[3][2] = 7;

	e_table->node_num[4] = 3;
	e_table->edge_type[4] = ELEM_LINE2;
	e_table->table[4][0] = 1;
	e_table->table[4][1] = 3;
	e_table->table[4][2] = 8;

	e_table->node_num[5] = 3;
	e_table->edge_type[5] = ELEM_LINE2;
	e_table->table[5][0] = 2;
	e_table->table[5][1] = 3;
	e_table->table[5][2] = 9;

	return(NORMAL_RC);
}


static RC
set_edge_table_penta1 (EDGE_TABLE *e_table)
{
	e_table->edge_num = 9;

	e_table->node_num[0] = 2;
	e_table->edge_type[0] = ELEM_LINE1;
	e_table->table[0][0] = 1;
	e_table->table[0][1] = 2;

	e_table->node_num[1] = 2;
	e_table->edge_type[1] = ELEM_LINE1;
	e_table->table[1][0] = 2;
	e_table->table[1][1] = 0;

	e_table->node_num[2] = 2;
	e_table->edge_type[2] = ELEM_LINE1;
	e_table->table[2][0] = 0;
	e_table->table[2][1] = 1;

	e_table->node_num[3] = 2;
	e_table->edge_type[3] = ELEM_LINE1;
	e_table->table[3][0] = 4;
	e_table->table[3][1] = 5;

	e_table->node_num[4] = 2;
	e_table->edge_type[4] = ELEM_LINE1;
	e_table->table[4][0] = 5;
	e_table->table[4][1] = 3;

	e_table->node_num[5] = 2;
	e_table->edge_type[5] = ELEM_LINE1;
	e_table->table[5][0] = 3;
	e_table->table[5][1] = 4;

	e_table->node_num[6] = 2;
	e_table->edge_type[6] = ELEM_LINE1;
	e_table->table[6][0] = 0;
	e_table->table[6][1] = 3;

	e_table->node_num[7] = 2;
	e_table->edge_type[7] = ELEM_LINE1;
	e_table->table[7][0] = 1;
	e_table->table[7][1] = 4;

	e_table->node_num[8] = 2;
	e_table->edge_type[8] = ELEM_LINE1;
	e_table->table[8][0] = 2;
	e_table->table[8][1] = 5;

	return(NORMAL_RC);
}


static RC
set_edge_table_penta2 (EDGE_TABLE *e_table)
{
	e_table->edge_num = 9;

	e_table->node_num[0] = 3;
	e_table->edge_type[0] = ELEM_LINE2;
	e_table->table[0][0] = 1;
	e_table->table[0][1] = 2;
	e_table->table[0][2] = 6;

	e_table->node_num[1] = 3;
	e_table->edge_type[1] = ELEM_LINE2;
	e_table->table[1][0] = 2;
	e_table->table[1][1] = 0;
	e_table->table[1][2] = 7;

	e_table->node_num[2] = 3;
	e_table->edge_type[2] = ELEM_LINE2;
	e_table->table[2][0] = 0;
	e_table->table[2][1] = 1;
	e_table->table[2][2] = 8;

	e_table->node_num[3] = 3;
	e_table->edge_type[3] = ELEM_LINE2;
	e_table->table[3][0] = 4;
	e_table->table[3][1] = 5;
	e_table->table[3][2] = 9;

	e_table->node_num[4] = 3;
	e_table->edge_type[4] = ELEM_LINE2;
	e_table->table[4][0] = 5;
	e_table->table[4][1] = 3;
	e_table->table[4][2] = 10;

	e_table->node_num[5] = 3;
	e_table->edge_type[5] = ELEM_LINE2;
	e_table->table[5][0] = 3;
	e_table->table[5][1] = 4;
	e_table->table[5][2] = 11;

	e_table->node_num[6] = 3;
	e_table->edge_type[6] = ELEM_LINE2;
	e_table->table[6][0] = 0;
	e_table->table[6][1] = 3;
	e_table->table[6][2] = 12;

	e_table->node_num[7] = 3;
	e_table->edge_type[7] = ELEM_LINE2;
	e_table->table[7][0] = 1;
	e_table->table[7][1] = 4;
	e_table->table[7][2] = 13;

	e_table->node_num[8] = 3;
	e_table->edge_type[8] = ELEM_LINE2;
	e_table->table[8][0] = 2;
	e_table->table[8][1] = 5;
	e_table->table[8][2] = 14;

	return(NORMAL_RC);
}


static RC
set_edge_table_hexa1 (EDGE_TABLE *e_table)
{
	e_table->edge_num = 12;

	e_table->node_num[0] = 2;
	e_table->edge_type[0] = ELEM_LINE1;
	e_table->table[0][0] = 0;
	e_table->table[0][1] = 1;

	e_table->node_num[1] = 2;
	e_table->edge_type[1] = ELEM_LINE1;
	e_table->table[1][0] = 1;
	e_table->table[1][1] = 2;

	e_table->node_num[2] = 2;
	e_table->edge_type[2] = ELEM_LINE1;
	e_table->table[2][0] = 2;
	e_table->table[2][1] = 3;

	e_table->node_num[3] = 2;
	e_table->edge_type[3] = ELEM_LINE1;
	e_table->table[3][0] = 3;
	e_table->table[3][1] = 0;

	e_table->node_num[4] = 2;
	e_table->edge_type[4] = ELEM_LINE1;
	e_table->table[4][0] = 4;
	e_table->table[4][1] = 5;

	e_table->node_num[5] = 2;
	e_table->edge_type[5] = ELEM_LINE1;
	e_table->table[5][0] = 5;
	e_table->table[5][1] = 6;

	e_table->node_num[6] = 2;
	e_table->edge_type[6] = ELEM_LINE1;
	e_table->table[6][0] = 6;
	e_table->table[6][1] = 7;

	e_table->node_num[7] = 2;
	e_table->edge_type[7] = ELEM_LINE1;
	e_table->table[7][0] = 7;
	e_table->table[7][1] = 4;

	e_table->node_num[8] = 2;
	e_table->edge_type[8] = ELEM_LINE1;
	e_table->table[8][0] = 0;
	e_table->table[8][1] = 4;

	e_table->node_num[9] = 2;
	e_table->edge_type[9] = ELEM_LINE1;
	e_table->table[9][0] = 1;
	e_table->table[9][1] = 5;

	e_table->node_num[10] = 2;
	e_table->edge_type[10] = ELEM_LINE1;
	e_table->table[10][0] = 2;
	e_table->table[10][1] = 6;

	e_table->node_num[11] = 2;
	e_table->edge_type[11] = ELEM_LINE1;
	e_table->table[11][0] = 3;
	e_table->table[11][1] = 7;

	return(NORMAL_RC);
}


static RC
set_edge_table_hexa2 (EDGE_TABLE *e_table)
{
	e_table->edge_num = 12;

	e_table->node_num[0] = 3;
	e_table->edge_type[0] = ELEM_LINE2;
	e_table->table[0][0] = 0;
	e_table->table[0][1] = 1;
	e_table->table[0][2] = 8;

	e_table->node_num[1] = 3;
	e_table->edge_type[1] = ELEM_LINE2;
	e_table->table[1][0] = 1;
	e_table->table[1][1] = 2;
	e_table->table[1][2] = 9;

	e_table->node_num[2] = 3;
	e_table->edge_type[2] = ELEM_LINE2;
	e_table->table[2][0] = 2;
	e_table->table[2][1] = 3;
	e_table->table[2][2] = 10;

	e_table->node_num[3] = 3;
	e_table->edge_type[3] = ELEM_LINE2;
	e_table->table[3][0] = 3;
	e_table->table[3][1] = 0;
	e_table->table[3][2] = 11;

	e_table->node_num[4] = 3;
	e_table->edge_type[4] = ELEM_LINE2;
	e_table->table[4][0] = 4;
	e_table->table[4][1] = 5;
	e_table->table[4][2] = 12;

	e_table->node_num[5] = 3;
	e_table->edge_type[5] = ELEM_LINE2;
	e_table->table[5][0] = 5;
	e_table->table[5][1] = 6;
	e_table->table[5][2] = 13;

	e_table->node_num[6] = 3;
	e_table->edge_type[6] = ELEM_LINE2;
	e_table->table[6][0] = 6;
	e_table->table[6][1] = 7;
	e_table->table[6][2] = 14;

	e_table->node_num[7] = 3;
	e_table->edge_type[7] = ELEM_LINE2;
	e_table->table[7][0] = 7;
	e_table->table[7][1] = 4;
	e_table->table[7][2] = 15;

	e_table->node_num[8] = 3;
	e_table->edge_type[8] = ELEM_LINE2;
	e_table->table[8][0] = 0;
	e_table->table[8][1] = 4;
	e_table->table[8][2] = 16;

	e_table->node_num[9] = 3;
	e_table->edge_type[9] = ELEM_LINE2;
	e_table->table[9][0] = 1;
	e_table->table[9][1] = 5;
	e_table->table[9][2] = 17;

	e_table->node_num[10] = 3;
	e_table->edge_type[10] = ELEM_LINE2;
	e_table->table[10][0] = 2;
	e_table->table[10][1] = 6;
	e_table->table[10][2] = 18;

	e_table->node_num[11] = 3;
	e_table->edge_type[11] = ELEM_LINE2;
	e_table->table[11][0] = 3;
	e_table->table[11][1] = 7;
	e_table->table[11][2] = 19;

	return(NORMAL_RC);
}

