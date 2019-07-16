/*********************************************************************
 * bmp_utl.c
 *
 * Copyright (C) 2001 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA>
 *
 * Thanks to following contributors.
 *  <Kenzen TAKEUCHI> <Tomoyuki TSUBATA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: bmp_utl.c,v 1.8 2003/11/26 08:10:15 sasaoka Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include "rc.h"
#include "bin_io.h"
#include "graphic_utl.h"

/*
 * $B05=L$7$?(BBMP$B%$%a!<%8$O%5%]!<%H$7$J$$(B
 * Write $B$O(B24bit$B$N$_BP1~(B
 * Read $B$O(B24bit,8bit[palette]$B$N$_BP1~(B
 * (4bit,2bit$B$N2hA|$O%Q%C%G%#%s%0$N=hM}$,ITL@$N$?$aL$BP1~(B)
 */

static RC read_bmp_palette(FILE *fp, GRAPHIC *gr, GRAPHIC_RGB *palette);
static RC read_bmp_image(FILE *fp, GRAPHIC *gr);

RC read_bmp(FILE *fp, GRAPHIC *gr)
{
	init_graphic(gr);
	RC_TRY( read_bmp_header(fp, gr) );
	RC_TRY( read_bmp_image(fp, gr) );
	return(NORMAL_RC);
}

RC read_bmp_header(FILE *fp, GRAPHIC *gr)
{
	unsigned int  ui, ui1, ui2;
	unsigned long ul;

	/* $B%X%C%@<1JLNN0h(B */
	if((ui1 = fgetc(fp)) == EOF) return(READ_ERROR_RC);
	if((ui2 = fgetc(fp)) == EOF) return(READ_ERROR_RC);
	if( (ui1 == 'B') && (ui2 == 'M') ){
		gr->source = GRAPHIC_BMP;
	}else{
		fprintf(stderr, "This is not a BMP file\n");
		return(READ_ERROR_RC);
	}

    /* $B%U%!%$%k%5%$%:(B */
	RC_TRY( lread_4bytes(fp, &ul) );

	/* $BM=LsNN0h(B1 & $BM=LsNN0h(B2 */
	RC_TRY( lread_2bytes(fp, &ui) );
	RC_TRY( lread_2bytes(fp, &ui) );

    /* $B%X%C%@%5%$%:(B */
	RC_TRY( lread_4bytes(fp, &ul) );

    /* $B>pJs%5%$%:(B */
	RC_TRY( lread_4bytes(fp, &ul) );

    /* $B2hA|(BX$B%5%$%:(B */
	RC_TRY( lread_4bytes(fp, &ul) );
	if(ul > 0){
		gr->size.x = ul;
	}else{
		fprintf(stderr, "Wrong image size [X]\n");
		return(READ_ERROR_RC);
	}

    /* $B2hA|(BY$B%5%$%:(B */
	RC_TRY( lread_4bytes(fp, &ul) );
	if(ul > 0){
		gr->size.y = ul;
	}else{
		fprintf(stderr, "Wrong image size [Y]\n");
		return(READ_ERROR_RC);
	}

    /* $B2hA|(BZ$B%5%$%:(B */
	gr->size.z = 1;

	/* $BLL?t(B */
	RC_TRY( lread_2bytes(fp, &ui) );

	/* $B?'%S%C%H?t(B */
	RC_TRY( lread_2bytes(fp, &ui) );
	if( (ui == 8) || (ui == 24) ){
		gr->type = GRAPHIC_RGB24;
		gr->i_info[0] = ui;
	}else if( (ui == 1) || (ui == 4) ){
		fprintf(stderr, "Non-Implement color bit\n");
		return(READ_ERROR_RC);
	}else{
		fprintf(stderr, "Unknown color bit\n");
		return(READ_ERROR_RC);
	}

	/* $B05=LJ}<0(B & $B05=L%5%$%:(B */
	RC_TRY( lread_4bytes(fp, &ul) );
	if(ul != 0){
		fprintf(stderr, "This is Compressed BMP file."
		                "We do not support compressed bmp\n");
		return(READ_ERROR_RC);
	}
	RC_TRY( lread_4bytes(fp, &ul) );

	/* $B?eJ?2rA|EY(B */
	RC_TRY( lread_4bytes(fp, &ul) );
	gr->resolution.x = ul/1000.0;

	/* $B?bD>2rA|EY(B */
	RC_TRY( lread_4bytes(fp, &ul) );
	gr->resolution.y = ul/1000.0;

	/* $B?'?t(B */
	RC_TRY( lread_4bytes(fp, &ul) );

	/* $B=EMW?'?t(B */
	RC_TRY( lread_4bytes(fp, &ul) );

	return(NORMAL_RC);
}


static RC read_bmp_palette(FILE *fp, GRAPHIC *gr, GRAPHIC_RGB *palette)
{
	int ii1, cl;
	for(ii1=0; ii1<(1<<(gr->i_info[0])); ii1++){
		if((cl = fgetc(fp)) == EOF) return(READ_ERROR_RC);
		  palette[ii1].b = (unsigned char)cl;
		if((cl = fgetc(fp)) == EOF) return(READ_ERROR_RC);
		  palette[ii1].g = (unsigned char)cl;
		if((cl = fgetc(fp)) == EOF) return(READ_ERROR_RC);
		  palette[ii1].r = (unsigned char)cl;
		/* padding$BBe$j$N(B4$B8DL\$,$"$k(B */
		if((cl = fgetc(fp)) == EOF) return(READ_ERROR_RC);
	}
	return(NORMAL_RC);
}


static RC read_bmp_image(FILE *fp, GRAPHIC *gr)
{
	int ii1, ii2, ii3;
	int padding = 0;
	int cl;
	
	GRAPHIC_SIZE_CHK_2D(gr->size);

	RC_TRY( allocate_graphic_data(gr) );

	if(gr->i_info[0] == 24){
		padding = (gr->size.x * 3) % 4;
		if(padding != 0){
			padding = 4 - padding; 
		}
		for(ii1=(gr->size.y-1); ; ii1--){
			for(ii2=0; (unsigned long)ii2<gr->size.x; ii2++){
				/* RGB$B$N=gHV$KCm0U!%(BBMP$B;EMM(B */
				if((cl = fgetc(fp)) == EOF) return(READ_ERROR_RC);
				  gr->data.rgb24[0][ii1][ii2].b = (unsigned char)cl;
				if((cl = fgetc(fp)) == EOF) return(READ_ERROR_RC);
				  gr->data.rgb24[0][ii1][ii2].g = (unsigned char)cl;
				if((cl = fgetc(fp)) == EOF) return(READ_ERROR_RC);
				  gr->data.rgb24[0][ii1][ii2].r = (unsigned char)cl;
			}
			for(ii3=0; ii3<padding; ii3++){
				if( (cl = fgetc(fp)) == EOF ) return(READ_ERROR_RC);
			}
			if(ii1 == 0) break;
		}
	}else if(gr->i_info[0] == 8){
		GRAPHIC_RGB palette[256];
		RC_TRY( read_bmp_palette(fp, gr, palette) );
		padding = gr->size.x % 4;
		if(padding != 0){
			padding = 4 - padding; 
		}
		for(ii1=(gr->size.y)-1; ; ii1--){
			for(ii2=0; (unsigned long)ii2<gr->size.x; ii2++){
				if((cl = fgetc(fp)) == EOF) return(READ_ERROR_RC);
				gr->data.rgb24[0][ii1][ii2].b = palette[cl].b;
				gr->data.rgb24[0][ii1][ii2].g = palette[cl].g;
				gr->data.rgb24[0][ii1][ii2].r = palette[cl].r;
			}
			for(ii3=0; ii3<padding; ii3++){
				if((cl = fgetc(fp)) == EOF) return(READ_ERROR_RC);
			}
			if(ii1 == 0) break;
		}
	}else{
		fprintf(stderr, "Unknown or non-support color bit.\n");
		return(READ_ERROR_RC);
	}

	return(NORMAL_RC);
}


RC write_bmp(FILE *fp, GRAPHIC gr)
{
	RC_TRY( write_bmp_header(fp, gr) );

	if(gr.type == GRAPHIC_RGB24){
		RC_TRY( write_bmp_image_rgb24(fp, gr) );
	}else{
		fprintf(stderr, "It is Non-Implementnon GRAPHIC_TYPE, or "
		                "unknown supported GRAPHIC_TYPE\n");
		return(WRITE_ERROR_RC);
	}
	return(NORMAL_RC);
}


RC write_bmp_header(FILE *fp, GRAPHIC gr)
{
	unsigned int color_num;
	unsigned long header_size;
	unsigned long file_size;
	int padding;
	
	GRAPHIC_SIZE_CHK_3D(gr.size);

	if(gr.type == GRAPHIC_RGB24){
		color_num = 0;  /* $B%Q%l%C%H$NM-L5(B */
		header_size = 14 + 40 + 4*color_num;
		padding = (gr.size.x * 3) % 4;
		if(padding != 0){
			padding = 4 - padding; 
		}
		file_size = header_size+(gr.size.x*3+padding)*(gr.size.y);
	}else{
		fprintf(stderr, "bad GRAPHIC_TYPE <");
		return(WRITE_ERROR_RC);
	}

	/* BMP$B%U%!%$%k%X%C%@(B */
	if(fputs("BM",fp) == EOF)return(WRITE_ERROR_RC);
	RC_TRY( lwrite_4bytes(fp, file_size) );   /* $B%U%!%$%k%5%$%:(B */
	RC_TRY( lwrite_2bytes(fp, 0) );           /* $BM=LsNN0h(B */
	RC_TRY( lwrite_2bytes(fp, 0) );           /* $BM=LsNN0h(B */
	RC_TRY( lwrite_4bytes(fp, header_size) ); /* $B%X%C%@%5%$%:(B */

	/* BMP$B>pJs%X%C%@(B */
	RC_TRY( lwrite_4bytes(fp, 40) );        /* $B>pJs%5%$%:(B */
	RC_TRY( lwrite_4bytes(fp, gr.size.x) ); /* $B2hA|(BX$B%5%$%:(B */
	RC_TRY( lwrite_4bytes(fp, gr.size.y) ); /* $B2hA|(BY$B%5%$%:(B */
	RC_TRY( lwrite_2bytes(fp, 1) );         /* $BLL?t(B */
	if(gr.type == GRAPHIC_RGB24){               /* $B?'%S%C%H?t(B */
		RC_TRY( lwrite_2bytes(fp, 24) );
	}else{
		fprintf(stderr, "bad GRAPHIC_TYPE <");
		return(WRITE_ERROR_RC);
	}
	RC_TRY( lwrite_4bytes(fp, 0) ); /* $B05=LJ}<0(B */
	RC_TRY( lwrite_4bytes(fp, 0) ); /* $B05=L%5%$%:(B */
	/*$B?eJ?2rA|EY(B[dot/m]*/
	RC_TRY( lwrite_4bytes(fp, (unsigned long)(gr.resolution.x*1000.0)) );
	/*$B?bD>2rA|EY(B[dot/m]*/
	RC_TRY( lwrite_4bytes(fp, (unsigned long)(gr.resolution.x*1000.0)) );
	if(gr.type == GRAPHIC_RGB24){       /* $B?'?t(B */
		RC_TRY( lwrite_4bytes(fp, 24) ); 
	}else{
		fprintf(stderr, "bad GRAPHIC_TYPE <");
		return(WRITE_ERROR_RC);
	}
	if(gr.type == GRAPHIC_RGB24){       /* $B=EMW?'?t(B */ 
		RC_TRY( lwrite_4bytes(fp, 24) ); 
	}else{
		fprintf(stderr, "bad GRAPHIC_TYPE <");
		return(WRITE_ERROR_RC);
	}

	return(NORMAL_RC);
}


RC write_bmp_image_rgb24(FILE *fp, GRAPHIC gr)
{
	unsigned int ii1, ii2, ii3;
	unsigned int padding;

	GRAPHIC_SIZE_CHK_3D(gr.size);

	if(gr.type != GRAPHIC_RGB24){
		fprintf(stderr, "This is not RGB 24 bit");
		return(ARG_ERROR_RC);
	}

	RC_NULL_CHK(gr.data.rgb24);

	padding = (gr.size.x * 3) % 4;
	if(padding != 0){
		padding = 4 - padding; 
	}
	for(ii1=(gr.size.y)-1; ; ii1--){
		for(ii2=0; (unsigned long)ii2<(gr.size.x); ii2++){
			GRAPHIC_RGB rgb24 = gr.data.rgb24[0][ii1][ii2];
			if(fputc(rgb24.b,fp) == EOF) return(WRITE_ERROR_RC);
			if(fputc(rgb24.g,fp) == EOF) return(WRITE_ERROR_RC);
			if(fputc(rgb24.r,fp) == EOF) return(WRITE_ERROR_RC);
		}
		for(ii3=0; ii3<padding; ii3++){
			if(fputc('0',fp) == EOF) return(WRITE_ERROR_RC);
		}
		if(ii1==0) break;
	}

	return(NORMAL_RC);
}

