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

/* $Id: bmp_utl.c 432 2005-08-24 03:00:31Z sasaoka $ */

#include <stdio.h>
#include <stdlib.h>
#include "rc.h"
#include "log_printf.h"
#include "bin_io.h"
#include "graphic_utl.h"

/*
 * 圧縮したBMPイメージはサポートしない
 * Write は24bitのみ対応
 * Read は24bit,8bit[palette]のみ対応
 * (4bit,2bitの画像はパッディングの処理が不明のため未対応)
 */

static RC read_bmp_palette(FILE *fp, GRAPHIC *gr, GRAPHIC_RGB *palette);
static RC read_bmp_image(FILE *fp, GRAPHIC *gr);

RC
read_bmp (FILE *fp, GRAPHIC *gr)
{
	init_graphic(gr);
	RC_TRY( read_bmp_header(fp, gr) );
	RC_TRY( read_bmp_image(fp, gr) );
	return(NORMAL_RC);
}

RC
read_bmp_header (FILE *fp, GRAPHIC *gr)
{
	unsigned int  ui, ui1, ui2;
	unsigned long ul;

	/* ヘッダ識別領域 */
	if((ui1 = fgetc(fp)) == EOF) return(READ_ERROR_RC);
	if((ui2 = fgetc(fp)) == EOF) return(READ_ERROR_RC);
	if( (ui1 == 'B') && (ui2 == 'M') ){
		gr->source = GRAPHIC_BMP;
	}else{
		RC_TRY( error_printf(1001, "BMP") );
		return(READ_ERROR_RC);
	}

    /* ファイルサイズ */
	RC_TRY( lread_4bytes(fp, &ul) );

	/* 予約領域1 & 予約領域2 */
	RC_TRY( lread_2bytes(fp, &ui) );
	RC_TRY( lread_2bytes(fp, &ui) );

    /* ヘッダサイズ */
	RC_TRY( lread_4bytes(fp, &ul) );

    /* 情報サイズ */
	RC_TRY( lread_4bytes(fp, &ul) );

    /* 画像Xサイズ */
	RC_TRY( lread_4bytes(fp, &ul) );
	if(ul > 0){
		gr->size.x = ul;
	}else{
		return(READ_ERROR_RC);
	}

    /* 画像Yサイズ */
	RC_TRY( lread_4bytes(fp, &ul) );
	if(ul > 0){
		gr->size.y = ul;
	}else{
		return(READ_ERROR_RC);
	}

    /* 画像Zサイズ */
	gr->size.z = 1;

	/* 面数 */
	RC_TRY( lread_2bytes(fp, &ui) );

	/* 色ビット数 */
	RC_TRY( lread_2bytes(fp, &ui) );
	if( (ui == 8) || (ui == 24) ){
		gr->type = GRAPHIC_RGB24;
		gr->i_info[0] = ui;
	}else if( (ui == 1) || (ui == 4) ){
		return(IMPLEMENT_ERROR_RC);
	}else{
		return(READ_ERROR_RC);
	}

	/* 圧縮方式 & 圧縮サイズ */
	RC_TRY( lread_4bytes(fp, &ul) );
	if(ul != 0){
		RC_TRY( error_printf(1002, "Compressed BMP") );
		return(READ_ERROR_RC);
	}
	RC_TRY( lread_4bytes(fp, &ul) );

	/* 水平解像度 */
	RC_TRY( lread_4bytes(fp, &ul) );
	gr->resolution.x = ul/1000.0;

	/* 垂直解像度 */
	RC_TRY( lread_4bytes(fp, &ul) );
	gr->resolution.y = ul/1000.0;

	/* 色数 */
	RC_TRY( lread_4bytes(fp, &ul) );

	/* 重要色数 */
	RC_TRY( lread_4bytes(fp, &ul) );

	return(NORMAL_RC);
}


static RC
read_bmp_palette (FILE *fp, GRAPHIC *gr, GRAPHIC_RGB *palette)
{
	int ii1, cl;
	for(ii1=0; ii1<(1<<(gr->i_info[0])); ii1++){
		if((cl = fgetc(fp)) == EOF) return(READ_ERROR_RC);
		  palette[ii1].b = (unsigned char)cl;
		if((cl = fgetc(fp)) == EOF) return(READ_ERROR_RC);
		  palette[ii1].g = (unsigned char)cl;
		if((cl = fgetc(fp)) == EOF) return(READ_ERROR_RC);
		  palette[ii1].r = (unsigned char)cl;
		/* padding代りの4個目がある */
		if((cl = fgetc(fp)) == EOF) return(READ_ERROR_RC);
	}
	return(NORMAL_RC);
}


static RC
read_bmp_image (FILE *fp, GRAPHIC *gr)
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
				/* RGBの順番に注意．BMP仕様 */
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
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


RC
write_bmp (FILE *fp, GRAPHIC gr)
{
	RC_TRY( write_bmp_header(fp, gr) );

	if(gr.type == GRAPHIC_RGB24){
		RC_TRY( write_bmp_image_rgb24(fp, gr) );
	}else{
		RC_TRY( error_printf(7001) );
		return(IMPLEMENT_ERROR_RC);
	}
	return(NORMAL_RC);
}


RC
write_bmp_header (FILE *fp, GRAPHIC gr)
{
	unsigned int color_num;
	unsigned long header_size;
	unsigned long file_size;
	int padding;
	
	GRAPHIC_SIZE_CHK_3D(gr.size);

	if(gr.type == GRAPHIC_RGB24){
		color_num = 0;  /* パレットの有無 */
		header_size = 14 + 40 + 4*color_num;
		padding = (gr.size.x * 3) % 4;
		if(padding != 0){
			padding = 4 - padding; 
		}
		file_size = header_size+(gr.size.x*3+padding)*(gr.size.y);
	}else{
		RC_TRY( error_printf(7001) );
		return(IMPLEMENT_ERROR_RC);
	}

	/* BMPファイルヘッダ */
	if(fputs("BM",fp) == EOF) return(WRITE_ERROR_RC);
	RC_TRY( lwrite_4bytes(fp, file_size) );   /* ファイルサイズ */
	RC_TRY( lwrite_2bytes(fp, 0) );           /* 予約領域 */
	RC_TRY( lwrite_2bytes(fp, 0) );           /* 予約領域 */
	RC_TRY( lwrite_4bytes(fp, header_size) ); /* ヘッダサイズ */

	/* BMP情報ヘッダ */
	RC_TRY( lwrite_4bytes(fp, 40) );        /* 情報サイズ */
	RC_TRY( lwrite_4bytes(fp, gr.size.x) ); /* 画像Xサイズ */
	RC_TRY( lwrite_4bytes(fp, gr.size.y) ); /* 画像Yサイズ */
	RC_TRY( lwrite_2bytes(fp, 1) );         /* 面数 */
	if(gr.type == GRAPHIC_RGB24){               /* 色ビット数 */
		RC_TRY( lwrite_2bytes(fp, 24) );
	}else{
		RC_TRY( error_printf(7001) );
		return(IMPLEMENT_ERROR_RC);
	}
	RC_TRY( lwrite_4bytes(fp, 0) ); /* 圧縮方式 */
	RC_TRY( lwrite_4bytes(fp, 0) ); /* 圧縮サイズ */
	/*水平解像度[dot/m]*/
	RC_TRY( lwrite_4bytes(fp, (unsigned long)(gr.resolution.x*1000.0)) );
	/*垂直解像度[dot/m]*/
	RC_TRY( lwrite_4bytes(fp, (unsigned long)(gr.resolution.x*1000.0)) );
	if(gr.type == GRAPHIC_RGB24){       /* 色数 */
		RC_TRY( lwrite_4bytes(fp, 24) ); 
	}else{
		RC_TRY( error_printf(7001) );
		return(IMPLEMENT_ERROR_RC);
	}
	if(gr.type == GRAPHIC_RGB24){       /* 重要色数 */ 
		RC_TRY( lwrite_4bytes(fp, 24) ); 
	}else{
		RC_TRY( error_printf(7001) );
		return(IMPLEMENT_ERROR_RC);
	}

	return(NORMAL_RC);
}


RC
write_bmp_image_rgb24 (FILE *fp, GRAPHIC gr)
{
	unsigned int ii1, ii2, ii3;
	unsigned int padding;

	GRAPHIC_SIZE_CHK_3D(gr.size);

	if(gr.type != GRAPHIC_RGB24) return(ARG_ERROR_RC);

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

