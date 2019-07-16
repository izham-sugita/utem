/*********************************************************************
 * tiff_utl.c
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

/* $Id: tiff_utl.c,v 1.9 2003/11/26 08:10:15 sasaoka Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include "rc.h"
#include "bin_io.h"
#include "graphic_utl.h"

/*
* 1つのファイル内に複数の画像を持つTIFFイメージ
* 及び、圧縮したTIFFイメージはサポートしない
* 解像度の単位に関する情報が不明のため、記述している値をそのまま入力している
* 現在対応しているのは "24bit","8bitパレットカラー","8Bitモノクロ","2bit"
*
* TIFF_UTLでのi_infoの内容
*  i_info[0] = 色ビット数[bit/dot]
*            = 1 4 8 24 32
*  i_info[1] = 画像開始位置
*  i_info[2] = パレット開始位置
*  i_info[3] = 色表現(0-8)
*/

RC static read_tiff_palette(FILE *fp, GRAPHIC *gr,
                            ENDIAN endian, GRAPHIC_RGB *palette);

RC read_tiff(FILE *fp, GRAPHIC *gr)
{
	ENDIAN endian;

	init_graphic(gr);
	RC_TRY( read_tiff_header(fp, gr, &endian) );
	RC_TRY( read_tiff_image(fp, gr, endian) );

	return(NORMAL_RC);
}


RC read_tiff_header(FILE *fp, GRAPHIC *gr, ENDIAN *endian)
{
	RC rc;
	int ii1;
	int itmp;
	unsigned long ltmp;
	unsigned ui, ui1, ui2;
	unsigned long ul;
	unsigned long cpos;
	int num_tag;
	int num_data;

	/* ヘッダ識別領域 */
	if((ui1 = fgetc(fp)) == EOF) return(READ_ERROR_RC);
	if((ui2 = fgetc(fp)) == EOF) return(READ_ERROR_RC);
	if( (ui1 == 'I') && (ui2 == 'I') ){
		*endian = LH;
	}else if( (ui1 == 'M') && (ui2 == 'M') ){
		*endian = HL;
	}else{
		return(UNKNOWN_ERROR_RC);
		fprintf(stderr, "Is it TIFF FILE?\n");
	}

	/* 版数 (42) */
	RC_TRY( read_2bytes(fp, *endian, &ui) );
	if(ui != 42){ fprintf(stderr, "edition number < ");
		return(READ_ERROR_RC);
	}

	/* IFD(Image File Header) 開始位置 */
	RC_TRY( read_4bytes(fp, *endian, &ul) );

	/* IFD Position に移動 */ 
	if(fseek(fp, ul, SEEK_SET) != 0){
		fprintf(stderr, "seek IFD position < ");
		return(SEEK_ERROR_RC);
	}

	/* IFD0 タグ数 */
	RC_TRY( read_2bytes(fp, *endian, &ui) );
	if(ui < 1){
		fprintf(stderr, "tag number < ");
		return(READ_ERROR_RC);
	}
	num_tag = ui;

	/* IFD0 */
	for(ii1=0; ii1<num_tag; ii1++){
		RC_TRY( read_2bytes(fp, *endian, &ui) );
		switch(ui){
			case 256: /* 画像 X size */
				RC_TRY( read_2bytes(fp, *endian, &ui) );
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				if(ul < 1){
					fprintf(stderr, "case 256 < ");
					return(READ_ERROR_RC);
				}else{
					gr->size.x = ul;
					gr->size.z = 1;
				}
				break;

			case 257: /* 画像 Y size */
				RC_TRY( read_2bytes(fp, *endian, &ui) );
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				if(ul < 1){
					fprintf(stderr, "case 257 < ");
					return(READ_ERROR_RC);
				}else{
					gr->size.y = ul;
				} break;

			case 258: /* 色ビット数[bit/dot] */
				RC_TRY( read_2bytes(fp, *endian, &ui) );
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				num_data = ul;
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				if(num_data == 1){
					if((ul==8)||(ul==1)){
						gr->i_info[0] = ul;
					}else{
						fprintf(stderr, "color bit = %ld\n", ul);
						fprintf(stderr, "unknown or non support color bit\n");
						fprintf(stderr, "case 258 < ");
						return(READ_ERROR_RC);
					}
				}else if(num_data > 1){
					cpos = ftell(fp);                 /* 現在の位置を記憶 */
					if(fseek(fp, ul, SEEK_SET) != 0){ /*色ビット数記述部に移動*/
						fprintf(stderr, "case 258 < ");
						return(SEEK_ERROR_RC);
					}
					itmp = 0;
					for(ii1=0;ii1<num_data;ii1++){
						RC_TRY( read_2bytes(fp, *endian, &ui) );
						itmp+=ui;
					}
					if(itmp == 24){
						gr->i_info[0] = itmp;
					}else{
						fprintf(stderr, "color bit = %d\n", itmp);
						fprintf(stderr, "unknown or non support color bit\n");
						fprintf(stderr, "case 258 < ");
						return(READ_ERROR_RC);
					}
					if(fseek(fp, cpos, SEEK_SET) != 0){ /* もとの位置に移動 */
						fprintf(stderr, "case 258 < ");
						return(SEEK_ERROR_RC);
					}
				}else{
					fprintf(stderr, "case 258 < ");
					return(READ_ERROR_RC);
				}
				break;

			case 259: /* 圧縮 */
				/*
				* 1: 無圧縮
				* 2: CCITT 1D
				* 3: CCITT Group3
				* 4: CCITT Group4
				* 5: LZW
				* 6: JPEG
				* 32771: 無圧縮
				* 32773: RLE, PackBits
				*/
				RC_TRY( read_2bytes(fp, *endian, &ui) );
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				if(ul==1||ul==32771){
					break;
				}else{
					fprintf(stderr, "Compression Method = %ld\n", ul);
					fprintf(stderr, "Don't support compression image !!\n");
					fprintf(stderr, "case 259 < ");
					return(READ_ERROR_RC);
				}
				break;

			case 262: /* 色表現 */
				/*
				* 0: 白が0のモノクロ(TIFF Version 4,5,6)
				* 1: 黒が0のモノクロ(Ver4,5,6)
				* 2: RGB(Ver4,5,6)
				* 3: RGBパレット(Ver5,6)
				* 4: 透過(Ver6)
				* 5: CMYK(Ver6)
				* 6: YCbCr(Ver6)
				* 8: CIELab(Ver6)
				*/
				RC_TRY( read_2bytes(fp, *endian, &ui) );
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				if(ul==0||ul==1||ul==2||ul==3){
					gr->i_info[3] = ul;
					break;
				}else{
					fprintf(stderr, "non support photometric interpretation\n");
					fprintf(stderr, "case 262 < ");
					return(READ_ERROR_RC);
				}
				break;

			case 273: /* 画像開始位置[byte] */
				RC_TRY( read_2bytes(fp, *endian, &ui) );
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				num_data = ul;
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				if(num_data == 1){
					gr->i_info[1] = ul;
				}else if(num_data > 1){
					cpos = ftell(fp); /* 現在の位置を記憶 */
					 /* 画像開始位置の記述部に移動 */
					if(fseek(fp, ul, SEEK_SET) != 0){
						fprintf(stderr, "case 273 < ");
						return(SEEK_ERROR_RC);
					}
					/* GIMPのTIFFには何故か4つの画像開始位置の記述がある */
					/* その最初のデータが画像開始位置 */
					RC_TRY( read_4bytes(fp, *endian, &ul) );
					if(ul < 1){
						fprintf(stderr, "unknown image data potion\n");
						fprintf(stderr, "case 273 < ");
						return(READ_ERROR_RC);
					}else{
						gr->i_info[1] = ul;
					}
					if(fseek(fp, cpos, SEEK_SET) != 0){ /* もとの位置に移動 */
						fprintf(stderr, "case 273 < ");
						return(SEEK_ERROR_RC);
					}
				}else{
					fprintf(stderr, "case 273 < ");
					return(READ_ERROR_RC);
				}
				break;
				
			case 282: /* 画像 X 解像度 */
				RC_TRY( read_2bytes(fp, *endian, &ui) );
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				num_data = ul;
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				if(num_data == 1){
					gr->resolution.x = ul;
				}else if(num_data > 1){
					cpos = ftell(fp); /* 現在の位置を記憶 */
					if(fseek(fp, ul, SEEK_SET) != 0){ /* 解像度記述部に移動 */
						fprintf(stderr, "case 282 < ");
						return(SEEK_ERROR_RC);
					}
					ltmp = 0;
					for(ii1=0;ii1<num_data;ii1++){
						RC_TRY( read_4bytes(fp, *endian, &ul) );
						ltmp+=ul;
					}
					gr->resolution.x = ltmp;
					if(fseek(fp, cpos, SEEK_SET) != 0){ /* もとの位置に移動 */
						fprintf(stderr, "case 282 < ");
						return(SEEK_ERROR_RC);
					}
				}else{
					fprintf(stderr, "case 282 < ");
					return(READ_ERROR_RC);
				}
				break;

			case 283: /* 画像 Y 解像度 */
				RC_TRY( read_2bytes(fp, *endian, &ui) );
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				num_data = ul;
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				if(num_data == 1){
					gr->resolution.y = ul;
				}else if(num_data > 1){
					cpos = ftell(fp); /* 現在の位置を記憶 */
					if(fseek(fp, ul, SEEK_SET) != 0){ /* 解像度記述部に移動 */
						fprintf(stderr, "case 283 < ");
						return(SEEK_ERROR_RC);
					}
					ltmp = 0;
					for(ii1=0;ii1<num_data;ii1++){
						RC_TRY( read_4bytes(fp, *endian, &ul) );
						ltmp+=ul;
					}
					gr->resolution.y = ltmp;
					if(fseek(fp, cpos, SEEK_SET) != 0){ /* もとの位置に移動 */
						fprintf(stderr, "case 283 < ");
						return(SEEK_ERROR_RC);
					}
				}else{
					fprintf(stderr, "case 283 < ");
					return(READ_ERROR_RC);
				}
				break;

			case 296: /* 解像度単位 (不明な点が多いので項目だけ) */
				RC_TRY( read_2bytes(fp, *endian, &ui) );
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				if(ul < 1){
					fprintf(stderr, "case 296 < ");
					return(READ_ERROR_RC);
				}
				break;

			/* パレット開始位置を示すタグ(320)のデータ数はあり得ないほど多い */
			/* よって、データ数を1として入力 */
			case 320: /* パレット開始位置[byte] */
				RC_TRY( read_2bytes(fp, *endian, &ui) );
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				gr->i_info[2] = ul;
				break;

			default:
				RC_TRY( read_2bytes(fp, *endian, &ui) );
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				break;
		}
	}

	rc = read_4bytes(fp, *endian, &ul);
	if(rc == READ_ERROR_RC ){
		/*
		 * IFD0 がファイルの終端にある時もある。よって、
		 * 次を読み込むとEOFである可能性が高いのでそのまま
		 * NORMAL_RC を返す。
		 */
		return(NORMAL_RC);
	}
	if((rc == NORMAL_RC) && (ul != 0)){
		/* IFD1以降が存在するため処理を中断 */
		fprintf(stderr, "next IFD = %ld\n", ul);
		fprintf(stderr, "Some IFD found.\n Don't support some image !!\n");
		return(READ_ERROR_RC);
	}

	return(NORMAL_RC);
}

RC read_tiff_image(FILE *fp, GRAPHIC *gr, ENDIAN endian)
{
	int ii1, ii2, ii3;
	int cl, cl1;
	GRAPHIC_RGB palette[256];
	

	if(gr->i_info[1] < 0){
		fprintf(stderr, "Bad Image File Position <");
		return(READ_ERROR_RC);
	}

	if(fseek(fp, gr->i_info[1], SEEK_SET) != 0) return(SEEK_ERROR_RC);

	GRAPHIC_SIZE_CHK_3D(gr->size);
	gr->type = GRAPHIC_RGB24;
	RC_TRY( allocate_graphic_data(gr) );

	if((gr->i_info[0] == 24)&&(gr->i_info[3] == 2)){
		GRAPHIC_SIZE_LOOP_2D(gr->size.y, gr->size.x){
			if( (cl=fgetc(fp)) == EOF ) return(READ_ERROR_RC);
			gr->data.rgb24[0][ii1][ii2].r = (unsigned char)cl;
			if( (cl=fgetc(fp)) == EOF ) return(READ_ERROR_RC);
			gr->data.rgb24[0][ii1][ii2].g = (unsigned char)cl;
			if( (cl=fgetc(fp)) == EOF ) return(READ_ERROR_RC);
			gr->data.rgb24[0][ii1][ii2].b = (unsigned char)cl;
		}GRAPHIC_SIZE_LOOP_2D_END;

	}else if((gr->i_info[0] == 8)&&(gr->i_info[3] == 0)){
		GRAPHIC_SIZE_LOOP_2D(gr->size.y, gr->size.x){
			if( (cl=fgetc(fp)) == EOF ) return(READ_ERROR_RC);
			cl = 255 - cl;
			gr->data.rgb24[0][ii1][ii2].r = cl;
			gr->data.rgb24[0][ii1][ii2].g = cl;
			gr->data.rgb24[0][ii1][ii2].b = cl;
		}GRAPHIC_SIZE_LOOP_2D_END;

	}else if((gr->i_info[0] == 8)&&(gr->i_info[3] == 1)){
		GRAPHIC_SIZE_LOOP_2D(gr->size.y, gr->size.x){
			if( (cl=fgetc(fp)) == EOF ) return(READ_ERROR_RC);
			gr->data.rgb24[0][ii1][ii2].r = (unsigned char)cl;
			gr->data.rgb24[0][ii1][ii2].g = (unsigned char)cl;
			gr->data.rgb24[0][ii1][ii2].b = (unsigned char)cl;
		}GRAPHIC_SIZE_LOOP_2D_END;

	}else if((gr->i_info[0] == 8)&&(gr->i_info[3] == 3)){
		RC_TRY( read_tiff_palette(fp, gr, endian, palette) );
		GRAPHIC_SIZE_LOOP_2D(gr->size.y, gr->size.x){
			if( (cl=fgetc(fp)) == EOF ) return(READ_ERROR_RC);
			gr->data.rgb24[0][ii1][ii2].r = palette[cl].r;
			gr->data.rgb24[0][ii1][ii2].g = palette[cl].g;
			gr->data.rgb24[0][ii1][ii2].b = palette[cl].b;
		}GRAPHIC_SIZE_LOOP_2D_END;

	}else if((gr->i_info[0] == 1)&&(gr->i_info[3] == 0)){
		for(ii1=0; (unsigned long)ii1<gr->size.y; ii1++){
			for(ii2=0; (unsigned long)ii2<gr->size.x; ii2+=8){
				if( (cl=fgetc(fp)) == EOF ) return(READ_ERROR_RC);
				for(ii3=0; ii3<8; ii3++){
					if((unsigned long)(ii2+ii3) > (gr->size.x - 1)) continue;
					cl1 = ((cl>>(7-ii3))&0x1);
					if(cl1 == 0) cl1 = 255;
					else         cl1 = 0;
					gr->data.rgb24[0][ii1][ii2+ii3].r = cl1;
					gr->data.rgb24[0][ii1][ii2+ii3].g = cl1;
					gr->data.rgb24[0][ii1][ii2+ii3].b = cl1;
				}
			}
		}

	}else if((gr->i_info[0] == 1)&&(gr->i_info[3] == 1)){
		for(ii1=0; (unsigned long)ii1<gr->size.y; ii1++){
			for(ii2=0; (unsigned long)ii2<gr->size.x; ii2+=8){
				if( (cl=fgetc(fp)) == EOF ) return(READ_ERROR_RC);
				for(ii3=0; ii3<8; ii3++){
					if((unsigned long)(ii2+ii3) > (gr->size.x - 1)) continue;
					cl1 = ((cl>>(7-ii3))&0x1);
					if(cl1 == 0) cl1 = 0;
					else         cl1 = 255;
					gr->data.rgb24[0][ii1][ii2+ii3].r = cl1;
					gr->data.rgb24[0][ii1][ii2+ii3].g = cl1;
					gr->data.rgb24[0][ii1][ii2+ii3].b = cl1;
				}
			}
		}

	}else{
		fprintf(stderr, "color bit = %d\n", gr->i_info[0]);
		fprintf(stderr, "photometric interpretation = %d\n", gr->i_info[3]);
		fprintf(stderr, "unknown or non-support color bit\n");
		return(READ_ERROR_RC);
	}

	return(NORMAL_RC);
}
	

RC static read_tiff_palette(FILE *fp, GRAPHIC *gr,
                            ENDIAN endian, GRAPHIC_RGB *palette)
{
	int ii1;
	unsigned int ui;


	if(gr->i_info[2] < 0){
		fprintf(stderr, "bad palette file position <");
		return(READ_ERROR_RC);
	}
	if(fseek(fp, gr->i_info[2], SEEK_SET) != 0) return(SEEK_ERROR_RC);

	/*
	* TIFFのパレットは16bitでRRR...GGG...BBB...の順番でで記述されている
	* RGB256色に納めるため、下位8bitは切り捨て
	*/
	for(ii1=0; ii1<(1<<(gr->i_info[0])); ii1++){
		RC_TRY( read_2bytes(fp, endian, &ui) );
		ui>>=8;
		palette[ii1].r = (unsigned char)ui;
	}
	for(ii1=0; ii1<(1<<(gr->i_info[0])); ii1++){
		RC_TRY( read_2bytes(fp, endian, &ui) );
		ui>>=8;
		palette[ii1].g = (unsigned char)ui;
	}
	for(ii1=0; ii1<(1<<(gr->i_info[0])); ii1++){
		RC_TRY( read_2bytes(fp, endian, &ui) );
		ui>>=8;
		palette[ii1].b = (unsigned char)ui;
	}

	return(NORMAL_RC);
}

