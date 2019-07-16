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
* 1$B$D$N%U%!%$%kFb$KJ#?t$N2hA|$r;}$D(BTIFF$B%$%a!<%8(B
* $B5Z$S!"05=L$7$?(BTIFF$B%$%a!<%8$O%5%]!<%H$7$J$$(B
* $B2rA|EY$NC10L$K4X$9$k>pJs$,ITL@$N$?$a!"5-=R$7$F$$$kCM$r$=$N$^$^F~NO$7$F$$$k(B
* $B8=:_BP1~$7$F$$$k$N$O(B "24bit","8bit$B%Q%l%C%H%+%i!<(B","8Bit$B%b%N%/%m(B","2bit"
*
* TIFF_UTL$B$G$N(Bi_info$B$NFbMF(B
*  i_info[0] = $B?'%S%C%H?t(B[bit/dot]
*            = 1 4 8 24 32
*  i_info[1] = $B2hA|3+;O0LCV(B
*  i_info[2] = $B%Q%l%C%H3+;O0LCV(B
*  i_info[3] = $B?'I=8=(B(0-8)
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

	/* $B%X%C%@<1JLNN0h(B */
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

	/* $BHG?t(B (42) */
	RC_TRY( read_2bytes(fp, *endian, &ui) );
	if(ui != 42){ fprintf(stderr, "edition number < ");
		return(READ_ERROR_RC);
	}

	/* IFD(Image File Header) $B3+;O0LCV(B */
	RC_TRY( read_4bytes(fp, *endian, &ul) );

	/* IFD Position $B$K0\F0(B */ 
	if(fseek(fp, ul, SEEK_SET) != 0){
		fprintf(stderr, "seek IFD position < ");
		return(SEEK_ERROR_RC);
	}

	/* IFD0 $B%?%0?t(B */
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
			case 256: /* $B2hA|(B X size */
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

			case 257: /* $B2hA|(B Y size */
				RC_TRY( read_2bytes(fp, *endian, &ui) );
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				if(ul < 1){
					fprintf(stderr, "case 257 < ");
					return(READ_ERROR_RC);
				}else{
					gr->size.y = ul;
				} break;

			case 258: /* $B?'%S%C%H?t(B[bit/dot] */
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
					cpos = ftell(fp);                 /* $B8=:_$N0LCV$r5-21(B */
					if(fseek(fp, ul, SEEK_SET) != 0){ /*$B?'%S%C%H?t5-=RIt$K0\F0(B*/
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
					if(fseek(fp, cpos, SEEK_SET) != 0){ /* $B$b$H$N0LCV$K0\F0(B */
						fprintf(stderr, "case 258 < ");
						return(SEEK_ERROR_RC);
					}
				}else{
					fprintf(stderr, "case 258 < ");
					return(READ_ERROR_RC);
				}
				break;

			case 259: /* $B05=L(B */
				/*
				* 1: $BL505=L(B
				* 2: CCITT 1D
				* 3: CCITT Group3
				* 4: CCITT Group4
				* 5: LZW
				* 6: JPEG
				* 32771: $BL505=L(B
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

			case 262: /* $B?'I=8=(B */
				/*
				* 0: $BGr$,(B0$B$N%b%N%/%m(B(TIFF Version 4,5,6)
				* 1: $B9u$,(B0$B$N%b%N%/%m(B(Ver4,5,6)
				* 2: RGB(Ver4,5,6)
				* 3: RGB$B%Q%l%C%H(B(Ver5,6)
				* 4: $BF)2a(B(Ver6)
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

			case 273: /* $B2hA|3+;O0LCV(B[byte] */
				RC_TRY( read_2bytes(fp, *endian, &ui) );
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				num_data = ul;
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				if(num_data == 1){
					gr->i_info[1] = ul;
				}else if(num_data > 1){
					cpos = ftell(fp); /* $B8=:_$N0LCV$r5-21(B */
					 /* $B2hA|3+;O0LCV$N5-=RIt$K0\F0(B */
					if(fseek(fp, ul, SEEK_SET) != 0){
						fprintf(stderr, "case 273 < ");
						return(SEEK_ERROR_RC);
					}
					/* GIMP$B$N(BTIFF$B$K$O2?8N$+(B4$B$D$N2hA|3+;O0LCV$N5-=R$,$"$k(B */
					/* $B$=$N:G=i$N%G!<%?$,2hA|3+;O0LCV(B */
					RC_TRY( read_4bytes(fp, *endian, &ul) );
					if(ul < 1){
						fprintf(stderr, "unknown image data potion\n");
						fprintf(stderr, "case 273 < ");
						return(READ_ERROR_RC);
					}else{
						gr->i_info[1] = ul;
					}
					if(fseek(fp, cpos, SEEK_SET) != 0){ /* $B$b$H$N0LCV$K0\F0(B */
						fprintf(stderr, "case 273 < ");
						return(SEEK_ERROR_RC);
					}
				}else{
					fprintf(stderr, "case 273 < ");
					return(READ_ERROR_RC);
				}
				break;
				
			case 282: /* $B2hA|(B X $B2rA|EY(B */
				RC_TRY( read_2bytes(fp, *endian, &ui) );
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				num_data = ul;
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				if(num_data == 1){
					gr->resolution.x = ul;
				}else if(num_data > 1){
					cpos = ftell(fp); /* $B8=:_$N0LCV$r5-21(B */
					if(fseek(fp, ul, SEEK_SET) != 0){ /* $B2rA|EY5-=RIt$K0\F0(B */
						fprintf(stderr, "case 282 < ");
						return(SEEK_ERROR_RC);
					}
					ltmp = 0;
					for(ii1=0;ii1<num_data;ii1++){
						RC_TRY( read_4bytes(fp, *endian, &ul) );
						ltmp+=ul;
					}
					gr->resolution.x = ltmp;
					if(fseek(fp, cpos, SEEK_SET) != 0){ /* $B$b$H$N0LCV$K0\F0(B */
						fprintf(stderr, "case 282 < ");
						return(SEEK_ERROR_RC);
					}
				}else{
					fprintf(stderr, "case 282 < ");
					return(READ_ERROR_RC);
				}
				break;

			case 283: /* $B2hA|(B Y $B2rA|EY(B */
				RC_TRY( read_2bytes(fp, *endian, &ui) );
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				num_data = ul;
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				if(num_data == 1){
					gr->resolution.y = ul;
				}else if(num_data > 1){
					cpos = ftell(fp); /* $B8=:_$N0LCV$r5-21(B */
					if(fseek(fp, ul, SEEK_SET) != 0){ /* $B2rA|EY5-=RIt$K0\F0(B */
						fprintf(stderr, "case 283 < ");
						return(SEEK_ERROR_RC);
					}
					ltmp = 0;
					for(ii1=0;ii1<num_data;ii1++){
						RC_TRY( read_4bytes(fp, *endian, &ul) );
						ltmp+=ul;
					}
					gr->resolution.y = ltmp;
					if(fseek(fp, cpos, SEEK_SET) != 0){ /* $B$b$H$N0LCV$K0\F0(B */
						fprintf(stderr, "case 283 < ");
						return(SEEK_ERROR_RC);
					}
				}else{
					fprintf(stderr, "case 283 < ");
					return(READ_ERROR_RC);
				}
				break;

			case 296: /* $B2rA|EYC10L(B ($BITL@$JE@$,B?$$$N$G9`L\$@$1(B) */
				RC_TRY( read_2bytes(fp, *endian, &ui) );
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				RC_TRY( read_4bytes(fp, *endian, &ul) );
				if(ul < 1){
					fprintf(stderr, "case 296 < ");
					return(READ_ERROR_RC);
				}
				break;

			/* $B%Q%l%C%H3+;O0LCV$r<($9%?%0(B(320)$B$N%G!<%??t$O$"$jF@$J$$$[$IB?$$(B */
			/* $B$h$C$F!"%G!<%??t$r(B1$B$H$7$FF~NO(B */
			case 320: /* $B%Q%l%C%H3+;O0LCV(B[byte] */
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
		 * IFD0 $B$,%U%!%$%k$N=*C<$K$"$k;~$b$"$k!#$h$C$F!"(B
		 * $B<!$rFI$_9~$`$H(BEOF$B$G$"$k2DG=@-$,9b$$$N$G$=$N$^$^(B
		 * NORMAL_RC $B$rJV$9!#(B
		 */
		return(NORMAL_RC);
	}
	if((rc == NORMAL_RC) && (ul != 0)){
		/* IFD1$B0J9_$,B8:_$9$k$?$a=hM}$rCfCG(B */
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
	* TIFF$B$N%Q%l%C%H$O(B16bit$B$G(BRRR...GGG...BBB...$B$N=gHV$G$G5-=R$5$l$F$$$k(B
	* RGB256$B?'$KG<$a$k$?$a!"2<0L(B8bit$B$O@Z$j<N$F(B
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

