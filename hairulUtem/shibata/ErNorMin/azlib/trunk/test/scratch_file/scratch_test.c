#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "base.h"
#include "mathematics.h"

int main (int argc, char **argv)
{
	SCRATCH_FILE *sf;
	FILE *fp;
	int ii1, ii2;
	int sc_num[9];

	if(argc < 2){
		fprintf(stderr, "Usage : %s [scratch directory(ex:/tmp)]\n", argv[0]);
		fprintf(stderr, "scratch directory need 10GB free space\n");
		return(EXIT_FAILURE);
	}

	/* Memory Manager */
	RC_TRY_MAIN( mm_init((size_t)200*MEGA_BYTE) );

	/* log Print */
	RC_TRY_MAIN( set_log_level(5) );
	RC_TRY_MAIN( set_log_file(0, stderr) );

	sf = open_scratch_file(argv[1]);
	if(sf == NULL){
		RC_TRY_MAIN( log_printf(0, "Error(Write): Scratch File Open\n") );
		return(EXIT_FAILURE);
	}

	RC_TRY_MAIN( log_printf(1,
	    "\nTry to write 9th scratch area in a scratch-file.\n"
	    "Each area is filled with the same character.\n"
	    "1st area : [aaa.....aaa] (Position: 0Byte -- 1G-1 Byte)\n"
	    "2nd area : [bbb.....bbb] (Position: 1GB   -- 2G-1 Byte)\n"
	    "...\n"
	    "9th area : [iii.....iii] (Position: 8GB   -- 9G-1 Byte)\n\n"
	) );

	/* 9GBの スクラッチファイル書き込みテスト開始 */
	RC_TRY_MAIN( log_printf(1, "[Scratch File Writting Test]\n") );

	/* 文字 a,b,c,d,e,f,g,h,i を1GB分づつ書き込む */
	for(ii1=0; ii1<9; ii1++){
		/* スクラッチファイル書き込み開始 */
		fp = add_scratch_file_start(sf);
		if(fp == NULL){
			RC_TRY_MAIN( log_printf(0, "Error: Add Scratch File\n") );
			return(EXIT_FAILURE);
		}

		/* 各文字で1GB分を埋め尽くす．*/
		for(ii2=0; ii2<GIGA_BYTE; ii2++){
			if(EOF == fputc('a'+ii1, fp))
				RC_TRY_MAIN( log_printf(1, "Write Error in %d GB\n", ii1+1) );
			if(ii2%MEGA_BYTE == 0){
				double percent = (double)ii2/GIGA_BYTE;
				RC_TRY_MAIN( progress_bar_printf(5, percent,
				             "  Write %2.2f GB. Now filling char \"%c\"",
				             percent+ii1, 'a'+ii1) );
			}
		}
		RC_TRY_MAIN( log_printf(5, "\n") );

		/* スクラッチファイル書き込み終了 */
		sc_num[ii1] = add_scratch_file_end(sf);
	}
	RC_TRY_MAIN( log_printf(5, "Complete !!\n\n") );

	/* スクラッチファイルの SEEK テスト */
	RC_TRY_MAIN( log_printf(1, "[Scratch File SEEK Test]\n") );
	for(ii1=0; ii1<9; ii1++){
		RC_TRY_MAIN( log_printf(1, "Write Char Check : \"%c\"\n", 'a'+ii1) );
		fp = read_scratch_file(sf, sc_num[ii1]);
		if(fp == NULL){
			RC_TRY_MAIN( log_printf(0, "Error(Seek): Add Scratch File\n") );
			return(EXIT_FAILURE);
		}

		/* WRAP64_INT の出力は環境依存なので要修正 */
		RC_TRY_MAIN( log_printf(1, " First : Potistion = %lld  %s\n",
		                        wrap64_ftell(fp),
		                        (fgetc(fp) == 'a'+ii1 ? "OK" : "NO")) );
		RC_TRY_MAIN( wrap64_fseek(fp, (WRAP64_INT)(GIGA_BYTE-2), SEEK_CUR) );
		RC_TRY_MAIN( log_printf(1, " Last  : Potistion = %lld  %s\n",
		                        wrap64_ftell(fp),
		                        (fgetc(fp) == 'a'+ii1 ? "OK" : "NO")) );
	}
	RC_TRY_MAIN( log_printf(1, "\n") );

	/* スクラッチファイルの消去 */
	RC_TRY_MAIN( close_scratch_file(sf) );

	/* Memory Manager */
	RC_TRY( mm_terminate() );

	return(EXIT_SUCCESS);
}

