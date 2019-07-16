#include <stdio.h>

/*
 *  [文字装飾]    [文字の色] [背景の色]
 *  0 : 初期状態   30 : 黒    40 : 黒
 *  1 : 太文字     31 : 赤    41 : 赤
 *  4 : 下線       32 : 緑    42 : 緑
 *  5 : 点滅       33 : 黄    43 : 黄
 *  7 : 色を反転   34 : 青    44 : 青
 *                 35 : 紫    45 : 紫
 *                 36 : 水    46 : 水
 *                 37 : 白    47 : 白
 *  (点滅は端末によっては太文字)
 */

int main(void);

int main(void)
{
	fprintf(stderr, "\x1b[30mBLACK\x1b[0m\n");
	fprintf(stderr, "\x1b[31mRED\x1b[0m\n");
	fprintf(stderr, "\x1b[32mGREEN\x1b[0m\n");
	fprintf(stderr, "\x1b[33mYELLOW\x1b[0m\n");
	fprintf(stderr, "\x1b[34mBLUE\x1b[0m\n");
	fprintf(stderr, "\x1b[35mMAGENTA\x1b[0m\n");
	fprintf(stderr, "\x1b[36mCYAN\x1b[0m\n");
	fprintf(stderr, "\x1b[37mWHITE\x1b[0m\n");
	
	fprintf(stderr, "\x1b[0;31mRED_NORMAL\x1b[0m\n");
	fprintf(stderr, "\x1b[1;31mRED_BOLD\x1b[0m\n");
	fprintf(stderr, "\x1b[4;31mRED_UNDER\x1b[0m\n");
	fprintf(stderr, "\x1b[5;31mRED_FLASH\x1b[0m\n");
	fprintf(stderr, "\x1b[7;31mRED_REV\x1b[0m\n");
	
	fprintf(stderr, "\x1b[0;31;40mRED_NORMAL_BLACK\x1b[0m\n");
	fprintf(stderr, "\x1b[1;31;41mRED_BOLD_RED\x1b[0m\n");
	fprintf(stderr, "\x1b[4;31;42mRED_UNDER_GREEN\x1b[0m\n");
	fprintf(stderr, "\x1b[5;31;43mRED_FLASH_YELLOW\x1b[0m\n");
	fprintf(stderr, "\x1b[7;31;44mRED_REV_BLUE\x1b[0m\n");
	fprintf(stderr, "\x1b[31;45mRED_MAGENTA\x1b[0m\n");
	fprintf(stderr, "\x1b[31;46mRED_CYAN\x1b[0m\n");
	fprintf(stderr, "\x1b[31;47mRED_WHITE\x1b[0m\n");

    return(0);
}
