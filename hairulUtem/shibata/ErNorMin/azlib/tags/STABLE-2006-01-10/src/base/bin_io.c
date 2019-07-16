/*********************************************************************
 * bin_io.c
 *
 * Copyright (C) 2001 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Ryu SASAOKA> 
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: bin_io.c 352 2004-07-26 01:52:25Z sasaoka $ */

#include <stdio.h>
#include "bin_io.h"
#include "rc.h"


RC
read_2bytes (FILE *fp, ENDIAN endian, unsigned int *ui)
{
	if(endian == LH){
		RC_TRY( lread_2bytes(fp, ui) );
	}else if(endian == HL){
		RC_TRY( hread_2bytes(fp, ui) );
	}else{
		return(ARG_ERROR_RC);
	}

	return(NORMAL_RC);
}


RC
read_4bytes (FILE *fp, ENDIAN endian, unsigned long *ul)
{
	if(endian == LH){
		RC_TRY( lread_4bytes(fp, ul) );
	}else if(endian == HL){
		RC_TRY( hread_4bytes(fp, ul) );
	}else{
		return(ARG_ERROR_RC);
	}

	return(NORMAL_RC);
}


/* read 2 bytes and convert to unsigned int value */
RC
hread_2bytes (FILE *fp,unsigned int *value)
{
	unsigned int b;

	*value = 0;
	if((b = fgetc(fp)) == EOF)return(READ_ERROR_RC);
	*value = (b << 8);
	if((b = fgetc(fp)) == EOF)return(READ_ERROR_RC);
	*value += b;

	return(NORMAL_RC);
}


/* read 4 bytes and convert to unsigned long value */
RC
hread_4bytes (FILE *fp,unsigned long *value)
{
	unsigned long b;

	*value = 0;
	if((b = (unsigned long)fgetc(fp)) == EOF)return(READ_ERROR_RC);
	*value = (b << 24);
	if((b = (unsigned long)fgetc(fp)) == EOF)return(READ_ERROR_RC);
	*value += (b << 16);
	if((b = (unsigned long)fgetc(fp)) == EOF)return(READ_ERROR_RC);
	*value += (b << 8);
	if((b = (unsigned long)fgetc(fp)) == EOF)return(READ_ERROR_RC);
	*value += b;

	return(NORMAL_RC);
}


/* write 2 bytes */
RC
hwrite_2bytes (FILE *fp,unsigned int value)
{
	int c;

	c = (0xff & (value >> 8));
	if(fputc(c,fp) == EOF)return(WRITE_ERROR_RC);
	c = (0xff & value);
	if(fputc(c,fp) == EOF)return(WRITE_ERROR_RC);

	return(NORMAL_RC);
}


/* write 4 bytes */
RC
hwrite_4bytes (FILE *fp,unsigned long value)
{
	int c;

	c = (0xff & (value >> 24));
	if(fputc(c,fp) == EOF)return(WRITE_ERROR_RC);
	c = (0xff & (value >> 16));
	if(fputc(c,fp) == EOF)return(WRITE_ERROR_RC);
	c = (0xff & (value >> 8));
	if(fputc(c,fp) == EOF)return(WRITE_ERROR_RC);
	c = (0xff & value);
	if(fputc(c,fp) == EOF)return(WRITE_ERROR_RC);

	return(NORMAL_RC);
}


/* read 2 bytes and convert to unsigned int value */
RC
lread_2bytes (FILE *fp,unsigned int *value)
{
	unsigned int b;

	*value = 0;
	if((b = fgetc(fp)) == EOF)return(READ_ERROR_RC);
	*value = b;
	if((b = fgetc(fp)) == EOF)return(READ_ERROR_RC);
	*value += (b << 8);

	return(NORMAL_RC);
}


/* read 4 bytes and convert to unsigned long value */
RC
lread_4bytes (FILE *fp,unsigned long *value)
{
	unsigned long b;

	*value = 0;
	if((b = (unsigned long)fgetc(fp)) == EOF)return(READ_ERROR_RC);
	*value = b;
	if((b = (unsigned long)fgetc(fp)) == EOF)return(READ_ERROR_RC);
	*value += (b << 8);
	if((b = (unsigned long)fgetc(fp)) == EOF)return(READ_ERROR_RC);
	*value += (b << 16);
	if((b = (unsigned long)fgetc(fp)) == EOF)return(READ_ERROR_RC);
	*value += (b << 24);

	return(NORMAL_RC);
}


/* write 2 bytes */
RC
lwrite_2bytes (FILE *fp,unsigned int value)
{
	int c;

	c = (0xff & value);
	if(fputc(c,fp) == EOF)return(WRITE_ERROR_RC);
	c = (0xff & (value >> 8));
	if(fputc(c,fp) == EOF)return(WRITE_ERROR_RC);

	return(NORMAL_RC);
}


/* write 4 bytes */
RC
lwrite_4bytes (FILE *fp,unsigned long value)
{
	int c;

	c = (0xff & value);
	if(fputc(c,fp) == EOF)return(WRITE_ERROR_RC);
	c = (0xff & (value >> 8));
	if(fputc(c,fp) == EOF)return(WRITE_ERROR_RC);
	c = (0xff & (value >> 16));
	if(fputc(c,fp) == EOF)return(WRITE_ERROR_RC);
	c = (0xff & (value >> 24));
	if(fputc(c,fp) == EOF)return(WRITE_ERROR_RC);

	return(NORMAL_RC);
}


