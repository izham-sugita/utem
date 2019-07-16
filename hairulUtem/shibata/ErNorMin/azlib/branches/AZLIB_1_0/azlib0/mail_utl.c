/*********************************************************************
 * mail_utl.c
 *
 * Copyright (C) 2001 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA> <Kenzen TAKEUCHI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: mail_utl.c,v 1.3 2003/07/22 10:21:27 adachi Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include "rc.h"


const char *FileName = "/tmp/aaa.txt";


RC send_mail(char *address, char *text)
{
	char com[2048];
	int out_size;
	FILE *fp;

	RC_NULL_CHK( fp = fopen(FileName, "w") );
	fprintf(fp, "%s\n", text);
	fclose(fp);

#ifndef WIN32
	out_size = snprintf(com, sizeof(com), "mail -s subj %s <%s",
	                    address, FileName);
#else  /* WIN32 */
	out_size = sprintf(com, "mail -s subj %s <%s", address, FileName);
#endif /* WIN32 */

	if((out_size < 1)||(out_size >= sizeof(com))){
		return(OVERFLOW_ERROR_RC);
	}
	system(com);

	if(remove(FileName) != 0){
		return(UNKNOWN_ERROR_RC);
	}

	return(NORMAL_RC);
}


