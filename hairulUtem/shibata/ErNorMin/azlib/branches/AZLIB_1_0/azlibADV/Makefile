#####################################################################
# Makefile(azlibADV)
#
# Copyright (C) 2001 AzLib Developers Group
#
# Written by
#  <Kenzen TAKEUCHI> <Yasuo ASAGA>
#
# Thanks to following contributors.
#  <Ryu SASAOKA>
#
# Refer to log documents for details.
#####################################################################

# $Id: Makefile,v 1.10 2004/02/05 10:28:57 sasaoka Exp $

## gcc on Linux ##
CC      = gcc
CFLAGS  = -O -Wall -fno-builtin
INC     = -I../azlib0 -I../azlibFEM -I/usr/local/Adventure/bin

## Intel C compiler on Linux ##
#CC      = icc
#CFLAGS  = -O -w1 -nolib_inline
#INC     = -I../azlib0 -I../azlibFEM -I/usr/local/Adventure/bin

ALL_OBJ = adv_component.o

all : $(ALL_OBJ)

.c.o :
	$(CC) $(CFLAGS) $(INC) `$(ADVSYSDIR)/advsys-config --cflags` -c $<

clean :
	rm -f *.obj *.o *.lib *.a *~ *.bak
