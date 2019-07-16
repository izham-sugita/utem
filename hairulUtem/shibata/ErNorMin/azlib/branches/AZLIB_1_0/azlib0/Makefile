#####################################################################
# Makefile(azlib0)
#
# Copyright (C) 2001 AzLib Developers Group
#
# Written by
#  <Kenzen TAKEUCHI>
#
# Thanks to following contributors.
#  <Ryu SASAOKA>
#
# Refer to log documents for details.
#####################################################################

# $Id: Makefile,v 1.16 2004/02/05 10:28:57 sasaoka Exp $

## gcc on Linux ##
CC      = gcc
CFLAGS  = -O -Wall -fno-builtin
AR      = ar
ARFLAGS = crsv
RANLIB  = ranlib
INC     =

## Intel C compiler on Linux ##
#CC      = icc
#CFLAGS  = -O -w1 -nolib_inline
#AR      = ar
#ARFLAGS = crsv
#RANLIB  = ranlib
#INC     =


ALL_OBJ = \
  bin_io.o \
  idc.o \
  integration.o \
  mail_utl.o \
  math_utl.o \
  pokecom.o \
  rc.o \
  string_utl.o

all : $(ALL_OBJ) azlib0.a

azlib0.a : $(ALL_OBJ)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

.c.o :
	$(CC) $(CFLAGS) $(INC) -c $<

clean :
	rm -f *.obj *.o *.lib *.a *~ *.bak

