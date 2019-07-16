#####################################################################
# Makefile(azlibPIX)
#
# Copyright (C) 2001 AzLib Developers Group
#
# Written by
#  <Kenzen TAKEUCHI> <Ryu SASAOKA>
#
# Thanks to following contributors.
#  <Yasuo ASAGA>
#
# Refer to log documents for details.
#####################################################################

# $Id: Makefile,v 1.11 2004/02/06 08:06:51 sasaoka Exp $

## gcc ##
CC      = gcc
CFLAGS  = -O -Wall -fno-builtin
AR      = ar
ARFLAGS = crsv
RANLIB  = ranlib
INC     = -I../azlib0 -I../azlibFEM

## icc for Linux ##
#CC      = icc
#CFLAGS  = -O -w1 -nolib_inline
#AR      = ar
#ARFLAGS = crsv
#RANLIB  = ranlib
#INC     = -I../azlib0 -I../azlibFEM


ALL_OBJ = \
  avs_field_utl.o \
  bmp_utl.o \
  dicom_utl.o \
  graphic_utl.o \
  tiff_utl.o

all : $(ALL_OBJ) azlibPIX.a

azlibPIX.a : $(ALL_OBJ)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

.c.o :
	$(CC) $(CFLAGS) $(INC) -c $<

clean :
	rm -f *.obj *.o *.lib *.a *~ *.bak

