#####################################################################
# Makefile(mathematics)
#
# Copyright (C) 2005 AzLib Developers Group
#
# Written by
#  <Ryu SASAOKA>
#
# Thanks to following contributors.
#  <Kenzen TAKEUCHI>
#
# Refer to log documents for details.
#####################################################################

# $Id: Makefile 1095 2016-12-05 07:29:35Z hayashi $

## gcc on Linux ##
CC      = gcc
CFLAGS  = -O -Wall -fno-builtin
AR      = ar
ARFLAGS = crs
RANLIB  = ranlib
INC     = -I../base

## Intel C compiler on Linux ##
#CC      = icc
#CFLAGS  = -O -w1 -fno-builtin
#AR      = ar
#ARFLAGS = crsv
#RANLIB  = ranlib
#INC     = -I../base


ALL_OBJ = \
  complex_utl.o \
  crs_complex_utl.o \
  geometry.o \
  lq.o \
  lu.o \
  math_utl.o \
  matrix.o \
  nonzero_cg.o \
  operation.o \
  pokecom2.o \
  qr.o \
  random.o \
  vect3d.o \
  vect3d_complex.o \
  vect_storage.c

all : $(ALL_OBJ) mathematics.a

mathematics.a : $(ALL_OBJ)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

.c.o :
	$(CC) $(CFLAGS) $(INC) -c $<

clean :
	rm -f *.obj *.o *.lib *.a *~ *.bak

