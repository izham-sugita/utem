#####################################################################
# Makefile(image)
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

# $Id: Makefile 928 2009-10-06 11:04:42Z aoyama $

## gcc ##
CC      = gcc
CFLAGS  = -O -Wall -fno-builtin
AR      = ar
ARFLAGS = crs
RANLIB  = ranlib
INC     = -I../base -I../mathematics -I../fem

## icc for Linux ##
#CC      = icc
#CFLAGS  = -O -w1 -fno-builtin
#AR      = ar
#ARFLAGS = crsv
#RANLIB  = ranlib
#INC     = -I../base -I../mathematics -I../fem


ALL_OBJ = \
  avs_field_utl.o \
  bmp_utl.o \
  dicom_utl.o \
  graphic_binio.o \
  graphic_utl.o \
  tiff_utl.o \
  voxel_utl.o

all : $(ALL_OBJ) image.a

image.a : $(ALL_OBJ)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

.c.o :
	$(CC) $(CFLAGS) $(INC) -c $<

clean :
	rm -f *.obj *.o *.lib *.a *~ *.bak

