#####################################################################
# Makefile(gui)
#
# Copyright (C) 2002 AzLib Developers Group
#
# Written by
#  <Kenzen TAKEUCHI> <Ryu SASAOKA>
#
# Thanks to following contributors.
#  <Yasuo ASAGA>
#
# Refer to log documents for details.
#####################################################################

# $Id: Makefile 420 2005-08-12 07:02:40Z sasaoka $

## gcc on Linux ##
CC      = gcc
CFLAGS  = -O -Wall -fno-builtin
AR      = ar
ARFLAGS = crs
RANLIB  = ranlib
INC     = -I../base -I../mathematics -I../fem -I../image

## Intel C compiler on Linux ##
#CC      = icc
#CFLAGS  = -O -w1 -nolib_inline
#AR      = ar
#ARFLAGS = crsv
#RANLIB  = ranlib
#INC     = -I../base -I../mathematics  -I../fem -I../image

## gtk+-2.0 & gtkglext-1.0 ##
GTK_CFLAGS = `/usr/bin/pkg-config --cflags gtk+-2.0 gtkglext-1.0`


ALL_OBJ = \
  glut_utl.o \
  gtk_model_view.o \
  gtk_utl.o \
  gtkgl_article.o \
  gtkgl_fem.o \
  gtkgl_image.o \
  gtkgl_utl.o

all : $(ALL_OBJ) gui.a

gui.a : $(ALL_OBJ)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

.c.o :
	$(CC) $(CFLAGS) $(INC) $(GTK_CFLAGS) -c $<

clean :
	rm -f *.obj *.o *.lib *.a *~ *.bak


