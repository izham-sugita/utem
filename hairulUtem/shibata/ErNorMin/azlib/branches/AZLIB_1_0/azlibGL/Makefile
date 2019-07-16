#####################################################################
# Makefile(azlibGL)
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

# $Id: Makefile,v 1.11 2004/02/05 10:28:57 sasaoka Exp $

## gcc on Linux ##
CC      = gcc
CFLAGS  = -O -Wall -fno-builtin
AR      = ar
ARFLAGS = crsv
RANLIB  = ranlib
INC     = -I../azlib0 -I../azlibFEM -I../azlibPIX

## Intel C compiler on Linux ##
#CC      = icc
#CFLAGS  = -O -w1 -nolib_inline
#AR      = ar
#ARFLAGS = crsv
#RANLIB  = ranlib
#INC     = -I../azlib0 -I../azlibFEM -I../azlibPIX

## gtk+-2.0 & gtkglext-1.0 ##
GTK_CFLAGS = `/usr/bin/pkg-config --cflags gtk+-2.0 gtkglext-1.0`


ALL_OBJ =  \
  glut_utl.o \
  gtkgl_utl.o \
  gtk_graph.o \
  gtk_model_view.o \
  gtk_utl.o 

all : $(ALL_OBJ) azlibGL.a

azlibGL.a : $(ALL_OBJ)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

.c.o :
	$(CC) $(CFLAGS) $(INC) $(GTK_CFLAGS) -c $<

clean :
	rm -f *.obj *.o *.lib *.a *~ *.bak

