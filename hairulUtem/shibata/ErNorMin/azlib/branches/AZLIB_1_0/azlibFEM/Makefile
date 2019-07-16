#####################################################################
# Makefile(azlibFEM)
#
# Copyright (C) 2001 AzLib Developers Group
#
# Written by
#  <Kenzen TAKEUCHI>
#
# Thanks to following contributors.
#  <Yasuo ASAGA> <Takaaki NAGATANI>
#
# Refer to log documents for details.
#####################################################################

# $Id: Makefile,v 1.15 2004/02/05 10:28:57 sasaoka Exp $

## gcc on Linux ##
CC      = gcc
CFLAGS  = -O -Wall -fno-builtin
AR      = ar
ARFLAGS = crsv
RANLIB  = ranlib
INC     = -I../azlib0

## Intel C compiler on Linux ##
#CC      = icc
#CFLAGS  = -O -w1 -nolib_inline
#AR      = ar
#ARFLAGS = crsv
#RANLIB  = ranlib
#INC     = -I../azlib0

ALL_OBJ = \
  ansys_component.o \
  avs_utl.o \
  dxf_utl.o \
  eigen_solver.o \
  element_volume.o \
  eps_utl.o \
  extract_surface.o \
  fem_renumber.o \
  fem_solver.o \
  fem_solver_6dof.o \
  nl_solver.o \
  fem_struct.o \
  fem_struct_binio.o \
  fem_struct_print.o \
  fem_utl.o \
  flow_solver.o \
  interpolation.o \
  nonzero_cg.o \
  nst_blk_component.o \
  nst_blk_llio.o \
  nst_f06_component.o \
  nst_pch_component.o \
  optimization_utl.o \
  pressure_force.o \
  sky_cholesky.o \
  sky_crout.o \
  stl_utl.o \
  sysnoise_component.o

all : $(ALL_OBJ) azlibFEM.a

azlibFEM.a : $(ALL_OBJ)
	$(AR) $(ARFLAGS) $@ $? 
	$(RANLIB) $@

.c.o :
	$(CC) $(CFLAGS) $(INC) -c $<

clean :
	rm -f *.obj *.o *.lib *.a *~ *.bak

