#####################################################################
# Makefile(fem)
#
# Copyright (C) 2005 AzLib Developers Group
#
# Written by
#  <Kenzen TAKEUCHI>
#
# Thanks to following contributors.
#  <Yasuo ASAGA> <Takaaki NAGATANI> <Ryu SASAOKA>
#
# Refer to log documents for details.
#####################################################################

# $Id: Makefile 1009 2012-03-03 13:48:01Z aoyama $

## gcc on Linux ##
CC      = gcc
CFLAGS  = -O -Wall -fno-builtin
AR      = ar
ARFLAGS = crs
RANLIB  = ranlib
INC     = -I../base -I../mathematics
## with Adoventure
#INC     = -I../base -I../mathematics -I/usr/local/Adventure/bin

## Intel C compiler on Linux ##
#CC      = icc
#CFLAGS  = -O -w1 -fno-builtin
#AR      = ar
#ARFLAGS = crsv
#RANLIB  = ranlib
#INC     = -I../base -I../mathematics 
## with Adoventure
#INC     = -I../base -I../mathematics -I/usr/local/Adventure/bin

ALL_OBJ = \
  ansys_component.o \
  avs_utl.o \
  contact.o \
  dxf_utl.o \
  eigen_solver.o \
  element_volume.o \
  eps_utl.o \
  extract_surface.o \
  face_edge_table.o \
  fem_geometry.o \
  fem_renumber.o \
  fem_solver.o \
  fem_solver_6dof.o \
  fem_solver_shell.o \
  fem_solver_plate_shell.o \
  fem_search_common.o \
  fem_struct.o \
  fem_struct_binio.o \
  fem_struct_print.o \
  fem_utl.o \
  flow_solver.o \
  integration.o \
  interpolation.o \
  matrix_db.o \
  matrix_db_amls.o \
  nl_solver.o \
  nl_solver_shell.o \
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

all : $(ALL_OBJ) fem.a

fem.a : $(ALL_OBJ)
	$(AR) $(ARFLAGS) $@ $? 
	$(RANLIB) $@

.c.o :
	$(CC) $(CFLAGS) $(INC) -c $<

clean :
	rm -f *.obj *.o *.lib *.a *~ *.bak

