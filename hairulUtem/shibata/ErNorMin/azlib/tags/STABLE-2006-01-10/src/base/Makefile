#####################################################################
# Makefile(base)
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

# $Id: Makefile 420 2005-08-12 07:02:40Z sasaoka $

## gcc on Linux ##
CC      = gcc
CFLAGS  = -O -Wall -fno-builtin
AR      = ar
ARFLAGS = crs
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
  data_model.o \
  idc.o \
  list.o \
  log_printf.o \
  mail_utl.o \
  memory_manager.o \
  rc.o \
  scratch_io.o \
  string_utl.o \
  time_utl.o \
  wrapper_64.o

all : $(ALL_OBJ) base.a

base.a : $(ALL_OBJ)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

.c.o :
	$(CC) $(CFLAGS) $(INC) -c $<

clean :
	rm -f *.obj *.o *.lib *.a *~ *.bak

