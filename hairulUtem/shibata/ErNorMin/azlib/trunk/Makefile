#####################################################################
# Makefile(azlib)
#
# Copyright (C) 2004 AzLib Developers Group
#
# Written by
#  <Ryu SASAOKA>
#
# Thanks to following contributors.
#  <Kenzen TAKEUCHI>
#
# Refer to log documents for details.
#####################################################################

# $Id: Makefile 1035 2016-01-21 12:55:58Z hayashi $

# gcc on Linux ##
CC      = gcc
CFLAGS  = -O -Wall -fno-builtin
#CFLAGS  = -O -Wall -fno-builtin  -pedantic
#CFLAGS  = -O -Wall -fno-builtin  -pedantic -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
RM      = rm -f
MAKE    = make "CC=$(CC)" "CFLAGS=$(CFLAGS)"

## Intel C compiler on Linux ##
#CC      = icc
#CFLAGS  = -O -w1 -nolib_inline
#RM      = rm -f
#MAKE    = make "CC=$(CC)" "CFLAGS=$(CFLAGS)"

## AZLIB ##
SRC    = src
BASE   = $(SRC)/base
MATH   = $(SRC)/mathematics
FEM    = $(SRC)/fem
IMAGE  = $(SRC)/image
GUI    = $(SRC)/gui


all : base_obj math_obj fem_obj image_obj gui_obj 

base_obj :
	cd $(BASE) ; $(MAKE) all

math_obj :
	cd $(MATH) ; $(MAKE) all

fem_obj :
	cd $(FEM) ; $(MAKE) all

image_obj :
	cd $(IMAGE) ; $(MAKE) all

gui_obj :
	cd $(GUI) ; $(MAKE) all

global :
	rm -rf GPATH GRTAGS GSYMS GTAGS HTML/
	gtags -v
	htags -vn

clean :
	cd $(BASE);  $(MAKE) clean
	cd $(MATH);  $(MAKE) clean
	cd $(FEM);   $(MAKE) clean
	cd $(IMAGE); $(MAKE) clean
	cd $(GUI);   $(MAKE) clean
	rm -f GPATH GRTAGS GSYMS GTAGS *.bak

