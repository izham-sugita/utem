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

# $Id: Makefile 455 2005-09-16 06:36:00Z sasaoka $

# gcc on Linux ##
CC      = gcc
CFLAGS  = -O -Wall -fno-builtin
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
GUI    = $(SRC)/gui
IMAGE  = $(SRC)/image


all : base_obj math_obj fem_obj gui_obj image_obj

base_obj :
	cd $(BASE) ; $(MAKE) all

math_obj :
	cd $(MATH) ; $(MAKE) all

fem_obj :
	cd $(FEM) ; $(MAKE) all

gui_obj :
	cd $(GUI) ; $(MAKE) all

image_obj :
	cd $(IMAGE) ; $(MAKE) all

global :
	rm -rf GPATH GRTAGS GSYMS GTAGS HTML/
	gtags -v
	htags -vn

clean :
	cd $(BASE);  $(MAKE) clean
	cd $(MATH);  $(MAKE) clean
	cd $(FEM);   $(MAKE) clean
	cd $(GUI);   $(MAKE) clean
	cd $(IMAGE); $(MAKE) clean
	rm -f GPATH GRTAGS GSYMS GTAGS *.bak

