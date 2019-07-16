################################################################################
## Makefile for fem_view
##
## $Id: Makefile 1072 2016-09-13 09:22:56Z hayashi $
##
################################################################################

## CC,CFLAGS,LIBS ##
## GCC ##
CC     = gcc
#CC     = gcc -pg  # with gprof
#CC      = gcc -g  # with GDB
CFLAGS  = -O -Wall -fno-builtin
GTK_CFLAGS = `/usr/bin/pkg-config --cflags gtk+-2.0 gtkglext-1.0`
#LIBS   = -lm -lpthread
LIBS    = -lm
GTK_LIBS = `/usr/bin/pkg-config --libs gtk+-2.0 gtkglext-1.0`
OPENGL_LIBS = -lGL -lGLU -lglut

## Intel C compiler on Linux ##
#CC      = icc
#CFLAGS  = -O -w1 -nolib_inline
#LIBS    = -lm

## SGI IRIX ##
#CC      = cc -64
#CFLAGS  = -O -ansi
#LIBS    = -lfastm

## Tru64 ##
#CC      = cc
#CFLAGS  = -O -fast -msg_disable longdoubleny1
#LIBS    = -lm

## Common ##
MAKE    = make "CC=$(CC)" "CFLAGS=$(CFLAGS)"
RM      = rm -f

## AzLib ##
AZLIB    = ../../src
BASE   = $(AZLIB)/base
FEM    = $(AZLIB)/fem
MATH   = $(AZLIB)/mathematics
IMAGE  = $(AZLIB)/image
GUI    = $(AZLIB)/gui


TARGET0 = main
OBJM0   = $(TARGET0).o
LIBAZ0  = $(GUI)/gui.a \
          $(IMAGE)/image.a \
          $(FEM)/fem.a \
          $(MATH)/mathematics.a \
          $(BASE)/base.a

ALL = $(TARGET0)
all : $(ALL)

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

$(TARGET0) : $(OBJM0) base_obj math_obj fem_obj image_obj gui_obj
	$(CC) -o $@ $(OBJM0) $(LIBAZ0) $(LIBS) $(GTK_LIBS) $(OPENGL_LIBS)

.c.o :
	$(CC) $(CFLAGS) $(GTK_CFLAGS) -I$(BASE) -I$(MATH) -I$(FEM) -I$(IMAGE) -I$(GUI) -c $<
 
#clean :
#	rm -f $(ALL) *.o *.obj *.bak *~ gmon.out core

clean :
	$(RM) $(ALL) *.o *.obj *.bak *~ gmon.out core

cleanall :
	cd $(BASE);  $(MAKE) clean
	cd $(MATH);  $(MAKE) clean
	cd $(FEM);   $(MAKE) clean
	cd $(IMAGE); $(MAKE) clean
	cd $(GUI);   $(MAKE) clean
	$(MAKE) clean
	
