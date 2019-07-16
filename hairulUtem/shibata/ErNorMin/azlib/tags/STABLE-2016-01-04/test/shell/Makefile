## CC,CFLAGS,LIBS ##
## GCC ##
CC     = gcc
#CC     = gcc -pg  # with gprof
#CC      = gcc -g  # with GDB
CFLAGS  = -O -Wall -fno-builtin -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
#LIBS   = -lm -lpthread
LIBS    = -lm

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
GUI    = $(AZLIB)/gui
IMAGE  = $(AZLIB)/image
MATH   = $(AZLIB)/mathematics


TARGET0 = fem_shell
OBJM0   = $(TARGET0).o
LIBAZ0  = $(FEM)/fem.a \
          $(MATH)/mathematics.a \
          $(BASE)/base.a

TARGET1 = amls_nonzero3U
OBJM1   = $(TARGET1).o
LIBAZ1  = $(FEM)/fem.a \
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

$(TARGET0) : $(OBJM0) fem_obj base_obj math_obj
	$(CC) -o $@ $(OBJM0) $(LIBAZ0) $(LIBS)

$(TARGET1) : $(OBJM1) fem_obj base_obj math_obj
	$(CC) -o $@ $(OBJM1) $(LIBAZ1) $(LIBS)

.c.o :
	$(CC) $(CFLAGS) -I$(BASE) -I$(MATH) -I$(FEM) -c $<
 
#clean :
#	rm -f $(ALL) *.o *.obj *.bak *~ gmon.out core

clean :
	cd $(BASE);  $(MAKE) clean
	cd $(FEM);   $(MAKE) clean
	cd $(GUI);   $(MAKE) clean
	cd $(IMAGE); $(MAKE) clean
	cd $(MATH); $(MAKE) clean
	$(RM) $(ALL) *.o *.obj *.bak *~ gmon.out core

