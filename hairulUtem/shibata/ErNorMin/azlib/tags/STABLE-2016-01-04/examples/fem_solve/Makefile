## CC,CFLAGS,LIBS ##
## GCC ##
#CC     = gcc
#CC     = gcc -pg  # with gprof
CC      = gcc -g  # with GDB
CFLAGS  = -O -Wall -fno-builtin
#LIBS   = -lm -lpthread
LIBS    = -lm

## Intel C compiler on Linux ##
#CC      = icc
#CFLAGS  = -O -w1 -nolib_inline

## SGI IRIX ##
#CC      = cc -64
#LIBS    = -lfastm
#CFLAGS  = -O -ansi

## Tru64 ##
#CC      = cc
#LIBS    = -lm
#CFLAGS  = -O -fast -msg_disable longdoubleny1


## Common ##
MAKE    = make "CC=$(CC)" "CFLAGS=$(CFLAGS)"
RM      = rm -f

## AzLib ##
AZLIB  = ../../src
BASE   = $(AZLIB)/base
FEM    = $(AZLIB)/fem
GUI    = $(AZLIB)/gui
IMAGE  = $(AZLIB)/image
MATH   = $(AZLIB)/mathematics

TARGET0 = fem_solve
OBJM0   = fem_solve.o\
          file_io.o
LIBAZ0  = $(FEM)/fem.a \
          $(MATH)/mathematics.a \
          $(BASE)/base.a

ALL = $(TARGET0)

all : $(ALL)

base_obj :
	cd $(BASE) ; $(MAKE) all

fem_obj :
	cd $(FEM) ; $(MAKE) all

math_obj :
	cd $(MATH) ; $(MAKE) all


$(TARGET0) : $(OBJM0) base_obj math_obj fem_obj
	$(CC) -o $@ $(OBJM0) $(LIBAZ0) $(LIBS)

.c.o :
	$(CC) $(CFLAGS) -I$(BASE) -I$(MATH) -I$(FEM) -c $<
 
clean :
	cd $(BASE);  $(MAKE) clean
	cd $(FEM);   $(MAKE) clean
	cd $(MATH); $(MAKE) clean
	$(RM) $(ALL) *.o *.obj *.bak gmon.out core

