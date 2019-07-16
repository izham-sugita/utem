## gcc ##
#CC      = gcc
#LIBS    = -L/home/aoyama/box/GotoBLAS -lgoto -lpthread -L/usr/lib/gcc/i686-momonga-linux/4.1.0 -lgfortranbegin -lgfortran -lm
#LIBS    = -L/home/aoyama/box/GotoBLAS -lgoto -lpthread -lgfortranbegin -lgfortran -lm
#CFLAGS  = -O -Wall -fno-builtin

## gcc and gprof ##
CC     = gcc -pg
#static link
#LIBS    = -lm

#dynamic link
#LIBS    = -lm -ldmumps-4.10.0 -lmumps_common-4.10.0 -lmpiseq_seq-4.10.0 -lpord-4.10.0 -lgfortranbegin -lgfortran
#LIBS    = -lm -ldmumps_seq -lmumps_common_seq -lpord_seq -lgfortranbegin -lgfortran -lpthread
LIBS    = -lm -ldmumps_seq -lmumps_common_seq -lpord_seq -lgfortran -lpthread
CFLAGS  = -O -Wall -fno-builtin -fopenmp

## icc ##
#CC      = icc
#LIBS    =  -lifcore -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_gnu_thread -lpthread 
#CFLAGS  = -O -w1 -nolib_inline -openmp -openmp-report1

## SGI IRIX ##
#CC      = cc -64
#LIBS    = -lfastm
#CFLAGS  = -O -ansi

## Tru64 ##
#CC      = cc
#LIBS    = -lm
#CFLAGS  = -O -fast -msg_disable longdoubleny1

AZLIB    = ./azlib/trunk/src/
BASE     = $(AZLIB)/base
FEM      = $(AZLIB)/fem
GUI      = $(AZLIB)/gui
IMAGE    = $(AZLIB)/image
MATH     = $(AZLIB)/mathematics
#ADVSYSDIR = /usr/local/Adventure/AdvIO-0.8b

#Need to modify location of the files
#MUMPS       = ./lib/MUMPS_4.10.0
#MUMPS = ./lib/mumps
#MUMPSLIB    = $(MUMPS)/lib

#my MUMPS build
#MUMPS = /home/sugita/cpp-codes/hairulUtem/shibata/ErNorMin/lib/mumps
#MUMPSLIB    = $(MUMPS)/lib
#MUMPSINC    = $(MUMPS)/include
#LIBSEQ      = $(MUMPS)/libseq
#PORDLIB     = $(MUMPS)/PORD/lib
#PORDINC     = $(MUMPS)/PORD/include

#Ubuntu system MUMPS build
MUMPS = /usr
MUMPSLIB    = $(MUMPS)/lib
MUMPSINC    = $(MUMPS)/include
LIBSEQ      = $(MUMPS)/libseq
PORDLIB     = $(MUMPS)/lib
PORDINC     = $(MUMPS)/include

MAKE    = make "CC=$(CC)" "CFLAGS=$(CFLAGS)"

TARGET0 = main
OBJM0   = $(TARGET0).o \
		  nl_solver_topology.o \
          ccs_format.o

LIBAZ0 = $(FEM)/fem.a \
$(BASE)/base.a \
$(MATH)/mathematics.a 
#$(MUMPSLIB)/libdmumps.a \
#$(MUMPSLIB)/libmumps_common.a \
#$(LIBSEQ)/libmpiseq.a \
#$(PORDLIB)/libpord.a


ALL =   $(TARGET0)

all : $(ALL)

base_obj :
	cd $(BASE) ; $(MAKE) all

fem_obj :
	cd $(FEM) ; $(MAKE) all

gui_obj :
	cd $(GUI) ; $(MAKE) all

image_obj :
	cd $(IMAGE) ; $(MAKE) all

math_obj :
	cd $(MATH) ; $(MAKE) all

#mumps_obj :
#	cd $(MUMPS) ; $(MAKE) 

#$(TARGET0) : $(OBJM0) base_obj fem_obj math_obj mumps_obj
#	$(CC) $(CFLAGS) -o $@ $(OBJM0) $(OBJS0) $(LIBAZ0) $(LIBS)

$(TARGET0) : $(OBJM0) base_obj fem_obj math_obj
	$(CC) $(CFLAGS) -o $@ $(OBJM0) $(OBJS0) $(LIBAZ0) $(LIBS)

#original make
#.c.o :
#	$(CC) $(CFLAGS) -I$(BASE) -I$(FEM) -I$(GUI) -I$(IMAGE) -I$(MATH) \
#	-I$(MUMPSINC) -I$(LIBSEQ) -I$(PORDINC) -c $<

.c.o :
	$(CC) $(CFLAGS) -I$(BASE) -I$(FEM) -I$(GUI) -I$(IMAGE) -I$(MATH) \
	-I$(MUMPSINC) -I$(LIBSEQ) -I$(PORDINC) -c $<

clean :
	rm -f $(ALL) *.o *.obj *.bak *.inp gmon.out core

cleanall :
	cd $(BASE);   $(MAKE) clean
	cd $(FEM); $(MAKE) clean
	cd $(GUI);  $(MAKE) clean
	cd $(IMAGE); $(MAKE) clean
	cd $(MATH); $(MAKE) clean
#	cd $(MUMPS); $(MAKE) clean
	rm -f $(ALL) *.o *.bak gmon.out core

