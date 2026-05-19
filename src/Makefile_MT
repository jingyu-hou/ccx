CFLAGS = -Wall -O3 -fopenmp -fno-pie -I ../SPOOLES.2.2 -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE -DUSE_MT=1 
FFLAGS = -Wall -O3 -fopenmp -fallow-argument-mismatch 

CC=cc
FC=gfortran

.c.o :
	$(CC) $(CFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<

include Makefile.inc

SCCXMAIN = WeICME.c

OCCXF = $(SCCXF:.f=.o)
OCCXC = $(SCCXC:.c=.o)
OCCXMAIN = $(SCCXMAIN:.c=.o)

DIR=../SPOOLES.2.2

LIBS = \
       $(DIR)/MT/src/spoolesMT.a \
       $(DIR)/spooles.a \
	../ARPACK/libarpack_REDHAT.a \
	../lapack-3.2.1/liblapackREDHAT.a \
	../lapack-3.2.1/libblasREDHAT.a \
       -lpthread -lm -lc

WeICME_MT: $(OCCXMAIN) WeICME_MT.a  $(LIBS)
	./date.pl; $(CC) $(CFLAGS) -c WeICME.c; $(FC) -fopenmp -Wall -O3 -no-pie -o $@ $(OCCXMAIN) WeICME_MT.a $(LIBS)

WeICME_MT.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?
                                                                               
clean:
	rm *.o *.a
