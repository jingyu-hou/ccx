# export We_HOME=$(shell pwd)
# ARDIR = ./ARPACK
# SPDIR = ./SPOOLES.2.2
# SPMTDIR = ./SPOOLES.2.2/MT/src
# CCXDIR = ./src
# MAKE = make
# MAKELIB = make lib
# subsystem:
# 	cd $(ARDIR) && $(MAKELIB)
# 	cd $(SPDIR) && $(MAKELIB)
# 	cd $(SPMTDIR) && $(MAKE)
# 	cd $(CCXDIR) && $(MAKE)


# Prool's makefile
all:
	cd SPOOLES.2.2/MT;make lib
	cp SPOOLES.2.2/MT/src/spoolesMT.a SPOOLES.2.2/spoolesMT.a
	cp src/Makefile_MT src/Makefile
	cd ARPACK;make home=`pwd` lib
	cd SPOOLES.2.2;make lib
	cd src;make
clean:
	cd ARPACK;make clean
	cd SPOOLES.2.2;make clean
	cd src;make clean
