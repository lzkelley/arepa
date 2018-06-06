EXEC   = B-HaloTrees
OBJS   = main.o mymalloc.o read_parameters.o allvars.o io.o mergertree.o hdf5_read.o
INCL   = allvars.h  proto.h  Makefile

OPTIONS  += -DLONGIDS
#OPTIONS  += -DSUBFIND_EXTRA      # read/include gas and GFM related group catalog fields

OPTIONS += -DSKIP_SNAP

SYSTYPE="odyssey-opteron"
SYSTYPE="odyssey2"
SYSTYPE="darwin"

CC       =   gcc                  # sets the C-compiler (default)
OPTIMIZE =   -g -Wall  -fopenmp   # optimization and warning flags (default)
HDF5INCL =  -DH5_USE_16_API
HDF5LIB  =  -lhdf5 -lz

ifeq ($(SYSTYPE),"darwin")
CC       = mpicc
OPTIMIZE = -g -Wall
HDF5INCL = -I/opt/local/inlcude -DH5_USE_16_API
HDF5LIB = -L/opt/local/lib -lhdf5 -lz
endif

ifeq ($(SYSTYPE),"odyssey2")
CC       =   mpicc
OPTIMIZE =   -g -Wall
endif

ifeq ($(SYSTYPE),"odyssey-opteron")
CC       =  icc                   # sets the C-compiler (default)
OPTIMIZE =  -O3 -Wall             # optimization and warning flags (default)
HDF5INCL =  -DH5_USE_16_API
HDF5LIB  =  -lhdf5
endif


LIBS   =   -lm $(HDF5LIB) -g
CFLAGS =   -g $(OPTIONS)  $(OPTIMIZE)  $(HDF5INCL)

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS)  $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 

clean:
	rm -f $(OBJS) $(EXEC)







