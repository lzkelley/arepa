EXEC   = B-BaseTree

OBJS   = main.o mymalloc.o read_parameters.o allvars.o io.o hdf5_read_routines.o

INCL   = allvars.h  proto.h  Makefile

#OPTIONS  += -DINCLUDE_STARS
OPTIONS  += -DLONGIDS
OPTIONS  += -DIDS_HAVE_GAPS
#OPTIONS  += -DSWAP_ENDIAN
#OPTIONS  += -DSKIP_BY_WEIGHT                          # Do A->C assigning also in case that substructure is assigned to be merged to larger structure in B
OPTIONS  += -DALPHA=0.8                                # Exponent in weighting function
#OPTIONS  += -DBACKWARD_CHECKING                       # Check also backwards (B->A) dependencies and does reassignement of decendent in case of halos appearing 
#OPTIONS  += -DHALO_SIZE_INCREASE_FOR_SWITCHING=1.5    # How much bigger a "newborn" Halo has to be to redirect in case of detected confusion
OPTIONS += -DSKIP_SNAP

SYSTYPE="odyssey-opteron"
SYSTYPE="odyssey2"
SYSTYPE="darwin"

CC       =   gcc    		    # sets the C-compiler (default)
OPTIMIZE =  -g -Wall -fopenmp       # optimization and warning flags (default)
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
CC       =   icc  
OPTIMIZE =   -g -Wall 
endif


LIBS   =   -lm $(HDF5LIB) -g

CFLAGS =   -g $(OPTIONS)  $(OPTIMIZE)  $(HDF5INCL)

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


clean:
	rm -f $(OBJS) $(EXEC)









