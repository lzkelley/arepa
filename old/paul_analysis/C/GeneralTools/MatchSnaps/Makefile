EXEC   = MatchSnaps 

OBJS   = main.o mymalloc.o read_parameters.o allvars.o io.o hdf5_read_routines.o

INCL   = allvars.h  proto.h  Makefile

OPTIONS  += -DSNAP_IN_HDF5
OPTIONS  += -DSUBFIND_IN_HDF5
OPTIONS  += -DLONGIDS
OPTIONS  += -DALPHA=0.0                


SYSTYPE="odyssey-opteron"
SYSTYPE="odyssey2"

CC       =   gcc    		    # sets the C-compiler (default)
OPTIMIZE =  -g -Wall -fopenmp       # optimization and warning flags (default)
HDF5INCL =  -DH5_USE_16_API
HDF5LIB  =  -lhdf5 -lz


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









