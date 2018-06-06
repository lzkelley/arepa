#OPT    += -DDIRECT_RAY_INTEGRATION

EXEC   = getnh.so

SRCS   = idlmakedustmap.c allocate.c printhead.c system.c setup_units.c \
	 localvalues.c cooling.c tree.c readsnap.c overhead.c \
	 ngbtree3d.c nrutil.c selectb.c no_tree.c

OBJS   = $(SRCS:.c=.o)
INCL   = overhead.h   \
	 nrutil.h

CFLAGS =  -g -Wall -fPIC #-O3 -g  #  -Wall
LNKCMD =  ld -L/usr/lib -L/usr/local/lib  -shared
#CFLAGS =  -shared -m32 $(OPT)
#CFLAGS =  -bundle $(OPT)

LIBS   = -lm  

CC     =  gcc 

#$(EXEC): $(OBJS) 
#	$(CC) $(CFLAGS) $(OBJS) $(LIBS)  -o $(EXEC)

#$(OBJS): $(SRCS) 

$(EXEC): $(OBJS) 
	$(CC) -L/usr/lib -L/usr/local/lib  -shared $(OBJS) $(LIBS) $(CFLAGS)  -o $(EXEC)

$(OBJS): $(SRCS)
#$(OBJS): $(INCL) 



clean:
	rm $(OBJS) $(EXEC)










