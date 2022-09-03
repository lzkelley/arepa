EXEC   = hop.so

SRCS   = EnzoHop.c hop_hop.c hop_kd.c hop_regroup.c hop_slice.c hop_smooth.c


OBJS   = $(SRCS:.c=.o)
INCL   = hop.h kd.h slice.h smooth.h


CFLAGS =  -shared -fPIC -O3 -g  #  -Wall
LNKCMD =  ld -L/usr/lib -L/usr/local/lib  -shared

LIBS   =  -lm 

CC     =  gcc 

$(EXEC): $(OBJS) 
	$(CC) -L/usr/lib -L/usr/local/lib -shared $(OBJS) $(LIBS) $(CFLAGS)  -o $(EXEC)


$(OBJS): $(INCL) 



clean:
	rm $(OBJS) $(EXEC)
