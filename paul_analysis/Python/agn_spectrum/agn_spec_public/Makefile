
EXEC   = agn_spectrum_shared.so

SRCS   = agn_spectrum_shared.c


OBJS   = $(SRCS:.c=.o)
INCL   = 


#CFLAGS =  -shared -m32 $(OPT)
#CFLAGS =  -bundle $(OPT)
# next lines for OS X compile only; otherwise, use one of the above two lines
SDK    = /Developer/SDKs/MacOSX10.4u.sdk
CFLAGS  =  -bundle $(OPT) -isysroot $(SDK) -arch ppc -arch ppc64 -arch i386
LDFLAGS = -isysroot $(SDK) -syslibroot,$(SDK)

LIBS   =   

CC     =  gcc 

$(EXEC): $(OBJS) 
	$(CC) $(CFLAGS) $(OBJS) $(LIBS)  -o $(EXEC)


$(OBJS): $(SRCS) 



clean:
	rm $(OBJS) $(EXEC)










