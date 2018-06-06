EXEC   = MakeGalaxy

OBJS   = main.o cooling.o read_param.o effmodel.o  globvars.o \
	halo.o haloset.o escape.o distfunc.o \
	disk.o diskset.o \
	bulge.o bulgeset.o \
	gas.o  gasset.o \
	save.o init.o structure.o \
	misc.o toomre.o force.o \
	forcetree.o \
	nrsrc/bessi0.o   nrsrc/bessj1.o   \
	nrsrc/erff.o     nrsrc/polint.o   nrsrc/spline.o \
	nrsrc/bessi1.o   nrsrc/bessk0.o   nrsrc/gasdev.o  \
	nrsrc/qromb.o    nrsrc/splint.o  \
	nrsrc/bessj0.o   nrsrc/bessk1.o   \
	nrsrc/nrutil.o   nrsrc/ran1.o     nrsrc/trapzd.o \
	nrsrc/gammp.o    nrsrc/gcf.o      nrsrc/gser.o  \
        nrsrc/gammln.o   nrsrc/indexx.o   nrsrc/zriddr.o


INCL   = prototypes.h globvars.h nrsrc/nrutil.h nrsrc/nr.h

OPTIMIZE =  -O3  -g

#OPT    +=  -DTESTING

#OPT    +=  -DREDSHIFT_SCALING             # Scale V200, M200, R200, and CC with redshift
#OPT    +=  -DV_SCALING                    # V200 is fixed with redshift, otherwise M200 is

#OPT    +=  -DHDF5                         # enable output to HDF 5 format

OPT    +=  -DMAXWELLIAN                   # use Hernquist 1993 Maxwellian assumption to set bulge/halo velocities
#OPT    +=  -DDF_EDDINGTON                 # use eddington's method to calcualte f_e
#OPT    +=  -DDF_H_MODEL                   # use analytic H90 distribution function

#OPT    +=  -DU4MOD
#OPT    +=  -DSCALEFIX
#OPT    +=  -DTJEOS
#OPT    +=  -DMOLECULAR_COOLING

CFLAGS = $(OPTIMIZE) $(OPT)

LIBS   =  -lm  -lgsl  -lgslcblas

#CC     =  icc
CC	=  cc

.c.o:
	$(CC) $(CFLAGS) -c -o $@ $<

$(EXEC): $(OBJS) 
	$(CC) $(CFLAGS) $(OBJS) -o $(EXEC)  $(LIBS)

$(OBJS): $(INCL) 


clean:
	rm -f $(OBJS) $(EXEC)












