## returns the attenuation at a given frequency nu (in Hz) for a 
##    given column density (log(NH/cm^-2)), and metallicity (in solar units)
##
##		i.e. for some luminosity L, observed = L * attenuate(nu,alog10(NH),Z/0.02)
##
##	flags: can change the reddening curve adopted 
##		SMC (default) = smc-like reddening curve (strong UV & mid-IR extinction)
##			observations find this most appropriate for quasars: Hopkins et al. 2004
##		MW  = MW-like (weaker UV & mid-IR; strong 2200 angstrom bump) 
##		LMC = LMC-like (basically halfway between the two)
##			(gas to dust ratio in each case is MW-like average scaled 
##				linearly with metallicity)
##
##	can call for specific bands (i.e. net attenuation integrated over the band):
##		BB = B-band
##      IR = mid-IR (15 microns), 
##		SX = soft X-ray (0.5-2 keV), 
##		HX = hard X-ray (2-10 keV)
##
##  the dust extinction curves are from Pei et al. 1992.; gas-to-dust ratios 
##	  from Bouchet et al. 1985. 
##
##  above the Lyman edge, photoionization cross sections computed following 
##   Morrison & McCammon 1983. the code follows Brant's decomposition of this into 
##   a metals-free component and a metals component that scales linearly with 
##   the input metallicity
##
##  compton scattering is treated as non-relativistic achromatic Thompson scattering
##	 at low energies (<~ 4 keV), and above this the transmission curves from 
##   Matt, Pompilio, & La Franca 1995 have been extracted to interpolate over a wide 
##   range of frequencies and NH. These are decomposed into a H-dominated continuum 
##   component that basically scales in a metal-free manner, and an iron flourescence 
##   component (between ~5.8-7.4 keV) that emits in a manner scaling with metallicity. 
##   the iron flourescence is fine for most applications, but can give strange 
##   (although fixed s.t. non-divergent and still monotonic) behavior when you're  
##   simultaneously at very low but non-zero metallicities (<~ 0.1 solar) 
##   and very high (logNH > 25-25.5) column densities -- in this regime the 
##   flourescence should really be calculated self-consistently. 
##
##

import numpy as np
import ctypes
import util

#def vdouble(x):
#    return x.ctypes.data_as(ctypes.POINTER(ctypes.c_double));


def opacity_per_solar_metallicity( f0 ):
    return opacity_per_solar_metallicity( f0 );

def attenuate( nu_in_Hz, log_NH, metallicity_in_solar, \
	SMC=0, LMC=0, MW=0, BB=0, IR=0, SX=0, HX=0):

    ## location of shared library
    exec_call=util.dir.atten_dir()+'/attenuate_py.so'
    lib=ctypes.cdll[exec_call];

    ## default to SMC-like reddening
    dust_key = 2
    if (MW==1):  dust_key = 0
    if (LMC==1): dust_key = 1
    if (SMC==1): dust_key = 2

    ## frequency
    if (BB==1): nu_in_Hz = -1.0
    if (IR==1): nu_in_Hz = -2.0
    if (SX==1): nu_in_Hz = -3.0
    if (HX==1): nu_in_Hz = -4.0	

    NH = 10.**np.array(log_NH,ndmin=1,dtype='d')
    nu_in_Hz=np.array(nu_in_Hz,ndmin=1,dtype='d');
    metallicity_in_solar=np.array(metallicity_in_solar,ndmin=1,dtype='d');
    N_nu = len(nu_in_Hz)
    N_NH = len(NH)
    N_metal = len(metallicity_in_solar)
    if (N_metal <= 1): metallicity_in_solar = 0.*NH + metallicity_in_solar[0]
    atten = np.zeros(N_nu*N_NH,dtype='d')
    NH = np.array(NH,ndmin=1,dtype='d')
    metallicity_in_solar = np.array(metallicity_in_solar,ndmin=1,dtype='d')
    out_cast=ctypes.c_double*(N_nu*N_NH); atten=out_cast();
    
    ## actual function call
    lib.main( ctypes.c_int(N_nu), vdouble(nu_in_Hz), \
        ctypes.c_int(N_NH), vdouble(NH), \
        vdouble(metallicity_in_solar), \
        ctypes.c_int(dust_key), \
        ctypes.byref(atten) );

    ## now put the output arrays into a useful format 
    atten = np.copy(np.ctypeslib.as_array(atten));

    atten_f = np.zeros((N_nu,N_NH),dtype='d')
    for i in range(N_nu):
        atten_f[i,:] = atten[i*N_NH:(i+1)*N_NH]
    atten=atten_f
    atten[atten==0]=1.0e-40; 
    atten[np.isnan(atten)]=1.0e-40;

    return atten



#-------------------------------------------------------------------------------------
# Pei (1992) MW,LMC,SMC dust laws used to calculate differential extinction
#   for dust at most frequencies
# (INPUT LAMBDA IN MICRONS)
#     NOTE: Pei's model is based on observations from lambda > 1000 Angstroms
#       (f < 3.0e15 Hz), but the grain models extrapolate these fits at least to
#		lambda = 0.001 micrometer (= 10 Angstroms, or f = 3.0e17 Hz), where
#		the extinction drops off rapidly. So it's probably safe to use for all freq.
#-------------------------------------------------------------------------------------
def pei_dustparam( lambda_microns, MW=0, LMC=0, SMC=1 ):
    
    if (MW==1):
        a = [165., 14., 0.045, 0.002, 0.002, 0.012]
        l = [0.047, 0.08, 0.22, 9.7, 18., 25.]
        b = [90., 4.00, -1.95, -1.95, -1.80, 0.00]
        n = [2.0, 6.5, 2.0, 2.0, 2.0, 2.0]
        R_V = 3.08
    if (LMC==1):
        a = [175., 19., 0.023, 0.005, 0.006, 0.020]
        l = [0.046, 0.08, 0.22, 9.7, 18., 25.]
        b = [90., 5.50, -1.95, -1.95, -1.80, 0.00]
        n = [2.0, 4.5, 2.0, 2.0, 2.0, 2.0]
        R_V = 3.16
    if (SMC==1):
        a = [185., 27., 0.005, 0.010, 0.012, 0.030]
        l = [0.042, 0.08, 0.22, 9.7, 18., 25.]
        b = [90., 5.50, -1.95, -1.95, -1.80, 0.00]
        n = [2.0, 4.0, 2.0, 2.0, 2.0, 2.0]
        R_V = 2.93
    
    xsi = 0.*lambda_microns;
    for i in range(len(a)):
        xsi += a[i] / ( (lambda_microns/l[i])**(n[i]) + (l[i]/lambda_microns)**(n[i]) + b[i] )
    R_lam = (1.0 + R_V) * xsi;

    return xsi;


#-------------------------------------------------------------------------------------
# For 0.03 keV < E < 10 keV
#   (7.2e15 < nu[Hz] < 2.4e18  or   1.2 < lambda[Angstroms] < 413)
#   we use the photoelectric absorption cross sections of
#   Morrison & McCammon (1983)
#     NOTE: these assume solar abundances and no ionization,
#             the appropriate number probably scales linearly with both
#				- keyword ABUNDANCE gives ratio to solar/MW metallicity
#   (this is all for the COMPTON THIN regime)
#-------------------------------------------------------------------------------------
def morrison_photoelec( frequency ):
    # (can restructure this to allow vector i/o if the time for the run is too long)
    f_003keV = 7.253e15
    f_10keV  = 2.418e18
    keV_per_Hz = 4.13608e-18
    fkeV = frequency * keV_per_Hz	# convert frequency from Hz to keV
    
    # Now set the appropriate polynomial terms from Table 2
    #   of Morrison & McCammon (for a given frequency range)
    if  (fkeV < 0.100):
        c = [	17.3,	608.1,	-2150.0	]
    elif(fkeV < 0.284):
        c = [	34.6,	267.9,	 -476.1	]
    elif(fkeV < 0.400):
        c = [	78.1,	 18.8,	    4.3	]
    elif(fkeV < 0.532):
        c = [	71.4,	 66.8,	  -51.4	]
    elif(fkeV < 0.707):
        c = [	95.5,	145.8,	  -61.1	]
    elif(fkeV < 0.867):
        c = [   308.9,  -380.6,	  294.0	]
    elif(fkeV < 1.303):
        c = [   120.6,	169.3,	  -47.7	]
    elif(fkeV < 1.840):
        c = [   141.3,	146.8,	  -31.5	]
    elif(fkeV < 2.471):
        c = [   202.7,	104.7,	  -17.0	]
    elif(fkeV < 3.210):
        c = [   342.7,	 18.7,	    0.0	]
    elif(fkeV < 4.038):
        c = [   352.2,	 18.7,	    0.0	]
    elif(fkeV < 7.111):
        c = [   433.9,	 -2.4,	  -0.75	]
    elif(fkeV < 8.331):
        c = [   629.0,	 30.9,      0.0	]
    elif(fkeV < 10.00):
        c = [   701.2,	 25.2,	    0.0	]
    else:
        c = [   701.2,	 25.2,	    0.0	]

    # Use these coefficients to calculate the cross section (in 10^-24 cm^2) per hydrogen atom
    return (c[0] + c[1]*fkeV + c[2]*(fkeV*fkeV)) / (fkeV*fkeV*fkeV)  * 1.0e-24;



#--------------------------------------------------------------------------
# Function to calculate the cross section per H atom at a given
#    frequency f (in Hz)
#--------------------------------------------------------------------------
def cross_section( f0, METALLICITY_OVER_SOLAR=1.0 ):
    f0=np.array(f0);
    SIGMA = 0.0*np.array(f0);
    
    # For 0.03 keV < E < 10 keV
    #   (7.2e15 < nu[Hz] < 2.4e18  or   1.2 < lambda[Angstroms] < 413)
    #   we use the photoelectric absorption cross sections of
    #   Morrison & McCammon (1983)
    #     NOTE: these assume solar abundances and no ionization,
    #             the appropriate number probably scales linearly with both
    #   (this is all for the COMPTON THIN regime)
    #
    f_003keV = 7.253e15
    f_H_edge = 1.362/3.0 * f_003keV
    f_10keV  = 2.418e18
    if (SIGMA.size > 1):
        for i in range(SIGMA.size):
            if (f_H_edge <= f0[i]): SIGMA[i] += morrison_photoelec(f0[i]);
    else:
        if (f_H_edge <= f0): SIGMA += morrison_photoelec(f0);
        SIGMA *= METALLICITY_OVER_SOLAR;
    # (make sure this uses the total column, not just the
    #   neutral column, as ions still contribute), otherwise divide by mean (~1/3,say)
    
    # For optical-IR regions, we use the Pei numerical approximations below.
    #
    # xsi = tau(lambda)/tau(B) is the ratio of extinction at lambda to the
    #    extinction in the B-band.
    # k = 10^21 (tau_B / NH)   (NH in cm^2) gives the dimensionless gas-to-dust
    #    ratio, with k=0.78 for MW, k=0.16 for LMC, k=0.08 for SMC.
    #    k is INDEPENDENT of the grain properties, and seems to scale rougly
    #    linearly with metallicity
    # so, for now, assume solar metallicity and k = k_MW = 0.78. we can rescale later.
    #
    # tau_B = ( NH / (10^21 cm^-2) ) * k --> SIGMA_B = k*10^-21  cm^2
    # tau_lambda = xsi * tau_B --> SIGMA = xsi * SIGMA_B
    #
    # k = 0.78 for the MW
    # k = 0.08 for the SMC, approximately in line with the MW/LMC/SMC metallicity
    #  sequence, so we take a k_MW then scaled by the metallicity
    #
    k = 0.78 * METALLICITY_OVER_SOLAR
    if (SIGMA.size > 1):
        ok=(f0 < f_003keV)
        lambda_microns = 2.998e14 / f0[ok];   # convert frequency [Hz] to wavelength [microns]
        xsi = pei_dustparam( lambda_microns, SMC=1 );
        SIGMA[ok] += xsi * k * 1.0e-21;
    else:
        if (f0 < f_003keV):
            lambda_microns = 2.998e14 / f0;   # convert frequency [Hz] to wavelength [microns]
            xsi = pei_dustparam( lambda_microns, SMC=1 );
            SIGMA += xsi * k * 1.0e-21;

# No double-counting, but I have checked it in detail, and there is really almost
#    absolutely no difference -- even up to NH~10^24-25, it's like a factor of 1.1
#    or so between including both these factors and not, and only in a very, very
#    narrow frequency range

# + Thompson scattering cross sections (NR Compton scattering)
#    full compton scattering dsigma/dOmega = (1/2)*r0^2*(1+cos^2(theta))
#    ( sigma_thompson = (8pi/3)*r0^2 )

    sigma_thompson = 6.65e-25	## cm^-2
    SIGMA += sigma_thompson * 0.5
    # b/c half scattered into viewing angle, half out (*very* roughly) (reflection)
    
    return SIGMA


## simple routine for opacity scaled to solar metallicity at frequencies f0
def opacity_per_solar_metallicity( f0 ):
    mu = 0.75 ## assume dense gas
    return cross_section( f0, METALLICITY_OVER_SOLAR=1.0 ) / ( mu * 1.67e-24 );
