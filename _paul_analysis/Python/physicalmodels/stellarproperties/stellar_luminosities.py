import numpy as np
import util
import ctypes
import units.springel_units
import physicalmodels.attenuation.attenuate as atten

def checklen(x):
    return len(np.array(x,ndmin=1));

def vfloat(x):
    return x.ctypes.data_as(ctypes.POINTER(ctypes.c_float));

def fcor(x):
    return np.array(x,dtype='f',ndmin=1)

def ok_scan(input,xmax=1.0e10,pos=0):
    if (pos==1):
        return (np.isnan(input)==False) & (np.fabs(input)<=xmax) & (input > 0.);
    else:
        return (np.isnan(input)==False) & (np.fabs(input)<=xmax);


def get_attenuated_stellar_luminosities( BAND_IDS, star_pos, gas_pos, bh_pos, \
                                        stellar_age, stellar_metallicity, stellar_mass, \
                                        gas_u, gas_rho, gas_hsml, gas_numh, gas_nume, gas_metallicity, gas_mass, \
                                        bh_luminosity, \
                                        xrange=0, yrange=0, zrange=0, \
                                        INCLUDE_BH=0, SKIP_ATTENUATION=0,
                                        IMF_SALPETER=0, IMF_CHABRIER=1, \
                                        MIN_CELL_SIZE=0.01, OUTER_RANGE_OF_INT=1200., \
                                        SCATTERED_FRACTION=0.0, \
                                        REDDENING_SMC=0, REDDENING_LMC=0, REDDENING_MW=0, \
                                        AGN_MARCONI=0, AGN_HRH=1, AGN_RICHARDS=0, AGN_SDSS=0 ):
    
    ## first some basic pre-processing to make sure the numbers are in order
    if ((checklen(gas_pos[0,:])==3) & (checklen(gas_pos[:,0]) !=3)): gas_pos=np.transpose(gas_pos);
    if ((checklen(star_pos[0,:])==3) & (checklen(star_pos[:,0]) !=3)): star_pos=np.transpose(star_pos);
    if (INCLUDE_BH==1):
        if ((checklen(bh_pos[0,:])==3) & (checklen(bh_pos[:,0]) !=3)): bh_pos=np.transpose(bh_pos);
    if checklen(stellar_metallicity.shape)>1: stellar_metallicity=stellar_metallicity[:,0];
    if checklen(gas_metallicity.shape)>1: gas_metallicity=gas_metallicity[:,0];
    
    gas_temp = units.springel_units.gas_code_to_temperature(gas_u,gas_nume);
    gas_metallicity[gas_temp > 1.0e6] = 0.0; ## don't allow hot gas to have dust
    
    
    ## now call the extinction calculation
    Nstar=checklen(star_pos[0,:]);
    if (SKIP_ATTENUATION==0):
        if (INCLUDE_BH==1):
            Nbh=checklen(bh_pos[0,:]);
            source_pos=np.zeros(3,Nstar+Nbh);
            for j in [0,1,2]:
                source_pos[j,0:Nstar]=star_pos[j,:];
                source_pos[j,Nstar:Nstar+Nbh]=bh_pos[j,:];
        else:
            source_pos=star_pos;
    
        LOS_NH, LOS_NH_HOT, LOS_Z = \
            return_columns_to_sources( source_pos, gas_pos, \
                                              gas_u, gas_rho, gas_hsml, gas_numh, gas_nume, gas_metallicity, gas_mass, \
                                              xrange=xrange, yrange=yrange, zrange=zrange, \
                                              MIN_CELL_SIZE=MIN_CELL_SIZE, OUTER_RANGE_OF_INT=OUTER_RANGE_OF_INT, \
                                              TRIM_PARTICLES=1 );

    else: ## SKIP_ATTENUATION==1
        N_sources=checklen(star_pos[0,:]);
        if(INCLUDE_BH==1): N_sources+=checklen(bh_pos[0,:]);
        NHmin=1.0e10; LOS_NH=np.zeros(N_sources)+NHmin; LOS_NH_HOT=np.copy(LOS_NH); LOS_Z=0.*LOS_NH+1.0;
    
    print '<LOS_NH> == ',np.median(LOS_NH),' <LOS_Z> == ',np.median(LOS_Z)


    ## alright now we're ready to get the (intrinsic) stellar luminosities
    nband=checklen(BAND_IDS); lums=np.zeros([nband,Nstar]); nu_eff_l=np.zeros([nband]);
    for i_band in range(nband):
        nu_eff_l[i_band] = colors_table(np.array([1.0]),np.array([1.0]), \
                                            BAND_ID=BAND_IDS[i_band],RETURN_NU_EFF=1);
        lums[i_band,:] = stellar_mass * colors_table( stellar_age, stellar_metallicity/0.02, \
                                            BAND_ID=BAND_IDS[i_band], CHABRIER_IMF=IMF_CHABRIER, SALPETER_IMF=IMF_SALPETER, CRUDE=1, \
                                            UNITS_SOLAR_IN_BAND=1); ## this is such that solar-type colors appear white

    ## if we're using the BH, also get its luminosities at the bands of interest
    if (INCLUDE_BH==1):
        Nbh=checklen(bh_pos[0,:]); Nbands=checklen(BAND_IDS); lums_bh=np.zeros([Nbands,Nbh]);
        for i_bh in range(Nbh):
            lums_bh[:,i_bh] = util.agn_spectrum( nu_eff_l, np.log10(bh_luminosity[i_bh]), \
                                                HRH=AGN_HRH,MARCONI=AGN_MARCONI,RICHARDS=AGN_RICHARDS,SDSS=AGN_SDSS );
        lums_new=np.zeros([Nbands,Nstar+Nbh]);
        for i_band in range(Nbands):
            lums_new[i_band,0:Nstar]=lums[i_band,:];
            lums_new[i_band,Nstar:Nstar+Nbh]=lums_bh[i_band,:];
        lums=lums_new


    ## call the attenuation routine to get the post-extinction luminosities
    lums_atten=1.0*lums;
    LOS_NH_TO_USE = LOS_NH;
    for i_band in range(checklen(BAND_IDS)):
        f_atten = attenuate( nu_eff_l[i_band], np.log10(LOS_NH), LOS_Z/0.02, \
                                 SMC=REDDENING_SMC, LMC=REDDENING_LMC, MW=REDDENING_MW );
        lums_atten[i_band,:] = lums[i_band,:] * \
                                     ((1.-SCATTERED_FRACTION)*f_atten + SCATTERED_FRACTION);
    
    return lums, lums_atten;




##
## return: los_NH_allgas, los_NH_hotphase, los_gas_metallicity
##
def return_columns_to_sources( source_pos, gas_pos, \
                              gas_u, gas_rho, gas_hsml, gas_numh, gas_nume, gas_metallicity, gas_mass, \
                              xrange=0, yrange=0, zrange=0, \
                              MIN_CELL_SIZE=0.01, OUTER_RANGE_OF_INT=1200., \
                              TRIM_PARTICLES=1 ):
    
    ## check the ordering of the position matrices:
    if ((checklen(gas_pos[0,:])==3) & (checklen(gas_pos[:,0]) !=3)): gas_pos=np.transpose(gas_pos);
    if ((checklen(source_pos[0,:])==3) & (checklen(source_pos[:,0]) !=3)): source_pos=np.transpose(source_pos);
    ## and that metallicities are a vector, not a matrix
    if (len(gas_metallicity.shape)>1): gas_metallicity=gas_metallicity[:,0]
    
    if ((checklen(gas_pos[:,0]) != 3) | (checklen(gas_pos[0,:]) <= 1)):
        print 'ERROR WILL OCCUR :: need pos to be (3,N)'
    
    x=source_pos[0,:] ; y=source_pos[1,:] ; z=source_pos[2,:]
    if(checklen(xrange)<=1): xrange=[np.min(x),np.max(x)];
    if(checklen(yrange)<=1): yrange=[np.min(y),np.max(y)];
    xr=xrange; yr=yrange;
    if(checklen(zrange)<=1):
        zrr=np.sqrt((xr[1]-xr[0])**2.+(yr[1]-yr[0])**2.)/np.sqrt(2.);
        zmin=np.median(z)-zrr; zmax=np.median(z)+zrr;
        if (np.min(z) > zmin): zmin=np.min(z);
        zrange=[zmin,zmax]; print 'z_range (calc) == ',zrange
    zr=zrange;
    x00=0.5*(xr[1]+xr[0]); y00=0.5*(yr[1]+yr[0]); z00=0.5*(zr[1]+zr[0]);
    tolfac = 1.0e10;
    if (TRIM_PARTICLES==1):
        tolfac = 0.05;
    #tolfac = -0.01;
    ## trim down the incoming list to only whats in the range plotted
    ##   (saves a ton of time and memory overflow crashes)
    
    dx=(0.5+tolfac)*(xr[1]-xr[0]); dy=(0.5+tolfac)*(yr[1]-yr[0]); dz=(0.5+tolfac)*(zr[1]-zr[0]);
    ok_sources=ok_scan(x-x00,xmax=dx) & ok_scan(y-y00,xmax=dy) & ok_scan(z-z00,xmax=dz);
    x=gas_pos[0,:] ; y=gas_pos[1,:] ; z=gas_pos[2,:]
    gw=gas_rho ; gh=gas_hsml ; gz=gas_metallicity ; gm=gas_mass
    ok_gas=ok_scan(x-x00,xmax=dx) & ok_scan(y-y00,xmax=dy) & ok_scan(z-z00,xmax=dz) & \
        ok_scan(gw,pos=1) & ok_scan(gh,pos=1) & ok_scan(gz,pos=1) & ok_scan(gm,pos=1,xmax=1.0e40);
    
    Ngas = checklen(gas_mass[ok_gas]);
    Nstars = checklen(source_pos[0,ok_sources]);
    if (Nstars<=1) or (Ngas<=1):
        print ' UH-OH: EXPECT ERROR NOW, there are no valid source/gas particles to send!'
        print 'Ngas=',Ngas,'Nstars=',Nstars,'dx=',dx,'dy=',dy,'dz=',dz,'x00=',x00,'y00=',y00,'z00=',z00
        return -1,-1,-1;
    
    dzmax=np.max(gas_pos[2,ok_gas])-z00;
    if(dzmax<OUTER_RANGE_OF_INT): OUTER_RANGE_OF_INT=dzmax;
    print 'PASSING: N_gas=',Ngas,'N_sources=',Nstars,'MaxDist=',OUTER_RANGE_OF_INT,'MinCell=',MIN_CELL_SIZE;
    Nbh=0; theta=1.0e-4; phi=1.0e-4;

    ## load the routine we need
    exec_call=util.dir.c_routines_dir()+'/LOS_column_singlePOV/getnh.so'
    NH_routine=ctypes.cdll[exec_call];



    ## cast the variables to store the results
    nh_out_cast=ctypes.c_float*Nstars;
    los_NH_out=nh_out_cast(); los_NH_hot_out=nh_out_cast(); los_Z_out=nh_out_cast();

    ## ok this is a bit arcane but the routine will read appropriately this block order
    Coord = np.zeros((Ngas+Nstars,10),dtype='f');
    Coord[0:Ngas,0] = gas_pos[0,ok_gas]-x00;
    Coord[0:Ngas,1] = gas_pos[1,ok_gas]-y00;
    Coord[0:Ngas,2] = gas_pos[2,ok_gas]-z00;
    Coord[0:Ngas,3] = gas_u[ok_gas]
    Coord[0:Ngas,4] = gas_rho[ok_gas]
    Coord[0:Ngas,5] = gas_hsml[ok_gas]
    Coord[0:Ngas,6] = gas_numh[ok_gas]
    Coord[0:Ngas,7] = gas_nume[ok_gas]
    Coord[0:Ngas,8] = gas_metallicity[ok_gas]
    Coord[0:Ngas,9] = gas_mass[ok_gas]
    Coord[Ngas:Nstars+Ngas,0] = source_pos[0,ok_sources]-x00;
    Coord[Ngas:Nstars+Ngas,1] = source_pos[1,ok_sources]-y00;
    Coord[Ngas:Nstars+Ngas,2] = source_pos[2,ok_sources]-z00;
    Coord=np.copy(np.transpose(Coord));

    ## main call to the NH-calculation routine
    NH_routine.getnh(   ctypes.c_int(Ngas),
                        ctypes.c_int(Nstars),
                        ctypes.c_int(Nbh),
                        ctypes.c_float(theta),
                        ctypes.c_float(phi),
                        vfloat(Coord),
                        ctypes.byref(los_NH_out),
                        ctypes.byref(los_NH_hot_out),
                        ctypes.byref(los_Z_out),
                        ctypes.c_float(OUTER_RANGE_OF_INT),
                        ctypes.c_float(MIN_CELL_SIZE) );
    ## now put the output arrays into a useful format
    print type(los_NH_out), los_NH_out
    los_NH = np.ctypeslib.as_array(los_NH_out);     # removed a np.copy() as below
    los_NH_hot = np.ctypeslib.as_array(np.copy(los_NH_hot_out));
    los_Z = np.ctypeslib.as_array(np.copy(los_Z_out));
     
    # trap for really low NH value and zero metallicity (make it small instead)
    low_NH = 1.0e10;
    los_NH[los_NH<low_NH]=low_NH; los_NH_hot[los_NH_hot<low_NH]=low_NH;
    los_Z[los_Z<=1.0e-5]=1.0e-5;
     
    ## assign strong attenuation to all 'off-grid' sources, then fill in calc. vals
    Nstarstot=checklen(source_pos[0,:]);
    los_NH_allgas=np.zeros(Nstarstot,dtype='f')+1.0e23;
    los_NH_hotgas=np.zeros(Nstarstot,dtype='f')+1.0e23;
    los_gas_metallicity=np.zeros(Nstarstot,dtype='f')+0.02;
    nok=checklen(los_NH_allgas[ok_sources])
    los_NH_allgas[ok_sources]=fcor(los_NH[0:Nstars]);
    los_NH_hotgas[ok_sources]=fcor(los_NH_hot[0:Nstars]);
    los_gas_metallicity[ok_sources]=fcor(los_Z[0:Nstars]);

    return los_NH_allgas, los_NH_hotgas, los_gas_metallicity;



## routines from colors_sps module
def colors_table( age_in_Gyr, metallicity_in_solar_units,
                 BAND_ID=0, SALPETER_IMF=0, CHABRIER_IMF=1, QUIET=0, CRUDE=0,
                 RETURN_NU_EFF=0, RETURN_LAMBDA_EFF=0, UNITS_SOLAR_IN_BAND=0 ):
    return colors_table( age_in_Gyr, metallicity_in_solar_units,
                             BAND_ID=BAND_ID, SALPETER_IMF=SALPETER_IMF, CHABRIER_IMF=CHABRIER_IMF, QUIET=QUIET, CRUDE=CRUDE,
                             RETURN_NU_EFF=RETURN_NU_EFF, RETURN_LAMBDA_EFF=RETURN_LAMBDA_EFF, UNITS_SOLAR_IN_BAND=UNITS_SOLAR_IN_BAND )


def colors_table( age_in_Gyr, metallicity_in_solar_units,
                 BAND_ID=0, SALPETER_IMF=0, CHABRIER_IMF=1, QUIET=0, CRUDE=0,
                 RETURN_NU_EFF=0, RETURN_LAMBDA_EFF=0, UNITS_SOLAR_IN_BAND=0 ):
    
    #import utilities as util
    import numpy as np
    import scipy.ndimage.interpolation as interpolate
    import struct
    
    age_in_Gyr=np.array(age_in_Gyr,ndmin=1);
    metallicity_in_solar_units=np.array(metallicity_in_solar_units,ndmin=1);
    
    band=BAND_ID; # default=bolometric
    j = [  0,  6,  7,  8,  9, 10, 11, 12, 13,  1,   2,   3,   4,   5] # ordering I'm used to
    i = [  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10,  11,  12,  13] # ordering of this
    band_standardordering = band
    band = j[band]
    if (band > 13):
        print 'BAND_ID must be < 13';
        return 0;
    
    b=['Bolometric', \
       'Sloan u','Sloan g','Sloan r','Sloan i','Sloan z', \
       'Johnsons U','Johnsons B', 'Johnsons V','Johnsons R','Johnsons I', \
       'Cousins J','Cousins H','Cousins K']
    if (QUIET==0): print 'Calculating M/L in '+str(b[band])+' ('+str(band)+','+str(band_standardordering)+')'
    
    if (RETURN_NU_EFF==1) or (RETURN_LAMBDA_EFF==1):
        lam_eff=np.array([1.e-5, 3541., 4653., 6147., 7461., 8904., 3600., 4400., \
                          5556., 6940., 8700., 12150., 16540., 21790.]);
        nu_eff = 2.998e18 / lam_eff;
        if (RETURN_NU_EFF==1): return nu_eff[band];
        if (RETURN_LAMBDA_EFF==1): return lam_eff[band];

#froot = util.return_python_routines_homedir()+'/colors_sps/'; # directory in which the data binaries are stored
    if (CHABRIER_IMF==1): fname=util.dir.sps_dir()+'colors.chabrier.dat'
    if (SALPETER_IMF==1): fname=util.dir.sps_dir()+'colors.salpeter.dat'

    lut = open(fname,'r');
    lut_dat = lut.read();
    Nl,Na,Nz = struct.unpack('3i',lut_dat[0:12])
    z_grid = np.array(struct.unpack(str(Nz)+'d',lut_dat[12:12+8*Nz]))
    age_grid = np.array(struct.unpack(str(Na)+'d',lut_dat[12+8*Nz:12+8*Nz+8*Na]))
    l_all_l = np.array(struct.unpack(str(Nl*Na*Nz)+'d',lut_dat[12+8*Nz+8*Na:12+8*Nz+8*Na+8*Nl*Na*Nz]))
    l_all = np.transpose(l_all_l.reshape(Nz,Na,Nl))
    lut.close()
    
    l_band = np.zeros((Na,Nz),dtype=np.float64);
    for iz in range(Nz): l_band[:,iz]=l_all[band,:,iz]
    
    # allow for extreme metallicities (extrapolate linearly past table)
    push_metals = 1;
    if (push_metals==1):
        Nz = Nz + 1;
        z_ext = [1000.0];
        z_grid = np.concatenate([z_grid,z_ext])
        lb1 = l_band[:,Nz-3]
        lb2 = l_band[:,Nz-2]
        lbx = np.zeros((Na,Nz),dtype=np.float64)
        lbx[:,0:Nz-1] = l_band
        lbx[:,Nz-1] = (lb2 - lb1) / (np.log10(z_grid[Nz-2]/z_grid[Nz-3])) * \
            np.log10(z_grid[Nz-1]/z_grid[Nz-2])
        l_band = lbx;

    # get the x-axis (age) locations of input points
    ia_pts=np.interp(np.log10(age_in_Gyr)+9.0,age_grid,np.arange(0,Na,1));
    # this returns the boundary values for points outside of them (no extrapolation)
    #f=interp.interp1d(age_grid,np.arange(0,Na,1),kind='linear');
    #ia_pts=f(np.log10(age_in_Gyr)+9.0);
    
    # get the y-axis (metallicity) locations of input points
    zsun = 0.02;
    iz_pts=np.interp(np.log10(metallicity_in_solar_units*zsun),np.log10(z_grid),np.arange(0,Nz,1));
    #f=interp.interp1d(np.log10(z_grid),np.arange(0,Nz,1),kind='linear');
    #iz_pts=f(np.log10(metallicity_in_solar_units*zsun));
    
    if (CRUDE==1):
        ia_pts=np.around(ia_pts).astype(int);
        iz_pts=np.around(iz_pts).astype(int);
        print ia_pts, iz_pts, ia_pts, iz_pts
        print np.min( ia_pts), np.min( iz_pts), np.min( ia_pts), np.min( iz_pts)
        print np.max( ia_pts), np.max( iz_pts), np.max( ia_pts), np.max( iz_pts)
        ia_pts[ia_pts < 0] = np.max(ia_pts) 
        iz_pts[iz_pts < 0] = np.max(iz_pts)
        l_b=l_band[ia_pts,iz_pts];
    else:
        l_b = interpolate.map_coordinates(l_band, (ia_pts,iz_pts), order=1);
    l_b = 10.**l_b
    
    # output is currently L/M in L_sun_IN_THE_BAND_OF_INTEREST/M_sun,
    # but we want our default to be L/M in units of L_bolometric/M_sun = 3.9e33/2.0e33, so
    #   need to get rid fo the L_sun_IN_THE_BAND_OF_INTEREST/L_bolometric
    
    # AB system solar luminosities used for determining L_sun in absolute units for each of these
    N_BANDS=14
    mag_sun_ab = np.zeros(N_BANDS,dtype=float)
    mag_sun_ab[0] = 4.74;
    l_bol_sun = 3.9e33; # bolometric solar in erg/s
    mag_sun_ab[1] = 6.34;  #U (BESSEL)
    mag_sun_ab[2] = 5.33;  #B (BESSEL)
    mag_sun_ab[3] = 4.81;  #V (BESSEL)
    mag_sun_ab[4] = 4.65;  #R (KPNO)
    mag_sun_ab[5] = 4.55;  #I (KPNO)
    mag_sun_ab[6] = 4.57;  #J (BESSEL)
    mag_sun_ab[7] = 4.71;  #H (BESSEL)
    mag_sun_ab[8] = 5.19;  #K (BESSEL)
    mag_sun_ab[9] = 6.75;  #SDSS u (unprimed AB)
    mag_sun_ab[10] = 5.33; #SDSS g (unprimed AB)
    mag_sun_ab[11] = 4.67; #SDSS r (unprimed AB)
    mag_sun_ab[12] = 4.48; #SDSS i (unprimed AB)
    mag_sun_ab[13] = 4.42; #SDSS z (unprimed AB)
    
    # Effective wavelengths of the bands (in Angstroms), to compute nuLnu<->Lnu
    # UBVRIJHK from http://cassfos02.ucsd.edu/physics/ph162/mags.html
    # SDSS ugriz from http://www.sdss.org/dr4/instruments/imager/index.html#filters
    lambda_eff = np.zeros(N_BANDS,dtype=float);
    lambda_eff[0] = 4243.93;  #bolometric, no nu
    lambda_eff[1] = 3600.0;  #U
    lambda_eff[2] = 4400.0;  #B
    lambda_eff[3] = 5556.0;  #V
    lambda_eff[4] = 6940.0;  #R
    lambda_eff[5] = 8700.0;  #I
    lambda_eff[6] = 12150.;  #J
    lambda_eff[7] = 16540.;  #H
    lambda_eff[8] = 21790.;  #K
    lambda_eff[9]  = 3551.;  #SDSS u
    lambda_eff[10] = 4686.;  #SDSS g
    lambda_eff[11] = 6165.;  #SDSS r
    lambda_eff[12] = 7481.;  #SDSS i
    lambda_eff[13] = 8931.;  #SDSS z
    c_light = 2.998e10; # speed of light in cm/s
    nu_eff  = c_light / (lambda_eff * 1.0e-8); # converts to nu_eff in Hz
    
    ten_pc   = 10.e0 * 3.086e18; # 10 pc in cm
    log_S_nu = -(mag_sun_ab + 48.6)/2.5; # zero point definition for ab magnitudes
    S_nu     = 10.**log_S_nu; # get the S_nu at 10 pc which defines M_AB
    lnu_sun_band = S_nu * (4.*3.14159*ten_pc*ten_pc); # multiply by distance modulus
    nulnu_sun_band = lnu_sun_band * nu_eff; # multiply by nu_eff to get nu*L_nu
    l_bol_sun = nulnu_sun_band[0];
    
    if (UNITS_SOLAR_IN_BAND==0):
        l_b *= nulnu_sun_band[band_standardordering] / l_bol_sun; 
    
    return l_b;




## routines from attenuation module
def attenuate( nu_in_Hz, log_NH, metallicity_in_solar, \
              SMC=0, LMC=0, MW=0, BB=0, IR=0, SX=0, HX=0):
    return atten.attenuate( nu_in_Hz, log_NH, metallicity_in_solar, \
                           SMC=SMC, LMC=LMC, MW=MW, BB=BB, IR=IR, SX=SX, HX=HX)
