import numpy as np
import math
import gadget_lib.gadget as gadget
import util.utilities as util
import visualization.return_columns_to_sources as vis_col

def checklen(x):
    return len(np.array(x,ndmin=1));

##
## this routine basically functions as a wrapper to get the stellar luminosities
##   and attenuation, to be sent to 'make_threeband_image' for final processing
##  units ::
##      masses in M_sun
##      distances in kpc
##      ages in Gyr
##      metallicity in absolute (solar=0.02)
##      gas: u, rho, in Gadget units (u in km/s^2)
##      luminosity in L_sun
##
def get_attenuated_stellar_luminosities( BAND_IDS, star_pos, gas_pos, bh_pos, \
        stellar_age, stellar_metallicity, stellar_mass, \
        gas_u, gas_rho, gas_hsml, gas_numh, gas_nume, gas_metallicity, gas_mass, \
        bh_luminosity, \
        xrange=0, yrange=0, zrange=0, \
        INCLUDE_BH=0, SKIP_ATTENUATION=0, 
        ADD_BASE_METALLICITY=0., ADD_BASE_AGE=0., 
        IMF_SALPETER=0, IMF_CHABRIER=1, \
        MIN_CELL_SIZE=0.01, OUTER_RANGE_OF_INT=1200., \
        SCATTERED_FRACTION=0.0, \
        REDDENING_SMC=0, REDDENING_LMC=0, REDDENING_MW=0, \
        AGN_MARCONI=0, AGN_HRH=1, AGN_RICHARDS=0, AGN_SDSS=0, \
        HOT_GAS_DTM_FACTOR=0.0 ):

    ## first some basic pre-processing to make sure the numbers are in order
    if ((checklen(gas_pos[0,:])==3) & (checklen(gas_pos[:,0]) !=3)): gas_pos=np.transpose(gas_pos);
    if ((checklen(star_pos[0,:])==3) & (checklen(star_pos[:,0]) !=3)): star_pos=np.transpose(star_pos);
    if (INCLUDE_BH==1): 
        if ((checklen(bh_pos[0,:])==3) & (checklen(bh_pos[:,0]) !=3)): bh_pos=np.transpose(bh_pos);

    stellar_metallicity += ADD_BASE_METALLICITY*0.02;
    gas_metallicity += ADD_BASE_METALLICITY*0.02;
    stellar_age += ADD_BASE_AGE;
    if checklen(stellar_metallicity.shape)>1: stellar_metallicity=stellar_metallicity[:,0];
    if checklen(gas_metallicity.shape)>1: gas_metallicity=gas_metallicity[:,0];
    gas_temp = gadget.gas_temperature(gas_u,gas_nume);
    ## multiply metallicity of >10^6 K gas by HOT_GAS_DTM_FACTOR to account for dust
    ## destruction in hot gas; set the factor to 0 if you want no dust in hot gas or 1
    ## if you want the same dust-to-metals ratio as in the ISM
    gas_metallicity[gas_temp > 1.0e6] *= HOT_GAS_DTM_FACTOR;


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
          vis_col.return_columns_to_sources( source_pos, gas_pos, \
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
        nu_eff_l[i_band] = util.colors_table(np.array([1.0]),np.array([1.0]), \
            BAND_ID=BAND_IDS[i_band],RETURN_NU_EFF=1);
        lums[i_band,:] = stellar_mass * util.colors_table( stellar_age, stellar_metallicity/0.02, \
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
        f_atten = util.attenuate( nu_eff_l[i_band], np.log10(LOS_NH), LOS_Z/0.02, \
          SMC=REDDENING_SMC, LMC=REDDENING_LMC, MW=REDDENING_MW );
        lums_atten[i_band,:] = lums[i_band,:] * \
            ((1.-SCATTERED_FRACTION)*f_atten + SCATTERED_FRACTION);

    return lums, lums_atten;
