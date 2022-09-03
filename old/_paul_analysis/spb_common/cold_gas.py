# -*- coding: utf-8 -*-
"""Module for finding the neutral hydrogen in a halo. Each class has a function get_reproc_rhoHI, which returns
the neutral hydrogen density in *physical* atoms / cm^3
    Contains:
        StarFormation - Partially implements the star formation model of Springel & Hernquist 2003.
        YajimaRT - implements a fit to the rt formula of Yajima 2011.
        RahmatiRT - implements a fit to the rt formula of Rahmati 2012.
                  get_code_rhoHI returns the rhoHI density as given by Arepo
    Method:
        get_reproc_HI - Gets a corrected neutral hydrogen density
"""

import numpy as np
import numexpr as nmex

class StarFormation:
    """Calculates the fraction of gas in cold clouds, following
    Springel & Hernquist 2003 (astro-ph/0206393) and
    Nagamine, Springel and Hernquist 2004 (astro-ph/0305409).

    Parameters (of the star formation model):
        hubble - hubble parameter in units of 100 km/s/Mpc
        t_0_star - star formation timescale at threshold density
             - (MaxSfrTimescale) 1.5 in internal time units ( 1 itu ~ 0.97 Gyr/h)
        beta - fraction of massive stars which form supernovae (FactorSN) 0.1 in SH03.
        T_SN - Temperature of the supernova in K- 10^8 K SH03. (TempSupernova) Used to calculate u_SN
        T_c  - Temperature of the cold clouds in K- 10^3 K SH03. (TempClouds) Used to calculate u_c.
        A_0  - Supernova evaporation parameter (FactorEVP = 1000).
        WARNING: cooling time for the cloud is hard-coded, as is rescaled beta.
        This uses values from the default GFM parameters. Do not try to change them!
    """
    def __init__(self,hubble=0.7,t_0_star=2.27,beta=0.1,T_SN=5.73e7,T_c = 1000, A_0=573):
        #Some constants and unit systems
        #Internal gadget mass unit: 1e10 M_sun/h in g/h
        self.UnitMass_in_g=1.989e43
        #Internal gadget length unit: 1 kpc/h in cm/h
        self.UnitLength_in_cm=3.085678e21
        #Internal velocity unit : 1 km/s in cm/s
        self.UnitVelocity_in_cm_per_s=1e5
        self.UnitTime_in_s = (self.UnitLength_in_cm/self.UnitVelocity_in_cm_per_s)
        #proton mass in g
        self.protonmass=1.67262178e-24
        self.hy_mass = 0.76 # Hydrogen massfrac
        self.gamma=5./3.
        #Boltzmann constant (cgs)
        self.boltzmann=1.38066e-16

        #Gravitational constant (cgs)
        #self.gravity = 6.672e-8

        #100 km/s/Mpc in s
        #self.h100 = 3.2407789e-18

        self.hubble=hubble

        #Supernova timescale in s
        self.t_0_star=t_0_star*self.UnitTime_in_s/self.hubble # Now in s

        #This is the value output by an IMF rescaling, if GFM_EVOLUTION is on.
        self.beta = 0.264089

        self.A_0 = A_0

        #u_c - thermal energy in the cold gas.
        meanweight = 4 / (1 + 3 * self.hy_mass)          #Assuming neutral gas for u_c
        self.u_c =  1. / meanweight * (1.0 / (self.gamma-1)) * (self.boltzmann / self.protonmass) *T_c

        #SN energy: u_SN = (1-beta)/beta epsilon_SN
        meanweight = 4 / (8 - 5 * (1 - self.hy_mass))    #Assuming FULL ionization for u_H
        self.u_SN =  1. / meanweight * (1.0 / (self.gamma -1)) * (self.boltzmann / self.protonmass) * T_SN

        self.UnitDensity_in_cgs = self.UnitMass_in_g/self.UnitLength_in_cm**3
        #This is the hard-coded "very high density" at which the cooling rate is computed
        #in the Gadget 2-phase SFR model.
        self.rhoinf = 277.476 * self.UnitDensity_in_cgs*self.hubble**2
        #The cooling time for cold clouds, as computed inside the Gadget star-forming model.
        self.tcool = 4.64419e-10 * self.UnitTime_in_s/self.hubble

    def cold_gas_frac(self,rho, tcool,rho_thresh):
        """Calculates the fraction of gas in cold clouds, following
        Springel & Hernquist 2003 (astro-ph/0206393) and
        Nagamine, Springel and Hernquist 2004 (astro-ph/0305409).

        Parameters:
            rho - Density of hot gas (hydrogen /cm^3)
            tcool - cooling time of the gas. Zero if gas is being heated (s)
            rho_thresh - SFR threshold density in hydrogen /cm^3
        Returns:
            The fraction of gas in cold clouds. In practice this is often 1.
        """
        #Star formation timescale
        t_star = self.t_0_star*np.sqrt(rho_thresh/rho)

        #Supernova evaporation parameter
        Arho = self.A_0 * (rho_thresh/rho)**0.8

        #Internal energy in the hot phase
        u_h = self.u_SN/ (1+Arho)+self.u_c

        # a parameter: y = t_star \Lambda_net(\rho_h,u_h) / \rho_h (\beta u_SN - (1-\beta) u_c) (SH03)
        # Or in Gadget: y = t_star /t_cool * u_SN / ( \beta u_SN - (1-\beta) u_c)
        y = u_h / (self.beta*self.u_SN - (1-self.beta)*self.u_c)
        times = tcool / t_star
        #The cold gas fraction
        f_c = 1.+ times/(2*y) - np.sqrt(times/y+times**2/(4*y**2))
        return f_c

    def get_rho_thresh(self):
        """
        This function calculates the physical density threshold for star formation.
        It can be specified in two ways: either as a physical density threshold
        (rho_phys_thresh ) in units of hydrogen atoms per cm^3
        or as a critical density threshold (rho_crit_thresh) which is in units of rho_baryon at z=0.
        DO NOT USE THIS FUNCTION as we cannot calculate the cooling rates here
        Returns:
                rho_thresh in units of H atoms/cm^3
        """
        u_h = self.u_SN / self.A_0

        #u_4 - thermal energy at 10^4K
        meanweight = 4 / (8 - 5 * (1 - self.hy_mass))    #Assuming FULL ionization for u_H
        u_4 =  1. / meanweight * (1.0 / (self.gamma-1)) * (self.boltzmann / self.protonmass) *1e4
        #Note: get_asymptotic_cool does not give the full answer, so do not use it.
        coolrate = u_h / self.tcool / self.rhoinf

        x = (u_h - u_4) / (u_h - self.u_c)
        physdens =  x / (1 - x)**2 * (self.beta * self.u_SN - (1 -self.beta) * self.u_c) /(self.t_0_star * coolrate)
        return physdens / self.protonmass *self.hy_mass

    def get_asmyptotic_cool(self,u_h):
        """
        Get the cooling time for the asymptotically dense limit of cooling,
        where ne = 1, at a temperature u_h.
        Neglect all cooling except free-free; Gadget includes Compton from the CMB,
        but this will be negligible for high temperatures.

        Assumes no heating.
        This function is not useful as it does not really do the right calculation,
        but I leave it here in case it is one day.
        """
        yhelium = (1 - self.hy_mass) / (4 * self.hy_mass)
        mu = 4 / (8 - 5 * (1 - self.hy_mass))    #Assuming FULL ionization for u_H
        temp = (self.gamma-1) *  mu * self.protonmass / self.boltzmann * u_h
        print "T=",temp
        #Very hot: H and He both fully ionized
        yhelium = (1 - self.hy_mass) / (4 * self.hy_mass)
        nHp = 1.0
        nHepp = yhelium
        ne = nHp + 2.0 * nHepp

        #Free-free cooling rate
        LambdaFF = 1.42e-27 * np.sqrt(temp) * (1.1 + 0.34 * np.exp(-(5.5 - np.log(temp))**2 / 3)) * (nHp + 4 * nHepp) * ne

	    # Inverse Compton cooling off the microwave background
	    #LambdaCmptn = 5.65e-36 * ne * (temp - 2.73 * (1. + self.redshift)) * pow(1. + redshift, 4.) / nH

        return LambdaFF


def get_tescari_rhoHI(star,bar,rho_phys_thresh=0.1):
    """Get a neutral hydrogen density in cm^-2
    applying the correction in eq. 1 of Tescari & Viel
    Parameters:
        bar = a baryon type from an HDF5 file.
        rho_phys_thresh - physical SFR threshold density in hydrogen atoms/cm^3
                        - 0.1 (Tornatore & Borgani 2007)
                        - 0.1289 (derived from the SH star formation model)
    Returns:
        nH0 - the density of neutral hydrogen in these particles in atoms/cm^3
    """
    inH0=np.array(bar["NeutralHydrogenAbundance"],dtype=np.float64)
    #Convert density to hydrogen atoms /cm^3: internal gadget density unit is h^2 (1e10 M_sun) / kpc^3
    irho=np.array(bar["Density"],dtype=np.float64)*(star.UnitMass_in_g/star.UnitLength_in_cm**3)*star.hubble**2/(star.protonmass/star.hy_mass)
    #Default density matches Tescari & Viel and Nagamine 2004
    dens_ind=np.where(irho > rho_phys_thresh)
    #UnitCoolingRate_in_cgs=UnitMass_in_g*(UnitVelocity_in_cm_per_s**3/UnitLength_in_cm**4)
    #Note: CoolingRate is really internal energy / cooling time = u / t_cool
    # HOWEVER, a CoolingRate of zero is really t_cool = 0, which is Lambda < 0, ie, heating.
    #For the star formation we are interested in y ~ t_star/t_cool,
    #So we want t_cool = InternalEnergy/CoolingRate,
    #except when CoolingRate==0, when we want t_cool = 0
    icool=np.array(bar["CoolingRate"],dtype=np.float64)
    ienergy=np.array(bar["InternalEnergy"],dtype=np.float64)
    cool=icool[dens_ind]
    ind=np.where(cool == 0)
    #Set cool to a very large number to avoid divide by zero
    cool[ind]=1e99
    tcool = ienergy[dens_ind]/cool
    #Convert from internal time units, normally 9.8x10^8 yr/h to s.
    tcool *= (star.UnitLength_in_cm/star.UnitVelocity_in_cm_per_s)/star.hubble # Now in s
    fcold=star.cold_gas_frac(irho[dens_ind],tcool,rho_phys_thresh)
    #Adjust the neutral hydrogen fraction
    inH0[dens_ind]=fcold
    #Calculate rho_HI
    nH0=irho*inH0
    #Now in atoms /cm^3
    return nH0

class YajimaRT:
    """Neutral hydrogen density with a self-shielding correction as suggested by Yajima Nagamine 2012 (1112.5691)
    This is just neutral over a certain density."""
    def __init__(self, redshift, hubble=0.71):
        self.redshift = redshift
        self.hubble=hubble
        #Internal gadget mass unit: 1e10 M_sun/h in g/h
        self.UnitMass_in_g=1.989e43
        #Internal gadget length unit: 1 kpc/h in cm/h
        self.UnitLength_in_cm=3.085678e21
        #Internal velocity unit : 1 km/s in cm/s
        self.UnitVelocity_in_cm_per_s=1e5
        #proton mass in g
        self.protonmass=1.67262178e-24
        self.hy_mass = 0.76 # Hydrogen massfrac

    def get_yajima_rhoHI(self,bar):
        """Get a neutral hydrogen density with a self-shielding correction as suggested by Yajima Nagamine 2012 (1112.5691)
        This is just neutral over a certain density."""
        inH0=np.array(bar["NeutralHydrogenAbundance"], dtype=np.float32)
        #Convert density to hydrogen atoms /cm^3: internal gadget density unit is h^2 (1e10 M_sun) / kpc^3
        irho=np.array(bar["Density"], dtype=np.float32)*(self.UnitMass_in_g/self.UnitLength_in_cm**3)*self.hubble**2*(self.hy_mass/self.protonmass)
        #Slightly less sharp cutoff power law fit to data
        r2 = 10**-2.3437
        r1 = 10**-1.81844
        dens_ind=np.where(irho > r1)
        inH0[dens_ind]=1.
        del dens_ind
        ind2 = np.where(np.logical_and(irho < r1,irho > r2))
        #Interpolate between r1 and r2
        n=2.6851
        inH0[ind2] = (inH0[ind2]*(r1-irho[ind2])**n+(irho[ind2]-r2)**n)/(r1-r2)**n
        del ind2
        return inH0

    def get_reproc_HI(self, bar):
        """Get a neutral hydrogen density using the fitting formula of Rahmati 2012"""
        return self.get_yajima_rhoHI(bar)

#Opacities for the FG09 UVB from Rahmati 2012.
#IMPORTANT: The values given for z > 5 are calculated by fitting a power law and extrapolating.
#Gray power law was: -1.12e-19*(zz-3.5)+2.1e-18 fit to z > 2.
#gamma_UVB was: -8.66e-14*(zz-3.5)+4.84e-13
#This is clearly wrong, but this model is equally a poor choice at these redshifts anyway.
gray_opac = [2.59e-18,2.37e-18,2.27e-18, 2.15e-18, 2.02e-18, 1.94e-18, 1.82e-18, 1.71e-18, 1.60e-18]
gamma_UVB = [3.99e-14, 3.03e-13, 6e-13, 5.53e-13, 4.31e-13, 3.52e-13, 2.678e-13,  1.81e-13, 9.43e-14]
zz = [0, 1, 2, 3, 4, 5, 6, 7,8]

import scipy.interpolate.interpolate as intp

class RahmatiRT:
    """Class implementing the neutral fraction ala Rahmati 2012"""
    def __init__(self, redshift,hubble = 0.71, fbar=0.17, molec = True):
        self.f_bar = fbar
        self.redshift = redshift
        self.molec = molec
        #Interpolate for opacity and gamma_UVB
        gamma_inter = intp.interp1d(zz,gamma_UVB)
        gray_inter = intp.interp1d(zz,gray_opac)
        print zz
	print gray_opac
	print redshift
        self.gray_opac = gray_inter(redshift)
        self.gamma_UVB = gamma_inter(redshift)
        #Some constants and unit systems
        #Internal gadget mass unit: 1e10 M_sun/h in g/h
        UnitMass_in_g=1.989e43
        #Internal gadget length unit: 1 kpc/h in cm/h
        self.UnitLength_in_cm=3.085678e21
        self.UnitDensity_in_cgs = UnitMass_in_g/self.UnitLength_in_cm**3
        #Internal velocity unit : 1 km/s in cm/s
        self.UnitVelocity_in_cm_per_s=1e5
        #proton mass in g
        self.protonmass=1.67262178e-24
        #self.hy_mass = 0.76 # Hydrogen massfrac
        self.gamma=5./3
        #Boltzmann constant (cgs)
        self.boltzmann=1.38066e-16

        self.hubble = hubble

        self.star = StarFormation(hubble)
        #Physical density threshold for star formation in H atoms / cm^3
        self.PhysDensThresh = self.star.get_rho_thresh()


    def photo_rate(self, nH, temp):
        """Photoionisation rate as  a function of density from Rahmati 2012, eq. 14.
        Calculates Gamma_{Phot}.
        Inputs: hydrogen density, temperature
            n_H

        The coefficients are their best-fit from appendix A."""
        nSSh = self.self_shield_dens(temp)
        photUVBratio= 0.98*(1+(nH/nSSh)**1.64)**-2.28+0.02*(1+nH/nSSh)**-0.84
        return photUVBratio * self.gamma_UVB

    def self_shield_dens(self,temp):
        """Calculate the critical self-shielding density. Rahmati 202 eq. 13.
        gray_opac and gamma_UVB are parameters of the UVB used.
        gray_opac is in cm^2 (2.49e-18 is HM01 at z=3)
        gamma_UVB in 1/s (1.16e-12 is HM01 at z=3)
        temp is particle temperature in K
        f_bar is the baryon fraction. 0.17 is roughly 0.045/0.265
        Returns density in atoms/cm^3"""
        T4 = temp/1e4
        G12 = self.gamma_UVB/1e-12
        return 6.73e-3 * (self.gray_opac / 2.49e-18)**(-2./3)*(T4)**0.17*(G12)**(2./3)*(self.f_bar/0.17)**(-1./3)

    def recomb_rate(self, temp):
        """The recombination rate from Rahmati eq A3, also Hui Gnedin 1997.
        Takes temperature in K, returns rate in cm^3 / s"""
        lamb = 315614./temp
        return 1.269e-13*lamb**1.503 / (1+(lamb/0.522)**0.47)**1.923

    def neutral_fraction(self, nH, temp):
        """The neutral fraction from Rahmati 2012 eq. A8"""
        alpha_A = self.recomb_rate(temp)
        #A6 from Theuns 98
        LambdaT = 1.17e-10*temp**0.5*np.exp(-157809./temp)/(1+np.sqrt(temp/1e5))
        A = alpha_A + LambdaT
        B = 2*alpha_A + self.photo_rate(nH, temp)/nH + LambdaT

        return (B - np.sqrt(B**2-4*A*alpha_A))/(2*A)

    def get_temp(self,bar):
        """Compute temperature (in K) from internal energy.
           Uses: internal energy
                 electron abundance
                 hydrogen mass fraction (0.76)
           Factor to convert U (J/kg) to T (K) : U = N k T / (γ - 1)
           T = U (γ-1) μ m_P / k_B
           where k_B is the Boltzmann constant
           γ is 5/3, the perfect gas constant
           m_P is the proton mass

           μ = 1 / (mean no. molecules per unit atomic weight)
             = 1 / (X + Y /4 + E)
             where E = Ne * X, and Y = (1-X).
             Can neglect metals as they are heavy.
             Leading contribution is from electrons, which is already included
             [+ Z / (12->16)] from metal species
             [+ Z/16*4 ] for OIV from electrons."""
        #convert U (J/kg) to T (K) : U = N k T / (γ - 1)
        #T = U (γ-1) μ m_P / k_B
        #where k_B is the Boltzmann constant
        #γ is 5/3, the perfect gas constant
        #m_P is the proton mass
        #μ is 1 / (mean no. molecules per unit atomic weight) calculated in loop.
        #Internal energy units are 10^-10 erg/g
        ienergy=np.array(bar["InternalEnergy"])*1e10
        #Calculate temperature from internal energy and electron abundance
        nelec=np.array(bar['ElectronAbundance'])
        try:
            hy_mass = np.array(bar["GFM_Metals"][:,0], dtype=np.float32)
        except KeyError:
            hy_mass = 0.76
        mu = nmex.evaluate("4 / (hy_mass * (3 + 4*nelec) + 1)")
        #So for T in K, boltzmann in erg/K, internal energy has units of erg/g
        temp = (self.gamma-1) * self.protonmass / self.boltzmann * nmex.evaluate("mu * ienergy")
        return temp

    def get_rahmati_HI(self, bar):
        """Get a neutral hydrogen density using the fitting formula of Rahmati 2012"""
        #Convert density to atoms /cm^3: internal gadget density unit is h^2 (1e10 M_sun) / kpc^3
        nH=self.get_code_rhoH(bar)
        temp = self.get_temp(bar)
        nH0 = self.neutral_fraction(nH, temp)
        return nH0

    def get_code_rhoH(self,bar):
        """Convert density to physical atoms /cm^3: internal gadget density unit is h^2 (1e10 M_sun) / kpc^3"""
        nH = np.array(bar["Density"])
        conv = np.float32(self.UnitDensity_in_cgs*self.hubble**2/(self.protonmass)*(1+self.redshift)**3)
        #Convert to physical
        nH=nmex.evaluate("conv*nH")
        return nH

    def code_neutral_fraction(self, bar):
        """Get the neutral fraction from the code"""
        return np.array(bar["NeutralHydrogenAbundance"])

    def get_reproc_HI(self, bar):
        """Get a neutral hydrogen *fraction* using values given by Arepo
        which are based on Rahmati 2012 if UVB_SELF_SHIELDING is on.
        Above the star formation density use the Rahmati fitting formula directly,
        as Arepo reports values for the eEOS. """
        nH0 = self.code_neutral_fraction(bar)
        nH = self.get_code_rhoH(bar)
        #Above star-formation threshold, we want a neutral fraction which includes
        #explicitly the amount of gas in cold clouds.
        #Ideally we should compute this fraction, and then do
        #  tcool = self.get_tcool(nH,bar)[ind]
        #  print np.median(tcool)
        #  cold_frac = self.star.cold_gas_frac(nH[ind], tcool,self.PhysDensThresh/0.76)
        #  print np.mean(cold_frac)
        #  ssnH0 = (1-cold_frac)*self.neutral_fraction(nH[ind], temp[ind])+ cold_frac*self.neutral_fraction(nH[ind], 1e4)
        #But the cooling time reported by the code is not quite what we want here,
        #because it uses the internal energy reported by U, whereas we really want
        #the cooling time using the energy for the hot phase, at a given density.
        #So just assume that at the threshold all gas is in cold clouds.
        #In reality, this should be about 90% of gas in cold clouds, so
        #we will overpredict the neutral fraction by a small amount.
        ind = np.where(nH > self.PhysDensThresh/0.76)
        ssnH0 = self.neutral_fraction(nH[ind], 1e4)
        nH0[ind] = ssnH0*(1-self.molec*self.get_H2_frac(nH[ind]))
        return nH0

    def get_H2_frac(self,nHI):
        """Get the molecular fraction for neutral gas from the ISM pressure:
           only meaningful when nH > 0.1."""
        fH2 = 1./(1+(0.1/nHI)**(0.92*5./3.)*35**0.92)
        return fH2

    def get_tcool(self, nH, bar):
        """
        Get the cooling time for a particle from the snapshot.

        UnitCoolingRate_in_cgs=UnitMass_in_g*(UnitVelocity_in_cm_per_s**3/UnitLength_in_cm**4)

        CoolingRate is really internal energy / cooling time = u / t_cool
         HOWEVER, a CoolingRate of zero is really t_cool = 0, which is Lambda < 0, ie, heating.
        For the star formation we are interested in y ~ t_star/t_cool,
        So we want t_cool = InternalEnergy/CoolingRate,
        except when CoolingRate==0, when we want t_cool = 0

        returns tcool in seconds
        """
        #GFM_Cooling rate is already in s/erg.
        cool=np.array(bar["GFM_CoolingRate"], dtype=np.float32)
        ienergy=np.array(bar["InternalEnergy"], dtype=np.float32)
        ind=np.where(cool < 0)
        #Set cool to a very large number to avoid divide by zero
        tcool = np.zeros_like(ienergy, dtype=np.float32)
        fact = -(0.76/self.protonmass)*nH[ind]
        tcool[ind] = ienergy[ind]*1e10/(cool[ind]*fact)
        return tcool
