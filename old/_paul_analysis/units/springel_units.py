import numpy as np
import units
from astropy.cosmology import LambdaCDM
import astropy.units as u


def hubble_function( a, H0=67.74, Omega_M0=0.3089, Omega_L=0.6911 ):
    return H0 * np.sqrt(  Omega_M0 / a**3 + Omega_L  )


def critical_density( a, H0=67.74, Omega_M0=0.3089, Omega_L=0.6911  ):
    gravitational_constant = 6.67e-8    # cm^3 / g / s^2
    H = hubble_function( a , H0=H0, Omega_M0 = Omega_M0, Omega_L = Omega_L) # km/s/Mpc
    H /= 3.086e19       # units of [s^-1]
    
    crit_val = 3.0 * H**2 / (8.0 * 3.14159 * gravitational_constant )   # units of g / cm^3
    crit_val /= 1.67e-24                                                # units of [#/cm^-3]

    return crit_val


def age_from_a(a, a0=1, H0=70.4, Om0=0.2726, Ode0=0.7274):
#a    cosmo = LambdaCDM( H0=H0, Om0=Om0, Ode0=Ode0)
#    return cosmo.age( 1/a - 1.0 ) 

      OmegaLambda = Ode0
      Omega0      = Om0

      factor1 = 2.0 / (3.0 * np.sqrt(OmegaLambda));

      term1 = np.sqrt(OmegaLambda / Omega0) * a0**1.5
      term2 = np.sqrt(1 + OmegaLambda / Omega0 * a0**3.0)
      factor2 = np.log(term1 + term2);

      t0 = factor1 * factor2;

      term1 = np.sqrt(OmegaLambda / Omega0) * a**1.5
      term2 = np.sqrt(1 + OmegaLambda / Omega0 * a**3.0)
      factor2 = np.log(term1 + term2);

      t1 = factor1 * factor2;

      result = t0 - t1;

      hubble            = 3.2407789e-18
      hubbleparam       = 0.7
      sec_per_megayear  = 3.15576e13

      ft = result / (hubble * hubbleparam);  	#/* now in seconds */
      ft /= sec_per_megayear * 1000;   		#a  /* now in gigayears */

      return ft 


   


def gas_mu(num_e):
    XH=0.76; # we track this with metal species now, could do better...
    yhelium=(1.-XH)/(4.*XH);
    return (1.+4.*yhelium)/(1.+yhelium+num_e);

def gas_code_to_temperature(gas_u,gas_nume, gamma=5.0/3.0):
#    mu = gas_mu(gas_nume);
#    MeanWeight= mu*units.constants["PROTONMASS"]
#    Temp= MeanWeight/units.constants["BoltzMann_ergs"] * (gamma-1.0) * gas_u * 1.e10
    return gas_mu(gas_nume) * units.constants["PROTONMASS"] /  units.constants["BoltzMann_ergs"] * (gamma-1.0) * gas_u * 1.e10
     #Temp

def gas_code_to_cgs_density(rho):
    return rho * 405.39

def gas_cs_effective_eos(u, q_eos=1.0):
    ## returns gas sound speed in km/s (for old effective gadget equation-of-state)
    ## u is in GADGET units = (kpc/Gyr)^2 = 0.957*(km/s)^2
    ## actually no, time unit in gadget isnt Gyr, exactly (velocity unit *is* km/s
    g_gamma = 5./3.
    u_min = 100. ## 10^4 K -- *roughly* 10 km/s, though exact conversion depends on mean molecular weight, etc.
    cs_gas = np.sqrt(g_gamma*(g_gamma-1.) * (q_eos*u + (1.-q_eos)*u_min))
    return cs_gas


def gas_xray_brems(mass_in_gadget_units, u_in_gadget_units, rho_in_gadget_units, num_e, num_h):
    ## returns gas x-ray bremstrahhlung luminosity (x-ray line cooling is separate)
    protonmass = 1.6726e-24;
    brem_normalization= 1.2e-24;
    m = mass_in_gadget_units * 1.989e43; ## convert to g
    u = u_in_gadget_units;
    
    MeanWeight = gas_mu(num_e) * protonmass;
    keV = gas_temperature(u, num_e, keV=1);
    density = rho_in_gadget_units * 6.76991e-22 ## convert to g/cm3
    
    # total bremstrahhlung luminosity (integrated over all wavelengths)
    brem_lum = brem_normalization * (m/MeanWeight) * (density/MeanWeight) * \
        np.sqrt(keV) * num_e * (1.-num_h);
    
    # thermal bremstrahhlung spectrum goes as I_nu~exp(-h*nu/kT):
    # crudely take Lx as the fraction radiated above > some freq in keV (nu_min_keV),
    # which gives:
    nu_min_keV = 0.5
    xray_lum = brem_lum * np.exp(-nu_min_keV/keV)
    
    return xray_lum
