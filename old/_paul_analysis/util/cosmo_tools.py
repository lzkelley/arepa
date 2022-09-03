#!/usr/bin/env python
"""Functions to calculate cosmological definitions in the Illustris cosmology."""

__author__ = "Ekta Patel and contributing authors"
__copyright__ = "Copyright 2015, The Authors"
__credits__ = ["Ekta Patel and contributing authors"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Ekta Patel"
__email__ = "ektapatel@email.arizona.edu"
__status__ = "Beta -- forever."

import math
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import norm
from scipy.integrate import simps

def time(z):
    '''
    Calculate lookback time for a flat cosmology
    '''
    if z == 0.:
        time = 0. 
    else: 
        h = 0.704
        H0 = h*100
        OM = 0.2726
        OL = 0.7274
        H0 = H0 * 3.241e-20 / 3.171e-17 # final units of[Gyr ^-1]

        def f(z):
            return 1/ (H0*(1+z)*np.sqrt(OM*(1+z)**3 + OL)) 

        zs = np.arange(0., z, 1e-4)
        y = f(zs)
        time = simps(y,zs)
    return time #in Gyrs

# delta_c taken from Bryan and Norman (1998)
def delta_vir(omegam):
    x = omegam - 1.
    deltac = 18*np.pi**2 + 82*x -39*x**2
    return deltac/omegam

# taken from van der Marel 2012, equation A1 ("THE M31 VELOCITY VECTOR. II. RADIAL ORBIT TOWARD THE MILKY WAY AND IMPLIED LOCAL GROUP MASS")
def r_vir(omegam, h, Mvir):
    a = 206./h
    b = delta_vir(omegam) * omegam / 97.2
    c = Mvir * h / 1e12
    return a * b**(-1./3.) * c**(1./3.)

def v_vir(omegam, h, Mvir):
    G = 4.302e-6 #kpc/Msun (km/s)^2
    rvir = r_vir(omegam, h, Mvir)
    return np.sqrt(G*Mvir/rvir)

# vdM 2012, equation A5
def q(omegam):
    return 200. / (delta_vir(omegam) * omegam)

# taken from Klypin 2011 on Bolshoi simulations, equation 10
def c_vir(h, Mvir):
    a = Mvir *h/ 1e12
    return 9.60*(a)**(-0.075)

# vdM 2012 definition for cvir = rvir/rs
def r_s(omegam, h, Mvir, c=False):
    rv = r_vir(omegam, h, Mvir)
    cv = c_vir(h, Mvir)
    if c: 
        cv = c 
    return rv/cv

# vdM 2012, equation A3
def f(x):
    a = np.log(1+x) 
    b = x / (1+x)
    return a - b

# from Battaner on NED IPAC
def NFW_Vmax(omegam, h, Mvir, cvir, r):
    x = r/r_vir(omegam, h, Mvir)
    top = np.log(1+cvir*x) - cvir*x/(1+cvir*x)
    v200 = h * c_200(cvir, omegam, h, Mvir) * r_s(omegam, h, Mvir, c=cvir)
    return v200 * (top/ (x * f(cvir)))**0.5

# from Klypin et al. (2011), equation 3
def NFW_Vmax_Klypin(omegam, h, Mvir, cvir):
    G = 4.302e-6 #kpc * (km/s)^2 / Msun
    xmax = 2.15
    rvir = r_vir(omegam, h, Mvir)
    a = G * f(xmax)/f(cvir) * cvir/xmax * Mvir**(1./3.)/rvir
    return a**(0.5) * Mvir**(1./3.)

# http://glowingpython.blogspot.com/2012/01/fixed-point-iteration.html
# fp iteration reference in link above; vdM 2012, equation A6
def c_200(cvir,omegam, h, Mvir):
    c0 = cvir
    qi = q(omegam)
    fvir = f(cvir)
    tol = 1e-7
    maxiter=200
    e = 1
    itr = 0
    cp = []
    while(e > tol and itr < maxiter):
        c = cvir * (f(c0)/(qi*fvir))**(1./3.)
        e = norm(c0 - c)
        c0 = c
        cp.append(c0)
        itr = itr + 1.
    return c

# vdM 2012 definition for c200=r200/rs
def r_200(omegam, h, Mvir):
    cvir = c_vir(h, Mvir)
    c200 = c_200(cvir, omegam, h, Mvir)
    rs = r_s(omegam, h, Mvir)
    return c200*rs

# vdM 2012, equation A7
def M200(omegam, h, Mvir):
    cvir = c_vir(h, Mvir)
    c200 = c_200(cvir, omegam, h, Mvir)
    f200 = f(c200)
    fvir = f(cvir)
    return f200 * Mvir / fvir

# vdM 2012, equation A7
def mass_ratio(omegam, h, Mvir):
    ''' Returns ratio of M200 to Mvir'''
    m200 = M200(omegam, h, Mvir)
    return m200/Mvir

# vdM 2012, equation A11
def ascale(omegam, h, Mvir, a200=False):
    ''' calculates virial hernquist scale radius from the NFW profile; 
    if 200=True calculated the 200 hernquist scale radius'''

    rs = r_s(omegam, h, Mvir)
    cvir = c_vir(h, Mvir)
    a = 1./np.sqrt(2*f(cvir))
    b = 1./cvir
    scale = rs/(a-b)

    if a200:
        c200 = c_200(cvir, omegam, h, Mvir)
        a = 1./np.sqrt(2*f(c200))
        b = 1./c200
        scale = rs/(a-b)
    return scale

# vdM 2012, equation A12
def Mh2Mvir(omegam, h, Mvir, a200=False):
    ''' Returns ratio of hernquist mass to virial mass '''

    avir = ascale(omegam, h, Mvir)
    rs = r_s(omegam, h, Mvir)
    cvir = c_vir(h, Mvir)
    ratio = (avir/rs)**2 / (2*f(cvir))
    
    if a200:
        a_200 = ascale(omegam, h, Mvir, a200=True)
        ratio = (a_200/rs)**2 / (2*f(cvir))
    return ratio

def snapnum2z(snapnum):
    snapnums = np.arange(0,136,1)
    zs =[4.67730000e+01,4.45620000e+01,4.24540000e+01,4.06400000e+01
    ,3.87120000e+01,3.68750000e+01,3.51220000e+01,3.36140000e+01
    ,3.20120000e+01,3.04840000e+01,2.90270000e+01,2.76380000e+01
    ,2.64420000e+01,2.51720000e+01,2.39610000e+01,2.28060000e+01
    ,2.18120000e+01,2.07560000e+01,1.97490000e+01,1.87890000e+01
    ,1.79630000e+01,1.70850000e+01,1.62480000e+01,1.54500000e+01
    ,1.47630000e+01,1.40340000e+01,1.33380000e+01,1.26750000e+01
    ,1.20420000e+01,1.14970000e+01,1.09190000e+01,1.03670000e+01
    ,9.99660000e+00,9.84140000e+00,9.38880000e+00,9.00230000e+00
    ,8.90800000e+00,8.44950000e+00,8.01220000e+00,7.59510000e+00
    ,7.23630000e+00,7.00540000e+00,6.85510000e+00,6.49160000e+00
    ,6.14490000e+00,6.01080000e+00,5.84660000e+00,5.52980000e+00
    ,5.22760000e+00,4.99590000e+00,4.93940000e+00,4.66450000e+00
    ,4.42800000e+00,4.17680000e+00,4.00790000e+00,3.93730000e+00
    ,3.70880000e+00,3.49090000e+00,3.28300000e+00,3.08480000e+00
    ,3.00810000e+00,2.89580000e+00,2.73310000e+00,2.57730000e+00
    ,2.44420000e+00,2.31610000e+00,2.20790000e+00,2.10330000e+00
    ,2.00200000e+00,1.90410000e+00,1.82270000e+00,1.74360000e+00
    ,1.66670000e+00,1.60420000e+00,1.53120000e+00,1.47200000e+00
    ,1.41410000e+00,1.35760000e+00,1.30240000e+00,1.24850000e+00
    ,1.20630000e+00,1.15460000e+00,1.11420000e+00,1.07450000e+00
    ,1.03550000e+00,9.97300000e-01,9.87850000e-01,9.50530000e-01
    ,9.23000000e-01,8.86900000e-01,8.51470000e-01,8.16710000e-01
    ,7.91070000e-01,7.57440000e-01,7.32640000e-01,7.00110000e-01
    ,6.76110000e-01,6.44640000e-01,6.21430000e-01,5.98540000e-01
    ,5.75980000e-01,5.46390000e-01,5.24570000e-01,5.03050000e-01
    ,4.81830000e-01,4.60920000e-01,4.40300000e-01,4.19970000e-01
    ,3.99930000e-01,3.80170000e-01,3.60690000e-01,3.47850000e-01
    ,3.28830000e-01,3.10070000e-01,2.91580000e-01,2.73350000e-01
    ,2.61340000e-01,2.43540000e-01,2.25990000e-01,2.14420000e-01
    ,1.97280000e-01,1.80390000e-01,1.69250000e-01,1.52750000e-01
    ,1.41880000e-01,1.25760000e-01,1.09870000e-01,9.94010000e-02
    ,8.38840000e-02,7.36620000e-02,5.85070000e-02,4.85240000e-02
    ,3.37240000e-02,2.39740000e-02,9.52180000e-03,0.00000000e+00]

    z = zs[np.where(snapnums == snapnum)[0][0]]
    return z

if __name__ == "__main__":

    print snapnum2z(42)
    zs = np.arange(0.,10., .01)
    plt.figure()
    plt.plot(zs, [time(z) for z in zs], '-')
    plt.xlabel('z')
    plt.ylabel('lookback time [Gyr]')
    plt.savefig('redshift.pdf')

    print [time(i) for i in [0.0,0.01,0.2,0.4,0.7,1.1,1.9,4.0]]
    h = 0.704
    om = 0.2726
    M = 1e12
    import nfw as nfw
    print NFW_Vmax(om,h,M,12.,2.15*r_s(om,h,M,c=12.)), NFW_Vmax_Klypin(om,h,M,12.)
    print nfw.NFW(M,2.15*r_s(om,h,M,c=12.), 12.).v_rot()

    print 'q:', q(om)
    print 'delta_vir:', delta_vir(om)
    print 'r_s:', r_s(om, h, M)
    print 'c200:', c_200(10., om, h, M)
    print 'cvir:', c_vir(h, M)
    print 'r200:', r_200(om, h, M), 'rvir:', r_vir(om, h, M)    
    print 'a200:', a(om, h, M, a200=True), 'avir:',ascale(om, h, M)

    print 'example for vdM 2012 paper appendix using cvir=10:'
    print 'c200:',c_200(10., om, h, M)
    print 'M200/Mvir:', mass_ratio(om ,h ,M)
    print 'a200/rs:', ascale(om, h, M, a200=True)/r_s(om, h, M)
    print 'avir/rs:', ascale(om, h, M)/r_s(om, h, M)
    print 'M_H,vir/ Mvir:', Mh2Mvir(om, h, M)
    print 'M_H,200/Mvir:', Mh2Mvir(om, h, M, a200=True)

