#!/usr/bin/env python
"""Routines to create representations of dark matter halo profiles including: Plummer sphere, Isothermal sphere, Hernquist bulge, NFW"""

__author__ = "Ekta Patel and contributing authors"
__copyright__ = "Copyright 2015, The Authors"
__credits__ = ["Ekta Patel and contributing authors"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Ekta Patel"
__email__ = "ektapatel@email.arizona.edu"
__status__ = "Beta -- forever."

import numpy as np
from util.cosmo_tools import *

class Plummer:
    def __init__(self, Mvir, r, a):
        """
        Inputs: Mvir: total mass (solar mass)
                a: Plummer length scale (kpc)
                r: radius (kpc)
        """
        self.Mvir = Mvir
        self.a = a
        self.r = r 
        self.G = 4.302e-6 # kpc/Msun (km/s)^2 

    def density(self):
        M = self.Mvir
        a = self.a 
        r = self.r 
        return 3*M/(4*np.pi*a**3) * (1+(r/a)**2)**(-2.5)

    def potential(self):
        M = self.Mvir
        a = self.a
        G = self.G
        r = self.r
        return - G*M/ np.sqrt(r**2 + a**2) 

    def v_esc(self):
        r = self.r
        phi = self.potential()
        return np.sqrt(-2*phi)
        
    def a_grav(self, position):
        ''' Output: 3D acceleration at a given 3D position, useful for orbit integrators'''

        M = self.Mvir
        a = self.a
        G = self.G
        x,y,z = position
        rr = np.sqrt(x**2 + y**2 + z**2)
        ax = (G * M * x**2)/(rr * (rr**2 + a**2)**(1.5) )
        ay = (G * M * y**2)/(rr * (rr**2 + a**2)**(1.5) )
        az = (G * M * z**2)/(rr * (rr**2 + a**2)**(1.5) )
        acc = [ax, ay, az]
        return acc

    def mass(self):
        ''' Mass enclosed within some input radius '''
        M = self.Mvir
        a = self.a
        r = self.r
        mass_enc = M*r**3/ (r**2 + a**2)**(1.5)
        return mass_enc

    def v_rot(self):
        r = self.r
        G = self.G
        M = self.mass()
        return np.sqrt(G*M/r)

    def Plummer_scale_length(self, M50):
        ''' If the Plummer scale length (a) is unknown, it can be determined by the virial mass and the mass enclosed within 50 kpc '''
        M = self.Mvir
        r = 50.
        a = np.sqrt(((M*(r**3))/M50)**(2/3.) - r**2)
        return a 

class Isotherm:
    def __init__(self, r, vc):
        """
        Inputs: r: radius (kpc)
                vc: circular velocity at a given position (i.e. solar circle) [km/s]
        """
        self.r = r
        self.G = 4.498768e-6 #(kpc^3/Gyr^2/Msol)
        self.vc = vc
        self.conv = 3.241e-17 / (3.171e-17) #km/s to kpc/Gyr

    def potential(self):
        r = self.r
        vc = self.vc
        return - vc**2. * np.log(r)

    def density(self):
        k = self.conv
        G = self.G
        r = self.r
        vc = self.vc
        return (vc*k)**2./ (4. * np.pi * G * r**2.)
    
    def mass(self):
        vc = self.vc
        G = self.G
        r = self.r
        k = self.conv
        return (vc*k)**2. * r /G

    def v_esc(self):
        phi = self.potential()
        return np.sqrt(-2.*phi)

    def a_grav(self, position):
        r = self.r
        vc = self.vc
        x,y,z = pos_vec
        rr = np.sqrt(x**2 + y**2 + z**2)
        ax = x * (vc * 3.24e-17)**2 /rr**2
        ay = y * (vc * 3.24e-17)**2 /rr**2
        az = z * (vc * 3.24e-17)**2 /rr**2
        acc = [ax, ay, az]
        return acc

class Hernquist:
    def __init__(self, Mvir, r, a):
        """
        Inputs: Mvir: total mass (solar mass)
                a: Hernquist length scale (kpc)
                r: radius (kpc)
        """
        self.Mvir = Mvir
        self.a = a
        self.r = r 
        self.G = 4.498768e-6 

    def density(self):
        M = self.Mvir
        r = self.r
        a = self.a
        return M*a / (2.*np.pi*r*(r+a)**3.)

    def potential(self):
        G = self.G
        M = self.Mvir
        a = self.a
        r = self.r
        return -G*M /(r+a)

    def mass(self):
        M = self.Mvir
        r = self.r
        a = self.a
        return M*r**2. / (r+a)**2.

    def v_esc(self):
        phi = self.potential()
        return np.sqrt(-2.*phi)

    def v_rot(self):
        G = self.G
        M = self.mass()
        r = self.r
        return np.sqrt(G*M/r)

    def a_grav(self, position):
        G = self.G
        M = self.Mvir
        a = self.a
        r = self.r
        x,y,z = position
        rr = np.sqrt(x**2 + y**2 + z**2.)
        ax = G*M* x**2 /(rr**2 * (rr+a)**2)
        ay = G*M* y**2 /(rr**2 * (rr+a)**2)
        az = G*M* z**2 /(rr**2 * (rr+a)**2)
        acc = [ax, ay, az]
        return acc
        
    def Hernquist_scale_length(self, M50):
        ''' If the Hernquist scale length (a) is unknown, it can be determined by the virial mass and the mass enclosed within 50 kpc '''
        M = self.Mvir
        r = 50.
        a = np.sqrt(Mvir * r**2 / M50kpc) - r
        return a 

class NFW:
    def __init__(self, Mvir, r, cvir):
        """
        Inputs: Mvir (solar mass)
        r: radius (kpc)
        c: r_vir/r_scale
        """
        self.Mvir = Mvir
        self.r = r
        self.G = 4.302e-6 # kpc/Msun (km/s)^2 
        self.cvir = cvir
        self.x = r/r_s(0.2726, 0.704, Mvir, c=cvir)

    def density(self):
        Mvir = self.Mvir
        cvir = self.cvir
        rs = r_s(0.2726, 0.704, Mvir)
        rhos = Mvir/ (4*np.pi*f(cvir)*rs**3)
        x = self.x       
        return rhos/(x*(1+x)**2)

    def mass(self):
        cvir = self.cvir
        Mvir = self.Mvir
        x = self.x
        return Mvir*f(x)/f(cvir)

    def potential(self):
        G = self.G
        Mvir = self.Mvir
        rs = r_s(0.2726, 0.704, Mvir)
        x = self.x
        cvir = self.cvir
        r = self.r
        return -G*Mvir/f(cvir) *np.log(1+x)/r

    def v_esc(self):
        phi = self.potential()
        return np.sqrt(-2*phi)

    def v_rot(self):
        G = self.G
        m = self.mass()
        r = self.r
        return np.sqrt(G*m/r)

    def a_grav(self, position):
        G = self.G
        Mvir = self.Mvir
        cvir = self.cvir
        r = self.r
        x = self.x
        i,j,k = position
        rr = np.sqrt(i**2 + j**2 + k**2)
        A = G * Mvir * f(x)
        B = f(cvir) * rr**4
        ax = A/B * i**2
        ay = A/B * j**2
        az = A/B * k**2
        acc = [ax, ay, az]
        return acc
        
