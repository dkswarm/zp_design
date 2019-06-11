import numpy as np
import matplotlib.pyplot as plt

import pdb

# By default, everything is in mm like common raytrace routines

hc = 1239.8 # in nm eV
mg_k_wave = hc/1253.7/10**6    # Wavelength of Mg K alpha_1 emission in mm
al_k_wave = hc/1486.7/10**6    # Wavelength of Al K alpha_1 emission in mm

mg_k_wave_width = mg_k_wave**2/hc*0.36*10**6          # Assuming a natural line width of 0.36 eV
al_k_wave_width = al_k_wave**2/hc*0.49*10**6          # Assuming a natural line width of 0.49 eV

arcsec = np.pi/180/3600

class zone_plate:
    def __init__(self, f, r_max, wave = mg_k_wave, dwave = mg_k_wave_width,order = 1.,r_min = 0., theta_seg = 2*np.pi):
        self.f = f
        self.r_max = r_max
        self.r_min = r_min          # Nonzero for a segmented zone plate -- defines how big an area is covered by the ZP.
        self.theta_seg = theta_seg  # Less than 2*pi for a segmented zone plate -- defines the angular coverage as measured from the +x axis.
        self.wave = wave
        self.dwave = dwave
        self.order = order

        self.material = 'Si'        # Setting the 'opaque' material
        self.sp_ratio = 0.5         # Setting the ratio of space (r_n+1 - r_n) to period (r_n+2 - r_n)

        self.N_zones = self.__compute_N_zones()
        self.r_first = self.__compute_r_first()
        self.dr_min = self.__compute_dr_min()

    def compute_N_zones(self):
        return np.floor(self.r_max**2/(self.wave*self.order)/self.f)

    def compute_r_first(self):
        return np.sqrt(self.f*self.wave*self.order)

    def compute_dr_min(self):
        return self.r_first/2/np.sqrt(self.N_zones)

    __compute_N_zones = compute_N_zones
    __compute_dr_min = compute_dr_min
    __compute_r_first = compute_r_first

    def trace_rays(self,ray_object,order):
        '''
        Defines the interaction with an ArcusRays-style ray object. Computes the 
        diffraction through the zone_plate as a function of order without regard for
        efficiency.
        Assumes rays are on zone plate and that the ZP lies in the XY plane.
        
        '''

        # Computing the relative orientation of the grating and the groove density at each ray point.
        phis = np.arctan2(ray_object.y,ray_object.x) - np.pi/2
        zone_num = np.ceil((ray_object.x**2 + ray_object.y**2)/self.r_first**2)
        ds = self.r_first / np.sqrt(zone_num)

        ray_object.vx = ray_object.vx + (order*ray_object.wave/ds)*np.sin(phis)
        ray_object.vy = ray_object.vy - (order*ray_object.wave/ds)*np.cos(phis)
        ray_object.vz = np.sqrt(1.0 - ray_object.vx**2 - ray_object.vy**2)

        ray_object.zone_num = zone_num