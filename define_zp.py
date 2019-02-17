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
    def __init__(self, f, r_max, wave = mg_k_wave, dwave = mg_k_wave_width,order = 1.):
        self.f = f
        self.r_max = r_max
        self.wave = wave
        self.dwave = dwave
        self.order = order

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