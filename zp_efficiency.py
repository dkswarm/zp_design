import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import root
import mendeleev as mv

from zp_design import define_zp as dzp

import pdb

ind_pointer = 'C:\Users\swarm\Software\python_repository\zp_design\cxro_lookup\z_index.npy'

'''
Look up the index of refraction of the material.
'''
def compute_index(wave,Z,database_fn = ind_pointer):
    ind_data = np.load(database_fn)     # Database is organized by (energy, delta, beta), Z, values
    mat_data = ind_data[:,Z - 1,:]      # Reduces to just the element we want.
    mat_energy_delta_func = interp1d(mat_data[0],mat_data[1],kind = 'cubic')
    mat_energy_beta_func = interp1d(mat_data[0],mat_data[2],kind = 'cubic')

    energy = dzp.hc/10**6/wave          # Convert the wavelength (supposed to be in mm) to eV.
    return mat_energy_delta_func(energy),mat_energy_beta_func(energy)       # Returns delta, beta


'''
Calculate the efficiency of a zone plate of order m.
'''
def m_eff(z_h,zp,m):
    '''
    Computes the diffraction efficiency of a zone plate 'zp' in order 'm' with opaque zone height 'z_h'.
    z_h -- zone height
    zp -- the zone plate
    m - order
    '''
    ratio = zp.sp_ratio
    Z = mv.element(zp.material).atomic_number
    wave,k = zp.wave,2*np.pi/zp.wave
    delta,beta = compute_index(wave,Z)

    return ((np.sin(m*np.pi*ratio)/(m*np.pi))**2)*(1 + np.exp(-2*k*beta*z_h) - 2*np.exp(-k*beta*z_h)*np.cos(k*delta*z_h))


def compute_z_opt(zp,m):
    '''
    Computes the zone height which maximizes the efficiency of a zone plate 'zp' in order 'm'.
    z_h -- zone height
    zp -- the zone plate
    m - order
    '''
    ratio = zp.sp_ratio
    Z = mv.element(zp.material).atomic_number
    wave,k = zp.wave,2*np.pi/zp.wave
    delta,beta = compute_index(wave,Z)

    zhs = np.logspace(-2,1.7,10001)/10**3      # Test heights, which range from 10 nm to 50 microns 

    data = m_eff(zhs,zp,m)
    return zhs[np.argmax(data)]                # Really crude way of getting at the maximum efficiency, but it works...