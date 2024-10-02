import numpy as np
import matplotlib.pyplot as plt
import zp_efficiency as zpe
import define_zp as dzp


# Zone plate measurements taken from 190611_ExampleZPTolerances.py

print('Setting values')
f = 40e3    # focal length for PSU beamline
r_max = 200 # maximum radius of zone plate
r_min = 150 # minimum radius of zone plate (segmented)
theta_seg = (24/180)*np.pi  # segmented zone plate arc in radians

duty_cycles = np.linspace(0.3,0.8,51)
zone_height_nominal = 2.15e-3
zone_heights = np.linspace(zone_height_nominal - 1e-3, zone_height_nominal + 1e-3, 201)

efficiencies = np.zeros((duty_cycles.size, zone_heights.size))

print('Calculating efficiencies')

for index, duty_cycle in enumerate(duty_cycles):
    print(duty_cycle)
    zp = dzp.zone_plate(f=f, r_max=r_max, r_min=r_min, theta_seg=theta_seg, sp_ratio=duty_cycle)
    efficiencies[index,:] = zpe.m_eff(z_h=zone_heights, zp=zp, m=1)

print('Saving data')

zhlayer, dclayer = np.meshgrid(zone_heights,duty_cycles)

# The meshgrid swaps the indices from what we established in the efficiencies array
np.save('fab_efficiency_tols.npy', np.asarray([dclayer,zhlayer,efficiencies]))
print('Efficiencies saved.')