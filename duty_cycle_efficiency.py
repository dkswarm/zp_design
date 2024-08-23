import numpy as np
import matplotlib.pyplot as plt
import define_zp as dzp
import zp_efficiency as zpe

'''
All measurements in units of millimeters
'''

materials = ['W','Si']
m = np.array([1,2,3]) # array of orders for efficiency calculation
duty_cycles = np.linspace(0,1,100)    # array of duty cyles for zone plates

# Zone plate measurements taken from 190611_ExampleZPTolerances.py

f = 40e3    # focal length for PSU beamline
r_max = 200 # maximum radius of zone plate
r_min = 150 # minimum radius of zone plate (segmented)
theta_seg = (24/180)*np.pi  # segmented zone plate arc in radians

# initializing array to store efficiencies for each zone plate in each order
efficiencies = np.zeros((m.size, duty_cycles.size))
m0efficiencies = np.zeros(duty_cycles.size)

z_h = 2.15e-3   # zone height from APRA proposal section 3.2, page 8, par 1
for material in materials:
    print(material)
    for index, duty_cycle in enumerate(duty_cycles):
        print('Initiating Zone Plate with duty cycle %f'%(duty_cycle))
        zoneplate = dzp.zone_plate(f=f, r_max=r_max, r_min=r_min, sp_ratio=duty_cycle, \
                                   material=material, theta_seg=theta_seg)
        print('Calculating Efficiencies')
        efficiencies[:,index] = zpe.m_eff(z_h=z_h, zp=zoneplate, m=m)
        m0efficiencies[index] = zpe.m0_eff(z_h=z_h, zp=zoneplate, m=0)

    efficiencies *= 100
    m0efficiencies *= 100

    plt.figure()
    plt.plot(duty_cycles, m0efficiencies, label='0th Order')
    plt.plot(duty_cycles, efficiencies[0], label='1st Order')
    plt.plot(duty_cycles, efficiencies[1], label='2nd Order')
    plt.plot(duty_cycles, efficiencies[2], label='3rd Order')
    plt.legend()
    plt.yscale('log')
    plt.ylim(1,100)
    plt.xlabel('Duty Cycle (S:P)')
    plt.ylabel('Efficiency (%)')
    plt.title('Efficiencies vs. Duty Cycle: ' + material)
    plt.tight_layout()
    plt.savefig('efficiencies_duty-cycle_' + material + '.png')

    print('Plot saved. Material processing complete.')
