import glob
import os
import numpy as np

# Other dependencies.
import n_lookup
from mendeleev import element

import pdb

index_dir = r'C:\Users\Casey\Software\python_repository\zp_design\cxro_lookup'

# An example of how to use the CXRO look-up: this finds the tabulated indices of refraction for Silicon Nitride from 30 eV - 10 keV.
#sn = n_lookup.index()
#data = sn.calculate_energy_scan(30.,10000.,100)

Zs = range(1,92)
# First index: energy, delta, beta
# Second index: atomic number
# Third index: values of the first index
material_database = np.empty((3,len(Zs),500))

material = n_lookup.index()

for i in range(len(Zs)):
    #pdb.set_trace()
    material.chemical_formula = element(Zs[i]).symbol 
    temp_data = np.asarray(material.calculate_energy_scan(1000.,1500.,500))
    material_database[:,i,:] = np.transpose(temp_data)

np.save(index_dir + '\z_index.npy',material_database)