import numpy as np 
import matplotlib.pyplot as plt

data = np.load('fab_efficiency_tols.npy')

duty_cycle = data[0,:,:]
zone_height = data[1,:,:]

zh_mm = zone_height * 1e3

efficiency = data[2,:,:]

max_efficiency = np.max(efficiency)

plt.figure(figsize=(8,8))
plt.pcolormesh(duty_cycle,zh_mm,efficiency/max_efficiency, cmap='binary')
plt.colorbar(label='Intensity / Max Intensity')
plt.contour(duty_cycle,zh_mm, efficiency/max_efficiency, levels=[0.8,0.9,0.95,0.99], \
            colors=['xkcd:bright red','xkcd:cyan', 'xkcd:bright green','xkcd:sun yellow'])
plt.scatter(0.5, 2.15, marker='+', color='xkcd:bright magenta', s=100, label='APRA ZP')

plt.legend(fontsize=12)
plt.xlim(0.3,0.7)
plt.xlabel('Duty Cycle (s:P)', fontsize=18)
plt.xticks(fontsize=12)
plt.ylabel(r'Zone Height ($\mu m$)', fontsize=18)
plt.yticks(fontsize=12)

plt.savefig('fabricaton_efficiency_tolerances.pdf',dpi=2000)
