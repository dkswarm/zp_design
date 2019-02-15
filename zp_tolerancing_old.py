from numpy import *
import matplotlib.pyplot as plt

'''
Zone Plate Tolerancing Script -- Making the Figure 4.10 plot from Menz Thesis
'''
#######################################################
arcsec = pi/180/3600

'''
Zone Plate Operating Characteristics -- all unit in meters
'''

# Operating wavelength  
lam = 1240./1487*10**-9    # 1240./1254*10**-9
# Natural linewidth of the line (Al Ka = 0.49 eV, Mg Ka = 0.36 eV)
dlam = 1240./(1487.)**2*0.49*10**-9
# Focal length of operation
f_op = 122.

# ZP Radii of interest
rad = logspace(-3,2,1001)

'''
Tolerancing functions outlined in Table 4.4
'''
# Limiting resolution
def res(r,lam = lam,f = f_op):
    return 1.22*lam*f/(2*r)/arcsec

# Focal depth
def focal_depth(r,df,f = f_op):
    return df*r/f**2/arcsec

# Chromatic aberration
def chrom(r,lam = lam,dlam = dlam, f = f_op):
    return dlam/lam*r/f/arcsec

# Radial offset from center
def rad_error(r,dr):
    return dr/r/arcsec

# Spherical aberration 
def spher_abb(r,f = f_op):
    return 0.5*r**3/f**3/arcsec

# Astigmatism and Field Curvature
def astig_fc(r,theta,phi = 0.,f = f_op):
    return 0.5*r/f*theta**2*sqrt(1 + sin(phi)**2)/arcsec

# Coma
def coma(r,theta,phi = 0.,f = f_op):
    return 0.25*r**2/f**2*theta*sqrt(1 + 8*cos(phi)**2)

'''
Making the Menz Fig. 4.10 plot for the equivalent ZP we're interested in.
'''

def plot_tols(rad):
    plt.ion()
    plt.figure(figsize = (10,10))
    # Plotting the limiting zone plate resolution.
    plt.loglog(rad,res(rad),label = 'ZP Res.',color = 'r')
    # Plotting the chromatic aberration expected for the wavelength of interest.
    plt.loglog(rad,chrom(rad),label = 'Chromatic',color = 'y')
    # Plotting the defocus expected for a fraction of the focal length
    plt.loglog(rad,focal_depth(rad,10**-5*f_op),label = 'Def. 10^-5',color = 'g',linestyle = 'dashed')
    plt.loglog(rad,focal_depth(rad,10**-4*f_op),label = 'Def. 10^-5',color = 'g',linestyle = 'dotted')
    # Plotting the impact of astigmatism
    plt.loglog(rad,astig_fc(rad,0.5*pi/180),label= 'Astig. 0.5 deg.',color = 'b',linestyle = 'dashed')
    plt.loglog(rad,astig_fc(rad,5.0*pi/180),label= 'Astig. 5.0 deg.',color = 'b',linestyle = 'dotted')
    # Plotting the impact of coma
    plt.loglog(rad,coma(rad,0.5*pi/180),label = 'Coma 0.5 deg.',color = 'm',linestyle = 'dashed')
    plt.loglog(rad,coma(rad,5.0*pi/180),label = 'Coma 5.0 deg.',color = 'm',linestyle = 'dotted')
    # Plotting the impact of spherical aberration
    plt.loglog(rad,spher_abb(rad),label = 'Sph. Ab.', color = 'k',linestyle = 'dotted')
    
    plt.xlabel('Zone Plate Radius (m)')
    plt.ylabel('Angular Resolution (arcseconds)')
    plt.title('Angular Resolution Impacts on ZP Performance for\nFocal Length = ' + str(f_op) + 'm ZPs')

    plt.xlim(10**-3,70.)
    plt.ylim(10**-4,100.)
    plt.legend()
    
