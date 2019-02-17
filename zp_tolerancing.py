from numpy import *
import matplotlib.pyplot as plt

import pdb

from zp_design import define_zp as dzp

class tol_allocation:
    '''
    This sets an error budget for zone plate operation. The Zone Plate (zp) is defined in the module zp_design, and is a class.
    dz -- the focal depth (displacement/misalignment along the optical axis)
    dr -- the average period error in dr. The Menz thesis calls this radial error.
    offaxis_ang -- the angle between the zone plate normal and the chief ray, denoted as theta in Menz thesis.
    sector_ang -- the sector angle, relevant only for partial ZPs of one sector
    '''
    def __init__(self,zp,dz,dr,offaxis_ang,sector_ang = 0.0):
        self.zp = zp
        self.dz = dz
        self.dr = dr
        self.offaxis_ang = offaxis_ang
        self.sector_ang = sector_ang

        self.__run_tol_assignment()
    
    '''
    Series of functions computing the angular resolution contribution of each of the errors.
    '''
    # Computes the diffraction limit for zone plate resolution.
    def compute_diffraction_limit_error(self):
        self.res_err = 1.22*self.zp.wave/(2*self.zp.r_max)

    # Computes the angular resolution error for defocus
    def compute_defocus_error(self):
        self.defocus_err = self.dz*self.zp.r_max/self.zp.f**2

    # Computes the chromatic aberration error.
    def compute_chromatic_aberration_error(self):
        self.chrom_err = (self.zp.dwave/self.zp.wave)*(self.zp.r_max/self.zp.f)

    # Computes the angular contribution from average period error.
    def compute_radial_error(self):
        self.rad_err = self.dr/self.zp.r_max

    # Errors due to spherical aberration.
    def compute_spherical_error(self):
        self.sph_err = 0.5*self.zp.r_max**3/self.zp.f**3

    # Compute the error due to astigmatism and field curvature.
    def compute_astigmatism_fc_error(self):
        self.astig_fc_err = 0.5*(self.zp.r_max/self.zp.f)*(self.offaxis_ang**2)*sqrt(1 + sin(self.sector_ang)**2)

    # Computes the error due to coma.
    def compute_coma_error(self):
        self.coma_err = 0.25*(self.zp.r_max/self.zp.f)**2*self.offaxis_ang*sqrt(1 + 8*cos(self.sector_ang)**2)

    '''
    Function that takes the assigned error budget, computes all of the relevant errors, then derives a angular resolution sum.
    We also then compute an effective focal length error 
    '''
    def run_tol_assignment(self):
        self.__compute_diffration_limit_error()
        self.__compute_defocus_error()
        self.__compute_chromatic_aberration_error()
        self.__compute_radial_error()
        self.__compute_spherical_error()
        self.__compute_astigmatism_fc_error()
        self.__compute_coma_error()

        self.total_err = sqrt(self.res_err**2 + self.defocus_err**2 + self.chrom_err**2 + self.rad_err**2 + self.sph_err**2 + self.astig_fc_err**2 + self.coma_err**2)

        # We next make the approximation that all of these errors behave like a defocus error when operated as a CZP, employing 
        # the idea that an optical system is identical when traced in reverse, and approximating that focusing to a point and being collimated 
        # from a point are roughly similar optical systems. We lastly compute that the source appears to at an effective distance of eff_source_dist, 
        # based on the Taylor expansion of the thin lens equation. The derivation for this simple expression is captured in the digital lab notebook.
        self.eff_f_err = self.zp.f*self.total_err           # Inverting 
        self.eff_source_dist = self.zp.r_max/self.total_err # Technically, this should have a minus sign in front of it.
      
    def print_err_contribs(self):
        print 'Diffraction Limit: ' + "{:3.2f}".format(self.res_err/dzp.arcsec) + ' arcsec'
        print 'Defocus Error: ' + "{:3.2f}".format(self.defocus_err/dzp.arcsec) + ' arcsec'
        print 'Chromatic Aberration: ' + "{:3.2f}".format(self.chrom_err/dzp.arcsec) + ' arcsec'
        print 'Period Error: ' + "{:3.2f}".format(self.rad_err/dzp.arcsec) + ' arcsec'
        print 'Spherical Error: ' + "{:3.2f}".format(self.sph_err/dzp.arcsec) + ' arcsec'
        print 'Astigmatism/FC Error: ' + "{:3.2f}".format(self.astig_fc_err/dzp.arcsec) + ' arcsec'
        print 'Coma Error: ' + "{:3.2f}".format(self.coma_err/dzp.arcsec) + ' arcsec'

    __run_tol_assignment = run_tol_assignment
    __compute_diffration_limit_error = compute_diffraction_limit_error
    __compute_defocus_error = compute_defocus_error
    __compute_chromatic_aberration_error = compute_chromatic_aberration_error
    __compute_radial_error = compute_radial_error
    __compute_spherical_error = compute_spherical_error
    __compute_astigmatism_fc_error = compute_astigmatism_fc_error
    __compute_coma_error = compute_coma_error

def test_plot():
    # Units of distance are in mm, just like the zone definitions.
    f_op = 122000.
    max_radii = logspace(0,5,1001)

    zps = [dzp.zone_plate(f = f_op,r_max = max_radii[i],wave = dzp.al_k_wave,dwave = dzp.al_k_wave_width) for i in range(len(max_radii))]

    plt.ion()
    plt.figure(figsize = (10,10))
    
    tols_perfect = [tol_allocation(zps[i],0.0,0.0,0.0,0.0) for i in range(len(zps))]
    tols_defocus_105 = [tol_allocation(zps[i],f_op*10**-5,0.0,0.0,0.0) for i in range(len(zps))]
    tols_defocus_104 = [tol_allocation(zps[i],f_op*10**-4,0.0,0.0,0.0) for i in range(len(zps))] 
    tols_offaxis_05deg = [tol_allocation(zps[i],0.0,0.0,0.5*pi/180,0.0) for i in range(len(zps))]
    tols_offaxis_50deg = [tol_allocation(zps[i],0.0,0.0,5.0*pi/180,0.0) for i in range(len(zps))]
    
    
    # Plotting the limiting zone plate resolution.
    res = asarray([tols_perfect[i].res_err for i in range(len(max_radii))])/dzp.arcsec
    plt.loglog(max_radii,res,label = 'ZP Res.',color = 'r')
    # Plotting the chromatic aberration expected for the wavelength of interest.
    chrom = asarray([tols_perfect[i].chrom_err for i in range(len(max_radii))])/dzp.arcsec
    plt.loglog(max_radii,chrom,label = 'Chromatic',color = 'y')
    # Plotting the defocus expected for a fraction of the focal length.
    defocus_104 = asarray([tols_defocus_104[i].defocus_err for i in range(len(max_radii))])/dzp.arcsec
    defocus_105 = asarray([tols_defocus_105[i].defocus_err for i in range(len(max_radii))])/dzp.arcsec
    plt.loglog(max_radii,defocus_105,label = 'Def. 10^-5',color = 'g',linestyle = 'dashed')
    plt.loglog(max_radii,defocus_104,label = 'Def. 10^-4',color = 'g',linestyle = 'dotted')
    # Plotting the impact of astigmatism.
    astig_fc_05 = asarray([tols_offaxis_05deg[i].astig_fc_err for i in range(len(max_radii))])/dzp.arcsec
    astig_fc_50 = asarray([tols_offaxis_50deg[i].astig_fc_err for i in range(len(max_radii))])/dzp.arcsec
    plt.loglog(max_radii,astig_fc_05,label= 'Astig. 0.5 deg.',color = 'b',linestyle = 'dashed')
    plt.loglog(max_radii,astig_fc_50,label= 'Astig. 5.0 deg.',color = 'b',linestyle = 'dotted')
    # Plotting the impact of coma.
    coma_05 = asarray([tols_offaxis_05deg[i].coma_err for i in range(len(max_radii))])/dzp.arcsec
    coma_50 = asarray([tols_offaxis_50deg[i].coma_err for i in range(len(max_radii))])/dzp.arcsec
    plt.loglog(max_radii,coma_05,label = 'Coma 0.5 deg.',color = 'm',linestyle = 'dashed')
    plt.loglog(max_radii,coma_50,label = 'Coma 5.0 deg.',color = 'm',linestyle = 'dotted')
    # Plotting the impact of spherical aberration.
    sph = asarray([tols_perfect[i].sph_err for i in range(len(max_radii))])/dzp.arcsec
    plt.loglog(max_radii,sph,label = 'Sph. Ab.', color = 'k',linestyle = 'dotted')
    
    plt.xlabel('Zone Plate Radius (mm)')
    plt.ylabel('Angular Resolution (arcseconds)')
    plt.title('Angular Resolution Impacts on ZP Performance for\nFocal Length = ' + str(f_op) + ' mm ZPs')

    plt.xlim(1,70000.)
    plt.ylim(10**-4,100.)
    plt.legend()
    pdb.set_trace()
