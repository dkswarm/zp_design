from numpy import *
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os
import pdb
from copy import deepcopy

# For heritage use of PyXFocus.
import PyXFocus.sources as src
import PyXFocus.transformations as trans
import PyXFocus.surfaces as surf
import PyXFocus.analyses as anal
import PyXFocus.conicsolve as conic

# For getting chromatic aberration implemented based on the actual line shape of Mg Kalpha.
import zp_design.emission_line_sources as emlines

# For the zone plates.
import zp_design.define_zp as dzp
import zp_design.zp_tolerancing as zpt

# Throwing in a zone plate so that we can test functionality.
n_rays = 10**5

#############################################################################################
# Defining the ray object class. This is identical to the ray class from arcusTrace.
class ray_obj:
    def __init__(self, PyXFocusRays, wave):
        self.opd = PyXFocusRays[0]
        self.x = PyXFocusRays[1]
        self.y = PyXFocusRays[2]
        self.z = PyXFocusRays[3]
        self.vx = PyXFocusRays[4]
        self.vy = PyXFocusRays[5]
        self.vz = PyXFocusRays[6]
        self.nx = PyXFocusRays[7]
        self.ny = PyXFocusRays[8]
        self.nz = PyXFocusRays[9]
        self.wave = wave
        self.index = arange(len(PyXFocusRays[0]))
        self.weight = ones(len(PyXFocusRays[0]))
            
    def set_prays(self,PyXFocusRays, ind = None):
        if ind is not None:
            self.opd[ind] = PyXFocusRays[0]
            self.x[ind] = PyXFocusRays[1]
            self.y[ind] = PyXFocusRays[2]
            self.z[ind] = PyXFocusRays[3]
            self.vx[ind] = PyXFocusRays[4]
            self.vy[ind] = PyXFocusRays[5]
            self.vz[ind] = PyXFocusRays[6]
            self.nx[ind] = PyXFocusRays[7]
            self.ny[ind] = PyXFocusRays[8]
            self.nz[ind] = PyXFocusRays[9]
        else:
            self.opd = PyXFocusRays[0]
            self.x = PyXFocusRays[1]
            self.y = PyXFocusRays[2]
            self.z = PyXFocusRays[3]
            self.vx = PyXFocusRays[4]
            self.vy = PyXFocusRays[5]
            self.vz = PyXFocusRays[6]
            self.nx = PyXFocusRays[7]
            self.ny = PyXFocusRays[8]
            self.nz = PyXFocusRays[9]
            
    def yield_prays(self, ind = None):
        if ind is not None:
            return [self.opd[ind],self.x[ind],self.y[ind],self.z[ind],self.vx[ind],self.vy[ind],self.vz[ind],self.nx[ind],self.ny[ind],self.nz[ind]]
        else:
            return [self.opd,self.x,self.y,self.z,self.vx,self.vy,self.vz,self.nx,self.ny,self.nz]
    
    def yield_object_indices(self, ind):
        new_object = deepcopy(self)
        for key in self.__dict__.keys():
            new_object.__dict__[key] = self.__dict__[key][ind]
        return new_object
    
    def pickle_me(self, pickle_file):
        new_object = copy.deepcopy(self)
        keys = new_object.__dict__.keys()
        attribs = []
        attribs.append(keys)
        attribs.append([new_object.__dict__[key] for key in keys])
        f = open(pickle_file,'wb')
        cPickle.dump(attribs,f)
        f.close()  

class wolterI:
    def __init__(self,r0,z0,z_ps = 2.5, z_ss = 2.5,mirror_length = 100.):
        # Radius of curvature of the Wolter-I optic.
        self.r0 = r0
        # Focal length of the Wolter-I optic.
        self.z0 = z0
        # Z separation between node point and primary.
        self.z_ps = z_ps
        # Z separation between node point and secondary.
        self.z_ss = z_ss
        # Z length of the mirror
        self.mirror_length = mirror_length
        # Distance in Z between the node point and the zone plate collimator. 
        # The first number means that the ZP is 50 mm away from the Wolter-I optic
        # nominally.
        self.zp_wolter_sep = 50. + mirror_length + z_ps
        self.primary_rmax = conic.primrad(self.z0 + self.z_ps + self.mirror_length,self.r0,self.z0)
        self.primary_rmin = conic.primrad(self.z0 + self.z_ps,self.r0,self.z0)

#############################################################################################
# Defining the functions needed to raytrace a zone plate. 

def illum_zp(xsource, ysource, zsource, rin, rout, tmin, tmax, num, source_size):
    '''
    Sub-apertured annulus beam confined to inner and outer radius 
    that converges to a location xsource, ysource, zsource.

    Parameters
    ----------
    xsource : int / float
        X-location of the source.
    ysource : int / float
        Y-location of the source.
    zsource : int / float
        Z-location of the source.
    rin : int / float
        Inner radius of sub-apertured annulus beam.
    rout : int / float
        Outer radius of sub-apertured annulus beam.
    tmin : int / float
        Minimum angular extent of sub-apertured annulus beam.
    tmax : int/ float
        Maximum angular extent of sub-apertured annulus beam.
    num : int
        Number of rays to create.
    source_size : int / float
        Radial dimension of the source.

    Returns
    -------
    rays : list
        List of ray parameters (opd, x, y, z, l, m, n, ux, uy, uz).

    '''
    rho = sqrt(rin**2+random.rand(num)*(rout**2-rin**2))
    theta = tmin + random.rand(num)*(tmax-tmin)
    x = rho*cos(theta)
    y = rho*sin(theta)
    z = zeros(num)

    rho_source, theta_source = random.rand(num)*source_size,random.rand(num)*2*pi
    vx = x - xsource + rho_source*cos(theta_source)
    vy = y - ysource + rho_source*sin(theta_source)
    vz = z - zsource
    
    norm = 1./sqrt(vx**2 + vy**2 + vz**2)
    vx,vy,vz = vx*norm,vy*norm,vz*norm 

    ux = repeat(0., num)
    uy = repeat(0., num)
    uz = repeat(0., num)
    opd = arange(num)

    return [opd, x, y, z, vx, vy, vz, ux, uy, uz]

def trace_zp(zp,src_dist,n_rays = 10**6,tols = None,order = 1,wave = emlines.mgk_limited,source_size = 0.0):
    '''
    Generalized raytracing function for a zone plate. We create a point source at the defined
    source distance, trace to the zone plate specified in zp, diffract based on order, and return the 
    rays immediately after they leave the zone plate. The zone plate is misaligned to the source based
    on the tolerances given in "tols".
    The zone plate is symmetrically placed about the +x axis.
    Parameters:
    -----------
    zp: zone plate object as defined from the zone plate design module.
    src_dist: distance between the zone plate and the source. This should be negative for forward propagation of the rays.
    n_rays: number of rays to create as propagating from the source
    tols: a list of the misalignment tolerances in the order: 
        [dz (defocus), offaxis_ang (theta in Menz dissertation), phi (mixing angle between misalignment in y vs. x), dd (error in groove density)].
        Note that dx and dy are resulted to the off-axis angle via offaxis_ang = arctan(sqrt(dx**2 + dy**2)/(source_dist + dz)) and phi = arcsin(dy/sqrt(dx**2 + dy**2))
    order: diffraction order
    wave: emission line instance from emlines to create the wavelength distribution of the source rays. If a single number, makes 
        the beam monochromatic.
    source_size: the size of the X-ray source. Note that due to implementation issue in PyXFocus, this doesn't work without a
        modified PyXFocus.sources package.
    '''
    if tols is not None:
        [dz,offaxis_ang,phi,dd] = tols
    else:
        [dz,offaxis_ang,phi,dd]= [0.,0.,0.,0.]

    # Calculating where to place the source for raytracing the ZP. These account for the alignment tolerances -- 
    # misplacement in z (defocus), off-axis angle and direction (offaxis_ang and phi).
    norm = -(src_dist - dz)*tan(offaxis_ang)
    xsource,ysource = cos(phi)*norm,sin(phi)*norm
    
    # Creates the rays already vignetted by the zone plate. A segmented zone plate is handled by having a nonzero 
    # r_min and a theta_seg less than 2*pi defined in the zp instance. 
    pyxrays = illum_zp(xsource,ysource,src_dist - dz,zp.r_min,zp.r_max,-zp.theta_seg/2.,zp.theta_seg/2.,n_rays,source_size)

    # We next draw on the wavelength distribution to give the rays each an individual wavelength 
    # and link the pyxrays to a ray object instance.
    if type(wave) is float:
        waves = ones(n_rays)*wave
    else:
        try:
            waves = wave.draw_waves(n_rays)
        except:
            raise TypeError('parameter "wave" is an incompatible format. Please specify a float or an emlines instance.' )

    # Preserving the idealized rays and creating the ray object for use with zp.trace_rays.
    ideal_pyxrays = deepcopy(pyxrays)
    rays = ray_obj(pyxrays,waves)

    # Finally, diffracting through the zone plate.
    zp.trace_rays(rays,order)
    return rays
    
def focus_zp_rays(rays):
    '''
    Quick handle function for diagnostics.
    '''
    pyxrays = rays.yield_prays()
    dz = anal.analyticImagePlane(pyxrays)
    trans.transform(pyxrays,0.,0.,dz,0.,0.,0.)
    surf.flat(pyxrays)
    return rays,pyxrays,dz

#############################################################################################
# This set of functions enables the tracing of a zone plate / Wolter-I optical system, as well 
# as facilitating a direct comparison between being directly illuminated by the beamline as 
# being illuminated by a collimating zone plate.

def wolterItrace(rays,wolterIoptic,scatter = False):
    '''
    A raytracing function taking a ray object centered on a Wolter-I optic defined 
    by the WolterI class and tracing through this optic handling vignetting. 

    Rays need to already be at the common Wolter-I focus.
    '''
    # First, getting the PyXFocus rays out of the ray object.    
    pyxrays = rays.yield_prays()

    # Find which rays intersect with primary mirror surface.
    surf.wolterprimary(pyxrays, wolterIoptic.r0, wolterIoptic.z0)
    mask1 = logical_and(pyxrays[3] < wolterIoptic.z_ps + wolterIoptic.mirror_length + wolterIoptic.z0 ,pyxrays[3] > wolterIoptic.z_ps + wolterIoptic.z0)

    # Filter out photons which would not strike primary surface and propagate good photons.
    pyxrays = [x[mask1] for x in pyxrays]
    new_waves = rays.wave[mask1]
    trans.reflect(pyxrays)

    # Creating a new ray object on the primary mirror.
    primary_pyxrays = deepcopy(pyxrays)
    primary_rays = ray_obj(primary_pyxrays,wave = new_waves)

    # Performing the same set of steps for the secondary mirror.
    surf.woltersecondary(pyxrays, wolterIoptic.r0, wolterIoptic.z0)
    mask2 = logical_and(pyxrays[3] > wolterIoptic.z0 - wolterIoptic.z_ss - wolterIoptic.mirror_length, wolterIoptic.z0, pyxrays[3] < wolterIoptic.z0 - wolterIoptic.z_ss)
    pyxrays = [x[mask2] for x in pyxrays]
    new_waves = new_waves[mask2]
    trans.reflect(pyxrays)

    # Creating a new ray object on the primary mirror.
    secondary_pyxrays = deepcopy(pyxrays)
    secondary_rays = ray_obj(secondary_pyxrays,wave = new_waves)

    # Add scatter.
    if scatter is True:
        pyxrays[4] = pyxrays[4] + random.normal(scale=1.5e-10, size=len(pyxrays[4]))
        pyxrays[5] = pyxrays[5] + random.normal(scale=1.5e-6, size=len(pyxrays[5]))
        pyxrays[6] = -sqrt(1.-pyxrays[5]**2-pyxrays[4]**2)

    # Propagating rays to this nominal focus location, computing the optimal focus location, 
    # and returning the computed rays.
    surf.flat(pyxrays)
    dz = anal.analyticImagePlane(pyxrays)
    pyxrays[3] = ones(len(pyxrays[3]))*dz   # Tracking the z-displacement relative to the nominal focus at zero.
    focused_rays = ray_obj(pyxrays,wave = new_waves)
    return primary_rays,secondary_rays,focused_rays,dz

def zp_before_wolter_trace(tol_budget,wolterIoptic,wave = emlines.mgk_limited,order = 1., zp_wolter_sep = 152.5):
    '''
    Calculates the radius of curvature resulting from a CZP with a given tolerance budget.
    '''
    tols = [tol_budget.dz,tol_budget.offaxis_ang,tol_budget.phi,tol_budget.dr]

    #############################
    # Tracing the zone plate first.
    rays = trace_zp(tol_budget.zp,-tol_budget.zp.f,tols = tols,order = order,wave = wave,source_size = tol_budget.source_size)
    pyxrays = rays.yield_prays()

    #############################
    # Handling the transforms necessary for the zone plate rays to interact with the WolterI optic correctly.

    # Next, translates the zone plate to be at the center of the ZP coordinates in XY, 
    # moves back to intersection point of the Wolter-I pair, and brings the rays there.
    trans.transform(pyxrays,0.,0.,zp_wolter_sep, 0.,0.,0.)    # mean(pyxrays[1]),mean(pyxrays[2])
    
    # Rotating the rays so that +z can point towards the mirrors and back towards the source.
    trans.transform(pyxrays,0.,0.,0.,0.,pi,pi/2)
    # Moving everything into the common coodinate system of the Wolter-I optic.
    trans.transform(pyxrays,0.,-wolterIoptic.r0,-wolterIoptic.z0,0.,0.,0.)

    # Putting rays on the zone plate into a new "zone_plate_rays" object so that everything is in the same coordinate system.
    zp_rays = deepcopy(pyxrays)
    zone_plate_rays = ray_obj(zp_rays,wave = rays.wave)

    ##############################
    # Next, adding an aperture mask to ensure that we're eliminating the zeroth order contribution. 
    # This happens immediately after (100 mm) the collimating zone plate.
    mask_czp_sep = 100.
    mask_tolerance = 0.5
    
    # Moving forward to CZP, then back by the mask separation.
    trans.transform(pyxrays,0.,0.,zp_wolter_sep + wolterIoptic.z0 - mask_czp_sep,0.,0.,0.)  
    surf.flat(pyxrays)

    # Calculating which rays make it through this mask.
    mask_condition = logical_and(sqrt(rays.x**2 + rays.y**2) < wolterIoptic.primary_rmax + mask_tolerance, sqrt(rays.x**2 + rays.y**2) > wolterIoptic.primary_rmin - mask_tolerance)
  
    # Now moving the rays back to the WolterI coordinates, and the rays should now be on the mask. 
    trans.transform(pyxrays,0.,0.,-zp_wolter_sep - wolterIoptic.z0 + mask_czp_sep,0.,0.,0.) 

    # And finally, applying the cut condition so that we get out only rays that pass through the mask.
    mask_rays = rays.yield_object_indices(ind = mask_condition)

    #if sum(mask_condition)
    ## Passing to the WolterI raytracing function, which handles everything else.
    #primary_rays,secondary_rays,focused_rays,dz = wolterItrace(mask_rays,wolterIoptic)

    return zone_plate_rays, mask_rays  # , primary_rays, secondary_rays, focused_rays, dz

def compare_beamline_vs_zp(tol_budget,wolterIoptic,wave = emlines.mgk_limited,source_size = 0.0):
    # Doing the zp trace with all the other functions that we built.
    zp_rays, zp_wolter_rays, zp_dz = zp_wolter_system_trace(tol_budget, wolterIoptic, wave = wave, source_size = source_size)

    # Doing the beamline-only trace by tracing to a ZP but setting the response to be that of zeroth order. 
    beamline_null_tols = zpt.tol_allocation(zp = tol_budget.zp, dz = 0., dr = 0.0, offaxis_ang = 0.0, phi = 0.0)
    beamline_rays,beamline_wolter_rays,beamline_zp = zp_wolter_system_trace(beamline_null_tols, wolterIoptic, order = 0.0, wave = wave, source_size = source_size)
    
    return zp_rays,beamline_rays,zp_wolter_rays,beamline_wolter_rays,zp_dz,beamline_zp

def phi_sample(zp,wolterIoptic,dzs,offaxis_angs,phis = linspace(0,9./5*pi,10)):
    zp_grid = zeros((len(dzs),len(offaxis_angs),len(phis)))
    beam_grid = zeros((len(dzs),len(offaxis_angs),len(phis)))
   
    for i in range(len(dzs)):
        for j in range(len(offaxis_angs)):
            for k in range(len(phis)):
                czp_tols = zpt.tol_allocation(zp,dz = dzs[i],dr = 0.0,offaxis_ang = offaxis_angs[j],phi = phis[k])
                zp_rays,beamline_rays,zp_wolter_rays,beamline_wolter_rays,zp_dz,beamline_zp = compare_beamline_vs_zp(czp_tols,wolterIoptic,wave = emlines.mgk_limited,source_size = 0.060)
                zp_grid[i,j,k] = zp_dz
                beam_grid[i,j,k] = beamline_zp
                print(dzs[i],offaxis_angs[j],phis[k])
    return zp_grid,beam_grid
