import json
import argparse
import matplotlib.pyplot as plt
import numpy as np
import time
from raysect.optical import World
from cherab.omfit import load_machine
from netCDF4 import Dataset
from collections import namedtuple
vertex = namedtuple("vertex", "x y")

def clamp(n, minn, maxn):
    if n < minn:
        return minn
    elif n > maxn:
        return maxn
    else:
        return n

class dms:
    """
    Divertor monitoring class
    """

    def __init__(self,world=None,config=None,plasma=None):
        if config is not None:
            self.world = world
            self.config = config
            self.plasma = plasma
            self.Nlos=100
            self.power = None
            self.spectra=None
            self.te_los=None
            self.ne_los=None
            self.ni_los=None
            self.nz_los=None
            self.d_los=None
            self.spec = None
            self.fibres=None
            
    def simulate(self):
        from cherab.omfit import load_dms_output,load_dms_fibres,load_dms_spectrometer        
        from raysect.core import Vector3D, Point3D
        # Load diagnostic geometry
        self.fibres = load_dms_fibres(self.config)
        # Load diagnostic settings
        self.spec = load_dms_spectrometer(self.config)
        # Simulate synthetic measurement
        self.power,self.spectra,self.te_los,self.ne_los,self.ni_los,self.nz_los,self.d_los = load_dms_output(self.config, self.world, self.plasma, self.spec, self.fibres, self.Nlos)        

    def write_cdf(self,ncfile='cherab.nc'):            
        # Output netCDF file
        dataset  = Dataset(ncfile,'a')
        # Geometry details
        dmsgroup = dataset.createGroup("DMS")

        nFibres = dmsgroup.createDimension('nFibres',self.fibres.numfibres)
        Fibres  = dmsgroup.createVariable('nFibres',np.float32,('nFibres'))

        Wavelength     = dmsgroup.createDimension('Wavelength',self.spec.pixels)
        Wavelength_arr = dmsgroup.createVariable('Wavelength',np.float32,('Wavelength'))

        nVec = dmsgroup.createDimension('nVec',3)

        uVec       = dmsgroup.createVariable('uVec',np.float32,('nVec','nFibres'))
        uVec.label = 'Fibre unit vectors'
        uVec.units = 'arb'

        distance = dmsgroup.createVariable('distance',np.float32,('nFibres'))
        distance.label = 'Line-of-sight distance'
        distance.units = 'm'
   
        orig       = dmsgroup.createVariable('origin',np.float32,('nVec','nFibres'))
        orig.label = 'Fibre origin'
        orig.units = 'm'

        power       = dmsgroup.createVariable('line-integrated_power',np.float32,('nFibres'))
        power.label = 'Line-integrated power'
        power.units = 'W'
        
        spectra       = dmsgroup.createVariable('line-integrated_spectrum',np.float32,('Wavelength','nFibres'))
        spectra.label = 'Line-integrated spectrum'
        spectra.units = 'W/m2'

        Fibres[:]          = range(1,self.fibres.numfibres+1)
        power[:]           = self.power
        Wavelength_arr[:]  = self.spec.wlngth

        LoS          = dmsgroup.createDimension('LoS',self.Nlos)

        te_los       = dmsgroup.createVariable('LoS_Te',np.float32,('LoS','nFibres'))
        te_los.label = 'Line-of-sight Te'
        te_los.units = 'eV'

        ne_los       = dmsgroup.createVariable('LoS_ne',np.float32,('LoS','nFibres'))
        ne_los.label = 'Line-of-sight ne'
        ne_los.units = 'm-3'

        Nions        = dmsgroup.createDimension('Nions',self.nz_los.shape[2])
        nz_los       = dmsgroup.createVariable('LoS_nz',np.float32,('LoS','nFibres','Nions'))
        nz_los.label = 'Line-of-sight nz'
        nz_los.units = 'm-3'

        Nd        = dmsgroup.createDimension('Nd',2)
        ni_los       = dmsgroup.createVariable('LoS_ni',np.float32,('LoS','nFibres','Nd'))
        ni_los.label = 'Line-of-sight ni'
        ni_los.units = 'm-3'

        d_los       = dmsgroup.createVariable('LoS_dist',np.float32,('LoS','nFibres'))
        d_los.label = 'Line-of-sight distance'
        d_los.units = 'm'

        for i in range(self.fibres.numfibres):
            self.fibres.set_fibre(number=int(i+1))        
            orig[:,i]    = [self.fibres.origin[0],self.fibres.origin[1],self.fibres.origin[2]]
            uVec[:,i]    = [self.fibres.xhat(),self.fibres.yhat(),self.fibres.zhat()]
            distance[i]  = self.fibres.fibre_distance_world(self.world)
            spectra[:,i] = self.spectra[:,i]
            te_los[:,i]  = self.te_los[:,i]
            ne_los[:,i]  = self.ne_los[:,i]
            ni_los[:,i,:]= self.ni_los[:,i,:]
            nz_los[:,i,:]= self.nz_los[:,i,:]
            d_los[:,i]   = self.d_los[:,i]
           
        dataset.close()
class camera:
    """
    Filtered camera details
    """
    def __init__(self,world=None,config=None,plasma=None):
        if config is not None:
            self.world  = world
            self.config = config
            self.plasma = plasma
            self.camera = None
    def simulate(self):
        from cherab.omfit import load_camera
        self.camera = load_camera(self.config, self.world)
        plt.ion()
        self.camera.observe()
        plt.ioff()
        plt.show()
        
    def write_cdf(self,ncfile='cherab.nc'):            
        # Output netCDF file
        dataset = Dataset(ncfile,'a')
        # Geometry details
        camgroup = dataset.createGroup("camera")
        xPixels = camgroup.createDimension('xPixels',self.camera.pipelines[0].frame.mean.shape[0])
        yPixels = camgroup.createDimension('yPixels',self.camera.pipelines[0].frame.mean.shape[1])
        image = camgroup.createVariable('image',np.float32,('xPixels','yPixels'))
        for i in range(self.camera.pipelines[0].frame.mean.shape[0]):      
            image[i,:] = self.camera.pipelines[0].frame.mean[:,i]

        dataset.close()
class simulation:
    """
    Plasma simulation details
    """
    def __init__(self,world=None,config=None):
        self.plasma     = None
        self.mesh       = None
        self.world      = world
        self.config     = config
        self.num        = 500
        Rrange          = np.array(config['plasma']['edge']['Rrange'])
        Zrange          = np.array(config['plasma']['edge']['Zrange'])
        self.xrange     = np.linspace(Rrange[0],Rrange[1], self.num)
        self.yrange     = np.linspace(Zrange[0],Zrange[1], self.num) 

    def load(self):
        if self.config['plasma']['edge']['present']:
            from cherab.omfit import load_emission, load_edge_simulation, load_profiles
            # Load the edge plasma solution if present	  
            self.plasma, self.mesh = load_edge_simulation(self.config, self.world)
            # Load all the emission line models and emission profiles
            self.em_plasma = load_emission(self.config, self.plasma, self.xrange, self.yrange, self.num)
            # Get profile data
            self.te_plasma, self.ne_plasma, self.ni_plasma, self.nz_plasma = load_profiles(self.config, self.plasma, self.xrange, self.yrange, self.num)

    def write_cdf(self,ncfile='cherab.nc'):            
        # Output netCDF file
        dataset = Dataset(ncfile,'a')
        # Mesh details
        if self.config['plasma']['edge']['mesh']:
            meshgroup       = dataset.createGroup("mesh")
            nDistributionX  = meshgroup.createDimension('nDistributionX',4)
            nDistributionsX = meshgroup.createVariable('nDistributionX',np.float32,('nDistributionX'))
            nDistributionY  = meshgroup.createDimension('nDistributionY',len(self.mesh.triangles))
            nDistributionsY = meshgroup.createVariable('nDistributionY',np.float32,('nDistributionY'))

            nDistributionsX[:]    = (0,1,2,3)
            nDistributionsX.units = '[-]'
            nDistributionsX.label = 'Vertice points'
            nDistributionsY[:]    = range(0,len(self.mesh.triangles))
            nDistributionsY.units = '[-]'
            nDistributionsY.label = 'Number of triangles'
            
            mesh_coords_x     = meshgroup.createVariable('mesh_coords_x',np.float32,('nDistributionX','nDistributionY'))
            mesh_coords_x .label = 'Mesh X vertice coordinates'
            mesh_coords_x .units = 'm'

            mesh_coords_y     = meshgroup.createVariable('mesh_coords_y',np.float32,('nDistributionX','nDistributionY'))
            mesh_coords_y.label = 'Mesh Y vertice coordinates'
            mesh_coords_y.units = 'm'

            i = 0 
            for triangle in self.mesh.triangles:
                i1, i2, i3 = triangle
                v1 = vertex(*self.mesh.vertex_coords[i1])
                v2 = vertex(*self.mesh.vertex_coords[i2])
                v3 = vertex(*self.mesh.vertex_coords[i3])
                mesh_coords_x[:,i] = [v1.x, v2.x, v3.x, v1.x]
                mesh_coords_y[:,i] = [v1.y, v2.y, v3.y, v1.y]
                i += 1
     
        # Plasma details
        plasmagroup     = dataset.createGroup("plasma")
        nDistributionX  = plasmagroup.createDimension('nDistributionX',self.num)
        nDistributionsX = plasmagroup.createVariable('nDistributionX',np.float32,('nDistributionX'))

        nDistributionY  = plasmagroup.createDimension('nDistributionY',self.num)
        nDistributionsY = plasmagroup.createVariable('nDistributionY',np.float32,('nDistributionY'))

        nDistributionsX[:] = self.xrange
        nDistributionsX.units = 'm'
        nDistributionsX.label = 'Plasma length'
        nDistributionsY[:] = self.yrange
        nDistributionsY.units = 'm'
        nDistributionsY.label = 'Plasma height'

        Te = plasmagroup.createVariable('Te',np.float32,('nDistributionX','nDistributionY'))
        Te.label='Plasma Te'
        Te.units='eV'            
        for i, x in enumerate(self.xrange):
            Te[i,:]=self.te_plasma[i,:]

        ne = plasmagroup.createVariable('ne',np.float32,('nDistributionX','nDistributionY'))
        ne.label='Plasma ne'
        ne.units='m-3'
        for i, x in enumerate(self.xrange):
            ne[i,:]=self.ne_plasma[i,:]

        ND = plasmagroup.createDimension('ND',2)
        ni = plasmagroup.createVariable('ni',np.float32,('nDistributionX','nDistributionY','ND'))
        ni.label='Plasma ni'
        ni.units='m-3'
        for i, x in enumerate(self.xrange):
            ni[i,:,:]=self.ni_plasma[i,:,:]

        emiss = plasmagroup.createVariable('emiss',np.float32,('nDistributionX','nDistributionY'))
        emiss.label='Plasma emission'
        emiss.units='ph/m-3/s'
        for i, x in enumerate(self.xrange):
                emiss[i,:]=self.em_plasma[i,:]

        Nions= plasmagroup.createDimension('Nions',self.nz_plasma.shape[2])
        nz = plasmagroup.createVariable('nz',np.float32,('nDistributionX','nDistributionY','Nions'))
        nz.label='Plasma impurity density'
        nz.units='m-3'
        for i, x in enumerate(self.xrange):
            nz[i,:,:]=self.nz_plasma[i,:,:]

        dataset.close()
  
        
if __name__ == '__main__':
       
    parser = argparse.ArgumentParser('Simulate a spectrum')
    parser.add_argument('configfile', type=str, help="A JSON format CHERAB configuration file.")
    args = parser.parse_args()
    # read/parse config file
    with open(args.configfile, 'r') as f:
        config = json.load(f)

    # Setup netCDF output
    ncfile = 'cherab.nc'
    dataset = Dataset(ncfile,'w')
    dataset.history = 'Created' + time.ctime(time.time())
    dataset.close()

    # Load the machine
    world = World()
    load_machine(config, world)

    # Load the plasma simulation
    sim = simulation(world=world,config=config)
    sim.load()
    sim.write_cdf()
    
    # Load each diagnostic
    if config['dms']['simulate']:
        diagDMS = dms(world=world,config=config,plasma=sim.plasma)
        diagDMS.simulate()
        diagDMS.write_cdf(ncfile=ncfile)

    if config['observer']['simulate']:
        diagCAM = camera(world=world,config=config,plasma=sim.plasma)
        diagCAM.simulate()
        diagCAM.write_cdf(ncfile=ncfile)

     
