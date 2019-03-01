import json
import argparse
import matplotlib.pyplot as plt
import numpy as np
import time
from raysect.optical import World
from cherab.omfit import load_machine
from netCDF4 import Dataset

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
            self.power = None
            self.spectra=None
            self.te_plasma=None
            self.ne_plasma=None
            self.spec = None
            self.fibres=None
            self.plasma = plasma
    def simulate(self):
        from cherab.omfit import load_dms_output,load_dms_fibres,load_dms_spectrometer        
        # Load diagnostic geometry
        self.fibres = load_dms_fibres(self.config)
        # Load diagnostic settings
        self.spec = load_dms_spectrometer(self.config)
        # Simulate synthetic measurement
        self.power,self.spectra = load_dms_output(self.config, self.world, self.plasma, self.spec, self.fibres)

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
   
        orig       = dmsgroup.createVariable('origin',np.float32,('nVec'))
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
        orig[:]            = [self.fibres.origin[0],self.fibres.origin[1],self.fibres.origin[2]]

        for i in range(self.fibres.numfibres):
            self.fibres.set_fibre(number=int(i+1))        
            uVec[:,i]   = [self.fibres.xhat(),self.fibres.yhat(),self.fibres.zhat()]
            distance[i] = self.fibres.fibre_distance_world(self.world)
            spectra[:,i]= self.spectra[:,i]

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
    Filtered camera details
    """
    def __init__(self,world=None,config=None):
        self.plasma    = None
        self.world     = world
        self.config    = config
        self.te_plasma = np.zeros((500, 500))
        self.ne_plasma = np.zeros((500, 500))
        self.xrange    = np.linspace(0.0, 2.5, 500)
        self.yrange    = np.linspace(-2.5, 2.5, 500) 
    def load(self):
        if self.config['plasma']['edge']['present']:
            from cherab.omfit import load_emission, load_edge_simulation
            # Load the edge plasma solution if present	  
            self.plasma = load_edge_simulation(self.config, self.world)
            # Load all the emission line models
            load_emission(self.config, self.plasma)

            for i, x in enumerate(self.xrange):
                for j, y in enumerate(self.yrange):
                    val = self.plasma.electron_distribution.effective_temperature(x, 0.0, y)
                    self.te_plasma[j, i] = clamp(val,0,50.0)
                    val = self.plasma.electron_distribution.density(x, 0.0, y)
                    self.ne_plasma[j, i] = val

    def write_cdf(self,ncfile='cherab.nc'):            
        # Output netCDF file
        dataset = Dataset(ncfile,'a')
        # Geometry details
        plasmagroup = dataset.createGroup("plasma")

        nDistributionX = plasmagroup.createDimension('nDistributionX',500)
        nDistributionsX = plasmagroup.createVariable('nDistributionX',np.float32,('nDistributionX'))

        nDistributionY = plasmagroup.createDimension('nDistributionY',500)
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

        Ne = plasmagroup.createVariable('Ne',np.float32,('nDistributionX','nDistributionY'))
        Ne.label='Plasma ne'
        Ne.units='m-3'

        if self.xrange is not None:
            for i, x in enumerate(self.xrange):
                    Te[i,:]=self.te_plasma[i,:]
                    Ne[i,:]=self.ne_plasma[i,:]
        else:
            Te[:,:] = self.te_plasma
            Ne[:,:] = self.ne_plasma
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
        diagDMS = dms(world=world,config=config,plasma=sim)
        diagDMS.simulate()
        diagDMS.write_cdf(ncfile=ncfile)

    if config['observer']['simulate']:
        diagCAM = camera(world=world,config=config,plasma=sim)
        diagCAM.simulate()
        diagCAM.write_cdf(ncfile=ncfile)


     
