import json
import argparse
import matplotlib.pyplot as plt
import numpy as np
import time
from raysect.optical import World
from cherab.omfit import load_machine, load_edge_simulation
from netCDF4 import Dataset

class dms:
    """
    Divertor monitoring class
    """

    def __init__(self,config=None):
        if config is not None:
            self.world = World()
            self.config = config
            # Load all the CAD file geometry for the selected machine
            load_machine(config, self.world)

    def simulate(self):
        from cherab.omfit import load_dms_power,load_dms_fibres,load_dms_spectrometer        
        # Load diagnostic geometry
        self.fibres = load_dms_fibres(self.config)
        # Load diagnostic geometry
        self.spec = load_dms_spectrometer(self.config)
        if self.config['plasma']['edge']['present'] or self.config['plasma']['core']['present']:
            # Load the edge plasma solution if present	  
            plasma = load_edge_simulation(self.config, self.world)
            # Load the specified filtered camera
            synth = load_dms_power(self.config, self.world, self.plasma,self.spec,self.fibgeom)
        else:
            print ("No SOL or core model specified")
    def write_cdf(self,ncfile='cherab.nc'):            
        # Output netCDF file
        dataset = Dataset(ncfile,'a')
        # Geometry details
        dmsgroup = dataset.createGroup("DMS")
        nFibres = dmsgroup.createDimension('nFibres',self.fibres.numfibres)
        nVec = dmsgroup.createDimension('nVec',3)
        uVec = dmsgroup.createVariable('uVec',np.float32,('nVec','nFibres'))
        uVec.label = 'Fibre unit vectors'
        uVec.units = 'arb'
        distance = dmsgroup.createVariable('distance',np.float32,('nFibres'))
        distance.label = 'Line-of-sight distance'
        distance.units = 'm'
        for i in range(self.fibres.numfibres):
            self.fibres.set_fibre(number=int(i+1))        
            uVec[:,i] = [self.fibres.xhat(),self.fibres.yhat(),self.fibres.zhat()]
            distance[i]=self.fibres.fibre_distance_world(self.world)
    
        orig = dmsgroup.createVariable('origin',np.float32,('nVec'))
        orig.label = 'Fibre origin'
        orig.units = 'm'
        orig[:]    = [self.fibres.origin[0],self.fibres.origin[1],self.fibres.origin[2]]
        dataset.close()
class camera:
    """
    Filtered camera
    """
    def __init__(self,config=None):
        if config is not None:
            self.world = World()
            self.config = config
            # Load all the CAD file geometry for the selected machine
            load_machine(config, self.world)

    def simulate(self):
        from cherab.omfit import load_camera, load_emission

        if self.config['plasma']['edge']['present'] or self.config['plasma']['core']['present']:
            # Load the edge plasma solution if present	  
            plasma = load_edge_simulation(self.config, self.world)
            # Load all the emission line models
            load_emission(self.config, self.plasma)
            # Load the specified filtered camera
        else:
            print ("No SOL or core model specified")

        self.camera = load_camera(self.config, self.world)
#        plt.ion()
        self.camera.observe()
#        plt.ioff()
#        plt.show()
    def write_cdf(self,ncfile='cherab.nc'):            
        # Output netCDF file
        print ("Not developed yet")        
        
if __name__ == '__main__':
       
    parser = argparse.ArgumentParser('Simulate a spectrum')
    parser.add_argument('configfile', type=str, help="A JSON format CHERAB configuration file.")
    args = parser.parse_args()
    # read/parse config file
    with open(args.configfile, 'r') as f:
        config = json.load(f)

    ncfile = 'cherab.nc'
    dataset = Dataset(ncfile,'w')
    dataset.history = 'Created' + time.ctime(time.time())
    dataset.close()
    
    if config['dms']['simulate']:
        sim = dms(config)
        sim.simulate()
        sim.write_cdf(ncfile=ncfile)

    if config['observer']['simulate']:
        sim = camera(config)
        sim.simulate()
        sim.write_cdf(ncfile=ncfile)


     
