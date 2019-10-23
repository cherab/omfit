import numpy as np
import matplotlib.pyplot as plt

def load_dms_output(config,world,plasma,spec,fibgeom):
    from raysect.optical.observer import FibreOptic, PowerPipeline0D, SpectralRadiancePipeline0D
    from raysect.optical import translate, rotate, rotate_basis
    from raysect.core import Vector3D, Point3D
    from raysect.core.ray import Ray as CoreRay

    power_arr = np.zeros(fibgeom.numfibres)
    spectra_arr = np.zeros((spec.pixels,fibgeom.numfibres))
    te_los = np.zeros((100,fibgeom.numfibres))
    ne_los = np.zeros((100,fibgeom.numfibres))
    d_los  = np.zeros((100,fibgeom.numfibres))
    if config['dms']['fibre_choice'] != -1:
        fib = [config['dms']['fibre_choice']]
    else:
        fib = np.linspace(1,fibgeom.numfibres,num=fibgeom.numfibres)
        
    for i, f in enumerate(fib):
        print ("Analysing fibre: ",int(f))
        fibgeom.set_fibre(number=int(f))        
        start_point    = Point3D(fibgeom.origin[0],fibgeom.origin[1],fibgeom.origin[2])
        forward_vector = Vector3D(fibgeom.xhat(),fibgeom.yhat(),fibgeom.zhat()).normalise()
        up_vector      = Vector3D(0, 0, 1.0)
        if config['dms']['power_pipeline']:
            power          = PowerPipeline0D()
            fibre  = FibreOptic([spectra], acceptance_angle=config['dms']['acceptance_angle'], radius=config['dms']['radius'], spectral_bins=spec.pixels, spectral_rays=config['dms']['spectral_rays'],
                                pixel_samples=config['dms']['pixel_samples'], transform=translate(*start_point)*rotate_basis(forward_vector, up_vector), parent=world)
            fibre.min_wavelength = spec.wlower
            fibre.max_wavelength = spec.wupper
            fibre.observe()
            power_arr[f] = power.value.mean
        else:
            power_arr[f] = None
        if config['dms']['radiance_pipeline']:
            spectra = SpectralRadiancePipeline0D(display_progress=False)
            fibre  = FibreOptic([spectra], acceptance_angle=config['dms']['acceptance_angle'], radius=config['dms']['radius'], spectral_bins=spec.pixels, spectral_rays=config['dms']['spectral_rays'],
                                pixel_samples=config['dms']['pixel_samples'], transform=translate(*start_point)*rotate_basis(forward_vector, up_vector), parent=world)
            fibre.min_wavelength = spec.wlower
            fibre.max_wavelength = spec.wupper
            fibre.observe()
            spectra_arr[:,f] = spectra.samples.mean
        else:
            spectra_arr[:,f] = None
        if config['dms']['los_profiles']:
            dist_var = np.linspace(0, fibgeom.fibre_distance_world(world), num=100, dtype=float)
            for j,t in enumerate(dist_var):
                x = start_point.x + fibgeom.xhat() * t
                y = start_point.y + fibgeom.yhat() * t
                z = start_point.z + fibgeom.zhat() * t
                te_los[j,f] = plasma.electron_distribution.effective_temperature(x,y,z)
                ne_los[j,f] = plasma.electron_distribution.density(x,y,z)
                d_los[j,f]  = t
        else:
            d_los[:,f]  = None
            te_los[:,f] = None
            ne_los[:,f] = None

    return power_arr,spectra_arr,te_los,ne_los,d_los

def load_dms_spectrometer(config):
    from cherab.mastu.div_spectrometer import spectrometer  
    spec=spectrometer()
    spec.set_range(setting=config['dms']['spectrometer']) 
    return spec	

def load_dms_fibres(config):
    from cherab.mastu.div_spectrometer import fibres
    fibgeom = fibres()
    fibgeom.set_bundle(group=config['dms']['fibres'])
    return fibgeom
	
