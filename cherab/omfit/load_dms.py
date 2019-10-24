import numpy as np
import matplotlib.pyplot as plt

def load_dms_output(config,world,plasma,spec,fibgeom,numlos):
    from cherab.core.atomic.elements import hydrogen, deuterium, carbon, helium, nitrogen, neon, argon, krypton, xenon
    from raysect.optical.observer import FibreOptic, PowerPipeline0D, SpectralRadiancePipeline0D
    from raysect.optical import translate, rotate, rotate_basis
    from raysect.core import Vector3D, Point3D
    from raysect.core.ray import Ray as CoreRay

    power_arr = np.zeros(fibgeom.numfibres)
    spectra_arr = np.zeros((spec.pixels,fibgeom.numfibres))
    te_los = np.zeros((numlos,fibgeom.numfibres))
    ne_los = np.zeros((numlos,fibgeom.numfibres))
    d_los  = np.zeros((numlos,fibgeom.numfibres))
    nN_los = np.zeros((numlos,fibgeom.numfibres,8))
    if config['dms']['fibre_choice'] != -1:
        fib = [config['dms']['fibre_choice']]
    else:
        fib = np.linspace(1,fibgeom.numfibres,num=fibgeom.numfibres)
        
    for i, f in enumerate(fib):
        arridx = int(f)-1
        print ("Analysing fibre: ",int(f))
        fibgeom.set_fibre(number=int(f))        
        start_point    = Point3D(fibgeom.origin[0],fibgeom.origin[1],fibgeom.origin[2])
        forward_vector = Vector3D(fibgeom.xhat(),fibgeom.yhat(),fibgeom.zhat()).normalise()
        up_vector      = Vector3D(1.0, 1.0, 1.0)
        if config['dms']['power_pipeline']:
            power          = PowerPipeline0D()
            fibre  = FibreOptic([spectra], acceptance_angle=config['dms']['acceptance_angle'], radius=config['dms']['radius'], spectral_bins=spec.pixels, spectral_rays=config['dms']['spectral_rays'],
                                pixel_samples=config['dms']['pixel_samples'], transform=translate(*start_point)*rotate_basis(forward_vector, up_vector), parent=world)
            fibre.min_wavelength = spec.wlower
            fibre.max_wavelength = spec.wupper
            fibre.observe()
            power_arr[arridx] = power.value.mean

        if config['dms']['radiance_pipeline']:
            spectra = SpectralRadiancePipeline0D(display_progress=False)
            fibre  = FibreOptic([spectra], acceptance_angle=config['dms']['acceptance_angle'], radius=config['dms']['radius'], spectral_bins=spec.pixels, spectral_rays=config['dms']['spectral_rays'],
                                pixel_samples=config['dms']['pixel_samples'], transform=translate(*start_point)*rotate_basis(forward_vector, up_vector), parent=world)
            fibre.min_wavelength = spec.wlower
            fibre.max_wavelength = spec.wupper
            fibre.observe()
            spectra_arr[:,arridx] = spectra.samples.mean

        if config['dms']['los_profiles']:
            n0 = plasma.composition.get(nitrogen, 0)
            n1 = plasma.composition.get(nitrogen, 1)
            n2 = plasma.composition.get(nitrogen, 2)
            n3 = plasma.composition.get(nitrogen, 3)
            n4 = plasma.composition.get(nitrogen, 4)
            n5 = plasma.composition.get(nitrogen, 5)
            n6 = plasma.composition.get(nitrogen, 6)
            n7 = plasma.composition.get(nitrogen, 7)
            dist_var = np.linspace(0, fibgeom.fibre_distance_world(world), num=numlos, dtype=float)
            for j,t in enumerate(dist_var):
                x = start_point.x + fibgeom.xhat() * t
                y = start_point.y + fibgeom.yhat() * t
                z = start_point.z + fibgeom.zhat() * t
                te_los[j,arridx] = plasma.electron_distribution.effective_temperature(x,y,z)
                ne_los[j,arridx] = plasma.electron_distribution.density(x,y,z)
                nN_los[j,arridx,0] = n0.distribution.density(x, y, z)
                nN_los[j,arridx,1] = n1.distribution.density(x, y, z)
                nN_los[j,arridx,2] = n2.distribution.density(x, y, z)
                nN_los[j,arridx,3] = n3.distribution.density(x, y, z)
                nN_los[j,arridx,4] = n4.distribution.density(x, y, z)
                nN_los[j,arridx,5] = n5.distribution.density(x, y, z)
                nN_los[j,arridx,6] = n6.distribution.density(x, y, z)
                nN_los[j,arridx,7] = n7.distribution.density(x, y, z)
                d_los[j,arridx]  = t

    return power_arr,spectra_arr,te_los,ne_los,nN_los,d_los

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
	
