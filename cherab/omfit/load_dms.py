import numpy as np
import matplotlib.pyplot as plt
from cherab.core.atomic.elements import hydrogen, deuterium, carbon, helium, nitrogen, neon, argon, krypton, xenon
_SPECIES_LOOKUP = {
    "hydrogen": hydrogen,
    "deuterium": deuterium,
    "carbon": carbon,
    'helium': helium,
    'nitrogen': nitrogen,
    'neon': neon,
    'argon': argon,
    'krypton': krypton,
    'xenon': xenon,
}

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
    ni_los = np.zeros((numlos,fibgeom.numfibres,2))
    d_los  = np.zeros((numlos,fibgeom.numfibres))
    if config['plasma']['edge']['nz2D']:
        try:
            species = _SPECIES_LOOKUP[config["plasma"]['edge']["nz2D_species"]]
        except KeyError:
            raise ValueError("The species emission specification for this config file is invalid.")
        nz_los = np.zeros((numlos,fibgeom.numfibres,species.atomic_number+1))
        ions = np.linspace(0,species.atomic_number,num=species.atomic_number+1)  
    else:
        nz_los=np.zeros((numlos,fibgeom.numfibres,1))

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
            fibre  = FibreOptic([power], acceptance_angle=config['dms']['acceptance_angle'], radius=config['dms']['radius'], spectral_bins=spec.pixels, spectral_rays=config['dms']['spectral_rays'],
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
            dist_var = np.linspace(0, fibgeom.fibre_distance_world(world), num=numlos, dtype=float)
            for j,t in enumerate(dist_var):
                x = start_point.x + fibgeom.xhat() * t
                y = start_point.y + fibgeom.yhat() * t
                z = start_point.z + fibgeom.zhat() * t
                if config['plasma']['edge']['nz2D']:
                    for k,zeff in enumerate(ions):
                        composition = plasma.composition.get(species, k)
                        nz_los[j, arridx, k]  = composition.distribution.density(x, y, z)
                if config['plasma']['edge']['ni2D']:
                    for k in range(2):
                        composition = plasma.composition.get(deuterium, k)
                        ni_los[j, arridx, k]  = composition.distribution.density(x, y, z)
                if config['plasma']['edge']['Te2D']:
                    te_los[j,arridx] = plasma.electron_distribution.effective_temperature(x,y,z)
                if config['plasma']['edge']['ne2D']:
                    ne_los[j,arridx] = plasma.electron_distribution.density(x,y,z)
                d_los[j,arridx]  = t

    return power_arr,spectra_arr,te_los,ne_los,ni_los,nz_los,d_los

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
	
