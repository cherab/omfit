from cherab.openadas import OpenADAS
from raysect.optical.observer import FibreOptic, PowerPipeline0D
from raysect.core import Vector3D, Point3D
from cherab.core.atomic import Line
from cherab.core.model import ExcitationLine, RecombinationLine, Bremsstrahlung
from cherab.core.model.lineshape import StarkBroadenedLine

def load_dms_power(config,world,plasma,spec,fibgeom):

    plasma.atomic_data = OpenADAS(permit_extrapolation=True)
    plasma.models= [d_alpha_excit, d_alpha_recom]

    power_arr = np.zeros(fibgeom.numfibres)

    d_alpha         = Line(deuterium, 0, (3, 2))
    d_alpha_excit   = ExcitationLine(d_alpha, lineshape=StarkBroadenedLine)
    d_alpha_recom   = RecombinationLine(d_alpha, lineshape=StarkBroadenedLine)

    for i, f in enumerate(power_arr):
        print ("Analysing fibre: ",int(i+1))
        fibgeom.set_fibre(number=int(i+1))        
        start_point    = Point3D(fibgeom.origin[0],fibgeom.origin[1],fibgeom.origin[2])
        forward_vector = Vector3D(fibgeom.xhat(),fibgeom.yhat(),fibgeom.zhat()).normalise()
        up_vector      = Vector3D(0, 0, 1.0)
        power          = PowerPipeline0D()
        fibre  = FibreOptic([power], acceptance_angle=1, radius=0.001, spectral_bins=spec.pixels, spectral_rays=1,
                            pixel_samples=5, transform=translate(*start_point)*rotate_basis(forward_vector, up_vector), parent=world)
        fibre.min_wavelength = spec.wlower
        fibre.max_wavelength = spec.wupper
        fibre.observe()
        power_arr[i] = power.value.mean

    return power_arr

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
	
