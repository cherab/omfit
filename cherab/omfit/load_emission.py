
from cherab.core.atomic import Line
from cherab.core.atomic.elements import hydrogen, deuterium, carbon, helium, nitrogen, neon, argon, krypton, xenon
from cherab.core.model import ExcitationLine, RecombinationLine, MultipletLineShape, StarkBroadenedLine, Bremsstrahlung
from raysect.core import Vector3D, Point3D
from raysect.optical import Spectrum
import numpy as np
from cherab.openadas import OpenADAS
from cherab.openadas.install import install_adf15
from cherab.openadas.repository import add_wavelength

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

_EMISSION_TYPE_LOOKUP = {
    "ExcitationLine": ExcitationLine,
    "RecombinationLine": RecombinationLine
}


def load_emission(config, plasma, n1_density,n1_emiss, xr, yr):

    models = []

    plasma.atomic_data = OpenADAS(permit_extrapolation=True)
    if config["plasma"]["installADF15"]:
        for adf15 in config["plasma"]["adf15"]:
            try:
                species = _SPECIES_LOOKUP[adf15["species"]]
            except KeyError:
                raise ValueError("The species emission specification for this config file is invalid.")                                   
            charge     = adf15["ionisation"]
            file_path  = adf15["file_path"]
            adas_path  = adf15["adas_path"]
            install_adf15(species,charge,file_path,adas_path=adas_path)
            
    if config["plasma"]["bremsstrahlung"]:
        models.append(Bremsstrahlung())

    for emission_instruction in config["plasma"]["emission_instructions"]:
       
        try:
            species = _SPECIES_LOOKUP[emission_instruction["species"]]
        except KeyError:
            raise ValueError("The species emission specification for this config file is invalid.")                                   
        ionisation = emission_instruction["ionisation"]
        upper      = emission_instruction["upper"]
        lower      = emission_instruction["lower"]
        transition =(upper,lower)
        wavelength = emission_instruction["wavelength"]
        if not wavelength == 0: 
            add_wavelength(species,ionisation,transition,wavelength)

        # Do a test call to see if the species exists on this plasma.
        _ = plasma.composition.get(species, ionisation)

        line = Line(species, ionisation, (upper, lower))
        if emission_instruction["multiplet"]:
            multipletWvlngths = emission_instruction["multipletWvlngths"]
            multipletRatios = emission_instruction["multipletRatios"]
            multiplet = [multipletWvlngths,multipletRatios]    
            models.append(_EMISSION_TYPE_LOOKUP[emission_instruction["type"]](line,lineshape=MultipletLineShape,lineshape_args=[multiplet]))
        elif emission_instruction['stark']:
            models.append(_EMISSION_TYPE_LOOKUP[emission_instruction["type"]](line,lineshape=StarkBroadenedLine))
        else:
            models.append(_EMISSION_TYPE_LOOKUP[emission_instruction["type"]](line))
    species  = _SPECIES_LOOKUP["nitrogen"]
    NII      = Line(species, 1, ("2s2 2p1 3p1 1D2.0", "2s2 2p1 3s1 1P1.0"))
    NII_line_exc = ExcitationLine(NII)
    NII_line_rec = RecombinationLine(NII)
    models.append(NII_line_exc)
    models.append(NII_line_rec)
    plasma.models = models
    n1 = plasma.composition.get(nitrogen, 1)

    for i, x in enumerate(xr):
        for j, y in enumerate(yr):
            n1_emiss[j, i]   = NII_line_exc.emission(Point3D(x, 0.0, y), Vector3D(1.0, 0, 0), Spectrum(395, 400, 8000)).total() +NII_line_rec.emission(Point3D(x, 0.0, y), Vector3D(1.0, 0, 0), Spectrum(395, 400, 8000)).total() 
            n1_density[j, i] = n1.distribution.density(x, 0.0, y)

    
