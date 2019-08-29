
from cherab.core.atomic import Line
from cherab.core.atomic.elements import hydrogen, deuterium, carbon, helium, nitrogen, neon, argon, krypton, xenon
from cherab.core.model import ExcitationLine, RecombinationLine, MultipletLineShape, StarkBroadenedLine, Bremsstrahlung

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


def load_emission(config, plasma):

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
    plasma.models = models
