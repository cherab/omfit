
from cherab.core.atomic import Line
from cherab.core.atomic.elements import hydrogen, deuterium, carbon, helium, nitrogen, neon, argon, krypton, xenon
from cherab.core.model import ExcitationLine, RecombinationLine

from cherab.openadas import OpenADAS


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

    for emission_instruction in config["plasma"]["emission_instructions"]:

        try:
            species = _SPECIES_LOOKUP[emission_instruction["species"]]
        except KeyError:
            raise ValueError("The species emission specification for this config file is invalid.")

        ionisation = emission_instruction["ionisation"]
        upper = emission_instruction["upper"]
        lower = emission_instruction["lower"]

        # Do a test call to see if the species exists on this plasma.
        _ = plasma.composition.get(species, ionisation)

        line = Line(species, ionisation, (upper, lower))

        models.append(_EMISSION_TYPE_LOOKUP[emission_instruction["type"]](line))

    plasma.models = models
