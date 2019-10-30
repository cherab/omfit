from cherab.core.atomic.elements import hydrogen, deuterium, carbon, helium, nitrogen, neon, argon, krypton, xenon
import numpy as np

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
def load_profiles(config, plasma, xr, yr, num):

    te_plasma = np.zeros((num,num))
    ne_plasma = np.zeros((num,num))
    ni_plasma = np.zeros((num,num,2))

    if config['plasma']['edge']['Te2D']:
        for i, x in enumerate(xr):
            for j, y in enumerate(yr):
                te_plasma[j, i]  = plasma.electron_distribution.effective_temperature(x, 0.0, y)
    
    if config['plasma']['edge']['ne2D']:
        for i, x in enumerate(xr):
            for j, y in enumerate(yr):
                ne_plasma[j, i]  = plasma.electron_distribution.density(x, 0.0, y)

    if config['plasma']['edge']['ni2D']:

        ni_plasma = np.zeros((num,num,2))

        for i, x in enumerate(xr):
            for j, y in enumerate(yr):
                for k in range(2):
                    composition = plasma.composition.get(deuterium, k)
                    ni_plasma[j, i, k] = composition.distribution.density(x, 0.0, y)

    if config['plasma']['edge']['nz2D']:
            try:
                species = _SPECIES_LOOKUP[config["plasma"]['edge']["nz2D_species"]]
            except KeyError:
                raise ValueError("The species emission specification for this config file is invalid.")

            nz_plasma = np.zeros((num,num,species.atomic_number+1))
            ions = np.linspace(0,species.atomic_number,num=species.atomic_number+1)  

            for i, x in enumerate(xr):
                for j, y in enumerate(yr):
                    for k,z in enumerate(ions):
                        composition = plasma.composition.get(species, k)
                        nz_plasma[j, i, k]  = composition.distribution.density(x, 0.0, y)
    else:
        nz_plasma = np.zeros((num,num,1))

    return te_plasma, ne_plasma, ni_plasma, nz_plasma
