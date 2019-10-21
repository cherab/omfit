
from raysect.optical.material import AbsorbingSurface


def load_machine(config, world):

    if config['machine']['reflecting']:
        override_material = AbsorbingSurface()
    else:
        override_material = None

    if config['machine']['name'] == "JET":
        from cherab.jet.machine import import_jet_mesh
        import_jet_mesh(world, override_material=override_material)

    elif config['machine']['name'] == "MAST-U":
        from cherab.mastu.machine import import_mastu_mesh
        import_mastu_mesh(world, override_material=override_material)

    elif config['machine']['name'] == "AUG":
        from cherab.aug.machine import import_aug_mesh
        import_aug_mesh(world, override_material=override_material)

    else:

        raise ValueError("The machine specified '{}' is not currently available in this package."
                         "".format(config['machine']['name']))
