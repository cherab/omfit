

def load_edge_simulation(config, world):

    # Only try to do something if edge plasma has been selected as True
    try:
        if not config['plasma']['edge']['present']:
            raise ValueError("This config file does not specify an edge plasma.")
    except KeyError:
        raise ValueError("This config file does not specify an edge plasma.")

    if config['plasma']['edge']['type'] == "SOLPS":
        return _load_solps_simulation(config, world)

    else:
        raise ValueError("Unrecognised simulation type.")


def _load_solps_simulation(config, world):

    if config['plasma']['edge']['SOLPS_format'] == 'MDSplus':

        print()
        print("loading MDPsplus SOLPS plasma")

        from cherab.solps import load_solps_from_mdsplus

        mds_server = config['plasma']['edge']['mds_server']
        mds_solps_reference = config['plasma']['edge']['mds_solps_reference']

        sim = load_solps_from_mdsplus(mds_server, mds_solps_reference)
        plasma = sim.create_plasma(parent=world)

    elif config['plasma']['edge']['SOLPS_format'] == 'Files':

        from cherab.solps import load_solps_from_raw_output

        solps_directory = config['plasma']['edge']['solps_directory']
        sim = load_solps_from_raw_output(solps_directory, debug=True)
        plasma = sim.create_plasma(parent=world)

    else:
        raise ValueError("Unrecognised SOLPS format '{}'.".format(config['plasma']['edge']['SOLPS_format']))

    return plasma

