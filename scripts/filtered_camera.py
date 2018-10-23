#!/bin/env python

import json
import argparse

from raysect.optical import World
from cherab.omfit import load_machine, load_edge_simulation, load_emission, load_camera


# handle arguments
parser = argparse.ArgumentParser('Simulate a filtered camera')
parser.add_argument('configfile', type=str, help="A JSON format CHERAB configuration file.")
args = parser.parse_args()

# read/parse config file
with open(args.configfile, 'r') as f:
    config = json.load(f)


world = World()

# Load all the CAD file geometry for the selected machine
load_machine(config, world)

# Load the edge plasma solution if present
plasma = load_edge_simulation(config, world)

# Load all the emission line models
load_emission(config, plasma)

# Load the specified filtered camera
camera = load_camera(config, world)


if config["observer"]["display_progress"]:
    import matplotlib.pyplot as plt
    plt.ion()

camera.observe()

if config["observer"]["display_progress"]:
    plt.ioff()
    plt.show()

