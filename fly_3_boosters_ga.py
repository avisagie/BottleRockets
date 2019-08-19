import numpy as np
from numpy import sin, cos, sqrt, pi
from util import *

from matplotlib import pylab

import rocket_architectures
import rocket
rocket.verbose = False

from genetic import Param, GeneticAlgorithm

radius = 0.084 / 2

center_length = 1.8 # m
center_volume = 1000 * pi * radius**2 * center_length # liters
booster_length = 1 # m
booster_volume = 1000 * pi * radius**2 * booster_length # liters

print(f"Rocket length:{center_length:0.01f}m, volume:{center_volume:0.01f}l")
print(f"Booster length:{booster_length:0.01f}m, volume:{booster_volume:0.01f}l")

params = {
        "radius": radius,
        "C_drag": 0.4, 
        "dry_mass": Param(0.4, 0.8),
        "volume": center_volume,
        "water_l": Param(center_volume / 8, center_volume / 2),
        "pressure": 10, # relative pressure
        "nozzle_radius": Param(0.0045, 0.0105),
        "launch_tube_length": 1.3, # m

        "booster_radius": radius,
        "booster_C_drag": 0.3,
        "booster_dry_mass": Param(0.3, 0.8),
        "booster_volume": booster_volume,
        "booster_water_l": Param(booster_volume / 8, booster_volume / 2),
        "booster_nozzle_radius": Param(0.0045, 0.0105),
        "booster_launch_tube_length": 0.3, # m

        "theta": Param(30, 65), # degrees
        "rail_length": 1.5, # m

        "timestep": 0.001,
}

def fitness(params):    
    try:
        traces = rocket_architectures.sim_3_boosters(**params)
    except Exception as ex:
        print(ex)
        return 0.0

    time, position, velocity, acceleration = traces
    speed = sqrt(np.sum(velocity * velocity, axis=1))
    accel = sqrt(np.sum(acceleration * acceleration, axis=1)) / 9.81 # in Gs
    max_speed = max(speed)
    max_acceleration = max(sqrt(np.sum(acceleration * acceleration, axis=1)))
    max_position = np.max(position[:,0])

    # print(f"Distance and speed: max height:{np.max(position[:,1]):0.1f}m, distance:{max_position:0.0f}m, "
    #     + f"max speed:{max_speed:0.0f}m/s, {ms2kmh(max_speed):0.0f}km/h), acceleration:{max_acceleration/9.81:0.0f}g")

    return max_position

ga = GeneticAlgorithm(params, fitness)
ga.run()
