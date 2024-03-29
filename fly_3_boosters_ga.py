import numpy as np
from numpy import sin, cos, sqrt, pi
from util import *

from matplotlib import pylab

import rocket_architectures
import rocket
rocket.verbose = False

from genetic import Param, GeneticAlgorithm

radius = 0.068 / 2
center_length = 1.5 # m
center_volume = 1000 * pi * radius**2 * center_length # liters
booster_radius = 0.086 / 2 # m
booster_length = 0.52 # m
booster_volume = 1000 * pi * booster_radius**2 * booster_length # liters

print(f"Rocket length:{center_length:0.01f}m, volume:{center_volume:0.01f}l")
print(f"Booster length:{booster_length:0.01f}m, volume:{booster_volume:0.01f}l")

number_of_boosters = 3

params = {
        "radius": radius,
        "C_drag": 0.32, 
        "dry_mass": Param(0.2, 1.0),
        "volume": center_volume,
        "water_l": Param(center_volume / 8, 0.75 * center_volume),
        "pressure": 9, # relative pressure
        "nozzle_radius": Param(0.005, 0.011),
        "launch_tube_length": 1.0, # Param(0.0, center_length*0.9), # m

        "booster_radius": booster_radius,
        "booster_C_drag": 0.6,
        "booster_dry_mass": Param(0.440/3, 0.8),
        "booster_volume": booster_volume,
        "booster_water_l": Param(booster_volume / 8, 0.75 * booster_volume),
        "booster_nozzle_radius": 0.011, # Param(0.005, 0.011), # Param(0.005, 0.011),
        "booster_launch_tube_length": 0.3, # Param(0.0, booster_length*0.9), # m

        "theta": Param(30, 75), # degrees
        "rail_length": 0.3, # m

        "timestep": 0.001,

        "number_of_boosters": 3,
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
    max_distance = np.max(position[:,0])

    # print(f"Distance and speed: max height:{np.max(position[:,1]):0.1f}m, distance:{max_distance:0.0f}m, "
    #     + f"max speed:{max_speed:0.0f}m/s, {ms2kmh(max_speed):0.0f}km/h), acceleration:{max_acceleration/9.81:0.0f}g")

    return max_distance

ga = GeneticAlgorithm(params, fitness, population_size=40, generations=30, temperature_factor=0.91)
ga.run()

best_params = ga.get_best_params()
traces = rocket_architectures.sim_3_boosters(**best_params)
rocket_architectures.plot_basic(traces)
