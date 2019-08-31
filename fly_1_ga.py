import numpy as np
from numpy import sin, cos, sqrt, pi
from util import *

from matplotlib import pylab

import rocket_architectures
import rocket
rocket.verbose = False

from genetic import Param, GeneticAlgorithm

radius = 0.084 / 2

length = 1.0 # m
volume = 1000 * pi * radius**2 * length # liters

params = {
    "radius": radius,
    "C_drag": 0.32, 
    "dry_mass": Param(0.1, 1.0),
    "volume": volume,
    "water_l": Param(0.125*volume, 0.75*volume),
    "pressure": 10.0, # relative pressure
    "nozzle_radius": Param(0.005, 0.0105),
    "launch_tube_length": 1.0, # Param(0.0, 0.9*length), # m

    "theta": Param(30.0, 75), # degrees
    "rail_length": 1.0, # m

    "timestep": 0.001,
}

def fitness(params):    
    traces = rocket_architectures.sim1(**params)

    time, position, velocity, acceleration = traces
    speed = sqrt(np.sum(velocity * velocity, axis=1))
    accel = sqrt(np.sum(acceleration * acceleration, axis=1)) / 9.81 # in Gs
    max_speed = max(speed)
    max_acceleration = max(sqrt(np.sum(acceleration * acceleration, axis=1)))
    max_distance = np.max(position[:,0])
    max_height = np.max(position[:,1])

    # print(f"Distance and speed: max height:{np.max(position[:,1]):0.1f}m, distance:{max_distance:0.0f}m, "
    #     + f"max speed:{max_speed:0.0f}m/s, {ms2kmh(max_speed):0.0f}km/h), acceleration:{max_acceleration/9.81:0.0f}g")

    return max_distance

ga = GeneticAlgorithm(params, fitness)
ga.run()

print(f"Rocket length:{length:0.01f}m, volume:{volume:0.01f}l")

best_params = ga.get_best_params()
traces = rocket_architectures.sim1(**best_params)
rocket_architectures.plot_basic(traces)

