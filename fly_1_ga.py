import numpy as np
from numpy import sin, cos, sqrt, pi
from util import *

from matplotlib import pylab

import rocket_architectures
import rocket
rocket.verbose = False

from genetic import Param, GeneticAlgorithm

radius = 0.09 / 2

length = 1.35 # m
volume = 1000 * pi * radius**2 * length # liters

print(f"Rocket length:{length:0.01f}m, volume:{volume:0.01f}l")

params = {
    "radius": radius,
    "C_drag": 0.3, 
    "dry_mass": Param(0.4, 0.8),
    "volume": volume,
    "water_l": Param(volume/8, volume/2),
    "pressure": 8, # relative pressure
    "nozzle_radius": Param(0.0045, 0.0105),
    "launch_tube_length": Param(0.0, 0.9*length), # m

    "theta": Param(30.0, 75), # degrees
    "rail_length": 2.0, # m

    "timestep": 0.001,
}

def fitness(params):    
    traces = rocket_architectures.sim1(**params)

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
