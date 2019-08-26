import math
import numpy as np
from numpy import sin, cos, sqrt, pi
from util import *

from matplotlib import pylab

import rocket_architectures
import rocket
rocket.verbose = False

from genetic import Param, GeneticAlgorithm

radius = 0.084 / 2

stage_length = 1.5 # m
stage_volume = 1000 * pi * radius**2 * stage_length # liters

print(f"Rocket length:{stage_length:0.01f}m, volume:{stage_length:0.01f}l")

params = {
    'pressure': 6, # relative pressure

    's1_radius': radius,
    's1_C_drag': 0.6,
    's1_dry_mass_base': 0.100,
    's1_volume': Param(1.0, 6.0),
    's1_water_l': Param(stage_volume/8, 0.75*stage_volume),
    's1_nozzle_radius': Param(0.0087/2, 0.0105),

    's2_radius': radius,
    's2_C_drag': 0.5,
    's2_dry_mass_base': 0.100,
    's2_volume': Param(1.0, 4.0),
    's2_water_l': Param(stage_volume/8, 0.75*stage_volume),
    's2_nozzle_radius': Param(0.0087/2, 0.0105),

    's3_radius': radius,
    's3_C_drag': 0.4,
    's3_dry_mass_base': Param(0.1, 1.0),
    's3_volume': Param(1.0, 2.0),
    's3_water_l': Param(stage_volume/8, 0.75*2.0),
    's3_nozzle_radius': Param(0.0087/2, 0.0105),

    'theta': Param(30.0, 90.0), # degrees
    'rail_length': 3, # m

    'timestep': 0.001,
}

def fitness(params):    
    try:
        traces = rocket_architectures.sim_3_stage(**params)
    except Exception as ex:
        print(ex)
        return 0.0

    time, position, velocity, acceleration = traces
    speed = sqrt(np.sum(velocity * velocity, axis=1))
    accel = sqrt(np.sum(acceleration * acceleration, axis=1)) / 9.81 # in Gs
    max_speed = max(speed)
    max_acceleration = max(sqrt(np.sum(acceleration * acceleration, axis=1)))
    max_distance = np.max(position[:,0])
    max_height = np.max(position[:,1])

    # print(f"Distance and speed: max height:{np.max(position[:,1]):0.1f}m, distance:{max_distance:0.0f}m, "
    #     + f"max speed:{max_speed:0.0f}m/s, {ms2kmh(max_speed):0.0f}km/h), acceleration:{max_acceleration/9.81:0.0f}g")

    if math.isnan(max_distance):
        return 0.0

    return max_distance

ga = GeneticAlgorithm(params, fitness, population_size=40, generations=100, temperature_factor=0.94)
ga.run()

best_params = ga.get_best_params()
traces = rocket_architectures.sim_3_stage(**best_params)
rocket_architectures.plot_basic(traces)
