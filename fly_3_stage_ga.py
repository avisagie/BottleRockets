import numpy as np
from numpy import sin, cos, sqrt, pi
from util import *

from matplotlib import pylab

import rocket_architectures
import rocket
rocket.verbose = False

from genetic import Param, GeneticAlgorithm

radius = 0.086 / 2

stage_length = 1.5 # m
stage_volume = 1000 * pi * radius**2 * stage_length # liters

print(f"Rocket length:{stage_length:0.01f}m, volume:{stage_length:0.01f}l")

params = {
    'pressure': 10, # relative pressure

    's1_radius': radius,
    's1_C_drag': 0.6,
    's1_dry_mass': 0.300,
    's1_volume': stage_volume,
    's1_water_l': Param(stage_volume/8, 0.75*stage_volume),
    's1_nozzle_radius': Param(0.0087/2, 0.0105),

    's2_radius': radius,
    's2_C_drag': 0.5,
    's2_dry_mass': 0.300,
    's2_volume': stage_volume,
    's2_water_l': Param(stage_volume/8, 0.75*stage_volume),
    's2_nozzle_radius': Param(0.0087/2, 0.0105),

    's3_radius': radius,
    's3_C_drag': 0.4,
    's3_dry_mass': Param(0.25, 1.5 - 2*0.38),
    's3_volume': stage_volume,
    's3_water_l': Param(stage_volume/8, 0.75*stage_volume),
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

    return max_height

ga = GeneticAlgorithm(params, fitness, population_size=30, generations=30, temperature_factor=0.90)
ga.run()

best_params = ga.get_best_params()
traces = rocket_architectures.sim_3_stage(**best_params)
rocket_architectures.plot_basic(traces)
