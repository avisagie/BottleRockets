import numpy as np
from numpy import sin, cos, sqrt, pi
from util import *

from matplotlib import pylab

from rocket_architectures import sim_3_boosters, sim_3_boosters_bullet, sim1, plot_basic

if __name__ == "__main__":

    radius = 0.084 / 2

    length = 1.1 # m
    volume = 5.38 # 1000 * pi * radius**2 * length # liters

    num_fins = 3
    fin_thickness = 0.002
    fin_length = 0.08

    traces = sim1(
        radius = radius,
        C_drag = 0.37, 
        dry_mass = 0.360,
        volume = volume,
        water_l = volume/3,
        pressure = 10, # relative pressure
        nozzle_radius = 0.0105,
        launch_tube_length = 1.0, # m

        theta = 43, # degrees
        rail_length = 0.5, # m

        extra_frontal_surface = num_fins * fin_length * fin_thickness,

        timestep = 0.001
    )

    time, position, velocity, acceleration = traces
    speed = sqrt(np.sum(velocity * velocity, axis=1))
    accel = sqrt(np.sum(acceleration * acceleration, axis=1)) / 9.81 # in Gs
    max_speed = max(speed)
    max_acceleration = max(sqrt(np.sum(acceleration * acceleration, axis=1)))
    print(f"Rocket length:{length:0.01f}m, volume:{volume:0.01f}l")
    print(f"Distance and speed: max height:{np.max(position[:,1]):0.1f}m, distance:{np.max(position[:,0]):0.0f}m, "
    + f"max speed:{max_speed:0.0f}m/s, {ms2kmh(max_speed):0.0f}km/h), acceleration:{max_acceleration/9.81:0.0f}g")

    plot_basic(traces)