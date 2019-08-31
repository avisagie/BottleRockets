import numpy as np
from numpy import sin, cos, sqrt, pi
from util import *

from matplotlib import pylab

from rocket_architectures import sim_3_boosters, sim1, plot_basic

if __name__ == "__main__":

    radius = 0.068 / 2

    center_length = 1.0 # m
    center_volume = 1500 * pi * radius**2 * center_length # liters
    booster_length = 1 # m
    booster_volume = 1000 * pi * radius**2 * booster_length # liters

    num_fins = 3
    fin_thickness = 0.003
    fin_length = 0.110

    traces = sim_3_boosters(
        C_drag = 0.320,
        booster_C_drag = 0.600,
        booster_dry_mass = 0.080, # optimised
        booster_launch_tube_length = 0.300,
        booster_nozzle_radius = 0.011, # optimised
        booster_radius = 0.034,
        booster_volume = 2.000,
        booster_water_l = 0.838, # optimised
        dry_mass = 0.400, # optimised
        launch_tube_length = 1.000,
        nozzle_radius = 0.005,
        pressure = 10.000,
        radius = 0.034,
        rail_length = 0.300,
        theta = 45.000,
        timestep = 0.001,
        volume = 3.632,
        water_l = 1.284, # optimised
    )

    time, position, velocity, acceleration = traces
    speed = sqrt(np.sum(velocity * velocity, axis=1))
    accel = sqrt(np.sum(acceleration * acceleration, axis=1)) / 9.81 # in Gs
    max_speed = max(speed)
    max_acceleration = max(sqrt(np.sum(acceleration * acceleration, axis=1)))
    # print(f"Center stage length:{center_length:0.01f}m, volume:{center_volume:0.01f}l")
    # print(f"Booster length:{booster_length:0.01f}m, volume:{booster_volume:0.01f}l")
    print(f"Distance and speed: max height:{np.max(position[:,1]):0.1f}m, distance:{np.max(position[:,0]):0.0f}m, "
    + f"max speed:{max_speed:0.0f}m/s, {ms2kmh(max_speed):0.0f}km/h), acceleration:{max_acceleration/9.81:0.0f}g")

    plot_basic(traces)
