import numpy as np
from numpy import sin, cos, sqrt, pi
from util import *

from matplotlib import pylab

from rocket_architectures import sim_3_boosters, sim1, plot_basic

if __name__ == "__main__":

    radius = 0.084 / 2

    center_length = 1.8 # m
    center_volume = 1000 * pi * radius**2 * center_length # liters
    booster_length = 1 # m
    booster_volume = 1000 * pi * radius**2 * booster_length # liters

    traces = sim_3_boosters(
     C_drag = 0.4,
     booster_C_drag = 0.3,
    booster_dry_mass = 0.25,
     booster_launch_tube_length = 1.0,
    booster_nozzle_radius = 0.011,
     booster_radius = 0.042,
     booster_volume = 5.541769440932395,
    booster_water_l = 2.815664357194821,
    dry_mass = 0.7379773338578499,
     launch_tube_length = 1.0,
    nozzle_radius = 0.00435,
     pressure = 10,
     radius = 0.042,
     rail_length = 1.5,
    theta = 46.36640198337145,
     timestep = 0.001,
     volume = 9.975184993678312,
    water_l = 4.609209988958819,    )

    time, position, velocity, acceleration = traces
    speed = sqrt(np.sum(velocity * velocity, axis=1))
    accel = sqrt(np.sum(acceleration * acceleration, axis=1)) / 9.81 # in Gs
    max_speed = max(speed)
    max_acceleration = max(sqrt(np.sum(acceleration * acceleration, axis=1)))
    print(f"Center stage length:{center_length:0.01f}m, volume:{center_volume:0.01f}l")
    print(f"Booster length:{booster_length:0.01f}m, volume:{booster_volume:0.01f}l")
    print(f"Distance and speed: max height:{np.max(position[:,1]):0.1f}m, distance:{np.max(position[:,0]):0.0f}m, "
    + f"max speed:{max_speed:0.0f}m/s, {ms2kmh(max_speed):0.0f}km/h), acceleration:{max_acceleration/9.81:0.0f}g")

    plot_basic(traces)