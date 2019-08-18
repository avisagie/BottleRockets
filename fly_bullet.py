import numpy as np
from numpy import sin, cos, sqrt, pi
from util import *

from matplotlib import pylab

from rocket_architectures import sim_3_boosters, sim_3_boosters_bullet, sim1, plot_basic

if __name__ == "__main__":

    radius = 0.084 / 2

    bullet_length = 1.5 # m
    bullet_volume = 1000 * pi * radius**2 * bullet_length # liters
    booster_length = 1.5 # m
    booster_volume = 1000 * pi * radius**2 * booster_length # liters

    traces = sim_3_boosters_bullet(
        radius = radius,
        C_drag = 0.3, 
        dry_mass = 0.5,
        volume = bullet_volume,
        water_l = bullet_volume / 3,
        pressure = 10, # relative pressure
        nozzle_radius = 0.0105,
        launch_tube_length = 0.0, # m

        booster_radius = radius,
        booster_C_drag = 0.4,
        booster_dry_mass = 0.5,
        booster_volume = booster_volume,
        booster_water_l = 3.0 / 3,
        booster_nozzle_radius = 0.0105,
        booster_launch_tube_length = 1.3, # m

        theta = 45, # degrees
        rail_length = 1.5, # m

        timestep = 0.001
    )

    time, position, velocity, acceleration = traces
    speed = sqrt(np.sum(velocity * velocity, axis=1))
    accel = sqrt(np.sum(acceleration * acceleration, axis=1)) / 9.81 # in Gs
    max_speed = max(speed)
    max_acceleration = max(sqrt(np.sum(acceleration * acceleration, axis=1)))
    print(f"Center stage length:{bullet_length:0.01f}m, volume:{bullet_volume:0.01f}l")
    print(f"Booster length:{booster_length:0.01f}m, volume:{booster_volume:0.01f}l")
    print(f"Distance and speed: max height:{np.max(position[:,1]):0.1f}m, distance:{np.max(position[:,0]):0.0f}m, "
    + f"max speed:{max_speed:0.0f}m/s, {ms2kmh(max_speed):0.0f}km/h), acceleration:{max_acceleration/9.81:0.0f}g")

    plot_basic(traces)