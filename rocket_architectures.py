import numpy as np
from numpy import sin, cos, sqrt, pi
from rocket import Stepper, Ballistic, BoosterScienceBits, RocketWithComponents
from util import *
from itertools import chain

from matplotlib import pylab

def sim_3_boosters(
    radius = 0.045,
    C_drag = 0.3,
    dry_mass = 0.5,
    volume = 8,
    water_l = 8.0 / 3,
    pressure = 10, # relative pressure
    nozzle_radius = 0.0105,
    launch_tube_length = 1.3, # m

    booster_radius = 0.045,
    booster_C_drag = 0.3,
    booster_dry_mass = 0.3,
    booster_volume = 3,
    booster_water_l = 3.0 / 3,
    booster_nozzle_radius = 0.0105,
    booster_launch_tube_length = 0.3, # m

    theta = 40, # degrees
    rail_length = 1.5, # m

    timestep = 0.001
    ):

    stepper = Stepper()

    position = np.array([0, 0.1])
    altitude = position[1]
    
    boosters = [BoosterScienceBits( t0=0, 
                                    water=booster_water_l,
                                    pressure=bar2pa(pressure), 
                                    dry_mass=booster_dry_mass, 
                                    volume=booster_volume, 
                                    C_drag=booster_C_drag, 
                                    A_cross_sectional_area=pi*booster_radius**2, 
                                    nozzle_radius=booster_nozzle_radius, 
                                    launch_tube_length=booster_launch_tube_length,
                                    timestep=timestep) for x in range(3)]

    center = BoosterScienceBits( t0=0, 
                                 water=water_l,
                                 pressure=bar2pa(pressure), 
                                 dry_mass=dry_mass, 
                                 volume=volume, 
                                 C_drag=C_drag, 
                                 A_cross_sectional_area=pi*radius**2, 
                                 nozzle_radius=nozzle_radius, 
                                 launch_tube_length=launch_tube_length,
                                 timestep=timestep )

    def validate():
        # TODO revisit, and make sure it gets called
        # Check that the center stage does nto pull away from the boosters
        a_center = center.F_thrust() / center.mass()
        for booster in boosters: # redundant because they're identical
            a_booster = booster.F_thrust() / booster.mass()
            if a_booster < a_center:
                return False
        return True

    phase = RocketWithComponents(position, position, 0.001*np.array([cos(deg2rad(theta)), sin(deg2rad(theta))]), 0.0, 
        components=list(chain(boosters, [center])), 
        rail_length=rail_length, 
        validate=validate, 
        timestep=timestep  
    )

    stepper.step(phase)

    phase = Ballistic(phase.position(), phase.velocity(), phase.t,
                    dry_mass=dry_mass, 
                    C_drag=C_drag,
                    A_cross_sectional_area=pi * (radius**2),
                    timestep=timestep)

    stepper.step(phase)

    return stepper.get_traces()


def sim1(
    radius = 0.045,
    C_drag = 0.3,
    dry_mass = 0.5,
    volume = 8,
    water_l = 8.0 / 3,
    pressure = 10, # relative pressure
    nozzle_radius = 0.0105,
    launch_tube_length = 0.0, # m

    theta = 40, # degrees
    rail_length = 1.5, # m

    timestep = 0.001
    ):

    stepper = Stepper()

    position = np.array([0, 0.1])
    
    center = BoosterScienceBits( t0=0, 
                                 water=water_l,
                                 pressure=bar2pa(pressure), 
                                 dry_mass=dry_mass, 
                                 volume=volume, 
                                 C_drag=C_drag, 
                                 A_cross_sectional_area=pi*radius**2, 
                                 nozzle_radius=nozzle_radius, 
                                 launch_tube_length=launch_tube_length,
                                 timestep=timestep )

    launcher_origin = np.array([0,0.3])
    phase = RocketWithComponents(position, position, 0.001*np.array([cos(deg2rad(theta)), sin(deg2rad(theta))]), 0.0, 
        components=[center], 
        rail_length=rail_length, 
        timestep=timestep  
    )

    stepper.step(phase)

    phase = Ballistic(phase.position(), phase.velocity(), phase.t,
                    dry_mass=dry_mass, 
                    C_drag=C_drag,
                    A_cross_sectional_area=pi * (radius**2),
                    timestep=timestep)

    stepper.step(phase)

    return stepper.get_traces()


if __name__ == "__main__":
    traces = sim1() # sim_3_boosters()
    time, position, velocity, acceleration = traces
    speed = sqrt(np.sum(velocity * velocity, axis=1))
    accel = sqrt(np.sum(acceleration * acceleration, axis=1)) / 9.81 # in Gs
    max_speed = max(speed)
    max_acceleration = max(sqrt(np.sum(acceleration * acceleration, axis=1)))

    print(position)

    ax1 = pylab.subplot(211)
    ax1.plot(position[:, 0], accel, 'b')
    ax1.set_ylabel("Acceleration (g)", color='b')
    ax2 = ax1.twinx()
    ax2.plot(position[:, 0], speed, 'r')
    ax2.set_ylabel("Speed (m/s)", color='r')
    ax1.grid()
    ax1.set_title(f"Distance:{np.max(position[:,0]):0.0f}m, flight time:{max(time):0.01f}s")
    print(f"Distance and speed (max height:{np.max(position[:,1]):0.1f}m, distance:{np.max(position[:,0]):0.0f}m, "
    + f"max speed:{max_speed:0.0f}m/s, {ms2kmh(max_speed):0.0f}km/h), acceleration:{max_acceleration/9.81:0.0f}g")
    ax1.set_xlabel("Distance (s)")

    ax1 = pylab.subplot(212)
    ax1.plot(position[:, 0], position[:, 1], 'b')
    ax1.set_ylabel("Height (m)", color='b')
    ax1.grid()
    ax1.set_xlabel("Distance (m)")

    pylab.show()
