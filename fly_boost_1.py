import numpy as np
from numpy import sin, cos, sqrt, pi
from rocket import Stepper, Ballistic, BoostScienceBits
from util import *

from matplotlib import pylab

timestep = 0.001
radius = 0.05
C_drag = 0.3
dry_mass = 0.2
volume = 2
water_l = volume / 5
pressure = 9 - 1.0 # relative pressure
nozzle_radius = 0.01
theta = 45 # degrees


def trajectory(theta, timestep, radius, C_drag, dry_mass):

    stepper = Stepper()

    position = np.array([0, 0.1])
    altitude = position[1]
    count = 0

    max_altitude = 0.0

    phase = BoostScienceBits(   position=position, 
                                velocity=0.001*np.array([cos(deg2rad(theta)), sin(deg2rad(theta))]), # almost still, just has a direction for the start
                                t0=0, 
                                water=water_l,
                                pressure=bar2pa(pressure), 
                                dry_mass=dry_mass, 
                                volume=volume, 
                                C_drag=C_drag, 
                                A_cross_sectional_area=pi*radius**2, 
                                nozzle_radius=nozzle_radius, 
                                rail_length=1.5,
                                timestep=timestep)

    stepper.step(phase)

    phase = Ballistic(phase.position(), phase.velocity(), phase.t,
                    dry_mass=dry_mass, 
                    C_drag=C_drag,
                    A_cross_sectional_area=pi * (radius**2),
                    timestep=timestep)

    stepper.step(phase)

    return stepper.get_traces()


traces = trajectory(theta, timestep=timestep, radius=radius, C_drag=C_drag, dry_mass=dry_mass)
time, position, velocity = traces
speed = sqrt(np.sum(velocity * velocity, axis=1))
max_speed = max(speed)

print(position)

ax1 = pylab.subplot(211)
ax1.plot(time, position[:, 1], 'b')
ax1.set_ylabel("Height (m)", color='b')
ax2 = ax1.twinx()
ax2.plot(time, speed, 'r')
ax2.set_ylabel("Speed (m/s)", color='r')
ax1.grid()
ax1.set_title(f"Distance and speed (max height:{np.max(position[:,1]):0.1f}m, distance:{np.max(position[:,0]):0.0f}m, "
                + f"max speed:{max_speed:0.0f}m/s, {ms2kmh(max_speed):0.0f}km/h)")
ax1.set_xlabel("Time (s)")

ax1 = pylab.subplot(212)
ax1.plot(position[:, 0], position[:, 1], 'b')
ax1.set_ylabel("Height (m)", color='b')
ax1.grid()
ax1.set_xlabel("Distance (m)")

pylab.show()


