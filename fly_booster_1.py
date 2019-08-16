import numpy as np
from numpy import sin, cos, sqrt, pi
from rocket import Stepper, Ballistic, BoosterScienceBits, RocketWithComponents
from util import *
from itertools import chain

from matplotlib import pylab

radius = 0.045
C_drag = 0.3
dry_mass = 0.5
volume = 8
water_l = volume / 3
pressure = 10 # relative pressure
nozzle_radius = 0.0075

booster_radius = 0.045
booster_C_drag = 0.2
booster_dry_mass = 0.3
booster_volume = 3
booster_water_l = booster_volume / 3
booster_nozzle_radius = 0.0105

theta = 45 # degrees

timestep = 0.001

def trajectory(theta, timestep, radius, C_drag, dry_mass):

    stepper = Stepper()

    position = np.array([0, 0.1])
    altitude = position[1]
    count = 0

    max_altitude = 0.0

    boosters = [BoosterScienceBits( t0=0, 
                                    water=booster_water_l,
                                    pressure=bar2pa(pressure), 
                                    dry_mass=booster_dry_mass, 
                                    volume=booster_volume, 
                                    C_drag=booster_C_drag, 
                                    A_cross_sectional_area=pi*booster_radius**2, 
                                    nozzle_radius=booster_nozzle_radius, 
                                    timestep=timestep) for x in range(3)]

    center = BoosterScienceBits( t0=0, 
                                 water=water_l,
                                 pressure=bar2pa(pressure), 
                                 dry_mass=dry_mass, 
                                 volume=volume, 
                                 C_drag=C_drag, 
                                 A_cross_sectional_area=pi*radius**2, 
                                 nozzle_radius=nozzle_radius, 
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

    phase = RocketWithComponents(np.array([0,0.3]), 0.001*np.array([cos(theta), sin(theta)]), 0.0, 
        components=list(chain(boosters, [center])), 
        rail_length=1.5, 
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


traces = trajectory(theta, timestep=timestep, radius=radius, C_drag=C_drag, dry_mass=dry_mass)
time, position, velocity, acceleration = traces
speed = sqrt(np.sum(velocity * velocity, axis=1))
max_speed = max(speed)
max_acceleration = max(sqrt(np.sum(acceleration * acceleration, axis=1)))

print(position)

ax1 = pylab.subplot(211)
ax1.plot(time, position[:, 1], 'b')
ax1.set_ylabel("Height (m)", color='b')
ax2 = ax1.twinx()
ax2.plot(time, speed, 'r')
ax2.set_ylabel("Speed (m/s)", color='r')
ax1.grid()
ax1.set_title(f"Distance:{np.max(position[:,0]):0.0f}m, flight time:{max(time):0.01f}s")
print(f"Distance and speed (max height:{np.max(position[:,1]):0.1f}m, distance:{np.max(position[:,0]):0.0f}m, "
+ f"max speed:{max_speed:0.0f}m/s, {ms2kmh(max_speed):0.0f}km/h), acceleration:{max_acceleration/9.81:0.0f}g")
ax1.set_xlabel("Time (s)")

ax1 = pylab.subplot(212)
ax1.plot(position[:, 0], position[:, 1], 'b')
ax1.set_ylabel("Height (m)", color='b')
ax1.grid()
ax1.set_xlabel("Distance (m)")

pylab.show()


