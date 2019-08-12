import numpy as np
from numpy import sin, cos, sqrt, pi
from rocket import Ballistic, BoostScienceBits
from util import *

from matplotlib import pylab

timestep = 0.001
radius = 0.05
C_drag = 1.0
dry_mass = 0.1
water_l = 0.66
pressure = 2.5331 - 1.0 # relative pressure
volume = 2.0
nozzle_radius = 0.01

def trajectory(theta, timestep, radius, C_drag, dry_mass):

    trace = []
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
                                timestep=timestep)

    while altitude > 0.0:
        # print(phase.state)

        ret = phase.step()
        if ret is None:
            # water's run out
            break

        time, position, velocity = ret
        altitude = position[1]

        t = np.append(position, [sqrt(sum(velocity * velocity)), time])
        trace.append(t)
        max_altitude = max(max_altitude, altitude)

        count += 1
        if count*timestep >= 0.001:
            print(f'{time:0.04f}: {position}, {velocity}')    
            count = 0


    phase = Ballistic(position, velocity, phase.t,
                    dry_mass=dry_mass, 
                    C_drag=C_drag,
                    A_cross_sectional_area=pi * (radius**2),
                    timestep=timestep)


    # while the y coordinate is above ground
    while altitude >= 0:
        time, position, velocity = phase.step()
        altitude = position[1]

        t = np.append(position, [sqrt(sum(velocity * velocity)), time])
        trace.append(t)
        max_altitude = max(max_altitude, altitude)

        count += 1
        if count*timestep > 0.1:
            print(f'{time:0.04f}: {position}, {velocity}')    
            count = 0

    flight_time = time
    print(f'{time:0.04f}: {position}, {velocity}')    
    print(f'Flight time: {flight_time}, Max altitude: {max_altitude}')

    return np.array(trace)

pylab.title(f"Trajectories at various launch angles")

angles = [35, 40, 45, 55, 65, 75, 90]
for theta in [90]:
    trace = trajectory(theta, timestep=timestep, radius=radius, C_drag=C_drag, dry_mass=dry_mass)
    pylab.plot(trace[:,3], trace[:,1])
    pylab.plot(trace[:,3], trace[:,2])

pylab.grid()
pylab.xlabel("Distance from origin (m)")
pylab.ylabel("Height (m)")
pylab.show()


