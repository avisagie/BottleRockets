import numpy as np
from numpy import sin, cos, sqrt, pi
from rocket import Ballistic
from util import *

from matplotlib import pylab

timestep = 0.001
radius = 0.045
speed = 150
C_drag = 0.3
dry_mass = 0.4

def trajectory(speed, theta, timestep, radius, C_drag, dry_mass):
    phase = Ballistic([0, 0], [speed*cos(deg2rad(theta)), speed*sin(deg2rad(theta))], 0.0,
                    dry_mass=dry_mass, 
                    C_drag=C_drag,
                    A_cross_sectional_area=pi * (radius**2),
                    timestep=timestep)


    trace = []
    position = phase.position()
    altitude = position[1]
    count = 0

    max_altitude = 0.0
    # while the y coordinate is above ground
    while altitude >= 0:
        time, position, velocity = phase.step()
        altitude = position[1]

        t = np.append(position, sqrt(sum(velocity * velocity)))
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
for theta in angles:
    trace = trajectory(kmh2ms(speed), theta, timestep=timestep, radius=radius, C_drag=C_drag, dry_mass=dry_mass)
    pylab.plot(trace[:,0], trace[:,1])

pylab.legend(angles)
pylab.grid()
pylab.xlabel("Distance from origin (m)")
pylab.ylabel("Height (m)")
pylab.show()


