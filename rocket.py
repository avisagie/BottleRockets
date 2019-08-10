import numpy as np
from numpy import sin, cos, sqrt, pi
import scipy.integrate

# Assume it's constant for the flight and use dry air, which is more dense than humid air.
# https://en.wikipedia.org/wiki/Density_of_air
# Can be refined later
air_density = 1.2 # kg/m^3 

class Ballistic:
    def __init__(self, position, velocity, t0, dry_mass, C_drag, A_cross_sectional_area, timestep = 0.001):
        """
        Starting conditions: 
            position: array (2,), displacement from the origin, m
            velocity: array (2,), m/s            
            time: t0 for this part of the flight
        """
        self.state = np.array([position, velocity])
        self.t = t0
        self.dry_mass = dry_mass
        self.C_drag = C_drag
        self.A_cross_sectional_area = A_cross_sectional_area
        self.timestep = timestep


    def position(self):
        return self.state[0]


    def velocity(self):
        return self.state[1]


    def fun(self):
        """
        The right hand side of the system of ODEs. 
        (https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.RK45.html#scipy.integrate.RK45)
        """
        velocity = self.velocity()

        speed2 = np.sum(velocity * velocity)
        direction = 1/sqrt(speed2) * velocity
        dvdt = -np.array([0, 9.81]) - 1/(2*self.dry_mass) * self.C_drag * self.A_cross_sectional_area * speed2 * direction

        position = self.position()
        dsdt = velocity

        return np.array([dsdt, dvdt])


    def step(self):
        """
        Step the system one step forward.
        Return: position, velocity, t
        """

        self.state = self.state + self.timestep * self.fun() 
        self.t += self.timestep

        return self.t, self.position(), self.velocity()
