import numpy as np
from numpy import sin, cos, sqrt, pi
import scipy.integrate

# Assume it's constant for the flight and use dry air, which is more dense than humid air.
# https://en.wikipedia.org/wiki/Density_of_air
# Can be refined later
air_density = 1.2 # kg / m^3 

# pretty close for our temperatures and height above sea level
water_density = 998.0 # kg / m^3

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
        a_drag = - 1/(2*self.dry_mass) * self.C_drag * self.A_cross_sectional_area * speed2 * direction

        a_grav = -np.array([0, 9.81])

        dvdt = a_grav + a_drag

        dsdt = velocity

        return np.array([dsdt, dvdt])


    def step(self):
        """
        Step the system one step forward.
        Return: t, position, velocity
        """

        self.state = self.state + self.timestep * self.fun() 
        self.t += self.timestep

        return self.t, self.position(), self.velocity()


class BoostScienceBits:

    """
    Boost phase from science bits: http://www.sciencebits.com/RocketEqs 
    """

    def __init__(self, position, velocity, t0, 
                 water, # in liter
                 pressure, # in kPa relative to outside pressure
                 dry_mass, volume, # in liter, total volume
                 C_drag, A_cross_sectional_area, nozzle_radius, # in meters
                 timestep = 0.001):
        # (mass and pressure are vectors in case I want to use scipy integrators)
        mass = dry_mass + water/1000 * water_density # 1l = 0.001m^3
        self.state = np.array([position, velocity, [mass, 0], [pressure, 0]]) 
        self.t = t0
        self.dry_mass = dry_mass
        self.mass_0 = mass
        self.volume = volume / 1000 # total volume in m^3, 1l = 0.001m^3
        self.volume_0 = self.volume - water/1000 # air volume at t0
        self.pressure_0 = pressure
        self.C_drag = C_drag
        self.A_cross_sectional_area = A_cross_sectional_area
        self.timestep = timestep
        self.nozzle_radius = nozzle_radius
        self.nozzle_area = pi*self.nozzle_radius**2
        self.gamma = 1.4 # adiabatic constant for dry air


    def position(self):
        return self.state[0]


    def velocity(self):
        return self.state[1]

    
    def fun(self):
        """
        The right hand side of the system of ODEs. 
        (https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.RK45.html#scipy.integrate.RK45)
        """

        pressure = self.state[3][0]
        mass = self.state[2][0]
        dmdt = - self.nozzle_area * water_density * sqrt(2 * pressure / water_density)

        velocity = self.velocity()
        speed2 = np.sum(velocity * velocity) # scalar
        direction = 1/sqrt(speed2) * velocity # vector

        a_drag = - 1/(2*mass) * self.C_drag * self.A_cross_sectional_area * speed2 * direction
        a_grav = -np.array([0, 9.81])
        # Assumption: It thrusts in the same direction as it's flying. 
        # It points in the direction that it's pointing. 
        # I.e. aerodynamically stable and perfectly responsive.
        a_thrust = 1/mass * 2*self.nozzle_area * pressure * direction
        dvdt = a_drag + a_grav + a_thrust
        # print(f"Mass: {mass}, Thrust force: {mass*a_thrust}, dmdt: {dmdt}")

        dsdt = velocity

        return np.array([dsdt, dvdt, [dmdt, 0], [0, 0]])


    def step(self):
        """
        Step the system one step forward.
        Return: t, position, velocity or None if the water's run out
        """

        # update pressure here.

        mass = self.state[2][0]
        self.state[3][0] = self.pressure_0 * ( (self.volume_0 + (self.mass_0 - mass)/water_density) / self.volume_0 ) ** -self.gamma

        next_state = self.state + self.timestep * self.fun() 
        next_mass = next_state[2][0]
        if next_mass < self.dry_mass:
            return None

        self.state = next_state
        self.t += self.timestep

        return self.t, self.position(), self.velocity()


class BoostWheeler:

    """
    Boost phase equations. From Dan Wheeler: https://www.et.byu.edu/~wheeler/benchtop/thrust.php 

    It's different in that it takes into account the acceleration force of the 
    rocket on the water. The pressure at the bottom of the water (i.e. the 
    nozzle) is also affected by the acceleration.

    INCOMPLETE!
    """

    def __init__(self, position, velocity, t0, water, pressure, dry_mass, volume, C_drag, A_cross_sectional_area, timestep = 0.001):
        # (water is vector in case I want to use scipy integrators)
        self.state = np.array([position, velocity, [water, 0], [pressure, 0]]) 
        self.t = t0
        self.dry_mass = dry_mass
        self.volume = volume
        self.C_drag = C_drag
        self.A_cross_sectional_area = A_cross_sectional_area
        self.timestep = timestep

        raise Exception("INCOMPLETE!")

    
        def fun(self):
            """
            The right hand side of the system of ODEs. 
            (https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.RK45.html#scipy.integrate.RK45)
            """
            velocity = self.velocity()

            speed2 = np.sum(velocity * velocity)
            direction = 1/sqrt(speed2) * velocity

            a_drag = - 1/(2*self.dry_mass) * self.C_drag * self.A_cross_sectional_area * speed2 * direction
            a_grav = -np.array([0, 9.81])
            a_int = 0 # TBD
            a_thrust = 0 # TBD
            dvdt = a_grav + a_thrust + a_drag
    
            dsdt = velocity

            dhdt = 0 # TBD

            dpdt = 0 # TBD
    
            return np.array([dsdt, dvdt, dhdt, dpdt])
    
                


        def step(self):
            """
            Step the system one step forward.
            Return: t, position, velocity, water (l), pressure(kPa)
            """
            pass
