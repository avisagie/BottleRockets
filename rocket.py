import numpy as np
from numpy import sin, cos, sqrt, pi
import scipy.integrate
from typing import Callable, List, Dict

# Assume it's constant for the flight and use dry air, which is more dense than humid air.
# https://en.wikipedia.org/wiki/Density_of_air
# Can be refined later
air_density = 1.2 # kg / m^3 

# pretty close for our temperatures and height above sea level
water_density = 1000.0 # kg / m^3

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
        self.acceleration = 0.0


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
        a_drag = - 1/(2*self.dry_mass) * self.C_drag * air_density * self.A_cross_sectional_area * speed2 * direction

        a_grav = -np.array([0, 9.81])

        dvdt = a_grav + a_drag
        self.acceleration = dvdt

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
                 rail_length = 0.0, # launch rails that keep it upright. gravity then work 
                 pipe_length = 0.0, # launch pipe for a piston effect
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
        self.rail_length = rail_length

        self.acceleration = 0.0


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
        speed2 = np.sum(velocity * velocity) # scalar
        direction = 1/sqrt(speed2) * velocity # vector

        F_drag = - 0.5 * self.C_drag * air_density * self.A_cross_sectional_area * speed2

        pressure = self.state[3][0]
        mass = self.state[2][0]
        dmdt = - self.nozzle_area * water_density * sqrt(2 * pressure / water_density)

        # Assumption: It thrusts in the same direction as it's flying. 
        # It points in the direction that it's pointing. 
        # I.e. aerodynamically stable and perfectly responsive.
        F_thrust = 2*self.nozzle_area * pressure

        # while it's on the rail, just get the component of gravity backwards along the rail
        distance_from_origin = sqrt(sum(self.position() * self.position()))
        if distance_from_origin > self.rail_length:
            a_grav = -np.array([0, 9.81])
        else:
            a_grav = np.dot(-np.array([0, 9.81]), direction) * direction

        a_thrust = 1/mass * F_thrust * direction
        a_drag = 1/mass * F_drag * -direction
        dvdt = a_drag + a_grav + a_thrust
        self.acceleration = dvdt
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


class RocketWithComponents:

    """
    Boost phase from science bits: http://www.sciencebits.com/RocketEqs 
    """

    def __init__(self, position, velocity, t0,
                 components,
                 rail_length = 0.0,
                 validate = Callable[[], bool], # tell if this system holds together
                 timestep = 0.001):
        # (mass and pressure are vectors in case I want to use scipy integrators)
        self.state = np.array([position, velocity]) 
        self.t = t0
        self.timestep = timestep
        self.acceleration = 0.0
        self.components = components
        self.rail_length = rail_length

        self.validate = validate


    def position(self):
        return self.state[0]


    def velocity(self):
        return self.state[1]

    
    def F_drag(self, speed2):
        return sum(c.F_drag(speed2) for c in self.components)


    def F_thrust(self):
        return sum(c.F_thrust() for c in self.components)


    def mass(self):
        return sum(c.mass() for c in self.components)


    def fun(self):
        """
        The right hand side of the system of ODEs. 
        (https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.RK45.html#scipy.integrate.RK45)
        """

        velocity = self.velocity()
        speed2 = np.sum(velocity * velocity) # scalar
        direction = 1/sqrt(speed2) * velocity # vector

        F_drag = self.F_drag(speed2)

        # Assumption: It thrusts in the same direction as it's flying. 
        # It points in the direction that it's pointing. 
        # I.e. aerodynamically stable and perfectly responsive.
        F_thrust = self.F_thrust()

        # while it's on the rail, just get the component of gravity backwards along the rail
        distance_from_origin = sqrt(sum(self.position() * self.position()))
        if distance_from_origin > self.rail_length:
            a_grav = -np.array([0, 9.81])
        else:
            a_grav = np.dot(-np.array([0, 9.81]), direction) * direction

        mass = self.mass()

        a_thrust = 1/mass * F_thrust * direction
        a_drag = 1/mass * F_drag * -direction
        dvdt = a_drag + a_grav + a_thrust
        self.acceleration = dvdt

        dsdt = velocity

        return np.array([dsdt, dvdt])


    def step(self):
        """
        Step the system one step forward.
        Return: t, position, velocity or None if the water's run out
        """

        self.state = self.state + self.timestep * self.fun() 
        self.t += self.timestep

        done = []
        for c in self.components:
            step = c.step()
            if step is None:
                done.append(c)

        for d in done:
            print(f"    removing {d} from a list of {len(self.components)} components")
            self.components.remove(d)

        if not self.components:
            return None

        return self.t, self.position(), self.velocity()


class BoosterScienceBits:

    """
    Boost phase from science bits: http://www.sciencebits.com/RocketEqs 
    Use with RocketWithComponents
    """

    def __init__(self, t0, 
                 water, # in liter
                 pressure, # in kPa relative to outside pressure
                 dry_mass, volume, # in liter, total volume
                 C_drag, A_cross_sectional_area, nozzle_radius, # in meters
                 timestep = 0.001):
        # (mass and pressure are vectors in case I want to use scipy integrators)
        mass = dry_mass + water/1000 * water_density # 1l = 0.001m^3
        self.state = np.array([[mass], [pressure]])
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

    
    def F_drag(self, speed2):
        return - 0.5 * self.C_drag * air_density * self.A_cross_sectional_area * speed2


    def F_thrust(self):
        pressure = self.state[1, 0]
        return 2*self.nozzle_area * pressure

    
    def mass(self):
        mass = self.state[0, 0]
        return mass


    def fun(self):
        """
        The right hand side of the system of ODEs. 
        (https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.RK45.html#scipy.integrate.RK45)
        """

        pressure = self.state[1, 0]
        mass = self.state[0, 0]
        dmdt = - self.nozzle_area * water_density * sqrt(2 * pressure / water_density)

        return np.array([[dmdt], [0]])


    def step(self):
        """
        Step the system one step forward.
        Return: t, position, velocity or None if the water's run out
        """

        # update pressure here.

        mass = self.mass()
        self.state[1, 0] = self.pressure_0 * ( (self.volume_0 + (self.mass_0 - mass)/water_density) / self.volume_0 ) ** -self.gamma

        next_state = self.state + self.timestep * self.fun()

        next_mass = next_state[0, 0]
        if next_mass < self.dry_mass:
            return None

        self.state = next_state
        self.t += self.timestep

        return self.t


class Stepper:

    def __init__(self, print_interval=0.01):
        self.t_time = []
        self.t_position = []
        self.t_velocity = []
        self.t_acceleration = []
        self.print_interval = print_interval
        self.cur_interval: int = 0
    

    def step(self, phase):
        """
        Step until the rocket hits the floor or the phase returns None.
        """

        altitude = phase.position()[1]

        while altitude > 0.0:

            ret = phase.step()
            if ret is None:
                # water's run out
                break

            time, position, velocity = ret
            altitude = position[1]

            interval = int(time/self.print_interval)
            if interval > self.cur_interval:
                print(f'{time:0.04f}: {position}, {velocity}')    
                self.cur_interval = interval

            self.t_time.append(time)
            self.t_position.append(position)
            self.t_velocity.append(velocity)
            self.t_acceleration.append(phase.acceleration)

        self.count = 0

    def get_traces(self):
        return np.array(self.t_time), np.array(self.t_position), np.array(self.t_velocity), np.array(self.t_acceleration)
