import numpy as np
from numpy import sin, cos, sqrt, pi
# import scipy.integrate
from typing import Tuple,Callable
from abc import ABC,abstractmethod
from collections import namedtuple
from physics import pressure_after_adiabatic_expansion
import scipy.integrate as spi
from scipy.optimize import brentq
from bottle_shapes import get_bottle_helper
from icecream import ic

# Assume it's constant for the flight and use dry air, which is more dense than humid air.
# https://en.wikipedia.org/wiki/Density_of_air
# Can be refined later
air_density = 1.2 # kg / m^3 

# pretty close for our temperatures and height above sea level
water_density = 1000.0 # kg / m^3

verbose = False

class Phase(ABC):
    @abstractmethod
    def step(self) -> Tuple[float,np.ndarray,np.ndarray]:
        """
        Step the system one step forward.
        Return: t, position, velocity
        """
        pass

    @abstractmethod
    def position(self) -> np.ndarray:
        pass

    @property
    def acceleration(self) -> np.ndarray:
        return self._acceleration
    
    @acceleration.setter
    def acceleration(self,value):
        self._acceleration = value

class WaterThruster():
    """
    Computes thrust and inertial forces due to the movement and expulsion of water.
    See "Water Motion from the Bernoulli Equation" from https://www.et.byu.edu/~wheeler/benchtop/pix/thrust_eqns.pdf
    Notation:
      z : distance from nozzle, 
      A(z) : cross-sectional area at z
      z_H : current "height" of the water. (distance of boundary to nozzle)
      u_out : speed of ejected water relative to the nozzle 
    """
    def __init__(self, nozzle_profile : Callable, starting_water_volume : float):
        """
        nozzle_profile: A callable A(z) that returns the cross sectional area in m^2 of the water
        vessel at distance z from the outlet/nozzle. This is important for calculating the
        flow of water at distance z from the nozzle.

        starting_volume: A float providing the amount of water in cubic meters.
        Used together with nozzle_profile (A(z)) to compute the starting height of the boundary
        area z_start.
        """

        self.A = nozzle_profile
        self.water_volume_0 = starting_water_volume
        self.z_H = self._distance_from_nozzle(nozzle_profile,starting_water_volume)
        self.A_out = nozzle_profile(0)

        self.u_out = 0

    @staticmethod
    def _distance_from_nozzle(A : Callable,starting_volume : float):

        def volume_difference(z):
            volume,_ = spi.quad(A,0,z)
            return volume - starting_volume
        
        z_0 = brentq(volume_difference, 0, 10)
        return z_0
    
    def _average_nozzle_to_body_ratio(self, z_H):
        "B(H)"
        ratio = lambda z: self.A_out/self.A(z)
        area,error = spi.quad(ratio,0,z_H)
        return area
    
    def _squared_boundary_ratio(self,z_H):
        "C(H)"
        ratio = self.A_out/self.A(z_H)
        return (ratio**2 - 1.0) / 2.0
    
    def step(self,P_gauge : float,d_t : float, acceleration : np.ndarray) -> Tuple[float,float]:
        """
        Returns thrust and mass of water ejected for one time step.
        Requires current gauge pressure (absolute - atmospheric) and acceleration of the thruster
        """
        if self.z_H <= 0:
            return 0,0

        B_H = self._average_nozzle_to_body_ratio(self.z_H)
        C_H = self._squared_boundary_ratio (self.z_H)
        a =  np.linalg.norm(acceleration)

        if numerically_stable := B_H > d_t*100:
            du_dt = - (C_H*self.u_out**2 + P_gauge/water_density + a*self.z_H)/B_H 

            du = du_dt*d_t
            self.u_out += du
            ic(numerically_stable)

        else:
            self.u_out = -(-(P_gauge/water_density + a*self.z_H)/C_H)**(0.5)
            du_dt = 0

        F_thrust = water_density * self.A_out * self.u_out**2
        F_internal = -water_density*self.A_out*(self.z_H*du_dt + self.A_out/self.A(self.z_H)*self.u_out**2)

        F_tot = F_thrust + F_internal
        # compute how far the water boundary moved
        A_H = self.A(self.z_H)
        d_H_dt = (self.A_out * self.u_out)/A_H
        d_H = d_H_dt*d_t
        self.z_H += d_H

        ejected_mass = abs(d_H*A_H*water_density)
        ic(du_dt,self.u_out,F_internal)

        return F_tot,ejected_mass
    
class NaiveWaterThruster:
    def __init__(self,nozzle_radius : float):
        self.A_out = np.pi*nozzle_radius**2

    def step(self,P_gauge : float,d_t : float, acceleration : np.ndarray):
        F_tot = 2.0*P_gauge*self.A_out

        dm_dt = self.A_out * water_density * sqrt(2 * P_gauge / water_density)
        ejected_mass = dm_dt*d_t
        return F_tot,ejected_mass


class Ballistic(Phase):
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


def always_happy(*args):
    return True

class BoosterScienceBits():

    """
    Boost phase from science bits: http://www.sciencebits.com/RocketEqs 
    Use with RocketWithComponents
    """

    def __init__(self, t0, 
                 water, # in liter
                 pressure, # in kPa relative to outside pressure
                 dry_mass, volume, # in liter, total volume
                 C_drag, A_cross_sectional_area, nozzle_radius, # in meters
                 launch_tube_length = 0.0, # m. Assumes that the launch tube fits snuggly and has no friction
                 removable = True,
                 timestep = 0.001,
                 bottle_shape = "simple"):
        # (mass and pressure are vectors in case I want to use scipy integrators)
        mass = dry_mass + water/1000 * water_density # 1l = 0.001m^3
        self.state = np.array([[mass], [pressure]])
        self.t = t0
        self.dry_mass = dry_mass
        self.mass_0 = mass
        self.nozzle_radius = nozzle_radius
        self.nozzle_area = pi*self.nozzle_radius**2
        self.volume = volume / 1000 # total volume in m^3, 1l = 0.001m^3
        self.volume_0 = self.volume - water/1000 - self.nozzle_area*launch_tube_length # air volume at t0
        self.pressure_0 = pressure
        self.C_drag = C_drag
        self.A_cross_sectional_area = A_cross_sectional_area
        self.timestep = timestep

        self.launch_tube_length = max(0.0, launch_tube_length)

        self.__removable = removable

        if bottle_shape == "naive":
            self.thruster = NaiveWaterThruster(nozzle_radius)
        else:
            bottle_radius = (A_cross_sectional_area/np.pi)**0.5
            nozzle_profile = get_bottle_helper(bottle_shape,nozzle_radius,bottle_radius)
            self.thruster = WaterThruster(nozzle_profile=nozzle_profile,starting_water_volume=water/1000)

        self.f_thrust = 0


    def position(self):
        return self.state[0]


    def velocity(self):
        return self.state[1]


    def removable(self):
        return self.__removable


    def F_drag(self, speed2):
        return - 0.5 * self.C_drag * air_density * self.A_cross_sectional_area * speed2


    def F_thrust(self):
        return self.f_thrust

    
    def mass(self):
        mass = self.state[0, 0]
        return mass

    def step(self,acceleration,distance_from_origin) -> bool:
        """
        Step the system one step forward.
        Return: True if the booster is done, else False
        """
        in_launch_tube_phase = distance_from_origin < self.launch_tube_length
        pressure = self.state[1, 0]
        if in_launch_tube_phase:
            self.f_thrust = self.nozzle_area * pressure # https://www.et.byu.edu/~wheeler/benchtop/pix/thrust_eqns.pdf Section III
            d_mass = 0
        else:
            F,d_mass = self.thruster.step(pressure,self.timestep,acceleration)
            self.state[0, 0] -= d_mass
            self.f_thrust = F
        
        mass = self.mass()
        water_volume_lost = (self.mass_0 - mass)/water_density
        launch_pipe_volume_lost = min(self.launch_tube_length, distance_from_origin) * self.nozzle_area
        self.state[1, 0] = pressure_after_adiabatic_expansion(
            starting_pressure=self.pressure_0,
            starting_gas_volume=self.volume_0,
            current_gas_volume=self.volume_0 + water_volume_lost + launch_pipe_volume_lost,
        )

        done = True if mass < self.dry_mass else False
        return done

class RocketWithComponents(Phase):

    """
    Boost phase that sheds components when they are spent (step returns None)
    """

    def __init__(self, position, origin, velocity, t0,
                 components,
                 rail_length = 0.0,
                 validate = always_happy, # tell if this system holds together
                 timestep = 0.001):
        # (mass and pressure are vectors in case I want to use scipy integrators)
        self.state = np.array([position, velocity]) 
        self.t = t0
        self.timestep = timestep
        self.acceleration = 0.0
        self.components = components
        self.rail_length = rail_length
        self.origin = origin

        self.validate = validate

    @property
    def components(self) -> BoosterScienceBits:
        return self._components

    @components.setter
    def components(self,values):
        if not isinstance(values,list):
            values = [values]

        if not all([isinstance(v,BoosterScienceBits) for v in values]):
            raise ValueError("Currently only components of type BoosterScienceBits are allowed")
        
        self._components = values

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

        distance_from_origin = np.linalg.norm(self.origin - self.position())

        # while it's on the rail, just get the component of gravity backwards along the rail
        if distance_from_origin > self.rail_length:
            a_grav = -np.array([0, 9.81])
        else:
            a_grav = np.dot(-np.array([0, 9.81]), direction) * direction

        # Assumption: It thrusts in the same direction as it's flying. 
        # It points in the direction that it's pointing. 
        # I.e. aerodynamically stable and perfectly responsive.
        F_thrust = self.F_thrust()

        mass = self.mass()

        # print(f'{self.t:0.03f}s: {F_thrust}, {F_drag}, {mass}kg')

        a_thrust = 1/mass * F_thrust * direction
        a_drag = 1/mass * F_drag * direction
        dvdt = a_drag + a_grav + a_thrust
        self.acceleration = dvdt

        dsdt = velocity

        if not self.validate(speed2):
            raise Exception("It broke apart ")

        return np.array([dsdt, dvdt])


    def step(self):
        """
        Step the system one step forward.
        Return: t, position, velocity or None if the water's run out
        """
        done = []
        for c in self.components:
            dist = np.linalg.norm(self.origin - self.position())
            if c.step(self.acceleration,distance_from_origin = dist):
                done.append(c)

        for d in done:            
            self.components.remove(d)
            if not d.removable():
                raise Exception("Non removable part wants to fall off!!!")
            if verbose: 
                print(f"Removing {d} from a list of {len(self.components)} components")

        if not self.components:
            return None
        
        self.state = self.state + self.timestep * self.fun() 
        self.t += self.timestep

        return self.t, self.position(), self.velocity()

Traces = namedtuple("Trace",["time","position","velocity","acceleration"])
class Stepper:

    def __init__(self, print_interval=0.025):
        self.t_time = []
        self.t_position = []
        self.t_velocity = []
        self.t_acceleration = []
        self.print_interval = print_interval
        self.cur_interval: int = 0
    

    def step(self, phase : Phase):
        """
        Step until the rocket hits the floor or the phase returns None.
        """

        altitude = phase.position()[1]

        while altitude > 0.0:

            ret = phase.step()
            if ret is None:
                # water's run out, or it's over somehow.
                break

            time, position, velocity = ret
            altitude = position[1]

            interval = int(time/self.print_interval)
            if interval > self.cur_interval:
                if verbose: print(f'{time:0.03f}s: position:{position}, velocity:{velocity}')    
                self.cur_interval = interval

            self.t_time.append(time)
            self.t_position.append(position)
            self.t_velocity.append(velocity)
            self.t_acceleration.append(phase.acceleration)

        self.count = 0

    def get_traces(self) -> Traces:
        trace = Traces(
            np.array(self.t_time),
            np.array(self.t_position),
            np.array(self.t_velocity),
            np.array(self.t_acceleration)
        )
        return trace 
