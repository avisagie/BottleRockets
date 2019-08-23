import numpy as np
from numpy import sin, cos, sqrt, pi
from rocket import Stepper, Ballistic, BoosterScienceBits, RocketWithComponents, water_density
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


def sim_3_boosters_bullet(
    radius = 0.045,
    C_drag = 0.3,
    dry_mass = 0.3,
    volume = 3,
    water_l = 3.0 / 3,
    pressure = 10, # relative pressure
    nozzle_radius = 0.0105,
    launch_tube_length = 0.0, # m, must add equations for pushing back on the first stage and the origin to not be fixed. For now make it 0.

    booster_radius = 0.045,
    booster_C_drag = 0.5,
    booster_dry_mass = 0.6,
    booster_volume = 8,
    booster_water_l = 8.0 / 3,
    booster_nozzle_radius = 0.0105,
    booster_launch_tube_length = 1.0, # m

    theta = 45, # degrees
    rail_length = 1.5, # m

    timestep = 0.001
    ):

    """
    Simulate a 3 booster first stage with a small streamlined "bullet". The bullet is a 
    smaller water rocket that starts the moment the first stage stops accelerating.
    """

    stepper = Stepper()

    position = np.array([0, 0.1])
    
    bullet_mass = dry_mass + water_l/1000 * water_density
    print(f'Bullet starting mass: {bullet_mass:0.01f}kg')

    boosters = [BoosterScienceBits( t0=0, 
                                    water=booster_water_l,
                                    pressure=bar2pa(pressure), 
                                    dry_mass=booster_dry_mass + 1.0/3.0*bullet_mass, 
                                    volume=booster_volume, 
                                    C_drag=booster_C_drag, 
                                    A_cross_sectional_area=pi*booster_radius**2, 
                                    nozzle_radius=booster_nozzle_radius, 
                                    launch_tube_length=booster_launch_tube_length,
                                    timestep=timestep) for x in range(3)]

    phase = RocketWithComponents(position, position, 0.001*np.array([cos(deg2rad(theta)), sin(deg2rad(theta))]), 0.0, 
        components=boosters, 
        rail_length=rail_length, 
        timestep=timestep  
    )

    stepper.step(phase)

    position = phase.position()
    velocity = phase.velocity()
    time = phase.t

    bullet = BoosterScienceBits( t0=0, 
                                 water=water_l,
                                 pressure=bar2pa(pressure), 
                                 dry_mass=dry_mass, 
                                 volume=volume, 
                                 C_drag=C_drag, 
                                 A_cross_sectional_area=pi*radius**2, 
                                 nozzle_radius=nozzle_radius, 
                                 launch_tube_length=launch_tube_length,
                                 timestep=timestep )

    phase = RocketWithComponents(position, position, velocity, time, 
        components=[bullet], 
        rail_length=0.0, 
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


def sim_3_stage(
    pressure = 10, # relative pressure

    s1_radius = 0.045,
    s1_C_drag = 0.6,
    s1_dry_mass = 0.3,
    s1_volume = 3,
    s1_water_l = 3.0 / 3,
    s1_nozzle_radius = 0.0105,

    s2_radius = 0.045,
    s2_C_drag = 0.6,
    s2_dry_mass = 0.3,
    s2_volume = 3,
    s2_water_l = 3.0 / 3,
    s2_nozzle_radius = 0.0105,

    s3_radius = 0.045,
    s3_C_drag = 0.6,
    s3_dry_mass = 0.3,
    s3_volume = 3,
    s3_water_l = 3.0 / 3,
    s3_nozzle_radius = 0.0105,

    theta = 45, # degrees
    rail_length = 3, # m

    timestep = 0.001
    ):

    """
    Simulate a 3 booster first stage with a small streamlined "bullet". The bullet is a 
    smaller water rocket that starts the moment the first stage stops accelerating.
    """

    stepper = Stepper()

    position = np.array([0, 0.1])
    velocity = 0.001 * np.array([cos(deg2rad(theta)), sin(deg2rad(theta))])
    time = 0.0

    stage1 = BoosterScienceBits( t0=0, 
                                 water=s1_water_l,
                                 pressure=bar2pa(pressure), 
                                 dry_mass=s1_dry_mass + s2_dry_mass + s2_water_l/1000.0*water_density + s3_dry_mass + s3_water_l/1000.0*water_density, 
                                 volume=s1_volume, 
                                 C_drag=s1_C_drag, 
                                 A_cross_sectional_area=pi*s1_radius**2, 
                                 nozzle_radius=s1_nozzle_radius, 
                                 launch_tube_length=0.5,
                                 timestep=timestep )

    phase = RocketWithComponents(position, position, velocity, time, 
        components=[stage1], 
        rail_length=0.0, 
        timestep=timestep  
    )

    print(f'Starting stage 1 at {phase.t:0.003}s')
    stepper.step(phase)

    position = phase.position()
    velocity = phase.velocity()
    time = phase.t

    stage2 = BoosterScienceBits( t0=0, 
                                 water=s2_water_l,
                                 pressure=bar2pa(pressure), 
                                 dry_mass=s2_dry_mass + s3_dry_mass + s3_water_l/1000.0*water_density, 
                                 volume=s2_volume, 
                                 C_drag=s2_C_drag, 
                                 A_cross_sectional_area=pi*s2_radius**2, 
                                 nozzle_radius=s2_nozzle_radius, 
                                 launch_tube_length=0.0,
                                 timestep=timestep )

    phase = RocketWithComponents(phase.position(), phase.position(), phase.velocity(), phase.t, 
        components=[stage2], 
        rail_length=0.0, 
        timestep=timestep  
    )

    print(f'Starting stage 2 at {phase.t:0.003}s')
    stepper.step(phase)

    stage3 = BoosterScienceBits( t0=0, 
                                 water=s3_water_l,
                                 pressure=bar2pa(pressure), 
                                 dry_mass=s3_dry_mass, 
                                 volume=s3_volume, 
                                 C_drag=s3_C_drag, 
                                 A_cross_sectional_area=pi*s3_radius**2, 
                                 nozzle_radius=s3_nozzle_radius, 
                                 launch_tube_length=0.0,
                                 timestep=timestep )

    phase = RocketWithComponents(phase.position(), phase.position(), phase.velocity(), phase.t, 
        components=[stage3], 
        rail_length=0.0, 
        timestep=timestep  
    )

    print(f'Starting stage 3 at {phase.t:0.003}s')
    stepper.step(phase)

    phase = Ballistic(phase.position(), phase.velocity(), phase.t,
                    dry_mass=s3_dry_mass, 
                    C_drag=s3_C_drag,
                    A_cross_sectional_area=pi * s3_radius**2,
                    timestep=timestep)

    print(f'Ballistic at {phase.t:0.003}s')
    stepper.step(phase)

    return stepper.get_traces()


def sim1(
    radius = 0.045,
    C_drag = 0.3,
    dry_mass = 0.5,
    volume = 8,
    water_l = 8.0 / 3,
    pressure = 6, # relative pressure
    nozzle_radius = 0.0105,
    launch_tube_length = 0.0, # m

    theta = 40, # degrees
    rail_length = 1.5, # m

    extra_frontal_surface = 0.0, # m^2, for things like fins.

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
                                 A_cross_sectional_area=pi*radius**2 + extra_frontal_surface, 
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


def plot_basic(traces):
    time, position, velocity, acceleration = traces
    speed = sqrt(np.sum(velocity * velocity, axis=1))
    accel = sqrt(np.sum(acceleration * acceleration, axis=1)) / 9.81 # in Gs
    max_speed = max(speed)
    max_acceleration = max(sqrt(np.sum(acceleration * acceleration, axis=1)))

    print(position)

    ax1 = pylab.subplot(211)
    ax1.plot(position[:,0], accel, 'b')
    ax1.set_ylabel("Acceleration (g)", color='b')
    ax2 = ax1.twinx()
    ax2.plot(position[:,0], speed, 'r')
    ax2.set_ylabel("Speed (m/s)", color='r')
    ax1.grid()
    ax1.set_title(f"Distance:{np.max(position[:,0]):0.0f}m, flight time:{max(time):0.01f}s")
    ax1.set_xlabel("Horizontal distance (m)")

    ax1 = pylab.subplot(212)
    ax1.plot(position[:, 0], position[:, 1], 'b')
    ax1.set_ylabel("Height (m)", color='b')
    ax1.grid()
    ax1.set_xlabel("Horizontal distance (m)")

    pylab.show()


if __name__ == "__main__":
    traces = sim_3_boosters_bullet()
    plot_basic(traces)