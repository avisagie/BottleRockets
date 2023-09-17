
from typing import Callable
import numpy as np
def simple_bottle(nozzle_radius,body_radius) -> Callable:
    """
    Return simple bottle profile:
    
    ▲  │                 │
    │  │                 │
    │  │                 │
    │  │                 │
    │  │                 │
    Z  │                 │
    │  │                 │
    │  │                 │
    │  │                 │
    │  │                 │
    │  │                 │
    │  │                 │
    │  │                 │
    │  │                 │
    0  └──────     ──────┘
    """
    nozzle_area = np.pi*nozzle_radius**2
    body_area = np.pi*body_radius**2
    def simple_profile(z):
        if z == 0:
            return nozzle_area
        if z > 0 :
            return body_area
        
    return simple_profile

def infinite_bottle(nozzle_radius):
    """
    Assumes the body radius is much larger than the nozzle radius.
    This means there is no water flow to the bottle except
    at the nozzle.

    <-- ──────     ────── -->
    """
    nozzle_area = np.pi*nozzle_radius**2
    def naive_profile(z):
        if z == 0:
            return nozzle_area
        if z > 0:
            return 1e6
    
    return naive_profile

def triangular_bottle(nozzle_radius,body_radius,taper_angle):
    """
    Return simple bottle profile:
    
    ▲  │                 │
    │  │                 │
    │  │                 │
    │  │                 │
    │  │                 │
    Z  │                 │
    │  │                 │
    │  │                 │
    │   \               /  
    │    \             /  
    │     \           /   
    │      \     |   /    
    │       \    |a /     
    │        \   | /      
    0         \  |/

    a: taper_angle
    """
    nozzle_area = np.pi*nozzle_radius**2
    body_area = np.pi*body_radius**2

    double_tan_a = 2.0*np.tan(taper_angle)
    taper_distance = (body_radius - nozzle_radius)/(double_tan_a)
    def taper_profile(z):
        if z == 0:
            return nozzle_area
        elif z > taper_distance:
            return body_area
        else:
            radius = nozzle_radius + z*double_tan_a
            taper_area = np.pi*radius**2
            return taper_area
        
    return taper_profile


def get_bottle_helper(type_of_bottle : str,nozzle_radius,body_radius) -> Callable:
    if type_of_bottle == "simple":
        return simple_bottle(nozzle_radius,body_radius)
    elif type_of_bottle == "infinite":
        return infinite_bottle(nozzle_radius)
    elif type_of_bottle == "taper_60":
        return triangular_bottle(nozzle_radius,body_radius,np.deg2rad(90.0 - 60.0))
    elif type_of_bottle == "taper_45":
        return triangular_bottle(nozzle_radius,body_radius,np.deg2rad(45.0))

