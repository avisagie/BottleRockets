
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