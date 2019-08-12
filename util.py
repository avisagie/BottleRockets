
from numpy import pi

def kmh2ms(kmh):
    return kmh  * (1000/3600)

def ms2kmh(ms):
    return ms  / (1000*3600)

def rad2deg(rad):
    return rad/(2*pi)*360

def deg2rad(deg):
    return deg / 360 * 2 * pi

def bar2pa(bar):
    return bar * 100000.0
