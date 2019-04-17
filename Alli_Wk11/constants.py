from astropy import units as u
import numpy as np

G = 6.67259*10**11 #N-m2/kg2, gravitational constant
#constants from http://www.braeunig.us/space/constant.htm
#equations from http://www.braeunig.us/space/orbmech.htm
class Sun():
    def __init__(self):
        Sun.mass = 	1.9891*10**30 *u.kg
        Sun.radius = 696000 *u.km
    def get_params(self):
        return Sun.mass, Sun.radius

class Earth():
    def __init__(self):
        Earth.mass = 	5.9737*10**24 *u.kg
        Earth.radius = 	6378.137 *u.km
        Earth.ecc = 0.2056
        Earth.i = 7.00 *u.deg
    def get_params(self):
        return Earth.mass, Earth.radius, Earth.ecc, Earth.i
