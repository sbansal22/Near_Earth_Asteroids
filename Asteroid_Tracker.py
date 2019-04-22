"""
Interactive Asteroid Visualization

Author: SPARSH BANSAL
This is a research project in the Nugent Lab, under Dr. Carrie Nugent

Data source: JPL Horizons small space-body dynamics
Orbital Mechanics references: 
A. Orbital Mechanics for Engineering Students- Howard D. Curtis
B. Orbital and Celestial Mechanics - Vinti
C. An Introduction to the Mathematics and Methods of Astrodynamics - Battin
D. Poliastro - Astrodynamics in Python

Research lab partner: Alli Busa

"""

import numpy as np

# Set up matplotlib and use a nicer set of plot parameters

import matplotlib
import matplotlib.pyplot as plt

# Configure Jupyter so figures appear in the notebook

# Configure Jupyter to display the assigned value after an assignment

# import functions from the modsim.py module
from modsim import *

# import math
import numpy as np
import matplotlib.pyplot as plt

#import time
import datetime as dt

from astropy.io import ascii

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

# Defines what the table is - JPL data on asteroids
tbl = ascii.read("ELEMENTS.NUMBR")
#print(tbl)

µ = 1.32712440041 * (10**20)

def display_params(asteroid_name):

    """
    Displays the values of orbital parameters

    Parameters
    ----------

    asteroid_name: Name of the asteroid

    """
    for i in range(0:len(tbl)-1):
        if tbl[i][1] == 'asteroid_name':
            ast_param = tbl[i]
            return ast_param

def assign_params(ast_param):

    """
    Assigns parameters to their respective symbol

    Parameters
    ----------

    ast_param: parameter list from display_params

    """
    
    a = ast_param[3] # Assigning the semi-major axis in AU
    e = ast_param[4] # Assigning the eccentricity 
    Ω = np.deg2rad(ast_param[7]) # Assigning the ascending node in radians
    ω = np.deg2rad(ast_param[6]) # Assigning the argument of perihelion in radians
    inc = np.deg2rad(ast_param[5]) # Assigning the inclination of the orbit w.r.t J2000 ecliptic plane (plane of Earth's orbit around the sun) in radians
    b = a*((1-(e**2))**(1/2)) # Computing the semi-minor axis in AU
    f = a*e # Calculating f, the distance from the center of the orbit to the Sun in AU
    P = (a)**(3/2) # Calculating P, the orbital period
    ma_4_27 = ast_param[8] # Assigning the mean anomaly as on the midnight of April 27, 2019

    return a, e, Ω, ω, inc, b, f, P

def current_epoch():
    """
    Finds the epoch past the Modified Julian Date (midnight on November 17, 1858)

    """
    # How should I make this dynamic so that the user can also input a time/date?

    MJD = dt.datetime(1858,11,17,00,00,00)
    NOW = dt.datetime.now()

    return (NOW-MJD).total_seconds()

def ma_mjd(ma_4_27):

    """
    Computes the mean anomaly of the asteroid at MJD

    Parameters
    ----------

    ma_4_27: mean anomaly as on the midnight of April 27, 2019

    """

    return ma_4_27 - (5063040000)*((µ/(a**3))**(1/2))

def ma_current():

    """
    Computes the current mean anomaly of the asteroid

    """

    return ma_mjd(ma_4_27) + (current_epoch())*((µ/(a**3))**(1/2))

def keplers_equation(ma, E, e):

    """
    Defines the Kepler's Equation for the eccentric anomaly

    Parameters
    ----------

    ma: mean anomaly of the asteroid

    E: eccentric anomaly of the asteroid

    e: eccentricity of the asteroid

    """

    return E - ecc * np.sin(E) - M

def ta_current(ma,e):

    """
    Computes the true anomaly of the asteroid

    Parameters
    ----------

    ma: mean anomaly of the asteroid

    e: eccentricity of the asteroid

    Note
    ----

    Used the angles.py script from poliastro

    """

    return np.deg2rad(M_to_nu(ma, e, delta=1e-2))

def ea_current(ma,e):

    """
    Computes the eccentric anomaly of the asteroid

    Parameters
    ----------

    ma: mean anomaly of the asteroid

    e: eccentricity of the asteroid

    Note
    ----

    Used the angles.py script from poliastro

    """
    
   return np.deg2rad(M_to_E(ma, e)

def distance_current(a, e, E):

    """
    Computes the distance of the asteroid from the Sun.

    Parameters
    ----------

    a: length of the semi-major axis of the orbit of the asteroid

    e: eccentricity of the asteroid's orbit

    E: current eccentric anomaly of the asteroid

    """

    return a * (1 - e*np.cos(E))

def position_current():

    """
    Computes the position vectors to the asteroid for current time in 3D space with respect to the Sun

    """
    p = distance_current(a, e, E))*(np.array([np.cos(ta_current(ma,e)),np.sin(ta_current(ma,e)),0])

    return p

def cartesian_coordinates_current(p):

    """
    Assigns the position vectors to the corresponding cartesian coordinates

    Parameters
    ----------

    p: position vectors of the asteroid at the current time

    """

    ccx = p[0] # x coordinate of the asteroid in 3D space
    ccy = p[1] # y coordinate of the asteroid in 3D space
    ccz = p[2] # z coordinate of the asteroid in 3D space

    return ccx, ccy, ccz

def transformed_cartesian_coordinates_current(ccx, ccy, ccz):

    """
    Transforms the cartesian coordinates of the asteroid's position based on the inclination, ascending node, and argument of perihelion

    Parameters
    ----------

    ccx: x coordinate of the asteroid in 3D space
    
    ccy: y coordinate of the asteroid in 3D space
    
    ccz: z coordinate of the asteroid in 3D space

    """

    tx = ccx*(np.cos(ω)*np.cos(Ω) - np.sin(ω)*np.cos(inc)*np.sin(Ω)) - ccy*(np.sin(ω)*np.cos(Ω) + np.cos(ω)*np.cos(inc)*np.sin(Ω))
    ty = ccx*(np.cos(ω)*np.sin(Ω) + np.sin(ω)*np.cos(inc)*np.cos(Ω)) + ccy*(np.cos(ω)*np.cos(inc)*np.cos(Ω) - np.sin(ω)*np.sin(Ω))
    tz = ccx*(np.sin(ω)*np.sin(inc)) + ccy*(np.cos(ω)*np.sin(inc))

    return tx, ty, tz

def plot(tx, ty, tz):

    """
    Plots the current location of the asteroid with respect to the Sun

    Parameters
    ----------

    tx: projected x coordinate

    ty: projected y coordinate

    tz: projected z coordinate

    """

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot(tx, ty, tz, label='Current location in space')
    ax.legend()
    plt.show()

    # How can I make this plot update with time?



# -----------------------------------------------------------------------------
# Run if called from the command line
# -----------------------------------------------------------------------------
if __name__ == "__main__":

    display_params(Ceres)

