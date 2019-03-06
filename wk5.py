##############
#Imports I used in Week 3 , Importing / Using MP data
import pandas as pd
import mpcorbget
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import datetime # for getting current time

#imports for plotting
import orbital
import poliastro
from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import
from poliastro.plotting.static import StaticOrbitPlotter
from astropy import units as u
###############################
#Classes and Functions

def returnOrbitalElements(asteroid):
    """
    Takes the orbital elements dictionary that MPCORBget creates, and returns each individually
    """
    epoch = asteroid.orbEl.get('epoch')
    mean_anomaly = asteroid.orbEl.get('ME')
    argument_of_perihelion = asteroid.orbEl.get('w')
    longitude_of_ascendingnode = asteroid.orbEl.get('0')
    inclination = asteroid.orbEl.get('i')
    eccentricity = asteroid.orbEl.get('e')
    semimajor_axis = asteroid.orbEl.get('a')
    mean_daily_motion = asteroid.orbEl.get('n')
    return epoch, mean_anomaly, argument_of_perihelion, longitude_of_ascendingnode, inclination, eccentricity, semimajor_axis,mean_daily_motion


def whatTimeIsIt():
    """
    returns current time and date, in order to plot orbitals at certain position
    """
    currentDT = datetime.datetime.now()
    year = int(currentDT.year)
    month = int(currentDT.month)
    day = int(currentDT.day)
    hour = int(currentDT.hour)
    minute = int(currentDT.minute)
    second = int(currentDT.second)

    return year, month, day, hour, minute, second


def ImportFromTxt(txtfile):
    """
    Reads data from eros.txt file
    """
    df = pd.read_csv(txtfile, sep="  ", header=None)
    #print(df.to_string())
    #epoch = dateUnpack(df.iloc[0,4])
    inclination = float(df.iloc[0,7])
    eccentricity = float(df.iloc[0,8])
    semimajor_axis = float(df.iloc[0,10])
    argument_of_perihelion = float(df.iloc[0,5])
    mean_anomaly = float(df.iloc[0,4]) #arugment _of_perihelion
    longitude_of_ascendingnode = float(df.iloc[0,6])
    return inclination, eccentricity, semimajor_axis, mean_anomaly, longitude_of_ascendingnode, argument_of_perihelion


def calculateTrueAnomaly(mean_anomaly, eccentricity):
    """
    Calculates true anomaly from mean anomaly. Equation taken from wikipedia.
    """
    tru_anom = mean_anomaly + (2*eccentricity - (1/4)*eccentricity**3)*np.sin(mean_anomaly) +  (5/4)*np.sin(2*mean_anomaly)*eccentricity**(2) + (13/12)*np.sin(3*mean_anomaly)*eccentricity**3 *3
    return tru_anom


def dateUnpack(self, packed):
        """A package from mpcorgbet"""
        yearcode = {"I":"18","J":"19","K":"20"}
        daycode = "123456789ABCDEFGHIJKLMNOPQRSTUV"
        year = yearcode[packed[0]]+packed[1:3]
        month = daycode.index(packed[3])+1
        day = daycode.index(packed[4])+1
        return "%s/%s/%s" % (month, day, year)
###############################
#Main Execution

#Commands for getting mpcorgbet data
#first = mpcorbget.MPCORB("162269")
#epoch, mean_anomaly, argument_of_perihelion, longitude_of_ascendingnode, inclination, eccentricity, semimajor_axis,mean_daily_motion = returnOrbitalElements(first)

#Commands for getting eros data from txt file
i, e, a, ME, long, peri = ImportFromTxt('/home/alli/NugentResearch/eros.txt')
tru_anom = calculateTrueAnomaly(ME, e)
print(i, e, a, ME, long, peri)
#Creating poliastro orbit object
ss = Orbit.from_classical(Earth, a* u.AU, e* u.one, i* u.deg, long * u.deg, peri * u.deg, tru_anom * u.deg)
#Plotting Eros
plotter = StaticOrbitPlotter() #creating plotter class
plotter.plot(ss, label="Eros") #plotting Eros
plt.show()
