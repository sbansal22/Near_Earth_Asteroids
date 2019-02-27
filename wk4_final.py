##############
#Imports I used in Week 3 , Importing / Using MP data
import pandas as pd
import mpcorbget
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import datetime # for getting current time

#imports suggested from Carrie
import orbital

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
    df = pd.read_csv(txtfile, sep="  ", header=None)
    print(float(df.iloc[0,9]),float(df.iloc[0,8]))
    #epoch = dateUnpack(df[20:25].strip())
    #mean_anomaly = float(df[26:35].str.strip())
    inclination = float(df.iloc[0,8])
    eccentricity = float(df.iloc[0,9])
    semimajor_axis = float(df.iloc[0,11])
    return inclination, eccentricity, semimajor_axis

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
i, e, a = ImportFromTxt('/home/alli/NugentResearch/eros.txt')
#first_to_plot = orbital.KeplerianElements(a=semimajor_axis, e=eccentricity, i=inclination, body=orbital.earth)
first_to_plot = orbital.KeplerianElements(a=a, e=e, i=i, body=orbital.earth)

#Imports to try to get animation to work
import matplotlib.animation as animation
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'

#Dummy Asteriod
orbit1 = orbital.KeplerianElements.with_period(90 * 60, body=orbital.earth)
orbital.plot(orbit1, title='Orbit Test')

#My Asteroid in 2D plot
orbital.plot(first_to_plot, title='MyOrbit 1' )
#line1_anim.save('MyOrbit-1.mp4', writer=writer)

#My Asteroid in 3D plot
line2_anim = orbital.plot3d(first_to_plot, title='MyOrbit 2', animate=True )
#line2_anim.save('MyOrbit-2.mp4', writer=writer)

plt.show()
