
##############
#Imports I used in Week 3 , Importing / Using MP data
import pandas as pd
import mpcorbget
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import ascii #for reading files from a csv file
import astropy.coordinates as coord
import astropy.units as u
##############
#Imports for Visualization
from mpl_toolkits import mplot3d
##############
#Importing Data
first = mpcorbget.MPCORB("162269")
print(first.objdata) #will do getMPC function, returns row of minor planet center data
print(first.geocentric('1998/08/20/00')) #returns geocentric coordinates for that specific time
print(first.target)
#Getting RA & dec
ra = first.target.a_ra
dec = first.target.a_dec
#Converting to Cartesian, based on stackoverflow answer
def convert_to_cartesian(ra,dec):
    A = (ra[1] * 15) + (ra[2] * 0.25) + (ra[3] * 0.004166)
#ra.to(u.hourangle)
np.savetxt('162269_data.csv', ("162269", str(ra), str(dec)), delimiter=',')
