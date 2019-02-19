

################################################################
#Imports which I find useful

import pandas as pd
import mpcorbget
################################################################
#Imports for Astropy data visualization

import numpy as np
# Set up matplotlib and use a nicer set of plot parameters
#%config InlineBackend.rc = {} - used in jupyter notebook, makes it pretty
import matplotlib
#matplotlib.rc_file("../../templates/matplotlibrc")
import matplotlib.pyplot as plt
#%matplotlib inline
from astropy.io import ascii #for reading files from a csv file
import astropy.coordinates as coord
import astropy.units as u

################################################################

#Block 1  -  import dat file
"""
#Method 1 for importing : Pandas
#Success : Not yet
#importing from pandas was causing issues b/c of different column numbers
near_asteroids = df = pd.read_csv('/home/alli/NugentResearch/nea_extended.dat', sep='\s+', header=None, skiprows=1)
print(near_asteroids.head())
"""

#Method 2 for importing : mpcorbget
#Success : Yes
#trying to import using mpcorbget
first = mpcorbget.MPCORB("162269")
print(first.objdata) #will do getMPC function, returns row of minor planet center data
print(first.geocentric('1998/08/20/00')) #returns geocentric coordinates for that specific time


"""
#Method 3 for importing : ascii
#Success : Not yet
data = ascii.read("/home/alli/NugentResearch/nea_extended.dat", header_start=1, data_start=2, guess=False, fast_reader=False)
print(data)
"""


#####################################################################
#Block 1  -  Visualize data using data viz tutorial

"""
Try 1: using ascii
first_row = data[0] # get the first (0th) row
first_row.xephem = self.rdes + ",e," + str(self.orbEl["i"]) + "," + str(self.orbEl["O"]) + "," + str(self.orbEl["w"]) + "," + str(self.orbEl["a"]) + "," + str(self.orbEl["n"]) + "," + str(self.orbEl["e"]) + "," + str(self.orbEl["ME"]) + "," + self.orbEl["epoch"] + ",2000,H" + str(self.H) + "," + str(self.G)
first_row.target = ephem.readdb(first_row.xephem)
print(first_row.target.a_ra, first_row.target.a_dec)
"""
ra = first.target.a_ra
#ra = ra.wrap_at(180*u.degree)
dec = first.target.a_dec
print(ra, dec)



fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection="mollweide")
ax.scatter(np.radians(ra), np.radians(dec))
ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'])
ax.grid(True)
plt.show()
