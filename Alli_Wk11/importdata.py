##############
#Imports I used in Week 3 , Importing / Using MP data
import pandas as pd
import mpcorbget
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import datetime # for getting current time

## I will add more functions later
#Now it just imports the text file I gave
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
