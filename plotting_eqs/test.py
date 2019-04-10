import constants
import equations
#import plotting2
import sys
import numpy as np
sys.path.insert(0, '/home/alli/NugentResearch')
from wk5 import ImportFromTxt
from IPython import sys_info
#from IPython.core.interactiveshell import InteractiveShell
#InteractiveShell.showtraceback()
print(sys_info())
#Import Eros Orbital Elements
i, ecc, a, M, long, peri = ImportFromTxt('/home/alli/NugentResearch/eros.txt')

nu = equations.calculateEllipticalTrueAnomaly(M, ecc)
print('True Anomaly of Given Asteroid: '+str(nu))
pos = equations.calculatePosition(nu, ecc, a)
angle = equations.calculatePeriapsisAngle(nu, ecc, a)


ex,ey,ez = equations.convertToCartesian2(angle, i, long, pos)

X,Y,Z = equations.createMotion(ecc, a, M, long)
print(X,Y,Z)

#plotting2.plotPoint(ex,ey,ez)
