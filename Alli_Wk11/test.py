import constants
import equations
import sys
import numpy as np
sys.path.insert(0, '/home/alli/NugentResearch')
from importdata import ImportFromTxt
from IPython import sys_info
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from IPython import get_ipython


#Import Eros Orbital Elements
i, ecc, a, M, long, peri = ImportFromTxt('/home/alli/NugentResearch/eros.txt')
print('i: '+str(i))
print('ecc: '+str(ecc))
print('a: '+str(a))
print('M: '+str(M))
print('long: '+str(long))

#Do calculations and print out results
inc, M, long, peri = equations.convert2rad(i, M, long, peri)
nu = equations.calculateEllipticalTrueAnomaly(M, ecc)
#print('True Anomaly of Given Asteroid: '+str(nu))
pos = equations.calculatePosition(nu, ecc, a)
#print('position: '+str(pos))
angle = equations.calculatePeriapsisAngle(nu, ecc, a)
#print('angle: '+str(angle))
ex,ey,ez = equations.convertToCartesian2(angle, i, long, pos)
#print('x: '+str(ex))
#print('y: '+str(ey))
#print('z: '+str(ez))

# Apply it to all mean anomalies
#This calculates all the x,y,z coordinates
X,Y,Z = equations.createMotion(ecc, a, M, long, inc)
#print(X,Y,Z) #printing all the x,y,z coordinates
#
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

#for i in zip(X,Y,Z):
    #print('printing')
    #print(i[0], i[1], i[2])
for i in zip(X,Y,Z):
    ax.scatter(i[0], i[1], i[2])

ax.set_xlim(-1,1)
ax.set_ylim(-1,1)
ax.set_zlim(-1,1)
plt.show()

#get_ipython().run_line_magic('matplotlib', 'qt')
