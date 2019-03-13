
import imageio
import glob, os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import fsolve
import math
import numpy as np
from datetime import datetime

#-----------------------------------------------
# System Conditions & Initial Orbital Elements
#-----------------------------------------------

dt = 8        # Time step is 10 minutes in seconds
G = 6.67e-11   # Universal Gravitational Constant
M = 5.972e24   # Mass of Earth in kg

# Create initial position and velocity vectors
r = np.array([6.771e6,0,0])                             # Initial position east
v = np.array([0,8.2e3,0])                              # Initial velocity north


#--------------------------------------------
# Main
#
# Calls program functions
#
# Input: none
# Output: none (images save to file system)
#--------------------------------------------
def main():
    # Start timing the script
    startTime = datetime.now()
    
    # Calculate the system coordinates
    newtonian_cors = newtonian_solve()
    orbital_elem_cors = orbital_elements_solve()
    
    # Graph each set of coordinates and save it as the frame
    # of a video (saves .png files to the filesystem)
    #make_frames([newtonian_cors,orbital_elem_cors])
    make_frames(newtonian_cors)

    # Turn the .png frames into an animated .gif file
    make_gif()
    
    err_list = calc_error(newtonian_cors, orbital_elem_cors)
    
    # Finish timing the script
    print("Script duration: " + str(datetime.now() - startTime))
    
    print("Returning err_data.")
    print(sum(err_list)/len(err_list))
    return [sum(err_list)/len(err_list),err_list]


#--------------------------------------------
# Orbital Elements Solve
#
# Solve positions using equations that derive position
# as a function of the system orbital elements
#
# Input: dt, G, M, r_0, v_0
# Output: [list(x cors), list(y cors), list(z cors)]
#--------------------------------------------
def orbital_elements_solve(dt=dt,G=G,M=M,r=r,v=v):
    print("Performing an orbital elements orbit generation...")
    # Initialize empty lists to hold orbital positions 
    xorb = [] 
    yorb = []          
    zorb = []

    # Initialize loop counter and theta sentinel variable
    i=0
    theta = 0
    
    # Loop to calculate orbital positions
    while 2*math.pi-theta > 0.01:
        # Recalculate orbital elements
        u=G*M
        h = np.cross(r,v)
        e = float(np.linalg.norm(np.cross(v,h)/u-r/np.linalg.norm(r)))

        Me = math.pow(u,2)/math.pow(np.linalg.norm(h),3)*math.pow(1-math.pow(e,2),3/2)*dt*i
        theta = Me + (2*e-math.pow(e,3)/4)*math.sin(Me)+5/4*math.pow(e,2)*math.sin(Me*2)+13/12*math.pow(e,3)*math.sin(Me*3)

        # Use theta in the orbit formula to determine the r distance
        dist = (math.pow(np.linalg.norm(h),2)/u) / (1+e*math.cos(theta))

        # Convert to Cartesian coordinates and save them in orbit arrays
        xorb.append(dist*math.cos(theta))
        yorb.append(dist*math.sin(theta))
        zorb.append(0)

        # Advance loop counter variable
        i = i+1

    print("Orbit generation complete.")
    return [xorb,yorb,zorb]

#--------------------------------------------
# Newtonian Solve
#
# Solve positions using discrete time steps using
# Newtonian gravity/acceleration/velocity vectors
#
# Input: floats dt, G, M, numpy vectors r_0, v_0
# Output: [list(x cors), list(y cors), list(z cors)]
#--------------------------------------------

def newtonian_solve(dt=dt,G=G,M=M,r=r,v=v):
    print("Performing a Newtonian orbit generation...")
    # Initialize empty lists to hold orbital positions
    xorb = []
    yorb = []
    zorb = []
    
    # Save a copy of the initial r vector
    r_0 = r
    
    # Initialize sentinel boolean
    sentinel = True
    
    # Loop to calculate orbital positions
    while sentinel:
        # Calculate changes in velocity and position
        a = r*(-G*M/(math.pow(np.linalg.norm(r),3)))
        dv = a*dt
        v = v + dv
        dr = v*dt
        # Update velocity and position
        r = r + dr
        # Write position vector to holding orbital position vectors
        xorb.append(float(r[0]))
        yorb.append(float(r[1]))
        zorb.append(float(r[2]))  
        
        # Check whether we need to toggle the sentinel
        if np.linalg.norm(r-r_0) < 200e3 and (r-r_0)[1]<0:
            sentinel = False
        
    print("Orbit generation complete.")
    return [xorb,yorb,zorb]

#--------------------------------------------
# Calc Error
#
# Calculate the error
# Newtonian gravity/acceleration/velocity vectors
#
# Input: floats dt, G, M, numpy vectors r_0, v_0
# Output: [list(x cors), list(y cors), list(z cors)]
#--------------------------------------------

def calc_error(nc,oc): #Shorthand for newtonian_cors and orbital_elem_cors
    
    error = []
    if len(nc[0])>len(oc[0]):  # If newtonian longer than orbitalelem
        for i in range(0,len(oc[0])-1):
            error.append(math.sqrt(math.pow(nc[0][i]-oc[0][i],2)+math.pow(nc[1][i]-oc[1][i],2)+math.pow(nc[2][i]-oc[2][i],2)))
            
    else:  # orbitalelem longer than newtonian
        for i in range(0,len(nc[0])-1):
            error.append(math.sqrt(math.pow(nc[0][i]-oc[0][i],2)+math.pow(nc[1][i]-oc[1][i],2)+math.pow(nc[2][i]-oc[2][i],2)))
            
    return error

#--------------------------------------------
# Make Frames
#
# Saves pngs of a series of plots of each orbital
# position which can be used as frames of a gif
#
# Input: [list(x cors), list(y cors), list(z cors)]
# Output: void.  Saves .png files to cwd
#--------------------------------------------
def make_frames(corMatrix):
    print("Generating frames...")
    # cd into current working directory
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)
    
    # Make a directory to store .png files in
    try:
        os.chdir('orbit_pngs')
    except:
        os.mkdir('orbit_pngs')
        os.chdir('orbit_pngs')
    
    # If there are two sets of coordinates to graph, unpack both
    if len(corMatrix) == 2:
        xorb = corMatrix[0][0]
        yorb = corMatrix[0][1]
        zorb = corMatrix[0][2]
        xorb2 = corMatrix[1][0]
        yorb2 = corMatrix[1][1]
        zorb2 = corMatrix[1][2]
    
    else:  # There's only one set of coordinates to graph
        # Unpack the x,y,z coordinates stored in coordinate matrix
        xorb = corMatrix[0]
        yorb = corMatrix[1]
        zorb = corMatrix[2]
    
    #--------------------------------------------
    # Produce video frames from orbit data
    #--------------------------------------------
    n = 50    # Graph every nth frame
    for timeinc in range(0, len(xorb),n):
        
        if timeinc%500<n:
            print(str(int(timeinc/n)) + " of " + str(int(len(xorb)/n)))
        #initialize plot
        fig=plt.figure(frameon=False)        
        ax=fig.gca(projection='3d')
        ax.set_facecolor('black')
        ax.set_aspect(1)        
        ax.view_init(azim=-90, elev=90)
    
        ax.set_xlim3d(-1.1*max(map(abs,xorb)),1.1*max(map(abs,xorb)))
        ax.set_ylim3d(-1.1*max(map(abs,yorb)),1.1*max(map(abs,yorb)))
        ax.set_zlim3d(-1.1*max(map(abs,zorb)),1.1*max(map(abs,zorb)))
        ax.set_axis_off()
        ax.scatter([0],[0],[0],marker='o',color='#B2B2FF',s=100)
            
        ax.plot(xorb,yorb,zorb,ls='solid', lw=1, color='#FF0000') #CCCCCC
        ax.scatter([xorb[timeinc]],[yorb[timeinc]],[zorb[timeinc]],marker='o',color='#FF0000', s=50)
        
        if len(corMatrix)==2:
            ax.plot(xorb2,yorb2,zorb2,ls='solid', lw=1, color='#00FF00')  #BBBBFF
            ax.scatter([xorb2[timeinc]],[yorb2[timeinc]],[zorb2[timeinc]],marker='o',color='#00FF00', s=25)

        timestr= "%04d" % timeinc
        pngname=str(timestr)+".png"
        plt.savefig(str(pngname), axisbg='black')
        plt.clf()
        plt.cla()
        plt.close()        

    print("Finished generating frames.")

#--------------------------------------------
# Make GIF
#
# String all .png files in the orbit_pngs directory
# into a single .gif file
#
# Input: void.  Strings together .png files from cwd
# Output: void.  Saves .gif file to cwd
#--------------------------------------------
def make_gif():
    print("Compiling gif file...")
    # cd into current working directory
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)
    # cd into the png directory
    os.chdir("orbit_pngs")

    # Make a list of filenames ending in .png
    filenames = []
    for file in glob.glob("*.png"):
        filenames.append(file)
    filenames.sort()
        
    # Read in the images from every .png filename
    images = []
    for filename in filenames:
        images.append(imageio.imread(filename))

    # Write the images into a gif
    imageio.mimsave('orbit.gif', images)

    print("Gif complete.")

#===============================================================
err_data = main()
