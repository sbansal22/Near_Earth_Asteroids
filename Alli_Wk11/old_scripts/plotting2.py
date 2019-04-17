from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from itertools import product, combinations
import PIL
from mpl_toolkits.mplot3d import proj3d
import constants
from constants import Earth



def plotPoint(ex,ey,ez):
    """

    """

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect("equal")



    # draw sphere
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = np.cos(u)*np.sin(v)
    y = np.sin(u)*np.sin(v)
    z = np.cos(v)
    r = 6378.137


    ax.plot_wireframe(r*x, r*y, r*z, color="r")
    ex = -1.6228595368879668
    ey = 0.7337559841691161
    ez =  0.037899361259070535

    # draw a point
    ax.scatter([ex], [ey], [ez], color="g", s=100)


    %config InlineBackend.figure_format = "png"
    plt.show()


def otherplot():


    # load bluemarble with PIL
    bm = PIL.Image.open('bluemarble.jpg')
    # it's big, so I'll rescale it, convert to array, and divide by 256 to get RGB values that matplotlib accept
    bm = np.array(bm.resize([int(d/5) for d in bm.size]))/256.

    # coordinates of the image - don't know if this is entirely accurate, but probably close
    lons = np.linspace(-180, 180, bm.shape[1]) * np.pi/180
    lats = np.linspace(-90, 90, bm.shape[0])[::-1] * np.pi/180

    # repeat code from one of the examples linked to in the question, except for specifying facecolors:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    x = np.outer(np.cos(lons), np.cos(lats)).T
    y = np.outer(np.sin(lons), np.cos(lats)).T
    z = np.outer(np.ones(np.size(lons)), np.sin(lats)).T
    ax.plot_surface(x, y, z, rstride=4, cstride=4, facecolors = bm)

    %config InlineBackend.figure_format = "png"
    plt.show()



%config InlineBackend.figure_format = "png"
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect("equal")



# draw sphere
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = np.cos(u)*np.sin(v)
y = np.sin(u)*np.sin(v)
z = np.cos(v)
r = 6378.137


ax.plot_wireframe(r*x, r*y, r*z, color="r")
ex = -1.6228595368879668
ey = 0.7337559841691161
ez =  0.037899361259070535

# draw a point
ax.scatter([ex], [ey], [ez], color="g", s=100)

plt.show()
