import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot3D(x,y,z):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_trisurf(x, y, z)
    plt.show()

def generate_sphere(radius, center, num=20):
    """
    Taken from poliastro
    """

    u1 = np.linspace(0, 2 * np.pi, num)
    v1 = u1.copy()
    uu, vv = np.meshgrid(u1, v1)

    x_center, y_center, z_center = center

    xx = x_center + radius * np.cos(uu) * np.sin(vv)
    yy = y_center + radius * np.sin(uu) * np.sin(vv)
    zz = z_center + radius * np.cos(vv)

    return xx, yy, zz
