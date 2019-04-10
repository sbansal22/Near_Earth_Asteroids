from astropy import units as u
import numpy as np



def calculateNewtonianVelocity(mass1, radius):
    """
    Calculating newtonian velocity
    M - mass of heavier object that your object is rotating around
    r - distance between objects
    """
    v = np.sqrt((G*mass1)/radius)
    return v



def periapsisDistanceapoapsisDistance(ecc, a):
    """
    Takes eccentricity and semimajor axis and returns periapsis distance and apoasis distance
    """
    Rp = a*(1-ecc)
    Ra = a*(1+ecc)
    return Rp, Ra

def _kepler_equation(E, M, ecc):
    """
    Taken from poliastro : https://github.com/poliastro/poliastro/blob/master/src/poliastro/core/angles.py
    """
    return E - ecc * np.sin(E) - M

def _kepler_equation_prime(E, M, ecc):
    return 1 - ecc * np.cos(E)

def newtonsmethodElliptical(x0, M, ecc, tol=1.48e-08, maxiter=50):
    """
    Calling Newton's Method to solve Kepler's equation
    Taken from poliastro : https://github.com/poliastro/poliastro/blob/master/src/poliastro/core/angles.py

    """
    p0 = 1.0 * x0
    for iter in range(maxiter):
        fval = _kepler_equation(p0, M, ecc)
        fder = _kepler_equation_prime(p0, M, ecc)

        newton_step = fval / fder
        p = p0 - newton_step
        if abs(p - p0) < tol:
            return p
        p0 = p
    return p0
def eccentricityAnomaly(M, ecc):
    """
    Eccentricity anomaly taken from mean anomaly
    """
    E = newtonsmethodElliptical( M, M, ecc)
    return E

def calculateEllipticalTrueAnomaly(M, ecc):
    """
    Taken from poliastro : https://github.com/poliastro/poliastro/blob/master/src/poliastro/core/angles.py
    """
    E = eccentricityAnomaly(M, ecc)
    beta = ecc / (1 + np.sqrt((1 - ecc) * (1 + ecc)))
    nu = E + 2 * np.arctan(beta * np.sin(E) / (1 - beta * np.cos(E)))
    return nu


def calculatePosition(nu, ecc, a):
    """
    Calculates the Distance of a secondary from a primary in an elliptical orbit
    """
    r = (a*(1-ecc**2))/ (1+ecc*np.cos(nu))
    return r

def calculatePeriapsisAngle(nu, ecc, a):
    """
    Calculates the angle in an orbit in an elliptical orbit
    """
    phi = np.arctan((ecc*np.sin(nu))/ (1+ecc*np.cos(nu)))
    return phi
def convertToCartesian(periangle, i, long):
    """
    Takes the orbital elements which are euler angles and returns cartesian coordinates
    """
    ex = np.array([(np.cos(periangle)*np.cos(long) - np.sin(long)*np.sin(periangle)*np.cos(i)) , (np.sin(long)*np.sin(periangle) + np.cos(long)*np.sin(periangle)*np.cos(i)) , np.sin(long)*np.sin(i)])
    ey = np.array([(-1*np.cos(long)*np.sin(long) - np.sin(long)*np.cos(periangle)*np.cos(i)),(-1*np.sin(long)*np.sin(periangle) + np.cos(periangle)*np.cos(long)*np.cos(i)), (np.cos(periangle)*np.sin(i)) ])
    ez = np.array([(np.sin(long)*np.sin(i)), (-np.cos(long)*np.sin(i)), (np.cos(i))])
    return ex,ey,ez

def convertToCartesian2(periangle, i, long, r):
    """
    Try 2
    """
    X = r*(np.cos(long)*np.cos(periangle) - np.sin(long)*np.sin(periangle)*np.cos(i))
    Y = r*(np.sin(long)*np.cos(periangle) + np.cos(long)*np.sin(periangle)*np.cos(i))
    Z = r*(np.sin(periangle)*np.sin(i))
    return X,Y,Z

def createRange(M):
    """
    creates list of all the mean anomalies to plot with
    """
    init_M = M
    range1 = list(range(int(init_M), 360,5))
    range2 = list(range(0, int(init_M),5))
    rangeT = range1 + range2
    return rangeT
def createMotion(ecc, a, M, long):
    """
    Iterates the functions to calculate cartesian coordinates
    for 5 degree increments of mean anomaly until one iteration
    has passed

    """
    M_range = createRange(M)

    nu_range = [calculateEllipticalTrueAnomaly(M, ecc) for M in M_range]

    pos_range = [calculatePosition(nu, ecc, a) for nu in nu_range]
    angle_range = [calculatePeriapsisAngle(nu, ecc, a) for nu in nu_range]
    cartesianx = []
    cartesiany = []
    cartesianz = []
    i = 1
    for angle, pos in zip(angle_range, pos_range):
        x,y,z = convertToCartesian2(angle, i, long, pos)
        cartesianx.append(x)
        cartesiany.append(y)
        cartesianz.append(z)
        i += 1
    return cartesianx, cartesiany, cartesianz
