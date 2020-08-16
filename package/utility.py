from typing import TypedDict, NamedTuple
import numpy as np
from datetime import datetime

print("Module Utility in Package")
class CONSTANT():
    __slots__ = ()

    def __new__(cls, *args, **kwargs):
        raise RuntimeError("CONSTANT class can't have instances")

    PI = 3.14159265358979323846
    G = 6.67408 * 10**-20  # Gravitional Constant (km^3 kg^-1 s^-2)
    M = 5.972 * 10**24  # Earth Mass (kg)
    R = 6367.5  # Earth Radius(km)
    GM = G * M  # Earth's standard gravitational parameter [km^3.s^-2]

class TypeOE(TypedDict, total=False):
    """ ISS (ZARYA)
    1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927
    2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537

    1   01–01   Line number 2
    2   03–07   Satellite Catalog number                        25544
    3   09–16   Inclination (degrees)                           51.6416
    4   18–25   Right Ascension of the Ascending Node (degrees) 247.4627
    5   27–33   Eccentricity (decimal point assumed)            0006703
    6   35–42   Argument of Perigee (degrees)                   130.5360
    7   44–51   Mean Anomaly (degrees)                          325.0288
    8   53–63   Mean Motion (revolutions per day)               15.72125391
    9   64–68   Revolution number at epoch (revolutions)        56353
    10  69–69   Checksum (modulo 10)                            7


    """
    eccentricity: float
    semimajor_axis: float
    inclination: float
    right_ascension: float
    argument_of_perigee: float
    mean_anomaly: float

    static = CONSTANT.GM**(1/3)

    @classmethod
    def calSemimajorAxis(cls, mean_motion:float)->float:
        """The TLE gives mean motion (n) in rev/day. This needs to be converted to rad/s which can be accomplished by multiplying the n TLE value by 2π/86400.
        Therefore, to go directly from n in TLE to the semi-major axis a. We can use the following formula:
                GM^(1/3)
        a =___________________________

            (2π*mean_motion *  /86400)^(2/3)
        For n=15.5918272revday, we get a=6768.16km.
        """

        a = cls.static
        a /= (2*CONSTANT.PI*mean_motion/86400)**(2/3)
        return a

class TypeLatLong(NamedTuple):
    """The "latitude" (abbreviation: Lat., φ, or phi) of a point on Earth's surface is the angle between the equatorial plane and the straight line that passes through that point and through (or close to) the center of the Earth.
    Lines joining points of the same latitude trace circles on the surface of Earth called parallels, as they are parallel to the Equator and to each other. 
    The North Pole is 90° N; the South Pole is 90° S. The 0° parallel of latitude is designated the Equator, the fundamental plane of all geographic coordinate systems. 
    The Equator divides the globe into Northern and Southern Hemispheres.

    The "longitude" (abbreviation: Long., λ, or lambda) of a point on Earth's surface is the angle east or west of a reference meridian to another meridian that passes through that point. 
    All meridians are halves of great ellipses (often called great circles), which converge at the North and South Poles. 
    The meridian of the British Royal Observatory in Greenwich, in southeast London, England, is the international prime meridian, 
    although some organizations—such as the French Institut Géographique National—continue to use other meridians for internal purposes. 
    The prime meridian determines the proper Eastern and Western Hemispheres, although maps often divide these hemispheres further west in order to keep the Old World on a single side. 
    The antipodal meridian of Greenwich is both 180°W and 180°E. This is not to be conflated with the International Date Line, 
    which diverges from it in several places for political and convenience reasons, including between far eastern Russia and the far western Aleutian Islands.
    """
    lattitude: float
    longtitude: float

    def __repr__(self):
        return f"LatLong(Lattitude={self.lattitude}, Longtitude={self.longtitude})"


class TypeXYZ(NamedTuple):
    """Every point that is expressed in ellipsoidal coordinates can be expressed as an rectilinear x y z (Cartesian) coordinate.
    Cartesian coordinates simplify many mathematical calculations. The Cartesian systems of different datums are not equivalent.
    """
    x: float
    y: float
    z: float

    def __repr__(self):
        return f"Coordinate(X={self.x:0.4f}, Y={self.y:0.4f}, Z={self.z:g})"


# CosinMatrix Transfrom ECI to PQW:
def getMatECItoPQW(omega: float=0, i: float=0, sigma: float=0) -> np.array:

    RzO = np.array([[np.cos(omega), np.sin(omega), 0], [-np.sin(omega), np.cos(omega), 0], [0, 0, 1]])
    RzS = np.array([[np.cos(sigma), np.sin(sigma), 0], [-np.sin(sigma), np.cos(sigma), 0], [0, 0, 1]])
    RxI = np.array([[1, 0, 0], [0, np.cos(i), np.sin(i)], [0, -np.sin(i), np.cos(i)]])

    return RzO @ RxI @ RzS


# Calculate Mean Abnomaly
def calMeanAnomaly(M0: float, a: float, t: float) -> float:
    """
    Calculate MeanAbnomaly of Satellite which Semimajor Axis a(km) after t(s) period
    Kepler III
    a^3      GM
   ----- = -------
    T^2     4 pi^2

    M = M0 + 2 pi*t/T
    """
    T = np.sqrt(4 * a**3 * CONSTANT.PI**2 / CONSTANT.G / CONSTANT.M)
    return M0 + 2 * CONSTANT.PI / T * t

# Calculate EccentricAbnomaly from MeanAbnomaly Using Newton Method


def calEccentricAnomaly(M, e, E0=1, repeated=10000, alias=0.1) -> float:
    """
    Using Kepler Equation to solve E(rad)
    M = E - e * sin(E)
    """
    counter = 0
    while(counter < repeated and np.abs(M + e * np.sin(E0) - E0) > alias):
        E0 -= (E0 - e * np.sin(E0) - M) / (1 - e * np.cos(E0))
        counter += 1

    return E0


def calVecInPQW(a: float, e: float, E: float) -> np.array:
    """
    p = a*cosE -c
    q = b*sinE
    w = 0
    c=a*e
    b=a*sqrt(1-e**2)
    """
    b = a * np.sqrt(1 - e**2)
    c = a * e
    return np.array([a * np.cos(E) - c, b * np.sin(E), 0], dtype='float64')


# Calculate GTMS in a arbitary Timepoint: datetime.datetime()
def getGMST(time: datetime = datetime.now()) -> float:
    # Calculate JND
    a = (14 - time.month) // 12
    y = time.year + 4800 - a
    m = time.month + 12 * a - 3
    JND = time.day + (153 * m + 2) // 5 + 365 * y + y // 4 - y // 100 + y // 400 - 32045

    # Calculate Century from J2000
    T0 = (JND - 2451545) / 36525

    # Calculate GTMS[deg] at 0:00 UT
    G0 = 100.4606184 + 36000.77004 * T0 + 0.000387933 * T0**2 - 0.00000002583 * T0**3

    t = (time.hour - 12) / 24 + time.minute / 1440 + time.second / 86400
    GMST = G0 + 360.98564724 * t
    return GMST


# Calculate Rotation Frame from ECEF to ECI
def getRie(time: datetime=datetime.now()) -> np.array:
    gmst = getGMST(time)
    rad = CONSTANT.PI * gmst / 180
    Rie = np.array([[np.cos(rad), np.sin(rad), 0], [-np.sin(rad), np.cos(rad), 0], [0, 0, 1]])
    return Rie


# Coordinates in Lad and Long in ECI
def toECIfromLatLong(location:TypeLatLong  , time: datetime=datetime.now()) -> np.array:
    x = CONSTANT.R * np.cos(location.lattitude) * np.cos(location.longtitude)
    y = CONSTANT.R * np.cos(location.lattitude) * np.sin(location.longtitude)
    z = CONSTANT.R * np.sin(location.lattitude)
    coord = np.array([x, y, z])
    coord = getRie(time).T @ coord
    return coord


def normalize(vec: np.array) -> np.array:
    """Normalize 1D Vector Array 
    """
    norm = np.linalg.norm(vec)
    if norm < 0.001:
        return vec
    return vec / norm


if __name__ == "__main__":
    print(CONSTANT.PI)
    print(CONSTANT.G)
    print(CONSTANT.GM)
    # print(toECIfromLatLong())
