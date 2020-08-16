import numpy as np
from .utility import CONSTANT, TypeOE
import utility


def calc_oe_from_sv(R: np.array, V: np.array) -> TypeOE:
    """ Calculate the orbital elements from the state vectors.
        Input:
            - R: satellite's position vector        [km]
            - V: satellite's velocity vector        [km/s]
        Output:
            - e: eccentricity                       [ ]
            - a: semimajor axis                     [km]
            - i: inclination                        [rad]
            - ra: longitude of the ascending node   [rad]
            - w: argument of periapsis              [rad]
            - MA: mean anomaly                      [rad]
    """
    OE: TypeOE = {}

    # Set a number below which the eccentricity is considered to be zero
    eps = 0.0000000001

    # Calculate the distance
    r = np.linalg.norm(R)

    # Calculate the speed
    v = np.linalg.norm(V)

    # Calculate the radial velocity
    vr = np.dot(R, V) / r

    # Calculate the specific angular momentum and its magnitude
    H = np.cross(R, V)
    h = np.linalg.norm(H)

    # Calculate the inclination
    OE['inclination'] = np.arccos(H[2] / h)

    # Calculate vector N that defines the node line and its magnitude
    N = np.cross(np.array([0, 0, 1]), H)
    n = np.linalg.norm(N)

    # Calculate the right ascension (longitude) of the ascending node
    if n != 0:
        ra = np.arccos(N[0] / n)
        if N[1] < 0:
            ra = 2 * np.pi - ra
    else:
        ra = 0
    OE['right_ascension'] = ra

    # Calculate the eccentricity
    E = 1 / CONSTANT.GM * ((v * v - CONSTANT.GM / r) * R - r * vr * V)
    e = np.linalg.norm(E)
    OE['eccentricity'] = e

    # Calculate the argument of perigee
    if n != 0:
        if e > eps:
            w = np.arccos(np.dot(N, E) / n / e)
            if E[2] < 0:
                w = 2 * np.pi - w
        else:
            w = 0
    else:
        w = 0
    OE['argument_of_perigee'] = w

    # Calculate the true anomaly
    if e > eps:
        TA = np.arccos(np.dot(E, R) / e / r)
        if vr < 0:
            TA = 2 * np.pi - TA
    else:
        cp = np.cross(N, R)
        TA = np.arccos(np.dot(N, R) / n / r)
        if cp[2] < 0:
            TA = 2 * np.pi - TA

    # Calculate the eccentric anomaly
    EA = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(TA / 2))

    # Calculate the mean anomaly
    MA = EA - e * np.sin(EA)
    OE['mean_anomaly'] = MA

    # Calculate the semimajor axis
    a = h * h / CONSTANT.GM / (1 - e * e)
    OE['semimajor_axis'] = a

    return OE


if __name__ == "__main__":
    print("Hello the world")
