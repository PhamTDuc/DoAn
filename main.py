
# Step 01
import numpy as np
from package.utility import CONSTANT, calMeanAnomaly, calEccentricAnomaly, calVecInPQW, getMatECItoPQW
from package.orbital_elements import calc_oe_from_sv
# Three Continious TimePoints
time_points = np.array([3900, 4000, 4100], dtype='float64')


# Locations of Observers
R1 = np.array([-2936.2922, -5654.01, 0], dtype='float64')
R2 = np.array([-389.5577, -6359.079, 0], dtype='float64')
R3 = np.array([2224.1957, -5970.1419, 0], dtype='float64')


# Directional Vectors from Observers to Satellite
p1 = np.array([-0.39868932, 0.91626844, -0.0387166], dtype='float64')
p2 = np.array([-0.56171108, 0.82642542, 0.03875162], dtype='float64')
p3 = np.array([-0.70363441, 0.70400468, 0.09631214], dtype='float64')

# Preparation for Calculating Satellite Orbit
t1 = time_points[0] - time_points[1]
t3 = time_points[2] - time_points[1]
t = t3 - t1


D0 = p1 @ np.cross(p2, p3)
D11 = R1 @ np.cross(p2, p3)
D21 = R2 @ np.cross(p2, p3)
D31 = R3 @ np.cross(p2, p3)

D12 = R1 @ np.cross(p1, p3)
D22 = R2 @ np.cross(p1, p3)
D32 = R3 @ np.cross(p1, p3)

D13 = R1 @ np.cross(p1, p2)
D23 = R2 @ np.cross(p1, p2)
D33 = R3 @ np.cross(p1, p2)


# Calculate Prameters for "r2" Equation
A = (-D12 * t3 / t + D22 + D32 * t1 / t) / D0
B = (-D12 * (t**2 - t3**2) * t3 / t + D32 * (t**2 - t1**2) * t1 / t) / 6 / D0
E = R2 @ p2
Muy = CONSTANT.G * CONSTANT.M


a = -(A**2 + 2 * A * E + np.linalg.norm(R2)**2)
b = -2 * Muy * B * (A + E)
c = -Muy**2 * B**2


# Solve |r2|
# Using Newton's method to Solve Equation of |r2|
# (x^8 + a*x^6 + b*x^3+  c = 0)
fx = np.poly1d([1, 0, a, 0, 0, b, 0, 0, c])
dx = fx.deriv()

r2 = 6000
alias = 0.1
counter = 0

while (counter < 10000 and np.abs(fx(r2)) > alias):
    # print(counter)
    r2 -= fx(r2) / dx(r2)
    counter += 1


# Find Other Constants (C1, C3)
C1 = t3 * (1 + Muy / 6 / r2**3 * (t**2 - t3**2)) / t
C3 = -t1 * (1 + Muy / 6 / r2**3 * (t**2 - t1**2)) / t


# Find Distances from Observations to Satellite (mul_p1, mul_p2, mul_p3)
mul_p1 = (-D11 + D21 / C1 - C3 * D31 / C1) / D0
mul_p2 = A + Muy * B / r2**3
mul_p3 = (-C3 * D13 / C1 + D23 / C3 - D33) / D0


# Find Position Vectors of Satellite (r1, r2, r3)
vec_r1 = R1 + mul_p1 * p1
vec_r2 = R2 + mul_p2 * p2
vec_r3 = R3 + mul_p3 * p3


# Find Velocity of Satellite at t2 (v2)
f1 = 1 - Muy * t1**2 / 2 / r2
g1 = t1 - Muy * t1**3 / 6 / r2**3
f3 = 1 - Muy * t3**2 / 2 / r2
g3 = t3 - Muy * t3**3 / 6 / r2**3

v2 = (f1 * vec_r3 - f3 * vec_r1) / (f1 * g3 - f3 * g1)  # Velocity v2
print("v2 = ", v2)


# Calculate OE from (R, V)
OE = calc_oe_from_sv(vec_r2, v2)
print("OE=", OE)

# Calculate MeanAbnomaly and EccentricAbnomaly
mean_anomaly = calMeanAnomaly(OE['mean_anomaly'], OE['semimajor_axis'], t1)
print("Mean Anomaly", mean_anomaly)
eccentric_anomaly = calEccentricAnomaly(mean_anomaly, OE['eccentricity'], 1)
print("EccentricAnomaly", eccentric_anomaly)


# Find Position in PQW
r1_PQW = calVecInPQW(OE['semimajor_axis'], OE['eccentricity'], eccentric_anomaly)
toECI = getMatECItoPQW(i=OE['inclination'], omega=OE['argument_of_perigee'], sigma=OE['right_ascension'])
print("r1 in PQW", r1_PQW)
print(vec_r1)
print("r1 in ECI", np.linalg.inv(toECI) @ r1_PQW)
