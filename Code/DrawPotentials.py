'''
Draw variation of WCA with distance to the 2d Cluster

'''


from matplotlib import pyplot as plt
from matplotlib import patches as patches
import numpy as np
import math


def WCA(r, rmin, rmax, epsilon, sigma):
    V = 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6) + epsilon
    F = 4 * epsilon / r * (12 * (sigma / r)**12 - 6 * (sigma / r)**6)
    return (V, F)


def WCA_corrected(r, rmin, rmax, epsilon, sigma, offset):
    print(offset)
    print(sigma)
    r = r + offset
    V = 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6) + epsilon
    F = 4 * epsilon / r * (12 * (sigma / r)**12 - 6 * (sigma / r)**6)
    return (V, F)


# Adjust number of halo spheres
N_spheres = 3

# #
# #
# #
# #
# Cluster data non adjusted
d = 10
theta_2d = 2*math.pi/N_spheres
radius_2d = d/math.sin(theta_2d)*math.sin(math.pi/2 - theta_2d/2)

# Plotting data

x1 = np.linspace(1, 2**(1/6)*10, 100)
y1, Force = WCA(x1, 0, 0, 1, 10)

x2 = np.linspace(1, 2**(1/6)*d, 100)

y2, Force = WCA(x2, 0, 0, 1, d)

x2 = x2 + radius_2d

# Plots
fig, ax = plt.subplots()
plt.grid(linestyle='--')
ax.set_aspect(1)

# Draw cluster non adjusted
for i in range(N_spheres):
    coord = (radius_2d*math.cos(theta_2d*i), radius_2d*math.sin(theta_2d*i))
    core_diam = (radius_2d-d/2)*2
    circle2 = plt.Circle(coord, color='g', radius=d/2, fill=True, alpha=0.5)
    ax.add_artist(circle2)

circle2 = plt.Circle((0, 0), color='g', radius=radius_2d -
                     d/2, fill=True, alpha=0.5)
ax.add_artist(circle2)

circle1 = plt.Circle((0, 0), 5, color='b', fill=False, alpha=0.5)
ax.add_artist(circle1)
plt.plot(x1, y1, 'b', label='disk')
plt.plot(x2, y2, 'g', label='ring')
plt.legend()

plt.hlines(0, -(radius_2d + d/2), 25, colors='k',
           linestyle='--', linewidth=0.75)
plt.vlines(0, -(radius_2d + d/2), 25, colors='k',
           linestyle='--', linewidth=0.75)
plt.vlines(radius_2d, -(radius_2d + d/2), 25, colors='k',
           linestyle='--', linewidth=0.75)
plt.ylim(-(radius_2d + d/2), 25)
plt.xlim(-(radius_2d + d/2), 25)

# #
# #
# #
# #
# Cluster data adjusted for cluster size
theta_2d = 2*math.pi/N_spheres
const = 1/math.sin(theta_2d)*math.sin(math.pi/2 - theta_2d/2)
radius_2d = 10*const/(1+2*const)
d = 10 - 2*radius_2d


# Plotting data
x1 = np.linspace(0.1, 2**(1/6)*10, 100)
y1, Force = WCA(x1, 0, 0, 1, 10)

x2 = np.linspace(0.1, 2**(1/6)*d, 100)
y2, Force = WCA(x2, 0, 0, 1, d)

x3 = np.linspace(0.1, 2**(1/6)*10 - radius_2d, 100)
y3, Force = WCA_corrected(x3, 0, 0, 1, 2*radius_2d+d, radius_2d)

x2 = x2 + radius_2d
x3 = x3 + radius_2d


# Potential data adjusted for cluster size

# Plots
fig, ax = plt.subplots()
plt.grid(linestyle='--', linewidth=0.75, color='lightgray')
ax.set_aspect(1)
plt.rcParams.update({'mathtext.default': 'regular'})

# Draw cluster adjusted
for i in range(N_spheres):
    coord = (radius_2d*math.cos(theta_2d*i), radius_2d*math.sin(theta_2d*i))
    core_diam = (radius_2d-d/2)*2
    circle2 = plt.Circle(coord, color='g', radius=d/2, fill=True, alpha=0.5)
    ax.add_artist(circle2)

circle2 = plt.Circle((0, 0), color='g', radius=radius_2d -
                     d/2, fill=True, alpha=0.5)
ax.add_artist(circle2)

circle1 = plt.Circle((0, 0), 5, color='b', fill=False, alpha=0.5)
ax.add_artist(circle1)

plt.plot(x1, y1, 'b', label='disk')
plt.plot(x2, y2, 'g', label='ring')
plt.plot(x3, y3, 'r', linestyle='--', linewidth=0.75, label='corrected')
plt.legend()

plt.hlines(0, -5, 15, colors='k',
           linestyle='--', linewidth=0.75)
plt.vlines(0, -5, 15, colors='k',
           linestyle='--', linewidth=0.75)
plt.vlines(radius_2d, -5, 15, colors='k',
           linestyle='--', linewidth=0.75)
plt.ylim(-5, 15)
plt.xlim(-5, 15)
plt.ylabel('$U_{WCA}$ / -')
plt.xlabel('r / -')

plt.show()
plt.clf()
