import math
from MarsonFunctions import core_properties, mom_inertia, settings
import numpy as np

N = 20
halo_diam = 0.5
halo_mass = 1
density = 0.60
dimensions = 2

outputInterval = 1
time_step = 0.005/2

therm_steps = 10000
npt_steps = 10000
equil_steps = 10000

kT = 1.0

poly_key = '2Dspheres'
N_cluster = 10


def get_ratios(poly=poly_key, d=halo_diam, N_cluster=10):
    sphere_diam = core_properties(poly='{}_sphere_diam'.format(
        poly_key), d=d, N_spheres=N_cluster)  # Diameter of sphere circunscribing PSC + halo_diam

    theta = 2*math.pi/N_cluster
    A = d**2 / (4*math.sin(theta)) * math.sin((math.pi-theta)/2) * \
        math.cos(theta/2) + d**2 / 4 * math.cos(math.pi/6) + d**2 / 8 * \
        (math.pi/6 + theta/2) - d**2 / 48 * math.pi
    area_psc = 2*N_cluster*A

    area_sphere = math.pi/4 * sphere_diam**2

    return [area_psc, area_sphere, area_psc/area_sphere]


ratios = []

for N_cluster in range(1, 20):
    ratios.append(get_ratios(N_cluster=N_cluster))

print(np.average(np.array(ratios)[2:-1, 2]))
#print([i[2] for i in ratios])
