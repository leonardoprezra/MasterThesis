'''
HARD SPHERES SIMULATION

Simulation of Platonic Polyhedral Sphere Clusters (PSCs)
based on the work of Marson et al
https://doi.org/10.1039/C9SM00664H

Clusters consist of a 'core' type particle in its center of mass,
and it is surrounded by 'halo' type particles

Regarding HOOMD-blue nomenclature of rigid bodies
Core particles = central
Halo particles = constituent

PSCs keys:
octa
tetra
ico
cube
dode
2Dspheres
3Dspheres
one

Integrators:
local
langevin
'''

import hoomd
import hoomd.hpmc
import random
import math
from MarsonFunctions import core_properties, mom_inertia, settings, PartCluster, create_snapshot, unif_pos, WCA_corrected

import os
import sys
import time

'''
# Attach VS Code debugger for debugging invoked files
import ptvsd
# 5678 is the default attach port in the VS Code debug configurations
print("Waiting for debugger attach")
ptvsd.enable_attach(address=('localhost', 3002), redirect_output=True)
ptvsd.wait_for_attach()
'''

start_time = time.time()

# Update settings from arguments
for i in range(1, len(sys.argv)):
    try:
        para, value = sys.argv[i].split("=")
        if para in settings:
            settings[para] = settings[para].__class__(value)
        else:
            print("Found an invalid Argument: "+sys.argv[i])
            sys.exit(5)
    except:
        print("Found an invalid Argument: "+sys.argv[i])
        exit()

# Updates parameters
N = settings['N']
halo_diam = settings['diameter']
poly_key = settings['poly']
halo_mass = settings['mass']
density = settings['density']
dimensions = settings['dimensions']

outputInterval_gsd = settings['outputInterval_gsd']
outputInterval_log = settings['outputInterval_log']
time_step = settings['time_step']


equil_steps = settings['equil_steps']
N_cluster = 0
settings['N_cluster'] = 0

sigma_u = settings['diameter']*settings['ratio']
epsilon_u = settings['epsilon']

# Print the correct number of clusters in the system
user_N = settings['N']
if settings['dimensions'] == 3:
    settings['N'] = settings['N']**3
elif settings['dimensions'] == 2:
    settings['N'] = settings['N']**2

nameString = "dataHighP_{}/".format(
    poly_key) + settings['nameString'].format(**settings)

settings['N'] = user_N

# Create directory to store simulation results
if(not os.path.exists("dataHighP_{}".format(poly_key))):
    os.mkdir("dataHighP_{}".format(poly_key))

# Print simulation information
print('Working on: {:s}\nUsing: {}'.format(nameString, __file__))

print("== using these settings ==")
for k, v in settings.items():
    print("  {:15s}  =  {:}".format(k, v))


print('halo_diam={}'.format(halo_diam))

# Initialize execution context
hoomd.context.initialize("")  # "--mode=gpu"
print('[I] Initialize . . . . done.')

# # Create snapshot for simulation
# Lattice as starting configuration
a1, a2, a3 = [[1.3*halo_diam, 0, 0],
              [0, 1.3*halo_diam, 0],
              [0, 0, 1.3*halo_diam]]
n = [N, N, N]

if dimensions == 2:
    a1, a2, a3 = [[1.3*halo_diam, 0, 0],
                  [0, 1.3*halo_diam, 0],
                  [0, 0, 1]]
    n = [N, N]

# Creates lattice

uc = hoomd.lattice.unitcell(N=1,
                            a1=a1,
                            a2=a2,
                            a3=a3,
                            dimensions=dimensions,
                            position=[(0, 0, 0)],
                            type_name=['core'],
                            mass=[halo_mass],
                            diameter=[halo_diam])


# Initialize sys configuration
system = hoomd.init.create_lattice(unitcell=uc, n=n)

gall = hoomd.group.all()

total_N = len(system.particles)  # total number of clusters


# Adjust density
if dimensions == 2:
    vol = math.pi/4*halo_diam**2 * total_N
elif dimensions == 3:
    vol = math.pi/6*halo_diam**3 * total_N


'''
dens = 0.3
if dimensions == 2:
    boxLen = math.sqrt(vol / dens)
    hoomd.update.box_resize(Lx=boxLen, Ly=boxLen,
                            period=None, scale_particles=True)
elif dimensions == 3:
    boxLen = math.pow(vol / dens, 1/3)
    hoomd.update.box_resize(L=boxLen, period=None, scale_particles=True)
'''

print('[II] Snapshot . . . . done.')

# # Equilibration
# Integrator selection
mc = hoomd.hpmc.integrate.sphere(d=halo_diam*0.1, seed=1)

mc.shape_param.set('core', diameter=halo_diam)

# Adjust density before actual run
print('!!!!!!!!!!!!!!!!!!!!!\nPre-Initial Compression')
print(vol/system.box.get_volume())

pre_dens = vol/system.box.get_volume()
dens = 0.49

if dimensions == 2:
    boxLen = math.sqrt(vol / dens)
    hoomd.update.box_resize(Lx=boxLen, Ly=boxLen,
                            period=None, scale_particles=True)

elif dimensions == 3:
    for dens in range(int(pre_dens*10000), int(dens*10000), 100):
        dens = dens / 10000
        boxLen = math.pow(vol / dens, 1/3)
        hoomd.update.box_resize(
            L=boxLen, period=None, scale_particles=True)

        hoomd.run(500, quiet=True)

# Store snapshot information

log = hoomd.analyze.log(filename='{:s}.log'.format(nameString),
                        quantities=['volume',
                                    'momentum',
                                    'time',
                                    'potential_energy',
                                    'kinetic_energy',
                                    'translational_kinetic_energy',
                                    'rotational_kinetic_energy',
                                    'temperature',
                                    'pressure'],
                        # 'pair_lj_energy'],
                        period=outputInterval_log,
                        overwrite=True)

hoomd.hpmc.analyze.sdf(mc=mc, filename='{:s}_sdf.log'.format(
    nameString), xmax=0.002, dx=1e-5, navg=100, period=outputInterval_log)

gsd = hoomd.dump.gsd(filename='{:s}.gsd'.format(nameString),
                     period=outputInterval_gsd,
                     group=hoomd.group.all(),
                     dynamic=['momentum'],
                     overwrite=True)

# Increase density

for dens in range(4900, 7410, 10):
    dens = dens / 10000
    if dimensions == 2:
        boxLen = math.sqrt(vol / dens)
        hoomd.update.box_resize(Lx=boxLen, Ly=boxLen,
                                period=None, scale_particles=True)
    elif dimensions == 3:
        boxLen = math.pow(vol / dens, 1/3)
        hoomd.update.box_resize(L=boxLen, period=None, scale_particles=True)

    hoomd.run(equil_steps, quiet=True)

print('!!!!!!!!!!!!!!!!!!!!!\nFinal Compression')
print(vol/system.box.get_volume())


end_time = time.time()

sim_time = end_time - start_time

print('TOTAL SIMULATION TIME [s] = \n{}'.format(sim_time))
