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
import hoomd.md
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

nameString = "data_{}/".format(
    poly_key) + settings['nameString'].format(**settings)

settings['N'] = user_N

# Create directory to store simulation results
if(not os.path.exists("data_{}".format(poly_key))):
    os.mkdir("data_{}".format(poly_key))

# Print simulation information
print('Working on: {:s}\nUsing: {}'.format(nameString, __file__))

print("== using these settings ==")
for k, v in settings.items():
    print("  {:15s}  =  {:}".format(k, v))


# Initialize execution context
hoomd.context.initialize("")  # "--mode=gpu"
print('[I] Initialize . . . . done.')

# # Create snapshot for simulation
# Lattice as starting configuration
a1, a2, a3 = [[1.3, 0, 0],
              [0, 1.3, 0],
              [0, 0, 1.3]]
n = [N, N, N]

if dimensions == 2:
    a1, a2, a3 = [[1.3, 0, 0],
                  [0, 1.3, 0],
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

'''
for p in gall:
    p.diameter = halo_diam
'''

total_N = len(system.particles)  # total number of clusters


# Adjust density
if dimensions == 2:
    vol = math.pi/4*halo_diam**2 * total_N
elif dimensions == 3:
    vol = math.pi/6*halo_diam**3 * total_N

dens = 0.55
if dimensions == 2:
    boxLen = math.sqrt(vol / dens)
    hoomd.update.box_resize(Lx=boxLen, Ly=boxLen,
                            period=None, scale_particles=True)
elif dimensions == 3:
    boxLen = math.pow(vol / dens, 1/3)
    hoomd.update.box_resize(L=boxLen, period=None, scale_particles=True)


print('[II] Snapshot . . . . done.')

# Neighbor list and Potential selection
nl = hoomd.md.nlist.cell()

# WCA is LJ with shift that causes potential to be zero a r_cut
# lj = hoomd.md.pair.lj(r_cut=2**(1/6)*sigma_u, nlist=nl)
dict_coeff = dict(epsilon=epsilon_u, sigma=sigma_u, offset=0)
table = hoomd.md.pair.table(width=100, nlist=nl)
table.pair_coeff.set('core', 'core', func=WCA_corrected,
                     rmin=0.1, rmax=(2 ** (1/6)*sigma_u), coeff=dict_coeff)

'''
# Shifts interaction potential, so that its value is zero at the r_cut
lj.set_params(mode='shift')
lj.pair_coeff.set('core', 'core', epsilon=epsilon_u, sigma=sigma_u)
'''

# # Equilibration
# Integrator selection
hoomd.md.integrate.mode_standard(dt=time_step)
# creates grooup consisting of central particles in rigid body
if settings['integrator'] == 'langevin':
    langevin = hoomd.md.integrate.langevin(
        group=hoomd.group.all(), kT=settings['kT_equil'], seed=settings['seed'])

    langevin.set_gamma(a='core', gamma=settings['fric_coeff'])

elif settings['integrator'] == 'nve':
    nve = hoomd.md.integrate.nve(group=hoomd.group.all())

    nve.randomize_velocities(kT=settings['kT_equil'], seed=settings['seed'])

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
                                    'pressure',
                                    'pair_table_energy'],
                        # 'pair_lj_energy'],
                        period=outputInterval_log,
                        overwrite=True)

g = hoomd.group.all()

logger = hoomd.analyze.log(filename='{:s}_special.log'.format(nameString),
                           quantities=['pair_table_force'],
                           period=outputInterval_log,
                           overwrite=True)

logger.register_callback(
    'pair_table_force', lambda timestep: math.sqrt(sum(i*i for i in table.get_net_force(g))))

gsd = hoomd.dump.gsd(filename='{:s}.gsd'.format(nameString),
                     period=outputInterval_gsd,
                     group=hoomd.group.all(),
                     dynamic=['momentum'],
                     overwrite=True)

# Increase density

for dens in range(5500, 8210, 10):
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


# Decrease density
for dens in range(8200, 5490, -10):
    dens = dens / 10000
    if dimensions == 2:
        boxLen = math.sqrt(vol / dens)
        hoomd.update.box_resize(Lx=boxLen, Ly=boxLen,
                                period=None, scale_particles=True)
    elif dimensions == 3:
        boxLen = math.pow(vol / dens, 1/3)
        hoomd.update.box_resize(L=boxLen, period=None, scale_particles=True)

    hoomd.run(equil_steps, quiet=True)

print('!!!!!!!!!!!!!!!!!!!!!\nFinal Expansion')
print(vol/system.box.get_volume())


end_time = time.time()

sim_time = end_time - start_time

print('TOTAL SIMULATION TIME [s] = \n{}'.format(sim_time))
