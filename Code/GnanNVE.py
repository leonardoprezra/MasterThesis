'''
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
import io
from contextlib import redirect_stdout

import hoomd
import hoomd.md
import random
import math
from GnanFunctions import create_snapshot_soft, hertzian
from MarsonFunctions import core_properties, mom_inertia, settings, PartCluster, create_snapshot, unif_pos

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
N_cluster = settings['N_cluster']

sigma_u = settings['diameter']*settings['ratio']
epsilon_u = settings['epsilon']

fene_k = settings['fene_k']
fene_r0 = settings['fene_r0']
harm_k = settings['harm_k']

# Print the correct number of clusters in the system
user_N = settings['N']
if settings['dimensions'] == 3:
    settings['N'] = settings['N']**3
elif settings['dimensions'] == 2:
    settings['N'] = settings['N']**2

nameString = "dataFLEX_{}/".format(
    poly_key) + settings['nameString'].format(**settings)

settings['N'] = user_N

# Create directory to store simulation results
if(not os.path.exists("dataFLEX_{}".format(poly_key))):
    os.mkdir("dataFLEX_{}".format(poly_key))

# Print simulation information
print('Working on: {:s}\nUsing: {}'.format(nameString, __file__))

print("== using these settings ==")
for k, v in settings.items():
    print("  {:15s}  =  {:}".format(k, v))

# Core particle properties
cluster = PartCluster(
    poly_key=poly_key, N_cluster=N_cluster, halo_diam=halo_diam, halo_mass=halo_mass, ratio=settings['ratio'])
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print('halo_diam={}'.format(cluster.halo_diam))
print('core_diam={}'.format(cluster.core_diam))
print('sphere_diam={}'.format(cluster.sphere_diam))

# Initialize execution context
hoomd.context.initialize("")  # "--mode=gpu"
print('[I] Initialize . . . . done.')

# # Create snapshot for simulation
# Lattice as starting configuration
system, total_N = create_snapshot_soft(
    cluster=cluster, dimensions=dimensions, N=N)


# Neighbor list and Potential selection
nl = hoomd.md.nlist.cell()

# WCA is LJ with shift that causes potential to be zero a r_cut
lj = hoomd.md.pair.lj(r_cut=2**(1/6)*cluster.halo_diam, nlist=nl)

# Shifts interaction potential, so that its value is zero at the r_cut
lj.set_params(mode='shift')

lj.pair_coeff.set('halo', 'halo', epsilon=epsilon_u, sigma=cluster.halo_diam)
# r_cut = False, excludes type pair interaction from neighbour list
# No interaction of halo-core or core-core
lj.pair_coeff.set(['halo', 'core'], 'core', r_cut=False, epsilon=0, sigma=0)

# Apply FENE bonds
FENE = hoomd.md.bond.fene()
FENE.bond_coeff.set('fene', k=fene_k, r0=fene_r0, sigma=cluster.halo_diam,
                    epsilon=epsilon_u)
FENE.bond_coeff.set('harmonic', k=0, r0=0, sigma=0, epsilon=0)
FENE.bond_coeff.set('fene_skip', k=0, r0=0, sigma=0, epsilon=0)

# Apply FENE SKIP bonds
FENE_SKIP = hoomd.md.bond.fene()
FENE_SKIP.bond_coeff.set('fene_skip', k=fene_k, r0=fene_r0*2, sigma=cluster.halo_diam*2,
                         epsilon=epsilon_u)
FENE_SKIP.bond_coeff.set('harmonic', k=0, r0=0, sigma=0, epsilon=0)
FENE_SKIP.bond_coeff.set('fene', k=0, r0=0, sigma=0, epsilon=0)


# Apply Harmonic bonds
HARMONIC = hoomd.md.bond.harmonic()
HARMONIC.bond_coeff.set('harmonic', k=harm_k, r0=(
    cluster.sphere_diam-halo_diam)/2)
HARMONIC.bond_coeff.set('fene', k=0, r0=0)
HARMONIC.bond_coeff.set('fene_skip', k=0, r0=0)

'''
# Apply Hertzian bonds
HERTZIAN = hoomd.md.bond.table(width=1000)
HERTZIAN.bond_coeff.set('hertzian', func=hertzian, rmin=0, rmax=cluster.sphere_diam*2, coeff=dict(U=1000, sigma_H=(cluster.sphere_diam-halo_diam)/2))
HERTZIAN.bond_coeff.set('fene', func=hertzian, rmin=0, rmax=cluster.sphere_diam*2, coeff=dict(U=0, sigma_H=(cluster.sphere_diam-halo_diam)/2))
'''
# # Equilibration
# Integrator selection
hoomd.md.integrate.mode_standard(dt=time_step)

# Creates grooup consisting of central particles in rigid body
all_group = hoomd.group.all()
langevin = hoomd.md.integrate.langevin(
    group=all_group, kT=settings['kT_therm'], seed=settings['seed'])
langevin.set_gamma(a='halo', gamma=settings['fric_coeff'])
langevin.set_gamma(a='core', gamma=0)


'''
halo_group = hoomd.group.type('halo')
core_group = hoomd.group.type('core')

langevin_halo = hoomd.md.integrate.langevin(
        group=halo_group, kT=settings['kT_therm'], seed=settings['seed'])
langevin_halo.set_gamma(a='halo', gamma=settings['fric_coeff'])
langevin_halo.disable()

langevin_core = hoomd.md.integrate.langevin(
    group=core_group, kT=settings['kT_therm'], seed=settings['seed'], noiseless_r=True, noiseless_t=True)
langevin_core.set_gamma(a='core', gamma=0)
langevin_core.disable()
'''
# Adjust density before actual run
vol = cluster.vol_cluster(dimensions) * total_N

print('!!!!!!!!!!!!!!!!!!!!!\nPre-Initial Compression')
print(vol/system.box.get_volume())

pre_dens = vol/system.box.get_volume()
dens = 0.55

if dimensions == 2:
    if dimensions == 2:
        boxLen = math.sqrt(vol / dens)
        hoomd.update.box_resize(Lx=boxLen, Ly=boxLen,
                                period=None, scale_particles=True)
    elif dimensions == 3:
        boxLen = math.pow(vol / dens, 1/3)
        hoomd.update.box_resize(L=boxLen, period=None, scale_particles=True)

elif dimensions == 3:
    for dens in range(int(pre_dens*10000), int(dens*10000), 100):
        dens = dens / 10000
        if dimensions == 2:
            boxLen = math.sqrt(vol / dens)
            hoomd.update.box_resize(Lx=boxLen, Ly=boxLen,
                                    period=None, scale_particles=True)
        elif dimensions == 3:
            boxLen = math.pow(vol / dens, 1/3)
            hoomd.update.box_resize(
                L=boxLen, period=None, scale_particles=True)

        hoomd.run(500, quiet=True)


# Store snapshot information
# hoomd.dump.gsd("final.gsd", group=hoomd.group.all(),
#                overwrite=True, period=None)
# Computes thermodynamical properties of halo particles
halo_group = hoomd.group.type(type='halo', name='halo')
core_group = hoomd.group.type(type='core', name='core')
halo_thermo = hoomd.compute.thermo(group=halo_group)
core_thermo = hoomd.compute.thermo(group=core_group)

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
                                    'pair_lj_energy',
                                    'potential_energy_core',
                                    'potential_energy_halo',
                                    'kinetic_energy_core',
                                    'kinetic_energy_halo'],
                        period=outputInterval_log,
                        overwrite=True)


gsd = hoomd.dump.gsd(filename='{:s}.gsd'.format(nameString),
                     period=outputInterval_gsd,
                     group=hoomd.group.all(),
                     dynamic=['momentum'],
                     overwrite=True)

'''
trap_stdout = io.StringIO()

with redirect_stdout(trap_stdout):
    for i in range(therm_steps):
        langevin_halo.enable()
        hoomd.run(1, quiet=True)
        langevin_halo.disable()

        langevin_core.enable()
        hoomd.run(1, quiet=True)
        langevin_core.disable()
'''

# Increase density

for dens in range(5500, 9910, 10):
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

for dens in range(9900, 5490, -10):
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
