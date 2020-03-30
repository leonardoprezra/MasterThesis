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
from GnanFunctions import core_properties, settings, PartCluster, create_snapshot_soft, unif_pos, hertzian

import os
import sys
import time

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

outputInterval = settings['outputInterval']
time_step = settings['time_step']

therm_steps = settings['therm_steps']
npt_steps = settings['npt_steps']
equil_steps = settings['equil_steps']
nve_steps = settings['nve_steps']
N_cluster = settings['N_cluster']

sigma_u = settings['sigma']
epsilon_u = settings['epsilon']

fene_k = settings['fene_k']
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
    poly_key=poly_key, N_cluster=N_cluster, halo_diam=halo_diam, halo_mass=halo_mass)

# Area of clusters
# 0.89 factor was chosen after study in 2D

# Initialize execution context
hoomd.context.initialize("")  # "--mode=gpu"
print('[I] Initialize . . . . done.')

# # Create snapshot for simulation
# Lattice as starting configuration
system, total_N = create_snapshot_soft(
    cluster=cluster, dimensions=dimensions, N=N, density=density)

if dimensions == 2:
    boxLen = math.sqrt(cluster.vol_cluster(dimensions) * total_N / density)
elif dimensions == 3:
    boxLen = math.pow(cluster.vol_cluster(
        dimensions) * total_N / density, 1/3)

# Neighbor list and Potential selection
nl = hoomd.md.nlist.cell()

# WCA is LJ with shift that causes potential to be zero a r_cut
lj = hoomd.md.pair.lj(r_cut=2**(1/6)*sigma_u, nlist=nl)

# Shifts interaction potential, so that its value is zero at the r_cut
lj.set_params(mode='shift')

lj.pair_coeff.set('halo', 'halo', epsilon=epsilon_u, sigma=sigma_u)
lj.pair_coeff.set(['halo', 'core'], 'core', epsilon=0, sigma=0) # No interaction of halo-core or core-core

# Apply FENE bonds
FENE = hoomd.md.bond.fene()
FENE.bond_coeff.set('fene', k=fene_k, r0=5.5, sigma=sigma_u,
                    epsilon=settings['epsilon'])
FENE.bond_coeff.set('harmonic', k=0, r0=5.5, sigma=0, epsilon=0)

# Apply Harmonic bonds
HARMONIC = hoomd.md.bond.harmonic()
HARMONIC.bond_coeff.set('harmonic', k=harm_k , r0=(cluster.sphere_diam-halo_diam)/2)
HARMONIC.bond_coeff.set('fene', k=0 , r0=0)

'''
# Apply Hertzian bonds
HERTZIAN = hoomd.md.bond.table(width=1000)
HERTZIAN.bond_coeff.set('hertzian', func=hertzian, rmin=0, rmax=cluster.sphere_diam*2, coeff=dict(U=1000, sigma_H=(cluster.sphere_diam-halo_diam)/2))
HERTZIAN.bond_coeff.set('fene', func=hertzian, rmin=0, rmax=cluster.sphere_diam*2, coeff=dict(U=0, sigma_H=(cluster.sphere_diam-halo_diam)/2))
'''
# # Thermalization
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

# Store snapshot information
# hoomd.dump.gsd("final.gsd", group=hoomd.group.all(),
#                overwrite=True, period=None)
# Computes thermodynamical properties of halo particles

#halo_thermo = hoomd.compute.thermo(group=group_halo)

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
                                    'pair_lj_energy'],
                        period=outputInterval,
                        overwrite=True)

gsd = hoomd.dump.gsd(filename='{:s}.gsd'.format(nameString),
                     period=outputInterval,
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

hoomd.run(therm_steps)
langevin.disable()

# # Compression

# Retrieve current pressure value
pressure = log.query('pressure')
print('PRESSURE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n{0:.8f}\n'.format(pressure))

npt = hoomd.md.integrate.npt(
    group=all_group, kT=settings['kT_npt'], tau=settings['tau'], P=settings['pressure'], tauP=settings['tauP'])

density_compression = cluster.vol_cluster(
    dimensions)*total_N/system.box.get_volume()
density_compression_control = density_compression


while density_compression < density:
    hoomd.run(10, quiet=True)

    density_compression = cluster.vol_cluster(
        dimensions)*total_N/system.box.get_volume()

npt.disable()
langevin.enable()

# Rescale simulation box
hoomd.update.box_resize(L=boxLen, period=None, scale_particles=False)

print('Final density')
print(cluster.vol_cluster(dimensions)*total_N/system.box.get_volume())

# # Equilibration at final volume fraction


for i in [0]:
    '''
    langevin_core.set_params(kT=i)
    langevin_halo.set_params(kT=i)
    '''
    langevin.set_params(kT=i)

    print('TEMPERATURE!!!!!!!!!!!!!!!!!!\n{}\n'.format(i))

    hoomd.run(equil_steps)
    '''
    for i in range(equil_steps):
        langevin_halo.enable()
        hoomd.run(1, quiet=True)
        langevin_halo.disable()

        langevin_core.enable()
        hoomd.run(1, quiet=True)
        langevin_core.disable()
    '''



end_time = time.time()

sim_time = end_time - start_time

print('TOTAL SIMULATION TIME [s] = \n{}'.format(sim_time))
