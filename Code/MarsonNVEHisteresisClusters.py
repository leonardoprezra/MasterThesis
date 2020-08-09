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
N_cluster = settings['N_cluster']

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

# Core particle properties
cluster = PartCluster(
    poly_key=poly_key, N_cluster=N_cluster, halo_diam=halo_diam, halo_mass=halo_mass, ratio=settings['ratio'], dimensions=dimensions)

print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print('halo_diam={}'.format(cluster.halo_diam))
print('core_diam={}'.format(cluster.core_diam))
print('sphere_diam={}'.format(cluster.sphere_diam))
print(cluster.core_coord)

# Updates parameters
settings['diameter'] = cluster.sphere_diam

if settings['pair'] == 'tabulated':
    settings['sigma'] = cluster.sphere_diam
    sigma_u = cluster.sphere_diam
elif settings['pair'] == 'LJ':
    settings['sigma'] = cluster.halo_diam
    sigma_u = cluster.halo_diam

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
system, rigid, group_core, group_halo, total_N = create_snapshot(
    cluster=cluster, dimensions=dimensions, N=N)

# Adjust density
vol = cluster.vol_cluster(dimensions) * total_N

dens = 0.55
if dimensions == 2:
    boxLen = math.sqrt(vol / dens)
    hoomd.update.box_resize(Lx=boxLen, Ly=boxLen,
                            period=None, scale_particles=True)
elif dimensions == 3:
    boxLen = math.pow(vol / dens, 1/3)
    hoomd.update.box_resize(L=boxLen, period=None, scale_particles=True)

# Neighbor list and Potential selection
nl = hoomd.md.nlist.cell()

# Choose pair potential

# Table potential
if settings['pair'] == 'tabulated':
    r_cut = (2 ** (1/6)*sigma_u)
    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    print('r_cut={}'.format(r_cut))
    # Table potential
    offset = (cluster.halo_diam/2+cluster.core_diam/2)
    r_vis = r_cut-offset

    print('offset={}'.format(offset))
    print('visible_radius={}'.format(r_cut-offset))

    #n_pos = math.pi / (2*math.asin(cluster.halo_diam / (2*r_vis)))
    #corr = int(n_pos)*N_cluster

    # print('n_pos={}'.format(n_pos))
    dict_coeff = dict(epsilon=epsilon_u, sigma=sigma_u,
                      offset=offset)

    table = hoomd.md.pair.table(width=100, nlist=nl)
    table.pair_coeff.set('halo', 'halo', func=WCA_corrected,
                         rmin=0, rmax=r_vis, coeff=dict_coeff)
    table.pair_coeff.set(['halo', 'core'], 'core',
                         func=WCA_corrected, rmin=0, rmax=0.002, coeff=dict_coeff)

# LJ shifted potential
elif settings['pair'] == 'LJ':
    # WCA is LJ with shift that causes potential to be zero a r_cut
    lj = hoomd.md.pair.lj(r_cut=2**(1/6)*sigma_u, nlist=nl)

    # Shifts interaction potential, so that its value is zero at the r_cut
    lj.set_params(mode='shift')
    lj.pair_coeff.set('halo', 'halo', epsilon=epsilon_u,
                      sigma=sigma_u)

    if settings['ratio'] == 1:
        lj.pair_coeff.set(['halo', 'core'], 'core',
                          r_cut=False, epsilon=0, sigma=0)
    elif settings['ratio'] < 1 and settings['ratio'] > 0:
        sigma_core_core = cluster.core_diam
        lj.pair_coeff.set('core', 'core', r_cut=2**(1/6) *
                          sigma_core_core, epsilon=epsilon_u, sigma=sigma_core_core)

        sigma_halo_core = cluster.core_diam/2 + cluster.halo_diam/2
        lj.pair_coeff.set('halo', 'core', r_cut=2**(1/6) *
                          sigma_halo_core, epsilon=epsilon_u, sigma=sigma_halo_core)


# # Equilibration
# Integrator selection
hoomd.md.integrate.mode_standard(dt=time_step)
# creates grooup consisting of central particles in rigid body
rigid = hoomd.group.rigid_center()

if settings['integrator'] == 'langevin':
    langevin = hoomd.md.integrate.langevin(
        group=rigid, kT=settings['kT_equil'], seed=settings['seed'])

    langevin.set_gamma(a='halo', gamma=settings['fric_coeff'])
    langevin.set_gamma(a='core', gamma=settings['fric_coeff'])

elif settings['integrator'] == 'nve':
    nve = hoomd.md.integrate.nve(group=rigid)

    nve.randomize_velocities(kT=settings['kT_equil'], seed=settings['seed'])


# Store snapshot information
if settings['pair'] == 'tabulated':
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
                            period=outputInterval_log,
                            overwrite=True)
elif settings['pair'] == 'LJ':
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
                            period=outputInterval_log,
                            overwrite=True)

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
