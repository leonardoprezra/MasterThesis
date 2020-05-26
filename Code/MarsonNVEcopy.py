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

import hoomd
import hoomd.md
import random
import math
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

therm_steps = settings['therm_steps']
equil_steps = settings['equil_steps']
N_cluster = settings['N_cluster']

sigma_u = settings['diameter']*settings['ratio']
epsilon_u = settings['epsilon']

quantities_log = ['volume',
                  'momentum',
                  'time',
                  'potential_energy',
                  'kinetic_energy',
                  'translational_kinetic_energy',
                  'rotational_kinetic_energy',
                  'temperature',
                  'pressure',
                  'pair_lj_energy']

restart_step = settings['outputInterval_log']

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

# Check if init file exists and print simulation information
if(os.path.exists('{:s}_INIT-EQUIL.gsd'.format(nameString))):
    init_file = True
    print('Reading init file: {:s}_INIT-EQUIL.gsd'.format(nameString))
    print('Working on: \n{:s}_EQUIL.log\nUsing: {}'.format(
        nameString, __file__))
elif(not os.path.exists('{:s}_INIT-EQUIL.gsd'.format(nameString))):
    init_file = False
    print('Init file does not exist')
    print('Working on: \n{:s}_INIT-EQUIL.gsd\n{:s}_EQUIL.log\nUsing: {}'.format(
        nameString, nameString, __file__))

print("== using these settings ==")
for k, v in settings.items():
    print("  {:15s}  =  {:}".format(k, v))

# Core particle properties
cluster = PartCluster(
    poly_key=poly_key, N_cluster=N_cluster, halo_diam=halo_diam, halo_mass=halo_mass, ratio=settings['ratio'])

# Initialize execution context
hoomd.context.initialize("")  # "--mode=gpu"
print('[I] Initialize . . . . done.')

# # # If init_file does not exist
if init_file == False:
    # # Create snapshot for simulation
    # Lattice as starting configuration
    system, rigid, group_core, group_halo, total_N = create_snapshot(
        cluster=cluster, dimensions=dimensions, N=N)

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

    # # Thermalization
    # Integrator selection
    hoomd.md.integrate.mode_standard(dt=time_step)
    # creates grooup consisting of central particles in rigid body
    rigid = hoomd.group.rigid_center()
    langevin = hoomd.md.integrate.langevin(
        group=rigid, kT=settings['kT_therm'], seed=settings['seed'])

    langevin.set_gamma(a='halo', gamma=settings['fric_coeff'])
    langevin.set_gamma(a='core', gamma=settings['fric_coeff'])

    # Run simuations to catch possible termination from exceeded wall time
    try:
        hoomd.run_upto(therm_steps, limit_multiple=restart_step)

    # # Compression
        langevin.disable()

        # Retrieve current pressure value
        pressure = log.query('pressure')
        print(
            'PRESSURE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n{0:.8f}\n'.format(pressure))

        npt = hoomd.md.integrate.npt(
            group=rigid, kT=settings['kT_npt'], tau=settings['tau'], P=settings['pressure'], tauP=settings['tauP'])

        density_compression = cluster.vol_cluster(
            dimensions)*total_N/system.box.get_volume()
        density_compression_control = density_compression

        while density_compression < density:
            hoomd.run_upto(10, quiet=True)
            density_compression = cluster.vol_cluster(
                dimensions)*total_N/system.box.get_volume()

        if poly_key == '2Dspheres' and N_cluster == 3:
            hoomd.update.box_resize(Ly=boxLen*boxLen/system.box.Lx,
                                    period=None, scale_particles=False)
        else:
            hoomd.update.box_resize(
                L=boxLen, period=None, scale_particles=False)

        print('Final density')
        print(cluster.vol_cluster(dimensions)*total_N/system.box.get_volume())

        gsd = hoomd.dump.gsd(filename='{:s}_INIT-EQUIL.gsd'.format(nameString),
                             period=None,
                             group=hoomd.group.all(),
                             dynamic=['momentum'],
                             truncate=True,
                             overwrite=True)
        npt.disable()

        # # Equilibration at final volume fraction
        langevin.enable()

        log = hoomd.analyze.log(filename='{:s}_EQUIL.log'.format(nameString),
                                quantities=quantities_log,
                                period=outputInterval_log,
                                overwrite=True)

        gsd = hoomd.dump.gsd(filename='{:s}_EQUIL.gsd'.format(nameString),
                             period=outputInterval_gsd,
                             group=hoomd.group.all(),
                             dynamic=['momentum'],
                             overwrite=True)

        langevin.set_params(kT=settings['kT_equil'])
        hoomd.run_upto(equil_steps, limit_multiple=restart_step)


# # # If init_file exist
elif init_file == True:

    # # Equilibration at final volume fraction
    langevin.enable()

    log = hoomd.analyze.log(filename='{:s}_EQUIL.log'.format(nameString),
                            quantities=quantities_log,
                            period=outputInterval_log,
                            overwrite=True)

    gsd = hoomd.dump.gsd(filename='{:s}_EQUIL.gsd'.format(nameString),
                         period=outputInterval_gsd,
                         group=hoomd.group.all(),
                         dynamic=['momentum'],
                         overwrite=True)

    langevin.set_params(kT=settings['kT_equil'])
    hoomd.run(equil_steps)


end_time = time.time()

sim_time = end_time - start_time

print('TOTAL SIMULATION TIME [s] = \n{}'.format(sim_time))
