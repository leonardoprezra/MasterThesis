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
2Dspheres

Integrators:
local
langevin
'''

from __future__ import division
import hoomd
import hoomd.md
import hoomd.hpmc
import random
import math
from MarsonFunctions import core_properties, mom_inertia, settings

import os
import sys

# Create directory to store simulation results
if(not os.path.exists("data_start2d")):
    os.mkdir("data_start2d")

# Parameters

N = 20
halo_diam = 0.5
halo_mass = 1
density = 0.60
dimensions = 2

outputInterval = 1
time_step = 0.005/2

therm_steps = 1000
npt_steps = 1000
equil_steps = 1000

kT = 1.0

poly_key = '2Dspheres'
N_cluster = 3


core_diameter = core_properties(poly='{}_diam'.format(
    poly_key), d=halo_diam, N_spheres=N_cluster)  # Diameter of core particles
core_type = core_properties(poly='{}_type'.format(
    poly_key), N_spheres=N_cluster)  # List for creation of clusters
core_coord = core_properties(poly='{}_coord'.format(
    poly_key), d=halo_diam, N_spheres=N_cluster)  # Coordinates of particles in cluster
core_mass = core_properties(poly='{}_mass'.format(
    poly_key), mass=halo_mass, N_spheres=N_cluster)  # Mass of core particles
N_cluster = core_properties(poly='{}_N'.format(
    poly_key), N_spheres=N_cluster)  # Number of particles in cluster
sphere_diam = core_properties(poly='{}_sphere_diam'.format(
    poly_key), d=halo_diam, N_spheres=N_cluster)  # Diameter of sphere circunscribing PSC + halo_diam

# boxLen = math.sqrt(area_particle*N*N/density)

area_particle = math.pi/4 * (N_cluster*halo_diam**2 + core_diameter**2)

boxLen = math.sqrt(
    (N*N*area_particle)/density)

print('boxLen')
print(boxLen*boxLen)

# Initialize execution context
hoomd.context.initialize("")
print('[I] Initialize . . . . done.')

# Creates lattice
'''
uc = hoomd.lattice.unitcell(N=1,
                            a1=[1, 0, 0],
                            a2=[0, 1, 0],
                            a3=[0, 0, 1],
                            dimensions=dimensions,
                            position=[[0, 0, 0]],
                            type_name='core',
                            mass=halo_mass,
                            diameter=halo_diam
                            )'''
uc = hoomd.lattice.unitcell(N=1,
                            a1=[1.5*sphere_diam, 0, 0],
                            a2=[0, 1.5*sphere_diam, 0],
                            a3=[0, 0, 1],
                            dimensions=dimensions,
                            position=[[0, 0, 0]],
                            type_name=['core'],
                            mass=[core_mass],
                            diameter=[core_diameter],
                            moment_inertia=[[1.0, 1.0, 1.0]],
                            orientation=[[0.707, 0, 0, 0.707]])

# Initialize system configuration
system = hoomd.init.create_lattice(unitcell=uc, n=[N, N])
system.particles.types.add('halo')
# Create rigid clusters
rigid = hoomd.md.constrain.rigid()

rigid.set_param('core',
                types=core_type,
                positions=core_coord)  # Set constituent particles of the rigid body

rigid.create_bodies()
rigid.validate_bodies()
# Set particle diameters
group_halo = hoomd.group.type('halo')
group_core = hoomd.group.type('core')

for p in group_halo:
    p.diameter = halo_diam

print('particles diameter')
print(system.particles[0].diameter)
print('paticles number!!!!!')
print(system.particles)
d = hoomd.dump.gsd(filename='./data_start2d/lattice.gsd',
                   period=None,
                   group=hoomd.group.all(),
                   overwrite=True)

print('[II] Snapshot . . . . done.')

# Neighbor list and Potential selection
nl = hoomd.md.nlist.cell()

# WCA is LJ with shift that causes potential to be zero a r_cut
sigma_u = halo_diam
epsilon_u = 1.0

lj = hoomd.md.pair.lj(r_cut=2**(1/6)*sigma_u, nlist=nl)

# Shifts interaction potential, so that its value is zero at the r_cut
lj.set_params(mode='shift')

lj.pair_coeff.set(['halo', 'core'], ['halo', 'core'],
                  epsilon=epsilon_u, sigma=sigma_u)
#lj.pair_coeff.set(['core'], 'core', epsilon=0, sigma=0)

# lj.pair_coeff.set('core', 'core', epsilon=epsilon_u, sigma=sigma_u)


# # Thermalization
# Integrator selection
standard_integration = hoomd.md.integrate.mode_standard(dt=time_step)
rigid = hoomd.group.rigid_center()
# creates grooup consisting of central particles in rigid body
all = hoomd.group.all()
langevin = hoomd.md.integrate.langevin(
    group=rigid, kT=kT, seed=123)


l = hoomd.analyze.log('./data_start2d/thermalization.log',
                      quantities=['potential_energy',
                                  'translational_kinetic_energy',
                                  'rotational_kinetic_energy',
                                  'volume',
                                  'pressure',
                                  'temperature'],
                      period=outputInterval,
                      overwrite=True)


d = hoomd.dump.gsd('./data_start2d/thermalization.gsd',
                   period=outputInterval,
                   group=hoomd.group.all(),
                   overwrite=False)


hoomd.run(therm_steps)

# # Compression
langevin.disable()

time_step_compression = time_step

standard_integration.set_params(dt=time_step_compression)

npt = hoomd.md.integrate.npt(
    group=rigid, kT=kT, tau=1.0, P=10, tauP=1.2)

l = hoomd.analyze.log('./data_start2d/compression.log',
                      quantities=['potential_energy',
                                  'translational_kinetic_energy',
                                  'rotational_kinetic_energy',
                                  'volume',
                                  'pressure',
                                  'temperature'],
                      period=outputInterval,
                      overwrite=True)


d = hoomd.dump.gsd('./data_start2d/compression.gsd',
                   period=outputInterval,
                   group=hoomd.group.all(),
                   overwrite=True)

# hoomd.run(npt_steps)

density_compression = (N*N*area_particle)/system.box.get_volume()
density_compression_control = density_compression

while density_compression < density:
    '''
    if density_compression > density_compression_control*2:
        time_step_compression *= 0.7
        density_compression_control = density_compression
        standard_integration.set_params(dt=time_step_compression)
    '''
    hoomd.run(500, quiet=True)
    density_compression = (N*N*area_particle)/system.box.get_volume()
    print('{} > {} || {} > {}'.format(density_compression,
                                      density, boxLen*boxLen, system.box.get_volume()))


hoomd.update.box_resize(Lx=boxLen, Ly=boxLen,
                        period=None, scale_particles=False)
print('Final density')
print((N*N*area_particle)/system.box.get_volume())

print('!!!!!ParticlesDiameterCORE!!!!!')
sadas = []
for p in group_core:
    sadas.append(p.diameter)

print(sadas)

print('!!!!!ParticlesDiameterHALO!!!!!')
sadas = []
for p in group_halo:
    sadas.append(p.diameter)

print(sadas)

# # Equilibration at final volume fraction
npt.disable()
langevin.enable()

standard_integration.set_params(dt=time_step)

l = hoomd.analyze.log('./data_start2d/equilibration.log',
                      quantities=['potential_energy',
                                  'translational_kinetic_energy',
                                  'rotational_kinetic_energy'],
                      period=outputInterval,
                      overwrite=True)


d = hoomd.dump.gsd('./data_start2d/equilibration.gsd',
                   period=outputInterval,
                   group=hoomd.group.all(),
                   overwrite=True)

hoomd.run(equil_steps)
