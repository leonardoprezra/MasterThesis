import math
import sys
import numpy as np
from numpy import linalg as LA
import hoomd
import hoomd.md
from scipy.spatial.transform import Rotation as R

# Create snapshot for system initialization


def create_snapshot_soft(cluster, N, dimensions=3):
    '''Create snapshot of soft particles and initializes it.

    Parameters
    ----------
    cluster : PartCluster
        Object that contains cluster information.
    dimensions : int
        Dimensions of the simulation 2D or 3D.
    N : int
        Total number of clusters is N**2 or N**3.

    Returns
    -------
    system
        System_data object of hoomd.
    rigid
        Rigid bodies.
    '''

    # Creates snapshot with first cluster
    snapshot = hoomd.data.make_snapshot(N=1+cluster.N_cluster,
                                        particle_types=['core', 'halo'],
                                        #bond_types=['fene', 'hertzian'],
                                        bond_types=['fene', 'harmonic'],
                                        box=hoomd.data.boxdim(
                                            L=cluster.sphere_diam * 1.1, dimensions=dimensions)
                                        )

    # Properties of core particles
    snapshot.particles.typeid[0] = 0
    snapshot.particles.mass[0] = cluster.halo_mass
    snapshot.particles.position[0] = (0, 0, 0)
    snapshot.particles.diameter[0] = cluster.core_diam*0.1
    snapshot.particles.body[0] = -2

    # Properties of halo particles
    snapshot.particles.typeid[1:] = 1

    for i in range(1, snapshot.particles.N):
        snapshot.particles.mass[i] = cluster.halo_mass
        # i goes one index higher than the lenght of coord
        snapshot.particles.position[i] = cluster.core_coord[i-1]
        snapshot.particles.diameter[i] = cluster.halo_diam
        snapshot.particles.body[i] = -2

    if dimensions == 2:
        # Set FENE bonds among halo particles
        snapshot.bonds.resize(cluster.N_cluster*2)
        snapshot.bonds.group[:cluster.N_cluster] = [
            [i, i+1] for i in range(1, snapshot.particles.N-1)] + [[1, snapshot.particles.N-1]]

        # Set Harmonic bonds among halo and core particles
        # Set Hertzian bonds among halo and core particles
        snapshot.bonds.group[cluster.N_cluster:] = [[0, i]
                                                    for i in range(1, snapshot.particles.N)]
        snapshot.bonds.typeid[cluster.N_cluster:] = 1

    elif dimensions == 3:
        # Set Harmonic bonds among halo and core particles
        # Set Hertzian bonds among halo and core particles
        snapshot.bonds.resize(cluster.N_cluster)
        snapshot.bonds.group[:cluster.N_cluster] = [[0, i]
                                                    for i in range(1, snapshot.particles.N)]
        snapshot.bonds.typeid[:cluster.N_cluster] = 1

        bonds_size = cluster.N_cluster
        # Set FENE bonds among halo particles
        for i1, p1 in enumerate(snapshot.particles.position):
            prev_bonds_size = bonds_size
            pairs = []
            for i2, p2 in enumerate(snapshot.particles.position):

                if i1 != i2 and i1 != 0 and i2 != 0:
                    dist = math.sqrt(sum([(p1[a]-p2[a])**2 for a in range(3)]))

                    if dist <= cluster.halo_diam * 2:
                        pairs.append(i2)

            bonds_size += len(pairs)

            if prev_bonds_size != bonds_size:
                snapshot.bonds.resize(bonds_size)
                snapshot.bonds.group[prev_bonds_size:] = [[i1, i]
                                                          for i in pairs]
                snapshot.bonds.typeid[prev_bonds_size:] = 0

    # Replicates cluster in snapshot

    if dimensions == 2:
        snapshot.replicate(N, N, 1)
    elif dimensions == 3:
        snapshot.replicate(N, N, N)

    # Initialize sys configuration
    system = hoomd.init.read_snapshot(snapshot)

    group_core = hoomd.group.type(type='core')
    total_N = len(group_core)  # total number of clusters
    body_flags = set([p.body for p in system.particles])

    # Set body flags to negative values (floppy bodies)
    i = -3
    for b in body_flags:
        for p in system.particles:
            if p.body == b:
                p.body = i
        i += -1

    print('[II] Snapshot . . . . done.')

    return (system, total_N)


# Hertzian potential

def hertzian(r, rmin, rmax, U, sigma_H):
    '''Create snapshot of soft particles and initializes it.

    Parameters
    ----------
    r : PartCluster
        Object that contains cluster information.
    dimensions : int
        Dimensions of the simulation 2D or 3D.
    N : int
        Total number of clusters is N**2 or N**3.
    density : float
        target density of the system.



    Returns
    -------
    system
        System_data object of hoomd.
    rigid
        Rigid bodies.
    '''
    V = U*(1-r/sigma_H)**(5/2)
    F = -5/2 * U*(1-r/sigma_H)**(3/2)*-1/sigma_H

    return (V, F)
