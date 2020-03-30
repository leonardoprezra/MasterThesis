import math
import sys
import numpy as np
from numpy import linalg as LA
import hoomd
import hoomd.md
from scipy.spatial.transform import Rotation as R

# Parameters
settings = {}

settings['N'] = 6  # N**2 or N**3 are the number of PSCs
settings['diameter'] = 1.0  # Diameter of halo particles
settings['poly'] = 'dode'  # Type of polyhedron
settings['mass'] = 1.0  # Mass of halo particles
settings['density'] = 0.70  # Volume fraction
settings['dimensions'] = 3  # 2d or 3d
settings['N_cluster'] = 9  # number of spheres in cluster

settings['integrator'] = 'nve'  # Integrator
settings['nameString'] = 'integrator-{integrator}_shape-{poly}_N-{N}_VF-{density:4.2f}_dim-{dimensions}_Nclus-{N_cluster}_tstep-{time_step}_feneK-{fene_k}_harmK-{harm_k}'
settings["initFile"] = "None"

settings['max_move'] = 0.002  # Maximum move displacement (HPMC)
settings['max_rot'] = 0.4  # Maximum move rotation (HPMC)
settings['seed'] = 42  # Random number seed (HPMC, LANGEVIN)

settings['sigma'] = settings['diameter']  # WCA-potential parameters (LANGEVIN)
settings['epsilon'] = 1.0  # WCA-potential parameters (LANGEVIN)
settings['kT_therm'] = 5.0  # Temperature of the simulation (LANGEVIN, NPT)
settings['kT_npt'] = 5.0  # Temperature of the simulation (LANGEVIN, NPT)
settings['kT_equil'] = 1.0  # Temperature of the simulation (LANGEVIN, NPT)
visc = 1/math.pi/6/settings['sigma']  # Solvent viscosity (LANGEVIN)
# Particle friction coefficient (LANGEVIN)
settings['fric_coeff'] = 6*math.pi*visc*settings['sigma']

settings['tau'] = 1.0  # Coupling constant for the thermostat (NPT)
settings['pressure'] = 5  # Isotropic pressure set point for barostat (NPT)
tauP = settings['tauP'] = 1.2  # Coupling constant for the barostat (NPT)

settings['hpmc_steps'] = 1000  # Number of time steps of hpmc simulation
settings['npt_steps'] = 1000  # Number of steps required during compression
settings['equil_steps'] = 10000  # Number of equilibration steps
settings['therm_steps'] = 10000  # Number of thermalization steps
settings['nve_steps'] = 1000  # Number of thermalization steps

settings['outputInterval'] = 10  # Number of time steps between data storage
a = math.sqrt(settings['mass']*settings['sigma']**2/settings['epsilon'])
# Time step of MD simulations
settings['time_step'] = 0.001 * \
    math.sqrt(settings['mass']*settings['sigma']**2/settings['epsilon'])

settings['fene_k'] = 15
settings['harm_k'] = 1000


# Core particle properties


def core_properties(d=0, poly=None, mass=0, N_spheres=0):
    '''Information of clusters.

    Parameters
    ----------
    d : float
        Diameter of halo sphere.
    poly : str
        Type of information requested (diam, type, coord, mass, N, sphere_diam).
    mass : float
        Mass of halo sphere.
    N_spheres : int
        Number of halo spheres in the cluster.

    Returns
    -------
    core_values[poly] : float or list
        Value of information requested .
    '''
    g_r = (1 + math.sqrt(5))/2  # Golden ratio

    if N_spheres != 0:
        theta_2d = 2*math.pi/N_spheres
        radius_2d = d/math.sin(theta_2d)*math.sin(math.pi/2 - theta_2d/2)
    else:
        theta_2d = 0
        radius_2d = 0

    if poly[:9] == '3Dspheres':
        coord, core_diam, sphere_diam = unif_pos(N=N_spheres, d=d)
    else:
        coord, core_diam, sphere_diam = [0, 0, 0]

    core_values = {
        'octa_diam': d * (-1+math.sqrt(2)),
        'octa_type': ['halo']*6,
        'octa_coord': [(-d/math.sqrt(2), 0, 0), (d/math.sqrt(2), 0, 0),
                       (0, -d/math.sqrt(2), 0), (0, d/math.sqrt(2), 0),
                       (0, 0, -d/math.sqrt(2)), (0, 0, d/math.sqrt(2))],
        'octa_mass': mass*6,
        'octa_N': 6,
        'octa_sphere_diam': math.sqrt(2)*d + d,

        'tetra_diam': d * math.sqrt(3/8) * (1-1/2),
        'tetra_type': ['halo']*4,
        'tetra_coord': [[d*math.sqrt(3/8) * i for i in (math.sqrt(8/9), 0, -1/3)], [d*math.sqrt(3/8)*i for i in (-math.sqrt(2/9), math.sqrt(2/3), -1/3)],
                        [d*math.sqrt(3/8)*i for i in (-math.sqrt(2/9), -math.sqrt(2/3), -1/3)], [d*math.sqrt(3/8)*i for i in (0, 0, 1)]],
        'tetra_mass': mass*4,
        'tetra_N': 4,
        'tetra_sphere_diam': math.sqrt(3/2)*d + d,

        'cube_diam':  d*(math.sqrt(3)-1),
        'cube_type': ['halo']*8,
        'cube_coord': [(d/2, d/2, d/2), (-d/2, -d/2, -d/2),
                       (d/2, -d/2, d/2), (-d/2, d/2, -d/2),
                       (d/2, d/2, -d/2), (-d/2, -d/2, d/2),
                       (-d/2, d/2, d/2), (d/2, -d/2, -d/2)],
        'cube_mass': mass*8,
        'cube_N': 8,
        'cube_sphere_diam': d*(math.sqrt(3)+1),

        'ico_diam': (math.sqrt(g_r+2)*d/2-d/2)*2,
        'ico_type': ['halo']*12,
        'ico_coord': [(0, d/2, d/2*g_r), (0, -d/2, -d/2*g_r),
                      (0, d/2, -d/2 * g_r), (0, -d/2, d/2*g_r),
                      (d/2*g_r, 0, d/2), (-d/2*g_r, 0, -d/2),
                      (d/2*g_r, 0, -d/2), (-d/2*g_r, 0, d/2),
                      (d/2, d/2*g_r, 0), (-d/2, -d/2*g_r, 0),
                      (d/2, -d/2*g_r, 0), (-d/2, d/2*g_r, 0)],
        'ico_mass': mass*12,
        'ico_N': 12,
        'ico_sphere_diam': (math.sqrt(g_r+2)*d/2+d/2)*2,

        'dode_diam': (math.sqrt(3)*d*g_r/2-d/2)*2,
        'dode_type': ['halo']*20,
        'dode_coord': [(i*d*g_r/2, j*d*g_r/2, k*d*g_r/2) for i in [1, -1] for j in [1, -1] for k in [1, -1]]
        + [(0, i*d*g_r/2, j*d*g_r/2) for i in [g_r, -g_r]
            for j in [1/g_r, -1/g_r]]
        + [(j*d*g_r/2, 0, i*d*g_r/2) for i in [g_r, -g_r]
            for j in [1/g_r, -1/g_r]]
        + [(i*d*g_r/2, j*d*g_r/2, 0) for i in [g_r, -g_r]
            for j in [1/g_r, -1/g_r]],
        'dode_mass': mass*20,
        'dode_N': 20,
        'dode_sphere_diam': (math.sqrt(3)*d*g_r/2+d/2)*2,

        '2Dspheres_type': ['halo']*N_spheres,
        '2Dspheres_coord': [(radius_2d*math.cos(theta_2d*i), radius_2d*math.sin(theta_2d*i), 0) for i in range(N_spheres)],
        '2Dspheres_mass': mass*N_spheres,
        '2Dspheres_N': N_spheres,
        '2Dspheres_diam': (radius_2d-d/2)*2,
        '2Dspheres_sphere_diam': radius_2d*2+d,

        '3Dspheres_type': ['halo']*N_spheres,
        '3Dspheres_coord': coord,
        '3Dspheres_mass': mass*N_spheres,
        '3Dspheres_N': N_spheres,
        '3Dspheres_diam':  core_diam,
        '3Dspheres_sphere_diam': sphere_diam,
    }
    try:
        return core_values[poly]
    except:
        sys.exit('Incorrect polyhedron-key')


# Cluster object
class PartCluster:
    '''Cluster of particles.

    Parameters
    ----------
    poly_key : str
        Type of cluster.
    N_cluster : int
        Number of halo spheres in the cluster.
    halo_diam : float
        Diameter of halo sphere.
    hal_mass: float
        Mass of halo sphere.

    Attributes:
    -----------
    core_diam : float
        Diameter of central sphere.
    core_type : list
        Type of central sphere.
    core_coord : list
        Local coordinates of halo spheres.
    core_mass : float
        Total mass of all halo speres.
    N_cluster : int
        Number of halo spheres in the cluster.
    sphere_diam : float
        Diameter of sphere that circumscribes the cluster.
    halo_diam : float
        Diameter of halo sphere/
    poly_key : str
        Type of cluster.
    rot_matrix : np.array
        Rotation matrix to rotate from local coordinates to coordinates of principal axis.
    halo_mass : float
        Mass of each halo particles

    Methods
    -------
    vol_cluster(dimensions)
        Volume (area in 2D) of the cluster
    '''

    def __init__(self, poly_key, N_cluster, halo_diam, halo_mass):
        self.core_diam = core_properties(poly='{}_diam'.format(
            poly_key), d=halo_diam, N_spheres=N_cluster)  # Diameter of core particles
        self.core_type = core_properties(poly='{}_type'.format(
            poly_key), N_spheres=N_cluster)  # List to create clusters
        self.core_coord = core_properties(poly='{}_coord'.format(
            poly_key), d=halo_diam, N_spheres=N_cluster)  # Coordinates of particles in cluster (global coordinates)
        self.core_mass = core_properties(poly='{}_mass'.format(
            poly_key), mass=halo_mass, N_spheres=N_cluster)  # Mass of core particles
        self.N_cluster = core_properties(poly='{}_N'.format(
            poly_key), N_spheres=N_cluster)  # Number of particles in cluster
        self.sphere_diam = core_properties(poly='{}_sphere_diam'.format(
            poly_key), d=halo_diam, N_spheres=N_cluster)  # Diameter of sphere circunscribing PSC + halo_diam
        self.halo_diam = halo_diam
        self.poly_key = poly_key
        self.halo_mass = halo_mass  # Mass of each halo particles

    def vol_cluster(self, dimensions):
        if dimensions == 2:

            return math.pi/4 * (self.core_diam**2 + self.N_cluster*self.halo_diam**2)
        elif dimensions == 3:
            return math.pi/6 * (self.core_diam**3 + self.N_cluster*self.halo_diam**3)


# Create snapshot for system initialization


def create_snapshot_soft(cluster, N, density, dimensions=3):
    '''Create snapshot of soft particles and initializes it.

    Parameters
    ----------
    cluster : PartCluster
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
    
    elif dimensions ==3:
        # Set Harmonic bonds among halo and core particles
        # Set Hertzian bonds among halo and core particles
        snapshot.bonds.resize(cluster.N_cluster)
        snapshot.bonds.group[:cluster.N_cluster] = [[0, i]
                                                    for i in range(1, snapshot.particles.N)]
        snapshot.bonds.typeid[:cluster.N_cluster] = 1

        bonds_size = cluster.N_cluster
        # Set FENE bonds among halo particles
        for i1,p1 in enumerate(snapshot.particles.position):
            prev_bonds_size = bonds_size
            pairs = []
            for i2,p2 in enumerate(snapshot.particles.position):

                if i1 != i2 and i1 != 0 and i2 !=0:
                    dist = math.sqrt(sum([(p1[a]-p2[a])**2 for a in range(3)]))

                    if dist <= cluster.halo_diam *2:
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


# Distribute points uniformly on a sphere


def unif_pos(N, d):
    """
    Generate almost uniformly distributed positions on the surface of a unit sphere.

    Parameters
    ----------
    N : int > 1
        Number of halo spheres in the cluster.
    d : float
        Diameter of halo sphere.

    Returns
    -------
    points : list
        Coordinates.
    core_diam : float
        Diameter of central sphere.
    sphere_diam : float
        Diameter of sphere that circumscribes the cluster.
    """
    if N == 1:
        return [[1, 0, 0]], [[1, 0, 0, 0]]
    points = []
    offset = 2/N
    increment = np.pi * (3 - np.sqrt(5))
    i = 0
    while i < N:
        y = ((i * offset) - 1) + (offset / 2)
        r = np.sqrt(1 - y**2)
        phi = (i % N) * increment
        x = np.cos(phi) * r
        z = np.sin(phi) * r
        points.append([x, y, z])
        i += 1

    # Find smallest distance between points
    min_dist = math.sqrt((points[0][0]-points[1][0])**2 + (points[0]
                                                           [1]-points[1][1])**2 + (points[0][2]-points[1][2])**2)

    for pi in points:
        for pj in points:
            dist = math.sqrt((pi[0]-pj[0])**2 +
                             (pi[1]-pj[1])**2 + (pi[2]-pj[2])**2)
            if min_dist > dist and dist != 0:
                min_dist = dist

    # Rescale unit sphere
    for i in range(len(points)):
        points[i] = [d/min_dist * b for b in points[i]]

    core_diam = (1*d/min_dist - d/2)*2*0.01
    sphere_diam = (1*d/min_dist + d/2)*2

    return points, core_diam, sphere_diam


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
