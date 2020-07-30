import math
import sys
import numpy as np
from numpy import linalg as LA
import hoomd
import hoomd.md
from scipy.spatial.transform import Rotation as R

# Parameters
settings = {}

settings['N'] = 2  # N**2 or N**3 are the number of PSCs
settings['diameter'] = 1  # Outer diameter of cluster
settings['poly'] = '3Dspheres'  # Type of polyhedron
settings['mass'] = 1.0  # Mass of halo particles
settings['density'] = 0.70  # Volume fraction
settings['dimensions'] = 3  # 2d or 3d
settings['N_cluster'] = 8  # number of spheres in cluster
settings['ratio'] = 1.0 # halo_diam/halo_edge


settings['integrator'] = 'langevin'  # Integrator
settings['nameString'] = 'integrator-{integrator}_shape-{poly}_N-{N:4d}_VF-{density:4.2f}_dim-{dimensions}_Nclus-{N_cluster}_tstep-{time_step:7.5f}_ratio-{ratio:5.3f}_tmult-{tstep_multiplier:5.3f}_pair-{pair}'
settings["initFile"] = "None"
settings['pair'] = 'LJ' # Pair potential used: tabulated or LJ

settings['seed'] = 42  # Random number seed (HPMC, LANGEVIN)

settings['sigma'] = settings['diameter']*settings['ratio']  # WCA-potential parameters (LANGEVIN)
settings['epsilon'] = 1.0  # WCA-potential parameters (LANGEVIN)
settings['kT_therm'] = 1.0  # Temperature of the simulation (LANGEVIN, NPT)
settings['kT_npt'] = 1.0  # Temperature of the simulation (LANGEVIN, NPT)
settings['kT_equil'] = 1.0  # Temperature of the simulation (LANGEVIN, NPT)
visc = 1/math.pi/6/settings['sigma'] # Solvent viscosity (LANGEVIN)
settings['fric_coeff'] = 6*math.pi*visc*settings['sigma'] # Particle friction coefficient (LANGEVIN)

settings['tau'] = 1.0  # Coupling constant for the thermostat (NPT)
settings['pressure'] = 50  # Isotropic pressure set point for barostat (NPT)
tauP = settings['tauP'] = 1.2  # Coupling constant for the barostat (NPT)

settings['equil_steps'] = 50  # Number of equilibration steps

settings['outputInterval_gsd'] = 50 # Number of time steps between data storage in gsd file
settings['outputInterval_log'] = 5 # Number of time steps between data storage in log file
settings['tstep_multiplier'] = 0.005
settings['time_step'] = settings['tstep_multiplier']*math.sqrt(settings['mass']*settings['sigma']**2/settings['epsilon'])  # Time step of MD simulations

settings['fene_k'] = 15
settings['fene_r0'] = 1.5
settings['harm_k'] = 500



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


    # Makes every cluster with sphere_diam = 1
    if poly.split('_')[0] == '2Dspheres':
        const = 1/math.sin(theta_2d)*math.sin(math.pi/2 - theta_2d/2)
        radius_2d = d*const/(1+2*const)
        d = d - 2*radius_2d
        core_values = {
            '2Dspheres_type': ['halo']*N_spheres,
            '2Dspheres_coord': [(radius_2d*math.cos(theta_2d*i), radius_2d*math.sin(theta_2d*i), 0) for i in range(N_spheres)],
            '2Dspheres_mass': mass*N_spheres,
            '2Dspheres_N': N_spheres,
            '2Dspheres_diam': (radius_2d-d/2)*2,
            '2Dspheres_sphere_diam': radius_2d*2+d,
            '2Dspheres_halo_diam': d
        }
    else :
        core_values = {
            'octa_diam': d * (-1+math.sqrt(2)),
            'octa_type': ['halo']*6,
            'octa_coord': [(-d/math.sqrt(2), 0, 0), (d/math.sqrt(2), 0, 0),
                        (0, -d/math.sqrt(2), 0), (0, d/math.sqrt(2), 0),
                        (0, 0, -d/math.sqrt(2)), (0, 0, d/math.sqrt(2))],
            'octa_mass': mass*6,
            'octa_N': 6,
            'octa_sphere_diam': math.sqrt(2)*d + d,
            'octa_halo_diam': d,

            'tetra_diam': d * math.sqrt(3/8) * (1-1/2),
            'tetra_type': ['halo']*4,
            'tetra_coord': [[d*math.sqrt(3/8) * i for i in (math.sqrt(8/9), 0, -1/3)], [d*math.sqrt(3/8)*i for i in (-math.sqrt(2/9), math.sqrt(2/3), -1/3)],
                            [d*math.sqrt(3/8)*i for i in (-math.sqrt(2/9), -math.sqrt(2/3), -1/3)], [d*math.sqrt(3/8)*i for i in (0, 0, 1)]],
            'tetra_mass': mass*4,
            'tetra_N': 4,
            'tetra_sphere_diam': math.sqrt(3/2)*d + d,
            'tetra_halo_diam': d,

            'cube_diam':  d*(math.sqrt(3)-1),
            'cube_type': ['halo']*8,
            'cube_coord': [(d/2, d/2, d/2), (-d/2, -d/2, -d/2),
                        (d/2, -d/2, d/2), (-d/2, d/2, -d/2),
                        (d/2, d/2, -d/2), (-d/2, -d/2, d/2),
                        (-d/2, d/2, d/2), (d/2, -d/2, -d/2)],
            'cube_mass': mass*8,
            'cube_N': 8,
            'cube_sphere_diam': d*(math.sqrt(3)+1),
            'cube_halo_diam': d,

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
            'ico_halo_diam': d,

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
            'dode_halo_diam': d,

            'Dspheres_type': ['halo']*N_spheres,
            'Dspheres_coord': np.array([(radius_2d*math.cos(theta_2d*i), radius_2d*math.sin(theta_2d*i), 0) for i in range(N_spheres)]),
            'Dspheres_mass': mass*N_spheres,
            'Dspheres_N': N_spheres,
            'Dspheres_diam': (radius_2d-d/2)*2,
            'Dspheres_sphere_diam': radius_2d*2+d,
            'Dspheres_halo_diam': d,

            '3Dspheres_type': ['halo']*N_spheres,
            '3Dspheres_coord': coord,
            '3Dspheres_mass': mass*N_spheres,
            '3Dspheres_N': N_spheres,
            '3Dspheres_diam':  core_diam,
            '3Dspheres_sphere_diam': sphere_diam,
            '3Dspheres_halo_diam': d,

            'one_type': ['core'],
            'one_coord': [(0,0,0)],
            'one_mass': mass,
            'one_N': 1,
            'one_diam':  d,
            'one_sphere_diam': d,
            'one_halo_diam': d,
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
    ratio : float
        halo_diam/distance-between-contiguous-particles

    Attributes:
    -----------
    core_diam : float
        Diameter of central sphere.
    core_type : list
        Type of central sphere.
    core_coord : list
        Local coordinates of halo spheres.
    core_mass : float
        Total mass of all halo spheres.
    N_cluster : int
        Number of halo spheres in the cluster.
    sphere_diam : float
        Diameter of sphere that circumscribes the cluster.
    halo_diam : float
        Diameter of halo sphere/
    poly_key : str
        Type of cluster.
    inertia: np.array
        Diagonal moment of inertia of the cluster.
    rot_matrix : np.array
        Rotation matrix to rotate from local coordinates to coordinates of principal axis.
    halo_mass: float
        Mass of a single halo sphere.

    Methods
    -------
    vol_cluster(dimensions)
        Volume (area in 2D) of the cluster
    '''

    def __init__(self, poly_key, N_cluster, halo_diam, halo_mass, ratio=1):
        self.core_type = core_properties(poly='{}_type'.format(
            poly_key), d=halo_diam, N_spheres=N_cluster)  # List to create clusters
        self.core_coord = core_properties(poly='{}_coord'.format(
            poly_key), d=halo_diam, N_spheres=N_cluster)  # Coordinates of particles in cluster (global coordinates)
        self.core_mass = core_properties(poly='{}_mass'.format(
            poly_key), mass=halo_mass, d=halo_diam, N_spheres=N_cluster)  # Mass of core particles
        self.N_cluster = core_properties(poly='{}_N'.format(
            poly_key), d=halo_diam, N_spheres=N_cluster)  # Number of particles in cluster
        self.halo_diam = core_properties(poly='{}_halo_diam'.format(
            poly_key), d=halo_diam, N_spheres=N_cluster) * ratio # Number of particles in cluster
        self.core_diam = core_properties(poly='{}_diam'.format(
            poly_key), d=halo_diam, N_spheres=N_cluster) + self.halo_diam*(1-ratio) # Diameter of core particles
        self.sphere_diam = self.core_diam + 2*self.halo_diam  # Diameter of sphere circunscribing PSC
        self.poly_key = poly_key
        self.inertia, self.rot_matrix, self.quaternion = mom_inertia(
            particles=self.core_coord, mass=halo_mass) # Moment of inertia (given as diagonal matrix),
                                                       # Rotation matrix to rotate from global coordinates to principal axes coordinates
        self.core_coord = quat_rotation(particles=self.core_coord, q=self.quaternion) # Coordinates of particles in cluster (principal axes coordinates)
        self.halo_mass = halo_mass
        
        if poly_key == '2Dspheres':
            # Corrects size to have sphere_diam = 1
            self.core_coord = [(coord[0]/self.sphere_diam,coord[1]/self.sphere_diam, coord[2]/self.sphere_diam ) for coord in self.core_coord]
            self.halo_diam = self.halo_diam/self.sphere_diam
            self.core_diam = self.core_diam/self.sphere_diam
            self.sphere_diam = self.sphere_diam/self.sphere_diam
            self.inertia, self.rot_matrix, self.quaternion = mom_inertia(
                particles=self.core_coord, mass=halo_mass)
        

    def vol_cluster(self, dimensions):
        if dimensions == 2:
            return math.pi/4 * (self.core_diam**2 + self.N_cluster*self.halo_diam**2)
        elif dimensions == 3:
            return math.pi/6 * (self.core_diam**3 + self.N_cluster*self.halo_diam**3)


# Moment of Inertia function


def mom_inertia(particles, mass):
    '''Calculates moment of inertia of rigid bodies.

    Parameters
    ----------
    particles : list
        Coordinates of particles in the cluster.
    mass : float
        Mass of halo sphere.

    Returns
    -------
    mom_inertia_princ : np.array
        Principal moment of inertia of the cluster.
    rot_matrix : np.array
        Rotation matrix to rotate from local coordinates to coordinates of principal axis.
    '''
    Ixx = 0.0
    Iyy = 0.0
    Izz = 0.0
    Ixy = 0.0
    Iyz = 0.0
    Ixz = 0.0

    for coord in particles:
        Ixx += (coord[1]**2 + coord[2]**2) * mass
        Iyy += (coord[0]**2 + coord[2]**2) * mass
        Izz += (coord[0]**2 + coord[1]**2) * mass
        Ixy += -coord[0]*coord[1]*mass
        Iyz += -coord[1]*coord[2]*mass
        Ixz += -coord[0]*coord[2]*mass

    inertia_matrix = np.array([[Ixx, Ixy, Ixz],
                               [Ixy, Iyy, Iyz],
                               [Ixz, Iyz, Izz]])

    mom_inertia_princ, princ_axis = LA.eig(inertia_matrix)

    rot_matrix = np.transpose(princ_axis) 
    #rot_matrix = princ_axis # Use to rotate on the opposite direction

    # Construct quaternion
    m = rot_matrix

    q0 = (1/2)*math.sqrt(m[0,0] + m[1,1] + m[2,2] +1)
    q1 = (m[1,2] - m[2,1])/(4*q0)
    q2 = (m[2,0] - m[0,2])/(4*q0)
    q3 = (m[0,1] - m[1,0])/(4*q0)

    quaternion = np.array([q0,q1,q2,q3])

    return mom_inertia_princ, rot_matrix, quaternion

# Create snapshot for system initialization


def create_snapshot(cluster, N, dimensions=3):
    '''Create snapshot from lattice and initializes it.

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
    part_per_cell = 1
    part_position = [[0, 0, 0]]
    part_orientation = [[1, 0, 0, 0]]

    n = [N, N, N]

    a1, a2, a3 = [[1.1*cluster.sphere_diam, 0, 0],
                  [0, 1.1*cluster.sphere_diam, 0],
                  [0, 0, 1.1*cluster.sphere_diam]]

    if dimensions == 2:
        a1, a2, a3 = [[1.1*cluster.sphere_diam, 0, 0],
                      [0, 1.1*cluster.sphere_diam, 0],
                      [0, 0, 1]]

        n = [N, N]
    '''
    if cluster.poly_key == '2Dspheres' and cluster.N_cluster == 3:
        a1, a2, a3 = [[3*cluster.halo_diam, 0, 0],
                      [0, math.sin(math.pi/3)*cluster.halo_diam*2, 0],
                      [0, 0, 1]]
        part_per_cell = 2
        part_position = [[-3*cluster.halo_diam/4, -cluster.halo_diam*math.sin(math.pi/3)/6, 0],
                         [3*cluster.halo_diam/4, cluster.halo_diam*math.sin(math.pi/3)/6, 0]]
        part_orientation = [[-math.cos(math.pi/6/2), 0, 0, math.sin(math.pi/6/2)],
                            [-math.cos(math.pi/2/2), 0, 0, math.sin(math.pi/2/2)]]
    '''
    # Creates lattice
    uc = hoomd.lattice.unitcell(N=part_per_cell,
                                a1=a1,
                                a2=a2,
                                a3=a3,
                                dimensions=dimensions,
                                position=part_position,
                                type_name=part_per_cell*['core'],
                                mass=part_per_cell*[cluster.core_mass],
                                diameter=part_per_cell*[cluster.core_diam],
                                moment_inertia=part_per_cell*[cluster.inertia],
                                orientation=part_orientation)

    # Initialize sys configuration
    system = hoomd.init.create_lattice(unitcell=uc, n=n)

    total_N = len(system.particles)  # total number of clusters

    system.particles.types.add('halo')

    
    # Create rigid clusters
    rigid = hoomd.md.constrain.rigid()

    rigid.set_param('core',
                    types=cluster.core_type,
                    positions=cluster.core_coord)  # Set constituent particles of the rigid body

    rigid.create_bodies()
    rigid.validate_bodies()

    # Set particle diameters
    group_halo = hoomd.group.type('halo')
    group_core = hoomd.group.type('core')

    for p in group_halo:
        p.diameter = cluster.halo_diam
    
    print('[II] Snapshot . . . . done.')

    return (system, rigid, group_core, group_halo, total_N)

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

# Rotate points with rotation matrix
def quat_rotation(particles, q):
    '''Rotate point coordinates according to given rotation matrix.

    Parameters
    ----------
    particles : list
        Coordinates of particles in the cluster.
    q : np.array
        Quaternion used for rotation.

    Returns
    -------
    new_particles : list
        Rotated coordinates.
    '''
    #qnorm = LA.norm(q)
    q0 = q[0]
    q = q[1:]

    particles = np.array(particles)
    new_particles = []

    for p in particles:
        v = p
        # rot_p = np.matmul(rot_matrix,p)
        qpq = ( 2*q0**2 - 1 )*v + 2*np.dot(q,v)*q + 2*q0*np.cross(q,v)

        rot_p = qpq

        new_particles.append(tuple(rot_p.tolist()))

    return new_particles

# WCA potential field corrected for position of halo particle
def WCA_corrected(r, rmin, rmax, epsilon, sigma, offset, corr=1):
    '''WCA potential corrected for offset of halo particles with regards to the center of the cluster.

    Parameters
    ----------
    r : double
        Positions relative to the center of the particle.
    rmin : double
        Below this distance potential is set to zero.
    rmax : double
        Above this distance potential is set to zero.
    epsilon : double
        WCA parameter.
    sigma : double
        WCA parameter.
    offset : double
        Distance from center of cluster to center of halo particle.
    corr : double
        Corrects for the increased interaction between particles.

    Returns
    -------
    V : double
        Potential evaluated at r.
    F : double
        Force evaluated at r.
    '''

    r = r + offset
    V = (4 * epsilon * ((sigma / r)**12 - (sigma / r)**6) + epsilon) / corr
    F = 4 * epsilon / r * (12 * (sigma / r)**12 - 6 * (sigma / r)**6) / corr
    return (V, F)



if __name__ == '__main__':
    cluster = PartCluster(
    poly_key='2Dspheres', N_cluster=3, halo_diam=1, halo_mass=1, ratio=0.6)
    print('\nCOORDINATES')
    print(cluster.core_coord)

    radius = [LA.norm(p) for p in cluster.core_coord]
    print('\nRADIUS CIRCLE')
    print(np.array(radius))

    print('\nCORE DIAMETER')
    print(cluster.core_diam)
    print('\nHALO DIAMETER')
    print(cluster.halo_diam)

