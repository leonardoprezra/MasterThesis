'''
Automated jobs
    - Create randomly placed clusters in simulation box in 2D and 3D
    - Create lattice of clusters in simulation box in 2D and 3D
    - Run MD simulation

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
nve

Arguments
    ----------
    N_cluster : int
        Number of halo spheres in the cluster.
        sys.argv[1]
    
'''

import os
import subprocess
import random
import math
import sys

import argparse

# Command line argument parsing
parser = argparse.ArgumentParser(
    description='Run simulations of rigid 2D clusters.')
parser.add_argument('-n', '--Number-Spheres', type=int, dest='N_cluster',
                    help='number of spheres in the cluster')

args = parser.parse_args()

# General simulation parameters
settings = {}
settings['N'] = 65  # N**2 or N**3 are the number of PSCs
settings['diameter'] = 10  # Diameter of halo particles
settings['epsilon'] = 1.0  # WCA-potential parameters
settings['mass'] = 1.0  # Mass of halo particles
settings['nameString'] = 'integrator-{integrator}_shape-{poly}_N-{N:4d}_VF-{density:4.2f}_dim-{dimensions}_Nclus-{N_cluster}_tstep-{time_step:7.5f}_ratio-{ratio:5.3f}_tmult-{tstep_multiplier:5.3f}_pair-{pair}'
settings["initFile"] = 'None'
# Number of time steps between data storage in gsd file
settings['outputInterval_gsd'] = 70000
# Number of time steps between data storage in log file
settings['outputInterval_log'] = 70
settings['equil_steps'] = 70000  # Number of equilibration steps
settings['ratio'] = 1
settings['tstep_multiplier'] = 0.005
settings['sigma'] = settings['diameter'] * \
    settings['ratio']  # WCA-potential parameters (LANGEVIN)

settings['density'] = 0.70

nameFormat = "data_{poly}/" + settings['nameString']


# Specific simulation parameters
parameterspace = []

tstep_multiplier = settings['tstep_multiplier']

parameterspace += [
    {**settings,
     'integrator': 'langevin',
     'poly': '2Dspheres',
     'dimensions': 2,
     'N_cluster': args.N_cluster,
     'ratio': 1,
     'time_step': tstep_multiplier*math.sqrt(settings['mass']*settings['sigma']**2/settings['epsilon']),
     'pair' : 'tabulated'
     # 'initFile': [nameFormat.format(**settings)+'_restart-000.gsd']
     }]


# Run Simulations
for initDict in parameterspace:

    # Print the correct number of clusters in the system
    user_N = initDict['N']
    if initDict['dimensions'] == 3:
        initDict['N'] = initDict['N']**3
    elif initDict['dimensions'] == 2:
        initDict['N'] = initDict['N']**2

    nameString = nameFormat.format(**initDict)

    initDict['N'] = user_N

    # Create directories
    try:
        if(not os.path.exists("data_{poly}/".format(**initDict))):
            os.mkdir("data_{poly}/".format(**initDict))
    except OSError as e:
        if e.errno != 17:
            raise
        pass

    # Print current working simulation
    if(os.path.exists(nameString+".outputs")):
        print("\nSkipping "+nameString)
        continue
    else:
        print("\nLauing "+nameString)

    # Create list with command line arguments
    initString = []

    for p, v in initDict.items():
        initString.append('{}={}'.format(str(p), str(v)))

    # Run simulations
    out = open(nameString+".outputs", "w")
    proc = subprocess.Popen(["python",  "-u",
                             # "MarsonNVE.py",  #
                             "/home/hpc/iwsp/iwsp023h/MasterThesis/Code/MarsonNVEHisteresisClusters.py",
                             # "/nishome/students/leonardo/Dokumente/Thesis/Code/MarsonNVEHisteresisClusters.py",
                             # "/home/leo/MasterThesis/Code/MarsonNVEHisteresisClusters.py",
                             *initString],
                            stdout=out,
                            stderr=out)
    proc.wait()
