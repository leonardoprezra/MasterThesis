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

Integrators:
local
langevin
nve
'''

import os
import subprocess
import random
import math
import sys
import argparse

# Command line argument parsing
parser = argparse.ArgumentParser(
    description='Run MD simulations of soft particles in 2D and 3D')
parser.add_argument('density', type=int,
                    help='volume fraction covered by clusters')
parser.add_argument('start_N_cluster', type=int,
                    help='smallest number of halo particles in a cluster')
parser.add_argument('end_N_cluster', type=int,
                    help='largest number of halo particles in a cluster')

args = parser.parse_args()


# General simulation parameters
settings = {}
settings['N'] = 5  # N**2 or N**3 are the number of PSCs
settings['diameter'] = 1  # Diameter of halo particles
settings['sigma'] = 1.0  # WCA-potential parameters
settings['epsilon'] = 1.0  # WCA-potential parameters
settings['mass'] = 1.0  # Mass of halo particles
settings['nameString'] = 'integrator-{integrator}_shape-{poly}_N-{N}_VF-{density:4.2f}_dim-{dimensions}_Nclus-{N_cluster}_tstep-{time_step}'
settings["initFile"] = 'None'
settings['outputInterval'] = 20  # Number of time steps between data storage
settings['therm_steps'] = 400  # Number of thermalization steps
settings['equil_steps'] = 700  # Number of equilibration steps


nameFormat = "dataFLEX_{poly}/" + settings['nameString']


# Specific simulation parameters
parameterspace = []

dens = int(args.density)/100
tstep_multiplier = 0.003

start_N_cluster = int(args.start_N_cluster)
end_N_cluster = int(args.end_N_cluster)

'''
for dens in dens_values:
    for tstep_multiplier in tstep_values:
        parameters = [
            {**settings,
             'integrator': 'nve',
             'density': dens,
             'poly': 'octa',
             'dimensions': 3,
             'N_cluster': 6,
             'time_step': tstep_multiplier*math.sqrt(settings['mass']*settings['sigma']**2/settings['epsilon']),
             # 'initFile': [nameFormat.format(**settings)+'_restart-000.gsd']
             },
            {**settings,
             'integrator': 'nve',
             'density': dens,
             'poly': 'tetra',
             'dimensions': 3,
             'N_cluster': 4,
             'time_step': tstep_multiplier*math.sqrt(settings['mass']*settings['sigma']**2/settings['epsilon']),
             # 'initFile': [nameFormat.format(**settings)+'_restart-000.gsd']
             },
            {**settings,
             'integrator': 'nve',
             'density': dens,
             'poly': 'dode',
             'dimensions': 3,
             'N_cluster': 20,
             'time_step': tstep_multiplier*math.sqrt(settings['mass']*settings['sigma']**2/settings['epsilon']),
             # 'initFile': [nameFormat.format(**settings)+'_restart-000.gsd']
             },
            {**settings,
             'integrator': 'nve',
             'density': dens,
             'poly': 'ico',
             'dimensions': 3,
             'N_cluster': 12,
             'time_step': tstep_multiplier*math.sqrt(settings['mass']*settings['sigma']**2/settings['epsilon']),
             # 'initFile': [nameFormat.format(**settings)+'_restart-000.gsd']
             },
            {**settings,
             'integrator': 'nve',
             'density': dens,
             'poly': 'cube',
             'dimensions': 3,
             'N_cluster': 8,
             'time_step': tstep_multiplier*math.sqrt(settings['mass']*settings['sigma']**2/settings['epsilon']),
             # 'initFile': [nameFormat.format(**settings)+'_restart-000.gsd']
             }
        ]

        parameterspace += parameters

for dens in dens_values:
    for tstep_multiplier in tstep_values:
        for i in range(3, 50, 1):
            parameterspace += [
                {**settings,
                 'integrator': 'nve',
                 'density': dens,
                 'poly': '3Dspheres',
                 'dimensions': 3,
                 'N_cluster': i,
                 'time_step': tstep_multiplier*math.sqrt(settings['mass']*settings['sigma']**2/settings['epsilon']),
                 # 'initFile': [nameFormat.format(**settings)+'_restart-000.gsd']
                 }]
'''

for i in range(start_N_cluster, end_N_cluster, 1):
    parameterspace += [
        {**settings,
            'integrator': 'nve',
            'density': dens,
            'poly': '2Dspheres',
            'dimensions': 2,
            'N_cluster': i,
            'time_step': tstep_multiplier*math.sqrt(settings['mass']*settings['sigma']**2/settings['epsilon']),
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
        if(not os.path.exists("dataFLeX_{poly}/".format(**initDict))):
            os.mkdir("dataFLEX_{poly}/".format(**initDict))
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
                             "GnanNVE.py",
                             *initString],
                            stdout=out,
                            stderr=out)
    proc.wait()
