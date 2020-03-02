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
from MarsonNVE import MarsonNVE

# General simulation parameters
settings = {}
settings['N'] = 10  # N**2 or N**3 are the number of PSCs
settings['diameter'] = 1  # Diameter of halo particles
settings['sigma'] = 1.0  # WCA-potential parameters
settings['epsilon'] = 1.0  # WCA-potential parameters
settings['mass'] = 1.0  # Mass of halo particles
settings['nameString'] = 'integrator-{integrator}_shape-{poly}_N-{N}_VF-{density:4.2f}_dim-{dimensions}_Nclus-{N_cluster}_tstep-{time_step}'
settings["initFile"] = 'None'
settings['outputInterval'] = 100  # Number of time steps between data storage

nameFormat = "data_{poly}/" + settings['nameString']


# Specific simulation parameters
parameterspace = []
dens_values = [0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8]
tstep_values = [0.001, 0.002, 0.003, 0.004, 0.005,
                0.006, 0.007, 0.008, 0.009, 0.01,
                0.02, 0.03, 0.04, 0.05, 0.06,
                0.07, 0.08, 0.09, 0.1]


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


# Run Simulations
for initDict in parameterspace:

    if initDict['dimensions'] == 3:
        initDict['N'] = initDict['N']**3
    elif initDict['dimensions'] == 2:
        initDict['N'] = initDict['N']**2

    nameString = nameFormat.format(**initDict)

    if(not os.path.exists("data_{poly}/".format(**initDict))):
        os.mkdir("data_{poly}/".format(**initDict))

    if(os.path.exists(nameString+".outputs")):
        print("Skipping "+nameString)
        continue
    else:
        print("Lauing "+nameString)

    normalstdout = sys.stdout
    normalstderr = sys.stderr
    '''
    out = open(nameString+".outputs", "w")
    sys.stderr = out
    sys.stdout = out
    '''
    try:
        MarsonNVE(**initDict)
    except:
        sys.stdout = normalstdout
        sys.stderr = normalstderr
        continue

    #sys.stdout = normalstdout
    #sys.stderr = normalstderr
