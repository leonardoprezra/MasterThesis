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

paraFormat = "{:s}={:s}"

# Create list to pass command line arguments to python scripts


def listOfParameterInits(inDict):
    myList = [[[], {}]]

    for parameter, valueList in inDict.items():
        newList = []
        for initEntry in myList:
            for newValue in valueList:
                newPara = paraFormat.format(str(parameter), str(newValue))
                newDict = {**initEntry[1], parameter: newValue}
                newList.append([[*initEntry[0], newPara], newDict])
        myList = newList
    return myList


# General simulation parameters
settings = {}
settings['N'] = [10]  # Number of PSCs
settings['diameter'] = [1]  # Diameter of halo particles
settings['poly'] = ['octa']  # Type of polyhedron
settings['sigma'] = [1.0]  # WCA-potential parameters
settings['epsilon'] = [1.0]  # WCA-potential parameters
settings['mass'] = [1.0]  # Mass of halo particles
settings['density'] = [0.3]  # Volume fraction
settings['dimensions'] = [3]  # 2d or 3d
settings['integrator'] = ['local']  # Integrator
settings['nameString'] = [
    'integrator-{integrator}_shape-{poly}_N-{N}_VF-{density:4.2f}_dim-{dimensions}_Nclus-{N_cluster}_tstep-{time_step}']
settings["initFile"] = ['None']
settings['outputInterval'] = [100]  # Number of time steps between data storage
settings['N_cluster'] = 17  # Number of halo spheres in cluster

nameFormat = "data_{poly}/" + settings['nameString'][0]



# Specific simulation parameters
parameterspace = []
dens_values = [0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8]
tstep_values = [0.001,0.002,0.003,0.004,0.005,
                0.006,0.007,0.008,0.009,0.01,
                0.02,0.03,0.04,0.05,0.06,
                0.07,0.08,0.09,0.1]


for dens in dens_values:
    for tstep_multiplier in tstep_values:
        parameters = [
            {**settings,
            'integrator': ['nve'],
            'density': [dens],
            'poly': ['octa'],
            'dimensions': [3],
            'N_cluster': [6],
            'time_step': [tstep_multiplier*math.sqrt(settings['mass'][0]*settings['sigma'][0]**2/settings['epsilon'][0])],
            # 'initFile': [nameFormat.format(**settings)+'_restart-000.gsd']
            },
            {**settings,
            'integrator': ['nve'],
            'density': [dens],
            'poly': ['tetra'],
            'dimensions': [3],
            'N_cluster': [4],
            'time_step': [tstep_multiplier*math.sqrt(settings['mass'][0]*settings['sigma'][0]**2/settings['epsilon'][0])],
            # 'initFile': [nameFormat.format(**settings)+'_restart-000.gsd']
            },
            {**settings,
            'integrator': ['nve'],
            'density': [dens],
            'poly': ['dode'],
            'dimensions': [3],
            'N_cluster': [20],
            'time_step': [tstep_multiplier*math.sqrt(settings['mass'][0]*settings['sigma'][0]**2/settings['epsilon'][0])],
            # 'initFile': [nameFormat.format(**settings)+'_restart-000.gsd']
            },
            {**settings,
            'integrator': ['nve'],
            'density': [dens],
            'poly': ['ico'],
            'dimensions': [3],
            'N_cluster': [12],
            'time_step': [tstep_multiplier*math.sqrt(settings['mass'][0]*settings['sigma'][0]**2/settings['epsilon'][0])],
            # 'initFile': [nameFormat.format(**settings)+'_restart-000.gsd']
            },
            {**settings,
            'integrator': ['nve'],
            'density': [dens],
            'poly': ['cube'],
            'dimensions': [3],
            'N_cluster': [8],
            'time_step': [tstep_multiplier*math.sqrt(settings['mass'][0]*settings['sigma'][0]**2/settings['epsilon'][0])],
            # 'initFile': [nameFormat.format(**settings)+'_restart-000.gsd']
            }
        ]

        parameterspace += parameters

for dens in dens_values:
    for tstep_multiplier in tstep_values:
        for i in range(3,50,1):    
            parameterspace += [
                {**settings,
                'integrator': ['3Dspheres'],
                'density': [dens],
                'poly': ['cube'],
                'dimensions': [3],
                'N_cluster': [i],
                'time_step': [tstep_multiplier*math.sqrt(settings['mass'][0]*settings['sigma'][0]**2/settings['epsilon'][0])],
                # 'initFile': [nameFormat.format(**settings)+'_restart-000.gsd']
                }]
            

# Run Simulations
for initString, initDict in [l for dictionary in parameterspace
                            for l in listOfParameterInits(dictionary)]:

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

    sys.stdout = open(nameString+".outputs", "w")
    sys.stderr = open(nameString+".err", "w")
    MarsonNVE(**initDict)

    sys.stdout = normalstdout
    sys.stderr = normalstderr
