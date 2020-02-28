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
from Marson22NVE import Marson2NVE

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
settings['N'] = [5]  # Number of PSCs
settings['diameter'] = [1]  # Diameter of halo particles
settings['poly'] = ['octa']  # Type of polyhedron
settings['sigma'] = [1.0]  # WCA-potential parameters
settings['epsilon'] = [1.0]  # WCA-potential parameters
settings['mass'] = [1.0]  # Mass of halo particles
settings['density'] = [0.3]  # Volume fraction
settings['dimensions'] = [3]  # 2d or 3d
settings['integrator'] = ['local']  # Integrator
settings['nameString'] = [
    'integrator_{integrator}_shape_{poly}_N_{N}_VF_{density:4.2f}_dim_{dimensions}_Nclus_{N_cluster}_tstep_{time_step}']
settings["initFile"] = ['None']
settings['outputInterval'] = [100]  # Number of time steps between data storage
settings['N_cluster'] = 17  # Number of halo spheres in cluster

nameFormat = "data_{poly}/" + settings['nameString'][0]

# Specific simulation parameters

parameterspace = [
    {**settings,
     'integrator': ['nve'],
     'density': [0.7],
     'poly': ['octa'],
     'dimensions': [3],
     'N': [5],
     'N_cluster': [6],
     'time_step': [0.007*math.sqrt(settings['mass'][0]*settings['sigma'][0]**2/settings['epsilon'][0])],
     # 'initFile': [nameFormat.format(**settings)+'_restart-000.gsd']
     },
    {**settings,
     'integrator': ['nve'],
     'density': [0.7],
     'poly': ['tetra'],
     'dimensions': [3],
     'N': [5],
     'N_cluster': [4],
     'time_step': [0.007*math.sqrt(settings['mass'][0]*settings['sigma'][0]**2/settings['epsilon'][0])],
     # 'initFile': [nameFormat.format(**settings)+'_restart-000.gsd']
     },
]

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
    Marson2NVE(**initDict)

    sys.stdout = normalstdout
    sys.stderr = normalstderr
