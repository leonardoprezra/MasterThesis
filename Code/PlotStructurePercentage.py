'''
Plots variables from log files in current directory.

Run in a directory that contain the data_* directories,
no sub directories should exist in data_*.


'''

from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import numpy as np
import os
import math
from MarsonFunctions import PartCluster


# Name of directory to store converted files
dir_name = 'StructurePercentage/'

# Create directories
try:
    if(not os.path.exists(dir_name)):
        os.mkdir(dir_name)
except OSError as e:
    if e.errno != 17:
        raise
    pass

data_ico = np.array([[000, 4.1, 94.9, 1, 0],
                     [300, 7.7, 83.8, 8.5, 0],
                     [400, 4.7, 94, 1.3, 0],
                     [500, 5.6, 75.8, 17.3, 1.3],
                     [600, 4.1, 94.4, 1.5, 0],
                     [700, 6.6, 89.4, 4, 0],
                     [800, 15.4, 58, 26.5, 0],
                     [900, 9.6, 77.3, 13, 0],
                     [1000, 9, 79.7, 11.2, 0]])

data_3DSoft = np.array([[000, 6.6, 45.5, 45.8, 2],
                        [100, 100, 0, 0, 0],
                        [200, 19.3, 59.8, 20.8, 0],
                        [300, 19.8, 56.4, 23.7, 0],
                        [400, 26.9, 45.5, 27.5, 0],
                        [500, 11.7, 71.7, 16.6, 0],
                        [600, 20.7, 46.9, 32.4, 0],
                        [700, 16.9, 52.2, 30.8, 0.1],
                        [800, 17.9, 34.5, 47.3, 0.2],
                        [900, 20.4, 54.9, 24.6, 0],
                        [1000, 24.7, 42, 33.3, 0]])

xlab = ['Other', 'FCC', 'HCP', 'BCC']

# #
# #
# #
# # Plot DATASOFT3D
fig5 = plt.figure(5)
ax5_1 = fig5.add_subplot(1, 1, 1)
ax5_1.plot(xlab, data_3DSoft[0, 1:], label='hard_particle')
ax5_1.set_ylabel('$Fraction$ / %')
ax5_1.set_xlabel('$Structure Type$')

e_t = int(data_3DSoft.shape[0]/2)

for i in range(1, e_t):
    ax5_1.plot(xlab, data_3DSoft[i, 1:],
               label='harmk-{}'.format(int(data_3DSoft[i, 0])))

ax5_1.legend()
fig5.savefig(dir_name + 'Structure_DATASOFT3D_1.pdf', format='pdf')
fig5.clear()

# #
# #
# #
# # Plot DATASOFT3D
fig5 = plt.figure(5)
ax5_1 = fig5.add_subplot(1, 1, 1)
ax5_1.plot(xlab, data_3DSoft[0, 1:], label='hard_particle')
ax5_1.set_ylabel('$Fraction$ / %')
ax5_1.set_xlabel('$Structure Type$')


for i in range(e_t, data_3DSoft.shape[0]):
    ax5_1.plot(xlab, data_3DSoft[i, 1:],
               label='harmk-{}'.format(int(data_3DSoft[i, 0])))
ax5_1.legend()
fig5.savefig(dir_name + 'Structure_DATASOFT3D_2.pdf', format='pdf')
fig5.clear()

# #
# #
# #
# # Plot ICO
fig5 = plt.figure(5)
ax5_1 = fig5.add_subplot(1, 1, 1)
ax5_1.plot(xlab, data_ico[0, 1:], label='ico_rigid')
ax5_1.set_ylabel('$Fraction$ / %')
ax5_1.set_xlabel('$Structure Type$')

e_t = int(data_ico.shape[0]/2)

for i in range(1, e_t):
    ax5_1.plot(xlab, data_ico[i, 1:],
               label='harmk-{}'.format(int(data_ico[i, 0])))

ax5_1.legend()
fig5.savefig(dir_name + 'Structure_ICO_1.pdf', format='pdf')
fig5.clear()

# #
# #
# #
# # Plot ICO
fig5 = plt.figure(5)
ax5_1 = fig5.add_subplot(1, 1, 1)
ax5_1.plot(xlab, data_ico[0, 1:], label='ico_rigid')
ax5_1.set_ylabel('$Fraction$ / %')
ax5_1.set_xlabel('$Structure Type$')

for i in range(e_t, data_ico.shape[0]):
    ax5_1.plot(xlab, data_ico[i, 1:],
               label='harmk-{}'.format(int(data_ico[i, 0])))

ax5_1.legend()
fig5.savefig(dir_name + 'Structure_ICO_2.pdf', format='pdf')
fig5.clear()
