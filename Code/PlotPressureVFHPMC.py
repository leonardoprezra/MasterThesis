'''
Plots P-VF from given _sdf.log files in current directory.
Simulations must follow HPMC.

!!!!!!!!!!!!!
dx and xmax must be read from hoomd.hpmc.analyze.sdf() from MarsonNVEHisteresisOneHPMC.py
dx = 1e-4
xmax=0.02

'''

from matplotlib import pyplot as plt
import numpy as np
import os
import math
from MarsonFunctions import PartCluster
import argparse
import sys
import scipy.constants

# Command line argument parsing
parser = argparse.ArgumentParser(
    description='Plots Volume Fraction vs Pressure from .log file, within the specified Volume Fractions.')
parser.add_argument('-a', '--frame-average', type=int, dest='frame_ave',
                    help='number of frames to average')
parser.add_argument('-j', '--frame-jump', type=int, dest='frame_jump',
                    help='number of frames at the same Volume Fraction')
parser.add_argument('-t', '--frame-total', type=int, dest='frame_total',
                    help='tota number of frames in the .log file')
parser.add_argument('-f', '--in-file', type=str,
                    dest='in_file', help='input .log file')
parser.add_argument('-s', '--sdf-file', type=str,
                    dest='sdf_file', help='input _sdf.log file')

args = parser.parse_args()


try:
    in_file = args.in_file
    d = args.in_file.split('/')[-1]
except:
    in_file = args.in_file
    d = args.in_file

# Name of directory to store converted files
dir_name = 'thermo_plots_hpmc/'

# Create directories
try:
    if(not os.path.exists(dir_name)):
        os.mkdir(dir_name)
except OSError as e:
    if e.errno != 17:
        raise
    pass

# #
# #
# #
# #
# Read and store data from files
# Read _sdf.log data
try:
    file_sdf = args.sdf_file
    d_sdf = args.sdf_file.split('/')[-1]
except:
    file_sdf = args.sdf_file
    d_sdf = args.sdf_file

print(file_sdf)
print(d_sdf)

# Read .log data
try:
    file_log = args.in_file
    d_log = args.in_file.split('/')[-1]
except:
    file_log = args.in_file
    d_log = args.in_file

print(file_log)
print(d_log)

# # Information from file name
# Label includes: (5) Nclus, (4) dim, (3) VF, (7) ratio
label = d_sdf.split('_')[4] + '_' + \
    d_sdf.split('_')[5] + '_' + d_sdf.split('_')[7]
# Number of clusters in simulation (2)
total_N = int(d_sdf.split('_')[2].split('-')[1])
# Number of particles per cluster (4)
N_cluster = int(d_sdf.split('_')[5].split('-')[1])
# Shape of clusters (1)
poly_key = d_sdf.split('_')[1].split('-')[1]
# Check dimensions
dimensions = int(d_sdf.split('_')[4].split('-')[1])

# Volume of particles in simulation
if poly_key == 'one':
    if dimensions == 3:
        vol = math.pi/6*1**3 * total_N
    if dimensions == 2:
        vol = math.pi/4*1**2 * total_N
else:
    sys.exit(
        "The file passed do not correspond to the results of HPMC simulation")

# Store _sdf.log data
data_sdf = np.genfromtxt(fname=file_sdf, skip_header=False)
# excludes data of last VF, change required because vol_log doesn't register the last VF
'''
# Use for contents of out700
data_sdf = data_sdf[:, 1:]

print(data_sdf.shape)
'''
data_sdf = data_sdf[:, 1:]
data_sdf = [np.mean(data_sdf[i-4:i+1, :], axis=0)
            for i in range((10-1), (args.frame_total+1)*10 - 1, 10)]

data_sdf = np.array(data_sdf)


# Store .log data
data_log = np.genfromtxt(fname=file_log, skip_header=True)

# Get VF from data from .log file
vol_log = np.array([[data_log[i-1, 0], vol / np.mean(data_log[i-args.frame_ave:i, 1])] for i in
                    range(args.frame_jump, (args.frame_total+1)*args.frame_jump, args.frame_jump)])

num_dens = np.array([[np.mean(data_log[i-args.frame_ave:i, 1])] for i in
                     range(args.frame_jump, (args.frame_total+1)*args.frame_jump, args.frame_jump)])
num_dens = total_N / num_dens

# print(vol_log[-1:])

# #
# #
# #
# #
# Calculate pressure
# Extrapolation of scale density function
# dx and xmax must be read from hoomd.hpmc.analyze.sdf() from MarsonNVEHisteresisOneHPMC.py
dx = 1e-5
xmax = 0.002


def extrapolate(s, dx, xmax, degree=5):
    # determine the number of values to fit
    n_fit = 0
    for i in range(int(math.ceil(xmax/dx))):

        if np.sum(s[0: i]*dx) > 0.5:
            n_fit = i
            break

    if n_fit == 0:
        n_fit = int(math.ceil(xmax/dx))  # uses the whole range of values

    s_fit = s[0: n_fit]

    # construct the x coordinates
    # x_fit = np.arange(0, xmax, dx) # uses the whole range of values
    x_fit = np.arange(0, dx*n_fit, dx)
    x_fit += dx/2
    # perform the fit and extrapolation
    p = np.polyfit(x_fit, s_fit, degree)
    '''
    pp = np.poly1d(p)

    # plot data and fit
    plt.plot(x_fit, data_sdf[0], label='sdf')
    plt.plot(x_fit, pp(x_fit), label='polyfit')
    plt.legend()
    plt.show()
    '''
    return np.polyval(p, 0.0)


# Values of scale density function at s(0+)
p_sdf = np.array([[extrapolate(i, dx=dx, xmax=xmax)] for i in data_sdf])

# P/kT * Area_single_particle
p_sdf = num_dens*(1+p_sdf/(2*dimensions)) * math.pi/4*1**2

# data = [timestep, VF, P/kT]
data_ave = np.append(vol_log, p_sdf, axis=1)

# #
# #
# #
# # Plot Pressure-VF
fig5 = plt.figure(5)
ax5_1 = fig5.add_subplot(1, 1, 1)
ax5_1.set_title('data-{}\nN-{}'.format(poly_key, total_N))
plt.rcParams.update({'mathtext.default': 'regular'})

mid_point = int(args.frame_total/2)
# Plot compresion
ax5_1.plot(data_ave[:mid_point, 1], data_ave[:mid_point, 2],
           label=label+'_compr', linewidth=0.5)

# Plot expansion
ax5_1.plot(data_ave[mid_point:, 1], data_ave[mid_point:, 2],
           label=label+'_exp', linewidth=0.5)


# Axis labels
ax5_1.set_ylabel('$\dfrac{PA_1}{kT}$ / -')
ax5_1.set_xlabel('$\phi$ / -')
ax5_1.legend()
fig5.tight_layout()
fig5.savefig(
    dir_name + 'PressureVF_data-{}_N-{}_Nclus-{}.pdf'.format(poly_key, total_N, N_cluster), format='pdf')
plt.clf()


# Array contains [VF,PRESSURE,STD]
#save_data = data_ave[:,1:]
save_data = np.concatenate((data_ave[:, 1:], np.zeros((182, 1))), axis=1)
print(save_data.shape)
# Saves data in .npy file
np.save(dir_name + d[:-4] + '_Pressure.npy', save_data)
