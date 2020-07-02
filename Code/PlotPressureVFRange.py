'''
Plots variables from log files in current directory.

Run in a directory that contain the data_* directories,
no sub directories should exist in data_*.


'''

from matplotlib import pyplot as plt
import numpy as np
import os
import math
from MarsonFunctions import PartCluster
import argparse

# Command line argument parsing
parser = argparse.ArgumentParser(
    description='Plots Volume Fraction vs Pressure from .log file, within the specified Volume Fractions.')
parser.add_argument('-a', '--frame-average', type=int, dest='frame_ave',
                    help='number of frames to average')
parser.add_argument('-j', '--frame-jump', type=int, dest='frame_jump',
                    help='number of frames at the same Volume Fraction')
parser.add_argument('-t', '--frame-total', type=int, dest='frame_total',
                    help='tota number of frames in the .log file')
parser.add_argument('-i', '--initial-density', type=float, dest='init_dens',
                    help='initial density of the range to plot')
parser.add_argument('-e', '--end-density', type=float, dest='end_dens',
                    help='end density of the range to plot')
parser.add_argument('-f', '--in-file', type=str, nargs='+',
                    dest='in_file', help='input .gsd file')

args = parser.parse_args()

# Name of directory to store converted files
dir_name = 'thermo_plots_range/'

# Create directories
try:
    if(not os.path.exists(dir_name)):
        os.mkdir(dir_name)
except OSError as e:
    if e.errno != 17:
        raise
    pass


# Read and plot data
for d in args.in_file:

    try:
        in_file = d
        d = d.split('/')[-1]
    except:
        in_file = d
        d = d

    print(in_file)
    print(d)

    # # Information from file name
    # Label includes: (5) Nclus, (4) dim, (3) VF, (7) ratio
    label = d.split('_')[4] + '_' + \
        d.split('_')[5] + '_' + d.split('_')[7]
    # Number of clusters in simulation (2)
    total_N = int(d.split('_')[2].split('-')[1])
    # Number of particles per cluster (4)
    N_cluster = int(d.split('_')[5].split('-')[1])
    # Shape of clusters (1)
    poly_key = d.split('_')[1].split('-')[1]
    # Check dimensions
    dimensions = int(d.split('_')[4].split('-')[1])

    # Volume of particles in simulation
    if poly_key == 'one':
        if dimensions == 3:
            vol = math.pi/6*1**3 * total_N
        if dimensions == 2:
            vol = math.pi/4*1**2 * total_N
    else:
        ratio = float(d.split('_')[7].split('-')[1])

        cluster = PartCluster(
            poly_key=poly_key, N_cluster=N_cluster, halo_diam=1, halo_mass=1, ratio=ratio)

        vol = cluster.vol_cluster(dimensions) * total_N

    # Read data
    data = np.genfromtxt(fname=in_file, skip_header=True)

    # # Plot Pressure-VF
    fig5 = plt.figure(5)
    ax5_1 = fig5.add_subplot(1, 1, 1)
    ax5_1.set_title('data-{}\nN-{}'.format(poly_key, total_N))

    # Averages the Pressure data points at each VF
    data_ave = [[vol / np.mean(
        data[i-args.frame_ave:i, 1]), np.mean(data[i-args.frame_ave:i, 9])] for i in
        range(args.frame_jump, args.frame_total*args.frame_jump, args.frame_jump)]
    data_ave = np.array(data_ave)

    # initial density frame of compresion section
    i_d = int((args.init_dens-data_ave[0, 0])*args.frame_jump)
    # end density frame of compresion section
    e_d = int((args.end_dens-data_ave[0, 0])*args.frame_jump)
    # initial density frame of expansion section
    i_d_r = -e_d - 1
    # end density frame of expansion section
    e_d_r = -i_d
    # initial density, correct frame
    i_d = int((args.init_dens-data_ave[0, 0])*args.frame_jump) + 1
    # end density, correct frame
    e_d = int((args.end_dens-data_ave[0, 0])*args.frame_jump) + 2

    # Error bars
    error = [np.std(data[i-args.frame_ave:i, 9]) for i in
             range(args.frame_jump, args.frame_total*args.frame_jump, args.frame_jump)]
    error = np.array([error]).T

    # Plot compresion
    ax5_1.errorbar(data_ave[i_d:e_d, 0], data_ave[i_d:e_d, 1],
                   yerr=error[i_d:e_d, 0], label=label+'_compr', linewidth=0.5)

    # Plot expansion
    ax5_1.errorbar(data_ave[i_d_r:e_d_r, 0], data_ave[i_d_r:e_d_r, 1],
                   yerr=error[i_d_r:e_d_r, 0], label=label+'_exp', linewidth=0.5)

    print('!!!!!!!!!!!PRESSURE!!!!!!!!!!!\n' +
          str(d + '\n' + str(data_ave.shape)))

    # Secondary axis, frames in gsd file
    secax = ax5_1.secondary_xaxis('top', functions=(lambda x: (
        x-data_ave[0, 0])*args.frame_jump, lambda x: x/args.frame_jump + data_ave[0, 0]))

    # Axis labels
    secax.set_xlabel('Frame number')
    ax5_1.set_ylabel('Pressure / -')
    ax5_1.set_xlabel('$\phi$ / -')
    # ax4_1.xaxis.set_minor_locator(plt.MultipleLocator(500))
    ax5_1.legend()
    fig5.tight_layout()
    fig5.savefig(
        dir_name + 'PressureVF_data-{}_N-{}_Nclus-{}_IDens-{}_EDens-{}.pdf'.format(poly_key, total_N, N_cluster, args.init_dens, args.end_dens), format='pdf')
    plt.clf()

    # Prepare array with only information in the range
    data_compr = data_ave[i_d:e_d, :]
    data_exp = data_ave[i_d_r:e_d_r, :]
    data_compl = np.concatenate((data_compr, data_exp), axis=0)
    error_compr = error[i_d:e_d]
    error_exp = error[i_d_r:e_d_r]
    err_compl = np.concatenate((error_compr, error_exp), axis=0)

    # Array contains [VF,PRESSURE,STD]
    save_data = np.concatenate((data_compl, err_compl), axis=1)

    # Saves data in .npy file
    np.save(dir_name + d[:-4] + '_Pressure.npy', save_data)

    # # Plot Energy-VF
    fig1 = plt.figure(1)
    ax1_1 = fig1.add_subplot(1, 1, 1)
    ax1_1.set_title('data-{}\nN-{}'.format(poly_key, total_N))

    # Averages the Energy data points at each VF
    data_ave = [[vol / np.mean(
        data[i-args.frame_ave:i, 1]), np.mean(data[i-args.frame_ave:i, 4] + data[i-args.frame_ave:i, 5])] for i in
        range(args.frame_jump, args.frame_total*args.frame_jump, args.frame_jump)]
    data_ave = np.array(data_ave)

    # Error bars
    error = [np.std(data[i-args.frame_ave:i, 4] + data[i-args.frame_ave:i, 5]) for i in
             range(args.frame_jump, args.frame_total*args.frame_jump, args.frame_jump)]
    error = np.array([error]).T

    # Plot compresion
    ax5_1.errorbar(
        data_ave[i_d:e_d, 0], data_ave[i_d:e_d, 1], yerr=error[i_d:e_d, 0], label=label+'_compr', linewidth=0.5)

    # Plot expansion
    ax5_1.errorbar(
        data_ave[i_d_r:e_d_r, 0], data_ave[i_d_r:e_d_r, 1], yerr=error[i_d_r:e_d_r, 0], label=label+'_exp', linewidth=0.5)

    print('!!!!!!!!!!!!ENERGY!!!!!!!!!!!!\n' +
          str(d + '\n' + str(data_ave.shape)))

    # Secondary axis, frames in gsd file
    secax = ax1_1.secondary_xaxis('top', functions=(lambda x: (
        x-data_ave[0, 0])*args.frame_jump, lambda x: x/args.frame_jump + data_ave[0, 0]))

    # Axis labels
    secax.set_xlabel('Frame number')
    ax1_1.set_ylabel('Energy / -')
    ax1_1.set_xlabel('$\phi$ / -')
    # ax4_1.xaxis.set_minor_locator(plt.MultipleLocator(500))
    ax1_1.legend()
    fig1.tight_layout()
    fig1.savefig(
        dir_name + 'EnergyVF_data-{}_N-{}_Nclus-{}_IDens-{}_EDens-{}.pdf'.format(poly_key, total_N, N_cluster, args.init_dens, args.end_dens), format='pdf')
    plt.clf()

    # Prepare array with only information in the range
    data_compr = data_ave[i_d:e_d, :]
    data_exp = data_ave[i_d_r:e_d_r, :]
    data_compl = np.concatenate((data_compr, data_exp), axis=0)
    error_compr = error[i_d:e_d]
    error_exp = error[i_d_r:e_d_r]
    err_compl = np.concatenate((error_compr, error_exp), axis=0)

    # Array contains [VF,ENERGY,STD]
    save_data = np.concatenate((data_compl, err_compl), axis=1)

    # Saves data in .npy file
    np.save(dir_name + d[:-4] + '_Energy.npy', save_data)
