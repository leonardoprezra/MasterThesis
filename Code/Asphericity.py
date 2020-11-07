'''
Extract the last <frame_num> of the <in_file> and saves them as .gsd and .pos.

'''
import sys
import argparse
import os

from MarsonFunctions import PartCluster

import gsd
import gsd.hoomd

import math
import numpy as np
from numpy import linalg as LA

from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

# Name of directory to store converted files
dir_name = 'asphericity/'

# Command line argument parsing
parser = argparse.ArgumentParser(
    description='Extract the last <frame_num> of the <in_file> and saves them as .gsd and .pos.')
parser.add_argument('-f', '--in-file', type=str, nargs='+',
                    dest='in_file', help='input .gsd file')
parser.add_argument('-j', '--frame-jump', type=int, dest='frame_jump',
                    help='number of frames at the same Volume Fraction')
parser.add_argument('-i', '--initial-density', type=float, dest='init_dens',
                    help='initial density of the range to plot')
parser.add_argument('-e', '--end-density', type=float, dest='end_dens',
                    help='end density of the range to plot')

args = parser.parse_args()

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
# Calculate distribution of radius of gyration
def gyr_r(snapshots, vol_tot, all_index):

    N = len(all_index[0])
    # print('size= {}'.format(N))
    # print('N={}'.format(N))

    b_ave = []
    dens_ave = []

    for snap in snapshots:
        # Initial density and frames of interest
        box_dim = snap.configuration.box[:3]

        if dimensions == 3:
            box_vol = box_dim[0]*box_dim[1]*box_dim[2]
        if dimensions == 2:
            box_vol = box_dim[0]*box_dim[1]

        dens = vol / box_vol

        r_cm_all = []
        b_all = []

        # Coordinates of center of mass
        for cluster_index in all_index:
            r_cm_prel = np.zeros(3)
            for indiv_index in cluster_index:
                r_cm_prel += snap.particles.position[indiv_index]
            r_cm_all.append(r_cm_prel/N)

        # S : Gyration tensor
        for i, cluster_index in enumerate(all_index):
            Sxx = 0.0
            Syy = 0.0
            Szz = 0.0
            Sxy = 0.0
            Syz = 0.0
            Sxz = 0.0

            # r_cm : center of mass
            r_cm = r_cm_all[i]

            # print('r_cm={}'.format(r_cm))

            for indiv_index in cluster_index:
                coord = snap.particles.position[indiv_index]

                Sxx += (coord[0] - r_cm[0])**2
                Syy += (coord[1] - r_cm[1])**2
                Szz += (coord[2] - r_cm[2])**2
                Sxy += (coord[0] - r_cm[0])*(coord[1] - r_cm[1])
                Syz += (coord[1] - r_cm[1])*(coord[2] - r_cm[2])
                Sxz += (coord[0] - r_cm[0])*(coord[2] - r_cm[2])

            S = np.array([[Sxx, Sxy, Sxz],
                          [Sxy, Syy, Syz],
                          [Sxz, Syz, Szz]])

            S = S/N

            eigen_values, princ_axis = LA.eig(S)

            lam = sorted(eigen_values)
            # Rg : Radius of gyration
            Rg = np.sqrt(np.sum(np.square(lam)))

            # b : asphericity
            b = ((lam[1] - lam[0])**2 + (lam[2] - lam[1]) ** 2) / \
                (2*(lam[0] + lam[1] + lam[2])**2)

            b_all.append(b)

            print('lam={}'.format(lam))
            print('b={}'.format(b))
            print('Rg={}'.format(Rg))
            print('dens={}'.format(dens)

        b_all=np.array(b_all)
        '''
        print('asphericity= {}'.format(b_all))
        print('size= {}'.format(len(b_all)))
        exit()
        '''
        b_ave.append(np.mean(b_all))
        dens_ave.append(dens)

        plt.hist(b_all, bins=1000)

        exit()

    return [dens_ave, b_ave]


# #
# #
# #
# #
# Loop over files
fig1=plt.figure(1)
ax1_1=fig1.add_subplot(1, 1, 1)


for d in args.in_file:
    # #
    # #
    # #
    # #
    # Store label data
    try:
        in_file=d
        d=d.split('/')[-1]
    except:
        in_file=d
        d=d

    print(in_file)
    print('Working on : \n' + d)

    # # Information from file name
    # Label includes: (5) Nclus, (4) dim, (3) VF, (7) ratio
    label='d-' + d.split('_')[4].split('-')[1] + '_' + \
        'n-' + d.split('_')[5].split('-')[1] + '_'  # + \
    # 'r-' + ''d.split('_')[7].split('-')[1]

    # Check for harmk
    if d.split('_')[-1].split('-')[0] == 'harmk':
        label='d-' + d.split('_')[4].split('-')[1] + '_' + \
            'n-' + d.split('_')[5].split('-')[1] + '_' + \
            'k-' + d.split('_')[-1].split('-')[1][:-4]

        harmk=d.split('_')[-1].split('-')[1][:-4]

        if int(d.split('_')[-1].split('-')[1][:-4]) == 0:
            label='d-' + d.split('_')[4].split('-')[1] + '_' + \
                'n-' + d.split('_')[5].split('-')[1] + '_' + \
                'k-' + 'one'

            harmk='once'

    # Number of clusters in simulation (2)
    total_N=int(d.split('_')[2].split('-')[1])
    # Number of particles per cluster (4)
    N_cluster=int(d.split('_')[5].split('-')[1])
    # Shape of clusters (1)
    poly_key=d.split('_')[1].split('-')[1]
    # Check dimensions
    dimensions=int(d.split('_')[4].split('-')[1])

    # Volume of particles in simulation
    if poly_key == 'one':
        if dimensions == 3:
            vol=math.pi/6*1**3 * total_N
        if dimensions == 2:
            vol=math.pi/4*1**2 * total_N
    else:
        ratio=float(d.split('_')[7].split('-')[1])

        cluster=PartCluster(
            poly_key=poly_key, N_cluster=N_cluster, halo_diam=1, halo_mass=1, ratio=ratio, dimensions=dimensions)

        vol=cluster.vol_cluster(dimensions) * total_N

    # #
    # #
    # #
    # #
    # Read initial snapshot
    t=gsd.hoomd.open(name=in_file, mode='rb')

    # Initial density and frames of interest
    box_dim=t[0].configuration.box[:3]

    if dimensions == 3:
        box_vol=box_dim[0]*box_dim[1]*box_dim[2]
    if dimensions == 2:
        box_vol=box_dim[0]*box_dim[1]

    dens=vol / box_vol

    # initial density frame of compresion section
    i_d=int((args.init_dens-dens)*args.frame_jump)
    # end density frame of compresion section
    e_d=int((args.end_dens-dens)*args.frame_jump)
    # initial density frame of expansion section
    i_d_r=-e_d - 2
    # end density frame of expansion section
    e_d_r=-i_d - 1
    # initial density, correct frame
    i_d=int((args.init_dens-dens)*args.frame_jump) + 1
    # end density, correct frame
    e_d=int((args.end_dens-dens)*args.frame_jump) + 2

    # Group particles in clusters
    body_flags=set()

    # Store unique body flags
    for i in t[0].particles.body[:]:
        s=set([i])
        body_flags.update(s)

    # List of indexes of particles of the same body, excluding core
    part_index=[]

    for body in body_flags:
        prelim_indexes=[]
        for index, b_f in enumerate(t[0].particles.body[:]):
            if b_f == body and t[0].particles.typeid[index] == 1:
                prelim_indexes.append(index)

        part_index.append(prelim_indexes)

    part_index=np.array(part_index)

    # #
    # #
    # #
    # #
    # Gyration radio and asphericity
    '''
    dens_c, b_c = gyr_r(snapshots=t[i_d:e_d],
                        vol_tot=vol, all_index=part_index)
    '''
    dens_e, b_e=gyr_r(snapshots=t[i_d_r:e_d_r],
                        vol_tot=vol, all_index=part_index)

    dens_ave=[]
    b_ave=[]
    for i in range(len(dens_c)):
        dens_ave.append((dens_c[i]+dens_e[i])/2)
        b_ave.append((b_c[i]+b_e[i])/2)

    ax1_1.set_title('N-{}'.format(total_N))
    ax1_1.plot(dens_c, b_c, label=label)

    print('\n')

ax1_1.legend()
ax1_1.set_ylabel('$b$ / -')
ax1_1.set_xlabel('$\phi$ / -')
# ax1_1.xaxis.set_major_locator(MultipleLocator(0.01))
plt.xticks(rotation=70)
fig1.tight_layout()
fig1.savefig(dir_name + 'asphericity.pdf', format='pdf')
plt.show()
plt.clf()
