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

# Get path to current working directory
path = os.getcwd()
data_names = []
data_names_ENERGY = []

# Create directories
try:
    if(not os.path.exists("thermo_plots")):
        os.mkdir("thermo_plots")
except OSError as e:
    if e.errno != 17:
        raise
    pass

# Save .log files paths
for root, dirs, files in os.walk(path, topdown=True):
    for name in files:
        if name[-10:] == 'ENERGY.log':
            val = (root.split('/')[-1], name.split('_ENERGY.log')
                   [0], os.path.join(root, name))
            data_names_ENERGY.append(val)

        elif name[-4:] == '.log':
            val = (root.split('/')[-1], name.split('.log')
                   [0], os.path.join(root, name))
            data_names.append(val)

print('Number of .log files = ', len(data_names))

print('Number of ENERGY.log files = ', len(data_names_ENERGY))

# Store combinations of type of cluster (title), number of clusters and particles per cluster (subtitle)
# i[0] shape, i[1](2) N, i[1](5) Nclus
titles_subtitle = set(
    [(i[0], i[1].split('_')[2], i[1].split('_')[5]) for i in data_names])

titles_subtitle_ENERGY = set(
    [(i[0], i[1].split('_')[2], i[1].split('_')[5]) for i in data_names_ENERGY])

# Read and plot data
for t_s in [titles_subtitle, titles_subtitle_ENERGY]:
    if t_s == titles_subtitle:
        print('Combination of .log files = ')
        print(titles_subtitle)
    elif t_s == titles_subtitle_ENERGY:
        print('Combination of ENERGY.log files = ')
        print(titles_subtitle_ENERGY)

    for t, st1, st2 in t_s:
        # Plot Pressure-VF
        fig5 = plt.figure(5)
        ax5_1 = fig5.add_subplot(1, 1, 1)
        ax5_1.set_title('{}\n{}'.format(t, st1))
        for d in data_names:
            # Check if file matches combination of type of cluster, number of clusters and particles per cluster
            if d[0] == t and d[1].split('_')[2] == st1 and d[1].split('_')[5] == st2:
                # Label includes: (5) Nclus, (4) dim, (3) VF, (7) ratio
                label = d[1].split('_')[4] + '_' + \
                    d[1].split('_')[5] + '_' + d[1].split('_')[7]
                data = np.genfromtxt(fname=d[2], skip_header=True)

                # Check dimensions
                dimensions = int(d[1].split('_')[4].split('-')[1])
                total_N = int(d[1].split('_')[2].split('-')[1])

                if d[1].split('_')[1].split('-')[1] == 'one':
                    if dimensions == 3:
                        vol = math.pi/6*1**3 * total_N
                    if dimensions == 2:
                        vol = math.pi/4*1**2 * total_N
                else:
                    poly_key = d[1].split('_')[1].split('-')[1]
                    N_cluster = int(d[1].split('_')[5].split('-')[1])
                    ratio = float(d[1].split('_')[7].split('-')[1])

                    cluster = PartCluster(
                        poly_key=poly_key, N_cluster=N_cluster, halo_diam=1, halo_mass=1, ratio=ratio)

                    vol = cluster.vol_cluster(dimensions) * total_N

                # Averages the last 100 Pressure data points at each VF
                data_ave = [[vol / np.mean(
                    data[i-10:i, 1]), np.mean(data[i-10:i, 9])] for i in
                    range(500, 2803*500, 500)]
                data_ave = np.array(data_ave)

                ax5_1.plot(data_ave[:, 0], data_ave[:, 1], label=label)

                print('!!!!!!!!!!!!!!!!!!!!!\n' +
                      str(d[1]) + '\n' + str(data_ave.shape))

        ax5_1.set_ylabel('Pressure / -')
        ax5_1.set_xlabel('VF / -')
        # ax4_1.xaxis.set_minor_locator(plt.MultipleLocator(500))
        ax5_1.legend()
        # fig4.tight_layout()
        fig5.savefig(
            'thermo_plots/PressureVF_{}_{}_{}.png'.format(t, st1, st2))
        plt.clf()
