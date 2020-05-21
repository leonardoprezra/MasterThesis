'''
Plots variables from log files in current directory.

Run in a directory that contain the data_* directories,
no sub directories should exist in data_*.


'''

from matplotlib import pyplot as plt
import numpy as np
import os

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

'''
# Store combinations of type of cluster (title) and number of clusters (subtitle)
titles_subtitle = set([(i[0], i[1].split('_')[-5]) for i in data_names])

# Read and plot data
print(titles_subtitle)
for t, st in titles_subtitle:
    # Plot Temperature
    fig1 = plt.figure(1)
    ax1_1 = fig1.add_subplot(1, 1, 1)
    ax1_1.set_title('{}\n{}'.format(t, st))
    for d in data_names:
        # Check if file matches combination of type of cluster and number of clusters
        if d[0] == t and d[1].split('_')[-5] == st:
            # Label includes: (-2) Nclus, (-3) dim, (-4) VF
            label = d[1].split('_')[-4]+'_' + \
                d[1].split('_')[-3]+'_'+d[1].split('_')[-2]
            data = np.genfromtxt(fname=d[2], skip_header=True)
            ax1_1.plot(data[:, 0], data[:, 8], label=label)

    ax1_1.set_ylabel('Temperature / -')
    ax1_1.set_xlabel('Time Step / -')
    ax1_1.legend()
    fig1.tight_layout()
    fig1.savefig('thermo_plots/Temperature_{}_{}.png'.format(t, st))
    plt.clf()

    # Plot Pressure
    fig2 = plt.figure(2)
    ax2_1 = fig2.add_subplot(1, 1, 1)
    ax2_1.set_title('{}\n{}'.format(t, st))
    for d in data_names:
        # Check if file matches combination of type of cluster and number of clusters
        if d[0] == t and d[1].split('_')[-5] == st:
            # Label includes: (-2) Nclus, (-3) dim, (-4) VF
            label = d[1].split('_')[-4]+'_' + \
                d[1].split('_')[-3]+'_'+d[1].split('_')[-2]
            data = np.genfromtxt(fname=d[2], skip_header=True)
            ax2_1.plot(data[:, 0], data[:, 9], label=label)

    ax2_1.set_ylabel('Pressure / -')
    ax2_1.set_xlabel('Time Step / -')
    ax2_1.legend()

    
    fig2.tight_layout()
    fig2.savefig('thermo_plots/Pressure_{}_{}.png'.format(t, st))
    plt.clf()
'''

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
        # Plot Temperature
        fig1 = plt.figure(1)
        ax1_1 = fig1.add_subplot(1, 1, 1)
        ax1_1.set_title('{}\n{}\n{}'.format(t, st1, st2))
        for d in data_names:
            # Check if file matches combination of type of cluster, number of clusters and particles per cluster
            if d[0] == t and d[1].split('_')[2] == st1 and d[1].split('_')[5] == st2:
                # Label includes: (5) Nclus, (4) dim, (3) VF, (7) ratio
                label = d[1].split('_')[3]+'_' + \
                    d[1].split('_')[4]+'_'+d[1].split('_')[5] + \
                    '_'+d[1].split('_')[7]
                data = np.genfromtxt(fname=d[2], skip_header=True)
                ax1_1.plot(data[:, 0], data[:, 8], label=label)

        ax1_1.set_ylabel('Temperature / -')
        ax1_1.set_xlabel('Time Step / -')
        # ax1_1.set_ylim(ymax=10)
        ax1_1.legend()
        # fig1.tight_layout()
        fig1.savefig(
            'thermo_plots/Temperature_{}_{}_{}.png'.format(t, st1, st2))
        plt.clf()

        # Plot Pressure
        fig2 = plt.figure(2)
        ax2_1 = fig2.add_subplot(1, 1, 1)
        ax2_1.set_title('{}\n{}\n{}'.format(t, st1, st2))
        for d in data_names:
            # Check if file matches combination of type of cluster, number of clusters and particles per cluster
            if d[0] == t and d[1].split('_')[2] == st1 and d[1].split('_')[5] == st2:
                # Label includes: (5) Nclus, (4) dim, (3) VF, (7) ratio
                label = d[1].split('_')[3]+'_' + \
                    d[1].split('_')[4]+'_'+d[1].split('_')[5] + \
                    '_'+d[1].split('_')[7]
                data = np.genfromtxt(fname=d[2], skip_header=True)
                ax2_1.plot(data[:, 0], data[:, 9], label=label)

        ax2_1.set_ylabel('Pressure / -')
        ax2_1.set_xlabel('Time Step / -')
        # ax2_1.set_ylim(ymax=10)
        ax2_1.legend()
        # fig2.tight_layout()
        fig2.savefig('thermo_plots/Pressure_{}_{}_{}.png'.format(t, st1, st2))
        plt.clf()

        # Plot Internal Energy
        fig3 = plt.figure(3)
        ax3_1 = fig3.add_subplot(1, 1, 1)
        ax3_1.set_title('{}\n{}\n{}'.format(t, st1, st2))
        for d in data_names:
            # Check if file matches combination of type of cluster, number of clusters and particles per cluster
            if d[0] == t and d[1].split('_')[2] == st1 and d[1].split('_')[5] == st2:
                # Label includes: (5) Nclus, (4) dim, (3) VF, (7) ratio
                label = d[1].split('_')[3]+'_' + \
                    d[1].split('_')[4]+'_'+d[1].split('_')[5] + \
                    '_'+d[1].split('_')[7]
                data = np.genfromtxt(fname=d[2], skip_header=True)
                ax3_1.plot(data[:, 0], data[:, 4]+data[:, 5], label=label)

        ax3_1.set_ylabel('Energy / -')
        ax3_1.set_xlabel('Time Step / -')
        # ax3_1.set_ylim(ymax=10)
        ax3_1.legend()
        # fig3.tight_layout()
        fig3.savefig('thermo_plots/Energy{}_{}_{}.png'.format(t, st1, st2))
        plt.clf()

        # Plot WCA energy
        fig4 = plt.figure(4)
        ax4_1 = fig4.add_subplot(1, 1, 1)
        ax4_1.set_title('{}\n{}\n{}'.format(t, st1, st2))
        for d in data_names:
            # Check if file matches combination of type of cluster, number of clusters and particles per cluster
            if d[0] == t and d[1].split('_')[2] == st1 and d[1].split('_')[5] == st2:
                # Label includes: (5) Nclus, (4) dim, (3) VF, (7) ratio
                label = d[1].split('_')[3]+'_' + \
                    d[1].split('_')[4]+'_'+d[1].split('_')[5] + \
                    '_'+d[1].split('_')[7]
                data = np.genfromtxt(fname=d[2], skip_header=True)
                ax4_1.plot(data[:, 0], data[:, 10], label=label)

        ax4_1.set_ylabel('Pair WCA energy / -')
        ax4_1.set_xlabel('Time Step / -')
        ax4_1.xaxis.set_minor_locator(plt.MultipleLocator(500))
        ax4_1.legend()
        # fig4.tight_layout()
        fig4.savefig(
            'thermo_plots/PairWCAEnergy_{}_{}_{}.png'.format(t, st1, st2))
        plt.clf()
