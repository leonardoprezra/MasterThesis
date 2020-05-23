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
        if name[-4:] == '.log':
            # val = [a for i in file.split('_') for a in i.split('-') if a[0]=='tstep']
            val = (root.split('/')[-1], name.split('.log')
                   [0], os.path.join(root, name))
            data_names.append(val)

print(len(data_names))
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
titles_subtitle = set(
    [(i[0], i[1].split('_')[-5], i[1].split('_')[-2]) for i in data_names])

# Read and plot data
print(titles_subtitle)
for t, st1, st2 in titles_subtitle:
    # Plot Temperature
    fig1 = plt.figure(1)
    ax1_1 = fig1.add_subplot(1, 1, 1)
    ax1_1.set_title('{}\n{}\n{}'.format(t, st1, st2))
    for d in data_names:
        # Check if file matches combination of type of cluster, number of clusters and particles per cluster
        if d[0] == t and d[1].split('_')[-5] == st1 and d[1].split('_')[-2] == st2:
            # Label includes: (-2) Nclus, (-3) dim, (-4) VF
            label = d[1].split('_')[-4]+'_' + \
                d[1].split('_')[-3]+'_'+d[1].split('_')[-2]
            data = np.genfromtxt(fname=d[2], skip_header=True)
            ax1_1.plot(data[:, 0], data[:, 8], label=label)

    ax1_1.set_ylabel('Temperature / -')
    ax1_1.set_xlabel('Time Step / -')
    #ax1_1.set_ylim(ymax=10)
    ax1_1.legend()
    fig1.tight_layout()
    fig1.savefig('thermo_plots/Temperature_{}_{}_{}.png'.format(t, st1, st2))
    plt.clf()

    # Plot Pressure
    fig2 = plt.figure(2)
    ax2_1 = fig2.add_subplot(1, 1, 1)
    ax2_1.set_title('{}\n{}\n{}'.format(t, st1, st2))
    for d in data_names:
        # Check if file matches combination of type of cluster, number of clusters and particles per cluster
        if d[0] == t and d[1].split('_')[-5] == st1 and d[1].split('_')[-2] == st2:
            # Label includes: (-2) Nclus, (-3) dim, (-4) VF
            label = d[1].split('_')[-4]+'_' + \
                d[1].split('_')[-3]+'_'+d[1].split('_')[-2]
            data = np.genfromtxt(fname=d[2], skip_header=True)
            ax2_1.plot(data[:, 0], data[:, 9], label=label)

    ax2_1.set_ylabel('Pressure / -')
    ax2_1.set_xlabel('Time Step / -')
    #ax2_1.set_ylim(ymax=10)
    ax2_1.legend()
    fig2.tight_layout()
    fig2.savefig('thermo_plots/Pressure_{}_{}_{}.png'.format(t, st1, st2))
    plt.clf()

    # Plot Internal Energy
    fig3 = plt.figure(3)
    ax3_1 = fig3.add_subplot(1, 1, 1)
    ax3_1.set_title('{}\n{}\n{}'.format(t, st1, st2))
    for d in data_names:
        # Check if file matches combination of type of cluster, number of clusters and particles per cluster
        if d[0] == t and d[1].split('_')[-5] == st1 and d[1].split('_')[-2] == st2:
            # Label includes: (-2) Nclus, (-3) dim, (-4) VF
            label = d[1].split('_')[-4]+'_' + \
                d[1].split('_')[-3]+'_'+d[1].split('_')[-2]
            data = np.genfromtxt(fname=d[2], skip_header=True)
            ax3_1.plot(data[:, 0], data[:, 4]+data[:, 5], label=label)

    ax3_1.set_ylabel('Energy / -')
    ax3_1.set_xlabel('Time Step / -')
    # ax3_1.set_ylim(ymax=10)
    ax3_1.legend()
    fig3.tight_layout()
    fig3.savefig('thermo_plots/Energy{}_{}_{}.png'.format(t, st1, st2))
    plt.clf()

    # Plot WCA energy
    fig1 = plt.figure(4)
    ax1_1 = fig1.add_subplot(1, 1, 1)
    ax1_1.set_title('{}\n{}\n{}'.format(t, st1, st2))
    for d in data_names:
        # Check if file matches combination of type of cluster, number of clusters and particles per cluster
        if d[0] == t and d[1].split('_')[-5] == st1 and d[1].split('_')[-2] == st2:
            # Label includes: (-2) Nclus, (-3) dim, (-4) VF
            label = d[1].split('_')[-4]+'_' + \
                d[1].split('_')[-3]+'_'+d[1].split('_')[-2]
            data = np.genfromtxt(fname=d[2], skip_header=True)
            ax1_1.plot(data[:, 0], data[:, 10], label=label)

    ax1_1.set_ylabel('Pair WCA energy / -')
    ax1_1.set_xlabel('Time Step / -')
    ax1_1.xaxis.set_minor_locator(plt.MultipleLocator(500))
    ax1_1.legend()
    fig1.tight_layout()
    fig1.savefig('thermo_plots/PairWCAEnergy_{}_{}_{}.png'.format(t, st1, st2))
    plt.clf()
