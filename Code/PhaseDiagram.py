'''
Plots Phase Diagram of 2D clusters.


'''

from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import os
import math
from MarsonFunctions import PartCluster
import argparse
import scipy.signal
from scipy.interpolate import interp1d
from scipy import optimize

# Command line argument parsing
parser = argparse.ArgumentParser(
    description='Plots Phase Diagram of 2D clusters.')
parser.add_argument('-f', '--in-file', type=str, nargs='+',
                    dest='in_file', help='input Plotting_data.npy file')

args = parser.parse_args()

# Loops over given files
for d in args.in_file:
    print(d)
    # Read data and store in array
    # [Nclus, Mean_Pressure, Mean_VF , Start_P, Start_VF, End_P, End_VF, Mean_std_comp, Std_std_comp, Mean_std_exp, Std_std_exp, Window_size]
    data = np.load(d)
    Nclus = data[1:, 0]
    mean_std_comp = data[:, 7]
    std_std_comp = data[:, 8]
    mean_std_exp = data[:, 9]
    std_std_exp = data[:, 10]
    start_p = data[1:, 3]
    start_vf = data[1:, 4]
    end_p = data[1:, 5]
    end_vf = data[1:, 6]

    # Reference pressure = Nclus-0
    p_ref = data[0, 1]
    vf_ref = data[0, 2]

    # Relative pressure at phase transition
    p_rel = data[1:, 1]/p_ref
    vf_rel = data[1:, 2]/vf_ref

    # Relative pressure at phase transition
    fig1 = plt.figure(1)
    ax1_1 = fig1.add_subplot(1, 1, 1)
    ax1_1.plot(Nclus, p_rel, label='Relative Pressure')
    ax1_1.set_ylabel('P/P(Nclus-0) / -')
    ax1_1.set_xlabel('Nclus / -')
    ax1_1.set_xticks(np.arange(3, 22))
    ax1_1.legend()

    # Relative density at phase transition
    fig2 = plt.figure(2)
    ax2_1 = fig2.add_subplot(1, 1, 1)
    ax2_1.plot(Nclus, vf_rel, label='Relative Density')
    ax2_1.set_ylabel('$\phi$/$\phi$(Nclus-0) / -')
    ax2_1.set_xlabel('Nclus / -')
    ax2_1.set_xticks(np.arange(3, 22))
    ax2_1.legend()

    # Pressure variation
    fig3 = plt.figure(3)
    ax3_1 = fig3.add_subplot(1, 1, 1)
    ax3_1.errorbar(np.concatenate(
        (np.array([0]), Nclus)), mean_std_comp, yerr=std_std_comp, label='Pressure Variation Compression')
    ax3_1.errorbar(np.concatenate(
        (np.array([0]), Nclus)), mean_std_exp, yerr=std_std_exp, label='Pressure Variation Expansion')
    ax3_1.set_ylabel('std / -')
    ax3_1.set_xlabel('Nclus / -')
    ax3_1.set_xticks(np.concatenate((np.array([0]), np.arange(3, 22))))
    ax3_1.legend()

    # #
    # #
    # #
    # Phase diagram
    fig4 = plt.figure(4)
    ax4_1 = fig4.add_subplot(1, 1, 1)
    # Liquid phase
    for ind, vf in enumerate(start_vf):
        n = Nclus[ind]
        vf_plot = []

        vf_range = np.arange(0.55, vf, 0.01)
        vf_range = list(vf_range)
        for a in vf_range:
            vf_plot.append(a)

        n_plot = [n]*len(vf_plot)
        ax4_1.plot(n_plot, vf_plot, 'bo', markersize=6,
                   label='phase transition')

    # Phase transition
    for ind, vf in enumerate(start_vf):
        n = Nclus[ind]
        vf_plot = []

        vf_range = np.arange(vf, end_vf[ind], 0.01)
        vf_range = list(vf_range)
        for a in vf_range:
            vf_plot.append(a)

        n_plot = [n]*len(vf_plot)
        ax4_1.plot(n_plot, vf_plot, 'y^', markersize=6,
                   label='phase transition')

    # Solid phase
    for ind, vf in enumerate(end_vf):
        n = Nclus[ind]
        vf_plot = []

        vf_range = np.arange(0.78, vf, -0.01)
        vf_range = list(vf_range)
        for a in vf_range:
            vf_plot.append(a)

        n_plot = [n]*len(vf_plot)
        ax4_1.plot(n_plot, vf_plot, 'rs', markersize=6, label='solid')

    ax4_1.set_ylabel('$\phi$ / -')
    ax4_1.set_xlabel('Nclus / -')
    ax4_1.set_xticks(np.arange(3, 22))
    # ax4_1.legend()

    with PdfPages("{}_PhaseTransition.pdf".format(d[:-4])) as pdf:
        for i in range(1, 5):
            pdf.savefig(figure=i)