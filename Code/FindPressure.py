'''
Finds Pressure of first order phase transformation.


'''

from matplotlib import pyplot as plt
import numpy as np
import os
import math
from MarsonFunctions import PartCluster
import argparse

# Command line argument parsing
parser = argparse.ArgumentParser(
    description='Finds Pressure of first order phase transformation.')
parser.add_argument('-f', '--in-file', type=str, nargs='+',
                    dest='in_file', help='input .gsd file')

args = parser.parse_args()

for d in args.in_file:
    # Read data and store in array
    data = np.load(d)
    x = data[:, 0]
    y = data[:, 1]
    std = data[:, 2]

    # Fit data to polynomial degree=3
    # Polynomial coefficients
    z = np.polyfit(x, y, 3)

    # Polynomial class
    # Returns the result of the value evaluated in the polynomial
    p = np.poly1d(z)
    py = p(x)

    # Find inflection point
    p_prime_prime = np.poly1d([6*z[0], 2*z[1]])
    x_infl = p_prime_prime.r  # -2*z[1]/(6*z[0])

    # # Plot Pressure-VF
    fig1 = plt.figure(1)
    ax1_1 = fig1.add_subplot(1, 1, 1)

    ax1_1.errorbar(x, y, yerr=std, linewidth=0.5, label='Original Data')
    ax1_1.plot(x, py, label='Fit', linewidth=0.5)
    ax1_1.plot(x_infl, p(x_infl), 'o', label='Inflection')

    ax1_1.set_ylabel('Pressure / -')
    ax1_1.set_xlabel('$\phi$ / -')
    ax1_1.legend()
    fig1.tight_layout()

    plt.show()

    # Area differences
    a_1 = []
    a_2 = []
    m_ = []
    diff_a = []

    # # Find intercepts
    for m in np.linspace(0, 2, 100):  # [0, 1, 2]:
        # Coefficients of the polynomial of the intercepts
        coeff_inter = [z[0], z[1], z[2]-m, z[3] - (p(x_infl)-m*x_infl)]

        # Equation of the intercepts
        p_inter = np.poly1d(coeff_inter)

        # Intercepts [x1,xp,x2]
        x_inter = p_inter.r

        # First section area
        coeff_a_1 = [1/4*z[0], 1/3*z[1], 1/2 *
                     (z[2]-m), (z[3] - (p(x_infl)-m*x_infl)), 0]
        p_a_1 = np.poly1d(coeff_a_1)

        # Second section area
        coeff_a_2 = [-1/4*z[0], -1/3*z[1], -1/2 *
                     (z[2]-m), -(z[3] - (p(x_infl)-m*x_infl)), 0]
        p_a_2 = np.poly1d(coeff_a_2)

        # Append values to lists
        a_1.append(p_a_1(x_inter[1]) - p_a_1(x_inter[0]))
        a_2.append(p_a_2(x_inter[2]) - p_a_2(x_inter[1]))
        print('Intercept='+str(x_inter))
        print('a_1={}       a_2={}'.format(a_1[-1], a_2[-1]))
        diff_a.append(a_1[-1] - a_2[-1])
        m_.append(m)

        '''

        # Plot intercepts
        print(x_inter)
        x_lin = np.poly1d([m, p(x_infl)-m*x_infl])

        fig1 = plt.figure(1)
        ax1_1 = fig1.add_subplot(1, 1, 1)

        ax1_1.errorbar(x, y, yerr=std, linewidth=0.5, label='Original Data')
        ax1_1.plot(x, py, label='Fit', linewidth=0.5)
        ax1_1.plot(x_infl, p(x_infl), 'o', label='Inflection')

        x_plot_inter = [x_inter[0], x_inter[2]]
        ax1_1.plot(x_plot_inter, p(x_plot_inter), 'o', label='Intercepts')
        ax1_1.plot(x, x_lin(x), label='m-{}'.format(m))

        ax1_1.set_ylabel('Pressure / -')
        ax1_1.set_xlabel('$\phi$ / -')
        ax1_1.legend()
        fig1.tight_layout()

        plt.show()
        '''

    # # Find start and end of co-existing phases
    # diff_a = np.array(diff_a)

    # Find where m yields diff=0
    index_diff_0 = [count for count, value in enumerate(
        diff_a) if value == 0 and value != 'numpy.complex128']

    start_index = 0
    co = 0
    for count, value in enumerate(diff_a):
        if type(value[0]) != complex:
            print('!!!!!!!!!!!!!!!!!!!Value')
            print(type(value[0]))
            start_index = count
            print('!!!!!!!!!!!!!!!start_index')
            print(start_index)
            break

        co += 1

    print("!!!!!!!!!!!!!!!!!co={}".format(co))

    a_1 = a_1[start_index:]
    a_2 = a_2[start_index:]
    m_ = m_[start_index:]
    diff_a = diff_a[start_index:]

    print('!!!!!!!!!!!!!!!!!!!!!a_1')
    print(a_1)
    print('!!!!!!!!!!!!!!!!!!!!!a_2')
    print(a_2)
    print('!!!!!!!!!!!!!!!!!!!!!m_')
    print(m_)
    print('!!!!!!!!!!!!!!!!!!!!!diff_a')
    print(diff_a)

    # Find where m yields diff=0
    index_diff_0 = [count for count, value in enumerate(
        diff_a) if value == 0 and value != complex]

    fig2 = plt.figure(2)
    ax2_1 = fig2.add_subplot(1, 1, 1)
    ax2_1.plot(m_, a_1, label='a_1', linewidth=0.5)
    ax2_1.plot(m_, a_2, label='a_2', linewidth=0.5)
    ax2_1.plot(m_, diff_a, label='diff', linewidth=0.5)
    ax2_1.plot(m_[index_diff_0[0]], diff_a[index_diff_0[0]],
               'o', label='First diff=0', linewidth=0.5)
    ax2_1.set_ylabel('m / -')
    ax2_1.set_xlabel('Area / -')
    ax2_1.legend()

    fig2.tight_layout()
    plt.show()
'''
    ax1_1.errorbar(x, y, yerr=std, linewidth=0.5, label='Original Data')
    ax1_1.plot(x, py, label='Fit', linewidth=0.5)
    ax1_1.plot(x_infl, p(x_infl), 'o', label='Inflection')

    print('!!!!!m_' + str(m_))
    x_lin = np.poly1d([m_[index_diff_0[0]], p(
        x_infl)-m_[index_diff_0[0]]*x_infl])
    print(x_lin(x))
    ax1_1.plot(x, x_lin(x), label='m-{}'.format(m_[index_diff_0[0]]))

    ax1_1.set_ylabel('Pressure / -')
    ax1_1.set_xlabel('$\phi$ / -')
    ax1_1.legend()
    fig1.tight_layout()

    plt.show()
'''
