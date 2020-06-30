'''
Finds Pressure of first order phase transformation.


'''

from matplotlib import pyplot as plt
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
    description='Finds Pressure of first order phase transformation.')
parser.add_argument('-f', '--in-file', type=str, nargs='+',
                    dest='in_file', help='input .gsd file')

args = parser.parse_args()

# Central difference


def cent_diff(x, y):
    # Forward difference for the first element of the compression
    dfCD = [(y[1]-y[0])/(x[1]-x[0])]
    print('size df={}'.format(len(dfCD)))

    # Central difference for elements on the compression
    for i in range(1, int(len(x))-1):
        dfCD.append((y[i+1]-y[i-1]) / (x[i+1]-x[i-1]))

    print('size df={}'.format(len(dfCD)))
    # Backward difference for the last element of the compression
    dfCD.append((y[int(len(x))-1]-y[int(len(x))-2]) /
                (x[int(len(x))-1]-x[int(len(x))-2]))

    print('size df={}'.format(len(dfCD)))
    '''
    # Forward difference for the first element of the expansion
    dfCD.append((y[int(len(x)/2)+1]-y[int(len(x)/2)]) /
                (x[int(len(x)/2)+1]-x[int(len(x)/2)]))

    print('size df={}'.format(len(dfCD)))
    # Central difference for elements on the expansion
    for i in range(int(len(x)/2)+1, int(len(x))-1):
        dfCD.append((y[i+1]-y[i-1]) / (x[i+1]-x[i-1]))

    print('size df={}'.format(len(dfCD)))
    # Backward difference for the last element of the expansion
    dfCD.append((y[int(len(x))-1]-y[int(len(x))-2]) /
                (x[int(len(x))-1]-x[int(len(x))-2]))

    print('size df={}'.format(len(dfCD)))
    '''

    return dfCD


for d in args.in_file:
    # Read data and store in array
    data = np.load(d)
    x = data[:, 0]
    y = data[:, 1]
    std = data[:, 2]

    # Divide data in compression and expansion runs
    x_comp = x[: int(len(x)/2)]
    y_comp = y[: int(len(y)/2)]
    std_comp = std[: int(len(std)/2)]

    x_exp = x[int(len(x)/2):]
    y_exp = y[int(len(y)/2):]
    std_exp = std[int(len(std)/2):]

    # Fit data to polynomial degree=3
    # Polynomial coefficients
    z = np.polyfit(np.concatenate((x_comp, x_exp), axis=0),
                   np.concatenate((y_comp, y_exp), axis=0), 3)

    # Polynomial class
    # Returns the result of the value evaluated in the polynomial
    p = np.poly1d(z)
    py = p(np.concatenate((x_comp, x_exp), axis=0))

    # Find inflection point of polynonmial fit
    p_prime_prime = np.poly1d([6*z[0], 2*z[1]])
    x_infl = p_prime_prime.r  # -2*z[1]/(6*z[0])

    # # Find inflection point using Central Differences
    # Filter data using Savitzky–Golay filter
    y_filtered_comp = scipy.signal.savgol_filter(
        y_comp, 29, 3, delta=x[1]-x[0])  # window size 10, polynomial order 3

    y_filtered_exp = scipy.signal.savgol_filter(
        y_exp, 29, 3, delta=x[1]-x[0])  # window size 10, polynomial order 3

    # Interpolate using spline on filtered data
    spline_y_filtered_comp = interp1d(x_comp, y_filtered_comp, kind='cubic')
    spline_y_filtered_exp = interp1d(x_exp, y_filtered_exp, kind='cubic')

    # First derivative of filtered signal
    df_comp = cent_diff(x_comp, y_filtered_comp)
    df_exp = cent_diff(x_exp, y_filtered_exp)

    # Second derivative of filtered signal
    df_df_comp = cent_diff(x_comp, df_comp)
    df_df_exp = cent_diff(x_exp, df_exp)

    # Spline interpolation of second derivative
    spline_df_df_comp = interp1d(x_comp, df_df_comp, kind='cubic')
    spline_df_df_exp = interp1d(x_exp, df_df_exp, kind='cubic')

    # Find inflection point using Newton-Rapson method
    root_comp = optimize.newton(spline_df_df_comp, x_infl)
    root_exp = optimize.newton(spline_df_df_exp, x_infl)

    # # Plot Pressure-VF
    fig1 = plt.figure(1)
    ax1_1 = fig1.add_subplot(1, 1, 1)

    ax1_1.errorbar(x_comp, y_comp, yerr=std_comp,
                   linewidth=0.5, label='Original Data COMPRESSION')
    ax1_1.errorbar(x_exp, y_exp, yerr=std_exp, linewidth=0.5,
                   label='Original Data EXPANSION')
    ax1_1.plot(x_comp, y_filtered_comp,
               label='Filtered COMPRESSION', linewidth=0.5)
    ax1_1.plot(x_exp, y_filtered_exp,
               label='Filtered EXPANSION', linewidth=0.5)
    ax1_1.plot(x, py, label='Fit', linewidth=0.5)
    ax1_1.plot(x_infl, p(x_infl), 'o', label='Inflection')
    # ax1_1.plot(x[:int(len(x)/2)], dfFFT.real, '--', label="FFT Derivative")
    #ax1_1.plot(x_comp, df_df_comp, '--', label='CD Derivative COMPRESSION')
    #ax1_1.plot(x_exp, df_df_exp, '--', label='CD Derivative EXPANSION')
    ax1_1.plot(root_comp, spline_y_filtered_comp(
        root_comp), 'o', label='Newton COMPRESSION')
    ax1_1.plot(root_exp, spline_y_filtered_exp(
        root_exp), 'o', label='Newton EXPANSION')

    ax1_1.set_ylabel('Pressure / -')
    ax1_1.set_xlabel('$\phi$ / -')
    ax1_1.legend()
    fig1.tight_layout()

    plt.show()

    # # Find area between splines and a line that crosses respective inflection point
    a_1_comp = []  # Area before inflection point
    a_2_comp = []  # Area after inflection point
    m_comp = []  # Slope of the curve
    diff_a_comp = []  # Difference between areas

    for m in np.linspace(0, 2, 100):  # [0, 1, 2]:
        # Data points of the line
        y_lin_comp = m*x_comp + \
            (spline_y_filtered_comp(root_comp) - m*root_comp)
        y_lin_exp = m*x_exp + (spline_y_filtered_exp(root_exp) - m*root_exp)

        # Find intercept of line and splines of comp and exp data
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        fig1 = plt.figure(1)
        ax1_1 = fig1.add_subplot(1, 1, 1)

        ax1_1.errorbar(x_comp, y_comp, yerr=std_comp,
                       linewidth=0.5, label='Original Data COMPRESSION')
        ax1_1.errorbar(x_exp, y_exp, yerr=std_exp, linewidth=0.5,
                       label='Original Data EXPANSION')
        ax1_1.plot(x_comp, y_filtered_comp,
                   label='Filtered COMPRESSION', linewidth=0.5)
        ax1_1.plot(x_exp, y_filtered_exp,
                   label='Filtered EXPANSION', linewidth=0.5)
        ax1_1.plot(x, py, label='Fit', linewidth=0.5)
        ax1_1.plot(x_infl, p(x_infl), 'o', label='Inflection')
        # ax1_1.plot(x[:int(len(x)/2)], dfFFT.real, '--', label="FFT Derivative")
        #ax1_1.plot(x_comp, df_df_comp, '--', label='CD Derivative COMPRESSION')
        #ax1_1.plot(x_exp, df_df_exp, '--', label='CD Derivative EXPANSION')
        ax1_1.plot(root_comp, spline_y_filtered_comp(
            root_comp), 'o', label='Newton COMPRESSION')
        ax1_1.plot(root_exp, spline_y_filtered_exp(
            root_exp), 'o', label='Newton EXPANSION')
        ax1_1.plot(x_comp, y_lin_comp, label='COMP', linewidth=0.5)
        ax1_1.plot(x_exp, y_lin_exp, label='EXP', linewidth=0.5)

        ax1_1.set_ylabel('Pressure / -')
        ax1_1.set_xlabel('$\phi$ / -')
        ax1_1.legend()
        fig1.tight_layout()

        plt.show()
        '''
        

        intercept

        # Append values to lists
        a_1.append(p_a_1(x_inter[1]) - p_a_1(x_inter[0]))
        a_2.append(p_a_2(x_inter[2]) - p_a_2(x_inter[1]))
        diff_a.append(a_1[-1] - a_2[-1])
        m_.append(m)
        '''

    '''
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
