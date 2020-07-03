'''
Finds Pressure of first order phase transformation.

Saves data in .npy file
[Nclus, Mean_Pressure, Mean_VF , Start_P, Start_VF, End_P, End_VF, Window_size]
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
    description='Finds Pressure of first order phase transformation.')
parser.add_argument('-f', '--in-file', type=str, nargs='+',
                    dest='in_file', help='input .gsd file')

args = parser.parse_args()

# Data fitting parameters
window_size = 67

# Array to store phase transition data from different files
plotting_data = np.empty((0, 7))

# Loops over given files
for d in args.in_file:
    print('Working on:')
    print(d)
    print('')

    # Read data
    data = np.load(d)
    x = data[:, 0]
    y = data[:, 1]
    std = data[:, 2]
    Nclus = int(d.split('_')[5].split('-')[1])

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

    # #
    # #
    # #
    # # Find inflection point using Central Differences and Newton-Rapson method
    # Filter data using Savitzky–Golay filter
    y_filtered_comp = scipy.signal.savgol_filter(
        y_comp, window_size, 3, delta=x[1]-x[0])  # window size 10, polynomial order 3

    y_filtered_exp = scipy.signal.savgol_filter(
        y_exp, window_size, 3, delta=x[1]-x[0])  # window size 10, polynomial order 3

    # Interpolate using spline on filtered data
    spline_y_filtered_comp = interp1d(x_comp, y_filtered_comp, kind='cubic')
    spline_y_filtered_exp = interp1d(x_exp, y_filtered_exp, kind='cubic')

    # First derivative of filtered signal
    df_comp = np.gradient(y_filtered_comp, x_comp)
    df_exp = np.gradient(y_filtered_exp, x_exp)

    spline_df_comp = interp1d(x_comp, df_comp, kind='cubic')
    spline_df_exp = interp1d(x_exp, df_exp, kind='cubic')

    # Second derivative of filtered signal
    df_df_comp = np.gradient(df_comp, x_comp)
    df_df_exp = np.gradient(df_exp, x_exp)

    spline_df_df_comp = interp1d(x_comp, df_df_comp, kind='cubic')
    spline_df_df_exp = interp1d(x_exp, df_df_exp, kind='cubic')

    # Third derivative of filtered signal
    df_df_df_comp = np.gradient(df_df_comp, x_comp)
    df_df_df_exp = np.gradient(df_df_exp, x_exp)

    spline_df_df_df_comp = interp1d(x_comp, df_df_df_comp, kind='cubic')
    spline_df_df_df_exp = interp1d(x_exp, df_df_df_exp, kind='cubic')

    # Find inflection point
    root_comp = optimize.newton(
        spline_df_df_comp, x_infl, fprime=spline_df_df_df_comp)
    root_exp = optimize.newton(
        spline_df_df_exp, x_infl, fprime=spline_df_df_df_exp)

    #
    #
    #
    #
    #
    #
    #

    # Find start and end of phase transformation
    def start_phase_change(x_data, y_data):
        phase_change = x_data[0]

        for i in range(len(x_data)):
            bef_grad = np.mean(y_data[i:i+10])
            aft_grad = np.mean(y_data[i+10:i+20])
            porc_diff = abs(aft_grad - bef_grad) / bef_grad

            if porc_diff > 0.3:
                phase_change = x_data[i]
                break

        return phase_change

    def end_phase_change(x_data, y_data):
        phase_change = x_data[-1]

        for i in range(len(x_data), 0, -1):
            bef_grad = np.mean(y_data[i-10:i])
            aft_grad = np.mean(y_data[i-20:i-10])
            porc_diff = abs(aft_grad - bef_grad) / bef_grad

            if porc_diff > 0.3:
                phase_change = x_data[i-1]
                break

        return phase_change

    x_start_comp = start_phase_change(x_comp, df_comp)
    x_start_exp = end_phase_change(x_exp, df_exp)
    x_end_comp = end_phase_change(x_comp, df_comp)
    x_end_exp = start_phase_change(x_exp, df_exp)
    #
    #
    #
    #
    #
    #
    #

    fig9 = plt.figure(2)
    ax9_1 = fig9.add_subplot(1, 1, 1)
    ax9_1.set_title(
        'Nclus-{}\nFIRST DERIVATIVE'.format(int(d.split('_')[5].split('-')[1])))
    ax9_1.plot(x_comp, df_comp,
               '--', label='CD Derivative COMPRESSION')
    ax9_1.plot(x_exp, df_exp, '--',
               label='CD Derivative EXPANSION')
    ax9_1.plot(x_start_comp, spline_comp(x_start_comp), 'x', x_end_comp,
               spline_comp(x_end_comp), 'x', label='Phase Transition COMPRESSION')
    ax9_1.plot(x_start_exp, spline_exp(x_start_exp), 'x', x_end_exp,
               spline_exp(x_end_exp), 'x', label='Phase Transition EXPANSION')
    ax9_1.legend()
    # plt.show()

    # # Find area between splines and a line that crosses respective inflection point
    a_1_comp_ = []  # Area before inflection point
    a_2_comp_ = []  # Area after inflection point
    diff_a_comp_ = []  # Difference between areas

    a_1_exp_ = []  # Area before inflection point
    a_2_exp_ = []  # Area after inflection point
    diff_a_exp_ = []  # Difference between areas

    m_ = []  # Slope of the curve

    m = 1

    # Data points of the line
    y_lin_comp = m*x_comp + \
        (spline_y_filtered_comp(root_comp) - m*root_comp)
    y_lin_exp = m*x_exp + (spline_y_filtered_exp(root_exp) - m*root_exp)

    # Difference of line and splines of comp and exp data
    diff_spline_y_comp = y_lin_comp - spline_y_filtered_comp(x_comp)

    diff_spline_y_exp = y_lin_exp - spline_y_filtered_exp(x_exp)

    # Interpolation function used to find zeros diff_spline_y
    spline_diff_spline_y_comp = interp1d(
        x_comp, diff_spline_y_comp, kind='cubic')
    spline_diff_spline_y_exp = interp1d(x_exp, diff_spline_y_exp, kind='cubic')

    # First intercept of line and splines of comp and exp data
    if np.sign(spline_diff_spline_y_comp(root_comp*0.999)) == np.sign(spline_diff_spline_y_comp(x_comp[0])):
        x_1_comp = root_comp
    else:
        #x_1_comp = optimize.newton(spline_diff_spline_y_comp, x_comp[0] + (root_comp-x_comp[0])/2)
        x_1_comp = optimize.brentq(
            spline_diff_spline_y_comp, x_comp[0], root_comp*0.999)

    if np.sign(spline_diff_spline_y_exp(root_exp*0.999)) == np.sign(spline_diff_spline_y_exp(x_exp[-1])):
        x_1_exp = root_exp
    else:
        #x_1_exp = optimize.newton(spline_diff_spline_y_exp, x_exp[-1] + (root_exp-x_exp[-1])/2)
        x_1_exp = optimize.brentq(
            spline_diff_spline_y_exp, x_exp[-1], root_exp*0.999)

    # Second intercept of line and splines of comp and exp data
    if np.sign(spline_diff_spline_y_comp(root_comp*1.001)) == np.sign(spline_diff_spline_y_comp(x_comp[-1])):
        x_2_comp = root_comp
    else:
        #x_2_comp = optimize.newton(spline_diff_spline_y_comp, x_comp[-1] - (x_comp[-1]-root_comp)/2)
        x_2_comp = optimize.brentq(
            spline_diff_spline_y_comp, root_comp*1.001, x_comp[-1])

    if np.sign(spline_diff_spline_y_exp(root_exp*1.001)) == np.sign(spline_diff_spline_y_exp(x_exp[0])):
        x_2_exp = root_exp
    else:
        #x_2_exp = optimize.newton(spline_diff_spline_y_exp, x_exp[0] - (x_exp[0]-root_exp)/2)
        x_2_exp = optimize.brentq(
            spline_diff_spline_y_exp, root_exp*1.001, x_exp[0])

    # Area in [x1,root] => line - spline
    ind_comp = np.where((x_comp >= x_1_comp) & (x_comp <= root_comp))
    ind_exp = np.where((x_exp >= x_1_exp) & (x_exp <= root_exp))

    a_1_comp = np.trapz(spline_y_filtered_comp(
        x_comp[ind_comp]), x_comp[ind_comp]) - np.trapz(y_lin_comp[ind_comp], x_comp[ind_comp])
    a_1_exp = np.trapz(spline_y_filtered_exp(np.flip(x_exp[ind_exp])), np.flip(
        x_exp[ind_exp])) - np.trapz(np.flip(y_lin_exp[ind_exp]), np.flip(x_exp[ind_exp]))

    # print('a_1_comp={}\na_1_exp={}'.format(a_1_comp,a_1_exp))

    # Area in [root,x2] => spline - line
    ind_comp = np.where((x_comp <= x_2_comp) & (x_comp >= root_comp))
    ind_exp = np.where((x_exp <= x_2_exp) & (x_exp >= root_exp))
    # print(root_exp)
    # print(x_2_exp)
    # print(x_exp[ind_exp])

    a_2_comp = -np.trapz(spline_y_filtered_comp(
        x_comp[ind_comp]), x_comp[ind_comp]) + np.trapz(y_lin_comp[ind_comp], x_comp[ind_comp])
    a_2_exp = -np.trapz(spline_y_filtered_exp(np.flip(x_exp[ind_exp])), np.flip(
        x_exp[ind_exp])) + np.trapz(np.flip(y_lin_exp[ind_exp]), np.flip(x_exp[ind_exp]))

    # Append values to lists
    a_1_comp_.append(a_1_comp)
    a_2_comp_.append(a_2_comp)
    diff_a_comp_.append(a_1_comp - a_2_comp)

    a_1_exp_.append(a_1_exp)
    a_2_exp_.append(a_2_exp)
    diff_a_exp_.append(a_1_exp - a_2_exp)

    m_.append(m)

    # print('a_2_comp={}\na_2_exp={}'.format(a_2_comp,a_2_exp))

    # print('x_1_comp={}\nroot_comp={}\nx_2_comp={}\nx_1_exp={}\nroot_exp={}\nx_2_exp{}'.format(x_1_comp,root_comp,x_2_comp,x_1_exp,root_exp,x_2_exp))

    # # Plot Diff between spline and y
    fig3 = plt.figure(3)
    ax3_1 = fig3.add_subplot(1, 1, 1)
    ax3_1.set_title('Nclus-{}'.format(int(d.split('_')[5].split('-')[1])))
    ax3_1.plot(x_comp, diff_spline_y_comp, label='COMPRESSION')
    ax3_1.plot(x_exp, diff_spline_y_exp, label='EXPANSION')
    ax3_1.plot(x_1_comp, spline_diff_spline_y_comp(x_1_comp), 'o', x_2_comp,
               spline_diff_spline_y_comp(x_2_comp), 'o', label='Intercept COMPRESSION')
    ax3_1.plot(x_1_exp, spline_diff_spline_y_exp(x_1_exp), 'o', x_2_exp,
               spline_diff_spline_y_exp(x_2_exp), 'o', label='Intercept EXPANSION')

    ax3_1.set_ylabel('Diff Splice - YLine / -')
    ax3_1.set_xlabel('$\phi$ / -')
    ax3_1.legend()
    fig3.tight_layout()

    # plt.show()

    # # Plot Pressure - VF
    fig1 = plt.figure(4)
    ax1_1 = fig1.add_subplot(1, 1, 1)
    ax1_1.set_title('Nclus-{}'.format(int(d.split('_')[5].split('-')[1])))

    '''
    ax1_1.errorbar(x_comp, y_comp, yerr=std_comp,
                   linewidth=0.5, label='Original Data COMPRESSION')
    ax1_1.errorbar(x_exp, y_exp, yerr=std_exp, linewidth=0.5,
                   label='Original Data EXPANSION')
    '''
    ax1_1.plot(x_comp, y_filtered_comp,
               label='Filtered COMPRESSION', linewidth=0.5)
    ax1_1.plot(x_exp, y_filtered_exp,
               label='Filtered EXPANSION', linewidth=0.5)
    ax1_1.plot(x, py, label='Cubic Polynomial Fit', linewidth=0.5)
    ax1_1.plot(x_infl, p(x_infl), 'o', label='Cubic Inflection')
    # ax1_1.plot(x[:int(len(x)/2)], dfFFT.real, '--', label="FFT Derivative")
    #ax1_1.plot(x_comp, df_df_comp, '--', label='CD Derivative COMPRESSION')
    #ax1_1.plot(x_exp, df_df_exp, '--', label='CD Derivative EXPANSION')
    ax1_1.plot(root_comp, spline_y_filtered_comp(
        root_comp), 'o', label='Inflection Newton COMPRESSION')
    ax1_1.plot(root_exp, spline_y_filtered_exp(
        root_exp), 'o', label='Inflection Newton EXPANSION')
    ax1_1.plot(x_comp, y_lin_comp, label='Line COMP', linewidth=0.5)
    ax1_1.plot(x_exp, y_lin_exp, label='Line EXP', linewidth=0.5)
    ax1_1.plot(x_1_comp, spline_y_filtered_comp(x_1_comp), 'o', x_2_comp,
               spline_y_filtered_comp(x_2_comp), 'o', label='Intercept COMPRESSION')
    ax1_1.plot(x_1_exp, spline_y_filtered_exp(x_1_exp), 'o', x_2_exp,
               spline_y_filtered_exp(x_2_exp), 'o', label='Intercept EXPANSION')
    ax1_1.plot(x_start_comp, spline_y_filtered_comp(x_start_comp), 'x', x_end_comp,
               spline_y_filtered_comp(x_end_comp), 'x', label='Phase Transition COMPRESSION')
    ax1_1.plot(x_start_exp, spline_y_filtered_exp(x_start_exp), 'x', x_end_exp,
               spline_y_filtered_exp(x_end_exp), 'x', label='Phase Transition EXPANSION')

    ax1_1.set_ylabel('Pressure / -')
    ax1_1.set_xlabel('$\phi$ / -')
    ax1_1.legend()
    fig1.tight_layout()

    # plt.show()

    # Plot second derivative
    fig9 = plt.figure(1)
    ax9_1 = fig9.add_subplot(1, 1, 1)
    ax9_1.set_title(
        'Nclus-{}\nSECOND DERIVATIVE'.format(int(d.split('_')[5].split('-')[1])))
    ax9_1.plot(x_comp, df_df_comp,
               '--', label='CD Derivative COMPRESSION')
    ax9_1.plot(x_exp, df_df_exp, '--',
               label='CD Derivative EXPANSION')
    ax9_1.legend()
    # plt.show()

    # save_data = [Nclus, Mean_Pressure, Mean_VF , Start_P, Start_VF, End_P, End_VF, Window_size]
    mean_P = (spline_y_filtered_comp(root_comp) +
              spline_y_filtered_exp(root_exp)) / 2
    mean_P = float(mean_P)

    mean_VF = (root_comp + root_exp) / 2
    mean_VF = float(mean_VF)

    start_P = (spline_y_filtered_comp(x_start_comp) +
               spline_y_filtered_exp(x_start_exp)) / 2

    start_VF = (x_start_comp + x_start_exp) / 2

    end_P = (spline_y_filtered_comp(x_end_comp) +
             spline_y_filtered_exp(x_end_exp)) / 2

    end_VF = (x_end_comp + x_end_exp) / 2

    save_data = np.array(
        [[Nclus, mean_P, mean_VF, start_P, start_VF, end_P, end_VF, window_size]])

    plotting_data = np.concatenate((plotting_data, save_data), axis=0)


# Create File to save figures
with PdfPages("Figures.pdf") as pdf:
    for i in range(1, 5):
        pdf.savefig(figure=i)


print('plotting_data')
print(plotting_data)

if(os.path.exists('Plotting_data.npy')):
    read_data = np.load('Plotting_data.npy')
    plotting_data = np.concatenate((read_data, plotting_data), axis=0)
    print('plotting_data')
    print(plotting_data)
    np.save('Plotting_data.npy', plotting_data)
else:
    np.save('Plotting_data.npy', plotting_data)
