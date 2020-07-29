'''
Finds Pressure of first order phase transformation.

Saves data in .npy file
[Nclus, Mean_Pressure, Mean_VF , Start_P, Start_VF, End_P, End_VF, Mean_std_comp, Std_std_comp, Mean_std_exp, Std_std_exp, Window_size]
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
                    dest='in_file', help='input .npy file obtained from PlotPressureRangeVF.py')

args = parser.parse_args()

# Storage file
store_file = 'Plotting_data.npy'
#store_file = 'Testing_data.npy'

# Data fitting parameters
window_size = 47

# Array to store phase transition data into .npy file
plotting_data = np.empty((0, 12))

# #
# #
# #
# Function definitions

# Find start (compression run) and end (expansion run) of phase transformation


def start_phase_change(x_data, y_data):
    phase_change = x_data[0]

    for i in range(len(x_data)):
        bef_grad = np.mean(y_data[i:i+5])
        aft_grad = np.mean(y_data[i+5:i+10])
        porc_diff = abs(aft_grad - bef_grad) / bef_grad

        if porc_diff > 0.1:
            phase_change = x_data[i+5]
            break

    return phase_change

# Find start (expansion run) and end (compression run) of phase transformation


def end_phase_change(x_data, y_data):
    phase_change = x_data[-1]

    for i in range(len(x_data), 0, -1):
        bef_grad = np.mean(y_data[i-5:i])
        aft_grad = np.mean(y_data[i-10:i-5])
        porc_diff = abs(aft_grad - bef_grad) / bef_grad

        if porc_diff > 0.1:
            phase_change = x_data[i-5]
            break

    return phase_change


# #
# #
# #
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

    if Nclus == 0:
        Nclus = 1

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
    x_infl_p3 = p_prime_prime.r  # -2*z[1]/(6*z[0])

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
    x_infl_comp = optimize.newton(
        spline_df_df_comp, x_infl_p3, fprime=spline_df_df_df_comp)
    x_infl_exp = optimize.newton(
        spline_df_df_exp, x_infl_p3, fprime=spline_df_df_df_exp)

    # #
    # #
    # #
    # Find start and end of phase transformation
    x_start_comp = start_phase_change(x_comp, df_comp)
    x_start_exp = end_phase_change(x_exp, df_exp)
    x_end_comp = end_phase_change(x_comp, df_comp)
    x_end_exp = start_phase_change(x_exp, df_exp)

    # #
    # #
    # #
    # Save data
    # Spread of pressure data
    mean_std_comp = np.mean(std_comp)
    mean_std_exp = np.mean(std_exp)
    std_std_comp = np.std(std_comp)
    std_std_exp = np.std(std_exp)

    # save_data = [Nclus, Mean_Pressure, Mean_VF , Start_P, Start_VF, End_P, End_VF, Window_size, Mean_std_comp, Std_std_comp, Mean_std_exp, Std_std_exp]
    mean_P = (spline_y_filtered_comp(x_infl_comp) +
              spline_y_filtered_exp(x_infl_exp)) / 2
    mean_P = float(mean_P)

    mean_VF = (x_infl_comp + x_infl_exp) / 2
    mean_VF = float(mean_VF)

    start_P = (spline_y_filtered_comp(x_start_comp) +
               spline_y_filtered_exp(x_start_exp)) / 2

    start_VF = (x_start_comp + x_start_exp) / 2

    end_P = (spline_y_filtered_comp(x_end_comp) +
             spline_y_filtered_exp(x_end_exp)) / 2

    end_VF = (x_end_comp + x_end_exp) / 2

    save_data = np.array(
        [[Nclus, mean_P, mean_VF, start_P, start_VF, end_P, end_VF, mean_std_comp, std_std_comp, mean_std_exp, std_std_exp, window_size]])

    # Array used to store data in .npy file
    plotting_data = np.concatenate((plotting_data, save_data), axis=0)

    # #
    # #
    # #
    # Plot data
    # Pressure - VF
    fig1 = plt.figure(1)
    ax1_1 = fig1.add_subplot(1, 1, 1)
    '''
    ax1_1.set_title('Nclus-{}'.format(Nclus))
    '''
    ax1_1.errorbar(x_comp, y_comp, yerr=std_comp,
                   linewidth=0.75, label='Original Data COMPRESSION')
    ax1_1.errorbar(x_exp, y_exp, yerr=std_exp, linewidth=0.75,
                   label='Original Data EXPANSION')
    ax1_1.plot(x_comp, y_filtered_comp,
               label='Filtered COMPRESSION', linewidth=0.5)
    ax1_1.plot(x_exp, y_filtered_exp,
               label='Filtered EXPANSION', linewidth=0.5)
    ax1_1.plot(x, py, label='Cubic Polynomial Fit', linewidth=0.5)
    ax1_1.plot(x_infl_p3, p(x_infl_p3), 'o')  # label='Cubic Inflection')
    ax1_1.plot(x_infl_comp, spline_y_filtered_comp(
        x_infl_comp), 'o')  # label='Inflection Newton COMPRESSION')
    ax1_1.plot(x_infl_exp, spline_y_filtered_exp(
        x_infl_exp), 'o')  # label='Inflection Newton EXPANSION')
    ax1_1.plot(x_start_comp, spline_y_filtered_comp(x_start_comp), 'x', x_end_comp,
               spline_y_filtered_comp(x_end_comp), 'x')  # label='Phase Transition COMPRESSION')
    ax1_1.plot(x_start_exp, spline_y_filtered_exp(x_start_exp), 'x', x_end_exp,
               spline_y_filtered_exp(x_end_exp), 'x')  # label='Phase Transition EXPANSION')
    ax1_1.set_ylabel('$\dfrac{PA_n}{kT}$ / -')
    ax1_1.set_xlabel('$\phi$ / -')
    ax1_1.legend()
    fig1.tight_layout()

    # plt.show()

    # First derivative
    fig2 = plt.figure(2)
    ax2_1 = fig2.add_subplot(1, 1, 1)
    ax2_1.set_title(
        'Nclus-{}\nFIRST DERIVATIVE'.format(Nclus))
    ax2_1.plot(x_comp, df_comp,
               '--', label='COMPRESSION')
    ax2_1.plot(x_exp, df_exp, '--',
               label='EXPANSION')
    ax2_1.plot(x_start_comp, spline_df_comp(x_start_comp), 'x', x_end_comp,
               spline_df_comp(x_end_comp), 'x', label='Phase Transition COMPRESSION')
    ax2_1.plot(x_start_exp, spline_df_exp(x_start_exp), 'x', x_end_exp,
               spline_df_exp(x_end_exp), 'x', label='Phase Transition EXPANSION')
    ax2_1.legend()
    fig2.tight_layout()
    # plt.show()

    # Second derivative
    fig3 = plt.figure(3)
    ax3_1 = fig3.add_subplot(1, 1, 1)
    ax3_1.set_title(
        'Nclus-{}\nSECOND DERIVATIVE'.format(Nclus))
    ax3_1.plot(x_comp, df_df_comp,
               '--', label='COMPRESSION')
    ax3_1.plot(x_exp, df_df_exp, '--',
               label='EXPANSION')
    ax3_1.legend()
    fig3.tight_layout()
    # plt.show()

    # Create File to save figures
    with PdfPages("{}_Figures.pdf".format(d[:-4])) as pdf:
        for i in range(1, 4):
            pdf.savefig(figure=i)

    fig1.clear()
    fig2.clear()
    fig3.clear()


if(os.path.exists(store_file)):
    read_data = np.load(store_file)
    plotting_data = np.concatenate((read_data, plotting_data), axis=0)
    print('Plotting_data')
    print(plotting_data)
    np.save(store_file, plotting_data)
else:
    print('Plotting_data')
    print(plotting_data)
    np.save(store_file, plotting_data)
