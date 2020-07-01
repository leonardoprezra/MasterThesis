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
        root_comp), 'o', label='Inflection Newton COMPRESSION')
    ax1_1.plot(root_exp, spline_y_filtered_exp(
        root_exp), 'o', label='Inflection Newton EXPANSION')

    ax1_1.set_ylabel('Pressure / -')
    ax1_1.set_xlabel('$\phi$ / -')
    ax1_1.legend()
    fig1.tight_layout()

    plt.show()

    # # Find area between splines and a line that crosses respective inflection point
    a_1_comp_ = []  # Area before inflection point
    a_2_comp_ = []  # Area after inflection point
    diff_a_comp_ = []  # Difference between areas

    a_1_exp_ = []  # Area before inflection point
    a_2_exp_ = []  # Area after inflection point
    diff_a_exp_ = []  # Difference between areas

    m_ = []  # Slope of the curve
    counter = 0

    for m in np.linspace(0, 20, 1000):  # [0, 1, 2]:
        # Data points of the line
        y_lin_comp = m*x_comp + \
            (spline_y_filtered_comp(root_comp) - m*root_comp)
        y_lin_exp = m*x_exp + (spline_y_filtered_exp(root_exp) - m*root_exp)

        # Difference of line and splines of comp and exp data
        diff_spline_y_comp = y_lin_comp - spline_y_filtered_comp(x_comp)

        diff_spline_y_exp = y_lin_exp - spline_y_filtered_exp(x_exp)

        # Interpolation function used to find zeros diff_spline_y
        spline_diff_spline_y_comp = interp1d(x_comp, diff_spline_y_comp, kind='cubic')
        spline_diff_spline_y_exp = interp1d(x_exp, diff_spline_y_exp, kind='cubic')
        
        # First intercept of line and splines of comp and exp data
        if np.sign(spline_diff_spline_y_comp(root_comp*0.999)) == np.sign(spline_diff_spline_y_comp(x_comp[0])):
            x_1_comp = root_comp
        else:
            #x_1_comp = optimize.newton(spline_diff_spline_y_comp, x_comp[0] + (root_comp-x_comp[0])/2)
            x_1_comp = optimize.brentq(spline_diff_spline_y_comp, x_comp[0], root_comp*0.999)

        if np.sign(spline_diff_spline_y_exp(root_exp*0.999)) == np.sign(spline_diff_spline_y_exp(x_exp[-1])):
            x_1_exp = root_exp
        else:
            #x_1_exp = optimize.newton(spline_diff_spline_y_exp, x_exp[-1] + (root_exp-x_exp[-1])/2)
            x_1_exp = optimize.brentq(spline_diff_spline_y_exp, x_exp[-1], root_exp*0.999)

        # Second intercept of line and splines of comp and exp data
        if np.sign(spline_diff_spline_y_comp(root_comp*1.001)) == np.sign(spline_diff_spline_y_comp(x_comp[-1])):
            x_2_comp = root_comp
        else:
            #x_2_comp = optimize.newton(spline_diff_spline_y_comp, x_comp[-1] - (x_comp[-1]-root_comp)/2)
            x_2_comp = optimize.brentq(spline_diff_spline_y_comp, root_comp*1.001, x_comp[-1])

        if np.sign(spline_diff_spline_y_exp(root_exp*1.001)) == np.sign(spline_diff_spline_y_exp(x_exp[0])):
            x_2_exp = root_exp
        else:
            #x_2_exp = optimize.newton(spline_diff_spline_y_exp, x_exp[0] - (x_exp[0]-root_exp)/2)
            x_2_exp = optimize.brentq(spline_diff_spline_y_exp, root_exp*1.001, x_exp[0])

        # Area in [x1,root] => line - spline
        ind_comp = np.where( (x_comp >= x_1_comp) & (x_comp <= root_comp) )
        ind_exp = np.where( (x_exp >= x_1_exp) & (x_exp <= root_exp) )

        a_1_comp = np.trapz(spline_y_filtered_comp(x_comp[ind_comp]), x_comp[ind_comp]) - np.trapz(y_lin_comp[ind_comp],x_comp[ind_comp])
        a_1_exp = np.trapz(spline_y_filtered_exp(np.flip(x_exp[ind_exp])), np.flip(x_exp[ind_exp])) - np.trapz(np.flip(y_lin_exp[ind_exp]),np.flip(x_exp[ind_exp]))

        #print('a_1_comp={}\na_1_exp={}'.format(a_1_comp,a_1_exp))

        # Area in [root,x2] => spline - line
        ind_comp = np.where( (x_comp <= x_2_comp) & (x_comp >= root_comp) )
        ind_exp = np.where( (x_exp <= x_2_exp) & (x_exp >= root_exp) )
        #print(root_exp)
        #print(x_2_exp)
        #print(x_exp[ind_exp])

        a_2_comp = -np.trapz(spline_y_filtered_comp(x_comp[ind_comp]), x_comp[ind_comp]) + np.trapz(y_lin_comp[ind_comp],x_comp[ind_comp])
        a_2_exp = -np.trapz(spline_y_filtered_exp(np.flip(x_exp[ind_exp])), np.flip(x_exp[ind_exp])) + np.trapz(np.flip(y_lin_exp[ind_exp]),np.flip(x_exp[ind_exp]))

        # Append values to lists
        a_1_comp_.append(a_1_comp)
        a_2_comp_.append(a_2_comp)
        diff_a_comp_.append(a_1_comp - a_2_comp)

        a_1_exp_.append(a_1_exp)
        a_2_exp_.append(a_2_exp)
        diff_a_exp_.append(a_1_exp - a_2_exp)

        m_.append(m)

        #print('a_2_comp={}\na_2_exp={}'.format(a_2_comp,a_2_exp))

        #print('x_1_comp={}\nroot_comp={}\nx_2_comp={}\nx_1_exp={}\nroot_exp={}\nx_2_exp{}'.format(x_1_comp,root_comp,x_2_comp,x_1_exp,root_exp,x_2_exp))
        counter += 1
        if False:
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n{}'.format(m))
            fig3 = plt.figure(3)
            ax3_1 = fig3.add_subplot(1,1,1)
            ax3_1.plot(x_comp, diff_spline_y_comp, label='COMPRESSION')
            ax3_1.plot(x_exp, diff_spline_y_exp, label='EXPANSION')
            ax3_1.plot(x_1_comp, spline_diff_spline_y_comp(x_1_comp), 'o', x_2_comp, spline_diff_spline_y_comp(x_2_comp), 'o', label='Intercept COMPRESSION')
            ax3_1.plot(x_1_exp, spline_diff_spline_y_exp(x_1_exp), 'o', x_2_exp, spline_diff_spline_y_exp(x_2_exp), 'o', label='Intercept EXPANSION')
            
            
            ax3_1.set_ylabel('Diff Splice - YLine / -')
            ax3_1.set_xlabel('$\phi$ / -')
            ax3_1.legend()
            fig3.tight_layout()

            plt.show()


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
                root_comp), 'o', label='Inflection Newton COMPRESSION')
            ax1_1.plot(root_exp, spline_y_filtered_exp(
                root_exp), 'o', label='Inflection Newton EXPANSION')
            ax1_1.plot(x_comp, y_lin_comp, label='COMP', linewidth=0.5)
            ax1_1.plot(x_exp, y_lin_exp, label='EXP', linewidth=0.5)
            ax1_1.plot(x_1_comp, spline_y_filtered_comp(x_1_comp), 'o', x_2_comp, spline_y_filtered_comp(x_2_comp), 'o', label='Intercept COMPRESSION')
            ax1_1.plot(x_1_exp, spline_y_filtered_exp(x_1_exp), 'o', x_2_exp, spline_y_filtered_exp(x_2_exp), 'o', label='Intercept EXPANSION')

            ax1_1.set_ylabel('Pressure / -')
            ax1_1.set_xlabel('$\phi$ / -')
            ax1_1.legend()
            fig1.tight_layout()

            plt.show()
            
        
    
    
    
    # Find where m yields diff=0
    
    fig2 = plt.figure(2)
    ax2_1 = fig2.add_subplot(1, 1, 1)
    ax2_1.plot(m_, diff_a_comp_, label='diff_a_comp', linewidth=0.5)
    ax2_1.plot(m_, diff_a_exp_, label='diff_a_exp', linewidth=0.5)
    ax2_1.set_ylabel('DIFF AREA / -')
    ax2_1.set_xlabel('m / -')
    ax2_1.legend()

    fig2.tight_layout()
    plt.show()

    fig4 = plt.figure(4)
    ax4_1 = fig4.add_subplot(1, 1, 1)
    ax4_1.plot(m_, a_1_comp_, label='a_1_comp', linewidth=0.5)
    ax4_1.plot(m_, a_2_comp_, label='a_2_comp', linewidth=0.5)
    ax4_1.plot(m_, a_1_exp_, label='a_1_exp', linewidth=0.5)
    ax4_1.plot(m_, a_2_exp_, label='a_2_exp', linewidth=0.5)
    ax4_1.set_ylabel('DIFF AREA / -')
    ax4_1.set_xlabel('m / -')
    ax4_1.legend()

    fig4.tight_layout()
    plt.show()