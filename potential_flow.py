''' 
Sokratis Anagnostopoulos - Potential Flow
'''

import numpy as np
import matplotlib.pyplot as plt
from math import sin, sinh, acos, cos, ceil

pi = acos(-1.0)

# Globar variables
dx = 0.01  # dx length
val = 120.0  # amplitude

# The tolerance value used throughout this study is 1.0e-4


def numerical(tolerance, dx, val):

    # Numerical solution

    size = int(1.0/dx)  # internal grid size = size - 2 (boundaries)
    ini = val/2  # initial values
    psi = np.zeros((size, size)) + ini  # initialisation

    # Boundary Conditions

    # Since the right and left sections of the initial 2m x 1m geometry are
    # symmetric, the geometry considered for this study is the right 1m x 1m square,
    # in favor of computation time.
    # A constant value of 120 can be set along the left boundary.
    # Symmetry boundary condition is usually applied but in this simple flow 
    # case it is not necessary, since the boundaries are known.

    psi[0, :] = 120  # top boundary
    psi[:, 0] = 120  # left boundary
    psi[:, -1] = 120*np.linspace(1, 0, size)  # right boundary
    psi[-1, :] = 120*np.linspace(1, 0, size)  # bottom boundary


    # Creation of i and j arrays that will be used in
    # the vectorised checkerboard method.
    # E.g. for a 3 x 3 psi array, the array i and j
    # are of the form: i = 1 1 1 2 2 2 3 3 3
    #                  j = 1 3 5 2 4 6 1 3 5
    #               or j = 2 4 6 1 3 5 2 4 6
    # Note that array j starts from even or odd sequences
    # depending on the iteration so that all the checkerboard
    # cells can be accessed.

    odd = np.arange(1, size-1, 2)
    even = np.arange(2, size-1, 2)
    even2_0 = np.hstack((odd,even))
    even2 = np.copy(even2_0)
    even3_0 = np.hstack((even,odd))
    even3 = np.copy(even3_0)
    odd2 = np.ones(ceil((size-2)/2))
    odd3 = np.copy(odd2)
    for _ in range(ceil((size-4)/2)):
        even2 = np.hstack((even2, even2_0))
        even3 = np.hstack((even3, even3_0))
    for i in range(1, size-2):
        odd3 = np.hstack((odd3, odd2 + i))
    odd3 = odd3.astype(int)

    i = odd3  # Sets array of all i. Note that it doesn't have to be 
              # inside the loop while j does.

    # Constants
    a = b = c = d = 1
    e = -4
    f = 0

    # tolerance = 1.0e-4  # convergence tolerance

    resid = np.zeros((size,size))  # Initialisation of residual array
    sum_resid = 1  # Sum of the residuals
    count = size*size  # counter
    omega = 1.5  # relaxation factor
    iteration = 0  # iterations counter

    # Numerical Calculations
    while (sum_resid > tolerance):
        sum_resid = 0

        j = even3
        resid[i, j] = a*psi[i+1, j] + b*psi[i-1, j] + c*psi[i, j+1] + d*psi[i, j-1] + e*psi[i, j] - f  # Residual
        psi[i, j] = psi[i, j] - omega*resid[i, j]/e  #Psi calculation, considering relaxation factor.
        sum_resid += np.sum(np.abs(resid))

        j = even2
        resid[i, j] = a*psi[i+1, j] + b*psi[i-1, j] + c*psi[i, j+1] + d*psi[i, j-1] + e*psi[i, j] - f  # Residual
        psi[i, j] = psi[i, j] - omega*resid[i, j]/e  #Psi calculation, considering relaxation factor.
        sum_resid += np.sum(np.abs(resid))
        
        sum_resid = sum_resid/count
        iteration += 1

    print(iteration)

    return(psi)


def analytical(n_max, dx, val):

    # Analytical solution

    size = int(1.0/dx)  # internal grid size = size - 2 (boundaries)
    # n_max = 200  # Max iterations
    gamma = np.zeros((4, n_max))  # Initialisation of gamma array
    psi_anal = np.zeros((size, size))  # Initialisation of psi

    # Analytical Calculations
    for n in range(1, n_max):
        gamma[0, n] = 2.0*val/(n*pi)*(-cos(n*pi) + 1.0)
        gamma[1, n] = 2.0*val/(n*pi)*(-cos(n*pi) + 1.0)
        gamma[2, n] = 2.0*val/(n*pi)*(-cos(n*pi))
        gamma[3, n] = 2.0*val/n/pi

    # for i in range(0, size):
    #     for j in range(0, size):
    
    ii = np.arange(0, size)
    jj = np.arange(0, size)
    i, j = np.meshgrid(ii, jj) 

    # Definition of x and y
    x = i*dx
    y = j*dx

    for n in range(n_max-1, 0, -1):  # Summing terms from small to large prevents small terms 
                                        # from being neglected from the summation as they are added first.

        # Definition of hyperbolic terms
        h10 = np.sinh(n*pi)
        a1 = n*pi*x
        a2 = n*pi*(1.0-y)
        a3 = n*pi*(1.0-x)
        a4 = n*pi*y
        beta = n*pi

        # The addition of hyperbolic terms was leading to large errors due to summation of large numbers.
        # E.g. when 100000000.4 and 100000000.0 are added the relative error is small (0.4 in 10e8), but compared 
        # to the values of the amplitude of psi in this study it is significant.
        # Thus the identity sinh(a)*cosh(b)-sinh(b)*cosh(a)=sinh(a-b)was used to simplify the terms 
        # so that the summation of hyperbolic functions wasreplaced by division. 
        # The division of two big numbers (here of the form sinh(beta-an)/h10) ) doesn't result in large computational errors.

        psi_anal[i, j] += gamma[0, n]*np.sin(n*pi*y)*(np.sinh(beta-a1)/h10)
        psi_anal[i, j] += gamma[1, n]*np.sin(n*pi*x)*(np.sinh(beta-a2)/h10)
        psi_anal[i, j] += gamma[2, n]*np.sin(n*pi*y)*(np.sinh(beta-a3)/h10)
        psi_anal[i, j] += gamma[3, n]*np.sin(n*pi*x)*(np.sinh(beta-a4)/h10)

    psi_anal = np.flipud(psi_anal.T)  # Set the correct orientation of the analytical solution

    return(psi_anal)


def plots(psi, psi_anal, dx, val, plotv = [1, 1, 1]):

    size = int(1.0/dx)  # internal grid size = size - 2 (boundaries)

    # Figure 1
    # Numerical pixel plot and contours (stream lines)
    if plotv[0] == 1:

        plt.figure(1)
        plt.imshow(psi)  # Plots psi values
        plt.colorbar()
        plt.xlabel('x gridpoints')
        plt.ylabel('y gridpoints')
        plt.title('Psi potential flow (numerical) (l/min)')

        plt.figure(2)
        plt.contour(np.arange(size), np.arange(size), np.flipud(psi), 20, origin = 'lower')  # Plots contour of psi
        plt.colorbar()
        plt.xlabel('x gridpoints')
        plt.ylabel('y gridpoints')
        plt.title('Psi stream lines (l/min)')

        plt.show()

    # Array plot initialisation for quiver
    vx = np.zeros((size, size))
    vy = np.zeros((size, size))

    # Arrow calculations for vy
    s = size -1
    for i in range(size):
        for j in range(1, size-1):
            vy[i, j] = (psi[i, j+1] - psi[i, j-1])/(2*100/size)  # central dif
    for i in range(size):
            vy[i, 0] = (psi[i, 1] - psi[i, 0])/(100/size)  # forward dif
            vy[i, s] = (psi[i, s] - psi[i, s-1])/(100/size)  # backward dif

    # Arrow calculations for vx
    for i in range(1, size-1):
        for j in range(size):
            vx[i, j] = -(psi[i+1, j] - psi[i-1, j])/(2*100/size)  # central
    for j in range(size):
            vx[0, j] = -(psi[1, j] - psi[0, j])/(100/size)  # forward
            vx[s, j] = -(psi[s, j] - psi[s-1, j])/(100/size) # backward

    # Divide by 60000 so that the flow velocity is in m/s
    vy = np.flipud(-vy)/60000
    vx = np.flipud(vx)/60000

    # Figure 2
    # Quiver plot
    if plotv[1] == 1:
        color = np.sqrt(np.square(vx) + np.square(vy))  # Colormap definition
        plt.figure(3)
        # plt.quiver(vx[1:-1, 1:-1], vy[1:-1, 1:-1], color, pivot = 'tip')
        plt.quiver(vx, vy, color, pivot = 'tip')  # Quiver arrow plot
        plt.xlabel('vx')
        plt.ylabel('vy')
        plt.title('Psi arrow velocity plot (m/s)')
        plt.colorbar()
        plt.show()

    # Figure 3
    # Analytical and numerical results without the boundaries
    if plotv[2] == 1:

        plt.figure(4)
        plt.imshow(psi[1:-1, 1:-1])
        plt.xlabel('x gridpoints')
        plt.ylabel('y gridpoints')
        plt.title('Psi numerical (l/min)')
        plt.colorbar()

        plt.figure(5)
        plt.imshow(psi_anal[1:-1, 1:-1])
        plt.xlabel('x gridpoints')
        plt.ylabel('y gridpoints')
        plt.title('Psi analytical (l/min)')
        plt.colorbar()

        plt.show()

            
def error_calc(psi, psi_anal):

    # Error bonus plot and calculations (numerical vs analytical)

    # Neglect the boundaries
    psi_anal = psi_anal[1:-1, 1:-1]
    psi = psi[1:-1, 1:-1]

    error_rel = abs(psi_anal - psi)/abs(psi_anal)*100  # relative error %
    error_abs = abs(psi_anal - psi) # absolute error (difference)

    plt.figure(1)
    plt.imshow(error_abs)
    plt.colorbar()
    plt.xlabel('x gridpoints')
    plt.ylabel('y gridpoints')
    plt.title('Absolute error (analytical - numerical) (l/min)')
    plt.show()

    # The following errors are calculated for the internal cells of the geometry

    # rel_max = np.amax(error_rel)  # Maximum relative error
    # abs_max = np.amax(error_abs)  # Maximum absolute error
    rel_av = np.average(error_rel[1:-1, 1:-1])  # Average relative error
    # abs_av = np.average(error_abs)  # Average absolute error
    return('Average relative error:', rel_av)


def analytical_conv(tolerance):

    # Convergence study to determine the numbers of terms to be
    # considered in the analytical solution.
    # The results of this convergence study show that a tolerance
    # of 1.0e-4 can be achieved after 170 summation terms of the
    # analytical expression.

    size = int(1.0/dx)  # internal grid size = size - 2 (boundaries)

    # Initialisation
    analytical_old = np.zeros((size-2, size-2))
    dif_av = 1
    count = (size-2)*(size-2)
    terms = 10
    dif_av_plot = np.array([])  # Average difference for plotting
    terms_x = np.array([])  # Term count for plotting
    while dif_av > tolerance:
        
        analytical_cur = analytical(terms, dx, val)[1:-1, 1:-1]
        dif = abs(analytical_old - analytical_cur)
        dif_av = np.sum(dif)/count
        analytical_old = analytical_cur

        dif_av_plot = np.append(dif_av_plot, dif_av)
        terms_x = np.append(terms_x, terms)
        print('terms:', terms, 'dif_av:', dif_av)
        terms += 10

    plt.plot(terms_x[1:], dif_av_plot[1:])  # Analytical convergence plot
    plt.xlabel('Summation terms')
    plt.ylabel('Average difference error (l/min)')
    plt.title('Psi analytical convergence')
    plt.show()
    return(terms)


def discrepancy(reductions):

    # Convergence analysis to examine the discrepancy between the
    # numerical and the analytical results as dx decreases.
    # The discrepancy converges asymptoticaly as the grid resolution
    # increases.
    # At this point ~ 0.2 l/min (dx ~1.0e-3), the solution could 
    # be slightly oscillating but continues to converge for further 
    # mesh refinement.

    # Initialisation
    dif_av = 1
    dif_av_plot = np.array([])  # Average difference for plotting
    terms_x = np.array([])  # dx for plotting

    dx_conv = 0.1
    # while dif_av > tolerance:
    for _ in range(reductions):

        size = int(1.0/dx_conv)  # internal grid size = size - 2 (boundaries)
        count = (size-2)*(size-2)

        analytical_cur = analytical(80, dx_conv, val)[1:-1, 1:-1]
        numerical_cur = numerical(1.0e-4, dx_conv, val)[1:-1, 1:-1]
        dif = abs(numerical_cur - analytical_cur)  # Difference calculation
        dif_av = np.sum(dif)/count  # Difference average

        dif_av_plot = np.append(dif_av_plot, dif_av)
        terms_x = np.append(terms_x, dx_conv)
        print('dx:', dx_conv, 'dif_av:', dif_av)
        dx_conv = round(float(dx_conv/2), 4)

    plt.plot(terms_x, dif_av_plot)  # Analytical convergence plot
    plt.gca().invert_xaxis()
    plt.xscale('log')
    plt.xlabel('dx (m)')
    plt.ylabel('Average difference error (l/min)')
    plt.title('Analytical/Numerical discrepancy logarithmic plot')
    plt.show()
    return(dx_conv)

# Function Calling
num = numerical(1.0e-4, dx, val)
an = analytical(100, dx, val)

plots(num, an, dx, val)
error_calc(num, an)
analytical_conv(1.0e-4)
discrepancy(6)
