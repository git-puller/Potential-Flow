''' 
Function Documentation
Sokratis Anagnostopoulos - Potential Flow
'''    


def numerical(tolerance, dx, val):

    '''

    Numerical solution of potential flow problem
    
    Parameters
    ----------
    tolerance : Convergence tolerance (1.0e-4)
    dx : Element length (dx = dy)
    val : Potential at inlet in l/min

    Returns
    ----------
    psi : Final state of potential psi.
    '''


def analytical(n_max, dx, val):

    '''

    Analytical solution of potential flow problem
    
    Parameters
    ----------
    n_max : Number of summation terms considered in the series
    dx : Element length (dx = dy)
    val : Potential at inlet in l/min

    Returns
    ----------
    psi_anal : Final state of potential psi analytical.
    '''


def plots(psi, psi_anal, dx, val):

    '''

    Plots of the following:

    Figure 1
    Numerical pixel plot and contours (stream lines) of psi numerical

    Figure 2
    Quiver plot of stream lines (numerical)

    Figure 3
    Analytical and numerical results without the boundaries
    
    Parameters
    ----------
    psi : Numerical solution of psi
    psi_anal: Analitical solution of psi
    dx : Element length (dx = dy)
    val : Potential at inlet in l/min
    plotv : Plots' array (figures 1-3 are plotted when their corresponding value is 1)

    Returns
    ----------
    
    '''


def error_calc(psi, psi_anal):

    '''

    Plots of numerical vs analitical results and calculations
    of max and average errors (relative and absolute).
    Only the absolute error is ploted but the rest can be easily
    plotted as well if un-commented.

    Parameters
    ----------
    psi : Numerical solution of psi
    psi_anal: Analitical solution of psi

    Returns
    ----------
    rel_av : Average relative error between numerical and analytical.
    '''


def analytical_conv(tolerance):

    '''

    Convergence study to determine the numbers of terms to be
    considered in the analytical solution.
    The results of this convergence study show that a tolerance
    of 1.0e-4 can be achieved after 170 summation terms of the
    analytical expression.

    Parameters
    ----------
    tolerance : Convergence tolerance (1.0e-4)

    Returns
    ----------
    terms : Number of terms required in the analytical solution 
    to achieve the desired tolerance
    '''


def discrepancy():

    '''

    Convergence analysis to examine the discrepancy between the
    numerical and the analytical results as dx decreases.
    The discrepancy converges asymptoticaly as the grid resolution
    increases.
    At this point ~ 0.2 l/min (dx ~1.0e-3), the solution could 
    be slightly oscillating but continues to converge for further 
    mesh refinement.

    Parameters
    ----------
    reductions : Defines how many times will dx be reduced by a factor of 2
    in order for the convergence plot to be created.

    Returns
    ----------
    dx_conv = Smallest dx used.
    '''
