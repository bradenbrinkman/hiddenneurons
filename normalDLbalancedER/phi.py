# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 13:08:02 2015

@author: gabeo (Original author)

Modified by Braden Brinkman.

This script defines the activation function for Poisson model, along with its derivatives.
"""

''' here, exponential nonlinearity '''
def phi(g,gain):

    import numpy as np

    g_calc = np.exp(g*1)

    r_out = g_calc

    return r_out

def phi_prime(g,gain):

    '''
    first derivative of phi wrt input
    '''

    import numpy as np

    phi_pr = np.exp(g*1)

    return phi_pr

def phi_prime2(g,gain):

    '''
    second derivative of phi wrt input
    '''

    import numpy as np
    
    phi_pr2 = np.exp(g*1)

    return phi_pr2