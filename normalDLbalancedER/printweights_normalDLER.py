# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 13:04:49 2015

@author: gabeo (original author)

Modified by Braden Brinkman.

This function produces and prints to file a weighted adjacency matrix with Erdos-Reyni 
connectivity and Gaussian weights, with Dale's law imposed.

"""
import os
import numpy as np
from phi import phi
from phi import phi_prime
from rates_linear_response_GaussianER import rates_ss
import params_Gaussian as params
from generate_adj_DLER import generate_adj_DLER as gen_adj
from correlation_functions import bin_pop_spiketrain, auto_covariance_pop


if __name__ == '__main__':

    par = params.params_Gaussian()
    
    N = par.N
    p = par.p
    J0 = 1.0  # typical stdev of Gaussian weights
    weights = J0*abs(np.random.randn(N,N))
    tau = par.tau
    b = par.b
    mu = par.mu
    gain = par.gain
    seed = 10

    W0 = gen_adj(N, p, False, seed) # generate adjacency matrix

    W = np.multiply(weights,W0) # generate weighted adjacency matrix
    
    np.savetxt('baseweights_normalDLER_N' + str(N) + '_p' + str(p) + '_mu' + str(mu)+ '_J0' + str(J0)+ '_seed' + str(seed) + '.txt',W,'%.6f',' ','\n')

#    r = rates_ss(W)
    
#    np.savetxt('MFTfull_normalDLER_N' + str(N) + '_p' + str(p) + '_mu' + str(mu)+ '_J0' + str(J0) + '_seed' + str(seed) + '.txt',r,'%.6f',' ','\n')

