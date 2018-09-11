# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 13:04:49 2015

@author: gabeo (original author)

Modified by Braden Brinkman.

This function computes the effective coupling Jeff(t) in the time domain and prints it to
file.

"""
import os
import math
import numpy as np
from phi import phi
from phi import phi_prime
from rates_linear_response_hiddenlinearrates_GaussianER import rates_ss
from rates_linear_response_hiddenlinearrates_GaussianER import hidrates_ss
from rates_linear_response_hiddenlinearrates_GaussianER import rates_timeavgapprox_linearhiddenrates
from rates_linear_response_hiddenlinearrates_GaussianER import linear_response_hiddenrate_fun
from rates_linear_response_hiddenlinearrates_GaussianER import empavg_pop
from rates_linear_response_hiddenlinearrates_GaussianER import Jeff_w
from rates_linear_response_hiddenlinearrates_GaussianER import Jeff_t
import params_hiddenlinearrates_Gaussian as params
from generate_adj_signedER import generate_adj_signedER as gen_adj
from correlation_functions import bin_pop_spiketrain, auto_covariance_pop

if __name__ == '__main__':

    par = params.params_Gaussian()
    
    N = par.N
    Nhid = par.Nhid
    Nrec = par.Nrec
    p = par.p
    J0 = 1.0 # typical stdev of Gaussian weights
    weights = J0*np.random.randn(N,N)/np.sqrt(p*N)
    tau = par.tau
    b = par.b
    mu = par.mu
    gain = par.gain
    seed = 9

    trans = 5. * tau  # simulation transient
    tstop = 4000. * tau + trans  # simulation time
#    dt = .2 * tau  # Euler step # original time step in Gabe's code, seems too long
    dt = 0.01 * tau # shorter time step to get convergence to MFT results
#    dt = .005 * tau  # Euler step -- this is better, but takes a while to run

    W0 = gen_adj(N, p, False, seed) # generate full adjacency matrix

    W = np.multiply(weights,W0) # generate full weighted adjacency matrix

    hidrange = np.arange(0,Nhid)
    recrange = np.arange(Nhid,N)
    
    Wrr = W[recrange[:,None],recrange[None,:]]
    Wrh = W[recrange[:,None],hidrange[None,:]]
    Whr = W[hidrange[:,None],recrange[None,:]]
    Whh = W[hidrange[:,None],hidrange[None,:]]
    
    mftrates_hid = hidrates_ss(Whh)

    Jefft = Jeff_t(Wrr, Wrh, Whr, Whh, mftrates_hid)
    Jefft = Jefft.real
        
    for ii in range(Nrec):
        for jj in range(Nrec):
            np.savetxt('Jefft_r1' + str(ii) + '_r2' + str(jj) + '_normalsignedbalancedER_N' + str(N) + '_Nhid' + str(Nhid) + '_p' + str(p) + '_mu' + str(mu)+ '_J0' + str(J0) + '_seed' + str(seed) + '.txt',Jefft[ii,jj,:],'%.6f','\n')

    np.savetxt('Wrr_normalsignedbalancedER_N' + str(N) + '_Nhid' + str(Nhid) + '_p' + str(p) + '_mu' + str(mu)+ '_J0' + str(J0) + '_seed' + str(seed) + '.txt',Wrr,'%.12f',' ','\n')




