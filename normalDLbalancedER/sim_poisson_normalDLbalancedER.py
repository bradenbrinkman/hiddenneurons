# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 13:04:49 2015

@author: gabeo (original author)

Modified by Braden Brinkman.

This function simulates a network of Poisson neurons with alpha-function synapses.

Network connectivity is a Dale's Law network with Erdos-Reyni statistics and scaled by
1/sqrt(N) (balanced scaling).

"""
import os
import numpy as np
from phi import phi
from phi import phi_prime
from rates_linear_response_GaussianER import rates_ss
import params_Gaussian as params
from generate_adj_DLER import generate_adj_DLER as gen_adj
from correlation_functions import bin_pop_spiketrain, auto_covariance_pop

def sim_poisson(W, tstop, trans, dt):

    '''
    :param W: weight matrix
    :param tstop: total simulation time (including initial transient)
    :param trans: initial transient (don't record spikes until trans milliseconds)
    :param dt: Euler step
    :return:
    '''
    
    # unpackage parameters
    par = params.params_Gaussian()
    N = par.N
    tau = par.tau
    b = par.b
    gain = par.gain

    # sim variables
    Nt = int(tstop/dt)
    t = 0
    numspikes = 0
    
    maxspikes = 500*N*tstop/1000  # 500 Hz / neuron
    spktimes = np.zeros((maxspikes, 2)) # store spike times and neuron labels
    g_vec = np.zeros((Nt, N))

    # alpha function synaptic variables
    s = np.zeros((N,))
    s0 = np.zeros((N,))
    s_dummy = np.zeros((N,))
    
    a = 1. / tau
    a2 = a**2
    
    for i in range(0,Nt,1):
        
        t += dt
        
        # update each neuron's output and plasticity traces
        s_dummy += dt*(-2*a*s_dummy - a2*s)    
        s += dt*s_dummy
        
        
        
        # compute each neuron's input
        g = np.dot(W, s) + b
        g_vec[i] = g

        # decide if each neuron spikes, update synaptic output of spiking neurons
        # each neurons's rate is phi(g)
        r = phi(g, gain)

        try:
            spiket = np.random.poisson(r*dt, size=(N,))
        except:
            break

        s_dummy += spiket*a2 # a for non-unit norm alpha function (& dimensionless Wij)
        # change to s_dummy += spiket*a2 for unit norm alpha function, but then Wij has units of time
    
        ### store spike times and counts
        if t > trans:
            for j in range(N):
                if spiket[j] >= 1 and numspikes < maxspikes:
                    spktimes[numspikes, 0] = t
                    spktimes[numspikes, 1] = j
                    numspikes += 1
    
    # truncate spike time array
    spktimes = spktimes[0:numspikes, :]

    return spktimes, g_vec


if __name__ == '__main__':

    par = params.params_Gaussian()
    
    N = par.N
    p = par.p
    J0 = 0.5 # typical stdev of Gaussian weights
    weights = J0*abs(np.random.randn(N,N))/np.sqrt(p*N)
    tau = par.tau
    b = par.b
    mu = par.mu
    gain = par.gain
    seed = 1

    trans = 5. * tau  # simulation transient
    tstop = 4000. * tau + trans  # simulation time
    dt = .2 * tau  # Euler step

    W0 = gen_adj(N, p, False, seed) # generate adjacency matrix

    W = np.multiply(weights,W0) # generate weighted adjacency matrix
    
    np.savetxt('weights_normalDLbalancedER_N' + str(N) + '_p' + str(p) + '_mu' + str(mu)+ '_J0' + str(J0)+ '_seed' + str(seed) + '.txt',W,'%.6f',' ','\n')

    spktimes, g_vec = sim_poisson(W, tstop, trans, dt)
    
    np.savetxt('spktimes_normalDLbalancedER_N' + str(N) + '_p' + str(p) + '_mu' + str(mu)+ '_J0' + str(J0)+ '_seed' + str(seed) + '.txt',spktimes,'%.6f',' ','\n')
    #np.savetxt('gvec.txt',g_vec,'%.6f',' ','\n')

    # compute some statistics
    mftratesfull = rates_ss(W)
    np.savetxt('mftratesfull_normalDLbalancedER_N' + str(N) + '_p' + str(p) + '_mu' + str(mu)+ '_J0' + str(J0) + '_seed' + str(seed) + '.txt',mftratesfull,'%.6f',' ','\n')
    # mftrateshidden = rates_ss(W[]) # need to sweep over different subsets of hidden units
    ind_include = range(N)  # indices of E neurons
    spk_Epop = bin_pop_spiketrain(spktimes, dt, 1, tstop, trans, ind_include)
    dt_ccg = 1.  # ms
    lags = np.arange(-10.*tau, 10.*tau, dt_ccg)
    numspikes = len(spktimes)
   # pop_2point = auto_covariance_pop(spktimes, ind_include, numspikes, dt, lags, tau, tstop, trans)



