# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 13:04:49 2015

@author: gabeo (original author)

Modified by Braden Brinkman.

This function simulates a network of Poisson neurons with alpha-function synapses.

Unlike the other simulation functions, which simulate the full network, this function
explicitly partitions the network into recorded and hidden neurons, and estimates the 
firing rates of both the full network and a network containing only the hidden neurons.
The rates are estimated both by simulating the networks and also using mean field 
approximations so that the two estimates can be compared.

"""
import os
import numpy as np
from phi import phi
from phi import phi_prime
from rates_linear_response_hiddenlinearrates_GaussianER import rates_ss
from rates_linear_response_hiddenlinearrates_GaussianER import hidrates_ss
from rates_linear_response_hiddenlinearrates_GaussianER import rates_timeavgapprox_linearhiddenrates
from rates_linear_response_hiddenlinearrates_GaussianER import linear_response_hiddenrate_fun
from rates_linear_response_hiddenlinearrates_GaussianER import empavg_pop
import params_hiddenlinearrates_Gaussian as params
from generate_adj_signedER import generate_adj_signedER as gen_adj
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
        # keep as s_dummy += spiket*a2 for unit norm alpha function, but then Wij has units of time
    
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
    
    
def hiddenspkfilter(W,recspikes,t,dt):

    '''
    :param W: weight matrix
    :param tstop: total simulation time (including initial transient)
    :param trans: initial transient (don't record spikes until trans milliseconds)
    :param dt: Euler step
    :return:
    '''


if __name__ == '__main__':

    par = params.params_Gaussian()
    
    N = par.N
    Nhid = par.Nhid
    p = par.p
    J0 = 1.0 # typical stdev of Gaussian weights
    weights = J0*np.random.randn(N,N)/np.sqrt(p*N)
    tau = par.tau
    b = par.b
    mu = par.mu
    gain = par.gain
    seed = 1

    trans = 5. * tau  # simulation transient
    tstop = 4000. * tau + trans  # simulation time
#    dt = .2 * tau  # Euler step # original time step in Gabe's code, seems too long
    dt = 0.01 * tau # shorter time step to get convergence to MFT results
#    dt = .005 * tau  # Euler step -- this is better, but takes a while to run

    W0 = gen_adj(N, p, False, seed) # generate full adjacency matrix

    W = np.multiply(weights,W0) # generate full weighted adjacency matrix
    
#    rperm = np.random.permutation(N)
#    hidrange = rperm[0:Nhid]
#    recrange = rperm[Nhid:]

    hidrange = np.arange(0,Nhid)
    recrange = np.arange(Nhid,N)
    
    Wrr = W[recrange[:,None],recrange[None,:]]
    Wrh = W[recrange[:,None],hidrange[None,:]]
    Whr = W[hidrange[:,None],recrange[None,:]]
    Whh = W[hidrange[:,None],hidrange[None,:]]
    
    mftrates_hid = hidrates_ss(Whh)
    Lhr0 = linear_response_hiddenrate_fun(0.0, Whh, Whr, mftrates_hid)
    Lhr0 = Lhr0.real
    
    numtrials = 5
    dt_ccg = 1.  # ms
    empavg = np.zeros((N,numtrials))
    
    for i in range(numtrials):
        spktimes, g_vec = sim_poisson(W, tstop, trans, dt)
        empavg[:,i] = empavg_pop(spktimes,N,dt,tstop,trans,range(N))
        # np.savetxt('spktimes_hiddenlinearrates_normalsignedbalancedER_N' + str(N) + '_p' + str(p) + '_mu' + str(mu)+ '_J0' + str(J0) + '_seed' + str(seed) + '_trial' + str(i) + '.txt',spktimes,'%.6f',' ','\n')

        
    empavgavg = np.mean(empavg,axis=1) #average over trials
#    np.savetxt('test1.txt',empavgavg,'%.6f','\n')

    empavgavg = empavgavg/(tstop-trans) #average spike rates
#    np.savetxt('test2.txt',empavgavg,'%.6f','\n')
    
    empavgavg_rec = empavgavg[recrange]
    empavgavg_hid = empavgavg[hidrange]
    
    
    Nhavgapprox = mftrates_hid + Lhr0.dot(empavgavg_rec)
    
    mftrates_full = rates_ss(W)
    
    np.savetxt('Avgrates_empirical_hiddenlinearrates_normalsignedbalancedER_N' + str(N) + '_Nhid' + str(Nhid) + '_p' + str(p) + '_mu' + str(mu)+ '_J0' + str(J0) + '_seed' + str(seed) + '.txt',empavgavg_hid,'%.6f',' ','\n')
    np.savetxt('Avgrates_approx_hiddenlinearrates_normalsignedbalancedER_N' + str(N) + '_Nhid' + str(Nhid) + '_p' + str(p) + '_mu' + str(mu)+ '_J0' + str(J0) + '_seed' + str(seed) + '.txt',Nhavgapprox,'%.6f',' ','\n')
    np.savetxt('Avgrates_MFT_hiddenlinearrates_normalsignedbalancedER_N' + str(N) + '_Nhid' + str(Nhid) + '_p' + str(p) + '_mu' + str(mu)+ '_J0' + str(J0) + '_seed' + str(seed) + '.txt',mftrates_hid,'%.6f',' ','\n')

    np.savetxt('Avgrates_MFTfull_hiddenlinearrates_normalsignedbalancedER_N' + str(N) + '_Nhid' + str(Nhid) + '_p' + str(p) + '_mu' + str(mu)+ '_J0' + str(J0) + '_seed' + str(seed) + '.txt',mftrates_full,'%.6f',' ','\n')
   
    
    #spktimesonly = spktimes[:,0]
    #neuronind = spktimes[:,1]
    #neuronindrange = recrange
    

    #Nhavg = rates_timeavgapprox_linearhiddenrates(Whh,Whr,mftrates_hid,spktimes,neuronindrange,dt,dt_ccg,tstop,trans)
    
#    np.savetxt('Avgrates_hiddenlinearrates_normalsignedbalancedER_N' + str(N) + '_Nhid' + str(Nhid) + '_p' + str(p) + '_mu' + str(mu)+ '_J0' + str(J0) + '_seed' + str(seed) + '.txt',Nhavg,'%.6f',' ','\n')
#    np.savetxt('MFTrates_hiddenlinearrates_normalsignedbalancedER_N' + str(N) + '_Nhid' + str(Nhid) + '_p' + str(p) + '_mu' + str(mu)+ '_J0' + str(J0) + '_seed' + str(seed) + '.txt',mftrates_hid,'%.6f',' ','\n')

    
    
#    np.savetxt('weights_hiddenlinearrates_normalsignedbalancedER_N' + str(N) + '_p' + str(p) + '_mu' + str(mu)+ '_J0' + str(J0) + '_seed' + str(seed) + '.txt',W,'%.6f',' ','\n')
    
#    np.savetxt('spktimes_hiddenlinearrates_normalsignedbalancedER_N' + str(N) + '_p' + str(p) + '_mu' + str(mu)+ '_J0' + str(J0) + '_seed' + str(seed) + '.txt',spktimes,'%.6f',' ','\n')


#    np.savetxt('mftratesfull_normalsignedbalancedER_N' + str(N) + '_p' + str(p) + '_mu' + str(mu)+ '_J0' + str(J0) + '_seed' + str(seed) + '.txt',mftratesfull,'%.6f',' ','\n')
#    ind_include = range(N)  # indices of E neurons
#    spk_Epop = bin_pop_spiketrain(spktimes, dt, 1, tstop, trans, ind_include)

  #  lags = np.arange(-10.*tau, 10.*tau, dt_ccg)
  #  numspikes = len(spktimes)
   # pop_2point = auto_covariance_pop(spktimes, ind_include, numspikes, dt, lags, tau, tstop, trans)



