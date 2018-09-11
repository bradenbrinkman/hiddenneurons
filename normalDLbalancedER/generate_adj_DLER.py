# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 13:16:21 2015

@author: gabeo (original author)

Modified by Braden Brinkman.

This function generates a non-symmetric, signed Erdos-Renyi adjacency matrix 
with Dale's law imposed.

Inputs are: 
N, the number of neurons
p, the expected fraction of non-zero connections
aut, a boolean that indicates whether to include autapses (non-zero self connections). Defaults to False if not provided.
rngseed, the seed for the random number generator. Defaults to 1 if not provided.

Outputs are:
W0, the weighted adjacency matrix.

dependency: numpy
"""
import numpy as np

def generate_adj_DLER(N,p,aut=False,rngseed=1):

    rng = np.random.RandomState(rngseed)
    
    W0 = np.random.rand(N,N)
    
    mask = W0<p
    
    W0[mask] = 1
    W0[np.logical_not(mask)] = 0
    
    if aut:
        np.fill_diagonal(W0,1) # autapses to allow for refractoriness or burstiness
    else:
        np.fill_diagonal(W0,0) # disallow autapses; set to -1 for refractory-like effects
        
    signs = 2*np.random.randint(0,2,size=(1,N))-1 #[(-1)**np.random.randint(0,1) for i in range(N)] 
    
    W0 = W0*signs #imposes a random sign to each column, enforcing Dale's law 
    
    
    return W0