# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 09:04:07 2015

@author: gabeo (Original author)

Modified by Braden Brinkman.

This function defines a class to hold parameters for an Erdos-Reyni network with Gaussian weights
"""

import numpy as np

class params_Gaussian:
    def __init__(self):
        self.N = 1000 # Total number of neurons
        self.Nhid = 997 # Number of hidden neurons
        self.Nrec = self.N-self.Nhid # Number of observed neurons

        self.p = 0.2 # fraction of connectivity in the network.

        self.tau = 10  # time constant
        self.gain = 1 # gain of threshold-linear input-output function
        self.mu = -1.0 # sets baseline firing rate of the network
        self.b = self.mu*np.ones((self.N,)) # array with each element equal to mu.