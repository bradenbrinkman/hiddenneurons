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
        self.N = 1000
        self.Nhid = 997
        self.Nrec = self.N-self.Nhid

        self.p = 0.2

        self.tau = 10  # time constant
        self.gain = 1 # gain of threshold-linear input-output function
        self.mu = -1.0
        self.b = self.mu*np.ones((self.N,))

        #  triplet plasticity parameters
        # self.b = .1*np.ones((self.N,)) # external input
       # self.A3plus = 6.5e-3
       # self.A2minus = 7.1e-3
       # self.tauplus = 17 # msec
       # self.tauminus = 34 # msec
       # self.taux = 101 # msec
       # self.tauy = 114 # msec
       # self.eta = 1*1e-2 # learning rate parameter