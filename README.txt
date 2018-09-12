README for code accompanying "Predicting How and When Hidden Neurons Skew Measured Synaptic Interactions.”

===========================================================================
=================MATLAB SCRIPTS FOR PLOTTING FIGURES======================
===========================================================================

MATLAB scripts were originally written for MATLAB version 2013b. The author makes no guarantee these scripts work for any older or newer version.

The following scripts necessary reproduce figures in the text. The necessary data to be plotted may be downloaded from the following links:

Jdata: https://www.dropbox.com/s/bhdrw0ht8gf0z7v/Jdata.zip?dl=0
Fig6: https://www.dropbox.com/s/6vzipad37t7d1qz/Fig6.zip?dl=0
SFigs: https://www.dropbox.com/s/w3q0e2xmml12zpu/SFigs.zip?dl=0

The extracted folders should be placed in the same folder as the MATLAB scripts. The functions of each script are:

1. ffwiplotter.m reproduces the inset plots in Figure 3.
2. ffwiplotter_4n.m reproduces the inset plots in Figure 4.
3. Janalysissubplots.m calls data in the subfolder Jdata to reproduce Figure 5.
4. Jefftplotter.m calls data in the subfolder Fig6 to reproduce Figure 6.
5. Janalysis_ERplots_dJ_sigmoidnonlin.m calls Jdata to plot Figures 7 and 8. 
6. Janalysis_ERplots_dJ_N100.m calls data in Jdata to plot Figure 9.
7. hiddenlinearrates_validationcheckplots.m call data in the subfolder SFigs to produce supplementary figures S1-3.

The scripts Janalysis_ERplots_dJ_sigmoidnonlin.m and hiddenlinearrates_validationcheckplots.m require the user to change some parameters and rerun to produce particular figures. These changes are noted in the comments of these scripts.

===========================================================================
=======================PYTHON SIMULATION CODE============================
===========================================================================

Python code for performing network simulations was originally written by Gabe Ocker to accompany the paper “Linking structure and activity in nonlinear spiking networks,” PLOS Computational Biology https://doi.org/10.1371/journal.pcbi.1005583. It has been modified by the lead author of this work (Braden Brinkman). The modified code is included in this repository.

There are four folders, each with a similar set of functions. The primary difference between the folders is the statistical structure of the networks used in each simulation. The four folders are named after the corresponding structure of the adjacency matrices used in the simulation:

normalDLbalancedER: Erdos-Reyni connectivity. Weights are drawn from a Gaussian distribution, but signs adjusted such that all of a neuron’s weights are the same sign (Dale’s Law). Weights are scaled by 1/sqrt(N), where N is the number of neurons, known as balanced scaling.

normalDLclassicalER: Erdos-Reyni connectivity. Weights are drawn from a Gaussian distribution, but signs adjusted such that all of a neuron’s weights are the same sign (Dale’s Law). Weights are scaled by 1/N, where N is the number of neurons, known as classical scaling.

normalsignedbalancedER: Erdos-Reyni connectivity. Weights are drawn from a Gaussian distribution. Weights are scaled by 1/sqrt(N), where N is the number of neurons, known as balanced scaling.

normalsignedclassicalER: Erdos-Reyni connectivity. Weights are drawn from a Gaussian distribution. Weights are scaled by 1/N, where N is the number of neurons, known as classical scaling.

The common functions in each folder are:

correlation_functions_hiddenlinearrates.py: defines a set of functions used to calculate the correlation functions of a network, both using the output of the network simulations and calculation techniques of Ocker et al. (2017).

generate_adj_X.py: Generates the appropriate adjacency matrix corresponding to statistics X (e.g., X = DLER for Dale’s Law Erdos-Reyni networks).

params_hiddenlinearrates_Gaussian.py: set the parameters for the network simulations.

phi.py: Sets the nonlinearity phi(x).

printweights_X.py: generates and prints to file an adjacency matrix with statistics X (e.g., X = DLER for Dale’s Law Erdos-Reyni networks).

rates_linear_response_hiddenlinearrates_GaussianER.py: defines a set of functions used to calculate the firing rates of a network, both using the output of the network simulations and calculation techniques of Ocker et al. (2017).

sim_poisson_X.py: Simulates the nonlinear Hawkes model with network statistics X. (e.g., X = normalDLbalancedER for a network with Dale’s Law, Erdos-Reyni connectivity, normal (Gaussian) weights, and balanced scaling).

In addition to the above functions, the folder normalsignedbalancedER contains the additional following functions:

sim_poisson_hiddenlinearrates_normalsignedbalancedER.py: unlike the normal sim_poisson_X.py, explicitly partitions the network into recorded and hidden neurons, and estimates the firing rates of both the full network and a network containing only the hidden neurons.

sim_Jeffw_normalsignedbalancedER.py: calculates the effective synaptic interaction in the frequency domain.

sim_Jefft_normalsignedbalancedER.py: calculates the effective synaptic interaction in the time domain.

