    '''

    A set of functions used by sim_poisson_hiddenlinearrates_normalsignedbalancedER.py,
    sim_Jefft_normalsignedbalancedER.py, and sim_Jeffw_normalsignedbalancedER.py

    '''

import math
import numpy as np
from phi import phi, phi_prime, phi_prime2
import params_hiddenlinearrates_Gaussian as params_Gaussian;

reload(params_Gaussian)


def rates_ss(W):  # set inputs here

    '''
    compute steady-state mean field rates through Euler step
    :param W: weight matrix
    :return: steady-state rates with transfer functions and cellular/synaptic parameters defined in params.py and phi.py
    '''

    par = params_Gaussian.params_Gaussian()
    b = par.b
    gain = par.gain
    tau = par.tau
    N = par.N

    dt = .02 * tau
    Tmax = int(50 * tau / dt)
    a = 1. / tau
    a2 = a ** 2

    r = np.zeros(N)
    s_dummy = np.zeros(N)
    s = np.zeros(N)

    r_vec = np.zeros((N, Tmax))
    for i in range(Tmax):
        s_dummy += dt * (-2 * a * s_dummy - a2 * s) + r * a2 * dt # r * a2 for unit norm alpha function
        s += dt * s_dummy

        g = W.dot(s) + b
        r = phi(g, gain)
        r_vec[:, i] = r

    return r
    
def hidrates_ss(W):  # set inputs here

    '''
    compute steady-state mean field rates through Euler step
    :param W: weight matrix
    :return: steady-state rates with transfer functions and cellular/synaptic parameters defined in params.py and phi.py
    '''

    par = params_Gaussian.params_Gaussian()
    b = par.b
    gain = par.gain
    tau = par.tau
    N = par.N
    Nhid = par.Nhid
    
    bhid = b[0:Nhid]

    dt = .02 * tau
    Tmax = int(50 * tau / dt)
    a = 1. / tau
    a2 = a ** 2

    r = np.zeros(Nhid)
    s_dummy = np.zeros(Nhid)
    s = np.zeros(Nhid)

    r_vec = np.zeros((Nhid, Tmax))
    for i in range(Tmax):
        s_dummy += dt * (-2 * a * s_dummy - a2 * s) + r * a2 * dt # r * a2 for unit norm alpha function
        s += dt * s_dummy

        g = W.dot(s) + bhid
        r = phi(g, gain)
        r_vec[:, i] = r

    return r


def rates_1loop(W):
    """
    inputs: weight matrix
    calculate one-loop correction for steady-state firing rates in fluctuation expansion
    """

    import params;
    reload(params)
    from phi import phi_prime
    from phi import phi_prime2

    par = params.params()
    N = par.N
    gain = par.gain
    b = par.b

    phi_r = rates_ss(W)
    Tmax = 100
    dt_ccg = 1.
    wmax = 1. / dt_ccg
    dw = 1. / Tmax

    w = np.arange(-wmax, wmax, dw) * math.pi
    dw = dw * math.pi
    Nw = w.size

    g0 = np.dot(W, phi_r) + b
    phi_1 = phi_prime(g0, gain)
    phi_1 = np.diag(phi_1)

    phi_2 = phi_prime2(g0, gain)
    Fbarsum = np.zeros((N, Nw), dtype=complex)

    for o in range(Nw):  # first compute Fbar over dummy frequency
        Fbar1 = np.dot(g_fun(w[o]) * W, linear_response_fun(w[o], np.dot(phi_1, W), phi_r))
        Fbarsum[:, o] = np.dot(Fbar1 * Fbar1.conj(), phi_r)  # sum over first inner vertex

    Fbarsum_int = np.sum(Fbarsum, axis=1) * dw  # integrate over dummy frequency

    F1 = linear_response_fun(0., np.dot(phi_1, W), phi_r)

    r_1loop = np.dot(F1, .5 * phi_2 * Fbarsum_int) / ((2 * math.pi) ** 1)  # sum over second inner vertex
    return r_1loop


def g_fun(w):
    import numpy as np
    import params;
    reload(params)
    par = params.params()
    tau = par.tau

    taud = 0.

    g = np.exp(-1j * w * taud) / ((1 + 1j * w * tau) ** 2)  # alpha function

    return g


def linear_response_fun(w, W, phi_r):
	
    import numpy as np
    par = params_Gaussian.params_Gaussian()
    N = par.N
    
    s = np.multiply(phi_r,W)

    Gamma = g_fun(w) * s 
    Delta1 = np.linalg.inv(np.eye(N) - Gamma)
    Delta = np.multiply(Delta1,phi_r)

    return Delta
    
def linear_response_fun_hid(w, Whh, phi_h):
	
    import numpy as np
    par = params_Gaussian.params_Gaussian()
    N = par.N
    Nhid = par.Nhid
    
    s = np.multiply(phi_h,Whh)

    Gamma = g_fun(w) * s
    Delta1 = np.linalg.inv(np.eye(Nhid) - Gamma)
#    Delta = np.linalg.inv(np.multiply(phi_h,np.eye(Nhid)) - g_fun(w) * Whh)
    Delta = np.multiply(Delta1,phi_h)

    return Delta
    
def Jeff_w(w, Wrr, Wrh, Whr, Whh, phi_h):
	
    import numpy as np
    par = params_Gaussian.params_Gaussian()
    N = par.N
    Nhid = par.Nhid
    
    Jdirect = g_fun(w) * Wrr
    
    Delta = linear_response_fun_hid(w, Whh, phi_h)
    
    dJ = g_fun(w)*g_fun(w)*Wrh.dot(Delta.dot(Whr))
    
    Jeff = Jdirect + dJ

    return Jeff
    
def Jeff_t(Wrr, Wrh, Whr, Whh, phi_h):
	
    import numpy as np
    par = params_Gaussian.params_Gaussian()
    N = par.N
    Nhid = par.Nhid
    Nrec = par.Nrec
    tau = par.tau
    
    Tmax = 20. * tau
    dt = 0.1*tau # doesn't have to be same dt as simulation dt.
    wmax = 10. / tau 
    dw = 1. / Tmax

    t = np.arange(0,Tmax,dt)
    w = np.arange(-wmax, wmax, dw) * math.pi
    dw = dw * math.pi
    Nw = w.size
    Nt = t.size

    Jsum = np.zeros((Nrec, Nrec, Nw), dtype=complex)
    Expiwt = np.zeros((Nw,Nt), dtype=complex)

    for o in range(Nw):  # first compute Fbar over dummy frequency
        Jsum[:, :, o] = Jeff_w(w[o], Wrr, Wrh, Whr, Whh, phi_h)
        for tt in range(Nt):
            Expiwt[o,tt] = np.exp(1j*w[o]*t[tt])

 #   Jeff_t = np.zeros((Nrec, Nrec, Nt), dtype=complex)

    Jeff_t = np.dot(Jsum,Expiwt) * dw / 2 / math.pi # integrate over frequency
    
    Jeff_t = Jeff_t.real

    return Jeff_t
    
def linear_response_hiddenrate_fun(w, Whh, Whr, phi_h):
	
    import numpy as np
    par = params_Gaussian.params_Gaussian()
    N = par.N
    Nhid = par.Nhid

    Delta = linear_response_fun_hid(w, Whh, phi_h)
    filter = g_fun(w) * Delta.dot(Whr)

    return filter
    
def empavg_single(spktimes,neuronind,dt,tstop,trans):
    
    import numpy as np
    
    # find spike times
    ind_t = np.where(spktimes[:,1]==neuronind)
  #  t_spk = spktimes[ind_t,0]-trans
  #  t_spk = np.transpose(t_spk)
    Nspk = np.size(ind_t)
    
    return Nspk
    
    
def empavg_pop(spktimes,N,dt,tstop,trans,ind_include):
    
    import numpy as np
    
    Nspk_pop = np.zeros((N,))    
    
    for i in ind_include:
        Nspk_pop[i] = empavg_single(spktimes, i, dt, tstop, trans)
    
    return Nspk_pop
    
def spiketrain_FT_hiddenrate(spktimes,neuronind,dt,dt_ccg,tstop,trans):

    import numpy as np
    par = params_Gaussian.params_Gaussian()
    N = par.N
    Nhid = par.Nhid
    
    spikes = bin_spiketrain(spktimes,neuronind,dt,dt_ccg,tstop,trans)

    spkfft_temp = np.fft.fft(spikes,axis=1)
    
    spkfft = np.fft.fftshift(spkfft_temp,axis=1)
    
    freqs = np.fft.fftfreq(spikes.size(1,),d=dt)

    return spkfft, freqs
    
def rates_timeavgapprox_linearhiddenrates(Whh,Whr,phi_h,spktimes,neuronind,dt,dt_ccg,tstop,trans): 

    import numpy as np
    from rates_linear_response_hiddenlinearrates_GaussianER import rates_ss
 #   from correlation_functions import bin_spiketrain
    from correlation_functions_hiddenlinearrates import bin_pop_array_spiketrain
    
    spikes = bin_pop_array_spiketrain(spktimes,dt,dt_ccg,tstop,trans,neuronind)
    
    spkfft = np.fft.rfft(spikes,axis=0)
    
    #freqs = np.fft.fftfreq(spikes.size(1,),d=dt)
    
    # Need only zero freq component for approx
    
    filter = linear_response_hiddenrate_fun(0.0, Whh, Whr, phi_h)
#        filteredspikes_freq[count] = np.multiply(filter,temp)
    
    spk0 = spkfft[0,:]
    L0 = filter.dot(spk0.real)
    
    
    vh = hidrates_ss(Whh)
    
    Nhavg = vh + 0.5*L0
    
    return spk0.real
    
#    return Nhavg
    
    
#def rates_linearhiddenrates(Whh,Whr,phi_h,spktimes,neuronind,dt,dt_ccg,tstop,trans): 

#   import numpy as np
    
#    spkfft, freqs = spiketrain_FT_hiddenrate(spktimes,neuronind,dt,dt_ccg,tstop,trans)
    
#    count = 0
#    filter = np.array(freqs.size)
#    for w in freqs:
#        filter = linear_response_hiddenrate_fun(w, Whh, Whr, phi_h)
#        temp = spkfft[:,count]
#        filteredspikes_freq[count] = np.multiply(filter,temp)
#        count += count
    
#    # inverse fourier transform
    
#    filteredspikes_temp = np.fft.ifftshift(filteredspikes_freq,axis=1)
#    
#    filteredspikes = np.fft.ifft(filteredspikes_temp,axis=1)
    
#    return filteredspikes
       


def linear_response_1loop(w, W, phi_r):
    '''
    calculate one-loop correction to the propagator around mean-field theory
    :param w: frequency
    :param W: weight matrix, weighted by postsynaptic gain
    :param phi_r: firing rates
    :return: propagator matrix
    '''

    par = params.params()

    b = par.b
    gain = par.gain
    N = par.N

    Tmax = 100
    dt_ccg = 1
    wmax = 1. / dt_ccg
    dw = 1. / Tmax

    w_calc = np.arange(-wmax, wmax, dw) * math.pi
    dw *= math.pi
    Nw = w_calc.size

    g0 = np.dot(W, phi_r) + b
    phi_1 = phi_prime(g0, gain)
    phi_1_diag = np.diag(phi_1)

    phi_2 = phi_prime2(g0, gain)
    phi_2_diag = np.diag(phi_2)

    F1 = linear_response_fun(w, np.dot(phi_1_diag, W), phi_r)
    Fbar = np.dot(g_fun(w) * W, F1)

    Fbar_int = np.zeros((N, N), dtype='complex128')

    for o in range(Nw):
        Fbar1 = np.dot(g_fun(w_calc[o]) * W, linear_response_fun(w_calc[o], np.dot(phi_1_diag, W), phi_r))
        Fbar2 = np.dot(g_fun(w - w_calc[0]) * W, linear_response_fun(w - w_calc[o], np.dot(phi_1_diag, W), phi_r))
        Fbar_int += np.dot(Fbar1 * Fbar2, np.dot(phi_1_diag, Fbar)) * dw

    linear_response_1loop = np.dot(np.dot(F1, phi_2_diag / 2.), Fbar_int) / (2 * math.pi ** 1)

    return linear_response_1loop
