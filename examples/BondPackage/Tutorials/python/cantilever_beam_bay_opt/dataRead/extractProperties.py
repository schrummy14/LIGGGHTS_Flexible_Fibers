
import numpy as np
from scipy.optimize import least_squares

def getFun(t,z):
    fun = lambda b: b[0]*np.exp(-b[1]*t)*np.cos(2*np.pi*b[2]*t)
    minFun = lambda b: fun(b) - z
    return fun, minFun

def diff(x):
    dx = np.zeros(len(x)-1)
    dx = x[1:]-x[0:-1]
    return dx

def findFreq(t,z,doGrad = True, count=0):
    if count > 1:
        return 0
    # Do the same as in matlab...
    z_org = np.copy(z)
    if (doGrad):
        z = np.gradient(z,t)
    z1 = z[0:-1]
    z2 = z[1:]
    ids = np.argwhere(z1*z2 < 0)
    t_id_1 = t[ids]
    t_id_2 = t[ids+1]  
    z_id_1 = z[ids]
    z_id_2 = z[ids+1]
    
    m = (z_id_2 - z_id_1) / (t_id_2 - t_id_1)
    t_zero = t_id_1 - z_id_1/m
    std_diff = np.std(diff(t_zero))
    if (std_diff > 0.001):
        return findFreq(t,z_org,False,count+1)
    T = 2.0*np.mean(diff(t_zero))

    return 1.0/T

def getDamp(t,z):
    fun, minFun = getFun(t,z)
    x0 = [1.0, 0.0, 100.0]
    lb = 0.0
    ub = np.inf
    res = least_squares(minFun, x0, bounds = (lb, ub))
    return res.x[-2]

def getProps(data):
    t = data[:,0]
    z = data[:,2]
    deflection_distance = 20.0e-3
    deflection_speed = 1.0
    wait_time = 25.0e-3
    t_release = 0.046415008 # deflection_distance/deflection_speed + wait_time
    
    z = z[t>=t_release]
    t = t[t>=t_release]
    t = t - t[0]

    freq = findFreq(t,z)

    z = z/z[0]
    t = freq*t

    damp = getDamp(t,z)
    
    if damp is None:
        damp = 0.0
    if freq is None:
        freq = 0.0

    return freq, damp
