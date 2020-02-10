
import numpy as np
import scipy as sp
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt

from glob import glob
from makeDOE import readDOEparams
from scipy.optimize import least_squares
from statsmodels.stats.anova import anova_lm

def stepwiselm(X,y,degree="quadratic"):
    # Get full model
    full_LM = fitlm(X,y,degree)
    X_quad = sm.add_constant(X)
    X_quad = getQuadratic(X_quad)

    # Find lowest performer in the model and remove it
    while True:
        red_LM, X_quad = getReducedModel(full_LM,X_quad,y)
        if lmNotSimilar(full_LM, red_LM):
            break
        full_LM = red_LM
    return full_LM

def lmNotSimilar(full, reduced):
    SSEr = reduced.mse_model*reduced.df_resid
    SSEf = full.mse_model*full.df_resid
    deltaSSE = SSEr-SSEf/(full.df_model - reduced.df_model)
    F_test = deltaSSE/full.mse_model
    pval = sp.stats.f.cdf(F_test,full.df_model - reduced.df_model, full.df_resid)
    if pval < 0.95:
        return False
    else:
        return True

def getReducedModel(LM,X,y):

    pVals = LM.pvalues
    k = 0
    maxVal = 0.0
    while True:
        if pVals.index[k] == "const":
            k += 1
            continue
        maxVal = pVals[k]
        maxID = k
        break
    for k in range(k+1,len(pVals)):
        if maxVal < pVals[k]:
            maxVal = pVals[k]
            maxID = k
    
    # Build new LM
    print("Attempting to Drop :", pVals.index[maxID])
    X = X.drop(pVals.index[maxID],axis=1)
    redLM = sm.OLS(y,X).fit()
    return redLM, X

def getInteractions(x):
    dfNames = x.columns
    for k1 in range(0, len(dfNames)):
        dfName1 = dfNames[k1]
        for k2 in range(k1+1, len(dfNames)):
            dfName2 = dfNames[k2]
            if dfName1 == dfName2 or dfName1 == "const" or dfName2 == "const":
                continue
            varName = dfName1 + ":" + dfName2
            varVale = x[dfName1]*x[dfName2]
            x[varName] = varVale
    return x

def getQuadratic(x):
    intTerms = []
    quadTerms = []
    dfNames = x.columns
    for k1 in range(0, len(dfNames)):
        dfName1 = dfNames[k1]
        if dfName1 == "const":
            continue
        for k2 in range(k1+1, len(dfNames)):
            dfName2 = dfNames[k2]
            if dfName1 == dfName2 or dfName2 == "const":
                continue
            varName = dfName1 + ":" + dfName2
            varVale = x[dfName1]*x[dfName2]
            intTerms.append([varName, varVale])
        varName = dfName1 + "^2"
        varVale = x[dfName1]*x[dfName1]
        quadTerms.append([varName, varVale])
    for intTerm in intTerms:
        x[intTerm[0]] = intTerm[1]
    for quadTerm in quadTerms:
        x[quadTerm[0]] = quadTerm[1]
    return x


def buildDF(x,degree):
    if degree == "linear":
        return x
    elif degree == "interactions":
        return getInteractions(x)
    elif degree == "quadratic":
        return getQuadratic(x)


def fitlm(x,y,degree = "linear"):
    
    x = sm.add_constant(x)
    x = buildDF(x,degree)
    
    lm = sm.OLS(y,x).fit()

    return lm

def getFun(t,z):
    fun = lambda b: b[0]*np.exp(-b[1]*t)*np.cos(2*np.pi*b[2]*t)
    minFun = lambda b: fun(b) - z
    return fun, minFun

def diff(x):
    dx = np.zeros(len(x)-1)
    dx = x[1:]-x[0:-1]
    return dx

def findFreq(t,z,doGrad = True):
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
        return findFreq(t,z_org,False)
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

    return freq, damp
    

if __name__ == "__main__":
    vars = readDOEparams()
    numFiles = len(glob("results/run*"))

    rid = np.array(vars[0][1])
    bld = np.array(vars[1][1])
    bym = np.array(vars[2][1])

    rho = 125.0
    diam = 0.00283
    lb = diam
    area = np.pi*(diam/2.0)*(diam/2.0)
    K = bym*area/lb
    vol = 4.0 * np.pi * (diam/2.0)*(diam/2.0)*(diam/2.0) / 3.0
    m = vol*rho
    w2 = K/m
    w = np.sqrt(w2)
    
    data = []
    for k in range(numFiles):
        fname = "results/run_%s/beam.csv" % str(k+1)
        curData = np.genfromtxt(fname = fname, delimiter = ",")
        data.append(curData[1:,:])
    
    frq = np.zeros(numFiles)
    bgd = np.zeros(numFiles)
    for k in range(numFiles):
        frq[k], bgd[k] = getProps(data[k])
    
    dfData = np.array([bld,bym,w,frq,bgd])
    dfColu = ["bld","bym","w","frq","bgd"]
    df = pd.DataFrame(
        data=dfData.transpose(),
        columns=dfColu
    )
    X1 = df["bgd"] # X = df[["frq","bgd"]]
    X2 = df[["frq","bgd"]] # X = df[["frq","bgd"]]
    y = df["bld"]

    md1 = fitlm(X1,y,"quadratic")
    # print(md1.summary())
    md2 = fitlm(X2,y,"interactions")
    # print(md2.summary())

    sLM = stepwiselm(X2,y)
    print(sLM.summary())

    # plt.plot(bld,bgd,'o')
    # plt.show()
