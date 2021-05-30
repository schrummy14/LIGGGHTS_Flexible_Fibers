'''
the black_box package is also needed
https://github.com/schrummy14/blackbox
This is a fork of the original package so if you follow the directions in the readme,
you will need to modify the search_min function yourself following the changes in
commit number 'b4527b6b5bb722944a5773f4853ff8d331749fc3'
'''
import os
import sys
import glob
import typing

import black_box as bb
import numpy as np
import bayOpt.bayOpt as bo
import dataRead.extractProperties as ep

from liggghts import liggghts

global _count
_count = 0

_lmp = liggghts("",["-screen","/dev/null","-log","/dev/null"])

def main():
    x0, y0 = getPreviousRuns()
    bounds = getBounds()
    xIter, yIter = doOpt(x0, y0, bounds, 2)
    if not yIter is None:
        np.savetxt("xIter.csv", xIter)
        np.savetxt("yIter.csv", yIter)
        yBest = -100
        for k in range(len(yIter)):
            if yIter[k] > yBest:
                yBest = yIter[k]
                xBest = xIter[k]
    else:
        xBest = xIter

def doOpt(x0, y0, bounds, algo=1):
    if algo == 1:
        x, y = bo.bayesian_optimisation(25, sample_loss, bounds, x0, y0, 5)
    elif algo == 2:
        x = bb.search_min(sample_loss2, bounds, 50, 1, "outfile.csv",executor=None,randseed=123457)
        y = None  

    return x, y
    
def runBest(x):
    global _count
    _count += 1
    print("Running Best Sim")
    params = {
        "bond_youngs_modulus": x[0],
        "bond_damp_val": x[1]
    }
    params.update({"run_num": _count})
    run_liggghts(params) # This creates a file that holds the simulation data
    frq, bgd = get_data() # Now we need to extract the data and make a conclusion
    e1 = (frq-60.0)/60.0
    e2 = (bgd-0.15)/0.15
    print("Error in x1 = %f, x2 = %f" % (100.0*e1, 100.0*e2))

def getBounds():
    return np.array([
    [9.0e8, 2.0e9],
    [0.0009, 0.01]
    ])

def getPreviousRuns():
    x0 = None
    y0 = None
    try:
        x0 = np.loadtxt("xIter.csv")
        y0 = np.loadtxt("yIter.csv")
        return x0, y0
    except:
        return None, None

def run_liggghts(params: typing.Dict[str, float]) -> None:
    _lmp.command("clear")
    set_liggghts_variables(params)
    _lmp.file("in.liggghts")

def sample_loss2(x: typing.List[float]) -> float:
    global _count
    _count += 1
    params = {
        "bond_youngs_modulus": x[0],
        "bond_damp_val": x[1]
    }
    params.update({"run_num": _count})
    # Need to check if we have already ran this iteration
    # If we have, simply return the result, if not then we need to run liggghts
    data = None
    if os.path.isfile('black_box_runs.csv'):
        with open('black_box_runs.csv','r') as f:
            data = f.readlines()
        for curData in data:
            d = curData.split(',')
            if int(d[0]) == _count:
                x0err = (x[0] - float(d[1]))/x[0]
                x1err = (x[1] - float(d[2]))/x[1]
                if abs(x0err) < 0.001 and abs(x1err) < 0.001:
                    return float(d[3])

    
    run_liggghts(params) # This creates a file that holds the simulation data
    frq, bgd = get_data() # Now we need to extract the data and make a conclusion
    e1 = (frq-60.0)/60.0
    e2 = (bgd-0.15)/0.15
    funErr = np.linalg.norm(np.array([e1, e2]))
    print("Error in x1 = %f, x2 = %f has function error = %f" % (100.0*e1, 100.0*e2, funErr))
    with open("black_box_runs.csv", 'a') as f:
        f.write('%d, %f, %f, %f\n' % (_count, x[0], x[1], funErr))
    return funErr

def sample_loss(x: typing.List[float]) -> float:
    global _count
    _count += 1
    print("Running Bay Opt Sim: %d" % (_count))
    params = {
        "bond_youngs_modulus": x[0],
        "bond_damp_val": x[1]
    }
    params.update({"run_num": _count})
    run_liggghts(params) # This creates a file that holds the simulation data
    frq, bgd = get_data() # Now we need to extract the data and make a conclusion
    e1 = (frq-60.0)/60.0
    e2 = (bgd-0.15)/0.15
    funErr = np.linalg.norm(np.array([e1, e2]))
    print("Error in x1 = %f, x2 = %f has function error = %f" % (100.0*e1, 100.0*e2, funErr))
    return -funErr

def get_data() -> typing.Tuple[float, float]:
    global _count
    fname = "results/run_%s/beam.csv" % str(_count)
    curData = np.genfromtxt(fname = fname, delimiter = ",")
    data = curData[1:,:]
    frq, bgd = ep.getProps(data)
    return frq, bgd

def set_liggghts_variables(params: typing.Dict[str, float]) -> None:
    # Set simulation variables
    for key in params.keys():
        curStr = "variable %s equal %e" % (key, params[key])
        _lmp.command(curStr)

if __name__ == "__main__":
    main()
