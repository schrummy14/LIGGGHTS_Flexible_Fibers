
import glob
import numpy as np
from os.path import sep as fsep
import matplotlib.pyplot as plt

def readLIGGGHTS(fName):
    data = []
    with open(fName,'r') as f:
        line = f.readline()
        line = f.readline()
        timeStep = int(line)
        line = f.readline()
        line = f.readline()
        numAtoms = int(line.strip())
        line = f.readline()
        line = f.readline()
        line = f.readline()
        line = f.readline()
        line = f.readline()
        fdata = []
        for k in range(numAtoms):
            line = f.readline()
            line.strip()
            fdata.append(eval('[' + line.replace(' ',',') + ']'))
    return timeStep, np.array(fdata)

def main():
    fileNames = glob.glob('post1'+fsep+'compression'+fsep+'dump*.liggghts')
    numFN = len(fileNames)
    numPosts = len(glob.glob('post*'))
    data = np.zeros([len(fileNames),numPosts-1])
    for k in range(1,numPosts):
        kk = -1
        tStep = []
        for f in fileNames:
            kk += 1
            print("Reading in:", f)
            tt, tData = readLIGGGHTS(f)
            tData = tData[tData[:,0].argsort()]
            cf = f
            cf = cf.replace('post1', 'post' + str(int(np.power(2.0,k))))
            print("Reading in:", cf)
            ct, cData = readLIGGGHTS(cf)
            cData = cData[cData[:,0].argsort()]
            res = cData-tData
            tStep.append(ct)
            data[kk,k-1] = np.linalg.norm(res)
        tStep = np.array(tStep)
        ids = tStep.argsort()
        data[:,k-1] = data[ids, k-1]
    tStep = tStep[ids]
    for k in range(numPosts-1):
        plt.plot(tStep,np.log10(data[:,k]+2.0e-16), label = 'Procs ' + str(k+2))
    plt.legend()
    print(np.linalg.norm(data))
    plt.show()

if __name__ == "__main__":
    main()
