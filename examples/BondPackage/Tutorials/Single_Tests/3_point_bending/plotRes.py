import os
import glob
import numpy as np
import matplotlib.pyplot as plt

def main():
    allData = []
    fname = 'post' + os.path.sep + '3-point.csv'
    fnames = glob.glob('post' + os.path.sep + '3-point_*.csv')
    is1 = 0
    k = -1
    for fname in fnames:
        k += 1
        numThreads = int(fname[-5])
        if numThreads == 1:
            is1 = k
        data = np.genfromtxt(fname,delimiter=',')
        data = data[1:,:]
        t = data[:,0]
        z = -data[:,1]
        fz = data[:,2]
        allData.append(np.array([z,fz]))

        # plt.plot(z,fz)
    z = allData[is1][:,0]
    fz = allData[is1][:,1]
    for k in range(len(allData)):
        if k == is1:
            continue
        ffz= allData[k][:,1]
        plt.plot(z,np.log10(abs(fz-ffz)))

    plt.show()

if __name__ == "__main__":
    main()  