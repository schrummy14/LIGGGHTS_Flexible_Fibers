import glob
import numpy as np
import matplotlib.pyplot as plt

def main():
    files = glob.glob("compression_*.csv")
    trueValues = np.genfromtxt("compression_1.csv",delimiter=',')

    for f in files:
        k = int(f[-5])
        if k == 1:
            continue
        vals = np.genfromtxt(f,delimiter=',')
        res = 100.0*np.linalg.norm(trueValues[1:,:]-vals[1:,:])/np.linalg.norm(trueValues[1:,:])
        plt.plot(k,res,'o')
        # plt.plot(trueValues[1:,0],trueValues[1:,2],trueValues[1:,0],vals[1:,2])
        # plt.show()
        print(k,res)
    plt.show()


if __name__ == "__main__":
    main()