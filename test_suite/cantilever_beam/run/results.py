
import os
import numpy as np

if __name__ == "__main__":
    d = np.genfromtxt('beam.csv', delimiter=',', skip_header=1)
    n = np.size(d)
    d = np.reshape(d,n)
    with open('..' + os.sep + 'results.txt','w') as f:
        for k in range(n):
            f.write('%e\n'%(d[k]))
