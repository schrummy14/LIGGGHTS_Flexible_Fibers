
import os
import numpy as np

if __name__ == "__main__":
  # Read in data
  d = np.genfromtxt('post' + os.sep + '3-point.csv', delimiter=',', skip_header=1)

  t = d[:,0]
  x = d[:,1]
  f = d[:,2]

  x *= -1 # Flip the sign of the displacement 

  # Get just the compression 
  dx = np.gradient(x,t)
  x = x[dx>0.002]
  f = f[dx>0.002]

  x = x[f>0]
  f = f[f>0]

  # Get linear fit mx+b = f
  n = len(x)
  A = np.zeros([n,2])
  A[:,0] = x
  A[:,1] = 1

  m, b = np.linalg.lstsq(A, f, rcond=None)[0]
  
  with open('..' + os.sep + 'results.txt','w') as f:
    f.write('%e\n%e\n'%(m,b))