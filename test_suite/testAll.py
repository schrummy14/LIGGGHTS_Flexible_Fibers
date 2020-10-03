import os
import sys
import subprocess
import numpy as np

LIGGGHTS = '../../../src/lmp_auto'

def singeCheck():
    process = subprocess.Popen(
                ['sh', 'run.sh', LIGGGHTS, sys.executable],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    return

if __name__ == "__main__":
    folders = next(os.walk('.'))[1]
    for fold in folders:
        print("Testing", fold)
        os.chdir(fold)
        try:
            d2 = np.genfromtxt('truth.txt')
        except:
            print("No truth.txt file found, skipping...")
            os.chdir('..')
            continue
        # os.chdir(fold)
        singeCheck()
        d1 = np.genfromtxt('results.txt')
        # d2 = np.genfromtxt('truth.txt')
        #d = 0
        #for k in range(len(d2)):
        #    d += np.abs(d1-d2)/np.max(np.array([1,d2[k]]))
        d = np.abs(d1-d2)/(np.abs(d2)+1.0e-10)
        d = np.sum(d)
        if d < 1.0e-10:
            print(fold, ': pass')
            os.remove('results.txt')
        else:
            print(fold, ' : fail')

        os.chdir('..')
