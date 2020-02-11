import numpy as np
import subprocess



Vi = 0.03
en = 0.5

Vout_analytical = en*Vi


filename = 'post/liggghts.dump'
line = subprocess.check_output(['tail', '-1', filename]).decode('UTF-8').split()
Vout_sim = -1.*float(line[9])

print("Simulated bounce velocity is {:.4f} m/s while analytical solution is {:.4f} m/s.\nError is {:.2f} %".
      format(Vout_sim,Vout_analytical,100*abs(Vout_analytical-Vout_sim)/Vout_analytical))


