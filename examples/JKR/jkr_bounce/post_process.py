import numpy as np
import subprocess


W = 0.2
nu = 0.3
E = 215e9
rho = 8300.
Rstar = 5e-6

Estar = 1/(2*(1-nu**2)/E)


K = 1.052440624204792

Vs = 1.84*((W/Rstar)**5./(rho**3.*Estar**2.))**(1./6.)


Vi = 0.03
en = (1-(Vs/Vi)**2.)**0.5

Vout_analytical = en*Vi


filename = 'post/liggghts.dump'
line = subprocess.check_output(['tail', '-1', filename]).decode('UTF-8').split()
Vout_sim = -1.*float(line[9])

print("Simulated bounce velocity is {:.4f} m/s while analytical solution is {:.4f} m/s.\nError is {:.2f} %".
      format(Vout_sim,Vout_analytical,100*abs(Vout_analytical-Vout_sim)/Vout_analytical))


