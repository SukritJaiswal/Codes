#python code for sigma-omega model

import numpy as np
from math import pi
from scipy.integrate import quad
import matplotlib.pyplot as plt

global hbarc,mn,ms,mw,gs,gw,gam

def scalarden(kF,mstar):
    EFstar = np.sqrt(kF**2+mstar**2)
    return gam*mstar*(kF*EFstar-mstar**2*np.log((kF+EFstar)/mstar))/(4*pi**2)

e_conv_fact = 1.7827e12 # [MeV/fm^3] to [g/cm^3]
p_conv_fact = 1.6022e33 # [MeV/fm^3] to [dyn/cm^2]

hbarc = 197.3269631
mn = 939
ms,mw = 550,783
gs,gw = 9.5727,11.6711

#specify g factor: 2 for neutron matter, 4 for nuclear matter
################################################################
gam = 2
################################################################

imax = 10000
accy = 1e-4

kmin,kmax,kstep=0,1.75,0.02                               #in units of fm^-1
kFvals = hbarc*np.arange(kmin,kmax+kstep,kstep)


rhoBvals,mnstarvals,epsvals,prsvals=[],[],[],[]
for kF in kFvals:                                   #unit conversion to MeV
    rhoB = gam*kF**3/(6*pi**2)

    mnstar,mnstar_old = mn,mn                         #starting point
    for i in range(0,imax):
        rhoS = scalarden(kF,mnstar)
        mnstar = mn - (gs/ms)**2*rhoS
        if abs(mnstar-mnstar_old) < accy:
            break
        mnstar_old = mnstar

    integrand_eps = lambda k,mnstar: k**2*np.sqrt(k**2+mnstar**2) 
    integrand_prs = lambda k,mnstar: k**4/np.sqrt(k**2+mnstar**2)

    integral_eps,error = quad(integrand_eps,0,kF,args=(mnstar))
    integral_prs,error = quad(integrand_prs,0,kF,args=(mnstar))

    eps = (gw*rhoB/mw)**2/2 + (ms*(mn-mnstar)/gs)**2/2 + gam*integral_eps/(2*pi**2)
    prs = (gw*rhoB/mw)**2/2 - (ms*(mn-mnstar)/gs)**2/2 + gam*integral_prs/(6*pi**2)

    print(f'{rhoB/hbarc**3:f}\t{eps/hbarc**3*e_conv_fact:f}\t{prs/hbarc**3*p_conv_fact:f}')
    
    rhoBvals.append(rhoB)
    mnstarvals.append(mnstar)
    epsvals.append(eps)
    prsvals.append(prs)

kFvals=kFvals/hbarc
rhoBvals=np.array(rhoBvals,float)/hbarc**3
mnstarvals=np.array(mnstarvals,float)/mn
epsvals=np.array(epsvals,float)/hbarc**3
prsvals=np.array(prsvals,float)/hbarc**3

fig=plt.figure(constrained_layout=True)
plt.plot(kFvals,epsvals/rhoBvals-mn)
plt.xlabel('$k_F$',fontsize = 16)
plt.ylabel('Binding Energy (MeV)',fontsize = 16) 

plt.show()

fig.savefig("EoS $\sigma - \omega$",dpi=600)

