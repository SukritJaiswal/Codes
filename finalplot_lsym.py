import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


def line(x,a,b):
     return a*x+b
     
# for lsym

file_50=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_lsym50.out")
nb_50=file_50[:,0] 
r50=file_50[:,1]
m50=file_50[:,2]
f50=file_50[:,3]
c50=m50*1.4766/r50 
d50=np.sqrt(m50/1.4/(r50/10)**3)
wm50=m50*f50*2.*np.pi*1.4766	#50 lsym

file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_lsym55.out")
nb_55=file_55[:,0] 
r55=file_55[:,1]
m55=file_55[:,2]
f55=file_55[:,3]
c55=m55*1.4766/r55 
d55=np.sqrt(m55/1.4/(r55/10)**3)
wm55=m55*f55*2.*np.pi*1.4766	#55 lsym

file_60=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_lsym60.out")
nb_60=file_60[:,0] 
r60=file_60[:,1]
m60=file_60[:,2]
f60=file_60[:,3]
c60=m60*1.4766/r60 
d60=np.sqrt(m60/1.4/(r60/10)**3)
wm60=m60*f60*2.*np.pi*1.4766	#60 lsym

m=np.concatenate( [m50,m55,m60], axis=0)	#mass
cf=np.concatenate( [c50,c55,c60], axis=0 )  #compactness
df=np.concatenate( [d50,d55,d60], axis=0 )	#density
f=np.concatenate( [f50,f55,f60], axis=0 )	#frequency
wmf=np.concatenate( [wm50,wm55,wm60], axis=0 )	#wM (kHz Km)

#freq vs mass

fig=plt.figure(constrained_layout=True)

plt.plot(m50,f50,'-',label=r'$L_{sym}=50$')
plt.plot(m55,f55,'-',label=r'$L_{sym}=55$')
plt.plot(m60,f60,'-',label=r'$L_{sym}=60$')


plt.xlabel('M ($M_{\odot}$)',fontsize = 18)
plt.ylabel('$\\nu$ (kHz)',fontsize = 18) 
plt.tick_params(axis = 'both', direction = 'in', which = 'both', right = True, top=True, labelsize=15)
plt.minorticks_on() 
plt.legend(fontsize = 16)
plt.show()
fig.savefig("massvsfreq_lsym.pdf",dpi=600)

#densityvsfreq

fig=plt.figure(constrained_layout=True)

plt.plot(d50,f50,'-',label=r'$L_{sym}=50$')
plt.plot(d55,f55,'-',label=r'$L_{sym}=55$')
plt.plot(d60,f60,'-',label=r'$L_{sym}=60$')


plt.xlabel('$(\\frac{M}{\\bar{R}^3})^{1/2}$',fontsize = 18)
plt.ylabel('$\\nu$ (kHz)',fontsize = 18) 
plt.tick_params(axis = 'both', direction = 'in', which = 'both', right = True, top=True, labelsize=15)
plt.minorticks_on() 
plt.legend(fontsize = 16)
plt.show()
fig.savefig("densityvsfreq_lsym.pdf",dpi=600)
