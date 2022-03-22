import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
                               
def line(x,a,b):
     return a*x+b
     
# for jsym

file_30=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_jsym30.out")
nb_30=file_30[:,0] 
r30=file_30[:,1]
m30=file_30[:,2]
f30=file_30[:,3]
c30=m30*1.4766/r30 
d30=np.sqrt(m30/1.4/(r30/10)**3)
wm30=m30*f30*2.*np.pi*1.4766	#30 jsym

file_32=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_jsym32.out")
nb_32=file_32[:,0] 
r32=file_32[:,1]
m32=file_32[:,2]
f32=file_32[:,3]
c32=m32*1.4766/r32 
d32=np.sqrt(m32/1.4/(r32/10)**3)
wm32=m32*f32*2.*np.pi*1.4766	#32 jsym

m=np.concatenate( [m30,m32], axis=0)	#mass
cf=np.concatenate( [c30,c32], axis=0 )  #compactness
df=np.concatenate( [d30,d32], axis=0 )	#density
f=np.concatenate( [f30,f32], axis=0 )	#frequency
wmf=np.concatenate( [wm30,wm32], axis=0 )	#wM (kHz Km)

#freq vs mass

fig=plt.figure(constrained_layout=True)

plt.plot(m30,f30,'-',label=r'$J_{sym}=30$')
plt.plot(m32,f32,'-',label=r'$J_{sym}=32$')


plt.xlabel('M ($M_{\odot}$)',fontsize = 18)
plt.ylabel('$\\nu$ (kHz)',fontsize = 18) 
plt.tick_params(axis = 'both', direction = 'in', which = 'both', right = True, top=True, labelsize=15)
plt.minorticks_on() 
plt.legend(fontsize = 16)
plt.show()
fig.savefig("massvsfreq_jsym.pdf",dpi=0)
