import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

def line(x,a,b):
     return a*x+b
     
# for esat

file_155=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_esat15.5.out")
nb_155=file_155[:,0] 
r155=file_155[:,1]
m155=file_155[:,2]
f155=file_155[:,3]
c155=m155*1.4766/r155 
d155=np.sqrt(m155/1.4/(r155/10)**3)
wm155=m155*f155*2.*np.pi*1.4766	#-15.5 esat

file_160=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_esat16.0.out")
nb_160=file_160[:,0] 
r160=file_160[:,1]
m160=file_160[:,2]
f160=file_160[:,3]
c160=m160*1.4766/r160 
d160=np.sqrt(m160/1.4/(r160/10)**3)
wm160=m160*f160*2.*np.pi*1.4766	#-16.0 esat

file_165=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_esat16.5.out")
nb_165=file_165[:,0] 
r165=file_165[:,1]
m165=file_165[:,2]
f165=file_165[:,3]
c165=m165*1.4766/r165 
d165=np.sqrt(m165/1.4/(r165/10)**3)
wm165=m165*f165*2.*np.pi*1.4766	#-16.5 esat

m=np.concatenate( [m155,m160,m165], axis=0)	#mass
cf=np.concatenate( [c155,c160,c165], axis=0 )  #compactness
df=np.concatenate( [d155,d160,d165], axis=0 )	#density
f=np.concatenate( [f155,f160,f165], axis=0 )	#frequency
wmf=np.concatenate( [wm155,wm160,wm165], axis=0 )	#wM (kHz Km)

#freq vs mass

fig=plt.figure(constrained_layout=True)

plt.plot(m155,f155,'-',label=r'$E_{sat}=15.5$')
plt.plot(m160,f160,'-',label=r'$E_{sat}=16.0$')
plt.plot(m165,f165,'-',label=r'$E_{sat}=16.5$')


plt.xlabel('M ($M_{\odot}$)',fontsize = 18)
plt.ylabel('$\\nu$ (kHz)',fontsize = 18)
plt.tick_params(axis = 'both', direction = 'in', which = 'both', right = True, top=True, labelsize=15)
plt.minorticks_on() 
plt.legend(fontsize = 16)
plt.show()
fig.savefig("massvsfreq_esat.pdf",dpi=600)
