import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

def line(x,a,b):
     return a*x+b
     
#for ksat

file_240=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_ksat240.out")
nb_240=file_240[:,0] 
r240=file_240[:,1]
m240=file_240[:,2]
f240=file_240[:,3]
c240=m240*1.4766/r240 
d240=np.sqrt(m240/1.4/(r240/10)**3)
wm240=m240*f240*2.*np.pi*1.4766	#0.240 effm

file_250=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_ksat250.out")
nb_250=file_250[:,0] 
r250=file_250[:,1]
m250=file_250[:,2]
f250=file_250[:,3]
c250=m250*1.4766/r250 
d250=np.sqrt(m250/1.4/(r250/10)**3)
wm250=m250*f250*2.*np.pi*1.4766	#0.250 effm

file_260=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_ksat260.out")
nb_260=file_260[:,0] 
r260=file_260[:,1]
m260=file_260[:,2]
f260=file_260[:,3]
c260=m260*1.4766/r260 
d260=np.sqrt(m260/1.4/(r260/10)**3)
wm260=m260*f260*2.*np.pi*1.4766	#0.260 effm

file_270=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_ksat270.out")
nb_270=file_270[:,0] 
r270=file_270[:,1]
m270=file_270[:,2]
f270=file_270[:,3]
c270=m270*1.4766/r270 
d270=np.sqrt(m270/1.4/(r270/10)**3)
wm270=m270*f270*2.*np.pi*1.4766	#0.270 effm

file_280=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_ksat280.out")
nb_280=file_280[:,0] 
r280=file_280[:,1]
m280=file_280[:,2]
f280=file_280[:,3]
c280=m280*1.4766/r280 
d280=np.sqrt(m280/1.4/(r280/10)**3)
wm280=m280*f280*2.*np.pi*1.4766	#0.280 effm



m=np.concatenate( [m240,m250,m260,m270,m280], axis=0)	#mass
cf=np.concatenate( [c240,c250,c260,c270,c280], axis=0 )  #compactness
df=np.concatenate( [d240,d250,d260,d270,d280], axis=0 )	#density
f=np.concatenate( [f240,f250,f260,f270,f280], axis=0 )	#frequency
wmf=np.concatenate( [wm240,wm250,wm260,wm270,wm280], axis=0 )	#wM (kHz Km)

#freq vs mass

fig=plt.figure(constrained_layout=True)

plt.plot(m240,f240,'-',label=r'$K_{sat}=240$')
plt.plot(m250,f250,'-',label=r'$K_{sat}=250$')
plt.plot(m260,f260,'-',label=r'$K_{sat}=260$')
plt.plot(m270,f270,'-',label=r'$K_{sat}=270$')
plt.plot(m280,f280,'-',label=r'$K_{sat}=280$')


plt.xlabel('M ($M_{\odot}$)',fontsize = 18)
plt.ylabel('$\\nu$ (kHz)',fontsize = 18) 
plt.tick_params(axis = 'both', direction = 'in', which = 'both', right = True, top=True, labelsize=15)
plt.minorticks_on() 
plt.legend(fontsize = 16)
plt.show()
fig.savefig("massvsfreq_ksat.pdf",dpi=600)
