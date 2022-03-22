import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def line(x,a,b):
     return a*x+b
     
# for GM1 nolamomega

file_gm = np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_gm1_nolamomega_fulleos.out")
nb_gm=file_gm[:,0] 
rgm=file_gm[:,1]
mgm=file_gm[:,2]
fgm=file_gm[:,3]
cgm=mgm*1.4766/rgm 
dgm=np.sqrt(mgm/1.4/(rgm/10)**3)
wmgm=mgm*fgm*2.*np.pi*1.4766

#freq vs mass

fig=plt.figure(constrained_layout=True)

plt.plot(mgm,fgm,'-',marker=".",markersize=5,label=r'GM1')

plt.xlabel('M ($M_{\odot}$)',fontsize = 16)
plt.ylabel('$\\nu$ (kHz)',fontsize = 16) 
plt.legend()
plt.show()
fig.savefig("massvsfreq_gm1_nolamomega_fulleos.pdf",dpi=600)
