import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def line(x,a,b):
     return a*x+b
     
# for gm1

file_gm = np.loadtxt("/home/sukrit/Desktop/Combine Code/fmode_lowmass_0.65.out")
nb_gm=file_gm[:,0] 
rgm=file_gm[:,1]
egm=file_gm[:,2]
pgm=file_gm[:,3]
#cgm=mgm*1.4766/rgm 
#dgm=np.sqrt(mgm/1.4/(rgm/10)**3)
#wmgm=mgm*fgm*2.*np.pi*1.4766   
     
#file_gm = np.loadtxt("/home/sukrit/Desktop/Combine Code/eosanmisovec_m0.65lsym60j32_1.out")
#nb_gm=file_gm[:,0] 
#rgm=file_gm[:,1]
#egm=file_gm[:,2]
#pgm=file_gm[:,3]
          
#freq vs mass

fig=plt.figure(constrained_layout=True)

plt.plot(egm,pgm,'-',marker=".",markersize=5,label=r'0.65 m*')

plt.xlabel('M $(M_{\odot})$',fontsize = 16)
plt.ylabel('$\\nu$ (kHz)' ,fontsize = 16) 
plt.legend()
plt.show()
fig.savefig("fmode_lowmass_m0.55.pdf",dpi=600)
