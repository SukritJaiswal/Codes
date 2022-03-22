import matplotlib.pyplot as plt
import numpy as np


file_59=np.loadtxt("/home/sukrit/Desktop/Combine Code/lowm_m0.55.out")
nb59=file_59[:,0] #prod nb/n0
rho59=file_59[:,1] #nb (fm^3)
r59=file_59[:,2] #radius (km)
m59=file_59[:,3] #mass (msol)
f59=file_59[:,4] #f mode freq 
p59=file_59[:,5] #p1 mode
p59_2=file_59[:,6] #p2 mode
dimtidal=file_59[:,7] #dimless tidal
tidal=file_59[:,8] #tidal

file_60=np.loadtxt("/home/sukrit/Desktop/Combine Code/lowm_m0.60.out")
nb60=file_60[:,0] #prod nb/n0
rho60=file_60[:,1] #nb (fm^3)
r60=file_60[:,2] #radius (km)
m60=file_60[:,3] #mass (msol)
f60=file_60[:,4] #f mode freq 
p60=file_60[:,5] #p1 mode
p60_2=file_60[:,6] #p2 mode
dimtidal=file_60[:,7] #dimless tidal
tidal=file_60[:,8] #tidal

file_61=np.loadtxt("/home/sukrit/Desktop/Combine Code/lowm_m0.65.out")
nb61=file_61[:,0] #prod nb/n0
rho61=file_61[:,1] #nb (fm^3)
r61=file_61[:,2] #radius (km)
m61=file_61[:,3] #mass (msol)
f61=file_61[:,4] #f mode freq 
p61=file_61[:,5] #p1 mode
p61_2=file_61[:,6] #p2 mode
dimtidal=file_61[:,7] #dimless tidal
tidal=file_61[:,8] #tidal

file_62=np.loadtxt("/home/sukrit/Desktop/Combine Code/lowm_m0.70.out")
nb62=file_62[:,0] #prod nb/n0
rho62=file_62[:,1] #nb (fm^3)
r62=file_62[:,2] #radius (km)
m62=file_62[:,3] #mass (msol)
f62=file_62[:,4] #f mode freq 
p62=file_62[:,5] #p1 mode
p62_2=file_62[:,6] #p2 mode
dimtidal=file_62[:,7] #dimless tidal
tidal=file_62[:,8] #tidal

file_63=np.loadtxt("/home/sukrit/Desktop/Combine Code/lowm_m0.75.out")
nb63=file_63[:,0] #prod nb/n0
rho63=file_63[:,1] #nb (fm^3)
r63=file_63[:,2] #radius (km)
m63=file_63[:,3] #mass (msol)
f63=file_63[:,4] #f mode freq 
p63=file_63[:,5] #p1 mode
p63_2=file_63[:,6] #p2 mode
dimtidal=file_63[:,7] #dimless tidal
tidal=file_63[:,8] #tidal

fig=plt.figure(constrained_layout=True)
plt.plot(m59,f59,'-',label=r'$m^*/m=0.55$')
plt.plot(m60,f60,'-',label=r'$m^*/m=0.60$')
plt.plot(m61,f61,'-',label=r'$m^*/m=0.65$')
plt.plot(m62,f62,'-',label=r'$m^*/m=0.70$')
plt.plot(m63,f63,'-',label=r'$m^*/m=0.75$')

plt.xlabel('$M (M_{\odot})$',fontsize = 18)
plt.ylabel('$\\nu$ (kHz)',fontsize = 18)
plt.xlim(0,0.5)
plt.ylim(0,3)
plt.tick_params(axis = 'both', direction = 'in', which = 'both', right = True, top=True, labelsize=15)
plt.minorticks_on() 
plt.legend(fontsize = 16)
fig.savefig("massvsfreqlowmasseffm.pdf",dpi=600)

plt.show()

'''
fig=plt.figure(constrained_layout=True)
plt.plot(nb1,f1,'-',label=r'$f-$mode')
plt.plot(nb1,p1,'-',label=r'$p_1-$mode')
plt.plot(nb1,p2,'-',label=r'$p_2-$mode')
     
#plt.colorbar(m1,label="M $(M_{\odot})$", orientation="Vertical")     
plt.xlabel('$\\rho/\\rho_{sat}$',fontsize = 18)
plt.ylabel('$\\nu$ (kHz)',fontsize = 18)
plt.tick_params(axis = 'both', direction = 'in', which = 'both', right = True, top=True, labelsize=15)
plt.minorticks_on() 
plt.legend(fontsize = 16)
fig.savefig("rho=0.15,esat=16,ksat=240,j32,lsym60,m=0.7;lowmass.pdf",dpi=600)

plt.show()
'''
