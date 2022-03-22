import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def line(x,a,b):
     return a*x+b
     
# for 

file_0 = np.loadtxt("/home/sukrit/Desktop/Combine Code/boom1.out")
nb0=file_0[:,0] 
r0=file_0[:,1]
m0=file_0[:,2]
f0=file_0[:,3]
c0=m0*1.4766/r0 
d0=np.sqrt(m0/1.4/(r0/10)**3)
wm0=m0*f0*2.*np.pi*1.4766   

file_1 = np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_sotani_4.out")
nb1=file_1[:,0] 
r1=file_1[:,1]
m1=file_1[:,2]
f1=file_1[:,3]
c1=m1*1.4766/r1 
d1=np.sqrt(m1/1.4/(r1/10)**3)
wm1=m1*f1*2.*np.pi*1.4766   

file_2 = np.loadtxt("/home/sukrit/Desktop/Combine Code/boom3.out")
nb2=file_2[:,0] 
r2=file_2[:,1]
m2=file_2[:,2]
f2=file_2[:,3]
c2=m2*1.4766/r2 
d2=np.sqrt(m2/1.4/(r2/10)**3)
wm2=m2*f2*2.*np.pi*1.4766   

#freq vs mass

fig=plt.figure(constrained_layout=True)

plt.plot(nb0,f0,'-',marker=".",markersize=5,label=r'$f$ mode')
plt.plot(nb1,f1,'-',marker=".",markersize=5,label=r'$p_1$ mode')
plt.plot(nb2,f2,'-',marker=".",markersize=5,label=r'$p_2$ mode')


plt.xlabel('M $(M_{\odot})$',fontsize = 16)
plt.ylabel('$\\nu$ (kHz)' ,fontsize = 16) 
#plt.xlim(0,0.5)
plt.legend()
plt.show()
fig.savefig("pmode_low_m0.7**.pdf",dpi=600)
