import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def line(x,a,b):
     return a*x+b

file_59=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmodeavoid")
prod1=file_59[:,0] 
rad1=file_59[:,1] 
mas1=file_59[:,2] 
f1=file_59[:,3]
p1=file_59[:,4]  
p2=file_59[:,5]

f11=f1*2.998*100/(6.28*6.28)
p11=p1*2.998*100/(6.28*6.28)
p22=p2*2.998*100/(6.28*6.28)


fig=plt.figure(constrained_layout=True)
plt.plot(prod1,f11,'-o',markersize=3,label=r'$f$')
plt.plot(prod1,p11,'-o',markersize=3,label=r'$p_1$')
plt.plot(prod1,p22,'-o',markersize=3,label=r'$p_2$')

     
plt.xlabel('$n/n_0$',fontsize = 16)
plt.ylabel('$\\nu$ (kHz)',fontsize = 16)
plt.xlim(0.7,1.6)
plt.legend()
fig.savefig("avoidedcrossing3.pdf",dpi=600)

plt.show()
