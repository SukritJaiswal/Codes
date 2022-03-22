import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def line(x,a,b):
     return a*x+b

file_59=np.loadtxt("/home/sukrit/Desktop/Combine Code/calceosanmisovec0.55sigma.out")
nb1=file_59[:,0] 
rho1=file_59[:,1] 
en1=file_59[:,2] 
pr1=file_59[:,3]
bn1=file_59[:,4]  
sigma1=file_59[:,5]

fig=plt.figure(constrained_layout=True)
plt.plot(nb1,sigma1,'-*',label=r'$m^*=0.55m_N$')
     
plt.xlabel('$n_b/n_0$',fontsize = 16)
plt.ylabel('$\\sigma$ value (MeV)',fontsize = 16)
plt.legend()
fig.savefig("m=0.55;sigma graph.pdf",dpi=600)

plt.show()


