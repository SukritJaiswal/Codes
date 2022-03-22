import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def line(x,a,b):
     return a*x+b

file_59=np.loadtxt("/home/sukrit/Desktop/Combine Code/TovAll")
m1=file_59[:,0] 
r1=file_59[:,1] 
m2=file_59[:,2] 
r2=file_59[:,3]
m3=file_59[:,4]  
r3=file_59[:,5]
m4=file_59[:,6]
r4=file_59[:,7]
m5=file_59[:,8]
r5=file_59[:,9]
m6=file_59[:,10]
r6=file_59[:,11]
m7=file_59[:,12]
r7=file_59[:,13]
m8=file_59[:,14]
r8=file_59[:,15]
m9=file_59[:,16]
r9=file_59[:,17]
m10=file_59[:,18]
r10=file_59[:,19]
m11=file_59[:,20]
r11=file_59[:,21]
m12=file_59[:,22]
r12=file_59[:,23]
m13=file_59[:,24]
r13=file_59[:,25]
m14=file_59[:,26]
r14=file_59[:,27]
m15=file_59[:,28]
r15=file_59[:,29]
m16=file_59[:,30]
r16=file_59[:,31]

fig=plt.figure(constrained_layout=True)
plt.plot(r1,m1,'--')
plt.plot(r2,m2,'--')
plt.plot(r3,m3,'--')
plt.plot(r4,m4,'--')
plt.plot(r5,m5,'--')
plt.plot(r6,m6,'--')
plt.plot(r7,m7,'--')
plt.plot(r8,m8,'--')
plt.plot(r9,m9,'--')
plt.plot(r10,m10,'--')
plt.plot(r11,m11,'--')
plt.plot(r12,m12,'--')
plt.plot(r13,m13,'--')
plt.plot(r14,m14,'--')
plt.plot(r15,m15,'--')
plt.plot(r16,m16,'--')
     
plt.xlabel('Radius ($km$)',fontsize = 16)
plt.ylabel('Mass ($M_{\odot}$)',fontsize = 16)
#plt.xlim(10,16)
#plt.ylim(0,3)
plt.legend()
fig.savefig("AllTov.pdf",dpi=600)

plt.show()


