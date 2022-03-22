import numpy as np
import matplotlib.pyplot as plt

filelist=[]

for i in range(1,6):
    filelist.append("pmode_m%s.out" %i)

plt.title("P Mode frequencies")
plt.xlabel('Mass ($M_{\odot}$)')
plt.ylabel('F (kHz)')
plt.yticks()

for fname in filelist:
    data=np.loadtxt(fname)
    X=data[:,2]
    Y=data[:,3]
    plt.legend(["purple", "red", "green" , "blue" , "magenta"], loc = "upper left", bbox_to_anchor=(0.75, 1.15), ncol=1)		
    plt.plot(X,Y,'o-r')

plt.show()

