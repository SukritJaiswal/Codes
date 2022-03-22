import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def line(x,a,b):
	return a*x+b
     
# for 

filelist=[]

for i in range(1,44):
	filelist.append("yy%s" %i)

l = len(filelist)

#for j in range(l):
#	x = []
#	y = []
#	for line in open(filelist[j], 'r'):
#		lines = [k for k in line.split()]
#		x.append(lines[0])
#		y.append(lines[1])
		

j=0		
for fname in filelist:
	data=np.loadtxt(fname)
	X=data[:,0]
	Y=data[:,1]
	fig=plt.figure(constrained_layout=True)
	plt.plot(X,Y,'-*',label=r'$m^*=0.7m_N,prod=%s$'%(j+1))
	plt.axhline(y=0, color='r', linestyle='-')     
	plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
	plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
	plt.xlim(0.0,0.3)
	plt.ylim(-100.0,100.0)
	plt.legend()
	fig.savefig("m=0.7;prod=%s.pdf"%(j+1),dpi=600)
	j=j+1
#	plt.show() 

	
	
