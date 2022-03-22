import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def line(x,a,b):
     return a*x+b
     
# for 


file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz1")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.119$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.119.pdf",dpi=600)

plt.show()  


file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz2")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.121$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.121.pdf",dpi=600)

plt.show()  



file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz3")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.124$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.124.pdf",dpi=600)

plt.show()  



file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz4")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.127$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.127.pdf",dpi=600)

plt.show()  



file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz5")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.130$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.130.pdf",dpi=600)

plt.show()  



file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz6")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.133$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.133.pdf",dpi=600)

plt.show() 



file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz7")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.136$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.136.pdf",dpi=600)

plt.show() 


file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz8")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.139$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.139.pdf",dpi=600)

plt.show() 


file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz9")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.142$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.142.pdf",dpi=600)

plt.show() 


file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz10")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.145$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.145.pdf",dpi=600)

plt.show() 


file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz11")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.148$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.148.pdf",dpi=600)

plt.show() 


file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz12")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.151$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.151.pdf",dpi=600)

plt.show() 


file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz13")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.154$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.154.pdf",dpi=600)

plt.show() 


file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz14")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.157$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.157.pdf",dpi=600)

plt.show() 


file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz15")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.160$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.160.pdf",dpi=600)

plt.show() 


file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz16")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.163$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.163.pdf",dpi=600)

plt.show() 


file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz17")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.167$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.167.pdf",dpi=600)

plt.show() 


file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz18")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.170$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.170.pdf",dpi=600)

plt.show() 


file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz19")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.173$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.173.pdf",dpi=600)

plt.show() 


file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz20")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.176$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.176.pdf",dpi=600)

plt.show() 


file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz21")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.180$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.180.pdf",dpi=600)

plt.show() 


file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz22")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.183$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.183.pdf",dpi=600)

plt.show() 


file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz23")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.186$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.186.pdf",dpi=600)

plt.show() 


file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz24")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.190$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.190.pdf",dpi=600)

plt.show() 


file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz25")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.193$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.193.pdf",dpi=600)

plt.show() 


file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz26")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.196$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.196.pdf",dpi=600)

plt.show() 


file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz27")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.200$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.200.pdf",dpi=600)

plt.show() 


file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz28")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.203$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.203.pdf",dpi=600)

plt.show() 


file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/zz29")
om1=file_55[:,0] 
p1=file_55[:,1]
x=0
fig=plt.figure(constrained_layout=True)
plt.plot(om1,p1,'-*',label=r'$m^*=0.7m_N,M=0.207$')
plt.axhline(y=0, color='r', linestyle='-')
     
plt.xlabel('$(\\omega)$ $(km^{-1})$',fontsize = 16)
plt.ylabel('$\\Delta P$ (surface) $(km^{-2})$',fontsize = 16)
plt.xlim(0.0,1.87)
plt.ylim(-1000.0,1000.0)
plt.legend()
fig.savefig("m=0.7;M=0.207.pdf",dpi=600)

plt.show() 




