import matplotlib.pyplot as plt
  
x = []
y = []
for line in open('pmode.out', 'r'):
    lines = [i for i in line.split()]
    x.append(lines[2])
    y.append(lines[3])
      
plt.title("P Mode frequencies")
plt.xlabel('Mass ($M_{\odot}$)')
plt.ylabel('F (kHz)')
plt.yticks()
plt.plot(x, y, marker = 'o', c = 'g')
  
plt.show()


####################################


