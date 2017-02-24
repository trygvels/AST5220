import matplotlib.pyplot as plt 
import numpy as np 
file = open('splint.dat')

lines = file.readlines()

eta = np.zeros(len(lines))
x = np.zeros(len(lines))

i= 0
for line in lines:
    line = line.split()
    eta[i] = float(line[1])
    x[i] = float(line[0])
    i += 1

file = open('eta.dat')

lines = file.readlines()

eta1 = np.zeros(len(lines))
x1 = np.zeros(len(lines))

i= 0
for line in lines:
    line = line.split()
    eta1[i] = float(line[1])
    x1[i] = float(line[0])
    i += 1



plt.plot(x,eta,x1,eta1)
plt.show()
