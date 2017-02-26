import matplotlib.pyplot as plt 
import numpy as np 
from matplotlib import rc


# Plotting splint
file = open('splint.dat')
lines = file.readlines()

eta = np.zeros(len(lines))
x = np.zeros(len(lines))

i = 0
for line in lines:
    line = line.split()
    x[i] = float(line[0])
    eta[i] = float(line[1])
    i += 1
"""
# Plotting omegas

file = open('omegas.dat')
lines = file.readlines()

Om = np.zeros(len(lines))
Ob = np.zeros(len(lines))
Or = np.zeros(len(lines))
Ol = np.zeros(len(lines))

i = 0
for line in lines:
    line = line.split()
    Om[i] = float(line[0])
    Ob[i] = float(line[1])
    Or[i] = float(line[2])
    Ol[i] = float(line[3])
    i += 1
"""
# Plotting Hx and Hz
file = open('HxHz.dat')
lines = file.readlines()

Hx = np.zeros(len(lines))
z  = np.zeros(len(lines))
Hz = np.zeros(len(lines))
i=0
for line in lines:
	line  = line.split()
	Hx[i] = float(line[0])
	z[i]  = float(line[1])
	Hz[i] = float(line[2])
	i+=1


# These are the "Tableau 20" colors as RGB.    
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]  
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.) 


#--------------------Configuration------------------------------
fig = plt.figure()
# Remove the plot frame lines. They are unnecessary chartjunk.    
ax = plt.subplot(111)  
ax.yaxis.grid() 
#ax.xaxis.grid()  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(False)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(False)   
plt.xticks(fontsize=10)    
plt.yticks(fontsize=10)    
# Remove the tick marks; they are unnecessary with the tick lines we just plotted.    
plt.tick_params(axis="both", which="both", bottom="off", top="off",    
                labelbottom="on", left="off", right="off", labelleft="on")

plt.grid(b=True, which='minor', alpha=0.2)

#------------------Plotting-----------------
"""
# Plotting omegas
plt.title(r"/Omega")
plt.xlabel(r"Percentage")
plt.ylabel(r"x = log(a)")
plt.plot(x,Om, tableau20[2])
plt.plot(x,Ob,tableau20[4])
plt.plot(x,Or,tableau20[6])
plt.plot(x,Ol,tableau20[8])
fig.savefig('Omegas.pdf', bbox_inches='tight',pad_inches=0.106)
plt.show()
"""

# Plotting conformal time
plt.title(r"$\eta$",fontsize=12)
plt.xlabel(r"$\eta$",fontsize=14)
plt.ylabel(r"x = log(a)",fontsize=14)
plt.plot(x,eta, color=tableau20[4])
fig.savefig('Hx.pdf', bbox_inches='tight',pad_inches=0.106)
plt.show()

# Plotting Hx
plt.title(r"H(x)",fontsize=12)
plt.xlabel(r"H",fontsize=14)
plt.ylabel(r"x = log(a)",fontsize=14)
plt.plot(x,Hx, color=tableau20[4])
fig.savefig('Hx.pdf', bbox_inches='tight',pad_inches=0.106)
plt.show()

# Plotting Hz
plt.title(r"Hz")
plt.xlabel(r"H")
plt.ylabel(r"z")
plt.plot(z, Hz, color=tableau20[4])
fig.savefig('Hz.pdf', bbox_inches='tight',pad_inches=0.106)
plt.show()

