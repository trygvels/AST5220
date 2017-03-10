import matplotlib.pyplot as plt
import numpy as np


# Plotting splint
file = open('etasplint.dat')
lines = file.readlines()

etas = np.zeros(len(lines)) # length - 500
x_t = np.zeros(len(lines)) # length - 500

i = 0
for line in lines:
    line = line.split()
    x_t[i] = float(line[1])
    etas[i] = float(line[0])
    i += 1

# Plotting splint
file = open('eta.dat')
lines = file.readlines()

eta = np.zeros(len(lines))
x = np.zeros(len(lines))

i = 0
for line in lines:
    line = line.split()
    x[i] = float(line[1]) #length - 1000
    eta[i] = float(line[0]) #length - 1000
    i += 1

# Plotting omegas

file = open('omegas1.dat')
lines = file.readlines()

Om = np.zeros(len(lines))
Ob = np.zeros(len(lines))

i = 0
for line in lines:
    line = line.split()
    Om[i] = float(line[0])
    Ob[i] = float(line[1])
    i += 1

file = open('omegas2.dat')
lines = file.readlines()

Or = np.zeros(len(lines))
Ol = np.zeros(len(lines))
i = 0
for line in lines:
    line = line.split()
    Or[i] = float(line[0])
    Ol[i] = float(line[1])
    i += 1

# Plotting Hx and Hz
file = open('HxHz.dat')
lines = file.readlines()

Hx = np.zeros(len(lines)) #Length 1000
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
t = np.linspace(0,13.77,1000)
# Plotting omegas
plt.title(r"Evolution of $\Omega$")
plt.ylabel(r"$\Omega$")
plt.xlabel(r"Gyrs")
plt.xlim(0,13.77)
plt.plot(np.exp(x)*13.77,Om,label=r"$\Omega_m$", color=tableau20[2])
plt.plot(np.exp(x)*13.77,Ob,label=r"$\Omega_b$", color=tableau20[4])
plt.plot(np.exp(x)*13.77,Or,label=r"$\Omega_r$", color=tableau20[6])
plt.plot(np.exp(x)*13.77,Ol,label=r"$\Omega_{\lambda}$", color=tableau20[8])
plt.legend(loc=6)
plt.show()
fig.savefig('Omegas.pdf', bbox_inches='tight',pad_inches=0.106)

# Plotting omegas
plt.title(r"Evolution of $\Omega$")
plt.ylabel(r"$\Omega$")
plt.xlabel(r"x = log(a)")
plt.xlim(-24,0)
plt.plot(x,Om,label=r"$\Omega_m$", color=tableau20[2])
plt.plot(x,Ob,label=r"$\Omega_b$", color=tableau20[4])
plt.plot(x,Or,label=r"$\Omega_r$", color=tableau20[6])
plt.plot(x,Ol,label=r"$\Omega_{\lambda}$", color=tableau20[8])
plt.legend(loc=6)
plt.show()
fig.savefig('OmegasLog.pdf', bbox_inches='tight',pad_inches=0.106)
"""
MPc= 3.085e22
# Plotting conformal time
plt.title(r"Conformal time $\eta$ interpolated",fontsize=12)
plt.ylabel(r"$\eta$ [MPc]",fontsize=14)
plt.xlabel(r"x = log(a)",fontsize=14)
plt.xlim(-24,0)
plt.semilogy(x_t,etas/MPc, "-o", label="Interpolated section", color=tableau20[4])
plt.semilogy(x,eta/MPc, label="Conformal time", color=tableau20[6])
plt.legend(loc=7)
plt.show()
fig.savefig('Hx.pdf', bbox_inches='tight',pad_inches=0.106)

# Plotting Hx
plt.title(r"Hubble factor H(x) as a function of log(a)",fontsize=12)
plt.ylabel(r"H(x)",fontsize=14)
plt.xlabel(r"x = log(a)",fontsize=14)
plt.xlim(-24,0)
plt.semilogy(x,Hx, color=tableau20[4])
plt.legend(loc=7)
plt.show()
fig.savefig('Hx.pdf', bbox_inches='tight',pad_inches=0.106)

# Plotting Hz
plt.title(r"Hubble factor H(z) as a function of redshift z")
plt.ylabel(r"H(z)")
plt.xlabel(r"z - redshift")
plt.semilogy(z, Hz, color=tableau20[4])
plt.legend(loc=7)
plt.show()
fig.savefig('Hz.pdf', bbox_inches='tight',pad_inches=0.106)
"""
