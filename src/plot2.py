import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

#-------------------- Colors --------------------
# These are the "Tableau 20" colors as RGB.
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

#--------------------Configuration--------------------
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
x, Xe = np.loadtxt("Xe.dat", unpack=True)

# Plotting X_e
plt.title(r"Electron fraction $X_e$")
plt.ylabel(r"$X_e$")
plt.xlabel(r"z")
plt.xlim(1800,0)
z = 1/np.exp(x)-1
plt.semilogy(z, Xe,label=r"$X_e$", color=tableau20[2])
plt.legend(fancybox=True,framealpha=0.5)
plt.show()
#fig.savefig('Xe.pdf', bbox_inches='tight',pad_inches=0.106)

#--------------------Configuration--------------------

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


tau, dtau,tau2 = np.loadtxt("tau.dat", unpack=True)
# Plotting Tau
plt.title(r"Optical depth $\tau$ and absolute value of derivatives")
plt.ylabel(r"$\tau$, $|\tau'|$, $|\tau''|$")
plt.xlabel(r"z")
plt.ylim(1e-8,1e4)
plt.xlim(-10,0)
plt.semilogy(x, abs(tau),  label=r"$\tau$",color=tableau20[8])
plt.semilogy(x, abs(dtau), label=r"$|\tau'|$",color=tableau20[2])
plt.semilogy(x, abs(tau2), label=r"$|\tau''|$",color=tableau20[4])
plt.legend(fancybox=True,framealpha=0.5)
plt.show()
#fig.savefig('tau.pdf', bbox_inches='tight',pad_inches=0.106)

#--------------------Configuration--------------------
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

g,dg,g2 = np.loadtxt("g.dat", unpack=True)
# Plotting g
plt.title(r"Visibility function $\tilde{g}$ and scaled derivatives")
plt.ylabel(r"$\tilde{g}$, $\tilde{g}'/10$, $\tilde{g}''/300$")
plt.xlabel(r"z")
plt.xlim(-7.4,-6.4)
plt.plot(x, g,label=r"$\tilde{g}$", color=tableau20[8])
plt.plot(x, dg/10.,label=r"$\tilde{g}'$",color=tableau20[4])
plt.plot(x, g2/300., label=r"$\tilde{g}''$", color=tableau20[2])
plt.legend(fancybox=True,framealpha=0.5)
plt.show()
#fig.savefig('g.pdf', bbox_inches='tight',pad_inches=0.106)
