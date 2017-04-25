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
phi1,phi2,phi3,phi4,phi5,phi6 = np.loadtxt("Phi.dat", unpack=True)
dritt, x = np.loadtxt("etasplint.dat", unpack=True)
print len(x), len(phi1)
plt.ylabel(r"$\Phi$")
plt.xlabel(r"x")
#plt.xlim(1800,0)
plt.plot(x, phi1[:-1],label=r"$\Phi_1$", color=tableau20[2])
plt.plot(x, phi2[:-1],label=r"$\Phi_2$", color=tableau20[3])
plt.plot(x, phi3[:-1],label=r"$\Phi_3$", color=tableau20[4])
plt.plot(x, phi4[:-1],label=r"$\Phi_4$", color=tableau20[5])
plt.plot(x, phi5[:-1],label=r"$\Phi_5$", color=tableau20[6])
plt.plot(x, phi6[:-1],label=r"$\Phi_6$", color=tableau20[7])
plt.legend(fancybox=True,framealpha=0.5)
plt.show()
#fig.savefig('Xe.pdf', bbox_inches='tight',pad_inches=0.106)