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

dritt, x = np.loadtxt("etasplint.dat", unpack=True) # x-grid

#--------------------PHI--------------------
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


phi1,phi2,phi3,phi4,phi5,phi6 = np.loadtxt("Phi.dat", unpack=True)
plt.ylabel(r"$\Phi$")
plt.xlabel(r"x")
plt.xlim(x[0],x[-1])

plt.plot(x, phi1,label=r"$\Phi, k = 1$", color=tableau20[0])
plt.plot(x, phi2,label=r"$\Phi, k = 5$", color=tableau20[2])
plt.plot(x, phi3,label=r"$\Phi, k = 10$", color=tableau20[4])
plt.plot(x, phi4,label=r"$\Phi, k = 40$", color=tableau20[6])
plt.plot(x, phi5,label=r"$\Phi, k = 60$", color=tableau20[16])
plt.plot(x, phi6,label=r"$\Phi, k = 100$", color=tableau20[18])
plt.legend(fancybox=True,framealpha=0.5,loc=3)
fig.savefig('phi.pdf', bbox_inches='tight',pad_inches=0.106)
plt.show()

#--------------------PSI--------------------
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

psi1,psi2,psi3,psi4,psi5,psi6 = np.loadtxt("Psi.dat", unpack=True)
plt.ylabel(r"$\Psi$")
plt.xlabel(r"x")
plt.xlim(x[0],x[-1])

plt.plot(x, psi1,label=r"$\Psi, k = 1$", color=tableau20[0])
plt.plot(x, psi2,label=r"$\Psi, k = 5$", color=tableau20[2])
plt.plot(x, psi3,label=r"$\Psi, k = 10$", color=tableau20[4])
plt.plot(x, psi4,label=r"$\Psi, k = 40$", color=tableau20[6])
plt.plot(x, psi5,label=r"$\Psi, k = 60$", color=tableau20[16])
plt.plot(x, psi6,label=r"$\Psi, k = 100$", color=tableau20[18])
plt.legend(fancybox=True,framealpha=0.5,loc=3)
fig.savefig('psi.pdf', bbox_inches='tight',pad_inches=0.106)
plt.show()

#--------------------DELTA--------------------
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

delta1,delta2,delta3,delta4,delta5,delta6 = np.loadtxt("delta.dat", unpack=True)
plt.ylabel(r"$\delta$")
plt.xlabel(r"x")
plt.xlim(x[0],x[-1])

plt.plot(x, delta1,label=r"$\delta, k = 1$", color=tableau20[0])
plt.plot(x, delta2,label=r"$\delta, k = 5$", color=tableau20[2])
plt.plot(x, delta3,label=r"$\delta, k = 10$", color=tableau20[4])
plt.plot(x, delta4,label=r"$\delta, k = 40$", color=tableau20[6])
plt.plot(x, delta5,label=r"$\delta, k = 60$", color=tableau20[16])
plt.plot(x, delta6,label=r"$\delta, k = 100$", color=tableau20[18])
plt.yscale('symlog')
plt.legend(fancybox=True,framealpha=0.5,loc=3)
fig.savefig('delta.pdf', bbox_inches='tight',pad_inches=0.106)
plt.show()

#--------------------DELTA_B--------------------
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

delta_b1,delta_b2,delta_b3,delta_b4,delta_b5,delta_b6 = np.loadtxt("delta_b.dat", unpack=True)
plt.ylabel(r"$\delta_b$")
plt.xlabel(r"x")
plt.xlim(x[0],x[-1])

plt.plot(x, delta_b1,label=r"$\delta_b, k = 1$", color=tableau20[0])
plt.plot(x, delta_b2,label=r"$\delta_b, k = 5$", color=tableau20[2])
plt.plot(x, delta_b3,label=r"$\delta_b, k = 10$", color=tableau20[4])
plt.plot(x, delta_b4,label=r"$\delta_b, k = 40$", color=tableau20[6])
plt.plot(x, delta_b5,label=r"$\delta_b, k = 60$", color=tableau20[16])
plt.plot(x, delta_b6,label=r"$\delta_b, k = 100$", color=tableau20[18])
plt.yscale('symlog')
plt.legend(fancybox=True,framealpha=0.5,loc=3)
fig.savefig('deltab.pdf', bbox_inches='tight',pad_inches=0.106)
plt.show()

#--------------------v--------------------
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

v1,v2,v3,v4,v5,v6 = np.loadtxt("v.dat", unpack=True)
plt.ylabel(r"$v$")
plt.xlabel(r"x")
plt.xlim(x[0],x[-1])

plt.plot(x, v1,label=r"$v, k = 1$", color=tableau20[0])
plt.plot(x, v2,label=r"$v, k = 5$", color=tableau20[2])
plt.plot(x, v3,label=r"$v, k = 10$", color=tableau20[4])
plt.plot(x, v4,label=r"$v, k = 40$", color=tableau20[6])
plt.plot(x, v5,label=r"$v, k = 60$", color=tableau20[16])
plt.plot(x, v6,label=r"$v, k = 100$", color=tableau20[18])
plt.yscale('symlog')
plt.legend(fancybox=True,framealpha=0.5,loc=3)
fig.savefig('v.pdf', bbox_inches='tight',pad_inches=0.106)
plt.show()

#--------------------v_b--------------------
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

vb1,vb2,vb3,vb4,vb5,vb6 = np.loadtxt("v_b.dat", unpack=True)
plt.ylabel(r"$v_b$")
plt.xlabel(r"x")
plt.xlim(x[0],x[-1])
plt.plot(x, vb1,label=r"$v_b, k = 1$", color=tableau20[0])
plt.plot(x, vb2,label=r"$v_b, k = 5$", color=tableau20[2])
plt.plot(x, vb3,label=r"$v_b, k = 10$", color=tableau20[4])
plt.plot(x, vb4,label=r"$v_b, k = 40$", color=tableau20[6])
plt.plot(x, vb5,label=r"$v_b, k = 60$", color=tableau20[16])
plt.plot(x, vb6,label=r"$v_b, k = 100$", color=tableau20[18])
plt.yscale('symlog')
plt.legend(fancybox=True,framealpha=0.5,loc=3)
fig.savefig('vb.pdf', bbox_inches='tight',pad_inches=0.106)
plt.show()

#--------------------theta0--------------------
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

theta1,theta2,theta3,theta4,theta5,theta6 = np.loadtxt("Theta0.dat", unpack=True)
plt.ylabel(r"$\theta_0$")
plt.xlabel(r"x")
#plt.xlim(1800,0)
plt.plot(x, theta1,label=r"$\theta_0, k = 1$", color=tableau20[0])
plt.plot(x, theta2,label=r"$\theta_0, k = 5$", color=tableau20[2])
plt.plot(x, theta3,label=r"$\theta_0, k = 10$", color=tableau20[4])
plt.plot(x, theta4,label=r"$\theta_0, k = 40$", color=tableau20[6])
plt.plot(x, theta5,label=r"$\theta_0, k = 60$", color=tableau20[16])
plt.plot(x, theta6,label=r"$\theta_0, k = 100$", color=tableau20[18])
plt.xlim(x[0],x[-1])
plt.yscale('symlog')
plt.legend(fancybox=True,framealpha=0.5,loc=3)
fig.savefig('theta0.pdf', bbox_inches='tight',pad_inches=0.106)
plt.show()
