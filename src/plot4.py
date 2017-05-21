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


#--------------------CL--------------------
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

planck1  = np.loadtxt("COM_PowerSpect_CMB-TT-loL-full_R2.02.txt",unpack=True,skiprows=3)
planck2  = np.loadtxt("COM_PowerSpect_CMB-TT-hiL-full_R2.02.txt",unpack=True,skiprows=3)
planck_l1 = planck1[0]
planck_l2 = planck2[0]

C_llow = planck1[1]
C_lhi  = planck2[1]

error1 = planck1[2]
error2 = planck2[2]

Clplanck = np.hstack([C_llow,C_lhi])
planck_l = np.hstack([planck_l1,planck_l2])
error = np.hstack([error1,error2])


plt.errorbar(planck_l,Clplanck,yerr=error,label='Planck data',color='grey')#,fmt='-o')

lhi,clhi = np.loadtxt("C_l.dat", unpack=True)
lhi,clhiH66 = np.loadtxt("C_lH66.dat", unpack=True)
lhi,clhiH74 = np.loadtxt("C_lH74.dat", unpack=True)
lhi,clhib42 = np.loadtxt("C_lb42.dat", unpack=True)
lhi,clhib50 = np.loadtxt("C_lb50.dat", unpack=True)
lhi,clhim200 = np.loadtxt("C_lm200.dat", unpack=True)
lhi,clhim248 = np.loadtxt("C_lm248.dat", unpack=True)

plt.xlabel(r'l')
plt.ylabel(r'$l(l+1)C_l/2\pi$')
plt.xlim(min(lhi),1200)
plt.ylim(0,6000)
clhi = clhi/max(clhi)*5775
clhiH66 = clhiH66/max(clhiH66)*5775
clhiH74 = clhiH74/max(clhiH74)*5775
clhib42 = clhib42/max(clhib42)*5775
clhib50 = clhib50/max(clhib50)*5775
clhim200 = clhim200/max(clhim200)*5775
clhim248 = clhim248/max(clhim248)*5775
plt.plot(lhi, clhi,label=r"$C_l default$", color=tableau20[0])
plt.plot(lhi, clhiH66,label=r"$C_l h = 66$", color=tableau20[1])
plt.plot(lhi, clhiH74,label=r"$C_l h = 74$", color=tableau20[2])
plt.plot(lhi, clhib42,label=r"$C_l b = 42$", color=tableau20[3])
plt.plot(lhi, clhib50,label=r"$C_l b = 50$", color=tableau20[4])
plt.plot(lhi, clhim200,label=r"$C_l m = 200$", color=tableau20[5])
plt.plot(lhi, clhim248,label=r"$C_l m = 248$", color=tableau20[6])
plt.legend(fancybox=True,framealpha=0.5,loc=7)
plt.show()
"""
#----------------INTEGRANDS----------------
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

khi,thetahi1 = np.loadtxt("integrand1.dat", unpack=True)
thetahi2 = np.loadtxt("integrand2.dat", unpack=True)
thetahi3 = np.loadtxt("integrand3.dat", unpack=True)
thetahi4 = np.loadtxt("integrand4.dat", unpack=True)
thetahi5 = np.loadtxt("integrand5.dat", unpack=True)
thetahi6 = np.loadtxt("integrand6.dat", unpack=True)

plt.xlabel(r'k*c/H_0')
plt.ylabel(r'$l(l+1)\Theta^2/(ck/H_0)$')
plt.xlim([0,500])
plt.ylim(0,0.08)
plt.plot(khi, thetahi1,label=r"$\Theta, l = 2$", color=tableau20[0])
plt.plot(khi, thetahi2,label=r"$\Theta, l = 50$", color=tableau20[1])
plt.plot(khi, thetahi3,label=r"$\Theta, l = 200$", color=tableau20[2])
plt.plot(khi, thetahi4,label=r"$\Theta, l = 500$", color=tableau20[3])
plt.plot(khi, thetahi5,label=r"$\Theta, l = 800$", color=tableau20[4])
plt.plot(khi, thetahi6,label=r"$\Theta, l = 1200$", color=tableau20[5])
plt.legend(fancybox=True,framealpha=0.5,loc=7)
plt.show()
"""
