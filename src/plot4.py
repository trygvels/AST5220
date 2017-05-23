
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
"""
plt.xkcd()

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
#plt.xticks([])
#plt.yticks([])
ax.set_ylim([0, 6000])
ax.set_xlim([0, 1200])

plt.annotate('ITS CORRECT,\nIF THE UNIVERSE\nIS COMPLETELY DIFFERENT.',
    xy=(640, 3000), arrowprops=dict(arrowstyle='->'), xytext=(250, 1000))
lhi,clhi = np.loadtxt("data/C_lM1.dat", unpack=True)
clhi = clhi/max(clhi[100:])*5775

plt.plot(clhi)

plt.xlabel('l')
plt.ylabel(r'l(l+1)C_l/2$\pi$')
fig.savefig('figs/xkcdcl.pdf', bbox_inches='tight',pad_inches=0.106)
plt.show()

"""
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


#--------------------CL change h --------------------
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

planck1  = np.loadtxt("data/COM_PowerSpect_CMB-TT-loL-full_R2.02.txt",unpack=True,skiprows=3)
planck2  = np.loadtxt("data/COM_PowerSpect_CMB-TT-hiL-full_R2.02.txt",unpack=True,skiprows=3)

planck_l1 = planck1[0]
planck_l2 = planck2[0]

C_llow = planck1[1]
C_lhi  = planck2[1]

error1 = planck1[2]
error2 = planck2[2]

Clplanck = np.hstack([C_llow,C_lhi])
planck_l = np.hstack([planck_l1,planck_l2])
error = np.hstack([error1,error2])



lhi,clhi = np.loadtxt("data/C_lM1.dat", unpack=True)
lhi,clhiH66 = np.loadtxt("data/C_lh67.dat", unpack=True)
lhi,clhiH74 = np.loadtxt("data/C_lh74.dat", unpack=True)

plt.xlabel(r'l')
plt.ylabel(r'$l(l+1)C_l/2\pi$')
plt.xlim(0,1200)
plt.ylim(0,6000)
K = 5775
clhi = clhi/max(clhi[100:])*K
clhiH66 = clhiH66/max(clhiH66)*K
clhiH74 = clhiH74/max(clhiH74)*K
plt.plot(lhi, clhi,label=r"$C_l$ Default", color=tableau20[0],zorder=4)
plt.plot(lhi, clhiH66,label=r"$C_l$ h = 67", color=tableau20[2],zorder=2)
plt.plot(lhi, clhiH74,label=r"$C_l$ h = 74", color=tableau20[4],zorder=3)
plt.errorbar(planck_l,Clplanck,yerr=error,label='Planck data',color=tableau20[15],zorder=1)#,fmt='-o')
plt.legend(fancybox=True,framealpha=0.5,loc=1)
fig.savefig('figs/clvarh.pdf', bbox_inches='tight',pad_inches=0.106)
plt.show()
#--------------------CL change b--------------------
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

planck1  = np.loadtxt("data/COM_PowerSpect_CMB-TT-loL-full_R2.02.txt",unpack=True,skiprows=3)
planck2  = np.loadtxt("data/COM_PowerSpect_CMB-TT-hiL-full_R2.02.txt",unpack=True,skiprows=3)

planck_l1 = planck1[0]
planck_l2 = planck2[0]

C_llow = planck1[1]
C_lhi  = planck2[1]

error1 = planck1[2]
error2 = planck2[2]

Clplanck = np.hstack([C_llow,C_lhi])
planck_l = np.hstack([planck_l1,planck_l2])
error = np.hstack([error1,error2])


plt.errorbar(planck_l,Clplanck,yerr=error,label='Planck data',color=tableau20[15],zorder=1)#,fmt='-o')

lhi,clhi = np.loadtxt("data/C_lM1.dat", unpack=True)
lhi,clhib42 = np.loadtxt("data/C_lb42.dat", unpack=True)
lhi,clhib50 = np.loadtxt("data/C_lb50.dat", unpack=True)

plt.xlabel(r'l')
plt.ylabel(r'$l(l+1)C_l/2\pi$')
plt.xlim(0,1200)
plt.ylim(0,6000)
K = 5775
clhi = clhi/max(clhi[100:])*K
clhib42 = clhib42/max(clhib42)*K
clhib50 = clhib50/max(clhib50)*K
plt.plot(lhi, clhi,label=r"$C_l$ Default", color=tableau20[0],zorder=4)
plt.plot(lhi, clhib42,label=r"$C_l$ b = 0.042", color=tableau20[2],zorder=2)
plt.plot(lhi, clhib50,label=r"$C_l$ b = 0.050", color=tableau20[4],zorder=3)
plt.legend(fancybox=True,framealpha=0.5,loc=1)
fig.savefig('figs/clvarb.pdf', bbox_inches='tight',pad_inches=0.106)
plt.show()
#--------------------CL Change m--------------------
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

planck1  = np.loadtxt("data/COM_PowerSpect_CMB-TT-loL-full_R2.02.txt",unpack=True,skiprows=3)
planck2  = np.loadtxt("data/COM_PowerSpect_CMB-TT-hiL-full_R2.02.txt",unpack=True,skiprows=3)

planck_l1 = planck1[0]
planck_l2 = planck2[0]

C_llow = planck1[1]
C_lhi  = planck2[1]

error1 = planck1[2]
error2 = planck2[2]

Clplanck = np.hstack([C_llow,C_lhi])
planck_l = np.hstack([planck_l1,planck_l2])
error = np.hstack([error1,error2])


plt.errorbar(planck_l,Clplanck,yerr=error,label='Planck data',color=tableau20[15],zorder=1)#,fmt='-o')

lhi,clhi = np.loadtxt("data/C_lM1.dat", unpack=True)
lhi,clhim200 = np.loadtxt("data/C_lm200.dat", unpack=True)
lhi,clhim248 = np.loadtxt("data/C_lm248.dat", unpack=True)

plt.xlabel(r'l')
plt.ylabel(r'$l(l+1)C_l/2\pi$')
plt.xlim(0,1200)
plt.ylim(0,6000)
K = 5775
clhi = clhi/max(clhi[100:])*K
clhim200 = clhim200/max(clhim200)*K
clhim248 = clhim248/max(clhim248)*K
plt.plot(lhi, clhi,label=r"$C_l$ Default", color=tableau20[0],zorder=4)
plt.plot(lhi, clhim200,label=r"$C_l$ m = 200", color=tableau20[2],zorder=2)
plt.plot(lhi, clhim248,label=r"$C_l$ m = 248", color=tableau20[4],zorder=3)
plt.legend(fancybox=True,framealpha=0.5,loc=1)
fig.savefig('figs/clvarm.pdf', bbox_inches='tight',pad_inches=0.106)
plt.show()

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

khi,thetahi1 = np.loadtxt("data/integrand1.dat", unpack=True)
thetahi2 = np.loadtxt("data/integrand2.dat", unpack=True)
thetahi3 = np.loadtxt("data/integrand3.dat", unpack=True)
thetahi4 = np.loadtxt("data/integrand4.dat", unpack=True)
thetahi5 = np.loadtxt("data/integrand5.dat", unpack=True)
thetahi6 = np.loadtxt("data/integrand6.dat", unpack=True)

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
fig.savefig('figs/integrands.pdf', bbox_inches='tight',pad_inches=0.106)
plt.show()
