import numpy as np
import math
import cmath
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib
import matplotlib.colors as mcolors


def gauss1D(x, mu, sigma):
    return np.exp(-(x-mu)**2/(2*sigma**2))

def FTgauss1D(p, mu, sigma):
    return sigma*np.exp(-1j*mu*p)*np.exp(-sigma**2*p**2/2)

def A(mu1, mu2, nu1, nu2, sigma, sign):
    #sign: plus =0, minus = 1
    normfactor = 2*np.pi*sigma**2*(1+(-1)**(sign)*np.exp(-((mu2-mu1)**2+(nu2-nu1)**2)/(4*sigma**2)))
    return 1/np.sqrt(normfactor)

def overlap1D(mu, nu, sigma1, sigma2):
    return sigma1*sigma2*np.sqrt(2*np.pi/(1+sigma1**2*sigma2**2))*np.exp(-(2j*mu*nu+nu**2*sigma1**2+mu**2*sigma2**2)/(2*(1+sigma1**2*sigma2**2)))

def psi_minus(x1, x2, mu1, mu2, nu1, nu2, sigma):
    return A(mu1, mu2, nu1, nu2, sigma, 1)*(gauss1D(x1,mu1,sigma)*gauss1D(x2,nu2,sigma) - gauss1D(x1,mu2,sigma)*gauss1D(x2,nu1, sigma))

def psi_plus(x1, x2, mu1, mu2, nu1, nu2, sigma):
    return A(mu1, mu2,nu1, nu2, sigma, 0)*(gauss1D(x1,mu1,sigma)*gauss1D(x2,nu2,sigma) + gauss1D(x1,mu2,sigma)*gauss1D(x2,nu1, sigma))

def phi_plus(x1, x2, mu1, mu2, nu1, nu2, sigma):
    return A(mu1, mu2,nu1, nu2, sigma, 0)*(gauss1D(x1,mu1,sigma)*gauss1D(x2,nu1,sigma) + gauss1D(x1,mu2,sigma)*gauss1D(x2,nu2, sigma))

def phi_minus(x1, x2, mu1, mu2, nu1, nu2, sigma):
    return A(mu1, mu2,nu1, nu2, sigma, 1)*(gauss1D(x1,mu1,sigma)*gauss1D(x2,nu1,sigma) - gauss1D(x1,mu2,sigma)*gauss1D(x2,nu2, sigma))

def FTpsi_minus(x1, x2, mu1, mu2, nu1, nu2, sigma):
    return A(mu1, mu2, nu1, nu2, sigma, 1)*(FTgauss1D(x1,mu1,sigma)*FTgauss1D(x2,nu2,sigma) - FTgauss1D(x1,mu2,sigma)*FTgauss1D(x2,nu1, sigma))

def AS_overlap_prob(Omega1, Omega2, Theta1, Theta2, sigma_w, tau1, tau2, pi1, pi2, sigma_t):
    """ function for overlap of psi- states only in conjugate bases """
    norm_factor = A(Omega1, Omega2, Theta1, Theta2, sigma_w,1)*A(tau1,tau2,pi1, pi2,sigma_t, 1)
    overlap_func = overlap1D(Omega1,tau1, sigma_w, sigma_t)*overlap1D(Theta2, pi2, sigma_w, sigma_t) + overlap1D(Omega2, tau2, sigma_w, sigma_t)*overlap1D(Theta1, pi1, sigma_w, sigma_t) - overlap1D(Omega2, tau1, sigma_w, sigma_t)*overlap1D(Theta1, pi2, sigma_w, sigma_t) - overlap1D(Omega1, tau2, sigma_w, sigma_t)*overlap1D(Theta2, pi1, sigma_w, sigma_t)
    return np.square(np.abs(norm_factor*overlap_func))



tau1 = 0 #ps
tau2 = 220 #ps
s_t = 17 #ps
offset = 2 * np.pi * (195.86)  # THz ~ 1550nm
# offset = 2 * np.pi * (195.86+0.019/2) # THz ~ 1550nm

# #paper Clementi 2023 Programmable frequency-bin quantum states in a nano-engineered silicon device
Omega1  = 2*np.pi*195.86 # THz
Omega2 = 2*np.pi*(195.86+0.019) # THz
s_w = 2*np.pi*0.0011 # THz (1.3GHz FWHM)

t_start = tau1-4*s_t  # ps
t_end = tau2+4*s_t  # ps
t_start2 = -400
t_end2 = 400
w_start = Omega1 - 4*s_w  # THz
w_end = Omega2 + 4*s_w  # THz
w_start2 = -2 * np.pi * 0.03  # THz
w_end2 = 2 * np.pi * 0.03  # THz
n_points = 1000

t1 = np.linspace(t_start, t_end, n_points)
t2 = np.linspace(t_start, t_end, n_points)
t1_1 = np.linspace(t_start2, t_end2, n_points)
t2_1 = np.linspace(t_start2, t_end2, n_points)
w1 = np.linspace(w_start2, w_end2, n_points)
w2 = np.linspace(w_start2, w_end2, n_points)
w1_1 = np.linspace(w_start, w_end, n_points)
w2_1 = np.linspace(w_start, w_end, n_points)

T1, T2 = np.meshgrid(t1, t2)
T1_1, T2_1 = np.meshgrid(t1_1, t2_1)
W1, W2 = np.meshgrid(w1, w2)
W1_1, W2_1 = np.meshgrid(w1_1, w2_1)

psi_minus1 = A(tau1, tau2, tau1, tau2, s_t, 1) * (
        gauss1D(T1, tau1, s_t) * gauss1D(T2, tau2, s_t) - gauss1D(T1, tau2, s_t) * gauss1D(T2, tau1, s_t))
psi_minus2 = psi_minus(T1, T2, tau1, tau2, tau1, tau2, s_t)
psi_minus3 = psi_minus(W1_1, W2_1, Omega1, Omega2, Omega1, Omega2, s_w)

psi_minusFT1 = FTpsi_minus(W1, W2, tau1, tau2, tau1, tau2, s_t)
psi_minusFT2 = FTpsi_minus(T1_1, T2_1, Omega1, Omega2, Omega1, Omega2, s_w)



fsize = 22  # fontsize for the figures
matplotlib.rcParams.update({'font.size': fsize})

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]

xvec = [0.10, 0.10, 0.8]
# Psi^-t in the time domain
fig = plt.figure(figsize=(8, 8))
#fig.patch.set_facecolor('#F9F9F9')
ax = fig.add_subplot(111)
bbox = ax.get_position()
new_bbox = (bbox.x0+xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
ax.set_position(new_bbox)
im1 = ax.imshow(psi_minus2, cmap='inferno', extent=[t_start, t_end, t_start, t_end],
                        interpolation='nearest', origin='lower')
ax.set_xlabel('Time of photon 1, $t_1$/ps')# $\sigma_t$ $\sigma_{\omega}$')
ax.set_ylabel('Time of photon 2, $t_2$/ps')
ax.set_xticks([0, 100, 200])
ax.set_xticklabels(['0', '100', '200'])
ax.set_yticks([0, 100, 200])
ax.set_yticklabels(['0', '100', '200'])
ax.set_title(r'$\Psi^-_t(t_1, t_2)$',fontweight='bold',pad=12)

axins = inset_axes(
    ax,
    width="5%",  # width: 5% of parent_bbox width
    height="100%",  # height: 50%
    loc="lower left",
    bbox_to_anchor=(1.05, 0., 1, 1),
    bbox_transform=ax.transAxes,
    borderpad=0,
)
fig.colorbar(im1, cax=axins,label='Normalised amplitude')

plt.show()

#Psi^-t in the freq domain
fig2 = plt.figure(figsize=(8, 8))
ax2 = fig2.add_subplot(111)
bbox = ax2.get_position()
new_bbox = (bbox.x0+xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
ax2.set_position(new_bbox)
im3 = ax2.imshow(np.abs(psi_minusFT1), cmap='inferno', extent=[(w_start2 + offset) / (2 * np.pi), (w_end2 + offset) / (2 * np.pi), (w_start2 + offset)/ (2 * np.pi), (w_end2 + offset) / (2 * np.pi)],
               interpolation='nearest', origin='lower')
ax2.set_xlabel(r'Freq. of photon 1, $(\omega_1/2\pi)$/THz')
ax2.set_ylabel(r'Freq. of photon 2, $(\omega_2/2\pi)$/THz')
ax2.set_xticks([195.84, 195.86, 195.88])
ax2.set_xticklabels(['195.84', '195.86', '195.88'])
ax2.set_yticks([195.84, 195.86, 195.88])
ax2.set_yticklabels(['195.84', '195.86', '195.88'])
ax2.set_title(r'$\Psi^-_t(\omega_1, \omega_2)$',fontweight='bold',pad=12)
axins = inset_axes(
    ax2,
    width="5%",  # width: 5% of parent_bbox width
    height="100%",  # height: 50%
    loc="lower left",
    bbox_to_anchor=(1.05, 0., 1, 1),
    bbox_transform=ax2.transAxes,
    borderpad=0,
)

fig.colorbar(im3, cax=axins,label='Absolute value of \nnormalised amplitude')

plt.show()

# Psi^-_f in the f domain
fig3 = plt.figure(figsize=(8, 8))
ax3 = fig3.add_subplot(111)
bbox = ax3.get_position()
new_bbox = (bbox.x0+xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
ax3.set_position(new_bbox)
im2 = ax3.imshow(psi_minus3, cmap='inferno',
                        extent=[w_start / (2 * np.pi), w_end / (2 * np.pi), w_start / (2 * np.pi),
                                w_end / (2 * np.pi)],
                        interpolation='nearest', origin='lower')
ax3.set_xlabel('Freq. of photon 1, $(\omega_1/2\pi)$/THz')
ax3.set_ylabel('Freq. of photon 2, $(\omega_2/2\pi)$/THz')
ax3.set_xticks([195.86, 195.87, 195.88])
ax3.set_xticklabels(['195.86', '195.87', '195.88'])
ax3.set_yticks([195.86, 195.87, 195.88])
ax3.set_yticklabels(['195.86', '195.87', '195.88'])
ax3.set_title(r'$\Psi^-_f(\omega_1, \omega_2)$',fontweight='bold',pad=12)

axins = inset_axes(
    ax3,
    width="5%",  # width: 5% of parent_bbox width
    height="100%",  # height: 50%
    loc="lower left",
    bbox_to_anchor=(1.05, 0., 1, 1),
    bbox_transform=ax3.transAxes,
    borderpad=0,
)
fig.colorbar(im2, cax=axins, ticks=[-50,-25,0,25,50], label='Normalised amplitude')

plt.show()

#Psi^-_f in the time domain
fig4 = plt.figure(figsize=(8, 8))
ax4 = fig4.add_subplot(111)
bbox = ax4.get_position()
new_bbox = (bbox.x0+xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
ax4.set_position(new_bbox)
im4 = ax4.imshow(np.abs(psi_minusFT2), cmap='inferno', extent=[t_start2, t_end2, t_start2, t_end2],
               interpolation='nearest', origin='lower')
ax4.set_xlabel('Time of photon 1, $t_1$/ps')
ax4.set_ylabel('Time of photon 2, $t_2$/ps')
ax4.set_xticks([-250, 0, 250])
ax4.set_xticklabels(['-250', '0', '250'])
ax4.set_yticks([-250, 0, 250])
ax4.set_yticklabels(['-250', '0', '250'])
ax4.set_title(r'$\Psi^-_f(t_1, t_2)$',fontweight='bold',pad=12)

axins = inset_axes(
    ax4,
    width="5%",  # width: 5% of parent_bbox width
    height="100%",  # height: 50%
    loc="lower left",
    bbox_to_anchor=(1.05, 0., 1, 1),
    bbox_transform=ax4.transAxes,
    borderpad=0,
)
cbar = fig.colorbar(im4, cax=axins, ticks=[0,0.001,0.002,0.003,0.004,0.005],  label='Absolute value of \nnormalised amplitude/10$^{-3}$')
cbar.ax.set_yticklabels(['0', '1', '2', '3', '4', '5'])  # vertically oriented colorbar

plt.show()




cmap = matplotlib.colormaps['inferno']
colors = cmap([0.3, 0.6, 0.9])
transparency = [1, 0.7, 0.55]
# colors = cmap([0.4, 0.7, 0.9])
# transparency = [1, 0.85, 0.7]

fsize = 20  # fontsize for the figures
matplotlib.rcParams.update({'font.size': fsize})
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]
#plt.rcParams["axes.unicode_minus"]
minus_sign = '\u2212'

xvec = [0.0, 0.0, 1]

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(projection='3d')
bbox = ax.get_position()
new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
ax.set_position(new_bbox)
ax.plot_surface(T1, T2, psi_minus2, cmap=plt.cm.get_cmap('PuOr', 256).reversed(), alpha=transparency[0], edgecolor='grey', lw=0.1, rstride=30, cstride=30)#, label=fr'$\sigma_\omega$ = {sigma_w[i]*1000} GHz')
ax.set_xlabel('Time of photon 1,\n $t_1$ (ps)',labelpad=23)
ax.set_ylabel('Time of photon 2,\n $t_2$ (ps)',labelpad=23)
ax.set_zlabel('Normalised \namplitude',labelpad=20)
ax.set_xticks([0, 220])
ax.set_xticklabels([r'$\tau_0$', r'$\tau_1$'])
ax.set_yticks([0, 220])
ax.set_yticklabels([r'$\tau_0$', r'$\tau_1$'])
ax.set_zticks([np.min(psi_minus2), 0, np.max(psi_minus2)])
ax.set_zticklabels([fr'${minus_sign}A_{{\tau}}$', '0', r'$A_{\tau}$'])
ax.tick_params(axis='z', pad=5)
ax.tick_params(axis='x', pad=0)
ax.tick_params(axis='y', pad=0)
ax.set_title(r'$\Psi^-_t(t_1, t_2)$',pad=-10) #,fontweight='bold'
ax.contour(T1, T2, psi_minus2, zdir='x', offset=-68, levels = [tau1, tau2], stride=0.01, cmap='PuOr')
ax.set_xlim([t_start, t_end])
ax.set_ylim([t_start, t_end])
ax.set_zlim([np.min(psi_minus2)-0.002, np.max(psi_minus2)+0.002])
#fig.fill_between(T1, T2, tau1-4*s_t, facecolor='lavender', interpolate=True)
#ax.set_zticks([0, 0.5, 1])
#ax.set_zticklabels(['0.0', '0.5', '1.0'])
# plt.legend(bbox_to_anchor=(-0.5, 0.5), loc='center left', fontsize=18)
#plt.legend(bbox_to_anchor=(1, 1.16), loc='upper right', fontsize=18)
ax.view_init(20, -50) #-70
plt.show()

print(np.max(psi_minus2))

cmap = matplotlib.colormaps['PuOr']
yvec = [0.0, 0.0, 0.6, 1.0]


fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(projection='3d')
bbox = ax.get_position()
new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
ax.set_position(new_bbox)
ax.plot_surface(W1_1, W2_1, psi_minus3, cmap=plt.cm.get_cmap('PuOr', 256).reversed(), alpha=transparency[0],edgecolor='grey', lw=0.1, rstride=30, cstride=30)#, label=fr'$\sigma_\omega$ = {sigma_w[i]*1000} GHz')
ax.set_xlabel('Freq. of photon 1,\n $\omega_1$ ($2\pi$ THz)',labelpad=23)
ax.set_ylabel('Freq. of photon 2,\n $\omega_2$ ($2\pi$ THz)',labelpad=23)
ax.set_zlabel('Normalised \namplitude',labelpad=20)
ax.set_xticks([Omega1, Omega2])
ax.set_xticklabels([r'$\Omega_0$', r'$\Omega_1$'])
ax.set_yticks([Omega1, Omega2])
ax.set_yticklabels([r'$\Omega_0$', r'$\Omega_1$'])
ax.set_zticks([np.min(psi_minus3), 0, np.max(psi_minus3)])
ax.set_zticklabels([f'${minus_sign}A_{{\Omega}}$', '0', r'$A_{\Omega}$'])
ax.tick_params(axis='z', pad=5)
ax.tick_params(axis='x', pad=0)
ax.tick_params(axis='y', pad=0)
ax.set_title(r'$\Psi^-_{f}(\omega_1, \omega_2)$',pad=-10) #,fontweight='bold'
ax.contour(W1_1, W2_1, psi_minus3, zdir='x', offset=Omega1 - 4*s_w, levels = [Omega1, Omega2], stride=1, cmap='PuOr')
ax.set_xlim([w_start, w_end])
ax.set_ylim([w_start, w_end])
ax.set_zlim([np.min(psi_minus3) - 3, np.max(psi_minus3) + 3])
#ax.set_zticks([0, 0.5, 1])
#ax.set_zticklabels(['0.0', '0.5', '1.0'])
# plt.legend(bbox_to_anchor=(-0.5, 0.5), loc='center left', fontsize=18)
#plt.legend(bbox_to_anchor=(1, 1.16), loc='upper right', fontsize=18)
ax.view_init(20, -50)
plt.show()

print(np.max(psi_minus3))


