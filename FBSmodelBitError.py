import numpy as np
import math
import cmath
import matplotlib.pyplot as plt
from scipy.special import erf
import matplotlib as mpl
import pandas as pd
import scipy.interpolate as interpolate
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#the definition of erf(x) used in the python package is 2/sqrt(pi)*integrate(exp(-z^2))z=0 to z=x, x can be a real number

def A(tau1, tau2, sigma):
    normfactor = 2*(1 - np.exp(-((tau2-tau1)**2)/(2*sigma**2)))
    return 1/np.sqrt(normfactor)

# def overlapshift(Omega, tau, sigma_w, sigma_t, shift, start, end):
#     coeff = np.sqrt(sigma_t*sigma_w/(2*(1+sigma_w**2*sigma_t**2)))
#     phase = np.exp(-1j*(Omega+shift)*tau/(1+sigma_w**2*sigma_t**2))
#     shiftdecay = np.exp(-sigma_t**2*shift**2*(1-sigma_w**2 + sigma_t**2*sigma_w**2)/(2*(1+sigma_w**2*sigma_t**2)))
#     bindecay = np.exp(-(tau**2*sigma_w**2 + Omega*2*sigma_t**2)/(2*(1+sigma_w**2*sigma_t**2)))
#     mixdecay = np.exp(-Omega*shift*sigma_t**2/(1+sigma_w**2*sigma_t**2))
#     scaling = np.sqrt((1+sigma_t**2*sigma_w**2)/(2*sigma_w**2))
#     offset = (1j*sigma_w**2*tau + shift*sigma_w**2*sigma_t**2 - Omega)/(1+sigma_w**2*sigma_t**2)
#     endval = erf(scaling * end + scaling*offset)
#     startval = erf(scaling * start + scaling*offset)
#     return coeff*phase*shiftdecay*bindecay*mixdecay*(endval - startval)

def overlapshift(Omega, tau, sigma_w, sigma_t, shift, start, end):
    coeff = np.sqrt(sigma_t*sigma_w/(2*(1+sigma_w**2*sigma_t**2)))
    phase = np.exp(-1j*(Omega+shift)*tau/(1+sigma_w**2*sigma_t**2))
    shiftdecay = np.exp(-sigma_t**2*shift**2/(2*(1+sigma_w**2*sigma_t**2)))
    bindecay = np.exp(-(tau**2*sigma_w**2 + Omega**2*sigma_t**2)/(2*(1+sigma_w**2*sigma_t**2)))
    mixdecay = np.exp(-Omega*shift*sigma_t**2/(1+sigma_w**2*sigma_t**2))
    scaling = np.sqrt((1+sigma_t**2*sigma_w**2)/(2*sigma_w**2))
    offset = (1j*sigma_w**2*tau + shift*sigma_w**2*sigma_t**2 - Omega)/(1+sigma_w**2*sigma_t**2)
    endval = erf(scaling * end + scaling*offset)
    startval = erf(scaling * start + scaling*offset)
    print(Omega, tau)
    print(coeff, phase, shiftdecay, bindecay, mixdecay)
    print(endval)
    print(startval)
    return coeff*phase*shiftdecay*bindecay*mixdecay*(endval - startval)




def Uoverlap1D(Omega1,tau2, sigma_w,sigma_t, Omega, epsilon, mu, theta=0, phi=0):
    print(Omega1, tau2)
    f1 = overlapshift(Omega1, tau2, sigma_w, sigma_t, 0, -math.inf, Omega-epsilon/2)
    f2 = np.cos(theta)*overlapshift(Omega1, tau2, sigma_w, sigma_t, 0, Omega-epsilon/2, Omega+epsilon/2) - np.exp(-1j*phi)*np.sin(theta)*overlapshift(Omega1, tau2, sigma_w, sigma_t, mu, Omega-epsilon/2, Omega+epsilon/2)
    f3 = overlapshift(Omega1, tau2, sigma_w, sigma_t, 0, Omega + epsilon/2, Omega + mu - epsilon/2)
    f4 = np.cos(theta)*overlapshift(Omega1, tau2, sigma_w, sigma_t, 0, Omega+mu-epsilon/2, Omega+mu+epsilon/2) + np.exp(1j*phi)*np.sin(theta)*overlapshift(Omega1, tau2, sigma_w, sigma_t, -mu, Omega+mu-epsilon/2, Omega+mu+epsilon/2)
    f5 = overlapshift(Omega1, tau2, sigma_w, sigma_t, 0, Omega+mu+epsilon/2, math.inf)
    return f1 + f2 + f3 + f4 + f5

def Uoverlap2D(tau1,tau2, sigma_t, Omega1, Omega2, sigma_w, Omega, epsilon, mu, theta=0, phi=0):
    return 2*A(tau1, tau2, sigma_t)*A(Omega1, Omega2, sigma_w)*(Uoverlap1D(Omega1, tau1, sigma_w, sigma_t, Omega, epsilon, mu, theta, phi)*Uoverlap1D(Omega2, tau2, sigma_w, sigma_t, Omega, epsilon, mu, theta, phi)-Uoverlap1D(Omega1, tau2, sigma_w, sigma_t, Omega, epsilon, mu, theta, phi)*Uoverlap1D(Omega2, tau1, sigma_w, sigma_t, Omega, epsilon, mu, theta, phi))


def Uprob2D(tau1,tau2, sigma_t, Omega1, Omega2, sigma_w, Omega, epsilon, mu, theta=0, phi=0):
    return np.abs(Uoverlap2D(tau1,tau2, sigma_t, Omega1, Omega2, sigma_w, Omega, epsilon, mu, theta, phi))**2

tau1 = 0
tau2 = 220  # ps
# sigma_t = 0.01876
sigma_t = 17  # ps

Omega1 = 0
Omega2 = 2 * np.pi * 0.019
# Omega2 = 2*np.pi*(1+(2*np.pi*0.0005)**2*sigma_t**2)/(tau2-tau1) + Omega1
# print(Omega2)
sigma_w = 2*np.pi*np.asarray([0.0001, 0.0005, 0.0011])


# Omega = -0.01*2*np.pi
# width = 2*np.pi*np.asarray([0.0005, 0.001])
# gap = 0.019 * 2 * np.pi  # THz #0.019*2*np.pi

#CALCULATING THE BIT ERROR RATE!

Omega = 0
width = 6*sigma_w
gap = 0.019 * 2 * np.pi  # THz #0.019*2*np.pi


theta = np.linspace(0, 2*np.pi, 251)
phi = np.linspace(0, 2*np.pi, 251)


# theta=np.linspace(0, 2*np.pi, 1200)
# phi=np.linspace(0, 2*np.pi,1200)
# delta=0

Theta, Phi = np.meshgrid(theta, phi)
probabilities = np.zeros((len(theta), len(phi), len(width)))

print('bit error, zero angles')
print(Uprob2D(tau1, tau2, sigma_t, Omega1, Omega2, sigma_w, Omega, width[1], gap, 0, 0))

#print(Uprob2D(tau1, tau2, sigma_t, Omega1, Omega2, sigma_, Omega1, width[0], gap, np.pi, 0, delta))

for i in range(len(width)):
    probabilities[:, :, i] = Uprob2D(tau1, tau2, sigma_t, Omega1, Omega2, sigma_w[i], Omega, width[i], gap, Theta, Phi)



cmap = mpl.colormaps['inferno']
colors = cmap([0.9, 0.7, 0.4])
transparency = [0.7, 0.85, 1]


fsize = 22  # fontsize for the figures
mpl.rcParams.update({'font.size': fsize})
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]

xvec = [0.0, -0.05, 1]

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(projection='3d')
bbox = ax.get_position()
new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
ax.set_position(new_bbox)
for i in range(len(width)):
    ax.plot_surface(Theta / np.pi, Phi / np.pi, probabilities[:, :, i], color=colors[i], alpha=transparency[i],
                    rstride=10, cstride=10, label=fr'$\epsilon$ = {round(width[i]/(2*np.pi) * 1000, 1)}GHz')
ax.set_xlabel(r'Rot. angle, $\theta$/$\pi$ rad', labelpad=15)
ax.set_ylabel('Deph. angle, \n$\phi$/$\pi$ rad', labelpad=25)
# ax.set_zlabel(r'$\left|\mathrm{1}\right\rangle$ meas. prob.', labelpad=10)
ax.set_zlabel(r"$|\left\langle\Psi'^-_t|\Psi^-_f\right\rangle|^2$/10$^{-4}$", labelpad=25)
# ax.zaxis.label.set_rotation(-90)
# ax.zaxis.set_rotate_label(False)
#ax.set_zlim([0, 1])
ax.set_xticks([0, 1, 2])
ax.set_xticklabels(['0', '1', '2'])
ax.set_yticks([0, 1, 2])
ax.set_yticklabels(['0', '1', '2'])
ax.set_zticks([0, 0.000025, 0.00005, 0.000075, 0.0001, 0.000125])
ax.set_zticklabels(['0.00', '0.25', '0.50', '0.75', '1.00', '1.25'])
ax.tick_params(axis='z', which='major', pad=10)
# plt.legend(bbox_to_anchor=(-0.5, 0.5), loc='center left', fontsize=18)
plt.legend(bbox_to_anchor=(1, 1.11), loc='upper right', fontsize=18)
ax.view_init(10, -60)
plt.show()



fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(projection='3d')
bbox = ax.get_position()
new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
ax.set_position(new_bbox)
ax.plot_surface(Theta / np.pi, Phi / np.pi, probabilities[:, :, 0], color=colors[0], alpha=transparency[0],
                    rstride=10, cstride=10, label=fr'$\epsilon$ = {round(width[0]/(2*np.pi) * 1000, 1)}GHz')
ax.set_xlabel(r'Rot. angle, $\theta$/$\pi$ rad', labelpad=15)
ax.set_ylabel('Deph. angle, \n$\phi$/$\pi$ rad', labelpad=25)
# ax.set_zlabel(r'$\left|\mathrm{1}\right\rangle$ meas. prob.', labelpad=10)
ax.set_zlabel(r"$|\left\langle\Psi'^-_t|\Psi^-_f\right\rangle|^2$", labelpad=10)
ax.set_xticks([0, 0.5, 1])
ax.set_xticklabels(['0.0', '0.5', '1.0'])
ax.set_yticks([0, 0.5, 1])
ax.set_yticklabels(['0.0', '0.5', '1.0'])
# ax.set_zticks([0, 0.5, 1])
# ax.set_zticklabels(['0.0', '0.5', '1.0'])
plt.legend(bbox_to_anchor=(1, 1.11), loc='upper right', fontsize=18)
ax.view_init(10, -60)
plt.show()



xvec = [0.10, 0.10, 0.8]

fsize = 30  # fontsize for the figures
mpl.rcParams.update({'font.size': fsize})

axlabels = [r"$|\left\langle\Psi'^-_t|\Psi^-_f\right\rangle|^2$/10$^{-6}$", r"$|\left\langle\Psi'^-_t|\Psi^-_f\right\rangle|^2$/10$^{-4}$", r"$|\left\langle\Psi'^-_t|\Psi^-_f\right\rangle|^2$/10$^{-5}$"]
axticks = [[8.2e-6, 8.25e-6, 8.3e-6], [1.1e-4, 1.15e-4, 1.2e-4, 1.25e-4], [5.15e-5, 5.25e-5, 5.35e-5]]
axticklabels = [['8.20', '8.25', '8.30'], ['1.10', '1.15', '1.20', '1.25'], ['5.15', '5.25', '5.35']]

for i in range(len(width)):
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    bbox = ax.get_position()
    new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    im1 = ax.imshow(probabilities[:, :, i], cmap='inferno',
                    extent=[theta[0] / np.pi, theta[-1] / np.pi, phi[0] / np.pi, phi[-1] / np.pi],
                    interpolation='nearest', origin='lower',
                    aspect=(theta[0] - theta[-1]) / (phi[0] - phi[-1]))
    ax.set_xlabel(r'Rot. angle $\theta$ ($\pi$ rad)')
    ax.set_ylabel(r'Deph. angle $\phi$ ($\pi$ rad)')
    ax.set_title(fr'$\epsilon$ = {round(width[i]/(2*np.pi)* 1000, 1)} GHz', pad=10)
    # ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    # ax.set_xticklabels(['0.00', '0.25', '0.50', '0.75', '1.00'])
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(['0', '1', '2'])
    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(['0', '1', '2'])
    ax.tick_params(axis='x', pad=9)
    axins = inset_axes(
        ax,
        width="5%",  # width: 5% of parent_bbox width
        height="100%",  # height: 50%
        loc="lower left",
        bbox_to_anchor=(1.05, 0., 1, 1),
        bbox_transform=ax.transAxes,
        borderpad=0,
    )
    cbar = fig.colorbar(im1, cax=axins, ticks=axticks[i])
    cbar.ax.set_yticklabels(axticklabels[i])
    cbar.set_label(axlabels[i], labelpad=10)
    plt.show()