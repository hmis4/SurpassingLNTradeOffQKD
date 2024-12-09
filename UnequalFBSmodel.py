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
    normfactor = 2*np.pi*sigma**2*(1 - np.exp(-((tau2-tau1)**2)/(2*sigma**2)))
    return 1/np.sqrt(normfactor)

def overlapshift(tau1, tau2, sigma_t, shift, start, end):
    coeff = (np.sqrt(np.pi)/2)*sigma_t
    phase = np.exp(-1j*shift*(tau1+tau2)/2)
    shiftdecay = np.exp(-sigma_t**2*shift**2/4)
    bindecay = np.exp(-(tau1-tau2)**2/(4*sigma_t**2))
    endval = erf(sigma_t * end + sigma_t * shift/2 + 1j*sigma_t * (tau2-tau1)/(2*sigma_t**2))
    startval = erf(sigma_t * start + sigma_t * shift / 2 + 1j * sigma_t * (tau2 - tau1) / (2 * sigma_t ** 2))
    return coeff*phase*shiftdecay*bindecay*(endval - startval)

def Uoverlap1D(tau1,tau2, sigma_t, Omega, epsilon, mu, theta=0, phi=0, delta=0):
    f1 = overlapshift(tau1, tau2, sigma_t, 0, -math.inf, Omega-epsilon/2)
    f2 = np.cos(theta)*overlapshift(tau1, tau2, sigma_t, 0, Omega-epsilon/2, Omega+epsilon/2) - np.exp(1j*delta)*np.exp(-1j*phi)*np.sin(theta)*overlapshift(tau1, tau2, sigma_t, mu, Omega-epsilon/2, Omega+epsilon/2)
    f3 = overlapshift(tau1, tau2, sigma_t, 0, Omega + epsilon/2, Omega + mu - epsilon/2)
    f4 = np.exp(1j*delta)*np.cos(theta)*overlapshift(tau1, tau2, sigma_t, 0, Omega+mu-epsilon/2, Omega+mu+epsilon/2) + np.exp(1j*phi)*np.sin(theta)*overlapshift(tau1, tau2, sigma_t, -mu, Omega+mu-epsilon/2, Omega+mu+epsilon/2)
    f5 = overlapshift(tau1, tau2, sigma_t, 0, Omega+mu+epsilon/2, math.inf)
    return f1 + f2 + f3 + f4 + f5


def Uoverlap2Deps(tau1,tau2, sigma_t, Omega, epsilon, mu, theta=0, phi=0, delta=0, diff=0):
    return 2*A(tau1, tau2, sigma_t)**2*(Uoverlap1D(tau1, tau1, sigma_t, Omega, epsilon, mu, theta, phi, delta)*Uoverlap1D(tau2, tau2, sigma_t, Omega, epsilon+diff, mu, theta, phi, delta)-Uoverlap1D(tau1, tau2, sigma_t, Omega, epsilon, mu, theta, phi, delta)*Uoverlap1D(tau2, tau1, sigma_t, Omega, epsilon+diff, mu, theta, phi, delta))

def Uprob2Deps(tau1, tau2, sigma_t, Omega, epsilon, mu, theta=0, phi=0, delta=0, diff=0):
    return np.abs(Uoverlap2Deps(tau1, tau2, sigma_t, Omega, epsilon, mu, theta, phi, delta, diff))**2


def Uoverlap2Dtheta(tau1,tau2, sigma_t, Omega, epsilon, mu, theta=0, phi=0, delta=0, diff=0):
    return 2*A(tau1, tau2, sigma_t)**2*(Uoverlap1D(tau1, tau1, sigma_t, Omega, epsilon, mu, theta, phi, delta)*Uoverlap1D(tau2, tau2, sigma_t, Omega, epsilon, mu, theta+diff, phi, delta)-Uoverlap1D(tau1, tau2, sigma_t, Omega, epsilon, mu, theta, phi, delta)*Uoverlap1D(tau2, tau1, sigma_t, Omega, epsilon, mu, theta+diff, phi, delta))

def Uprob2Dtheta(tau1, tau2, sigma_t, Omega, epsilon, mu, theta=0, phi=0, delta=0, diff=0):
    return np.abs(Uoverlap2Dtheta(tau1, tau2, sigma_t, Omega, epsilon, mu, theta, phi, delta, diff))**2


def Uoverlap2Dphi(tau1,tau2, sigma_t, Omega, epsilon, mu, theta=0, phi=0, delta=0, diff=0):
    return 2*A(tau1, tau2, sigma_t)**2*(Uoverlap1D(tau1, tau1, sigma_t, Omega, epsilon, mu, theta, phi, delta)*Uoverlap1D(tau2, tau2, sigma_t, Omega, epsilon, mu, theta, phi+diff, delta)-Uoverlap1D(tau1, tau2, sigma_t, Omega, epsilon, mu, theta, phi, delta)*Uoverlap1D(tau2, tau1, sigma_t, Omega, epsilon, mu, theta, phi+diff, delta))

def Uprob2Dphi(tau1, tau2, sigma_t, Omega, epsilon, mu, theta=0, phi=0, delta=0, diff=0):
    return np.abs(Uoverlap2Dphi(tau1, tau2, sigma_t, Omega, epsilon, mu, theta, phi, delta, diff))**2


def Uoverlap2Domega(tau1,tau2, sigma_t, Omega, epsilon, mu, theta=0, phi=0, delta=0, diff=0):
    return 2*A(tau1, tau2, sigma_t)**2*(Uoverlap1D(tau1, tau1, sigma_t, Omega, epsilon, mu, theta, phi, delta)*Uoverlap1D(tau2, tau2, sigma_t, Omega+diff, epsilon, mu, theta, phi, delta)-Uoverlap1D(tau1, tau2, sigma_t, Omega, epsilon, mu, theta, phi, delta)*Uoverlap1D(tau2, tau1, sigma_t, Omega+diff, epsilon, mu, theta, phi, delta))

def Uprob2Domega(tau1, tau2, sigma_t, Omega, epsilon, mu, theta=0, phi=0, delta=0, diff=0):
    return np.abs(Uoverlap2Domega(tau1, tau2, sigma_t, Omega, epsilon, mu, theta, phi, delta, diff))**2


def Uoverlap2Dmu(tau1,tau2, sigma_t, Omega, epsilon, mu, theta=0, phi=0, delta=0, diff=0):
    return 2*A(tau1, tau2, sigma_t)**2*(Uoverlap1D(tau1, tau1, sigma_t, Omega, epsilon, mu, theta, phi, delta)*Uoverlap1D(tau2, tau2, sigma_t, Omega, epsilon, mu+diff, theta, phi, delta)-Uoverlap1D(tau1, tau2, sigma_t, Omega, epsilon, mu, theta, phi, delta)*Uoverlap1D(tau2, tau1, sigma_t, Omega, epsilon, mu+diff, theta, phi, delta))

def Uprob2Dmu(tau1, tau2, sigma_t, Omega, epsilon, mu, theta=0, phi=0, delta=0, diff=0):
    return np.abs(Uoverlap2Dmu(tau1, tau2, sigma_t, Omega, epsilon, mu, theta, phi, delta, diff))**2


def conjoverlapshift(Omega, tau, sigma_w, sigma_t, shift, start, end):
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




def conjUoverlap1D(Omega1,tau2, sigma_w,sigma_t, Omega, epsilon, mu, theta=0, phi=0):
    print(Omega1, tau2)
    f1 = conjoverlapshift(Omega1, tau2, sigma_w, sigma_t, 0, -math.inf, Omega-epsilon/2)
    f2 = np.cos(theta)*conjoverlapshift(Omega1, tau2, sigma_w, sigma_t, 0, Omega-epsilon/2, Omega+epsilon/2) + np.exp(1j*phi)*np.sin(theta)*conjoverlapshift(Omega1, tau2, sigma_w, sigma_t, mu, Omega-epsilon/2, Omega+epsilon/2)
    f3 = conjoverlapshift(Omega1, tau2, sigma_w, sigma_t, 0, Omega + epsilon/2, Omega + mu - epsilon/2)
    f4 = np.cos(theta)*conjoverlapshift(Omega1, tau2, sigma_w, sigma_t, 0, Omega+mu-epsilon/2, Omega+mu+epsilon/2) - np.exp(-1j*phi)*np.sin(theta)*conjoverlapshift(Omega1, tau2, sigma_w, sigma_t, -mu, Omega+mu-epsilon/2, Omega+mu+epsilon/2)
    f5 = conjoverlapshift(Omega1, tau2, sigma_w, sigma_t, 0, Omega+mu+epsilon/2, math.inf)
    return f1 + f2 + f3 + f4 + f5

def conjUoverlap2D(tau1,tau2, sigma_t, Omega1, Omega2, sigma_w, Omega, epsilon, mu, theta=0, phi=0,diff=0):
    return 2*A(tau1, tau2, sigma_t)*A(Omega1, Omega2, sigma_w)*(conjUoverlap1D(Omega1, tau1, sigma_w, sigma_t, Omega, epsilon, mu, theta, phi)*conjUoverlap1D(Omega2, tau2, sigma_w, sigma_t, Omega, epsilon+diff, mu, theta, phi)-conjUoverlap1D(Omega1, tau2, sigma_w, sigma_t, Omega, epsilon, mu, theta, phi)*conjUoverlap1D(Omega2, tau1, sigma_w, sigma_t, Omega, epsilon+diff, mu, theta, phi))


def conjUprob2D(tau1,tau2, sigma_t, Omega1, Omega2, sigma_w, Omega, epsilon, mu, theta=0, phi=0,diff=0):
    return np.abs(conjUoverlap2D(tau1,tau2, sigma_t, Omega1, Omega2, sigma_w, Omega, epsilon, mu, theta, phi, diff))**2

# def conjUoverlap2D(tau1,tau2, sigma_t, Omega1, Omega2, sigma_w, Omega, epsilon, mu, theta=0, phi=0,diff=0):
#     return 2*A(tau1, tau2, sigma_t)*A(Omega1, Omega2, sigma_w)*(conjUoverlap1D(Omega1, tau1, sigma_w, sigma_t, Omega, epsilon, mu, theta, phi)*conjUoverlap1D(Omega1, tau2, sigma_w, sigma_t, Omega, epsilon+diff, mu, theta, phi)+ conjUoverlap1D(Omega2, tau1, sigma_w, sigma_t, Omega, epsilon, mu, theta, phi)*conjUoverlap1D(Omega2, tau2, sigma_w, sigma_t, Omega, epsilon+diff, mu, theta, phi)-conjUoverlap1D(Omega1, tau2, sigma_w, sigma_t, Omega, epsilon, mu, theta, phi)*conjUoverlap1D(Omega1, tau1, sigma_w, sigma_t, Omega, epsilon+diff, mu, theta, phi) - conjUoverlap1D(Omega2, tau2, sigma_w, sigma_t, Omega, epsilon, mu, theta, phi)*conjUoverlap1D(Omega2, tau1, sigma_w, sigma_t, Omega, epsilon+diff, mu, theta, phi))
#
#
# def conjUprob2D(tau1,tau2, sigma_t, Omega1, Omega2, sigma_w, Omega, epsilon, mu, theta=0, phi=0,diff=0):
#     return np.abs(conjUoverlap2D(tau1,tau2, sigma_t, Omega1, Omega2, sigma_w, Omega, epsilon, mu, theta, phi, diff))**2


#krcompphase = True
difftheta = True
diffepsilon = True
diffphi = True
diffmu = True
diffOmega = True


if difftheta:
    tau1 = 0
    tau2 = 220 # ps
    sigma_t = 17 # ps
    gap = 0.019*2*np.pi # THz #0.019*2*np.pi
    Omega1 = 0
    sigma_w = 0.0005 # THz
    width = 2*np.pi*6*sigma_w

    theta=np.linspace(0,np.pi, 101)
    phi=np.linspace(0,np.pi,101)
    delta=0

    difference = 2*np.pi*np.asarray([0, 0.001, 0.01, 0.1])

    #print(Uprob2D(tau1,tau2, sigma_t, Omega1, width, gap,0,0,0))

    Theta, Phi = np.meshgrid(theta, phi)
    probabilities = np.zeros((len(theta), len(phi), len(difference)))
    biterror = np.zeros((len(theta), len(phi), len(difference)))

    for i in range(len(difference)):
        probabilities[:,:,i] = Uprob2Dtheta(tau1, tau2, sigma_t, Omega1, width, gap, Theta, Phi, delta, difference[i])
        biterror[:,:,i] = conjUprob2D(tau1, tau2, sigma_t, Omega1, width, gap, Theta, Phi, delta, difference[i])

    cmap = mpl.colormaps['inferno']
    colors = cmap([0.9, 0.7, 0.5, 0.3])
    transparency = [0.55, 0.7, 0.85, 1]
    # colors = cmap([0.4, 0.7, 0.9])
    # transparency = [1, 0.85, 0.7]

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
    for i in range(len(difference)):
        ax.plot_surface(Theta / np.pi, Phi / np.pi, probabilities[:, :, i], color=colors[i], alpha=transparency[i],
                        rstride=10, cstride=10, label=fr'$\delta_\theta$ = {difference[i]/(np.pi)}$\pi$ rad')
    ax.set_xlabel(r"Rot. angle, "
           "\n" 
           r"$\theta$ ($\pi$ rad)", labelpad=25)
    ax.set_ylabel('Deph. angle, \n$\phi$ ($\pi$ rad)', labelpad=25)
    ax.set_zlabel(r"$|\left\langle\Psi'^-_t|\Psi^-_t\right\rangle|^2$", labelpad=25)
    ax.set_zlim([0, 1])
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(['0.0', '0.5', '1.0'])
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(['0.0', '0.5', '1.0'])
    ax.set_zticks([0, 0.5, 1])
    ax.set_zticklabels(['0.0', '0.5', '1.0'])
    ax.tick_params(axis='z', pad=10)
    ax.tick_params(axis='x', pad=-1)
    ax.tick_params(axis='y', pad=5)
    # plt.legend(bbox_to_anchor=(-0.5, 0.5), loc='center left', fontsize=18)
    plt.legend(bbox_to_anchor=(1, 1.16), loc='upper right', fontsize=18)
    #ax.view_init(10, -60)
    ax.view_init(10, -83)
    plt.show()


    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(projection='3d')
    bbox = ax.get_position()
    new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    for i in range(len(difference)):
        ax.plot_surface(Theta / np.pi, Phi / np.pi, biterror[:, :, i], color=colors[i], alpha=transparency[i],
                        rstride=10, cstride=10, label=fr'$\delta_\theta$ = {difference[i]/(np.pi)}$\pi$ rad')
    ax.set_xlabel(r"Rot. angle, "
           "\n" 
           r"$\theta$ ($\pi$ rad)", labelpad=25)
    ax.set_ylabel('Deph. angle, \n$\phi$ ($\pi$ rad)', labelpad=25)
    ax.set_zlabel(r"$|\left\langle\Psi'^-_t|\Psi^-_t\right\rangle|^2$", labelpad=25)
    # ax.set_zlim([0, 1])
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(['0.0', '0.5', '1.0'])
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(['0.0', '0.5', '1.0'])
    # ax.set_zticks([0, 0.5, 1])
    # ax.set_zticklabels(['0.0', '0.5', '1.0'])
    ax.tick_params(axis='z', pad=10)
    ax.tick_params(axis='x', pad=-1)
    ax.tick_params(axis='y', pad=5)
    # plt.legend(bbox_to_anchor=(-0.5, 0.5), loc='center left', fontsize=18)
    plt.legend(bbox_to_anchor=(1, 1.16), loc='upper right', fontsize=18)
    #ax.view_init(10, -60)
    ax.view_init(10, -60)
    plt.show()




    overlaps = np.zeros((len(theta), len(phi), len(difference)), dtype=complex)
    mags = np.zeros((len(theta), len(phi), len(difference)))
    phases = np.zeros((len(theta), len(phi), len(difference)))

    for i in range(len(difference)):
        overlaps[:,:,i] = Uoverlap2Dtheta(tau1, tau2, sigma_t, Omega1, width, gap, Theta, Phi, delta, difference[i])


    #mags = np.abs(overlaps)**2
    phases = np.arctan(np.divide(np.imag(overlaps),np.real(overlaps)))

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(projection='3d')
    bbox = ax.get_position()
    new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    for i in range(len(difference)):
        ax.plot_surface(Theta / np.pi, Phi / np.pi, phases[:, :, i]/np.pi, color=colors[i], alpha=transparency[i],
                        rstride=10, cstride=10, label=fr'$\delta_\theta$ = {difference[i]/(np.pi)}$\pi$ rad')
    ax.set_xlabel(r"Rot. angle, "
           "\n" 
           r"$\theta$ ($\pi$ rad)", labelpad=25)
    ax.set_ylabel('Deph. angle, \n$\phi$ ($\pi$ rad)', labelpad=25)
    ax.set_zlabel(r"Global phase/$\pi$ rad", labelpad=25)
    # ax.set_zlim([0, 1])
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(['0.0', '0.5', '1.0'])
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(['0.0', '0.5', '1.0'])
    # ax.set_zticks([0, 0.5, 1])
    # ax.set_zticklabels(['0.0', '0.5', '1.0'])
    ax.tick_params(axis='z', pad=10)
    ax.tick_params(axis='x', pad=-1)
    ax.tick_params(axis='y', pad=5)
    # plt.legend(bbox_to_anchor=(-0.5, 0.5), loc='center left', fontsize=18)
    plt.legend(bbox_to_anchor=(1, 1.16), loc='upper right', fontsize=18)
    #ax.view_init(10, -60)
    ax.view_init(10, -60)
    plt.show()



if diffepsilon:
    tau1 = 0
    tau2 = 220 # ps
    sigma_t = 17 # ps
    gap = 0.019*2*np.pi # THz #0.019*2*np.pi
    Omega1 = 0
    sigma_w = 0.0005 # THz
    width = 2*np.pi*6*sigma_w

    theta=np.linspace(0,np.pi, 101)
    phi=np.linspace(0,np.pi,101)
    delta=0

    difference = width*np.asarray([0, 0.01, 0.1, 1])

    #print(Uprob2D(tau1,tau2, sigma_t, Omega1, width, gap,0,0,0))

    Theta, Phi = np.meshgrid(theta, phi)
    probabilities = np.zeros((len(theta), len(phi), len(difference)))
    biterror = np.zeros((len(theta), len(phi), len(difference)))

    for i in range(len(difference)):
        probabilities[:,:,i] = Uprob2Deps(tau1, tau2, sigma_t, Omega1, width, gap, Theta, Phi, delta, difference[i])
        biterror[:, :, i] = conjUprob2D(tau1, tau2, sigma_t, Omega1, width, gap, Theta, Phi, delta, difference[i])

    cmap = mpl.colormaps['inferno']
    colors = cmap([0.9, 0.7, 0.5, 0.3])
    transparency = [0.55, 0.7, 0.85, 1]
    # colors = cmap([0.4, 0.7, 0.9])
    # transparency = [1, 0.85, 0.7]

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
    for i in range(len(difference)):
        ax.plot_surface(Theta / np.pi, Phi / np.pi, probabilities[:, :, i], color=colors[i], alpha=transparency[i],
                        rstride=10, cstride=10, label=fr'$\delta_\epsilon$ = {difference[i]/width}$\epsilon$')
    ax.set_xlabel(r"Rot. angle, "
           "\n" 
           r"$\theta$ ($\pi$ rad)", labelpad=25)
    ax.set_ylabel('Deph. angle, \n$\phi$ ($\pi$ rad)', labelpad=25)
    ax.set_zlabel(r"$|\left\langle\Psi'^-_t|\Psi^-_t\right\rangle|^2$", labelpad=10)
    ax.set_zlim([0, 1])
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(['0.0', '0.5', '1.0'])
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(['0.0', '0.5', '1.0'])
    ax.set_zticks([0, 0.5, 1])
    ax.set_zticklabels(['0.0', '0.5', '1.0'])
    # plt.legend(bbox_to_anchor=(-0.5, 0.5), loc='center left', fontsize=18)
    plt.legend(bbox_to_anchor=(1, 1.16), loc='upper right', fontsize=18)
    ax.view_init(10, -60)
    plt.show()


    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(projection='3d')
    bbox = ax.get_position()
    new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    for i in range(len(difference)):
        ax.plot_surface(Theta / np.pi, Phi / np.pi, biterror[:, :, i], color=colors[i], alpha=transparency[i],
                        rstride=10, cstride=10, label=fr'$\delta_\epsilon$ = {difference[i]/width}$\epsilon$')
    ax.set_xlabel(r"Rot. angle, "
           "\n" 
           r"$\theta$ ($\pi$ rad)", labelpad=25)
    ax.set_ylabel('Deph. angle, \n$\phi$ ($\pi$ rad)', labelpad=25)
    ax.set_zlabel(r"$|\left\langle\Psi'^-_t|\Psi^-_t\right\rangle|^2$", labelpad=10)
    # ax.set_zlim([0, 1])
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(['0.0', '0.5', '1.0'])
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(['0.0', '0.5', '1.0'])
    # ax.set_zticks([0, 0.5, 1])
    # ax.set_zticklabels(['0.0', '0.5', '1.0'])
    # plt.legend(bbox_to_anchor=(-0.5, 0.5), loc='center left', fontsize=18)
    plt.legend(bbox_to_anchor=(1, 1.16), loc='upper right', fontsize=18)
    ax.view_init(10, -60)
    plt.show()


    overlaps = np.zeros((len(theta), len(phi), len(difference)), dtype=complex)
    mags = np.zeros((len(theta), len(phi), len(difference)))
    phases = np.zeros((len(theta), len(phi), len(difference)))

    for i in range(len(difference)):
        overlaps[:,:,i] = Uoverlap2Deps(tau1, tau2, sigma_t, Omega1, width, gap, Theta, Phi, delta, difference[i])


    #mags = np.abs(overlaps)**2
    phases = np.arctan(np.divide(np.imag(overlaps),np.real(overlaps)))

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(projection='3d')
    bbox = ax.get_position()
    new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    for i in range(len(difference)):
        ax.plot_surface(Theta / np.pi, Phi / np.pi, phases[:, :, i]/np.pi, color=colors[i], alpha=transparency[i],
                        rstride=10, cstride=10, label=fr'$\delta_\epsilon$ = {difference[i]/width}$\epsilon$')
    ax.set_xlabel(r"Rot. angle, "
           "\n" 
           r"$\theta$ ($\pi$ rad)", labelpad=25)
    ax.set_ylabel('Deph. angle, \n$\phi$ ($\pi$ rad)', labelpad=25)
    ax.set_zlabel(r"Global phase/$\pi$ rad", labelpad=25)
    # ax.set_zlim([0, 1])
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(['0.0', '0.5', '1.0'])
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(['0.0', '0.5', '1.0'])
    # ax.set_zticks([0, 0.5, 1])
    # ax.set_zticklabels(['0.0', '0.5', '1.0'])
    ax.tick_params(axis='z', pad=10)
    ax.tick_params(axis='x', pad=-1)
    ax.tick_params(axis='y', pad=5)
    # plt.legend(bbox_to_anchor=(-0.5, 0.5), loc='center left', fontsize=18)
    plt.legend(bbox_to_anchor=(1, 1.16), loc='upper right', fontsize=18)
    #ax.view_init(10, -60)
    ax.view_init(10, -60)
    plt.show()




if diffphi:
    tau1 = 0
    tau2 = 220 # ps
    sigma_t = 17 # ps
    gap = 0.019*2*np.pi # THz #0.019*2*np.pi
    Omega1 = 0
    sigma_w = 0.0005 # THz
    width = 2*np.pi*6*sigma_w

    theta=np.linspace(0,np.pi, 101)
    phi=np.linspace(0,np.pi,101)
    delta=0

    difference = 2*np.pi*np.asarray([0, 0.001, 0.01, 0.1])

    #print(Uprob2D(tau1,tau2, sigma_t, Omega1, width, gap,0,0,0))

    Theta, Phi = np.meshgrid(theta, phi)
    probabilities = np.zeros((len(theta), len(phi), len(difference)))

    for i in range(len(difference)):
        probabilities[:,:,i] = Uprob2Dphi(tau1, tau2, sigma_t, Omega1, width, gap, Theta, Phi, delta, difference[i])

    cmap = mpl.colormaps['inferno']
    colors = cmap([0.9, 0.7, 0.5, 0.3])
    transparency = [0.55, 0.7, 0.85, 1]
    # colors = cmap([0.4, 0.7, 0.9])
    # transparency = [1, 0.85, 0.7]

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
    for i in range(len(difference)):
        ax.plot_surface(Theta / np.pi, Phi / np.pi, probabilities[:, :, i], color=colors[i], alpha=transparency[i],
                        rstride=10, cstride=10, label=fr'$\delta_\phi$ = {difference[i]/(np.pi)}$\pi$ rad')
    ax.set_xlabel(r"Rot. angle, "
           "\n" 
           r"$\theta$ ($\pi$ rad)", labelpad=25)
    ax.set_ylabel('Deph. angle, \n$\phi$ ($\pi$ rad)', labelpad=25)
    ax.set_zlabel(r"$|\left\langle\Psi'^-_t|\Psi^-_t\right\rangle|^2$", labelpad=10)
    ax.set_zlim([0, 1])
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(['0.0', '0.5', '1.0'])
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(['0.0', '0.5', '1.0'])
    ax.set_zticks([0, 0.5, 1])
    ax.set_zticklabels(['0.0', '0.5', '1.0'])
    # plt.legend(bbox_to_anchor=(-0.5, 0.5), loc='center left', fontsize=18)
    plt.legend(bbox_to_anchor=(1, 1.16), loc='upper right', fontsize=18)
    ax.view_init(10, -60)
    plt.show()

    overlaps = np.zeros((len(theta), len(phi), len(difference)), dtype=complex)
    mags = np.zeros((len(theta), len(phi), len(difference)))
    phases = np.zeros((len(theta), len(phi), len(difference)))

    for i in range(len(difference)):
        overlaps[:,:,i] = Uoverlap2Deps(tau1, tau2, sigma_t, Omega1, width, gap, Theta, Phi, delta, difference[i])


    #mags = np.abs(overlaps)**2
    phases = np.arctan(np.divide(np.imag(overlaps),np.real(overlaps)))

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(projection='3d')
    bbox = ax.get_position()
    new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    for i in range(len(difference)):
        ax.plot_surface(Theta / np.pi, Phi / np.pi, phases[:, :, i]/np.pi, color=colors[i], alpha=transparency[i],
                        rstride=10, cstride=10, label=fr'$\delta_\phi$ = {difference[i]/(np.pi)}$\pi$ rad')
    ax.set_xlabel(r"Rot. angle, "
           "\n" 
           r"$\theta$ ($\pi$ rad)", labelpad=25)
    ax.set_ylabel('Deph. angle, \n$\phi$ ($\pi$ rad)', labelpad=25)
    ax.set_zlabel(r"Global phase/$\pi$ rad", labelpad=25)
    # ax.set_zlim([0, 1])
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(['0.0', '0.5', '1.0'])
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(['0.0', '0.5', '1.0'])
    # ax.set_zticks([0, 0.5, 1])
    # ax.set_zticklabels(['0.0', '0.5', '1.0'])
    ax.tick_params(axis='z', pad=10)
    ax.tick_params(axis='x', pad=-1)
    ax.tick_params(axis='y', pad=5)
    # plt.legend(bbox_to_anchor=(-0.5, 0.5), loc='center left', fontsize=18)
    plt.legend(bbox_to_anchor=(1, 1.16), loc='upper right', fontsize=18)
    #ax.view_init(10, -60)
    ax.view_init(10, -60)
    plt.show()



if diffmu:
    tau1 = 0
    tau2 = 220 # ps
    sigma_t = 17 # ps
    gap = 0.019*2*np.pi # THz #0.019*2*np.pi
    Omega1 = 0
    sigma_w = 0.0005 # THz
    width = 2*np.pi*6*sigma_w

    theta=np.linspace(0,np.pi, 101)
    phi=np.linspace(0,np.pi,101)
    delta=0

    difference = gap*np.asarray([0, 0.01, 0.1, 1])

    #print(Uprob2D(tau1,tau2, sigma_t, Omega1, width, gap,0,0,0))

    Theta, Phi = np.meshgrid(theta, phi)
    probabilities = np.zeros((len(theta), len(phi), len(difference)))

    for i in range(len(difference)):
        probabilities[:,:,i] = Uprob2Dmu(tau1, tau2, sigma_t, Omega1, width, gap, Theta, Phi, delta, difference[i])

    cmap = mpl.colormaps['inferno']
    colors = cmap([0.3, 0.5, 0.7, 0.9])
    transparency = [1, 0.85, 0.7, 0.55]
    # colors = cmap([0.4, 0.7, 0.9])
    # transparency = [1, 0.85, 0.7]

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
    for i in range(len(difference)):
        ax.plot_surface(Theta / np.pi, Phi / np.pi, probabilities[:, :, i], color=colors[i], alpha=transparency[i],
                        rstride=10, cstride=10, label=fr'$\delta_\mu$ = {difference[i]/gap}$\mu$')
    ax.set_xlabel(r"Rot. angle, "
           "\n" 
           r"$\theta$ ($\pi$ rad)", labelpad=25)
    ax.set_ylabel('Deph. angle, \n$\phi$ ($\pi$ rad)', labelpad=25)
    ax.set_zlabel(r"$|\left\langle\Psi'^-_t|\Psi^-_t\right\rangle|^2$", labelpad=10)
    ax.set_zlim([0, 1])
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(['0.0', '0.5', '1.0'])
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(['0.0', '0.5', '1.0'])
    ax.set_zticks([0, 0.5, 1])
    ax.set_zticklabels(['0.0', '0.5', '1.0'])
    # plt.legend(bbox_to_anchor=(-0.5, 0.5), loc='center left', fontsize=18)
    plt.legend(bbox_to_anchor=(1, 1.16), loc='upper right', fontsize=18)
    ax.view_init(10, -60)
    plt.show()


    overlaps = np.zeros((len(theta), len(phi), len(difference)), dtype=complex)
    mags = np.zeros((len(theta), len(phi), len(difference)))
    phases = np.zeros((len(theta), len(phi), len(difference)))

    for i in range(len(difference)):
        overlaps[:,:,i] = Uoverlap2Dmu(tau1, tau2, sigma_t, Omega1, width, gap, Theta, Phi, delta, difference[i])


    #mags = np.abs(overlaps)**2
    phases = np.arctan(np.divide(np.imag(overlaps),np.real(overlaps)))

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(projection='3d')
    bbox = ax.get_position()
    new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    for i in range(len(difference)):
        ax.plot_surface(Theta / np.pi, Phi / np.pi, phases[:, :, i]/np.pi, color=colors[i], alpha=transparency[i],
                        rstride=10, cstride=10, label=fr'$\delta_\mu$ = {difference[i]/gap}$\mu$')
    ax.set_xlabel(r"Rot. angle, "
           "\n" 
           r"$\theta$ ($\pi$ rad)", labelpad=25)
    ax.set_ylabel('Deph. angle, \n$\phi$ ($\pi$ rad)', labelpad=25)
    ax.set_zlabel(r"Global phase/$\pi$ rad", labelpad=25)
    # ax.set_zlim([0, 1])
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(['0.0', '0.5', '1.0'])
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(['0.0', '0.5', '1.0'])
    # ax.set_zticks([0, 0.5, 1])
    # ax.set_zticklabels(['0.0', '0.5', '1.0'])
    ax.tick_params(axis='z', pad=10)
    ax.tick_params(axis='x', pad=-1)
    ax.tick_params(axis='y', pad=5)
    # plt.legend(bbox_to_anchor=(-0.5, 0.5), loc='center left', fontsize=18)
    plt.legend(bbox_to_anchor=(1, 1.16), loc='upper right', fontsize=18)
    #ax.view_init(10, -60)
    ax.view_init(10, -60)
    plt.show()


if diffOmega:
    tau1 = 0
    tau2 = 220 # ps
    sigma_t = 17 # ps
    gap = 0.019*2*np.pi # THz #0.019*2*np.pi
    Omega1 = 0
    sigma_w = 0.0005 # THz
    width = 2*np.pi*6*sigma_w

    theta=np.linspace(0,np.pi, 101)
    phi=np.linspace(0,np.pi,101)
    delta=0

    difference = 2*np.pi*0.01*np.asarray([0, 0.01, 0.1, 1])

    #print(Uprob2D(tau1,tau2, sigma_t, Omega1, width, gap,0,0,0))

    Theta, Phi = np.meshgrid(theta, phi)
    probabilities = np.zeros((len(theta), len(phi), len(difference)))

    for i in range(len(difference)):
        probabilities[:,:,i] = Uprob2Domega(tau1, tau2, sigma_t, Omega1, width, gap, Theta, Phi, delta, difference[i])

    cmap = mpl.colormaps['inferno']
    colors = cmap([0.3, 0.5, 0.7, 0.9])
    transparency = [1, 0.7, 0.85, 0.55]
    # colors = cmap([0.4, 0.7, 0.9])
    # transparency = [1, 0.85, 0.7]

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
    for i in range(len(difference)):
        ax.plot_surface(Theta / np.pi, Phi / np.pi, probabilities[:, :, i], color=colors[i], alpha=transparency[i],
                        rstride=10, cstride=10, label=fr'$\delta_\Omega$ = {round(difference[i]*1000/(2*np.pi),2)} GHz')
    ax.set_xlabel(r"Rot. angle, "
           "\n" 
           r"$\theta$ ($\pi$ rad)", labelpad=25)
    ax.set_ylabel('Deph. angle, \n$\phi$ ($\pi$ rad)', labelpad=25)
    ax.set_zlabel(r"$|\left\langle\Psi'^-_t|\Psi^-_t\right\rangle|^2$", labelpad=10)
    ax.set_zlim([0, 1])
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(['0.0', '0.5', '1.0'])
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(['0.0', '0.5', '1.0'])
    ax.set_zticks([0, 0.5, 1])
    ax.set_zticklabels(['0.0', '0.5', '1.0'])
    # plt.legend(bbox_to_anchor=(-0.5, 0.5), loc='center left', fontsize=18)
    plt.legend(bbox_to_anchor=(1, 1.16), loc='upper right', fontsize=18)
    ax.view_init(10, -60)
    plt.show()


    overlaps = np.zeros((len(theta), len(phi), len(difference)), dtype=complex)
    mags = np.zeros((len(theta), len(phi), len(difference)))
    phases = np.zeros((len(theta), len(phi), len(difference)))

    for i in range(len(difference)):
        overlaps[:,:,i] = Uoverlap2Domega(tau1, tau2, sigma_t, Omega1, width, gap, Theta, Phi, delta, difference[i])


    #mags = np.abs(overlaps)**2
    phases = np.arctan(np.divide(np.imag(overlaps),np.real(overlaps)))

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(projection='3d')
    bbox = ax.get_position()
    new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    for i in range(len(difference)):
        ax.plot_surface(Theta / np.pi, Phi / np.pi, phases[:, :, i]/np.pi, color=colors[i], alpha=transparency[i],
                        rstride=10, cstride=10, label=fr'$\delta_\Omega$ = {round(difference[i]*1000/(2*np.pi),2)} GHz')
    ax.set_xlabel(r"Rot. angle, "
           "\n" 
           r"$\theta$ ($\pi$ rad)", labelpad=25)
    ax.set_ylabel('Deph. angle, \n$\phi$ ($\pi$ rad)', labelpad=25)
    ax.set_zlabel(r"Global phase/$\pi$ rad", labelpad=25)
    # ax.set_zlim([0, 1])
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(['0.0', '0.5', '1.0'])
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(['0.0', '0.5', '1.0'])
    # ax.set_zticks([0, 0.5, 1])
    # ax.set_zticklabels(['0.0', '0.5', '1.0'])
    ax.tick_params(axis='z', pad=10)
    ax.tick_params(axis='x', pad=-1)
    ax.tick_params(axis='y', pad=5)
    # plt.legend(bbox_to_anchor=(-0.5, 0.5), loc='center left', fontsize=18)
    plt.legend(bbox_to_anchor=(1, 1.16), loc='upper right', fontsize=18)
    #ax.view_init(10, -60)
    ax.view_init(10, -60)
    plt.show()


