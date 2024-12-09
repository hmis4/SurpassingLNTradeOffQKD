import numpy as np
import math
import cmath
import matplotlib.pyplot as plt
from scipy.special import erf
import matplotlib as mpl
import pandas as pd
import scipy.interpolate as interpolate
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.font_manager import FontProperties
from matplotlib import rc_context

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

def Uoverlap2D(tau1,tau2, sigma_t, Omega, epsilon, mu, theta=0, phi=0, delta=0):
    return 2*A(tau1, tau2, sigma_t)**2*(Uoverlap1D(tau1, tau1, sigma_t, Omega, epsilon, mu, theta, phi, delta)*Uoverlap1D(tau2, tau2, sigma_t, Omega, epsilon, mu, theta, phi, delta)-Uoverlap1D(tau1, tau2, sigma_t, Omega, epsilon, mu, theta, phi, delta)*Uoverlap1D(tau2, tau1, sigma_t, Omega, epsilon, mu, theta, phi, delta))

def Uprob2D(tau1, tau2, sigma_t, Omega, epsilon, mu, theta=0, phi=0, delta=0):
    return np.abs(Uoverlap2D(tau1, tau2, sigma_t, Omega, epsilon, mu, theta, phi, delta))**2


findlim = False
comptau2 = True
comps_t = False
compOmega = False
compeps = False
compmu = False
compsigma = False
savedat= False
phasefac = True
krcomp= True


#krcompphase = True


if findlim:
    tau1 = 0
    tau2 = np.asarray([250, 300, 350]) # ps
    gap = 0.019*2*np.pi # THz #0.019*2*np.pi
    Omega1 = 0#-0.03*2*np.pi
    theta=np.pi
    phi=0
    delta=0

    sigma_t = np.linspace(0.01, 40, 1000)
    sigma_w = np.linspace(0.00001, 0.01,1000)

    ST, SW = np.meshgrid(sigma_t, sigma_w)
    probabilities = np.zeros((len(sigma_t), len(sigma_w), len(tau2)))

    for i in range(len(tau2)):
        probabilities[:,:,i] = Uprob2D(tau1, tau2[i], ST, Omega1, 12*np.pi*SW, gap, theta, phi, delta)


    cmap =  mpl.colormaps['inferno']
    colors = cmap([0.9, 0.7, 0.4])
    transparency = [0.7,0.85,1]

    fsize = 22  # fontsize for the figures
    mpl.rcParams.update({'font.size': fsize})
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial"]

    xvec = [0.0, -0.05, 1]

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(projection='3d')
    bbox = ax.get_position()
    new_bbox = (bbox.x0+xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    for i in range(len(tau2)):
        ax.plot_surface(ST, SW, probabilities[:,:,i], color=colors[i], alpha=transparency[i], rstride=10, cstride=10, label=fr'$\Omega_2$ = {tau2[i]} THz')
    ax.set_xlabel(r'sigma_t/ps', labelpad=15)
    ax.set_ylabel('sigma_w/THz', labelpad=25)
    ax.set_zlabel(r'$\left|\mathrm{1}\right\rangle$ meas. prob.', labelpad=10)
    ax.set_zlim([0, 1])
    ax.set_zticks([0,0.5,1])
    ax.set_zticklabels(['0.0', '0.5', '1.0'])
    plt.legend(bbox_to_anchor=(1, 1.11), loc='upper right', fontsize=18)
    ax.view_init(10, -60)
    plt.show()


    xvec = [0.10, 0.10, 0.8]

    fsize = 30  # fontsize for the figures
    mpl.rcParams.update({'font.size': fsize})

    for i in range(len(tau2)):
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111)
        bbox = ax.get_position()
        new_bbox = (bbox.x0+xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
        ax.set_position(new_bbox)
        im1 = ax.imshow(probabilities[:,:,i], cmap='inferno', vmin=0, vmax=1,
                       extent=[sigma_t[0], sigma_t[-1], sigma_w[0], sigma_w[-1]],
                       interpolation='nearest', origin='lower',
                       aspect=(sigma_t[0] - sigma_t[-1]) / (sigma_w[0] - sigma_w[-1]))
        ax.set_xlabel(r'sigma_t/ps', labelpad=15)
        ax.set_ylabel('sigma_w/THz', labelpad=25)
        ax.set_title(fr'$\Omega_2$ = {tau2[i]} THz', pad=10)
        axins = inset_axes(
            ax,
            width="5%",  # width: 5% of parent_bbox width
            height="100%",  # height: 50%
            loc="lower left",
            bbox_to_anchor=(1.05, 0., 1, 1),
            bbox_transform=ax.transAxes,
            borderpad=0,
        )
        cbar = fig.colorbar(im1, cax=axins)
        cbar.set_label(label=r'$\left|\mathrm{1}\right\rangle$ measurement prob.', labelpad=10)
        plt.show()


if comptau2:
    tau1 = 0
    tau2 = np.asarray([120, 220, 320]) # ps
    sigma_t = 17 # ps
    gap = 0.019*2*np.pi # THz #0.019*2*np.pi
    Omega1 = 0
    sigma_w = 0.0005 # THz
    width = 2*np.pi*6*sigma_w

    theta=np.linspace(0,np.pi, 101)
    phi=np.linspace(0,np.pi,101)
    delta=0

    Theta, Phi = np.meshgrid(theta, phi)
    probabilities = np.zeros((len(theta), len(phi), len(tau2)))

    for i in range(len(tau2)):
        probabilities[:,:,i] = Uprob2D(tau1, tau2[i], sigma_t, Omega1, width, gap, Theta, Phi, delta)

    cmap = mpl.colormaps['inferno']
    colors = cmap([0.4, 0.7, 0.9])
    transparency = [1, 0.85, 0.7]

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
    for i in range(len(tau2)):
        ax.plot_surface(Theta / np.pi, Phi / np.pi, probabilities[:, :, i], color=colors[i], alpha=transparency[i],
                        rstride=10, cstride=10, label=fr'$\tau_1$ = {tau2[i]} ps')
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
    plt.legend(bbox_to_anchor=(1, 1.11), loc='upper right', fontsize=18)
    ax.view_init(10, -60)
    plt.show()



    xvec = [0.10, 0.10, 0.8]

    fsize = 30  # fontsize for the figures
    mpl.rcParams.update({'font.size': fsize})

    for i in range(len(tau2)):
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111)
        bbox = ax.get_position()
        new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
        ax.set_position(new_bbox)
        im1 = ax.imshow(probabilities[:, :, i], cmap='inferno', vmin=0, vmax=1,
                        extent=[theta[0] / np.pi, theta[-1] / np.pi, phi[0] / np.pi, phi[-1] / np.pi],
                        interpolation='nearest', origin='lower',
                        aspect=(theta[0] - theta[-1]) / (phi[0] - phi[-1]))
        ax.set_xlabel(r'Rot. angle $\theta$/$\pi$ rad')
        ax.set_ylabel(r'Deph. angle $\phi$/$\pi$ rad')
        ax.set_title(fr'$\tau_2$ = {tau2[i]} ps', pad=10)
        # ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
        # ax.set_xticklabels(['0.00', '0.25', '0.50', '0.75', '1.00'])
        ax.set_xticks([0, 0.5, 1])
        ax.set_xticklabels(['0.00', '0.50', '1.00'])
        ax.set_yticks([0, 0.5, 1])
        ax.set_yticklabels(['0.00', '0.50', '1.00'])
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
        cbar = fig.colorbar(im1, cax=axins)
        cbar.set_label(label=r'$\left|\mathrm{1}\right\rangle$ measurement prob.', labelpad=10)
        plt.show()

if comps_t:
    tau1 = 0
    tau2 = 220 # ps
    #sigma_t = 0.01876
    sigma_t = np.asarray([10, 17, 30]) # ps
    gap = 0.019*2*np.pi # THz #0.019*2*np.pi
    Omega1 = 0
    sigma_w = 0.0005 # THz
    width = 2*np.pi*6*sigma_w

    theta=np.linspace(0,np.pi, 101)
    phi=np.linspace(0,np.pi,101)
    delta=0

    Theta, Phi = np.meshgrid(theta, phi)
    probabilities = np.zeros((len(theta), len(phi), len(sigma_t)))

    for i in range(len(sigma_t)):
        probabilities[:,:,i] = Uprob2D(tau1, tau2, sigma_t[i], Omega1, width, gap, Theta, Phi, delta)

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
    for i in range(len(sigma_t)):
        ax.plot_surface(Theta / np.pi, Phi / np.pi, probabilities[:, :, i], color=colors[i], alpha=transparency[i],
                        rstride=10, cstride=10, label=fr'$\sigma_t$ = {sigma_t[i]} ps')
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
    plt.legend(bbox_to_anchor=(1, 1.11), loc='upper right', fontsize=18)
    ax.view_init(10, -60)
    plt.show()

    xvec = [0.10, 0.10, 0.8]

    fsize = 30  # fontsize for the figures
    mpl.rcParams.update({'font.size': fsize})

    for i in range(len(sigma_t)):
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111)
        bbox = ax.get_position()
        new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
        ax.set_position(new_bbox)
        im1 = ax.imshow(probabilities[:, :, i], cmap='inferno', vmin=0, vmax=1,
                        extent=[theta[0] / np.pi, theta[-1] / np.pi, phi[0] / np.pi, phi[-1] / np.pi],
                        interpolation='nearest', origin='lower',
                        aspect=(theta[0] - theta[-1]) / (phi[0] - phi[-1]))
        ax.set_xlabel(r'Rot. angle $\theta$/$\pi$ rad')
        ax.set_ylabel(r'Deph. angle $\phi$/$\pi$ rad')
        ax.set_title(fr'$\sigma_t$ = {sigma_t[i]} ps', pad=10)
        # ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
        # ax.set_xticklabels(['0.00', '0.25', '0.50', '0.75', '1.00'])
        ax.set_xticks([0, 0.5, 1])
        ax.set_xticklabels(['0.00', '0.50', '1.00'])
        ax.set_yticks([0, 0.5, 1])
        ax.set_yticklabels(['0.00', '0.50', '1.00'])
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
        cbar = fig.colorbar(im1, cax=axins)
        cbar.set_label(label=r'$\left|\mathrm{1}\right\rangle$ measurement prob.', labelpad=10)
        plt.show()

if compOmega:
    tau1 = 0
    tau2 = 220 # ps
    #sigma_t = 0.01876
    sigma_t = 17 # ps
    gap = 0.019*2*np.pi # THz #0.019*2*np.pi
    Omega1 = 2*np.pi*np.asarray([0, 0.01, 0.015, 0.02])

    sigma_w = 0.0011 # THz
    width = 2*np.pi*6*sigma_w

    theta=np.linspace(0,np.pi, 101)
    phi=np.linspace(0,np.pi,101)
    delta=0

    Theta, Phi = np.meshgrid(theta, phi)
    probabilities = np.zeros((len(theta), len(phi), len(Omega1)))

    for i in range(len(Omega1)):
        probabilities[:,:,i] = Uprob2D(tau1, tau2, sigma_t, Omega1[i], width, gap, Theta, Phi, delta)

    cmap = mpl.colormaps['inferno']
    colors = cmap([0.3, 0.5, 0.7, 0.9])
    transparency = [1, 0.85, 0.7, 0.55]

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
    for i in range(len(Omega1)):
        ax.plot_surface(Theta / np.pi, Phi / np.pi, probabilities[:, :, i], color=colors[i], alpha=transparency[i],
                        rstride=10, cstride=10, label=fr'$\Omega$ = {Omega1[i]/(2*np.pi)} THz')
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


    xvec = [0.10, 0.10, 0.8]

    fsize = 30  # fontsize for the figures
    mpl.rcParams.update({'font.size': fsize})

    for i in range(len(Omega1)):
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111)
        bbox = ax.get_position()
        new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
        ax.set_position(new_bbox)
        im1 = ax.imshow(probabilities[:, :, i], cmap='inferno', vmin=0, vmax=1,
                        extent=[theta[0] / np.pi, theta[-1] / np.pi, phi[0] / np.pi, phi[-1] / np.pi],
                        interpolation='nearest', origin='lower',
                        aspect=(theta[0] - theta[-1]) / (phi[0] - phi[-1]))
        ax.set_xlabel(r'Rot. angle $\theta$/$\pi$ rad')
        ax.set_ylabel(r'Deph. angle $\phi$/$\pi$ rad')
        ax.set_title(fr'$\Omega_1$ = {Omega1[i]/(2*np.pi)} THz', pad=10)
        # ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
        # ax.set_xticklabels(['0.00', '0.25', '0.50', '0.75', '1.00'])
        ax.set_xticks([0, 0.5, 1])
        ax.set_xticklabels(['0.00', '0.50', '1.00'])
        ax.set_yticks([0, 0.5, 1])
        ax.set_yticklabels(['0.00', '0.50', '1.00'])
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
        cbar = fig.colorbar(im1, cax=axins)
        cbar.set_label(label=r'$\left|\mathrm{1}\right\rangle$ measurement prob.', labelpad=10)
        plt.show()

if compeps:
    tau1 = 0
    tau2 = 220 # ps
    sigma_t = 17 # ps
    gap = 0.019*2*np.pi # THz
    Omega1 = 0
    sigma_w = np.asarray([0.0001, 0.0005, 0.0011]) # THz
    width = 2*np.pi*6*sigma_w

    theta=np.linspace(0, np.pi, 101)
    phi=np.linspace(0, np.pi,101)
    delta=0

    Theta, Phi = np.meshgrid(theta, phi)
    probabilities = np.zeros((len(theta), len(phi), len(width)))

    for i in range(len(width)):
        probabilities[:,:,i] = Uprob2D(tau1, tau2, sigma_t, Omega1, width[i], gap, Theta, Phi, delta)


    cmap =  mpl.colormaps['inferno']
    colors = cmap([0.9, 0.7, 0.4])
    transparency = [0.7,0.85,1]

    fsize = 14  # fontsize for the figures
    mpl.rcParams.update({'font.size': fsize})
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial"]
    # mpl.rcParams['mathtext.fontset'] = 'custom'
    # mpl.rcParams['mathtext.rm'] = 'Arial'
    # mpl.rcParams['mathtext.it'] = 'Arial:italic'
    # mpl.rcParams['mathtext.bf'] = 'Arial:bold'
    # arial_font = FontProperties(family='Arial')

    xvec = [0.0, -0.05, 1]

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(projection='3d')
    bbox = ax.get_position()
    new_bbox = (bbox.x0+xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    for i in range(len(width)):
        ax.plot_surface(Theta/np.pi, Phi/np.pi, probabilities[:,:,i], color=colors[i], alpha=transparency[i], rstride=10, cstride=10, label=fr'$\epsilon$ = {round(6*sigma_w[i]*sigma_t,2)} $\frac{{1}}{{\sigma_t}}$') #label=fr'$\epsilon$ = {round(6*sigma_w[i]*1000,1)}GHz')
    ax.set_xlabel(r"Rot. angle, "
           "\n" 
           r"$\theta$ ($\pi$ rad)", labelpad=10)
    ax.set_ylabel('Deph. angle, \n$\phi$ ($\pi$ rad)', labelpad=10)
    #ax.set_zlabel(r'$\left|\mathrm{1}\right\rangle$ meas. prob.', labelpad=10)
    ax.set_zlabel(r"$|\left\langle\Psi'^-_t|\Psi^-_t\right\rangle|^2$", labelpad=5)
    #ax.zaxis.label.set_rotation(-90)
    #ax.zaxis.set_rotate_label(False)
    ax.set_zlim([0, 1])
    ax.set_xticks([0,0.5,1])
    ax.set_xticklabels(['0.0', '0.5', '1.0'])
    ax.set_yticks([0,0.5,1])
    ax.set_yticklabels(['0.0', '0.5', '1.0'])
    ax.set_zticks([0,0.5,1])
    ax.set_zticklabels(['0.0', '0.5', '1.0'])
    ax.tick_params(axis='x', pad=-1)
    ax.tick_params(axis='y', pad=-1)
    with rc_context({
        'mathtext.fontset': 'custom',
        'mathtext.rm': 'Arial',
        'mathtext.it': 'Arial:italic',
        'mathtext.bf': 'Arial:bold'
    }):
        plt.legend(bbox_to_anchor=(0.95, 1.11), loc='upper right', fontsize=14)
    ax.view_init(10, -60)
    plt.show()



    xvec = [0.10, 0.10, 0.8]

    fsize = 30  # fontsize for the figures
    mpl.rcParams.update({'font.size': fsize})

    for i in range(len(width)):
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111)
        bbox = ax.get_position()
        new_bbox = (bbox.x0+xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
        ax.set_position(new_bbox)
        im1 = ax.imshow(probabilities[:,:,i], cmap='inferno', vmin=0, vmax=1,
                       extent=[theta[0] / np.pi, theta[-1] / np.pi, phi[0] / np.pi, phi[-1] / np.pi],
                       interpolation='nearest', origin='lower',
                       aspect=(theta[0] - theta[-1]) / (phi[0] - phi[-1]))
        ax.set_xlabel(r'Rot. angle $\theta$/$\pi$ rad')
        ax.set_ylabel(r'Deph. angle $\phi$/$\pi$ rad')
        ax.set_title(fr'$\epsilon$ = {round(6*sigma_w[i]*1000,1)} GHz', pad=10)
        # ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
        # ax.set_xticklabels(['0.00', '0.25', '0.50', '0.75', '1.00'])
        ax.set_xticks([0, 0.5, 1])
        ax.set_xticklabels(['0.00', '0.50', '1.00'])
        ax.set_yticks([0, 0.5, 1])
        ax.set_yticklabels(['0.00', '0.50', '1.00'])
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
        cbar = fig.colorbar(im1, cax=axins)
        cbar.set_label(r"$|\left\langle\Psi'^-_t|\Psi^-_t\right\rangle|^2$", labelpad=10)
        plt.show()

if compmu:
    tau1 = 0
    tau2 = 220 # ps
    sigma_t = 17 # ps
    gap = 2*np.pi*np.asarray([0.003, 0.01, 0.019, 1]) # THz #0.019*2*np.pi
    Omega1 = 0
    sigma_w = 0.0005 # THz
    width = 2*np.pi*6*sigma_w

    theta=np.linspace(0,np.pi, 101)
    phi=np.linspace(0,np.pi,101)
    delta=0

    Theta, Phi = np.meshgrid(theta, phi)
    probabilities = np.zeros((len(theta), len(phi), len(gap)))

    for i in range(len(gap)):
        probabilities[:,:,i] = Uprob2D(tau1, tau2, sigma_t, Omega1, width, gap[i], Theta, Phi, delta)


    cmap =  mpl.colormaps['inferno']
    colors = cmap([0.3, 0.5, 0.7, 0.9])
    transparency = [1, 0.85, 0.7, 0.55]

    fsize = 22  # fontsize for the figures
    mpl.rcParams.update({'font.size': fsize})
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial"]

    xvec = [0.0, -0.05, 1]

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(projection='3d')
    bbox = ax.get_position()
    new_bbox = (bbox.x0+xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    for i in range(len(gap)):
        ax.plot_surface(Theta/np.pi, Phi/np.pi, probabilities[:,:,i], color=colors[i], alpha=transparency[i], rstride=10, cstride=10, label=fr'$\mu$= {gap[i]/(2*np.pi)} THz')
    ax.set_xlabel(r"Rot. angle, "
           "\n" 
           r"$\theta$ ($\pi$ rad)", labelpad=25)
    ax.set_ylabel('Deph. angle, \n$\phi$ ($\pi$ rad)', labelpad=25)
    ax.set_zlabel(r"$|\left\langle\Psi'^-_t|\Psi^-_t\right\rangle|^2$", labelpad=10)
    ax.set_zlim([0, 1])
    ax.set_xticks([0,0.5,1])
    ax.set_xticklabels(['0.0', '0.5', '1.0'])
    ax.set_yticks([0,0.5,1])
    ax.set_yticklabels(['0.0', '0.5', '1.0'])
    ax.set_zticks([0,0.5,1])
    ax.set_zticklabels(['0.0', '0.5', '1.0'])
    #plt.legend(bbox_to_anchor=(-0.5, 0.5), loc='center left', fontsize=18)
    plt.legend(bbox_to_anchor=(1, 1.16), loc='upper right', fontsize=18)
    ax.view_init(10, -60)
    plt.show()

    xvec = [0.10, 0.10, 0.8]

    fsize = 30  # fontsize for the figures
    mpl.rcParams.update({'font.size': fsize})

    for i in range(len(gap)):
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111)
        bbox = ax.get_position()
        new_bbox = (bbox.x0+xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
        ax.set_position(new_bbox)
        im1 = ax.imshow(probabilities[:,:,i], cmap='inferno', vmin=0, vmax=1,
                       extent=[theta[0] / np.pi, theta[-1] / np.pi, phi[0] / np.pi, phi[-1] / np.pi],
                       interpolation='nearest', origin='lower',
                       aspect=(theta[0] - theta[-1]) / (phi[0] - phi[-1]))
        ax.set_xlabel(r'Rot. angle $\theta$/$\pi$ rad')
        ax.set_ylabel(r'Deph. angle $\phi$/$\pi$ rad')
        ax.set_title(fr'$\mu$ = {gap[i]/(2*np.pi)} THz', pad=10)
        # ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
        # ax.set_xticklabels(['0.00', '0.25', '0.50', '0.75', '1.00'])
        ax.set_xticks([0, 0.5, 1])
        ax.set_xticklabels(['0.00', '0.50', '1.00'])
        ax.set_yticks([0, 0.5, 1])
        ax.set_yticklabels(['0.00', '0.50', '1.00'])
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
        cbar = fig.colorbar(im1, cax=axins)
        cbar.set_label(label=r'$\left|\mathrm{1}\right\rangle$ measurement prob.', labelpad=10)
        plt.show()

if compsigma:
    tau1 = 0
    tau2 = 220 # ps
    sigma_t = np.asarray([10, 20, 30, 40]) # ps
    gap = 0.019*2*np.pi # THz
    Omega1 = 0

    theta=np.linspace(0, np.pi, 101)
    phi=0
    delta=0
    width = 2*np.pi*np.linspace(0.0005, 0.002,101)

    Theta, Width = np.meshgrid(theta, width)
    probabilities = np.zeros((len(theta), len(width), len(sigma_t)))
    keyrates = np.zeros((len(theta), len(width), len(sigma_t)))

    for i in range(len(sigma_t)):
        probabilities[:,:,i] = Uprob2D(tau1, tau2, sigma_t[i], Omega1, Width, gap, Theta, phi, delta)

    cmap =  mpl.colormaps['inferno']

    colors = cmap([0.9, 0.7, 0.5, 0.3])
    transparency = [0.55, 0.7, 0.85, 1]

    fsize = 22  # fontsize for the figures
    mpl.rcParams.update({'font.size': fsize})
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial"]

    xvec = [0.0, -0.05, 1]

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(projection='3d')
    bbox = ax.get_position()
    new_bbox = (bbox.x0+xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    for i in range(len(sigma_t)):
        ax.plot_surface(Theta/np.pi, Width*1000/(2*np.pi), probabilities[:,:,i], color=colors[i], alpha=transparency[i], rstride=10, cstride=10, label=fr'$\sigma_t$ = {sigma_t[i]} ps')
    ax.set_xlabel(r'Rot. angle, $\theta$/$\pi$ rad', labelpad=15)
    ax.set_ylabel('Noise width, \n$\epsilon$/GHz', labelpad=25)
    #ax.set_zlabel(r'$\left|\mathrm{1}\right\rangle$ meas. prob.', labelpad=10)
    ax.set_zlabel(r"$|\left\langle\Psi'^-_t|\Psi^-_t\right\rangle|^2$", labelpad=10)
    #ax.zaxis.label.set_rotation(-90)
    #ax.zaxis.set_rotate_label(False)
    ax.set_zlim([0, 1])
    ax.set_xticks([0,0.5,1])
    ax.set_xticklabels(['0.0', '0.5', '1.0'])
    ax.set_yticks([0.5, 1, 1.5, 2])
    ax.set_yticklabels(['0.5', '1.0', '1.5', '2.0'])
    ax.set_zticks([0,0.5,1])
    ax.set_zticklabels(['0.0', '0.5', '1.0'])
    #plt.legend(bbox_to_anchor=(-0.5, 0.5), loc='center left', fontsize=18)
    plt.legend(bbox_to_anchor=(1, 1.11), loc='upper right', fontsize=18)
    ax.view_init(10, -60)
    plt.show()


    fsize = 14  # fontsize for the figures
    mpl.rcParams.update({'font.size': fsize})
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial"]

    xvec = [0.0, -0.05, 1]

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(projection='3d')
    bbox = ax.get_position()
    new_bbox = (bbox.x0+xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    for i in range(len(sigma_t)):
        ax.plot_surface(Theta / np.pi, Width * 1000 / (2 * np.pi), probabilities[:, :, i], color=colors[i],
                        alpha=transparency[i], rstride=10, cstride=10, label=fr'$\sigma_t$ = {sigma_t[i]} ps')
    ax.set_xlabel(r"Rot. angle, "
           "\n" 
           r"$\theta$ ($\pi$ rad)", labelpad=10)
    ax.set_ylabel('Noise width, \n$\epsilon$/GHz', labelpad=10)
    #ax.set_zlabel(r'$\left|\mathrm{1}\right\rangle$ meas. prob.', labelpad=10)
    ax.set_zlabel(r"$|\left\langle\Psi'^-_t|\Psi^-_t\right\rangle|^2$", labelpad=5)
    #ax.zaxis.label.set_rotation(-90)
    #ax.zaxis.set_rotate_label(False)
    ax.set_zlim([0, 1])
    ax.set_xticks([0,0.5,1])
    ax.set_xticklabels(['0.0', '0.5', '1.0'])
    ax.set_yticks([0.5, 1, 1.5, 2])
    ax.set_yticklabels(['0.5', '1.0', '1.5', '2.0'])
    ax.set_zticks([0,0.5,1])
    ax.set_zticklabels(['0.0', '0.5', '1.0'])
    ax.tick_params(axis='x', pad=-1)
    ax.tick_params(axis='y', pad=-1)
    #plt.legend(bbox_to_anchor=(-0.5, 0.5), loc='center left', fontsize=18)
    plt.legend(bbox_to_anchor=(0.95, 1.07), loc='upper right', fontsize=12)
    ax.view_init(10, 40)
    plt.show()



    tau1 = 0
    tau2 = 220 # ps
    #sigma_t = np.asarray([10, 20, 30, 40]) # ps
    gap = 0.019*2*np.pi # THz #0.019*2*np.pi
    Omega1 = 0

    # theta=np.linspace(0, np.pi, 101)
    theta=np.pi/2
    phi=0
    delta=0
    width = 2*np.pi*np.linspace(0.0005, 0.002,101)
    sigma_t = np.linspace(10, 40, 101)


    Width, Sigma = np.meshgrid(width, sigma_t)
    probabilities = np.zeros((len(width), len(sigma_t)))
    keyrates = np.zeros((len(width), len(sigma_t)))


    probabilities = Uprob2D(tau1, tau2, Sigma, Omega1, Width, gap, theta, phi, delta)

    cmap =  mpl.colormaps['inferno']
    colors = cmap([0.9, 0.7, 0.4])

    colors = cmap([0.9, 0.7, 0.5, 0.3])
    transparency = [0.55, 0.7, 0.85, 1]

    fsize = 22  # fontsize for the figures
    mpl.rcParams.update({'font.size': fsize})
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial"]

    xvec = [0.0, -0.05, 1]

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(projection='3d')
    bbox = ax.get_position()
    new_bbox = (bbox.x0+xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    ax.plot_surface(Width*1000/(2*np.pi), Sigma, probabilities, color=colors[1], alpha=transparency[3], rstride=10, cstride=10)
    ax.set_xlabel('Noise width, \n$\epsilon$/GHz', labelpad=25)
    ax.set_ylabel(r'Time width, $\sigma_t$/ps', labelpad=15)
    #ax.set_zlabel(r'$\left|\mathrm{1}\right\rangle$ meas. prob.', labelpad=10)
    ax.set_zlabel(r"$|\left\langle\Psi'^-_t|\Psi^-_t\right\rangle|^2$", labelpad=10)
    #ax.zaxis.label.set_rotation(-90)
    #ax.zaxis.set_rotate_label(False)
    ax.set_zlim([0, 1])
    ax.set_xticks([0.5, 1, 1.5, 2])
    ax.set_xticklabels(['0.5', '1.0', '1.5', '2.0'])
    ax.set_yticks([10,20,30,40])
    ax.set_yticklabels(['10', '20', '30', '40'])
    ax.set_zticks([0,0.5,1])
    ax.set_zticklabels(['0.0', '0.5', '1.0'])
    #plt.legend(bbox_to_anchor=(-0.5, 0.5), loc='center left', fontsize=18)
    plt.legend(bbox_to_anchor=(1, 1.11), loc='upper right', fontsize=18)
    ax.view_init(10, -60)
    plt.show()

if savedat:
    tau1 = 0
    tau2 = 220  # ps
    sigma_t = 17  # ps
    gap = 0.019 * 2 * np.pi  # THz #0.019*2*np.pi
    Omega1 = 0
    sigma_w = 0.0001 # THz
    width = 2 * np.pi * 6 * sigma_w

    theta = np.linspace(0, 2 * np.pi, 1200)
    phi = 0
    delta = 0

    overlaps = np.zeros((len(theta)), dtype=complex)
    mags = np.zeros((len(theta)))
    phases = np.zeros((len(theta)))

    overlaps = Uoverlap2D(tau1, tau2, sigma_t, Omega1, width, gap, theta, phi, delta)


    probabilities = np.abs(overlaps)**2
    phases = np.arctan(np.divide(np.imag(overlaps),np.real(overlaps)))



    np.savetxt(
        'meas_probs_500MHz_phi0.csv',
        probabilities, delimiter=',')
    print("saved")

    np.savetxt(
        'meas_phases_500MHz_phi0.csv',
        phases, delimiter=',')
    print("saved")

    np.savetxt(
        '500MHz_phi0_thetas.csv',
        theta, delimiter=',')
    print("saved")

    fsize = 22  # fontsize for the figures
    mpl.rcParams.update({'font.size': fsize})
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial"]
    cmap = mpl.colormaps['inferno']

    xvec = [0.0, -0.05, 1]

    plt.plot(theta, probabilities)
    plt.xlabel('theta')
    plt.ylabel('eta0')
    plt.show()

    plt.plot(theta, phases)
    plt.show()

if phasefac:
    tau1 = 0
    tau2 = 220 # ps
    sigma_t = 17 # ps
    gap = 0.019*2*np.pi # THz #0.019*2*np.pi
    Omega1 = 0
    sigma_w = np.asarray([0.0001, 0.0005, 0.0011]) # THz
    width = 2*np.pi*6*sigma_w

    theta=np.linspace(0,2*np.pi, 501)
    phi=np.linspace(0,2*np.pi,501)
    delta=0

    Theta, Phi = np.meshgrid(theta, phi)

    overlaps = np.zeros((len(theta), len(phi), len(width)), dtype=complex)
    mags = np.zeros((len(theta), len(phi), len(width)))
    phases = np.zeros((len(theta), len(phi), len(width)))

    for i in range(len(width)):
        overlaps[:,:,i] = Uoverlap2D(tau1, tau2, sigma_t, Omega1, width[i], gap, Theta, Phi, delta)


    mags = np.abs(overlaps)**2
    phases = np.arctan(np.divide(np.imag(overlaps),np.real(overlaps)))

    if not krcomp:
        cmap = mpl.colormaps['inferno']
        colors = cmap([0.3, 0.6, 0.9])
        transparency = [1, 0.7, 0.55]

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
            ax.plot_surface(Theta / np.pi, Phi / np.pi, phases[:,:,i]/np.pi, color=colors[i], alpha=transparency[i],rstride=10, cstride=10, label=fr'$\sigma_\omega$ = {sigma_w[i]*1000} GHz')
            ax.set_xlabel(r'Rot. angle, $\theta$/$\pi$ rad', labelpad=15)
        ax.set_ylabel('Deph. angle, \n$\phi$/$\pi$ rad', labelpad=25)
        ax.set_zlabel(r"$\Psi'^-_t$ global phase/$\pi$ rad", labelpad=10)
        #ax.set_zlim([0, 1])
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(['0.0', '1.0', '2.0'])
        ax.set_yticks([0, 1, 2])
        ax.set_yticklabels(['0.0', '1.0', '2.0'])
        plt.legend(bbox_to_anchor=(1, 1.16), loc='upper right', fontsize=18)
        ax.view_init(10, -60)
        plt.show()

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(projection='3d')
        bbox = ax.get_position()
        new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
        ax.set_position(new_bbox)
        for i in range(len(width)):
            ax.plot_surface(Theta / np.pi, Phi / np.pi, mags[:,:,i]/np.pi, color=colors[i], alpha=transparency[i],rstride=10, cstride=10, label=fr'$\sigma_\omega$ = {sigma_w[i]*1000} GHz')
        ax.set_xlabel(r'Rot. angle, $\theta$/$\pi$ rad', labelpad=15)
        ax.set_ylabel('Deph. angle, \n$\phi$/$\pi$ rad', labelpad=25)
        ax.set_zlabel(r"$\Psi'^-_t$ global phase/$\pi$ rad", labelpad=10)
        #ax.set_zlim([0, 1])
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(['0.0', '1.0', '2.0'])
        ax.set_yticks([0, 1, 2])
        ax.set_yticklabels(['0.0', '1.0', '2.0'])
        plt.legend(bbox_to_anchor=(1, 1.16), loc='upper right', fontsize=18)
        ax.view_init(10, -60)
        plt.show()


        xvec = [0.10, 0.10, 0.8]

        fsize = 30  # fontsize for the figures
        mpl.rcParams.update({'font.size': fsize})

        for i in range(len(width)):
            fig = plt.figure(figsize=(8, 8))
            ax = fig.add_subplot(111)
            bbox = ax.get_position()
            new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
            ax.set_position(new_bbox)
            im1 = ax.imshow(phases[:,:,i]/np.pi, cmap='inferno',
                            extent=[theta[0] / np.pi, theta[-1] / np.pi, phi[0] / np.pi, phi[-1] / np.pi],
                            interpolation='nearest', origin='lower',
                            aspect=(theta[0] - theta[-1]) / (phi[0] - phi[-1]))
            ax.set_xlabel(r'Rot. angle $\theta$ ($\pi$ rad)')
            ax.set_ylabel(r'Deph. angle $\phi$ ($\pi$ rad)')
            ax.set_title(fr'$\epsilon$ = {round(6*sigma_w[i]*1000,1)} GHz', pad=10)
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
            cbar = fig.colorbar(im1, cax=axins)
            cbar.set_label(label=r"$\Psi'^-_t$ global phase ($\pi$ rad)", labelpad=10)
            plt.show()


if (krcomp and phasefac):
    #probability to key rate conversion
    file = pd.ExcelFile('AllKeyRates_SimulatedData.xlsx')
    ###INTERPOLATION FUCNTION###
    # Asymmetric losses keyrate as a function of asymmetric loss (for zero loss channel)
    Asym_loss = pd.read_excel(file, sheet_name='OurProt_AmpDamp_fine_smallgap')
    params = pd.read_excel(file, sheet_name='OurProt_AmpDamp_fine_params')
    eta0 = params['eta0'].to_numpy()
    delta = params['delta'].to_numpy()
    # we want to interpolate the data with a given key rate to get a given eta0

    Eta, Delta = np.meshgrid(eta0, delta)

    num_pts=126

    asymm_loss = np.zeros((num_pts, num_pts))
    for i in range(num_pts):
        asymm_loss[:,i] = Asym_loss[f'column {i + 1}']

    print(eta0, delta)
    print(np.shape(asymm_loss) == (len(eta0), len(delta)))
    print(mags[1,1,0], phases[1,1,0])


    interp_func_spline_xy = interpolate.RectBivariateSpline(eta0, delta, asymm_loss)
    print("test2")
    print(interp_func_spline_xy.ev(mags[1,1,0], np.abs(phases[1,1,0])))


    epsilon_keyrate = np.zeros((len(theta), len(phi), len(width)))
    epsilon_keyrate_mags = np.zeros((len(theta), len(phi), len(width)))

    for i in range(len(width)):
        #Mag, Pha = np.meshgrid(mags[:,:,i], np.abs(phases[:,:,i]), indexing='ij')
        epsilon_keyrate[:,:,i] = interp_func_spline_xy.ev(mags[:,:,i],phases[:,:,i])  # use np.abs(phases) because the keyrate is calculated for 0 to pi/2 phase, but it is symmetric about phi=0
        epsilon_keyrate_mags[:,:,i] = interp_func_spline_xy.ev(mags[:,:,i], np.zeros(np.shape(mags[:,:,i])))
    #keyrate for our protocol

    difference =  epsilon_keyrate_mags[0, :, 1]-epsilon_keyrate[0,:,1] #calculating the difference between the keyrates when the global phase are considered vs when we neglect them
    print("diff = ", difference)

    print("KR=", epsilon_keyrate)

    xvec = [0.10, 0.10, 0.8]


    #### PLOTTING THE KEY RATE COMPARISON AGAINST THETA
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial"]

    fsize = 30  # fontsize for the figures
    mpl.rcParams.update({'font.size': fsize})


    vmin_vec = [0.416, 0.127, 0.0]
    tickvals = [[0.425, 0.45, 0.475, 0.5], [0.2, 0.3, 0.4, 0.5], [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]]

    for i in range(len(width)):
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111)
        bbox = ax.get_position()
        new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
        ax.set_position(new_bbox)
        im1 = ax.imshow(epsilon_keyrate[:,:,i], cmap='inferno', vmin=vmin_vec[i], vmax=0.5,
                        extent=[theta[0] / np.pi, theta[-1] / np.pi, phi[0] / np.pi, phi[-1] / np.pi],
                        interpolation='nearest', origin='lower',
                        aspect=(theta[0] - theta[-1]) / (phi[0] - phi[-1]))
        ax.set_xlabel(r'Rot. angle, $\theta$ ($\pi$ rad)')
        ax.set_ylabel(r'Deph. angle, $\phi$ ($\pi$ rad)')
        ax.set_title(fr'$\epsilon$ = {round(6*sigma_w[i]*1000,1)} GHz', pad=10)
        # ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
        # ax.set_xticklabels(['0.00', '0.25', '0.50', '0.75', '1.00'])
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(['0.0', '1.0', '2.0'])
        ax.set_yticks([0, 1, 2])
        ax.set_yticklabels(['0.0', '1.0', '2.0'])
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
        cbar = fig.colorbar(im1, cax=axins)
        cbar.set_label(label=r'Key rate (bits/pulse)', labelpad=10)
        cbar.set_ticks(tickvals[i])
        # cbar.set_ticklabels(['0.00', '0.25', '0.50'])
        plt.show()


    cmap =  mpl.colormaps['inferno']
    colors = cmap([0.9, 0.7, 0.4])
    transparency = [0.7,0.85,1]

    fsize = 22  # fontsize for the figures
    mpl.rcParams.update({'font.size': fsize})
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial"]

    xvec = [0.0, -0.05, 1]

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(projection='3d')
    bbox = ax.get_position()
    new_bbox = (bbox.x0+xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    for i in range(len(width)):
        ax.plot_surface(Theta / np.pi, Phi / np.pi, epsilon_keyrate[:,:,i], color=colors[i], alpha=transparency[i], rstride=10,cstride=10, label=fr'$\sigma_\omega$ = {sigma_w[i]*1000} GHz')
    ax.set_xlabel(r'Rot. angle, $\theta$/$\pi$ rad', labelpad=15)
    ax.set_ylabel('Deph. angle, \n$\phi$/$\pi$ rad', labelpad=25)
    ax.set_zlabel(r'Key rate', labelpad=10)
    # ax.set_zlim([0, 0.25])
    # ax.set_xticks([0,0.5, 1])
    # ax.set_xticklabels(['0.0', '0.5', '1.0'])
    # ax.set_yticks([0,0.5, 1])
    # ax.set_yticklabels(['0.0', '0.5', '1.0'])
    # ax.set_zticks([0,0.25])
    # ax.set_zticklabels(['0.0', '0.25'])
    plt.legend(bbox_to_anchor=(1, 1.11), loc='upper right', fontsize=18)
    ax.view_init(10, -60)
    plt.show()

    bb84 = pd.read_excel(file, sheet_name='BB84_theta')
    wang = pd.read_excel(file, sheet_name='Wang_theta')
    boileau4 = pd.read_excel(file, sheet_name='Boileau4_theta')
    li_deph = pd.read_excel(file, sheet_name='LiDeph_theta')

    cmap = mpl.colormaps['inferno']
    transparency = [1, 0.85, 0.7, 0.55]

    colors = cmap([0.1, 0.3, 0.5, 0.9])

    fsize = 18  # fontsize for the figures

    #### PLOTTING THE KEY RATE COMPARISON AGAINST THETA
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial"]

    xvec = [0.1, 0.03, 0.82]
    #xvec = [0.1, -0.02, 0.95]

    # colors = cmap([0.8, 0.6, 0.4, 0.0])
    colors = cmap([0.85, 0.6, 0.4, 0.0])

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    bbox = ax.get_position()
    new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    ax.plot(bb84['theta']/np.pi, bb84['keyrate'], color='#ffba7c', linewidth=2, label='BB84', linestyle='-')
    ax.plot(wang['theta']/np.pi, wang['keyrate'], color='tab:orange', linewidth=2, label='Wang', linestyle='-')
    #ax.plot(li_deph['theta']/np.pi, li_deph['keyrate'], color=colors[2], linewidth=2, label='Li')
    ax.plot(boileau4['theta'] / np.pi, boileau4['keyrate'], color = '#863f00', linewidth = 2, label = 'Boileau et al.', linestyle='-')
    ax.plot(theta / (np.pi), epsilon_keyrate[0, :, 0], color= 'tab:blue', linewidth=3,linestyle='-',
            label=fr'This work, $\epsilon$ = {round(6*sigma_w[0]*sigma_t,2)} $\frac{{1}}{{\sigma_t}}$')
    ax.plot(theta / (np.pi), epsilon_keyrate[0, :, 1], color= '#4da4e0', linewidth=3, linestyle='--',
            label=fr'This work, $\epsilon$ = {round(6*sigma_w[1]*sigma_t,2)} $\frac{{1}}{{\sigma_t}}$')
    ax.plot(theta / (np.pi), epsilon_keyrate[0, :, 2], color= '#89c3eb', linewidth=3, linestyle=':',
            label=fr'This work, $\epsilon$ = {round(6*sigma_w[2]*sigma_t,2)} $\frac{{1}}{{\sigma_t}}$')
    ax.set_xlabel(r'Rotation angle, $\theta$ ($\pi$ rad)')
    ax.set_ylabel('Key rate (bits/pulse)')
    ax.set_aspect(1.2)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(['0.00', '0.25', '0.50', '0.75', '1.00'])
    ax.tick_params(axis='x', pad=7)
    with rc_context({
        'mathtext.fontset': 'custom',
        'mathtext.rm': 'Arial',
        'mathtext.it': 'Arial:italic',
        'mathtext.bf': 'Arial:bold'
    }):
        ll = ax.legend(ncol=2, bbox_to_anchor=(0.5, 1.38), loc='upper center', fontsize=16, handlelength=1.5, labelspacing=0, handleheight=2.1)
        #ll = ax.legend(ncol=2, bbox_to_anchor=(0.5, 1.25), loc='upper center', fontsize=20, handlelength=1.5, labelspacing=0, handleheight=2.1)
        for text in ll.texts:
            # text.set_va('bottom')
            text_x, text_y = text.get_position()
            text.set_position((text_x, text_y + 3))
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 0.53])
    plt.show()

    # #9e4a00



