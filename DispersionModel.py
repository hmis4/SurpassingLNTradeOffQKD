import numpy as np
import math
import cmath
import matplotlib.pyplot as plt
from scipy.special import erf
import matplotlib as mpl
import pandas as pd
import scipy.interpolate as interpolate
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def probff1(a, sigma):
    return np.exp(-a**2*sigma**2)

def probtt1(a, sigma):
    return np.exp(-a**2/sigma**2)

def probtf1(a, Omega1, Omega2, sigma_w, tau1, tau2, sigma_t):
    coeff = 8*sigma_t**2*sigma_w**2/(1+sigma_w**2*sigma_t**2)**2
    normt = 1/(1-np.exp(-(tau1-tau2)**2/(2*sigma_t**2)))
    normf = 1 / (1 - np.exp(-(Omega1 - Omega2) ** 2 / (2 * sigma_w ** 2)))
    freqdecay = np.exp(-sigma_t**2*(Omega1**2 + Omega2**2)/(1+sigma_w**2*sigma_t**2))
    timedecay = np.exp(-sigma_w**2*((a+tau1)**2 + (a+tau2)**2)/(1+sigma_w**2*sigma_t**2))
    osc = 1 - np.cos((Omega2 - Omega1)*(tau2-tau1)/(1+sigma_w**2*sigma_t**2))
    return coeff*normt*normf*freqdecay*timedecay*osc

# def overff1(a, sigma, Omega1, Omega2):
#     return np.exp(-a**2*sigma**2/2)*np.exp(1j*a*(Omega1+Omega2))

# def overtt1(a, sigma):
#     return np.exp(-a**2/(2*sigma**2))

def phaseff1(a, Omega1, Omega2):
    return a*(Omega1+Omega2)

def probff2(b, Omega1, Omega2, sigma_w):
    normf = 1 / (1 - np.exp(-(Omega1 - Omega2) ** 2 / (2 * sigma_w ** 2)))**2
    coeff = 1/(1+b**2*sigma_w**4)
    freqdecay = np.exp(-2*b**2*sigma_w**2*(Omega1**2+Omega2**2)/(1+b**2*sigma_w**4))
    oscdecay = 1 + np.exp(-(Omega1-Omega2)**2/((1+b**2*sigma_w**4)*sigma_w**2)) - 2*np.exp(-(Omega1-Omega2)**2/(2*(1+b**2*sigma_w**4)*sigma_w**2))*np.cos(b*(Omega1-Omega2)**2/(2*(1+b**2*sigma_w**4)))
    return normf*coeff*freqdecay*oscdecay

def probtt2(b, tau1, tau2, sigma_t):
    normt = 1 / (1 - np.exp(-(tau1 - tau2) ** 2 / (2 * sigma_t ** 2)))**2
    coeff = sigma_t**4/(sigma_t**4 + b**2)
    oscdecay = 1 + np.exp(-sigma_t**2*(tau1-tau2)**2/(sigma_t**4 + b**2)) - 2*np.exp(-sigma_t**2*(tau1-tau2)**2/(2*(sigma_t**4 + b**2)))*np.cos(b*(tau1-tau2)**2/(2*(sigma_t**4+b**2)))
    return normt*coeff*oscdecay

def probtf2(b, Omega1, Omega2, sigma_w, tau1, tau2, sigma_t):
    coeff = 4*sigma_t**2*sigma_w**2/((1+sigma_w**2*sigma_t**2)**2 + 4*sigma_w**4*b**2)
    normt = 1/(1-np.exp(-(tau1-tau2)**2/(2*sigma_t**2)))
    normf = 1 / (1 - np.exp(-(Omega1 - Omega2) ** 2 / (2 * sigma_w ** 2)))
    freqdecay = np.exp(-(sigma_t**2*(1+sigma_t**2*sigma_w**2)+4*sigma_w**2*b**2)*(Omega1**2 + Omega2**2)/((1+sigma_w**2*sigma_t**2)**2 + 4*sigma_w**4*b**2))
    timedecay = np.exp(-sigma_w**2*(1+sigma_w**2*sigma_t**2)*(tau1**2 + tau2**2)/((1+sigma_w**2*sigma_t**2)**2+4*sigma_w**4*b**2))
    oscdecay = np.exp(-4*sigma_w**2*b*(Omega1*tau1+Omega2*tau2)/((1+sigma_w**2*sigma_t**2)**2 + 4*sigma_w**4*b**2)) + np.exp(-4*sigma_w**2*b*(Omega1*tau2+Omega2*tau1)/((1+sigma_w**2*sigma_t**2)**2 + 4*sigma_w**4*b**2)) - 2*np.exp(-2*sigma_w**2*b*(Omega1+Omega2)*(tau1+tau2)/((1+sigma_w**2*sigma_t**2)**2 + 4*sigma_w**4*b**2))*np.cos((1+sigma_w**2*sigma_t**2)*(Omega2 - Omega1)*(tau2-tau1)/((1+sigma_w**2*sigma_t**2)**2 + 4*sigma_w**4*b**2))
    return coeff*normt*normf*freqdecay*timedecay*oscdecay

def A(mu1, mu2, sigma):
    return 1/np.sqrt(2*(1-np.exp(-(mu1-mu2)**2/(2*sigma**2))))

# def overlapff1d(b, Omega1, Omega2, sigma_w):
#     return np.sqrt(1/(1-1j*b*sigma_w**2))*np.exp(-(Omega1**2+Omega2**2)/(2*sigma_w**2))*np.exp((Omega1+Omega2)**2/(4*(1-1j*b*sigma_w**2)*sigma_w**2))
#
# def overlapff2d (b, Omega1, Omega2, sigma_w):
#     return 2*A(Omega1, Omega2, sigma_w)**2*(overlapff1d(b, Omega1, Omega1, sigma_w)*overlapff1d(b, Omega2, Omega2, sigma_w) - overlapff1d(b, Omega1, Omega2, sigma_w)*overlapff1d(b, Omega2, Omega1, sigma_w))


def overlapff2d(b, Omega1, Omega2, sigma_w):
    return 2*A(Omega1, Omega2, sigma_w)**2*(1/(1-1j*b*sigma_w**2))*np.exp(-(Omega1**2 + Omega2**2)/sigma_w**2)*(np.exp((Omega1**2 + Omega2**2)/((1-1j*b*sigma_w**2)*sigma_w**2))-np.exp((Omega1+Omega2)**2/(2*(1-1j*b*sigma_w**2)*sigma_w**2)))


def overlaptt2d(b, tau1, tau2, sigma_t):
    return 2*A(tau1, tau2, sigma_t)**2*(sigma_t**2/(sigma_t**2-1j*b))*(1 - np.exp(-(tau1-tau2)**2/(2*(sigma_t**2-1j*b))))

def overlap1d(a, tau0, tau1, sigma):
    return sigma*1/np.sqrt(sigma**2-1j*a)*np.exp(-(tau0-tau1)**2/(4*(sigma**2-1j*a)))

def prob1d(a, tau0, tau1, sigma):
    return np.abs(overlap1d(a, tau0, tau1, sigma))**2






example = False
compwidth = False
compOmega2 = False
compOmega1 = False
keyrate1 = False
keyrate2 = False

# probability to key rate conversion
file = pd.ExcelFile(
    'AllKeyRates_SimulatedData.xlsx')
###INTERPOLATION FUCNTION###
# Asymmetric losses keyrate as a function of asymmetric loss (for zero loss channel)
Asym_loss = pd.read_excel(file, sheet_name='OurProt_AmpDamp_fine_smallgap')
params = pd.read_excel(file, sheet_name='OurProt_AmpDamp_fine_params')
eta0 = params['eta0'].to_numpy()
delta = params['delta'].to_numpy()
# we want to interpolate the data with a given key rate to get a given eta0

Eta, Delta = np.meshgrid(eta0, delta)

num_pts = 126

asymm_loss = np.zeros((num_pts, num_pts))
for i in range(num_pts):
    asymm_loss[:, i] = Asym_loss[f'column {i + 1}']

interp_func_spline_xy = interpolate.RectBivariateSpline(eta0, delta, asymm_loss)


fsize = 28  # fontsize for the figures
mpl.rcParams.update({'font.size': fsize})
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]

if example:
    tau1 = 0
    tau2 = 220 # ps
    sigma_t = 17 # ps
    Omega1 = 0 # THz - offset from central frequency of tau1
    Omega2 = 0.019*2*np.pi # THz #0.019*2*np.pi
    sigma_w = 2*np.pi*0.0011 # THz

    # tau1 = 0
    # tau2 = 800 # ps
    # sigma_t = 120 # ps
    # Omega1 = 0 # THz - offset from central frequency of tau1
    # Omega2 = 0.019*2*np.pi # THz #0.019*2*np.pi
    # sigma_w = 0.0011*2*np.pi


    a = np.linspace(0, 100, 1001)

    prob1f = probff1(a,sigma_w)
    prob1t = probtt1(a,sigma_t)
    prob1tf = probtf1(a, Omega1, Omega2, sigma_w, tau1, tau2, sigma_t)
    print(prob1tf)

    b = np.linspace(0, 2000, 1001)

    prob2f = probff2(b, Omega1, Omega2, sigma_w)
    prob2t = probtt2(b, tau1, tau2, sigma_t)
    prob2tf = probtf2(b, Omega1, Omega2, sigma_w, tau1, tau2, sigma_t)

    print('prod tf 2')
    print(prob2tf)


    fsize=28 # fontsize for the figures
    mpl.rcParams.update({'font.size': fsize})
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial"]

    cmap =  mpl.colormaps['inferno']
    #colors = cmap([0.7, 0.5, 0.0])
    colors = cmap([0.0, 0.55, 0.7])

    xvec = [0.05, 0.07, 0.95]

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    bbox = ax.get_position()
    new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    ax.plot(a, prob1f, color = colors[0], linewidth=3, label=r"$|\left\langle\Psi'^-_f|\Psi^-_f\right\rangle|^2$")
    ax.plot(a, prob1t, color = colors[1], linewidth=3, label=r"$|\left\langle\Psi'^-_t|\Psi^-_t\right\rangle|^2$")
    ax.plot(a, prob1tf, color = colors[2], linewidth=3, label=r"$|\left\langle\Psi'^-_t|\Psi^-_f\right\rangle|^2$")
    ax.set_xlabel(r'$\alpha_1$ (ps)')
    ax.set_ylabel('Fidelity')
    # plt.title('Title')
    ax.legend(bbox_to_anchor=(0.99, 0.08), loc='lower right', fontsize=18)
    ax.set_aspect(60)
    #ax.set_aspect(1000000)
    ax.set_xlim([a[0],a[-1]])
    #ax.set_ylim([0,1])
    plt.show()


    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    bbox = ax.get_position()
    new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    ax.plot(b, prob2f, color = colors[0], linewidth=3, label=r"$|\left\langle\Psi'^-_f|\Psi^-_f\right\rangle|^2$")
    ax.plot(b, prob2t, color = colors[1], linewidth=3, label=r"$|\left\langle\Psi'^-_t|\Psi^-_t\right\rangle|^2$")
    ax.plot(b, prob2tf, color = colors[2], linewidth=3, label=r"$|\left\langle\Psi'^-_t|\Psi^-_f\right\rangle|^2$")
    ax.set_xlabel(r'$\alpha_2$ (ps$^2$)')
    ax.set_ylabel(r"Fidelity")
    # plt.title('Title')
    ax.legend(bbox_to_anchor=(0.98, 0.98), loc='upper right', fontsize=18)
    ax.set_aspect(1200)
    #ax.set_aspect(1000000)
    ax.set_xlim([b[0],b[-1]])
    #ax.set_ylim([0,1])
    plt.show()


if compwidth:
    tau1 = 0
    tau2 = 320 # ps
    sigma_t = np.asarray([10,20,30,40]) # ps
    Omega1 = 0 # THz - offset from central frequency of tau1
    Omega2 = 0.08*2*np.pi # THz #0.019*2*np.pi
    sigma_w = np.asarray([0.0001,0.0005, 0.001, 0.005])*2*np.pi

    #a = np.linspace(0, 100, 1001)
    b = np.linspace(0, 2000, 2001)

    # A1, SW1 = np.meshgrid(a, sigma_w)
    # A2, ST1 = np.meshgrid(a, sigma_t)
    B1, SW2 = np.meshgrid(b, sigma_w)
    B2, ST2 = np.meshgrid(b, sigma_t)

    #
    # prob1f = probff1(A1, SW1)
    # prob1t = probtt1(A2, ST1)
    prob2f = probff2(B1, Omega1, Omega2, SW2)
    prob2t = probtt2(B2, tau1, tau2, ST2)


    fsize=28 # fontsize for the figures
    mpl.rcParams.update({'font.size': fsize})
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial"]

    cmap =  mpl.colormaps['inferno']
    #colors = cmap([0.7, 0.5, 0.0])
    colors = cmap([0.9, 0.7, 0.5, 0.0])

    #colors = ['#ffba7c', '#ff9537', '#d16200', '#9e4a00']
    colors = ['#E18FAD','#B46783','#864059', '#59182F']

    xvec = [0.05, 0.07, 0.95]

    # fig = plt.figure(figsize=(8, 8))
    # ax = fig.add_subplot(111)
    # bbox = ax.get_position()
    # new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    # ax.set_position(new_bbox)
    # for i in range(len(sigma_w)):
    #     ax.plot(a, prob1f[i, :], color = colors[i], linewidth=2, label=f"$P(\Psi^-_f|\Psi'^-_f)$, $\sigma_\omega =${sigma_w[i]}")
    # ax.set_xlabel(r'$\alpha_1$/ps')
    # ax.set_ylabel('Measurement probability')
    # ax.legend(bbox_to_anchor=(0.95, 0.05), loc='lower right', fontsize=14)
    # ax.set_aspect(60)
    # ax.set_xlim([a[0],a[-1]])
    # #ax.set_ylim([0,1])
    # plt.show()
    #
    # fig = plt.figure(figsize=(8, 8))
    # ax = fig.add_subplot(111)
    # bbox = ax.get_position()
    # new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    # ax.set_position(new_bbox)
    # for i in range(len(sigma_t)):
    #     ax.plot(a, prob1t[i,:], color = colors[i], label=f"$P(\Psi^-_f|\Psi'^-_f)$, $\sigma_t =${sigma_t[i]}")
    # ax.set_xlabel(r'$\alpha_1$/ps')
    # ax.set_ylabel('Measurement probability')
    # ax.legend(bbox_to_anchor=(0.95, 0.05), loc='lower right', fontsize=14)
    # ax.set_aspect(60)
    # ax.set_xlim([a[0],a[-1]])
    # #ax.set_ylim([0,1])
    # plt.show()

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    bbox = ax.get_position()
    new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    for i in range(len(sigma_w)):
        ax.plot(b, prob2f[i,:], color = colors[i], linewidth=3, label=f"$\sigma_\omega =${sigma_w[i]/(2*np.pi)*1000} GHz")
    ax.set_xlabel(r'$\alpha_2$ (ps$^2$)')
    ax.set_ylabel(r"$|\left\langle\Psi'^-_f|\Psi^-_f\right\rangle|^2$")
    ax.legend(bbox_to_anchor=(0.98, 0.46), loc='upper right', fontsize=18)
    ax.set_aspect(1200)
    ax.set_xlim([b[0],b[-1]])
    #ax.set_ylim([0,1])
    plt.show()

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    bbox = ax.get_position()
    new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    for i in range(len(sigma_t)):
        ax.plot(b, prob2t[i,:], color = colors[i], linewidth=3, label=f"$\sigma_t =${sigma_t[i]} ps")
    ax.set_xlabel(r'$\alpha_2$ (ps$^2$)')
    ax.set_ylabel(r"$|\left\langle\Psi'^-_t|\Psi^-_t\right\rangle|^2$")
    ax.legend(bbox_to_anchor=(0.98, 0.98), loc='upper right', fontsize=18)
    ax.set_aspect(1200)
    ax.set_xlim([b[0],b[-1]])
    #ax.set_ylim([0,1])
    plt.show()

if compOmega2:
    tau1 = 0
    tau2 = np.asarray([200, 600, 1000]) # ps
    sigma_t = 20 # ps
    Omega1 = 0 # THz - offset from central frequency of tau1
    Omega2 = np.asarray([0.01, 0.05, 0.1])*2*np.pi # THz #0.019*2*np.pi
    sigma_w = 0.001*2*np.pi

    b = np.linspace(0, 2000, 1001)

    B1, W2 = np.meshgrid(b, Omega2)
    B2, T2 = np.meshgrid(b, tau2)

    prob2f = probff2(B1, Omega1, W2, sigma_w)
    prob2t = probtt2(B2, tau1, T2, sigma_t)

    cmap =  mpl.colormaps['inferno']
    #colors = cmap([0.7, 0.5, 0.0])
    colors = cmap([0.9, 0.7, 0.5, 0.0])

    colors = ['#E18FAD', '#9D546E', '#59182F']

    xvec = [0.05, 0.07, 0.95]

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    bbox = ax.get_position()
    new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    for i in range(len(Omega2)):
        ax.plot(b, prob2f[i,:], color = colors[i], linewidth=3, label=f"$\Omega_1 =${Omega2[i]/(2*np.pi)} THz")
    ax.set_xlabel(r'$\alpha_2$ (ps$^2$)')
    ax.set_ylabel(r"$|\left\langle\Psi'^-_f|\Psi^-_f\right\rangle|^2$")
    ax.legend(bbox_to_anchor=(0.98, 0.98), loc='upper right', fontsize=18)
    ax.set_aspect(1200)
    ax.set_xlim([b[0],b[-1]])
    #ax.set_ylim([0,1])
    plt.show()

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    bbox = ax.get_position()
    new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    for i in range(len(tau2)):
        ax.plot(b, prob2t[i,:], color = colors[i], linewidth=3, label=fr"$\tau_1=${tau2[i]} ps")
    ax.set_xlabel(r'$\alpha_2$ (ps$^2$)')
    ax.set_ylabel(r"$|\left\langle\Psi'^-_t|\Psi^-_t\right\rangle|^2$")
    ax.legend(bbox_to_anchor=(0.98, 0.98), loc='upper right', fontsize=18)
    ax.set_aspect(1200)
    ax.set_xlim([b[0],b[-1]])
    #ax.set_ylim([0,1])
    plt.show()

if compOmega1:
    tau1 = np.asarray([0, 50, 100])
    tau2 = 300 # ps
    sigma_t = 20 # ps
    Omega1 = np.asarray([0, 0.005, 0.01])*2*np.pi # THz - offset from central frequency of tau1
    Omega2 = 0.02*2*np.pi # THz #0.019*2*np.pi
    sigma_w = 0.001*2*np.pi

    b = np.linspace(0, 2000, 1001)

    B1, W1 = np.meshgrid(b, Omega1)
    B2, T1 = np.meshgrid(b, tau1)

    prob2f = probff2(B1, W1, Omega2, sigma_w)
    prob2t = probtt2(B2, T1, tau2, sigma_t)

    cmap =  mpl.colormaps['inferno']
    #colors = cmap([0.7, 0.5, 0.0])
    colors = cmap([0.9, 0.7, 0.5, 0.0])

    colors = ['#E18FAD', '#9D546E', '#59182F']

    xvec = [0.05, 0.07, 0.95]

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    bbox = ax.get_position()
    new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    for i in range(len(Omega1)):
        ax.plot(b, prob2f[i,:], color = colors[i], linewidth=3, label=f"$\Omega_0 =${Omega1[i]/(2*np.pi)} THz")
    ax.set_xlabel(r'$\alpha_2$ (ps$^2$)')
    ax.set_ylabel(r"$|\left\langle\Psi'^-_f|\Psi^-_f\right\rangle|^2$")
    ax.legend(bbox_to_anchor=(0.98, 0.98), loc='upper right', fontsize=18)
    ax.set_aspect(1200)
    ax.set_xlim([b[0],b[-1]])
    #ax.set_ylim([0,1])
    plt.show()

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    bbox = ax.get_position()
    new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    for i in range(len(tau1)):
        ax.plot(b, prob2t[i,:], color = colors[i], linewidth=3, label=fr"$\tau_0=${tau1[i]} ps")
    ax.set_xlabel(r'$\alpha_2$ (ps$^2$)')
    ax.set_ylabel(r"$|\left\langle\Psi'^-_t|\Psi^-_f\right\rangle|^2$")
    ax.legend(bbox_to_anchor=(0.98, 0.98), loc='upper right', fontsize=18)
    ax.set_aspect(1200)
    ax.set_xlim([b[0],b[-1]])
    #ax.set_ylim([0,1])
    plt.show()


if keyrate1:
    tau1 = 0
    tau2 = 220 # ps
    sigma_t = np.asarray([10,17,30]) # ps
    Omega1 = 0 # THz - offset from central frequency of tau1
    Omega2 = 0.019*2*np.pi # THz #0.019*2*np.pi
    sigma_w = 0.0011*2*np.pi

    a = np.linspace(0, 100, 1001)
    #b = np.linspace(0, 2000, 2001)

    prob1f = probff1(a, sigma_w)
    prob1t = np.zeros((len(a),len(sigma_t)))
    for i in range(len(sigma_t)):
        prob1t[:,i] = probtt1(a, sigma_t[i])
    phase1f = phaseff1(a, Omega1, Omega2)
    # prob2f = probff2(B1, Omega1, Omega2, SW2)
    # prob2t = probtt2(B2, tau1, tau2, ST2)

    for i in range(len(phase1f)):
        if (phase1f[i] > np.pi/2) and (phase1f[i] < np.pi):
            phase1f[i] = np.pi - phase1f[i]

        elif (phase1f[i] > np.pi) and (phase1f[i] < 3*np.pi/2):
            phase1f[i] = phase1f[i] - np.pi

        elif (phase1f[i] > 3*np.pi/2) and (phase1f[i] < 2*np.pi):
            phase1f[i] = 2*np.pi - phase1f[i]

        elif (phase1f[i] > 2*np.pi) and (phase1f[i] < 5*np.pi/2):
            phase1f[i] = phase1f[i] - 2*np.pi

        elif (phase1f[i] > 5*np.pi/2) and (phase1f[i] < 3*np.pi):
            phase1f[i] = 3*np.pi - phase1f[i]

        elif (phase1f[i] > 3*np.pi) and (phase1f[i] < 7*np.pi/2):
            phase1f[i] = phase1f[i] - 3*np.pi


        elif (phase1f[i] > 7* np.pi / 2) and (phase1f[i] < 4* np.pi):
            phase1f[i] = 4* np.pi - phase1f[i]


    # epsilon_keyrate = interp_func_spline_xy.ev(prob1t, phase1f)
    #
    # epsilon_keyrate_2 = interp_func_spline_xy.ev(np.divide(prob1t, prob1f), phase1f)

    epsilon_keyrate_3 = np.zeros((len(a), len(sigma_t)))

    for i in range(len(sigma_t)):
        epsilon_keyrate_3[:,i] = prob1f*interp_func_spline_xy.ev(np.divide(prob1t[:,i], prob1f), phase1f)

    fsize=22 # fontsize for the figures
    mpl.rcParams.update({'font.size': fsize})
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial"]

    cmap =  mpl.colormaps['inferno']
    #colors = cmap([0.7, 0.5, 0.0])
    colors = cmap([0.9, 0.7, 0.5, 0.0])

    colors = ['#E18FAD', '#9D546E', '#59182F']

    xvec = [0.05, 0.05, 0.95]


    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    bbox = ax.get_position()
    new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    for i in range(len(sigma_t)):
        ax.plot(a, epsilon_keyrate_3[:,i], color=colors[i], linewidth=3, label=f'$\sigma_t$ = {sigma_t[i]} ps')
    #ax.plot(a, phase1f/np.pi, color=colors[3], linewidth=2)
    ax.set_xlabel(r'Lin. dispersion parameter, $\alpha_1$ (ps)')
    ax.set_ylabel("Key rate (bits/pulse)")
    ax.legend(bbox_to_anchor=(0.98, 0.98), loc='upper right', fontsize=20)
    ax.set_aspect(125)
    ax.set_xlim([a[0],a[-1]])
    ax.set_ylim([-0.005, 0.505])
    plt.show()




    tau1 = 0
    tau2 = 220 # ps
    #sigma_t = 17 #ps
    sigma_t = 30 # ps
    Omega1 = 0 # THz - offset from central frequency of tau1
    Omega2 = 2*np.pi*np.asarray([0.01, 0.015, 0.02]) # THz #0.019*2*np.pi
    sigma_w = 0.0011*2*np.pi

    a = np.linspace(0, 100, 1001)
    # b = np.linspace(0, 2000, 2001)

    prob1f = probff1(a, sigma_w)
    prob1t = probtt1(a, sigma_t)
    phase1f = np.zeros((len(a),len(Omega2)))
    for i in range(len(Omega2)):
        phase1f[:, i] = phaseff1(a, Omega1, Omega2[i])
    # prob2f = probff2(B1, Omega1, Omega2, SW2)
    # prob2t = probtt2(B2, tau1, tau2, ST2)

    for j in range(len(Omega2)):
        for i in range(len(a)):
            if (phase1f[i][j] > np.pi / 2) and (phase1f[i][j] < np.pi):
                phase1f[i][j] = np.pi - phase1f[i][j]

            elif (phase1f[i][j] > np.pi) and (phase1f[i][j] < 3 * np.pi / 2):
                phase1f[i][j] = phase1f[i][j] - np.pi

            elif (phase1f[i][j] > 3 * np.pi / 2) and (phase1f[i][j] < 2 * np.pi):
                phase1f[i][j] = 2 * np.pi - phase1f[i][j]

            elif (phase1f[i][j] > 2 * np.pi) and (phase1f[i][j] < 5 * np.pi / 2):
                phase1f[i][j] = phase1f[i][j] - 2 * np.pi

            elif (phase1f[i][j] > 5 * np.pi / 2) and (phase1f[i][j] < 3 * np.pi):
                phase1f[i][j] = 3 * np.pi - phase1f[i][j]

            elif (phase1f[i][j] > 3 * np.pi) and (phase1f[i][j] < 7 * np.pi / 2):
                phase1f[i][j] = phase1f[i][j] - 3 * np.pi


            elif (phase1f[i][j] > 7 * np.pi / 2) and (phase1f[i][j] < 4 * np.pi):
                phase1f[i][j] = 4 * np.pi - phase1f[i][j]

    # epsilon_keyrate = interp_func_spline_xy.ev(prob1t, phase1f)
    #
    # epsilon_keyrate_2 = interp_func_spline_xy.ev(np.divide(prob1t, prob1f), phase1f)

    epsilon_keyrate_3 = np.zeros((len(a), len(Omega2)))

    for i in range(len(Omega2)):
        epsilon_keyrate_3[:, i] = prob1f * interp_func_spline_xy.ev(np.divide(prob1t, prob1f), phase1f[:,i])

    fsize = 22  # fontsize for the figures
    mpl.rcParams.update({'font.size': fsize})
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial"]

    cmap = mpl.colormaps['inferno']
    # colors = cmap([0.7, 0.5, 0.0])
    colors = cmap([0.9, 0.7, 0.5, 0.0])

    colors = ['#E18FAD', '#9D546E', '#59182F']

    xvec = [0.05, 0.05, 0.95]

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    bbox = ax.get_position()
    new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    for i in range(len(Omega2)):
        ax.plot(a, epsilon_keyrate_3[:, i], color=colors[i], linewidth=3, label=f'$\Omega_1$ = {round(Omega2[i]/(2*np.pi),3)} THz')
    # ax.plot(a, phase1f/np.pi, color=colors[3], linewidth=2)
    ax.set_xlabel(r'Lin. dispersion parameter, $\alpha_1$ (ps)')
    ax.set_ylabel("Key rate (bits/pulse)")
    ax.legend(bbox_to_anchor=(0.98, 0.98), loc='upper right', fontsize=20)
    ax.set_aspect(125)
    ax.set_xlim([a[0], a[-1]])
    ax.set_ylim([-0.005, 0.505])
    plt.show()



if keyrate2:
    tau1 = 0
    tau2 = 220 # ps
    #sigma_t = 17 #ps
    sigma_t = np.asarray([10,17,30]) # ps
    Omega1 = 0 # THz - offset from central frequency of tau1
    Omega2 = 0.019*2*np.pi # THz #0.019*2*np.pi
    sigma_w = 0.0011*2*np.pi


    #a = np.linspace(0, 100, 1001)
    b = np.linspace(0, 2000, 20001)
    #
    # prob1f = probff1(a, sigma_w)
    # prob1t = np.zeros((len(a),len(sigma_t)))
    # for i in range(len(sigma_t)):
    #     prob1t[:,i] = probtt1(a, sigma_t[i])
    # phase1f = phaseff1(a, Omega1, Omega2)



    overlap2f = np.zeros(len(b), dtype=complex)
    mags2f = np.zeros(len(b))
    phase2f = np.zeros(len(b))
    overlap2t = np.zeros((len(b),len(sigma_t)), dtype=complex)
    mags2t = np.zeros((len(b),len(sigma_t)))
    phase2t = np.zeros((len(b),len(sigma_t)))


    overlap2f = overlapff2d(b, Omega1, Omega2, sigma_w)

    for i in range(len(sigma_t)):
        overlap2t[:,i] = overlaptt2d(b, tau1, tau2, sigma_t[i])

    print('overlap')
    print(overlap2f)

    mags2f = np.abs(overlap2f)**2
    phase2f = np.arctan(np.divide(np.imag(overlap2f),np.real(overlap2f)))
    mags2t = np.abs(overlap2t)**2
    phase2t = np.arctan(np.divide(np.imag(overlap2t),np.real(overlap2t)))

    print(mags2f)

    phasediff = np.zeros((len(b), len(sigma_t)))

    for i in range(len(sigma_t)):
        phasediff[:,i] = phase2f - phase2t[:,i]

    plt.plot(b, mags2f, label='mag 2f')
    plt.plot(b, phase2f, label='phase 2f')
    plt.plot(b, mags2t[:,1], label='mag 2t')
    plt.plot(b, phase2t[:,1], label='phase 2t')
    plt.plot(b, phasediff[:,1], label='phase diff')
    # plt.plot(b, prob2f, label='orig f')
    # plt.plot(b, prob2t, label='orig t')
    plt.legend()
    plt.show()


    for j in range(len(sigma_t)):
        for i in range(len(b)):
            if (phasediff[i,j]<0):
                phasediff[i,j] = -phasediff[i,j]

            if (phasediff[i,j] > np.pi/2) and (phasediff[i,j] < np.pi):
                phasediff[i,j] = np.pi - phasediff[i,j]

            elif (phasediff[i,j] > np.pi) and (phasediff[i,j] < 3*np.pi/2):
                phasediff[i,j] = phasediff[i,j] - np.pi

            elif (phasediff[i,j] > 3*np.pi/2) and (phasediff[i,j] < 2*np.pi):
                phasediff[i,j] = 2*np.pi - phasediff[i,j]

            elif (phasediff[i,j] > 2*np.pi) and (phasediff[i,j] < 5*np.pi/2):
                phasediff[i,j] = phasediff[i,j] - 2*np.pi

            elif (phasediff[i,j] > 5*np.pi/2) and (phasediff[i,j] < 3*np.pi):
                phasediff[i,j] = 3*np.pi - phasediff[i,j]

            elif (phasediff[i,j] > 3*np.pi) and (phasediff[i,j] < 7*np.pi/2):
                phasediff[i,j] = phasediff[i,j] - 3*np.pi

            elif (phasediff[i,j] > 7* np.pi / 2) and (phasediff[i,j] < 4* np.pi):
                phasediff[i,j] = 4* np.pi - phasediff[i,j]


    # epsilon_keyrate = interp_func_spline_xy.ev(prob1t, phase1f)
    #
    # epsilon_keyrate_2 = interp_func_spline_xy.ev(np.divide(prob1t, prob1f), phase1f)

    # epsilon_keyrate_3 = np.zeros((len(a), len(sigma_t)))

    epsilon_keyrate = np.zeros((len(b),len(sigma_t)))

    for j in range(len(sigma_t)):
        for i in range(len(b)):
            if mags2f[i] >= mags2t[i,j]:
                epsilon_keyrate[i,j] = mags2f[i]*interp_func_spline_xy.ev(np.divide(mags2t[i,j], mags2f[i]), phasediff[i,j])

            else:
                epsilon_keyrate[i,j] = mags2t[i,j] * interp_func_spline_xy.ev(np.divide(mags2f[i], mags2t[i,j]), phasediff[i,j])


    fsize=22 # fontsize for the figures
    mpl.rcParams.update({'font.size': fsize})
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial"]

    cmap =  mpl.colormaps['inferno']
    #colors = cmap([0.7, 0.5, 0.0])
    colors = cmap([0.9, 0.7, 0.5, 0.0])

    colors = ['#E18FAD', '#9D546E', '#59182F']

    xvec = [0.05, 0.05, 0.95]


    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    bbox = ax.get_position()
    new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    for i in range(len(sigma_t)):
        ax.plot(b, epsilon_keyrate[:,i], color=colors[i], linewidth=3,  label=f'$\sigma_t$ = {sigma_t[i]} ps')
    #ax.plot(a, phase1f/np.pi, color=colors[3], linewidth=2)
    #ax.plot(b, prob1d(b, tau1, tau2, sigma_t[1]))
    ax.set_xlabel(r'Quad. dispersion parameter, $\alpha_2$ (ps$^2$)')
    ax.set_ylabel("Key rate (bits/pulse)")
    ax.legend(bbox_to_anchor=(0.98, 0.98), loc='upper right', fontsize=20)
    ax.set_aspect(2500)
    ax.set_xlim([b[0],b[-1]])
    ax.set_ylim([-0.005, 0.505])
    plt.show()


    tau1 = 0
    tau2 = 220 # ps
    #sigma_t = 17 #ps
    sigma_t = 17 # ps
    Omega1 = 0 # THz - offset from central frequency of tau1
    Omega2 = 2*np.pi*np.asarray([0.01, 0.015, 0.02]) # THz #0.019*2*np.pi
    sigma_w = 0.0011*2*np.pi


    #a = np.linspace(0, 100, 1001)
    b = np.linspace(0, 2000, 20001)
    #
    # prob1f = probff1(a, sigma_w)
    # prob1t = np.zeros((len(a),len(sigma_t)))
    # for i in range(len(sigma_t)):
    #     prob1t[:,i] = probtt1(a, sigma_t[i])
    # phase1f = phaseff1(a, Omega1, Omega2)



    overlap2t = np.zeros(len(b), dtype=complex)
    mags2t = np.zeros(len(b))
    phase2t = np.zeros(len(b))
    overlap2f = np.zeros((len(b),len(Omega2)), dtype=complex)
    mags2f = np.zeros((len(b),len(Omega2)))
    phase2f = np.zeros((len(b),len(Omega2)))



    overlap2t = overlaptt2d(b, tau1, tau2, sigma_t)

    for i in range(len(Omega2)):
        overlap2f[:,i] = overlapff2d(b, Omega1, Omega2[i], sigma_w)

    print('overlap')
    print(overlap2f)

    mags2f = np.abs(overlap2f)**2
    phase2f = np.arctan(np.divide(np.imag(overlap2f),np.real(overlap2f)))
    mags2t = np.abs(overlap2t)**2
    phase2t = np.arctan(np.divide(np.imag(overlap2t),np.real(overlap2t)))

    print(mags2f)

    phasediff = np.zeros((len(b), len(Omega2)))

    for i in range(len(Omega2)):
        phasediff[:,i] = phase2f[:,i] - phase2t



    for j in range(len(Omega2)):
        for i in range(len(b)):
            if (phasediff[i,j]<0):
                phasediff[i,j] = -phasediff[i,j]

            if (phasediff[i,j] > np.pi/2) and (phasediff[i,j] < np.pi):
                phasediff[i,j] = np.pi - phasediff[i,j]

            elif (phasediff[i,j] > np.pi) and (phasediff[i,j] < 3*np.pi/2):
                phasediff[i,j] = phasediff[i,j] - np.pi

            elif (phasediff[i,j] > 3*np.pi/2) and (phasediff[i,j] < 2*np.pi):
                phasediff[i,j] = 2*np.pi - phasediff[i,j]

            elif (phasediff[i,j] > 2*np.pi) and (phasediff[i,j] < 5*np.pi/2):
                phasediff[i,j] = phasediff[i,j] - 2*np.pi

            elif (phasediff[i,j] > 5*np.pi/2) and (phasediff[i,j] < 3*np.pi):
                phasediff[i,j] = 3*np.pi - phasediff[i,j]

            elif (phasediff[i,j] > 3*np.pi) and (phasediff[i,j] < 7*np.pi/2):
                phasediff[i,j] = phasediff[i,j] - 3*np.pi

            elif (phasediff[i,j] > 7* np.pi / 2) and (phasediff[i,j] < 4* np.pi):
                phasediff[i,j] = 4* np.pi - phasediff[i,j]


    # epsilon_keyrate = interp_func_spline_xy.ev(prob1t, phase1f)
    #
    # epsilon_keyrate_2 = interp_func_spline_xy.ev(np.divide(prob1t, prob1f), phase1f)

    # epsilon_keyrate_3 = np.zeros((len(a), len(sigma_t)))

    epsilon_keyrate = np.zeros((len(b),len(Omega2)))

    for j in range(len(Omega2)):
        for i in range(len(b)):
            if mags2f[i,j] >= mags2t[i]:
                epsilon_keyrate[i,j] = mags2f[i,j]*interp_func_spline_xy.ev(np.divide(mags2t[i], mags2f[i,j]), phasediff[i,j])

            else:
                epsilon_keyrate[i,j] = mags2t[i] * interp_func_spline_xy.ev(np.divide(mags2f[i,j], mags2t[i]), phasediff[i,j])


    fsize=22 # fontsize for the figures
    mpl.rcParams.update({'font.size': fsize})
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial"]

    cmap =  mpl.colormaps['inferno']
    #colors = cmap([0.7, 0.5, 0.0])
    colors = cmap([0.9, 0.7, 0.5, 0.0])

    colors = ['#E18FAD', '#9D546E', '#59182F']

    xvec = [0.05, 0.05, 0.95]


    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    bbox = ax.get_position()
    new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    for i in range(len(Omega2)):
        ax.plot(b, epsilon_keyrate[:,i], color=colors[i], linewidth=3,  label=f'$\Omega_1$ = {round(Omega2[i]/(2*np.pi),3)} THz')
    #ax.plot(a, phase1f/np.pi, color=colors[3], linewidth=2)
    #ax.plot(b, prob1d(b, tau1, tau2, sigma_t[1]))
    ax.set_xlabel(r'Quad. dispersion parameter, $\alpha_2$ (ps$^2$)')
    ax.set_ylabel("Key rate (bits/pulse)")
    ax.legend(bbox_to_anchor=(0.98, 0.98), loc='upper right', fontsize=20)
    ax.set_aspect(2500)
    ax.set_xlim([b[0],b[-1]])
    ax.set_ylim([-0.005, 0.505])
    plt.show()



