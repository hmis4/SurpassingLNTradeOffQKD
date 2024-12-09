import numpy as np
import math
import cmath
import matplotlib.pyplot as plt
from scipy.special import erf
import matplotlib as mpl
import pandas as pd
import scipy.interpolate as interpolate
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def A(mu0, mu1, sigma):
    normfactor = 2*np.pi*sigma**2*(1 - np.exp(-((mu0-mu1)**2)/(2*sigma**2)))
    return 1/np.sqrt(normfactor)

def lin1Doverlaptt(a, tau0, tau1, sigma):
    return np.exp(-(a+tau0-tau1)**2/(4*sigma**2))

def lin1Doverlapff(a, Omega0, Omega1, sigma):
    return np.exp(-(Omega1-Omega0)**2/(4*sigma**2))*np.exp(-(a**2*sigma**2)/4)*np.exp(1j*a*(Omega0+Omega1)/2)

def overlaptt1(a, delta, tau0, tau1, sigma):
    return (lin1Doverlaptt(a, tau0, tau0, sigma)*lin1Doverlaptt(a+delta, tau1, tau1, sigma)-lin1Doverlaptt(a, tau0, tau1, sigma)*lin1Doverlaptt(a+delta, tau1, tau0, sigma))

def overlapff1(a, delta, Omega0, Omega1, sigma):
    return (lin1Doverlapff(a, Omega0, Omega0, sigma)*lin1Doverlapff(a+delta, Omega1, Omega1, sigma)-lin1Doverlapff(a, Omega0, Omega1, sigma)*lin1Doverlapff(a+delta, Omega1, Omega0, sigma))


def quad1Doverlaptt(a, tau0, tau1, sigma):
    norm = sigma*np.sqrt(1/(sigma**2-1j*a))
    expfacs = np.exp(-(tau0-tau1)**2/(4*(sigma**2-1j*a)))
    return norm*expfacs

def quad1Doverlapff(a, Omega0, Omega1, sigma):
    norm = np.sqrt(1/(1-1j*a*sigma**2))
    expfacs = np.exp(-(Omega0**2+Omega1**2)/(2*sigma**2))*np.exp((Omega0+Omega1)**2/(4*sigma**2*(1-1j*a*sigma**2)))
    return norm*expfacs

def overlaptt2(a, delta, tau0, tau1, sigma):
    return (quad1Doverlaptt(a, tau0, tau0, sigma)*quad1Doverlaptt(a+delta, tau1, tau1, sigma)-quad1Doverlaptt(a, tau0, tau1, sigma)*quad1Doverlaptt(a+delta, tau1, tau0, sigma))

def overlapff2(a, delta, Omega0, Omega1, sigma):
    return (quad1Doverlapff(a, Omega0, Omega0, sigma)*quad1Doverlapff(a+delta, Omega1, Omega1, sigma)-quad1Doverlapff(a, Omega0, Omega1, sigma)*quad1Doverlapff(a+delta, Omega1, Omega0, sigma))

keyrate1 = True
keyrate2 = True

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


if keyrate1:
    tau1 = 0
    tau2 = 220 # ps
    sigma_t = 30 # ps
    Omega1 = 0 # THz - offset from central frequency of tau1
    Omega2 = 0.019*2*np.pi # THz #0.019*2*np.pi
    sigma_w = 0.0011*2*np.pi
    delta = np.asarray([0, 1, 5])
    #delta = np.asarray([0])

    print(lin1Doverlapff(0, Omega2, Omega2, sigma_w))
    print(overlapff1(0, 0, Omega1, Omega2, sigma_w))

    a = np.linspace(0, 100, 1001)
    #b = np.linspace(0, 2000, 2001)

    overlap1f = np.zeros((len(a),len(delta)), dtype=complex)
    overlap1t = np.zeros((len(a), len(delta)),dtype=complex)
    prob1f = np.zeros((len(a), len(delta)))
    phase1f = np.zeros((len(a), len(delta)))
    prob1t = np.zeros((len(a), len(delta)))
    phase1t = np.zeros((len(a), len(delta)))


    for i in range(len(delta)):
        overlap1f[:,i] = overlapff1(a, delta[i], Omega1, Omega2, sigma_w)
        overlap1t[:, i] = overlaptt1(a, delta[i], tau1, tau2, sigma_t)


    prob1f = np.abs(overlap1f)**2
    prob1t = np.abs(overlap1t)**2
    phase1f = np.arctan(np.divide(np.imag(overlap1f),np.real(overlap1f)))
    phase1t = np.arctan(np.divide(np.imag(overlap1t),np.real(overlap1t)))

    print(prob1f)
    print(prob1t)

    phasediff = np.zeros((len(a), len(delta)))
    phasediff = phase1f - phase1t

    print(phasediff)


    for j in range(len(delta)):
        for i in range(len(a)):
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



    epsilon_keyrate = np.zeros((len(a), len(delta)))

    for j in range(len(delta)):
        for i in range(len(a)):
            if prob1f[i,j] >= prob1t[i,j]:
                epsilon_keyrate[i,j] = prob1f[i,j]*interp_func_spline_xy.ev(np.divide(prob1t[i,j], prob1f[i,j]), phasediff[i,j])

            else:
                epsilon_keyrate[i,j] = prob1t[i,j] * interp_func_spline_xy.ev(np.divide(prob1f[i,j], prob1t[i,j]), phasediff[i,j])



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
    for i in range(len(delta)):
        ax.plot(a, epsilon_keyrate[:,i], color=colors[i], linewidth=3, label=f'$\delta$ = {delta[i]} ps')
    #ax.plot(a, phase1f/np.pi, color=colors[3], linewidth=2)
    ax.set_xlabel(r'Lin. dispersion parameter, $\alpha_1$ (ps)')
    ax.set_ylabel("Key rate (bits/pulse)")
    ax.legend(bbox_to_anchor=(0.98, 0.98), loc='upper right', fontsize=20)
    ax.set_aspect(125)
    ax.set_xlim([a[0],a[-1]])
    ax.set_ylim([-0.005, 0.505])
    plt.show()



if keyrate2:
    tau1 = 0
    tau2 = 220 # ps
    sigma_t = 30 # ps
    Omega1 = 0 # THz - offset from central frequency of tau1
    Omega2 = 0.019*2*np.pi # THz #0.019*2*np.pi
    sigma_w = 0.0011*2*np.pi
    delta = np.asarray([0, 10, 50])
    #delta = np.asarray([0])

    print(lin1Doverlapff(0, Omega2, Omega2, sigma_w))
    print(overlapff1(0, 0, Omega1, Omega2, sigma_w))

    b = np.linspace(0, 2000, 2001)
    #b = np.linspace(0, 2000, 2001)

    overlap1f = np.zeros((len(b),len(delta)), dtype=complex)
    overlap1t = np.zeros((len(b), len(delta)),dtype=complex)
    prob1f = np.zeros((len(b), len(delta)))
    phase1f = np.zeros((len(b), len(delta)))
    prob1t = np.zeros((len(b), len(delta)))
    phase1t = np.zeros((len(b), len(delta)))


    for i in range(len(delta)):
        overlap1f[:,i] = overlapff2(b, delta[i], Omega1, Omega2, sigma_w)
        overlap1t[:, i] = overlaptt2(b, delta[i], tau1, tau2, sigma_t)


    prob1f = np.abs(overlap1f)**2
    prob1t = np.abs(overlap1t)**2
    phase1f = np.arctan(np.divide(np.imag(overlap1f),np.real(overlap1f)))
    phase1t = np.arctan(np.divide(np.imag(overlap1t),np.real(overlap1t)))

    print(prob1f)
    print(prob1t)

    phasediff = np.zeros((len(a), len(delta)))
    phasediff = phase1f - phase1t

    print(phasediff)


    for j in range(len(delta)):
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



    epsilon_keyrate = np.zeros((len(b), len(delta)))

    for j in range(len(delta)):
        for i in range(len(b)):
            if prob1f[i,j] >= prob1t[i,j]:
                epsilon_keyrate[i,j] = prob1f[i,j]*interp_func_spline_xy.ev(np.divide(prob1t[i,j], prob1f[i,j]), phasediff[i,j])

            else:
                epsilon_keyrate[i,j] = prob1t[i,j] * interp_func_spline_xy.ev(np.divide(prob1f[i,j], prob1t[i,j]), phasediff[i,j])



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
    for i in range(len(delta)):
        ax.plot(b, epsilon_keyrate[:,i], color=colors[i], linewidth=3, label=f'$\delta$ = {delta[i]} ps')
    #ax.plot(a, phase1f/np.pi, color=colors[3], linewidth=2)
    ax.set_xlabel(r'Quad. dispersion parameter, $\alpha_2$ (ps$^2$)')
    ax.set_ylabel("Key rate (bits/pulse)")
    ax.legend(bbox_to_anchor=(0.98, 0.98), loc='upper right', fontsize=20)
    ax.set_aspect(2500)
    ax.set_xlim([b[0],b[-1]])
    ax.set_ylim([-0.005, 0.505])
    plt.show()

