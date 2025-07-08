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


def overlapshift(tau1, tau2, sigma_t, shift, start, end):
    coeff = 0.5
    phase = np.exp(-1j*shift*(tau1+tau2)/2)
    shiftdecay = np.exp(-sigma_t**2*shift**2/4)
    bindecay = np.exp(-(tau1-tau2)**2/(4*sigma_t**2))
    endval = erf(sigma_t * end + sigma_t * shift/2 + 1j*sigma_t * (tau2-tau1)/(2*sigma_t**2))
    startval = erf(sigma_t * start + sigma_t * shift / 2 + 1j * sigma_t * (tau2 - tau1) / (2 * sigma_t ** 2))
    return coeff*phase*shiftdecay*bindecay*(endval - startval)

def Uoverlap1D(tau1,tau2, sigma_t, Omega, epsilon, mu, theta=0, phi=0):
    f1 = overlapshift(tau1, tau2, sigma_t, 0, -math.inf, Omega-epsilon/2)
    f2 = np.cos(theta)*overlapshift(tau1, tau2, sigma_t, 0, Omega-epsilon/2, Omega+epsilon/2) - np.exp(-1j*phi)*np.sin(theta)*overlapshift(tau1, tau2, sigma_t, mu, Omega-epsilon/2, Omega+epsilon/2)
    f3 = overlapshift(tau1, tau2, sigma_t, 0, Omega + epsilon/2, Omega + mu - epsilon/2)
    f4 = np.cos(theta)*overlapshift(tau1, tau2, sigma_t, 0, Omega+mu-epsilon/2, Omega+mu+epsilon/2) + np.exp(1j*phi)*np.sin(theta)*overlapshift(tau1, tau2, sigma_t, -mu, Omega+mu-epsilon/2, Omega+mu+epsilon/2)
    f5 = overlapshift(tau1, tau2, sigma_t, 0, Omega+mu+epsilon/2, math.inf)
    return f1 + f2 + f3 + f4 + f5

def conjoverlapshift(Omega, tau, sigma_w, sigma_t, shift, start, end):
    coeff = np.sqrt(sigma_t*sigma_w/(2*(1+sigma_w**2*sigma_t**2))) #yes
    phase = np.exp(-1j*(Omega+shift)*tau/(1+sigma_w**2*sigma_t**2)) #yes
    shiftdecay = np.exp(-sigma_t**2*shift**2/(2*(1+sigma_w**2*sigma_t**2))) #yes
    bindecay = np.exp(-(tau**2*sigma_w**2 + Omega**2*sigma_t**2)/(2*(1+sigma_w**2*sigma_t**2))) #yes
    mixdecay = np.exp(-Omega*shift*sigma_t**2/(1+sigma_w**2*sigma_t**2)) #yes
    scaling = np.sqrt((1+sigma_t**2*sigma_w**2)/(2*sigma_w**2))
    offset = (1j*sigma_w**2*tau + shift*sigma_w**2*sigma_t**2 - Omega)/(1+sigma_w**2*sigma_t**2)
    endval = erf(scaling * end + scaling*offset)
    startval = erf(scaling * start + scaling*offset)
    return coeff*phase*shiftdecay*bindecay*mixdecay*(endval - startval)


def conjUoverlap1D(Omega1,tau2, sigma_w,sigma_t, Omega, epsilon, mu, theta=0, phi=0):
    print(Omega1, tau2)
    f1 = conjoverlapshift(Omega1, tau2, sigma_w, sigma_t, 0, -math.inf, Omega-epsilon/2)
    f2 = np.cos(theta)*conjoverlapshift(Omega1, tau2, sigma_w, sigma_t, 0, Omega-epsilon/2, Omega+epsilon/2) - np.exp(-1j*phi)*np.sin(theta)*conjoverlapshift(Omega1, tau2, sigma_w, sigma_t, mu, Omega-epsilon/2, Omega+epsilon/2)
    f3 = conjoverlapshift(Omega1, tau2, sigma_w, sigma_t, 0, Omega + epsilon/2, Omega + mu - epsilon/2)
    f4 = np.cos(theta)*conjoverlapshift(Omega1, tau2, sigma_w, sigma_t, 0, Omega+mu-epsilon/2, Omega+mu+epsilon/2) + np.exp(1j*phi)*np.sin(theta)*conjoverlapshift(Omega1, tau2, sigma_w, sigma_t, -mu, Omega+mu-epsilon/2, Omega+mu+epsilon/2)
    f5 = conjoverlapshift(Omega1, tau2, sigma_w, sigma_t, 0, Omega+mu+epsilon/2, math.inf)
    return f1 + f2 + f3 + f4 + f5

Omega0 = -0.0095
Omega1 = 0.0095
sigma_w = 0.0011
tau0 = -110
tau1 = 110
sigma_t = 17

theta = np.linspace(0, 2*np.pi, 11)
phi =0

print(conjUoverlap1D(Omega0, tau0, sigma_w, sigma_t, Omega0, 6*sigma_w, Omega1-Omega0, theta, phi))

print(erf(1+2j))

print(Uoverlap1D(tau0, tau0, sigma_t, Omega0, 6*sigma_w, Omega1-Omega0, theta, phi))