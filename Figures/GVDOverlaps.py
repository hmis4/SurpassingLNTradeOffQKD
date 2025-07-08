import numpy as np
import math
import cmath
import matplotlib.pyplot as plt
from scipy.special import erf
import matplotlib as mpl
import pandas as pd
import scipy.interpolate as interpolate
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def probff1(a, Omega0, Omega1, sigma):
    return np.exp(-a**2*sigma**2/4)*np.exp(-(Omega1 - Omega0)**2/(4*sigma**2))*np.exp(1j*a*(Omega0 + Omega1))

def probtt1(a, tau0, tau1, sigma):
    return np.exp(-(a + tau1 - tau0)**2/(4*sigma**2))

def probtf1(a, Omega, sigma_w, tau, sigma_t):
    coeff = np.sqrt(2*sigma_t*sigma_w/(1+sigma_w**2*sigma_t**2))
    freqdecay = np.exp(-sigma_t**2*Omega**2/(2*(1+sigma_w**2*sigma_t**2)))

    timedecay = np.exp(-sigma_w**2*(a+tau)**2/(2*(1+sigma_w**2*sigma_t**2)))
    phase = np.exp(1j*Omega*(a+tau)/(1+sigma_w**2*sigma_t**2))
    return coeff*freqdecay*timedecay*phase


def ftoverlap(Omega, sigma_w, tau, sigma_t):
    denom = 1 + sigma_w**2*sigma_t**2
    norm = np.sqrt(2*sigma_w*sigma_t/denom)
    decay = np.exp(-(Omega**2*sigma_t**2 + tau**2*sigma_w**2)/(2*denom))
    phase = np.exp(1j*Omega*tau/denom)
    return norm*decay*phase


print(erf(1+1j))

Omega0 = -0.0095
Omega1 = 0.0095
sigma_w = 0.0011
tau0 = -110
tau1 = 110
sigma_t = 17



####KEYRATE N=1

file2 = pd.ExcelFile(
    'C:/Users/ir22317/OneDrive - University of Bristol/PhD/Noiseless QKD Project/Security Analysis/KeyRates_DispersionChannel_v2.xlsx')
lindisp_30ps = pd.read_excel(file2, sheet_name='LinDisp_1.1GHz_30ps_0.19GHz')
lindisp_17ps = pd.read_excel(file2, sheet_name='LinDisp_1.1GHz_17ps')
lindisp_10ps = pd.read_excel(file2, sheet_name='LinDisp_1.1GHz_10ps')

lindisp_10GHz = pd.read_excel(file2, sheet_name='LinDisp_Omega1_0.01THz')
lindisp_15GHz = pd.read_excel(file2, sheet_name='LinDisp_Omega1_0.015THz')
lindisp_20GHz = pd.read_excel(file2, sheet_name='LinDisp_Omega1_0.02THz')

fsize = 22  # fontsize for the figures
mpl.rcParams.update({'font.size': fsize})
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]

cmap = mpl.colormaps['inferno']
# colors = cmap([0.7, 0.5, 0.0])
colors = cmap([0.9, 0.7, 0.5, 0.0])

colors = ['#E18FAD', '#9D546E', '#59182F']

xvec = [0.05, 0.05, 0.95]


#Vary Sigma_t
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
bbox = ax.get_position()
new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
ax.set_position(new_bbox)
ax.plot(lindisp_10ps['alpha1'], lindisp_10ps['keyrate'], color=colors[0], linewidth=3, label=f'$\sigma_t$ = 10 ps')
ax.plot(lindisp_17ps['alpha1'], lindisp_17ps['keyrate'], color=colors[1], linewidth=3, label=f'$\sigma_t$ = 17 ps')
ax.plot(lindisp_30ps['alpha1'], lindisp_30ps['keyrate'], color=colors[2], linewidth=3, label=f'$\sigma_t$ = 30 ps')
ax.set_xlabel(r'Lin. dispersion parameter, $\alpha_1$ (ps)')
ax.set_ylabel("Key rate (bits/pulse)")
ax.legend(bbox_to_anchor=(0.98, 0.98), loc='upper right', fontsize=20)
ax.set_aspect(125)
ax.set_xlim([0, 100])
ax.set_ylim([-0.005, 0.505])
plt.show()


### Vary Omega1
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
bbox = ax.get_position()
new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
ax.set_position(new_bbox)
ax.plot(lindisp_10GHz['alpha1'], lindisp_10GHz['keyrate'], color=colors[0], linewidth=3, label=f'$\Omega_1$ = 0.01 THz')
ax.plot(lindisp_15GHz['alpha1'], lindisp_15GHz['keyrate'], color=colors[1], linewidth=3, label=f'$\Omega_1$ = 0.015 THz')
ax.plot(lindisp_20GHz['alpha1'], lindisp_20GHz['keyrate'], color=colors[2], linewidth=3, label=f'$\Omega_1$ = 0.02 THz')
ax.set_xlabel(r'Lin. dispersion parameter, $\alpha_1$ (ps)')
ax.set_ylabel("Key rate (bits/pulse)")
ax.legend(bbox_to_anchor=(0.98, 0.98), loc='upper right', fontsize=20)
ax.set_aspect(125)
ax.set_xlim([0, 100])
ax.set_ylim([-0.005, 0.505])
plt.show()



####KEYRATE N=2

quaddisp_30ps = pd.read_excel(file2, sheet_name='QuadDisp_1.1GHz_30ps')
quaddisp_17ps = pd.read_excel(file2, sheet_name='QuadDisp_1.1GHz_17ps')
quaddisp_10ps = pd.read_excel(file2, sheet_name='QuadDisp_1.1GHz_10ps')

quaddisp_10GHz_v1 = pd.read_excel(file2, sheet_name='QuadDisp_Omega1_0.01THz')

quaddisp_10GHz = pd.read_excel(file2, sheet_name='QuadDisp_Omega1_0.01THz_v2')
quaddisp_15GHz = pd.read_excel(file2, sheet_name='QuadDisp_Omega1_0.015THz')
quaddisp_20GHz = pd.read_excel(file2, sheet_name='QuadDisp_Omega1_0.02THz')
quaddisp_13GHz = pd.read_excel(file2, sheet_name='QuadDisp_0.013THz')
quaddisp_12GHz = pd.read_excel(file2, sheet_name='QuadDisp_0.012THz')
quaddisp_10GHz_v2 = pd.read_excel(file2, sheet_name='QuadDisp_0.01THz')

fsize = 22  # fontsize for the figures
mpl.rcParams.update({'font.size': fsize})
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]

cmap = mpl.colormaps['inferno']
# colors = cmap([0.7, 0.5, 0.0])
colors = cmap([0.9, 0.7, 0.5, 0.0])

colors = ['#E18FAD', '#9D546E', '#59182F']

xvec = [0.05, 0.05, 0.95]

#Vary sigma_t
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
bbox = ax.get_position()
new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
ax.set_position(new_bbox)
ax.plot(quaddisp_10ps['alpha2'], quaddisp_10ps['keyrate'], color=colors[0], linewidth=3,label=f'$\sigma_t$ = 10 ps')
ax.plot(quaddisp_17ps['alpha2'], quaddisp_17ps['keyrate'], color=colors[1], linewidth=3,label=f'$\sigma_t$ = 17 ps')
ax.plot(quaddisp_30ps['alpha2'], quaddisp_30ps['keyrate'], color=colors[2], linewidth=3,label=f'$\sigma_t$ = 30 ps')
ax.set_xlabel(r'Quad. dispersion parameter, $\alpha_2$ (ps$^2$)')
ax.set_ylabel("Key rate (bits/pulse)")
ax.legend(bbox_to_anchor=(0.98, 0.98), loc='upper right', fontsize=20)
ax.set_aspect(2500)
ax.set_xlim([0,2000])
#ax.set_xlim([b[0], b[-1]])
ax.set_ylim([-0.005, 0.505])
plt.show()



#Vary Omega_1
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
bbox = ax.get_position()
new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
ax.set_position(new_bbox)
ax.plot(quaddisp_10GHz['alpha2'], quaddisp_10GHz['keyrate'], color=colors[0], linewidth=3,label=f'$\Omega_1$ = 0.01 THz')
# ax.plot(quaddisp_10GHz_v1['alpha2'], quaddisp_10GHz_v1['keyrate'], color='tab:blue', linewidth=2,label=f'$\Omega_1$ = 0.01 THz', linestyle='--')
# ax.scatter(1120, 0.0503048904004472, color='tab:blue')
# ax.scatter(1120, 0.050305078279107, color='tab:orange')
# ax.scatter(1120, 0.050335042074785, color='tab:green', label='UB')
# ax.plot(quaddisp_13GHz['alpha2'], quaddisp_13GHz['keyrate'], color='tab:orange', linewidth=3,label=f'$\Omega_1$ = 0.013 THz')
# ax.plot(quaddisp_12GHz['alpha2'], quaddisp_12GHz['keyrate'], color='tab:blue', linewidth=3,label=f'$\Omega_1$ = 0.012 THz')
#ax.scatter(quaddisp_10GHz_v2['alpha2'], quaddisp_10GHz_v2['keyrate'], color='tab:orange', linewidth=3,label=f'$\Omega_1$ = 0.013 THz')
ax.plot(quaddisp_15GHz['alpha2'], quaddisp_15GHz['keyrate'], color=colors[1], linewidth=3,label=f'$\Omega_1$ = 0.015 THz')
ax.plot(quaddisp_20GHz['alpha2'], quaddisp_20GHz['keyrate'], color=colors[2], linewidth=3,label=f'$\Omega_1$ = 0.02 THz')
ax.set_xlabel(r'Quad. dispersion parameter, $\alpha_2$ (ps$^2$)')
ax.set_ylabel("Key rate (bits/pulse)")
ax.legend(bbox_to_anchor=(0.98, 0.98), loc='upper right', fontsize=20)
ax.set_aspect(2500)
ax.set_xlim([0,2000])
#ax.set_xlim([b[0], b[-1]])
ax.set_ylim([-0.005, 0.505])
plt.show()






####LIN DISP UNEQUAL GVD
lindisp_delta_1ps = pd.read_excel(file2, sheet_name='LinDisp_Delta_1ps_sigma_30ps')
lindisp_delta_5ps = pd.read_excel(file2, sheet_name='LinDisp_Delta_5ps_sigma_30ps')

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
ax.plot(lindisp_30ps['alpha1'], lindisp_30ps['keyrate'], color=colors[0], linewidth=3, label=f'$\delta$ = 0 ps')
ax.plot(lindisp_delta_1ps['alpha1'], lindisp_delta_1ps['keyrate'], color=colors[1], linewidth=3, label=f'$\delta$ = 1 ps')
ax.plot(lindisp_delta_5ps['alpha1'], lindisp_delta_5ps['keyrate'], color=colors[2], linewidth=3, label=f'$\delta$ = 5 ps')
# ax.plot(a, phase1f/np.pi, color=colors[3], linewidth=2)
ax.set_xlabel(r'Lin. dispersion parameter, $\alpha_1$ (ps)')
ax.set_ylabel("Key rate (bits/pulse)")
ax.legend(bbox_to_anchor=(0.98, 0.98), loc='upper right', fontsize=20)
ax.set_aspect(125)
ax.set_xlim([0,100])
ax.set_ylim([-0.005, 0.505])
plt.show()


####QUAD DISP UNEQUAL GVD
quaddisp_delta_10ps = pd.read_excel(file2, sheet_name='QuadDisp_Delta_10ps2_sigma_30ps')
quaddisp_delta_50ps = pd.read_excel(file2, sheet_name='QuadDisp_Delta_50ps2_sigma_30ps')

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
ax.plot(quaddisp_30ps['alpha2'], quaddisp_30ps['keyrate'], color=colors[0], linewidth=3, label=f'$\delta$ = 0 ps$^2$')
ax.plot(quaddisp_delta_10ps['alpha2'], quaddisp_delta_10ps['keyrate'], color=colors[1], linewidth=3, label=f'$\delta$ = 10 ps$^2$')
ax.plot(quaddisp_delta_50ps['alpha2'], quaddisp_delta_50ps['keyrate'], color=colors[2], linewidth=3, label=f'$\delta$ = 50 ps$^2$')
ax.set_xlabel(r'Quad. dispersion parameter, $\alpha_2$ (ps$^2$)')
ax.set_ylabel("Key rate (bits/pulse)")
ax.legend(bbox_to_anchor=(0.98, 0.98), loc='upper right', fontsize=20)
ax.set_aspect(2500)
ax.set_xlim([0,2000])
ax.set_ylim([-0.005, 0.505])
plt.show()
