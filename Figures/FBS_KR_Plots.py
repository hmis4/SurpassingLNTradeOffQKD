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

file = pd.ExcelFile(
    'C:/Users/ir22317/OneDrive - University of Bristol/PhD/Noiseless QKD Project/Security Analysis/AllKeyRates_ForPlots_v3.xlsx')
file2 = pd.ExcelFile(
    'C:/Users/ir22317/OneDrive - University of Bristol/PhD/Noiseless QKD Project/Security Analysis/KeyRates_FBSChannel.xlsx')

bb84 = pd.read_excel(file, sheet_name='BB84_theta')
wang = pd.read_excel(file, sheet_name='Wang_theta')
boileau4 = pd.read_excel(file, sheet_name='Boileau4_theta')
li_deph = pd.read_excel(file, sheet_name='LiDeph_theta')
ourprot_1100MHz = pd.read_excel(file2, sheet_name='FBS_theta_sigma_1.1GHz')
ourprot_500MHz = pd.read_excel(file2, sheet_name='FBS_theta_sigma_0.5GHz')
ourprot_100MHz = pd.read_excel(file2, sheet_name='FBS_theta_sigma_0.1GHz')
#ourprot_v2 = pd.read_excel(file, sheet_name='OurProt_v2_theta_alpha_0.1')
#ourprot_v2_alpha2 = pd.read_excel(file, sheet_name='OurProt_v2_theta_alpha_0.2')

cmap = mpl.colormaps['inferno']
transparency = [1, 0.85, 0.7, 0.55]

colors = cmap([0.1, 0.3, 0.5, 0.9])



fsize = 20  # fontsize for the figures
mpl.rcParams.update({'font.size': fsize})

#### PLOTTING THE KEY RATE COMPARISON AGAINST THETA
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]

xvec = [0.1, 0.03, 0.82]
# xvec = [0.1, -0.02, 0.95]

# colors = cmap([0.8, 0.6, 0.4, 0.0])
colors = cmap([0.85, 0.6, 0.4, 0.0])

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
bbox = ax.get_position()
new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
ax.set_position(new_bbox)
ax.plot(bb84['theta'] / np.pi, bb84['keyrate'], color='#ffba7c', linewidth=2, label='BB84', linestyle='-')
ax.plot(wang['theta'] / np.pi, wang['keyrate'], color='tab:orange', linewidth=2, label='Wang', linestyle='-')
# ax.plot(li_deph['theta']/np.pi, li_deph['keyrate'], color=colors[2], linewidth=2, label='Li')
ax.plot(boileau4['theta'] / np.pi, boileau4['keyrate'], color='#863f00', linewidth=2, label='Boileau et al.',
        linestyle='-')
# ax.plot(theta / (np.pi), epsilon_keyrate[0, :, 0], color='tab:blue', linewidth=3, linestyle='-',
#         label=fr'This work, $\epsilon$ = {round(6 * sigma_w[0] * sigma_t, 2)} $\frac{{1}}{{\sigma_t}}$')
# ax.plot(theta / (np.pi), epsilon_keyrate[0, :, 1], color='#4da4e0', linewidth=3, linestyle='--',
#         label=fr'This work, $\epsilon$ = {round(6 * sigma_w[1] * sigma_t, 2)} $\frac{{1}}{{\sigma_t}}$')
# ax.plot(theta / (np.pi), epsilon_keyrate[0, :, 2], color='#89c3eb', linewidth=3, linestyle=':',
#         label=fr'This work, $\epsilon$ = {round(6 * sigma_w[2] * sigma_t, 2)} $\frac{{1}}{{\sigma_t}}$')
ax.plot(ourprot_100MHz['theta'] / np.pi, ourprot_100MHz['keyrate'], color='tab:blue', linewidth=3, linestyle='-',
         label=fr'This work, $\epsilon$ = 0.01 $\frac{{1}}{{\sigma_t}}$')
ax.plot(ourprot_500MHz['theta'] / np.pi, ourprot_500MHz['keyrate'], color='#4da4e0', linewidth=3, linestyle='--',
         label=fr'This work, $\epsilon$ = 0.05 $\frac{{1}}{{\sigma_t}}$')
ax.plot(ourprot_1100MHz['theta'] / np.pi, ourprot_1100MHz['keyrate'], color='#89c3eb', linewidth=3, linestyle=':',
        label=r'This work, $\epsilon$ = 0.11 $\frac{{1}}{{\sigma_t}}$')
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
    ll = ax.legend(ncol=2, bbox_to_anchor=(0.5, 1.38), loc='upper center', fontsize=16, handlelength=1.5,
                   labelspacing=0, handleheight=2.1)
    # ll = ax.legend(ncol=2, bbox_to_anchor=(0.5, 1.25), loc='upper center', fontsize=20, handlelength=1.5, labelspacing=0, handleheight=2.1)
    for text in ll.texts:
        # text.set_va('bottom')
        text_x, text_y = text.get_position()
        text.set_position((text_x, text_y + 3))
ax.set_xlim([0, 1])
ax.set_ylim([0, 0.53])
plt.show()