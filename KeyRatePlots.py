import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import math
import numpy as np
import os.path
import matplotlib.tri as tri
import pyexcel
from matplotlib import cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

Uheatmap = False
losscomp = True
lossphasechan = False

fsize=22 # fontsize for the figures
matplotlib.rcParams.update({'font.size': fsize})
file = pd.ExcelFile('AllKeyRates_SimulatedData.xlsx')
cmap = matplotlib.colormaps['inferno']
#plt.rcParams.update({"text.usetex": True})

if Uheatmap:
    num_pts = 201

    BB84_u = pd.read_excel(file, sheet_name='BB84_U_coarse_smallgap')
    #BB84_u = pd.read_excel(file2, sheet_name='BB84_heatmap_26pts_shaped')
    #print(BB84_u)


    bb84_u = np.zeros((num_pts,num_pts))
    for i in range(num_pts):
        bb84_u[i,:] = BB84_u[f'column {i+1}']

    print(bb84_u)


    xvec = [0.10, 0.10, 0.8]

    fsize = 30  # fontsize for the figures
    matplotlib.rcParams.update({'font.size': fsize})
    #    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial"]



    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    bbox = ax.get_position()
    new_bbox = (bbox.x0+xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    im1 = ax.imshow(bb84_u, cmap='inferno',
                   extent=[0, 2, 0, 2],
                   interpolation='nearest', origin='lower')
    ax.set_xlabel(r'Rot. angle, $\theta$ ($\pi$ rad)')
    ax.set_ylabel(r'Deph. angle, $\phi$ ($\pi$ rad)')
    #ax.set_title(fr'BB84', pad=10)
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
    cbar.set_ticks([0, 0.25, 0.5])
    cbar.set_ticklabels(['0.00', '0.25', '0.50'])
    plt.show()

    Wang_u = pd.read_excel(file, sheet_name='Wang_U_coarse')
    #Wang_u = pd.read_excel(file2, sheet_name='Wang_U_updated_126pts_0.5pi')
    #Wang_u = pd.read_excel(file2, sheet_name='Wang_heatmap_26pts_shaped')
    wang_u = np.zeros((num_pts,num_pts))
    for j in range(num_pts):
        print(j)
        wang_u[j,:] = Wang_u[f'column {j+1}']

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    bbox = ax.get_position()
    new_bbox = (bbox.x0+xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    im1 = ax.imshow(wang_u, cmap='inferno',
                   extent=[0, 2, 0, 2],
                   interpolation='nearest', origin='lower')
    ax.set_xlabel(r'Rot. angle, $\theta$ ($\pi$ rad)')
    ax.set_ylabel(r'Deph. angle, $\phi$ ($\pi$ rad)')
    #ax.set_title(fr'BB84', pad=10)
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
    cbar.set_ticks([0, 0.25, 0.5])
    cbar.set_ticklabels(['0.00', '0.25', '0.50'])
    plt.show()

    num_pts = 126

    BB84_u_fine = pd.read_excel(file, sheet_name='BB84_U_fine_smallgap')
    #BB84_u = pd.read_excel(file2, sheet_name='BB84_heatmap_26pts_shaped')
    #print(BB84_u)


    bb84_u_fine = np.zeros((num_pts,num_pts))
    for i in range(num_pts):
        bb84_u_fine[i,:] = BB84_u_fine[f'column {i+1}']


    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    bbox = ax.get_position()
    new_bbox = (bbox.x0+xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    im1 = ax.imshow(bb84_u_fine, cmap='inferno',
                   extent=[0, 0.5, 0, 0.5],
                   interpolation='nearest', origin='lower')
    ax.set_xlabel(r'Rot. angle, $\theta$ ($\pi$ rad)')
    ax.set_ylabel(r'Deph. angle, $\phi$ ($\pi$ rad)')
    #ax.set_title(fr'BB84', pad=10)
    # ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    # ax.set_xticklabels(['0.00', '0.25', '0.50', '0.75', '1.00'])
    ax.set_xticks([0, 0.25, 0.5])
    ax.set_xticklabels(['0.00', '0.25', '0.50'])
    ax.set_yticks([0, 0.25, 0.5])
    ax.set_yticklabels(['0.00', '0.25', '0.50'])
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
    cbar.set_ticks([0, 0.25, 0.5])
    cbar.set_ticklabels(['0.00', '0.25', '0.50'])
    plt.show()



    Wang_u_fine = pd.read_excel(file, sheet_name='Wang_U_fine')
    #BB84_u = pd.read_excel(file2, sheet_name='BB84_heatmap_26pts_shaped')
    #print(BB84_u)


    wang_u_fine = np.zeros((num_pts,num_pts))
    for i in range(num_pts):
        wang_u_fine[i,:] = Wang_u_fine[f'column {i+1}']


    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    bbox = ax.get_position()
    new_bbox = (bbox.x0+xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    im1 = ax.imshow(wang_u_fine, cmap='inferno',
                   extent=[0, 0.5, 0, 0.5],
                   interpolation='nearest', origin='lower')
    ax.set_xlabel(r'Rot. angle, $\theta$ ($\pi$ rad)')
    ax.set_ylabel(r'Deph. angle, $\phi$ ($\pi$ rad)')
    #ax.set_title(fr'BB84', pad=10)
    # ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    # ax.set_xticklabels(['0.00', '0.25', '0.50', '0.75', '1.00'])
    ax.set_xticks([0, 0.25, 0.5])
    ax.set_xticklabels(['0.00', '0.25', '0.50'])
    ax.set_yticks([0, 0.25, 0.5])
    ax.set_yticklabels(['0.00', '0.25', '0.50'])
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
    cbar.set_ticks([0, 0.25, 0.5])
    cbar.set_ticklabels(['0.00', '0.25', '0.50'])
    plt.show()



if losscomp:
    bb84 = pd.read_excel(file, sheet_name='BB84_loss')
    print(bb84)

    # bb84_tenth = pd.read_excel(file, sheet_name='BB84_sqrt(0.1)_loss')

    sr24 = pd.read_excel(file, sheet_name='OurProt_loss')
    print(sr24)

    boileau4 = pd.read_excel(file, sheet_name='Boileau4_loss')
    # print(boileau4)

    boileau3 = pd.read_excel(file, sheet_name='Boileau3_loss')
    # print(boileau3)

    wang = pd.read_excel(file, sheet_name='Wang_loss')
    # print(wang)

    li = pd.read_excel(file, sheet_name='LiRot_loss')
    # print(li)

    tf = pd.read_excel('C:/Users/ir22317/OneDrive - University of Bristol/PhD/Noiseless QKD Project/Security Analysis/TFQKD_KeyRates.xlsx', sheet_name='TF_QKD_muopt_0.1146')


    colors = cmap([0, 0.45, 0.5, 0.67, 0.87])

    #colors= ['k', 'k', 'k', 'tab:red', 'tab:red', 'tab:blue', 'tab:blue']
    #linesty= ['-', '-', '--', '-', '--', '--', '-']
    #linesty = ['-', '-', '-', '-', '-', '-', '-']
    linesty = ['-', '-', '--', '-', '--', '--', '-']

    fsize = 18  # fontsize for the figures
    matplotlib.rcParams.update({'font.size': fsize})

    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial"]

    xvec = [0.1, -0.005, 0.88]

    boileau4_loss_kr = 0.5*10**(-0.4*boileau4['loss'])
    li_loss_kr = 10 ** (-0.4 * li['loss'])
    bb84_loss_kr = 0.5*10 ** (-0.1 * bb84['loss'])

    print('li @96 dB = ', li_loss_kr[95])
    print('boileau @96 dB = ', boileau4_loss_kr[95])

    colors = cmap([0.85, 0.6, 0.57, 0.0])

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    bbox = ax.get_position()
    new_bbox = (bbox.x0+xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    ax.plot(bb84['loss'], bb84['keyrate'], linestyle='-', color='tab:orange',linewidth=2, label='BB84, theory')
    ax.plot(bb84['loss'], bb84['keyrate']/10**3, linestyle=':', color='tab:orange', linewidth=2, label='BB84, expt.')
    ax.plot(tf['loss'], tf['keyrate'], linestyle='-', color='#ffba7c', linewidth=2, label='TF-QKD, theory')
    ax.plot(tf['loss'], tf['keyrate']/10**2, linestyle=':', color='#ffba7c', linewidth=2, label='TF-QKD, expt.')
    #ax.plot(boileau4['loss'][95],  5.88268440666069e-05, 'x')
    #plt.plot(sr24_XZ['loss'], sr24_XZ['keyrate'], linestyle=':', color='k', label='Experimental implementation')
    ax.plot(boileau3['loss'], boileau3['keyrate'], linestyle='-', color='#863f00', linewidth=2, label='Boileau 3-photon')
    ax.plot(boileau4['loss'], boileau4['keyrate'], linestyle='-.', color='#863f00', linewidth=2, label='Boileau 4-photon')
    #ax.plot(boileau4['loss'], boileau4_loss_kr, linestyle='-', color=colors[4], linewidth=2, label='Boileau 4-photon')
    #ax.plot(li['loss'], li['keyrate'], linestyle=linesty[5], color=colors[5], linewidth=2, label='Li Rotation-only')
    #ax.plot(li['loss'], li_loss_kr, linestyle='-', color=colors[5], linewidth=2, label='Li Rotation-only')
    ax.plot(wang['loss'], wang['keyrate'], linestyle='-', color='tab:blue',linewidth=2, label='Wang')
    ax.plot(sr24['loss'], sr24['keyrate'], linestyle='-', color='tab:blue', linewidth=3, label='This work')
    #plt.plot(cv['loss'], cv['keyrate'], linestyle='-', color='b', label='CV-QKD protocol')
    ax.set_xlabel('Total loss (dB)')
    ax.set_ylabel('Key rate (bits/pulse)')
    ax.set_yscale('log')
    #plt.title('Title')
    #ax.legend(loc='lower left', fontsize=14)
    ax.legend(ncol=2, bbox_to_anchor=(0.5, 1.32), loc='upper center', fontsize=16, handlelength=1.3, labelspacing=0.2)
    ax.set_aspect(1.0)
    ax.set_xlim([0, 8])
    plt.show()

    fsize = 18 # fontsize for the figures
    matplotlib.rcParams.update({'font.size': fsize})

    xvec = [0.1, 0.01, 0.9]

    bb84 = pd.read_excel(file, sheet_name='BB84_wide_loss')

    wang = pd.read_excel(file, sheet_name='Wang_wide_loss')

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    bbox = ax.get_position()
    new_bbox = (bbox.x0+xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    ax.plot(bb84['loss'], bb84['keyrate'], linestyle='-', color='#ffba7c',linewidth=2, label='BB84, theory')
    ax.plot(wang['loss'], wang['keyrate'], linestyle='-', color='#ff9537',linewidth=2, label='Wang')
    #ax.plot(boileau3['loss'], boileau3['keyrate'], linestyle='-', color='#863f00', linewidth=2, label='Boileau 3-photon')
    ax.plot(boileau4['loss'], boileau4['keyrate'], linestyle='-', color='#d16200', linewidth=2, label='Boileau 4-photon')
    ax.plot(li['loss'], li['keyrate'], linestyle='-', color='#9e4a00', linewidth=2, label='Li Rotation-only')
    ax.set_xlabel('Total loss (dB)')
    ax.set_ylabel('Key rate (bits/pulse)')
    ax.set_yscale('log')
    #plt.title('Title')
    #ax.legend(loc='lower left', fontsize=14)
    ax.legend(bbox_to_anchor=(1.0, 1.0), loc='upper right', fontsize=16, handlelength=1.3, labelspacing=0.2)
    ax.set_aspect(5.0)
    ax.set_xlim([0, 65])
    plt.show()









if lossphasechan:
    #file = pd.ExcelFile('C:/Users/ir22317/OneDrive - University of Bristol/PhD/Noiseless QKD Project/Security Analysis/Asymm_Losses_Interpolate.xlsx')
    ###INTERPOLATION FUCNTION###
    # Asymmetric losses keyrate as a function of asymmetric loss (for zero loss channel)
    Asym_loss = pd.read_excel(file, sheet_name='OurProt_AmpDamp_fine_smallgap')
    params = pd.read_excel(file, sheet_name='OurProt_AmpDamp_fine_params')
    eta0 = params['eta0']
    phase = params['delta']
    # we want to interpolate the data with a given key rate to get a given eta0

    Eta, Phase = np.meshgrid(eta0, phase)

    num_pts=126

    asymm_loss = np.zeros((num_pts, num_pts))
    for i in range(num_pts):
        asymm_loss[i, :] = Asym_loss[f'column {i + 1}']


    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial"]


    xvec = [0.10, 0.10, 0.8]

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    bbox = ax.get_position()
    new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    im1 = ax.imshow(asymm_loss, cmap='inferno',
                    extent=[0, 1, 0, 0.5],
                    interpolation='nearest', origin='lower',
                    aspect=2)
    ax.set_xlabel(r'State $|1\rangle$ trans. prob., $\eta_0$')
    ax.set_ylabel(r'Decoh. angle $\Delta$/$\pi$ rad')
    #ax.set_title(fr'$\Omega_1$ = {Omega1[i]/(2*np.pi)} THz', pad=10)
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(['0.0', '0.5', '1.0'])
    # ax.set_xticks([0, 1, 2])
    # ax.set_xticklabels(['0.0', '1.0', '2.0'])
    ax.set_yticks([0, 0.25, 0.5])
    ax.set_yticklabels(['0.00', '0.25', '0.50'])
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
    cbar.ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
    cbar.ax.set_yticklabels(['0.0', '0.1', '0.2', '0.3', '0.4', '0.5'])
    plt.show()

    transparency = [1, 0.85, 0.7, 0.55]
    colors = cmap([0.1, 0.3, 0.6, 0.9])

    xvec = [0.0, 0.05, 1]

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(projection='3d')
    bbox = ax.get_position()
    new_bbox = (bbox.x0+xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    ax.plot_surface(Eta, Phase / np.pi, asymm_loss, color=colors[2], alpha=0.9, rstride=10,cstride=10, label=fr'Our protocol')
    ax.set_xlabel(r"Asymm. trans."
           "\n"
           r"prob., $\eta_0$", labelpad=25)
    ax.set_ylabel('Decoh. angle, \n $\Delta$ ($\pi$ rad)', labelpad=25)
    ax.set_zlabel('Key rate\n (bits/pulse)', labelpad=25)
    # plt.legend(bbox_to_anchor=(1, 1.11), loc='upper right', fontsize=18)
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(['0.0', '0.5', '1.0'])
    ax.set_yticks([0, 0.25, 0.5])
    ax.set_yticklabels(['0.0', '0.25', '0.50'])
    ax.set_zticks([0, 0.25, 0.5])
    ax.set_zticklabels(['0.00', '0.25', '0.50'])
    ax.zaxis.label.set_rotation(90)
    ax.zaxis.set_rotate_label(False)
    ax.set_xlim([0,1])
    ax.set_ylim([0, 0.5])
    ax.tick_params(axis='z', pad=9)
    ax.tick_params(axis='x', pad=7)
    ax.tick_params(axis='y', pad=7)
    ax.view_init(20, -150)
    plt.show()

    xvec = [0.05, 0.05, 0.95]

    colors = cmap([0.8, 0.6, 0.4, 0.0])

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    bbox = ax.get_position()
    new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
    ax.set_position(new_bbox)
    ax.plot(eta0, asymm_loss[0, :], color=cmap(0), linewidth=3)#, label='BB84')
    ax.set_xlabel(r'Asymmetric transmission probability, $\eta_0$',labelpad=7)
    ax.set_ylabel('Key rate (bits/pulse)',labelpad=7)
    #ax.set_title(r'Key rate vs unitary rotation $\theta$, with $\phi$=0', pad=10)
    ax.set_aspect(1.4)
    #ax.legend(bbox_to_anchor=(0.5, 0.05), loc='lower center', fontsize=20, handlelength=1, labelspacing=0.3)
    ax.set_xlim([0, 1])
    ax.set_ylim([0,0.5])
    ax.tick_params(axis='x', pad=5)
    ax.tick_params(axis='y', pad=5)
    plt.show()

