import numpy as np
import math
import cmath
import matplotlib.pyplot as plt
from scipy.special import erf
import matplotlib as mpl
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

def Uarboverlap1D(tau1,tau2, sigma_t, Omega1, epsilon1, mu1, Omega2, epsilon2, mu2, Omega3, mu3, epsilon3, theta1, phi1, theta2, phi2, theta3, phi3):
    a = Omega1 - epsilon1/2
    b = Omega1 + epsilon1 / 2
    c = Omega1 - epsilon1/2 + mu1
    d = Omega1 + epsilon1 / 2 + mu1
    e = Omega2 - epsilon2/2
    f = Omega2 + epsilon2 / 2
    g = Omega2 - epsilon2/2 + mu2
    h = Omega2 + epsilon2 / 2 + mu2
    i = Omega3 - epsilon3/2
    j = Omega3 + epsilon3 / 2
    k = Omega3 - epsilon3/2 + mu3
    l = Omega3 + epsilon3 / 2 + mu3
    startings = np.asarray([a,b,c,d,e,f,g,h,i,j,k,l])
    angles = np.asarray([-theta1, -phi1, theta1, phi1, -theta2, -phi2, theta2, phi2, -theta3, -phi3, theta3, phi3])
    shifts = np.asarray([mu1, 0, -mu1, 0, mu2, 0, -mu2, 0, mu3, 0, -mu3, 0])
    # num_in_pair = np.asarray([0,0,1,1,0,0,1,1,0,0,1,1])
    # pairings = np.asarray([0,0,0,0,1,1,1,1,2,2,2,2])

    #print(startings)

    #ordered_startings = np.sort(startings)
    ord_ind = np.argsort(startings)
    ord_start = startings[ord_ind]
    ord_angles = angles[ord_ind]
    ord_shifts = shifts[ord_ind]
    # ord_pairs = pairings[ord_ind]
    # ord_num = num_in_pair[ord_ind]

    # print(ord_ind)
    # print(ord_start)
    # print(ord_angles)
    # print(ord_shifts)

    # print(ord_pairs)
    # print(ord_num)
    # print(startings[np.where(startings==ord_start[0])[0]+2])
    # np.where(ord_ind==0)[0]


    f1 = overlapshift(tau1, tau2, sigma_t, 0, -math.inf, ord_start[0])
    f2 = np.cos(ord_angles[0])*overlapshift(tau1, tau2, sigma_t, 0, ord_start[0], ord_start[1]) + np.exp(1j*ord_angles[1])*np.sin(ord_angles[0])*overlapshift(tau1, tau2, sigma_t, ord_shifts[0], ord_start[0],ord_start[1])
    f3 = overlapshift(tau1, tau2, sigma_t, 0, ord_start[1], ord_start[2])
    f4 = np.cos(ord_angles[2])*overlapshift(tau1, tau2, sigma_t, 0, ord_start[2], ord_start[3]) + np.exp(1j*ord_angles[3])*np.sin(ord_angles[2])*overlapshift(tau1, tau2, sigma_t, ord_shifts[2], ord_start[2], ord_start[3])
    f5 = overlapshift(tau1, tau2, sigma_t, 0, ord_start[3], ord_start[4])
    f6 = np.cos(ord_angles[4]) * overlapshift(tau1, tau2, sigma_t, 0, ord_start[4], ord_start[5]) + np.exp(1j * ord_angles[5]) * np.sin(ord_angles[4]) * overlapshift(tau1, tau2, sigma_t, ord_shifts[4], ord_start[4], ord_start[5])
    f7 = overlapshift(tau1, tau2, sigma_t, 0, ord_start[5], ord_start[6])
    f8 = np.cos(ord_angles[6])*overlapshift(tau1, tau2, sigma_t, 0, ord_start[6], ord_start[7]) + np.exp(1j*ord_angles[7])*np.sin(ord_angles[6])*overlapshift(tau1, tau2, sigma_t, ord_shifts[6], ord_start[6], ord_start[7])
    f9 = overlapshift(tau1, tau2, sigma_t, 0, ord_start[7], ord_start[8])
    f10 = np.cos(ord_angles[8])*overlapshift(tau1, tau2, sigma_t, 0, ord_start[8], ord_start[9]) + np.exp(1j*ord_angles[9])*np.sin(ord_angles[8])*overlapshift(tau1, tau2, sigma_t, ord_shifts[8], ord_start[8], ord_start[9])
    f11 = overlapshift(tau1, tau2, sigma_t, 0, ord_start[9], ord_start[10])
    f12 = np.cos(ord_angles[10])*overlapshift(tau1, tau2, sigma_t, 0, ord_start[10], ord_start[11]) + np.exp(1j*ord_angles[11])*np.sin(ord_angles[10])*overlapshift(tau1, tau2, sigma_t, ord_shifts[10], ord_start[10], ord_start[11])
    f13 = overlapshift(tau1, tau2, sigma_t, 0, ord_start[11], math.inf)
    return f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8 + f9 + f10 + f11 + f12 + f13


def Uarboverlap2D(tau1,tau2, sigma_t, Omega1, epsilon1, mu1, Omega2, epsilon2, mu2, Omega3, mu3, epsilon3, theta1, phi1, theta2, phi2, theta3, phi3):
    return 2*A(tau1, tau2, sigma_t)**2*(Uarboverlap1D(tau1, tau1, sigma_t, Omega1, epsilon1, mu1, Omega2, epsilon2, mu2, Omega3, mu3, epsilon3, theta1, phi1, theta2, phi2, theta3, phi3)*Uarboverlap1D(tau2, tau2, sigma_t, Omega1, epsilon1, mu1, Omega2, epsilon2, mu2, Omega3, mu3, epsilon3, theta1, phi1, theta2, phi2, theta3, phi3)-Uarboverlap1D(tau1, tau2, sigma_t, Omega1, epsilon1, mu1, Omega2, epsilon2, mu2, Omega3, mu3, epsilon3, theta1, phi1, theta2, phi2, theta3, phi3)*Uarboverlap1D(tau2, tau1, sigma_t, Omega1, epsilon1, mu1, Omega2, epsilon2, mu2, Omega3, mu3, epsilon3, theta1, phi1, theta2, phi2, theta3, phi3))


def Uarbprob2D(tau1, tau2, sigma_t, Omega1, epsilon1, mu1, Omega2, epsilon2, mu2, Omega3, mu3, epsilon3, theta1, phi1, theta2, phi2, theta3, phi3):
    return np.abs(Uarboverlap2D(tau1, tau2, sigma_t, Omega1, epsilon1, mu1, Omega2, epsilon2, mu2, Omega3, mu3, epsilon3, theta1, phi1, theta2, phi2, theta3, phi3))**2



fsize=20 # fontsize for the figures
mpl.rcParams.update({'font.size': fsize})

# tau1 = 0
# tau2 = 6
# sigma_t = 1
tau1 = 0
tau2 = 220 # ps
sigma_t = 17 #ps

# Omega = [0, 0.25, 0.5]
# epsilon = [0.1, 0.1, 0.1]
# mu = [0.1, 0.1, 0.1]

# Omega = 2*np.pi*np.asarray([0, -0.026, 0.007])
# epsilon = 2*np.pi*6*np.asarray([0.0011, 0.0011, 0.0011])
# mu = 2*np.pi*np.asarray([0.019, 0.019, 0.019])

Omega = 2*np.pi*np.asarray([0, -0.026, 0.007])
epsilon = 2*np.pi*6*np.asarray([0.0005, 0.0005, 0.0005])
mu = 2*np.pi*np.asarray([0.019, 0.019, 0.019])

theta_vec = [0, np.pi/2, np.pi/2]
phi_vec = [0, 0, 0]
# theta_vec = [0, 2,3]
# phi_vec=[0, 2.5, 3.5]

# tau1 = 0
# tau2 = 1000
# sigma_t = 30
#
# Omega = np.asarray([0, 0.012, 0.024])*2*np.pi
# epsilon = np.asarray([0.0013, 0.0013, 0.0013])*6*2*np.pi
# mu = np.asarray([0.01, 0.01, 0.01])*2*np.pi
#
# theta_vec = [0, np.pi/2, np.pi/2]
# phi_vec = [0, 0, 0]


theta = np.linspace(0, np.pi, 101)
phi = np.linspace(0, np.pi, 101)

Theta, Phi = np.meshgrid(theta, phi)
probs1 = np.zeros((len(theta),len(phi)))
probs2 = np.zeros((len(theta),len(phi)))
probs3 = np.zeros((len(theta),len(phi)))

# #TESTING TESTING
# probs1 = Uarbprob2D(tau1, tau2, sigma_t, Omega[0], epsilon[0], mu[0], Omega[1], epsilon[1], mu[1], Omega[2], mu[2], epsilon[2], Theta, Phi, 0, 0, 0, 0)
# probs2 = Uarbprob2D(tau1, tau2, sigma_t, Omega[0], epsilon[0], mu[0], Omega[1], epsilon[1], mu[1], Omega[2], mu[2], epsilon[2], Theta, Phi, theta_vec[1], phi_vec[1],0,0)
# probs3 = Uarbprob2D(tau1, tau2, sigma_t, Omega[0], epsilon[0], mu[0], Omega[1], epsilon[1], mu[1], Omega[2], mu[2], epsilon[2], Theta, Phi, theta_vec[1], phi_vec[1], theta_vec[2], phi_vec[2])

for i in range(len(theta)):
    for j in range(len(phi)):
        print(i,j)
        probs1[j,i] = Uarbprob2D(tau1, tau2, sigma_t, Omega[0], epsilon[0], mu[0], Omega[1], epsilon[1], mu[1], Omega[2], mu[2], epsilon[2], theta[i], phi[j], 0, 0, 0, 0)
        probs2[j,i] = Uarbprob2D(tau1, tau2, sigma_t, Omega[0], epsilon[0], mu[0], Omega[1], epsilon[1], mu[1], Omega[2], mu[2], epsilon[2], theta[i], phi[j], theta_vec[1], phi_vec[1],0,0)
        probs3[j,i] = Uarbprob2D(tau1, tau2, sigma_t, Omega[0], epsilon[0], mu[0], Omega[1], epsilon[1], mu[1], Omega[2], mu[2], epsilon[2], theta[i], phi[j],theta_vec[1], phi_vec[1], theta_vec[2], phi_vec[2])


# probs1 = Uarbprob2D(tau1, tau2, sigma_t, Omega[0], epsilon[0], mu[0], Omega[1], epsilon[1], mu[1], Omega[2], mu[2], epsilon[2], Theta, Phi, 0, 0, 0, 0)
# probs2 = Uarbprob2D(tau1, tau2, sigma_t, Omega[0], epsilon[0], mu[0], Omega[1], epsilon[1], mu[1], Omega[2], mu[2], epsilon[2], Theta, Phi, theta_vec[1], phi_vec[1],0,0)
# probs3 = Uarbprob2D(tau1, tau2, sigma_t, Omega[0], epsilon[0], mu[0], Omega[1], epsilon[1], mu[1], Omega[2], mu[2], epsilon[2], Theta, Phi, theta_vec[1], phi_vec[1], theta_vec[2], phi_vec[2])
#



# c = plt.imshow(probs3, cmap='viridis',
#                extent=[theta[0], theta[-1], phi[0], phi[-1]],
#                interpolation='nearest', origin='lower',
#                aspect=(theta[0] - theta[-1]) / (phi[0] - phi[-1]))
# plt.colorbar(c)
# # plt.title('Overlap with shifting')
# plt.xlabel('theta')
# plt.ylabel('phi')
# plt.show()

cmap = mpl.colormaps['inferno']
colors = cmap([0.9, 0.7, 0.4])
transparency = [0.7,0.85,1]

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
ax.plot_surface(Theta/np.pi, Phi/np.pi, probs1, rstride=10, cstride=10, label='1 pair', alpha=transparency[0],color=colors[0])
ax.plot_surface(Theta/np.pi, Phi/np.pi, probs2, rstride=10, cstride=10, label='2 pairs', alpha=transparency[1], color=colors[1])
ax.plot_surface(Theta/np.pi, Phi/np.pi, probs3, rstride=10, cstride=10, label='3 pairs', alpha = transparency[2], color=colors[2])
ax.set_xlabel(r"Rot. angle, "
              "\n"
              r"$\theta$ ($\pi$ rad)", labelpad=25)
ax.set_ylabel('Deph. angle, \n$\phi$ ($\pi$ rad)', labelpad=25)
#ax.set_zlabel(r'$\left|\mathrm{1}\right\rangle$ meas. prob.', labelpad=10)
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


# c = plt.imshow(probs1, cmap='inferno',
#                extent=[theta[0]/np.pi, theta[-1]/np.pi, phi[0]/np.pi, phi[-1]/np.pi],
#                interpolation='nearest', origin='lower',
#                aspect=(theta[0] - theta[-1]) / (phi[0] - phi[-1]))
# plt.colorbar(c)
# # plt.title('Overlap with shifting')
# plt.xlabel(r'Rotation angle $\theta$/$\pi$')
# plt.ylabel(r'Dephasing angle $\phi$/$\pi$')
# plt.title('One swap')
# plt.show()
#
# c = plt.imshow(probs2, cmap='inferno',
#                extent=[theta[0]/np.pi, theta[-1]/np.pi, phi[0]/np.pi, phi[-1]/np.pi],
#                interpolation='nearest', origin='lower',
#                aspect=(theta[0] - theta[-1]) / (phi[0] - phi[-1]))
# plt.colorbar(c)
# # plt.title('Overlap with shifting')
# plt.xlabel(r'Rotation angle $\theta$/$\pi$')
# plt.ylabel(r'Dephasing angle $\phi$/$\pi$')
# plt.title('Two swaps')
# plt.show()
#
# c = plt.imshow(probs3, cmap='inferno',
#                extent=[theta[0]/np.pi, theta[-1]/np.pi, phi[0]/np.pi, phi[-1]/np.pi],
#                interpolation='nearest', origin='lower',
#                aspect=(theta[0] - theta[-1]) / (phi[0] - phi[-1]))
# plt.colorbar(c)
# # plt.title('Overlap with shifting')
# plt.xlabel(r'Rotation angle $\theta$/$\pi$')
# plt.ylabel(r'Dephasing angle $\phi$/$\pi$')
# plt.title('Three swaps')
# plt.show()


xvec = [0.10, 0.10, 0.8]

fsize = 30  # fontsize for the figures
mpl.rcParams.update({'font.size': fsize})

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
bbox = ax.get_position()
new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
ax.set_position(new_bbox)
c = ax.imshow(probs1, cmap='inferno', vmin=0, vmax=1,
               extent=[theta[0]/np.pi, theta[-1]/np.pi, phi[0]/np.pi, phi[-1]/np.pi],
               interpolation='nearest', origin='lower',
               aspect=(theta[0] - theta[-1]) / (phi[0] - phi[-1]))
ax.set_xlabel(r'Rot. angle $\theta$/$\pi$ rad')
ax.set_ylabel(r'Deph. angle $\phi$/$\pi$ rad')
ax.set_title('1 $U_f$ operation', pad=10)
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
cbar = fig.colorbar(c, cax=axins)
cbar.set_label(label=r'$\left|\mathrm{1}\right\rangle$ measurement prob.', labelpad=10)
plt.show()


fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
bbox = ax.get_position()
new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
ax.set_position(new_bbox)
c = ax.imshow(probs2, cmap='inferno', vmin=0, vmax=1,
               extent=[theta[0]/np.pi, theta[-1]/np.pi, phi[0]/np.pi, phi[-1]/np.pi],
               interpolation='nearest', origin='lower',
               aspect=(theta[0] - theta[-1]) / (phi[0] - phi[-1]))
ax.set_xlabel(r'Rot. angle $\theta$/$\pi$ rad')
ax.set_ylabel(r'Deph. angle $\phi$/$\pi$ rad')
ax.set_title('2 $U_f$ operations', pad=10)
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
cbar = fig.colorbar(c, cax=axins)
cbar.set_label(label=r'$\left|\mathrm{1}\right\rangle$ measurement prob.', labelpad=10)
plt.show()


fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
bbox = ax.get_position()
new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
ax.set_position(new_bbox)
c = plt.imshow(probs3, cmap='inferno',vmin=0, vmax=1,
               extent=[theta[0]/np.pi, theta[-1]/np.pi, phi[0]/np.pi, phi[-1]/np.pi],
               interpolation='nearest', origin='lower',
               aspect=(theta[0] - theta[-1]) / (phi[0] - phi[-1]))
ax.set_xlabel(r'Rot. angle $\theta$/$\pi$ rad')
ax.set_ylabel(r'Deph. angle $\phi$/$\pi$ rad')
ax.set_title('3 $U_f$ operations', pad=10)
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
cbar = fig.colorbar(c, cax=axins)
cbar.set_label(label=r'$\left|\mathrm{1}\right\rangle$ measurement prob.', labelpad=10)
plt.show()

