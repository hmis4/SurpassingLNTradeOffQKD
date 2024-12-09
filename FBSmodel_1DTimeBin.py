import numpy as np
import math
import cmath
import matplotlib.pyplot as plt
from scipy.special import erf
from scipy import integrate
import matplotlib as mpl
#the definition of erf(x) used in the python package is 2/sqrt(pi)*integrate(exp(-z^2))z=0 to z=x, x can be a real number

def gauss1D(t, tau, sigma):
    return (np.pi)**(-0.25)*sigma**(-0.5)*np.exp(-(t-tau)**2/(2*sigma**2))

def FTgauss1D(w, tau, sigma):
    return (np.pi)**(-0.25)*sigma**(0.5)*sigma*np.exp(-1j*tau*w)*np.exp(-sigma**2*w**2/2)

def inverseFTshift(t, tau, sigma, shift, a, b):
    return 1/2*(np.pi)**(-0.25)*sigma**(-0.5)*np.exp(-1j*shift*t)*np.exp(-(t-tau)**2/(2*sigma**2))*(erf(sigma/np.sqrt(2)*b+sigma/np.sqrt(2)*shift-sigma/np.sqrt(2)*1j*(t-tau)/sigma**2) - erf(sigma/np.sqrt(2)*a+sigma/np.sqrt(2)*shift-sigma/np.sqrt(2)*1j*(t-tau)/sigma**2))

def swapgauss1D(t,tau,sigma,Omega,mu,epsilon,theta,phi):
    return inverseFTshift(t,tau,sigma,0,-np.inf,Omega-epsilon/2) + np.cos(theta)*inverseFTshift(t,tau,sigma,0,Omega-epsilon/2, Omega+epsilon/2) - np.exp(-1j*phi)*np.sin(theta)*inverseFTshift(t,tau,sigma,mu,Omega-epsilon/2,Omega+epsilon/2) + inverseFTshift(t,tau,sigma,0,Omega+epsilon/2, Omega+mu-epsilon/2) + np.cos(theta)*inverseFTshift(t,tau,sigma,0,Omega+mu-epsilon/2,Omega+mu+epsilon/2) + np.exp(1j*phi)*np.sin(theta)*inverseFTshift(t,tau,sigma,-mu,Omega+mu-epsilon/2,Omega+mu+epsilon/2) + inverseFTshift(t,tau,sigma,0,Omega+mu+epsilon/2, np.inf)
###THIS HAS BEEN EDITED TO HAVE THE CORRECT U CV MODEL

tau=0
tau2=220
sigma=17

Omega = 0
epsilon=0.019
mu=0.0005

theta=np.pi
phi=np.pi/2

t = np.linspace(-400,400,10001)
w = (-20,20, 10001)
# errorfunc=erf(0.5*(-np.inf) + 0.5*1j*t)
# print(errorfunc)
# plt.plot(t,np.real(errorfunc))
# plt.show()
#
# #timebinswapsec3=1/2*np.exp(-1j*mu*t)*np.exp(-(t-tau)**2/(2*sigma**2))*(erf(sigma/np.sqrt(2)*(Omega+epsilon+mu-1j*(t-tau)/sigma**2)) - erf(sigma/np.sqrt(2)*(Omega+mu-1j*(t-tau)/sigma**2)))
# timebinswapsec3=1/2*np.exp(-1j*mu*t)*np.exp(-(t-tau)**2/(2*sigma**2))*(erf(sigma/np.sqrt(2)*(-np.inf)+sigma/np.sqrt(2)*epsilon+sigma/np.sqrt(2)*mu-sigma/np.sqrt(2)*1j*(t-tau)/sigma**2) - erf(sigma/np.sqrt(2)*Omega+sigma/np.sqrt(2)*epsilon+sigma/np.sqrt(2)*mu-sigma/np.sqrt(2)*1j*(t-tau)/sigma**2))
# print(timebinswapsec3)
#
# plt.plot(t,timebinswapsec3)
# plt.show()
#
# timebinswapsec=inverseFTshift(t,tau,sigma,mu,-np.inf,Omega)
#
# plt.plot(t,timebinswapsec)
# plt.show()



timebin_noswap = gauss1D(t,tau,sigma)
timebin2_noswap = gauss1D(t,tau2,sigma)
timebin_swap = swapgauss1D(t,tau,sigma,Omega,mu,epsilon,theta,phi)
timebin2_swap = swapgauss1D(t,tau2,sigma,Omega,mu,epsilon, theta,phi)
freqbin1 = FTgauss1D(t,Omega,epsilon/3)
I1 = integrate.simpson(np.conjugate(timebin_noswap)*timebin_noswap, x=t)
I2 = integrate.simpson(np.conjugate(timebin_noswap)*timebin_swap, x=t)
I3 = integrate.simpson(np.conjugate(timebin_noswap)*freqbin1,x=t)
I4 = integrate.simpson(np.conjugate(timebin_swap)*freqbin1,x=t)

J1 = integrate.simpson(np.conjugate(timebin_noswap)*timebin2_noswap,x=t)
J2 = integrate.simpson(np.conjugate(timebin_swap)*timebin2_swap,x=t)

# timebin_f_noswap = np.zeros(len(w))
# timebin_f_swap =  np.zeros(len(w))
# for i in range(len(timebin_f_swap)):
#     timebin_f_noswap[i] = integrate.simpson(np.conjugate(timebin_f_noswap * np.exp(-1j * w[i] * t), x=t))
#     timebin_f_swap[i] = integrate.simpson(np.conjugate(timebin_swap*np.exp(-1j*w[i]*t), x=t))

print(I1,I2)
print(np.abs(I1)**2, np.abs(I2)**2, np.abs(I3)**2, np.abs(I4)**2)

print(J1,J2)
print(np.abs(J1)**2,np.abs(J2)**2)

cmap = mpl.colormaps['inferno']
transparency = [1, 0.85, 0.7, 0.55]
colors = cmap([0.0, 0.7, 0.5, 0.4])


cmap = mpl.colormaps['inferno']
transparency = [1, 0.85, 0.7, 0.55]
colors = cmap([0.0, 0.7, 0.5, 0.4])

xvec = [0.05, 0.05, 0.95]

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]
fsize=22 # fontsize for the figures
mpl.rcParams.update({'font.size': fsize})


fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
bbox = ax.get_position()
new_bbox = (bbox.x0 + xvec[0], bbox.y0 + xvec[1], bbox.width * xvec[2], bbox.height * xvec[2])
ax.set_position(new_bbox)
ax.plot(t,timebin_noswap, label='Original time-bin state', color=colors[0])
#plt.plot(t,freqbin1,label='freqbin')
#plt.plot(t,np.abs(timebin_swap), label='abs swap 1')
ax.plot(t,np.real(timebin_swap), label='Re(Transformed time-bin state)', color=colors[1])
ax.plot(t,np.imag(timebin_swap), label='Im(Transformed time-bin state)', color=colors[2])
#plt.plot(t,timebin2_noswap, label='noswap 2')
#plt.plot(t,np.abs(timebin2_swap), label='re swap 2')
#plt.plot(t,np.imag(timebin2_swap), label='im swap 2')
ax.set_xlabel('Time (ps)')
ax.set_ylabel('Wavefunction amplitude')
plt.legend(bbox_to_anchor=(0, 1), loc='upper left', fontsize=18)
plt.show()

# plt.plot(w,np.abs(timebin_f_noswap), label='Original time-bin state', color=colors[0])
# #plt.plot(t,freqbin1,label='freqbin')
# #plt.plot(t,np.abs(timebin_swap), label='abs swap 1')
# plt.plot(w,np.real(timebin_f_swap), label='Re(Transformed time-bin state)', color=colors[1])
# plt.plot(w,np.imag(timebin_f_swap), label='Im(Transformed time-bin state)', color=colors[2])
# #plt.plot(t,timebin2_noswap, label='noswap 2')
# #plt.plot(t,np.abs(timebin2_swap), label='re swap 2')
# #plt.plot(t,np.imag(timebin2_swap), label='im swap 2')
# plt.xlabel('Time/ps')
# plt.ylabel('Wavefunction amplitude')
# plt.legend()
# plt.show()

sentstate1=swapgauss1D(t,tau,sigma,Omega,mu,epsilon,np.pi/2,np.pi)
sentstate2=swapgauss1D(t,tau2,sigma,Omega,mu,epsilon,np.pi/2,np.pi)
measurestate1 = swapgauss1D(t,tau,sigma,Omega,mu,epsilon,0,0)
measurestate2 = swapgauss1D(t,tau2,sigma,Omega,mu,epsilon,0,0)

J1 = integrate.simpson(np.conjugate(sentstate1)*measurestate1,x=t)
J2 = integrate.simpson(np.conjugate(sentstate2)*measurestate2,x=t)
J3 = integrate.simpson(np.conjugate(sentstate2)*measurestate1,x=t)
J4 = integrate.simpson(np.conjugate(sentstate1)*measurestate2,x=t)
norm11 = integrate.simpson(np.conjugate(sentstate1)*sentstate1,x=t)
norm12 = integrate.simpson(np.conjugate(sentstate2)*sentstate2,x=t)
norm13 = integrate.simpson(np.conjugate(sentstate1)*sentstate2,x=t)
norm14 = integrate.simpson(np.conjugate(sentstate2)*sentstate1,x=t)
norm21 = integrate.simpson(np.conjugate(measurestate1)*measurestate1,x=t)
norm22 = integrate.simpson(np.conjugate(measurestate2)*measurestate2,x=t)
norm23 = integrate.simpson(np.conjugate(measurestate1)*measurestate2,x=t)
norm24 = integrate.simpson(np.conjugate(measurestate2)*measurestate1,x=t)

print(J1,J2,J3,J4)
print(norm11,norm12, norm13,norm14)
print(norm21,norm22, norm23, norm24)
norm1 = 1/np.sqrt(2*np.abs(norm11*norm12-norm13*norm14))
norm2 = 1/np.sqrt(2*np.abs(norm21*norm22-norm23*norm24))
print(np.abs(norm11)**2,np.abs(norm12)**2)
overlap=2*(J1*J2-J3*J4)*norm1*norm2
prob=np.abs(overlap)**2
print(prob)



# thetavec = np.linspace(0,np.pi,100)
# phivec = np.linspace(0,np.pi,100)
# overlap_prob = np.zeros((len(thetavec),len(phivec)))
#
# for i in range(len(thetavec)):
#     for j in range(len(phivec)):
#         bin1 = swapgauss1D(t, tau, sigma, Omega, mu, epsilon, theta, phi)
#         bin2 = swapgauss1D(t, tau2, sigma, Omega, mu, epsilon, theta, phi)
#         overlap = integrate.simpson(np.conjugate(bin1)*bin2,x=t)
#         overlap_prob[i][j] = np.abs(overlap)**2
#
# c = plt.imshow(np.abs(overlap_prob), cmap='viridis',
#                extent=[thetavec[0], thetavec[-1], phivec[0], phivec[-1]],
#                interpolation='nearest', origin='lower', aspect=(thetavec[0]-thetavec[-1])/(phivec[0]-phivec[-1]))
# plt.colorbar(c)
# plt.title('Abs value')
# plt.xlabel('Frequency bin 1/2pi*THz')
# plt.ylabel('Frequency bin width/2pi*THz')
# plt.show()
