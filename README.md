# SurpassingLNTradeOffQKD
This repository provides code to reproduce the figures and simulations presented in the paper Surpassing the loss-noise robustness trade-off in quantum key distribution by Hannah Seabrook, Emilien Lavie, Teodor Stroemberg, Matthew P Stafford and Giulia Rubino. This README file provides detailed instructions on using the code to reproduce the data and figures in the paper.


-----
## Figures
-----

### BellStates.py 
Produces figures of the antisymmetric frequency-bin and time-bin Bell-states used for the encoding, i.e. Fig 2 in main text (3D plot) and Fig 1 in the SM (2D spectra)

---

### DispersionModel.py
Produces figures relating to the dispersion model.
- example - makes Fig 12 and Fig 14 in SM (Bell-state fidelities for n=1 and n=2 dispersion)
- compwidth - makes Fig 15 in SM (Bell-state fidelities for varying sigma_omega and sigma_t for n=2 dispersion)
- compOmega2 - makes Fig 16b and 17b in SM (Bell-state fidelities for varying Omega1 and tau1 for n=2 dispersion)
- compOmega1 - makes Fig 16a and 17a in SM (Bell-state fidelities for varying Omega0 and tau0 for n=2 dispersion)
- keyrate1 - makes Fig 13 in SM (key rate for n=1 dispersion with varying sigma_t and Omega1)
- keyrate2 - makes Fig 18 in SM (key rate for n=2 dispersion with varying sigma_t and Omega1)

---

### FBSmodel.py
Produces figures relating to the FBS noise model. 
- If comptau2 - makes Fig 3c in SM (FBS with varying tau_1)
- If comps_t - makes Fig 3d in SM (FBS with varying sigma_t)
- If compOmega - makes Fig 3a in SM (FBS with varying Omega)
- If compeps - makes Fig 4 in the main text (FBS with varying epsilon)
- If compmu - makes Fig 3b in SM (FBS with varying mu)
- If compsigma - makes Fig 8 in Appendix C (FBS with varying sigma_t, for optimising over noise fluctuations)
- If phasefac - makes Fig 4 in SM (the global phase acquired by the time-bin state under the FBS channel)
- If phasefac and krcomp - makes Fig 7 in SM and Fig 5 in the main text (the key rate of our protocol under the FBS channel)
- If savedat - generates and saves data used for the probabilistic FBS channel analysis outlined in Sec III A of the main text

---

### FBSmodel_1DTimeBin.py
Produces illustration of the effect of FBS noise one a 1D time-bin state (Fig 2 in SM).

---

### FBSmodelBitError.py
Produces figures of the bit error fidelity of our protocol (Fig 5 in SM).

---

### KeyRatePlots.py
Produces figures relating to key rate.
- Uheatmap - makes Figs 8 and 9 in SM (key rate of BB84 and Wang under collective unitary noise, with dependence on theta and phi)
- losscomp - makes Fig 6 in the main text and Figs 10 and 11 in SM (key rate as a function of loss; key rate for Li et al (2008) as a function of loss; and key rate for higher losses where instabilities take effect in the numerical simulation)
- lossphasechan - makes Fig 9 in Appendix E (key rate of BB84 under an amplitude damping and decoherence channel)

---

### MultiFBSmodel.py
Produces the figure of the time-bin state fidelity under multiple pairs of FBS operations (Fig 6 in SM).

---

### UnequalDispersionModel.py
Produces the figure on collective noise with different noise parameters between the two photons for the dispersion model with the amplitude damping approximation.
- keyrate1 - makes Fig 20a in SM (key rate under varying deviations in dispersion noise parameters for n=1)
- keyrate2 - makes Fig 20b in SM (key rate under varying deviations in dispersion noise parameters for n=2)

---

### UnequalFBSmodel.py
Produces the figure on collective noise with different noise parameters between the two photons for the FBS model with the amplitude damping approximation.
- difftheta - makes Fig 19d in SM (fidelity for deviating theta between the two photons)
- diffepsilon - makes Fig 19a in SM (fidelity for deviating epsilon between the two photons)
- diffphi - makes Fig 19e in SM (fidelity for deviating phi between the two photons)
- diffmu - makes Fig 19c in SM (fidelity for deviating mu between the two photons)
- diffOmega - makes Fig 19b in SM (fidelity for deviating Omega between the two photons)


---

### FBS_KR_Plots.py
Produces the key rate plots for our protocol under the FBS channel with the two-photon ququart security analysis simulations. (Fig 5)


---

### GVDOverlaps.py
Produces the key rate plots for our protocol under the dispersion channel with the two-photon ququart security analysis simulations, for both linear (Fig S6) and quadratic dispersion (Fig S8), and with unequal dispersion parameters on the two photons (S20).

---

### AllKeyRates_SimulatedData.xlsx
Excel file containing all the data simulated by the security analysis for varying parameters.
- BB84_loss - BB84 key rate as a function of loss from 0 to 10dB
- OurProt_loss - our key rate as a function of loss from 0 to 10dB with amplitude damping model
- Wang_loss - Wang key rate as a function of loss from 0 to 10dB
- Boileau3_loss - Boileau 3-photon protocol key rate as a function of loss from 0 to 10dB
- Boileau4_loss - Boileau 4-photon key rate as a function of loss from 0 to 10dB
- LiRot_loss - Li rotation-only protocol key rate as a function of loss from 0 to 10dB
- LiDeph_loss - Li dephasing-only key rate as a function of loss from 0 to 10dB
- BB84_wide_loss - BB84 key rate as a function of loss from 0 to 70dB
- Wang_wide_loss - BB84 key rate as a function of loss from 0 to 35dB
- BB84_theta - BB84 key rate as a function of unitary parameter theta from 0 to pi, with phi=0
- Wang_theta - Wang key rate as a function of unitary parameter theta from 0 to pi, with phi=0
- Boileau4_theta - Boileau 4-photon protocol key rate as a function of theta from 0 to pi, with phi=0
- LiRot_theta - Li rotation-only key rate as a function of theta from 0 to pi, with phi=0
- LiRot_phi - Li rotation-only key rate as a function of phi from 0 to pi, for different fixed values of theta
- LiDeph_theta - Li dephasing-only key rate as a function of theta from 0 to pi, with phi=0
- OurProt_AmpDamp_fine_smallgap - BB84 under an amplitude damping and decoherence channel, arranged as a heatmap array, corresponds to the action of the FBS on our protocol, with amplitude damping model
- OurProt_AmpDAmp_fine_params - the parameter vectors used for the amplitude damping and decoherence phase, with amplitude damping model
- BB84_U_coarse_smallgap - BB84 under a unitary operation, arranged as a heatmap array, over the range 0 to 2 pi
- BB84_U_coarse_params - the parameter vectors used for the theta and phi values
- BB84_U_fine_smallgap - BB84 under a unitary operation, arranged as a heatmap array, over the range 0 to pi/2
- BB84_U_fine_params - the parameter vectors used for the theta and phi values
- Wang_U_coarse - Wang under a unitary operation, arranged as a heatmap array, over the range 0 to 2 pi
- Wang_U_coarse_params - the parameter vectors used for the theta and phi values
- Wang_U_fine - Wang under a unitary operation, arranged as a heatmap array, over the range 0 to pi/2
- Wang_U_fine_params - the parameter vectors used for the theta and phi values

---

### KeyRates_FBSChannel.xlsx
Excel file containing all the data simulated by the security analysis for varying parameters for our protocol with the FBS model.
- FBS_theta_sigma_1.1GHz - scan with theta of the FBS channel for epsilon = 6sigma_w and sigma_w = 1.1GHz
- FBS_theta_sigma_0.5GHz - scan with theta of the FBS channel for epsilon = 6sigma_w and sigma_w = 0.5GHz
- FBS_theta_sigma_0.1GHz - scan with theta of the FBS channel for epsilon = 6sigma_w and sigma_w = 0.1GHz
- Loss_only - ideal channel with only loss
- LossMixed_1.1GHz - scan with loss for the mixed FBS channel for epsilon = 6sigma_w and sigma_w = 1.1GHz
- LossMixed_0.5GHz - scan with loss for the mixed FBS channel for epsilon = 6sigma_w and sigma_w = 0.5GHz
- LossMixed_0.1GHz - scan with loss for the mixed FBS channel for epsilon = 6sigma_w and sigma_w = 0.1GHz
- FBS_1.1GHz_theta_phi_0_pi - 2D scan of theta and phi between 0 and pi for the FBS channel for epsilon = 6sigma_w and sigma_w = 1.1GHz
- FBS_1.1GHz_theta_phi_0_2pi - 2D scan of theta and phi between 0 and 2pi for the FBS channel for epsilon = 6sigma_w and sigma_w = 1.1GHz
- FBS_0.5GHz_theta_phi_0_2pi - 2D scan of theta and phi between 0 and 2pi for the FBS channel for epsilon = 6sigma_w and sigma_w = 0.5GHz
- FBS_0.1GHz_theta_phi_0_2pi - 2D scan of theta and phi between 0 and 2pi for the FBS channel for epsilon = 6sigma_w and sigma_w = 0.1GHz

---

### KeyRates_DispersionChannel_v2.xlsx
Excel file containing all the data simulated by the security analysis for varying parameters for our protocol with the dispersion model (linear and quadratic).
- LinDisp_1.1GHz_30ps_0.19GHz - scan with alpha_1 linear dispersion parameter for sigma_w = 1.1GHz, Omega1 = 0.019THz, sigma_t = 30ps
- LinDisp_1.1GHz_17ps - scan with alpha_1 linear dispersion parameter for sigma_w = 1.1GHz, Omega1 = 0.019THz, sigma_t = 17ps
- LinDisp_1.1GHz_10ps - scan with alpha_1 linear dispersion parameter for sigma_w = 1.1GHz, Omega1 = 0.019THz, sigma_t = 10ps
- LinDisp_Omega1_0.01THz - scan with alpha_1 linear dispersion parameter for sigma_w = 1.1GHz, Omega1 = 0.01THz, sigma_t = 17ps
- LinDisp_Omega1_0.015THz - scan with alpha_1 linear dispersion parameter for sigma_w = 1.1GHz, Omega1 = 0.015THz, sigma_t = 17ps
- LinDisp_Omega1_0.02THz - scan with alpha_1 linear dispersion parameter for sigma_w = 1.1GHz, Omega1 = 0.02THz, sigma_t = 17ps
- LinDisp_Delta_1ps_sigma_30ps - scan with alpha_1 linear dispersion parameter for sigma_w = 1.1GHz, Omega1 = 0.01GHz, sigma_t = 30ps, for unequal dispersion on the two photons, separated by alpha_1 = 1ps
- LinDisp_Delta_5ps_sigma_30ps - scan with alpha_1 linear dispersion parameter for sigma_w = 1.1GHz, Omega1 = 0.01GHz, sigma_t = 30ps, for unequal dispersion on the two photons, separated by alpha_1 = 5ps
- QuadDisp_1.1GHz_30ps - scan with alpha_12 quadratic dispersion parameter for sigma_w = 1.1GHz, Omega1 = 0.019THz, sigma_t = 30ps
- QuadDisp_1.1GHz_17ps - scan with alpha_2 quadratic dispersion parameter for sigma_w = 1.1GHz, Omega1 = 0.019THz, sigma_t = 17ps
- QuadDisp_1.1GHz_10ps - scan with alpha_2 quadratic dispersion parameter for sigma_w = 1.1GHz, Omega1 = 0.019THz, sigma_t = 10ps
- QuadDisp_Omega1_0.01THz - scan with alpha_2 quadratic dispersion parameter for sigma_w = 1.1GHz, Omega1 = 0.01THz, sigma_t = 17ps
- QuadDisp_Omega1_0.015THz - scan with alpha_2 quadratic dispersion parameter for sigma_w = 1.1GHz, Omega1 = 0.015THz, sigma_t = 17ps
- QuadDisp_Omega1_0.02THz - scan with alpha_2 quadratic dispersion parameter for sigma_w = 1.1GHz, Omega1 = 0.02THz, sigma_t = 17ps
- QuadDisp_Delta_10ps2_sigma_30ps - scan with alpha_1 linear dispersion parameter for sigma_w = 1.1GHz, Omega1 = 0.01GHz, sigma_t = 30ps, for unequal dispersion on the two photons, separated by alpha_2 = 10ps^2
- QuadDisp_Delta_50ps2_sigma_30ps - scan with alpha_1 linear dispersion parameter for sigma_w = 1.1GHz, Omega1 = 0.01GHz, sigma_t = 30ps, for unequal dispersion on the two photons, separated by alpha_2 = 50ps^2

---

### TFQKD_KeyRates.xlsx
Excel file containing the simulated points of the analytical expression for the key rate of TF-QKD in a lossy channel. Used to add the TF-QKD plot to the loss comparison.

-----

## Key Rate Simulations

Our key rate analysis used the Version 1 of the Open QKD Security package, developed by A. Winick, et al., available at https://github.com/Optical-Quantum-Communication-Theory/ImprovedDecoyStateAndFlagStateSquashingMethods/tree/main/openQKDsecurityV1. Here, we provide the "Presets" and "Protocols" files that we developed for the analysis. These sit on top of and must be used in conjuction with the Open QKD Security package hyperlinked.

-----

### Presets

- BB84_asymptotic.m - file for BB84 key rate analysis. The BB84 key rate is analysed using the parameters theta and phi, with loss scaling eta^1. The same file is used for our protocol with the amplitude damping model, using only the parameters eta0 and delta for the amplitude damping and decoherence channel, and loss scaling eta^2.
- Boileau3_asymptotic.m - file for the key rate analysis of the 3-photon protocol from Boileau et al. (2004).
- Boileau4_asymptotic.m - file for the key rate analysis of the 4-photon protocol from Boileau et al. (2004).
- Li08Dep_asymptotic.m - file for the key rate analysis of the dephasing-only protocol from Li et al. (2008).
- Li08Rot_asymptotic.m - file for the key rate analysis of the rotation-only protocol from Li et al. (2008).
- Wang05_asymptotic.m - file for the key rate analysis of the protocol from Wang (2005).
- OurProt_asymptotic.m - file for the key rate analysis of our protocol using the ququart method of analysing the security.

---

### Protocols

For the BB84 code (and our protocol):
- BB84Channel1
- BB84LossyDescription
- BB84LossyNoisyChannel

For the Boileau 3-photon code:
- Boileau3LossyChannel
- Boileau3LossyDescription

For the Boileau 4-photon code:
- Boileau4LossyChannel
- Boileau4LossyDescription

For the Li dephasing-only code:
- Li08DepLossyChannel
- Li08DepLossyDescription

For the Li rotation-only code:
- Li08RotLossyChannel
- Li08RotLossyDescription

For the Wang code:
- WangLossyNoisyChannel
- WangLossyDescription

For our protocol with ququart model:
- OurProt_FBSChannel - FBS channel model
- OurProt_LinDispChannel - linear dispersion channel model
- OurProt_LossMixedFBSChannel - lossy channel with mixed FBS channel
- OurProt_LossyChannel - loss only channel
- OurProt_LossyDescription
- OurProt_MixedFBSChannel - mixed channel (Eq(8)) with FBS channel
- OurProt_MixedLinDispChannel - mixed channel (Eq(8)) with linear dispersion channel
- OurProt_MixedQuadDispChannel - mixed channel (Eq(8)) with quadratic dispersion channel
- OurProt_NonCollectiveLinDispChannel - linear dispersion with unequal effect on the two photons
- OurProt_NonCollectiveQuadDispChannel - quadratic dispersion with unequal effect on the two photons
- OurProt_QuadDispChannel - quadratic dispersion channel model

