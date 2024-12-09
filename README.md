# SurpassingLNTradeOffQKD
This repository provides code to reproduce the figures and simulations presented in the paper Surpassing the loss-noise robustness trade-off in quantum key distribution by Hannah Seabrook, Emilien Lavie, Teodor Stroemberg, Matthew P Stafford and Giulia Rubino. This README file provides detailed instructions on using the code to reproduce the data and figures in the paper.

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
Produces figures of the bit error rate of our protocol (Fig 5 in SM).

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
Produces the figure on collective noise with different noise parameters between the two photons for the dispersion model.
- keyrate1 - makes Fig 20a in SM (key rate under varying deviations in dispersion noise parameters for n=1)
- keyrate2 - makes Fig 20b in SM (key rate under varying deviations in dispersion noise parameters for n=2)

---

### UnequalFBSmodel.py
Produces the figure on collective noise with different noise parameters between the two photons for the FBS model.
- difftheta - makes Fig 19d in SM (fidelity for deviating theta between the two photons)
- diffepsilon - makes Fig 19a in SM (fidelity for deviating epsilon between the two photons)
- diffphi - makes Fig 19e in SM (fidelity for deviating phi between the two photons)
- diffmu - makes Fig 19c in SM (fidelity for deviating mu between the two photons)
- diffOmega - makes Fig 19b in SM (fidelity for deviating Omega between the two photons)


