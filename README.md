# SurpassingLNTradeOffQKD
This repository provides code to reproduce the figures and simulations presented in the paper Surpassing the loss-noise robustness trade-off in quantum key distribution by Hannah Seabrook, Emilien Lavie, Teodor Stroemberg, Matthew P Stafford and Giulia Rubino. This README file provides detailed instructions on using the code to reproduce the data and figures in the paper.

BellStates.py - produces figures of the antisymmetric frequency-bin and time-bin Bell-states used for the encoding, i.e. Fig 2 in main text (3D plot) and Fig 1 in the SM (2D spectra)

FBSmodel.py - produces figures relating to the FBS noise model. 
              If comptau2 - makes Fig 3c in SM
              If comps_t - makes Fig 3d in SM
              If compOmega - makes Fig 3a in SM
              If compeps - makes Fig 4 in the main text
              If compmu - makes Fig 3b in SM
              If compsigma - makes Fig 8 in Appendix C
              If phasefac - makes Fig 4 in SM
              If phasefac and krcomp - makes Fig 7 in SM and Fig 5 in the main text
              If savedat - generates and saves data used for the probabilistic FBS channel analysis outlined in Sec III A of the main text


