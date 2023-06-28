# Super-sample covariance of lensing power spectra

This code was used in the numerical analysis of the paper:

- Barreira, Krause & Schmidt 2018a, https://arxiv.org/abs/1711.07467, *Complete super-sample lensing covariance in the response approach*

It evaluates the lensing power spectrum covariance matrix including all physical contributions: Gaussian (G), connected non-Gaussian (cNG, up to 1-loop terms) and super-sample covariance (SSC). For the G and cNG terms, it assumes the Limber approximation, but for the SSC term it does also a beyond-Limber calculation.

This code was subsequently incorporated into the cosmological likelihood analysis code [CosmoLike](https://github.com/CosmoLike/CosmoCov). This was later used in

- Barreira, Krause & Schmidt 2018b, https://arxiv.org/abs/1807.04266, *Accurate cosmic shear errors: do we need ensembles of simulations?*

to demonstrate the accuracy of analytical approaches to the covariance matrix for analyses of current and future large imaging galaxy surveys like DES, Euclid and Vera Rubin.


### Dependencies


- python: numpy, scipy, matplotlib
- python: healpy (pip install --user healpy), for Healpix-format mask operations

### Code overview

- The files parameters.py and functions.py define global parameters, variables and functions
- The scripts in compute_cov/ execute the covariance calculation
- The scripts in plots/ make plots (figures are stored here too).

### Gallery

Summary of the kinematic regimes in $k_1-k_2$ space

<img src="plots/fig_regimes_v2.png" width="1000" height=auto/>

Stitching of the tree-level standard perturbation theory (SPT) and response results

<img src="plots/fig_tree_transition_mono.png" width="1000" height=auto/>

Contributions along the diagonal

<img src="plots/fig_covariance_diagonal_mono.png" width="600" height=auto/>

Response functions used in the calculation

<img src="plots/fig_responses_eul.png" width="1000" height=auto/>
