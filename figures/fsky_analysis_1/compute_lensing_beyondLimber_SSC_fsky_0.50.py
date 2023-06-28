import sys; sys.path.append('../../'); from prepare_for_lenscov import *; from params_here import *; Omega_S = 6.28318530718

# ================================================================================ #
# Draw mask and compute its angular power spectrum 
# ================================================================================ #
print 'Computing power spectrum of mask ... '
# Choose parameter nside, which sets number of pixels in the healpix map
nside = 128
npix  = hp.nside2npix(nside)

# array mask stores the value of the map in each pixel
fsky  = Omega_S/4./pi
mask = zeros(npix)

mask[0:int(round(npix*fsky))] = 1.0 # Circular polar cap

cl_mask  = hp.anafast(mask,lmax = 4*nside)/(4.*pi*fsky*fsky)
L_mask   = arange(len(cl_mask))

# ================================================================================ #
# Load f(l,L,p) (see notes and/or paper on lensing SSC)
# ================================================================================ #
print 'Loading f(l,L,p) table and creating interpolator ...  '

f_lLp_table         = zeros([len(l1_array), len(L_array), len(log10p_array_coarse)])
f_lLp_table_2deriv  = zeros([len(l1_array), len(L_array), len(log10p_array_coarse)])

for i1 in range(len(l1_array)):
    f_lLp_table[i1]        = loadtxt(lpath+'compute_lensing_beyondLimber_flLp_table/data_f_lLp_lindex_'        + str(i1) + '.dat')
    f_lLp_table_2deriv[i1] = loadtxt(lpath+'compute_lensing_beyondLimber_flLp_table/data_f_lLp_2deriv_lindex_' + str(i1) + '.dat')

f_lLp_int        = RegularGridInterpolator((l1_array, L_array, log10p_array_coarse), f_lLp_table       , method = 'nearest') #call like f_lLp_int([x,y,z])
f_lLp_2deriv_int = RegularGridInterpolator((l1_array, L_array, log10p_array_coarse), f_lLp_table_2deriv, method = 'nearest') 

# ================================================================================ #
# Define functions to compute sigma(ell_1, ell_2, L) 
# ================================================================================ #

def sigma_l1l2_L_Limber(l1,l2,L): #sigma_l1l2_L, but after Limber's approximation; used below for high-L speed-ups
    def sigma_l1l2_L_integrand_Limber(chi,  l1,l2,L): 
        z  = zofchi_int(chi)
        k1 = (l1+0.5)/chi
        k2 = (l2+0.5)/chi
        p  = (L+0.5)/chi
        weights   = gkernel_int(chi)**4. / chi**6.
        spectra   = Pnl_int(z, k1)*Pnl_int(z, k2)*Plin_int(z, p)
        responses = (R_1_int(z, k1) + (1./6)*R_K_int(z, k1)) * (R_1_int(z, k2) + (1./6)*R_K_int(z, k2))
        return weights * spectra * responses
    return integrate.quad(sigma_l1l2_L_integrand_Limber, chi_array[0], chi_array[-1], args=(l1,l2,L), epsabs = 0.0, epsrel = epsrel_sigma, limit = 1000)[0]

def sigma_l1l2_L(l1,l2,L): #returns the value of sigma(l1,l2,L)
    def sigma_l1l2_L_integrand(logp,  l1,l2,L): 
        p            = exp(logp)
        f_lLp_term_1 = f_lLp_int([l1, L, log10(p)])[0] + f_lLp_2deriv_int([l1, L, log10(p)])[0]
        f_lLp_term_2 = f_lLp_int([l2, L, log10(p)])[0] + f_lLp_2deriv_int([l2, L, log10(p)])[0]
        return (2./pi) * p**3. * f_lLp_term_1 * f_lLp_term_2 * Plin_int(0.0, p)[0]
    # Use full or Limber-approximated result depending on L-value
    if(L > L_max_insigma):
        out = sigma_l1l2_L_Limber(l1,l2,L)
    else:
        logp_min = log(10.**min(log10p_array_coarse))
        logp_max = log(10.**max(log10p_array_coarse))
        out = integrate.quad(sigma_l1l2_L_integrand, logp_min, logp_max, args=(l1,l2,L), epsabs = 0.0, epsrel = epsrel_sigma, limit = 1000)[0]
    return out

# ================================================================================ #
# Compute covariance matrix 
# ================================================================================ #
print 'Looping over l1, l2 ... '
# Loop over l1 and l2 values
cov_l1l2_beyondLimber_mono = zeros([len(l1_array), len(l2_array)])

for m in range(len(l1_array)):
    print m, 'of', len(l1_array), '; it gets faster though ... '
    for n in range(m, len(l2_array)):

        # Sum over L
        for iL in range(len(L_insum)):
            cov_l1l2_beyondLimber_mono[m,n] += cl_mask[iL] * sigma_l1l2_L(l1_array[m], l2_array[n], L_insum[iL]) * (2.*L_insum[iL] + 1.) / (4.*pi)

# ================================================================================ #
# Symmetrize and write matrices
# ================================================================================ #
print 'Symmetrizing and writing matrices and data vector ... '

cov_l1l2_beyondLimber_mono = symmetrize_matrix(cov_l1l2_beyondLimber_mono)

covwriter_matrix(cov_l1l2_beyondLimber_mono, 'data_ssc_cov_l1l2_beyondLimber_mono_fsky_0.50.dat')

