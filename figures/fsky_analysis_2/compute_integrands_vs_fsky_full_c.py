import sys; sys.path.append('../../'); from prepare_for_lenscov import *; 
from params_here import *

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

def sigma_l1l2_L_full(l1,l2,L): #returns the value of sigma(l1,l2,L)
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

def sigma_l1l2_L_node(l1,l2,L): #returns the value of sigma(l1,l2,L) without the derivative term contributions
    def sigma_l1l2_L_integrand(logp,  l1,l2,L):
        p            = exp(logp)
        f_lLp_term_1 = f_lLp_int([l1, L, log10(p)])[0]
        f_lLp_term_2 = f_lLp_int([l2, L, log10(p)])[0]
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
# Compute cl_mask for the various f_sky and the sigma_full and sigma_node
# ================================================================================ #

print 'NOT EVALUATING MASK .. let that be done by other script ... '
#ccl_masks  = zeros([len(L_insum), len(ffsky)])
#for i in range(len(ffsky)):
#    print i
#    npix = hp.nside2npix(nside)
#    mask = zeros(npix)
#    mask[0:int(round(npix*ffsky[i]))] = 1.0 # Circular polar cap
#    cl_mask        = hp.anafast(mask,lmax = 4*nside)/(4.*pi*ffsky[i]*ffsky[i])
#    ccl_masks[:,i] = cl_mask[0:len(L_insum)]


print 'Computing sigma for l1 = l2 = ', l1use_c

sigma_flat = zeros(len(L_insum))
sigma_full = zeros(len(L_insum))
sigma_node = zeros(len(L_insum))
for i in range(len(L_insum)):
    print i
    sigma_flat[i] = sigma_l1l2_L_Limber(l1use_c, l2use_c, L_insum[i])
    sigma_full[i] = sigma_l1l2_L_full(l1use_c, l2use_c, L_insum[i])
    sigma_node[i] = sigma_l1l2_L_node(l1use_c, l2use_c, L_insum[i])

# ================================================================================ #
# Write to files
# ================================================================================ #
print 'Writing ... '

fout = open('data_sigma_flat_c.dat', 'w')
for i in range(len(L_insum)):
    fout.write(str(sigma_flat[i])); fout.write('\n')
fout.close()
fout = open('data_sigma_full_c.dat', 'w')
for i in range(len(L_insum)):
    fout.write(str(sigma_full[i])); fout.write('\n')
fout.close()
fout = open('data_sigma_node_c.dat', 'w')
for i in range(len(L_insum)):
    fout.write(str(sigma_node[i])); fout.write('\n')
fout.close()

