import sys; sys.path.append('../../'); from prepare_for_lenscov import *; from params_here import *; Omega_S = 2.51327412287

# ================================================================================ #
# Do variance integrals 
# ================================================================================ #
print 'Doing variance integrals ... '
pperp_array   = 10.**linspace(log10(0.2*100./c_light), 3.0, 10000)

# Do variance integrals
varint_A_chi = zeros(len(chi_array))
for i in range(len(chi_array)):
    theta_S = sqrt(Omega_S/pi)
    xnow    = pperp_array*chi_array[i]*theta_S
    z_now   = zofchi_int(chi_array[i])

    integrand_A_now = pperp_array**2. * Plin_int(z_now, pperp_array)[:,0] * (2.*special.j1(xnow)/xnow)**2.

    varint_A_chi[i] = integrate.trapz(integrand_A_now, log(pperp_array))/2./pi

# For isotropic window functions and averaging also over theta_12, all variance integrals are related in this way (see mathematica notebook)
varint_B_chi = (1./36) * varint_A_chi
varint_C_chi =  (1./6) * varint_A_chi
varint_D_chi =  (1./6) * varint_A_chi

# ================================================================================ #
# Do chi integrals 
# ================================================================================ #
print 'Doing chi integrals ... '

cov_l1l2_A = zeros([len(l1_array), len(l2_array)])
cov_l1l2_B = zeros([len(l1_array), len(l2_array)])
cov_l1l2_C = zeros([len(l1_array), len(l2_array)])
cov_l1l2_D = zeros([len(l1_array), len(l2_array)])

# Loop over l1 and l2 values
for m in range(len(l1_array)):
    print m, 'of', len(l1_array), '; it gets faster though ... '
    for n in range(m, len(l2_array)):
        k1_now = (l1_array[m]+0.5)/chi_array
        k2_now = (l2_array[n]+0.5)/chi_array
        R_1_l1 = diagonal(flipud(array(R_1_int(z_array, k1_now))))
        R_1_l2 = diagonal(flipud(array(R_1_int(z_array, k2_now))))
        R_K_l1 = diagonal(flipud(array(R_K_int(z_array, k1_now))))
        R_K_l2 = diagonal(flipud(array(R_K_int(z_array, k2_now))))
        Pnl_l1 = diagonal(flipud(array(Pnl_int(z_array, k1_now))))
        Pnl_l2 = diagonal(flipud(array(Pnl_int(z_array, k2_now))))

        ssc_integrand_A =  varint_A_chi * R_1_l1 * R_1_l2 * Pnl_l1 * Pnl_l2 * gkernel_array**4. / chi_array**4.
        ssc_integrand_B =  varint_B_chi * R_K_l1 * R_K_l2 * Pnl_l1 * Pnl_l2 * gkernel_array**4. / chi_array**4.
        ssc_integrand_C =  varint_C_chi * R_1_l1 * R_K_l2 * Pnl_l1 * Pnl_l2 * gkernel_array**4. / chi_array**4.
        ssc_integrand_D =  varint_D_chi * R_K_l1 * R_1_l2 * Pnl_l1 * Pnl_l2 * gkernel_array**4. / chi_array**4.

        cov_l1l2_A[m,n] = integrate.trapz(ssc_integrand_A, chi_array)
        cov_l1l2_B[m,n] = integrate.trapz(ssc_integrand_B, chi_array)
        cov_l1l2_C[m,n] = integrate.trapz(ssc_integrand_C, chi_array)
        cov_l1l2_D[m,n] = integrate.trapz(ssc_integrand_D, chi_array)

# ================================================================================ #
# Symmetrize and write matrices
# ================================================================================ #
print 'Symmetrizing and writing matrices and data vector ... '

cov_l1l2_A = symmetrize_matrix(cov_l1l2_A)
cov_l1l2_B = symmetrize_matrix(cov_l1l2_B)
cov_l1l2_C = symmetrize_matrix(cov_l1l2_C)
cov_l1l2_D = symmetrize_matrix(cov_l1l2_D)

cov_l1l2 = cov_l1l2_A + cov_l1l2_B + cov_l1l2_C + cov_l1l2_D

covwriter_matrix(cov_l1l2, 'data_ssc_cov_l1l2_fsky_0.20.dat')
