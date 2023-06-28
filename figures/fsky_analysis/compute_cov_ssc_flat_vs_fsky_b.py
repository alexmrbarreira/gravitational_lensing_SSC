import sys; sys.path.append('../../'); 
from prepare_for_lenscov import *; 
from params_here import *

# ================================================================================ #
# Do variance integrals 
# ================================================================================ #
print ('Doing variance integrals for various fsky... ')

logpmin = log(0.02*100./c_light)
logpmax = log(1.0e5)

def var_integrand(logp,  chi,fsky):
    p       = exp(logp)
    z       = zofchi_int(chi)
    Omega_S = fsky * 4.*pi
    theta_S = sqrt(Omega_S/pi)
    x       = p * chi * theta_S
    out     = p**2. * Plin_int(z, p)[0] * (2.*special.j1(x)/x)**2.
    return out / (2.*pi)

variance_integrals = zeros([len(chi_array), len(ffsky)])
for j in range(len(ffsky)):
    fsky_now = ffsky[j]
    print (j, len(ffsky))
    for i in range(len(chi_array)):
        chi_now = chi_array[i]
        
        variance_integrals[i,j] = integrate.quad(var_integrand, logpmin, logpmax, args=(chi_now, fsky_now),  epsabs = 0.0, epsrel = epsrel_varint, limit = 1000)[0]


# ================================================================================ #
# Do chi integral
# ================================================================================ #
print ('Doing chi integrals ... ')

cov_ssc_vs_fsky = zeros(len(ffsky))
for i in range(len(ffsky)):
    print (i, len(ffsky))
    k1_now = (l1use_b+0.5)/chi_array
    k2_now = (l2use_b+0.5)/chi_array

    R_1_l1 = diagonal(flipud(array(R_1_int(z_array, k1_now))))
    R_1_l2 = diagonal(flipud(array(R_1_int(z_array, k2_now))))
    R_K_l1 = diagonal(flipud(array(R_K_int(z_array, k1_now))))
    R_K_l2 = diagonal(flipud(array(R_K_int(z_array, k2_now))))
    Pnl_l1 = diagonal(flipud(array(Pnl_int(z_array, k1_now))))
    Pnl_l2 = diagonal(flipud(array(Pnl_int(z_array, k2_now))))

    ssc_integrand_A      =   1.0   * variance_integrals[:,i] * R_1_l1 * R_1_l2 * Pnl_l1 * Pnl_l2 * gkernel_array**4. / chi_array**4.
    ssc_integrand_B      = (1./36) * variance_integrals[:,i] * R_K_l1 * R_K_l2 * Pnl_l1 * Pnl_l2 * gkernel_array**4. / chi_array**4.
    ssc_integrand_C      = (1./6)  * variance_integrals[:,i] * R_1_l1 * R_K_l2 * Pnl_l1 * Pnl_l2 * gkernel_array**4. / chi_array**4.
    ssc_integrand_D      = (1./6)  * variance_integrals[:,i] * R_K_l1 * R_1_l2 * Pnl_l1 * Pnl_l2 * gkernel_array**4. / chi_array**4.

    ssc_integral_A       = integrate.trapz(ssc_integrand_A, chi_array) 
    ssc_integral_B       = integrate.trapz(ssc_integrand_B, chi_array) 
    ssc_integral_C       = integrate.trapz(ssc_integrand_C, chi_array) 
    ssc_integral_D       = integrate.trapz(ssc_integrand_D, chi_array) 

    cov_ssc_vs_fsky[i] = ssc_integral_A + ssc_integral_B + ssc_integral_C + ssc_integral_D


print ('Writing ... ')
fout = open('data_cov_ssc_flat_vs_fsky_b.dat', 'w')
for i in range(len(ffsky)):
    fout.write(str(cov_ssc_vs_fsky[i])); fout.write('\n')
fout.close()

