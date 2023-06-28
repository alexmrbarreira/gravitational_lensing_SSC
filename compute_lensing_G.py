from prepare_for_lenscov import *

# ================================================================================ #
# Compute lensing power spectrum 
# ================================================================================ #

print ('Computing power spectrum and Gaussian covariance ... ')
C_l1 = zeros(len(l1_array))
for m in range(len(l1_array)):
    z_array       = zofchi_int(chi_array)
    k1_now        = (l1_array[m]+0.5)/chi_array
    Pnl_l1        = diagonal(flipud(array(Pnl_int(z_array, k1_now))))
    C_l_integrand = Pnl_l1 * gkernel_array**2. / chi_array**2.
    C_l1[m]       = integrate.trapz(C_l_integrand, chi_array)

# Shape noise
Clnoise  = sigma_eps**2./sounumden/2.

cov_l1l2        = zeros([len(l1_array), len(l1_array)])
cov_l1l2_wnoise = zeros([len(l1_array), len(l1_array)])
for i in range(len(l1_array)):
    deltal_now           = deltal_array[i]
    l_now                = l1_array[i]
    cov_l1l2[i,i]        = 2. * (4.*pi/Omega_S/(2.*l_now+1)/deltal_now) *  C_l1[i]**2.
    cov_l1l2_wnoise[i,i] = 2. * (4.*pi/Omega_S/(2.*l_now+1)/deltal_now) * (C_l1[i] + Clnoise)**2.

# ================================================================================ #
# Symmetrize and write matrices
# ================================================================================ #

print ('Writing matrices and data vector ... ')

covwriter_matrix(cov_l1l2       , 'data_store/data_g_cov_l1l2.dat')
covwriter_matrix(cov_l1l2_wnoise, 'data_store/data_g_cov_l1l2_wnoise.dat')

fout = open('data_store/data_C_l1.dat', 'w')
for i in range(len(l1_array)):
    fout.write(str(l1_array[i])); fout.write(' '); fout.write(str(C_l1[i])); fout.write('\n');
fout.close()

