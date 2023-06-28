import sys; sys.path.append('../'); from prepare_for_lenscov import *


# Set here index of l1_array array that this script will deal with ... 
index_l1 = 12
l1_touse = l1_array[index_l1]

# ================================================================================ #
# Compute table for f(l,L,p) (see notes and/or paper on lensing SSC)
# ================================================================================ #
print 'Constructing f(l,L,p) table ... '
print ' ...... doing ell = ', l1_touse, 'which is index', index_l1

# Get D(z) interpolator ; D(z=0)=1
k_dummy  = 0.1
Dofz     = sqrt(Plin_int(z_array, k_dummy) / Plin_int(0.0, k_dummy)[0])
Dofz_int = interpolate.interp1d(z_array, Dofz)

def jell(ell, x): #spherical Bessel function 
    return special.spherical_jn(int(ell), x, derivative=False)

def jell_primeprime(ell, x): #2nd derivative of the spherical Bessel function w.r.t. the argument
    L = int(ell)
    term0 = ((L*(L-1) - x**2.)/x**2.) * special.spherical_jn(int(ell)  , x, derivative=False)
    term1 = (2./x)                    * special.spherical_jn(int(ell+1), x, derivative=False)
    return term0 + term1

def f_lLp_integrand(chi,  ell,L,p): #chi-integrand of the f(l,L,p) function (without derivative terms)
    k_here        = (ell+0.5) / chi
    z_here        = zofchi_int(chi)
    R_perp        = R_1_int(z_here, k_here) + (1./6)*R_K_int(z_here, k_here)
    Pnl_here      = Pnl_int(z_here, k_here)
    sphBess       = jell(L, p*chi)
    return (gkernel_int(chi)**2. / chi**2.) * Pnl_here * Dofz_int(z_here) * R_perp * sphBess

def f_lLp_integrand_2deriv(chi,  ell,L,p): #chi-integrand of the f(l,L,p) function (the derivative terms)
    k_here        = (ell+0.5) / chi
    z_here        = zofchi_int(chi)
    R_perp        = R_1_int(z_here, k_here) + (1./6)*R_K_int(z_here, k_here)
    R_para        = R_1_int(z_here, k_here) - (1./3)*R_K_int(z_here, k_here)
    Pnl_here      = Pnl_int(z_here, k_here)
    sphBss_2deriv = jell_primeprime(L, p*chi)
    return (gkernel_int(chi)**2. / chi**2.) * Pnl_here * Dofz_int(z_here) * (R_perp-R_para) * sphBss_2deriv

def f_lLp(ell,L,log10p): #returns the value of f(l,L,p) (without the derivative terms)
    p = 10**log10p
    return integrate.quad(f_lLp_integrand, chi_array[0], chi_array[-1], args=(ell,L,p), epsabs = 0.0, epsrel = epsrel_flLp, limit = 1000)[0]

def f_lLp_2deriv(ell,L,log10p): #returns the value of f(l,L,p) (the derivative terms)
    p = 10**log10p
    return integrate.quad(f_lLp_integrand_2deriv, chi_array[0], chi_array[-1], args=(ell,L,p), epsabs = 0.0, epsrel = epsrel_flLp, limit = 1000)[0]

# Construct 3D lookup table for f(l,L,p)
f_lLp_table         = zeros([len(L_array), len(log10p_array_coarse)])
f_lLp_table_2deriv  = zeros([len(L_array), len(log10p_array_coarse)])

for i2 in range(len(L_array)):
    t0 = time.clock()
    print i2, 'of', len(L_array)
    for i3 in range(len(log10p_array_coarse)):

        f_lLp_table[i2,i3]        =        f_lLp(l1_touse, L_array[i2], log10p_array_coarse[i3])
        f_lLp_table_2deriv[i2,i3] = f_lLp_2deriv(l1_touse, L_array[i2], log10p_array_coarse[i3])

    t1 = time.clock()
    print 'This took ', t1-t0, 'sec'

print 'Writing table to file ... '
# Write table without the derivative terms
fout = open('data_f_lLp_lindex_' + str(index_l1) + '.dat', 'w')
for i2 in range(len(L_array)):
    for i3 in range(len(log10p_array_coarse)):
        fout.write(str(f_lLp_table[i2,i3])); fout.write(' ');
    fout.write('\n')
fout.close()

# Write table with the derivative terms only
fout = open('data_f_lLp_2deriv_lindex_' + str(index_l1) + '.dat', 'w')
for i2 in range(len(L_array)):
    for i3 in range(len(log10p_array_coarse)):
        fout.write(str(f_lLp_table_2deriv[i2,i3])); fout.write(' ');
    fout.write('\n')
fout.close()



