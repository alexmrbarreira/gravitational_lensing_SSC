from numpy import *
from scipy import *
from pylab import *
from scipy import interpolate, integrate, special
from scipy.interpolate import RegularGridInterpolator, UnivariateSpline
import pylab, time, warnings
import healpy as hp
import os
rcParams.update({'text.usetex': False, 'mathtext.fontset': 'stix'})
warnings.filterwarnings("ignore", category=DeprecationWarning)

lpath = '/home/barreira/a.Other/d-github_workrepos/b-gravitational_lensing_SSC/' 

# ==================================================== #
# Compute z(chi) relation
# ==================================================== #
print ('Computing z(chi) relation ... ')
Om0 = 0.30
h   = 0.7
c_light = 299792.458 #km/s

if (Om0 != 0.3):
    print ('')
    print ('Om0 =', Om0, 'but should be equal to 0.3 to match cosmology of tables in lookup_tables/')
    print ('Should you wish to modify the cosmology, you need to regenerate the spectra tables in lookup_tables/')
    print ('')
    quit()
elif (h!= 0.7):
    print ('')
    print ('h =', h, 'but should be equal to 0.7 to match cosmology of tables in lookup_tables/')
    print ('Should you wish to modify the cosmology, you need to regenerate the spectra tables in lookup_tables/')
    print ('')
    quit()

radian_to_arcmin = 180.*60./pi          #Multiply by this amount to convert radians to arcmin
radian_to_arcsec = radian_to_arcmin*60. #Multiply by this amount to convert radians to arcsec

# Source redshift
z_s      = 1.0
zz_array = linspace(0., 3.*z_s, 10000)

# Background cosmology functions
def E_lcdm(a, Om0): # This is E(a) = H(a)/H0
    return sqrt(Om0*a**(-3.) + (1-Om0))
def d_ang(z_i, z_f, Nz, E, Om0): # This gives the angular diameter distance in Mpc/h between z_i and z_f for given E = H(a)/H0
    array_int_z = linspace(z_i, z_f, Nz)
    d_ang_int   = 1./E(1./(1 + array_int_z), Om0)
    return (c_light/100)*integrate.trapz(d_ang_int, array_int_z)/(1 + z_f)

# Compute the chi(z)
chiofz = zeros(len(zz_array))
for i in range(len(zz_array)):
    chiofz[i] = d_ang(0.0, zz_array[i], 1000, E_lcdm, Om0)*(1. + zz_array[i])

chiofz_int = interpolate.interp1d(sort(zz_array), sort(chiofz))
zofchi_int = interpolate.interp1d(sort(chiofz), sort(zz_array))

# ================================================================================ #
# Load spectra and growth-only response tables and create 2D interpolators
# ================================================================================ #
print ('Making interpolators ... ')
# Create the interpolators that are needed; it does trick of setting to zero (or linear result for responses) for k = 0 (to prevent crashes/problems when k<k_min)
# This needs to be adapted if the format of the files changes
def make_interpolator_forP(path1, path2, value_at_k0):
    filetmp  = loadtxt(path1, skiprows = 1)
    ztmp     = loadtxt(path2)
    ktmp     = append(0., filetmp[:,0])
    Ptmp     = vstack([value_at_k0*ones(len(ztmp)), filetmp[:,1:]])
    return interpolate.interp2d(ztmp, ktmp, Ptmp, kind='cubic')

def make_interpolator_forG(path1, path2, value_at_k0, index_highk, asymptote, doextrapolation):
    filetmp  = loadtxt(path1, skiprows = 1)
    ztmp     = loadtxt(path2)
    # Set to linear theory for low-k
    ktmp     = append(0., filetmp[:,0])
    Gtmp     = vstack([value_at_k0*ones(len(ztmp)), filetmp[:,1:]])
    if(doextrapolation):
        # Set correction for high-k
        khig = 10.**linspace(log10(max(ktmp+0.000001)), 3, 500)
        kout = append(ktmp, khig)
        Ghig = zeros([len(kout), len(ztmp)])
        Ghig[0:len(ktmp), :] = Gtmp
        for i in range(len(ztmp)):
            G_at_kmax = Ghig[len(ktmp)-1, i]
            kmax      = ktmp[-1]
            Ghig[len(ktmp)::, i] = (kmax/khig)**index_highk * (G_at_kmax - asymptote) + asymptote
        return interpolate.interp2d(ztmp, kout, Ghig, kind='cubic')
    else:
        return interpolate.interp2d(ztmp, ktmp, Gtmp, kind='cubic')

# P_L(z, k)
Plin_int   = make_interpolator_forP(lpath+'lookup_tables/P_L_camb.dat', lpath+'lookup_tables/P_L_zvalues.dat', 0.0)
# dP_L/dk
dPlin_int  = make_interpolator_forP(lpath+'lookup_tables/dP_L_postcamb.dat', lpath+'lookup_tables/P_L_zvalues.dat', 0.0)
# d^2P_L/dk^2
ddPlin_int = make_interpolator_forP(lpath+'lookup_tables/ddP_L_postcamb.dat', lpath+'lookup_tables/P_L_zvalues.dat', 0.0)

## P_m(z, k) w/ Camb halofit
#Pnl_int   = make_interpolator_forP(lpath+'lookup_tables/P_m_camb.dat', lpath+'lookup_tables/P_m_zvalues.dat', 0.0)
## dP_m/dk w/ CLASS
#dPnl_int  = make_interpolator_forP(lpath+'lookup_tables/dP_m_postcamb.dat', lpath+'lookup_tables/P_m_zvalues.dat', 0.0)
## d^2P_L/dk^2 w/ CLASS
#ddPnl_int = make_interpolator_forP(lpath+'lookup_tables/ddP_m_postcamb.dat', lpath+'lookup_tables/P_m_zvalues.dat', 0.0)

# P_m(z, k) w/ Coyote
Pnl_int   = make_interpolator_forP(lpath+'lookup_tables/P_m_coyote.dat', lpath+'lookup_tables/P_m_zvalues.dat', 0.0)
# dP_m/dk w/ Coyote
dPnl_int  = make_interpolator_forP(lpath+'lookup_tables/dP_m_postcoyote.dat', lpath+'lookup_tables/P_m_zvalues.dat', 0.0)
# d^2P_L/dk^2 w/ Coyote
ddPnl_int = make_interpolator_forP(lpath+'lookup_tables/ddP_m_postcoyote.dat', lpath+'lookup_tables/P_m_zvalues.dat', 0.0)

# G_1(z, k)
G1_int   = make_interpolator_forG(lpath+'lookup_tables/Resp_G1_fromsims.dat', lpath+'lookup_tables/Resp_zvalues_fromsims.dat', 26./21, 0.5, -3./4, True)
# G_2(z, k)
G2_int   = make_interpolator_forG(lpath+'lookup_tables/Resp_G2_fromsims.dat', lpath+'lookup_tables/Resp_zvalues_fromsims.dat', 3002./1323, 000, 000, False)
# G_K(z, k)
GK_int   = make_interpolator_forG(lpath+'lookup_tables/Resp_GK_fromsims_Andreas.dat', lpath+'lookup_tables/Resp_zvalues_fromsims_Andreas.dat', 8./7, 0.5, -9./4, True)

# R_1(z, k)
def R_1_int(z, k):
    return 1. - (1./3)*k*dPnl_int(z, k)/Pnl_int(z, k) + G1_int(z, k)
# R_K(z, k)
def R_K_int(z, k):
    return GK_int(z, k) - k*dPnl_int(z, k)/Pnl_int(z, k) 

# ============================================================== #
# Specify lensing kernel, Omega_S, shape noise  and l1, l2 arrays
# ============================================================== #
chi_s         = chiofz_int(z_s)
chi_array     = linspace(0.01, chi_s, 1000)
z_array       = zofchi_int(chi_array)
a_array       = 1./(1.+z_array)
gkernel_array = (3.*(100.)**2*Om0/2./c_light**2.) * (chi_s - chi_array)*chi_array / chi_s / a_array
gkernel_int   = interpolate.interp1d(chi_array, gkernel_array) 

Omega_S       = 0.363610*4.*pi

sigma_eps = 0.37
sounumden = 30.*radian_to_arcmin**2. #in /radian^2

lbinedges    = 10.**linspace(log10(20.), log10(5000.), 21)
l1_array     = (lbinedges[:-1] + lbinedges[1:])/2.
l2_array     = (lbinedges[:-1] + lbinedges[1:])/2.
deltal_array = lbinedges[1:] - lbinedges[:-1]

# ============================================================== #
# Specify parameters of beyond Limber SSC calculation
# ============================================================== #

log10p_array_coarse = linspace(-4.0, log10(0.3), 100) # p values of f_lLp table
L_max_flLp          = 99
L_array             = array(range(int(L_max_flLp)+1))   # L values of f_lLp table 

L_max_insigma       = 15 # L up to which to use the beyond Limber result. Note L_max_insigma <= L_max_flLp must hold
if(L_max_insigma > L_max_flLp): print ('L_max_insigma > L_max_flLp! Reduce L_max_insigma. STOP!'); quit()
L_insum             = array(range(200))

epsrel_flLp  = 1.0e-2 # Relative accuracy of flLp integrals
epsrel_sigma = 1.0e-2 # Relative accuracy of sigma integrals

# ============================================================== #
# Other functions
# ============================================================== #

# Function that writes covariances to files
def covwriter_matrix(cov_2d_array, filepath):
    fout = open(filepath, 'w')
    leny = len(cov_2d_array[:,0])
    lenx = len(cov_2d_array[0,:])
    for i in range(leny):
        for j in range(lenx):
            fout.write(str(cov_2d_array[i,j])); fout.write(' ');
        fout.write('\n')
    fout.close()
    return 0

# Function that symmetrizes matrices (after computing only one of its triangles)
def symmetrize_matrix(matrix):
    if(len(matrix[:,0]) != len(matrix[0,:])):
        print ('Not square matrix!'); quit()
    n = len(matrix[:,0])
    new_matrix = copy(matrix)
    for j in range(n):
        for i in range(1+j, n):
            new_matrix[i, j] = matrix[j, i]
    return new_matrix


