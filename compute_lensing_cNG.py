from prepare_for_lenscov import *

# ======================================================== #
# Define functions that'll be used throughout 
# ======================================================== #

epsilon_treelevel = 1.0e-9 #using this tiny margin to avoid have cosine angles = 1 (-1); I find this dirty trick works sufficiently well

# The integrand of the standard SPT tree level trispectrum 
def alpha(q,k,mu):
    return 1. + k*mu/q
def beta(q,k,mu):
    return (mu/2.)*(k/q + q/k) + mu**2.
def F2(q1,q2,mu12):
    return (1./7)*(5.*alpha(q1,q2,mu12) + 2.*beta(q1,q2,mu12))
def G2(q1,q2,mu12):
    return (1./7)*(3.*alpha(q1,q2,mu12) + 4.*beta(q1,q2,mu12))
def F2s(q1,q2,mu12):
    return (1./2)*(F2(q1,q2,mu12) + F2(q2,q1,mu12))
def G2s(q1,q2,mu12):
    return (1./2)*(G2(q1,q2,mu12) + G2(q2,q1,mu12))
def F3(q1,q2,q3,mu12,mu13,mu23):
    term1 = 2.* beta(q1,sqrt(q2**2.+q3**2.+2.*q2*q3*mu23),(q2*mu12+q3*mu13)/sqrt(q2**2.+q3**2.+2.*q2*q3*mu23)) * G2(q2,q3,mu23)
    term2 = 7.*alpha(q1,sqrt(q2**2.+q3**2.+2.*q2*q3*mu23),(q2*mu12+q3*mu13)/sqrt(q2**2.+q3**2.+2.*q2*q3*mu23)) * F2(q2,q3,mu23)
    term3 = 2.* beta(sqrt(q1**2.+q2**2.+2.*q1*q2*mu12),q3,(q1*mu13+q2*mu23)/sqrt(q1**2.+q2**2.+2.*q1*q2*mu12)) * G2(q1,q2,mu12)
    term4 = 7.*alpha(sqrt(q1**2.+q2**2.+2.*q1*q2*mu12),q3,(q1*mu13+q2*mu23)/sqrt(q1**2.+q2**2.+2.*q1*q2*mu12)) * G2(q1,q2,mu12)
    return (1./18)*(term1 + term2 + term3 + term4)
def F3s(q1,q2,q3,mu12,mu13,mu23):
    out = F3(q1,q2,q3,mu12,mu13,mu23) + F3(q1,q3,q2,mu13,mu12,mu23) + F3(q2,q1,q3,mu12,mu23,mu13) + F3(q2,q3,q1,mu23,mu12,mu13) + F3(q3,q1,q2,mu13,mu23,mu12) + F3(q3,q2,q1,mu23,mu13,mu12)
    return out/6.
def cov_stdtree_angle(theta, k1, k2, z):
    mu12 = cos(theta)
    # Shorthand notation of vector diferences and angles
    k1mk2          = sqrt(k1**2. + k2**2. - 2.*k1*k2*mu12) # |\vec{k1}-\vec{k2}|
    k1pk2          = sqrt(k1**2. + k2**2. + 2.*k1*k2*mu12) # |\vec{k1}+\vec{k2}|
    ang_k1mk2andk2 = (k1*mu12-k2)/k1mk2                    # Cosine angle between \vec{k1}-\vec{k2} and \vec{k2}
    ang_k1pk2andk2 = (k1*mu12+k2)/k1pk2                    # Cosine angle between \vec{k1}+\vec{k2} and \vec{k2}
    ang_k2mk1andk1 = (k2*mu12-k1)/k1mk2                    # Cosine angle between \vec{k2}-\vec{k1} and \vec{k1}
    ang_k2pk1andk1 = (k2*mu12+k1)/k1pk2                    # Cosine angle between \vec{k2}+\vec{k1} and \vec{k1}
    # Interpolate power spectra
    P1   = Plin_int(z, k1)[0]
    P2   = Plin_int(z, k2)[0]
    P1m2 = flipud(Plin_int(z, flipud(k1mk2))[:,0])
    P1p2 = Plin_int(z, k1pk2)[:,0]
    # Define terms
    term1 = 12.*F3s(k1,k1,k2,-1.+epsilon_treelevel,mu12,-mu12)             *  P1**2. * P2
    term2 = 12.*F3s(k2,k2,k1,-1.+epsilon_treelevel,mu12,-mu12)             *  P2**2. * P1
    term3 = 4.*F2s(k1mk2,k2, ang_k1mk2andk2)**2. * P2**2. * P1m2
    term4 = 4.*F2s(k1pk2,k2,-ang_k1pk2andk2)**2. * P2**2. * P1p2
    term5 = 4.*F2s(k1mk2,k1, ang_k2mk1andk1)**2. * P1**2. * P1m2
    term6 = 4.*F2s(k1pk2,k1,-ang_k2pk1andk1)**2. * P1**2. * P1p2
    term7 = 8.*F2s(k1mk2,k2, ang_k1mk2andk2)*F2s(k1mk2,k1, ang_k2mk1andk1)  *  P1 * P2 * P1m2
    term8 = 8.*F2s(k1pk2,k2,-ang_k1pk2andk2)*F2s(k1pk2,k1,-ang_k2pk1andk1)  *  P1 * P2 * P1p2
    return term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8

# ======================================================================================= #
# Define function covNG_resp(k1, k2, z) function -- actually, it is angle-averaged T(k1, -k2, k2, -k2), i.e., no 1/V factor
# Interpolators above are treated as global, but can be passed on as arguments if needed be
# ======================================================================================= #
def covNG_resp(k1, k2, z):
    fsq      = 0.5
    k_hard = max([k1, k2])
    k_soft = min([k1, k2])
    G1_hard = G1_int(z, k_hard)
    G1_soft = G1_int(z, k_soft)
    G2_hard = G2_int(z, k_hard)
    G2_soft = G2_int(z, k_soft)
    Plin_soft =  Plin_int(z, k_soft)
    Pnl_soft  =  Pnl_int(z, k_soft)
    Pnl_hard  =  Pnl_int(z, k_hard)
    # k*P'/P and k^2*P''/P
    fk1_hard = k_hard     *  dPnl_int(z, k_hard)/Pnl_int(z, k_hard)
    fk1_soft = k_soft     *  dPnl_int(z, k_soft)/Pnl_int(z, k_soft)
    fk2_hard = k_hard**2. * ddPnl_int(z, k_hard)/Pnl_int(z, k_hard)
    fk2_soft = k_soft**2. * ddPnl_int(z, k_soft)/Pnl_int(z, k_soft)
    # kG1'
    fractk = 0.02
    kG1prime_hard = k_hard * ( G1_int(z, (1.+fractk)*k_hard) - G1_int(z, (1.-fractk)*k_hard) ) / (2.*fractk*k_hard)
    kG1prime_soft = k_soft * ( G1_int(z, (1.+fractk)*k_soft) - G1_int(z, (1.-fractk)*k_soft) ) / (2.*fractk*k_soft)
    # Eulerian responses that are needed: R_2, R_Kdelta, R_K^2, R_K.K, R_KK
    R3e_hard = (8./21)*G1_hard + G2_hard + (-(2./9) - (2./3)*G1_hard)*fk1_hard + (1./9)*fk2_hard - (2./3)*kG1prime_hard
    R3e_soft = (8./21)*G1_soft + G2_soft + (-(2./9) - (2./3)*G1_soft)*fk1_soft + (1./9)*fk2_soft - (2./3)*kG1prime_soft
    R4e_hard = (1518./1813)*((8./21)*G1_hard + G2_hard) + (41./22)*(-2./9 - (2./3)*G1_hard)*fk1_hard + (1./3)*fk2_hard
    R4e_soft = (1518./1813)*((8./21)*G1_soft + G2_soft) + (41./22)*(-2./9 - (2./3)*G1_soft)*fk1_soft + (1./3)*fk2_soft
    R5e_hard = (1./21)*G1_hard - (1./6)*fk1_hard
    R5e_soft = (1./21)*G1_soft - (1./6)*fk1_soft
    R6e_hard = -(22./13)*G1_hard + (3./2)*fk1_hard
    R6e_soft = -(22./13)*G1_soft + (3./2)*fk1_soft
    R7e_hard = (1476./1813)*((8./21)*G1_hard + G2_hard) + (69./44)*(-2./9 - (2./3)*G1_hard)*fk1_hard + (1./2)*fk2_hard
    R7e_soft = (1476./1813)*((8./21)*G1_soft + G2_soft) + (69./44)*(-2./9 - (2./3)*G1_soft)*fk1_soft + (1./2)*fk2_soft
    # \mathcal{A, B, C} functions in Eq.(2.4) of arXiv:1705.01092
    A_hard = (1./2)*R3e_hard + (2./3)*R5e_hard + (2./9)*R6e_hard
    A_soft = (1./2)*R3e_soft + (2./3)*R5e_soft + (2./9)*R6e_soft
    B_hard = (2./3)*R4e_hard + (2./9)*R6e_hard
    B_soft = (2./3)*R4e_soft + (2./9)*R6e_soft
    C_hard = (4./9)*R7e_hard
    C_soft = (4./9)*R7e_soft
    # Compute stitched tree level piece
    theta_range = linspace(pi, 2.*pi, 1000)+1.0e-6 #Small correction avoids singularities; half the range suffices by symmetry (choice pi...2pi is to ensure current use of interpolators)
    if (k_soft/k_hard <= fsq): # Squeezed, so response
        covNG_stitched_tree = (2. * ( (A_hard + B_hard/4. + (11./32)*C_hard ) * Pnl_hard * Plin_soft**2.))[0]
    else: # Non-squeezed, so standard
        integrand_stdtree   = cov_stdtree_angle(theta_range, k1, k2, z)
        covNG_stitched_tree = integrate.trapz(integrand_stdtree, theta_range) / (pi)
    # Compute 1-loop response piece: the p-integral
    # Get knl
    p_array    = 10.**linspace(-3, 2, 1000)
    knl        = p_array[where(abs(p_array**3.*Plin_int(z, p_array)[:,0]/2./pi**2.-1.) == min(abs(p_array**3.*Plin_int(z, p_array)[:,0]/2./pi**2.-1.)))][0]
    p_array    = 10.**linspace(log10(0.01*min([fsq*k_soft, knl])), log10(min([fsq*k_soft, knl])), 1000)
    logp_array = log(p_array)
    p_integral = integrate.trapz(p_array**3.*Plin_int(z, p_array)[:,0]**2., logp_array)
    # Compute 1-loop response piece: include the angular part with responses
    covNG_resp1loop = (2.*A_hard*A_soft + B_hard*B_soft/10. + (2./5)*(A_hard*C_soft + A_soft*C_hard) + (B_hard*C_soft + B_soft*C_hard)/35.  + (27./280)*C_hard*C_soft)
    covNG_resp1loop *= (2. * (Pnl_hard*Pnl_soft/(2.*pi)**2.) * p_integral)[0]
    # Return the total result 
    return covNG_stitched_tree + covNG_resp1loop

# ================================================================================ #
# Do chi integrals 
# ================================================================================ #
print ('Doing chi integral ... ')

cov_l1l2 = zeros([len(l1_array), len(l2_array)])

# Loop over l1 and l2 values
for m in range(len(l1_array)):
    print (m, 'of', len(l1_array), '; it gets faster though ... ')
    for n in range(m, len(l2_array)):
        integrand_now = zeros(len(chi_array))
        for i in range(1,len(chi_array)):
            k1_now = (l1_array[m]+0.5)/chi_array[i]
            k2_now = (l2_array[n]+0.5)/chi_array[i]
            z_now  = z_array[i]

            integrand_now[i] = covNG_resp(k1_now, k2_now, min(z_now, 2.)) * gkernel_array[i]**4. / chi_array[i]**6.

        cov_l1l2[m,n] = integrate.trapz(integrand_now[1::], chi_array[1::]) / Omega_S

# ================================================================================ #
# Symmetrize and write matrices
# ================================================================================ #
print ('Symmetrizing and writing matrices ... ')

cov_l1l2 = symmetrize_matrix(cov_l1l2)

covwriter_matrix(cov_l1l2, 'data_store/data_cng_cov_l1l2.dat')

