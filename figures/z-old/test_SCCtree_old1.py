import sys; sys.path.append('../'); from prepare_for_lenscov import *
from numpy import random

## ==============================================================
## Define relevant functions
## ==============================================================
def alpha(q,k,mu):
    return 1. + k*mu/q
def beta(q,k,mu):
    return (mu/2.)*(k/q + q/k) + mu**2.
def F2(q1,q2,mu12):
    return (1./7)*(5.*alpha(q1,q2,mu12) + 2.*beta(q1,q2,mu12))
def F2s(q1,q2,mu12):
    return (1./2)*(F2(q1,q2,mu12) + F2(q2,q1,mu12))

## The SSC trispectrum contribution at tree level in perturbation theory
def T_ssc_tree_pt(k1, k2, mu12, p, mu, phi, z):
    # Define angles and vector amplitudes
    mu2         = mu*mu12 + sqrt(1.-mu**2.)*sqrt(1.-mu12**2.)*cos(phi) #cosine angle between \vec{k2} and \vec{p}
    mod_k1mp    = sqrt(k1**2. + p**2. - 2.*k1*p*mu ) #|\vec{k1} - \vec{p}|
    mod_k2pp    = sqrt(k2**2. + p**2. + 2.*k2*p*mu2) #|\vec{k2} - \vec{p}|
    ang_p_k1mp  = ( k1*mu  - p)/mod_k1mp # cosine angle between \vec{p} and  \vec{k1}-\vec{p}
    ang_p_mk2mp = (-k2*mu2 - p)/mod_k2pp # cosine angle between \vec{p} and -\vec{k2}-\vec{p}
    # Define the spectra that are needed
    Pp    = Plin_int(z, p)[0]
    Pk1   = Plin_int(z, k1)[0]
    Pk2   = Plin_int(z, k2)[0]
    Pk1mp = Plin_int(z, mod_k1mp)[0]
    Pk2pp = Plin_int(z, mod_k2pp)[0]
    # Define the four diagrams
    term1 = 4. * F2s(p, mod_k1mp, ang_p_k1mp) * F2s(p, mod_k2pp, ang_p_mk2mp) * Pp * Pk1mp * Pk2pp
    term2 = 4. * F2s(p,       k1,        -mu) * F2s(p, mod_k2pp, ang_p_mk2mp) * Pp * Pk1   * Pk2pp
    term3 = 4. * F2s(p, mod_k1mp, ang_p_k1mp) * F2s(p,       k2,         mu2) * Pp * Pk1mp * Pk2
    term4 = 4. * F2s(p,       k1,        -mu) * F2s(p,       k2,         mu2) * Pp * Pk1   * Pk2
    # That's it, return
    return term1 + term2 + term3 + term4

# The response based result at tree level
def T_ssc_tree_response(k1, k2, mu12, p, mu, phi, z):
    mu2         = mu*mu12 + sqrt(1.-mu**2.)*sqrt(1.-mu12**2.)*cos(phi) #cosine angle between \vec{k2} and \vec{p}
    # Define power spectrum terms
    Pp    = Plin_int(z, p)[0]
    Pk1   = Plin_int(z, k1)[0]
    Pk2   = Plin_int(z, k2)[0]
    # Define power spectrum derivative terms
    fk1_k1 = k1 * dPlin_int(z, k1)[0]/Plin_int(z, k1)[0]
    fk1_k2 = k2 * dPlin_int(z, k2)[0]/Plin_int(z, k2)[0]
    # Define the response terms
    curlyR1_k1 = (47./21 - fk1_k1/3.) + (8./7 - fk1_k1)*(mu**2. - 1./3)
    curlyR1_k2 = (47./21 - fk1_k2/3.) + (8./7 - fk1_k2)*(mu2**2. - 1./3)
    # That's it, return
    return curlyR1_k1 * curlyR1_k2 * Pp * Pk1 * Pk2 

## ====================================================================================================
## Fix all parameters, vary p, plot
## ====================================================================================================
# Choose parameters
drawrandom = True
if (drawrandom):
    print ''
    print 'Drawing them randomnly'
    print ''
    k1use = 10.**((2. + 1.)*random.rand() - 2.)
    k2use = 10.**((2. + 1.)*random.rand() - 2.)
    mu12use = 2.*random.rand() - 1.
    muuse   = 2.*random.rand() - 1.
    phiuse  = 2.*pi*random.rand()
    zuse    = 2.*random.rand()
else:
    k1use   = 0.1
    k2use   = 0.01
    mu12use = 0.98
    muuse   = -0.4
    phiuse  = pi*1.5
    zuse    = 0.234

print 'k1   = ', k1use
print 'k2   = ', k2use
print 'mu12 = ', mu12use
print 'mu   = ', muuse
print 'phi  = ', phiuse
print 'zuse = ', zuse

# Compute PT and response based results
p_array   = 10.**linspace(-5, 1, 1000)
result_pt_tree = zeros(len(p_array))
result_re_tree = zeros(len(p_array))
for i in range(len(p_array)):
    result_pt_tree[i] = T_ssc_tree_pt(      k1use, k2use, mu12use, p_array[i], muuse, phiuse, zuse)
    result_re_tree[i] = T_ssc_tree_response(k1use, k2use, mu12use, p_array[i], muuse, phiuse, zuse)

# Plot
labelsize = 22
ticksize  = 20
textsize  = 24
titlesize = 24
text_font = 20
legend_font = 20

fig1 = plt.figure(1, figsize=(8., 10.))
fig1.subplots_adjust(left=0.17, right=0.95, top=0.95, bottom=-0.20, hspace = 0.15)#, space = 0.2)
fig1.add_subplot(2,1,1)
plot(p_array, result_pt_tree / max(result_pt_tree), linestyle = 'solid', linewidth = 2, c = 'b', label = r'$Tree\ PT$')
plot(p_array, result_re_tree / max(result_re_tree), linestyle = 'solid', linewidth = 2, c = 'r', label = r'$Tree\ Response$')
axvline(k1use, linestyle = 'dashed', c = 'k')
axvline(k2use, linestyle = 'dashed', c = 'k')
ylabel(r'$T_m^{\rm SSC}/({\rm its\ max})$'       , fontsize = labelsize)
xlabel(r'$p\ \left[h/{\rm Mpc}\right]$'       , fontsize = labelsize)
xscale('log')
#yscale('log')
xlim(min(p_array), max(p_array))
#ylim(1.0e0, 1.0e4)
xticks(size = ticksize)
yticks(size = ticksize)
params = {'legend.fontsize': legend_font}; pylab.rcParams.update(params); legend(loc = 'upper right', ncol = 1)
fig1.add_subplot(4,1,3)
plot(p_array, result_re_tree / result_pt_tree - 1., linestyle = 'solid', linewidth = 2, c = 'r')
axhline(0.0  , linestyle = 'solid', c = 'b')
axvline(k1use, linestyle = 'dashed', c = 'k')
axvline(k2use, linestyle = 'dashed', c = 'k')
ylabel(r'$Rel.\ dif.\ to\ PT$'       , fontsize = labelsize)
xlabel(r'$p\ \left[h/{\rm Mpc}\right]$'       , fontsize = labelsize)
xscale('log')
#yscale('log')
xlim(min(p_array), max(p_array))
ylim(-0.50, 0.50)
xticks(size = ticksize)
yticks(size = ticksize)

## ====================================================================================================
## Generate many parameters, store result as function of p, plot ratios
## ====================================================================================================
Nrandom = 500
result_pt_tree_store = zeros([len(p_array), Nrandom])
result_re_tree_store = zeros([len(p_array), Nrandom])
mink_store = zeros(Nrandom)
for i in range(Nrandom):
    # Get values
    k1use = 10.**((2. + 1.)*random.rand() - 2.)
    k2use = 10.**((2. + 1.)*random.rand() - 2.)
    mu12use = 2.*random.rand() - 1.
    muuse   = 2.*random.rand() - 1.
    phiuse  = 2.*pi*random.rand()
    zuse    = 2.*random.rand()
    for j in range(len(p_array)):
        result_pt_tree_store[j, i] = T_ssc_tree_pt(      k1use, k2use, mu12use, p_array[j], muuse, phiuse, zuse)
        result_re_tree_store[j, i] = T_ssc_tree_response(k1use, k2use, mu12use, p_array[j], muuse, phiuse, zuse)
    mink_store[i] = min([k1use, k2use])

fig2 = plt.figure(2, figsize=(8., 6.))
fig2.subplots_adjust(left=0.17, right=0.95, top=0.95, bottom=0.13, hspace = 0.15)#, space = 0.2)
fig2.add_subplot(1,1,1)
for i in range(Nrandom):
    plot(p_array/mink_store[i], result_pt_tree_store[:,i]/result_re_tree_store[:,i] - 1., c = 'grey')
annotate(r"$From\ "+str(Nrandom)+"\ random\ \{k_1, k_2, \mu_{12}, \mu, \phi\}\ pairs$", xy = (0.04, 0.90), xycoords='axes fraction',
xytext = (0.04, 0.90), textcoords='axes fraction', bbox=dict(boxstyle="round4", fc = 'w', alpha = 0.6), color = 'k', fontsize = textsize-2)
ylabel(r'$T_m^{\rm SSC, PT}/T_m^{\rm SSC, Response}-1$'       , fontsize = labelsize)
xlabel(r'$p/k_{\rm soft}$'       , fontsize = labelsize)
xlim(0., 1.0)
ylim(-0.50, 0.50)
xticks(size = ticksize)
yticks(size = ticksize)

fig2.savefig('fig_SSC_pamp_test.png')

show()




