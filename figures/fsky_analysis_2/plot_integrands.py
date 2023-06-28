import sys; sys.path.append('../../'); from prepare_for_lenscov import *; 
from params_here import *

# ================================================================================ #
# Define flat window function and compute values to plot 
# ================================================================================ #
def Wl2_flat(l, fsky):
    OmegaW = fsky * 4.*pi
    thetaW = sqrt(OmegaW/pi)
    x      = l * thetaW    
    out    = OmegaW**2. * (2.*special.j1(x)/x)**2.
    return out 

cls_masks_flat = zeros([len(L_insum), len(ffsky)])
for i in range(len(ffsky)):
    cls_masks_flat[:,i] = Wl2_flat(L_insum+1.0e-3, ffsky[i])

# ================================================================================ #
# Load files and compute stuff to plot
# ================================================================================ #

cov_ssc_vs_ffsky_flat_a = loadtxt('data_cov_ssc_flat_vs_fsky_a.dat')
cov_ssc_vs_ffsky_flat_b = loadtxt('data_cov_ssc_flat_vs_fsky_b.dat')
cov_ssc_vs_ffsky_flat_c = loadtxt('data_cov_ssc_flat_vs_fsky_c.dat')

sigma_flat_a = loadtxt('data_sigma_flat_a.dat')
sigma_full_a = loadtxt('data_sigma_full_a.dat')
sigma_node_a = loadtxt('data_sigma_node_a.dat')

sigma_flat_b = loadtxt('data_sigma_flat_b.dat')
sigma_full_b = loadtxt('data_sigma_full_b.dat')
sigma_node_b = loadtxt('data_sigma_node_b.dat')

sigma_flat_c = loadtxt('data_sigma_flat_c.dat')
sigma_full_c = loadtxt('data_sigma_full_c.dat')
sigma_node_c = loadtxt('data_sigma_node_c.dat')

cls_masks  = loadtxt('data_cls_mask.dat')

cov_ssc_vs_ffsky_full_a = zeros(len(ffsky))
cov_ssc_vs_ffsky_node_a = zeros(len(ffsky))
cov_ssc_vs_ffsky_full_b = zeros(len(ffsky))
cov_ssc_vs_ffsky_node_b = zeros(len(ffsky))
cov_ssc_vs_ffsky_full_c = zeros(len(ffsky))
cov_ssc_vs_ffsky_node_c = zeros(len(ffsky))

for i in range(len(ffsky)):
    cov_ssc_vs_ffsky_full_a[i] = sum( cls_masks[:,i] * sigma_full_a * (2.*L_insum + 1.) / (4.*pi) )
    cov_ssc_vs_ffsky_node_a[i] = sum( cls_masks[:,i] * sigma_node_a * (2.*L_insum + 1.) / (4.*pi) )

    cov_ssc_vs_ffsky_full_b[i] = sum( cls_masks[:,i] * sigma_full_b * (2.*L_insum + 1.) / (4.*pi) )
    cov_ssc_vs_ffsky_node_b[i] = sum( cls_masks[:,i] * sigma_node_b * (2.*L_insum + 1.) / (4.*pi) )

    cov_ssc_vs_ffsky_full_c[i] = sum( cls_masks[:,i] * sigma_full_c * (2.*L_insum + 1.) / (4.*pi) )
    cov_ssc_vs_ffsky_node_c[i] = sum( cls_masks[:,i] * sigma_node_c * (2.*L_insum + 1.) / (4.*pi) )

# ================================================================================ #
# Plot stuff ...  
# ================================================================================ #
labelsize = 24
ticksize  = 22
textsize  = 30
titlesize = 26
text_font = 24
legend_font = 20
v_min = 0.95
v_max = 1.15

fig0 = plt.figure(0, figsize=(17.5, 7.))
fig0.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom= 0.14, hspace = 0.20)#, space = 0.2)

ax = fig0.add_subplot(1,2,1)

plot(ffsky, cov_ssc_vs_ffsky_full_a/cov_ssc_vs_ffsky_flat_a, linestyle = 'solid', linewidth = 2., c = 'b', label = r"$\ell_1=\ell_2 = "+str("%.0f" % l1use_a)+"$")
plot(ffsky, cov_ssc_vs_ffsky_full_b/cov_ssc_vs_ffsky_flat_b, linestyle = 'solid', linewidth = 2., c = 'g', label = r"$\ell_1=\ell_2 = "+str("%.0f" % l1use_b)+"$")
plot(ffsky, cov_ssc_vs_ffsky_full_c/cov_ssc_vs_ffsky_flat_c, linestyle = 'solid', linewidth = 2., c = 'r', label = r"$\ell_1=\ell_2 = "+str("%.0f" % l1use_c)+"$")

#plot(ffsky, cov_ssc_vs_ffsky_node/cov_ssc_vs_ffsky_flat, linestyle = 'solid', linewidth = 2., c = 'g', label = r'$Without\ \partial^2/\partial x^2$')

fill_between(ffsky, ones(len(ffsky))*1.01, ones(len(ffsky))*0.99, color = 'grey', alpha = 0.5)

annotate(r"$Mask\ :\ Spherical\ cap$", xy = (0.06,0.60), xycoords='axes fraction', color = 'k', fontsize = text_font)
xlabel(r'$f_{\rm sky}$' , fontsize = labelsize)
ylabel(r'${\rm Cov}^{SSC}:\ \ {\rm curved\ sky}/{\rm flat\ sky}$' , fontsize = labelsize)
xscale('log')
xlim(min(ffsky), max(ffsky))
ylim(0.98, 1.30)
xticks(size = ticksize)
yticks(size = ticksize)
params = {'legend.fontsize': legend_font}; pylab.rcParams.update(params); legend(loc = 'upper left', ncol = 1)

ax = fig0.add_subplot(2,2,2)

#plot(L_insum+1, sigma_flat_a*1.0e21, linestyle = 'dashed', linewidth = 2., c = 'k', label = r'$\sigma_{\ell_1\ell_2}^{L, Limber}$')
#plot(L_insum+1, sigma_full_a*1.0e21, linestyle = 'solid', linewidth = 2., c = 'k', label = r'$\sigma_{\ell_1\ell_2}^{L}$')
#annotate(r"$\ell_1=\ell_2 = "+str("%.2f" % l1use_a)+"$", xy = (0.10,0.15), xycoords='axes fraction', color = 'k', fontsize = text_font)
#ylabel(r'$\times 10^{-21}$' , fontsize = labelsize)
#ylim(0.0, 2.2)

plot(L_insum+1, sigma_flat_b*1.0e23, linestyle = 'dashed' , linewidth = 2., c = 'k', label = r'$\sigma_{\ell_1\ell_2}^{L, {\rm flat}}$')
plot(L_insum+1, sigma_full_b*1.0e23, linestyle = 'solid', linewidth = 2., c = 'k', label = r'$\sigma_{\ell_1\ell_2}^{L}$')
annotate(r"$\ell_1=\ell_2 = "+str("%.0f" % l1use_b)+"$", xy = (0.10,0.15), xycoords='axes fraction', color = 'k', fontsize = text_font)
ylabel(r'$\times 10^{-23}$' , fontsize = labelsize)
ylim(0.0, 16.)

#plot(L_insum+1, sigma_flat_c*1.0e23, linestyle = 'dashed' , linewidth = 2., c = 'k', label = r'$\sigma_{\ell_1\ell_2}^{L, Limber}$')
#plot(L_insum+1, sigma_full_c*1.0e23, linestyle = 'solid', linewidth = 2., c = 'k', label = r'$\sigma_{\ell_1\ell_2}^{L}$')
#annotate(r"$\ell_1=\ell_2 = "+str("%.2f" % l1use_c)+"$", xy = (0.10,0.15), xycoords='axes fraction', color = 'k', fontsize = text_font)
#ylabel(r'$\times 10^{-23}$' , fontsize = labelsize)
#ylim(0.0, 2.7)

#ylabel(r'$\times 10^{-21}$' , fontsize = labelsize)
xscale('log')
xlim(min(L_insum), max(L_insum))
xticks(size = ticksize)
yticks(size = ticksize)
params = {'legend.fontsize': legend_font}; pylab.rcParams.update(params); legend(loc = 'upper right', ncol = 1)
ax.yaxis.set_major_locator(MaxNLocator(nbins=5))

ax = fig0.add_subplot(2,2,4)

clnorm = 1.0
ifsky  = [0, 50, -1]

plot(L_insum+1, cls_masks[:,[ifsky[0]]]*clnorm, linestyle = 'solid', linewidth = 2., c = 'm')
plot(L_insum+1, cls_masks[:,[ifsky[1]]]*clnorm, linestyle = 'solid', linewidth = 2., c = 'darkorange')
plot(L_insum+1, cls_masks[:,[ifsky[2]]]*clnorm, linestyle = 'solid', linewidth = 2., c = 'r')

plot(L_insum+1, cls_masks_flat[:,[ifsky[0]]] * clnorm / max(cls_masks_flat[:,[ifsky[0]]]) , linestyle = 'dashed', linewidth = 2., c = 'm')
plot(L_insum+1, cls_masks_flat[:,[ifsky[1]]] * clnorm / max(cls_masks_flat[:,[ifsky[1]]]) , linestyle = 'dashed', linewidth = 2., c = 'darkorange')
plot(L_insum+1, cls_masks_flat[:,[ifsky[2]]] * clnorm / max(cls_masks_flat[:,[ifsky[2]]]) , linestyle = 'dashed', linewidth = 2., c = 'r')

plot(L_insum+1, cls_masks[:,[ifsky[0]]]*clnorm*1.0e200, linestyle = 'dashed', linewidth = 2., c = 'k', label = r"$\propto \left|\mathcal{W}(L)\right|^2$")
plot(L_insum+1, cls_masks[:,[ifsky[0]]]*clnorm*1.0e200, linestyle = 'solid', linewidth = 2., c = 'k', label = r"$\propto \sum_{M}\left|b_{LM}\right|^2$")

annotate(r"$f_{\rm sky} = "+str("%.3f" % ffsky[ifsky[0]])+"$", xy = (0.70,0.20), xycoords='axes fraction', color = 'm'         , fontsize = text_font-2)
annotate(r"$f_{\rm sky} = "+str("%.2f" % ffsky[ifsky[1]])+"$", xy = (0.39,0.20), xycoords='axes fraction', color = 'darkorange', fontsize = text_font-2)
annotate(r"$f_{\rm sky} = "+str("%.1f" % ffsky[ifsky[2]])+"$", xy = (0.17,0.20), xycoords='axes fraction', color = 'r'         , fontsize = text_font-2)

xlabel(r'$L + 1$' , fontsize = labelsize)
ylabel(r'$Mask\ spectra$' , fontsize = labelsize)
xscale('log')
#xlim(min(L_insum), max(L_insum))
xlim(1., 200.)
ylim(-0.01, 1.1)
xticks(size = ticksize)
yticks(size = ticksize)
params = {'legend.fontsize': legend_font-2}; pylab.rcParams.update(params); legend(loc = 'upper right', ncol = 1)

ax.yaxis.set_major_locator(MaxNLocator(nbins=5))

fig0.savefig('fig_integrands.png')
fig0.savefig('fig_integrands.pdf')

show()
