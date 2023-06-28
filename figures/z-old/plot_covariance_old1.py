import sys; sys.path.append('../'); from prepare_for_lenscov import *

## ========================================================
## Prepare what to plot
## ========================================================

data_vector = loadtxt('../data_store/data_C_l1.dat')[:,1]
matrix_C1C2 = outer(data_vector, data_vector)

cov_ssc_a = loadtxt('../data_store/data_ssc_cov_l1l2_A.dat')
cov_ssc_b = loadtxt('../data_store/data_ssc_cov_l1l2_B.dat')
cov_ssc_c = loadtxt('../data_store/data_ssc_cov_l1l2_C.dat')
cov_ssc_d = loadtxt('../data_store/data_ssc_cov_l1l2_D.dat')

cov_ssc     = cov_ssc_a + cov_ssc_b + cov_ssc_c + cov_ssc_d
cov_ssc_bcd = cov_ssc_b + cov_ssc_c + cov_ssc_d
cov_cng   = loadtxt('../data_store/data_cng_cov_l1l2.dat')
#cov_g     = loadtxt('../data_store/data_g_cov_l1l2.dat') 
cov_g     = loadtxt('../data_store/data_g_cov_l1l2_wnoise.dat') 
#cov_gwn   = loadtxt('../data_store/data_g_cov_l1l2_wnoise.dat') 
cov_tot   = cov_ssc + cov_cng + cov_g

invcov_g           = linalg.inv(cov_g)
invcov_g_ssc_a     = linalg.inv(cov_g + cov_ssc_a)
invcov_g_ssc_a_cng = linalg.inv(cov_g + cov_ssc_a + cov_cng)
invcov_g_ssc_abcd  = linalg.inv(cov_g + cov_ssc)
invcov_g_cng       = linalg.inv(cov_g + cov_cng)
invcov_g_ssc_bcd   = linalg.inv(cov_g + cov_ssc_b + cov_ssc_c + cov_ssc_d)
invcov_tot         = linalg.inv(cov_tot)

## ========================================================
## Compute Signal-to-Noise ratio as a function of ell 
## ========================================================
snratio_g           = zeros(len(l1_array))
snratio_g_ssc_a     = zeros(len(l1_array))
snratio_g_ssc_a_cng = zeros(len(l1_array))
snratio_g_ssc_abcd  = zeros(len(l1_array))
snratio_g_cng       = zeros(len(l1_array))
snratio_g_ssc_bcd   = zeros(len(l1_array))
snratio_tot         = zeros(len(l1_array))

for i in range(len(l1_array)):
    signal_now             = data_vector[:i+1]
    invcov_g_now           = invcov_g[:i+1, :i+1]
    invcov_g_ssc_a_now     = invcov_g_ssc_a[:i+1, :i+1]
    invcov_g_ssc_a_cng_now = invcov_g_ssc_a_cng[:i+1, :i+1]
    invcov_g_ssc_abcd_now  = invcov_g_ssc_abcd[:i+1, :i+1]
    invcov_g_cng_now       = invcov_g_cng[:i+1, :i+1]
    invcov_g_ssc_bcd_now   = invcov_g_ssc_bcd[:i+1, :i+1]
    invcov_tot_now         = invcov_tot[:i+1, :i+1]

    snratio_g[i]           = sqrt(dot(dot(signal_now, invcov_g_now)            , signal_now))
    snratio_g_ssc_a[i]     = sqrt(dot(dot(signal_now, invcov_g_ssc_a_now)      , signal_now))
    snratio_g_ssc_a_cng[i] = sqrt(dot(dot(signal_now, invcov_g_ssc_a_cng_now)  , signal_now))
    snratio_g_ssc_abcd[i]  = sqrt(dot(dot(signal_now, invcov_g_ssc_abcd_now)   , signal_now))
    snratio_g_cng[i]       = sqrt(dot(dot(signal_now, invcov_g_cng_now)        , signal_now))
    snratio_g_ssc_bcd[i]   = sqrt(dot(dot(signal_now, invcov_g_ssc_bcd_now)    , signal_now))
    snratio_tot[i]         = sqrt(dot(dot(signal_now, invcov_tot_now)          , signal_now))

## ========================================================
## Plot
## ========================================================
labelsize = 22
ticksize  = 20
textsize  = 24
titlesize = 24
text_font = 20
legend_font = 20
v_min = 0.01
v_max = 0.60

fig0 = plt.figure(0, figsize=(17., 7.))
fig0.subplots_adjust(left=0.05, right=1.00, top=0.93, bottom=0.14, wspace = 0.15, hspace = 0.15)

panel = fig0.add_subplot(1,2,1)
#pcolor(l1_array, l2_array, 1.0e4*cov_ssc/matrix_C1C2)
pcolor(l1_array, l2_array, 1.0e4*cov_ssc/matrix_C1C2, vmin=v_min, vmax=v_max)
panel.set_aspect('equal')
colorbar().ax.tick_params(labelsize=ticksize-2)
annotate(r"${\rm SSC}$", xy = (0.04, 0.90), xycoords='axes fraction',
xytext = (0.04, 0.90), textcoords='axes fraction', bbox=dict(boxstyle="round4", fc = 'w', alpha = 0.6), color = 'k', fontsize = textsize+4)
title(r'${\rm Cov}(\ell_1, \ell_2)/(C(\ell_1)C(\ell_2)) \times 10^4$'   , fontsize = labelsize-2)
xlabel(r'$\ell_1$' , fontsize = labelsize)
ylabel(r'$\ell_2$' , fontsize = labelsize)
xscale('log')
yscale('log')
xlim(min(l1_array), max(l1_array))
ylim(min(l2_array), max(l2_array))
xticks(size = ticksize)
yticks(size = ticksize)
panel = fig0.add_subplot(1,2,2)
#pcolor(l1_array, l2_array, 1.0e4*cov_cng/matrix_C1C2)
pcolor(l1_array, l2_array, 1.0e4*(cov_cng+cov_g)/matrix_C1C2, vmin=v_min, vmax=v_max)
panel.set_aspect('equal')
colorbar().ax.tick_params(labelsize=ticksize-2)
annotate(r"${\rm G\ +\ cNG}$", xy = (0.04, 0.90), xycoords='axes fraction',
xytext = (0.04, 0.90), textcoords='axes fraction', bbox=dict(boxstyle="round4", fc = 'w', alpha = 0.6), color = 'k', fontsize = textsize+4)
title(r'${\rm Cov}(\ell_1, \ell_2)/(C(\ell_1)C(\ell_2)) \times 10^4$'   , fontsize = labelsize-2)
xlabel(r'$\ell_1$' , fontsize = labelsize)
ylabel(r'$\ell_2$' , fontsize = labelsize-2)
xscale('log')
yscale('log')
xlim(min(l1_array), max(l1_array))
ylim(min(l2_array), max(l2_array))
xticks(size = ticksize)
yticks(size = ticksize)

fig1 = plt.figure(1, figsize=(17.5, 10.))
fig1.subplots_adjust(left=0.07, right=0.99, top=0.96, bottom=0.10, wspace = 0.16, hspace = 0.13)
l2indices = [int(len(l1_array)*0.), int(len(l1_array)*0.33), int(len(l1_array)*0.66), len(l1_array)-1] 
for i in range(4):
    panel = fig1.add_subplot(2,2,i+1)
    l2index = l2indices[i]
    plot(l1_array, 1.0e4*cov_ssc_a[:,l2index]  /data_vector/data_vector[l2index], linestyle = 'solid' , linewidth = 2, c = 'g', label = r'${\rm Cov}^{\rm SSC}_{\kappa, A}$')
    plot(l1_array, 1.0e4*cov_ssc_bcd[:,l2index]/data_vector/data_vector[l2index], linestyle = 'dashed', linewidth = 2, c = 'g', label = r'${\rm Cov}^{\rm SSC}_{\kappa, BCD}$')
    #plot(l1_array, 1.0e4*cov_ssc_b[:,l2index]  /data_vector/data_vector[l2index], linestyle = '-.'    , linewidth = 2, c = 'g', label = r'${\rm Cov}^{\rm SSC}_{\kappa, B}$')
    #plot(l1_array, 1.0e4*cov_ssc_c[:,l2index]  /data_vector/data_vector[l2index], linestyle = 'dashed', linewidth = 2, c = 'g', label = r'${\rm Cov}^{\rm SSC}_{\kappa, C}$')
    #plot(l1_array, 1.0e4*cov_ssc_d[:,l2index]  /data_vector/data_vector[l2index], linestyle = 'dotted', linewidth = 2, c = 'g', label = r'${\rm Cov}^{\rm SSC}_{\kappa, D}$')
    plot(l1_array, 1.0e4*cov_ssc[:,l2index]    /data_vector/data_vector[l2index], linestyle = 'solid' , linewidth = 2, c = 'b', label = r'${\rm Cov}^{\rm SSC}_{\kappa}$')
    plot(l1_array, 1.0e4*cov_cng[:,l2index]    /data_vector/data_vector[l2index], linestyle = 'solid' , linewidth = 2, c = 'r', label = r'${\rm Cov}^{\rm cNG}_{\kappa}$')
    plot(l1_array, 1.0e4*cov_tot[:,l2index]    /data_vector/data_vector[l2index], linestyle = 'solid' , linewidth = 2, c = 'k', label = r'${\rm Cov}^{\rm Total}_{\kappa}$')
    #plot(l1_array,1.0e4*(cov_tot-cov_ssc_bcd)[:,l2index]/data_vector/data_vector[l2index],linestyle='dashed',linewidth=1,c='k',label=r'${\rm Cov}^{\rm Total*}_{\kappa}$')
    annotate(r'$\ell_2\ =\ ' +  str("%.0f" % l2_array[l2index]) + '$', xy = (0.80,0.86), xycoords='axes fraction', color = 'k', fontsize = text_font+1)
    ylabel(r'${\rm Cov}_{12}/(C_1C_2) \times 10^4$'        , fontsize = labelsize)
    xscale('log')
    xlim(min(l1_array), max(l1_array))
    ylim(-0.05, 0.65)
    xticks(size = ticksize)
    yticks(size = ticksize)
    if(i>1):
        xlabel(r'$\ell_1$' , fontsize = labelsize)
    if(i==3):
        params = {'legend.fontsize': legend_font}; pylab.rcParams.update(params); legend(loc = 'upper left', ncol = 2)

fig2 = plt.figure(2, figsize=(17.5, 10.))
fig2.subplots_adjust(left=0.09, right=0.98, top=0.95, bottom=-0.20, hspace = 0.15)#, space = 0.2)
fig2.add_subplot(2,2,1)
plot(l1_array, 1.0e4*diagonal(cov_ssc_a  /matrix_C1C2), linestyle = 'solid' , linewidth = 2, c = 'g', label = r'${\rm Cov}^{\rm SSC}_{\kappa, A}$')
plot(l1_array, 1.0e4*diagonal(cov_ssc_bcd/matrix_C1C2), linestyle = 'dashed', linewidth = 2, c = 'g', label = r'${\rm Cov}^{\rm SSC}_{\kappa, BCD}$')
#plot(l1_array,-1.0e4*diagonal(cov_ssc_b  /matrix_C1C2), linestyle = '-.'    , linewidth = 2, c = 'g', label = r'$-{\rm Cov}^{\rm SSC}_{\kappa, B}$')
#plot(l1_array, 1.0e4*diagonal(cov_ssc_c  /matrix_C1C2), linestyle = 'dashed', linewidth = 2, c = 'g', label = r'${\rm Cov}^{\rm SSC}_{\kappa, C}$')
#plot(l1_array, 1.0e4*diagonal(cov_ssc_d  /matrix_C1C2), linestyle = 'dotted', linewidth = 2, c = 'g', label = r'${\rm Cov}^{\rm SSC}_{\kappa, D}$')
plot(l1_array, 1.0e4*diagonal(cov_ssc    /matrix_C1C2), linestyle = 'solid' , linewidth = 2, c = 'b', label = r'${\rm Cov}^{\rm SSC}_{\kappa}$')
plot(l1_array, 1.0e4*diagonal(cov_cng    /matrix_C1C2), linestyle = 'solid' , linewidth = 2, c = 'r', label = r'${\rm Cov}^{\rm cNG}_{\kappa}$')
plot(l1_array, 1.0e4*diagonal(cov_g      /matrix_C1C2), linestyle = 'solid' , linewidth = 2, c = 'c', label = r'${\rm Cov}^{\rm G}_{\kappa}$')
plot(l1_array, 1.0e4*diagonal(cov_tot    /matrix_C1C2), linestyle = 'solid' , linewidth = 2, c = 'k', label = r'${\rm Cov}^{\rm Total}_{\kappa}$')
#plot(l1_array, 1.0e4*diagonal((cov_tot-cov_ssc_bcd)  /matrix_C1C2), linestyle = 'dashed' , linewidth = 1, c = 'k', label = r'${\rm Cov}^{\rm Total*}_{\kappa}$')
ylabel(r'${\rm Cov}(\ell, \ell)/[C(\ell)]^2 \times 10^4$'        , fontsize = labelsize)
xscale('log')
yscale('log')
xlim(min(l1_array), max(l1_array))
ylim(1.0e-3, 1.0e3)
xticks(size = ticksize)
yticks(size = ticksize)
params = {'legend.fontsize': legend_font}; pylab.rcParams.update(params); legend(loc = 'upper right', ncol = 3)
fig2.add_subplot(4,2,5)
plot(l1_array, diagonal(cov_g)       / diagonal(cov_tot), linestyle = 'solid'  , linewidth = 2, c = 'c')
plot(l1_array, diagonal(cov_ssc)     / diagonal(cov_tot), linestyle = 'solid'  , linewidth = 2, c = 'b')
plot(l1_array, diagonal(cov_ssc_a)   / diagonal(cov_tot), linestyle = 'solid'  , linewidth = 2, c = 'g')
plot(l1_array, diagonal(cov_ssc_bcd) / diagonal(cov_tot), linestyle = 'dashed' , linewidth = 2, c = 'g')
plot(l1_array, diagonal(cov_cng)     / diagonal(cov_tot), linestyle = 'solid'  , linewidth = 2, c = 'r')
axhline(1.0, linestyle = 'solid' , linewidth = 2, c = 'k',)
xlabel(r'$\ell$' , fontsize = labelsize)
ylabel(r'$Fractional\ contribution$'        , fontsize = labelsize)
xscale('log')
xlim(min(l1_array), max(l1_array))
ylim(0., 1.03)
xticks(size = ticksize)
yticks(size = ticksize)

fig2 = plt.figure(2, figsize=(17.5, 10.))
fig2.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=-0.20, hspace = 0.15)#, space = 0.2)
fig2.add_subplot(2,2,2)
plot(l1_array, snratio_g          , linestyle = 'dashed', linewidth = 2, c = 'k', label = r'$G$')
plot(l1_array, snratio_g_cng      , linestyle = 'solid' , linewidth = 2, c = 'orange', label = r'$G+cNG$')
plot(l1_array, snratio_g_ssc_bcd  , linestyle = 'dashed', linewidth = 2, c = 'orange', label = r'$G+SSC_{BCD}$')
plot(l1_array, snratio_g_ssc_a    , linestyle = 'solid' , linewidth = 2, c = 'r', label = r'$G+SSC_A$')
plot(l1_array, snratio_g_ssc_a_cng, linestyle = 'solid' , linewidth = 2, c = 'm', label = r'$G+SSC_A+cNG$')
plot(l1_array, snratio_g_ssc_abcd , linestyle = 'dashed', linewidth = 2, c = 'm', label = r'$G+SSC$')
plot(l1_array, snratio_tot        , linestyle = 'solid' , linewidth = 2, c = 'k', label = r'$G+cNG+SSC$')
ylabel(r'$Signal-to-Noise(\ell<\ell_{\rm max})$'       , fontsize = labelsize)
xscale('log')
yscale('log')
xlim(min(l1_array), max(l1_array))
#ylim(1.0e0, 1.0e4)
xticks(size = ticksize)
yticks(size = ticksize)
params = {'legend.fontsize': legend_font}; pylab.rcParams.update(params); legend(loc = 'upper left', ncol = 1)
fig2.add_subplot(4,2,6)
plot(l1_array, snratio_g           / snratio_g - 1., linestyle = 'dashed', linewidth = 2, c = 'k')
plot(l1_array, snratio_g_cng       / snratio_g - 1., linestyle = 'solid' , linewidth = 2, c = 'orange')
plot(l1_array, snratio_g_ssc_bcd   / snratio_g - 1., linestyle = 'dashed', linewidth = 2, c = 'orange')
plot(l1_array, snratio_g_ssc_a     / snratio_g - 1., linestyle = 'solid' , linewidth = 2, c = 'r')
plot(l1_array, snratio_g_ssc_a_cng / snratio_g - 1., linestyle = 'solid' , linewidth = 2, c = 'm')
plot(l1_array, snratio_g_ssc_abcd  / snratio_g - 1., linestyle = 'dashed', linewidth = 2, c = 'm')
plot(l1_array, snratio_tot         / snratio_g - 1., linestyle = 'solid' , linewidth = 2, c = 'k')
axhline(0.0, linestyle = 'dashed' , linewidth = 1, c = 'k',)
xlabel(r'$\ell_{\rm max}$' , fontsize = labelsize)
ylabel(r'$Rel.\ diff.\ to\ G$'        , fontsize = labelsize)
xscale('log')
xlim(min(l1_array), max(l1_array))
ylim(-0.70, 0.03)
xticks(size = ticksize)
yticks(size = ticksize)

fig0.savefig('fig_cov_l1l2_maps.png')
fig0.savefig('../../lensing_response/paper/figures/fig_cov_l1l2_maps.png')
fig1.savefig('fig_cov_l1l2_slices.png')
fig1.savefig('../../lensing_response/paper/figures/fig_cov_l1l2_slices.pdf')
fig2.savefig('fig_cov_l1l2_diagonal.png')
fig2.savefig('../../lensing_response/paper/figures/fig_cov_l1l2_diagonal.pdf')

#show()
