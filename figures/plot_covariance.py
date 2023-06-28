import sys; sys.path.append('../'); from prepare_for_lenscov import *

## ========================================================
## Prepare what to plot
## ========================================================

data_vector = loadtxt('../data_store/data_C_l1.dat')[:,1]
matrix_C1C2 = outer(data_vector, data_vector)

# Load Limber SSC
cov_ssc_a = loadtxt('../data_store/data_ssc_cov_l1l2_A.dat')
cov_ssc_b = loadtxt('../data_store/data_ssc_cov_l1l2_B.dat')
cov_ssc_c = loadtxt('../data_store/data_ssc_cov_l1l2_C.dat')
cov_ssc_d = loadtxt('../data_store/data_ssc_cov_l1l2_D.dat')

cov_ssc     = cov_ssc_a + cov_ssc_b + cov_ssc_c + cov_ssc_d
cov_ssc_bcd = cov_ssc_b + cov_ssc_c + cov_ssc_d

# Load beyondLimber SSC
cov_ssc_bLimber       = loadtxt('../data_store/data_ssc_cov_l1l2_beyondLimber_mono.dat')
cov_ssc_bLimber_noRK  = loadtxt('../data_store/data_ssc_cov_l1l2_beyondLimber_mono_noRK.dat')
cov_ssc_bLimber_tidal = cov_ssc_bLimber - cov_ssc_bLimber_noRK

# Load cNG
cov_cng   = loadtxt('../data_store/data_cng_cov_l1l2.dat')

# Load Gaussian
#cov_g     = loadtxt('../data_store/data_g_cov_l1l2.dat') 
cov_g     = loadtxt('../data_store/data_g_cov_l1l2_wnoise.dat') 

# Compute total
#cov_tot   = cov_g + cov_cng + cov_ssc
cov_tot   = cov_g + cov_cng + cov_ssc_bLimber



cov_ssc_bLimber       = cov_ssc_a + cov_ssc_b + cov_ssc_c + cov_ssc_d 
cov_ssc_bLimber_noRK  = cov_ssc_a
cov_ssc_bLimber_tidal = cov_ssc_bLimber - cov_ssc_bLimber_noRK


cov_tot   = cov_g + cov_cng + cov_ssc_bLimber


## ========================================================
## Compute Signal-to-Noise ratio as a function of ell 
## ========================================================
def compute_SN(n_tomo, n_lbin, l_predic, d_predic, covmat):
    n_data = len(d_predic)
    n_dvec = int(n_tomo*(n_tomo+1.)/2.)

    indices_to_delete = [] # will store indices to tell which colums/lines/entries to remove

    ellmax_array = zeros(n_lbin)
    sn_array     = zeros(n_lbin)
    # Loop over number of ell values
    for k in range(n_lbin):
        # Delete the wanted columns
        covmat_now = delete(covmat    , indices_to_delete, 1)
        covmat_now = delete(covmat_now, indices_to_delete, 0)

        l_predic_now   = delete(l_predic, indices_to_delete)
        d_predic_now   = delete(d_predic, indices_to_delete)

        # Invert the resulting reduced covariance matrix
        invcov_now = linalg.inv(covmat_now)

        # Compute SN and fill the arrays backwards
        ellmax_array[n_lbin-1-k] = max(l_predic_now)
        sn_array[n_lbin-1-k]     = sqrt(dot(d_predic_now, dot(invcov_now  , d_predic_now)))

        # Increment number of indices to delete
        indices_to_delete  = indices_to_delete + [i*n_lbin+(n_lbin-1)-k for i in range(n_dvec)]
    # That's it; return
    return ellmax_array, sn_array

ellmax, snratio_g           = compute_SN(1, len(l1_array), l1_array, data_vector, cov_g) # ellmax is the same for all so load only here; it is equal to l1_array anyway
snratio_g_cng               = compute_SN(1, len(l1_array), l1_array, data_vector, cov_g + cov_cng              )[1]
snratio_g_ssc_bLimber_noRK  = compute_SN(1, len(l1_array), l1_array, data_vector, cov_g + cov_ssc_bLimber_noRK )[1]
snratio_g_ssc_bLimber_tidal = compute_SN(1, len(l1_array), l1_array, data_vector, cov_g + cov_ssc_bLimber_tidal)[1]
snratio_g_ssc_bLimber       = compute_SN(1, len(l1_array), l1_array, data_vector, cov_g + cov_ssc_bLimber      )[1]
snratio_tot                 = compute_SN(1, len(l1_array), l1_array, data_vector, cov_tot                      )[1]

## ========================================================
## Plot
## ========================================================
labelsize = 26
ticksize  = 24
textsize  = 30
titlesize = 26
text_font = 22
legend_font = 22
v_min = 0.01
v_max = 0.60

fig0 = plt.figure(0, figsize=(17., 7.))
fig0.subplots_adjust(left=0.05, right=1.00, top=0.93, bottom=0.14, wspace = 0.15, hspace = 0.15)

panel = fig0.add_subplot(1,2,1)
#pcolor(l1_array, l2_array, 1.0e4*cov_ssc_beyondLimber_mono/matrix_C1C2)
pcolor(l1_array, l2_array, 1.0e4*cov_ssc_bLimber/matrix_C1C2, vmin=v_min, vmax=v_max)
panel.set_aspect('equal')
colorbar().ax.tick_params(labelsize=ticksize)
annotate(r"${\rm SSC}$", xy = (0.04, 0.90), xycoords='axes fraction',
xytext = (0.04, 0.90), textcoords='axes fraction', bbox=dict(boxstyle="round4", fc = 'w', alpha = 0.6), color = 'k', fontsize = textsize)
title(r'${\rm Cov}(\ell_1, \ell_2)/(C(\ell_1)C(\ell_2)) \times 10^4$'   , fontsize = titlesize)
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
colorbar().ax.tick_params(labelsize=ticksize)
annotate(r"${\rm G\ +\ cNG}$", xy = (0.04, 0.90), xycoords='axes fraction',
xytext = (0.04, 0.90), textcoords='axes fraction', bbox=dict(boxstyle="round4", fc = 'w', alpha = 0.6), color = 'k', fontsize = textsize)
title(r'${\rm Cov}(\ell_1, \ell_2)/(C(\ell_1)C(\ell_2)) \times 10^4$'   , fontsize = titlesize)
xlabel(r'$\ell_1$' , fontsize = labelsize)
ylabel(r'$\ell_2$' , fontsize = labelsize)
xscale('log')
yscale('log')
xlim(min(l1_array), max(l1_array))
ylim(min(l2_array), max(l2_array))
xticks(size = ticksize)
yticks(size = ticksize)

fig1 = plt.figure(1, figsize=(17.5, 10.))
fig1.subplots_adjust(left=0.07, right=0.99, top=0.96, bottom=0.10, wspace = 0.18, hspace = 0.15)
l2indices = [int(len(l1_array)*0.)+0, int(len(l1_array)*0.33), int(len(l1_array)*0.66), len(l1_array)-1] 
for i in range(4):
    panel = fig1.add_subplot(2,2,i+1)
    l2index = l2indices[i]
    plot(l1_array, 1.0e4*cov_ssc_bLimber_noRK[:,l2index]  /data_vector/data_vector[l2index], linestyle = 'solid' , linewidth = 2, c = 'g', label = r'${\rm Cov}^{\rm SSC}_{\kappa, Density}$')
    plot(l1_array, 1.0e4*cov_ssc_bLimber_tidal[:,l2index] /data_vector/data_vector[l2index], linestyle = 'dashed', linewidth = 2, c = 'g', label = r'${\rm Cov}^{\rm SSC}_{\kappa, Tidal}$')
    plot(l1_array, 1.0e4*cov_ssc_bLimber[:,l2index]       /data_vector/data_vector[l2index], linestyle = 'solid' , linewidth = 2, c = 'b', label = r'${\rm Cov}^{\rm SSC}_{\kappa}$')
    plot(l1_array, 1.0e4*cov_cng[:,l2index]               /data_vector/data_vector[l2index], linestyle = 'solid' , linewidth = 2, c = 'r', label = r'${\rm Cov}^{\rm cNG}_{\kappa}$')
    plot(l1_array, 1.0e4*cov_tot[:,l2index]               /data_vector/data_vector[l2index], linestyle = 'solid' , linewidth = 2, c = 'k', label = r'${\rm Cov}^{\rm Total}_{\kappa}$')
    annotate(r'$\ell_2\ =\ ' +  str("%.0f" % l2_array[l2index]) + '$', xy = (0.76,0.86), xycoords='axes fraction', color = 'k', fontsize = text_font+4)
    ylabel(r'${\rm Cov}_{12}/(C_1C_2) \times 10^4$'        , fontsize = labelsize+2)
    xscale('log')
    xlim(min(l1_array), max(l1_array))
    ylim(0.0, 0.80)
    xticks(size = ticksize+2)
    yticks(size = ticksize+2)
    if(i>1):
        xlabel(r'$\ell_1$' , fontsize = labelsize+2)
    if(i==3):
        params = {'legend.fontsize': legend_font+2}; pylab.rcParams.update(params); legend(loc = 'upper left', ncol = 2)

fig2 = plt.figure(2, figsize=(17.5, 10.))
fig2.subplots_adjust(left=0.09, right=0.98, top=0.95, bottom=-0.20, hspace = 0.15)#, space = 0.2)
fig2.add_subplot(2,2,1)
plot(l1_array, 1.0e4*diagonal(cov_ssc_bLimber_noRK /matrix_C1C2), linestyle = 'solid'  , linewidth = 2, c = 'g', label = r'${\rm Cov}^{\rm SSC}_{\kappa, Density}$')
plot(l1_array, 1.0e4*diagonal(cov_ssc_bLimber_tidal/matrix_C1C2), linestyle = 'dashed' , linewidth = 2, c = 'g', label = r'${\rm Cov}^{\rm SSC}_{\kappa, Tidal}$')
plot(l1_array, 1.0e4*diagonal(cov_ssc_bLimber      /matrix_C1C2), linestyle = 'solid'  , linewidth = 2, c = 'b', label = r'${\rm Cov}^{\rm SSC}_{\kappa}$')
plot(l1_array, 1.0e4*diagonal(cov_cng              /matrix_C1C2), linestyle = 'solid'  , linewidth = 2, c = 'r', label = r'${\rm Cov}^{\rm cNG}_{\kappa}$')
plot(l1_array, 1.0e4*diagonal(cov_g                /matrix_C1C2), linestyle = 'solid'  , linewidth = 2, c = 'c', label = r'${\rm Cov}^{\rm G}_{\kappa}$')
plot(l1_array, 1.0e4*diagonal(cov_tot              /matrix_C1C2), linestyle = 'solid'  , linewidth = 2, c = 'k', label = r'${\rm Cov}^{\rm Total}_{\kappa}$')
ylabel(r'${\rm Cov}(\ell, \ell)/[C(\ell)]^2 \times 10^4$'        , fontsize = labelsize+2)
xscale('log')
yscale('log')
xlim(min(l1_array), max(l1_array))
ylim(1.0e-3, 1.0e3)
xticks(size = ticksize+1)
yticks(size = ticksize+1)
params = {'legend.fontsize': legend_font-1.5}; pylab.rcParams.update(params); legend(loc = 'upper right', ncol = 3)
fig2.add_subplot(4,2,5)
plot(l1_array, diagonal(cov_g)                 / diagonal(cov_tot), linestyle = 'solid'  , linewidth = 2, c = 'c')
plot(l1_array, diagonal(cov_ssc_bLimber)       / diagonal(cov_tot), linestyle = 'solid'  , linewidth = 2, c = 'b')
plot(l1_array, diagonal(cov_ssc_bLimber_noRK)  / diagonal(cov_tot), linestyle = 'solid'  , linewidth = 2, c = 'g')
plot(l1_array, diagonal(cov_ssc_bLimber_tidal) / diagonal(cov_tot), linestyle = 'dashed' , linewidth = 2, c = 'g')
plot(l1_array, diagonal(cov_cng)               / diagonal(cov_tot), linestyle = 'solid'  , linewidth = 2, c = 'r')
axhline(1.0, linestyle = 'solid' , linewidth = 2, c = 'k',)
xlabel(r'$\ell$' , fontsize = labelsize+2)
ylabel(r'$Frac.\ contribution$'        , fontsize = labelsize+2)
xscale('log')
xlim(min(l1_array), max(l1_array))
ylim(0., 1.03)
xticks(size = ticksize)
yticks(size = ticksize)

fig2.add_subplot(2,2,2)

plot(ellmax, snratio_g                     , linestyle = 'dashed' , linewidth = 2, c = 'k'     , label = r'$G$')
plot(ellmax, snratio_g_cng                 , linestyle = 'solid'  , linewidth = 2, c = 'orange', label = r'$G+cNG$')
plot(ellmax, snratio_g_ssc_bLimber_tidal   , linestyle = 'dashed' , linewidth = 2, c = 'orange', label = r'$G+SSC_{Tidal}$')
plot(ellmax, snratio_g_ssc_bLimber_noRK    , linestyle = 'dashed' , linewidth = 2, c = 'm'     , label = r'$G+SSC_{Density}$')
plot(ellmax, snratio_g_ssc_bLimber         , linestyle = 'solid'  , linewidth = 2, c = 'm'     , label = r'$G+SSC$')
plot(ellmax, snratio_tot                   , linestyle = 'solid'  , linewidth = 2, c = 'k'     , label = r'$G+cNG+SSC$')

ylabel(r'$Signal-to-Noise(\ell<\ell_{\rm max})$'       , fontsize = labelsize+2)
xscale('log')
yscale('log')
xlim(min(ellmax), max(ellmax))
ylim(7.0e0, 1.0e3)
xticks(size = ticksize+1)
yticks(size = ticksize+1)
params = {'legend.fontsize': legend_font}; pylab.rcParams.update(params); legend(loc = 'lower right', ncol = 2)
fig2.add_subplot(4,2,6)
plot(ellmax, snratio_g                    / snratio_g - 1., linestyle = 'dashed', linewidth = 2, c = 'k')
plot(ellmax, snratio_g_cng                / snratio_g - 1., linestyle = 'solid' , linewidth = 2, c = 'orange')
plot(ellmax, snratio_g_ssc_bLimber_tidal  / snratio_g - 1., linestyle = 'dashed', linewidth = 2, c = 'orange')
plot(ellmax, snratio_g_ssc_bLimber_noRK   / snratio_g - 1., linestyle = 'dashed', linewidth = 2, c = 'm')
plot(ellmax, snratio_g_ssc_bLimber        / snratio_g - 1., linestyle = 'solid' , linewidth = 2, c = 'm')
plot(ellmax, snratio_tot                  / snratio_g - 1., linestyle = 'solid' , linewidth = 2, c = 'k')

axhline(0.0, linestyle = 'dashed' , linewidth = 1, c = 'k',)
xlabel(r'$\ell_{\rm max}$' , fontsize = labelsize+2)
ylabel(r'$Rel.\ diff.\ to\ G$'        , fontsize = labelsize+2)
xscale('log')
xlim(min(ellmax), max(ellmax))
ylim(-0.70, 0.03)
xticks(size = ticksize+1)
yticks(size = ticksize+1)

fig0.savefig('fig_cov_l1l2_maps.png')
fig1.savefig('fig_cov_l1l2_slices.png')
fig1.savefig('fig_cov_l1l2_slices.pdf')
fig2.savefig('fig_cov_l1l2_diagonal.png')
fig2.savefig('fig_cov_l1l2_diagonal.pdf')

show()
