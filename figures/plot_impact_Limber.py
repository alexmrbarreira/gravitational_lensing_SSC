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

cov_ssc_beyondLimber_mono      = loadtxt('../data_store/data_ssc_cov_l1l2_beyondLimber_mono.dat')
cov_ssc_beyondLimber_mono_noRK = loadtxt('../data_store/data_ssc_cov_l1l2_beyondLimber_mono_noRK.dat')

cov_ssc     = cov_ssc_a + cov_ssc_b + cov_ssc_c + cov_ssc_d
cov_ssc_bcd = cov_ssc_b + cov_ssc_c + cov_ssc_d
cov_cng   = loadtxt('../data_store/data_cng_cov_l1l2.dat')
#cov_g     = loadtxt('../data_store/data_g_cov_l1l2.dat') 
cov_g     = loadtxt('../data_store/data_g_cov_l1l2_wnoise.dat') 
cov_tot   = cov_ssc + cov_cng + cov_g

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
pcolor(l1_array, l2_array, 1.0e4*cov_ssc/matrix_C1C2, vmin=v_min, vmax=v_max)
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
pcolor(l1_array, l2_array, 1.0e4*(cov_ssc_beyondLimber_mono)/matrix_C1C2, vmin=v_min, vmax=v_max)
panel.set_aspect('equal')
colorbar().ax.tick_params(labelsize=ticksize)
annotate(r"$beyond\ Limber$", xy = (0.04, 0.90), xycoords='axes fraction',
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
    plot(l1_array, 1.0e4*cov_ssc[:,l2index]                   /data_vector/data_vector[l2index], linestyle = 'solid' , linewidth = 2, c = 'b', label = r'${\rm Cov}^{\rm SSC}_{\kappa}$')
    plot(l1_array, 1.0e4*cov_ssc_beyondLimber_mono[:,l2index] /data_vector/data_vector[l2index], linestyle = 'solid' , linewidth = 2, c = 'r', label = r'$beyond\ Limber$')


    annotate(r'$\ell_2\ =\ ' +  str("%.0f" % l2_array[l2index]) + '$', xy = (0.76,0.86), xycoords='axes fraction', color = 'k', fontsize = text_font+4)
    ylabel(r'${\rm Cov}_{12}/(C_1C_2) \times 10^4$'        , fontsize = labelsize+2)
    xscale('log')
    xlim(min(l1_array), max(l1_array))
    #ylim(0.0, 0.80)
    xticks(size = ticksize+2)
    yticks(size = ticksize+2)
    if(i>1):
        xlabel(r'$\ell_1$' , fontsize = labelsize+2)
    if(i==3):
        params = {'legend.fontsize': legend_font+3}; pylab.rcParams.update(params); legend(loc = 'upper left', ncol = 2)

fig2 = plt.figure(2, figsize=(9., 6.))
fig2.subplots_adjust(left=0.17, right=0.98, top=0.94, bottom= 0.16, hspace = 0.15)#, space = 0.2)
fig2.add_subplot(1,1,1)

plot(l1_array, 1.0e4*diagonal(cov_ssc                   /matrix_C1C2), linestyle = 'solid' , linewidth = 2, c = 'b', label = r'${\rm Cov}^{\rm SSC}_{\kappa}$')
plot(l1_array, 1.0e4*diagonal(cov_ssc_beyondLimber_mono /matrix_C1C2), linestyle = 'solid' , linewidth = 2, c = 'r', label = r'$beyond\ Limber$')

plot(l1_array, 1.0e4*diagonal(cov_ssc_a                      /matrix_C1C2), linestyle = 'dashed' , linewidth = 2, c = 'b')
plot(l1_array, 1.0e4*diagonal(cov_ssc_beyondLimber_mono_noRK /matrix_C1C2), linestyle = 'dashed' , linewidth = 2, c = 'r')

annotate(r'$Solid\  :\ Full$'     , xy = (0.05,0.15), xycoords='axes fraction', color = 'k', fontsize = text_font+4)
annotate(r'$Dashed\ :\ Isotropic$', xy = (0.05,0.08), xycoords='axes fraction', color = 'k', fontsize = text_font+4)

xlabel(r'$\ell$' , fontsize = labelsize)
ylabel(r'${\rm Cov}(\ell, \ell)/[C(\ell)]^2 \times 10^4$'        , fontsize = labelsize+2)
xscale('log')
yscale('log')
xlim(min(l1_array), max(l1_array))
ylim(1.0e-3, 1.0e3)
xticks(size = ticksize+1)
yticks(size = ticksize+1)
params = {'legend.fontsize': legend_font-1}; pylab.rcParams.update(params); legend(loc = 'upper right', ncol = 3)

show()
