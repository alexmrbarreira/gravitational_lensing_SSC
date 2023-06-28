import sys; sys.path.append('../../'); from prepare_for_lenscov import *
from params_here import *

## ======================================================================= 
## Load the data
## ======================================================================= 

data_vector = loadtxt('../../data_store/data_C_l1.dat')[:,1]
matrix_C1C2 = outer(data_vector, data_vector)

cov_ssc_fsky_0005 = loadtxt('data_ssc_cov_l1l2_fsky_0.005.dat')
cov_ssc_fsky_001  = loadtxt('data_ssc_cov_l1l2_fsky_0.01.dat')
cov_ssc_fsky_005  = loadtxt('data_ssc_cov_l1l2_fsky_0.05.dat')
cov_ssc_fsky_012  = loadtxt('data_ssc_cov_l1l2_fsky_0.12.dat')
cov_ssc_fsky_020  = loadtxt('data_ssc_cov_l1l2_fsky_0.20.dat')
cov_ssc_fsky_036  = loadtxt('data_ssc_cov_l1l2_fsky_0.36.dat')
cov_ssc_fsky_050  = loadtxt('data_ssc_cov_l1l2_fsky_0.50.dat')

cov_ssc_beyondLimber_fsky_0005 = loadtxt('data_ssc_cov_l1l2_beyondLimber_mono_fsky_0.005.dat')
cov_ssc_beyondLimber_fsky_001  = loadtxt('data_ssc_cov_l1l2_beyondLimber_mono_fsky_0.01.dat')
cov_ssc_beyondLimber_fsky_005  = loadtxt('data_ssc_cov_l1l2_beyondLimber_mono_fsky_0.05.dat')
cov_ssc_beyondLimber_fsky_012  = loadtxt('data_ssc_cov_l1l2_beyondLimber_mono_fsky_0.12.dat')
cov_ssc_beyondLimber_fsky_020  = loadtxt('data_ssc_cov_l1l2_beyondLimber_mono_fsky_0.20.dat')
cov_ssc_beyondLimber_fsky_036  = loadtxt('data_ssc_cov_l1l2_beyondLimber_mono_fsky_0.36.dat')
cov_ssc_beyondLimber_fsky_050  = loadtxt('data_ssc_cov_l1l2_beyondLimber_mono_fsky_0.50.dat')

#cov_ssc_beyondLimber_fsky_0005 = loadtxt('data_ssc_cov_l1l2_beyondLimber_mono_fsky_0.005_noderiv.dat')
#cov_ssc_beyondLimber_fsky_001  = loadtxt('data_ssc_cov_l1l2_beyondLimber_mono_fsky_0.01_noderiv.dat')
#cov_ssc_beyondLimber_fsky_005  = loadtxt('data_ssc_cov_l1l2_beyondLimber_mono_fsky_0.05_noderiv.dat')
#cov_ssc_beyondLimber_fsky_012  = loadtxt('data_ssc_cov_l1l2_beyondLimber_mono_fsky_0.12_noderiv.dat')
#cov_ssc_beyondLimber_fsky_020  = loadtxt('data_ssc_cov_l1l2_beyondLimber_mono_fsky_0.20_noderiv.dat')
#cov_ssc_beyondLimber_fsky_036  = loadtxt('data_ssc_cov_l1l2_beyondLimber_mono_fsky_0.36_noderiv.dat')
#cov_ssc_beyondLimber_fsky_050  = loadtxt('data_ssc_cov_l1l2_beyondLimber_mono_fsky_0.50_noderiv.dat')

## ========================================================
## Plot
## ========================================================
labelsize = 22
ticksize  = 22
textsize  = 30
titlesize = 26
text_font = 22
legend_font = 22
v_min = 0.95
v_max = 1.10
cbticks = [v_min, 1.0, 1.05, 1.10]

fig0 = plt.figure(0, figsize=(17.5, 10.))
fig0.subplots_adjust(left=0.075, right=0.98, top=0.96, bottom=-0.20, hspace = 0.15)#, space = 0.2)
fig0.add_subplot(2,2,1)


panel = fig0.add_subplot(2,2,1)
plot(l1_array, 1.0e4*diagonal(cov_ssc_beyondLimber_fsky_0005/matrix_C1C2), linestyle = 'dashed' , linewidth = 2, c = 'b', label = r'$f_{\rm sky} = 0.005$')
plot(l1_array, 1.0e4*diagonal(cov_ssc_beyondLimber_fsky_001 /matrix_C1C2) , linestyle = 'dashed' , linewidth = 2, c = 'g', label = r'$f_{\rm sky} = 0.01$')
plot(l1_array, 1.0e4*diagonal(cov_ssc_beyondLimber_fsky_005 /matrix_C1C2) , linestyle = 'dashed' , linewidth = 2, c = 'r', label = r'$f_{\rm sky} = 0.05$')
plot(l1_array, 1.0e4*diagonal(cov_ssc_beyondLimber_fsky_012 /matrix_C1C2) , linestyle = 'dashed' , linewidth = 2, c = 'c', label = r'$f_{\rm sky} = 0.12$')
plot(l1_array, 1.0e4*diagonal(cov_ssc_beyondLimber_fsky_020 /matrix_C1C2) , linestyle = 'dashed' , linewidth = 2, c = 'm', label = r'$f_{\rm sky} = 0.20$')
plot(l1_array, 1.0e4*diagonal(cov_ssc_beyondLimber_fsky_036 /matrix_C1C2) , linestyle = 'dashed' , linewidth = 2, c = 'k', label = r'$f_{\rm sky} = 0.36$')
#plot(l1_array, 1.0e4*diagonal(cov_ssc_beyondLimber_fsky_050 /matrix_C1C2), linestyle = 'dashed' , linewidth = 2, c = 'k', label = r'$f_{\rm sky} = 0.50$')
plot(l1_array, 1.0e4*diagonal(cov_ssc_fsky_0005/matrix_C1C2), linestyle = 'solid' , linewidth = 2, c = 'b')
plot(l1_array, 1.0e4*diagonal(cov_ssc_fsky_001 /matrix_C1C2), linestyle = 'solid' , linewidth = 2, c = 'g')
plot(l1_array, 1.0e4*diagonal(cov_ssc_fsky_005 /matrix_C1C2), linestyle = 'solid' , linewidth = 2, c = 'r')
plot(l1_array, 1.0e4*diagonal(cov_ssc_fsky_012 /matrix_C1C2), linestyle = 'solid' , linewidth = 2, c = 'c')
plot(l1_array, 1.0e4*diagonal(cov_ssc_fsky_020 /matrix_C1C2), linestyle = 'solid' , linewidth = 2, c = 'm')
plot(l1_array, 1.0e4*diagonal(cov_ssc_fsky_036 /matrix_C1C2), linestyle = 'solid' , linewidth = 2, c = 'k')
#plot(l1_array, 1.0e4*diagonal(cov_ssc_fsky_050 /matrix_C1C2), linestyle = 'solid' , linewidth = 2, c = 'k')
annotate(r'$Solid\ \ \ \ \   :\ Limber$', xy = (0.02,0.13), xycoords='axes fraction', color = 'k', fontsize = text_font)
annotate(r'$Dashed\        :\ Beyond\ Limber$', xy = (0.02,0.05), xycoords='axes fraction', color = 'k', fontsize = text_font)
#xlabel(r'$\ell$' , fontsize = labelsize)
ylabel(r'${\rm Cov}^{SSC}(\ell, \ell)/[C(\ell)]^2 \times 10^4$' , fontsize = labelsize)
xscale('log')
yscale('log')
xlim(min(l1_array), max(l1_array))
ylim(2.0e-2, 1.0e3)
xticks(size = ticksize)
yticks(size = ticksize)
params = {'legend.fontsize': legend_font-2}; pylab.rcParams.update(params); legend(loc = 'upper right', ncol = 2)

fig0.add_subplot(4,2,5)
plot(l1_array, diagonal(cov_ssc_beyondLimber_fsky_0005/cov_ssc_fsky_0005), linestyle = 'solid' , linewidth = 2, c = 'b')
plot(l1_array, diagonal(cov_ssc_beyondLimber_fsky_001 /cov_ssc_fsky_001) , linestyle = 'solid' , linewidth = 2, c = 'g')
plot(l1_array, diagonal(cov_ssc_beyondLimber_fsky_005 /cov_ssc_fsky_005) , linestyle = 'solid' , linewidth = 2, c = 'r')
plot(l1_array, diagonal(cov_ssc_beyondLimber_fsky_012 /cov_ssc_fsky_012) , linestyle = 'solid' , linewidth = 2, c = 'c')
plot(l1_array, diagonal(cov_ssc_beyondLimber_fsky_020 /cov_ssc_fsky_020) , linestyle = 'solid' , linewidth = 2, c = 'm')
plot(l1_array, diagonal(cov_ssc_beyondLimber_fsky_036 /cov_ssc_fsky_036) , linestyle = 'solid' , linewidth = 2, c = 'k')
#plot(l1_array, diagonal(cov_ssc_beyondLimber_fsky_050 /cov_ssc_fsky_050) , linestyle = 'solid' , linewidth = 2, c = 'k')
fill_between(l1_array, ones(len(l1_array))*1.01, ones(len(l1_array))*0.99, color = 'grey', alpha = 0.5)
xlabel(r'$\ell$' , fontsize = labelsize+2)
ylabel(r'$Beyond\ Limber\ /\ Limber$'        , fontsize = labelsize+2)
xscale('log')
xlim(min(l1_array), max(l1_array))
ylim(v_min, v_max)
xticks(size = ticksize)
yticks(size = ticksize)

ii = [3,4,7,8,11,12]

for i in ii:
    panel = fig0.add_subplot(4,4,i)

    if(i== 3): pcolor(l1_array, l2_array, cov_ssc_beyondLimber_fsky_036 /cov_ssc_fsky_036 , vmin = v_min, vmax = v_max)
    if(i== 4): pcolor(l1_array, l2_array, cov_ssc_beyondLimber_fsky_020 /cov_ssc_fsky_020 , vmin = v_min, vmax = v_max)
    if(i== 7): pcolor(l1_array, l2_array, cov_ssc_beyondLimber_fsky_012 /cov_ssc_fsky_012 , vmin = v_min, vmax = v_max)
    if(i== 8): pcolor(l1_array, l2_array, cov_ssc_beyondLimber_fsky_005 /cov_ssc_fsky_005 , vmin = v_min, vmax = v_max)
    if(i==11): pcolor(l1_array, l2_array, cov_ssc_beyondLimber_fsky_001 /cov_ssc_fsky_001 , vmin = v_min, vmax = v_max)
    if(i==12): pcolor(l1_array, l2_array, cov_ssc_beyondLimber_fsky_0005/cov_ssc_fsky_0005, vmin = v_min, vmax = v_max)

    if(i== 3):
        annotate(r"$f_{\rm sky} = 0.36$",xy = (0.06, 0.85), xycoords='axes fraction', xytext = (0.06, 0.85), textcoords='axes fraction', bbox=dict(boxstyle="round4", fc = 'w', alpha = 0.6), color = 'k', fontsize = textsize-6)
        title(r'$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ Beyond\ Limber\ /\ Limber$'   , fontsize = titlesize)
        ylabel(r'$\ell_2$' , fontsize = labelsize)
        panel.axes.get_xaxis().set_visible(False)
    if(i== 4):
        annotate(r"$f_{\rm sky} = 0.20$",xy = (0.06, 0.85), xycoords='axes fraction', xytext = (0.06, 0.85), textcoords='axes fraction', bbox=dict(boxstyle="round4", fc = 'w', alpha = 0.6), color = 'k', fontsize = textsize-6)
        panel.axes.get_xaxis().set_visible(False)
        panel.axes.get_yaxis().set_visible(False)
    if(i== 7):
        annotate(r"$f_{\rm sky} = 0.12$",xy = (0.06, 0.85), xycoords='axes fraction', xytext = (0.06, 0.85), textcoords='axes fraction', bbox=dict(boxstyle="round4", fc = 'w', alpha = 0.6), color = 'k', fontsize = textsize-6)
        ylabel(r'$\ell_2$' , fontsize = labelsize)
        panel.axes.get_xaxis().set_visible(False)
    if(i== 8):
        annotate(r"$f_{\rm sky} = 0.05$",xy = (0.06, 0.85), xycoords='axes fraction', xytext = (0.06, 0.85), textcoords='axes fraction', bbox=dict(boxstyle="round4", fc = 'w', alpha = 0.6), color = 'k', fontsize = textsize-6)
        panel.axes.get_xaxis().set_visible(False)
        panel.axes.get_yaxis().set_visible(False)
    if(i==11):
        annotate(r"$f_{\rm sky} = 0.01$",xy = (0.06, 0.85), xycoords='axes fraction', xytext = (0.06, 0.85), textcoords='axes fraction', bbox=dict(boxstyle="round4", fc = 'w', alpha = 0.6), color = 'k', fontsize = textsize-6)
        xlabel(r'$\ell_1$' , fontsize = labelsize)
        ylabel(r'$\ell_2$' , fontsize = labelsize)
    if(i==12):
        annotate(r"$f_{\rm sky} = 0.005$",xy = (0.06, 0.85), xycoords='axes fraction', xytext = (0.06, 0.85), textcoords='axes fraction', bbox=dict(boxstyle="round4", fc = 'w', alpha = 0.6), color = 'k', fontsize = textsize-6)
        xlabel(r'$\ell_1$' , fontsize = labelsize)
        panel.axes.get_yaxis().set_visible(False)

    panel.set_aspect('equal')
    
    colorbar(ticks=cbticks).ax.tick_params(labelsize=labelsize)

    xscale('log')
    yscale('log')
    xlim(min(l1_array), max(l1_array))
    ylim(min(l2_array), max(l2_array))
    xticks(size = ticksize)
    yticks(size = ticksize)

fig0.savefig('fig_Limber_vs_fky.png')

show()

