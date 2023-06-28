import sys; sys.path.append('../'); 
from prepare_for_lenscov import *
from matplotlib.collections import LineCollection

## ========================================================
## Prepare what to plot
## ========================================================
zzhere = linspace(0., 1.8, 100)
kkhere = 10.**linspace(-2, 1, 500)

R1_store = [] 
RK_store = [] 
for i in range(len(zzhere)):
    R1_tmp = zeros(len(kkhere))
    RK_tmp = zeros(len(kkhere))
    for j in range(len(kkhere)):
        R1_tmp[j] = R_1_int(zzhere[i], kkhere[j])
        RK_tmp[j] = R_K_int(zzhere[i], kkhere[j])
    R1_store.append(R1_tmp)
    RK_store.append(RK_tmp)

line_segments_R1 = LineCollection([list(zip(kkhere,y)) for y in R1_store], linewidths= 2.0, linestyles = 'solid')
line_segments_RK = LineCollection([list(zip(kkhere,y)) for y in RK_store], linewidths= 2.0, linestyles = 'solid')

line_segments_R1.set_array(zzhere)
line_segments_RK.set_array(zzhere)

## ========================================================
## Plot parameters
## ========================================================
labelsize = 26
ticksize  = 22
textsize  = 20
titlesize = 24
text_font = 24

## ========================================================
## Plot total response functions
## ========================================================
fig0 = plt.figure(0, figsize=(17., 6.))
fig0.subplots_adjust(left=0.065, right=1.0, top=0.96, bottom=0.16, wspace = 0.15, hspace = 0.15)

# R_1(k,z)
panel = fig0.add_subplot(1,2,1)
panel.add_collection(line_segments_R1)
colorbar(line_segments_R1).ax.tick_params(labelsize=ticksize)
axvline(pi*512./500., linestyle = 'dashed', c = 'k')
annotate(r"$Colors\ indicate\ redshift$", xy = (0.05,0.80), xycoords='axes fraction', color = 'k', fontsize = text_font)
xlabel(r'$k\ \left[h/{\rm Mpc}\right]$'        , fontsize = labelsize)
ylabel(r'$R_1(k, z)$'        , fontsize = labelsize)
xscale('log')
xlim(min(kkhere), max(kkhere))
ylim(0.5, 4.0)
xticks(size = ticksize)
yticks(size = ticksize)

# R_K(k,z)
panel = fig0.add_subplot(1,2,2)
panel.add_collection(line_segments_RK)
colorbar(line_segments_RK).ax.tick_params(labelsize=ticksize)
axvline(2.0, linestyle = 'dashed', c = 'k')
annotate(r"$Colors\ indicate\ redshift$", xy = (0.05,0.80), xycoords='axes fraction', color = 'k', fontsize = text_font)
xlabel(r'$k\ \left[h/{\rm Mpc}\right]$'        , fontsize = labelsize)
ylabel(r'$R_K(k, z)$'        , fontsize = labelsize)
xscale('log')
xlim(min(kkhere), max(kkhere))
ylim(0.5, 4.0)
xticks(size = ticksize)
yticks(size = ticksize)

fig0.savefig('fig_R1RK.png')

show()
