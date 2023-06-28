import sys; sys.path.append('../../'); from prepare_for_lenscov import *
from params_here import *


print ''
print ''
print 'Did not run! Uncomment quit() if you are sure you want to rewrite the compute_* files; they will have to be further adjusted by hand.'
print ''
print ''
quit()

os.system('rm -f compute_lensing_*')

for i in range(len(fsky_values)):
    os.system("cp ../../compute_lensing_SSC.py compute_lensing_SSC_fsky_"+str('%.2f'%fsky_values[i])+".py")

    os.system("cp ../../compute_lensing_beyondLimber_SSC.py compute_lensing_beyondLimber_SSC_fsky_"+str('%.2f'%fsky_values[i])+".py")
  
    os.system(r"sed -i '/^from/c\import sys; sys.path.append('../../'); from prepare_for_lenscov import *; from params_here import *; Omega_S = "+str(Omega_S_values[i])+"' compute_lensing_SSC_fsky_"+str('%.2f'%fsky_values[i])+".py")

    os.system(r"sed -i '/^from/c\import sys; sys.path.append('../../'); from prepare_for_lenscov import *; from params_here import *; Omega_S = "+str(Omega_S_values[i])+"' compute_lensing_beyondLimber_SSC_fsky_"+str('%.2f'%fsky_values[i])+".py")


