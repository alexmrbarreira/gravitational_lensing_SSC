import sys; sys.path.append('../'); from prepare_for_lenscov import *

i_to_submit = [4,5,6,7]

print ''
print 'Submitting ', i_to_submit, 'as background jobs ... '
print ''

for i in (i_to_submit):
    os.system(r"python compute_lensing_beyondLimber_flLp_table_index_"+str(i)+".py &")

