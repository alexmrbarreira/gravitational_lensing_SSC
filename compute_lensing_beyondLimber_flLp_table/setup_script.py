import sys; sys.path.append('../'); from prepare_for_lenscov import *

print ''
print 'Preparing scripts for l1_array with size', len(l1_array)
print ''

#Delete prexisting files
os.system("rm -f compute_lensing_beyondLimber_flLp_table_index*")

for i in range(len(l1_array)):
    # Make copy of files
    os.system("cp compute_lensing_beyondLimber_flLp_table_masterfiletocopy.py compute_lensing_beyondLimber_flLp_table_index_"+str(i)+".py")
    # Adjust the index of ell each copy will deal with
    os.system(r"sed -i '/^index_l1/c\index_l1 = " + str(i) +  "' compute_lensing_beyondLimber_flLp_table_index_"+str(i)+".py")

