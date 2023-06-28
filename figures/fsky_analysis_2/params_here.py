import sys; sys.path.append('../../'); from prepare_for_lenscov import *;


#ffsky = 10**(linspace(log10(0.001), log10(0.5), 100))
#ffsky = 10**(linspace(log10(0.0001), log10(0.001), 10))
ffsky = 10**(linspace(log10(0.001), log10(1.), 100))

l1use_a = l1_array[8]
l2use_a = l1use_a

l1use_b = l1_array[12]
l2use_b = l1use_b

l1use_c = l1_array[16]
l2use_c = l1use_b

L_max_insigma = 15
L_insum       = array(range(2000))

epsrel_sigma  = 1.0e-2
epsrel_varint = 1.0e-3

nside = 1024
