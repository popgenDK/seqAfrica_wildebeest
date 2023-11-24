##run_loter.py ceu.vcf yri.vcf mex.vcf test 30  ## test is only prefix for output
import sys
import allel
import numpy as np

print(sys.argv[1])
print(sys.argv[2])
print(sys.argv[3])
print(sys.argv[4])
print(int(sys.argv[5]))
r1=sys.argv[1]
r2=sys.argv[2]
mix=sys.argv[3]
output=sys.argv[4]
threads=int(sys.argv[5])

workdir='/home/users/xiaodong/Documents/Project/Wildebeest/loter'
#print(r1)
#print(threads)

def vcf2npy(vcfpath):
    callset = allel.read_vcf(vcfpath)
    haplotypes_1 = callset['calldata/GT'][:,:,0]
    haplotypes_2 = callset['calldata/GT'][:,:,1]

    m, n = haplotypes_1.shape
    mat_haplo = np.empty((2*n, m))
    mat_haplo[::2] = haplotypes_1.T
    mat_haplo[1::2] = haplotypes_2.T

    return mat_haplo.astype(np.uint8)


import os

popR1 = vcf2npy(os.path.join(workdir,'data',r1 ))
popR2 = vcf2npy(os.path.join(workdir,'data',r2 ))
popMix = vcf2npy(os.path.join(workdir,'data',mix ))

print ("Preparation done!")

import loter.locanc.local_ancestry as lc

res_loter = lc.loter_smooth(l_H=[popR1, popR2], h_adm=popMix, num_threads=threads) ## set the number of threads

import matplotlib.pyplot as plt

#plt.imshow(res_loter, interpolation='nearest', aspect='auto')
#plt.colorbar()
#plt.savefig(os.path.join(workdir,"results",output+'.png'))

np.savetxt(os.path.join(workdir,"results",output+"_anc.txt"), res_loter, fmt="%i")
