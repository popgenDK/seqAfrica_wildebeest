#!/usr/bin/python

import argparse
import sys
import allel
import numpy as np
import csv
import multiprocessing
from itertools import repeat

#print(allel.__version__)

def window_fst(chr,vcffile,pop1,pop2,windowsize,stepratio):
    stepsize = windowsize * stepratio
    callset = allel.read_vcf(vcffile,region=chr)
    pop1_list = [];
    pop2_list = [];
    
    for i,j in enumerate(callset['samples']):
        if j in pop1:
            pop1_list.append(i)
        elif j in pop2:
            pop2_list.append(i)


    subpops = [pop1_list,pop2_list]
    max_pos =  int(max(callset['variants/POS']))
    start = 1;
    end = start + windowsize -1
    output = list()
    while  (start+ windowsize) <= (max_pos+1) :
        index = np.argwhere( (callset['variants/POS'] < end ) & (callset['variants/POS'] >= start) )
        if (index.size == 0) :
            output.append([chr,start,end,start+windowsize/2,0,"NA"])


        else:
            gt = allel.GenotypeArray( callset['calldata/GT'][index.flatten().tolist()] )
            ac1 = gt.count_alleles(subpop=subpops[0])
            ac2 = gt.count_alleles(subpop=subpops[1])

            num, den = allel.hudson_fst(ac1, ac2)
            fst = np.sum(num) / np.sum(den)
            output.append( [chr,start,end,start+windowsize/2,index.size,fst] )
            #print(sca,"\t",start,"\t",end,"\t",start+windowsize/2,"\t",index.size,"\t",fst)

        start += stepsize
        end += stepsize

    return output

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', action='store', dest='vcffile', help="VCF file")
    parser.add_argument('-c', '--chr', action='store', dest='chrfile', help="CHR list")
    parser.add_argument('-p', '--pop', action='store', dest='popfile', help="Individual and population file")
    parser.add_argument('-p1','--pop1',action='store', dest='pop1', help="Population1")
    parser.add_argument('-p2','--pop2',action='store', dest='pop2', help="Population2")
    parser.add_argument('-w','--window',action='store', dest='windowsize', type=int,help="Size of sliding window")
    parser.add_argument('-f', '--factor', action='store', dest='factor', type=float,help="Proportion of step size to window size")
    parser.add_argument('-t', '--thread', action='store', dest='ncore', type=int, help="Number of threads to use")

    if len(sys.argv)==1:
         parser.print_help()
         parser.exit()
    
    args = parser.parse_args()
    #print(vars(args))
    
    vcffile = args.vcffile
    popfile = args.popfile
    chrfile = args.chrfile
    windowsize = args.windowsize
    factor = args.factor
    ncore = args.ncore
    p1 = args.pop1
    p2 = args.pop2
    
    with open(popfile, mode='r') as infile1:
        reader = csv.reader(infile1,delimiter="\t")
        mydict = {rows[0]:rows[1] for rows in reader}
    infile1.close()

    with open(chrfile,mode='r') as infile2:
        chrs = [line.rstrip('\n') for line in infile2]
    infile2.close()


    pop1 = [k for k,v in mydict.items() if v ==p1]
    pop2 = [k for k,v in mydict.items() if v ==p2]

    with multiprocessing.Pool(processes = ncore) as pool:
        result = pool.starmap(window_fst, zip(chrs,repeat(vcffile), repeat(pop1),repeat(pop2),repeat(windowsize),repeat(factor)))
    #result0 = window_fst("HiC_scaffold_29",vcffile,pop1,pop2,100000,0.5)

    for  i in range(len(chrs)) :
        for row in result[i]:
            print('\t'.join(map(str,row)))




