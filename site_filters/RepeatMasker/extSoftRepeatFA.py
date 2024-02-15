#!/usr/bin/python

import re
import textwrap
import sys

def ranges(nums):
    nums = sorted(set(nums))
    gaps = [[s, e] for s, e in zip(nums, nums[1:]) if s+1 < e]
    edges = iter(nums[:1] + sum(gaps, []) + nums[-1:])
    return list(zip(edges, edges))


infile = sys.argv[1]
with open(infile, 'r') as input:
    block= [] # save sequence
    chr = ""
    for line in input:
        if line.startswith('>'):
            if block:# deal with the previous sequence
                sequence = ''.join(block)
                positions = [i for i,a in enumerate(sequence) if a.islower()]
                #for pos in positions:
                   # print(chr+"\t"+str(pos)+"\t"+str(pos+1)+"\t"+str(sequence[pos]))
                for i in ranges(positions):
                    startend = list(i)
                    print( chr+"\t"+str(startend[0])+"\t"+str(startend[1]+1) )
                block = [] # reset sequence
                
            
            header = line.rstrip("\n")[1:]
            chr = header.partition(' ')[0]
            #print(chr)
        else:
            block.append(line.rstrip("\n"))

    if block: # the very last sequence
        sequence = ''.join(block)
        positions = [i for i,a in enumerate(sequence) if a.islower()]
        #for pos in positions:
        #    print(chr+"\t"+str(pos)+"\t"+str(pos+1))
        for i in ranges(positions):
            startend = list(i)
            print( chr+"\t"+str(startend[0])+"\t"+str(startend[1]+1) )
        block = [] # reset sequence         
