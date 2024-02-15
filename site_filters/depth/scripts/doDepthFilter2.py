# adpated from frederik ffs script. this one creates a bed list with sites to remove!!
# python3 doDepthFilter.py indir min thres max thres outprefix
#      indir: folder with input .pos.gz files
#      medianfile: file with a single number which is median depth
#      chromfile: file with list of chromosomes depth has been stimated for
#      outprefix: prefix of outpuf file. will create outprefix.bed and outprefix.log

import gzip
import sys
import glob
import os

def is_within_thresholds(x, min_thres, max_thres):
    global nout
    global nin
    
    if x == 'NA\n':
        nout += 1
        return False

    x = int(x)
    if x < min_thres or x > max_thres:
        nout += 1
        return False
    
    nin += 1
    return True


def get_bad_ranges(infile, out, min_thres, max_thres):

    is_out = False
    prev_chrr = None
    prev_pos = None
    with open(infile, 'rt') as f:
        for line in f:
            chrr, pos, d = line.split('\t')
            
            # when chromosome ends write to bed if is out
            if not chrr == prev_chrr and is_out:
                end = int(prev_pos)
                out.write(f'{prev_chrr}\t{start}\t{end}\n')
                start = int(pos) - 1 # update to new position
            prev_chrr = chrr
            prev_pos = pos

            if not is_out:
                if is_within_thresholds(d, min_thres, max_thres):
                    continue
                elif not is_within_thresholds(d, min_thres, max_thres):
                    start = int(pos) - 1
                    is_out = True
            elif is_out:
                if is_within_thresholds(d, min_thres, max_thres):
                    end = int(pos) - 1
                    out.write(f'{chrr}\t{start}\t{end}\n')
                    is_out = False

                elif not is_within_thresholds(d, min_thres, max_thres):
                    continue

        if is_out:
             end = int(pos)
             out.write(f'{chrr}\t{start}\t{end}\n')



infile = sys.argv[1] # sys.argv[0] is the name of the script!   
outpre = sys.argv[2]

minl = int(sys.argv[3])
maxl = int(sys.argv[4])

outfile = outpre + ".bed"
logfile = outpre + ".log"

out = open(outfile, 'w+')
log = open(logfile, 'w+')


log.write(f"# Will create bed list with regions to remove {outfile}.\n")
log.write(f"# Will only keep sites with depth between {minl} and {maxl}.\n\n")
log.write("#Summary of filter\n")
         
nout = 0
nin = 0

get_bad_ranges(infile, out, minl, maxl)
 

log.write(f"Total Kept: {nin}\nTotal Removed: {nout}\nTotal: {nin+nout}\nProportion Kept: {nin/(nin+nout)}\nProportion Removed: {nout/(nin+nout)}\n")


out.close()
log.close()

