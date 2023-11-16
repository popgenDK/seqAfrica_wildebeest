//Parameters for the coalescence simulation program
3 sampels to simulate :
//Population effective size (number of genes)
NPOPB
NPOPT
NPOPW
//Haploid sample sizes
18
10
56
//Growth rates: negative growth implies population expansion
0
0
0
//Number of migration matrices: 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
3 historical event
TADM2        1  0       0.12    1       0       0
TDIV1        1  2       1       RES1    0       0
TDIV2        2  0       1       RES2    0       0
//Number of independent loci (chromosome)
1        0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ  1          0     1.45e-8   OUTEXP
