"""Calculate between-population divergence metrics.  Input is a tab delimited file with no header, containing Scaff, pos, ac, an, dp, and genotypes coded as 0-4.
Input file custom: Filename should end with _XXX_raw.table.recode.txt, where XXX is three-letter population abbreviation"""

import argparse
from random import randint
import numpy as np
import math
import itertools

# export PYTHONPATH="$PYTHONPATH:/Users/monnahap/Documents/Research/code/GenomeScan/"

# Add within population varP...could just use average SSI


def generateFSC2input(input_file, output, outname, numpops, window_size, minimum_snps, num_bootstraps):

    scaff_lengths = [33684748, 19642879, 24872290, 23717143, 21575646, 25532148, 25060017, 23333815]
    num_winds = sum([float(x) / float(window_size) for x in scaff_lengths])
    per_scaff = [float(x) / float(window_size) for x in scaff_lengths]
    wind_counts = []
    for j in range(0, num_bootstraps):
        wind_counts.append([0 for x in range(0, int(num_winds))])
        print(len(wind_counts[j]))
        for x in range(0, int(num_winds)):
            print(j, len(wind_counts[j]), randint(0, int(num_winds)), x)
            wind_counts[j][randint(0, int(num_winds) - 1)] += 1

    # Prepare output file
    outfile = output + outname + '_DSFS.obs'
    out = open(outfile, 'w')
    out.write("1 observations.  No. of demes and sample sizes are on next line\n")

    # Sort intput data
    data = open(input_file, 'r')
    data = [j.strip("\n").strip("\t").split("\t") for j in data]
    data = sorted(data, key=lambda k: (int(k[2].split("_")[1]), int(k[3]), k[0]))  # Sorts by scaffold then position, then population

    # Begin loop over data file
    snp_count = 0
    winexclcount = 0
    num_wind = 0
    get_pops = True
    num_inds = []
    ploidies = []
    num_alleles = []
    for i, line in enumerate(data):

        pop, ploidy, scaff, pos, ac, an, dp = line[:7]
        pos = float(pos)

        if i % 100000 == 0:
            print(i)
        if i == 0:
            old_pos = pos
            Locus = []  # Hold information of single locus across populations
            ac_i = []

        if get_pops is True and len(Locus) == numpops:
            string = str(numpops) + "\t\t\t"
            string1 = ""
            for pop_site in Locus:
                num_inds.append(len(pop_site[7:]))
                ploidies.append(int(pop_site[1]))
                num_alleles.append(len(pop_site[7:]) * int(pop_site[1]))
                string += str(len(pop_site[7:]) * int(pop_site[1])) + "\t"
                string1 += str(len(pop_site[7:]) * int(pop_site[1])) + ","
            string1.strip(",")
            string.strip("\t")
            string += "\n"
            get_pops = False
            for rep in range(0, num_bootstraps):
                exec("out%d = open('%s%s.rep%d_DSFS.obs', 'w')" % (rep + 1, output, outname, rep + 1))
                exec("DSFS%d = np.zeros(shape=(%s), dtype=np.int)" % (rep, string1))
                exec("out%d.write(%s)" % (rep + 1, string))

        if pos == old_pos:  # Accruing information from multiple populations but same locus
            Locus.append(line)
            try:
                line.remove("-9")
            except ValueError:
                ac_i.append(int(ac))
        elif len(Locus) == numpops:  # Skip locus calc if data not present from all populations
            snp_count += 1
            scaff_num = int(Locus[0][2].split("_")[1])
            cur_pos = float(Locus[0][3])
            wind_num = sum(per_scaff[:scaff_num - 1]) + int(math.ceil(cur_pos / float(window_size))) - 1
            for rep in wind_counts:
                index = ""
                for aa in ac_i:
                    index += '[' + aa + ']'
                for samp_num in range(0, rep[wind_num]):
                    exec("DSFS%d%s += 1" % (rep, index))
            ac_i = []
            Locus = []
            Locus.append(line)
            old_pos = pos

        else:  # Previous site contained data from only one population, so skip calculations
            ac_i = []
            Locus = []
            Locus.append(line)
            old_pos = pos

    ff = []
    for i, pop in enumerate(num_alleles):
        ff.append([jj for jj in range(0, pop + 1)])
    states = list(itertools.product(*ff))
    for state in states:
        dstring = ""
        for st in range(0, len(state)):
            dstring += "[" + state[st] + "]"
        for rep in range(0, num_bootstraps):
            exec("entry = DSFS%d%s" % (rep, dstring))
            exec("out%d.write(%s)" % (rep, str(entry) + "\t"))


    return num_wind, winexclcount


if __name__ == '__main__':  # Used to run code from command line

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, metavar='input_file', required=True, help='input file created with recode012.py')
    parser.add_argument('-o', type=str, metavar='output_directory', required=True, help='Output Directory')
    parser.add_argument('-np', type=int, metavar='number_of_populations', required=True)
    parser.add_argument('-prefix', type=str, metavar='output_file_prefix', required=True, help='Name indicating populations in input file')
    parser.add_argument('-ws', type=float, metavar='window_size', required=False, default='10000.0', help='Size of windows in bp')
    parser.add_argument('-ms', type=int, metavar='minimum_snps', required=False, default='2', help='minimum number of snps in a window')

    args = parser.parse_args()

    j1, j2, j3 = generateFSC2input(args.i, args.o, args.prefix, args.ws, args.ms)
