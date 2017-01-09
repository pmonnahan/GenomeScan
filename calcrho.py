"""Calculate within-population metrics for diversity and selection.  Input is a tab delimited file with no header, containing Scaff, pos, ac, an, dp, and genotypes coded as 0-4.
Input file custom: Filename should end with _XXX_raw.table.recode.txt, where XXX is three-letter population abbreviation"""

import argparse
import numpy

# export PYTHONPATH="$PYTHONPATH:/Users/monnahap/Documents/Research/code/GenomeScan/"

# concatenated, recoded files with  sort <filename> -k 3,4

# Substitute missingness for number of individuals to downsample to.  That way you can enter one number for all populations to be downsampled to.
def calcrho(input_file, output, popnum, sampind=5, window_size=50000, minimum_snps=2):

    snp_count = 0
    start = 0.0
    end = window_size
    winexclcount = 0
    num_wind = 0
    outfile = output + "_RHO.txt"
    out1 = open(outfile, 'w')
    out1.write("scaff\tstart\tend\twin_size\tnum_snps\tsampind\tRho\n")

    # To downsample or not to downsample?? 
    # Cat files prior running script and sort by scaffold and position
    with open(input_file, 'r') as infile:
        for i, line in enumerate(infile):

            line = line.strip("\n")
            line = line.strip("\t")
            line = line.split("\t")

            pop, ploidy, scaff, pos, ac, an, dp = line[:7]

            gt = line[7:]

            if i % 100000 == 0:
                print(i)
            if i == 0:
                old_pos = pos
                Locus = [] # Hold information of single locus across populations
                oldscaff = scaff
                AN = int(sampind * ploidy)
                n = float(AN)
                p = []
                Wd = [] # S1 - r, where r is number of subpops
                Ww = [] # r - 1
                SSs = []
                SSi = []
                Wa = [] # (S1 - S2)/S1, where S1 is sum(ni) and S2 is sum(ni^2) and ni is number individuals
                Num = 0.0
                Den = 0.0

            if int(pos) > start and int(pos) <= end and scaff == oldscaff:
                if pos == old_pos: # Accruing information from multiple populations but same locus
                    Locus.append(line)
                elif len(Locus) == popnum: # Within current window but have moved on from previous locus
                    if snp_count == 0:
                        max_ind = max([len(x) for x in Locus])
                    X_i = []
                    sss = 0.0
                    ssi = 0.0
                    num = 0.0
                    den = 0.0
                    s1 = 0.0
                    s2 = 0.0
                    wd = 0.0
                    X_bar = 0.0
                    Tot_alleles = 0.0
                    for pop_site in Locus:
                        ploidy = pop_site[1]
                        num_ind = len(pop_site[7:])  # THIS MUST BE CHANGED IF DOWNSAMPLING IS IMPLEMENTED
                        x_i = sum([int(geno) for geno in pop_site[7:]]) / (num_ind * ploidy)
                        X_i.append(x_i)
                        s1 += num_ind
                        s2 += num_ind**2
                        for ind in pop_site[7:]:
                            ssi += (float(ind) / float(ploidy)) - x_i
                            Tot_alleles += float(ploidy)
                            X_bar += float(ind)
                        # Do rho calcs for site
                    X_bar = X_bar / Tot_alleles
                    wa = (s1 - s2) / s1
                    wd = s1 - len(Locus)
                    ww = len(Locus) - 1
                    for x in X_i:
                        sss += x - X_bar

                    num += (wd * sss) - (ww * ssi)
                    den += (wa * ssi)
                    Num += (wd * sss) - (ww * ssi)
                    Den += (wa * ssi)

                    # RESET SUMS OF SQUARES
                    Locus = []
                    Locus.append(line)
                    old_pos = pos
                    snp_count += 1
                    # sgt = numpy.random.choice(gt, size=sampind, replace=False)
                    # sac = sum([int(x) for x in sgt])
                    p1 = float(ac) / float(AN)
                    p.append(p1)
                else:
                    Locus = []
                    Locus.append(line)
                    old_pos = pos
                    snp_count += 1
                    # sgt = numpy.random.choice(gt, size=sampind, replace=False)
                    # sac = sum([int(x) for x in sgt])
                    p1 = float(ac) / float(AN)
                    p.append(p1)


            elif int(pos) > end or scaff != oldscaff:

                n = float(AN)

                if snp_count >= minimum_snps:  # Use or exclude window from 
                    num_wind += 1

                    fac = num / den

                    rho = fac / (1 + fac)

                    out1.write(scaff + '\t' +
                               str(start) + '\t' +
                               str(end) + '\t' +
                               str(window_size) + '\t' +
                               str(snp_count) + '\t' +
                               str(sampind) + '\t' +
                               str(numpy.mean(p)) + '\t' +
                               str(rho))

                else:
                    winexclcount += 1

                snp_count = 0
                p = []
                num = 0.0
                den = 0.0

                if float(pos) > end:
                    while float(pos) > end:
                        end += window_size
                    start = end - window_size
                elif scaff != oldscaff:
                    oldscaff = scaff

                    start = 0.0
                    end = window_size

                if int(pos) > start and int(pos) <= end and int(an) >= AN and scaff == oldscaff:
                    snp_count += 1
                    sgt = numpy.random.choice(gt, size=sampind, replace=False)
                    sac = sum([int(x) for x in sgt])
                    p1 = float(sac) / float(AN)
                    p.append(p1)


        if snp_count >= minimum_snps: #accomodate final window
            n = float(AN)
            num_wind += 1

            out1.write(scaff + '\t' +
               str(start) + '\t' +
               str(end) + '\t' +
               str(window_size) + '\t' +
               str(snp_count) + '\t' +
               str(sampind) + '\t' +
               str(numpy.mean(p)) + '\t' +
               str(rho))

    n = float(AN)

    out1.write(scaff + '\t' +
         str(start) + '\t' +
         str(end) + '\t' +
         str(window_size) + '\t' +
         str(snp_count) + '\t' +
         str(sampind) + '\t' +
         str(numpy.mean(p)) + '\t' +
         str(rho))

    out1.close()

    return num_wind, winexclcount, rho


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, metavar='input_file', required=True, help='input file created with recode012.py')
    parser.add_argument('-o', type=str, metavar='output_directory', required=True, help='Output Directory')
    parser.add_argument('-popnum', type=str, metavar='number_of_populations', required=True, help='')
    parser.add_argument('-sampind', type=int, metavar='DownSampled_individuals', required=False, default='5', help='Number of individuals to downsample the data to')
    parser.add_argument('-ws', type=float, metavar='window_size', required=False, default='10000.0', help='Size of windows in bp')
    parser.add_argument('-ms', type=int, metavar='minimum_snps', required=False, default='2', help='minimum number of snps in a window')

    args = parser.parse_args()

    j1, j2, j3 = calcrho(args.i, args.o, args.sampind, args.ws, args.ms)
