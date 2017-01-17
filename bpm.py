"""Calculate between-population divergence metrics.  Input is a tab delimited file with no header, containing Scaff, pos, ac, an, dp, and genotypes coded as 0-4.
Input file custom: Filename should end with _XXX_raw.table.recode.txt, where XXX is three-letter population abbreviation"""

import argparse
import numpy

# export PYTHONPATH="$PYTHONPATH:/Users/monnahap/Documents/Research/code/GenomeScan/"

# Use 'sort <filename> -k3,3 -k4,4n > rho_test3.sort.txt' to sort by scaffold and then position (numerically)


# Substitute missingness for number of individuals to downsample to.  That way you can enter one number for all populations to be downsampled to.
def calcbpm(input_file1, input_file2, output, popnum, window_size, minimum_snps, minimum_individuals):

    def calcRholocus(locus_list):
        locus = []
        for l in locus_list:
            try:
                locus.append(l.remove("-9"))
            except ValueError:
                locus.append(l)
        min_ind = min([len(x[7:]) for x in locus])
        # Do rho calcs for site
        X_i = []
        sss = 0.0
        ssi = 0.0
        ssg = 0.0
        s1 = 0.0
        s2 = 0.0
        ni =[]
        ploidy_list = []
        for pop_site in locus:
            ploidy = int(pop_site[1])
            num_ind = len(pop_site[7:])  # THIS MUST BE CHANGED IF DOWNSAMPLING IS IMPLEMENTED
            
            x_i = sum([float(geno) for geno in pop_site[7:]]) / (num_ind * ploidy)
            ni.append(num_ind)
            X_i.append(x_i)
            ploidy_list.append(ploidy)
            s1 += num_ind
            s2 += num_ind**2
            for ind in pop_site[7:]:
                x_g = float(ind)/float(ploidy)
                ssi += ((x_g - x_i)**2)*float(ploidy)
        X_bar = numpy.mean(X_i)
        wa = s1 - (s2 / s1)
        wd = s1 - len(locus)
        ww = len(locus) - 1
        sst = 0.0
        for i, x in enumerate(X_i):
            sss += ((x - X_bar)**2) * ni[i] * ploidy_list[i]
        num = (wd * sss) - (ww * ssi)
        den = (wa * ssi)

        return num, den, jj


    def calcFst(locus_list):
        locus = []
        for l in locus_list:
            try:
                locus.append(l.remove("-9"))
            except ValueError:
                locus.append(l)
        r = float(len(locus))
        n_bar = float(sum([len(x[7:]) for x in locus]) / r)
        j = float(sum([len(x[7:])**2 for x in locus]))
        n_C = (r * n_bar - j / (r * n_bar)) / (r - 1.0)
        p_bar_i = []
        n_i = []
        h_i = []

        for pop_site in locus:
            ploidy = float(pop_site[1])
            nnn = float(len(pop_site[7:]))
            an = ploidy * nnn
            ac = sum([float(geno) for geno in pop_site[7:]])
            n_i.append(nnn)
            p_bar_i.append(ac / an)
            h = 0.0
            if ploidy == 4.0
                for ind in pop_site[7:]:
                    if float(ind) == 1.0 or float(ind) == 3.0:
                        h += 0.5
                    if float(ind) == 2.0:
                        h += 0.667
            if ploidy == 2.0
                for ind in pop_site[7:]
                    if float(ind) == 2.0:
                        h += 1.0

            h_i.append(h / nnn)

        p_bar = sum([j*n_i[i] for i,j in enumerate(p_bar_i)]) / (r * n_bar)
        h_bar = sum([j*n_i[i] for i,j in enumerate(p_bar_i)]) / (r * n_bar)

        var_p = 0.0
        for i,p in enumerate(p_bar_i):
            var_p += (n_i[i] * (p - p_bar)**2) / ((r - 1) * n_bar)

        a = (n_bar * n_C) * (var_p - ((1.0/(n_bar -1.0)) * ((p_bar * (1.0 - p_bar)) - (((r - 1.0) / r) * var_p) - (0.25 * h_bar))))
        b = (n_bar / (n_bar - 1.0)) * ((p_bar * (1.0 - p_bar)) - (((r - 1.0) / r) * var_p) - (((2.0 * (n_bar - 1.0)) / (4.0 * n_bar)) * h_bar))
        c = 0.5 * h_bar

        return a, b, c


    def calcDxy(locus_info):
        locus = []
        for l in locus_list:
            try:
                locus.append(l.remove("-9"))
            except ValueError:
                locus.append(l)

        for i, pop_site in enumerate(locus):
            ploidy = float(pop_site[1])
            nnn = float(len(pop_site[7:]))
            an = ploidy * nnn
            ac = sum([float(geno) for geno in pop_site[7:]])
            if i == 0:
                p1 = (ac / an)
            if i == 1:
                p2 = (ac / an)

        dxy = (p1 * (1.0 - p2)) + (p2 * (1.0 - p1))

        return dxy


    # Concatenate input files and sort them
    in1 = open(input_file1,'r').readlines()
    in2 = open(input_file2,'r').readlines()
    in1 = [j.split("\t") for j in in1]
    in2 = [j.split("\t") for j in in2]
    data = in1 + in2
    data = sorted(data, key=lambda k: (int(k[2].split("_")[1]), int(k[3]))) # Sorts by scaffold then position


    snp_count = 0
    Snp_count = 0
    start = 0.0
    end = window_size
    winexclcount = 0
    num_wind = 0
    infile = input_file.split("/")[-1]
    outfile = output + infile.replace('.txt','_RHO.txt')
    out1 = open(outfile, 'w')
    out1.write("pop1\tpop2\tscaff\tstart\tend\twin_size\tnum_snps\tRho\tFst\tdxy\n")
    Dxy = 0.0

    # To downsample or not to downsample?? 
    # Cat files prior running script and sort by scaffold and position
    # Does not handle missing data currently
    for i, line in enumerate(data):

        line = line.strip("\n")
        line = line.strip("\t")
        line = line.split("\t")

        pop, ploidy, scaff, pos, ac, an, dp = line[:7]
        pos = float(pos)
        gt = line[7:]

        if i % 100000 == 0:
            print(i)
        if i == 0:
            pop1 = pop
            old_pos = pos
            Locus = [] # Hold information of single locus across populations
            oldscaff = scaff
            Wd = [] # S1 - r, where r is number of subpops
            Ww = [] # r - 1
            SSs = []
            SSi = []
            Wa = [] # (S1 - S2)/S1, where S1 is sum(ni) and S2 is sum(ni^2) and ni is number individuals
            FNum = 0.0  # Numerator for entire genome
            FDen = 0.0  # Denominator for entire genome
            Fnum = 0.0  # Numerator for the window
            Fden = 0.0  # Denominator for the window
            RNum = 0.0  # Numerator for entire genome
            RDen = 0.0  # Denominator for entire genome
            Rnum = 0.0  # Numerator for the window
            Rden = 0.0  # Denominator for the window
            dxy = 0.0  # Dxy for window

        if pos > start and pos <= end and scaff == oldscaff:
            if pos == old_pos: # Accruing information from multiple populations but same locus
                Locus.append(line)
                if Snp_count == 0.0:
                    pop2 = pop
                    assert pop1 != pop2
            elif len(Locus) == popnum: # Within current window but have moved on from previous locus
                n, d, j = calclocus(Locus)
                a, b, c = calcFst(Locus)
                dxy += calcDxy(Locus)
                Dxy += dxy
                snp_count += j
                Fnum += a
                Fden += (a + b + c)
                FNum += a
                FDen += (a + b + c)
                Rnum += n
                Rden += d
                RNum += n
                RDen += d

                # RESET SUMS OF SQUARES
                Locus = []
                Locus.append(line)
                old_pos = pos
                
            else:
                Locus = []
                Locus.append(line)
                old_pos = pos


        elif int(pos) > end or scaff != oldscaff:
            if len(Locus) == popnum: # Current snp is onto next window, but must do calculation for previous locus before moving on
                snp_count += 1
                n, d, j = calclocus(Locus)
                a, b, c = calcFst(Locus)
                dxy += calcDxy(Locus)
                Dxy += dxy
                Fnum += a
                Fden += (a + b + c)
                FNum += a
                FDen += (a + b + c)
                Rnum += n
                Rden += d
                RNum += n
                RDen += d

            if snp_count >= minimum_snps:  # Use or exclude window from
                Snp_count += snp_count 
                num_wind += 1

                fst = Fnum / Fden

                fac = num / den

                rho = fac / (1 + fac)

                dxy = dxy / float(snp_count)

                out1.write(pop1 + '\t' + pop2 + '\t' + scaff + '\t' +
                           str(start) + '\t' +
                           str(end) + '\t' +
                           str(window_size) + '\t' +
                           str(snp_count) + '\t' +
                           str(rho) + '\t' +
                           str(fst) + '\t' +
                           str(dxy) + '\n')

            else:
                winexclcount += 1

            snp_count = 0
            p = []
            Fnum = 0.0
            Fden = 0.0
            Rnum = 0.0
            Rden = 0.0
            dxy = 0.0

            if float(pos) > end:
                while float(pos) > end:
                    end += window_size
                start = end - window_size
            elif scaff != oldscaff:
                oldscaff = scaff

                start = 0.0
                end = window_size

            if int(pos) > start and int(pos) <= end and scaff == oldscaff:
                Locus = []
                Locus.append(line)
                old_pos = pos

    if len(Locus) == popnum:  # Final window calculations
        snp_count += 1
        n, d = calcRholocus(Locus)
        a, b, c = calcFst(Locus)
        dxy += calcDxy(Locus)
        Dxy += dxy
        Fnum += a
        Fden += (a + b + c)
        FNum += a
        FDen += (a + b + c)
        Rnum += n
        Rden += d
        RNum += n
        RDen += d

    if snp_count >= minimum_snps:  # Use or exclude window from 
        num_wind += 1

        fac = num / den

        rho = fac / (1 + fac)

        dxy = dxy / snp_count

        out1.write(pop1 + '\t' + pop2 + '\t' + scaff + '\t' +
                       str(start) + '\t' +
                       str(end) + '\t' +
                       str(window_size) + '\t' +
                       str(snp_count) + '\t' +
                       str(rho) + '\t' +
                       str(fst) + '\t' +
                       str(dxy) + '\n')

    else:
        winexclcount += 1

    FAC = Num / Den
    rho = FAC / (1 + FAC)
    Fst = FNum / FDen
    out1.write(pop1 + '\t' + pop2 + '\t' + "Genome" + '\t' +
         "-99" + '\t' +
         "-99" + '\t' +
         str(window_size) + '\t' +
         str(Snp_count) + '\t' +
         str(rho) + '\t' +
         str(Fst) + '\t' +
         str(dxy) + '\n')

    out1.close()

    return num_wind, winexclcount, rho


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, metavar='input_file', required=True, help='input file created with recode012.py')
    parser.add_argument('-o', type=str, metavar='output_directory', required=True, help='Output Directory')
    parser.add_argument('-popnum', type=int, metavar='number_of_populations', required=True, help='')
    parser.add_argument('-sampind', type=int, metavar='DownSampled_individuals', required=False, default='5', help='Number of individuals to downsample the data to')
    parser.add_argument('-ws', type=float, metavar='window_size', required=False, default='10000.0', help='Size of windows in bp')
    parser.add_argument('-ms', type=int, metavar='minimum_snps', required=False, default='2', help='minimum number of snps in a window')

    args = parser.parse_args()

    j1, j2, j3 = calcrho(args.i, args.o, args.popnum, args.ws, args.ms, args.sampind)
