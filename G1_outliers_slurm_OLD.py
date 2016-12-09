## genome scan pipeline, written by Jeff DaCosta, streamlined by Christian Sailer
## September 2016, updated 11 November 2016

import os, sys, subprocess, statistics, argparse
from natsort import natsorted
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser(description='Start of genome scan pipeline. This script combines the two allele counte tables of one contrast, '+
                                 'removes fixed non-variant sites and calculates Dxy, Fst, pi and allele frequency difference.')

parser.add_argument('-i', type=str, metavar='inputdir_path', required=True, help='REQUIRED: Full or relative path to input AC table directory')
parser.add_argument('-coh1', type=str, metavar='cohort_1', required=True, help='REQUIRED: Name of first cohort')
parser.add_argument('-coh2', type=str, metavar='cohort_2', required=True, help='REQUIRED: Name of second cohort')
parser.add_argument('-snps', type=int, metavar='snps_per_window', default='100', help='Number of SNPs per window [100]')
parser.add_argument('-o', type=str, metavar='outputdir_path', default='na', help='Optional output directory relative path, contrast directory will be a subdir of this one')
parser.add_argument('-ci', type=int, metavar='stdev_for_CI', default='4', help='Number of standard deviation for confidence interval: 3 for 99percent and 4 for 99.9percent CI [4]')
parser.add_argument('-per', type=float, metavar='percent_missing_data', default='0.25', help='Ratio of missing data allowed, eg 0.2 allows 80percent missing data [0.25]')

args = parser.parse_args()

contrast = args.i
contrast = contrast.split("/")[-2]

if args.o is 'na':
    outputdir = "output/"
else:
    outputdir = str(args.o+"output/")
    if os.path.exists(args.o) == False:
        os.mkdir(args.o)

# check if output directory exists, create it if necessary
if os.path.exists(outputdir) == False:
    os.mkdir(outputdir)
    print('\nCreated directory '+outputdir)

# before repeating SNP window calculation, test if that has been done and only the ci has to be changed
if os.path.isfile(outputdir+'/'+contrast+'_WG_'+str(args.snps)+'SNPs_3metrics.txt') == False:
    if os.path.isfile(outputdir+'/'+contrast+'_scaf_8_AFs.table') == False:
        

        ###### STEP 1 ######
        ## Prepare input data for analysis
        print('\nSearching input directory for *_raw.table')
        incoh1 = []
        incoh2 = []
        for dirName, subdirList, fileList in os.walk(args.i, topdown = False):
            for fname in fileList:
                if fname.endswith('_raw.table') and args.coh1 in fname:
                    incoh1.append(dirName+'/'+fname)
                elif fname.endswith('_raw.table') and args.coh2 in fname:
                    incoh2.append(dirName+'/'+fname)
        incoh1=natsorted(incoh1)
        incoh2=natsorted(incoh2)
                                    

        # paste cohort 1 & cohort 2 into a table next to each other & remove CHROM POS from second cohort in joint table
        # yields: CHROM POS AC AN DP AC AN DP
        print('\n\tSTEP 1: Paste contrast AC tables\n')
        for i in range(len(incoh1)):
            print('Processing '+incoh1[i])
            pastecmd = open('paste_'+args.coh1+'_'+args.coh2+'.unix', 'w')
            pastecmd.write('paste ')
            pastecmd.write(incoh1[i]+' '+incoh2[i]+' | cut -f -5,8,9,10 ')
            pastecmd.write('> '+outputdir+'/'+args.coh1+args.coh2+'_scaf_'+str(i+1)+'_temp.table')
            pastecmd.close()

            # run in unix
            cmd = (open('paste_'+args.coh1+'_'+args.coh2+'.unix', 'r'))
            p = subprocess.Popen(cmd, shell=True)
            sts = os.waitpid(p.pid, 0)[1]

        os.remove('paste_'+args.coh1+'_'+args.coh2+'.unix') # removes unix script file from folder


        ###### STEP 2 ######
        ## This section is Jeff's FixedDerivedAlleleCheck.py
        # search directory for output of above step to make new input list
        intable = []
        for dirName, subdirList, fileList in os.walk(outputdir):
            for table in fileList:
                if table.startswith(args.coh1+args.coh2) and table.endswith('_temp.table'):
                    intable.append(dirName+'/'+table)
        intable=natsorted(intable)            

        print('\n\tSTEP 2: Remove sites that are fixed in both cohorts and are missing in one but not the other cohort\n')
        print('Found '+str(len(intable))+' input tables')

        # obtain ANmax for each cohort
        print('\nGathering allele count info:\n')
        for i in range(len(intable)):
            infile = open(intable[i],'r')
            infile.readline()
            tANmax_1 = []
            tANmax_2 = []
            for line in infile:
                data = line.split()
                tANmax_1.append(int(data[3]))
                tANmax_2.append(int(data[6]))
            ANmax_1 = max(tANmax_1)
            ANmax_2 = max(tANmax_2)
            print('For '+intable[i]+': ')
            print('ANmax = '+str(ANmax_1)+' in cohort '+args.coh1)
            print('ANmax = '+str(ANmax_2)+' in cohort '+args.coh2)
            infile.close()

            # filter
            infile = open(intable[i],'r')
            outfile = open(outputdir+'/'+args.coh1+args.coh2+'_scaf_'+str(i+1)+'.tab','w')
            header = infile.readline()
            outfile.write(header)

            count = 0
            target = 100000
            Zero_count = 0
            All_count = 0

            for line in infile:
                data = line.split()
                if int(data[3])>=int(round((1-args.per)*ANmax_1,0)) and int(data[6])>=int(round((1-args.per)*ANmax_2,0)):
                    total_allele_count = int(data[2])+int(data[5])
                    total_alleles_sampled = int(data[3])+int(data[6])
                    if total_allele_count == total_alleles_sampled:
                        All_count +=1
                    elif total_allele_count == 0:
                        Zero_count +=1
                    else:
                        outfile.write(line)
                    count += 1
            infile.close()
            outfile.close()

            print('Total variants in file: '+str(count))
            print('Variants with total allele count = 0: '+str(Zero_count))
            print('Variants fixed in both cohorts: '+str(All_count))
            print('Retained '+str(count-Zero_count-All_count)+' variants\n')
            os.remove(outputdir+'/'+args.coh1+args.coh2+'_scaf_'+str(i+1)+'_temp.table') # removes output files of part 1, as they are no longer needed


        ###### STEP 3 ######
        ## jump highdepthFilter and remove DP column, Jeff's script adjusted
        print('\n\tSTEP 3: Remove DP column')
        intable = []
        for dirName, subdirList, fileList in os.walk(outputdir):
            for table in fileList:
                if table.startswith(args.coh1+args.coh2) and table.endswith('.tab'):
                    intable.append(dirName+'/'+table)
        intable=natsorted(intable)

        for i in range(len(intable)):
            infile = open(intable[i],'r')
            infile.readline()
            outfile = open(outputdir+'/'+contrast+'_scaf_'+str(i+1)+'_AFs.table','w')
            outfile.write('CHROM\tPOS\tAC1\tAN1\tAC2\tAN2\n')

            for line in infile:
                data = line.split()
                outfile.write(data[0]+'\t'+data[1]+'\t'+data[2]+'\t'+data[3]+'\t'+data[5]+'\t'+data[6]+'\n')
            infile.close()
            outfile.close()
            os.remove(outputdir+'/'+args.coh1+args.coh2+'_scaf_'+str(i+1)+'.tab') # removes output files of part 2, as they are no longer needed


        ###### STEP 4 ######
        # calculate Dxy, Fst, pi, allele frequency and allele frequncy difference
        # this is adjusted from Jeff's AFs_2_DD_Fst_Dxy.py script and should be able to handle missing data
        # description='This script calculates differentiation (AFD), diversity (pi), and Fst statistics across non-overlapping SNP windows in an
        # allele frequency table. The input file should have a header and the following columns:
        # scaffold, position, cohort1_allele_count, cohort1_all_alleles_count, cohort2_allele_count, cohort2_all_alleles_count

        print('\n\tSTEP 4: Calculate Dxy, Fst, DD')
        header = 'scaffold\tstart\tend\tmidpoint\tlength\t'+args.coh1+'_freq\t'+args.coh2+'_freq\tabsdiff\t'+args.coh1+'-'+args.coh2+'\t'+args.coh2+'-'+args.coh1+'\t'+args.coh1+'_pi\t'+args.coh2+'_pi\tFst\tDxy\n'
        outfile = open(outputdir+'/'+contrast+'_metrics_WG_'+str(args.snps)+'SNPs.txt','w')
        outfile.write(header)

        AFs = []
        for dirName, subdirList, fileList in os.walk(outputdir):
            for file in fileList:
                if file.startswith(contrast) and file.endswith('_AFs.table') == True:
                    AFs.append(dirName+'/'+file)
        AFs_sort = natsorted(AFs)

        print('\nFound '+str(len(AFs_sort))+' files starting with '+contrast+' and ending with _AFs.table\n')

        file_count = 0
        win_count = 0
        winexclcount = 0

        for AF in AFs_sort:
            count = -1 # to remove header from count
            infile = open(AF,'r')
            for line in infile:
                count += 1
            infile.close()
            num_win = int(count/args.snps)

            print('Processing '+str(num_win)+' windows in '+AF)

            if num_win > 0:
                file_count += 1
                win_count += num_win

                infile = open(AF,'r')
                infile.readline() # to read the header and start reading the numbers for calcluation from the second line of the file
                for i in range(num_win):
                    scaf = []
                    pos = []
                    c1_freq = []
                    c2_freq = []
                    absdiff = []
                    c1_minus_c2 = []
                    c2_minus_c1 = []
                    c1_pi = []
                    c2_pi = []
                    fst = []
                    dxy = []
                    for i in range(args.snps):
                        data = infile.readline()
                        data = data.split()
                        scaf.append(data[0])
                        pos.append(int(data[1]))
                        allele_c1_count = int(data[2])
                        allele_c2_count = int(data[4])
                        allele_c1_freq = allele_c1_count/int(data[3])
                        allele_c2_freq = allele_c2_count/int(data[5])
                        allele_count = allele_c1_count + allele_c2_count
                        allele_freq = allele_count/(int(data[3]) + int(data[5]))
                        c1_freq.append(allele_c1_freq)
                        c2_freq.append(allele_c2_freq)
                        c1_minus_c2.append(allele_c1_freq - allele_c2_freq)
                        c2_minus_c1.append(allele_c2_freq - allele_c1_freq)
                        absdiff.append(max([allele_c1_freq-allele_c2_freq, allele_c2_freq-allele_c1_freq]))
                        c1_pi.append((2*allele_c1_freq*(1-allele_c1_freq)))    
                        c2_pi.append((2*allele_c2_freq*(1-allele_c2_freq)))
                        ht = 2*allele_freq*(1-allele_freq)
                        h1 = 2*allele_c1_freq*(1-allele_c1_freq)
                        h2 = 2*allele_c2_freq*(1-allele_c2_freq)
                        h12 = ((h1*int(data[3]))+(h2*int(data[5])))/(int(data[3]) + int(data[5]))
                        allele_fst = abs(ht-h12)/ht
                        fst.append(allele_fst)
                        dxy.append((allele_c1_freq*(1-allele_c2_freq))+(allele_c2_freq*(1-allele_c1_freq)))

                    wstart = min(pos)
                    wend = max(pos)
                    wlength = wend-wstart
                    wmid = wstart+(wlength/2)
                    mean_c1_freq = statistics.mean(c1_freq)
                    mean_c2_freq = statistics.mean(c2_freq)
                    mean_absdiff = statistics.mean(absdiff)
                    mean_c1_minus_c2 = statistics.mean(c1_minus_c2)
                    mean_c2_minus_c1 = statistics.mean(c2_minus_c1)
                    mean_c1_pi = statistics.mean(c1_pi)
                    mean_c2_pi = statistics.mean(c2_pi)
                    mean_fst = statistics.mean(fst)
                    sum_dxy = sum(dxy)
                    win_dxy = (1/args.snps)*sum_dxy

                    if wlength <= 26560:
                        outfile.write(scaf[0]+'\t'+
                                      str(wstart)+'\t'+
                                      str(wend)+'\t'+
                                      str(wmid)+'\t'+
                                      str(wlength)+'\t'+
                                      str(mean_c1_freq)+'\t'+
                                      str(mean_c2_freq)+'\t'+
                                      str(mean_absdiff)+'\t'+
                                      str(mean_c1_minus_c2)+'\t'+
                                      str(mean_c2_minus_c1)+'\t'+
                                      str(mean_c1_pi)+'\t'+
                                      str(mean_c2_pi)+'\t'+
                                      str(mean_fst)+'\t'+
                                      str(win_dxy)+'\n')
                    else:
                        winexclcount +=1

        outfile.close()

        print('\nAnalyzed '+str(win_count)+' windows in '+str(file_count)+' files')
        print('Excluded '+str(winexclcount)+' windows longer than 26560bp\n')

        ###### STEP 4b ######
        # Add DD residuals using R script
        inlist = []
        for dirName, subdirList, fileList in os.walk(outputdir):
            for file in fileList:
                if file.startswith(contrast) and file.endswith(str(args.snps)+'SNPs.txt'):
                    inlist.append(dirName+'/'+file)
        print('Add DDresiduals')

        count=0
        for file in inlist:

            snp_file=pd.read_table(file,header=0)

            slope, intercept, r_value, p_value, std_err = stats.linregress(snp_file['absdiff'],snp_file[args.coh2+'_pi'])

            snp_file['prediction'] = intercept + (slope*snp_file['absdiff'])

            snp_file['DD'] = snp_file[args.coh2+'_pi'] - snp_file['prediction']

            slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(snp_file['absdiff'],snp_file[args.coh1+'_pi'])

            print(slope, intercept)
            print(slope1,intercept1)
            snp_file['prediction1'] = intercept1 + (slope1*snp_file['absdiff'])

            snp_file['DD1'] = snp_file[args.coh1+'_pi'] - snp_file['prediction1']

            snp_file.to_csv(outputdir+'/'+contrast+'_WG_'+str(args.snps)+'SNPs_3metrics.txt',sep="\t",index=False)

            count += 1

        print('Processed '+str(count)+' files for DDresiduals')
        os.remove(outputdir+'/'+contrast+'_metrics_WG_'+str(args.snps)+'SNPs.txt')


        ###### STEP 5 ######
        # define outliers using confidence intervall assuming a normal distribution
        print('\n\tSTEP 5: Define and select outliers')
        out = []
        for dirName, subdirList, fileList in os.walk(outputdir):
            for file in fileList:
                if file.startswith(contrast) and file.endswith(str(args.snps)+'SNPs_3metrics.txt') == True:
                    out.append(dirName+'/'+file)
        out_sort = natsorted(out)

        for out in out_sort:
            infile = open(out, 'r')
            header = infile.readline()
            dxy = []
            fst = []
            dd = []
            length = []
            
            for line in infile:
                data = line.split()
                dxy.append(float(data[13]))
                fst.append(float(data[12]))
                dd.append(float(data[15]))
                length.append(int(data[4]))
            mean_Dxy = statistics.mean(dxy)
            median_Dxy = statistics.median(dxy)
            sd_Dxy = statistics.stdev(dxy)
            CV_Dxy = sd_Dxy/mean_Dxy
            out_dxy = mean_Dxy+args.ci*sd_Dxy
            mean_Fst = statistics.mean(fst)
            median_Fst = statistics.median(fst)
            sd_Fst = statistics.stdev(fst)
            CV_Fst = sd_Fst/mean_Fst
            out_fst = mean_Fst+args.ci*sd_Fst
            mean_DD = statistics.mean(dd)
            median_DD = statistics.median(dd)
            sd_DD = statistics.stdev(dd)
            CV_DD = sd_DD/mean_DD
            out_dd = mean_DD-args.ci*sd_DD
            mean_win = statistics.mean(length)
            median_win = statistics.median(length)
            sd_win = statistics.stdev(length)
            CV_win = sd_win/mean_win
            out_win = mean_win+args.ci*sd_win
            # for graphs
            out3_dxy = mean_Dxy+3*sd_Dxy
            out4_dxy = mean_Dxy+4*sd_Dxy
            out3_fst = mean_Fst+3*sd_Fst
            out4_fst = mean_Fst+4*sd_Fst
            out3_dd = mean_DD-3*sd_DD
            out4_dd = mean_DD-4*sd_DD
            
            
            outlier_values = open(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_outlier_descriptive_stats.txt', 'w')
            outlier_values.write('stat\tDxy\tFst\tDD\twin_length\n'+
                                 'min\t'+str(round(min(dxy),4))+'\t'+str(round(min(fst),4))+'\t'+str(round(min(dd),4))+'\t'+str(round(min(length),0))+'\n'+
                                 'mean\t'+str(round(mean_Dxy,4))+'\t'+str(round(mean_Fst,4))+'\t'+str(round(mean_DD,4))+'\t'+str(round(mean_win,0))+'\n'+
                                 'median\t'+str(round(median_Dxy,4))+'\t'+str(round(median_Fst,4))+'\t'+str(round(median_DD,4))+'\t'+str(round(median_win,0))+'\n'+
                                 'sd\t'+str(round(sd_Dxy,4))+'\t'+str(round(sd_Fst,4))+'\t'+str(round(sd_DD,4))+'\t'+str(round(sd_win,0))+'\n'+
                                 'CV\t'+str(round(CV_Dxy,4))+'\t'+str(round(CV_Fst,4))+'\t'+str(round(CV_DD,4))+'\t'+str(round(CV_win,0))+'\n'+
                                 'max\t'+str(round(max(dxy),4))+'\t'+str(round(max(fst),4))+'\t'+str(round(max(dd),4))+'\t'+str(round(max(length),0))+'\n'+
                                 'CI\t'+str(round(out_dxy,4))+'\t'+str(round(out_fst,4))+'\t'+str(round(out_dd,4))+'\t'+str(round(out_win,0))+'\n')
            
            infile.close()

            # use pandas to read infiles as dataframe and select outlier combinations
            df_dat = pd.read_table(out)
            df_dat['dxyout'] = np.where(df_dat.Dxy >= out_dxy, 1, 0)
            df_dat['fstout'] = np.where(df_dat.Fst >= out_fst, 1, 0)
            df_dat['ddout'] = np.where(df_dat.DD <= out_dd, 1, 0)
            df_dat['dxyfst'] = df_dat.dxyout+df_dat.fstout
            df_dat['dxydd'] = df_dat.dxyout+df_dat.ddout
            df_dat['fstdd'] = df_dat.fstout+df_dat.ddout
            df_dat['dxyfstdd'] = df_dat.dxyout+df_dat.fstout+df_dat.ddout
            df_dat['bedstart'] = df_dat.start-1
            print('\nFor confidence intervals of '+str(args.ci)+' standard deviations we find:\n')
            print(str(sum(1 for x in df_dat.dxyout if x == 1))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4)))
            print(str(sum(1 for x in df_dat.fstout if x == 1))+'\toutlier windows for Fst values of >='+str(round(out_fst, 4)))
            print(str(sum(1 for x in df_dat.ddout if x == 1))+'\toutlier windows for DD values of <='+str(round(out_dd, 4)))
            print(str(sum(1 for x in df_dat.dxyfst if x == 2))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4))+' and Fst values of >='+str(round(out_fst, 4)))
            print(str(sum(1 for x in df_dat.dxydd if x == 2))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4))+' and DD values of <='+str(round(out_dd, 4)))
            print(str(sum(1 for x in df_dat.fstdd if x == 2))+'\toutlier windows for Fst values of >='+str(round(out_fst, 4))+' and DD values of <='+str(round(out_dd, 4)))
            print(str(sum(1 for x in df_dat.dxyfstdd if x == 3))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4))+', Fst values of >='+str(round(out_fst, 4))+' and DD values of <='+str(round(out_dd, 4))+'\n')
            print('Select oulier windows\n')
            
            df_dat.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_allsites.csv', index=False)
            # select all windows that are outliers for at least one metric
            df_outlier = df_dat[((df_dat.dxyout != 0) | (df_dat.fstout != 0) | (df_dat.ddout != 0))]
            df_outlier.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_ALL_outliers.csv', index=False)
            # write bedfile for gene retrieval and orthologous gene onthology match
            header = ["scaffold", "bedstart", "end"]
            df_outlier.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_ALL_outliers.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

            # select Dxy only outliers only
            df_dxy = df_dat[((df_dat.dxyout == 1) & (df_dat.fstout == 0) & (df_dat.ddout == 0))]
            if df_dxy.empty:
                print('No single Dxy outlier windows')
            else:
                df_dxy.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_Dxy_outliers.csv', index=False)
                # write bedfile for gene retrieval and orthologous gene onthology match
                header = ["scaffold", "bedstart", "end"]
                df_dxy.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_Dxy_outliers.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

            # select Fst only outliers only
            df_fst = df_dat[((df_dat.dxyout == 0) & (df_dat.fstout == 1) & (df_dat.ddout == 0))]
            if df_fst.empty:
                print('No single Fst outlier windows')
            else:
                df_fst.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_Fst_outliers.csv', index=False)
                # write bedfile for gene retrieval and orthologous gene onthology match
                header = ["scaffold", "bedstart", "end"]
                df_fst.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_Fst_outliers.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)
                
            # select DD only outliers only
            df_dd = df_dat[((df_dat.dxyout == 0) & (df_dat.fstout == 0) & (df_dat.ddout == 1))]
            if df_dd.empty:
                print('No single DD outlier windows')
            else:
                df_dd.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_DD_outliers.csv', index=False)
                # write bedfile for gene retrieval and orthologous gene onthology match
                header = ["scaffold", "bedstart", "end"]
                df_dd.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_DD_outliers.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

            # select Dxy and Fst double outliers only
            df_dxyfst = df_dat[((df_dat.dxyfst == 2) & (df_dat.dxyfstdd != 3))]
            if df_dxyfst.empty:
                print('No double Dxy Fst outlier windows')
            else:
                df_dxyfst.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_DxyFst_outliers.csv', index=False)
                # write bedfile for gene retrieval and orthologous gene onthology match
                header = ["scaffold", "bedstart", "end"]
                df_dxyfst.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_DxyFst_outliers.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

            # select Dxy and DD double outliers only
            df_dxydd = df_dat[((df_dat.dxydd == 2) & (df_dat.dxyfstdd != 3))]
            if df_dxydd.empty:
                print('No double Dxy DD outlier windows')
            else:
                df_dxydd.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_DxyDD_outliers.csv', index=False)
                # write bedfile for gene retrieval and orthologous gene onthology match
                header = ["scaffold", "bedstart", "end"]
                df_dxydd.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_DxyDD_outliers.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

            # select Fst and DD double outliers only
            df_fstdd = df_dat[((df_dat.fstdd == 2) & (df_dat.dxyfstdd != 3))]
            if df_fstdd.empty:
                print('No Fst DD double outlier windows')
            else:
                df_fstdd.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_FstDD_outliers.csv', index=False)
                # write bedfile for gene retrieval and orthologous gene onthology match
                header = ["scaffold", "bedstart", "end"]
                df_fstdd.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_FstDD_outliers.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

            # select Dxy, Fst, DD triple outliers only
            df_dxyfstdd = df_dat[(df_dat.dxyfstdd == 3)]
            if df_dxyfstdd.empty:
                print('No Dxy Fst DD triple outlier windows')
            else:
                df_dxyfstdd.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_DxyFstDD_outliers.csv', index=False)
                # write bedfile for gene retrieval and orthologous gene onthology match
                header = ["scaffold", "bedstart", "end"]
                df_dxyfstdd.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_DxyFstDD_outliers.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)


        ###### STEP 6 ######
        ## create histograms with the CI cutoff value for Dxy, Fst, DD and window length
        print('\n\tSTEP 6: Create population genetics histograms')
        graphs = []
        for dirName, subdirList, fileList in os.walk(outputdir):
            for file in fileList:
                if file.startswith(contrast) and file.endswith(str(args.snps)+'SNPs_'+str(args.ci)+'sd_allsites.csv') == True:
                    graphs.append(dirName+'/'+file)

        print('Creating histograms for Dxy, Fst, DD and window length')

        count=1
        num_bins = 50
        for file in graphs:

            data=pd.read_csv(file,header=0)
            f=plt.figure(1)
            plt.subplot(311)
            plt.hist(data['Fst'], num_bins, facecolor='green', alpha=0.5)
            plt.xlabel('Fst')
            plt.ylabel('Frequency')
            plt.subplot(312)
            plt.hist(data['Dxy'], num_bins, facecolor='blue', alpha=0.5)
            plt.xlabel('Dxy')
            plt.ylabel('Frequency')
            plt.subplot(313)
            plt.hist(data['DD'], num_bins, facecolor='red', alpha=0.5)
            plt.xlabel('DD')
            plt.ylabel('Frequency')

            plt.tight_layout()

            plt.show()

            f.savefig(outputdir+'histograms.pdf')

    
            # the histogram of the data

            rfile = open(outputdir+'graphs/'+contrast+'_graphs.r', 'w')
            rfile.write('# load table\n'+
                        'test <- read.csv("'+file+'", header=T)\n'+
                        '#str(test)\n'+
                        '#load library\n'+
                        'library(ggplot2)\n'+
                        'library(methods)\n'+
                        'p.dxy <- qplot(Dxy, data=test, binwidth=0.001, xlim=c(0, (max(test$Dxy)+0.01))) + geom_vline(xintercept='+str(out3_dxy)+', color="grey30", linetype="dashed") + geom_vline(xintercept='+str(out4_dxy)+', color="grey30", linetype="dotted") + theme_bw()\n'+
                        'p.fst <- qplot(Fst, data=test, binwidth=0.001, xlim=c(0, (max(test$Fst)+0.01))) + geom_vline(xintercept='+str(out3_fst)+', color="grey30", linetype=2) + geom_vline(xintercept='+str(out4_fst)+', color="grey30", linetype="dotted") + theme_bw()\n'+
                        'p.dd <- qplot(DD, data=test, binwidth=0.001) + geom_vline(xintercept='+str(out3_dd)+', color="grey30", linetype=2) + geom_vline(xintercept='+str(out4_dd)+', color="grey30", linetype="dotted") + theme_bw()\n'+
                        'p.length <- qplot(length, data=test, binwidth=100, xlim=c(0,26560)) + theme_bw()\n\n'+
                        '# export as pdf\n'+
                        'pdf(file="'+outputdir+'graphs/'+contrast+'_'+str(args.snps)+'SNPs_Dxy_histogram.pdf")\n'+
                        'p.dxy\n'+
                        'dev.off()\n'+
                        'pdf(file="'+outputdir+'graphs/'+contrast+'_'+str(args.snps)+'SNPs_Fst_histogram.pdf")\n'+
                        'p.fst\n'+
                        'dev.off()\n'+
                        'pdf(file="'+outputdir+'graphs/'+contrast+'_'+str(args.snps)+'SNPs_DD_histogram.pdf")\n'+
                        'p.dd\n'+
                        'dev.off()\n'+
                        'pdf(file="'+outputdir+'graphs/'+contrast+'_'+str(args.snps)+'SNPs_length_histogram.pdf")\n'+
                        'p.length\n'+
                        'dev.off()')
            rfile.close()

            cmd = ('Rscript '+outputdir+'/graphs/'+contrast+'_graphs.r')
            p = subprocess.Popen(cmd, shell=True)
            sts = os.waitpid(p.pid, 0)[1]

            count += 1
            #os.remove(outputdir+'/graphs/'+args.coh1+args.coh2+'_graphs.r')
        print('\n\tDONE\n')    

    else:
        print('\n\t SNPs have already been filtered for '+str(args.coh1)+str(args.coh2)+':')
        print('-> jump to')
        ###### STEP 4 ######
        print('\tSTEP 4: Calculate Dxy, Fst, DD')
        header = 'scaffold\tstart\tend\tmidpoint\tlength\t'+args.coh1+'_freq\t'+args.coh2+'_freq\tabsdiff\t'+args.coh1+'-'+args.coh2+'\t'+args.coh2+'-'+args.coh1+'\t'+args.coh1+'_pi\t'+args.coh2+'_pi\tFst\tDxy\n'
        outfile = open(outputdir+'/'+contrast+'_metrics_WG_'+str(args.snps)+'SNPs.txt','w')
        outfile.write(header)

        AFs = []
        for dirName, subdirList, fileList in os.walk(outputdir):
            for file in fileList:
                if file.startswith(contrast) and file.endswith('_AFs.table') == True:
                    AFs.append(dirName+'/'+file)
        AFs_sort = natsorted(AFs)

        print('\nFound '+str(len(AFs_sort))+' files starting with '+contrast+' and ending with _AFs.table\n')

        file_count = 0
        win_count = 0
        winexclcount = 0

        for AF in AFs_sort:
            count = -1 # to remove header from count
            infile = open(AF,'r')
            for line in infile:
                count += 1
            infile.close()
            num_win = int(count/args.snps)

            print('Processing '+str(num_win)+' windows in '+AF)

            if num_win > 0:
                file_count += 1
                win_count += num_win

                infile = open(AF,'r')
                infile.readline() # to read the header and start reading the numbers for calcluation from the second line of the file
                for i in range(num_win):
                    scaf = []
                    pos = []
                    c1_freq = []
                    c2_freq = []
                    absdiff = []
                    c1_minus_c2 = []
                    c2_minus_c1 = []
                    c1_pi = []
                    c2_pi = []
                    fst = []
                    dxy = []
                    for i in range(args.snps):
                        data = infile.readline()
                        data = data.split()
                        scaf.append(data[0])
                        pos.append(int(data[1]))
                        allele_c1_count = int(data[2])
                        allele_c2_count = int(data[4])
                        allele_c1_freq = allele_c1_count/int(data[3])
                        allele_c2_freq = allele_c2_count/int(data[5])
                        allele_count = allele_c1_count + allele_c2_count
                        allele_freq = allele_count/(int(data[3]) + int(data[5]))
                        c1_freq.append(allele_c1_freq)
                        c2_freq.append(allele_c2_freq)
                        c1_minus_c2.append(allele_c1_freq - allele_c2_freq)
                        c2_minus_c1.append(allele_c2_freq - allele_c1_freq)
                        absdiff.append(max([allele_c1_freq-allele_c2_freq, allele_c2_freq-allele_c1_freq]))
                        c1_pi.append((2*allele_c1_freq*(1-allele_c1_freq)))    
                        c2_pi.append((2*allele_c2_freq*(1-allele_c2_freq)))
                        ht = 2*allele_freq*(1-allele_freq)
                        h1 = 2*allele_c1_freq*(1-allele_c1_freq)
                        h2 = 2*allele_c2_freq*(1-allele_c2_freq)
                        h12 = ((h1*int(data[3]))+(h2*int(data[5])))/(int(data[3]) + int(data[5]))
                        allele_fst = abs(ht-h12)/ht
                        fst.append(allele_fst)
                        dxy.append((allele_c1_freq*(1-allele_c2_freq))+(allele_c2_freq*(1-allele_c1_freq)))

                    wstart = min(pos)
                    wend = max(pos)
                    wlength = wend-wstart
                    wmid = wstart+(wlength/2)
                    mean_c1_freq = statistics.mean(c1_freq)
                    mean_c2_freq = statistics.mean(c2_freq)
                    mean_absdiff = statistics.mean(absdiff)
                    mean_c1_minus_c2 = statistics.mean(c1_minus_c2)
                    mean_c2_minus_c1 = statistics.mean(c2_minus_c1)
                    mean_c1_pi = statistics.mean(c1_pi)
                    mean_c2_pi = statistics.mean(c2_pi)
                    mean_fst = statistics.mean(fst)
                    sum_dxy = sum(dxy)
                    win_dxy = (1/args.snps)*sum_dxy

                    if wlength <= 26560:
                        outfile.write(scaf[0]+'\t'+
                                      str(wstart)+'\t'+
                                      str(wend)+'\t'+
                                      str(wmid)+'\t'+
                                      str(wlength)+'\t'+
                                      str(mean_c1_freq)+'\t'+
                                      str(mean_c2_freq)+'\t'+
                                      str(mean_absdiff)+'\t'+
                                      str(mean_c1_minus_c2)+'\t'+
                                      str(mean_c2_minus_c1)+'\t'+
                                      str(mean_c1_pi)+'\t'+
                                      str(mean_c2_pi)+'\t'+
                                      str(mean_fst)+'\t'+
                                      str(win_dxy)+'\n')
                    else:
                        winexclcount +=1

        outfile.close()

        print('\nAnalyzed '+str(win_count)+' windows in '+str(file_count)+' files')
        print('Excluded '+str(winexclcount)+' windows longer than 26560bp\n')

        ###### STEP 4b ######
        # Add DD residuals using R script
        inlist = []
        for dirName, subdirList, fileList in os.walk(outputdir):
            for file in fileList:
                if file.startswith(contrast) and file.endswith(str(args.snps)+'SNPs.txt'):
                    inlist.append(dirName+'/'+file)
        print('Add DDresiduals')

        count = 0
        for file in inlist:
            
            snp_file=pd.read_table(file,header=0)

            slope, intercept, r_value, p_value, std_err = stats.linregress(snp_file['absdiff'],snp_file[args.coh2+'_pi'])

            snp_file['prediction'] = intercept + (slope*snp_file['absdiff'])

            snp_file['DD'] = snp_file[args.coh2+'_pi'] - snp_file['prediction']

            snp_file.to_csv(outputdir+'/'+contrast+'_WG_'+str(args.snps)+'SNPs_3metrics.txt',sep="\t",index=False)

            count += 1

        print('Processed '+str(count)+' files for DDresiduals')
        os.remove(outputdir+'/'+contrast+'_metrics_WG_'+str(args.snps)+'SNPs.txt')


        ###### STEP 5 ######
        # define outliers using confidence intervall assuming a normal distribution
        print('\n\tSTEP 5: Define and select outliers')
        out = []
        for dirName, subdirList, fileList in os.walk(outputdir):
            for file in fileList:
                if file.startswith(contrast) and file.endswith(str(args.snps)+'SNPs_3metrics.txt') == True:
                    out.append(dirName+'/'+file)
        out_sort = natsorted(out)

        for out in out_sort:
            infile = open(out, 'r')
            header = infile.readline()
            dxy = []
            fst = []
            dd = []
            length = []
            
            for line in infile:
                data = line.split()
                dxy.append(float(data[13]))
                fst.append(float(data[12]))
                dd.append(float(data[15]))
                length.append(int(data[4]))
            mean_Dxy = statistics.mean(dxy)
            median_Dxy = statistics.median(dxy)
            sd_Dxy = statistics.stdev(dxy)
            CV_Dxy = sd_Dxy/mean_Dxy
            out_dxy = mean_Dxy+args.ci*sd_Dxy
            mean_Fst = statistics.mean(fst)
            median_Fst = statistics.median(fst)
            sd_Fst = statistics.stdev(fst)
            CV_Fst = sd_Fst/mean_Fst
            out_fst = mean_Fst+args.ci*sd_Fst
            mean_DD = statistics.mean(dd)
            median_DD = statistics.median(dd)
            sd_DD = statistics.stdev(dd)
            CV_DD = sd_DD/mean_DD
            out_dd = mean_DD-args.ci*sd_DD
            mean_win = statistics.mean(length)
            median_win = statistics.median(length)
            sd_win = statistics.stdev(length)
            CV_win = sd_win/mean_win
            out_win = mean_win+args.ci*sd_win
            # for graphs
            out3_dxy = mean_Dxy+3*sd_Dxy
            out4_dxy = mean_Dxy+4*sd_Dxy
            out3_fst = mean_Fst+3*sd_Fst
            out4_fst = mean_Fst+4*sd_Fst
            out3_dd = mean_DD-3*sd_DD
            out4_dd = mean_DD-4*sd_DD
            
            
            outlier_values = open(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_outlier_descriptive_stats.txt', 'w')
            outlier_values.write('stat\tDxy\tFst\tDD\twin_length\n'+
                                 'min\t'+str(round(min(dxy),4))+'\t'+str(round(min(fst),4))+'\t'+str(round(min(dd),4))+'\t'+str(round(min(length),0))+'\n'+
                                 'mean\t'+str(round(mean_Dxy,4))+'\t'+str(round(mean_Fst,4))+'\t'+str(round(mean_DD,4))+'\t'+str(round(mean_win,0))+'\n'+
                                 'median\t'+str(round(median_Dxy,4))+'\t'+str(round(median_Fst,4))+'\t'+str(round(median_DD,4))+'\t'+str(round(median_win,0))+'\n'+
                                 'sd\t'+str(round(sd_Dxy,4))+'\t'+str(round(sd_Fst,4))+'\t'+str(round(sd_DD,4))+'\t'+str(round(sd_win,0))+'\n'+
                                 'CV\t'+str(round(CV_Dxy,4))+'\t'+str(round(CV_Fst,4))+'\t'+str(round(CV_DD,4))+'\t'+str(round(CV_win,0))+'\n'+
                                 'max\t'+str(round(max(dxy),4))+'\t'+str(round(max(fst),4))+'\t'+str(round(max(dd),4))+'\t'+str(round(max(length),0))+'\n'+
                                 'CI\t'+str(round(out_dxy,4))+'\t'+str(round(out_fst,4))+'\t'+str(round(out_dd,4))+'\t'+str(round(out_win,0))+'\n')
            
            infile.close()

            # use pandas to read infiles as dataframe and select outlier combinations
            df_dat = pd.read_table(out)
            df_dat['dxyout'] = np.where(df_dat.Dxy >= out_dxy, 1, 0)
            df_dat['fstout'] = np.where(df_dat.Fst >= out_fst, 1, 0)
            df_dat['ddout'] = np.where(df_dat.DD <= out_dd, 1, 0)
            df_dat['dxyfst'] = df_dat.dxyout+df_dat.fstout
            df_dat['dxydd'] = df_dat.dxyout+df_dat.ddout
            df_dat['fstdd'] = df_dat.fstout+df_dat.ddout
            df_dat['dxyfstdd'] = df_dat.dxyout+df_dat.fstout+df_dat.ddout
            df_dat['bedstart'] = df_dat.start-1
            print('\nFor confidence intervals of '+str(args.ci)+' standard deviations we find:\n')
            print(str(sum(1 for x in df_dat.dxyout if x == 1))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4)))
            print(str(sum(1 for x in df_dat.fstout if x == 1))+'\toutlier windows for Fst values of >='+str(round(out_fst, 4)))
            print(str(sum(1 for x in df_dat.ddout if x == 1))+'\toutlier windows for DD values of <='+str(round(out_dd, 4)))
            print(str(sum(1 for x in df_dat.dxyfst if x == 2))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4))+' and Fst values of >='+str(round(out_fst, 4)))
            print(str(sum(1 for x in df_dat.dxydd if x == 2))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4))+' and DD values of <='+str(round(out_dd, 4)))
            print(str(sum(1 for x in df_dat.fstdd if x == 2))+'\toutlier windows for Fst values of >='+str(round(out_fst, 4))+' and DD values of <='+str(round(out_dd, 4)))
            print(str(sum(1 for x in df_dat.dxyfstdd if x == 3))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4))+', Fst values of >='+str(round(out_fst, 4))+' and DD values of <='+str(round(out_dd, 4))+'\n')
            print('Select outlier windows\n')
            
            df_dat.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_allsites.csv', index=False)
            # select all windows that are outliers for at least one metric
            df_outlier = df_dat[((df_dat.dxyout != 0) | (df_dat.fstout != 0) | (df_dat.ddout != 0))]
            df_outlier.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_ALL_outliers.csv', index=False)
            # write bedfile for gene retrieval and orthologous gene onthology match
            header = ["scaffold", "bedstart", "end"]
            df_outlier.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_ALL_outliers.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

            # select Dxy only outliers only
            df_dxy = df_dat[((df_dat.dxyout == 1) & (df_dat.fstout == 0) & (df_dat.ddout == 0))]
            if df_dxy.empty:
                print('No single Dxy outlier windows')
            else:
                df_dxy.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_Dxy_outliers.csv', index=False)
                # write bedfile for gene retrieval and orthologous gene onthology match
                header = ["scaffold", "bedstart", "end"]
                df_dxy.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_Dxy_outliers.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

            # select Fst only outliers only
            df_fst = df_dat[((df_dat.dxyout == 0) & (df_dat.fstout == 1) & (df_dat.ddout == 0))]
            if df_fst.empty:
                print('No single Fst outlier windows')
            else:
                df_fst.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_Fst_outliers.csv', index=False)
                # write bedfile for gene retrieval and orthologous gene onthology match
                header = ["scaffold", "bedstart", "end"]
                df_fst.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_Fst_outliers.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)
                
            # select DD only outliers only
            df_dd = df_dat[((df_dat.dxyout == 0) & (df_dat.fstout == 0) & (df_dat.ddout == 1))]
            if df_dd.empty:
                print('No single DD outlier windows')
            else:
                df_dd.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_DD_outliers.csv', index=False)
                # write bedfile for gene retrieval and orthologous gene onthology match
                header = ["scaffold", "bedstart", "end"]
                df_dd.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_DD_outliers.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

            # select Dxy and Fst double outliers only
            df_dxyfst = df_dat[((df_dat.dxyfst == 2) & (df_dat.dxyfstdd != 3))]
            if df_dxyfst.empty:
                print('No double Dxy Fst outlier windows')
            else:
                df_dxyfst.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_DxyFst_outliers.csv', index=False)
                # write bedfile for gene retrieval and orthologous gene onthology match
                header = ["scaffold", "bedstart", "end"]
                df_dxyfst.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_DxyFst_outliers.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

            # select Dxy and DD double outliers only
            df_dxydd = df_dat[((df_dat.dxydd == 2) & (df_dat.dxyfstdd != 3))]
            if df_dxydd.empty:
                print('No double Dxy DD outlier windows')
            else:
                df_dxydd.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_DxyDD_outliers.csv', index=False)
                # write bedfile for gene retrieval and orthologous gene onthology match
                header = ["scaffold", "bedstart", "end"]
                df_dxydd.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_DxyDD_outliers.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

            # select Fst and DD double outliers only
            df_fstdd = df_dat[((df_dat.fstdd == 2) & (df_dat.dxyfstdd != 3))]
            if df_fstdd.empty:
                print('No Fst DD double outlier windows')
            else:
                df_fstdd.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_FstDD_outliers.csv', index=False)
                # write bedfile for gene retrieval and orthologous gene onthology match
                header = ["scaffold", "bedstart", "end"]
                df_fstdd.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_FstDD_outliers.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

            # select Dxy, Fst, DD triple outliers only
            df_dxyfstdd = df_dat[(df_dat.dxyfstdd == 3)]
            if df_dxyfstdd.empty:
                print('No Dxy Fst DD triple outlier windows')
            else:
                df_dxyfstdd.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_DxyFstDD_outliers.csv', index=False)
                # write bedfile for gene retrieval and orthologous gene onthology match
                header = ["scaffold", "bedstart", "end"]
                df_dxyfstdd.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_DxyFstDD_outliers.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)


else: #this is the else argument connecting to the top of the script               
    ###### STEP 5 ######
    # define outliers using confidence intervall assuming a normal distribution
    print('\n\t Popgen metrics have already been computed for '+str(args.snps)+'SNPs windows:')
    print('-> jump to')
    print('\tSTEP 5: Define and select outliers')
    out = []
    for dirName, subdirList, fileList in os.walk(outputdir):
        for file in fileList:
            if file.startswith(contrast) and file.endswith(str(args.snps)+'SNPs_3metrics.txt') == True:
                out.append(dirName+'/'+file)
    out_sort = natsorted(out)

    for out in out_sort:
        infile = open(out, 'r')
        header = infile.readline()
        dxy = []
        fst = []
        dd = []
        length = []
        
        for line in infile:
            data = line.split()
            dxy.append(float(data[13]))
            fst.append(float(data[12]))
            dd.append(float(data[15]))
            length.append(int(data[4]))
        mean_Dxy = statistics.mean(dxy)
        median_Dxy = statistics.median(dxy)
        sd_Dxy = statistics.stdev(dxy)
        CV_Dxy = sd_Dxy/mean_Dxy
        out_dxy = mean_Dxy+args.ci*sd_Dxy
        mean_Fst = statistics.mean(fst)
        median_Fst = statistics.median(fst)
        sd_Fst = statistics.stdev(fst)
        CV_Fst = sd_Fst/mean_Fst
        out_fst = mean_Fst+args.ci*sd_Fst
        mean_DD = statistics.mean(dd)
        median_DD = statistics.median(dd)
        sd_DD = statistics.stdev(dd)
        CV_DD = sd_DD/mean_DD
        out_dd = mean_DD-args.ci*sd_DD
        mean_win = statistics.mean(length)
        median_win = statistics.median(length)
        sd_win = statistics.stdev(length)
        CV_win = sd_win/mean_win
        out_win = mean_win+args.ci*sd_win
        
        outlier_values = open(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_outlier_descriptive_stats.txt', 'w')
        outlier_values.write('stat\tDxy\tFst\tDD\twin_length\n'+
                             'mean\t'+str(round(mean_Dxy,4))+'\t'+str(round(mean_Fst,4))+'\t'+str(round(mean_DD,4))+'\t'+str(round(mean_win,0))+'\n'+
                             'median\t'+str(round(median_Dxy,4))+'\t'+str(round(median_Fst,4))+'\t'+str(round(median_DD,4))+'\t'+str(round(median_win,0))+'\n'+
                             'sd\t'+str(round(sd_Dxy,4))+'\t'+str(round(sd_Fst,4))+'\t'+str(round(sd_DD,4))+'\t'+str(round(sd_win,0))+'\n'+
                             'CV\t'+str(round(CV_Dxy,4))+'\t'+str(round(CV_Fst,4))+'\t'+str(round(CV_DD,4))+'\t'+str(round(CV_win,0))+'\n'+
                             'CI\t'+str(round(out_dxy,4))+'\t'+str(round(out_fst,4))+'\t'+str(round(out_dd,4))+'\t'+str(round(out_win,0))+'\n')
        
        infile.close()

        # use pandas to read infiles as dataframe and select outlier combinations
        df_dat = pd.read_table(out)
        df_dat['dxyout'] = np.where(df_dat.Dxy >= out_dxy, 1, 0)
        df_dat['fstout'] = np.where(df_dat.Fst >= out_fst, 1, 0)
        df_dat['ddout'] = np.where(df_dat.DD <= out_dd, 1, 0)
        df_dat['dxyfst'] = df_dat.dxyout+df_dat.fstout
        df_dat['dxydd'] = df_dat.dxyout+df_dat.ddout
        df_dat['fstdd'] = df_dat.fstout+df_dat.ddout
        df_dat['dxyfstdd'] = df_dat.dxyout+df_dat.fstout+df_dat.ddout
        df_dat['bedstart'] = df_dat.start-1
        print('\nFor confidence intervals of '+str(args.ci)+' standard deviations we find:\n')
        print(str(sum(1 for x in df_dat.dxyout if x == 1))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4)))
        print(str(sum(1 for x in df_dat.fstout if x == 1))+'\toutlier windows for Fst values of >='+str(round(out_fst, 4)))
        print(str(sum(1 for x in df_dat.ddout if x == 1))+'\toutlier windows for DD values of <='+str(round(out_dd, 4)))
        print(str(sum(1 for x in df_dat.dxyfst if x == 2))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4))+' and Fst values of >='+str(round(out_fst, 4)))
        print(str(sum(1 for x in df_dat.dxydd if x == 2))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4))+' and DD values of <='+str(round(out_dd, 4)))
        print(str(sum(1 for x in df_dat.fstdd if x == 2))+'\toutlier windows for Fst values of >='+str(round(out_fst, 4))+' and DD values of <='+str(round(out_dd, 4)))
        print(str(sum(1 for x in df_dat.dxyfstdd if x == 3))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4))+', Fst values of >='+str(round(out_fst, 4))+' and DD values of <='+str(round(out_dd, 4))+'\n')
        print('Select oulier windows\n')
        
        df_dat.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_allsites.csv', index=False)
        # select all windows that are outliers for at least one metric
        df_outlier = df_dat[((df_dat.dxyout != 0) | (df_dat.fstout != 0) | (df_dat.ddout != 0))]
        df_outlier.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_ALL_outliers.csv', index=False)
        # write bedfile for gene retrieval and orthologous gene onthology match
        header = ["scaffold", "bedstart", "end"]
        df_outlier.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_ALL_outliers.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

        # select Dxy only outliers only
        df_dxy = df_dat[((df_dat.dxyout == 1) & (df_dat.fstout == 0) & (df_dat.ddout == 0))]
        if df_dxy.empty:
            print('No single Dxy outlier windows')
        else:
            df_dxy.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_Dxy_outliers.csv', index=False)
            # write bedfile for gene retrieval and orthologous gene onthology match
            header = ["scaffold", "bedstart", "end"]
            df_dxy.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_Dxy_outliers.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

        # select Fst only outliers only
        df_fst = df_dat[((df_dat.dxyout == 0) & (df_dat.fstout == 1) & (df_dat.ddout == 0))]
        if df_fst.empty:
            print('No single Fst outlier windows')
        else:
            df_fst.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_Fst_outliers.csv', index=False)
            # write bedfile for gene retrieval and orthologous gene onthology match
            header = ["scaffold", "bedstart", "end"]
            df_fst.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_Fst_outliers.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)
            
        # select DD only outliers only
        df_dd = df_dat[((df_dat.dxyout == 0) & (df_dat.fstout == 0) & (df_dat.ddout == 1))]
        if df_dd.empty:
            print('No single DD outlier windows')
        else:
            df_dd.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_DD_outliers.csv', index=False)
            # write bedfile for gene retrieval and orthologous gene onthology match
            header = ["scaffold", "bedstart", "end"]
            df_dd.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_DD_outliers.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

        # select Dxy and Fst double outliers only
        df_dxyfst = df_dat[((df_dat.dxyfst == 2) & (df_dat.dxyfstdd != 3))]
        if df_dxyfst.empty:
            print('No double Dxy Fst outlier windows')
        else:
            df_dxyfst.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_DxyFst_outliers.csv', index=False)
            # write bedfile for gene retrieval and orthologous gene onthology match
            header = ["scaffold", "bedstart", "end"]
            df_dxyfst.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_DxyFst_outliers.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

        # select Dxy and DD double outliers only
        df_dxydd = df_dat[((df_dat.dxydd == 2) & (df_dat.dxyfstdd != 3))]
        if df_dxydd.empty:
            print('No double Dxy DD outlier windows')
        else:
            df_dxydd.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_DxyDD_outliers.csv', index=False)
            # write bedfile for gene retrieval and orthologous gene onthology match
            header = ["scaffold", "bedstart", "end"]
            df_dxydd.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_DxyDD_outliers.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

        # select Fst and DD double outliers only
        df_fstdd = df_dat[((df_dat.fstdd == 2) & (df_dat.dxyfstdd != 3))]
        if df_fstdd.empty:
            print('No Fst DD double outlier windows')
        else:
            df_fstdd.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_FstDD_outliers.csv', index=False)
            # write bedfile for gene retrieval and orthologous gene onthology match
            header = ["scaffold", "bedstart", "end"]
            df_fstdd.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_FstDD_outliers.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

        # select Dxy, Fst, DD triple outliers only
        df_dxyfstdd = df_dat[(df_dat.dxyfstdd == 3)]
        if df_dxyfstdd.empty:
            print('No Dxy Fst DD triple outlier windows')
        else:
            df_dxyfstdd.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_DxyFstDD_outliers.csv', index=False)
            # write bedfile for gene retrieval and orthologous gene onthology match
            header = ["scaffold", "bedstart", "end"]
            df_dxyfstdd.to_csv(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(args.ci)+'sd_DxyFstDD_outliers.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

    print('\n Succesfully found selective sweep candidates for '+args.coh1+'_'+args.coh2+' contrast and '+str(args.snps)+' SNP windows\n\n')
    print('\n\tDONE\n')
