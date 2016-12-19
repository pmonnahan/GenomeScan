## genome scan pipeline class file
## has to be loaded within the Gx_x.py genome scan scripts
## by Christian Sailer with help from Katie Barr, 5 December 2016
import os, sys, subprocess, statistics, argparse
from natsort import natsorted
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import math


## create Class 
class G1():
    def __init__(self, AFs_table=None, three_metrics=None, args=None, natsorted=natsorted):
        self.AFs_table, self.three_metrics = None, None
        self.args = args
        self.natsorted = natsorted

    def step1(self):
        ###### STEP 1 ######
        ## Prepare input data for analysis
        # test where to place the outputdirectory
        args = self.args
        cwd = os.getcwd()
        contrast = args.i
        
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
        contrast = args.coh1+args.coh2
        
        print('\nSearching input directory for *_raw.table')
        tincoh1 = []
        tincoh2 = []
        for dirName, subdirList, fileList in os.walk(args.i, topdown = False):
            for fname in fileList:
                if fname.endswith('_raw.table') and args.coh1 in fname:
                    tincoh1.append(dirName+'/'+fname)
                elif fname.endswith('_raw.table') and args.coh2 in fname:
                    tincoh2.append(dirName+'/'+fname)
        incoh1 = natsorted(tincoh1)
        incoh2 = natsorted(tincoh2)
                                    
        # paste cohort 1 & cohort 2 into a table next to each other & remove CHROM POS from second cohort in joint table
        # yields: CHROM POS AC AN DP AC AN DP
        print('\n\tSTEP 1: Paste contrast AC tables\n')
        for i in range(len(incoh1)):
            print('Processing '+incoh1[i])
            pastecmd = open('paste_'+args.coh1+'_'+args.coh2+'.unix', 'w')
            pastecmd.write('paste ')
            pastecmd.write(incoh1[i]+' '+incoh2[i]+' | cut -f -5,8,9,10 ')
            pastecmd.write('> '+outputdir+'/'+contrast+'_scaf_'+str(i+1)+'_temp.table')
            pastecmd.close()

            # run in unix
            cmd = (open('paste_'+args.coh1+'_'+args.coh2+'.unix', 'r'))
            p = subprocess.Popen(cmd, shell=True)
            sts = os.waitpid(p.pid, 0)[1]

        os.remove('paste_'+args.coh1+'_'+args.coh2+'.unix') # removes unix script file from folder

        # variables to be inherited
        self.outputdir = outputdir
        self.contrast = contrast

    def step2(self):
        ###### STEP 2 ######
        ## This section is Jeff's FixedDerivedAlleleCheck.py
        # search directory for output of above step to make new input list
        args = self.args
        outputdir = self.outputdir
        contrast = self.contrast
        tintable = []
        for dirName, subdirList, fileList in os.walk(outputdir):
            for table in fileList:
                if table.startswith(contrast) and table.endswith('_temp.table'):
                    tintable.append(dirName+'/'+table)
        intable = natsorted(tintable)
        
        print('\n\tSTEP 2: Remove sites that are fixed in both cohorts and are missing in one but not the other cohort\n')
        print('Found '+str(len(intable))+' input tables')

        # obtain ANmax for each cohort
        print('\nGathering allele count info:\n')
        for i in range(len(intable)):
            with open(intable[i], 'r') as infile:
                infile.readline()
                tANmax_1 = []
                tANmax_2 = []
                for line in infile:
                    scaffold, position, AC_1, AN_1, DP_1, AC_2, AN_2, DP_2 = line.split()
                    tANmax_1.append(int(AN_1))
                    tANmax_2.append(int(AN_2))
                ANmax_1 = max(tANmax_1)
                ANmax_2 = max(tANmax_2)
                print('For '+intable[i]+': ')
                print('ANmax = '+str(ANmax_1)+' in cohort '+args.coh1)
                print('ANmax = '+str(ANmax_2)+' in cohort '+args.coh2)

                # filter
                infile.seek(0)
                with open(outputdir+'/'+contrast+'_scaf_'+str(i+1)+'.tab','w') as outfile:
                    header = infile.readline()
                    outfile.write(header)

                    count = 0
                    target = 100000
                    Zero_count = 0
                    All_count = 0

                    for line in infile:
                        scaffold, position, AC_1, AN_1, DP_1, AC_2, AN_2, DP_2 = line.split()
                        if int(AN_1)>=int(round((1-args.per)*ANmax_1,0)) and int(AN_2)>=int(round((1-args.per)*ANmax_2,0)):
                            total_allele_count = int(AC_1)+int(AC_2)
                            total_alleles_sampled = int(AN_1)+int(AN_2)
                            if total_allele_count == total_alleles_sampled:
                                All_count +=1
                            elif total_allele_count == 0:
                                Zero_count +=1
                            else:
                                outfile.write(line)
                            count += 1

                    print('Total variants in file: '+str(count))
                    print('Variants with total allele count = 0: '+str(Zero_count))
                    print('Variants fixed in both cohorts: '+str(All_count))
                    print('Retained '+str(count-Zero_count-All_count)+' variants\n')
                    os.remove(outputdir+'/'+contrast+'_scaf_'+str(i+1)+'_temp.table') # removes output files of part 1, as they are no longer needed


    def step3(self):
        ###### STEP 3 ######
        ## jump highdepthFilter and remove DP column, Jeff's script adjusted
        args = self.args
        outputdir = self.outputdir
        contrast = self.contrast

        print('\n\tSTEP 3: Remove DP column')
        tintable = []
        for dirName, subdirList, fileList in os.walk(outputdir):
            for table in fileList:
                if table.startswith(contrast) and table.endswith('.tab'):
                    tintable.append(dirName+'/'+table)
        intable = natsorted(tintable)

        for i in range(len(intable)):
            with open(intable[i],'r') as infile:
                infile.readline()
                with open(outputdir+'/'+contrast+'_scaf_'+str(i+1)+'_AFs'+args.suf+'.table','w') as outfile:
                    outfile.write('CHROM\tPOS\tAC1\tAN1\tAC2\tAN2\n')

                    for line in infile:
                        scaffold, position, AC_1, AN_1, DP_1, AC_2, AN_2, DP_2 = line.split()
                        outfile.write(scaffold+'\t'+position+'\t'+AC_1+'\t'+AN_1+'\t'+AC_2+'\t'+AN_2+'\n')
                    os.remove(outputdir+'/'+contrast+'_scaf_'+str(i+1)+'.tab') # removes output files of part 2, as they are no longer needed
            

    def step4(self):
        ###### STEP 4 ######
        # calculate Dxy, Fst, pi, allele frequency and allele frequncy difference. Input required as:
        # scaffold, position, cohort1_allele_count, cohort1_all_alleles_count, cohort2_allele_count, cohort2_all_alleles_count
        args = self.args
        if args.o is 'na':
            outputdir = str(args.coh1+args.coh2)
        else:
            outputdir = str(args.o+args.coh1+args.coh2)

        contrast = args.coh1 + args.coh2
        self.outputdir = outputdir
        self.contrast = contrast
        print('\n\tSTEP 4: Calculate Dxy, Fst, DD')
        header = 'scaffold\tstart\tend\tmidpoint\tlength\tnum_snps\t'+args.coh1+'_freq\t'+args.coh2+'_freq\tabsdiff\t'+args.coh1+'_'+args.coh2+'\t'+args.coh2+'_'+args.coh1+'\tvaru\t'+args.coh1+'_pi\t'+args.coh2+'_pi\tFst\tDxy\n'
        with open(outputdir+'/'+contrast+'_metrics_WG_'+str(args.snps)+'SNPs'+args.suf+'.txt','w') as outfile:
            outfile.write(header)

            AFs = []
            for dirName, subdirList, fileList in os.walk(outputdir):
                for file in fileList:
                    if file.startswith(contrast) and file.endswith('_AFs'+args.suf+'.table') == True:
                        AFs.append(dirName+'/'+file)
            AFs_sort = natsorted(AFs)

            print('\nFound '+str(len(AFs_sort))+' files starting with '+contrast+' and ending with _AFs'+args.suf+'table\n')

            file_count = 0
            win_count = 0
            winexclcount = 0

            for AF in AFs_sort:
                count = -1 # to remove header from count
                with open(AF,'r') as infile:
                    for line in infile:
                        print(line)
                        count += 1
                    num_win = int(count/args.snps)

                    print('Processing '+str(num_win)+' windows in '+AF)

                if num_win > 0:
                    file_count += 1
                    win_count += num_win

                    with open(AF,'r') as infile:
                        infile.readline() # to read the header and start reading the numbers for calcluation from the second line of the file
                        if args.win_type=='snp':
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
                                varu = []
                                for i in range(args.snps):
                                    data = infile.readline()
                                    scaffold, position, AC_1, AN_1, AC_2, AN_2 = data.split()
                                    scaf.append(scaffold)
                                    pos.append(int(position))
                                    allele_c1_count = int(AC_1)
                                    allele_c2_count = int(AC_2)
                                    allele_c1_freq = allele_c1_count/int(AN_1)
                                    allele_c2_freq = allele_c2_count/int(AN_2)
                                    allele_count = allele_c1_count + allele_c2_count
                                    allele_freq = allele_count/(int(AN_1) + int(AN_2))
                                    c1_freq.append(allele_c1_freq)
                                    c2_freq.append(allele_c2_freq)
                                    c1_minus_c2.append(allele_c1_freq - allele_c2_freq)
                                    c2_minus_c1.append(allele_c2_freq - allele_c1_freq)
                                    absdiff.append(abs(allele_c1_freq-allele_c2_freq))
                                    c1_pi.append((2*allele_c1_freq*(1-allele_c1_freq)))    
                                    c2_pi.append((2*allele_c2_freq*(1-allele_c2_freq)))
                                    ht = 2*allele_freq*(1-allele_freq)
                                    h1 = 2*allele_c1_freq*(1-allele_c1_freq)
                                    h2 = 2*allele_c2_freq*(1-allele_c2_freq)
                                    h12 = ((h1*int(AN_1))+(h2*int(AN_2)))/(int(AN_1) + int(AN_2))
                                    allele_fst = abs(ht-h12)/ht
                                    fst.append(allele_fst)
                                    dxy.append((allele_c1_freq*(1-allele_c2_freq))+(allele_c2_freq*(1-allele_c1_freq)))
                                    varu1 = (1/(int(AN_1)*int(AN_1)))*allele_c1_freq*(1-allele_c1_freq)
                                    varu2 = (1/(int(AN_2)*int(AN_2)))*allele_c2_freq*(1-allele_c2_freq)
                                    varu.append(varu1 + varu2)

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
                                mean_varu = statistics.mean(varu)

                                if wlength <= 26560:
                        #             string = "scaf[0]\t %d \t %d \t %d" % (wstart, wend, wmid)
                        #             outfile.write(string)
                                    outfile.write(scaf[0]+'\t'+
                                                  str(wstart)+'\t'+
                                                  str(wend)+'\t'+
                                                  str(wmid)+'\t'+
                                                  str(wlength)+'\t'+
                                                  str(args.snps)+'\t'+
                                                  str(mean_c1_freq)+'\t'+
                                                  str(mean_c2_freq)+'\t'+
                                                  str(mean_absdiff)+'\t'+
                                                  str(mean_c1_minus_c2)+'\t'+
                                                  str(mean_c2_minus_c1)+'\t'+
                                                  str(mean_varu)+'\t'+
                                                  str(mean_c1_pi)+'\t'+
                                                  str(mean_c2_pi)+'\t'+
                                                  str(mean_fst)+'\t'+
                                                  str(win_dxy)+'\n')
                                else:
                                    winexclcount +=1
                        #GEOGRAPHIC WINDOW ANALYSIS
                        else:
                            AN1=math.ceil(ANmax_1*args.m)
                            AN2=math.ceil(ANmax_2*args.m)
                            AFS1=[0 for cat in range(0,AN1)]
                            AFS2=[0 for cat in range(0,AN2)]
                            AFS12=np.zeros((AN1,AN2),dtype=int)
                            # for AF in AFs_sort:
                            #     count = -1 # to remove header from count
                            #     infile = open(AF,'r')
                            #     for line in infile:
                            #         count += 1
                            #         line=line.split("\t")
                            #         pos=int(line[1])
                            #     maxpos = pos + args.ws
                            #     infile.close()
                            #     num_win = int(maxpos/args.ws)
                            start=0
                            end=args.ws
                            idx=0
                            snp_count = 0
                            scafs = []
                            Start = []
                            End = []
                            Snp_counts = [0 for jjj in range(0,num_win)]
                            c1_freq = []
                            c2_freq = []
                            absdiff=[]
                            c1_minus_c2 = []
                            c2_minus_c1 = []
                            c1_pi = []
                            c2_pi = []
                            fst = []
                            dxy = []

                            win_afs1=[0 for cat in range(0,AN1)]
                            win_afs2=[0 for cat in range(0,AN2)]
                            win_afs12=np.zeros((AN1,AN2),dtype=int)

                            for line in infile:
                                data=line.split("\t")
                                scaf = data[0]
                                pos = int(data[1])
                                allele_c1_count = int(data[2])
                                allele_c2_count = int(data[4])
                                allele_c1_freq = allele_c1_count/int(data[3])
                                allele_c2_freq = allele_c2_count/int(data[5])
                                allele_count = allele_c1_count + allele_c2_count
                                allele_freq = allele_count/(int(data[3]) + int(data[5]))

                                if pos > start and pos <=end:
                                    snp_count+=1
                                    win_afs1[allele_c1_count] += 1
                                    win_afs2[allele_c2_count] += 1
                                    win_afs12[allele_c1_count,allele_c2_count]+=1
                                    AFS1[allele_c1_count] += 1
                                    AFS2[allele_c2_count] += 1
                                    AFS12[allele_c1_count,allele_c2_count]+=1
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
                                    print(start, end, pos, allele_c1_count,allele_c2_count)
                                    print(win_afs1)
                                    print(win_afs2)
                                    print(AFS1)
                                    print(AFS2)


                                elif pos > end:
                                    Snp_counts[idx]=snp_count
                                    if snp_count >= args.ms:
                                            outfile.write(scaf[0]+'\t'+
                                            str(start)+'\t'+
                                            str(end)+'\t'+
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

                                    #SETUP NEW WINDOW
                                    win_afs1=[0 for cat in range(0,AN1)]
                                    win_afs2=[0 for cat in range(0,AN2)]
                                    win_afs12=np.zeros((AN1,AN2),dtype=int)
                                    c1_freq = []
                                    c2_freq = []
                                    absdiff=[]
                                    c1_minus_c2 = []
                                    c2_minus_c1 = []
                                    c1_pi = []
                                    c2_pi = []
                                    fst = []
                                    dxy = []

                                    while pos > end:
                                        idx+=1
                                        end+=args.ws/2

                                    start=end-args.ws

                                    assert (pos < end and pos > start)

                                    snp_count=1
                                    win_afs1[allele_c1_count] += 1
                                    win_afs2[allele_c2_count] += 1
                                    win_afs12[allele_c1_count,allele_c2_count]+=1
                                    AFS1[allele_c1_count] += 1
                                    AFS2[allele_c2_count] += 1
                                    AFS12[allele_c1_count,allele_c2_count]+=1
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

                                else:
                                    print("SNP NOT CAUGHT BY WINDOW")

        print('\nAnalyzed '+str(win_count)+' windows in '+str(file_count)+' files')
        print('Excluded '+str(winexclcount)+' windows longer than 26560bp\n')

        ###### STEP 4b ######
        # Add DD residuals using R script
        inlist = []
        for dirName, subdirList, fileList in os.walk(outputdir):
            for file in fileList:
                if file.startswith(contrast) and file.endswith(str(args.snps)+'SNPs'+args.suf+'.txt'):
                    inlist.append(dirName+'/'+file)
        inlist_sort = natsorted(inlist)

        print('Add DDresiduals')

        for file in inlist_sort:
            count = 0

            snp_file = pd.read_table(file, header=0)
            slope, intercept, r_value, p_value, std_err = stats.linregress(snp_file['absdiff'], snp_file[args.coh1+'_pi'])
            snp_file['prediction'] = intercept + (slope*snp_file['absdiff'])
            print('Uncorrected slope '+str(slope)+' uncorrected intercept '+str(intercept))

        # calculate slope correction using var(u)
            varu = snp_file.varu
            mean_absdiff = statistics.mean(snp_file.absdiff)
            mean_pi1 = statistics.mean(snp_file[args.coh1+'_pi'])
            var_absdiff = statistics.variance(snp_file.absdiff)
            cor_slope = statistics.mean((var_absdiff/(var_absdiff - varu)))*slope
            cor_intercept = mean_pi1-(mean_absdiff*cor_slope)

            print('\nVar(u) of '+str(file)+' is '+str(statistics.mean(varu)))
            print('Corrected slope '+str(cor_slope)+' corrected intercept '+str(cor_intercept))
            
            snp_file['cor_prediction'] = cor_intercept + (cor_slope*snp_file.absdiff)
            snp_file['DD'] = snp_file[args.coh1+'_pi'] - snp_file.cor_prediction

            snp_file.to_csv(outputdir+'/'+contrast+'_WG_'+str(args.snps)+'SNPs_3metrics'+args.suf+'.txt',sep="\t", index=False)

            count += 1
        print('Processed '+str(count)+' files for DDresiduals')
        os.remove(outputdir+'/'+contrast+'_metrics_WG_'+str(args.snps)+'SNPs'+args.suf+'.txt')


    def step5(self):
        ###### STEP 5 ######
        # define outliers using top percentile cut-off
        print('\n\tSTEP 5: Select outlier windows')
        args = self.args
        cwd = os.getcwd()
        if args.o is 'na':
            outputdir = str(args.coh1+args.coh2)
        else:
            outputdir = str(args.o+args.coh1+args.coh2)

        contrast = args.coh1 + args.coh2
        self.outputdir = outputdir
        self.contrast = contrast
        self.cwd = cwd

        out = []
        for dirName, subdirList, fileList in os.walk(outputdir):
            for file in fileList:
                if file.startswith(contrast) and file.endswith(str(args.snps)+'SNPs_3metrics'+args.suf+'.txt') == True:
                    out.append(dirName+'/'+file)
        out_sort = natsorted(out)

        for out in out_sort:
            with open(out, 'r') as infile:
                header = infile.readline()
                afd_c1_c2 = []
                afd_c2_c1 = []
                dxy = []
                fst = []
                dd = []
                length = []
                abs_afd = []
                
                for line in infile:
                    data = line.split()
                    afd_c1_c2.append(float(data[8]))
                    dxy.append(float(data[14]))
                    fst.append(float(data[13]))
                    dd.append(float(data[17]))
                    length.append(int(data[4]))

#                    dxy.append(float(Dxy))
#                    fst.append(float(Fst))
#                    dd.append(float(DD))
#                    varu.append(float(varu))
#                    length.append(int(length))
                    abs_afd.append(float(data[9]))
                mean_afd = statistics.mean(afd_c1_c2)
                median_afd = statistics.mean(afd_c1_c2)
                sd_afd = statistics.stdev(afd_c1_c2)
                CV_afd = sd_afd/mean_afd
                out_afd_up = np.percentile(afd_c1_c2, (100-args.cut))
                out_afd_low = np.percentile(afd_c1_c2, args.cut)
                mean_Dxy = statistics.mean(dxy)
                median_Dxy = statistics.median(dxy)
                sd_Dxy = statistics.stdev(dxy)
                CV_Dxy = sd_Dxy/mean_Dxy
                out_dxy = np.percentile(dxy, (100-args.cut))
                mean_Fst = statistics.mean(fst)
                median_Fst = statistics.median(fst)
                sd_Fst = statistics.stdev(fst)
                CV_Fst = sd_Fst/mean_Fst
                out_fst = np.percentile(fst, (100-args.cut))
                mean_DD = statistics.mean(dd)
                median_DD = statistics.median(dd)
                sd_DD = statistics.stdev(dd)
                CV_DD = sd_DD/mean_DD
                out_dd = np.percentile(dd, args.cut)
                mean_win = statistics.mean(length)
                median_win = statistics.median(length)
                sd_win = statistics.stdev(length)
                CV_win = sd_win/mean_win
                out_win = np.percentile(length, (100-args.cut))
                
                outlier_values = open(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(int(10000*args.cut))+'ppm_descriptive_stats'+args.suf+'.txt', 'w')
                outlier_values.write('stat\tAFD_1_2\tDxy\tFst\tDD\twin_length\n'+
                                     'min\t'+str(round(min(afd_c1_c2),4))+'\t'+str(round(min(dxy),4))+'\t'+str(round(min(fst),4))+'\t'+str(round(min(dd),4))+'\t'+str(round(min(length),0))+'\n'+
                                     'mean\t'+str(round(mean_afd,4))+'\t'+str(round(mean_Dxy,4))+'\t'+str(round(mean_Fst,4))+'\t'+str(round(mean_DD,4))+'\t'+str(round(mean_win,0))+'\n'+
                                     'median\t'+str(round(median_afd,4))+'\t'+str(round(median_Dxy,4))+'\t'+str(round(median_Fst,4))+'\t'+str(round(median_DD,4))+'\t'+str(round(median_win,0))+'\n'+
                                     'sd\t'+str(round(sd_afd,4))+'\t'+str(round(sd_Dxy,4))+'\t'+str(round(sd_Fst,4))+'\t'+str(round(sd_DD,4))+'\t'+str(round(sd_win,0))+'\n'+
                                     'CV\t'+str(round(CV_afd,4))+'\t'+str(round(CV_Dxy,4))+'\t'+str(round(CV_Fst,4))+'\t'+str(round(CV_DD,4))+'\t'+str(round(CV_win,0))+'\n'+
                                     'max\t'+str(round(max(afd_c1_c2),4))+'\t'+str(round(max(dxy),4))+'\t'+str(round(max(fst),4))+'\t'+str(round(max(dd),4))+'\t'+str(round(max(length),0))+'\n'+
                                     'Cutoff\t'+str(round(out_afd_up,4))+'\t'+str(round(out_dxy,4))+'\t'+str(round(out_fst,4))+'\t'+str(round(out_dd,4))+'\t'+str(round(out_win,0))+'\n'+
                                     'Cutoff2\t'+str(round(out_afd_low,4))+'\t'+str(round(out_dxy,4))+'\t'+str(round(out_fst,4))+'\t'+str(round(out_dd,4))+'\t'+str(round(out_win,0))+'\n')
                

                # use pandas to read infiles as dataframe and select outlier combinations
                c1_c2 = str(args.coh1+'_'+args.coh2)
                df_dat = pd.read_table(out)
                df_dat['afdout'] = np.where(df_dat[c1_c2] >= out_afd_up, 1, 0)
                df_dat['afd21out'] = np.where(df_dat[c1_c2] <= out_afd_low, 1, 0)
                df_dat['dxyout'] = np.where(df_dat.Dxy >= out_dxy, 1, 0)
                df_dat['fstout'] = np.where(df_dat.Fst >= out_fst, 1, 0)
                df_dat['ddout'] = np.where(df_dat.DD <= out_dd, 1, 0)
                # doubles
                df_dat['afddxy'] = df_dat.afdout+df_dat.dxyout
                df_dat['afdfst'] = df_dat.afdout+df_dat.fstout
                df_dat['afddd'] = df_dat.afdout+df_dat.ddout
                df_dat['afd21dxy'] = df_dat.afd21out+df_dat.dxyout
                df_dat['afd21fst'] = df_dat.afd21out+df_dat.fstout
                df_dat['afd21dd'] = df_dat.afd21out+df_dat.ddout
                df_dat['dxyfst'] = df_dat.dxyout+df_dat.fstout
                df_dat['dxydd'] = df_dat.dxyout+df_dat.ddout
                df_dat['fstdd'] = df_dat.fstout+df_dat.ddout
                # triples
                df_dat['afddxyfst'] = df_dat.afdout+df_dat.dxyout+df_dat.fstout
                df_dat['afddxydd'] = df_dat.afdout+df_dat.dxyout+df_dat.ddout
                df_dat['afdfstdd'] = df_dat.afdout+df_dat.fstout+df_dat.ddout
                df_dat['afd21dxyfst'] = df_dat.afd21out+df_dat.dxyout+df_dat.fstout
                df_dat['afd21dxydd'] = df_dat.afd21out+df_dat.dxyout+df_dat.ddout
                df_dat['afd21fstdd'] = df_dat.afd21out+df_dat.fstout+df_dat.ddout
                df_dat['dxyfstdd'] = df_dat.dxyout+df_dat.fstout+df_dat.ddout
                # quadruple
                df_dat['afddxyfstdd'] = df_dat.afdout+df_dat.dxyout+df_dat.fstout+df_dat.ddout
                df_dat['afd21dxyfstdd'] = df_dat.afd21out+df_dat.dxyout+df_dat.fstout+df_dat.ddout

                df_dat['bedstart'] = df_dat.start-1
                print('\nFor top '+str(args.cut)+' percent we find:\n')
                print(str(sum(1 for x in df_dat.afdout if x == 1))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4)))
                print(str(sum(1 for x in df_dat.afd21out if x == 1))+'\toutlier windows for AFD('+args.coh2+'-'+args.coh1+') values of <='+str(round(out_afd_low, 4)))
                print(str(sum(1 for x in df_dat.dxyout if x == 1))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4)))
                print(str(sum(1 for x in df_dat.fstout if x == 1))+'\toutlier windows for Fst values of >='+str(round(out_fst, 4)))
                print(str(sum(1 for x in df_dat.ddout if x == 1))+'\toutlier windows for DD values of <='+str(round(out_dd, 4)))
                print('\tof which are double outliers for:')
                print(str(sum(1 for x in df_dat.afddxy if x == 2))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+' and Dxy values of >='+str(round(out_dxy, 4)))
                print(str(sum(1 for x in df_dat.afdfst if x == 2))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+' and Fst values of >='+str(round(out_fst, 4)))
                print(str(sum(1 for x in df_dat.afddd if x == 2))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+' and DD values of >='+str(round(out_dd, 4)))
                print(str(sum(1 for x in df_dat.afd21dxy if x == 2))+'\toutlier windows for AFD('+args.coh2+'-'+args.coh1+') values of <='+str(round(out_afd_low, 4))+' and Dxy values of >='+str(round(out_dxy, 4)))
                print(str(sum(1 for x in df_dat.afd21fst if x == 2))+'\toutlier windows for AFD('+args.coh2+'-'+args.coh1+') values of <='+str(round(out_afd_low, 4))+' and Fst values of >='+str(round(out_fst, 4)))
#                print(str(sum(1 for x in df_dat.afd21dd if x == 2))+'\toutlier windows for AFD_2_1 values of <='+str(round(out_afd_low, 4))+' and DD values of >='+str(round(out_dd, 4)))
                print(str(sum(1 for x in df_dat.dxyfst if x == 2))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4))+' and Fst values of >='+str(round(out_fst, 4)))
                print(str(sum(1 for x in df_dat.dxydd if x == 2))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4))+' and DD values of <='+str(round(out_dd, 4)))
                print(str(sum(1 for x in df_dat.fstdd if x == 2))+'\toutlier windows for Fst values of >='+str(round(out_fst, 4))+' and DD values of <='+str(round(out_dd, 4)))
                print('\tof which are triple outliers for:')
                print(str(sum(1 for x in df_dat.afddxyfst if x == 3))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+', Dxy values of >='+str(round(out_dxy, 4))+' and Fst values of <='+str(round(out_fst, 4)))
                print(str(sum(1 for x in df_dat.afddxydd if x == 3))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+', Dxy values of >='+str(round(out_dxy, 4))+' and DD values of <='+str(round(out_dd, 4)))
                print(str(sum(1 for x in df_dat.afdfstdd if x == 3))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+', Fst values of >='+str(round(out_fst, 4))+' and DD values of <='+str(round(out_dd, 4)))
                print(str(sum(1 for x in df_dat.afd21dxyfst if x == 3))+'\toutlier windows for AFD('+args.coh2+'-'+args.coh1+') values of <='+str(round(out_afd_low, 4))+', Dxy values of >='+str(round(out_dxy, 4))+' and Fst values of <='+str(round(out_fst, 4)))
#                print(str(sum(1 for x in df_dat.afd21dxydd if x == 3))+'\toutlier windows for AFD_2_1 values of <='+str(round(out_afd_low, 4))+', Dxy values of >='+str(round(out_dxy, 4))+' and DD values of <='+str(round(out_dd, 4)))
#                print(str(sum(1 for x in df_dat.afd21fstdd if x == 3))+'\toutlier windows for AFD_2_1 values of <='+str(round(out_afd_low, 4))+', Fst values of >='+str(round(out_fst, 4))+' and DD values of <='+str(round(out_dd, 4)))
                print(str(sum(1 for x in df_dat.dxyfstdd if x == 3))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4))+', Fst values of >='+str(round(out_fst, 4))+' and DD values of <='+str(round(out_dd, 4)))
                print('\tof which are quadruple outliers for:')
                print(str(sum(1 for x in df_dat.afddxyfstdd if x == 4))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+', Dxy values of >='+str(round(out_dxy, 4))+', Fst values of >='+str(round(out_fst, 4))+' and DD values of <='+str(round(out_dd, 4))+'\n')
#                print(str(sum(1 for x in df_dat.afd21dxyfstdd if x == 4))+'\toutlier windows for AFD_2_1 values of <='+str(round(out_afd_low, 4))+', Dxy values of >='+str(round(out_dxy, 4))+', Fst values of >='+str(round(out_fst, 4))+' and DD values of <='+str(round(out_dd, 4))+'\n')

                print('Select oulier windows\n')

                file_basename = cwd+'/'+outputdir+'/'+args.coh1+args.coh2+'_'+str(args.snps)+'SNPs_'+str(int(10000*args.cut))+'ppm_'
                df_dat.to_csv(file_basename+'allsites'+args.suf+'.csv', index=False)
                # select all windows that are outliers for at least one metric
                df_outlier = df_dat[((df_dat.dxyout != 0) | (df_dat.fstout != 0) | (df_dat.ddout != 0))]
                df_outlier.to_csv(file_basename+'ALL_outliers'+args.suf+'.csv', index=False)
                # write bedfile for gene retrieval and orthologous gene onthology match
                header = ["scaffold", "bedstart", "end"]
                df_outlier.to_csv(file_basename+'ALL_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

                # create inclusive outlier window lists, that is windows that are outliers for 2 or more metrics are listed in the 2 or more list as well as in the single metric lists
                # select AFD_1_2 only outliers
                df_afd12 = df_dat[(df_dat.afdout == 1)]
                if df_afd12.empty:
                    print('No single AFD('+args.coh1+'-'+args.coh2+') outlier windows')
                else:
                    df_afd12.to_csv(file_basename+'afd12_outliers'+args.suf+'.csv', index=False)
                    # write bedfile for gene retrieval and orthologous gene onthology match
                    header = ["scaffold", "bedstart", "end"]
                    df_afd12.to_csv(file_basename+'afd12_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

                # select AFD_2_1 only outliers
                df_afd21 = df_dat[(df_dat.afd21out == 1)]
                if df_afd21.empty:
                    print('No single AFD('+args.coh2+'-'+args.coh1+') outlier windows')
                else:
                    df_afd21.to_csv(file_basename+'afd21_outliers'+args.suf+'.csv', index=False)
                    header = ["scaffold", "bedstart", "end"]
                    df_afd21.to_csv(file_basename+'afd21_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

                # select Dxy only outliers
                df_dxy = df_dat[(df_dat.dxyout == 1)] # exclusive df_dat[((df_dat.dxyout == 1) & (df_dat.fstout == 0) & (df_dat.ddout == 0))]
                if df_dxy.empty:
                    print('No single Dxy outlier windows')
                else:
                    df_dxy.to_csv(file_basename+'Dxy_outliers'+args.suf+'.csv', index=False)
                    # write bedfile for gene retrieval and orthologous gene onthology match
                    header = ["scaffold", "bedstart", "end"]
                    df_dxy.to_csv(file_basename+'Dxy_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

                # select Fst only outliers
                df_fst = df_dat[(df_dat.fstout == 1)]
                if df_fst.empty:
                    print('No single Fst outlier windows')
                else:
                    df_fst.to_csv(file_basename+'Fst_outliers'+args.suf+'.csv', index=False)
                    # write bedfile for gene retrieval and orthologous gene onthology match
                    header = ["scaffold", "bedstart", "end"]
                    df_fst.to_csv(file_basename+'Fst_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)
                    
                # select DD only outliers
                df_dd = df_dat[(df_dat.ddout == 1)]
                if df_dd.empty:
                    print('No single DD outlier windows')
                else:
                    df_dd.to_csv(file_basename+'DD_outliers'+args.suf+'.csv', index=False)
                    # write bedfile for gene retrieval and orthologous gene onthology match
                    header = ["scaffold", "bedstart", "end"]
                    df_dd.to_csv(file_basename+'DD_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

                # select AFD_1_2 and Dxy double outliers
                df_afddxy = df_dat[(df_dat.afddxy == 2) ]
                if df_afddxy.empty:
                    print('No double AFD('+args.coh1+'-'+args.coh2+') Dxy outlier windows')
                else:
                    df_afddxy.to_csv(file_basename+'afd12Dxy_outliers'+args.suf+'.csv', index=False)
                    # write bedfile for gene retrieval and orthologous gene onthology match
                    header = ["scaffold", "bedstart", "end"]
                    df_afddxy.to_csv(file_basename+'afd12Dxy_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)
                    
                # select AFD_1_2 and Fst double outliers
                df_afdfst = df_dat[(df_dat.afdfst == 2) ]
                if df_afdfst.empty:
                    print('No double AFD('+args.coh1+'-'+args.coh2+') Fst outlier windows')
                else:
                    df_afdfst.to_csv(file_basename+'afd12Fst_outliers'+args.suf+'.csv', index=False)
                    # write bedfile for gene retrieval and orthologous gene onthology match
                    header = ["scaffold", "bedstart", "end"]
                    df_afdfst.to_csv(file_basename+'afd12Fst_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)
                    
                # select AFD_1_2 and DD double outliers
                df_afddd = df_dat[(df_dat.afddd == 2) ]
                if df_afddd.empty:
                    print('No double AFD('+args.coh1+'-'+args.coh2+') outlier windows')
                else:
                    df_afddd.to_csv(file_basename+'afd12DD_outliers'+args.suf+'.csv', index=False)
                    # write bedfile for gene retrieval and orthologous gene onthology match
                    header = ["scaffold", "bedstart", "end"]
                    df_afddd.to_csv(file_basename+'afd12DD_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

                # select AFD_2_1 and Dxy double outliers
                df_afd21dxy = df_dat[(df_dat.afd21dxy == 2) ]
                if df_afd21dxy.empty:
                    print('No double AFD('+args.coh2+'-'+args.coh1+') Dxy outlier windows')
                else:
                    df_afd21dxy.to_csv(file_basename+'afd21Dxy_outliers'+args.suf+'.csv', index=False)
                    header = ["scaffold", "bedstart", "end"]
                    df_afd21dxy.to_csv(file_basename+'afd21Dxy_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)
                    
                # select AFD_2_1 and Fst double outliers
                df_afd21fst = df_dat[(df_dat.afd21fst == 2) ]
                if df_afd21fst.empty:
                    print('No double AFD('+args.coh2+'-'+args.coh1+') Fst outlier windows')
                else:
                    df_afd21fst.to_csv(file_basename+'afd21Fst_outliers'+args.suf+'.csv', index=False)
                    header = ["scaffold", "bedstart", "end"]
                    df_afd21fst.to_csv(file_basename+'afd21Fst_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)
                    
                # select Dxy and Fst double outliers
                df_dxyfst = df_dat[(df_dat.dxyfst == 2) ]
                if df_dxyfst.empty:
                    print('No double Dxy Fst outlier windows')
                else:
                    df_dxyfst.to_csv(file_basename+'DxyFst_outliers'+args.suf+'.csv', index=False)
                    # write bedfile for gene retrieval and orthologous gene onthology match
                    header = ["scaffold", "bedstart", "end"]
                    df_dxyfst.to_csv(file_basename+'DxyFst_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

                # select Dxy and DD double outliers
                df_dxydd = df_dat[(df_dat.dxydd == 2)]
                if df_dxydd.empty:
                    print('No double Dxy DD outlier windows')
                else:
                    df_dxydd.to_csv(file_basename+'DxyDD_outliers'+args.suf+'.csv', index=False)
                    # write bedfile for gene retrieval and orthologous gene onthology match
                    header = ["scaffold", "bedstart", "end"]
                    df_dxydd.to_csv(file_basename+'DxyDD_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

                # select Fst and DD double outliers 
                df_fstdd = df_dat[(df_dat.fstdd == 2)]
                if df_fstdd.empty:
                    print('No Fst DD double outlier windows')
                else:
                    df_fstdd.to_csv(file_basename+'FstDD_outliers'+args.suf+'.csv', index=False)
                    # write bedfile for gene retrieval and orthologous gene onthology match
                    header = ["scaffold", "bedstart", "end"]
                    df_fstdd.to_csv(file_basename+'FstDD_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

                # select AFD_1_2, Fst, DD triple outliers 
                df_afddxyfst = df_dat[(df_dat.afddxyfst == 3)]
                if df_afddxyfst.empty:
                    print('No AFD('+args.coh1+'-'+args.coh2+') Dxy Fst triple outlier windows')
                else:
                    df_afddxyfst.to_csv(file_basename+'afd12DxyFst_outliers'+args.suf+'.csv', index=False)
                    # write bedfile for gene retrieval and orthologous gene onthology match
                    header = ["scaffold", "bedstart", "end"]
                    df_afddxyfst.to_csv(file_basename+'afd12DxyFst_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)
                    
                # select AFD_1_2, Dxy, DD triple outliers 
                df_afddxydd = df_dat[(df_dat.afddxydd == 3)]
                if df_afddxydd.empty:
                    print('No AFD('+args.coh1+'-'+args.coh2+') Dxy DD triple outlier windows')
                else:
                    df_afddxydd.to_csv(file_basename+'afd12DxyDD_outliers'+args.suf+'.csv', index=False)
                    # write bedfile for gene retrieval and orthologous gene onthology match
                    header = ["scaffold", "bedstart", "end"]
                    df_afddxydd.to_csv(file_basename+'afd12DxyDD_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)
                    
                # select AFD_1_2, Fst, DD triple outliers 
                df_afdfstdd = df_dat[(df_dat.afdfstdd == 3)]
                if df_afdfstdd.empty:
                    print('No AFD('+args.coh1+'-'+args.coh2+') Fst DD triple outlier windows')
                else:
                    df_afdfstdd.to_csv(file_basename+'afd12FstDD_outliers'+args.suf+'.csv', index=False)
                    # write bedfile for gene retrieval and orthologous gene onthology match
                    header = ["scaffold", "bedstart", "end"]
                    df_afdfstdd.to_csv(file_basename+'afd12FstDD_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

                # select AFD_2_1, Dxy, DD triple outliers 
                df_afd21dxyfst = df_dat[(df_dat.afd21dxyfst == 3)]
                if df_afd21dxyfst.empty:
                    print('No AFD('+args.coh2+'-'+args.coh1+') Dxy Fst triple outlier windows')
                else:
                    df_afd21dxyfst.to_csv(file_basename+'afd21DxyFst_outliers'+args.suf+'.csv', index=False)
                    header = ["scaffold", "bedstart", "end"]
                    df_afd21dxyfst.to_csv(file_basename+'afd21DxyFst_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)
                                        
                # select Dxy, Fst, DD triple outliers 
                df_dxyfstdd = df_dat[(df_dat.dxyfstdd == 3)]
                if df_dxyfstdd.empty:
                    print('No Dxy Fst DD triple outlier windows')
                else:
                    df_dxyfstdd.to_csv(file_basename+'DxyFstDD_outliers'+args.suf+'.csv', index=False)
                    # write bedfile for gene retrieval and orthologous gene onthology match
                    header = ["scaffold", "bedstart", "end"]
                    df_dxyfstdd.to_csv(file_basename+'DxyFstDD_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

                # select AFD12, Dxy, Fst, DD quadruple outliers 
                df_afddxyfstdd = df_dat[(df_dat.afddxyfstdd == 4)]
                if df_afddxyfstdd.empty:
                    print('No AFD('+args.coh1+'-'+args.coh2+') Dxy Fst DD quadruple outlier windows')
                else:
                    df_afddxyfstdd.to_csv(file_basename+'afd12DxyFstDD_outliers'+args.suf+'.csv', index=False)
                    # write bedfile for gene retrieval and orthologous gene onthology match
                    header = ["scaffold", "bedstart", "end"]
                    df_afddxyfstdd.to_csv(file_basename+'afd12DxyFstDD_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)


                # set outliers to be inherited to step 6 to generate graphs
                self.out_afd_up = out_afd_up
                self.out_afd_low = out_afd_low
                self.out_dxy = out_dxy
                self.out_fst = out_fst
                self.out_dd = out_dd
                self.c1_c2 = c1_c2


    def step6(self):
        ###### STEP 6 ######
        ## create histograms with the CI cutoff value for Dxy, Fst, DD and window length
        args = self.args
        if args.o is 'na':
            outputdir = str(args.coh1+args.coh2)
        else:
            outputdir = str(args.o+args.coh1+args.coh2)

        contrast = args.coh1 + args.coh2
        args = self.args
        self.outputdir = outputdir
        self.contrast = contrast
        cwd = self.cwd
        print('\n\tSTEP 6: Create population genetics histograms')
        graphs = []
        for dirName, subdirList, fileList in os.walk(outputdir):
            for file in fileList:
                if file.startswith(contrast) and file.endswith(str(args.snps)+'SNPs_'+str(int(10000*args.cut))+'ppm_allsites'+args.suf+'.csv') == True:
                    graphs.append(dirName+'/'+file)
        graphs_sorted = natsorted(graphs)

        if os.path.exists(outputdir+'/graphs') == False:
            os.mkdir(outputdir+'/graphs')

        print('Creating histograms for Dxy, Fst, DD and window length')
        for file in graphs_sorted:
            count = 0
            rfile = open(outputdir+'/graphs/'+contrast+'_graphs.r', 'w')
            rfile.write('# load table\n'+
                        'test <- read.csv("'+cwd+'/'+file+'", header=T)\n'+
                        '#str(test)\n'+
                        '#load library\n'+
                        'library(ggplot2)\n'+
                        'library(methods)\n'+
                        'p.afd <- qplot('+str(self.c1_c2)+', data=test, binwidth=0.001) + geom_vline(xintercept='+str(self.out_afd_up)+', color="grey30", linetype="dashed") + geom_vline(xintercept='+str(self.out_afd_low)+', color="grey30", linetype="dotted") + theme_bw()\n'+
                        'p.dxy <- qplot(Dxy, data=test, binwidth=0.001, xlim=c(0, (max(test$Dxy)+0.01))) + geom_vline(xintercept='+str(self.out_dxy)+', color="grey30", linetype="dashed") + theme_bw()\n'+
                        'p.fst <- qplot(Fst, data=test, binwidth=0.001, xlim=c(0, (max(test$Fst)+0.01))) + geom_vline(xintercept='+str(self.out_fst)+', color="grey30", linetype=2) + theme_bw()\n'+
                        'p.dd <- qplot(DD, data=test, binwidth=0.001) + geom_vline(xintercept='+str(self.out_dd)+', color="grey30", linetype=2) + theme_bw()\n'+
                        'p.length <- qplot(length, data=test, binwidth=100, xlim=c(0,26560)) + theme_bw()\n\n'+
                        '# export as pdf\n'+
                        'pdf(file="'+outputdir+'/graphs/'+contrast+'_'+str(args.snps)+'SNPs_'+str(int(10000*args.cut))+'ppm_AFD_histogram'+args.suf+'.pdf")\n'+
                        'p.afd\n'+
                        'dev.off()\n'+
                        'pdf(file="'+outputdir+'/graphs/'+contrast+'_'+str(args.snps)+'SNPs_'+str(int(10000*args.cut))+'ppm_Dxy_histogram'+args.suf+'.pdf")\n'+
                        'p.dxy\n'+
                        'dev.off()\n'+
                        'pdf(file="'+outputdir+'/graphs/'+contrast+'_'+str(args.snps)+'SNPs_'+str(int(10000*args.cut))+'ppm_Fst_histogram'+args.suf+'.pdf")\n'+
                        'p.fst\n'+
                        'dev.off()\n'+
                        'pdf(file="'+outputdir+'/graphs/'+contrast+'_'+str(args.snps)+'SNPs_'+str(int(10000*args.cut))+'ppm_DD_histogram'+args.suf+'.pdf")\n'+
                        'p.dd\n'+
                        'dev.off()\n'+
                        'pdf(file="'+outputdir+'/graphs/'+contrast+'_'+str(args.snps)+'SNPs_length_histogram.pdf")\n'+
                        'p.length\n'+
                        'dev.off()')
            rfile.close()

            cmd = ('Rscript '+outputdir+'/graphs/'+args.coh1+args.coh2+'_graphs.r')
            p = subprocess.Popen(cmd, shell=True)
            sts = os.waitpid(p.pid, 0)[1]

            count += 1
            os.remove(outputdir+'/graphs/'+args.coh1+args.coh2+'_graphs.r')
        print('\n\tDONE\n')    
    

    # define the sequence of steps that should be run
    def step1to6(self):
        self.step1()
        self.step2()
        self.step3()
        self.step4()
        self.step5()
        self.step6()
    def step4to6(self):
        self.step4()
        self.step5()
        self.step6()
    def step5to6(self):
        self.step5()
        self.step6()

if __name__ == '__main__':
    import os, sys
