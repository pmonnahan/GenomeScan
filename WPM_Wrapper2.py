#!/usr/bin/env python35

import os
import subprocess
import argparse
import pandas
import math
import statistics
import numpy as np

# Directions:
# import WPM_Wrapper2 into python console
# and pass path of Individual_Key file as object to class PopGen

# USE BELOW COMMANDS TO CONCAT POP FILES after all are finished running
# head -1 BEL_WPM.txt > WPM_All.txt
# tail -n+2 -q *_WPM.txt >> WPM_All.txt

# create variables that can be entered in the command line


class PopGen:

    def __init__(self, WorkingDir):
        if WorkingDir.endswith("/") is False:
            WorkingDir += "/"
        if os.path.exists(WorkingDir) is False:
            os.mkdir(WorkingDir)
            os.mkdir(WorkingDir + "OandE/")
        if os.path.exists(WorkingDir + "OandE/") is False:
            os.mkdir(WorkingDir + "OandE/")
        POP_file = pandas.read_csv(WorkingDir + "PopKey.csv", header=0)
        POP_names = list(POP_file.columns.values)[1:]
        sample_names = list(POP_file['Samples'])
        samps = {}
        samp_nums = {}
        for pop in POP_names:
            pop_list = []
            include_index = list(POP_file[pop])
            for i, sample in enumerate(sample_names):
                if include_index[i] == 1:
                    pop_list.append(sample)
            samps[pop] = pop_list
            samp_nums[pop] = len(pop_list)


        # Determine number of individuals to downsample all populations to
        min_ind = min([sum(list(POP_file[pop])) for pop in POP_names])
        self.pops = POP_names
        self.samps = samps
        self.samp_nums = samp_nums
        self.min_ind = min_ind
        self.dir = WorkingDir
        self.oande = WorkingDir + "OandE/"
        self.split_dir = self.dir + "VCFs/"


    def removePop(self, popname):
        popname = str(popname)
        if popname in self.pops:
            self.pops.remove(popname)
            self.samps.pop(popname, None)
            self.samp_nums.pop(popname, None)
            # Recalculate min_ind
            min_ind = min([self.samp_nums[pop] for pop in self.pops])
            self.min_ind = min_ind
        else:
            print("Population does not exist")


    def splitVCFs(self, vcf_dir, ref_path, mem=16000, numcores=1, print1=False, overwrite=False):
        if vcf_dir.endswith("/") is False:
            vcf_dir += "/"
        outdir = self.split_dir
        self.vcf_dir = vcf_dir

        mem1 = int(mem / 1000)

        if os.path.exists(outdir) is False:
            os.mkdir(outdir)
        elif overwrite is False:
            print("VCF directory already exists.  Set 'overwrite = True' if you want to overwrite existing files")
        else:
            print("Overwriting files in existing VCF directory")

        for pop in self.pops:
            # Add samples to list for each population according to PF file
            sample_string1 = ""
            for samp in self.samps[pop]:
                sample_string1 += " -sn " + samp
            joblist = []

            vcf_list = []
            vcf_basenames = []
            for file in os.listdir(vcf_dir):
                if file[-6:] == 'vcf.gz':
                    vcf_list.append(file)
                    vcf_basenames.append(file[:-7])
                elif file[-3:] == 'vcf':
                    vcf_list.append(file)
                    vcf_basenames.append(file[:-4])
            for v, vcf in enumerate(vcf_list):
                # Select single population and biallelic SNPs for each scaffold and convert to variants table
                shfile1 = open(pop + '.sh', 'w')
                shfile1.write('#!/bin/bash\n' +
                              '#SBATCH -J ' + pop + '.sh' + '\n' +
                              '#SBATCH -e ' + self.oande + pop + vcf + '.gatk.err' + '\n' +
                              '#SBATCH -o ' + self.oande + pop + vcf + '.gatk.out' + '\n' +
                              '#SBATCH -p nbi-medium\n' +
                              '#SBATCH -n ' + str(numcores) + '\n' +
                              '#SBATCH -t 0-4:00\n' +
                              '#SBATCH --mem=' + str(mem) + '\n' +
                              'source GATK-3.6.0\n' +
                              'java -Xmx' + str(mem1) + 'g -jar /nbi/software/testing/GATK/3.6.0/src/GenomeAnalysisTK.jar -T SelectVariants -R ' + ref_path + ' -V ' + vcf_dir + vcf + sample_string1 + ' -o ' + outdir + vcf_basenames[v] + '.' + pop + '.vcf\n' +
                              'java -Xmx' + str(mem1) + 'g -jar /nbi/software/testing/GATK/3.6.0/src/GenomeAnalysisTK.jar -T SelectVariants -R ' + ref_path + ' -V ' + outdir + vcf_basenames[v] + '.' + pop + '.vcf --restrictAllelesTo BIALLELIC -env -o ' + outdir + vcf_basenames[v] + '.' + pop + '.bi.vcf\n' +
                              'java -Xmx' + str(mem1) + 'g -jar /nbi/software/testing/GATK/3.6.0/src/GenomeAnalysisTK.jar -T VariantsToTable -R ' + ref_path + ' -V ' + outdir + vcf_basenames[v] + '.' + pop + '.bi.vcf -F CHROM -F POS -F AC -F AN -F DP -GF GT -o ' + outdir + vcf_basenames[v] + '.' + pop + '_raw.table\n'
                              'gzip ' + outdir + vcf_basenames[v] + '.' + pop + '.vcf')
                shfile1.close()

                if print1 is False:  # send slurm job to NBI SLURM cluster
                    cmd1 = ('sbatch ' + pop + '.sh')
                    p1 = subprocess.Popen(cmd1, shell=True)
                    sts1 = os.waitpid(p1.pid, 0)[1]
                    joblist.append(p1.pid)

                else:
                    file1 = open(pop + '.sh', 'r')
                    data1 = file1.read()
                    print(data1)

                os.remove(pop + '.sh')

            # combine all variants table for each scaffold within a population
            shfile3 = open(pop + '.sh', 'w')

            shfile3.write('#!/bin/bash\n' +
                          '#SBATCH -J ' + pop + '.sh' + '\n' +
                          '#SBATCH -e ' + self.oande + pop + '.cat.err' + '\n' +
                          '#SBATCH -o ' + self.oande + pop + '.cat.out' + '\n' +
                          '#SBATCH -p nbi-medium\n' +
                          '#SBATCH -n ' + str(numcores) + '\n' +
                          '#SBATCH -t 0-12:00\n' +
                          '#SBATCH --mem=' + str(mem) + '\n' +
                          'cat ' + outdir + '*' + pop + '_raw.table | tail -n+2 > ' + outdir + pop + '.table\n')
            shfile3.close()

            if print1 is False:
                cmd3 = ('sbatch -d singleton ' + pop + '.sh')
                p3 = subprocess.Popen(cmd3, shell=True)
                sts3 = os.waitpid(p3.pid, 0)[1]
            else:
                file3 = open(pop + '.sh', 'r')
                data3 = file3.read()
                print(data3)

            os.remove(pop + '.sh')


    def recode(self, min_avg_dp, missingness, print1=False, mem=16000, numcores=1):

        sampind = int(math.ceil(self.min_ind * (1.0 - missingness)))
        if sampind == self.min_ind and missingness != 0.0:
            sampind = self.min_ind - 1
        self.samp_ind = sampind
        # Determine if the recoded vcf files already exist and if so, set VCF_Parse to False
        recode_dir = self.dir + "Recoded.DP" + str(min_avg_dp) + ".M" + str(missingness) + "/"
        if os.path.exists(self.split_dir) is True:
            existing_files = []
            if os.path.exists(recode_dir) is False:
                os.mkdir(recode_dir)
            else:
                for file in os.listdir(recode_dir):
                    if file.endswith('.table.recode.txt'):
                        existing_files.append(file.split('.')[0])
            if set(self.pops).issuperset(set(existing_files)) is True and len(existing_files) != 0:
                print("Recoded vcf files already exist.  Delete folder or change parameters")
            # Look for '.recode' table files corresponding to POP_names...if they exist...skip all but wpm step and print message to output.

            else:
                for pop in self.pops:
                    # FORMAT TABLE FOR WPM FILE.

                    shfile3 = open(pop + '.sh', 'w')
                    shfile3.write('#!/bin/bash\n' +
                                  '#SBATCH -J ' + pop + '.sh' + '\n' +
                                  '#SBATCH -e ' + self.oande + pop + '.recode012.err' + '\n' +
                                  '#SBATCH -o ' + self.oande + pop + '.recode012.out' + '\n' +
                                  '#SBATCH -p nbi-medium\n' +
                                  '#SBATCH -n ' + str(numcores) + '\n' +
                                  '#SBATCH -t 1-00:00\n' +
                                  '#SBATCH --mem=' + str(mem) + '\n' +
                                  'source python-3.5.1\n' +
                                  'source env/bin/activate\n' +
                                  'python3 /usr/users/JIC_c1/monnahap/GenomeScan/recode012.py -i ' + self.split_dir + pop + '.table -pop ' + pop + ' -mf ' + str(1.0 - missingness) + ' -dp ' + str(min_avg_dp) + ' -o ' + recode_dir + '\n')
                    shfile3.close()

                    if print1 is False:
                        cmd3 = ('sbatch -d singleton ' + pop + '.sh')
                        p3 = subprocess.Popen(cmd3, shell=True)
                        sts3 = os.waitpid(p3.pid, 0)[1]
                    else:
                        file3 = open(pop + '.sh', 'r')
                        data3 = file3.read()
                        print(data3)

                    os.remove(pop + '.sh')
        else:
            print("Must run splitVCFs prior to running recode")

    def getPloidies(self, recode_dir):

        print("Be sure that 'recode' scripts have all finished")
        if recode_dir.endswith("/") is False:
            recode_dir += "/"
        ploidies = {}
        dips = []
        tets = []
        if os.path.exists(recode_dir) is True:
            for pop in self.pops:
                try:
                    tmp = open(recode_dir + pop + '.table.recode.txt', 'r')
                    line = tmp.readline()
                    ploidy = line.split("\t")[1]
                    ploidies[pop] = ploidy
                    if ploidy == "4.0":
                        tets.append(pop)
                    elif ploidy == "2.0":
                        dips.append(pop)
                    else:
                        print("Ploidy level not recognized")
                except (FileNotFoundError, IndexError):
                    print("Error determining ploidy for population: ", pop)
            self.ploidies = ploidies
            self.dips = dips
            self.tets = tets
        else:
            print("recode_dir does not exist")

    # CALCULATE WITHIN POPULATION METRICS
    def calcwpm(self, recode_dir, window_size, min_snps, population="all", print1=False, mem=16000, numcores=1, sampind="-99"):
        if sampind == "-99":
            sind = self.samp_ind
        else:
            sind = sampind
        if population == "all":
            pops = self.pops
        else:
            pops = [population]
        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        if os.path.exists(recode_dir) is True:
            summary = open(recode_dir + "WPMInputSummary.txt", 'w')
            summary.write('VCF Directory = ' + self.vcf_dir + '\n' +
                          'Window size in bp = ' + str(window_size) + '\n' +
                          "Minimum number of SNPs per window  = " + str(min_snps) + "\n")
            summary.close()

            for pop in pops:
                shfile3 = open(pop + '.sh', 'w')

                shfile3.write('#!/bin/bash\n' +
                              '#SBATCH -J ' + pop + '.sh' + '\n' +
                              '#SBATCH -e ' + self.oande + pop + '.wpm.err' + '\n' +
                              '#SBATCH -o ' + self.oande + pop + '.wpm.out' + '\n' +
                              '#SBATCH -p nbi-medium\n' +
                              '#SBATCH -n ' + str(numcores) + '\n' +
                              '#SBATCH -t 1-00:00\n' +
                              '#SBATCH --mem=' + str(mem) + '\n' +
                              'source python-3.5.1\n' +
                              'source env/bin/activate\n' +
                              'python3 /usr/users/JIC_c1/monnahap/GenomeScan/wpm.py -i ' + recode_dir + pop + '.table.recode.txt -o ' + recode_dir + ' -sampind ' + str(sind) + ' -ws ' + str(window_size) + ' -ms ' + str(min_snps) + '\n')
                shfile3.close()

                if print1 is False:
                    cmd3 = ('sbatch -d singleton ' + pop + '.sh')
                    p3 = subprocess.Popen(cmd3, shell=True)
                    sts3 = os.waitpid(p3.pid, 0)[1]
                else:
                    file3 = open(pop + '.sh', 'r')
                    data3 = file3.read()
                    print(data3)

                os.remove(pop + '.sh')
        else:
            print("Did not find recode_dir.  Must run splitVCFs followed by recode before able to calculate within population metrics")

    def calcbpm(self, recode_dir, pops, output_name, window_size, minimum_snps, print1=False, mem=16000, numcores=1):

        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        if os.path.exists(recode_dir) is True:

            # Concatenate input files and sort them
            print("Concatenating input files")
            concat_file = open(recode_dir + output_name + '.table.recode.txt', 'w')
            for pop in pops:  # Add data from all populations to single, huge list
                try:
                    with open(recode_dir + pop + '.table.recode.txt', 'r') as in1:
                        for line in in1:
                            concat_file.write(line)
                except IOError:
                    print("Did not find input file for pop ", pop)
            print("Finished preparing input data")

            shfile3 = open(output_name + '.bpm.sh', 'w')

            shfile3.write('#!/bin/bash\n' +
                          '#SBATCH -J ' + output_name + '.bpm.sh' + '\n' +
                          '#SBATCH -e ' + self.oande + output_name + '.bpm.err' + '\n' +
                          '#SBATCH -o ' + self.oande + output_name + '.bpm.out' + '\n' +
                          '#SBATCH -p nbi-medium\n' +
                          '#SBATCH -n ' + str(numcores) + '\n' +
                          '#SBATCH -t 1-00:00\n' +
                          '#SBATCH --mem=' + str(mem) + '\n' +
                          'source python-3.5.1\n' +
                          'source env/bin/activate\n' +
                          'python3 /usr/users/JIC_c1/monnahap/GenomeScan/bpm.py -i ' + recode_dir + output_name + '.table.recode.txt -o ' + recode_dir + ' -prefix ' + output_name + ' -ws ' + str(window_size) + ' -ms ' + str(minimum_snps) + '\n')
            shfile3.close()

            if print1 is False:
                cmd3 = ('sbatch -d singleton ' + output_name + '.sh')
                p3 = subprocess.Popen(cmd3, shell=True)
                sts3 = os.waitpid(p3.pid, 0)[1]
            else:
                file3 = open(output_name + '.sh', 'r')
                data3 = file3.read()
                print(data3)

            os.remove(output_name + '.sh')

        else:
            print("Did not find recode_dir.  Must run splitVCFs followed by recode before able to calculate between population metrics")


    def findOutliers(self, recode_dir, in_file, column_index_list, percentile, metrics, tails='upper'):

        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        if os.path.exists(recode_dir) is True:

            data = pandas.read_table(recode_dir + in_file, header=0)
            metrics = []
            for i in column_index_list:
                metrics.append(list(data.columns.values)[i])
            data.sort(metrics, ascending=[1 for x in range(0, len(metrics))])

            for metric in metrics:
                data[metric + '.out'] = 0
                data = data.sort([metric], ascending=[0])
                if tails == 'both':
                    data[metric + '.out'].loc[(data[metric] > data.quantile(q=percentile, axis=1))] = 1
                    data[metric + '.out'].loc[(data[metric] < data.quantile(q=1.0 - percentile, axis=1))] = 1
                elif tails == 'lower':
                    data[metric + '.out'].loc[(data[metric] < data.quantile(q=1.0 - percentile, axis=1))] = 1
                elif tails == 'upper':
                    data[metric + '.out'].loc[(data[metric] > data.quantile(q=percentile, axis=1))] = 1
                else:
                    print("Did not specify tails option correctly.  Options are: both, upper, and lower")
            data['num_outliers'] = data.iloc[:, -len(metrics):].sum(1)
            data.to_csv(recode_dir + in_file.replace(".txt", "") + '_' + str(percentile) + 'tile_OutLabelled.csv', index=False)
            # select all windows that are outliers for at least one metric
            df_outlier = data[(data.num_outliers != 0)]
            df_outlier.to_csv(recode_dir + in_file.replace(".txt", "") + '_' + str(percentile) + 'tile_OutOnly.csv', index=False)
            df_outlier.to_csv(recode_dir + in_file.replace(".txt", "") + '_' + str(percentile) + 'tile_OutOnly.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)


    def annotateOutliers(self, recode_dir, in_file, basename, annotation_file, overlap_proportion=0.000001):

        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        if os.path.exists(recode_dir) is True:
            shfile1 = open(recode_dir + 'bedtools_gff.sh', 'w')
            shfile1.write('#!/bin/bash\n' +
                          '#SBATCH -J GS.bedtools.sh' + '\n' +
                          '#SBATCH -e GS.bedtools.err\n' +
                          '#SBATCH -o GS.bedtools.out\n' +
                          '#SBATCH -p nbi-short\n' +
                          '#SBATCH -n 1\n' +
                          '#SBATCH -t 0-12:00\n' +
                          '#SBATCH --mem=16000\n' +
                          'source bedtools-2.17.0\n' +
                          'bedtools intersect -a ' + recode_dir + in_file + ' -b ' + annotation_file + ' -f ' + str(overlap_proportion) + ' -wb | ' +
                          """awk '{$1=$2=$3=""; print $4,$5,$6,$7,$8,$9,$10,$11,$12}'""" +
                          '| grep transcript | grep -v transcription | sort -u | ' +
                          """tr ' ' '\t' """
                          '> ' + recode_dir + basename + '_' + str(overlap_proportion * 100) + 'ol_genes.gff')
            shfile1.close()

            cmd1 = ('sbatch ' + recode_dir + 'bedtools_gff.sh')
            p1 = subprocess.Popen(cmd1, shell=True)
            sts1 = os.waitpid(p1.pid, 0)[1]

        else:
            print("recode_dir not found")


    def mergeAnnotation(self, recode_dir, outlier_file, annotated_outlier_file):
        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        if os.path.exists(recode_dir) is True:
            try:
                outliers = pandas.read_table(recode_dir + outlier_file, header=0)
                annotation = pandas.read_table(recode_dir + annotated_outlier_file, names=["scaffold", "start", "end", "info1", "info2", "info3", "info4", "info5", "info6", "info7", "info8", "info9"])
            except IOError:
                print("Did not find either original outlier file or the annotated outlier file")
            merged = pandas.merge(outliers, annotation, ["scaffold", "start", "end"],)
            merged.to_csv(recode_dir + outlier_file.replace(".txt", "") + '_OutAnnot.csv', index=False)

