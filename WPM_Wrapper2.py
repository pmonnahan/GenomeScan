#!/usr/bin/env python35

import os
import subprocess
import argparse
import pandas
import math

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


    def splitVCFs(self, vcf_dir, ref_path, mem=16000, numcores=1,print1=False, overwrite=False):
        if vcf_dir.endswith("/") is False:
            vcf_dir += "/"
        outdir = self.split_dir

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

    # SPLIT VCF, CONVERT TO TABLES, REMOVE VCF
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

                    if print1 is True:
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

    # CALCULATE WITHIN POPULATION METRICS
    def calcwpm(self, recode_dir, window_size, min_snps, print1=False, mem=16000, numcores=1):

        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        if os.path.exists(recode_dir) is True:
            summary = open(recode_dir + "WPMInputSummary.txt", 'w')
            summary.write('VCF Directory = ' + self.vcf_dir + '\n' +
                          'Window size in bp = ' + window_size + '\n' +
                          "Minimum number of SNPs per window  = " + min_snps + "\n")
            summary.close()

            for pop in self.pops:
                shfile3 = open(pop + '.sh', 'w')

                shfile3.write('#!/bin/bash\n' +
                              '#SBATCH -J ' + pop + '.sh' + '\n' +
                              '#SBATCH -e ' + self.oande + pop + '.wpm.err' + '\n' +
                              '#SBATCH -o ' + self.oande + pop + '.wpm.out' + '\n' +
                              '#SBATCH -p nbi-medium\n' +
                              '#SBATCH -n ' + str(numcores) + '\n' +
                              '#SBATCH -t 1-00:00\n' +
                              '#SBATCH --mem=' +str(mem) + '\n' +
                              'source python-3.5.1\n' +
                              'source env/bin/activate\n' +
                              'python3 /usr/users/JIC_c1/monnahap/GenomeScan/wpm.py -i ' + recode_dir + pop + '.table.recode.txt -o ' + recode_dir + ' -sampind ' + str(self.samp_ind) + ' -ws ' + str(window_size) + ' -ms ' + str(min_snps) + '\n')
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
            print("Must run splitVCFs followed by recode before able to calculate within population metrics")

    def calcpairwisebpm(self, recode_dir, pop1, pop2, window_size, minimum_snps, print1=False, mem=16000, numcores=1):

        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        if os.path.exists(recode_dir) is True:

            shfile3 = open(pop1 + 'v' + pop2 + '.sh', 'w')

            shfile3.write('#!/bin/bash\n' +
                          '#SBATCH -J ' + pop1 + 'v' + pop2 + '.bpm.sh' + '\n' +
                          '#SBATCH -e ' + self.oande + pop1 + 'v' + pop2 + '.bpm.err' + '\n' +
                          '#SBATCH -o ' + self.oande + pop1 + 'v' + pop2 + '.bpm.out' + '\n' +
                          '#SBATCH -p nbi-medium\n' +
                          '#SBATCH -n ' + str(numcores) + '\n' +
                          '#SBATCH -t 1-00:00\n' +
                          '#SBATCH --mem=' + str(mem) + '\n' +
                          'source python-3.5.1\n' +
                          'source env/bin/activate\n' +
                          'python3 /usr/users/JIC_c1/monnahap/GenomeScan/bpm.py -i1 ' + recode_dir + pop1 + '.table.recode.txt -i1 ' + recode_dir + pop2 + '.table.recode.txt -o ' + recode_dir + ' -ws ' + str(window_size) + ' -ms ' + str(minimum_snps) + '\n')
            shfile3.close()

            if print1 is False:
                cmd3 = ('sbatch -d singleton ' + pop1 + 'v' + pop2 + '.sh')
                p3 = subprocess.Popen(cmd3, shell=True)
                sts3 = os.waitpid(p3.pid, 0)[1]
            else:
                file3 = open(pop1 + 'v' + pop2 + '.sh', 'r')
                data3 = file3.read()
                print(data3)

            os.remove(pop1 + 'v' + pop2 + '.sh')

        else:
            print("Must run splitVCFs followed by recode before able to calculate between population metrics")
