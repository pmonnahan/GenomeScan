## genome scan pipeline, written by Jeff DaCosta, streamlined by Christian Sailer
## September 2016, updated 11 November 2016

import os, sys, argparse, subprocess, statistics
from natsort import natsorted
import pandas as pd
import numpy as np


#create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser(description='Part 2 of genome scan pipeline. This script identifies the genes that lie '+
                                 'within the given interval (-win) using the annotation gff file. In the second step '+
                                 'it greps the gene functions from the ortholog gene function list.')

parser.add_argument('-i', type=str, metavar='inputdir_path', required=True, help='REQUIRED: Full or relative path to contrast directory')
parser.add_argument('-an', type=str, metavar='gff_annotation_file', default='LyV2.gff', help='Full path to annotation gff file')
parser.add_argument('-gf', type=str, metavar='gene_function_file', default='LyV2_TAIR10orth_des_20150927.txt', help='Full path to lyrate - thaliana orthologs file')
parser.add_argument('-ovlp', type=float, metavar='proportion_overlap', default='0.000001', help='Proportion of base pairs that overlap the search pattern as percentage [0.0001]')
args = parser.parse_args()

print('\nSearching directories...')

if args.ovlp == 0.000001:
    if os.path.exists(args.i+'lax_genes/') == False:
        os.mkdir(args.i+'lax_genes/')
    outputdir = str(args.i+'lax_genes/')
else:
    if os.path.exists(args.i+'genes/') == False:
        os.mkdir(args.i+'genes/')
    outputdir = str(args.i+'genes/')


###### STEP 1 ######
# obtain outlier gene annotation list using bedtools
print('\n\tStep 1: Obtain outlier gene annotation list, '+str(args.ovlp)+' proportion overlap\n')
a = []
#files = []
for dirName, subdirList, fileList in os.walk(args.i):
    for file in fileList:
        if file.endswith('_outliers.bed') == True:
            a.append(file)
a.sort()

print('Found '+str(len(a))+' bedfiles to select gene annotation:')
count = 0
for file in a:
    print(' Processing '+str(file))
    basename = file.replace('_outliers.bed', '')
    # gcmd = open(outputdir+'bedtools_gff.unix', 'w')
    # gcmd.write('bedtools intersect -a '+args.i+'/'+str(file)+' -b '+str(args.an)+' -f '+str(args.ovlp)+' -wb | ')
    # gcmd.write('bedtools intersect -a '+args.i+'/'+str(file)+' -b '+str(args.an)+' -f '+str(args.ovlp)+' -wb | ')
    # gcmd.write("""awk '{$1=$2=$3=""; print $4,$5,$6,$7,$8,$9,$10,$11,$12}' """)
    # gcmd.write('| grep transcript | grep -v transcription | sort -u | ')
    # gcmd.write("""tr ' ' '\t' """)
    # gcmd.write('> '+outputdir+basename+'_'+str(args.ovlp)+'ol_genes.gff')
    # gcmd.close()
    shfile1 = open(outputdir+'bedtools_gff.sh','w')
    shfile1.write('#!/bin/bash\n'+
                '#SBATCH -J GS.bedtools.sh'+'\n'+
                '#SBATCH -e GS.bedtools.err\n'+
                '#SBATCH -o GS.bedtools.out\n'+
                '#SBATCH -p nbi-short\n'+
                '#SBATCH -n 1\n'+
                '#SBATCH -t 2-5:00\n'+
                '#SBATCH --mem=16000\n'+
                'source bedtools-2.17.0\n'+
                'bedtools intersect -a '+args.i+str(file)+' -b '+str(args.an)+' -f '+str(args.ovlp)+' -wb | '+
                """awk '{$1=$2=$3=""; print $4,$5,$6,$7,$8,$9,$10,$11,$12}'""" +
                '| grep transcript | grep -v transcription | sort -u | '+
                """tr ' ' '\t' """
                '> '+outputdir+basename+'_'+str(args.ovlp)+'ol_genes.gff')                        
    shfile1.close()

    cmd1 = ('sbatch '+outputdir+'bedtools_gff.sh')
    p1 = subprocess.Popen(cmd1, shell=True)
    sts1 = os.waitpid(p1.pid, 0)[1]

    #run in unix
    # cmd = (open(outputdir+'bedtools_gff.unix', 'r'))
    # p = subprocess.Popen(cmd, shell=True)
    # sts = os.waitpid(p.pid, 0)[1]

    count +=1


                
print('\nExtracted gene annotation lists from '+str(count)+' bedfiles\n')
#os.remove(outputdir+'bedtools_gff.sh')


###### STEP 2 ######
# obtain outlier gene list
print('\n\tStep 2: Obtain outlier gene list, '+str(args.ovlp)+' percent overlap\n')
gff = []
for dirName, subdirList, fileList in os.walk(outputdir):
    for file in fileList:
        if file.endswith(str(args.ovlp)+'ol_genes.gff') == True:
            gff.append(file)
gff.sort()

print(' Found '+str(len(gff))+' gff candidate files to select genes:')
count = 0
for file in gff:
    print('Processing '+str(file))
    basename = file.replace('.gff', '')
    gcmd = open(outputdir+'sed.unix', 'w')
    gcmd.write("""sed -n -e 's/^.*;Parent=//p' """)
    gcmd.write(outputdir+file+' | sort -u > '+outputdir+basename+'.txt')
    gcmd.close()

    #run in unix
    cmd = (open(outputdir+'sed.unix', 'r'))
    p = subprocess.Popen(cmd, shell=True)
    sts = os.waitpid(p.pid, 0)[1]

    count +=1

print('\nExtracted gene lists from '+str(count)+' bedfiles\n')
os.remove(outputdir+'sed.unix')


###### STEP 3 ######
# create interval list file in bed format for candidate genes
print('\n\tStep 3: Create interval bedfiles for candidate genes\n')
inlist = []
for dirName, subdirList, fileList in os.walk(outputdir):
    for file in fileList:
        if file.endswith(str(args.ovlp)+'ol_genes.gff'):
            inlist.append(file)
inlist.sort()
print('Found '+str(len(inlist))+' gff candidate files to obtain intervals:')

for file in inlist:
    print(' Processing '+str(file))
    basename = file.replace('_genes.gff','')
    infile = open(outputdir+file,'r')
    outfile = open(outputdir+basename+'_intervals.bed', 'w')
    for line in infile:
        data = line.split()
        bedstart = int(data[3])-1
        outfile.write(data[0]+'\t'+str(bedstart)+'\t'+data[4]+'\t'+data[8]+'\n')
    infile.close()
    outfile.close()

print('\nCreated '+str(count)+' candidate interval bedfiles\n')


###### Step 4 ######
# obtain gene onthologies from A. thaliana orthologs
print('\n\tSTEP 4: Obtain GO terms for gene lists\n')

a = []
for dirName, subdirList, fileList in os.walk(outputdir):
    for file in fileList:
        if file.endswith(str(args.ovlp)+'ol_genes.txt') == True:
            a.append(file)
a.sort()

print('Found '+str(len(a))+' genelists to select genes:')

count = 0
for file in a:
    print(' Processing '+str(file))
    basename = file.replace('.txt', '')
    outfile = open(outputdir+basename+'_GF.txt', 'w')
    query = open(outputdir+file, 'r')
    for tline in query:
        line = tline.replace('\n', '')
        test = open(args.gf, 'r')
        for testline in test:
            data = testline.split('\t')
            if line in data[0]:
                outfile.write(data[0]+'\t'+data[1]+'\t'+data[3]+'\t'+data[4]+'\t'+data[5]+'\t'+data[6]+'\n')
    query.close()
    outfile.close()
    count +=1

print('\nObtained Arabidopsis thaliana ortholog gene function for '+str(count)+' outlier gene lists\n')
#os.remove(outputdir+'GF_grep.unix')
print('\n\tDONE\n')
