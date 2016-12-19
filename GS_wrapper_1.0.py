#!/usr/bin/env python35

import os, sys, subprocess, argparse,pandas

# create variables that can be entered in the command line
parser = argparse.ArgumentParser()
parser.add_argument('-CF', type=str, metavar='contrast_file', required=True, help='path to csv file containing information on what samples to include')
parser.add_argument('-V', type=str, metavar='vcf_directory', required=False, default="-99", help='directory to vcf files to be split and converted to table')
parser.add_argument('-T', type=str, metavar='table_directory', required=False, default="-99",help='directory to vcf files to be split and converted to table')
parser.add_argument('-GZ', type=str, metavar='gzipped',required=False,default="true",help='are vcfs gzipped (true) or not (false)')
parser.add_argument('-GS', type=str, metavar='path_to_GS_scripts',required=False,default='-99',help='full path to folder containing GenomeScan scripts; required if args.PS is true (default).  Must also contain the annotation and gene function files')
parser.add_argument('-PS', type=str,metavar='perform_GS_scan',required=False,default='true',help='true or false to perform genome scan; if false it just splits vcfs and converts to tables')
parser.add_argument('-OVLP', type=str, metavar='percent_overlap',required=False,default='0.00001',help='Percentage of base pairs that should overlap in the search, default is 0.00001')
parser.add_argument('-SNPS', type=str, metavar='Num_SNPs',required=False,default='100',help='Number of SNPs per analysed window, default is 100SNPs')
parser.add_argument('-WIN_TYPE', type=str, metavar='Window_Type',required=False,default='geo',help='geo for geographic windows; snp for snp-based windows')
parser.add_argument('-WS', type=str, metavar='Window_Size_BP',required=False,default='10000',help='size of windows in terms of bp; for geo window_type')
parser.add_argument('-MS', type=str, metavar='Minimum_SNPs',required=False,default='10',help='Minimum number of SNPs in a window')
parser.add_argument('-M', type=str, metavar='Missingness',required=False,default='0.0',help='Proportion of missing individuals allowed in each group')
parser.add_argument('-P1', type=str, metavar='Print_VCF_Commands', required=False, default='false', help='if true then print shell script to screen')
parser.add_argument('-P2', type=str, metavar='Print_GS_Commands', required=False, default='false', help='if true then print shell script to screen')
parser.add_argument('-metrics', type=str, metavar='metrics',required=True,help='Outlier/Outlier combinations for which to produce graphs, possiblities are: Dxy, Fst, DD, DxyFst, DxyDD, FstDD, DxyFstDD')
parser.add_argument('-PM', type=str, metavar='Percentage_missing_allowed',required=False, default = '0.25',help='Percentage of allowed missing data specifed as a ratio, eg. for 20% type 0.2, default is 0.25')
parser.add_argument('-MEM', type=str, metavar='Memory_Requested',required=True,help='number of Gb of memory to request for all jobs')
parser.add_argument('-K', type=str, metavar='Keep_vcf_files',required=False,default='false',help='Do you want to keep the split vcf files?')
parser.add_argument('-CUT', type=str, metavar='Percentile_cutoff',required=False,default='99.95',help='Percentile for defining outliers; alternative to ci option')
parser.add_argument('-O', type=str, metavar='output_directory',required=True,help='Output Directory')

args = parser.parse_args()

contrast_file=pandas.read_csv(args.CF,header=0)
contrast_names=list(contrast_file.columns.values)[1:]
sample_names=list(contrast_file['Samples'])



mem = args.MEM+'000'

AN ="-99"
for file in os.listdir(args.GS):
	if file == 'LyV2.gff':
		AN = 'LyV2.gff'
if AN == "-99":
	print("Annotation File not found!!")
	args.GZ=-99

GF ="-99"
for file in os.listdir(args.GS):
	if file == 'LyV2_TAIR10orth_des_20150927.txt':
		GF = 'LyV2_TAIR10orth_des_20150927.txt'
if GF == "-99":
	print("Gene Function File not found!!")
	args.GZ=-99

if args.O.endswith("/")==False:
	args.O=args.O+"/"
if args.V!="-99" and args.V.endswith("/")==False:
	args.V=args.V+"/"
if args.T!="-99" and args.T.endswith("/")==False:
	args.T=args.T+"/"

if os.path.exists(args.O) == False:
	os.mkdir(args.O)


for contrast in contrast_names:
	include_index=list(contrast_file[contrast])
	group1=[]
	group2=[]
	sample_string1=""
	sample_string2=""
	ss1=""
	ss2=""
	ss3=""
	if len(contrast.split("_"))==2:
		group1name,group2name=contrast.split("_")
	else:
		group1name="group1"
		group2name="group2"

	for i,sample in enumerate(sample_names):
		if include_index[i]==1:
			group1.append(sample)
			sample_string1 += ' -sn ' + sample
			ss1+=sample+','
		elif include_index[i]==2:
			group2.append(sample)
			sample_string2 += ' -sn ' + sample
			ss2+=sample+','
		elif include_index[i]==0:
			ss3+=sample+','

	if args.T!="-99":
		indir=args.T
	elif args.V!="-99":
		indir=args.O+contrast
	elif args.V != "-99" and args.T != "-99":
		print("Must Provide EITHER directory to VCF's or .table files...not both")
		break
	else:
		print("Must Provide directory to VCF's or .table files")
		break

	if os.path.exists(args.O+contrast) == False:
		os.mkdir(args.O+contrast)
	if os.path.exists(args.O+contrast+"/OandE") == False:
		os.mkdir(args.O+contrast+"/OandE")
	outdir=args.O+contrast+"/"
	OandE=args.O+contrast+"/OandE/"


	summary=open(outdir+"InputSummary_"+contrast+".txt",'w')
	summary.write('Contrast Name = '+contrast+'\n'+
				'VCF Directory = '+args.V+'\n'+
				'Window Type = '+args.WIN_TYPE+'\n')
	if args.WIN_TYPE=="snp":
		summary.write('Number of SNPs per Window = '+args.SNPS+'\n')
	elif args.WIN_TYPE=='geo':
		summary.write('Window size in bp = '+args.WS+'\n')
		summary.write("Minimum number of SNPs per window  = "+args.MS+"\n")

	summary.write('Was GS scan performed? '+args.PS+'\n'+
				'Annotation file = '+AN+'\n'+
				'Gene Function file = '+GF+'\n'+
				'Percentage base pairs overlapping (-ovlp) = '+args.OVLP+'\n'
				'Percentile Cutoff? =' +args.CUT+'\n'
				'Proportion missing data allowed = '+args.PM+'\n'
				'Metrics for which graphs requested = '+args.metrics+'\n'
				'Group 1 = '+ss1+'\n'
				'Group 2 = '+ss2+'\n'
				'Excluded Samples = '+ss3+'\n')
	summary.close()



	joblist=[]
	if args.V!="-99":

		vcf_list = []
		vcf_basenames = []
		for file in os.listdir(args.V):
			if file[-6:] == 'vcf.gz':
				vcf_list.append(file)
				vcf_basenames.append(file[:-7])
			elif file[-3:] == 'vcf':
				vcf_list.append(file)
				vcf_basenames.append(file[:-4])

		for v,vcf in enumerate(vcf_list):

			#THIS MAY HAVE MESSED UP THE FINAL STEP BY RENAMING SHELL SCRIPTS...SBATCH AT END NEEDS TO WAIT FOR ALL JOBS NAMED CONTRAST.SH TO END.
			shfile1 = open('vcf'+str(v)+'.group1.sh','w')
			shfile1.write('#!/bin/bash\n'+
						'#SBATCH -J GS.'+contrast+'.sh'+'\n'+
						'#SBATCH -e '+OandE+contrast+'.'+vcf+'.err\n'+
						'#SBATCH -o '+OandE+contrast+'.'+vcf+'.out\n'+
						'#SBATCH -p nbi-long\n'+
						'#SBATCH -n 1\n'+
						'#SBATCH -t 2-5:00\n'+
						'#SBATCH --mem='+mem+'\n'+
						'source GATK-3.6.0\n'+
						'java -XX:ParallelGCThreads=2 -Xmx'+args.MEM+'g -jar /nbi/software/testing/GATK/3.6.0/src/GenomeAnalysisTK.jar -T SelectVariants -R /nbi/group-data/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta -V ' + args.V + vcf + sample_string1 +' -o ' + outdir + vcf_basenames[v] + '.' + group1name + '.vcf\n'+
						'java -XX:ParallelGCThreads=2 -Xmx'+args.MEM+'g -jar /nbi/software/testing/GATK/3.6.0/src/GenomeAnalysisTK.jar -T VariantsToTable -R /nbi/group-data/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta -V ' + outdir + vcf_basenames[v] + '.' + group1name+'.vcf -F CHROM -F POS -F AC -F AN -F DP -o ' + outdir + vcf_basenames[v] + '.' + group1name + '_raw.table\n')


			shfile2 = open('vcf'+str(v)+'.group2.sh','w')
			shfile2.write('#!/bin/bash\n'+
						'#SBATCH -J GS.'+contrast+'.sh'+'\n'+
						'#SBATCH -e '+OandE+contrast+'.'+vcf+'.err\n'+
						'#SBATCH -o '+OandE+contrast+'.'+vcf+'.out\n'+
						'#SBATCH -p nbi-long\n'+
						'#SBATCH -n 1\n'+
						'#SBATCH -t 2-5:00\n'+
						'#SBATCH --mem='+mem+'\n'+
						'source GATK-3.6.0\n'+
						'java -XX:ParallelGCThreads=2 -Xmx'+args.MEM+'g -jar /nbi/software/testing/GATK/3.6.0/src/GenomeAnalysisTK.jar -T SelectVariants -R /nbi/group-data/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta -V ' + args.V + vcf + sample_string2 +' -o ' + outdir + vcf_basenames[v] + '.' + group2name + '.vcf\n'+
						'java -XX:ParallelGCThreads=2 -Xmx'+args.MEM+'g -jar /nbi/software/testing/GATK/3.6.0/src/GenomeAnalysisTK.jar -T VariantsToTable -R /nbi/group-data/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta -V ' + outdir + vcf_basenames[v] + '.' + group2name + '.vcf -F CHROM -F POS -F AC -F AN -F DP -o ' + outdir + vcf_basenames[v] + '.' + group2name + '_raw.table\n')
						

			if args.K == 'false':
				shfile1.write('rm '+ outdir + vcf_basenames[v] + '.' + group1name + '.vcf\n'+
							'rm '+ outdir + vcf_basenames[v] + '.' + group1name + '.vcf.idx')						 
				shfile2.write('rm '+ outdir + vcf_basenames[v] + '.' + group2name + '.vcf\n'+
							'rm '+ outdir + vcf_basenames[v] + '.' + group2name + '.vcf.idx')
				shfile1.close()
				shfile2.close()
			else:
				shfile1.write('gzip '+ outdir + vcf_basenames[v] + '.' + group1name + '.vcf')
				shfile2.write('gzip '+ outdir + vcf_basenames[v] + '.' + group2name + '.vcf')
				shfile1.close()
				shfile2.close()


			if args.P1=='false': #send slurm job to NBI SLURM cluster
				cmd1 = ('sbatch vcf'+str(v)+'.group1.sh')
				p1 = subprocess.Popen(cmd1, shell=True)
				sts1 = os.waitpid(p1.pid, 0)[1]
				joblist.append(p1.pid)

				cmd2 = ('sbatch vcf'+str(v)+'.group2.sh')
				p2 = subprocess.Popen(cmd2, shell=True)
				sts2 = os.waitpid(p2.pid, 0)[1]
				joblist.append(p2.pid)

			else:
				print('vcf'+str(v)+'.group1.sh')
				file1 = open('vcf'+str(v)+'.group1.sh','r')
				data1 = file1.read()
				print(data1)
				file2 = open('vcf'+str(v)+'.group2.sh','r')
				data2 = file2.read()
				print(data2)

			os.remove('vcf'+str(v)+'.group1.sh')
			os.remove('vcf'+str(v)+'.group2.sh')

	#PERFORM GENOME SCAN
	if args.PS == 'true':
		if args.GS!='-99':
			shfile3 = open('contrast.sh','w')

			shfile3.write('#!/bin/bash\n'+
							'#SBATCH -J GS.'+contrast+'.sh'+'\n'+
							'#SBATCH -e '+OandE+'GS.'+contrast+'.err'+'\n'+
							'#SBATCH -o '+OandE+'GS.'+contrast+'.out'+'\n'+
							'#SBATCH -p nbi-long\n'+
							'#SBATCH -n 1\n'+
							'#SBATCH -t 2-5:00\n'+
							'#SBATCH --mem='+mem+'\n'+
							'source python-3.5.1\n'+
							'source R-3.2.3\n'+
							'virtualenv --system-site-packages env\n'+
							'source env/bin/activate\n'+
							'python3 ' + args.GS + 'G1_outliers_slurm_new_PJM.py -i '+indir+ ' -coh1 ' + group1name + ' -coh2 ' + group2name + ' -snps '+args.SNPS+' -per '+ args.PM+ ' -win_type ' + args.WIN_TYPE+ ' -m ' +args.M+ ' -ms ' +args.MS+ ' -ws '+ args.WS+ ' -cut ' +args.CUT+ ' -o '+outdir)
			shfile3.close()

			if args.P2 == 'false':
				cmd3 = ('sbatch -d singleton contrast.sh')
				p3 = subprocess.Popen(cmd3, shell=True)
				sts3 = os.waitpid(p3.pid, 0)[1]
			elif args.P2 == 'true':
				file3 = open('contrast.sh','r')
				data3 = file3.read()
				print(data3)

			os.remove('contrast.sh')

		else:
			print("Must provide full path to location of GenomeScan scripts")
		




		



