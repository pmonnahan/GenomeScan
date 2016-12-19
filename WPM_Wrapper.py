#!/usr/bin/env python35

import os, sys, subprocess, argparse,pandas

# create variables that can be entered in the command line
parser = argparse.ArgumentParser()
parser.add_argument('-PF', type=str, metavar='Individual_key', required=True, help='path to csv file containing information on what samples to include')
parser.add_argument('-V', type=str, metavar='vcf_directory', required=True, default="-99", help='directory to vcf files to be split and converted to table')
parser.add_argument('-WS', type=str, metavar='Window_Size_BP',required=False,default='10000',help='size of windows in terms of bp; for geo window_type')
parser.add_argument('-MS', type=str, metavar='Minimum_SNPs',required=False,default='10',help='Minimum number of SNPs in a window')
parser.add_argument('-DP', type=str, metavar='Minimum_Depth',required=False,default='10',help='minimum average depth per individual')
parser.add_argument('-M', type=float, metavar='Missingness',required=False,default='0.0',help='Proportion of missing individuals allowed in each group')
parser.add_argument('-P1', type=str, metavar='Print_VCF_Commands', required=False, default='false', help='if true then print shell script to screen')
parser.add_argument('-P2', type=str, metavar='Print_GS_Commands', required=False, default='false', help='if true then print shell script to screen')
parser.add_argument('-K', type=str, metavar='Keep_vcf_files',required=False,default='false',help='Do you want to keep the split vcf files?')
parser.add_argument('-O', type=str, metavar='output_directory',required=True,help='Output Directory')

args = parser.parse_args()

POP_file=pandas.read_csv(args.PF,header=0)
POP_names=list(POP_file.columns.values)[1:]
sample_names=list(POP_file['Samples'])


if args.O.endswith("/")==False:
	args.O=args.O+"/"
if args.V!="-99" and args.V.endswith("/")==False:
	args.V=args.V+"/"

if os.path.exists(args.O) == False:
	os.mkdir(args.O)
if os.path.exists(args.O+"Within-Population-Metrics/") == False:
	os.mkdir(args.O+"Within-Population-Metrics/")
if os.path.exists(args.O+"Within-Population-Metrics/OandE/") == False:
	os.mkdir(args.O+"Within-Population-Metrics/OandE/")

for pop in POP_names:
	include_index=list(POP_file[pop])
	group1=[]
	sample_string1=""
	ss1=""

	for i,sample in enumerate(sample_names):
		if include_index[i]==1:
			group1.append(sample)
			sample_string1 += ' -sn ' + sample
			ss1+=sample+','

	if os.path.exists(args.O+"Within-Population-Metrics/"+pop+"/") == False:
		os.mkdir(args.O+"Within-Population-Metrics/"+pop+"/")

	outdir1=args.O+"Within-Population-Metrics/"
	outdir=args.O+"Within-Population-Metrics/"+pop+"/"
	oande=args.O+"Within-Population-Metrics/OandE/"

	summary=open(outdir1+"InputSummary.txt",'w')
	summary.write('VCF Directory = '+args.V+'\n'+
				'Window size in bp = '+args.WS+'\n'+
				"Minimum number of SNPs per window  = "+args.MS+"\n"+
				"Minimum average depth per individual  = "+args.DP+"\n"+
				'Proportion missing data allowed = '+str(args.M)+'\n')
	summary.close()

	joblist=[]

	vcf_list = []
	vcf_basenames = []
	for file in os.listdir(args.V):
		if file[-6:] == 'vcf.gz':
			vcf_list.append(file)
			vcf_basenames.append(file[:-7])
		elif file[-3:] == 'vcf':
			vcf_list.append(file)
			vcf_basenames.append(file[:-4])


#SPLIT VCF, CONVERT TO TABLES, REMOVE VCF
	for v,vcf in enumerate(vcf_list):

		shfile1 = open(pop+'.sh','w')
		shfile1.write('#!/bin/bash\n'+
					'#SBATCH -J ' +pop +'.sh'+'\n'+
					'#SBATCH -e '+oande+pop+'.gatk.err'+'\n'+
					'#SBATCH -o '+oande+pop+'.gatk.out'+'\n'+
					'#SBATCH -p nbi-long\n'+
					'#SBATCH -n 1\n'+
					'#SBATCH -t 0-4:00\n'+
					'#SBATCH --mem=16000\n'+
					'source GATK-3.6.0\n'+
					'java -Xmx16g -jar /nbi/software/testing/GATK/3.6.0/src/GenomeAnalysisTK.jar -T SelectVariants -R /nbi/group-data/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta -V ' + args.V + vcf + sample_string1 +' -o ' + outdir + vcf_basenames[v] + '.' + pop + '.vcf\n'+
					'java -Xmx16g -jar /nbi/software/testing/GATK/3.6.0/src/GenomeAnalysisTK.jar -T SelectVariants -R /nbi/group-data/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta -V ' + outdir + vcf_basenames[v] + '.' + pop + '.vcf --restrictAllelesTo BIALLELIC -env -o ' + outdir + vcf_basenames[v] + '.' + pop + '.bi.vcf\n'+
					'java -Xmx16g -jar /nbi/software/testing/GATK/3.6.0/src/GenomeAnalysisTK.jar -T VariantsToTable -R /nbi/group-data/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta -V ' + outdir + vcf_basenames[v] + '.' + pop+'.bi.vcf -F CHROM -F POS -F AC -F AN -F DP -GF GT -o ' + outdir + vcf_basenames[v] + '.' + pop + '_raw.table\n')


		if args.K == 'false':
			shfile1.write('rm '+ outdir + vcf_basenames[v] + '.' + pop + '.vcf\n'+
						'rm '+ outdir + vcf_basenames[v] + '.' + pop + '.bi.vcf\n'+
						'rm '+ outdir + vcf_basenames[v] + '.' + pop + '.vcf.idx\n'+
						'rm '+ outdir + vcf_basenames[v] + '.' + pop + '.bi.vcf.idx')						 
			shfile1.close()
		else:
			shfile1.write('gzip '+ outdir + vcf_basenames[v] + '.' + pop + '.vcf')
			shfile1.close()


		if args.P1=='false': #send slurm job to NBI SLURM cluster
			cmd1 = ('sbatch ' +pop+ '.sh')
			p1 = subprocess.Popen(cmd1, shell=True)
			sts1 = os.waitpid(p1.pid, 0)[1]
			joblist.append(p1.pid)

		else:
			file1 = open(pop+'.sh','r')
			data1 = file1.read()
			print(data1)

		os.remove(pop+'.sh')

#CONCATENATE THE INDIVIDUAL SCAFFOLD TABLES TO A SINGLE TABLE IN THE PARENT DIRECTORY
#THEN, REMOVE THE DIRECTORY CONTAINING SCAFFOLD FILES		

	shfile3 = open(pop+'.sh','w')

	shfile3.write('#!/bin/bash\n'+
					'#SBATCH -J '+pop+'.sh'+'\n'+
					'#SBATCH -e '+oande+pop+'.cat.err'+'\n'+
					'#SBATCH -o '+oande+pop+'.cat.out'+'\n'+
					'#SBATCH -p nbi-long\n'+
					'#SBATCH -n 1\n'+
					'#SBATCH -t 0-12:00\n'+
					'#SBATCH --mem=32000\n'+
					'cat ' + outdir + '*'+ pop + '_raw.table > ' + outdir1 + pop + '.table\n'+
					'rm -r ' + outdir + '\n')
	shfile3.close()

	if args.P2 == 'false':
		cmd3 = ('sbatch -d singleton ' + pop+ '.sh')
		p3 = subprocess.Popen(cmd3, shell=True)
		sts3 = os.waitpid(p3.pid, 0)[1]
	elif args.P2 == 'true':
		file3 = open(pop+'.sh','r')
		data3 = file3.read()
		print(data3)

	os.remove(pop+'.sh')

#FORMAT TABLE FOR WPM FILE.
#not done

	shfile3 = open(pop+'.sh','w')

	shfile3.write('#!/bin/bash\n'+
					'#SBATCH -J '+pop+'.sh'+'\n'+
					'#SBATCH -e '+oande+pop+'.recode012.err'+'\n'+
					'#SBATCH -o '+oande+pop+'.recode012.out'+'\n'+
					'#SBATCH -p nbi-long\n'+
					'#SBATCH -n 1\n'+
					'#SBATCH -t 1-00:00\n'+
					'#SBATCH --mem=32000\n'+
					'source python-3.5.1\n'+
					'source env/bin/activate\n'+
					'python3 /nbi/Research-Groups/JIC/Levi-Yant/Patrick/code/recode012.py -i ' + outdir1 + pop + '.table -mf ' + str(1.0 - args.M) + ' -dp ' + args.DP + ' -o ' + outdir1 +'\n'
					'rm ' + outdir1 + pop + '.table')
	shfile3.close()

	if args.P2 == 'false':
		cmd3 = ('sbatch -d singleton ' + pop+ '.sh')
		p3 = subprocess.Popen(cmd3, shell=True)
		sts3 = os.waitpid(p3.pid, 0)[1]
	elif args.P2 == 'true':
		file3 = open(pop+'.sh','r')
		data3 = file3.read()
		print(data3)

	os.remove(pop+'.sh')

#CALCULATE WITHIN POPULATION METRICS
#not done
	shfile3 = open(pop+'.sh','w')

	shfile3.write('#!/bin/bash\n'+
					'#SBATCH -J '+pop+'.sh'+'\n'+
					'#SBATCH -e '+oande+pop+'.wpm.err'+'\n'+
					'#SBATCH -o '+oande+pop+'.wpm.out'+'\n'+
					'#SBATCH -p nbi-long\n'+
					'#SBATCH -n 1\n'+
					'#SBATCH -t 1-00:00\n'+
					'#SBATCH --mem=32000\n'+
					'source python-3.5.1\n'+
					'source env/bin/activate\n'+
					'python3 /nbi/Research-Groups/JIC/Levi-Yant/Patrick/code/wpm.py -i '+ outdir1 + pop + '.recode.txt -pop ' + pop + ' -o ' + outdir1 + ' -m ' + str(args.M) + ' -ws ' + args.WS + ' -ms ' + args.MS + '\n'+
					'rm ' + outdir1 + pop + '.recode.txt')
	shfile3.close()

	if args.P2 == 'false':
		cmd3 = ('sbatch -d singleton ' + pop + '.sh')
		p3 = subprocess.Popen(cmd3, shell=True)
		sts3 = os.waitpid(p3.pid, 0)[1]
	elif args.P2 == 'true':
		file3 = open(pop+'.sh','r')
		data3 = file3.read()
		print(data3)

	os.remove(pop+'.sh')

#CONCATENATE ALL WPM FILES BUT HAVE TO DEAL WITH HEADERS IN EACH.
		




		



