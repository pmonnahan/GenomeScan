"""Calculate within-population metrics for diversity and selection.  Input is a tab delimited file with no header, containing Scaff, pos, ac, an, dp, and genotypes coded as 0-4. 
Input file custom: Filename should end with _XXX_raw.table.recode.txt, where XXX is three-letter population abbreviation"""

import os, sys, subprocess, argparse,pandas,math,numpy
from scipy import stats

# create variables that can be entered in the command line
parser = argparse.ArgumentParser()
parser.add_argument('-i', type=str, metavar='input_file', required=True, help='input file created with recode012.py')
parser.add_argument('-o', type=str, metavar='output_directory',required=True,help='Output Directory')
parser.add_argument('-m', type=float, metavar='Missingness',required=False,default='0.0',help='Proportion of missing individuals allowed in each group')
parser.add_argument('-ws', type=float, metavar='window_size',required=False,default='10000.0',help='Size of windows in bp')
parser.add_argument('-ms', type=int, metavar='minimum_snps',required=False,default='2',help='minimum number of snps in a window')

args = parser.parse_args()

if args.i.endswith("_raw.table.recode.txt"):
	outfile=args.i.split("/")[-1].strip("_raw.table.recode.txt")+"_Within-Pop-Metrics.txt"
	print outfile

if args.o.endswith("/")!=True:
	args.o += "/"

out1=open(args.o+outfile,'w')
out1.write("scaff\tstart\tend\twin_size\tnum_snps\tavg_freq\tavg_Ehet\tThetaW\tPi\tThetaH\tThetaL\tD\tH\tE\n")

snp_count = 0
start = 0.0
end = args.ws
winexclcount=0

with open(args.i, 'r') as infile:
	for i,line in enumerate(infile):

		line=line.strip("\n")
		line=line.strip("\t")
		line=line.split("\t")

		scaff,pos,ac,an,dp = line[:5]

		gt = line[5:]

		if i%100000==0:
			print i	
		if i == 0:
			totind = len(line[5:])
			ploidy = int(an)/totind
			sampind =int(math.ceil(totind * (1.0-args.m))) 
			AN = sampind*4
			n = float(AN)
			p = []
			Ehet =[]
			afs = [0 for cat in range(0,AN+1)]
			AFS = [0 for cat in range(0,AN+1)]
			aw = 0.0
			bw = 0.0 #This is b sub n+1 in Zeng. a just goes to n but b goes to n+1
			a2 = 0.0
			for j in range(1,AN):
				aw += 1.0/float(j) #a1 in Tajima 1989 and an in Zeng 2006
				a2 += 1.0/float(j**2) #This is bn according to Zeng 2006
				bw += 1.0/float(j**2)
			bw += 1.0/float(n**2) #This is the n+1 part
			b1 = (n + 1)/(3*(n - 1)) #From Tajima 1989
			b2 = (2*n**2 + n + 3)/(9*n*(n - 1))
			c1 = b1 - (1/aw)
			c2 = b2 - (n + 2)/(aw*n) + a2/aw**2
			e1 = c1/aw
			e2 = c2/(aw**2 + a2)

		if int(pos) > start and int(pos) <=end and an>=AN:
			snp_count+=1
			sgt=numpy.random.choice(gt,size=sampind,replace=False)
			sac=sum([int(x) for x in sgt])
			p1 = float(sac)/float(AN)
			p.append(p1)
			Ehet.append(p1*(1-p1))
			afs[sac]+=1
			AFS[sac]+=1

		elif int(pos) > end:
			Pi = 0.0
			h = 0.0
			L = 0.0
			n = float(AN)
			if snp_count >= args.ms:
				S = float(sum(afs[1:-1]))
				W = S/aw
				W2 = S*(S-1)/((aw**2)+bw)
				for j in range(1,AN):

					Pi += afs[j]*j*(AN - j)
					h += afs[j]*(j**2)
					L += afs[j]*j

				Pi = 2*Pi/(n*(n-1))
				h = 2*h/(n*(n-1))
				L = L/(n - 1)
				
				#When calculating variance of D, H, and E should we use W value for window or for entire genome.  Probably for window??

				varPi_W = e1*S + e2*S*(S-1)
				varPi_L = (((n - 2)/(6*n - 6))*W) + (((18*n**2*(3*n + 2)*bw) - (88*n**3 + 9*n**2 - 13*n + 6))/(9*n*(n - 1)**2)*W2)
				varL_W = (((n/(2*n-2)) - 1/aw)*W) + ((a2/aw**2 + (2*(n/(n - 1))**2)*a2 - 2*(n*a2 - n + 1)/((n - 1)*aw) - (3*n + 1)/(n - 1))*W2)
				print "1",varPi_W,varPi_L, varL_W,W,W2
				print "2",((n - 2)/(6*n - 6)),((18*(n**2)*(3*n + 2)*bw) - (88*n**3 + 9*n**2 - 13*n + 6))/(9*n*(n - 1)**2)
				print "3",((n/(2*n-2)) - 1/aw),(a2/aw**2 + (2*(n/(n - 1))**2)*a2 - 2*(n*a2 - n + 1)/(n - 1)*aw - (3*n + 1)/(n - 1))
				print "4",aw,a2,b1,b2,c1,c2,e1,S,e2,S-1

				D = (Pi - W)/(varPi_W**0.5)
				H = (Pi - L)/(varPi_L**0.5) #Normalized H according to Zeng 2006
				E = (L - W)/(varL_W**0.5)

				out1.write(scaff+'\t'+
				str(start)+'\t'+
				str(end)+'\t'+
				str(args.ws)+'\t'+
				str(snp_count)+'\t'+
				str(numpy.mean(p1))+'\t'+
				str(numpy.mean(Ehet))+'\t'+
				str(W)+'\t'+
				str(Pi)+'\t'+
				str(h)+'\t'+
				str(L)+'\t'+
				str(D)+'\t'+
				str(H)+'\t'+
				str(E)+'\n')

				

			else:
				winexclcount +=1

			snp_count = 0
			p = []
			Ehet =[]
			afs = [0 for cat in range(0,AN+1)]
			AFS = [0 for cat in range(0,AN+1)]

			while float(pos) > end:
				end+=args.ws/2

			start=end-args.ws

out1.close()




			


		