import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-i', type=str, metavar='input_file', required=True, help='path to csv file containing information on what samples to include')
parser.add_argument('-o', type=str, metavar='output_dir', required=True, help='Output Directory')
parser.add_argument('-length', type=float, metavar='genome_length', required=True, help='Total size of the genome to be simulated')
parser.add_argument('-ploidy', type=int, metavar='ploidy', required=True, help='ploidy of populations to be simulated')

args = parser.parse_args()
if os.path.exists(args.o + 'MSrecodePloidy' + str(args.ploidy) + '/') is False:
    os.mkdir(args.o + 'MSrecodePloidy' + str(args.ploidy) + '/')
outdir = args.o + 'MSrecodePloidy' + str(args.ploidy) + '/'
with open(args.i) as inp:
    pop = 0
    for line in inp:
        line = line.strip('\n')
        line = line.strip(' ')
        line = line.split(' ')

        if line[0] == '':
            pass
        elif line[0] == '//':
            if pop != 0:  # Write final individual's genotype info from previous population
                outline = ''
                for j in ind:
                    outline += str(j) + '\t'
                outline.strip('\t')
                exec("out%d.write('%s')" % (pop, outline))
                exec('out%d.write("""\n""")' % (pop))
                exec("Out%d = open('%sMS%d.table.recode.txt','w')" % (pop, outdir, pop))
                with open(outdir + 'MS' + str(pop) + '.table.recode.trans.txt','r') as f:
                    lis = [x.strip('\n').split('\t') for x in f]
                    for x in zip(*lis):
                        for y in x:
                            exec("Out%d.write('%s\t')" % (pop, y))
                        exec('Out%d.write("""\n""")' % (pop))
            
                os.remove(outdir + 'MS' + str(pop) + '.table.recode.trans.txt')
            hap = 0
            pop += 1
            print(pop)
            exec("out%d = open('%sMS%d.table.recode.trans.txt','w')" % (pop, outdir, pop))
        elif line[0] == 'segsites:':
            numsites = int(line[1])
            outline = ''
            for j in range(0, numsites):
                outline += 'MS' + str(pop)+ '\t'
            outline.strip('\t')
            exec("out%d.write('%s')" % (pop, outline))
            exec('out%d.write("""\n""")' % (pop))
            outline = ''
            for j in range(0, numsites):
                outline += str(args.ploidy) + '\t'
            outline.strip('\t')
            exec("out%d.write('%s')" % (pop, outline))
            exec('out%d.write("""\n""")' % (pop))
            outline = ''
            for j in range(0, numsites):
                outline += 'scaffold_1\t'
            outline.strip('\t')
            exec("out%d.write('%s')" % (pop, outline))
            exec('out%d.write("""\n""")' % (pop))
        elif line[0] == 'positions:':
            positions = line[1:]
            positions = [int(round(args.length * float(x))) for x in positions]
            outline = ''
            for j in positions:
                outline += str(j) + '\t'
            outline.strip('\t')
            exec("out%d.write('%s')" % (pop, outline))
            exec('out%d.write("""\n""")' % (pop))
            outline = ''
            for j in range(0, numsites):
                outline += '99\t'
            outline.strip('\t')
            exec("out%d.write('%s')" % (pop, outline))
            exec('out%d.write("""\n""")' % (pop))
            utline = ''
            for j in range(0, numsites):
                outline += '99\t'
            outline.strip('\t')
            exec("out%d.write('%s')" % (pop, outline))
            exec('out%d.write("""\n""")' % (pop))
            utline = ''
            for j in range(0, numsites):
                outline += '99\t'
            outline.strip('\t')
            exec("out%d.write('%s')" % (pop, outline))
            exec('out%d.write("""\n""")' % (pop))
            ind = [0 for kk in range(0, numsites)]

        elif line[0][:2] == '01' or line[0][:2] == '10' or line[0][:2] == '11' or line[0][:2] == '00':
            hap += 1
            if hap <= args.ploidy:
                for a, allele in enumerate(line[0]):
                    ind[a] += int(allele)

            else:
                outline = ''
                for j in ind:
                    outline += str(j) + '\t'
                outline.strip('\t')
                exec("out%d.write('%s')" % (pop, outline))
                exec('out%d.write("""\n""")' % (pop))
                hap = 1
                ind = [0 for kk in range(0, numsites)]
                for a, allele in enumerate(line[0]):
                    ind[a] += int(allele)