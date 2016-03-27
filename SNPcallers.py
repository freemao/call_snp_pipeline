#!/usr/lib/python
#-*- coding:utf-8 -*-
import datetime
from optparse import OptionParser
from call_snp_pipeline import FreebayesPipe
from subprocess import call

msg_usage = 'usage: %prog [-C] callers [-T] DANorRNA [-R] reference [-O] output'
descr = '''This script gather the different snp callers:freebyes, GATK, and
samtools.By pointing -c, You can choose any tool you like. You should choose
the data you have, such as fastq? sam? bam?. You also should point the type
of your data,DNA? or RNA?
'''
optparser = OptionParser(usage = msg_usage, description = descr)
optparser.add_option('-C', '--caller', dest = 'snpcaller',
                     help = 'point the caller you want to use. FB, SB, or GATK?')
optparser.add_option('-T', '--type', dest = 'datatype',
                     help = 'your data is DNA or RNA ?')
optparser.add_option('-R', '--ref', dest = 'refpath',
                     help = 'where is your reference seq ?')
optparser.add_option('-O', '--output', dest = 'output',
                     help = 'point your outcome file name.')
options, args = optparser.parse_args()

DNA_step = FreebayesPipe('.')
DNA_step.getrmpfilelist()
print 'files used to call DNA SNP %s'%DNA_step.namelist

RNA_step = FreebayesPipe('.')
RNA_step.getreadyfilelist()
print 'files used to call RNA SNP: %s'%RNA_step.namelist

def get_seqs(ref):
    '''in order to run SNPcallers simultaneously'''
    f0 = open(ref)
    contiglist = []
    for i in f0:
        if i.startswith('>'):
            seq = i.strip().split()[0][1:]
            contiglist.append(seq)
    print 'there are %s scaffolds in reference \
file:\n%s'%(len(contiglist),contiglist)
    return contiglist

def combine_vcf(yournamelist, final_name):
    nf = open(final_name, 'w')
    firfile = yournamelist[0]
    f1 = open(firfile)
    for i in f1:
        nf.write(i)
    f1.close()
    for j in yournamelist[1:]:
        f = open(j)
        for line in f:
            if line.startswith('#'):
                pass
            else:
                nf.write(line)
        f.close()
    nf.close()

def split_namelist(contiglist, n):
    '''if there are so many contig in reference file. I have to split
    seq namelist into smaller file because parallel many process may
    lead to cpu crazy~ so run parallel one by one'''
    newseqlist = []
    lens = len(contiglist)
    a = 0
    while a+n <= lens-1:
        nl = contiglist[a:a+n]
        newseqlist.append(nl)
        a += n
    else:
        nl = contiglist[a:]
        newseqlist.append(nl)
    return newseqlist

def runfreebayesfile(fileslist, ref, outputname):
    ''' no mnp and no complex!!!!'''
    seqlist = get_seqs(ref)
    f0 = open('bamfile.fb.list', 'w')
    for i in fileslist:
        f0.write(i + '\n')
    f0.close()
    splitnamelist = []
    f1 = open('run_freebayes.sh', 'w')
    for seq in seqlist:
        f1.write('freebayes -r %s -f %s -C 2 -F 0.01 -L bamfile.fb.list > %s\n'%(seq,\
ref, seq+'.'+outputname))
        splitnamelist.append(seq+'.'+outputname)
    f1.close()
    if len(seqlist) <= 15:
        call('parallel < run_freebayes.sh', shell = True)
        combine_vcf(splitnamelist, outputname)
    else:
        call('parallel -j 40 < run_freebayes.sh', shell = True)
        combine_vcf(splitnamelist, outputname)

def runsambcffile(fileslist, ref, outputname):
    seqlist = get_seqs(ref)
    f0 = open('bamfile.sb.list', 'w')
    for i in fileslist:
        f0.write(i + '\n')
    f0.close()
    splitnamelist = []
    f1 = open('run_samtools.sh', 'w')
    for seq in seqlist:
        f1.write('samtools mpileup -r %s -f %s -ugDV -b \
bamfile.sb.list | bcftools view -cNegv - > %s\n'%(seq, ref, seq+'.'+outputname))
        splitnamelist.append(seq+'.'+outputname)
    f1.close()
    if len(seqlist) <= 15:
        call('parallel < run_samtools.sh', shell = True)
        combine_vcf(splitnamelist, outputname)
    else:
        call('parallel -j 15 < run_samtools.sh', shell = True)
        combine_vcf(splitnamelist, outputname)

def runGATKfile(fileslist,ref, outputname):
    '''I dont now wether and how GATK can call a population's SNP.
    this version is just for one sample'''
    seqlist = get_seqs(ref)
    filenames = []
    for i in fileslist:
        filenames.append(i)
    addI = []
    for i in filenames:
        addI.append(' -I ' + i)
    filearg = ''.join(addI)
    f1 = open('run_gatk.sh', 'w')
    splitnamelist = []
    for seq in seqlist:
        f1.write('java -jar /share/Public/cmiao/GATK_tools/\
GenomeAnalysisTK.jar -nct 50 -allowPotentiallyMisencodedQuals \
-T HaplotypeCaller -L %s -R %s %s -o %s\n'%(seq, ref, filearg, seq+'.'+outputname))
        splitnamelist.append(seq+'.'+outputname)
    f1.close()
    call('chmod 777 run_gatk.sh', shell = True)
    call('./run_gatk.sh', shell = True)
    combine_vcf(splitnamelist, outputname)

if __name__ == '__main__':
    import sys
    C = options.snpcaller
    T = options.datatype
    O = options.output
    R = options.refpath

    if C == 'FB' and T == 'DNA':
        runfreebayesfile(DNA_step.namelist,R,O)
    elif C == 'GATK' and T == 'DNA':
        runGATKfile(DNA_step.namelist,R, O)
    elif C == 'SB' and T == 'DNA':
        runsambcffile(DNA_step.namelist,R,O)
    elif C == 'FB' and T == 'RNA':
        runfreebayesfile(RNA_step.namelist,R,O)
    elif C == 'GATK' and T == 'RNA':
        runGATKfile(RNA_step.namelist,R, O)
    elif C == 'SB' and T == 'RNA':
        runsambcffile(RNA_step.namelist,R, O)
    else:
        print 'Please choose the caller.'

