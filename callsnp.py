#!/usr/lib/python
#-*- coding:utf-8 -*-
from optparse import OptionParser
from call_snp_pipeline import FreebayesPipe
from subprocess import call

msg_usage = 'usage: %prog [-D] [-C] [-T] [-R]'
descr = '''This script gather the different snp callers:freebyes, GATK, and
samtools.By pointing -c, You can choose any tool you like. You should choose
the data you have, such as fastq? sam? bam?. You also should point the type
of your data,DNA? or RNA?
'''
optparser = OptionParser(usage = msg_usage, description = descr)
optparser.add_option('-D', '--data', dest = 'dataformat',
                     help = 'point the data you have at present.')
optparser.add_option('-C', '--caller', dest = 'snpcaller',
                     help = 'point the caller you want to use.')
optparser.add_option('-T', '--type', dest = 'datatype',
                     help = 'your data is DNA or RNA ?')
optparser.add_option('-R', '--reference', dest = 'referencename',
                     help = 'point your reference sequence name.')
options, args = optparser.parse_args()

step = FreebayesPipe('.')
step.getrmpfilelist()
print step.namelist

def runfreebayesfile():
    f0 = open('bamfile.fb.list', 'w')
    vcfname = set()
    for i in step.namelist:
        j = i.split('-')[0][1:]
        vcfname.add(j)
        f0.write(i + '\n')
    f0.close()
    f1 = open('run_freebayes.sh', 'w')
    if len(vcfname) == 1:
        f1.write('freebayes -f Osativa_204.fa -F 0.1 -L bamfile.fb.list > ' \
+ j + '.fb.vcf')
    else:
        print 'the vcf file name is not unique! check your files please.'
    f1.close()
    call('chmod 777 run_freebayes.sh', shell = True)
    call('./run_freebayes.sh', shell = True)

def runsambcffile():
    f0 = open('bamfile.sb.list', 'w')
    vcfname = set()
    for i in step.namelist:
        j = i.split('-')[0][1:]
        vcfname.add(j)
        f0.write(i + '\n')
    f0.close()
    f1 = open('run_samtools1.sh', 'w')
    f2 = open('run_bcftools2.sh', 'w')
    if len(vcfname) == 1:
        f1.write('samtools mpileup -f Osativa_204.fa -P ILLUMINA -EgD -b \
bamfile.sb.list > ' + j + '.sb.bcf')
        f2.write('bcftools view -cNegv ' + j + '.sb.bcf' + ' > ' + j + \
'.sb.vcf')
    else:
        print 'the vcf file name is not unique! check your files please.'
    f1.close()
    call('chmod 777 run_samtools1.sh', shell = True)
    call('chmod 777 run_bcftools2.sh', shell = True)
    call('./run_samtools1.sh', shell = True)
    call('./run_bcftools2.sh', shell = True)

def runGATKfile():
    filenames = []
    vcfname = set()
    for i in step.namelist:
        j = i.split('-')[0][1:]
        vcfname.add(j)
        filenames.append(i)
    addI = []
    for i in filenames:
        addI.append(' -I ' + i)
    filearg = ''.join(addI)
    f1 = open('run_gatk.sh', 'w')
    if len(vcfname) == 1:
        f1.write('java -jar /share/Public/cmiao/GATK_tools/\
GenomeAnalysisTK.jar -nct 30 -T HaplotypeCaller -R Osativa_204.fa' + \
filearg + ' -o ' + j +'.gatk.vcf')
    else:
        print 'the vcf file name is not unique! check your files please.'
    f1.close()
    call('chmod 777 run_gatk.sh', shell = True)
    call('./run_gatk.sh', shell = True)

if __name__ == '__main__':
    import sys
    C = options.snpcaller
    if C == 'FB':
        runfreebayesfile()
    elif C == 'GATK':
        runGATKfile()
    elif C == 'SB':
        runsambcffile()
    else:
        print 'Please choose the caller.'



