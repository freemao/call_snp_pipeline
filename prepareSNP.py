#!/usr/lib/python
#-*- coding:utf-8 -*-
from optparse import OptionParser
from call_snp_pipeline import FreebayesPipe
from subprocess import call

msg_usage = 'usage: %prog [-D] [-T] [-R]'
descr = '''This script mainly create the datas which meet the snp callers'
needs.You should choose the data you have, such as fastq? sam? bam?. You
also should point the type of your data,DNA? or RNA?
'''
optparser = OptionParser(usage = msg_usage, description = descr)
optparser.add_option('-D', '--data', dest = 'dataformat',
                     help = 'point the data you have at present.')
optparser.add_option('-T', '--type', dest = 'datatype',
                     help = 'your data is DNA or RNA ?')
#optparser.add_option('-R', '--reference', dest = 'referencename',
#                     help = 'point your reference sequence name.')
options, args = optparser.parse_args()

def pre_DNA_fqgz():
    ref = FreebayesPipe('.')
    ref.pre_ref()
    step1 =  FreebayesPipe('.')
    step1.getgzfilelist()
    print step1.namelist
    step1.pre_bwa()
    step1.run_bwa()

    step2 = FreebayesPipe('.')
    step2.getsamfilelist()
    print step2.namelist
    step2.runsam2bamfile()

    step3 = FreebayesPipe('.')
    step3.getbamfilelist()
    print step3.namelist
    step3.runsortfile()

    step4 = FreebayesPipe('.')
    step4.getsortfilelist()
    print step4.namelist
    step4.rmdupfile()
    step4.runbaifile()

def pre_DNA_bam():
    step3 = FreebayesPipe('.')
    step3.getbamfilelist()
    print step3.namelist
    step3.runsortfile()

    step4 = FreebayesPipe('.')
    step4.getsortfilelist()
    print step4.namelist
    step4.rmdupfile()
    step4.runbaifile()

def pre_RNA_fqgz():
    ref = FreebayesPipe('.')
    ref.pre_ref()
    step1 =  FreebayesPipe('.')
    step1.getgzfilelist()
    print step1.namelist
    step1.pre_tophat2()
    step1.runtophat2file()

    step2 = FreebayesPipe('.')
    step2.getbamfilelist()
    print step2.namelist
    step2.runsortfile()

    step3 = FreebayesPipe('.')
    step3.getsortfilelist()
    print step3.namelist
    step3.runrmdupfile()

    step4 = FreebayesPipe('.')
    step4.getrmpfilelist()
    print step4.namelist
    step4.runaddrgfile()

    step5 = FreebayesPipe('.')
    step5.getrgfilelist()
    print step5.namelist
    step5.runorderfile()

    step6 = FreebayesPipe('.')
    step6.getorderedfilelist()
    print step6.namelist
    step6.runsplitNtrimfile()

def pre_RNA_bam():
    step2 = FreebayesPipe('.')
    step2.getbamfilelist()
    print step2.namelist
    step2.runsortfile()

    step3 = FreebayesPipe('.')
    step3.getsortfilelist()
    print step3.namelist
    step3.runrmdupfile()

    step4 = FreebayesPipe('.')
    step4.getrmpfilelist()
    print step4.namelist
    step4.runaddrgfile()

    step5 = FreebayesPipe('.')
    step5.getrgfilelist()
    print step5.namelist
    step5.runorderfile()

    step6 = FreebayesPipe('.')
    step6.getorderedfilelist()
    print step6.namelist
    step6.runbaifile()
    step6.runsplitNtrimfile()


def pre_RNA_ordered_bam():
    step6 = FreebayesPipe('.')
    step6.getorderedfilelist()
    print step6.namelist
    step6.runbaifile()
    step6.runsplitNtrimfile()

if __name__ == '__main__':
    T = options.datatype
    D = options.dataformat
    if T == 'DNA' and D == 'fqgz':
        pre_DNA_fqgz()
    elif T == 'DNA' and D == 'bam':
        pre_DNA_bam()
    elif T == 'RNA' and D == 'fqgz':
        pre_RNA_fqgz()
    elif T == 'RNA' and D == 'bam':
        pre_RNA_bam()
    elif T == 'RNA' and D == 'orderedbam':
        pre_RNA_ordered_bam()
    else:
        print 'Please point your data format and data type.'
