#!/usr/lib/python
#-*- coding:utf-8 -*-
from optparse import OptionParser
from call_snp_pipeline import FreebayesPipe
from subprocess import call

msg_usage = 'usage: %prog [-D] fqgz, fq, sam, bam... [-T] DNA,RNA [-R] refseq'
descr = '''This script mainly generates what SNP callers'
needs.You should choose the data you have, such as fastq? fq? sam? bam?. You
also should point the type of your data,DNA? or RNA?
Caveat: the fasq.gz file should follow the name format:
'F12-1_CTTGA_L006_R?.fastq.gz'
F: The material is from which part of your sample? one letter, for example,
L means leaf, F means flower, R means root.
12: Your sample name.Can not contain any symbols.
-: dash. 1: repeats,one digit. _: underline. CTTGA: barcode.
L006: lane?
'''
optparser = OptionParser(usage = msg_usage, description = descr)
optparser.add_option('-D', '--data', dest = 'dataformat',
                     help = 'point the data you have at present.\
                     fqgz? fq? sam? bam? rg? orderbam?')
optparser.add_option('-T', '--type', dest = 'datatype',
                     help = 'your data is DNA or RNA ?')
optparser.add_option('-R', '--reference', dest = 'reference',
                     help = 'point your reference sequence absolute path.')
optparser.add_option('-G', '--gtf', dest = 'gtf',
                     help = 'point your gtf file.')
options, args = optparser.parse_args()

def pre_DNA_fqgz(refseq):
    ref = FreebayesPipe()
    ref.pre_ref(refseq)
    step1 =  FreebayesPipe()
    step1.getgzfilelist()
    print 'fast.gz list: %s'%step1.namelist
    step1.pre_bwa(refseq)
    step1.runbwafile(refseq)

    step2 = FreebayesPipe()
    step2.getsamfilelist()
    print step2.namelist
    step2.runsam2bamfile()

    step3 = FreebayesPipe()
    step3.getbamfilelist()
    print step3.namelist
    step3.runsortfile()

    step4 = FreebayesPipe()
    step4.getsortfilelist()
    print step4.namelist
    step4.runrmdupfile()

    step5 = FreebayesPipe()
    step5.getrmpfilelist()
    print step5.namelist
    step5.runbaifile()
    print 'done~'

def pre_DNA_fq(refseq):
    ref = FreebayesPipe()
    ref.pre_ref(refseq)
    step1 =  FreebayesPipe()
    step1.getfqfilelist()
    print step1.namelist
    step1.pre_bwa(refseq)
    step1.runbwafile(refseq)

    step2 = FreebayesPipe()
    step2.getsamfilelist()
    print step2.namelist
    step2.runsam2bamfile()

    step3 = FreebayesPipe()
    step3.getbamfilelist()
    print step3.namelist
    step3.runsortfile()

    step4 = FreebayesPipe()
    step4.getsortfilelist()
    print step4.namelist
    step4.runrmdupfile()

    step5 = FreebayesPipe()
    step5.getrmpfilelist()
    print step5.namelist
    step5.runbaifile()
    print 'done~'

def pre_DNA_sam():
    step2 = FreebayesPipe()
    step2.getsamfilelist()
    print step2.namelist
    step2.runsam2bamfile()

    step3 = FreebayesPipe()
    step3.getbamfilelist()
    print step3.namelist
    step3.runsortfile()

    step4 = FreebayesPipe()
    step4.getsortfilelist()
    print step4.namelist
    step4.runrmdupfile()

    step5 = FreebayesPipe()
    step5.getrmpfilelist()
    print step5.namelist
    step5.runbaifile()
    print 'done~'

def pre_DNA_bam():
    step3 = FreebayesPipe()
    step3.getbamfilelist()
    print step3.namelist
    step3.runsortfile()

    step4 = FreebayesPipe()
    step4.getsortfilelist()
    print step4.namelist
    step4.runrmdupfile()

    step5 = FreebayesPipe()
    step5.getrmpfilelist()
    print step5.namelist
    step5.runbaifile()
    print 'done~'

def pre_DNA_sortedbam():
    step4 = FreebayesPipe()
    step4.getsortfilelist()
    print step4.namelist
    step4.runrmdupfile()

    step5 = FreebayesPipe()
    step5.getrmpfilelist()
    print step5.namelist
    step5.runbaifile()
    print 'done~'

def pre_DNA_rmpbam():
    step5 = FreebayesPipe()
    step5.getrmpfilelist()
    print step5.namelist
    step5.runbaifile()
    print 'done~'

def pre_RNA_fqgz(refseq, gtf):
    ref = FreebayesPipe('.')
    ref.pre_ref(refseq)
    step1 =  FreebayesPipe('.')
    step1.getgzfilelist()
    print 'your fastq.gz file list: %s'%step1.namelist
    print 'total %s files.'%len(step1.namelist)
    print 'next step: align using tophat2.'
    step1.pre_tophat2(refseq)
    step1.runtophat2file(refseq, gtf)

    step2 = FreebayesPipe('.')
    step2.getbamfilelist()
    print 'your bam file list: %s'%step2.namelist
    print 'total %s files.'%len(step2.namelist)
    print 'next step: remove duplicate reads in  your bam files.'
    step2.runrmdupfile()

    step4 = FreebayesPipe('.')
    step4.getrmpfilelist()
    print 'your removed bam files list: %s'%step4.namelist
    print 'total %s files.'%len(step4.namelist)
    print 'next step: add the read group info on each sample'
    step4.runaddrgfile()

    step5 = FreebayesPipe('.')
    step5.getrgfilelist()
    print 'your bam files with rg info list: %s'%step5.namelist
    print 'total %s files.'%len(step5.namelist)
    print 'next step: keep the order consistent between bam and ref.'
    step5.runorderfile(refseq)

    step6 = FreebayesPipe('.')
    step6.getorderfilelist()
    print 'your order bam files list: %s'%step6.namelist
    print 'total %s files.'%len(step6.namelist)
    step6.runbaifile()
    print 'next step: split N cigar'
    step6.runsplitNfile(refseq)

#    step7 = FreebayesPipe()
#    step7.getreadyfilelist()
#    print 'your splited N cigar files list: %s'%step7.namelist
#    print 'total %s files.'%len(step7.namelist)
#    print 'the next step: index your files used to call snp'
#    step7.runbaifile()
#    print 'done~'

def pre_RNA_bam(refseq):
    step2 = FreebayesPipe('.')
    step2.getbamfilelist()
    print 'your bam file list: %s'%step2.namelist
    print 'total %s files.'%len(step2.namelist)
    print 'next step: remove duplicate.'
    step2.runrmdupfile()

    step4 = FreebayesPipe('.')
    step4.getrmpfilelist()
    print 'your removed bam files list: %s'%step4.namelist
    print 'total %s files.'%len(step4.namelist)
    print 'next step: add the read group info on each sample'
    step4.runaddrgfile()

    step5 = FreebayesPipe('.')
    step5.getrgfilelist()
    print 'your bam files with rg info list: %s'%step5.namelist
    print 'total %s files.'%len(step5.namelist)
    print 'next step: keep the order consistent between bam and ref.'
    step5.runorderfile(refseq)

    step6 = FreebayesPipe('.')
    step6.getorderfilelist()
    print 'your order bam files list: %s'%step6.namelist
    print 'total %s files.'%len(step6.namelist)
    step6.runbaifile()
    print 'next step: split N cigar'
    step6.runsplitNfile(refseq)

    step7 = FreebayesPipe()
    step7.getreadyfilelist()
    print 'your splited N cigar files list: %s'%step7.namelist
    print 'total %s files.'%len(step7.namelist)
    print 'the next step: index your files used to call snp'
    step7.runbaifile()
    print 'done~'

def pre_RNA_sortedbam(refseq):
    step3 = FreebayesPipe('.')
    step3.getsortfilelist()
    print 'your sorted bam files list: %s'%step3.namelist
    print 'total %s files.'%len(step3.namelist)
    print 'next step: remove duplicate.'
    step3.runrmdupfile()

    step4 = FreebayesPipe('.')
    step4.getrmpfilelist()
    print 'your removed bam files list: %s'%step4.namelist
    print 'total %s files.'%len(step4.namelist)
    print 'next step: add the read group info on each sample'
    step4.runaddrgfile()

    step5 = FreebayesPipe('.')
    step5.getrgfilelist()
    print 'your bam files with rg info list: %s'%step5.namelist
    print 'total %s files.'%len(step5.namelist)
    print 'next step: keep the order consistent between bam and ref.'
    step5.runorderfile(refseq)

    step6 = FreebayesPipe('.')
    step6.getorderfilelist()
    print 'your order bam files list: %s'%step6.namelist
    print 'total %s files.'%len(step6.namelist)
    step6.runbaifile()
    print 'next step: split N cigar'
    step6.runsplitNfile(refseq)

    step7 = FreebayesPipe()
    step7.getreadyfilelist()
    print 'your splited N cigar files list: %s'%step7.namelist
    print 'total %s files.'%len(step7.namelist)
    print 'the next step: index your files used to call snp'
    step7.runbaifile()
    print 'done~'

def pre_RNA_rg_bam(refseq):
    '''if your file is bam files which have been added rg info'''
    step5 = FreebayesPipe('.')
    step5.getrgfilelist()
    print 'your bam files with rg info list: %s'%step5.namelist
    print 'total %s files.'%len(step5.namelist)
    print 'next step: keep the order consistent between bam and ref.'
    step5.runorderfile(refseq)

    step6 = FreebayesPipe('.')
    step6.getorderfilelist()
    print 'your order bam files list: %s'%step6.namelist
    print 'total %s files.'%len(step6.namelist)
    step6.runbaifile()
    print 'next step: split N cigar'
    step6.runsplitNfile(refseq)

    step7 = FreebayesPipe()
    step7.getreadyfilelist()
    print 'your splited N cigar files list: %s'%step7.namelist
    print 'total %s files.'%len(step7.namelist)
    print 'the next step: index your files used to call snp'
    step7.runbaifile()
    print 'done~'

def pre_RNA_order_bam(refseq):
    step6 = FreebayesPipe('.')
    step6.getorderfilelist()
    print 'your order bam files list: %s'%step6.namelist
    print 'total %s files.'%len(step6.namelist)
    step6.runbaifile()
    print 'next step: split N cigar'
    step6.runsplitNfile(refseq)

    step7 = FreebayesPipe()
    step7.getreadyfilelist()
    print 'your splited N cigar files list: %s'%step7.namelist
    print 'total %s files.'%len(step7.namelist)
    print 'the next step: index your files used to call snp'
    step7.runbaifile()
    print 'done~'

def pre_RNA_Nsplited_bam():
    step7 = FreebayesPipe()
    step7.getreadyfilelist()
    print 'your splited N cigar files list: %s'%step7.namelist
    print 'total %s files.'%len(step7.namelist)
    print 'the next step: index your files used to call snp'
    step7.runbaifile()
    print 'done~'


if __name__ == '__main__':
    T = options.datatype
    D = options.dataformat
    R = options.reference
    G = options.gtf
    if T == 'DNA' and D == 'fqgz':
        pre_DNA_fqgz(R)
    elif T == 'DNA' and D == 'fq':
        pre_DNA_fq(R)
    elif T == 'DNA' and D == 'sam':
        pre_DNA_sam()
    elif T == 'DNA' and D == 'bam':
        pre_DNA_bam()
    elif T == 'DNA' and D == 'sortedbam':
        pre_DNA_sortedbam()
    elif T == 'DNA' and D == 'rmpbam':
        pre_DNA_rmpbam()
    elif T == 'RNA' and D == 'fqgz':
        pre_RNA_fqgz(R, G)
    elif T == 'RNA' and D == 'bam':
        pre_RNA_bam(R)
    elif T == 'RNA' and D == 'sortedbam':
        pre_RNA_sortedbam(R)
    elif T == 'RNA' and D == 'rgbam':
        pre_RNA_rg_bam(R)
    elif T == 'RNA' and D == 'orderbam':
        pre_RNA_order_bam(R)
    elif T == 'RNA' and D == 'Nsplitedbam':
        pre_RNA_Nsplited_bam()
    else:
        print 'Please point your data format and data type.'
