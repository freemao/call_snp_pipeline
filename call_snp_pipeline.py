#!/usr/lib/python
#-*- coding:utf-8 -*-
import os
from subprocess import call

class FreebayesPipe:
    '''This FreebayesPipe class contains all the  steps to prepare to call
SNP whether for DNA data or RNA data. For DNA data, the synopsis is 1,fq.gz.
2,sam. 3,bam. 4,sorted.bam. 5,sorted.rmp.bam.6,sorted.rmp.bam.bai. While for
the RNA data, the synopsis is 1,fq.gz. 2,bam. 3,sorted.bam. 4,sorted.rmp.bam.
5,sorted.rmp.rg.bam. 6,ordered.bam. 7,Nsplited.ready.bam.'''
    def __init__(self, dirname='.'):
        self.dirname = dirname
        self.allnamelist = os.listdir(dirname)
        self.namelist = []
        self.arg1 = []
        self.arg2 = []
        self.arg3 = []
        self.arg4 = []

    def pre_ref(self, refseq):
        '''designate refseq path'''
        abspathrefseq = os.path.abspath(refseq)
        if os.path.isfile(abspathrefseq):
            print 'the reference: %s'%abspathrefseq
        else:
            return 'rhe reference %s does not exists!!'%abspathrefseq
        dirseq = os.path.dirname(abspathrefseq)
        all_suffixes = [i.split('.')[-1] for i in os.listdir(dirseq)]
        if 'fai' not in all_suffixes:
            print 'reference index file does not exist, building...'
            cmd1 = 'samtools faidx %s'%abspathrefseq
            call(cmd1, shell=True)
            print 'reference index is built.'
        if 'dict' not in all_suffixes:
            dictname = '.'.join(abspathrefseq.split('.')[0:-1])+'.dict'
            cmd2 = 'java -jar /share/Public/cmiao/picard-tools-1.112/\
CreateSequenceDictionary.jar R=%s O=%s'%(abspathrefseq, dictname)
            call(cmd2, shell=True)
            print 'reference dict is built.'

    def getgzfilelist(self):
        for fn in self.allnamelist:
            temp = os.path.join(self.dirname, fn)
            if (os.path.isfile(temp) and temp.split('.')[-1] == 'gz'):
                self.namelist.append(fn)
                self.namelist.sort()

    def rungzfile(self):
        f0 = open('run_gz.txt', 'w')
        for fn in self.namelist:
            f0.write('gzip -d ' + fn + '\n')
        f0.close()
        call('parallel < run_gz.txt', shell = True)

    def getfqfilelist(self):
        for fn in self.allnamelist:
            temp = os.path.join(self.dirname, fn)
            if (os.path.isfile(temp)
                and temp.split('.')[-1] in ['fastq', 'fq']):
                self.namelist.append(fn)
                self.namelist.sort()

    def pre_bwa(self, refseq):
        abspathrefseq = os.path.abspath(refseq)
        dirseq = os.path.dirname(abspathrefseq)
        index_suffixes = ('amb', 'ann', 'bwt', 'pac', 'sa')
        all_suffixes = [i.split('.')[-1] for i in os.listdir(dirseq)]
        if set(index_suffixes).issubset(all_suffixes):
            print 'The bwa index has been built.'
        else:
            print "The reference sequence don't build bwa index yet. \
building..."
            cmd = 'bwa index %s'%abspathrefseq
            call(cmd, shell=True)

    def runbwafile(self, refseq):
        print 'running bwa_mem...'
        L = range(0, len(self.namelist), 2)
        for i in L:
            SM = self.namelist[i].split('-')[0][1:]
            R = r"'@RG\tID:%s\tSM:%s'"%(SM, SM)
            sam_prefix = '_'.join(self.namelist[i].split('_')[:-1])
            self.arg1.append(self.namelist[i])
            self.arg2.append(self.namelist[i+1])
            self.arg3.append(sam_prefix+'.sam')
            self.arg4.append(R)
        f0 = open('run_bwa.txt', 'w')
        for x, y, z, w in zip(self.arg1, self.arg2, self.arg3, self.arg4):
            cmd1 = 'bwa mem -t 1 -R %s %s %s %s > %s \n'%(w, refseq, x, y, z)
            f0.write(cmd1)
        f0.close()
        call('parallel < run_bwa.txt', shell = True)

    def runbwa_aln_file(self, refseq):
        print 'run bwa_aln...'
        f1 = open('run_sampe.txt', 'w')
        L = range(0, len(self.namelist), 2)
        for i in L:
            f0 = open('run_sai.txt','w')
            sai1 = self.namelist[i].split('.')[0]+'.sai'
            sai2 = self.namelist[i+1].split('.')[0]+'.sai'
            cmd1 = 'bwa aln -n 0.04 -t 20 %s %s > %s\n'%(refseq,\
self.namelist[i],sai1)
            f0.write(cmd1)
            cmd2 = 'bwa aln -n 0.04 -t 20 %s %s > %s\n'%(refseq,\
self.namelist[i+1],sai2)
            f0.write(cmd2)
            f0.close()
            call('parallel < run_sai.txt', shell = True)
            sam = '_'.join(self.namelist[i].split('_')[:-1])+'.sam'
            SM = self.namelist[i].split('-')[0][1:]
            R = r"'@RG\tID:%s\tSM:%s'"%(SM, SM)
            cmd3 = 'bwa sampe -r %s -f %s %s %s %s %s %s\n'%(R, sam, refseq,\
sai1, sai2, self.namelist[i], self.namelist[i+1])
            f1.write(cmd3)
        f1.close()
        call('parallel -j 20 < run_sampe.txt', shell = True)

    def pre_tophat2(self, refseq):
        abspathrefseq = os.path.abspath(refseq)
        dirseq = os.path.dirname(abspathrefseq)
        all_suffixes = [i.split('.')[-1] for i in os.listdir(dirseq)]
        index_suffixes = ('bt2',)
        if set(index_suffixes).issubset(all_suffixes):
            print 'The tophat2 reference index has been built.'
        else:
            print "The reference sequence don't build bowtie2 index yet. \
building..."
            prefix_seq = '.'.join(abspathrefseq.split('.')[0:-1])
            print 'the prefix of refseq name: %s'%(prefix_seq)
            cmd = 'bowtie2-build -f %s %s'%(abspathrefseq, prefix_seq)
            print cmd
            call(cmd, shell=True)

    def runtophat2file(self, refseq, gtf):
        '''if you have gtf or gtf, you'd better provide it.'''
        abspathrefseq = os.path.abspath(refseq)
        abspathgtf = os.path.abspath(gtf)
        faIdxName = '.'.join(abspathrefseq.split('.')[0:-1])
        L = range(0, len(self.namelist), 2)
        for i in L:
            self.arg1.append(self.namelist[i])
            self.arg2.append(self.namelist[i+1])
        f0 = open('run_tophat2.sh', 'w')
        for x, y in zip(self.arg1, self.arg2):
            cmd = 'tophat2 -g 1 -o . -p 40 --no-discordant -G %s \
-r 100 --mate-std-dev 50  %s %s %s \n'%(abspathgtf,faIdxName, x, y)
            f0.write(cmd)
            bamname = '_'.join(x.split('_')[:-1]) + '.bam'
            unmapedBamName = '_'.join(x.split('_')[:-1]) + '.unmap.bam'
            f0.write('mv ./accepted_hits.bam ' + bamname + '\n')
            f0.write('mv ./unmapped.bam ' + unmapedBamName + '\n')
        f0.close()
        call('chmod 777 run_tophat2.sh', shell = True)
        call('./run_tophat2.sh', shell = True)

    def getsamfilelist(self):
        for fn in self.allnamelist:
            temp = os.path.join(self.dirname,fn)
            if (os.path.isfile(temp) and temp.split('.')[-1] == 'sam') :
                self.namelist.append(fn)
                self.namelist.sort()

    def runsam2bamfile(self):
        print 'running sam to bam...'
        for i in self.namelist:
            self.arg1.append(i)
            self.arg2.append(i.split('.')[0] + '.bam')
        f0 = open('run_sam2bam.txt', 'w')
        for x, y in zip(self.arg1, self.arg2):
            f0.write('samtools view -bS ' + x + ' > ' + y + '\n')
        f0.close()
        call('parallel < run_sam2bam.txt', shell = True)

    def getbamfilelist(self):
        for fn in self.allnamelist:
            temp = os.path.join(self.dirname, fn)
            if (os.path.isfile(temp) and fn.split('.')[-1] == 'bam'
                and not set(['unmap','sorted', 'rmp', 'rg', 'reorder','splitN']) & set(fn.split('.')) ):
                self.namelist.append(fn)
                self.namelist.sort()

    def runsortfile(self):
        print 'running sort bam...'
        for i in self.namelist:
            self.arg1.append(i)
            self.arg2.append(i.split('.')[0] + '.sorted')
        f0 = open('run_sort.txt', 'w')
        for x, y in zip(self.arg1, self.arg2):
            f0.write('samtools sort ' + x + ' ' + y + '\n')
        f0.close()
        call('parallel < run_sort.txt', shell =True)

    def getsortfilelist(self):
        for fn in self.allnamelist:
            temp = os.path.join(self.dirname, fn)
            if (os.path.isfile(temp)
                and fn.split('.')[-2:] == ['sorted', 'bam']
                and not set(['rmp', 'rg']) & set(fn.split('.')) ):
                self.namelist.append(fn)
                self.namelist.sort()

    def runrmdupfile(self):
        print 'running remove duplicate reads...'
        for i in self.namelist:
            self.arg1.append(i)
            self.arg2.append(i.split('.')[0] + '.rmp.bam')
        f0 = open('run_rmp.txt', 'w')
        for x, y in zip(self.arg1, self.arg2):
            f0.write('samtools_0.1.18 rmdup ' + x + ' ' + y + '\n')
        f0.close()
        call('parallel < run_rmp.txt', shell = True)

    def getrmpfilelist(self):
        for fn in self.allnamelist:
            temp = os.path.join(self.dirname, fn)
            if (os.path.isfile(temp)
                and '.'.join(fn.split('.')[-2:]) == 'rmp.bam'):
                self.namelist.append(fn)
                self.namelist.sort()

    def runaddrgfile(self):
        f0 = open('run_addrg.txt', 'w')
        for i in self.namelist:
            self.arg1.append(i)
            rg_name = i.split('_')[0]
            self.arg2.append(rg_name)
            self.arg3.append(i.split('.')[0] + '.rmp.rg.bam')
        for x, z, y in zip(self.arg1, self.arg2, self.arg3):
            f0.write('bamaddrg -b %s -s %s -r %s > %s\n'%(x, z, z, y))
        f0.close()
        call('parallel < run_addrg.txt', shell = True)

    def getrgfilelist(self):
        for fn in self.allnamelist:
            seg = fn.split('.')
            temp = os.path.join(self.dirname, fn)
            if (os.path.isfile(temp)
                and '.'.join(fn.split('.')[-3:]) == 'rmp.rg.bam'):
                self.namelist.append(fn)
                self.namelist.sort()

    def getrecalfilelist(self):
        for fn in self.allnamelist:
            seg = fn.split('.')
            temp = os.path.join(self.dirname, fn)
            if (os.path.isfile(temp)
                and seg[1:] == ['sorted','rmp','rg','recal', 'bam']):
                self.namelist.append(fn)
                self.namelist.sort()

    def runbaifile(self):
        print 'buiding index of bam file...'
        for i in self.namelist:
            self.arg1.append(i)
        f0 = open('run_bai.txt', 'w')
        for x in self.arg1:
            cmd = 'samtools index %s'%x
            f0.write(cmd + '\n')
        f0.close()
        call('parallel < run_bai.txt', shell = True)

    def runorderfile(self, refseq):
        '''let the bam and refseq have same seqname order'''
        for i in self.namelist:
            self.arg1.append(i)
            self.arg2.append(i.split('.')[0] + '.rmp.rg.reorder.bam')
        f0 = open('run_order.txt', 'w')
        for x, y in zip(self.arg1, self.arg2):
            cmd = 'java -jar /share/Public/cmiao/picard-tools-1.112/\
ReorderSam.jar I=%s O=%s REFERENCE=%s \n'%(x, y, refseq)
            f0.write(cmd)
        f0.close()
        call('parallel < run_order.txt', shell = True)

    def getorderfilelist(self):
        for fn in self.allnamelist:
            seg = fn.split('.')
            temp = os.path.join(self.dirname, fn)
            if (os.path.isfile(temp)
                and '.'.join(fn.split('.')[-4:]) == 'rmp.rg.reorder.bam'):
                self.namelist.append(fn)
                self.namelist.sort()

    def runsplitNfile(self, refseq):
        for i in self.namelist:
            self.arg1.append(i)
            self.arg2.append(i.split('.')[0] + '.rmp.rg.reorder.splitN.bam')
        f0 = open('run_splitN.txt', 'w')
        for x, y in zip(self.arg1, self.arg2):
            cmd = 'java -jar /share/bioinfo/miaochenyong/ASE-PROJECT/GenomeAnalysisTK.jar \
-T SplitNCigarReads -I %s -U ALLOW_N_CIGAR_READS -o %s -R %s \n'%(x, y, refseq)
            f0.write(cmd)
        f0.close()
        call('parallel < run_splitN.txt', shell = True)

    def getreadyfilelist(self):
        for fn in self.allnamelist:
            seg = fn.split('.')
            temp = os.path.join(self.dirname, fn)
            if (os.path.isfile(temp)
                and '.'.join(fn.split('.')[-5:]) == 'rmp.rg.reorder.splitN.bam'):
                self.namelist.append(fn)
                self.namelist.sort()

if __name__ == '__main__':
    print 'Feel free to use this class!'
    print FreebayesPipe.__doc__
