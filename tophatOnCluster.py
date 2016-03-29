import os.path
def getjobs(refseq,gtf):
    header = '#!/bin/sh\n\n#PBS -N ASE_project\n#PBS -l nodes=1:ppn=10\n#PBS -o \
/share/workplace/home/cmiao/lastjob-pbs_out.$PBS_JOBID\n#PBS -e \
/share/workplace/home/cmiao/lastjob-pbs_err.$PBS_JOBID\n#PBS -l \
walltime=10000:00:00\n#PBS -q batch\n'
    namelist = []
    allfiles = os.listdir('.')
    for fn in allfiles:
        temp = os.path.join('.', fn)
        if (os.path.isfile(temp) and temp.split('.')[-1] == 'gz'):
            namelist.append(fn)
    namelist.sort()
    abspathrefseq = os.path.abspath(refseq)
    abspathgtf = os.path.abspath(gtf)
    faIdxName = '.'.join(abspathrefseq.split('.')[0:-1])
    L = range(0, len(namelist), 2)
    arg1 = []
    arg2 = []
    for i in L:
        arg1.append(namelist[i])
        arg2.append(namelist[i+1])
    for x, y in zip(arg1, arg2):
        oDir = x.split('_')[0]
        abspathoDir = os.path.abspath(oDir)
        abspathx = os.path.abspath(x)
        abspathy = os.path.abspath(y)
        cmd1 = 'tophat -g 1 -o %s -p 10 --no-discordant -G %s \
-r 100 --mate-std-dev 50  %s %s %s \n'%(abspathoDir, abspathgtf,faIdxName,\
abspathx, abspathy)
        bamname = '_'.join(x.split('_')[:-1]) + '.bam'
        unmapedBamName = '_'.join(x.split('_')[:-1]) + '.unmap.bam'
        cmd2 = 'ln -s %s/accepted_hits.bam %s'%(abspathoDir,bamname)
        cmd3 = 'ln -s %s/unmapped.bam %s'%(abspathoDir,unmapedBamName)
        f = open('job_'+oDir, 'w')
        f.write(header+'\n'+cmd1+'\n'+cmd2+'\n'+cmd3)
        f.close()
if __name__ == "__main__":
    import sys
    if len(sys.argv) == 3:
        getjobs(sys.argv[1], sys.argv[2])
    else:
        print 'Usage:\npython tophatOnCluster.py refseq gtf'

