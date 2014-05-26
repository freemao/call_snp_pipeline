#!/sur/lib/python
#-*- cdding:utf-8 -*-
f0 = open('C18.vcf', 'r')
f1 = open('C18.filteredtest.vcf', 'w')

for i in f0:
    if i.startswith('#'):
        print i
        f1.write(i)
    else:
        j = i.split()
        if (float(j[5]) > 30 and j[7].split('=')[-1] == 'snp'
            and int(j[9].split(':')[1]) > 10) :
            f1.write(i)

f0.close()
f1.close()


