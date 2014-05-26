#!/usr/lib/pyhton
#-*- coding:utf-8 -*-

from freebayes_pipeline import FreebayesPipe
#ob1 = FreebayesPipe('.')
#ob1.getrgfilelist()
#print ob1.namelist
#ob1.runfreebayesfile()

ob2 = FreebayesPipe('.')
ob2.getvcffilelist()
ob2.runfilterfile()
