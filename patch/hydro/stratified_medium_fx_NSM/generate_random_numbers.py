from __future__ import division
import sys
import math
import numpy
import scipy
import csv

fileout = 'list_of_random.dat'
fout = open(fileout, 'a')

for i in range(0,500000):
    r=numpy.random.ranf()
    x=numpy.random.ranf()
    y=numpy.random.ranf()
    z=0.4+0.2*numpy.random.ranf()
    r_NSM=numpy.random.ranf()
    r_Ia=numpy.random.ranf()
    dice=numpy.random.rand()
    #print r, x, y, z
    format = "%e\t %e\t %e\t %e\t %e\t %e\t %e\n"
    vals = (r,x,y,z,r_NSM,r_Ia,dice)
    string = format % vals
    fout.write(string)

fout.close()
