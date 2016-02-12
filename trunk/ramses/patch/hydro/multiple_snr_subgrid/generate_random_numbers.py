from __future__ import division
import sys
import math
import numpy
import scipy
import csv

fileout = 'list_of_random.dat'
fout = open(fileout, 'a')

for i in range(0,10000):
    r=numpy.random.ranf()
    x=numpy.random.ranf()
    y=numpy.random.ranf()
    z=numpy.random.ranf()
    #print r, x, y, z
    format = "%e\t %e\t %e\t %e\n"
    vals = (r,x,y,z)
    string = format % vals
    fout.write(string)

fout.close()
