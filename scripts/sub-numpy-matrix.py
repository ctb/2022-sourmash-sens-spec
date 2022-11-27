#! /usr/bin/env python
import numpy
import sys

with open(sys.argv[1], 'rb') as fp:
    mat1 = numpy.load(fp)

with open(sys.argv[2], 'rb') as fp:
    mat2 = numpy.load(fp)

numpy.set_printoptions(precision=2, suppress=True)
print(numpy.abs(mat1 - mat2))
