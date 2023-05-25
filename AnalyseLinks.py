import math
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import minimize

numtwists = 1.0
twist = numtwists
a,b,c = 10,-2,2.0#b is negative.a
Lx, Ly,Lz = 101, 101, 51
tolerance, a = 0.95, 0.1
# if not os.path.exists('results'):
#     os.makedirs('results')
# if not os.path.exists('results/a%fb%fc%ftwist%d' % (a,b,c,twist)):
#     os.makedirs('results/a%fb%fc%ftwist%d' % (a,b,c,twist))
# os.chdir('results/a%fb%fc%ftwist%d' % (a,b,c,twist))

director = np.loadtxt('idirectorfield.dat')
ndotzfile = open('results/a%fb%fc%ftwist%d/ndotz.dat', % (a,b,c,twist), 'w')

print(len(director))
i = 0
for z in range(0,Lz):
    for x in range(0,Lx):
        for y in range(0,Ly):
            #for i in range(0,len(director)):
            ndotz = (director[i,2])
            ndotzfile.write("%f %f  %f  %f\n" % (x, y, z, ndotz))
            i+=1
            print(i)
ndotzfile.close()
