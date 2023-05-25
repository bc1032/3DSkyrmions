import math
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import minimize
from numba import jit
import pandas as pd
import glob
import re

directory = 'newInitialise/'
#directory = 'startingfromnewintialised/'
#os.chdir('/home/bc1032/Desktop/Work/Cholesterics/Skyrmions/DefectInitialisation/')
#os.chdir('/media/bc1032/Seagate Portable Drive/Skyrmions/')

# all_dirs = glob.glob(os.path.join(directory, '*/'))

#exclude_file = 'gapvstime.dat'

# Use glob to get a list of all directories in the directory
all_dirs = glob.glob(os.path.join(directory, '*dz1.0*/'))

# Use a list comprehension to filter out the directories that contain the exclude_file
#all_dirs = [d for d in all_dirsinitial if not os.path.isfile(os.path.join(d, exclude_file))]

minen = []
b = []
dz = 1.0
Lx,Ly,Lz=51,51,int(51/dz)
cutoff = 0.0
lastfile = 10000

def countgap():
    for f in range (0,25*lastfile,25):
        file = pd.read_csv('finaldirectorfield%d.csv' % f)

        i = 0
        for z in range(0,Lz):
            for x in range(26,27):
                for y in range(26,27):
                    vz = file.vz[z*Lx*Ly + x*Ly + y]
                    if vz < cutoff:
                        i+=1
    return(i)

for dir in range(0,len(all_dirs)):
    #energy = pd.read_csv('%sEnergy.csv' % all_dirs[dir])
    string = all_dirs[dir]
    # Define the regular expression pattern for float values
    pattern = r'\d+\.\d+'

    # Use re.findall to obtain all float values in the string
    float_values = re.findall(pattern, string)
    ks = float(float_values[0])
    kb = float(float_values[1])
    kt = float(float_values[2])
    b = -float(float_values[3])
    q0 = float(float_values[4])
    # dx = float(float_values[5])
    # dy = float(float_values[6])
    # dz = float(float_values[7])
    #dz = 1
    print(dz)
    #pd.read_csv('%s/Energy.csv' % all_dirs[dir])
    path = ('%s/finaldirectorfield' % all_dirs[dir])

    lst = os.listdir(path) # your directory path
    lastfile = len(lst)
    print(lastfile*25)
    #os.chdir('results/ks1.000000kb1.000000kt0.500000b-0.400000c0.000000twist3.000000/finaldirectorfield')
    gapfile = open('%sgapvstime.dat' % all_dirs[dir], 'w')
    gapfile.write("gap\n")

    #file = pandas.read_csv('finaldirectorfield3975.csv')
    #countgap()
    for f in range (0,(25*lastfile),25):
        if len(lst) > 1:
            file = pd.read_csv('%sfinaldirectorfield/finaldirectorfield%d.csv' % (all_dirs[dir],f))

        else:
            file = pd.read_csv('%sfinaldirectorfield/finaldirectorfield.csv' % all_dirs[dir])

        i = 0
        for z in range(0,Lz):
            for x in range(26,27):
                for y in range(26,27):
                    vz = file.vz[z*Lx*Ly + x*Ly + y]
                    if vz < cutoff:
                        i+=1
        gapfile.write("%f\n" % (i*dz))
    print(i,dz)
    print(f)

    gapfile.close()
