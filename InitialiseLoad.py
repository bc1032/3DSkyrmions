import math
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
import pandas
from scipy.optimize import minimize
from numba import jit

#@jit
def load(Lx, Ly, Lz):
    filephi = pandas.read_csv('/home/bc1032/Desktop/Work/Cholesterics/Skyrmions/DefectInitialisation/ResultsPositive/a1.000000bminus0.200000c0.000000twist3.000000/phi/phi4300.csv')
    filetheta = pandas.read_csv('/home/bc1032/Desktop/Work/Cholesterics/Skyrmions/DefectInitialisation/ResultsPositive/a1.000000bminus0.200000c0.000000twist3.000000/theta/theta4300.csv')

    theta = np.zeros([Lx,Ly,Lz])
    phi = np.zeros([Lx,Ly,Lz])
    #print(filephi.phi)
    i = 0
    for z in range(0,int(Lz)):
        for x in range(0,Lx):
            for y in range(0,Ly):
                phi[x,y,z] = filephi.phi[i]
                theta[x,y,z] = filetheta.theta[i]
                i+=1
                #print(filephi.phi[i])
    #print(theta[0,0,0],theta[25,8,5])
    return(phi,theta)

def loadnew(Lx, Ly, Lz):
    #filephi = pandas.read_csv('/home/bc1032/Desktop/Work/Cholesterics/Skyrmions/DefectInitialisation/newInitialise/ks1.000000kb1.000000kt1.000000b0.300000q00.150000dx1.000000dy1.000000dz0.500000twist3.000000/phi/phi7700.csv')
    #filetheta = pandas.read_csv('/home/bc1032/Desktop/Work/Cholesterics/Skyrmions/DefectInitialisation/newInitialise/ks1.000000kb1.000000kt1.000000b0.300000q00.150000dx1.000000dy1.000000dz0.500000twist3.000000/theta/theta7700.csv')
    filephi = pandas.read_csv('/home/bc1032/Desktop/Work/Cholesterics/Skyrmions/DefectInitialisation/newInitialise/ks1.000000kb1.000000kt1.000000b0.300000q00.150000dx1.000000dy1.000000dz1.000000twist3.000000/phi/phi1875.csv')
    filetheta = pandas.read_csv('/home/bc1032/Desktop/Work/Cholesterics/Skyrmions/DefectInitialisation/newInitialise/ks1.000000kb1.000000kt1.000000b0.300000q00.150000dx1.000000dy1.000000dz1.000000twist3.000000/theta/theta1875.csv')

    theta = np.zeros([Lx,Ly,Lz])
    phi = np.zeros([Lx,Ly,Lz])
    #print(filephi.phi)
    i = 0
    for z in range(0,int(Lz)):
        for x in range(0,Lx):
            for y in range(0,Ly):
                phi[x,y,z] = filephi.phi[i]
                theta[x,y,z] = filetheta.theta[i]
                i+=1
                #print(filephi.phi[i])
    #print(theta[0,0,0],theta[25,8,5])
    return(phi,theta)

def loadrestart(Lx, Ly, Lz):
    filephi = pandas.read_csv('/home/bc1032/Desktop/Work/Cholesterics/Skyrmions/DefectInitialisation/resultsnosaddlesplay/ks1.000000kb1.000000kt0.300000b-0.600000c0.000000twist3.000000/phi/phi9975.csv')
    filetheta = pandas.read_csv('/home/bc1032/Desktop/Work/Cholesterics/Skyrmions/DefectInitialisation/resultsnosaddlesplay/ks1.000000kb1.000000kt0.300000b-0.600000c0.000000twist3.000000/theta/theta9975.csv')

    theta = np.zeros([Lx,Ly,Lz])
    phi = np.zeros([Lx,Ly,Lz])
    #print(filephi.phi)
    i = 0
    for z in range(0,int(Lz)):
        for x in range(0,Lx):
            for y in range(0,Ly):
                phi[x,y,z] = filephi.phi[i]
                theta[x,y,z] = filetheta.theta[i]
                i+=1
                #print(filephi.phi[i])
    #print(theta[0,0,0],theta[25,8,5])
    return(phi,theta)

#load(Lx,Ly,Lz)
