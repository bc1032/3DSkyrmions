import math
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import minimize

if not os.path.exists('results'):
    os.makedirs('results')

os.chdir('results')


Lx, Ly,Lz = 101, 101, 51
R = 40
numtwists,numtwiststwo = 1.0, 2.0
twist = 1.0/numtwists
twisttwo = 1.0/numtwists

def initialise(Lx, Ly, R,twist):
    theta = np.zeros([Lx, Ly])
    phi = np.zeros([Lx,Ly])
    centrex = int(Lx/2)
    centrey = int(Ly/2)
    file = open('directorfield.dat', 'w')
    filei = open('i.dat', 'w')
    filej = open('j.dat', 'w')
    filek = open('k.dat', 'w')
    for z in range(0,Lz):
        for x in range(0,Lx):
            for y in range(0,Ly):
                dx = abs(centrex - x)
                dy = abs(centrey - y)
                if z <= int(Lz/3.0) or z >= int(2.0*Lz/3.0):
                # If inside circle of Radius, R, change elements to desired quantities.
                    if math.sqrt(dx**2 + dy**2) <= R:
                        phi[x,y] = math.pi/2 + math.atan2(y-centrey,x-centrex)
                        theta[x,y] =  math.pi*math.sqrt((x-centrex)**2 + (y-centrey)**2)/(twist*R)
                    else:
                        phi[x,y] = math.pi/2
                        theta[x,y] = math.pi
                else:
                    if math.sqrt(dx**2 + dy**2) <= R:
                        phi[x,y] = math.pi/2 + math.atan2(y-centrey,x-centrex)
                        theta[x,y] =  math.pi*math.sqrt((x-centrex)**2 + (y-centrey)**2)/(twisttwo*R)
                    else:
                        phi[x,y] = math.pi/2
                        theta[x,y] = math.pi
                i = math.sin(theta[x,y])*math.cos(phi[x,y])
                j = math.sin(theta[x,y])*math.sin(phi[x,y])
                k = math.cos(theta[x,y])
                file.write("%f  %f  %f\n" % (i,j,k))
                filei.write("%f\n" % (i))
                filej.write("%f\n" % (j))
                filek.write("%f\n" % (k))


    #print(phi,theta)
    np.savetxt('phi.dat', phi)
    np.savetxt('theta.dat', theta)
    file.close()
    filei.close()
    filej.close()
    filek.close()

    return(phi,theta)

initialise(Lx,Ly,R,twist)
