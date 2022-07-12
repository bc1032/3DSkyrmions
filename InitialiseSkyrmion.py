import math
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import minimize
# a,b,c = 10,-10,5e-2#b is negative.a
#
# if not os.path.exists('results'):
#     os.makedirs('results')
# if not os.path.exists('results/c%f' % c):
#     os.makedirs('results/c%f' % c)
# os.chdir('results/c%f' % c)


def initialise(Lx, Ly,Lz, R,twist):
    theta = np.zeros([Lx, Ly,Lz])
    phi = np.zeros([Lx,Ly,Lz])
    centrex = int(Lx/2)
    centrey = int(Ly/2)
    file = open('directorfield.dat', 'w')
    filei = open('i.dat', 'w')
    filej = open('j.dat', 'w')
    filek = open('k.dat', 'w')
    for z in range(0,int(Lz/3)):
        for x in range(0,Lx):
            for y in range(0,Ly):
                dx = abs(centrex - x)
                dy = abs(centrey - y)

                # If inside circle of Radius, R, change elements to desired quantities.
                if math.sqrt(dx**2 + dy**2) <= R:
                    phi[x,y,z] = math.pi/2 + math.atan2(y-centrey,x-centrex)
                    theta[x,y,z] =  twist*math.pi*math.sqrt((x-centrex)**2 + (y-centrey)**2)/(R)
                else:
                    phi[x,y,z] = math.pi/2 #+ math.atan2(y-centrey,x-centrex)
                    theta[x,y,z] = math.pi

                i = math.sin(theta[x,y,z])*math.cos(phi[x,y,z])
                j = math.sin(theta[x,y,z])*math.sin(phi[x,y,z])
                k = math.cos(theta[x,y,z])
                file.write("%f  %f  %f\n" % (i,j,k))
                filei.write("%f\n" % (i))
                filej.write("%f\n" % (j))
                filek.write("%f\n" % (k))
    for z in range(int(Lz/3),int(2.*Lz/3)):
        for x in range(0,Lx):
            for y in range(0,Ly):
                dx = abs(centrex - x)
                dy = abs(centrey - y)

                #No skyrmion in middle slice
                phi[x,y,z] = math.pi/2 #+ math.atan2(y-centrey,x-centrex)
                theta[x,y,z] = math.pi

                i = math.sin(theta[x,y,z])*math.cos(phi[x,y,z])
                j = math.sin(theta[x,y,z])*math.sin(phi[x,y,z])
                k = math.cos(theta[x,y,z])
                file.write("%f  %f  %f\n" % (i,j,k))
                filei.write("%f\n" % (i))
                filej.write("%f\n" % (j))
                filek.write("%f\n" % (k))

    for z in range(int(2.*Lz/3),int(Lz)):
        for x in range(0,Lx):
            for y in range(0,Ly):
                dx = abs(centrex - x)
                dy = abs(centrey - y)

                # If inside circle of Radius, R, change elements to desired quantities.
                if math.sqrt(dx**2 + dy**2) <= R:
                    phi[x,y,z] = math.pi/2 + math.atan2(y-centrey,x-centrex)
                    theta[x,y,z] =  twist*math.pi*math.sqrt((x-centrex)**2 + (y-centrey)**2)/(R)
                else:
                    phi[x,y,z] = math.pi/2 #+ math.atan2(y-centrey,x-centrex)
                    theta[x,y,z] = math.pi

                i = math.sin(theta[x,y,z])*math.cos(phi[x,y,z])
                j = math.sin(theta[x,y,z])*math.sin(phi[x,y,z])
                k = math.cos(theta[x,y,z])
                file.write("%f  %f  %f\n" % (i,j,k))
                filei.write("%f\n" % (i))
                filej.write("%f\n" % (j))
                filek.write("%f\n" % (k))

    #print(phi,theta)
    #np.savetxt('phi.dat', phi)
    #np.savetxt('theta.dat', theta)
    #np.savetxt('iphi.dat', phi)
    #np.savetxt('itheta.dat', theta)
    file.close()
    filei.close()
    filej.close()
    filek.close()

    return(phi,theta)
