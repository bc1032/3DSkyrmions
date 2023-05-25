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
# if not os.path.exists('results/c%f' % c):import math
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
    split = 4
    phi = np.zeros([Lx,Ly,Lz])
    centrex = int(Lx/2)
    centrey = int(Ly/2)
    file = open('idirectorfield.csv', 'w')
    file2 = open('idirectorfield.dat', 'w')
    file.write('x,y,z,vx,vy,vz\n')
    filei = open('i.dat', 'w')
    filej = open('j.dat', 'w')
    filek = open('k.dat', 'w')
    count = 0
    for z in range(0,int(Lz/3)):
        for x in range(0,Lx):
            for y in range(0,Ly):
                dx = abs(centrex - x)
                dy = abs(centrey - y)

                # If inside circle of Radius, R, change elements to desired quantities.
                if math.sqrt(dx**2 + dy**2) <= R:
                    phi[x,y,z] = math.pi/2 + math.atan2(y-centrey,x-centrex)# + 2.*math.pi*z/(Lz-1)#  + math.pi/2 (f(Z)) instead
                    theta[x,y,z] =  twist*(math.pi/2.0)*math.sqrt((x-centrex)**2 + (y-centrey)**2)/(R)
                else:
                    phi[x,y,z] = math.pi/2 + math.atan2(y-centrey,x-centrex)
                    theta[x,y,z] =  twist*math.pi/2

                i = math.sin(theta[x,y,z])*math.cos(phi[x,y,z])
                j = math.sin(theta[x,y,z])*math.sin(phi[x,y,z])
                k = math.cos(theta[x,y,z])
                file.write("%d,%d,%d,%f,%f,%f\n" % (x,y,z,i,j,k))
                file2.write("%f  %f  %f\n" % (i,j,k))
                filei.write("%f\n" % (i))
                filej.write("%f\n" % (j))
                filek.write("%f\n" % (k))
        count+=1
    count = 0
    for z in range(int(Lz/3),int(2*Lz/3.0)):
        for x in range(0,Lx):
            for y in range(0,Ly):
                dx = abs(centrex - x)
                dy = abs(centrey - y)
                if math.sqrt(dx**2 + dy**2) <= R*(twist-1)/twist:
                    phi[x,y,z] = math.pi/2 + 2*math.pi*count/(Lz/3) + math.atan2(y-centrey,x-centrex) #  + math.pi/2 (f(Z)) instead
                    theta[x,y,z] =  twist*(math.pi/2.0)*math.sqrt((x-centrex)**2 + (y-centrey)**2)/(R)

                elif math.sqrt(dx**2  + dy**2) <= R and math.sqrt(dx**2 + dy**2) > R*(twist-1)/twist:
                    phi[x,y,z] = math.pi/2 + math.atan2(y-centrey,x-centrex) #  + math.pi/2 (f(Z)) instead
                    theta[x,y,z] =  twist*(math.pi/2.0)*math.sqrt((x-centrex)**2 + (y-centrey)**2)/(R)

                else:
                    phi[x,y,z] = math.pi/2 + math.atan2(y-centrey,x-centrex)# - 2*math.pi#math.pi/2 - 2*math.pi*count/((Lz)/3) + math.atan2(y-centrey,x-centrex)
                    theta[x,y,z] =  twist*math.pi/2

                i = math.sin(theta[x,y,z])*math.cos(phi[x,y,z])
                j = math.sin(theta[x,y,z])*math.sin(phi[x,y,z])
                k = math.cos(theta[x,y,z])
                file.write("%d,%d,%d,%f,%f,%f\n" % (x,y,z,i,j,k))
                file2.write("%f  %f  %f\n" % (i,j,k))
                filei.write("%f\n" % (i))
                filej.write("%f\n" % (j))
                filek.write("%f\n" % (k))
        count+=1

    for z in range(int(2*Lz/3),Lz):
        for x in range(0,Lx):
            for y in range(0,Ly):
                dx = abs(centrex - x)
                dy = abs(centrey - y)

                # If inside circle of Radius, R, change elements to desired quantities.
                if math.sqrt(dx**2 + dy**2) <= R:
                    phi[x,y,z] = math.pi/2 + math.atan2(y-centrey,x-centrex)# + 2.*math.pi#*z/(Lz-1)#  + math.pi/2 (f(Z)) instead
                    theta[x,y,z] =  twist*(math.pi/2.0)*math.sqrt((x-centrex)**2 + (y-centrey)**2)/(R)
                else:
                    phi[x,y,z] = math.pi/2 + math.atan2(y-centrey,x-centrex)
                    theta[x,y,z] =  twist*math.pi/2


                i = math.sin(theta[x,y,z])*math.cos(phi[x,y,z])
                j = math.sin(theta[x,y,z])*math.sin(phi[x,y,z])
                k = math.cos(theta[x,y,z])
                file.write("%d,%d,%d,%f,%f,%f\n" % (x,y,z,i,j,k))
                file2.write("%f  %f  %f\n" % (i,j,k))
                filei.write("%f\n" % (i))
                filej.write("%f\n" % (j))
                filek.write("%f\n" % (k))
    file.close()
    file2.close()
    filei.close()
    filej.close()
    filek.close()

    return(phi,theta)
