import math
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
import pandas
from scipy.optimize import minimize

Lx, Ly,Lz = 51, 51, 51
R = 20.0
numtwists = 3
twist = numtwists
dx,dy,dz = 1,1,1

ks,kb,kt,b = .3,1.0,1.0,-0.2

filephi = pandas.read_csv('/home/bc1032/Desktop/Work/Cholesterics/Skyrmions/3DCodeFull/results/ks%fkb%fkt%fb%fc0.000000twist3.000000/phi/phi1075.csv' % (ks,kb,kt,b))
filetheta = pandas.read_csv('/home/bc1032/Desktop/Work/Cholesterics/Skyrmions/3DCodeFull/results/ks%fkb%fkt%fb%fc0.000000twist3.000000/theta/theta1075.csv' % (ks,kb,kt,b))
theta = np.zeros([Lx,Ly,Lz])
phi = np.zeros([Lx,Ly,Lz])
twist = np.zeros([Lx,Ly,Lz])

#print(filephi.phi)
i = 0
for z in range(0,int(Lz)):
    for x in range(0,Lx):
        for y in range(0,Ly):
            phi[x,y,z] = filephi.phi[i]
            theta[x,y,z] = filetheta.theta[i]
            i+=1

for z in range(0,Lz):

    if z == 0:
            zb = Lz-1
            za = 1
    elif z == Lz-1:
        zb = Lz-2
        za = 0
    else:
        zb = z - 1
        za = z + 1

    for x in range(0,Lx):

        if x == 0:
            xl = Lx-1
            xr = 1
        elif x == Lx-1:
            xl = Lx-2
            xr = 0
        else:
            xl = x - 1
            xr = x + 1

        for y in range(0,Ly):

            if y == 0:
                ya = y+1
                yb = Ly-1
            elif y == Ly-1:
                ya = 0
                yb = Ly-2
            else:
                ya = y+1
                yb = y-1

            dthetadz = (theta[x,y,za]-theta[x,y,zb])/(2.*dx)
            dthetadzsq = (theta[x,y,za]+theta[x,y,zb] - 2.*theta[x,y,z])/(dx)
            dthetadx = (theta[xr,y,z]-theta[xl,y,z])/(2.*dx)
            dthetady = (theta[x,ya,z]-theta[x,yb,z])/(2.*dx)
            dthetadxsq = (theta[xr,y,z]+theta[xl,y,z] - 2.*theta[x,y,z])/(dx)
            dthetadysq = (theta[x,ya,z]+theta[x,yb,z] - 2.*theta[x,y,z])/(dx)
            dphidz = (phi[x,y,za]-phi[x,y,zb])/(2.*dx)
            dphidzsq = (phi[x,y,za]+phi[x,y,zb] - 2.*phi[x,y,z])/(dx)

            dthetadxdy = (theta[xr,ya,z] + theta[xl,yb,z] - theta[xl,ya,z] - theta[xr,yb,z])/(4*dx*dy)
            dthetadydz = (theta[x,ya,za] + theta[x,yb,zb] - theta[x,ya,zb] - theta[x,yb,za])/(4*dy*dz)
            dthetadxdz = (theta[xr,y,za] + theta[xl,y,zb] - theta[xl,y,za] - theta[xr,y,zb])/(4*dx*dz)

            dphidx1 = (phi[xr,y,z]-phi[xl,y,z])/(2.*dx)
            dphidxm = (phi[xr,y,z]-phi[xl,y,z] - 2*math.pi)/(2.*dx)
            dphidxp = (phi[xr,y,z]-phi[xl,y,z] + 2*math.pi)/(2.*dx)

            if abs(dphidx1) < abs(dphidxm) and abs(dphidx1) < abs(dphidxp):
                dphidx = dphidx1
            elif abs(dphidxm) < abs(dphidx1) and abs(dphidxm) < abs(dphidxp):
                dphidx = dphidxm
            else:
                dphidx = dphidxp
            #dphidx = min(abs(dphidx1),abs(dphidxm),abs(dphidxp))

            dphidy1 = (phi[x,ya,z]-phi[x,yb,z])/(2.*dx)
            dphidym = (phi[x,ya,z]-phi[x,yb,z] - 2*math.pi)/(2.*dx)
            dphidyp = (phi[x,ya,z]-phi[x,yb,z] + 2*math.pi)/(2.*dx)

            if abs(dphidy1) < abs(dphidym) and abs(dphidy1) < abs(dphidyp):
                dphidy = dphidy1
            elif abs(dphidym) < abs(dphidy1) and abs(dphidym) < abs(dphidyp):
                dphidy = dphidym
            else:
                dphidy = dphidyp
            #dphidy = min(abs(dphidy1),abs(dphidym),abs(dphidyp))

            dphidxdy1 = (phi[xr,ya,z] + phi[xl,yb,z] - phi[xl,ya,z] - phi[xr,yb,z])/(4*dx*dy)
            dphidxdyp = (phi[xr,ya,z] + phi[xl,yb,z] - phi[xl,ya,z] - phi[xr,yb,z] + 2*math.pi)/(4*dx*dy)
            dphidxdym = (phi[xr,ya,z] + phi[xl,yb,z] - phi[xl,ya,z] - phi[xr,yb,z] - 2*math.pi)/(4*dx*dy)

            if abs(dphidxdy1) < abs(dphidxdym) and abs(dphidxdy1) < abs(dphidxdyp):
                dphidxdy = dphidxdy1
            elif abs(dphidxdym) < abs(dphidxdy1) and abs(dphidxdym) < abs(dphidxdyp):
                dphidxdy = dphidxdym
            else:
                dphidxdy = dphidxdyp

            dphidxdz1 = (phi[xr,y,za] + phi[xl,y,zb] - phi[xl,y,za] - phi[xr,y,zb])/(4*dx*dz)
            dphidxdzp = (phi[xr,y,za] + phi[xl,y,zb] - phi[xl,y,za] - phi[xr,y,zb] + 2*math.pi)/(4*dx*dz)
            dphidxdzm = (phi[xr,y,za] + phi[xl,y,zb] - phi[xl,y,za] - phi[xr,y,zb] - 2*math.pi)/(4*dx*dz)

            if abs(dphidxdz1) < abs(dphidxdzm) and abs(dphidxdz1) < abs(dphidxdzp):
                dphidxdz = dphidxdz1
            elif abs(dphidxdzm) < abs(dphidxdz1) and abs(dphidxdzm) < abs(dphidxdzp):
                dphidxdz = dphidxdzm
            else:
                dphidxdz = dphidxdzp

            dphidydz1 = (phi[x,ya,za] + phi[x,yb,zb] - phi[x,yb,za] - phi[x,ya,zb])/(4*dy*dz)
            dphidydzp = (phi[x,ya,za] + phi[x,yb,zb] - phi[x,yb,za] - phi[x,ya,zb] + 2*math.pi)/(4*dy*dz)
            dphidydzm = (phi[x,ya,za] + phi[x,yb,zb] - phi[x,yb,za] - phi[x,ya,zb] - 2*math.pi)/(4*dy*dz)

            if abs(dphidydz1) < abs(dphidydzm) and abs(dphidydz1) < abs(dphidydzp):
                dphidydz = dphidydz1
            elif abs(dphidydzm) < abs(dphidydz1) and abs(dphidydzm) < abs(dphidydzp):
                dphidydz = dphidydzm
            else:
                dphidydz = dphidydzp

            dphidxsq1 = (phi[xr,y,z]+phi[xl,y,z] - 2*phi[x,y,z])/(dx*dx)
            dphidxsqm = (phi[xr,y,z]+phi[xl,y,z] - 2*phi[x,y,z] - 2*math.pi)/(dx*dx)
            dphidxsqp = (phi[xr,y,z]+phi[xl,y,z] - 2*phi[x,y,z] + 2*math.pi)/(dx*dx)

            if abs(dphidxsq1) < abs(dphidxsqm) and abs(dphidxsq1) < abs(dphidxsqp):
                dphidxsq = dphidxsq1
            elif abs(dphidxsqm) < abs(dphidxsq1) and abs(dphidxsqm) < abs(dphidxsqp):
                dphidxsq = dphidxsqm
            else:
                dphidxsq = dphidxsqp

            #dphidxsq = min(abs(dphidxsq1),abs(dphidxsqm),abs(dphidxsqp))

            dphidysq1 = (phi[x,ya,z]+phi[x,yb,z] - 2*phi[x,y,z])/(dy*dy)
            dphidysqm = (phi[x,ya,z]+phi[x,yb,z] - 2*phi[x,y,z] - 2*math.pi)/(dy*dy)
            dphidysqp = (phi[x,ya,z]+phi[x,yb,z] - 2*phi[x,y,z] + 2*math.pi)/(dy*dy)

            if abs(dphidysq1) < abs(dphidysqm) and abs(dphidysq1) < abs(dphidysqp):
                dphidysq = dphidysq1
            elif abs(dphidysqm) < abs(dphidysq1) and abs(dphidysqm) < abs(dphidysqp):
                dphidysq = dphidysqm
            else:
                dphidysq = dphidysqp


            twist[x,y,z] = (-2*dphidz*math.sin(theta[x,y,z])**2 + math.cos(phi[x,y,z])*(-2*dthetady + dphidx*math.sin(2*theta[x,y,z]))\
                + (2*dthetadx + dphidy*math.sin(2*theta[x,y,z]))*math.sin(phi[x,y,z]))/2.

file = open('/home/bc1032/Desktop/Work/Cholesterics/Skyrmions/3DCodeFull/results/ks%fkb%fkt%fb%fc0.000000twist3.000000/twisten.csv' % (ks,kb,kt,b), 'w')
file.write('x,y,z,twisten\n')


for z in range(0,Lz):
    for x in range(0,Lx):
        for y in range(0,Ly):
            file.write("%d,%d,%d,%f\n" % (x,y,z,twist[x,y,z]))
file.close()
