import math
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import minimize
from numba import jit
import pandas
import glob
import re

import InitialiseSkyrmion

#os.chdir('/media/bc1032/Seagate Portable Drive/Skyrmions/')
dz = 1.0
Lx, Ly,Lz = 51, 51, int(51/dz)
R = 20.0
numtwists = 3
twist = numtwists
numits = 10000

tolerance = 0.001
gamma = 1e-2
#ks,kb,kt,b,c = 1.0,1.0,1.0,-0.6,0.0#b is negative
n = 0

# Set the directory you want to search in
directory = 'startingfromnewintialised/'
#directory = 'newInitialise/'
#directory = 'resultsoneconstant'

# Use glob to get a list of all directories in the directory
all_dirs = glob.glob(os.path.join(directory, '*b0.4*dz1.0*/'))

# Set the file name you want to exclude from the search
#exclude_file = 'Energy.csv'

# Use a list comprehension to filter out the directories that contain the exclude_file
#all_dirs = [d for d in all_dirsinitial if not os.path.isfile(os.path.join(d, exclude_file))]

def load(Lx, Ly, Lz,i):

    if i == 0:
        os.chdir('../phi')
        path=os.getcwd()
        lst = os.listdir(path) # your directory path
        lastfile = len(lst)
        if len(lst) > 1:
            filephi = pandas.read_csv('phi%d.csv' % (i))
        else:
            filephi = pandas.read_csv('phi.csv')

        os.chdir('../theta')
        path=os.getcwd()
        lst = os.listdir(path) # your directory path
        lastfile = len(lst)
        if len(lst) > 1:
            filetheta = pandas.read_csv('theta%d.csv' % (i))
        else:
            filetheta = pandas.read_csv('theta.csv')

    else:
        os.chdir('../phi')

        path=os.getcwd()
        lst = os.listdir(path) # your directory path
        lastfile = len(lst)
        if len(lst) > 1:
            filephi = pandas.read_csv('phi%d.csv' % (i))
        else:
            filephi = pandas.read_csv('phi.csv')

        os.chdir('../theta')

        path=os.getcwd()
        lst = os.listdir(path) # your directory path
        lastfile = len(lst)
        if len(lst) > 1:
            filetheta = pandas.read_csv('theta%d.csv' % (i))
        else:
            filetheta = pandas.read_csv('theta.csv')


    path=os.getcwd()
    lst = os.listdir(path) # your directory path
    lastfile = len(lst)

    theta = np.zeros([Lx,Ly,Lz])
    phi = np.zeros([Lx,Ly,Lz])
    #print(filephi.phi)
    d = 0
    for z in range(0,int(Lz)):
        for x in range(0,Lx):
            for y in range(0,Ly):
                phi[x,y,z] = filephi.phi[d]
                theta[x,y,z] = filetheta.theta[d]
                d+=1
    return(phi,theta)

@jit
def calcen(ks,kb,kt,b,c, phi,theta,Lx,Ly,Lz,n,i,dx,dy,dz):

    a = 1
    E,splay,bend,twist2,twist1,sadsplay,bendtwist,twist0= 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    kss = 1
    centrex = int(Lx/2)
    centrey = int(Ly/2)
    for z in range(0,Lz):

    #    for z in range(0,Lz):
        if z == 0:
            za = 1
            zb = Lz-1
        elif z == Lz-1:
            za = 0
            zb = Lz-2
        else:
            za = z+1
            zb = z-1

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

                dddx = abs(centrex - x)
                dddy = abs(centrey - y)

                dthetadz = (theta[x,y,za]-theta[x,y,zb])/(2.*dz)
                dthetadzsq = (theta[x,y,za]+theta[x,y,zb] - 2.*theta[x,y,z])/(dz**2)
                dthetadx = (theta[xr,y,z]-theta[xl,y,z])/(2.*dx)
                dthetady = (theta[x,ya,z]-theta[x,yb,z])/(2.*dy)
                dthetadxsq = (theta[xr,y,z]+theta[xl,y,z] - 2.*theta[x,y,z])/(dx*dx)
                dthetadysq = (theta[x,ya,z]+theta[x,yb,z] - 2.*theta[x,y,z])/(dy*dy)
                dphidz = (phi[x,y,za]-phi[x,y,zb])/(2.*dz)
                dphidzsq = (phi[x,y,za]+phi[x,y,zb] - 2.*phi[x,y,z])/(dz*dz)

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

                dphidy1 = (phi[x,ya,z]-phi[x,yb,z])/(2.*dy)
                dphidym = (phi[x,ya,z]-phi[x,yb,z] - 2*math.pi)/(2.*dy)
                dphidyp = (phi[x,ya,z]-phi[x,yb,z] + 2*math.pi)/(2.*dy)

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

                if math.sqrt(dddx**2 + dddy**2) <= R:

    # Energy with no saddle splay.
                    splay += ks*(np.sin(theta[x,y,z])*(dthetadz - dphidy*np.cos(phi[x,y,z]) + dphidx*np.sin(phi[x,y,z])) - np.cos(theta[x,y,z])*(dthetadx*np.cos(phi[x,y,z]) + dthetady*np.sin(phi[x,y,z])))**2

                    bend += kb*(np.sin(theta[x,y,z])**2*(dthetadz*np.cos(theta[x,y,z]) + np.sin(theta[x,y,z])*(dthetadx*np.cos(phi[x,y,z]) + dthetady*np.sin(phi[x,y,z])))**2 + \
                        (dthetadz*np.cos(theta[x,y,z])**2*np.sin(phi[x,y,z]) + np.cos(phi[x,y,z])*np.sin(theta[x,y,z])**2*(dphidx*np.cos(phi[x,y,z]) + dphidy*np.sin(phi[x,y,z])) + \
                        (np.sin(2*theta[x,y,z])*(dthetady*np.sin(phi[x,y,z])**2 + np.cos(phi[x,y,z])*(dphidz + dthetadx*np.sin(phi[x,y,z]))))/2.)**2 + \
                        (dthetadz*np.cos(theta[x,y,z])**2*np.cos(phi[x,y,z]) - np.sin(theta[x,y,z])**2*np.sin(phi[x,y,z])*(dphidx*np.cos(phi[x,y,z]) + dphidy*np.sin(phi[x,y,z])) + \
                        (np.sin(2*theta[x,y,z])*(dthetadx + dthetadx*np.cos(2*phi[x,y,z]) - 2*dphidz*np.sin(phi[x,y,z]) + dthetady*np.sin(2*phi[x,y,z])))/4.)**2)

                    sadsplay += 2*kss*np.sin(theta[x,y,z])*((-(dphidy*dthetadx) + dphidx*dthetady)*np.cos(theta[x,y,z]) + \
                        np.sin(theta[x,y,z])*((-(dphidz*dthetady) + dphidy*dthetadz)*np.cos(phi[x,y,z]) + (dphidz*dthetadx - dphidx*dthetadz)*np.sin(phi[x,y,z])))

                    twist1 += (b*(-2*dphidz*np.sin(theta[x,y,z])**2 + np.cos(phi[x,y,z])*(-2*dthetady + dphidx*np.sin(2*theta[x,y,z])) + (2*dthetadx + dphidy*np.sin(2*theta[x,y,z]))*np.sin(phi[x,y,z])))/2.

                    twist2 += (kt*(-2*dphidz*np.sin(theta[x,y,z])**2 + np.cos(phi[x,y,z])*(-2*dthetady + dphidx*np.sin(2*theta[x,y,z])) + (2*dthetadx + dphidy*np.sin(2*theta[x,y,z]))*np.sin(phi[x,y,z]))**2)/4.

                    bendtwist += (-4*dphidz*dthetadx*np.cos(theta[x,y,z])*np.sin(theta[x,y,z]) - 2*dphidxdz*np.sin(theta[x,y,z])**2 +\
                        np.cos(phi[x,y,z])*(-2*dthetadxdy + 2*dphidx*dthetadx*np.cos(2*theta[x,y,z]) + dphidxsq*np.sin(2*theta[x,y,z])) + dphidx*np.cos(phi[x,y,z])*(2*dthetadx + dphidy*np.sin(2*theta[x,y,z])) -\
                        dphidx*(-2*dthetady + dphidx*np.sin(2*theta[x,y,z]))*np.sin(phi[x,y,z]) + (2*dthetadxsq + 2*dphidy*dthetadx*np.cos(2*theta[x,y,z]) + dphidxdy*np.sin(2*theta[x,y,z]))*np.sin(phi[x,y,z]))**2/4.\
                        + (-4*dphidz*dthetadz*np.cos(theta[x,y,z])*np.sin(theta[x,y,z]) - 2*dphidzsq*np.sin(theta[x,y,z])**2 +\
                        np.cos(phi[x,y,z])*(-2*dthetadydz + 2*dphidx*dthetadz*np.cos(2*theta[x,y,z]) + dphidxdz*np.sin(2*theta[x,y,z])) + dphidz*np.cos(phi[x,y,z])*(2*dthetadx + dphidy*np.sin(2*theta[x,y,z])) -\
                        dphidz*(-2*dthetady + dphidx*np.sin(2*theta[x,y,z]))*np.sin(phi[x,y,z]) + (2*dthetadxdz + 2*dphidy*dthetadz*np.cos(2*theta[x,y,z]) + dphidydz*np.sin(2*theta[x,y,z]))*np.sin(phi[x,y,z]))**2/4.\
                        + (-4*dphidz*dthetady*np.cos(theta[x,y,z])*np.sin(theta[x,y,z]) - 2*dphidydz*np.sin(theta[x,y,z])**2 +\
                        np.cos(phi[x,y,z])*(-2*dthetadysq + 2*dphidx*dthetady*np.cos(2*theta[x,y,z]) + dphidxdy*np.sin(2*theta[x,y,z])) + dphidy*np.cos(phi[x,y,z])*(2*dthetadx + dphidy*np.sin(2*theta[x,y,z])) -\
                        dphidy*(-2*dthetady + dphidx*np.sin(2*theta[x,y,z]))*np.sin(phi[x,y,z]) + (2*dthetadxdy + 2*dphidy*dthetady*np.cos(2*theta[x,y,z]) + dphidysq*np.sin(2*theta[x,y,z]))*np.sin(phi[x,y,z]))**2/4.

                    twist0 += b**2/(4.0*kt)

    E = splay + bend + twist1 + twist2 + twist0

    print(E, splay, bend, sadsplay, twist1,twist2,bendtwist,twist0)
    return(E, splay, bend, sadsplay, twist1,twist2,bendtwist,twist0)

@jit
def calcendensity(ks,kb,kt,b,c, phi,theta,Lx,Ly,Lz,n,i,dx,dy,dz):
    a = 1
    E,splay,bend,twist2,twist1,sadsplay,bendtwist,twist0 = np.zeros((Lx,Ly,Lz)),np.zeros((Lx,Ly,Lz)),np.zeros((Lx,Ly,Lz)),np.zeros((Lx,Ly,Lz)),np.zeros((Lx,Ly,Lz)),np.zeros((Lx,Ly,Lz)),np.zeros((Lx,Ly,Lz)),np.zeros((Lx,Ly,Lz))
    kss = 1
    centrex = int(Lx/2)
    centrey = int(Ly/2)
    for z in range(0,Lz):

        #for z in range(0,Lz):
        if z == 0:
            za = 1
            zb = Lz-1
        elif z == Lz-1:
            za = 0
            zb = Lz-2
        else:
            za = z+1
            zb = z-1

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

                dddx = abs(centrex - x)
                dddy = abs(centrey - y)

                dthetadz = (theta[x,y,za]-theta[x,y,zb])/(2.*dz)
                dthetadzsq = (theta[x,y,za]+theta[x,y,zb] - 2.*theta[x,y,z])/(dz**2)
                dthetadx = (theta[xr,y,z]-theta[xl,y,z])/(2.*dx)
                dthetady = (theta[x,ya,z]-theta[x,yb,z])/(2.*dy)
                dthetadxsq = (theta[xr,y,z]+theta[xl,y,z] - 2.*theta[x,y,z])/(dx*dx)
                dthetadysq = (theta[x,ya,z]+theta[x,yb,z] - 2.*theta[x,y,z])/(dy*dy)
                dphidz = (phi[x,y,za]-phi[x,y,zb])/(2.*dz)
                dphidzsq = (phi[x,y,za]+phi[x,y,zb] - 2.*phi[x,y,z])/(dz*dz)

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

                dphidy1 = (phi[x,ya,z]-phi[x,yb,z])/(2.*dy)
                dphidym = (phi[x,ya,z]-phi[x,yb,z] - 2*math.pi)/(2.*dy)
                dphidyp = (phi[x,ya,z]-phi[x,yb,z] + 2*math.pi)/(2.*dy)

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
                if math.sqrt(dddx**2 + dddy**2) <= R:

    # Energy with no saddle splay.
                    splay[x,y,z] = ks*(np.sin(theta[x,y,z])*(dthetadz - dphidy*np.cos(phi[x,y,z]) + dphidx*np.sin(phi[x,y,z])) - np.cos(theta[x,y,z])*(dthetadx*np.cos(phi[x,y,z]) + dthetady*np.sin(phi[x,y,z])))**2
                #    print(splay[x,y,z])
                    bend[x,y,z] = kb*(np.sin(theta[x,y,z])**2*(dthetadz*np.cos(theta[x,y,z]) + np.sin(theta[x,y,z])*(dthetadx*np.cos(phi[x,y,z]) + dthetady*np.sin(phi[x,y,z])))**2 + \
                        (dthetadz*np.cos(theta[x,y,z])**2*np.sin(phi[x,y,z]) + np.cos(phi[x,y,z])*np.sin(theta[x,y,z])**2*(dphidx*np.cos(phi[x,y,z]) + dphidy*np.sin(phi[x,y,z])) + \
                        (np.sin(2*theta[x,y,z])*(dthetady*np.sin(phi[x,y,z])**2 + np.cos(phi[x,y,z])*(dphidz + dthetadx*np.sin(phi[x,y,z]))))/2.)**2 + \
                        (dthetadz*np.cos(theta[x,y,z])**2*np.cos(phi[x,y,z]) - np.sin(theta[x,y,z])**2*np.sin(phi[x,y,z])*(dphidx*np.cos(phi[x,y,z]) + dphidy*np.sin(phi[x,y,z])) + \
                        (np.sin(2*theta[x,y,z])*(dthetadx + dthetadx*np.cos(2*phi[x,y,z]) - 2*dphidz*np.sin(phi[x,y,z]) + dthetady*np.sin(2*phi[x,y,z])))/4.)**2)

                    sadsplay[x,y,z] = 2*kss*np.sin(theta[x,y,z])*((-(dphidy*dthetadx) + dphidx*dthetady)*np.cos(theta[x,y,z]) + \
                        np.sin(theta[x,y,z])*((-(dphidz*dthetady) + dphidy*dthetadz)*np.cos(phi[x,y,z]) + (dphidz*dthetadx - dphidx*dthetadz)*np.sin(phi[x,y,z])))

                    twist1[x,y,z] = (b*(-2*dphidz*np.sin(theta[x,y,z])**2 + np.cos(phi[x,y,z])*(-2*dthetady + dphidx*np.sin(2*theta[x,y,z])) + (2*dthetadx + dphidy*np.sin(2*theta[x,y,z]))*np.sin(phi[x,y,z])))/2.

                    twist2[x,y,z] = (kt*(-2*dphidz*np.sin(theta[x,y,z])**2 + np.cos(phi[x,y,z])*(-2*dthetady + dphidx*np.sin(2*theta[x,y,z])) + (2*dthetadx + dphidy*np.sin(2*theta[x,y,z]))*np.sin(phi[x,y,z]))**2)/4.

                    twist0[x,y,z] = (b**2/(4.0*kt))

                    bendtwist[x,y,z] = (-4*dphidz*dthetadx*np.cos(theta[x,y,z])*np.sin(theta[x,y,z]) - 2*dphidxdz*np.sin(theta[x,y,z])**2 +\
                        np.cos(phi[x,y,z])*(-2*dthetadxdy + 2*dphidx*dthetadx*np.cos(2*theta[x,y,z]) + dphidxsq*np.sin(2*theta[x,y,z])) + dphidx*np.cos(phi[x,y,z])*(2*dthetadx + dphidy*np.sin(2*theta[x,y,z])) -\
                        dphidx*(-2*dthetady + dphidx*np.sin(2*theta[x,y,z]))*np.sin(phi[x,y,z]) + (2*dthetadxsq + 2*dphidy*dthetadx*np.cos(2*theta[x,y,z]) + dphidxdy*np.sin(2*theta[x,y,z]))*np.sin(phi[x,y,z]))**2/4.\
                        + (-4*dphidz*dthetadz*np.cos(theta[x,y,z])*np.sin(theta[x,y,z]) - 2*dphidzsq*np.sin(theta[x,y,z])**2 +\
                        np.cos(phi[x,y,z])*(-2*dthetadydz + 2*dphidx*dthetadz*np.cos(2*theta[x,y,z]) + dphidxdz*np.sin(2*theta[x,y,z])) + dphidz*np.cos(phi[x,y,z])*(2*dthetadx + dphidy*np.sin(2*theta[x,y,z])) -\
                        dphidz*(-2*dthetady + dphidx*np.sin(2*theta[x,y,z]))*np.sin(phi[x,y,z]) + (2*dthetadxdz + 2*dphidy*dthetadz*np.cos(2*theta[x,y,z]) + dphidydz*np.sin(2*theta[x,y,z]))*np.sin(phi[x,y,z]))**2/4.\
                        + (-4*dphidz*dthetady*np.cos(theta[x,y,z])*np.sin(theta[x,y,z]) - 2*dphidydz*np.sin(theta[x,y,z])**2 +\
                        np.cos(phi[x,y,z])*(-2*dthetadysq + 2*dphidx*dthetady*np.cos(2*theta[x,y,z]) + dphidxdy*np.sin(2*theta[x,y,z])) + dphidy*np.cos(phi[x,y,z])*(2*dthetadx + dphidy*np.sin(2*theta[x,y,z])) -\
                        dphidy*(-2*dthetady + dphidx*np.sin(2*theta[x,y,z]))*np.sin(phi[x,y,z]) + (2*dthetadxdy + 2*dphidy*dthetady*np.cos(2*theta[x,y,z]) + dphidysq*np.sin(2*theta[x,y,z]))*np.sin(phi[x,y,z]))**2/4.

                    E[x,y,z] = splay[x,y,z] + bend[x,y,z] + twist1[x,y,z] + twist2[x,y,z] + twist0[x,y,z]

                #endensfile.write('%d, %d, %d, %f, %f, %f, %f, %f, %f, %f, %f\n' % (x,y,z, E[x,y,z], splay[x,y,z],\
                                    #bend[x,y,z], sadsplay[x,y,z], twist1[x,y,z],twist2[x,y,z],bendtwist[x,y,z],twist0[x,y,z]))

            #    print(E[x,y,z])
    #print(E, splay, bend, sadsplay, twist1,twist2,bendtwist,twist0)
    return(E, splay, bend, sadsplay, twist1,twist2,bendtwist,twist0)


#@jit
def analyse(ks,kb,kt,b,c,q0,dir,dx,dy,dz):
    #for n in range(0,len(dirs_without_file)):

    if dir == 0:
        os.chdir('%sks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f/phi' % (directory,ks,kb,kt,abs(b),q0,dx,dy,dz,twist))
        #os.chdir('%sks%fkb%fkt%fb%fq0%ftwist3.000000/phi' % (directory,ks,kb,kt,abs(b),q0))
        path=os.getcwd()
        lst = os.listdir(path) # your directory path
        lastfile = len(lst)
    else:
        os.chdir('../../ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f/phi' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist))
        path=os.getcwd()
        lst = os.listdir(path) # your directory path
        lastfile = len(lst)
        print(path)
        #phi, theta = load(Lx,Ly,Lz,i)
    enfile = open('../Energy.csv', 'w')
    enfile.write('E,s,b,ss,t1,t2,bt,t0\n')

    if not os.path.exists('../energydensity'):
        os.makedirs('../energydensity')

    for i in range(0,(lastfile)*25,25):
        phi, theta = load(Lx,Ly,Lz,i)

        E, splay, bend, sadsplay, twist1,twist2,bendtwist,twist0 = calcen(ks,kb,kt,b,c, phi,theta,Lx,Ly,Lz,n,i,dx,dy,dz)
        enfile.write('%f, %f, %f, %f, %f, %f, %f, %f\n' % (E, splay, bend, sadsplay, twist1,twist2,bendtwist,twist0))

        endensfile = open('../energydensity/finalenergydensity%d.csv' % i,'w')
        endensfile.write('x,y,z,E,s,b,ss,t1,t2,bt,t0\n')
        E, splay, bend, sadsplay, twist1,twist2,bendtwist,twist0 = calcendensity(ks,kb,kt,b,c, phi,theta,Lx,Ly,Lz,n,i,dx,dy,dz)

        for z in range(0,Lz):
            for x in range(0,Lx):
                for y in range(0,Ly):
        #            print(E[x,y,z])
                    endensfile.write('%d, %d, %d, %f, %f, %f, %f, %f, %f, %f, %f\n' % (x,y,z, E[x,y,z], splay[x,y,z], bend[x,y,z], sadsplay[x,y,z], twist1[x,y,z],twist2[x,y,z],bendtwist[x,y,z],twist0[x,y,z]))
        endensfile.close()

        print(i)

#    enfile.close()

for dir in range(0,len(all_dirs)):
    # Print the list of directories without the exclude_file
    print(all_dirs[dir])
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
    dx = float(float_values[5])
    dy = float(float_values[6])
    dz = float(float_values[7])
    Lz = int(51/dz)
    c = 0#float(float_values[4])
    # Print the list of float values
    print(ks,kb,kt,b, dx,dy,dz)
    analyse(ks,kb,kt,b,c,q0,dir,dx,dy,dz)
