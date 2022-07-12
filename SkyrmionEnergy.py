import math
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import minimize

import InitialiseSkyrmion

# if not os.path.exists('results'):
#     os.makedirs('results')
#
# os.chdir('results')

def sphere(a,b,c,phi,theta,Lx,Ly,Lz):
    dx = 1
    dy = 1
    E=0.0
    for z in range(0,Lz):

        for z in range(0,Lz):
            if z == 0:
                za = z+1
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

                dthetadz = (theta[x,y,za]-theta[x,y,zb])/(2.*dx)
                dthetadx = (theta[xr,y,z]-theta[xl,y,z])/(2.*dx)
                dthetady = (theta[x,ya,z]-theta[x,yb,z])/(2.*dx)
                dthetadxsq = (theta[xr,y,z]+theta[xl,y,z] - 2.*theta[x,y])/(dx)
                dthetadysq = (theta[x,ya,z]+theta[x,yb,z] - 2.*theta[x,y])/(dx)

                dphidz = (phi[x,y,za]-phi[x,y,zb])/(2.*dx)
                dphidx1 = (phi[xr,y,z]-phi[xl,y,z])/(2.*dx)
                dphidxm = (phi[xr,y,z]-phi[xl,y,z] - 2*math.pi)/(2.*dx)
                dphidxp = (phi[xr,y,z]-phi[xl,y,z] + 2*math.pi)/(2.*dx)

                if abs(dphidx1) < abs(dphidxm) and abs(dphidx1) < abs(dphidxp):
                    dphidx = dphidx1
                elif abs(dphidxm) < abs(dphidx1) and abs(dphidxm) < abs(dphidxp):
                    dphidx = dphidxm
                else:
                    dphidx = dphidxp
                #dphidx = min(dphidx1,dphidxm,dphidxp)

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

                E += c*math.sin(theta[x,y,z])**4 + a*(dthetadx**2 + dthetady**2 + dthetadz**2 + (dphidx**2 + dphidy**2 + dphidz**2)*math.sin(theta[x,y,z])**2) +\
                    (b*(-2*dphidz*math.sin(theta[x,y,z])**2 + math.cos(phi[x,y,z])*(-2*dthetady + dphidx*math.sin(2*theta[x,y,z])) +\
                    (2*dthetadx + dphidy*math.sin(2*theta[x,y,z]))*math.sin(phi[x,y,z])))/2.

            #print(E)
    print(E)
    return(E)
