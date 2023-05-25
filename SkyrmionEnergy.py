import math
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import minimize
from numba import jit

import InitialiseSkyrmion

# if not os.path.exists('results'):
#     os.makedirs('results')
#
# os.chdir('results')
@jit
def sphere(ks,kb,kt,b,c, phi,theta,Lx,Ly,Lz,dx,dy,dz,R):
    E = 0.0
    centrex = int(Lx/2)
    centrey = int(Ly/2)
    for z in range(0,Lz):

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
                dthetadzsq = (theta[x,y,za]+theta[x,y,zb] - 2.*theta[x,y,z])/(dz*dz)
                dthetadx = (theta[xr,y,z]-theta[xl,y,z])/(2.*dx)
                dthetady = (theta[x,ya,z]-theta[x,yb,z])/(2.*dy)
                dthetadxsq = (theta[xr,y,z]+theta[xl,y,z] - 2.*theta[x,y,z])/(dx*dx)
                dthetadysq = (theta[x,ya,z]+theta[x,yb,z] - 2.*theta[x,y,z])/(dy*dy)
                dphidz = (phi[x,y,za]-phi[x,y,zb])/(2.*dz)
                dphidzsq = (phi[x,y,za]+phi[x,y,zb] - 2.*phi[x,y,z])/(dz**2)

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

                    E += (ks*(np.sin(theta[x,y,z])*(dthetadz - dphidy*np.cos(phi[x,y,z]) + dphidx*np.sin(phi[x,y,z])) - \
                        np.cos(theta[x,y,z])*(dthetadx*np.cos(phi[x,y,z]) + dthetady*np.sin(phi[x,y,z])))**2) + \
                        (kb*(np.sin(theta[x,y,z])**2*(dthetadz*np.cos(theta[x,y,z]) + np.sin(theta[x,y,z])*(dthetadx*np.cos(phi[x,y,z]) + dthetady*np.sin(phi[x,y,z])))**2 + \
                        (dthetadz*np.cos(theta[x,y,z])**2*np.sin(phi[x,y,z]) + np.cos(phi[x,y,z])*np.sin(theta[x,y,z])**2*(dphidx*np.cos(phi[x,y,z]) + dphidy*np.sin(phi[x,y,z])) + \
                        (np.sin(2*theta[x,y,z])*(dthetady*np.sin(phi[x,y,z])**2 + np.cos(phi[x,y,z])*(dphidz + dthetadx*np.sin(phi[x,y,z]))))/2.)**2 + \
                        (dthetadz*np.cos(theta[x,y,z])**2*np.cos(phi[x,y,z]) - np.sin(theta[x,y,z])**2*np.sin(phi[x,y,z])*(dphidx*np.cos(phi[x,y,z]) + dphidy*np.sin(phi[x,y,z])) + \
                        (np.sin(2*theta[x,y,z])*(dthetadx + dthetadx*np.cos(2*phi[x,y,z]) - 2*dphidz*np.sin(phi[x,y,z]) + dthetady*np.sin(2*phi[x,y,z])))/4.)**2)) + \
                        (b*(-2*dphidz*np.sin(theta[x,y,z])**2 + np.cos(phi[x,y,z])*(-2*dthetady + dphidx*np.sin(2*theta[x,y,z])) + (2*dthetadx + dphidy*np.sin(2*theta[x,y,z]))*np.sin(phi[x,y,z])))/2. + \
                        (kt*(-2*dphidz*np.sin(theta[x,y,z])**2 + np.cos(phi[x,y,z])*(-2*dthetady + dphidx*np.sin(2*theta[x,y,z])) + (2*dthetadx + dphidy*np.sin(2*theta[x,y,z]))*np.sin(phi[x,y,z]))**2)/4.

                #dphidy = min(abs(dphidy1),abs(dphidym),abs(dphidyp))
                #One constant energy simplification.
                # E += c*math.sin(theta[x,y,z])**4 + a*(dthetadx**2 + dthetady**2 + dthetadz**2 + (dphidx**2 + dphidy**2 + dphidz**2)*math.sin(theta[x,y,z])**2) +\
                #     (b*(-2*dphidz*math.sin(theta[x,y,z])**2 + math.cos(phi[x,y,z])*(-2*dthetady + dphidx*math.sin(2*theta[x,y,z])) +\
                #     (2*dthetadx + dphidy*math.sin(2*theta[x,y,z]))*math.sin(phi[x,y,z])))/2.

                # Energy with all constants included.
  #               E += (dthetadz**2*(kb + ks + (kb - ks)*math.cos(2*theta[x,y,z])))/2. - b*dthetady*math.cos(phi[x,y,z]) + dthetady**2*kt*math.cos(phi[x,y,z])**2 + \
  # dthetadx**2*ks*math.cos(theta[x,y,z])**2*math.cos(phi[x,y,z])**2 + (dphidz**2*(kb + kt + (kb - kt)*math.cos(2*theta[x,y,z]))*math.sin(theta[x,y,z])**2)/2. + \
  # dphidy**2*ks*math.cos(phi[x,y,z])**2*math.sin(theta[x,y,z])**2 + dthetadx**2*kb*math.cos(theta[x,y,z])**2*math.cos(phi[x,y,z])**4*math.sin(theta[x,y,z])**2 + \
  # (dphidx**2*(kb + 2*ks + kt + 2*(-kb + kt)*math.cos(2*theta[x,y,z])*math.cos(phi[x,y,z])**2 + (kb - 2*ks + kt)*math.cos(2*phi[x,y,z]))*math.sin(theta[x,y,z])**2)/4. + c*math.sin(theta[x,y,z])**4 +\
  # dthetadx**2*kb*math.cos(phi[x,y,z])**2*math.sin(theta[x,y,z])**4 - a*dphidy*dthetadx*math.cos(phi[x,y,z])**2*math.sin(2*theta[x,y,z]) + dphidy*dthetadx*ks*math.cos(phi[x,y,z])**2*math.sin(2*theta[x,y,z]) +\
  # b*dthetadx*math.sin(phi[x,y,z]) + 2*dthetadx*dthetady*kb*math.cos(theta[x,y,z])**2*math.cos(phi[x,y,z])**3*math.sin(theta[x,y,z])**2*math.sin(phi[x,y,z]) + (b*dphidy*math.sin(2*theta[x,y,z])*math.sin(phi[x,y,z]))/2. +\
  # dthetadx**2*kt*math.sin(phi[x,y,z])**2 + dthetady**2*ks*math.cos(theta[x,y,z])**2*math.sin(phi[x,y,z])**2 + dthetady**2*kb*math.sin(theta[x,y,z])**4*math.sin(phi[x,y,z])**2 -\
  # a*dphidy*dthetadx*math.sin(2*theta[x,y,z])*math.sin(phi[x,y,z])**2 + dphidy*dthetadx*kt*math.sin(2*theta[x,y,z])*math.sin(phi[x,y,z])**2 + (dphidy**2*kt*math.sin(2*theta[x,y,z])**2*math.sin(phi[x,y,z])**2)/4. +\
  # (dthetadx*dthetady*kb*math.cos(phi[x,y,z])*math.sin(2*theta[x,y,z])**2*math.sin(phi[x,y,z])**3)/2. + dphidy**2*kb*math.sin(theta[x,y,z])**4*math.sin(phi[x,y,z])**4 +\
  # (dthetady**2*kb*math.sin(2*theta[x,y,z])**2*math.sin(phi[x,y,z])**4)/4. + dphidz*math.sin(theta[x,y,z])**2*\
  #  (-b + 2*dthetady*(-a + kt)*math.cos(phi[x,y,z]) + 2*dthetadx*(a - kt)*math.sin(phi[x,y,z]) + (kb - kt)*math.sin(2*theta[x,y,z])*(dphidx*math.cos(phi[x,y,z]) + dphidy*math.sin(phi[x,y,z]))) +\
  # 2*dthetadz*math.sin(theta[x,y,z])*((a - ks)*math.sin(theta[x,y,z])*(dphidy*math.cos(phi[x,y,z]) - dphidx*math.sin(phi[x,y,z])) +\
  #    (kb - ks)*math.cos(theta[x,y,z])*(dthetadx*math.cos(phi[x,y,z]) + dthetady*math.sin(phi[x,y,z]))) - dthetadx*dthetady*kt*math.sin(2*phi[x,y,z]) +\
  # dthetadx*dthetady*ks*math.cos(theta[x,y,z])**2*math.sin(2*phi[x,y,z]) + dthetadx*dthetady*kb*math.sin(theta[x,y,z])**4*math.sin(2*phi[x,y,z]) +\
  # (dphidy*dthetady*ks*math.sin(2*theta[x,y,z])*math.sin(2*phi[x,y,z]))/2. - (dphidy*dthetady*kt*math.sin(2*theta[x,y,z])*math.sin(2*phi[x,y,z]))/2. +\
  # (dphidy**2*kb*math.sin(theta[x,y,z])**4*math.sin(2*phi[x,y,z])**2)/4. + (dthetadx**2*kb*math.sin(2*theta[x,y,z])**2*math.sin(2*phi[x,y,z])**2)/16. +\
  # (dthetady**2*kb*math.sin(2*theta[x,y,z])**2*math.sin(2*phi[x,y,z])**2)/16. + (dphidx*\
  #    (dthetady*(2*a - ks - kt + (ks - kt)*math.cos(2*phi[x,y,z]))*math.sin(2*theta[x,y,z]) + math.cos(phi[x,y,z])*math.sin(2*theta[x,y,z])*(b + 2*dthetadx*(-ks + kt)*math.sin(phi[x,y,z])) +\
  #      dphidy*(kb - 2*ks + kt + (-kb + kt)*math.cos(2*theta[x,y,z]))*math.sin(theta[x,y,z])**2*math.sin(2*phi[x,y,z])))/2.


# Old/Legacy Energy with no saddle splay.
  #               E += c*math.sin(theta[x,y,z])**4 + ks*(-(dthetadz*math.sin(theta[x,y,z])) + math.cos(phi[x,y,z])*(dthetadx*math.cos(theta[x,y,z]) + dphidy*math.sin(theta[x,y,z])) + \
  #     (dthetady*math.cos(theta[x,y,z]) - dphidx*math.sin(theta[x,y,z]))*math.sin(phi[x,y,z]))**2 + \
  # (b*(-2*dphidz*math.sin(theta[x,y,z])**2 + math.cos(phi[x,y,z])*(-2*dthetady + dphidx*math.sin(2*theta[x,y,z])) + (2*dthetadx + dphidy*math.sin(2*theta[x,y,z]))*math.sin(phi[x,y,z])))/2. + \
  # (kt*(-2*dphidz*math.sin(theta[x,y,z])**2 + math.cos(phi[x,y,z])*(-2*dthetady + dphidx*math.sin(2*theta[x,y,z])) + (2*dthetadx + dphidy*math.sin(2*theta[x,y,z]))*math.sin(phi[x,y,z]))**2)/4. + \
  # kb*(dthetadz**2*math.cos(theta[x,y,z])**2 + dthetadz*math.sin(2*theta[x,y,z])*(dthetadx*math.cos(phi[x,y,z]) + dthetady*math.sin(phi[x,y,z])) + \
  #    math.sin(theta[x,y,z])**2*((dthetadx*math.cos(phi[x,y,z]) + dthetady*math.sin(phi[x,y,z]))**2 + (dphidz*math.cos(theta[x,y,z]) + math.sin(theta[x,y,z])*(dphidx*math.cos(phi[x,y,z]) + dphidy*math.sin(phi[x,y,z])))**2))


    #print(E)
    return(E)
