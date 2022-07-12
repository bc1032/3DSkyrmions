import math
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import minimize

import InitialiseSkyrmion
import SkyrmionEnergy
Lx, Ly,Lz = 101, 101, 21
R = 30.0
numtwists = 1.0
twist = numtwists
dx,dy = 1,1
numits = 50000

tolerance = -1e-6
gamma = 1e-3
a,b,c = 10,-1,1#b is negative.a

if not os.path.exists('results'):
    os.makedirs('results')
if not os.path.exists('results/c%f' % c):
    os.makedirs('results/c%f' % c)
os.chdir('results/c%f' % c)

def derivativesphere(gamma,phi,theta):
    i=0
    Energy = []
    #phi = np.loadtxt('phi.dat')
    #theta = np.loadtxt('theta.dat')
    E = SkyrmionEnergy.sphere(a,b,c, phi,theta,Lx,Ly,Lz)
    phi,theta = minimise(gamma,phi,theta,a,b,c,i)
    #SkyrmionEnergy.sphere(a,b,c, phi,theta)
    n=0
    Energy.append(E)
    for i in range(0,numits):
        print(i)
        Eprev = E
        E = SkyrmionEnergy.sphere(a,b,c, phi, theta,Lx,Ly,Lz)
        Energy.append(E)
        print(-(Eprev - E))
        if (E - Eprev) > tolerance:
            np.savetxt('phifinal',phi)
            np.savetxt('thetafinal',theta)
            return(phi,theta)
        else:
            phi,theta = minimise(gamma,phi,theta,a,b,c,i)
            phiprev, thetaprev = phi,theta

    Energy = np.asarray(Energy)
    np.savetxt('Energyevolution.dat',Energy)
    np.savetxt('phifinal',phi)
    np.savetxt('thetafinal',theta)
    return(phi,theta)

def minimise(gamma,phi,theta,a,b,c,i):
    dphi = np.zeros((Lx,Ly,Lz))
    dtheta = np.zeros((Lx,Ly,Lz))
    #gamma = -1e-3
    for z in range(0,Lz):

        if z == 0:
                zb = Lz-1
                za = 1
        elif z == Lz-1:
            zb = Lz-2
            za = 0
        else:
            xl = x - 1
            xr = x + 1

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

                #dphidxdy = min(abs(dphidxdy1),abs(dphidxdym),abs(dphidxdyp))

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

                #dphidysq = min(abs(dphidysq1),abs(dphidysqm),abs(dphidysqp))

                dphi[x,y,z] = gamma*(2*math.sin(theta[x,y,z])*(-2*a*(dphidx*dthetadx + dphidy*dthetady)*math.cos(theta[x,y,z]) + (b - 2*a*dphidz)*dthetadz*math.cos(theta[x,y,z]) - \
                                a*(dphidxsq + dphidysq + dphidzsq)*math.sin(theta[x,y,z]) + b*dthetadx*math.cos(phi[x,y,z])*math.sin(theta[x,y,z]) + \
                                b*dthetady*math.sin(theta[x,y,z])*math.sin(phi[x,y,z])))

                dtheta[x,y,z] = gamma*((-2.0)*a*(dthetadxsq + dthetadysq + dthetadzsq)\
                                + a*(dphidx**2 + dphidy**2)*math.sin(2*theta[x,y,z])\
                                - b*dphidz*math.sin(2*theta[x,y,z]) + a*dphidz**2*math.sin(2*theta[x,y,z]) + 2*math.sin(theta[x,y,z])**2*(c*math.sin(2*theta[x,y,z])\
                                - b*(dphidx*math.cos(phi[x,y,z]) + dphidy*math.sin(phi[x,y,z]))))

    phi -= dphi
    theta -= dtheta
    file = open('finaldirectorfield.dat', 'w')
    phifile = open('phi.dat', 'w')
    thetafile = open('theta.dat', 'w')

    for z in range(0,Lz):
        for x in range(0,Lx):
            for y in range(0,Ly):
                i = math.sin(theta[x,y,z])*math.cos(phi[x,y,z])
                j = math.sin(theta[x,y,z])*math.sin(phi[x,y,z])
                k = math.cos(theta[x,y,z])
                phifile.write("%f\n" % phi[x,y,z])
                thetafile.write("%f\n" % theta[x,y,z])
                file.write("%f  %f  %f\n" % (i,j,k))

    file.close()
    phifile.close()
    thetafile.close()
    # np.savetxt('dtheta.dat', dtheta)
    # np.savetxt('dphi.dat', dphi)
    # np.savetxt('theta.dat', theta)
    # np.savetxt('phi.dat', phi)

    return(phi,theta)
phi,theta = InitialiseSkyrmion.initialise(Lx,Ly,Lz,R,twist)

phi,theta = derivativesphere(gamma,phi,theta)
file = open('finaldirectorfield.dat', 'w')

for x in range(0,Lx):
    for y in range(0,Ly):
        i = math.sin(theta[x,y,z])*math.cos(phi[x,y,z])
        j = math.sin(theta[x,y,z])*math.sin(phi[x,y,z])
        k = math.cos(theta[x,y,z])
        file.write("%f  %f  %f\n" % (i,j,k))
file.close()
