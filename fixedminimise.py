import math
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import minimize
from numba import jit


import InitialiseSkyrmion
import InitialiseLoad
import SkyrmionEnergy

dz = 0.5
dx = 1.0
dy = dx

Lx, Ly,Lz = 75, 75, int(51/dz)
R = 19
numtwists = 3
twist = numtwists
numits = 10000

tolerance = 0.001
gamma = 1e-2
#q0 = math.pi/13.0
q0 = math.pi/(2*R/numtwists)
ks,kb,kt = 1.0,1.0,1.0
b,c = -2.0*kt*q0,0.0#b is negative.a

b = -0.3#b is negative.

q0 = -b/(2.0*kt)
a=1
#0 = one constant, 1 = Initialising defects, 2 = full with no saddlesplay, 3 = varydz, 4 = startfromnew initials
mode = 1

if mode == 0:
    if not os.path.exists('newresultsoneconstant'):
        os.makedirs('newresultsoneconstant')
    if not os.path.exists('newresultsoneconstant/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist)):
        os.makedirs('newresultsoneconstant/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist))
    if not os.path.exists('newresultsoneconstant/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f/finaldirectorfield' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist)):
        os.makedirs('newresultsoneconstant/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f/finaldirectorfield' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist))
    if not os.path.exists('newresultsoneconstant/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f/phi' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist)):
        os.makedirs('newresultsoneconstant/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f/phi' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist))
    if not os.path.exists('newresultsoneconstant/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f/theta' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist)):
        os.makedirs('newresultsoneconstant/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f/theta' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist))
    os.chdir('newresultsoneconstant/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist))

elif mode == 1:
    if not os.path.exists('newInitialise'):
        os.makedirs('newInitialise')
    if not os.path.exists('newInitialise/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%fR%f' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist,R)):
        os.makedirs('newInitialise/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%fR%f' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist,R))
    if not os.path.exists('newInitialise/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%fR%f/finaldirectorfield' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist,R)):
        os.makedirs('newInitialise/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%fR%f/finaldirectorfield' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist,R))
    if not os.path.exists('newInitialise/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%fR%f/phi' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist,R)):
        os.makedirs('newInitialise/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%fR%f/phi' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist,R))
    if not os.path.exists('newInitialise/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%fR%f/theta' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist,R)):
        os.makedirs('newInitialise/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%fR%f/theta' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist,R))
    os.chdir('newInitialise/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%fR%f' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist,R))

elif mode == 2:
    if not os.path.exists('newresultsconstantq0'):
        os.makedirs('newresultsconstantq0')
    if not os.path.exists('newresultsconstantq0/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist)):
        os.makedirs('newresultsconstantq0/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist))
    if not os.path.exists('newresultsconstantq0/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f/finaldirectorfield' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist)):
        os.makedirs('newresultsconstantq0/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f/finaldirectorfield' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist))
    if not os.path.exists('newresultsconstantq0/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f/phi' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist)):
        os.makedirs('newresultsconstantq0/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f/phi' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist))
    if not os.path.exists('newresultsconstantq0/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f/theta' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist)):
        os.makedirs('newresultsconstantq0/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f/theta' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist))
    os.chdir('newresultsconstantq0/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist))

elif mode == 3:
    if not os.path.exists('varydz'):
        os.makedirs('varydz')
    if not os.path.exists('varydz/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist)):
        os.makedirs('varydz/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist))
    if not os.path.exists('varydz/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f/finaldirectorfield' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist)):
        os.makedirs('varydz/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f/finaldirectorfield' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist))
    if not os.path.exists('varydz/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f/phi' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist)):
        os.makedirs('varydz/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f/phi' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist))
    if not os.path.exists('varydz/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f/theta' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist)):
        os.makedirs('varydz/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f/theta' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist))
    os.chdir('varydz/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%f' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist))

elif mode == 4:
    if not os.path.exists('startingfromnewintialised'):
        os.makedirs('startingfromnewintialised')
    if not os.path.exists('startingfromnewintialised/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%fR%f' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist,R)):
        os.makedirs('startingfromnewintialised/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%fR%f' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist,R))
    if not os.path.exists('startingfromnewintialised/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%fR%f/finaldirectorfield' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist,R)):
        os.makedirs('startingfromnewintialised/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%fR%f/finaldirectorfield' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist,R))
    if not os.path.exists('startingfromnewintialised/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%fR%f/phi' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist,R)):
        os.makedirs('startingfromnewintialised/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%fR%f/phi' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist,R))
    if not os.path.exists('startingfromnewintialised/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%fR%f/theta' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist,R)):
        os.makedirs('startingfromnewintialised/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%fR%f/theta' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist,R))
    os.chdir('startingfromnewintialised/ks%fkb%fkt%fb%fq0%fdx%fdy%fdz%ftwist%fR%f' % (ks,kb,kt,abs(b),q0,dx,dy,dz,twist,R))

#@jit
def derivativesphere(gamma,phi,theta):
    i=0
    #phi = np.loadtxt('phi.csv')
    #theta = np.loadtxt('theta.csv')
    #E = SkyrmionEnergy.sphere(ks,kb,kt,b,c, phi,theta,Lx,Ly,Lz,dx,dy,dz,R)
    #SkyrmionEnergy.sphere(a,b,c, phi,theta)
    n=0
    E = 0.0
    #phi,theta = minimise(gamma,phi,theta,ks,kb,kt,b,c,i)
    for i in range(0,numits+1):
        Eprev = E

        E = SkyrmionEnergy.sphere(ks,kb,kt,b,c, phi,theta,Lx,Ly,Lz,dx,dy,dz,R)

        if i % 25 == 0:
            print(i, E)
            print(-(Eprev - E))
        if abs((E - Eprev)) < tolerance:
            return(phi,theta)
        else:
            phi,theta = minimise(gamma,phi,theta,ks,kb,kt,b,c,i)
            phiprev, thetaprev = phi,theta
    return(phi,theta)


#@jit
def minimise(gamma,phi,theta,ks,kb,kt,b,c,i):
    dphi = np.zeros((Lx,Ly,Lz))
    dtheta = np.zeros((Lx,Ly,Lz))
    #gamma = -1e-3

    centrex = int(Lx/2)
    centrey = int(Ly/2)

    if i % 25 == 0:
        file = open('finaldirectorfield/finaldirectorfield%d.csv' % i, 'w')
        file.write('x,y,z,vx,vy,vz\n')
        phifile = open('phi/phi%d.csv' % i, 'w')
        thetafile = open('theta/theta%d.csv' % i, 'w')
        phifile.write('x,y,z,phi\n')
        thetafile.write('x,y,z,theta\n')

        for z in range(0,Lz):
            for x in range(0,Lx):
                for y in range(0,Ly):
                    i = math.sin(theta[x,y,z])*math.cos(phi[x,y,z])
                    j = math.sin(theta[x,y,z])*math.sin(phi[x,y,z])
                    k = math.cos(theta[x,y,z])
                    file.write("%d,%d,%d,%f,%f,%f\n" % (x,y,z,i,j,k))
                    phifile.write("%d,%d,%d,%f\n" % (x,y,z,phi[x,y,z]))
                    thetafile.write("%d,%d,%d,%f\n" % (x,y,z,theta[x,y,z]))
        file.close()
        phifile.close()
        thetafile.close()

    phi, theta = minimisecalculation(dphi,dtheta,centrex, centrey,gamma,phi,theta,ks,kb,kt,b,c,i)

    return(phi,theta)

@jit
def minimisecalculation(dphi,dtheta,centrex, centrey,gamma,phi,theta,ks,kb,kt,b,c,i):

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

                ddx = abs(centrex - x)
                ddy = abs(centrey - y)

                if y == 0:
                    ya = y+1
                    yb = Ly-1
                elif y == Ly-1:
                    ya = 0
                    yb = Ly-2
                else:
                    ya = y+1
                    yb = y-1

                dthetadz = (theta[x,y,za]-theta[x,y,zb])/(2.*dz)
                dthetadzsq = (theta[x,y,za]+theta[x,y,zb] - 2.*theta[x,y,z])/(dz**2)
                dthetadx = (theta[xr,y,z]-theta[xl,y,z])/(2.*dx)
                dthetady = (theta[x,ya,z]-theta[x,yb,z])/(2.*dy)
                dthetadxsq = (theta[xr,y,z]+theta[xl,y,z] - 2.*theta[x,y,z])/(dx**2)
                dthetadysq = (theta[x,ya,z]+theta[x,yb,z] - 2.*theta[x,y,z])/(dy**2)
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

                #dphidysq = min(abs(dphidysq1),abs(dphidysqm),abs(dphidysqp))
                #if math.sqrt(ddx**2 + ddy**2) <= R:
                if math.sqrt(ddx**2 + ddy**2) <= R*(twist-1)/twist:
                    #E one constant Lagrangian
                    if mode == 0:
                        dphi[x,y,z] =  gamma*(2*math.sin(theta[x,y,z])*(-2*a*(dphidx*dthetadx + dphidy*dthetady)*math.cos(theta[x,y,z]) + (b - 2*a*dphidz)*dthetadz*math.cos(theta[x,y,z]) - \
                                        a*(dphidxsq + dphidysq + dphidzsq)*math.sin(theta[x,y,z]) + b*dthetadx*math.cos(phi[x,y,z])*math.sin(theta[x,y,z]) + \
                                        b*dthetady*math.sin(theta[x,y,z])*math.sin(phi[x,y,z])))

                        dtheta[x,y,z] = gamma*((-2.0)*a*(dthetadxsq + dthetadysq + dthetadzsq)\
                                        + a*(dphidx**2 + dphidy**2)*math.sin(2*theta[x,y,z])\
                                        - b*dphidz*math.sin(2*theta[x,y,z]) + a*dphidz**2*math.sin(2*theta[x,y,z]) + 2*math.sin(theta[x,y,z])**2*(c*math.sin(2*theta[x,y,z])\
                                        - b*(dphidx*math.cos(phi[x,y,z]) + dphidy*math.sin(phi[x,y,z]))))

                    #E_full Lagrangian.
     #                elif mode == 1:
     #                    dphi[x,y,z] =  gamma*((math.sin(theta[x,y,z])*(-4*dphidx*dthetadx*kb*math.cos(theta[x,y,z]) - 4*dphidy*dthetady*kb*math.cos(theta[x,y,z]) - 8*dphidx*dthetadx*ks*math.cos(theta[x,y,z]) - \
     #       8*dphidy*dthetady*ks*math.cos(theta[x,y,z]) - 4*dphidx*dthetadx*kt*math.cos(theta[x,y,z]) - 4*dphidy*dthetady*kt*math.cos(theta[x,y,z]) + 4*dphidx*dthetadx*kb*math.cos(3*theta[x,y,z]) + \
     #       4*dphidy*dthetady*kb*math.cos(3*theta[x,y,z]) - 4*dphidx*dthetadx*kt*math.cos(3*theta[x,y,z]) - 4*dphidy*dthetady*kt*math.cos(3*theta[x,y,z]) - 3*dphidxsq*kb*math.sin(theta[x,y,z]) - \
     #       3*dphidysq*kb*math.sin(theta[x,y,z]) - 2*dphidzsq*kb*math.sin(theta[x,y,z]) - 4*dphidxsq*ks*math.sin(theta[x,y,z]) - 4*dphidysq*ks*math.sin(theta[x,y,z]) - dphidxsq*kt*math.sin(theta[x,y,z]) - \
     #       dphidysq*kt*math.sin(theta[x,y,z]) - 6*dphidzsq*kt*math.sin(theta[x,y,z]) + \
     #       4*math.cos(phi[x,y,z])*math.sin(theta[x,y,z])*(-(dphidz*(kb - kt)*(2*dthetadx*(1 + 2*math.cos(2*theta[x,y,z])) + dphidy*math.sin(2*theta[x,y,z]))) + \
     #          2*(b*dthetadx + dthetadydz*(ks - kt) + dphidxdz*(-kb + kt)*math.sin(2*theta[x,y,z]))) + dphidxsq*kb*math.sin(3*theta[x,y,z]) + dphidysq*kb*math.sin(3*theta[x,y,z]) - \
     #       2*dphidzsq*kb*math.sin(3*theta[x,y,z]) - dphidxsq*kt*math.sin(3*theta[x,y,z]) - dphidysq*kt*math.sin(3*theta[x,y,z]) + 2*dphidzsq*kt*math.sin(3*theta[x,y,z]) + \
     #       math.cos(2*phi[x,y,z])*(4*(-(dphidx*dthetadx*(kb - 2*ks + kt)) + 2*dthetadxdy*(-ks + kt))*math.cos(theta[x,y,z]) + 4*dphidx*dthetadx*(kb - kt)*math.cos(3*theta[x,y,z]) + \
     #          (-dphidxsq - 2*dphidx*dphidy + dphidysq)*(3*kb - 4*ks + kt)*math.sin(theta[x,y,z]) + \
     #          8*dthetady*(dphidy*math.cos(theta[x,y,z])*(kb - ks + (-kb + kt)*math.cos(2*theta[x,y,z])) + dthetadx*(kb + ks - 2*kt)*math.sin(theta[x,y,z])) - \
     #          (-dphidxsq - 2*dphidx*dphidy + dphidysq)*(kb - kt)*math.sin(3*theta[x,y,z])) + \
     #       4*math.sin(theta[x,y,z])*(2*dthetadxdz*(-ks + kt) + 2*dthetady*(b - dphidz*(kb - kt)*(1 + 2*math.cos(2*theta[x,y,z]))) - (2*dphidydz - dphidx*dphidz)*(kb - kt)*math.sin(2*theta[x,y,z]))*\
     #        math.sin(phi[x,y,z]) + 8*dthetadz*(b*math.cos(theta[x,y,z]) + 2*dphidz*math.cos(theta[x,y,z])*(-kt + (-kb + kt)*math.cos(2*theta[x,y,z])) - \
     #          (kb - kt)*math.sin(3*theta[x,y,z])*(dphidx*math.cos(phi[x,y,z]) + dphidy*math.sin(phi[x,y,z])) + (kb + ks - 2*kt)*math.cos(theta[x,y,z])*(dthetady*math.cos(phi[x,y,z]) - dthetadx*math.sin(phi[x,y,z])))\
     #        + (4*(-dthetadxsq + dthetadysq)*(-ks + kt)*math.cos(theta[x,y,z]) + 8*dphidy*dthetadx*math.cos(theta[x,y,z])*(-kb + ks + (kb - kt)*math.cos(2*theta[x,y,z])) + \
     #          8*dphidx*dthetady*math.cos(theta[x,y,z])*(-kb + ks + (kb - kt)*math.cos(2*theta[x,y,z])) + 4*dthetady**2*(kb + ks - 2*kt)*math.sin(theta[x,y,z]) + \
     #          (-4*dthetadx**2*(kb + ks - 2*kt) + (dphidx**2 - 2*dphidxdy)*(3*kb - 4*ks + kt))*math.sin(theta[x,y,z]) + \
     #          2*dphidy**2*(-kb + 2*ks - kt + (kb - kt)*math.cos(2*theta[x,y,z]))*math.sin(theta[x,y,z]) - (dphidx**2 - 2*dphidxdy)*(kb - kt)*math.sin(3*theta[x,y,z]))*math.sin(2*phi[x,y,z])))/4.)
     #
     #                    dtheta[x,y,z] = gamma*((-2*dthetadxsq*kb - 2*dthetadysq*kb - 4*dthetadzsq*kb - 2*dthetadxsq*ks - 2*dthetadysq*ks - 4*dthetadzsq*ks - 4*dthetadxsq*kt - 4*dthetadysq*kt + \
     # 2*dthetadxsq*kb*math.cos(2*theta[x,y,z]) + 2*dthetadysq*kb*math.cos(2*theta[x,y,z]) - 4*dthetadzsq*kb*math.cos(2*theta[x,y,z]) - 2*dthetadxsq*ks*math.cos(2*theta[x,y,z]) - \
     # 2*dthetadysq*ks*math.cos(2*theta[x,y,z]) + 4*dthetadzsq*ks*math.cos(2*theta[x,y,z]) + 16*c*math.cos(theta[x,y,z])*math.sin(theta[x,y,z])**3 - \
     # 8*math.cos(phi[x,y,z])*math.sin(theta[x,y,z])*(2*dthetadxdz*(kb - ks)*math.cos(theta[x,y,z]) + (b*dphidx + dphidydz*(-ks + kt))*math.sin(theta[x,y,z]) + \
     #    dphidz*(dthetady*(kb - ks)*math.cos(theta[x,y,z]) + dphidx*(-kb + ks + 2*(-kb + kt)*math.cos(2*theta[x,y,z]))*math.sin(theta[x,y,z]))) - 4*b*dphidz*math.sin(2*theta[x,y,z]) + \
     # 2*dphidx**2*kb*math.sin(2*theta[x,y,z]) + 2*dphidy**2*kb*math.sin(2*theta[x,y,z]) - 2*dthetadx**2*kb*math.sin(2*theta[x,y,z]) - 2*dthetady**2*kb*math.sin(2*theta[x,y,z]) + \
     # 4*dthetadz**2*(kb - ks)*math.sin(2*theta[x,y,z]) + 2*dphidx**2*ks*math.sin(2*theta[x,y,z]) + 2*dphidy**2*ks*math.sin(2*theta[x,y,z]) + 2*dthetadx**2*ks*math.sin(2*theta[x,y,z]) + \
     # 2*dthetady**2*ks*math.sin(2*theta[x,y,z]) + 4*dphidz**2*kt*math.sin(2*theta[x,y,z]) - dphidx**2*kb*math.sin(4*theta[x,y,z]) - dphidy**2*kb*math.sin(4*theta[x,y,z]) + \
     # 2*dphidz**2*kb*math.sin(4*theta[x,y,z]) + dphidx**2*kt*math.sin(4*theta[x,y,z]) + dphidy**2*kt*math.sin(4*theta[x,y,z]) - 2*dphidz**2*kt*math.sin(4*theta[x,y,z]) + \
     # math.cos(2*phi[x,y,z])*(2*(-2*dphidy*dthetadx - dthetadxsq - 2*dphidx*dthetady + dthetadysq)*(kb + ks - 2*kt) + \
     #    2*(2*dphidy*dthetadx + dthetadxsq + 2*dphidx*dthetady - dthetadysq)*(kb - ks)*math.cos(2*theta[x,y,z]) + \
     #    2*(dphidx**2*kb - dthetadx**2*kb + dthetady**2*(kb - ks) + dphidx**2*ks + dthetadx**2*ks - dphidy**2*(kb + ks - 2*kt) - 2*dphidx**2*kt + \
     #       2*dphidxdy*(-ks + kt))*math.sin(2*theta[x,y,z]) + (-dphidx**2 + dphidy**2)*(kb - kt)*math.sin(4*theta[x,y,z])) + \
     # 8*math.sin(theta[x,y,z])*(-((-(dphidz*dthetadx) + 2*dthetadydz)*(kb - ks)*math.cos(theta[x,y,z])) + dphidxdz*(-ks + kt)*math.sin(theta[x,y,z]) + \
     #    dphidy*(-b + dphidz*(kb - ks + 2*(kb - kt)*math.cos(2*theta[x,y,z])))*math.sin(theta[x,y,z]))*math.sin(phi[x,y,z]) - \
     # 4*dthetadz*(kb - ks)*(math.sin(2*theta[x,y,z])*(dphidy*math.cos(phi[x,y,z]) - dphidx*math.sin(phi[x,y,z])) + 2*math.cos(2*theta[x,y,z])*(dthetadx*math.cos(phi[x,y,z]) + dthetady*math.sin(phi[x,y,z]))) + \
     # 2*(-2*(-(dphidx*dthetadx) + dthetadxdy + dphidy*dthetady)*(kb + ks - 2*kt) + 2*(-(dphidx*dthetadx) + dthetadxdy + dphidy*dthetady)*(kb - ks)*math.cos(2*theta[x,y,z]) + \
     #    (2*dthetadx*dthetady*(-kb + ks) + 2*dphidx*dphidy*(kb + ks - 2*kt) + dphidxsq*(ks - kt) + dphidysq*(-ks + kt))*math.sin(2*theta[x,y,z]) + \
     #    dphidx*dphidy*(-kb + kt)*math.sin(4*theta[x,y,z]))*math.sin(2*phi[x,y,z]))/4.)

                    elif mode == 1 or mode == 2 or mode == 3 or mode == 4:
                            #E Lagrangian no saddle splay.
                        dphi[x,y,z] =  gamma*((math.sin(theta[x,y,z])*(-4*dphidx*dthetadx*kb*math.cos(theta[x,y,z]) - 4*dphidy*dthetady*kb*math.cos(theta[x,y,z]) - 8*dphidx*dthetadx*ks*math.cos(theta[x,y,z]) - \
       8*dphidy*dthetady*ks*math.cos(theta[x,y,z]) - 4*dphidx*dthetadx*kt*math.cos(theta[x,y,z]) - 4*dphidy*dthetady*kt*math.cos(theta[x,y,z]) + 4*dphidx*dthetadx*kb*math.cos(3*theta[x,y,z]) + \
       4*dphidy*dthetady*kb*math.cos(3*theta[x,y,z]) - 4*dphidx*dthetadx*kt*math.cos(3*theta[x,y,z]) - 4*dphidy*dthetady*kt*math.cos(3*theta[x,y,z]) - 3*dphidxsq*kb*math.sin(theta[x,y,z]) - \
       3*dphidysq*kb*math.sin(theta[x,y,z]) - 2*dphidzsq*kb*math.sin(theta[x,y,z]) - 4*dphidxsq*ks*math.sin(theta[x,y,z]) - 4*dphidysq*ks*math.sin(theta[x,y,z]) - dphidxsq*kt*math.sin(theta[x,y,z]) - \
       dphidysq*kt*math.sin(theta[x,y,z]) - 6*dphidzsq*kt*math.sin(theta[x,y,z]) + \
       4*math.cos(phi[x,y,z])*math.sin(theta[x,y,z])*(-(dphidz*(kb - kt)*(2*dthetadx*(1 + 2*math.cos(2*theta[x,y,z])) + dphidy*math.sin(2*theta[x,y,z]))) + \
          2*(b*dthetadx + dthetadydz*(ks - kt) + dphidxdz*(-kb + kt)*math.sin(2*theta[x,y,z]))) + dphidxsq*kb*math.sin(3*theta[x,y,z]) + dphidysq*kb*math.sin(3*theta[x,y,z]) - \
       2*dphidzsq*kb*math.sin(3*theta[x,y,z]) - dphidxsq*kt*math.sin(3*theta[x,y,z]) - dphidysq*kt*math.sin(3*theta[x,y,z]) + 2*dphidzsq*kt*math.sin(3*theta[x,y,z]) + \
       math.cos(2*phi[x,y,z])*(4*(-(dphidx*dthetadx*(kb - 2*ks + kt)) + 2*dthetadxdy*(-ks + kt))*math.cos(theta[x,y,z]) + 4*dphidx*dthetadx*(kb - kt)*math.cos(3*theta[x,y,z]) + \
          (-dphidxsq - 2*dphidx*dphidy + dphidysq)*(3*kb - 4*ks + kt)*math.sin(theta[x,y,z]) + \
          8*dthetady*(dphidy*math.cos(theta[x,y,z])*(kb - ks + (-kb + kt)*math.cos(2*theta[x,y,z])) + dthetadx*(kb + ks - 2*kt)*math.sin(theta[x,y,z])) - \
          (-dphidxsq - 2*dphidx*dphidy + dphidysq)*(kb - kt)*math.sin(3*theta[x,y,z])) + \
       4*math.sin(theta[x,y,z])*(2*dthetadxdz*(-ks + kt) + 2*dthetady*(b - dphidz*(kb - kt)*(1 + 2*math.cos(2*theta[x,y,z]))) - (2*dphidydz - dphidx*dphidz)*(kb - kt)*math.sin(2*theta[x,y,z]))*\
        math.sin(phi[x,y,z]) + 8*dthetadz*(b*math.cos(theta[x,y,z]) + 2*dphidz*math.cos(theta[x,y,z])*(-kt + (-kb + kt)*math.cos(2*theta[x,y,z])) - \
          (kb - kt)*math.sin(3*theta[x,y,z])*(dphidx*math.cos(phi[x,y,z]) + dphidy*math.sin(phi[x,y,z])) + (kb + ks - 2*kt)*math.cos(theta[x,y,z])*(dthetady*math.cos(phi[x,y,z]) - dthetadx*math.sin(phi[x,y,z])))\
        + (4*(-dthetadxsq + dthetadysq)*(-ks + kt)*math.cos(theta[x,y,z]) + 8*dphidy*dthetadx*math.cos(theta[x,y,z])*(-kb + ks + (kb - kt)*math.cos(2*theta[x,y,z])) + \
          8*dphidx*dthetady*math.cos(theta[x,y,z])*(-kb + ks + (kb - kt)*math.cos(2*theta[x,y,z])) + 4*dthetady**2*(kb + ks - 2*kt)*math.sin(theta[x,y,z]) + \
          (-4*dthetadx**2*(kb + ks - 2*kt) + (dphidx**2 - 2*dphidxdy)*(3*kb - 4*ks + kt))*math.sin(theta[x,y,z]) + \
          2*dphidy**2*(-kb + 2*ks - kt + (kb - kt)*math.cos(2*theta[x,y,z]))*math.sin(theta[x,y,z]) - (dphidx**2 - 2*dphidxdy)*(kb - kt)*math.sin(3*theta[x,y,z]))*math.sin(2*phi[x,y,z])))/4.)

                        dtheta[x,y,z] = gamma*(4*c*math.cos(theta[x,y,z])*math.sin(theta[x,y,z])**3 + (-2*dthetadxsq*kb - 2*dthetadysq*kb - 4*dthetadzsq*kb - 2*dthetadxsq*ks - 2*dthetadysq*ks - 4*dthetadzsq*ks - 4*dthetadxsq*kt -\
      4*dthetadysq*kt + 2*dthetadxsq*kb*math.cos(2*theta[x,y,z]) + 2*dthetadysq*kb*math.cos(2*theta[x,y,z]) - 4*dthetadzsq*kb*math.cos(2*theta[x,y,z]) - 2*dthetadxsq*ks*math.cos(2*theta[x,y,z]) - \
      2*dthetadysq*ks*math.cos(2*theta[x,y,z]) + 4*dthetadzsq*ks*math.cos(2*theta[x,y,z]) - \
      8*math.cos(phi[x,y,z])*math.sin(theta[x,y,z])*(2*dthetadxdz*(kb - ks)*math.cos(theta[x,y,z]) + (b*dphidx + dphidydz*(-ks + kt))*math.sin(theta[x,y,z]) + \
         dphidz*(dthetady*(kb - ks)*math.cos(theta[x,y,z]) + dphidx*(-kb + ks + 2*(-kb + kt)*math.cos(2*theta[x,y,z]))*math.sin(theta[x,y,z]))) - 4*b*dphidz*math.sin(2*theta[x,y,z]) + \
      2*dphidx**2*kb*math.sin(2*theta[x,y,z]) + 2*dphidy**2*kb*math.sin(2*theta[x,y,z]) - 2*dthetadx**2*kb*math.sin(2*theta[x,y,z]) - 2*dthetady**2*kb*math.sin(2*theta[x,y,z]) + \
      4*dthetadz**2*(kb - ks)*math.sin(2*theta[x,y,z]) + 2*dphidx**2*ks*math.sin(2*theta[x,y,z]) + 2*dphidy**2*ks*math.sin(2*theta[x,y,z]) + 2*dthetadx**2*ks*math.sin(2*theta[x,y,z]) + \
      2*dthetady**2*ks*math.sin(2*theta[x,y,z]) + 4*dphidz**2*kt*math.sin(2*theta[x,y,z]) - dphidx**2*kb*math.sin(4*theta[x,y,z]) - dphidy**2*kb*math.sin(4*theta[x,y,z]) + \
      2*dphidz**2*kb*math.sin(4*theta[x,y,z]) + dphidx**2*kt*math.sin(4*theta[x,y,z]) + dphidy**2*kt*math.sin(4*theta[x,y,z]) - 2*dphidz**2*kt*math.sin(4*theta[x,y,z]) + \
      math.cos(2*phi[x,y,z])*(2*(-2*dphidy*dthetadx - dthetadxsq - 2*dphidx*dthetady + dthetadysq)*(kb + ks - 2*kt) + \
         2*(2*dphidy*dthetadx + dthetadxsq + 2*dphidx*dthetady - dthetadysq)*(kb - ks)*math.cos(2*theta[x,y,z]) + \
         2*(dphidx**2*kb - dthetadx**2*kb + dthetady**2*(kb - ks) + dphidx**2*ks + dthetadx**2*ks - dphidy**2*(kb + ks - 2*kt) - 2*dphidx**2*kt + \
            2*dphidxdy*(-ks + kt))*math.sin(2*theta[x,y,z]) + (-dphidx**2 + dphidy**2)*(kb - kt)*math.sin(4*theta[x,y,z])) + \
      8*math.sin(theta[x,y,z])*(-((-(dphidz*dthetadx) + 2*dthetadydz)*(kb - ks)*math.cos(theta[x,y,z])) + dphidxdz*(-ks + kt)*math.sin(theta[x,y,z]) + \
         dphidy*(-b + dphidz*(kb - ks + 2*(kb - kt)*math.cos(2*theta[x,y,z])))*math.sin(theta[x,y,z]))*math.sin(phi[x,y,z]) - \
      4*dthetadz*(kb - ks)*(math.sin(2*theta[x,y,z])*(dphidy*math.cos(phi[x,y,z]) - dphidx*math.sin(phi[x,y,z])) + 2*math.cos(2*theta[x,y,z])*(dthetadx*math.cos(phi[x,y,z]) + dthetady*math.sin(phi[x,y,z]))) + \
      2*(-2*(-(dphidx*dthetadx) + dthetadxdy + dphidy*dthetady)*(kb + ks - 2*kt) + 2*(-(dphidx*dthetadx) + dthetadxdy + dphidy*dthetady)*(kb - ks)*math.cos(2*theta[x,y,z]) + \
         (2*dthetadx*dthetady*(-kb + ks) + 2*dphidx*dphidy*(kb + ks - 2*kt) + dphidxsq*(ks - kt) + dphidysq*(-ks + kt))*math.sin(2*theta[x,y,z]) + \
         dphidx*dphidy*(-kb + kt)*math.sin(4*theta[x,y,z]))*math.sin(2*phi[x,y,z]))/4.)

                    else:
                        dphi[x,y,z] = 0.0
                        dtheta[x,y,z] = 0.0
    phi -= dphi
    theta -= dtheta

    return(phi,theta)
if mode == 1:
    phi,theta = InitialiseSkyrmion.initialise(Lx,Ly,Lz, R, twist)
elif mode == 4:
    phi,theta = InitialiseLoad.loadnew(Lx,Ly,Lz)
else:
    phi,theta = InitialiseLoad.load(Lx,Ly,Lz)
phi,theta = derivativesphere(gamma,phi,theta)
#InitialiseSkyrmion.initialiseidealfinal(Lx,Ly,Lz, R, twist)
