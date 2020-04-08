import os
import sys
from math import pi
from numpy import *
from matplotlib.pyplot import *
import pyPLUTO as pp

plutodir = os.environ['PLUTO_DIR']
wdir = './'
print(wdir)
nlinf = pp.nlast_info(w_dir=wdir)
print(nlinf)
for timestep in range(1,2):
    D = pp.pload (timestep, w_dir=wdir) # Loading the data into a pload object D.
    x = D.x1
    y = D.x2
    z = D.x3
    vx = D.vx1
    vy = D.vx2
    vz = D.vx3
    rho = D.rho

    nx = size(x)
    ny = size(y)
    nz = size(z)

    r = zeros((nx,ny,nz))

    for ix in range (nx):
        for iy in range (ny):
            for iz in range (nz):
                r[ix,iy,iz] = sqrt(x[ix]*x[ix] + y[iy]*y[iy] + z[iz]*z[iz])

                vr = zeros((nx,ny,nz))

    for ix in range (nx):
        for iy in range (ny):
            for iz in range (nz):
                vr[ix,iy,iz] = (x[ix]*vx[ix,iy,iz] + y[iy]*vy[ix,iy,iz] + z[iz]*vz[ix,iy,iz]) / r[ix,iy,iz]

    l = 2.49597871e+13
    rho_0 = 1.64e-24
    v_0 = 1e5
    mdot = 4*pi*r*r*l*l*rho*rho_0*vr*v_0
    print(mdot)

    for ix in range (nx):
        for iy in range (ny):
            for iz in range (nz):
                if timestep == 1:
                    plot(r[ix,iy,iz] ,mdot[ix,iy,iz],'ro')
                else:
                    plot(r[ix,iy,iz] ,mdot[ix,iy,iz],'k+')
                    
xlabel('radius')
ylabel('accretion rate')
title('bondi problem accrection rate 1 black hole')
savefig('bondi1_009.png')
