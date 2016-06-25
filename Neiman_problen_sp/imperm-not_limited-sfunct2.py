#!/usr/bin/python
import numpy as np
from sys import path
path.append("./../")
from smodule.sfunct2 import *

Kx=1
Py=1
Phi1NL=-0.25/sp.pi*sp.log((x0-y0)**2+(x1-y1)**2)
omgNL=[ Py*sp.diff(Phi1NL,sy[j]) for j in range(2) ]
Phi2NL=-0.5/sp.pi*sp.atan((x1-y1)/(x0-y0))
V2NL=[sp.ratsimp(Kx*sp.diff(Phi2NL,sx[j])) for j in range(2)]
VSNL=[sp.ratsimp(Kx*sp.diff(Phi1NL,sx[j])) for j in range(2)]

nI=800
qSrc=np.array([np.pi])
zSrc=np.array([[0.0],[2.0]],dtype='d')
minDist=0.02
hv=0.2 # for velocity field grid

y_crcl=0.0
x_crcl=0.0
R=1.0
thet=np.linspace(2*np.pi,0,nI+1,endpoint=True)
LI=np.array([x_crcl+R*np.cos(thet),y_crcl+R*np.sin(thet)])

xIc=(LI[:,:nI]+LI[:,1:])/2.
Nxc=Normal(LI)
M=MVDL(x=xIc,b=LI,sV2=V2NL)
A=M[0]*Nxc[0]+M[1]*Nxc[1]
WS=VSources(q=qSrc,x=xIc,z=zSrc,sVS=VSNL)
f=-WS[0]*Nxc[0]-WS[1]*Nxc[1]
#print 'cond(A)', np.linalg.cond(A)
#print 'det(A)', np.linalg.det(A)
A=np.r_[np.c_[np.ones((nI,1)),A],np.c_[0,np.ones((1,nI))]]
#print 'cond(A)', np.linalg.cond(A)
#print 'det(A)', np.linalg.det(A)
f=np.r_[f,0]
g=np.linalg.solve(A,f)
gI=g[1:]
#print 'gamma=',g[0]
xs=np.mgrid[-2.5:2.5:hv,-2.:3:hv]
xs=xs[:,((xs[0]-zSrc[0][0])**2+(xs[1]-zSrc[1][0])**2)>0.2**2]
xs=xs[:,np.logical_or((xs[0]**2+xs[1]**2)>(1.+minDist)**2,(xs[0]**2+xs[1]**2)<(1.-minDist)**2)]
Ws=VSources(qSrc,xs,zSrc,sVS=VSNL)+VDL(x=xs,b=LI,g=gI,sV2=V2NL)

# Testing on the exact solution
qOut=np.array([np.pi,np.pi,-np.pi])
zOut=np.array([[0.0,0.0,0.0],[2.0,0.5,0.0]],dtype='d')
xOut=xs[:,(xs[0]**2+xs[1]**2)>(1+minDist)**2]
WOut=VSources(qOut,xOut,zOut,VSNL)
WNOut=VSources(qSrc,xOut,zSrc,VSNL)+VDL(xOut,LI,g[1:],V2NL)
etaOut=(1.-np.sqrt(WNOut[0]**2+WNOut[1]**2)/np.sqrt(WOut[0]**2+WOut[1]**2))*100
etaOutMax=np.max(etaOut)


print "nI = ", nI
print "etaOutMax = ",etaOutMax

np.savetxt('bI-not_limited.dat',np.c_[LI[0],LI[1]],fmt='%6.4f',delimiter=' ')
np.savetxt('WI-not_limited.dat',np.c_[xs[0],xs[1],Ws[0],Ws[1]],fmt='%6.4f',delimiter=' ')

from pylab import quiver,show, gca,plot,text, axis, grid
plot(LI[0],LI[1],color='black',lw=2)
quiver(xs[0],xs[1],Ws[0],Ws[1],color='r')
quiver(xOut[0],xOut[1],WOut[0],WOut[1],color='g')
text(-2.5,4.1,ur"$\eta_{Out,max}=%6.4f$per." % etaOutMax)
grid(True)
axis('equal')
show()
