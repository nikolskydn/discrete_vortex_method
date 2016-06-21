#!/usr/bin/python

# agreements:
#     border has the format plus point
#     global symbol's variables:
#         x=(x0,x1) --- point collocation  ( row )
#         y=(y0,y1) --- point with the singularity ( column )
#     z --- point z contains a source (sink)
#
# a,b,c,d --- vector; A,B,C,D --- matrix

import numpy as np
import sympy as sp
myeps=0.0000000000001
myeps2=myeps**2

x0,x1,y0,y1=sp.symbols('x0 x1 y0 y1',real=True)
_Phi1_=-0.25/sp.pi*sp.log((x0-y0)**2+(x1-y1)**2)
_Psi2_=-_Phi1_
_Kx_=1
_Hx_=1
_Px_=1
_Ky_=1
_Hy_=1
_Py_=1

def withoutZero(sFunct=None,minDist=myeps):
    def vectorFunct(A,B,C,D):
        gd=((A-C)**2+(B-D)**2)>minDist**2
        AA=A.copy(); AA[np.logical_not(gd)]+=1
        result=np.zeros_like(A)
        result[gd]=sp.lambdify((x0,x1,y0,y1),sFunct,modules='numpy')(AA,B,C,D)[gd]
        return result
    return vectorFunct

def lmbd(sFunct=None):
    def vectorFunct(A,B,C,D):
        return sp.lambdify((x0,x1,y0,y1),sFunct,modules='numpy')(A,B,C,D)
    return vectorFunct

def varphiSources(q,x,z,Phi1=_Phi1_):
    Z0,X0=np.meshgrid(z[0].ravel(),x[0].ravel())
    Z1,X1=np.meshgrid(z[1].ravel(),x[1].ravel())
    return -np.sum(q*lmbd(Phi1)(X0,X1,Z0,Z1),axis=1)

def VSources(q,x,z,Kx=_Kx_,Phi1=_Phi1_):
    Z0,X0=np.meshgrid(z[0].ravel(),x[0].ravel())
    Z1,X1=np.meshgrid(z[1].ravel(),x[1].ravel())
    V0=lmbd(Kx*sp.diff(Phi1,x0))(X0,X1,Z0,Z1)
    V1=lmbd(Kx*sp.diff(Phi1,x1))(X0,X1,Z0,Z1)
    return -np.array([np.sum(q*V0,axis=1), np.sum(q*V1,axis=1)]) 

def Normal(x):
    n=x[0].shape[0]-1
    return np.array([-x[1,1:]+x[1,:n],x[0,1:]-x[0,:n]])

def normal(x):
    nr=Normal(x)
    return nr/np.sqrt(nr[0]**2+nr[1]**2)

def ieFredgolm2dk(xrw,xcl,lmd,alph,q,z,vphFun,Py=_Py_,Phi1=_Phi1_,isDiag=1,K=1,minDist=myeps):
    nrw=xrw[0].shape[0]-1; ncl=xcl[0].shape[0]-1
    x=np.array([(xrw[0,:nrw]+xrw[0,1:])/2,(xrw[1,:nrw]+xrw[1,1:])/2])
    y=np.array([(xcl[0,:ncl]+xcl[0,1:])/2,(xcl[1,:ncl]+xcl[1,1:])/2])
    YC0,XC0=np.meshgrid(y[0,:],x[0,:])
    YC1,XC1=np.meshgrid(y[1,:],x[1,:])
    Nrml=Normal(xcl)
    KF0=withoutZero(Py*sp.diff(Phi1,y0),minDist)(XC0,XC1,YC0,YC1)
    KF1=withoutZero(Py*sp.diff(Phi1,y1),minDist)(XC0,XC1,YC0,YC1)
    MtrA=KF0*Nrml[0]+KF1*Nrml[1]
    MtrA=-2*lmd*MtrA + ( np.eye(nrw) if isDiag else 0 )
    vtrf=2.*lmd*vphFun(q,x,z,Phi1)-2*K*alph*x[1]
    return MtrA,vtrf

def Theta(X0,X1,Y0,Y1,rEps=myeps):
    R2=(X0-Y0)**2+(X1-Y1)**2
    THETA=np.ones_like(R2)
    rEps2=rEps**2
    R2=R2/rEps2; R=R2**.5; R3=R2*R; R5=R3*R2; R7=R5*R2; R9=R7*R2
    THETA[R2<rEps]=.125*(63.*R5-90.*R7+35.*R9)[R2<rEps]
    return THETA 

def V2(X0,X1,Y0,Y1,Hx=_Hx_,Psi2=_Psi2_,rEps=myeps):
    return Theta(X0,X1,Y0,Y1,rEps)*np.array([withoutZero(sp.diff(Psi2,x1)/Hx)(X0,X1,Y0,Y1),withoutZero(-sp.diff(Psi2,x0)/Hx)(X0,X1,Y0,Y1)])

def MatrixVDL(xrw,xcl,Hx=_Hx_,Psi2=_Psi2_,rEps=myeps):
    ncl=xcl[0].shape[0]-1
    Y0,X0=np.meshgrid(xcl[0,:ncl],xrw[0].ravel())
    Y1,X1=np.meshgrid(xcl[1,:ncl],xrw[1].ravel())
    YNEXT0,X0=np.meshgrid(xcl[0,1:],xrw[0].ravel())
    YNEXT1,X1=np.meshgrid(xcl[1,1:],xrw[1].ravel())
    return V2(X0,X1,Y0,Y1,Hx,Psi2,rEps)-V2(X0,X1,YNEXT0,YNEXT1,Hx,Psi2,rEps)

def VDL(xrw,xcl,gg,Hx=_Hx_,Psi2=_Psi2_,rEps=myeps):
    W=MatrixVDL(xrw,xcl,Hx,Psi2,rEps)*gg
    return np.array([np.sum(W[0],axis=1),np.sum(W[1],axis=1)])


def ieSingular(xrw,xcl,q,z,WFun,Kx=_Kx_,Hx=_Hx_,Phi1=_Phi1_,Psi2=_Psi2_): 
    nrw=xrw[0].shape[0]-1
    ncl=xcl[0].shape[0]-1
    xc=np.array([.5*(xrw[0,:nrw]+xrw[0,1:]),.5*(xrw[1,:nrw]+xrw[1,1:])])
    W=MatrixVDL(xc,xcl,Hx,Psi2,myeps)
    nrml=normal(xrw)
    A=W[0]*nrml[0].reshape(-1,1)+W[1]*nrml[1].reshape(-1,1)
    W0=WFun(q,xc,z,Kx,Phi1)
    f=-W0[0]*nrml[0]-W0[1]*nrml[1]
    return A,f

