#!/usr/bin/python
import numpy as np
import sympy as sp

myeps=0.000001
x0,x1,y0,y1=sp.symbols('x0 x1 y0 y1',real=True)
sx=[x0,x1]
sy=[y0,y1]
_Kx_=1
_Py_=1
_Phi1_=-0.25/sp.pi*sp.log((x0-y0)**2+(x1-y1)**2)
_sOmega_=[_Py_*sp.diff(_Phi1_,sy[j]) for j in range(2)]
_Phi2_=-0.5/sp.pi*sp.atan((x1-y1)/(x0-y0))
_sV2_=[sp.ratsimp(_Kx_*sp.diff(_Phi2_,sx[j])) for j in range(2)]
_sVS_=[_Kx_*sp.diff(_Phi1_,sx[j]) for j in range(2)]

def Normal(x):
    n=x[0].shape[0]-1
    return np.array([-x[1,1:]+x[1,:n],x[0,1:]-x[0,:n]])

def normal(x):
    nr=Normal(x)
    return nr/np.sqrt(nr[0]**2+nr[1]**2)

def withoutZero(sFunct):
    def wrapped(sExpress,M,minDist=myeps):
        gd=np.sum([(M[i,1]-M[i,0])**2 for i in range(2)],axis=0)>minDist**2
        AA=M[0,0].copy()
        AA[np.logical_not(gd)]+=1
        result=np.zeros_like(M[0,0])
        result[gd]=sFunct(sExpress,np.array([[AA,M[0,1]],[M[1,0],M[1,1]]]),minDist)[gd]
        return result
    wrapped.orig=sFunct
    return wrapped

@withoutZero
def sF(sExpress,M,minDist=myeps):
    return sp.lambdify((y0,x0,y1,x1),sExpress,modules='numpy')(M[0,0],M[0,1],M[1,0],M[1,1])

def MPDL(x,y,Ny,sOmega=_sOmega_,minDist=myeps):
    M=np.array([np.meshgrid(y[j],x[j]) for j in range(2)])
    return np.sum([sF(sOmega[j],M,minDist)*Ny[j] for j in range(2)],axis=0)

def PDL(x,y,g,Ny,sOmega=_sOmega_,minDist=myeps):
    return np.sum(MPDL(x,y,Ny,sOmega,minDist)*g,axis=1)

def Theta(M,rEps=myeps):
    R2=np.sum([(M[i,1]-M[i,0])**2 for i in range(2)],axis=0)
    THETA=np.ones_like(R2)
    rEps2=rEps**2
    R2=R2/rEps2; R=R2**.5; R3=R2*R; R5=R3*R2; R7=R5*R2; R9=R7*R2
    THETA[R2<rEps]=.125*(63.*R5-90.*R7+35.*R9)[R2<rEps]
    return THETA

def V2(x,y,sV2=_sV2_,rEps=myeps):
    M=np.array([np.meshgrid(y[j],x[j]) for j in range(2)])
    return Theta(M,rEps)*np.array([sF(sV2[j],M) for j in range(2)])

def MVDL(x,b,sV2=_sV2_,rEps=myeps):
    return -V2(x,b[:,1:],sV2,rEps)+V2(x,b[:,:b.shape[1]-1],sV2,rEps)

def VDL(x,b,g,sV2=_sV2_,rEps=myeps):
    return np.sum(MVDL(x,b,sV2,rEps)*g,axis=2)

def varphiSources(q,x,z,sPhi1=_Phi1_):
    M=np.array([np.meshgrid(z[j],x[j]) for j in range(2)])
    return -np.sum(q*sF.orig(sPhi1,M),axis=1)

def VSources(q,x,z,sVS=_sVS_):
    M=np.array([np.meshgrid(z[j],x[j]) for j in range(2)])
    return -np.array([np.sum(q*sF.orig(sVS[j],M),axis=1) for j in range(2)])

#if __name__ == "__main__":
