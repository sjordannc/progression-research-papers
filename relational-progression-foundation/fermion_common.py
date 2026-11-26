"""Doubler-free radial lattice Dirac fermion (staggered u on sites j, v on half-sites j+1/2) with the covariant radial
derivative B = e^{-Phi} d/dr e^{Phi}, Phi' = kappa/f.  Geometry dl^2 = dr^2 + f(r)^2 dOmega^2, lapse N(r): H_N = sqrt(N) H sqrt(N)."""
import numpy as np
from scipy.linalg import eigvalsh_tridiagonal
import os
JCUT=float(os.environ.get('JCUT','5'))
PV_c=np.array([1.,-3.,3.,-1.]); PV_M2=np.array([0.,1.,2.,3.])
_gx,_gw=np.polynomial.legendre.leggauss(6)
def int_inv_f(fn,a,b):
    mid=(a+b)/2; half=(b-a)/2
    return half*sum(wt/fn(mid+half*x) for x,wt in zip(_gx,_gw))
def sites(kappa,Nr,xmin=0.0):
    """positions (step 1/2) and types (0=u,1=v). Termination chosen so the regular solution is kept and the
    lattice-only zero modes (u~r^-kappa at the origin for kappa>0, v~r^-|kappa| for kappa<0) are excluded."""
    if kappa<0: x=np.arange(1,Nr+0.75,0.5)          # u1 v1.5 ... uNr vNr+1/2
    else:       x=np.arange(0.5,Nr+0.25,0.5)        # v0.5 u1 ... uNr
    t=((2*x)%2==1).astype(int)                       # half-integer -> v
    keep=x>=xmin
    x=x[keep]; t=t[keep]
    want=0 if kappa<0 else 1                          # first site type
    if t[0]!=want: x=x[1:]; t=t[1:]
    return x,t
class Geo:
    def __init__(s,fn,Nfn,Nr): s.fn=fn; s.Nfn=Nfn; s.Nr=Nr; s.cache={}
def chain(geo,kappa,M,xmin=None):
    if xmin is None: xmin=max(0.0,abs(kappa)/JCUT-1)
    key=(kappa,xmin)
    if key not in geo.cache:
        x,t=sites(kappa,geo.Nr,xmin); I=int_inv_f(geo.fn,x[:-1],x[1:]); N=geo.Nfn(x)
        e=np.where(t[:-1]==0,-np.exp(-kappa*I),np.exp(kappa*I))*np.sqrt(N[:-1]*N[1:])
        geo.cache[key]=(x,t,N,e)
    x,t,N,e=geo.cache[key]
    d=np.where(t==0,M*N,-M*N)
    return d,e
def E0_PV(fn,Nfn,Nr,m,Kmax=3600):
    geo=Geo(fn,Nfn,Nr); per=np.zeros(Kmax)
    for K in range(1,Kmax+1):
        if K/JCUT-1>=Nr-2: break
        tot=0.0
        for kappa in (K,-K):
            for c,M2 in zip(PV_c,PV_M2*m*m):
                d,e=chain(geo,kappa,np.sqrt(M2))
                tot+=c*(-0.5)*np.sum(np.abs(eigvalsh_tridiagonal(d,e)))
        per[K-1]=2*K*tot
    return per
def H_dense(geo,kappa,M,xmin=None):
    d,e=chain(geo,kappa,M,xmin); x=geo.cache[(kappa,xmin if xmin is not None else max(0.0,abs(kappa)/JCUT-1))][0]
    return np.diag(d)+np.diag(e,1)+np.diag(e,-1),x
