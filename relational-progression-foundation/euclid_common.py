"""Euclidean-time lattice version (A_T = time spacing; 0 = continuous time): single-particle energies E -> (2/a_t) arcsinh(a_t E/2)."""
"""Planar geometry dl^2 = dx^2 + e^{2 sigma(x)}(dy^2+dz^2); one Dirac fermion reduced to 1D staggered chains per transverse k.
Smooth physical transverse cutoff: dispersion mu(q) = q exp(q^2/Lam^2), q = k e^{-sigma(x)}; sites removed where q > qcut.
The 4-component 1D problem is H = H_chain(q) (x) tau_z + M (1 (x) tau_x), so H^2 = H_chain^2 + M^2 exactly and
|E| = sqrt(eps_n^2 + M^2) with eps_n the eigenvalues of the massless chain (this is how Pauli-Villars masses enter)."""
import numpy as np
from scipy.linalg import eigvalsh_tridiagonal
import os
POW=float(os.environ.get('REG_POW','2'))
AT=float(os.environ.get('A_T','1'))
def feps(E):
    return E if AT==0 else (2/AT)*np.arcsinh(AT*E/2)
def mu_of(q,Lam):
    # REG_POW=0: transverse lattice dispersion 2 sin(q/2) (spacing 1, same as x and t), cut at the zone edge q=pi (disk approximation)
    return 2*np.sin(q/2) if POW==0 else q*np.exp((q/Lam)**POW)
def chain(k,su,sv,Lam,qcut):
    n=len(su); q=np.empty(2*n); q[0::2]=k*np.exp(-su); q[1::2]=k*np.exp(-sv)
    keep=q<qcut; mu=mu_of(q,Lam); d=np.where(np.arange(2*n)%2==0,mu,-mu); e=np.where(np.arange(2*n-1)%2==0,-1.0,1.0)
    idx=np.where(keep)[0]
    if len(idx)==0: return []
    dd=d[idx]; ee=np.zeros(len(idx)-1); adj=(idx[1:]-idx[:-1])==1; ee[adj]=e[idx[:-1][adj]]
    segs=[]; cuts=np.where(ee==0)[0]; starts=np.concatenate([[0],cuts+1]); ends=np.concatenate([cuts+1,[len(dd)]])
    for a,b in zip(starts,ends):
        ev=np.array([dd[a]]) if b-a==1 else eigvalsh_tridiagonal(dd[a:b],ee[a:b-1]); segs.append((ev,np.abs(dd[a:b])))
    return segs
def chain_E(k,su,sv,Lam,qcut,Ms):
    """returns, for each mass M in Ms: 2*sum sqrt(eps^2+M^2) - 2*sum_sites sqrt(mu^2+M^2)  (= sum|E| of the 4-comp system minus exact volume term)"""
    out=np.zeros(len(Ms))
    for ev,dabs in chain(k,su,sv,Lam,qcut):
        for i,M in enumerate(Ms):
            out[i]+=2*(np.sum(feps(np.sqrt(ev**2+M*M)))-np.sum(feps(np.sqrt(dabs**2+M*M))))
    return out
def gl(a,b,n):
    x,w=np.polynomial.legendre.leggauss(n); return (b-a)/2*x+(a+b)/2,(b-a)/2*w
def E_vac(sig_fn,L,Lam,qcut,eps_max,Ms=(0.0,),nbulk=8,nseg=40,npiece=6):
    """vacuum energy per unit coordinate transverse area, -(1/2) sum|E| integrated with d^2k/(2pi)^2, for each mass"""
    xu=np.arange(1,L+1,dtype=float); xv=xu+0.5; su=sig_fn(xu); sv=sig_fn(xv)
    K0=qcut*np.exp(-eps_max)*0.999
    edges=np.concatenate([[0],K0*np.geomspace(1e-4,1,nseg)]); tot=np.zeros(len(Ms))
    for a,b in zip(edges[:-1],edges[1:]):
        ks,ws=gl(a,b,nbulk)
        for k,w in zip(ks,ws): tot+=w*k*chain_E(k,su,sv,Lam,qcut,Ms)
    bp=np.unique(np.round(qcut*np.exp(np.concatenate([su,sv])),13)); bp=np.concatenate([[K0],bp[bp>K0]])
    for a,b in zip(bp[:-1],bp[1:]):
        if b-a<1e-14: continue
        ks,ws=gl(a,b,npiece)
        for k,w in zip(ks,ws): tot+=w*k*chain_E(k,su,sv,Lam,qcut,Ms)
    return -0.5*tot/(2*np.pi)
def S_plane(L,Lam,qcut,M=0.0,nbulk=8,nseg=40):
    """entanglement per unit area of the plane cut between v_{L+1/2} and u_{L+1} in a flat chain of 2L cells (4-comp)"""
    n=2*L; tz=np.diag([1.,-1.]); tx=np.array([[0,1.],[1,0]])
    def S_k(k):
        m=mu_of(k,Lam); d=np.where(np.arange(2*n)%2==0,m,-m); e=np.where(np.arange(2*n-1)%2==0,-1.0,1.0)
        Hc=np.diag(d)+np.diag(e,1)+np.diag(e,-1)
        if M==0.0:
            w,v=np.linalg.eigh(Hc); Vm=v[:,w<0]; CA=Vm[:2*L]@Vm[:2*L].T; c=np.linalg.eigvalsh(CA); f=2
        else:
            H=np.kron(Hc,tz)+M*np.kron(np.eye(2*n),tx); w,v=np.linalg.eigh(H); Vm=v[:,w<0]; CA=Vm[:4*L]@Vm[:4*L].T; c=np.linalg.eigvalsh(CA); f=1
        c=c[(c>1e-15)&(c<1-1e-15)]; return f*float(np.sum(-c*np.log(c)-(1-c)*np.log(1-c)))
    edges=np.concatenate([[0],qcut*np.geomspace(1e-4,1,nseg)]); tot=0.0
    for a,b in zip(edges[:-1],edges[1:]):
        ks,ws=gl(a,b,nbulk)
        for k,w in zip(ks,ws): tot+=w*k*S_k(k)
    return tot/(2*np.pi)
