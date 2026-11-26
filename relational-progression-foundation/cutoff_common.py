import numpy as np
from scipy.linalg import eigvalsh_tridiagonal
def chain_E(k,sig_u,sig_v,xu,Lam):
    """sum|E| of the staggered 1D Dirac chain(s) for coordinate transverse momentum k; sites kept where k e^{-sigma}<Lam.
    order u1 v1.5 u2 ... ; hopping -1 (u->v) and +1 (v->u) between neighbouring kept sites only."""
    n=len(xu); mu=np.empty(2*n); mu[0::2]=k*np.exp(-sig_u); mu[1::2]=k*np.exp(-sig_v)
    keep=mu<Lam
    d=np.where(np.arange(2*n)%2==0,mu,-mu)
    e=np.where(np.arange(2*n-1)%2==0,-1.0,1.0)
    e=e*(keep[:-1]&keep[1:])
    d=d[keep]; idx=np.where(keep)[0]
    # off-diagonals between consecutive kept sites (zero if not adjacent)
    ee=np.zeros(len(idx)-1)
    adj=(idx[1:]-idx[:-1])==1
    ee[adj]=e[idx[:-1][adj]]
    if len(d)==0: return 0.0
    if len(d)==1: return abs(d[0])
    return np.sum(np.abs(eigvalsh_tridiagonal(d,ee)))
def gl(a,b,n):
    x,w=np.polynomial.legendre.leggauss(n); return (b-a)/2*x+(a+b)/2,(b-a)/2*w
def E_vac(sig_fn,L,Lam,eps_max,nbulk=6,nseg=24,npiece=6):
    xu=np.arange(1,L+1,dtype=float); xv=xu+0.5
    su=sig_fn(xu); sv=sig_fn(xv)
    K0=Lam*np.exp(-eps_max)*0.999          # below every breakpoint, for every config: identical nodes
    # bulk [0,K0]: composite Gauss, refined near 0
    edges=np.concatenate([[0],K0*np.geomspace(1e-4,1,nseg)])
    tot=0.0
    for a,b in zip(edges[:-1],edges[1:]):
        ks,ws=gl(a,b,nbulk)
        for k,w in zip(ks,ws): tot+=w*k*chain_E(k,su,sv,xu,Lam)
    bp=np.unique(np.round(Lam*np.exp(np.concatenate([su,sv])),14))
    bp=np.concatenate([[K0],bp[bp>K0]])
    for a,b in zip(bp[:-1],bp[1:]):
        if b-a<1e-15: continue
        ks,ws=gl(a,b,npiece)
        for k,w in zip(ks,ws): tot+=w*k*chain_E(k,su,sv,xu,Lam)
    return -tot/(2*np.pi)          # two 2-component chains per k: -2*(1/2)sum|E| * (1/2pi) int k dk
def S_plane(L,Lam,nbulk=8,nseg=30):
    """eta: entanglement per unit area of a plane cut (between v_{L/2+1/2} and u_{L/2+1}) in the flat chain"""
    xu=np.arange(1,L+1,dtype=float); n=L; cut=L  # sites 0..2n-1; region = first L sites (u1..v_{L/2+1/2})
    def S_k(k):
        d=np.where(np.arange(2*n)%2==0,k,-k).astype(float); e=np.where(np.arange(2*n-1)%2==0,-1.0,1.0)
        H=np.diag(d)+np.diag(e,1)+np.diag(e,-1); w,v=np.linalg.eigh(H); Vm=v[:,w<0]
        CA=Vm[:cut]@Vm[:cut].T; c=np.linalg.eigvalsh(CA); c=c[(c>1e-15)&(c<1-1e-15)]
        return float(np.sum(-c*np.log(c)-(1-c)*np.log(1-c)))
    edges=np.concatenate([[0],Lam*np.geomspace(1e-4,1,nseg)]); tot=0.0
    for a,b in zip(edges[:-1],edges[1:]):
        ks,ws=gl(a,b,nbulk)
        for k,w in zip(ks,ws): tot+=w*k*2*S_k(k)
    return tot/(2*np.pi)
