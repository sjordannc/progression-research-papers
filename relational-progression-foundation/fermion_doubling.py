"""P0: doubling check. Vacuum-energy density Lambda (first-order response to a volume bump) for
(a) the staggered covariant fermion, (b) the naive unstaggered fermion (u,v on the same sites, central differences),
compared with the real scalar (note 5): covariant prediction Lambda_D/Lambda_S = -4; one doubler pair -> -8."""
import numpy as np,sys
from fermion_common import *
from scipy.linalg import eigvalsh
m=0.2; Nr=360; r0=150.; w=24.; eps=0.02; one=lambda r:np.ones_like(r)
g=lambda r:np.exp(-(r-r0)**2/(2*w*w)); rr=np.linspace(0,Nr,200001); dr=rr[1]-rr[0]
Vol1=8*np.pi*np.sum(rr**2*g(rr))*dr
Kmax=int(sys.argv[1]) if len(sys.argv)>1 else 3600
a=[]
for s in (1,-1):
    a.append(E0_PV(lambda r,s=s:r*(1+s*eps*g(r)),one,Nr,m,Kmax).sum())
LamD=(a[0]-a[1])/(2*eps)/Vol1
LamS=-4.02503e-06   # scalar, note 5 (m=0.2)
print(f"staggered covariant fermion: Lambda_D={LamD:.6e}  Lambda_D/Lambda_S={LamD/LamS:.4f}  (prediction -4)",flush=True)
# naive: u_j, v_j on the same sites; d/dr -> central difference; kappa/f on site
def naive_E(fn):
    j=np.arange(1,Nr+1,dtype=float); f=fn(j); tot=0.0
    for K in range(1,Kmax+1):
        j0=max(0,K//JCUT-1)
        if j0>=Nr-2: break
        for kappa in (K,-K):
            n=Nr-j0; D=(np.eye(n,k=1)-np.eye(n,k=-1))/2; B=D+np.diag(kappa/f[j0:])
            for c,M2 in zip(PV_c,PV_M2*m*m):
                ev=eigvalsh(B.T@B+M2*np.eye(n)); tot+=c*2*K*(-0.5)*2*np.sum(np.sqrt(np.abs(ev)))
    return tot
if len(sys.argv)>2:
    b=[naive_E(lambda r,s=s:r*(1+s*eps*g(r))) for s in (1,-1)]
    LamN=(b[0]-b[1])/(2*eps)/Vol1
    print(f"naive fermion: Lambda_N={LamN:.6e}  Lambda_N/Lambda_S={LamN/LamS:.4f}  (one doubler pair -> -8)",flush=True)
