"""Two-clock test (amended pre-registration): insert the linearised solution into the lattice stationarity equations."""
import numpy as np,sys,json,time
from scipy.integrate import cumulative_trapezoid as ctz
from clock_common import *
use_ct=int(sys.argv[1]) if len(sys.argv)>1 else 1
R0=float(sys.argv[2]) if len(sys.argv)>2 else 100.0; W=float(sys.argv[3]) if len(sys.argv)>3 else 12.0
etaD=4.6404e-4; G=1/(4*etaD); C=-1/(16*np.pi*G); M=1e-3*R0/G
rr=np.linspace(1e-6,400,400001); dr=rr[1]-rr[0]
mu=M*np.exp(-(rr-R0)**2/(2*W*W))/(np.sqrt(2*np.pi)*W); m=ctz(mu,rr,initial=0)
g=G*m/rr**2
nu=-(np.trapezoid(g,rr)-ctz(g,rr,initial=0))-G*M/rr[-1]   # Newtonian lapse, nu(infinity)=0
I1=ctz(np.where(rr>0,m/rr,0),rr,initial=0); s=-G*I1/rr
nuF=Tab(rr,nu); sF=Tab(rr,s); zero=Tab(rr,0*rr)
nuM=Tab(rr,-nu); sM=Tab(rr,-s)
T,sup=test_splines(0.0,220.0,10.0); nT=len(T)
J=np.array([np.sum(mu*t(rr))*dr for t in T])
Lam1=1.560709e-5; dL2=-2.47e-4*Lam1 if use_ct else 0.0
def ct(nu_v,s_v):
    N=np.exp(nu_v); f=rr*(1+s_v)
    gn=np.array([4*np.pi*Lam1*np.sum(N*f*f*t(rr))*dr for t in T])
    gs=np.array([4*np.pi*np.sum((Lam1*N*2*f*rr+dL2*2*rr*rr*s_v)*t(rr))*dr for t in T])
    return gn,gs
t0=time.time(); out={}
def sym(nuP,sP,nuMm,sMm,nuv,sv):
    _,an,as_=grad_fn(nuP,sP,T,T); _,bn,bs=grad_fn(nuMm,sMm,T,T)
    cpn,cps=ct(nuv,sv); cmn,cms=ct(-nuv,-sv)
    return (an-bn)/2-(cpn-cmn)/2,(as_-bs)/2-(cps-cms)/2
Hn_n,Hn_s=sym(nuF,zero,nuM,zero,nu,0*s); print('nu response done',f'{time.time()-t0:.0f}s',flush=True)
Hs_n,Hs_s=sym(zero,sF,zero,sM,0*nu,s); print('s response done',f'{time.time()-t0:.0f}s',flush=True)
# rows: [nu rows (energy constraint) ; s rows (lapse equation)]
A=np.vstack([np.concatenate([Hn_n,Hn_s]),np.concatenate([Hs_n,Hs_s])]).T; rhs=-np.concatenate([J,np.zeros(nT)])
# GR expectation for comparison: residual of exact solution should vanish; scale of each row group
lam=np.linalg.lstsq(A,rhs,rcond=None)[0]; res=A@lam-rhs
scaleN=np.linalg.norm(J); 
r_unit=A@np.ones(2)-rhs
inside=[k for k,(a,b) in enumerate(sup) if b<=80]
print(f"lam_nu = {lam[0]:.4f}  lam_s = {lam[1]:.4f}")
print(f"residual of exact-G solution (lam=1): |R|/|J| = {np.linalg.norm(r_unit)/scaleN:.4f}  (constraint rows {np.linalg.norm(r_unit[:nT])/scaleN:.4f}, lapse rows {np.linalg.norm(r_unit[nT:])/scaleN:.4f})")
print(f"residual after fit: |R|/|J| = {np.linalg.norm(res)/scaleN:.4f};  inside-shell rows: {np.linalg.norm(np.concatenate([res[:nT][inside],res[nT:][inside]]))/scaleN:.4f}")
pts=dict(A1=30.,A2=60.,B1=160.,B2=210.)
Phi={k:float(nuF(np.array([r]))[0]) for k,r in pts.items()}
for a,b in [('B1','B2'),('A2','B2'),('A1','A2')]:
    print(f"clock ratio {a}/{b}: GR 1{Phi[a]-Phi[b]:+.4e}  lattice 1{lam[0]*(Phi[a]-Phi[b]):+.4e}")
out=dict(lam=lam.tolist(),res_unit=float(np.linalg.norm(r_unit)/scaleN),res_fit=float(np.linalg.norm(res)/scaleN),
 res_inside=float(np.linalg.norm(np.concatenate([res[:nT][inside],res[nT:][inside]]))/scaleN),Phi=Phi,G=G,M=M,use_ct=use_ct,
 A=A.tolist(),rhs=rhs.tolist(),sup=sup,s_pts={k:float(sF(np.array([r]))[0]) for k,r in pts.items()})
json.dump(out,open(f"clock_test_ct{use_ct}_R{R0:g}_W{W:g}.json","w"),indent=1)
