"""Mechanism check (R3): a second clock built from a different mechanism -- the bound level of a Dirac mass well on the
radial lattice (kappa=1) -- placed at several radii in the solved background (lattice-fitted lapse and ruler),
compared with a two-level clock whose rate is N(r) Delta E. Reports both relative to the far reference."""
import numpy as np,json,sys
from scipy.integrate import cumulative_trapezoid as ctz
from scipy.linalg import eigh_tridiagonal
from fermion_common import sites
d=json.load(open(sys.argv[1])); lam_nu,lam_s=d['lam']; G=d['G']; M=d['M']; R0,W=100.0,float(sys.argv[2])
rr=np.linspace(1e-6,400,400001); mu=M*np.exp(-(rr-R0)**2/(2*W*W))/(np.sqrt(2*np.pi)*W); m=ctz(mu,rr,initial=0); g=G*m/rr**2
nu=lam_nu*(-(np.trapezoid(g,rr)-ctz(g,rr,initial=0))-G*M/rr[-1]); s=lam_s*(-G*ctz(np.where(rr>0,m/rr,0),rr,initial=0)/rr)
amp=float(sys.argv[3]) if len(sys.argv)>3 else 1.0   # optional amplification of the field (tests tidal corrections)
nuF=lambda x:np.interp(x,rr,amp*nu); sF=lambda x:np.interp(x,rr,amp*s)
_gx,_gw=np.polynomial.legendre.leggauss(6)
def level(rc,curved,kappa=1,Nr=300,Mout=2.0,Min=0.3,hw=4.0):
    x,t=sites(kappa,Nr,0.0); N=np.exp(nuF(x)) if curved else np.ones_like(x)
    xa,xb=x[:-1],x[1:]; mid=(xa+xb)/2; half=(xb-xa)/2; q=mid[:,None]+half[:,None]*_gx; qw=half[:,None]*_gw
    I=np.sum(qw/(q*(1+(sF(q) if curved else 0))),axis=1)
    tau=np.where(t[:-1]==0,-1.0,1.0); e=tau*np.exp(tau*kappa*I)*np.sqrt(N[:-1]*N[1:])
    Mx=np.where(np.abs(x-rc)<hw,Min,Mout); dd=np.where(t==0,1.0,-1.0)*Mx*N
    w=eigh_tridiagonal(dd,e,eigvals_only=True,select='v',select_range=(1e-9,Mout*0.999*(N.max())))
    return w.min()
pts=dict(A1=30.,A2=60.,B1=160.,B2=210.,REF=280.)
two={k:float(np.exp(nuF(np.array([r]))[0])) for k,r in pts.items()}
lev={k:level(r,True)/level(r,False) for k,r in pts.items()}
print("radius  two-level rate/flat   bound-level rate/flat   difference")
for k,r in pts.items(): print(f"{k} r={r:5.0f}  {two[k]-1:+.6e}   {lev[k]-1:+.6e}   {lev[k]-two[k]:+.2e}")
for a,b in [('B1','B2'),('A2','B2'),('A1','A2')]:
    rt=two[a]/two[b]-1; rl=lev[a]/lev[b]-1
    # reference independence: same ratio computed via REF as intermediate reference
    rref=(two[a]/two['REF'])/(two[b]/two['REF'])-1
    print(f"{a}/{b}: two-level {rt:+.6e}  bound-level {rl:+.6e}  (mechanism difference {abs(rl-rt)/abs(rt):.2%} of the shift)  via reference clock: {rref:+.6e}")
