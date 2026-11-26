"""Two-clock test: self-consistent static solve of the induced-gravity equations with the lattice fermion vacuum energy.
usage: clock_solve.py [iters] [ct: 1 = with the pre-registered k^0 artifact counterterm, 0 = without]"""
import numpy as np,sys,json,time
from clock_common import *
from gr_model import GR, rr, dr
iters=int(sys.argv[1]) if len(sys.argv)>1 else 3; use_ct=int(sys.argv[2]) if len(sys.argv)>2 else 1
b=Basis(); n=b.n
etaD=4.6404e-4; G=1/(4*etaD); C=-1/(16*np.pi*G)
R0,W=100.0,12.0; Phi_target=1e-3; M=Phi_target*R0/G
mu=lambda r: M*np.exp(-(r-R0)**2/(2*W*W))/(np.sqrt(2*np.pi)*W)
gr=GR(b,C,mu); H=gr.hess(); J=gr.Jsrc()
rows=[k for k,bb in enumerate(b.nu) if bb[2]<=220+1e-9]+[b.nn+k for k,bb in enumerate(b.s) if bb[2]<=220+1e-9]
Hr=H[rows]; print('unknowns',n,'equations',len(rows),'rank',np.linalg.matrix_rank(Hr),flush=True)
P=np.linalg.pinv(Hr)
u_gr=-P@J[rows]
Lam1=1.560709e-5; dL2=-2.47e-4*Lam1 if use_ct else 0.0
def ct_grad(a):
    """gradient of Lam1*4pi int N f^2 dr + dL2*4pi int r^2 s^2 dr (fine grid)"""
    nu,S,_,_=gr.fields(a); N=np.exp(nu); f=rr*(1+S)
    gn=[4*np.pi*Lam1*np.sum(N*f*f*B[0])*dr for B in gr.Bn]
    gs=[4*np.pi*np.sum((Lam1*N*2*f*rr+dL2*2*rr*rr*S)*B[0])*dr for B in gr.Bs]
    return np.array(gn+gs)
def R_lat(u):
    _,gp=grad_E0(u,b); _,gm=grad_E0(-u,b)
    return ((gp-gm)/2-(ct_grad(u)-ct_grad(-u))/2+J)[rows]
t=time.time(); u=u_gr.copy(); hist=[]
for it in range(iters):
    R=R_lat(u); rel=np.linalg.norm(R)/np.linalg.norm(J[rows]); hist.append(rel)
    print(f"iter {it}: |residual|/|source| = {rel:.3e}  {time.time()-t:.0f}s",flush=True)
    u=u-P@R
R=R_lat(u); rel=np.linalg.norm(R)/np.linalg.norm(J[rows]); hist.append(rel); print(f"final: |residual|/|source| = {rel:.3e}",flush=True)
pts=dict(A1=30.,A2=60.,B1=160.,B2=210.)
def at(a,r): return float(b.nu_of(a,np.array([r]))[0]),float(b.s_of(a,np.array([r]))[0])
out=dict(u_lat=u.tolist(),u_gr=u_gr.tolist(),hist=hist,G=G,M=M,use_ct=use_ct)
for k,r in pts.items(): out[k]=dict(lat=at(u,r),gr=at(u_gr,r))
json.dump(out,open(f"clock_solution_ct{use_ct}.json","w"),indent=1)
for k in pts: print(k,'nu lat %.6e gr %.6e | s lat %.6e gr %.6e'%(out[k]['lat'][0],out[k]['gr'][0],out[k]['lat'][1],out[k]['gr'][1]))
