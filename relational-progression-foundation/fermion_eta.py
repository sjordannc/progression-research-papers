"""Entanglement side for the lattice Dirac fermion: PV-combined entanglement entropy of a ball (flat).
Usage: fermion_eta.py m"""
import numpy as np,sys,json,time
from fermion_common import *
m=float(sys.argv[1]); ns=[20,30,40,50]; Kmax=1500; one=lambda r:np.ones_like(r)
def S_cut(H,x,Rc):
    w,v=np.linalg.eigh(H); Vm=v[:,w<0]; sel=x<Rc
    if not sel.any(): return 0.0
    CA=Vm[sel]@Vm[sel].T
    c=np.linalg.eigvalsh(CA); c=c[(c>1e-15)&(c<1-1e-15)]
    return float(np.sum(-c*np.log(c)-(1-c)*np.log(1-c)))
out={}; t=time.time()
for n in ns:
    Nr=n+110; geo=Geo(lambda r:r,one,Nr); per=np.zeros(Kmax)
    for K in range(1,Kmax+1):
        if K/JCUT-1>=n: break
        for kappa in (K,-K):
            for c,M2 in zip(PV_c,PV_M2*m*m):
                H,x=H_dense(geo,kappa,np.sqrt(M2)); per[K-1]+=c*2*K*S_cut(H,x,n+0.75)
    out[n]=dict(S=[per[:L].sum() for L in (100,200,300,500,1500)])
    print(n,out[n],f"{time.time()-t:.0f}s",flush=True)
json.dump(out,open(f"feta_m{m}.json","w"),indent=1)
