"""Cutoff-scale induced G for one Dirac fermion with a physical cutoff (lattice in x, sharp physical transverse cutoff).
Usage: cutoff_G.py Lambda_perp L widths"""
import numpy as np,sys,json,time
from cutoff_common import *
Lam=float(sys.argv[1]); L=int(sys.argv[2]); W=[float(w) for w in sys.argv[3].split(',')]
x0=L/2; eps=np.array([0,0.03,-0.03,0.06,-0.06,0.09,-0.09]); out={}; t=time.time()
E0=E_vac(lambda x:0*x,L,Lam,eps.max())
for w in W:
    g=lambda x:np.exp(-(x-x0)**2/(2*w*w))
    Es=[E0 if e==0 else E_vac(lambda x,e=e:e*g(x),L,Lam,eps.max()) for e in eps]
    P=np.polyfit(eps,np.array(Es)-E0,6); a1=P[-2]; a2=P[-3]
    P4=np.polyfit(eps,np.array(Es)-E0,4); 
    xx=np.linspace(0,L,400001); dx=xx[1]-xx[0]; gg=g(xx); gp=np.gradient(gg,dx)
    Vol1=2*np.sum(gg)*dx; Vol2=2*np.sum(gg**2)*dx; I2=2*np.sum(gp**2)*dx
    out[w]=dict(a1=a1,a2=a2,a1_d4=P4[-2],a2_d4=P4[-3],Vol1=Vol1,Vol2=Vol2,I2=I2)
    Lam1=a1/Vol1; C=(a2-Lam1*Vol2)/I2
    print(f"Lam_perp={Lam:.4f} L={L} w={w}: Lambda_vac={Lam1:.8f}  C(Lambda2=Lambda1)={C:.6f}  (deg4 fit: {(P4[-3]-P4[-2]/Vol1*Vol2)/I2:.6f})  {time.time()-t:.0f}s",flush=True)
json.dump({str(k):v for k,v in out.items()},open(f"cutG_L{Lam:.3f}_N{L}.json","w"),indent=1)
