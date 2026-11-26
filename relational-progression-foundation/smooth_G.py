"""Cutoff-scale induced G, smooth physical cutoff. Usage: smooth_G.py Lambda L widths qcut_factor"""
import numpy as np,sys,json,time
from smooth_common import *
Lam=float(sys.argv[1]); L=int(sys.argv[2]); W=[float(w) for w in sys.argv[3].split(',')]; qf=float(sys.argv[4]) if len(sys.argv)>4 else 3.7
qc=qf*Lam; x0=L/2; eps=np.array([0,0.03,-0.03,0.06,-0.06,0.09,-0.09]); out={}; t=time.time()
E0=E_vac(lambda x:0*x,L,Lam,qc,eps.max())[0]
for w in W:
    g=lambda x:np.exp(-(x-x0)**2/(2*w*w))
    Es=np.array([E0 if e==0 else E_vac(lambda x,e=e:e*g(x),L,Lam,qc,eps.max())[0] for e in eps])-E0
    P6=np.polyfit(eps,Es,6); P4=np.polyfit(eps,Es,4)
    xx=np.linspace(0,L,400001); dx=xx[1]-xx[0]; gg=g(xx); gp=np.gradient(gg,dx); gpp=np.gradient(gp,dx)
    v=dict(a1=P6[-2],a2=P6[-3],a1_4=P4[-2],a2_4=P4[-3],Vol1=2*np.sum(gg)*dx,Vol2=2*np.sum(gg**2)*dx,I2=2*np.sum(gp**2)*dx,I4=np.sum(gpp**2)*dx)
    out[w]=v; Lam1=v['a1']/v['Vol1']
    print(f"Lam={Lam} L={L} qf={qf} w={w}: Lambda_vac={Lam1:.10f}  C(L2=L1)={(v['a2']-Lam1*v['Vol2'])/v['I2']:+.6f}  deg4: {(v['a2_4']-v['a1_4']/v['Vol1']*v['Vol2'])/v['I2']:+.6f}  {time.time()-t:.0f}s",flush=True)
json.dump({str(k):v for k,v in out.items()},open(f"sG_L{Lam}_N{L}_q{qf}.json","w"),indent=1)
