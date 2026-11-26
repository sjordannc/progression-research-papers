"""Validation of the planar code against the covariant Pauli-Villars predictions C_D=-A m^2/(96 pi^2), eta_D=A m^2/(24 pi)."""
import numpy as np,sys,time
from euclid_common import *
m=float(sys.argv[4]) if len(sys.argv)>4 else 0.2; Lam=1.5; qc=(float(sys.argv[5]) if len(sys.argv)>5 else 3.7)*Lam; L=int(sys.argv[1]); W=[float(w) for w in sys.argv[2].split(',')]; mode=sys.argv[3]
c=np.array([1,-3,3,-1]); Ms=np.sqrt(np.array([0,1,2,3])*m*m); A=6*np.log(2)-3*np.log(3); t=time.time()
if mode=='eta':
    e=sum(ci*S_plane(L,Lam,qc,M=M) for ci,M in zip(c,Ms)); print(f"m={m} eta_PV={e:.6e}  pred {A*m*m/(24*np.pi):.6e} ratio {e/(A*m*m/(24*np.pi)):.4f}",flush=True); sys.exit()
x0=L/2; eps=np.array([0,0.03,-0.03,0.06,-0.06,0.09,-0.09])
def E(f): return c@E_vac(f,L,Lam,qc,eps.max(),Ms=Ms)
E0=E(lambda x:0*x); Cp=-A*m*m/(96*np.pi**2)
for w in W:
    g=lambda x:np.exp(-(x-x0)**2/(2*w*w))
    Es=np.array([E0 if e==0 else E(lambda x,e=e:e*g(x)) for e in eps])-E0; P=np.polyfit(eps,Es,6)
    xx=np.linspace(0,L,400001); dx=xx[1]-xx[0]; gg=g(xx); gp=np.gradient(gg,dx)
    Vol1=2*np.sum(gg)*dx; Vol2=2*np.sum(gg**2)*dx; I2=2*np.sum(gp**2)*dx; L1=P[-2]/Vol1
    print(f"m={m} PV w={w}: C/C_pred={(P[-3]-L1*Vol2)/I2/Cp:+.4f}  {time.time()-t:.0f}s",flush=True)
