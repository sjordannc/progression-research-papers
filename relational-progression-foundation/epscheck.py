import numpy as np,sys
from induced_common import *
m=float(sys.argv[1]); w=float(sys.argv[2]); r0=float(sys.argv[3]); Nr=int(sys.argv[4]); lmax=int(sys.argv[5]) if len(sys.argv)>5 else 3600
one=lambda r:np.ones_like(r); g=lambda r:np.exp(-(r-r0)**2/(2*w*w))
rr=np.linspace(0,Nr,400001); dr=rr[1]-rr[0]; I=lambda y:np.sum(y)*dr
Vol1=8*np.pi*I(rr**2*g(rr)); Vol2=4*np.pi*I(rr**2*g(rr)**2); I2=8*np.pi*I(np.gradient(rr*g(rr),dr)**2)
E0=E0_PV(lambda r:r,one,Nr,m,lmax=lmax)
A=6*np.log(2)-3*np.log(3); Cp=-A*m*m/(192*np.pi**2)
for eps in [0.005,0.01,0.02]:
    Ep=E0_PV(lambda r:r*(1+eps*g(r)),one,Nr,m,lmax=lmax); Em=E0_PV(lambda r:r*(1-eps*g(r)),one,Nr,m,lmax=lmax)
    a1=(Ep-Em)/(2*eps); a2=(Ep+Em-2*E0)/(2*eps**2); Lam=a1/Vol1
    print(f"eps={eps}: a1={a1:.8f} a2={a2:.8f}  Lam={Lam:.8e}  C/Cp={(a2-Lam*Vol2)/I2/Cp:.4f}",flush=True)
