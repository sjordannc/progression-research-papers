"""Geometry side: induced Einstein coefficient C from the PV-combined vacuum (Killing) energy E0[N,f] of the lattice field.
E0 = Lam*int N sqrt(h) + C*int N sqrt(h) R3 + (4-derivative), C = -1/(16 pi G_ind).
Usage: induced_G.py m xi mode   (mode f: spatial bumps; mode n: lapse-lapse and lapse-f cross terms)"""
import numpy as np,sys,time,json
from induced_common import *
m=float(sys.argv[1]); xi=float(sys.argv[2]); mode=sys.argv[3]
WS=[int(x) for x in sys.argv[4].split(',')] if len(sys.argv)>4 else None
Nr=360; r0=150.; eps=0.02; LS=[2400,3000,3600]; lmax=LS[-1]
one=lambda r:np.ones_like(r); rr=np.linspace(0,Nr,200001); dr=rr[1]-rr[0]
def integ(y): return np.sum(y)*dr
def tailsum(per):
    S=[per[:L].sum() for L in LS]; d1=S[1]-S[0]; d2=S[2]-S[1]
    # tail ~ L^-k : fit k from successive differences, extrapolate
    if d1*d2<=0 or abs(d2)>=abs(d1): return S[2],0.0
    q=d2/d1; tail=d2*q/(1-q); return S[2]+tail,abs(tail)
def R4(fn,Nfn,r,h=1e-3):
    f=fn(r); f1=(fn(r+h)-fn(r-h))/(2*h); N=Nfn(r); N1=(Nfn(r+h)-Nfn(r-h))/(2*h); N2=(Nfn(r+h)-2*N+Nfn(r-h))/h**2
    return Ricci(fn,r)-2*(N2+2*f1*N1/f)/N
def E(fn,Nfn):
    return E0_PV(fn,Nfn,Nr,m,xi=xi,lmax=lmax,Rfn=(lambda r:R4(fn,Nfn,r)) if xi!=0 else None,return_l=True)[1]
t=time.time(); base=E(lambda r:r,one); out={}
for w in (WS if WS else ([12,16,24,32] if mode=='f' else [16,24])):
    g=lambda r,w=w:np.exp(-(r-r0)**2/(2*w*w))
    rg=rr*g(rr); rgp=np.gradient(rg,dr); rgpp=np.gradient(rgp,dr)
    Vol1=8*np.pi*integ(rr**2*g(rr)); Vol2=4*np.pi*integ(rr**2*g(rr)**2); I2=8*np.pi*integ(rgp**2); I4=8*np.pi*integ(rgpp**2)
    if mode=='f':
        Ep=E(lambda r:r*(1+eps*g(r)),one); Em=E(lambda r:r*(1-eps*g(r)),one)
        a1,e1=tailsum((Ep-Em)/(2*eps)); a2,e2=tailsum((Ep+Em-2*base)/(2*eps**2))
        out[w]=dict(a1=a1,a2=a2,e1=e1,e2=e2,Vol1=Vol1,Vol2=Vol2,I2=I2,I4=I4)
    else:
        Nb=lambda r,s,w=w:1+s*eps*np.exp(-(r-r0)**2/(2*w*w))
        Enp=E(lambda r:r,lambda r:Nb(r,1)); Enm=E(lambda r:r,lambda r:Nb(r,-1))
        ann,enn=tailsum((Enp+Enm-2*base)/(2*eps**2)); an1,_=tailsum((Enp-Enm)/(2*eps))
        X={}
        for s1 in (1,-1):
            for s2 in (1,-1):
                X[(s1,s2)]=E(lambda r:r*(1+s2*eps*g(r)),lambda r:Nb(r,s1))
        ax,ex=tailsum((X[(1,1)]+X[(-1,-1)]-X[(1,-1)]-X[(-1,1)])/(4*eps*eps))
        nn=np.exp(-(rr-r0)**2/(2*w*w)); nnp=np.gradient(nn,dr)
        Jn1=4*np.pi*integ(rr**2*nn); Jx0=8*np.pi*integ(rr**2*nn*g(rr)); Jx2=16*np.pi*integ(rr*nnp*rgp); Jnn2=4*np.pi*integ(rr**2*nnp**2)
        out[w]=dict(an1=an1,ann=ann,enn=enn,ax=ax,ex=ex,Jn1=Jn1,Jx0=Jx0,Jx2=Jx2,Jnn2=Jnn2)
    print(w,out[w],f"{time.time()-t:.0f}s",flush=True)
json.dump({str(k):v for k,v in out.items()},open(f"G_m{m}_xi{xi:.4f}_{mode}_{sys.argv[4] if len(sys.argv)>4 else 'd'}.json","w"),indent=1)
