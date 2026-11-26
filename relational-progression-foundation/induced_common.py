import numpy as np
from scipy.linalg import eigh_tridiagonal, eigvalsh_tridiagonal
JCUT=10
PV_c=np.array([1.,-3.,3.,-1.]); PV_M2=np.array([0.,1.,2.,3.])   # times m^2
def chain(l,f,fh,Nl,Nlh,M2,xiR):
    """tridiagonal of D^{1/2} V D^{1/2}; arrays indexed j=1..Nr stored at [j]"""
    Nr=len(f)-2; j=np.arange(1,Nr+1)
    d=(Nlh[j]*fh[j]**2+Nlh[j-1]*fh[j-1]**2*(j>1))/f[j]**2 + Nl[j]*(l*(l+1)/f[j]**2+M2+xiR[j])
    e=-Nlh[j[:-1]]*fh[j[:-1]]**2/(f[j[:-1]]*f[j[:-1]+1])
    s=np.sqrt(Nl[j]); return d*s*s, e*s[:-1]*s[1:]
def E0_PV(fn,Nfn,Nr,m,xi=0.0,lmax=800,Rfn=None,return_l=False):
    r=np.arange(0,Nr+2,dtype=float); f=fn(r); fh=fn(r+.5); Nl=Nfn(r); Nlh=Nfn(r+.5)
    xiR=xi*Rfn(r) if Rfn is not None else np.zeros_like(r)
    per=np.zeros(lmax)
    for l in range(lmax):
        tot=0.0
        for c,M2 in zip(PV_c,PV_M2*m*m):
            d,e=chain(l,f,fh,Nl,Nlh,M2,xiR)
            j0=max(0,l//JCUT-1)   # drop sites with angular gap l/r > JCUT (PV-suppressed; removes rounding noise)
            if j0>=len(d)-2: continue
            d=d[j0:]; e=e[j0:]
            w=eigvalsh_tridiagonal(d,e); tot+=c*0.5*np.sum(np.sqrt(np.abs(w)))
        per[l]=(2*l+1)*tot
    return (per.sum(),per) if return_l else per.sum()
def Ricci(fn,r,h=1e-3):
    f=fn(r); f1=(fn(r+h)-fn(r-h))/(2*h); f2=(fn(r+h)-2*f+fn(r-h))/h**2
    return 2/f**2*(1-f1**2-2*f*f2)
