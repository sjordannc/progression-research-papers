import numpy as np
from scipy.interpolate import BSpline
rr=np.linspace(0.25,300,60000); dr=rr[1]-rr[0]
def bder(basis_list,x):
    """values, first, second derivatives of each basis function on grid x"""
    out=[]
    for f,lo,hi in basis_list:
        if not f.__defaults__: out.append((np.ones_like(x),np.zeros_like(x),np.zeros_like(x))); continue
        b=f.__defaults__[0]; out.append((np.nan_to_num(b(x)),np.nan_to_num(b.derivative(1)(x)),np.nan_to_num(b.derivative(2)(x))))
    return out
class GR:
    def __init__(s,basis,C,mu):
        s.b=basis; s.C=C; s.mu=mu(rr); s.Bn=bder(basis.nu,rr); s.Bs=bder(basis.s,rr)
    def fields(s,a):
        nu=sum(c*B[0] for c,B in zip(a[:s.b.nn],s.Bn)); S=sum(c*B[0] for c,B in zip(a[s.b.nn:],s.Bs))
        S1=sum(c*B[1] for c,B in zip(a[s.b.nn:],s.Bs)); S2=sum(c*B[2] for c,B in zip(a[s.b.nn:],s.Bs))
        return nu,S,S1,S2
    def E(s,a):
        nu,S,S1,S2=s.fields(a); f=rr*(1+S); fp=1+S+rr*S1; fpp=2*S1+rr*S2; N=np.exp(nu)
        return s.C*4*np.pi*np.sum(N*(-4*f*fpp+2*(1-fp**2)))*dr+np.sum(N*s.mu)*dr
    def hess(s,h=1e-4):
        n=s.b.n; H=np.zeros((n,n)); E0=s.E(np.zeros(n))
        for i in range(n):
            for j in range(i,n):
                a=np.zeros(n); 
                def Ev(x,y): a2=np.zeros(n); a2[i]+=x; a2[j]+=y; return s.E(a2)
                H[i,j]=H[j,i]=(Ev(h,h)-Ev(h,-h)-Ev(-h,h)+Ev(-h,-h))/(4*h*h)
        return H
    def Jsrc(s):
        return np.array([np.sum(s.mu*B[0])*dr for B in s.Bn]+[0.0]*s.b.ns)
