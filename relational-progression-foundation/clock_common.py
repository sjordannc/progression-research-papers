"""Two-clock test: gradients of the Pauli-Villars Dirac vacuum energy E0[N,f] on the note-6 radial lattice by
Hellmann-Feynman, for N = exp(nu), f = r(1+s), nu and s expanded in cubic B-splines."""
import numpy as np, os
from scipy.linalg import eigh_tridiagonal
from scipy.interpolate import BSpline
from fermion_common import sites, JCUT, PV_c, PV_M2
_gx,_gw=np.polynomial.legendre.leggauss(6)
NR=int(os.environ.get('NR','300')); MPV=0.2
def bspline_basis(lo,hi,h):
    """cubic B-splines with uniform knots spacing h whose support lies in [lo-?, hi]; returns list of callables"""
    kn=np.arange(lo-3*h,hi+h/2,h); out=[]
    for i in range(len(kn)-4):
        t=kn[i:i+5]
        if t[-1]>hi+1e-9 or t[-1]<=max(lo,0)+1e-9: continue
        b=BSpline.basis_element(t,extrapolate=False); out.append((lambda x,b=b:np.nan_to_num(b(x)), t[0], t[-1]))
    return out
class Basis:
    def __init__(s,h=10.0,nu_hi=270.,s_lo=0.,s_hi=270.):
        # nu: a constant (reference-clock normalisation) plus splines vanishing with two derivatives at r=0 (regular lapse);
        # s: splines with support in r>=0, so f = r(1+s) is regular at the origin
        sp=[b for b in bspline_basis(0.0,nu_hi,h) if b[1]>=-1e-9]
        s.nu=[(lambda x:np.ones_like(np.asarray(x,float)),0.0,np.inf)]+sp; s.s=[b for b in bspline_basis(s_lo,s_hi,h) if b[1]>=s_lo-1e-9]
        s.nn=len(s.nu); s.ns=len(s.s); s.n=s.nn+s.ns
    def nu_of(s,a,x): return sum(c*b[0](x) for c,b in zip(a[:s.nn],s.nu))
    def s_of(s,a,x): return sum(c*b[0](x) for c,b in zip(a[s.nn:],s.s))
def chain_grad(K_list,a,basis,masses_c=None):
    """returns (E0 contribution, gradient wrt basis coeffs) summed over the given K (both signs, PV masses)."""
    g=np.zeros(basis.n); E=0.0
    for K in K_list:
        if K/JCUT-1>=NR-2: continue
        for kappa in (K,-K):
            xmin=max(0.0,abs(kappa)/JCUT-1); x,t=sites(kappa,NR,xmin)
            N=np.exp(basis.nu_of(a,x)); xa=x[:-1]; xb=x[1:]; mid=(xa+xb)/2; half=(xb-xa)/2
            Q=[(mid[:,None]+half[:,None]*_gx[None,:]),half[:,None]*_gw[None,:]]
            sq=basis.s_of(a,Q[0]); I=np.sum(Q[1]/(Q[0]*(1+sq)),axis=1)
            tau=np.where(t[:-1]==0,-1.0,1.0); e=tau*np.exp(tau*kappa*I)*np.sqrt(N[:-1]*N[1:])
            sig=np.where(t==0,1.0,-1.0)
            dEdN=np.zeros(len(x)); dEdI=np.zeros(len(I))
            for c,M2 in zip(PV_c,PV_M2*MPV*MPV):
                M=np.sqrt(M2); d=sig*M*N
                w,v=eigh_tridiagonal(d,e); sg=np.sign(w)
                Qd=(v*v)@sg; Qo=(v[:-1]*v[1:])@sg
                W=2*K*c; E+=W*(-0.5)*np.sum(np.abs(w))
                dEd=W*(-0.5)*Qd; dEe=W*(-1.0)*Qo
                dEdN+=dEd*sig*M; dEdN[:-1]+=dEe*e/(2*N[:-1]); dEdN[1:]+=dEe*e/(2*N[1:]); dEdI+=dEe*tau*kappa*e
            for k,b in enumerate(basis.nu): g[k]+=np.sum(dEdN*N*b[0](x))
            for k,b in enumerate(basis.s):
                dI=-np.sum(Q[1]*b[0](Q[0])/(Q[0]*(1+sq)**2),axis=1); g[basis.nn+k]+=np.sum(dEdI*dI)
    return E,g
def grad_E0(a,basis,nproc=4,Kmax=None):
    from multiprocessing import Pool
    Kmax=Kmax or int(JCUT*(NR-1)); Ks=np.arange(1,Kmax+1)
    chunks=[Ks[i::nproc*4] for i in range(nproc*4)]
    with Pool(nproc) as p: res=p.starmap(chain_grad,[(list(c),a,basis) for c in chunks])
    return sum(r[0] for r in res), sum(r[1] for r in res)

def chain_grad_fn(K_list,nu_fn,s_fn,tests_nu,tests_s):
    """E0 and its directional derivatives along test functions, for fields given as callables."""
    gn=np.zeros(len(tests_nu)); gs=np.zeros(len(tests_s)); E=0.0
    for K in K_list:
        if K/JCUT-1>=NR-2: continue
        for kappa in (K,-K):
            xmin=max(0.0,abs(kappa)/JCUT-1); x,t=sites(kappa,NR,xmin)
            N=np.exp(nu_fn(x)); xa=x[:-1]; xb=x[1:]; mid=(xa+xb)/2; half=(xb-xa)/2
            q=mid[:,None]+half[:,None]*_gx[None,:]; qw=half[:,None]*_gw[None,:]
            sq=s_fn(q); I=np.sum(qw/(q*(1+sq)),axis=1)
            tau=np.where(t[:-1]==0,-1.0,1.0); e=tau*np.exp(tau*kappa*I)*np.sqrt(N[:-1]*N[1:])
            sig=np.where(t==0,1.0,-1.0); dEdN=np.zeros(len(x)); dEdI=np.zeros(len(I))
            for c,M2 in zip(PV_c,PV_M2*MPV*MPV):
                M=np.sqrt(M2); d=sig*M*N; w,v=eigh_tridiagonal(d,e); sg=np.sign(w)
                Qd=(v*v)@sg; Qo=(v[:-1]*v[1:])@sg; W=2*K*c; E+=W*(-0.5)*np.sum(np.abs(w))
                dEd=W*(-0.5)*Qd; dEe=W*(-1.0)*Qo
                dEdN+=dEd*sig*M; dEdN[:-1]+=dEe*e/(2*N[:-1]); dEdN[1:]+=dEe*e/(2*N[1:]); dEdI+=dEe*tau*kappa*e
            for k,T in enumerate(tests_nu): gn[k]+=np.sum(dEdN*N*T(x))
            w2=qw/(q*(1+sq)**2)
            for k,T in enumerate(tests_s): gs[k]+=np.sum(dEdI*(-np.sum(w2*T(q),axis=1)))
    return E,gn,gs
def grad_fn(nu_fn,s_fn,tests_nu,tests_s,nproc=4):
    from multiprocessing import Pool
    Ks=np.arange(1,int(JCUT*(NR-1))+1); chunks=[list(Ks[i::nproc*6]) for i in range(nproc*6)]
    with Pool(nproc) as p: res=p.starmap(chain_grad_fn,[(c,nu_fn,s_fn,tests_nu,tests_s) for c in chunks])
    return sum(r[0] for r in res),sum(r[1] for r in res),sum(r[2] for r in res)
class Tab:
    """picklable tabulated function (linear interpolation on a fine grid; constant extrapolation)"""
    def __init__(s,x,y): s.x=np.asarray(x,float); s.y=np.asarray(y,float)
    def __call__(s,r): return np.interp(r,s.x,s.y)
class Spl:
    def __init__(s,t): s.b=BSpline.basis_element(np.asarray(t,float),extrapolate=False)
    def __call__(s,r): return np.nan_to_num(s.b(r))
def test_splines(lo,hi,h):
    kn=np.arange(lo,hi+h/2,h); return [Spl(kn[i:i+5]) for i in range(len(kn)-4)],[ (kn[i],kn[i+4]) for i in range(len(kn)-4)]
