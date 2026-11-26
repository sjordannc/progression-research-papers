"""Entanglement side: PV-combined entanglement entropy of a sphere, S_PV(R)=sum_i c_i S(M_i;R), on the flat radial lattice,
and PV-combined <phi^2> at a point. Usage: induced_eta.py m"""
import numpy as np,sys,json,time
from induced_common import PV_c,PV_M2
m=float(sys.argv[1]); Rs=[20,30,40,50]; lmax=1500
def K(l,Nr,M2):
    j=np.arange(1,Nr+1,dtype=float); Kd=l*(l+1)/j**2+((j+.5)**2+(j-.5)**2*(j>1))/j**2+M2
    Ko=-(j[:-1]+.5)**2/(j[:-1]*(j[:-1]+1)); return np.diag(Kd)+np.diag(Ko,1)+np.diag(Ko,-1)
def S_and_phi2(Kmat,n,jstar):
    w,v=np.linalg.eigh(Kmat); Om=(v*np.sqrt(w))@v.T; Xjj=0.5*np.sum(v[jstar]**2/np.sqrt(w))
    A=Om[:n,:n];B=Om[:n,n:];C=Om[n:,n:]
    beta=0.5*B.T@np.linalg.solve(A,B); gam=C-beta
    xi=np.linalg.eigvals(np.linalg.solve(gam,beta)).real; xi=xi[(xi>1e-17)&(xi<1)]
    xp=xi/(1+np.sqrt(1-xi**2)); return float(np.sum(-np.log(1-xp)-xp/(1-xp)*np.log(xp))),Xjj
out={}; t=time.time()
for n in Rs:
    Nr=n+110; jstar=n-1  # <phi^2> evaluated at r=n (inside the ball, flat space)
    perS=np.zeros(lmax); perP=np.zeros(lmax)
    for l in range(lmax):
        for c,M2 in zip(PV_c,PV_M2*m*m):
            s,x=S_and_phi2(K(l,Nr,M2),n,jstar); perS[l]+=c*(2*l+1)*s; perP[l]+=c*(2*l+1)*x/(4*np.pi*n**2)
    out[n]=dict(S=[perS[:L].sum() for L in (500,1000,1500)],phi2=[perP[:L].sum() for L in (500,1000,1500)])
    print(n,out[n],f"{time.time()-t:.0f}s",flush=True)
json.dump(out,open(f"eta_m{m}.json","w"),indent=1)
