"""Test A: does the entanglement entropy of a ball respond to a change of PHYSICAL area with the
flat-space area-law coefficient?  Scalar field on the spatial metric dl^2 = dr^2 + f(r)^2 dOmega^2,
Srednicki radial lattice generalised: H_l = 1/2 sum[pi_j^2 + f_{j+1/2}^2 (phi_j/f_j - phi_{j+1}/f_{j+1})^2 + l(l+1)/f_j^2 phi_j^2]
Prediction (pre-registered): dS = kappa * (f(R)^2 - R^2), kappa = flat coefficient (0.295).
Null alternatives: dS ~ 0 (entropy tracks coordinate area) or a different coefficient."""
import numpy as np, sys
N=int(sys.argv[1]) if len(sys.argv)>1 else 60
def K_l(l,f,fh):
    K=np.zeros((N,N))
    for j in range(1,N+1):
        i=j-1
        K[i,i]+= l*(l+1)/f[j]**2 + fh[j]**2/f[j]**2
        if j>1: K[i,i]+= fh[j-1]**2/f[j]**2
        if j<N: K[i,i+1]=K[i+1,i]=-fh[j]**2/(f[j]*f[j+1])
    return K
def S_cut(K,n):
    w,v=np.linalg.eigh(K); Om=(v*np.sqrt(w))@v.T
    A=Om[:n,:n];B=Om[:n,n:];C=Om[n:,n:]
    beta=0.5*B.T@np.linalg.solve(A,B); gam=C-beta
    xi=np.linalg.eigvals(np.linalg.solve(gam,beta)).real; xi=xi[(xi>1e-15)&(xi<1)]
    xp=xi/(1+np.sqrt(1-xi**2)); return float(np.sum(-np.log(1-xp)-xp/(1-xp)*np.log(xp)))
def total_S(fn,n,lmax=3000):
    r=np.arange(0,N+2,dtype=float); f=fn(r); fh=fn(r+0.5)   # f[j] at site j, fh[j] at j+1/2
    S=0.0
    for l in range(lmax):
        s=S_cut(K_l(l,f,fh),n); S+=(2*l+1)*s
        if l>30 and (2*l+1)*s<1e-9*max(S,1): break
    return S,l
kappa=0.295
print(" R    eps   width   S_flat     S_deformed   dS        kappa*dA/a^2   ratio   (lmax)")
for n in [10,15,20,25]:
    R=n+0.5
    S0,l0=total_S(lambda r: r, n)
    for eps,wid in [(0.02,3.0),(0.02,6.0),(-0.02,3.0),(0.05,3.0)]:
        fn=lambda r,eps=eps,wid=wid,R=R: r*(1+eps*np.exp(-((r-R)/wid)**2))
        S1,l1=total_S(fn,n)
        dA=fn(np.array([R]))[0]**2-R**2
        print(f"{R:5.1f} {eps:+.2f} {wid:5.1f}  {S0:9.3f}  {S1:9.3f}  {S1-S0:+8.4f}  {kappa*dA:+9.4f}   {(S1-S0)/(kappa*dA):.3f}   ({l0},{l1})",flush=True)
