"""Test B: is the exact lattice modular Hamiltonian of a ball the CHM operator 2 pi int w T00, w=(R^2-r^2)/(2R)?
And with WHICH stress tensor: minimal T00, or improved T00 - xi lap(phi^2), xi=1/6 (adds (pi/R) int psi^2 dV)?
3D massless scalar, Srednicki radial lattice, sector l (u_j = r_j * psi_l(r_j)), Dirichlet at r=N+1.
Probe: classical-noise deposit X -> X + eps v v^T with v_j = j*psi(j), psi a bump strictly inside the ball.
Pre-registered: Q_exact = 1/2 v^T M_x v  should equal Q_imp (support) ; Q_min deviating by the predicted
(pi/R) sum v^2 => minimal ruled out. Contradiction: Q_exact matches neither within finite-size error.
Control: momentum deposit, 1/2 p^T M_p p vs 2 pi sum w_j p_j^2/2 (no improvement term enters)."""
import mpmath as mp, sys, time
mp.mp.dps=int(sys.argv[1]) if len(sys.argv)>1 else 120
N=int(sys.argv[2]) if len(sys.argv)>2 else 80
n=int(sys.argv[3]) if len(sys.argv)>3 else 24
ls=[int(s) for s in (sys.argv[4] if len(sys.argv)>4 else "0,1,2").split(",")]
R=mp.mpf(n)+mp.mpf(1)/2; t0=time.time()
def Kmat(l):
    K=mp.matrix(N,N)
    for j in range(1,N+1):
        i=j-1; K[i,i]=mp.mpf(l*(l+1))/j**2+(j+mp.mpf(1)/2)**2/mp.mpf(j)**2
        if j>1: K[i,i]+=(j-mp.mpf(1)/2)**2/mp.mpf(j)**2
        if j<N: K[i,i+1]=K[i+1,i]=-(j+mp.mpf(1)/2)**2/(mp.mpf(j)*(j+1))
    return K
def fun(A,f):
    E,Q=mp.eigsy(A); return Q*mp.diag([f(e) for e in E])*Q.T
w=lambda r:(R**2-mp.mpf(r)**2)/(2*R)
sc=n/24.0
bumps=[(round(a*sc),max(1,round(b*sc))) for a,b in [(0,12),(0,24),(12,8),(18,6),(24,4),(24,8)]]
for l in ls:
    K=Kmat(l)
    X=fun(K,lambda e:1/(2*mp.sqrt(e))); P=fun(K,lambda e:mp.sqrt(e)/2)
    XA=X[0:n,0:n]; PA=P[0:n,0:n]
    E,Qx=mp.eigsy(XA); Xh=Qx*mp.diag([mp.sqrt(e) for e in E])*Qx.T; Xhi=Qx*mp.diag([1/mp.sqrt(e) for e in E])*Qx.T
    E2,Q=mp.eigsy(Xh*PA*Xh); nu=[mp.sqrt(e) for e in E2]
    ep=[mp.log((2*v+1)/(2*v-1)) for v in nu]
    Mx=Xhi*Q*mp.diag([v*e for v,e in zip(nu,ep)])*Q.T*Xhi
    Mp=Xh*Q*mp.diag([e/v for v,e in zip(nu,ep)])*Q.T*Xh
    print(f"\n=== l={l}  N={N} n={n} R={float(R)}  min(nu-1/2)={mp.nstr(min(nu)-0.5,3)}  max eps={mp.nstr(max(ep),4)}  ({time.time()-t0:.0f}s)",flush=True)
    rows=[]
    print(" bump(r0,sig) |  Q_exact    Q_min     Q_imp   | exact/min  exact/imp | predicted imp/min")
    for r0,s in bumps:
        v=[j*mp.e**(-((j-r0)/mp.mpf(s))**2/2) for j in range(1,n+1)]
        Qe=mp.fsum(v[a]*Mx[a,b]*v[b] for a in range(n) for b in range(n))/2
        # minimal CHM: gradient bonds at r=j+1/2 (bond j,j+1 inside ball for j<n; v_{n+1}=0 bond at R has w=0),
        # centrifugal term at sites
        g=mp.fsum(w(j+mp.mpf(1)/2)*(j+mp.mpf(1)/2)**2*(v[j-1]/j-(v[j]/(j+1) if j<n else 0))**2/2 for j in range(1,n+1))
        c=mp.fsum(w(j)*l*(l+1)*v[j-1]**2/(2*mp.mpf(j)**2) for j in range(1,n+1))
        Qm=2*mp.pi*(g+c); Qi=Qm+mp.pi/R*mp.fsum(x**2 for x in v); uR=R*mp.e**(-((R-r0)/mp.mpf(s))**2/2)
        Qb=Qi-mp.pi/3*uR**2
        rows.append((float(Qe),float(Qm),float(Qi-Qm)));print(f'   edge: Q_ex/Q_imp_bulk={float(Qe/Qi):.4f}  Q_ex/Q_imp_with_edge={float(Qe/Qb):.4f}  edge share={float((Qi-Qb)/Qi):.3f}')
        print(f" ({r0:2d},{s})      | {float(Qe):9.4f} {float(Qm):9.4f} {float(Qi):9.4f} |  {float(Qe/Qm):.4f}    {float(Qe/Qi):.4f}   |  {float(Qi/Qm):.4f}",flush=True)
    import numpy as np
    A=np.array([[r[1],r[2]] for r in rows]); y=np.array([r[0] for r in rows])
    (a,b),*_=np.linalg.lstsq(A,y,rcond=None); res=np.max(np.abs(A@[a,b]-y)/y)
    print(f" FIT Q_exact = a*Q_min + b*Q_impterm:  a={a:.4f}  b={b:.4f}  b/a={b/a:.4f}  max rel resid={res:.2e}  (improved: b/a=1, minimal: b/a=0)")
    print(" momentum control: 1/2 p^T M_p p  vs 2pi sum w p^2/2")
    for r0,s in bumps[2:5]:
        p=[mp.e**(-((j-r0)/mp.mpf(s))**2/2) for j in range(1,n+1)]
        Qe=mp.fsum(p[a]*Mp[a,b]*p[b] for a in range(n) for b in range(n))/2
        Qc=2*mp.pi*mp.fsum(w(j)*p[j-1]**2/2 for j in range(1,n+1))
        print(f"   ({r0:2d},{s})  {float(Qe):9.4f} {float(Qc):9.4f} ratio {float(Qe/Qc):.4f}",flush=True)
    if False:
        r0,s=bumps[2]; v=mp.matrix([j*mp.e**(-((j-r0)/mp.mpf(s))**2/2) for j in range(1,n+1)])
        S=lambda ns: mp.fsum((x+mp.mpf(.5))*mp.log(x+mp.mpf(.5))-(x-mp.mpf(.5))*mp.log(x-mp.mpf(.5)) for x in ns)
        S0=S(nu); eps=mp.mpf(10)**(-mp.mp.dps//3)
        Xe=XA+eps*v*v.T; Ee,Qe_=mp.eigsy(Xe); Xeh=Qe_*mp.diag([mp.sqrt(e) for e in Ee])*Qe_.T
        E3,_=mp.eigsy(Xeh*PA*Xeh); dS=(S([mp.sqrt(e) for e in E3])-S0)/eps
        print(f" exact dS/eps for bump {bumps[2]} = {mp.nstr(dS,8)}   (S0={mp.nstr(S0,6)})")
print(f"total {time.time()-t0:.0f}s")
