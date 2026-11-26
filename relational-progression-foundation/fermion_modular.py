"""P4F: exact modular Hamiltonian of a ball for the massless lattice Dirac fermion (sector kappa), 100-digit precision,
compared with CHM 2 pi (1/2){W,H} and with candidate extra terms (pi/R) psi^dag psi and (pi/R) psibar psi.
Usage: fermion_modular.py dps Nr n kappa_list"""
import mpmath as mp, numpy as np, sys, time
mp.mp.dps=int(sys.argv[1]); Nr=int(sys.argv[2]); n=int(sys.argv[3]); kaps=[int(x) for x in sys.argv[4].split(',')]
t0=time.time()
def layout(kappa):
    if kappa<0: x=[mp.mpf(k)/2 for k in range(2,2*Nr+2)]       # u1 v1.5 ... uNr vNr+1/2
    else:       x=[mp.mpf(k)/2 for k in range(1,2*Nr+1)]       # v0.5 u1 ... uNr
    t=[0 if (2*xx)%2==0 else 1 for xx in x]; return x,t
def Bmat(kappa):
    x,t=layout(kappa); ui=[i for i in range(len(x)) if t[i]==0]; vi=[i for i in range(len(x)) if t[i]==1]
    upos={i:k for k,i in enumerate(ui)}; vpos={i:k for k,i in enumerate(vi)}
    B=mp.matrix(len(vi),len(ui))
    for i in range(len(x)-1):
        if t[i]==0: B[vpos[i+1],upos[i]]=-(x[i]/x[i+1])**kappa
        else:       B[vpos[i],upos[i+1]]=(x[i+1]/x[i])**kappa
    return B,[x[i] for i in ui],[x[i] for i in vi]
for Rshift in [mp.mpf(3)/4]:
  R=n+Rshift
  for kappa in kaps:
    B,xu,xv=Bmat(kappa); BtB=B.T*B; E,Q=mp.eigsy(BtB)
    S1=Q*mp.diag([1/mp.sqrt(e) for e in E])*Q.T; O=B*S1          # P_- = 1/2 [[I,-O^T],[-O,I]]
    nu=sum(1 for q in xu if q<R); nv=sum(1 for q in xv if q<R); n2=nu+nv
    CA=mp.matrix(n2,n2)
    for i in range(nu): CA[i,i]=mp.mpf(1)/2
    for i in range(nv): CA[nu+i,nu+i]=mp.mpf(1)/2
    for i in range(nv):
        for k in range(nu): CA[nu+i,k]=-O[i,k]/2; CA[k,nu+i]=-O[i,k]/2
    c,V=mp.eigsy(CA)
    h=V*mp.diag([mp.log((1-x)/x) for x in c])*V.T
    S=mp.fsum(-x*mp.log(x)-(1-x)*mp.log(1-x) for x in c)
    # CHM operator on region (u_1..u_n, v_1..v_n at r=j+1/2)
    rs=xu[:nu]+xv[:nv]
    w=[(R**2-q**2)/(2*R) for q in rs]
    HA=mp.matrix(n2,n2)
    for i in range(nv):
        for k in range(nu): HA[nu+i,k]=B[i,k]; HA[k,nu+i]=B[i,k]
    W=mp.diag(w); hC=mp.pi*(W*HA+HA*W)
    print(f"\n=== kappa={kappa} Nr={Nr} n={n} R={float(R)}  S={mp.nstr(S,8)}  min c={mp.nstr(min(c),3)} ({time.time()-t0:.0f}s)",flush=True)
    rows=[]
    probes=[]
    for r0,sg in [(0,8),(6,4),(10,5),(14,4),(18,3),(24,4),(24,8)]:
        for k in [0.3,0.6,1.0]:
            for ph in [0,1]:
                G=[mp.e**(-((r-r0)/sg)**2/2) for r in rs]
                vec=[G[i]*(mp.cos(k*rs[i]) if (i<nu)!=(ph==1) else mp.sin(k*rs[i])) for i in range(n2)]
                probes.append(((r0,sg,k,ph),mp.matrix(vec)))
        G=[mp.e**(-((r-r0)/sg)**2/2) for r in rs]
        probes.append(((r0,sg,'u-only'),mp.matrix([G[i] if i<nu else 0 for i in range(n2)])))
        probes.append(((r0,sg,'v-only'),mp.matrix([G[i] if i>=nu else 0 for i in range(n2)])))
    for lab,v in probes:
        Qe=(v.T*h*v)[0]; Qc=(v.T*hC*v)[0]; X1=mp.pi/R*mp.fsum(x**2 for x in v); X2=mp.pi/R*(mp.fsum(v[i]**2 for i in range(nu))-mp.fsum(v[nu+i]**2 for i in range(nv)))
        if len(lab)==3: print(f"   {lab}: Q_exact={float(Qe):+.3e}  Q_CHM={float(Qc):+.5f}  Q_exact/[(pi/R)sum psi^2]={float(Qe/X1):+.5f}  (scalar analog: 1)",flush=True)
        rows.append([float(Qe),float(Qc),float(X1),float(X2)])
    A=np.array(rows); y=A[:,0]
    for lab2,cols in [("CHM only",[1]),("CHM + psi^dag psi + psibar psi",[1,2,3])]:
        co,*_=np.linalg.lstsq(A[:,cols],y,rcond=None); res=np.abs(A[:,cols]@co-y)/np.abs(y)
        extra=""
        print(f"  fit [{lab2}]: coeffs={np.round(co,5)}  median rel resid={np.median(res):.2e} max={np.max(res):.2e}{extra}",flush=True)
    edge=np.array([lab[0]==24 for lab,_ in probes])
    sel=(np.abs(A[:,1])>0.05*np.max(np.abs(A[:,1])))&~edge; r=y[sel]/A[sel,1]; print(f"  Q_exact/Q_CHM over {sel.sum()} probes with non-negligible CHM: min {r.min():.4f} median {np.median(r):.4f} max {r.max():.4f}",flush=True)
    se=(np.abs(A[:,1])>0.05*np.max(np.abs(A[:,1])))&edge
    if se.any(): re=y[se]/A[se,1]; print(f"  EDGE probes (bump at R): Q_exact/Q_CHM min {re.min():.4f} median {np.median(re):.4f} max {re.max():.4f}",flush=True)
    print(f"  scalar comparison: the scalar's (pi/R)phi^2 term was 4-77% of Q for comparable bumps",flush=True)
print(f"total {time.time()-t0:.0f}s")
