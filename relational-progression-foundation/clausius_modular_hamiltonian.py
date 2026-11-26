"""High-precision Clausius toy model.
Environment: infinite harmonic chain (discretised 1+1 massless scalar, IR mass m), vacuum.
Region A: L sites to the right of a cut. Clocks: static observers at distance x from the cut,
whose Rindler (boost-time) lapse is Pi(x)=a*x; Unruh/Tolman local temperature T(x)=1/(2 pi x).
"""
import mpmath as mp, numpy as np, time, sys
mp.mp.dps=int(sys.argv[1]) if len(sys.argv)>1 else 40
L=int(sys.argv[2]) if len(sys.argv)>2 else 120
m=mp.mpf('1e-4')
t0=time.time()
def corr(d):
    fx=lambda k: mp.cos(k*d)/(2*mp.sqrt(m**2+4*mp.sin(k/2)**2))
    fp=lambda k: mp.cos(k*d)*mp.sqrt(m**2+4*mp.sin(k/2)**2)/2
    return mp.quad(fx,[0,mp.mpf('1e-3'),0.1,mp.pi])/mp.pi, mp.quad(fp,[0,mp.pi])/mp.pi
cx=[];cp=[]
for d in range(L):
    a,b=corr(d); cx.append(a); cp.append(b)
X=mp.matrix(L,L); P=mp.matrix(L,L)
for i in range(L):
    for j in range(L):
        X[i,j]=cx[abs(i-j)]; P[i,j]=cp[abs(i-j)]
print("correlators done %.0fs"%(time.time()-t0),flush=True)
def sqrtm_sym(A):
    E,Q=mp.eigsy(A); return Q*mp.diag([mp.sqrt(e) for e in E])*Q.T, Q, E
Xh,Qx,Ex=sqrtm_sym(X)
Xhi=Qx*mp.diag([1/mp.sqrt(e) for e in Ex])*Qx.T
W=Xh*P*Xh; E2,Q=mp.eigsy(W)
nu=[mp.sqrt(e) for e in E2]
eps=[mp.log((2*v+1)/(2*v-1)) for v in nu]
def S_of(nus): return mp.fsum((v+mp.mpf(.5))*mp.log(v+mp.mpf(.5))-(v-mp.mpf(.5))*mp.log(v-mp.mpf(.5)) for v in nus)
S0=S_of(nu)
Mx=Xhi*Q*mp.diag([v*e for v,e in zip(nu,eps)])*Q.T*Xhi
Mp=Xh*Q*mp.diag([e/v for v,e in zip(nu,eps)])*Q.T*Xh
print("modular Hamiltonian done %.0fs, S0=%s"%(time.time()-t0,mp.nstr(S0,8)),flush=True)
x=[i+0.5 for i in range(L)]
print("\n== Modular Hamiltonian profile: coefficient of pi_i^2/2 vs CFT interval prediction 2 pi x (L-x)/L")
for i in [0,1,2,4,8,16,32,48,64,90]:
    if i<L:
        pred=2*np.pi*x[i]*(L-x[i])/L
        print(f" x={x[i]:5.1f}  M_pi={float(Mp[i,i]):9.4f}  2pi x(L-x)/L={pred:9.4f}  ratio={float(Mp[i,i])/pred:.4f}")
print(" bond couplings -M_phi[i,i+1] vs 2 pi xb (L-xb)/L:")
for i in [k for k in [0,2,8,16,32,64] if k<L-1]:
    xb=i+1.0; pred=2*np.pi*xb*(L-xb)/L
    print(f"  xb={xb:5.1f}  {float(-Mx[i,i+1]):9.4f}  {pred:9.4f}  ratio={float(-Mx[i,i+1])/pred:.4f}")
np.save('./Mx.npy',np.array([[float(Mx[i,j]) for j in range(L)] for i in range(L)]))
np.save('./Mp.npy',np.array([[float(Mp[i,j]) for j in range(L)] for i in range(L)]))
# ---- exact first-law / Clausius check for small static heat deposits (phi-noise bumps)
print("\n== Static heat deposits: classical-noise bump of width w at distance x0 (epsilon -> 0)")
print(" x0   dS/eps (exact entropy)   d<K>/eps (first law)   2pi*sum x dE/eps (Clausius, T=1/2pi x)   CFT-interval weighting")
w=3.0
for x0 in [6,12,24,40]:
    u=[mp.e**(-((i+0.5-x0)/w)**2/2) for i in range(L)]
    # energy of the deposit per unit eps: gradient energy of u (assign bond i,i+1 to xb=i+1) + mass term
    eb=[ (u[i+1]-u[i])**2/2 for i in range(L-1)]
    dK=mp.fsum(u[i]*Mx[i,j]*u[j] for i in range(L) for j in range(L))/2
    BW=2*mp.pi*mp.fsum((i+1)*eb[i] for i in range(L-1))
    BWi=2*mp.pi*mp.fsum((i+1)*(L-i-1)/mp.mpf(L)*eb[i] for i in range(L-1))
    res=[]
    for e in [mp.mpf('1e-12'),mp.mpf('1e-16')]:
        Xe=X+e*mp.matrix(u)*mp.matrix(u).T
        Xeh,_,_=sqrtm_sym(Xe); Ee,_=mp.eigsy(Xeh*P*Xeh)
        res.append((S_of([mp.sqrt(z) for z in Ee])-S0)/e)
    print(f" {x0:3d}  {mp.nstr(res[0],7):>12} / {mp.nstr(res[1],7):>12}   {mp.nstr(dK,7):>12}   {mp.nstr(BW,7):>12}   {mp.nstr(BWi,7):>12}",flush=True)
print("total %.0fs"%(time.time()-t0))
