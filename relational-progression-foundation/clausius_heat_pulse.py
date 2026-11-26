"""Heat pulse crossing the cut: first-law entropy change d<K_A>(t) (modular Hamiltonian from
clausius_mp.py, 200-digit precision) vs the Clausius sum  sum_x dE(x,t)/T(x)  with T(x)=1/beta(x),
beta(x)=2 pi x (L-x)/L (interval form of the Unruh/Bisognano-Wichmann weight)."""
import numpy as np
Mx=np.load('./Mx.npy'); Mp=np.load('./Mp.npy'); L=Mx.shape[0]
N=1200; c=600; m=1e-4
V=np.diag(np.full(N,2.0+m*m)); V-=np.diag(np.ones(N-1),1)+np.diag(np.ones(N-1),-1)
w2,O=np.linalg.eigh(V); w=np.sqrt(w2)
k=0.25; om=2*np.sin(k/2); vg=np.cos(k/2); width=5.0; x0=c-25
j=np.arange(N); g=np.exp(-((j-x0)/width)**2/2)
u=g*np.cos(k*j); v=om*g*np.sin(k*j)
uh=O.T@u; vh=O.T@v
A=slice(c,c+L); xs=np.arange(L)+0.5; xb=np.arange(1,L)
beta_s=2*np.pi*xs*(L-xs)/L; beta_b=2*np.pi*xb*(L-xb)/L
print(" t    energy in A   d<K_A>/eps (first law)   sum dE/T /eps (Clausius)   ratio")
for t in range(0,61,5):
    C=np.cos(w*t); S=np.sin(w*t)
    ut=O@(C*uh+S/w*vh); vt=O@(-w*S*uh+C*vh)
    ua=ut[A]; va=vt[A]
    dK=0.5*(ua@Mx@ua+va@Mp@va)
    es=0.5*va**2+0.5*m*m*ua**2; eb=0.5*np.diff(ua)**2
    EA=es.sum()+eb.sum()
    Cl=np.sum(beta_s*es)+np.sum(beta_b*eb)
    print(f"{t:3d}  {EA:11.4e}  {dK:14.5e}  {Cl:14.5e}   {dK/Cl if Cl>1e-12 else float('nan'):.4f}")
