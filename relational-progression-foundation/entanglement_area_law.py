import numpy as np
N=60
def K_l(l):
    K=np.zeros((N,N))
    for j in range(1,N+1):
        i=j-1
        K[i,i]+= l*(l+1)/j**2
        K[i,i]+= (j+0.5)**2/j**2
        if j>1: K[i,i]+= (j-0.5)**2/j**2
        if j<N:
            K[i,i+1]=K[i+1,i]=-(j+0.5)**2/(j*(j+1))
    return K
def S_of(K,n):
    w,v=np.linalg.eigh(K); Om=(v*np.sqrt(w))@v.T
    A=Om[:n,:n];Bm=Om[:n,n:];C=Om[n:,n:]
    beta=0.5*Bm.T@np.linalg.solve(A,Bm); gam=C-beta
    xi=np.linalg.eigvals(np.linalg.solve(gam,beta)).real
    xi=xi[(xi>1e-14)&(xi<1)]
    xp=xi/(1+np.sqrt(1-xi**2))
    return float(np.sum(-np.log(1-xp)-xp/(1-xp)*np.log(xp)))
ns=list(range(6,31,4)); Ls=range(0,1500)
S=np.zeros(len(ns))
for l in Ls:
    K=K_l(l); sl=np.array([S_of(K,n) for n in ns]); S+= (2*l+1)*sl
    if l>50 and sl.max()*(2*l+1)<1e-10: break
R=np.array(ns)+0.5
c=np.polyfit(R**2,S,1)
for r,s in zip(R,S): print(f"R/a={r:5.1f}  S={s:8.2f}  S/(R/a)^2={s/r**2:.3f}")
print("fit S = %.3f (R/a)^2 + %.2f, last l used %d"%(c[0],c[1],l))
