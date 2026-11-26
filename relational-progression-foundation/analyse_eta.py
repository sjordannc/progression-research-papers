import json,numpy as np
A=6*np.log(2)-3*np.log(3)
for m in [0.2,0.14]:
    d=json.load(open(f"eta_m{m}.json")); n=np.array([int(k) for k in d]); R=n+0.5
    S=np.array([d[k]['S'][-1] for k in d]); P=np.array([d[k]['phi2'][-1] for k in d])
    etap=A*m*m/(48*np.pi); php=A*m*m/(16*np.pi**2)
    for lab,X in [("area+log+const",np.vstack([4*np.pi*R**2,np.log(R),np.ones_like(R)]).T),("area+const",np.vstack([4*np.pi*R**2,np.ones_like(R)]).T),("area+R+const",np.vstack([4*np.pi*R**2,R,np.ones_like(R)]).T)]:
        co,res,*_=np.linalg.lstsq(X,S,rcond=None); print(f"m={m} {lab}: eta={co[0]:.5e}  eta/pred={co[0]/etap:.4f}  other={co[1:]}  maxresid={np.max(np.abs(X@co-S)):.1e}")
    print(f"m={m} <phi^2>_PV at r=n: {P}  pred {php:.5e}  ratio {P/php}")
    print(f"m={m} covariant relation eta = pi<phi^2>/3 : eta_fit/(pi*phi2[-1]/3)=",end="")
    co,*_=np.linalg.lstsq(np.vstack([4*np.pi*R**2,np.log(R),np.ones_like(R)]).T,S,rcond=None); print(co[0]/(np.pi*P[-1]/3))
