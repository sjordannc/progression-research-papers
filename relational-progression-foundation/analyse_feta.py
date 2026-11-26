import json,numpy as np
A=6*np.log(2)-3*np.log(3)
etaS={0.2:2.26395e-4,0.14:1.10883e-4}   # scalar, note 5 (area+log+const)
for m in [0.2,0.14]:
    d=json.load(open(f"feta_m{m}.json")); n=np.array([int(k) for k in d]); R=n+0.75; S=np.array([d[k]['S'][-1] for k in d])
    etap=A*m*m/(24*np.pi)
    for lab,X in [("area+log+const",np.vstack([4*np.pi*R**2,np.log(R),np.ones_like(R)]).T),("area+R+const",np.vstack([4*np.pi*R**2,R,np.ones_like(R)]).T),("area+const",np.vstack([4*np.pi*R**2,np.ones_like(R)]).T)]:
        co=np.linalg.lstsq(X,S,rcond=None)[0]; print(f"m={m} {lab}: eta_D={co[0]:.5e} eta_D/pred={co[0]/etap:.4f}  eta_D/eta_S={co[0]/etaS[m]:.4f} (pred 2)  maxres={np.max(np.abs(X@co-S)):.1e}")
