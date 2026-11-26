import json,numpy as np,glob
A=6*np.log(2)-3*np.log(3)
phi2={0.2:2.18122e-4,0.14:1.06235e-4}; eta={0.2:2.26395e-4,0.14:1.10883e-4}
R={}
for fn in sorted(glob.glob("G_m*_f_*.json")):
    m=float(fn.split('_m')[1].split('_')[0]); xi=round(float(fn.split('_xi')[1].split('_')[0]),3)
    d=json.load(open(fn)); Cp=-A*m*m/(192*np.pi**2)
    w=np.array([float(k) for k in d]); y=np.array([v['a2']/v['Vol2'] for v in d.values()]); x=np.array([v['I2']/v['Vol2'] for v in d.values()])
    L1=np.mean([v['a1']/v['Vol1'] for v in d.values()])
    out=[]
    for lab,sel,cols in [("3 widest, (L2,C)",w>=24,lambda x:[np.ones_like(x),x]),("all, (L2,C,D x^2)",w>0,lambda x:[np.ones_like(x),x,x**2]),("w>=16, (L2,C,D x^2)",w>=16,lambda x:[np.ones_like(x),x,x**2])]:
        X=np.vstack(cols(x[sel])).T; co=np.linalg.lstsq(X,y[sel],rcond=None)[0]
        out.append(co[1]); print(f"m={m} xi={xi} {lab}: C/C_pred={co[1]/Cp:+.4f}  (L2-L1)/L1={(co[0]-L1)/L1:+.2e}")
    R[(m,xi)]=(np.array(out),Cp)
print()
for m in [0.2,0.14]:
    C0,Cp=R[(m,0.0)]; C6,_=R[(m,0.167)]
    print(f"m={m}: C(0)/C_pred {np.round(C0/Cp,3)};  C(1/6)/|C(0)| {np.round(C6/abs(C0),3)};  [C(0)-C(1/6)]/(-<phi2>/12) {np.round((C0-C6)/(-phi2[m]/12),3)}")
    print(f"        P2: -4pi C(0)/eta = {np.round(-4*np.pi*C0/eta[m],3)}")
C0,Cp=R[(0.28,0.0)]; print("m=0.28: C(0)/C_pred",np.round(C0/Cp,3))
