import json,numpy as np,glob
A=6*np.log(2)-3*np.log(3)
etaD={}
for m in (0.2,0.14,0.28):
    try:
        d=json.load(open(f"feta_m{m}.json")); n=np.array([int(k) for k in d]); R=n+0.75; S=np.array([d[k]['S'][-1] for k in d])
        X=np.vstack([4*np.pi*R**2,np.log(R),np.ones_like(R)]).T; etaD[m]=np.linalg.lstsq(X,S,rcond=None)[0][0]
    except FileNotFoundError: pass
print("eta_D:",{k:f"{v:.5e} ({v/(A*k*k/(24*np.pi)):.4f} of pred)" for k,v in etaD.items()})
for fn in sorted(glob.glob("FG_m*_f_*.json")):
    m=float(fn.split('_m')[1].split('_')[0]); d=json.load(open(fn)); Cp=-A*m*m/(96*np.pi**2); Lp=-4*(12*np.log(2)-9*np.log(3))*m**4/(64*np.pi**2)
    w=np.array([float(k) for k in d]); y=np.array([v['a2']/v['Vol2'] for v in d.values()]); x=np.array([v['I2']/v['Vol2'] for v in d.values()])
    L1=np.array([v['a1']/v['Vol1'] for v in d.values()])
    print(f"{fn}: Lambda1/Lambda_pred={np.mean(L1)/Lp:.4f} (spread {np.ptp(L1)/abs(np.mean(L1)):.1e})")
    for k in range(len(w)): print(f"   w={w[k]:.0f}: naive C(w)/C_pred (Lambda2=Lambda1) = {((y[k]-L1[k])/x[k])/Cp:+.4f}")
    for lab,sel,cols in [("3 widest (L2,C)",w>=24,lambda x:[np.ones_like(x),x]),("all (L2,C,Dx^2)",w>0,lambda x:[np.ones_like(x),x,x**2]),("w>=16 (L2,C,Dx^2)",w>=16,lambda x:[np.ones_like(x),x,x**2]),("w>=16 (L2,C)",w>=16,lambda x:[np.ones_like(x),x])]:
        if sel.sum()<len(cols(x[:1]))+0: continue
        co=np.linalg.lstsq(np.vstack(cols(x[sel])).T,y[sel],rcond=None)[0]
        ratio=""
        if m in etaD: ratio=f"   P2F: -4pi C/eta_D = {-4*np.pi*co[1]/etaD[m]:.4f}"
        print(f"   {lab}: C/C_pred={co[1]/Cp:+.4f}  (L2-L1)/L1={(co[0]-np.mean(L1))/np.mean(L1):+.2e}{ratio}")
for fn in sorted(glob.glob("FGn_merged.json")):
    m=0.2; d=json.load(open(fn)); Cp=-A*m*m/(96*np.pi**2); print(fn)
    for w,v in sorted(d.items(),key=lambda kv:float(kv[0])):
        Lam=v['an1']/v['Jn1']
        print(f"   w={w}: Lambda_n={Lam:.6e}  nn/(|C_pred| int(grad n)^2)={v['ann']/(abs(Cp)*v['Jnn2']):+.4f}  cross C_x/C_pred={(v['ax']-Lam*v['Jx0'])/v['Jx2']/Cp:+.4f}  cross/eta: {-4*np.pi*(v['ax']-Lam*v['Jx0'])/v['Jx2']/etaD[0.2]:.4f}")
