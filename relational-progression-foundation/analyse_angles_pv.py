import numpy as np,json,os
sp=np.sqrt(np.pi); A=6*np.log(2)-3*np.log(3); m=0.2; Cp=-A*m*m/(96*np.pi**2); Lam=1.5893e-5
for o in ['100','110','111']:
    r={}
    for k in 'SG':
        f=f'G_{o}_{k}_PV0.2_L240_n24.json'
        if not os.path.exists(f): continue
        d=json.load(open(f)); ws=np.array(sorted(float(x) for x in d)); a2=np.array([d[str(w)]['a2'] for w in ws])
        if len(ws)<3: continue
        X=np.vstack([ws*sp,sp/(2*ws)]).T; r[k]=np.linalg.lstsq(X,a2,rcond=None)[0]
    if 'S' in r: print(o,'C_PV/C_pred=%.4f'%(r['S'][1]/2/Cp),'volume A_S/(2 Lam)=%.4f'%(r['S'][0]/(2*Lam)),end=' ')
    if 'S' in r and 'G' in r: print('Q_G/Q_S=%+.4f'%(r['G'][1]/r['S'][1]),'A_G/(Lam/2)=%.4f'%(r['G'][0]/(Lam/2)),end='')
    print()
