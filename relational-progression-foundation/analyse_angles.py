import numpy as np,json,glob,sys
sp=np.sqrt(np.pi)
def fit(d,ws=None,three=True):
    ws=sorted(float(k) for k in d) if ws is None else ws
    a2=np.array([d[str(w)]['a2'] for w in ws]); w=np.array(ws)
    V=w*sp; I2=sp/(2*w); I4=3*sp/(4*w**3)
    X=np.vstack([V,I2,I4]).T if three else np.vstack([V,I2]).T
    c,res,_,_=np.linalg.lstsq(X,a2,rcond=None); return c
eta={}
for o in ['100','110','111']:
    f=sorted(glob.glob(f'eta_{o}_full0.0_L40_n24.json')); eta[o]=json.load(open(f[0]))['eta']
for mode,tag in [('full','full0.0_L140_n24')]:
    for o in ['100','110','111']:
        S=json.load(open(f'G_{o}_S_{tag}.json')); G=json.load(open(f'G_{o}_G_{tag}.json'))
        cS3=fit(S); cS2=fit(S,[8.0,12.0,16.0],False); cG3=fit(G); cG2=fit(G,[8.0,12.0,16.0],False)
        C=cS3[1]/2; rho=-4*np.pi*C/eta[o]
        print(f"{o}: A_S={cS3[0]:.6f} Q_S={cS3[1]:+.5f} (2-par w>=8: {cS2[1]:+.5f}) D_S={cS3[2]:+.4f} | Q_G={cG3[1]:+.5f} (2-par {cG2[1]:+.5f}) A_G={cG3[0]:.6f} | C={C:+.5f} eta={eta[o]:.5f} rho={rho:+.3f}  Q_G/Q_S={cG3[1]/cS3[1]:+.3f}")
