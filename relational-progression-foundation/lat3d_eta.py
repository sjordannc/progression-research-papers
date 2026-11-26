"""usage: lat3d_eta.py orient Nside_phys nth mode(full|PV) [m]"""
import numpy as np,sys,json
from lat3d_common import *
name=sys.argv[1]; Lp=float(sys.argv[2]); nth=int(sys.argv[3]); mode=sys.argv[4]; m=float(sys.argv[5]) if len(sys.argv)>5 else 0.0
o=Orient(name); Ns=int(round(Lp*o.dn))
if mode=='PV': Ms=np.sqrt(np.array([0,1,2,3]))*m; cw=np.array([1,-3,3,-1.])
else: Ms=np.array([0.0]); cw=np.array([1.])
S=S_cut(o,Ns,Ms,nth); eta=float(cw@S)
A=6*np.log(2)-3*np.log(3); extra=f" pred {A*m*m/(24*np.pi):.6e} ratio {eta/(A*m*m/(24*np.pi)):.4f}" if mode=='PV' else ''
print(f"{name} eta {mode} m={m} Nside={Ns} nth={nth}: eta={eta:.8e}{extra}",flush=True)
json.dump(dict(eta=eta,S=S.tolist()),open(f"eta_{name}_{mode}{m}_L{Lp:g}_n{nth}.json","w"))
