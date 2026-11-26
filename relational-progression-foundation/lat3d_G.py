"""usage: lat3d_G.py orient kind(S|G) widths Lphys nth mode(full|PV) [m]
writes a2(w) (second-order coefficient of E per area) for each width; for PV mode the PV combination."""
import numpy as np,sys,json,time
from lat3d_common import *
name,kind=sys.argv[1],sys.argv[2]; W=[float(x) for x in sys.argv[3].split(',')]; Lp=float(sys.argv[4]); nth=int(sys.argv[5]); mode=sys.argv[6]
m=float(sys.argv[7]) if len(sys.argv)>7 else 0.0
o=Orient(name); N=int(round(Lp*o.dn)); eps=np.array([0,0.03,-0.03,0.06,-0.06,0.09,-0.09]); t=time.time()
if mode=='PV': Ms=np.sqrt(np.array([0,1,2,3]))*m; cw=np.array([1,-3,3,-1.])
else: Ms=np.array([0.0]); cw=np.array([1.])
E0=cw@E_vac(o,N,lambda u:0*u,kind,Ms,nth); out={}
for w in W:
    Es=np.array([E0 if e==0 else cw@E_vac(o,N,lambda u,e=e:e*np.exp(-u*u/(2*w*w)),kind,Ms,nth) for e in eps])-E0
    P6=np.polyfit(eps,Es,6); P4=np.polyfit(eps,Es,4); out[w]=dict(a1=P6[-2],a2=P6[-3],a2_4=P4[-3])
    print(f"{name} {kind} {mode} m={m} nth={nth} N={N} w={w}: a1={P6[-2]:.10f} a2={P6[-3]:.8f} (deg4 {P4[-3]:.8f}) {time.time()-t:.0f}s",flush=True)
json.dump({str(k):v for k,v in out.items()},open(f"G_{name}_{kind}_{mode}{m}_L{Lp:g}_n{nth}.json","w"),indent=1)
