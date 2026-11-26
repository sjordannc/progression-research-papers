import numpy as np,sys
import fermion_common as fc
fc.JCUT=float(sys.argv[3])
m=0.2; w=float(sys.argv[1]); Nr=360; r0=150.; eps=float(sys.argv[2]); fl=lambda r:r
P={}
for s in (0,1,-1):
    P[s]=fc.E0_PV(fl,lambda r,s=s:1+s*eps*np.exp(-(r-r0)**2/(2*w*w)),Nr,m,1600)
q=(P[1]+P[-1]-2*P[0])/(2*eps**2)
print("JCUT",fc.JCUT,"eps",eps,[(L,round(q[:L].sum(),6)) for L in [100,200,400,600,800,1000,1200,1600]],flush=True)
