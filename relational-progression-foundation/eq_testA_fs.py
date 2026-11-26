# Finite-size / angular-cutoff check for Test A: dS at checkpoints of lmax, for given N and n
import numpy as np, sys
exec(open('eq_testA.py').read().split("def total_S")[0].replace("N=int(sys.argv[1]) if len(sys.argv)>1 else 60","N=int(sys.argv[1])"))
n=int(sys.argv[2]); R=n+0.5; kappa=0.295; cps=[750,1500,3000,6000]
def partial(fn):
    r=np.arange(0,N+2,dtype=float); f=fn(r); fh=fn(r+0.5); S=0.0; out=[]
    for l in range(max(cps)):
        S+=(2*l+1)*S_cut(K_l(l,f,fh),n)
        if l+1 in cps: out.append(S)
    return np.array(out)
S0=partial(lambda r:r)
for eps,wid in [(0.02,3.0),(0.02,6.0)]:
    fn=lambda r: r*(1+eps*np.exp(-((r-R)/wid)**2)); dA=fn(np.array([R]))[0]**2-R**2
    d=partial(fn)-S0
    print(f"N={N} R={R} eps={eps} w={wid}: ratio at lmax {cps}: "+"  ".join(f"{x/(kappa*dA):.4f}" for x in d),flush=True)
