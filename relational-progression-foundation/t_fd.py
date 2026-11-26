# finite-difference check of the Hellmann-Feynman gradient on a few K
import numpy as np
from clock_common import *
b=Basis(); print('basis',b.nn,b.ns)
rng=np.random.default_rng(1); a=1e-3*rng.standard_normal(b.n)
Ks=[1,7,40,300]
E,g=chain_grad(Ks,a,b)
for k in [3,10,b.nn+4,b.nn+12]:
    da=np.zeros(b.n); h=1e-5; da[k]=h
    Ep,_=chain_grad(Ks,a+da,b); Em,_=chain_grad(Ks,a-da,b); print(k,g[k],(Ep-Em)/(2*h))
