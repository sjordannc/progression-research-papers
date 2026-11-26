import numpy as np,sys,json
from smooth_common import *
Lam=float(sys.argv[1]); L=int(sys.argv[2]); qf=float(sys.argv[3]) if len(sys.argv)>3 else 3.7
eta=S_plane(L,Lam,qf*Lam); print(f"Lam={Lam} L={L} qf={qf}: eta={eta:.8f}",flush=True)
json.dump(dict(eta=eta),open(f"seta_L{Lam}_N{L}_q{qf}.json","w"))
