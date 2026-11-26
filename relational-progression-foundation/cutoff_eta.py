import numpy as np,sys,json
from cutoff_common import *
Lam=float(sys.argv[1]); L=int(sys.argv[2])
eta=S_plane(L,Lam); print(f"Lam_perp={Lam:.4f} L={L}: eta={eta:.6f}",flush=True)
json.dump(dict(eta=eta),open(f"cuteta_L{Lam:.3f}_N{L}.json","w"))
