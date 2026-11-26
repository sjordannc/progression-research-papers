"""Masses are twisted (M X5, X5 = alpha1 alpha2 alpha3 beta), chirally equivalent to a Dirac mass in the continuum, so that H^2 = H0^2 + M^2 and Pauli-Villars cancels on the lattice too.
Full 3D cubic Wilson-Dirac lattice (r=1), continuous time, in a planar metric profile along an arbitrary lattice direction d.
Metric enters through the inverse triad M on links (kinetic coefficient A_j = sum_i alpha^i M_ij at the link midpoint,
Wilson weight (M M^T)_jj), in variables chi = h^{1/4} psi with flat site measure.
Translations perpendicular to d are exact symmetries (basis t1, t2); Bloch phases theta1, theta2 reduce the problem to a
block-tridiagonal chain in the layer index n = d_int . x (layer spacing 1/|d_int| along d, T-cell area |d_int|)."""
import numpy as np,os
from scipy.linalg import eig_banded, eigh
s1=np.array([[0,1],[1,0]],complex); s2=np.array([[0,-1j],[1j,0]]); s3=np.array([[1,0],[0,-1]],complex); Z=np.zeros((2,2))
ALPHA=[np.block([[Z,s],[s,Z]]) for s in (s1,s2,s3)]; BETA=np.diag([1.,1.,-1.,-1.]).astype(complex)
X5=ALPHA[0]@ALPHA[1]@ALPHA[2]@BETA  # Hermitian, X5^2=1, anticommutes with all alpha and beta: twisted mass, H^2=H0^2+M^2
ORIENT={'100':((1,0,0),(0,1,0),(0,0,1)),'110':((1,1,0),(1,-1,0),(0,0,1)),'111':((1,1,1),(1,-1,0),(0,1,-1))}
class Orient:
    def __init__(self,name):
        d,t1,t2=[np.array(v,float) for v in ORIENT[name]]; self.dint=d.astype(int); self.dn=np.linalg.norm(d); self.dhat=d/self.dn
        self.T=np.array([t1,t2]); self.area=np.linalg.norm(np.cross(t1,t2)); self.Ms=np.vstack([t1,t2,d])
    def k_of(self,th1,th2): return np.linalg.solve(self.Ms,np.array([th1,th2,0.0]))
def Mtriad(dhat,val,kind):
    P=np.outer(dhat,dhat); return (P+np.exp(-val)*(np.eye(3)-P)) if kind=='S' else (np.exp(-val)*P+(np.eye(3)-P))
def link_data(dhat,val,kind):
    M=Mtriad(dhat,val,kind); A=[sum(ALPHA[i]*M[i,j] for i in range(3)) for j in range(3)]; w=np.sum(M*M,axis=1)  # (M M^T)_jj
    return A,w
def blocks(o,k,prof,kind,m,r=1.0):
    """prof(u) metric function; returns diag blocks D[n], upper blocks U[n] (n->n+1)."""
    N=prof.N; D=np.zeros((N,4,4),complex); U=np.zeros((N-1,4,4),complex)
    for n in range(N):
        A,w=link_data(o.dhat,prof.at_layer[n],kind); D[n]+=m*BETA
        for j in range(3):
            if o.dint[j]==0: D[n]+=A[j]*np.sin(k[j])+0.5*r*w[j]*BETA*(2-2*np.cos(k[j]))
    for n in range(N-1):
        A,w=link_data(o.dhat,prof.at_mid[n],kind)
        for j in range(3):
            if o.dint[j]==1:
                ph=np.exp(1j*k[j]); U[n]+=(-0.5j)*A[j]*ph-0.5*r*w[j]*BETA*ph; D[n]+=0.5*r*w[j]*BETA; D[n+1]+=0.5*r*w[j]*BETA
    return D,U
def band_eigs(D,U):
    N=D.shape[0]; n=4*N; ab=np.zeros((8,n),complex)
    H_lo=lambda i,j: None
    for b in range(N):
        for a in range(4):
            for c in range(a,4):  # lower part of diag block: row b*4+c, col b*4+a
                ab[c-a,b*4+a]=D[b][c,a]
    for b in range(N-1):
        L=U[b].conj().T  # block (b+1,b)
        for a in range(4):
            for c in range(4):
                ab[4+c-a,b*4+a]=L[c,a]
    return eig_banded(ab,lower=True,eigvals_only=True,check_finite=False)
class Profile:
    def __init__(self,o,N,fn):
        self.N=N; u=(np.arange(N)-N/2)/o.dn; self.at_layer=fn(u); self.at_mid=fn(u+0.5/o.dn)
def theta_grid(nth):
    """midpoint grid on [-pi,pi)^2, keeping half by theta -> -theta symmetry (weights doubled)"""
    t=-np.pi+(np.arange(nth)+0.5)*2*np.pi/nth; T1,T2=np.meshgrid(t,t,indexing='ij'); T1=T1.ravel(); T2=T2.ravel()
    keep=(T1>0)|((T1==0)&(T2>0)) if False else (T1>0); return T1[keep],T2[keep],2.0/(nth*nth)
def theta_grid_refined(nth,ncore=2,levels=7,ng=3):
    """coarse midpoint grid (half, by symmetry) outside the central (2 ncore)^2 coarse cells; inside, nested square rings
    halving toward theta=0 with ng x ng Gauss points per sub-square (resolves the Dirac point and the PV mass scale)."""
    h=2*np.pi/nth; t=-np.pi+(np.arange(nth)+0.5)*h; T1,T2=np.meshgrid(t,t,indexing='ij'); T1=T1.ravel(); T2=T2.ravel()
    K=ncore*h; out=(np.abs(T1)>K)|(np.abs(T2)>K); keep=out&(T1>0)
    P1=list(T1[keep]); P2=list(T2[keep]); Wt=[2.0/(nth*nth)]*int(keep.sum())
    x,wg=np.polynomial.legendre.leggauss(ng)
    def square(cx,cy,half,f):
        for a,wa in zip(x,wg):
            for b,wb in zip(x,wg):
                P1.append(cx+half*a); P2.append(cy+half*b); Wt.append(f*wa*wb*half*half/(4*np.pi**2))
    Kj=K
    for lev in range(levels):
        q=Kj/2  # ring between Kj and Kj/2: 4x4 squares of side q, minus central 2x2
        for i in range(4):
            for j in range(4):
                if 1<=i<=2 and 1<=j<=2: continue
                cx=-Kj+(i+0.5)*q; cy=-Kj+(j+0.5)*q
                if cx>0: square(cx,cy,q/2,2.0)  # use theta->-theta symmetry: keep cx>0 squares, double
        Kj=q
    for i in range(2):
        for j in range(2):
            cx=-Kj+(i+0.5)*Kj; cy=-Kj+(j+0.5)*Kj
            if cx>0: square(cx,cy,Kj/2,2.0)
    return np.array(P1),np.array(P2),np.array(Wt)
GRID=os.environ.get('GRID','refined')
def get_grid(nth):
    if GRID=='uniform':
        T1,T2,w=theta_grid(nth); return T1,T2,np.full(len(T1),w)
    return theta_grid_refined(nth)
def E_vac(o,N,fn,kind,masses,nth):
    """vacuum energy per unit physical transverse area, -1/2 sum|E|, for each mass. returns array over masses"""
    prof=Profile(o,N,fn); T1,T2,WW=get_grid(nth); out=np.zeros(len(masses))
    for a,b,wgt in zip(T1,T2,WW):
        k=o.k_of(a,b)
        D,U=blocks(o,k,prof,kind,0.0); e2=band_eigs(D,U)**2
        for i,m in enumerate(masses): out[i]+=wgt*np.sum(np.sqrt(e2+m*m))
    return -0.5*out/o.area
def S_cut(o,Nside,masses,nth):
    """entanglement per unit physical area across the lattice plane between layers Nside-1 and Nside (flat), each mass"""
    prof=Profile(o,2*Nside,lambda u:0*u); T1,T2,WW=get_grid(nth); out=np.zeros(len(masses))
    for a,b,wgt in zip(T1,T2,WW):
        k=o.k_of(a,b)
        for i,m in enumerate(masses):
            D,U=blocks(o,k,prof,'S',0.0); N=2*Nside; H=np.zeros((4*N,4*N),complex)
            for n in range(N): H[4*n:4*n+4,4*n:4*n+4]=D[n]+m*X5
            for n in range(N-1): H[4*n:4*n+4,4*n+4:4*n+8]=U[n]; H[4*n+4:4*n+8,4*n:4*n+4]=U[n].conj().T
            e,v=eigh(H); V=v[:,e<0][:4*Nside]; c=np.linalg.eigvalsh(V@V.conj().T); c=c[(c>1e-14)&(c<1-1e-14)]
            out[i]+=wgt*np.sum(-c*np.log(c)-(1-c)*np.log(1-c))
    return out/o.area
