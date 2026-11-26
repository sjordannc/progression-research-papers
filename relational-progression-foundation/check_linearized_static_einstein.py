import sympy as sp
t,x,y,z,e=sp.symbols('t x y z epsilon')
X=[t,x,y,z]
A=sp.Function('A')(x,y,z); B=sp.Function('B')(x,y,z)
g=sp.diag(-(1+2*e*A),1+2*e*B,1+2*e*B,1+2*e*B)  # c=1, static
gi=g.inv()
def Gam(a,b,c): return sum(gi[a,d]*(sp.diff(g[d,b],X[c])+sp.diff(g[d,c],X[b])-sp.diff(g[b,c],X[d])) for d in range(4))/2
G=[[[sp.series(Gam(a,b,c),e,0,2).removeO() for c in range(4)] for b in range(4)] for a in range(4)]
def Ric(b,c):
    r=0
    for a in range(4):
        r+=sp.diff(G[a][b][c],X[a])-sp.diff(G[a][b][a],X[c])
        for d in range(4): r+=G[a][a][d]*G[d][b][c]-G[a][c][d]*G[d][b][a]
    return sp.expand(sp.series(r,e,0,2).removeO()).coeff(e)
R=sp.Matrix(4,4,lambda i,j:Ric(i,j))
eta=sp.diag(-1,1,1,1)
Rs=sum(eta[i,i]*R[i,i] for i in range(4))
Gt=sp.simplify(R-eta*Rs/2)
lap=lambda f: sum(sp.diff(f,v,2) for v in (x,y,z))
S=A+B
print("G00 + 2 lap B =",sp.simplify(Gt[0,0]+2*lap(B)))
for (i,j) in [(1,1),(1,2),(2,3)]:
    expect=(lap(S) if i==j else 0)-sp.diff(S,X[i],X[j])
    print(f"G{i}{j} - (delta lap - d d)(A+B) =",sp.simplify(Gt[i,j]-expect))
