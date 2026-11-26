import numpy as np
# paraxial: dpsi/dz = (i/2k) d2psi/dx2 + i k dn(x) psi
lam=1.0; k=2*np.pi/lam; L=400.0; N=4096
x=np.linspace(-L/2,L/2,N,endpoint=False); dx=x[1]-x[0]
kx=2*np.pi*np.fft.fftfreq(N,dx)
w=10.0; Z=2000.0; nz=2000; dz=Z/nz
def run(psi0,dn):
    psi=psi0.astype(complex); cents=[]
    lin=np.exp(-1j*kx**2/(2*k)*dz/2)
    for s in range(nz+1):
        I=abs(psi)**2; cents.append((x*I).sum()/I.sum())
        if s==nz: break
        psi=np.fft.ifft(lin*np.fft.fft(psi)); psi*=np.exp(1j*k*dn*dz); psi=np.fft.ifft(lin*np.fft.fft(psi))
    return np.array(cents)
g=np.exp(-x**2/(2*w**2))
G=1e-6  # index gradient per unit length
cases={
 "symmetric beam, uniform medium":(g,0*x),
 "amplitude skewed toward -x, uniform medium":(g*(1+0.8*np.tanh(-x/w)),0*x),
 "symmetric beam, index gradient (higher n at -x)":(g,-G*x),
 "skewed beam, index gradient":(g*(1+0.8*np.tanh(-x/w)),-G*x),
}
pred=-0.5*G*Z**2
for name,(p,dn) in cases.items():
    c=run(p,dn); print(f"{name:50s} centroid shift over z: {c[-1]-c[0]:+.4f}")
print("ray-optics prediction for gradient case:",pred)
