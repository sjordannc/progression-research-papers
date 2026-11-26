# Pre-registration, fermion test (written Tue Nov 24 02:20:20 UTC 2026, before any fermion run)
Field: free Dirac fermion in 3D, radial lattice, angular quantum number kappa = ±1, ±2, ... (degeneracy 2|kappa|).
Doubling: staggered radial lattice (upper radial component u on sites j, lower v on half-sites j+1/2) with the
"gauge-covariant" radial derivative B = e^{-Phi} d e^{Phi}, Phi' = kappa/f, so that B annihilates the continuum zero mode
and the angular gap is >= 2 sinh(kappa/2f) at all lattice momenta (no low-energy mode at k = pi).
Geometry coupling: radial operator d/dr + kappa/f(r) (spin connection absorbed by psi = chi/f); lapse as H_N = sqrt(N) H sqrt(N).
Vacuum energy -1/2 sum |E|. Same PV scheme as note 5: M^2 = (0,1,2,3) m^2, weights (1,-3,3,-1). A = 6 ln2 - 3 ln3 = 0.8630.
Covariant (heat-kernel) predictions for one Dirac fermion, relative to one minimally coupled real scalar:
  Lambda_D = -4 Lambda_S          (P0, doubling check: a doubler would give -8; the naive unstaggered lattice is run to show it)
  C_D = 2 C_S(xi=0) = -A m^2/(96 pi^2) = -9.11e-4 m^2
  eta_D = 2 eta_S = A m^2/(24 pi) = 1.145e-2 m^2      (no contact term for spin 1/2: Kabat 1995)
P1F: induced action Einstein-like: C_D width-independent after the free k^0 constant; lapse-lapse O(k^2) term absent (<10% of EH scale, falling with w).
P2F: -4 pi C_D / eta_D = 1 within 10%.  FAIL otherwise.
P4F (local source): ball modular Hamiltonian of the massless fermion equals 2 pi * (1/2){W, H} (symmetric/gravitating T00),
   with NO extra bulk term: in a fit h = a*CHM + b1*(pi/R) psi^dag psi + b2*(pi/R) psibar psi, b1 and b2 consistent with 0
   (no dimension-1 operator exists; psibar psi has dimension 3 -> R^6, negligible). FAIL if b1 or b2 is comparable to the
   scalar's (pi/R) phi^2 term, i.e. |b|/a > 0.1.
Consequence logic: P2F and P4F pass => in this model entanglement equilibrium holds for all states of the fermion with
G_ind = 1/(4 eta) (mass-dependent part). Cutoff-scale part of 1/G still not computed.
