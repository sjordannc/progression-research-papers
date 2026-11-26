# Pre-registration, curvature at different angles to the lattice axes (written Tue Nov 24 12:15:50 UTC 2026, before any run)
Question (Scott): test curvature at different angles to the lattice axes (the sharper test proposed in note 9).
Regulator: full 3D cubic lattice, one Wilson (r=1) Dirac fermion, continuous time. A regular 3D lattice cannot keep a
fixed PHYSICAL spacing in a curved space, so the lattice is a fixed graph and the metric enters through the link
couplings (inverse triad M = e^{-1}, kinetic link coefficient sum_i alpha^i M_ij at the link midpoint; Wilson weight
(M M^T)_jj). Variables chi = h^{1/4} psi, so the site measure is flat. This is the standard way to put a lattice field on a
curved background; it is a different fixed-lattice regulator from notes 7-8.
Geometry: planar profile along unit direction d, d in {[100], [110], [111]} (0, 45, 54.7 degrees to the nearest axis),
  (S) dl^2 = du^2 + e^{2 sigma(u)} dl_perp^2  (int sqrt(h) R = 2 int sigma'^2, the notes 7-8 geometry),
  (G) dl^2 = e^{2 alpha(u)} du^2 + dl_perp^2   (pure coordinate change: R = 0; a covariant regulator gives no two-derivative term).
Reduction: translations perpendicular to d are exact lattice symmetries; Bloch momenta for them leave 1D block chains.
Extraction: E(eps) for sigma = eps g, g Gaussian of width w; second-order coefficient a2(w) = A Vol2 + Q I2 + D I4
  (A, the zero-derivative term, is fitted, not fixed, because a coordinate lattice gives a non-sqrt(h) volume term).
  C(d) = Q_S(d)/2, rho(d) = -4 pi C(d)/eta(d), eta(d) = entanglement per area across a lattice plane normal to d.
Caveat stated in advance: with a full cubic triad, the Einstein term is itself expressible through frame (torsion) terms
(teleparallel identity), so C_R is not separately defined for this regulator. The test therefore asks whether frame terms
are present and how large, not for C_R itself.
Pre-registered readings:
  FRAME TERMS PRESENT (supports note 9): rho(d) differs between orientations by more than 3x its numerical error, or the
    pure-coordinate-change coefficient |Q_G| exceeds 5% of |Q_S| (at the same orientation).
  NO FRAME TERMS (contradicts note 9's account): rho(d) orientation-independent within errors AND |Q_G| < 5% |Q_S|.
  Separately: PASS of equilibrium for this regulator if rho(d) within 10% of 1 for all three d.
Validation (must hold or the run is void): flat vacuum energy density equal for all d (to 1e-6 relative); Pauli-Villars
(m = 0.2, masses sqrt(0,1,2,3) m): rho_PV(d) within 10% of 1 and |Q_G,PV| < 10% of |Q_S,PV| for each d.
