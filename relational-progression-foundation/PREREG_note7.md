# Pre-registration, cutoff-scale part of G (written Tue Nov 24 03:29:18 UTC 2026, before any run)
Goal: compare the FULL (cutoff-scale, no Pauli-Villars) induced 1/G with the FULL entanglement area coefficient eta,
same regulator on both sides, for one Dirac fermion.
Why a new geometry: on the radial lattice the angular directions are a continuum, so the cutoff-scale vacuum energy is
unregulated. We use planar geometry dl^2 = dx^2 + e^{2 sigma(x)} (dy^2+dz^2):
  x: staggered lattice, physical spacing a=1 (fixed, since g_xx=1);
  transverse: plane waves, PHYSICAL momentum q = k e^{-sigma(x)} cut sharply at q < Lambda_perp: a chain site at x
  exists for coordinate momentum k only if k e^{-sigma(x)} < Lambda_perp. So the number of modes per physical volume is
  the same everywhere: the cutoff is physical in all three directions (time is continuous).
Dirac reduction: psi = e^{-sigma} chi gives 1D Dirac chains with position-dependent mass mu(x) = k e^{-sigma(x)}, two per k.
Vacuum energy per unit coordinate transverse area: E = -(1/2pi) int k dk sum|E_n(k)| (k-integral done piecewise between
the moving cutoff breakpoints, so E(eps) is smooth).
Expansion: E = Lam*int e^{2sigma} dx + C*int sqrt(h) R + ... ; int sqrt(h) R = 2 int sigma'^2 (2nd order). C = -1/(16 pi G_ind).
Entanglement: eta = (1/2pi) int_0^Lambda k dk * 2 * S_1D(k), plane cut in the flat chain, same cutoff.
Covariant expectation (Susskind-Uglum/Jacobson): -4 pi C = eta. Our regulator is NOT covariant (lattice in x, sharp disk
in the transverse plane, continuous time), so the outcome is not guaranteed.
Runs: Lambda_perp in {2, pi, 2 sqrt(pi)} (2 = max lattice frequency in x; pi = zone edge; 2 sqrt(pi) = isotropic mode density).
PASS: rho = -4 pi C / eta within 10% of 1 for all three Lambda_perp.
PARTIAL: within 10% for some but not all (=> the cutoff part matches only for a tuned regulator; G regulator-dependent).
FAIL: no Lambda_perp within 10%.
Also: k^0 check (second-order volume coefficient = first-order Lambda within 1e-3) and width/box/precision checks.
Reported: a/l_P = sqrt(16 pi N_D |C|) for each regulator (G traded for a: a consistency statement, not a prediction).

## Amendment (written Tue Nov 24 03:33:15 UTC 2026, after the sharp-cutoff runs, before any smooth-cutoff run)
Sharp-cutoff result: C(w) grows linearly with bump width (35 -> 140 for w = 6 -> 24 at Lambda_perp=2), i.e. the
second-order energy has a w-independent term: the cutoff surface in (x, k) space costs energy, giving a non-analytic
int|grad sigma| term. With a sharp physical cutoff the induced action is NOT Einstein-like, so G is undefined for it.
Smooth physical cutoff (second regulator): transverse dispersion mu(q) = q exp(q^2/Lambda^2), q = k e^{-sigma(x)};
modes removed only where mu > mu_cut (q_cut = 3.7 Lambda, mu_cut ~ 1e6), where the removal edge is suppressed by ~1/mu_cut.
The purely diagonal (sum over sites of mu) part of the vacuum energy is exactly a lattice volume term and is subtracted
analytically for precision. Same pass/fail conditions as above (rho within 10% of 1), for Lambda in {1, 1.5, 2};
additional checks: insensitivity to q_cut (3.4 vs 3.7 Lambda) and to box length.
