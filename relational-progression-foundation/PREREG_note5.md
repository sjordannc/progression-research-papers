# Pre-registration (written Tue Nov 24 00:21:21 UTC 2026 before any run in this step)
Model: semiclassical induced gravity (Sakharov) on the 3D radial lattice. The geometry (spatial metric f(r), lapse N(r))
has no bare action; its dynamics is the field's own vacuum (Killing) energy E0[N,f]. Equilibrium test: is the geometry's
induced Einstein coefficient C (E0 ⊃ C ∫√h R, C = -1/(16πG_ind)) equal to -η/(4π), η the field's entanglement area coefficient?
Regulator-independent part: Pauli–Villars combination over masses M^2 = (0,1,2,3)m^2, weights (1,-3,3,-1)
(a third finite difference in M^2 of physical quantities).
Analytic covariant predictions (m=lattice-unit mass scale; A = 6ln2-3ln3 = 0.8630):
  eta_PV  = A m^2/(48π) = 5.72e-3 m^2
  C_PV    = -eta_PV/(4π) = -4.55e-4 m^2   (xi=0)
  <phi^2>_PV = 3 eta_PV/π = 5.47e-3 m^2
  C_PV(xi=1/6) = C_PV(0) + <phi^2>/12 = 0
P1 (geometry dynamics is Einstein-like): C extracted from f-deformations is width-independent after removing
   the 4-derivative tail; lapse-lapse Hessian has no O(k^2) term at the level of 10% of the EH scale.
P2 (equilibrium = dynamics, xi=0): -4π C_PV / eta_PV = 1 within 10%.  FAIL if outside 10% (after lattice/finite-size errors).
P3 (xi): C_PV(1/6) consistent with 0 => plain entanglement equilibrium fails for conformal coupling (contact term).
B' (edge term): for bumps reaching the entangling sphere, Q_exact = Q_min + (π/R)Σv^2 - (π/3) u(R)^2 (improved, with
   its boundary term). FAIL if the boundary term is absent or of different size.
Consequence logic: if P2 and B' pass, then at xi=0 the induced geometry is sourced by T^min while the entropy balance
requires T^imp -> equilibrium FAILS locally by the improvement tensor (∂∂-g∂²)<φ²>/6.
