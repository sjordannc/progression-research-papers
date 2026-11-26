# Pre-registration, two-clock test (written Wed Nov 25 18:48:44 UTC 2026, before any run)
Question (Scott): do the two-clock test, i.e. P6: recover d tau_A / d tau_B = 1 + (Phi_A - Phi_B)/c^2 and the ruler half
(sigma = -theta) from clocks plus an environment, without inserting Phi.
Model (the draft's induced-gravity sector, notes 5-6): classical static geometry (lapse N(r), spatial metric
dl^2 = dr^2 + f(r)^2 dOmega^2, r = proper radial distance) with NO gravitational action. Its only energy is the vacuum
(Killing) energy of the environment, one free Dirac fermion on the note-6 staggered radial lattice, Pauli-Villars
combined (m = 0.2; the part for which notes 6-7 verified -4 pi C/eta = 1), plus a static dust shell of rest energy
E_src = int N mu dr (universal coupling, P3; no Newtonian coupling, no potential-dependent weight: R1).
Field equations = stationarity of E_tot[N,f] (dE/dN: local energy constraint; dE/df: static lapse equation).
Counterterms, fixed in advance from flat-space measurements in note 6 (no source involved): the vacuum term
Lambda1 * int N sqrt(h) (Lambda1 = 1.560709e-5) and the lattice k^0 artifact (L2-L1) * 4 pi int r^2 s^2,
(L2-L1)/L1 = -2.47e-4, with f = r(1+s). (Cosmological-constant tuning is not addressed.)
Solve: N = e^nu, nu and s on cubic B-splines (spacing 10; nu on [0,270], s on [70,270]); stationarity imposed for test
functions supported in r <= 220; lattice gradients by Hellmann-Feynman; linear response by symmetric differences;
defect-correction iteration preconditioned by the Einstein functional (converged when the lattice residual < 1e-3 of the
source). Shell: Gaussian, radius 100, width 8, weak field (|nu| ~ 1e-3).
Prediction (no fit): the same discretised equations with the Einstein functional C int N sqrt(h) R, C = -1/(16 pi G),
G = 1/(4 eta_D), eta_D = 4.6404e-4 = the fermion's own entanglement coefficient measured in note 6.
Clocks: A1 = 30, A2 = 60 (inside the shell), B1 = 160, B2 = 210 (outside). Clock phase rates relative to the reference
are N(x) * Delta E (P3); checked for a second clock mechanism (a localized level of a separate Dirac chain).
PASS (all):
 (P1) outside: [nu(B1)-nu(B2)]_lattice / [same]_prediction within 10% of 1;
 (P2) inside vs outside: [nu(A2)-nu(B2)] lattice/prediction within 10% of 1;
 (P3) inside the shell: |nu(A1)-nu(A2)| < 10% of |nu(A2)-nu(B2)| (no potential difference inside a shell);
 (P4) ruler half: s(B1), s(B2) lattice/prediction within 15% of 1;
 (P5) the inside-outside clock ratio matches 1 + Phi_A - Phi_B with Phi the Newtonian potential of the shell for
      G = 1/(4 eta_D), within 10%.
FAIL if any fails; PARTIAL if P1-P3 pass and P4 fails.
Scope stated in advance: tests the regulator-independent (Pauli-Villars) part only; the cutoff part of G failed on all
fixed lattices (notes 7, 8, 10) and is assumed cancelled. Weak field, static, spherical, one species.

## Amendment (written Wed Nov 25 18:49:38 UTC 2026, before any field-equation run; only the code check t_fd.py had been run)
Shell width 8 -> 12: the Pauli-Villars Compton length is 1/m = 5, so a width-8 shell would put higher-derivative terms
at the (1/(m w))^2 ~ 40% level near the shell. With w = 12 they are ~17% near the shell and decay away from it; the clock
criteria use points >= 40 from the shell.

## Amendment 2 (written Wed Nov 25 18:53:06 UTC 2026, before any lattice field-equation run; only GR-discretisation tests had been run)
A finite-box Galerkin solve is ill-posed for an isolated positive mass: the spatial metric outside the shell cannot
vanish at a finite boundary (GR test: inside-shell lapse not constant, rank-deficient rows). Replaced by an equivalent,
box-independent test. The exact linearised solution in this gauge is (derived from the same Einstein functional):
  nu'(r) = G m(r)/r^2  (Newton), s(r) = -(G/r) int_0^r m(r')/r' dr',  m(r) = int_0^r mu.
Test: insert nu = lam_nu * nu_GR, s = lam_s * s_GR (G = 1/(4 eta_D)) into the LATTICE stationarity equations, tested
against B-splines supported in [0,220] (nu rows: local energy constraint; s rows: lapse equation). Fit lam_nu, lam_s by
least squares; the lattice clock ratio is then 1 + lam_nu (Phi_A - Phi_B).
PASS (all): (Q1) lam_nu within 10% of 1 (clock half: gives P1, P2, P5); (Q2) lam_s within 15% of 1 (ruler half, P4);
(Q3) residual after the fit < 10% of the source projection, overall and on the rows inside the shell (r<80), which
tests the inside-shell constancy P3 and the 1/r shape. Same counterterms as before. Without the k^0 artifact counterterm
also reported.
