# Pre-registration, cutoff-scale G with a space-time symmetric cutoff (written Tue Nov 24 04:36:30 UTC 2026, before any run)
Question (Scott): repeat the note-7 test with a cutoff that treats space and time the same way.
What matters for the horizon argument (Susskind-Uglum, Jacobson) is the boost plane (t, x); in Euclidean signature,
rotations of the (tau, x) plane. Note 7 discretised x (staggered, cell a=1, dispersion p_hat = 2 sin(p/2)) but kept time
continuous (omega unbounded). Here Euclidean time is a lattice with the SAME spacing a_t = a = 1 and a Wilson (r=1) time
derivative, so the free Dirac operator obeys D^dagger D = omega_hat^2 + p_hat^2 + mu(q)^2 with omega_hat = 2 sin(omega/2):
the (tau, x) plane is a symmetric square lattice (no time doublers). Transverse directions keep the note-7 smooth physical
cutoff mu(q) = q exp((q/Lam)^p), p=2 (transverse rotations are exact; they are not the boost).
Vacuum energy (static background): E0 = -(1/2) sum_n eps(E_n), eps(E) = (2/a_t) arcsinh(a_t |E| / 2)
  (from int_{-pi}^{pi} d omega/2pi ln(omega_hat^2 + E^2) = 2 arcsinh(|E|/2); the transfer-matrix energy, cosh eps = 1 + E^2/2).
Entanglement: the propagator D^{-1} = D^dagger/(omega_hat^2 + K^2) depends on the spatial Dirac operator K only through
functions of K, so the transfer-matrix vacuum is the same Dirac sea as before; eta (same regulator) equals the note-7 eta.
Extraction, fits, geometry, widths, boxes: identical to note 7.
Test: rho = -4 pi C / eta at a_t = 1, Lam in {1, 1.5, 2}.
PASS: rho within 10% of 1 for all three. PARTIAL: for some. FAIL: none. Also report the sign of C.
Diagnostics (not part of pass/fail): a_t in {0.5, 0.25} (a_t -> 0 must recover note 7), steeper cutoff p=4 at Lam=1.5,
Pauli-Villars validation of the new code (m = 0.2: expect C/C_pred near 1, since PV parts are UV-insensitive),
width (w = 6..16) and box (L = 300 vs 600) stability.
Reported: a/l_P from 1/(4G) = N eta and, if C has the right sign, from C directly.
