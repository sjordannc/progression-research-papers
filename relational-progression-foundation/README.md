# Relational Progression as an Operational Foundation for General Relativity

Consolidated draft (November 2026) with companion notes 2–12, an earlier unnumbered note, their scripts, pre-registrations and results. The draft and notes are working research documents, not peer-reviewed results. Each claim in the draft carries a status tag (Standard, Derived, Conditional, Open, Fails), summarised in the draft's status ledger.

**Draft:** [PDF](relational_progression_foundation.pdf) · [LaTeX](relational_progression_foundation.tex)

## Companion notes

| Note | Topic | Outcome | Files |
| --- | --- | --- | --- |
| Earlier note | Energy imbalance as the mechanism of free fall and refraction | Analysis (not cited as a numbered note) | [PDF](energy_imbalance_mechanism_note.pdf) |
| 2 | Making the hypersurface-deformation route rigorous | Assumptions traced to postulates | [PDF](hda_route_rigor_note.pdf) |
| 3 | Lattice toy model of the Clausius relation | Derived (numerical toy model) | [PDF](clausius_toy_model_note.pdf) |
| 4 | Inputs of entanglement equilibrium on a 3D lattice | Inputs derived (numerical) | [PDF](entanglement_equilibrium_lattice_note.pdf) |
| 5 | Entanglement equilibrium with a responding geometry (scalar) | Conditional: G = 1/(4η) for the regulator-independent part; fails for φ² states at minimal coupling | [PDF](induced_gravity_equilibrium_note.pdf) |
| 6 | Fermion environment | Derived (numerical): equilibrium holds, G = 1/(4η) to 3% | [PDF](fermion_equilibrium_note.pdf) |
| 7 | Cutoff-scale part of G, physical lattice cutoff | Fails (wrong sign) | [PDF](cutoff_scale_G_note.pdf) |
| 8 | Cutoff treating space and time alike | Fails (wrong sign, smaller) | [PDF](euclid_cutoff_G_note.pdf) |
| 9 | Why fixed lattices fail and covariant cutoffs pass | Analysis | [PDF](why_covariant_note.pdf) |
| 10 | Curvature at different angles to the lattice axes | Frame terms present; cutoff part fails at every angle | [PDF](angles_note.pdf) |
| 11 | Can the pseudotensor close G? | No (analysis) | [PDF](pseudotensor_note.pdf) |
| 12 | Two-clock test: weak-field clock ratio from clocks plus an environment | Derived (numerical), conditional; pre-registered run failed narrowly, passes after a convergence fix | [PDF](two_clock_note.pdf) |

The LaTeX source of each note sits beside its PDF. Pre-registrations are in `PREREG_note*.md` and raw outputs (logs, JSON) in `note*_results/`; runs that were discarded are kept under `void/` where applicable.

## Open problem

The cutoff-scale part of G is not derived: every fixed-spacing lattice regulator tried gives it the wrong sign (notes 7, 8, 10). The candidates left are frame-free discreteness (for example random sprinkling) or a species sum rule.

## Running the scripts

Python 3 with NumPy, SciPy and mpmath. Set `OMP_NUM_THREADS=1` when running several jobs in parallel. Each note lists the scripts it uses and their arguments.
