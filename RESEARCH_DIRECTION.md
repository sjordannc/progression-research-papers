# Relational Progression and Weighted Conservation

Scott A. Jordan · Research direction · 8 October 2026

This page presents the current research direction and connects it to the [manuscript collection](README.md). It is a web adaptation of the research-direction note dated 8 October 2026.

## The question

I propose treating progression as a relational property of physical evolution whose comparative accumulation is constrained by conservation. Time is the operational measure of those relationships. The question is whether energy expressed relative to a local clock can be related to a conserved description through a weight derived from physical dynamics.

An energy measurement can be accurate in its local clock convention while requiring a specified conversion when compared with quantities defined using another reference. The physical task is to derive that conversion and identify what is conserved. The formulation leaves open whether a progression field is useful as an effective description; it does not require such a field as a fundamental starting point.

## Clocks and comparative accumulation

A clock comparison must specify preparation, calibration, histories, and endpoints or a signal protocol. Reunited clocks allow direct comparison of accumulated readings. Separated clocks require an operational pairing of events. No single universal progression value is assumed for all histories.

For two calibrated clocks, let θ denote accumulated unwrapped relative phase and ω₀ the angular frequency under the stated reference calibration. A comparison is

$$
P_{AB}=\frac{\Delta\theta_A/\omega_{A0}}{\Delta\theta_B/\omega_{B0}}.
$$

For identical clocks the calibration factors cancel. For different mechanisms they are essential. A single energy eigenstate's overall phase is not a clock reading; relative phases, interference, transitions, or other measurable correlations are required.

Standard quantum mechanics supplies the benchmark

$$
\Delta\theta_C=\frac{1}{\hbar}\int \Delta E_C\,d\tau_C,
$$

with phase orientation chosen to make the accumulation positive. Because this expression already uses proper time, it motivates an operational comparison but does not derive emergent time.

## Energy relative to a clock

For one phase history described using a local clock variable τ_C and a differentiable reference variable T, the chain rule gives

$$
E_T=\hbar\frac{d\theta}{dT}
=\frac{d\tau_C}{dT}\,\hbar\frac{d\theta}{d\tau_C}
=w_C E_C.
$$

Here E is the energy gap or phase generator relevant to the history. T is a specified reference, not a universal background time. At this stage w_C is a conversion factor, not an independently propagating field.

This relation does not establish conservation of E_T. A conserved quantity requires an appropriate symmetry or a complete balance law, including the interaction, environment, and physical reference system. Relabeling a clock changes no observation. A new result would either derive the relation from more primitive dynamics or predict a measurable departure.

## Weighted conservation

A candidate local formulation distinguishes a measured current j from the current entering a conservation law:

$$
\mathcal J^\mu=w j^\mu,
\qquad \partial_\mu\mathcal J^\mu=0.
$$

For positive differentiable w, this entails

$$
\partial_\mu j^\mu=-j^\mu\partial_\mu\ln w.
$$

The identity explains how an apparent source in one description can coexist with weighted conservation. It does not identify j as energy, determine w, or establish a gravitational mechanism. If w=Π^q, the exponent and the physical meaning of Π must follow from the model rather than be chosen after observing the current.

The charge also requires an integration measure and boundary conditions. Once a metric is supplied, the covariant expression is

$$
\nabla_\mu\mathcal J^\mu
=\frac{1}{\sqrt{-g}}\partial_\mu\!\left(\sqrt{-g}\,\mathcal J^\mu\right)=0.
$$

Before deriving a metric, the model must supply a suitable measure or vector density directly. A vector-current conservation law is not the same object as stress-energy balance. In a theory with extra sectors, all exchanges must be included in the total balance. The existing matter and optical papers accordingly distinguish weighted subsystem transport from conservation of total stress-energy.

## The stationary gravitational benchmark

Established general relativity already distinguishes locally measured energy from energy associated with a stationary symmetry. For a static metric,

$$
ds^2=-N(\mathbf x)^2c^2dT^2+h_{ij}dx^i dx^j,
\qquad d\tau=N\,dT
$$

for observers at fixed spatial coordinates. With the stationary Killing vector normalized relative to T, a freely propagating particle or photon has

$$
\mathcal E=N E_{\rm loc}
$$

in units c=1. The Killing energy is conserved along a geodesic; local energy can differ between locations. For a photon exchanged between static observers,

$$
\frac{E_{{\rm loc},B}}{E_{{\rm loc},A}}=\frac{N_A}{N_B},
\qquad
\frac{d\tau_A}{d\tau_B}=\frac{N_A}{N_B},
$$

where clock intervals are paired by the same stationary coordinate interval. The first relation concerns a transported signal. Locally identical stationary clocks retain their calibrated local transition gaps.

At weak field,

$$
\frac{d\tau_A}{d\tau_B}\simeq 1+\frac{\Phi_A-\Phi_B}{c^2}.
$$

These are targets to recover. Setting w equal to a known lapse reproduces the benchmark by construction. A relational derivation would calculate the observable ratio first and identify its geometric representation afterward. The stationary construction does not provide a global conserved energy in every dynamical or cosmological spacetime. See [Carroll's general-relativity notes](https://arxiv.org/abs/gr-qc/9712019).

## A working mathematical setting

A reparameterization-invariant constrained action can include physical clocks and their interactions without assigning physical status to the history label λ:

$$
S=\int d\lambda\left[\sum_a p_a\frac{dq^a}{d\lambda}-n(\lambda)\mathcal K(q,p)\right],
\qquad\mathcal K=0.
$$

Where a clock observable C is monotonic and its bracket with the constraint is nonzero,

$$
\frac{dO}{dC}=\frac{\{O,\mathcal K\}}{\{C,\mathcal K\}}.
$$

The arbitrary multiplier cancels. This removes the privilege of the history label; it does not independently derive energy conservation or gravity. Internal-clock quantum descriptions have an established precedent in [Page and Wootters](https://doi.org/10.1103/PhysRevD.27.2885).

A schematic candidate is

$$
\mathcal K=\sum_i w_i(q)H_i(q_i,p_i)+H_{\rm int}(q,p)-Q_0.
$$

This is an ansatz, not a completed physical model. It assumes a setting where the candidate fixed charge Q₀ is meaningful. The weights and interaction term must arise from a specified parent model.

## How the existing results guide the next step

The [connected Jacobi model](relational-clock-transport/manuscript.pdf) provides a finite example of transported phase clocks, conserved phase charges, and relational observables. Its stable fast-mode construction induces clock inertias rather than putting their organization dependence directly into the bare kinetic coefficients. The mechanism remains model dependent and does not derive gravitational redshift.

Its cancellation result is particularly important: a common fractional response is absent from local pairwise clock comparisons. Testing universality therefore requires additional operational structure. Unequal responses of model clocks establish a differential effect, not universal temporal geometry.

The [HDA manuscript](progression_weighted_conservation.pdf) supplies a complementary conditional route toward continuum geometry. Its downstream canonical and matter calculations help identify the structure a more microscopic theory must reproduce. Its local-evolution, integrability, locality, and calibration assumptions remain part of the foundational problem.

The field-based papers explore concrete effective completions. Their results remain useful within their assumptions even when the current foundational formulation does not postulate a progression field. The [overview](README.md#what-has-been-obtained-and-what-remains-open) records the results and limits of each approach.

## Criteria for progress

The first calculation should contain two calibrated clocks, an interaction or signal linking their records, and a dynamical environment. It should define the parent symmetry and charge or balance law, solve the clock correlations, account for exchanged energy in every sector, and verify reference-clock independence.

The first gravitational target is the stationary weak-field clock ratio derived without inserting its potential dependence in advance. Different clock mechanisms must then be compared to distinguish a universal relation from altered clock physics. Recovering gravity more broadly also requires spatial and causal relations, light propagation, free-fall trajectories, and field dynamics; a scalar clock weight alone does not determine them.

Three outcomes should be distinguished:

- A reformulation reproduces known observations through a change of description.
- An emergent derivation recovers those observations from a physically specified model whose assumptions do not already contain the target relation.
- A new prediction changes a dimensionless observable and supplies a quantitative test.

The central open question is whether comparative phase accumulation and conserved exchange can be derived together in a relational system, so that gravitational clock comparisons follow without being prescribed.
