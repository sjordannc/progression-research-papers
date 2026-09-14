"""Reproduce the figures and numerical checks in the manuscript.

Run from the submission directory with:

    python Online_Resource_1.py

Requirements: numpy, scipy, matplotlib.
Outputs: Fig1.pdf, Fig1.eps, Fig2.pdf, Fig2.eps, plus a numerical
summary printed to standard output.
"""

from pathlib import Path
from shutil import which
import subprocess

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid, quad
from scipy.optimize import root


OUT = Path(__file__).resolve().parent

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 9,
        "axes.labelsize": 9,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

# Model parameters used in Sec. 8.
m = 1.0
K = 1.0
mu0 = 0.4
BV = 0.2
kappa = 0.5
M1 = 1.2
M2 = 1.0
g1 = 0.12
g2 = 0.25
E = 5.0
target = np.array([0.5, 0.8, 0.55])  # Delta z, Delta theta_1, Delta theta_2


def organization(u):
    return 2.0 + (1.0 + kappa) * u**2


def coefficients(u, mu):
    cval = organization(u)
    az = m + mu * (1.0 + 2.0 * u**2)
    i1 = M1 + g1 * cval
    i2 = M2 + g2 * cval
    return np.asarray([az, i1, i2])


def potential(u):
    return BV * u**2 * (1.0 - u) ** 2


def endpoint_derivatives(u, charges, mu):
    avec = coefficients(u, mu)
    qc = np.sum(charges**2 / (2.0 * (E - potential(u)) * avec))
    if qc >= 1.0:
        raise ValueError("The trial charges leave the regular monotonic branch.")
    return charges / avec * np.sqrt(K / (2.0 * (E - potential(u)) * (1.0 - qc)))


def endpoint_residual(charges, mu):
    values = np.array(
        [
            quad(
                lambda x, j=j: endpoint_derivatives(x, charges, mu)[j],
                0.0,
                1.0,
                epsabs=1.0e-11,
                epsrel=1.0e-11,
            )[0]
            for j in range(3)
        ]
    )
    return values - target


def solve_charges(mu, guess):
    solution = root(lambda c: endpoint_residual(c, mu), guess, tol=1.0e-11)
    if not solution.success or np.linalg.norm(endpoint_residual(solution.x, mu)) > 1.0e-9:
        raise RuntimeError(f"Boundary-value solve failed: {solution.message}")
    return solution.x


def stationary_action(mu, charges):
    def integrand(u):
        avec = coefficients(u, mu)
        qc = np.sum(charges**2 / (2.0 * (E - potential(u)) * avec))
        return np.sqrt(2.0 * K * (E - potential(u)) / (1.0 - qc))

    return quad(integrand, 0.0, 1.0, epsabs=1.0e-11, epsrel=1.0e-11)[0]


def sensitivity_integral(charges):
    def integrand(u):
        avec = coefficients(u, mu0)
        zprime, theta1prime, theta2prime = endpoint_derivatives(u, charges, mu0)
        metric_speed_sq = (
            K
            + avec[0] * zprime**2
            + avec[1] * theta1prime**2
            + avec[2] * theta2prime**2
        )
        omega = E - potential(u)
        active_link_speed_sq = (1.0 + 2.0 * u**2) * zprime**2
        return omega * active_link_speed_sq / np.sqrt(2.0 * omega * metric_speed_sq)

    return quad(integrand, 0.0, 1.0, epsabs=1.0e-11, epsrel=1.0e-11)[0]


charges = solve_charges(mu0, np.array([1.5, 2.2, 1.65]))
u = np.linspace(0.0, 1.0, 700)
az, i1, i2 = coefficients(u, mu0)
relative_ratio = (i2 / i1) / (i2[0] / i1[0])
veff = potential(u) + charges[0] ** 2 / (2.0 * az) + charges[1] ** 2 / (
    2.0 * i1
) + charges[2] ** 2 / (2.0 * i2)

# Figure 1: connected link activation and two-clock response.
fig = plt.figure(figsize=(7.2, 2.35), constrained_layout=True)
gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.0, 1.0])

ax = fig.add_subplot(gs[0, 0])
layouts = [
    (0.0, np.array([[0.00, 0.55], [0.30, 0.05], [0.60, 0.55]]), (0, (2, 2))),
    (0.5, np.array([[0.88, 0.55], [1.18, 0.05], [1.48, 0.55]]), (0, (5, 2))),
    (1.0, np.array([[1.76, 0.55], [2.06, 0.05], [2.36, 0.55]]), "solid"),
]
for strength, pts, closure_style in layouts:
    ax.plot(pts[:2, 0], pts[:2, 1], color="black", linewidth=1.1)
    ax.plot(pts[1:, 0], pts[1:, 1], color="black", linewidth=1.1)
    ax.plot(
        [pts[2, 0], pts[0, 0]],
        [pts[2, 1], pts[0, 1]],
        color="black",
        linewidth=0.8 + 0.7 * strength,
        linestyle=closure_style,
    )
    ax.scatter(pts[:, 0], pts[:, 1], s=17, color="black", zorder=3)
    ax.text(pts[:, 0].mean(), -0.13, rf"$u={strength:g}$", ha="center", fontsize=8)
ax.text(0.01, 0.98, "(a)", transform=ax.transAxes, va="top", fontweight="bold")
ax.set_xlim(-0.13, 2.49)
ax.set_ylim(-0.24, 0.77)
ax.axis("off")

ax = fig.add_subplot(gs[0, 1])
ax.plot(u, i1, color="#0072B2", linewidth=1.7, label=r"$I_1(u)$")
ax.plot(u, i2, color="#D55E00", linewidth=1.7, linestyle="--", label=r"$I_2(u)$")
ax.text(0.02, 0.98, "(b)", transform=ax.transAxes, va="top", fontweight="bold")
ax.set_xlabel(r"link coordinate $u$")
ax.set_ylabel("clock-fiber inertia")
ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.08, 1.0))
ax.spines[["top", "right"]].set_visible(False)

ax = fig.add_subplot(gs[0, 2])
ax.plot(u, relative_ratio, color="#0072B2", linewidth=1.8)
ax.text(0.02, 0.98, "(c)", transform=ax.transAxes, va="top", fontweight="bold")
ax.set_xlabel(r"link coordinate $u$")
ax.set_ylabel(r"$R_{12}(u)/R_{12}(0)$")
ax.set_ylim(0.995, 1.12)
ax.spines[["top", "right"]].set_visible(False)

for suffix in ("pdf", "eps"):
    fig.savefig(OUT / f"Fig1.{suffix}", bbox_inches="tight")
plt.close(fig)

# Figure 2: regular branch and normalized affine parameter.
dsigma_du = np.sqrt(K / 2.0) / np.sqrt(E - veff)
sigma = cumulative_trapezoid(dsigma_du, u, initial=0.0)
sigma /= sigma[-1]

fig, (ax1, ax2) = plt.subplots(
    1, 2, figsize=(7.2, 2.35), constrained_layout=True
)
ax1.plot(u, veff, color="#0072B2", linewidth=1.7, label=r"$V_{\rm eff}(u)$")
ax1.axhline(E, color="black", linewidth=1.1, linestyle="--", label=r"$E_u$")
ax1.set_xlabel(r"link coordinate $u$")
ax1.set_ylabel("effective potential / energy")
ax1.set_ylim(2.65, 5.15)
ax1.legend(frameon=False, loc="lower left")
ax1.spines[["top", "right"]].set_visible(False)

ax2.plot(sigma, u, color="#0072B2", linewidth=1.7)
ax2.set_xlabel(r"normalized affine label $\sigma/\sigma_+$")
ax2.set_ylabel(r"link coordinate $u$")
ax2.set_xlim(0.0, 1.0)
ax2.set_ylim(0.0, 1.0)
ax2.spines[["top", "right"]].set_visible(False)

for suffix in ("pdf", "eps"):
    fig.savefig(OUT / f"Fig2.{suffix}", bbox_inches="tight")
plt.close(fig)

# Numerical checks quoted in Sec. 8.
step = 1.0e-4
charges_plus = solve_charges(mu0 + step, charges)
charges_minus = solve_charges(mu0 - step, charges)
finite_difference = (
    stationary_action(mu0 + step, charges_plus)
    - stationary_action(mu0 - step, charges_minus)
) / (2.0 * step)
direct_sensitivity = sensitivity_integral(charges)
xi = relative_ratio[-1]

print("Conserved charges (p_z, pi_1, pi_2):")
print("  " + ", ".join(f"{value:.10f}" for value in charges))
print(f"Endpoint transport ratio Xi_+-: {xi:.10f}")
print(f"Fixed-endpoint finite-difference derivative: {finite_difference:.10f}")
print(f"Direct sensitivity integral: {direct_sensitivity:.10f}")
print(f"Absolute sensitivity discrepancy: {abs(finite_difference-direct_sensitivity):.3e}")

# Optional figure integrity information when Poppler is available.
pdfinfo = which("pdfinfo")
if pdfinfo is not None:
    for name in ("Fig1.pdf", "Fig2.pdf"):
        subprocess.run([pdfinfo, str(OUT / name)], check=True, stdout=subprocess.DEVNULL)
