"""
Paper V — Yukawa modification of gravity from self-consistent feedback,
and its Eöt-Wash exclusion bound.

Within the mobility framework of Papers I-II, a non-linear self-consistent
treatment of the deficit field (deficit acts as its own secondary source,
not just the externally placed mass) gives a Yukawa-type correction to the
Newtonian potential at second order:

    V(r) = -(G m1 m2 / r) * [1 + alpha_Y * exp(-r / lambda_Y)]

with Yukawa length scale (in lattice units)
    lambda_Y_lat = sqrt(alpha_D / |beta_fb|)
where alpha_D is the drainage coupling of Paper I and beta_fb is the
feedback strength (the coefficient of the quadratic self-source term in
the rule). Physical scale:
    lambda_Y_phys = lambda_Y_lat * d_min
with d_min the lattice spacing (typically identified with the Planck
length under Planck-scale anchoring of Paper II).

Eöt-Wash precision tests of the inverse-square law (Adelberger 2009)
constrain alpha_Y as a function of lambda over the range
[~10 micrometres, ~10^5 metres]. For alpha_Y of order unity (the natural
strength expected from full gravitational coupling), this excludes
Yukawa scales in that range.

The output of this script is the allowed region of (alpha_D, beta_fb)
parameter space under Planck-scale anchoring.
"""
import json
import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
DATA = HERE / "data"; DATA.mkdir(exist_ok=True)

# Physical constants
ELL_PLANCK = 1.616e-35  # m

# Eöt-Wash exclusion (schematic, Adelberger 2009 95% CL bounds)
# lambda in metres, alpha_Y_max is the max Yukawa strength NOT excluded
LAMBDA_BOUNDS_M = np.array([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10, 100, 1e3, 1e4, 1e5])
ALPHA_Y_MAX     = np.array([1e2, 1e0, 1e-1, 1e-2, 5e-3, 5e-4, 1e-3, 1e-2, 1e-1, 1.0, 10])


def lambda_phys_from_lattice(alpha_D, beta_fb, ell=ELL_PLANCK):
    """Yukawa length in metres from lattice parameters."""
    if beta_fb <= 0:
        return float('inf')
    return np.sqrt(alpha_D / beta_fb) * ell


def excluded(lambda_phys, alpha_Y=1.0):
    """Return True if alpha_Y > Eöt-Wash bound at this lambda."""
    if lambda_phys < LAMBDA_BOUNDS_M.min() or lambda_phys > LAMBDA_BOUNDS_M.max():
        return False  # outside the experimental sensitivity range
    aY_max = float(np.interp(np.log10(lambda_phys),
                             np.log10(LAMBDA_BOUNDS_M),
                             np.log10(ALPHA_Y_MAX)))
    return alpha_Y > 10 ** aY_max


def main():
    print("Paper V — Yukawa scale from self-consistent feedback")
    print("=" * 60)
    print(f"Planck length: {ELL_PLANCK:.3e} m")
    print(f"Eot-Wash forbidden range: lambda_phys in [{LAMBDA_BOUNDS_M.min():.0e}, "
          f"{LAMBDA_BOUNDS_M.max():.0e}] m for alpha_Y ~ 1")
    print()

    # === Test case 1: typical lattice parameters ===
    cases = [
        # (alpha_D, beta_fb, label)
        (1.0,   1.0,   "natural lattice (alpha_D = beta_fb = 1)"),
        (0.15,  0.01,  "weak feedback"),
        (1.0,   1e-30, "very weak feedback"),
        (1.0,   1e-60, "ultra-weak feedback"),
        (1.0,   1e60,  "ultra-strong feedback"),
        (4*np.pi, 1.0, "Planck-anchored alpha_D / R_0 = 4 pi"),
    ]

    print(f"{'case':>50} {'lambda_phys (m)':>20} {'excluded?':>12}")
    rows = []
    for alpha_D, beta_fb, label in cases:
        lam = lambda_phys_from_lattice(alpha_D, beta_fb)
        exc = excluded(lam, alpha_Y=1.0)
        marker = "EXCLUDED" if exc else "allowed"
        print(f"{label:>50} {lam:20.3e} {marker:>12}")
        rows.append({
            "label": label,
            "alpha_D": alpha_D,
            "beta_fb": beta_fb,
            "lambda_phys_m": float(lam),
            "excluded_at_alpha_Y_1": exc,
        })
    print()

    # === Required range of beta_fb ===
    # For Yukawa scale safely below 10 mu m: beta_fb > alpha_D / (1e-5/ell_P)^2
    # = alpha_D / (6.2e29)^2 = alpha_D / 3.8e59
    print("Required allowed regions (assuming alpha_Y ~ 1):")
    print("  Either lambda_phys < 10^-5 m  (below Eot-Wash sensitivity)")
    print("  OR     lambda_phys > 10^5 m   (above Eot-Wash sensitivity)")
    print()
    for alpha_D in [1.0, 4*np.pi, 0.1]:
        beta_lo_below = alpha_D / (1e-5 / ELL_PLANCK) ** 2
        beta_hi_above = alpha_D / (1e5 / ELL_PLANCK) ** 2
        print(f"  alpha_D = {alpha_D}: beta_fb > {beta_lo_below:.2e} (below 10 mu m)")
        print(f"                       OR beta_fb < {beta_hi_above:.2e} (above 10^5 m)")
        rows.append({
            "alpha_D": alpha_D,
            "beta_fb_min_for_short_yukawa": float(beta_lo_below),
            "beta_fb_max_for_long_yukawa": float(beta_hi_above),
        })

    out = {
        "ell_planck_m": ELL_PLANCK,
        "lambda_bounds_m": LAMBDA_BOUNDS_M.tolist(),
        "alpha_Y_max": ALPHA_Y_MAX.tolist(),
        "results": rows,
    }
    with open(DATA / "01_yukawa_bounds.json", "w") as f:
        json.dump(out, f, indent=2, default=lambda o: float(o) if isinstance(o, (np.floating,)) else o)
    print()
    print(f"Saved {DATA / '01_yukawa_bounds.json'}")


if __name__ == "__main__":
    main()
