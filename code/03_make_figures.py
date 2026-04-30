"""Paper V — figure generation."""
import json, numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path("/sessions/sweet-determined-tesla/mnt/Physique/paperV_falsifiers")
DATA = HERE / "data"
FIG = HERE / "figures"; FIG.mkdir(exist_ok=True)

# ============================================================
# Figure 1: Yukawa exclusion plot (alpha_Y vs lambda_Y)
# ============================================================
print("Fig 1: Yukawa exclusion...")
ELL_P = 1.616e-35

# Eot-Wash schematic (from script 01)
lam_bound = np.array([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10, 100, 1e3, 1e4, 1e5])
aY_max = np.array([1e2, 1e0, 1e-1, 1e-2, 5e-3, 5e-4, 1e-3, 1e-2, 1e-1, 1.0, 10])

fig, ax = plt.subplots(figsize=(8, 5.5))
# Excluded region (above the curve)
lam_grid = np.logspace(-6, 6, 400)
aY_grid = 10**np.interp(np.log10(lam_grid), np.log10(lam_bound), np.log10(aY_max))
# Cap edges
aY_grid = np.where(lam_grid < 1e-5, 1e3, aY_grid)
aY_grid = np.where(lam_grid > 1e5, 1e3, aY_grid)
ax.fill_between(lam_grid, aY_grid, 1e4, color='red', alpha=0.2,
                label='excluded by Eöt-Wash')
ax.plot(lam_bound, aY_max, 'r-', lw=2, label='Eöt-Wash 95% CL bound')

# DDD natural prediction: lambda ~ ell_P, alpha_Y ~ 1
ax.plot([ELL_P], [1.0], 'b*', ms=18, mec='black', mew=0.8,
        label=r'DDD natural lattice ($\lambda_Y \sim \ell_P$, $\alpha_Y \sim 1$)')

# Annotate the safe regions
ax.annotate('safe (below sensitivity)',
            xy=(1e-7, 0.3), color='green', fontsize=9, ha='left')
ax.annotate('safe (above sensitivity)',
            xy=(3e5, 0.3), color='green', fontsize=9, ha='left')

ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(1e-8, 1e6); ax.set_ylim(1e-5, 1e3)
ax.set_xlabel(r'Yukawa length scale $\lambda_Y$  [m]')
ax.set_ylabel(r'Yukawa strength $\alpha_Y$')
ax.set_title(r'Eöt-Wash exclusion plot and the DDD natural prediction')
ax.legend(loc='upper right', fontsize=10)
ax.grid(alpha=0.3, which='both')
fig.tight_layout()
fig.savefig(FIG / "fig1_yukawa_exclusion.pdf", bbox_inches='tight')
fig.savefig(FIG / "fig1_yukawa_exclusion.png", dpi=150, bbox_inches='tight')
plt.close(fig)


# ============================================================
# Figure 2: strong-field deflection coefficient
# ============================================================
print("Fig 2: strong-field deflection...")
with open(DATA / "02_strong_field_deflection.json") as f:
    d = json.load(f)
pts = d["data_points"]
b_over_beta = np.array([p["b_over_beta"] for p in pts])
rel_dev = np.array([p["rel_dev"] for p in pts])
inv_b = np.array([p["beta_over_b"] for p in pts])
C_DDD = d["C_DDD_fitted"]
C_PPN = d["C_PPN_GR"]

fig, ax = plt.subplots(figsize=(8, 5.5))
mask = b_over_beta >= 50
ax.loglog(inv_b, rel_dev, 'ro', ms=8, label='DDD eikonal (Paper IV data)')
# DDD fit line
inv_b_fit = np.logspace(-3.2, -1.2, 50)
ax.loglog(inv_b_fit, C_DDD * inv_b_fit, 'r-', lw=1.5,
          label=rf'DDD fit: $C_{{\rm DDD}} = {C_DDD:.3f}$')
# PPN line
ax.loglog(inv_b_fit, C_PPN * inv_b_fit, 'b--', lw=1.5,
          label=rf'Schwarzschild PPN: $C_{{\rm PPN}} = 15\pi/32 = {C_PPN:.3f}$')

ax.set_xlabel(r'$\beta/b$')
ax.set_ylabel(r'$(\Delta\theta - 4GM/c^2 b) / (4GM/c^2 b)$')
ax.set_title(r'Strong-field deflection: DDD vs Schwarzschild PPN')
ax.legend(loc='lower right', fontsize=10)
ax.grid(alpha=0.3, which='both')
fig.tight_layout()
fig.savefig(FIG / "fig2_strong_field_C.pdf", bbox_inches='tight')
fig.savefig(FIG / "fig2_strong_field_C.png", dpi=150, bbox_inches='tight')
plt.close(fig)


# ============================================================
# Figure 3: Lorentz deviation vs gamma, with cutoff scale
# ============================================================
print("Fig 3: Lorentz deviation vs gamma...")
# Pull data from Paper III
with open("/sessions/sweet-determined-tesla/mnt/Physique/paperIII_inertia_kinematic/data/05_comoving_clock.json") as f:
    p3 = json.load(f)
res = p3["results"]
gammas = np.array([r["gamma_floq"] for r in res])
err = np.array([100 * r["rel_err"] for r in res])  # rel_err is signed; multiply by 100 for %

fig, ax = plt.subplots(figsize=(8, 5.5))
ax.semilogx(gammas, err, 'rs-', ms=7, lw=1.2, label='measured (Paper III data)')
ax.axhline(0, color='k', ls=':', alpha=0.5)
ax.fill_between([0.5, 1e23], -1, 1, color='green', alpha=0.1, label=r'$|err| < 1\%$ band')

# Mark "current observable gamma" for cosmic rays (UHECR ~ 10^11)
ax.axvline(1e11, color='blue', ls='--', alpha=0.5,
           label=r'UHECR $\gamma \sim 10^{11}$ (current limit)')
# Mark "Planck cutoff" for electron mass
ax.axvline(1e22, color='purple', ls='--', alpha=0.5,
           label=r'Planck cutoff $\gamma \sim 10^{22}$ (electron)')

ax.set_xlabel(r'$\gamma$')
ax.set_ylabel(r'$(\omega_{\rm comov} - m/\gamma) / (m/\gamma)$  [\%]')
ax.set_title(r'Sub-Lorentz deviation: measured ($\gamma \leq 5$) and conjectural extension')
ax.set_xlim(0.5, 1e23)
ax.set_ylim(-30, 5)
ax.legend(loc='lower left', fontsize=9)
ax.grid(alpha=0.3, which='both')
fig.tight_layout()
fig.savefig(FIG / "fig3_lorentz_deviation.pdf", bbox_inches='tight')
fig.savefig(FIG / "fig3_lorentz_deviation.png", dpi=150, bbox_inches='tight')
plt.close(fig)


# ============================================================
# Figure 4: falsification overview / regimes
# ============================================================
print("Fig 4: falsification overview...")
fig, ax = plt.subplots(figsize=(10, 5))

# Four predictions, each as a horizontal bar showing testability range
predictions = [
    {"label": "Yukawa scale", "x_low": 1e-12, "x_high": 1e-30,
     "obs_low": 1e-5, "obs_high": 1e5, "y": 4, "color": "#3366cc"},
    {"label": "deflection $C$", "x_low": 1e-3, "x_high": 1.0,
     "obs_low": 1e-7, "obs_high": 1e-4, "y": 3, "color": "#cc3366"},
    {"label": r"Lorentz dev.\ at $\gamma$", "x_low": 1.0, "x_high": 1e22,
     "obs_low": 1.0, "obs_high": 1e11, "y": 2, "color": "#cc9933"},
    {"label": "ringdown deviation", "x_low": 0.01, "x_high": 1.0,
     "obs_low": 0.001, "obs_high": 0.1, "y": 1, "color": "#33cc99"},
]
for p in predictions:
    # DDD predicted range
    ax.barh(p["y"], np.log10(p["x_high"] / p["x_low"]),
            left=np.log10(p["x_low"]),
            color=p["color"], alpha=0.4, height=0.5,
            label="DDD prediction" if p["y"] == 4 else None)
    # Observable range (where current data probe)
    ax.barh(p["y"] + 0.25, np.log10(p["obs_high"] / p["obs_low"]),
            left=np.log10(p["obs_low"]),
            color='gray', alpha=0.3, height=0.2,
            label="current sensitivity" if p["y"] == 4 else None)
    ax.text(np.log10(p["x_low"]) - 0.5, p["y"], p["label"],
            ha='right', va='center', fontsize=10, fontweight='bold')

ax.set_xlim(-32, 23)
ax.set_ylim(0.3, 4.7)
ax.set_xticks(np.arange(-30, 24, 5))
ax.set_xticklabels([rf'$10^{{{x}}}$' for x in np.arange(-30, 24, 5)])
ax.set_xlabel(r'characteristic scale (units depending on observable)')
ax.set_yticks([])
ax.set_title(r'Falsification regime map: DDD predicted scale vs current experimental envelope')
ax.legend(loc='upper right', fontsize=10)
ax.grid(axis='x', alpha=0.3)
fig.tight_layout()
fig.savefig(FIG / "fig4_falsification_map.pdf", bbox_inches='tight')
fig.savefig(FIG / "fig4_falsification_map.png", dpi=150, bbox_inches='tight')
plt.close(fig)

print("All four figures saved to:", FIG)
