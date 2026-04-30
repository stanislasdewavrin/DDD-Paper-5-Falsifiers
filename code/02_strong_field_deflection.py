"""
Paper V — strong-field deflection: DDD eikonal vs PPN.

The eikonal calculation of Paper IV gives a finite-impact-parameter
deflection of the form
    Delta_theta_DDD = (4 G M / c^2 b) * [1 + C_DDD * (beta/b) + ...]
with C_DDD ~ 1.2 measured numerically (the 2.4*beta/b correction to
the relative dev = 2 * C_DDD * beta/b in the leading 4GM/c^2 b form).

Standard post-post-Newtonian (Schwarzschild geodesic) gives
    Delta_theta_GR_PPN = (4 G M / c^2 b) * [1 + (15 pi / 16) * (G M / c^2 b) + ...]
                      = (4 G M / c^2 b) * [1 + (15 pi / 32) * (beta / b) + ...]
i.e. C_PPN = 15*pi/32 ~ 1.47.

This script reads Paper IV deflection data, fits the leading correction
coefficient, and compares to PPN. The result is a structural prediction:
DDD predicts ~80% of the PPN second-order correction, distinguishable in
strong-field gravitational lensing (e.g. close-passage starlight near
neutron stars or near-light-cone trajectories around black holes).
"""
import json
import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
DATA = HERE / "data"; DATA.mkdir(exist_ok=True)

# Load Paper IV deflection data (combine original + extended-b)
PAPER_IV_DATA = Path("/sessions/sweet-determined-tesla/mnt/Physique"
                     "/paperIV_gravitational_phenomenology/data")

with open(PAPER_IV_DATA / "01_photon_deflection.json") as f:
    d_orig = json.load(f)
with open(PAPER_IV_DATA / "01b_deflection_high_b.json") as f:
    d_ext = json.load(f)

# Combine
seen = set()
b_all = []
defl_all = []
gr_all = []
for r in d_orig["results"]:
    bb = r["b_over_beta"]
    if bb in seen: continue
    seen.add(bb)
    b_all.append(bb)
    defl_all.append(r["deflection_alpha2"])
    gr_all.append(r["deflection_GR_pred"])
for r in d_ext["results"]:
    bb = r["b_over_beta"]
    if bb in seen: continue
    seen.add(bb)
    b_all.append(bb)
    defl_all.append(r["deflection_alpha2"])
    gr_all.append(r["deflection_GR_pred"])

idx = np.argsort(b_all)
b_all = np.array(b_all)[idx]
defl_all = np.array(defl_all)[idx]
gr_all = np.array(gr_all)[idx]

# Fit: relative deviation = (deflection - GR) / GR = C * (beta/b) + O(beta^2/b^2)
# We have beta = 1, so beta/b = 1/b_all
beta = 1.0
rel_dev = (defl_all - gr_all) / gr_all
inv_b = beta / b_all

# Linear fit on the asymptotic regime b/beta >= 50
mask = b_all >= 50
slope, intercept = np.polyfit(inv_b[mask], rel_dev[mask], 1)
C_DDD = float(slope)
C_PPN = 15 * np.pi / 32

print("Paper V — strong-field deflection coefficient")
print("=" * 60)
print(f"{'b/beta':>8} {'rel_dev':>10} {'beta/b':>10}")
for bb, rd, ib in zip(b_all, rel_dev, inv_b):
    print(f"{bb:8.1f} {rd:10.4f} {ib:10.4f}")
print()
print(f"Linear fit on b/beta >= 50:")
print(f"  rel_dev = {slope:.4f} * (beta/b) + {intercept:.4e}")
print(f"  C_DDD (fitted)         = {C_DDD:.3f}")
print(f"  C_PPN (Schwarzschild)  = 15*pi/32 = {C_PPN:.3f}")
print(f"  Ratio C_DDD / C_PPN    = {C_DDD/C_PPN:.3f}")
print()
print(f"DDD predicts ~{100*C_DDD/C_PPN:.0f}% of the PPN second-order coefficient.")
print(f"The (C_DDD - C_PPN) excess in strong-field deflection is")
print(f"a structural prediction of the present photon-substrate ansatz.")

out = {
    "C_DDD_fitted": C_DDD,
    "C_PPN_GR": C_PPN,
    "ratio_DDD_to_PPN": float(C_DDD / C_PPN),
    "fit_window_b_over_beta": "[50, infinity]",
    "data_points": [
        {"b_over_beta": float(b), "rel_dev": float(r), "beta_over_b": float(ib)}
        for b, r, ib in zip(b_all, rel_dev, inv_b)
    ],
}
with open(DATA / "02_strong_field_deflection.json", "w") as f:
    json.dump(out, f, indent=2)
print()
print(f"Saved {DATA / '02_strong_field_deflection.json'}")
