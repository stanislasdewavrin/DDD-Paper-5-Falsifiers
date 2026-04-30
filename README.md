# Paper V — Predictions and Falsifiers of DDD Gravity

**Title:** Discrete Drainage Dynamics V: Predictions and Falsifiers
of the Gravitational Sector

## What this paper does

Lists where DDD makes distinct testable predictions vs Newton/GR and
quantifies the experimental envelopes that would falsify them. Four
predictions:

1. **Yukawa modification** from self-consistent feedback in the scalar
   sector; Eöt-Wash bounds constrain the feedback strength `β_fb`.
2. **Strong-field deflection coefficient** `C_DDD ≈ 2.34` vs Schwarzschild
   PPN `C_PPN ≈ 1.47` — a 60% excess at second order, distinguishable in
   close-passage gravitational lensing.
3. **Ultra-relativistic Lorentz deviations** from Paper III's lattice
   corrections; sub-Lorentz residual (clock runs slower than 1/γ).
4. **Compact-object phenomenology** near the bandwidth horizon
   `r = β`; qualitative for now.

A falsification table summarises the four against current experimental
bounds: none currently excluded, each testable in principle.

## Code

- `code/01_yukawa_bounds.py` — Yukawa scale from feedback strength;
  Eöt-Wash exclusion regions
- `code/02_strong_field_deflection.py` — fit `C_DDD` from Paper IV
  ray-tracing data; compare to Schwarzschild PPN

## Build

```
make
```
