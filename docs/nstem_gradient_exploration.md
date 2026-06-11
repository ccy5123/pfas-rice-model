# Multi-height stem & the vertical PFAS gradient — exploration record

> Kept so the reasoning isn't re-discovered later. This documents an avenue we
> investigated for GAP B (root→shoot loading `f_xy`) and **deliberately set aside**
> in favour of the aggregate-BAF absolute fit. Code lives in
> `src/pfas_rice_plant_module_nstem.py`; the check is `validation/nstem_gradient_check.py`.

## Motivation
Yamazaki 2023 (SI Tables S18/S19) resolves rice stem PFAS by height (0–20 / 20–40 /
40–60 / >60 cm). The **vertical gradient is congener-dependent**: short chains
concentrate strongly **upward** (PFBA top/bottom ≈ 7.4, transpiration-driven),
long chains are flat/**down** (PFUnDA ≈ 0.66). A single well-mixed "straw"
compartment cannot host this, and it was offered as the reason the saturated W2
transport fit had to inflate long-chain `f_xy`. We tested whether resolving the
stem into N serial segments lets a **monotone** `f_xy` reproduce the gradient.

## What we built
`pfas_rice_plant_module_nstem.py`:
- **`NStemModel`** — N serial stem(+leaf) segments; xylem flows root→s1→…→sN→grain
  with transpiration draw-off `τ_s` (water leaves, the non-volatile anion stays →
  upward concentration) + radial **equilibrium** exchange + growth dilution.
  Mass-conserving (sole source `M_root·j_R`).
- **`NStemKineticModel`** — decouples xylem and tissue with a finite radial
  conductance `k_rad` (quasi-steady xylem; `dC_s/dt = k_rad·(Cw_xyl−Cw_s) − μ_s·C_s`);
  recovers `NStemModel` as `k_rad→∞`. Mass-conserving.

Driven by the **measured** forcings (`forcing_rice.Q_TP` = Kumari/Nay Htoon;
`growth_rice.M_s` = ORYZA IR72), i.e. no hand-tuning.

## What we found
1. **Short/mid-chain upward gradient: reproduced** (correct direction) — transpiration
   water-removal concentrates the mobile fraction upward.
2. **Long-chain reversal: NOT reproduced.** Equilibrium segments give a
   ~chain-independent gradient at realistic biomass, because the crossover
   `B* ~ Q_s/(M_s·μ_s) ≈ 100+` sits far above the congener range (`B` 1–60). An
   earlier apparent flip was an artifact of **inflated placeholder biomass**.
3. **Kinetic `k_rad` did not rescue it.** Sweeping `k_rad` (0.05–5 L/(day·kg)) the
   predicted top/bottom stays ~2.4–4.9 for all congeners vs the observed 7.4→0.66.

## Why (the diagnosis)
In both variants the tissue free conc stays in **reversible** balance with the
upward-concentrating xylem, so the **bound** fraction concentrates upward too. The
observed long-chain **down** gradient needs the high-`B` bound fraction to be
**delivery-limited and slow to re-release** — i.e. **irreversible / hysteretic
sorption** acting as a sink that strips the ascending stream low in the stem.

## The idea kept for later (if the within-stem profile is ever needed)
Add **asymmetric (hysteretic) radial kinetics**: a fast sorption rate `k_sorb`
and a much slower desorption `k_desorb ≪ k_sorb` (or an irreversibly-bound pool).
High-`B`/long-chain congeners then accumulate low and cannot re-supply the upper
segments → flat/down; short chains stay reversible → up. This is one extra
parameter (the sorption reversibility) and is the natural next step.

## Decision
The within-stem vertical profile is a **fine detail** relative to the model's
purpose (compartment BAFs), and chasing it needs increasingly elaborate sorption
physics. We therefore **pin the absolute `f_xy` via the aggregate root/straw/grain
BAF** (with the measured `Q_TP`/`M_s`) instead — that meets the actual task-2 goal
— and keep the multi-height stem + this note as a documented, mass-conserving
option to revisit with hysteretic sorption.
