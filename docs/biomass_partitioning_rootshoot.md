# Rice biomass partitioning & the root:shoot anchor

> Status: literature-grounded for the **shoot split**; the **root fraction** is a
> documented data gap (anchored to a literature *range*, not a site-specific value).
> Scope: explains the organ-biomass `M_k(t)` that feeds the burden = `C_k·M_k`
> question and the calibration coupling it creates. No default behaviour changed.

## 1. Why this matters

"Where does PFAS accumulate?" is a **burden** (`C_k·M_k`, µg) question, and burden
depends directly on the per-organ biomass `M_k`. A wrong partitioning silently
biases the burden ranking — and, as shown in §4, the model's transport calibration
turns out to be **entangled** with the biomass assumption.

## 2. What the literature says (maturity, lowland rice)

**Aboveground split — well established, and the model matches it:**

| organ | % of aboveground DM at maturity | source |
|---|---|---|
| panicle / grain | ~45–58% (HI 0.45–0.55, modern) | Ntanos & Koutroubas 2002; Amanullah & Inamullah 2016 |
| stem / culm + sheath | ~24–33% | Amanullah & Inamullah 2016 |
| leaf blade | ~18–24% | Amanullah & Inamullah 2016 |

**Root fraction — the weak link:**

- Root mass fraction declines from ~0.2 (seedling) to ~0.1 (heading) of total plant;
  at **maturity it is ~0.07–0.13** (root:shoot ~0.08–0.15), pushed lower in high-yield
  cultivars by grain-fill dilution (root:shoot negatively correlates with yield).
- A concrete anchor: heading shoot 3 Mg/ha ↔ root ~0.33 Mg/ha → root:shoot ≈ 0.11.
- **Combine, don't co-locate**: the full per-organ distribution is built from the
  shoot split (above) × an *independent* root:shoot ratio — the two need not come from
  the same paper. Combining gives (total-plant %): root 9–13, grain 46–48, stem 25–26,
  leaf 16–17 for root:shoot 0.10–0.15.

### Citations (DOI)
- Ntanos & Koutroubas 2002, *Field Crops Res.* 74:93–101 — **10.1016/S0378-4290(01)00203-9** (maturity HI 0.47–0.61; aboveground only).
- Amanullah & Inamullah 2016, *Rice Science* 23(2):78–87 — maturity panicle/culm/leaf % (aboveground only).
- Nada & Abogadallah 2018, *Acta Physiol. Plant.* 40:123 — **10.1007/s11738-018-2697-5** (root:shoot is genotype/condition dependent; pot/greenhouse).
- "Response of Grain Yield and Root … to Nitrogen Levels in Paddy Rice", *Front. Plant Sci.* 2021 — **10.3389/fpls.2021.713814** (japonica paddy, tillering→maturity SDW/RDW/root:shoot).
- Yoshida 1981, *Fundamentals of Rice Crop Science*, IRRI (classic reference for root fraction ~5–10%).
- Model partitioning source: Bouman & van Laar 2006, *Agric. Syst.* 87:249–273 (**10.1016/j.agsy.2004.09.011**); Li et al. 2017, *Agric. For. Meteorol.* 237–238:246–256 (ORYZA v3).

> Honest caveat: none of these is the model's exact target system (e.g. Yamazaki
> Andosol greenhouse). Partitioning is genotype × environment dependent, so the
> literature gives a defensible **range**, not a site truth. A maturity, root-
> inclusive measurement for the target system would close the gap.

## 3. What the model assumes vs literature

| biomass driver | root % of total | root:shoot | HI |
|---|---|---|---|
| `reproduce_demo` placeholder logistic | 23 | **0.30** | **0.07** |
| `growth_rice` (app/`simulate` default) | 3.4 | **0.035** | 0.51 |
| literature (maturity) | 7–13 | 0.08–0.15 | 0.45–0.55 |

The **shoot split of `growth_rice` matches the literature**; only its **root
fraction is too low** (DVS partitioning drives root growth to ~0 after flowering,
so the final root:shoot collapses to 0.035). The `reproduce_demo` placeholder is the
opposite extreme — root *too high* (0.30) and HI *non-physical* (0.07).

`growth_rice.organ_biomass(..., root_shoot=<r>)` now rescales the root trajectory to
a chosen final root:shoot (default `None` = unchanged), so a literature value can be
plugged in for sensitivity work.

## 4. Calibration coupling — the key finding

`validation/root_shoot_biomass_sensitivity.py` holds the **W2 transport fit**
(`params/parameters.json`) fixed and swaps the biomass driver, scoring the Yamazaki
root/straw/grain BAF (log10 RMSE):

| biomass | root:shoot | HI | W2 RMSE |
|---|---|---|---|
| reproduce placeholder | 0.30 | 0.07 | **0.029** |
| growth_rice default | 0.035 | 0.51 | 0.305 |
| growth_rice R/S=0.10 (lit) | 0.10 | 0.48 | 0.255 |
| growth_rice R/S=0.15 (lit) | 0.15 | 0.46 | 0.251 |
| growth_rice R/S=0.30 | 0.30 | 0.40 | 0.258 |

**The celebrated RMSE 0.029 reproduction is attained only with the non-physical
placeholder biomass** (root:shoot 0.30 *and* HI 0.07). With any realistic HI (~0.5)
the same W2 parameters give RMSE ~0.25–0.31 **regardless of root:shoot** — i.e. it is
the unrealistically low HI (almost no grain), not just the root, that the fit leans
on. The transport fit (`f_xy_W2fit`, `L_Ph_W2fit`, `kappa_d_W2fit`) is therefore
**entangled with a non-physical biomass** and is not valid on a literature-consistent
one.

## 5. Implications for the burden / "leaf dominates" question

- The leaf-dominant burden seen in the default `simulate` (4-pool) is partly the
  documented **leaf-sink runaway** (`nstem_leaf` redistributes it to grain≈leaf≈stem)
  and partly the **root mass being too low** (so root burden is under-counted,
  especially for root-dominated PFOS / long chains). See
  `validation/root_shoot_biomass_sensitivity.py` and the burden sensitivity in chat.
- Correcting root:shoot to ~0.10 roughly triples root burden share (PFOS root
  54%→78%) and ~doubles total uptake (bigger root = more uptake surface).

## 6. Follow-up (open)

1. **Re-fit the transport parameters** (`f_xy`, `L_Ph`, `kappa_d`) on a literature-
   consistent biomass (root:shoot ~0.10, HI ~0.5) — the current W2 fit does not
   transfer (§4). This is the substantive next step (model task #7).
2. **Close the root data gap**: a maturity, root-inclusive biomass for the target
   system to pin root:shoot (currently a literature range 0.08–0.15).
3. Decide whether to make a literature root:shoot the `growth_rice` default (kept
   `None`/unchanged for now to preserve reproducibility until the re-fit in #1).
