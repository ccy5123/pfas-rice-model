#!/usr/bin/env python3
# =============================================================================
# validation/hwang2017_lettuce.py
# -----------------------------------------------------------------------------
# The NEUTRAL path against Hwang, Lee & Kim 2017 (PLOS ONE 12(2):e0172254) --
# chlorpyrifos taken up by lettuce from treated soil.  Handoff task A3.
#
# WHY THIS DATASET. It is the only TIME-RESOLVED PER-ORGAN neutral dataset to
# hand (3 sampling times x 2 soil levels), it is LIPOPHILIC (log Kow 4.01, where
# Ge's difenoconazole sits, on the falling limb of the Briggs bell), and -- unlike
# every other candidate -- it supplies a MEASURED Kd, so a soil concentration
# converts into the pore-water exposure the model actually wants instead of being
# an unusable soil-basis number.
#
# WHAT IT CAN AND CANNOT SETTLE -- read this before quoting any number below.
# Four limits stack, and two of them are not resolvable from the article:
#
#   1. LETTUCE, NOT RICE, and ONE compound.  A single point cannot test a QSPR;
#      it tests the dynamics at one log Kow.
#   2. FRESH-vs-DRY BASIS OF TABLE 1 IS NOT STATED IN THE ARTICLE.  For a tissue
#      at ~95 % water this is a ~20x lever, and it SPANS THE VERDICT (section 3):
#      read one way the model over-predicts the root, read the other it cannot
#      reach it.  This is the same trap that produced the Tang fresh/dry artifact
#      and the lipid-basis error, so it is quantified rather than assumed away.
#   3. THE GROWTH CURVE'S FUNCTIONAL FORM IS NOT IN THE TRANSCRIPTION.  Table 2
#      gives `Ig` ("log initial plant weight") and `Kg` ("plant growth constant")
#      but not the equation they parameterise.  Leaf concentration scales with
#      Q_TP/M_leaf, so this propagates directly -- section 5 scans it.
#   4. THE ROOTS GREW IN SOIL.  Chlorpyrifos has Koc ~2219, so soil adhering to
#      the rhizoplane carries residue; the root endpoint is contact-confounded in
#      exactly the way Li 2025 is (docs/twopool_root_exploration.md). Section 4
#      quantifies how much of the root signal that can explain.
#
#   5. The in-plant half-life `Tp` is NOT the authors' measurement -- they say so
#      explicitly (section 6). It must NOT be presented as independent
#      confirmation of the ~7 d the Ge scan predicts.
#
# So this is a SECONDARY, CONDITION-LIMITED CHECK, not an a-priori validation on
# the footing of Liu 2023 / Ge 2017. What it does deliver is a second, independent
# sighting of an ALREADY-OPEN problem: the Briggs root partition under-predicts
# measured root concentrations of lipophilic compounds in soil-grown plants (the
# Brunetti 2021 `K_RW` = 13.3 vs Briggs ~1.0 disagreement, docs/neutral_dpu_
# validation.md section 5).
#
#   python validation/hwang2017_lettuce.py
# =============================================================================
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT_DIR, "src"))

from pfas_rice_plant_module_4pool_surf import binding_factors  # noqa: E402
import neutral_dpu as ND  # noqa: E402

# ---------------------------------------------------------------------------
# The data, transcribed from the article (docs/HANDOFF_neutral_next.md section 6)
# ---------------------------------------------------------------------------
LOG_KOW = 4.01               # Table 2: Kow = 1.02e4
KD = 82.1                    # Table 2, mL/g == L/kg, MEASURED (adsorption isotherm)
QW_L_PER_DAY = 0.0468        # Table 2: 46.8 mL/d, MEASURED (potometer)

# Table 2, by nominal soil treatment [mg/kg]
SOIL = {
    10: dict(C0=15.2, T_soil=17.2, Tp=8.7, Ig=0.3062, Kg=1.1020),
    20: dict(C0=24.9, T_soil=7.9, Tp=8.0, Ig=0.3092, Kg=1.2031),
}

# Table 1, measured residues [mg/kg], mean of n=3. BASIS NOT STATED (limit 2).
#   (level, day) -> (leaf, root, whole)
TABLE1 = {
    (10, 21): (0.5, 7.5, 0.9), (10, 30): (0.4, 5.6, 0.6), (10, 40): (0.1, 1.8, 0.2),
    (20, 21): (0.8, 2.1, 0.8), (20, 30): (0.3, 3.6, 0.5), (20, 40): (0.1, 0.4, 0.1),
}

# Tissue water contents. Lettuce leaf is ~95-96 % water (a head lettuce is
# proverbially water); the ROOT value is NOT measured anywhere in reach and is
# carried over from the repo's rice root -- an assumption, scanned in section 5.
LETTUCE_WATER = {"root": 0.90, "stem": 0.92, "leaf": 0.95, "grain": 0.95}
ROOT_MASS_FRACTION = 0.054   # DERIVED FROM TABLE 1 ITSELF -- see check_internal_consistency


# ---------------------------------------------------------------------------
# exposure and growth
# ---------------------------------------------------------------------------
def exposure(level: int, t):
    """Pore-water concentration [ug/L] from the authors' own exposure model.

        Ce(t) = C0 * (1/2)^(t/T) / Kd

    i.e. the soil residue decays first order and the measured solid-water
    distribution coefficient converts it to a soil-solution concentration. Every
    input here is measured (C0 by analysis, Kd by adsorption isotherm) except the
    soil half-life T, which the authors flag as taken from the literature
    (section 6) -- so the exposure TRAJECTORY carries that borrowed rate even
    though its level does not.
    """
    p = SOIL[level]
    return p["C0"] * 0.5 ** (np.asarray(t, float) / p["T_soil"]) / KD * 1e3


def growth(level: int, t, form: str = "log10"):
    """Total plant fresh mass [kg] over time.

    Table 2 reports `Ig` = "log initial plant weight" [g] and `Kg` = "plant growth
    constant" for a curve fitted at r > 0.99, but NOT the equation. Of the forms
    that fit those two symbols, only a LOG-LOG one gives lettuce-scale masses:

        log10  : log10 M = Ig + Kg*log10(t)   -> M(1 d) = 2.0 g, M(40 d) = 118 g
        ln     : ln M    = Ig + Kg*ln(t)      -> M(1 d) = 1.4 g, M(40 d) =  75 g
        exp    : M = 10^Ig * exp(Kg*t)        -> M(40 d) = 10^19 g  (rejected)

    Both surviving forms are plausible for a lettuce head 40 days from transplant,
    so `log10` is used as the primary and `ln` as the sensitivity (section 5).
    This is a RECONSTRUCTION, not a transcription: treat absolute leaf
    concentrations as carrying its uncertainty.
    """
    p = SOIL[level]
    t = np.maximum(np.asarray(t, float), 1.0)     # the fit is undefined at t=0
    if form == "log10":
        g = 10.0 ** (p["Ig"] + p["Kg"] * np.log10(t))
    elif form == "ln":
        g = np.exp(p["Ig"] + p["Kg"] * np.log(t))
    elif form == "exp":
        g = 10.0 ** p["Ig"] * np.exp(p["Kg"] * t)
    else:
        raise ValueError(f"unknown growth form {form!r}")
    return g / 1e3                                 # g -> kg


def drivers(level: int, season: float = 40.0, n_t: int = 161,
            growth_form: str = "log10", root_frac: float = ROOT_MASS_FRACTION):
    """`simulate_neutral` drivers for one soil treatment.

    Lettuce is a rosette: it is root + leaf, with no meaningful stem and no fruit.
    The 4-compartment core is filled accordingly -- a token stem (1 % of mass, so
    the xylem has something to pass through) and the grain pinned at the mass
    floor, which the formation gate then keeps unloaded. The endpoints that mean
    anything here are ROOT and LEAF.
    """
    t = np.linspace(0.0, float(season), int(n_t))
    total = growth(level, t, growth_form)
    root = root_frac * total
    stem = 0.01 * total
    leaf = np.maximum(total - root - stem, 1e-9)
    M = np.column_stack([root, stem, leaf, np.full_like(t, 1e-9)])
    return dict(t=t, Cwo=exposure(level, t), Qtp=np.full_like(t, QW_L_PER_DAY), M=M)


def lettuce_compartments(half_life: float | None = None, waters: dict | None = None):
    """Briggs composition on lettuce water contents, Trapp 1994 lipids."""
    gam = 0.0 if not half_life else float(np.log(2.0) / half_life)
    return ND.rice_compartments(waters=dict(LETTUCE_WATER, **(waters or {})),
                                gammas=dict.fromkeys(LETTUCE_WATER, gam))


def predict(level: int, half_life: float | None = None, growth_form: str = "log10",
            waters: dict | None = None, root_frac: float = ROOT_MASS_FRACTION):
    """Run the neutral path for one treatment; returns the simulate_neutral dict."""
    c = ND.NeutralCompound("chlorpyrifos", LOG_KOW, MW=350.6, K_AW=4.8e-4)
    return ND.simulate_neutral(
        c, drivers(level, growth_form=growth_form, root_frac=root_frac),
        comps=lettuce_compartments(half_life, waters))


def _conc_at(res, tissue, day):
    return float(np.interp(day, res["t"], res["conc"][tissue]))


# ---------------------------------------------------------------------------
# 1. does Table 1 hang together?
# ---------------------------------------------------------------------------
def check_internal_consistency():
    """Is `whole` the mass-weighted mean of leaf and root -- and if so, what root
    mass fraction does Table 1 imply?

    This matters for two reasons. It tests that the three columns share ONE basis
    (a mass-weighted mean only closes if they do), and the implied root mass
    fraction is itself a weak read on WHICH basis, because a lettuce root is a
    much larger share of the plant on dry weight (~11 %) than on fresh (~5 %).
    It also supplies the root:shoot split the model run needs, from the data
    rather than from an assumption.
    """
    print("=" * 84)
    print("1. INTERNAL CONSISTENCY — is `whole` the mass-weighted mean of leaf+root?")
    print("=" * 84)
    print("   whole = leaf*(1-x) + root*x  solved for the root mass fraction x.")
    print("   If the three columns share one basis and one mass split, x must come")
    print("   out the SAME at every sampling — which is a real test, not a fit.")
    print(f"\n{'soil':>6}{'day':>5}{'leaf':>8}{'root':>8}{'whole':>8}{'x':>10}")
    xs = []
    for (lvl, day), (lf, rt, wh) in sorted(TABLE1.items()):
        x = (wh - lf) / (rt - lf) if rt != lf else float("nan")
        note = ""
        if x <= 0.005:
            note = "  <- forced to 0 by 1-dp rounding"
        else:
            xs.append(x)
        print(f"{lvl:>6}{day:>5}{lf:>8.1f}{rt:>8.1f}{wh:>8.1f}{x:>10.3f}{note}")
    m, s = float(np.mean(xs)), float(np.std(xs))
    print(f"\n   x = {m:.3f} +/- {s:.3f} over the {len(xs)} rows the rounding leaves usable.")
    print("   => the columns ARE mutually consistent on one basis and one mass split.")
    print(f"   The implied root mass fraction {m:.1%} is characteristic of FRESH weight")
    print("   for lettuce; the same plant on dry weight would put ~11 % in the root")
    print("   (leaf ~95 % water vs root ~90 %). Weak evidence, not proof — lettuce")
    print("   root:shoot varies with pot size and harvest date — but it points the")
    print("   same way as the convention for produce residues, which are reported on")
    print("   an as-eaten FRESH basis.")
    return m, s


# ---------------------------------------------------------------------------
# 2/3. the a-priori comparison, and the basis ambiguity
# ---------------------------------------------------------------------------
def observed_baf(level, day, tissue, basis="fw"):
    """Measured tissue conc / the modelled pore-water exposure at that time [L/kg].

    `basis='dw'` converts the measurement onto the model's fresh-weight basis with
    the run's own water contents, exactly as the shared obs harness does.
    """
    lf, rt, _ = TABLE1[(level, day)]
    val = {"leaf": lf, "root": rt}[tissue] * 1e3        # mg/kg -> ug/kg
    if basis == "dw":
        val *= (1.0 - LETTUCE_WATER[tissue])
    return val / float(exposure(level, day))


def _rmse(pairs):
    if not pairs:
        return float("nan")
    e = [(np.log10(max(p, 1e-9)) - np.log10(max(o, 1e-9))) ** 2 for p, o in pairs]
    return float(np.sqrt(np.mean(e)))


def compare(half_life=None, growth_form="log10", basis="fw", quiet=False,
            by_tissue=False):
    """Model vs Hwang for both treatments and all three samplings.

    `by_tissue` also returns the per-tissue errors, which is where this dataset
    turns out to be most informative -- the two readings of the basis fail on
    OPPOSITE organs (section 3).
    """
    if not quiet:
        hl = "none (gamma=0)" if not half_life else f"{half_life:.1f} d"
        print(f"\n   basis={basis}, half-life={hl}, growth={growth_form}")
        print(f"{'soil':>6}{'day':>5}{'Cw ug/L':>10}{'tissue':>8}"
              f"{'obs':>10}{'model':>10}{'ratio':>9}")
    per = {"root": [], "leaf": []}
    for level in sorted(SOIL):
        res = predict(level, half_life, growth_form)
        for day in (21, 30, 40):
            for tissue in ("root", "leaf"):
                o = observed_baf(level, day, tissue, basis)
                p = _conc_at(res, tissue, day) / float(exposure(level, day))
                per[tissue].append((p, o))
                if not quiet:
                    print(f"{level:>6}{day:>5}{float(exposure(level, day)):>10.1f}"
                          f"{tissue:>8}{o:>10.2f}{p:>10.2f}{p / o:>9.2f}")
    pairs = per["root"] + per["leaf"]
    rmse = _rmse(pairs)
    if not quiet:
        print(f"   log10 RMSE (n={len(pairs)}) = {rmse:.3f}"
              f"   [root {_rmse(per['root']):.3f}, leaf {_rmse(per['leaf']):.3f}]")
    if by_tissue:
        return rmse, {k: _rmse(v) for k, v in per.items()}
    return rmse


def check_basis():
    """THE decisive limitation: the unstated basis spans the verdict.

    The model's root BAF has a hard ceiling -- the root cannot exceed its
    equilibrium partition K_PW, because that is what the compartment IS. So the
    question "does the model reach the measured root?" has a yes/no answer that
    flips with the basis, and no amount of parameter freedom changes it.
    """
    print("\n" + "=" * 84)
    print("2. THE UNSTATED BASIS SPANS THE VERDICT (the reason A3 cannot conclude)")
    print("=" * 84)
    comps = lettuce_compartments()
    core = ND.neutral_compound(ND.NeutralCompound("chlorpyrifos", LOG_KOW))
    K = binding_factors(comps, core)
    names = ("root", "stem", "leaf", "grain")
    kpw = {n: float(K[i]) for i, n in enumerate(names)}
    print(f"   Briggs K_PW at log Kow {LOG_KOW}: root {kpw['root']:.1f}, "
          f"leaf {kpw['leaf']:.1f} L/kg fw")
    print(f"   TSCF (Briggs bell) = {ND.briggs_tscf(LOG_KOW):.3f} — the falling limb:")
    print("   a lipophilic compound is predicted to stay in the root.")
    print("\n   The ROOT is an equilibrium partition, so K_PW is a structural CEILING")
    print("   on the modelled root BAF. Where does the measurement sit relative to it?")
    print(f"\n{'soil':>6}{'day':>5}{'root BAF (fw)':>15}{'root BAF (dw)':>15}"
          f"{'ceiling':>10}")
    for (lvl, day) in sorted(TABLE1):
        print(f"{lvl:>6}{day:>5}{observed_baf(lvl, day, 'root', 'fw'):>15.1f}"
              f"{observed_baf(lvl, day, 'root', 'dw'):>15.1f}{kpw['root']:>10.1f}")
    fw = [observed_baf(l, d, "root", "fw") for (l, d) in TABLE1]
    dw = [observed_baf(l, d, "root", "dw") for (l, d) in TABLE1]
    print(f"\n   FRESH-weight reading: root BAF {min(fw):.0f}-{max(fw):.0f}, i.e. "
          f"{min(fw) / kpw['root']:.1f}-{max(fw) / kpw['root']:.1f}x ABOVE the ceiling")
    print("     => the model structurally CANNOT reach the measurement.")
    print(f"   DRY-weight reading:   root BAF {min(dw):.1f}-{max(dw):.1f}, i.e. "
          f"{min(dw) / kpw['root']:.2f}-{max(dw) / kpw['root']:.2f}x the ceiling")
    print("     => the measurement sits UNDER it and the model over-predicts.")
    print("\n   One unstated table footnote therefore decides whether this dataset")
    print("   REFUTES the partition core or CONFIRMS it. That is the finding: Hwang")
    print("   cannot adjudicate anything until the basis is known. Section 1 leans")
    print("   fresh weight, which is the direction that matters — see section 6.")
    return kpw


# ---------------------------------------------------------------------------
# 4. how much of the root signal can adhering soil explain?
# ---------------------------------------------------------------------------
def check_soil_contact():
    """The roots grew in soil at a known concentration, so contact is a confound.

    Chlorpyrifos Koc ~2219 keeps it on the solids, so soil left on the rhizoplane
    reads as root residue. The bound is easy: soil at C0 is the MOST concentrated
    thing the root can be contaminated with, so the fraction of root mass that
    would have to be soil to explain the measurement is a hard, assumption-free
    statement about how much of the signal contact can account for.
    """
    print("\n" + "=" * 84)
    print("4. ROOT SOIL-CONTACT CONFOUND — how much can it explain?")
    print("=" * 84)
    print("   If a fraction phi of the reported root mass were adhering soil at the")
    print("   soil concentration C0, the apparent residue would be phi*C0. Solving")
    print("   for the phi that would explain the WHOLE measurement:")
    print(f"\n{'soil':>6}{'day':>5}{'C0 mg/kg':>10}{'root mg/kg':>12}{'phi needed':>12}")
    for (lvl, day), (_, rt, _) in sorted(TABLE1.items()):
        phi = rt / SOIL[lvl]["C0"]
        print(f"{lvl:>6}{day:>5}{SOIL[lvl]['C0']:>10.1f}{rt:>12.1f}{phi:>12.1%}")
    print("\n   The 10 mg/kg treatment would need 12-49 % of the washed root mass to")
    print("   be soil, which is not credible for washed roots. So contact is a REAL")
    print("   confound but CANNOT explain a 3-10x exceedance on its own — the fresh-")
    print("   weight discrepancy in section 2 survives it. (The same reasoning is")
    print("   why Li 2025's root TF is called surface-confounded but not dismissed.)")


# ---------------------------------------------------------------------------
# 5. what the unresolved inputs are worth
# ---------------------------------------------------------------------------
def sensitivity():
    """Span the three things this dataset does not pin down.

    None of these is a fitted parameter -- they are scanned to show how much of
    the comparison rests on them, which is the only honest way to report a result
    that depends on a reconstruction.
    """
    print("\n" + "=" * 84)
    print("5. SENSITIVITY to the inputs the article does not supply")
    print("=" * 84)
    print("   (a) in-plant half-life. Tp = 8.7 / 8.0 d is in Table 2 but is NOT the")
    print("       authors' measurement (section 6) — scanned, not adopted.")
    print(f"\n{'half-life':>14}{'RMSE (fw)':>12}{'RMSE (dw)':>12}")
    for hl in (None, 30.0, 8.4, 4.0):
        lbl = "none" if hl is None else f"{hl:.1f} d"
        print(f"{lbl:>14}{compare(hl, quiet=True, basis='fw'):>12.3f}"
              f"{compare(hl, quiet=True, basis='dw'):>12.3f}")
    print("\n   (b) growth-curve form (Ig/Kg parameterise an equation not in the")
    print("       transcription; `exp` is rejected on magnitude, see growth()).")
    print(f"\n{'form':>14}{'M(40 d) g':>12}{'RMSE (fw)':>12}{'RMSE (dw)':>12}")
    for form in ("log10", "ln"):
        m40 = float(growth(10, 40.0, form)) * 1e3
        print(f"{form:>14}{m40:>12.0f}"
              f"{compare(None, form, 'fw', quiet=True):>12.3f}"
              f"{compare(None, form, 'dw', quiet=True):>12.3f}")
    print("\n   (c) lettuce ROOT water content (assumed 0.90; not measured anywhere")
    print("       in reach). It enters the dw reading twice — once converting the")
    print("       measurement, once through K_PW — so it is scanned on both.")
    print(f"\n{'root water':>14}{'RMSE (fw)':>12}{'RMSE (dw)':>12}")
    old = LETTUCE_WATER["root"]
    try:
        for w in (0.85, 0.90, 0.95):
            LETTUCE_WATER["root"] = w
            print(f"{w:>14.2f}{compare(None, 'log10', 'fw', quiet=True):>12.3f}"
                  f"{compare(None, 'log10', 'dw', quiet=True):>12.3f}")
    finally:
        LETTUCE_WATER["root"] = old


# ---------------------------------------------------------------------------
# 6. verdict
# ---------------------------------------------------------------------------
def main():
    print("HWANG 2017 (lettuce / chlorpyrifos) — NEUTRAL path, handoff task A3")
    print("PLOS ONE 12(2):e0172254 · log Kow 4.01 · Kd 82.1 L/kg (measured) ·")
    print("Q_w 46.8 mL/d (measured) · exposure Ce(t) = C0*(1/2)^(t/T)/Kd\n")
    check_internal_consistency()
    kpw = check_basis()
    print("\n" + "=" * 84)
    print("3. MODEL vs MEASUREMENT, both readings of the basis")
    print("=" * 84)
    rmse, per = {}, {}
    for basis in ("fw", "dw"):
        rmse[basis], per[basis] = compare(None, "log10", basis, by_tissue=True)
    print("\n   THE TWO READINGS FAIL ON OPPOSITE ORGANS:")
    print(f"\n{'basis':>8}{'root':>10}{'leaf':>10}   what it looks like")
    print(f"{'fw':>8}{per['fw']['root']:>10.3f}{per['fw']['leaf']:>10.3f}"
          "   leaf predicted well, root unreachable")
    print(f"{'dw':>8}{per['dw']['root']:>10.3f}{per['dw']['leaf']:>10.3f}"
          "   root bracketed, leaf over by 3-24x")
    print("\n   Neither reading makes the plant coherent — one organ is always wrong,")
    print("   and swapping the basis just moves which. That is a stronger statement")
    print("   than either RMSE: it says the discrepancy is NOT a units artifact that")
    print("   the right footnote would dissolve. Something about a lipophilic")
    print("   compound in a soil-grown plant is genuinely missing, and the basis")
    print("   decides only WHERE the model is wrong, not WHETHER.")
    check_soil_contact()
    sensitivity()

    print("\n" + "=" * 84)
    print("VERDICT")
    print("=" * 84)
    print(f"  a-priori log10 RMSE: {rmse['fw']:.3f} reading Table 1 as FRESH weight,")
    print(f"                       {rmse['dw']:.3f} reading it as DRY weight.")
    print("  NEITHER NUMBER IS A VALIDATION RESULT and neither should be quoted as")
    print("  one. Four limits stack (see the module header), two of them — the")
    print("  unstated basis and the reconstructed growth curve — unresolvable from")
    print("  the article.")
    print("\n  What this run DOES establish:")
    print("    * Table 1 is internally consistent on ONE basis (section 1): `whole`")
    print("      is the mass-weighted mean of leaf and root at a root mass fraction")
    print("      of 5.4 +/- 0.9 %, the same at every sampling. So the columns can be")
    print("      read together, and the implied fraction points to FRESH weight.")
    print("    * THE TWO READINGS FAIL ON OPPOSITE ORGANS (section 3). Fresh weight")
    print("      predicts the leaf well and cannot reach the root; dry weight")
    print("      brackets the root and over-predicts the leaf 3-24x. The basis")
    print("      therefore decides WHERE the model is wrong, not WHETHER — so the")
    print("      discrepancy is not a units artifact awaiting a footnote.")
    print(f"    * On the fresh reading the measured root exceeds the model's")
    print(f"      structural partition ceiling K_PW = {kpw['root']:.1f} L/kg by 3-10x, and")
    print("      soil contact cannot account for it (section 4). That is the SAME")
    print("      DIRECTION and a comparable magnitude to the already-open Brunetti")
    print("      2021 disagreement (calibrated pea root K_RW = 13.3 vs a Briggs")
    print("      K_PW of ~1.0 for carbamazepine) — a second, independent sighting")
    print("      of one problem: for LIPOPHILIC compounds in SOIL-GROWN plants the")
    print("      Briggs root partition looks too low.")
    print("\n  A TRAP TO NAME EXPLICITLY. The half-life scan (section 5a) improves the")
    print("  DRY reading a lot (0.73 -> 0.30) and makes the FRESH one steadily worse,")
    print("  so the fit 'prefers' dry weight. That is NOT evidence about the basis:")
    print("  choosing the reading that lets the model fit is circular, and it would")
    print("  be choosing against the only non-circular evidence available (section 1")
    print("  and the convention that produce residues are reported as-eaten). The")
    print("  basis is an open question about the paper, not about the model.")
    print("\n  Cheapest way to close it: ask the authors which basis Table 1 uses.")
    print("  Until then this stays a secondary check — recorded, not quoted.")
    return rmse


if __name__ == "__main__":
    main()
