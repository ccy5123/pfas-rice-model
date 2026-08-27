"""Neutral / weak-electrolyte extension -- structural probe (phase 1).

Checks the four structural claims of the DPU-base extension documented in
`docs/NEUTRAL_DPU_EXTENSION_DESIGN_KR.md`:

  1. SPECIATION LIMITS -- (f_n, f_d) is the switch, not the valence.  With f_d=0 the
     GHK branch dies and the result is independent of z; with f_n=0 the neutral
     branch dies and the PFAS behaviour is recovered exactly.
  2. TSCF BELL -- root/xylem loading follows Briggs, but the resulting straw
     CONCENTRATION peaks at a HIGHER log K_ow than the TSCF bell does, because
     conc ~ TSCF * B and B rises monotonically with K_ow.
  3. ROOT EQUILIBRIUM -- with the default fast-exchange P_n the root approaches the
     Briggs/DPU equilibrium C_root = B_root * C_w^o.
  4. BRIGGS RCF COMPARISON -- the mechanistic B_k built from MEASURED rice lipid
     reproduces Briggs' K_ow slope (0.77/log unit) exactly but sits ~12x below its
     lipid coefficient, because Briggs' 0.03 is a fresh-weight "octanol equivalent",
     not an analytically measured lipid fraction.

This is STRUCTURAL verification against theory, NOT validation against data: there
is no per-organ neutral-compound time series for rice in docs/literature_db/.

Run:  python validation/neutral_probe.py
"""
from __future__ import annotations
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

import literature_params as lp                                        # noqa: E402
import model_api as api                                              # noqa: E402
from pfas_rice_plant_module_4pool_surf import (                      # noqa: E402
    Compound, Environment, root_uptake)

LOGKOW_GRID = np.arange(0.0, 6.01, 0.25)


def speciation_limits():
    """Claim 1: (f_n, f_d) selects the pathway; the valence is inert when f_d=0."""
    print("=" * 78)
    print("1. SPECIATION LIMITS  --  the switch is (f_n, f_d), not z")
    print("=" * 78)
    env = Environment()
    kw = dict(K_prot=0.0, K_PL=0.0, K_cw=0.0, Km_in=1.0, Km_out=1.0, L_Ph=1.0, f_xy=1.0)

    # neutral: fd=0 -> GHK term is dead -> the valence cannot matter
    vals = []
    for z in (-2, -1, 0, 1, 2):
        c = Compound("neutral", kappa_d=0.5, Vmax_in=0.0, Vmax_out=0.0,
                     fn=1.0, fd=0.0, P_n=10.0, z=z, **kw)
        vals.append(root_uptake(2.0, 0.5, c, env))
    same = all(v == vals[0] for v in vals)
    print(f"  neutral (fn=1, fd=0): j_R over z=-2..+2 -> {vals[0]:.6f}  identical={same}")
    assert same, "a neutral compound must not feel the membrane potential"
    assert abs(vals[0] - 10.0 * (2.0 - 0.5)) < 1e-12, "neutral term must be Fickian"

    # PFAS: fn=0 -> neutral term is dead -> P_n cannot matter
    vals = []
    for P_n in (0.0, 1.0, 1e4):
        c = Compound("pfas", kappa_d=0.5, Vmax_in=20.0, Vmax_out=8.0,
                     fn=0.0, fd=1.0, P_n=P_n, **kw)
        vals.append(root_uptake(2.0, 0.5, c, env))
    same = all(v == vals[0] for v in vals)
    print(f"  PFAS    (fn=0, fd=1): j_R over P_n=0..1e4 -> {vals[0]:.6f}  identical={same}")
    assert same, "a permanent anion must not feel the neutral permeability"

    # weak acid: both branches alive, and the split follows Henderson-Hasselbalch
    for pKa in (2.0, 4.0, 6.0, 8.0):
        fn, fd = lp.speciation(pKa, lp.PADDY_PH)
        print(f"  weak acid pKa={pKa:4.1f} at pH {lp.PADDY_PH}: f_n={fn:.4f}  f_d={fd:.4f}")
    fn, fd = lp.speciation(-3.0, lp.PADDY_PH)                 # PFSA-like
    print(f"  PFSA-like pKa=-3.0            : f_n={fn:.2e}  f_d={fd:.6f}  (-> the PFAS limit)")
    print()


def ion_trap():
    """Claim 1b: the phloem ion trap switches itself off for a permanent anion."""
    print("=" * 78)
    print("2. PHLOEM ION TRAP  --  leaf cytosol pH 7.2 -> phloem pH 8.0")
    print("=" * 78)
    print("  Lambda is the EQUILIBRIUM enrichment; P_n f_n / P_d f_d is whether the")
    print("  neutral pathway can actually deliver it.  The trap turns off for a")
    print("  permanent anion KINETICALLY (f_n -> 0), not thermodynamically.\n")
    print(f"{'pKa':>6} {'Lambda':>8} {'P_n f_n / P_d f_d':>19}")
    for pKa in (-3.0, 2.0, 4.0, 6.0, 8.0, 10.0):
        lam = lp.ion_trap_factor(pKa, 7.2, 8.0)
        dom = lp.neutral_pathway_ratio(pKa, 7.2)
        tag = "   <- PFSA-like" if pKa < 0 else ("   <- acidic herbicide" if pKa == 4.0 else "")
        print(f"{pKa:6.1f} {lam:8.3f} {dom:19.3e}{tag}")

    lam_strong = lp.ion_trap_factor(-3.0, 7.2, 8.0)
    lam_limit = 10.0 ** (8.0 - 7.2)
    print(f"\n  Strong-acid limit: Lambda -> 10**(8.0-7.2) = {lam_limit:.3f} "
          f"(got {lam_strong:.3f}), NOT 1.")
    assert abs(lam_strong - lam_limit) < 1e-3, "strong-acid Lambda must tend to 10**dpH"
    assert lp.ion_trap_factor(4.0, 7.2, 8.0) > 6.0, "a pKa-4 acid should enrich ~6x"

    dom_acid, dom_pfas = lp.neutral_pathway_ratio(4.0, 7.2), lp.neutral_pathway_ratio(-3.0, 7.2)
    print(f"  Kinetic switch-off: the neutral pathway drops {dom_acid / dom_pfas:.1e}x "
          f"from a pKa-4 acid to a PFSA.")
    assert dom_acid > 1.0, "a pKa-4 acid should be carried mainly by the neutral species"
    assert dom_pfas < 1e-5, "PFAS must have a negligible neutral pathway"
    print()


def kow_sweep():
    """Claims 2 & 3: TSCF bell vs the tissue-concentration peak, and root equilibrium."""
    print("=" * 78)
    print("3. log K_ow SWEEP  --  TSCF bell vs tissue concentration")
    print("=" * 78)
    print(f"{'logKow':>7} {'K_lip':>10} {'TSCF':>6} {'B_root':>9} | "
          f"{'root':>9} {'straw':>8} {'grain':>9} | root/B_root")
    rows = []
    for lk in LOGKOW_GRID:
        r = api.simulate_neutral(float(lk), n_t=121)
        root = float(r["baf_final"]["root"])
        b_root = float(r["B_k"]["root"])
        rows.append((float(lk), root, float(r["straw_baf"]), b_root))
        if abs(lk * 4 - round(lk * 4)) < 1e-9 and (round(lk * 4) % 4 == 0 or abs(lk - 3.75) < 1e-9
                                                   or abs(lk - 1.75) < 1e-9):
            print(f"{lk:7.2f} {r['briggs']['K_lip']:10.1f} {r['briggs']['tscf']:6.3f} "
                  f"{b_root:9.3f} | {root:9.3f} {r['straw_baf']:8.3f} "
                  f"{r['baf_final']['grain']:9.2f} | {root / b_root:.4f}")

    lk_arr = np.array([x[0] for x in rows])
    straw = np.array([x[2] for x in rows])
    eq = np.array([x[1] / x[3] for x in rows])
    peak = float(lk_arr[int(np.argmax(straw))])
    # analytic optimum of  TSCF(x) * K_ow^b  (lipid-dominated binding)
    analytic = lp.TSCF_LOGKOW_OPT + lp.BRIGGS_B * np.log(10.0) * lp.TSCF_WIDTH / 2.0
    print(f"\n  TSCF bell peaks at log K_ow = {lp.TSCF_LOGKOW_OPT}")
    print(f"  straw CONC peaks at log K_ow = {peak}   (analytic ~{analytic:.2f})")
    print("  -> the tissue peak is shifted RIGHT because conc ~ TSCF * B and B rises with K_ow.")
    assert 3.0 <= peak <= 4.5, f"straw peak {peak} outside the expected [3.0, 4.5]"
    print(f"  root equilibration C_root/B_root: min={eq.min():.4f} over the whole sweep"
          f"  (fast-exchange limit, P_n={lp.PN_DEFAULT:g})")
    assert eq.min() > 0.95, "the root should sit near the Briggs/DPU equilibrium"
    print()
    return rows


def briggs_rcf_comparison():
    """Claim 4: the slope matches Briggs; the coefficient does not, and why."""
    print("=" * 78)
    print("4. MECHANISTIC B_root vs BRIGGS RCF")
    print("=" * 78)
    theta = 0.90
    f_lip = api.TISSUE_TOTAL_LIPID_DW["root"]
    coef_model = (1.0 - theta) * f_lip * lp.BRIGGS_A
    print(f"{'logKow':>7} {'B_root':>10} {'RCF':>10} {'ratio':>7}")
    prev = None
    slopes = []
    for lk in (1.0, 2.0, 3.0, 4.0, 5.0):
        b = theta + (1.0 - theta) * f_lip * lp.briggs_klip(lk)
        rcf = lp.briggs_rcf(lk)
        print(f"{lk:7.1f} {b:10.3f} {rcf:10.2f} {rcf / b:7.2f}")
        if prev:
            slopes.append(np.log10((b - theta) / prev))
        prev = b - theta
    print(f"\n  Briggs lipid coefficient : 0.03      * K_ow^0.77   (FRESH weight, "
          f"'octanol equivalent')")
    print(f"  mechanistic coefficient  : {coef_model:.5f} * K_ow^0.77   "
          f"((1-theta) * f_lip_dw * a)")
    print(f"  ratio                    : {0.03 / coef_model:.1f}x")
    print(f"  K_ow slope, model        : {np.mean(slopes):.3f} per log unit "
          f"(Briggs b = {lp.BRIGGS_B})")
    need_dw = 0.03 / ((1.0 - theta) * lp.BRIGGS_A)
    print(f"\n  To match Briggs the root would need f_lip = {need_dw:.2f} dw "
          f"({need_dw * 100:.0f}% of dry weight as lipid) -- no root has that.")
    print(f"  Equivalently, Trapp's default L=0.025 is FRESH weight; at theta={theta} "
          f"that is {lp.f_lip_from_fresh_weight(0.025, theta):.2f} dw.")
    print("  => The STRUCTURE (K_ow slope) is reproduced; the empirical intercept is a")
    print("     lumped sorption capacity, not measured lipid.  Measured lipid is the")
    print("     default here; lp.LIPID_OCT_EQUIV_FW holds the Briggs-consistent anchor.")
    assert abs(np.mean(slopes) - lp.BRIGGS_B) < 1e-6, "the K_ow slope must match Briggs"
    print()


def main():
    speciation_limits()
    ion_trap()
    kow_sweep()
    briggs_rcf_comparison()
    print("=" * 78)
    print("ALL STRUCTURAL CHECKS PASSED")
    print("Reminder: structural verification against theory, NOT validation against")
    print("data -- there is no per-organ neutral time series for rice in the DB.")
    print("=" * 78)


if __name__ == "__main__":
    main()
