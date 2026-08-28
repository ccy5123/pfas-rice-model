#!/usr/bin/env python3
# =============================================================================
# validation/neutral_dpu_validation.py
# -----------------------------------------------------------------------------
# Tests the NEUTRAL-organic path (`src/neutral_dpu.py`) -- the Briggs/Kow base of
# the framework, which had been derived (docs/dpu_model_summary_corrected.tex) but
# never implemented, and never checked against anything.
#
# WHY THIS IS THE INTERESTING TEST. Every PFAS result in this repo is entangled
# with PFAS-specific parameters that had to be FIT (f_xy, k_seq, the lipid
# conductances), which is why the honest verdicts read "reproduction, not
# prediction". For a neutral compound there is nothing to fit: partitioning
# (K_PW) and root->shoot loading (TSCF) are both fixed from outside by published
# QSPRs on log Kow. So this is the one setting where the DPU BACKBONE -- four
# compartments, xylem advection, growth dilution, terminal accumulation -- is
# exposed on its own.
#
# WHAT THIS SCRIPT DOES AND DOES NOT ESTABLISH. Sections 1-4 are checks against
# published QSPRs and against the model's own structure; they can falsify the
# implementation and they quantify its scope, but they are NOT validation against
# measured plant data. Section 5 is the harness that does that -- and it is inert
# until a measured table is supplied (`data_obs/neutral_obs_template.csv` gives the
# schema). The environment this was written in blocks all publisher/PMC/Crossref
# access, so no measured neutral dataset could be obtained; see
# docs/neutral_dpu_validation.md for the shortlist of datasets to drop in.
#
#   python validation/neutral_dpu_validation.py
#   python validation/neutral_dpu_validation.py --obs data_obs/my_neutral_obs.csv
# =============================================================================
from __future__ import annotations
import csv, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT_DIR, "src"))

import forcing_rice as fr                                    # noqa: E402
from growth_rice import organ_biomass                        # noqa: E402
from pfas_rice_plant_module_4pool_surf import binding_factors  # noqa: E402
import neutral_dpu as ND                                     # noqa: E402

SEASON, NT = 120.0, 241
TISSUES = ("root", "stem", "leaf", "grain")


def drivers(season=SEASON, n_t=NT, Cwo=1.0):
    """Standard measured forcings (forcing_rice Q_TP + growth_rice organ biomass)."""
    t = np.linspace(0.0, float(season), int(n_t))
    b = organ_biomass(t, season)
    M = np.maximum(np.column_stack([b["root"], b["stem"], b["leaf"], b["grain"]]), 1e-4)
    return dict(t=t, Cwo=np.full(len(t), float(Cwo)), Qtp=fr.Q_TP(t, season), M=M)


# ---------------------------------------------------------------------------
# 1. the partition adapter reproduces Briggs exactly
# ---------------------------------------------------------------------------
def check_partitioning():
    print("=" * 84)
    print("1. PARTITION ADAPTER — does binding_factors() evaluate the Briggs K_PW?")
    print("=" * 84)
    print("   The core computes B_k = theta + (1-theta)*(f_prot*K_prot + f_PL*K_PL +")
    print("   f_cw*K_cw). With f_prot=f_cw=0 and K_PL=a*Kow^b that must collapse to")
    print("   Briggs' RCF = 0.82 + 10^(0.77 logKow - 1.52) on his barley root.")
    root = ND.briggs_root_compartment()
    print(f"\n{'logKow':>8}{'Briggs RCF':>13}{'model K_PW':>13}{'rel.err':>11}")
    worst = 0.0
    for lk in (-1.0, 0.0, 1.0, 1.78, 3.0, 4.5, 6.0):
        c = ND.neutral_compound(ND.NeutralCompound("probe", lk))
        model = float(binding_factors([root], c)[0])
        briggs = ND.briggs_rcf(lk)
        err = abs(model - briggs) / briggs
        worst = max(worst, err)
        print(f"{lk:>8.2f}{briggs:>13.4f}{model:>13.4f}{err:>11.2e}")
    ok = worst < 1e-12
    print(f"\n   max relative error {worst:.2e}  ->  {'PASS' if ok else 'FAIL'}")
    print("   (this is a test of the ADAPTER, not of Briggs: it shows the neutral")
    print("    path really is the published partition core and not a re-fit of it)")
    return ok


# ---------------------------------------------------------------------------
# 2. the Kow signature: does the backbone reproduce the known ordering?
# ---------------------------------------------------------------------------
def check_kow_signature(drv=None):
    print("\n" + "=" * 84)
    print("2. Kow SIGNATURE — shoot-dominated (polar) to root-retained (lipophilic)")
    print("=" * 84)
    print("   NOTHING is fitted here: K_PW and TSCF both come from log Kow. The")
    print("   qualitative law the model must reproduce is the one every uptake study")
    print("   reports -- polar compounds translocate to the shoot, lipophilic ones stay")
    print("   in the root -- with the turnover set by the Briggs bell (peak logKow 1.78).")
    drv = drv or drivers()
    print(f"\n{'logKow':>8}{'TSCF':>8}{'K_PW root':>11}{'root':>10}{'straw':>10}"
          f"{'grain':>10}{'straw/root':>12}")
    rows = []
    for lk in (-1.0, 0.0, 1.0, 1.78, 2.5, 3.5, 4.5, 5.5):
        r = ND.simulate_neutral(ND.NeutralCompound(f"logKow{lk}", lk), drv)
        tf = r["straw_baf"] / max(r["baf_final"]["root"], 1e-12)
        rows.append((lk, r["TSCF"], r["K_PW"]["root"], r["baf_final"]["root"],
                     r["straw_baf"], r["baf_final"]["grain"], tf))
        print(f"{lk:>8.2f}{r['TSCF']:>8.3f}{r['K_PW']['root']:>11.2f}"
              f"{r['baf_final']['root']:>10.2f}{r['straw_baf']:>10.2f}"
              f"{r['baf_final']['grain']:>10.2f}{tf:>12.2f}")
    tfs = [x[6] for x in rows]
    lks = [x[0] for x in rows]
    cross = next((lks[i] for i in range(1, len(tfs)) if tfs[i] < 1.0 <= tfs[i - 1]), None)
    peak_lk = lks[int(np.argmax(tfs))]
    monotone_tail = all(tfs[i] >= tfs[i + 1] for i in range(int(np.argmax(tfs)), len(tfs) - 1))
    print(f"\n   straw/root peaks at logKow {peak_lk:.2f} (Briggs TSCF peak 1.78);")
    print(f"   crosses 1 (root becomes the dominant compartment) near logKow "
          f"{cross if cross is not None else float('nan'):.2f};")
    print(f"   declines monotonically past the peak: {monotone_tail}")
    ok = monotone_tail and cross is not None and 1.0 <= peak_lk <= 3.0
    print(f"   -> {'PASS' if ok else 'FAIL'} (qualitative law reproduced with zero fitted parameters)")
    return ok, rows


# ---------------------------------------------------------------------------
# 3. scope: the leaf is UNBOUNDED without metabolism or volatilisation
# ---------------------------------------------------------------------------
def check_terminal_sink(drv=None):
    print("\n" + "=" * 84)
    print("3. SCOPE — for neutrals, metabolism/volatilisation are LOAD-BEARING")
    print("=" * 84)
    print("   With the phloem off (the neutral base has none) the leaf is the sole")
    print("   xylem terminal and its only sink is growth dilution, which -> 0 at")
    print("   maturity. So a recalcitrant, non-volatile neutral compound MUST run")
    print("   away: the leaf integrates the whole transpiration stream. That is")
    print("   arithmetic, not a bug -- but it means gamma=0 (fine for PFAS) is NOT a")
    print("   safe default here, and the air terms this module omits are not optional")
    print("   for a volatile compound. Quantifying it:")
    drv = drv or drivers()
    print(f"\n{'half-life (d)':>14}{'gamma (1/d)':>13}{'leaf BAF':>11}{'root BAF':>11}")
    base = None
    for hl in (np.inf, 60.0, 21.0, 7.0, 2.0):
        g = 0.0 if not np.isfinite(hl) else float(np.log(2.0) / hl)
        comps = ND.rice_compartments(gammas={k: g for k in TISSUES})
        r = ND.simulate_neutral(ND.NeutralCompound("probe", 1.78), drv, comps=comps)
        base = base if base is not None else r["baf_final"]["leaf"]
        lbl = "inf (recalcitrant)" if not np.isfinite(hl) else f"{hl:.0f}"
        print(f"{lbl:>14}{g:>13.4f}{r['baf_final']['leaf']:>11.2f}"
              f"{r['baf_final']['root']:>11.2f}")
    print("\n   => report neutral runs WITH a measured half-life, or state the result")
    print("      as an upper bound. `neutral_dpu.k_aw_warning` flags the volatile case.")
    return True


# ---------------------------------------------------------------------------
# 4. the ionic machinery really is off
# ---------------------------------------------------------------------------
def check_neutral_switch(drv=None):
    print("\n" + "=" * 84)
    print("4. SWITCH — z=0 removes anion exclusion and the carrier, exactly")
    print("=" * 84)
    drv = drv or drivers()
    r = ND.simulate_neutral(ND.NeutralCompound("probe", 2.0), drv)
    core = ND.neutral_compound(ND.NeutralCompound("probe", 2.0))
    ok = (abs(r["N"]) < 1e-12 and abs(r["eN"] - 1.0) < 1e-12
          and core.Vmax_in == 0.0 and core.Vmax_out == 0.0)
    print(f"   N = {r['N']:.3e} (PFAS: 4.67),  e^N = {r['eN']:.6f} (PFAS: ~107),")
    print(f"   Vmax_in = {core.Vmax_in}, Vmax_out = {core.Vmax_out}")
    print("   => the membrane term degenerates to passive Fickian diffusion, i.e. the")
    print("      neutral path is the SAME ODE with the physics of z=0, not a special case.")
    print(f"   -> {'PASS' if ok else 'FAIL'}")
    return ok


# ---------------------------------------------------------------------------
# 5. measured-data harness (inert until an observation table is supplied)
# ---------------------------------------------------------------------------
OBS_SCHEMA = ("compound", "log_kow", "tissue", "value", "basis", "endpoint")


def load_neutral_obs(path):
    """Read a measured neutral-organic table.

    Columns (see data_obs/neutral_obs_template.csv):
      compound, log_kow, tissue (root|stem|leaf|grain|straw), value,
      basis (fw|dw), endpoint (baf|conc|tf), and optionally exposure_ugL,
      half_life_d, K_AW, season_d, note.
    `basis=dw` values are converted to the model's fresh-weight basis with the
    tissue water contents actually used in the run, since (1-theta) differs by
    tissue and therefore does NOT cancel in a tissue/root ratio.
    """
    rows, skipped = [], 0
    with open(path, newline="") as f:
        for r in csv.DictReader(x for x in f if not x.lstrip().startswith("#")):
            if not r.get("compound"):
                continue
            missing = [c for c in OBS_SCHEMA if c not in r]
            if missing:
                raise ValueError(f"{path}: missing column(s) {missing}")
            # the template ships schema-example rows with fake values; refuse to
            # treat them as measurements, however the file is passed around
            if r["compound"].upper().startswith("EXAMPLE"):
                skipped += 1
                continue
            rows.append(dict(r, log_kow=float(r["log_kow"]), value=float(r["value"])))
    if skipped:
        print(f"   [skipped {skipped} EXAMPLE schema row(s) -- these are NOT data]")
    return rows


def compare_to_obs(path, drv=None):
    print("\n" + "=" * 84)
    print(f"5. MEASURED-DATA COMPARISON — {path}")
    print("=" * 84)
    obs = load_neutral_obs(path)
    if not obs:
        print("   (no rows) — supply a measured table to turn this into validation.")
        return None
    drv = drv or drivers()
    waters = ND.RICE_WATER
    by_cmpd = {}
    for r in obs:
        by_cmpd.setdefault((r["compound"], r["log_kow"]), []).append(r)

    pairs = []
    print(f"{'compound':16}{'logKow':>8}{'tissue':>8}{'obs':>10}{'model':>10}{'ratio':>9}")
    for (name, lk), rs in sorted(by_cmpd.items(), key=lambda x: x[0][1]):
        hl = next((float(r["half_life_d"]) for r in rs
                   if r.get("half_life_d") not in (None, "")), None)
        gam = float(np.log(2.0) / hl) if hl else 0.0
        comps = ND.rice_compartments(gammas={k: gam for k in TISSUES})
        m = ND.simulate_neutral(ND.NeutralCompound(name, lk), drv, comps=comps)
        for r in rs:
            tis, ep, basis = r["tissue"], r["endpoint"], r["basis"]
            if tis == "straw":
                pred = m["straw_baf"] if ep != "tf" else m["straw_baf"] / m["baf_final"]["root"]
            elif ep == "tf":
                pred = m["tf_final"][tis]
            else:
                pred = m["baf_final"][tis]
            o = r["value"]
            if basis == "dw":     # dw -> fw on the model's own water contents
                if ep == "tf":
                    wt = waters["grain"] if tis == "grain" else waters.get(tis, 0.80)
                    o = o * (1.0 - wt) / (1.0 - waters["root"])
                else:
                    wt = waters.get(tis, 0.80)
                    o = o * (1.0 - wt)
            pairs.append((pred, o))
            ratio = pred / o if o > 0 else float("inf")
            print(f"{name:16}{lk:>8.2f}{tis:>8}{o:>10.3f}{pred:>10.3f}"
                  + (f"{ratio:>9.2f}" if np.isfinite(ratio) else f"{'n/a':>9}"))
    e = [(np.log10(max(p, 1e-6)) - np.log10(max(o, 1e-6))) ** 2 for p, o in pairs]
    rmse = float(np.sqrt(np.mean(e)))
    print(f"\n   log10 RMSE (n={len(pairs)}) = {rmse:.3f}")
    print("   NOTE: zero parameters were fitted to this table -- K_PW and TSCF come")
    print("   from log Kow alone, so this is a genuine a-priori prediction, and it is")
    print("   the only such test in this repo. Interpret it against the PFAS side's")
    print("   a-priori error (log10 RMSE ~0.84-0.95, CLAUDE.md section 6).")
    return rmse


def main(obs_path=None):
    drv = drivers()
    print("NEUTRAL-ORGANIC DPU PATH — validation")
    print(f"forcings: measured Q_TP (peak {drv['Qtp'].max():.3f} L/d/hill), "
          f"growth_rice biomass, season {SEASON:.0f} d, Cwo = 1 ug/L\n")
    ok1 = check_partitioning()
    ok2, _ = check_kow_signature(drv)
    check_terminal_sink(drv)
    ok4 = check_neutral_switch(drv)
    rmse = compare_to_obs(obs_path, drv) if obs_path else None

    print("\n" + "=" * 84)
    print("STATUS")
    print("=" * 84)
    print(f"  structural checks: partition {'PASS' if ok1 else 'FAIL'}, "
          f"Kow signature {'PASS' if ok2 else 'FAIL'}, neutral switch "
          f"{'PASS' if ok4 else 'FAIL'}")
    if rmse is None:
        print("  measured-data comparison: NOT RUN — no observation table supplied.")
        print("  => The neutral path is IMPLEMENTED and internally consistent with the")
        print("     published Briggs QSPRs, and it reproduces the qualitative Kow law")
        print("     with zero fitted parameters. It is NOT YET VALIDATED against measured")
        print("     plant data. That step needs a measured table (schema:")
        print("     data_obs/neutral_obs_template.csv); candidate datasets and the access")
        print("     problem are documented in docs/neutral_dpu_validation.md.")
    else:
        print(f"  measured-data comparison: log10 RMSE {rmse:.3f} (a-priori, nothing fitted)")
    return ok1 and ok2 and ok4


if __name__ == "__main__":
    p = None
    if "--obs" in sys.argv:
        p = sys.argv[sys.argv.index("--obs") + 1]
    main(p)
