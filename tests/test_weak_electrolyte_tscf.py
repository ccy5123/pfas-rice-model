"""Guards for the first EMPIRICAL test of the weak-electrolyte path.

`tests/test_weak_electrolyte.py` pins the port STRUCTURALLY (bit-identity,
continuity, acid/base asymmetry). This file pins what the DATA said, in
`validation/weak_electrolyte_tscf.py`: direction supported, magnitude refuted.

The trap guard is the important one. Reading the table's pKa column instead of
its logD column manufactures a spectacular false counterexample out of the eight
pKa-1.62 barley rows, which are NOT ionised at the test pH. A later session that
"discovers" it is going backwards, so it is pinned here rather than only warned
about in a comment.
"""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "validation"))

import weak_electrolyte_tscf as WE                             # noqa: E402


@pytest.fixture(scope="module")
def rows():
    return WE.load()


# ---------------------------------------------------------------------------
# the data and how f_n is obtained
# ---------------------------------------------------------------------------
def test_the_table_splits_into_the_held_out_ionisable_rows(rows):
    """30 un-ionised rows are what every prior score used; 67 were held out."""
    assert len(rows) == 97
    assert sum(r["neutral"] for r in rows) == 30
    assert sum(not r["neutral"] for r in rows) == 67


def test_fn_comes_from_logD_and_never_exceeds_one(rows):
    """f_n = 10**(logD - logP), clipped. The clip absorbs the 20 rows where the
    two source columns disagree in the impossible direction (logD > logP)."""
    for r in rows:
        assert 0.0 < r["fn"] <= 1.0
        raw = 10.0 ** (r["logD"] - r["logP"])
        assert r["fn"] == pytest.approx(min(1.0, raw), rel=1e-12)


def test_the_pKa_162_rows_are_NOT_ionisable_the_false_counterexample(rows):
    """THE TRAP. pKa 1.62 read as an acid is fully dissociated at any test pH,
    and those eight barley rows have TSCF 0.63-0.98 -- which looks like a
    devastating refutation and is not one. Their own logD sits ABOVE their logP,
    so the table says they are un-ionised (a basic centre), and the shipped
    `neutral_at_test_pH` flag already classes them neutral."""
    fam = [r for r in rows if r["pKa"] == 1.62]
    assert len(fam) == 8
    assert all(r["neutral"] for r in fam), "pKa 1.62 rows must stay OUT of the test set"
    assert all(r["fn"] == 1.0 for r in fam)
    assert max(r["TSCF"] for r in fam) > 0.9        # the tempting-looking part


def test_the_pKa_column_and_the_logD_column_disagree(rows):
    """Why f_n is not taken from the pKa: Henderson-Hasselbalch on the pKa column
    does not reproduce the table's own logD, so the two columns cannot both be
    describing the same ionisation. Pinned so the derivation is not 'simplified'
    back onto the pKa."""
    import literature_params as LP
    resid = []
    for r in rows:
        if r["pKa"] == WE.MISSING:
            continue
        fn = LP.speciation(r["pKa"], r["pH"], True)[0]
        resid.append(abs(r["logD"] - (r["logP"] + np.log10(max(fn, 1e-12)))))
    assert np.median(resid) > 1.0        # recorded as 1.41 on the acid reading


# ---------------------------------------------------------------------------
# what the port predicts
# ---------------------------------------------------------------------------
def test_influx_ratio_is_exactly_one_in_the_neutral_limit():
    """The continuity property the whole port rests on: an un-ionised compound
    must feel no speciation penalty at all."""
    assert WE.influx_ratio(1.0, True) == pytest.approx(1.0, rel=1e-12)
    assert WE.influx_ratio(1.0, False) == pytest.approx(1.0, rel=1e-12)


def test_influx_ratio_spans_orders_of_magnitude_over_this_table(rows):
    """The magnitude claim being tested: Phi moves ~4 orders of magnitude across
    the f_n range these 67 rows cover."""
    fn = [r["fn"] for r in rows if not r["neutral"]]
    span = WE.influx_ratio(1.0) / WE.influx_ratio(min(fn))
    assert span > 1e3


def test_a_cation_is_less_penalised_than_an_anion():
    """The port's asymmetry, at the level this script uses it: the inside-negative
    membrane excludes the anion and attracts the cation."""
    assert WE.influx_ratio(1e-4, False) > WE.influx_ratio(1e-4, True)


def test_pKa_for_inverts_the_speciation(rows):
    import literature_params as LP
    for fn in (0.9, 0.1, 1e-3, 1e-5):
        for is_acid in (True, False):
            pKa = WE.pKa_for(fn, 6.5, is_acid)
            assert LP.speciation(pKa, 6.5, is_acid)[0] == pytest.approx(fn, rel=1e-9)


# ---------------------------------------------------------------------------
# what the data said
# ---------------------------------------------------------------------------
def test_direction_is_supported_measured_transfer_falls_as_it_ionises(rows):
    """The positive half of the verdict, and the port's first empirical support:
    measured TSCF really does rise with the neutral fraction."""
    ion = [r for r in rows if not r["neutral"]]
    rho = WE.spearman([r["fn"] for r in ion], [r["TSCF"] for r in ion])
    assert rho > 0.4, f"direction result lost (rho={rho:.3f})"


def test_magnitude_is_refuted_strongly_ionised_rows_still_transfer(rows):
    """The negative half: where the model predicts effectively nothing, the data
    show a mean TSCF around 0.13."""
    deep = [r for r in rows if not r["neutral"] and r["fn"] < 1e-3]
    assert len(deep) >= 10
    assert np.mean([r["TSCF"] for r in deep]) > 0.05
    # the model's own influx ratio there is ~1e-3 or below
    assert max(WE.influx_ratio(r["fn"]) for r in deep) < 1e-2


@pytest.mark.parametrize("is_acid", [True, False])
def test_speciation_on_orders_better_and_under_delivers(rows, is_acid):
    """The load-bearing half of the two-metric result, pinned on a deterministic
    subsample so the guard stays fast.

    ONLY the robust claims are asserted. The RMSE degradation is deliberately NOT
    asserted here: it holds on the full table (0.272 -> 0.304) but survives only
    ~82% of bootstrap resamples, and this very subsample flips it. Writing it as
    a test would pin an artefact -- see `bootstrap_wins` and section 4 of the
    validation script, which is where that comparison belongs."""
    ion = [r for r in rows if not r["neutral"]]
    sub = ion[::6]                          # 12 rows, stable under row order
    obs = np.array([r["TSCF"] for r in sub])
    off = np.array([WE.model_tscf(r["logP"], 1.0, r["pH"]) for r in sub])
    on = np.array([WE.model_tscf(r["logP"], r["fn"], r["pH"], is_acid) for r in sub])

    assert WE.spearman(on, obs) > WE.spearman(off, obs)          # better ordering
    assert np.mean(on - obs) < -0.1                              # and under-delivers


# ---------------------------------------------------------------------------
# the apoplastic bypass the result motivated
# ---------------------------------------------------------------------------
def test_g_apo_defaults_to_zero_everywhere_it_could_leak():
    """The bypass must be structurally absent unless asked for, on both paths."""
    import literature_params as LP
    import neutral_dpu as ND
    assert LP.literature_compound("PFOA").g_apo == 0.0
    assert ND.NeutralCompound("x", 2.45).g_apo == 0.0
    assert ND.neutral_compound(ND.NeutralCompound("x", 2.45)).g_apo == 0.0
    assert ND.neutral_compound(ND.NeutralCompound("x", 2.45, pKa=4.0)).g_apo == 0.0


def test_g_apo_is_gated_by_neither_speciation_nor_the_membrane_potential():
    """What makes it a separate term rather than a larger kappa_d. The bypass
    contribution to j_R must be identical for a neutral and for a fully ionised
    compound -- a route around the membrane cannot feel (fn, fd) or the GHK
    factor. This is the property the whole repair rests on."""
    from dataclasses import replace
    import literature_params as LP
    import pfas_rice_plant_module_4pool_surf as P4

    env = LP.literature_environment()
    base = LP.literature_compound("PFOA")
    Cwo, Cw = 1.0, 0.01
    contribs = []
    for fn, fd in ((1.0, 0.0), (0.5, 0.5), (0.0, 1.0)):
        off = replace(base, fn=fn, fd=fd, g_apo=0.0)
        on = replace(base, fn=fn, fd=fd, g_apo=0.7)
        contribs.append(P4.root_uptake(Cwo, Cw, on, env)
                        - P4.root_uptake(Cwo, Cw, off, env))
    assert contribs[0] == pytest.approx(contribs[1], rel=1e-12)
    assert contribs[1] == pytest.approx(contribs[2], rel=1e-12)
    assert contribs[0] == pytest.approx(0.7 * (Cwo - Cw), rel=1e-12)


def test_g_apo_is_self_targeting():
    """Why it does not simply inflate everything: uptake conductance sets how fast
    the root equilibrates, not the level it equilibrates to, so the bypass does
    almost nothing where the membrane is already fast (an un-ionised compound) and
    a great deal where speciation has collapsed it. That asymmetry is what lets one
    parameter address the ionisable rows without disturbing the 30 neutral ones."""
    lk, pH, g = 2.28, 6.5, 2.0
    neutral_off = WE.model_tscf(lk, 1.0, pH)
    neutral_on = WE.model_tscf(lk, 1.0, pH, g_apo=g)
    ion_off = WE.model_tscf(lk, 1e-4, pH)
    ion_on = WE.model_tscf(lk, 1e-4, pH, g_apo=g)

    assert neutral_on / neutral_off < 1.2, "bypass must barely move an un-ionised compound"
    assert ion_on / ion_off > 100.0, "bypass must rescue a strongly ionised one"


def test_a_small_bypass_helps_and_a_large_one_absorbs(rows):
    """The shape of the g_apo curve IS the result (docs section 4m), so it is the
    thing to pin: the ordering peaks at a SMALL bypass and the RMSE optimum lies
    past that peak. A later session that fits g_apo on RMSE and reports the 0.245
    without the +0.450 beside it has reintroduced exactly the absorption this was
    pre-registered against.

    Run on a deterministic subsample to stay fast; the full-table curve is in the
    docs and reproduced by validation/weak_electrolyte_tscf.py section 5.
    """
    ion = [r for r in rows if not r["neutral"]][::5]        # 14 rows
    obs = np.array([r["TSCF"] for r in ion])

    def preds(g):
        return np.array([WE.model_tscf(r["logP"], r["fn"], r["pH"], True, g_apo=g)
                         for r in ion])

    p0, p_small, p_big = preds(0.0), preds(0.5), preds(5.0)
    rmse = lambda p: float(np.sqrt(np.mean((p - obs) ** 2)))    # noqa: E731

    # a small bypass helps BOTH ways -- the Pareto point the docs quote
    assert WE.spearman(p_small, obs) > WE.spearman(p0, obs)
    assert rmse(p_small) < rmse(p0)
    # the large one buys RMSE by giving up the ordering
    assert rmse(p_big) < rmse(p_small)
    assert WE.spearman(p_big, obs) < WE.spearman(p_small, obs)


def test_the_rank_gain_is_robust_and_the_rmse_loss_is_not():
    """Pins the asymmetry itself, on synthetic arrays so it costs no ODE solves:
    `bootstrap_wins` must be able to report the two frequencies separately. The
    measured values (rank 0.94, RMSE 0.18) are recorded in the docs; what is
    pinned here is that the comparison is made at all and is direction-aware."""
    obs = np.linspace(0.05, 0.9, 40)
    off = np.random.default_rng(7).permutation(obs)  # right scale, ordering lost
    on = 0.01 * obs                                  # right ordering, wrong scale
    b = WE.bootstrap_wins(off, on, obs, n=200, seed=1)
    assert b["rank"] > 0.9                           # ON orders better
    assert b["rmse"] < 0.1                           # ON scales worse
