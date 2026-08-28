"""
Tests for the Hwang 2017 lettuce/chlorpyrifos check -- validation/hwang2017_lettuce.py.

This dataset is deliberately NOT wired into the shared `--obs` harness (see the
module header there): that harness runs every row on the rice drivers, and Hwang
is lettuce over 40 days under a DECAYING exposure, so a number from it would be
silently meaningless. The tests below therefore pin the dedicated script instead,
and in particular they pin the CAVEATS -- because the value of this dataset is
entirely in what it refuses to conclude, and a later session must not quietly
promote either RMSE into a validation claim.
"""
import os
import sys

import numpy as np
import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "src"))
sys.path.insert(0, os.path.join(_ROOT, "validation"))

import hwang2017_lettuce as HW  # noqa: E402
import neutral_dpu as ND  # noqa: E402
from pfas_rice_plant_module_4pool_surf import binding_factors  # noqa: E402


# --------------------------------------------------------------------------
# the transcribed data and the authors' own exposure model
# --------------------------------------------------------------------------
def test_exposure_is_the_authors_equation():
    """Ce(t) = C0 * (1/2)^(t/T) / Kd, in ug/L. Every input is measured except the
    soil half-life, which the authors flag as literature-derived."""
    for level in (10, 20):
        p = HW.SOIL[level]
        for day in (0.0, 21.0, 40.0):
            assert float(HW.exposure(level, day)) == pytest.approx(
                p["C0"] * 0.5 ** (day / p["T_soil"]) / HW.KD * 1e3, rel=1e-12)
    # exposure declines, and faster in the treatment with the shorter soil half-life
    assert HW.exposure(10, 40.0) < HW.exposure(10, 0.0)
    assert HW.exposure(20, 40.0) / HW.exposure(20, 0.0) < \
        HW.exposure(10, 40.0) / HW.exposure(10, 0.0)


def test_table1_is_internally_consistent_on_one_basis():
    """`whole` must be the mass-weighted mean of leaf and root at ONE root mass
    fraction, the same at every sampling. That is a real check on the
    transcription and on the columns sharing a basis -- not a fit."""
    x = []
    for (lvl, day), (lf, rt, wh) in HW.TABLE1.items():
        f = (wh - lf) / (rt - lf)
        if f > 0.005:                      # two rows are forced to 0 by 1-dp rounding
            x.append(f)
    assert len(x) == 4
    assert np.mean(x) == pytest.approx(0.054, abs=0.005)
    assert np.std(x) < 0.02                # tight: one mass split explains them all
    # and it is the value the model run actually uses
    assert HW.ROOT_MASS_FRACTION == pytest.approx(float(np.mean(x)), abs=0.005)


def test_growth_reconstruction_rejects_the_exponential_form():
    """`Ig`/`Kg` parameterise an equation the transcription does not carry. Only
    the log-log forms give lettuce-scale masses; the exponential reading is
    rejected on magnitude, and that rejection must stay visible."""
    for form in ("log10", "ln"):
        m40 = float(HW.growth(10, 40.0, form)) * 1e3      # g
        assert 50.0 < m40 < 200.0                          # a plausible lettuce head
        assert HW.growth(10, 40.0, form) > HW.growth(10, 21.0, form)
    assert float(HW.growth(10, 40.0, "exp")) * 1e3 > 1e6   # absurd -> rejected
    with pytest.raises(ValueError):
        HW.growth(10, 10.0, "nope")


# --------------------------------------------------------------------------
# the finding: the unstated basis spans the verdict
# --------------------------------------------------------------------------
def test_root_partition_ceiling_is_structural():
    """The modelled root cannot exceed its equilibrium partition K_PW -- that is
    what the compartment IS -- so no parameter choice can reach a measurement
    above it. Pinned because the whole section-2 argument rests on it."""
    comps = HW.lettuce_compartments()
    core = ND.neutral_compound(ND.NeutralCompound("chlorpyrifos", HW.LOG_KOW))
    kpw_root = float(binding_factors(comps, core)[0])
    assert kpw_root == pytest.approx(15.8, abs=0.3)
    res = HW.predict(10)
    assert res["conc"]["root"].max() / float(HW.exposure(10, 0.0)) < kpw_root * 1.01


def test_the_basis_flips_the_verdict():
    """Read fresh, the measured root is far ABOVE the model's ceiling; read dry it
    is BELOW. One unstated table footnote therefore decides whether this dataset
    refutes or confirms the partition core -- which is why A3 concludes nothing
    about it."""
    comps = HW.lettuce_compartments()
    core = ND.neutral_compound(ND.NeutralCompound("chlorpyrifos", HW.LOG_KOW))
    ceiling = float(binding_factors(comps, core)[0])
    fw = [HW.observed_baf(l, d, "root", "fw") for (l, d) in HW.TABLE1]
    dw = [HW.observed_baf(l, d, "root", "dw") for (l, d) in HW.TABLE1]
    assert min(fw) > 2.5 * ceiling          # unreachable on every sampling
    assert max(dw) <= 1.1 * ceiling         # bracketed, never far above
    assert all(f > d for f, d in zip(fw, dw))


def test_the_two_readings_fail_on_opposite_organs():
    """The sharpest result here: fresh weight predicts the leaf well and misses the
    root, dry weight brackets the root and over-predicts the leaf. So the basis
    decides WHERE the model is wrong, not WHETHER -- the discrepancy is not a
    units artifact that the right footnote would dissolve."""
    _, fw = HW.compare(basis="fw", quiet=True, by_tissue=True)
    _, dw = HW.compare(basis="dw", quiet=True, by_tissue=True)
    assert fw["leaf"] < fw["root"]          # fresh: leaf is the good organ
    assert dw["root"] < dw["leaf"]          # dry:   root is the good organ
    assert dw["root"] < fw["root"] and fw["leaf"] < dw["leaf"]


def test_soil_contact_cannot_explain_the_root_exceedance():
    """Adhering soil is a real confound, but bounded: soil at C0 is the most
    concentrated thing a root can carry, so the mass fraction that would have to
    be soil is a hard number. It is implausibly large for the 10 mg/kg rows."""
    phi = {k: rt / HW.SOIL[k[0]]["C0"] for k, (_, rt, _) in HW.TABLE1.items()}
    assert phi[(10, 21)] > 0.4              # ~49 % of washed root mass: not credible
    assert phi[(10, 30)] > 0.3
    assert all(0.0 < v < 1.0 for v in phi.values())


def test_half_life_scan_is_not_evidence_about_the_basis():
    """Metabolism can only LOWER modelled concentrations, so it helps the reading
    where the model is too high (dry) and hurts the one where it is too low
    (fresh). The fit therefore 'prefers' dry weight -- and using that to choose
    the basis would be circular. Pinned so the asymmetry stays on the record."""
    fw0 = HW.compare(None, basis="fw", quiet=True)
    dw0 = HW.compare(None, basis="dw", quiet=True)
    fw1 = HW.compare(8.4, basis="fw", quiet=True)
    dw1 = HW.compare(8.4, basis="dw", quiet=True)
    assert fw1 > fw0                        # metabolism makes the fresh reading worse
    assert dw1 < dw0                        # and the dry reading better
