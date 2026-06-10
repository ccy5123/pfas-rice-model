"""
Tests for the literature-derived parameter module.

These lock in the *verified* QSPR relationships transcribed from the curated
literature database (docs/literature_db) -- the per-CF2 slopes, the sulfonate
offset, the measured Koc anchor, the f_d speciation, and the U-shaped protein
affinity -- and check that the builders produce objects that drop straight into
the plant/soil modules and preserve the structural root > straw > grain result.

Absolute placeholder intercepts (K_PL/K_prot/K_cw magnitude) are NOT asserted as
physical truths -- only their chain-length scaling, which IS from the literature.
"""
import numpy as np
import pytest

import literature_params as lp
from pfas_rice_plant_module import (
    Compound, Environment, Compartment, RiceUptakeModel, binding_factors,
    _logistic, ROOT, STEM, LEAF, FRUIT,
)
from soil_paddy import FreundlichSoil, inputs_from_soil


# ---------------------------------------------------------------------------
# C6 -- speciation: f_d
# ---------------------------------------------------------------------------
def test_fd_is_essentially_one_at_paddy_ph():
    # permanently-dissociated anion across the paddy porewater pH range
    for pH in (5.0, 6.0, 7.0):
        assert float(lp.f_d(lp.PKA["carboxylate"], pH)) > 0.999
        assert float(lp.f_d(lp.PKA["sulfonate"], pH)) > 0.9999


def test_fd_formula_and_monotonicity():
    assert float(lp.f_d(6.5, 6.5)) == pytest.approx(0.5)         # pKa == pH
    # robust even to the contested high PFOA pKa (Burns 2008)
    assert float(lp.f_d(3.8, 5.0)) >= 0.94
    # f_d decreases as pKa rises (less dissociation)
    assert float(lp.f_d(1.0, 6.0)) > float(lp.f_d(3.0, 6.0))


# ---------------------------------------------------------------------------
# C3 -- soil Koc QSPR
# ---------------------------------------------------------------------------
def test_koc_reproduces_measured_pfoa_anchor():
    _, npfc, hg = lp.SPECIES["PFOA"]
    assert float(lp.koc(npfc, hg)) == pytest.approx(96.0, rel=1e-9)


def test_koc_per_cf2_slope():
    # adjacent carboxylates differ by exactly the per-CF2 slope
    ratio = float(lp.koc(8, "carboxylate")) / float(lp.koc(7, "carboxylate"))
    assert ratio == pytest.approx(10.0 ** lp.KOC_PER_CF2)


def test_koc_sulfonate_offset():
    ratio = float(lp.koc(8, "sulfonate")) / float(lp.koc(8, "carboxylate"))
    assert ratio == pytest.approx(10.0 ** lp.KOC_SULFONATE_OFFSET)


def test_koc_ordering_and_log_option():
    assert float(lp.koc(8, "sulfonate")) > float(lp.koc(7, "carboxylate")) > float(lp.koc(4, "sulfonate"))
    assert float(lp.koc(7, "carboxylate", log10=True)) == pytest.approx(np.log10(96.0))


def test_koc_rejects_unknown_head_group():
    with pytest.raises(ValueError):
        lp.koc(7, "phosphonate")


def test_koc_to_KF_linear_is_exact():
    assert lp.koc_to_KF(96.0, 0.02, n=1.0) == pytest.approx(96.0 * 0.02)


# ---------------------------------------------------------------------------
# C4 -- binding factor QSPR
# ---------------------------------------------------------------------------
def test_kpl_anchor_and_slope():
    assert lp.k_pl(lp.KPL_ANCHOR_NPFC, "carboxylate") == pytest.approx(lp.KPL_ANCHOR_LKG)
    r_c = lp.k_pl(8, "carboxylate") / lp.k_pl(7, "carboxylate")
    r_s = lp.k_pl(8, "sulfonate") / lp.k_pl(7, "sulfonate")
    assert r_c == pytest.approx(10.0 ** lp.KPL_PER_CF2["carboxylate"])
    assert r_s == pytest.approx(10.0 ** lp.KPL_PER_CF2["sulfonate"])


def test_kprot_is_u_shaped_and_plant_weaker():
    # plateau optimum C6-C10 (equal); shorter/longer bind weaker
    plateau = [lp.k_prot(n, plant=False) for n in (6, 7, 8, 9, 10)]
    assert all(p == pytest.approx(plateau[0]) for p in plateau)
    assert lp.k_prot(4, plant=False) < plateau[0]
    assert lp.k_prot(12, plant=False) < plateau[0]
    # plant storage protein binds weaker than serum albumin (Zhou 2025)
    assert lp.k_prot(7, plant=True) == pytest.approx(lp.k_prot(7, plant=False) * lp.PLANT_PROTEIN_SCALE)
    assert lp.k_prot(7, plant=True) < lp.k_prot(7, plant=False)


# ---------------------------------------------------------------------------
# C4 -- MEASURED per-congener values (Chen 2025 SI / Zhou 2025 SI)
# ---------------------------------------------------------------------------
def test_kpl_uses_measured_chen_when_named():
    # named congener -> measured K_MW (10**logK_MW); PFOA logK_MW = 3.28
    assert lp.k_pl(name="PFOA") == pytest.approx(10.0 ** 3.28, rel=1e-6)
    assert lp.k_pl(name="PFOS") > lp.k_pl(name="PFOA") > lp.k_pl(name="PFBA")
    # unknown name falls back to the slope rule (needs nPFC)
    assert lp.k_pl(7, "carboxylate", name="not-a-pfas") == pytest.approx(lp.KPL_ANCHOR_LKG)


def test_kprot_albumin_from_hsa_kd():
    # K_prot = 1/(K_D[mol/L] * MW_HSA[kg/mol]); PFOA K_D = 2.57 umol/L
    expect = 1.0 / (2.57e-6 * lp.MW_HSA_KG_MOL)
    assert lp.k_prot_albumin("PFOA") == pytest.approx(expect, rel=1e-6)
    assert lp.k_prot_albumin("PFOS") > lp.k_prot_albumin("PFOA")   # PFOS binds HSA stronger
    assert lp.k_prot_albumin("not-a-pfas") is None
    # named k_prot uses the measured albumin value, plant-scaled
    assert lp.k_prot(name="PFOA", plant=False) == pytest.approx(expect, rel=1e-6)
    assert lp.k_prot(name="PFOA", plant=True) == pytest.approx(expect * lp.PLANT_PROTEIN_SCALE, rel=1e-6)


# ---------------------------------------------------------------------------
# C1 -- MEASURED calibration data (Kim et al. 2019)
# ---------------------------------------------------------------------------
def test_kim2019_grain_baf():
    baf = lp.kim2019_grain_baf("porewater")
    # PFOA: brown rice 0.349 ng/g, porewater 78.7 ng/L -> 0.349/(78.7/1000) = 4.43 L/kg
    assert baf["PFOA"] == pytest.approx(0.349 / (78.7 / 1000.0), rel=1e-6)
    assert all(v > 0 for v in baf.values())
    soil = lp.kim2019_grain_baf("soil")
    assert soil["PFOA"] == pytest.approx(0.349 / 0.160, rel=1e-6)


# ---------------------------------------------------------------------------
# C1 -- grain BAF chain-length trend
# ---------------------------------------------------------------------------
def test_grain_baf_factor_decreases_with_chain_length():
    assert lp.grain_baf_chain_factor(7) == pytest.approx(1.0)          # reference = PFOA
    assert lp.grain_baf_chain_factor(3) > lp.grain_baf_chain_factor(7) > lp.grain_baf_chain_factor(11)
    # -0.5 log per CF2 -> factor 0.1 per 2 CF2
    assert lp.grain_baf_chain_factor(9) / lp.grain_baf_chain_factor(7) == pytest.approx(0.1)


# ---------------------------------------------------------------------------
# species table
# ---------------------------------------------------------------------------
def test_species_info_known_and_unknown():
    assert lp.species_info("PFOA") == (8, 7, "carboxylate")
    assert lp.species_info("PFOS") == (8, 8, "sulfonate")
    with pytest.raises(KeyError):
        lp.species_info("not-a-pfas")


# ---------------------------------------------------------------------------
# C5 -- environment / membrane potential
# ---------------------------------------------------------------------------
def test_literature_environment_is_anion_excluding():
    env = lp.literature_environment()
    assert isinstance(env, Environment)
    lo, hi = lp.EM_RICE_ROOT_RANGE_V
    assert lo <= env.E <= hi              # within the rice-specific range (Wang 1994)
    assert env.N > 0                      # inside-negative membrane excludes the anion


# ---------------------------------------------------------------------------
# builders -> Compound / FreundlichSoil
# ---------------------------------------------------------------------------
def test_literature_compound_binding_scales_with_chain_length():
    cb, c8, cd = (lp.literature_compound("PFBA"),
                  lp.literature_compound("PFOA"),
                  lp.literature_compound("PFDA"))
    assert isinstance(c8, Compound)
    assert cd.K_PL > c8.K_PL > cb.K_PL    # longer chain -> stronger membrane binding
    assert c8.fd == pytest.approx(1.0, abs=1e-3)   # fully dissociated


def test_literature_compound_requires_identifiable_species():
    with pytest.raises(ValueError):
        lp.literature_compound("unknown-pfas")      # not in SPECIES, no nPFC/head_group


def test_literature_paddy_soil_prefers_measured_anchor():
    soil = lp.literature_paddy_soil("PFOA", f_oc=0.02, n=1.0)
    assert isinstance(soil, FreundlichSoil)
    assert soil.K_F == pytest.approx(lp.KOC_ANCHORS_LKG["PFOA"] * 0.02)
    # a species without a measured anchor falls back to the QSPR
    soil2 = lp.literature_paddy_soil("PFNA", f_oc=0.02, n=1.0)
    assert soil2.K_F == pytest.approx(float(lp.koc(8, "carboxylate")) * 0.02)


# ---------------------------------------------------------------------------
# end-to-end: the builders plug into the model and preserve the structure
# ---------------------------------------------------------------------------
def _full_model(species="PFOA"):
    t = np.linspace(0.0, 120.0, 481)
    C_total = np.full_like(t, 5.0)
    Qtp = 0.05 + 0.35 * np.exp(-((t - 75.0) ** 2) / (2 * 25.0 ** 2))
    M = np.column_stack([
        _logistic(t, 1e-3, 0.030, 0.10, 20.0), _logistic(t, 1e-3, 0.040, 0.10, 25.0),
        _logistic(t, 1e-3, 0.050, 0.12, 30.0), _logistic(t, 1e-5, 0.025, 0.18, 80.0)])
    soil = lp.literature_paddy_soil(species, f_oc=0.02, n=0.85)
    inputs = inputs_from_soil(t, C_total, Qtp, M, soil)
    comps = [Compartment("root",  0.70, 0.05, 0.010, 0.30),
             Compartment("stem",  0.80, 0.01, 0.005, 0.08),
             Compartment("leaf",  0.80, 0.03, 0.020, 0.04, S=20.0),
             Compartment("grain", 0.15, 0.08, 0.010, 0.10, S=2.0)]
    model = RiceUptakeModel(env=lp.literature_environment(),
                            cmpd=lp.literature_compound(species),
                            comps=comps, inputs=inputs)
    return t, model


def test_literature_parametrised_model_runs_and_orders_root_straw_grain():
    t, model = _full_model("PFOA")
    sol = model.solve(t)
    assert sol.success and np.all(np.isfinite(sol.y))
    C = sol.y[:, -1]
    Mf = model.inputs.M_(t[-1])
    straw = (C[STEM] * Mf[STEM] + C[LEAF] * Mf[LEAF]) / (Mf[STEM] + Mf[LEAF])
    assert C[ROOT] > straw > C[FRUIT]


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
