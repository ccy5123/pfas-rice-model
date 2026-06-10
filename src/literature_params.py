"""
Literature-derived parameters for the PFAS rice uptake model
============================================================

This module turns the curated literature database
(``docs/literature_db/PFAS_rice_parameter_database.xlsx`` and the per-sheet CSV
exports next to it) into machine-usable parameters and *builders* that plug
straight into :mod:`pfas_rice_plant_module` and :mod:`soil_paddy`.

Design / honesty notes
----------------------
The database supplies, for each model term, **verified QSPR slopes** and, in a
few cases, **verified absolute anchors**.  We encode them faithfully and flag
what is solid vs. what is still a placeholder:

* FULLY literature-grounded (verified DOIs):
    - ``f_d`` from pKa (C6; Goss 2008) -- the permanently-anionic assumption.
    - soil ``Koc`` chain-length QSPR (C3; Higgins & Luthy 2006 slope
      +0.50..0.60 log/CF2 and +0.23 sulfonate, anchored on the Milinovic 2015
      measured PFOA/PFOS/PFBS Koc values).
    - root membrane potential ``E_m`` for the GHK term (C5; Wang & Glass 1994,
      rice -116..-140 mV -- DOI not confirmed this session, see ``DOI_status``).
    - QSPR *slopes* for the membrane/protein binding that feeds ``B_k``
      (C4; Chen 2025 K_MW +0.36/CF2 (PFCA), +0.37/CF2 (PFSA)).

* PLACEHOLDER intercepts (flagged ``GAP`` in the database):
    - absolute ``K_PL``, ``K_prot``, ``K_cw`` per congener: the database gives
      the chain-length *slopes/shape* but the absolute values must be read from
      the cited SI tables (not done this session) or fitted.  We therefore scale
      a NOMINAL anchor (kept equal to the illustrative demo value) by the
      verified slope, so the model runs with chain-length-resolved binding while
      remaining explicit that the intercept is uncalibrated.

* NOT BAF-identifiable (fitted, not from the database):
    - ``f_xy`` (TSCF), ``L_Ph`` (phloem loading), ``kappa_d``, ``Vmax/Km``.
      The database confirms the *mechanisms* (carrier M-M, low anion TSCF,
      phloem mobility) and the chain-length *ordering*, but the numeric values
      are Tier-1/2 calibration targets.  Builders accept them as arguments with
      documented defaults matching the demo.

Symbol map (report / code -> this module)
-----------------------------------------
    f_d   (speciation)         -> f_d()
    Koc   (soil sub-model)     -> koc(), koc_to_KF(), literature_paddy_soil()
    K_PL  (binding, B_k)       -> k_pl()
    K_prot(binding, B_k)       -> k_prot()
    E, N  (GHK root uptake)    -> literature_environment()
    Compound / Environment     -> literature_compound(), literature_environment()

``DOI_status`` mirrors the database: ``verified`` = DOI seen verbatim during the
literature search; ``UNVERIFIED`` = lead not confirmed -- confirm before citing.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pfas_rice_plant_module import Compound, Environment

# ---------------------------------------------------------------------------
# references (mirrors docs/literature_db Source_Shortlist; status per database)
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Citation:
    key: str
    authors: str
    year: int
    journal: str
    doi: str
    doi_status: str          # "verified" | "UNVERIFIED"
    category: str            # database tab (C1..C6)


REFERENCES: dict[str, Citation] = {
    "goss2008": Citation("goss2008", "Goss", 2008, "Environ. Sci. Technol. 42:456-458",
                         "10.1021/es702192c", "verified", "C6"),
    "cheng2009": Citation("cheng2009", "Cheng et al.", 2009, "J. Phys. Chem. A 113:8152-8156",
                          "10.1021/jp9051352", "verified", "C6"),
    "torralba2023": Citation("torralba2023", "Torralba-Sanchez et al.", 2023,
                             "Environ. Toxicol. Chem. 42:2317", "10.1002/etc.5716", "verified", "C6"),
    "ebert2020": Citation("ebert2020", "Ebert, Allendorf, Berger, Goss, Ulrich", 2020,
                          "Environ. Sci. Technol. 54:5051-5061", "10.1021/acs.est.0c00175",
                          "verified", "C6/C4"),
    "higgins2006": Citation("higgins2006", "Higgins & Luthy", 2006,
                            "Environ. Sci. Technol. 40:7251-7256", "10.1021/es061000n",
                            "verified", "C3"),
    "milinovic2015": Citation("milinovic2015", "Milinovic, Lacorte, Vidal, Rigol", 2015,
                              "Sci. Total Environ. 511:63-71", "10.1016/j.scitotenv.2014.12.017",
                              "verified", "C3"),
    "fabregat2021": Citation("fabregat2021", "Fabregat-Palau, Vidal, Rigol", 2021,
                             "Sci. Total Environ. 801:149343", "10.1016/j.scitotenv.2021.149343",
                             "verified", "C3"),
    "jakobsen2026": Citation("jakobsen2026", "Jakobsen et al.", 2026, "Vadose Zone J.",
                             "10.1002/vzj2.70114", "verified", "C3"),
    "chen2025": Citation("chen2025", "Chen et al.", 2025, "Environ. Sci. Technol.",
                         "10.1021/acs.est.4c06734", "verified", "C4"),
    "droge2019": Citation("droge2019", "Droge", 2019, "Environ. Sci. Technol.",
                          "10.1021/acs.est.8b05052", "verified", "C4"),
    "zhou2025": Citation("zhou2025", "Zhou, Liu, Yuan, Ruan, Chen", 2025,
                         "Ecotox. Environ. Saf. 291:117902", "10.1016/j.ecoenv.2025.117902",
                         "verified", "C4"),
    "orosea2025": Citation("orosea2025", "O. rosea hyperaccumulator study", 2025,
                           "Nat. Commun.", "10.1038/s41467-025-65191-3", "verified", "C4"),
    "zhang2019": Citation("zhang2019", "Zhang, Sun et al.", 2019, "Sci. Total Environ. 654:19-27",
                          "10.1016/j.scitotenv.2018.10.443", "verified", "C5"),
    "wang1994": Citation("wang1994", "Wang, Glass et al.", 1994, "Plant Physiol. 104(3):899",
                         "10.1104/pp.104.3.899", "UNVERIFIED", "C5"),
    "liu2019": Citation("liu2019", "Liu, Lu, Song et al.", 2019, "Environ. Int. 127:671-684",
                        "10.1016/j.envint.2019.04.008", "verified", "C1"),
    "kim2019": Citation("kim2019", "Kim, Ekpe, Lee, Kim, Oh", 2019, "Sci. Total Environ. 671:714-721",
                        "10.1016/j.scitotenv.2019.03.240", "verified", "C1"),
}


# ---------------------------------------------------------------------------
# species table  (chain-length proxy + acid head group)
#   n_perfluoroC = number of perfluorinated carbons (database convention:
#   PFOA = C8 carboxylate has 7 perfluoro-C; PFOS = C8 sulfonate has 8).
# ---------------------------------------------------------------------------
SPECIES: dict[str, tuple[int, int, str]] = {
    # name      (total_C, n_perfluoroC, head_group)
    "PFBA":   (4, 3, "carboxylate"),
    "PFPeA":  (5, 4, "carboxylate"),
    "PFHxA":  (6, 5, "carboxylate"),
    "PFHpA":  (7, 6, "carboxylate"),
    "PFOA":   (8, 7, "carboxylate"),
    "PFNA":   (9, 8, "carboxylate"),
    "PFDA":   (10, 9, "carboxylate"),
    "PFUnDA": (11, 10, "carboxylate"),
    "PFDoDA": (12, 11, "carboxylate"),
    "PFBS":   (4, 4, "sulfonate"),
    "PFHxS":  (6, 6, "sulfonate"),
    "PFOS":   (8, 8, "sulfonate"),
}


def species_info(name: str) -> tuple[int, int, str]:
    """(total_C, n_perfluoroC, head_group) for a known PFAS name."""
    try:
        return SPECIES[name]
    except KeyError:
        raise KeyError(f"unknown PFAS {name!r}; known: {sorted(SPECIES)}") from None


# ===========================================================================
# C6 -- speciation: pKa -> f_d   (justifies the permanently-anionic anion)
# ===========================================================================
# PFCA pKa ~ -0.5..1 (Goss 2008; consensus near 0); PFSA pKa < 0 (strong acid).
# f_d is robust to the pKa controversy: even pKa=3.8 gives f_d>=0.94 at pH 5-7.
PKA: dict[str, float] = {"carboxylate": 0.5, "sulfonate": -3.0}   # goss2008 / torralba2023
PADDY_PH = 6.5                                                    # paddy porewater (range 5-7)


def f_d(pKa: float, pH: float = PADDY_PH):
    """Fraction dissociated (anion)  f_d = 1 / (1 + 10**(pKa - pH)).  [C6, verified]"""
    return 1.0 / (1.0 + np.power(10.0, np.asarray(pKa, float) - np.asarray(pH, float)))


# ===========================================================================
# C3 -- soil sorption: Koc(chain length, head group)
# ===========================================================================
KOC_PER_CF2 = 0.55             # log10 L/kg per CF2 (higgins2006; range 0.50-0.60)  [verified]
KOC_SULFONATE_OFFSET = 0.23    # log10 L/kg, sulfonate vs CF2-matched carboxylate    [verified]
# measured Koc anchors (milinovic2015), L/kg == mL/g  [verified]
KOC_ANCHORS_LKG: dict[str, float] = {"PFBS": 17.0, "PFOA": 96.0, "PFOS": 710.0}
_KOC_QSPR_ANCHOR = ("PFOA", 7, 96.0)   # (name, n_perfluoroC, Koc) carboxylate anchor


def koc(n_perfluoroC: float, head_group: str = "carboxylate", *, log10: bool = False):
    """Predicted organic-carbon-normalized partition Koc [L/kg].  [C3]

    QSPR anchored on the measured PFOA Koc (96 L/kg, milinovic2015) with the
    Higgins & Luthy per-CF2 slope and sulfonate offset:

        log Koc = log10(96) + 0.55*(nPFC - 7) + (0.23 if sulfonate else 0)

    Caveat (jakobsen2026): pure Koc UNDER-predicts short-chain (C4-C6) sorption;
    add a clay term there.  Where a measured anchor exists
    (``KOC_ANCHORS_LKG``) prefer it over this extrapolation.
    """
    _, anch_npfc, anch_koc = _KOC_QSPR_ANCHOR
    logK = np.log10(anch_koc) + KOC_PER_CF2 * (np.asarray(n_perfluoroC, float) - anch_npfc)
    if head_group == "sulfonate":
        logK = logK + KOC_SULFONATE_OFFSET
    elif head_group != "carboxylate":
        raise ValueError("head_group must be 'carboxylate' or 'sulfonate'")
    return logK if log10 else np.power(10.0, logK)


def koc_to_KF(Koc_LkG: float, f_oc: float, n: float = 1.0) -> float:
    """Freundlich capacity K_F from Koc and organic-carbon fraction f_oc.

    For the linear case (n=1) this is exact: K_d = Koc * f_oc [L/kg].  For
    n != 1 the Freundlich K_F carries concentration units that depend on n, so
    ``K_F ~= Koc * f_oc`` is a documented first-pass approximation (use a
    paddy-specific multi-point isotherm -- e.g. Qian/Wang 2020 -- to refine n).
    """
    return float(Koc_LkG) * float(f_oc)


# ===========================================================================
# C4 -- tissue binding factor B_k: K_PL, K_prot, K_cw
# ===========================================================================
# membrane(phospholipid)-water partition slope, per CF2 (chen2025)  [verified]
KPL_PER_CF2: dict[str, float] = {"carboxylate": 0.36, "sulfonate": 0.37}
KPL_PER_CF2_DROGE = 0.53       # alt PFSA slope (droge2019)                          [verified]
# NOMINAL anchors (absolute intercept = GAP; kept equal to the demo values so the
# model runs unchanged in magnitude while scaling correctly with chain length).
KPL_ANCHOR_LKG = 100.0         # K_PL at nPFC=7 carboxylate  [PLACEHOLDER intercept]
KPL_ANCHOR_NPFC = 7
KPROT_ANCHOR_LKG = 50.0        # K_prot serum-albumin scale  [PLACEHOLDER intercept]
KCW_ANCHOR_LKG = 20.0          # K_cw cell-wall              [PLACEHOLDER -- no coefficient in lit.]
# plant proteins (soy protein isolate, a grain storage-protein analog) bind
# WEAKER than serum albumin (zhou2025) -> grain B_k via storage protein is lower
# than animal-albumin estimates.  Ratio is illustrative (no precise value given).
PLANT_PROTEIN_SCALE = 0.3      # [PLACEHOLDER ratio < 1]


def k_pl(n_perfluoroC: float, head_group: str = "carboxylate",
         anchor: float = KPL_ANCHOR_LKG, anchor_npfc: int = KPL_ANCHOR_NPFC) -> float:
    """Phospholipid membrane-water partition K_PL [L/kg].  [C4]

    Scales a NOMINAL anchor by the verified per-CF2 slope (chen2025):
        log K_PL = log(anchor) + slope*(nPFC - anchor_npfc).
    """
    slope = KPL_PER_CF2[head_group]
    return float(anchor) * float(np.power(10.0, slope * (n_perfluoroC - anchor_npfc)))


def _kprot_chain_factor(n_perfluoroC: float) -> float:
    """U-shaped (inverted-V) chain-length factor for serum-albumin affinity.

    Affinity is NON-monotonic with an optimum plateau at nPFC 6-10 (chen2025;
    cross-species albumin study); shorter and longer chains bind weaker.  The
    plateau is from the literature; the off-plateau decline rate is illustrative.
    """
    if 6 <= n_perfluoroC <= 10:
        return 1.0
    d = (6 - n_perfluoroC) if n_perfluoroC < 6 else (n_perfluoroC - 10)
    return float(np.power(10.0, -0.25 * d))     # ~ -0.25 log/CF2 off the optimum [illustrative]


def k_prot(n_perfluoroC: float, plant: bool = True, anchor: float = KPROT_ANCHOR_LKG) -> float:
    """Protein-water partition K_prot [L/kg].  [C4]

    U-shaped chain-length factor (optimum C6-C10) on a NOMINAL anchor; when
    ``plant=True`` it is reduced by ``PLANT_PROTEIN_SCALE`` (plant storage
    proteins bind weaker than serum albumin, zhou2025).
    """
    val = float(anchor) * _kprot_chain_factor(n_perfluoroC)
    return val * PLANT_PROTEIN_SCALE if plant else val


# ===========================================================================
# C1/C5 -- chain-length translocation trend (diagnostic; informs f_xy ordering)
# ===========================================================================
GRAIN_BAF_PER_CF2 = -0.5       # log10 grain-BAF per CF2 (liu2019)  [verified trend]


def grain_baf_chain_factor(n_perfluoroC: float, ref_n_perfluoroC: int = 7) -> float:
    """Relative grain-BAF vs chain length (1.0 at the reference, PFOA).  [C1]

    log(grain BAF) falls ~0.5 per CF2 (short chains preferentially translocated;
    long chains retained in root).  Use to order ``f_xy`` across a chain-length
    series; it is a diagnostic scaling, not a direct model parameter.
    """
    return float(np.power(10.0, GRAIN_BAF_PER_CF2 * (n_perfluoroC - ref_n_perfluoroC)))


# ===========================================================================
# C5 -- root membrane potential for the GHK exclusion term
# ===========================================================================
EM_RICE_ROOT_V = -0.120                 # central value (wang1994: -0.116..-0.140 V)
EM_RICE_ROOT_RANGE_V = (-0.140, -0.116) # rice-specific bounds (wang1994)
EM_PLANT_GENERIC_RANGE_V = (-0.160, -0.090)  # generic plant sweep (sensitivity)
# NOTE: NH4+ (dominant N in flooded paddy) DEPOLARIZES the membrane -> less
# negative E_m -> weaker anion exclusion.  In-situ paddy E_m is a database GAP.


# ===========================================================================
# builders -- literature-parametrised model objects
# ===========================================================================
def literature_environment(E_V: float = EM_RICE_ROOT_V, T: float = 298.15,
                           z: int = -1) -> Environment:
    """:class:`Environment` with the rice root membrane potential (C5, wang1994)."""
    return Environment(T=T, E=E_V, z=z)


def literature_compound(name: str | None = None, *,
                        n_perfluoroC: int | None = None, head_group: str | None = None,
                        pH: float = PADDY_PH, plant_protein: bool = True,
                        # Tier-1/2 transport params: NOT in the database (fitted);
                        # defaults match the demo so the model runs out of the box.
                        f_xy: float = 0.02, L_Ph: float = 0.005, kappa_d: float = 0.5,
                        Vmax_in: float = 20.0, Km_in: float = 5.0,
                        Vmax_out: float = 8.0, Km_out: float = 5.0,
                        kpl_anchor: float = KPL_ANCHOR_LKG,
                        kprot_anchor: float = KPROT_ANCHOR_LKG,
                        kcw_anchor: float = KCW_ANCHOR_LKG) -> Compound:
    """Build a :class:`Compound` with literature-derived binding & speciation.

    Binding (``K_prot``, ``K_PL``) scales with chain length via the verified C4
    QSPR slopes; ``f_d`` comes from the C6 pKa.  ``K_cw`` and the transport
    parameters are placeholders / fitted (see module docstring).
    """
    if name is not None and name in SPECIES:
        _, n_perfluoroC, head_group = SPECIES[name]
    if n_perfluoroC is None or head_group is None:
        raise ValueError("provide a known `name` (in SPECIES) or both "
                         "`n_perfluoroC` and `head_group`")
    label = name or f"PF{n_perfluoroC}{'S' if head_group == 'sulfonate' else 'A'}"
    return Compound(
        name=label,
        K_prot=k_prot(n_perfluoroC, plant=plant_protein, anchor=kprot_anchor),
        K_PL=k_pl(n_perfluoroC, head_group, anchor=kpl_anchor),
        K_cw=kcw_anchor,
        kappa_d=kappa_d, Vmax_in=Vmax_in, Km_in=Km_in,
        Vmax_out=Vmax_out, Km_out=Km_out,
        L_Ph=L_Ph, f_xy=f_xy,
        fd=float(f_d(PKA[head_group], pH)), fn=0.0,
    )


def literature_paddy_soil(name: str | None = None, *,
                          n_perfluoroC: int | None = None, head_group: str | None = None,
                          f_oc: float = 0.02, n: float = 1.0, theta_g: float = 0.45,
                          Koc: float | None = None, redox: str = "aerobic"):
    """Build a :class:`soil_paddy.FreundlichSoil` from the C3 Koc QSPR.

    ``Koc`` precedence: explicit arg > measured anchor (``KOC_ANCHORS_LKG``) >
    QSPR ``koc()``.  ``K_F = Koc * f_oc`` (see :func:`koc_to_KF`).  Anoxic/flooded
    sorption is a database GAP, so for a redox switch wrap two states (one with a
    reduced ``f_oc``/``Koc``) in :class:`soil_paddy.PaddyRedox` and treat the
    flooded state as illustrative.
    """
    from soil_paddy import FreundlichSoil
    if Koc is None:
        if name is not None and name in KOC_ANCHORS_LKG:
            Koc = KOC_ANCHORS_LKG[name]
        else:
            if name is not None and name in SPECIES:
                _, n_perfluoroC, head_group = SPECIES[name]
            if n_perfluoroC is None or head_group is None:
                raise ValueError("need a measured-anchor `name`, an explicit `Koc`, "
                                 "or (n_perfluoroC, head_group) for the QSPR")
            Koc = float(koc(n_perfluoroC, head_group))
    K_F = koc_to_KF(Koc, f_oc, n)
    return FreundlichSoil(K_F=K_F, n=n, theta_g=theta_g,
                          name=f"{name or 'PFAS'}/{redox}")


# ===========================================================================
# demo: show the QSPRs and an end-to-end literature-parametrised run
# ===========================================================================
def _demo():
    from pfas_rice_plant_module import (
        RiceUptakeModel, Compartment, binding_factors, _logistic,
        ROOT, STEM, LEAF, FRUIT,
    )
    from soil_paddy import inputs_from_soil

    print("== C6  f_d (fraction dissociated) at paddy pH  [goss2008, verified] ==")
    for hg in ("carboxylate", "sulfonate"):
        print(f"  {hg:11s} pKa={PKA[hg]:+.1f}  f_d(pH5)={float(f_d(PKA[hg],5)):.5f}  "
              f"f_d(pH7)={float(f_d(PKA[hg],7)):.5f}")

    print("\n== C3  Koc QSPR vs measured anchors  [higgins2006 + milinovic2015] ==")
    for nm in ("PFBS", "PFBA", "PFHxA", "PFOA", "PFNA", "PFDA", "PFOS"):
        _, npfc, hg = SPECIES[nm]
        meas = KOC_ANCHORS_LKG.get(nm)
        ms = f"   measured={meas:6.1f}" if meas else ""
        print(f"  {nm:7s} nPFC={npfc:2d} {hg:11s} Koc_pred={float(koc(npfc,hg)):7.1f} L/kg{ms}")

    print("\n== C4  K_PL chain-length scaling  [chen2025 slope; intercept PLACEHOLDER] ==")
    for nm in ("PFBA", "PFHxA", "PFOA", "PFDA", "PFDoDA"):
        _, npfc, hg = SPECIES[nm]
        print(f"  {nm:7s} nPFC={npfc:2d}  K_PL={k_pl(npfc,hg):8.2f}  "
              f"K_prot(plant)={k_prot(npfc):7.2f} L/kg")

    print("\n== C1  grain-BAF chain-length factor (rel. PFOA)  [liu2019, verified trend] ==")
    for nm in ("PFBA", "PFHxA", "PFOA", "PFDA", "PFDoDA"):
        _, npfc, _ = SPECIES[nm]
        print(f"  {nm:7s} nPFC={npfc:2d}  grain-BAF factor = {grain_baf_chain_factor(npfc):.3f}")

    print("\n== end-to-end: literature-parametrised PFOA in paddy rice ==")
    t = np.linspace(0.0, 120.0, 481)
    C_total = np.full_like(t, 5.0)                       # ug/kg dry soil
    Qtp = 0.05 + 0.35 * np.exp(-((t - 75.0) ** 2) / (2 * 25.0 ** 2))
    M = np.column_stack([
        _logistic(t, 1e-3, 0.030, 0.10, 20.0), _logistic(t, 1e-3, 0.040, 0.10, 25.0),
        _logistic(t, 1e-3, 0.050, 0.12, 30.0), _logistic(t, 1e-5, 0.025, 0.18, 80.0)])
    soil = literature_paddy_soil("PFOA", f_oc=0.02, n=0.85, theta_g=0.45)
    inputs = inputs_from_soil(t, C_total, Qtp, M, soil)
    env = literature_environment()
    cmpd = literature_compound("PFOA")
    comps = [Compartment("root",  0.70, 0.05, 0.010, 0.30),
             Compartment("stem",  0.80, 0.01, 0.005, 0.08),
             Compartment("leaf",  0.80, 0.03, 0.020, 0.04, S=20.0),
             Compartment("grain", 0.15, 0.08, 0.010, 0.10, S=2.0)]
    model = RiceUptakeModel(env=env, cmpd=cmpd, comps=comps, inputs=inputs)
    sol = model.solve(t)
    Cend = sol.y[:, -1]
    B = binding_factors(comps, cmpd)
    print(f"  E_m={env.E*1e3:.0f} mV  N={env.N:.3f}  e^N={np.exp(env.N):.1f}  f_d={cmpd.fd:.5f}")
    print(f"  K_prot={cmpd.K_prot:.2f}  K_PL={cmpd.K_PL:.2f}  K_cw={cmpd.K_cw:.2f}  [L/kg]")
    print(f"  soil: Koc={KOC_ANCHORS_LKG['PFOA']:.0f} L/kg -> K_F={soil.K_F:.3f} "
          f"(f_oc=0.02), n={soil.n}; C_w^o(day45)={inputs.Cwo_(45):.4f} ug/L")
    print(f"  B_k [L/kg]: " + ", ".join(f"{c.name}={b:.2f}" for c, b in zip(comps, B)))
    Mf = inputs.M_(t[-1])
    straw = (Cend[STEM] * Mf[STEM] + Cend[LEAF] * Mf[LEAF]) / (Mf[STEM] + Mf[LEAF])
    print(f"  final conc [ug/kg]: root={Cend[ROOT]:.3f}  straw={straw:.3f}  grain={Cend[FRUIT]:.3f}")
    ok = Cend[ROOT] > straw > Cend[FRUIT]
    print(f"  ordering root > straw > grain: {'OK' if ok else 'VIOLATED'}")


if __name__ == "__main__":
    _demo()
