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

* MEASURED per-congener values, extracted from the cited papers (this session;
  see ``docs/literature_db/raw_si/``):
    - ``K_PL`` (= membrane-water K_MW): Chen 2025 Table S5, log L/kg lipid,
      cross-checked against Droge 2019 SSLM (same method, within ~0.2 log).
    - ``K_prot`` (= protein-water K_prow): Zhou 2025 Table 1 dialysis values for
      soy protein isolate (the PLANT/grain storage-protein analog) and BSA
      (animal reference).  [The Chen 2025 HSA K_D is reference-only -- the
      single-site binding-constant route overestimates the partition ~50x.]
    - calibration target: Kim 2019 per-congener brown-rice (grain) BAF, paired
      with paddy pore water.

* STILL PLACEHOLDER / GAP:
    - ``K_cw`` (cell-wall): no partition coefficient exists in the literature
      (only the binding components -- pectin/hemicellulose -- are identified);
      kept as a nominal anchor.

* NOT BAF-identifiable (fitted, not from the database):
    - ``f_xy`` (TSCF), ``L_Ph`` (phloem loading), ``kappa_d``, ``Vmax/Km``.
      The database confirms the *mechanisms* and chain-length *ordering*, but the
      numeric values are Tier-1/2 calibration targets.  ``_demo()`` shows fitting
      ``L_Ph`` to the Kim 2019 PFOA grain BAF.

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
    "tang2026": Citation("tang2026", "Tang, Xiao, Wu, Wang, Ge, Zhu, Chu, Chen", 2026,
                         "J. Hazard. Mater.", "10.1016/j.jhazmat.2025.141017", "verified", "C1/C5"),
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
# --- MEASURED membrane-water partition K_MW (= K_PL), log10 L/kg lipid, pH 7.0 ---
# chen2025 Table S5 (SSLM/TRANSIL kit).  Cross-checked against droge2019 SSLM
# (same method) for the PFCAs -> agree within ~0.2 log, confirming L/kg lipid.
# (see docs/literature_db/raw_si/chen2025_kmw_hsa.csv, droge2019_kmw.csv)
KMW_CHEN2025_LOG: dict[str, float] = {
    "PFBA": 1.63, "PFPeA": 2.02, "PFHxA": 2.42, "PFHpA": 2.85, "PFOA": 3.28,
    "PFNA": 3.75, "PFDA": 4.18, "PFUnDA": 4.50, "PFDoDA": 4.82,
    "PFBS": 2.72, "PFHxS": 3.65, "PFOS": 4.50,
}
# --- MEASURED HSA binding: equilibrium dissociation constant K_D [umol/L], pH 7.4 ---
# chen2025 Table S5.  Converted to a protein-water partition by k_prot_albumin().
KD_HSA_CHEN2025_UMOL: dict[str, float] = {
    "PFBA": 84.9, "PFPeA": 16.27, "PFHxA": 13.98, "PFHpA": 6.2, "PFOA": 2.57,
    "PFNA": 3.75, "PFDA": 3.64, "PFUnDA": 4.73, "PFDoDA": 15.76,
    "PFBS": 15.6, "PFHxS": 1.31, "PFOS": 0.94,
}
MW_HSA_KG_MOL = 66.5           # human serum albumin molecular weight
# BSA association constants K_A [L/mol] @300 K (zhou2025 Table S4) -- a BINDING
# constant from fluorescence quenching (NOT a partition coefficient); see note below.
KA_BSA_ZHOU2025: dict[str, float] = {
    "PFHxA": 3.03e4, "PFOA": 2.726e5, "PFDA": 7.31e5,
    "PFBS": 47.6, "PFHxS": 2.64e4, "PFOS": 4.12e5,
}
# --- MEASURED protein-water partition K_prow [log10 L/kg], dialysis (zhou2025 Table 1) ---
# This is the PARTITION coefficient the model needs (C_protein/C_water at low conc),
# measured directly for FOUR proteins.  PREFERRED over the K_A/K_D route, which
# (single-site) overestimates by ~50x (BSA PFOA: dialysis 110 vs K_D-derived ~5850).
# "soy" = soy protein isolate = grain storage-protein analog -> the PLANT K_prot.
# (see docs/literature_db/raw_si/zhou2025_kprow.csv)
KPROW_ZHOU2025_LOG: dict[str, dict[str, float]] = {
    "bsa":         {"PFBA": 2.31, "PFHxA": 2.16, "PFOA": 2.04, "PFDA": 1.67,
                    "PFBS": 2.28, "PFHxS": 2.18, "PFOS": 1.88},
    "faf_bsa":     {"PFBA": 2.33, "PFHxA": 2.16, "PFOA": 2.00, "PFDA": 1.64,
                    "PFBS": 2.37, "PFHxS": 2.08, "PFOS": 1.94},
    "phycocyanin": {"PFBA": 1.08, "PFHxA": 1.11, "PFOA": 1.28, "PFDA": 1.54,
                    "PFBS": 1.10, "PFHxS": 1.68, "PFOS": 1.69},
    "soy":         {"PFBA": 1.03, "PFHxA": 1.15, "PFOA": 1.09, "PFDA": 1.37,
                    "PFBS": 1.56, "PFHxS": 1.63, "PFOS": 1.54},
}

# membrane(phospholipid)-water partition slope, per CF2 (chen2025)  [verified] --
# used only as a FALLBACK for congeners absent from the measured dict.
KPL_PER_CF2: dict[str, float] = {"carboxylate": 0.36, "sulfonate": 0.37}
KPL_PER_CF2_DROGE = 0.53       # alt PFSA slope (droge2019)                          [verified]
KPL_ANCHOR_LKG = 10.0 ** KMW_CHEN2025_LOG["PFOA"]   # measured PFOA anchor (~1905 L/kg)
KPL_ANCHOR_NPFC = 7
KPROT_ANCHOR_LKG = 12.0        # K_prot last-resort fallback [L/kg] (~ measured soy PFOA),
                               # used only if a congener is absent from KPROW_ZHOU2025_LOG
KCW_ANCHOR_LKG = 20.0          # K_cw cell-wall  [PLACEHOLDER -- no coefficient in literature]
# Plant storage proteins (soy protein isolate) bind PFAS WEAKER than serum
# albumin -- now captured DIRECTLY by the measured K_prow (zhou2025 Table 1), so
# NO scale factor is needed: plant/grain tissues use protein="soy", animal "bsa".


def k_pl(n_perfluoroC: float | None = None, head_group: str = "carboxylate",
         anchor: float = KPL_ANCHOR_LKG, anchor_npfc: int = KPL_ANCHOR_NPFC,
         *, name: str | None = None) -> float:
    """Phospholipid membrane-water partition K_PL [L/kg lipid].  [C4]

    Prefers the MEASURED chen2025 K_MW when ``name`` is a known congener; otherwise
    falls back to the per-CF2 slope from the measured PFOA anchor:
        log K_PL = log(anchor) + slope*(nPFC - anchor_npfc).
    """
    if name is not None and name in KMW_CHEN2025_LOG:
        return float(np.power(10.0, KMW_CHEN2025_LOG[name]))
    if n_perfluoroC is None:
        raise ValueError("k_pl needs a measured `name` or an `n_perfluoroC` (fallback)")
    slope = KPL_PER_CF2[head_group]
    return float(anchor) * float(np.power(10.0, slope * (n_perfluoroC - anchor_npfc)))


def _kprot_chain_factor(n_perfluoroC: float) -> float:
    """U-shaped (inverted-V) fallback chain-length factor (optimum nPFC 6-10).
    Used ONLY when no measured K_prow exists for the congener or its family."""
    if 6 <= n_perfluoroC <= 10:
        return 1.0
    d = (6 - n_perfluoroC) if n_perfluoroC < 6 else (n_perfluoroC - 10)
    return float(np.power(10.0, -0.25 * d))


def k_prot_measured(name: str, protein: str = "soy") -> float | None:
    """Dialysis-measured protein-water partition K_prow [L/kg] (zhou2025 Table 1).

    ``protein`` in {"soy", "bsa", "faf_bsa", "phycocyanin"}.  "soy" (soy protein
    isolate) is the grain storage-protein analog.  Returns None if the congener
    was not measured for that protein.
    """
    log = KPROW_ZHOU2025_LOG.get(protein, {}).get(name)
    return None if log is None else float(np.power(10.0, log))


def _kprow_interp(name: str, protein: str) -> float | None:
    """Interpolate log K_prow within the same head-group family (zhou2025) for a
    congener not directly measured.  ``np.interp`` clamps beyond the measured range."""
    if name not in SPECIES:
        return None
    _, npfc, head = SPECIES[name]
    pts = sorted((n, KPROW_ZHOU2025_LOG[protein][sp])
                 for sp, (_, n, h) in SPECIES.items()
                 if h == head and sp in KPROW_ZHOU2025_LOG.get(protein, {}))
    if not pts:
        return None
    return float(np.power(10.0, np.interp(npfc, [p[0] for p in pts], [p[1] for p in pts])))


def k_prot_albumin(name: str) -> float | None:
    """Serum-albumin partition [L/kg] from the HSA K_D (chen2025), single-site:
    ``K_prot = 1/(K_D[mol/L]*MW_HSA[kg/mol])``.  NOTE: this binding-constant route
    OVERESTIMATES the partition by ~50x vs the dialysis K_prow (zhou2025, BSA PFOA:
    110 vs ~5850) -- kept for reference only; prefer :func:`k_prot`."""
    kd = KD_HSA_CHEN2025_UMOL.get(name)
    return None if kd is None else 1.0 / ((kd * 1e-6) * MW_HSA_KG_MOL)


def k_prot(n_perfluoroC: float | None = None, plant: bool = True,
           anchor: float = KPROT_ANCHOR_LKG, *, name: str | None = None,
           protein: str | None = None) -> float:
    """Protein-water partition K_prot [L/kg protein].  [C4, zhou2025]

    Uses the dialysis-MEASURED K_prow: ``protein="soy"`` (soy protein isolate, the
    grain storage-protein analog) for plant tissues, ``"bsa"`` for an animal
    reference; ``plant`` selects the default.  Congeners not measured are
    interpolated within their head-group family; absent that, a U-shaped fallback
    on a nominal anchor is used.
    """
    prot = protein or ("soy" if plant else "bsa")
    base = k_prot_measured(name, prot) if name is not None else None
    if base is None and name is not None:
        base = _kprow_interp(name, prot)
    if base is None:
        if n_perfluoroC is None:
            raise ValueError("k_prot needs a measured `name` or an `n_perfluoroC` (fallback)")
        base = float(anchor) * _kprot_chain_factor(n_perfluoroC)
    return base


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


# --- MEASURED field calibration data: Kim et al. 2019 (Korean paddy) ------------
# Table 4 averages (porewater n=27; soil & brown rice n=30), per congener:
#   (porewater [ng/L], soil [ng/g dw], brown rice [ng/g], rice detection freq [%]).
# Only congeners with detectable brown rice are kept.  These are field ENSEMBLE
# averages (NOT same-site paired), so the BAFs below are approximate anchors.
# (see docs/literature_db/raw_si/kim2019_field_conc.csv, kim2019_grain_baf.csv)
KIM2019_FIELD: dict[str, tuple[float, float, float, float]] = {
    "PFHpA":  (24.7,   0.499,  0.00974, 13.0),
    "PFOA":   (78.7,   0.160,  0.349,   57.0),
    "PFNA":   (13.0,   0.272,  0.0287,  20.0),
    "PFDA":   (3.54,   0.217,  0.00551, 6.7),
    "PFUnDA": (0.244,  0.149,  0.00807, 13.0),
    "PFDoDA": (0.0488, 0.0802, 0.00172, 3.3),
}


def kim2019_grain_baf(basis: str = "porewater") -> dict[str, float]:
    """Per-congener brown-rice (grain) BAF from Kim et al. 2019.  [C1, kim2019]

    ``basis="porewater"`` -> C_grain/C_porewater [L/kg] (the model's BAF, since
    C_w^o is pore water; ng/g == ug/kg, ng/L -> ug/L by /1000).  ``basis="soil"``
    -> C_grain/C_soil [kg/kg].  Field ensemble averages -> treat as approximate.
    """
    out = {}
    for sp, (pw, soil, rice, _df) in KIM2019_FIELD.items():
        out[sp] = rice / (pw / 1000.0) if basis == "porewater" else rice / soil
    return out


# ===========================================================================
# C5 -- root membrane potential for the GHK exclusion term
# ===========================================================================
EM_RICE_ROOT_V = -0.120                 # central value (wang1994: -0.116..-0.140 V)
EM_RICE_ROOT_RANGE_V = (-0.140, -0.116) # rice-specific bounds (wang1994)
EM_PLANT_GENERIC_RANGE_V = (-0.160, -0.090)  # generic plant sweep (sensitivity)
# NOTE: NH4+ (dominant N in flooded paddy) DEPOLARIZES the membrane -> less
# negative E_m -> weaker anion exclusion.  In-situ paddy E_m is a database GAP.


# ===========================================================================
# GAP-B headgroup-specific f_xy offset (task 8: PFSA / ether transport term)
# ===========================================================================
# Root->shoot loading differs by head group at matched chain length.  Expressed
# as a multiplicative offset on the carboxylate f_xy:  f_xy = f_xy_PFCA * exp(off).
# SIGN NOW CONFIRMED (was "uncertain"): PFSA translocates LESS than the
# CF2-comparable PFCA -- the transfer factor TF=tissue/root is consistently lower:
#   * Tang 2026 (paddy rice, low dose): PFOS/PFOA TF_stalk = 0.57/2.22 = 0.26
#   * Yamazaki 2023 (lysimeter): PFOS/PFOA TF(straw/root) = 0.73/1.69 = 0.43
#   geomean ratio ~0.33 -> ln offset ~ -1.1 (refines the placeholder -1.5).
# Ether (GenX/PFECA): Tang GenX/PFOA TF ~0.4-1.2 (variable, single study) -> ~-0.7;
#   provisional, and not among the 12 core congeners.
# (see docs/literature_db/raw_si/tang2026_tf_bcf.csv)
FXY_HEADGROUP_LN_OFFSET: dict[str, float] = {
    "carboxylate": 0.0,
    "sulfonate": -1.1,      # Tang 2026 + Yamazaki 2023 (TF PFOS/PFOA ~0.26-0.43)
    "ether": -0.7,          # Tang 2026 GenX (provisional; ether/replacement)
}


def f_xy_headgroup(f_xy_carboxylate: float, head_group: str = "carboxylate") -> float:
    """Apply the head-group offset to a carboxylate f_xy:  f_xy * exp(offset).  [task 8]

    PFSA (sulfonate) translocates ~3x less than the CF2-matched PFCA (Tang 2026,
    Yamazaki 2023); ether (GenX) is intermediate (provisional).
    """
    return float(f_xy_carboxylate) * float(np.exp(FXY_HEADGROUP_LN_OFFSET[head_group]))


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
        K_prot=k_prot(n_perfluoroC, plant=plant_protein, anchor=kprot_anchor, name=name),
        K_PL=k_pl(n_perfluoroC, head_group, anchor=kpl_anchor, name=name),
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

    print("\n== C4  MEASURED K_PL (chen2025 K_MW) & K_prot (zhou2025 dialysis K_prow) [L/kg] ==")
    for nm in ("PFBA", "PFHxA", "PFOA", "PFDA", "PFDoDA", "PFOS"):
        star = "" if nm in KPROW_ZHOU2025_LOG["soy"] else " (interp)"
        print(f"  {nm:7s}  K_PL={k_pl(name=nm):9.1f}   "
              f"K_prot(soy/plant)={k_prot(name=nm, plant=True):7.1f}   "
              f"K_prot(BSA)={k_prot(name=nm, plant=False):7.1f}{star}")

    print("\n== C1  MEASURED grain BAF (Kim 2019, porewater basis) [L/kg]  [kim2019] ==")
    for nm, b in kim2019_grain_baf("porewater").items():
        print(f"  {nm:7s}  grain BAF = {b:7.2f}")

    print("\n== task 8  head-group f_xy offset (Tang2026 + Yamazaki2023 TF) ==")
    for hg in ("carboxylate", "sulfonate", "ether"):
        print(f"  {hg:11s}  f_xy x exp({FXY_HEADGROUP_LN_OFFSET[hg]:+.1f}) = "
              f"x{np.exp(FXY_HEADGROUP_LN_OFFSET[hg]):.2f}  (PFCA f_xy=0.04 -> {f_xy_headgroup(0.04, hg):.4f})")

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
    Cwo_end = inputs.Cwo_(t[-1])
    Mf = inputs.M_(t[-1])
    straw = (Cend[STEM] * Mf[STEM] + Cend[LEAF] * Mf[LEAF]) / (Mf[STEM] + Mf[LEAF])
    print(f"  final conc [ug/kg]: root={Cend[ROOT]:.3f}  straw={straw:.3f}  grain={Cend[FRUIT]:.3f}")
    print(f"  model BAF  [L/kg] : root={Cend[ROOT]/Cwo_end:.2f}  straw={straw/Cwo_end:.2f}  "
          f"grain={Cend[FRUIT]/Cwo_end:.3f}")
    print(f"  Kim 2019 PFOA grain BAF (porewater) = {kim2019_grain_baf()['PFOA']:.2f} L/kg "
          f"<- Tier-1 calibration target")
    ok = Cend[ROOT] / Cwo_end > straw / Cwo_end > Cend[FRUIT] / Cwo_end
    print(f"  ordering root > straw > grain: {'OK' if ok else 'VIOLATED (transport params need calibration)'}")

    # Tier-1 calibration: fit phloem loading L_Ph to the Kim 2019 PFOA grain BAF.
    try:
        from calibration import calibrate, Param, ObservedBAF
        tgt = kim2019_grain_baf()["PFOA"]
        res = calibrate(model, [Param("L_Ph", 1e-4, 1.0)],
                        [ObservedBAF("grain", tgt, sigma=0.3)], global_search=False)
        print(f"\n== Tier-1 calibration to Kim 2019 PFOA grain BAF ({tgt:.2f} L/kg) ==")
        print(f"  fitted L_Ph = {res.values['L_Ph']:.4f}  ->  grain BAF = "
              f"{res.predicted['grain']:.3f}  (was {Cend[FRUIT]/Cwo_end:.3f})")
        print("  NOTE: grain-only constraint -> fits the phloem knob; f_xy (root->shoot)"
              " needs root/straw data (Kim is grain-only) -- a DB gap.")
    except Exception as e:
        print(f"  (calibration step skipped: {e})")


if __name__ == "__main__":
    _demo()
