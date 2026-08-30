"""
Neutral-organic (non-ionised) DPU path -- the Briggs/Kow base of the framework
==============================================================================

Everything else in this repo models PFAS: a **permanently dissociated anion**
(`f_d ~ 1`) for which the neutral-compound Briggs/Kow partition core explicitly
does **not** apply (CLAUDE.md section 2).  The neutral base itself
(`docs/dpu_model_summary_corrected.tex`) was derived but never implemented: the
compartment dataclass carries an `fn` field that is fixed at 0.0 and appears in
no equation, the neutral passive-uptake term exists only as a comment, and there
is no neutral observation anywhere in `data_obs/`.  This module supplies the
missing path.

Why it is worth having
----------------------
Every PFAS validation in this repo is entangled with PFAS-specific parameters
that had to be FIT (`f_xy`, `k_seq`, the lipid conductances) -- which is exactly
why the honest verdicts read "reproduction, not prediction".  A neutral compound
is different: its partitioning and its root->shoot loading are both fixed from
OUTSIDE the model by published QSPRs on log Kow, with nothing left to tune.  So
the neutral path is the one setting in which the **DPU backbone itself** --
four compartments, xylem advection, growth dilution, terminal accumulation --
can be tested independently of the ionic extension.  A failure here would be a
failure of the skeleton; a failure on the PFAS side alone is not.

The key structural insight: **no new ODE is needed.**  A neutral compound is the
existing 4-compartment model with the ionic machinery switched off by physics
rather than by special-casing:

  * **valence z = 0** => the electrochemical driving force `N = zEF/RT` is 0, so
    the GHK factor `N/(e^N - 1)` -> 1 and the exclusion factor `e^N` -> 1.  The
    membrane term `kappa_d*g*(Cwo - e^N*Cw_root)` degenerates *exactly* to
    passive Fickian diffusion `kappa_d*(Cwo - Cw_root)` -- the "optional neutral
    passive term" the module headers describe.  Anion exclusion (`e^N ~ 107`)
    simply is not there for an uncharged molecule.
  * **Vmax = 0** => no carrier.  PFAS need one to overcome exclusion; neutral
    organics cross the membrane unaided.
  * **binding by Briggs instead of basis-A**: with `f_prot = f_cw = 0` and
    `K_PL = a*Kow^b`, the existing `binding_factors` returns
    `theta + (1-theta)*f_PL*a*Kow^b = W + L*a*Kow^b` -- the Trapp/Briggs
    plant-water partition coefficient K_PW, term for term.
  * **f_xy = TSCF(log Kow)**, the Briggs bell, *computed* rather than fitted.

Is TSCF double-counted?  No, and it is worth being explicit.  Briggs defines TSCF
as C_xylem / C_external.  Here the membrane term is diffusive, so it drives the
root free concentration toward the external one (`Cw_root -> Cwo`), and the xylem
is then loaded at `Cw_xyl = TSCF * Cw_root ~ TSCF * Cwo` -- Briggs' definition,
recovered rather than applied twice.  (`kappa_d` therefore sets how fast the root
equilibrates, not the level it equilibrates to.)

Anchors -- VERIFIED AT SOURCE (Briggs, Bromilow & Evans 1982, Pestic. Sci.
13:495-504, eqs. 2 and 3 of that paper; they also appear in this repo's
`docs/theory_anchor.tex` and are reproduced verbatim as eq. 4 of Schriever &
Lamshoeft 2020, Sci. Total Environ. 713:136667):

    TSCF = 0.784 * exp[ -(log Kow - 1.78)^2 / 2.44 ]        (bell, peak 0.784)
    RCF  = 0.82 + 10^(0.77*log Kow - 1.52)                  (barley roots)

The RCF relation is the same K_PW form: `W = 0.82`, `L*a = 10^-1.52 = 0.0302`,
`b = 0.77`.  With the conventional lipid-vs-octanol correction `a = 1.22` that is
a barley-root lipid content `L = 0.0247`; `a` and `L` are only identifiable as
their product, so `a` is a convention and `L*a` is the anchored quantity.
`briggs_root_compartment()` returns that anchor exactly, which makes the
"does the model reproduce Briggs" check a real test rather than a tautology.

Briggs derives the 0.82 floor by subtracting the fitted lipid term from the
measured RCF of nine polar compounds, and attributes it to "equilibration of the
chemical between the external solution and the water contained within the roots"
-- noting that roots at 90% water would give ~0.9.  That is an independent
confirmation that the first term of K_PW is the tissue WATER content, which is
exactly how this module maps it onto `Compartment.theta`.

Schriever & Lamshoeft 2020 refit the same Gaussian form to a larger, modern
compilation (97 TSCF values from intact-plant hydroponic tests, 42 compounds),
in log D rather than log Kow: A = 0.746, B = 2.160, C = 7.230 (two outliers
removed; A = 0.713, B = 2.209, C = 7.862 with them).  The peak HEIGHT is
essentially Briggs' (0.75 vs 0.78) but the bell is much BROADER (C 7.2 vs 2.4)
and shifted right.  `schriever_tscf()` exposes it, and `simulate_neutral(
tscf_model="schriever")` runs with it -- the difference between the two is a
fair estimate of how well TSCF is actually known, which matters because on this
path TSCF is not a fitted parameter but an input.

Scope and honesty
-----------------
* Air exchange (volatilisation / gaseous uptake, sections 6.3-6.5 of the tex) is
  implemented in `plant_air` and is **opt-in**: `simulate_neutral(air=True)`.  It
  is off by default because it needs two inputs the Kow-only a-priori run does not
  have -- `K_AW` and a molar mass -- so switching it on silently would quietly
  change what "a-priori from log Kow alone" means.  With it off this module is
  valid for **non-volatile** neutral organics (low `K_AW`: most pesticides,
  carbamazepine, neonicotinoids, triazoles) and `k_aw_warning()` flags a compound
  where that assumption is likely violated; with it on, that warning is replaced
  by an actual run.  The flux scales with the tissue's specific SURFACE AREA,
  which this repo has only ever used as a leaf/grain ratio -- see the `plant_air`
  header before trusting an absolute volatilisation magnitude.
* Tissue lipid contents are taken from Trapp, McFarlane & Matthies 1994 (the
  canonical validation of this framework): root 1%, stem and leaf 3% of FRESH
  weight -- the same basis as the water contents, as `K_PW = W + L*a*Kow^b`
  requires (see TRAPP1994_LIPID_FW).  Those are SOYBEAN values, not rice -- an
  organ-resolved total-lipid measurement for rice remains a genuine gap -- but
  they are a cited anchor and they are what the DPU base itself was exercised
  with.  Note `params/parameters.json` cannot supply them: it carries
  *phospholipid* fractions (membrane binding of anions), a lower bound on the
  total lipid Briggs partitioning refers to.
  The ROOT value is the one live disagreement in this module: 1% fw sits 2.48x
  below the lipid term Briggs' own regression implies, and which of the two is
  right is NOT settled -- so it is a named mode, `lipid_source="measured"`
  (default) vs `"briggs_anchor"`, not a constant.  See LIPID_SOURCES, and run
  `validation/neutral_dpu_validation.py --lipid-source both` for the evidence.
  Every published number from this module is on the default.
* Metabolism `gamma` is genuinely non-zero for most neutral organics (unlike
  PFAS), so it is exposed per compartment and defaults to 0 only so that a run
  without a measured half-life is obviously an upper bound.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, replace

import numpy as np

from pfas_rice_plant_module_4pool_surf import (
    Compartment, Compound, Environment, PlantInputs, RiceUptakeModel,
    binding_factors, ROOT, STEM, LEAF, FRUIT,
    _speciation, P_N_OVER_P_D, LEAF_CYTOSOL_PH,
)

# --- Briggs 1982 coefficients (docs/theory_anchor.tex eqs. briggsT / briggsR) --
TSCF_MAX = 0.784          # bell maximum [-]
TSCF_LOGKOW_PEAK = 1.78   # log Kow at the peak [-]
TSCF_WIDTH = 2.44         # Gaussian width parameter [-]
RCF_FLOOR = 0.82          # water / free-space floor of the root [L/kg]
RCF_SLOPE = 0.77          # b: lipophilic sorption slope [-]
RCF_INTERCEPT = -1.52     # log10 of (L*a) for macerated barley roots [-]
# lipid-vs-octanol correction; only the PRODUCT L*a is identifiable (see header)
LIPID_OCTANOL_A = 1.22


def briggs_tscf(log_kow: float) -> float:
    """Transpiration stream concentration factor (Briggs 1982 eq. 3).

    The bell: uptake into the transpiration stream is most efficient at
    intermediate lipophilicity (peak 0.784 at log Kow 1.78). Very polar compounds
    do not cross the membrane; very lipophilic ones are retained by the root.
    Fitted to 17 NON-IONISED compounds -- do not apply it to a dissociating species
    (that is precisely why the PFAS side needs its own `f_xy`).
    """
    return TSCF_MAX * float(np.exp(-((log_kow - TSCF_LOGKOW_PEAK) ** 2) / TSCF_WIDTH))


def briggs_rcf(log_kow: float) -> float:
    """Root concentration factor for barley roots (Briggs 1982 eq. 2):
    log(RCF - 0.82) = 0.77 log Kow - 1.52."""
    return RCF_FLOOR + 10.0 ** (RCF_SLOPE * log_kow + RCF_INTERCEPT)


# Schriever & Lamshoeft 2020 (STOTEN 713:136667) Table 3, Gaussian refit of the
# same functional form to 97 intact-plant TSCF values, in log D. Their "n=78"
# row (two outliers removed) is the reported best fit.
SCHRIEVER_TSCF = dict(A=0.746, B=2.160, C=7.230)
SCHRIEVER_TSCF_ALL = dict(A=0.713, B=2.209, C=7.862)


def schriever_tscf(log_d: float, params: dict | None = None) -> float:
    """TSCF from the Schriever & Lamshoeft 2020 Gaussian refit.

    Same form as Briggs, fitted to a larger modern compilation and expressed in
    log D (which equals log Kow for a non-ionisable compound -- the only kind this
    module handles, so the two descriptors coincide here). Nearly the same peak
    height as Briggs but a much broader bell: use it to see how much of a result
    rests on TSCF being known precisely.
    """
    p = params or SCHRIEVER_TSCF
    return p["A"] * float(np.exp(-((log_d - p["B"]) ** 2) / p["C"]))


def tscf(log_kow: float, model: str = "briggs") -> float:
    """TSCF from the selected published QSPR ('briggs' default, or 'schriever')."""
    if model == "briggs":
        return briggs_tscf(log_kow)
    if model == "schriever":
        return schriever_tscf(log_kow)
    raise ValueError(f"unknown TSCF model {model!r}; use 'briggs' or 'schriever'")


# --- Briggs 1983 STEM coefficients (Pestic. Sci. 14:492-500, their eqs. 2 and 3;
# VERIFIED AT SOURCE) --------------------------------------------------------
# The companion paper to Briggs 1982, and the only published anchor for the STEM
# compartment. It has exactly the K_PW form -- a water floor plus a lipophilic
# term -- but fitted to barley shoot bases rather than roots:
#
#     log(K_stem/xylem_sap - 0.82) = 0.95*log Kow - 2.05                  (eq. 2)
#     SCF = K_stem/xylem_sap * TSCF                                       (eq. 3)
#
# so the anchored stem lipid term is L*a = 10^-2.05 = 0.0089 with b = 0.95.
# Note what the repo runs instead: L = 0.03 (Trapp 1994 soybean) with the
# conventional a = 1.22 and Briggs' ROOT exponent b = 0.77, i.e. L*a = 0.0366 --
# 4.1x ABOVE this anchor, and with a flatter slope. That is the same class of
# mismatch as the root's (BRIGGS_ANCHORED_LIPID_FW) but in the OPPOSITE
# direction, and it is likewise left as-is rather than silently corrected.
# `validation/briggs1983_stem.py` measures what it costs.
STEM_SLOPE = 0.95         # b for the stem/xylem-sap partition [-]
STEM_INTERCEPT = -2.05    # log10 of (L*a) for barley shoot bases [-]
STEM_FLOOR = 0.82         # water floor, the same value Briggs found for roots


def briggs_stem_xylem_partition(log_kow: float) -> float:
    """K_stem/xylem sap (Briggs 1983 eq. 2) [L/kg]."""
    return STEM_FLOOR + 10.0 ** (STEM_SLOPE * log_kow + STEM_INTERCEPT)


def briggs_scf(log_kow: float, tscf_model: str = "briggs") -> float:
    """Stem concentration factor (Briggs 1983 eq. 3): the stem/xylem-sap
    partition times the fraction of the external solution reaching the xylem.

    Reproduces the paper's own stated result -- a maximum of "about 6 ... at
    about log Kow = 4.5" -- which is what makes this a transcription check
    rather than an assertion.
    """
    return briggs_stem_xylem_partition(log_kow) * tscf(log_kow, model=tscf_model)


def k_pw(log_kow: float, W: float, L: float, a: float = LIPID_OCTANOL_A,
         b: float = RCF_SLOPE) -> float:
    """Plant-water partition coefficient K_PW = W + L*a*Kow^b  [L/kg fw].

    W = fresh-weight water content [L/kg], L = fresh-weight lipid fraction [kg/kg].
    This is the neutral analogue of the PFAS `B_k`; both map a free aqueous
    concentration to a tissue concentration, so the ODE is indifferent to which
    produced the number.
    """
    return W + L * a * 10.0 ** (b * log_kow)


# ---------------------------------------------------------------------------
# compound / compartment construction
# ---------------------------------------------------------------------------
@dataclass
class NeutralCompound:
    """A non-ionised organic: everything the model needs comes from log Kow.

    kappa_d : root-membrane permeability-area product [L/(day kg)]. For a neutral
        molecule this is plain passive diffusion (no GHK exclusion). It sets how
        fast the root equilibrates with the pore water, NOT the equilibrium level
        (that is K_PW), so results are insensitive to it once it is fast relative
        to the season -- unlike the PFAS case where it fights anion exclusion.
    K_AW : air-water partition coefficient [-]; only used to warn, since air
        exchange is not implemented (see the module header).
    gamma : first-order metabolism [1/day] -- per compartment, 0 = recalcitrant.
    """
    name: str
    log_kow: float
    MW: float = float("nan")          # g/mol (reporting only; no air terms)
    K_AW: float = 0.0
    kappa_d: float = 20.0
    gamma: float = 0.0
    tscf: float | None = None         # override the QSPR if TSCF was measured
    tscf_model: str = "briggs"        # 'briggs' (default) or 'schriever'
    # --- weak electrolyte (opt-in; None = the strictly-neutral path above) -----
    # With a pKa the compound stops being purely neutral and becomes a weak
    # acid/base: part neutral molecule, part ion. That case cannot be expressed by
    # this module's z=0 trick -- see the speciation block in the plant module --
    # so it is routed through the (fn, fd) pair instead. `pKa=None` keeps every
    # published neutral number bit-identical.
    pKa: float | None = None
    is_acid: bool = True              # False -> weak BASE (its ion is a CATION)
    # APOPLASTIC BYPASS conductance [L/(day kg)] -- entry that never crosses a
    # membrane, so it feels neither the (fn, fd) weighting nor the GHK factor.
    # 0 (default) removes the term identically. It exists for the weak-electrolyte
    # case: section 4l of docs/neutral_dpu_validation.md showed the membrane-only
    # model orders ionisable compounds well and under-delivers by orders of
    # magnitude, and this is the one lever that can raise the floor without
    # flattening that ordering. EXPLORATORY -- see `fit_g_apo` in
    # validation/weak_electrolyte_tscf.py for what the data say the value is.
    g_apo: float = 0.0
    pH: float = 6.5                   # root-zone pH setting (fn, fd); paddy pore water

    @property
    def TSCF(self) -> float:
        if self.tscf is not None:
            return float(self.tscf)
        return tscf(self.log_kow, self.tscf_model)


def neutral_environment() -> Environment:
    """Environment with valence z = 0: no electrochemical driving force, so the
    GHK term degenerates to passive diffusion and `e^N = 1` (no anion exclusion).
    This is the single switch that turns the ionic model into the neutral one."""
    return Environment(z=0)


def neutral_compound(c: NeutralCompound, a: float = LIPID_OCTANOL_A,
                     b: float = RCF_SLOPE) -> Compound:
    """Adapt a NeutralCompound onto the core `Compound` container.

    K_prot = K_cw = 0 and K_PL = a*Kow^b, so `binding_factors` evaluates the
    Briggs/Trapp K_PW (see header). Vmax = 0 removes the carrier; f_xy is the
    Briggs TSCF, computed rather than fitted.
    """
    kw = dict(
        name=c.name, K_prot=0.0, K_PL=a * 10.0 ** (b * c.log_kow), K_cw=0.0,
        kappa_d=float(c.kappa_d),
        Vmax_in=0.0, Km_in=1.0, Vmax_out=0.0, Km_out=1.0,
        L_Ph=0.0,                      # no phloem in the neutral base (see simulate_neutral)
        f_xy=c.TSCF,
        # NOTE on fd: in the ionic core `fd` multiplies the membrane term as the
        # DISSOCIATED fraction. Here it is 1.0 not because the compound is
        # dissociated -- it is not -- but because with z=0 that term has already
        # become the neutral passive-diffusion flux kappa_d*(Cwo - Cw_root), and
        # fd is simply its (unit) prefactor. `fn=1` records the actual speciation.
        # Setting fd=0 would zero the uptake entirely, which is the wrong reading.
        fd=1.0,
        fn=1.0,
        # apoplastic bypass; 0 by default, so the term is structurally absent
        g_apo=float(c.g_apo),
    )
    kw.update(_weak_electrolyte_kw(c))    # {} unless a pKa was supplied
    return Compound(**kw)


def _weak_electrolyte_kw(c: NeutralCompound) -> dict:
    """The extra `Compound` fields a WEAK ELECTROLYTE needs; `{}` when `pKa` is
    None, so the strictly-neutral path is untouched.

    Two choices are made here and both are calibrated to stay CONTINUOUS with the
    neutral limit above, which is the property that makes this safe to add:

    `P_n = kappa_d`. Above, `kappa_d` is the neutral molecule's passive
    conductance (with z=0 the GHK term has already degenerated to Fickian), so the
    neutral species keeps exactly the conductance it had. A weak acid with pKa far
    ABOVE the pH is then numerically the same compound as the same molecule run
    with `pKa=None` -- no discontinuity at the boundary between the two modes.

    `kappa_d -> kappa_d * 10^-3.5` for the ION. That is Trapp's relation
    (`P_N_OVER_P_D`): the charged species is ~10^3.5 times less membrane-permeable
    than its neutral form. So as pKa falls the compound does not merely lose the
    neutral route, it becomes genuinely harder to take up -- which is the physics,
    and is why a permanent anion needs a carrier at all.

    NOTE this is NOT how a PFAS is run: `simulate()` fits `kappa_d` directly as the
    lumped ionic conductance and does not derive it from a neutral counterpart.
    """
    if c.pKa is None:
        return {}
    fn, fd = _speciation(c.pKa, c.pH, c.is_acid)
    return dict(
        pKa=float(c.pKa), is_acid=bool(c.is_acid),
        z=(-1 if c.is_acid else +1),
        P_n=float(c.kappa_d),
        kappa_d=float(c.kappa_d) / P_N_OVER_P_D,
        fn=float(fn), fd=float(fd),
    )


def neutral_compartment(name: str, W: float, L: float, S: float = 0.0,
                        gamma: float = 0.0) -> Compartment:
    """Compartment carrying Briggs composition: water content W [L/kg fw] and
    fresh-weight lipid fraction L [kg/kg fw].

    `binding_factors` rescales dry-weight fractions by (1 - theta), so the stored
    `f_PL` is L/(1-W); the product returned is then exactly W + L*a*Kow^b.
    """
    if not 0.0 <= W < 1.0:
        raise ValueError(f"water content W must be in [0,1): got {W}")
    return Compartment(name=name, theta=W, f_prot=0.0, f_PL=L / (1.0 - W),
                       f_cw=0.0, S=S, gamma=gamma)


def briggs_root_compartment() -> Compartment:
    """The barley root Briggs actually measured: W = 0.82, L*a = 10^-1.52.

    Reproduces `briggs_rcf` exactly through `binding_factors`, which is what makes
    `validation/neutral_dpu_validation.py`'s partition check a real test of the
    adapter rather than a restatement of the formula.
    """
    L = 10.0 ** RCF_INTERCEPT / LIPID_OCTANOL_A
    return neutral_compartment("root", W=RCF_FLOOR, L=L)


# Tissue composition for NEUTRAL partitioning.
#   WATER: the repo's own measured rice `theta_fw` (params/parameters.json).
#   LIPID: from Trapp, McFarlane & Matthies 1994 (Environ. Toxicol. Chem.
#     13:413-422), the canonical validation of this very framework, which states
#     "root lipid content l_R = 1%; stem and leaf lipid content l_St = 3%,
#     l_L = 3%" for its soybean runs. These are the values the DPU base itself was
#     exercised with. They are NOT rice measurements -- rice-specific total lipid
#     resolved by organ is still a genuine gap -- but they are a cited anchor
#     rather than the guesses that stood here before.
# BASIS (important, and easy to get wrong): these are FRESH-weight fractions, the
# same basis as W. K_PW = W + L*a*Kow^b only adds up if both terms share a basis,
# and Briggs' own barley anchor is L*a = 10^-1.52 => L = 0.0247, i.e. 2.5% of
# FRESH weight -- next to his W = 0.82. Trapp likewise lists l_R = 1% alongside
# "water content root W_R = 94.2%". Treating them as dry-weight fractions (as an
# earlier revision of this file did) under-states K_PW ~10x for lipophilic
# compounds; the Liu 2023 root data caught it (log10 RMSE 0.605 -> 0.198).
RICE_WATER = {"root": 0.90, "stem": 0.83, "leaf": 0.78, "grain": 0.14}
TRAPP1994_LIPID_FW = {"root": 0.01, "stem": 0.03, "leaf": 0.03, "grain": 0.03}

# CORROBORATION of the root value, added when Acquisition_Queue.csv row A1
# arrived. Li, Chiou, Li & Schnoor 2019 (Environ. Int. 126:46-53) tabulate root
# lipid for 13 crops on an explicitly FRESH-weight basis (they convert dry-basis
# reports at 90% root water) inside an RCF expression of this same form, so their
# values are commensurable with L here. Their cereals -- barley 1.00%, wheat
# 1.10-1.14%, maize 0.53% -- bracket the 1% above. RICE IS ABSENT from their
# table, so the rice-specific gap is NOT closed; what improves is the provenance
# of the number, from a soybean model run to measured cereal roots.
#
# THE ANCHOR DISCREPANCY -- read before changing anything here.
# Briggs' fitted lipid term is L*a = 10^-1.52 = 0.0302. Substituting a measured
# L = 0.01 while keeping the conventional a = 1.22 changes that product to
# 0.0122, i.e. the rice root runs 2.5x BELOW the anchor this module claims. The
# header says L and a are identifiable only as their product; this dict replaces
# one factor and leaves the other, which is not a neutral substitution.
# Measured on Briggs' own barley rows the cost is log10 RMSE 0.266 against the
# anchor's 0.111 (validation/li2019_rcf_apriori.py section 3).
# It is left as-is DELIBERATELY. Restoring the anchor improves Li 2019's
# hydroponic rows and Ge 2017 but DEGRADES three tables -- Liu 2023 (the only
# rice one), Kodesova 2019, and Li 2019's own 376-row SOIL half, which is 12x
# larger than the hydroponic half arguing for it. Picking an intermediate L
# would make the path's headline claim -- nothing is fitted -- false.
# The physical reading is that Briggs' coefficient carries non-lipid
# sorption (cell wall, lignin: this repo's f_cw*K_cw, GAP A) that the neutral
# composition zeroes, so the fix is a measured cell-wall coefficient, not a
# bigger lipid fraction. Run the anchor with the named mode instead:
#     rice_compartments(lipid_source="briggs_anchor")
BRIGGS_ANCHORED_LIPID_FW = dict(TRAPP1994_LIPID_FW,
                                root=10.0 ** RCF_INTERCEPT / LIPID_OCTANOL_A)

# The two readings above as a NAMED MODE, in the repo's own idiom (`f_xy_source`,
# `cwo_profile`, `biomass`, `tscf_model`). The point is that the question this
# selects between is genuinely unsettled -- see the discrepancy note above and
# docs/neutral_dpu_validation.md section 5 -- so the alternative is kept runnable
# and testable rather than decided by whichever constant got written down.
#
#   "measured"      L_root = 1% fw. Measured cereal root lipid (Li 2019 A1:
#                   barley 1.00, wheat 1.10-1.14, maize 0.53%), adopted from
#                   Trapp 1994. DEFAULT.
#   "briggs_anchor" L_root = 2.47% fw, i.e. L*a = 10^-1.52 exactly. NOT a
#                   measurement -- the lipid term Briggs' 1982 barley regression
#                   implies. Briggs measured no lipid at all.
#
# The evidence is 3 tables to 1 AGAINST restoring the anchor, and the one for it
# is the hydroponic half of a paper whose own soil half (12x larger) says the
# opposite; `validation/neutral_dpu_validation.py --lipid-source both` prints the
# comparison on any table. The default is "measured" because the evidence no
# longer supports moving -- NOT because 1% is validated for rice, which it is
# not. Only the PRODUCT L*a is identifiable, so raising L to 0.0247 and raising
# `a` to 3.02 are the same model; note `a = 1.22` itself has no citation here.
LIPID_SOURCES = {
    "measured": TRAPP1994_LIPID_FW,
    "briggs_anchor": BRIGGS_ANCHORED_LIPID_FW,
}
DEFAULT_LIPID_SOURCE = "measured"
RICE_SURFACE = {"leaf": 20.0, "grain": 2.0}
# Briggs' sorption exponent b is 0.77 for roots and stems; Trapp 1994 uses 0.95
# for LEAVES (their eq. 1 discussion), i.e. leaf lipid is closer to octanol.
LEAF_SORPTION_EXPONENT = 0.95


def rice_compartments(lipids: dict | None = None, waters: dict | None = None,
                      gammas: dict | None = None,
                      lipid_source: str = DEFAULT_LIPID_SOURCE) -> list[Compartment]:
    """[root, stem, leaf, grain] with neutral (Briggs) composition.

    `lipids` are FRESH-weight fractions, the same basis as the water contents
    (see the note on TRAPP1994_LIPID_FW -- mixing the bases is the classic error
    here and understates K_PW by ~10x for lipophilic compounds).

    `lipid_source` selects the tissue-lipid reading -- "measured" (default) or
    "briggs_anchor"; see LIPID_SOURCES for what each one is and why the choice is
    still open. `lipids` overrides individual tissues on top of whichever is
    selected, so an explicit dict still wins.
    """
    if lipid_source not in LIPID_SOURCES:
        raise ValueError(f"unknown lipid_source {lipid_source!r}; "
                         f"expected one of {sorted(LIPID_SOURCES)}")
    W = dict(RICE_WATER, **(waters or {}))
    L = dict(LIPID_SOURCES[lipid_source], **(lipids or {}))
    g = dict.fromkeys(W, 0.0)
    g.update(gammas or {})
    return [neutral_compartment(k, W=W[k], L=L[k],
                                S=RICE_SURFACE.get(k, 0.0), gamma=g[k])
            for k in ("root", "stem", "leaf", "grain")]


K_AW_VOLATILE = 1e-4      # above this, volatilisation rivals the other leaf sinks


def k_aw_warning(c: NeutralCompound) -> str | None:
    """Flag a compound whose result is an upper bound because air exchange is OFF.

    K_AW above ~1e-4 puts volatilisation on a par with the other loss terms for a
    leaf, so a run without air terms over-predicts the shoot. Returns None when
    the non-volatile assumption is defensible.

    This is the warning for a run with `air` disabled. Air exchange is now
    IMPLEMENTED (`plant_air`, `simulate_neutral(air=True)`), so the remedy is to
    run it rather than to caveat it -- which is what the message says.
    """
    if c.K_AW and c.K_AW > K_AW_VOLATILE:
        return (f"{c.name}: K_AW={c.K_AW:.2e} -- volatilisation is not included in "
                "THIS run, so shoot concentrations are an UPPER bound. Re-run with "
                "simulate_neutral(..., air=True) to model it (needs MW); see "
                "src/plant_air.py and docs/dpu_model_summary_corrected.tex sec:permeability.")
    return None


def air_exchange(c: NeutralCompound, **kw):
    """Build a `plant_air.AirExchange` for a NeutralCompound.

    Pulls `K_AW`, `MW` and `log_kow` off the compound; `kw` overrides anything on
    `AirExchange` (ambient `C_air`, `rh`, `T`, the specific areas `S`, ...).
    Raises if the compound carries no molar mass, since every air-side equation
    needs one -- better than silently running with a wrong permeability.
    """
    from plant_air import AirExchange
    return AirExchange(K_AW=float(c.K_AW), MW=float(c.MW),
                       log_kow=float(c.log_kow), **kw)


# ---------------------------------------------------------------------------
# forward run
# ---------------------------------------------------------------------------
def simulate_neutral(cmpd: NeutralCompound, drivers: dict, comps=None,
                     phloem=False, phi: float = 0.1, T_C_Ph: float = 10.0,
                     L_Ph: float = 1.0, tscf_model=None, C0=None,
                     air=False, air_kw=None,
                     lipid_source: str = DEFAULT_LIPID_SOURCE):
    """Run the 4-compartment DPU for a NEUTRAL organic.

    drivers : {t, Cwo, Qtp, M} on a common grid -- the same driver contract as
        `model_api.simulate(drivers=...)`; M is (n_t, 4) organ fresh mass [kg].

    lipid_source : which tissue-lipid reading to build the compartments from,
        "measured" (default) or "briggs_anchor" -- see LIPID_SOURCES. Ignored
        when `comps` is passed, since those are already built.

    air : plant-air exchange (volatilisation + gaseous uptake), OFF by default.
        `True` builds a `plant_air.AirExchange` from the compound's `K_AW`, `MW`
        and `log_kow` (`air_kw` passes through ambient `C_air`, `rh`, `T`, the
        specific areas `S`, ...); an `AirExchange` instance may be given directly.
        It is off by default because it needs `K_AW` and `MW`, which the strict
        Kow-only a-priori run does not use -- turning it on silently would change
        what the published RMSEs mean. With `K_AW = 0` the terms are identically
        zero, so enabling it cannot perturb a non-volatile result (tested).

    phloem : OFF by default, faithfully to the neutral base -- "no dissociation,
        pH-dependent speciation, membrane electrical potential, ion-trap, or
        phloem transport is included" (dpu_model_summary_corrected.tex section 2).
        The phloem in the PFAS core is an addition of the ionisable extension
        (the grain is phloem-fed there). With phloem off the grain is reached only
        by its share of the xylem, split by surface area, as the neutral
        derivation specifies. Turning it on (with `L_Ph` the loading factor) is
        an explicit departure from the base, not the default -- and it matters:
        an unrestricted neutral phloem (L_Ph=1) drives the grain, a terminal
        accumulator of small mass, to BAF ~20-60, which is a statement about
        assuming unrestricted loading, not a prediction of the neutral base.

    Returns a dict with the usual tissue series plus the neutral diagnostics
    (K_PW per tissue, the Briggs TSCF actually used, RCF/TF at maturity).
    """
    t = np.asarray(drivers["t"], dtype=float)
    Cwo = np.asarray(drivers["Cwo"], dtype=float)
    Qtp = np.asarray(drivers["Qtp"], dtype=float)
    M = np.asarray(drivers["M"], dtype=float)
    comps = (comps if comps is not None
             else rice_compartments(lipid_source=lipid_source))
    if cmpd.pKa is not None and comps[LEAF].pH is None:
        # give the leaf a cytosolic pH so the phloem ION TRAP is available. It only
        # does anything when `phloem=True`; with the neutral base's default
        # (phloem off, L_Ph=0) this is inert. Copy rather than mutate a
        # caller-supplied compartment list.
        comps = list(comps)
        comps[LEAF] = replace(comps[LEAF], pH=LEAF_CYTOSOL_PH)

    if tscf_model is not None:
        cmpd = replace(cmpd, tscf_model=tscf_model)
    core = neutral_compound(cmpd)
    if phloem:
        core.L_Ph = float(L_Ph)
    env = neutral_environment()
    inputs = PlantInputs(t=t, Cwo=Cwo, Qtp=Qtp, M=M,
                         leaf_loss=drivers.get("leaf_loss"))
    air_obj = None
    if air is not False and air is not None:
        air_obj = air_exchange(cmpd, **(air_kw or {})) if air is True else air
    model = RiceUptakeModel(env=env, cmpd=core, comps=comps, inputs=inputs,
                            phi=(phi if phloem else 0.0),
                            T_C_Ph=(T_C_Ph if phloem else 0.0),
                            air=air_obj)
    sol = model.solve(t, C0=C0)
    Y = sol.y
    K = binding_factors(comps, core)
    names = ("root", "stem", "leaf", "grain")
    conc = {k: Y[i] for i, k in enumerate(names)}
    ref = float(Cwo[-1]) if Cwo[-1] > 0 else float(np.nanmax(Cwo)) or 1.0
    baf_final = {k: float(v[-1] / ref) for k, v in conc.items()}
    Mf = M[-1]
    straw = (conc["stem"] * Mf[STEM] + conc["leaf"] * Mf[LEAF]) / (Mf[STEM] + Mf[LEAF])
    return dict(
        t=t, compound=cmpd.name, success=bool(sol.success), conc=conc,
        Cwo=Cwo, Qtp=Qtp, M=M, straw=straw,
        baf={k: v / ref for k, v in conc.items()}, baf_final=baf_final,
        straw_baf=float(straw[-1] / ref),
        tf_final={k: baf_final[k] / max(baf_final["root"], 1e-12) for k in names},
        K_PW={k: float(K[i]) for i, k in enumerate(names)},
        phloem=bool(phloem),
        TSCF=float(core.f_xy), log_kow=float(cmpd.log_kow),
        rcf_briggs=briggs_rcf(cmpd.log_kow),
        # air exchange: the diagnostics (per-organ permeabilities and the implied
        # volatilisation half-life) only when it is actually on; the "this is an
        # upper bound" warning only when it is not.
        air=air_obj is not None,
        air_summary=(air_obj.summary(comps, K, M=M[-1], Qtp=float(Qtp[-1]))
                     if air_obj is not None else None),
        warning=(None if air_obj is not None else k_aw_warning(cmpd)),
        N=float(env.N), eN=float(np.exp(env.N)),      # 0 and 1: no exclusion
    )


if __name__ == "__main__":
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import forcing_rice as fr
    from growth_rice import organ_biomass

    season, n = 120.0, 241
    t = np.linspace(0.0, season, n)
    b = organ_biomass(t, season)
    drivers = dict(t=t, Cwo=np.full(n, 1.0), Qtp=fr.Q_TP(t, season),
                   M=np.maximum(np.column_stack(
                       [b["root"], b["stem"], b["leaf"], b["grain"]]), 1e-4))
    print(f"{'compound':16}{'logKow':>8}{'TSCF':>8}{'K_PW root':>11}"
          f"{'root':>9}{'straw':>9}{'grain':>9}")
    for c in (NeutralCompound("thiamethoxam", -0.13),
              NeutralCompound("imidacloprid", 0.57),
              NeutralCompound("carbamazepine", 2.45),
              NeutralCompound("isoprothiolane", 3.30),
              NeutralCompound("difenoconazole", 4.36)):
        r = simulate_neutral(c, drivers)
        print(f"{c.name:16}{c.log_kow:>8.2f}{r['TSCF']:>8.3f}{r['K_PW']['root']:>11.2f}"
              f"{r['baf_final']['root']:>9.2f}{r['straw_baf']:>9.2f}"
              f"{r['baf_final']['grain']:>9.2f}")
