"""Request schemas for the web API.

A ``Scenario`` is the whole sidebar of the Streamlit app as one JSON document:
WHICH compound (curated congener / SMILES / neutral organic), HOW the exposure is
supplied (the six data-source modes), and the PFAS mechanism switches. The
service (``server/service.py``) turns it into exactly the ``model_api`` calls the
Streamlit app makes, so the two front ends run the same model on the same inputs.

Defaults reproduce the Streamlit defaults (PFOA, Cwᵒ = 1 µg/L, 120-day season,
measured forcings, ORYZA2000 biomass, the shipped carrier uptake).
"""
from __future__ import annotations

from typing import Dict, List, Literal, Optional

from pydantic import BaseModel, ConfigDict, Field

Lang = Literal["ko", "en"]
Mode = Literal["parametric", "tables", "csv", "hydrus", "inventory", "biomon"]


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid")


class NeutralSpec(_Strict):
    """A NEUTRAL organic (or weak electrolyte): the ``simulate_neutral`` inputs.

    Mirrors ``ui/sidebar._neutral_panel``: log Kow is the one required input; the
    air term, the weak-electrolyte speciation and the apoplastic bypass are opt-in
    and stay structurally absent unless switched on.
    """
    log_kow: float = Field(2.45, ge=-2.0, le=8.0)
    name: str = "carbamazepine"
    half_life: float = Field(7.0, ge=0.0, le=365.0)          # 0 = no metabolism (γ = 0)
    tscf_model: Literal["briggs", "schriever"] = "briggs"
    lipid_source: Literal["measured", "briggs_anchor"] = "measured"
    phloem: bool = False
    air: bool = False
    MW: float = Field(236.3, ge=1.0, le=1e5)
    K_AW: float = Field(1e-3, ge=0.0, le=10.0)
    C_air: float = Field(0.0, ge=0.0)
    weak_electrolyte: bool = False
    pKa: float = Field(4.5, ge=-5.0, le=14.0)
    pH: float = Field(6.5, ge=3.0, le=10.0)
    is_acid: bool = True
    g_apo: float = Field(0.0, ge=0.0, le=50.0)
    Koc: Optional[float] = Field(None, ge=0.0, le=1e9)      # None -> Karickhoff 1981


class CompoundSpec(_Strict):
    kind: Literal["congener", "smiles", "neutral"] = "congener"
    congener: str = "PFOA"
    smiles: Optional[str] = None
    neutral: NeutralSpec = Field(default_factory=NeutralSpec)


class MechanismSpec(_Strict):
    """PFAS levers. ``biomass`` applies to every compound class; the rest are PFAS-only
    and are ignored on the neutral path (exactly as the Streamlit sidebar hides them)."""
    E_m: float = Field(-120.0, ge=-200.0, le=0.0)
    fxy_source: Literal["recommended", "W2fit"] = "recommended"
    biomass: Literal["oryza", "growth_rice"] = "oryza"
    uptake: Literal["carrier", "bypass"] = "carrier"
    lipid_loading: bool = False
    vmax_scale: Optional[float] = Field(None, ge=0.0, le=100.0)  # None -> the mode's own value
    g_apo: Optional[float] = Field(None, ge=0.0, le=100.0)
    km_scale: float = Field(1.0, ge=0.01, le=1000.0)


class TablesSpec(_Strict):
    """Editable growth + pore-water tables, column-major ({column: [values]})."""
    growth: Optional[Dict[str, List[Optional[float]]]] = None   # day, root, stem, leaf, grain
    growth_units: Literal["g/hill", "kg/hill", "g/m2", "kg/ha", "t/ha"] = "g/hill"
    cwo: Optional[Dict[str, List[Optional[float]]]] = None      # day, Cwo
    density: Dict[str, float] = Field(default_factory=lambda: {
        "root": 1.0, "stem": 0.30, "leaf": 0.30, "grain": 1.20})


class BiomonSpec(_Strict):
    conc: Dict[str, float] = Field(default_factory=lambda: {
        "root": 0.49, "straw": 0.83, "grain": 0.46})
    Cwo: Optional[float] = Field(1.0, gt=0.0)


class ExposureSpec(_Strict):
    mode: Mode = "parametric"
    # parametric (and the tables-mode fallback constant)
    Cwo: float = Field(1.0, ge=0.0, le=1e6)
    season: float = Field(120.0, ge=30.0, le=200.0)
    measured_forcing: bool = True
    cwo_profile: Literal["constant", "flooded"] = "constant"
    k_leach: Optional[float] = Field(None, ge=0.0, le=1.0)   # None -> HYDRUS-calibrated default
    # custom tables
    tables: TablesSpec = Field(default_factory=TablesSpec)
    # HYDRUS / CSV drivers (None -> the bundled example)
    csv_text: Optional[str] = Field(None, max_length=5_000_000)
    # live HYDRUS-1D
    f_oc: float = Field(0.02, ge=0.001, le=0.2)
    flood_until: float = Field(90.0, ge=0.0, le=200.0)
    percolation: float = Field(0.30, ge=0.0, le=5.0)
    # soil inventory -> pore water (Freundlich)
    C_total: float = Field(5.0, ge=0.0, le=1e6)
    K_F: float = Field(2.0, gt=0.0, le=1e4)
    n_F: float = Field(0.85, gt=0.0, le=2.0)
    theta_g: float = Field(0.35, gt=0.0, le=1.0)
    inv_flood: bool = True
    inv_k_leach: Optional[float] = Field(None, ge=0.0, le=1.0)
    # biomonitoring: manual values, or a CSV (None + no manual -> bundled example)
    biomon: Optional[BiomonSpec] = None
    biomon_csv_text: Optional[str] = Field(None, max_length=1_000_000)


class Scenario(_Strict):
    compound: CompoundSpec = Field(default_factory=CompoundSpec)
    mechanism: MechanismSpec = Field(default_factory=MechanismSpec)
    exposure: ExposureSpec = Field(default_factory=ExposureSpec)
    lang: Lang = "en"                                      # language of the server messages


class FigureRequest(_Strict):
    scenario: Scenario = Field(default_factory=Scenario)
    name: str
    lang: Lang = "en"
    opts: Dict[str, object] = Field(default_factory=dict)


class InverseRequest(_Strict):
    scenario: Scenario = Field(default_factory=Scenario)
    measured: Dict[str, float]
    sigma_log10: float = Field(0.15, gt=0.0, le=2.0)
    lang: Lang = "en"


class DownloadRequest(_Strict):
    scenario: Scenario = Field(default_factory=Scenario)
    lang: Lang = "en"
    scenario_label: Optional[str] = None


class CwoPreviewRequest(_Strict):
    compound: CompoundSpec = Field(default_factory=CompoundSpec)
    level: float = Field(1.0, ge=0.0)
    season: float = Field(120.0, ge=30.0, le=200.0)
    k_leach: Optional[float] = Field(None, ge=0.0, le=1.0)
