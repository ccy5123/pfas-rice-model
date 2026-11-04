"""Request schemas for the web API.

A ``Scenario`` is the whole sidebar of the Streamlit app as one JSON document:
WHICH compound (curated congener / SMILES / neutral organic), HOW the exposure is
supplied (the six data-source modes), and the PFAS mechanism switches. The
service (``server/service.py``) turns it into exactly the ``model_api`` calls the
Streamlit app makes (``ui/sidebar.build`` + ``ui/common.run_model``), so the two
front ends run the same model on the same inputs.

Defaults are the model's own (``model_api``): PFOA, Cwᵒ = 1 µg/L, the calibration
site-year crop drivers (Yamazaki 2020, Tsukuba; ORYZA2000 port), the redistributed
shoot with f_xy fitted per congener and the free anion (lipid loading off). A field
left ``None`` is filled server-side from the same source the Streamlit widget is
seeded from (e.g. the soil-inventory defaults), never from a number typed here.
"""
from __future__ import annotations

from typing import Dict, List, Literal, Optional

from pydantic import BaseModel, ConfigDict, Field

import model_api as api

Lang = Literal["ko", "en"]
Mode = Literal["parametric", "tables", "csv", "hydrus", "inventory", "biomon"]
# ui.common.CROPS: the calibration site-year (default) or a generic-climatology biomass driver
Crop = Literal["site2020", "oryza2000", "oryza", "growth_rice"]


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid")


class NeutralSpec(_Strict):
    """A NEUTRAL organic (or weak electrolyte): the ``simulate_neutral`` inputs.

    Mirrors ``ui/sidebar._neutral_panel``: log Kow is the one required input; the
    air term, the weak-electrolyte speciation and the apoplastic bypass are opt-in
    and stay structurally absent unless switched on. ``grain_feed`` is
    ``simulate_neutral``'s own switch for the neutral brown rice ("xylem", the model
    default: a sap stream at the xylem-sap concentration, PREREG 9.15; None = not
    predicted)."""
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
    grain_feed: Optional[Literal["xylem"]] = "xylem"


class CompoundSpec(_Strict):
    kind: Literal["congener", "smiles", "neutral"] = "congener"
    congener: str = "PFOA"
    smiles: Optional[str] = None
    neutral: NeutralSpec = Field(default_factory=NeutralSpec)


class MechanismSpec(_Strict):
    """The crop drivers (every compound class) and the PFAS transport variant.

    ``crop`` is the Streamlit "Crop drivers" choice (``ui.common.CROPS``). ``shoot``,
    ``fxy_source`` and ``lipid_loading`` pick the fitted variant (``model_api.variant_key``);
    a combination the model has no fit for (f_xy 'fit' with lipid loading or with the
    4-pool shoot) runs 'recommended' with a message, as the Streamlit sidebar does. They are
    PFAS-only and ignored on the neutral path."""
    crop: Crop = "site2020"
    shoot: Literal["redistributed", "4pool"] = api.DEFAULT_SHOOT
    fxy_source: Literal["fit", "recommended"] = api.DEFAULT_F_XY_SOURCE
    lipid_loading: bool = api.DEFAULT_LIPID_LOADING


class TablesSpec(_Strict):
    """Editable growth + pore-water tables, column-major ({column: [values]})."""
    growth: Optional[Dict[str, List[Optional[float]]]] = None   # day, root, stem, leaf, grain
    growth_units: Literal["g/hill", "kg/hill", "g/m2", "kg/ha", "t/ha"] = "g/hill"
    cwo: Optional[Dict[str, List[Optional[float]]]] = None      # day, Cwo


class BiomonSpec(_Strict):
    conc: Dict[str, float] = Field(default_factory=lambda: {
        "root": 0.49, "straw": 0.83, "grain": 0.46})
    Cwo: Optional[float] = Field(1.0, gt=0.0)


class ExposureSpec(_Strict):
    mode: Mode = "parametric"
    # parametric (and the tables-mode fallback constant)
    Cwo: float = Field(1.0, ge=0.0, le=1e6)
    season: float = Field(120.0, ge=90.0, le=160.0)          # generic crops only
    cwo_profile: Literal["constant", "flooded"] = "constant"
    # custom tables
    tables: TablesSpec = Field(default_factory=TablesSpec)
    # HYDRUS / CSV drivers (None -> the bundled example)
    csv_text: Optional[str] = Field(None, max_length=5_000_000)
    # live HYDRUS-1D: the soil organic carbon (None -> the site's, soil_hydrus.FOC)
    f_oc: Optional[float] = Field(None, ge=0.005, le=0.05)
    # soil inventory -> pore water (None -> the Streamlit seeds: the Yamazaki pre-season
    # soil of SI Table S16 and the site Kd; `inventory_defaults`)
    C_total: Optional[float] = Field(None, ge=0.0, le=1e4)
    K_F: Optional[float] = Field(None, ge=1e-4, le=1e5)
    n_F: float = Field(1.0, ge=0.6, le=1.0)
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
    """``rsd_percent`` is the user's analytical relative SD (0 = not given); the model's own
    transfer error is always added inside ``estimate_exposure_bayesian``. ``view`` picks the
    mechanism the inverse runs under, as the Streamlit tabs do: "expert" passes the shoot and
    lipid switches, "simple" runs the crop and f_xy source only."""
    scenario: Scenario = Field(default_factory=Scenario)
    measured: Dict[str, float]
    rsd_percent: float = Field(0.0, ge=0.0, le=300.0)
    view: Literal["simple", "expert"] = "expert"
    lang: Lang = "en"


class DownloadRequest(_Strict):
    scenario: Scenario = Field(default_factory=Scenario)
    lang: Lang = "en"
    scenario_label: Optional[str] = None


class CwoPreviewRequest(_Strict):
    compound: CompoundSpec = Field(default_factory=CompoundSpec)
    crop: Crop = "site2020"
    level: float = Field(1.0, ge=0.0)
    season: float = Field(120.0, ge=90.0, le=160.0)        # generic crops only
