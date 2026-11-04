// Shapes shared with the FastAPI backend (server/schemas.py + server/service.py).

export type Lang = "ko" | "en";
export type Tissue = "root" | "stem" | "leaf" | "grain";
export const TISSUES: Tissue[] = ["root", "stem", "leaf", "grain"];
export type Mode = "parametric" | "tables" | "csv" | "hydrus" | "inventory" | "biomon";
export type CompoundKind = "congener" | "smiles" | "neutral";
export type Table = Record<string, (number | null)[]>;

export interface NeutralSpec {
  log_kow: number;
  name: string;
  half_life: number;
  tscf_model: "briggs" | "schriever";
  lipid_source: "measured" | "briggs_anchor";
  phloem: boolean;
  air: boolean;
  MW: number;
  K_AW: number;
  C_air: number;
  weak_electrolyte: boolean;
  pKa: number;
  pH: number;
  is_acid: boolean;
  g_apo: number;
  Koc: number | null;
  grain_feed: "xylem" | null;     // simulate_neutral's brown-rice feed (null = not predicted)
}

export interface CompoundSpec {
  kind: CompoundKind;
  congener: string;
  smiles: string | null;
  neutral: NeutralSpec;
}

// ui.common.CROPS keys: the calibration site-year (default) or a generic-climatology driver
export type Crop = "site2020" | "oryza2000" | "oryza" | "growth_rice";
export type Shoot = "redistributed" | "4pool";
export type FxySource = "fit" | "recommended";

export interface MechanismSpec {
  crop: Crop;
  shoot: Shoot;
  fxy_source: FxySource;
  lipid_loading: boolean;
}

export interface TablesSpec {
  growth: Table | null;
  growth_units: "g/hill" | "kg/hill" | "g/m2" | "kg/ha" | "t/ha";
  cwo: Table | null;
}

export interface ExposureSpec {
  mode: Mode;
  Cwo: number;
  season: number;                  // generic crops only (the site-year sets its own)
  cwo_profile: "constant" | "flooded";
  tables: TablesSpec;
  csv_text: string | null;
  f_oc: number | null;             // null -> the site's (soil_hydrus.FOC)
  C_total: number | null;          // null -> the Yamazaki pre-season soil (server seed)
  K_F: number | null;              // null -> the site Kd (server seed)
  n_F: number;
  biomon: { conc: Record<string, number>; Cwo: number | null } | null;
  biomon_csv_text: string | null;
}

export interface Scenario {
  compound: CompoundSpec;
  mechanism: MechanismSpec;
  exposure: ExposureSpec;
  lang: Lang;
}

export const DEFAULT_NEUTRAL: NeutralSpec = {
  log_kow: 2.45, name: "carbamazepine", half_life: 7.0, tscf_model: "briggs",
  lipid_source: "measured", phloem: false, air: false, MW: 236.3, K_AW: 1e-3, C_air: 0,
  weak_electrolyte: false, pKa: 4.5, pH: 6.5, is_acid: true, g_apo: 0, Koc: null, grain_feed: "xylem",
};

/** The default scenario. The transport variant and the crop come from the server's
 *  `model_api` defaults (`/api/meta`), never from a value typed here. */
export function defaultScenario(lang: Lang, meta: Meta): Scenario {
  return {
    compound: { kind: "congener", congener: "PFOA", smiles: null, neutral: { ...DEFAULT_NEUTRAL } },
    mechanism: {
      crop: meta.default_crop, shoot: meta.default_shoot, fxy_source: meta.default_f_xy_source,
      lipid_loading: meta.default_lipid_loading,
    },
    exposure: {
      mode: "parametric", Cwo: 1, season: meta.generic_season, cwo_profile: "constant",
      tables: { growth: null, growth_units: "g/hill", cwo: null },
      csv_text: null, f_oc: null, C_total: null, K_F: null, n_F: 1.0,
      biomon: null, biomon_csv_text: null,
    },
    lang,
  };
}

export interface Band { value: number; lo: number | null; hi: number | null; factor: number }
export interface Intake {
  grain_conc: number; rice_intake_g_day: number; body_weight_kg: number;
  weekly_intake_ng_per_kg_bw: number | null; twi: number; fraction: number | null;
  percent: number | null; in_group: boolean | null; signal: "green" | "amber" | "red" | null;
}
export interface Message { level: "info" | "success" | "warning" | "error"; text: string }

export interface Descriptors {
  formula: string; mol_weight: number; n_perfluoroC: number; head_group: string;
  n_ether_O: number; n_CF3: number; branched: boolean; is_linear: boolean;
  matched_name: string | null; notes: string[];
}

export interface RunResult {
  congener: string;
  kind: CompoundKind;
  mode: Mode;
  success: boolean;
  t: number[];
  conc: Record<Tissue, (number | null)[]>;
  baf: Record<Tissue, (number | null)[]>;
  straw: number[];
  Cwo: number[];
  Qtp: number[];
  M: Record<Tissue, number[]>;
  baf_final: Record<Tissue, number>;
  straw_baf: number;
  cwo_ref: number;
  baf_final_dw: Record<string, number> | null;     // dry basis (the Yamazaki comparison)
  conc_dw_final: Record<string, number>;
  B_k: Record<Tissue, number> | null;
  params: Record<string, any>;
  variant: { key: string; withdrawn: string | null } | null;
  finals: { root: number; straw: number; grain: number; white: number };
  bands: { root: Band; straw: Band; grain: Band; white: Band };
  input_range: Record<"root" | "straw" | "grain", [number, number] | null>;
  uncertainty_factor: number;
  intake: Intake;
  interp_spread: Record<"root" | "straw" | "grain", number> | null;
  schematic: Record<"conc" | "baf", { cmin: number; cmax: number; label: string }>;
  obs: Record<string, number>;
  bio_baf: Record<string, number> | null;
  measured_bio: { conc: Record<string, number>; Cwo: number | null } | null;
  descriptors: Descriptors | null;
  provisional: boolean;
  messages: Message[];
  has_isotherm: boolean;
  inventory: { C_total: number; K_F: number; n_F: number } | null;
  effective: { Cwo: number; season: number; cwo_profile: string; drivers: boolean;
               sim_kw: Record<string, any> | null };
  // neutral only
  TSCF?: number; K_PW?: Record<Tissue, number>; rcf_briggs?: number; half_life?: number | null;
  warning?: string | null; air_summary?: any; koc?: number; grain_predicted?: boolean;
}

export type Strings = Record<string, Record<string, string>>;

export interface Meta {
  app_version: string;
  repo_url: string;
  docs_url: string;
  congeners: string[];
  congener_info: Record<string, { n_C: number; group: string }>;
  labels: Record<Lang, Record<string, string>>;
  presets: { en: [string, number][]; ko: [string, number, string][] };
  scenarios_ko: [string, string, number, string, string][];
  strings: Strings;
  crops: [Crop, string][];
  default_crop: Crop;
  site_season: number;
  generic_season: number;
  shoots: Shoot[];
  default_shoot: Shoot;
  default_f_xy_source: FxySource;
  default_lipid_loading: boolean;
  default_variant: string;
  withdrawn_variants: Record<string, string>;
  tang_congeners: string[];
  example_smiles: [string, string][];
  intake: { rice_g: number; rice_range: [number, number]; bw: number; twi: number; congeners: string[] };
  uncertainty_factor: number;
  apriori_log10_rmse: number;
  neutral_scores: Record<string, number>;
  neutral_equilibrium_scores: Record<string, number>;
  neutral_lipid_fw: { measured: number; briggs_anchor: number };
  karickhoff: { slope: number; intercept: number };
  koc_max: number;
  f_n_tested_floor: number;
  soil_foc: number;
  hydrus_foc: number;
  available: { rdkit: boolean; hydrus: boolean; ctx: boolean };
}

export interface FigurePayload { figure: { data: any[]; layout: any }; data: Record<string, any> }

export interface Estimate {
  congener: string; median: number; ci68: [number | null, number | null];
  ci95: [number | null, number | null]; sd_log10: number | null; used_tissues: string[];
  measured: Record<string, number>; model_fit: Record<string, number>;
}

export interface CtxResult {
  ok: boolean; note: string; dtxsid: string; name: string; casrn: string; smiles: string;
  kwargs: Record<string, any>; badges: Record<string, string>; sources: Record<string, string>;
  outside_ad: Record<string, boolean>; pka_measured: boolean;
}
