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
}

export interface CompoundSpec {
  kind: CompoundKind;
  congener: string;
  smiles: string | null;
  neutral: NeutralSpec;
}

export interface MechanismSpec {
  E_m: number;
  fxy_source: "recommended" | "W2fit";
  biomass: "oryza" | "growth_rice";
  uptake: "carrier" | "bypass";
  lipid_loading: boolean;
  vmax_scale: number | null;
  g_apo: number | null;
  km_scale: number;
}

export interface TablesSpec {
  growth: Table | null;
  growth_units: "g/hill" | "kg/hill" | "g/m2" | "kg/ha" | "t/ha";
  cwo: Table | null;
  density: Record<Tissue, number>;
}

export interface ExposureSpec {
  mode: Mode;
  Cwo: number;
  season: number;
  measured_forcing: boolean;
  cwo_profile: "constant" | "flooded";
  k_leach: number | null;
  tables: TablesSpec;
  csv_text: string | null;
  f_oc: number;
  flood_until: number;
  percolation: number;
  C_total: number;
  K_F: number;
  n_F: number;
  theta_g: number;
  inv_flood: boolean;
  inv_k_leach: number | null;
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
  weak_electrolyte: false, pKa: 4.5, pH: 6.5, is_acid: true, g_apo: 0, Koc: null,
};

export function defaultScenario(lang: Lang): Scenario {
  return {
    compound: { kind: "congener", congener: "PFOA", smiles: null, neutral: { ...DEFAULT_NEUTRAL } },
    mechanism: {
      E_m: -120, fxy_source: "recommended", biomass: "oryza", uptake: "carrier",
      lipid_loading: false, vmax_scale: null, g_apo: null, km_scale: 1,
    },
    exposure: {
      mode: "parametric", Cwo: 1, season: 120, measured_forcing: true, cwo_profile: "constant",
      k_leach: null,
      tables: { growth: null, growth_units: "g/hill", cwo: null,
                density: { root: 1.0, stem: 0.3, leaf: 0.3, grain: 1.2 } },
      csv_text: null, f_oc: 0.02, flood_until: 90, percolation: 0.3,
      C_total: 5, K_F: 2, n_F: 0.85, theta_g: 0.35, inv_flood: true, inv_k_leach: null,
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
  B_k: Record<Tissue, number> | null;
  N: number | null;
  eN: number | null;
  params: Record<string, any>;
  finals: { root: number; straw: number; grain: number };
  bands: { root: Band; straw: Band; grain: Band };
  uncertainty_factor: number;
  intake: Intake;
  schematic: Record<"conc" | "baf", { cmin: number; cmax: number; label: string }>;
  obs: Record<string, number>;
  bio_baf: Record<string, number> | null;
  measured_bio: { conc: Record<string, number>; Cwo: number | null } | null;
  descriptors: Descriptors | null;
  provisional: boolean;
  messages: Message[];
  implied_volume_ml: Record<Tissue, number> | null;
  has_isotherm: boolean;
  effective: { Cwo: number; season: number; measured: boolean; cwo_profile: string;
               k_leach: number | null; drivers: boolean };
  // neutral only
  TSCF?: number; K_PW?: Record<Tissue, number>; rcf_briggs?: number; half_life?: number | null;
  warning?: string | null; air_summary?: any; koc?: number;
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
  organ_labels: Record<Lang, Record<Tissue, string>>;
  uptake_modes: Record<string, { vmax_scale: number; g_apo: number }>;
  default_uptake: string;
  tang_congeners: string[];
  tang_refit_fxy: Record<string, number>;
  example_smiles: [string, string][];
  tissue_density: Record<Tissue, number>;
  intake: { rice_g: number; bw: number; twi: number; congeners: string[] };
  uncertainty_factor: number;
  default_k_leach: Record<string, number>;
  karickhoff: { slope: number; intercept: number };
  koc_max: number;
  f_n_tested_floor: number;
  inverse_precision: Record<Lang, [string, number][]>;
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
