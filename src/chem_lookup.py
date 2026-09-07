"""Chemical identity + physicochemical property lookup (EPA CompTox / CTX API).

WHAT THIS IS FOR
----------------
The NEUTRAL-organic path needs a **log Kow** and, for its optional terms, a molar
mass, a Henry's-law constant, a pKa and a soil Koc. Those are properties of the
compound, not of this model, and typing them by hand is where a user's error
enters. This module resolves a **name, CAS-RN or SMILES** to an EPA DTXSID and
pulls those values from the CompTox Chemicals Dashboard's CTX API.

THE ONE RULE THIS MODULE ENFORCES: PROVENANCE TRAVELS WITH THE VALUE
-------------------------------------------------------------------
CompTox serves EXPERIMENTAL and PREDICTED (OPERA) values side by side, and the
neutral path's only claim is that *nothing in it is fitted* -- every published
a-priori number (Liu 0.206/0.281, Ge 0.783, Briggs-stem 0.299) was computed on a
measured/reported log Kow. Silently substituting an OPERA prediction would leave
those numbers describing something the app no longer does. So every value here is
a `Prop` carrying `source` ("experimental" / "predicted") and the reporting lab or
model, `experimental` is PREFERRED when both exist, and the UI is expected to show
the badge. Nothing is silently upgraded to fact.

NOT AVAILABLE HERE, on purpose: the **in-planta half-life**. No dashboard field
corresponds to it -- and Kodesova 2019 measured the parent fraction varying 4.8x
BETWEEN SPECIES for one compound, so it is not a compound property to look up
(docs/neutral_dpu_validation.md section 4i).

NETWORK / KEY
-------------
The CTX API needs a free EPA key sent as `x-api-key`. It is read from the
`CTX_API_KEY` environment variable or, inside Streamlit, `st.secrets["ctx_api_key"]`
-- never committed. With no key, no network, or an unknown compound this module
returns an empty `Lookup` with a `note` saying which of those happened; it never
raises into the caller and the app stays fully usable with manual entry.

CLI (run it on a machine that has the key and outbound HTTPS):

    python src/chem_lookup.py carbamazepine
    python src/chem_lookup.py 298-46-4
    python src/chem_lookup.py "NC(=O)N1c2ccccc2C=Cc2ccccc21"     # SMILES (needs RDKit)
    python src/chem_lookup.py carbamazepine --raw                # dump raw JSON

`--raw` dumps every endpoint's payload, because EPA revises paths and field names
and a mismatch produces EMPTY COLUMNS, not an error: probe one chemical and read
the raw JSON before trusting a batch (the working notes' first rule). The parsers
match property names case-insensitively, take alias lists for every field, and
accept a bare list or a `{"data": [...]}` envelope, so a schema change degrades to
a MISSING property rather than a wrong one.

Verify the whole pipeline on a compound with well-established properties before
trusting it -- benzoic acid and benzyl alcohol are the working notes' suggestions.

BASE URL: `https://comptox.epa.gov/ctx-api`. The older `api-ccte.epa.gov` no
longer resolves -- do not use it. `CTX_BASE_URL` overrides.
"""
from __future__ import annotations

import json
import os
import re
import sys
from dataclasses import dataclass, field, asdict

DEFAULT_BASE_URL = "https://comptox.epa.gov/ctx-api"
TIMEOUT_S = 20.0
RETRY_STATUS = (429, 502, 503, 504)      # transient; anything else is fatal for the record
RETRIES = 3

# The five CTX traps this module is written against (working notes, verified
# against the live API by the repo owner -- every one of them fails SILENTLY):
#  1. Batch bodies are newline-separated PLAIN TEXT, not JSON arrays (a JSON array
#     returns HTTP 200 with every field null). Not used here -- single lookups only.
#  2. "NaN"/"inf"/"null"/"N/A" arrive as STRINGS, and float("NaN") SUCCEEDS, so a
#     missing experimental value would outrank and displace a good prediction.
#     `_num` rejects them by name before converting (`_NULL_STRINGS`).
#  3. Field naming differs BETWEEN endpoints: the property endpoints are camelCase
#     (`propName`/`propValue`/`propUnit`/`sourceName`), records nested in the fate
#     endpoint are snake_case (`prop_name`/`prop_value`/...). Every reader below
#     takes an ALIAS LIST rather than one spelling.
#  4. NOTHING in a property payload says whether a value is measured or predicted
#     -- only the URL you called does. So provenance is tagged AT CALL TIME from
#     the endpoint, never inferred from a field that may be absent. (The fate
#     endpoint is the one exception: it mixes both and carries `propType`.)
#  5. Units are not what you expect -- Henry's law comes as atm-m3/mol. Reading
#     Pa-m3/mol as atm-m3/mol is a clean log10(101325) = 5.006 log-unit offset that
#     still looks like a plausible number, so `henry_to_kaw` converts explicitly and
#     REFUSES an unrecognised unit.

# Henry's law: the dashboard reports H in atm-m3/mol; the model wants the
# DIMENSIONLESS air-water partition K_AW = H / (R*T)  (R in atm-m3/(mol*K)).
R_ATM_M3 = 8.20573660809596e-5
T_REF_K = 298.15
_RT = R_ATM_M3 * T_REF_K                       # 0.024465 atm-m3/mol at 25 C

# Values that must never become numbers even though float() would accept them
# (trap 2). Compared lower-cased against the raw string.
_NULL_STRINGS = {"nan", "inf", "-inf", "infinity", "-infinity", "none", "null", "n/a", "na", ""}

# Field-name ALIASES (trap 3): camelCase on the property endpoints, snake_case
# inside the fate endpoint's nested records.
_NAME_KEYS = ("propName", "prop_name", "propertyId", "propertyName", "name", "property")
_VALUE_KEYS = ("propValue", "prop_value", "value", "resultValue", "medianValue", "meanValue")
_UNIT_KEYS = ("propUnit", "prop_unit", "unit", "units")
_SOURCE_KEYS = ("sourceName", "source_name", "source", "modelName", "model_name", "dataSource")
_TYPE_KEYS = ("propType", "prop_type", "propertyType", "type")
# OPERA applicability domain: the per-model fields come back null in practice; the
# populated ones are the Global variants. Outside the AD a prediction is not just
# uncertain, it is out of scope -- so it reaches the badge.
_AD_KEYS = ("adConclusionGlobal", "ad_conclusion_global", "adConclusion", "ad_conclusion")

# What each model input is called on the dashboard. Matched case-insensitively as
# SUBSTRINGS of the returned property name, most specific first, because the CTX
# property vocabulary is not stable enough to hard-code exact ids.
_PROPERTY_PATTERNS = {
    "log_kow": ("octanol-water partition", "logkow", "log kow", "logp", "log p"),
    "henry": ("henry",),
    "pka_acidic": ("pka_a", "acidic pka", "pka (acidic", "strongest acidic"),
    "pka_basic": ("pka_b", "basic pka", "pka (basic", "strongest basic"),
    "log_koc": ("koc",),
    "water_solubility": ("water solubility",),
    "melting_point": ("melting point",),
}


@dataclass
class Prop:
    """One looked-up property value, with where it came from.

    `source` is "experimental" or "predicted"; `origin` names the lab/model when
    the API reports one (e.g. "OPERA"). `raw_unit` is kept so a caller can see
    what was converted.
    """
    value: float
    source: str = "unknown"
    origin: str = ""
    raw_value: float | None = None
    raw_unit: str = ""
    ad: str = ""                    # OPERA applicability domain (Global conclusion)

    @property
    def is_experimental(self) -> bool:
        return self.source == "experimental"

    @property
    def outside_ad(self) -> bool:
        """A prediction the model itself says is outside its applicability domain."""
        a = (self.ad or "").strip().lower()
        return bool(a) and ("outside" in a or a.startswith("no"))

    def badge(self) -> str:
        """Short provenance label for a UI ('experimental (PHYSPROP)',
        'predicted (OPERA)', 'predicted (OPERA, OUTSIDE the applicability domain)')."""
        s = self.source if self.source != "unknown" else "source unknown"
        bits = [b for b in (self.origin,
                            "OUTSIDE the applicability domain" if self.outside_ad else "") if b]
        return f"{s} ({', '.join(bits)})" if bits else s


@dataclass
class Lookup:
    """Everything a query resolved to. `ok` is False when nothing was found."""
    query: str = ""
    ok: bool = False
    note: str = ""
    dtxsid: str = ""
    preferred_name: str = ""
    casrn: str = ""
    smiles: str = ""
    inchikey: str = ""
    mol_weight: float | None = None
    props: dict = field(default_factory=dict)          # name -> Prop

    def get(self, name):
        return self.props.get(name)

    def value(self, name, default=None):
        p = self.props.get(name)
        return default if p is None else p.value

    def to_dict(self):
        d = asdict(self)
        d["props"] = {k: asdict(v) for k, v in self.props.items()}
        return d

    # -- the model-facing view -------------------------------------------------
    def neutral_kwargs(self):
        """The subset that maps onto `model_api.simulate_neutral(...)` arguments.

        `K_AW` is converted from the Henry's-law constant; `Koc` feeds the soil
        modes; `pKa`/`is_acid` are filled only when a pKa was actually found (a
        strictly neutral compound must keep `pKa=None`, which is the bit-identical
        neutral path). The in-planta HALF-LIFE is deliberately absent -- it is not
        a dashboard property and not a compound constant.
        """
        out = {}
        if "log_kow" in self.props:
            out["log_kow"] = self.props["log_kow"].value
        if self.mol_weight:
            out["MW"] = float(self.mol_weight)
        if "K_AW" in self.props:
            out["K_AW"] = self.props["K_AW"].value
        if "log_koc" in self.props:
            out["Koc"] = float(10.0 ** self.props["log_koc"].value)
        acid, base = self.props.get("pka_acidic"), self.props.get("pka_basic")
        if acid is not None or base is not None:
            # Which centre governs is a judgement the DATA cannot make: report the
            # one that is defined, and prefer the ACIDIC one when both are, since
            # that is the case the ionisable-organic extension was built for.
            out["pKa"] = (acid or base).value
            out["is_acid"] = acid is not None
        return out


# ---------------------------------------------------------------------------
# key / transport
# ---------------------------------------------------------------------------
def api_key(explicit=None):
    """CTX API key from (1) the argument, (2) $CTX_API_KEY, (3) st.secrets.

    Returns "" when unset -- callers degrade to manual entry rather than failing.
    """
    if explicit:
        return str(explicit).strip()
    env = os.environ.get("CTX_API_KEY", "").strip()
    if env:
        return env
    try:                                                # only inside Streamlit
        import streamlit as st
        return str(st.secrets.get("ctx_api_key", "")).strip()
    except Exception:                                   # noqa: BLE001
        return ""


def base_url(explicit=None):
    """CTX base URL: the argument, else $CTX_BASE_URL, else EPA's.

    The override exists so the client can be pointed at a mirror or a local stub
    (that is how the success path is exercised where EPA is unreachable), not
    because the endpoint is expected to move."""
    return (explicit or os.environ.get("CTX_BASE_URL") or DEFAULT_BASE_URL).rstrip("/")


def _get(path, key, base_url=DEFAULT_BASE_URL, params=None):
    """One GET, with a short retry on TRANSIENT statuses. Returns (json, error);
    never raises. 429/502/503/504 are retried with exponential backoff; any other
    non-200 is fatal for this record rather than retried."""
    import time
    import requests
    url = f"{base_url.rstrip('/')}{path}"
    headers = {"x-api-key": key, "accept": "application/json",
               "content-type": "application/json"}
    last = "unknown error"
    for attempt in range(RETRIES):
        try:
            r = requests.get(url, headers=headers, params=params, timeout=TIMEOUT_S)
        except Exception as e:                          # noqa: BLE001 (network is optional)
            last = f"{type(e).__name__}: {e}"
            if attempt + 1 < RETRIES:
                time.sleep(2.0 ** attempt)
                continue
            return None, last
        if r.status_code in RETRY_STATUS and attempt + 1 < RETRIES:
            time.sleep(2.0 ** attempt)
            last = f"HTTP {r.status_code}"
            continue
        if r.status_code in (401, 403):
            return None, f"HTTP {r.status_code} — the CTX API rejected the key"
        if r.status_code == 404:
            return None, "not found (HTTP 404)"
        if r.status_code >= 400:
            return None, f"HTTP {r.status_code}"
        try:
            return r.json(), None
        except Exception as e:                          # noqa: BLE001
            return None, f"unparseable response: {type(e).__name__}: {e}"
    return None, last


def _rows(payload):
    """CTX returns either a bare list or an object wrapping one; take the rows."""
    if payload is None:
        return []
    if isinstance(payload, list):
        return [r for r in payload if isinstance(r, dict)]
    if isinstance(payload, dict):
        for k in ("data", "results", "content", "items"):
            v = payload.get(k)
            if isinstance(v, list):
                return [r for r in v if isinstance(r, dict)]
        return [payload]
    return []


# ---------------------------------------------------------------------------
# identity
# ---------------------------------------------------------------------------
_CASRN_RE = re.compile(r"^\d{2,7}-\d{2}-\d$")
_DTXSID_RE = re.compile(r"^DTXSID\d+$", re.I)
_INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")


def identifier_kind(query):
    """Classify a query as 'dtxsid' | 'casrn' | 'inchikey' | 'smiles' | 'name'.

    SMILES detection is deliberately conservative -- a string is only called a
    SMILES when RDKit parses it AND it does not look like a plain word, because a
    chemical NAME must never be mistaken for a structure.
    """
    q = (query or "").strip()
    if not q:
        return "name"
    if _DTXSID_RE.match(q):
        return "dtxsid"
    if _CASRN_RE.match(q):
        return "casrn"
    if _INCHIKEY_RE.match(q.upper()):
        return "inchikey"
    if re.fullmatch(r"[A-Za-z0-9 ,'\-\(\)\[\]\.]+", q) and not re.search(r"[=#@\\/]|\d\)", q):
        # letters/spaces only, no bond or ring syntax -> a name (e.g. "carbamazepine")
        looks_like_word = bool(re.fullmatch(r"[A-Za-z][A-Za-z \-',\.]*", q))
        if looks_like_word:
            return "name"
    try:
        from rdkit import Chem
        if Chem.MolFromSmiles(q) is not None:
            return "smiles"
    except Exception:                                   # noqa: BLE001 (RDKit optional)
        pass
    return "name"


def smiles_to_inchikey(smiles):
    """SMILES -> InChIKey (RDKit). The CTX chemical search has no SMILES route, so
    a structure query is resolved through its InChIKey. Returns "" if unavailable."""
    try:
        from rdkit import Chem
        m = Chem.MolFromSmiles(smiles)
        return Chem.MolToInchiKey(m) if m is not None else ""
    except Exception:                                   # noqa: BLE001
        return ""


def resolve(query, key=None, base=None, _get_fn=None):
    """Resolve a name / CAS-RN / InChIKey / SMILES / DTXSID to a `Lookup` identity.

    Returns a Lookup with `ok=False` and a `note` when the key is missing, the
    network is unreachable, or nothing matched -- the caller keeps working.
    """
    get = _get_fn or _get
    url = base_url(base)
    q = (query or "").strip()
    out = Lookup(query=q)
    if not q:
        out.note = "empty query"
        return out
    k = api_key(key)
    if not k:
        out.note = ("no CTX API key — set CTX_API_KEY or streamlit secrets "
                    "`ctx_api_key` (a free key comes from EPA)")
        return out

    kind = identifier_kind(q)
    term = q
    if kind == "smiles":
        term = smiles_to_inchikey(q)
        if not term:
            out.note = "could not turn that SMILES into an InChIKey (RDKit needed)"
            return out
    payload, err = get(f"/chemical/search/equal/{term}", k, url)
    rows = _rows(payload)
    if not rows and kind in ("smiles", "inchikey") and "-" in term:
        # fall back to the InChIKey SKELETON (first block): the dashboard entry may
        # differ in stereo/protonation from the structure the user drew
        payload, err = get(f"/chemical/search/start-with/{term.split('-')[0]}", k, url)
        rows = _rows(payload)
    if not rows:
        out.note = err or f"no CompTox match for {q!r}"
        return out

    row = rows[0]
    out.dtxsid = str(row.get("dtxsid") or row.get("dtxsId") or "")
    out.preferred_name = str(row.get("preferredName") or row.get("searchName") or "")
    out.casrn = str(row.get("casrn") or "")
    out.smiles = str(row.get("smiles") or "")
    out.inchikey = str(row.get("inchikey") or row.get("inchiKey") or "")
    if not out.dtxsid:
        out.note = "a match came back without a DTXSID"
        return out
    out.ok = True
    if len(rows) > 1:
        out.note = f"{len(rows)} matches; using {out.preferred_name or out.dtxsid}"
    return out


# ---------------------------------------------------------------------------
# properties
# ---------------------------------------------------------------------------
def _num(x):
    """Float, or None -- rejecting the STRING nulls first (trap 2).

    `float("NaN")` and `float("inf")` both SUCCEED, so a `propValue` of "NaN" would
    otherwise become a real number and, being on the experimental endpoint, would
    outrank and displace a perfectly good prediction."""
    if x is None:
        return None
    if isinstance(x, str) and x.strip().lower() in _NULL_STRINGS:
        return None
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if v == v and abs(v) != float("inf") else None      # belt and braces


def _first(row, keys):
    for k in keys:
        if row.get(k) not in (None, ""):
            return row[k]
    return None


def _row_source(row, default=""):
    """Provenance from the PAYLOAD -- only valid for the fate endpoint (trap 4).

    Everywhere else the payload says nothing about measured-vs-predicted and the
    caller must pass the bucket the ENDPOINT implies; `default` is that bucket."""
    t = str(_first(row, _TYPE_KEYS) or "").lower()
    if "exp" in t:
        return "experimental"
    if "pred" in t or "opera" in t:
        return "predicted"
    return default


def _match_property(row, patterns):
    name = " ".join(str(row.get(k, "")) for k in _NAME_KEYS).lower()
    return any(p in name for p in patterns)


def _pick(tagged, patterns):
    """Best (row, source) for a property. EXPERIMENTAL wins over predicted; a row
    whose value is missing or a string null is not a candidate at all (trap 2), so a
    blank measurement can never displace a usable prediction. Among predictions, one
    INSIDE its applicability domain beats one outside."""
    hits = [(r, src) for r, src in tagged
            if _match_property(r, patterns) and _num(_first(r, _VALUE_KEYS)) is not None]
    if not hits:
        return None, ""
    def rank(item):
        r, src = item
        ad_bad = 1 if str(_first(r, _AD_KEYS) or "").lower().startswith(("outside", "no")) else 0
        return (0 if src == "experimental" else 1, ad_bad)
    hits.sort(key=rank)
    return hits[0]


def _value_of(row):
    return _first(row, _VALUE_KEYS)


def _unit_of(row):
    return str(_first(row, _UNIT_KEYS) or "")


def henry_to_kaw(value, unit=""):
    """Henry's-law constant -> the DIMENSIONLESS K_AW the model uses.

    atm-m3/mol (the dashboard's usual unit) -> H/(R*T); Pa-m3/mol -> /101325 first;
    an already-dimensionless value passes through. Returns None for a unit this
    does not recognise, rather than guessing -- a wrong K_AW silently switches the
    leaf's volatilisation sink on or off.
    """
    v = _num(value)
    if v is None:
        return None
    u = (unit or "").lower().replace(" ", "")
    if not u or "atm" in u:                       # atm-m3/mol (default reading)
        return v / _RT
    if "pa" in u:                                 # Pa-m3/mol
        return (v / 101325.0) / _RT
    if "dimensionless" in u or u in ("-", "unitless"):
        return v
    return None


def properties(dtxsid, key=None, base=None, _get_fn=None):
    """Physicochemical properties for a DTXSID -> ({name: Prop}, note).

    PROVENANCE COMES FROM THE ENDPOINT, NOT THE PAYLOAD (trap 4). A property record
    carries nothing that says whether it was measured or modelled, so this calls the
    experimental and predicted paths SEPARATELY and tags each batch of rows with the
    bucket its URL implies. The fate endpoint is the one exception -- it mixes both
    and does carry `propType`/`prop_type` -- so its rows are read with that field and
    fall back to the endpoint tag only when it is absent.

    `K_AW` is derived from the Henry's-law row (see `henry_to_kaw`); everything else
    is passed through with its unit recorded.
    """
    get = _get_fn or _get
    url = base_url(base)
    k = api_key(key)
    if not k or not dtxsid:
        return {}, ("no CTX API key" if not k else "no DTXSID")

    tagged, errs = [], []                      # [(row, source_from_the_url), ...]
    for path, tag in (
            (f"/chemical/property/experimental/search/by-dtxsid/{dtxsid}", "experimental"),
            (f"/chemical/property/predicted/search/by-dtxsid/{dtxsid}", "predicted"),
            (f"/chemical/fate/search/by-dtxsid/{dtxsid}", "")):
        payload, err = get(path, k, url)
        if err:
            errs.append(f"{path.rsplit('/', 2)[0].rsplit('/', 1)[-1]}: {err}")
            continue
        for row in _rows(payload):
            # the fate payload nests its records under experimental/predicted keys
            nested = False
            for nk, ntag in (("experimentalFateData", "experimental"),
                             ("experimental_fate_data", "experimental"),
                             ("predictedFateData", "predicted"),
                             ("predicted_fate_data", "predicted")):
                inner = row.get(nk)
                if isinstance(inner, list):
                    nested = True
                    tagged += [(r, _row_source(r, ntag)) for r in inner if isinstance(r, dict)]
            if not nested:
                tagged.append((row, _row_source(row, tag)))

    props = {}
    for name, patterns in _PROPERTY_PATTERNS.items():
        row, src = _pick(tagged, patterns)
        if row is None:
            continue
        v = _num(_first(row, _VALUE_KEYS))
        origin = str(_first(row, _SOURCE_KEYS) or ("OPERA" if src == "predicted" else ""))
        props[name] = Prop(value=v, source=src or "unknown", origin=origin,
                           raw_value=v, raw_unit=_unit_of(row),
                           ad=str(_first(row, _AD_KEYS) or ""))

    h = props.get("henry")
    if h is not None:
        kaw = henry_to_kaw(h.raw_value, h.raw_unit)
        if kaw is not None:
            props["K_AW"] = Prop(value=kaw, source=h.source, origin=h.origin,
                                 raw_value=h.raw_value, raw_unit=h.raw_unit or "atm-m3/mol",
                                 ad=h.ad)
    note = "" if props else ("; ".join(errs) or "no properties returned")
    return props, note


def lookup(query, key=None, base=None, _get_fn=None):
    """Resolve `query` and fetch its properties in one call. Never raises."""
    out = resolve(query, key=key, base=base, _get_fn=_get_fn)
    if not out.ok:
        return out
    props, note = properties(out.dtxsid, key=key, base=base, _get_fn=_get_fn)
    out.props = props
    detail, _ = (_get_fn or _get)(f"/chemical/detail/search/by-dtxsid/{out.dtxsid}",
                                  api_key(key), base_url(base))
    for row in _rows(detail):
        mw = _num(row.get("averageMass") or row.get("molWeight") or row.get("monoisotopicMass"))
        if mw:
            out.mol_weight = mw
        out.smiles = out.smiles or str(row.get("smiles") or "")
        break
    if note and not out.note:
        out.note = note
    return out


# ---------------------------------------------------------------------------
def _cli(argv):
    args = [a for a in argv if not a.startswith("--")]
    raw = "--raw" in argv
    if not args:
        print(__doc__.strip().split("CLI (")[-1])
        return 2
    q = " ".join(args)
    k = api_key()
    if not k:
        print("No CTX API key. Set it first:\n"
              "  export CTX_API_KEY='...'        # or streamlit secrets: ctx_api_key\n"
              "A free key is issued by EPA for the CompTox (CTX) APIs.")
        return 1
    if raw:
        base = base_url()
        print(f"# base {base}")
        payload, err = _get(f"/chemical/search/equal/{q}", k, base)
        print(f"--- /chemical/search/equal/{q}\n{err or json.dumps(payload, indent=2)[:3000]}")
        rows = _rows(payload)
        if rows:
            sid = rows[0].get("dtxsid") or rows[0].get("dtxsId")
            for path in (f"/chemical/detail/search/by-dtxsid/{sid}",
                         f"/chemical/property/experimental/search/by-dtxsid/{sid}",
                         f"/chemical/property/predicted/search/by-dtxsid/{sid}",
                         f"/chemical/fate/search/by-dtxsid/{sid}"):
                pl, e = _get(path, k, base)
                print(f"\n--- {path}\n{e or json.dumps(pl, indent=2)[:4000]}")
        return 0
    r = lookup(q)
    if not r.ok:
        print(f"✗ {r.note}")
        return 1
    print(f"{r.preferred_name or q}   {r.dtxsid}   CAS {r.casrn or '—'}")
    print(f"  SMILES {r.smiles or '—'}")
    print(f"  MW     {r.mol_weight if r.mol_weight else '—'}")
    for name, p in sorted(r.props.items()):
        unit = f" [{p.raw_unit}]" if p.raw_unit and name != "K_AW" else ""
        print(f"  {name:16} {p.value:<14.6g}{unit:14} {p.badge()}")
    nk = r.neutral_kwargs()
    print("\n  -> simulate_neutral(" + ", ".join(f"{k_}={v!r}" for k_, v in nk.items()) + ")")
    print("  NOTE: the in-planta half-life is NOT a dashboard property — set it yourself.")
    print("  Sanity-check the pipeline on a well-characterised compound before trusting a")
    print("  batch: 'benzoic acid' and 'benzyl alcohol' both have solid measured values.")
    if r.note:
        print(f"  note: {r.note}")
    return 0


if __name__ == "__main__":                              # pragma: no cover
    sys.exit(_cli(sys.argv[1:]))
