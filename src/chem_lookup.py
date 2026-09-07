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

`--raw` exists because EPA's response SHAPE is what this module has to guess at:
the parsers below match property names case-insensitively and accept both a bare
list and a `{"data": [...]}` envelope, so a schema change degrades to "not found"
rather than a crash. If a lookup comes back empty for a compound you know is in
the dashboard, `--raw` shows what actually arrived.
"""
from __future__ import annotations

import json
import os
import re
import sys
from dataclasses import dataclass, field, asdict

DEFAULT_BASE_URL = "https://api-ccte.epa.gov"
TIMEOUT_S = 20.0

# Henry's law: the dashboard reports H in atm-m3/mol; the model wants the
# DIMENSIONLESS air-water partition K_AW = H / (R*T)  (R in atm-m3/(mol*K)).
R_ATM_M3 = 8.20573660809596e-5
T_REF_K = 298.15
_RT = R_ATM_M3 * T_REF_K                       # 0.024465 atm-m3/mol at 25 C

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

    @property
    def is_experimental(self) -> bool:
        return self.source == "experimental"

    def badge(self) -> str:
        """Short provenance label for a UI ('experimental' / 'predicted (OPERA)')."""
        s = self.source if self.source != "unknown" else "source unknown"
        return f"{s} ({self.origin})" if self.origin else s


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
    """One GET. Returns (parsed_json, error_string); never raises."""
    import requests
    url = f"{base_url.rstrip('/')}{path}"
    try:
        r = requests.get(url, headers={"x-api-key": key, "accept": "application/json"},
                         params=params, timeout=TIMEOUT_S)
    except Exception as e:                              # noqa: BLE001 (network is optional)
        return None, f"{type(e).__name__}: {e}"
    if r.status_code == 401 or r.status_code == 403:
        return None, f"HTTP {r.status_code} — the CTX API rejected the key"
    if r.status_code == 404:
        return None, "not found (HTTP 404)"
    if r.status_code >= 400:
        return None, f"HTTP {r.status_code}"
    try:
        return r.json(), None
    except Exception as e:                              # noqa: BLE001
        return None, f"unparseable response: {type(e).__name__}: {e}"


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
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if v == v and abs(v) != float("inf") else None


def _classify(row):
    """(source, origin) for a property row: experimental vs predicted, and by whom."""
    blob = " ".join(str(row.get(k, "")) for k in
                    ("propType", "propertyType", "type", "source", "modelName", "dataSource")).lower()
    if "exp" in blob:
        return "experimental", str(row.get("source") or row.get("dataSource") or "")
    if "pred" in blob or "opera" in blob or "model" in blob:
        return "predicted", str(row.get("modelName") or row.get("source") or
                                ("OPERA" if "opera" in blob else ""))
    return "unknown", str(row.get("source") or "")


def _match_property(row, patterns):
    name = " ".join(str(row.get(k, "")) for k in
                    ("propertyId", "name", "propertyName", "property")).lower()
    return any(p in name for p in patterns)


def _pick(rows, patterns):
    """Best row for a property: EXPERIMENTAL wins over predicted, else first found."""
    hits = [r for r in rows if _match_property(r, patterns) and _num(_value_of(r)) is not None]
    if not hits:
        return None
    hits.sort(key=lambda r: 0 if _classify(r)[0] == "experimental" else 1)
    return hits[0]


def _value_of(row):
    for k in ("value", "propValue", "resultValue", "medianValue", "meanValue"):
        if k in row:
            return row[k]
    return None


def _unit_of(row):
    for k in ("unit", "units", "propUnit"):
        if row.get(k):
            return str(row[k])
    return ""


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
    """Physicochemical properties for a DTXSID -> {name: Prop} (+ 'K_AW' derived).

    Both the experimental and the predicted endpoints are read (the combined one
    first); whichever answers, the rows are matched by NAME, so a vocabulary change
    degrades to a missing property instead of a wrong one.
    """
    get = _get_fn or _get
    url = base_url(base)
    k = api_key(key)
    if not k or not dtxsid:
        return {}, ("no CTX API key" if not k else "no DTXSID")

    rows, errs = [], []
    for path in (f"/chemical/property/search/by-dtxsid/{dtxsid}",
                 f"/chemical/property/experimental/search/by-dtxsid/{dtxsid}",
                 f"/chemical/property/predicted/search/by-dtxsid/{dtxsid}"):
        payload, err = get(path, k, url)
        if err:
            errs.append(err)
        rows.extend(_rows(payload))                     # all three; _pick prefers experimental

    props = {}
    for name, patterns in _PROPERTY_PATTERNS.items():
        row = _pick(rows, patterns)
        if row is None:
            continue
        src, origin = _classify(row)
        props[name] = Prop(value=_num(_value_of(row)), source=src, origin=origin,
                           raw_value=_num(_value_of(row)), raw_unit=_unit_of(row))

    h = props.get("henry")
    if h is not None:
        kaw = henry_to_kaw(h.raw_value, h.raw_unit)
        if kaw is not None:
            props["K_AW"] = Prop(value=kaw, source=h.source, origin=h.origin,
                                 raw_value=h.raw_value, raw_unit=h.raw_unit or "atm-m3/mol")
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
        payload, err = _get(f"/chemical/search/equal/{q}", k)
        print(f"--- search/equal/{q}\n{err or json.dumps(payload, indent=2)[:4000]}")
        rows = _rows(payload)
        if rows:
            sid = rows[0].get("dtxsid")
            for path in (f"/chemical/detail/search/by-dtxsid/{sid}",
                         f"/chemical/property/search/by-dtxsid/{sid}"):
                p, e = _get(path, k)
                print(f"\n--- {path}\n{e or json.dumps(p, indent=2)[:6000]}")
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
    if r.note:
        print(f"  note: {r.note}")
    return 0


if __name__ == "__main__":                              # pragma: no cover
    sys.exit(_cli(sys.argv[1:]))
