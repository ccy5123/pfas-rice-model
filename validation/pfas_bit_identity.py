"""PFAS bit-identity check for the neutral / weak-electrolyte extension.

The extension's top invariant is that the PFAS path does not move a single bit.
Every field the neutral branch adds defaults to the value at which its term vanishes
identically -- not approximately -- so the correct test is exact ``==``, NOT
``np.allclose``: a tolerance would hide precisely the kind of drift this is guarding
against (a reordered sum, an added-then-subtracted zero, a changed default).

Usage::

    # extract the checkout to compare against, anywhere
    mkdir -p /tmp/base && git archive <commit> | tar -x -C /tmp/base
    python validation/pfas_bit_identity.py /tmp/base

    # or compare two explicit trees
    python validation/pfas_bit_identity.py /tmp/base /path/to/new

Exits 0 when identical, 1 on any difference (so it can gate a commit).  It runs the
grid in a SUBPROCESS per checkout, because both trees define the same module names
and cannot be imported into one interpreter.

This script was reconstructed each session before it was committed; it lives here so
the method is reproducible rather than remembered.  See
docs/HANDOFF_neutral_extension.md section 2.
"""
import json
import os
import subprocess
import sys

# The grid: 4 congeners spanning the chain-length range (short / mid / sulfonate /
# long) x every parameter source and both opt-in mechanisms, so a change confined to
# one branch still shows up somewhere.
CONGENERS = ("PFBA", "PFOA", "PFOS", "PFDoDA")
F_XY_SOURCES = ("recommended", "W2fit", "oryza")

_RUNNER = r'''
import json, os, sys
sys.path.insert(0, os.path.join(sys.argv[1], "src"))
import model_api as api
CONGENERS = %(cong)r
F_XY_SOURCES = %(srcs)r
out = {}
for c in CONGENERS:
    for src in F_XY_SOURCES:
        r = api.simulate(c, f_xy_source=src)
        out["%%s|%%s" %% (c, src)] = (
            [r["baf_final"][t] for t in ("root", "stem", "leaf", "grain")]
            + [r["straw_baf"], r["N"], r["eN"]]
            + [r["B_k"][t] for t in ("root", "stem", "leaf", "grain")])
    for tag, kw in (("lipid", dict(lipid_loading=True)),
                    ("flooded", dict(cwo_profile="flooded"))):
        r = api.simulate(c, **kw)
        out["%%s|%%s" %% (c, tag)] = (
            [r["baf_final"][t] for t in ("root", "stem", "leaf", "grain")] + [r["straw_baf"]])
print(json.dumps(out))
''' % {"cong": CONGENERS, "srcs": F_XY_SOURCES}


def _run(root):
    """Run the grid inside `root`'s own checkout and return {scenario: [floats]}."""
    p = subprocess.run([sys.executable, "-c", _RUNNER, root],
                       capture_output=True, text=True, cwd=root)
    if p.returncode != 0:
        sys.exit(f"runner failed in {root}:\n{p.stderr[-3000:]}")
    return json.loads(p.stdout)


def compare(base, new, verbose=True):
    """Return the number of differing floats between two checkouts (0 == identical)."""
    a, b = _run(base), _run(new)
    if set(a) != set(b):
        sys.exit(f"scenario sets differ: {set(a) ^ set(b)}")
    ndiff = nfloat = 0
    for k in sorted(a):
        for i, (x, y) in enumerate(zip(a[k], b[k])):
            nfloat += 1
            if x != y:                      # EXACT equality, deliberately not allclose
                ndiff += 1
                if verbose:
                    print(f"  DIFF {k}[{i}]  base={x!r}  new={y!r}")
    if verbose:
        print(f"scenarios={len(a)}  floats={nfloat}  differing={ndiff}")
        print("PFAS BIT-IDENTICAL ✓" if ndiff == 0 else "PFAS CHANGED ✗")
    return ndiff


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    base = sys.argv[1]
    new = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))
    sys.exit(1 if compare(base, new) else 0)
