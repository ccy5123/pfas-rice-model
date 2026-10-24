"""A small thread-safe LRU memo with in-flight de-duplication.

The web front end asks for one scenario's result from several places at once (the
headline, then whichever tab is open), so two requests for the SAME solve often
arrive together. The first caller computes; the others wait for it instead of
starting a second identical ODE solve. Keys are canonical strings (``key_of``),
numpy arrays are hashed rather than serialised.
"""
from __future__ import annotations

import hashlib
import json
import math
import threading
from collections import OrderedDict

import numpy as np


def _canon(o):
    if isinstance(o, np.ndarray):
        a = np.ascontiguousarray(o, dtype=float)
        return {"__nd__": hashlib.sha1(a.tobytes()).hexdigest(), "shape": list(a.shape)}
    if isinstance(o, dict):
        return {str(k): _canon(v) for k, v in sorted(o.items(), key=lambda kv: str(kv[0]))}
    if isinstance(o, (list, tuple)):
        return [_canon(v) for v in o]
    if isinstance(o, (np.floating, float)):
        f = float(o)
        return f if math.isfinite(f) else repr(f)
    if isinstance(o, np.integer):
        return int(o)
    if hasattr(o, "model_dump"):                       # pydantic models
        return _canon(o.model_dump())
    return o


def key_of(*parts, **kw):
    """Canonical string key for any mix of scalars, dicts, lists and arrays."""
    return json.dumps(_canon([list(parts), kw]), sort_keys=True, default=str)


class Memo:
    def __init__(self, maxsize=128):
        self.maxsize = int(maxsize)
        self._d = OrderedDict()
        self._lock = threading.Lock()
        self._inflight = {}
        self.hits = 0
        self.misses = 0

    def get(self, key, fn):
        with self._lock:
            if key in self._d:
                self._d.move_to_end(key)
                self.hits += 1
                return self._d[key]
            ev = self._inflight.get(key)
            owner = ev is None
            if owner:
                ev = threading.Event()
                self._inflight[key] = ev
        if not owner:
            ev.wait()
            with self._lock:
                if key in self._d:
                    self.hits += 1
                    return self._d[key]
            return fn()                                # the owner failed: compute ourselves
        try:
            val = fn()
            with self._lock:
                self.misses += 1
                self._d[key] = val
                while len(self._d) > self.maxsize:
                    self._d.popitem(last=False)
            return val
        finally:
            with self._lock:
                self._inflight.pop(key, None)
            ev.set()

    def peek(self, key):
        """The cached value for `key`, or None (never computes)."""
        with self._lock:
            return self._d.get(key)

    def clear(self):
        with self._lock:
            self._d.clear()
            self.hits = self.misses = 0

    def __len__(self):
        return len(self._d)
