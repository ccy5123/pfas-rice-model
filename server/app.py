"""FastAPI app: JSON endpoints over ``server/service.py`` + the built React front end.

    pip install -r requirements-web.txt
    (cd frontend && npm ci && npm run build)      # once; produces frontend/dist
    uvicorn server.app:app --port 8000             # from the repo root

Open http://localhost:8000. During front-end development run ``npm run dev`` in
``frontend/`` instead (Vite proxies ``/api`` to port 8000). Interactive API docs
are at ``/api/docs``.
"""
from __future__ import annotations

import os
import threading
from contextlib import asynccontextmanager
from urllib.parse import quote

from fastapi import FastAPI, Query, Request
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.responses import FileResponse, JSONResponse, Response, StreamingResponse

from . import ROOT
from . import service as S
from .schemas import (CwoPreviewRequest, DownloadRequest, FigureRequest, InverseRequest,
                      Scenario)

DIST = os.environ.get("PFAS_WEB_DIST", os.path.join(ROOT, "frontend", "dist"))

@asynccontextmanager
async def _lifespan(_app):
    # Pre-solve the Simple view's default scenarios in the background (best-effort).
    if os.environ.get("PFAS_WEB_WARMUP", "1") != "0":
        threading.Thread(target=S.warmup, daemon=True).start()
    yield


app = FastAPI(title="PFAS–Rice Uptake Model API", version=S.APP_VERSION, lifespan=_lifespan,
              docs_url="/api/docs", openapi_url="/api/openapi.json", redoc_url=None)
app.add_middleware(GZipMiddleware, minimum_size=1024)


def _json(text_or_obj, status=200):
    """Pre-serialised JSON text (or an object) -> response, without re-encoding."""
    body = text_or_obj if isinstance(text_or_obj, str) else S.dumps(text_or_obj)
    return Response(content=body, media_type="application/json", status_code=status)


@app.exception_handler(S.ApiError)
async def _api_error(_: Request, ex: S.ApiError):
    return JSONResponse({"detail": ex.detail}, status_code=ex.status)


# ---------------------------------------------------------------- API
@app.get("/api/health")
def health():
    return {"ok": True, "version": S.APP_VERSION}


@app.get("/api/meta")
def meta():
    return _json(S.meta())


@app.post("/api/run")
def run(scn: Scenario):
    return _json(S.payload(S.run(scn)))


@app.post("/api/figure")
def figure(req: FigureRequest):
    return _json(S.figure(req))


@app.post("/api/cwo-preview")
def cwo_preview(req: CwoPreviewRequest):
    return _json(S.cwo_preview(req))


@app.post("/api/inverse")
def inverse(req: InverseRequest):
    return StreamingResponse(S.inverse_stream(req), media_type="application/x-ndjson")


@app.post("/api/download/{kind}")
def download(kind: str, req: DownloadRequest):
    body, media, fname, fname_utf8 = S.download(kind, req)
    disp = f'attachment; filename="{fname}"'
    if fname_utf8:
        disp += f"; filename*=UTF-8''{quote(fname_utf8)}"
    return Response(content=body, media_type=media + "; charset=utf-8",
                    headers={"Content-Disposition": disp})


@app.get("/api/tables/defaults")
def tables_defaults(season: float = Query(120.0, ge=30, le=200),
                    biomass: str = Query("oryza", pattern="^(oryza|growth_rice)$"),
                    Cwo: float = Query(1.0, ge=0)):
    return _json(S.tables_defaults(season, biomass, Cwo))


@app.get("/api/k-leach")
def k_leach(congener: str | None = None, log_kow: float | None = None, koc: float | None = None):
    """The HYDRUS-calibrated flooded-shape default (per congener; from Koc for a neutral)."""
    from .schemas import CompoundSpec, NeutralSpec
    if log_kow is not None or koc is not None:
        comp = CompoundSpec(kind="neutral", neutral=NeutralSpec(
            log_kow=2.45 if log_kow is None else log_kow, Koc=koc))
    elif congener:
        comp = CompoundSpec(kind="congener", congener=congener)
    else:
        comp = CompoundSpec(kind="smiles")
    return {"k_leach": S.default_k_leach(comp)}


@app.get("/api/smiles/svg")
def smiles_svg(smiles: str = Query(..., max_length=2000)):
    return S.mol_svg(smiles)


@app.get("/api/ctx")
def ctx(q: str = Query(..., min_length=1, max_length=500)):
    return _json(S.ctx_lookup(q.strip()))


@app.get("/api/hydrus/status")
def hydrus_status():
    import model_api as api
    return {"available": bool(api.hydrus_available())}


@app.post("/api/hydrus/build")
def hydrus_build():
    import model_api as api
    ok, log = api.build_hydrus_engine()
    S.FIG.clear()
    S.RUN.clear()
    return {"ok": bool(ok), "log": list(log)}


@app.get("/api/cache")
def cache_stats():
    return {n: {"size": len(m), "hits": m.hits, "misses": m.misses}
            for n, m in (("sim", S.SIM), ("run", S.RUN), ("drivers", S.DRV),
                         ("figures", S.FIG), ("aux", S.AUX))}


# ---------------------------------------------------------------- the React app
@app.get("/{path:path}", include_in_schema=False)
def spa(path: str):
    if path.startswith("api/"):
        return JSONResponse({"detail": "not found"}, status_code=404)
    index = os.path.join(DIST, "index.html")
    if not os.path.isfile(index):
        return JSONResponse({"detail": "front end not built — run `npm ci && npm run build` "
                                       "in frontend/ (the API itself is live: /api/docs)"},
                            status_code=503)
    target = os.path.normpath(os.path.join(DIST, path))
    if path and target.startswith(os.path.normpath(DIST) + os.sep) and os.path.isfile(target):
        headers = ({"Cache-Control": "public, max-age=31536000, immutable"}
                   if "/assets/" in target.replace(os.sep, "/") else None)
        return FileResponse(target, headers=headers)
    return FileResponse(index, headers={"Cache-Control": "no-cache"})
