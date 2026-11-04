"""FastAPI web backend for the PFAS–rice model (the React front end's API).

The Streamlit app (``app.py``) re-runs its whole script on every widget touch;
this backend instead exposes the SAME ``model_api`` computations as JSON
endpoints so a JavaScript front end can keep all interaction (sliders, tabs,
the day scrubber) in the browser and only ask the server when a model input
changes. Run with ``uvicorn server.app:app`` from the repo root.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src")
for _p in (SRC, ROOT):
    if _p not in sys.path:
        sys.path.insert(0, _p)
