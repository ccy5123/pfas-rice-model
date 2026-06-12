# Deploying / demoing the dashboard on another computer

The showcase artifact is the Streamlit app (`app.py`) — the plant/soil accumulation
colormap + interactive time series. Three ways to show it elsewhere, simplest first.

## A. Local run (most reliable; works offline)
Any computer with Python 3.9+:
```bash
git clone https://github.com/ccy5123/pfas-rice-model.git
cd pfas-rice-model
pip install -r requirements.txt        # full app stack
streamlit run app.py                    # opens http://localhost:8501 in the browser
```
4 of the 5 exposure modes (model / HYDRUS-CSV / soil-inventory / biomonitoring) work out
of the box; the **live HYDRUS-1D** mode needs the compiled engine (gfortran) and `phydrus`
and is hidden automatically when absent. The **SMILES** structure input is an API feature
(`model_api.simulate_from_smiles`), not an app mode, so RDKit is not needed to run the app.

## B. Streamlit Community Cloud (free hosting; share by URL — no install on the demo machine)
1. Make sure the repo is on GitHub (it is) and `requirements.txt` has the app stack (it does).
2. Go to **[share.streamlit.io](https://share.streamlit.io)** and sign in with GitHub.
3. **New app → From existing repo**:
   - Repository: `ccy5123/pfas-rice-model`
   - Branch: `main`
   - Main file path: `app.py`
4. **Deploy**. First build takes ~2–3 min (installs `requirements.txt`). You get a public
   URL like `https://<name>.streamlit.app` — open it on any computer or phone.

Notes:
- Streamlit Cloud installs **`requirements.txt`** (that is why the app deps live there, not
  only in `requirements-app.txt`). It deploys as-is.
- The **live HYDRUS-1D** mode will not run on Cloud (it needs the compiled FORTRAN engine);
  the app detects this and shows the build steps instead — the other four modes demo fully.
- RDKit / phydrus are **not** installed on Cloud (not in `requirements.txt`); nothing in
  `app.py` imports them, so the app boots cleanly.
- To redeploy after changes: just `git push` to the deployed branch — Cloud auto-rebuilds.

## C. No install at all (non-interactive backup)
- Show `docs/OVERVIEW_KR.md` (model schematic + validation tables) and the PNGs in
  `validation/figures/` as slides/PDF.
- Or run the app locally once and **screen-record / screenshot** it into a slide deck —
  zero risk of a live-demo failure.

## Which to pick
| situation | use |
|---|---|
| you control the demo laptop, want it bulletproof / offline | **A. local run** |
| share a link, others open on their own devices | **B. Streamlit Cloud** |
| no internet / no install / safest backup | **C. figures or a recording** |
