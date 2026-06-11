# Example inputs for the visualization tool

Ready-to-load CSVs for `app.py` (`streamlit run app.py`). Each loads automatically when the
matching mode is selected and no file is uploaded. See `docs/visualization_tool.md`.

## `hydrus_drivers_example.csv` — HYDRUS / CSV drivers mode
A synthetic HYDRUS-1D/Phydrus-style hand-off (Method A, one-way coupling).

| column | meaning | HYDRUS-1D source |
|---|---|---|
| `t` | day after transplant | output times |
| `Cwo` | pore-water free anion [µg/L] | `Conc` at the root-zone node (`Obs_Node.out`/`solute1.out`) |
| `Qtp` | transpiration stream [L/day] | `vRoot` / `T_act` (`T_Level.out`) |
| `M_root,M_stem,M_leaf,M_grain` | organ fresh mass [kg] | a plant growth sub-model |

Here `Cwo(t)` was produced by inverting a 5 µg/kg-dry soil inventory through the Freundlich
paddy soil with early-season flooding + leaching (`model_api.pore_water_from_inventory`); the
forcings are the measured FAO-56 transpiration + ORYZA IR72 biomass. Replace with your own
HYDRUS run — `Qtp` and the `M_*` columns are optional (omit them to fall back to the measured
crop forcings).

## `biomonitoring_example.csv` — Biomonitoring mode
Measured tissue concentrations + the pore-water/soil-solution concentration; no soil model
needed (`BAF = conc / Cwo`).

| column | meaning |
|---|---|
| `tissue` | `root` / `straw` / `stem` / `leaf` / `grain` (free text; these map onto the plant map) |
| `conc` | measured tissue concentration [µg/kg] |
| `Cwo` | measured pore-water / soil-solution concentration [µg/L] (first row used) |

The bundled values are the Yamazaki 2023 PFOA brown-rice field BAFs at `Cwo = 1` (so the
concentration equals the BAF).
