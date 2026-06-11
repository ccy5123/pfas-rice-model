# Canonical alias: the original soil_paddy.py / calibration.py import the name
# `pfas_rice_plant_module`. Point it at the basis-A 4-pool (+surf) branch.
# (H7: the legacy naive `pfas_rice_plant_module.py` is DEPRECATED — do not use.)
from pfas_rice_plant_module_4pool_surf import (
    PlantInputs, Environment, Compound, Compartment, RiceUptakeModel,
    binding_factors, root_uptake, _ghk_factor, _logistic,
    ROOT, STEM, LEAF, FRUIT,
)
