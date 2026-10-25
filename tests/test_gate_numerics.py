"""Guard for validation/gate_numerics_sensitivity.py: its 'shipped' row must be the code's
actual organ-mass floor and gate ramp, or the register's gate row describes another model."""
import inspect
import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "src"))
sys.path.insert(0, os.path.join(_ROOT, "validation"))


def test_shipped_variant_is_the_code():
    import gate_numerics_sensitivity as G
    import organ_mass as om
    import pfas_rice_plant_module_nstem_leaf as NL
    floor, ramp = next(iter(G.VARIANTS.values()))
    assert inspect.signature(om.organ_matrix).parameters["floor"].default == floor
    # the gate ramp: fully open at (1 + 0.5) x the floor
    t = np.linspace(0.0, 10.0, 11)
    gm = np.r_[np.full(5, floor), floor * np.array([ramp, 2.0, 4.0, 8.0, 16.0, 32.0])]
    M = np.column_stack([np.ones(11)] * 6 + [gm])
    inp = NL.PlantInputsNL(t=t, Cwo=np.ones(11), Qtp=np.ones(11), M=M)
    assert inp.grain_gate_(4.0) == 0.0
    assert inp.grain_gate_(5.0) == 1.0
    assert 0.0 < inp.grain_gate_(4.5) < 1.0
