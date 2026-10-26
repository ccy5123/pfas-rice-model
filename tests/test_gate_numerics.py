"""Guards for the numerical floor and the grain formation gate (validation/numerics_convergence.py).

The floor is admissible only because the results do not depend on it; these tests pin that the
convergence script's 'shipped' row is the code, that there is ONE floor, and that the gate is the
step it is documented to be."""
import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "src"))
sys.path.insert(0, os.path.join(_ROOT, "validation"))


def test_convergence_script_shipped_row_is_the_code():
    import inspect
    import numerics_convergence as N
    import organ_mass as om
    import pfas_rice_plant_module_4pool_surf as P4
    assert N.SHIPPED["floor"] == om.MASS_FLOOR == P4.MASS_FLOOR
    assert inspect.signature(om.organ_matrix).parameters["floor"].default == om.MASS_FLOOR
    assert inspect.signature(P4.grain_formed).parameters["floor"].default == om.MASS_FLOOR


def test_gate_is_a_step_off_the_floor():
    import organ_mass as om
    import pfas_rice_plant_module_nstem_leaf as NL
    f = om.MASS_FLOOR
    gm = np.r_[np.full(5, f), 1.0001 * f, 2 * f, 1e-6, 1e-4, 1e-3, 1e-2]
    np.testing.assert_array_equal(NL.grain_formed(gm), (gm > f).astype(float))
    np.testing.assert_array_equal(NL.grain_formed(np.full(4, 0.05)), np.ones(4))   # never on the floor
    t = np.linspace(0.0, 10.0, 11)
    M = np.column_stack([np.ones(11)] * 6 + [gm])
    inp = NL.PlantInputsNL(t=t, Cwo=np.ones(11), Qtp=np.ones(11), M=M)
    assert inp.grain_gate_(4.0) == 0.0 and inp.grain_gate_(5.0) == 1.0


def test_crop_drivers_carry_no_floor_of_their_own():
    """The crop model's grain before panicle formation is 0; only organ_matrix floors it."""
    import model_api as api
    run = api.site_run()
    assert float(np.min(run["grain"])) == 0.0
    import organ_mass as om
    d = api.scenario_drivers()
    assert float(np.min(np.asarray(d["M"])[:, 3])) == om.MASS_FLOOR
