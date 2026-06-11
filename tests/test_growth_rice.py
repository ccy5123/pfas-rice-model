"""Tests for the ORYZA-based biomass forcing M_s(t) (task 2, part 2)."""
import numpy as np
import growth_rice as gr


def test_final_biomass_magnitude_and_partition():
    t = np.linspace(0, 120.0, 481)
    b = gr.organ_biomass(t, 120.0)
    shoot = b["stem"][-1] + b["leaf"][-1] + b["grain"][-1]
    # renormalized to ~1740 g/m^2 shoot (Bouman 2006 IR72) -> kg/hill at 25 hills/m^2
    assert shoot == __import__("pytest").approx(1740.0 * gr.G_M2_TO_KG_HILL, rel=1e-3)
    hi = b["grain"][-1] / shoot
    assert 0.40 < hi < 0.60               # harvest index, modern rice
    assert b["root"][-1] > 0              # roots accrue during vegetative growth


def test_masses_monotone_and_grain_after_flowering():
    t = np.linspace(0, 120.0, 481)
    b = gr.organ_biomass(t, 120.0)
    for k in ("root", "stem", "leaf", "grain"):
        assert np.all(np.diff(b[k]) >= -1e-12)        # cumulative biomass non-decreasing
    # grain is ~0 well before flowering (DVS=1 at day 65), present at harvest
    assert b["grain"][np.argmin(abs(t - 30))] < 0.02 * b["grain"][-1]
    assert b["grain"][-1] > b["leaf"][-1] * 0.5


def test_M_for_nstem_shape():
    t = np.linspace(0, 120.0, 200)
    M = gr.M_for_nstem(t, N=4, season=120.0)
    assert M.shape == (200, 6)            # root + 4 stem segments + grain
    assert np.all(M > 0)                  # floored, strictly positive for the ODE


if __name__ == "__main__":
    import sys, pytest
    sys.exit(pytest.main([__file__, "-v"]))
