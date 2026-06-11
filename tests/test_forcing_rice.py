"""Tests for the rice transpiration forcing Q_TP(t) (task 2, part 1).

Locks in that the measured-data-anchored transpiration curve (Kumari 2022 Kc/ET0
+ Nay Htoon 2018 T/ET partition + FAO-56 stages) has the right magnitude and the
right transpiration/ET fraction -- the drivers that set the nstem f_xy scale."""
import numpy as np
import forcing_rice as fr


def test_seasonal_transpiration_fraction_matches_nayhtoon():
    # flooded paddy: seasonal T/ET ~ 0.42 (Nay Htoon 2018; evaporation ~58%)
    assert fr.seasonal_T_over_ET(120.0) == __import__("pytest").approx(0.42, abs=0.03)


def test_peak_QTP_is_realistic_and_mid_late_season():
    t = np.linspace(0, 120.0, 481)
    Q = fr.Q_TP(t, 120.0)
    # realistic per-hill transpiration peak (mm/day * ~0.04 m^2): ~0.05-0.20 L/d
    assert 0.05 < Q.max() < 0.20
    # peak at canopy closure (mid-season), not at transplant or harvest
    assert 0.3 < t[Q.argmax()] / 120.0 < 0.85
    # rises from a low base (bare/flooded -> evaporation-dominated) and ends lower
    assert Q[0] < 0.2 * Q.max()
    assert Q[-1] < Q.max()


def test_kc_and_ft_within_bounds():
    t = np.linspace(0, 120.0, 200)
    kc = fr.kc_curve(t, 120.0)
    ft = fr.ft_curve(t, 120.0)
    assert np.all((kc >= 0.9) & (kc <= 1.3))           # Kumari actual Kc range
    assert np.all((ft >= fr.F_T_INI - 1e-9) & (ft <= fr.F_T_FULL + 1e-9))
    assert ft[0] < ft[-1]                               # transpiration fraction rises with canopy


if __name__ == "__main__":
    import sys, pytest
    sys.exit(pytest.main([__file__, "-v"]))
