"""Checks on the transcribed tables themselves.

The tables were typed in by hand from the printed JANAF pages, so the main risk is a
typing error. Cp, H - H298 and S are not independent (dH = Cp dT, dS = Cp/T dT), so a
wrong digit in one column shows up as a mismatch with the other two.
"""
import numpy as np
import pytest
from scipy.integrate import quad

from janaf_thermo import CAL, raw_table, species, species_keys

KEYS = species_keys()

# Intervals (upper T of the 100 K step) that hold a real lambda or alpha-beta transition.
# The Cp peak is narrower than the 100 K grid, so H cannot be rebuilt from Cp there.
TRANSITIONS = {
    "Fe2O3": {1000},   # hematite magnetic transition near 950 K
    "Fe3O4": {900},    # magnetite Curie point near 848 K
    "SiO2": {900},     # alpha-beta quartz near 848 K
}


def test_fifteen_species():
    assert len(KEYS) == 15


@pytest.mark.parametrize("key", KEYS)
def test_table_shape(key):
    d = raw_table(key)
    n = len(d["T"])
    assert len(d["Cp"]) == len(d["S"]) == len(d["HmH298"]) == n
    assert d["T"][0] == 298 and d["HmH298"][0] == 0
    assert d["T"][1] == 400 and np.all(np.diff(d["T"][1:]) == 100)  # 298, 400, 500, ...
    assert d["source"].startswith("JANAF") and isinstance(d["page"], int)


@pytest.mark.parametrize("key", KEYS)
def test_columns_increase(key):
    d = raw_table(key)
    assert np.all(np.diff(d["HmH298"]) > 0)
    assert np.all(np.diff(d["S"]) > 0)
    assert np.all(np.asarray(d["Cp"]) > 0)


@pytest.mark.parametrize("key", KEYS)
def test_entropy_agrees_with_enthalpy(key):
    """dS from the tabulated dH (mean Cp over the step) matches tabulated dS within 1 %."""
    d = raw_table(key)
    T = np.asarray(d["T"], float)
    H = np.asarray(d["HmH298"]) * 1000
    S = np.asarray(d["S"])
    dS_pred = np.diff(H) / np.diff(T) * np.log(T[1:] / T[:-1])
    assert np.max(np.abs(dS_pred - np.diff(S)) / np.diff(S)) < 0.01


@pytest.mark.parametrize("key", KEYS)
def test_enthalpy_agrees_with_cp(key):
    """The integral of the Cp column over each step matches dH within 1 %, except at transitions."""
    d = raw_table(key)
    sp = species(key)
    T = d["T"]
    H = np.asarray(d["HmH298"]) * 1000 * CAL
    for a, b, dH in zip(T[:-1], T[1:], np.diff(H)):
        if b in TRANSITIONS.get(key, set()):
            continue
        integral = quad(lambda t: float(sp.cp(t)), a, b, limit=200)[0]
        assert abs(integral - dH) / dH < 0.01, (key, b)


# Standard values at 298.15 K from the current NIST-JANAF tables (4th ed., 1998).
# The 2nd edition differs a little; these checks only catch gross typing errors.
REF_4TH = {           # dHf298 kJ/mol, S298 J/(mol K)
    "CO2": (-393.522, 213.785),
    "CO": (-110.527, 197.653),
    "H2O": (-241.826, 188.834),
    "CH4": (-74.873, 186.251),
    "O2": (0.0, 205.147),
    "N2": (0.0, 191.609),
}


@pytest.mark.parametrize("key", REF_4TH)
def test_standard_values_close_to_4th_edition(key):
    sp = species(key)
    dHf, S = REF_4TH[key]
    assert sp.dHf298 / 1000 == pytest.approx(dHf, abs=0.2)
    assert float(sp.s(298.15)) == pytest.approx(S, abs=0.2)


@pytest.mark.parametrize("key,M", [("CO2", 44.0095), ("H2O", 18.0153), ("Fe2O3", 159.69), ("Fe3O4", 231.53)])
def test_molar_mass(key, M):
    assert species(key).M * 1000 == pytest.approx(M, rel=1e-3)
