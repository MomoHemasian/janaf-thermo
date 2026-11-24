"""Reaction values checked against numbers that are well known in the literature."""
import numpy as np
import pytest

from janaf_thermo import EXAMPLES, R_GAS, reaction, species


@pytest.mark.parametrize("name", EXAMPLES)
def test_examples_balance_elements(name):
    """Every example reaction conserves Fe, C, O, H, Mg (counts per formula unit)."""
    atoms = {
        "Fe2O3": {"Fe": 2, "O": 3}, "Fe3O4": {"Fe": 3, "O": 4}, "O2": {"O": 2}, "C": {"C": 1},
        "CO2": {"C": 1, "O": 2}, "CO": {"C": 1, "O": 1}, "CH4": {"C": 1, "H": 4},
        "H2O": {"H": 2, "O": 1}, "MgCO3": {"Mg": 1, "C": 1, "O": 3}, "MgO": {"Mg": 1, "O": 1},
    }
    total = {}
    for sp, c in EXAMPLES[name].items():
        for el, n in atoms[sp].items():
            total[el] = total.get(el, 0) + c * n
    assert all(abs(v) < 1e-12 for v in total.values()), total


def test_methane_lhv():
    """Lower heating value of methane, about 50.0 MJ/kg."""
    dH = reaction("methane_combustion", 298.15)["dH"]
    assert -dH / species("CH4").M / 1e6 == pytest.approx(50.0, abs=0.1)


def test_boudouard_k_equals_one_near_700_C():
    T = np.arange(900.0, 1100.0, 0.5)
    lnK = np.log(reaction("boudouard", T)["K"])
    T_eq = T[np.argmin(np.abs(lnK))]
    assert T_eq - 273.15 == pytest.approx(700, abs=10)


def test_hematite_dissociates_in_air_near_1385_C():
    """6 Fe2O3 -> 4 Fe3O4 + O2: K = pO2 = 0.21 atm at about 1385 C."""
    T = np.arange(1550.0, 1750.0, 0.5)
    K = reaction("hematite_dissociation", T)["K"]
    T_air = T[np.argmin(np.abs(np.log(K / 0.21)))]
    assert T_air - 273.15 == pytest.approx(1385, abs=15)


def test_magnetite_oxidation_heat():
    """Magnetite oxidation releases about 0.5 MJ per kg Fe3O4."""
    dH = reaction("magnetite_oxidation", 298.15)["dH"]
    assert dH / species("Fe3O4").M / 1e3 == pytest.approx(-507, abs=10)


def test_consistency_of_outputs():
    T = np.array([600.0, 900.0, 1200.0])
    r = reaction({"CO": -1, "O2": -0.5, "CO2": 1}, T)
    assert np.allclose(r["dG"], r["dH"] - T * r["dS"])
    assert np.allclose(r["K"], np.exp(-r["dG"] / (R_GAS * T)))
    assert r["extrapolated"] is False


def test_reaction_flags_extrapolation():
    assert reaction("magnesite_decomposition", 1200.0)["extrapolated"] is True
