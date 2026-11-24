"""Checks on the property functions built from the tables."""
import numpy as np
import pytest

from janaf_thermo import CAL, raw_table, species, species_keys

KEYS = species_keys()


@pytest.mark.parametrize("key", KEYS)
def test_interpolation_hits_table_points(key):
    d = raw_table(key)
    sp = species(key)
    T = np.asarray(d["T"], float)
    assert np.allclose(sp.h(T) - sp.dHf298, np.asarray(d["HmH298"]) * 1000 * CAL, atol=1e-6)
    assert np.allclose(sp.s(T), np.asarray(d["S"]) * CAL, atol=1e-9)
    assert np.allclose(sp.cp(T), np.asarray(d["Cp"]) * CAL, atol=1e-9)


@pytest.mark.parametrize("key", KEYS)
def test_gibbs_helmholtz(key):
    """dG/dT = -S, checked with a central difference.

    H and S are interpolated separately, so between table points this holds only
    approximately (measured worst case 3.1 % for MgO, 2.5 J/(mol K) for Fe3O4).
    The limit below is set just above that, to catch a regression, not to hide it.
    """
    sp = species(key)
    T = np.linspace(sp.Tmin + 10, sp.Tmax - 10, 400)
    e = 1e-3
    dGdT = (sp.g(T + e) - sp.g(T - e)) / (2 * e)
    err = np.abs(dGdT + sp.s(T))
    assert np.max(err) < 3.0
    assert np.max(err / sp.s(T)) < 0.035


@pytest.mark.parametrize("key", ["O2", "CO2", "Fe2O3"])
def test_extrapolation_is_constant_cp(key):
    """Above Tmax, H grows with the last Cp and S with Cp ln(T/Tmax), so they stay consistent."""
    sp = species(key)
    T1, T2 = sp.Tmax, sp.Tmax + 300
    cp_end = float(sp.cp(T1))
    assert float(sp.h(T2) - sp.h(T1)) == pytest.approx(cp_end * 300, rel=1e-12)
    assert float(sp.s(T2) - sp.s(T1)) == pytest.approx(cp_end * np.log(T2 / T1), rel=1e-12)
    assert not sp.in_range(T2) and sp.in_range(T1)


def test_vectorised_and_scalar_agree():
    sp = species("Fe3O4")
    T = np.array([500.0, 1000.0, 1500.0])
    assert np.allclose(sp.h(T), [float(sp.h(t)) for t in T])


def test_unknown_species_message():
    with pytest.raises(KeyError, match="Known"):
        species("unobtainium")


def test_raw_table_is_a_copy():
    d = raw_table("O2")
    d["S"][0] = -1
    assert raw_table("O2")["S"][0] != -1
