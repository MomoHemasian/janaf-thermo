"""Reaction thermodynamics from the species tables.

A reaction is a dict species -> stoichiometric coefficient (products +, reactants -).
Gases are referred to 1 atm (the JANAF standard state), condensed phases to activity 1.
"""
from __future__ import annotations

import numpy as np

from .species import species

R_GAS = 8.314462618  # J/(mol K)

# A few textbook reactions that use only species in this data set.
EXAMPLES = {
    "magnetite_oxidation": {"Fe3O4": -1, "O2": -0.25, "Fe2O3": 1.5},       # 4 Fe3O4 + O2 -> 6 Fe2O3
    "hematite_dissociation": {"Fe2O3": -6, "Fe3O4": 4, "O2": 1},           # 6 Fe2O3 -> 4 Fe3O4 + O2
    "carbon_combustion": {"C": -1, "O2": -1, "CO2": 1},                    # C + O2 -> CO2
    "boudouard": {"C": -1, "CO2": -1, "CO": 2},                            # C + CO2 -> 2 CO
    "co_oxidation": {"CO": -1, "O2": -0.5, "CO2": 1},                      # CO + 1/2 O2 -> CO2
    "methane_combustion": {"CH4": -1, "O2": -2, "CO2": 1, "H2O": 2},       # CH4 + 2 O2 -> CO2 + 2 H2O(g)
    "magnesite_decomposition": {"MgCO3": -1, "MgO": 1, "CO2": 1},          # MgCO3 -> MgO + CO2
}


def reaction(nu, T) -> dict:
    """dH, dS, dG (per mol of reaction as written) and K at temperature(s) T in K.

    nu: dict species -> coefficient, or the name of one of EXAMPLES.
    Returns a dict with 'dH' (J), 'dS' (J/K), 'dG' (J), 'K' (-) and 'extrapolated'
    (True if any species was used outside its tabulated range).
    """
    if isinstance(nu, str):
        nu = EXAMPLES[nu]
    T = np.asarray(T, float)
    dH = np.zeros_like(T)
    dG = np.zeros_like(T)
    extrapolated = False
    for key, c in nu.items():
        sp = species(key)
        dH = dH + c * sp.h(T)
        dG = dG + c * sp.g(T)
        extrapolated = extrapolated or not sp.in_range(T)
    return {"dH": dH, "dS": (dH - dG) / T, "dG": dG, "K": np.exp(-dG / (R_GAS * T)),
            "extrapolated": extrapolated}
