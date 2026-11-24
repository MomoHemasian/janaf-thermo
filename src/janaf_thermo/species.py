"""Properties of one species from the JANAF tables, in SI units.

The tables are stored as printed in JANAF (cal, kcal, 100 K grid). Here they are
converted to J and interpolated with a shape-preserving PCHIP spline on H - H298,
S and Cp. Outside the tabulated range the end value of Cp is held constant, so
enthalpy and entropy stay consistent with each other.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from importlib import resources
from typing import Callable

import numpy as np
from scipy.interpolate import PchipInterpolator

CAL = 4.184  # J per thermochemical calorie
DATA_FILE = "janaf_2nd_ed.json"


@lru_cache(maxsize=1)
def _database() -> dict:
    text = resources.files("janaf_thermo").joinpath("data", DATA_FILE).read_text(encoding="utf-8")
    return json.loads(text)


def species_keys() -> list[str]:
    """Keys of all species in the data set, for example 'Fe3O4' or 'CO2'."""
    return list(_database().keys())


def raw_table(key: str) -> dict:
    """The table of one species as transcribed (cal and kcal units), with its source page."""
    db = _database()
    if key not in db:
        raise KeyError(f'Unknown species "{key}". Known: {", ".join(db)}')
    return json.loads(json.dumps(db[key]))  # a copy, so the cache cannot be changed


@dataclass(frozen=True)
class Species:
    key: str
    name: str
    phase: str
    M: float          # molar mass, kg/mol
    dHf298: float     # standard enthalpy of formation at 298.15 K, J/mol
    Tmin: float       # K
    Tmax: float       # K
    source: str
    cp: Callable      # Cp(T), J/(mol K)
    h: Callable       # H(T) on the formation basis, H(298) = dHf298, J/mol
    s: Callable       # S(T), J/(mol K)
    g: Callable       # G(T) = H - T S, J/mol

    def in_range(self, T) -> bool:
        """True if every T lies inside the tabulated range (no extrapolation)."""
        T = np.asarray(T, float)
        return bool(np.all((T >= self.Tmin) & (T <= self.Tmax)))


@lru_cache(maxsize=None)
def species(key: str) -> Species:
    """Build the property functions of one species. All inputs in K, outputs in SI."""
    d = raw_table(key)
    T = np.asarray(d["T"], float)
    H = np.asarray(d["HmH298"], float) * 1000 * CAL
    S = np.asarray(d["S"], float) * CAL
    Cp = np.asarray(d["Cp"], float) * CAL
    Tmin, Tmax = float(T[0]), float(T[-1])
    ppH, ppS, ppC = PchipInterpolator(T, H), PchipInterpolator(T, S), PchipInterpolator(T, Cp)
    dHf = float(d["dHf298"]) * 1000 * CAL

    def cp(t):
        return ppC(np.clip(np.asarray(t, float), Tmin, Tmax))

    def dh(t):
        t = np.asarray(t, float)
        return ppH(np.clip(t, Tmin, Tmax)) + Cp[-1] * np.maximum(t - Tmax, 0) + Cp[0] * np.minimum(t - Tmin, 0)

    def s(t):
        t = np.asarray(t, float)
        return (ppS(np.clip(t, Tmin, Tmax))
                + Cp[-1] * np.log(np.maximum(t, Tmax) / Tmax)
                + Cp[0] * np.log(np.minimum(t, Tmin) / Tmin))

    def h(t):
        return dHf + dh(t)

    def g(t):
        return h(t) - np.asarray(t, float) * s(t)

    return Species(key=key, name=d["name"], phase=d["phase"], M=float(d["M"]) / 1000,
                   dHf298=dHf, Tmin=Tmin, Tmax=Tmax, source=f'{d["source"]}, p. {d["page"]}',
                   cp=cp, h=h, s=s, g=g)
