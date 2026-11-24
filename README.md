# janaf-thermo

A small Python library with JANAF thermochemical tables for 15 species that come up in
iron ore processing, combustion and carbonate calcination. It gives Cp, H, S and G in SI
units at any temperature, and dH, dG and the equilibrium constant K of any reaction you
write with these species.

I wrote it as the thermodynamics layer of my own model of iron ore pellet induration.
I am publishing it on its own because the data took time to type in and check, and
someone else may find it useful.

## What is inside

| Key | Species | Phase | T range (K) |
|---|---|---|---|
| Fe2O3 | hematite | cr | 298 to 2000 |
| Fe3O4 | magnetite | cr | 298 to 2000 |
| FeO | wustite, Fe0.947O | cr | 298 to 1600 |
| SiO2 | quartz | cr | 298 to 2000 |
| MgO | periclase | cr | 298 to 1600 |
| MgCO3 | magnesite | cr | 298 to 1000 |
| Al2O3 | corundum | cr | 298 to 2000 |
| C | graphite | cr | 298 to 2000 |
| O2, N2, CO, CO2, H2O, CH4 | gases | g | 298 to 2000 |
| SO2 | gas | g | 298 to 1500 |

All tables come from the JANAF Thermochemical Tables, 2nd edition (1971). The page of
each table is stored next to the data. See [DATA_SOURCES.md](DATA_SOURCES.md).

## Install

```
pip install git+https://github.com/MomoHemasian/janaf-thermo.git
```

or clone the repository and run `pip install -e .[test]`. It needs Python 3.9 or newer,
NumPy and SciPy.

## Use

```python
from janaf_thermo import species, reaction

fe3o4 = species("Fe3O4")
fe3o4.cp(1000.0)      # J/(mol K)
fe3o4.h(1000.0)       # J/mol, formation basis: h(298.15) = standard enthalpy of formation
fe3o4.s(1000.0)       # J/(mol K)
fe3o4.g(1000.0)       # J/mol, G = H - T S
fe3o4.source          # 'JANAF Thermochemical Tables, 2nd ed., NSRDS-NBS 37 (1971), p. 791'

# any reaction as a dict: products positive, reactants negative
r = reaction({"C": -1, "CO2": -1, "CO": 2}, [900.0, 1000.0, 1100.0])
r["dH"], r["dG"], r["dS"], r["K"]   # J, J, J/K, K with gases in atm
r["extrapolated"]                   # True if a species was used outside its table
```

All functions take a number or a NumPy array of temperatures in K. A few ready reactions
are in `janaf_thermo.EXAMPLES` (magnetite oxidation, hematite dissociation, carbon and CO
combustion, Boudouard, methane combustion, magnesite decomposition).
`examples/quick_start.py` shows a short run.

## How it works

The tables are stored as printed (cal and kcal, 100 K steps) and converted to J when
loaded. H - H298, S and Cp are each interpolated with a shape-preserving PCHIP spline, so
the functions return the table values exactly at the table temperatures. Outside the
table the last Cp is held constant, and H and S are extended with that Cp, so they stay
consistent with each other.

## Tests

```
pytest
```

There are 120 test cases. They check three things:

1. The typed tables against themselves. Cp, H - H298 and S are linked (dH = Cp dT and
   dS = Cp/T dT), so a wrong digit in one column shows up as a mismatch with the others.
   Every 100 K step agrees within 1 %, except the three steps that contain a real
   transition (hematite near 950 K, magnetite near 848 K, quartz near 848 K), where the
   Cp peak is narrower than the grid.
2. Standard values at 298.15 K against the current NIST-JANAF tables (4th edition).
3. Reaction numbers that are well known: the Boudouard equilibrium near 700 °C, hematite
   losing oxygen in air near 1385 °C, about 50 MJ/kg for the lower heating value of methane
   and about -507 kJ/kg Fe3O4 for magnetite oxidation.

## Limits

- 15 species only, and only the 2nd edition values. Newer data (JANAF 4th edition, 1998)
  differs a little for some species. For most gases the difference is about 0.1 J/(mol K)
  in S at 298 K.
- H and S are interpolated separately. Between table points dG/dT = -S holds only to
  about 3 % in the worst case (MgO). At the table points everything is exact. A version
  where Cp, H and S come from one consistent fit is on my list.
- Near the three transitions above, `cp(T)` is a smooth curve and does not show the peak.
  `h(T)` still follows the table, so heat balances are correct over those steps.
- No phase changes beyond what the tables contain (for example wustite above 1600 K or
  liquid water are not included).

## Licence

MIT, see [LICENSE](LICENSE). The data source and its status are described in
[DATA_SOURCES.md](DATA_SOURCES.md).

## Author

Mohamad Hosein (Momo) Hemasian Etefagh, chemical engineer, Padova.
