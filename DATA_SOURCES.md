# Data sources

## JANAF Thermochemical Tables, 2nd edition

D. R. Stull and H. Prophet (project directors), *JANAF Thermochemical Tables*, 2nd ed.,
NSRDS-NBS 37, U.S. National Bureau of Standards, Washington DC, 1971.

All 15 tables in `src/janaf_thermo/data/janaf_2nd_ed.json` come from this publication.
I typed them in by hand from the printed pages. The page number of each table is stored
in the `page` field.

What was copied for each species: molar mass, standard enthalpy of formation at 298 K,
and the columns T, Cp, S and H - H298 at 298 K and every 100 K up to the end of the
table (or 2000 K). Units are as printed: Cp and S in cal/(mol K), H - H298 and dHf298 in
kcal/mol (thermochemical calorie, 4.184 J).

NSRDS-NBS 37 is a publication of the U.S. National Bureau of Standards (now NIST).
Works of the U.S. federal government are not under copyright in the United States.
The numbers themselves are measured and evaluated physical data.

The tests check the typed values for internal consistency and against the current
NIST-JANAF tables. If you find a value that does not match the printed page, please
open an issue with the species, temperature and page.

## Reference values used only in the tests

M. W. Chase Jr., *NIST-JANAF Thermochemical Tables*, 4th ed., J. Phys. Chem. Ref. Data,
Monograph 9, 1998. Six standard values at 298.15 K (dHf and S of CO2, CO, H2O, CH4, O2,
N2) are written in `tests/test_data.py` as a check. No tables from the 4th edition are
included.

## Not included

No data from other compilations (for example Barin, *Thermochemical Data of Pure
Substances*) is in this repository.
