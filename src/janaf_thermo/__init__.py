"""janaf_thermo: JANAF (2nd ed., 1971) thermochemical tables for 15 species, in Python.

    from janaf_thermo import species, reaction
    sp = species("Fe3O4")
    sp.h(1000.0), sp.s(1000.0), sp.g(1000.0), sp.cp(1000.0)   # SI units
    reaction("boudouard", 1000.0)["K"]
"""
from .reactions import EXAMPLES, R_GAS, reaction
from .species import CAL, Species, raw_table, species, species_keys

__version__ = "0.1.0"
__all__ = ["species", "species_keys", "raw_table", "Species", "reaction", "EXAMPLES", "R_GAS", "CAL"]
