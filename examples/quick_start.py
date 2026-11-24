"""Short tour: species properties and two reaction equilibria."""
import numpy as np

from janaf_thermo import reaction, species

fe3o4 = species("Fe3O4")
print(fe3o4.name, fe3o4.source)
for T in (298.15, 800.0, 1200.0):
    print(f"  T = {T:7.2f} K   Cp = {fe3o4.cp(T):7.2f} J/(mol K)   H = {fe3o4.h(T) / 1000:9.2f} kJ/mol")

# Temperature where hematite starts to lose oxygen in air (pO2 = 0.21 atm)
T = np.arange(1500.0, 1800.0, 0.5)
K = reaction("hematite_dissociation", T)["K"]
print(f"Fe2O3 -> Fe3O4 in air at about {T[np.argmin(abs(np.log(K / 0.21)))] - 273.15:.0f} C")

# Heat of magnetite oxidation per kg Fe3O4
dH = reaction("magnetite_oxidation", 298.15)["dH"]
print(f"Magnetite oxidation: {dH / fe3o4.M / 1000:.0f} kJ per kg Fe3O4")
