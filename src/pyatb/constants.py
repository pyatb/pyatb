"""
Convert
"""
# convert Rydberg Constant to Electron-volt
Ry_to_eV = 13.605698066

# convert Angstrom to Bohr Radius
Ang_to_Bohr = 1.8897259886

# convert Joule to Electron-volt
J_to_eV =  6.241506363094e18


"""
Constant
"""
# Boltzmann constant
k_B_SI = 1.380649e-23    # J/K
k_B_eV = k_B_SI * J_to_eV  # eV/K
k_B_Ry = k_B_SI * J_to_eV / Ry_to_eV # Ry/K

# elementary charge
elem_charge_SI = 1.60217662e-19  # C

# electron mass
electron_mass_SI = 9.10938188e-31 # kg

# Planck constant
h_SI = 6.62607004e-34  # m^2 kg / s
hbar_SI = 1.054571628e-34  # m^2 kg / s

# von Klitzing constant = h / e^2 
R_k = 25812.80745  # \omega