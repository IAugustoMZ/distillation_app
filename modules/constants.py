import os

# project root path
PROJECT_ROOT = os.path.dirname(
    os.path.dirname(
        os.path.realpath(__file__)
    )
)
# path for the physical and chemical properties
# database
PHYSICAL_CHEMICAL_SOURCE = os.path.join(
    PROJECT_ROOT, 'data', 'physical_chemical_properties.xlsx'
)

# path for the file that stores the
# vle equilibrium data
EQUILIBRIUM_DATA = os.path.join(
    PROJECT_ROOT, 'data', 'equilibrium_curve_data.csv'
)

# path for the thermodynamic module
THERMODYNAMICS_PATH = os.path.join(
    PROJECT_ROOT, 'modules', 'thermodynamic'
)

# information about the feed
Cp = 159        # heat capacity [J/mol.K]
Hvap = 32.099   # enthalpy of vaporization [kJ/mol]
Tref = 298.15   # reference state temperature [K]