import os

PROJECT_ROOT = os.path.dirname(
    os.path.dirname(
        os.path.realpath(__file__)
    )
)
PHYSICAL_CHEMICAL_SOURCE = os.path.join(
    PROJECT_ROOT, 'data', 'physical_chemical_properties.xlsx'
)
EQUILIBRIUM_DATA = os.path.join(
    PROJECT_ROOT, 'data', 'equilibrium_curve_data.csv'
)