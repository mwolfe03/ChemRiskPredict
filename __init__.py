# ChemRiskPredict/__init__.py

# Import functions and classes from your modules
from .Data_Collection_from_Pubchem import *
from .Hazard_Prediction_Model import *
from .Hazard_and_Subgroup_Column_Initialization import *
from .ring_handling_w_rdkit import *
from .sub_groups_from_smiles import *

# For accessing the CSV files
import os

def get_data_file(filename):
    """
    Get the full path to one of the data files in the package.
    """
    package_directory = os.path.dirname(__file__)
    return os.path.join(package_directory, filename)


__all__ = [
    'determine_potential_hazards_from_dataframe', 'hazard_model_testing', 'create_dataframe_from_ids',
    'update_existing_dataframe_from_dataframe', 'update_existing_ids_dataframe_from_cids'
]
