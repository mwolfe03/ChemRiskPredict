# ChemRiskPredict/__init__.py

# Import functions and classes from your modules
from ChemRiskPredict.dataframe_and_model_creation.Data_Collection_from_Pubchem import *
from ChemRiskPredict.dataframe_and_model_creation.Hazard_Prediction_Model import *
from ChemRiskPredict.dataframe_and_model_creation.Hazard_and_Subgroup_Column_Initialization import *
from ChemRiskPredict.smiles_to_groupings.ring_handling_w_rdkit import *
from ChemRiskPredict.smiles_to_groupings.sub_groups_from_smiles import *

# For accessing the CSV files
import os

def get_data_file(filename):
    """
    Get the full path to one of the data files in the package.
    """
    package_directory = os.path.dirname(__file__)
    return os.path.join(package_directory, filename)


__all__ = [
    'determine_potential_hazards_from_dataframe', 'hazard_model_testing',
    'update_existing_dataframe_from_dataframe', 'update_existing_ids_dataframe_from_cids'
]
