from rdkit import Chem
from rdkit.Chem import rdmolops



def ring_handling_rdkit(SMILES_string: str) -> list | bool:
    """
    Input: SMILES_string is a string of the SMILES representation
    Output: list of lists of ring indices, returns False if there is an issue
    """

    molecule = Chem.MolFromSmiles(SMILES_string, sanitize=False)

    if molecule is None:
        print("Invalid SMILES string.")
        return False
    else:
        # Detect all rings in the molecule
        ring_vector = rdmolops.GetSymmSSSR(molecule, False)
        ring_idx_list = []

        for ring in ring_vector:
            current_ring = []
            for idx in ring:
                current_ring.append(idx)
            ring_idx_list.append(current_ring)

    idx_to_index_dict = map_idx_to_string_index(SMILES_string)
    ring_indices_list = []
    for ring_idxs in ring_idx_list:
        ring_indices = idx_to_index_of_rings(ring_idxs, idx_to_index_dict, SMILES_string)

        ring_indices_list.append(ring_indices)

    return ring_indices_list



def map_idx_to_string_index(SMILES_string: str) -> dict:
    """
    Input: SMILES_string is a string of the SMILES representation
    Output: returns a dictionary that maps IDX values to indices in the SMILES_string
    """

    idx = 0
    idx_to_index_dict = {}
    for index in range(len(SMILES_string)):
        character = SMILES_string[index]
        if character.isalpha():
            idx_to_index_dict[idx] = index
            idx += 1

    return idx_to_index_dict



def idx_to_index_of_rings(ring_idx_list: list, idx_to_index_dict: dict, SMILES_string: str) -> list:
    """
    Input: ring_idx_list list of one ring's idx values. idx_to_index_dict dictionary
           that maps idx values to indices within the smiles string
    Output: list of indices of the ring
    """

    index_list = []
    for idx in ring_idx_list:
        index = idx_to_index_dict[idx]
        index_list.append(index)

        if index < len(SMILES_string) and SMILES_string[index + 1] == "(":
            index_list.append(index + 1)
        if index - 1 >= 0 and SMILES_string[index - 1] == ")":
            index_list.append(index - 1)

    return index_list



def create_dict_of_ring_indices(ring_indices_list: list) -> dict:
    """
    Input: ring_indices_list created by ring_handling_rdkit
    Output: dict that maps indices of rings to their ring's respective ring_num
    """

    ring_index_dict = {}
    ring_num = 1
    for ring_index_list in ring_indices_list:

        for i in ring_index_list:
            if i in ring_index_dict:
                ring_index_dict[i].append(ring_num)
            else:
                ring_index_dict[i] = [ring_num]

        ring_num += 1
    return ring_index_dict