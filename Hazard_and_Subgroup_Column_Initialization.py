import pandas as pd
import re


import sub_groups_from_smiles as subgroups
import Data_Collection_from_Pubchem as pbchem_coll



def create_dataframe_from_ids(compound_IDs: list,
                              save_to_csv: bool=True, wait_time: float = 1.5) -> pd.DataFrame|bool:
    """
    Input: compound_IDs list of ints representing compound IDs. save_to_csv bool to determine whether to
           save dataframe to csv or not, default is True. wait_time is a float value of the desired wait time between
           requests in seconds
    Output: pd.DataFrame of compound data including a column for each type of hazard, and columns of group
            and group combination hash keys
    """

    try:
        # List of compound identifiers
        all_data = []

        # Loop through each compound ID and fetch data
        for compound_ID in compound_IDs:
            compound_data = pbchem_coll.fetch_compound_data(compound_ID, wait_time)

            all_data.append(compound_data)

        df = pd.DataFrame(all_data)

        # Collect grouping data
        initialize_grouping_data(df)

        # Save for later use
        if save_to_csv:
            df.to_csv('compound_data.csv', index=False)

        return df
    except ValueError:
        return False



def update_existing_ids_dataframe(main_df: pd.DataFrame, new_compound_IDs: list,
                                   save_to_csv: bool=True, overwrite_old_data: bool=False) -> pd.DataFrame|bool:
    """
    Input: main_df is a pd.DataFrame that will be mutated to include new data. new_compound_IDS is list of ints
           representing compound IDs.  save_to_csv bool to determine whether to
           save dataframe to csv or not, default is True. overwrite_old_data is a bool that determine whether data for
           old compound IDs can
           be overwritten by new data. overwrite_old_data is False by default
    Output: The merged pd.DataFrame if parent_dataframe was successfully mutated to include new compound IDs, False
            otherwise
    """

    try:
        new_df = create_dataframe_from_ids(new_compound_IDs, False)
    except ValueError:
        print("Issue with new compound IDs")
        return False

    try:

        if overwrite_old_data:
            main_df = main_df[~main_df['Compound ID'].isin(new_df['Compound ID'])]
        else:
            new_df = new_df[~new_df['Compound ID'].isin(main_df['Compound ID'])]

        updated_df = pd.concat([main_df, new_df], axis=0, ignore_index=True).fillna(0)


        if save_to_csv:
            updated_df.to_csv('compound_data.csv', index=False)

        return updated_df

    except ValueError:
        print("Issue combining dataframes")
        return False





def generate_grouping_hash_lists(df: pd.DataFrame, Canonical_SMILES_column_name: str="Canonical SMILES",
                                 create_columns: bool=True) -> bool:
    """
    Input: df is a pd.DataFrame including a column of Canonical SMILES strings. Canonical_SMILES_column_name is the
           column name of the column containing the Canonical SMILES strings, it is "Canonical SMILES" by default.
           create_columns is a bool that determines if the columns "smiles_group_hash_list",
           "smiles_group_combination_hash_list", and "could_not_collect_grouping_data" need to be initialized.
           create_columns is True by default
    Output: bool that is True if columns where successfully created and or updated to include hash lists, False
            otherwise
    """

    try:
        if create_columns:
            df["smiles_group_hash_list"] = [[] for _ in range(len(df))]
            df["smiles_group_combination_hash_list"] = [[] for _ in range(len(df))]
            df["could_not_collect_grouping_data"] = 0


        for index, row in df.iterrows():
            canonical_smiles = row[Canonical_SMILES_column_name]
            group_list_info = subgroups.smiles_to_sub_groups(canonical_smiles)
            if group_list_info == 1:
                df.at[index, "could_not_collect_grouping_data"] = 1
                continue

            group_list = group_list_info[0]
            group_combination_list = group_list_info[1]

            smiles_group_hash_list = []
            smiles_group_combination_hash_list = []

            # get data for individual groups
            for group in group_list:
                group_hash_key = hash_smiles_group(group)
                smiles_group_hash_list.append(group_hash_key)

            for group_combo in group_combination_list:
                group_combo_hash_key = hash_smiles_group_combination(group_combo)
                smiles_group_combination_hash_list.append(group_combo_hash_key)

            df.at[index, "smiles_group_hash_list"] = smiles_group_combination_hash_list
            df.at[index, "smiles_group_combination_hash_list"] = smiles_group_hash_list

            return True
    except ValueError:
        return False



def hash_smiles_group(subgroup: str,
                      rotation_matters: bool=False) -> int|bool:
    """
    Input: subgroup is string representation of a subgroup. These should have been
           created by subgroups.smiles_to_sub_groups(). rotation_matters is a bool that represents whether
           groups other than rings that maintain the same order but are shifted can map to the same group, is False by
           default
    Output: hash representing the input subgroup

    Ex rotation_matters: if False, "ABCDE" and "DEABC" will map to the same hashcode and thus be counted as the same
                        group. If True, they will map to different hashcodes. Regardless of rotation_matters,
                        "0ABCDE" and "0DEABC" will map to the same hashcode as they both represent a ring of the same
                        components in the same order. The "0" at index 0 marks the subgroup as a ring.
    """

    # check if subgroup is empty
    if subgroup == "":
        return False

    is_ring = False
    if subgroup[0] == "0":
        is_ring = True
        subgroup = subgroup[1:]

    if not rotation_matters or is_ring:
        rotations = [subgroup[i:] + subgroup[:i] for i in range(len(subgroup))]
        pre_hashed_key = is_ring + min(rotations)

    else:
        pre_hashed_key = str(is_ring) + subgroup

    return hash(pre_hashed_key)



def hash_smiles_group_combination(group_combo: tuple,
                                  rotation_matters: bool=False) -> int:
    """
    Input: group_combo is a tuple of length two containing two subgroups of the same molecule. These should have been
           created by subgroups.smiles_to_sub_groups(). rotation_matters is a bool that is False by default.
           rotation_matters determines whether non ring subgroups with the same elements and the same order, but are
           shifted, represent the same subgroup.
    Output: Hash code to represent the subgroup combination
    """

    group_1_hash = hash_smiles_group(group_combo[0], rotation_matters)
    group_2_hash = hash_smiles_group(group_combo[1], rotation_matters)

    pre_hashed_key = min((group_1_hash + group_2_hash), (group_2_hash + group_1_hash))

    return hash(pre_hashed_key)



def create_grouping_columns(df: pd.DataFrame,
                            group_hashs_list_column_name: str="smiles_group_hash_list") -> bool:
    """
    Input: df is a pd.DataFrame containing a column that contains a list of subgroup hash's
    Output: bool of True if df was successfully mutated to include and update subgroup hashcode columns, False otherwise
    """

    try:
        for index, row in df.iterrows():
            hash_keys = row[group_hashs_list_column_name]

            # For each hash key in the list
            for hash_key in hash_keys:
                # Check if the column exists

                if hash_key in df.columns:
                    # Increment the value in the existing column
                    df.at[index, hash_key] += 1

                else:
                    # Create a new column and initialize it
                    df[hash_key] = 0  # Initialize the new column with zeros
                    df.at[index, hash_key] = 1

        return True

    except ValueError:
        return False



def initialize_grouping_data(df: pd.DataFrame,
                             Canonical_SMILES_column_name: str = "Canonical SMILES") -> bool:
    """
    Input: DataFrame containing a column with Canonical SMILES strings. Canonical_SMILES_column_name is the name
           of the column containing the Canonical SMILES representations.
    Output: True if df was successfully mutated to include hazard columns and subgroup columns, False otherwise
    """

    try:
        split_hazard_data(df)
        generate_grouping_hash_lists(df, Canonical_SMILES_column_name)
        create_grouping_columns(df)
        return True
    except ValueError:
        return False



def split_hazard_data(df: pd.DataFrame,
                      hazard_column_name: str="Hazards") -> bool:
    """
    Input: pd.DataFrame containing a column that contains compound hazards
    Output: True if df was successfully mutated to include a column for each type of hazard, False otherwise
    """

    try:
        for index, row in df.iterrows():
            hazards_string = row[hazard_column_name]
            if type(hazards_string) != str:
                continue

            hazard_list = re.sub(r'\s+', '', hazards_string).split(",")

            # For each hazard in the list
            for hazard in hazard_list:

                # Check if the column exists
                if hazard in df.columns:
                    # Set row value to True
                    df.at[index, hazard] = True

                else:
                    # Create a new column and initialize it
                    df[hazard] = 0  # Initialize the new column with 0 to represent False
                    df.at[index, hazard] = 1   # set this row's value to 1 to represent True

        return True

    except ValueError:
        return False