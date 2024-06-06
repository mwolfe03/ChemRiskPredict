import pandas as pd
import re

import sub_groups_from_smiles as subgroups
import Data_Collection_from_Pubchem as pubchem_coll

import warnings

# Suppress the specific PerformanceWarning from pandas
# These performance errors will be looked into in the future, but for now it works
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
warnings.simplefilter('ignore', category=pd.errors.SettingWithCopyWarning)
# this gets raised by .applymap(), but the version of pandas on google collab crashes using .map()
warnings.simplefilter('ignore', category=FutureWarning)



def create_dataframe_from_cids(compound_IDs: list,
                              save_to_csv: bool=True, csv_name: str= 'compound_data.csv', wait_time: float=1,
                              drop_empty_hazard_rows: bool=True) -> pd.DataFrame|bool:
    """
    Input: compound_IDs list of ints representing compound IDs. save_to_csv bool to determine whether to
           save dataframe to csv or not, default is True. wait_time is a float value of the desired wait time between
           requests in seconds. drop_empty_hazard_rows is a bool that is True by default. drops rows in the columns that
           do not have any hazard data associated with them. This is due to the inability to differentiate between
           chemicals that have no hazards because they are safe and chemicals that have no hazards because there simply
           is not enough data.
    Output: pd.DataFrame of compound data including a column for each type of hazard, and columns of group
            and group combination hash keys
    """

    try:
        # List of compound identifiers
        all_data = []

        # Loop through each compound ID and fetch data
        for compound_ID in compound_IDs:
            compound_data = pubchem_coll.fetch_compound_data(compound_ID, wait_time)

            all_data.append(compound_data)

        df = pd.DataFrame(all_data)

        # Collect grouping data
        initialize_grouping_data(df)

        # Removes rows that don't have any hazards attached. Currently, there is no way to differentiate
        # between that chemicals that are safe, and chemicals that simply do not have any hazard data.
        if drop_empty_hazard_rows:
            df = df[df['Hazards'] != '']

        # Save for later use
        if save_to_csv:
            df.to_csv(csv_name, index=False)
        columns_to_remove = ["smiles_group_hash_list", "smiles_group_combination_hash_list"]
        df = df.drop(columns=columns_to_remove, errors='ignore')

        return df
    except ValueError:
        return False



def generate_grouping_hash_lists_from_df(df: pd.DataFrame,
                                         Canonical_SMILES_column_name: str="Canonical SMILES",
                                         create_columns: bool=True) -> bool:
    """
    Input: df is a pd.DataFrame including a column of Canonical SMILES strings. Canonical_SMILES_column_name is the
           column name of the column containing the Canonical SMILES strings, it is "Canonical SMILES" by default.
           create_columns is a bool that determines if the columns "smiles_group_hash_list",
           "smiles_group_combination_hash_list", and "could_not_collect_grouping_data" need to be initialized.
           create_columns is True by default
    Output: bool that is True if columns were successfully created and/or updated to include hash lists, False
            otherwise
    """

    try:
        if create_columns:
            df["smiles_group_hash_list"] = [[] for _ in range(len(df))]
            df["smiles_group_combination_hash_list"] = [[] for _ in range(len(df))]
            df["could_not_collect_grouping_data"] = 0

        for index, row in df.iterrows():
            canonical_smiles = row[Canonical_SMILES_column_name]

            group_hash_lists = generate_grouping_hash_lists_from_SMILES_string(canonical_smiles)

            if group_hash_lists == False:
                df.at[index, "could_not_collect_grouping_data"] = 1
                continue

            smiles_group_hash_list = group_hash_lists[0]
            smiles_group_combination_hash_list = group_hash_lists[1]

            df.at[index, "smiles_group_hash_list"] = smiles_group_hash_list
            df.at[index, "smiles_group_combination_hash_list"] = smiles_group_combination_hash_list

        return True
    except ValueError:
        return False



def generate_grouping_hash_lists_from_SMILES_string(canonical_smiles: str) -> tuple|bool:
    """
    Input: canonical_smiles of type str
    Output: tuple containing a hashlist of the different subgroups at index 0, and a hashlist of the different subgroup
            combinations at index 1. False if no subgroups could be pulled
    """
    group_list_info = subgroups.smiles_to_sub_groups(canonical_smiles)

    # group data could not be pulled
    if not group_list_info:
        return False

    group_list = group_list_info[0]
    group_combination_list = group_list_info[1]

    smiles_group_hash_list = []
    smiles_group_combination_hash_list = []

    # get data for individual groups
    for group in group_list:
        group_hash_key = hash_smiles_group(group)
        smiles_group_hash_list.append(group_hash_key)

    # get data for group combinations
    for group_combo in group_combination_list:
        group_combo_hash_key = hash_smiles_group_combination(group_combo)
        smiles_group_combination_hash_list.append(group_combo_hash_key)

    return smiles_group_hash_list, smiles_group_combination_hash_list



def convert_smiles_to_dataframe(canonical_smiles: str) -> pd.DataFrame:
    """
    Input: canonical_smiles is a string.
    Output: A pd.DataFrame with only one row referring to canonical_smiles. pd.Dataframe contains columns of the
            hashed subgroups and subgroup combinations
    """

    hash_dict = {}
    hash_lists = generate_grouping_hash_lists_from_SMILES_string(canonical_smiles)
    hash_lists_combined = hash_lists[0] + hash_lists[1]
    for hash in hash_lists_combined:
        if hash in hash_dict:
            hash_dict[hash][0] += 1
        else:
            hash_dict[hash] = [1]

    df = pd.DataFrame(hash_dict)

    return df



def fit_dataframes(this_smiles_df: pd.DataFrame, cleaned_main_df: pd.DataFrame) -> tuple:
    """
    Input: this_smiles_df is a pd.DataFrame. cleaned_main_df is a pd.DataFrame
    Output: tuple containing this_smiles_df at index 0 and cleaned_main_df where both dataframes have each other's
            columns added and set to zero. Both DataFrames have the same columns and are in the same order
    """

    # removes all columns not in this_smiles_df that are sum to less than 2
    cleaned_filtered_main_df = filter_columns_by_sum_and_input(cleaned_main_df, this_smiles_df)

    # Get the common columns between this_smiles_df and cleaned_main_df
    common_columns = this_smiles_df.columns.intersection(cleaned_main_df.columns)

    # Only keep the common columns in this_smiles_df
    this_smiles_df = this_smiles_df[common_columns]

    # Add missing columns from cleaned_filtered_main_df to this_smiles_df and set their values to zero
    missing_columns = [col for col in cleaned_filtered_main_df.columns if col not in this_smiles_df.columns]
    for col in missing_columns:
        this_smiles_df[col] = 0

    this_smiles_df = this_smiles_df[cleaned_filtered_main_df.columns]
    # Reorder the columns of this_smiles_df to match cleaned_main_df

    return (this_smiles_df, cleaned_filtered_main_df)



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
        pre_hashed_key = str(is_ring) + min(rotations)

    else:
        pre_hashed_key = str(is_ring) + subgroup

    hashed = hash(pre_hashed_key)
    return hashed



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
        generate_grouping_hash_lists_from_df(df, Canonical_SMILES_column_name)
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
                    df.at[index, hazard] = 1

                else:
                    # Create a new column and initialize it
                    df[hazard] = 0  # Initialize the new column with 0 to represent False
                    df.at[index, hazard] = 1   # set this row's value to 1 to represent True

        return True
    except ValueError:
        return False



def update_existing_ids_dataframe_from_cids(main_df: pd.DataFrame, new_compound_IDs: list,
                                            save_to_csv: bool=True, csv_name: str="compound_data.csv",
                                            overwrite_old_data: bool=False) -> pd.DataFrame|bool:
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
        new_df = create_dataframe_from_cids(new_compound_IDs, save_to_csv=False)

    except ValueError:
        print("Issue with new compound IDs")
        return False

    try:
        return update_existing_dataframe_from_dataframe(main_df, new_df, save_to_csv=save_to_csv, csv_name=csv_name, overwrite_old_data=overwrite_old_data)

    except ValueError:
        print("Issue combining dataframes")
        return False



def update_existing_dataframe_from_dataframe(main_df: pd.DataFrame, second_df: pd.DataFrame,
                                             save_to_csv: bool=True, csv_name: str="compound_data.csv",
                                             overwrite_old_data: bool=False) -> pd.DataFrame|bool:
    """
    Input: main_df is a pd.DataFrame that will be mutated to include new data. second_df is a pd.DataFrame that will be
           used to update main_df. save_to_csv bool to determine whether to save dataframe to csv or not, default is
           True. overwrite_old_data is a bool that determine whether data for old compound IDs can be overwritten by
           new data. overwrite_old_data is False by default
    Output: The merged pd.DataFrame if parent_dataframe was successfully mutated to include new compound IDs, False
            otherwise
    """

    try:
        if "Compound ID" not in main_df.columns or "Compound ID" not in second_df.columns:
            raise ValueError("Both dataframes must have a 'Compound ID' column")

        if overwrite_old_data:

            matching_ids = second_df["Compound ID"].isin(main_df["Compound ID"])
            matching_rows = second_df[matching_ids]

            # this is what will be iterated over to add missing columns

            second_df = second_df[~matching_ids]
            main_df = main_df.set_index("Compound ID")
            matching_rows = matching_rows.set_index("Compound ID")
            main_df.update(matching_rows)
            main_df = main_df.reset_index()

        else:
            # filters out rows in second_df that have matching "Compound ID" with values in main_df
            second_df = second_df[~second_df["Compound ID"].isin(main_df["Compound ID"])]

        # Step 1: Add missing columns from main_df to second_df and set them to 0
        for col in main_df.columns:
            if col not in second_df.columns:
                second_df.loc[:, col] = 0

        # Step 2: Add missing columns from second_df to main_df and set them to 0
        for col in second_df.columns:
            if col not in main_df.columns:
                main_df.loc[:, col] = 0

        # Step 3: Ensure both dataframes have the same columns and order
        main_df = main_df[sorted(main_df.columns, key=str)]
        second_df = second_df[sorted(second_df.columns, key=str)]

        # Step 4: Concatenate the dataframes
        combined_df = pd.concat([main_df, second_df], ignore_index=True)

        if save_to_csv:
            combined_df.to_csv(csv_name, index=False)

        return combined_df

    except ValueError:
        print("Issue combining dataframes")
        return False



def filter_columns_by_sum_and_input(main_df: pd.DataFrame, second_df: pd.DataFrame) -> pd.DataFrame:
    """
    Input: main_df is a pd.DataFrame, second_df pd.DataFrame
    Output: updated main_df that contains only columns that have more than one "1" or are present in second_df
    """
    columns_to_keep = []
    for col in main_df.columns:
        if main_df[col].sum() >= 2 or col in second_df.columns:
            columns_to_keep.append(col)
    return main_df[columns_to_keep]