import pandas as pd
from sklearn.neighbors import NearestNeighbors

import Hazard_and_Subgroup_Column_Initialization as data_init

import warnings

# Suppress the specific PerformanceWarning from pandas
# These performance errors will be looked into in the future, but for now it works
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

# Initializes default_dataframe that comes with the package
try:
    default_data_df = pd.read_csv("default_compound_data.csv")
except FileNotFoundError:
    default_data_df = None

# Initializes default_testing_dataframe that comes with the package
try:
    default_testing_data_df = pd.read_csv("default_testing_compound_data.csv")
except FileNotFoundError:
    default_data_df = None



def determine_potential_hazards_from_dataframe(canonical_smiles: str,
                                               main_df: pd.DataFrame=default_data_df,
                                               min_probability: float=.5, num_neighbors: int=10, algorithm: str= "auto") \
                                               -> pd.DataFrame:
    """
    Input: canonical_smiles is a str. main_df is the dataframe that will be used to fit the model that predicts the
           most similar compounds. min_probability is a float value between 0 and 1 that determines what proportion of
           nearest neighbors must have a hazard for the model to predict that the input molecule will have said hazard.
           num_neighbors is the number of nearest neighbors the model pulls to determine the potential hazards.
           algorithm is a str that represents the algorithm the model will use.
    Output: pd.DataFrame containing columns for each hazard for the input compound. 1 indicates that the model predicts
            that compound will have said hazard. 0 predicts the opposite.
    """

    main_df_cloned = main_df.copy()
    this_smiles_df = data_init.convert_smiles_to_dataframe(canonical_smiles)

    cleaned_main_df = clean_data_frame(main_df, limit_num_groups=False)

    df_tuple = data_init.fit_dataframes(this_smiles_df, cleaned_main_df)

    this_smiles_finished_df = df_tuple[0]
    finished_main_df = df_tuple[1]

    # convert column names to string so that scikit can read them
    this_smiles_finished_df.columns = this_smiles_finished_df.columns.astype(str)
    finished_main_df.columns = finished_main_df.columns.astype(str)

    nn_model = NearestNeighbors(n_neighbors=num_neighbors, algorithm=algorithm).fit(finished_main_df)
    neighbors_indices = nn_model.kneighbors(this_smiles_finished_df)[1]

    hazard_probability_dict = pull_hazards_from_dataframe(neighbors_indices[0], main_df_cloned)

    # set to 1 if the probability is high enough
    for hazard in hazard_probability_dict:
        if hazard_probability_dict[hazard] >= min_probability:
            hazard_probability_dict[hazard] = [1]
        else:
            hazard_probability_dict[hazard] = [0]


    this_smiles_hazard_df = pd.DataFrame(hazard_probability_dict)

    return this_smiles_hazard_df



def clean_data_frame(main_df: pd.DataFrame,
                     limit_num_groups: bool=True, num_groups_to_include: int=200,
                     group_counts_matters: bool=False, columns_to_drop=None) -> pd.DataFrame:
    """
    Input: main_df is the dataframe that will be used to fit the model that predicts the
           most similar compounds. limit_num_groups is a bool that determines whether to limit the number of
           subgroups and subgroup combinations hash columns. if limit_num_groups = True, then num_groups_to_include
           is an int that determines how many groups to include. group_counts_matters determines whether the
           number of occurrences a subgroup/subgroup combination has in a molecule has any effect on the model.
           columns_to_drop is a list of names of columns to drop that are unneeded (this is for all columns other than
           hash columns)
    Output: pd.DataFrame containing only the desired columns to be fit with the model
    """

    # drop molecules that did not have any hazard data associated with them
    if columns_to_drop is None:
        columns_to_drop = ["IUPAC Name", "Canonical SMILES", "Hazards", "AcuteToxic",
                           "HealthHazard", "Irritant", "Corrosive", "EnvironmentalHazard", "CompressedGas",
                           "could_not_collect_grouping_data", "no_hazard_data_to_collect"]
    training_df = main_df[main_df["Hazards"] != ""]

    # drop molecules that were unsuccessfully parsed for grouping data
    training_df = training_df[training_df["could_not_collect_grouping_data"] != 1]

    # drop unwanted columns
    training_df = training_df.drop(columns=columns_to_drop, errors='ignore')

    # converts all non 0 values to 1 if group_count_matters
    # Ex: if group_counts_matters is True, then a molecule containing 100 =O groups and a molecule
    # containing 1 =O group will have both be set to 1 for  the column corresponding to =O
    if not group_counts_matters:
                                    # .applymap() used inplace of .map() because of google collab having issues with .map()
            training_df = training_df.applymap(lambda x: 1 if x != 0 else 0)


    # Calculate the number of non-zero values for each column
    non_zero_counts = (training_df != 0).sum()

    # Sort columns by the number of non-zero values in descending order
    sorted_columns = non_zero_counts.sort_values(ascending=False)

    if limit_num_groups:
        top_columns = sorted_columns.head(num_groups_to_include).index
    else:
        # does not limit the num of subgroup hash columns
        top_columns = sorted_columns.index

    finished_training_df = training_df[top_columns]

    return finished_training_df



def pull_hazards_from_dataframe(index_list: list, main_df: pd.DataFrame,
                                hazard_list=None) -> dict:
    """
    Input: index_list is a list of the indices of the nearest neighbors. main_df is the pd.DataFrame that was used in
           fitting the model.  hazard_list is a list of potential hazard names.
    Output: dict with hazard names mapping to the proportion of neighbors who have that hazard
    """

    if hazard_list is None:
        hazard_list = ["AcuteToxic", "Flammable", "HealthHazard", "Irritant", "Corrosive",
                       "EnvironmentalHazard", "CompressedGas"]

    filtered_df = main_df.loc[index_list]
    len_index_list = len(index_list) # this is the number of neighbors we are using
    hazard_probability_dict = {}

    for hazard in hazard_list:
        if hazard in main_df.columns:
            hazard_count = filtered_df[hazard].sum()
            hazard_prob = hazard_count / len_index_list  # of a neighbor to have this hazard

        else:
            hazard_prob = 0

        hazard_probability_dict[hazard] = hazard_prob

    return hazard_probability_dict



def hazard_model_testing(testing_df: pd.DataFrame=default_testing_data_df,
                         main_df: pd.DataFrame=default_data_df, min_probability: float=.5,
                         num_neighbors: int=10, algorithm: str= "auto") -> pd.DataFrame:
    """
    Input: testing_df is the pd.DataFrame that will be used to test a parameter set's effectiveness (The rest of these
           optional inputs are said parameters that are being tested). main_df is the dataframe that will be used to fit
           the model that predicts the most similar compounds. min_probability is a float value between 0 and 1 that
           determines what proportion of nearest neighbors must have a hazard for the model to predict that the input
           molecule will have said hazard. num_neighbors is the number of nearest neighbors the model pulls to determine
           the potential hazards. algorithm is a str that represents the algorithm the model will use.
    Output: pd.DataFrame displaying correct prediction proportion, false negative proportion, and false positive
            proportion
    """

    len_testing_data = len(testing_df)  # Calculate the length of the testing DataFrame

    hazard_prediction_accuracy_count_dict = {}
    hazard_prediction_false_pos_count_dict = {}
    hazard_prediction_false_neg_count_dict = {}

    for index, row in testing_df.iterrows():

        this_canonical_smiles = row['Canonical SMILES']

        this_predicted_hazard_df = determine_potential_hazards_from_dataframe(this_canonical_smiles, main_df=main_df,
                                                                               min_probability=min_probability,
                                                                               num_neighbors=num_neighbors,
                                                                               algorithm=algorithm)

        predicted_hazard_column_name_list = this_predicted_hazard_df.columns

        for hazard in predicted_hazard_column_name_list:
            # prediction is either 1 or 0
            prediction = this_predicted_hazard_df[hazard].iloc[0]
            # actual is either 1 or 0
            actual = row.get(hazard)
            if type(actual) != int:
                actual = 0

            # False negative
            if prediction < actual:
                if hazard not in hazard_prediction_false_neg_count_dict:
                    hazard_prediction_false_neg_count_dict[hazard] = 1
                else:
                    hazard_prediction_false_neg_count_dict[hazard] += 1

            # Correct prediction
            elif prediction == actual:
                # add 1 for a correct prediction
                if hazard not in hazard_prediction_accuracy_count_dict:
                    hazard_prediction_accuracy_count_dict[hazard] = 1
                else:
                    hazard_prediction_accuracy_count_dict[hazard] += 1

            # False positive
            elif prediction > actual:
                if hazard not in hazard_prediction_false_pos_count_dict:
                    hazard_prediction_false_pos_count_dict[hazard] = 1
                else:
                    hazard_prediction_false_pos_count_dict[hazard] += 1

    # Calculate proportions based on the total number of rows in the testing DataFrame
    hazard_accuracy_dict = {}
    hazard_accuracy_dict["Statistic"] = "Correct Prediction Proportion"
    for hazard in hazard_prediction_accuracy_count_dict:
        num_correct = hazard_prediction_accuracy_count_dict[hazard]
        correctness = num_correct / len_testing_data  # Using len_testing_data here
        hazard_accuracy_dict[hazard] = correctness

    hazard_false_pos_dict = {}
    hazard_false_pos_dict["Statistic"] = "False Positive Proportion"
    for hazard in hazard_prediction_false_pos_count_dict:

        num_correct = hazard_prediction_false_pos_count_dict[hazard]
        correctness = num_correct / len_testing_data  # Using len_testing_data here
        hazard_false_pos_dict[hazard] = correctness

    hazard_false_neg_dict = {}
    hazard_false_neg_dict["Statistic"] = "False Negative Proportion"
    for hazard in hazard_prediction_false_neg_count_dict:

        num_correct = hazard_prediction_false_neg_count_dict[hazard]
        correctness = num_correct / len_testing_data  # Using len_testing_data here
        hazard_false_neg_dict[hazard] = correctness

    hazard_probability_df = pd.DataFrame([hazard_accuracy_dict, hazard_false_pos_dict, hazard_false_neg_dict])
    return hazard_probability_df