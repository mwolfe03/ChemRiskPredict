import pandas as pd
from sklearn.neighbors import NearestNeighbors
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline

from dataframe_and_model_creation import Hazard_and_Subgroup_Column_Initialization as data_init



# Initializes default_dataframe that comes with the package
try:
    default_data_df = pd.read_csv("default_data/training_compound_data.csv")
except FileNotFoundError:
    print("failed to load default training data")
    default_data_df = None

# Initializes default_testing_dataframe that comes with the package
try:
    default_testing_data_df = pd.read_csv("default_data/testing_compound_data.csv")
except FileNotFoundError:
    print("failed to load default testing data")
    default_testing_data_df = None



def determine_potential_hazards_from_dataframe(canonical_smiles: str,
                                               main_df: pd.DataFrame=default_data_df,
                                               min_probability: float=.3, num_neighbors: int=5, algorithm: str= "auto",
                                               metric: str="cosine") \
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

    # DATA INITIALIZATION
    #start_time = time.time()
    # function returning a dict of all the different subgroups and subgroup combinations with a 1. Turn them into strings
    this_smiles_dict = data_init.convert_smiles_to_dict(canonical_smiles)
    cleaned_cloned_main_df = main_df.filter(like='key', axis=1)  # include only key columns
    # Convert the input dictionary to a DataFrame and reindex to include all columns
    # function removing any columns in the dict that are not IN the main df

    this_smiles_dict_cleaned = data_init.clean_smiles_dict(this_smiles_dict, cleaned_cloned_main_df)

    this_smiles_df = pd.DataFrame([this_smiles_dict_cleaned])
    this_smiles_df = this_smiles_df.reindex(columns=cleaned_cloned_main_df.columns, fill_value=0)

    # MODEL BUILDING
    imputer = SimpleImputer(strategy='constant', fill_value=0)
    knn = NearestNeighbors(n_neighbors=num_neighbors, metric=metric)  # NearestNeighbors instead of KNeighborsClassifier
    pipeline = Pipeline(steps=[('imputer', imputer), ('knn', knn)])
    pipeline.fit(cleaned_cloned_main_df)
    #print("--- %s end Model ---" % (time.time() - start_time))
    # MODEL PREDICTIONS
    neighbor_distances, neighbor_indices = pipeline.named_steps['knn'].kneighbors(this_smiles_df.iloc[0].to_numpy().reshape(1, -1))
    #print(neighbor_distances)
    #print(neighbor_indices)
    #print("--- %s predict neighbors ---" % (time.time() - start_time))
    hazard_probability_dict = pull_hazards_from_dataframe(neighbor_indices[0], main_df)

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
           subgroups and subgroup combinations group columns. if limit_num_groups = True, then num_groups_to_include
           is an int that determines how many groups to include. group_counts_matters determines whether the
           number of occurrences a subgroup/subgroup combination has in a molecule has any effect on the model.
           columns_to_drop is a list of names of columns to drop that are unneeded (this is for all columns other than
           group columns)
    Output: pd.DataFrame containing only the desired columns to be fit with the model
    """
    # drop molecules that did not have any hazard data associated with them
    if columns_to_drop is None:
        columns_to_drop = ["IUPAC Name", "Canonical SMILES", "Hazards",
                           "could_not_collect_grouping_data", "no_hazard_data_to_collect", ]
    training_df = main_df[main_df["Hazards"] != ""] #check, might not need


    # drop molecules that were unsuccessfully parsed for grouping data
    training_df = training_df[training_df["could_not_collect_grouping_data"] != 1]
    # drop unwanted columns
    training_df = training_df.drop(columns=columns_to_drop, errors='ignore')
    # converts all non 0 values to 1 if group_count_matters
    # Ex: if group_counts_matters is True, then a molecule containing 100 =O groups and a molecule
    # containing 1 =O group will have both be set to 1 for  the column corresponding to =O
    #if not group_counts_matters:
                                    # .applymap() used inplace of .map() because of google collab having issues with .map()
           # training_df = training_df.applymap(lambda x: 1 if x != 0 else 0)
    #print("--- %s BC4 ---" % (time.time() - start_time))

    # Calculate the number of non-zero values for each column
    non_zero_counts = (training_df != 0).sum()
    #print("--- %s BC5 ---" % (time.time() - start_time))
    # Sort columns by the number of non-zero values in descending order
    sorted_columns = non_zero_counts.sort_values(ascending=False)
    #print("--- %s BC6 ---" % (time.time() - start_time))
    if limit_num_groups:
        top_columns = sorted_columns.head(num_groups_to_include).index
    else:
        # does not limit the num of subgroup group columns
        top_columns = sorted_columns.index
    #print("--- %s BC7 ---" % (time.time() - start_time))
    finished_training_df = training_df[top_columns]
    #print("--- %s BC8 ---" % (time.time() - start_time))
    return finished_training_df



def pull_hazards_from_dataframe(index_list: list, main_df: pd.DataFrame,
                                hazard_list=None) -> dict:
    """
    Input: index_list is a list of the indices of the nearest neighbors. main_df is the pd.DataFrame that was used in
           fitting the model.  hazard_list is a list of potential hazard names.
    Output: dict with hazard names mapping to the proportion of neighbors who have that hazard
    """

    if hazard_list is None:

        hazard_list = ["AcuteToxic", "Flammable", "HealthHazard", "Corrosive",
                       "EnvironmentalHazard"]

    filtered_df = main_df.iloc[index_list]
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
                         main_df: pd.DataFrame=default_data_df, min_probability: float=.3,
                         num_neighbors: int=5, algorithm: str= "auto", metric: str="cosine") -> pd.DataFrame:
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
    hazard_prediction_true_pos_count_dict = {}
    hazard_prediction_true_neg_count_dict = {}

    hazard_prediction_false_pos_count_dict = {}
    hazard_prediction_false_neg_count_dict = {}

    for index, row in testing_df.iterrows():

        this_canonical_smiles = row['Canonical SMILES']

        this_predicted_hazard_df = determine_potential_hazards_from_dataframe(this_canonical_smiles, main_df=main_df,
                                                                               min_probability=min_probability,
                                                                               num_neighbors=num_neighbors,
                                                                               algorithm=algorithm, metric=metric)

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

                # True negative
                if prediction == 0:
                    if hazard not in hazard_prediction_true_neg_count_dict:
                        hazard_prediction_true_neg_count_dict[hazard] = 1
                    else:
                        hazard_prediction_true_neg_count_dict[hazard] += 1
                # True positive
                if prediction == 1:
                    if hazard not in hazard_prediction_true_pos_count_dict:
                        hazard_prediction_true_pos_count_dict[hazard] = 1
                    else:
                        hazard_prediction_true_pos_count_dict[hazard] += 1


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

    hazard_prediction_true_neg_dict = {}
    hazard_prediction_true_neg_dict["Statistic"] = "True Negative Proportion"
    for hazard in hazard_prediction_true_neg_count_dict:
        num_correct = hazard_prediction_true_neg_count_dict[hazard]
        correctness = num_correct / len_testing_data  # Using len_testing_data here
        hazard_prediction_true_neg_dict[hazard] = correctness

    hazard_prediction_true_pos_dict = {}
    hazard_prediction_true_pos_dict["Statistic"] = "True Positive Proportion"
    for hazard in hazard_prediction_true_pos_count_dict:
        num_correct = hazard_prediction_true_pos_count_dict[hazard]
        correctness = num_correct / len_testing_data  # Using len_testing_data here
        hazard_prediction_true_pos_dict[hazard] = correctness


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

    hazard_probability_df = pd.DataFrame([hazard_accuracy_dict, hazard_prediction_true_pos_dict,
                                          hazard_prediction_true_neg_dict, hazard_false_pos_dict, hazard_false_neg_dict])
    return hazard_probability_df