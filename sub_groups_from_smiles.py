import ring_handling_w_rdkit as ring_handle



def smiles_to_sub_groups(SMILES_string : str) -> tuple|bool:
    """
    Input: SMILES_string smiles representation of a molecule in string format
    Output: tuple containing tuple of all subgroups of the given molecule at index 0,
            and tuple of all combinations of adjacent subgroups at index 1
    Ex Input: "CC1=C(C(=CC=C1)[N+](=O)[O-])N"
    Ex Output: "(('0CCCCCC', 'C', '[N+]', '=O', '[O-]', 'N'), (('CCCCCC', 'C'), ('[N+]', 'CCCCCC'), ('[N+]', '=O'), ('[O-]', '[N+]'), ('N', 'CCCCCC')))"
    Note: subgroups that represent rings will have "0" at index 0 as a marker
    """

    if type(SMILES_string) != str:
        print("SMILES_string is of wrong type")
        return False

    try:
        sub_group_index_list = subgroups_from_SMILES(SMILES_string)
        test_neighbors_dict = determine_group_neighbors(SMILES_string, sub_group_index_list)
        cleaned_dict = clean_subgroup_dict(test_neighbors_dict, SMILES_string)
        group_list = get_group_list(cleaned_dict)
        combination_groups = create_group_combinations(cleaned_dict)

        return group_list, combination_groups

    except ValueError:
        print("Invalid SMILES_string")
        return False



def subgroups_from_SMILES(SMILES_string: str) -> list:
    subgroup_index_list = []
    index_to_group_dict = {}
    if "1" in SMILES_string: # checks if there are any rings
        ring_indices_lists = ring_handle.ring_handling_rdkit(SMILES_string)
        index_to_group_dict = ring_handle.create_dict_of_ring_indices(ring_indices_lists)
        ring_num = 1
        for ring_indices in ring_indices_lists:
            subgroup_index_list.append([ring_num, ring_indices])
            ring_num += 1

    bond_type_tuple = ("=", "#")
    grouping_type_tuple = ("(", ")")

    current_sub_group = []
    current_group_ascii = 97
    current_group_name = chr(current_group_ascii)

    SMILES_length = len(SMILES_string)
    for i in range(SMILES_length):
        char = SMILES_string[i]
        if i in index_to_group_dict:
            # With the exception of rings, no atom/bond index should appear in more than one group.
            # An atom/bond index can be in 1 group OR in 1 or more rings
            # everything before this should be it's own group
            if len(current_sub_group) > 0:
                current_group_ascii += 1
                current_group_name = chr(current_group_ascii)
                subgroup_index_list.append((current_group_name, current_sub_group))
            current_sub_group = []
            continue

        elif char.isdigit():
            # no reason to include digits
            # everything before this should be its own group
            # in theory, current_sub_group should always be a list of len 0 here
            if len(current_sub_group) > 0:
                subgroup_index_list.append(current_sub_group)
                current_group_name = chr(current_group_ascii)
                subgroup_index_list.append((current_group_name, current_sub_group))
            current_sub_group = []
            continue

        elif char in grouping_type_tuple:
            if char == "(":
                current_sub_group.append(i)     # including opening now helps with finding neighbors later
            # leaves as list for now
            subgroup_index_list.append((current_group_name, current_sub_group))
            previous_sub_group = current_sub_group
            current_sub_group = []

            if char == ")":
                current_sub_group.append(i)

            current_group_ascii += 1
            current_group_name = chr(current_group_ascii)
        else:    # wanted character types to represent a subgroup
            current_sub_group.append(i)  # add index of char to list
            index_to_group_dict[i] = [current_group_name]

        if i+1 == SMILES_length: # for items at the very end of the string
            # leave as a list for now
            subgroup_index_list.append((current_group_name, current_sub_group))

    subgroup_index_list = clean_subgroup_index_list(subgroup_index_list, SMILES_string)
    return subgroup_index_list



def clean_subgroup_index_list(subgroup_index_list: list, SMILES_string: str) -> list:
    """
    Input: subgroup_index_list is a list of lists of subgroup indices
    Output: subgroup_index_list with empty index tuples and single item tuples of non atoms removed
    """

    groups_to_remove = []
    for subgroup_tuple in subgroup_index_list:
        tuple_of_indices = subgroup_tuple[1] # index 0 is the group name
        length_tuple_of_indices = len(tuple_of_indices)

        if length_tuple_of_indices == 0:
            groups_to_remove.append((subgroup_tuple))

        elif length_tuple_of_indices == 1 and (not SMILES_string[tuple_of_indices[0]].isalpha()):   # remove any unhelpful single item tuples
            groups_to_remove.append((subgroup_tuple))

    for bad_group in groups_to_remove:
        subgroup_index_list.remove(bad_group)
    return subgroup_index_list



def create_group_dict(subgroup_index_list: list) -> dict:
    """
    Input: subgroup_index_list is a list of lists representing all groups. Each interior list contains the
           group's symbol at [0] and indices list at [1]
    Output: dictionary in the format {group_symbol: [[group_symbol, [group_indices]], [empty neighbor_list]]}
    """

    subgroup_dict = {}
    for subgroup in subgroup_index_list:
        subgroup_dict[subgroup[0]] = [subgroup, []] # empty second list will hold neighbors
    return subgroup_dict



def grouping_bridges_dict(SMILES_string: str) -> dict:
    """
    Input: Smiles string
    Output: list of the index BEFORE starting parens, and index AFTER after corresponding end parens
    Ex Input: "Cc(c(o)ccccc)C"
    Ex Output: {2: 12, 4: 6}
    """
    len_smiles_string = len(SMILES_string)
    group_bridge_dict = {}
    for i in range(len_smiles_string):
        char = SMILES_string[i]

        if char == "(":

            paren_num = 1
            start_paren_index = i
            for index in range(i + 1, len_smiles_string):
                current_char = SMILES_string[index]

                if current_char == "(":
                    paren_num += 1
                elif current_char == ")":
                    paren_num -= 1

                if paren_num == 0 and current_char == ")":
                    end_paren_index = index
                    group_bridge_dict[start_paren_index] = end_paren_index
                    break

    return group_bridge_dict



def determine_group_neighbors(SMILES_string: str, subgroup_index_list: list) -> dict:
    """
    Input: SMILES_string is a string of the SMILES representation. subgroup_index_list is a list of lists of subgroup indices.
           subgroup_index_list should have been created by subgroups_from_SMILES()
    Output: dict mapping subgroups to adjacent/ overlapping subgroups. Connections are one directional. dict follows
            the format of {group_symbol: [[group_symbol, [group_indices]], [neighbor_list]]}
    """
    group_bridge_dict = grouping_bridges_dict(SMILES_string)
    subgroup_dict = create_group_dict(subgroup_index_list)
    for this_subgroup in subgroup_index_list:
        this_subgroup_symbol = this_subgroup[0]
        this_subgroup_indices = this_subgroup[1]

        for i in this_subgroup_indices:
            if i in group_bridge_dict:
                # if i is a key in this dict, then it is a start paren "(" index
                # "(" index maps to corresponding ")" index
                bridged_index = group_bridge_dict[i]

                # add to list so that this group, and the group originally containing the ")" index are marked as neighbors
                this_subgroup_indices.append(bridged_index)


        for other_subgroup in subgroup_index_list:

            for index in this_subgroup_indices:
                neighbor_index_back = index - 1
                neighbor_index_front = index + 1

                other_subgroup_index_list = other_subgroup[1]
                other_subgroup_symbol = other_subgroup[0]
                other_subgroup_neighbors = subgroup_dict[other_subgroup_symbol][1]

                if SMILES_string[index] == ")":
                    continue

                if neighbor_index_front in other_subgroup_index_list and SMILES_string[neighbor_index_front] == ")":
                    continue

                if this_subgroup_symbol == other_subgroup_symbol:
                    break

                elif  this_subgroup_symbol in other_subgroup_neighbors:
                    # We want each neighbor connection to be one directional
                    # so that there are no duplicate tuples later on.

                    break

                elif (index in other_subgroup_index_list
                        or neighbor_index_back in other_subgroup_index_list
                        or neighbor_index_front in other_subgroup_index_list):
                    # Then they are neighbors
                    # Example subgroup: key -> [(key, [index_list]), [neighbor_symbol_list]

                    subgroup_dict[this_subgroup_symbol][1].append(other_subgroup_symbol)
                    break # stop checking indices for current other_subgroup

    return subgroup_dict



def clean_subgroup_dict(subgroup_dict: dict, SMILES_string: str) -> dict:
    """
    Input: subgroup_dict should have been created by determine_group_neighbors() and follows the format of
           {group_symbol: [[group_symbol, [group_indices]], [neighbor_list]]}, SMILES_string is a string
    Output: dict of the format {group_symbol: [group_string, [neighbor_symbols]]}
    """

    cleaned_subgroup_dict = {}
    for group_key in subgroup_dict:
        subgroup_info = subgroup_dict[group_key]
        group_neighbor_list = subgroup_info[1]
        group_index_list = subgroup_info[0][1]

        # convert to characters
        group_string = ""
        is_ring = False # Special case for rings, we need a way to identify structures that are meant to represent rings
        for index in group_index_list:
            this_char = SMILES_string[index]

            if this_char in ["(", ")"]:
                continue
            elif this_char.isdigit():   # Don't want to include numbers in the end result, but we do need to include some sort of marker
                is_ring = True
                continue
            else:
                group_string += this_char

        if is_ring:
            # 0 at index 0 will represent
            group_string = "0" + this_char

        cleaned_subgroup_dict[group_key] = [group_string, group_neighbor_list]

    return cleaned_subgroup_dict



def get_group_list(cleaned_subgroup_dict: dict) -> tuple:
    """
    Input: cleaned_subgroup_dict should have been created by clean_subgroup_dict() and follow the format
           {group_symbol: [group_string, [neighbor_symbols]]
    Output: list of all subgroups of the format [group_string_1, group_string_2, ... ]
    Note: subgroups that represent rings will have "0" at index 0 as a marker
    """

    group_list = []
    for group_key in cleaned_subgroup_dict:
        group_info = cleaned_subgroup_dict[group_key]
        group_string = group_info[0]

        if len(group_string) == 0:
            continue

        # add a "0" to the front of the string if it is a ring as a marker
        if type(group_key) == int:
            group_string = "0" + group_string

        group_list.append(group_string)

    return tuple(group_list)



def create_group_combinations(cleaned_subgroup_dict: dict) -> tuple:
    """
    Input: cleaned_subgroup_dict this dictionary should have already been processed by clean_subgroup_dict().
            Format should match {group_symbol: [group_string, [neighbor_symbols]]}
    Output: List of all possible combinations of adjacent subgroups
            Format follows tuple ((group_string_1, group_string_2, ... ]

    Ex Output: (('CCCCCC', 'C'), ('[N+]', 'CCCCCC'), ('[N+]', '=O'), ('[O-]', '[N+]'), ('N', 'CCCCCC'))
    """

    group_combinations = []

    for group_key in cleaned_subgroup_dict:
        this_group_info = cleaned_subgroup_dict[group_key]
        this_group_string = this_group_info[0]
        this_group_neighbors = this_group_info[1]

        if this_group_string in ["", "=", "#"]:     # Unwanted groups
            continue

        for neighbor_key in this_group_neighbors:
            other_group_info = cleaned_subgroup_dict[neighbor_key]

            other_group_string = other_group_info[0]

            # don't include these neighbors that escaped cleaning
            if other_group_string in ["", "=", "#"]:
                continue

            group_combinations.append((this_group_string, other_group_string))

    return tuple(group_combinations)