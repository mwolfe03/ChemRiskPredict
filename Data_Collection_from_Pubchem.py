import requests
import time



def fetch_compound_data(compound_ID: int, wait_time: float) -> dict:
    """
    Input: compound_ID is an int representation of the compound ID. wait_time is a float value of the desired
    wait time between requests in seconds
    Output: dictionary containing compound info including Canonical Smiles, and Hazards
    """

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{compound_ID}/JSON/"
    response = requests.get(url)

    # cool down inbetween calls
    time.sleep(wait_time)

    data = response.json()

    iupac_name = None
    canonical_smiles = None
    hazards = []

    # Extracting IUPAC name and canonical SMILES
    for section in data['Record']['Section']:
        if section['TOCHeading'] == 'Names and Identifiers':
            for subsection in section['Section']:
                if subsection['TOCHeading'] == 'Computed Descriptors':
                    for descriptor in subsection['Section']:
                        if descriptor['TOCHeading'] == 'IUPAC Name':
                            iupac_name = descriptor['Information'][0]['Value']['StringWithMarkup'][0]['String']
                        elif descriptor['TOCHeading'] == 'Canonical SMILES':
                            canonical_smiles = descriptor['Information'][0]['Value']['StringWithMarkup'][0]['String']

    # Extracting hazards
    for section in data['Record']['Section']:
        if section['TOCHeading'] == 'Safety and Hazards':
            for subsection in section['Section']:
                if subsection['TOCHeading'] == 'Hazards Identification':
                    for subsubsection in subsection['Section']:
                        if subsubsection['TOCHeading'] == 'GHS Classification':
                            for info in subsubsection['Information']:
                                if info['Name'] == 'Pictogram(s)':
                                    for pictogram in info['Value']['StringWithMarkup'][0]['Markup']:
                                        hazards.append(pictogram['Extra'])

    return {
        "Compound ID": compound_ID,
        "IUPAC Name": iupac_name,
        "Canonical SMILES": canonical_smiles,
        "Hazards": ', '.join(hazards)
    }