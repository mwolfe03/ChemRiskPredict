<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ChemRiskPredict README</title>
</head>
<body>
    <h1>ChemRiskPredict</h1>
    <p><strong>ChemRiskPredict</strong> is a Python package designed to predict the potential hazards of chemical molecules based on their Canonical SMILES structures. This package leverages cheminformatics and machine learning techniques to provide accurate hazard predictions. The data used for training the models is collected from PubChem's PUG-REST API.</p>

    <h2>Features</h2>
    <ul>
        <li>Predict various hazards from Canonical SMILES structures.</li>
        <li>Update existing dataframes with new hazard predictions.</li>
        <li>Retrieve data directly from PubChem using their PUG-REST API.</li>
        <li>Handle molecular rings using RDKit.</li>
    </ul>

    <h2>Installation</h2>
    <p>To install ChemRiskPredict, you can use pip:</p>
    <pre><code>pip install git+https://github.com/mwolfe03/ChemRiskPredict.git</code></pre>

    <h2>Dependencies</h2>
    <p>The following Python packages are required:</p>
    <ul>
        <li>pandas</li>
        <li>scikit-learn</li>
        <li>rdkit</li>
        <li>requests</li>
    </ul>
    <p>These dependencies will be installed automatically when you install ChemRiskPredict.</p>

    <h2>Usage</h2>
    <p>Here is an example of how to use ChemRiskPredict:</p>
    <h3>Importing the Package</h3>
    <pre><code>import ChemRiskPredict as crp</code></pre>

    <h3>Creating a DataFrame from CIDs</h3>
    <pre><code>cids = [2244, 3672, 5281]  # Example PubChem Compound IDs
dataframe = crp.create_dataframe_from_cids(cids)
print(dataframe)</code></pre>

    <h3>Predicting Hazards</h3>
    <pre><code># Assuming you have a dataframe with a 'SMILES' column
hazard_predictions = crp.determine_potential_hazards_from_dataframe(dataframe)
print(hazard_predictions)</code></pre>

    <h3>Updating an Existing DataFrame</h3>
    <pre><code>new_cids = [5357, 6256]
updated_dataframe = crp.update_existing_ids_dataframe_from_cids(dataframe, new_cids)
print(updated_dataframe)</code></pre>

    <h2>Functions</h2>
    <h3>create_dataframe_from_cids(cids)</h3>
    <p>Create a pandas DataFrame from a list of PubChem Compound IDs (CIDs).</p>

    <h3>determine_potential_hazards_from_dataframe(dataframe)</h3>
    <p>Determine potential hazards from a DataFrame containing SMILES structures.</p>

    <h3>update_existing_dataframe_from_dataframe(existing_df, new_df)</h3>
    <p>Update an existing DataFrame with data from another DataFrame.</p>

    <h3>update_existing_ids_dataframe_from_cids(existing_df, cids)</h3>
    <p>Update an existing DataFrame with new data retrieved using PubChem CIDs.</p>

    <h2>Data Files</h2>
    <p>You can access additional data files included in the package using the <code>get_data_file</code> function:</p>
    <pre><code>file_path = crp.get_data_file('data_file.csv')</code></pre>

    <h2>License</h2>
    <p>This project is licensed under the MIT License - see the <a href="LICENSE">LICENSE</a> file for details.</p>

    <h2>Contact</h2>
    <p>For any inquiries or issues, please contact the author:</p>
    <ul>
        <li><strong>Author:</strong> Mateo Wolfe</li>
        <li><strong>Email:</strong> mwolfe08@calpoly.edu</li>
    </ul>

    <h2>Contributing</h2>
    <p>Contributions are welcome! Please fork the repository and submit a pull request with your changes. Make sure to include tests for any new features or bug fixes.</p>

    <h2>Acknowledgments</h2>
    <p>Data is collected from PubChem's PUG-REST API: <a href="https://pubchem.ncbi.nlm.nih.gov/rest/pug_view">https://pubchem.ncbi.nlm.nih.gov/rest/pug_view</a></p>
</body>
</html>
