from setuptools import setup, find_packages

# Package metadata
NAME = 'ChemRiskPredict'
VERSION = '1.0.0'
DESCRIPTION = 'A package for predicting chemical hazards.'
AUTHOR = 'Mateo Wolfe'
EMAIL = 'mwolfe08@calpoly.edu'
URL = 'https://github.com/mwolfe03/ChemRiskPredict'
LICENSE = 'MIT'
KEYWORDS = ['chemistry', 'risk prediction', 'hazard prediction, cheminformatics']


# Define dependencies
INSTALL_REQUIRES = [
    'pandas',
    'scikit-learn',
    'rdkit',
    'requests',
]

setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description= "Predicts varius hazards from Canonical SMILES structures",
    long_description_content_type='text/markdown',
    author=AUTHOR,
    author_email=EMAIL,
    url=URL,
    license=LICENSE,
    keywords=KEYWORDS,
    packages=find_packages(),
    install_requires=INSTALL_REQUIRES,
    python_requires='>=3.6',
)
