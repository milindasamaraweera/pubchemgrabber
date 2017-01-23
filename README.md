#pubchemgrabber
pubchemgrabber is an automated-LC-MS based work-flow written in python, intended to be used in untargated-metabolomics. 

pubchemgrabber allows...
* Retrieval of candidate sets using exact mass from PubChem remotely or a local SD file
* Apply a set of pre-filters to reduce the retrieved candidate sets, based on the chemical knowledge of the unknown (i.e charge)
* Apply an Artifical-Neural-Network (ANN) based predictive model to calculate RI values of each candidate and match it with the experimental RI value to reduce the number of candidates
* Prepare input files to be used in in-silico fragmenters such as: magmaplus and cfm-id
