#pubchemgrabber
pubchemgrabber is an automated-LC-MS based work-flow written in python (version 3.5), intended to be used in untargated-metabolomics. 

pubchemgrabber allows...
* Retrieval of candidate sets using exact mass from PubChem remotely or a local SD file
* A pre-screening of the candidate sets (based on the knowledge of the unknown, i.e charge)
* RI based screening of candidate sets using an Artifical-Neural-Network (ANN) based predictive model
* Prepare input files to be used in in-silico fragmenters such as: magmaplus and cfm-id
