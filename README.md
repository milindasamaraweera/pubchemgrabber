#pubchemgrabber
pubchemgrabber is an automated-LC-MS based work-flow written in python, intended to be used in untargated-metabolomics. 

pubchemgrabber allows...
* Retrieval of candidate sets using exact mass from PubChem remotely or from a local SD file
* Apply a set of pre-filters to reduce the size of retrieved candidate sets, based on the prior chemical knowledge of the unknown (i.e charge, mode of operation of the mass-spec)
* Apply Artifical-Neural-Network (ANN) based predictive model to calculate a RI value of each candidate and match it with the experimental RI value, finally removing candidates that have RI values that greatly vary compared to the experimental RI value
* Prepare input files to be used in in-silico fragmenters such as: magmaplus and cfm-id
