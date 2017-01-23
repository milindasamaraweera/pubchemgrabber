#pubchemgrabber
pubchemgrabber is an automated-LC-MS based work-flow written in python, intended to be used in untargated-metabolomics. 

pubchemgrabber allows...
* Retrieval of candidate using exact mass from PubChem or a local SD file
* Apply a set of pre-filters to reduce the candidate sets, based on the chemical knowledge of the unknown (i.e charge)
* Apply a ANN based RI Filter to reduce the candidate sets further based on experimental RI value of the unknown
