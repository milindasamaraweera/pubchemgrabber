#pubchemgrabber
An automated work-flow to be using in untargated-metabolomics to retrieve and filter candidate sets based on experimental RI values obtained using LC-MS experiments

pubchemgrabber allows...
* Retrieval of candidate using exact mass from PubChem or a local SD file
* Apply a set of pre-filters to reduce the candidate sets, based on the chemical knowledge of the unknown (i.e charge)
* Apply a ANN based RI Filter to reduce the candidate sets further based on experimental RI value of the unknown
