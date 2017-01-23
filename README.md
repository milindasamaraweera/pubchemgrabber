#pubchemgrabber
An automated work-flow intended to be used in untargated-metabolomics to retrieve and filter candidate by matching experimental RI values (obtained using LC-MS experiments) with predicted RI values based on an ANN based RI predictor.

pubchemgrabber allows...
* Retrieval of candidate using exact mass from PubChem or a local SD file
* Apply a set of pre-filters to reduce the candidate sets, based on the chemical knowledge of the unknown (i.e charge)
* Apply a ANN based RI Filter to reduce the candidate sets further based on experimental RI value of the unknown
