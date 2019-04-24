# MutationExplorer

MutationExplorer is a lightweight web application for exploring somatic mutations in cancer datasets.  It is based on R and the Shiny framework.

MutationExplorer is currently most usable for breast cancer samples, since the relevant biomarkers for patient selection (e.g. clinical status of ER and PgR) are hardcoded.


## Features

- Filtering of mutations by type (nonsense, missense, etc) and [COSMIC](https://cancer.sanger.ac.uk/cosmic/) status
- Filtering of the patient sample cohort by histological/genomic biomarker status
- Mutational landscape waterfall plots
- Kaplan-Meier survival analysis based on gene and pathway mutation status, as well as tumor mutation burden (TMB)
- Protein plots similar to cBioPortal oncoprints or lollipop plots


## Authors

- Christian Brueffer ([@cbrueffer](http://github.com/cbrueffer/))
- Sergii Gladchuk ([@sergiigladchuk](http://github.com/sergiigladchuk/))

This software was developed as part of a PhD research project in the
[laboratory of Lao H. Saal, Translational Oncogenomics Unit, Department of Oncology and Pathology, Lund University, Sweden](https://www.med.lu.se/saalgroup).


## Availability

A MutationExplorer instance with the data from a large breast cancer patient cohort together with an accompanying pre-print will be available shortly.
