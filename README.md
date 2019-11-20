# MutationExplorer

MutationExplorer is a lightweight web application for exploring somatic mutations in cancer datasets.  It is based on R and the Shiny framework.

MutationExplorer is currently most usable for breast cancer samples, since the relevant biomarkers for patient selection (e.g. clinical status of ER and PgR) are hardcoded.

<p align="center">
<img src="https://github.com/cbrueffer/MutationExplorer/blob/master/doc/images/mutationexplorer-main.png" width="600" height="330" />
</p>


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

A MutationExplorer instance with the data from 3,217 primary breast tumor transcriptomes from the [Sweden Cancerome Analysis Network--Breast (SCAN-B)](https://portal.research.lu.se/portal/en/projects/sweden-cancerome-analysis-network--breast-scanb-a-largescale-multicenter-infrastructure-towards-implementation-of-breast-cancer-genomic-analyses-in-the-clinical-routine%2898b7c09f-85d5-4f98-bcae-a43bffb6f869%29/publications.html) can be accessed at https://oncogenomics.bmc.lu.se/MutationExplorer.


## Citation

A preprint describing this tool and the dataset mentioned above will be available shortly. For now, please cite the SCAN-B MutationExplorer website: https://oncogenomics.bmc.lu.se/MutationExplorer
