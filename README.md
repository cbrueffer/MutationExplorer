# ShinyMutationExplorer

ShinyMutationExplorer is a lightweight web application for exploring somatic mutations in cancer datasets.  It is based on R and the Shiny framework.

The application is currently most usable for breast cancer samples, since the relevant biomarkers for patient selection (clinical status ER, PgR, ...) are hardcoded.

This software was developed as part of a PhD research project in the
[laboratory of Lao H. Saal, Translational Oncogenomics Unit, Department of Oncology and Pathology, Lund University, Sweden](https://www.med.lu.se/saalgroup).


## TODO

- don't allow plot download when no plot available (e.g. no gene selected)
- only redrawing when more genes are being selected, see https://code.i-harness.com/en/q/1db7c1c for an example
- Add color palette selection from https://github.com/gabrifc/raincloud-shiny/blob/master/source/paletteColours.R
