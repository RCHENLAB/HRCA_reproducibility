# HRCA_reproducibility

This repository comprises code and analyses for the Human Retina Cell Atlas (HRCA) project within the HCA framework, emphasizing the [Eye Biological Network](https://www.humancellatlas.org/biological-networks/).

## File organization for figures and tables

The main code can be found in the [analysis](./analysis) directory of this repository. The subdirectories are organized according to the main figures corresponding to the manuscript. Each figure folder contains a document file that explains the running script, input, and output files. The [utility](./utility) folder comprises standalone scripts that have been used in multiple places.

The repositories below are also developed and utilized in the HRCA.

1. The pipeline to process the unpublished and collected public datasets is accessible at https://github.com/lijinbio/cellqc

2. Scripts related to the benchmark study, integration pipeline, and label transfer using scArches are available at https://github.com/theislab/HRCA-reproducibility

## HRCA Reference Model with scArches

Please refer to [scArches](./scArches) for a brief tutorial on training the reference model for HRCA v1.0 and performing label transfer using the trained reference model with scArches.

The HRCA v1.0 reference model with scArches is available for download on Zenodo at https://doi.org/10.5281/zenodo.14014720.

## HRCA on Interactive Browsers

HRCA can be accessed through several interactive web browsers, including:

- [HCA Data Portal](https://data.humancellatlas.dev.clevercanary.com/hca-bio-networks/eye)
- [CELLxGENE](https://cellxgene.cziscience.com/collections/4c6eaf5c-6d57-4c76-b1e9-60df8c655f1e)
- [UCSC Cell Browser](https://retina.cells.ucsc.edu)
- [Broad Single Cell Portal](https://singlecell.broadinstitute.org/single_cell/study/SCP2805)
- [Cell Annotation Platform](https://celltype.info/project/381)

## References

The HRCA paper: [Li et al.](https://doi.org/10.1101/2023.11.07.566105), bioRxiv 2023.

## Questions

If you have any questions, please submit an issue to this repository.

