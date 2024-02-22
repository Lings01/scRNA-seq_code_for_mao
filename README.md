# Single-Cell Analysis for CD276-dependent Efferocytosis in Bladder Cancer

This repository contains the single-cell RNA sequencing analysis code used in the research paper "CD276-dependent efferocytosis by tumor-associated macrophages promotes immune evasion in bladder cancer". The code herein is intended to process, analyze, and visualize single-cell transcriptomics data to elucidate the role of CD276 in the efferocytosis process conducted by tumor-associated macrophages (TAMs) and its implications in immune evasion within the bladder cancer microenvironment.

## Repository Contents

- `cellchat.R`: Code for performing cell-cell communication analysis using CellChat.
- `cellphoneDB.R`: Scripts for statistical analysis of cell-cell communication via CellPhoneDB.
- `gsva.R`: Gene set variation analysis (GSVA) for pathway enrichment analysis.
- `harmony.R`: Harmony algorithm for integrating multiple single-cell datasets.
- `scenic.R`: SCENIC pipeline for reconstructing gene regulatory networks.

## Usage

Each script is named according to the analysis or tool used. To reproduce the results from the paper, run each script in the R environment. Ensure that all dependencies are installed before executing the scripts.

## Dependencies

- R (version 3.6.0 or higher)
- CellChat
- CellPhoneDB
- GSVA
- Harmony
- SCENIC
- Additional R packages as required by the above

## Installation

To run the analysis scripts, you will need to install the necessary R packages. Installation instructions for each package can be found at their respective websites or repositories.

## Data

The data used for this analysis is not included in this repository and must be obtained separately.

## Contact

For questions or support, please open an issue in this repository or contact the corresponding author of the research paper.

## Citation

If you use the code or methodologies in this repository in your research, please cite our paper:

> [Author(s)], "CD276-dependent efferocytosis by tumor-associated macrophages promotes immune evasion in bladder cancer," [Journal], [Year].
