# Cell-Cell Ratio in IMC Data of LUAD Patients

## Introduction
This study is based on research published by Sorin et al. titled "Single-cell Spatial Landscapes of the Lung Tumor Immune Microenvironment" (available at [PubMed](https://pubmed.ncbi.nlm.nih.gov/36725934/)), which provides data from 416 LUAD patients, including IMC (Imaging Mass Cytometry) data and clinical information. The downloading and processing of raw IMC data are detailed in the original paper.

The contents of this repository comprise:

1. Cell density, cell fraction, rescaled immune cell fraction, and cell-cell ratio calculated for each IMC sample in this study.
2. R code for survival analysis and down-sampling analysis.

## File Descriptions

### Cell Density/Fraction/Ratio Supplementary Tables
1. **Cell_Fraction.csv**: This file contains the proportions of each of the 16 types of cells within each IMC sample, relative to the total number of cells. The sum of cell fractions in each sample equals 1.

2. **Cell_Fraction_rescaled_after_removed_cancer_endothelial.csv**: Similar to the previous file, this one presents the proportions of each cell type within the remaining cell population after removing Endothelial and cancer cells. The sum of the immune cell fractions in each sample equals 1.

3. **Cell_Density.csv**: This file provides the counts of each type of cell per unit area of the image within the IMC sample. The unit area is calculated from the length * width divided by 1e6 (1,000,000).

4. **CellCell_Ratio.csv**: In this file, the fractions of the 16 types of cells are presented as ratios calculated by combining the numerator and denominator, resulting in a total of 240 permutations.

### R Data and Code
1. **IMC_LUAD.RData**: This data file contains two components: "data," which encompasses all the tables mentioned above, and "info," which contains clinical data for the 416 patients.

2. **CellCell_Ratio.csv**: This file includes R code for conducting survival analysis and down-sampling validation using the aforementioned data.

Please feel free to explore and utilize these resources for your research and analysis.
