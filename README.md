# attractiveness_penalties

----------------------------------------------------------------------------------

Title: Leveraging Genomic Data to Document Within-Race Attractiveness Penalties Among Black Americans

Journal: Sociological Science

Authors: Taddess, Zhang, and Trejo

*Note, Taddess and Zhang are co-first authors and Trejo is the corresponding author

Date: May 2026

----------------------------------------------------------------------------------

This replication package contains code to reproduce the analyses, tables, and figures included in the main text and supplementary material. The Add Health survey and genetic data used in this project are restricted-use and cannot be redistributed with the replication code. Researchers must apply separately for access through Add Health/CPC and dbGaP.

Code used to construct the genetic similarity proportions (GSPs) used in this project is available at: https://github.com/luyin-z/estimating_admixture

----------------------------------------------------------------------------------

Description of the script files: 


1. 01-prepare_pheno.R: prepare survey data and merge them with GSPs
2. 02-analysis.R: run the main analyses and generate tables and figures
3. 03-ref_map.R: create the map showing the geographic locations of reference populations used to estimate GSPs

----------------------------------------------------------------------------------

Description of the data files: 


1. igsr_populations.tsv: coordinates information for the reference panel, downloaded from
https://www.internationalgenome.org/data-portal/population
2. reference_composition.txt: composition of the reference panel
