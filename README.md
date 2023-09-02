# Predicting Inflation with Generalized Random Forests

## Overview of scripts, data files and folders

### Scripts

| File name    | Description |
| :-------- | :------- |
| main_script.R  | runs the full pseudo-out-of-sample procedure and exports csv files containing results   |
| main_script_functions.R | includes a number of support functions, including functions cleaning/preparing the datasets and dimensionality reduction functions  |
| pca.R   | contains the principal components estimation function  |
| analysis_main_script.R  | turns the output from main_script.R into plots and summary statistics |
| get_grf.R  | includes the GRF forest functions (standard RF and LLF)  |
| get_rf_gls.R   | includes the RF-GLS forest function  |
| parallel_get_grf.R  | includes an “add-on” function utilizing parallel computing  |
| parallel_get_rf_gls.R  | includes an “add-on” function utilizing parallel computing  |
| data_transformations.R   | obtains the results from the ADF-KPSS testing algorithm  |
| data_transformation_functions.R  | contains the functions implementing the ADFKPSS testing algorithm|
| simulation_study.R   | runs simulation study including predictions and confidence intervals produced for Friedman’s DGP and exports csv files with results  |
| simulation_studies_analysis.R   | obtains summary statistics based on the output from simulation_study.R  |

### Data files

| File name    | Description |
| :-------- | :------- |
| dataset_friedmans_d10_s5_scalingfalse.xlsx  | contains data simulated from Friedman’s DGP |

### Folders

| File name    | Description |
| :-------- | :------- |
| main_results   | results from main_script.R are exported into this folder |
| simulation_results_ci_regression_forest  | results from simulation_study.R are exported into this folder |
| simulation_results_ci_llf | results from simulation_study.R are exported into this folder |
| simulation_results_ci_rf_gls  | results from simulation_study.R are exported into this folder |

## Instructions to run simulation study 
The file named simulation_study.R runs a simulation study comparing self-written forest functions to established packages. To run the simulation study, complete the following steps. 

* Step 1: Save all files of the repository in a folder on your computer.
* Step 2: Install needed packages by running the following command in R:

```R
install.packages(c("openxlsx", "powerplus", "Matrix", "doParallel", "foreach", "readxl", "urca", "grf", "randomForest"))
```
* Step 3: Open the file named simulation_study.R and set your working directory in line 4.
* Step 4: Run the full script. 
