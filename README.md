# Predicting Inflation with Generalized Random Forests

This repository contains R code accompanying my master thesis titled Predicting Inflation with Generalized Random Forests. Beside various testing procedures, the repository contains self-written implementations from scratch of the standard Random Forest (Breiman, 2001) written in the Generalized Random Forest framework (Athey et. al., 2019), the Local Linear Forest (Friedberg et. al., 2021) and the RF-GLS for spatially dependent data (Saha et. al., 2021). 

Instructions for how to run a simulation study comparing the self-written forest functions to established forest packages are found at the bottom of this README file, along with instructions for how to download and use the self-written forest functions separately. The script needed to obtain the main results of my thesis, main_script.R, is included in the repository, but is not runnable. That is because the macro-financial dataset used is currently not included in the repository, as it is an internal ECB dataset, which is not publicly available. 

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

## References 
Athey, S., Tibshirani, J., and Wager, S. (2019). Generalized Random Forests. The Annals of Statistics, 47:1148-1178.

Breiman, L. (2001) Random Forests. Machine Learning, 45:5-32.

Friedberg, R., Tibshirani, J., Athey, S., and Wager, S. (2021). Local Linear Forests. Journal of Computational and Graphical Statistics, 30(2):503-517.

Saha, A., Basu, S., and Datta, A. (2021). Random Forests for Spatially Dependent Data. Journal of the American Statistical Association. 
