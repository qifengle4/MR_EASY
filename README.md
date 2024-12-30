
# MR_EASY
This package provides a set of custom functions for conducting Mendelian Randomization (MR) analysis, including process_gwas_data, perform_mr_analysis, and perform_re_mr_analysis. It primarily achieves this by processing genome-wide association study (GWAS) data. The code focuses on two main functionalities: data preprocessing and MR analysis.
#Overview of Main Functionality
1. process_gwas_data Function
1.1 Functionality
This function is responsible for preprocessing GWAS raw data, which mainly includes reading, formatting, filtering, local clustering, and F-statistic calculation.
1.2 Input Parameters
(1) rawdata_path: The path pointing to the original GWAS data file.
(2) output_no_filter_path: The output path for storing the unfiltered GWAS data file.
(3) output_clumped_path: The output path for storing results after local clustering.
(4) output_f_filter_path: The output path for storing results after filtering based on F-statistics.
(5) work_path: Sets the working directory.
(6) file_ids: A vector of GWAS file IDs to be processed.
(7) file_suffix: A custom suffix for the input file, generally used to indicate the file format.
(8) pFilter: The p-value threshold used for filtering SNPs.
(9) plink_path: The path to the PLINK executable, used for SNP clustering.
(10) bfile_path: The path to the reference genome data used in the clustering process.
(11) metadata_file: The file path to store metadata associated with GWAS (e.g., sample size).
(12) metadata_id_col, metadata_sample_size_col: The column names specifying the sample identifier and sample size in the metadata.
(13) type: The data type (e.g., 'exposure' or 'outcome').
(14) Other parameters include names of relevant columns such as SNP, beta, standard error, effect, etc.
1.3 Key Steps
(1) Directory Check and Creation: Checks whether the input/output directories exist and creates them if they do not.
(2) Data Reading and Formatting:
Uses the fread function to read the original GWAS data.
Adds sample size and phenotype columns.
Formats the data to comply with the requirements of the TwoSampleMR package.
(3) Data Output: Writes the unfiltered data to the specified output directory.
(4) P-value Filtering and Clustering:
Filters the data according to the defined p-value to ensure the validity of selected SNPs.
Uses ieugwasr::ld_clump to cluster SNPs, employing PLINK to remove highly correlated SNPs.
(5) F-statistic Calculation: Calculates the F-statistic for the clustered SNPs and applies further filtering (e.g., using thresholds for MAF and F values).
(6) The final data output is stored in a specified path.

2. perform_mr_analysis Function
2.1 Functionality
Conducts Mendelian randomization analysis by combining exposure and outcome data for statistical testing and generating a series of related visualizations.
2.2 Input Parameters
(1) exposure_data_path: The path to the processed exposure data.
(2) outcome_data_path: The path to the outcome data.
(3) result_output_path: The directory to store results.
(4) file_list: A list of file IDs to analyze.
(5) outcome_name: A vector specifying the outcome name.
(6) pFilter: The p-value threshold for filtering.
(7) exposure_file_suffix, outcome_file_suffix: Custom file suffixes to help determine the file types being read.
2.3 Key Steps
(1) Directory Check: Ensures that the exposure data, outcome data, and output paths exist.
(2) Loop Processing Each Result:
Reads the exposure and outcome data one by one and formats them according to specified columns.
Merges exposure and outcome data to ensure data integrity.
(3) MR Analysis:
Uses harmonise_data to align the output data.
Applies MR-PRESSO for outlier detection and visualization.
(4) Statistical Testing:
Calculates the results of Mendelian randomization analysis, including methods like IVW (Inverse Variance Weighted), BWMR (Bayesian Weighted MR), MR-Egger, weighted median methods, Simple mode, and Weighted mode.
Conducts heterogeneity and pleiotropy tests.
FDR correction for positive results from the IVW method.
(5) Visualization of Results:
Generates scatter plots, forest plots, and funnel plots to visualize the analysis results, facilitating data and result comprehension.
(6) Output Results:
Stores the analysis result and relevant test results in CSV format at the specified path, ensuring easy access for subsequent use and review.

3. perform_re_mr_analysis Function
3.1 Functionality
Conducts retrospective analysis of previous Mendelian randomization results, mainly utilizing existing IVW results for further verification and refinement.
3.2 Input Parameters
(1) exposure_data_path: The path to the exposure data.
(2) outcome_data_path: The path to the outcome data.
(3) result_output_path: The output path for storing results.
(4) mr_results_path: The path to existing MR results.
(5) outcome_name: A list of outcome names.
(6) pFilter: The p-value for filtering.
(7) exposure_file_suffix, outcome_file_suffix: The suffixes to identify files.
3.3 Key Steps
(1) Results Directory Check: Ensures that all necessary folders and paths exist.
(2) File List Generation: Extracts the file list from existing IVW results, ensuring the continuity and consistency of the analysis.
(3) Reading Exposure and Outcome Data: Loops through and reads specified exposure and outcome data, performing data integration and validation.
(4) MR Analysis:
Similar to perform_mr_analysis, performs data merging and randomization analysis, but focuses more on examining existing IVW results.
(5) Testing and Output: Cleans up data frames after outputting analysis results to ensure clear state management for each iteration.
(6) Final Results Storage: Stores results from all different analysis methods and cleaned data for easier follow-up analyses.
